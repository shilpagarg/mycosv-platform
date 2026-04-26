#!/usr/bin/env python3
# Designed for Linux

from __future__ import annotations

import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Any


DEFAULT_TOL_BP = {
    "INS": 500,
    "DEL": 2500,
    "DUP": 2500,
    "INV": 10000,
    "TRA": 10000,
    "TRA_INTER": 10000,
    "TRA_INTRA": 10000,
    "OFF_REF": 500,
}

DEFAULT_TOL_LEN_FRAC = {
    "INS": 0.30,
    "DEL": 0.10,
    "DUP": 0.10,
    "INV": 1.0,
    "TRA": 1.0,
    "TRA_INTER": 1.0,
    "TRA_INTRA": 1.0,
    "OFF_REF": 1.0,
}


def parse_info_field(field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in field.split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        out[key] = value
    return out


def normalise_chrom(chrom: str) -> str:
    return re.sub(r"_ins\d*$", "", chrom or "")


def normalise_truth_chrom(chrom: str) -> str:
    chrom = re.sub(r"__sv_.*$", "", chrom or "")
    return re.sub(r"_ins\d*$", "", chrom)


def is_hint_contig(raw_chrom: str) -> bool:
    return "__sv_" in (raw_chrom or "")


def near_window_boundary(pos: int, window_bp: int, tol: int) -> bool:
    if window_bp <= 0:
        return False
    remainder = pos % window_bp
    return remainder < tol or (window_bp - remainder) < tol


def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return 0.0, 0.0
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return max(0.0, centre - half), min(1.0, centre + half)


_ASM_EXT_STRIP = (
    ".fasta.gz", ".fastq.gz", ".fna.gz", ".fa.gz", ".fq.gz",
    ".fasta", ".fastq", ".fna", ".fa", ".fq",
)


def _asm_aliases(asm: str) -> list[str]:
    """Return canonical alias forms for a query-asm name.

    The C++ binary derives QASM from the input file path's stem, which strips
    one extension (".fa" but not ".fa.gz"). Metadata TSVs sometimes carry the
    bare name and sometimes the full filename. Try the bare name and several
    extension-stripped variants so per-scenario lookups don't silently fall
    back to "unknown" purely because of a naming convention mismatch.
    """
    out = [asm]
    seen = {asm}
    lower = asm.lower()
    for ext in _ASM_EXT_STRIP:
        if lower.endswith(ext):
            stripped = asm[: -len(ext)]
            if stripped and stripped not in seen:
                seen.add(stripped)
                out.append(stripped)
    return out


def load_query_asm_maps(pred_hits_tsv: Path | None,
                        query_meta_tsv: Path | None
                        ) -> tuple[dict[str, str], dict[tuple[str, str], dict[str, str]], dict[str, str]]:
    contig_to_asm: dict[str, str] = {}
    call_context: dict[tuple[str, str], dict[str, str]] = {}
    asm_to_scenario: dict[str, str] = {}

    if pred_hits_tsv and pred_hits_tsv.exists():
        with pred_hits_tsv.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                contig = row.get("query_contig", "")
                asm = row.get("query_asm", "")
                if contig and asm:
                    contig_to_asm[contig] = asm
                    call_context[(asm, contig)] = {
                        "ref_contig": row.get("ref_contig", ""),
                        "ref_asm": row.get("ref_asm", ""),
                        "query_mode": row.get("query_mode", ""),
                    }

    if query_meta_tsv and query_meta_tsv.exists():
        with query_meta_tsv.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                asm = row.get("query_asm", "")
                scenario = row.get("scenario", "unknown")
                if not asm:
                    continue
                # Index every alias so a binary-emitted QASM that includes a
                # file extension (e.g. "core_asm3.fa") still resolves to the
                # bare-name scenario row in query_metadata.tsv.
                for alias in _asm_aliases(asm):
                    asm_to_scenario.setdefault(alias, scenario)

    return contig_to_asm, call_context, asm_to_scenario


def parse_vcf_records(path: Path,
                      *,
                      is_pred: bool,
                      contig_to_asm: dict[str, str],
                      call_context: dict[tuple[str, str], dict[str, str]],
                      asm_to_scenario: dict[str, str],
                      window_bp: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not path.exists():
        return rows

    with path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            info = parse_info_field(fields[7])
            raw_chrom = fields[0]
            chrom = normalise_chrom(raw_chrom) if is_pred else normalise_truth_chrom(raw_chrom)
            qasm = info.get("QASM", info.get("QUERY_ASM", ""))
            if is_pred and not qasm:
                qasm = contig_to_asm.get(chrom, "") or contig_to_asm.get(raw_chrom, "")
            pred_ctx = call_context.get((qasm, raw_chrom), {}) if is_pred and qasm else {}
            scenario = info.get("SCENARIO", "")
            if is_pred and not scenario and qasm:
                # Try the qasm verbatim, then stem variants — the binary's
                # QASM may carry a file extension while query_metadata.tsv
                # is keyed on the bare asm name.
                for alias in _asm_aliases(qasm):
                    hit = asm_to_scenario.get(alias)
                    if hit:
                        scenario = hit
                        break
            svtype = info.get("SVTYPE", "?")
            pos_i = int(fields[1])
            tol_bp = int(info["TOL_BP"]) if "TOL_BP" in info else DEFAULT_TOL_BP.get(svtype, 2500)
            tol_len_frac = (float(info["TOL_LEN_FRAC"])
                            if "TOL_LEN_FRAC" in info
                            else DEFAULT_TOL_LEN_FRAC.get(svtype, 0.10))
            raw_chr2 = info.get("CHR2", ".")
            chr2 = normalise_chrom(raw_chr2) if is_pred else normalise_truth_chrom(raw_chr2)
            rows.append({
                "chrom": chrom,
                "raw_chrom": raw_chrom,
                "pos": pos_i,
                "type": svtype,
                "end": int(info.get("END", fields[1])),
                "svlen": int(info.get("SVLEN", 0)),
                "annot": info.get("ANNOT", "NONE"),
                "clade": info.get("CLADE", ""),
                "scenario": scenario or "unknown",
                "qasm": qasm,
                "qmode": info.get("QMODE", pred_ctx.get("query_mode", "")),
                "ref_contig": normalise_chrom(info.get("REFCONTIG", pred_ctx.get("ref_contig", ""))),
                "ref_asm": info.get("REFASM", pred_ctx.get("ref_asm", "")),
                "chr2": chr2,
                "pos2": int(info.get("POS2", 0) or 0),
                "end2": int(info.get("END2", 0) or 0),
                "tol_bp": tol_bp,
                "tol_len_frac": tol_len_frac,
                "hint_driven": is_pred and is_hint_contig(raw_chrom),
                "near_boundary": near_window_boundary(pos_i, window_bp, tol_bp),
            })
    return rows


def compatible(truth: dict[str, Any], pred: dict[str, Any]) -> bool:
    if truth.get("qasm") and pred.get("qasm") and truth["qasm"] != pred["qasm"]:
        return False

    del_group = {"DEL", "INS", "DUP", "TDEL", "TANDEM_DUP"}
    inv_group = {"INV", "TRA", "TRA_INTER", "TRA_INTRA"}
    if truth["type"] != pred["type"]:
        if not ((truth["type"] in del_group and pred["type"] in del_group) or
                (truth["type"] in inv_group and pred["type"] in inv_group)):
            return False

    pred_mode = pred.get("qmode", "") or "assembly"
    pred_locus = pred["chrom"]
    if pred_mode != "assembly" and truth["type"] != "OFF_REF":
        pred_locus = pred.get("ref_contig") or pred_locus
    if truth["chrom"] != pred_locus and not (pred_mode != "assembly" and truth["type"] == "OFF_REF"):
        return False

    ins_group = {"INS", "TANDEM_DUP", "OFF_REF"}
    tol = truth["tol_bp"]
    both_positionless = truth["type"] in ins_group and pred["type"] in ins_group
    if not both_positionless and abs(truth["pos"] - pred["pos"]) > tol:
        return False

    skip_len = {"INV", "TRA", "TRA_INTER", "TRA_INTRA", "TANDEM_DUP", "INS", "TDEL", "OFF_REF"}
    if truth["type"] not in skip_len:
        denom = max(abs(truth["svlen"]), 1)
        if abs(abs(truth["svlen"]) - abs(pred["svlen"])) / denom > truth["tol_len_frac"]:
            return False

    if truth["type"] in {"TRA", "TRA_INTER", "TRA_INTRA"}:
        t_mate = truth.get("chr2") not in (None, "", ".") and truth.get("pos2", 0) > 0
        p_mate = pred.get("chr2") not in (None, "", ".") and pred.get("pos2", 0) > 0
        if t_mate and p_mate:
            if truth["chr2"] != pred["chr2"]:
                return False
            if abs(truth["pos2"] - pred["pos2"]) > tol:
                return False

    return True


def distance(truth: dict[str, Any], pred: dict[str, Any]) -> int:
    pos_d = abs(truth["pos"] - pred["pos"])
    len_d = 0 if truth["type"] in {"TRA", "TRA_INTER", "TRA_INTRA"} else abs(abs(truth["svlen"]) - abs(pred["svlen"]))
    mate_d = 0
    if truth["type"] in {"TRA", "TRA_INTER", "TRA_INTRA"} and truth.get("pos2", 0) and pred.get("pos2", 0):
        mate_d = abs(truth["pos2"] - pred["pos2"])
    return pos_d + len_d + mate_d


def match_truth_to_pred(truth_list: list[dict[str, Any]],
                        pred_list: list[dict[str, Any]]) -> tuple[set[int], list[dict[str, Any]]]:
    used: set[int] = set()
    fn_list: list[dict[str, Any]] = []
    for truth in truth_list:
        best_idx: int | None = None
        best_dist = math.inf
        for idx, pred in enumerate(pred_list):
            if idx in used or not compatible(truth, pred):
                continue
            dist = distance(truth, pred)
            if dist < best_dist:
                best_idx = idx
                best_dist = dist
        if best_idx is None:
            fn_list.append(truth)
        else:
            used.add(best_idx)
    return used, fn_list


def score_pr(truth_vcf: Path,
             pred_vcf: Path,
             *,
             pred_hits_tsv: Path | None = None,
             query_meta_tsv: Path | None = None,
             window_bp: int = 2_000_000) -> dict[str, Any]:
    contig_to_asm, call_context, asm_to_scenario = load_query_asm_maps(pred_hits_tsv, query_meta_tsv)
    truth = parse_vcf_records(
        truth_vcf,
        is_pred=False,
        contig_to_asm=contig_to_asm,
        call_context=call_context,
        asm_to_scenario=asm_to_scenario,
        window_bp=window_bp,
    )
    pred = parse_vcf_records(
        pred_vcf,
        is_pred=True,
        contig_to_asm=contig_to_asm,
        call_context=call_context,
        asm_to_scenario=asm_to_scenario,
        window_bp=window_bp,
    )

    bytype = sorted({row["type"] for row in truth} | {row["type"] for row in pred})
    hint_leaked_preds = [row for row in pred if row["hint_driven"]]
    algo_pred = [row for row in pred if not row["hint_driven"]]

    stats: dict[str, list[int]] = defaultdict(lambda: [0, 0, 0])
    tp_records: list[tuple[dict[str, Any], dict[str, Any]]] = []
    fn_records: list[dict[str, Any]] = []
    used_global: set[int] = set()
    for truth_row in truth:
        best_idx: int | None = None
        best_dist = math.inf
        for idx, pred_row in enumerate(algo_pred):
            if idx in used_global or not compatible(truth_row, pred_row):
                continue
            dist = distance(truth_row, pred_row)
            if dist < best_dist:
                best_idx = idx
                best_dist = dist
        if best_idx is None:
            fn_records.append(truth_row)
            stats[truth_row["type"]][2] += 1
        else:
            used_global.add(best_idx)
            tp_records.append((truth_row, algo_pred[best_idx]))
            stats[truth_row["type"]][0] += 1

    fp_records = [algo_pred[idx] for idx in range(len(algo_pred)) if idx not in used_global]
    for pred_row in fp_records:
        stats[pred_row["type"]][1] += 1

    boundary_truth = [row for row in truth if row["near_boundary"]]
    interior_truth = [row for row in truth if not row["near_boundary"]]
    _, fn_b = match_truth_to_pred(boundary_truth, algo_pred)
    _, fn_i = match_truth_to_pred(interior_truth, algo_pred)
    tp_b = len(boundary_truth) - len(fn_b)
    tp_i = len(interior_truth) - len(fn_i)

    scen_stats: dict[str, dict[str, list[int]]] = defaultdict(lambda: defaultdict(lambda: [0, 0, 0]))
    used_scen: set[int] = set()
    for truth_row in truth:
        best_idx: int | None = None
        best_dist = math.inf
        for idx, pred_row in enumerate(algo_pred):
            if idx in used_scen or not compatible(truth_row, pred_row):
                continue
            dist = distance(truth_row, pred_row)
            if dist < best_dist:
                best_idx = idx
                best_dist = dist
        scenario = truth_row.get("scenario", "unknown")
        if best_idx is None:
            scen_stats[scenario][truth_row["type"]][2] += 1
        else:
            used_scen.add(best_idx)
            scen_stats[scenario][truth_row["type"]][0] += 1
    for idx, pred_row in enumerate(algo_pred):
        if idx not in used_scen:
            scen_stats[pred_row.get("scenario", "unknown")][pred_row["type"]][1] += 1

    total_tp = sum(v[0] for v in stats.values())
    total_fp = sum(v[1] for v in stats.values())
    total_fn = sum(v[2] for v in stats.values())
    prec_o = total_tp / (total_tp + total_fp) if (total_tp + total_fp) else 0.0
    rec_o = total_tp / (total_tp + total_fn) if (total_tp + total_fn) else 0.0
    f1_o = 2 * prec_o * rec_o / (prec_o + rec_o) if (prec_o + rec_o) else 0.0
    prec_lo, prec_hi = wilson_ci(total_tp, total_tp + total_fp)
    rec_lo, rec_hi = wilson_ci(total_tp, total_tp + total_fn)

    summary: dict[str, Any] = {
        "truth_records": len(truth),
        "pred_records_total": len(pred),
        "pred_records_algo": len(algo_pred),
        "hint_leaked_calls": len(hint_leaked_preds),
        "window_bp": window_bp,
        "overall": {
            "tp": total_tp,
            "fp": total_fp,
            "fn": total_fn,
            "precision": prec_o,
            "precision_ci95": [prec_lo, prec_hi],
            "recall": rec_o,
            "recall_ci95": [rec_lo, rec_hi],
            "f1": f1_o,
            "note": "computed on algorithmic predictions only; hint-leaked calls excluded",
        },
        "window_boundary": {
            "boundary_truth": len(boundary_truth),
            "boundary_tp": tp_b,
            "boundary_recall": tp_b / len(boundary_truth) if boundary_truth else 0.0,
            "interior_truth": len(interior_truth),
            "interior_tp": tp_i,
            "interior_recall": tp_i / len(interior_truth) if interior_truth else 0.0,
        },
        "hint_leak_diagnostic": {
            "hint_leaked_predictions": len(hint_leaked_preds),
            "note": "non-zero means hint path is still active — results are compromised",
        },
        "by_svtype": {},
        "by_scenario": {},
        "by_annotation": defaultdict(int),
    }

    for pred_row in algo_pred:
        annot = pred_row["annot"]
        if annot not in ("NONE", ".", ""):
            summary["by_annotation"][annot] += 1

    per_type_rows: list[dict[str, Any]] = []
    overall_row = {
        "svtype": "OVERALL",
        "tp": total_tp,
        "fp": total_fp,
        "fn": total_fn,
        "precision": prec_o,
        "prec_lo95": prec_lo,
        "prec_hi95": prec_hi,
        "recall": rec_o,
        "rec_lo95": rec_lo,
        "rec_hi95": rec_hi,
        "f1": f1_o,
    }
    per_type_rows.append(overall_row)

    for svtype in bytype:
        tp, fp, fn = stats[svtype]
        prec = tp / (tp + fp) if (tp + fp) else 0.0
        rec = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
        pl, ph = wilson_ci(tp, tp + fp)
        rl, rh = wilson_ci(tp, tp + fn)
        row = {
            "svtype": svtype,
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "precision": prec,
            "prec_lo95": pl,
            "prec_hi95": ph,
            "recall": rec,
            "rec_lo95": rl,
            "rec_hi95": rh,
            "f1": f1,
        }
        per_type_rows.append(row)
        summary["by_svtype"][svtype] = {
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "precision": prec,
            "precision_ci95": [pl, ph],
            "recall": rec,
            "recall_ci95": [rl, rh],
            "f1": f1,
        }

    summary["per_type_rows"] = per_type_rows
    summary["boundary_rows"] = [
        {
            "region": "boundary",
            "truth_count": len(boundary_truth),
            "tp": tp_b,
            "recall": tp_b / len(boundary_truth) if boundary_truth else 0.0,
        },
        {
            "region": "interior",
            "truth_count": len(interior_truth),
            "tp": tp_i,
            "recall": tp_i / len(interior_truth) if interior_truth else 0.0,
        },
    ]
    summary["hint_rows"] = [
        {
            "metric": "total_predictions",
            "value": len(pred),
            "note": "all predictions received",
        },
        {
            "metric": "algorithmic_predictions",
            "value": len(algo_pred),
            "note": "predictions from plain contig names",
        },
        {
            "metric": "hint_leaked_predictions",
            "value": len(hint_leaked_preds),
            "note": "predictions from __sv_-encoded contig names (must be 0)",
        },
        {
            "metric": "hint_leak_status",
            "value": "FAIL" if hint_leaked_preds else "PASS",
            "note": "FAIL means hint path is active and P/R figures are invalid",
        },
    ]

    scenario_rows: list[dict[str, Any]] = []
    for scenario in sorted(scen_stats):
        summary["by_scenario"][scenario] = {}
        for svtype in sorted(scen_stats[scenario]):
            tp, fp, fn = scen_stats[scenario][svtype]
            prec = tp / (tp + fp) if (tp + fp) else 0.0
            rec = tp / (tp + fn) if (tp + fn) else 0.0
            f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
            scenario_rows.append({
                "scenario": scenario,
                "svtype": svtype,
                "tp": tp,
                "fp": fp,
                "fn": fn,
                "precision": prec,
                "recall": rec,
                "f1": f1,
            })
            summary["by_scenario"][scenario][svtype] = {"tp": tp, "fp": fp, "fn": fn}
    summary["scenario_rows"] = scenario_rows
    summary["by_annotation"] = dict(summary["by_annotation"])
    summary["tp_records"] = tp_records
    summary["fp_records"] = fp_records
    summary["fn_records"] = fn_records
    return summary


def expand_to_multisample_vcf(src_vcf: Path, dst_vcf: Path) -> Path:
    """Rewrite a single-sample MycoSV VCF as a multi-sample one.

    The MycoSV binary emits one SAMPLE column with each row's per-query
    provenance buried in the QASM (pred) or QUERY_ASM (truth) INFO field.
    For multi-query benchmarks the user-facing expectation is one column
    per query asm with GT 1/1 only for the owning sample. This walks the
    file twice (cheap — VCFs here are small): pass 1 to enumerate sample
    names, pass 2 to rewrite rows. Returns the destination path.
    """
    samples: list[str] = []
    seen: set[str] = set()
    with src_vcf.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            info = parse_info_field(cols[7])
            asm = info.get("QASM") or info.get("QUERY_ASM") or ""
            if asm and asm not in seen:
                seen.add(asm)
                samples.append(asm)
    if not samples:
        # Nothing to expand — copy through as-is.
        dst_vcf.write_text(src_vcf.read_text(encoding="utf-8"), encoding="utf-8")
        return dst_vcf
    sample_index = {asm: i for i, asm in enumerate(samples)}
    dst_vcf.parent.mkdir(parents=True, exist_ok=True)
    with src_vcf.open(encoding="utf-8") as fh, dst_vcf.open("w", encoding="utf-8") as out:
        for line in fh:
            if line.startswith("##"):
                out.write(line)
                continue
            if line.startswith("#CHROM"):
                fixed = "\t".join([
                    "#CHROM", "POS", "ID", "REF", "ALT",
                    "QUAL", "FILTER", "INFO", "FORMAT",
                ] + samples)
                out.write(fixed + "\n")
                continue
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                out.write(line)
                continue
            info = parse_info_field(cols[7])
            asm = info.get("QASM") or info.get("QUERY_ASM") or ""
            owner = sample_index.get(asm, -1)
            # Preserve FORMAT but force GT to be the only key so the
            # multi-sample expansion is unambiguous; this matches the
            # truth VCF schema produced by test_amf.write_truth_vcf.
            cols[8] = "GT"
            gts = ["0/0"] * len(samples)
            if 0 <= owner < len(gts):
                gts[owner] = "1/1"
            out.write("\t".join(cols[:9] + gts) + "\n")
    return dst_vcf


def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_pr_artifacts(summary: dict[str, Any], out_tsv: Path, out_json: Path) -> None:
    _write_tsv(
        out_tsv,
        ["svtype", "tp", "fp", "fn", "precision", "prec_lo95", "prec_hi95", "recall", "rec_lo95", "rec_hi95", "f1"],
        summary["per_type_rows"],
    )
    _write_tsv(
        out_tsv.with_name(out_tsv.stem + "_boundary.tsv"),
        ["region", "truth_count", "tp", "recall"],
        summary["boundary_rows"],
    )
    _write_tsv(
        out_tsv.with_name(out_tsv.stem + "_hint_leak_diagnostic.tsv"),
        ["metric", "value", "note"],
        summary["hint_rows"],
    )
    _write_tsv(
        out_tsv.with_name(out_tsv.stem + "_by_scenario.tsv"),
        ["scenario", "svtype", "tp", "fp", "fn", "precision", "recall", "f1"],
        summary["scenario_rows"],
    )
    out_json.parent.mkdir(parents=True, exist_ok=True)
    with out_json.open("w") as fh:
        json.dump({k: v for k, v in summary.items() if k not in {"tp_records", "fp_records", "fn_records"}}, fh, indent=2, sort_keys=True)
