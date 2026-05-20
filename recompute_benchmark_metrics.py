#!/usr/bin/env python3
"""Recompute exact_benchmark_summary.tsv from existing benchmark outputs.

Reads each panel's prepared/query_manifest.tsv + mycosv/calls.vcf +
comparators/{tool}/{query}/* and re-derives the TSV with the new "Fix B"
(any-clade row) and comparator-agreement row logic, plus
OFF_REF/misrouted-call diagnostics in the JSON. Tools are NOT re-run.

Usage:
    python3 recompute_benchmark_metrics.py PANEL_DIR [PANEL_DIR ...]

Each PANEL_DIR is expected to contain prepared/ and benchmark_*/ subdirs.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

from run_real_fungal_benchmark import (
    DEFAULT_DATA_CACHE,
    NormalizedCall,
    build_consensus_truth,
    fasta_contig_names,
    load_mycosv_query_calls,
    load_mycosv_reference_calls,
    load_minigraph_bubble_calls,
    load_reference_vcf_calls,
    load_syri_query_calls,
    score_callsets,
    validation_basis_for_label,
    _emit_per_svtype_rows,
    _lift_calls_to_benchmark_ref,
    tool_path,
    write_agreement_summary,
)


def load_query_manifest(prepared_dir: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    path = prepared_dir / "query_manifest.tsv"
    if not path.exists():
        return rows
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows.append(row)
    return rows


def discover_truth_sets(
    out_dir: Path, mode: str, query_manifest: list[dict[str, str]]
) -> dict[str, dict[tuple[str, str], list[NormalizedCall]]]:
    """Walk benchmark_*/comparators/{tool}/{query}/ and load comparator callsets."""
    truth: dict[str, dict[tuple[str, str], list[NormalizedCall]]] = defaultdict(dict)
    comp_root = out_dir / "comparators"
    if not comp_root.exists():
        return truth

    for query_row in query_manifest:
        if query_row.get("query_mode") and query_row["query_mode"] != mode:
            continue
        qasm = query_row["query_asm"]

        # Assembly-mode comparators
        if mode == "assembly":
            mg_dir = comp_root / "minigraph" / qasm
            if mg_dir.exists():
                bubble_bed = mg_dir / "bubbles.bed"
                sample_bed = mg_dir / "sample.bed"
                if bubble_bed.exists() and sample_bed.exists():
                    calls = load_minigraph_bubble_calls(bubble_bed, sample_bed, qasm)
                    if calls:
                        truth[qasm][("reference", "minigraph")] = calls

            syri_tsv = comp_root / "syri" / qasm / "normalized.tsv"
            if syri_tsv.exists():
                calls = load_syri_query_calls(syri_tsv, qasm)
                if calls:
                    truth[qasm][("query", "syri")] = calls

            for label, vcf_rel in (
                ("cactus", "pangenome.vcf.gz"),
                ("svim_asm", "svim_asm_out/variants.vcf"),
                ("anchorwave", "anchorwave.vcf"),
                ("pggb", "pggb.smoothxg.fix.vcf"),
            ):
                vcf_path = comp_root / label / qasm / vcf_rel
                if not vcf_path.exists() and label == "pggb":
                    # pggb output filename varies by version; pick first .vcf
                    pggb_dir = comp_root / "pggb" / qasm
                    if pggb_dir.exists():
                        vcfs = sorted(pggb_dir.glob("*.vcf*"))
                        vcf_path = vcfs[0] if vcfs else vcf_path
                if vcf_path.exists() and vcf_path.stat().st_size > 0:
                    calls = load_reference_vcf_calls(vcf_path, label, qasm)
                    if calls:
                        truth[qasm][("reference", label)] = calls

        elif mode == "long-reads":
            for label in ("svim", "sniffles", "cutesv"):
                d = comp_root / label / qasm
                if not d.exists():
                    continue
                # svim writes to <d>/svim_out/variants.vcf, not <d>/*.vcf —
                # the previous glob missed it and silently dropped svim from
                # the comparator callsets for every long-reads recompute.
                candidates: list[Path] = []
                if label == "svim":
                    candidates.append(d / "svim_out" / "variants.vcf")
                candidates.extend(d.glob("*.vcf"))
                candidates.extend(d.glob("*.vcf.gz"))
                vcfs = sorted({p for p in candidates if p.exists() and p.stat().st_size > 0})
                if vcfs:
                    calls = load_reference_vcf_calls(vcfs[0], label, qasm)
                    if calls:
                        truth[qasm][("reference", label)] = calls

        elif mode == "short-reads":
            for label in ("delly", "manta"):
                d = comp_root / label / qasm
                if not d.exists():
                    continue
                # Manta writes to <d>/results/variants/diploidSV.vcf.gz; delly
                # writes a single <d>/delly.vcf next to the BAM. The flat glob
                # alone misses Manta on real runs.
                candidates: list[Path] = []
                if label == "manta":
                    candidates.extend([
                        d / "results" / "variants" / "diploidSV.vcf.gz",
                        d / "results" / "variants" / "candidateSV.vcf.gz",
                    ])
                candidates.extend(d.glob("*.vcf"))
                candidates.extend(d.glob("*.vcf.gz"))
                vcfs = sorted({p for p in candidates if p.exists() and p.stat().st_size > 0})
                if vcfs:
                    calls = load_reference_vcf_calls(vcfs[0], label, qasm)
                    if calls:
                        truth[qasm][("reference", label)] = calls

    return truth


def recompute_for_mode(panel_dir: Path, mode: str) -> bool:
    out_dir = panel_dir / f"benchmark_{mode}"
    if not out_dir.exists():
        return False
    prepared_dir = panel_dir / "prepared"
    if not prepared_dir.exists():
        prepared_dir = panel_dir
    mycosv_vcf = out_dir / "mycosv" / "calls.vcf"
    if not mycosv_vcf.exists():
        sys.stderr.write(f"[skip] no mycosv VCF at {mycosv_vcf}\n")
        return False
    if mycosv_vcf.stat().st_size == 0:
        sys.stderr.write(
            f"[warn] mycosv VCF is empty for {panel_dir.name}/{mode}; "
            f"recompute will produce zero-pred rows\n"
        )

    query_manifest = load_query_manifest(prepared_dir)
    if not query_manifest:
        sys.stderr.write(f"[skip] no query manifest at {prepared_dir}\n")
        return False
    query_manifest = [q for q in query_manifest if q.get("query_mode", mode) == mode]

    mycosv_calls_by_query: dict[str, dict[str, list[NormalizedCall]]] = {}
    for q in query_manifest:
        qasm = q["query_asm"]
        mycosv_calls_by_query[qasm] = {
            "query": load_mycosv_query_calls(mycosv_vcf, qasm),
            "reference": load_mycosv_reference_calls(mycosv_vcf, qasm),
        }

    truth_sets = discover_truth_sets(out_dir, mode, query_manifest)

    agreement_rows: list[dict] = []
    summary_json: dict = {"mode": mode, "queries": {}}
    lift_cache_dir = out_dir / "lift_cache"
    for q in query_manifest:
        qasm = q["query_asm"]
        mq = mycosv_calls_by_query.get(qasm, {}).get("query", [])
        mr_all = mycosv_calls_by_query.get(qasm, {}).get("reference", [])
        bench_ref = q.get("benchmark_ref_fasta") or "."
        bench_contigs = (
            fasta_contig_names(Path(bench_ref))
            if bench_ref not in {"", "."}
            else frozenset()
        )
        # Lift mycosv sibling-clade refs onto benchmark_ref_fasta for reads
        # modes so per-comparator PR scoring is apples to apples. Mirrors the
        # online benchmark() path; the cached PAF is reused if present.
        if (
            mode in {"long-reads", "short-reads"}
            and bench_contigs
            and bench_ref not in {"", "."}
            and tool_path("minimap2") is not None
        ):
            mr_all = _lift_calls_to_benchmark_ref(
                mr_all,
                Path(bench_ref),
                DEFAULT_DATA_CACHE,
                lift_cache_dir / qasm,
                threads=4,
            )
        mr = (
            [c for c in mr_all if c.ref_contig in bench_contigs]
            if bench_contigs
            else mr_all
        )
        off_ref_dropped = sum(1 for c in mq if c.svtype == "OFF_REF")
        misrouted = max(0, len(mr_all) - len(mr))
        summary_json["queries"][qasm] = {
            "mycosv_calls": {
                "query": len(mq),
                "reference": len(mr),
                "reference_total": len(mr_all),
                "benchmark_ref_contigs": len(bench_contigs),
                "off_ref_dropped": off_ref_dropped,
                "misrouted_to_sibling_clade": misrouted,
            },
            "exact_benchmarks": {
                "query": {},
                "reference": {},
                "reference_any_clade": {},
            },
        }

        truth_for_query = truth_sets.get(qasm, {})

        for (coord_space, label), truth_calls in truth_for_query.items():
            pred_calls = mr if coord_space == "reference" else mq
            m = score_callsets(truth_calls, pred_calls)
            agreement_rows.extend(_emit_per_svtype_rows(
                query_asm=qasm,
                coord_space=coord_space,
                truth_label=label,
                method="mycosv",
                truth_calls=truth_calls,
                pred_calls=pred_calls,
            ))
            summary_json["queries"][qasm]["exact_benchmarks"][coord_space][label] = m

            supported_preds = [
                c for c in pred_calls
                if c.read_support is not None and c.read_support >= 1
            ]
            if supported_preds:
                agreement_rows.extend(_emit_per_svtype_rows(
                    query_asm=qasm,
                    coord_space=coord_space,
                    truth_label=label,
                    method="mycosv_read_supported",
                    truth_calls=truth_calls,
                    pred_calls=supported_preds,
                ))

        for (coord_space, label), truth_calls in truth_for_query.items():
            if coord_space != "reference":
                continue
            m_any = score_callsets(truth_calls, mr_all)
            agreement_rows.extend(_emit_per_svtype_rows(
                query_asm=qasm,
                coord_space="reference_any_clade",
                truth_label=label,
                method="mycosv",
                truth_calls=truth_calls,
                pred_calls=mr_all,
            ))
            summary_json["queries"][qasm]["exact_benchmarks"]["reference_any_clade"][
                label
            ] = m_any

        for coord_space in ("query", "reference"):
            ref_labels = [
                lbl
                for (cs, lbl) in truth_for_query.keys()
                if cs == coord_space and not lbl.startswith("consensus_")
            ]
            if len(ref_labels) < 2:
                continue
            consensus = build_consensus_truth(
                [truth_for_query[(coord_space, lbl)] for lbl in ref_labels],
                min_support=2,
            )
            pred_calls = mr if coord_space == "reference" else mq
            m_c = score_callsets(consensus, pred_calls)
            label = f"consensus_2of_{len(ref_labels)}"
            agreement_rows.extend(_emit_per_svtype_rows(
                query_asm=qasm,
                coord_space=coord_space,
                truth_label=label,
                method="mycosv",
                truth_calls=consensus,
                pred_calls=pred_calls,
            ))
            summary_json["queries"][qasm]["exact_benchmarks"][coord_space][label] = m_c

    if not agreement_rows:
        for q in query_manifest:
            qasm = q["query_asm"]
            for coord_space, preds in (
                ("query", mycosv_calls_by_query.get(qasm, {}).get("query", [])),
                ("reference", mycosv_calls_by_query.get(qasm, {}).get("reference", [])),
            ):
                agreement_rows.append({
                    "query_asm": qasm,
                    "coordinate_space": coord_space,
                    "truth_label": "no_comparator",
                    "validation_basis": validation_basis_for_label("no_comparator"),
                    "svtype": "ALL",
                    "method": "mycosv",
                    "truth_calls": float("nan"),
                    "pred_calls": len(preds),
                    "tp": float("nan"),
                    "fp": float("nan"),
                    "fn": float("nan"),
                    "precision": float("nan"),
                    "recall": float("nan"),
                    "f1": float("nan"),
                    "prec_lo95": float("nan"),
                    "prec_hi95": float("nan"),
                    "rec_lo95": float("nan"),
                    "rec_hi95": float("nan"),
                    "status": "no_truth",
                })

    tsv_path = out_dir / "exact_benchmark_summary.tsv"
    write_agreement_summary(tsv_path, agreement_rows)
    json_path = out_dir / "benchmark_summary.json"
    json_path.write_text(json.dumps(summary_json, indent=2))
    sys.stderr.write(
        f"[ok] {panel_dir.name}/{mode}: "
        f"{len(agreement_rows)} rows -> {tsv_path}\n"
    )
    return True


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("panel_dirs", nargs="+", type=Path)
    p.add_argument("--modes", default="assembly,short-reads,long-reads")
    args = p.parse_args()
    modes = [m.strip() for m in args.modes.split(",") if m.strip()]
    any_done = False
    for panel_dir in args.panel_dirs:
        for mode in modes:
            if recompute_for_mode(panel_dir, mode):
                any_done = True
    return 0 if any_done else 1


if __name__ == "__main__":
    sys.exit(main())
