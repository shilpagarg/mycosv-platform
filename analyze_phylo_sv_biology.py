#!/usr/bin/env python3
"""
Biological analysis of fungal SV calls across phylogeny.

Addresses three questions:
  1. How SVs distribute across phylogeny (by phylum/class/order)
  2. How MGEs shape genome architecture globally
  3. How HGT propagates across clades (TRA/OFF_REF as HGT proxies)

Inputs come from:
  --vcf-dirs  : one or more directories containing calls.vcf (real benchmark output)
  --biology-dirs : one or more biology/ subdirs with biology_candidates.tsv
  --taxonomy  : data_cache/taxonomy_cache.json (built during prepare)
  --phenotype : data_cache/phenotypic_metadata.json (built during prepare)
  --query-manifest : query_manifest.tsv (from prepare --out-dir)
  --out-dir   : where to write analysis outputs

Outputs (all in --out-dir):
  phylo_sv_distribution.tsv   -- SV counts by phylum/class/order/svtype
  mge_architecture.tsv        -- MGE type × SV type co-occurrence per clade
  hgt_propagation.tsv         -- TRA/OFF_REF events grouped by clade pair
  summary.json                -- machine-readable aggregation of all three
  plots/                      -- PNG figures for each question
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# VCF parsing
# ---------------------------------------------------------------------------

SV_TYPES = {"DEL", "INS", "DUP", "INV", "TRA", "OFF_REF"}

MGE_CLASSES = {
    "TE", "TE_LTR", "TE_TIR", "TE_LINE", "TE_SINE",
    "LTR_GYPSY", "LTR_COPIA", "HELITRON", "MITE",
    "STARSHIP", "HGT", "REPEAT", "LINE", "SINE", "DNA_TIR", "RIP",
}

MGE_INTEGRATIVE = {"STARSHIP", "HGT"}
MGE_TRANSPOSABLE = {"TE", "TE_LTR", "TE_TIR", "TE_LINE", "TE_SINE",
                    "LTR_GYPSY", "LTR_COPIA", "HELITRON", "MITE",
                    "LINE", "SINE", "DNA_TIR"}
MGE_REPEAT_BASED = {"REPEAT", "RIP"}

HGT_SV_TYPES = {"TRA", "OFF_REF"}


def _info_val(info_str: str, key: str) -> str:
    for field in info_str.split(";"):
        if field.startswith(key + "="):
            return field[len(key) + 1:]
        if field == key:
            return "1"
    return ""


def _first_info_val(info_str: str, *keys: str) -> str:
    for key in keys:
        val = _info_val(info_str, key)
        if val:
            return val
    return ""


def parse_vcf(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    try:
        with path.open() as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    continue
                info = parts[7]
                svtype = _info_val(info, "SVTYPE") or parts[4].lstrip("<").rstrip(">")
                if svtype not in SV_TYPES:
                    continue
                rows.append({
                    "chrom": parts[0],
                    "pos": parts[1],
                    "svtype": svtype,
                    "svlen": _info_val(info, "SVLEN"),
                    "query_asm": _first_info_val(info, "QUERY_ASM", "QUERY", "SAMPLE"),
                    "element_class": _first_info_val(info, "ELEMENT_CLASS", "EC", "TE_CLASS"),
                    "chr2": _info_val(info, "CHR2"),
                    "end2": _info_val(info, "END2"),
                    "hgt_flag": _info_val(info, "HGT"),
                    "score": _info_val(info, "SCORE"),
                })
    except (OSError, UnicodeDecodeError):
        pass
    return rows


def parse_biology_tsv(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    try:
        with path.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                rows.append(dict(row))
    except (OSError, UnicodeDecodeError):
        pass
    return rows


def parse_query_manifest(path: Path) -> dict[str, dict[str, str]]:
    meta: dict[str, dict[str, str]] = {}
    try:
        with path.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                asm = row.get("query_asm", "").strip()
                if asm:
                    meta[asm] = dict(row)
    except (OSError, UnicodeDecodeError):
        pass
    return meta


# ---------------------------------------------------------------------------
# Analysis 1: SV distribution across phylogeny
# ---------------------------------------------------------------------------

def analyze_phylo_sv_distribution(
    sv_rows: list[dict[str, str]],
    asm_meta: dict[str, dict[str, str]],
    taxonomy: dict[str, dict[str, str]],
) -> list[dict[str, Any]]:
    """SV counts aggregated by taxonomic rank."""
    # Build taxid→lineage lookup keyed by species name too
    species_to_lineage: dict[str, dict[str, str]] = {}
    for lineage in taxonomy.values():
        sp = lineage.get("species", "")
        if sp and sp != ".":
            species_to_lineage[sp] = lineage

    counts: dict[tuple[str, str, str, str], int] = defaultdict(int)

    for sv in sv_rows:
        asm = sv.get("query_asm", "")
        svtype = sv.get("svtype", ".")
        meta = asm_meta.get(asm, {})
        species = meta.get("species", "")
        lineage = species_to_lineage.get(species, {})
        phylum = lineage.get("phylum", meta.get("phylum", ".")) or "."
        class_ = lineage.get("class", meta.get("class", ".")) or "."
        order = lineage.get("order", meta.get("order", ".")) or "."
        counts[(phylum, class_, order, svtype)] += 1

    rows = []
    for (phylum, class_, order, svtype), count in sorted(counts.items()):
        rows.append({
            "phylum": phylum,
            "class": class_,
            "order": order,
            "svtype": svtype,
            "count": count,
        })
    return rows


# ---------------------------------------------------------------------------
# Analysis 2: MGE shaping of genome architecture
# ---------------------------------------------------------------------------

def analyze_mge_architecture(
    sv_rows: list[dict[str, str]],
    bio_rows: list[dict[str, str]],
    asm_meta: dict[str, dict[str, str]],
    taxonomy: dict[str, dict[str, str]],
) -> list[dict[str, Any]]:
    """MGE type × SV type co-occurrence; genome-scale architecture effects."""
    species_to_lineage: dict[str, dict[str, str]] = {}
    for lineage in taxonomy.values():
        sp = lineage.get("species", "")
        if sp and sp != ".":
            species_to_lineage[sp] = lineage

    # Use biology candidates as primary source (richer element_class data)
    counts: dict[tuple[str, str, str, str], int] = defaultdict(int)

    for row in bio_rows:
        asm = row.get("query_asm", "") or row.get("sample", "")
        ec = row.get("element_class", "") or row.get("te_class", "") or row.get("effect", "") or "OTHER"
        svtype = row.get("svtype", "") or row.get("sv_type", "") or row.get("type", ".") or "."
        meta = asm_meta.get(asm, {})
        species = meta.get("species", "") or row.get("phylum", "")
        lineage = species_to_lineage.get(species, {})
        phylum = lineage.get("phylum", row.get("phylum", ".")) or "."

        mge_cat = (
            "integrative" if ec in MGE_INTEGRATIVE
            else "transposable" if ec in MGE_TRANSPOSABLE
            else "repeat" if ec in MGE_REPEAT_BASED
            else "none"
        )
        if mge_cat == "none":
            continue
        counts[(phylum, mge_cat, ec, svtype)] += 1

    # Also fold in raw SV rows that have element_class set
    for sv in sv_rows:
        ec = sv.get("element_class", "") or ""
        if not ec or ec == ".":
            continue
        asm = sv.get("query_asm", "")
        svtype = sv.get("svtype", ".") or "."
        meta = asm_meta.get(asm, {})
        species = meta.get("species", "")
        lineage = species_to_lineage.get(species, {})
        phylum = lineage.get("phylum", meta.get("phylum", ".")) or "."
        mge_cat = (
            "integrative" if ec in MGE_INTEGRATIVE
            else "transposable" if ec in MGE_TRANSPOSABLE
            else "repeat" if ec in MGE_REPEAT_BASED
            else "none"
        )
        if mge_cat != "none":
            counts[(phylum, mge_cat, ec, svtype)] += 1

    rows = []
    for (phylum, mge_cat, ec, svtype), count in sorted(counts.items()):
        rows.append({
            "phylum": phylum,
            "mge_category": mge_cat,
            "element_class": ec,
            "svtype": svtype,
            "count": count,
        })
    return rows


# ---------------------------------------------------------------------------
# Analysis 3: HGT propagation across clades
# ---------------------------------------------------------------------------

def analyze_hgt_propagation(
    sv_rows: list[dict[str, str]],
    bio_rows: list[dict[str, str]],
    asm_meta: dict[str, dict[str, str]],
    taxonomy: dict[str, dict[str, str]],
) -> list[dict[str, Any]]:
    """TRA and OFF_REF events as candidate proxies for HGT between clades.

    A TRA call where chr (source contig) and chr2 (mate contig) are in different
    species/orders is treated as screening evidence for possible inter-clade HGT
    propagation, not as proof without phylogenetic/contamination follow-up.
    """
    species_to_lineage: dict[str, dict[str, str]] = {}
    for lineage in taxonomy.values():
        sp = lineage.get("species", "")
        if sp and sp != ".":
            species_to_lineage[sp] = lineage

    # Use bio_rows for HGT candidates first (already scored and classified)
    hgt_rows = [
        r for r in bio_rows
        if r.get("hgt_flag") == "1"
        or r.get("candidate_type") == "hgt_candidate"
        or (r.get("element_class") or r.get("te_class") or r.get("effect")) in MGE_INTEGRATIVE
    ]
    # Supplement with raw TRA/OFF_REF from VCF
    vcf_hgt = [sv for sv in sv_rows if sv.get("svtype") in HGT_SV_TYPES]

    counts: dict[tuple[str, str, str], int] = defaultdict(int)
    event_rows: list[dict[str, Any]] = []

    def get_lineage_for_asm(asm: str) -> dict[str, str]:
        meta = asm_meta.get(asm, {})
        sp = meta.get("species", "")
        return species_to_lineage.get(sp, {
            "phylum": meta.get("phylum", "."),
            "class": meta.get("class", "."),
            "order": meta.get("order", "."),
            "species": sp,
        })

    for row in hgt_rows:
        asm = row.get("query_asm", "") or row.get("sample", "")
        lin = get_lineage_for_asm(asm)
        svtype = row.get("svtype", "") or row.get("sv_type", "") or row.get("type", "TRA")
        ec = row.get("element_class", "") or row.get("te_class", "") or row.get("effect", "") or "HGT"
        phylum = lin.get("phylum", ".")
        order = lin.get("order", ".")
        counts[(phylum, order, svtype)] += 1
        event_rows.append({
            "source_asm": asm,
            "source_phylum": phylum,
            "source_order": order,
            "svtype": svtype,
            "element_class": ec,
            "pos": row.get("pos", "."),
            "evidence": "biology_candidate",
            "hgt_flag": row.get("hgt_flag", "0"),
            "candidate_type": row.get("candidate_type", "."),
            "rationale": row.get("rationale", "."),
        })

    for sv in vcf_hgt:
        asm = sv.get("query_asm", "")
        lin = get_lineage_for_asm(asm)
        svtype = sv.get("svtype", "TRA")
        phylum = lin.get("phylum", ".")
        order = lin.get("order", ".")
        counts[(phylum, order, svtype)] += 1
        event_rows.append({
            "source_asm": asm,
            "source_phylum": phylum,
            "source_order": order,
            "svtype": svtype,
            "element_class": sv.get("element_class", "."),
            "pos": sv.get("pos", "."),
            "evidence": "vcf_call",
            "hgt_flag": sv.get("hgt_flag", "0"),
            "candidate_type": ".",
            "rationale": ".",
        })

    summary_rows: list[dict[str, Any]] = []
    for (phylum, order, svtype), count in sorted(counts.items(), key=lambda x: -x[1]):
        summary_rows.append({
            "phylum": phylum,
            "order": order,
            "svtype": svtype,
            "event_count": count,
        })

    return summary_rows, event_rows


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def make_plots(
    phylo_rows: list[dict[str, Any]],
    mge_rows: list[dict[str, Any]],
    hgt_summary: list[dict[str, Any]],
    out_dir: Path,
) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("[warn] matplotlib not available; skipping plots")
        return

    plots_dir = out_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    # --- Plot 1: SV distribution across phyla ---
    phyla = sorted({r["phylum"] for r in phylo_rows if r["phylum"] != "."})
    sv_types = ["DEL", "INS", "DUP", "INV", "TRA", "OFF_REF"]
    matrix: dict[str, dict[str, int]] = {ph: {sv: 0 for sv in sv_types} for ph in phyla}
    for r in phylo_rows:
        ph = r["phylum"]
        sv = r["svtype"]
        if ph in matrix and sv in sv_types:
            matrix[ph][sv] += r["count"]

    if phyla and any(sum(matrix[ph].values()) > 0 for ph in phyla):
        x = np.arange(len(phyla))
        width = 0.12
        fig, ax = plt.subplots(figsize=(max(8, len(phyla) * 1.5), 5))
        colors = ["#d62728", "#1f77b4", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b"]
        for i, (sv, color) in enumerate(zip(sv_types, colors)):
            vals = [matrix[ph][sv] for ph in phyla]
            ax.bar(x + i * width, vals, width, label=sv, color=color, alpha=0.85)
        ax.set_xlabel("Phylum")
        ax.set_ylabel("SV count")
        ax.set_title("SV Distribution Across Fungal Phylogeny")
        ax.set_xticks(x + width * (len(sv_types) - 1) / 2)
        ax.set_xticklabels(phyla, rotation=30, ha="right", fontsize=8)
        ax.legend(title="SV type", bbox_to_anchor=(1, 1))
        fig.tight_layout()
        fig.savefig(plots_dir / "fig1_phylo_sv_distribution.png", dpi=150)
        fig.savefig(plots_dir / "fig1_phylo_sv_distribution.svg")
        plt.close(fig)

    # --- Plot 2: MGE architecture — stacked bar by phylum and MGE category ---
    mge_phyla = sorted({r["phylum"] for r in mge_rows if r["phylum"] != "."})
    mge_cats = ["integrative", "transposable", "repeat"]
    mge_matrix: dict[str, dict[str, int]] = {ph: {c: 0 for c in mge_cats} for ph in mge_phyla}
    for r in mge_rows:
        ph = r["phylum"]
        cat = r["mge_category"]
        if ph in mge_matrix and cat in mge_cats:
            mge_matrix[ph][cat] += r["count"]

    if mge_phyla and any(sum(mge_matrix[ph].values()) > 0 for ph in mge_phyla):
        fig, ax = plt.subplots(figsize=(max(7, len(mge_phyla) * 1.2), 5))
        x = np.arange(len(mge_phyla))
        bottoms = np.zeros(len(mge_phyla))
        cat_colors = {"integrative": "#e41a1c", "transposable": "#377eb8", "repeat": "#4daf4a"}
        for cat in mge_cats:
            vals = np.array([mge_matrix[ph][cat] for ph in mge_phyla], dtype=float)
            ax.bar(x, vals, bottom=bottoms, label=cat, color=cat_colors[cat], alpha=0.85)
            bottoms += vals
        ax.set_xlabel("Phylum")
        ax.set_ylabel("MGE-linked SV count")
        ax.set_title("MGE Impact on Genome Architecture by Phylum")
        ax.set_xticks(x)
        ax.set_xticklabels(mge_phyla, rotation=30, ha="right", fontsize=8)
        ax.legend(title="MGE category")
        fig.tight_layout()
        fig.savefig(plots_dir / "fig2_mge_architecture.png", dpi=150)
        fig.savefig(plots_dir / "fig2_mge_architecture.svg")
        plt.close(fig)

    # --- Plot 3: HGT propagation by order ---
    hgt_orders = sorted({r["order"] for r in hgt_summary if r["order"] != "."}, key=lambda o: -sum(r["event_count"] for r in hgt_summary if r["order"] == o))[:15]
    if hgt_orders:
        fig, ax = plt.subplots(figsize=(10, 5))
        hgt_sv_types = ["TRA", "OFF_REF"]
        x = np.arange(len(hgt_orders))
        width = 0.35
        hgt_mat: dict[str, dict[str, int]] = {o: {sv: 0 for sv in hgt_sv_types} for o in hgt_orders}
        for r in hgt_summary:
            o = r["order"]
            sv = r["svtype"]
            if o in hgt_mat and sv in hgt_sv_types:
                hgt_mat[o][sv] += r["event_count"]
        tra_vals = [hgt_mat[o]["TRA"] for o in hgt_orders]
        offref_vals = [hgt_mat[o]["OFF_REF"] for o in hgt_orders]
        ax.bar(x - width / 2, tra_vals, width, label="TRA", color="#e41a1c", alpha=0.85)
        ax.bar(x + width / 2, offref_vals, width, label="OFF_REF", color="#ff7f0e", alpha=0.85)
        ax.set_xlabel("Order")
        ax.set_ylabel("HGT-proxy event count")
        ax.set_title("HGT Propagation Across Fungal Orders (top 15)")
        ax.set_xticks(x)
        ax.set_xticklabels(hgt_orders, rotation=40, ha="right", fontsize=8)
        ax.legend()
        fig.tight_layout()
        fig.savefig(plots_dir / "fig3_hgt_propagation.png", dpi=150)
        fig.savefig(plots_dir / "fig3_hgt_propagation.svg")
        plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf-dirs", nargs="+", type=Path, default=[],
                    help="Directories containing calls.vcf from real benchmark runs")
    ap.add_argument("--biology-dirs", nargs="+", type=Path, default=[],
                    help="Directories containing biology_candidates.tsv")
    ap.add_argument("--query-manifests", nargs="+", type=Path, default=[],
                    help="query_manifest.tsv files from prepare runs")
    ap.add_argument("--taxonomy", type=Path,
                    default=Path("data_cache/taxonomy_cache.json"),
                    help="Taxonomy lineage cache JSON (from prepare)")
    ap.add_argument("--phenotype", type=Path,
                    default=Path("data_cache/phenotypic_metadata.json"),
                    help="BioSample phenotype cache JSON (from prepare)")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    # Load taxonomy
    taxonomy: dict[str, dict[str, str]] = {}
    if args.taxonomy.exists():
        try:
            taxonomy = json.loads(args.taxonomy.read_text())
            print(f"[taxonomy] loaded {len(taxonomy)} lineages from {args.taxonomy}")
        except Exception as e:
            print(f"[warn] could not load taxonomy: {e}")

    # Load phenotype metadata
    phenotype: dict[str, dict[str, str]] = {}
    if args.phenotype.exists():
        try:
            phenotype = json.loads(args.phenotype.read_text())
            print(f"[phenotype] loaded {len(phenotype)} BioSample records from {args.phenotype}")
        except Exception as e:
            print(f"[warn] could not load phenotype: {e}")

    # Build assembly metadata from query manifests
    asm_meta: dict[str, dict[str, str]] = {}
    for qm in args.query_manifests:
        asm_meta.update(parse_query_manifest(qm))
    print(f"[manifest] {len(asm_meta)} assemblies from {len(args.query_manifests)} manifest(s)")

    # Collect SV calls
    all_sv_rows: list[dict[str, str]] = []
    for vcf_dir in args.vcf_dirs:
        vcf = vcf_dir / "calls.vcf"
        if not vcf.exists():
            vcf = next(vcf_dir.glob("*.vcf"), None)
        if vcf:
            rows = parse_vcf(vcf)
            print(f"[vcf] {len(rows)} SVs from {vcf}")
            all_sv_rows.extend(rows)

    # Collect biology candidates
    all_bio_rows: list[dict[str, str]] = []
    for bio_dir in args.biology_dirs:
        tsv = bio_dir / "biology_candidates.tsv"
        if tsv.exists():
            rows = parse_biology_tsv(tsv)
            print(f"[bio] {len(rows)} biology candidates from {tsv}")
            all_bio_rows.extend(rows)

    print(f"[total] {len(all_sv_rows)} SV calls, {len(all_bio_rows)} biology candidates")

    # Run the three analyses
    phylo_rows = analyze_phylo_sv_distribution(all_sv_rows, asm_meta, taxonomy)
    mge_rows = analyze_mge_architecture(all_sv_rows, all_bio_rows, asm_meta, taxonomy)
    hgt_summary, hgt_events = analyze_hgt_propagation(all_sv_rows, all_bio_rows, asm_meta, taxonomy)

    # Write TSVs
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out_dir / "phylo_sv_distribution.tsv", phylo_rows,
              ["phylum", "class", "order", "svtype", "count"])
    write_tsv(args.out_dir / "mge_architecture.tsv", mge_rows,
              ["phylum", "mge_category", "element_class", "svtype", "count"])
    write_tsv(args.out_dir / "hgt_propagation.tsv", hgt_summary,
              ["phylum", "order", "svtype", "event_count"])
    write_tsv(args.out_dir / "hgt_events.tsv", hgt_events,
              ["source_asm", "source_phylum", "source_order", "svtype",
               "element_class", "pos", "evidence", "hgt_flag", "candidate_type", "rationale"])

    # Write summary JSON
    total_by_svtype: dict[str, int] = {}
    for r in phylo_rows:
        total_by_svtype[r["svtype"]] = total_by_svtype.get(r["svtype"], 0) + r["count"]

    mge_totals: dict[str, int] = {}
    for r in mge_rows:
        mge_totals[r["mge_category"]] = mge_totals.get(r["mge_category"], 0) + r["count"]

    summary = {
        "total_sv_calls": len(all_sv_rows),
        "total_biology_candidates": len(all_bio_rows),
        "assemblies_with_metadata": len(asm_meta),
        "q1_sv_across_phylogeny": {
            "phyla_covered": len({r["phylum"] for r in phylo_rows}),
            "by_svtype": total_by_svtype,
            "rows": len(phylo_rows),
        },
        "q2_mge_architecture": {
            "mge_linked_events": sum(mge_totals.values()),
            "by_mge_category": mge_totals,
            "rows": len(mge_rows),
        },
        "q3_hgt_propagation": {
            "hgt_proxy_events": sum(r["event_count"] for r in hgt_summary),
            "clades_implicated": len({r["order"] for r in hgt_summary if r["order"] != "."}),
            "by_svtype": {sv: sum(r["event_count"] for r in hgt_summary if r["svtype"] == sv) for sv in HGT_SV_TYPES},
        },
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))

    print(f"\n=== Biological Analysis Summary ===")
    print(f"Q1 Phylo SV distribution: {summary['q1_sv_across_phylogeny']['phyla_covered']} phyla, SVs: {total_by_svtype}")
    print(f"Q2 MGE architecture:      {summary['q2_mge_architecture']['mge_linked_events']} MGE-linked events, categories: {mge_totals}")
    print(f"Q3 HGT propagation:       {summary['q3_hgt_propagation']['hgt_proxy_events']} events across {summary['q3_hgt_propagation']['clades_implicated']} orders")

    make_plots(phylo_rows, mge_rows, hgt_summary, args.out_dir)

    print(f"\nOutputs written to {args.out_dir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
