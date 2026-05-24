#!/usr/bin/env python3
"""Plot MycoSV pangenome-call summaries for fungal biology reports.

The main visualization report is broad: benchmarking, biology, and novel-SV
questions all live together. This script makes a smaller, community-facing set
of figures focused on the central MycoSV claim: pangenome routing recovers SVs
that a single-reference benchmark cannot represent.
"""

from __future__ import annotations

import argparse
import csv
import html
import math
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    fallback_python = Path(
        os.environ.get(
            "MYCOSV_ENV_PYTHON",
            "/mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/tools/envs/envs/fungi_graph_sv/bin/python",
        )
    )
    if fallback_python.exists() and Path(sys.executable).resolve() != fallback_python.resolve():
        os.execv(str(fallback_python), [str(fallback_python), *sys.argv])
    plt = None


PALETTE = {
    "pangenome_only": "#0072B2",
    "single_reference_equivalent": "#D55E00",
    "read_supported": "#009E73",
    "intrinsic_only": "#CC79A7",
    "HGT": "#E69F00",
    "RIP": "#56B4E9",
    "TE": "#009E73",
    "NONE": "#7F7F7F",
}


def read_tsv(path: Path | None) -> list[dict[str, str]]:
    if path is None or not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def as_int(value: object, default: int = 0) -> int:
    try:
        if value in {"", ".", None}:
            return default
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def as_float(value: object, default: float = 0.0) -> float:
    try:
        if value in {"", ".", None}:
            return default
        return float(str(value))
    except (TypeError, ValueError):
        return default


def norm_token(value: str | None, fallback: str = "unknown") -> str:
    v = (value or "").strip()
    return v if v and v != "." else fallback


def call_scope(row: dict[str, str]) -> str:
    return (
        "single_reference_equivalent"
        if row.get("single_reference_equivalent", "").lower() == "yes"
        else "pangenome_only"
    )


def element_group(value: str | None) -> str:
    v = norm_token(value, "NONE").upper()
    if v.startswith("TE_") or v in {"REPEAT", "LTR", "LINE", "SINE", "TIR"}:
        return "TE / repeat"
    if v == "HGT":
        return "HGT / Starship"
    if v == "RIP":
        return "RIP / genome defense"
    if v == "NONE":
        return "Structural (other)"
    return v


def top_items(counter: Counter[str], n: int = 12) -> list[tuple[str, int]]:
    items = [(k, v) for k, v in counter.items() if k not in {"", "."} and v > 0]
    return sorted(items, key=lambda kv: (-kv[1], kv[0]))[:n]


def save_bar(
    path: Path,
    labels: list[str],
    values: list[int | float],
    *,
    title: str,
    ylabel: str = "Calls",
    colors: list[str] | None = None,
    rotate: int = 25,
) -> None:
    if plt is None or not labels:
        return
    fig_w = max(7.0, min(14.0, 0.55 * len(labels) + 4.5))
    fig, ax = plt.subplots(figsize=(fig_w, 4.6))
    ax.bar(labels, values, color=colors or "#4C78A8")
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=rotate)
    for tick in ax.get_xticklabels():
        tick.set_ha("right")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def save_stacked_bar(
    path: Path,
    table: dict[str, Counter[str]],
    *,
    title: str,
    ylabel: str = "Calls",
    normalize: bool = False,
) -> None:
    if plt is None or not table:
        return
    x_labels = sorted(table)
    stack_labels = sorted({k for c in table.values() for k in c})
    if not x_labels or not stack_labels:
        return
    bottoms = [0.0 for _ in x_labels]
    fig_w = max(7.0, min(14.0, 0.65 * len(x_labels) + 4.5))
    fig, ax = plt.subplots(figsize=(fig_w, 4.8))
    for stack in stack_labels:
        vals: list[float] = []
        for x in x_labels:
            total = sum(table[x].values())
            raw = table[x].get(stack, 0)
            vals.append((raw / total) if normalize and total else raw)
        ax.bar(x_labels, vals, bottom=bottoms, label=stack)
        bottoms = [b + v for b, v in zip(bottoms, vals)]
    ax.set_title(title)
    ax.set_ylabel("Fraction" if normalize else ylabel)
    ax.tick_params(axis="x", rotation=25)
    for tick in ax.get_xticklabels():
        tick.set_ha("right")
    ax.legend(fontsize=8, frameon=False, ncols=2)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def save_heatmap(
    path: Path,
    table: dict[str, Counter[str]],
    *,
    title: str,
) -> None:
    if plt is None or not table:
        return
    rows = sorted(table)
    cols = sorted({k for c in table.values() for k in c})
    if not rows or not cols:
        return
    matrix = [[table[row].get(col, 0) for col in cols] for row in rows]
    fig, ax = plt.subplots(figsize=(max(6, 0.7 * len(cols) + 3), max(4, 0.35 * len(rows) + 2)))
    image = ax.imshow(matrix, cmap="YlGnBu", aspect="auto")
    ax.set_title(title)
    ax.set_xticks(range(len(cols)), cols, rotation=30, ha="right")
    ax.set_yticks(range(len(rows)), rows)
    for i, row in enumerate(matrix):
        for j, val in enumerate(row):
            if val:
                ax.text(j, i, str(val), ha="center", va="center", fontsize=8)
    fig.colorbar(image, ax=ax, label="Calls")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def pct(num: int | float, den: int | float) -> str:
    return "0.0" if not den else f"{100.0 * num / den:.1f}"


def build_summary(
    novel_rows: list[dict[str, str]],
    layer_rows: list[dict[str, str]],
    followup_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    total = len(novel_rows)
    pangenome_only = sum(1 for r in novel_rows if call_scope(r) == "pangenome_only")
    single_ref = total - pangenome_only
    read_supported = sum(1 for r in novel_rows if r.get("read_supported", "").lower() == "yes")
    hgt = sum(1 for r in novel_rows if element_group(r.get("element_class")) == "HGT / Starship")
    two_speed = sum(
        1 for r in novel_rows
        if element_group(r.get("element_class")) in {"RIP / genome defense", "TE / repeat"}
    )
    layer_all = next((r for r in layer_rows if r.get("query_asm") == "ALL"), layer_rows[0] if layer_rows else {})
    return [{
        "metric": "novel_mycosv_calls",
        "value": total,
        "percent": "100.0",
        "note": "Rows in novel_mycosv_calls.tsv",
    }, {
        "metric": "pangenome_only_calls",
        "value": pangenome_only,
        "percent": pct(pangenome_only, total),
        "note": "Not representable as single-reference-equivalent calls",
    }, {
        "metric": "single_reference_equivalent_calls",
        "value": single_ref,
        "percent": pct(single_ref, total),
        "note": "Projected to the benchmark reference",
    }, {
        "metric": "read_supported_calls",
        "value": read_supported,
        "percent": pct(read_supported, total),
        "note": "Marked read_supported=yes",
    }, {
        "metric": "hgt_or_starship_candidate_calls",
        "value": hgt,
        "percent": pct(hgt, total),
        "note": "element_class grouped as HGT / Starship",
    }, {
        "metric": "te_or_two_speed_candidate_calls",
        "value": two_speed,
        "percent": pct(two_speed, total),
        "note": "RIP / genome-defense or TE/repeat classes",
    }, {
        "metric": "raw_pairwise_pangenome_observations",
        "value": as_int(layer_all.get("raw_pairwise_pangenome_observations")),
        "percent": "",
        "note": "From pangenome_call_layers.tsv",
    }, {
        "metric": "deduplicated_pangenome_loci",
        "value": as_int(layer_all.get("deduplicated_pangenome_loci")),
        "percent": "",
        "note": "From pangenome_call_layers.tsv",
    }, {
        "metric": "pangenome_powered_loci",
        "value": as_int(layer_all.get("deduplicated_pangenome_loci")),
        "percent": "",
        "note": "Deduplicated graph-native SV loci discovered across the pangenome",
    }, {
        "metric": "pangenome_only_read_supported_loci",
        "value": as_int(layer_all.get("pangenome_only_read_supported")),
        "percent": "",
        "note": "Graph loci not recovered by the single-reference projection but supported by validation evidence",
    }, {
        "metric": "single_ref_fraction_of_raw",
        "value": layer_all.get("single_ref_fraction_of_raw", ""),
        "percent": "",
        "note": "Fraction of raw pangenome observations representable on the benchmark reference",
    }, {
        "metric": "followup_candidates",
        "value": len(followup_rows),
        "percent": "",
        "note": "Rows in mycosv_validation_followup.tsv",
    }]


def make_html(outdir: Path, title: str, summary_rows: list[dict[str, object]], figures: Iterable[Path]) -> None:
    fig_blocks = []
    for fig in figures:
        if fig.exists():
            fig_blocks.append(
                f"<section><h2>{html.escape(fig.stem.replace('_', ' ').title())}</h2>"
                f"<img src='{html.escape(fig.name)}' alt='{html.escape(fig.stem)}'></section>"
            )
    summary_html = "\n".join(
        "<tr>"
        + "".join(f"<td>{html.escape(str(row.get(col, '')))}</td>" for col in ["metric", "value", "percent", "note"])
        + "</tr>"
        for row in summary_rows
    )
    page = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html.escape(title)}</title>
  <style>
    body {{ font-family: system-ui, -apple-system, Segoe UI, sans-serif; margin: 2rem; color: #222; }}
    h1 {{ margin-bottom: 0.25rem; }}
    table {{ border-collapse: collapse; margin: 1rem 0 2rem; }}
    th, td {{ border-bottom: 1px solid #ddd; padding: 0.45rem 0.65rem; text-align: left; }}
    img {{ max-width: 100%; height: auto; border: 1px solid #ddd; }}
    section {{ margin: 2rem 0; }}
  </style>
</head>
<body>
  <h1>{html.escape(title)}</h1>
  <p>Focused plots for MycoSV pangenome-only, single-reference-equivalent, evidence, and fungal biology call layers.</p>
  <h2>Summary</h2>
  <table>
    <thead><tr><th>Metric</th><th>Value</th><th>Percent</th><th>Note</th></tr></thead>
    <tbody>{summary_html}</tbody>
  </table>
  {''.join(fig_blocks)}
</body>
</html>
"""
    (outdir / "mycosv_pangenome_calls_report.html").write_text(page, encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate focused MycoSV pangenome-call plots from benchmark TSV outputs.",
    )
    parser.add_argument("--benchmark-dir", type=Path, required=True,
                        help="Directory containing novel_mycosv_calls.tsv and pangenome_call_layers.tsv.")
    parser.add_argument("--outdir", type=Path, default=None,
                        help="Output directory. Defaults to <benchmark-dir>/pangenome_plots.")
    parser.add_argument("--title", default="MycoSV pangenome call biology plots")
    args = parser.parse_args(argv)

    bench = args.benchmark_dir
    outdir = args.outdir or (bench / "pangenome_plots")
    outdir.mkdir(parents=True, exist_ok=True)

    novel = read_tsv(bench / "novel_mycosv_calls.tsv")
    layers = read_tsv(bench / "pangenome_call_layers.tsv")
    tiers = read_tsv(bench / "mycosv_evidence_tiers.tsv")
    biology = read_tsv(bench / "biology_findings.tsv")
    followup = read_tsv(bench / "mycosv_validation_followup.tsv")

    if not novel and not layers:
        raise SystemExit(f"No pangenome plotting inputs found under {bench}")

    summary_rows = build_summary(novel, layers, followup)
    write_tsv(outdir / "pangenome_summary.tsv", summary_rows, ["metric", "value", "percent", "note"])

    figures: list[Path] = []

    if layers:
        layer = next((r for r in layers if r.get("query_asm") == "ALL"), layers[0])
        labels = [
            "raw observations",
            "deduplicated loci",
            "single-ref equivalent",
            "pangenome-only",
            "pangenome-only read-supported",
        ]
        values = [
            as_int(layer.get("raw_pairwise_pangenome_observations")),
            as_int(layer.get("deduplicated_pangenome_loci")),
            as_int(layer.get("single_reference_equivalent_calls")),
            as_int(layer.get("pangenome_only_calls")),
            as_int(layer.get("pangenome_only_read_supported")),
        ]
        fig = outdir / "pangenome_call_layers.png"
        save_bar(fig, labels, values, title="MycoSV pangenome call layers", rotate=20)
        figures.append(fig)

    if novel:
        scope_counts = Counter(call_scope(r) for r in novel)
        fig = outdir / "single_reference_vs_pangenome_only.png"
        save_bar(
            fig,
            ["single-reference equivalent", "pangenome-only"],
            [scope_counts.get("single_reference_equivalent", 0), scope_counts.get("pangenome_only", 0)],
            title="Single-reference-equivalent vs pangenome-only MycoSV calls",
            colors=[PALETTE["single_reference_equivalent"], PALETTE["pangenome_only"]],
            rotate=10,
        )
        figures.append(fig)

        by_scope_svtype: dict[str, Counter[str]] = defaultdict(Counter)
        by_scope_element: dict[str, Counter[str]] = defaultdict(Counter)
        by_bucket_element: dict[str, Counter[str]] = defaultdict(Counter)
        contig_counts: Counter[str] = Counter()
        for row in novel:
            scope = call_scope(row).replace("_", " ")
            by_scope_svtype[scope][norm_token(row.get("svtype"))] += 1
            by_scope_element[scope][element_group(row.get("element_class"))] += 1
            by_bucket_element[norm_token(row.get("discovery_bucket"))][element_group(row.get("element_class"))] += 1
            if call_scope(row) == "pangenome_only":
                contig_counts[norm_token(row.get("query_contig"))] += 1

        fig = outdir / "svtype_by_pangenome_scope.png"
        save_stacked_bar(fig, by_scope_svtype, title="SV types by pangenome scope")
        figures.append(fig)

        fig = outdir / "biology_class_by_pangenome_scope.png"
        save_stacked_bar(fig, by_scope_element, title="Fungal biology classes by pangenome scope")
        figures.append(fig)

        fig = outdir / "biology_class_by_discovery_bucket.png"
        save_stacked_bar(fig, by_bucket_element, title="Biology classes by MycoSV discovery bucket", normalize=True)
        figures.append(fig)

        top_contigs = top_items(contig_counts, 12)
        fig = outdir / "top_pangenome_only_contigs.png"
        save_bar(
            fig,
            [k for k, _ in top_contigs],
            [v for _, v in top_contigs],
            title="Top contigs by pangenome-only MycoSV calls",
        )
        figures.append(fig)

        rows = [
            {"scope": scope, "category": cat, "n_calls": n}
            for scope, counts in by_scope_element.items()
            for cat, n in sorted(counts.items())
        ]
        write_tsv(outdir / "biology_class_by_pangenome_scope.tsv", rows, ["scope", "category", "n_calls"])

    if tiers:
        tier_table: dict[str, Counter[str]] = defaultdict(Counter)
        for row in tiers:
            svtype = norm_token(row.get("svtype"))
            tier = norm_token(row.get("tier"))
            tier_table[svtype][tier] += as_int(row.get("n_calls"))
        fig = outdir / "evidence_tier_by_svtype.png"
        save_stacked_bar(fig, tier_table, title="MycoSV evidence tiers by SV type")
        figures.append(fig)

    if followup:
        bucket_counts = Counter(norm_token(r.get("validation_bucket")) for r in followup)
        top_buckets = top_items(bucket_counts, 10)
        fig = outdir / "validation_followup_buckets.png"
        save_bar(
            fig,
            [k for k, _ in top_buckets],
            [v for _, v in top_buckets],
            title="Validation follow-up buckets for MycoSV novel calls",
        )
        figures.append(fig)

    if biology:
        ecology_table: dict[str, Counter[str]] = defaultdict(Counter)
        candidate_counts: Counter[str] = Counter()
        for row in biology:
            trophic = norm_token(row.get("trophic_mode") or row.get("ecological_trait"))
            candidate = norm_token(row.get("candidate_type"))
            group = element_group(row.get("element_class"))
            ecology_table[trophic][group] += 1
            candidate_counts[candidate] += 1
        fig = outdir / "ecology_by_genome_biology_heatmap.png"
        save_heatmap(fig, ecology_table, title="Ecological mode by MycoSV genome-biology class")
        figures.append(fig)
        top_candidates = top_items(candidate_counts, 12)
        fig = outdir / "biology_candidate_types.png"
        save_bar(
            fig,
            [k for k, _ in top_candidates],
            [v for _, v in top_candidates],
            title="Top MycoSV fungal biology candidate types",
        )
        figures.append(fig)

    make_html(outdir, args.title, summary_rows, figures)
    print(f"Wrote MycoSV pangenome plots to: {outdir}")
    print(f"HTML report: {outdir / 'mycosv_pangenome_calls_report.html'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
