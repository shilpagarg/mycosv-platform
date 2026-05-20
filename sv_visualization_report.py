#!/usr/bin/env python3
"""
Create an integrated visualization report for structural variant (SV) analyses
across simulated and real datasets.

What this script covers
-----------------------
1. Simulated data benchmarking
   - Precision / recall / F1 by caller, dataset, and SV type
   - Breakpoint error distributions
   - Size-stratified performance
   - Support counts and truth overlap summaries

2. Real data summaries
   - SV burden by sample / cohort
   - SV type composition
   - Length distributions
   - Chromosome-level SV landscapes
   - Caller overlap on real data

3. Biological findings
   - Recurrent genes affected by SVs
   - Pathway / annotation summaries
   - Sample-level burden vs phenotype associations
   - Optional gene expression or copy-number linked summaries

Outputs
-------
- A self-contained HTML report
- PNG figures in an output directory
- Summary TSV files for downstream use

Expected input
--------------
The script is intentionally permissive. It accepts CSV/TSV tables and tries to
harmonize common column names.

Recommended files:
- simulated_benchmark.tsv
- real_sv_calls.tsv
- biological_findings.tsv
- sample_metadata.tsv (optional)

Example usage
-------------
python sv_visualization_report.py \
  --simulated results/simulated_benchmark.tsv \
  --real results/real_sv_calls.tsv \
  --biology results/biological_findings.tsv \
  --metadata results/sample_metadata.tsv \
  --outdir report_output \
  --title "SV analysis report"

Notes
-----
- Uses matplotlib only for plotting.
- Avoids hard-coding a single schema by normalizing synonymous column names.
- If some sections cannot be generated because inputs are missing, the report is
  still created with the available content.
"""

from __future__ import annotations

import argparse
import base64
import io
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

try:
    import numpy as np
    import pandas as pd
except ModuleNotFoundError:
    fallback_python = Path(
        os.environ.get(
            "MYCOSV_ENV_PYTHON",
            "/mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/tools/envs/envs/fungi_graph_sv/bin/python",
        )
    )
    if fallback_python.exists() and Path(sys.executable).resolve() != fallback_python.resolve():
        os.execv(str(fallback_python), [str(fallback_python), *sys.argv])
    raise

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None


# -----------------------------
# Utilities
# -----------------------------


def read_table(path: Optional[str]) -> Optional[pd.DataFrame]:
    if path is None:
        return None
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Input file not found: {p}")
    if p.stat().st_size == 0:
        return pd.DataFrame()
    sep = "\t" if p.suffix.lower() in {".tsv", ".txt"} else ","
    # Tolerate ragged rows that slip through upstream merging — pandas' default
    # C parser bails on the first row whose field count differs from the
    # header. Falling back to the python engine with on_bad_lines='skip'
    # drops only the offending rows instead of aborting the whole report.
    # low_memory=False reads each column in a single pass so pandas can pick a
    # consistent dtype; the previous chunked read emitted a DtypeWarning on the
    # 30+ mixed-type columns produced by the comparator merge.
    try:
        return pd.read_csv(p, sep=sep, low_memory=False)
    except pd.errors.ParserError:
        sys.stderr.write(
            f"[read_table] Falling back to lenient parse for {p} "
            "(skipping malformed rows)\n"
        )
        return pd.read_csv(p, sep=sep, engine="python", on_bad_lines="skip")


COLUMN_ALIASES: Dict[str, Sequence[str]] = {
    "sample": ["sample", "sample_id", "tumor_sample", "specimen", "case", "query_asm"],
    "dataset": ["dataset", "dataset_name", "cohort", "group"],
    "caller": ["caller", "tool", "method", "sv_caller", "alignment_mode"],
    "sv_type": ["sv_type", "type", "svclass", "sv_class", "event_type", "svtype"],
    "chrom": ["chrom", "chr", "chromosome", "query_contig"],
    "start": ["start", "pos", "position", "breakpoint1", "bp1"],
    "end": ["end", "stop", "breakpoint2", "bp2"],
    "sv_len": ["sv_len", "length", "size", "event_size", "span", "ancestral_segment_bp"],
    "clade": ["clade", "phylum", "clade_rank", "cladeRank", "taxonomy"],
    "te_class": ["te_class", "element_class", "elementClass", "repeat_class", "te_type"],
    "novelty_tier": ["novelty_tier", "novelty", "hgt_tier", "annotation_tier"],
    "precision": ["precision", "prec"],
    "recall": ["recall", "sensitivity", "tpr"],
    "f1": ["f1", "f1_score"],
    "tp": ["tp", "true_positive", "true_positives"],
    "fp": ["fp", "false_positive", "false_positives"],
    "fn": ["fn", "false_negative", "false_negatives"],
    "breakpoint_error": ["breakpoint_error", "bp_error", "distance_to_truth", "expression_distance_bp"],
    "support": ["support", "read_support", "split_read_support", "evidence", "comparator_support_count"],
    "truth_overlap": ["truth_overlap", "matched_truth", "is_matched", "overlap_truth"],
    # Biology candidates columns from analyze_new_biology_candidates.py
    "gene": ["gene", "gene_symbol", "affected_gene", "nearest_gene",
             "expression_gene", "gene_name", "gene_id"],
    "pathway": ["pathway", "term", "annotation", "hallmark",
                "evidence_axis", "functional_example", "candidate_type"],
    "effect": ["effect", "impact", "functional_effect", "consequence",
               "element_class", "novelty", "rationale"],
    "phenotype": ["phenotype", "status", "response", "subtype",
                  "scenario", "lifestyle", "architecture"],
    "expression": ["expression", "expr", "log2_expression", "gene_expression",
                   "expression_log2_fc", "log2_fc"],
    "cn": ["cn", "copy_number", "copy_number_state", "priority"],
}


def harmonize_columns(df: pd.DataFrame) -> pd.DataFrame:
    if df is None:
        return df
    out = df.copy()
    lower_map = {c.lower(): c for c in df.columns}
    for canonical, aliases in COLUMN_ALIASES.items():
        if canonical in out.columns:
            continue
        for alias in aliases:
            if alias.lower() in lower_map:
                out[canonical] = out[lower_map[alias.lower()]]
                break

    if "sv_len" not in out.columns and {"start", "end"}.issubset(out.columns):
        out["sv_len"] = (pd.to_numeric(out["end"], errors="coerce") - pd.to_numeric(out["start"], errors="coerce")).abs()

    if "f1" not in out.columns and {"precision", "recall"}.issubset(out.columns):
        p = pd.to_numeric(out["precision"], errors="coerce")
        r = pd.to_numeric(out["recall"], errors="coerce")
        out["f1"] = np.where((p + r) > 0, 2 * p * r / (p + r), np.nan)

    for col in ["precision", "recall", "f1", "tp", "fp", "fn", "breakpoint_error", "support", "sv_len", "expression", "cn"]:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")

    return out



def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)



def save_df(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False)


def plotting_available() -> bool:
    return plt is not None


def plotting_unavailable_html() -> str:
    return (
        "<p><i>Plotting skipped because matplotlib is not installed in this "
        "Python environment. Summary tables and the HTML report were still "
        "generated.</i></p>"
    )



def fig_to_base64(fig: plt.Figure) -> str:
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    buffer.seek(0)
    return base64.b64encode(buffer.read()).decode("utf-8")



def write_figure(fig: plt.Figure, path: Path) -> str:
    fig.savefig(path, dpi=180, bbox_inches="tight")
    encoded = fig_to_base64(fig)
    return encoded



def safe_groupby_mean(df: pd.DataFrame, group_cols: List[str], value_cols: List[str]) -> pd.DataFrame:
    cols = [c for c in value_cols if c in df.columns]
    if not cols:
        return pd.DataFrame()
    return df.groupby(group_cols, dropna=False)[cols].mean(numeric_only=True).reset_index()



NULL_LABELS = frozenset({"", ".", "NA", "N/A", "NONE", "NULL", "NAN", "NONE_SET", "UNKNOWN"})


def non_null_label_mask(series: pd.Series) -> pd.Series:
    labels = series.fillna("").astype(str).str.strip()
    return ~labels.str.upper().isin(NULL_LABELS)


def clean_label_series(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def first_non_null_label(df: pd.DataFrame, columns: Sequence[str], default: str = "") -> pd.Series:
    out = pd.Series(default, index=df.index, dtype="object")
    unresolved = pd.Series(True, index=df.index)
    for col in columns:
        if col not in df.columns:
            continue
        values = clean_label_series(df[col])
        use = unresolved & non_null_label_mask(values)
        out.loc[use] = values.loc[use]
        unresolved &= ~use
    return out


def top_n(df: pd.DataFrame, col: str, n: int = 20, *, drop_null_like: bool = True) -> pd.DataFrame:
    values = clean_label_series(df[col])
    if drop_null_like:
        values = values[non_null_label_mask(values)]
    vc = values.value_counts().head(n)
    return vc.rename_axis(col).reset_index(name="count")


def biology_effect_plot_table(bio_df: pd.DataFrame, top_n_labels: int = 15) -> pd.DataFrame:
    """Rows for the biology_gene_effects plot.

    Rank by informative repeat/mobile-element effects, not by total gene
    counts. Rows such as TE_LINE can be locus-level biology candidates with
    no affected gene assignment, so fall back to affected_locus/query_contig
    instead of dropping them from the effect-composition panel.
    """
    needed = {"effect"}
    if bio_df is None or bio_df.empty or not needed.issubset(bio_df.columns):
        return pd.DataFrame()
    out = bio_df.copy()
    out["_effect_label"] = first_non_null_label(
        out,
        ["gene", "affected_gene", "expression_gene", "affected_locus", "query_contig"],
    )
    out["_effect_class"] = clean_label_series(out["effect"])
    out = out[
        non_null_label_mask(out["_effect_label"])
        & non_null_label_mask(out["_effect_class"])
    ].copy()
    if out.empty:
        return pd.DataFrame()
    label_counts = out["_effect_label"].value_counts()
    selected = list(label_counts.head(top_n_labels).index)

    # Preserve rare informative classes (for example TE_LINE) even when many
    # HGT/RIP loci tie ahead of them by count.
    for effect in sorted(out["_effect_class"].unique()):
        best_for_effect = (
            out[out["_effect_class"].eq(effect)]["_effect_label"]
            .value_counts()
            .index
        )
        if len(best_for_effect) and best_for_effect[0] not in selected:
            selected.append(best_for_effect[0])

    return out[out["_effect_label"].isin(selected)].copy()


def numeric_pair_table(df: pd.DataFrame, x: str, y: str) -> pd.DataFrame:
    if df is None or df.empty or not {x, y}.issubset(df.columns):
        return pd.DataFrame(columns=[x, y])
    out = df[[x, y]].copy()
    out[x] = pd.to_numeric(out[x], errors="coerce")
    out[y] = pd.to_numeric(out[y], errors="coerce")
    return out.replace([np.inf, -np.inf], np.nan).dropna(subset=[x, y])



def add_size_bin(df: pd.DataFrame, source_col: str = "sv_len") -> pd.DataFrame:
    if source_col not in df.columns:
        return df
    bins = [-np.inf, 50, 100, 500, 1_000, 10_000, 100_000, np.inf]
    labels = ["<50", "50-100", "100-500", "500-1k", "1k-10k", "10k-100k", ">100k"]
    out = df.copy()
    out["size_bin"] = pd.cut(out[source_col], bins=bins, labels=labels)
    return out


SVTYPE_ORDER: tuple[str, ...] = ("INS", "DEL", "INV", "DUP", "TRA", "OFF_REF")


def short_sample_label(value: object, max_len: int = 30) -> str:
    s = str(value)
    if len(s) <= max_len:
        return s
    return s[: max_len - 1] + "…"


def is_benchmark_summary_table(df: pd.DataFrame) -> bool:
    return {"query_asm", "method", "svtype", "pred_calls", "truth_label"}.issubset(df.columns)


def benchmark_pred_call_counts(df: pd.DataFrame, *, include_all: bool = False) -> pd.DataFrame:
    """Deduplicate exact_benchmark_summary-style rows into prediction counts.

    exact_benchmark_summary.tsv repeats the same prediction counts once per
    truth label. For landscape plots we want one count per sample/caller/SV
    type, so max(pred_calls) is the stable de-duplicated value.
    """
    if not is_benchmark_summary_table(df):
        return pd.DataFrame()
    tmp = df.copy()
    tmp["pred_calls"] = pd.to_numeric(tmp["pred_calls"], errors="coerce").fillna(0)
    tmp["svtype"] = tmp["svtype"].astype(str)
    if "coordinate_space" in tmp.columns:
        coord = tmp["coordinate_space"].astype(str).str.lower()
        tmp = tmp[coord.ne("reference_any_clade")].copy()
        coord = tmp["coordinate_space"].astype(str).str.lower()
        if coord.eq("query").any():
            tmp = tmp[coord.eq("query")].copy()
    if not include_all:
        tmp = tmp[tmp["svtype"] != "ALL"]
    keep_methods = tmp["method"].astype(str).ne("")
    tmp = tmp[keep_methods]
    if tmp.empty:
        return pd.DataFrame()
    out = (
        tmp.groupby(["query_asm", "method", "svtype"], dropna=False)["pred_calls"]
        .max()
        .reset_index(name="count")
        .rename(columns={"query_asm": "sample", "method": "caller", "svtype": "sv_type"})
    )
    out["count"] = pd.to_numeric(out["count"], errors="coerce").fillna(0)
    return out


@dataclass
class FigureRecord:
    title: str
    filename: str
    encoded_png: str
    caption: str


# -----------------------------
# Plotting helpers
# -----------------------------


def plot_bar(df: pd.DataFrame, x: str, y: str, title: str, xlabel: str, ylabel: str, rotate_xticks: bool = False) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(max(9, 0.34 * len(df) + 3), 5.5))
    labels = [short_sample_label(v) for v in df[x]]
    ax.bar(labels, df[y])
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if rotate_xticks:
        ax.tick_params(axis="x", rotation=35, labelsize=8)
        for label in ax.get_xticklabels():
            label.set_ha("right")
    fig.tight_layout()
    return fig



def plot_grouped_metric(df: pd.DataFrame, category: str, metric: str, series: str, title: str) -> plt.Figure:
    pivot = df.pivot_table(index=category, columns=series, values=metric, aggfunc="mean")
    pivot.index = [short_sample_label(v) for v in pivot.index]
    fig, ax = plt.subplots(figsize=(max(10, 0.34 * len(pivot.index) + 3), 5.5))
    pivot.plot(kind="bar", ax=ax)
    ax.set_title(title)
    ax.set_xlabel(category)
    ax.set_ylabel(metric)
    ax.legend(title=series, bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.tick_params(axis="x", rotation=35, labelsize=8)
    for label in ax.get_xticklabels():
        label.set_ha("right")
    fig.tight_layout()
    return fig


def plot_caller_count_proxy(df: pd.DataFrame, title: str, *, top_n: int = 25) -> plt.Figure:
    totals = df.groupby("sample", dropna=False)["count"].sum().sort_values(ascending=False)
    top_samples = list(totals.head(top_n).index)
    plot_df = df[df["sample"].isin(top_samples)].copy()
    plot_df["sample"] = pd.Categorical(plot_df["sample"], categories=top_samples, ordered=True)
    pivot = (
        plot_df.pivot_table(index="sample", columns="caller", values="count", aggfunc="sum", observed=False)
        .fillna(0)
        .reindex(top_samples)
    )
    pivot.index = [short_sample_label(v) for v in pivot.index]
    fig, ax = plt.subplots(figsize=(max(11, 0.42 * len(pivot.index) + 4), 5.8))
    pivot.plot(kind="bar", ax=ax, width=0.84)
    ax.set_yscale("symlog", linthresh=5)
    ax.set_title(title)
    ax.set_xlabel("Sample")
    ax.set_ylabel("SV count (symlog scale)")
    ax.legend(title="caller", bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.tick_params(axis="x", rotation=35, labelsize=8)
    for label in ax.get_xticklabels():
        label.set_ha("right")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    return fig



def plot_hist(values: pd.Series, title: str, xlabel: str, bins: int = 30, log_x: bool = False) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(8, 5))
    vals = values.dropna()
    if log_x:
        vals = vals[vals > 0]
        if len(vals) == 0:
            vals = pd.Series([1])
        bins_edges = np.logspace(np.log10(vals.min()), np.log10(vals.max()), bins)
        ax.hist(vals, bins=bins_edges)
        ax.set_xscale("log")
    else:
        ax.hist(vals, bins=bins)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    fig.tight_layout()
    return fig


def boxplot_with_labels(ax: plt.Axes, data: list, labels: list, **kwargs) -> None:
    try:
        ax.boxplot(data, tick_labels=labels, **kwargs)
    except TypeError:
        ax.boxplot(data, labels=labels, **kwargs)



def plot_box(df: pd.DataFrame, x: str, y: str, title: str, rotate_xticks: bool = False) -> plt.Figure:
    categories = [g[y].dropna().values for _, g in df.groupby(x)]
    labels = [str(k) for k, _ in df.groupby(x)]
    fig, ax = plt.subplots(figsize=(10, 5))
    boxplot_with_labels(ax, categories, labels, showfliers=False)
    ax.set_title(title)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    if rotate_xticks:
        ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    return fig



def plot_stacked_counts(df: pd.DataFrame, index_col: str, stack_col: str, title: str, normalize: bool = False) -> plt.Figure:
    counts = pd.crosstab(df[index_col].fillna("NA"), df[stack_col].fillna("NA"), normalize="index" if normalize else False)
    counts = counts[[c for c in SVTYPE_ORDER if c in counts.columns] + [c for c in counts.columns if c not in SVTYPE_ORDER]]
    counts.index = [short_sample_label(v) for v in counts.index]
    fig, ax = plt.subplots(figsize=(max(10, 0.34 * len(counts.index) + 3), 5.5))
    counts.plot(kind="bar", stacked=True, ax=ax)
    ax.set_title(title)
    ax.set_xlabel(index_col)
    ax.set_ylabel("Fraction" if normalize else "Count")
    ax.legend(title=stack_col, bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.tick_params(axis="x", rotation=35, labelsize=8)
    for label in ax.get_xticklabels():
        label.set_ha("right")
    fig.tight_layout()
    return fig



def plot_scatter(df: pd.DataFrame, x: str, y: str, title: str, xlabel: Optional[str] = None, ylabel: Optional[str] = None) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(7, 5))
    clean = numeric_pair_table(df, x, y)
    ax.scatter(clean[x], clean[y], alpha=0.7)
    ax.set_title(title)
    ax.set_xlabel(xlabel or x)
    ax.set_ylabel(ylabel or y)
    fig.tight_layout()
    return fig


# -----------------------------
# Reporting sections
# -----------------------------


def build_simulated_section(sim_df: Optional[pd.DataFrame], outdir: Path) -> Tuple[str, List[FigureRecord], List[pd.DataFrame]]:
    if sim_df is None or sim_df.empty:
        return "<p>No simulated benchmark table provided.</p>", [], []

    sim_df = harmonize_columns(sim_df)
    sim_df = add_size_bin(sim_df)
    figs: List[FigureRecord] = []
    tables: List[pd.DataFrame] = []
    blocks: List[str] = []

    group_cols = [c for c in ["caller", "dataset", "sv_type"] if c in sim_df.columns]
    metric_cols = [c for c in ["precision", "recall", "f1", "tp", "fp", "fn"] if c in sim_df.columns]
    if group_cols and metric_cols:
        summary = safe_groupby_mean(sim_df, group_cols, metric_cols)
        if not summary.empty:
            tables.append(summary)
            save_df(summary, outdir / "simulated_summary.tsv")
            preview = summary.head(20).to_html(index=False, border=0)
            blocks.append("<h3>Simulated benchmark summary</h3>" + preview)

    if not plotting_available():
        blocks.append(plotting_unavailable_html())
        blocks.append("<p>This section summarizes caller performance on simulated truth sets, including accuracy, breakpoint precision, and size-stratified sensitivity.</p>")
        return "\n".join(blocks), figs, tables

    if {"sv_type", "f1", "caller"}.issubset(sim_df.columns):
        fig = plot_grouped_metric(sim_df, "sv_type", "f1", "caller", "F1 by SV type and caller (simulated)")
        encoded = write_figure(fig, outdir / "sim_f1_by_svtype_caller.png")
        figs.append(FigureRecord(
            title="F1 by SV type and caller",
            filename="sim_f1_by_svtype_caller.png",
            encoded_png=encoded,
            caption="Average F1 across simulated truth sets stratified by SV type and caller."
        ))

    if {"size_bin", "recall", "caller"}.issubset(sim_df.columns):
        size_perf = sim_df.dropna(subset=["size_bin"]).groupby(["size_bin", "caller"], observed=False)["recall"].mean().reset_index()
        fig = plot_grouped_metric(size_perf, "size_bin", "recall", "caller", "Recall by SV size bin (simulated)")
        encoded = write_figure(fig, outdir / "sim_recall_by_size.png")
        figs.append(FigureRecord(
            title="Recall by SV size bin",
            filename="sim_recall_by_size.png",
            encoded_png=encoded,
            caption="Size-stratified recall highlights caller sensitivity across small and large events."
        ))

    if "breakpoint_error" in sim_df.columns:
        fig = plot_hist(sim_df["breakpoint_error"], "Breakpoint error distribution (simulated)", "Breakpoint error (bp)")
        encoded = write_figure(fig, outdir / "sim_breakpoint_error_hist.png")
        figs.append(FigureRecord(
            title="Breakpoint error distribution",
            filename="sim_breakpoint_error_hist.png",
            encoded_png=encoded,
            caption="Distribution of deviation between called and truth breakpoints."
        ))

    if {"caller", "support"}.issubset(sim_df.columns):
        fig = plot_box(sim_df, "caller", "support", "Support by caller for matched simulated SVs", rotate_xticks=True)
        encoded = write_figure(fig, outdir / "sim_support_by_caller.png")
        figs.append(FigureRecord(
            title="Support by caller",
            filename="sim_support_by_caller.png",
            encoded_png=encoded,
            caption="Distribution of supporting evidence for simulated SV calls across callers."
        ))

    blocks.append("<p>This section summarizes caller performance on simulated truth sets, including accuracy, breakpoint precision, and size-stratified sensitivity.</p>")
    return "\n".join(blocks), figs, tables



def build_real_section(real_df: Optional[pd.DataFrame], metadata_df: Optional[pd.DataFrame], outdir: Path) -> Tuple[str, List[FigureRecord], List[pd.DataFrame]]:
    if real_df is None or real_df.empty:
        return "<p>No real-data SV table provided.</p>", [], []

    raw_real_df = real_df.copy()
    real_df = harmonize_columns(real_df)
    real_df = add_size_bin(real_df)
    benchmark_counts = benchmark_pred_call_counts(raw_real_df)
    figs: List[FigureRecord] = []
    tables: List[pd.DataFrame] = []
    blocks: List[str] = []

    if metadata_df is not None and not metadata_df.empty:
        metadata_df = harmonize_columns(metadata_df)
        if "sample" in real_df.columns and "sample" in metadata_df.columns:
            real_df = real_df.merge(metadata_df, on="sample", how="left", suffixes=("", "_meta"))

    if not plotting_available():
        if "sample" in real_df.columns:
            burden = real_df.groupby("sample", dropna=False).size().reset_index(name="sv_count").sort_values("sv_count", ascending=False)
            tables.append(burden)
            save_df(burden, outdir / "real_sv_burden_by_sample.tsv")
        if {"chrom", "sv_type"}.issubset(real_df.columns):
            chrom_counts = real_df.groupby(["chrom", "sv_type"], dropna=False).size().reset_index(name="count")
            tables.append(chrom_counts)
            save_df(chrom_counts, outdir / "real_sv_by_chrom_and_type.tsv")
        blocks.append(plotting_unavailable_html())
        blocks.append("<p>This section summarizes real-data structural variant landscapes, sample burden, class composition, size distributions, and cohort-level heterogeneity.</p>")
        return "\n".join(blocks), figs, tables

    if not benchmark_counts.empty:
        mycosv_counts = benchmark_counts[benchmark_counts["caller"].eq("mycosv")]
        burden_source = mycosv_counts if not mycosv_counts.empty else benchmark_counts
        burden = (
            burden_source.groupby("sample", dropna=False)["count"]
            .sum()
            .reset_index(name="sv_count")
            .sort_values("sv_count", ascending=False)
        )
        tables.append(burden)
        save_df(burden, outdir / "real_sv_burden_by_sample.tsv")
        fig = plot_bar(burden.head(30), "sample", "sv_count", "MycoSV SV burden by sample (top 30)", "Sample", "SV count", rotate_xticks=True)
        encoded = write_figure(fig, outdir / "real_sv_burden_by_sample.png")
        figs.append(FigureRecord(
            title="SV burden by sample",
            filename="real_sv_burden_by_sample.png",
            encoded_png=encoded,
            caption="De-duplicated MycoSV prediction burden across real samples, using pred_calls from benchmark summaries."
        ))
    elif "sample" in real_df.columns:
        burden = real_df.groupby("sample", dropna=False).size().reset_index(name="sv_count").sort_values("sv_count", ascending=False)
        tables.append(burden)
        save_df(burden, outdir / "real_sv_burden_by_sample.tsv")
        fig = plot_bar(burden.head(30), "sample", "sv_count", "SV burden by sample (top 30)", "Sample", "SV count", rotate_xticks=True)
        encoded = write_figure(fig, outdir / "real_sv_burden_by_sample.png")
        figs.append(FigureRecord(
            title="SV burden by sample",
            filename="real_sv_burden_by_sample.png",
            encoded_png=encoded,
            caption="SV burden across real samples, showing the most variant-rich samples."
        ))

    if not benchmark_counts.empty:
        mycosv_counts = benchmark_counts[benchmark_counts["caller"].eq("mycosv")]
        comp_source = mycosv_counts if not mycosv_counts.empty else benchmark_counts
        comp_source = comp_source[comp_source["count"] > 0]
        if not comp_source.empty:
            fig = plot_stacked_counts(
                comp_source.loc[comp_source.index.repeat(comp_source["count"].astype(int).clip(lower=0))],
                "sample", "sv_type", "MycoSV SV composition by sample", normalize=True,
            )
            encoded = write_figure(fig, outdir / "real_sv_composition_by_sample.png")
            figs.append(FigureRecord(
                title="SV composition by sample",
                filename="real_sv_composition_by_sample.png",
                encoded_png=encoded,
                caption="Relative MycoSV SV class composition per sample; ALL rows and repeated truth-label rows are excluded."
            ))
    elif {"sample", "sv_type"}.issubset(real_df.columns):
        real_sv_only = real_df[real_df["sv_type"].astype(str).ne("ALL")]
        fig = plot_stacked_counts(real_sv_only, "sample", "sv_type", "SV composition by sample", normalize=True)
        encoded = write_figure(fig, outdir / "real_sv_composition_by_sample.png")
        figs.append(FigureRecord(
            title="SV composition by sample",
            filename="real_sv_composition_by_sample.png",
            encoded_png=encoded,
            caption="Relative composition of SV classes for each sample."
        ))

    if "sv_len" in real_df.columns:
        fig = plot_hist(real_df["sv_len"], "SV length distribution (real data)", "SV length (bp)", log_x=True)
        encoded = write_figure(fig, outdir / "real_sv_length_distribution.png")
        figs.append(FigureRecord(
            title="SV length distribution",
            filename="real_sv_length_distribution.png",
            encoded_png=encoded,
            caption="Length distribution of real-data SVs on a logarithmic x-axis."
        ))

    if {"chrom", "sv_type"}.issubset(real_df.columns):
        chrom_counts = real_df.groupby(["chrom", "sv_type"], dropna=False).size().reset_index(name="count")
        tables.append(chrom_counts)
        save_df(chrom_counts, outdir / "real_sv_by_chrom_and_type.tsv")
        fig = plot_grouped_metric(chrom_counts, "chrom", "count", "sv_type", "SV landscape by chromosome")
        encoded = write_figure(fig, outdir / "real_sv_by_chromosome.png")
        figs.append(FigureRecord(
            title="SV landscape by chromosome",
            filename="real_sv_by_chromosome.png",
            encoded_png=encoded,
            caption="Chromosome-level distribution of SV classes in real data."
        ))

    if not benchmark_counts.empty:
        caller_summary = (
            benchmark_counts.groupby(["sample", "caller"], dropna=False)["count"]
            .sum()
            .reset_index(name="count")
        )
        fig = plot_caller_count_proxy(caller_summary, "Caller-specific SV counts in real data (top 25 samples)")
        encoded = write_figure(fig, outdir / "real_caller_overlap_proxy.png")
        figs.append(FigureRecord(
            title="Caller-specific counts",
            filename="real_caller_overlap_proxy.png",
            encoded_png=encoded,
            caption="De-duplicated caller-specific prediction counts from pred_calls; y-axis uses symlog scaling so large comparator callsets do not hide smaller callsets."
        ))
    elif {"caller", "sample"}.issubset(real_df.columns):
        caller_summary = real_df.groupby(["sample", "caller"], dropna=False).size().reset_index(name="count")
        fig = plot_grouped_metric(caller_summary, "sample", "count", "caller", "Caller-specific SV counts in real data")
        encoded = write_figure(fig, outdir / "real_caller_overlap_proxy.png")
        figs.append(FigureRecord(
            title="Caller-specific counts",
            filename="real_caller_overlap_proxy.png",
            encoded_png=encoded,
            caption="A practical proxy for caller concordance in real data when truth labels are unavailable."
        ))

    if {"phenotype", "sample"}.issubset(real_df.columns):
        sample_burden = real_df.groupby(["sample", "phenotype"], dropna=False).size().reset_index(name="sv_count")
        fig = plot_box(sample_burden, "phenotype", "sv_count", "SV burden by phenotype", rotate_xticks=True)
        encoded = write_figure(fig, outdir / "real_sv_burden_by_phenotype.png")
        figs.append(FigureRecord(
            title="SV burden by phenotype",
            filename="real_sv_burden_by_phenotype.png",
            encoded_png=encoded,
            caption="Distribution of overall SV burden across phenotype-defined groups."
        ))

    blocks.append("<p>This section summarizes real-data structural variant landscapes, sample burden, class composition, size distributions, and cohort-level heterogeneity.</p>")
    return "\n".join(blocks), figs, tables



def build_biology_section(bio_df: Optional[pd.DataFrame], outdir: Path) -> Tuple[str, List[FigureRecord], List[pd.DataFrame]]:
    if bio_df is None or bio_df.empty:
        return "<p>No biological findings table provided.</p>", [], []

    bio_df = harmonize_columns(bio_df)
    figs: List[FigureRecord] = []
    tables: List[pd.DataFrame] = []
    blocks: List[str] = []

    if not plotting_available():
        if "gene" in bio_df.columns:
            gene_counts = top_n(bio_df, "gene", 20)
            tables.append(gene_counts)
            save_df(gene_counts, outdir / "biology_top_genes.tsv")
        if "pathway" in bio_df.columns:
            pathway_counts = top_n(bio_df, "pathway", 15)
            tables.append(pathway_counts)
            save_df(pathway_counts, outdir / "biology_top_pathways.tsv")
        blocks.append(plotting_unavailable_html())
        blocks.append("<p>This section summarizes recurrent genes, pathways, functional effects, and expression/copy-number links associated with SVs.</p>")
        return "\n".join(blocks), figs, tables

    if "gene" in bio_df.columns:
        gene_counts = top_n(bio_df, "gene", 20)
        tables.append(gene_counts)
        save_df(gene_counts, outdir / "biology_top_genes.tsv")
        if not gene_counts.empty:
            fig = plot_bar(gene_counts, "gene", "count", "Top recurrent genes affected by SVs", "Gene", "Event count", rotate_xticks=True)
            encoded = write_figure(fig, outdir / "biology_top_genes.png")
            figs.append(FigureRecord(
                title="Top recurrent genes",
                filename="biology_top_genes.png",
                encoded_png=encoded,
                caption="Most recurrent genes linked to called structural variants; placeholder labels such as NA or . are excluded."
            ))

    if "pathway" in bio_df.columns:
        pathway_counts = top_n(bio_df, "pathway", 15)
        tables.append(pathway_counts)
        save_df(pathway_counts, outdir / "biology_top_pathways.tsv")
        fig = plot_bar(pathway_counts, "pathway", "count", "Top affected pathways / annotations", "Pathway", "Count", rotate_xticks=True)
        encoded = write_figure(fig, outdir / "biology_top_pathways.png")
        figs.append(FigureRecord(
            title="Top pathways / annotations",
            filename="biology_top_pathways.png",
            encoded_png=encoded,
            caption="Most frequently implicated pathways or annotations among biological findings."
        ))

    if "effect" in bio_df.columns:
        subset = biology_effect_plot_table(bio_df, 15)
        if not subset.empty:
            effect_summary = (
                subset.groupby(["_effect_label", "_effect_class"], dropna=False)
                .size()
                .reset_index(name="count")
                .rename(columns={"_effect_label": "gene_or_locus", "_effect_class": "effect"})
            )
            save_df(effect_summary, outdir / "biology_gene_effects.tsv")
            fig = plot_stacked_counts(
                subset,
                "_effect_label",
                "_effect_class",
                "Functional effect composition of top genes/loci",
                normalize=True,
            )
            encoded = write_figure(fig, outdir / "biology_gene_effects.png")
            figs.append(FigureRecord(
                title="Functional effect composition",
                filename="biology_gene_effects.png",
                encoded_png=encoded,
                caption=(
                    "Relative composition of informative repeat/mobile-element "
                    "effects among recurrent genes or loci. Rows with placeholder "
                    "effects such as NONE are excluded; locus labels are used "
                    "when no affected gene was assigned."
                )
            ))

    if {"expression", "cn"}.issubset(bio_df.columns):
        cn_expr = numeric_pair_table(bio_df, "cn", "expression")
        save_df(cn_expr, outdir / "biology_cn_vs_expression.tsv")
        if not cn_expr.empty:
            fig = plot_scatter(cn_expr, "cn", "expression", "Copy number vs expression for SV-associated loci")
            encoded = write_figure(fig, outdir / "biology_cn_vs_expression.png")
            figs.append(FigureRecord(
                title="Copy number vs expression",
                filename="biology_cn_vs_expression.png",
                encoded_png=encoded,
                caption="Association between copy-number state and expression at SV-associated loci."
            ))
        else:
            blocks.append(
                "<p>Copy number vs expression plot skipped: no rows have both "
                "numeric copy-number/priority and expression log2 fold-change "
                "values.</p>"
            )

    if {"phenotype", "gene"}.issubset(bio_df.columns):
        gene_df = bio_df[non_null_label_mask(bio_df["gene"])].copy()
        top_genes = clean_label_series(gene_df["gene"]).value_counts().head(10).index
        subset = gene_df[clean_label_series(gene_df["gene"]).isin(top_genes)].copy()
        if not subset.empty:
            subset["gene"] = clean_label_series(subset["gene"])
            phen = pd.crosstab(subset["phenotype"].fillna("NA"), subset["gene"])
            phen = phen.reset_index().melt(id_vars="phenotype", var_name="gene", value_name="count")
            fig = plot_grouped_metric(phen, "phenotype", "count", "gene", "Top recurrent genes across phenotypes")
            encoded = write_figure(fig, outdir / "biology_genes_by_phenotype.png")
            figs.append(FigureRecord(
                title="Genes across phenotypes",
                filename="biology_genes_by_phenotype.png",
                encoded_png=encoded,
                caption="Distribution of recurrent SV-linked genes across phenotype groups."
            ))

    blocks.append("<p>This section focuses on biological interpretation, highlighting recurrent genes, pathways, functional effects, and optional multi-omic associations.</p>")
    return "\n".join(blocks), figs, tables


# -----------------------------
# New biology plots: clade-SV, TE-architecture, HGT-propagation
# -----------------------------


def plot_clade_sv(df: pd.DataFrame, outdir: Path) -> Optional[FigureRecord]:
    """SV type counts stratified by fungal clade/phylum."""
    needed = {"clade", "sv_type"}
    if not needed.issubset(df.columns):
        return None
    counts = df.groupby(["clade", "sv_type"], dropna=False).size().reset_index(name="count")
    if counts.empty:
        return None
    fig, ax = plt.subplots(figsize=(11, 5))
    pivot = counts.pivot_table(index="clade", columns="sv_type", values="count",
                               aggfunc="sum", fill_value=0)
    pivot.plot(kind="bar", stacked=True, ax=ax)
    ax.set_title("SV burden by clade and type")
    ax.set_xlabel("Clade / phylum")
    ax.set_ylabel("SV count")
    ax.legend(title="SV type", bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fname = "clade_sv_burden.png"
    encoded = write_figure(fig, outdir / fname)
    return FigureRecord(
        title="SV burden by clade",
        filename=fname,
        encoded_png=encoded,
        caption="Stacked SV counts per clade/phylum, coloured by SV class.",
    )


def plot_te_architecture(df: pd.DataFrame, outdir: Path) -> Optional[FigureRecord]:
    """TE class counts and median size for each TE family."""
    te_col = "te_class" if "te_class" in df.columns else None
    if te_col is None:
        return None
    te_df = df[df[te_col].notna() & (df[te_col] != "NONE")].copy()
    if te_df.empty:
        return None
    counts = te_df[te_col].value_counts().reset_index()
    counts.columns = [te_col, "count"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    axes[0].barh(counts[te_col].astype(str), counts["count"])
    axes[0].set_title("TE class counts")
    axes[0].set_xlabel("Count")
    axes[0].set_ylabel("TE class")
    if "sv_len" in te_df.columns:
        order = counts[te_col].tolist()
        data = [te_df.loc[te_df[te_col] == cls, "sv_len"].dropna().values for cls in order]
        boxplot_with_labels(axes[1], data, order, showfliers=False, vert=False)
        axes[1].set_title("TE size distribution (bp)")
        axes[1].set_xlabel("SV length (bp)")
    else:
        axes[1].axis("off")
    fig.tight_layout()
    fname = "te_architecture.png"
    encoded = write_figure(fig, outdir / fname)
    return FigureRecord(
        title="TE architecture",
        filename=fname,
        encoded_png=encoded,
        caption="TE class abundance (left) and size distribution by class (right).",
    )


def plot_hgt_propagation(df: pd.DataFrame, outdir: Path) -> Optional[FigureRecord]:
    """Novelty-tier distribution for HGT candidate loci across clades."""
    nt_col = "novelty_tier" if "novelty_tier" in df.columns else None
    if nt_col is None:
        return None
    hgt_df = df.copy()
    if nt_col not in hgt_df.columns:
        return None
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    tier_counts = hgt_df[nt_col].fillna("UNKNOWN").value_counts().reset_index()
    tier_counts.columns = [nt_col, "count"]
    axes[0].bar(tier_counts[nt_col].astype(str), tier_counts["count"])
    axes[0].set_title("HGT novelty tier distribution")
    axes[0].set_xlabel("Novelty tier")
    axes[0].set_ylabel("Locus count")
    axes[0].tick_params(axis="x", rotation=30)
    if "clade" in hgt_df.columns:
        novel = hgt_df[hgt_df[nt_col] == "NOVEL"]
        if not novel.empty:
            clade_counts = novel["clade"].fillna("unknown").value_counts().head(15).reset_index()
            clade_counts.columns = ["clade", "count"]
            axes[1].barh(clade_counts["clade"].astype(str), clade_counts["count"])
            axes[1].set_title("NOVEL loci by clade (top 15)")
            axes[1].set_xlabel("Count")
        else:
            axes[1].axis("off")
    else:
        axes[1].axis("off")
    fig.tight_layout()
    fname = "hgt_propagation.png"
    encoded = write_figure(fig, outdir / fname)
    return FigureRecord(
        title="HGT propagation",
        filename=fname,
        encoded_png=encoded,
        caption="Novelty tier distribution (left) and NOVEL loci per clade (right).",
    )


def _select_wins_subset(real_df: pd.DataFrame) -> pd.DataFrame:
    """Slice real_merged.tsv to the rows that drive the wins-matrix panel.

    We want per-(query, label, svtype, method) F1 scored against the most
    independent validation basis available. Prefer raw-read validated
    `read_level_union` rows. Fall back to read-supported comparator agreement,
    then to unvalidated comparator consensus for older/partial outputs.
    Deliberately exclude per-comparator labels such as minigraph_read_supported;
    those are "mycosv vs one comparator baseline" diagnostics, not the
    apples-to-apples wins view.
    """
    if real_df is None or real_df.empty:
        return pd.DataFrame()
    needed = {"truth_label", "method", "svtype", "f1"}
    if not needed.issubset(real_df.columns):
        return pd.DataFrame()
    out = real_df.copy()
    label = out["truth_label"].astype(str)
    keep = label.eq("read_level_union") | label.str.startswith("consensus_")
    if "truth_calls" in out.columns:
        keep &= pd.to_numeric(out["truth_calls"], errors="coerce").fillna(0) > 0
    if "status" in out.columns:
        keep &= ~out["status"].astype(str).eq("no_truth")
    out = out[keep].copy()
    label = out["truth_label"].astype(str)
    raw_label = label.eq("read_level_union")
    if raw_label.any():
        out = out[raw_label].copy()
    else:
        rs_label = label.str.contains("_read_supported", regex=False)
        if rs_label.any():
            out = out[rs_label].copy()
    out["f1"] = pd.to_numeric(out["f1"], errors="coerce")
    out = out.dropna(subset=["f1"])
    return out


def plot_wins_matrix(real_df: pd.DataFrame, outdir: Path) -> List[FigureRecord]:
    """Render the per-SV-type wins matrix: for each comparator-as-method row,
    is MycoSV's F1 ≥ that comparator's F1 on the same validation basis?

    Outputs three figure cards:
      A. F1 heatmap, rows=method (mycosv + each comparator), cols=svtype.
         Cell value = mean F1 across queries.
      B. Wins bar: per (svtype, comparator), fraction of queries where
         mycosv F1 ≥ comparator F1 on the same validation label.
      C. Per-query F1 scatter: mycosv vs each comparator, one panel per
         comparator, dot per (query × svtype).
    Prefer raw-read validated labels so we never promote a comparator-only
    baseline to biological truth.
    """
    figs: List[FigureRecord] = []
    subset = _select_wins_subset(real_df)
    if subset.empty:
        return figs
    rs = subset.copy()
    if rs.empty:
        return figs
    rs["svtype"] = rs["svtype"].astype(str)
    rs["method"] = rs["method"].astype(str)

    # ── (A) F1 heatmap, mean across queries ─────────────────────────────
    pivot = rs.pivot_table(
        index="method", columns="svtype", values="f1", aggfunc="mean",
    )
    if not pivot.empty:
        # Put mycosv at the top so the eye reads it as "the method we are
        # testing"; sort the rest alphabetically so panels are stable across
        # runs even when comparator availability changes.
        ordered = ["mycosv"] + sorted(m for m in pivot.index if m != "mycosv")
        pivot = pivot.reindex([m for m in ordered if m in pivot.index])
        sv_order = [s for s in ("ALL", "INS", "DEL", "INV", "DUP", "TRA", "OFF_REF") if s in pivot.columns]
        if sv_order:
            pivot = pivot[sv_order]
        fig, ax = plt.subplots(figsize=(max(6, 0.9 * len(pivot.columns) + 4), max(3, 0.5 * len(pivot.index) + 2)))
        masked = np.ma.masked_invalid(pivot.values.astype(float))
        cmap = plt.get_cmap("viridis").copy()
        cmap.set_bad(color="#f2f2f2")
        im = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=0.0, vmax=1.0)
        ax.set_xticks(range(len(pivot.columns)))
        ax.set_xticklabels(pivot.columns, rotation=0)
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels(pivot.index)
        for i in range(pivot.shape[0]):
            for j in range(pivot.shape[1]):
                v = float(pivot.values[i, j])
                label = "NA" if np.isnan(v) else f"{v:.2f}"
                ax.text(j, i, label, ha="center", va="center",
                        color=("white" if not np.isnan(v) and v < 0.55 else "black"), fontsize=9)
        ax.set_title("Mean F1 by validation basis (method × SV type)")
        fig.colorbar(im, ax=ax, label="F1")
        fig.tight_layout()
        encoded = write_figure(fig, outdir / "wins_f1_heatmap.png")
        save_df(pivot.reset_index(), outdir / "wins_f1_heatmap.tsv")
        figs.append(FigureRecord(
            title="Mean F1 by method × SV type",
            filename="wins_f1_heatmap.png",
            encoded_png=encoded,
            caption=(
                "Each cell is the mean F1 across queries when scoring a "
                "method against the strongest available validation basis for "
                "that SV type: raw-read validated rows first, then "
                "read-supported comparator agreement, then comparator "
                "consensus for older outputs. mycosv row is pinned to the top. "
                "wins_f1_heatmap.tsv carries the underlying numbers."
            ),
        ))

    # ── (B) Wins bar: % queries where mycosv ≥ comparator ───────────────
    if "query_asm" in rs.columns:
        wide = rs.pivot_table(
            index=["query_asm", "truth_label", "svtype"], columns="method",
            values="f1", aggfunc="mean",
        )
        if "mycosv" in wide.columns:
            wins_rows = []
            for cmp in [c for c in wide.columns if c != "mycosv"]:
                pair = wide[["mycosv", cmp]].dropna()
                if pair.empty:
                    continue
                for sv in sorted(pair.index.get_level_values("svtype").unique()):
                    sub = pair.xs(sv, level="svtype", drop_level=False)
                    if sub.empty:
                        continue
                    n = len(sub)
                    wins = int((sub["mycosv"] >= sub[cmp]).sum())
                    ties = int((sub["mycosv"] == sub[cmp]).sum())
                    zero_ties = int(((sub["mycosv"] == 0) & (sub[cmp] == 0)).sum())
                    wins_rows.append({
                        "comparator": cmp,
                        "svtype": sv,
                        "queries": n,
                        "mycosv_wins": wins,
                        "mycosv_ties": ties,
                        "zero_ties": zero_ties,
                        "win_rate": wins / n if n else 0.0,
                    })
            if wins_rows:
                wins_df = pd.DataFrame(wins_rows)
                pivot_wins = wins_df.pivot_table(
                    index="svtype", columns="comparator", values="win_rate",
                )
                sv_order = [s for s in ("ALL", "INS", "DEL", "INV", "DUP", "TRA", "OFF_REF") if s in pivot_wins.index]
                if sv_order:
                    pivot_wins = pivot_wins.reindex(sv_order)
                fig, ax = plt.subplots(figsize=(max(7, 0.85 * len(pivot_wins.columns) + 3), max(3.5, 0.55 * len(pivot_wins.index) + 2.2)))
                masked = np.ma.masked_invalid(pivot_wins.values.astype(float))
                cmap = plt.get_cmap("magma").copy()
                cmap.set_bad(color="#f2f2f2")
                im = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=0.0, vmax=1.0)
                ax.set_xticks(range(len(pivot_wins.columns)))
                ax.set_xticklabels(pivot_wins.columns, rotation=35, ha="right")
                ax.set_yticks(range(len(pivot_wins.index)))
                ax.set_yticklabels(pivot_wins.index)
                for i in range(pivot_wins.shape[0]):
                    for j in range(pivot_wins.shape[1]):
                        v = float(pivot_wins.values[i, j])
                        label = "NA" if np.isnan(v) else f"{v:.2f}"
                        ax.text(j, i, label, ha="center", va="center",
                                color=("white" if not np.isnan(v) and v < 0.65 else "black"), fontsize=9)
                ax.set_xlabel("Comparator")
                ax.set_ylabel("SV type")
                ax.set_title("Per-(SV type, comparator) win/tie rate of MycoSV")
                fig.colorbar(im, ax=ax, label="MycoSV win/tie rate (F1 >= comparator)")
                fig.tight_layout()
                encoded = write_figure(fig, outdir / "wins_rate.png")
                save_df(wins_df, outdir / "wins_rate.tsv")
                figs.append(FigureRecord(
                    title="MycoSV win rate vs each comparator",
                    filename="wins_rate.png",
                    encoded_png=encoded,
                    caption=(
                        "Per (SV type, comparator), the fraction of comparable "
                        "query/truth-label cells where MycoSV F1 ≥ that "
                        "comparator's F1 on the same validation label. "
                        "Ties, including 0-vs-0 ties with non-empty truth, "
                        "count because the panel is a win/tie rate. Grey "
                        "cells are unavailable comparisons. wins_rate.tsv "
                        "also reports tie counts."
                    ),
                ))

            # ── (C) Per-query F1 scatter: mycosv vs each comparator ────
            cmps = [c for c in wide.columns if c != "mycosv"]
            if cmps:
                ncols = min(3, len(cmps))
                nrows = math.ceil(len(cmps) / ncols)
                fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows), squeeze=False)
                for idx, cmp in enumerate(cmps):
                    ax = axes[idx // ncols][idx % ncols]
                    pair = wide[["mycosv", cmp]].dropna()
                    pair = pair[(pair["mycosv"] > 0) | (pair[cmp] > 0)].reset_index()
                    if pair.empty:
                        ax.axis("off")
                        continue
                    for sv, group in pair.groupby("svtype"):
                        ax.scatter(group[cmp], group["mycosv"], label=sv, alpha=0.7)
                    ax.plot([0, 1], [0, 1], color="grey", linestyle="--", linewidth=0.8)
                    ax.set_xlim(0, 1.02); ax.set_ylim(0, 1.02)
                    ax.set_xlabel(f"{cmp} F1")
                    ax.set_ylabel("mycosv F1")
                    ax.set_title(f"mycosv vs {cmp}")
                    ax.legend(loc="lower right", fontsize=7)
                # Hide unused panes.
                for j in range(len(cmps), nrows * ncols):
                    axes[j // ncols][j % ncols].axis("off")
                fig.tight_layout()
                encoded = write_figure(fig, outdir / "wins_scatter.png")
                figs.append(FigureRecord(
                    title="Per-query F1 scatter: MycoSV vs each comparator",
                    filename="wins_scatter.png",
                    encoded_png=encoded,
                    caption=(
                        "Each panel: dot per (query × SV type), x-axis is "
                        "the comparator F1, y-axis is MycoSV F1, grey "
                        "diagonal is parity. Dots above the diagonal are "
                        "MycoSV wins."
                    ),
                ))
    return figs


def build_wins_matrix_section(
    real_df: Optional[pd.DataFrame],
    outdir: Path,
) -> Tuple[str, List[FigureRecord], List[pd.DataFrame]]:
    if real_df is None or real_df.empty:
        return ("<p>No real-data PR table provided — wins matrix skipped.</p>", [], [])
    if not plotting_available():
        return (plotting_unavailable_html(), [], [])
    df = real_df.copy()
    # Don't run real_df through harmonize_columns: the wins matrix needs the
    # raw `truth_label`, `svtype`, `method`, `f1`, and `query_asm` columns
    # written by run_real_fungal_benchmark.write_agreement_summary, and
    # harmonize_columns aliases `svtype` → `sv_type` which would break the
    # downstream filter.
    figs = plot_wins_matrix(df, outdir)
    if not figs:
        return ("<p>No raw-read or comparator-agreement PR rows found in the real-data table — "
                "wins matrix skipped.</p>", [], [])
    html = (
        "<p>Apples-to-apples comparison: every method (MycoSV + each "
        "algorithmic comparator) is scored as a predictor against the same "
        "validation basis. The report prefers <code>read_level_union</code> "
        "rows, which are independently anchored in raw FASTQ/read evidence. "
        "Comparator consensus rows are retained as baseline agreement only, "
        "not treated as fungal truth.</p>"
    )
    return html, figs, []


def build_clade_te_hgt_section(
    real_df: Optional[pd.DataFrame],
    bio_df: Optional[pd.DataFrame],
    outdir: Path,
) -> Tuple[str, List[FigureRecord], List[pd.DataFrame]]:
    figs: List[FigureRecord] = []
    tables: List[pd.DataFrame] = []
    if not plotting_available():
        return plotting_unavailable_html(), [], []

    for src_df in [real_df, bio_df]:
        if src_df is None or src_df.empty:
            continue
        df = harmonize_columns(src_df)

        rec = plot_clade_sv(df, outdir)
        if rec is not None:
            figs.append(rec)

        rec = plot_te_architecture(df, outdir)
        if rec is not None:
            figs.append(rec)

        rec = plot_hgt_propagation(df, outdir)
        if rec is not None:
            figs.append(rec)

        if figs:
            break  # plots generated from first non-empty source

    if not figs:
        return "<p>No clade/TE/HGT columns found in provided tables.</p>", [], []

    html = "<p>Clade-stratified SV landscape, TE family architecture, and HGT novelty propagation across fungal lineages.</p>"
    return html, figs, tables


# -----------------------------
# Novel-SV biological-question section
#
# Three questions, scored over MycoSV-only (mycosv_unique=yes) calls joined
# to the biology candidate annotations:
#   Q1. Which novel SVs are HGT/Starship cargo crossing clade boundaries?
#   Q2. Which novel SVs sit in two-speed accessory / TE-rich architecture?
#   Q3. Which novel SVs have direct expression evidence at a nearby gene?
# -----------------------------


_TE_ELEMENT_CLASSES = {
    "TE", "TE_LTR", "TE_TIR", "TE_LINE", "TE_SINE",
    "LTR_GYPSY", "LTR_COPIA", "LINE", "SINE",
    "DNA_TIR", "HELITRON", "MITE", "RIP", "REPEAT",
}
_HGT_ELEMENT_CLASSES = {"HGT", "STARSHIP"}


def _join_novel_to_biology(
    novel_df: pd.DataFrame,
    bio_df: Optional[pd.DataFrame],
) -> pd.DataFrame:
    if novel_df is None or novel_df.empty:
        return pd.DataFrame()
    novel_only = novel_df.copy()
    if "mycosv_unique" in novel_only.columns:
        novel_only = novel_only[novel_only["mycosv_unique"].astype(str).str.lower() == "yes"]
    if bio_df is None or bio_df.empty:
        return novel_only.copy()
    bio = bio_df.copy()
    # harmonize_columns() may have renamed the join keys: query_asm → sample,
    # query_contig → chrom, pos → start, svtype → sv_type. Try both the
    # original and the harmonized names so the join works in both code paths.
    candidate_keys = [
        "query_asm", "sample",
        "query_contig", "chrom",
        "pos", "start",
        "end",
        "svtype", "sv_type",
    ]
    join_cols = [
        c for c in candidate_keys
        if c in novel_only.columns and c in bio.columns
    ]
    # De-dupe while preserving order (preserves "sample, chrom, start, end, sv_type"
    # ordering when both column-name conventions match the same logical key).
    seen: set = set()
    join_cols = [c for c in join_cols if not (c in seen or seen.add(c))]
    if not join_cols:
        return novel_only
    # Coerce numeric join keys to consistent types so the merge is order-stable.
    for col in ("pos", "start", "end"):
        if col in novel_only.columns:
            novel_only[col] = pd.to_numeric(novel_only[col], errors="coerce").astype("Int64")
        if col in bio.columns:
            bio[col] = pd.to_numeric(bio[col], errors="coerce").astype("Int64")
    return novel_only.merge(bio, on=join_cols, how="left", suffixes=("", "_bio"))


def _novel_q1_hgt(joined: pd.DataFrame, outdir: Path) -> Optional[FigureRecord]:
    """Q1: novel SVs with HGT/Starship cargo across clade or phylum boundaries."""
    if joined.empty:
        return None
    # element_class is aliased to "effect" by harmonize_columns; check both.
    ec_col = next(
        (c for c in ("element_class", "effect", "element_class_bio", "effect_bio") if c in joined.columns),
        None,
    )
    if ec_col is None:
        return None
    df = joined[joined[ec_col].astype(str).isin(_HGT_ELEMENT_CLASSES)].copy()
    # Translocations / off-reference also surface HGT-style breakpoints when
    # the element_class is missing on the row but the call-type indicates a
    # cross-locus event flagged NOVEL by the routing layer. svtype and
    # annotation are aliased to sv_type and pathway by harmonize_columns().
    sv_col = "sv_type" if "sv_type" in joined.columns else ("svtype" if "svtype" in joined.columns else None)
    annot_col = "pathway" if "pathway" in joined.columns else ("annotation" if "annotation" in joined.columns else None)
    if sv_col is not None and annot_col is not None:
        annot_lower = joined[annot_col].astype(str).str.upper()
        df_extra = joined[
            joined[sv_col].astype(str).isin({"TRA", "OFF_REF"})
            & annot_lower.isin({"NOVEL", "NOVEL_WEAK"})
        ]
        df = pd.concat([df, df_extra], ignore_index=True).drop_duplicates()
    if df.empty:
        return None
    # phylum is aliased to "clade" by harmonize_columns; pick whichever exists.
    phylum_col = next(
        (c for c in ("phylum", "clade", "query_asm") if c in df.columns),
        None,
    )
    by_phylum = (
        df[phylum_col].fillna("unknown").astype(str).value_counts().head(15)
        if phylum_col is not None
        else pd.Series(dtype="int64")
    )
    if by_phylum.empty:
        return None
    fig, ax = plt.subplots(figsize=(10, 5))
    by_phylum.plot(kind="bar", ax=ax, color="#b04a3a")
    ax.set_title("Q1. Novel HGT / Starship-cargo SV candidates per clade")
    ax.set_xlabel("Clade / phylum / query_asm")
    ax.set_ylabel("Novel HGT-class SV count")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    encoded = write_figure(fig, outdir / "novel_q1_hgt_cargo.png")
    df.to_csv(outdir / "novel_q1_hgt_cargo.tsv", sep="\t", index=False)
    return FigureRecord(
        title="Q1. HGT / Starship cargo (novel)",
        filename="novel_q1_hgt_cargo.png",
        encoded_png=encoded,
        caption=(
            "Novel MycoSV-only SVs whose element_class is HGT or STARSHIP, "
            "or cross-locus TRA/OFF_REF flagged NOVEL by the routing layer. "
            "Bars show count per clade. Per-row table: "
            "novel_q1_hgt_cargo.tsv."
        ),
    )


def _novel_q2_two_speed(joined: pd.DataFrame, outdir: Path) -> Optional[FigureRecord]:
    """Q2: novel SVs in two-speed / TE-rich accessory architecture."""
    if joined.empty:
        return None
    ec_col = next(
        (c for c in ("element_class", "effect", "element_class_bio", "effect_bio") if c in joined.columns),
        None,
    )
    if ec_col is None:
        return None
    df = joined[joined[ec_col].astype(str).isin(_TE_ELEMENT_CLASSES)].copy()
    if df.empty:
        return None
    # Stratify by architecture (two_speed, te_rich, smut_pathogen, …) when
    # the metadata column is present; otherwise fall back to scenario or to
    # the harmonized "phenotype" column (architecture/scenario/lifestyle all
    # map to phenotype via COLUMN_ALIASES).
    arch_col = next(
        (c for c in ("architecture", "scenario", "lifestyle", "phenotype") if c in df.columns),
        None,
    )
    if arch_col is None:
        return None
    by_arch = (
        df[arch_col]
        .fillna("unknown")
        .astype(str)
        .value_counts()
        .head(12)
    )
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    by_arch.plot(kind="bar", ax=axes[0], color="#3a6db0")
    axes[0].set_title(f"Q2. Novel TE-class SVs by {arch_col}")
    axes[0].set_xlabel(arch_col)
    axes[0].set_ylabel("Novel TE-class SV count")
    axes[0].tick_params(axis="x", rotation=45)
    by_class = df[ec_col].fillna("NONE").astype(str).value_counts().head(12)
    by_class.plot(kind="barh", ax=axes[1], color="#3a6db0")
    axes[1].set_title("Q2. TE element-class composition")
    axes[1].set_xlabel("Count")
    axes[1].set_ylabel("element_class")
    fig.tight_layout()
    encoded = write_figure(fig, outdir / "novel_q2_two_speed.png")
    df.to_csv(outdir / "novel_q2_two_speed.tsv", sep="\t", index=False)
    return FigureRecord(
        title="Q2. Two-speed / TE-rich accessory architecture (novel)",
        filename="novel_q2_two_speed.png",
        encoded_png=encoded,
        caption=(
            "Novel MycoSV-only SVs sitting in TE-class sequence (LTRs, DNA "
            "transposons, helitrons, repeat-rich regions, RIP). Left: count "
            "per architecture/scenario, exposing two-speed / TE-rich "
            "compartments. Right: TE element_class composition. Per-row "
            "table: novel_q2_two_speed.tsv."
        ),
    )


def _novel_q3_expression(joined: pd.DataFrame, outdir: Path) -> Optional[FigureRecord]:
    """Q3: novel SVs with direct gene-expression support at a nearby gene."""
    if joined.empty:
        return None
    es_col = next((c for c in ("expression_supported", "expression_supported_bio") if c in joined.columns), None)
    if es_col is None:
        return None
    df = joined[joined[es_col].astype(str).str.lower() == "yes"].copy()
    if df.empty:
        return None
    # Plot the log2_fc / -log10(padj) "volcano" of the expression-supported
    # nearby genes for these novel SVs, plus a top-gene bar.
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    log2_col = next((c for c in ("expression_log2_fc", "log2_fc", "expression") if c in df.columns), None)
    padj_col = next((c for c in ("expression_padj", "padj") if c in df.columns), None)
    if log2_col and padj_col:
        x = pd.to_numeric(df[log2_col], errors="coerce")
        p = pd.to_numeric(df[padj_col], errors="coerce").clip(lower=1e-300)
        y = -np.log10(p)
        axes[0].scatter(x, y, alpha=0.7, color="#2e8b57")
        axes[0].axhline(-math.log10(0.05), color="grey", linestyle="--", linewidth=0.8)
        axes[0].set_xlabel(f"{log2_col}")
        axes[0].set_ylabel(f"-log10({padj_col})")
        axes[0].set_title("Q3. Volcano of expression-supported novel SVs")
    else:
        axes[0].axis("off")
    gene_col = next((c for c in ("expression_gene", "gene") if c in df.columns), None)
    if gene_col is not None:
        top = df[gene_col].fillna("unknown").astype(str).value_counts().head(15)
        top.plot(kind="barh", ax=axes[1], color="#2e8b57")
        axes[1].set_title("Q3. Top genes near novel expression-supported SVs")
        axes[1].set_xlabel("Novel SV count")
        axes[1].set_ylabel("Gene")
    else:
        axes[1].axis("off")
    fig.tight_layout()
    encoded = write_figure(fig, outdir / "novel_q3_expression.png")
    df.to_csv(outdir / "novel_q3_expression.tsv", sep="\t", index=False)
    return FigureRecord(
        title="Q3. Direct gene-expression-supported novel SVs",
        filename="novel_q3_expression.png",
        encoded_png=encoded,
        caption=(
            "Novel MycoSV-only SVs where the analyzer flagged a "
            "differentially-expressed gene within the expression window. "
            "Left: volcano of nearby genes. Right: top recurrent genes "
            "across novel events. Per-row table: novel_q3_expression.tsv."
        ),
    )


def build_evidence_panorama_section(
    tiers_df: Optional[pd.DataFrame],
    outdir: Path,
) -> Tuple[str, List[FigureRecord], List[pd.DataFrame]]:
    """Per-query stacked-bar panorama of MycoSV calls by evidence tier.

    Reads `mycosv_evidence_tiers.tsv` (query_asm, svtype, tier, n_calls) and
    plots one stacked bar per query (ALL svtypes summed) plus a per-svtype
    grid. Tier ordering is fixed at strong > moderate > intrinsic_only > weak
    so the legend and color scheme are stable across panels. The point of this
    panel is to surface MycoSV calls that are real-but-unvalidatable
    (intrinsic_only) — these would otherwise be invisible in the F1 panels
    because validation drops them from the independently validated set.
    """
    if tiers_df is None or tiers_df.empty:
        return ("<p>No evidence-tier table provided "
                "(--evidence-tiers mycosv_evidence_tiers.tsv).</p>", [], [])
    if not plotting_available():
        return (plotting_unavailable_html(), [], [])

    df = tiers_df.copy()
    needed = {"query_asm", "svtype", "tier", "n_calls"}
    if not needed.issubset(df.columns):
        missing = ", ".join(sorted(needed - set(df.columns)))
        return (f"<p>evidence-tier table is missing columns: {missing}</p>", [], [])
    df["n_calls"] = pd.to_numeric(df["n_calls"], errors="coerce").fillna(0).astype(int)
    df = df[df["n_calls"] > 0]
    if df.empty:
        return ("<p>No MycoSV calls to tier in this run.</p>", [], [])

    tier_order = ["strong", "moderate", "intrinsic_only", "weak"]
    tier_colors = {
        "strong":         "#1b7837",  # comparator + read evidence
        "moderate":       "#7fbf7b",  # one of the two
        "intrinsic_only": "#f1a340",  # MycoSV cluster only — the real-but-unvalidatable bucket
        "weak":           "#998ec3",  # very low intrinsic, rare
    }

    figs: List[FigureRecord] = []
    # ── (A) per query, all SV types summed ───────────────────────────────
    by_query = df.groupby(["query_asm", "tier"], as_index=False)["n_calls"].sum()
    pivot_q = by_query.pivot(index="query_asm", columns="tier", values="n_calls").fillna(0)
    for tier in tier_order:
        if tier not in pivot_q.columns:
            pivot_q[tier] = 0
    pivot_q = pivot_q[tier_order]
    fig, ax = plt.subplots(figsize=(max(7, 0.6 * len(pivot_q) + 4), 4.5))
    pivot_q.plot(kind="bar", stacked=True, ax=ax,
                 color=[tier_colors[t] for t in tier_order])
    ax.set_title("MycoSV call panorama by evidence tier (per query)")
    ax.set_xlabel("query")
    ax.set_ylabel("MycoSV calls")
    ax.legend(title="evidence", loc="upper right")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    fig.tight_layout()
    encoded = write_figure(fig, outdir / "mycosv_evidence_panorama.png")
    save_df(pivot_q.reset_index(), outdir / "mycosv_evidence_panorama.tsv")
    figs.append(FigureRecord(
        title="MycoSV call panorama — strong / moderate / intrinsic_only / weak",
        filename="mycosv_evidence_panorama.png",
        encoded_png=encoded,
        caption=(
            "Stacked counts of every MycoSV call in this run, tiered by the "
            "evidence backing it. <b>strong</b> = matched by ≥1 comparator AND "
            "read-validated. <b>moderate</b> = matched by one of the two. "
            "<b>intrinsic_only</b> = no external evidence, but the MycoSV "
            "C++ caller clustered the call at SUPPORT≥2 (sibling-clade "
            "contig absent from the validation BAM, low query coverage, …). "
            "These are real-but-unvalidatable calls that the F1 plots "
            "otherwise penalize as FP. <b>weak</b> = residual intrinsic≤1. "
            "Source: mycosv_evidence_tiers.tsv."
        ),
    ))

    # ── (B) per (query × svtype) grid ────────────────────────────────────
    by_q_svt = df.groupby(["query_asm", "svtype", "tier"], as_index=False)["n_calls"].sum()
    pivot_sv = by_q_svt.pivot_table(
        index=["query_asm", "svtype"], columns="tier",
        values="n_calls", aggfunc="sum",
    ).fillna(0)
    for tier in tier_order:
        if tier not in pivot_sv.columns:
            pivot_sv[tier] = 0
    pivot_sv = pivot_sv[tier_order]
    # Render each query as a row of svtype subplots so the reader can compare
    # tier mix across SV classes within the same query.
    queries = sorted(df["query_asm"].unique())
    sv_order = [s for s in ("INS", "DEL", "DUP", "INV", "TRA", "OFF_REF") if s in df["svtype"].unique()]
    if queries and sv_order:
        fig2, axes = plt.subplots(
            len(queries), 1,
            figsize=(max(7, 0.8 * len(sv_order) + 3), max(3, 2.3 * len(queries))),
            sharex=True, squeeze=False,
        )
        for i, qa in enumerate(queries):
            ax = axes[i, 0]
            sub = pivot_sv.loc[qa] if qa in pivot_sv.index.get_level_values(0) else pd.DataFrame()
            if sub.empty:
                ax.set_visible(False)
                continue
            sub = sub.reindex(sv_order).fillna(0)[tier_order]
            sub.plot(kind="bar", stacked=True, ax=ax,
                     color=[tier_colors[t] for t in tier_order], legend=(i == 0))
            ax.set_title(short_sample_label(qa))
            ax.set_xlabel("")
            ax.set_ylabel("calls")
            plt.setp(ax.get_xticklabels(), rotation=0)
        fig2.suptitle("MycoSV evidence tier × SV type (per query)")
        fig2.tight_layout()
        encoded2 = write_figure(fig2, outdir / "mycosv_evidence_by_svtype.png")
        save_df(pivot_sv.reset_index(), outdir / "mycosv_evidence_by_svtype.tsv")
        figs.append(FigureRecord(
            title="MycoSV evidence tier × SV type (per query)",
            filename="mycosv_evidence_by_svtype.png",
            encoded_png=encoded2,
            caption=(
                "Per-SV-type breakdown of the panorama. Tall <b>intrinsic_only</b> "
                "bars on a given SV type usually indicate one of: (a) sibling-clade "
                "REFCONTIG that the validation BAM never indexed, (b) DUP/INV/TRA "
                "lacking a clean split-read signal in the assembly-mode "
                "contig-vs-ref BAM, or (c) low query coverage. None of these "
                "imply the call is wrong; they just mean the bench's external "
                "scaffolding can't reach it."
            ),
        ))

    body = (
        "<p>This panel makes <b>real-but-unvalidatable</b> MycoSV calls "
        "first-class. F1 plots above use a validation set that has been shrunk "
        "by read-validation; calls in the <b>intrinsic_only</b> tier here are "
        "MycoSV predictions that the external pipeline could not confirm "
        "even though MycoSV's own clustering passed (SUPPORT≥2). Counts come "
        "from <code>mycosv_evidence_tiers.tsv</code>; per-call tiers are in "
        "<code>novel_mycosv_calls.tsv</code> and <code>biology_findings.tsv</code>.</p>"
    )
    return body, figs, [pivot_q.reset_index(), pivot_sv.reset_index()]


def build_novel_questions_section(
    novel_df: Optional[pd.DataFrame],
    bio_df: Optional[pd.DataFrame],
    outdir: Path,
) -> Tuple[str, List[FigureRecord], List[pd.DataFrame]]:
    if novel_df is None or novel_df.empty:
        return "<p>No novel-SV table provided (--novel novel_mycosv_calls.tsv).</p>", [], []
    novel_df = harmonize_columns(novel_df)
    bio_df_h = harmonize_columns(bio_df) if bio_df is not None else None
    joined = _join_novel_to_biology(novel_df, bio_df_h)
    if joined.empty:
        return "<p>No MycoSV-unique novel SVs to highlight.</p>", [], []

    save_df(joined, outdir / "novel_questions_joined.tsv")
    if not plotting_available():
        return (
            "<p>Novel-SV biological-question table was generated, but plots "
            "were skipped because matplotlib is not installed.</p>",
            [],
            [],
        )

    figs: List[FigureRecord] = []
    for builder in (_novel_q1_hgt, _novel_q2_two_speed, _novel_q3_expression):
        try:
            rec = builder(joined, outdir)
        except Exception as exc:
            sys.stderr.write(
                f"[novel-questions] {builder.__name__} failed: "
                f"{type(exc).__name__}: {exc}\n"
            )
            rec = None
        if rec is not None:
            figs.append(rec)

    html_intro = (
        "<p>Three biological questions over the MycoSV-unique (novel) SV set, "
        "joined to <code>biology_findings.tsv</code> for element_class, "
        "phylum / scenario / architecture metadata, and expression evidence:</p>"
        "<ol>"
        "<li><b>Q1.</b> Which novel SVs are HGT / Starship cargo crossing clade boundaries?</li>"
        "<li><b>Q2.</b> Which novel SVs sit in two-speed accessory / TE-rich architecture?</li>"
        "<li><b>Q3.</b> Which novel SVs have direct nearby-gene expression support?</li>"
        "</ol>"
    )
    if not figs:
        html_intro += "<p><i>No rows matched any of the three questions; nothing to plot.</i></p>"
    return html_intro, figs, []


# -----------------------------
# HTML rendering
# -----------------------------


def render_figure_cards(figs: List[FigureRecord]) -> str:
    chunks = []
    for fig in figs:
        chunks.append(
            f"""
            <div class=\"card\">
              <h3>{fig.title}</h3>
              <img src=\"data:image/png;base64,{fig.encoded_png}\" alt=\"{fig.title}\" />
              <p class=\"caption\">{fig.caption}</p>
              <p class=\"small\">Saved file: <code>{fig.filename}</code></p>
            </div>
            """
        )
    return "\n".join(chunks)



def summary_stats(sim_df: Optional[pd.DataFrame], real_df: Optional[pd.DataFrame], bio_df: Optional[pd.DataFrame]) -> Dict[str, str]:
    stats: Dict[str, str] = {}
    if sim_df is not None and not sim_df.empty:
        sdf = harmonize_columns(sim_df)
        stats["Simulated rows"] = f"{len(sdf):,}"
        if "caller" in sdf.columns:
            stats["Simulated callers"] = str(sdf["caller"].nunique(dropna=True))
        if "f1" in sdf.columns:
            stats["Mean simulated F1"] = f"{sdf['f1'].mean():.3f}"
    if real_df is not None and not real_df.empty:
        rdf = harmonize_columns(real_df)
        stats["Real SV calls"] = f"{len(rdf):,}"
        if "sample" in rdf.columns:
            stats["Real samples"] = str(rdf["sample"].nunique(dropna=True))
        if "sv_type" in rdf.columns:
            stats["SV classes observed"] = str(rdf["sv_type"].nunique(dropna=True))
    if bio_df is not None and not bio_df.empty:
        bdf = harmonize_columns(bio_df)
        stats["Biology rows"] = f"{len(bdf):,}"
        if "gene" in bdf.columns:
            stats["Unique genes"] = str(bdf["gene"].nunique(dropna=True))
        if "pathway" in bdf.columns:
            stats["Unique pathways"] = str(bdf["pathway"].nunique(dropna=True))
    return stats



def render_summary_tiles(stats: Dict[str, str]) -> str:
    tiles = []
    for k, v in stats.items():
        tiles.append(f"<div class=\"tile\"><div class=\"tile-value\">{v}</div><div class=\"tile-label\">{k}</div></div>")
    return "\n".join(tiles)



def build_html_report(
    title: str,
    stats: Dict[str, str],
    sim_html: str,
    sim_figs: List[FigureRecord],
    real_html: str,
    real_figs: List[FigureRecord],
    bio_html: str,
    bio_figs: List[FigureRecord],
    clade_html: str = "",
    clade_figs: Optional[List[FigureRecord]] = None,
    novel_html: str = "",
    novel_figs: Optional[List[FigureRecord]] = None,
    wins_html: str = "",
    wins_figs: Optional[List[FigureRecord]] = None,
    panorama_html: str = "",
    panorama_figs: Optional[List[FigureRecord]] = None,
) -> str:
    if clade_figs is None:
        clade_figs = []
    if novel_figs is None:
        novel_figs = []
    if wins_figs is None:
        wins_figs = []
    if panorama_figs is None:
        panorama_figs = []
    return f"""
<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\" />
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\" />
  <title>{title}</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 0; padding: 0; background: #fafafa; color: #222; }}
    .container {{ max-width: 1300px; margin: 0 auto; padding: 28px; }}
    .hero {{ padding: 24px 0 10px 0; border-bottom: 1px solid #ddd; margin-bottom: 18px; }}
    h1, h2, h3 {{ margin-top: 0; }}
    .tiles {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 12px; margin: 20px 0; }}
    .tile {{ background: white; border: 1px solid #e3e3e3; border-radius: 14px; padding: 16px; }}
    .tile-value {{ font-size: 1.7rem; font-weight: 700; }}
    .tile-label {{ font-size: 0.92rem; color: #555; margin-top: 6px; }}
    .section {{ margin: 28px 0; }}
    .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap: 18px; }}
    .card {{ background: white; border: 1px solid #e3e3e3; border-radius: 16px; padding: 16px; box-shadow: 0 2px 8px rgba(0,0,0,0.03); }}
    img {{ width: 100%; height: auto; border-radius: 10px; }}
    table {{ border-collapse: collapse; width: 100%; background: white; }}
    th, td {{ border-bottom: 1px solid #ececec; padding: 8px 10px; text-align: left; font-size: 0.92rem; }}
    .caption {{ color: #444; }}
    .small {{ color: #666; font-size: 0.88rem; }}
    code {{ background: #f1f1f1; padding: 2px 6px; border-radius: 6px; }}
  </style>
</head>
<body>
  <div class=\"container\">
    <div class=\"hero\">
      <h1>{title}</h1>
      <p>Automated visualization report for simulated benchmarking, real-data SV analyses, and biological interpretation.</p>
    </div>

    <div class=\"tiles\">{render_summary_tiles(stats)}</div>

    <div class=\"section\">
      <h2>1. Simulated data benchmarking</h2>
      {sim_html}
      <div class=\"grid\">{render_figure_cards(sim_figs)}</div>
    </div>

    <div class=\"section\">
      <h2>2. Real data structural variant analyses</h2>
      {real_html}
      <div class=\"grid\">{render_figure_cards(real_figs)}</div>
    </div>

    <div class=\"section\">
      <h2>2b. MycoSV vs comparators per SV type &mdash; validation-basis wins matrix</h2>
      {wins_html}
      <div class=\"grid\">{render_figure_cards(wins_figs)}</div>
    </div>

    <div class=\"section\">
      <h2>2c. MycoSV call panorama &mdash; evidence tiers</h2>
      {panorama_html}
      <div class=\"grid\">{render_figure_cards(panorama_figs)}</div>
    </div>

    <div class=\"section\">
      <h2>3. Biological findings</h2>
      {bio_html}
      <div class=\"grid\">{render_figure_cards(bio_figs)}</div>
    </div>

    <div class=\"section\">
      <h2>4. Clade-SV landscape, TE architecture &amp; HGT propagation</h2>
      {clade_html}
      <div class=\"grid\">{render_figure_cards(clade_figs)}</div>
    </div>

    <div class=\"section\">
      <h2>5. Novel SVs &mdash; three biological questions</h2>
      {novel_html}
      <div class=\"grid\">{render_figure_cards(novel_figs)}</div>
    </div>
  </div>
</body>
</html>
"""


# -----------------------------
# Main
# -----------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate an SV visualization report from simulated and real datasets.")
    parser.add_argument("--simulated", type=str, default=None, help="Simulated benchmark CSV/TSV")
    parser.add_argument("--real", type=str, default=None, help="Real-data SV calls CSV/TSV")
    parser.add_argument("--biology", type=str, default=None, help="Biological findings CSV/TSV")
    parser.add_argument("--metadata", type=str, default=None, help="Optional sample metadata CSV/TSV")
    parser.add_argument("--novel", type=str, default=None,
                        help="novel_mycosv_calls.tsv from a benchmark run; "
                             "powers the three biological-question section.")
    parser.add_argument("--evidence-tiers", type=str, default=None,
                        help="mycosv_evidence_tiers.tsv from a benchmark run; "
                             "powers the evidence-tier panorama panel that "
                             "surfaces real-but-unvalidatable MycoSV calls.")
    parser.add_argument("--outdir", type=str, required=True, help="Output directory")
    parser.add_argument("--title", type=str, default="SV visualization report", help="Report title")
    return parser.parse_args()



def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    ensure_dir(outdir)

    sim_df = read_table(args.simulated)
    real_df = read_table(args.real)
    bio_df = read_table(args.biology)
    metadata_df = read_table(args.metadata)
    novel_df = read_table(args.novel)
    tiers_df = read_table(args.evidence_tiers)

    sim_html, sim_figs, _ = build_simulated_section(sim_df, outdir)
    real_html, real_figs, _ = build_real_section(real_df, metadata_df, outdir)
    bio_html, bio_figs, _ = build_biology_section(bio_df, outdir)
    clade_html, clade_figs, _ = build_clade_te_hgt_section(real_df, bio_df, outdir)
    novel_html, novel_figs, _ = build_novel_questions_section(novel_df, bio_df, outdir)
    panorama_html, panorama_figs, _ = build_evidence_panorama_section(tiers_df, outdir)
    # Wins matrix runs on the *raw* real_df (pre-harmonization) so the
    # `truth_label`, `svtype`, `method` columns survive.
    wins_html, wins_figs, _ = build_wins_matrix_section(real_df, outdir)

    stats = summary_stats(sim_df, real_df, bio_df)
    html = build_html_report(
        args.title, stats,
        sim_html, sim_figs,
        real_html, real_figs,
        bio_html, bio_figs,
        clade_html, clade_figs,
        novel_html, novel_figs,
        wins_html, wins_figs,
        panorama_html=panorama_html, panorama_figs=panorama_figs,
    )

    report_path = outdir / "sv_visualization_report.html"
    report_path.write_text(html, encoding="utf-8")

    print(f"Report written to: {report_path}")
    print(f"Figures written under: {outdir}")


if __name__ == "__main__":
    main()
