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
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# -----------------------------
# Utilities
# -----------------------------


def read_table(path: Optional[str]) -> Optional[pd.DataFrame]:
    if path is None:
        return None
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Input file not found: {p}")
    sep = "\t" if p.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(p, sep=sep)


COLUMN_ALIASES: Dict[str, Sequence[str]] = {
    "sample": ["sample", "sample_id", "tumor_sample", "specimen", "case", "query_asm"],
    "dataset": ["dataset", "dataset_name", "cohort", "group"],
    "caller": ["caller", "tool", "method", "sv_caller", "alignment_mode"],
    "sv_type": ["sv_type", "type", "svclass", "sv_class", "event_type", "svtype"],
    "chrom": ["chrom", "chr", "chromosome", "query_contig"],
    "start": ["start", "pos", "position", "breakpoint1", "bp1"],
    "end": ["end", "stop", "breakpoint2", "bp2"],
    "sv_len": ["sv_len", "length", "size", "event_size", "span", "ancestral_segment_bp"],
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
    renamed = {}
    lower_map = {c.lower(): c for c in df.columns}
    for canonical, aliases in COLUMN_ALIASES.items():
        for alias in aliases:
            if alias.lower() in lower_map:
                renamed[lower_map[alias.lower()]] = canonical
                break
    out = df.rename(columns=renamed).copy()

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



def top_n(df: pd.DataFrame, col: str, n: int = 20) -> pd.DataFrame:
    vc = df[col].fillna("NA").value_counts().head(n)
    return vc.rename_axis(col).reset_index(name="count")



def add_size_bin(df: pd.DataFrame, source_col: str = "sv_len") -> pd.DataFrame:
    if source_col not in df.columns:
        return df
    bins = [-np.inf, 50, 100, 500, 1_000, 10_000, 100_000, np.inf]
    labels = ["<50", "50-100", "100-500", "500-1k", "1k-10k", "10k-100k", ">100k"]
    out = df.copy()
    out["size_bin"] = pd.cut(out[source_col], bins=bins, labels=labels)
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
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(df[x].astype(str), df[y])
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if rotate_xticks:
        ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    return fig



def plot_grouped_metric(df: pd.DataFrame, category: str, metric: str, series: str, title: str) -> plt.Figure:
    pivot = df.pivot_table(index=category, columns=series, values=metric, aggfunc="mean")
    fig, ax = plt.subplots(figsize=(10, 5))
    pivot.plot(kind="bar", ax=ax)
    ax.set_title(title)
    ax.set_xlabel(category)
    ax.set_ylabel(metric)
    ax.legend(title=series, bbox_to_anchor=(1.02, 1), loc="upper left")
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



def plot_box(df: pd.DataFrame, x: str, y: str, title: str, rotate_xticks: bool = False) -> plt.Figure:
    categories = [g[y].dropna().values for _, g in df.groupby(x)]
    labels = [str(k) for k, _ in df.groupby(x)]
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.boxplot(categories, labels=labels, showfliers=False)
    ax.set_title(title)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    if rotate_xticks:
        ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    return fig



def plot_stacked_counts(df: pd.DataFrame, index_col: str, stack_col: str, title: str, normalize: bool = False) -> plt.Figure:
    counts = pd.crosstab(df[index_col].fillna("NA"), df[stack_col].fillna("NA"), normalize="index" if normalize else False)
    fig, ax = plt.subplots(figsize=(10, 5))
    counts.plot(kind="bar", stacked=True, ax=ax)
    ax.set_title(title)
    ax.set_xlabel(index_col)
    ax.set_ylabel("Fraction" if normalize else "Count")
    ax.legend(title=stack_col, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    return fig



def plot_scatter(df: pd.DataFrame, x: str, y: str, title: str, xlabel: Optional[str] = None, ylabel: Optional[str] = None) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(7, 5))
    clean = df[[x, y]].dropna()
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

    real_df = harmonize_columns(real_df)
    real_df = add_size_bin(real_df)
    figs: List[FigureRecord] = []
    tables: List[pd.DataFrame] = []
    blocks: List[str] = []

    if metadata_df is not None and not metadata_df.empty:
        metadata_df = harmonize_columns(metadata_df)
        if "sample" in real_df.columns and "sample" in metadata_df.columns:
            real_df = real_df.merge(metadata_df, on="sample", how="left", suffixes=("", "_meta"))

    if "sample" in real_df.columns:
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

    if {"sample", "sv_type"}.issubset(real_df.columns):
        fig = plot_stacked_counts(real_df, "sample", "sv_type", "SV composition by sample", normalize=True)
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

    if {"caller", "sample"}.issubset(real_df.columns):
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

    if "gene" in bio_df.columns:
        gene_counts = top_n(bio_df, "gene", 20)
        tables.append(gene_counts)
        save_df(gene_counts, outdir / "biology_top_genes.tsv")
        fig = plot_bar(gene_counts, "gene", "count", "Top recurrent genes affected by SVs", "Gene", "Event count", rotate_xticks=True)
        encoded = write_figure(fig, outdir / "biology_top_genes.png")
        figs.append(FigureRecord(
            title="Top recurrent genes",
            filename="biology_top_genes.png",
            encoded_png=encoded,
            caption="Most recurrent genes linked to called structural variants."
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

    if {"effect", "gene"}.issubset(bio_df.columns):
        effect_counts = pd.crosstab(bio_df["gene"].fillna("NA"), bio_df["effect"].fillna("NA"))
        effect_counts = effect_counts.sum(axis=1).sort_values(ascending=False).head(15).index
        subset = bio_df[bio_df["gene"].isin(effect_counts)]
        fig = plot_stacked_counts(subset, "gene", "effect", "Functional effect composition of top genes", normalize=True)
        encoded = write_figure(fig, outdir / "biology_gene_effects.png")
        figs.append(FigureRecord(
            title="Functional effect composition",
            filename="biology_gene_effects.png",
            encoded_png=encoded,
            caption="Relative composition of annotated functional effects among recurrent genes."
        ))

    if {"expression", "cn"}.issubset(bio_df.columns):
        fig = plot_scatter(bio_df, "cn", "expression", "Copy number vs expression for SV-associated loci")
        encoded = write_figure(fig, outdir / "biology_cn_vs_expression.png")
        figs.append(FigureRecord(
            title="Copy number vs expression",
            filename="biology_cn_vs_expression.png",
            encoded_png=encoded,
            caption="Association between copy-number state and expression at SV-associated loci."
        ))

    if {"phenotype", "gene"}.issubset(bio_df.columns):
        top_genes = bio_df["gene"].fillna("NA").value_counts().head(10).index
        subset = bio_df[bio_df["gene"].isin(top_genes)]
        phen = pd.crosstab(subset["phenotype"].fillna("NA"), subset["gene"].fillna("NA"))
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
) -> str:
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
      <h2>3. Biological findings</h2>
      {bio_html}
      <div class=\"grid\">{render_figure_cards(bio_figs)}</div>
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

    sim_html, sim_figs, _ = build_simulated_section(sim_df, outdir)
    real_html, real_figs, _ = build_real_section(real_df, metadata_df, outdir)
    bio_html, bio_figs, _ = build_biology_section(bio_df, outdir)

    stats = summary_stats(sim_df, real_df, bio_df)
    html = build_html_report(args.title, stats, sim_html, sim_figs, real_html, real_figs, bio_html, bio_figs)

    report_path = outdir / "sv_visualization_report.html"
    report_path.write_text(html, encoding="utf-8")

    print(f"Report written to: {report_path}")
    print(f"Figures written under: {outdir}")


if __name__ == "__main__":
    main()
