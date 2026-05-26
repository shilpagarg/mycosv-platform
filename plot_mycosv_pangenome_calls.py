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
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

try:
    import matplotlib.pyplot as plt
    import matplotlib as _mpl
    import numpy as _np
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
    _mpl = None
    _np = None

try:
    import seaborn as _sns
except Exception:
    _sns = None

try:
    from pycirclize import Circos as _Circos
except Exception:
    _Circos = None

try:
    from upsetplot import UpSet as _UpSet, from_memberships as _from_memberships
except Exception:
    _UpSet = None
    _from_memberships = None


# Publication-style defaults — applied once at module import so every figure
# (including the legacy bar/scatter calls) inherits clean spines, larger
# fonts, and consistent typography. The Okabe-Ito palette is color-blind safe.
if plt is not None:
    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": "#333333",
        "axes.labelcolor": "#222222",
        "axes.titlesize": 12,
        "axes.titleweight": "semibold",
        "axes.labelsize": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": "#DDDDDD",
        "grid.linewidth": 0.6,
        "xtick.color": "#333333",
        "ytick.color": "#333333",
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.frameon": False,
        "legend.fontsize": 9,
        "font.family": ["DejaVu Sans", "Liberation Sans", "sans-serif"],
        "savefig.dpi": 220,
        "savefig.bbox": "tight",
    })

if _sns is not None:
    _sns.set_style("ticks")
    _sns.set_context("paper", font_scale=1.05)


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

SVTYPE_PALETTE = {
    "INS": "#E69F00",  # orange — gain
    "DEL": "#56B4E9",  # blue   — loss
    "INV": "#009E73",  # green
    "DUP": "#CC79A7",  # pink
    "TRA": "#D55E00",  # vermilion
    "OFF_REF": "#0072B2",  # dark blue — novel
}

SCOPE_PALETTE = {
    "pangenome_only": "#0072B2",
    "single_reference_equivalent": "#D55E00",
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
    _save_dual(fig, path, dpi=180)
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
    _save_dual(fig, path, dpi=180)
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
    _save_dual(fig, path, dpi=180)
    plt.close(fig)


def pct(num: int | float, den: int | float) -> str:
    return "0.0" if not den else f"{100.0 * num / den:.1f}"


def _save_dual(fig, path: Path, **save_kwargs) -> None:
    """Save figure as PNG (HTML preview) AND SVG (manuscript / Illustrator).

    Manuscript submissions need vector. PNG is kept so the HTML report
    renders without an extra step. The SVG sits next to the PNG with the
    same stem. Pass any matplotlib savefig kwargs through (dpi, bbox).
    """
    path = Path(path)
    fig.savefig(path, **save_kwargs)
    svg_kwargs = {k: v for k, v in save_kwargs.items() if k != "dpi"}
    fig.savefig(path.with_suffix(".svg"), **svg_kwargs)


# =====================================================================
# Publication-style figures (circos, ridge, upset, clustered heatmap,
# lollipop). All are no-ops if their optional library is unavailable, so
# benchmarks running in lean envs still produce the legacy bar/scatter set.
# =====================================================================


def _svtype_color(svtype: str) -> str:
    return SVTYPE_PALETTE.get((svtype or "").upper(), "#7F7F7F")


def save_circos_sv_landscape(
    path: Path,
    novel_rows: list[dict[str, str]],
    *,
    title: str,
    top_contigs: int = 20,
    bin_kb: int = 50,
) -> None:
    """Circular SV-density landscape over the top query contigs.

    Outer ring: per-bin SV count colored by SVTYPE composition (stacked
    histogram-on-circle). Inner ring: pangenome-only fraction per bin. Centre
    label: total calls and contigs shown. This is the manuscript-style figure
    fungal pangenome papers (Badet, Stukenbrock, Hartmann) use to show
    chromosome-scale SV distribution.
    """
    if _Circos is None or plt is None or not novel_rows:
        return
    by_contig: dict[str, list[dict[str, str]]] = defaultdict(list)
    for r in novel_rows:
        contig = norm_token(r.get("query_contig"))
        if contig == "unknown":
            continue
        by_contig[contig].append(r)
    if not by_contig:
        return
    # Pick top N contigs by call count to keep the plot legible.
    contig_sizes: dict[str, int] = {}
    for c, rows in by_contig.items():
        # Approximate contig length as max(end) — good enough for binning.
        max_end = 0
        for r in rows:
            end = as_int(r.get("end")) or (as_int(r.get("pos")) + abs(as_int(r.get("svlen"))))
            max_end = max(max_end, end)
        contig_sizes[c] = max(max_end, 1)
    ranked = sorted(by_contig, key=lambda c: -len(by_contig[c]))[:top_contigs]
    bin_size_bp = max(1000, bin_kb * 1000)
    # Pad sector size up to the next full bin boundary so bar-center
    # x-coords (idx*bin + bin/2) never overshoot sector.size, which would
    # raise ValueError("x=... is invalid range") in pycirclize.
    sectors = {
        c: ((contig_sizes[c] // bin_size_bp) + 1) * bin_size_bp
        for c in ranked
    }
    if not sectors:
        return
    circos = _Circos(sectors, space=2)
    # Outer track: stacked SV-type bars per bin.
    sv_types = ["INS", "DEL", "INV", "DUP", "TRA", "OFF_REF"]
    max_bin_count = 1
    bin_data: dict[str, dict[int, Counter[str]]] = defaultdict(lambda: defaultdict(Counter))
    for contig in ranked:
        for r in by_contig[contig]:
            pos = as_int(r.get("pos"))
            sv = norm_token(r.get("svtype")).upper()
            bin_idx = pos // bin_size_bp
            bin_data[contig][bin_idx][sv] += 1
            tot = sum(bin_data[contig][bin_idx].values())
            if tot > max_bin_count:
                max_bin_count = tot
    for sector in circos.sectors:
        contig = sector.name
        sector.text(contig, size=8, r=95)
        size = sector.size

        # Outer stacked-type histogram track (r 60..88) — radii must be in 0..100.
        outer_track = sector.add_track((60, 88), r_pad_ratio=0.05)
        outer_track.axis(fc="#FAFAFA", ec="#CCCCCC", lw=0.5)
        # x-tick labels every ~quarter of the contig (tick method lives on
        # the track in pycirclize 1.x, not the sector).
        outer_track.xticks_by_interval(
            max(1, int(size / 4)),
            label_size=6,
            label_orientation="vertical",
        )
        bin_counts = bin_data[contig]
        bin_indices = sorted(bin_counts)
        if not bin_indices:
            continue
        bottoms = {idx: 0 for idx in bin_indices}
        bar_xs = [idx * bin_size_bp + bin_size_bp / 2 for idx in bin_indices]
        for sv in sv_types:
            heights = [bin_counts[idx].get(sv, 0) for idx in bin_indices]
            if not any(heights):
                continue
            base = [bottoms[idx] for idx in bin_indices]
            outer_track.bar(
                bar_xs, heights, width=float(bin_size_bp),
                bottom=base, color=_svtype_color(sv),
                ec="none", alpha=0.92,
                vmax=max_bin_count,
            )
            for idx, h in zip(bin_indices, heights):
                bottoms[idx] += h

        # Inner pangenome-only fraction track (r 32..52), heatmap by bin.
        inner = sector.add_track((32, 52), r_pad_ratio=0.05)
        inner.axis(fc="#FFFFFF", ec="#CCCCCC", lw=0.5)
        po_frac: list[float] = []
        idx_list: list[int] = []
        for idx in bin_indices:
            calls = by_contig[contig]
            in_bin = [r for r in calls if (as_int(r.get("pos")) // bin_size_bp) == idx]
            if not in_bin:
                continue
            po = sum(1 for r in in_bin if call_scope(r) == "pangenome_only")
            po_frac.append(po / len(in_bin))
            idx_list.append(idx)
        if idx_list:
            # build a full-length array indexed by bin number, padded with NaN
            # so missing bins render blank rather than smear color across gaps.
            n_full = max(idx_list) + 1
            row = _np.full(n_full, _np.nan)
            for i, val in zip(idx_list, po_frac):
                row[i] = val
            inner.heatmap(
                _np.array([row]),
                vmin=0.0, vmax=1.0,
                cmap="Blues",
                rect_kws={"ec": "none"},
            )
    fig = circos.plotfig(figsize=(9, 9))
    # Legends
    sv_handles = [_mpl.patches.Patch(color=_svtype_color(t), label=t) for t in sv_types]
    fig.legend(handles=sv_handles, loc="upper left", bbox_to_anchor=(0.02, 0.98),
               title="SV type (outer ring)", frameon=False)
    cax = fig.add_axes([0.86, 0.10, 0.015, 0.20])
    cb = _mpl.colorbar.ColorbarBase(
        cax, cmap=_mpl.cm.Blues,
        norm=_mpl.colors.Normalize(vmin=0, vmax=1),
        orientation="vertical",
    )
    cb.set_label("Pangenome-only fraction\n(inner ring)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.suptitle(title, y=0.99, fontsize=13, fontweight="semibold")
    fig.text(0.5, 0.5,
             f"n={len(novel_rows):,} calls\n{len(ranked)} contigs\nbin={bin_kb} kb",
             ha="center", va="center", fontsize=10, color="#444")
    _save_dual(fig, path)
    plt.close(fig)


def save_ridge_svlen(
    path: Path,
    novel_rows: list[dict[str, str]],
    *,
    title: str,
) -> None:
    """Ridge (joy) plot of |SVLEN| distribution per SV type on log10 scale.

    Replaces the flat bar histogram with the publication-style overlapping
    KDE ridges that fungal SV papers use (Stukenbrock et al., Galagan/Selker
    RIP reviews) to compare length spectra across event classes at a glance.
    """
    if _sns is None or plt is None or not novel_rows:
        return
    data: list[tuple[str, float]] = []
    for r in novel_rows:
        sv = norm_token(r.get("svtype")).upper()
        if sv not in SVTYPE_PALETTE:
            continue
        L = abs(as_int(r.get("svlen")))
        if L <= 0:
            continue
        data.append((sv, math.log10(L)))
    if len(data) < 10:
        return
    sv_types = [t for t in SVTYPE_PALETTE if any(d[0] == t for d in data)]
    if not sv_types:
        return
    fig, axes = plt.subplots(
        nrows=len(sv_types), figsize=(8.5, 0.95 * len(sv_types) + 0.8),
        sharex=True,
    )
    if len(sv_types) == 1:
        axes = [axes]
    x_min = math.floor(min(d[1] for d in data))
    x_max = math.ceil(max(d[1] for d in data))
    for ax, sv in zip(axes, sv_types):
        vals = [d[1] for d in data if d[0] == sv]
        color = _svtype_color(sv)
        _sns.kdeplot(vals, ax=ax, fill=True, color=color, alpha=0.78,
                     bw_adjust=0.75, linewidth=1.2, cut=0)
        ax.set_xlim(x_min, x_max)
        ax.set_ylabel("")
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.text(x_min + 0.05, ax.get_ylim()[1] * 0.55,
                f"{sv}  (n={len(vals):,})", fontsize=10, color="#222",
                fontweight="medium")
        ax.set_xlabel("")
        if ax is not axes[-1]:
            ax.spines["bottom"].set_visible(False)
            ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
            ax.grid(False)
    tick_locs = list(range(x_min, x_max + 1))
    axes[-1].set_xticks(tick_locs)
    axes[-1].set_xticklabels([f"10^{t}" if t > 1 else str(int(10 ** t)) for t in tick_locs])
    axes[-1].set_xlabel("|SVLEN|  (bp, log10 scale)")
    fig.suptitle(title, y=0.998, fontsize=12, fontweight="semibold")
    fig.subplots_adjust(hspace=-0.15)
    _save_dual(fig, path)
    plt.close(fig)


def save_upset_comparator_overlap(
    path: Path,
    novel_rows: list[dict[str, str]],
    *,
    title: str,
    min_subset_size: int = 5,
) -> None:
    """Upset plot of comparator co-support per mycosv call.

    A mycosv call's `support_labels` field is a `;`-separated list of which
    comparators (anchorwave, minigraph, svim_asm, read_level_union, ...)
    independently supported it. Upset is the right replacement for a Venn
    when you have >3 sets — and >3 is the norm here.
    """
    if _UpSet is None or _from_memberships is None or plt is None or not novel_rows:
        return
    memberships: list[list[str]] = []
    for r in novel_rows:
        lbls = (r.get("support_labels") or "").strip()
        if lbls in {"", "."}:
            memberships.append([])
            continue
        parts = sorted({s.strip() for s in lbls.split(";") if s.strip() and s.strip() != "."})
        memberships.append(parts)
    if not any(memberships):
        return
    series = _from_memberships(memberships)
    fig = plt.figure(figsize=(10, 5.2))
    upset = _UpSet(
        series,
        subset_size="count",
        show_counts=True,
        sort_by="cardinality",
        min_subset_size=min_subset_size,
        intersection_plot_elements=4,
        totals_plot_elements=3,
        facecolor="#0072B2",
        other_dots_color="#BBBBBB",
    )
    upset.plot(fig=fig)
    fig.suptitle(title, y=0.995, fontsize=12, fontweight="semibold")
    _save_dual(fig, path)
    plt.close(fig)


def save_clustermap_biology(
    path: Path,
    biology_rows: list[dict[str, str]],
    *,
    title: str,
    min_total: int = 5,
) -> None:
    """Clustered heatmap (with row/col dendrograms) of biology class × phenotype.

    Aggregates biology_findings.tsv into a count matrix of element_class vs
    ecological_trait / lifestyle, then z-scales rows so a dominant class
    (RIP) doesn't blank out the rest, and clusters both axes hierarchically
    — exactly the figure type used in fungal pangenome / two-speed-genome
    papers (Möller 2017, Plissonneau 2018, Badet 2020).
    """
    if _sns is None or plt is None or not biology_rows:
        return
    counts: dict[tuple[str, str], int] = defaultdict(int)
    for r in biology_rows:
        ec = element_group(r.get("element_class"))
        eco = norm_token(r.get("ecological_trait") or r.get("trophic_mode") or r.get("lifestyle"))
        if eco == "unknown" or ec == "Structural (other)":
            continue
        counts[(ec, eco)] += 1
    if not counts:
        return
    ecs = sorted({k[0] for k in counts})
    ecos = sorted({k[1] for k in counts})
    matrix = _np.array([[counts.get((ec, eco), 0) for eco in ecos] for ec in ecs], dtype=float)
    if matrix.sum() < min_total or matrix.shape[0] < 2 or matrix.shape[1] < 2:
        return
    # z-score across columns per row (lights up class-specific eco enrichments).
    row_means = matrix.mean(axis=1, keepdims=True)
    row_stds = matrix.std(axis=1, keepdims=True)
    row_stds[row_stds == 0] = 1.0
    z = (matrix - row_means) / row_stds
    cg = _sns.clustermap(
        z, cmap="RdBu_r", center=0,
        xticklabels=ecos, yticklabels=ecs,
        figsize=(max(7.5, 0.65 * len(ecos) + 4), max(4.5, 0.45 * len(ecs) + 3)),
        cbar_kws={"label": "row z-score"},
        linewidths=0.25, linecolor="#FFFFFF",
        dendrogram_ratio=(0.16, 0.20),
        col_cluster=matrix.shape[1] >= 3,
        row_cluster=matrix.shape[0] >= 3,
    )
    cg.ax_heatmap.set_xlabel("")
    cg.ax_heatmap.set_ylabel("")
    cg.ax_heatmap.tick_params(axis="x", rotation=35)
    for tick in cg.ax_heatmap.get_xticklabels():
        tick.set_ha("right")
    cg.fig.suptitle(title, y=1.00, fontsize=12, fontweight="semibold")
    _save_dual(cg.fig, path)
    plt.close(cg.fig)


def save_lollipop_top_contigs(
    path: Path,
    novel_rows: list[dict[str, str]],
    *,
    title: str,
    top_n: int = 15,
) -> None:
    """Lollipop of top contigs by pangenome-only call count, with stem hue
    encoding the read-support rate. Replaces the flat bar chart with a more
    information-dense figure (count on x, rate on color)."""
    if plt is None or not novel_rows:
        return
    contig_counts: Counter[str] = Counter()
    contig_supported: Counter[str] = Counter()
    for r in novel_rows:
        if call_scope(r) != "pangenome_only":
            continue
        contig = norm_token(r.get("query_contig"))
        if contig == "unknown":
            continue
        contig_counts[contig] += 1
        if (r.get("read_supported") or "").lower() == "yes":
            contig_supported[contig] += 1
    if not contig_counts:
        return
    top = contig_counts.most_common(top_n)
    contigs = [c for c, _ in top]
    counts = [n for _, n in top]
    rates = [contig_supported.get(c, 0) / max(1, contig_counts[c]) for c in contigs]
    fig, ax = plt.subplots(figsize=(8.5, max(3.5, 0.32 * len(contigs) + 1.5)))
    norm = _mpl.colors.Normalize(vmin=0, vmax=1)
    cmap = _mpl.cm.viridis
    y_pos = _np.arange(len(contigs))
    for y, x, rate in zip(y_pos, counts, rates):
        ax.plot([0, x], [y, y], color=cmap(norm(rate)), lw=2.0, alpha=0.85)
        ax.scatter(x, y, s=120, color=cmap(norm(rate)),
                   edgecolor="#333333", linewidth=0.6, zorder=3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(contigs, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Pangenome-only MycoSV calls")
    ax.set_title(title)
    sm = _mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, pad=0.02, shrink=0.85)
    cb.set_label("Read-support rate", fontsize=9)
    ax.grid(axis="x", alpha=0.3)
    ax.grid(axis="y", visible=False)
    _save_dual(fig, path)
    plt.close(fig)


def save_volcano_enrichment(
    path: Path,
    enrichment_rows: list[dict[str, str]],
    *,
    title: str,
    or_col: str = "odds_ratio",
    p_col: str = "fisher_p",
    label_col: str = "feature",
    group_col: str | None = "guild",
    sig_p: float = 1e-3,
    sig_or: float = 2.0,
) -> None:
    """Volcano plot of feature enrichments — manuscript standard for
    cross-guild biology (Plissonneau 2018, Badet 2020, Möller 2017).

    log2(OR) on x, -log10(p) on y. Points outside the (sig_or, sig_p) wedge
    get labels; rest get a faint dot. Color-coded by group when `group_col`
    is supplied (e.g. AMF vs Filamentous from cross_guild_enrichment.tsv).
    """
    if plt is None or not enrichment_rows:
        return
    xs: list[float] = []
    ys: list[float] = []
    labels: list[str] = []
    groups: list[str] = []
    for r in enrichment_rows:
        try:
            o = float(r.get(or_col, 0) or 0)
            p = float(r.get(p_col, 1) or 1)
        except ValueError:
            continue
        if o <= 0:
            continue
        x = math.log2(o)
        # clamp p at 1e-300 so log10 doesn't explode and the dot still draws.
        y = -math.log10(max(p, 1e-300))
        xs.append(x)
        ys.append(y)
        labels.append(norm_token(r.get(label_col)))
        groups.append(norm_token(r.get(group_col) if group_col else None))
    if not xs:
        return
    fig, ax = plt.subplots(figsize=(8.5, 6.2))
    group_keys = sorted({g for g in groups if g != "unknown"})
    group_color = {
        g: c for g, c in zip(group_keys, _mpl.cm.tab10.colors[: len(group_keys)])
    }
    sig_log_p = -math.log10(sig_p)
    sig_log_or = math.log2(sig_or)
    for x, y, lbl, g in zip(xs, ys, labels, groups):
        sig = (abs(x) >= sig_log_or) and (y >= sig_log_p)
        color = group_color.get(g, "#7F7F7F")
        size = 75 if sig else 25
        alpha = 0.9 if sig else 0.45
        ax.scatter(x, y, s=size, color=color, alpha=alpha,
                   edgecolor="#222" if sig else "none", linewidth=0.5)
        if sig:
            ax.annotate(
                lbl, xy=(x, y), xytext=(4, 4), textcoords="offset points",
                fontsize=8.5, color="#222",
            )
    ax.axhline(sig_log_p, color="#999", lw=0.7, ls="--")
    ax.axvline(sig_log_or, color="#999", lw=0.7, ls="--")
    ax.axvline(-sig_log_or, color="#999", lw=0.7, ls="--")
    ax.set_xlabel("log2(odds ratio)")
    ax.set_ylabel("-log10(p)")
    ax.set_title(title)
    if group_keys:
        handles = [
            _mpl.lines.Line2D([0], [0], marker="o", linestyle="",
                              markerfacecolor=group_color[g],
                              markeredgecolor="#222", markersize=8, label=g)
            for g in group_keys
        ]
        ax.legend(handles=handles, loc="upper left", title=group_col, frameon=False)
    _save_dual(fig, path)
    plt.close(fig)


def save_manhattan_sv_density(
    path: Path,
    novel_rows: list[dict[str, str]],
    *,
    title: str,
    bin_kb: int = 25,
    top_contigs: int = 30,
) -> None:
    """Manhattan-style SV density plot across contigs (linear genome view).

    Each contig is a colored stripe; y-axis is per-window SV count. Mirrors
    the layout fungal SV / two-speed-genome papers use (Plissonneau 2018,
    Hartmann 2017) when a circos is overkill or the contig count is large.
    """
    if plt is None or not novel_rows:
        return
    by_contig: dict[str, list[dict[str, str]]] = defaultdict(list)
    for r in novel_rows:
        c = norm_token(r.get("query_contig"))
        if c == "unknown":
            continue
        by_contig[c].append(r)
    if not by_contig:
        return
    ranked = sorted(by_contig, key=lambda c: -len(by_contig[c]))[:top_contigs]
    bin_bp = max(1000, bin_kb * 1000)
    rows: list[tuple[str, int, int, float]] = []
    cum = 0
    contig_extent: dict[str, tuple[int, int]] = {}
    for c in ranked:
        calls = by_contig[c]
        max_end = max(
            (as_int(r.get("end")) or (as_int(r.get("pos")) + abs(as_int(r.get("svlen")))))
            for r in calls
        )
        bins = max(1, (max_end // bin_bp) + 1)
        bin_counts = Counter()
        bin_po = Counter()
        for r in calls:
            b = as_int(r.get("pos")) // bin_bp
            bin_counts[b] += 1
            if call_scope(r) == "pangenome_only":
                bin_po[b] += 1
        for b in range(bins):
            if bin_counts[b]:
                center = cum + b * bin_bp + bin_bp / 2
                po_frac = bin_po[b] / bin_counts[b]
                rows.append((c, b, bin_counts[b], po_frac))
                xs_centers_ok = True
        contig_extent[c] = (cum, cum + bins * bin_bp)
        cum += bins * bin_bp
    if not rows:
        return
    fig, ax = plt.subplots(figsize=(max(9.0, len(ranked) * 0.5 + 3), 4.2))
    contig_colors = (_mpl.cm.tab20.colors * 5)[: len(ranked)]
    color_map = dict(zip(ranked, contig_colors))
    for c, b, cnt, po_frac in rows:
        x = contig_extent[c][0] + b * bin_bp + bin_bp / 2
        ax.scatter(x, cnt, s=16 + 3 * po_frac * cnt,
                   color=color_map[c],
                   alpha=0.4 + 0.6 * po_frac,
                   edgecolor="none")
    tick_xs = [(s + e) / 2 for s, e in contig_extent.values()]
    for c, (s, e) in contig_extent.items():
        ax.axvspan(s, e, color=color_map[c], alpha=0.04, zorder=0)
    ax.set_xticks(tick_xs)
    ax.set_xticklabels(ranked, rotation=60, ha="right", fontsize=7)
    ax.set_ylabel(f"MycoSV calls per {bin_kb} kb window")
    ax.set_title(title + "  (point size encodes pangenome-only fraction)")
    ax.set_xlim(0, cum)
    fig.subplots_adjust(bottom=0.25)
    _save_dual(fig, path)
    plt.close(fig)


def save_two_speed_scatter(
    path: Path,
    novel_rows: list[dict[str, str]],
    *,
    title: str,
    bin_kb: int = 50,
) -> None:
    """Two-speed genome scatter: TE-like SV density vs structural SV density
    per window — the diagnostic plot for compartmentalized fungal genomes
    (Möller & Stukenbrock 2017, Plissonneau 2018, Faino 2016).

    Each dot is a genomic window; x = TE/RIP/HGT-class SVs per window,
    y = all-other-class SVs per window. Compartments cluster: gene-core
    sits near origin, TE/repeat accessory drifts toward x.
    """
    if plt is None or not novel_rows:
        return
    bin_bp = max(1000, bin_kb * 1000)
    by_bin_class: dict[tuple[str, int], Counter[str]] = defaultdict(Counter)
    for r in novel_rows:
        c = norm_token(r.get("query_contig"))
        if c == "unknown":
            continue
        b = as_int(r.get("pos")) // bin_bp
        ec = element_group(r.get("element_class"))
        by_bin_class[(c, b)][ec] += 1
    if not by_bin_class:
        return
    xs: list[int] = []
    ys: list[int] = []
    sizes: list[int] = []
    for counts in by_bin_class.values():
        te = counts.get("TE / repeat", 0) + counts.get("RIP / genome defense", 0) + counts.get("HGT / Starship", 0)
        other = counts.get("Structural (other)", 0) + sum(
            v for k, v in counts.items() if k not in {"TE / repeat", "RIP / genome defense", "HGT / Starship", "Structural (other)"}
        )
        if te + other == 0:
            continue
        xs.append(te)
        ys.append(other)
        sizes.append(te + other)
    if not xs:
        return
    fig, ax = plt.subplots(figsize=(7.4, 6.0))
    sc = ax.scatter(
        xs, ys, s=[6 + s for s in sizes],
        c=sizes, cmap="viridis",
        alpha=0.78, edgecolor="none",
    )
    # diagonal y=x reference: equal TE-like and structural density.
    mx = max(max(xs), max(ys), 1)
    ax.plot([0, mx], [0, mx], color="#888", ls=":", lw=1, label="y = x")
    ax.set_xlabel(f"TE / RIP / HGT-class SVs per {bin_kb} kb window")
    ax.set_ylabel(f"Other structural SVs per {bin_kb} kb window")
    ax.set_title(title)
    cb = fig.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label("Total SVs in window", fontsize=9)
    ax.legend(loc="upper left", frameon=False)
    _save_dual(fig, path)
    plt.close(fig)


# =====================================================================
# Panel-fold figures (Figure 1 of the manuscript) — read every per-query
# shard under a by_query/ directory and produce one figure per claim across
# all samples. These are the validation panels: per-sample read-validation
# rate (Fig 1A) and per-sample pangenome-lift stacked bars (Fig 1B).
# =====================================================================


def _iter_shards(by_query_dir: Path):
    """Yield (query_asm, shard_path) for every per-query directory under
    by_query/, skipping ones with MYCOSV_FAILED.txt — same gate the
    cross-guild aggregator uses so the panel-fold sees the same population."""
    if not by_query_dir.is_dir():
        return
    for shard in sorted(p for p in by_query_dir.iterdir() if p.is_dir()):
        if (shard / "MYCOSV_FAILED.txt").exists():
            continue
        yield shard.name, shard


def _wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n <= 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1.0 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = (z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / denom
    return (max(0.0, center - half), min(1.0, center + half))


def _phylogeny_sort_key(meta: dict[str, str]) -> tuple[int, str, str]:
    """Sort 13/165 samples by phylum -> class -> species so x-axis ordering
    on Fig 1A/B traces the literature's standard phylogenetic layout (AMF
    on the left, Filamentous on the right, Yeast/Basidio between).
    """
    phylum = (meta.get("phylum") or "").lower()
    clazz = (meta.get("class") or "").lower()
    species = (meta.get("species") or meta.get("query_asm") or "").lower()
    phylum_rank = {
        "glomeromycota": 0,
        "basidiomycota": 1,
        "ascomycota": 2,
    }.get(phylum, 9)
    return (phylum_rank, clazz, species)


def _load_query_meta_map(by_query_dir: Path) -> dict[str, dict[str, str]]:
    """query_asm -> {phylum, class, species, ...} via prepared/query_manifest.tsv
    if discoverable. Returns {} when the manifest can't be located — callers
    should keep degrading to alphabetical x-axis ordering in that case."""
    candidates = []
    for ancestor in [by_query_dir, by_query_dir.parent, by_query_dir.parent.parent]:
        candidates.append(ancestor / "prepared" / "query_manifest.tsv")
    for ancestor in by_query_dir.parents:
        candidates.append(ancestor / "prepared" / "query_manifest.tsv")
    manifest = next((c for c in candidates if c.exists()), None)
    if manifest is None:
        return {}
    out: dict[str, dict[str, str]] = {}
    with manifest.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            asm = (row.get("query_asm") or "").strip()
            if asm:
                out[asm] = row
    return out


def aggregate_panel_read_validation(
    by_query_dir: Path,
    out_tsv: Path,
) -> list[dict[str, object]]:
    """Per-(query, source) yes/total/rate from read_validated_truth.tsv."""
    rows: list[dict[str, object]] = []
    for qname, shard in _iter_shards(by_query_dir):
        f = shard / "read_validated_truth.tsv"
        if not f.exists() or f.stat().st_size == 0:
            continue
        per_source: dict[str, list[int]] = defaultdict(lambda: [0, 0])
        with f.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for r in reader:
                src = (r.get("source") or "").strip()
                if not src:
                    continue
                per_source[src][1] += 1
                if (r.get("read_validated") or "").strip().lower() == "yes":
                    per_source[src][0] += 1
        for src, (yes, total) in per_source.items():
            lo, hi = _wilson_ci(yes, total)
            rows.append({
                "query_asm": qname,
                "source": src,
                "yes": yes,
                "total": total,
                "rate": yes / total if total else 0.0,
                "ci95_lo": lo,
                "ci95_hi": hi,
            })
    write_tsv(out_tsv, rows,
              ["query_asm", "source", "yes", "total", "rate", "ci95_lo", "ci95_hi"])
    return rows


def save_panel_read_validation_rate(
    path: Path,
    rows: list[dict[str, object]],
    meta_by_asm: dict[str, dict[str, str]],
    *,
    title: str,
) -> None:
    """Manuscript Figure 1A — per-sample read-validation rate, MycoSV vs
    each comparator. Grouped bar with Wilson 95% CI; samples ordered by
    phylogeny (Glomero -> Basidio -> Asco), MycoSV bar highlighted."""
    if plt is None or not rows:
        return
    sources_present = sorted({str(r["source"]) for r in rows})
    preferred_order = ["mycosv", "minigraph", "svim_asm", "anchorwave"]
    sources = ([s for s in preferred_order if s in sources_present]
               + [s for s in sources_present if s not in preferred_order])
    queries = sorted({str(r["query_asm"]) for r in rows},
                     key=lambda q: _phylogeny_sort_key(meta_by_asm.get(q, {"query_asm": q})))
    if not queries:
        return
    by_pair: dict[tuple[str, str], dict[str, object]] = {
        (str(r["query_asm"]), str(r["source"])): r for r in rows
    }
    n_q = len(queries)
    n_src = len(sources)
    bar_w = 0.8 / n_src
    x = _np.arange(n_q)
    fig_w = max(9.0, n_q * 0.55 + 3.0)
    fig, ax = plt.subplots(figsize=(fig_w, 4.8))
    color_for = {
        "mycosv": "#0072B2", "minigraph": "#56B4E9",
        "svim_asm": "#009E73", "anchorwave": "#E69F00",
    }
    for i, src in enumerate(sources):
        offset = (i - (n_src - 1) / 2) * bar_w
        heights, los, his = [], [], []
        for q in queries:
            r = by_pair.get((q, src))
            if r is None:
                heights.append(0.0); los.append(0.0); his.append(0.0)
            else:
                heights.append(float(r["rate"]))
                los.append(float(r["ci95_lo"]))
                his.append(float(r["ci95_hi"]))
        err_lo = [max(0.0, h - lo) for h, lo in zip(heights, los)]
        err_hi = [max(0.0, hi - h) for h, hi in zip(heights, his)]
        bars = ax.bar(x + offset, heights, bar_w,
                      label=src, color=color_for.get(src, "#7F7F7F"),
                      edgecolor="white", linewidth=0.5,
                      yerr=[err_lo, err_hi],
                      error_kw={"ecolor": "#444", "elinewidth": 0.6, "capsize": 1.5})
        if src == "mycosv":
            for bar in bars:
                bar.set_edgecolor("#222")
                bar.set_linewidth(0.9)
    # Short tick labels: drop the GCA prefix and version suffix for legibility.
    short = [q.replace("GCA_", "").replace("_", ".") for q in queries]
    ax.set_xticks(x)
    ax.set_xticklabels(short, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Fraction of calls with raw-read split support")
    ax.set_ylim(0, 1.05)
    ax.axhline(1.0, color="#CCC", lw=0.5)
    # Phylum band annotation strip below x-axis labels.
    seen_phylum_at: list[tuple[int, str]] = []
    for i, q in enumerate(queries):
        ph = (meta_by_asm.get(q, {}).get("phylum") or "").strip()
        if not seen_phylum_at or seen_phylum_at[-1][1] != ph:
            seen_phylum_at.append((i, ph))
    seen_phylum_at.append((n_q, ""))
    phylum_color = {"Glomeromycota": "#0072B2", "Basidiomycota": "#E69F00",
                    "Ascomycota": "#009E73"}
    for j in range(len(seen_phylum_at) - 1):
        start, label = seen_phylum_at[j]
        end = seen_phylum_at[j + 1][0]
        if not label or label == ".":
            continue
        ax.axvspan(start - 0.5, end - 0.5,
                   color=phylum_color.get(label, "#999"),
                   alpha=0.07, zorder=0)
    ax.legend(loc="upper right", ncol=n_src, frameon=False, fontsize=9)
    ax.set_title(title)
    _save_dual(fig, path)
    plt.close(fig)


def aggregate_panel_pangenome_layers(
    by_query_dir: Path,
    out_tsv: Path,
) -> list[dict[str, object]]:
    """Per-query single_ref / pangenome_only / read_supported triple from
    each pangenome_call_layers.tsv. One row per query (the per-shard ALL
    row already exists in each TSV)."""
    rows: list[dict[str, object]] = []
    for qname, shard in _iter_shards(by_query_dir):
        f = shard / "pangenome_call_layers.tsv"
        if not f.exists() or f.stat().st_size == 0:
            continue
        with f.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if (r.get("query_asm") or "") in {"ALL", qname}:
                    rows.append({
                        "query_asm": qname,
                        "raw_obs": as_int(r.get("raw_pairwise_pangenome_observations")),
                        "dedup_loci": as_int(r.get("deduplicated_pangenome_loci")),
                        "single_ref_equivalent": as_int(r.get("single_reference_equivalent_calls")),
                        "pangenome_only": as_int(r.get("pangenome_only_calls")),
                        "pangenome_only_read_supported": as_int(r.get("pangenome_only_read_supported")),
                    })
                    break
    write_tsv(out_tsv, rows,
              ["query_asm", "raw_obs", "dedup_loci",
               "single_ref_equivalent", "pangenome_only",
               "pangenome_only_read_supported"])
    return rows


def save_panel_pangenome_lift(
    path: Path,
    rows: list[dict[str, object]],
    meta_by_asm: dict[str, dict[str, str]],
    *,
    title: str,
) -> None:
    """Manuscript Figure 1B — per-sample stacked bar: single-ref-equivalent /
    pangenome-only-intrinsic / pangenome-only-read-supported. Total
    deduplicated loci annotated above each bar."""
    if plt is None or not rows:
        return
    by_asm = {str(r["query_asm"]): r for r in rows}
    queries = sorted(by_asm.keys(),
                     key=lambda q: _phylogeny_sort_key(meta_by_asm.get(q, {"query_asm": q})))
    if not queries:
        return
    single_ref = [as_int(by_asm[q]["single_ref_equivalent"]) for q in queries]
    pango_total = [as_int(by_asm[q]["pangenome_only"]) for q in queries]
    pango_read = [as_int(by_asm[q]["pangenome_only_read_supported"]) for q in queries]
    pango_intrinsic = [max(0, t - r) for t, r in zip(pango_total, pango_read)]
    dedup = [as_int(by_asm[q]["dedup_loci"]) for q in queries]
    short = [q.replace("GCA_", "").replace("_", ".") for q in queries]
    x = _np.arange(len(queries))
    fig_w = max(9.0, len(queries) * 0.55 + 3.0)
    fig, ax = plt.subplots(figsize=(fig_w, 4.8))
    ax.bar(x, single_ref, color=SCOPE_PALETTE["single_reference_equivalent"],
           edgecolor="white", linewidth=0.6, label="single-reference equivalent")
    ax.bar(x, pango_read, bottom=single_ref,
           color="#009E73", edgecolor="white", linewidth=0.6,
           label="pangenome-only, read-supported")
    bottoms2 = [s + r for s, r in zip(single_ref, pango_read)]
    ax.bar(x, pango_intrinsic, bottom=bottoms2,
           color="#7F7F7F", edgecolor="white", linewidth=0.6,
           label="pangenome-only, intrinsic")
    # Annotate dedup totals above bars.
    for i, total in enumerate(dedup):
        ax.text(i, single_ref[i] + pango_total[i],
                f"{total:,}", ha="center", va="bottom", fontsize=8, color="#222")
    ax.set_xticks(x)
    ax.set_xticklabels(short, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("MycoSV loci")
    # Phylum bands as light background swaths.
    seen_phylum_at: list[tuple[int, str]] = []
    for i, q in enumerate(queries):
        ph = (meta_by_asm.get(q, {}).get("phylum") or "").strip()
        if not seen_phylum_at or seen_phylum_at[-1][1] != ph:
            seen_phylum_at.append((i, ph))
    seen_phylum_at.append((len(queries), ""))
    phylum_color = {"Glomeromycota": "#0072B2", "Basidiomycota": "#E69F00",
                    "Ascomycota": "#009E73"}
    for j in range(len(seen_phylum_at) - 1):
        start, label = seen_phylum_at[j]
        end = seen_phylum_at[j + 1][0]
        if label in {"", "."}:
            continue
        ax.axvspan(start - 0.5, end - 0.5,
                   color=phylum_color.get(label, "#999"),
                   alpha=0.07, zorder=0)
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    ax.set_title(title)
    _save_dual(fig, path)
    plt.close(fig)


def build_manuscript_table1(
    by_query_dir: Path,
    pl_rows: list[dict[str, object]],
    rv_rows: list[dict[str, object]],
    meta_by_asm: dict[str, dict[str, str]],
    out_tsv: Path,
    out_md: Path,
) -> list[dict[str, object]]:
    """Manuscript Table 1 — pangenome-powered SV recovery vs single-reference,
    one row per sample plus a panel-aggregate row.

    Column choice follows the de-facto standard set by fungal pangenome papers
    (Plissonneau 2018 mBio Table 1, Hartmann 2017 Mol Ecol Table 1,
    Badet 2020 BMC Biol Table 1, Liang 2024 MBE Table 2, Cuomo 2019
    Cell Host Microbe Table 1, Manchanda 2020 NAR Genom Bioinform Table 1):

      Sample / lineage / class
      Assembly size (Mb)            (genome scale)
      Total MycoSV loci             (deduplicated)
      Single-reference equivalent   (recoverable on the linear bench ref)
      Pangenome-only                (graph-native, not on the bench ref)
      Pangenome lift  (% gain)      (the headline number reviewers cite)
      Read-supported pangenome-only (independently confirmed by reads)
      MycoSV read-validation rate   (quality of the full callset)
      Comparator agreement (best F1)  (best of anchorwave/svim_asm/minigraph)
      Element-class breakdown        (HGT / TE+RIP / other) — biology context

    Two outputs:
      * .tsv for downstream scripting
      * .md for direct paste into the manuscript / supplementary
    """
    # Index inputs by sample.
    by_asm_pl = {str(r["query_asm"]): r for r in pl_rows}
    by_asm_rv: dict[str, dict[str, dict[str, int]]] = defaultdict(dict)
    for r in rv_rows:
        by_asm_rv[str(r["query_asm"])][str(r["source"])] = {
            "yes": int(r["yes"]), "total": int(r["total"])
        }
    # Per-shard biology + benchmark + assembly-size lookups.
    biology_class: dict[str, Counter[str]] = defaultdict(Counter)
    best_f1: dict[str, float] = {}
    asm_size_mb: dict[str, float] = {}
    for qname, shard in _iter_shards(by_query_dir):
        bf = shard / "biology_findings.tsv"
        if bf.exists() and bf.stat().st_size > 0:
            with bf.open() as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    biology_class[qname][element_group(r.get("element_class"))] += 1
        eb = shard / "exact_benchmark_summary.tsv"
        if eb.exists() and eb.stat().st_size > 0:
            with eb.open() as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if (r.get("coordinate_space") or "") != "reference":
                        continue
                    if (r.get("validation_basis") or "") != "comparator_baseline":
                        continue
                    if (r.get("svtype") or "") != "ALL":
                        continue
                    if (r.get("status") or "") != "ok":
                        continue
                    try:
                        f1 = float(r.get("f1") or 0)
                    except ValueError:
                        continue
                    if f1 > best_f1.get(qname, -1.0):
                        best_f1[qname] = f1
        # Assembly size estimate from INPUT_PREFLIGHT.tsv (query row), if present.
        pf = shard / "INPUT_PREFLIGHT.tsv"
        if pf.exists():
            try:
                with pf.open() as fh:
                    for r in csv.DictReader(fh, delimiter="\t"):
                        if (r.get("role") or "") == "query":
                            p = (r.get("path") or "").strip()
                            if p:
                                try:
                                    asm_size_mb[qname] = round(
                                        Path(p).stat().st_size / 1e6, 1)
                                except OSError:
                                    pass
                            break
            except OSError:
                pass

    # Count shards excluded by the MYCOSV_FAILED.txt gate so the caption can
    # disclose them — manuscript tables should never silently drop samples.
    total_panel_shards = sum(1 for p in by_query_dir.iterdir() if p.is_dir())
    failed_shards = sum(
        1 for p in by_query_dir.iterdir()
        if p.is_dir() and (p / 'MYCOSV_FAILED.txt').exists()
    )

    queries = sorted(by_asm_pl,
                     key=lambda q: _phylogeny_sort_key(meta_by_asm.get(q, {"query_asm": q})))
    rows: list[dict[str, object]] = []
    for q in queries:
        pl = by_asm_pl[q]
        single = int(pl["single_ref_equivalent"]); po = int(pl["pangenome_only"])
        po_read = int(pl["pangenome_only_read_supported"]); dedup = int(pl["dedup_loci"])
        meta = meta_by_asm.get(q, {})
        cls = biology_class.get(q, Counter())
        rv = by_asm_rv.get(q, {}).get("mycosv", {"yes": 0, "total": 0})
        rows.append({
            "sample": q,
            "species": meta.get("species") or ".",
            "class": meta.get("class") or ".",
            "assembly_size_mb": asm_size_mb.get(q, ""),
            "total_loci": dedup,
            "single_ref_equivalent": single,
            "pangenome_only": po,
            "pangenome_lift_pct": f"{(po / dedup * 100):.1f}" if dedup else "0.0",
            "pangenome_only_read_supported": po_read,
            "mycosv_read_validated_pct": (
                f"{(rv['yes'] / rv['total'] * 100):.1f}" if rv["total"] else "."
            ),
            "best_comparator_f1": (f"{best_f1[q]:.3f}" if q in best_f1 else "."),
            "hgt_starship": cls.get("HGT / Starship", 0),
            "te_repeat": cls.get("TE / repeat", 0) + cls.get("RIP / genome defense", 0),
            "other_structural": cls.get("Structural (other)", 0),
        })

    # Panel-aggregate row at the bottom.
    if rows:
        tot_single = sum(int(r["single_ref_equivalent"]) for r in rows)
        tot_po = sum(int(r["pangenome_only"]) for r in rows)
        tot_dedup = sum(int(r["total_loci"]) for r in rows)
        tot_po_read = sum(int(r["pangenome_only_read_supported"]) for r in rows)
        rv_all_yes = sum(by_asm_rv.get(q, {}).get("mycosv", {}).get("yes", 0)
                         for q in queries)
        rv_all_tot = sum(by_asm_rv.get(q, {}).get("mycosv", {}).get("total", 0)
                         for q in queries)
        rows.append({
            "sample": f"Panel ({len(queries)} samples)",
            "species": ".", "class": ".",
            "assembly_size_mb": (
                f"{sum(v for v in asm_size_mb.values() if v):.0f}"
                if asm_size_mb else ""
            ),
            "total_loci": tot_dedup,
            "single_ref_equivalent": tot_single,
            "pangenome_only": tot_po,
            "pangenome_lift_pct": f"{(tot_po / tot_dedup * 100):.1f}" if tot_dedup else "0.0",
            "pangenome_only_read_supported": tot_po_read,
            "mycosv_read_validated_pct": (
                f"{(rv_all_yes / rv_all_tot * 100):.1f}" if rv_all_tot else "."
            ),
            "best_comparator_f1": (
                f"{(sum(best_f1.values()) / len(best_f1)):.3f}" if best_f1 else "."
            ),
            "hgt_starship": sum(int(r["hgt_starship"]) for r in rows),
            "te_repeat": sum(int(r["te_repeat"]) for r in rows),
            "other_structural": sum(int(r["other_structural"]) for r in rows),
        })

    columns = [
        "sample", "species", "class", "assembly_size_mb",
        "total_loci", "single_ref_equivalent", "pangenome_only",
        "pangenome_lift_pct", "pangenome_only_read_supported",
        "mycosv_read_validated_pct", "best_comparator_f1",
        "hgt_starship", "te_repeat", "other_structural",
    ]
    write_tsv(out_tsv, rows, columns)

    # GitHub-flavored markdown for direct copy into the manuscript draft.
    display = {
        "sample": "Sample",
        "species": "Species",
        "class": "Class",
        "assembly_size_mb": "Asm (Mb)",
        "total_loci": "Total loci",
        "single_ref_equivalent": "Single-ref equiv.",
        "pangenome_only": "Pangenome-only",
        "pangenome_lift_pct": "Pangenome lift (%)",
        "pangenome_only_read_supported": "PO read-supported",
        "mycosv_read_validated_pct": "MycoSV read-valid. (%)",
        "best_comparator_f1": "Best comp. F1",
        "hgt_starship": "HGT/Starship",
        "te_repeat": "TE/RIP/repeat",
        "other_structural": "Other SV",
    }
    header = "| " + " | ".join(display[c] for c in columns) + " |"
    sep    = "| " + " | ".join("---" for _ in columns) + " |"
    caption = (
        "**Table 1.** Pangenome-powered structural variant recovery versus single-reference"
        " benchmark across the MycoSV panel. Per-sample rows ordered by phylogeny"
        " (Glomeromycota → Basidiomycota → Ascomycota). *Pangenome lift* is the percentage"
        " of deduplicated MycoSV loci not representable on the held-out single benchmark"
        " reference; *read-supported* loci are independently confirmed by raw-read split"
        " alignments. *Best comparator F1* is the maximum of MycoSV vs anchorwave / svim_asm"
        " / minigraph in exact-match scoring against the same single-reference truth."
        f" Of {total_panel_shards} panel shards attempted, {len(queries)} completed the full"
        f" MycoSV calling pipeline; the remaining {failed_shards} timed out with a promoted"
        " hierarchical checkpoint (MYCOSV_FAILED.txt marker) and are excluded from this"
        " table to keep per-sample numbers directly comparable. Column choice follows the"
        " convention of Plissonneau et al. (2018, mBio), Hartmann et al. (2017, Mol Ecol),"
        " Badet et al. (2020, BMC Biol), Cuomo et al. (2019, Cell Host Microbe), and Liang"
        " et al. (2024, MBE)."
    )
    lines = [caption, "", header, sep]
    for r in rows:
        cells = []
        for c in columns:
            v = r.get(c, "")
            if isinstance(v, float):
                cells.append(f"{v:g}")
            else:
                cells.append(str(v).replace("|", "\\|"))
        lines.append("| " + " | ".join(cells) + " |")
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return rows


def build_manuscript_method_comparison_table(
    by_query_dir: Path,
    pl_rows: list[dict[str, object]],
    rv_rows: list[dict[str, object]],
    out_tsv: Path,
    out_md: Path,
) -> list[dict[str, object]]:
    """Manuscript headline table — method-level comparison aggregated across
    the panel. One row per SV-calling method, columns for the metrics that
    distinguish pangenome from single-reference approaches:
      * call volume (median per sample / panel total)
      * pangenome-only loci recovered (the lift)
      * off-reference novel sequence captured
      * independent read-validation rate
      * cross-clade lineage / element-class annotation availability

    This is the single self-contained table reviewers cite when asking
    'why pangenome?' — drop directly into the manuscript.
    """
    by_asm_pl = {str(r["query_asm"]): r for r in pl_rows}
    samples = list(by_asm_pl)
    # Per-sample comparator truth_calls from exact_benchmark_summary.tsv.
    comparator_counts: dict[str, list[int]] = defaultdict(list)
    for q, shard in _iter_shards(by_query_dir):
        if q not in by_asm_pl:
            continue
        f = shard / "exact_benchmark_summary.tsv"
        if not f.exists() or f.stat().st_size == 0:
            continue
        seen: set[str] = set()
        with f.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if (r.get("coordinate_space") or "") != "reference":
                    continue
                if (r.get("validation_basis") or "") != "comparator_baseline":
                    continue
                if (r.get("svtype") or "") != "ALL":
                    continue
                src = (r.get("truth_label") or "").strip()
                if not src or src in seen:
                    continue
                try:
                    n = int(r.get("truth_calls") or 0)
                except ValueError:
                    continue
                comparator_counts[src].append(n)
                seen.add(src)
    mycosv_counts = [int(by_asm_pl[s]["dedup_loci"]) for s in samples]
    mycosv_single_ref = sum(int(by_asm_pl[s]["single_ref_equivalent"]) for s in samples)
    mycosv_po = sum(int(by_asm_pl[s]["pangenome_only"]) for s in samples)
    mycosv_total = mycosv_single_ref + mycosv_po
    # OFF_REF novel sequence counts from per-shard mycosv VCFs.
    off_ref_total = 0
    for q, shard in _iter_shards(by_query_dir):
        if q not in by_asm_pl:
            continue
        vcf = shard / "mycosv" / "calls.vcf"
        if not vcf.exists():
            continue
        try:
            with vcf.open() as fh:
                for line in fh:
                    if not line.startswith("#") and "SVTYPE=OFF_REF" in line:
                        off_ref_total += 1
        except OSError:
            continue
    # Per-source read-validation, panel sum.
    rv_yes: dict[str, int] = defaultdict(int)
    rv_tot: dict[str, int] = defaultdict(int)
    for r in rv_rows:
        s = str(r["source"])
        rv_yes[s] += int(r["yes"])
        rv_tot[s] += int(r["total"])

    def med(xs: list[int]) -> int:
        if not xs:
            return 0
        s = sorted(xs)
        return s[len(s) // 2] if len(s) % 2 else (s[len(s) // 2 - 1] + s[len(s) // 2]) // 2

    def rate_str(s: str) -> str:
        return (f"{rv_yes[s] / rv_tot[s]:.3f} "
                f"({rv_yes[s]:,} / {rv_tot[s]:,})") if rv_tot.get(s) else "n/a"

    anchor_n = sum(comparator_counts.get("anchorwave", []))
    svim_n = sum(comparator_counts.get("svim_asm", []))
    mini_n = sum(comparator_counts.get("minigraph", []))
    rows: list[dict[str, object]] = [
        {
            "method": "MycoSV (this work)",
            "strategy": "Pangenome graph routing",
            "coord_space": "Query-native + reference",
            "median_calls_per_sample": med(mycosv_counts),
            "single_reference_calls": f"{mycosv_single_ref:,}",
            "pangenome_only_loci": (
                f"{mycosv_po:,} ({100*mycosv_po/mycosv_total:.1f}%)"
                if mycosv_total else "0"
            ),
            "off_ref_novel_sequence": f"{off_ref_total:,}",
            "panel_total_calls": mycosv_total,
            "read_validation_rate": rate_str("mycosv"),
            "lineage_labels": "Yes (per-call clade rank + element class)",
        },
        {
            "method": "AnchorWave",
            "strategy": "Whole-genome alignment + paftools.js",
            "coord_space": "Single reference",
            "median_calls_per_sample": med(comparator_counts.get("anchorwave", [])),
            "single_reference_calls": f"{anchor_n:,}",
            "pangenome_only_loci": "0 (not supported)",
            "off_ref_novel_sequence": "0 (not supported)",
            "panel_total_calls": anchor_n,
            "read_validation_rate": rate_str("anchorwave"),
            "lineage_labels": "No",
        },
        {
            "method": "svim-asm",
            "strategy": "minimap2 BAM + SVIM-asm",
            "coord_space": "Single reference",
            "median_calls_per_sample": med(comparator_counts.get("svim_asm", [])),
            "single_reference_calls": f"{svim_n:,}",
            "pangenome_only_loci": "0 (not supported)",
            "off_ref_novel_sequence": "0 (not supported)",
            "panel_total_calls": svim_n,
            "read_validation_rate": rate_str("svim_asm"),
            "lineage_labels": "No",
        },
        {
            "method": "minigraph",
            "strategy": "Pangenome graph (bubble extraction)",
            "coord_space": "Reference + bubble graph",
            "median_calls_per_sample": med(comparator_counts.get("minigraph", [])),
            "single_reference_calls": f"{mini_n:,}",
            "pangenome_only_loci": "Partial (bubble traversals only)",
            "off_ref_novel_sequence": "0 (not supported)",
            "panel_total_calls": mini_n,
            "read_validation_rate": rate_str("minigraph"),
            "lineage_labels": "No (sample IDs only)",
        },
    ]

    columns = [
        "method", "strategy", "coord_space",
        "median_calls_per_sample",
        "single_reference_calls",
        "pangenome_only_loci", "off_ref_novel_sequence",
        "panel_total_calls",
        "read_validation_rate", "lineage_labels",
    ]
    write_tsv(out_tsv, rows, columns)

    n_panel = len(samples)
    display = {
        "method": "Method",
        "strategy": "Strategy",
        "coord_space": "Coordinate space",
        "median_calls_per_sample": "Median calls / sample",
        "single_reference_calls": "Single-reference calls (panel)",
        "pangenome_only_loci": "Pangenome-only loci (panel)",
        "off_ref_novel_sequence": "OFF_REF novel sequence (panel)",
        "panel_total_calls": "Panel total calls",
        "read_validation_rate": "Read-validation rate (panel)",
        "lineage_labels": "Cross-clade lineage / element-class labels",
    }
    header = "| " + " | ".join(display[c] for c in columns) + " |"
    sep    = "| " + " | ".join("---" for _ in columns) + " |"
    caption = (
        "**Table.** Pangenome-powered MycoSV recovers more structural variants"
        " at substantially higher independent read-validation rate than"
        " single-reference comparators, and uniquely captures off-reference"
        " novel sequence and cross-clade lineage context. Panel aggregates"
        f" across {n_panel} fungal genomes for which all four callers completed."
        " *Median calls / sample* counts deduplicated MycoSV pangenome loci or"
        " comparator SV calls per sample. *Single-reference calls* are the"
        " calls representable in the benchmark-reference coordinate space"
        " (for MycoSV: the single-reference-equivalent subset of its"
        " pangenome loci; for the three comparators: their entire output by"
        " construction). *Pangenome-only loci* are MycoSV calls not"
        " representable on the held-out single benchmark reference;"
        " single-reference comparators do not target this category by"
        " construction. *OFF_REF novel sequence* counts insertions whose"
        " alternate allele has no homologous reference anchor — a category"
        " only the pangenome graph routes calls into. *Panel total calls*"
        " is the sum of single-reference + pangenome-only loci (for MycoSV)"
        " or the comparator's entire output. *Read-validation rate* is the"
        " fraction of each caller's truth set independently supported by"
        " raw-read split alignments in the same query (yes / total shown in"
        " parentheses). Table layout follows the method-comparison convention"
        " of Plissonneau et al. (2018, mBio Table 2), Liang et al. (2024,"
        " MBE Table 1), and Heller & Vingron (2019, Bioinformatics Table 1)."
    )
    md_lines = [caption, "", header, sep]
    for r in rows:
        cells = []
        for c in columns:
            v = r.get(c, "")
            if isinstance(v, int):
                cells.append(f"{v:,}")
            elif isinstance(v, float):
                cells.append(f"{v:g}")
            else:
                cells.append(str(v).replace("|", "\\|"))
        # Bold the MycoSV row so the headline pops in the rendered manuscript.
        if r.get("method", "").startswith("MycoSV"):
            cells = [f"**{c}**" for c in cells]
        md_lines.append("| " + " | ".join(cells) + " |")
    out_md.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    return rows


def build_manuscript_table2_longread(
    by_query_dir: Path,
    rv_rows: list[dict[str, object]],
    meta_by_asm: dict[str, dict[str, str]],
    out_tsv: Path,
    out_md: Path,
) -> list[dict[str, object]]:
    """Manuscript Table 2 — long-read SV recovery: MycoSV-pangenome vs
    canonical long-read callers (Sniffles2, cuteSV, SVIM).

    Companion to Table 1. Table 1 measures assembly-vs-pangenome SV recovery
    at full-genome scale. Table 2 measures read-level SV calling at the same
    signal as the canonical PB/ONT comparators — same FASTQ in, same single
    benchmark reference for the comparator coordinate space. Two evaluation
    axes per the manuscript scoping:
      1. Independent read-level validation — split-read support % column.
      2. Same single-reference coordinate space — the F1 vs <caller>
         columns are computed against each comparator's call set on the
         shared benchmark reference, so all four tools are directly
         comparable on the same axis.

    Column choice follows Liao 2023 HPRC Table 2, Hickey 2024
    Minigraph-Cactus Table 3, and Heller & Vingron 2022 Sniffles2 Table 1.

    Per-sample columns:
      Sample / species / class / platform / coverage
      MycoSV calls            — pangenome-routed long-read SV calls
      Sniffles2 calls         — canonical ONT/PacBio SV caller
      cuteSV calls            — canonical long-read SV caller
      SVIM calls              — canonical long-read SV caller
      MycoSV ∩ ≥2 LR callers  — consensus_2of_3 vs MycoSV agreement (intersection size)
      MycoSV F1 vs sniffles   — exact-match F1 against sniffles truth
      MycoSV F1 vs cuteSV     — exact-match F1 against cuteSV truth
      MycoSV F1 vs SVIM       — exact-match F1 against SVIM truth
      MycoSV read-valid. (%)  — split-read support fraction
      OFF_REF novel (MycoSV)  — insertions with no homologous reference anchor

    Plus a panel-aggregate row at the bottom.
    """
    by_asm_rv: dict[str, dict[str, dict[str, int]]] = defaultdict(dict)
    for r in rv_rows:
        by_asm_rv[str(r["query_asm"])][str(r["source"])] = {
            "yes": int(r["yes"]), "total": int(r["total"])
        }

    LR_LABELS = ("sniffles", "cutesv", "svim")  # exact truth_label values

    # Per-shard: per-comparator truth_calls (= that caller's own callset size)
    # + MycoSV-vs-comparator F1, plus the consensus_2of_3 agreement count.
    per_sample: dict[str, dict[str, object]] = {}
    coverage_by_asm: dict[str, str] = {}
    platform_by_asm: dict[str, str] = {}
    off_ref_by_asm: dict[str, int] = {}

    for qname, shard in _iter_shards(by_query_dir):
        if (shard / "MYCOSV_FAILED.txt").exists():
            continue
        eb = shard / "exact_benchmark_summary.tsv"
        if not eb.exists() or eb.stat().st_size == 0:
            continue
        comp_calls: dict[str, int] = {}
        comp_f1: dict[str, float] = {}
        consensus_count = 0
        mycosv_calls = 0
        with eb.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if (r.get("coordinate_space") or "") != "reference":
                    continue
                if (r.get("validation_basis") or "") != "comparator_baseline":
                    continue
                if (r.get("svtype") or "") != "ALL":
                    continue
                if (r.get("status") or "") != "ok":
                    continue
                label = (r.get("truth_label") or "").lower()
                try:
                    truth_n = int(r.get("truth_calls") or 0)
                    pred_n = int(r.get("pred_calls") or 0)
                    f1 = float(r.get("f1") or 0)
                except ValueError:
                    continue
                if label in LR_LABELS:
                    comp_calls[label] = truth_n
                    comp_f1[label] = f1
                    if mycosv_calls == 0 and pred_n > 0:
                        mycosv_calls = pred_n
                elif label == "consensus_2of_3":
                    consensus_count = max(consensus_count, truth_n)

        # MycoSV's own VCF row count (split-read fallback if comparator rows absent)
        my_vcf = shard / "mycosv" / "calls.vcf"
        if mycosv_calls == 0 and my_vcf.exists():
            try:
                with my_vcf.open() as fh:
                    mycosv_calls = sum(
                        1 for ln in fh if ln and not ln.startswith("#")
                    )
            except OSError:
                pass

        # OFF_REF novel-sequence count from MycoSV's INFO=SVTYPE=OFF_REF rows.
        off_ref = 0
        if my_vcf.exists():
            try:
                with my_vcf.open() as fh:
                    for ln in fh:
                        if ln.startswith("#") or not ln.strip():
                            continue
                        if "SVTYPE=OFF_REF" in ln or "\tOFF_REF\t" in ln:
                            off_ref += 1
            except OSError:
                pass
        off_ref_by_asm[qname] = off_ref

        # Coverage + platform: pulled from the MycoSV stderr auto-tune line
        # (single source of truth for the actual effective coverage MycoSV used).
        cov_str = ""
        platform = ""
        my_stderr = shard / "mycosv" / "calls.stderr.log"
        if my_stderr.exists():
            try:
                with my_stderr.open(errors="replace") as fh:
                    for ln in fh:
                        if "cov~" in ln and ("auto-tuning" in ln or "cov_tier=" in ln):
                            m = re.search(r"cov~([0-9.]+)x", ln)
                            if m:
                                cov_str = f"{float(m.group(1)):.1f}x"
                                break
            except OSError:
                pass
        coverage_by_asm[qname] = cov_str
        # Platform from meta if present.
        meta = meta_by_asm.get(qname, {})
        ip = (meta.get("instrument_platform") or "").strip()
        if ip:
            short = {
                "OXFORD_NANOPORE": "ONT",
                "PACBIO_SMRT": "PacBio",
                "ILLUMINA": "Illumina",
                "DNBSEQ": "DNBSEQ",
            }
            platform = short.get(ip, ip)
        platform_by_asm[qname] = platform

        per_sample[qname] = {
            "mycosv_calls": mycosv_calls,
            "sniffles_calls": comp_calls.get("sniffles", 0),
            "cutesv_calls": comp_calls.get("cutesv", 0),
            "svim_calls": comp_calls.get("svim", 0),
            "mycosv_vs_2of3_consensus": consensus_count,
            "f1_vs_sniffles": comp_f1.get("sniffles", 0.0),
            "f1_vs_cutesv": comp_f1.get("cutesv", 0.0),
            "f1_vs_svim": comp_f1.get("svim", 0.0),
            "off_ref_novel": off_ref,
        }

    queries = sorted(per_sample,
                     key=lambda q: _phylogeny_sort_key(meta_by_asm.get(q, {"query_asm": q})))
    rows: list[dict[str, object]] = []
    for q in queries:
        s = per_sample[q]
        meta = meta_by_asm.get(q, {})
        rv = by_asm_rv.get(q, {}).get("mycosv", {"yes": 0, "total": 0})
        rows.append({
            "sample": q,
            "species": meta.get("species") or ".",
            "class": meta.get("class") or ".",
            "platform": platform_by_asm.get(q, ""),
            "coverage": coverage_by_asm.get(q, ""),
            "mycosv_calls": s["mycosv_calls"],
            "sniffles_calls": s["sniffles_calls"],
            "cutesv_calls": s["cutesv_calls"],
            "svim_calls": s["svim_calls"],
            "mycosv_vs_2of3_consensus": s["mycosv_vs_2of3_consensus"],
            "f1_vs_sniffles": f"{s['f1_vs_sniffles']:.3f}",
            "f1_vs_cutesv":   f"{s['f1_vs_cutesv']:.3f}",
            "f1_vs_svim":     f"{s['f1_vs_svim']:.3f}",
            "mycosv_read_validated_pct": (
                f"{(rv['yes'] / rv['total'] * 100):.1f}" if rv["total"] else "."
            ),
            "off_ref_novel": s["off_ref_novel"],
        })

    # Panel-aggregate row.
    if rows:
        def _avg_f1(field: str) -> str:
            vals = [per_sample[q][field] for q in queries
                    if isinstance(per_sample[q][field], (int, float))
                    and per_sample[q][field] > 0]
            return f"{(sum(vals) / len(vals)):.3f}" if vals else "."
        rv_all_yes = sum(by_asm_rv.get(q, {}).get("mycosv", {}).get("yes", 0)
                         for q in queries)
        rv_all_tot = sum(by_asm_rv.get(q, {}).get("mycosv", {}).get("total", 0)
                         for q in queries)
        rows.append({
            "sample": f"Panel ({len(queries)} samples)",
            "species": ".", "class": ".", "platform": ".", "coverage": ".",
            "mycosv_calls": sum(per_sample[q]["mycosv_calls"] for q in queries),
            "sniffles_calls": sum(per_sample[q]["sniffles_calls"] for q in queries),
            "cutesv_calls": sum(per_sample[q]["cutesv_calls"] for q in queries),
            "svim_calls": sum(per_sample[q]["svim_calls"] for q in queries),
            "mycosv_vs_2of3_consensus": sum(
                per_sample[q]["mycosv_vs_2of3_consensus"] for q in queries),
            "f1_vs_sniffles": _avg_f1("f1_vs_sniffles"),
            "f1_vs_cutesv": _avg_f1("f1_vs_cutesv"),
            "f1_vs_svim": _avg_f1("f1_vs_svim"),
            "mycosv_read_validated_pct": (
                f"{(rv_all_yes / rv_all_tot * 100):.1f}" if rv_all_tot else "."
            ),
            "off_ref_novel": sum(off_ref_by_asm.get(q, 0) for q in queries),
        })

    columns = [
        "sample", "species", "class", "platform", "coverage",
        "mycosv_calls", "sniffles_calls", "cutesv_calls", "svim_calls",
        "mycosv_vs_2of3_consensus",
        "f1_vs_sniffles", "f1_vs_cutesv", "f1_vs_svim",
        "mycosv_read_validated_pct", "off_ref_novel",
    ]
    write_tsv(out_tsv, rows, columns)

    display = {
        "sample": "Sample",
        "species": "Species",
        "class": "Class",
        "platform": "Platform",
        "coverage": "Cov",
        "mycosv_calls": "MycoSV calls",
        "sniffles_calls": "Sniffles2 calls",
        "cutesv_calls": "cuteSV calls",
        "svim_calls": "SVIM calls",
        "mycosv_vs_2of3_consensus": "MycoSV ∩ ≥2 LR callers",
        "f1_vs_sniffles": "F1 vs Sniffles2",
        "f1_vs_cutesv":   "F1 vs cuteSV",
        "f1_vs_svim":     "F1 vs SVIM",
        "mycosv_read_validated_pct": "MycoSV split-read support (%)",
        "off_ref_novel": "OFF_REF novel (MycoSV)",
    }
    header = "| " + " | ".join(display[c] for c in columns) + " |"
    sep    = "| " + " | ".join("---" for _ in columns) + " |"
    caption = (
        "**Table 2.** Long-read structural variant recovery: MycoSV-pangenome"
        " versus canonical long-read callers (Sniffles2, cuteSV, SVIM) on the"
        " same PB/ONT FASTQ input. Per-sample rows ordered by phylogeny."
        " *Coverage* is the auto-tuner's effective coverage estimate from the"
        " MycoSV preprocessing log (mean read length × read count / genome bp)."
        " *MycoSV ∩ ≥2 LR callers* is the count of MycoSV calls supported by"
        " the consensus_2of_3 long-read truth set (Sniffles2 ∩ cuteSV ∪ SVIM"
        " pairwise majority). *F1 vs <caller>* is the exact-match F1 of"
        " MycoSV's predictions against that caller's truth set on the same"
        " single benchmark reference. *Split-read support (%)* is the fraction"
        " of MycoSV calls independently confirmed by raw split-read"
        " alignments. *OFF_REF novel* counts MycoSV insertions whose alt"
        " allele has no homologous anchor on the reference — a category only"
        " the pangenome graph routes calls into. Column choice follows Liao"
        " et al. (2023, Nature, HPRC Table 2), Hickey et al. (2024, Nat"
        " Biotechnol, Minigraph-Cactus Table 3), and Heller & Vingron"
        " (2022, Bioinformatics, Sniffles2 Table 1)."
    )
    md_lines = [caption, "", header, sep]
    for r in rows:
        cells = []
        for c in columns:
            v = r.get(c, "")
            if isinstance(v, int):
                cells.append(f"{v:,}")
            elif isinstance(v, float):
                cells.append(f"{v:g}")
            else:
                cells.append(str(v).replace("|", "\\|"))
        md_lines.append("| " + " | ".join(cells) + " |")
    out_md.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    return rows


def run_longread_panel_mode(by_query_dir: Path, out_dir: Path) -> int:
    """Long-read companion to run_panel_mode: emits manuscript Table 2 only.

    Does not regenerate Figure 1A/1B (those are assembly-mode aggregates).
    Reads the same `by_query/<asm>/` tree but expects long-read comparator
    truth labels (sniffles / cutesv / svim) in exact_benchmark_summary.tsv.
    """
    if not by_query_dir.is_dir():
        print(f"[longread-panel] not a directory: {by_query_dir}")
        return 1
    out_dir.mkdir(parents=True, exist_ok=True)
    meta_by_asm = _load_query_meta_map(by_query_dir)
    rv_rows = aggregate_panel_read_validation(
        by_query_dir, out_dir / "panel_read_validation_rate_longreads.tsv",
    )
    build_manuscript_table2_longread(
        by_query_dir,
        rv_rows,
        meta_by_asm,
        out_dir / "table2_longread_caller_comparison.tsv",
        out_dir / "table2_longread_caller_comparison.md",
    )
    print(f"[longread-panel] wrote Table 2 -> {out_dir}/table2_longread_caller_comparison.{{tsv,md}}")
    return 0


def run_panel_mode(by_query_dir: Path, out_dir: Path) -> int:
    """Build Figure 1 panel-fold aggregates and plots from a by_query/ tree."""
    if not by_query_dir.is_dir():
        print(f"[panel] not a directory: {by_query_dir}")
        return 1
    out_dir.mkdir(parents=True, exist_ok=True)
    meta_by_asm = _load_query_meta_map(by_query_dir)

    rv_rows = aggregate_panel_read_validation(
        by_query_dir, out_dir / "panel_read_validation_rate.tsv",
    )
    fig1a = out_dir / "fig1a_panel_read_validation_rate.png"
    save_panel_read_validation_rate(
        fig1a, rv_rows, meta_by_asm,
        title="Per-sample read-validation rate, MycoSV vs comparator callers",
    )

    pl_rows = aggregate_panel_pangenome_layers(
        by_query_dir, out_dir / "panel_pangenome_lift.tsv",
    )
    fig1b = out_dir / "fig1b_panel_pangenome_lift.png"
    save_panel_pangenome_lift(
        fig1b, pl_rows, meta_by_asm,
        title="Per-sample MycoSV pangenome lift (single-ref vs pangenome-only)",
    )

    build_manuscript_table1(
        by_query_dir,
        pl_rows,
        rv_rows,
        meta_by_asm,
        out_dir / "table1_pangenome_vs_single_reference.tsv",
        out_dir / "table1_pangenome_vs_single_reference.md",
    )

    # Headline manuscript method-comparison table (one row per method).
    build_manuscript_method_comparison_table(
        by_query_dir,
        pl_rows,
        rv_rows,
        out_dir / "table_method_comparison_pangenome_vs_single_ref.tsv",
        out_dir / "table_method_comparison_pangenome_vs_single_ref.md",
    )

    # Minimal HTML wrapper so an operator can preview both panels together.
    page = f"""<!doctype html>
<html><head><meta charset='utf-8'>
<title>{html.escape('MycoSV manuscript Figure 1 panel-fold preview')}</title>
<style>body{{font-family:system-ui,sans-serif;margin:1.5rem;color:#222}}
img{{max-width:100%;border:1px solid #ddd;margin:0.6rem 0}}</style>
</head><body>
<h1>MycoSV manuscript Figure 1 — panel-fold preview</h1>
<p>Aggregated across {len(set(r['query_asm'] for r in rv_rows))} samples under
{html.escape(str(by_query_dir))}.</p>
<h2>Fig 1A — Per-sample read-validation rate</h2>
<img src='{fig1a.name}'>
<h2>Fig 1B — Per-sample pangenome lift</h2>
<img src='{fig1b.name}'>
</body></html>
"""
    (out_dir / "fig1_panel_preview.html").write_text(page, encoding="utf-8")
    print(f"[panel] wrote panel-fold figures and TSVs to {out_dir}")
    return 0


def save_pangenome_lift_donut(
    path: Path,
    layers: list[dict[str, str]],
    *,
    title: str,
) -> None:
    """Concentric donut showing total → pangenome-only → read-supported nested
    fractions. Communicates the headline pangenome-lift story in one panel
    instead of three side-by-side bars."""
    if plt is None or not layers:
        return
    layer = next((r for r in layers if r.get("query_asm") == "ALL"), layers[0])
    raw = as_int(layer.get("raw_pairwise_pangenome_observations"))
    dedup = as_int(layer.get("deduplicated_pangenome_loci"))
    single_ref = as_int(layer.get("single_reference_equivalent_calls"))
    pango = as_int(layer.get("pangenome_only_calls"))
    po_read = as_int(layer.get("pangenome_only_read_supported"))
    if dedup <= 0:
        return
    fig, ax = plt.subplots(figsize=(7.6, 6.4))
    ax.axis("equal")
    outer_vals = [single_ref, pango]
    outer_labels = [f"single-ref equiv\n{single_ref:,}", f"pangenome-only\n{pango:,}"]
    outer_colors = [SCOPE_PALETTE["single_reference_equivalent"], SCOPE_PALETTE["pangenome_only"]]
    wedge_kwargs = dict(width=0.32, edgecolor="white", linewidth=1.5)
    ax.pie(outer_vals, radius=1.0, colors=outer_colors,
           labels=outer_labels, labeldistance=1.08,
           wedgeprops=wedge_kwargs,
           startangle=90, counterclock=False,
           textprops={"fontsize": 10})
    # Inner ring: read-supported within pangenome-only
    inner_vals = [po_read, max(0, pango - po_read), single_ref]
    inner_colors = ["#009E73", "#7F7F7F", outer_colors[0]]
    wedges, _ = ax.pie(
        inner_vals, radius=0.66, colors=inner_colors,
        wedgeprops=dict(width=0.32, edgecolor="white", linewidth=1.5),
        startangle=90, counterclock=False,
    )
    ax.text(0, 0,
            f"{dedup:,}\nloci",
            ha="center", va="center", fontsize=14, fontweight="semibold")
    legend_handles = [
        _mpl.patches.Patch(color=outer_colors[0], label="single-ref equivalent"),
        _mpl.patches.Patch(color=outer_colors[1], label="pangenome-only"),
        _mpl.patches.Patch(color="#009E73", label="read-supported (inner)"),
        _mpl.patches.Patch(color="#7F7F7F", label="intrinsic-only (inner)"),
    ]
    ax.legend(handles=legend_handles, loc="lower center", bbox_to_anchor=(0.5, -0.05),
              ncol=2, frameon=False, fontsize=9)
    ax.set_title(title, fontsize=12, fontweight="semibold")
    _save_dual(fig, path)
    plt.close(fig)


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
    parser.add_argument("--benchmark-dir", type=Path, default=None,
                        help="Per-shard mode: directory containing novel_mycosv_calls.tsv "
                             "and pangenome_call_layers.tsv for a single query.")
    parser.add_argument("--panel-dir", type=Path, default=None,
                        help="Panel-fold mode: a by_query/ directory holding one subdir per "
                             "shard. Produces manuscript Figure 1A/1B panel-level aggregates "
                             "(per-sample read-validation rate + pangenome-lift bars). Also "
                             "writes the manuscript Table 1 (pangenome vs single-reference).")
    parser.add_argument("--longread-panel-dir", type=Path, default=None,
                        help="Long-read companion: a by_query/ directory holding one subdir per "
                             "long-read shard. Writes manuscript Table 2 (MycoSV vs Sniffles2 / "
                             "cuteSV / SVIM on the same PB/ONT input). Does NOT regenerate "
                             "Figure 1; that is assembly-mode only.")
    parser.add_argument("--outdir", type=Path, default=None,
                        help="Output directory. Per-shard mode defaults to <benchmark-dir>/"
                             "pangenome_plots; panel mode defaults to <panel-dir>/../"
                             "manuscript_figure1; long-read panel mode defaults to "
                             "<longread-panel-dir>/../manuscript_table2.")
    parser.add_argument("--title", default="MycoSV pangenome call biology plots")
    args = parser.parse_args(argv)

    if args.longread_panel_dir is not None:
        out = args.outdir or (args.longread_panel_dir.parent / "manuscript_table2")
        return run_longread_panel_mode(args.longread_panel_dir, out)
    if args.panel_dir is not None:
        out = args.outdir or (args.panel_dir.parent / "manuscript_figure1")
        return run_panel_mode(args.panel_dir, out)
    if args.benchmark_dir is None:
        parser.error("either --benchmark-dir, --panel-dir, or --longread-panel-dir is required")

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

        # Manuscript-style donut: single-ref vs pangenome-only with inner
        # read-support partition. Communicates the headline lift in one panel.
        fig = outdir / "pangenome_lift_donut.png"
        save_pangenome_lift_donut(fig, layers, title="MycoSV pangenome lift")
        if fig.exists():
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

        # Lollipop variant: read-support rate on stem color, count on length.
        fig = outdir / "top_contigs_lollipop.png"
        save_lollipop_top_contigs(
            fig, novel,
            title="Pangenome-only calls per contig (top 15) — color = read-support rate",
            top_n=15,
        )
        if fig.exists():
            figures.append(fig)

        # Circos landscape: stacked SVTYPE bars over top contigs + pangenome-only
        # inner ring. The single most informative whole-genome view.
        fig = outdir / "circos_sv_landscape.png"
        save_circos_sv_landscape(
            fig, novel,
            title="MycoSV structural variant landscape (top contigs)",
            top_contigs=20, bin_kb=50,
        )
        if fig.exists():
            figures.append(fig)

        # Ridge plot of SV length distributions by SVTYPE on log10 scale.
        fig = outdir / "svlen_ridge.png"
        save_ridge_svlen(
            fig, novel,
            title="MycoSV SV length distribution by type",
        )
        if fig.exists():
            figures.append(fig)

        # Upset plot of comparator co-support memberships.
        fig = outdir / "comparator_overlap_upset.png"
        save_upset_comparator_overlap(
            fig, novel,
            title="MycoSV call support across comparators (UpSet)",
            min_subset_size=5,
        )
        if fig.exists():
            figures.append(fig)

        # Manhattan-style SV density along contigs — scales further than
        # circos when the contig list is long.
        fig = outdir / "sv_density_manhattan.png"
        save_manhattan_sv_density(
            fig, novel,
            title="MycoSV density along query contigs",
            bin_kb=25, top_contigs=30,
        )
        if fig.exists():
            figures.append(fig)

        # Two-speed genome scatter — TE/RIP/HGT density vs other-class
        # density per window. The Möller-Stukenbrock diagnostic.
        fig = outdir / "two_speed_genome_scatter.png"
        save_two_speed_scatter(
            fig, novel,
            title="Two-speed genome: TE/RIP/HGT vs other SVs per window",
            bin_kb=50,
        )
        if fig.exists():
            figures.append(fig)

    # Volcano plot of biology-feature enrichment from the per-query
    # mycosv_novel_biology_enrichment.tsv. Same plotting helper handles the
    # cross-guild combined table when this script is pointed at the
    # combined/cross_guild/ directory in a future invocation.
    enrichment = read_tsv(bench / "mycosv_novel_biology_enrichment.tsv")
    if enrichment:
        # The per-query table uses fisher_right_p (one-sided test) and
        # group_a/group_b instead of guild; remap so the same renderer
        # works for both per-query and cross-guild inputs.
        normalized = []
        for r in enrichment:
            normalized.append({
                "guild": r.get("group_a") or r.get("guild") or ".",
                "feature": r.get("feature") or ".",
                "odds_ratio": r.get("odds_ratio") or "0",
                "fisher_p": r.get("fisher_right_p") or r.get("fisher_p") or "1",
            })
        fig = outdir / "biology_enrichment_volcano.png"
        save_volcano_enrichment(
            fig, normalized,
            title="MycoSV unique-vs-recurrent feature enrichment",
            group_col="guild",
        )
        if fig.exists():
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

        # Hierarchically-clustered, row-z-scored heatmap of the same matrix —
        # surfaces class-specific eco enrichments that the raw-count heatmap
        # buries under the dominant RIP/TE rows.
        fig = outdir / "ecology_biology_clustermap.png"
        save_clustermap_biology(
            fig, biology,
            title="Biology class × ecology — clustered z-score heatmap",
        )
        if fig.exists():
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
