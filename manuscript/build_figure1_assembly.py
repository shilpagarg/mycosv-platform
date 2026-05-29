#!/usr/bin/env python3
"""Build a Nature Biotechnology-style assembly-mode Figure 1.

This composite is intentionally assembly-only. Long-read data are reserved for
Table 2, and the panel-200 run is reserved for Figure 2.
"""

from __future__ import annotations

import argparse
import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch


CALLER_ORDER = ["mycosv", "minigraph", "anchorwave", "svim_asm", "pggb", "cactus"]
CALLER_LABEL = {
    "mycosv": "MycoSV",
    "minigraph": "minigraph",
    "anchorwave": "AnchorWave",
    "svim_asm": "svim-asm",
    "pggb": "PGGB",
    "cactus": "Cactus",
}
CALLER_COLOR = {
    "mycosv": "#0072B2",
    "minigraph": "#009E73",
    "anchorwave": "#E69F00",
    "svim_asm": "#CC79A7",
    "pggb": "#D55E00",
    "cactus": "#56B4E9",
}
SVTYPE_COLOR = {
    "INS": "#E69F00",
    "DEL": "#56B4E9",
    "DUP": "#CC79A7",
    "INV": "#009E73",
    "TRA": "#D55E00",
    "OFF_REF": "#0072B2",
    "OTHER": "#7F7F7F",
}
SCOPE_COLOR = {
    "single": "#9E9E9E",
    "pangenome": "#0072B2",
    "supported": "#009E73",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        import csv
        return list(csv.DictReader(fh, delimiter="\t"))


def as_float(value: str | None, default: float = 0.0) -> float:
    try:
        if value in {None, "", "."}:
            return default
        return float(value)
    except ValueError:
        return default


def as_int(value: str | None, default: int = 0) -> int:
    return int(round(as_float(value, default)))


def info_dict(info: str) -> dict[str, str | bool]:
    out: dict[str, str | bool] = {}
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            key, val = item.split("=", 1)
            out[key] = val
        else:
            out[item] = True
    return out


def read_vcf_calls(path: Path, sample: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    if not path.exists():
        return rows
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, _vid, _ref, _alt, _qual, _flt, info = fields[:8]
            meta = info_dict(info)
            svtype = str(meta.get("SVTYPE", "OTHER")).upper()
            svlen = abs(as_int(str(meta.get("SVLEN", "0")), 0))
            if svlen == 0 and "END" in meta:
                svlen = abs(as_int(str(meta.get("END")), 0) - as_int(pos, 0))
            rows.append({
                "sample": sample,
                "chrom": chrom,
                "pos": as_int(pos, 0),
                "svtype": svtype if svtype in SVTYPE_COLOR else "OTHER",
                "svlen": svlen,
                "offref": "OFFREF" in meta or str(meta.get("OFF_REF_TIER", "")).startswith("NOVEL"),
            })
    return rows


def short_species(name: str) -> str:
    if not name or name == ".":
        return ""
    clean = name.replace("[", "").replace("]", "")
    parts = clean.split()
    if len(parts) >= 2:
        return f"{parts[0][0]}. {parts[1]}"
    return clean


def add_panel_label(ax, label: str) -> None:
    ax.text(-0.12, 1.08, label, transform=ax.transAxes, fontsize=13,
            fontweight="bold", va="bottom", ha="left")


def plot_validation(ax, rv_rows, sample_order, sample_label) -> None:
    by_sample = defaultdict(dict)
    for row in rv_rows:
        by_sample[row["query_asm"]][row["source"]] = row
    offsets = {"mycosv": -0.32, "minigraph": -0.19, "anchorwave": -0.06,
               "svim_asm": 0.06, "pggb": 0.19, "cactus": 0.32}
    for i, sample in enumerate(sample_order):
        for caller in CALLER_ORDER:
            row = by_sample.get(sample, {}).get(caller)
            if not row:
                continue
            y = as_float(row.get("rate"))
            lo = as_float(row.get("ci95_lo"), y)
            hi = as_float(row.get("ci95_hi"), y)
            x = i + offsets[caller]
            ax.errorbar(
                x, y, yerr=[[max(0, y - lo)], [max(0, hi - y)]],
                fmt="o", ms=5.2 if caller == "mycosv" else 4.4,
                color=CALLER_COLOR[caller], mec="white", mew=0.5,
                elinewidth=1.0, capsize=2.2, alpha=0.95, zorder=3,
            )
    ax.set_ylim(-0.04, 1.18)
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.set_ylabel("Read-validation rate")
    ax.set_xticks(range(len(sample_order)))
    ax.set_xticklabels([sample_label[s] for s in sample_order], rotation=30, ha="right",
                       fontstyle="italic")
    ax.tick_params(axis="x", pad=2)
    ax.axhline(1.0, color="#CCCCCC", lw=0.6, zorder=1)
    ax.grid(axis="y", color="#EAEAEA", lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title("Independent support across five assemblies", loc="left", pad=22)
    add_panel_label(ax, "a")


def plot_pangenome_lift(ax, lift_rows, sample_order, sample_label) -> None:
    rows = {r["query_asm"]: r for r in lift_rows}
    y_positions = list(range(len(sample_order)))
    max_total = 0
    pct_labels: list[tuple[int, float]] = []
    for y, sample in zip(y_positions, sample_order):
        row = rows[sample]
        single = as_int(row.get("single_ref_equivalent"))
        po = as_int(row.get("pangenome_only"))
        supported = as_int(row.get("pangenome_only_read_supported"))
        unsupported_po = max(0, po - supported)
        dedup = as_int(row.get("dedup_loci"))
        ax.barh(y, single, color=SCOPE_COLOR["single"], height=0.62)
        ax.barh(y, unsupported_po, left=single, color="#BFD7EA", height=0.62)
        ax.barh(y, supported, left=single + unsupported_po, color=SCOPE_COLOR["supported"], height=0.62)
        stacked_total = max(1, single + po)
        pct = 100 * po / max(1, dedup)
        max_total = max(max_total, stacked_total)
        pct_labels.append((y, pct))
    pct_x = max_total * 1.08
    for y, pct in pct_labels:
        ax.text(pct_x, y, f"{pct:.0f}%", va="center", ha="left", fontsize=8.2, color="#222222")
    ax.text(pct_x, -0.74, "pangenome-\nonly", va="bottom", ha="left", fontsize=7.6, color="#555555")
    ax.set_xlim(0, max_total * 1.22)
    ax.set_yticks(y_positions)
    ax.set_yticklabels([sample_label[s] for s in sample_order], fontstyle="italic")
    ax.invert_yaxis()
    ax.set_xlabel("MycoSV loci by representability class")
    ax.set_title("Pangenome-only loci dominate assembly calls", loc="left")
    ax.grid(axis="x", color="#EAEAEA", lw=0.6)
    ax.set_axisbelow(True)
    add_panel_label(ax, "b")


def add_clean_legends(fig, ax_a, ax_b) -> None:
    caller_handles = [
        mpl.lines.Line2D([0], [0], marker="o", linestyle="", color=CALLER_COLOR[c],
                         label=CALLER_LABEL[c], markersize=5.6)
        for c in CALLER_ORDER
    ]
    ax_a.legend(
        handles=caller_handles,
        ncols=len(CALLER_ORDER),
        loc="lower left",
        bbox_to_anchor=(0.0, 1.04),
        frameon=False,
        fontsize=8,
        columnspacing=1.1,
        handletextpad=0.32,
        borderaxespad=0,
    )

    scope_handles = [
        Patch(color=SCOPE_COLOR["single"], label="single-ref equivalent"),
        Patch(color=SCOPE_COLOR["supported"], label="pangenome-only, read-supported"),
        Patch(color="#BFD7EA", label="pangenome-only, intrinsic"),
    ]
    bbox = ax_b.get_position()
    fig.legend(
        handles=scope_handles,
        ncols=3,
        loc="upper left",
        bbox_to_anchor=(bbox.x0, bbox.y0 - 0.035),
        frameon=False,
        fontsize=7.8,
        handlelength=1.2,
        columnspacing=1.25,
        handletextpad=0.4,
        borderaxespad=0,
    )


def plot_svlen_violin(ax, calls, sample_order, sample_label) -> None:
    values = []
    for sample in sample_order:
        vals = [math.log10(max(50, int(r["svlen"]))) for r in calls if r["sample"] == sample and int(r["svlen"]) > 0]
        values.append(vals[:12000])
    parts = ax.violinplot(values, positions=range(len(sample_order)), widths=0.74,
                          showmeans=False, showmedians=True, showextrema=False)
    for body in parts["bodies"]:
        body.set_facecolor("#B7C9E2")
        body.set_edgecolor("#335C81")
        body.set_alpha(0.85)
        body.set_linewidth(0.7)
    parts["cmedians"].set_color("#1F1F1F")
    parts["cmedians"].set_linewidth(1.2)
    for i, vals in enumerate(values):
        if not vals:
            continue
        q1, q3 = sorted(vals)[int(0.25 * (len(vals)-1))], sorted(vals)[int(0.75 * (len(vals)-1))]
        ax.vlines(i, q1, q3, color="#333333", lw=2.1, zorder=4)
    ax.set_xticks(range(len(sample_order)))
    ax.set_xticklabels([sample_label[s] for s in sample_order], rotation=30, ha="right",
                       fontstyle="italic")
    ax.set_ylabel("SV length")
    ticks = [2, 3, 4, 5, 6]
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"$10^{t}$" for t in ticks])
    ax.set_ylim(1.7, 6.25)
    ax.grid(axis="y", color="#E5E5E5", lw=0.7)
    ax.set_title("SV-size spectra span focal and chromosome-scale events", loc="left")
    add_panel_label(ax, "c")


def plot_density(ax, calls, sample_order, sample_label, bin_size: int = 100_000) -> None:
    per_sample_bins: dict[str, Counter[tuple[str, int]]] = {}
    max_count = 1
    for sample in sample_order:
        counter: Counter[tuple[str, int]] = Counter()
        for r in calls:
            if r["sample"] != sample:
                continue
            counter[(str(r["chrom"]), int(r["pos"]) // bin_size)] += 1
        per_sample_bins[sample] = counter
        if counter:
            max_count = max(max_count, max(counter.values()))

    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "svdensity", ["#F7F7F7", "#B7D7E8", "#56B4E9", "#005A8D"]
    )
    norm = mpl.colors.Normalize(vmin=0, vmax=max_count)
    for y, sample in enumerate(sample_order):
        bins = per_sample_bins[sample]
        ranked = sorted(bins.items(), key=lambda kv: kv[1], reverse=True)[:220]
        ranked = sorted(ranked, key=lambda kv: (kv[0][0], kv[0][1]))
        for x, (_key, count) in enumerate(ranked):
            ax.scatter(x, y, s=10 + 18 * (count / max_count), color=cmap(norm(count)),
                       edgecolor="none", alpha=0.9)
        ax.text(-8, y, sample_label[sample], va="center", ha="right", fontsize=8,
                fontstyle="italic")
    ax.set_ylim(-0.6, len(sample_order) - 0.4)
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_xlabel("Contig-ordered 100 kb windows (top occupied windows per sample)")
    ax.set_title("Genome-wide SV-density hotspots", loc="left")
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = plt.colorbar(sm, ax=ax, fraction=0.032, pad=0.015)
    cbar.set_label("SVs / 100 kb", fontsize=8, labelpad=4)
    cbar.ax.tick_params(labelsize=8)
    add_panel_label(ax, "d")


def plot_biology_heatmap(ax, table1_rows, sample_order, sample_label) -> None:
    rows = {r["sample"]: r for r in table1_rows if not r["sample"].startswith("Panel")}
    cols = [
        ("hgt_starship", "HGT/\nStarship"),
        ("te_repeat", "TE/\nrepeat"),
        ("other_structural", "Other\nannotated"),
        ("unclassified", "Not classed\nhere"),
    ]
    matrix: list[list[float]] = []
    for sample in sample_order:
        row = rows[sample]
        total = max(1, as_int(row.get("total_svs") or row.get("total_loci")))
        hgt = as_int(row.get("hgt_starship"))
        te = as_int(row.get("te_repeat"))
        other = as_int(row.get("other_sv") or row.get("other_structural"))
        unclassified = max(0, total - hgt - te - other)
        vals = {
            "hgt_starship": hgt,
            "te_repeat": te,
            "other_structural": other,
            "unclassified": unclassified,
        }
        matrix.append([100 * vals[key] / total for key, _label in cols])
    im = ax.imshow(matrix, aspect="auto", cmap="YlGnBu", vmin=0, vmax=max(max(r) for r in matrix))
    ax.set_xticks([x - 0.5 for x in range(1, len(cols))], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, len(sample_order))], minor=True)
    ax.grid(which="minor", color="white", linewidth=1.3)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_yticks(range(len(sample_order)))
    ax.set_yticklabels([sample_label[s] for s in sample_order], fontstyle="italic")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([label for _key, label in cols])
    ax.set_title("Routed biological annotation is sample-specific", loc="left")
    cbar = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cbar.set_label("% of MycoSV loci", fontsize=8)
    cbar.ax.tick_params(labelsize=8)
    add_panel_label(ax, "e")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--by-query-dir", type=Path,
                    default=Path("experiments/million_real/full_fungal_assembly_20260522_070100/assembly/by_query"))
    ap.add_argument("--manuscript-dir", type=Path, default=Path("manuscript"))
    ap.add_argument("--out-prefix", type=Path, default=Path("manuscript/figures/fig1_assembly_five_sample"))
    args = ap.parse_args()

    table1 = read_tsv(args.manuscript_dir / "tables/table1_pangenome_vs_single_reference.tsv")
    table1_samples = [r for r in table1 if not r["sample"].startswith("Panel")]
    for r in table1_samples:
        r["sample"] = r["sample"].replace(".", "_")
    sample_order = [r["sample"] for r in table1_samples]
    sample_label = {r["sample"]: short_species(r["species"]) for r in table1_samples}
    rv_rows = read_tsv(args.manuscript_dir / "data/panel_read_validation_rate.tsv")
    lift_rows = read_tsv(args.manuscript_dir / "data/panel_pangenome_lift.tsv")

    calls: list[dict[str, object]] = []
    for sample in sample_order:
        calls.extend(read_vcf_calls(args.by_query_dir / sample / "mycosv/calls.vcf", sample))

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update({
        "font.family": ["DejaVu Sans", "Liberation Sans", "sans-serif"],
        "font.size": 8.5,
        "axes.titlesize": 10,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "savefig.dpi": 300,
    })
    fig = plt.figure(figsize=(13.8, 10.2), constrained_layout=False)
    gs = GridSpec(3, 4, figure=fig, height_ratios=[1.0, 1.05, 1.0], hspace=0.78, wspace=0.62)
    ax_a = fig.add_subplot(gs[0, :2])
    ax_b = fig.add_subplot(gs[0, 2:])
    ax_c = fig.add_subplot(gs[1, :2])
    ax_d = fig.add_subplot(gs[1, 2:])
    ax_e = fig.add_subplot(gs[2, 1:3])

    plot_validation(ax_a, rv_rows, sample_order, sample_label)
    plot_pangenome_lift(ax_b, lift_rows, sample_order, sample_label)
    plot_svlen_violin(ax_c, calls, sample_order, sample_label)
    plot_density(ax_d, calls, sample_order, sample_label)
    plot_biology_heatmap(ax_e, table1_samples, sample_order, sample_label)

    fig.suptitle("Assembly-mode pangenome SV discovery across five fungal genomes",
                 x=0.02, ha="left", y=0.985, fontsize=13, fontweight="bold")
    fig.text(0.02, 0.955,
             "Figure 1 uses assembly-mode calls only; long-read benchmarking is reported separately in Table 2.",
             ha="left", va="top", fontsize=8.5, color="#555555")
    fig.subplots_adjust(top=0.875, left=0.07, right=0.965, bottom=0.08)
    add_clean_legends(fig, ax_a, ax_b)
    for ext in ("png", "svg"):
        fig.savefig(args.out_prefix.with_suffix(f".{ext}"))
    plt.close(fig)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
