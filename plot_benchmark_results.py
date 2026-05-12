#!/usr/bin/env python3
"""
Visualize precision/recall benchmark results for the fungal SV calling pipeline.

Reads pr_metrics.tsv and pr_metrics_by_scenario.tsv from each mode directory
under a simulated experiment, and optionally from real data benchmark output.

Usage
-----
  python plot_benchmark_results.py --exp-dir experiments/simulated/20260424_194837/benchmarks
  python plot_benchmark_results.py            # auto-discovers the latest experiment

Outputs
-------
  <exp-dir>/plots/fig1_f1_by_svtype.png
  <exp-dir>/plots/fig2_pr_scatter.png
  <exp-dir>/plots/fig3_scenario_heatmap_<mode>.png
  <exp-dir>/plots/fig4_mode_comparison.png
  <exp-dir>/plots/fig5_biology_summary.png   (if biology_candidates.tsv exists)
  <exp-dir>/plots/summary_table.tsv
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

SV_TYPE_ORDER = ["DEL", "INS", "DUP", "INV", "TRA", "OFF_REF"]
SV_COLORS = {
    "DEL": "#e07b39",
    "INS": "#4e9abe",
    "DUP": "#6ab187",
    "INV": "#c0495a",
    "TRA": "#9b59b6",
    "OFF_REF": "#7f8c8d",
    "OVERALL": "#2c3e50",
}

MODE_LABELS = {
    "assembly":         "Assembly",
    "short-reads":      "Short-reads",
    "long-reads_hifi":  "Long-reads (HiFi)",
    "long-reads_ont-r10": "Long-reads (ONT R10)",
    "long-reads_ont-r9":  "Long-reads (ONT R9)",
    "long-reads":       "Long-reads",
}
MODE_COLORS = ["#2c3e50", "#2980b9", "#27ae60", "#8e44ad", "#e67e22"]


def _label(mode: str) -> str:
    return MODE_LABELS.get(mode, mode)


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def float_or(val: str | None, default: float = 0.0) -> float:
    try:
        out = float(val)
        if np.isnan(out):
            return default
        return out
    except (TypeError, ValueError):
        return default


def metric_or_nan(val: str | None) -> float:
    return float_or(val, np.nan)


def int_or(val: str | None, default: int = 0) -> int:
    try:
        return int(float(val))
    except (TypeError, ValueError):
        return default


def nonnegative_err(center: float, lo: float, hi: float) -> tuple[float, float]:
    if any(np.isnan(v) for v in (center, lo, hi)):
        return 0.0, 0.0
    return max(0.0, center - lo), max(0.0, hi - center)


# ─────────────────────────────────────────────────────────────────────────────
# Discovery
# ─────────────────────────────────────────────────────────────────────────────

def discover_exp_dir(root: Path) -> Path:
    sim_root = root / "experiments" / "simulated"
    if not sim_root.exists():
        raise FileNotFoundError(f"No experiments/simulated/ under {root}")
    candidates = sorted(
        [d for d in sim_root.iterdir() if d.is_dir() and (d / "benchmarks").exists()],
        reverse=True,
    )
    if not candidates:
        raise FileNotFoundError(f"No timestamped experiment dirs found in {sim_root}")
    bench = candidates[0] / "benchmarks"
    print(f"[auto-discover] using {bench}", flush=True)
    return bench


def discover_modes(exp_dir: Path) -> list[str]:
    modes: list[str] = []
    for d in sorted(exp_dir.iterdir()):
        if d.is_dir() and (d / "pr_metrics.tsv").exists():
            modes.append(d.name)
    return modes


# ─────────────────────────────────────────────────────────────────────────────
# Data loading
# ─────────────────────────────────────────────────────────────────────────────

def load_mode_metrics(exp_dir: Path, modes: list[str]) -> dict[str, list[dict]]:
    """Returns {mode: [row_dict, ...]} from pr_metrics.tsv per mode."""
    data: dict[str, list[dict]] = {}
    for mode in modes:
        rows = read_tsv(exp_dir / mode / "pr_metrics.tsv")
        if rows:
            data[mode] = rows
    return data


def load_scenario_metrics(exp_dir: Path, modes: list[str]) -> dict[str, list[dict]]:
    """Returns {mode: [row_dict, ...]} from pr_metrics_by_scenario.tsv per mode."""
    data: dict[str, list[dict]] = {}
    for mode in modes:
        rows = read_tsv(exp_dir / mode / "pr_metrics_by_scenario.tsv")
        if rows:
            data[mode] = rows
    return data


def load_biology(exp_dir: Path, modes: list[str]) -> dict[str, list[dict]]:
    data: dict[str, list[dict]] = {}
    for mode in modes:
        bio_tsv = exp_dir / mode / "biology" / "biology_candidates.tsv"
        rows = read_tsv(bio_tsv)
        if rows:
            data[mode] = rows
    return data


# ─────────────────────────────────────────────────────────────────────────────
# Figure 1: F1 by SV type, grouped by mode
# ─────────────────────────────────────────────────────────────────────────────

def fig_f1_by_svtype(mode_metrics: dict[str, list[dict]], out: Path) -> None:
    modes = list(mode_metrics.keys())
    svtypes = [t for t in SV_TYPE_ORDER if
               any(any(r["svtype"] == t for r in mode_metrics[m]) for m in modes)]

    x = np.arange(len(svtypes))
    width = 0.8 / max(len(modes), 1)

    fig, ax = plt.subplots(figsize=(max(8, 2 * len(svtypes)), 5))
    for i, mode in enumerate(modes):
        by_type = {r["svtype"]: r for r in mode_metrics[mode]}
        f1s = [metric_or_nan(by_type.get(t, {}).get("f1")) for t in svtypes]
        tps = [int_or(by_type.get(t, {}).get("tp")) for t in svtypes]
        offset = (i - len(modes) / 2 + 0.5) * width
        heights = [0.0 if np.isnan(v) else v for v in f1s]
        bars = ax.bar(x + offset, heights, width * 0.9,
                      label=_label(mode), color=MODE_COLORS[i % len(MODE_COLORS)],
                      alpha=0.85)
        for bar, f1, tp in zip(bars, f1s, tps):
            if not np.isnan(f1) and f1 > 0.02:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                        f"{f1:.2f}", ha="center", va="bottom", fontsize=6.5)

    ax.set_xticks(x)
    ax.set_xticklabels(svtypes, fontsize=10)
    ax.set_ylim(0, 1.12)
    ax.set_ylabel("F1 score", fontsize=11)
    ax.set_title("F1 score by SV type and query mode", fontsize=13)
    ax.legend(fontsize=9, loc="upper right")
    ax.axhline(1.0, color="gray", linewidth=0.5, linestyle="--")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  saved {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 2: Precision vs Recall scatter (with CI error bars)
# ─────────────────────────────────────────────────────────────────────────────

def fig_pr_scatter(mode_metrics: dict[str, list[dict]], out: Path) -> None:
    modes = list(mode_metrics.keys())
    svtypes = [t for t in SV_TYPE_ORDER if
               any(any(r["svtype"] == t for r in mode_metrics[m]) for m in modes)]

    fig, ax = plt.subplots(figsize=(7, 7))
    for i, mode in enumerate(modes):
        by_type = {r["svtype"]: r for r in mode_metrics[mode]}
        color = MODE_COLORS[i % len(MODE_COLORS)]
        for j, t in enumerate(svtypes):
            row = by_type.get(t, {})
            prec = metric_or_nan(row.get("precision"))
            rec = metric_or_nan(row.get("recall"))
            if np.isnan(prec) or np.isnan(rec) or (prec == 0 and rec == 0):
                continue
            p_lo = metric_or_nan(row.get("prec_lo95"))
            p_hi = metric_or_nan(row.get("prec_hi95"))
            r_lo = metric_or_nan(row.get("rec_lo95"))
            r_hi = metric_or_nan(row.get("rec_hi95"))
            r_err_lo, r_err_hi = nonnegative_err(rec, r_lo, r_hi)
            p_err_lo, p_err_hi = nonnegative_err(prec, p_lo, p_hi)
            xerr = [[r_err_lo], [r_err_hi]]
            yerr = [[p_err_lo], [p_err_hi]]
            marker = ["o", "s", "^", "D", "v", "P"][j % 6]
            ax.errorbar(rec, prec, xerr=xerr, yerr=yerr,
                        fmt=marker, color=color, markersize=8,
                        capsize=3, alpha=0.85,
                        label=f"{_label(mode)} / {t}" if i == 0 else None)
            ax.annotate(t, (rec, prec), textcoords="offset points",
                        xytext=(4, 4), fontsize=7, color=color)

    # Legend for modes
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], marker="o", color="w",
                               markerfacecolor=MODE_COLORS[i % len(MODE_COLORS)],
                               markersize=10, label=_label(m))
                       for i, m in enumerate(modes)]
    ax.legend(handles=legend_elements, fontsize=9, loc="lower left")
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("Recall", fontsize=11)
    ax.set_ylabel("Precision", fontsize=11)
    ax.set_title("Precision vs Recall by SV type and mode", fontsize=13)
    ax.axhline(1.0, color="gray", linewidth=0.5, linestyle="--")
    ax.axvline(1.0, color="gray", linewidth=0.5, linestyle="--")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  saved {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 3: Per-scenario F1 heatmap (one per mode)
# ─────────────────────────────────────────────────────────────────────────────

def fig_scenario_heatmap(scenario_metrics: dict[str, list[dict]],
                         plots_dir: Path) -> None:
    for mode, rows in scenario_metrics.items():
        scenarios = sorted({r["scenario"] for r in rows})
        svtypes = [t for t in SV_TYPE_ORDER
                   if any(r["svtype"] == t for r in rows)]
        if not scenarios or not svtypes:
            continue

        matrix = np.full((len(scenarios), len(svtypes)), np.nan)
        for r in rows:
            if r["svtype"] not in svtypes:
                continue
            si = scenarios.index(r["scenario"])
            ti = svtypes.index(r["svtype"])
            matrix[si, ti] = metric_or_nan(r.get("f1"))

        fig, ax = plt.subplots(figsize=(max(5, len(svtypes) * 1.2),
                                         max(3, len(scenarios) * 1.0)))
        masked = np.ma.masked_invalid(matrix)
        cmap = cm.RdYlGn.copy()
        cmap.set_bad("lightgrey")
        im = ax.imshow(masked, cmap=cmap, vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(len(svtypes)))
        ax.set_xticklabels(svtypes, fontsize=10)
        ax.set_yticks(range(len(scenarios)))
        ax.set_yticklabels(scenarios, fontsize=9)
        for si in range(len(scenarios)):
            for ti in range(len(svtypes)):
                val = matrix[si, ti]
                if not np.isnan(val):
                    ax.text(ti, si, f"{val:.2f}", ha="center", va="center",
                            fontsize=8, color="black" if val > 0.5 else "white")
        plt.colorbar(im, ax=ax, label="F1 score")
        ax.set_title(f"F1 by scenario × SV type  [{_label(mode)}]", fontsize=12)
        fig.tight_layout()
        out = plots_dir / f"fig3_scenario_heatmap_{mode.replace('/', '_')}.png"
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(f"  saved {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 4: Mode comparison (overall precision / recall / F1)
# ─────────────────────────────────────────────────────────────────────────────

def fig_mode_comparison(mode_metrics: dict[str, list[dict]], out: Path) -> None:
    modes = list(mode_metrics.keys())
    metrics = ["precision", "recall", "f1"]
    metric_labels = ["Precision", "Recall", "F1"]

    x = np.arange(len(modes))
    width = 0.25
    colors = ["#2980b9", "#27ae60", "#e74c3c"]

    fig, ax = plt.subplots(figsize=(max(6, len(modes) * 1.8), 5))
    for i, (metric, mlabel, color) in enumerate(zip(metrics, metric_labels, colors)):
        vals = []
        for mode in modes:
            overall = next((r for r in mode_metrics[mode] if r["svtype"] == "OVERALL"), {})
            vals.append(metric_or_nan(overall.get(metric)))
        heights = [0.0 if np.isnan(v) else v for v in vals]
        bars = ax.bar(x + (i - 1) * width, heights, width * 0.9,
                      label=mlabel, color=color, alpha=0.85)
        for bar, val in zip(bars, vals):
            if np.isnan(val):
                continue
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
                    f"{val:.3f}", ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels([_label(m) for m in modes], fontsize=9, rotation=15, ha="right")
    ax.set_ylim(0, 1.12)
    ax.set_ylabel("Score", fontsize=11)
    ax.set_title("Overall precision / recall / F1 by query mode", fontsize=13)
    ax.legend(fontsize=10)
    ax.axhline(1.0, color="gray", linewidth=0.5, linestyle="--")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  saved {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 5: Biology candidates summary
# ─────────────────────────────────────────────────────────────────────────────

def fig_biology(biology_data: dict[str, list[dict]], out: Path) -> None:
    if not biology_data:
        return
    # Aggregate across modes
    type_counts: dict[str, dict[str, int]] = {}
    for mode, rows in biology_data.items():
        cnt: dict[str, int] = {}
        for row in rows:
            ct = row.get("candidate_type", "other")
            cnt[ct] = cnt.get(ct, 0) + 1
        type_counts[_label(mode)] = cnt

    all_types = sorted({t for cnts in type_counts.values() for t in cnts})
    if not all_types:
        return

    modes_list = list(type_counts.keys())
    x = np.arange(len(all_types))
    width = 0.8 / max(len(modes_list), 1)

    fig, ax = plt.subplots(figsize=(max(8, len(all_types) * 1.5), 5))
    for i, mode in enumerate(modes_list):
        vals = [type_counts[mode].get(t, 0) for t in all_types]
        offset = (i - len(modes_list) / 2 + 0.5) * width
        ax.bar(x + offset, vals, width * 0.9,
               label=mode, color=MODE_COLORS[i % len(MODE_COLORS)], alpha=0.85)

    ax.set_xticks(x)
    ax.set_xticklabels(all_types, fontsize=8, rotation=20, ha="right")
    ax.set_ylabel("Count", fontsize=11)
    ax.set_title("Biology candidate types by mode", fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  saved {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Summary TSV
# ─────────────────────────────────────────────────────────────────────────────

def write_summary_table(mode_metrics: dict[str, list[dict]], out: Path) -> None:
    rows: list[dict] = []
    for mode, metrics in mode_metrics.items():
        by_type = {r["svtype"]: r for r in metrics}
        for svtype in ["OVERALL"] + SV_TYPE_ORDER:
            r = by_type.get(svtype, {})
            if not r:
                continue
            rows.append({
                "mode": mode,
                "svtype": svtype,
                "tp": r.get("tp", ""),
                "fp": r.get("fp", ""),
                "fn": r.get("fn", ""),
                "precision": r.get("precision", ""),
                "recall": r.get("recall", ""),
                "f1": r.get("f1", ""),
            })
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["mode", "svtype", "tp", "fp", "fn", "precision", "recall", "f1"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"  saved {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exp-dir", type=Path, default=None,
                    help="Path to the benchmarks/ directory of a simulated experiment. "
                         "Auto-discovered from the project root if omitted.")
    ap.add_argument("--root", type=Path,
                    default=Path(__file__).resolve().parent,
                    help="Project root (default: directory of this script).")
    ap.add_argument("--modes", default="",
                    help="Comma-separated list of mode sub-directories to include. "
                         "Default: all directories that contain pr_metrics.tsv.")
    ap.add_argument("--plots-dir", type=Path, default=None,
                    help="Output directory for plots (default: <exp-dir>/plots/).")
    ap.add_argument("--real-data-dirs", nargs="*", type=Path, default=[],
                    help="One or more real data benchmark output directories containing calls.vcf.")
    ap.add_argument("--query-manifests", nargs="*", type=Path, default=[],
                    help="query_manifest.tsv files from real data prepare runs (for biological analysis).")
    ap.add_argument("--run-bio-analysis", action="store_true",
                    help="Run analyze_phylo_sv_biology.py on collected real data SV calls.")
    args = ap.parse_args()

    exp_dir = args.exp_dir
    if exp_dir is None:
        exp_dir = discover_exp_dir(args.root)
    exp_dir = exp_dir.resolve()
    if not exp_dir.exists():
        print(f"ERROR: exp-dir not found: {exp_dir}", file=sys.stderr)
        return 1

    plots_dir = args.plots_dir or exp_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    modes: list[str] = []
    if args.modes:
        modes = [m.strip() for m in args.modes.split(",") if m.strip()]
    else:
        modes = discover_modes(exp_dir)
    if not modes:
        print(f"ERROR: no mode directories with pr_metrics.tsv found in {exp_dir}",
              file=sys.stderr)
        return 1

    print(f"Modes found: {modes}")

    mode_metrics = load_mode_metrics(exp_dir, modes)
    scenario_metrics = load_scenario_metrics(exp_dir, modes)
    biology_data = load_biology(exp_dir, modes)

    if not mode_metrics:
        print("ERROR: no pr_metrics.tsv data loaded.", file=sys.stderr)
        return 1

    print("Generating plots…")
    fig_f1_by_svtype(mode_metrics, plots_dir / "fig1_f1_by_svtype.png")
    fig_pr_scatter(mode_metrics, plots_dir / "fig2_pr_scatter.png")
    fig_scenario_heatmap(scenario_metrics, plots_dir)
    fig_mode_comparison(mode_metrics, plots_dir / "fig4_mode_comparison.png")
    if biology_data:
        fig_biology(biology_data, plots_dir / "fig5_biology_summary.png")
    write_summary_table(mode_metrics, plots_dir / "summary_table.tsv")

    # Run biological analysis if real data dirs provided
    if args.run_bio_analysis and (args.real_data_dirs or args.query_manifests):
        bio_script = Path(__file__).resolve().parent / "analyze_phylo_sv_biology.py"
        if bio_script.exists():
            bio_out = plots_dir.parent / "bio_analysis"
            cmd = [sys.executable, str(bio_script), "--out-dir", str(bio_out)]
            for d in (args.real_data_dirs or []):
                cmd += ["--vcf-dirs", str(d)]
                bio_dir = d / "biology"
                if bio_dir.exists():
                    cmd += ["--biology-dirs", str(bio_dir)]
            for m in (args.query_manifests or []):
                cmd += ["--query-manifests", str(m)]
            print(f"\nRunning biological analysis → {bio_out}/")
            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError as e:
                print(f"[warn] biological analysis failed: {e}", file=sys.stderr)

    print(f"\nAll plots written to {plots_dir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
