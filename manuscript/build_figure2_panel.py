#!/usr/bin/env python3
"""Build a Nature Biotechnology-style panel-200 Figure 2.

Figure 2 demonstrates the power of pangenome SV discovery across the 200-genome
fungal panel for biology, ecology, and scale. Panels answer:

  a  Q1  Cross-clade HGT / Starship-type cargo (hgt_flag is set when a novel
         segment matches a reference outside the query's own clade).
  b  Q2  Two-speed / TE-rich accessory architecture (TE + RIP + repeat element
         load per genome; RIP is the canonical genome-defense signature of
         fungal accessory compartments).
  c  Q3  Gene-proximal (genic) context of novel SVs, from the gene-model
         annotations. NOTE: direct RNA-seq expression support is NOT available
         in this run (no transcriptome was joined); this panel reports the
         structural nearby-gene basis for functional impact, which is the
         prerequisite layer for expression follow-up.
  d      Novel-SV biology partitions by fungal trophic mode (ecology).
  e      Biologically salient novel SVs are disproportionately pangenome-only
         (MycoSV-unique) -> the central "power of the pangenome" result.
  f      Discovery scales across the 200-genome panel (accumulation curve).

Style follows build_figure1_assembly.py (Okabe-Ito palette, panel labels,
left-aligned descriptive titles, png+svg @300 dpi) and the clean conceptual
clarity of classic fungal-genome-biology figures (two-speed genome, RIP, HGT,
Starships, trophic guilds).
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import random
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch

# ---- shared palette (refined Okabe-Ito, matching Figure 1 family) ---------
INK = "#22272B"           # near-black for text/axes
SCOPE_COLOR = {"single": "#DCDFE3", "pangenome": "#1F6F8B", "supported": "#2E8B6B"}
GUILD_COLOR = {
    "Saprotroph": "#1F6F8B", "Symbiotroph": "#2E8B6B",
    "Pathotroph": "#D1603D", "Mixed guild": "#8E6FAF", "Unassigned": "#C7CDD2",
}
ELEMENT_GROUP_COLOR = {
    "HGT / Starship": "#D1603D", "TE / RIP": "#E0A53A",
    "Repeat": "#B47AA8", "Other / unclassified": "#D9DDE1",
}
# sequential map for the ecology heatmap (white -> teal -> deep ink-teal)
HEAT_CMAP = LinearSegmentedColormap.from_list(
    "mycoteal", ["#F7FAFB", "#CDE3E6", "#7FB6BD", "#2F7E89", "#15454C"])
TE_RIP_CLASSES = {"RIP", "REPEAT", "TE_TIR", "TE_LINE", "TE_SINE", "TE_LTR"}
SVTYPES = {"INS", "DEL", "DUP", "INV", "TRA", "OFF_REF", "OTHER"}


def as_int(v, default=0):
    try:
        if v in (None, "", "."):
            return default
        return int(round(float(v)))
    except (ValueError, TypeError):
        return default


def truthy(v) -> bool:
    return str(v).strip().lower() in {"yes", "true", "1"}


def guild_of(trophic: str) -> str:
    t = (trophic or "").strip()
    if t in {"", "."}:
        return "Unassigned"
    if "-" in t:                       # e.g. Pathotroph-Saprotroph
        return "Mixed guild"
    base = t.split("_")[0]             # Saprotroph_wood_decay -> Saprotroph
    if base in {"Saprotroph", "Symbiotroph", "Pathotroph"}:
        return base
    return "Unassigned"


def element_group(ec: str) -> str:
    ec = (ec or "").strip().upper()
    if ec == "HGT":
        return "HGT / Starship"
    if ec in TE_RIP_CLASSES:
        return "TE / RIP" if ec in {"RIP"} or ec.startswith("TE_") else "Repeat"
    return "Other / unclassified"


# ---------------------------------------------------------------------------
def aggregate(by_query: Path):
    """Single streaming pass over all shards -> per-genome + global aggregates."""
    novel_files = sorted(glob.glob(str(by_query / "*" / "novel_mycosv_calls.tsv")))
    bio_files = sorted(glob.glob(str(by_query / "*" / "biology_findings.tsv")))

    per_genome = defaultdict(lambda: {
        "novel_total": 0, "pangenome_only": 0,
        "te_rip": 0, "hgt": 0, "genic": 0,
        "trophic": Counter(), "phylum": Counter(), "species": "",
    })
    # global tallies
    hgt_by_svtype = defaultdict(lambda: Counter())     # svtype -> {pangenome_only, single}
    discovery_bucket = Counter()
    panonly_within = {"all": [0, 0], "HGT": [0, 0], "TE_RIP": [0, 0], "genic": [0, 0]}  # [panonly, total]
    gene_dist = []
    biotype = Counter()
    guild_element = defaultdict(lambda: Counter())     # guild -> element_group counts

    # ---- novel_mycosv_calls.tsv: structural + pangenome-scope per genome ----
    for f in novel_files:
        g = Path(f).parent.name
        d = per_genome[g]
        try:
            rows = csv.DictReader(open(f, encoding="utf-8", errors="replace"), delimiter="\t")
        except OSError:
            continue
        for r in rows:
            d["novel_total"] += 1
            pan = truthy(r.get("mycosv_unique"))
            if pan:
                d["pangenome_only"] += 1
            ec = (r.get("element_class") or "").strip().upper()
            if ec in TE_RIP_CLASSES:
                d["te_rip"] += 1
            if ec == "HGT":
                d["hgt"] += 1
            discovery_bucket[(r.get("discovery_bucket") or ".").strip()] += 1
            panonly_within["all"][1] += 1
            panonly_within["all"][0] += int(pan)
            grp = element_group(ec)
            if grp == "HGT / Starship":
                panonly_within["HGT"][1] += 1
                panonly_within["HGT"][0] += int(pan)
            if grp == "TE / RIP" or grp == "Repeat":
                panonly_within["TE_RIP"][1] += 1
                panonly_within["TE_RIP"][0] += int(pan)

    # ---- biology_findings.tsv: HGT cross-clade, gene context, ecology ------
    for f in bio_files:
        g = Path(f).parent.name
        d = per_genome[g]
        try:
            rows = csv.DictReader(open(f, encoding="utf-8", errors="replace"), delimiter="\t")
        except OSError:
            continue
        for r in rows:
            tm = (r.get("trophic_mode") or ".").strip()
            if tm not in {"", "."}:
                d["trophic"][tm] += 1
            ph = (r.get("phylum") or ".").strip()
            if ph not in {"", "."}:
                d["phylum"][ph] += 1
            if not d["species"]:
                d["species"] = (r.get("species") or "").strip()
            # HGT cross-clade cargo by SV type
            if truthy(r.get("hgt_flag")):
                sv = (r.get("svtype") or "OTHER").strip().upper()
                sv = sv if sv in SVTYPES else "OTHER"
                pan = truthy(r.get("mycosv_unique"))
                hgt_by_svtype[sv]["pangenome_only" if pan else "single"] += 1
            # gene-proximal context
            gene = (r.get("affected_gene") or "").strip()
            dist = (r.get("affected_gene_distance_bp") or "").strip()
            if gene not in {"", ".", "NONE"} and dist not in {"", "."}:
                d["genic"] += 1
                gene_dist.append(as_int(dist))
                bt = (r.get("affected_gene_biotype") or ".").strip()
                biotype[bt if bt not in {"", "."} else "unspecified"] += 1
                panonly_within["genic"][1] += 1
                panonly_within["genic"][0] += int(truthy(r.get("mycosv_unique")))
            # ecology x element-group
            grp = element_group(r.get("element_class"))
            guild_element[guild_of(tm)][grp] += 1

    # ---- calls.vcf: cross-clade TRA junctions (cargo-free, HGT-associated) ---
    # TRA records carry no variant sequence, so they are never element-classified
    # as HGT cargo (panel a). But a TRA whose breakpoint anchors to a reference
    # at a SUPRA-GENUS clade rank is a cross-clade rearrangement junction -- the
    # adjacency signature of HGT / Starship integration. We surface it from the
    # VCF CLADE_RANK as a complementary, cargo-free track.
    import re as _re
    supra = {"family", "order", "class", "phylum"}
    tra_total = 0
    tra_crossclade = 0
    tra_rank = Counter()
    rank_re = _re.compile(r"CLADE_RANK=([^;\t]+)")
    for f in sorted(glob.glob(str(by_query / "*" / "mycosv" / "calls.vcf"))):
        try:
            fh = open(f, encoding="utf-8", errors="replace")
        except OSError:
            continue
        with fh:
            for line in fh:
                if line.startswith("#") or "SVTYPE=TRA" not in line:
                    continue
                tra_total += 1
                m = rank_re.search(line)
                rank = m.group(1) if m else "."
                tra_rank[rank] += 1
                if rank in supra:
                    tra_crossclade += 1

    # finalize per-genome labels
    for g, d in per_genome.items():
        d["guild"] = guild_of(d["trophic"].most_common(1)[0][0]) if d["trophic"] else "Unassigned"
        d["phylum_label"] = d["phylum"].most_common(1)[0][0] if d["phylum"] else "."

    return {
        "per_genome": dict(per_genome),
        "hgt_by_svtype": hgt_by_svtype,
        "discovery_bucket": discovery_bucket,
        "panonly_within": panonly_within,
        "gene_dist": gene_dist,
        "biotype": biotype,
        "guild_element": guild_element,
        "tra_total": tra_total,
        "tra_crossclade": tra_crossclade,
        "tra_rank": tra_rank,
        "n_genomes": len(novel_files),
    }


# ---- styling helpers -------------------------------------------------------
def panel_label(ax, lab):
    ax.text(-0.16, 1.05, lab, transform=ax.transAxes, fontsize=14,
            fontweight="bold", va="bottom", ha="left", color=INK)


def bold_legend(leg):
    """Bold every legend entry + its title (per user request)."""
    if leg is None:
        return leg
    for t in leg.get_texts():
        t.set_fontweight("bold")
    if leg.get_title() is not None:
        leg.get_title().set_fontweight("bold")
    return leg


def style_axes(ax):
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#9AA1A7")
        ax.spines[s].set_linewidth(0.9)
    ax.tick_params(colors=INK, length=3, width=0.8)
    ax.set_axisbelow(True)


def title_left(ax, text, pad=8):
    ax.set_title(text, loc="left", pad=pad, fontsize=9.6, color=INK, fontweight="bold")


# ---- panels ----------------------------------------------------------------
def panel_a_hgt(ax, agg):
    hbs = agg["hgt_by_svtype"]
    order = [s for s in ["INS", "DEL", "DUP", "INV", "OFF_REF", "OTHER"] if s in hbs]
    pan = [hbs[s]["pangenome_only"] for s in order]
    sgl = [hbs[s]["single"] for s in order]
    tra_x = agg.get("tra_crossclade", 0)         # cross-clade TRA junctions
    tra_all = agg.get("tra_total", 0)
    y = list(range(len(order)))
    ax.barh(y, pan, color=SCOPE_COLOR["pangenome"], height=0.64, zorder=3,
            label="pangenome-only cargo (MycoSV-unique)")
    ax.barh(y, sgl, left=pan, color=SCOPE_COLOR["single"], height=0.64, zorder=3,
            label="single-ref equivalent cargo")
    xmax = max(1, max([pan[i] + sgl[i] for i in range(len(order))] + [tra_x]))
    for i in range(len(order)):
        tot = pan[i] + sgl[i]
        if tot:
            ax.text(tot + xmax * 0.012, i, f"{tot:,}  ·  {100*pan[i]/tot:.0f}% pan",
                    va="center", ha="left", fontsize=7.3, color=INK, fontweight="bold")
    labels = list(order)
    # cross-clade TRA junction row (distinct: hatched grey, "junction" not cargo)
    ti = len(order)
    if tra_x:
        ax.axhline(ti - 0.5, color="#C7CDD2", lw=0.8, ls=(0, (2, 2)), zorder=1)
        ax.barh([ti], [tra_x], color="#E3E7EA", height=0.52, zorder=3,
                hatch="////", edgecolor="#7A848B", linewidth=0.5,
                label="cross-clade TRA junction (cargo-free)")
        ax.text(tra_x + xmax * 0.012, ti,
                f"{tra_x:,} junctions (class-rank anchored)",
                va="center", ha="left", fontsize=7.0, color="#5B656C", fontweight="bold")
        labels.append("TRA$^\\dagger$")
    ax.set_yticks(list(range(len(labels))))
    ax.set_yticklabels(labels, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel("Cross-clade HGT / Starship events (novel SV calls)")
    title_left(ax, "a   Cross-clade HGT / Starship cargo (+ junctions)")
    ax.grid(axis="x", color="#ECEFF1", lw=0.7)
    ax.set_xlim(0, xmax * 1.34)
    style_axes(ax)
    leg = ax.legend(loc="lower right", frameon=True, fontsize=7.1, framealpha=0.95,
                    edgecolor="#D5DADE", borderpad=0.6)
    leg.get_frame().set_linewidth(0.6)
    bold_legend(leg)
    ax.text(0.0, -0.205,
            "$\\dagger$ TRA = cross-clade rearrangement junction (no cargo sequence to class as HGT); "
            "shown as the adjacency signature of integration.",
            transform=ax.transAxes, fontsize=6.4, color="#7A848B", style="italic")


def panel_b_twospeed(ax, ax_top, ax_right, ax_leg, agg):
    xs, ys, cs, ss = [], [], [], []
    for d in agg["per_genome"].values():
        if d["novel_total"] < 30:
            continue
        xs.append(d["novel_total"])
        ys.append(100 * d["te_rip"] / max(1, d["novel_total"]))
        cs.append(GUILD_COLOR.get(d["guild"], GUILD_COLOR["Unassigned"]))
        ss.append(18 + 90 * (d["hgt"] / max(1, d["novel_total"])))
    ax.scatter(xs, ys, c=cs, s=ss, alpha=0.85, edgecolor="white", linewidth=0.5, zorder=3)
    med = sorted(ys)[len(ys) // 2] if ys else 0
    ax.axhline(med, color="#7A848B", lw=1.0, ls=(0, (4, 3)), zorder=1)
    ax.text(0.015, med, f" median {med:.0f}%", transform=ax.get_yaxis_transform(),
            va="bottom", ha="left", fontsize=7.0, color="#5B656C", fontweight="bold")
    ax.set_xscale("log")
    ax.set_xlabel("Novel SVs per genome (log scale)")
    ax.set_ylabel("TE + RIP + repeat share of novel SVs (%)")
    ax.grid(color="#ECEFF1", lw=0.7)
    style_axes(ax)
    # marginal distributions (joint-plot style)
    lo, hi = min(xs), max(xs)
    log_bins = [lo * (hi / lo) ** (i / 24) for i in range(25)]
    ax_top.hist(xs, bins=log_bins, color="#B9C2C8", edgecolor="white", linewidth=0.3)
    ax_top.set_xscale("log")
    ax_right.hist(ys, bins=18, orientation="horizontal", color="#B9C2C8",
                  edgecolor="white", linewidth=0.3)
    for a in (ax_top, ax_right):
        a.axis("off")
    ax_top.set_title("b   Two-speed signature: per-genome TE/RIP accessory load",
                     loc="left", pad=6, fontsize=9.6, color=INK, fontweight="bold")
    # guild legend in the empty top-right corner axis (never overlaps points)
    ax_leg.axis("off")
    handles = [mpl.lines.Line2D([0], [0], marker="o", linestyle="", color=c, label=g,
                                markersize=7, markeredgecolor="white")
               for g, c in GUILD_COLOR.items() if g != "Unassigned"]
    leg = ax_leg.legend(handles=handles, loc="center", frameon=False, fontsize=6.8,
                        title="trophic guild\n(size ∝ HGT fraction)", title_fontsize=6.8,
                        labelspacing=0.45, handletextpad=0.3, borderpad=0.1)
    bold_legend(leg)
    ax.text(-0.20, 1.30, "b", transform=ax.transAxes, fontsize=14,
            fontweight="bold", va="bottom", ha="left", color=INK)


def panel_c_gene(ax, agg):
    dists = [d for d in agg["gene_dist"] if d >= 0]
    logd = [math.log10(d + 1) for d in dists]
    ax.axvspan(0, 3, color="#E0A53A", alpha=0.10, zorder=0)   # <=1 kb band
    ax.hist(logd, bins=44, color=SCOPE_COLOR["pangenome"], alpha=0.92,
            edgecolor="white", linewidth=0.3, zorder=2)
    n1 = sum(1 for d in dists if d <= 1000)
    ax.axvline(3, color="#D1603D", lw=1.2, ls=(0, (4, 3)), zorder=3)
    ax.text(2.92, ax.get_ylim()[1] * 0.96,
            f"≤1 kb of a gene\n{n1:,} SVs ({100*n1/max(1,len(dists)):.0f}%)",
            va="top", ha="right", fontsize=7.4, color=INK, fontweight="bold")
    ax.set_xlabel("Distance from novel SV to nearest annotated gene (bp)")
    ax.set_ylabel("Novel SVs")
    ax.set_xticks([0, 1, 2, 3, 4, 5])
    ax.set_xticklabels(["0", "10", "$10^2$", "$10^3$", "$10^4$", "$10^5$"])
    title_left(ax, "c   Novel SVs are gene-proximal (genic-impact context)")
    ax.grid(axis="y", color="#ECEFF1", lw=0.7)
    style_axes(ax)
    # inset: nearest-gene biotype mix
    bt = agg["biotype"]
    top = [(k, v) for k, v in bt.most_common() if k not in {"unspecified"}][:4]
    if top:
        ins = ax.inset_axes([0.55, 0.46, 0.42, 0.46])
        labels = [k for k, _ in top][::-1]
        vals = [v for _, v in top][::-1]
        ins.barh(range(len(vals)), vals, color="#2F7E89", edgecolor="white", linewidth=0.4)
        ins.set_yticks(range(len(labels)))
        ins.set_yticklabels(labels, fontsize=6.2, fontweight="bold")
        ins.set_title("nearest-gene biotype", fontsize=6.6, fontweight="bold", color=INK, pad=2)
        ins.tick_params(labelsize=6.0, length=2)
        for sp in ("top", "right"):
            ins.spines[sp].set_visible(False)
        ins.set_xticks([])
        for j, v in enumerate(vals):
            ins.text(v, j, f" {v:,}", va="center", ha="left", fontsize=5.8, color=INK)
        ins.margins(x=0.18)
    panel_label(ax, "c")


def panel_d_ecology(ax, agg):
    guilds = [g for g in ["Saprotroph", "Symbiotroph", "Pathotroph", "Mixed guild"]
              if g in agg["guild_element"]]
    cols = ["HGT / Starship", "TE / RIP", "Repeat"]
    matrix = []
    for g in guilds:
        tot = sum(agg["guild_element"][g].values()) or 1
        matrix.append([100 * agg["guild_element"][g][c] / tot for c in cols])
    vmax = max(max(r) for r in matrix) if matrix else 1
    im = ax.imshow(matrix, aspect="auto", cmap=HEAT_CMAP, vmin=0, vmax=vmax)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(["HGT /\nStarship", "TE / RIP", "Repeat"], fontweight="bold")
    ax.set_yticks(range(len(guilds)))
    ax.set_yticklabels(guilds, fontweight="bold")
    ax.set_xticks([x - 0.5 for x in range(1, len(cols))], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, len(guilds))], minor=True)
    ax.grid(which="minor", color="white", lw=1.6)
    ax.tick_params(which="minor", length=0)
    ax.tick_params(which="major", length=0)
    for i in range(len(guilds)):
        for j in range(len(cols)):
            v = matrix[i][j]
            ax.text(j, i, f"{v:.1f}%", ha="center", va="center", fontsize=8.4,
                    fontweight="bold", color="white" if v > vmax * 0.55 else INK)
    title_left(ax, "d   Novel-SV biology partitions by fungal trophic mode")
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.set_label("% of novel-SV findings", fontsize=7.6, fontweight="bold")
    cb.ax.tick_params(labelsize=7)
    cb.outline.set_linewidth(0.5)
    panel_label(ax, "d")


def panel_e_pangenome(ax, agg):
    cats = [("all", "All\nnovel SVs"), ("HGT", "HGT /\nStarship"),
            ("TE_RIP", "TE / RIP"), ("genic", "Gene-\nproximal")]
    colors = ["#9AA7AE", ELEMENT_GROUP_COLOR["HGT / Starship"],
              ELEMENT_GROUP_COLOR["TE / RIP"], SCOPE_COLOR["pangenome"]]
    x = list(range(len(cats)))
    pct, ns = [], []
    for key, _ in cats:
        pan, tot = agg["panonly_within"][key]
        pct.append(100 * pan / max(1, tot)); ns.append(tot)
    ax.bar(x, pct, color=colors, width=0.62, edgecolor="white", linewidth=0.6, zorder=3)
    for i in range(len(cats)):
        ax.text(i, pct[i] + max(pct) * 0.02, f"{pct[i]:.0f}%", ha="center", va="bottom",
                fontsize=9.0, fontweight="bold", color=INK)
        ax.text(i, pct[i] + max(pct) * 0.02, f"\nn={ns[i]:,}", ha="center", va="top",
                fontsize=6.6, color="#5B656C")
    ax.axhline(pct[0], color="#9AA7AE", lw=1.0, ls=(0, (4, 3)), zorder=1)
    ax.text(len(cats) - 0.5, pct[0], " baseline", va="bottom", ha="right",
            fontsize=6.8, color="#5B656C", fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([lab for _, lab in cats], fontweight="bold")
    ax.set_ylabel("Pangenome-only (MycoSV-unique) share (%)")
    ax.set_ylim(0, max(pct) * 1.30)
    title_left(ax, "e   A pangenome-only SV layer (HGT, TE/RIP) single refs miss")
    ax.grid(axis="y", color="#ECEFF1", lw=0.7)
    style_axes(ax)
    panel_label(ax, "e")


def panel_f_scalability(ax, agg):
    per = agg["per_genome"]
    genomes = list(per.keys())
    n = len(genomes)
    rng = random.Random(1)
    acc_all, acc_pan = [], []
    for _ in range(60):
        order = genomes[:]
        rng.shuffle(order)
        ca = cp = 0
        a_run, p_run = [], []
        for gname in order:
            ca += per[gname]["novel_total"]; cp += per[gname]["pangenome_only"]
            a_run.append(ca); p_run.append(cp)
        acc_all.append(a_run); acc_pan.append(p_run)
    mean_all = [sum(c) / len(c) for c in zip(*acc_all)]
    mean_pan = [sum(c) / len(c) for c in zip(*acc_pan)]
    lo_all = [min(c) for c in zip(*acc_all)]; hi_all = [max(c) for c in zip(*acc_all)]
    xs = list(range(1, n + 1))
    ax.fill_between(xs, lo_all, hi_all, color=SCOPE_COLOR["pangenome"], alpha=0.12, lw=0)
    ax.fill_between(xs, 0, mean_pan, color=SCOPE_COLOR["supported"], alpha=0.16, lw=0)
    ax.plot(xs, mean_all, color=SCOPE_COLOR["pangenome"], lw=2.4, label="all novel SVs", zorder=4)
    ax.plot(xs, mean_pan, color=SCOPE_COLOR["supported"], lw=2.4,
            label="pangenome-only SVs", zorder=4)
    ax.scatter([n, n], [mean_all[-1], mean_pan[-1]], s=22, zorder=5,
               color=[SCOPE_COLOR["pangenome"], SCOPE_COLOR["supported"]],
               edgecolor="white", linewidth=0.6)
    ax.text(n, mean_all[-1], f"  {int(mean_all[-1]):,}", va="center", ha="left",
            fontsize=7.6, fontweight="bold", color=SCOPE_COLOR["pangenome"])
    ax.text(n, mean_pan[-1], f"  {int(mean_pan[-1]):,}", va="center", ha="left",
            fontsize=7.6, fontweight="bold", color=SCOPE_COLOR["supported"])
    ax.set_xlim(1, n * 1.12)
    ax.set_xlabel("Genomes sampled")
    ax.set_ylabel("Cumulative novel SV calls")
    title_left(ax, f"f   Discovery scales across the {n}-genome panel (no saturation)")
    ax.grid(color="#ECEFF1", lw=0.7)
    style_axes(ax)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    leg = ax.legend(loc="upper left", frameon=True, fontsize=7.8, framealpha=0.95,
                    edgecolor="#D5DADE", borderpad=0.6)
    leg.get_frame().set_linewidth(0.6)
    bold_legend(leg)
    panel_label(ax, "f")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--by-query-dir", type=Path,
                    default=Path("experiments/million_real/full_fungal_assembly_panel200_20260526_053633/assembly/by_query"))
    ap.add_argument("--out-prefix", type=Path,
                    default=Path("manuscript/figures/fig2_pangenome_panel200"))
    args = ap.parse_args()

    agg = aggregate(args.by_query_dir)
    print(f"[fig2] aggregated {agg['n_genomes']} genomes; "
          f"HGT svtypes={dict({k: sum(v.values()) for k, v in agg['hgt_by_svtype'].items()})}; "
          f"genic={len(agg['gene_dist'])}")

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update({
        "font.family": ["DejaVu Sans", "Liberation Sans", "sans-serif"],
        "font.size": 8.6, "axes.titlesize": 9.6, "axes.labelsize": 8.6,
        "axes.labelcolor": INK, "axes.edgecolor": "#9AA1A7",
        "text.color": INK, "xtick.color": INK, "ytick.color": INK,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.facecolor": "white", "savefig.dpi": 300, "savefig.bbox": "tight",
    })
    fig = plt.figure(figsize=(13.8, 14.6))
    gs = GridSpec(3, 2, figure=fig, hspace=0.52, wspace=0.26)
    panel_a_hgt(fig.add_subplot(gs[0, 0]), agg)
    # panel b is a joint plot: main + top/right marginals + legend corner
    gb = GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[0, 1],
                                 width_ratios=[4.0, 1.15], height_ratios=[1.15, 4.0],
                                 hspace=0.06, wspace=0.06)
    ax_b_main = fig.add_subplot(gb[1, 0])
    ax_b_top = fig.add_subplot(gb[0, 0], sharex=ax_b_main)
    ax_b_right = fig.add_subplot(gb[1, 1], sharey=ax_b_main)
    ax_b_leg = fig.add_subplot(gb[0, 1])
    panel_b_twospeed(ax_b_main, ax_b_top, ax_b_right, ax_b_leg, agg)
    panel_c_gene(fig.add_subplot(gs[1, 0]), agg)
    panel_d_ecology(fig.add_subplot(gs[1, 1]), agg)
    panel_e_pangenome(fig.add_subplot(gs[2, 0]), agg)
    panel_f_scalability(fig.add_subplot(gs[2, 1]), agg)

    fig.suptitle("Pangenome structural-variant discovery across 200 fungal genomes: biology, ecology and scale",
                 x=0.015, ha="left", y=0.992, fontsize=14, fontweight="bold", color=INK)
    fig.text(0.015, 0.965,
             "MycoSV-only novel calls (a, b, e) integrated with gene-model annotations and ecological traits (c, d). "
             "Panel c shows gene proximity — the structural basis for functional impact; direct RNA-seq expression was "
             "not joined in this run.",
             ha="left", va="top", fontsize=8.0, color="#5B656C")
    fig.subplots_adjust(top=0.905, left=0.07, right=0.965, bottom=0.05)
    for ext in ("png", "svg"):
        fig.savefig(args.out_prefix.with_suffix(f".{ext}"))
    plt.close(fig)
    print(f"[fig2] wrote {args.out_prefix}.png / .svg")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
