#!/usr/bin/env python3
"""Intra-species population SV analysis for the A. flavus strain panel.

Every strain was called against a common held-out reference (NRRL 3357,
GCA_055697395.1), so calls share a coordinate system. This script merges SV
loci across strains into a strain-by-locus occupancy matrix and reports the
site-frequency spectrum: how many SV loci are private to one strain, segregate
at intermediate frequency, or are near-fixed across the panel.

Occupancy is presence of a call, not a genotype: MycoSV assigns per-strain
support but does not jointly genotype the population, so frequencies are
call-frequency estimates and are not corrected for per-strain detection power.

Usage: build_population_sfs.py <by_query_dir> <out_prefix>
"""
import csv
import glob
import os
import sys
from collections import defaultdict

# Per-type merge tolerances (bp), matching the benchmark's truth-anchored gates.
TOL = {"INS": 500, "DEL": 2500, "DUP": 2500, "INV": 10000, "TRA": 10000,
       "OFF_REF": 500}
LEN_FRAC = 0.30  # fractional length tolerance for merging same-type events


def canon(t):
    t = (t or "").upper()
    if t in ("INS", "DUP"):
        return "INSDUP"
    if t in ("DEL",):
        return "DEL"
    if t in ("INV",):
        return "INV"
    if t in ("TRA", "BND"):
        return "TRA"
    if "OFF" in t:
        return "OFF_REF"
    return t


def parse_info(info):
    d = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


def load_calls(vcf):
    """Yield (ref_contig, pos, end, canon_type, svlen) for one strain's VCF."""
    out = []
    if not os.path.exists(vcf):
        return out
    for line in open(vcf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 8:
            continue
        info = parse_info(f[7])
        # Prefer reference coordinates when present so loci are comparable
        # across strains; fall back to the VCF CHROM/POS otherwise.
        contig = info.get("REFCONTIG") or f[0]
        try:
            pos = int(info.get("REFPOS") or f[1])
            end = int(info.get("REFEND") or info.get("END") or pos)
        except ValueError:
            continue
        if end < pos:
            pos, end = end, pos
        try:
            svlen = abs(int(info.get("SVLEN") or (end - pos)))
        except ValueError:
            svlen = end - pos
        out.append((contig, pos, end, canon(info.get("SVTYPE")), svlen))
    return out


def main():
    by_query = sys.argv[1]
    out_prefix = sys.argv[2]

    strains = sorted(
        d for d in os.listdir(by_query)
        if os.path.isdir(os.path.join(by_query, d))
    )
    # Only strains whose call set is a completed (non-censored) run.
    usable = []
    censored = []
    per_strain = {}
    for s in strains:
        vcf = os.path.join(by_query, s, "mycosv", "calls.vcf")
        if not os.path.exists(vcf):
            continue
        if os.path.exists(os.path.join(by_query, s, "MYCOSV_FAILED.txt")):
            censored.append(s)
        calls = load_calls(vcf)
        per_strain[s] = calls
        usable.append(s)

    # Merge loci across strains, per (contig, canon_type), greedily by position.
    # A locus accumulates the set of strains whose calls fall within tolerance.
    loci = []  # each: dict(contig,type,pos,end,strains=set)
    index = defaultdict(list)  # (contig,type) -> list of locus dicts
    for s in usable:
        for contig, pos, end, ctype, svlen in per_strain[s]:
            tol = TOL.get(ctype, TOL.get("INS", 500))
            hit = None
            for loc in index[(contig, ctype)]:
                if abs(loc["pos"] - pos) <= tol and abs(loc["end"] - end) <= tol:
                    lo = min(loc["svlen"], svlen)
                    hi = max(loc["svlen"], svlen)
                    if hi == 0 or (hi - lo) <= LEN_FRAC * hi + tol:
                        hit = loc
                        break
            if hit is None:
                hit = {"contig": contig, "type": ctype, "pos": pos,
                       "end": end, "svlen": svlen, "strains": set()}
                loci.append(hit)
                index[(contig, ctype)].append(hit)
            hit["strains"].add(s)

    n = len(usable)
    # Site-frequency spectrum over occupancy count.
    sfs = defaultdict(int)
    for loc in loci:
        sfs[len(loc["strains"])] += 1
    singleton = sfs.get(1, 0)
    near_fixed = sum(v for k, v in sfs.items() if k >= max(1, int(0.9 * n)))
    intermediate = sum(v for k, v in sfs.items() if 1 < k < max(1, int(0.9 * n)))
    total = len(loci)

    # Write occupancy matrix (locus x strain) and the SFS.
    with open(f"{out_prefix}.occupancy.tsv", "w") as fh:
        fh.write("contig\ttype\tpos\tend\tn_strains\tfreq\t" + "\t".join(usable) + "\n")
        for loc in sorted(loci, key=lambda x: (-len(x["strains"]), x["contig"], x["pos"])):
            row = [loc["contig"], loc["type"], str(loc["pos"]), str(loc["end"]),
                   str(len(loc["strains"])), f"{len(loc['strains'])/n:.4f}"]
            row += ["1" if s in loc["strains"] else "0" for s in usable]
            fh.write("\t".join(row) + "\n")

    with open(f"{out_prefix}.sfs.tsv", "w") as fh:
        fh.write("occupancy_count\tallele_frequency\tn_loci\n")
        for k in sorted(sfs):
            fh.write(f"{k}\t{k/n:.4f}\t{sfs[k]}\n")

    with open(f"{out_prefix}.summary.tsv", "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"n_strains\t{n}\n")
        fh.write(f"n_strains_censored_included\t{len(censored)}\n")
        fh.write(f"total_sv_loci\t{total}\n")
        fh.write(f"singleton_loci\t{singleton}\n")
        fh.write(f"singleton_fraction\t{singleton/total:.4f}\n" if total else "singleton_fraction\tNA\n")
        fh.write(f"intermediate_freq_loci\t{intermediate}\n")
        fh.write(f"near_fixed_loci\t{near_fixed}\n")
        med = sorted(len(per_strain[s]) for s in usable)
        fh.write(f"median_calls_per_strain\t{med[len(med)//2] if med else 0}\n")

    print(f"strains={n} (censored included={len(censored)})  total_loci={total}")
    print(f"singleton={singleton} ({singleton/total:.1%})  "
          f"intermediate={intermediate}  near_fixed={near_fixed}"
          if total else "no loci")
    print(f"wrote {out_prefix}.{{occupancy,sfs,summary}}.tsv")


if __name__ == "__main__":
    main()
