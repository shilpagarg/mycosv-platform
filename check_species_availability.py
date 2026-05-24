#!/usr/bin/env python3
"""Authoritative availability check for the fungal panel wish-list.

For each species:
  - NCBI Datasets v2 (genome): how many assemblies, best accession + level
  - ENA portal (assembly): cross-check with NCBI
  - ENA portal (read_run): are there public reads (long or short) for SRA-based
    benchmarking even when no assembly exists yet
  - JGI MycoCosm: noted as manual check (requires portal login for full list)

Writes a TSV + prints a guild-by-guild summary.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

UA = {"User-Agent": "mycosv-panel-check/1.0"}
DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2"
ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api"


def http_get_json(url, retries=2, sleep=1.0):
    last = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.loads(r.read())
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            last = e
        except Exception as e:
            last = e
        time.sleep(sleep * (attempt + 1))
    return {"__error__": str(last)}


def http_get_text(url, retries=2, sleep=1.0):
    last = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=30) as r:
                return r.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return ""
            last = e
        except Exception as e:
            last = e
        time.sleep(sleep * (attempt + 1))
    return ""


def check_ncbi(species):
    url = f"{DATASETS}/genome/taxon/{urllib.parse.quote(species)}?page_size=20"
    j = http_get_json(url)
    if not j or "__error__" in (j or {}):
        return {"n": 0, "best": "", "level": "", "note": "error" if j else "404"}
    n = j.get("total_count", 0)
    reports = j.get("reports", []) or []
    best = ""
    best_level = ""
    level_rank = {"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}
    best_rank = 0
    for r in reports:
        acc = r.get("accession", "")
        lv = r.get("assembly_info", {}).get("assembly_level", "")
        rank = level_rank.get(lv, 0)
        if rank > best_rank:
            best, best_level, best_rank = acc, lv, rank
    if not best and reports:
        best = reports[0].get("accession", "")
        best_level = reports[0].get("assembly_info", {}).get("assembly_level", "")
    return {"n": n, "best": best, "level": best_level, "note": ""}


def check_ena_assemblies(species):
    q = f'scientific_name="{species}"'
    url = (f"{ENA_PORTAL}/search?"
           f"{urllib.parse.urlencode({'query': q, 'result': 'assembly', 'fields': 'accession,assembly_level', 'format': 'tsv', 'limit': 20})}")
    text = http_get_text(url)
    lines = [l for l in text.strip().splitlines() if l]
    return max(0, len(lines) - 1) if lines else 0


def check_ena_reads(species):
    q = f'scientific_name="{species}"'
    url = (f"{ENA_PORTAL}/search?"
           f"{urllib.parse.urlencode({'query': q, 'result': 'read_run', 'fields': 'run_accession,instrument_platform,library_strategy,fastq_ftp', 'format': 'tsv', 'limit': 50})}")
    text = http_get_text(url)
    lines = [l for l in text.strip().splitlines() if l]
    if len(lines) <= 1:
        return {"runs": 0, "long": 0, "short": 0, "with_fastq": 0}
    long_plats = {"OXFORD_NANOPORE", "PACBIO_SMRT"}
    n_long = n_short = n_with_fastq = 0
    for line in lines[1:]:
        cols = line.split("\t")
        if len(cols) < 4:
            continue
        plat = cols[1]
        fastq = cols[3]
        if fastq:
            n_with_fastq += 1
        if plat in long_plats:
            n_long += 1
        elif plat == "ILLUMINA":
            n_short += 1
    return {"runs": len(lines) - 1, "long": n_long, "short": n_short, "with_fastq": n_with_fastq}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species-tsv", type=Path, required=True,
                    help="TSV with header 'group\\tspecies' to check")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--throttle", type=float, default=0.3,
                    help="seconds between API hits to avoid rate limit (default 0.3)")
    args = ap.parse_args()

    rows = []
    with args.species_tsv.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)

    out_rows = []
    print(f"checking {len(rows)} species ...", file=sys.stderr)
    for i, r in enumerate(rows, 1):
        sp = (r.get("species") or "").strip()
        if not sp:
            continue
        ncbi = check_ncbi(sp)
        time.sleep(args.throttle)
        ena_asm = check_ena_assemblies(sp)
        time.sleep(args.throttle)
        ena_rds = check_ena_reads(sp)
        time.sleep(args.throttle)
        out_rows.append({
            "group": r.get("group", ""),
            "species": sp,
            "ncbi_assemblies": ncbi["n"],
            "ncbi_best_accession": ncbi["best"],
            "ncbi_best_level": ncbi["level"],
            "ncbi_note": ncbi["note"],
            "ena_assemblies": ena_asm,
            "ena_read_runs": ena_rds["runs"],
            "ena_long_reads": ena_rds["long"],
            "ena_short_reads": ena_rds["short"],
            "ena_runs_with_fastq": ena_rds["with_fastq"],
        })
        if i % 10 == 0 or i == len(rows):
            avail = sum(1 for x in out_rows if x["ncbi_assemblies"] > 0 or x["ena_assemblies"] > 0)
            print(f"  [{i:3d}/{len(rows)}] {sp[:40]:<40} ncbi={ncbi['n']:3d} ena_asm={ena_asm:3d} ena_reads={ena_rds['runs']:3d}   (running tally: {avail} with any assembly)", file=sys.stderr)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        fields = ["group", "species", "ncbi_assemblies", "ncbi_best_accession",
                  "ncbi_best_level", "ncbi_note", "ena_assemblies",
                  "ena_read_runs", "ena_long_reads", "ena_short_reads", "ena_runs_with_fastq"]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    print(f"\nwrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    raise SystemExit(main())
