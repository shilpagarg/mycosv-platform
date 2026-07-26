#!/usr/bin/env python3
"""Re-annotate the 200-genome panel element_class with the fixed classifier.

For each held-out query: run the compiled reannotate tool (which applies the
fixed classify_repeat_element: structural-first, Pezizomycotina Starship gate,
RIP filamentous-only gate, adaptive HGT threshold) over the query's VCF VARSEQ,
passing the query's taxonomic class (from the manifest) and background GC
(computed once per benchmark reference). Aggregate element_class counts across
all 200 into the manuscript categories.
"""
import csv
import gzip
import os
import subprocess
import sys
from collections import defaultdict

ROOT = os.environ.get("MYCOSV_ROOT", os.getcwd())
PANEL = f"{ROOT}/experiments/million_real/full_fungal_assembly_panel200_20260526_053633"
REANN = "/tmp/reannotate"
OUT = f"{ROOT}/experiments/biology_benchmark/panel200_reann"
os.makedirs(OUT, exist_ok=True)


def gc_of(path):
    if not path or not os.path.exists(path):
        return 0.45
    op = gzip.open if path.endswith(".gz") else open
    g = a = 0
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            s = line.upper()
            g += s.count("G") + s.count("C")
            a += s.count("A") + s.count("T")
    return g / (g + a) if (g + a) else 0.45


rows = list(csv.DictReader(open(f"{PANEL}/prepared/query_manifest.tsv"), delimiter="\t"))
gc_cache = {}
agg = defaultdict(int)
per_genome = []
done = 0
for r in rows:
    q = r["query_asm"]
    vcf = f"{PANEL}/assembly/by_query/{q}/mycosv/calls.vcf"
    if not os.path.exists(vcf):
        continue
    phylum = (r.get("phylum") or ".").strip() or "."
    cls = (r.get("class") or ".").strip() or "."
    ref = (r.get("benchmark_ref_fasta") or "").strip()
    if ref not in gc_cache:
        gc_cache[ref] = gc_of(ref)
    cg = gc_cache[ref]
    try:
        out = subprocess.run([REANN, vcf, f"{cg:.3f}", phylum, cls],
                             capture_output=True, text=True, timeout=600).stdout
    except Exception as e:
        sys.stderr.write(f"{q}: {e}\n")
        continue
    cnt = defaultdict(int)
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) >= 5:
            cnt[f[4]] += 1
    hs = cnt["HGT"] + cnt["STARSHIP"]
    te = sum(v for k, v in cnt.items() if k.startswith("TE_") or k == "REPEAT")
    rip = cnt["RIP"]
    for k, v in cnt.items():
        agg[k] += v
    per_genome.append((q, cls, cg, rip, hs, te, cnt.get("STARSHIP", 0)))
    done += 1
    if done % 25 == 0:
        sys.stderr.write(f"...{done} genomes\n")

with open(f"{OUT}/per_genome.tsv", "w") as fh:
    fh.write("query\tclass\tcladeGc\tRIP\tHGT_Starship\tTE_repeat\tSTARSHIP\n")
    for row in per_genome:
        fh.write("\t".join(str(x) for x in row) + "\n")

RIP = agg["RIP"]
HS = agg["HGT"] + agg["STARSHIP"]
TE = sum(v for k, v in agg.items() if k.startswith("TE_") or k == "REPEAT")
with open(f"{OUT}/aggregate.tsv", "w") as fh:
    fh.write("category\tcount\n")
    fh.write(f"n_genomes\t{done}\n")
    fh.write(f"RIP-like\t{RIP}\n")
    fh.write(f"HGT-/Starship-like\t{HS}\n")
    fh.write(f"TE/repeat-like\t{TE}\n")
    fh.write(f"STARSHIP_strict\t{agg['STARSHIP']}\n")
    for k in sorted(agg):
        fh.write(f"raw::{k}\t{agg[k]}\n")

print(f"genomes={done}")
print(f"RIP-like={RIP}  (old 164,314)")
print(f"HGT-/Starship-like={HS}  (old 58,540 cargo)")
print(f"TE/repeat-like={TE}")
print(f"STARSHIP strict={agg['STARSHIP']}")
print(f"wrote {OUT}/aggregate.tsv and per_genome.tsv")
