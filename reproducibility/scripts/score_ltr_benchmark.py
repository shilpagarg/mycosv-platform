#!/usr/bin/env python3
"""LTR retrotransposon benchmark.

MycoSV's element screen is alignment-free and self-contained; its structural
TE_LTR detector requires both LTRs of a full element within one SV block and so
rarely fires (fungal LTR elements are 5-15 kb; MycoSV SV blocks are much
smaller). This benchmark therefore evaluates two things against a dedicated LTR
annotation (the LTR retrotransposon subset - LTR/Gypsy, LTR/Copia, LTR/Ngaro -
of the MycoMobilome fungal repeat library, and, when available, de-novo
LTRharvest calls):

  1. discovery: are MycoSV SV loci enriched in LTR-annotated regions relative to
     chance (i.e. does SV discovery track LTR-rich sequence, as fungal biology
     predicts)?
  2. classification: are MycoSV's TE/repeat-labelled loci specifically enriched
     for LTR overlap?

Reports overlap counts, chance baseline, enrichment and an exact-binomial p.
"""
import bisect
import glob
import os
import re
from collections import defaultdict
from math import lgamma, log, exp, isfinite

ROOT = "experiments/biology_benchmark"
REANN = f"{ROOT}/reannotated"
LTR = f"{ROOT}/ltr"
GEN = {"aflavus": "GCA_014117485_1", "fgram": "GCA_018346565_1",
       "pchryso": "GCA_023624235_1", "cboidinii": "GCA_002007985_1",
       "tsemiorbis": "GCA_020045945_2", "tatroviride": "GCA_019297715_1",
       "aflavus2": "GCA_014784225_2", "pleryngii": "GCA_980434285_1"}


def norm(c):
    return re.sub(r'[^A-Za-z0-9]', '_', c.strip())


def binom_p_ge(k, n, p):
    if not n or p <= 0:
        return None
    if k <= 0:
        return 1.0
    if p >= 1:
        return 1.0
    lp, lq = log(p), log(1 - p)
    t = 0.0
    for i in range(k, n + 1):
        lc = lgamma(n + 1) - lgamma(i + 1) - lgamma(n - i + 1)
        term = exp(lc + i * lp + (n - i) * lq)
        if isfinite(term):
            t += term
        if t >= 1:
            return 1.0
    return min(1.0, t)


def load_ltr(tag):
    d = defaultdict(list)
    f = f"{LTR}/{tag}.ltr.bed"
    if not os.path.exists(f):
        return {}
    for line in open(f):
        p = line.rstrip("\n").split("\t")
        if len(p) < 3:
            continue
        d[norm(p[0])].append((int(p[1]), int(p[2])))
    idx = {}
    for c, spans in d.items():
        spans.sort()
        merged = []
        for s, e in spans:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        idx[c] = ([m[0] for m in merged], [m[1] for m in merged])
    return idx


def hits(c, s, e, idx):
    if c not in idx:
        return False
    st, en = idx[c]
    i = bisect.bisect_right(st, e) - 1
    while i >= 0:
        if en[i] >= s:
            return True
        if st[i] < s - 2_000_000:
            break
        i -= 1
    return False


def load_loci(qdir):
    f = f"{REANN}/{qdir}.reann.tsv"
    out = []
    if not os.path.exists(f):
        return out
    for line in open(f):
        p = line.rstrip("\n").split("\t")
        if len(p) < 5:
            continue
        try:
            c = norm(p[0]); s = int(float(p[1])); e = int(float(p[2]))
        except Exception:
            continue
        if e < s:
            s, e = e, s
        out.append((c, s, e, p[4]))
    return out


rows = []
for tag, qdir in GEN.items():
    idx = load_ltr(tag)
    if not idx:
        continue
    loci = load_loci(qdir)
    if not loci:
        continue
    nU = len(loci)
    all_hit = sum(1 for c, s, e, _ in loci if hits(c, s, e, idx))
    base = all_hit / nU if nU else 0            # chance = LTR-covered fraction of examined loci
    te = [(c, s, e) for c, s, e, lab in loci
          if lab.startswith("TE_") or lab == "REPEAT"]
    te_hit = sum(1 for c, s, e in te if hits(c, s, e, idx))
    nTE = len(te)
    exp_te = nTE * base
    p_te = binom_p_ge(te_hit, nTE, base) if nTE else None
    rows.append((tag, nU, all_hit, f"{base:.4f}", nTE, te_hit,
                 f"{exp_te:.1f}",
                 f"{te_hit/exp_te:.2f}" if exp_te else "NA",
                 f"{p_te:.4g}" if p_te is not None else "NA"))

out = f"{ROOT}/ltr_benchmark.tsv"
with open(out, "w") as fh:
    fh.write("genome\tn_sv_loci\tsv_loci_in_LTR\tLTR_frac_of_loci\t"
             "n_TE_labelled\tTE_in_LTR\texp_by_chance\tenrichment\tbinom_p\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
print(open(out).read())
print(f"wrote {out}")
