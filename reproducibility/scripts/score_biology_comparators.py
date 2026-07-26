import glob, re, os, csv, bisect
from collections import defaultdict

ROOT = "experiments/biology_benchmark"
PRED = "experiments/million_real/full_fungal_assembly_20260522_070100/assembly/by_query"
GEN = {"aflavus":    ("GCA_014117485_1", "GCA_014117485.1_ASM1411748v1"),
       "fgram":      ("GCA_018346565_1", "GCA_018346565.1_ASM1834656v1"),
       "pchryso":    ("GCA_023624235_1", "GCA_023624235.1_ASM2362423v1"),
       "cboidinii":  ("GCA_002007985_1", "GCA_002007985.1_ASM200798v1"),
       "tsemiorbis": ("GCA_020045945_2", "GCA_020045945.2_ASM2004594v2"),
       # Read-validation substitution assemblies (assembly-matched PacBio reads
       # exist for both); analysed as full panel members, reported separately
       # from the original five-genome panel.
       "tatroviride": ("GCA_019297715_1", "GCA_019297715.1_ASM1929771v1"),
       "aflavus2":    ("GCA_014784225_2", "GCA_014784225.2_ASM1478422v2"),
       # DToL Basidiomycota — natural Starship negative control (Ascomycota-only element)
       "pleryngii":   ("GCA_980434285_1", "GCA_980434285.1_gfPleEryn1.hap1.1")}

def comparator_ran(compname, tag, stem):
    """Did this comparator actually produce output for this genome?

    Detected from the filesystem rather than hardcoded, so a genome whose
    comparator job has not finished is reported as not_run instead of being
    scored against an empty interval set (which would read as precision 0).
    """
    if compname.startswith("RIPCAL"):
        return os.path.exists(os.path.join(ROOT, "ripcal", f"{tag}_ripcal_scan.gff"))
    if compname == "TEsorter":
        return bool(glob.glob(os.path.join(ROOT, "tesorter", f"{stem}*.gff3")) or
                    glob.glob(os.path.join(ROOT, "tesorter", f"{stem}*.cls.tsv")))
    if compname == "starfish":
        return bool(glob.glob(os.path.join(ROOT, "starfish", tag, "elementFinder", "*.bed")) or
                    glob.glob(os.path.join(ROOT, "starfish", tag, "geneFinder", f"{tag}.tyr.bed")))
    if compname == "OcculterCut":
        return os.path.exists(os.path.join(ROOT, "occultercut", tag, "regions.gff3"))
    if compname == "MycoMobilome":
        d = os.path.join(ROOT, "mycomobilome")
        return bool(glob.glob(os.path.join(d, f"{tag}.te.tsv")) or
                    glob.glob(os.path.join(d, f"{tag}.te.part*.tsv")))
    if compname == "starbase":
        d = os.path.join(ROOT, "starbase")
        return bool(glob.glob(os.path.join(d, f"{tag}.ship.tsv")) or
                    glob.glob(os.path.join(d, f"{tag}.ship.part*.tsv")))
    return False

def norm(c):
    return re.sub(r'[^A-Za-z0-9]', '_', c.strip())

# ---------------------------------------------------------------- MycoSV loci
def mycosv_loci(qdir):
    """Every SV locus MycoSV emitted for this query, with its element label.

    The universe is ALL rows, not just the annotated ones: precision/recall for
    an element class is only interpretable on the shared denominator of loci
    that MycoSV actually examined.
    """
    f = os.path.join(PRED, qdir, "biology_candidates.tsv")
    loci = {}
    if not os.path.exists(f):
        return loci
    with open(f) as fh:
        r = csv.reader(fh, delimiter='\t')
        hdr = next(r)
        ci = {h: i for i, h in enumerate(hdr)}
        for row in r:
            try:
                contig = norm(row[ci['query_contig']])
                s = int(float(row[ci['pos']])); e = int(float(row[ci['end']]))
                ec = row[ci['element_class']].strip().upper()
                sub = row[ci['mge_subtype']].strip().lower()
                ct = row[ci['candidate_type']].strip().lower()
            except Exception:
                continue
            if e < s:
                s, e = e, s
            key = (contig, s, e)
            prev = loci.get(key)
            labels = prev[0] if prev else set()
            if ec == 'RIP':
                labels.add('RIP')
            if ec.startswith('TE_') or ec == 'REPEAT':
                labels.add('TE')
            if ec == 'STARSHIP':
                labels.add('STARSHIP_STRICT'); labels.add('STARSHIP_WIDE')
            # The caller labels most Starship-like mega-elements as HGT with an
            # integrative subtype; the strict STARSHIP class fires only on the
            # full AT-hull + GC-cargo signature. Score both.
            if ec == 'HGT' and (sub == 'integrative' or ct == 'hgt_candidate'):
                labels.add('STARSHIP_WIDE')
            loci[key] = (labels,)
    return loci

REANN_DIR = os.path.join(ROOT, "reannotated")

def mycosv_loci_reann(qdir):
    """Post-fix element labels from the reannotated calls (RIP-reorder +
    Pezizomycotina Starship gate). element_class is authoritative here, so the
    strict STARSHIP class is used directly rather than the HGT proxy."""
    f = os.path.join(REANN_DIR, f"{qdir}.reann.tsv")
    loci = {}
    if not os.path.exists(f):
        return None
    for line in open(f):
        p = line.rstrip("\n").split("\t")
        if len(p) < 5:
            continue
        try:
            contig = norm(p[0]); s = int(float(p[1])); e = int(float(p[2]))
        except Exception:
            continue
        ec = p[4].strip().upper()
        if e < s:
            s, e = e, s
        key = (contig, s, e)
        labels = loci.get(key, set())
        if ec == "RIP":
            labels.add("RIP")
        if ec.startswith("TE_") or ec == "REPEAT":
            labels.add("TE")
        if ec == "STARSHIP":
            labels.add("STARSHIP_STRICT"); labels.add("STARSHIP_WIDE")
        loci[key] = labels
    return {k: (v,) for k, v in loci.items()}


# ------------------------------------------------------------ comparator sets
def read_fasta(path):
    seqs = {}
    name = None
    buf = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if name:
                    seqs[name] = "".join(buf).upper()
                name = norm(line[1:].split()[0])
                buf = []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs

def rip_indices(s):
    """Selker RIP product and substrate indices for one window."""
    d = defaultdict(int)
    for i in range(len(s) - 1):
        d[s[i:i + 2]] += 1
    prod_n, prod_d = d['TA'], d['AT']
    sub_n, sub_d = d['CA'] + d['TG'], d['AC'] + d['GT']
    return (prod_n / prod_d if prod_d else 0.0,
            sub_n / sub_d if sub_d else 0.0)

def ripcal_intervals(tag, stem, stringent=False):
    """RIPCAL index-scan regions.

    This RIPCAL build ignores the -q/-w threshold flags (a scan run with the
    stringent thresholds returns a byte-identical region set), so the stringent
    tier is applied post hoc: RIPCAL's own regions are re-scored with the
    standard Selker indices and kept only at TpA/ApT >= 1.1 and
    (CpA+TpG)/(ApC+GpT) <= 0.9. The comparator call set is still RIPCAL's; only
    the stringency cut is ours, and both tiers are reported.
    """
    f = os.path.join(ROOT, "ripcal", f"{tag}_ripcal_scan.gff")
    iv = []
    if not os.path.exists(f):
        return iv
    for line in open(f):
        if line.startswith('#'):
            continue
        p = line.rstrip('\n').split('\t')
        if len(p) < 5 or p[2] != 'region':
            continue
        try:
            iv.append((norm(p[0]), int(p[3]), int(p[4])))
        except Exception:
            pass
    if not stringent:
        return iv
    seqs = read_fasta(os.path.join(ROOT, "genomes", f"{stem}.fna"))
    keep = []
    for c, s, e in iv:
        sub = seqs.get(c, "")[s - 1:e]
        if len(sub) < 100:
            continue
        prod, subst = rip_indices(sub)
        if prod >= 1.1 and subst <= 0.9:
            keep.append((c, s, e))
    return keep

def tesorter_intervals(stem):
    iv = []
    for f in sorted(set(glob.glob(os.path.join(ROOT, "tesorter", f"{stem}*.gff3")))):
        for line in open(f):
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 5:
                continue
            try:
                iv.append((norm(p[0]), int(p[3]), int(p[4])))
            except Exception:
                pass
    if not iv:
        # genome-mode fallback: TEsorter encodes coords in the sequence id
        for f in glob.glob(os.path.join(ROOT, "tesorter", f"{stem}*.cls.tsv")):
            for line in open(f):
                m = re.match(r'([^\s:]+):(\d+)\.\.(\d+)', line)
                if m:
                    try:
                        iv.append((norm(m.group(1)), int(m.group(2)), int(m.group(3))))
                    except Exception:
                        pass
    return iv

def occultercut_intervals(tag, stem):
    """AT-rich regions from OcculterCut's genome segmentation.

    OcculterCut segments the genome by local GC content and reports the GC
    cutoff separating the AT-rich mode from the GC-equilibrated mode. In most
    fungal taxa the AT-rich mode is the signature of RIP. The comparator set is
    the segments falling below that run's own reported cutoff.
    """
    d = os.path.join(ROOT, "occultercut", tag)
    gff = os.path.join(d, "regions.gff3")
    log = os.path.join(d, "run.log")
    if not (os.path.exists(gff) and os.path.exists(log)):
        return []
    m = re.search(r'cutoff:\s*([0-9.]+)', open(log).read())
    if not m:
        return []
    cutoff = float(m.group(1))
    seqs = read_fasta(os.path.join(ROOT, "genomes", f"{stem}.fna"))
    keep = []
    for line in open(gff):
        p = line.rstrip('\n').split('\t')
        if len(p) < 5 or p[2] != 'region':
            continue
        try:
            c = norm(p[0]); s = int(p[3]); e = int(p[4])
        except Exception:
            continue
        sub = seqs.get(c, "")[s - 1:e]
        if len(sub) < 100:
            continue
        gc = sum(sub.count(x) for x in "GC")
        at = sum(sub.count(x) for x in "AT")
        if not (gc + at):
            continue
        if 100.0 * gc / (gc + at) < cutoff:
            keep.append((c, s, e))
    return keep

def _blast_intervals(path, min_len, min_pident):
    """`path` may be a single file or the stem of chunked `.partN.` files."""
    """Genome-coordinate intervals from a blastn -outfmt 6 run with the genome
    as query: qseqid qstart qend sseqid pident length evalue bitscore."""
    iv = []
    files = [path] if os.path.exists(path) else []
    stem, ext = os.path.splitext(path)
    files += sorted(glob.glob(f"{stem}.part*{ext}"))
    if not files:
        return iv
    for line in (ln for f in files for ln in open(f)):
        p = line.rstrip('\n').split('\t')
        if len(p) < 6:
            continue
        try:
            c = norm(p[0]); s = int(p[1]); e = int(p[2])
            pid = float(p[4]); ln = int(p[5])
        except Exception:
            continue
        if ln < min_len or pid < min_pident:
            continue
        if e < s:
            s, e = e, s
        iv.append((c, s, e))
    return iv

def mycomobilome_intervals(tag):
    """TE copies called by homology to the MycoMobilome v1.1 fungal TE library."""
    return _blast_intervals(os.path.join(ROOT, "mycomobilome", f"{tag}.te.tsv"),
                            min_len=100, min_pident=70.0)

def starbase_intervals(tag):
    """Starship-homologous regions from the starbase curated Starship set.
    Starships are 20-700 kb; a 5 kb minimum keeps short cargo/TE homology out."""
    return _blast_intervals(os.path.join(ROOT, "starbase", f"{tag}.ship.tsv"),
                            min_len=5000, min_pident=70.0)

def starfish_intervals(tag):
    """starfish Starship evidence for a single genome.

    starfish's `insert`/`flank` element-finding compares a genome carrying an
    element against sister genomes lacking it, so on single-genome runs it
    exits with "0 tyr captains with candidate insertion sites" and writes no
    elementFinder BED. We therefore fall back to the tyrosine recombinase
    ("captain") genes from the `annotate`+`sketch` stage, which are the
    diagnostic Starship feature and are genuinely starfish-derived. Captain
    counts track the expected phylogenetic distribution (0 in the
    Saccharomycotina yeast, 2-10 in the Pezizomycotina genomes).
    """
    files = glob.glob(os.path.join(ROOT, "starfish", tag, "elementFinder", "*.bed"))
    if not files:
        files = glob.glob(os.path.join(ROOT, "starfish", tag, "geneFinder", f"{tag}.tyr.bed"))
    iv = []
    for f in files:
        for line in open(f):
            if line.startswith('#') or line.startswith('track'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 3:
                continue
            try:
                iv.append((norm(p[0]), int(p[1]), int(p[2])))
            except Exception:
                pass
    return iv

# ---------------------------------------------------------------- interval ops
def merge_by_contig(iv):
    d = defaultdict(list)
    for c, s, e in iv:
        if e < s:
            s, e = e, s
        d[c].append((s, e))
    out = {}
    for c, spans in d.items():
        spans.sort()
        merged = []
        for s, e in spans:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[c] = ([m[0] for m in merged], [m[1] for m in merged])
    return out

def hits(locus, idx):
    c, s, e = locus
    if c not in idx:
        return False
    starts, ends = idx[c]
    i = bisect.bisect_right(starts, e) - 1
    while i >= 0:
        if ends[i] >= s:
            return True
        if starts[i] < s - 5_000_000:
            break
        i -= 1
    return False

def binom_p_ge(k, n, p):
    """P(X >= k) for X~Bin(n,p): is observed agreement above the chance rate?

    Computed in log space; comb(n, i) overflows a float for the label counts
    seen here (n up to ~3000).
    """
    from math import lgamma, log, exp, isfinite
    if not n or p is None:
        return None
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0
    lp, lq = log(p), log(1.0 - p)
    total = 0.0
    for i in range(k, n + 1):
        lc = lgamma(n + 1) - lgamma(i + 1) - lgamma(n - i + 1)
        term = exp(lc + i * lp + (n - i) * lq)
        if not isfinite(term):
            continue
        total += term
        if total >= 1.0:
            return 1.0
    return min(1.0, total)

def wilson(k, n):
    if not n:
        return ("NA", "NA")
    from math import sqrt
    z = 1.959963985
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return ("%.4f" % max(0.0, c - h), "%.4f" % min(1.0, c + h))

# ---------------------------------------------------------------------- score
CLASSES = [
    ('RIP',             'RIPCAL',            'RIP'),
    ('RIP',             'RIPCAL_stringent',  'RIP_STRINGENT'),
    ('RIP',             'OcculterCut',       'RIP_ATRICH'),
    ('TE',              'TEsorter',          'TE'),
    ('TE',              'MycoMobilome',      'TE_LIB'),
    ('STARSHIP_STRICT', 'starfish',          'STARSHIP'),
    ('STARSHIP_WIDE',   'starfish',          'STARSHIP'),
    ('STARSHIP_STRICT', 'starbase',          'SHIP_LIB'),
    ('STARSHIP_WIDE',   'starbase',          'SHIP_LIB'),
]

rows = []
for tag, (qdir, stem) in sorted(GEN.items()):
    loci = mycosv_loci_reann(qdir) or mycosv_loci(qdir)
    universe = sorted(loci.keys())
    comp_raw = {'RIP': ripcal_intervals(tag, stem),
                'RIP_STRINGENT': ripcal_intervals(tag, stem, stringent=True),
                'RIP_ATRICH': occultercut_intervals(tag, stem),
                'TE': tesorter_intervals(stem),
                'TE_LIB': mycomobilome_intervals(tag),
                'STARSHIP': starfish_intervals(tag),
                'SHIP_LIB': starbase_intervals(tag)}
    comp_idx = {k: merge_by_contig(v) for k, v in comp_raw.items()}

    for grp, compname, compkey in CLASSES:
        if not comparator_ran(compname, tag, stem):
            rows.append(dict(
                genome=tag, element_class=grp, comparator=compname,
                n_comparator_intervals=0, n_sv_loci_examined=len(universe),
                n_mycosv_labelled=sum(1 for k in universe if grp in loci[k][0]),
                n_comparator_supported_loci="NA", tp="NA",
                precision="not_run", precision_ci95="NA",
                recall="not_run", recall_ci95="NA", f1="NA",
                chance_precision="NA", enrichment_over_chance="NA",
                expected_tp_by_chance="NA", binom_p="NA", power_verdict="NA"))
            continue
        idx = comp_idx[compkey]
        n_comp_raw = len(comp_raw[compkey])
        M = [k for k in universe if grp in loci[k][0]]
        C = [k for k in universe if hits(k, idx)]
        tp = sum(1 for k in M if hits(k, idx))
        nM, nC, nU = len(M), len(C), len(universe)
        prec = tp / nM if nM else None
        rec = tp / nC if nC else None
        f1 = (2 * prec * rec / (prec + rec)) if (prec and rec) else (0.0 if (prec is not None and rec is not None) else None)
        # chance precision: fraction of examined loci the comparator covers
        base = nC / nU if nU else None
        enr = (prec / base) if (prec is not None and base) else None
        # Power guard: with few labels and a low chance rate, a precision of 0
        # is uninformative rather than evidence of disagreement.
        exp_tp = nM * base if base is not None else None
        pval = binom_p_ge(tp, nM, base)
        # The test is an exact binomial, so a small expected count does not
        # invalidate a positive result; it only makes a NULL result
        # uninformative. Significance is therefore checked first.
        if exp_tp is None or pval is None:
            power = "NA"
        elif pval < 0.05:
            power = "enriched_p<0.05"
        elif exp_tp < 3:
            power = "underpowered"
        else:
            power = "not_distinguishable"
        pl, ph = wilson(tp, nM)
        rl, rh = wilson(tp, nC)
        rows.append(dict(
            genome=tag, element_class=grp, comparator=compname,
            n_comparator_intervals=n_comp_raw, n_sv_loci_examined=nU,
            n_mycosv_labelled=nM, n_comparator_supported_loci=nC, tp=tp,
            precision=("%.4f" % prec) if prec is not None else "NA",
            precision_ci95=f"{pl}-{ph}",
            recall=("%.4f" % rec) if rec is not None else "NA",
            recall_ci95=f"{rl}-{rh}",
            f1=("%.4f" % f1) if f1 is not None else "NA",
            chance_precision=("%.4f" % base) if base is not None else "NA",
            enrichment_over_chance=("%.2f" % enr) if enr is not None else "NA",
            expected_tp_by_chance=("%.1f" % exp_tp) if exp_tp is not None else "NA",
            binom_p=("%.4g" % pval) if pval is not None else "NA",
            power_verdict=power,
        ))

os.makedirs(ROOT, exist_ok=True)
out = os.path.join(ROOT, "biology_comparator_pr.tsv")
cols = ["genome", "element_class", "comparator", "n_comparator_intervals",
        "n_sv_loci_examined", "n_mycosv_labelled", "n_comparator_supported_loci",
        "tp", "precision", "precision_ci95", "recall", "recall_ci95", "f1",
        "chance_precision", "enrichment_over_chance", "expected_tp_by_chance",
        "binom_p", "power_verdict"]
with open(out, 'w') as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")
print("WROTE", out)
print(open(out).read())
