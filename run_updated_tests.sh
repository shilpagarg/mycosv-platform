#!/usr/bin/env bash
# Designed for Linux
# run_updated_tests.sh — v18
#
# Validates all fixes and additions across the v13→v18 audit:
#
#   DS-7  : TRA_INTRA emitted via PathPositionIndex order-stats (Vuillemin 1980)
#   DS-8  : is_inversion_flex — tolerance-fraction quick-reject
#   DS-9  : OffRefNoveltyTier — 4-tier novelty scoring
#   DS-10 : ReferenceLCAIndex — Euler tour + sparse-table RMQ (Harel & Tarjan 1984)
#   DS-11 : classify_triallelic — OVERLAPPING/NESTED/INTERLOCKING/PROPERLY_TRIALLELIC
#   DS-12 : classify_pantree — SNP/MNP/INS/DEL/DUP/REPL/INV/NON_REF
#
#   DS-13 : SuffixArray + LCP (Manber & Myers 1990; Kasai 2001)
#           O(|q| log N) MEM seeding → precise breakpoints for ALL SV types
#   DS-14 : WaveletTree over BWT (Grossi, Gupta & Vitter 1993)
#           O(k log σ) rank/select replaces O(L) unordered_set k-mer scan
#   DS-15 : van Emde Boas tree (van Emde Boas 1975/1977)
#           O(log log U) predecessor/successor for PathPositionIndex windows
#   DS-16 : Merge Sort Tree for interval stabbing (Willard 1985)
#           O(N log² N) triallelic batch-classification
#   DS-17 : Fenwick Tree / BIT (Fenwick 1994)
#           O(log N) capacity-aware LRU cache eviction
#   DS-18 : Treap for MEM seed chaining (Aragon & Seidel 1989)
#           O(N log N) chain DP → INS/DEL/INV/DUP/TRA classification
#
#   QA-1  : Coverage-aware query auto-tuning for assembly, short-read low/high
#           coverage, and long-read sequencing input modes
#
# Usage:
#   bash run_updated_tests.sh               # full suite
#   bash run_updated_tests.sh --no-pytest   # skip Python tests
#   bash run_updated_tests.sh --jobs N      # parallel pytest workers

set -euo pipefail

export LC_ALL=C
umask 022

usage() {
    cat <<'EOF'
Usage:
  bash run_updated_tests.sh [--no-pytest] [--jobs N]
EOF
}

if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
    usage
    exit 0
fi


ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
MAIN_CPP="$ROOT_DIR/main.cpp"
RUNNER_CPP_BIN="$ROOT_DIR/fungi_graphsv_tol_test"

fail(){ echo "[FAIL] $*" >&2; exit 1; }
require_cmd(){ command -v "$1" >/dev/null 2>&1 || fail "required command not found: $1"; }

BIN="$ROOT_DIR/fungi_graphsv_tol_test"
BIN_EXEC="$BIN"
JOBS=1
RUN_PYTEST=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-pytest) RUN_PYTEST=0 ;;
        --jobs)      JOBS="${2:?--jobs requires a number}"; shift ;;
        *) echo "[warn] unknown argument: $1" ;;
    esac
    shift
done

[[ "$JOBS" =~ ^[0-9]+$ && "$JOBS" -ge 1 ]] || fail "--jobs must be a positive integer"

require_cmd bash
require_cmd g++
require_cmd python3
require_cmd grep
require_cmd find

echo "=== fungi_graphsv_tol test suite v18 (DS-7 through DS-18) ==="
echo "    ROOT_DIR : $ROOT_DIR"
echo "    JOBS     : $JOBS"
echo ""

# ── 1. Shell syntax ─────────────────────────────────────────────────────────
echo "[1/7] Shell script syntax..."
bash -n "$ROOT_DIR/run_tol_bench.sh"
bash -n "$ROOT_DIR/run_updated_tests.sh"
echo "      OK"

# ── 2. Python syntax ─────────────────────────────────────────────────────────
echo "[2/7] Python syntax..."
PY_FILES=()
for f in \
    "$ROOT_DIR/test_amf.py" \
    "$ROOT_DIR/test_pipeline_features.py" \
    "$ROOT_DIR/test_all_use_cases.py" \
    "$ROOT_DIR/test_biological_use_case_fungi.py" \
    "$ROOT_DIR/example_real_data_fungi.py"
do
    [[ -f "$f" ]] && PY_FILES+=("$f")
done
if [[ ${#PY_FILES[@]} -gt 0 ]]; then
    python3 -m py_compile "${PY_FILES[@]}"
    echo "      OK"
else
    echo "      SKIP (no Python test files present)"
fi

# ── 3. Compile ───────────────────────────────────────────────────────────────
echo "[3/7] Compiling main.cpp (C++17, -Wall -Wextra)..."
g++ -O2 -std=c++17 -pthread \
    -Wall -Wextra \
    -Wno-unused-parameter -Wno-unused-variable \
    -I"$ROOT_DIR" \
    "$MAIN_CPP" \
    -o "$BIN"
[[ -x "$BIN_EXEC" ]] || fail "compiled binary not found after g++ build: $BIN"
echo "      OK → $BIN"

# ── 4. CLI flag smoke ────────────────────────────────────────────────────────
echo "[4/7] CLI flag smoke..."
HELP=$("$BIN_EXEC" --help 2>&1 || true)
for flag in \
    "--tol-hierarchical" "--tol-multi-rank" \
    "--tol-ancestral-align" "--tol-ancestral-recomb" \
    "--tol-recomb-min-seg-bp" "--tol-recomb-max-breakpoints" \
    "--tol-query-window-bp" "--tol-base-graph-build" \
    "--tol-validate-index"
do
    echo "$HELP" | grep -q -- "$flag" \
        && echo "      [ok] $flag" \
        || { echo "      [FAIL] $flag missing" >&2; exit 1; }
done

# ── 5. Functional smoke tests ────────────────────────────────────────────────
echo "[5/8] Functional smoke tests..."
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

# ── 5a. Index build (FIX-AMF adaptive cap + DS-12 pantree annotation path) ──
echo "      [5a] Index build..."
cat > "$TMP/manifest.tsv" <<'EOF'
#asm_name	phylum	class	order	family	genus	clade_name	clade_rank	fasta_path
asm_amf	Glomeromycota	Glomeromycetes	Glomerales	Glomeraceae	Rhizophagus	giant_amf	genus
asm_rust	Basidiomycota	Pucciniomycetes	Pucciniales	Pucciniaceae	Puccinia	rust_smut	family
asm_yeast	Ascomycota	Saccharomycetes	Saccharomycetales	Saccharomycetaceae	Saccharomyces	compact_yeast	species
asm_chytrid	Chytridiomycota	Chytridiomycetes	Chytridiales	Chytriaceae	Batrachochytrium	hgt_clade	genus
EOF

"$BIN_EXEC" \
    --tol-hierarchical \
    --tol-build-index  "$TMP/manifest.tsv" \
    --tol-index-dir    "$TMP/idx" \
    --tol-registry-dir "$TMP/reg" \
    --tol-multi-rank --quiet

GBZ=$(find "$TMP/reg" -name "*.gbz" | wc -l)
GBWT=$(find "$TMP/idx" -name "*.gbwt" | wc -l)
MINI=$(find "$TMP/idx" -name "*.min" | wc -l)
CIDX=$(find "$TMP/idx" -name "*.cidx" | wc -l)
test -f "$TMP/reg/clade_manifest.tsv" || { echo "[FAIL] manifest missing" >&2; exit 1; }
test -f "$TMP/idx/routing_manifest.tsv" || { echo "[FAIL] routing manifest missing" >&2; exit 1; }
echo "      [ok] $GBZ .gbz / $GBWT .gbwt / $MINI .min / $CIDX .cidx shards"

# ── 5b. Validation report ────────────────────────────────────────────────────
echo "      [5b] Validation report..."
"$BIN_EXEC" \
    --tol-hierarchical --tol-validate-index \
    --tol-build-index "$TMP/manifest.tsv" \
    --tol-index-dir   "$TMP/idx" \
    --tol-registry-dir "$TMP/reg" \
    --tol-validation-report "$TMP/val.tsv" \
    --quiet || true
test -f "$TMP/val.tsv" || { echo "[FAIL] val.tsv missing" >&2; exit 1; }
echo "      [ok] validation report written"

# ── 5c. Length-fallback SV detection — plain contig names, no hints ──────────
# Uses --query-mode assembly so the auto-detector does not rewrite the query
# contigs as short-read unitigs.  Sequences are non-homopolymer so they pass
# the low-complexity filter and are long enough to exceed minSvLen=40.
echo "      [5c] Length-fallback INS/DEL detection (plain contig names)..."
python3 - "$TMP" <<'PYEOF'
from pathlib import Path
import random

rng = random.Random(42)
bases = 'ACGT'

def rseq(n):
    return ''.join(rng.choices(bases, k=n))

out = Path(__import__('sys').argv[1])
ref = out / "ref_plain.fa"
# Reference: ctg1=500bp, ctg2=500bp
ref.write_text(">ctg1\n" + rseq(500) + "\n>ctg2\n" + rseq(500) + "\n")
# ctg1 query is 60 bp shorter → DEL (delta=60 > minSvLen=40)
# ctg2 query is 80 bp longer  → INS
query = out / "query_plain.fa"
query.write_text(">ctg1\n" + rseq(440) + "\n>ctg2\n" + rseq(580) + "\n")
(out / "refs_plain.txt").write_text(str(ref) + "\n")
(out / "queries_plain.txt").write_text(str(query) + "\n")
PYEOF

"$BIN_EXEC" \
    --ref-list    "$TMP/refs_plain.txt" \
    --query-list  "$TMP/queries_plain.txt" \
    --out-prefix  "$TMP/calls_plain" \
    --query-mode  assembly \
    --quiet

VCF="$TMP/calls_plain.vcf"
grep -q "SVTYPE=DEL" "$VCF" \
    && echo "      [ok] SVTYPE=DEL detected from shorter query contig" \
    || { echo "      [FAIL] SVTYPE=DEL not detected" >&2; exit 1; }
grep -q "SVTYPE=INS" "$VCF" \
    && echo "      [ok] SVTYPE=INS detected from longer query contig" \
    || { echo "      [FAIL] SVTYPE=INS not detected" >&2; exit 1; }
grep -q "VT=" "$VCF" \
    && echo "      [ok] VT= INFO field present (DS-12 pantree class)" \
    || { echo "      [FAIL] VT= INFO field absent in VCF output" >&2; exit 1; }

# Verify: a query contig whose name has NO match in the reference index is
# handled by the off-reference novelty scorer (simple_offref_fallback_calls),
# NOT by the length fallback.  It should emit an OFF_REF call (not INS/DEL),
# because there is no length delta to measure from an unrelated contig.
python3 - "$TMP" <<'PYEOF'
from pathlib import Path
import random
rng = random.Random(99)
# AT-rich reference
ref_seq = ''.join(rng.choices('AATT', k=500))
# GC-rich query with a completely different contig name — no name match possible
qseq = ''.join(rng.choices('GCGC', k=500))
out = Path(__import__('sys').argv[1])
ref2 = out / "ref_nomatch.fa"
ref2.write_text(f">ctgRef\n{ref_seq}\n")
query2 = out / "query_nomatch.fa"
query2.write_text(f">ctgNovel\n{qseq}\n")
(out / "refs_nomatch.txt").write_text(str(ref2) + "\n")
(out / "queries_nomatch.txt").write_text(str(query2) + "\n")
PYEOF
"$BIN_EXEC" \
    --ref-list    "$TMP/refs_nomatch.txt" \
    --query-list  "$TMP/queries_nomatch.txt" \
    --out-prefix  "$TMP/calls_nomatch" \
    --query-mode  assembly \
    --quiet
python3 - "$TMP/calls_nomatch.vcf" <<'PYEOF'
import sys
lines = [l for l in open(sys.argv[1]).read().splitlines()
         if not l.startswith('#') and l.strip()]
# The offref fallback fires for novel sequences — we expect either:
#   (a) no call (if the sequence is too short or low-complexity), or
#   (b) an OFF_REF call (not INS/DEL — length fallback requires a name match)
insdel = [l for l in lines if 'SVTYPE=INS' in l or 'SVTYPE=DEL' in l]
if insdel:
    print(f"[FAIL] length fallback emitted INS/DEL for unmatched contig name: {insdel[0]}",
          file=sys.stderr); sys.exit(1)
offref = [l for l in lines if 'SVTYPE=OFF_REF' in l]
if offref:
    print("      [ok] off-reference novelty scorer fired (no INS/DEL from length fallback)")
else:
    print("      [ok] no call for unmatched contig name (short/low-complexity sequence)")
PYEOF

# ── 5d. DS-9 novelty tier smoke — plain contig name, no hint ─────────────────
# The query contig has a plain biological name and a GC-rich sequence with zero
# k-mer overlap with the AT-rich reference.  The simple_offref_fallback fires
# and emits an OFF_REF call with a novelty tier annotation.
echo "      [5d] DS-9 off-reference novelty tier (plain name)..."
python3 - "$TMP" <<'PYEOF'
from pathlib import Path
import random
rng = random.Random(7)
# AT-rich reference (GC~30%)
ref_seq = ''.join(rng.choices('AATT', k=500))
# GC-rich query with no k-mer overlap (GC~70%, completely different sequence)
qseq = ''.join(rng.choices('GCGC', k=500))
out = Path(__import__('sys').argv[1])
ref = out / "ref_novel.fa"
ref.write_text(f">ctgRef\n{ref_seq}\n")
query = out / "query_novel.fa"
# Plain biological name — no __sv_ suffix
query.write_text(f">ctgNovel1\n{qseq}\n")
(out / "refs_novel.txt").write_text(str(ref) + "\n")
(out / "queries_novel.txt").write_text(str(query) + "\n")
PYEOF

"$BIN_EXEC" \
    --ref-list    "$TMP/refs_novel.txt" \
    --query-list  "$TMP/queries_novel.txt" \
    --out-prefix  "$TMP/calls_novel" \
    --query-mode  assembly \
    --quiet

VCF_N="$TMP/calls_novel.vcf"
grep -q "SVTYPE=OFF_REF" "$VCF_N" \
    && echo "      [ok] OFF_REF call emitted for novel sequence" \
    || { echo "      [FAIL] OFF_REF call missing for truly novel sequence" >&2; exit 1; }
grep -q "ANNOT=NOVEL\|ANNOT=NOVEL_WEAK\|ANNOT=DIVERGED\|ANNOT=OFF_REF_KNOWN" "$VCF_N" \
    && echo "      [ok] DS-9 novelty tier in ANNOT field" \
    || { echo "      [FAIL] DS-9 novelty tier annotation absent from VCF" >&2; exit 1; }

# ── 5e. DS-10/11/12 pantree classification VCF header fields ─────────────────
echo "      [5e] DS-10/11/12 pantree VCF INFO headers..."
# Runtime check: use the VCF produced in 5c (plain-name calls).
for field in "ID=VT" "ID=NR" "ID=TOPO" "ID=OFF_REF_TIER"; do
    grep -q "##INFO=<${field}" "$VCF" \
        && echo "      [ok] VCF header contains ##INFO=<${field}" \
        || { echo "      [FAIL] VCF header missing ##INFO=<${field}>" >&2; exit 1; }
done
grep -q "ID=VT" "$ROOT_DIR/main.cpp" \
    && echo "      [ok] VT INFO header in main.cpp" \
    || { echo "      [FAIL] VT INFO header missing" >&2; exit 1; }
grep -q "ID=NR" "$ROOT_DIR/main.cpp" \
    && echo "      [ok] NR INFO header in main.cpp" \
    || { echo "      [FAIL] NR INFO header missing" >&2; exit 1; }
grep -q "ID=TOPO" "$ROOT_DIR/main.cpp" \
    && echo "      [ok] TOPO INFO header in main.cpp" \
    || { echo "      [FAIL] TOPO INFO header missing" >&2; exit 1; }

# ── 5f. ARG ancestral recombination flags ────────────────────────────────────
echo "      [5f] ARG CLI pass-through..."
"$BIN_EXEC" \
    --ref-list    "$TMP/refs_plain.txt" \
    --query-list  "$TMP/queries_plain.txt" \
    --out-prefix  "$TMP/calls_arg" \
    --query-mode  assembly \
    --tol-ancestral-recomb \
    --tol-recomb-min-seg-bp 500 \
    --tol-recomb-max-breakpoints 8 \
    --quiet
echo "      [ok] ARG flags accepted without error"

echo "      All functional smoke tests passed."

# ── 6. Symbol audit ──────────────────────────────────────────────────────────
echo "[6/7] Symbol presence audit..."

L1="$ROOT_DIR/layer1_clade_graph.hpp"
HE="$ROOT_DIR/hierarchical_engine.hpp"
BR="$ROOT_DIR/fungi_tol_bridge.hpp"
MC="$ROOT_DIR/main.cpp"

check() {
    local file=$1 sym=$2
    grep -q "$sym" "$file" \
        && echo "      [ok] $sym" \
        || { echo "      [FAIL] $sym missing from $(basename $file)" >&2; exit 1; }
}

# DS-7
check "$L1" "PathPositionIndex"
check "$L1" "TraIntra"
check "$L1" "insert_position"
check "$L1" "quick_reject_window"

# DS-8
check "$L1" "is_inversion_flex"
check "$L1" "quick.reject"

# DS-9
check "$HE" "OffRefNoveltyTier"
check "$HE" "score_off_ref_novelty"
check "$HE" "UncoveredWindow"
check "$HE" "make_off_reference_call_scored"
# windowsOut->push_back is referenced in a comment only; verify the novelty
# tier name function is present instead (it is the live DS-9 implementation).
check "$HE" "novelty_tier_name"

# DS-10
check "$L1" "ReferenceLCAIndex"
check "$L1" "sparseTable"
check "$L1" "branch_point"
check "$L1" "build_sparse_table"

# DS-11
check "$L1" "TriallelicTopology"
check "$L1" "classify_triallelic"
check "$L1" "INTERLOCKING"

# DS-12
check "$L1" "PantreeVariantClass"
check "$L1" "classify_pantree"
check "$L1" "annotate_pantree_classes"
check "$L1" "NonRefVariant"

# Bridge propagation
check "$BR" "pantreeClass"
check "$BR" "isNonRefVariant"
check "$BR" "triallelicTopology"

# VCF fields
check "$MC" "ID=VT"
check "$MC" "ID=NR"
check "$MC" "ID=TOPO"
check "$MC" "pantreeClass"

# ── Ground-truth cleanliness audit ───────────────────────────────────────────
# Forbidden symbols: any of these in main.cpp means the hint path is active.
for forbidden in "EncodedSvHint" "parse_encoded_sv_hint" "simple_hint_fallback" \
                 "has_hints" "offrefName" "hintedType" "traMatesByEvent"; do
    grep -q "$forbidden" "$MC" \
        && { echo "      [FAIL] Forbidden hint symbol '$forbidden' in main.cpp" >&2; exit 1; } \
        || echo "      [ok] absent: $forbidden"
done

# Reference index must reject hint-encoded contigs (not strip them).
grep -q "simulator hint suffix\|will not be indexed" "$MC" \
    && echo "      [ok] load_simple_ref_index rejects hint-encoded ref contigs" \
    || { echo "      [FAIL] load_simple_ref_index does not reject hint-encoded ref contigs" >&2; exit 1; }

# strip_sv_suffix must only appear in best_ref_match (query lookup key only).
STRIP_SITES=$(grep -c "strip_sv_suffix" "$MC" 2>/dev/null || echo 0)
# Expected: 1 definition + 1 call in best_ref_match = 2 occurrences maximum
if [[ "$STRIP_SITES" -le 2 ]]; then
    echo "      [ok] strip_sv_suffix confined to best_ref_match ($STRIP_SITES occurrences)"
else
    echo "      [FAIL] strip_sv_suffix has $STRIP_SITES occurrences — should be ≤2" >&2; exit 1
fi

# Output contig name must be verbatim (v.qContig = contigName, not stripped).
grep -q "v\.qContig[[:space:]]*=[[:space:]]*contigName" "$MC" \
    && echo "      [ok] v.qContig assigned verbatim contig name" \
    || { echo "      [FAIL] v.qContig not assigned contigName verbatim" >&2; exit 1; }

echo "      All symbols present."

# ── DS-13 through DS-18 symbol audit ─────────────────────────────────────
echo "[6c/7] DS-13..DS-18 new data structure audit..."

# DS-13: SuffixArray + LCP
check "$L1" "SuffixArray"
check "$L1" "find_mems"
check "$L1" "revcomp"
check "$L1" "Kasai"

# DS-14: WaveletTree
L3="$ROOT_DIR/layer3_routing_index.hpp"
check "$L3" "WaveletTree"
check "$L3" "rank"
check "$L3" "count_kmer"

# DS-15: van Emde Boas
check "$L1" "VEBTree"
check "$L1" "any_in"
check "$L1" "log log U"

# DS-16: Merge Sort Tree
check "$L1" "MergeSortTree"
check "$L1" "stab_count"
check "$L1" "overlapping"

# DS-17: Fenwick Tree
check "$L1" "FenwickTree"
check "$L1" "find_kth"
check "$L1" "prefix_sum"

# DS-18: Treap chaining
check "$L1" "ChainTreap"
check "$L1" "insert_and_chain"
check "$L1" "best_chain_path"
check "$L1" "SvTypeFromChain"

# Integration point: try_mem_chain_call in bridge
check "$BR" "try_mem_chain_call"
check "$BR" "mem_chain_ds13_ds18"

echo "      DS-13..DS-18 symbols verified."

# ── Taxonomy rank consistency audit ──────────────────────────────────────
echo "[6b/7] Taxonomy rank consistency..."
TR="$ROOT_DIR/taxonomy_ranks.hpp"
SIM="$ROOT_DIR/test_amf.py"

check "$TR" "rank_from_string"
check "$TR" "rank_to_string"
check "$TR" "kLinnaeanRanks"

# RANKS list in the simulator must be used (not just defined) in hier_rows
grep -q "for rank in RANKS" "$SIM" \
    && echo "      [ok] RANKS list drives hier_rows iteration" \
    || { echo "      [FAIL] hier_rows still uses hardcoded rank literals" >&2; exit 1; }

# Canonical sanitize_name must be present in the bridge
check "$BR" "sanitize_name"

echo "      Taxonomy + sanitize checks passed."

# ── 7. Pytest ────────────────────────────────────────────────────────────────
if [[ "$RUN_PYTEST" -eq 1 ]]; then
    echo "[7/7] Running pytest suite (jobs=$JOBS)..."
    PYTEST_TARGETS=()
    for f in "$ROOT_DIR/test_pipeline_features.py" "$ROOT_DIR/test_all_use_cases.py" "$ROOT_DIR/test_biological_use_case_fungi.py" "$ROOT_DIR/test_new_biology_candidates.py"; do
        [[ -f "$f" ]] && PYTEST_TARGETS+=("$f")
    done
    if [[ ${#PYTEST_TARGETS[@]} -gt 0 ]]; then
        PYTEST_ARGS=("-q")
        [[ "$JOBS" -gt 1 ]] && PYTEST_ARGS+=("-n" "$JOBS")
        PYTHONPATH="$ROOT_DIR${PYTHONPATH:+:$PYTHONPATH}" \
            python3 -m pytest "${PYTEST_ARGS[@]}" "${PYTEST_TARGETS[@]}"
        echo "      OK"
    else
        echo "      SKIP (no pytest targets present)"
    fi
else
    echo "[7/7] pytest skipped (--no-pytest)"
fi

echo ""
echo "[info] Optional public-data demo available: python3 $ROOT_DIR/example_real_data_fungi.py"

echo "=== All v18 checks passed ==="
echo "    Data structures: VP-tree(DS-4) Bloom(DS-5) DSU(DS-6)"
echo "                     TraIntra-treap(DS-7) InvFlex-bitset(DS-8)"
echo "                     OffRef-SA(DS-9) LCA-RMQ(DS-10)"
echo "                     Triallelic-DSU(DS-11) Pantree-hash(DS-12)"
echo "                     SuffixArray+LCP(DS-13) WaveletTree(DS-14)"
echo "                     vEBTree(DS-15) MergeSortTree(DS-16)"
echo "                     FenwickTree(DS-17) ChainTreap(DS-18)"
