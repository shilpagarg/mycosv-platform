#!/usr/bin/env bash
#SBATCH --job-name=mycosv-debug64
#SBATCH -p multicore
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=0:30:00
#SBATCH --output=slurm-debug-%j.out

set -u
set -o pipefail

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "${PROJECT_DIR}"

PREPARED_DIR="${PREPARED_DIR:-experiments/million_real/20260518_105922}"
OUT_DIR="${OUT_DIR:-${PREPARED_DIR}/benchmark_assembly_smoke30}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
MAX_BENCHMARK_QUERIES="${MAX_BENCHMARK_QUERIES:-1}"
BENCHMARK_REF_CAP="${BENCHMARK_REF_CAP:-8}"
MAX_CALLS_PER_CONTIG="${MAX_CALLS_PER_CONTIG:-1000}"
MIN_BLOCK_SCORE="${MIN_BLOCK_SCORE:-4.0}"
TOL_MIN_CHAIN_ANCHORS="${TOL_MIN_CHAIN_ANCHORS:-2}"
MYCOSV_TOOL_TIMEOUT="${MYCOSV_TOOL_TIMEOUT:-1500}"
# Bug 4 fix: cap unique ref contigs the flat fallback iterates. 8 fungal ref
# FASTAs can carry 6 000+ contigs; without this the fallback never finishes
# inside MYCOSV_TOOL_TIMEOUT on a smoke run.
MAX_FLAT_REF_CONTIGS="${MAX_FLAT_REF_CONTIGS:-256}"
# Bug 2 fix: skip the flat-MEM-chain fallback when hierarchical already
# produced this many calls for a given query. 0 disables the gate.
SKIP_FLAT_IF_HIER_CALLS="${SKIP_FLAT_IF_HIER_CALLS:-5}"

echo "MycoSV bounded assembly debug"
echo "start: $(date)"
echo "cwd: $(pwd)"
echo "prepared_dir: ${PREPARED_DIR}"
echo "out_dir: ${OUT_DIR}"
echo "threads: ${THREADS}"
echo "max_benchmark_queries: ${MAX_BENCHMARK_QUERIES}"
echo "benchmark_ref_cap: ${BENCHMARK_REF_CAP}"
echo "max_calls_per_contig: ${MAX_CALLS_PER_CONTIG}"
echo "max_flat_ref_contigs: ${MAX_FLAT_REF_CONTIGS}"
echo "skip_flat_if_hier_calls: ${SKIP_FLAT_IF_HIER_CALLS}"
echo "mycosv_tool_timeout: ${MYCOSV_TOOL_TIMEOUT}"

if [[ ! -f "${PREPARED_DIR}/query_manifest.tsv" ]]; then
  echo "[error] missing ${PREPARED_DIR}/query_manifest.tsv" >&2
  exit 2
fi
if [[ ! -d "${PREPARED_DIR}/index" || ! -d "${PREPARED_DIR}/registry" ]]; then
  echo "[error] missing prepared index/registry under ${PREPARED_DIR}" >&2
  exit 2
fi

# For this diagnostic we intentionally keep one query and a tiny benchmark ref
# set. This tells us quickly whether MycoSV can emit assembly SVs before
# spending time on comparators/read validation.
export MYCOSV_FORCE_FLAT_REF_FALLBACK=1
export MYCOSV_TOOL_TIMEOUT

python3 -u run_real_fungal_benchmark.py benchmark \
  --prepared-dir "${PREPARED_DIR}" \
  --out-dir "${OUT_DIR}" \
  --mode assembly \
  --threads "${THREADS}" \
  --max-clade-genomes 32 \
  --max-benchmark-queries "${MAX_BENCHMARK_QUERIES}" \
  --mycosv-only \
  --no-validate-with-reads \
  --reuse-index-dir "${PREPARED_DIR}/index" \
  --reuse-registry-dir "${PREPARED_DIR}/registry" \
  --benchmark-ref-cap "${BENCHMARK_REF_CAP}" \
  --mycosv-arg=--max-calls-per-contig \
  --mycosv-arg="${MAX_CALLS_PER_CONTIG}" \
  --mycosv-arg=--min-block-score \
  --mycosv-arg="${MIN_BLOCK_SCORE}" \
  --mycosv-arg=--tol-min-chain-anchors \
  --mycosv-arg="${TOL_MIN_CHAIN_ANCHORS}" \
  --mycosv-arg=--max-ref-memory-mb \
  --mycosv-arg=8192 \
  --mycosv-arg=--max-flat-ref-contigs \
  --mycosv-arg="${MAX_FLAT_REF_CONTIGS}" \
  --mycosv-arg=--skip-flat-if-hier-calls \
  --mycosv-arg="${SKIP_FLAT_IF_HIER_CALLS}" \
  --mycosv-arg=--no-gfa

echo
echo "Debug counts"
if [[ -f "${OUT_DIR}/mycosv/calls.vcf" ]]; then
  printf "raw_vcf_calls\t"
  grep -vc '^#' "${OUT_DIR}/mycosv/calls.vcf" || true
else
  echo "raw_vcf_calls	MISSING"
fi

for f in \
  "${OUT_DIR}/pangenome_call_layers.tsv" \
  "${OUT_DIR}/sv_volume_audit.tsv" \
  "${OUT_DIR}/exact_benchmark_summary.tsv" \
  "${OUT_DIR}/novel_mycosv_calls.tsv" \
  "${OUT_DIR}/biology_findings.tsv"
do
  if [[ -f "${f}" ]]; then
    echo
    echo "== ${f} =="
    sed -n '1,12p' "${f}"
  else
    echo
    echo "== ${f} MISSING =="
  fi
done

echo
echo "end: $(date)"
