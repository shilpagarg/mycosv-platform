#!/usr/bin/env bash
#SBATCH --job-name=mycosv-dbg-hier
#SBATCH -p multicore
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=1:00:00
#SBATCH --output=slurm-dbg-hier-%j.out

# Single-query diagnostic: emit per-contig breadcrumbs through hierarchical
# Path A/B/C so we can see why the multisample VCF is OFF_REF-only.

set -u
set -o pipefail

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "${PROJECT_DIR}"

COMPARATOR_ENV="/mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/tools/envs/envs/fungi_graph_sv"
if [[ -d "${COMPARATOR_ENV}/bin" ]]; then
  export PATH="${COMPARATOR_ENV}/bin:${PATH}"
fi

PREPARED_DIR="experiments/million_real/20260518_105922"
OUT_DIR="${PREPARED_DIR}/benchmark_dbg_hier"
THREADS=8

export MYCOSV_FORCE_FLAT_REF_FALLBACK=1
export MYCOSV_TOOL_TIMEOUT=3000
export MYCOSV_DEBUG_HIER=1

echo "start: $(date)"
echo "binary: $(realpath ./fungi_graphsv_tol_bin.dbg)"

python3 -u run_real_fungal_benchmark.py benchmark \
  --prepared-dir "${PREPARED_DIR}" \
  --out-dir "${OUT_DIR}" \
  --mode assembly \
  --threads "${THREADS}" \
  --max-clade-genomes 32 \
  --max-benchmark-queries 1 \
  --mycosv-only \
  --no-validate-with-reads \
  --reuse-index-dir "${PREPARED_DIR}/index" \
  --reuse-registry-dir "${PREPARED_DIR}/registry" \
  --benchmark-ref-cap 8 \
  --binary-path "$(realpath ./fungi_graphsv_tol_bin.dbg)" \
  --mycosv-arg=--max-calls-per-contig --mycosv-arg=1000 \
  --mycosv-arg=--min-block-score    --mycosv-arg=4.0 \
  --mycosv-arg=--tol-min-chain-anchors --mycosv-arg=2 \
  --mycosv-arg=--max-ref-memory-mb  --mycosv-arg=8192 \
  --mycosv-arg=--max-flat-ref-contigs --mycosv-arg=256 \
  --mycosv-arg=--skip-flat-if-hier-calls --mycosv-arg=5 \
  --mycosv-arg=--no-gfa

echo "end: $(date)"
echo
echo "=== HIER-DBG lines ==="
grep -E "\[hier-dbg\]" "${OUT_DIR}/mycosv/calls.stderr.log" | head -80
echo
echo "=== SVTYPE summary ==="
if [[ -f "${OUT_DIR}/mycosv/calls.vcf" ]]; then
  grep -v "^#" "${OUT_DIR}/mycosv/calls.vcf" | grep -oE "SVTYPE=[A-Z_]+" | sort | uniq -c
fi
