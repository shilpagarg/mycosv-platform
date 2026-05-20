#!/usr/bin/env bash
#SBATCH --job-name=mycosv-fusarium
#SBATCH -p multicore
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --output=slurm-fusarium-%A-%a.out
#SBATCH --array=0-2

# Replaces the ad-hoc benchmark_mycosv_{asm,lr,lr2k}_fusarium runs
# (--mycosv-only, 30 min walltime, MYCOSV_TOOL_TIMEOUT=900) that produced:
#   - 0 DUP, 0 TRA on the F. falciforme vs F. oxysporum asm comparison
#   - 3,584 phantom NOVEL_WEAK 250 bp tile stubs filling kMaxOffRefWindows
#   - header-only read_validated_truth.tsv (Python killed before validation)
#   - exact_benchmark_summary.tsv rows pinned at no_comparator / NaN F1
#
# What's different here:
#   - 24 h SLURM wall (was 30 min) so the post-MycoSV validation + biology
#     pipeline can complete after the binary finishes.
#   - MYCOSV_TOOL_TIMEOUT=10800 s (3 h) so the long-read 23.4× run does not
#     get SIGKILL'd mid-call the way the 900 s budget did.
#   - Real comparators enabled per mode (was --mycosv-only):
#       asm        : --run-minigraph --run-svim-asm --run-anchorwave
#       long-reads : --run-svim --run-sniffles --run-cutesv
#   - Full long-read coverage (no 2k subsample) for the lr task — the 2k
#     subsample collapsed coverage to 0.06× and just hid the timeout.

set -u
set -o pipefail

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "${PROJECT_DIR}"

# Comparator binaries (minigraph / svim-asm / sniffles / svim / cuteSV /
# samtools / minimap2) live in the conda env the wrapper expects on PATH.
COMPARATOR_ENV="/mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/tools/envs/envs/fungi_graph_sv"
if [[ -d "${COMPARATOR_ENV}/bin" ]]; then
  export PATH="${COMPARATOR_ENV}/bin:${PATH}"
fi

PREPARED_DIR="${PREPARED_DIR:-experiments/million_real/20260518_105922}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"

# Per-array-task config: 0=asm, 1=long-reads full, 2=long-reads quick.
MODES=("assembly" "long-reads" "long-reads")
OUTNAMES=("benchmark_fusarium_asm" "benchmark_fusarium_lr" "benchmark_fusarium_lr_quick")
LR_READ_CAPS=(0 0 20000)   # 0 = no cap (full coverage)

idx="${SLURM_ARRAY_TASK_ID:-0}"
if (( idx < 0 || idx > 2 )); then
  echo "[error] SLURM_ARRAY_TASK_ID=${idx} out of range 0..2" >&2
  exit 2
fi
MODE="${MODES[$idx]}"
OUT_NAME="${OUTNAMES[$idx]}"
LR_READ_CAP="${LR_READ_CAPS[$idx]}"
OUT_DIR="${PREPARED_DIR}/${OUT_NAME}"

# Per-mode comparator flag set. Mirrors run_all_experiments.sh:485-501 but
# tied to the Fusarium one-genus panel instead of the full million_real
# matrix.
COMPARATOR_FLAGS=()
case "${MODE}" in
  assembly)
    COMPARATOR_FLAGS+=(--run-minigraph --run-svim-asm)
    if [[ "${FUSARIUM_RUN_ANCHORWAVE:-1}" == "1" ]]; then
      COMPARATOR_FLAGS+=(--run-anchorwave)
    fi
    READVAL_SUPPORT="${FUSARIUM_READVAL_SUPPORT_ASM:-1}"
    ;;
  long-reads)
    COMPARATOR_FLAGS+=(--run-svim --run-sniffles --run-cutesv)
    READVAL_SUPPORT="${FUSARIUM_READVAL_SUPPORT_LR:-2}"
    ;;
  *)
    echo "[error] unsupported MODE=${MODE}" >&2
    exit 2
    ;;
esac

# Long-read budget: 0 = default cap (200 000 reads / ~23× cov on the
# Fusarium SRR21444561 sample), >0 = explicit subsample for the quick task.
LR_READ_CAP_FLAGS=()
if [[ "${MODE}" == "long-reads" && "${LR_READ_CAP}" -gt 0 ]]; then
  LR_READ_CAP_FLAGS+=(--max-comparator-long-reads "${LR_READ_CAP}")
fi

# Total per-call timeout for the C++ MycoSV binary. 900 s killed the long-
# read full-coverage run mid-call last time; 10800 s leaves headroom for
# all 14 F. falciforme chromosomes plus three sibling-clade neighbors.
export MYCOSV_TOOL_TIMEOUT="${MYCOSV_TOOL_TIMEOUT:-10800}"
# Hierarchical-only is fine for these comparisons; the flat-MEM-chain
# fallback is what blew the 900 s budget previously. Leave the env var
# untouched unless explicitly debugging that path.
unset MYCOSV_FORCE_FLAT_REF_FALLBACK || true

# Memory caps mirror submit.sh defaults (CLAUDE memory note: cgroup is
# 12 G on interactive shells; SLURM job sees the full 128 G allocation).
export MILLION_REAL_SINGLE_REF_CACHE_MB="${MILLION_REAL_SINGLE_REF_CACHE_MB:-4096}"
export MILLION_REAL_MAX_REF_MEMORY_MB="${MILLION_REAL_MAX_REF_MEMORY_MB:-8192}"

# Force a clean rebuild before the array starts — header-only fixes
# (the OFF_REF NOVEL_WEAK tile filter in fungi_tol_bridge.hpp, the
# benchmark-ref self-exclusion in run_real_fungal_benchmark.py) only
# take effect with a fresh binary. Re-builds are guarded so only the
# first array task does it.
MYCOSV_BIN="${PROJECT_DIR}/fungi_graphsv_tol_bin"
REBUILD_LOCK="${PROJECT_DIR}/.fusarium_rebuild.lock"
if [[ "${idx}" == "0" ]] && [[ "${FUSARIUM_FORCE_REBUILD:-1}" == "1" ]]; then
  rm -f "${MYCOSV_BIN}" "${REBUILD_LOCK}"
elif [[ -f "${REBUILD_LOCK}" ]]; then
  : # other tasks wait until task 0 has rebuilt
fi

if [[ ! -f "${PREPARED_DIR}/query_manifest.tsv" ]]; then
  echo "[error] missing ${PREPARED_DIR}/query_manifest.tsv" >&2
  exit 2
fi
if [[ ! -d "${PREPARED_DIR}/index" || ! -d "${PREPARED_DIR}/registry" ]]; then
  echo "[error] missing prepared index/registry under ${PREPARED_DIR}" >&2
  exit 2
fi

echo "=== Fusarium benchmark task ${idx} ==="
echo "start:      $(date)"
echo "mode:       ${MODE}"
echo "out_dir:    ${OUT_DIR}"
echo "threads:    ${THREADS}"
echo "comparators:${COMPARATOR_FLAGS[*]}"
echo "lr_cap:     ${LR_READ_CAP}"
echo "tool_timeout:${MYCOSV_TOOL_TIMEOUT}s"
echo "minigraph:  $(command -v minigraph  || echo MISSING)"
echo "svim-asm:   $(command -v svim-asm   || echo MISSING)"
echo "anchorwave: $(command -v anchorwave || echo MISSING)"
echo "svim:       $(command -v svim       || echo MISSING)"
echo "sniffles:   $(command -v sniffles   || echo MISSING)"
echo "cuteSV:     $(command -v cuteSV     || command -v cutesv || echo MISSING)"
echo "samtools:   $(command -v samtools   || echo MISSING)"
echo "minimap2:   $(command -v minimap2   || echo MISSING)"

python3 -u run_real_fungal_benchmark.py benchmark \
  --prepared-dir "${PREPARED_DIR}" \
  --out-dir "${OUT_DIR}" \
  --mode "${MODE}" \
  --threads "${THREADS}" \
  --max-clade-genomes 32 \
  --benchmark-query-genera Fusarium \
  --max-benchmark-queries 1 \
  --reuse-index-dir "${PREPARED_DIR}/index" \
  --reuse-registry-dir "${PREPARED_DIR}/registry" \
  --benchmark-ref-cap 8 \
  --read-validation-min-support "${READVAL_SUPPORT}" \
  "${COMPARATOR_FLAGS[@]}" \
  "${LR_READ_CAP_FLAGS[@]}" \
  --mycosv-arg=--max-calls-per-contig --mycosv-arg=2000 \
  --mycosv-arg=--min-block-score --mycosv-arg=4.0 \
  --mycosv-arg=--tol-min-chain-anchors --mycosv-arg=2 \
  --mycosv-arg=--max-ref-memory-mb --mycosv-arg=8192 \
  --mycosv-arg=--max-flat-ref-contigs --mycosv-arg=256 \
  --mycosv-arg=--skip-flat-if-hier-calls --mycosv-arg=5 \
  --mycosv-arg=--no-gfa

RC=$?

echo
echo "=== Result for ${OUT_NAME} ==="
for f in mycosv/calls.vcf exact_benchmark_summary.tsv read_validated_truth.tsv \
         biology_findings.tsv pangenome_call_layers.tsv sv_volume_audit.tsv \
         benchmark_summary.json; do
  p="${OUT_DIR}/${f}"
  if [[ -f "${p}" ]]; then
    n=$(grep -vc '^#' "${p}" 2>/dev/null || echo "?")
    printf "  %-32s %s lines\n" "${f}" "${n}"
  else
    printf "  %-32s MISSING\n" "${f}"
  fi
done

VCF="${OUT_DIR}/mycosv/calls.vcf"
if [[ -f "${VCF}" ]]; then
  echo
  echo "  SVTYPE mix:"
  grep -v '^#' "${VCF}" | grep -oE 'SVTYPE=[A-Z_]+' | sort | uniq -c \
    | awk '{printf "    %-18s %s\n", $2, $1}'
fi

echo
echo "end: $(date)"
exit "${RC}"
