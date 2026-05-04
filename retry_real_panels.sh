#!/usr/bin/env bash
# retry_real_panels.sh
#
# Re-run the panels of a real-data benchmark that died from cgroup OOM or
# timeout, reusing the existing prepared/ directory of an earlier run.
# Run under srun --mem=32G to escape the 12 GiB user-slice cap.
#
# Usage:
#   sbatch -p <queue> --mem=32G --cpus-per-task=32 --time=12:00:00 \
#       retry_real_panels.sh <RUN_TIMESTAMP_DIR> [panel ...]
#
# Or interactive:
#   srun --mem=32G --cpus-per-task=32 --time=12:00:00 \
#       bash retry_real_panels.sh experiments/real_data/20260428_182921 \
#            te_rich_pathogen two_speed_pathogen amf_large cross_phylum_hgt
#
# If no panels are given, all panels under RUN_DIR are retried. Each panel/mode
# writes to benchmark_<mode>_retry/ so the original outputs are preserved.

set -u
set -o pipefail

WORK_DIR="/mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/AMF/scale"
RUN_DIR="${1:?usage: $0 <RUN_TIMESTAMP_DIR> [panel ...]}"
shift || true
RUN_DIR="$(readlink -f "$RUN_DIR")"

THREADS="${THREADS:-32}"
# Raised vs run_all_experiments.sh: under srun --mem=32G we can pack more
# clade genomes into the routing index, which improves both the
# routing-accuracy and the chance mycosv finds the right reference per region.
REAL_MAX_CLADE_GENOMES="${REAL_MAX_CLADE_GENOMES:-16}"
# 4-hour per-tool ceiling so the slow-but-correct cactus / pggb runs on
# 100-500 Mb fungal genomes don't get killed mid-run like in 20260428.
TOOL_TIMEOUT_S="${TOOL_TIMEOUT_S:-14400}"
export FUNGI_BENCH_TOOL_TIMEOUT="${TOOL_TIMEOUT_S}"

if [[ ! -d "${RUN_DIR}" ]]; then
  echo "[retry] RUN_DIR does not exist: ${RUN_DIR}" >&2
  exit 1
fi

if [[ $# -eq 0 ]]; then
  PANELS=()
  for d in "${RUN_DIR}"/*/; do
    [[ -d "${d}prepared" ]] && PANELS+=("$(basename "${d%/}")")
  done
else
  PANELS=("$@")
fi

if [[ ${#PANELS[@]} -eq 0 ]]; then
  # Bootstrap missing panels (e.g. two_speed_pathogen) from a fresh prepare.
  PANELS=(compact_yeast amf_large cross_phylum_hgt te_rich_pathogen two_speed_pathogen)
fi

# Shared across retries and fresh prepares so large FASTA/GFF/FASTQ downloads
# and metadata caches are not repeated for each timestamped run.
DATA_CACHE_DIR="${DATA_CACHE_DIR:-${WORK_DIR}/data_cache}"

cd "${WORK_DIR}"

retry_one_panel() {
  local panel="$1"
  local panel_dir="${RUN_DIR}/${panel}"
  mkdir -p "${panel_dir}"

  # Ensure prepared/ exists. If not (panel never ran prepare), do a fresh
  # prepare; otherwise reuse the existing one.
  if [[ ! -f "${panel_dir}/prepared/query_manifest.tsv" ]]; then
    echo "[retry] ${panel}: no prepared dir, running prepare..." >&2
    REAL_MAX_REF_DOWNLOADS="${REAL_MAX_REF_DOWNLOADS:-6}"
    REAL_MAX_ASMS_PER_SPECIES="${REAL_MAX_ASMS_PER_SPECIES:-3}"
    REAL_QUERIES_PER_SPECIES="${REAL_QUERIES_PER_SPECIES:-3}"
    REAL_MAX_QUERY_DOWNLOADS="${REAL_MAX_QUERY_DOWNLOADS:-6}"
    python3 run_real_fungal_benchmark.py prepare \
        --out-dir "${panel_dir}/prepared" \
        --panel "${panel}" \
        --source ncbi-genbank \
        --max-assemblies-per-species "${REAL_MAX_ASMS_PER_SPECIES}" \
        --querys-per-species "${REAL_QUERIES_PER_SPECIES}" \
        --max-ref-downloads "${REAL_MAX_REF_DOWNLOADS}" \
        --max-query-downloads "${REAL_MAX_QUERY_DOWNLOADS}" \
        --query-mode mixed \
        --read-accessions-per-species 2 \
        --allow-no-queries \
        --data-cache-dir "${DATA_CACHE_DIR}" \
        2>&1 | tee "${panel_dir}/prepare.retry.log"
  fi

  if [[ ! -f "${panel_dir}/prepared/query_manifest.tsv" ]]; then
    echo "[retry] ${panel}: prepare did not produce a query manifest, skipping" >&2
    return 1
  fi

  # Validate query inputs upfront — the 2026-04-28 run had a 3.9 KB Rhizophagus
  # FASTQ that gave mycosv "no sequences read" with no clean failure.
  awk -F'\t' 'NR>1 {print $3}' "${panel_dir}/prepared/query_manifest.tsv" | while read -r path; do
    if [[ -n "${path}" && -f "${path}" ]]; then
      bytes=$(stat -c %s "${path}")
      if [[ ${bytes} -lt 100000 ]]; then
        echo "[retry][warn] ${panel}: query file looks truncated (${bytes} bytes): ${path}" >&2
      fi
    fi
  done

  for mode in assembly short-reads long-reads; do
    local out="${panel_dir}/benchmark_${mode}_retry"
    local read_validation_min_support
    if [[ "${mode}" == "assembly" ]]; then
      read_validation_min_support="${REAL_READ_VALIDATION_MIN_SUPPORT_ASSEMBLY:-1}"
    else
      read_validation_min_support="${REAL_READ_VALIDATION_MIN_SUPPORT_READS:-3}"
    fi
    mkdir -p "${out}"
    echo "[retry] ${panel} / ${mode} -> ${out}" >&2
    python3 run_real_fungal_benchmark.py benchmark \
        --prepared-dir "${panel_dir}/prepared" \
        --mode "${mode}" \
        --out-dir "${out}" \
        --threads "${THREADS}" \
        --max-clade-genomes "${REAL_MAX_CLADE_GENOMES}" \
        --read-validation-min-support "${read_validation_min_support}" \
        --run-all-comparators \
        2>&1 | tee "${panel_dir}/benchmark_${mode}_retry.log"
  done
}

for panel in "${PANELS[@]}"; do
  retry_one_panel "${panel}" || echo "[retry] panel ${panel} failed" >&2
done

# Post-pass: recompute exact_benchmark_summary.tsv for every panel/mode using
# the new reference_any_clade + consensus_2of_n columns.
python3 recompute_benchmark_metrics.py "${RUN_DIR}"/*/

echo "[retry] done. New outputs under benchmark_<mode>_retry/ in each panel."
