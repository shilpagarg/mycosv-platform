#!/usr/bin/env bash
#SBATCH --job-name=mycosv-scale
#SBATCH -p multicore
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00

set -u
set -o pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"
export THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"

# Force a clean rebuild of the MycoSV binary before the matrix runs.
# Header-only edits (e.g. the SvTypeFromChain DUP per-pair detector in
# layer1_clade_graph.hpp) only take effect if the binary is recompiled.
# compile_binary_if_needed inside run_real_fungal_benchmark.py uses mtime
# comparison, but a stale binary whose mtime was touched after the headers
# would silently keep the old object code. Removing the binary here
# guarantees the next invocation rebuilds from current source.
MYCOSV_BIN="$(dirname "${BASH_SOURCE[0]}")/fungi_graphsv_tol_bin"
if [[ "${MYCOSV_FORCE_REBUILD:-1}" == "1" && -f "${MYCOSV_BIN}" ]]; then
  echo "[submit] forcing rebuild: removing ${MYCOSV_BIN}"
  rm -f "${MYCOSV_BIN}"
fi

# Million-real prepare/index build is the stage that has historically lost
# read_validated_truth.tsv rows to mid-loop OOM kills. The benchmark step
# now flushes the TSV after every query, but the underlying memory pressure
# is still there. Bumping the per-thread SA cache budget keeps the binary
# from spilling to std::bad_alloc under the SLURM cgroup. Override on
# tighter nodes or when running with --cpus-per-task < 16.
export MILLION_REAL_SINGLE_REF_CACHE_MB="${MILLION_REAL_SINGLE_REF_CACHE_MB:-4096}"
export MILLION_REAL_MAX_REF_MEMORY_MB="${MILLION_REAL_MAX_REF_MEMORY_MB:-1024}"

exec bash /mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/AMF/scale/run_all_experiments.sh "$@"
