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

# Keep million-real benchmark memory bounded by default. Large benchmark ref
# lists now run hierarchical-only unless MYCOSV_FORCE_FLAT_REF_FALLBACK=1, but
# these caps still protect explicit flat-fallback debug runs.
export MILLION_REAL_SINGLE_REF_CACHE_MB="${MILLION_REAL_SINGLE_REF_CACHE_MB:-4096}"
export MILLION_REAL_MAX_REF_MEMORY_MB="${MILLION_REAL_MAX_REF_MEMORY_MB:-4096}"

exec bash /mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/AMF/scale/run_all_experiments.sh --million-real "$@"
