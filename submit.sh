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
exec bash /mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/AMF/scale/run_all_experiments.sh "$@"
