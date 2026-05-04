#!/usr/bin/env bash
#SBATCH -p multicore
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00

set -u
set -o pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"
exec bash /mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/AMF/scale/run_all_experiments.sh "$@"
