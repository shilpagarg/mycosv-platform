#!/usr/bin/env bash
#SBATCH -p serial
#SBATCH --job-name=te-sorter
#SBATCH --cpus-per-task=1
#SBATCH --mem=5300M
#SBATCH --time=12:00:00
#SBATCH --array=1-5
#SBATCH --output=logs/slurm/slurm-tesorter-%A_%a.out
set -uo pipefail
R=/mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/AMF/scale
cd $R; export PATH="$R/tools_src/elembench/bin:$PATH"
read -r tag stem < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $R/experiments/biology_benchmark/genomes.tsv)
cd $R/experiments/biology_benchmark; mkdir -p tesorter; cd tesorter
echo "=== TEsorter $tag ($stem) $(date) ==="
[ -s ${stem}.rexdb.cls.tsv ] && { echo "already done"; exit 0; }
TEsorter ../genomes/${stem}.fna -genome -db rexdb -p 1 -pre ${stem}.rexdb
echo "=== done $(date) rc=$? ==="; ls -la ${stem}* 2>/dev/null
