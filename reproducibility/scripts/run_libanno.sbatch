#!/usr/bin/env bash
#SBATCH -p serial
#SBATCH --job-name=libanno
#SBATCH --cpus-per-task=1
#SBATCH --mem=5300M
#SBATCH --time=06:00:00
#SBATCH --array=1-78%72
#SBATCH --output=logs/slurm/slurm-libanno-%A_%a.out
set -uo pipefail
R=/mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/AMF/scale
cd $R; export PATH="$R/tools_src/elembench/bin:$PATH"
R2=$R/experiments/biology_benchmark
read -r tag k chunk lib < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $R2/chunks/blast_tasks.tsv)
if [ "$lib" = "te" ]; then DB=$R2/mycomobilome/telib; OUT=$R2/mycomobilome/${tag}.te.part${k}.tsv; EV=1e-10
else                        DB=$R2/starbase/shipdb;   OUT=$R2/starbase/${tag}.ship.part${k}.tsv; EV=1e-20; fi
echo "=== $tag chunk$k $lib $(date) ==="
[ -s "$OUT" ] && { echo "already done"; exit 0; }
blastn -task dc-megablast -query "$chunk" -db "$DB" \
  -outfmt '6 qseqid qstart qend sseqid pident length evalue bitscore' \
  -evalue $EV -num_threads 1 -max_target_seqs 5 -out "$OUT.tmp" && mv "$OUT.tmp" "$OUT"
echo "rows: $(wc -l < $OUT 2>/dev/null || echo 0) $(date)"
