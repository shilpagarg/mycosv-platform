#!/usr/bin/env bash
#SBATCH -p serial
#SBATCH --job-name=starfish
#SBATCH --cpus-per-task=1
#SBATCH --mem=5300M
#SBATCH --time=12:00:00
#SBATCH --array=1-5
#SBATCH --output=logs/slurm/slurm-starfish-%A_%a.out
set -uo pipefail
R="${MYCOSV_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
cd $R; export PATH="$R/tools_src/starfish/bin:$R/tools_src/elembench/bin:$PATH"
DB="$R/tools_src/starfish/db"; GEN="$R/experiments/biology_benchmark/genomes"
read -r tag stem < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $R/experiments/biology_benchmark/genomes.tsv)
wd="$R/experiments/biology_benchmark/starfish/$tag"
mkdir -p "$wd/geneFinder" "$wd/elementFinder" "$wd/blastdb"; cd "$wd"
echo "===== STARFISH $tag ($stem) $(date) ====="
printf "%s\t%s\n" "$tag" "$GEN/$stem.fna" > ome2assembly.txt
makeblastdb -in "$GEN/$stem.fna" -out blastdb/asm -parse_seqids -dbtype nucl >/dev/null 2>&1 || echo "[warn] makeblastdb"
starfish annotate -T 1 -x $tag -a ome2assembly.txt -p "$DB/YRsuperfams.p1-512.hmm" -P "$DB/YRsuperfamRefs.faa" -i tyr -o geneFinder/ || echo "[warn] annotate rc=$?"
gff=$(ls geneFinder/${tag}*.filt_intersect.gff geneFinder/${tag}*.filt.gff 2>/dev/null | head -1)
ids=$(ls geneFinder/${tag}*.filt_intersect.ids geneFinder/${tag}*.filt.ids 2>/dev/null | head -1)
if [[ -s "$gff" && -s "$ids" ]]; then
  printf "%s\t%s\n" "$tag" "$PWD/$gff" > ome2gff.txt
  starfish sketch -m 10000 -q "$ids" -g ome2gff.txt -i s -x $tag -o geneFinder/ || echo "[warn] sketch rc=$?"
  grep -P '\ttyr\t' geneFinder/${tag}.bed > geneFinder/${tag}.tyr.bed 2>/dev/null || echo "[warn] no tyr.bed"
  if [[ -s geneFinder/${tag}.tyr.bed ]]; then
    starfish insert -T 1 -a ome2assembly.txt -d blastdb/asm -b geneFinder/${tag}.tyr.bed -i tyr -x $tag -o elementFinder/ || echo "[warn] insert rc=$?"
    starfish flank -a ome2assembly.txt -b elementFinder/${tag}.insert.bed -x $tag -o elementFinder/ || echo "[warn] flank rc=$?"
  fi
fi
echo "--- outputs ---"; find geneFinder elementFinder -maxdepth 1 -type f | sort
echo "===== done $(date) ====="
