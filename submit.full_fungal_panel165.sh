#!/usr/bin/env bash
# Stage-3 launcher: comprehensive 165-species fungal panel.
#
# Thin wrapper around submit.full_fungal_assembly.sh that bumps QUERY_COUNT
# from 15 to 165 and broadens QUERY_GROUPS to cover every genus in the
# AMF / EMF / endophyte / filamentous / yeast wish-list. Launches the same
# 3-phase chain (bootstrap -> array -> combine) into a fresh OUT_ROOT so the
# 15-query results stay intact for comparison.
#
# NOTE: QUERY_GROUPS contains commas, which sbatch's --export= flag treats
# as a variable separator. We dodge that by `export`-ing everything in this
# shell and passing only `--export=ALL` to sbatch.
#
# Usage (after the in-flight 15-query run completes):
#   ./submit.full_fungal_panel165.sh                 # submits the chain
#   ./submit.full_fungal_panel165.sh --dry-run       # prints commands only
#
# Wall-time estimate: ~1.5 days with default array concurrency.

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${PROJECT_DIR}"

DRY_RUN=0
if [[ "${1:-}" == "--dry-run" ]]; then
  DRY_RUN=1
fi

# Fresh timestamped run dir so the in-flight 15-query results stay immutable.
TIMESTAMP="${TIMESTAMP:-panel165_$(date +%Y%m%d_%H%M%S)}"

# Export every var that submit.full_fungal_assembly.sh + run_real_fungal_benchmark.py
# rely on. With --export=ALL the SLURM job inherits these intact, commas and all.
export OUT_ROOT="experiments/million_real/full_fungal_assembly_${TIMESTAMP}"
export QUERY_COUNT=165
export QUERY_GROUPS="Acaulospora,Acremonium,Agaricus,Albahypha,Alternaria,Amanita,Ambispora,Archaeospora,Aspergillus,Austroboletus,Beauveria,Blaszkowskia,Boletus,Butyriboletus,Candida,Cantharellus,Cenococcum,Cetraspora,Chroogomphus,Cladosporium,Claroideoglomus,Clavulina,Colletotrichum,Cortinarius,Dentiscutata,Diversispora,Dominikia,Entrophospora,Eremothecium,Funneliformis,Fusarium,Geastrum,Geopora,Gigaspora,Glomus,Gomphidius,Gyroporus,Hebeloma,Helvella,Hydnellum,Hydnum,Hygrophorus,Inocybe,Laccaria,Lactarius,Lanmaoa,Leccinum,Lentinula,Meliniomyces,Monascus,mycorrhiza,Neurospora,Oidiodendron,Pacispora,Paecilomyces,Paraglomus,Paxillus,Penicillium,Pervetustus,Piriformospora,Pisolithus,Pleurotus,Pseudosperma,Racocetra,Redeckera,Rhizoglomus,Rhizophagus,Rhizopogon,Rhizopus,Rhizoscyphus,Russula,Saccharomyces,Sclerocystis,Scleroderma,Scutellospora,Septoglomus,Serendipita,Sieverdingia,Simiglomus,Sparassis,Sphaerosporella,Suillellus,Suillus,Talaromyces,Thelephora,Trichoderma,Tricholoma,Trichophaea,Tuber,Tulasnella,Tylopilus,Ustilago,Xerocomus,Yarrowia"
export MIN_ASSEMBLY_LEVEL=contig
export MYCOSV_BIOLOGY_TOP_N=5000

ARRAY_RANGE="0-$((QUERY_COUNT - 1))"
# Throttle concurrent array tasks so we don't starve the cluster.
ARRAY_CONCURRENCY="${ARRAY_CONCURRENCY:-30}"

echo "=== panel-165 launch ==="
echo "  out_root:         ${OUT_ROOT}"
echo "  query_count:      ${QUERY_COUNT}"
echo "  array_range:      ${ARRAY_RANGE}  (concurrency=${ARRAY_CONCURRENCY})"
echo "  query_groups:     ${#QUERY_GROUPS} chars, 94 keywords"
echo "  min_assembly_lvl: ${MIN_ASSEMBLY_LEVEL} (so conspecific contig assemblies are eligible)"
echo "  biology_top_n:    ${MYCOSV_BIOLOGY_TOP_N}"
echo "  dry_run:          ${DRY_RUN}"
echo

# Wrapper that prints (dry) or submits and returns just the job ID.
submit() {
  if (( DRY_RUN )); then
    echo "  [dry] sbatch $*" >&2
    echo "9999999"  # fake id so the dependency strings still look sane
  else
    sbatch --parsable "$@"
  fi
}

# Phase A — bootstrap: runs prepare-million-real (downloads + indexes 165
# queries, writes prepared/query_manifest.tsv + routing index + read-
# validation manifest), then BOOTSTRAP_ONLY=1 exits before the benchmark
# launches. The submit.full_fungal_assembly.sh script was reordered so the
# BOOTSTRAP_ONLY exit happens AFTER prepare, not before — that fixes the
# earlier "bootstrap exited in 4s without preparing anything" failure.
BOOTSTRAP_OVERRIDES="BOOTSTRAP_ONLY=1,SKIP_PREPARE=0"
BS=$(submit \
  --job-name=panel165-bootstrap \
  -p multicore -c 8 --mem=32G --time=24:00:00 \
  --output=slurm-panel165-bootstrap-%j.out \
  --export=ALL,${BOOTSTRAP_OVERRIDES} \
  submit.full_fungal_assembly.sh)
echo "bootstrap_jobid=${BS}"

# Phase B — 165-way array, throttled to ARRAY_CONCURRENCY parallel tasks.
# SKIP_PREPARE=1 reuses prepared/ from bootstrap; FORCE_RERUN_SHARDS=1 starts
# each shard from scratch (matters when bug fixes change call semantics).
ARRAY_OVERRIDES="FORCE_RERUN_SHARDS=1,SKIP_PREPARE=1"
AR=$(submit \
  --dependency=afterok:${BS} \
  --array=${ARRAY_RANGE}%${ARRAY_CONCURRENCY} \
  --export=ALL,${ARRAY_OVERRIDES} \
  submit.full_fungal_assembly.sh)
echo "array_jobid=${AR}"

# Phase C — combine all 165 shards into combined/ + master reports.
COMBINE_OVERRIDES="FULL_ASSEMBLY_SHARDS=1,RESUME_SHARDS=1,SKIP_PREPARE=1"
CB=$(submit \
  --dependency=afterany:${AR} \
  --job-name=panel165-combine \
  --export=ALL,${COMBINE_OVERRIDES} \
  submit.full_fungal_assembly.sh)
echo "combine_jobid=${CB}"

echo
if (( DRY_RUN )); then
  echo "[dry-run] no jobs submitted; rerun without --dry-run to launch"
else
  echo "[submitted] chain: ${BS} -> ${AR} -> ${CB}"
  echo "[monitor]   squeue -u \$USER | grep panel165"
  echo "[monitor]   tail -f slurm-panel165-bootstrap-${BS}.out"
fi
