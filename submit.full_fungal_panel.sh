#!/usr/bin/env bash
# Stage-3 launcher: comprehensive fungal panel benchmark (parameterized).
#
# Single launcher that supersedes submit.full_fungal_panel165.sh and
# submit.full_fungal_panel200.sh, which were 90% identical and differed only
# in QUERY_COUNT and the trailing genus list. Pass --panel <size> to pick
# a configuration; add a new case in select_panel() to onboard a future
# scale (e.g. --panel 500).
#
# Panel definitions:
#   --panel 165 : baseline 94-genus wish-list (AMF, EMF, endophyte,
#                 filamentous, yeast) selecting 165 query samples.
#   --panel 200 : 165 panel + 20 high-impact medical / agricultural /
#                 model-organism genera (Pyricularia rice blast,
#                 Puccinia wheat rust, Zymoseptoria two-speed wheat
#                 pathogen, Botrytis gray mold, Sclerotinia white mold,
#                 Verticillium wilt, Macrophomina charcoal rot,
#                 Cryphonectria chestnut blight, Cryptococcus meningitis,
#                 Histoplasma dimorphic, Trichophyton dermatophytes,
#                 Schizosaccharomyces fission yeast, Pichia industrial,
#                 Schizophyllum mating model, Metarhizium / Cordyceps
#                 entomopathogens, Aureobasidium black yeast, Microbotryum
#                 anther smut, Mucor + Neocallimastix early-diverging).
#                 Default panel size (manuscript-grade benchmark coverage).
#
# All three phases (bootstrap -> array -> combine) land in a fresh
# OUT_ROOT so prior runs stay immutable for comparison.
#
# NOTE: QUERY_GROUPS contains commas, which sbatch's --export= flag treats
# as a variable separator. We dodge that by `export`-ing everything in this
# shell and passing only `--export=ALL` to sbatch.
#
# Usage:
#   ./submit.full_fungal_panel.sh                        # default panel=200
#   ./submit.full_fungal_panel.sh --panel 165            # smaller panel
#   ./submit.full_fungal_panel.sh --panel 200 --dry-run  # print commands only
#
# Wall-time estimate: ~1.5-2 days with default array concurrency.

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${PROJECT_DIR}"

# ------------------------------------------------------------------- args
PANEL_SIZE=200
DRY_RUN=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --panel)   PANEL_SIZE="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1;       shift ;;
    -h|--help)
      grep '^#' "$0" | sed 's/^# \{0,1\}//' | head -40
      exit 0 ;;
    *) echo "unknown argument: $1 (use --panel <size> [--dry-run])" >&2; exit 2 ;;
  esac
done

# ----------------------------------------------------- panel configuration
# Shared base wish-list — AMF (Glomeromycota), EMF (Basidio + a few Asco),
# endophyte, filamentous Pezizomycotina, yeasts, and model-organism genera
# already used by the 13-sample / 165-sample validations.
BASE_GENERA="Acaulospora,Acremonium,Agaricus,Albahypha,Alternaria,Amanita,Ambispora,Archaeospora,Aspergillus,Austroboletus,Beauveria,Blaszkowskia,Boletus,Butyriboletus,Candida,Cantharellus,Cenococcum,Cetraspora,Chroogomphus,Cladosporium,Claroideoglomus,Clavulina,Colletotrichum,Cortinarius,Dentiscutata,Diversispora,Dominikia,Entrophospora,Eremothecium,Funneliformis,Fusarium,Geastrum,Geopora,Gigaspora,Glomus,Gomphidius,Gyroporus,Hebeloma,Helvella,Hydnellum,Hydnum,Hygrophorus,Inocybe,Laccaria,Lactarius,Lanmaoa,Leccinum,Lentinula,Meliniomyces,Monascus,mycorrhiza,Neurospora,Oidiodendron,Pacispora,Paecilomyces,Paraglomus,Paxillus,Penicillium,Pervetustus,Piriformospora,Pisolithus,Pleurotus,Pseudosperma,Racocetra,Redeckera,Rhizoglomus,Rhizophagus,Rhizopogon,Rhizopus,Rhizoscyphus,Russula,Saccharomyces,Sclerocystis,Scleroderma,Scutellospora,Septoglomus,Serendipita,Sieverdingia,Simiglomus,Sparassis,Sphaerosporella,Suillellus,Suillus,Talaromyces,Thelephora,Trichoderma,Tricholoma,Trichophaea,Tuber,Tulasnella,Tylopilus,Ustilago,Xerocomus,Yarrowia"

# Panel-200 additions: 20 medical / agricultural / model-organism genera.
PANEL200_EXTRAS="Aureobasidium,Botrytis,Cordyceps,Cryphonectria,Cryptococcus,Histoplasma,Macrophomina,Metarhizium,Microbotryum,Mucor,Neocallimastix,Pichia,Puccinia,Pyricularia,Schizophyllum,Schizosaccharomyces,Sclerotinia,Trichophyton,Verticillium,Zymoseptoria"

case "${PANEL_SIZE}" in
  165) QUERY_GROUPS="${BASE_GENERA}" ;;
  200) QUERY_GROUPS="${BASE_GENERA},${PANEL200_EXTRAS}" ;;
  *)   echo "unsupported panel size: ${PANEL_SIZE} (known: 165, 200)" >&2; exit 2 ;;
esac

# ----------------------------------------------------------------- exports
TIMESTAMP="${TIMESTAMP:-panel${PANEL_SIZE}_$(date +%Y%m%d_%H%M%S)}"
export OUT_ROOT="experiments/million_real/full_fungal_assembly_${TIMESTAMP}"
export QUERY_COUNT="${PANEL_SIZE}"
export QUERY_GROUPS
export MIN_ASSEMBLY_LEVEL=contig
export MYCOSV_BIOLOGY_TOP_N=5000

ARRAY_RANGE="0-$((QUERY_COUNT - 1))"
ARRAY_CONCURRENCY="${ARRAY_CONCURRENCY:-30}"
KEYWORD_COUNT=$(awk -F',' '{print NF}' <<< "${QUERY_GROUPS}")

JOB_PREFIX="panel${PANEL_SIZE}"
LOG_PREFIX="slurm-${JOB_PREFIX}"

echo "=== ${JOB_PREFIX} launch ==="
echo "  out_root:         ${OUT_ROOT}"
echo "  query_count:      ${QUERY_COUNT}"
echo "  array_range:      ${ARRAY_RANGE}  (concurrency=${ARRAY_CONCURRENCY})"
echo "  query_groups:     ${#QUERY_GROUPS} chars, ${KEYWORD_COUNT} keywords"
echo "  min_assembly_lvl: ${MIN_ASSEMBLY_LEVEL}"
echo "  biology_top_n:    ${MYCOSV_BIOLOGY_TOP_N}"
echo "  dry_run:          ${DRY_RUN}"
echo

# Wrapper that prints (dry) or submits and returns just the job ID.
submit() {
  if (( DRY_RUN )); then
    echo "  [dry] sbatch $*" >&2
    echo "9999999"
  else
    sbatch --parsable "$@"
  fi
}

# Phase A — bootstrap: runs prepare-million-real (downloads + indexes
# QUERY_COUNT queries, writes prepared/query_manifest.tsv + routing index
# + read-validation manifest), then BOOTSTRAP_ONLY=1 exits before the
# benchmark launches.
#
# Memory: 128G. panel-165 attempt 1 crashed at 44.15G MaxRSS on 32G;
# panel-200 adds ~21% scope, so 128G keeps the same ~3x safety margin.
# The C++ index builder is resumption-safe (skips clades whose .gbz
# already exists). multicore caps memory at 8192 MB/core, so 128G needs
# >=16 cores — the 8 cores beyond the C++ build's 8 threads are SLURM
# memory-budget headroom, not parallelism.
BOOTSTRAP_OVERRIDES="BOOTSTRAP_ONLY=1,SKIP_PREPARE=0"
BS=$(submit \
  --job-name="${JOB_PREFIX}-bootstrap" \
  -p multicore -c 16 --mem=128G --time=24:00:00 \
  --output="${LOG_PREFIX}-bootstrap-%j.out" \
  --export=ALL,${BOOTSTRAP_OVERRIDES} \
  submit.full_fungal_assembly.sh)
echo "bootstrap_jobid=${BS}"

# Phase B — N-way array, throttled to ARRAY_CONCURRENCY parallel tasks.
# SKIP_PREPARE=1 reuses prepared/ from bootstrap; FORCE_RERUN_SHARDS=1
# starts each shard from scratch (matters when bug fixes change call
# semantics). Comparator set (minigraph / svim-asm / anchorwave / PGGB /
# Cactus) is selected by submit.full_fungal_assembly.sh's RUN_PGGB /
# RUN_CACTUS defaults — both on by default for manuscript-grade coverage.
ARRAY_OVERRIDES="FORCE_RERUN_SHARDS=1,SKIP_PREPARE=1"
AR=$(submit \
  --dependency=afterok:${BS} \
  --array=${ARRAY_RANGE}%${ARRAY_CONCURRENCY} \
  --export=ALL,${ARRAY_OVERRIDES} \
  submit.full_fungal_assembly.sh)
echo "array_jobid=${AR}"

# Phase C — combine all shards into combined/ + master reports.
COMBINE_OVERRIDES="FULL_ASSEMBLY_SHARDS=1,RESUME_SHARDS=1,SKIP_PREPARE=1"
CB=$(submit \
  --dependency=afterany:${AR} \
  --job-name="${JOB_PREFIX}-combine" \
  --export=ALL,${COMBINE_OVERRIDES} \
  submit.full_fungal_assembly.sh)
echo "combine_jobid=${CB}"

echo
if (( DRY_RUN )); then
  echo "[dry-run] no jobs submitted; rerun without --dry-run to launch"
else
  echo "[submitted] chain: ${BS} -> ${AR} -> ${CB}"
  echo "[monitor]   squeue -u \$USER | grep ${JOB_PREFIX}"
  echo "[monitor]   tail -f ${LOG_PREFIX}-bootstrap-${BS}.out"
fi
