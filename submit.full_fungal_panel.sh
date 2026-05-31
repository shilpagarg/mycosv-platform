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
# Wall-time estimate: --panel 200 is comparator-free by default
# (ARRAY_MYCOSV_ONLY=1), so per-shard cost is MycoSV (seconds) + read-level
# validation + biology, ~5-15 min/shard. At the multicore_small QOS cap of ~16
# concurrent tasks (128 CPU / 8 CPU per task), 200 queries finish in ~3-6 h
# after bootstrap. --panel 165 keeps comparators ON by default (as before).
# Toggle either panel with ARRAY_MYCOSV_ONLY=0/1; with comparators on each
# shard adds up to ~4 h (90-min timeout each for minigraph/svim-asm/anchorwave)
# and the panel runs past a day.

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
ARRAY_CONCURRENCY="${ARRAY_CONCURRENCY:-60}"
ARRAY_CPUS="${ARRAY_CPUS:-8}"
ARRAY_MEM="${ARRAY_MEM:-40G}"
ARRAY_TIME="${ARRAY_TIME:-23:30:00}"
ARRAY_PARTITION="${ARRAY_PARTITION:-multicore_small,multicore}"
ARRAY_RUN_PGGB="${ARRAY_RUN_PGGB:-0}"
ARRAY_RUN_CACTUS="${ARRAY_RUN_CACTUS:-0}"
# Comparator-free is scoped to the 200-panel ONLY. That manuscript-grade panel
# rests on MycoSV graph-native calls + independent read-level validation;
# minigraph/svim-asm/anchorwave each carry a 90-min timeout and were the
# dominant per-shard wall-time cost over 200 queries, so they are disabled.
#
# Every other entry point keeps comparators and runs end-to-end as before:
#   - --panel 165 (and any future smaller panel) defaults to comparators ON.
#   - The 5-sample table1/table2 assembly + long-read runs invoke
#     submit.full_fungal_assembly.sh directly without MYCOSV_ONLY, so the
#     worker's MYCOSV_ONLY:-0 default leaves their comparator path untouched
#     (--run-minigraph/svim-asm/anchorwave[/pggb/cactus], sniffles/cutesv/svim).
#
# Override either way with ARRAY_MYCOSV_ONLY=0 (force comparators on) or =1
# (force comparator-free); when unset the panel size picks the default below.
if [[ -z "${ARRAY_MYCOSV_ONLY:-}" ]]; then
  if [[ "${PANEL_SIZE}" == "200" ]]; then
    ARRAY_MYCOSV_ONLY=1
  else
    ARRAY_MYCOSV_ONLY=0
  fi
fi
# Cap each assembly query to its N longest contigs. MycoSV's hierarchical
# assembly caller stalls for hours above a few hundred query contigs, and 85 of
# the 200 panel queries are fragmented public drafts (up to 181k contigs). The
# cap keeps every query in the panel as its N-longest-contig subset so the run
# finishes inside the wall-time budget; well-assembled queries (<=N contigs) are
# untouched. 0 = original behaviour (process every contig). Override with
# ARRAY_MAX_QUERY_CONTIGS=<N>.
ARRAY_MAX_QUERY_CONTIGS="${ARRAY_MAX_QUERY_CONTIGS:-150}"
KEYWORD_COUNT=$(awk -F',' '{print NF}' <<< "${QUERY_GROUPS}")

JOB_PREFIX="panel${PANEL_SIZE}"
LOG_PREFIX="slurm-${JOB_PREFIX}"

echo "=== ${JOB_PREFIX} launch ==="
echo "  out_root:         ${OUT_ROOT}"
echo "  query_count:      ${QUERY_COUNT}"
echo "  array_range:      ${ARRAY_RANGE}  (concurrency=${ARRAY_CONCURRENCY})"
echo "  array_resources:  ${ARRAY_PARTITION}, ${ARRAY_CPUS} CPUs, ${ARRAY_MEM}, ${ARRAY_TIME}"
echo "  heavy_graphs:     pggb=${ARRAY_RUN_PGGB}, cactus=${ARRAY_RUN_CACTUS}"
echo "  mycosv_only:      ${ARRAY_MYCOSV_ONLY} (1=no comparators, read-level validation only)"
echo "  max_query_contigs:${ARRAY_MAX_QUERY_CONTIGS} (keep N longest contigs/query; 0=all)"
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
# Memory: 256G. Panel-200 job 15448995 crashed at 116G MaxRSS / 128G cap
# during routing-index build (std::bad_alloc) because the C++ builder
# launched 16 parallel workers each holding up to ~1.5 GB of full-base
# graph state, plus glibc per-arena fragmentation across ~9k full-base
# shards. 256G fits well inside `multicore` (1.5 TB/node, MaxMemPerNode
# UNLIMITED) — no himem node needed. Tiered worker pools, streaming FASTA
# reads, and centroid roll-up in fungi_tol_bridge.hpp cut peak RSS further.
#
# MALLOC_ARENA_MAX / MALLOC_TRIM_THRESHOLD_: tame glibc malloc per-thread
# arenas (default = 8 × CPU_COUNT, so 16 cores × 8 = 128 arenas, each
# retaining freed pages). Capping at 2 arenas and a tight trim threshold
# returns memory to the kernel more aggressively — cheap substitute for
# jemalloc on a node that doesn't have it installed.
#
# Cores=32 not 16: multicore enforces 8192 MB/CPU, so 256G needs ≥32
# cores. The C++ index builder still uses --tol-index-threads=16 (see
# submit.full_fungal_assembly.sh); the extra 16 cores are pure memory
# headroom, not extra parallelism.
BOOTSTRAP_OVERRIDES="BOOTSTRAP_ONLY=1,SKIP_PREPARE=0,MALLOC_ARENA_MAX=2,MALLOC_TRIM_THRESHOLD_=131072"
BS=$(submit \
  --job-name="${JOB_PREFIX}-bootstrap" \
  -p multicore -c 32 --mem=256G --time=24:00:00 \
  --output="${LOG_PREFIX}-bootstrap-%j.out" \
  --export=ALL,${BOOTSTRAP_OVERRIDES} \
  submit.full_fungal_assembly.sh)
echo "bootstrap_jobid=${BS}"

# Phase B — N-way array, throttled to ARRAY_CONCURRENCY parallel tasks.
# SKIP_PREPARE=1 reuses prepared/ from bootstrap; FORCE_RERUN_SHARDS=1
# starts each shard from scratch (matters when bug fixes change call
# semantics). The default catch-up profile keeps the always-on assembly
# comparators (minigraph / svim-asm / anchorwave) and disables PGGB/Cactus,
# which were the main wall-time and scheduling pressure for the 200-sample
# panel. Set ARRAY_RUN_PGGB=1 ARRAY_RUN_CACTUS=1 for a slower, full graph
# comparator rerun.
ARRAY_OVERRIDES="FORCE_RERUN_SHARDS=1,SKIP_PREPARE=1,RUN_PGGB=${ARRAY_RUN_PGGB},RUN_CACTUS=${ARRAY_RUN_CACTUS},MYCOSV_ONLY=${ARRAY_MYCOSV_ONLY},MAX_ASSEMBLY_QUERY_CONTIGS_KEEP=${ARRAY_MAX_QUERY_CONTIGS}"
AR=$(submit \
  --dependency=afterok:${BS} \
  --job-name="${JOB_PREFIX}-array" \
  -p "${ARRAY_PARTITION}" -c "${ARRAY_CPUS}" --mem="${ARRAY_MEM}" --time="${ARRAY_TIME}" \
  --array=${ARRAY_RANGE}%${ARRAY_CONCURRENCY} \
  --export=ALL,${ARRAY_OVERRIDES} \
  submit.full_fungal_assembly.sh)
echo "array_jobid=${AR}"

# Phase C — combine all shards into combined/ + master reports.
COMBINE_OVERRIDES="FULL_ASSEMBLY_SHARDS=1,RESUME_SHARDS=1,SKIP_PREPARE=1,MYCOSV_ONLY=${ARRAY_MYCOSV_ONLY},MAX_ASSEMBLY_QUERY_CONTIGS_KEEP=${ARRAY_MAX_QUERY_CONTIGS}"
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
