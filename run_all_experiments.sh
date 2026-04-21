#!/usr/bin/env bash
# run_all_experiments.sh
# Master script to run all experiments (small-scale, million-scale, and real fungal data)
# with complete intermediate file preservation.
#
# Usage:
#   bash run_all_experiments.sh                    # Run all experiments
#   bash run_all_experiments.sh --small            # Small-scale only
#   bash run_all_experiments.sh --large            # Large-scale/million-scale only
#   bash run_all_experiments.sh --real             # Real panels only (per-species NCBI + ENA)
#   bash run_all_experiments.sh --million-real     # Only the real million-scale NCBI index build
#
# Environment overrides for --million-real:
#   MILLION_REAL_MAX_ASSEMBLIES=1000    # how many NCBI fungal assemblies to download
#   MILLION_REAL_TARGET_CENTROIDS=1000000  # total centroids after decoy padding (0 = no padding)
#
# IMPORTANT: we intentionally do NOT use `set -e` at the script level for the
# real-data stage. A single panel/network/tool failure must not abort the rest
# of the matrix — each panel is guarded with `|| true` and its log is kept for
# inspection. Fatal errors per stage are still reported in the summary.

set -u
set -o pipefail

WORK_DIR="/mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/AMF/scale"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
EXPERIMENT_TYPE="${1:-all}"  # all, small, large, real, million-real

# Strip leading dashes for convenience: --small and small are both accepted.
EXPERIMENT_TYPE="${EXPERIMENT_TYPE#--}"

# How many real NCBI fungal assemblies to download when building the real
# million-scale routing index. Override with env var when needed.
MILLION_REAL_MAX_ASSEMBLIES="${MILLION_REAL_MAX_ASSEMBLIES:-1000}"
MILLION_REAL_TARGET_CENTROIDS="${MILLION_REAL_TARGET_CENTROIDS:-1000000}"

# Create experiment directories
SMALL_DIR="${WORK_DIR}/experiments/small_tests/${TIMESTAMP}"
LARGE_DIR="${WORK_DIR}/experiments/large_scale/${TIMESTAMP}"
REAL_DIR="${WORK_DIR}/experiments/real_data/${TIMESTAMP}"
MILLION_REAL_DIR="${WORK_DIR}/experiments/million_real/${TIMESTAMP}"

mkdir -p "${SMALL_DIR}" "${LARGE_DIR}" "${REAL_DIR}" "${MILLION_REAL_DIR}" "${MILLION_REAL_DIR}"

cd "${WORK_DIR}"

# Color codes for output
GREEN='\033[0;32m'
RED='\033[0;31m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Track stage outcomes so we can still print a meaningful summary even if
# individual panels or modes fail.
FAILED_STAGES=()
SUCCESS_STAGES=()

mark_success() { SUCCESS_STAGES+=("$1"); }
mark_failure() { FAILED_STAGES+=("$1"); echo -e "${RED}✗ $1${NC}"; }

echo -e "${BLUE}========================================"
echo "MycoSV Comprehensive Experiments"
echo "========================================${NC}"
echo "Start time: $(date)"
echo "Timestamp: ${TIMESTAMP}"
echo "Working directory: ${WORK_DIR}"
echo "Experiment type: ${EXPERIMENT_TYPE}"
echo ""

# Scenario sets that guarantee all 5 SV types (INS, DEL, DUP, INV, TRA).
# compact_yeast alone only generates DEL+INS, so we always include at least
# one of {te_rich_pathogen, cross_phylum_hgt} where all five types are emitted.
SMALL_SCENARIO_SET="compact_yeast,te_rich_pathogen,cross_phylum_hgt"
LARGE_SCENARIO_SET="te_rich_pathogen,cross_phylum_hgt"

# ============================================================================
# 1. SMALL-SCALE SIMULATED DATA TESTS
# ============================================================================

if [[ "$EXPERIMENT_TYPE" == "all" || "$EXPERIMENT_TYPE" == "small" ]]; then
  echo -e "${YELLOW}[1/5] Running small-scale simulated data tests...${NC}"
  echo "      Output: ${SMALL_DIR}/simulated"
  mkdir -p "${SMALL_DIR}/simulated"

  if python3 -m pytest \
      test_pipeline_features.py \
      test_amf.py \
      test_all_use_cases.py \
      -v \
      --tb=short \
      --basetemp="${SMALL_DIR}/simulated/.pytest_tmp" \
      2>&1 | tee "${SMALL_DIR}/simulated/pytest_output.log"; then
    mark_success "small.simulated_pytest"
    echo -e "${GREEN}✓ Small-scale simulated tests complete${NC}"
  else
    mark_failure "small.simulated_pytest"
  fi
  echo ""
fi

# ============================================================================
# 2. SMALL-SCALE REAL FUNGAL DATA TESTS + PR METRICS
# ============================================================================

if [[ "$EXPERIMENT_TYPE" == "all" || "$EXPERIMENT_TYPE" == "small" ]]; then
  echo -e "${YELLOW}[2/5] Running small-scale real fungal data tests...${NC}"
  echo "      Output: ${SMALL_DIR}/real_data"
  mkdir -p "${SMALL_DIR}/real_data"

  if python3 -m pytest \
      test_real_fungal_benchmark.py \
      test_new_biology_candidates.py \
      -v \
      --tb=short \
      --basetemp="${SMALL_DIR}/real_data/.pytest_tmp" \
      2>&1 | tee "${SMALL_DIR}/real_data/pytest_output.log"; then
    mark_success "small.real_pytest"
    echo -e "${GREEN}✓ Small-scale real data tests complete${NC}"
  else
    mark_failure "small.real_pytest"
  fi
  echo ""

  # Run small-scale benchmark to actually produce precision/recall TSVs.
  # Pytest alone leaves only VCFs; the benchmark driver produces pr_metrics.*.
  echo -e "${YELLOW}[2b/5] Generating small-scale precision/recall metrics...${NC}"
  echo "      Output: ${SMALL_DIR}/benchmarks"
  mkdir -p "${SMALL_DIR}/benchmarks"

  if python3 run_million_mode_query_benchmark.py \
      --out-dir "${SMALL_DIR}/benchmarks" \
      --modes assembly,short-reads,long-reads \
      --scenario-set "${SMALL_SCENARIO_SET}" \
      --n-centroids 10000 \
      --n-genomes 4 \
      --n-reps 2 \
      --seed 42 \
      2>&1 | tee "${SMALL_DIR}/benchmarks/benchmark.log"; then
    mark_success "small.pr_metrics_benchmark"
    echo -e "${GREEN}✓ Small-scale benchmark metrics complete${NC}"
  else
    mark_failure "small.pr_metrics_benchmark"
  fi

  # Post-process pytest VCFs -> precision/recall TSVs so small_tests/*/metrics
  # has the same shape as the large-scale output regardless of whether the
  # million-mode driver succeeded above.
  echo "      Post-processing pytest VCFs into PR TSVs..."
  if bash generate_small_test_metrics.sh "${TIMESTAMP}" \
      2>&1 | tee "${SMALL_DIR}/metrics_postprocess.log"; then
    mark_success "small.pr_metrics_postprocess"
  else
    mark_failure "small.pr_metrics_postprocess"
  fi
  echo ""
fi

# ============================================================================
# 3. MILLION-SCALE SIMULATED DATA BENCHMARK (all 5 SV types)
# ============================================================================

if [[ "$EXPERIMENT_TYPE" == "all" || "$EXPERIMENT_TYPE" == "large" ]]; then
  echo -e "${YELLOW}[3/5] Running million-scale simulated data benchmark...${NC}"
  echo "      Modes: assembly, short-reads, long-reads"
  echo "      Scenarios: ${LARGE_SCENARIO_SET} (covers INS/DEL/DUP/INV/TRA)"
  echo "      Output: ${LARGE_DIR}/million_scale_simulated"
  mkdir -p "${LARGE_DIR}/million_scale_simulated"

  if python3 run_million_mode_query_benchmark.py \
      --out-dir "${LARGE_DIR}/million_scale_simulated" \
      --modes assembly,short-reads,long-reads \
      --scenario-set "${LARGE_SCENARIO_SET}" \
      --n-centroids 1000000 \
      --n-genomes 8 \
      --n-reps 3 \
      --seed 42 \
      2>&1 | tee "${LARGE_DIR}/million_scale_simulated/benchmark.log"; then
    mark_success "large.million_scale"
    echo -e "${GREEN}✓ Million-scale simulated benchmark complete${NC}"
  else
    mark_failure "large.million_scale"
  fi
  echo ""
fi

# ============================================================================
# 4. MILLION-SCALE MODE PRECISION-RECALL BENCHMARK
# ============================================================================

if [[ "$EXPERIMENT_TYPE" == "all" || "$EXPERIMENT_TYPE" == "large" ]]; then
  echo -e "${YELLOW}[4/5] Running mode precision-recall benchmarks...${NC}"
  echo "      Output: ${LARGE_DIR}/mode_pr_benchmark"
  mkdir -p "${LARGE_DIR}/mode_pr_benchmark"

  for mode in assembly short-reads long-reads; do
    echo "  - Benchmarking mode: ${mode}"
    mkdir -p "${LARGE_DIR}/mode_pr_benchmark/${mode}"

    if python3 run_mode_pr_benchmark.py \
        --modes "${mode}" \
        --out-dir "${LARGE_DIR}/mode_pr_benchmark/${mode}" \
        --scenario-set "${LARGE_SCENARIO_SET}" \
        --n-refs 500 \
        --n-queries 20 \
        --seed 42 \
        2>&1 | tee "${LARGE_DIR}/mode_pr_benchmark/${mode}/benchmark.log"; then
      mark_success "large.mode_pr.${mode}"
    else
      mark_failure "large.mode_pr.${mode}"
    fi
  done

  echo -e "${GREEN}✓ Mode precision-recall benchmarks complete${NC}"
  echo ""
fi

# ============================================================================
# 4b. MILLION-SCALE *REAL* FUNGAL INDEX (downloads NCBI assemblies)
# ============================================================================
#
# This stage downloads up to N real fungal assemblies from NCBI RefSeq, builds
# a real MycoSV routing index over them, then pads the routing store with
# synthetic decoys up to --target-centroids. This is the bridge between the
# small per-panel real-data downloads and the previously-all-synthetic
# million-scale benchmark. Disk/bandwidth heavy — opt-in via `million-real`
# or `all`.
# ============================================================================

if [[ "$EXPERIMENT_TYPE" == "all" || "$EXPERIMENT_TYPE" == "million-real" ]]; then
  echo -e "${YELLOW}[4b] Building real million-scale fungal routing index...${NC}"
  echo "      Downloading up to ${MILLION_REAL_MAX_ASSEMBLIES} NCBI RefSeq assemblies"
  echo "      Target centroids (real+decoys): ${MILLION_REAL_TARGET_CENTROIDS}"
  echo "      Output: ${MILLION_REAL_DIR}"

  if python3 run_real_fungal_benchmark.py prepare-million-real \
      --out-dir "${MILLION_REAL_DIR}" \
      --source ncbi-refseq \
      --max-assemblies "${MILLION_REAL_MAX_ASSEMBLIES}" \
      --target-centroids "${MILLION_REAL_TARGET_CENTROIDS}" \
      --min-assembly-level scaffold \
      --latest-only \
      --threads 4 \
      --seed 42 \
      2>&1 | tee "${MILLION_REAL_DIR}/prepare_million_real.log"; then
    mark_success "million_real.index_build"
    echo -e "${GREEN}✓ Real million-scale index ready${NC}"
  else
    mark_failure "million_real.index_build"
  fi
  echo ""
fi

# ============================================================================
# 5. REAL FUNGAL DATA BENCHMARKS
# ============================================================================
#
# Each panel is fully isolated: a failure in one panel (download error,
# missing tool, empty selection) does not stop the rest.
#
# State-of-the-art comparators are enabled per mode:
#   assembly    : SyRI, minigraph+gfatools, PGGB, Minigraph-Cactus,
#                 SVIM-asm, AnchorWave
#   short-reads : Delly, Manta
#   long-reads  : SVIM, Sniffles, cuteSV
# ============================================================================

if [[ "$EXPERIMENT_TYPE" == "all" || "$EXPERIMENT_TYPE" == "real" ]]; then
  echo -e "${YELLOW}[5/5] Running real fungal data benchmarks...${NC}"
  echo "      Panels: compact_yeast, amf_large, cross_phylum_hgt, te_rich_pathogen, two_speed_pathogen"
  echo "      Output: ${REAL_DIR}"

  PANELS=("compact_yeast" "amf_large" "cross_phylum_hgt" "te_rich_pathogen" "two_speed_pathogen")

  for panel in "${PANELS[@]}"; do
    echo ""
    echo "  ── Processing panel: ${panel} ──"
    PANEL_DIR="${REAL_DIR}/${panel}"
    mkdir -p "${PANEL_DIR}"

    # Prepare real data — guarded so one panel's failure does not abort the rest.
    # --query-mode mixed asks the preparer to also pull public ENA read runs
    # for each panel species (short-reads + long-reads), which is what makes
    # benchmark_short-reads/ and benchmark_long-reads/ produce actual results.
    # Previously those folders were silently empty because the NCBI panel
    # preparer only wrote assembly-mode query rows.
    echo "    - Preparing real data (mixed: assembly + short-reads + long-reads)..."
    if python3 run_real_fungal_benchmark.py prepare \
        --out-dir "${PANEL_DIR}/prepared" \
        --panel "${panel}" \
        --max-assemblies-per-species 8 \
        --querys-per-species 5 \
        --max-ref-downloads 20 \
        --max-query-downloads 10 \
        --query-mode mixed \
        --read-accessions-per-species 2 \
        --allow-no-queries \
        2>&1 | tee "${PANEL_DIR}/prepare.log"; then
      mark_success "real.${panel}.prepare"
    else
      mark_failure "real.${panel}.prepare"
      echo "    - Skipping benchmarks for ${panel} (prepare failed)"
      continue
    fi

    # Must have at least one manifest + queries produced by prepare.
    if [[ ! -f "${PANEL_DIR}/prepared/selected_catalog.tsv" \
       && ! -f "${PANEL_DIR}/prepared/reference_catalog.tsv" ]]; then
      echo "    - Skipping benchmarks for ${panel} (no catalog produced)"
      mark_failure "real.${panel}.no_catalog"
      continue
    fi
    if [[ ! -f "${PANEL_DIR}/prepared/query_manifest.tsv" ]]; then
      echo "    - Skipping benchmarks for ${panel} (no query manifest produced)"
      mark_failure "real.${panel}.no_queries"
      continue
    fi

    for mode in assembly short-reads long-reads; do
      echo "    - Benchmarking mode: ${mode}..."
      mkdir -p "${PANEL_DIR}/benchmark_${mode}"

      # Mode-specific SOTA comparator flags.
      comparator_flags=()
      case "${mode}" in
        assembly)
          # Assembly-to-assembly and pangenome-graph comparators:
          #   SyRI + minigraph + PGGB are the original trio;
          #   Minigraph-Cactus (cactus-pangenome), SVIM-asm, and AnchorWave
          #   are fungi-oriented / pangenome-oriented additions that run
          #   whenever the binary is on $PATH (each adapter no-ops if not).
          comparator_flags+=(--run-syri --run-minigraph --run-pggb)
          comparator_flags+=(--run-cactus --run-svim-asm --run-anchorwave)
          ;;
        short-reads)
          comparator_flags+=(--run-delly --run-manta)
          ;;
        long-reads)
          comparator_flags+=(--run-svim --run-sniffles --run-cutesv)
          ;;
      esac

      if python3 run_real_fungal_benchmark.py benchmark \
          --prepared-dir "${PANEL_DIR}/prepared" \
          --mode "${mode}" \
          --out-dir "${PANEL_DIR}/benchmark_${mode}" \
          --threads 4 \
          "${comparator_flags[@]}" \
          2>&1 | tee "${PANEL_DIR}/benchmark_${mode}.log"; then
        mark_success "real.${panel}.${mode}"
      else
        mark_failure "real.${panel}.${mode}"
      fi
    done
  done

  echo ""
  echo -e "${GREEN}✓ Real fungal data benchmarks stage complete${NC}"
  echo ""
fi

# ============================================================================
# SUMMARY REPORT
# ============================================================================

echo -e "${BLUE}========================================"
echo "Experiment Summary"
echo "========================================${NC}"
echo ""
echo "All experiment outputs saved to:"
echo "  Small-scale:   ${SMALL_DIR}"
echo "  Large-scale:   ${LARGE_DIR}"
echo "  Real data:     ${REAL_DIR}"
echo "  Million-real:  ${MILLION_REAL_DIR}"
echo ""
echo "Stage outcomes:"
echo "  Succeeded: ${#SUCCESS_STAGES[@]}"
echo "  Failed:    ${#FAILED_STAGES[@]}"
if [[ ${#FAILED_STAGES[@]} -gt 0 ]]; then
  echo "  Failed stages:"
  for s in "${FAILED_STAGES[@]}"; do echo "    - $s"; done
fi
echo ""
echo "Intermediate files preserved:"
find "${SMALL_DIR}" "${LARGE_DIR}" "${REAL_DIR}" "${MILLION_REAL_DIR}" -type f \
  \( -name "*.vcf" -o -name "*.vcf.gz" -o -name "*.tsv" -o -name "*.fasta" -o -name "*.fastq" \) \
  2>/dev/null | wc -l | xargs echo "  Total:"
echo ""
echo "Disk usage:"
du -sh "${SMALL_DIR}" "${LARGE_DIR}" "${REAL_DIR}" "${MILLION_REAL_DIR}" 2>/dev/null | awk '{print "  " $0}'
echo ""
echo "Log files:"
find "${SMALL_DIR}" "${LARGE_DIR}" "${REAL_DIR}" "${MILLION_REAL_DIR}" -name "*.log" 2>/dev/null | wc -l | xargs echo "  Total:"
echo ""
echo -e "${GREEN}✓ All experiments complete! End time: $(date)${NC}"
echo ""
echo "Next steps:"
echo "  1. Review logs for errors: grep -r 'ERROR' ${SMALL_DIR} ${LARGE_DIR} ${REAL_DIR} ${MILLION_REAL_DIR}"
echo "  2. Analyze results: python3 analyze_results.py --input-dir ${LARGE_DIR}"
echo "  3. Generate report: bash run_comprehensive_experiments.sh"
echo ""

# Exit non-zero only if *everything* in a requested stage failed, so callers
# (CI, cron) can distinguish "some panels flaky" from "nothing ran".
if [[ ${#SUCCESS_STAGES[@]} -eq 0 ]]; then
  exit 1
fi
exit 0
