#!/usr/bin/env bash
# run_all_experiments.sh
# Master script to run all experiments (simulated benchmark and real fungal data).
#
# Usage:
#   bash run_all_experiments.sh                    # Run all experiments
#   bash run_all_experiments.sh --simulated        # Simulated benchmark only
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
EXPERIMENT_TYPE="${1:-all}"  # all, simulated, real, million-real

# Shared download cache: FASTA files downloaded here are reused across runs.
DATA_CACHE_DIR="${DATA_CACHE_DIR:-${WORK_DIR}/data_cache}"
mkdir -p "${DATA_CACHE_DIR}"

# Strip leading dashes for convenience: --simulated and simulated are both accepted.
EXPERIMENT_TYPE="${EXPERIMENT_TYPE#--}"

# How many real NCBI fungal assemblies to download when building the real
# million-scale routing index. Override with env var when needed.
MILLION_REAL_MAX_ASSEMBLIES="${MILLION_REAL_MAX_ASSEMBLIES:-10000}"
MILLION_REAL_TARGET_CENTROIDS="${MILLION_REAL_TARGET_CENTROIDS:-1000000}"

# Worker threads passed to every tool that accepts a thread count:
#   MycoSV binary (--threads / --tol-index-threads), minimap2 (-t),
#   samtools sort (-@), sniffles2 (--threads), cuteSV (--threads),
#   Delly (OMP_NUM_THREADS), Manta (-j), minigraph (-t), pggb (-t).
THREADS="${THREADS:-32}"

# Long-read platform used for both simulated and real experiments.
#   hifi    PacBio HiFi CCS (Revio / Sequel IIe) — 15 kb reads, ≥Q20 accuracy.
#           minimap2 map-hifi → sniffles2 / cuteSV (HiFi params) / SVIM.
#   ont-r10 ONT R10.4.1 standard simplex — 10 kb, ~Q20 on PromethION/GridION.
#           minimap2 map-ont → sniffles2 --long-read-model ont_r10_q20 (v2.2+).
#           WhatsHap phase+haplotag: applicable for dikaryotic fungi such as
#           Puccinia spp., Leptosphaeria maculans, Zymoseptoria tritici.
#   ont-r9  ONT R9.4.1 legacy — 8 kb, ~Q15.  Still prevalent in public ENA data.
LONG_READ_PLATFORM="${LONG_READ_PLATFORM:-ont-r10}"

# Create experiment directories
SIM_DIR="${WORK_DIR}/experiments/simulated/${TIMESTAMP}"
REAL_DIR="${WORK_DIR}/experiments/real_data/${TIMESTAMP}"
MILLION_REAL_DIR="${WORK_DIR}/experiments/million_real/${TIMESTAMP}"

mkdir -p "${SIM_DIR}" "${REAL_DIR}" "${MILLION_REAL_DIR}"

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
# compact_yeast:              DEL + INS
# two_speed_pathogen_extreme: INV + TRA + INS
# arbuscular_mf:              DUP + INS
SIM_SCENARIO_SET="compact_yeast,two_speed_pathogen_extreme,arbuscular_mf"

# ============================================================================
# 1. SIMULATED BENCHMARK (precision/recall at million scale, hundreds of SVs)
#
# Parameters chosen to produce ≥270 truth SVs:
#   n_genomes=30, n_reps=3, n_contigs=10, total_len=200000
#   → (30-3)=27 query genomes × 10 contigs = 270 SVs across 3 scenarios
# ============================================================================

if [[ "$EXPERIMENT_TYPE" == "all" || "$EXPERIMENT_TYPE" == "simulated" ]]; then
  echo -e "${YELLOW}[1/4] Running simulated benchmark (million-scale routing, 270 SVs)...${NC}"
  echo "      Modes: assembly, short-reads, long-reads"
  echo "      Scenarios: ${SIM_SCENARIO_SET} (covers INS/DEL/DUP/INV/TRA)"
  echo "      Genomes: 30 total (3 refs + 27 queries), 10 contigs each"
  echo "      Output: ${SIM_DIR}/benchmarks"
  mkdir -p "${SIM_DIR}/benchmarks"

  if python3 run_million_mode_query_benchmark.py \
      --out-dir "${SIM_DIR}/benchmarks" \
      --modes assembly,short-reads,long-reads \
      --scenario-set "${SIM_SCENARIO_SET}" \
      --n-centroids 1000000 \
      --n-genomes 30 \
      --n-reps 3 \
      --n-contigs 10 \
      --total-len 200000 \
      --seed 42 \
      --target-svs-per-scenario 3000 \
      --min-contig-bp 12000 \
      --long-read-platform "${LONG_READ_PLATFORM}" \
      --threads "${THREADS}" \
      2>&1 | tee "${SIM_DIR}/benchmarks/benchmark.log"; then
    mark_success "simulated.pr_metrics_benchmark"
    echo -e "${GREEN}✓ Simulated benchmark complete${NC}"
  else
    mark_failure "simulated.pr_metrics_benchmark"
  fi
  echo ""
fi

# ============================================================================
# 2. MILLION-SCALE *REAL* FUNGAL INDEX (downloads NCBI assemblies)
# ============================================================================

if [[ "$EXPERIMENT_TYPE" == "all" || "$EXPERIMENT_TYPE" == "million-real" ]]; then
  echo -e "${YELLOW}[2/4] Building real million-scale fungal routing index...${NC}"
  echo "      Downloading up to ${MILLION_REAL_MAX_ASSEMBLIES} NCBI GenBank assemblies (contig level or better)"
  echo "      Target centroids (real+decoys): ${MILLION_REAL_TARGET_CENTROIDS}"
  echo "      Output: ${MILLION_REAL_DIR}"

  if python3 run_real_fungal_benchmark.py prepare-million-real \
      --out-dir "${MILLION_REAL_DIR}" \
      --source ncbi-genbank \
      --max-assemblies "${MILLION_REAL_MAX_ASSEMBLIES}" \
      --target-centroids "${MILLION_REAL_TARGET_CENTROIDS}" \
      --min-assembly-level contig \
      --threads "${THREADS}" \
      --seed 42 \
      --data-cache-dir "${DATA_CACHE_DIR}" \
      2>&1 | tee "${MILLION_REAL_DIR}/prepare_million_real.log"; then
    mark_success "million_real.index_build"
    echo -e "${GREEN}✓ Real million-scale index ready${NC}"
  else
    mark_failure "million_real.index_build"
  fi
  echo ""
fi

# ============================================================================
# 3. REAL FUNGAL DATA BENCHMARKS
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
  echo -e "${YELLOW}[3/4] Running real fungal data benchmarks...${NC}"
  echo "      Panels: compact_yeast, amf_large, cross_phylum_hgt, te_rich_pathogen, two_speed_pathogen"
  echo "      Output: ${REAL_DIR}"

  PANELS=("compact_yeast" "amf_large" "cross_phylum_hgt" "te_rich_pathogen" "two_speed_pathogen")

  for panel in "${PANELS[@]}"; do
    echo ""
    echo "  ── Processing panel: ${panel} ──"
    PANEL_DIR="${REAL_DIR}/${panel}"
    mkdir -p "${PANEL_DIR}"

    echo "    - Preparing real data (mixed: assembly + short-reads + long-reads)..."
    if python3 run_real_fungal_benchmark.py prepare \
        --out-dir "${PANEL_DIR}/prepared" \
        --panel "${panel}" \
        --source ncbi-genbank \
        --max-assemblies-per-species 20 \
        --querys-per-species 5 \
        --max-ref-downloads 20 \
        --max-query-downloads 10 \
        --query-mode mixed \
        --read-accessions-per-species 2 \
        --allow-no-queries \
        --data-cache-dir "${DATA_CACHE_DIR}" \
        2>&1 | tee "${PANEL_DIR}/prepare.log"; then
      mark_success "real.${panel}.prepare"
    else
      mark_failure "real.${panel}.prepare"
      echo "    - Skipping benchmarks for ${panel} (prepare failed)"
      continue
    fi

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
    if [[ $(wc -l < "${PANEL_DIR}/prepared/query_manifest.tsv") -le 1 ]]; then
      echo "    - Skipping benchmarks for ${panel} (query manifest is empty)"
      mark_failure "real.${panel}.empty_queries"
      continue
    fi

    for mode in assembly short-reads long-reads; do
      echo "    - Benchmarking mode: ${mode}..."
      mkdir -p "${PANEL_DIR}/benchmark_${mode}"

      comparator_flags=()
      case "${mode}" in
        assembly)
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
          --threads "${THREADS}" \
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
# 4. VISUALIZATION REPORT
# ============================================================================
#
# Builds an integrated HTML report across simulated benchmarks, real-data SV
# results, and biological findings when available. Missing inputs are handled
# gracefully so this stage never blocks the rest of the experiment matrix.
# ============================================================================

REPORT_DIR="${WORK_DIR}/experiments/reports/${TIMESTAMP}"
mkdir -p "${REPORT_DIR}"

if [[ "$EXPERIMENT_TYPE" == "all" || "$EXPERIMENT_TYPE" == "simulated" || "$EXPERIMENT_TYPE" == "real" ]]; then
  echo -e "${YELLOW}[4/4] Generating visualization report...${NC}"
  echo "      Output: ${REPORT_DIR}"

  SIM_RESULTS=""
  if [[ -f "${SIM_DIR}/benchmarks/million_mode_summary.tsv" ]]; then
    SIM_RESULTS="${SIM_DIR}/benchmarks/million_mode_summary.tsv"
  elif [[ -f "${SIM_DIR}/benchmarks/pr_metrics_simulated_summary.tsv" ]]; then
    SIM_RESULTS="${SIM_DIR}/benchmarks/pr_metrics_simulated_summary.tsv"
  fi

  REAL_RESULTS="${REPORT_DIR}/real_merged.tsv"
  BIO_RESULTS="${REPORT_DIR}/biology_merged.tsv"
  : > "${REAL_RESULTS}"
  : > "${BIO_RESULTS}"

  merge_tsv_group() {
    local out_file="$1"
    shift
    local first_written=0
    local f
    for f in "$@"; do
      [[ -f "$f" ]] || continue
      if [[ $first_written -eq 0 ]]; then
        cat "$f" >> "$out_file"
        first_written=1
      else
        tail -n +2 "$f" >> "$out_file"
      fi
    done
    return 0
  }

  mapfile -t REAL_TSVS < <(find "${REAL_DIR}" -type f \( \
      -name "*summary*.tsv" -o \
      -name "*pr_metrics*.tsv" -o \
      -name "*normalized_calls*.tsv" -o \
      -name "*score*.tsv" \
    \) 2>/dev/null | sort)

  mapfile -t BIO_TSVS < <(find "${REAL_DIR}" -type f \( \
      -name "*biology*.tsv" -o \
      -name "*candidate*.tsv" -o \
      -name "*annotation*.tsv" -o \
      -name "*pathway*.tsv" \
    \) 2>/dev/null | sort)

  if [[ ${#REAL_TSVS[@]} -gt 0 ]]; then
    merge_tsv_group "${REAL_RESULTS}" "${REAL_TSVS[@]}"
  fi
  if [[ ${#BIO_TSVS[@]} -gt 0 ]]; then
    merge_tsv_group "${BIO_RESULTS}" "${BIO_TSVS[@]}"
  fi

  report_cmd=(python3 sv_visualization_report.py --outdir "${REPORT_DIR}" --title "MycoSV comprehensive report (${TIMESTAMP})")
  [[ -n "${SIM_RESULTS}" ]] && report_cmd+=(--simulated "${SIM_RESULTS}")
  [[ -s "${REAL_RESULTS}" ]] && report_cmd+=(--real "${REAL_RESULTS}")
  [[ -s "${BIO_RESULTS}" ]] && report_cmd+=(--biology "${BIO_RESULTS}")

  if [[ -f "${WORK_DIR}/sv_visualization_report.py" ]]; then
    if "${report_cmd[@]}" 2>&1 | tee "${REPORT_DIR}/report.log"; then
      if [[ -f "${REPORT_DIR}/sv_visualization_report.html" ]]; then
        mark_success "report.visualization"
        echo -e "${GREEN}✓ Visualization report generated${NC}"
      else
        mark_failure "report.visualization_missing_output"
      fi
    else
      mark_failure "report.visualization"
    fi
  else
    echo "      Report script not found: ${WORK_DIR}/sv_visualization_report.py"
    mark_failure "report.script_missing"
  fi
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
echo "  Simulated:     ${SIM_DIR}"
echo "  Real data:     ${REAL_DIR}"
echo "  Million-real:  ${MILLION_REAL_DIR}"
echo "  Report:        ${REPORT_DIR}"
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
find "${SIM_DIR}" "${REAL_DIR}" "${MILLION_REAL_DIR}" "${REPORT_DIR}" -type f \
  \( -name "*.vcf" -o -name "*.vcf.gz" -o -name "*.tsv" -o -name "*.fasta" -o -name "*.fastq" -o -name "*.html" -o -name "*.png" \) \
  2>/dev/null | wc -l | xargs echo "  Total:"
echo ""
echo "Disk usage:"
du -sh "${SIM_DIR}" "${REAL_DIR}" "${MILLION_REAL_DIR}" "${REPORT_DIR}" 2>/dev/null | awk '{print "  " $0}'
echo ""
echo "Log files:"
find "${SIM_DIR}" "${REAL_DIR}" "${MILLION_REAL_DIR}" "${REPORT_DIR}" -name "*.log" 2>/dev/null | wc -l | xargs echo "  Total:"
echo ""
echo -e "${GREEN}✓ All experiments complete! End time: $(date)${NC}"
echo ""
echo "Next steps:"
echo "  1. Review logs for errors: grep -r 'ERROR' ${SIM_DIR} ${REAL_DIR} ${MILLION_REAL_DIR} ${REPORT_DIR}"
echo "  2. Open report: ${REPORT_DIR}/sv_visualization_report.html"
echo ""

if [[ ${#SUCCESS_STAGES[@]} -eq 0 ]]; then
  exit 1
fi
exit 0
