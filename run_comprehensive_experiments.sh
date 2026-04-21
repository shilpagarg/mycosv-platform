#!/usr/bin/env bash
# run_comprehensive_experiments.sh
# Designed for Linux.
# Comprehensive fungal genome experiments: million-scale + real-data benchmarks
# with accuracy/efficiency metrics for MycoSV vs. state-of-the-art SV callers.
#
# Note: the previous revision of this script used `--panels` (plural) when
# invoking run_real_fungal_benchmark.py, but the CLI actually accepts `--panel`
# (singular, repeatable via action="append"). That mismatch caused every panel
# to fail silently (no `set -e`), leaving real_data/* empty. Fixed below.

set -u
set -o pipefail

WORK_DIR="/mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/AMF/scale"
OUT_DIR="${WORK_DIR}/comprehensive_experiments"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Shared download cache: FASTA files downloaded here are reused across runs.
DATA_CACHE_DIR="${DATA_CACHE_DIR:-${WORK_DIR}/data_cache}"
mkdir -p "${DATA_CACHE_DIR}"

# Create output directories for organizing intermediate files
mkdir -p "${OUT_DIR}"
mkdir -p "${WORK_DIR}/experiments/small_tests/${TIMESTAMP}"
mkdir -p "${WORK_DIR}/experiments/large_scale/${TIMESTAMP}"
cd "${WORK_DIR}"

# Scenarios covering all 5 SV types (INS/DEL/DUP/INV/TRA).
# compact_yeast: DEL+INS  two_speed_pathogen_extreme: INV+TRA+INS  arbuscular_mf: DUP+INS
LARGE_SCENARIO_SET="compact_yeast,two_speed_pathogen_extreme,arbuscular_mf"

echo "========================================"
echo "Fungal Genome Comprehensive Experiments"
echo "========================================"
echo "Start time: $(date)"
echo "Working directory: ${WORK_DIR}"
echo "Output directory: ${OUT_DIR}"
echo ""

FAILED_STAGES=()
SUCCESS_STAGES=()
mark_success() { SUCCESS_STAGES+=("$1"); }
mark_failure() { FAILED_STAGES+=("$1"); echo "  ✗ Stage failed: $1"; }

# ============================================================================
# 1. Million-scale simulated data experiments
# ============================================================================

echo "[1/3] Running million-scale simulated data benchmarks..."
echo "      Modes: assembly, short-reads, long-reads"
echo "      Scenarios: ${LARGE_SCENARIO_SET} (all 5 SV types)"
echo "      Genomes: 1M centroids with 8 query genomes, 3 replicates"

MILLION_OUT="${OUT_DIR}/million_scale_expanded_${TIMESTAMP}"
MILLION_LOG="${OUT_DIR}/million_scale_expanded_${TIMESTAMP}.log"

if python3 run_million_mode_query_benchmark.py \
    --out-dir "${MILLION_OUT}" \
    --modes assembly,short-reads,long-reads \
    --scenario-set "${LARGE_SCENARIO_SET}" \
    --n-centroids 1000000 \
    --n-genomes 8 \
    --n-reps 3 \
    --seed 42 \
    2>&1 | tee "${MILLION_LOG}"; then
  mark_success "million_scale"
else
  mark_failure "million_scale"
fi

# Collect million-scale results
echo ""
echo "[1/3] Million-scale benchmarks complete. Results:"
MILLION_SUMMARY="${MILLION_OUT}/million_mode_summary.tsv"
if [ -f "${MILLION_SUMMARY}" ]; then
  echo "      $(cat "${MILLION_SUMMARY}")"
else
  echo "      (Summary file not found)"
fi

# ----------------------------------------------------------------------------
# [1b/3] Real million-scale fungal routing index
# ----------------------------------------------------------------------------
# Downloads up to MILLION_REAL_MAX_ASSEMBLIES real NCBI RefSeq assemblies,
# builds a MycoSV routing index over them, and pads to MILLION_REAL_TARGET_CENTROIDS
# with synthetic decoys. This is the real-data counterpart to the simulated
# million-scale index above. Skip with SKIP_MILLION_REAL=1 if you just want
# the panel benchmarks.
# ----------------------------------------------------------------------------
MILLION_REAL_MAX_ASSEMBLIES="${MILLION_REAL_MAX_ASSEMBLIES:-1000}"
MILLION_REAL_TARGET_CENTROIDS="${MILLION_REAL_TARGET_CENTROIDS:-1000000}"
MILLION_REAL_OUT="${OUT_DIR}/million_real_${TIMESTAMP}"

if [[ "${SKIP_MILLION_REAL:-0}" != "1" ]]; then
  echo ""
  echo "[1b/3] Building real million-scale fungal routing index..."
  echo "       Downloading up to ${MILLION_REAL_MAX_ASSEMBLIES} NCBI RefSeq assemblies"
  echo "       Target centroids: ${MILLION_REAL_TARGET_CENTROIDS}"
  echo "       Output: ${MILLION_REAL_OUT}"
  mkdir -p "${MILLION_REAL_OUT}"
  if python3 run_real_fungal_benchmark.py prepare-million-real \
      --out-dir "${MILLION_REAL_OUT}" \
      --source ncbi-refseq \
      --max-assemblies "${MILLION_REAL_MAX_ASSEMBLIES}" \
      --target-centroids "${MILLION_REAL_TARGET_CENTROIDS}" \
      --min-assembly-level scaffold \
      --latest-only \
      --threads 4 \
      --seed 42 \
      --data-cache-dir "${DATA_CACHE_DIR}" \
      2>&1 | tee "${MILLION_REAL_OUT}/prepare_million_real.log"; then
    mark_success "million_real"
  else
    mark_failure "million_real"
  fi
else
  echo ""
  echo "[1b/3] Skipping real million-scale index build (SKIP_MILLION_REAL=1)"
fi

echo ""
echo "[2/3] Running real fungal data benchmarks..."
echo "      Panels: compact_yeast, amf_large, cross_phylum_hgt, te_rich_pathogen, two_speed_pathogen"
echo "      Comparators by mode:"
echo "         assembly    -> SyRI + minigraph + PGGB + Minigraph-Cactus + SVIM-asm + AnchorWave"
echo "         short-reads -> Delly + Manta"
echo "         long-reads  -> SVIM + Sniffles + cuteSV"

# Real data benchmark panels to test
PANELS=("compact_yeast" "amf_large" "cross_phylum_hgt" "te_rich_pathogen" "two_speed_pathogen")

for PANEL in "${PANELS[@]}"; do
  echo ""
  echo "  ── Panel: ${PANEL} ──"
  PANEL_OUT="${OUT_DIR}/real_data_${PANEL}_${TIMESTAMP}"
  mkdir -p "${PANEL_OUT}"

  # Prepare real data. Note the CLI flag is --panel (singular, repeatable).
  # --query-mode mixed also pulls public ENA read runs so benchmark_short-reads/
  # and benchmark_long-reads/ actually produce VCFs instead of sitting empty.
  echo "    - Preparing real data (mixed: assembly + short-reads + long-reads)..."
  if python3 run_real_fungal_benchmark.py prepare \
      --out-dir "${PANEL_OUT}/prepared" \
      --panel "${PANEL}" \
      --max-assemblies-per-species 8 \
      --querys-per-species 5 \
      --max-ref-downloads 20 \
      --max-query-downloads 10 \
      --query-mode mixed \
      --read-accessions-per-species 2 \
      --allow-no-queries \
      --data-cache-dir "${DATA_CACHE_DIR}" \
      2>&1 | tee "${PANEL_OUT}/prepare_${PANEL}.log"; then
    mark_success "prepare.${PANEL}"
  else
    mark_failure "prepare.${PANEL}"
    echo "    - Skipping benchmarks for ${PANEL} (prepare failed)"
    continue
  fi

  # Require a catalog AND a query manifest before benchmarking.
  if [[ ! -f "${PANEL_OUT}/prepared/selected_catalog.tsv" \
     && ! -f "${PANEL_OUT}/prepared/reference_catalog.tsv" ]]; then
    echo "    - Skipping (no catalog produced)"
    mark_failure "prepare.${PANEL}.no_catalog"
    continue
  fi
  if [[ ! -f "${PANEL_OUT}/prepared/query_manifest.tsv" ]]; then
    echo "    - Skipping (no query manifest produced)"
    mark_failure "prepare.${PANEL}.no_queries"
    continue
  fi

  for MODE in assembly short-reads long-reads; do
    echo "    - Benchmarking mode: ${MODE}..."

    comparator_flags=()
    case "${MODE}" in
      assembly)
        # SyRI + minigraph + PGGB baseline, plus fungi/pangenome-oriented
        # additions (Minigraph-Cactus, SVIM-asm, AnchorWave) that no-op
        # gracefully when the binaries aren't installed.
        comparator_flags+=(--run-syri --run-minigraph --run-pggb)
        comparator_flags+=(--run-cactus --run-svim-asm --run-anchorwave)
        ;;
      short-reads) comparator_flags+=(--run-delly --run-manta) ;;
      long-reads)  comparator_flags+=(--run-svim --run-sniffles --run-cutesv) ;;
    esac

    if python3 run_real_fungal_benchmark.py benchmark \
        --prepared-dir "${PANEL_OUT}/prepared" \
        --mode "${MODE}" \
        --out-dir "${PANEL_OUT}/benchmark_${MODE}" \
        --threads 4 \
        "${comparator_flags[@]}" \
        2>&1 | tee "${PANEL_OUT}/benchmark_${MODE}.log"; then
      mark_success "benchmark.${PANEL}.${MODE}"
    else
      mark_failure "benchmark.${PANEL}.${MODE}"
    fi
  done
done

echo ""
echo "[3/3] Generating comprehensive summary report..."

# Compile results
REPORT="${OUT_DIR}/comprehensive_experiment_report_${TIMESTAMP}.md"
cat > "${REPORT}" << 'EOF'
# Comprehensive Fungal Genome Experiments Report

## Executive Summary

This report compiles accuracy and efficiency metrics for MycoSV across:
- **Simulated data**: 1M-centroid catalog with expanded query samples
- **Real data**: Multiple fungal panels with diverse evolutionary scenarios
- **Comparators**: SyRI, minigraph, PGGB, Minigraph-Cactus, SVIM-asm,
  AnchorWave (assembly); Delly, Manta (short-reads);
  SVIM, Sniffles, cuteSV (long-reads).

## Simulated Data Results

### Million-Scale Configuration
- Catalog size: 1,000,000 fungal genomes (centroids)
- Query samples: 8 test genomes, 3 replicates each
- Seed: 42 (reproducible)
- Scenarios: te_rich_pathogen + cross_phylum_hgt (all 5 SV types)

### Query Accuracy & Efficiency by Mode

EOF

# Append million-scale summary
if [ -f "${MILLION_SUMMARY}" ]; then
  echo "#### Table: Million-Scale Results" >> "${REPORT}"
  echo '```' >> "${REPORT}"
  cat "${MILLION_SUMMARY}" >> "${REPORT}"
  echo '```' >> "${REPORT}"
else
  echo "(Results not available — see ${MILLION_LOG})" >> "${REPORT}"
fi

# Add real data summary
{
  echo ""
  echo "## Real Fungal Data Results"
  echo ""
  echo "### Panels Tested"
  echo ""
  for PANEL in "${PANELS[@]}"; do
    PANEL_OUT="${OUT_DIR}/real_data_${PANEL}_${TIMESTAMP}"
    echo "- **${PANEL}**: \`real_data_${PANEL}_${TIMESTAMP}/\`"
    for MODE in assembly short-reads long-reads; do
      SUM="${PANEL_OUT}/benchmark_${MODE}/exact_benchmark_summary.tsv"
      if [ -f "${SUM}" ]; then
        rows=$(( $(wc -l < "${SUM}") - 1 ))
        echo "  - ${MODE}: ${rows} benchmark rows written"
      else
        echo "  - ${MODE}: (no exact_benchmark_summary.tsv — see log)"
      fi
    done
  done

  echo ""
  echo "## Summary Statistics"
  echo ""
  echo "Experiment completed at: $(date)"
  echo ""
  echo "### Stage Outcomes"
  echo "- Succeeded: ${#SUCCESS_STAGES[@]}"
  echo "- Failed:    ${#FAILED_STAGES[@]}"
  if [[ ${#FAILED_STAGES[@]} -gt 0 ]]; then
    echo ""
    echo "### Failed Stages"
    for s in "${FAILED_STAGES[@]}"; do echo "- ${s}"; done
  fi
  echo ""
  echo "### Key Metrics"
  echo "- **Million-scale**: See summary above"
  echo "- **Real data panels**: Check individual benchmark directories"
} >> "${REPORT}"

echo ""
echo "Report saved to: ${REPORT}"
echo ""
echo "========================================"
echo "Experiments complete!"
echo "End time: $(date)"
echo "Succeeded: ${#SUCCESS_STAGES[@]}  Failed: ${#FAILED_STAGES[@]}"
echo "========================================"

if [[ ${#SUCCESS_STAGES[@]} -eq 0 ]]; then
  exit 1
fi
exit 0
