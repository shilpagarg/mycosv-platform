#!/usr/bin/env bash
# run_quick_tests.sh
# Quick test suite without real data downloads (for local validation)
# Useful for testing bugs in small environment

set -e

WORK_DIR="/mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/AMF/scale"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Create test directory
TEST_DIR="${WORK_DIR}/experiments/quick_tests/${TIMESTAMP}"
mkdir -p "${TEST_DIR}"

cd "${WORK_DIR}"

echo "========================================"
echo "MycoSV Quick Test Suite"
echo "========================================"
echo "Start time: $(date)"
echo "Timestamp: ${TIMESTAMP}"
echo "Test directory: ${TEST_DIR}"
echo ""

# ============================================================================
# 1. SIMULATED DATA TESTS (All SV types coverage)
# ============================================================================

echo "[1/2] Running simulated data tests (all SV types)..."
mkdir -p "${TEST_DIR}/simulated"

python3 -m pytest \
  test_pipeline_features.py::test_assembly_mode_recovers_inv_dup_and_tra \
  test_pipeline_features.py::test_reads_modes_recover_inv_dup_and_tra \
  test_pipeline_features.py::test_reads_modes_recover_insertion_and_offref_end_to_end \
  -v \
  --tb=short \
  --basetemp="${TEST_DIR}/simulated/.pytest_tmp" \
  2>&1 | tee "${TEST_DIR}/simulated/pytest_output.log"

echo ""
echo "✓ Simulated data tests complete"
echo ""

# ============================================================================
# 2. MILLION-SCALE BENCHMARKS (All SV types coverage)
# ============================================================================

echo "[2/2] Running million-scale benchmarks (all 5 SV types)..."
mkdir -p "${TEST_DIR}/million_all_types"

# Test with te_rich_pathogen scenario (generates all 5 SV types)
python3 run_million_mode_query_benchmark.py \
  --out-dir "${TEST_DIR}/million_all_types" \
  --modes assembly \
  --scenario-set "te_rich_pathogen" \
  --n-centroids 100000 \
  --n-genomes 4 \
  --n-reps 2 \
  --seed 42 \
  2>&1 | tee "${TEST_DIR}/million_all_types/benchmark.log"

echo ""
echo "✓ Million-scale benchmark complete"
echo ""

# ============================================================================
# SUMMARY
# ============================================================================

echo "========================================"
echo "Quick Test Summary"
echo "========================================" 
echo "Results directory: ${TEST_DIR}"
echo ""
echo "Files generated:"
find "${TEST_DIR}" -type f \( -name "*.vcf" -o -name "*.tsv" -o -name "*.log" \) | wc -l | xargs echo "  Total:"
echo ""
echo "Test status:"
if grep -q "failed\|ERROR" "${TEST_DIR}/simulated/pytest_output.log"; then
  echo "  ✗ Simulated tests: FAILED (check log)"
else
  echo "  ✓ Simulated tests: PASSED"
fi

if grep -q "ERROR\|Traceback" "${TEST_DIR}/million_all_types/benchmark.log"; then
  echo "  ✗ Benchmark: FAILED (check log)"
else
  echo "  ✓ Benchmark: PASSED"
fi

echo ""
echo "✓ Quick tests complete! End time: $(date)"
echo ""
echo "Next steps:"
echo "  1. Review test logs: cat ${TEST_DIR}/simulated/pytest_output.log"
echo "  2. Check benchmark results: grep -E 'svtype|precision|recall|f1' ${TEST_DIR}/million_all_types/benchmark.log"
echo "  3. Verify all 5 SV types present: cut -f2 ${TEST_DIR}/million_all_types/*/pr_metrics_by_scenario.tsv | sort -u"
echo "  4. Run full suite: bash run_all_experiments.sh"
echo ""
