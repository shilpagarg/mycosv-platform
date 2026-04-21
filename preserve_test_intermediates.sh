#!/usr/bin/env bash
# preserve_test_intermediates.sh
# Wrapper to preserve VCF, TSV, FASTA, FASTQ intermediates from test runs
# Usage: bash preserve_test_intermediates.sh [small|large|all]

set -e

WORK_DIR="/mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/AMF/scale"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
TEST_TYPE="${1:-small}"  # small, large, or all

case "$TEST_TYPE" in
  small)
    TARGET_DIR="${WORK_DIR}/experiments/small_tests/${TIMESTAMP}"
    echo "Running small-scale tests with intermediate preservation..."
    echo "Output: ${TARGET_DIR}"
    ;;
  large)
    TARGET_DIR="${WORK_DIR}/experiments/large_scale/${TIMESTAMP}"
    echo "Running large-scale tests with intermediate preservation..."
    echo "Output: ${TARGET_DIR}"
    ;;
  all)
    TARGET_DIR="${WORK_DIR}/experiments/all_tests/${TIMESTAMP}"
    echo "Running all-scale tests with intermediate preservation..."
    echo "Output: ${TARGET_DIR}"
    ;;
  *)
    echo "Usage: $0 [small|large|all]"
    exit 1
    ;;
esac

mkdir -p "${TARGET_DIR}"

# Set pytest output/cache options to preserve test artifacts
export PYTEST_CACHE_DIR="${TARGET_DIR}/.pytest_cache"
export TMPDIR="${TARGET_DIR}/tmp"
mkdir -p "${TMPDIR}"

echo ""
echo "Starting test run at $(date)..."
echo "Test intermediates will be preserved in: ${TARGET_DIR}"
echo ""

# Run tests and capture output
cd "${WORK_DIR}"

# Run pytest with basetemp to preserve temp directory
python3 -m pytest \
  --basetemp="${TARGET_DIR}/pytest_tmp" \
  -v \
  test_amf.py \
  test_pipeline_features.py \
  test_all_use_cases.py \
  2>&1 | tee "${TARGET_DIR}/test_run_${TIMESTAMP}.log"

# Organize the preserved files
echo ""
echo "Test run complete. Organizing intermediate files..."
echo ""

# Check what was generated
if [[ -d "${TARGET_DIR}/pytest_tmp" ]]; then
  echo "Preserved pytest temporary files:"
  find "${TARGET_DIR}/pytest_tmp" -type f \( -name "*.vcf" -o -name "*.tsv" -o -name "*.fasta" -o -name "*.fastq" -o -name "*.fq" \) | head -20
  echo ""
  echo "Total preserved files: $(find "${TARGET_DIR}/pytest_tmp" -type f | wc -l)"
  
  # Create summary
  echo ""
  echo "Summary of preserved intermediate files:"
  echo "=========================================="
  find "${TARGET_DIR}/pytest_tmp" -type f \( -name "*.vcf" -o -name "*.tsv" -o -name "*.fasta" -o -name "*.fastq" -o -name "*.fq" \) | wc -l | xargs echo "  - VCF/TSV/FASTA/FASTQ files:"
  find "${TARGET_DIR}/pytest_tmp" -type d | wc -l | xargs echo "  - Directories:"
  
  du -sh "${TARGET_DIR}" | awk '{print "  - Total size: " $1}'
fi

echo ""
echo "Test run completed at $(date)"
echo "Intermediate files location: ${TARGET_DIR}"
