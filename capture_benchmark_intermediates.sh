#!/usr/bin/env bash
# capture_benchmark_intermediates.sh
# Runs benchmarks and organizes intermediate files by query mode
# Usage: bash capture_benchmark_intermediates.sh [assembly|short-reads|long-reads|all]

set -e

WORK_DIR="/mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/AMF/scale"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
QUERY_MODE="${1:-all}"

CAPTURE_DIR="${WORK_DIR}/experiments/simulated/${TIMESTAMP}"
mkdir -p "${CAPTURE_DIR}"

echo "========================================"
echo "Capturing Benchmark Intermediates"
echo "========================================"
echo "Timestamp: ${TIMESTAMP}"
echo "Query mode: ${QUERY_MODE}"
echo "Capture directory: ${CAPTURE_DIR}"
echo ""

cd "${WORK_DIR}"

# Helper function to run a benchmark and capture outputs
run_benchmark_with_capture() {
  local mode=$1
  local mode_dir="${CAPTURE_DIR}/${mode}"
  mkdir -p "${mode_dir}"
  
  echo "[$(date +%H:%M:%S)] Running ${mode} mode benchmark..."
  
  # Run benchmark and capture VCF/TSV outputs
  python3 run_mode_pr_benchmark.py \
    --modes "${mode}" \
    --out-dir "${mode_dir}/output" \
    --n-refs 500 \
    --n-queries 20 \
    --seed 42 \
    2>&1 | tee "${mode_dir}/benchmark_${mode}.log"
  
  # Collect intermediate files
  echo ""
  echo "Collecting VCF and TSV files from ${mode}..."
  
  if [[ -d "${mode_dir}/output" ]]; then
    # Count and list VCF files
    vcf_count=$(find "${mode_dir}/output" -name "*.vcf" -o -name "*.vcf.gz" | wc -l)
    if [[ $vcf_count -gt 0 ]]; then
      echo "  ✓ Found ${vcf_count} VCF files"
      find "${mode_dir}/output" -name "*.vcf" -o -name "*.vcf.gz" | head -5
    fi
    
    # Count and list TSV files
    tsv_count=$(find "${mode_dir}/output" -name "*.tsv" | wc -l)
    if [[ $tsv_count -gt 0 ]]; then
      echo "  ✓ Found ${tsv_count} TSV files"
      find "${mode_dir}/output" -name "*.tsv" | head -5
    fi
    
    # Count and list FASTA/FASTQ
    seq_count=$(find "${mode_dir}/output" \( -name "*.fasta" -o -name "*.fastq" -o -name "*.fa" -o -name "*.fq" \) | wc -l)
    if [[ $seq_count -gt 0 ]]; then
      echo "  ✓ Found ${seq_count} sequence files"
    fi
    
    # Show directory structure
    echo ""
    echo "  Directory structure:"
    du -sh "${mode_dir}/output" | awk '{print "    Total size: " $1}'
  fi
  
  echo ""
}

# Run based on query mode
case "$QUERY_MODE" in
  assembly)
    run_benchmark_with_capture "assembly"
    ;;
  short-reads)
    run_benchmark_with_capture "short-reads"
    ;;
  long-reads)
    run_benchmark_with_capture "long-reads"
    ;;
  all)
    for mode in assembly short-reads long-reads; do
      run_benchmark_with_capture "$mode"
    done
    ;;
  *)
    echo "Usage: $0 [assembly|short-reads|long-reads|all]"
    exit 1
    ;;
esac

# Generate summary
echo ""
echo "========================================"
echo "Benchmark Capture Complete"
echo "========================================"
echo ""
echo "Summary of captured intermediates:"
echo ""

for mode_dir in "${CAPTURE_DIR}"/*; do
  if [[ -d "$mode_dir" && "$(basename "$mode_dir")" != .pytest_cache ]]; then
    mode=$(basename "$mode_dir")
    echo "  ${mode}:"
    
    vcf_count=$(find "$mode_dir" -name "*.vcf" -o -name "*.vcf.gz" 2>/dev/null | wc -l)
    tsv_count=$(find "$mode_dir" -name "*.tsv" 2>/dev/null | wc -l)
    size=$(du -sh "$mode_dir" 2>/dev/null | awk '{print $1}')
    
    [[ $vcf_count -gt 0 ]] && echo "    - VCF files: $vcf_count"
    [[ $tsv_count -gt 0 ]] && echo "    - TSV files: $tsv_count"
    echo "    - Total size: ${size}"
  fi
done

echo ""
echo "All intermediates: ${CAPTURE_DIR}"
