#!/usr/bin/env bash
# quick_start_intermediates.sh
# Quick reference: One-liner commands to preserve test outputs

cat << 'EOF'
╔════════════════════════════════════════════════════════════════════════════╗
║                  PRESERVE TEST INTERMEDIATES - QUICK START                 ║
╚════════════════════════════════════════════════════════════════════════════╝

OVERVIEW
────────
Intermediate files (VCF, TSV, FASTA, FASTQ) are captured automatically.
Use these commands to preserve outputs from simulated benchmarks.

QUICK COMMANDS
──────────────

1. PRESERVE TEST INTERMEDIATES (unit tests)
   ──────────────────────────────────────────
   bash preserve_test_intermediates.sh

   Output: experiments/simulated/YYYYMMDD_HHMMSS/
   Includes: VCF truth calls, TSV metadata, FASTA refs, FASTQ queries

2. CAPTURE BENCHMARK INTERMEDIATES (all query modes)
   ───────────────────────────────────────────────────
   bash capture_benchmark_intermediates.sh all

   Output: experiments/simulated/YYYYMMDD_HHMMSS/{assembly,short-reads,long-reads}/
   Includes: Precision/recall TSV, VCF calls, benchmark logs

3. SINGLE QUERY MODE
   ──────────────────
   bash capture_benchmark_intermediates.sh assembly
   bash capture_benchmark_intermediates.sh short-reads
   bash capture_benchmark_intermediates.sh long-reads

   Output: experiments/simulated/YYYYMMDD_HHMMSS/{mode}/
   Includes: Mode-specific VCF/TSV outputs

EXAMPLE WORKFLOW
────────────────

# Step 1: Run the full simulated benchmark
bash run_all_experiments.sh --simulated

# Step 2: Find the output directory
ls -d experiments/simulated/20* | tail -1

# Step 3: Inspect VCF truth calls
cat experiments/simulated/20*/benchmarks/*/truth/*.vcf | head -20

# Step 4: View precision/recall metrics
ls experiments/simulated/20*/benchmarks/

# Step 5: Debug with raw VCF files
vcftools --vcf experiments/simulated/20*/benchmarks/*/truth/*.vcf --freq

DIRECTORY STRUCTURE
───────────────────

experiments/
  simulated/
    YYYYMMDD_HHMMSS/
      benchmarks/
        assembly/
          pr_metrics.tsv        ← Overall precision/recall
          pr_metrics_by_scenario.tsv
        short-reads/
          pr_metrics.tsv
        long-reads/
          pr_metrics.tsv
      metrics/
        pr_metrics_simulated_summary.tsv
        pr_metrics_simulated_by_svtype.tsv

  real_data/
    YYYYMMDD_HHMMSS/
      {panel}/
        prepared/               ← Reference catalog + query manifest
        benchmark_{mode}/       ← Mode-specific VCF/TSV outputs

  million_real/
    YYYYMMDD_HHMMSS/
      routing_index/            ← Real NCBI-based index

  test_results/
    *.md                        ← Final analysis reports
    *.log                       ← Test logs

KEY FILE TYPES
──────────────

VCF    - Variant Call Format (truth + calls)
TSV    - Tab-separated: metadata, metrics, manifests
FASTA  - Reference assemblies
FASTQ  - Query reads (short or long)
LOG    - Execution logs

SIZE MANAGEMENT
───────────────

Check disk usage:
  du -sh experiments/simulated/ experiments/real_data/

Remove old runs (>7 days):
  find experiments/simulated -maxdepth 1 -type d -mtime +7 -exec rm -rf {} \;

Compress old runs:
  tar -czf experiments/backup_20240420.tar.gz experiments/simulated/YYYYMMDD*

COMMON ISSUES
─────────────

❌ "Intermediates not appearing?"
   → Check log: tail experiments/simulated/YYYYMMDD_HHMMSS/benchmarks/benchmark.log
   → Verify test passed: grep FAILED benchmarks/benchmark.log
   → Check disk: df -h

❌ "Disk space exhausted?"
   → Remove old runs: find experiments -mtime +30 -type d -exec rm -rf {} \;
   → Compress and move: tar -czf backup.tar.gz experiments/ && rm -rf experiments/

❌ "Want to run specific test only?"
   → python3 -m pytest test_amf.py::test_name --basetemp=./experiments/simulated/custom

DOCUMENTATION
──────────────

Full guide: INTERMEDIATE_FILES_GUIDE.md
This file: $(basename $0)

═════════════════════════════════════════════════════════════════════════════
EOF
