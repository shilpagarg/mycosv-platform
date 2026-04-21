#!/usr/bin/env bash
# quick_start_intermediates.sh
# Quick reference: One-liner commands to preserve test outputs

cat << 'EOF'
╔════════════════════════════════════════════════════════════════════════════╗
║                  PRESERVE TEST INTERMEDIATES - QUICK START                 ║
╚════════════════════════════════════════════════════════════════════════════╝

OVERVIEW
────────
Intermediate files (VCF, TSV, FASTA, FASTQ) are now captured automatically.
Use these commands to preserve outputs by test type and query mode.

QUICK COMMANDS
──────────────

1. SMALL-SCALE TESTS (simulated data, all modes)
   ────────────────────────────────────────────────
   bash preserve_test_intermediates.sh small

   Output: experiments/small_tests/YYYYMMDD_HHMMSS/
   Includes: VCF truth calls, TSV metadata, FASTA refs, FASTQ queries

2. LARGE-SCALE BENCHMARKS (all query modes)
   ─────────────────────────────────────────
   bash capture_benchmark_intermediates.sh all

   Output: experiments/large_scale/YYYYMMDD_HHMMSS/{assembly,short-reads,long-reads}/
   Includes: Precision/recall TSV, VCF calls, benchmark logs

3. SINGLE QUERY MODE
   ──────────────────
   bash capture_benchmark_intermediates.sh assembly
   bash capture_benchmark_intermediates.sh short-reads
   bash capture_benchmark_intermediates.sh long-reads

   Output: experiments/large_scale/YYYYMMDD_HHMMSS/{mode}/
   Includes: Mode-specific VCF/TSV outputs

EXAMPLE WORKFLOW
────────────────

# Step 1: Capture small-scale test intermediates
bash preserve_test_intermediates.sh small

# Step 2: Find the output directory
ls -d experiments/small_tests/20* | tail -1

# Step 3: Inspect VCF truth calls
cat experiments/small_tests/20*/pytest_tmp/test_*/truth/*.vcf | head -20

# Step 4: View precision/recall metrics
ls experiments/small_tests/20*/pytest_tmp/*/

# Step 5: Debug with raw VCF files
vcftools --vcf experiments/small_tests/20*/pytest_tmp/*/truth/*.vcf --freq

DIRECTORY STRUCTURE
───────────────────

experiments/
  small_tests/
    YYYYMMDD_HHMMSS/
      pytest_tmp/           ← All test outputs
        test_amf0/
          query_truth.tsv   ← Truth calls per query
          truth/*.vcf       ← VCF format truth
          *.fasta           ← Refs
          *.fastq           ← Queries
  
  large_scale/
    YYYYMMDD_HHMMSS/
      assembly/
        output/             ← Assembly mode outputs
          *.vcf             ← VCF calls
          *.tsv             ← Precision/recall
      short-reads/
        output/             ← SR mode outputs
      long-reads/
        output/             ← LR mode outputs

  test_results/
    *.md                    ← Final analysis reports
    *.log                   ← Test logs

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
  du -sh experiments/small_tests/ experiments/large_scale/

Remove old runs (>7 days):
  find experiments/small_tests -maxdepth 1 -type d -mtime +7 -exec rm -rf {} \;

Compress old runs:
  tar -czf experiments/backup_20240420.tar.gz experiments/small_tests/YYYYMMDD*

COMMON ISSUES
─────────────

❌ "Intermediates not appearing?"
   → Check log: tail experiments/small_tests/YYYYMMDD_HHMMSS/test_run_*.log
   → Verify test passed: grep FAILED test_run_*.log
   → Check disk: df -h

❌ "Disk space exhausted?"
   → Remove old runs: find experiments -mtime +30 -type d -exec rm -rf {} \;
   → Compress and move: tar -czf backup.tar.gz experiments/ && rm -rf experiments/

❌ "Want to run specific test only?"
   → python3 -m pytest test_amf.py::test_name --basetemp=./experiments/small_tests/custom

DOCUMENTATION
──────────────

Full guide: INTERMEDIATE_FILES_GUIDE.md
This file: $(basename $0)

═════════════════════════════════════════════════════════════════════════════
EOF
