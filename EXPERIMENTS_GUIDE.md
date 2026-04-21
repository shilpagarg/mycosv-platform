# Running All MycoSV Experiments

Quick reference for executing all experiments from shell.

## Quick Start

```bash
cd /mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/AMF/scale

# Run ALL experiments (small-scale + large-scale + real data)
bash run_all_experiments.sh

# Run only small-scale tests
bash run_all_experiments.sh --small

# Run only large-scale/million-scale benchmarks
bash run_all_experiments.sh --large

# Run only real fungal data benchmarks
bash run_all_experiments.sh --real
```

## What Gets Executed

### 1. Small-Scale Simulated Data Tests (~10-30 minutes)
```bash
python3 -m pytest test_pipeline_features.py test_amf.py test_all_use_cases.py -v
```
- Comprehensive validation of all SV types (INS, DEL, DUP, INV, TRA)
- Tests all ecological scenarios (8 total)
- Output: `experiments/small_tests/YYYYMMDD_HHMMSS/simulated/`

### 2. Small-Scale Real Fungal Data Tests (~5-10 minutes)
```bash
python3 -m pytest test_real_fungal_benchmark.py test_new_biology_candidates.py -v
```
- Real fungal genome validation
- Breakpoint precision verification
- Output: `experiments/small_tests/YYYYMMDD_HHMMSS/real_data/`

### 3. Million-Scale Simulated Data Benchmark (~20-40 minutes)
```bash
python3 run_million_mode_query_benchmark.py \
  --out-dir experiments/large_scale/YYYYMMDD_HHMMSS/million_scale_simulated \
  --modes assembly,short-reads,long-reads \
  --n-centroids 1000000 \
  --n-genomes 8 \
  --n-reps 3 \
  --seed 42
```
- Tests with 1M-centroid catalog
- All query modes (assembly, short-reads, long-reads)
- 8 query genomes × 3 replicates
- Output: `experiments/large_scale/YYYYMMDD_HHMMSS/million_scale_simulated/`

### 4. Mode Precision-Recall Benchmarks (~15-25 minutes per mode)
```bash
for mode in assembly short-reads long-reads; do
  python3 run_mode_pr_benchmark.py \
    --modes "$mode" \
    --out-dir experiments/large_scale/YYYYMMDD_HHMMSS/mode_pr_benchmark/$mode \
    --n-refs 500 \
    --n-queries 20 \
    --seed 42
done
```
- Precision/recall metrics for each query mode
- 500 references × 20 queries per mode
- Output: `experiments/large_scale/YYYYMMDD_HHMMSS/mode_pr_benchmark/{assembly,short-reads,long-reads}/`

### 5. Real Fungal Data Benchmarks (~30-60 minutes per panel)
```bash
for panel in compact_yeast amf_large cross_phylum_hgt te_rich_pathogen two_speed_pathogen; do
  # Prepare
  python3 run_real_fungal_benchmark.py prepare \
    --out-dir experiments/real_data/YYYYMMDD_HHMMSS/$panel/prepared \
    --panels "$panel"
  
  # Benchmark all modes
  for mode in assembly short-reads long-reads; do
    python3 run_real_fungal_benchmark.py benchmark \
      --prepared-dir experiments/real_data/YYYYMMDD_HHMMSS/$panel/prepared \
      --mode "$mode" \
      --out-dir experiments/real_data/YYYYMMDD_HHMMSS/$panel/benchmark_$mode
  done
done
```
- Tests 5 distinct fungal genome panels
- Each panel has unique evolutionary characteristics
- All 3 query modes tested for each panel
- Output: `experiments/real_data/YYYYMMDD_HHMMSS/{panel}/`

**Panels covered:**
- `compact_yeast`: Model organism (S. cerevisiae)
- `amf_large`: Arbuscular mycorrhizal fungi (symbiotic)
- `cross_phylum_hgt`: Horizontal gene transfer scenarios
- `te_rich_pathogen`: Transposable element diversity
- `two_speed_pathogen`: Evolutionary rate variation

## Intermediate Files Preserved

Each experiment automatically saves:
- **VCF files**: Structural variant calls
- **TSV files**: Benchmark metrics and summaries
- **FASTA files**: Assembly/reference sequences
- **FASTQ files**: Read data
- **Log files**: Full execution logs for reproducibility

## Output Directory Structure

```
experiments/
├── small_tests/
│   └── YYYYMMDD_HHMMSS/
│       ├── simulated/
│       │   ├── pytest_output.log
│       │   ├── .pytest_tmp/
│       │   └── [VCF, TSV, FASTA, FASTQ files]
│       └── real_data/
│           ├── pytest_output.log
│           ├── .pytest_tmp/
│           └── [VCF, TSV, FASTA, FASTQ files]
├── large_scale/
│   └── YYYYMMDD_HHMMSS/
│       ├── million_scale_simulated/
│       │   ├── benchmark.log
│       │   └── [results, intermediates]
│       └── mode_pr_benchmark/
│           ├── assembly/
│           ├── short-reads/
│           └── long-reads/
└── real_data/
    └── YYYYMMDD_HHMMSS/
        ├── compact_yeast/
        ├── amf_large/
        ├── cross_phylum_hgt/
        ├── te_rich_pathogen/
        └── two_speed_pathogen/
            ├── prepared/
            ├── benchmark_assembly/
            ├── benchmark_short-reads/
            └── benchmark_long-reads/
```

## Analyzing Results

After experiments complete:

### View logs
```bash
TIMESTAMP="YYYYMMDD_HHMMSS"
# Check for errors
grep -r "ERROR" experiments/small_tests/$TIMESTAMP
grep -r "ERROR" experiments/large_scale/$TIMESTAMP
grep -r "ERROR" experiments/real_data/$TIMESTAMP

# View test output
cat experiments/small_tests/$TIMESTAMP/simulated/pytest_output.log
cat experiments/small_tests/$TIMESTAMP/real_data/pytest_output.log
```

### Analyze results
```bash
python3 analyze_results.py --input-dir experiments/large_scale/$TIMESTAMP
```

### Generate comprehensive report
```bash
bash run_comprehensive_experiments.sh
```

### Check intermediate files
```bash
TIMESTAMP="YYYYMMDD_HHMMSS"
echo "=== VCF files ==="
find experiments -name "*.vcf" | head -10

echo "=== TSV files ==="
find experiments -name "*.tsv" | head -10

echo "=== Total intermediate files ==="
find experiments -type f \( -name "*.vcf" -o -name "*.vcf.gz" -o -name "*.tsv" -o -name "*.fasta" -o -name "*.fastq" \) | wc -l

echo "=== Disk usage ==="
du -sh experiments/small_tests/$TIMESTAMP
du -sh experiments/large_scale/$TIMESTAMP
du -sh experiments/real_data/$TIMESTAMP
```

## System Requirements

- **Python 3.9+** with pytest, numpy, scipy
- **g++ 9+** (C++17 compiler)
- **Threading support** (pthread)
- **Disk space**: ~10-50 GB for all experiments
- **Time**: ~2-4 hours for complete execution

## Customizing Experiments

Edit parameters in the master script or run individual commands:

```bash
# Small configuration (quick testing)
python3 run_million_mode_query_benchmark.py \
  --n-centroids 100000 \
  --n-genomes 4 \
  --n-reps 2

# Large configuration (comprehensive)
python3 run_million_mode_query_benchmark.py \
  --n-centroids 5000000 \
  --n-genomes 16 \
  --n-reps 5

# Run specific tests only
python3 -m pytest test_pipeline_features.py::test_tandem_repeat -v

# Run with specific random seed
python3 run_real_fungal_benchmark.py prepare --seed 12345
```

## Troubleshooting

### Memory issues
```bash
# Reduce dataset size
python3 run_million_mode_query_benchmark.py --n-centroids 500000 --n-genomes 4
```

### Tests hanging
```bash
# Run with timeout
timeout 3600 bash run_all_experiments.sh --small
```

### Missing intermediates
```bash
# Check if test directories exist
ls -la experiments/
ls -la experiments/small_tests/
ls -la experiments/large_scale/
ls -la experiments/real_data/
```

### Re-run specific panel
```bash
TIMESTAMP="YYYYMMDD_HHMMSS"
python3 run_real_fungal_benchmark.py prepare \
  --out-dir "experiments/real_data/${TIMESTAMP}/amf_large/prepared" \
  --panels "amf_large"
```

## Performance Notes

- **Small-scale**: ~10-30 min (good for quick validation)
- **Large-scale**: ~40-80 min (comprehensive benchmarking)
- **Real data**: ~2-3 hours (all 5 panels with all modes)
- **Complete suite**: ~4-5 hours total

Run subset experiments during development, full suite for publication/release.
