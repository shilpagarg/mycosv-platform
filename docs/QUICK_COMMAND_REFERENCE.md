# Quick Command Reference

Set `MYCOSV_ROOT` to your clone before running any command below:

```bash
export MYCOSV_ROOT=/path/to/mycosv-platform
```

## One-Liner Commands to Run All Experiments

### Complete Suite (All experiments with all outputs)
```bash
cd $MYCOSV_ROOT && bash run_all_experiments.sh
```

### Small-Scale Only (Quick validation, ~30 min)
```bash
cd $MYCOSV_ROOT && bash run_all_experiments.sh --small
```

### Large-Scale Only (Million-scale benchmarks, ~90 min)
```bash
cd $MYCOSV_ROOT && bash run_all_experiments.sh --large
```

### Real Data Only (Fungal genome benchmarks, ~2 hours)
```bash
cd $MYCOSV_ROOT && bash run_all_experiments.sh --real
```

---

## Individual Experiment Commands

### Small-Scale Simulated Tests Only
```bash
cd $MYCOSV_ROOT
mkdir -p experiments/small_tests/$(date +%Y%m%d_%H%M%S)
python3 -m pytest test_pipeline_features.py test_amf.py test_all_use_cases.py -v --tb=short
```

### Small-Scale Real Data Tests Only
```bash
cd $MYCOSV_ROOT
python3 -m pytest test_real_fungal_benchmark.py test_new_biology_candidates.py -v --tb=short
```

### Million-Scale Benchmark Only (All modes)
```bash
cd $MYCOSV_ROOT
mkdir -p experiments/large_scale/$(date +%Y%m%d_%H%M%S)
python3 run_million_mode_query_benchmark.py \
  --out-dir experiments/large_scale/$(date +%Y%m%d_%H%M%S)/million_scale \
  --modes assembly,short-reads,long-reads \
  --n-centroids 1000000 \
  --n-genomes 8 \
  --n-reps 3 \
  --seed 42
```

### Million-Scale Per-Mode Benchmark
```bash
cd $MYCOSV_ROOT
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

for mode in assembly short-reads long-reads; do
  python3 run_million_mode_query_benchmark.py \
    --modes "$mode" \
    --out-dir experiments/large_scale/$TIMESTAMP/mode_$mode \
    --n-centroids 1000000 \
    --n-genomes 8 \
    --seed 42
done
```

### Visualization Report (HTML + PNG)
```bash
cd $MYCOSV_ROOT
python3 sv_visualization_report.py \
  --sv-tsv experiments/large_scale/$TIMESTAMP/results.sv.tsv \
  --out-dir reports/$TIMESTAMP
# Output: reports/$TIMESTAMP/mycosv_report.html
#   Section 1: SV type breakdown
#   Section 2: Precision/recall per mode
#   Section 3: TE classification
#   Section 4: Clade-SV, TE-architecture, HGT-propagation plots
```

### Real Fungal Data - Single Panel (e.g., compact_yeast)
```bash
cd $MYCOSV_ROOT
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
PANEL="compact_yeast"
mkdir -p experiments/real_data/$TIMESTAMP/$PANEL/{prepared,benchmark_assembly,benchmark_short-reads,benchmark_long-reads}

# Prepare panel
python3 run_real_fungal_benchmark.py prepare \
  --out-dir experiments/real_data/$TIMESTAMP/$PANEL/prepared \
  --panel $PANEL

# Benchmark all modes
python3 run_real_fungal_benchmark.py benchmark \
  --prepared-dir experiments/real_data/$TIMESTAMP/$PANEL/prepared \
  --mode assembly \
  --out-dir experiments/real_data/$TIMESTAMP/$PANEL/benchmark_assembly

python3 run_real_fungal_benchmark.py benchmark \
  --prepared-dir experiments/real_data/$TIMESTAMP/$PANEL/prepared \
  --mode short-reads \
  --out-dir experiments/real_data/$TIMESTAMP/$PANEL/benchmark_short-reads

python3 run_real_fungal_benchmark.py benchmark \
  --prepared-dir experiments/real_data/$TIMESTAMP/$PANEL/prepared \
  --mode long-reads \
  --out-dir experiments/real_data/$TIMESTAMP/$PANEL/benchmark_long-reads
```

### Real Fungal Data - Expanded Default Panels
```bash
cd $MYCOSV_ROOT
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
mkdir -p experiments/real_data/$TIMESTAMP

for PANEL in compact_yeast amf_large te_rich_pathogen two_speed_pathogen; do
  echo "Processing: $PANEL"
  mkdir -p experiments/real_data/$TIMESTAMP/$PANEL/{prepared,benchmark_assembly,benchmark_short-reads,benchmark_long-reads}
  
  python3 run_real_fungal_benchmark.py prepare \
    --out-dir experiments/real_data/$TIMESTAMP/$PANEL/prepared \
    --panel $PANEL \
    --query-mode mixed \
    --read-accessions-per-species ${REAL_READ_ACCESSIONS_PER_SPECIES:-1}
  
  if [ -f "experiments/real_data/$TIMESTAMP/$PANEL/prepared/reference_catalog.tsv" ]; then
    for MODE in assembly short-reads long-reads; do
      python3 run_real_fungal_benchmark.py benchmark \
        --prepared-dir experiments/real_data/$TIMESTAMP/$PANEL/prepared \
        --mode $MODE \
        --out-dir experiments/real_data/$TIMESTAMP/$PANEL/benchmark_$MODE
    done
  fi
done
```

---

## Checking Results

### List all experiment outputs
```bash
ls -lah $MYCOSV_ROOT/experiments/
```

### See intermediate files for a specific experiment
```bash
TIMESTAMP="YYYYMMDD_HHMMSS"  # Replace with actual timestamp
find $MYCOSV_ROOT/experiments -type f \( -name "*.vcf" -o -name "*.tsv" -o -name "*.fasta" -o -name "*.fastq" \) | head -20
```

### Check disk usage
```bash
du -sh $MYCOSV_ROOT/experiments/small_tests
du -sh $MYCOSV_ROOT/experiments/large_scale
du -sh $MYCOSV_ROOT/experiments/real_data
```

### View test logs
```bash
TIMESTAMP="YYYYMMDD_HHMMSS"
cat $MYCOSV_ROOT/experiments/small_tests/$TIMESTAMP/simulated/pytest_output.log
cat $MYCOSV_ROOT/experiments/small_tests/$TIMESTAMP/real_data/pytest_output.log
```

### Check for errors in all logs
```bash
grep -r "ERROR\|FAIL\|Exception" $MYCOSV_ROOT/experiments/ 2>/dev/null | head -20
```

### Count total intermediate files saved
```bash
find $MYCOSV_ROOT/experiments -type f \( -name "*.vcf" -o -name "*.vcf.gz" -o -name "*.tsv" -o -name "*.fasta" -o -name "*.fastq" \) | wc -l
```

---

## Environment Setup

### Activate virtual environment (if needed)
```bash
source $MYCOSV_ROOT/.venv/bin/activate
```

### Check Python version
```bash
python3 --version
```

### Verify required packages
```bash
python3 -c "import pytest; import numpy; import scipy; import pandas; print('✓ All dependencies installed')"
```

### Install missing packages
```bash
pip install pytest numpy scipy pandas
```

---

## Monitoring Progress

### Watch experiment in progress
```bash
# In one terminal
cd $MYCOSV_ROOT && bash run_all_experiments.sh

# In another terminal, monitor
watch -n 5 "du -sh $MYCOSV_ROOT/experiments"
```

### See test progress
```bash
# Watch for new log files being created
watch -n 5 "find $MYCOSV_ROOT/experiments -name '*.log' | wc -l"
```

### Monitor resource usage
```bash
# CPU and memory
top -u $USER

# Disk I/O
iostat -x 1 5

# File count growth
watch -n 5 "find $MYCOSV_ROOT/experiments -type f | wc -l"
```

---

## Estimated Runtimes

| Experiment | Time | Output |
|-----------|------|--------|
| Small simulated | 10-30 min | ~100 MB |
| Small real data | 5-10 min | ~50 MB |
| Million-scale | 20-40 min | ~500 MB |
| Mode PR benchmarks | 15-25 min each | ~300 MB per mode |
| Real data panel | 30-60 min | ~200-400 MB per panel |
| **Complete suite** | **2-4 hours** | **2-3 GB** |
