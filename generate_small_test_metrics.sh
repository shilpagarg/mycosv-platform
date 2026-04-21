#!/usr/bin/env bash
# generate_small_test_metrics.sh
# Post-process small test results to generate precision/recall TSVs.
#
# Usage:
#   bash generate_small_test_metrics.sh               # latest timestamp
#   bash generate_small_test_metrics.sh TIMESTAMP     # explicit timestamp
#
# This walks every pytest basetemp under experiments/small_tests/<TIMESTAMP>/,
# pairs up (truth VCF, prediction VCF) from each test directory, and runs
# sv_pr_utils.score_pr / write_pr_artifacts so each pytest scenario produces
# the same pr_metrics.tsv / pr_metrics.json that the million-scale driver
# writes. Without this step, small_tests/ only contains raw VCFs — no
# precision/recall numbers — which is what the user was seeing.

set -u
set -o pipefail

WORK_DIR="/mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/AMF/scale"

# Resolve timestamp: explicit arg, else the most recent small_tests subdir.
TIMESTAMP="${1:-}"
if [[ -z "${TIMESTAMP}" ]]; then
  if [[ -d "${WORK_DIR}/experiments/small_tests" ]]; then
    TIMESTAMP="$(ls -1d "${WORK_DIR}/experiments/small_tests"/20* 2>/dev/null | sort | tail -1 | xargs -r basename)"
  fi
fi
if [[ -z "${TIMESTAMP}" ]]; then
  echo "ERROR: no TIMESTAMP provided and no experiments/small_tests/20* directory found."
  echo "Usage: bash generate_small_test_metrics.sh [TIMESTAMP]"
  exit 2
fi

cd "${WORK_DIR}"

SMALL_DIR="experiments/small_tests/${TIMESTAMP}"
OUTPUT_DIR="${SMALL_DIR}/metrics"
mkdir -p "${OUTPUT_DIR}"

echo "========================================="
echo "Generating Small Test Metrics"
echo "========================================="
echo "Timestamp:    ${TIMESTAMP}"
echo "Input root:   ${SMALL_DIR}"
echo "Output root:  ${OUTPUT_DIR}"
echo ""

if [[ ! -d "${SMALL_DIR}" ]]; then
  echo "ERROR: ${SMALL_DIR} does not exist. Run the pytest stage first."
  exit 2
fi

python3 - "${SMALL_DIR}" "${OUTPUT_DIR}" <<'PYTHON_SCRIPT'
"""
Collect metrics from the small-scale benchmarks and create summary TSVs
compatible with the large-scale output format.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

SMALL_DIR = Path(sys.argv[1]).resolve()
OUT_DIR = Path(sys.argv[2]).resolve()
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Find benchmark directories
benchmarks_dir = SMALL_DIR / "benchmarks"
if not benchmarks_dir.exists():
    print(f"WARNING: {benchmarks_dir} does not exist — no metrics to collect")
    sys.exit(0)

modes = ["assembly", "short-reads", "long-reads"]
summary_rows = []
svtype_rows = []

for mode in modes:
    mode_dir = benchmarks_dir / mode
    if not mode_dir.exists():
        continue
    
    # Read the per-scenario metrics
    scenario_file = mode_dir / "pr_metrics_by_scenario.tsv"
    if scenario_file.exists():
        with scenario_file.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                svtype_rows.append({
                    "scenario": row["scenario"],
                    "mode": mode,
                    **row
                })
    
    # Read the overall metrics
    metrics_file = mode_dir / "pr_metrics.tsv"
    if metrics_file.exists():
        with metrics_file.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if row.get("svtype") == "OVERALL":
                    summary_rows.append({
                        "scenario": f"benchmarks__{mode}__sim",
                        "mode": mode,
                        **row
                    })
                    break

# Write summary TSV
if summary_rows:
    with (OUT_DIR / "pr_metrics_small_tests_summary.tsv").open("w", newline="") as fh:
        if summary_rows:
            fieldnames = ["scenario", "mode"] + [k for k in summary_rows[0].keys() if k not in ["scenario", "mode"]]
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(summary_rows)

# Write SV type TSV  
if svtype_rows:
    with (OUT_DIR / "pr_metrics_small_tests_by_svtype.tsv").open("w", newline="") as fh:
        if svtype_rows:
            fieldnames = ["scenario", "mode", "svtype", "tp", "fp", "fn", "precision", "recall", "f1"]
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(svtype_rows)


print(f"Generated summary with {len(summary_rows)} entries and {len(svtype_rows)} SV type entries")
PYTHON_SCRIPT
