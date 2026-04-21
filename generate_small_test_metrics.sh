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
Walk every pytest basetemp under the given small-test directory, find pairs
of (truth VCF, prediction VCF), and write per-scenario pr_metrics.tsv /
pr_metrics.json using sv_pr_utils.score_pr / write_pr_artifacts.

Also emit a rolled-up summary TSV (pr_metrics_small_tests_summary.tsv) so
results for the small-tests stage are comparable to large_scale/.
"""
from __future__ import annotations

import csv
import json
import sys
import traceback
from pathlib import Path

# Make the repo importable regardless of cwd.
REPO = Path.cwd().resolve()
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

try:
    from sv_pr_utils import score_pr, write_pr_artifacts
except Exception as exc:  # pragma: no cover - defensive
    print(f"ERROR: failed to import sv_pr_utils: {exc}")
    sys.exit(3)

SMALL_DIR = Path(sys.argv[1]).resolve()
OUT_DIR = Path(sys.argv[2]).resolve()
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------

def find_truth_vcfs(root: Path) -> list[Path]:
    """Every test_amf-style scenario writes truth/all_queries.truth.ref.vcf."""
    return sorted(root.rglob("truth/all_queries.truth.ref.vcf"))


def find_pred_vcf_for(truth_vcf: Path) -> Path | None:
    """
    Given a truth VCF at <scenario_dir>/truth/all_queries.truth.ref.vcf,
    prefer a sibling pipeline output VCF. Search order:
      1. <scenario_dir>/calls.vcf         (pipeline convention)
      2. <scenario_dir>/all.vcf
      3. <scenario_dir>/**/calls.vcf      (nested cli output)
      4. <scenario_dir>/**/*.vcf          (first non-truth VCF found)
    """
    scenario_dir = truth_vcf.parent.parent
    for name in ("calls.vcf", "all.vcf"):
        cand = scenario_dir / name
        if cand.exists():
            return cand
    for cand in sorted(scenario_dir.rglob("calls.vcf")):
        return cand
    for cand in sorted(scenario_dir.rglob("*.vcf")):
        # Skip the truth tree to avoid truth==pred.
        if "truth" in cand.parts:
            continue
        return cand
    return None


def scenario_name_for(truth_vcf: Path) -> str:
    """Use the scenario directory path relative to SMALL_DIR as the label."""
    scenario_dir = truth_vcf.parent.parent
    try:
        rel = scenario_dir.relative_to(SMALL_DIR)
    except ValueError:
        rel = scenario_dir
    return str(rel).replace("/", "__").replace(" ", "_")


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

truth_vcfs = find_truth_vcfs(SMALL_DIR)
print(f"Found {len(truth_vcfs)} truth VCF(s) under {SMALL_DIR}")
if not truth_vcfs:
    print("No truth VCFs found — small tests may not have produced simulations.")
    print("  Hint: ensure the pytest stage ran (test_amf / test_pipeline_features).")
    # Still create an empty summary so downstream tooling sees a file.
    (OUT_DIR / "pr_metrics_small_tests_summary.tsv").write_text(
        "scenario\tmode\ttruth_records\tpred_records_algo\ttp\tfp\tfn\tprecision\trecall\tf1\tnotes\n",
        encoding="utf-8",
    )
    sys.exit(0)

summary_rows: list[dict[str, object]] = []
failures: list[str] = []

for truth_vcf in truth_vcfs:
    pred_vcf = find_pred_vcf_for(truth_vcf)
    scenario = scenario_name_for(truth_vcf)
    scen_out = OUT_DIR / scenario
    scen_out.mkdir(parents=True, exist_ok=True)

    if pred_vcf is None:
        note = "no_prediction_vcf_found"
        print(f"  [skip] {scenario}: {note} (truth={truth_vcf})")
        summary_rows.append({
            "scenario": scenario, "mode": "?", "truth_records": 0,
            "pred_records_algo": 0, "tp": 0, "fp": 0, "fn": 0,
            "precision": 0.0, "recall": 0.0, "f1": 0.0, "notes": note,
        })
        continue

    # Optional sidecars if the pipeline left them alongside the prediction.
    scenario_dir = truth_vcf.parent.parent
    hits_tsv = None
    for cand in (scenario_dir / "calls.hits.tsv", scenario_dir / "all.hits.tsv"):
        if cand.exists():
            hits_tsv = cand
            break
    if hits_tsv is None:
        for cand in sorted(scenario_dir.rglob("*.hits.tsv")):
            hits_tsv = cand
            break
    meta_tsv = scenario_dir / "query_metadata.tsv"
    if not meta_tsv.exists():
        # Look one level up (sim dir sometimes sits beside the scenario dir).
        for cand in sorted(scenario_dir.rglob("query_metadata.tsv")):
            meta_tsv = cand
            break

    # Pull --query-mode from the metadata if we can.
    mode = "?"
    if meta_tsv.exists():
        try:
            with meta_tsv.open() as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                row = next(reader, None)
                if row and row.get("query_mode"):
                    mode = row["query_mode"]
        except Exception:
            pass

    try:
        summary = score_pr(
            truth_vcf,
            pred_vcf,
            pred_hits_tsv=hits_tsv if hits_tsv and hits_tsv.exists() else None,
            query_meta_tsv=meta_tsv if meta_tsv and meta_tsv.exists() else None,
        )
        write_pr_artifacts(
            summary,
            scen_out / "pr_metrics.tsv",
            scen_out / "pr_metrics.json",
        )
        overall = summary.get("overall", {})
        tp = int(overall.get("tp", 0))
        fp = int(overall.get("fp", 0))
        fn = int(overall.get("fn", 0))
        prec = float(overall.get("precision", 0.0) or 0.0)
        rec = float(overall.get("recall", 0.0) or 0.0)
        f1 = float(overall.get("f1", 0.0) or 0.0)
        summary_rows.append({
            "scenario": scenario,
            "mode": mode,
            "truth_records": int(summary.get("truth_records", 0)),
            "pred_records_algo": int(summary.get("pred_records_algo", 0)),
            "tp": tp, "fp": fp, "fn": fn,
            "precision": round(prec, 6),
            "recall": round(rec, 6),
            "f1": round(f1, 6),
            "notes": "",
        })
        print(
            f"  [ok]   {scenario}: truth={truth_vcf.name} pred={pred_vcf.name} "
            f"TP={tp} FP={fp} FN={fn} P={prec:.3f} R={rec:.3f} F1={f1:.3f}"
        )
    except Exception as exc:  # keep going on a single failure
        err = f"{type(exc).__name__}: {exc}"
        failures.append(f"{scenario}: {err}")
        print(f"  [fail] {scenario}: {err}")
        traceback.print_exc(limit=1)
        summary_rows.append({
            "scenario": scenario, "mode": mode, "truth_records": 0,
            "pred_records_algo": 0, "tp": 0, "fp": 0, "fn": 0,
            "precision": 0.0, "recall": 0.0, "f1": 0.0,
            "notes": err[:200],
        })

# ---------------------------------------------------------------------------
# Rolled-up summary
# ---------------------------------------------------------------------------

summary_tsv = OUT_DIR / "pr_metrics_small_tests_summary.tsv"
fieldnames = [
    "scenario", "mode", "truth_records", "pred_records_algo",
    "tp", "fp", "fn", "precision", "recall", "f1", "notes",
]
with summary_tsv.open("w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    for row in summary_rows:
        writer.writerow({k: row.get(k, "") for k in fieldnames})

# Per-SV-type roll-up (sum tp/fp/fn across scenarios via the per-type JSONs).
per_type_counts: dict[str, dict[str, int]] = {}
for truth_vcf in truth_vcfs:
    scenario = scenario_name_for(truth_vcf)
    pr_json = OUT_DIR / scenario / "pr_metrics.json"
    if not pr_json.exists():
        continue
    try:
        data = json.loads(pr_json.read_text(encoding="utf-8"))
    except Exception:
        continue
    for row in data.get("per_type_rows", []):
        svtype = row.get("svtype", "?")
        bucket = per_type_counts.setdefault(svtype, {"tp": 0, "fp": 0, "fn": 0})
        bucket["tp"] += int(row.get("tp", 0) or 0)
        bucket["fp"] += int(row.get("fp", 0) or 0)
        bucket["fn"] += int(row.get("fn", 0) or 0)

per_type_tsv = OUT_DIR / "pr_metrics_small_tests_by_svtype.tsv"
with per_type_tsv.open("w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=["svtype", "tp", "fp", "fn", "precision", "recall", "f1"],
        delimiter="\t",
    )
    writer.writeheader()
    for svtype, counts in sorted(per_type_counts.items()):
        tp = counts["tp"]; fp = counts["fp"]; fn = counts["fn"]
        prec = tp / (tp + fp) if (tp + fp) else 0.0
        rec = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = (2 * prec * rec / (prec + rec)) if (prec + rec) else 0.0
        writer.writerow({
            "svtype": svtype, "tp": tp, "fp": fp, "fn": fn,
            "precision": round(prec, 6),
            "recall": round(rec, 6),
            "f1": round(f1, 6),
        })

print("")
print(f"Summary TSV:          {summary_tsv}")
print(f"By-SV-type TSV:       {per_type_tsv}")
print(f"Scenarios scored:     {sum(1 for r in summary_rows if not r['notes'])}")
print(f"Scenarios skipped:    {sum(1 for r in summary_rows if r['notes'])}")

if failures:
    print("")
    print("Failures during scoring:")
    for f in failures:
        print(f"  - {f}")

# Non-zero exit only if nothing at all could be scored — lets downstream
# shells keep going on partial success.
if not any(not r["notes"] for r in summary_rows):
    sys.exit(4)
PYTHON_SCRIPT

rc=$?

echo ""
if [[ $rc -eq 0 ]]; then
  echo "✓ Metrics generation complete"
else
  echo "✗ Metrics generation finished with rc=${rc} — see messages above."
fi
exit $rc
