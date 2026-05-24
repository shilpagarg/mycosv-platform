#!/usr/bin/env bash
#SBATCH --job-name=mycosv-fus-debug
#SBATCH -p multicore
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --array=0-1
#SBATCH --output=slurm-debug-fusarium-modes-%A-%a.out

set -euo pipefail

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "${PROJECT_DIR}"

COMPARATOR_ENV="/mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/tools/envs/envs/fungi_graph_sv"
if [[ -d "${COMPARATOR_ENV}/bin" ]]; then
  export PATH="${COMPARATOR_ENV}/bin:${PATH}"
fi

PREPARED_DIR="${PREPARED_DIR:-experiments/million_real/20260518_105922}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
DEBUG_TAG="${DEBUG_TAG:-${SLURM_ARRAY_JOB_ID:-manual}}"
OUT_ROOT="${OUT_ROOT:-${PREPARED_DIR}/debug_fusarium_modes_${DEBUG_TAG}}"

MODES=("assembly" "long-reads")
OUTNAMES=("assembly" "long_reads")

idx="${SLURM_ARRAY_TASK_ID:-0}"
if (( idx < 0 || idx > 1 )); then
  echo "[error] SLURM_ARRAY_TASK_ID=${idx} out of range 0..1" >&2
  exit 2
fi

MODE="${MODES[$idx]}"
OUT_NAME="${OUTNAMES[$idx]}"
OUT_DIR="${OUT_ROOT}/${OUT_NAME}"

case "${MODE}" in
  assembly)
    COMPARATOR_FLAGS=(--run-minigraph --run-svim-asm --run-anchorwave)
    READVAL_SUPPORT="${DEBUG_READVAL_SUPPORT_ASM:-1}"
    ;;
  long-reads)
    COMPARATOR_FLAGS=(--run-svim --run-sniffles --run-cutesv)
    READVAL_SUPPORT="${DEBUG_READVAL_SUPPORT_LR:-2}"
    ;;
  *)
    echo "[error] unsupported MODE=${MODE}" >&2
    exit 2
    ;;
esac

LR_CAP_FLAGS=()
if [[ "${MODE}" == "long-reads" && "${DEBUG_MAX_LONG_READS:-0}" -gt 0 ]]; then
  LR_CAP_FLAGS=(--max-comparator-long-reads "${DEBUG_MAX_LONG_READS}")
fi

export MYCOSV_TOOL_TIMEOUT="${MYCOSV_TOOL_TIMEOUT:-10800}"
export MILLION_REAL_SINGLE_REF_CACHE_MB="${MILLION_REAL_SINGLE_REF_CACHE_MB:-4096}"
export MILLION_REAL_MAX_REF_MEMORY_MB="${MILLION_REAL_MAX_REF_MEMORY_MB:-8192}"
unset MYCOSV_FORCE_FLAT_REF_FALLBACK || true

mkdir -p "${OUT_ROOT}"

echo "=== MycoSV Fusarium debug ${MODE} ==="
echo "start:        $(date)"
echo "project:      ${PROJECT_DIR}"
echo "prepared:     ${PREPARED_DIR}"
echo "out_root:     ${OUT_ROOT}"
echo "out_dir:      ${OUT_DIR}"
echo "threads:      ${THREADS}"
echo "comparators:  ${COMPARATOR_FLAGS[*]}"
echo "readval_min:  ${READVAL_SUPPORT}"
echo "lr_cap:       ${DEBUG_MAX_LONG_READS:-default}"
echo "binary:       ${PROJECT_DIR}/fungi_graphsv_tol_bin"
stat -c "binary_mtime: %y" "${PROJECT_DIR}/fungi_graphsv_tol_bin"
echo "samtools:     $(command -v samtools || echo MISSING)"
echo "minimap2:     $(command -v minimap2 || echo MISSING)"
echo "minigraph:    $(command -v minigraph || echo MISSING)"
echo "svim-asm:     $(command -v svim-asm || echo MISSING)"
echo "anchorwave:   $(command -v anchorwave || echo MISSING)"
echo "svim:         $(command -v svim || echo MISSING)"
echo "sniffles:     $(command -v sniffles || echo MISSING)"
echo "cuteSV:       $(command -v cuteSV || command -v cutesv || echo MISSING)"

python3 -u run_real_fungal_benchmark.py benchmark \
  --prepared-dir "${PREPARED_DIR}" \
  --out-dir "${OUT_DIR}" \
  --mode "${MODE}" \
  --threads "${THREADS}" \
  --max-clade-genomes 32 \
  --benchmark-query-genera Fusarium \
  --max-benchmark-queries 1 \
  --reuse-index-dir "${PREPARED_DIR}/index" \
  --reuse-registry-dir "${PREPARED_DIR}/registry" \
  --benchmark-ref-cap 8 \
  --read-validation-min-support "${READVAL_SUPPORT}" \
  "${COMPARATOR_FLAGS[@]}" \
  "${LR_CAP_FLAGS[@]}" \
  --mycosv-arg=--max-calls-per-contig --mycosv-arg=2000 \
  --mycosv-arg=--min-block-score --mycosv-arg=4.0 \
  --mycosv-arg=--tol-min-chain-anchors --mycosv-arg=2 \
  --mycosv-arg=--max-ref-memory-mb --mycosv-arg=8192 \
  --mycosv-arg=--max-flat-ref-contigs --mycosv-arg=256 \
  --mycosv-arg=--skip-flat-if-hier-calls --mycosv-arg=5 \
  --mycosv-arg=--no-gfa

python3 - <<'PY' "${OUT_DIR}" "${MODE}"
import csv
import json
import sys
from collections import Counter
from pathlib import Path

out_dir = Path(sys.argv[1])
mode = sys.argv[2]
rows = []

def count_tsv(rel):
    path = out_dir / rel
    if not path.exists():
        return {"file": rel, "status": "missing", "rows": "."}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        data = list(reader)
    status_counts = Counter(row.get("status", ".") for row in data)
    return {
        "file": rel,
        "status": "ok",
        "rows": str(len(data)),
        "summary": ";".join(f"{k}:{v}" for k, v in status_counts.most_common(5)) or ".",
    }

def count_vcf(rel):
    path = out_dir / rel
    if not path.exists():
        return {"file": rel, "status": "missing", "rows": "."}
    n = 0
    sv = Counter()
    qmode = Counter()
    support = Counter()
    with path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            n += 1
            fields = line.rstrip("\n").split("\t")
            info = {}
            for item in fields[7].split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info[k] = v
            sv[info.get("SVTYPE", ".")] += 1
            qmode[info.get("QMODE", ".")] += 1
            support[info.get("SUPPORT", ".")] += 1
    return {
        "file": rel,
        "status": "ok",
        "rows": str(n),
        "summary": (
            "sv=" + ",".join(f"{k}:{v}" for k, v in sv.most_common()) +
            "|qmode=" + ",".join(f"{k}:{v}" for k, v in qmode.most_common()) +
            "|support_top=" + ",".join(f"{k}:{v}" for k, v in support.most_common(3))
        ),
    }

for rel in (
    "mycosv/calls.vcf",
    "mycosv/calls.multisample.vcf",
    "mycosv/calls.hierarchical.vcf",
):
    rows.append(count_vcf(rel))

for rel in (
    "read_validated_truth.tsv",
    "exact_benchmark_summary.tsv",
    "loo_consensus_summary.tsv",
    "match_failures.tsv",
    "biology_findings.tsv",
    "biology_candidates.tsv",
    "novel_mycosv_calls.tsv",
    "mycosv_validation_followup.tsv",
    "pangenome_call_layers.tsv",
    "sv_volume_audit.tsv",
):
    rows.append(count_tsv(rel))

summary = out_dir / "benchmark_summary.json"
if summary.exists():
    payload = json.loads(summary.read_text())
    rows.append({
        "file": "benchmark_summary.json",
        "status": "ok",
        "rows": ".",
        "summary": f"mode={payload.get('mode', mode)};queries={len(payload.get('queries', {}))}",
    })

audit = out_dir / "debug_step_audit.tsv"
with audit.open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=["file", "status", "rows", "summary"], delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
print(f"[audit] wrote {audit}")
PY

REPORT_OUT="${OUT_DIR}/report"
if [[ -f "${OUT_DIR}/exact_benchmark_summary.tsv" && -f "${OUT_DIR}/novel_mycosv_calls.tsv" ]]; then
  echo
  echo "=== Visualization report for ${MODE} ==="
  mkdir -p "${REPORT_OUT}"
  REPORT_ARGS=(
    --real "${OUT_DIR}/exact_benchmark_summary.tsv"
    --novel "${OUT_DIR}/novel_mycosv_calls.tsv"
    --outdir "${REPORT_OUT}"
    --title "MycoSV Fusarium ${MODE} debug report"
  )
  [[ -f "${OUT_DIR}/biology_findings.tsv" ]] \
    && REPORT_ARGS+=(--biology "${OUT_DIR}/biology_findings.tsv")
  [[ -f "${OUT_DIR}/mycosv_evidence_tiers.tsv" ]] \
    && REPORT_ARGS+=(--evidence-tiers "${OUT_DIR}/mycosv_evidence_tiers.tsv")
  if python3 -u sv_visualization_report.py "${REPORT_ARGS[@]}" \
      > "${REPORT_OUT}/report.log" 2>&1; then
    cat "${REPORT_OUT}/report.log"
  else
    rc=$?
    cat "${REPORT_OUT}/report.log" >&2 || true
    echo "[warn] report generation failed (rc=${rc}); benchmark outputs remain available under ${OUT_DIR}" >&2
  fi
else
  echo "[skip] report: exact_benchmark_summary.tsv / novel_mycosv_calls.tsv missing under ${OUT_DIR}"
fi

if [[ -f "${OUT_DIR}/novel_mycosv_calls.tsv" && -f "${OUT_DIR}/pangenome_call_layers.tsv" ]]; then
  echo
  echo "=== MycoSV pangenome-call plots for ${MODE} ==="
  python3 -u plot_mycosv_pangenome_calls.py \
    --benchmark-dir "${OUT_DIR}" \
    --outdir "${OUT_DIR}/pangenome_plots" \
    --title "MycoSV Fusarium ${MODE} pangenome-call biology" \
    || echo "[warn] pangenome-call plotting failed (rc=$?)" >&2
else
  echo "[skip] pangenome plots: novel_mycosv_calls.tsv / pangenome_call_layers.tsv missing under ${OUT_DIR}"
fi

echo "finish:       $(date)"
