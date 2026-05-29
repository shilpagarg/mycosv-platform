#!/usr/bin/env bash
#SBATCH --job-name=mycosv-full-fungi-asm
#SBATCH -p multicore
#SBATCH --mem=256G
#SBATCH --cpus-per-task=32
#SBATCH --time=72:00:00
#SBATCH --output=slurm-full-fungal-assembly-%j.out

set -euo pipefail

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "${PROJECT_DIR}"

COMPARATOR_ENV="/mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/tools/envs/envs/fungi_graph_sv"
if [[ -d "${COMPARATOR_ENV}/bin" ]]; then
  export PATH="${COMPARATOR_ENV}/bin:${PATH}"
fi

THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-32}}"
TIMESTAMP="${TIMESTAMP:-$(date +%Y%m%d_%H%M%S)}"
OUT_ROOT="${OUT_ROOT:-experiments/million_real/full_fungal_assembly_${TIMESTAMP}}"
PREPARED_DIR="${PREPARED_DIR:-${OUT_ROOT}/prepared}"
ASM_OUT="${ASM_OUT:-${OUT_ROOT}/assembly}"
DATA_CACHE_DIR="${DATA_CACHE_DIR:-${PROJECT_DIR}/data_cache}"

# Held-out query groups requested for the fungal-biology community panel.
# Misspellings such as mychorrzia/penchillium are normalized by the Python
# driver, but keep the defaults canonical for readable logs.
QUERY_GROUPS="${QUERY_GROUPS:-Aspergillus,Candida,mycorrhiza,Penicillium,Trichoderma,Neurospora crassa,Fusarium}"
QUERY_COUNT="${QUERY_COUNT:-15}"

# max-assemblies=0 means all currently discoverable public fungal assemblies
# from the NCBI-best catalog are cached/downloaded. The routing store is then
# padded to TARGET_CENTROIDS, so the experiment remains the million-scale
# benchmark even though NCBI currently exposes far fewer than 1,000,000 real
# fungal assemblies.
MAX_ASSEMBLIES="${MAX_ASSEMBLIES:-0}"
TARGET_CENTROIDS="${TARGET_CENTROIDS:-1000000}"
MAX_CLADE_GENOMES="${MAX_CLADE_GENOMES:-32}"
# Match the successful debug mode by default: each full-panel query is run as
# its own bounded shard, so 8 nearby benchmark refs is enough and avoids the
# 64-ref monolithic MycoSV timeout seen in the first full run.
BENCHMARK_REF_CAP="${BENCHMARK_REF_CAP:-8}"
READVAL_SUPPORT="${READVAL_SUPPORT:-1}"
# Cap each assembly query to its N longest contigs. MycoSV's hierarchical
# assembly caller scales with query contig count and stalls for hours above a
# few hundred contigs (fragmented public drafts run 2.7k-180k contigs). 0
# disables (every contig processed, original behaviour). When set, oversized
# queries are kept in the panel as their N-longest-contig subset rather than
# hanging or being skipped; see TRUNCATED_ASSEMBLY_QUERIES.tsv.
MAX_ASSEMBLY_QUERY_CONTIGS_KEEP="${MAX_ASSEMBLY_QUERY_CONTIGS_KEEP:-0}"
SEED="${SEED:-1}"
FULL_ASSEMBLY_SHARDS="${FULL_ASSEMBLY_SHARDS:-1}"
# Default to a generated read-validation manifest under prepared/ so the
# read_validated_truth.tsv stops being filled with `validation_unavailable`
# stubs. Override with explicit RAW_READ_VALIDATION_TSV or set to "" to
# disable read-validation entirely.
RAW_READ_VALIDATION_TSV="${RAW_READ_VALIDATION_TSV-${OUT_ROOT}/prepared/read_validation_manifest.tsv}"
RAW_READ_VALIDATION_MAX_READS="${RAW_READ_VALIDATION_MAX_READS:-200000}"
# Cap on rows kept in biology_candidates.tsv / biology_findings.tsv. The
# previous hard-coded 200 throttled per-query biology to a token sample; lift
# to 5000 so HGT/Starship/RIP/TE candidates across all SV types come through.
BIOLOGY_TOP_N="${BIOLOGY_TOP_N:-5000}"
export MYCOSV_BIOLOGY_TOP_N="${BIOLOGY_TOP_N}"
RESUME_SHARDS="${RESUME_SHARDS:-1}"
FORCE_RERUN_SHARDS="${FORCE_RERUN_SHARDS:-0}"
INDEX_CACHE_ID="${INDEX_CACHE_ID:-ncbi_best_q${QUERY_COUNT}_seed${SEED}_centroids${TARGET_CENTROIDS}_clade${MAX_CLADE_GENOMES}}"
INDEX_CACHE_DIR="${INDEX_CACHE_DIR:-${DATA_CACHE_DIR}/mycosv_indices/${INDEX_CACHE_ID}/index}"
REGISTRY_CACHE_DIR="${REGISTRY_CACHE_DIR:-${DATA_CACHE_DIR}/mycosv_indices/${INDEX_CACHE_ID}/registry}"
GENE_ANNOTATION_CACHE="${GENE_ANNOTATION_CACHE:-${DATA_CACHE_DIR}/gene_annotations/${INDEX_CACHE_ID}.gene_annotations.tsv}"

export MYCOSV_TOOL_TIMEOUT="${MYCOSV_TOOL_TIMEOUT:-21600}"
export MYCOSV_COMPARATOR_TIMEOUT="${MYCOSV_COMPARATOR_TIMEOUT:-5400}"
export MYCOSV_READ_VALIDATION_VIEW_TIMEOUT="${MYCOSV_READ_VALIDATION_VIEW_TIMEOUT:-30}"
export MILLION_REAL_DOWNLOAD_WORKERS="${MILLION_REAL_DOWNLOAD_WORKERS:-6}"
export MILLION_REAL_SINGLE_REF_CACHE_MB="${MILLION_REAL_SINGLE_REF_CACHE_MB:-8192}"
export MILLION_REAL_MAX_REF_MEMORY_MB="${MILLION_REAL_MAX_REF_MEMORY_MB:-16384}"

mkdir -p "${OUT_ROOT}"

echo "=== MycoSV full fungal assembly experiment ==="
echo "start:             $(date)"
echo "project:           ${PROJECT_DIR}"
echo "out_root:          ${OUT_ROOT}"
echo "prepared_dir:      ${PREPARED_DIR}"
echo "assembly_out:      ${ASM_OUT}"
echo "data_cache_dir:    ${DATA_CACHE_DIR}"
echo "threads:           ${THREADS}"
echo "max_assemblies:    ${MAX_ASSEMBLIES}"
echo "target_centroids:  ${TARGET_CENTROIDS}"
echo "query_groups:      ${QUERY_GROUPS}"
echo "query_count:       ${QUERY_COUNT}"
echo "full_shards:       ${FULL_ASSEMBLY_SHARDS}"
echo "benchmark_ref_cap: ${BENCHMARK_REF_CAP}"
echo "raw_read_val_tsv:  ${RAW_READ_VALIDATION_TSV:-none}"
echo "raw_read_val_cap:  ${RAW_READ_VALIDATION_MAX_READS}"
echo "resume_shards:     ${RESUME_SHARDS}"
echo "force_rerun:       ${FORCE_RERUN_SHARDS}"
echo "index_cache_dir:   ${INDEX_CACHE_DIR}"
echo "registry_cache:    ${REGISTRY_CACHE_DIR}"
echo "gene_annot_cache:  ${GENE_ANNOTATION_CACHE}"
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "array_task_id:     ${SLURM_ARRAY_TASK_ID}"
fi
echo "minimap2:          $(command -v minimap2 || echo MISSING)"
echo "samtools:          $(command -v samtools || echo MISSING)"
echo "minigraph:         $(command -v minigraph || echo MISSING)"
echo "svim-asm:          $(command -v svim-asm || echo MISSING)"
echo "anchorwave:        $(command -v anchorwave || echo MISSING)"

# Prebuild the MycoSV binary before any per-shard work. SLURM array jobs
# would otherwise have N tasks racing to compile the same output path; the
# Python compile_binary_if_needed() also takes an flock, but a single
# serialized build here keeps the per-task path a fast no-op.
#
# The lock file lives on node-local TMPDIR (or /tmp) because flock(2) on the
# bmh01-rds NFS share returns ENOLCK ("No locks available") under array
# contention — the prior run lost 8/15 shards to that exact failure. With the
# lock local to the node, array tasks only contend with siblings on the same
# node (rare for 32-CPU tasks) and a fresh binary mtime lets them no-op the
# lock entirely.
MYCOSV_BIN="${MYCOSV_BIN:-${PROJECT_DIR}/fungi_graphsv_tol_bin}"
BIN_LOCK="${BIN_LOCK:-${TMPDIR:-/tmp}/$(basename "${MYCOSV_BIN}").${USER}.lock}"
if [[ -x "${MYCOSV_BIN}" && "${FORCE_PREBUILD:-0}" != "1" ]]; then
  BIN_MTIME=$(stat -c %Y "${MYCOSV_BIN}" 2>/dev/null || echo 0)
  # Newest of main.cpp AND every *.hpp — must match the Python's _needs_build()
  # which also checks headers. Comparing only main.cpp meant a newer .hpp left
  # the binary "stale" by the Python's reckoning, so each array task tried to
  # lock+rebuild on the NFS share and ~half hit ENOLCK ("No locks available").
  SRC_MTIME=$(stat -c %Y "${PROJECT_DIR}/main.cpp" ${PROJECT_DIR}/*.hpp 2>/dev/null | sort -nr | head -1)
  SRC_MTIME="${SRC_MTIME:-0}"
  if (( BIN_MTIME >= SRC_MTIME && BIN_MTIME > 0 )); then
    echo "[prebuild] ${MYCOSV_BIN} mtime ${BIN_MTIME} >= newest source ${SRC_MTIME} (main.cpp + *.hpp); skipping lock+rebuild"
    SKIP_PREBUILD=1
  fi
fi
if [[ "${SKIP_PREBUILD:-0}" != "1" ]]; then
  echo "[prebuild] ensuring ${MYCOSV_BIN} is up to date (lock=${BIN_LOCK})"
  mkdir -p "$(dirname "${BIN_LOCK}")"
  if ! (
    flock -w 300 -x 9 || exit 75
    python3 -u - <<PY
from pathlib import Path
import sys
sys.path.insert(0, "${PROJECT_DIR}")
from run_real_fungal_benchmark import compile_binary_if_needed
compile_binary_if_needed(Path("${MYCOSV_BIN}").resolve())
print("[prebuild] binary is ready:", "${MYCOSV_BIN}")
PY
  ) 9>"${BIN_LOCK}"; then
    rc=$?
    if (( rc == 75 )) && [[ -x "${MYCOSV_BIN}" ]]; then
      echo "[prebuild] flock unavailable (rc=${rc}); binary exists at ${MYCOSV_BIN}, proceeding without lock"
    else
      echo "[prebuild] flock or compile failed (rc=${rc}) and no binary present at ${MYCOSV_BIN}" >&2
      exit 2
    fi
  fi
fi
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" && "${FORCE_ARRAY_PREPARE:-0}" != "1" && -s "${PREPARED_DIR}/prepare_million_real_summary.json" ]]; then
  echo "[array] using existing prepared/index; set FORCE_ARRAY_PREPARE=1 only if you really want each task to prepare"
elif [[ -n "${SLURM_ARRAY_TASK_ID:-}" && "${FORCE_ARRAY_PREPARE:-0}" != "1" ]]; then
  echo "[array] missing ${PREPARED_DIR}/prepare_million_real_summary.json; run one prepare job first, or set FORCE_ARRAY_PREPARE=1 for a single task only" >&2
  exit 2
elif [[ "${SKIP_PREPARE:-0}" != "1" || "${FORCE_ARRAY_PREPARE:-0}" == "1" ]]; then
  python3 -u run_real_fungal_benchmark.py prepare-million-real \
    --out-dir "${PREPARED_DIR}" \
    --data-cache-dir "${DATA_CACHE_DIR}" \
    --source ncbi-best \
    --max-assemblies "${MAX_ASSEMBLIES}" \
    --min-assembly-level "${MIN_ASSEMBLY_LEVEL:-contig}" \
    --target-centroids "${TARGET_CENTROIDS}" \
    --million-real-queries "${QUERY_COUNT}" \
    --million-real-query-genera "${QUERY_GROUPS}" \
    --threads "${THREADS}" \
    --max-clade-genomes "${MAX_CLADE_GENOMES}" \
    --seed "${SEED}" \
    --million-real-index-cache-dir "${INDEX_CACHE_DIR}" \
    --million-real-registry-cache-dir "${REGISTRY_CACHE_DIR}" \
    --million-real-gene-annotations-cache "${GENE_ANNOTATION_CACHE}" \
    --million-real-download-gff \
    --million-real-phenotypes
else
  echo "[skip] prepare-million-real because SKIP_PREPARE=1"
fi

if [[ -s "${INDEX_CACHE_DIR}/routing_manifest.tsv" ]]; then
  BENCHMARK_INDEX_DIR="${INDEX_CACHE_DIR}"
  BENCHMARK_REGISTRY_DIR="${REGISTRY_CACHE_DIR}"
  echo "[index] benchmark will reuse data-cache index: ${BENCHMARK_INDEX_DIR}"
elif [[ -s "${PREPARED_DIR}/index/routing_manifest.tsv" ]]; then
  BENCHMARK_INDEX_DIR="${PREPARED_DIR}/index"
  BENCHMARK_REGISTRY_DIR="${PREPARED_DIR}/registry"
  echo "[index] benchmark will reuse prepared-dir index: ${BENCHMARK_INDEX_DIR}"
else
  echo "[index] no reusable routing index found in ${INDEX_CACHE_DIR} or ${PREPARED_DIR}/index" >&2
  exit 2
fi

python3 - <<'PY' "${PREPARED_DIR}" "${QUERY_COUNT}" "${QUERY_GROUPS}" "${MODE:-assembly}"
import csv
import sys
from pathlib import Path

prepared = Path(sys.argv[1])
expected = int(sys.argv[2])
groups = sys.argv[3]
mode = sys.argv[4]
manifest = prepared / "query_manifest.tsv"
rows = list(csv.DictReader(manifest.open(), delimiter="\t"))
sel = [r for r in rows if (r.get("query_mode") or "assembly") == mode]
print(f"[query-check] mode={mode} queries={len(sel)} expected={expected} groups={groups}")
if len(sel) < expected:
    raise SystemExit(f"expected at least {expected} {mode} queries, found {len(sel)}")
PY

# BOOTSTRAP_ONLY=1 short-circuit: the prepare + index/manifest-check above
# is everything a bootstrap-phase job needs to do. The read-validation
# manifest is then built (idempotent if already populated by prepare) and the
# script exits before launching the benchmark, so the bootstrap slot doesn't
# accidentally start running mycosv on the head node of a 165-task array.
if [[ "${BOOTSTRAP_ONLY:-0}" == "1" ]]; then
  READ_MANIFEST_OUT="${READ_VALIDATION_MANIFEST:-${PROJECT_DIR}/${OUT_ROOT}/prepared/read_validation_manifest.tsv}"
  if [[ ! -s "${READ_MANIFEST_OUT}" ]]; then
    mkdir -p "$(dirname "${READ_MANIFEST_OUT}")"
    python3 -u "${PROJECT_DIR}/build_read_validation_manifest.py" \
      --query-manifest "${PROJECT_DIR}/${PREPARED_DIR}/query_manifest.tsv" \
      --reads-cache "${DATA_CACHE_DIR}/raw_reads" \
      --max-bases "${READ_VALIDATION_MAX_BASES:-300000000}" \
      --out-tsv "${READ_MANIFEST_OUT}" \
      || echo "[bootstrap] read-manifest build returned non-zero; manifest may have partial entries" >&2
  else
    echo "[bootstrap] read-validation manifest already present at ${READ_MANIFEST_OUT}"
  fi
  echo "[bootstrap] done at $(date) — prepared/ ready for array tasks"
  exit 0
fi

run_benchmark_dir() {
  local out_dir="$1"
  local query_asm="${2:-}"
  local query_args=()
  local raw_read_args=()
  if [[ -n "${query_asm}" ]]; then
    query_args=(--benchmark-query-asm "${query_asm}")
  fi
  if [[ -n "${RAW_READ_VALIDATION_TSV}" ]]; then
    raw_read_args=(
      --raw-read-validation-tsv "${RAW_READ_VALIDATION_TSV}"
      --raw-read-validation-max-reads "${RAW_READ_VALIDATION_MAX_READS}"
    )
  fi
  # Comparator selection per benchmark mode.
  #   assembly    : minigraph / svim-asm / anchorwave (always on); PGGB + Cactus
  #                 optional via RUN_PGGB / RUN_CACTUS. Reviewers (NM, NB)
  #                 expect head-to-head against modern pangenome graph callers
  #                 (Liao 2023 HPRC, Hickey 2024 Minigraph-Cactus, Garrison 2018 vg).
  #   long-reads  : Sniffles2 / cuteSV / SVIM — the canonical PB/ONT SV callers
  #                 (Heller 2022 Genome Biol, Jiang 2020 GB, Heller 2019 Bioinf).
  #   short-reads : Delly / Manta — canonical Illumina SV callers (Rausch 2012,
  #                 Chen 2016 Bioinf). Off by default; opt in with MODE=short-reads.
  # MODE is an env var so the same launcher handles assembly + reads modes
  # without forking the script — per the project policy of editing existing
  # pipeline + heredocs rather than spawning new shell scripts.
  local mode="${MODE:-assembly}"
  local comparator_flags=()
  # MYCOSV_ONLY=1 short-circuits every comparator (minigraph/svim-asm/
  # anchorwave/pggb/cactus/syri and the reads-mode baselines). The run then
  # rests entirely on MycoSV's graph-native calls + independent read-level
  # validation. This is the comparator-free full-fungal experiment: each
  # comparator carries a 90-min timeout and was the dominant per-shard
  # wall-time cost, so disabling them is what brings the 200-query panel
  # inside a single-digit-hour budget.
  if [[ "${MYCOSV_ONLY:-0}" == "1" ]]; then
    comparator_flags=(--mycosv-only)
    echo "[mode] ${mode} (mycosv-only: comparators disabled)"
    echo "[comparators] enabled: none (--mycosv-only)"
  else
    case "${mode}" in
      assembly)
        comparator_flags=(--run-minigraph --run-svim-asm --run-anchorwave)
        if [[ "${RUN_PGGB:-1}" == "1" ]]; then
          comparator_flags+=(--run-pggb)
        fi
        if [[ "${RUN_CACTUS:-1}" == "1" ]]; then
          comparator_flags+=(--run-cactus)
        fi
        if [[ "${RUN_SYRI:-0}" == "1" ]]; then
          comparator_flags+=(--run-syri)
        fi
        ;;
      long-reads)
        comparator_flags=(--run-sniffles --run-cutesv --run-svim)
        ;;
      short-reads)
        comparator_flags=(--run-delly --run-manta)
        ;;
      *)
        echo "[error] unsupported MODE=${mode} (expected: assembly|long-reads|short-reads)" >&2
        return 2
        ;;
    esac
    echo "[mode] ${mode}"
    echo "[comparators] enabled: ${comparator_flags[*]}"
  fi
  python3 -u run_real_fungal_benchmark.py benchmark \
    --prepared-dir "${PREPARED_DIR}" \
    --out-dir "${out_dir}" \
    --mode "${mode}" \
    --threads "${THREADS}" \
    --max-clade-genomes "${MAX_CLADE_GENOMES}" \
    --reuse-index-dir "${BENCHMARK_INDEX_DIR}" \
    --reuse-registry-dir "${BENCHMARK_REGISTRY_DIR}" \
    --benchmark-ref-cap "${BENCHMARK_REF_CAP}" \
    --read-validation-min-support "${READVAL_SUPPORT}" \
    --max-assembly-query-contigs-keep-largest "${MAX_ASSEMBLY_QUERY_CONTIGS_KEEP}" \
    "${comparator_flags[@]}" \
    "${query_args[@]}" \
    "${raw_read_args[@]}" \
    --mycosv-arg=--max-calls-per-contig --mycosv-arg="${MAX_CALLS_PER_CONTIG:-2000}" \
    --mycosv-arg=--min-block-score --mycosv-arg="${MIN_BLOCK_SCORE:-4.0}" \
    --mycosv-arg=--tol-min-chain-anchors --mycosv-arg="${TOL_MIN_CHAIN_ANCHORS:-2}" \
    --mycosv-arg=--graph-native-mode \
    --mycosv-arg=--tol-base-graph-build \
    --mycosv-arg=--max-ref-memory-mb --mycosv-arg="${MYCOSV_MAX_REF_MEMORY_MB:-8192}" \
    --mycosv-arg=--max-flat-ref-contigs --mycosv-arg="${MAX_FLAT_REF_CONTIGS:-256}" \
    --mycosv-arg=--skip-flat-if-hier-calls --mycosv-arg="${SKIP_FLAT_IF_HIER_CALLS:-5}"
}

write_debug_audit() {
  local out_dir="$1"
  python3 - <<'PY' "${out_dir}"
import csv
import json
import sys
from collections import Counter
from pathlib import Path

out_dir = Path(sys.argv[1])
rows = []

def count_tsv(rel):
    path = out_dir / rel
    if not path.exists():
        return {"file": rel, "status": "missing", "rows": ".", "summary": "."}
    with path.open() as fh:
        data = list(csv.DictReader(fh, delimiter="\t"))
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
        return {"file": rel, "status": "missing", "rows": ".", "summary": "."}
    n = 0
    sv = Counter()
    qasm = Counter()
    with path.open(errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            n += 1
            fields = line.rstrip("\n").split("\t")
            info = {}
            if len(fields) > 7:
                for item in fields[7].split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        info[k] = v
            sv[info.get("SVTYPE", ".")] += 1
            qasm[info.get("QASM", ".")] += 1
    return {
        "file": rel,
        "status": "ok",
        "rows": str(n),
        "summary": (
            "sv=" + ",".join(f"{k}:{v}" for k, v in sv.most_common()) +
            "|qasm=" + ",".join(f"{k}:{v}" for k, v in qasm.most_common(5))
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
    "read_validation_summary.tsv",
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
        "summary": f"mode={payload.get('mode', 'assembly')};queries={len(payload.get('queries', {}))}",
    })

audit = out_dir / "debug_step_audit.tsv"
with audit.open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=["file", "status", "rows", "summary"], delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
print(f"[audit] wrote {audit}")
PY
}

build_reports() {
  local out_dir="$1"
  local title="$2"
  local report_out="${out_dir}/report"
  if [[ -f "${out_dir}/exact_benchmark_summary.tsv" && -f "${out_dir}/novel_mycosv_calls.tsv" ]]; then
    mkdir -p "${report_out}"
    local report_args=(
      --real "${out_dir}/exact_benchmark_summary.tsv"
      --novel "${out_dir}/novel_mycosv_calls.tsv"
      --outdir "${report_out}"
      --title "${title}"
    )
    [[ -f "${out_dir}/biology_findings.tsv" ]] \
      && report_args+=(--biology "${out_dir}/biology_findings.tsv")
    [[ -f "${out_dir}/mycosv_evidence_tiers.tsv" ]] \
      && report_args+=(--evidence-tiers "${out_dir}/mycosv_evidence_tiers.tsv")
    python3 -u sv_visualization_report.py "${report_args[@]}" \
      > "${report_out}/report.log" 2>&1 || cat "${report_out}/report.log" >&2
  else
    echo "[skip] report for ${out_dir}: exact_benchmark_summary.tsv / novel_mycosv_calls.tsv missing"
  fi

  if [[ -f "${out_dir}/novel_mycosv_calls.tsv" && -f "${out_dir}/pangenome_call_layers.tsv" ]]; then
    python3 -u plot_mycosv_pangenome_calls.py \
      --benchmark-dir "${out_dir}" \
      --outdir "${out_dir}/pangenome_plots" \
      --title "${title} pangenome-call biology" \
      || echo "[warn] pangenome-call plotting failed for ${out_dir}" >&2
  else
    echo "[skip] pangenome plots for ${out_dir}: novel_mycosv_calls.tsv / pangenome_call_layers.tsv missing"
  fi
}

combine_shard_tsvs() {
  local combined_dir="$1"
  shift
  mkdir -p "${combined_dir}"
  local files=(
    read_validated_truth.tsv
    read_validation_summary.tsv
    exact_benchmark_summary.tsv
    loo_consensus_summary.tsv
    match_failures.tsv
    biology_findings.tsv
    biology_candidates.tsv
    novel_mycosv_calls.tsv
    mycosv_validation_followup.tsv
    mycosv_evidence_tiers.tsv
    pangenome_call_layers.tsv
    sv_volume_audit.tsv
  )
  local rel shard out first
  for rel in "${files[@]}"; do
    out="${combined_dir}/${rel}"
    first=1
    : > "${out}"
    for shard in "$@"; do
      [[ -f "${shard}/${rel}" ]] || continue
      if (( first )); then
        cat "${shard}/${rel}" >> "${out}"
        first=0
      else
        tail -n +2 "${shard}/${rel}" >> "${out}"
      fi
    done
    if (( first )); then
      rm -f "${out}"
    fi
  done
}

shard_complete() {
  local shard_dir="$1"
  [[ "${RESUME_SHARDS}" == "1" && "${FORCE_RERUN_SHARDS}" != "1" ]] || return 1
  [[ -s "${shard_dir}/benchmark_summary.json" ]] || return 1
  [[ -s "${shard_dir}/mycosv/calls.multisample.vcf" || -s "${shard_dir}/mycosv/calls.vcf" ]] || return 1
  [[ -s "${shard_dir}/pangenome_call_layers.tsv" ]] || return 1
  [[ -s "${shard_dir}/sv_volume_audit.tsv" ]] || return 1
  [[ -s "${shard_dir}/read_validation_summary.tsv" ]] || return 1
  return 0
}

if [[ "${FULL_ASSEMBLY_SHARDS}" == "1" ]]; then
  mkdir -p "${ASM_OUT}"
  SHARD_ROOT="${ASM_OUT}/by_query"
  mkdir -p "${SHARD_ROOT}"
  mapfile -t QUERY_ASMS < <(python3 - <<'PY' "${PREPARED_DIR}" "${MODE:-assembly}"
import csv
import sys
from pathlib import Path

prepared = Path(sys.argv[1])
mode = sys.argv[2]
rows = list(csv.DictReader((prepared / "query_manifest.tsv").open(), delimiter="\t"))
for row in rows:
    if (row.get("query_mode") or "assembly") == mode:
        print(row["query_asm"])
PY
)
  if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    if (( SLURM_ARRAY_TASK_ID < 0 || SLURM_ARRAY_TASK_ID >= ${#QUERY_ASMS[@]} )); then
      echo "[array] SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} outside query range 0..$((${#QUERY_ASMS[@]} - 1))"
      exit 2
    fi
    QUERY_ASMS=("${QUERY_ASMS[${SLURM_ARRAY_TASK_ID}]}")
    SHARD_STATUS="${ASM_OUT}/full_assembly_shards.array_${SLURM_ARRAY_TASK_ID}.tsv"
  else
    SHARD_STATUS="${ASM_OUT}/full_assembly_shards.tsv"
  fi
  printf "query_asm\tout_dir\tstatus\n" > "${SHARD_STATUS}"
  SHARD_DIRS=()
  SHARD_FAILURES=0
  for query_asm in "${QUERY_ASMS[@]}"; do
    shard_dir="${SHARD_ROOT}/${query_asm}"
    SHARD_DIRS+=("${shard_dir}")
    echo
    echo "=== Full fungal assembly shard: ${query_asm} ==="
    if shard_complete "${shard_dir}"; then
      echo "[resume] ${query_asm} already complete at ${shard_dir}; skipping"
      write_debug_audit "${shard_dir}"
      build_reports "${shard_dir}" "MycoSV full fungal assembly shard ${query_asm}"
      printf "%s\t%s\tresume_ok\n" "${query_asm}" "${shard_dir}" >> "${SHARD_STATUS}"
    elif run_benchmark_dir "${shard_dir}" "${query_asm}"; then
      write_debug_audit "${shard_dir}"
      build_reports "${shard_dir}" "MycoSV full fungal assembly shard ${query_asm}"
      printf "%s\t%s\tok\n" "${query_asm}" "${shard_dir}" >> "${SHARD_STATUS}"
    else
      rc=$?
      printf "%s\t%s\tfailed_rc_%s\n" "${query_asm}" "${shard_dir}" "${rc}" >> "${SHARD_STATUS}"
      SHARD_FAILURES=$((SHARD_FAILURES + 1))
      echo "[shards] ${query_asm} failed with rc=${rc}; continuing with remaining queries"
    fi
  done
  echo "[shards] wrote ${SHARD_STATUS}"
  if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    combine_shard_tsvs "${ASM_OUT}/combined" "${SHARD_DIRS[@]}"
    write_debug_audit "${ASM_OUT}/combined"
    build_reports "${ASM_OUT}/combined" "MycoSV full fungal assembly combined shards"
  else
    echo "[array] task ${SLURM_ARRAY_TASK_ID} finished ${QUERY_ASMS[0]}; combine after all array tasks complete"
  fi
  if (( SHARD_FAILURES > 0 )); then
    echo "[shards] completed with ${SHARD_FAILURES} failed shard(s); see ${SHARD_STATUS}"
    exit 1
  fi
else
  run_benchmark_dir "${ASM_OUT}"
  write_debug_audit "${ASM_OUT}"
  build_reports "${ASM_OUT}" "MycoSV full fungal assembly benchmark"
fi

echo "finish:            $(date)"
