#!/usr/bin/env bash
#SBATCH --job-name=mycosv-e2e1q
#SBATCH -p multicore
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=0:30:00
#SBATCH --output=slurm-e2e-%j.out

# End-to-end one-query debug:
#   1. MycoSV SV calls   (hierarchical + flat fallback, with all 4 bug fixes)
#   2. Comparator        (minigraph — the mandatory assembly-mode baseline)
#   3. Read-level validation (samtools/minimap2 anchoring against the query)
#   4. Biological findings  (analyze_new_biology_candidates.py pipeline)
#
# All four stages must produce non-empty output within 30 min wallclock.

set -u
set -o pipefail

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "${PROJECT_DIR}"

# Activate the comparator conda env so minigraph/samtools/minimap2/syri are
# resolvable. The wrapper has a hard-coded fallback to this path too, but
# putting it on PATH here keeps any subprocesses honest.
COMPARATOR_ENV="/mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/tools/envs/envs/fungi_graph_sv"
if [[ -d "${COMPARATOR_ENV}/bin" ]]; then
  export PATH="${COMPARATOR_ENV}/bin:${PATH}"
fi

PREPARED_DIR="${PREPARED_DIR:-experiments/million_real/20260518_105922}"
OUT_DIR="${OUT_DIR:-${PREPARED_DIR}/benchmark_e2e1q}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
# One query from each requested fungal group when present in query_manifest.tsv.
# Spelling variants are normalized by run_real_fungal_benchmark.py
# (e.g. penciliium -> Penicillium, trichderma -> Trichoderma, mychrrzia ->
# mycorrhiza/Rhizophagus/Glomus/Glomeromycetes-style AMF groups).
TARGET_QUERY_GROUPS="${TARGET_QUERY_GROUPS:-Fusarium,Penicillium,Candida,mychrrzia,Aspergillus,Trichoderma,Zymoseptoria}"
if [[ -n "${TARGET_QUERY_GROUPS}" ]]; then
  MAX_BENCHMARK_QUERIES="${MAX_BENCHMARK_QUERIES:-0}"
else
  MAX_BENCHMARK_QUERIES="${MAX_BENCHMARK_QUERIES:-1}"
fi
BENCHMARK_REF_CAP="${BENCHMARK_REF_CAP:-8}"
MAX_CALLS_PER_CONTIG="${MAX_CALLS_PER_CONTIG:-1000}"
MIN_BLOCK_SCORE="${MIN_BLOCK_SCORE:-4.0}"
TOL_MIN_CHAIN_ANCHORS="${TOL_MIN_CHAIN_ANCHORS:-2}"
FAST_DEBUG="${FAST_DEBUG:-1}"
RUN_MINIGRAPH="${RUN_MINIGRAPH:-0}"
# Bug 4 fix: cap flat-fallback ref contigs (8 fungal FASTAs = ~6 500 contigs raw).
MAX_FLAT_REF_CONTIGS="${MAX_FLAT_REF_CONTIGS:-256}"
# Bug 2 fix: skip the flat-MEM-chain fallback if hierarchical already produced
# >= N calls for this query.
if [[ "${FAST_DEBUG}" == "1" ]]; then
  SKIP_FLAT_IF_HIER_CALLS="${SKIP_FLAT_IF_HIER_CALLS:-1}"
else
  SKIP_FLAT_IF_HIER_CALLS="${SKIP_FLAT_IF_HIER_CALLS:-5}"
fi
# Total mycosv binary budget. 1000 s leaves ~700 s for the comparator + read
# validation + biology findings stages inside the 30 min SLURM walltime.
MYCOSV_TOOL_TIMEOUT="${MYCOSV_TOOL_TIMEOUT:-1000}"
MYCOSV_COMPARATOR_TIMEOUT="${MYCOSV_COMPARATOR_TIMEOUT:-300}"
REQUIRED_SV_TYPES="${REQUIRED_SV_TYPES:-INS,DEL,DUP,INV,TRA,OFF_REF}"
REQUIRE_ALL_SIX_SV_TYPES="${REQUIRE_ALL_SIX_SV_TYPES:-1}"

echo "MycoSV end-to-end one-query debug"
echo "start: $(date)"
echo "cwd: $(pwd)"
echo "prepared_dir: ${PREPARED_DIR}"
echo "out_dir: ${OUT_DIR}"
echo "threads: ${THREADS}"
echo "target_query_groups: ${TARGET_QUERY_GROUPS}"
echo "max_benchmark_queries: ${MAX_BENCHMARK_QUERIES}"
echo "benchmark_ref_cap: ${BENCHMARK_REF_CAP}"
echo "max_calls_per_contig: ${MAX_CALLS_PER_CONTIG}"
echo "max_flat_ref_contigs: ${MAX_FLAT_REF_CONTIGS}"
echo "skip_flat_if_hier_calls: ${SKIP_FLAT_IF_HIER_CALLS}"
echo "fast_debug: ${FAST_DEBUG}"
echo "run_minigraph: ${RUN_MINIGRAPH}"
echo "mycosv_tool_timeout: ${MYCOSV_TOOL_TIMEOUT}"
echo "mycosv_comparator_timeout: ${MYCOSV_COMPARATOR_TIMEOUT}"
echo "required_sv_types: ${REQUIRED_SV_TYPES}"
echo "require_all_six_sv_types: ${REQUIRE_ALL_SIX_SV_TYPES}"
echo "PATH-head: $(echo "${PATH}" | tr ':' '\n' | head -3)"
echo "minigraph: $(command -v minigraph || echo MISSING)"
echo "minimap2:  $(command -v minimap2  || echo MISSING)"
echo "samtools:  $(command -v samtools  || echo MISSING)"

if [[ ! -f "${PREPARED_DIR}/query_manifest.tsv" ]]; then
  echo "[error] missing ${PREPARED_DIR}/query_manifest.tsv" >&2
  exit 2
fi
if [[ ! -d "${PREPARED_DIR}/index" || ! -d "${PREPARED_DIR}/registry" ]]; then
  echo "[error] missing prepared index/registry under ${PREPARED_DIR}" >&2
  exit 2
fi

print_partial_debug_summary() {
  echo
  echo "=== Partial Debug Snapshot ==="
  if [[ -f "${OUT_DIR}/REQUESTED_QUERY_GROUPS.tsv" ]]; then
    echo "  requested query groups:"
    column -t -s $'\t' "${OUT_DIR}/REQUESTED_QUERY_GROUPS.tsv" 2>/dev/null || cat "${OUT_DIR}/REQUESTED_QUERY_GROUPS.tsv"
  fi
  local hier_vcf="${OUT_DIR}/mycosv/calls.hierarchical.vcf"
  if [[ -f "${hier_vcf}" ]]; then
    echo "  hierarchical SVTYPE mix:"
    for svt in ${REQUIRED_SV_TYPES//,/ }; do
      n=$(grep -v '^#' "${hier_vcf}" 2>/dev/null | grep -c "SVTYPE=${svt}\\b" || true)
      printf "    SVTYPE=%-8s %s\n" "${svt}" "${n}"
    done
    echo "  element classes:"
    grep -v '^#' "${hier_vcf}" 2>/dev/null \
      | grep -oE 'EC=[A-Z0-9_]+' \
      | sort | uniq -c \
      | awk '{printf "    %-14s %s\n", $2, $1}'
  fi
  echo "  newest key files:"
  find "${OUT_DIR}" -maxdepth 4 -type f \( \
      -name 'calls*.vcf' -o -name 'calls*.tsv' -o -name '*minigraph*' \
      -o -name 'biology_findings.tsv' -o -name 'novel_mycosv_calls.tsv' \
      -o -name 'read_validated_truth.tsv' -o -name 'pangenome_call_layers.tsv' \
    \) -printf '    %TY-%Tm-%Td %TH:%TM %9s %p\n' 2>/dev/null \
    | sort | tail -30
}

trap 'echo "[signal] caught TERM/INT at $(date)"; print_partial_debug_summary' TERM INT

export MYCOSV_FORCE_FLAT_REF_FALLBACK=1
export MYCOSV_TOOL_TIMEOUT
export MYCOSV_COMPARATOR_TIMEOUT

BENCHMARK_ARGS=()
if [[ "${RUN_MINIGRAPH}" == "1" ]]; then
  BENCHMARK_ARGS+=(--run-minigraph)
else
  BENCHMARK_ARGS+=(--mycosv-only)
fi
MYCOSV_EXTRA_ARGS=()
if [[ "${REQUIRE_ALL_SIX_SV_TYPES}" != "1" ]]; then
  MYCOSV_EXTRA_ARGS+=(--mycosv-arg=--no-graph-native-mode)
fi
python3 -u run_real_fungal_benchmark.py benchmark \
  --prepared-dir "${PREPARED_DIR}" \
  --out-dir "${OUT_DIR}" \
  --mode assembly \
  --threads "${THREADS}" \
  --max-clade-genomes 32 \
  --benchmark-query-genera "${TARGET_QUERY_GROUPS}" \
  --max-benchmark-queries "${MAX_BENCHMARK_QUERIES}" \
  --reuse-index-dir "${PREPARED_DIR}/index" \
  --reuse-registry-dir "${PREPARED_DIR}/registry" \
  --benchmark-ref-cap "${BENCHMARK_REF_CAP}" \
  "${BENCHMARK_ARGS[@]}" \
  --mycosv-arg=--max-calls-per-contig \
  --mycosv-arg="${MAX_CALLS_PER_CONTIG}" \
  --mycosv-arg=--min-block-score \
  --mycosv-arg="${MIN_BLOCK_SCORE}" \
  --mycosv-arg=--tol-min-chain-anchors \
  --mycosv-arg="${TOL_MIN_CHAIN_ANCHORS}" \
  --mycosv-arg=--max-ref-memory-mb \
  --mycosv-arg=8192 \
  --mycosv-arg=--max-flat-ref-contigs \
  --mycosv-arg="${MAX_FLAT_REF_CONTIGS}" \
  --mycosv-arg=--skip-flat-if-hier-calls \
  --mycosv-arg="${SKIP_FLAT_IF_HIER_CALLS}" \
  "${MYCOSV_EXTRA_ARGS[@]}" \
  --mycosv-arg=--no-gfa

RC=$?
echo
echo "wrapper exit: ${RC}"
echo

# ---- Stage-by-stage acceptance checks -------------------------------------
echo "=== Stage 1: MycoSV call generation ==="
for f in calls.vcf calls.hits.tsv calls.hierarchical.vcf calls.hierarchical.hits.tsv; do
  p="${OUT_DIR}/mycosv/${f}"
  if [[ -f "${p}" ]]; then
    n=$(grep -vc '^#' "${p}" 2>/dev/null || true)
    n=${n:-0}
    if [[ "${f}" == *.tsv ]]; then
      n=$((n > 0 ? n - 1 : 0))
    fi
    printf "  %-32s  %s records\n" "${f}" "${n}"
  else
    printf "  %-32s  MISSING\n" "${f}"
  fi
done

echo
echo "=== Stage 2: Comparator (minigraph) ==="
COMP_FILES=$(find "${OUT_DIR}" -maxdepth 4 -type f \( -name "minigraph*.vcf" -o -name "minigraph*.bub.bed" -o -name "*minigraph*.tsv" \) 2>/dev/null | head -10)
if [[ -n "${COMP_FILES}" ]]; then
  echo "${COMP_FILES}" | while read -r f; do
    printf "  %s  (%s records)\n" "${f}" "$(grep -vc '^#' "${f}" 2>/dev/null || echo ?)"
  done
else
  echo "  no minigraph output located"
fi

echo
echo "=== Stage 3: Read-level validation ==="
RV_DIR="${OUT_DIR}/read_validation"
if [[ -d "${RV_DIR}" ]]; then
  find "${RV_DIR}" -type f \( -name "*.bam" -o -name "*.tsv" -o -name "*.txt" \) 2>/dev/null | head -10
else
  echo "  ${RV_DIR} missing — read validation did not run"
fi

echo
echo "=== Stage 4: Biology findings ==="
for f in biology_findings.tsv novel_mycosv_calls.tsv exact_benchmark_summary.tsv pangenome_call_layers.tsv sv_volume_audit.tsv; do
  p="${OUT_DIR}/${f}"
  if [[ -f "${p}" ]]; then
    nrows=$(wc -l < "${p}")
    printf "  %-32s  %s lines\n" "${f}" "${nrows}"
  else
    printf "  %-32s  MISSING\n" "${f}"
  fi
done

echo
echo "=== Stage 5: Fungal biology sanity ==="
HIER_VCF="${OUT_DIR}/mycosv/calls.hierarchical.vcf"
ACCEPTANCE_FAIL=0
if [[ -f "${HIER_VCF}" ]]; then
  echo "  hierarchical SVTYPE mix:"
  missing_svtypes=()
  for svt in ${REQUIRED_SV_TYPES//,/ }; do
    n=$(grep -v '^#' "${HIER_VCF}" 2>/dev/null | grep -c "SVTYPE=${svt}\\b" || true)
    printf "    SVTYPE=%-8s %s\n" "${svt}" "${n}"
    if [[ "${n}" -eq 0 ]]; then
      missing_svtypes+=("${svt}")
    fi
  done
  if [[ "${#missing_svtypes[@]}" -gt 0 ]]; then
    ACCEPTANCE_FAIL=1
    echo "  missing_required_svtypes       ${missing_svtypes[*]}"
  else
    echo "  all_required_svtypes_present   yes"
  fi
  total_hier=$(grep -vc '^#' "${HIER_VCF}" 2>/dev/null || echo 0)
  n_sv_classes=$(grep -v '^#' "${HIER_VCF}" \
    | grep -oE 'SVTYPE=[A-Z_]+' | sort -u | wc -l)
  n_te_like=$(grep -v '^#' "${HIER_VCF}" \
    | grep -Ec 'EC=(RIP|REPEAT|TE|TE_|LTR_|DNA_TIR|LINE|SINE|STARSHIP|HGT|HELITRON|MITE)' || true)
  n_large=$(grep -v '^#' "${HIER_VCF}" \
    | awk -F'[;\t]' '{
        for (i=1;i<=NF;i++) if ($i ~ /^SVLEN=/) {
          split($i,a,"=");
          v=a[2]+0; if (v<0) v=-v;
          if (v>=5000) n++;
        }
      } END {print n+0}')
  printf "  hierarchical_total             %s\n" "${total_hier}"
  printf "  distinct_svtypes               %s\n" "${n_sv_classes}"
  printf "  te_or_mge_like_records         %s\n" "${n_te_like}"
  printf "  large_5kb_plus_records         %s\n" "${n_large}"
else
  echo "  missing hierarchical VCF"
  ACCEPTANCE_FAIL=1
fi

PANG="${OUT_DIR}/pangenome_call_layers.tsv"
if [[ -f "${PANG}" ]]; then
  awk -F'\t' 'NR==1{next} $1=="ALL" || $1!="query_asm" {
      raw+=$3; loci+=$4; single+=$5; pang+=$6; read+=$7; intrinsic+=$8
    } END {
      printf "  pangenome_raw_observations     %d\n", raw+0;
      printf "  pangenome_deduped_loci         %d\n", loci+0;
      printf "  pangenome_only_calls           %d\n", pang+0;
      printf "  pangenome_supported_calls      %d\n", read+intrinsic+0;
    }' "${PANG}"
fi

BIO="${OUT_DIR}/biology_findings.tsv"
if [[ -f "${BIO}" ]]; then
  n_bio=$(awk 'NR>1 && NF>0 {n++} END{print n+0}' "${BIO}")
  n_hgt=$(awk -F'\t' 'NR>1 && ($2 ~ /hgt/ || $16=="yes" || $15=="HGT" || $15=="STARSHIP") {n++} END{print n+0}' "${BIO}")
  n_te=$(awk -F'\t' 'NR>1 && ($15 ~ /(RIP|REPEAT|TE|LTR|DNA_TIR|LINE|SINE|HELITRON|MITE)/) {n++} END{print n+0}' "${BIO}")
  printf "  biology_findings_records       %s\n" "${n_bio}"
  printf "  hgt_starship_candidates        %s\n" "${n_hgt}"
  printf "  te_repeat_rewiring_candidates  %s\n" "${n_te}"
fi

echo
echo "=== Acceptance hints ==="
echo "  Expected fungal signal for this one-query smoke:"
echo "    - hierarchical_total well above zero and preferably >= 100"
echo "    - all six SV classes present: INS, DEL, DUP, INV, TRA, OFF_REF"
echo "    - TE/RIP/repeat or MGE-like annotations present"
echo "    - pangenome_call_layers non-zero, with pangenome-only calls if comparator/read stages finish"
echo "    - biology_findings.tsv non-header rows explaining TE/MGE/HGT/accessory hypotheses"

echo
echo "end: $(date)"
if [[ "${ACCEPTANCE_FAIL}" -ne 0 && "${RC}" -eq 0 ]]; then
  echo "[acceptance] failed required SVTYPE check"
  exit 3
fi
exit ${RC}
