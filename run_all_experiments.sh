#!/usr/bin/env bash
# run_all_experiments.sh
# Master script to run all experiments (simulated benchmark and real fungal data).
#
# Usage:
#   bash run_all_experiments.sh                    # Run all experiments
#   bash run_all_experiments.sh --simulated        # Simulated benchmark only
#   bash run_all_experiments.sh --real             # Real panels only (per-species NCBI + ENA)
#   bash run_all_experiments.sh --million-real     # Only the real million-scale NCBI index build
#
# Environment overrides for --million-real:
#   MILLION_REAL_MAX_ASSEMBLIES=1000    # how many NCBI fungal assemblies to download
#   MILLION_REAL_TARGET_CENTROIDS=1000000  # total centroids after decoy padding (0 = no padding)
#
# IMPORTANT: we intentionally do NOT use `set -e` at the script level for the
# real-data stage. A single panel/network/tool failure must not abort the rest
# of the matrix — each panel is guarded with `|| true` and its log is kept for
# inspection. Fatal errors per stage are still reported in the summary.

set -u
set -o pipefail

WORK_DIR="/mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/AMF/scale"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
EXPERIMENT_TYPE="${1:-all}"  # all, simulated, real, million-real

# Shared download cache: FASTA/GFF/FASTQ files and metadata fetched here are
# reused across timestamped experiment runs.
DATA_CACHE_DIR="${DATA_CACHE_DIR:-${WORK_DIR}/data_cache}"
mkdir -p "${DATA_CACHE_DIR}"

# Strip leading dashes for convenience: --simulated and simulated are both accepted.
EXPERIMENT_TYPE="${EXPERIMENT_TYPE#--}"

# How many real NCBI fungal assemblies to download when building the real
# million-scale routing index. A full 10000-assembly pull is appropriate for a
# dedicated --million-real run, but it does not fit comfortably inside the
# comprehensive all-panels matrix on a 10-24h Slurm allocation. The all-mode
# default still pads the routing store to one million centroids via decoys.
if [[ -z "${MILLION_REAL_MAX_ASSEMBLIES:-}" ]]; then
  if [[ "${EXPERIMENT_TYPE}" == "million-real" ]]; then
    MILLION_REAL_MAX_ASSEMBLIES=10000
  else
    MILLION_REAL_MAX_ASSEMBLIES=2000
  fi
fi
MILLION_REAL_TARGET_CENTROIDS="${MILLION_REAL_TARGET_CENTROIDS:-1000000}"
# How many of the downloaded assemblies to hold out as MycoSV-only benchmark
# queries (excluded from the index, then run end-to-end through indexing,
# alignment, SV calling, TE classification, biology candidates, and
# visualization in step 2b). 5 keeps wall time bounded on large indexes;
# raise on long-walltime nodes.
MILLION_REAL_QUERIES="${MILLION_REAL_QUERIES:-5}"
# Per-clade RAM cap on the MycoSV binary's hierarchical graph build; the
# default fits a 12 GiB cgroup but can be raised on bigger nodes.
MILLION_REAL_MAX_CLADE_GENOMES="${MILLION_REAL_MAX_CLADE_GENOMES:-8}"

# Per-query memory caps for the MycoSV binary in the million-real bench.
# Without a per-thread cache cap, the SingleRefMemCache (suffix arrays for
# every ref contig touched by a query) grows unbounded; with --threads
# parallel queries each holding their own cache the binary previously died
# with std::bad_alloc on multi-Gbp routed sub-graphs.  Default 4096 MB per
# thread fits comfortably in a 128 GiB SLURM job at 16 threads (~64 GiB
# total cache budget); shrink on tight cgroups or raise on bigger nodes.
MILLION_REAL_SA_MAX_CONTIG_MB="${MILLION_REAL_SA_MAX_CONTIG_MB:-25}"
MILLION_REAL_SINGLE_REF_CACHE_MB="${MILLION_REAL_SINGLE_REF_CACHE_MB:-4096}"
MILLION_REAL_MAX_REF_MEMORY_MB="${MILLION_REAL_MAX_REF_MEMORY_MB:-1024}"
MILLION_REAL_MAX_ASSEMBLY_QUERY_CONTIGS="${MILLION_REAL_MAX_ASSEMBLY_QUERY_CONTIGS:-5000}"
MILLION_REAL_MAX_ASSEMBLY_QUERY_BP="${MILLION_REAL_MAX_ASSEMBLY_QUERY_BP:-350000000}"

# Effective cgroup memory limit, when Slurm/cgroup v2 exposes one. This is
# more useful than host RAM on shared HPC nodes: the failed 14535674 run had a
# 64 GiB cgroup despite a larger physical node.
detect_cgroup_memory_gib() {
  local line rel cur candidate value depth i
  local IFS='/'
  local -a parts
  [[ -r /proc/self/cgroup ]] || return 0
  line="$(head -n 1 /proc/self/cgroup)"
  [[ "${line}" == *"::"* ]] || return 0
  rel="${line#*::}"
  rel="${rel#/}"
  cur="/sys/fs/cgroup"
  read -r -a parts <<< "${rel}"
  for (( depth=${#parts[@]}; depth>=0; depth-- )); do
    candidate="${cur}"
    for (( i=0; i<depth; i++ )); do
      [[ -n "${parts[$i]}" ]] && candidate+="/${parts[$i]}"
    done
    candidate+="/memory.max"
    [[ -r "${candidate}" ]] || continue
    value="$(<"${candidate}")"
    if [[ -n "${value}" && "${value}" != "max" ]]; then
      echo $(( value / 1024 / 1024 / 1024 ))
      return 0
    fi
  done
}
CGROUP_MEMORY_GIB="$(detect_cgroup_memory_gib)"

# Worker threads passed to every tool that accepts a thread count:
#   MycoSV binary (--threads / --tol-index-threads), minimap2 (-t),
#   samtools sort (-@), sniffles2 (--threads), cuteSV (--threads),
#   Delly (OMP_NUM_THREADS), Manta (-j), minigraph (-t), pggb (-t).
if [[ -z "${THREADS+x}" ]]; then
  THREADS="${SLURM_CPUS_PER_TASK:-32}"
  if [[ -n "${CGROUP_MEMORY_GIB}" && "${CGROUP_MEMORY_GIB}" -le 80 && "${THREADS}" -gt 16 ]]; then
    THREADS=16
  fi
fi

# AMF / repeat-rich assembly mode can still spike memory inside whole-contig
# matching. Run those panel/mode invocations with a much smaller thread budget
# by default; override if you are on a genuinely larger cgroup.
REAL_HEAVY_ASSEMBLY_THREADS="${REAL_HEAVY_ASSEMBLY_THREADS:-1}"
REAL_HEAVY_SA_MAX_CONTIG_MB="${REAL_HEAVY_SA_MAX_CONTIG_MB:-25}"
REAL_HEAVY_SINGLE_REF_CACHE_MB="${REAL_HEAVY_SINGLE_REF_CACHE_MB:-1024}"
REAL_HEAVY_MAX_REF_MEMORY_MB="${REAL_HEAVY_MAX_REF_MEMORY_MB:-1024}"
REAL_HEAVY_MAX_ASSEMBLY_QUERY_CONTIGS="${REAL_HEAVY_MAX_ASSEMBLY_QUERY_CONTIGS:-5000}"
REAL_HEAVY_MAX_ASSEMBLY_QUERY_BP="${REAL_HEAVY_MAX_ASSEMBLY_QUERY_BP:-350000000}"
# AMF (Gigaspora, Rhizophagus) assemblies are inherently fragmented — public
# Gigaspora drafts ship with 50–100k+ contigs because the genomes are large
# (~700 Mbp) and repeat-rich. The default 5000-contig cap was rejecting every
# AMF query, leaving the panel with a single Rhizophagus assembly. Use a much
# higher AMF-specific cap so the panel actually exercises >1 query.
REAL_HEAVY_MAX_ASSEMBLY_QUERY_CONTIGS_AMF="${REAL_HEAVY_MAX_ASSEMBLY_QUERY_CONTIGS_AMF:-200000}"
REAL_HEAVY_MAX_ASSEMBLY_QUERY_BP_AMF="${REAL_HEAVY_MAX_ASSEMBLY_QUERY_BP_AMF:-1500000000}"

# Timeout used inside run_real_fungal_benchmark.py for individual external
# tool calls. Keep it at least as large as the assembly benchmark cap so AMF
# MycoSV runs are governed by the stage timeout, not an older 2h hard stop.
export MYCOSV_TOOL_TIMEOUT="${MYCOSV_TOOL_TIMEOUT:-14400}"

# Stage-level wall-clock budgets. These are fail-forward guards: a stage that
# exceeds its cap is marked failed/timeout, then the matrix moves on to the
# next independent stage where possible.
SIM_BENCHMARK_TIMEOUT="${SIM_BENCHMARK_TIMEOUT:-3h}"
REAL_PREPARE_TIMEOUT="${REAL_PREPARE_TIMEOUT:-2h}"
REPORT_TIMEOUT="${REPORT_TIMEOUT:-30m}"

# Per-panel/per-mode wall-clock budget. A single benchmark invocation that
# overruns this limit is killed and the next mode/panel still runs — without
# it, one slow panel (typically AMF assembly mode where cactus chews on
# ~1 Gbp Rhizophagus genomes) consumes the SLURM time budget and starves
# every subsequent panel in the matrix. Override via env var when running
# on long-walltime partitions. Set to 0 to disable the cap.
BENCHMARK_TIMEOUT_ASSEMBLY="${BENCHMARK_TIMEOUT_ASSEMBLY:-3h}"
BENCHMARK_TIMEOUT_READS="${BENCHMARK_TIMEOUT_READS:-2h}"

# Per-panel overrides. Puccinia genomes in te_rich_pathogen are ~80 Mbp and
# ~85% TE content; even with cactus/pggb/anchorwave skipped, the surviving
# minigraph + svim_asm + MycoSV combination consistently spills past 3h on a
# single thread. Five hours empirically clears the matrix without reaching
# the SLURM walltime.
BENCHMARK_TIMEOUT_TE_RICH_PATHOGEN_ASSEMBLY="${BENCHMARK_TIMEOUT_TE_RICH_PATHOGEN_ASSEMBLY:-5h}"

# Bound public-read comparator inputs. MycoSV caps read consumption internally,
# but tools such as SVIM/Sniffles/cuteSV/Delly/Manta align the FASTQ path they
# are given. These caps keep public ENA runs from dominating wall time.
MAX_COMPARATOR_SHORT_READS="${MAX_COMPARATOR_SHORT_READS:-150000}"
MAX_COMPARATOR_LONG_READS="${MAX_COMPARATOR_LONG_READS:-20000}"

# Comma-separated list of panel names whose assembly-mode benchmark should
# skip cactus / pggb / anchorwave. These three are the heaviest comparators
# on multi-Gbp inputs (AMF panel, te_rich_pathogen for Puccinia ~80 Mbp,
# two_speed_pathogen for Fusarium oxysporum ~60 Mbp + Zymoseptoria tritici);
# skipping them cuts wall time enough that the matrix completes within
# typical SLURM time limits while still leaving syri / minigraph / svim_asm
# as truth comparators. two_speed_pathogen was added 2026-05-06 after a SLURM
# run timed out mid-cactus-on-Fusarium and starved the panel's short-reads
# and long-reads modes. Override or extend via env var.
HEAVY_COMPARATOR_SKIP_PANELS="${HEAVY_COMPARATOR_SKIP_PANELS:-amf_large,te_rich_pathogen,two_speed_pathogen}"

# Real-data panels to run by default. compact_yeast + amf_large cover the
# current baseline; te_rich_pathogen adds TE-rich full-SV stress; and
# two_speed_pathogen adds plant pathogen / effector-region biology. Other
# PANEL_PRESETS, such as cross_phylum_hgt, remain available via REAL_PANELS=...
REAL_PANELS_DEFAULT="${REAL_PANELS_DEFAULT:-compact_yeast,amf_large,te_rich_pathogen,two_speed_pathogen}"

# Long-read platform used for both simulated and real experiments.
#   hifi    PacBio HiFi CCS (Revio / Sequel IIe) — 15 kb reads, ≥Q20 accuracy.
#           minimap2 map-hifi → sniffles2 / cuteSV (HiFi params) / SVIM.
#   ont-r10 ONT R10.4.1 standard simplex — 10 kb, ~Q20 on PromethION/GridION.
#           minimap2 map-ont → sniffles2 --long-read-model ont_r10_q20 (v2.2+).
#           WhatsHap phase+haplotag: applicable for dikaryotic fungi such as
#           Puccinia spp., Leptosphaeria maculans, Zymoseptoria tritici.
#   ont-r9  ONT R9.4.1 legacy — 8 kb, ~Q15.  Still prevalent in public ENA data.
LONG_READ_PLATFORM="${LONG_READ_PLATFORM:-ont-r10}"

# Create experiment directories
SIM_DIR="${WORK_DIR}/experiments/simulated/${TIMESTAMP}"
REAL_DIR="${WORK_DIR}/experiments/real_data/${TIMESTAMP}"
MILLION_REAL_DIR="${WORK_DIR}/experiments/million_real/${TIMESTAMP}"

mkdir -p "${SIM_DIR}" "${REAL_DIR}" "${MILLION_REAL_DIR}"

cd "${WORK_DIR}"

# Color codes for output
GREEN='\033[0;32m'
RED='\033[0;31m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Track stage outcomes so we can still print a meaningful summary even if
# individual panels or modes fail.
FAILED_STAGES=()
SUCCESS_STAGES=()

mark_success() { SUCCESS_STAGES+=("$1"); }
mark_failure() { FAILED_STAGES+=("$1"); echo -e "${RED}✗ $1${NC}"; }

echo -e "${BLUE}========================================"
echo "MycoSV Comprehensive Experiments"
echo "========================================${NC}"
echo "Start time: $(date)"
echo "Timestamp: ${TIMESTAMP}"
echo "Working directory: ${WORK_DIR}"
echo "Experiment type: ${EXPERIMENT_TYPE}"
echo ""

# Scenario sets that guarantee all 5 SV types (INS, DEL, DUP, INV, TRA).
# compact_yeast:              DEL + INS
# two_speed_pathogen_extreme: INV + TRA + INS
# arbuscular_mf:              DUP + INS
SIM_SCENARIO_SET="compact_yeast,two_speed_pathogen_extreme,arbuscular_mf"

# ============================================================================
# 1. SIMULATED BENCHMARK (precision/recall at million scale)
#
# Per-scenario query budget controls runtime. With QUERIES_PER_SCENARIO=20 and
# 3 scenarios → 60 query genomes (+ refs). At n_contigs=10 that is ~600 truth
# SVs per mode — fast turnaround for testing. Override with the env var
# SIM_QUERIES_PER_SCENARIO to scale up (e.g. 200 for ~6000 SVs).
# ============================================================================

SIM_QUERIES_PER_SCENARIO="${SIM_QUERIES_PER_SCENARIO:-20}"

if [[ "$EXPERIMENT_TYPE" == "all" || "$EXPERIMENT_TYPE" == "simulated" ]]; then
  N_SCENARIOS=$(awk -F',' '{print NF}' <<< "${SIM_SCENARIO_SET}")
  TOTAL_QUERIES=$(( SIM_QUERIES_PER_SCENARIO * N_SCENARIOS ))
  echo -e "${YELLOW}[1/4] Running simulated benchmark (million-scale routing)...${NC}"
  echo "      Modes: assembly, short-reads, long-reads"
  echo "      Scenarios: ${SIM_SCENARIO_SET} (covers INS/DEL/DUP/INV/TRA)"
  echo "      Queries per scenario: ${SIM_QUERIES_PER_SCENARIO} → ${TOTAL_QUERIES} query genomes (+ 3 refs), 10 contigs each"
  echo "      Output: ${SIM_DIR}/benchmarks"
  mkdir -p "${SIM_DIR}/benchmarks"

  sim_cmd=(python3 -u run_million_mode_query_benchmark.py
      --out-dir "${SIM_DIR}/benchmarks" \
      --modes assembly,short-reads,long-reads \
      --scenario-set "${SIM_SCENARIO_SET}" \
      --n-centroids 1000000 \
      --n-reps 3 \
      --n-contigs 10 \
      --total-len 200000 \
      --seed 42 \
      --queries-per-scenario "${SIM_QUERIES_PER_SCENARIO}" \
      --min-contig-bp 12000 \
      --long-read-platform "${LONG_READ_PLATFORM}" \
      --threads "${THREADS}")
  if [[ -n "${SIM_BENCHMARK_TIMEOUT}" && "${SIM_BENCHMARK_TIMEOUT}" != "0" ]] && command -v timeout >/dev/null; then
    sim_cmd=(timeout --signal=TERM --kill-after=60 "${SIM_BENCHMARK_TIMEOUT}" "${sim_cmd[@]}")
  fi
  if "${sim_cmd[@]}" 2>&1 | tee "${SIM_DIR}/benchmarks/benchmark.log"; then
    mark_success "simulated.pr_metrics_benchmark"
    echo -e "${GREEN}✓ Simulated benchmark complete${NC}"
  else
    rc=${PIPESTATUS[0]}
    if [[ "${rc}" == "124" || "${rc}" == "137" ]]; then
      mark_failure "simulated.pr_metrics_benchmark.timeout(${SIM_BENCHMARK_TIMEOUT})"
    else
      mark_failure "simulated.pr_metrics_benchmark"
    fi
  fi
  echo ""
fi

# ============================================================================
# 2. MILLION-SCALE *REAL* FUNGAL INDEX (downloads NCBI assemblies)
# ============================================================================

if [[ "$EXPERIMENT_TYPE" == "all" || "$EXPERIMENT_TYPE" == "million-real" ]]; then
  echo -e "${YELLOW}[2/4] Building real million-scale fungal routing index + MycoSV-only benchmark...${NC}"
  echo "      Downloading up to ${MILLION_REAL_MAX_ASSEMBLIES} NCBI GenBank assemblies (contig level or better)"
  echo "      Target centroids (real+decoys): ${MILLION_REAL_TARGET_CENTROIDS}"
  echo "      Holding out ${MILLION_REAL_QUERIES} assemblies as MycoSV-only benchmark queries"
  echo "      Output: ${MILLION_REAL_DIR}"

  # Wall-clock caps per sub-step so a slow NCBI download burst or runaway
  # binary cannot starve the rest of the matrix.  Override / disable with
  # MILLION_REAL_PREPARE_TIMEOUT=0 / MILLION_REAL_BENCH_TIMEOUT=0.
  MILLION_REAL_PREPARE_TIMEOUT="${MILLION_REAL_PREPARE_TIMEOUT:-${MILLION_REAL_TIMEOUT:-8h}}"
  MILLION_REAL_BENCH_TIMEOUT="${MILLION_REAL_BENCH_TIMEOUT:-4h}"

  prepare_million_succeeded=0
  # Reads-mode coverage in the million-real bench. When 1, prepare-million-real
  # also resolves public ENA reads for each held-out query species so 2b can
  # exercise short-reads (Illumina) and long-reads (PacBio HiFi / ONT R10.4)
  # in addition to assembly mode.  Default ON so benchmark_short-reads/ and
  # benchmark_long-reads/ get populated; set MILLION_REAL_INCLUDE_READS=0 to
  # skip the ENA hops on time-constrained runs.
  MILLION_REAL_INCLUDE_READS="${MILLION_REAL_INCLUDE_READS:-1}"
  MILLION_REAL_READ_MODES="${MILLION_REAL_READ_MODES:-both}"
  MILLION_REAL_READ_RUNS_PER_QUERY="${MILLION_REAL_READ_RUNS_PER_QUERY:-1}"
  MILLION_REAL_NCBI_SOURCE="${MILLION_REAL_NCBI_SOURCE:-ncbi-best}"

  # ── 2a) prepare: download assemblies, build the routing index, hold out
  # ────── MILLION_REAL_QUERIES assemblies for the MycoSV-only benchmark.
  # python3 -u keeps stdout line-buffered under tee so progress is visible.
  prepare_cmd=(python3 -u run_real_fungal_benchmark.py prepare-million-real
      --out-dir "${MILLION_REAL_DIR}"
      --source "${MILLION_REAL_NCBI_SOURCE}"
      --max-assemblies "${MILLION_REAL_MAX_ASSEMBLIES}"
      --target-centroids "${MILLION_REAL_TARGET_CENTROIDS}"
      --min-assembly-level contig
      --threads "${THREADS}"
      --seed 42
      --max-clade-genomes "${MILLION_REAL_MAX_CLADE_GENOMES}"
      --million-real-queries "${MILLION_REAL_QUERIES}"
      --million-real-read-modes "${MILLION_REAL_READ_MODES}"
      --million-real-read-runs-per-query "${MILLION_REAL_READ_RUNS_PER_QUERY}"
      --data-cache-dir "${DATA_CACHE_DIR}")
  if [[ "${MILLION_REAL_INCLUDE_READS}" == "1" ]]; then
    prepare_cmd+=(--million-real-include-reads)
  fi
  if [[ -n "${MILLION_REAL_PREPARE_TIMEOUT}" && "${MILLION_REAL_PREPARE_TIMEOUT}" != "0" ]] && command -v timeout >/dev/null; then
    prepare_cmd=(timeout --signal=TERM --kill-after=60 "${MILLION_REAL_PREPARE_TIMEOUT}" "${prepare_cmd[@]}")
  fi
  if "${prepare_cmd[@]}" 2>&1 | tee "${MILLION_REAL_DIR}/prepare_million_real.log"; then
    mark_success "million_real.index_build"
    echo -e "${GREEN}✓ Real million-scale index ready${NC}"
    prepare_million_succeeded=1
  else
    rc=${PIPESTATUS[0]}
    if [[ "${rc}" == "124" || "${rc}" == "137" ]]; then
      mark_failure "million_real.index_build.timeout(${MILLION_REAL_PREPARE_TIMEOUT})"
      echo -e "${RED}✗ million-real prepare exceeded ${MILLION_REAL_PREPARE_TIMEOUT} — skipping 2b${NC}"
    else
      mark_failure "million_real.index_build"
    fi
  fi
  echo ""

  # ── 2b) Benchmark on the held-out queries: indexing was done in 2a, so
  # ────── reuse that index (no rebuild). Per-mode comparator subset chosen
  # ────── for bounded wall time (assembly: minigraph+svim_asm, long-reads:
  # ────── sniffles+cutesv, short-reads: delly). Set MILLION_REAL_MYCOSV_ONLY=1
  # ────── to skip every comparator and only exercise indexing → alignment →
  # ────── SV calling → TE classification → biology candidates →
  # ────── biology_findings.tsv → novel_mycosv_calls.tsv. Visualization
  # ────── (step 4) picks all of these up via MILLION_REAL_DIR scan.
  if [[ "${prepare_million_succeeded}" == "1" \
        && -f "${MILLION_REAL_DIR}/query_manifest.tsv" \
        && $(wc -l < "${MILLION_REAL_DIR}/query_manifest.tsv") -gt 1 ]]; then
    # Read-level validation needs a single per-query benchmark_ref_fasta;
    # in the million-real flow that ref is a sibling-genus assembly chosen
    # by prepare-million-real. It is ON by default so million-real produces
    # read_validated_truth.tsv / benchmark_summary.json support counts; set
    # MILLION_REAL_VALIDATE_WITH_READS=0 to skip on very tight walltime runs.
    val_flag="--validate-with-reads"
    if [[ "${MILLION_REAL_VALIDATE_WITH_READS:-1}" == "0" ]]; then
      val_flag="--no-validate-with-reads"
    fi

    # ── Modes to bench. assembly is always on; reads modes are added if
    # ────── prepare-million-real materialised matching query rows. The
    # ────── benchmark step itself filters by mode and writes a
    # ────── NO_QUERIES_FOR_MODE.txt marker when none exist, so it's safe
    # ────── to invoke unconditionally — but skipping the call avoids a
    # ────── spurious comparator-pre-flight burst per missing mode.
    declare -A million_real_modes_present
    million_real_modes_present[assembly]=0
    million_real_modes_present[short-reads]=0
    million_real_modes_present[long-reads]=0
    while IFS=$'\t' read -r _qa qmode _rest; do
      case "${qmode}" in
        assembly|short-reads|long-reads)
          million_real_modes_present[${qmode}]=$(( million_real_modes_present[${qmode}] + 1 ))
          ;;
      esac
    done < <(tail -n +2 "${MILLION_REAL_DIR}/query_manifest.tsv")

    for million_real_mode in assembly short-reads long-reads; do
      if [[ "${million_real_modes_present[${million_real_mode}]}" -eq 0 ]]; then
        if [[ "${million_real_mode}" != "assembly" ]]; then
          echo "      [skip] mode=${million_real_mode}: no query rows in manifest"
          echo "             (re-run prep with MILLION_REAL_INCLUDE_READS=1 to materialise reads queries)"
        fi
        continue
      fi
      echo -e "${YELLOW}[2b/4] Running MycoSV-only benchmark mode=${million_real_mode}...${NC}"
      bench_dir="${MILLION_REAL_DIR}/benchmark_${million_real_mode}"
      mkdir -p "${bench_dir}"
      # Per-mode read-level validation support: assembly truth is
      # high-confidence so 1 supporting read is enough; reads-mode truth
      # is read-derived already so we want a stricter ≥3 to reject
      # alignment artefacts.  Match the per-panel real-data defaults.
      if [[ "${million_real_mode}" == "assembly" ]]; then
        million_real_readval_support="${MILLION_REAL_READ_VALIDATION_MIN_SUPPORT:-1}"
      else
        million_real_readval_support="${MILLION_REAL_READ_VALIDATION_MIN_SUPPORT_READS:-3}"
      fi
      # Per-mode comparator selection: keep wall time bounded by skipping
      # the heaviest tool in each mode while still producing real F1 numbers
      # in exact_benchmark_summary.tsv. Override with MILLION_REAL_MYCOSV_ONLY=1
      # to fall back to the previous mycosv-only behaviour.
      million_real_comparator_flags=()
      if [[ "${MILLION_REAL_MYCOSV_ONLY:-0}" == "1" ]]; then
        million_real_comparator_flags=(--mycosv-only)
      else
        case "${million_real_mode}" in
          assembly)    million_real_comparator_flags=(--run-minigraph --run-svim-asm) ;;
          long-reads)  million_real_comparator_flags=(--run-sniffles --run-cutesv) ;;
          short-reads) million_real_comparator_flags=(--run-delly) ;;
        esac
      fi
      bench_cmd=(python3 -u run_real_fungal_benchmark.py benchmark
          --prepared-dir "${MILLION_REAL_DIR}"
          --out-dir "${bench_dir}"
          --mode "${million_real_mode}"
          --threads "${THREADS}"
          --max-clade-genomes "${MILLION_REAL_MAX_CLADE_GENOMES}"
          "${million_real_comparator_flags[@]}"
          --reuse-index-dir "${MILLION_REAL_DIR}/index"
          --reuse-registry-dir "${MILLION_REAL_DIR}/registry"
          --read-validation-min-support "${million_real_readval_support}"
          --mycosv-arg=--sa-max-contig-mb
          --mycosv-arg="${MILLION_REAL_SA_MAX_CONTIG_MB}"
          --mycosv-arg=--single-ref-cache-mb
          --mycosv-arg="${MILLION_REAL_SINGLE_REF_CACHE_MB}"
          --mycosv-arg=--max-ref-memory-mb
          --mycosv-arg="${MILLION_REAL_MAX_REF_MEMORY_MB}"
          --mycosv-arg=--no-flat-ref-fallback
          --mycosv-arg=--no-gfa
          --mycosv-arg=--tol-min-chain-anchors
          --mycosv-arg=2
          "${val_flag}")
      if [[ "${million_real_mode}" == "assembly" ]]; then
        bench_cmd+=(
          --max-assembly-query-contigs "${MILLION_REAL_MAX_ASSEMBLY_QUERY_CONTIGS}"
          --max-assembly-query-bp "${MILLION_REAL_MAX_ASSEMBLY_QUERY_BP}"
        )
      fi
      if [[ -n "${MILLION_REAL_BENCH_TIMEOUT}" && "${MILLION_REAL_BENCH_TIMEOUT}" != "0" ]] && command -v timeout >/dev/null; then
        bench_cmd=(timeout --signal=TERM --kill-after=60 "${MILLION_REAL_BENCH_TIMEOUT}" "${bench_cmd[@]}")
      fi
      if "${bench_cmd[@]}" 2>&1 | tee "${bench_dir}/benchmark.log"; then
        mark_success "million_real.mycosv_benchmark.${million_real_mode}"
        echo -e "${GREEN}✓ MycoSV-only million-real benchmark complete (mode=${million_real_mode})${NC}"
      else
        rc=${PIPESTATUS[0]}
        if [[ "${rc}" == "124" || "${rc}" == "137" ]]; then
          mark_failure "million_real.mycosv_benchmark.${million_real_mode}.timeout(${MILLION_REAL_BENCH_TIMEOUT})"
        else
          mark_failure "million_real.mycosv_benchmark.${million_real_mode}"
        fi
      fi
      echo ""
    done
  elif [[ "${prepare_million_succeeded}" == "1" ]]; then
    echo "      [skip] no held-out queries written (MILLION_REAL_QUERIES=${MILLION_REAL_QUERIES})"
  fi
fi

# ============================================================================
# 3. REAL FUNGAL DATA BENCHMARKS
# ============================================================================
#
# Each panel is fully isolated: a failure in one panel (download error,
# missing tool, empty selection) does not stop the rest.
#
# State-of-the-art comparators are enabled per mode:
#   assembly    : SyRI, minigraph+gfatools, PGGB, Minigraph-Cactus,
#                 SVIM-asm, AnchorWave
#   short-reads : Delly, Manta
#   long-reads  : SVIM, Sniffles, cuteSV
# ============================================================================

if [[ "$EXPERIMENT_TYPE" == "all" || "$EXPERIMENT_TYPE" == "real" ]]; then
  IFS=',' read -r -a PANELS <<< "${REAL_PANELS:-${REAL_PANELS_DEFAULT}}"
  IFS=',' read -r -a REAL_MODES <<< "${REAL_BENCHMARK_MODES:-assembly,short-reads,long-reads}"

  echo -e "${YELLOW}[3/4] Running real fungal data benchmarks...${NC}"
  echo "      Panels: ${PANELS[*]}"
  echo "      Modes: ${REAL_MODES[*]}"
  echo "      Output: ${REAL_DIR}"

  for panel in "${PANELS[@]}"; do
    echo ""
    echo "  ── Processing panel: ${panel} ──"
    PANEL_DIR="${REAL_DIR}/${panel}"
    mkdir -p "${PANEL_DIR}"

    echo "    - Preparing real data (mixed: assembly + short-reads + long-reads)..."
    # Reference and per-species caps are sized for a 12 GiB cgroup. The
    # MycoSV binary loads multiple references into RAM when building the
    # routing index; with REAL_MAX_REF_DOWNLOADS=6 / max-assemblies-per-species=3
    # the index build peaks well below 12 GiB. Override with env vars when
    # running on a larger node (e.g. REAL_MAX_REF_DOWNLOADS=20 on >=24 GiB).
    REAL_MAX_REF_DOWNLOADS="${REAL_MAX_REF_DOWNLOADS:-6}"
    REAL_MAX_ASMS_PER_SPECIES="${REAL_MAX_ASMS_PER_SPECIES:-3}"
    REAL_QUERIES_PER_SPECIES="${REAL_QUERIES_PER_SPECIES:-3}"
    REAL_MAX_QUERY_DOWNLOADS="${REAL_MAX_QUERY_DOWNLOADS:-6}"
    REAL_READ_ACCESSIONS_PER_SPECIES="${REAL_READ_ACCESSIONS_PER_SPECIES:-1}"
    REAL_NCBI_SOURCE="${REAL_NCBI_SOURCE:-ncbi-best}"
    prepare_panel_cmd=(python3 -u run_real_fungal_benchmark.py prepare
        --out-dir "${PANEL_DIR}/prepared" \
        --panel "${panel}" \
        --source "${REAL_NCBI_SOURCE}" \
        --max-assemblies-per-species "${REAL_MAX_ASMS_PER_SPECIES}" \
        --querys-per-species "${REAL_QUERIES_PER_SPECIES}" \
        --max-ref-downloads "${REAL_MAX_REF_DOWNLOADS}" \
        --max-query-downloads "${REAL_MAX_QUERY_DOWNLOADS}" \
        --query-mode mixed \
        --read-accessions-per-species "${REAL_READ_ACCESSIONS_PER_SPECIES}" \
        --allow-no-queries \
        --data-cache-dir "${DATA_CACHE_DIR}")
    if [[ -n "${REAL_PREPARE_TIMEOUT}" && "${REAL_PREPARE_TIMEOUT}" != "0" ]] && command -v timeout >/dev/null; then
      prepare_panel_cmd=(timeout --signal=TERM --kill-after=60 "${REAL_PREPARE_TIMEOUT}" "${prepare_panel_cmd[@]}")
    fi
    if "${prepare_panel_cmd[@]}" 2>&1 | tee "${PANEL_DIR}/prepare.log"; then
      mark_success "real.${panel}.prepare"
    else
      rc=${PIPESTATUS[0]}
      if [[ "${rc}" == "124" || "${rc}" == "137" ]]; then
        mark_failure "real.${panel}.prepare.timeout(${REAL_PREPARE_TIMEOUT})"
      else
        mark_failure "real.${panel}.prepare"
      fi
      echo "    - Skipping benchmarks for ${panel} (prepare failed)"
      continue
    fi

    if [[ ! -f "${PANEL_DIR}/prepared/selected_catalog.tsv" \
       && ! -f "${PANEL_DIR}/prepared/reference_catalog.tsv" ]]; then
      echo "    - Skipping benchmarks for ${panel} (no catalog produced)"
      mark_failure "real.${panel}.no_catalog"
      continue
    fi
    if [[ ! -f "${PANEL_DIR}/prepared/query_manifest.tsv" ]]; then
      echo "    - Skipping benchmarks for ${panel} (no query manifest produced)"
      mark_failure "real.${panel}.no_queries"
      continue
    fi
    if [[ $(wc -l < "${PANEL_DIR}/prepared/query_manifest.tsv") -le 1 ]]; then
      echo "    - Skipping benchmarks for ${panel} (query manifest is empty)"
      mark_failure "real.${panel}.empty_queries"
      continue
    fi

    for mode in "${REAL_MODES[@]}"; do
      [[ -z "${mode}" ]] && continue
      echo "    - Benchmarking mode: ${mode}..."
      mkdir -p "${PANEL_DIR}/benchmark_${mode}"

      # Auto-enable every comparator whose binaries are detected for this mode.
      # Prefer this over hard-listing flags so a missing tool is reported and
      # skipped instead of crashing the whole panel/mode. The benchmark itself
      # prepends the project conda env bin to PATH at startup, so the operator
      # does not need to `conda activate` first.
      comparator_flags=(--run-all-comparators)

      # Override for assembly mode on heavy panels: skip cactus/pggb/anchorwave
      # (the multi-Gbp killers) and instead enable only the lighter comparators
      # explicitly. Without this carve-out a single AMF assembly run consumed
      # the entire SLURM time budget in the previous matrix run.
      if [[ "${mode}" == "assembly" ]] && [[ ",${HEAVY_COMPARATOR_SKIP_PANELS}," == *",${panel},"* ]]; then
        comparator_flags=(--run-syri --run-minigraph --run-svim-asm)
        echo "      [skip] heavy comparators (cactus/pggb/anchorwave) for ${panel} (HEAVY_COMPARATOR_SKIP_PANELS)"
      fi

      # Per-clade RAM cap for the MycoSV index build. Default 8 fits
      # comfortably in a 12 GiB cgroup; raise on larger nodes.
      REAL_MAX_CLADE_GENOMES="${REAL_MAX_CLADE_GENOMES:-8}"
      benchmark_threads="${THREADS}"
      mycosv_extra_args=()
      if [[ "${mode}" == "assembly" ]] && [[ ",${HEAVY_COMPARATOR_SKIP_PANELS}," == *",${panel},"* ]]; then
        benchmark_threads="${REAL_HEAVY_ASSEMBLY_THREADS}"
        # AMF panels need a much looser fragmentation cap because Gigaspora /
        # Rhizophagus assemblies legitimately have 50–100k+ contigs. Apply the
        # AMF-specific cap only for amf_large; the other heavy panels (te_rich,
        # two_speed) keep the tighter default to stop noisy MAG/contaminant
        # drafts sneaking through.
        if [[ "${panel}" == "amf_large" ]]; then
          panel_max_contigs="${REAL_HEAVY_MAX_ASSEMBLY_QUERY_CONTIGS_AMF}"
          panel_max_bp="${REAL_HEAVY_MAX_ASSEMBLY_QUERY_BP_AMF}"
        else
          panel_max_contigs="${REAL_HEAVY_MAX_ASSEMBLY_QUERY_CONTIGS}"
          panel_max_bp="${REAL_HEAVY_MAX_ASSEMBLY_QUERY_BP}"
        fi
        mycosv_extra_args=(
          --mycosv-arg=--no-gfa
          --mycosv-arg=--sa-max-contig-mb
          --mycosv-arg="${REAL_HEAVY_SA_MAX_CONTIG_MB}"
          --mycosv-arg=--single-ref-cache-mb
          --mycosv-arg="${REAL_HEAVY_SINGLE_REF_CACHE_MB}"
          --mycosv-arg=--max-ref-memory-mb
          --mycosv-arg="${REAL_HEAVY_MAX_REF_MEMORY_MB}"
          --mycosv-arg=--no-flat-ref-fallback
          --mycosv-arg=--tol-min-chain-anchors
          --mycosv-arg=2
          --max-assembly-query-contigs
          "${panel_max_contigs}"
          --max-assembly-query-bp
          "${panel_max_bp}"
        )
        echo "      [mem] heavy assembly settings: threads=${benchmark_threads}, MycoSV SA contig cap=${REAL_HEAVY_SA_MAX_CONTIG_MB} MB, SA cache=${REAL_HEAVY_SINGLE_REF_CACHE_MB} MB"
        echo "      [filter] heavy assembly query caps: contigs<=${panel_max_contigs}, bp<=${panel_max_bp}"
      fi

      # Wrap the benchmark invocation in `timeout` so a single slow panel
      # (e.g. cactus on AMF-scale genomes) cannot consume the SLURM time
      # budget and starve every subsequent panel. Long-form duration suffixes
      # (e.g. "4h", "30m") are passed straight through to GNU `timeout`.
      if [[ "${mode}" == "assembly" ]]; then
        bench_timeout="${BENCHMARK_TIMEOUT_ASSEMBLY}"
        read_validation_min_support="${REAL_READ_VALIDATION_MIN_SUPPORT_ASSEMBLY:-1}"
      else
        bench_timeout="${BENCHMARK_TIMEOUT_READS}"
        read_validation_min_support="${REAL_READ_VALIDATION_MIN_SUPPORT_READS:-3}"
      fi
      # Per-panel/per-mode timeout overrides. te_rich_pathogen (Puccinia ~80 Mbp,
      # ~85% TE content) consistently overruns the global 3h assembly cap even
      # after we drop cactus/pggb/anchorwave — minigraph + svim_asm + MycoSV's
      # own SA walker against multiple references is the bottleneck. Allow each
      # panel/mode to claim its own budget without bumping the global default
      # for fast panels like compact_yeast.
      panel_mode_upper="$(echo "${panel}_${mode}" | tr '[:lower:]-' '[:upper:]_')"
      panel_timeout_var="BENCHMARK_TIMEOUT_${panel_mode_upper}"
      if [[ -n "${!panel_timeout_var:-}" ]]; then
        bench_timeout="${!panel_timeout_var}"
      fi
      bench_cmd=(python3 -u run_real_fungal_benchmark.py benchmark
        --prepared-dir "${PANEL_DIR}/prepared"
        --mode "${mode}"
        --out-dir "${PANEL_DIR}/benchmark_${mode}"
        --threads "${benchmark_threads}"
        --max-clade-genomes "${REAL_MAX_CLADE_GENOMES}"
        --read-validation-min-support "${read_validation_min_support}"
        --max-comparator-short-reads "${MAX_COMPARATOR_SHORT_READS}"
        --max-comparator-long-reads "${MAX_COMPARATOR_LONG_READS}"
        "${mycosv_extra_args[@]}"
        "${comparator_flags[@]}")
      if [[ -n "${bench_timeout}" && "${bench_timeout}" != "0" ]] && command -v timeout >/dev/null; then
        bench_cmd=(timeout --signal=TERM --kill-after=60 "${bench_timeout}" "${bench_cmd[@]}")
      fi
      if "${bench_cmd[@]}" 2>&1 | tee "${PANEL_DIR}/benchmark_${mode}.log"; then
        mark_success "real.${panel}.${mode}"
      else
        rc=${PIPESTATUS[0]}
        # GNU `timeout` exits 124 on TERM and 137 on KILL — surface that as
        # a structured failure tag so the summary makes the bottleneck
        # obvious instead of just printing "stage failed".
        if [[ "${rc}" == "124" || "${rc}" == "137" ]]; then
          mark_failure "real.${panel}.${mode}.timeout(${bench_timeout})"
        else
          mark_failure "real.${panel}.${mode}"
        fi
      fi
    done
  done

  echo ""
  echo -e "${GREEN}✓ Real fungal data benchmarks stage complete${NC}"
  echo ""
fi


# ============================================================================
# 4. VISUALIZATION REPORT
# ============================================================================
#
# Builds an integrated HTML report across simulated benchmarks, real-data SV
# results, and biological findings when available. Missing inputs are handled
# gracefully so this stage never blocks the rest of the experiment matrix.
# ============================================================================

REPORT_DIR="${WORK_DIR}/experiments/reports/${TIMESTAMP}"
mkdir -p "${REPORT_DIR}"

if [[ "$EXPERIMENT_TYPE" == "all" || "$EXPERIMENT_TYPE" == "simulated" || "$EXPERIMENT_TYPE" == "real" ]]; then
  echo -e "${YELLOW}[4/4] Generating visualization report...${NC}"
  echo "      Output: ${REPORT_DIR}"

  SIM_RESULTS=""
  if [[ -f "${SIM_DIR}/benchmarks/million_mode_summary.tsv" ]]; then
    SIM_RESULTS="${SIM_DIR}/benchmarks/million_mode_summary.tsv"
  elif [[ -f "${SIM_DIR}/benchmarks/pr_metrics_simulated_summary.tsv" ]]; then
    SIM_RESULTS="${SIM_DIR}/benchmarks/pr_metrics_simulated_summary.tsv"
  fi

  REAL_RESULTS="${REPORT_DIR}/real_merged.tsv"
  BIO_RESULTS="${REPORT_DIR}/biology_merged.tsv"
  NOVEL_RESULTS="${REPORT_DIR}/novel_merged.tsv"
  : > "${REAL_RESULTS}"
  : > "${BIO_RESULTS}"
  : > "${NOVEL_RESULTS}"

  # Schema-aware merge: input TSVs may have different (but overlapping) headers.
  # We union all columns, then emit one merged TSV with empty fills for columns
  # missing in any source file. Falls back to a no-op when no inputs are given.
  merge_tsv_group() {
    local out_file="$1"
    shift
    [[ $# -eq 0 ]] && return 0
    python3 - "$out_file" "$@" <<'PY'
import csv, sys
out_path = sys.argv[1]
inputs = sys.argv[2:]
rows = []
columns = []
seen = set()
for path in inputs:
    try:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if not reader.fieldnames:
                continue
            for col in reader.fieldnames:
                if col not in seen:
                    seen.add(col)
                    columns.append(col)
            for row in reader:
                rows.append(row)
    except OSError:
        continue
if not columns:
    open(out_path, "w").close()
    sys.exit(0)
with open(out_path, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                            extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow({c: row.get(c, "") for c in columns})
PY
    return 0
  }

  # Scan REAL_DIR (per-panel benchmarks) AND MILLION_REAL_DIR (the
  # MycoSV-only million-scale benchmark from step 2b). Without the latter,
  # the held-out queries' SV calls / TE classifications / biology candidates
  # never reach the report and the million-real flow looks empty.
  mapfile -t REAL_TSVS < <(find "${REAL_DIR}" "${MILLION_REAL_DIR}" -type f \( \
      -name "*summary*.tsv" -o \
      -name "*pr_metrics*.tsv" -o \
      -name "*normalized_calls*.tsv" -o \
      -name "*score*.tsv" \
    \) 2>/dev/null | sort)

  mapfile -t BIO_TSVS < <(find "${REAL_DIR}" "${MILLION_REAL_DIR}" -type f \( \
      -name "*biology*.tsv" -o \
      -name "*candidate*.tsv" -o \
      -name "*annotation*.tsv" -o \
      -name "*pathway*.tsv" \
    \) 2>/dev/null | sort)

  mapfile -t NOVEL_TSVS < <(find "${REAL_DIR}" "${MILLION_REAL_DIR}" -type f \
      -name "novel_mycosv_calls.tsv" 2>/dev/null | sort)

  if [[ ${#REAL_TSVS[@]} -gt 0 ]]; then
    merge_tsv_group "${REAL_RESULTS}" "${REAL_TSVS[@]}"
  fi
  if [[ ${#BIO_TSVS[@]} -gt 0 ]]; then
    merge_tsv_group "${BIO_RESULTS}" "${BIO_TSVS[@]}"
  fi
  if [[ ${#NOVEL_TSVS[@]} -gt 0 ]]; then
    merge_tsv_group "${NOVEL_RESULTS}" "${NOVEL_TSVS[@]}"
  fi

  report_cmd=(python3 sv_visualization_report.py --outdir "${REPORT_DIR}" --title "MycoSV comprehensive report (${TIMESTAMP})")
  [[ -n "${SIM_RESULTS}" ]] && report_cmd+=(--simulated "${SIM_RESULTS}")
  [[ -s "${REAL_RESULTS}" ]] && report_cmd+=(--real "${REAL_RESULTS}")
  [[ -s "${BIO_RESULTS}" ]] && report_cmd+=(--biology "${BIO_RESULTS}")
  [[ -s "${NOVEL_RESULTS}" ]] && report_cmd+=(--novel "${NOVEL_RESULTS}")
  if [[ -n "${REPORT_TIMEOUT}" && "${REPORT_TIMEOUT}" != "0" ]] && command -v timeout >/dev/null; then
    report_cmd=(timeout --signal=TERM --kill-after=60 "${REPORT_TIMEOUT}" "${report_cmd[@]}")
  fi

  if [[ -f "${WORK_DIR}/sv_visualization_report.py" ]]; then
    if "${report_cmd[@]}" 2>&1 | tee "${REPORT_DIR}/report.log"; then
      if [[ -f "${REPORT_DIR}/sv_visualization_report.html" ]]; then
        mark_success "report.visualization"
        echo -e "${GREEN}✓ Visualization report generated${NC}"
      else
        mark_failure "report.visualization_missing_output"
      fi
    else
      rc=${PIPESTATUS[0]}
      if [[ "${rc}" == "124" || "${rc}" == "137" ]]; then
        mark_failure "report.visualization.timeout(${REPORT_TIMEOUT})"
      else
        mark_failure "report.visualization"
      fi
    fi
  else
    echo "      Report script not found: ${WORK_DIR}/sv_visualization_report.py"
    mark_failure "report.script_missing"
  fi
  echo ""
fi

# ============================================================================
# SUMMARY REPORT
# ============================================================================

echo -e "${BLUE}========================================"
echo "Experiment Summary"
echo "========================================${NC}"
echo ""
echo "All experiment outputs saved to:"
echo "  Simulated:     ${SIM_DIR}"
echo "  Real data:     ${REAL_DIR}"
echo "  Million-real:  ${MILLION_REAL_DIR}"
echo "  Report:        ${REPORT_DIR}"
echo ""
echo "Stage outcomes:"
echo "  Succeeded: ${#SUCCESS_STAGES[@]}"
echo "  Failed:    ${#FAILED_STAGES[@]}"
if [[ ${#FAILED_STAGES[@]} -gt 0 ]]; then
  echo "  Failed stages:"
  for s in "${FAILED_STAGES[@]}"; do echo "    - $s"; done
fi
echo ""
echo "Intermediate files preserved:"
find "${SIM_DIR}" "${REAL_DIR}" "${MILLION_REAL_DIR}" "${REPORT_DIR}" -type f \
  \( -name "*.vcf" -o -name "*.vcf.gz" -o -name "*.tsv" -o -name "*.fasta" -o -name "*.fastq" -o -name "*.html" -o -name "*.png" \) \
  2>/dev/null | wc -l | xargs echo "  Total:"
echo ""
echo "Disk usage:"
du -sh "${SIM_DIR}" "${REAL_DIR}" "${MILLION_REAL_DIR}" "${REPORT_DIR}" 2>/dev/null | awk '{print "  " $0}'
echo ""
echo "Log files:"
find "${SIM_DIR}" "${REAL_DIR}" "${MILLION_REAL_DIR}" "${REPORT_DIR}" -name "*.log" 2>/dev/null | wc -l | xargs echo "  Total:"
echo ""
echo -e "${GREEN}✓ All experiments complete! End time: $(date)${NC}"
echo ""
echo "Next steps:"
echo "  1. Review logs for errors: grep -r 'ERROR' ${SIM_DIR} ${REAL_DIR} ${MILLION_REAL_DIR} ${REPORT_DIR}"
echo "  2. Open report: ${REPORT_DIR}/sv_visualization_report.html"
echo ""

if [[ ${#SUCCESS_STAGES[@]} -eq 0 ]]; then
  exit 1
fi
exit 0
