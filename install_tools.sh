#!/usr/bin/env bash
# install_tools.sh - install MycoSV comparator tools without Bioconda.
#
# Installs the comparator pre-flight dependency closure from GitHub source trees,
# GitHub-hosted release artifacts, or official GitHub/OCI containers for very
# large pangenome stacks.
#
# Comparator closure covered:
#   syri        : minimap2 + syri
#   minigraph  : minigraph + gfatools
#   pggb       : official GHCR pggb image wrapper, source cloned for provenance
#   cactus     : cactus-pangenome official release/container wrapper
#   svim_asm   : svim-asm + minimap2 + samtools
#   anchorwave : anchorwave + minimap2 + samtools + paftools.js
#   delly/manta helpers: bcftools for BCF/VCF normalization
#
# Usage:
#   bash install_tools.sh              # install comparator closure
#   bash install_tools.sh --check      # report tool availability only
#   bash install_tools.sh python-deps  # install Python reporting/test deps
#   bash install_tools.sh all          # install broad SV benchmark stack
#   bash install_tools.sh syri pggb    # install selected tools
#
# Environment overrides:
#   ENV_PATH=/path/to/conda/env
#   CONDA_INIT=/path/to/conda.sh
#   WORK_DIR=/path/to/project/workdir
#   SRC_DIR=/path/to/source/cache
#   THREADS=8
#   FORCE=1
#   USE_APPTAINER=1
#   CACTUS_MODE=container|release     # default: container

set -u
set -o pipefail

readonly DEFAULT_ENV_PATH="/mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/tools/envs/envs/fungi_graph_sv"
readonly DEFAULT_CONDA_INIT="/opt/apps/apps/binapps/conda/miniforge3/25.9.1-0/etc/profile.d/conda.sh"

ENV_PATH="${ENV_PATH:-$DEFAULT_ENV_PATH}"
CONDA_INIT="${CONDA_INIT:-$DEFAULT_CONDA_INIT}"
WORK_DIR="${WORK_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
SRC_DIR="${SRC_DIR:-$WORK_DIR/tools_src}"
THREADS="${THREADS:-4}"
FORCE="${FORCE:-0}"
USE_APPTAINER="${USE_APPTAINER:-1}"
CACTUS_MODE="${CACTUS_MODE:-container}"

GREEN='\033[0;32m'; RED='\033[0;31m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'
ok()   { printf '%b  OK    %s%b\n' "$GREEN" "$*" "$NC"; }
fail() { printf '%b  FAIL  %s%b\n' "$RED" "$*" "$NC"; }
info() { printf '%b  ->    %s%b\n' "$BLUE" "$*" "$NC"; }
hdr()  { printf '\n%b== %s ==%b\n' "$YELLOW" "$*" "$NC"; }

have() { command -v "$1" >/dev/null 2>&1; }
need() { ! have "$1" || [[ "$FORCE" == "1" ]]; }

require_cmd() {
    local missing=0 cmd
    for cmd in "$@"; do
        if ! have "$cmd"; then
            fail "required command missing: $cmd"
            missing=1
        fi
    done
    [[ "$missing" == 0 ]]
}

activate_env() {
    if [[ ! -f "$CONDA_INIT" ]]; then
        echo "[error] conda init script not found: $CONDA_INIT" >&2
        return 1
    fi
    # shellcheck disable=SC1090
    source "$CONDA_INIT"
    conda activate "$ENV_PATH"
    ENV_BIN="$CONDA_PREFIX/bin"
    mkdir -p "$SRC_DIR" "$ENV_BIN"
}

pip_install() { python -m pip install "$@"; }
clone_fresh() { local repo="$1" dest="$2"; rm -rf "$dest"; git clone --depth 1 "$repo" "$dest"; }

install_minimap2() {
    if ! need minimap2; then ok "minimap2 already on PATH"; return 0; fi
    require_cmd git make install || return 1
    clone_fresh https://github.com/lh3/minimap2.git "$SRC_DIR/minimap2" || return 1
    (cd "$SRC_DIR/minimap2" && make -j"$THREADS") || return 1
    install -m 755 "$SRC_DIR/minimap2/minimap2" "$ENV_BIN/minimap2"
    have minimap2
}

install_k8() {
    if ! need k8; then ok "k8 already on PATH"; return 0; fi
    require_cmd git make install || return 1
    clone_fresh https://github.com/lh3/k8.git "$SRC_DIR/k8" || return 1
    (cd "$SRC_DIR/k8" && make -j"$THREADS") || return 1
    install -m 755 "$SRC_DIR/k8/k8" "$ENV_BIN/k8"
    have k8
}

install_paftools_js() {
    if ! need paftools.js; then ok "paftools.js already on PATH"; return 0; fi
    require_cmd git install || return 1
    install_k8 || return 1
    if [[ ! -f "$SRC_DIR/minimap2/misc/paftools.js" || "$FORCE" == "1" ]]; then
        clone_fresh https://github.com/lh3/minimap2.git "$SRC_DIR/minimap2" || return 1
    fi
    install -m 755 "$SRC_DIR/minimap2/misc/paftools.js" "$ENV_BIN/paftools.js"
    have paftools.js
}

install_minigraph() {
    if ! need minigraph; then ok "minigraph already on PATH"; return 0; fi
    require_cmd git make install || return 1
    clone_fresh https://github.com/lh3/minigraph.git "$SRC_DIR/minigraph" || return 1
    (cd "$SRC_DIR/minigraph" && make -j"$THREADS") || return 1
    install -m 755 "$SRC_DIR/minigraph/minigraph" "$ENV_BIN/minigraph"
    have minigraph
}

install_gfatools() {
    if ! need gfatools; then ok "gfatools already on PATH"; return 0; fi
    require_cmd git make install || return 1
    clone_fresh https://github.com/lh3/gfatools.git "$SRC_DIR/gfatools" || return 1
    (cd "$SRC_DIR/gfatools" && make -j"$THREADS") || return 1
    install -m 755 "$SRC_DIR/gfatools/gfatools" "$ENV_BIN/gfatools"
    have gfatools
}

install_htslib() {
    if have bgzip && have tabix && [[ "$FORCE" != "1" ]]; then ok "htslib tools already on PATH"; return 0; fi
    require_cmd git make install autoreconf || return 1
    clone_fresh https://github.com/samtools/htslib.git "$SRC_DIR/htslib" || return 1
    (cd "$SRC_DIR/htslib" && autoreconf -i && ./configure --prefix="$CONDA_PREFIX" && make -j"$THREADS" && make install) || return 1
    have bgzip && have tabix
}

install_samtools() {
    if ! need samtools; then ok "samtools already on PATH"; return 0; fi
    require_cmd git make install autoreconf || return 1
    install_htslib || return 1
    clone_fresh https://github.com/samtools/samtools.git "$SRC_DIR/samtools" || return 1
    (cd "$SRC_DIR/samtools" && autoreconf -i && ./configure --prefix="$CONDA_PREFIX" --with-htslib="$CONDA_PREFIX" && make -j"$THREADS" && make install) || return 1
    have samtools
}

install_bcftools() {
    if ! need bcftools; then ok "bcftools already on PATH"; return 0; fi
    require_cmd git make install autoreconf || return 1
    install_htslib || return 1
    clone_fresh https://github.com/samtools/bcftools.git "$SRC_DIR/bcftools" || return 1
    (cd "$SRC_DIR/bcftools" && autoreconf -i && ./configure --prefix="$CONDA_PREFIX" --with-htslib="$CONDA_PREFIX" && make -j"$THREADS" && make install) || return 1
    have bcftools
}

install_syri() {
    if ! need syri; then ok "syri already on PATH"; return 0; fi
    require_cmd git || return 1
    install_minimap2 || return 1
    install_python-deps || return 1
    clone_fresh https://github.com/schneebergerlab/syri.git "$SRC_DIR/syri-src" || return 1
    (cd "$SRC_DIR/syri-src" && pip_install --no-build-isolation .) 2>&1 | tail -20 || return 1
    have syri
}

install_python-deps() {
    local missing=() pkg
    for pkg in Cython numpy psutil pandas pysam scipy pytest matplotlib; do
        python - "$pkg" <<'PY' >/dev/null 2>&1 || missing+=("$pkg")
import importlib.util
import sys
pkg = sys.argv[1]
module = {"Cython": "Cython"}.get(pkg, pkg)
raise SystemExit(0 if importlib.util.find_spec(module) else 1)
PY
    done
    if [[ ${#missing[@]} -eq 0 && "$FORCE" != "1" ]]; then
        ok "Python deps already installed"
        return 0
    fi
    info "Installing Python deps: ${missing[*]:-Cython numpy psutil pandas pysam scipy pytest matplotlib}"
    pip_install -U pip setuptools wheel 2>&1 | tail -5 || true
    if [[ "$FORCE" == "1" ]]; then
        pip_install -U Cython numpy psutil pandas pysam scipy pytest matplotlib 2>&1 | tail -20 || return 1
    else
        pip_install "${missing[@]}" 2>&1 | tail -20 || return 1
    fi
    python - <<'PY'
import importlib
for pkg in ["Cython", "numpy", "psutil", "pandas", "pysam", "scipy", "pytest", "matplotlib"]:
    importlib.import_module(pkg)
print("Python deps OK")
PY
}

install_svim() {
    if ! need svim; then ok "svim already on PATH"; return 0; fi
    require_cmd git || return 1
    pip_install "git+https://github.com/eldariont/svim.git" 2>&1 | tail -20 || return 1
    have svim
}

install_svim-asm() {
    if ! need svim-asm; then ok "svim-asm already on PATH"; return 0; fi
    require_cmd git || return 1
    install_minimap2 || return 1
    install_samtools || return 1
    pip_install "git+https://github.com/eldariont/svim-asm.git" 2>&1 | tail -20 || return 1
    have svim-asm
}

install_sniffles() {
    if ! need sniffles; then ok "sniffles already on PATH"; return 0; fi
    require_cmd git || return 1
    pip_install "git+https://github.com/fritzsedlazeck/Sniffles.git" 2>&1 | tail -20 || return 1
    have sniffles
}

install_cuteSV() {
    if ! need cuteSV; then ok "cuteSV already on PATH"; return 0; fi
    require_cmd git || return 1
    pip_install --no-build-isolation "git+https://github.com/brentp/cigar.git" 2>&1 | tail -5 || true
    clone_fresh https://github.com/tjiangHIT/cuteSV.git "$SRC_DIR/cuteSV-src" || return 1
    (cd "$SRC_DIR/cuteSV-src" && pip_install --no-build-isolation .) 2>&1 | tail -20 || return 1
    have cuteSV
}

install_anchorwave() {
    if ! need anchorwave; then ok "anchorwave already on PATH"; return 0; fi
    require_cmd git cmake make install || return 1
    install_minimap2 || return 1
    install_samtools || return 1
    install_paftools_js || return 1
    clone_fresh https://github.com/baoxingsong/AnchorWave.git "$SRC_DIR/AnchorWave" || return 1
    (cd "$SRC_DIR/AnchorWave" && cmake -DCMAKE_BUILD_TYPE=Release . && make -j"$THREADS") || return 1
    install -m 755 "$SRC_DIR/AnchorWave/anchorwave" "$ENV_BIN/anchorwave"
    have anchorwave
}

install_delly() {
    if ! need delly; then ok "delly already on PATH"; return 0; fi
    install_bcftools || return 1
    require_cmd curl chmod || return 1
    local version="v1.7.3"
    curl -fsSL -o "$ENV_BIN/delly" "https://github.com/dellytools/delly/releases/download/${version}/delly-${version}-linux-amd64" || return 1
    chmod +x "$ENV_BIN/delly"
    have delly
}

install_manta() {
    if ! need configManta.py; then ok "manta already on PATH"; return 0; fi
    require_cmd curl tar ln || return 1
    local version="1.6.0" tarball="$SRC_DIR/manta-${version}.tar.bz2" root="$SRC_DIR/manta-${version}.centos6_x86_64"
    rm -rf "$root" "$tarball"
    curl -fsSL -o "$tarball" "https://github.com/Illumina/manta/releases/download/v${version}/manta-${version}.centos6_x86_64.tar.bz2" || return 1
    tar -xjf "$tarball" -C "$SRC_DIR" || return 1
    [[ -x "$root/bin/configManta.py" ]] || return 1
    ln -sf "$root/bin/configManta.py" "$ENV_BIN/configManta.py"
    ln -sf "$root/bin/configManta.py" "$ENV_BIN/manta"
    have configManta.py
}

write_apptainer_wrapper() {
    local sif="$1"; shift
    local tool
    for tool in "$@"; do
        if [[ "$tool" == "pggb" ]]; then
            cat > "$ENV_BIN/$tool" <<EOF2
#!/usr/bin/env bash
# pggb's container can inherit host shell aliases/functions that make its
# internal \`which time\` probe call debianutils \`which\` with unsupported
# flags ("Illegal option --"). Start from a clean environment and a non-login
# shell so the resolver uses the plain executable lookup path.
exec apptainer exec --cleanenv --bind /mnt --bind /tmp "$sif" \\
    env -i HOME=/tmp PATH=/usr/local/bin:/usr/bin:/bin LC_ALL=C \\
    bash --noprofile --norc /usr/local/bin/pggb "\$@"
EOF2
        else
            cat > "$ENV_BIN/$tool" <<EOF2
#!/usr/bin/env bash
exec apptainer exec --bind /mnt --bind /tmp "$sif" "$tool" "\$@"
EOF2
        fi
        chmod +x "$ENV_BIN/$tool"
    done
}

pull_apptainer_image() {
    local sif="$1" image="$2"
    require_cmd apptainer || return 1
    if [[ ! -s "$sif" || "$FORCE" == "1" ]]; then
        apptainer pull --force "$sif" "$image" || return 1
    fi
}

install_pggb() {
    if ! need pggb; then ok "pggb already on PATH"; return 0; fi
    [[ "$USE_APPTAINER" == "1" ]] || { fail "pggb dependency stack requires Apptainer; set USE_APPTAINER=1"; return 1; }
    require_cmd git || return 1
    clone_fresh https://github.com/pangenome/pggb.git "$SRC_DIR/pggb-src" || info "pggb source clone failed (provenance only); continuing with apptainer image"
    local sif="$SRC_DIR/pggb.sif"
    pull_apptainer_image "$sif" docker://ghcr.io/pangenome/pggb:latest || return 1
    write_apptainer_wrapper "$sif" pggb
    have pggb
}

install_cactus_release() {
    require_cmd curl tar ln || return 1
    local version="3.1.4" tarball="$SRC_DIR/cactus-bin-v${version}.tar.gz" root="$SRC_DIR/cactus-bin-v${version}"
    rm -rf "$root" "$tarball"
    curl -fsSL -o "$tarball" "https://github.com/ComparativeGenomicsToolkit/cactus/releases/download/v${version}/cactus-bin-v${version}.tar.gz" || return 1
    tar -xzf "$tarball" -C "$SRC_DIR" || return 1
    if [[ -f "$root/setup.py" ]]; then
        pip_install --no-build-isolation -U setuptools wheel 2>&1 | tail -5 || true
        pip_install --no-build-isolation "$root" 2>&1 | tail -20 || true
        [[ -f "$root/toil-requirement.txt" ]] && pip_install --no-build-isolation -r "$root/toil-requirement.txt" 2>&1 | tail -10 || true
    fi
    local f b
    for f in "$root"/bin/*; do
        [[ -e "$f" ]] || continue
        b="$(basename "$f")"
        ln -sf "$f" "$ENV_BIN/$b"
    done
    have cactus-pangenome
}

install_cactus-pangenome() {
    if ! need cactus-pangenome; then ok "cactus-pangenome already on PATH"; return 0; fi
    if [[ "$CACTUS_MODE" == "release" ]]; then
        install_cactus_release && return 0
        return 1
    fi
    [[ "$USE_APPTAINER" == "1" ]] || { fail "cactus container install requires USE_APPTAINER=1; or set CACTUS_MODE=release"; return 1; }
    local sif="$SRC_DIR/cactus.sif"
    pull_apptainer_image "$sif" docker://quay.io/comparative-genomics-toolkit/cactus:v3.1.4 || return 1
    write_apptainer_wrapper "$sif" cactus cactus-pangenome cactus-prepare cactus-graphmap cactus-update-prepare
    have cactus-pangenome
}

readonly COMPARATOR_TOOLS=(python-deps minimap2 paftools_js syri minigraph gfatools samtools bcftools svim-asm anchorwave pggb cactus-pangenome)
readonly ALL_TOOLS=(python-deps minimap2 k8 paftools_js minigraph gfatools htslib samtools bcftools syri svim-asm anchorwave pggb cactus-pangenome svim sniffles cuteSV delly manta)
readonly CHECK_TOOLS=(minimap2 paftools.js syri minigraph gfatools pggb cactus-pangenome svim-asm samtools bcftools anchorwave svim sniffles cuteSV delly configManta.py)

normalize_tool() {
    case "$1" in
        cactus) echo cactus-pangenome ;;
        python|python_deps|python-deps|deps|pythondeps) echo python-deps ;;
        paftools|paftools.js|paftools_js) echo paftools_js ;;
        sniffles2) echo sniffles ;;
        cutesv) echo cuteSV ;;
        svim_asm) echo svim-asm ;;
        all) echo __ALL__ ;;
        comparators|benchmarking|preflight) echo __COMPARATORS__ ;;
        *) echo "$1" ;;
    esac
}

check_tools() {
    hdr "Comparator/tool availability"
    local missing=0 tool path
    for tool in "${CHECK_TOOLS[@]}"; do
        if have "$tool"; then
            path="$(command -v "$tool")"
            ok "$(printf '%-18s %s' "$tool" "$path")"
        else
            fail "$(printf '%-18s missing' "$tool")"
            missing=$((missing + 1))
        fi
    done
    echo
    hdr "Python package availability"
    local pkg module
    for pkg in Cython numpy psutil pandas pysam scipy pytest matplotlib; do
        module="$pkg"
        python - "$module" <<'PY' >/dev/null 2>&1
import importlib.util
import sys
raise SystemExit(0 if importlib.util.find_spec(sys.argv[1]) else 1)
PY
        if [[ "$?" == "0" ]]; then
            ok "$(printf '%-18s present' "$pkg")"
        else
            fail "$(printf '%-18s missing' "$pkg")"
            missing=$((missing + 1))
        fi
    done
    echo
    [[ "$missing" == 0 ]] && echo "All checked tools are present." || echo "Missing tools: $missing"
}

run_one() {
    local tool fn
    tool="$(normalize_tool "$1")"
    fn="install_${tool}"
    if ! declare -F "$fn" >/dev/null; then
        fail "unknown tool: $1"
        return 1
    fi
    hdr "$tool"
    if "$fn"; then ok "$tool installed"; return 0; fi
    fail "$tool install failed"
    return 1
}

print_usage() { sed -n '1,31p' "$0" | sed 's/^# \{0,1\}//'; }

main() {
    if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then print_usage; exit 0; fi
    activate_env || exit 1
    if [[ "${1:-}" == "--check" || "${1:-}" == "check" ]]; then check_tools; exit 0; fi

    local raw targets=() norm failed=0 t
    if [[ "$#" -eq 0 ]]; then
        targets=("${COMPARATOR_TOOLS[@]}")
    else
        for raw in "$@"; do
            norm="$(normalize_tool "$raw")"
            if [[ "$norm" == "__ALL__" ]]; then targets=("${ALL_TOOLS[@]}"); break; fi
            if [[ "$norm" == "__COMPARATORS__" ]]; then targets=("${COMPARATOR_TOOLS[@]}"); break; fi
            targets+=("$norm")
        done
    fi

    hdr "Install plan"
    echo "Env         : $ENV_PATH"
    echo "Bin         : $ENV_BIN"
    echo "Source/cache: $SRC_DIR"
    echo "Threads     : $THREADS"
    echo "Force       : $FORCE"
    echo "Cactus mode : $CACTUS_MODE"
    echo "Tools       : ${targets[*]}"

    for t in "${targets[@]}"; do
        run_one "$t" || failed=$((failed + 1))
    done
    check_tools
    if [[ "$failed" -gt 0 ]]; then fail "$failed install stage(s) failed"; exit 1; fi
}

main "$@"
