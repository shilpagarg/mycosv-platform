#!/usr/bin/env bash
# install_tools.sh — Install all tools required by the MycoSV benchmark environment.
#
# Sets up a conda environment named "mycosv" containing:
#   Core build tools   : g++, cmake, make
#   SV callers         : SyRI+minimap2, minigraph, PGGB, Minigraph-Cactus,
#                        SVIM-asm, AnchorWave, Delly, Manta,
#                        SVIM, Sniffles2, cuteSV, samtools, bcftools
#   TE classifiers     : DeepTE, NeuralTE, TERL, Terrier, ClassifyTE, CREATE, TEClass2
#   Python packages    : biopython, pandas, scikit-learn, pytorch (CPU), tensorflow (CPU)
#   MycoSV binary      : compiled from main.cpp after environment activation
#
# Usage:
#   bash install_tools.sh                  # full install (creates/updates mycosv env)
#   bash install_tools.sh --check          # only check what is/isn't installed
#   bash install_tools.sh --sv-only        # SV tools only (skip TE classifiers)
#   bash install_tools.sh --te-only        # TE tools only (skip SV callers)
#   bash install_tools.sh --mycosv-only    # build MycoSV binary only
#
# Environment variables:
#   CONDA_ENV_NAME=mycosv      Override the conda environment name
#   MAMBA=1                    Use mamba instead of conda for speed
#   SKIP_CONDA=1               Skip conda steps (assume env already active)
#
# Requirements:
#   - conda or mamba on PATH
#   - Internet access for package downloads
#   - ~15 GB free disk space for full install

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="${CONDA_ENV_NAME:-mycosv}"
CONDA_CMD="${MAMBA:-0}" ; [[ "${MAMBA:-0}" == "1" ]] && CONDA_CMD="mamba" || CONDA_CMD="conda"
MODE="${1:-all}"
MODE="${MODE#--}"  # strip leading --

# ---------------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------------
GREEN='\033[0;32m'; RED='\033[0;31m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'
ok()   { echo -e "${GREEN}  ✓ $*${NC}"; }
warn() { echo -e "${YELLOW}  ⚠ $*${NC}"; }
fail() { echo -e "${RED}  ✗ $*${NC}"; }
info() { echo -e "${BLUE}  → $*${NC}"; }

# ---------------------------------------------------------------------------
# Tool check helpers
# ---------------------------------------------------------------------------
MISSING_TOOLS=()
AVAILABLE_TOOLS=()

check_tool() {
    local name="$1"
    if command -v "$name" &>/dev/null; then
        local ver
        ver=$(command -v "$name")
        AVAILABLE_TOOLS+=("$name")
        ok "$name  →  $ver"
        return 0
    else
        MISSING_TOOLS+=("$name")
        fail "$name  (not found)"
        return 1
    fi
}

check_python_module() {
    local mod="$1"
    if python3 -c "import $mod" 2>/dev/null; then
        ok "python:$mod"
    else
        MISSING_TOOLS+=("python:$mod")
        fail "python:$mod  (not importable)"
    fi
}

# ---------------------------------------------------------------------------
# --check mode: just report what is/isn't available
# ---------------------------------------------------------------------------
if [[ "$MODE" == "check" ]]; then
    echo -e "${BLUE}========================================"
    echo "MycoSV Tool Availability Check"
    echo "========================================${NC}"
    echo ""

    echo "=== MycoSV binary ==="
    check_tool "fungi_graphsv_tol" || true
    echo ""

    echo "=== Compiler ==="
    check_tool "g++" || true
    echo ""

    echo "=== Assembly-mode SV comparators ==="
    check_tool "minimap2"        || true
    check_tool "syri"            || true
    check_tool "minigraph"       || true
    check_tool "gfatools"        || true
    check_tool "pggb"            || true
    check_tool "cactus-pangenome"|| true
    check_tool "svim-asm"        || true
    check_tool "anchorwave"      || true
    echo ""

    echo "=== Short-reads SV comparators ==="
    check_tool "delly"           || true
    check_tool "configManta.py"  || true
    check_tool "samtools"        || true
    check_tool "bcftools"        || true
    echo ""

    echo "=== Long-reads SV comparators ==="
    check_tool "svim"            || true
    check_tool "sniffles"        || true
    check_tool "cuteSV"          || true
    echo ""

    echo "=== TE classifiers ==="
    check_tool "DeepTE.py"       || true
    check_tool "NeuralTE.py"     || true
    check_tool "TERL.py"         || true
    check_tool "terrier.py"      || true
    check_tool "ClassifyTE.py"   || true
    check_tool "CREATE.py"       || true
    check_tool "TEClass2"        || true
    echo ""

    echo "=== Python packages ==="
    check_python_module "Bio"    || true
    check_python_module "pandas" || true
    check_python_module "sklearn"|| true
    check_python_module "torch"  || true
    check_python_module "tensorflow" || true
    echo ""

    echo "========================================"
    echo "Summary:"
    echo "  Available : ${#AVAILABLE_TOOLS[@]}"
    echo "  Missing   : ${#MISSING_TOOLS[@]}"
    if [[ ${#MISSING_TOOLS[@]} -gt 0 ]]; then
        echo "  Missing tools:"
        for t in "${MISSING_TOOLS[@]}"; do echo "    - $t"; done
        echo ""
        echo "To install missing tools, run:"
        echo "  bash install_tools.sh"
    fi
    echo "========================================"
    exit 0
fi

# ---------------------------------------------------------------------------
# Require conda
# ---------------------------------------------------------------------------
if [[ "${SKIP_CONDA:-0}" != "1" ]]; then
    if ! command -v conda &>/dev/null && ! command -v mamba &>/dev/null; then
        echo -e "${RED}[error] conda or mamba not found on PATH.${NC}"
        echo "Install Miniconda from: https://docs.conda.io/en/latest/miniconda.html"
        echo "Then re-run this script."
        exit 1
    fi
fi

echo -e "${BLUE}========================================"
echo "MycoSV Environment Installer"
echo "========================================${NC}"
echo "Environment : ${ENV_NAME}"
echo "Installer   : ${CONDA_CMD}"
echo "Mode        : ${MODE}"
echo "Script dir  : ${SCRIPT_DIR}"
echo ""

# ---------------------------------------------------------------------------
# Create / activate environment
# ---------------------------------------------------------------------------
if [[ "${SKIP_CONDA:-0}" != "1" && "$MODE" != "mycosv-only" ]]; then
    if conda env list 2>/dev/null | grep -qE "^${ENV_NAME}\s"; then
        info "Conda environment '${ENV_NAME}' already exists — updating"
    else
        info "Creating conda environment '${ENV_NAME}' (Python 3.11)"
        ${CONDA_CMD} create -y -n "${ENV_NAME}" python=3.11
    fi

    # Activate for subsequent conda install calls
    # shellcheck disable=SC1091
    eval "$(conda shell.bash hook)"
    conda activate "${ENV_NAME}"
    info "Activated environment: ${ENV_NAME}"
fi

# ---------------------------------------------------------------------------
# Helper: install a conda package, report result
# ---------------------------------------------------------------------------
install_conda() {
    local pkg="$1"
    local channel="${2:--c conda-forge}"
    info "Installing conda package: $pkg"
    if ${CONDA_CMD} install -y ${channel} "$pkg" 2>&1 | tail -3; then
        ok "Installed: $pkg"
    else
        warn "Failed to install $pkg via conda — may need manual install"
    fi
}

install_pip() {
    local pkg="$1"
    info "Installing pip package: $pkg"
    if python3 -m pip install -q "$pkg"; then
        ok "Installed pip: $pkg"
    else
        warn "Failed to install pip: $pkg"
    fi
}

# ---------------------------------------------------------------------------
# 1. Build tools (always needed)
# ---------------------------------------------------------------------------
if [[ "$MODE" =~ ^(all|sv-only|mycosv-only)$ ]]; then
    echo ""
    echo -e "${YELLOW}[1] Build tools${NC}"
    install_conda "gxx_linux-64"  "-c conda-forge"
    install_conda "cmake"         "-c conda-forge"
    install_conda "make"          "-c conda-forge"
fi

# ---------------------------------------------------------------------------
# 2. Assembly-mode SV comparators
# ---------------------------------------------------------------------------
if [[ "$MODE" =~ ^(all|sv-only)$ ]]; then
    echo ""
    echo -e "${YELLOW}[2] Assembly-mode SV comparators${NC}"

    # minimap2 + samtools + bcftools (shared dependency)
    install_conda "minimap2"    "-c bioconda -c conda-forge"
    install_conda "samtools"    "-c bioconda -c conda-forge"
    install_conda "bcftools"    "-c bioconda -c conda-forge"
    install_conda "htslib"      "-c bioconda -c conda-forge"

    # SyRI — whole-genome synteny + SV detection (assembly vs reference)
    # Paper: Goel et al. (2019); doi:10.1186/s13059-019-1793-5
    install_conda "syri"        "-c bioconda -c conda-forge"

    # minigraph — pangenome graph builder + GFA-based SV calling
    # Paper: Li (2020); doi:10.1186/s13059-020-02168-z
    install_conda "minigraph"   "-c bioconda -c conda-forge"
    install_conda "gfatools"    "-c bioconda -c conda-forge"

    # PGGB — PanGenome Graph Builder (sequence-to-graph, odgi)
    # Paper: Garrison et al. (2023); doi:10.1038/s41592-023-02014-7
    install_conda "pggb"        "-c bioconda -c conda-forge"
    install_conda "odgi"        "-c bioconda -c conda-forge"
    install_conda "vg"          "-c bioconda -c conda-forge"

    # SVIM-asm — SV calling from haplotype-resolved assemblies via minimap2
    # Paper: Heller & Vingron (2021); doi:10.1093/bioinformatics/btab705
    install_conda "svim-asm"    "-c bioconda -c conda-forge"

    # AnchorWave — sensitive whole-genome alignment for plant/fungal genomes
    # Paper: Song et al. (2022); doi:10.1073/pnas.2106652119
    install_conda "anchorwave"  "-c bioconda -c conda-forge"
    # paftools.js (part of minimap2 package, but ensure it is on PATH)
    install_conda "k8"          "-c bioconda -c conda-forge"

    # Minigraph-Cactus (cactus-pangenome) — sequence-to-pangenome graph
    # Paper: Armstrong et al. (2020); doi:10.1038/s41592-020-0939-8
    # Note: large package; may take several minutes
    info "Installing cactus (Minigraph-Cactus) — this may take several minutes..."
    if ${CONDA_CMD} install -y -c bioconda -c conda-forge cactus 2>&1 | tail -3; then
        ok "Installed: cactus (cactus-pangenome)"
    else
        warn "cactus install failed — try: pip install progressiveCactus"
        warn "Or download a release binary: https://github.com/ComparativeGenomicsToolkit/cactus/releases"
    fi
fi

# ---------------------------------------------------------------------------
# 3. Short-reads SV comparators
# ---------------------------------------------------------------------------
if [[ "$MODE" =~ ^(all|sv-only)$ ]]; then
    echo ""
    echo -e "${YELLOW}[3] Short-reads SV comparators${NC}"

    # Delly — structural variant detection from paired-end reads
    # Paper: Rausch et al. (2012); doi:10.1093/bioinformatics/bts378
    install_conda "delly"       "-c bioconda -c conda-forge"

    # Manta — structural variant and indel caller from paired-end reads
    # Paper: Chen et al. (2016); doi:10.1093/bioinformatics/btv710
    install_conda "manta"       "-c bioconda -c conda-forge"
fi

# ---------------------------------------------------------------------------
# 4. Long-reads SV comparators
# ---------------------------------------------------------------------------
if [[ "$MODE" =~ ^(all|sv-only)$ ]]; then
    echo ""
    echo -e "${YELLOW}[4] Long-reads SV comparators${NC}"

    # SVIM — SV detection from long reads (ONT/PacBio)
    # Paper: Heller & Vingron (2019); doi:10.1093/bioinformatics/btz041
    install_conda "svim"        "-c bioconda -c conda-forge"

    # Sniffles2 — SV detection from long reads with population support
    # Paper: Sedlazeck et al. (2018); doi:10.1038/s41592-018-0001-7
    install_conda "sniffles"    "-c bioconda -c conda-forge"

    # cuteSV — sensitive SV detection from long reads
    # Paper: Jiang et al. (2020); doi:10.1186/s13059-020-02107-y
    install_conda "cutesv"      "-c bioconda -c conda-forge"
fi

# ---------------------------------------------------------------------------
# 5. TE classification tools
# ---------------------------------------------------------------------------
if [[ "$MODE" =~ ^(all|te-only)$ ]]; then
    echo ""
    echo -e "${YELLOW}[5] TE classification tools${NC}"

    # Python ML/DL dependencies shared by TE tools
    install_pip "biopython>=1.81"
    install_pip "pandas>=2.0"
    install_pip "scikit-learn>=1.3"
    install_pip "numpy>=1.24"

    # DeepTE — CNN-based TE classifier (fungi/animal/plant)
    # Paper: Yan et al. (2020); doi:10.1093/nar/gkaa323
    # Requires TensorFlow
    info "Installing TensorFlow (CPU) for DeepTE/CREATE..."
    install_pip "tensorflow-cpu>=2.12"
    if command -v pip3 &>/dev/null; then
        info "Installing DeepTE via pip (GitHub)..."
        python3 -m pip install -q \
            "git+https://github.com/LiLabAtVT/DeepTE.git" 2>/dev/null \
            || warn "DeepTE GitHub install failed — install manually from https://github.com/LiLabAtVT/DeepTE"
    fi

    # NeuralTE — transformer-based TE classifier (best fungal F1 in PanTEon)
    # Paper: Han et al. (2024); doi:10.1093/nar/gkae009
    # Requires PyTorch
    info "Installing PyTorch (CPU) for NeuralTE/TERL/Terrier..."
    install_pip "torch>=2.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu"
    python3 -m pip install -q \
        "git+https://github.com/CSU-KangHu/NeuralTE.git" 2>/dev/null \
        || warn "NeuralTE GitHub install failed — install from https://github.com/CSU-KangHu/NeuralTE"

    # TERL — LSTM-based TE classifier
    # Paper: Lopez-Escamilla et al. (2021); doi:10.1093/molbev/msab143
    python3 -m pip install -q \
        "git+https://github.com/simonorozcoarias/TERL.git" 2>/dev/null \
        || warn "TERL GitHub install failed — install from https://github.com/simonorozcoarias/TERL"

    # Terrier — k-mer-based TE classifier
    # Paper: Ruiz et al. (2022); doi:10.1093/bioinformatics/btac274
    python3 -m pip install -q terrier-te 2>/dev/null \
        || warn "Terrier pip install failed — try: pip install terrier-te"

    # ClassifyTE — SVM-based TE classifier
    # Paper: Orozco-Arias et al. (2021); doi:10.1093/bioinformatics/btaa1101
    python3 -m pip install -q \
        "git+https://github.com/simonorozcoarias/ClassifyTE.git" 2>/dev/null \
        || warn "ClassifyTE GitHub install failed — install from https://github.com/simonorozcoarias/ClassifyTE"

    # CREATE — Random-forest TE classifier
    # Paper: Orozco-Arias et al. (2022); doi:10.3390/genes13020239
    python3 -m pip install -q \
        "git+https://github.com/simonorozcoarias/create.git" 2>/dev/null \
        || warn "CREATE GitHub install failed — install from https://github.com/simonorozcoarias/CREATE"

    # TEClass2 — kNN/SVM hierarchical TE classifier
    # Paper: Nicolás et al. (2021)
    install_conda "teclass2"    "-c bioconda -c conda-forge" \
        || python3 -m pip install -q teclass2 2>/dev/null \
        || warn "TEClass2 install failed — install from https://github.com/HelixPG/TEClass2"
fi

# ---------------------------------------------------------------------------
# 6. Build MycoSV binary
# ---------------------------------------------------------------------------
if [[ "$MODE" =~ ^(all|mycosv-only)$ ]]; then
    echo ""
    echo -e "${YELLOW}[6] Building MycoSV binary (fungi_graphsv_tol)${NC}"
    BINARY="${SCRIPT_DIR}/fungi_graphsv_tol"

    # Detect g++ from environment
    GPP="$(command -v g++ 2>/dev/null || command -v c++ 2>/dev/null || echo "")"
    if [[ -z "$GPP" ]]; then
        fail "g++ not found — install via: conda install -c conda-forge gxx_linux-64"
        exit 1
    fi

    info "Compiling with $GPP ..."
    "${GPP}" -O2 -DNDEBUG -std=c++17 -pthread \
        -I"${SCRIPT_DIR}" \
        "${SCRIPT_DIR}/main.cpp" \
        -o "${BINARY}" \
        && ok "Built: ${BINARY}" \
        || { fail "Build failed — check compiler errors above"; exit 1; }
fi

# ---------------------------------------------------------------------------
# 7. Verification
# ---------------------------------------------------------------------------
echo ""
echo -e "${YELLOW}[7] Verification${NC}"
echo ""
echo "=== MycoSV binary ==="
check_tool "fungi_graphsv_tol" || check_tool "${SCRIPT_DIR}/fungi_graphsv_tol" || true
echo ""

if [[ "$MODE" =~ ^(all|sv-only)$ ]]; then
    echo "=== Assembly SV tools ==="
    check_tool minimap2 || true
    check_tool syri     || true
    check_tool minigraph|| true
    check_tool gfatools  || true
    check_tool pggb     || true
    check_tool svim-asm || true
    check_tool anchorwave|| true
    echo ""
    echo "=== Short-reads SV tools ==="
    check_tool delly          || true
    check_tool configManta.py || true
    echo ""
    echo "=== Long-reads SV tools ==="
    check_tool svim     || true
    check_tool sniffles || true
    check_tool cuteSV   || true
    echo ""
fi

if [[ "$MODE" =~ ^(all|te-only)$ ]]; then
    echo "=== TE classification tools ==="
    check_tool DeepTE.py   || true
    check_tool NeuralTE.py || true
    check_tool TERL.py     || true
    check_tool terrier.py  || true
    check_tool ClassifyTE.py|| true
    check_tool CREATE.py   || true
    check_tool TEClass2    || true
    check_python_module tensorflow || true
    check_python_module torch      || true
    echo ""
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo -e "${BLUE}========================================"
echo "Install Summary"
echo "========================================${NC}"
echo "  Available: ${#AVAILABLE_TOOLS[@]}"
echo "  Missing  : ${#MISSING_TOOLS[@]}"
if [[ ${#MISSING_TOOLS[@]} -gt 0 ]]; then
    echo ""
    warn "The following tools could not be installed automatically:"
    for t in "${MISSING_TOOLS[@]}"; do
        echo "    - $t"
    done
    echo ""
    echo "To activate the environment and retry:"
    echo "  conda activate ${ENV_NAME}"
    echo "  bash install_tools.sh --check"
fi
echo ""
echo "To activate the environment:"
echo "  conda activate ${ENV_NAME}"
echo ""
echo "To run experiments:"
echo "  bash run_all_experiments.sh --real"
echo "  python3 run_te_benchmark.py --download-fungi-demo --out-dir te_results/"
echo ""
echo -e "${GREEN}Done.${NC}"
