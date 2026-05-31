#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${ROOT_DIR}"

ARCHIVE_DIR="cleanup_archive"
STAMP="$(date +%Y%m%d)"
mkdir -p "${ARCHIVE_DIR}"

shopt -s nullglob

remove_path() {
    local path="$1"
    local attempt

    for attempt in 1 2 3; do
        [[ -e "${path}" ]] || return 0
        rm -rf "${path}" 2>/dev/null && return 0
        sleep 1
    done

    rm -rf "${path}"
}

slurm_logs=(slurm-*.out)
if ((${#slurm_logs[@]})); then
    archive="${ARCHIVE_DIR}/slurm_logs_${STAMP}.tar.gz"
    suffix=1
    while [[ -e "${archive}" ]]; do
        archive="${ARCHIVE_DIR}/slurm_logs_${STAMP}_${suffix}.tar.gz"
        ((suffix += 1))
    done

    tar -czf "${archive}" "${slurm_logs[@]}"
    rm -f "${slurm_logs[@]}"
    echo "Archived ${#slurm_logs[@]} Slurm log(s) to ${archive}"
fi

remove_path __pycache__
remove_path .pytest_cache
rm -f fungi_graphsv_tol_bin fungi_graphsv_tol_bin.b1 fungi_graphsv_tol_bin.lock

echo "Workspace cleanup complete."
