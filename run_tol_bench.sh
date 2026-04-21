#!/usr/bin/env bash
set -euo pipefail

RUN_ALL_MODES="${RUN_ALL_MODES:-0}"
QUERY_MODE="${QUERY_MODE:-assembly}"

run_all_modes(){
    local root_base="$1"
    local modes=(assembly short-reads long-reads)
    for mode in "${modes[@]}"; do
        RUN_ALL_MODES=0 QUERY_MODE="$mode" bash "$0" "$TEST_AMF" "$MAIN_CPP" "${root_base}_${mode}"
    done
}

TEST_AMF="${1:-test_amf.py}"
MAIN_CPP="${2:-main.cpp}"
ROOT_BASE="${3:-bench}"

if [[ "$RUN_ALL_MODES" == "1" ]]; then
    run_all_modes "$ROOT_BASE"
    exit 0
fi

echo "run_tol_bench.sh placeholder mode=$QUERY_MODE root=$ROOT_BASE"
