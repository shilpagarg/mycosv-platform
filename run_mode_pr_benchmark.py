#!/usr/bin/env python3
# Designed for Linux

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

from mycosv_cli_runner import run_mycosv_command
from sv_pr_utils import score_pr, write_pr_artifacts


ROOT = Path(__file__).resolve().parent
DEFAULT_SIM = ROOT / "test_amf.py"
DEFAULT_MAIN = ROOT / "main.cpp"
DEFAULT_BIN = ROOT / "fungi_graphsv_tol_bin"


def run(cmd: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, cwd=str(cwd) if cwd else None, text=True, capture_output=True, check=True)


def compile_binary(main_cpp: Path, binary_path: Path) -> None:
    binary_path.parent.mkdir(parents=True, exist_ok=True)
    run(["g++", "-O2", "-std=c++17", "-pthread", str(main_cpp), "-o", str(binary_path)], cwd=ROOT)


def mode_specific_caller_args(mode: str, genome_size_hint: int) -> list[str]:
    args: list[str] = []
    if mode != "assembly" and genome_size_hint > 0:
        args.extend(["--genome-size-hint", str(genome_size_hint)])
    if mode == "short-reads":
        args.extend(["--sr-min-unitig-len", "100", "--sr-min-kmer-freq", "2"])
    return args




def ensure_sv_volume(n_genomes: int, n_reps: int, n_contigs: int, total_len: int,
                     scenario_set: str, target_svs_per_scenario: int,
                     min_contig_bp: int) -> tuple[int, int, int]:
    scenarios = [s.strip() for s in scenario_set.split(",") if s.strip()]
    n_scen = max(1, len(scenarios))
    n_reps = max(n_reps, n_scen)
    query_genomes = max(0, n_genomes - n_reps)
    per_query = max(1, n_contigs)
    current = max(0, (query_genomes // n_scen) * per_query)
    if current < target_svs_per_scenario:
        needed_queries_per_scenario = (target_svs_per_scenario + per_query - 1) // per_query
        n_genomes = max(n_genomes, n_reps + needed_queries_per_scenario * n_scen)
    per_contig_bp = max(1, total_len // max(1, n_contigs))
    if per_contig_bp < min_contig_bp:
        total_len = min_contig_bp * max(1, n_contigs)
    return n_genomes, n_reps, total_len
def parse_modes(raw: str) -> list[str]:
    modes = [m.strip() for m in raw.split(",") if m.strip()]
    valid = {"assembly", "short-reads", "long-reads"}
    for mode in modes:
        if mode not in valid:
            raise ValueError(f"Unsupported mode {mode!r}; valid={sorted(valid)}")
    if not modes:
        raise ValueError("At least one mode is required")
    return modes


def write_mode_summary(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "mode", "truth_records", "pred_records_algo", "hint_leaked_calls",
                "tp", "fp", "fn", "precision", "recall", "f1",
                "prec_lo95", "prec_hi95", "rec_lo95", "rec_hi95",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_svtype_summary(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "mode", "svtype", "tp", "fp", "fn",
                "precision", "prec_lo95", "prec_hi95",
                "recall", "rec_lo95", "rec_hi95", "f1",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description="One-command assembly/short-read/long-read SV precision/recall benchmark runner.")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--simulator", type=Path, default=DEFAULT_SIM)
    ap.add_argument("--main-cpp", type=Path, default=DEFAULT_MAIN)
    ap.add_argument("--binary-path", type=Path, default=DEFAULT_BIN)
    ap.add_argument("--skip-build", action="store_true")
    ap.add_argument("--modes", default="assembly,short-reads,long-reads")
    ap.add_argument("--phylum", default="Ascomycota")
    ap.add_argument("--scenario-set", default="compact_yeast,pathogenic,arbuscular_mf,cross_phylum_hgt_stress")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--n-genomes", type=int, default=4,
                    help="Number of simulated genomes per run (legacy alias: --n-refs).")
    ap.add_argument("--n-reps", type=int, default=2,
                    help="Number of simulated replicates / queries per run (legacy alias: --n-queries).")
    ap.add_argument("--n-refs", dest="n_genomes", type=int,
                    help="Backward-compatible alias for --n-genomes.")
    ap.add_argument("--n-queries", dest="n_reps", type=int,
                    help="Backward-compatible alias for --n-reps.")
    ap.add_argument("--total-len", type=int, default=12000)
    ap.add_argument("--n-contigs", type=int, default=2)
    ap.add_argument("--divergence", type=float, default=0.01)
    ap.add_argument("--routing-top-n", type=int, default=4)
    ap.add_argument("--min-svlen", type=int, default=40)
    ap.add_argument("--window-bp", type=int, default=2_000_000)
    ap.add_argument("--target-svs-per-scenario", type=int, default=3000)
    ap.add_argument("--min-contig-bp", type=int, default=12000)
    args = ap.parse_args()

    modes = parse_modes(args.modes)
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    args.n_genomes, args.n_reps, args.total_len = ensure_sv_volume(
        args.n_genomes,
        args.n_reps,
        args.n_contigs,
        args.total_len,
        args.scenario_set,
        args.target_svs_per_scenario,
        args.min_contig_bp,
    )

    if not args.skip_build or not args.binary_path.exists():
        compile_binary(args.main_cpp.resolve(), args.binary_path.resolve())

    mode_rows: list[dict[str, object]] = []
    svtype_rows: list[dict[str, object]] = []
    summary_json: dict[str, object] = {
        "modes": {},
        "config": {
            "phylum": args.phylum,
            "scenario_set": args.scenario_set,
            "n_genomes": args.n_genomes,
            "n_reps": args.n_reps,
            "total_len": args.total_len,
            "n_contigs": args.n_contigs,
            "divergence": args.divergence,
            "target_svs_per_scenario": args.target_svs_per_scenario,
            "min_contig_bp": args.min_contig_bp,
            "routing_top_n": args.routing_top_n,
            "min_svlen": args.min_svlen,
            "window_bp": args.window_bp,
        },
    }

    for mode in modes:
        mode_dir = out_dir / mode
        sim_dir = mode_dir / "sim"
        sim_dir.mkdir(parents=True, exist_ok=True)
        sim_cmd = [
            sys.executable, str(args.simulator.resolve()),
            "--phylum", args.phylum,
            "--n-genomes", str(args.n_genomes),
            "--n-reps", str(args.n_reps),
            "--total-len", str(args.total_len),
            "--n-contigs", str(args.n_contigs),
            "--out-dir", str(sim_dir),
            "--scenario-set", args.scenario_set,
            "--write-extended-manifest",
            "--query-mode", mode,
            "--divergence", str(args.divergence),
            "--seed", str(args.seed),
        ]
        sim_run = run(sim_cmd, cwd=ROOT)

        out_prefix = mode_dir / "calls"
        caller_cmd = [
            str(args.binary_path.resolve()),
            "--ref-list", str(sim_dir / "ref_list.txt"),
            "--query-list", str(sim_dir / "query_list.txt"),
            "--out-prefix", str(out_prefix),
            "--query-mode", mode,
            "--routing-top-n", str(args.routing_top_n),
            "--min-svlen", str(args.min_svlen),
            *mode_specific_caller_args(mode, args.total_len),
        ]
        caller_run = run_mycosv_command(caller_cmd, cwd=ROOT)

        summary = score_pr(
            sim_dir / "truth" / "all_queries.truth.ref.vcf",
            out_prefix.with_suffix(".vcf"),
            pred_hits_tsv=out_prefix.with_suffix(".hits.tsv"),
            query_meta_tsv=sim_dir / "query_metadata.tsv",
            window_bp=args.window_bp,
        )
        write_pr_artifacts(summary, mode_dir / "pr_metrics.tsv", mode_dir / "pr_metrics.json")
        overall = summary["overall"]
        mode_rows.append({
            "mode": mode,
            "truth_records": summary["truth_records"],
            "pred_records_algo": summary["pred_records_algo"],
            "hint_leaked_calls": summary["hint_leaked_calls"],
            "tp": overall["tp"],
            "fp": overall["fp"],
            "fn": overall["fn"],
            "precision": overall["precision"],
            "recall": overall["recall"],
            "f1": overall["f1"],
            "prec_lo95": overall["precision_ci95"][0],
            "prec_hi95": overall["precision_ci95"][1],
            "rec_lo95": overall["recall_ci95"][0],
            "rec_hi95": overall["recall_ci95"][1],
        })
        for row in summary["per_type_rows"]:
            svtype_rows.append({"mode": mode, **row})

        summary_json["modes"][mode] = {
            "simulation_stdout": sim_run.stdout,
            "caller_stdout": caller_run.stdout,
            "caller_stderr": caller_run.stderr,
            "metrics": {k: v for k, v in summary.items() if k not in {"tp_records", "fp_records", "fn_records"}},
        }

    write_mode_summary(out_dir / "mode_pr_summary.tsv", mode_rows)
    write_svtype_summary(out_dir / "mode_svtype_pr_summary.tsv", svtype_rows)
    with (out_dir / "mode_pr_summary.json").open("w") as fh:
        json.dump(summary_json, fh, indent=2, sort_keys=True)

    for row in mode_rows:
        print(
            f"{row['mode']}\tTP={row['tp']}\tFP={row['fp']}\tFN={row['fn']}\t"
            f"precision={float(row['precision']):.4f}\trecall={float(row['recall']):.4f}\tf1={float(row['f1']):.4f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
