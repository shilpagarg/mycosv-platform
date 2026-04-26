#!/usr/bin/env python3
# Designed for Linux

from __future__ import annotations

import argparse
import csv
import json
import struct
import subprocess
import sys
import time
from pathlib import Path

from mycosv_cli_runner import run_mycosv_command
from sv_pr_utils import expand_to_multisample_vcf, score_pr, write_pr_artifacts


ROOT = Path(__file__).resolve().parent
DEFAULT_SIM = ROOT / "test_amf.py"
DEFAULT_MAIN = ROOT / "main.cpp"
DEFAULT_BIN = ROOT / "fungi_graphsv_tol_bin"
DEFAULT_ANALYZE = ROOT / "analyze_new_biology_candidates.py"


def run(cmd: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, cwd=str(cwd) if cwd else None, text=True, capture_output=True, check=True, timeout=7200)


def compile_binary(main_cpp: Path, binary_path: Path) -> None:
    binary_path.parent.mkdir(parents=True, exist_ok=True)
    run(["g++", "-O2", "-std=c++17", "-pthread", str(main_cpp), "-o", str(binary_path)], cwd=ROOT)


def parse_modes(raw: str) -> list[str]:
    modes = [m.strip() for m in raw.split(",") if m.strip()]
    valid = {"assembly", "short-reads", "long-reads"}
    for mode in modes:
        if mode not in valid:
            raise ValueError(f"Unsupported mode {mode!r}; valid={sorted(valid)}")
    if not modes:
        raise ValueError("At least one mode is required")
    return modes


def mode_specific_caller_args(mode: str, genome_size_hint: int) -> list[str]:
    args: list[str] = []
    if mode != "assembly" and genome_size_hint > 0:
        args.extend(["--genome-size-hint", str(genome_size_hint)])
    if mode == "short-reads":
        args.extend([
            "--sr-kmer-size", "21",
            "--sr-min-kmer-freq", "2",
            "--min-anchors-per-block", "3",
        ])
    if mode == "long-reads":
        args.extend([
            "--lr-anchor-k", "15",       # shorter k for long-read error tolerance
            "--chain-gap-band", "15000",  # wider band for large read-length gaps
            "--min-anchors-per-block", "1",
            "--min-block-score", "3.0",
        ])
    return args




def ensure_sv_volume(n_genomes: int, n_reps: int, n_contigs: int, total_len: int,
                     scenario_set: str, target_svs_per_scenario: int,
                     min_contig_bp: int,
                     queries_per_scenario: int | None = None) -> tuple[int, int, int]:
    scenarios = [s.strip() for s in scenario_set.split(",") if s.strip()]
    n_scen = max(1, len(scenarios))
    n_reps = max(n_reps, n_scen)
    if queries_per_scenario is not None and queries_per_scenario > 0:
        # Explicit per-scenario query budget overrides the SV-volume heuristic.
        # Useful for fast testing: e.g. 20 queries × 3 scenarios = 60 query genomes.
        n_genomes = n_reps + queries_per_scenario * n_scen
    else:
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
def splitmix64(x: int) -> int:
    x = (x + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
    z = x
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
    return z ^ (z >> 31)


def make_hashes(n: int, seed: int) -> list[int]:
    hashes: list[int] = []
    seen: set[int] = set()
    cur = seed & 0xFFFFFFFFFFFFFFFF
    while len(hashes) < n:
        cur = splitmix64(cur)
        if cur in seen:
            continue
        seen.add(cur)
        hashes.append(cur)
    hashes.sort()
    return hashes


def parse_routing_manifest(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        hashes = [int(tok) for tok in cols[10].split(",") if tok]
        rows.append(
            {
                "clade_name": cols[0],
                "clade_rank": cols[1],
                "phylum": cols[2],
                "hashes": hashes,
            }
        )
    return rows


def write_disk_record(fh, clade_name: str, phylum: str, clade_rank: str, hashes: list[int]) -> None:
    def write_string(text: str) -> None:
        raw = text.encode("utf-8")
        fh.write(struct.pack("<I", len(raw)))
        fh.write(raw)

    write_string(clade_name)
    write_string(phylum)
    write_string(clade_rank)
    fh.write(struct.pack("<I", len(hashes)))
    if hashes:
        fh.write(struct.pack(f"<{len(hashes)}Q", *hashes))


def augment_routing_store(index_dir: Path, target_centroids: int, seed: int) -> dict[str, int]:
    routing_manifest = index_dir / "routing_manifest.tsv"
    store_path = index_dir / "routing_centroids.bin"
    real_rows = parse_routing_manifest(routing_manifest)
    if not real_rows:
        raise RuntimeError(f"No routing rows found in {routing_manifest}")
    if target_centroids < len(real_rows):
        raise ValueError(
            f"target_centroids={target_centroids} is smaller than real routing rows={len(real_rows)}"
        )

    hashes_per_centroid = max(16, max(len(row["hashes"]) for row in real_rows if row["hashes"]))
    decoys = target_centroids - len(real_rows)
    with store_path.open("wb") as fh:
        fh.write(struct.pack("<Q", target_centroids))
        for row in real_rows:
            write_disk_record(
                fh,
                str(row["clade_name"]),
                str(row["phylum"]),
                str(row["clade_rank"]),
                list(row["hashes"]),
            )
        for i in range(decoys):
            phylum = f"DecoyPhylum_{i % 16}"
            rank_cycle = ("phylum", "class", "order", "family", "genus", "species")
            rank = rank_cycle[i % len(rank_cycle)]
            hashes = make_hashes(hashes_per_centroid, seed ^ (i * 0x9E3779B97F4A7C15))
            write_disk_record(fh, f"decoy_clade_{i}", phylum, rank, hashes)

    skip_path = Path(str(store_path) + ".skip")
    if skip_path.exists():
        skip_path.unlink()

    return {
        "real_centroids": len(real_rows),
        "decoy_centroids": decoys,
        "total_centroids": target_centroids,
        "hashes_per_centroid": hashes_per_centroid,
    }


def write_summary(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "mode",
                "n_centroids",
                "truth_records",
                "pred_records_algo",
                "tp",
                "fp",
                "fn",
                "precision",
                "recall",
                "f1",
                "query_seconds",
                "skip_index_bytes",
                "store_bytes",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Run assembly/short-read/long-read query benchmarks against a million-scale hierarchical routing catalog."
    )
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--simulator", type=Path, default=DEFAULT_SIM)
    ap.add_argument("--main-cpp", type=Path, default=DEFAULT_MAIN)
    ap.add_argument("--binary-path", type=Path, default=DEFAULT_BIN)
    ap.add_argument("--skip-build", action="store_true")
    ap.add_argument("--modes", default="assembly,short-reads,long-reads")
    ap.add_argument("--n-centroids", type=int, default=1_000_000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--phylum", default="Ascomycota")
    ap.add_argument("--scenario-set", default="compact_yeast,two_speed_pathogen_extreme,arbuscular_mf")
    ap.add_argument("--n-genomes", type=int, default=4)
    ap.add_argument("--n-reps", type=int, default=2)
    ap.add_argument("--total-len", type=int, default=8000)
    ap.add_argument("--n-contigs", type=int, default=2)
    ap.add_argument("--divergence", type=float, default=0.01)
    ap.add_argument("--routing-top-n", type=int, default=4)
    ap.add_argument("--min-svlen", type=int, default=40)
    ap.add_argument("--window-bp", type=int, default=2_000_000)
    ap.add_argument("--target-svs-per-scenario", type=int, default=3000)
    ap.add_argument("--queries-per-scenario", type=int, default=0,
                    help="If >0, set the number of query genomes per scenario directly "
                         "(n_genomes = n_reps + queries_per_scenario * n_scenarios). "
                         "Overrides --target-svs-per-scenario. Useful for fast testing, "
                         "e.g. --queries-per-scenario 20.")
    ap.add_argument("--min-contig-bp", type=int, default=12000)
    ap.add_argument("--threads", type=int, default=32,
                    help="Parallel worker threads passed to the MycoSV binary (--threads) "
                         "and index build (--tol-index-threads). Default: 32.")
    ap.add_argument("--resume", action="store_true",
                    help="Skip simulation/build/query for any mode where calls.vcf already "
                         "exists in the output directory, jumping straight to PR scoring. "
                         "Useful after a session disconnect that killed Python mid-run.")
    ap.add_argument("--long-read-platform", default="hifi,ont-r10",
                    help="Comma-separated list of long-read platform presets to run. "
                         "Each generates a separate long-reads_<platform>/ output directory. "
                         "hifi=PacBio HiFi CCS (15 kb, ≥Q20), "
                         "ont-r10=ONT R10.4.1 simplex (10 kb, ~Q20), "
                         "ont-r9=ONT R9.4.1 (8 kb, ~Q15), "
                         "generic=simulator defaults. "
                         "Default: hifi,ont-r10 (both platforms).")
    args = ap.parse_args()

    modes = parse_modes(args.modes)
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Validate long-read platforms.
    _valid_lr_platforms = {"hifi", "ont-r10", "ont-r9", "generic"}
    lr_platforms = [p.strip() for p in args.long_read_platform.split(",") if p.strip()]
    for p in lr_platforms:
        if p not in _valid_lr_platforms:
            raise ValueError(f"Unknown long-read platform {p!r}; valid: {sorted(_valid_lr_platforms)}")

    # Expand modes: long-reads → one iteration per platform; others run once.
    # label      = output directory name and summary key (e.g. "long-reads_hifi")
    # actual_mode = query-mode flag passed to simulator/binary ("long-reads")
    # lr_plat    = platform preset forwarded to the simulator (empty for non-lr)
    run_tasks: list[tuple[str, str, str]] = []
    for mode in modes:
        if mode == "long-reads":
            for plat in lr_platforms:
                run_tasks.append((f"long-reads_{plat}", "long-reads", plat))
        else:
            run_tasks.append((mode, mode, ""))

    args.n_genomes, args.n_reps, args.total_len = ensure_sv_volume(
        args.n_genomes,
        args.n_reps,
        args.n_contigs,
        args.total_len,
        args.scenario_set,
        args.target_svs_per_scenario,
        args.min_contig_bp,
        queries_per_scenario=(args.queries_per_scenario if args.queries_per_scenario > 0 else None),
    )

    if not args.skip_build or not args.binary_path.exists():
        compile_binary(args.main_cpp.resolve(), args.binary_path.resolve())

    summary_rows: list[dict[str, object]] = []
    summary_json: dict[str, object] = {
        "config": {
            "n_centroids": args.n_centroids,
            "phylum": args.phylum,
            "scenario_set": args.scenario_set,
            "n_genomes": args.n_genomes,
            "n_reps": args.n_reps,
            "total_len": args.total_len,
            "n_contigs": args.n_contigs,
            "divergence": args.divergence,
            "target_svs_per_scenario": args.target_svs_per_scenario,
            "queries_per_scenario": args.queries_per_scenario,
            "min_contig_bp": args.min_contig_bp,
            "long_read_platforms": lr_platforms,
        },
        "modes": {},
    }

    for label, actual_mode, lr_plat in run_tasks:
        mode_dir = out_dir / label
        sim_dir = mode_dir / "sim"
        idx_dir = mode_dir / "idx"
        reg_dir = mode_dir / "reg"
        sim_dir.mkdir(parents=True, exist_ok=True)
        idx_dir.mkdir(parents=True, exist_ok=True)
        reg_dir.mkdir(parents=True, exist_ok=True)

        out_prefix = mode_dir / "calls"
        resuming = args.resume and out_prefix.with_suffix(".vcf").exists()

        if resuming:
            # Session was disconnected mid-run; the binary already wrote its output.
            # Reconstruct the minimal state needed for PR scoring and summary.
            print(f"[resume] {label}: calls.vcf exists — skipping sim/build/query", flush=True)
            store_path = idx_dir / "routing_centroids.bin"
            store_info = {
                "real_centroids": 0,
                "decoy_centroids": 0,
                "total_centroids": args.n_centroids,
                "hashes_per_centroid": 0,
            }
            if (idx_dir / "routing_manifest.tsv").exists():
                try:
                    real_rows = parse_routing_manifest(idx_dir / "routing_manifest.tsv")
                    store_info["real_centroids"] = len(real_rows)
                    store_info["decoy_centroids"] = args.n_centroids - len(real_rows)
                except Exception:
                    pass
            sim_run = build_run = query_run = None  # type: ignore[assignment]
            query_seconds = 0.0
        else:
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
                "--query-mode", actual_mode,
                "--divergence", str(args.divergence),
                # Platform preset for long-reads: controls read length, error rate,
                # and coverage in the simulator (see _LR_PLATFORM_PRESETS in test_amf.py).
                # hifi    → 15 kb, ≥Q20 — minimap2 map-hifi downstream
                # ont-r10 → 10 kb, ~Q20 — minimap2 map-ont; WhatsHap applicable for
                #           dikaryotic fungi (Puccinia, Leptosphaeria, Zymoseptoria)
                # ont-r9  → 8 kb, ~Q15
                *( ["--long-read-platform", lr_plat] if lr_plat else [] ),
            ]
            sim_run = run(sim_cmd, cwd=ROOT)

            build_cmd = [
                str(args.binary_path.resolve()),
                "--tol-hierarchical",
                "--tol-build-index", str(sim_dir / "hierarchy_manifest.tsv"),
                "--tol-index-dir", str(idx_dir),
                "--tol-registry-dir", str(reg_dir),
                "--tol-multi-rank",
                "--tol-base-graph-build",
                "--tol-max-clade-genomes", "32",
                "--tol-index-threads", str(args.threads),
            ]
            build_run = run_mycosv_command(build_cmd, cwd=ROOT)
            store_info = augment_routing_store(idx_dir, args.n_centroids, args.seed)

            query_cmd = [
                str(args.binary_path.resolve()),
                "--tol-hierarchical",
                "--tol-index-dir", str(idx_dir),
                "--tol-registry-dir", str(reg_dir),
                "--ref-list", str(sim_dir / "ref_list.txt"),
                "--query-list", str(sim_dir / "query_list.txt"),
                "--out-prefix", str(out_prefix),
                "--query-mode", actual_mode,
                "--routing-top-n", str(args.routing_top_n),
                "--min-svlen", str(args.min_svlen),
                "--threads", str(args.threads),
                *mode_specific_caller_args(actual_mode, args.total_len),
            ]
            q_start = time.perf_counter()
            query_run = run_mycosv_command(query_cmd, cwd=ROOT)
            query_seconds = time.perf_counter() - q_start

        # The C++ binary emits calls.vcf as a single-sample VCF with per-query
        # provenance only in INFO (QASM=). Produce a sibling multi-sample VCF
        # that materializes one column per query asm so downstream tools (and
        # spot-checks) see all queries explicitly.
        try:
            expand_to_multisample_vcf(
                out_prefix.with_suffix(".vcf"),
                out_prefix.parent / (out_prefix.name + ".multisample.vcf"),
            )
        except Exception as exc:
            sys.stderr.write(f"[multisample-vcf] expand failed: {exc}\n")

        summary = score_pr(
            sim_dir / "truth" / "all_queries.truth.ref.vcf",
            out_prefix.with_suffix(".vcf"),
            pred_hits_tsv=out_prefix.with_suffix(".hits.tsv"),
            query_meta_tsv=sim_dir / "query_metadata.tsv",
            window_bp=args.window_bp,
        )
        write_pr_artifacts(summary, mode_dir / "pr_metrics.tsv", mode_dir / "pr_metrics.json")
        overall = summary["overall"]

        # Biology candidate analysis — classifies TE-linked / novel SVs and
        # produces biology_candidates.tsv picked up by sv_visualization_report.py.
        if DEFAULT_ANALYZE.exists():
            bio_dir = mode_dir / "biology"
            bio_dir.mkdir(parents=True, exist_ok=True)
            try:
                run([
                    sys.executable, str(DEFAULT_ANALYZE),
                    "--vcf",            str(out_prefix.with_suffix(".vcf")),
                    "--hits",           str(out_prefix.with_suffix(".hits.tsv")),
                    "--query-metadata", str(sim_dir / "query_metadata.tsv"),
                    "--phylum",         args.phylum,
                    "--out-tsv",        str(bio_dir / "biology_candidates.tsv"),
                    "--summary-json",   str(bio_dir / "biology_candidates.json"),
                    "--top-n",          "200",
                ], cwd=ROOT)
            except subprocess.CalledProcessError:
                pass  # biology analysis is non-fatal

        store_path = idx_dir / "routing_centroids.bin"
        skip_path = Path(str(store_path) + ".skip")
        summary_rows.append(
            {
                "mode": label,
                "n_centroids": store_info["total_centroids"],
                "truth_records": summary["truth_records"],
                "pred_records_algo": summary["pred_records_algo"],
                "tp": overall["tp"],
                "fp": overall["fp"],
                "fn": overall["fn"],
                "precision": overall["precision"],
                "recall": overall["recall"],
                "f1": overall["f1"],
                "query_seconds": query_seconds,
                "skip_index_bytes": skip_path.stat().st_size if skip_path.exists() else 0,
                "store_bytes": store_path.stat().st_size if store_path.exists() else 0,
            }
        )
        summary_json["modes"][label] = {
            "simulation_stdout": sim_run.stdout if sim_run else "(resumed)",
            "build_stdout": build_run.stdout if build_run else "(resumed)",
            "query_stdout": query_run.stdout if query_run else "(resumed)",
            "query_stderr": query_run.stderr if query_run else "(resumed)",
            "store_info": store_info,
            "metrics": {k: v for k, v in summary.items() if k not in {"tp_records", "fp_records", "fn_records"}},
            "query_seconds": query_seconds,
        }

    write_summary(out_dir / "million_mode_summary.tsv", summary_rows)
    with (out_dir / "million_mode_summary.json").open("w", encoding="utf-8") as fh:
        json.dump(summary_json, fh, indent=2, sort_keys=True)

    for row in summary_rows:
        print(
            f"{row['mode']}\tn={row['n_centroids']}\tTP={row['tp']}\tFP={row['fp']}\tFN={row['fn']}\t"
            f"precision={float(row['precision']):.4f}\trecall={float(row['recall']):.4f}\t"
            f"f1={float(row['f1']):.4f}\tquery_seconds={float(row['query_seconds']):.2f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
