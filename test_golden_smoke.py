#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
BIN = ROOT / "fungi_graphsv_tol_bin"
SIM = ROOT / "test_amf.py"
MAIN = ROOT / "main.cpp"


def run(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, text=True, capture_output=True, check=True)


def ensure_binary() -> None:
    if not BIN.exists():
        run(["g++", "-O2", "-std=c++17", "-pthread", str(MAIN), "-o", str(BIN)])


def test_seeded_simulation_is_reproducible(tmp_path: Path) -> None:
    out1 = tmp_path / "run1"
    out2 = tmp_path / "run2"

    base = [
        sys.executable, str(SIM),
        "--phylum", "Ascomycota",
        "--n-genomes", "4",
        "--n-reps", "2",
        "--total-len", "12000",
        "--n-contigs", "2",
        "--scenario-set", "compact_yeast",
        "--seed", "42",
        "--write-extended-manifest",
    ]

    run(base + ["--out-dir", str(out1)])
    run(base + ["--out-dir", str(out2)])

    assert (out1 / "query_truth.tsv").read_text() == (out2 / "query_truth.tsv").read_text()
    assert (out1 / "query_metadata.tsv").read_text() == (out2 / "query_metadata.tsv").read_text()


def test_cross_mode_smoke_produces_vcf(tmp_path: Path) -> None:
    ensure_binary()

    ref = tmp_path / "ref.fa"
    ref.write_text(">ctg1\n" + "ACGT" * 150 + "\n")
    refs = tmp_path / "refs.txt"
    refs.write_text(str(ref) + "\n")

    asm = tmp_path / "query.fa"
    asm.write_text(">ctg1\n" + "ACGT" * 120 + "\n")
    asm_q = tmp_path / "asm.txt"
    asm_q.write_text(str(asm) + "\n")

    for mode, qlist in [("assembly", asm_q)]:
        out = tmp_path / f"{mode}_out"
        run([str(BIN), "--ref-list", str(refs), "--query-list", str(qlist), "--out-prefix", str(out), "--query-mode", mode])
        assert Path(str(out) + ".vcf").exists()
