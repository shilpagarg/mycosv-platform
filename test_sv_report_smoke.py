#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPORT = ROOT / "sv_visualization_report.py"


def run(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, text=True, capture_output=True, check=True)


def test_report_builds_from_minimal_inputs(tmp_path: Path) -> None:
    sim = tmp_path / "sim.tsv"
    sim.write_text(
        "caller\tdataset\tsv_type\tprecision\trecall\tf1\n"
        "mycosv\tsim1\tDEL\t1.0\t0.9\t0.947\n",
        encoding="utf-8",
    )

    real = tmp_path / "real.tsv"
    real.write_text(
        "sample\tsv_type\tchrom\tsv_len\n"
        "s1\tDEL\tchr1\t120\n"
        "s1\tINS\tchr2\t80\n",
        encoding="utf-8",
    )

    bio = tmp_path / "bio.tsv"
    bio.write_text(
        "gene\tpathway\teffect\n"
        "ERG11\tsterol\tgain\n",
        encoding="utf-8",
    )

    outdir = tmp_path / "report"
    run([
        sys.executable, str(REPORT),
        "--simulated", str(sim),
        "--real", str(real),
        "--biology", str(bio),
        "--outdir", str(outdir),
        "--title", "Smoke report",
    ])

    html = outdir / "sv_visualization_report.html"
    assert html.exists()
    text = html.read_text(encoding="utf-8")
    assert "Simulated data benchmarking" in text
    assert "Real data structural variant analyses" in text
    assert "Biological findings" in text
