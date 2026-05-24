#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "plot_mycosv_pangenome_calls.py"


def test_pangenome_plot_script_builds_summary_and_html(tmp_path: Path) -> None:
    bench = tmp_path / "bench"
    bench.mkdir()
    (bench / "novel_mycosv_calls.tsv").write_text(
        "query_asm\tquery_contig\tpos\tend\tsvtype\tsvlen\telement_class\t"
        "single_reference_equivalent\tread_supported\tdiscovery_bucket\n"
        "q1\tctg1\t10\t20\tINS\t10\tHGT\tno\tyes\tpangenome_only_read_supported\n"
        "q1\tctg1\t30\t40\tDEL\t-10\tRIP\tyes\tyes\tsingle_reference_equivalent\n",
        encoding="utf-8",
    )
    (bench / "pangenome_call_layers.tsv").write_text(
        "query_asm\tquery_mode\traw_pairwise_pangenome_observations\t"
        "deduplicated_pangenome_loci\tsingle_reference_equivalent_calls\t"
        "pangenome_only_calls\tpangenome_only_read_supported\t"
        "pangenome_only_intrinsic_supported\traw_to_deduplicated_ratio\t"
        "single_ref_fraction_of_raw\n"
        "ALL\tall\t2\t2\t1\t1\t1\t0\t1.0\t0.5\n",
        encoding="utf-8",
    )
    (bench / "mycosv_evidence_tiers.tsv").write_text(
        "query_asm\tsvtype\ttier\tn_calls\nq1\tINS\tstrong\t1\n",
        encoding="utf-8",
    )

    out = tmp_path / "plots"
    subprocess.run(
        [sys.executable, str(SCRIPT), "--benchmark-dir", str(bench), "--outdir", str(out)],
        check=True,
        text=True,
        capture_output=True,
    )

    summary = out / "pangenome_summary.tsv"
    html = out / "mycosv_pangenome_calls_report.html"
    assert summary.exists()
    assert html.exists()
    text = summary.read_text(encoding="utf-8")
    assert "pangenome_only_calls\t1\t50.0" in text
    assert "hgt_or_starship_candidate_calls\t1\t50.0" in text
