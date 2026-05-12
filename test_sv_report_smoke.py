#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd

from analyze_phylo_sv_biology import analyze_mge_architecture, parse_vcf
from sv_visualization_report import benchmark_pred_call_counts, harmonize_columns


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


def test_harmonize_columns_preserves_original_alias_columns() -> None:
    df = pd.DataFrame({
        "query_asm": ["asm1"],
        "svtype": ["INS"],
        "element_class": ["LTR_GYPSY"],
    })

    out = harmonize_columns(df)

    assert out.loc[0, "query_asm"] == "asm1"
    assert out.loc[0, "sample"] == "asm1"
    assert out.loc[0, "svtype"] == "INS"
    assert out.loc[0, "sv_type"] == "INS"
    assert out.loc[0, "element_class"] == "LTR_GYPSY"
    assert out.loc[0, "te_class"] == "LTR_GYPSY"


def test_benchmark_summary_counts_deduplicate_truth_labels_and_drop_all() -> None:
    df = pd.DataFrame({
        "query_asm": ["asm1", "asm1", "asm1"],
        "method": ["mycosv", "mycosv", "mycosv"],
        "svtype": ["INS", "INS", "ALL"],
        "truth_label": [
            "consensus_2of_N_read_supported",
            "minigraph_read_supported",
            "consensus_2of_N_read_supported",
        ],
        "pred_calls": [7, 7, 99],
    })

    out = benchmark_pred_call_counts(df)

    assert out.to_dict("records") == [{
        "sample": "asm1",
        "caller": "mycosv",
        "sv_type": "INS",
        "count": 7,
    }]


def test_phylo_vcf_parser_accepts_mycosv_info_aliases(tmp_path: Path) -> None:
    vcf = tmp_path / "calls.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "ctgA\t100\t.\tN\t<OFF_REF>\t60\tPASS\tSVTYPE=OFF_REF;QUERY=asm1;EC=STARSHIP\n",
        encoding="utf-8",
    )

    rows = parse_vcf(vcf)

    assert rows[0]["query_asm"] == "asm1"
    assert rows[0]["element_class"] == "STARSHIP"


def test_mge_architecture_excludes_non_mge_category() -> None:
    rows = analyze_mge_architecture(
        sv_rows=[],
        bio_rows=[
            {"query_asm": "asm1", "element_class": "NONE", "svtype": "DEL", "phylum": "Ascomycota"},
            {"query_asm": "asm1", "element_class": "LTR_GYPSY", "svtype": "INS", "phylum": "Ascomycota"},
        ],
        asm_meta={},
        taxonomy={},
    )

    assert {row["mge_category"] for row in rows} == {"transposable"}
