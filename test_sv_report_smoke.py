#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

pd = pytest.importorskip("pandas")

from analyze_phylo_sv_biology import analyze_mge_architecture, parse_vcf
from sv_visualization_report import (
    _select_wins_subset,
    biology_effect_plot_table,
    benchmark_pred_call_counts,
    build_biology_section,
    harmonize_columns,
    numeric_pair_table,
    plot_wins_matrix,
    top_n,
)


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
        "query_asm": ["asm1", "asm1", "asm1", "asm1"],
        "method": ["mycosv", "mycosv", "mycosv", "mycosv"],
        "svtype": ["INS", "INS", "ALL", "INS"],
        "truth_label": [
            "consensus_2of_N_read_supported",
            "minigraph_read_supported",
            "consensus_2of_N_read_supported",
            "consensus_2of_N_read_supported",
        ],
        "coordinate_space": ["query", "query", "query", "reference_any_clade"],
        "pred_calls": [7, 7, 99, 123],
    })

    out = benchmark_pred_call_counts(df)

    assert out.to_dict("records") == [{
        "sample": "asm1",
        "caller": "mycosv",
        "sv_type": "INS",
        "count": 7,
    }]


def test_wins_subset_keeps_consensus_rows_and_excludes_single_tool_truth() -> None:
    df = pd.DataFrame({
        "query_asm": ["asm1", "asm1", "asm1"],
        "truth_label": [
            "consensus_2of_2",
            "consensus_2of_2_read_supported",
            "minigraph_read_supported",
        ],
        "method": ["mycosv", "minigraph", "mycosv"],
        "svtype": ["INS", "INS", "INS"],
        "truth_calls": [3, 2, 5],
        "f1": [0.25, 0.75, 0.90],
        "status": ["ok", "ok", "ok"],
    })

    out = _select_wins_subset(df)

    assert out["truth_label"].tolist() == ["consensus_2of_2_read_supported"]
    assert out["f1"].tolist() == [0.75]


def test_wins_rate_counts_zero_ties_as_win_ties(tmp_path: Path) -> None:
    df = pd.DataFrame({
        "query_asm": ["asm1", "asm1"],
        "truth_label": ["consensus_2of_2_read_supported"] * 2,
        "method": ["mycosv", "sniffles"],
        "svtype": ["INV", "INV"],
        "truth_calls": [1, 1],
        "f1": [0.0, 0.0],
        "status": ["ok", "ok"],
    })

    plot_wins_matrix(df, tmp_path)
    wins = pd.read_csv(tmp_path / "wins_rate.tsv", sep="\t")

    assert wins.loc[0, "comparator"] == "sniffles"
    assert wins.loc[0, "mycosv_wins"] == 1
    assert wins.loc[0, "mycosv_ties"] == 1
    assert wins.loc[0, "zero_ties"] == 1
    assert wins.loc[0, "win_rate"] == 1.0


def test_top_n_excludes_null_like_gene_labels() -> None:
    df = pd.DataFrame({"gene": ["NA", ".", None, "", "VPS27", "VPS27", "THI12"]})

    out = top_n(df, "gene", 20)

    assert out.to_dict("records") == [
        {"gene": "VPS27", "count": 2},
        {"gene": "THI12", "count": 1},
    ]


def test_biology_effect_plot_table_uses_informative_locus_effects() -> None:
    df = pd.DataFrame({
        "gene": ["VPS27", ".", "", "HXT4"],
        "affected_locus": ["loc_none", "loc_line", "loc_rip", "loc_hgt"],
        "query_contig": ["ctg1", "ctg2", "ctg3", "ctg4"],
        "effect": ["NONE", "TE_LINE", "RIP", "HGT"],
    })

    out = biology_effect_plot_table(df, top_n_labels=1)

    assert out[["_effect_label", "_effect_class"]].to_dict("records") == [
        {"_effect_label": "loc_line", "_effect_class": "TE_LINE"},
        {"_effect_label": "loc_rip", "_effect_class": "RIP"},
        {"_effect_label": "HXT4", "_effect_class": "HGT"},
    ]


def test_numeric_pair_table_drops_non_numeric_expression_pairs() -> None:
    df = pd.DataFrame({
        "cn": ["1", "2", "bad"],
        "expression": [".", "1.5", "3.0"],
    })

    out = numeric_pair_table(df, "cn", "expression")

    assert out.to_dict("records") == [{"cn": 2.0, "expression": 1.5}]


def test_biology_section_skips_cn_expression_when_no_numeric_pairs(tmp_path: Path) -> None:
    df = pd.DataFrame({
        "gene": ["VPS27"],
        "effect": ["HGT"],
        "cn": ["6"],
        "expression": ["."],
    })

    html, figs, _ = build_biology_section(df, tmp_path)

    assert "Copy number vs expression plot skipped" in html
    assert "biology_cn_vs_expression.png" not in [fig.filename for fig in figs]
    assert (tmp_path / "biology_cn_vs_expression.tsv").exists()


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
