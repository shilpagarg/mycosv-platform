#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "plot_mycosv_pangenome_calls.py"


def _load_plot_module():
    spec = importlib.util.spec_from_file_location("plot_mycosv_pangenome_calls", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


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


def test_evaluate_panel_claims_gates_both_manuscript_claims(tmp_path: Path) -> None:
    plot = _load_plot_module()
    # rv_rows mimics aggregate_panel_read_validation output. Three queries:
    #   q_high — 90 % MycoSV read-validation rate -> in 78-100 band
    #   q_low  — 50 % rate                         -> below band
    #   q_top  — 100 % rate                        -> in band
    rv_rows = [
        {"query_asm": "q_high", "source": "mycosv",
         "yes": 90, "total": 100, "rate": 0.90,
         "ci95_lo": 0.82, "ci95_hi": 0.95},
        {"query_asm": "q_low", "source": "mycosv",
         "yes": 50, "total": 100, "rate": 0.50,
         "ci95_lo": 0.40, "ci95_hi": 0.60},
        {"query_asm": "q_top", "source": "mycosv",
         "yes": 40, "total": 40, "rate": 1.00,
         "ci95_lo": 0.91, "ci95_hi": 1.00},
        # An unrelated comparator source must be ignored by evaluate_panel_claims.
        {"query_asm": "q_high", "source": "minigraph",
         "yes": 5, "total": 50, "rate": 0.10,
         "ci95_lo": 0.04, "ci95_hi": 0.21},
    ]
    pl_rows = [
        # 75 % pangenome-only-read-supported (in 60-90 band)
        {"query_asm": "q_high", "dedup_loci": 100,
         "single_ref_equivalent": 25, "pangenome_only": 75,
         "pangenome_only_read_supported": 75},
        # 30 % pangenome-only-read-supported (below band)
        {"query_asm": "q_low", "dedup_loci": 100,
         "single_ref_equivalent": 70, "pangenome_only": 30,
         "pangenome_only_read_supported": 30},
        # 95 % pangenome-only-read-supported (above band)
        {"query_asm": "q_top", "dedup_loci": 40,
         "single_ref_equivalent": 2, "pangenome_only": 38,
         "pangenome_only_read_supported": 38},
    ]
    out_tsv = tmp_path / "panel_claim_validation.tsv"
    out_md = tmp_path / "panel_claim_validation.md"
    rows = plot.evaluate_panel_claims(
        rv_rows, pl_rows, out_tsv, out_md, panel_label="5sample",
    )
    by_q = {r["query_asm"]: r for r in rows}
    assert by_q["q_high"]["claim1_read_validation_status"] == "in_range"
    assert by_q["q_low"]["claim1_read_validation_status"] == "below_range"
    assert by_q["q_top"]["claim1_read_validation_status"] == "in_range"
    assert by_q["q_high"]["claim2_pangenome_only_status"] == "in_range"
    assert by_q["q_low"]["claim2_pangenome_only_status"] == "below_range"
    assert by_q["q_top"]["claim2_pangenome_only_status"] == "above_range"
    # Pooled ALL row: 90+50+40=180 / 240=0.75; PO_read 75+30+38=143 / 240=~0.596
    all_row = by_q["ALL"]
    assert all_row["mycosv_calls_total"] == 240
    assert all_row["mycosv_calls_read_validated"] == 180
    assert all_row["pangenome_only_read_supported_loci"] == 143
    # In-range counters denominator is the row count (3 queries + ALL).
    assert "in_range_queries=2/" in all_row["claim1_read_validation_status"]
    assert "in_range_queries=1/" in all_row["claim2_pangenome_only_status"]
    assert out_tsv.exists()
    assert out_md.exists()
    md_text = out_md.read_text(encoding="utf-8")
    assert "Claim 1" in md_text and "Claim 2" in md_text


def test_aggregate_panel_read_validation_unions_bam_and_refree(tmp_path: Path) -> None:
    plot = _load_plot_module()
    by_query = tmp_path / "by_query"
    shard = by_query / "q1"
    shard.mkdir(parents=True)
    # Same call (ctgA:100-200 DEL) marked NO in the BAM-anchored file and
    # YES in the reference-free file: the union must count it as validated.
    # A second call (ctgB:50-50 INS) is YES only in BAM. A third
    # (ctgC:300-400 INV) is YES only in refree. Total = 3 calls, all
    # validated under the union -> rate = 1.0.
    header = (
        "query_asm\tref_contig\tpos\tend\tsvtype\tsource\tcoord_space\t"
        "read_support\tvalidation_support\tsupport_source\tread_validated\tstatus\n"
    )
    (shard / "read_validated_truth.tsv").write_text(
        header
        + "q1\tctgA\t100\t200\tDEL\tmycosv\tquery\t-1\t0\texternal_validation\tno\tnot_validated\n"
        + "q1\tctgB\t50\t50\tINS\tmycosv\tquery\t-1\t5\texternal_validation\tyes\tvalidated\n",
        encoding="utf-8",
    )
    (shard / "read_validated_truth_refree.tsv").write_text(
        header
        + "q1\tctgA\t100\t200\tDEL\tmycosv\tquery\t-1\t6\treference_free_junction_kmer_alignment\tyes\tvalidated\n"
        + "q1\tctgC\t300\t400\tINV\tmycosv\tquery\t-1\t4\treference_free_junction_kmer_alignment\tyes\tvalidated\n",
        encoding="utf-8",
    )
    out = tmp_path / "panel_read_validation_rate.tsv"
    rows = plot.aggregate_panel_read_validation(by_query, out)
    mycosv = [r for r in rows if r["source"] == "mycosv"]
    assert len(mycosv) == 1
    rec = mycosv[0]
    assert rec["total"] == 3
    assert rec["yes"] == 3
    assert rec["rate"] == 1.0
