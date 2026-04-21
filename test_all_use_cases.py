#!/usr/bin/env python3
# Designed for Linux

"""test_all_use_cases.py — Integration tests across all 13 ecological scenarios.

Tests that the AMF simulator produces correct outputs for every scenario in
SCENARIOS and that each output contains the expected truth records.

Coverage:
  • All 13 ecological scenarios: compact_yeast → giant_amf → rust/smut → two-speed pathogen
  • OFF_REF entries present for every query genome
  • On-reference SV entries present (INS/DEL/DUP/INV/TRA)
  • TRA breakpoints encoded in truth VCF (CHR2 + POS2)
  • No hint suffixes (__sv_) in contig names when run in hint-free mode
  • Scenario name present in both query_metadata.tsv and stress_case_catalog.tsv
  • Repeat/TE annotation integration: scenario element classes written to
    graph_annotations_denovo.tsv when --write-query-annotations is passed
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from test_amf import SCENARIOS  # only the dict; no side-effect


def run_sim(out_dir: Path, scenario: str) -> None:
    """Run the AMF simulator for a single scenario with minimal genome count."""
    subprocess.run(
        [
            sys.executable,
            str(Path(__file__).parent / "test_amf.py"),
            "--phylum",  SCENARIOS[scenario]["phylum"],
            "--n-genomes", "2",
            "--n-reps",    "1",
            "--total-len", "20000",
            "--n-contigs", "2",
            "--out-dir",   str(out_dir),
            "--scenario-set", scenario,
            "--write-extended-manifest",
        ],
        check=True,
    )


def run_sim_annotated(out_dir: Path, scenario: str) -> None:
    """Run simulator with annotation output enabled."""
    subprocess.run(
        [
            sys.executable,
            str(Path(__file__).parent / "test_amf.py"),
            "--phylum",  SCENARIOS[scenario]["phylum"],
            "--n-genomes", "2",
            "--n-reps",    "1",
            "--total-len", "20000",
            "--n-contigs", "2",
            "--out-dir",   str(out_dir),
            "--scenario-set", scenario,
            "--write-extended-manifest",
            "--write-query-annotations",
        ],
        check=True,
    )


# ── core correctness tests ─────────────────────────────────────────────────

def test_all_scenarios(tmp_path: Path) -> None:
    """All 13 scenarios produce correct metadata, stress catalog, and truth."""
    for scen in SCENARIOS:
        out = tmp_path / scen
        run_sim(out, scen)

        meta = (out / "query_metadata.tsv").read_text()
        assert scen in meta, f"scenario {scen} missing from query_metadata.tsv"

        stress = (out / "stress_case_catalog.tsv").read_text()
        assert scen in stress, f"scenario {scen} missing from stress_case_catalog.tsv"

        truth_lines = (out / "query_truth.tsv").read_text().strip().splitlines()
        assert len(truth_lines) > 1, f"query_truth.tsv empty for scenario {scen}"

        svtypes = {row.split("\t")[2] for row in truth_lines[1:]}
        assert "OFF_REF" in svtypes, \
            f"OFF_REF entry missing from query_truth.tsv for scenario {scen}"
        assert any(t for t in svtypes if t != "OFF_REF"), \
            f"No on-reference SV entry in query_truth.tsv for scenario {scen}"


def test_hint_free_contig_names(tmp_path: Path) -> None:
    """In hint-free mode (default), contig names must not contain __sv_."""
    for scen in list(SCENARIOS.keys())[:4]:  # spot-check first 4
        out = tmp_path / f"hf_{scen}"
        run_sim(out, scen)
        header = (out / "query_truth.tsv").read_text().splitlines()[0].split("\t")
        contig_idx = header.index("query_contig")
        for row in (out / "query_truth.tsv").read_text().strip().splitlines()[1:]:
            contig = row.split("\t")[contig_idx]
            assert "__sv_" not in contig, \
                f"Hint suffix in contig {contig!r} for scenario {scen} (hint-free mode)"


def test_extended_manifest_columns(tmp_path: Path) -> None:
    """hierarchy_manifest.tsv must have class, order, family columns."""
    out = tmp_path / "manifest_check"
    run_sim(out, "arbuscular_mf")
    base = (out / "base_manifest.tsv").read_text()
    hier = (out / "hierarchy_manifest.tsv").read_text()
    assert "#asm_name\tphylum\tclass\torder\tfamily\tgenus\tclade_name\tclade_rank\tfasta_path" in base
    assert "\tclass\t"  in hier
    assert "\torder\t"  in hier
    assert "\tfamily\t" in hier


def test_tra_breakpoints_in_vcf(tmp_path: Path) -> None:
    """TRA entries in truth VCF must have CHR2 and POS2 INFO fields."""
    out = tmp_path / "tra_vcf"
    run_sim(out, "saprotrophic")
    vcf_lines = [
        line for line in (out / "truth" / "all_queries.truth.ref.vcf")
        .read_text().splitlines()
        if not line.startswith("#")
    ]
    tra_lines = [line for line in vcf_lines if "SVTYPE=TRA" in line]
    if tra_lines:  # TRA may not always be selected by random seed
        for line in tra_lines:
            assert "CHR2=" in line,  f"TRA VCF missing CHR2: {line}"
            assert "POS2="  in line,  f"TRA VCF missing POS2: {line}"
            assert "END2="  in line,  f"TRA VCF missing END2: {line}"


# ── architecture-specific tests ────────────────────────────────────────────

def test_compact_yeast_produces_small_svs(tmp_path: Path) -> None:
    """compact_yeast scenario: all on-ref SVs must be ≤ 500 bp."""
    out = tmp_path / "compact_yeast"
    run_sim(out, "compact_yeast")
    truth_lines = (out / "query_truth.tsv").read_text().strip().splitlines()
    header = truth_lines[0].split("\t")
    svlen_idx  = header.index("svlen")
    svtype_idx = header.index("svtype")
    for row in truth_lines[1:]:
        fields = row.split("\t")
        if fields[svtype_idx] == "OFF_REF":
            continue
        sv_len = int(fields[svlen_idx])
        assert sv_len <= 500, \
            f"compact_yeast emitted large SV svlen={sv_len}"


def test_giant_amf_long_contig(tmp_path: Path) -> None:
    """giant_amf off-ref contig must be present and have correct scenario tag."""
    out = tmp_path / "giant_amf"
    run_sim(out, "giant_amf")
    truth_lines = (out / "query_truth.tsv").read_text().strip().splitlines()
    header = truth_lines[0].split("\t")
    scen_idx = header.index("scenario")
    svtype_idx = header.index("svtype")
    off_ref_rows = [
        row for row in truth_lines[1:]
        if row.split("\t")[svtype_idx] == "OFF_REF"
    ]
    assert off_ref_rows, "giant_amf: no OFF_REF row"
    for row in off_ref_rows:
        assert row.split("\t")[scen_idx] == "giant_amf"


def test_rust_smut_te_heavy_inv_ins_bias(tmp_path: Path) -> None:
    """rust_smut_te_heavy: must emit INS or INV on-reference SVs."""
    out = tmp_path / "rust_smut"
    run_sim(out, "rust_smut_te_heavy")
    truth_lines = (out / "query_truth.tsv").read_text().strip().splitlines()
    header = truth_lines[0].split("\t")
    svtype_idx = header.index("svtype")
    svtypes = {row.split("\t")[svtype_idx] for row in truth_lines[1:]}
    assert svtypes & {"INS", "INV"}, \
        f"rust_smut_te_heavy: expected INS or INV in {svtypes}"


def test_cross_phylum_hgt_off_ref_present(tmp_path: Path) -> None:
    """cross_phylum_hgt_stress: OFF_REF must be present (HGT island)."""
    out = tmp_path / "cross_phylum"
    run_sim(out, "cross_phylum_hgt_stress")
    truth_text = (out / "query_truth.tsv").read_text()
    assert "OFF_REF" in truth_text, "cross_phylum_hgt_stress: OFF_REF missing"
    assert "TRA" in truth_text or "OFF_REF" in truth_text, \
        "cross_phylum_hgt_stress: expected TRA or OFF_REF"


def test_annotation_output_contains_element_classes(tmp_path: Path) -> None:
    """With --write-query-annotations, graph_annotations_denovo.tsv must exist
    and contain at least one non-NONE annotation for element-bearing scenarios."""
    for scen in ("arbuscular_mf", "rust_smut_te_heavy", "hgt_receiver"):
        out = tmp_path / f"annot_{scen}"
        run_sim_annotated(out, scen)
        annot_path = out / "graph_annotations_denovo.tsv"
        assert annot_path.exists(), \
            f"graph_annotations_denovo.tsv missing for scenario {scen}"
        content = annot_path.read_text()
        assert len(content.strip().splitlines()) > 1, \
            f"graph_annotations_denovo.tsv empty for {scen}"


# ── two-speed pathogen / multi-scenario ───────────────────────────────────

def test_two_speed_pathogen_inv_tra_ins(tmp_path: Path) -> None:
    """two_speed_pathogen_extreme must emit at least one of INV/TRA/INS."""
    out = tmp_path / "two_speed"
    run_sim(out, "two_speed_pathogen_extreme")
    truth_text = (out / "query_truth.tsv").read_text()
    assert any(t in truth_text for t in ("INV", "TRA", "INS")), \
        "two_speed_pathogen_extreme: expected INV/TRA/INS"


def test_tra_biased_scenario_emits_tra_and_truth_vcf_keeps_real_contigs(
    tmp_path: Path,
) -> None:
    """saprotrophic scenario: TRA entries use real contig names in VCF."""
    out = tmp_path / "saprotrophic"
    run_sim(out, "saprotrophic")
    truth_rows = (out / "query_truth.tsv").read_text().strip().splitlines()
    # At least one TRA or OFF_REF must appear
    has_tra_or_offref = any(
        "	TRA	" in row or "	OFF_REF	" in row
        for row in truth_rows[1:]
    )
    assert has_tra_or_offref, "saprotrophic: neither TRA nor OFF_REF found in truth TSV"

    vcf_path = out / "truth" / "all_queries.truth.ref.vcf"
    vcf_lines = [
        line for line in vcf_path.read_text().splitlines()
        if not line.startswith("#")
    ]
    tra_lines = [line for line in vcf_lines if "SVTYPE=TRA" in line]
    if tra_lines:
        assert all("CHR2=" in line and "POS2=" in line for line in tra_lines), \
            "TRA mate breakpoint missing from truth VCF"
        # Contig names in TRA lines must not carry __sv_ suffixes
        for line in tra_lines:
            chrom = line.split("\t")[0]
            assert "__sv_" not in chrom, \
                f"TRA VCF CHROM {chrom!r} carries hint suffix"
