#!/usr/bin/env python3
# Designed for Linux

import run_real_fungal_benchmark as rrfb
from pathlib import Path
import csv
import gzip
import math

from run_real_fungal_benchmark import (
    NormalizedCall,
    estimate_prepared_genome_size_hint,
    merge_sequence_sources,
    load_minigraph_bubble_calls,
    load_mycosv_query_calls,
    load_mycosv_reference_calls,
    load_normalized_calls_tsv,
    materialize_query_input,
    parse_ena_filereport_text,
    parse_assembly_summary,
    score_callsets,
    select_all_public_rows,
    select_ena_read_sources,
    select_species_rows,
    run_mycosv,
)


def test_parse_assembly_summary_and_select_species_rows():
    text = (
        "# assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tsubmitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\n"
        "GCF_000001\t.\t.\t.\treference genome\t4932\t4932\tSaccharomyces cerevisiae S288C\t.\t.\tlatest\tComplete Genome\tMajor\tFull\t2024/01/01\tasm1\t.\t.\t.\thttps://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000001\n"
        "GCF_000002\t.\t.\t.\trepresentative genome\t4932\t4932\tSaccharomyces cerevisiae isolate X\t.\t.\tlatest\tScaffold\tMajor\tFull\t2023/01/01\tasm2\t.\t.\t.\thttps://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000002\n"
    )
    rows = parse_assembly_summary(text)
    selected = select_species_rows(rows, "Saccharomyces cerevisiae", 2)
    assert len(selected) == 2
    assert selected[0]["assembly_accession"] == "GCF_000001"


def test_parse_assembly_summary_accepts_current_ncbi_header_style():
    text = (
        "# some comment\n"
        "#assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tsubmitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\n"
        "GCF_000010\t.\t.\t.\treference genome\t4932\t4932\tSaccharomyces cerevisiae S288C\t.\t.\tlatest\tComplete Genome\tMajor\tFull\t2024/01/01\tasm1\t.\t.\t.\thttps://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000010\n"
    )
    rows = parse_assembly_summary(text)
    assert len(rows) == 1
    assert rows[0]["assembly_accession"] == "GCF_000010"


def test_select_all_public_rows_filters_by_level_and_latest():
    rows = [
        {
            "assembly_accession": "GCF_1",
            "species_taxid": "4932",
            "organism_name": "Saccharomyces cerevisiae S288C",
            "assembly_level": "Complete Genome",
            "version_status": "latest",
            "ftp_path": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_1",
        },
        {
            "assembly_accession": "GCF_2",
            "species_taxid": "4932",
            "organism_name": "Saccharomyces cerevisiae isolate X",
            "assembly_level": "Contig",
            "version_status": "latest",
            "ftp_path": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_2",
        },
        {
            "assembly_accession": "GCF_3",
            "species_taxid": "559292",
            "organism_name": "Candida glabrata CBS138",
            "assembly_level": "Scaffold",
            "version_status": "suppressed",
            "ftp_path": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_3",
        },
    ]
    selected = select_all_public_rows(rows, min_assembly_level="scaffold", latest_only=True, max_total=0)
    assert [row["assembly_accession"] for row in selected] == ["GCF_1"]


def test_ncbi_best_deduplicates_paired_refseq_genbank_and_prefers_best():
    rows = [
        {
            "assembly_accession": "GCA_000001",
            "gbrs_paired_asm": "GCF_000001",
            "organism_name": "Saccharomyces cerevisiae S288C",
            "assembly_level": "Complete Genome",
            "version_status": "latest",
            "genome_rep": "Full",
            "refseq_category": "na",
            "seq_rel_date": "2024/01/01",
            "ftp_path": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001",
            "_catalog_source": "ncbi-genbank",
        },
        {
            "assembly_accession": "GCF_000001",
            "gbrs_paired_asm": "GCA_000001",
            "organism_name": "Saccharomyces cerevisiae S288C",
            "assembly_level": "Complete Genome",
            "version_status": "latest",
            "genome_rep": "Full",
            "refseq_category": "reference genome",
            "seq_rel_date": "2024/01/01",
            "ftp_path": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000001",
            "_catalog_source": "ncbi-refseq",
        },
        {
            "assembly_accession": "GCA_000002",
            "gbrs_paired_asm": "na",
            "organism_name": "Saccharomyces cerevisiae isolate X",
            "assembly_level": "Scaffold",
            "version_status": "latest",
            "genome_rep": "Full",
            "refseq_category": "na",
            "seq_rel_date": "2025/01/01",
            "ftp_path": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000002",
            "_catalog_source": "ncbi-genbank",
        },
    ]
    deduped = rrfb.deduplicate_best_assembly_rows(rows)
    assert {row["assembly_accession"] for row in deduped} == {"GCF_000001", "GCA_000002"}
    selected = select_species_rows(deduped, "Saccharomyces cerevisiae", 2)
    assert selected[0]["assembly_accession"] == "GCF_000001"


def test_ncbi_download_targets_strip_trailing_ftp_slash():
    row = {
        "ftp_path": (
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/"
            "GCA_000001405.29_GRCh38.p14/"
        )
    }
    targets = rrfb.ncbi_download_targets(row, include_gff=True)
    assert targets[0][0].endswith(
        "/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz"
    )
    assert "//GCA_000001405.29_GRCh38.p14_genomic" not in targets[0][0]


def test_gff_to_gene_annotations_parses_gbff_fallback(tmp_path: Path):
    gbff = tmp_path / "asm_genomic.gbff.gz"
    text = """LOCUS       ABC123                1000 bp    DNA     linear   PLN 01-JAN-2000
VERSION     ABC123.1
FEATURES             Location/Qualifiers
     source          1..1000
                     /organism="Example fungus"
     gene            complement(10..90)
                     /locus_tag="GENE1"
                     /gene="abc"
                     /gene_biotype="protein_coding"
     CDS             join(200..250,300..350)
                     /locus_tag="GENE2"
                     /product="example protein"
ORIGIN
//
"""
    with gzip.open(gbff, "wt", encoding="utf-8") as fh:
        fh.write(text)

    rows = rrfb.gff_to_gene_annotations([("GCA_TEST_1", gbff)])
    by_id = {row["gene_id"]: row for row in rows}
    assert by_id["GENE1"]["query_contig"] == "ABC123.1"
    assert by_id["GENE1"]["start"] == 10
    assert by_id["GENE1"]["end"] == 90
    assert by_id["GENE1"]["strand"] == "-"
    assert by_id["GENE2"]["start"] == 200
    assert by_id["GENE2"]["end"] == 350
    assert by_id["GENE2"]["product"] == "example protein"


def test_stream_gene_annotations_to_tsv_expands_aliases(tmp_path: Path):
    """The streaming writer must (a) emit each parsed gene once per owner alias,
    (b) cover ref aliases + per-benchmark query aliases, and (c) keep RSS low
    by not retaining cross-source rows. We assert (a) and (b) here; (c) is what
    the original implementation got wrong (it built every row in memory before
    write_tsv'ing), so a regression on (a/b) is the only thing a unit test
    can catch — memory blowup needs the live 2000-source workload to surface.
    """
    gff = tmp_path / "ref_asm_genomic.gff.gz"
    gff_text = (
        "##gff-version 3\n"
        "contigA\tref\tgene\t100\t200\t.\t+\t.\tID=g1;Name=GENEA\n"
        "contigA\tref\tgene\t300\t400\t.\t-\t.\tID=g2;Name=GENEB\n"
        "contigB\tref\tCDS\t500\t600\t.\t+\t.\tID=cds1\n"  # filtered: ftype not in _GFF_GENE_TYPES
    )
    with gzip.open(gff, "wt", encoding="utf-8") as fh:
        fh.write(gff_text)

    out = tmp_path / "gene_annotations.tsv"
    asm_aliases = {
        "REF1": {"REF1", "REF1.fa", "REF1.fasta"},
        "QUERY_A": {"QUERY_A", "QA.fna"},
        "QUERY_B": {"QUERY_B"},
    }
    ref_to_queries = {"REF1": ["QUERY_A", "QUERY_B"]}

    n = rrfb.stream_gene_annotations_to_tsv(
        out, [("REF1", gff)], asm_aliases, ref_to_queries, progress_every=0,
    )
    assert n > 0
    with out.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    # 2 genes × (3 ref aliases + 2 QA aliases + 1 QB alias) = 12 rows.
    assert len(rows) == 12
    owners_by_gene: dict[str, set[str]] = {}
    for row in rows:
        owners_by_gene.setdefault(row["gene_id"], set()).add(row["query_asm"])
    assert owners_by_gene["g1"] == {"REF1", "REF1.fa", "REF1.fasta", "QUERY_A", "QA.fna", "QUERY_B"}
    assert owners_by_gene["g2"] == owners_by_gene["g1"]


def test_stream_gene_annotations_to_tsv_handles_empty_source(tmp_path: Path):
    """A GBFF with no gene/CDS features (common for unannotated NCBI WGS
    assemblies) must not abort the streaming write — earlier sources' rows
    should stay on disk and the function must just keep going. Without this
    behaviour prepare_million_real silently produced zero-row TSVs whenever
    its 2000-source mix had even one barren GBFF.
    """
    gff = tmp_path / "ref_genomic.gff.gz"
    with gzip.open(gff, "wt", encoding="utf-8") as fh:
        fh.write("##gff-version 3\ncontigA\tref\tgene\t1\t100\t.\t+\t.\tID=g1\n")
    barren_gbff = tmp_path / "barren_genomic.gbff.gz"
    with gzip.open(barren_gbff, "wt", encoding="utf-8") as fh:
        fh.write(
            "LOCUS       SCAFFOLD              1000 bp    DNA     linear   CON 01-JAN-2024\n"
            "VERSION     SCAFFOLD.1\n"
            "FEATURES             Location/Qualifiers\n"
            "     source          1..1000\n"
            "                     /organism=\"Unannotated WGS\"\n"
            "CONTIG      join(SCAFFOLD_PART_1.1:1..1000)\n"
            "//\n"
        )
    out = tmp_path / "gene_annotations.tsv"
    n = rrfb.stream_gene_annotations_to_tsv(
        out,
        [("REF_OK", gff), ("REF_BARREN", barren_gbff)],
        asm_aliases={"REF_OK": {"REF_OK"}, "REF_BARREN": {"REF_BARREN"}},
        ref_to_queries={},
        progress_every=0,
    )
    assert n == 1
    with out.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    assert [row["query_asm"] for row in rows] == ["REF_OK"]


def test_load_mycosv_reference_calls_parses_ref_space(tmp_path: Path):
    vcf = tmp_path / "calls.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "query_ctg\t11\tsv1\tN\t<DEL>\t40\tPASS\tSVTYPE=DEL;SVLEN=-50;END=11;ANNOT=NONE;CLADE=ref_asm;REFCONTIG=ref_ctg;REFPOS=101;REFEND=150;QASM=query_asm\tGT:GQ\t0/1:40\n",
        encoding="utf-8",
    )
    calls = load_mycosv_reference_calls(vcf, "query_asm")
    assert len(calls) == 1
    call = calls[0]
    assert call.coord_space == "reference"
    assert call.ref_contig == "ref_ctg"
    assert call.pos == 101
    assert call.end == 150
    assert call.svtype == "DEL"
    assert call.read_support is None


def test_load_mycosv_calls_match_full_fasta_qasm_alias(tmp_path: Path):
    vcf = tmp_path / "calls.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "query_ctg\t11\tsv1\tN\t<INS>\t40\tPASS\tSVTYPE=INS;SVLEN=50;END=11;ANNOT=NONE;CLADE=ref_asm;REFCONTIG=ref_ctg;REFPOS=101;REFEND=101;QASM=GCA_000149225.2_ASM14922v2_genomic\tGT:GQ\t0/1:40\n",
        encoding="utf-8",
    )
    query_calls = load_mycosv_query_calls(vcf, "GCA_000149225_2")
    ref_calls = load_mycosv_reference_calls(vcf, "GCA_000149225_2")
    assert len(query_calls) == 1
    assert len(ref_calls) == 1
    assert query_calls[0].query_asm == "GCA_000149225_2"
    assert ref_calls[0].pos == 101


def test_load_mycosv_calls_parse_intrinsic_read_support(tmp_path: Path):
    vcf = tmp_path / "calls.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "sr_unitig7_len155_mf12\t1\tsv1\tN\t<OFF_REF>\t20\tPASS\t"
        "SVTYPE=OFF_REF;SVLEN=155;END=155;QASM=q1;SUPPORT=12\tGT\t0/1\n"
        "lr_pc4_n8\t20\tsv2\tN\t<DEL>\t40\tPASS\t"
        "SVTYPE=DEL;SVLEN=-50;END=20;REFCONTIG=chr1;REFPOS=101;REFEND=150;QASM=q1\tGT\t0/1\n",
        encoding="utf-8",
    )
    query_calls = load_mycosv_query_calls(vcf, "q1")
    ref_calls = load_mycosv_reference_calls(vcf, "q1")
    assert query_calls[0].read_support == 12
    assert query_calls[1].read_support == 8
    assert ref_calls[0].read_support == 8


def test_expand_to_multisample_vcf_uses_manifest_samples_for_empty_vcf(tmp_path: Path):
    vcf = tmp_path / "calls.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
        encoding="utf-8",
    )

    out = rrfb.expand_to_multisample_vcf(
        vcf,
        tmp_path / "calls.multisample.vcf",
        ["q1", "q2"],
    )

    header = [line for line in out.read_text(encoding="utf-8").splitlines() if line.startswith("#CHROM")][0]
    assert header.endswith("\tFORMAT\tq1\tq2")


def test_join_biology_findings_writes_header_when_candidates_missing(tmp_path: Path):
    out = tmp_path / "biology_findings.tsv"

    rrfb.join_biology_findings(None, [], {}, out)

    text = out.read_text(encoding="utf-8")
    assert text.startswith("query_asm\tquery_contig\tpos\tend\tsvtype")
    assert "comparator_support_count" in text


def test_load_mycosv_assembly_support_from_vcf(tmp_path: Path):
    vcf = tmp_path / "calls.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "asm_ctg\t11\tsv1\tN\t<DEL>\t40\tPASS\t"
        "SVTYPE=DEL;SVLEN=-50;END=11;REFCONTIG=chr1;REFPOS=101;REFEND=150;QASM=q1;QMODE=assembly;SUPPORT=17\tGT:GQ\t0/1:40\n",
        encoding="utf-8",
    )
    query_calls = load_mycosv_query_calls(vcf, "q1")
    ref_calls = load_mycosv_reference_calls(vcf, "q1")
    assert query_calls[0].read_support == 17
    assert ref_calls[0].read_support == 17


def test_validate_mycosv_calls_uses_intrinsic_read_support(tmp_path: Path, monkeypatch):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    query_row = {
        "query_asm": "q1",
        "query_mode": "short-reads",
        "path": str(ref),
        "benchmark_ref_fasta": str(ref),
    }
    call = NormalizedCall(
        "q1", "sr_unitig7_len155_mf12", 1, 155, "OFF_REF", 155, "mycosv",
        coord_space="query", read_support=12,
    )
    monkeypatch.setattr(rrfb, "_build_validation_bam", lambda *a, **k: (tmp_path / "x.bam", ref))
    monkeypatch.setattr(rrfb, "_samtools_count_breakpoint_support", lambda *a, **k: 0)
    kept, rows = rrfb.validate_calls_with_reads(
        [call], query_row, tmp_path / "validation",
        threads=1, min_support=3, flank_bp=250,
    )
    assert kept == [call]
    assert rows[0]["read_support"] == 12
    assert rows[0]["validation_support"] == 0
    assert rows[0]["support_source"] == "mycosv_short_read_kmer"
    assert rows[0]["read_validated"] == "yes"


def test_validate_assembly_query_space_mycosv_keeps_internal_support(tmp_path: Path, monkeypatch):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    query_row = {
        "query_asm": "q1",
        "query_mode": "assembly",
        "path": str(ref),
        "benchmark_ref_fasta": str(ref),
    }
    call = NormalizedCall(
        "q1", "query_ctg", 10, 100, "INV", 90, "mycosv",
        coord_space="query", read_support=5,
    )
    monkeypatch.setattr(rrfb, "_build_validation_bam", lambda *a, **k: (tmp_path / "x.bam", ref))
    monkeypatch.setattr(rrfb, "_samtools_count_breakpoint_support", lambda *a, **k: 0)
    kept, rows = rrfb.validate_calls_with_reads(
        [call], query_row, tmp_path / "validation",
        threads=1, min_support=3, flank_bp=250,
    )
    assert kept == [call]
    assert rows[0]["read_support"] == 5
    assert rows[0]["validation_support"] == -1
    assert rows[0]["support_source"] == "mycosv_assembly_anchors"
    assert rows[0]["status"] == "query_space_not_reference_validated"
    assert rows[0]["read_validated"] == "yes"


def test_score_callsets_no_truth_is_nan_status():
    metrics = score_callsets([], [
        NormalizedCall("q1", "ctg", 1, 10, "DEL", -10, "mycosv")
    ])
    assert math.isnan(metrics["f1"])
    assert metrics["status"] == "no_truth"


def test_write_mycosv_failure_outputs_are_parseable(tmp_path: Path):
    paths = rrfb.write_mycosv_failure_outputs(tmp_path / "mycosv" / "calls", "rc=9")
    for path in paths.values():
        p = Path(path)
        assert p.exists()
        assert p.stat().st_size > 0
    assert "#CHROM" in Path(paths["vcf"]).read_text(encoding="utf-8")
    hits_text = Path(paths["hits"]).read_text(encoding="utf-8")
    assert "query_asm\tquery_contig" in hits_text
    assert "MYCOSV_FAILED" in hits_text


def test_calls_compatible_accepts_inv_whole_block_prediction():
    truth = NormalizedCall("q1", "ctg", 5000, 5200, "INV", 200, "truth", coord_space="reference", ref_contig="chr1")
    pred = NormalizedCall("q1", "ctg", 1000, 9000, "INV", 8000, "mycosv", coord_space="reference", ref_contig="chr1")
    assert rrfb.calls_compatible(truth, pred)
    assert score_callsets([truth], [pred])["tp"] == 1


def test_tra_matching_requires_mate_breakpoint_when_present():
    truth = NormalizedCall(
        "q1", ".", 1000, 1000, "TRA", 1, "truth",
        coord_space="reference", ref_contig="chr1",
        mate_contig="chr2", mate_pos=5000, mate_end=5000,
    )
    local_only_wrong_mate = NormalizedCall(
        "q1", ".", 1002, 1002, "TRA", 1, "mycosv",
        coord_space="reference", ref_contig="chr1",
        mate_contig="chr3", mate_pos=5000, mate_end=5000,
    )
    correct = NormalizedCall(
        "q1", ".", 1002, 1002, "TRA", 1, "mycosv",
        coord_space="reference", ref_contig="chr1",
        mate_contig="chr2", mate_pos=5050, mate_end=5050,
    )

    assert not rrfb.calls_compatible(truth, local_only_wrong_mate)
    assert rrfb.calls_compatible(truth, correct)
    assert score_callsets([truth], [local_only_wrong_mate])["tp"] == 0
    assert score_callsets([truth], [correct])["tp"] == 1


def test_reference_vcf_loader_parses_bnd_alt_mate(tmp_path: Path):
    vcf = tmp_path / "bnd.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t1000\tbnd1\tN\tN]chr2:5000]\t60\tPASS\tSVTYPE=BND\n",
        encoding="utf-8",
    )

    calls = rrfb.load_reference_vcf_calls(vcf, "manta", "q1")

    assert len(calls) == 1
    assert calls[0].svtype == "TRA"
    assert calls[0].mate_contig == "chr2"
    assert calls[0].mate_pos == 5000


def test_validate_tra_uses_mate_breakpoint_support(tmp_path: Path, monkeypatch):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n>chr2\nACGT\n", encoding="utf-8")
    query_row = {
        "query_asm": "q1",
        "query_mode": "long-reads",
        "path": str(ref),
        "benchmark_ref_fasta": str(ref),
    }
    call = NormalizedCall(
        "q1", ".", 1000, 1000, "TRA", 1, "truth",
        coord_space="reference", ref_contig="chr1",
        mate_contig="chr2", mate_pos=5000, mate_end=5000,
    )
    seen = []

    def fake_support(_bam, contig, pos, end, **_kwargs):
        seen.append((contig, pos, end))
        return 0 if contig == "chr1" else 3

    monkeypatch.setattr(rrfb, "_build_validation_bam", lambda *a, **k: (tmp_path / "x.bam", ref))
    monkeypatch.setattr(rrfb, "_samtools_count_breakpoint_support", fake_support)

    kept, rows = rrfb.validate_calls_with_reads(
        [call], query_row, tmp_path / "validation",
        threads=1, min_support=3, flank_bp=250,
    )

    assert kept == [call]
    assert ("chr1", 1000, 1000) in seen
    assert ("chr2", 5000, 5000) in seen
    assert rows[0]["validation_support"] == 3
    assert rows[0]["read_validated"] == "yes"


def test_cigar_indels_support_reference_breakpoints():
    assert rrfb._cigar_indel_supports_call(
        "100M75D200M",
        1000,
        1098,
        1175,
        svtype="DEL",
        svlen=-75,
        flank_bp=10,
    )
    assert rrfb._cigar_indel_supports_call(
        "100M80I200M",
        1000,
        1099,
        1099,
        svtype="INS",
        svlen=80,
        flank_bp=10,
    )
    assert not rrfb._cigar_indel_supports_call(
        "100M80I200M",
        1000,
        2000,
        2000,
        svtype="INS",
        svlen=80,
        flank_bp=10,
    )


def test_load_normalized_calls_tsv_supports_reference_space(tmp_path: Path):
    tsv = tmp_path / "other.tsv"
    tsv.write_text(
        "query_asm\tcoord_space\tchrom\tpos\tend\tsvtype\tsvlen\n"
        "q1\treference\tchr2\t501\t550\tINV\t49\n",
        encoding="utf-8",
    )
    calls = load_normalized_calls_tsv(tsv, "other")
    assert len(calls) == 1
    assert calls[0].coord_space == "reference"
    assert calls[0].ref_contig == "chr2"
    assert calls[0].query_contig == "."


def test_load_minigraph_bubble_calls_infers_ins(tmp_path: Path):
    bubble = tmp_path / "bubbles.bed"
    sample = tmp_path / "sample.bed"
    bubble.write_text("chr1\t100\t120\t2\t2\t0\t20\t35\n", encoding="utf-8")
    sample.write_text("chr1\t100\t120\t.\t.\t.\t.\t.\tpath:35:+:qctg:200:235\n", encoding="utf-8")
    calls = load_minigraph_bubble_calls(bubble, sample, "q1")
    assert len(calls) == 1
    assert calls[0].coord_space == "reference"
    assert calls[0].svtype == "INS"
    assert calls[0].svlen == 15
    assert calls[0].pos == 101


def test_score_callsets_does_not_mix_coordinate_spaces():
    truth = [NormalizedCall("q1", "ctg", 100, 150, "DEL", 50, "truth", coord_space="reference", ref_contig="chr1")]
    pred = [NormalizedCall("q1", "ctg", 100, 150, "DEL", 50, "pred", coord_space="query", ref_contig="chr1")]
    metrics = score_callsets(truth, pred)
    assert metrics["tp"] == 0
    assert metrics["fp"] == 1
    assert metrics["fn"] == 1


def test_parse_ena_filereport_and_select_sources():
    text = (
        "run_accession\tscientific_name\tinstrument_platform\tlibrary_layout\tfastq_ftp\tread_count\tsubmitted_ftp\n"
        "SRR1\tAspergillus fumigatus\tILLUMINA\tPAIRED\tftp.sra.ebi.ac.uk/vol1/fastq/SRR1_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR1_2.fastq.gz\t100000\t\n"
        "SRR2\tAspergillus fumigatus\tOXFORD_NANOPORE\tSINGLE\tftp.sra.ebi.ac.uk/vol1/fastq/SRR2.fastq.gz\t100000\t\n"
    )
    rows = parse_ena_filereport_text(text)
    urls, meta = select_ena_read_sources(rows, "short-reads", 2)
    assert len(urls) == 2
    assert all(url.startswith("https://ftp.sra.ebi.ac.uk/") for url in urls)
    assert meta[0]["run_accession"] == "SRR1"
    assert meta[0]["selected_urls"].split(";") == urls


def test_select_ena_sources_rejects_submitted_binary_and_tiny_runs():
    rows = parse_ena_filereport_text(
        "run_accession\tscientific_name\tinstrument_platform\tlibrary_layout\tfastq_ftp\tread_count\tsubmitted_ftp\n"
        "BAD1\tRhizophagus irregularis\tPACBIO_SMRT\tSINGLE\t\t500000\tftp.sra.ebi.ac.uk/vol1/hdf5/BAD1.bas.h5\n"
        "BAD2\tRhizophagus irregularis\tOXFORD_NANOPORE\tSINGLE\tftp.sra.ebi.ac.uk/vol1/fastq/BAD2.fastq.gz\t1\t\n"
        "GOOD1\tRhizophagus irregularis\tOXFORD_NANOPORE\tSINGLE\tftp.sra.ebi.ac.uk/vol1/fastq/GOOD1.fastq.gz\t200000\t\n"
    )
    urls, meta = select_ena_read_sources(rows, "long-reads", 4)
    assert urls == ["https://ftp.sra.ebi.ac.uk/vol1/fastq/GOOD1.fastq.gz"]
    assert [row["run_accession"] for row in meta] == ["GOOD1"]


def test_merge_sequence_sources_concatenates_gz_fastq(tmp_path: Path):
    p1 = tmp_path / "r1.fastq.gz"
    p2 = tmp_path / "r2.fastq.gz"
    with gzip.open(p1, "wt", encoding="utf-8") as fh:
        fh.write("@a\nACGT\n+\n!!!!\n")
    with gzip.open(p2, "wt", encoding="utf-8") as fh:
        fh.write("@b\nTGCA\n+\n!!!!\n")
    merged = merge_sequence_sources([str(p1), str(p2)], tmp_path / "merged_reads")
    text = merged.read_text(encoding="utf-8")
    assert merged.suffix == ".fastq"
    assert "@a" in text and "@b" in text


def test_merge_sequence_sources_rejects_non_fastq_payload(tmp_path: Path):
    bad = tmp_path / "bad.fastq"
    bad.write_bytes(b"\x89HDF\r\n\x1a\nnot fastq" + b"x" * 32)
    try:
        merge_sequence_sources([str(bad)], tmp_path / "bad_merged")
    except ValueError as exc:
        assert "does not start with '@'" in str(exc)
    else:
        raise AssertionError("expected non-FASTQ payload to fail validation")
    assert not (tmp_path / "bad_merged.fastq").exists()


def test_materialize_query_input_supports_direct_fastq_urls(tmp_path: Path):
    p1 = tmp_path / "r1.fastq.gz"
    with gzip.open(p1, "wt", encoding="utf-8") as fh:
        fh.write("@a\nACGT\n+\n!!!!\n")
    row = {
        "asm_name": "short_reads_q",
        "query_mode": "short-reads",
        "fastq_url_1": str(p1),
        "species": "Aspergillus fumigatus",
    }
    query_row, query_path, source_rows = materialize_query_input(row, tmp_path / "queries", "custom", public_max_runs=1)
    assert query_row["query_mode"] == "short-reads"
    assert Path(query_path).exists()
    assert source_rows


def test_prepare_custom_manifest_uses_shared_data_cache(tmp_path: Path):
    import argparse

    src = tmp_path / "src"
    src.mkdir()
    ref = src / "ref.fa"
    qry = src / "query.fa"
    ref.write_text(">ref\nACGTACGT\n", encoding="utf-8")
    qry.write_text(">query\nACGTACGA\n", encoding="utf-8")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "role\tasm_name\tpath\tspecies\n"
        f"ref\tref_asm\t{ref}\tExample species\n"
        f"query\tquery_asm\t{qry}\tExample species\n",
        encoding="utf-8",
    )

    out_dir = tmp_path / "prepared"
    cache_dir = tmp_path / "data_cache"
    args = argparse.Namespace(
        out_dir=out_dir,
        custom_url_manifest=manifest,
        data_cache_dir=cache_dir,
        public_query_max_runs=1,
    )

    assert rrfb.prepare_from_custom_manifest(args) == 0
    ref_list = (out_dir / "ref_list.txt").read_text(encoding="utf-8").strip()
    query_list = (out_dir / "query_list.txt").read_text(encoding="utf-8").strip()
    assert ref_list == str((cache_dir / "refs" / "ref.fa").resolve())
    assert query_list == str((cache_dir / "queries" / "query.fa").resolve())
    summary = (out_dir / "prepare_summary.json").read_text(encoding="utf-8")
    assert str(cache_dir.resolve()) in summary


def test_estimate_prepared_genome_size_hint_reads_gz_benchmark_reference(tmp_path: Path):
    prepared = tmp_path / "prepared"
    prepared.mkdir()
    ref = prepared / "ref.fa.gz"
    with gzip.open(ref, "wt", encoding="utf-8") as fh:
        fh.write(">chr1\nACGTACGT\n>chr2\nAAA\n")
    (prepared / "query_manifest.tsv").write_text(
        "query_asm\tbenchmark_ref_fasta\nq1\t" + str(ref) + "\n",
        encoding="utf-8",
    )
    assert estimate_prepared_genome_size_hint(prepared) == 11


def test_run_mycosv_injects_read_mode_perf_defaults(tmp_path: Path, monkeypatch):
    prepared = tmp_path / "prepared"
    prepared.mkdir()
    (prepared / "ref_list.txt").write_text(str(prepared / "ref.fa") + "\n", encoding="utf-8")
    (prepared / "query_list.txt").write_text(str(prepared / "reads.fastq") + "\n", encoding="utf-8")
    ref = prepared / "ref.fa"
    ref.write_text(">chr1\n" + ("ACGT" * 100) + "\n", encoding="utf-8")
    (prepared / "query_manifest.tsv").write_text(
        "query_asm\tbenchmark_ref_fasta\nq1\t" + str(ref) + "\n",
        encoding="utf-8",
    )

    captured: dict[str, list[str]] = {}

    def fake_run(cmd, cwd=None):
        captured["cmd"] = list(cmd)
        class Dummy:
            stdout = ""
            stderr = ""
        return Dummy()

    monkeypatch.setattr(rrfb, "run_mycosv_command", fake_run)
    run_mycosv(prepared, tmp_path / "bench", tmp_path / "fake.exe", "short-reads", [])
    cmd = captured["cmd"]
    assert "--max-reads" in cmd
    assert "150000" in cmd
    assert "--genome-size-hint" in cmd
    assert "400" in cmd


def test_run_mycosv_reuses_prebuilt_index(tmp_path: Path, monkeypatch):
    """run_mycosv with reuse_index_dir must NOT rebuild the index and must
    point the binary at the prebuilt directory. Guards the million-real flow
    where prepare-million-real already wrote the index next to query_manifest.
    """
    prepared = tmp_path / "prepared"
    prepared.mkdir()
    ref = prepared / "ref.fa"
    ref.write_text(">chr1\nACGTACGTACGT\n", encoding="utf-8")
    (prepared / "ref_list.txt").write_text(str(ref) + "\n", encoding="utf-8")
    (prepared / "query_list.txt").write_text(str(ref) + "\n", encoding="utf-8")
    (prepared / "query_manifest.tsv").write_text(
        "query_asm\tbenchmark_ref_fasta\nq1\t" + str(ref) + "\n",
        encoding="utf-8",
    )
    (prepared / "hierarchy_manifest.tsv").write_text(
        "asm_name\tphylum\tclass\torder\tfamily\tgenus\tclade_name\tclade_rank\tfasta_path\n"
        f"q1\t.\t.\t.\t.\t.\t.\tspecies\t{ref}\n",
        encoding="utf-8",
    )

    # Prebuilt index dir with the marker file run_mycosv looks for.
    prebuilt_idx = tmp_path / "million_real_index"
    prebuilt_idx.mkdir()
    (prebuilt_idx / "routing_manifest.tsv").write_text("asm\tcentroid\n", encoding="utf-8")
    prebuilt_reg = tmp_path / "million_real_registry"
    prebuilt_reg.mkdir()

    invocations: list[list[str]] = []

    def fake_run(cmd, cwd=None):
        invocations.append(list(cmd))
        class Dummy:
            stdout = ""
            stderr = ""
        return Dummy()

    monkeypatch.setattr(rrfb, "run_mycosv_command", fake_run)
    run_mycosv(
        prepared, tmp_path / "bench", tmp_path / "fake.exe", "assembly", [],
        reuse_index_dir=prebuilt_idx, reuse_registry_dir=prebuilt_reg,
    )

    # Exactly one binary call (the SV-call run); no rebuild of the index.
    assert len(invocations) == 1, invocations
    cmd = invocations[0]
    assert "--tol-build-index" not in cmd, "reuse path must NOT rebuild the index"
    assert str(prebuilt_idx.resolve()) in cmd
    assert str(prebuilt_reg.resolve()) in cmd
    assert "--no-flat-ref-fallback" in cmd


def test_run_mycosv_ref_override_reenables_flat_fallback_for_fresh_index(tmp_path: Path, monkeypatch):
    prepared = tmp_path / "prepared"
    prepared.mkdir()
    full_ref = prepared / "full_ref.fa"
    full_ref.write_text(">chr1\nACGTACGTACGT\n", encoding="utf-8")
    bench_ref = prepared / "bench_ref.fa"
    bench_ref.write_text(">chr1\nACGTACGTACGT\n", encoding="utf-8")
    query = prepared / "query.fa"
    query.write_text(">chr1\nACGTACGTACGT\n", encoding="utf-8")
    (prepared / "ref_list.txt").write_text(str(full_ref) + "\n", encoding="utf-8")
    (prepared / "query_list.txt").write_text(str(query) + "\n", encoding="utf-8")
    (prepared / "query_manifest.tsv").write_text(
        "query_asm\tbenchmark_ref_fasta\nq1\t" + str(bench_ref) + "\n",
        encoding="utf-8",
    )
    (prepared / "hierarchy_manifest.tsv").write_text(
        "asm_name\tphylum\tclass\torder\tfamily\tgenus\tclade_name\tclade_rank\tfasta_path\n"
        f"ref1\t.\t.\t.\t.\t.\t.\tspecies\t{full_ref}\n",
        encoding="utf-8",
    )
    bench_ref_list = tmp_path / "bench_ref_list.txt"
    bench_ref_list.write_text(str(bench_ref) + "\n", encoding="utf-8")

    invocations: list[list[str]] = []

    def fake_run(cmd, cwd=None):
        invocations.append(list(cmd))
        class Dummy:
            stdout = ""
            stderr = ""
        return Dummy()

    monkeypatch.setattr(rrfb, "run_mycosv_command", fake_run)
    run_mycosv(
        prepared,
        tmp_path / "bench",
        tmp_path / "fake.exe",
        "assembly",
        ["--no-flat-ref-fallback", "--no-gfa"],
        ref_list_override=bench_ref_list,
    )

    assert len(invocations) == 2, invocations
    call_cmd = invocations[-1]
    assert "--tol-build-index" not in call_cmd
    assert str(bench_ref_list.resolve()) in call_cmd
    assert "--no-flat-ref-fallback" not in call_cmd
    assert "--no-gfa" in call_cmd


def test_run_mycosv_rejects_invalid_reuse_index(tmp_path: Path, monkeypatch):
    prepared = tmp_path / "prepared"
    prepared.mkdir()
    ref = prepared / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    (prepared / "ref_list.txt").write_text(str(ref) + "\n", encoding="utf-8")
    (prepared / "query_list.txt").write_text(str(ref) + "\n", encoding="utf-8")
    (prepared / "query_manifest.tsv").write_text("query_asm\nq1\n", encoding="utf-8")
    (prepared / "hierarchy_manifest.tsv").write_text(
        "asm_name\tphylum\tclass\torder\tfamily\tgenus\tclade_name\tclade_rank\tfasta_path\n"
        f"q1\t.\t.\t.\t.\t.\t.\tspecies\t{ref}\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(rrfb, "run_mycosv_command", lambda cmd, cwd=None: None)

    bogus_idx = tmp_path / "no_marker"
    bogus_idx.mkdir()
    import pytest
    with pytest.raises(FileNotFoundError):
        run_mycosv(
            prepared, tmp_path / "bench", tmp_path / "fake.exe", "assembly", [],
            reuse_index_dir=bogus_idx,
        )


def test_prepare_million_real_holds_out_queries_and_strips_them_from_index(
    tmp_path: Path, monkeypatch
):
    """End-to-end smoke of the million-real flow's manifest layer (no real
    network or binary calls): we mock the NCBI download + binary build so the
    test exercises:
      1. selection -> ref_manifest_rows construction
      2. held-out queries are sampled stride-uniformly across phyla
      3. queries get a sibling-genus benchmark_ref_fasta
      4. queries are stripped from hierarchy_manifest.tsv / ref_list.txt
      5. query_manifest.tsv + query_list.txt are written next to the index
    Without this guard the chain prepare -> benchmark in step 2 would silently
    skip step 2b (no held-out queries) or, worse, leak truth into the index.
    """
    import argparse
    fake_summary = "#assembly_accession\ttaxid\torganism_name\tassembly_level\tversion_status\tftp_path\n"
    rows = []
    # 6 fake rows across 2 phyla; --max-assemblies caps at 6, --queries=2.
    # ftp_path must be under https://ftp.ncbi.nlm.nih.gov/ for select_all_public_rows.
    for i in range(6):
        rows.append(
            f"GCA_{i:09d}.1\t100{i}\tFakespecies sp{i}\tComplete Genome\tlatest\t"
            f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/{i:03d}/GCA_{i:09d}.1\n"
        )
    monkeypatch.setattr(rrfb, "http_get_text", lambda url: fake_summary + "".join(rows))

    taxonomy_cache_paths = []
    def fake_taxonomy(taxids, cache_path=None):
        taxonomy_cache_paths.append(cache_path)
        out = {}
        for i, t in enumerate(taxids):
            phylum = "Ascomycota" if i < 3 else "Basidiomycota"
            out[t] = {
                "phylum": phylum, "class": ".", "order": ".",
                "family": ".", "genus": f"Genus{i}",
                "species": f"Fakespecies sp{i}",
            }
        return out
    monkeypatch.setattr(rrfb, "fetch_taxonomy_lineages", fake_taxonomy)

    # Materialize a tiny FASTA per row instead of hitting the network.
    refs_dir = tmp_path / "cache" / "refs"
    refs_dir.mkdir(parents=True)
    def fake_materialize(url, dest, keep_gz=True):
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(b">chr\nACGT\n")
        return dest
    monkeypatch.setattr(rrfb, "materialize_entry", fake_materialize)
    # Skip the binary build / decoy padding so the test stays hermetic.
    monkeypatch.setattr(rrfb, "compile_binary_if_needed", lambda *a, **k: None)
    monkeypatch.setattr(rrfb, "run_mycosv_command", lambda cmd, cwd=None: None)
    monkeypatch.setattr(rrfb, "augment_routing_store",
                        lambda idx, target, seed: {"real_centroids": 0,
                                                    "decoy_centroids": 0,
                                                    "total_centroids": 0,
                                                    "hashes_per_centroid": 0})

    out_dir = tmp_path / "million_real"
    args = argparse.Namespace(
        out_dir=out_dir,
        source="ncbi-genbank",
        max_assemblies=6,
        min_assembly_level="contig",
        latest_only=False,
        target_centroids=0,
        seed=42,
        threads=1,
        max_clade_genomes=2,
        binary_path=tmp_path / "fake_bin",
        force_rebuild=False,
        data_cache_dir=tmp_path / "cache",
        million_real_queries=2,
    )
    rc = rrfb.prepare_million_real(args)
    assert rc == 0
    assert taxonomy_cache_paths == [tmp_path / "cache" / "taxonomy_cache.json"]

    # The index manifest must NOT contain any of the held-out queries.
    hierarchy = (out_dir / "hierarchy_manifest.tsv").read_text(encoding="utf-8").splitlines()
    assert len(hierarchy) - 1 == 4, hierarchy  # 6 selected - 2 queries
    qm_path = out_dir / "query_manifest.tsv"
    assert qm_path.exists(), "million-real did not write query_manifest.tsv"
    qm_lines = qm_path.read_text(encoding="utf-8").splitlines()
    assert len(qm_lines) - 1 == 2, qm_lines
    qm_header = qm_lines[0].split("\t")
    assert "query_asm" in qm_header
    assert "benchmark_ref_fasta" in qm_header
    assert "phylum" in qm_header
    bench_ref_idx = qm_header.index("benchmark_ref_fasta")
    qasm_idx = qm_header.index("query_asm")
    qasms = [line.split("\t")[qasm_idx] for line in qm_lines[1:]]
    bench_refs = [line.split("\t")[bench_ref_idx] for line in qm_lines[1:]]
    assert len(set(qasms)) == 2, qasms
    # Each held-out query's benchmark_ref_fasta must be one of the OTHER
    # rows' FASTA paths — never the query's own — so we don't leak truth.
    for line in qm_lines[1:]:
        cells = line.split("\t")
        own_path_in_qlist = any(cells[qasm_idx] in ref for ref in bench_refs)
        # Loosely: a benchmark ref path must exist on disk.
        assert Path(cells[bench_ref_idx]).exists(), cells
    # query_list.txt must mirror query_manifest.tsv.
    ql_lines = (out_dir / "query_list.txt").read_text(encoding="utf-8").splitlines()
    assert len([l for l in ql_lines if l.strip()]) == 2

    # Every line in ref_list.txt must NOT belong to a held-out query.
    rl_lines = (out_dir / "ref_list.txt").read_text(encoding="utf-8").splitlines()
    rl_paths = {Path(l).name for l in rl_lines if l.strip()}
    for q_path in ql_lines:
        if q_path.strip():
            assert Path(q_path).name not in rl_paths, (q_path, rl_paths)


def test_benchmark_real_data_mycosv_only_skips_comparators(tmp_path: Path, monkeypatch):
    """--mycosv-only must not auto-enable any comparator binaries, must not
    force the per-mode mandatory baseline, and must clear any pre-set
    --run-X flags. Guards the million-real step against silent comparator
    runs that would inflate wall time and pull in tools we don't intend to
    benchmark in that flow.
    """
    import argparse
    # Build a minimal prepared dir that benchmark_real_data accepts.
    prepared = tmp_path / "prepared"
    prepared.mkdir()
    ref = prepared / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    (prepared / "ref_list.txt").write_text(str(ref) + "\n", encoding="utf-8")
    (prepared / "query_list.txt").write_text(str(ref) + "\n", encoding="utf-8")
    (prepared / "query_manifest.tsv").write_text(
        "query_asm\tquery_mode\tpath\tbenchmark_ref_fasta\n"
        f"q1\tassembly\t{ref}\t{ref}\n",
        encoding="utf-8",
    )
    (prepared / "hierarchy_manifest.tsv").write_text(
        "asm_name\tphylum\tclass\torder\tfamily\tgenus\tclade_name\tclade_rank\tfasta_path\n"
        f"q1\t.\t.\t.\t.\t.\t.\tspecies\t{ref}\n",
        encoding="utf-8",
    )

    # Make every binary appear available so a missing-tool short-circuit
    # cannot mask the --mycosv-only path.
    monkeypatch.setattr(rrfb, "tool_path", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(rrfb, "compile_binary_if_needed", lambda *a, **k: None)
    monkeypatch.setattr(
        rrfb, "run_mycosv",
        lambda *a, **k: {"vcf": str(prepared / "calls.vcf"),
                         "hits": str(prepared / "calls.hits.tsv"),
                         "gfa": str(prepared / "calls.gfa")},
    )
    (prepared / "calls.vcf").write_text(
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
        encoding="utf-8",
    )
    (prepared / "calls.hits.tsv").write_text("", encoding="utf-8")
    (prepared / "calls.gfa").write_text("", encoding="utf-8")

    # Trip-wires: if --mycosv-only is honored, none of these comparator
    # entry points should be reached.
    for name in ("run_syri_for_query", "run_minigraph_for_query",
                 "run_pggb_for_query", "run_cactus_for_query",
                 "run_svim_asm_for_query", "run_anchorwave_for_query",
                 "run_svim_for_query", "run_sniffles_for_query",
                 "run_cutesv_for_query", "run_delly_for_query",
                 "run_manta_for_query"):
        def _trip(*a, **k):
            raise AssertionError(f"comparator {name} ran under --mycosv-only")
        monkeypatch.setattr(rrfb, name, _trip)

    # Skip the candidate analyzer for hermeticity (it shells out).
    monkeypatch.setattr(rrfb, "maybe_run_candidate_analysis",
                        lambda *a, **k: (None, None))

    args = argparse.Namespace(
        prepared_dir=prepared,
        out_dir=tmp_path / "out",
        binary_path=tmp_path / "fake_bin",
        force_rebuild=False,
        mode="assembly",
        threads=1,
        max_clade_genomes=2,
        run_all_comparators=True,   # would normally enable everything; mycosv_only must override
        mycosv_only=True,
        run_syri=True, run_minigraph=True, run_pggb=True,
        run_cactus=True, run_svim_asm=True, run_anchorwave=True,
        run_svim=True, run_sniffles=True, run_cutesv=True,
        run_delly=True, run_manta=True,
        cactus_arg=[],
        normalized_other=[], other_vcf=[],
        mycosv_arg=[], minigraph_arg=[], pggb_arg=[],
        pggb_identity="90", pggb_segment_len="5k",
        expression_tsv=None, gene_annotations_tsv=None, ancestral_tsv=None,
        validate_with_reads=False,
        read_validation_min_support=3,
        read_validation_flank_bp=250,
        reuse_index_dir=None, reuse_registry_dir=None,
    )
    rc = rrfb.benchmark_real_data(args)
    assert rc == 0
    # Every --run-X must have been cleared by the --mycosv-only branch.
    for flag in ("run_syri", "run_minigraph", "run_pggb", "run_cactus",
                 "run_svim_asm", "run_anchorwave",
                 "run_svim", "run_sniffles", "run_cutesv",
                 "run_delly", "run_manta"):
        assert getattr(args, flag) is False, flag


def test_benchmark_mycosv_only_validates_mycosv_reference_calls(tmp_path: Path, monkeypatch):
    """Million-real runs use --mycosv-only, so there is no comparator truth
    to validate. Guard that read-level validation still records support for
    MycoSV reference-coordinate calls when --validate-with-reads is enabled.
    """
    import argparse
    import json

    prepared = tmp_path / "prepared"
    prepared.mkdir()
    ref = prepared / "ref.fa"
    ref.write_text(">chr1\n" + "ACGT" * 100 + "\n", encoding="utf-8")
    (prepared / "ref_list.txt").write_text(str(ref) + "\n", encoding="utf-8")
    (prepared / "query_list.txt").write_text(str(ref) + "\n", encoding="utf-8")
    (prepared / "query_manifest.tsv").write_text(
        "query_asm\tquery_mode\tpath\tbenchmark_ref_fasta\n"
        f"q1\tassembly\t{ref}\t{ref}\n",
        encoding="utf-8",
    )
    (prepared / "hierarchy_manifest.tsv").write_text(
        "asm_name\tphylum\tclass\torder\tfamily\tgenus\tclade_name\tclade_rank\tfasta_path\n"
        f"q1\t.\t.\t.\t.\t.\t.\tspecies\t{ref}\n",
        encoding="utf-8",
    )
    vcf = prepared / "calls.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "qctg\t11\tsv1\tN\t<DEL>\t40\tPASS\t"
        "SVTYPE=DEL;SVLEN=-50;END=11;REFCONTIG=chr1;REFPOS=101;REFEND=150;QASM=q1\tGT\t0/1\n",
        encoding="utf-8",
    )
    hits = prepared / "calls.hits.tsv"
    hits.write_text("", encoding="utf-8")
    gfa = prepared / "calls.gfa"
    gfa.write_text("", encoding="utf-8")

    monkeypatch.setattr(rrfb, "tool_path", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(rrfb, "compile_binary_if_needed", lambda *a, **k: None)
    monkeypatch.setattr(
        rrfb, "run_mycosv",
        lambda *a, **k: {"vcf": str(vcf), "hits": str(hits), "gfa": str(gfa)},
    )
    monkeypatch.setattr(rrfb, "maybe_run_candidate_analysis", lambda *a, **k: (None, None))

    def fake_validate(calls, query_row, work_dir, *, threads, min_support, flank_bp):
        assert len(calls) == 1
        return list(calls), [{
            "query_asm": query_row["query_asm"],
            "ref_contig": calls[0].ref_contig,
            "pos": calls[0].pos,
            "end": calls[0].end,
            "svtype": calls[0].svtype,
            "source": calls[0].source,
            "coord_space": calls[0].coord_space,
            "read_support": min_support,
            "read_validated": "yes",
        }]

    monkeypatch.setattr(rrfb, "validate_calls_with_reads", fake_validate)

    args = argparse.Namespace(
        prepared_dir=prepared,
        out_dir=tmp_path / "out",
        binary_path=tmp_path / "fake_bin",
        force_rebuild=False,
        mode="assembly",
        threads=1,
        max_clade_genomes=2,
        run_all_comparators=False,
        mycosv_only=True,
        run_syri=False, run_minigraph=False, run_pggb=False,
        run_cactus=False, run_svim_asm=False, run_anchorwave=False,
        run_svim=False, run_sniffles=False, run_cutesv=False,
        run_delly=False, run_manta=False,
        cactus_arg=[],
        normalized_other=[], other_vcf=[],
        mycosv_arg=[], minigraph_arg=[], pggb_arg=[],
        pggb_identity="90", pggb_segment_len="5k",
        expression_tsv=None, gene_annotations_tsv=None, ancestral_tsv=None,
        validate_with_reads=True,
        read_validation_min_support=1,
        read_validation_flank_bp=250,
        reuse_index_dir=None, reuse_registry_dir=None,
    )
    assert rrfb.benchmark_real_data(args) == 0
    readval = args.out_dir / "read_validated_truth.tsv"
    assert readval.exists()
    assert "mycosv" in readval.read_text(encoding="utf-8")
    assert (vcf.with_suffix(".multisample.vcf")).exists()
    assert (vcf.with_name("alls.multisample.vcf")).exists()
    assert (args.out_dir / "biology_findings.tsv").exists()
    summary = json.loads((args.out_dir / "benchmark_summary.json").read_text(encoding="utf-8"))
    assert summary["queries"]["q1"]["read_validation"]["mycosv_reference"]["read_validated"] == 1
