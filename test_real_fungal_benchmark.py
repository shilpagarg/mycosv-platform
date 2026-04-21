#!/usr/bin/env python3
# Designed for Linux

import run_real_fungal_benchmark as rrfb
from pathlib import Path
import gzip

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
        "run_accession\tscientific_name\tinstrument_platform\tlibrary_layout\tfastq_ftp\tsubmitted_ftp\n"
        "SRR1\tAspergillus fumigatus\tILLUMINA\tPAIRED\tftp.sra.ebi.ac.uk/vol1/fastq/SRR1_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR1_2.fastq.gz\t\n"
        "SRR2\tAspergillus fumigatus\tOXFORD_NANOPORE\tSINGLE\tftp.sra.ebi.ac.uk/vol1/fastq/SRR2.fastq.gz\t\n"
    )
    rows = parse_ena_filereport_text(text)
    urls, meta = select_ena_read_sources(rows, "short-reads", 2)
    assert len(urls) == 2
    assert all(url.startswith("https://ftp.sra.ebi.ac.uk/") for url in urls)
    assert meta[0]["run_accession"] == "SRR1"


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
