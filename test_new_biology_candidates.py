#!/usr/bin/env python3
# Designed for Linux

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / 'analyze_new_biology_candidates.py'


def run(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    if cmd and cmd[0] == 'python3':
        cmd = cmd.copy()
        cmd[0] = sys.executable
    return subprocess.run(cmd, text=True, capture_output=True, check=True)


def test_new_biology_candidate_report_smoke(tmp_path: Path) -> None:
    vcf = tmp_path / 'calls.vcf'
    vcf.write_text(
        '##fileformat=VCFv4.3\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
        'ctgA\t100\t.\tN\t<OFF_REF>\t60\tPASS\tSVTYPE=OFF_REF;END=400;ANNOT=NOVEL;EC=LTR_GYPSY\tGT\t0/1\n'
        'ctgB\t250\t.\tN\t<INV>\t60\tPASS\tSVTYPE=INV;END=650;ANNOT=NONE;EC=NONE\tGT\t0/1\n'
    )
    hits = tmp_path / 'calls.hits.tsv'
    hits.write_text(
        'query_asm\tquery_contig\ttype\tref_asm\tref_contig\tpos\tend\tsvlen\tblock_score\tanchors\tgenotype\tgq\tannotation\talignment_mode\tquery_mode\n'
        'asm1\tctgA\tOFF_REF\tOFF_REFERENCE\t.\t100\t400\t300\t8\t0\t0/1\t20\tNOVEL\tsimple_offref_fallback\tassembly\n'
        'asm2\tctgB\tINV\tref1\tchr1\t250\t650\t400\t15\t4\t0/1\t30\tNONE\tmem_chain\tassembly\n'
    )
    anc = tmp_path / 'calls.ancestral.tsv'
    anc.write_text(
        'query_asm\tquery_contig\tclade\tclade_rank\tphylum\tvariant_type\tbreakpoints\tsegment_bp\n'
        'asm1\tctgA\tClade1\tgenus\tAscomycota\tOFF_REF\t.\t300\n'
        'asm1\tctgA\tClade2\tfamily\tAscomycota\tOFF_REF\t100-220\t150\n'
    )
    meta = tmp_path / 'query_metadata.tsv'
    meta.write_text(
        'query_asm\tscenario\tlifestyle\tarchitecture\n'
        'asm1\thgt_receiver\thgt_receiver\tgc_shifted_hgt_receiver\n'
        'asm2\tcore\tbaseline_control\tcompact_baseline\n'
    )
    out_tsv = tmp_path / 'new_biology.tsv'
    out_json = tmp_path / 'new_biology.json'
    run([
        'python3', str(SCRIPT),
        '--vcf', str(vcf),
        '--hits', str(hits),
        '--ancestral', str(anc),
        '--query-metadata', str(meta),
        '--out-tsv', str(out_tsv),
        '--summary-json', str(out_json),
        '--phylum', 'Ascomycota',
        '--top-n', '10',
    ])
    rows = list(csv.DictReader(out_tsv.open(), delimiter='\t'))
    assert rows, 'expected at least one candidate row'
    assert rows[0]['candidate_type'] == 'novel_te_architecture'
    assert rows[0]['element_class'] == 'LTR_GYPSY'
    assert rows[0]['scenario'] == 'hgt_receiver'
    assert rows[0]['functional_example'] == 'methylation-silenced Hop insertion at the b1 locus'
    assert rows[0]['evidence_axis'] == 'expression_epigenetic'
    assert 'allele-specific silencing or activation' in rows[0]['real_data_signal']
    summary = json.loads(out_json.read_text())
    assert summary['candidate_count'] >= 1
    assert 'novel_te_architecture' in summary['by_candidate_type']
    assert summary['by_evidence_axis']['expression_epigenetic'] >= 1


def test_new_biology_candidate_report_accepts_caller_te_aliases(tmp_path: Path) -> None:
    vcf = tmp_path / 'calls.vcf'
    vcf.write_text(
        '##fileformat=VCFv4.3\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
        'ctgA\t100\t.\tN\t<OFF_REF>\t60\tPASS\tSVTYPE=OFF_REF;END=400;ANNOT=NOVEL;EC=TE_LTR\tGT\t0/1\n'
    )
    hits = tmp_path / 'calls.hits.tsv'
    hits.write_text(
        'query_asm\tquery_contig\ttype\tref_asm\tref_contig\tpos\tend\tsvlen\tblock_score\tanchors\tgenotype\tgq\tannotation\talignment_mode\tquery_mode\n'
        'asm1\tctgA\tOFF_REF\tOFF_REFERENCE\t.\t100\t400\t300\t8\t0\t0/1\t20\tNOVEL\tsimple_offref_fallback\tassembly\n'
    )
    meta = tmp_path / 'query_metadata.tsv'
    meta.write_text(
        'query_asm\tscenario\tlifestyle\tarchitecture\n'
        'asm1\thgt_receiver\thgt_receiver\tgc_shifted_hgt_receiver\n'
    )
    out_tsv = tmp_path / 'new_biology.tsv'
    out_json = tmp_path / 'new_biology.json'
    run([
        'python3', str(SCRIPT),
        '--vcf', str(vcf),
        '--hits', str(hits),
        '--query-metadata', str(meta),
        '--out-tsv', str(out_tsv),
        '--summary-json', str(out_json),
        '--phylum', 'Ascomycota',
        '--top-n', '10',
    ])
    rows = list(csv.DictReader(out_tsv.open(), delimiter='\t'))
    assert rows, 'expected TE_LTR to be recognized as a TE-linked candidate'
    assert rows[0]['candidate_type'] == 'novel_te_architecture'
    assert rows[0]['element_class'] == 'TE_LTR'
    assert rows[0]['functional_example'] == 'methylation-silenced Hop insertion at the b1 locus'


def test_new_biology_candidate_report_resolves_qasm_and_downsample_alias(tmp_path: Path) -> None:
    vcf = tmp_path / 'calls.vcf'
    vcf.write_text(
        '##fileformat=VCFv4.3\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
        'ctgA\t100\t.\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=250;ANNOT=NONE;EC=NONE;QASM=asm1.20000\tGT\t0/1\n'
    )
    hits = tmp_path / 'calls.hits.tsv'
    hits.write_text(
        'query_asm\tquery_contig\ttype\tref_asm\tref_contig\tpos\tend\tsvlen\tblock_score\tanchors\tgenotype\tgq\tannotation\talignment_mode\tquery_mode\n'
    )
    meta = tmp_path / 'query_metadata.tsv'
    meta.write_text(
        'query_asm\tscenario\tlifestyle\tarchitecture\n'
        'asm1\tte_rich_pathogen\tplant_pathogen\tte_rich\n'
    )
    out_tsv = tmp_path / 'new_biology.tsv'
    out_json = tmp_path / 'new_biology.json'
    run([
        'python3', str(SCRIPT),
        '--vcf', str(vcf),
        '--hits', str(hits),
        '--query-metadata', str(meta),
        '--out-tsv', str(out_tsv),
        '--summary-json', str(out_json),
        '--phylum', 'Basidiomycota',
        '--top-n', '10',
    ])
    rows = list(csv.DictReader(out_tsv.open(), delimiter='\t'))
    assert rows
    assert rows[0]['query_asm'] == 'asm1.20000'
    assert rows[0]['scenario'] == 'te_rich_pathogen'
    assert rows[0]['candidate_type'] == 'structural_sv_signal'


def test_new_biology_candidate_report_keeps_svtype_diversity(tmp_path: Path) -> None:
    vcf = tmp_path / 'calls.vcf'
    lines = [
        '##fileformat=VCFv4.3\n',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n',
    ]
    for i in range(20):
        lines.append(
            f'off{i}\t1\t.\tN\t<OFF_REF>\t60\tPASS\tSVTYPE=OFF_REF;END=300;ANNOT=NOVEL;EC=RIP;QASM=asm1\tGT\t0/1\n'
        )
    lines.append('del1\t500\t.\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=700;ANNOT=NONE;EC=NONE;QASM=asm1\tGT\t0/1\n')
    lines.append('inv1\t900\t.\tN\t<INV>\t60\tPASS\tSVTYPE=INV;END=1200;ANNOT=NONE;EC=NONE;QASM=asm1\tGT\t0/1\n')
    vcf.write_text(''.join(lines))
    hits = tmp_path / 'calls.hits.tsv'
    hits.write_text(
        'query_asm\tquery_contig\ttype\tref_asm\tref_contig\tpos\tend\tsvlen\tblock_score\tanchors\tgenotype\tgq\tannotation\talignment_mode\tquery_mode\n'
    )
    meta = tmp_path / 'query_metadata.tsv'
    meta.write_text('query_asm\tscenario\tlifestyle\tarchitecture\nasm1\tpathogenic\tplant_pathogen\tte_rich\n')
    out_tsv = tmp_path / 'new_biology.tsv'
    out_json = tmp_path / 'new_biology.json'
    run([
        'python3', str(SCRIPT),
        '--vcf', str(vcf),
        '--hits', str(hits),
        '--query-metadata', str(meta),
        '--out-tsv', str(out_tsv),
        '--summary-json', str(out_json),
        '--phylum', 'Ascomycota',
        '--top-n', '6',
    ])
    svtypes = {row['svtype'] for row in csv.DictReader(out_tsv.open(), delimiter='\t')}
    assert {'OFF_REF', 'DEL', 'INV'} <= svtypes


def test_new_biology_candidate_report_uses_direct_expression_support(tmp_path: Path) -> None:
    vcf = tmp_path / 'calls.vcf'
    vcf.write_text(
        '##fileformat=VCFv4.3\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
        'ctgA\t100\t.\tN\t<INS>\t60\tPASS\tSVTYPE=INS;END=140;ANNOT=NONE;EC=TE_LTR\tGT\t0/1\n'
    )
    hits = tmp_path / 'calls.hits.tsv'
    hits.write_text(
        'query_asm\tquery_contig\ttype\tref_asm\tref_contig\tpos\tend\tsvlen\tblock_score\tanchors\tgenotype\tgq\tannotation\talignment_mode\tquery_mode\n'
        'asm1\tctgA\tINS\tref1\tchr1\t100\t140\t40\t18\t4\t0/1\t40\tNONE\tmem_chain\tassembly\n'
    )
    meta = tmp_path / 'query_metadata.tsv'
    meta.write_text(
        'query_asm\tscenario\tlifestyle\tarchitecture\n'
        'asm1\tpathogenic\tplant_pathogen\ttwo_speed_pathogen_extreme\n'
    )
    expr = tmp_path / 'expression.tsv'
    expr.write_text(
        'query_asm\tquery_contig\tpos\tend\tsvtype\tgene_id\tgene_name\tdistance_bp\tlog2_fc\tpadj\tcondition\n'
        'asm1\tctgA\t100\t140\tINS\tgene_1\tEFF1\t350\t2.4\t0.003\tstress\n'
    )
    out_tsv = tmp_path / 'new_biology.tsv'
    out_json = tmp_path / 'new_biology.json'
    run([
        'python3', str(SCRIPT),
        '--vcf', str(vcf),
        '--hits', str(hits),
        '--query-metadata', str(meta),
        '--expression-tsv', str(expr),
        '--out-tsv', str(out_tsv),
        '--summary-json', str(out_json),
        '--phylum', 'Ascomycota',
        '--top-n', '10',
    ])
    rows = list(csv.DictReader(out_tsv.open(), delimiter='\t'))
    assert rows, 'expected expression-supported candidate row'
    assert rows[0]['candidate_type'] == 'te_expression_link'
    assert rows[0]['evidence_axis'] == 'expression_direct'
    assert rows[0]['expression_supported'] == 'yes'
    assert rows[0]['expression_gene'] == 'EFF1'
    summary = json.loads(out_json.read_text())
    assert summary['by_expression_support']['yes'] >= 1


def test_new_biology_candidate_report_quantifies_expression_from_sample_level_table(tmp_path: Path) -> None:
    vcf = tmp_path / 'calls.vcf'
    vcf.write_text(
        '##fileformat=VCFv4.3\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
        'ctgA\t100\t.\tN\t<INS>\t60\tPASS\tSVTYPE=INS;END=140;ANNOT=NONE;EC=TE_LTR\tGT\t0/1\n'
    )
    hits = tmp_path / 'calls.hits.tsv'
    hits.write_text(
        'query_asm\tquery_contig\ttype\tref_asm\tref_contig\tpos\tend\tsvlen\tblock_score\tanchors\tgenotype\tgq\tannotation\talignment_mode\tquery_mode\n'
        'asm1\tctgA\tINS\tref1\tchr1\t100\t140\t40\t18\t4\t0/1\t40\tNONE\tmem_chain\tassembly\n'
    )
    meta = tmp_path / 'query_metadata.tsv'
    meta.write_text(
        'query_asm\tscenario\tlifestyle\tarchitecture\n'
        'asm1\tpathogenic\tplant_pathogen\ttwo_speed_pathogen_extreme\n'
    )
    genes = tmp_path / 'genes.tsv'
    genes.write_text(
        'query_asm\tquery_contig\tgene_id\tgene_name\tstart\tend\n'
        'asm1\tctgA\tgene_1\tEFF1\t260\t420\n'
    )
    expr_long = tmp_path / 'expression_long.tsv'
    expr_long.write_text(
        'sample\tquery_asm\tgene_id\texpression\tcondition\n'
        'c1\tasm1\tgene_1\t1.0\tcontrol\n'
        'c2\tasm1\tgene_1\t1.2\tcontrol\n'
        's1\tasm1\tgene_1\t24.0\tstress\n'
        's2\tasm1\tgene_1\t20.0\tstress\n'
    )
    derived = tmp_path / 'derived_expression.tsv'
    out_tsv = tmp_path / 'new_biology.tsv'
    out_json = tmp_path / 'new_biology.json'
    run([
        'python3', str(SCRIPT),
        '--vcf', str(vcf),
        '--hits', str(hits),
        '--query-metadata', str(meta),
        '--expression-long-tsv', str(expr_long),
        '--gene-annotations', str(genes),
        '--derived-expression-out', str(derived),
        '--expression-group-a', 'control',
        '--expression-group-b', 'stress',
        '--out-tsv', str(out_tsv),
        '--summary-json', str(out_json),
        '--phylum', 'Ascomycota',
        '--top-n', '10',
    ])
    rows = list(csv.DictReader(out_tsv.open(), delimiter='\t'))
    assert rows, 'expected direct-quant candidate row'
    assert rows[0]['candidate_type'] == 'te_expression_link'
    assert rows[0]['expression_supported'] == 'yes'
    assert rows[0]['expression_gene'] == 'EFF1'
    assert rows[0]['expression_condition'] == 'stress_vs_control'
    derived_rows = list(csv.DictReader(derived.open(), delimiter='\t'))
    assert derived_rows, 'expected derived expression support output'
    assert derived_rows[0]['gene_name'] == 'EFF1'
    assert derived_rows[0]['condition'] == 'stress_vs_control'
