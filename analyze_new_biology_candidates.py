#!/usr/bin/env python3
# Designed for Linux

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

NOVEL_TIERS = {"NOVEL", "NOVEL_WEAK", "DIVERGED", "OFF_REF_KNOWN"}

# The C++ caller emits coarse element classes such as TE_LTR / TE_TIR / TE_LINE,
# while some curated examples in this report are keyed by finer biological labels
# such as LTR_GYPSY, DNA_TIR, or the generic TE bucket.  Normalizing the emitted
# labels keeps the report biologically consistent with the caller instead of
# silently downgrading true TE-linked calls to non-TE "other" rows.
ELEMENT_CLASS_ALIASES = {
    "TE_LTR": "TE",
    "TE_TIR": "DNA_TIR",
    "TE_LINE": "LINE",
    "TE_SINE": "SINE",
}
TE_CLASSES = {"TE", "LTR_GYPSY", "LTR_COPIA", "LINE", "SINE", "DNA_TIR", "HELITRON", "MITE", "RIP", "STARSHIP", "HGT", "REPEAT", "TE_LTR", "TE_TIR", "TE_LINE", "TE_SINE"}

# MGE subtypes for Mobile Genetic Element breakdown.
# Separated from the TE_CLASSES superset to allow independent MGE reporting.
MGE_INTEGRATIVE = {"STARSHIP", "HGT"}    # integrative island-type MGEs
MGE_TRANSPOSABLE = {"TE", "LTR_GYPSY", "LTR_COPIA", "DNA_TIR", "HELITRON", "MITE",
                    "TE_LTR", "TE_TIR", "LINE", "SINE"}
MGE_REPEAT_BASED = {"REPEAT", "RIP"}     # repeat/RIP elements (not strictly MGE)

# Curated functional exemplars used to make the report more actionable.
# These are short, concrete analogies that tie a structural signal to a real-data
# biological pattern reviewers will recognize: expression changes, chromatin
# shifts, or adaptive cargo movement.
FUNCTIONAL_EXAMPLES: list[dict[str, object]] = [
    {
        'name': 'methylation-silenced Hop insertion at the b1 locus',
        'match_ec': {'LTR_GYPSY', 'LTR_COPIA', 'TE'},
        'match_candidate_types': {'novel_te_architecture', 'te_architecture_rewiring'},
        'match_svtypes': {'INS', 'OFF_REF', 'DUP'},
        'priority': 5,
        'evidence_axis': 'expression_epigenetic',
        'system': 'maize b1 / Hopscotch-style TE promoter rewiring analogy',
        'real_data_signal': 'TE insertion or nearby repeat amplification can create allele-specific silencing or activation of adjacent genes.',
        'why_relevant': 'Best fit for novel TE insertions or TE-linked duplications near genes where RNA-seq and methylation can test cis-regulatory effects.',
        'suggested_readout': 'RNA-seq/qPCR across conditions plus methylation or chromatin profiling at the adjacent locus.',
    },
    {
        'name': 'stress-responsive LTR activation next to effector-like loci',
        'match_ec': {'LTR_GYPSY', 'LTR_COPIA', 'RIP', 'TE'},
        'match_candidate_types': {'novel_te_architecture', 'te_architecture_rewiring'},
        'match_scenarios': {'pathogenic', 'two_speed_pathogen_extreme', 'necrotrophic', 'saprotrophic'},
        'priority': 4,
        'evidence_axis': 'expression_stress_response',
        'system': 'fungal two-speed genome / effector-proximal LTR analogy',
        'real_data_signal': 'TE-rich compartments in fungal pathogens are often associated with stress-inducible genes and altered transcript output during host interaction.',
        'why_relevant': 'Fits TE-linked insertions, inversions, or duplications in pathogen-like scenarios where repeats can expose genes to inducible chromatin states.',
        'suggested_readout': 'Condition-matched RNA-seq under host-mimic or stress exposure, then test whether nearby genes are differentially expressed.',
    },
    {
        'name': 'mobile Starship-style accessory cargo movement',
        'match_ec': {'STARSHIP', 'HGT'},
        'match_candidate_types': {'novel_te_architecture', 'mosaic_te_lineage_switch'},
        'match_svtypes': {'OFF_REF', 'INS', 'TRA'},
        'priority': 5,
        'evidence_axis': 'adaptive_cargo_transfer',
        'system': 'fungal Starship / accessory gene island analogy',
        'real_data_signal': 'Large mobile cargo-bearing elements can shuttle metabolic or niche-adaptive genes between genomic backgrounds and create lineage-restricted accessory modules.',
        'why_relevant': 'Strong fit for off-reference cargo, translocations, or multi-clade ancestry around HGT- or Starship-like sequence.',
        'suggested_readout': 'Inspect cargo genes, test presence/absence across clades, and profile expression of linked metabolic or stress-response genes.',
    },
    {
        'name': 'TE-seeded recombination breakpoint that rewires regulatory neighborhoods',
        'match_ec': {'LTR_GYPSY', 'LTR_COPIA', 'DNA_TIR', 'HELITRON', 'MITE', 'TE', 'REPEAT'},
        'match_candidate_types': {'mosaic_te_lineage_switch', 'mosaic_lineage_switch', 'te_architecture_rewiring'},
        'match_svtypes': {'INV', 'TRA', 'DUP', 'DEL'},
        'require_ancestral_breakpoints': True,
        'priority': 5,
        'evidence_axis': 'structural_regulatory_rewiring',
        'system': 'repeat-mediated ectopic recombination analogy',
        'real_data_signal': 'Homologous or microhomology-mediated recombination between repeats can reshape gene neighborhoods, copy number, and local regulatory context.',
        'why_relevant': 'Best fit when the call already has breakpoint-resolved ancestry or mixed-clade support, suggesting repeat-guided rearrangement.',
        'suggested_readout': 'Validate breakpoints with long reads, inspect neighboring genes, and test for copy-number or expression shifts across the breakpoint.',
    },
    {
        'name': 'direct local expression shift linked to a TE-associated structural event',
        'match_ec': {'LTR_GYPSY', 'LTR_COPIA', 'DNA_TIR', 'HELITRON', 'MITE', 'TE', 'LINE', 'SINE', 'STARSHIP', 'HGT'},
        'match_candidate_types': {'te_expression_link', 'mosaic_te_expression_link'},
        'priority': 6,
        'require_expression_support': True,
        'evidence_axis': 'expression_direct',
        'system': 'direct RNA-seq/qPCR-supported cis-regulatory SV follow-up',
        'real_data_signal': 'A nearby gene already shows a significant expression shift, so the structural event becomes a direct mechanistic candidate rather than only a novelty hypothesis.',
        'why_relevant': 'Best fit when the caller and expression data agree on a local TE-linked or ancestry-shifted structural event near a significantly changed gene.',
        'suggested_readout': 'Validate the breakpoint and quantify allele- or condition-specific expression around the affected gene set.',
    },
    {
        'name': 'inter-phylum HGT translocation breakpoint',
        'match_ec': {'HGT', 'STARSHIP'},
        'match_candidate_types': {'hgt_candidate'},
        'match_svtypes': {'TRA', 'OFF_REF'},
        'priority': 7,
        'evidence_axis': 'horizontal_gene_transfer',
        'system': 'cross-phylum HGT translocation (Rhizophagus / Puccinia accessory island analogy)',
        'real_data_signal': 'Translocation breakpoints with HGT-class sequence and low same-clade overlap indicate donor-recipient boundaries of an horizontally transferred genomic island.',
        'why_relevant': 'Highest-priority hit when a TRA or off-reference segment is both phylogenetically novel within its clade and carried in a sequence class (HGT / Starship) associated with lateral transfer.',
        'suggested_readout': 'Confirm breakpoint with long reads, survey presence/absence across species panel, BLAST cargo genes against fungal + prokaryotic databases, and test expression of cargo under relevant conditions.',
    },
]


def parse_info(field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in field.split(';'):
        if not item:
            continue
        if '=' in item:
            k, v = item.split('=', 1)
            out[k] = v
        else:
            out[item] = '1'
    return out


def load_vcf_records(path: Path) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            chrom, pos, vid, ref, alt, qual, flt, info_field, *_ = line.rstrip('\n').split('\t')
            info = parse_info(info_field)
            svtype = info.get('SVTYPE', alt.strip('<>'))
            end = info.get('END', pos)
            records.append({
                'chrom': chrom,
                'pos': pos,
                'pos_int': int(pos),
                'end': end,
                'end_int': int(end),
                'info': info,
                'svtype': svtype,
            })
    return records


def normalize_element_class(ec: str) -> str:
    ec = (ec or 'NONE').strip()
    return ELEMENT_CLASS_ALIASES.get(ec, ec)


def load_hits(path: Path | None) -> dict[tuple[str, str, str, str], dict[str, str]]:
    out: dict[tuple[str, str, str, str], dict[str, str]] = {}
    if path is None or not path.exists():
        return out
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            key = (
                row.get('query_contig', ''),
                row.get('pos', ''),
                row.get('end', ''),
                row.get('type', ''),
            )
            out[key] = row
    return out



def load_query_meta(path: Path | None) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    if path is None or not path.exists():
        return out
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            asm = row.get('query_asm')
            if asm:
                out[asm] = row
    return out



def load_ancestral(path: Path | None) -> dict[tuple[str, str], dict[str, object]]:
    summary: dict[tuple[str, str], dict[str, object]] = {}
    if path is None or not path.exists():
        return summary
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            key = (row.get('query_asm', ''), row.get('query_contig', ''))
            ent = summary.setdefault(key, {
                'clades': set(),
                'ranks': set(),
                'has_breakpoints': False,
                'segment_bp': 0,
                'variant_types': set(),
            })
            clade = row.get('clade', '')
            if clade and clade != '.':
                ent['clades'].add(clade)
            rank = row.get('clade_rank', '')
            if rank and rank != '.':
                ent['ranks'].add(rank)
            bp = row.get('breakpoints', '')
            if bp and bp != '.':
                ent['has_breakpoints'] = True
            try:
                ent['segment_bp'] = max(int(row.get('segment_bp', '0') or 0), ent['segment_bp'])
            except ValueError:
                pass
            vt = row.get('variant_type', '')
            if vt:
                ent['variant_types'].add(vt)
    return summary


def parse_float(value: str | None, default: float | None = None) -> float | None:
    if value is None or value == '':
        return default
    try:
        return float(value)
    except ValueError:
        return default


def load_expression(path: Path | None) -> dict[str, object]:
    exact: dict[tuple[str, str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    by_contig: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    if path is None or not path.exists():
        return {'exact': exact, 'by_contig': by_contig}

    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            qasm = row.get('query_asm', '')
            contig = row.get('query_contig', '')
            if not qasm or not contig:
                continue
            parsed = {
                'gene_id': row.get('gene_id', '.'),
                'gene_name': row.get('gene_name') or row.get('gene_id', '.'),
                'distance_bp': parse_float(row.get('distance_bp'), 10**9),
                'log2_fc': parse_float(row.get('log2_fc'), 0.0),
                'padj': parse_float(row.get('padj'), 1.0),
                'condition': row.get('condition', '.'),
            }
            pos = row.get('pos') or row.get('sv_pos') or ''
            end = row.get('end') or row.get('sv_end') or ''
            svtype = row.get('svtype') or row.get('type') or ''
            if pos and end and svtype:
                exact[(qasm, contig, pos, end, svtype)].append(parsed)
            by_contig[(qasm, contig)].append(parsed)
    return {'exact': exact, 'by_contig': by_contig}


def merge_expression_support(*datasets: dict[str, object]) -> dict[str, object]:
    exact: dict[tuple[str, str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    by_contig: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for ds in datasets:
        if not ds:
            continue
        for key, rows in ds.get('exact', {}).items():
            exact[key].extend(rows)
        for key, rows in ds.get('by_contig', {}).items():
            by_contig[key].extend(rows)
    return {'exact': exact, 'by_contig': by_contig}


def load_gene_annotations(path: Path | None) -> dict[tuple[str, str], list[dict[str, object]]]:
    genes: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    if path is None or not path.exists():
        return genes
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            qasm = row.get('query_asm') or row.get('asm') or '.'
            contig = row.get('query_contig') or row.get('contig') or row.get('chrom') or row.get('seqid') or ''
            gene_id = row.get('gene_id') or row.get('feature_id') or row.get('id') or ''
            if not contig or not gene_id:
                continue
            start = int(float(row.get('start') or row.get('gene_start') or row.get('pos') or 0))
            end = int(float(row.get('end') or row.get('gene_end') or row.get('stop') or start))
            if end < start:
                start, end = end, start
            genes[(qasm, contig)].append({
                'gene_id': gene_id,
                'gene_name': row.get('gene_name') or row.get('name') or gene_id,
                'start': start,
                'end': end,
            })
    for rows in genes.values():
        rows.sort(key=lambda r: (int(r['start']), int(r['end'])))
    return genes


def load_expression_long(path: Path | None) -> dict[tuple[str, str], dict[str, object]]:
    out: dict[tuple[str, str], dict[str, object]] = {}
    if path is None or not path.exists():
        return out
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            qasm = row.get('query_asm') or row.get('asm') or '.'
            gene_id = row.get('gene_id') or row.get('feature_id') or row.get('id') or ''
            group = row.get('condition') or row.get('group') or row.get('state') or row.get('contrast') or ''
            expr = None
            for field in ('expression', 'expr', 'tpm', 'fpkm', 'count', 'value'):
                expr = parse_float(row.get(field))
                if expr is not None:
                    break
            if not gene_id or not group or expr is None:
                continue
            key = (qasm, gene_id)
            ent = out.setdefault(key, {
                'gene_name': row.get('gene_name') or row.get('name') or gene_id,
                'groups': defaultdict(list),
            })
            ent['groups'][group].append(max(0.0, float(expr)))
    return out


def interval_distance(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    if end_a < start_b:
        return start_b - end_a
    if end_b < start_a:
        return start_a - end_b
    return 0


def log2_transform(values: list[float], pseudocount: float) -> list[float]:
    return [math.log2(max(0.0, v) + pseudocount) for v in values]


def mean_and_variance(values: list[float]) -> tuple[float, float]:
    if not values:
        return 0.0, 0.0
    mean = sum(values) / float(len(values))
    if len(values) < 2:
        return mean, 0.0
    var = sum((v - mean) ** 2 for v in values) / float(len(values) - 1)
    return mean, var


def welch_normal_pvalue(group_a: list[float], group_b: list[float]) -> float:
    mean_a, var_a = mean_and_variance(group_a)
    mean_b, var_b = mean_and_variance(group_b)
    denom = math.sqrt(var_a / max(1, len(group_a)) + var_b / max(1, len(group_b)))
    if denom == 0.0:
        return 1.0 if abs(mean_b - mean_a) < 1e-12 else 0.0
    z = abs(mean_b - mean_a) / denom
    return math.erfc(z / math.sqrt(2.0))


def quantify_gene_support(
    groups: dict[str, list[float]],
    gene_name: str,
    group_a: str | None,
    group_b: str | None,
    min_reps: int,
    pseudocount: float,
) -> dict[str, object] | None:
    eligible: dict[str, list[float]] = {
        name: log2_transform(vals, pseudocount)
        for name, vals in groups.items()
        if len(vals) >= min_reps
    }
    if len(eligible) < 2:
        return None

    chosen_a = group_a
    chosen_b = group_b
    if chosen_a and chosen_b:
        if chosen_a not in eligible or chosen_b not in eligible:
            return None
    else:
        names = list(eligible)
        if len(names) == 2:
            chosen_a, chosen_b = names[0], names[1]
        else:
            best_pair: tuple[str, str] | None = None
            best_delta = -1.0
            for i, name_a in enumerate(names):
                mean_a = mean_and_variance(eligible[name_a])[0]
                for name_b in names[i + 1:]:
                    mean_b = mean_and_variance(eligible[name_b])[0]
                    delta = abs(mean_b - mean_a)
                    if delta > best_delta:
                        best_delta = delta
                        best_pair = (name_a, name_b)
            if best_pair is None:
                return None
            chosen_a, chosen_b = best_pair

    vals_a = eligible[chosen_a]
    vals_b = eligible[chosen_b]
    mean_a = mean_and_variance(vals_a)[0]
    mean_b = mean_and_variance(vals_b)[0]
    if group_a is None and group_b is None and mean_b < mean_a:
        chosen_a, chosen_b = chosen_b, chosen_a
        vals_a, vals_b = vals_b, vals_a
        mean_a, mean_b = mean_b, mean_a

    return {
        'gene_name': gene_name,
        'log2_fc': mean_b - mean_a,
        'pvalue': welch_normal_pvalue(vals_a, vals_b),
        'condition': f'{chosen_b}_vs_{chosen_a}',
        'group_a': chosen_a,
        'group_b': chosen_b,
        'n_group_a': len(vals_a),
        'n_group_b': len(vals_b),
    }


def benjamini_hochberg(rows: list[dict[str, object]], pkey: str, outkey: str) -> None:
    if not rows:
        return
    ranked = sorted(
        [(idx, max(0.0, min(1.0, float(row.get(pkey, 1.0) or 1.0)))) for idx, row in enumerate(rows)],
        key=lambda item: item[1],
    )
    m = len(ranked)
    adjusted = [1.0] * m
    running = 1.0
    for rev_idx in range(m - 1, -1, -1):
        rank = rev_idx + 1
        _, pval = ranked[rev_idx]
        running = min(running, pval * m / rank)
        adjusted[rev_idx] = min(1.0, running)
    for adj, (row_idx, _) in zip(adjusted, ranked):
        rows[row_idx][outkey] = adj


def derive_expression_support_from_quant(
    records: list[dict[str, object]],
    hits: dict[tuple[str, str, str, str], dict[str, str]],
    gene_annotations: dict[tuple[str, str], list[dict[str, object]]],
    expression_long: dict[tuple[str, str], dict[str, object]],
    window_bp: int,
    group_a: str | None,
    group_b: str | None,
    min_reps: int,
    pseudocount: float,
    out_path: Path | None = None,
) -> dict[str, object]:
    exact: dict[tuple[str, str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    by_contig: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    quantified_rows: list[dict[str, object]] = []

    for rec in records:
        chrom = str(rec['chrom'])
        pos = str(rec['pos'])
        end = str(rec['end'])
        svtype = str(rec['svtype'])
        info = rec.get('info', {})
        hit = hits.get((chrom, pos, end, svtype), {})
        qasm = hit.get('query_asm') or info.get('QUERY_ASM', '.')
        if not qasm or qasm == '.':
            continue
        genes = gene_annotations.get((qasm, chrom)) or gene_annotations.get(('.', chrom)) or []
        if not genes:
            continue
        pos_int = int(rec['pos_int'])
        end_int = int(rec['end_int'])
        lo = pos_int - max(0, window_bp)
        hi = end_int + max(0, window_bp)
        candidate_key = (qasm, chrom, pos, end, svtype)
        for gene in genes:
            gstart = int(gene['start'])
            gend = int(gene['end'])
            if gend < lo:
                continue
            if gstart > hi:
                break
            gene_id = str(gene['gene_id'])
            expr = expression_long.get((qasm, gene_id)) or expression_long.get(('.', gene_id))
            if not expr:
                continue
            quant = quantify_gene_support(
                expr.get('groups', {}),
                str(gene.get('gene_name', gene_id)),
                group_a,
                group_b,
                min_reps,
                pseudocount,
            )
            if quant is None:
                continue
            quantified_rows.append({
                'query_asm': qasm,
                'query_contig': chrom,
                'pos': pos,
                'end': end,
                'svtype': svtype,
                'gene_id': gene_id,
                'gene_name': quant['gene_name'],
                'distance_bp': interval_distance(pos_int, end_int, gstart, gend),
                'log2_fc': quant['log2_fc'],
                'pvalue': quant['pvalue'],
                'condition': quant['condition'],
                'group_a': quant['group_a'],
                'group_b': quant['group_b'],
                'n_group_a': quant['n_group_a'],
                'n_group_b': quant['n_group_b'],
                'candidate_key': candidate_key,
            })

    benjamini_hochberg(quantified_rows, 'pvalue', 'padj')

    for row in quantified_rows:
        parsed = {
            'gene_id': row['gene_id'],
            'gene_name': row['gene_name'],
            'distance_bp': row['distance_bp'],
            'log2_fc': row['log2_fc'],
            'padj': row.get('padj', 1.0),
            'condition': row['condition'],
        }
        exact[row['candidate_key']].append(parsed)
        by_contig[(row['query_asm'], row['query_contig'])].append(parsed)

    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open('w', newline='') as fh:
            fieldnames = [
                'query_asm', 'query_contig', 'pos', 'end', 'svtype',
                'gene_id', 'gene_name', 'distance_bp', 'log2_fc', 'pvalue', 'padj',
                'condition', 'group_a', 'group_b', 'n_group_a', 'n_group_b',
            ]
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for row in quantified_rows:
                writer.writerow({k: row.get(k, '.') for k in fieldnames})

    return {'exact': exact, 'by_contig': by_contig}


def summarize_expression(records: list[dict[str, object]]) -> dict[str, object] | None:
    if not records:
        return None
    ordered = sorted(
        records,
        key=lambda r: (
            0 if (r.get('padj') is not None and float(r['padj']) <= 0.05) else 1,
            -(abs(float(r.get('log2_fc') or 0.0))),
            float(r.get('distance_bp') or 10**9),
        ),
    )
    sig = [
        r for r in records
        if float(r.get('padj') or 1.0) <= 0.05 and abs(float(r.get('log2_fc') or 0.0)) >= 1.0
    ]
    best = ordered[0]
    return {
        'supported': bool(sig),
        'best_gene': best.get('gene_name') or best.get('gene_id') or '.',
        'distance_bp': int(float(best.get('distance_bp') or 10**9)) if best.get('distance_bp') is not None else 10**9,
        'log2_fc': float(best.get('log2_fc') or 0.0),
        'padj': float(best.get('padj') or 1.0),
        'condition': best.get('condition', '.'),
        'significant_genes': len(sig),
        'max_abs_log2_fc': max(abs(float(r.get('log2_fc') or 0.0)) for r in records),
    }


def expression_for_candidate(
    expr: dict[str, object],
    qasm: str,
    contig: str,
    pos: str,
    end: str,
    svtype: str,
) -> dict[str, object] | None:
    exact = expr.get('exact', {})
    by_contig = expr.get('by_contig', {})
    if (qasm, contig, pos, end, svtype) in exact:
        return summarize_expression(exact[(qasm, contig, pos, end, svtype)])
    return summarize_expression(by_contig.get((qasm, contig), []))



def is_hgt_candidate(ec: str, svtype: str, annot: str) -> bool:
    """True when the call has hallmarks of horizontal gene transfer:
    HGT element class OR a translocation with cross-clade novelty signal."""
    if ec == 'HGT':
        return True
    if svtype == 'TRA' and annot in {'NOVEL', 'NOVEL_WEAK'}:
        return True
    if svtype == 'OFF_REF' and annot in {'NOVEL', 'NOVEL_WEAK'} and ec in MGE_INTEGRATIVE:
        return True
    return False


def mge_subtype(ec: str) -> str:
    """Coarser MGE category for reporting: integrative, transposable, or repeat."""
    if ec in MGE_INTEGRATIVE:
        return 'integrative'
    if ec in MGE_TRANSPOSABLE:
        return 'transposable'
    if ec in MGE_REPEAT_BASED:
        return 'repeat'
    return 'none'


def classify_candidate(
    svtype: str,
    ec: str,
    annot: str,
    anc: dict[str, object] | None,
    expr: dict[str, object] | None,
) -> tuple[str, int, str]:
    novelty = annot in NOVEL_TIERS or svtype == 'OFF_REF'
    te_like = ec in TE_CLASSES and ec != 'NONE'
    hgt = is_hgt_candidate(ec, svtype, annot)
    ancestry_switch = bool(anc and (len(anc.get('clades', set())) > 1 or anc.get('has_breakpoints')))
    expr_supported = bool(expr and expr.get('supported'))

    # HGT candidates rank highest: cross-clade or Starship/integrative-island events
    # that are novel relative to the same-clade references.
    if hgt and ancestry_switch:
        return 'hgt_candidate', 7, 'Cross-clade TRA or HGT-class insertion with multi-clade ancestry — strong HGT signal.'
    if hgt:
        return 'hgt_candidate', 6, 'HGT element class or cross-clade translocation with novel sequence tier.'
    if expr_supported and te_like and ancestry_switch:
        return 'mosaic_te_expression_link', 6, 'TE-associated ancestry switch near a gene with direct expression support.'
    if expr_supported and te_like:
        return 'te_expression_link', 6, 'TE-associated structural event near a significantly shifted gene.'
    if te_like and novelty:
        return 'novel_te_architecture', 5, 'Novel or strongly diverged sequence with TE-like content.'
    if ancestry_switch and te_like:
        return 'mosaic_te_lineage_switch', 4, 'TE-associated event with multi-clade or breakpoint-resolved ancestry.'
    if ancestry_switch:
        return 'mosaic_lineage_switch', 3, 'Contig shows mixed ancestry or recombination-style clade switches.'
    if te_like and svtype in {'INS', 'DUP', 'INV', 'TRA', 'DEL'}:
        return 'te_architecture_rewiring', 3, 'TE-like structural event may rewire genome architecture.'
    if novelty:
        return 'sequence_novelty', 2, 'Sequence appears novel or strongly diverged relative to indexed references.'
    return 'other', 1, 'Interesting structural event, but without clear novelty or TE evidence.'



def choose_functional_example(
    candidate_type: str,
    svtype: str,
    ec: str,
    scenario: str,
    anc: dict[str, object] | None,
    expr: dict[str, object] | None,
) -> dict[str, str]:
    best_score = -1
    best: dict[str, object] | None = None
    for example in FUNCTIONAL_EXAMPLES:
        score = int(example.get('priority', 0))
        match_ec = set(example.get('match_ec', set()))
        if match_ec:
            if ec not in match_ec:
                continue
            score += 3
        match_ct = set(example.get('match_candidate_types', set()))
        if match_ct:
            if candidate_type not in match_ct:
                continue
            score += 3
        match_sv = set(example.get('match_svtypes', set()))
        if match_sv:
            if svtype not in match_sv:
                continue
            score += 2
        match_scenarios = set(example.get('match_scenarios', set()))
        if match_scenarios:
            if scenario not in match_scenarios:
                continue
            score += 2
        if example.get('require_expression_support') and not (expr and expr.get('supported')):
            continue
        if example.get('require_ancestral_breakpoints') and not (anc and anc.get('has_breakpoints')):
            continue
        if anc and anc.get('has_breakpoints'):
            score += 1
        if expr and expr.get('supported'):
            score += 2
        if score > best_score:
            best_score = score
            best = example

    if best is None:
        return {
            'functional_example': 'Nearby-gene regulatory follow-up',
            'evidence_axis': 'expression_screen',
            'example_system': 'generic TE-near-gene cis-regulatory screen',
            'real_data_signal': 'Novel sequence or structural change may alter nearby transcription even when the specific mechanism is not yet known.',
            'functional_hypothesis': 'Test whether genes near the breakpoint or insertion shift expression across relevant conditions.',
            'suggested_assay': 'RNA-seq or qPCR on genes flanking the locus, then inspect local copy number and TE context.',
        }

    return {
        'functional_example': str(best['name']),
        'evidence_axis': str(best['evidence_axis']),
        'example_system': str(best['system']),
        'real_data_signal': str(best['real_data_signal']),
        'functional_hypothesis': str(best['why_relevant']),
        'suggested_assay': str(best['suggested_readout']),
    }



def main() -> int:
    ap = argparse.ArgumentParser(description='Summarize new-biology candidate signals from fungi_graphsv_tol outputs.')
    ap.add_argument('--vcf', type=Path, required=True)
    ap.add_argument('--hits', type=Path)
    ap.add_argument('--ancestral', type=Path)
    ap.add_argument('--expression-tsv', type=Path)
    ap.add_argument('--expression-long-tsv', type=Path,
                    help='Sample-level long-form expression table with query_asm, gene_id, expression, and condition/group columns.')
    ap.add_argument('--gene-annotations', type=Path,
                    help='Gene coordinate TSV with query_asm, query_contig, gene_id, gene_name, start, end.')
    ap.add_argument('--derived-expression-out', type=Path,
                    help='Optional TSV path to write candidate-level expression quantification derived from --expression-long-tsv.')
    ap.add_argument('--expression-window-bp', type=int, default=5000)
    ap.add_argument('--expression-group-a')
    ap.add_argument('--expression-group-b')
    ap.add_argument('--expression-min-reps', type=int, default=2)
    ap.add_argument('--expression-pseudocount', type=float, default=1.0)
    ap.add_argument('--query-metadata', type=Path)
    ap.add_argument('--out-tsv', type=Path, required=True)
    ap.add_argument('--summary-json', type=Path, required=True)
    ap.add_argument('--phylum', default='unknown')
    ap.add_argument('--top-n', type=int, default=50)
    args = ap.parse_args()

    if bool(args.expression_long_tsv) != bool(args.gene_annotations):
        ap.error('--expression-long-tsv and --gene-annotations must be provided together')

    hits = load_hits(args.hits)
    meta = load_query_meta(args.query_metadata)
    ancestral = load_ancestral(args.ancestral)
    records = load_vcf_records(args.vcf)

    expression = load_expression(args.expression_tsv)
    if args.expression_long_tsv and args.gene_annotations:
        derived_expression = derive_expression_support_from_quant(
            records,
            hits,
            load_gene_annotations(args.gene_annotations),
            load_expression_long(args.expression_long_tsv),
            args.expression_window_bp,
            args.expression_group_a,
            args.expression_group_b,
            args.expression_min_reps,
            args.expression_pseudocount,
            args.derived_expression_out,
        )
        expression = merge_expression_support(expression, derived_expression)

    rows: list[dict[str, object]] = []
    for rec in records:
        chrom = str(rec['chrom'])
        pos = str(rec['pos'])
        end = str(rec['end'])
        info = rec['info']
        svtype = str(rec['svtype'])
        raw_ec = info.get('EC', 'NONE')
        ec = normalize_element_class(raw_ec)
        annot = info.get('ANNOT', '.')
        hit = hits.get((chrom, pos, end, svtype), {})
        qasm = hit.get('query_asm') or info.get('QUERY_ASM', '.')
        scenario = meta.get(qasm, {}).get('scenario', '.') if qasm in meta else '.'
        architecture = meta.get(qasm, {}).get('architecture', '.') if qasm in meta else '.'
        lifestyle = meta.get(qasm, {}).get('lifestyle', '.') if qasm in meta else '.'
        anc = ancestral.get((qasm, chrom)) or ancestral.get(('.', chrom))
        expr = expression_for_candidate(expression, qasm, chrom, pos, end, svtype)
        candidate_type, priority, rationale = classify_candidate(svtype, ec, annot, anc, expr)
        if candidate_type == 'other':
            continue
        clades = ','.join(sorted(anc['clades'])) if anc else '.'
        clade_ranks = ','.join(sorted(anc['ranks'])) if anc else '.'
        example = choose_functional_example(candidate_type, svtype, ec, scenario, anc, expr)
        expr_supported = 'yes' if expr and expr.get('supported') else 'no'
        expr_gene = expr.get('best_gene', '.') if expr else '.'
        expr_distance = expr.get('distance_bp', '.') if expr else '.'
        expr_log2_fc = expr.get('log2_fc', '.') if expr else '.'
        expr_padj = expr.get('padj', '.') if expr else '.'
        expr_condition = expr.get('condition', '.') if expr else '.'
        rows.append({
            'priority': priority,
            'candidate_type': candidate_type,
            'phylum': args.phylum,
            'query_asm': qasm,
            'query_contig': chrom,
            'scenario': scenario,
            'lifestyle': lifestyle,
            'architecture': architecture,
            'svtype': svtype,
            'element_class': raw_ec,
            'mge_subtype': mge_subtype(ec),
            'hgt_flag': 'yes' if is_hgt_candidate(ec, svtype, annot) else 'no',
            'novelty': annot,
            'pos': int(pos),
            'end': int(end),
            'ref_target': hit.get('ref_asm', info.get('CL', '.')),
            'alignment_mode': hit.get('alignment_mode', info.get('ALIGNMENT_MODE', '.')),
            'ancestral_clades': clades,
            'ancestral_ranks': clade_ranks,
            'ancestral_breakpoints': 'yes' if anc and anc.get('has_breakpoints') else 'no',
            'ancestral_segment_bp': anc.get('segment_bp', 0) if anc else 0,
            'expression_supported': expr_supported,
            'expression_gene': expr_gene,
            'expression_distance_bp': expr_distance,
            'expression_log2_fc': expr_log2_fc,
            'expression_padj': expr_padj,
            'expression_condition': expr_condition,
            'rationale': rationale,
            'functional_example': example['functional_example'],
            'evidence_axis': example['evidence_axis'],
            'example_system': example['example_system'],
            'real_data_signal': example['real_data_signal'],
            'functional_hypothesis': example['functional_hypothesis'],
            'suggested_assay': example['suggested_assay'],
            'follow_up': 'Prioritize for long-read validation, locus inspection, and RNA-seq/qPCR if a nearby gene is biologically relevant.',
        })

    rows.sort(key=lambda r: (-int(r['priority']), str(r['candidate_type']), str(r['query_asm']), str(r['query_contig']), int(r['pos'])))
    rows = rows[: max(1, args.top_n)]

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_tsv.open('w', newline='') as fh:
        fieldnames = [
            'priority', 'candidate_type', 'phylum', 'query_asm', 'query_contig', 'scenario', 'lifestyle',
            'architecture', 'svtype', 'element_class', 'mge_subtype', 'hgt_flag', 'novelty', 'pos', 'end',
            'ref_target', 'alignment_mode', 'ancestral_clades', 'ancestral_ranks', 'ancestral_breakpoints',
            'ancestral_segment_bp', 'expression_supported', 'expression_gene', 'expression_distance_bp',
            'expression_log2_fc', 'expression_padj', 'expression_condition', 'rationale',
            'functional_example', 'evidence_axis', 'example_system',
            'real_data_signal', 'functional_hypothesis', 'suggested_assay', 'follow_up'
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

    counts = Counter(r['candidate_type'] for r in rows)
    by_sv = Counter(r['svtype'] for r in rows)
    by_ec = Counter(r['element_class'] for r in rows)
    by_axis = Counter(r['evidence_axis'] for r in rows)
    by_expr = Counter(r['expression_supported'] for r in rows)
    example_names = Counter(r['functional_example'] for r in rows)

    # SV phylogeny: count SVs per phylum / scenario / lifestyle for phylogenetic landscape.
    phylo_dist: dict[str, dict[str, int]] = {}
    for r in rows:
        ph = str(r.get('phylum', '.') or '.')
        sc = str(r.get('scenario', '.') or '.')
        sv = str(r.get('svtype', '.') or '.')
        phylo_dist.setdefault(ph, {})
        phylo_dist[ph][sv] = phylo_dist[ph].get(sv, 0) + 1
        phylo_dist[ph].setdefault('_scenario', sc)

    # MGE breakdown: separate integrative islands (HGT/Starship), transposable elements,
    # and repeat-based elements for downstream MGE-specific reporting.
    mge_breakdown: dict[str, int] = {'integrative': 0, 'transposable': 0, 'repeat': 0, 'none': 0}
    for r in rows:
        mge_breakdown[mge_subtype(str(r.get('element_class', 'NONE')))] += 1

    # HGT-specific summary: candidates classified as hgt_candidate with TRA or OFF_REF type.
    hgt_rows = [r for r in rows if r.get('candidate_type') == 'hgt_candidate']
    hgt_by_svtype = dict(Counter(str(r.get('svtype')) for r in hgt_rows))
    hgt_by_phylum = dict(Counter(str(r.get('phylum')) for r in hgt_rows))

    with args.summary_json.open('w') as fh:
        json.dump({
            'phylum': args.phylum,
            'candidate_count': len(rows),
            'by_candidate_type': dict(counts),
            'by_svtype': dict(by_sv),
            'by_element_class': dict(by_ec),
            'by_evidence_axis': dict(by_axis),
            'by_expression_support': dict(by_expr),
            'functional_examples': dict(example_names),
            'phylo_sv_distribution': phylo_dist,
            'mge_breakdown': mge_breakdown,
            'hgt_summary': {
                'count': len(hgt_rows),
                'by_svtype': hgt_by_svtype,
                'by_phylum': hgt_by_phylum,
            },
            'top_priorities': rows[:10],
        }, fh, indent=2, sort_keys=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
