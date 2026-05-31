#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import re
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

# Pezizomycotina classes - the only fungal subphylum where RIP (Repeat-Induced
# Point mutation) operates. Outside this set the sequence-only RIP detector in
# layer1_clade_graph.hpp can still fire on AT-rich / TpA-enriched repeats, but
# the RIP label is biologically wrong (Saccharomycotina yeasts, Glomeromycota
# AMF, Basidiomycota, early-diverging fungi lack the meiotic RIP machinery).
RIP_SUPPORTING_CLASSES = {
    'sordariomycetes', 'eurotiomycetes', 'dothideomycetes', 'leotiomycetes',
    'pezizomycetes', 'orbiliomycetes', 'lecanoromycetes', 'xylonomycetes',
    'lichinomycetes', 'coniocybomycetes', 'geoglossomycetes', 'arthoniomycetes',
}


def rip_supported_by_clade(meta_row: dict[str, str]) -> bool:
    """True if the host clade is in Pezizomycotina (RIP-supporting).

    Returns True when class is unknown so we don't silently drop labels on
    panels without taxonomic metadata; the per-shard caller is responsible
    for supplying class via query_manifest.tsv.
    """
    clazz = (meta_row.get('class') or '').strip().lower()
    if not clazz or clazz == '.':
        return True
    return clazz in RIP_SUPPORTING_CLASSES


ASM_EXT_STRIP = (
    ".fasta.gz", ".fastq.gz", ".fna.gz", ".fa.gz", ".fq.gz",
    ".fasta", ".fastq", ".fna", ".fa", ".fq",
)

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
        'real_data_signal': 'Translocation breakpoints with HGT-class sequence and low same-clade overlap are consistent with candidate donor-recipient boundaries of a horizontally transferred genomic island.',
        'why_relevant': 'Highest-priority screening hit when a TRA or off-reference segment is both phylogenetically novel within its clade and carried in a sequence class (HGT / Starship) associated with lateral transfer.',
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



def query_asm_aliases(asm: str) -> list[str]:
    out = [asm]
    seen = {asm}
    lower = asm.lower()
    for ext in ASM_EXT_STRIP:
        if lower.endswith(ext):
            stripped = asm[: -len(ext)]
            if stripped and stripped not in seen:
                out.append(stripped)
                seen.add(stripped)
    downsample_stripped = re.sub(r"\.\d+$", "", asm)
    if downsample_stripped and downsample_stripped not in seen:
        out.append(downsample_stripped)
        seen.add(downsample_stripped)
    acc_match = re.search(r"(GC[AF]_\d+)\.(\d+)", asm)
    if acc_match:
        accession_alias = f"{acc_match.group(1)}_{acc_match.group(2)}"
        if accession_alias not in seen:
            out.append(accession_alias)
            seen.add(accession_alias)
    base = Path(asm).name
    if base and base != asm:
        for alias in query_asm_aliases(base):
            if alias not in seen:
                out.append(alias)
                seen.add(alias)
    return out


def load_ecological_traits(path: Path | None) -> dict[str, dict[str, str]]:
    traits: dict[str, dict[str, str]] = {}
    if path is None or not path.exists():
        return traits
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            qasm = row.get('query_asm') or row.get('asm') or ''
            if not qasm:
                continue
            for alias in query_asm_aliases(qasm):
                traits.setdefault(alias, row)
    return traits


def load_fungaltraits_csv(path: Path | None) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    by_species: dict[str, dict[str, str]] = {}
    by_genus: dict[str, dict[str, str]] = {}
    if path is None or not path.exists():
        return by_species, by_genus
    with path.open(encoding='utf-8', errors='replace') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            species = (row.get('speciesMatched') or row.get('species') or row.get('SPECIES') or '').replace('_', ' ').strip().lower()
            genus = (row.get('GENUS') or row.get('Genus') or '').strip().lower()
            if not genus and species:
                genus = species.split()[0]
            trait_name = (row.get('trait_name') or '').strip().lower()
            value = (row.get('value') or '').strip()
            if trait_name and value:
                slim = {'species': species}
                if trait_name == 'trophic_mode_fg':
                    slim['primary_lifestyle'] = value
                    slim['trophic_mode'] = value
                elif trait_name == 'substrate':
                    slim['substrate_or_host'] = value
                elif trait_name in {'growth_form', 'fruitbody_type', 'guild'}:
                    slim[trait_name] = value
                else:
                    continue
                row = slim
            if species:
                by_species.setdefault(species, {}).update({k: v for k, v in row.items() if v})
            if genus:
                by_genus.setdefault(genus, {}).update({k: v for k, v in row.items() if v})
    return by_species, by_genus


def ecological_context(
    qasm: str,
    meta_row: dict[str, str],
    traits_by_qasm: dict[str, dict[str, str]],
    fungal_by_species: dict[str, dict[str, str]],
    fungal_by_genus: dict[str, dict[str, str]],
) -> dict[str, str]:
    def pick(*values: str | None) -> str:
        for value in values:
            if value not in (None, '', '.'):
                return str(value)
        return '.'

    def inferred_trait() -> str:
        species_text = (meta_row.get('species') or '').lower()
        genus = (meta_row.get('genus') or '').lower()
        order = (meta_row.get('order') or '').lower()
        family = (meta_row.get('family') or '').lower()
        phylum = (meta_row.get('phylum') or '').lower()
        clazz = (meta_row.get('class') or '').lower()
        scenario = (meta_row.get('scenario') or '').lower()
        # All Glomeromycota / Glomeromycotina are obligate arbuscular mycorrhizal
        # symbionts. The phylum/class check ensures Diversispora, Glomus,
        # Funneliformis, Claroideoglomus, Septoglomus, Acaulospora, Scutellospora,
        # Paraglomus, Archaeospora, etc. are not silently misclassified as
        # filamentous saprotrophs when FungalTraits has no species-level row.
        amf_genera = {
            'rhizophagus', 'gigaspora', 'glomus', 'diversispora',
            'funneliformis', 'claroideoglomus', 'septoglomus',
            'acaulospora', 'scutellospora', 'paraglomus', 'archaeospora',
            'racocetra', 'cetraspora', 'dentiscutata', 'pacispora',
            'redeckera', 'sclerocystis', 'entrophospora', 'ambispora',
        }
        amf_higher = {
            'glomeromycota', 'glomeromycotina', 'glomeromycetes',
            'glomerales', 'diversisporales', 'archaeosporales', 'paraglomerales',
            'glomeraceae', 'diversisporaceae', 'gigasporaceae',
            'claroideoglomeraceae', 'acaulosporaceae', 'paraglomeraceae',
            'archaeosporaceae', 'ambisporaceae',
        }
        if (
            'arbuscular' in scenario
            or genus in amf_genera
            or phylum in amf_higher
            or clazz in amf_higher
            or order in amf_higher
            or family in amf_higher
        ):
            return 'Symbiotroph_arbuscular_mycorrhizal'
        if any(x in species_text or x == genus for x in (
            'puccinia', 'ustilago', 'mycosarcoma', 'pyricularia', 'fusarium',
            'leptosphaeria', 'zymoseptoria'
        )):
            return 'Pathotroph_plant_pathogen'
        if genus in {'candida', 'cryptococcus', 'nakaseomyces'} or 'candida' in species_text:
            return 'Yeast_opportunistic_pathogen'
        if genus in {
            'saccharomyces', 'lachancea', 'kluyveromyces',
            'ogataea', 'pichia', 'komagataella', 'yarrowia',
            'hanseniaspora', 'eremothecium', 'torulaspora', 'zygosaccharomyces',
        }:
            return 'Saprotroph_yeast'
        if order in {'polyporales', 'agaricales', 'russulales'}:
            return 'Saprotroph_wood_decay'
        if 'pathogen' in scenario:
            return 'Pathotroph'
        return '.'

    trait_row: dict[str, str] = {}
    for alias in query_asm_aliases(qasm):
        if alias in traits_by_qasm:
            trait_row = traits_by_qasm[alias]
            break
    species = pick(trait_row.get('species'), meta_row.get('species'))
    species_low = species.replace('_', ' ').strip().lower()
    fungal_row = (
        fungal_by_species.get(species_low)
        or fungal_by_genus.get(species_low.split()[0] if species_low else '')
        or {}
    )
    primary = pick(
        trait_row.get('primary_lifestyle'),
        trait_row.get('trophic_mode'),
        fungal_row.get('primary_lifestyle'),
        fungal_row.get('trophic_mode'),
        meta_row.get('lifestyle'),
        inferred_trait(),
    )
    secondary = pick(trait_row.get('secondary_lifestyle'), fungal_row.get('secondary_lifestyle'))
    trophic = pick(trait_row.get('trophic_mode'), fungal_row.get('trophic_mode'), primary)
    substrate = pick(trait_row.get('substrate_or_host'), fungal_row.get('substrate_or_host'), fungal_row.get('Substrate'))
    return {
        'species': species,
        'ecological_trait': primary,
        'secondary_lifestyle': secondary,
        'trophic_mode': trophic,
        'substrate_or_host': substrate,
    }


def load_query_meta(path: Path | None) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    if path is None or not path.exists():
        return out
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            asm = row.get('query_asm')
            if asm:
                for alias in query_asm_aliases(asm):
                    out.setdefault(alias, row)
            for path_field in ('path', 'benchmark_ref_fasta'):
                path_val = row.get(path_field) or ''
                if path_val:
                    for alias in query_asm_aliases(Path(path_val).name):
                        out.setdefault(alias, row)
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
                'strand': row.get('strand') or '.',
                'biotype': row.get('biotype') or '.',
                'product': row.get('product') or '.',
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


def nearest_gene_for_candidate(
    gene_annotations: dict[tuple[str, str], list[dict[str, object]]],
    qasm: str,
    contig: str,
    pos: str,
    end: str,
) -> dict[str, object] | None:
    """Return {'gene_id', 'gene_name', 'distance_bp'} for the nearest annotated
    gene to the (pos, end) breakpoint on (qasm, contig). Falls back to the
    wildcard ('.', contig) bucket if no asm-keyed match exists. Returns None
    when no gene_annotations were loaded for this contig.

    This is the "no expression matrix" fallback: it lets biology_candidates.tsv
    show *which* gene a breakpoint is sitting on or near, even when nobody has
    measured its expression in a public DB. Without this, expression_gene was
    always '.' and the analyzer's follow-up suggestion stayed generic.
    """
    genes = gene_annotations.get((qasm, contig)) or gene_annotations.get(('.', contig))
    if not genes:
        return None
    try:
        sv_lo = int(pos)
        sv_hi = int(end)
    except ValueError:
        return None
    if sv_hi < sv_lo:
        sv_lo, sv_hi = sv_hi, sv_lo
    best: dict[str, object] | None = None
    best_dist: int | None = None
    for gene in genes:
        gstart = int(gene['start'])
        gend = int(gene['end'])
        if gend < gstart:
            gstart, gend = gend, gstart
        if gend < sv_lo:
            dist = sv_lo - gend
        elif gstart > sv_hi:
            dist = gstart - sv_hi
        else:
            dist = 0
        if best_dist is None or dist < best_dist:
            best_dist = dist
            best = {
                'gene_id': str(gene['gene_id']),
                'gene_name': str(gene.get('gene_name') or gene['gene_id']),
                'distance_bp': dist,
                'product': str(gene.get('product') or '.'),
                'biotype': str(gene.get('biotype') or '.'),
            }
            if dist == 0:
                break
    return best


def nearest_gene_for_locus(
    gene_annotations: dict[tuple[str, str], list[dict[str, object]]],
    qasm_candidates: list[str],
    contig: str,
    pos: str,
    end: str,
) -> dict[str, object] | None:
    for qasm in qasm_candidates:
        if not qasm or qasm == '.':
            continue
        for alias in query_asm_aliases(qasm):
            hit = nearest_gene_for_candidate(gene_annotations, alias, contig, pos, end)
            if hit is not None:
                return hit
    return nearest_gene_for_candidate(gene_annotations, '.', contig, pos, end)



def is_hgt_candidate(ec: str, svtype: str, annot: str) -> bool:
    """True when the call has candidate hallmarks of horizontal gene transfer:
    HGT-like sequence OR a translocation with cross-clade novelty signal."""
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
        return 'hgt_candidate', 7, 'Cross-clade TRA or HGT-class insertion with multi-clade ancestry - strong candidate HGT signal.'
    if hgt:
        return 'hgt_candidate', 6, 'HGT-like element class or cross-clade translocation with novel sequence tier.'
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
    if svtype in {'INS', 'DEL', 'DUP', 'INV', 'TRA'}:
        return 'structural_sv_signal', 1, 'Breakpoint-resolved structural variant with potential genome architecture impact.'
    return 'other', 1, 'Interesting structural event, but without clear novelty or TE evidence.'


def select_diverse_rows(rows: list[dict[str, object]], top_n: int) -> list[dict[str, object]]:
    """Return a ranked top-N while preserving representation across SV classes.

    Fungal panels with many novel/off-reference unitigs can otherwise fill the
    entire top-N with OFF_REF candidates and hide biologically relevant DEL,
    INS, DUP, INV, or TRA rows.  Keep the global priority order, but reserve a
    small quota for each observed SV type before filling the remaining slots.
    """
    if top_n <= 0 or len(rows) <= top_n:
        return rows

    svtypes = sorted({str(r.get('svtype', '.')) for r in rows if r.get('svtype') not in {None, '', '.'}})
    if len(svtypes) <= 1:
        return rows[:top_n]

    quota = max(3, min(25, top_n // (2 * len(svtypes))))
    selected: list[dict[str, object]] = []
    selected_ids: set[int] = set()

    for svtype in svtypes:
        taken = 0
        for row in rows:
            if id(row) in selected_ids or str(row.get('svtype', '.')) != svtype:
                continue
            selected.append(row)
            selected_ids.add(id(row))
            taken += 1
            if taken >= quota or len(selected) >= top_n:
                break
        if len(selected) >= top_n:
            break

    for row in rows:
        if len(selected) >= top_n:
            break
        if id(row) not in selected_ids:
            selected.append(row)
            selected_ids.add(id(row))

    selected.sort(key=lambda r: (-int(r['priority']), str(r['candidate_type']), str(r['query_asm']), str(r['query_contig']), int(r['pos'])))
    return selected



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



# --- Cross-guild aggregation -------------------------------------------------
# These functions read the per-shard biology_findings.tsv files produced by
# the per-query analyzer and synthesize a guild-level view (AMF / EMF / yeast /
# saprotrophic Pezizomycotina / endophyte / mushroom / truffle / etc.). They
# explicitly stratify by *pangenome scope* - single-reference-equivalent vs
# pangenome-only - so the "what would have been missed without a pangenome"
# fraction is computed per guild, per candidate axis, and tied back to ecology
# (lifestyle / ecological_trait / substrate_or_host) and gene-expression links
# (expression_supported / affected_gene). Triggered by --cross-guild-shard-dir;
# left out of the per-shard code path so existing single-shard invocations
# stay unchanged.

# Genus -> guild mapping. Same source-of-truth used by the augment + 165 panel
# helpers. Genera that legitimately span guilds (e.g. Serendipita = orchid-MF
# *or* endophyte; Lactarius = EMF) are mapped to their most common guild here;
# fine-grained re-classification can be done downstream by reading lifestyle.
GUILD_BY_GENUS = {
    # AMF (Glomeromycotina)
    **{g: 'AMF' for g in (
        'Rhizophagus','Rhizoglomus','Diversispora','Acaulospora','Funneliformis',
        'Gigaspora','Glomus','Paraglomus','Scutellospora','Claroideoglomus',
        'Entrophospora','Ambispora','Septoglomus','Sieverdingia','Simiglomus',
        'Archaeospora','Pacispora','Racocetra','Dentiscutata','Cetraspora',
        'Albahypha','Dominikia','Blaszkowskia','Sclerocystis','Redeckera',
        'Pervetustus',
    )},
    # EMF (Basidiomycota mostly, plus a few Ascomycota)
    **{g: 'EMF' for g in (
        'Suillus','Pisolithus','Laccaria','Cenococcum','Paxillus','Trichophaea',
        'Sphaerosporella','Geastrum','Scleroderma','Chroogomphus','Gomphidius',
        'Hebeloma','Helvella','Amanita','Boletus','Lactarius','Xerocomus',
        'Tylopilus','Leccinum','Clavulina','Austroboletus','Sparassis','Hydnum',
        'Cantharellus','Rhizopogon','Russula','Thelephora','Cortinarius',
        'Hydnellum','Gyroporus','Geopora','Hygrophorus','Tricholoma','Lanmaoa',
        'Butyriboletus','Pseudosperma','Inocybe','Suillellus',
    )},
    # Yeast (model)
    **{g: 'Yeast' for g in ('Saccharomyces','Yarrowia')},
    # Edible mushrooms
    **{g: 'Mushroom' for g in ('Pleurotus','Lentinula','Agaricus')},
    # Truffles
    **{g: 'Truffle' for g in ('Tuber',)},
    # Endophytes (incl. orchid- and ericoid-MF that aren't single-guild)
    **{g: 'Endophyte' for g in (
        'Serendipita','Piriformospora','Beauveria','Paecilomyces',
        'Colletotrichum','Oidiodendron','Meliniomyces','Rhizoscyphus',
        'Tulasnella','Alternaria',
    )},
    # Filamentous Ascomycota saprotrophs / pathogens
    **{g: 'Filamentous' for g in (
        'Aspergillus','Penicillium','Talaromyces','Fusarium','Trichoderma',
        'Neurospora','Monascus','Cladosporium','Acremonium','Rhizoctonia',
        'Rhizopus','Mucor','Eremothecium','Ustilago','Mycosarcoma','Candida',
        'Blumeria',
    )},
}


def classify_guild(genus: str, species: str = '', lifestyle: str = '') -> str:
    """Return the guild label for a species. Defaults to 'Other' so callers
    can surface unmapped genera without dropping their rows.

    Kept for backward compatibility with code that doesn't have phylum/class
    in hand. New callers should prefer classify_guild_taxonomic() which uses
    phylum + class for a phylogeny-aware 4-way bucketing
    (AMF / Filamentous / Basidio / Yeast).
    """
    if genus:
        g = genus.split()[0].lstrip('[').rstrip(']')
        if g in GUILD_BY_GENUS:
            return GUILD_BY_GENUS[g]
    if species:
        first = species.split()[0].lstrip('[').rstrip(']')
        if first in GUILD_BY_GENUS:
            return GUILD_BY_GENUS[first]
    if lifestyle and lifestyle != '.':
        ll = lifestyle.lower()
        if 'mycorrhiz' in ll:
            return 'AMF' if 'arbusc' in ll else 'EMF'
        if 'endophyt' in ll:
            return 'Endophyte'
        if 'sapro' in ll:
            return 'Filamentous'
        if 'pathotroph' in ll or 'patho' in ll:
            return 'Filamentous'
    return 'Other'


# Class-based 4-way taxonomy used by the cross-guild biology figure (manuscript
# Fig 2A). Resolution gates on subphylum/class, not phylum alone, because
# Ascomycota contains both Pezizomycotina (filamentous, RIP-active) and
# Saccharomycotina (yeasts, no RIP) - collapsing them loses the architectural
# distinction the volcano is designed to surface.
YEAST_CLASSES = {
    'saccharomycetes', 'pichiomycetes', 'dipodascomycetes',
    'schizosaccharomycetes', 'taphrinomycetes',
}
# Pezizomycotina classes - same set used by RIP_SUPPORTING_CLASSES, kept
# separate so the two gates can evolve independently if needed.
FILAMENTOUS_ASCO_CLASSES = {
    'sordariomycetes', 'eurotiomycetes', 'dothideomycetes', 'leotiomycetes',
    'pezizomycetes', 'orbiliomycetes', 'lecanoromycetes', 'xylonomycetes',
    'lichinomycetes', 'coniocybomycetes', 'geoglossomycetes', 'arthoniomycetes',
}


def classify_guild_taxonomic(
    phylum: str = '',
    clazz: str = '',
    genus: str = '',
    species: str = '',
    lifestyle: str = '',
) -> str:
    """4-way phylogeny-aware guild label: AMF / Filamentous / Basidio / Yeast.

    Resolution priority:
      1. Glomeromycota phylum (or AMF genus list) -> AMF
      2. Saccharomycotina/Schizosaccharomycotina classes -> Yeast
      3. Pezizomycotina classes -> Filamentous
      4. Basidiomycota phylum -> Basidio
      5. Fall through to genus/lifestyle fallback for 'Other'/edge cases.

    The 13-sample panel resolves cleanly into AMF + Filamentous + Yeast; the
    165-sample panel adds the dominant Basidio (74/165) bucket.
    """
    p = (phylum or '').strip().lower()
    c = (clazz or '').strip().lower()
    if p == 'glomeromycota':
        return 'AMF'
    if c in YEAST_CLASSES:
        return 'Yeast'
    if c in FILAMENTOUS_ASCO_CLASSES:
        return 'Filamentous'
    if p == 'basidiomycota':
        return 'Basidio'
    # Genus list still catches AMF genera even when phylum/class were not
    # populated upstream (some hierarchy_manifest rows have '.' for class).
    if genus:
        g = genus.split()[0].lstrip('[').rstrip(']')
        if GUILD_BY_GENUS.get(g) == 'AMF':
            return 'AMF'
    return classify_guild(genus, species, lifestyle)


def _fisher_one_tailed_right(a: int, b: int, c: int, d: int) -> float:
    """Right-tailed Fisher's exact P(X>=a) for a 2x2 contingency table:
        [[a, b], [c, d]]. Uses log-space hypergeometric tail sum to stay
    well-behaved at large N."""
    n = a + b + c + d
    if n == 0:
        return 1.0
    row1 = a + b
    col1 = a + c
    a_min = max(0, col1 - (n - row1))
    a_max = min(row1, col1)

    def logC(n_, k_):
        return math.lgamma(n_ + 1) - math.lgamma(k_ + 1) - math.lgamma(n_ - k_ + 1)

    log_denom = logC(n, col1)
    log_tail = []
    for k in range(a, a_max + 1):
        logp = logC(row1, k) + logC(n - row1, col1 - k) - log_denom
        log_tail.append(logp)
    if not log_tail:
        return 1.0
    m = max(log_tail)
    p = math.exp(m) * sum(math.exp(x - m) for x in log_tail)
    return min(1.0, max(0.0, p))


def _truthy(value: str) -> bool:
    return value.strip().lower() in {'yes', 'true', '1', 'y'}


def _is_te_class(ec: str) -> bool:
    return ec.upper() in TE_CLASSES


def _shard_iter(shard_dir: Path):
    """Yield (query_asm, biology_findings_path) for shards that have one.

    Skips shards where the mycosv binary failed mid-run (MYCOSV_FAILED.txt
    marker) - those produce a near-empty checkpoint that pollutes panel
    aggregates with no_truth/zero-precision rows.
    """
    for q in sorted(p for p in shard_dir.iterdir() if p.is_dir()):
        if (q / 'MYCOSV_FAILED.txt').exists():
            print(f'[cross-guild] skipping failed shard {q.name} (MYCOSV_FAILED.txt present)',
                  flush=True)
            continue
        bf = q / 'biology_findings.tsv'
        if bf.exists() and bf.stat().st_size > 0:
            yield q.name, bf


def _locate_query_manifest(shard_dir: Path) -> Path | None:
    """Find the query_manifest.tsv that produced this by_query/ directory.

    Conventional layout:
        <run>/prepared/query_manifest.tsv
        <run>/<mode>/by_query/      <-- shard_dir is here
    so the prepared/ sibling lives at shard_dir.parent.parent / 'prepared'.
    Falls back to scanning ancestors for prepared/query_manifest.tsv.
    """
    for ancestor in [shard_dir, shard_dir.parent, shard_dir.parent.parent]:
        candidate = ancestor / 'prepared' / 'query_manifest.tsv'
        if candidate.exists():
            return candidate
    for ancestor in shard_dir.parents:
        candidate = ancestor / 'prepared' / 'query_manifest.tsv'
        if candidate.exists():
            return candidate
    return None


def _load_query_class_map(manifest_path: Path | None) -> dict[str, dict[str, str]]:
    """query_asm -> {phylum, class, ...} from prepared/query_manifest.tsv.

    Falls back to {} when the manifest is missing; classify_guild_taxonomic
    then degrades to its genus-list / lifestyle fallback. With the manifest in
    place the 165-sample run resolves into
    4 guilds (AMF / Filamentous / Basidio / Yeast) instead of 2.
    """
    out: dict[str, dict[str, str]] = {}
    if manifest_path is None or not manifest_path.exists():
        return out
    with manifest_path.open() as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            asm = (row.get('query_asm') or '').strip()
            if not asm:
                continue
            for alias in query_asm_aliases(asm):
                out.setdefault(alias, row)
    return out


def run_cross_guild(shard_dir: Path, out_dir: Path) -> int:
    """Aggregate biology_findings.tsv across shards, stratify by guild AND by
    pangenome scope, and emit cross-guild summary + enrichment + pangenome-lift
    tables that connect SV biology to ecology and gene-expression context."""
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = _locate_query_manifest(shard_dir)
    class_map = _load_query_class_map(manifest_path)
    if class_map:
        print(f'[cross-guild] taxonomy lookup loaded from {manifest_path} '
              f'({len(class_map)} query rows)', flush=True)
    else:
        print('[cross-guild] no query_manifest.tsv found; guild classifier '
              'will degrade to genus-list fallback (2-way bucketing)', flush=True)

    rows: list[dict[str, str]] = []
    n_shards = 0
    for qname, bf in _shard_iter(shard_dir):
        n_shards += 1
        with bf.open(newline='') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for r in reader:
                r['_query_asm'] = qname
                rows.append(r)
    if not rows:
        print(f'[cross-guild] no biology_findings.tsv shards found under {shard_dir}')
        return 1
    print(f'[cross-guild] read {len(rows)} biology rows across {n_shards} shards', flush=True)

    # Classify each row by guild + pangenome scope.
    for r in rows:
        qasm = r.get('_query_asm') or r.get('query_asm') or ''
        manifest_row = {}
        for alias in query_asm_aliases(qasm):
            if alias in class_map:
                manifest_row = class_map[alias]
                break
        r['_guild'] = classify_guild_taxonomic(
            phylum=manifest_row.get('phylum') or r.get('phylum') or '',
            clazz=manifest_row.get('class') or '',
            genus=r.get('genus', '') or manifest_row.get('genus') or '',
            species=r.get('species', '') or manifest_row.get('species') or '',
            lifestyle=r.get('lifestyle', '') or manifest_row.get('lifestyle') or '',
        )
        # pangenome scope: mycosv_unique=yes means the call was NOT recovered by
        # any single-reference comparator; single_reference_equivalent=yes means
        # it would be (so single-ref methods would also find it).
        r['_pangenome_only'] = _truthy(r.get('mycosv_unique', '') or '')
        r['_single_ref_eq'] = _truthy(r.get('single_reference_equivalent', '') or '')
        r['_has_eco'] = any(
            (r.get(k, '') or '').strip() not in ('.', '')
            for k in ('lifestyle','ecological_trait','substrate_or_host','trophic_mode')
        )
        r['_has_expr'] = _truthy(r.get('expression_supported', '') or '')
        r['_has_gene'] = (r.get('affected_gene', '.') or '.') not in ('.', '')
        r['_is_hgt'] = (r.get('candidate_type', '') or '').strip() == 'hgt_candidate' \
            or (r.get('element_class', '').upper() in ('HGT','STARSHIP'))
        r['_is_rip'] = r.get('element_class', '').upper() == 'RIP'
        r['_is_te'] = _is_te_class(r.get('element_class', '') or '')

    # ---- Per-guild summary -------------------------------------------------
    per_guild = defaultdict(lambda: {
        'n_shards': set(), 'n_candidates': 0,
        'n_hgt': 0, 'n_rip': 0, 'n_te': 0,
        'pangenome_only': 0, 'single_ref_eq': 0,
        'hgt_pangenome_only': 0, 'rip_pangenome_only': 0, 'te_pangenome_only': 0,
        'with_ecology': 0, 'with_expression': 0, 'with_gene': 0,
        'genera': Counter(),
        'svtypes': Counter(),
        'priorities': [],
    })
    for r in rows:
        g = r['_guild']
        b = per_guild[g]
        b['n_shards'].add(r['_query_asm'])
        b['n_candidates'] += 1
        if r['_is_hgt']: b['n_hgt'] += 1
        if r['_is_rip']: b['n_rip'] += 1
        if r['_is_te']: b['n_te'] += 1
        if r['_pangenome_only']: b['pangenome_only'] += 1
        if r['_single_ref_eq']: b['single_ref_eq'] += 1
        if r['_is_hgt'] and r['_pangenome_only']: b['hgt_pangenome_only'] += 1
        if r['_is_rip'] and r['_pangenome_only']: b['rip_pangenome_only'] += 1
        if r['_is_te'] and r['_pangenome_only']: b['te_pangenome_only'] += 1
        if r['_has_eco']: b['with_ecology'] += 1
        if r['_has_expr']: b['with_expression'] += 1
        if r['_has_gene']: b['with_gene'] += 1
        g_first = (r.get('genus','') or r.get('species','')).split()[0:1]
        if g_first:
            b['genera'][g_first[0]] += 1
        b['svtypes'][r.get('svtype','.')] += 1
        try:
            b['priorities'].append(int(r.get('priority','0')))
        except ValueError:
            pass

    summary_path = out_dir / 'cross_guild_summary.tsv'
    with summary_path.open('w', newline='') as fh:
        fields = ['guild','n_shards','n_candidates','n_hgt','n_rip','n_te',
                  'pangenome_only','single_ref_eq','pangenome_only_fraction',
                  'hgt_pangenome_only','rip_pangenome_only','te_pangenome_only',
                  'hgt_pangenome_only_fraction','rip_pangenome_only_fraction',
                  'with_ecology','with_expression','with_gene',
                  'mean_priority','top_svtypes','top_genera']
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader()
        for g, b in sorted(per_guild.items()):
            n = max(1, b['n_candidates'])
            w.writerow({
                'guild': g,
                'n_shards': len(b['n_shards']),
                'n_candidates': b['n_candidates'],
                'n_hgt': b['n_hgt'],
                'n_rip': b['n_rip'],
                'n_te': b['n_te'],
                'pangenome_only': b['pangenome_only'],
                'single_ref_eq': b['single_ref_eq'],
                'pangenome_only_fraction': f'{b["pangenome_only"]/n:.3f}',
                'hgt_pangenome_only': b['hgt_pangenome_only'],
                'rip_pangenome_only': b['rip_pangenome_only'],
                'te_pangenome_only': b['te_pangenome_only'],
                'hgt_pangenome_only_fraction': (f'{b["hgt_pangenome_only"]/b["n_hgt"]:.3f}' if b['n_hgt'] else '.'),
                'rip_pangenome_only_fraction': (f'{b["rip_pangenome_only"]/b["n_rip"]:.3f}' if b['n_rip'] else '.'),
                'with_ecology': b['with_ecology'],
                'with_expression': b['with_expression'],
                'with_gene': b['with_gene'],
                'mean_priority': (f'{sum(b["priorities"])/len(b["priorities"]):.2f}' if b['priorities'] else '.'),
                'top_svtypes': ';'.join(f'{k}:{v}' for k,v in b['svtypes'].most_common(4)),
                'top_genera': ';'.join(f'{k}:{v}' for k,v in b['genera'].most_common(4)),
            })
    print(f'[cross-guild] wrote {summary_path}', flush=True)

    # ---- Per-guild x feature Fisher enrichment -----------------------------
    # For each (guild, feature), test whether the feature is enriched in this
    # guild vs the rest of the panel (2x2 right-tailed Fisher).
    features = [
        ('hgt_candidate',           lambda r: r['_is_hgt']),
        ('rip_candidate',           lambda r: r['_is_rip']),
        ('te_or_repeat',            lambda r: r['_is_te']),
        ('pangenome_only',          lambda r: r['_pangenome_only']),
        ('single_ref_equiv',        lambda r: r['_single_ref_eq']),
        ('ecology_linked',          lambda r: r['_has_eco']),
        ('expression_supported',    lambda r: r['_has_expr']),
        ('gene_proximal',           lambda r: r['_has_gene']),
        ('hgt_AND_pangenome_only',  lambda r: r['_is_hgt'] and r['_pangenome_only']),
        ('hgt_AND_ecology_linked',  lambda r: r['_is_hgt'] and r['_has_eco']),
        ('hgt_AND_expression',      lambda r: r['_is_hgt'] and r['_has_expr']),
        ('pangenome_AND_expression',lambda r: r['_pangenome_only'] and r['_has_expr']),
        ('pangenome_AND_ecology',   lambda r: r['_pangenome_only'] and r['_has_eco']),
    ]
    enrichment_path = out_dir / 'cross_guild_enrichment.tsv'
    n_tests = len(features) * len(per_guild)
    with enrichment_path.open('w', newline='') as fh:
        fields = ['guild','feature','in_guild_yes','in_guild_no','out_guild_yes','out_guild_no','odds_ratio','fisher_p','bonferroni_p','interpretation']
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader()
        for guild_name in sorted(per_guild.keys()):
            for feat_name, pred in features:
                in_yes = sum(1 for r in rows if r['_guild'] == guild_name and pred(r))
                in_no  = sum(1 for r in rows if r['_guild'] == guild_name and not pred(r))
                out_yes= sum(1 for r in rows if r['_guild'] != guild_name and pred(r))
                out_no = sum(1 for r in rows if r['_guild'] != guild_name and not pred(r))
                # odds ratio (Haldane-Anscombe corrected)
                a, b, c, d = in_yes + 0.5, in_no + 0.5, out_yes + 0.5, out_no + 0.5
                or_val = (a * d) / (b * c)
                p = _fisher_one_tailed_right(in_yes, in_no, out_yes, out_no)
                bonf = min(1.0, p * n_tests)
                interp = '.'
                if in_yes > 0 and bonf < 0.05 and or_val > 1.0:
                    interp = f'enriched_in_{guild_name}'
                elif in_yes > 0 and bonf < 0.05 and or_val < 1.0:
                    interp = f'depleted_in_{guild_name}'
                w.writerow({
                    'guild': guild_name, 'feature': feat_name,
                    'in_guild_yes': in_yes, 'in_guild_no': in_no,
                    'out_guild_yes': out_yes, 'out_guild_no': out_no,
                    'odds_ratio': f'{or_val:.4g}',
                    'fisher_p': f'{p:.4g}',
                    'bonferroni_p': f'{bonf:.4g}',
                    'interpretation': interp,
                })
    print(f'[cross-guild] wrote {enrichment_path}', flush=True)

    # ---- Pangenome lift per guild x candidate axis -------------------------
    # The headline number for the writeup: "X% of HGT (RIP, TE, off-ref) in
    # guild Y would have been MISSED by single-reference methods."
    lift_path = out_dir / 'cross_guild_pangenome_lift.tsv'
    with lift_path.open('w', newline='') as fh:
        fields = ['guild','candidate_axis','pangenome_only','single_ref_eq','total','pangenome_lift_fraction','interpretation']
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader()
        axes = [
            ('all_candidates', lambda r: True),
            ('hgt_candidate',  lambda r: r['_is_hgt']),
            ('rip_candidate',  lambda r: r['_is_rip']),
            ('te_or_repeat',   lambda r: r['_is_te']),
            ('ecology_linked', lambda r: r['_has_eco']),
            ('expression_supported', lambda r: r['_has_expr']),
            ('gene_proximal',  lambda r: r['_has_gene']),
        ]
        for guild_name in sorted(per_guild.keys()):
            for ax_name, pred in axes:
                pan = sum(1 for r in rows if r['_guild']==guild_name and pred(r) and r['_pangenome_only'])
                sre = sum(1 for r in rows if r['_guild']==guild_name and pred(r) and r['_single_ref_eq'])
                tot = sum(1 for r in rows if r['_guild']==guild_name and pred(r))
                lift = (pan / tot) if tot else 0.0
                interp = '.'
                if tot >= 20 and lift >= 0.5:
                    interp = f'{guild_name}/{ax_name}_{lift*100:.0f}pct_pangenome_only'
                w.writerow({
                    'guild': guild_name, 'candidate_axis': ax_name,
                    'pangenome_only': pan, 'single_ref_eq': sre, 'total': tot,
                    'pangenome_lift_fraction': f'{lift:.3f}',
                    'interpretation': interp,
                })
    print(f'[cross-guild] wrote {lift_path}', flush=True)

    # ---- Brief textual summary JSON ----------------------------------------
    summary_json_path = out_dir / 'cross_guild_summary.json'
    with summary_json_path.open('w') as fh:
        json.dump({
            'n_shards': n_shards,
            'n_rows': len(rows),
            'guilds': {
                g: {
                    'n_candidates': b['n_candidates'],
                    'pangenome_only_fraction': (b['pangenome_only']/max(1,b['n_candidates'])),
                    'hgt_n': b['n_hgt'], 'rip_n': b['n_rip'], 'te_n': b['n_te'],
                    'with_ecology_fraction': (b['with_ecology']/max(1,b['n_candidates'])),
                    'with_expression_fraction': (b['with_expression']/max(1,b['n_candidates'])),
                    'top_svtypes': dict(b['svtypes'].most_common(5)),
                    'top_genera': dict(b['genera'].most_common(5)),
                }
                for g, b in per_guild.items()
            },
        }, fh, indent=2, sort_keys=True, default=str)
    print(f'[cross-guild] wrote {summary_json_path}', flush=True)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description='Summarize new-biology candidate signals from fungi_graphsv_tol outputs.')
    # Cross-guild aggregation mode (replaces --vcf / --hits / --out-tsv when set).
    ap.add_argument('--cross-guild-shard-dir', type=Path,
                    help='Aggregate biology_findings.tsv across all per-query shards under this directory '
                         '(e.g. .../assembly/by_query/). Stratifies by guild AND by pangenome scope, '
                         'and ties results back to ecology + gene-expression context. Writes '
                         'cross_guild_{summary,enrichment,pangenome_lift}.tsv into --cross-guild-out-dir.')
    ap.add_argument('--cross-guild-out-dir', type=Path,
                    help='Output directory for the cross-guild tables. Defaults to '
                         '<cross-guild-shard-dir>/../combined/cross_guild/.')
    ap.add_argument('--vcf', type=Path)
    ap.add_argument('--hits', type=Path)
    ap.add_argument('--ancestral', type=Path)
    ap.add_argument('--expression-tsv', type=Path)
    ap.add_argument('--expression-long-tsv', type=Path,
                    help='Sample-level long-form expression table with query_asm, gene_id, expression, and condition/group columns.')
    ap.add_argument('--gene-annotations', type=Path,
                    help='Gene coordinate TSV with query_asm, query_contig, gene_id, gene_name, start, end.')
    ap.add_argument('--ecological-traits', type=Path,
                    help='Ecological trait TSV keyed by query_asm, typically prepared/ecological_traits.tsv.')
    ap.add_argument('--fungaltraits-csv', type=Path,
                    help='Cached FungalTraits CSV used as a species/genus fallback for ecological_trait columns.')
    ap.add_argument('--derived-expression-out', type=Path,
                    help='Optional TSV path to write candidate-level expression quantification derived from --expression-long-tsv.')
    ap.add_argument('--expression-window-bp', type=int, default=5000)
    ap.add_argument('--expression-group-a')
    ap.add_argument('--expression-group-b')
    ap.add_argument('--expression-min-reps', type=int, default=2)
    ap.add_argument('--expression-pseudocount', type=float, default=1.0)
    ap.add_argument('--query-metadata', type=Path)
    ap.add_argument('--out-tsv', type=Path)
    ap.add_argument('--summary-json', type=Path)
    ap.add_argument('--phylum', default='unknown')
    ap.add_argument('--top-n', type=int, default=50)
    args = ap.parse_args()

    # Cross-guild mode short-circuits the per-shard pipeline. Used by the
    # combine phase to synthesize the panel-wide story; works for any number
    # of shards (15-query and 165-query runs alike - just point it at the
    # by_query/ directory).
    if args.cross_guild_shard_dir is not None:
        shard_dir = args.cross_guild_shard_dir.resolve()
        if not shard_dir.is_dir():
            ap.error(f'--cross-guild-shard-dir not a directory: {shard_dir}')
        out_dir = (args.cross_guild_out_dir
                   or (shard_dir.parent / 'combined' / 'cross_guild')).resolve()
        return run_cross_guild(shard_dir, out_dir)

    # Per-shard mode: --vcf / --out-tsv / --summary-json are all required.
    missing = [name for name, val in (
        ('--vcf', args.vcf), ('--out-tsv', args.out_tsv), ('--summary-json', args.summary_json),
    ) if val is None]
    if missing:
        ap.error(f'per-shard mode requires: {", ".join(missing)}')

    # gene_annotations alone is now valid: it powers the nearest-gene fallback
    # so expression_gene / expression_distance_bp populate without an RNA-seq
    # matrix. The reverse - expression_long without gene_annotations - still
    # cannot resolve gene-coord lookups and is an error.
    if args.expression_long_tsv and not args.gene_annotations:
        ap.error('--expression-long-tsv requires --gene-annotations to map gene_id -> contig coordinates')

    hits = load_hits(args.hits)
    meta = load_query_meta(args.query_metadata)
    ancestral = load_ancestral(args.ancestral)
    records = load_vcf_records(args.vcf)
    ecological_traits = load_ecological_traits(args.ecological_traits)
    fungal_by_species, fungal_by_genus = load_fungaltraits_csv(args.fungaltraits_csv)

    expression = load_expression(args.expression_tsv)
    # Always load gene_annotations even when no expression_long is supplied -
    # it lets the per-candidate fallback below populate expression_gene /
    # expression_distance_bp from the nearest annotated gene, so the operator
    # at least sees which gene is closest to each breakpoint without needing
    # an RNA-seq matrix (which most fungal panels don't have publicly).
    gene_annotations_lookup = load_gene_annotations(args.gene_annotations)
    if args.expression_long_tsv and args.gene_annotations:
        derived_expression = derive_expression_support_from_quant(
            records,
            hits,
            gene_annotations_lookup,
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
        qasm = hit.get('query_asm') or info.get('QASM') or info.get('QUERY_ASM', '.')
        meta_row = {}
        for alias in query_asm_aliases(qasm):
            if alias in meta:
                meta_row = meta[alias]
                break
        # RIP is a Pezizomycotina-restricted defense mechanism. The cpp
        # detector fires on TpA-enriched / CpA-depleted sequence which can
        # also arise in yeast / AMF / Basidiomycota for non-RIP reasons, so
        # relabel to REPEAT outside RIP-supporting classes.
        if ec == 'RIP' and not rip_supported_by_clade(meta_row):
            raw_ec = 'REPEAT'
            ec = 'REPEAT'
        scenario = meta_row.get('scenario', '.') if meta_row else '.'
        architecture = meta_row.get('architecture', '.') if meta_row else '.'
        lifestyle = meta_row.get('lifestyle', '.') if meta_row else '.'
        eco = ecological_context(qasm, meta_row, ecological_traits, fungal_by_species, fungal_by_genus)
        anc = ancestral.get((qasm, chrom)) or ancestral.get(('.', chrom))
        expr = expression_for_candidate(expression, qasm, chrom, pos, end, svtype)
        # Even when no expression matrix was supplied, the gene_annotations.tsv
        # alone tells us which annotated gene the breakpoint is closest to.
        # Surface that as expression_gene + expression_distance_bp; leave
        # expression_log2_fc / expression_padj as '.' since we have no measurement.
        nearest = None
        affected_locus = chrom
        ref_contig = hit.get('ref_contig') or info.get('REFCONTIG') or '.'
        ref_pos = hit.get('ref_pos') or info.get('REFPOS') or ''
        ref_end = hit.get('ref_end') or info.get('REFEND') or ref_pos
        ref_asm = hit.get('ref_asm') or info.get('CLADE') or info.get('CL') or '.'
        if expr is None or not expr.get('best_gene'):
            nearest = nearest_gene_for_locus(
                gene_annotations_lookup,
                [qasm],
                chrom,
                pos,
                end,
            )
        if nearest is None and ref_contig not in {'', '.'} and ref_pos not in {'', '.', '0'}:
            affected_locus = f"{ref_contig}:{ref_pos}-{ref_end or ref_pos}"
            nearest = nearest_gene_for_locus(
                gene_annotations_lookup,
                [qasm, meta_row.get('benchmark_ref_asm', ''), ref_asm, meta_row.get('benchmark_ref_fasta', '')],
                ref_contig,
                ref_pos,
                ref_end or ref_pos,
            )
        candidate_type, priority, rationale = classify_candidate(svtype, ec, annot, anc, expr)
        if candidate_type == 'other':
            continue
        clades = ','.join(sorted(anc['clades'])) if anc else '.'
        clade_ranks = ','.join(sorted(anc['ranks'])) if anc else '.'
        example = choose_functional_example(candidate_type, svtype, ec, scenario, anc, expr)
        expr_supported = 'yes' if expr and expr.get('supported') else 'no'
        expr_gene = (expr.get('best_gene') if expr and expr.get('best_gene') else
                     (nearest['gene_name'] if nearest else '.'))
        expr_distance = (expr.get('distance_bp') if expr and expr.get('best_gene') else
                         (nearest['distance_bp'] if nearest else '.'))
        expr_log2_fc = expr.get('log2_fc', '.') if expr else '.'
        expr_padj = expr.get('padj', '.') if expr else '.'
        expr_condition = expr.get('condition', '.') if expr else '.'
        affected_gene_id = nearest.get('gene_id') if nearest else '.'
        affected_gene = nearest.get('gene_name') if nearest else expr_gene
        affected_distance = nearest.get('distance_bp') if nearest else expr_distance
        affected_product = nearest.get('product') if nearest else '.'
        affected_biotype = nearest.get('biotype') if nearest else '.'
        rows.append({
            'priority': priority,
            'candidate_type': candidate_type,
            'phylum': args.phylum,
            'query_asm': qasm,
            'query_contig': chrom,
            'scenario': scenario,
            'species': eco['species'],
            'lifestyle': lifestyle,
            'architecture': architecture,
            'ecological_trait': eco['ecological_trait'],
            'secondary_lifestyle': eco['secondary_lifestyle'],
            'trophic_mode': eco['trophic_mode'],
            'substrate_or_host': eco['substrate_or_host'],
            'svtype': svtype,
            'element_class': raw_ec,
            'mge_subtype': mge_subtype(ec),
            'hgt_flag': 'yes' if is_hgt_candidate(ec, svtype, annot) else 'no',
            'novelty': annot,
            'pos': int(pos),
            'end': int(end),
            'ref_target': hit.get('ref_asm', info.get('CL', '.')),
            'affected_locus': affected_locus,
            'affected_gene_id': affected_gene_id,
            'affected_gene': affected_gene,
            'affected_gene_distance_bp': affected_distance,
            'affected_gene_biotype': affected_biotype,
            'affected_gene_product': affected_product,
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
    rows = select_diverse_rows(rows, max(1, args.top_n))

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_tsv.open('w', newline='') as fh:
        fieldnames = [
            'priority', 'candidate_type', 'phylum', 'query_asm', 'query_contig', 'scenario', 'species',
            'lifestyle', 'architecture', 'ecological_trait', 'secondary_lifestyle', 'trophic_mode',
            'substrate_or_host', 'svtype', 'element_class', 'mge_subtype', 'hgt_flag', 'novelty', 'pos', 'end',
            'ref_target', 'affected_locus', 'affected_gene_id', 'affected_gene', 'affected_gene_distance_bp',
            'affected_gene_biotype', 'affected_gene_product', 'alignment_mode', 'ancestral_clades',
            'ancestral_ranks', 'ancestral_breakpoints',
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
