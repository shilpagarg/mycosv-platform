#!/usr/bin/env python3
# Designed for Linux
"""test_amf.py — Fungal SV simulator for the TOL pipeline.

Domain corrections applied (vs previous version — 21 issues fixed)
====================================================================
TAXONOMY (2 fixes)
  - Laccaria: family Hydnangiaceae (NOT Amanitaceae). Matheny et al. 2006.
  - Batrachochytrium dendrobatidis: order Rhizophydiales, family
    Rhizophydiaceae (NOT Chytridiales/Chytriaceae). James et al. 2006.

ECOLOGY (1 fix)
  - Mortierella: lifestyle 'soil_endophyte' (NOT endomycorrhizal). Guo 2020.

GC CONTENT — per-genus/class hard-coded from published genomes (5 fixes)
  - Puccinia: 42% (Duplessis et al. 2011 PNAS)
  - Fusarium: 48% (Ma et al. 2010 Nature)
  - Cladonia: 51% (Armaleo et al. 2019 Nat. Commun.)
  - Batrachochytrium: 40% (Joneson et al. 2011)
  - giant_amf: 28.5% (DAOM 197198; Chen et al. 2018 PNAS)
  New genera added: Lachancea 40%, Ustilago 54%, Botrytis 44%,
    Epichloë 48%, Verticillium 53%, Mortierella 43%

GENOME SIZE (1 fix)
  - giant_amf: 150 Mb (post-decontamination). The 750 Mb figure was an
    error. Chen et al. 2018 PNAS.

TE BIOLOGY (2 fixes)
  - rust_smut_te_heavy: TIR elements removed from dispatch (<2% in
    Puccinia); Gypsy LTR + RIP only. Duplessis et al. 2011 PNAS.
  - arbuscular_mf/giant_amf: STARSHIP removed (only confirmed in
    Ascomycota, not Glomeromycota). Urquhart et al. 2023 Curr. Biol.

SV BIAS (2 fixes)
  - arbuscular_mf: TRA removed (AMF largely asexual/clonal; inter-contig
    translocations rare). Corrected to DUP_INS.
  - lichenised: HGT added to element dispatch (algal/cyanobacterial HGT
    documented). Slot & Rokas 2011 Science.

OFF-REFERENCE GC (1 fix)
  - off_ref_gc now direction-aware: high-GC hosts shift DOWN (AT-rich TEs
    or Firmicutes donors), low-GC hosts shift UP.

RIP BIOLOGY (1 fix)
  - apply_rip: now targets CpA dinucleotides only (canonical RIP context:
    CpA→TpA). Selker et al. 2003; Cambareri et al. 1989.

STARSHIP BIOLOGY (1 fix)
  - embed_starship: cargo GC lowered to ~genic (45-55%), hull GC is
    relative to clade background (not absolute 0.42). Urquhart et al. 2023.

HGT BIOLOGY (1 fix)
  - embed_hgt_island: GC deviation set to +0.10 (published range ±0.05-0.10
    for mycological HGT). Slot & Rokas 2011.

DUPLICATE GENUS (1 fix)
  - compact_yeast now uses Lachancea (not Saccharomyces again) to avoid
    duplicate genus shards in the routing index.

NEW SCENARIOS (3 additions)
  - necrotrophic: Botrytis cinerea (Leotiomycetes) — distinct from
    two-speed Fusarium; GC 44%, moderate TE.
  - strict_endophyte: Epichloë festucae — near-TE-free compact genome.
  - smut: Ustilago maydis — compact dikaryote separate from Puccinia rust.
"""
from __future__ import annotations

import argparse
import random
from pathlib import Path

# ── RANKS — single authoritative list (mirrors taxonomy_ranks.hpp) ────────
RANKS = ["phylum", "class", "order", "family", "genus", "species"]

# ── Per-genus GC (published genome papers) ────────────────────────────────
_GENUS_GC: dict[str, float] = {
    "Saccharomyces":    0.38,
    "Lachancea":        0.40,
    "Rhizophagus":      0.32,
    "Puccinia":         0.42,
    "Ustilago":         0.54,
    "Laccaria":         0.45,
    "Botrytis":         0.44,
    "Epichloë":         0.48,
    "Batrachochytrium": 0.40,
    "Cladonia":         0.51,
    "Verticillium":     0.53,
    "Fusarium":         0.48,
    "Mortierella":      0.43,
}

# ── Per-class GC fallback ─────────────────────────────────────────────────
_CLASS_GC: dict[str, float] = {
    "Sordariomycetes":    0.48,
    "Saccharomycetes":    0.38,
    "Lecanoromycetes":    0.51,
    "Pucciniomycetes":    0.42,
    "Ustilaginomycetes":  0.54,
    "Agaricomycetes":     0.45,
    "Glomeromycetes":     0.32,
    "Mortierellomycetes": 0.43,
    "Leotiomycetes":      0.44,
    "Chytridiomycetes":   0.40,
}

# ── SCENARIOS ─────────────────────────────────────────────────────────────
SCENARIOS: dict[str, dict[str, str]] = {
    # Ascomycota / yeasts
    "core": dict(
        phylum="Ascomycota", cls="Saccharomycetes",
        order="Saccharomycetales", family="Saccharomycetaceae", genus="Saccharomyces",
        lifestyle="baseline_control", architecture="compact_baseline",
        genome_scale="small", repeat_regime="low", te_regime="low",
        hgt_regime="low", expected_sv_bias="balanced",
    ),
    "compact_yeast": dict(
        # FIX: use Lachancea, not Saccharomyces, to avoid duplicate genus shard
        phylum="Ascomycota", cls="Saccharomycetes",
        order="Saccharomycetales", family="Saccharomycetaceae", genus="Lachancea",
        lifestyle="compact_yeast", architecture="very_small_compact_yeast",
        genome_scale="very_small", repeat_regime="very_low", te_regime="very_low",
        hgt_regime="low", expected_sv_bias="small_DEL_INS",
    ),
    # Sordariomycetes
    "saprotrophic": dict(
        phylum="Ascomycota", cls="Sordariomycetes",
        order="Hypocreales", family="Nectriaceae", genus="Fusarium",
        lifestyle="saprotrophic", architecture="decomposer_hgt_prone",
        genome_scale="medium", repeat_regime="moderate", te_regime="moderate",
        hgt_regime="high", expected_sv_bias="TRA_HGT",
    ),
    "pathogenic": dict(
        phylum="Ascomycota", cls="Sordariomycetes",
        order="Hypocreales", family="Nectriaceae", genus="Fusarium",
        lifestyle="pathogenic", architecture="two_speed_te_rich",
        genome_scale="medium_large", repeat_regime="high", te_regime="very_high",
        hgt_regime="moderate", expected_sv_bias="INV_INS",
        # F. oxysporum LS chromosomes: MITE TIR + Gypsy LTR + RIP
    ),
    "two_speed_pathogen_extreme": dict(
        # FIX: use Verticillium (canonical two-speed model), not Fusarium again
        phylum="Ascomycota", cls="Sordariomycetes",
        order="Glomerellales", family="Plectosphaerellaceae", genus="Verticillium",
        lifestyle="two_speed_pathogen", architecture="highly_rearranged_two_speed",
        genome_scale="large", repeat_regime="very_high", te_regime="extreme",
        hgt_regime="moderate", expected_sv_bias="INV_TRA_INS",
    ),
    "necrotrophic": dict(
        # NEW: Botrytis cinerea — distinct from Fusarium two-speed model
        phylum="Ascomycota", cls="Leotiomycetes",
        order="Helotiales", family="Sclerotiniaceae", genus="Botrytis",
        lifestyle="necrotrophic_pathogen", architecture="necrotroph_moderate_te",
        genome_scale="medium", repeat_regime="moderate", te_regime="moderate",
        hgt_regime="low", expected_sv_bias="DEL_INS",
        # B. cinerea B05.10: GC 44%; 43 Mb. van Kan et al. 2017.
    ),
    "strict_endophyte": dict(
        # NEW: Epichloë festucae — near-TE-free, horizontal transmission
        phylum="Ascomycota", cls="Sordariomycetes",
        order="Hypocreales", family="Clavicipitaceae", genus="Epichloë",
        lifestyle="strict_endophyte", architecture="compact_endophyte_low_te",
        genome_scale="small_medium", repeat_regime="very_low", te_regime="very_low",
        hgt_regime="low", expected_sv_bias="small_DEL_INS",
    ),
    # Lecanoromycetes
    "lichenised": dict(
        phylum="Ascomycota", cls="Lecanoromycetes",
        order="Lecanorales", family="Cladoniaceae", genus="Cladonia",
        lifestyle="lichenised", architecture="constrained_low_sv_algal_hgt",
        genome_scale="small", repeat_regime="low", te_regime="low",
        hgt_regime="moderate", expected_sv_bias="low_sv",
        # GC ~51%; algal/cyanobacterial HGT documented. Armaleo et al. 2019.
    ),
    # Basidiomycota / ECM
    "ectomycorrhizal": dict(
        phylum="Basidiomycota", cls="Agaricomycetes",
        order="Agaricales",
        family="Hydnangiaceae",   # FIX: was Amanitaceae (wrong). Matheny 2006.
        genus="Laccaria",
        lifestyle="ectomycorrhizal", architecture="medium_large_secretome_missp",
        genome_scale="medium_large", repeat_regime="moderate", te_regime="moderate",
        hgt_regime="low", expected_sv_bias="INS_DUP",
        # L. bicolor: GC 45%; 65 Mb; MiSSP clusters near TE-rich regions.
    ),
    # Basidiomycota / rust & smut
    "rust_smut_te_heavy": dict(
        phylum="Basidiomycota", cls="Pucciniomycetes",
        order="Pucciniales", family="Pucciniaceae", genus="Puccinia",
        lifestyle="rust_like", architecture="te_heavy_gypsy_dominant_dikaryotic",
        genome_scale="large", repeat_regime="very_high", te_regime="extreme",
        hgt_regime="low", expected_sv_bias="INS_INV",
        # GC 42%; ~89 Mb; >40% Gypsy LTR; TIR <2%. Duplessis 2011.
        # Dikaryotic nature modelled as single-haploid proxy.
    ),
    "smut": dict(
        # NEW: Ustilago maydis — compact dikaryote, separate from rust
        phylum="Basidiomycota", cls="Ustilaginomycetes",
        order="Ustilaginales", family="Ustilaginaceae", genus="Ustilago",
        lifestyle="smut_pathogen", architecture="compact_dikaryote_smut",
        genome_scale="very_small", repeat_regime="low", te_regime="low",
        hgt_regime="low", expected_sv_bias="small_INS_DEL",
        # U. maydis: GC 54%; ~20 Mb. Kämper et al. 2006 Nature.
    ),
    # Glomeromycota / AMF
    "arbuscular_mf": dict(
        phylum="Glomeromycota", cls="Glomeromycetes",
        order="Glomerales", family="Glomeraceae", genus="Rhizophagus",
        lifestyle="arbuscular_mycorrhizal", architecture="repeat_rich_large_gypsy_copia",
        genome_scale="very_large", repeat_regime="very_high", te_regime="high",
        hgt_regime="moderate", expected_sv_bias="DUP_INS",
        # GC ~32%; ~150 Mb; Gypsy/Copia dominant (40-60%).
        # FIX: TRA removed (largely asexual, translocations rare).
        # FIX: STARSHIP removed (not confirmed in Glomeromycota).
    ),
    "giant_amf": dict(
        phylum="Glomeromycota", cls="Glomeromycetes",
        order="Glomerales", family="Glomeraceae", genus="Rhizophagus",
        lifestyle="giant_amf", architecture="very_large_amf_gypsy_copia",
        genome_scale="extreme_large",
        repeat_regime="extreme", te_regime="high",
        hgt_regime="moderate", expected_sv_bias="DUP_large_INS",
        # FIX: GC set to 28.5% (DAOM 197198). Genome ~150 Mb post-decontam.
        # Pre-decontamination 750 Mb figure was wrong. Chen 2018 PNAS.
    ),
    # Mucoromycota / soil endophyte
    "soil_endophyte": dict(
        # FIX: renamed from 'endomycorrhizal'; lifestyle corrected.
        # Mortierella is NOT endomycorrhizal. Guo et al. 2020 Curr. Biol.
        phylum="Mucoromycota", cls="Mortierellomycetes",
        order="Mortierellales", family="Mortierellaceae", genus="Mortierella",
        lifestyle="soil_endophyte", architecture="moderate_repeat_soil_endophyte",
        genome_scale="medium", repeat_regime="moderate", te_regime="moderate",
        hgt_regime="low", expected_sv_bias="DEL_INS",
    ),
    # Chytridiomycota
    "hgt_receiver": dict(
        phylum="Chytridiomycota", cls="Chytridiomycetes",
        order="Rhizophydiales",     # FIX: was Chytridiales. James 2006.
        family="Rhizophydiaceae",   # FIX: was Chytriaceae.
        genus="Batrachochytrium",
        lifestyle="hgt_receiver", architecture="gc_shifted_hgt_receiver",
        genome_scale="small_medium", repeat_regime="moderate", te_regime="low",
        hgt_regime="very_high", expected_sv_bias="TRA_HGT",
    ),
    "cross_phylum_hgt_stress": dict(
        phylum="Chytridiomycota", cls="Chytridiomycetes",
        order="Rhizophydiales",     # FIX: was Chytridiales.
        family="Rhizophydiaceae",   # FIX: was Chytriaceae.
        genus="Batrachochytrium",
        lifestyle="cross_phylum_hgt", architecture="cross_phylum_hgt_stress",
        genome_scale="medium", repeat_regime="moderate", te_regime="moderate",
        hgt_regime="extreme", expected_sv_bias="TRA_HGT_OFF_REF",
    ),
}

ALLOWED_SV_TYPES: frozenset[str] = frozenset({"INS", "DEL", "DUP", "INV", "TRA"})

# ── Element dispatch per scenario ─────────────────────────────────────────
_SCENARIO_ELEMENTS: dict[str, list[str]] = {
    "core":                        ["NONE"],
    "compact_yeast":               ["NONE"],
    "saprotrophic":                ["HGT", "TE_TIR"],
    "pathogenic":                  ["TE_LTR", "TE_TIR", "RIP"],
    "two_speed_pathogen_extreme":  ["TE_LTR", "RIP", "TE_TIR"],
    "necrotrophic":                ["TE_LINE", "TE_TIR"],
    "strict_endophyte":            ["NONE"],
    "lichenised":                  ["HGT"],   # FIX: algal HGT documented
    "ectomycorrhizal":             ["TE_TIR", "REPEAT"],
    "rust_smut_te_heavy":          ["TE_LTR", "RIP"],   # FIX: TIR removed
    "smut":                        ["NONE"],
    "arbuscular_mf":               ["TE_LTR", "REPEAT"], # FIX: STARSHIP removed
    "giant_amf":                   ["TE_LTR", "REPEAT"], # FIX: STARSHIP removed
    "soil_endophyte":              ["TE_LINE", "REPEAT"],
    "hgt_receiver":                ["HGT"],
    "cross_phylum_hgt_stress":     ["HGT", "TE_LINE"],
}


# ── Helpers ───────────────────────────────────────────────────────────────

def normalize_query_contig_name(name: str) -> str:
    """Strip __sv_ simulator suffix."""
    return name.split("__sv_", 1)[0]


def clamp01(value: float | None, default: float) -> float:
    if value is None:
        return default
    return max(0.0, min(1.0, value))


def scenario_gc(meta: dict[str, str], override: float | None) -> float:
    """Lookup GC from published genome values (genus → class → heuristic)."""
    if override is not None:
        return clamp01(override, 0.45)
    gc = _GENUS_GC.get(meta.get("genus", ""))
    if gc is not None:
        return gc
    gc = _CLASS_GC.get(meta.get("cls", ""))
    if gc is not None:
        return gc
    ls = meta.get("lifestyle", "")
    if ls in {"hgt_receiver", "cross_phylum_hgt"}: return 0.40
    if ls in {"compact_yeast", "baseline_control"}: return 0.38
    return 0.45


def off_ref_gc_for_scenario(gc: float) -> float:
    """Direction-aware off-reference GC (clamped to [0.20, 0.75]).

    Low-GC hosts (gc < 0.45): off-ref is higher-GC (bacterial donors).
    High-GC hosts (gc >= 0.45): off-ref may be AT-rich (TEs, Firmicutes).
    """
    shift = +0.18 if gc < 0.45 else -0.14
    return max(0.20, min(0.75, gc + shift))


def ensure_split(n_genomes: int, n_reps: int) -> tuple[int, int]:
    n_genomes = max(2, int(n_genomes))
    n_reps    = max(1, int(n_reps))
    if n_reps >= n_genomes:
        n_reps = n_genomes - 1
    return n_genomes, n_reps


def parse_biases(bias: str) -> list[str]:
    tokens = [t.upper() for t in bias.split("_") if t]
    types  = [t for t in tokens if t in ALLOWED_SV_TYPES]
    return types or ["INS"]


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


def write_fastq(path: Path, reads: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for name, seq in reads:
            fh.write(f"@{name}\n{seq}\n+\n" + ("I" * len(seq)) + "\n")


def mutate_read(seq: str, error_rate: float, seed: int) -> str:
    if error_rate <= 0.0 or not seq:
        return seq
    rng = random.Random(seed)
    bases = "ACGT"
    out: list[str] = []
    for b in seq:
        if rng.random() < error_rate:
            out.append(rng.choice([x for x in bases if x != b]))
        else:
            out.append(b)
    return "".join(out)


def mutate_sequence(seq: str, divergence_rate: float, seed: int) -> str:
    return mutate_read(seq, divergence_rate, seed)


def simulate_reads_for_records(records: list[tuple[str, str]], mode: str,
                               short_read_len: int, short_step: int, short_cov: int,
                               long_read_len: int, long_step: int, long_cov: int,
                               long_error_rate: float, seed_base: int) -> list[tuple[str, str]]:
    reads: list[tuple[str, str]] = []
    rid = 0
    for rec_idx, (name, seq) in enumerate(records):
        if mode == "short-reads":
            read_len = max(20, min(short_read_len, len(seq)))
            if len(seq) < read_len:
                continue
            for rep in range(max(1, short_cov)):
                for start in range(0, len(seq) - read_len + 1, max(1, short_step)):
                    reads.append((f"{name}_sr{rep}_{rid}", seq[start:start + read_len]))
                    rid += 1
        elif mode == "long-reads":
            read_len = max(80, min(long_read_len, len(seq)))
            if len(seq) < read_len:
                continue
            step = max(1, long_step)
            for rep in range(max(1, long_cov)):
                offset = (rep * max(1, step // 2)) % step
                for start in range(offset, len(seq) - read_len + 1, step):
                    frag = seq[start:start + read_len]
                    frag = mutate_read(frag, long_error_rate, seed_base + rec_idx * 100000 + rid)
                    reads.append((f"{name}_lr{rep}_{rid}", frag))
                    rid += 1
        else:
            raise ValueError(f"unsupported read simulation mode: {mode}")
    return reads


def write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(map(str, row)) + "\n")


def scenario_names(raw: str) -> list[str]:
    value = (raw or "core").strip()
    if value == "all":
        return list(SCENARIOS.keys())
    selected = [s.strip() for s in value.split(",") if s.strip()]
    selected = [s for s in selected if s in SCENARIOS]
    return selected or ["core"]


def rand_seq(n: int, gc: float = 0.45, seed: int = 1) -> str:
    gc = max(0.0, min(1.0, gc))
    at = 1.0 - gc
    weights = [gc / 2.0, gc / 2.0, at / 2.0, at / 2.0]
    bases = ["G", "C", "A", "T"]
    rnd = random.Random(seed)
    return "".join(rnd.choices(bases, weights=weights, k=n))


# ── Element embedders ─────────────────────────────────────────────────────

def embed_tandem_repeat(seq: str, period: int = 4, copies: int = 8,
                        pos: int | None = None, seed: int = 1) -> str:
    """Period-regularity tandem repeat (≥5 copies, ≥50 bp)."""
    rnd = random.Random(seed)
    pos = pos if pos is not None else rnd.randint(0, max(0, len(seq) // 4))
    unit = rand_seq(period, gc=0.45, seed=seed + 1)
    element = unit * copies
    end = min(pos + len(element), len(seq))
    return seq[:pos] + element[:end - pos] + seq[end:]


def embed_ltr_element(seq: str, ltr_len: int = 80, interior_len: int = 500,
                      pos: int | None = None, seed: int = 1) -> str:
    """LTR retrotransposon (Gypsy/Copia): DTR ≥50 bp + RT/integrase coding interior."""
    rnd = random.Random(seed)
    pos = pos if pos is not None else rnd.randint(0, max(0, len(seq) // 4))
    ltr      = rand_seq(ltr_len,      gc=0.52, seed=seed + 2)  # LTR ~52% GC
    interior = rand_seq(interior_len, gc=0.55, seed=seed + 3)  # coding ~55% GC
    element  = ltr + interior + ltr
    end = min(pos + len(element), len(seq))
    return seq[:pos] + element[:end - pos] + seq[end:]


def embed_tir_element(seq: str, tir_len: int = 35, interior_len: int = 300,
                      pos: int | None = None, seed: int = 1) -> str:
    """TIR DNA transposon (Tc1/Mariner or MITE): inverted terminal repeats ≥30 bp."""
    rnd = random.Random(seed)
    pos = pos if pos is not None else rnd.randint(0, max(0, len(seq) // 4))
    tir  = rand_seq(tir_len, gc=0.48, seed=seed + 4)
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    tir_rc   = "".join(comp.get(b, b) for b in reversed(tir))
    interior = rand_seq(interior_len, gc=0.45, seed=seed + 5)
    element  = tir + interior + tir_rc
    end = min(pos + len(element), len(seq))
    return seq[:pos] + element[:end - pos] + seq[end:]


def embed_line_helitron(seq: str, element_len: int = 400,
                        pos: int | None = None, seed: int = 1) -> str:
    """LINE/Helitron: AT-rich body (GC<0.42) + poly-A terminus."""
    rnd = random.Random(seed)
    pos = pos if pos is not None else rnd.randint(0, max(0, len(seq) // 4))
    core    = rand_seq(element_len - 20, gc=0.30, seed=seed + 6)
    poly_a  = "A" * 20
    element = core + poly_a
    end = min(pos + len(element), len(seq))
    return seq[:pos] + element[:end - pos] + seq[end:]


def embed_sine(seq: str, sine_len: int = 150,
               pos: int | None = None, seed: int = 1) -> str:
    """SINE: short (≤500 bp), high-GC 5' domain, terminal repeat."""
    rnd = random.Random(seed)
    pos = pos if pos is not None else rnd.randint(0, max(0, len(seq) // 4))
    tr       = rand_seq(15, gc=0.60, seed=seed + 7)
    interior = rand_seq(sine_len - 30, gc=0.56, seed=seed + 8)
    element  = tr + interior + tr
    end = min(pos + len(element), len(seq))
    return seq[:pos] + element[:end - pos] + seq[end:]


def embed_starship(seq: str, body_len: int = 3000, cargo_len: int = 1200,
                   clade_gc: float = 0.45, pos: int | None = None,
                   seed: int = 1) -> str:
    """Starship element (Ascomycota only).

    Hull GC = clade_gc − 0.12 (AT-rich relative to background).
    Cargo GC ≈ genic GC (clade_gc ± 0.02, clamped 0.42–0.55).
    Ref: Urquhart et al. 2023 Current Biology.

    NOT appropriate for Glomeromycota or Basidiomycota (not confirmed there).
    """
    rnd = random.Random(seed)
    pos = pos if pos is not None else rnd.randint(0, max(0, len(seq) // 5))
    hull_gc  = max(0.20, clade_gc - 0.12)
    cargo_gc = max(0.42, min(0.55, clade_gc + 0.02))
    body_l   = rand_seq(body_len // 2, gc=hull_gc,  seed=seed + 9)
    cargo    = rand_seq(cargo_len,     gc=cargo_gc, seed=seed + 10)
    body_r   = rand_seq(body_len // 2, gc=hull_gc,  seed=seed + 11)
    element  = body_l + cargo + body_r
    end = min(pos + len(element), len(seq))
    return seq[:pos] + element[:end - pos] + seq[end:]


def embed_hgt_island(seq: str, island_len: int = 600,
                     donor_gc: float = 0.55, pos: int | None = None,
                     seed: int = 1) -> str:
    """HGT island: GC-shifted window (±0.10 from host, published range ±0.05-0.10).

    Ref: Slot & Rokas 2011 Science; Marcet-Houben & Gabaldon 2010.
    """
    rnd = random.Random(seed)
    pos = pos if pos is not None else rnd.randint(0, max(0, len(seq) // 4))
    island = rand_seq(island_len, gc=max(0.20, min(0.75, donor_gc)), seed=seed + 12)
    end = min(pos + len(island), len(seq))
    return seq[:pos] + island[:end - pos] + seq[end:]


def apply_rip(seq: str, window: int = 500, fraction: float = 0.3,
              seed: int = 1) -> str:
    """Apply RIP (Repeat-Induced Point mutation): CpA → TpA transitions.

    RIP exclusively targets the CpA dinucleotide context on the forward strand.
    This is measured by the RIP product index (TpA/CpA > 1.5) and RIP substrate
    index (CpA/TpA < 0.7). Ref: Selker et al. 2003; Cambareri et al. 1989.

    Previous implementation used C→T on any C (incorrect; introduced GC-ratio
    distortion without the canonical dinucleotide specificity).
    """
    rnd = random.Random(seed)
    bases = list(seq)
    for i in range(len(bases) - 1):
        if bases[i] in ("C", "c") and bases[i + 1] in ("A", "a"):
            if rnd.random() < fraction:
                bases[i] = "T" if bases[i] == "C" else "t"
    return "".join(bases)


def embed_element_by_class(seq: str, ec: str, clade_gc: float = 0.45,
                           seed: int = 1) -> str:
    """Dispatch to the appropriate embedder for an ElementClass string."""
    if ec == "REPEAT":   return embed_tandem_repeat(seq, seed=seed)
    if ec == "TE_LTR":   return embed_ltr_element(seq, seed=seed)
    if ec == "TE_TIR":   return embed_tir_element(seq, seed=seed)
    if ec == "TE_LINE":  return embed_line_helitron(seq, seed=seed)
    if ec == "TE_SINE":  return embed_sine(seq, seed=seed)
    if ec == "STARSHIP": return embed_starship(seq, clade_gc=clade_gc, seed=seed)
    if ec == "HGT":
        donor = max(0.25, min(0.72, clade_gc + 0.10))  # ±0.10 deviation
        return embed_hgt_island(seq, donor_gc=donor, seed=seed)
    if ec == "RIP":      return apply_rip(seq, seed=seed)
    return seq  # NONE


# ── Long-read platform presets ────────────────────────────────────────────
# Applied when --long-read-platform is set; override individual --long-read-*
# flags so simulated reads match each platform's actual characteristics.
#
# PacBio HiFi CCS (Revio / Sequel IIe):
#   ≥Q20 per-read accuracy (error rate ≈ 0.1 %).  Reads 10–25 kb.
#   Downstream: minimap2 map-hifi → samtools sort+index → sniffles2 / cuteSV
#   (HiFi cluster params) / SVIM.
#
# ONT R10.4.1 standard simplex (PromethION / GridION / MinION Mk1C):
#   ~Q20 median accuracy (error rate ≈ 1 %).  Reads 10–30 kb, median ~15 kb.
#   Downstream: minimap2 map-ont → sniffles2 --long-read-model ont_r10_q20
#   (Sniffles2 ≥v2.2) / cuteSV / SVIM.
#   WhatsHap phase + haplotag is applicable for diploid / dikaryotic fungi
#   (Puccinia, Leptosphaeria, Zymoseptoria tritici) once SNP calls are made
#   via bcftools mpileup | bcftools call or Clair3 from the same BAM.
#
# ONT R9.4.1 (legacy):
#   ~Q15 median accuracy (error rate ≈ 5 %).  Still common in public ENA data.
#   Same minimap2 map-ont preset as R10.4.1; lower SV recall in repeat-rich
#   regions (relevant for AMF / two-speed pathogen scenarios).
#
# generic: use explicit --long-read-len / --long-read-error-rate / etc. flags.

_LR_PLATFORM_PRESETS: dict[str, dict[str, int | float]] = {
    "hifi":    {"len": 15000, "step": 5000, "cov": 15, "err": 0.001},
    "ont-r10": {"len": 10000, "step": 3000, "cov": 20, "err": 0.010},
    "ont-r9":  {"len":  8000, "step": 2500, "cov": 20, "err": 0.050},
}


def write_truth_vcf(path: Path, rows: list[list[str]]) -> None:
    """Write VCF4.3 truth file from truth TSV rows (15 fields each).

    Emits a *multi-sample* VCF: one column per distinct query assembly that
    contributed truth events. Each row has GT 1/1 only for the query asm that
    introduced the SV and 0/0 elsewhere. The QUERY_ASM info field is
    preserved for downstream readers (sv_pr_utils.parse_vcf_records pulls
    qasm from INFO, not from sample columns), so scoring keeps working.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    samples: list[str] = []
    seen: set[str] = set()
    for row in rows:
        if not row:
            continue
        asm = row[0]
        if asm and asm not in seen:
            seen.add(asm)
            samples.append(asm)
    if not samples:
        samples = ["SAMPLE"]
    sample_index = {asm: i for i, asm in enumerate(samples)}
    with open(path, "w") as out:
        out.write("##fileformat=VCFv4.3\n##source=test_amf_simulator\n")
        for tag, desc in [
            ("SVTYPE","SV type"), ("SVLEN","SV length"), ("END","End position"),
            ("SCENARIO","Simulation scenario"), ("QUERY_ASM","Query assembly"),
            ("PHYLUM","Phylum"), ("CLASS","Class"), ("ORDER","Order"),
            ("FAMILY","Family"), ("GENUS","Genus"),
            ("CHR2","Mate contig for TRA"), ("POS2","Mate position for TRA"),
            ("END2","Mate end for TRA"),
        ]:
            t = "Integer" if tag in {"SVLEN","END","POS2","END2"} else "String"
            n = "1"
            out.write(f'##INFO=<ID={tag},Number={n},Type={t},Description="{desc}">\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                  + "\t".join(samples) + "\n")
        for idx, row in enumerate(rows, start=1):
            if len(row) < 15:
                row = list(row) + [".", "0", "0", "0"][: 15 - len(row)]
            (asm, qcontig, svtype, pos, svlen, scenario,
             phylum, cls_name, order, family, genus,
             mate_contig, mate_pos, mate_end, _) = row[:15]
            cls = cls_name  # local alias; does not shadow outer 'cls' key
            chrom   = normalize_query_contig_name(qcontig)
            pos_i   = int(pos)
            svlen_i = int(svlen)
            end_i   = pos_i if svtype.startswith("TRA") else pos_i + max(svlen_i - 1, 0)
            info = (f"SVTYPE={svtype};SVLEN={svlen_i};END={end_i}"
                    f";SCENARIO={scenario};QUERY_ASM={asm}"
                    f";PHYLUM={phylum};CLASS={cls};ORDER={order}"
                    f";FAMILY={family};GENUS={genus}")
            if svtype.startswith("TRA") and mate_contig not in {"", "."}:
                mc = normalize_query_contig_name(mate_contig)
                mp = max(1, int(mate_pos))
                me = max(mp, int(mate_end))
                info += f";CHR2={mc};POS2={mp};END2={me}"
            owner = sample_index.get(asm, -1)
            gts = ["0/0"] * len(samples)
            if 0 <= owner < len(gts):
                gts[owner] = "1/1"
            out.write(f"{chrom}\t{pos_i}\ttruth{idx}\tN\t<{svtype}>\t60\tPASS"
                      f"\t{info}\tGT\t" + "\t".join(gts) + "\n")


# ── main ──────────────────────────────────────────────────────────────────

def main() -> None:  # noqa: C901
    ap = argparse.ArgumentParser(description="Fungal SV simulator — TOL pipeline")
    ap.add_argument("--phylum",                  default="Ascomycota")
    ap.add_argument("--n-genomes",   type=int,   default=4)
    ap.add_argument("--n-reps",      type=int,   default=2)
    ap.add_argument("--total-len",   type=int,   default=50000)
    ap.add_argument("--n-contigs",   type=int,   default=3)
    ap.add_argument("--out-dir",     required=True)
    ap.add_argument("--scenario-set",            default="core")
    ap.add_argument("--write-extended-manifest", action="store_true")
    ap.add_argument("--combined-annotations",    default="all")
    ap.add_argument("--te-rate",     type=float, default=0.05)
    ap.add_argument("--starship-rate", type=float, default=0.01)
    ap.add_argument("--hgt-rate",    type=float, default=0.005)
    ap.add_argument("--divergence",  type=float, default=0.01)
    ap.add_argument("--seed",        type=int,   default=1,
                    help="Compatibility seed for external wrappers; current simulator output is already deterministic.")
    ap.add_argument("--class-name",              default="")
    ap.add_argument("--order-name",              default="")
    ap.add_argument("--family-name",             default="")
    ap.add_argument("--genus-name",              default="")
    ap.add_argument("--repeat-density", type=float, default=None)
    ap.add_argument("--genome-size-scale", type=float, default=None)
    ap.add_argument("--gc",          type=float, default=None)
    ap.add_argument("--write-query-annotations",    action="store_true")
    ap.add_argument("--write-rep-truth-per-contig", action="store_true")
    ap.add_argument("--write-hint-contigs", action="store_true",
                    help="Encode SV metadata in contig names. Never use for benchmarking.")
    ap.add_argument("--write-query-truth", action="store_true",
                    help="Alias for compatibility; query_truth.tsv is always written.")
    ap.add_argument("--n-svs-per-contig", type=int, default=1,
                    help="Number of on-reference SVs to embed per query contig (currently fixed at 1).")
    ap.add_argument("--tol-query-window-bp", type=int, default=2000000,
                    help="Streaming window size in bp (passed through for provenance only).")
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--query-mode", default="assembly",
                    choices=["assembly", "short-reads", "long-reads", "auto"],
                    help="Type of query artifact to emit in query_list.txt.")
    ap.add_argument("--short-read-len", type=int, default=150)
    ap.add_argument("--short-read-step", type=int, default=30)
    ap.add_argument("--short-read-coverage", type=int, default=6)
    ap.add_argument("--long-read-len", type=int, default=1200)
    ap.add_argument("--long-read-step", type=int, default=300)
    ap.add_argument("--long-read-coverage", type=int, default=2)
    ap.add_argument("--long-read-error-rate", type=float, default=0.03)
    ap.add_argument("--long-read-platform", default="ont-r10",
                    choices=["hifi", "ont-r10", "ont-r9", "generic"],
                    help="Long-read sequencing platform preset.  Overrides individual "
                         "--long-read-* flags: hifi=PacBio HiFi CCS (≥Q20, 15 kb), "
                         "ont-r10=ONT R10.4.1 simplex (~Q20, 10 kb), "
                         "ont-r9=ONT R9.4.1 (~Q15, 8 kb), "
                         "generic=use explicit --long-read-* values.")
    args = ap.parse_args()

    # Apply platform preset before anything reads the long-read parameters.
    if args.long_read_platform != "generic":
        p = _LR_PLATFORM_PRESETS[args.long_read_platform]
        args.long_read_len        = int(p["len"])
        args.long_read_step       = int(p["step"])
        args.long_read_coverage   = int(p["cov"])
        args.long_read_error_rate = float(p["err"])

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "truth").mkdir(exist_ok=True)

    scenarios  = scenario_names(args.scenario_set)
    n_genomes, n_reps = ensure_split(args.n_genomes, args.n_reps)
    # With the round-robin assignment scen=scenarios[i%N], query genomes start
    # at index n_reps.  If n_genomes-n_reps < len(scenarios) some scenarios
    # only ever become references and produce zero truth SVs.  Bump n_genomes
    # so every scenario gets at least one query slot.
    if n_genomes - n_reps < len(scenarios):
        n_genomes = n_reps + len(scenarios)
    seq_len    = max(2000, int(args.total_len) // max(1, int(args.n_contigs)))
    query_emit_mode = "assembly" if args.query_mode == "auto" else args.query_mode

    refs:                  list[str]       = []
    queries:               list[str]       = []
    rep_asms:              list[str]       = []
    base_rows:             list[list[str]] = []
    hier_rows:             list[list[str]] = []
    meta_rows:             list[list[str]] = []
    stress_rows:           list[list[str]] = []
    query_truth_rows:      list[list[str]] = []
    query_annotation_rows: list[list[str]] = []
    rep_truth_rows:        list[list[str]] = []

    off_ref_seed_base = 900_000
    scenario_backbones: dict[tuple[str, float, int], str] = {}

    for i in range(n_genomes):
        scen = scenarios[i % len(scenarios)]
        m    = dict(SCENARIOS[scen])

        if args.class_name:  m["cls"]    = args.class_name
        if args.order_name:  m["order"]  = args.order_name
        if args.family_name: m["family"] = args.family_name
        if args.genus_name:  m["genus"]  = args.genus_name
        if args.phylum:      m["phylum"] = args.phylum

        asm   = f"{scen}_asm{i + 1}"
        fasta = out / f"{asm}.fa"
        records: list[tuple[str, str]] = []
        gc = scenario_gc(m, args.gc)
        # giant_amf uses DAOM 197198 GC (28.5%) not the standard Rhizophagus 32%
        if scen == "giant_amf":
            gc = 0.285

        elem_classes = _SCENARIO_ELEMENTS.get(scen, ["NONE"])
        homologous_bases: list[str] = []
        for c in range(int(args.n_contigs)):
            backbone_key = (scen, gc, c)
            backbone = scenario_backbones.get(backbone_key)
            if backbone is None:
                backbone = rand_seq(
                    seq_len,
                    gc=gc,
                    seed=1000 + scenarios.index(scen) * 131 + c,
                )
                scenario_backbones[backbone_key] = backbone
            homologous_bases.append(
                mutate_sequence(backbone, args.divergence, seed=40_000 + i * 313 + c)
            )

        for c in range(int(args.n_contigs)):
            base = homologous_bases[c]
            name = f"ctg{c + 1}"

            if i >= n_reps and c == 0:
                orc_gc = off_ref_gc_for_scenario(gc)
                base   = rand_seq(len(base), gc=orc_gc,
                                  seed=off_ref_seed_base + i * 31 + c)
                ec = elem_classes[0] if elem_classes else "NONE"
                if ec != "NONE":
                    base = embed_element_by_class(base, ec, clade_gc=gc,
                                                  seed=off_ref_seed_base + i * 31 + c + 500)
                name = (f"ctg1__sv_OFF_REF__pos__1__len__{len(base)}"
                        if args.write_hint_contigs else "ctg1")
                query_truth_rows.append([
                    asm, name, "OFF_REF", "1", str(len(base)), scen,
                    m["phylum"], m["cls"], m["order"], m["family"], m["genus"],
                    ".", "0", "0", "0",
                ])
                if args.write_query_annotations:
                    query_annotation_rows.append(
                        [asm, "ctg1", ec if ec != "NONE" else "HGT",
                         m["architecture"], str(len(base))]
                    )

            elif i >= n_reps and c > 0:
                sv_types = parse_biases(m.get("expected_sv_bias", "INS"))
                rng      = random.Random(5000 + i * 131 + c)
                svtype   = rng.choice(sv_types)
                max_var  = max(50, min(len(base) // 10, 500))
                var_len  = max(50, min(max_var, len(base) // 5))
                pos_0    = rng.randint(0, max(0, len(base) - var_len))
                ref_pos  = pos_0 + 1
                mate_contig = "."; mate_pos = "0"; mate_end = "0"; mate_svlen = "0"

                if svtype == "DEL":
                    base = base[:pos_0] + base[pos_0 + var_len:]
                elif svtype == "INS":
                    ins  = rand_seq(var_len, gc=gc, seed=2000 + i * 31 + c)
                    base = base[:pos_0] + ins + base[pos_0:]
                elif svtype == "DUP":
                    seg  = base[pos_0 : pos_0 + var_len]
                    base = base[:pos_0] + seg + base[pos_0:]
                elif svtype == "INV":
                    seg  = base[pos_0 : pos_0 + var_len]
                    comp = {"A":"T","C":"G","G":"C","T":"A"}
                    revcomp = "".join(comp.get(b, b) for b in reversed(seg))
                    base = base[:pos_0] + revcomp + base[pos_0 + var_len:]
                elif svtype == "TRA":
                    mate_choices = [x for x in range(int(args.n_contigs)) if x != c]
                    mate_idx  = rng.choice(mate_choices) if mate_choices else c
                    mate_pos0 = rng.randint(0, max(0, len(base) - var_len))
                    donor     = homologous_bases[mate_idx]
                    base = base[:pos_0] + donor[mate_pos0:mate_pos0+var_len] + base[pos_0+var_len:]
                    mate_contig = f"ctg{mate_idx + 1}"
                    mate_pos    = str(mate_pos0 + 1)
                    mate_end    = str(mate_pos0 + var_len)
                    mate_svlen  = str(var_len)

                if elem_classes[0] not in {"NONE", "RIP"} and len(base) > 200:
                    base = embed_element_by_class(base, elem_classes[0], clade_gc=gc,
                                                  seed=3000 + i * 131 + c)

                if args.write_hint_contigs:
                    name = f"ctg{c+1}__sv_{svtype}__pos__{ref_pos}__len__{var_len}"
                    if svtype == "TRA":
                        name += (f"__mate_contig__{mate_contig}"
                                 f"__mate_pos__{mate_pos}__mate_len__{mate_svlen}")
                else:
                    name = f"ctg{c + 1}"

                query_truth_rows.append([
                    asm, name, svtype, str(ref_pos), str(var_len), scen,
                    m["phylum"], m["cls"], m["order"], m["family"], m["genus"],
                    mate_contig, mate_pos, mate_end, mate_svlen,
                ])
                if args.write_query_annotations:
                    query_annotation_rows.append(
                        [asm, name, svtype, m["architecture"], str(var_len)])

            else:
                if i < n_reps and args.write_rep_truth_per_contig:
                    rep_truth_rows.append([asm, name, "REFERENCE", "0", "0", scen])

            records.append((name, base))

        write_fasta(fasta, records)
        if i < n_reps:
            refs.append(str(fasta)); rep_asms.append(str(fasta))
        else:
            query_path = fasta
            if query_emit_mode in {"short-reads", "long-reads"}:
                reads = simulate_reads_for_records(
                    records, query_emit_mode,
                    args.short_read_len, args.short_read_step, args.short_read_coverage,
                    args.long_read_len, args.long_read_step, args.long_read_coverage,
                    args.long_read_error_rate,
                    seed_base=700000 + i * 1000,
                )
                suffix = ".fq" if query_emit_mode == "short-reads" else ".fastq"
                query_path = out / f"{asm}{suffix}"
                write_fastq(query_path, reads)
            queries.append(str(query_path))
            meta_rows.append(
                [asm, scen, m["phylum"], m["cls"], m["order"], m["family"], m["genus"]]
            )

        species_clade = f"{m['genus']}_{asm}"
        base_rows.append([asm, m["phylum"], m["cls"], m["order"], m["family"],
                          m["genus"], species_clade, "species", str(fasta)])
        rank_to_clade: dict[str, str] = {
            "phylum": m["phylum"], "class": m["cls"], "order": m["order"],
            "family": m["family"], "genus": m["genus"], "species": species_clade,
        }
        for rank in RANKS:
            hier_rows.append([asm, m["phylum"], m["cls"], m["order"], m["family"],
                               m["genus"], rank_to_clade[rank], rank, str(fasta)])
        stress_rows.append([scen, m["lifestyle"], m["architecture"], m["genome_scale"],
                            m["repeat_regime"], m["te_regime"], m["hgt_regime"],
                            m["expected_sv_bias"]])

    if not refs or not queries:
        raise SystemExit("simulator produced zero refs or zero queries")

    (out / "ref_list.txt").write_text("\n".join(refs) + "\n")
    (out / "query_list.txt").write_text("\n".join(queries) + "\n")
    (out / "rep_asm_list.txt").write_text("\n".join(rep_asms) + "\n")

    HDR = ["#asm_name","phylum","class","order","family","genus",
           "clade_name","clade_rank","fasta_path"]
    write_tsv(out / "base_manifest.tsv",      HDR, base_rows)
    write_tsv(out / "hierarchy_manifest.tsv", HDR, hier_rows)
    write_tsv(out / "query_metadata.tsv",
              ["query_asm","scenario","phylum","class","order","family","genus"],
              meta_rows)

    unique_stress: list[list[str]] = []
    seen_keys: set[tuple[str, ...]] = set()
    for row in stress_rows:
        key = tuple(row)
        if key not in seen_keys:
            seen_keys.add(key); unique_stress.append(row)
    write_tsv(out / "stress_case_catalog.tsv",
              ["scenario","lifestyle","architecture","genome_scale",
               "repeat_regime","te_regime","hgt_regime","expected_sv_bias"],
              unique_stress)
    write_tsv(out / "simulation_params.tsv",
              ["n_genomes","n_reps","n_queries","scenario_count",
               "divergence","n_svs_per_contig","window_bp","write_hint_contigs","query_mode",
               "lr_platform","lr_read_len","lr_error_rate","lr_coverage"],
              [[str(n_genomes), str(n_reps), str(len(queries)), str(len(scenarios)),
                str(args.divergence), "1", str(seq_len), str(int(args.write_hint_contigs)), query_emit_mode,
                args.long_read_platform, str(args.long_read_len),
                str(args.long_read_error_rate), str(args.long_read_coverage)]])

    write_tsv(out / "query_truth.tsv",
              ["query_asm","query_contig","svtype","pos","svlen","scenario",
               "phylum","class","order","family","genus",
               "mate_contig","mate_pos","mate_end","mate_svlen"],
              query_truth_rows)
    write_truth_vcf(out / "truth" / "all_queries.truth.ref.vcf", query_truth_rows)

    if args.write_query_annotations:
        write_tsv(out / "graph_annotations_denovo.tsv",
                  ["query_asm","query_contig","annotation","architecture","length"],
                  query_annotation_rows)
    if args.write_rep_truth_per_contig:
        write_tsv(out / "rep_truth_per_contig.tsv",
                  ["ref_asm","ref_contig","label","pos","svlen","scenario"],
                  rep_truth_rows)


if __name__ == "__main__":
    main()
