#!/usr/bin/env python3
# Designed for Linux

from __future__ import annotations

import argparse
import csv
import ctypes
import gzip
import io
import json
import os
import shutil
import subprocess
import sys
import time
import hashlib
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from run_million_mode_query_benchmark import augment_routing_store
from sv_pr_utils import DEFAULT_TOL_BP, DEFAULT_TOL_LEN_FRAC, wilson_ci


ROOT = Path(__file__).resolve().parent
DEFAULT_BIN = ROOT / "fungi_graphsv_tol_bin"
DEFAULT_ANALYZE = ROOT / "analyze_new_biology_candidates.py"
MYCOSV_BRIDGE_CPP = ROOT / "mycosv_cli_bridge.cpp"

NCBI_ASSEMBLY_SUMMARY = {
    "ncbi-refseq": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt",
    "ncbi-genbank": "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/assembly_summary.txt",
}

PUBLIC_RESOURCE_LINKS: list[dict[str, str]] = [
    {
        "label": "ncbi_refseq_fungi_assembly_summary",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt",
        "description": "NCBI RefSeq fungal assembly catalog",
    },
    {
        "label": "ncbi_genbank_fungi_assembly_summary",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/assembly_summary.txt",
        "description": "NCBI GenBank fungal assembly catalog",
    },
    {
        "label": "ena_filereport_docs",
        "url": "https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/file-reports.html",
        "description": "ENA filereport API documentation",
    },
    {
        "label": "ena_portal_api",
        "url": "https://www.ebi.ac.uk/ena/portal/api/",
        "description": "ENA portal API base",
    },
    {
        "label": "ena_sra_ftp_docs",
        "url": "https://ena-docs.readthedocs.io/en/latest/retrieval/file-download/sra-ftp-structure.html",
        "description": "ENA SRA FTP structure for public reads",
    },
    {
        "label": "ensembl_fungi",
        "url": "https://fungi.ensembl.org/index.html",
        "description": "Ensembl Fungi portal and public dumps",
    },
    {
        "label": "mycocosm",
        "url": "https://mycocosm.jgi.doe.gov/",
        "description": "JGI MycoCosm fungal genome portal",
    },
]

ENA_FILEREPORT_FIELDS = [
    "run_accession",
    "scientific_name",
    "study_accession",
    "sample_accession",
    "experiment_accession",
    "instrument_platform",
    "instrument_model",
    "library_layout",
    "library_source",
    "library_strategy",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
    "submitted_ftp",
    "submitted_md5",
    "submitted_bytes",
]

PANEL_PRESETS: dict[str, list[dict[str, str]]] = {
    "compact_yeast": [
        {"species": "Saccharomyces cerevisiae", "scenario": "compact_yeast", "lifestyle": "yeast", "architecture": "compact"},
        {"species": "Candida glabrata", "scenario": "compact_yeast", "lifestyle": "yeast_pathogen", "architecture": "compact"},
        {"species": "Lachancea kluyveri", "scenario": "compact_yeast", "lifestyle": "yeast", "architecture": "compact"},
    ],
    "amf_large": [
        {"species": "Rhizophagus irregularis", "scenario": "amf_large", "lifestyle": "arbuscular_mycorrhizal", "architecture": "large_repeat_rich"},
        {"species": "Gigaspora rosea", "scenario": "amf_large", "lifestyle": "arbuscular_mycorrhizal", "architecture": "large_repeat_rich"},
    ],
    "te_rich_pathogen": [
        {"species": "Puccinia graminis", "scenario": "te_rich_pathogen", "lifestyle": "plant_pathogen", "architecture": "te_rich"},
        {"species": "Puccinia striiformis", "scenario": "te_rich_pathogen", "lifestyle": "plant_pathogen", "architecture": "te_rich"},
        {"species": "Ustilago maydis", "scenario": "te_rich_pathogen", "lifestyle": "plant_pathogen", "architecture": "smut_pathogen"},
    ],
    "two_speed_pathogen": [
        {"species": "Leptosphaeria maculans", "scenario": "two_speed_pathogen", "lifestyle": "plant_pathogen", "architecture": "two_speed"},
        {"species": "Zymoseptoria tritici", "scenario": "two_speed_pathogen", "lifestyle": "plant_pathogen", "architecture": "two_speed"},
        {"species": "Fusarium oxysporum", "scenario": "two_speed_pathogen", "lifestyle": "plant_pathogen", "architecture": "two_speed"},
    ],
    "cross_phylum_hgt": [
        {"species": "Aspergillus fumigatus", "scenario": "cross_phylum_hgt", "lifestyle": "saprotroph_pathogen", "architecture": "aspergillus"},
        {"species": "Cryptococcus neoformans", "scenario": "cross_phylum_hgt", "lifestyle": "yeast_pathogen", "architecture": "basidiomycete"},
        {"species": "Rhizophagus irregularis", "scenario": "cross_phylum_hgt", "lifestyle": "symbiont", "architecture": "amf_large"},
    ],
}

TYPE_CANON = {
    "INS": "INS",
    "DEL": "DEL",
    "INV": "INV",
    "TRA": "TRA",
    "TRANS": "TRA",
    "TRANSLOCATION": "TRA",
    "INVTR": "TRA",
    "BND": "TRA",          # SVIM/Sniffles/Delly/Manta emit BND for translocations
    "DUP": "DUP",
    "DUP:TANDEM": "DUP",
    "DUP:INT": "DUP",
    "INVDP": "DUP",
    "CNV": "DUP",          # cuteSV sometimes emits CNV for duplications
    "OFF_REF": "OFF_REF",
}


@dataclass
class NormalizedCall:
    query_asm: str
    query_contig: str
    pos: int
    end: int
    svtype: str
    svlen: int
    source: str
    coord_space: str = "query"
    annotation: str = "."
    element_class: str = "NONE"
    ref_asm: str = "."
    ref_contig: str = "."


def run(cmd: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        text=True,
        capture_output=True,
        check=True,
    )


def maybe_add_gpp_dll_dir() -> Any | None:
    if os.name != "nt":
        return None
    gpp = shutil.which("g++")
    if not gpp:
        return None
    try:
        return os.add_dll_directory(str(Path(gpp).resolve().parent))
    except (AttributeError, FileNotFoundError, OSError):
        return None


def current_bridge_dll_path() -> Path:
    digest = hashlib.sha1()
    source_paths = sorted({
        *ROOT.glob("*.cpp"),
        *ROOT.glob("*.hpp"),
    })
    for path in source_paths:
        if path.exists():
            digest.update(str(path.resolve()).encode("utf-8"))
            digest.update(str(path.stat().st_mtime_ns).encode("utf-8"))
    return ROOT / f"mycosv_cli_bridge_{digest.hexdigest()[:12]}.dll"


def ensure_mycosv_bridge_dll(force: bool = False) -> Path:
    dll_path = current_bridge_dll_path()
    need_build = force or not dll_path.exists()
    if not need_build:
        return dll_path
    run(
        [
            "g++",
            "-O2",
            "-std=c++17",
            "-pthread",
            "-shared",
            "-I",
            str(ROOT),
            str(MYCOSV_BRIDGE_CPP),
            "-o",
            str(dll_path),
        ],
        cwd=ROOT,
    )
    return dll_path


def run_mycosv_via_dll(argv: list[str]) -> subprocess.CompletedProcess[str]:
    dll_path = ensure_mycosv_bridge_dll()
    handle = maybe_add_gpp_dll_dir()
    try:
        lib = ctypes.CDLL(str(dll_path))
        func = lib.run_mycosv_cli
        func.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_char_p)]
        func.restype = ctypes.c_int
        encoded = [arg.encode("utf-8") for arg in argv]
        arr = (ctypes.c_char_p * len(encoded))(*encoded)
        rc = func(len(encoded), arr)
    finally:
        if handle is not None:
            handle.close()
    if rc != 0:
        raise subprocess.CalledProcessError(rc, argv)
    return subprocess.CompletedProcess(argv, rc, "", "")


def run_mycosv_command(cmd: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    return run(cmd, cwd=cwd)


def list_panels_text() -> str:
    lines = []
    for name, entries in sorted(PANEL_PRESETS.items()):
        species = ", ".join(item["species"] for item in entries)
        lines.append(f"{name}\t{species}")
    return "\n".join(lines)


def normalize_name(raw: str) -> str:
    out = []
    for ch in raw.strip():
        if ch.isalnum():
            out.append(ch)
        else:
            out.append("_")
    cleaned = "".join(out).strip("_")
    while "__" in cleaned:
        cleaned = cleaned.replace("__", "_")
    return cleaned or "unknown"


def tool_path(name: str) -> str | None:
    return shutil.which(name)


def split_values(raw: str) -> list[str]:
    norm = raw.replace(";", ",").replace("|", ",")
    return [item.strip() for item in norm.split(",") if item.strip()]


def looks_like_url(raw: str) -> bool:
    return raw.startswith(("http://", "https://", "ftp://"))


def normalise_download_url(raw: str) -> str:
    if raw.startswith("ftp.sra.ebi.ac.uk/"):
        return "https://" + raw
    if raw.startswith("ftp://ftp.sra.ebi.ac.uk/"):
        return "https://" + raw[len("ftp://"):]
    if raw.startswith("ftp://ftp.ncbi.nlm.nih.gov/"):
        return "https://" + raw[len("ftp://"):]
    return raw


def http_get_text(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "MycoSV-real-benchmark/1.0"})
    with urllib.request.urlopen(req) as resp:
        return resp.read().decode("utf-8")


def http_download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        return dest
    req = urllib.request.Request(url, headers={"User-Agent": "MycoSV-real-benchmark/1.0"})
    with urllib.request.urlopen(req) as resp, dest.open("wb") as out:
        shutil.copyfileobj(resp, out)
    return dest


def maybe_gunzip(path: Path, keep_gz: bool = True) -> Path:
    if path.suffix != ".gz":
        return path
    out = path.with_suffix("")
    if out.exists():
        return out
    with gzip.open(path, "rb") as src, out.open("wb") as dst:
        shutil.copyfileobj(src, dst)
    if not keep_gz:
        path.unlink()
    return out


def open_text_auto(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def write_public_resource_links(path: Path) -> None:
    rows = [dict(row) for row in PUBLIC_RESOURCE_LINKS]
    write_tsv(path, rows, ["label", "url", "description"])


def parse_assembly_summary(text: str) -> list[dict[str, str]]:
    header: list[str] | None = None
    rows: list[dict[str, str]] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        if line.startswith("#"):
            stripped = line.lstrip("#").strip()
            if stripped.startswith("assembly_accession"):
                header = stripped.split("\t")
            continue
        if header is None:
            raise ValueError("assembly_summary header not found")
        fields = line.split("\t")
        if len(fields) < len(header):
            fields.extend([""] * (len(header) - len(fields)))
        rows.append(dict(zip(header, fields)))
    return rows


def assembly_level_rank(level: str) -> int:
    level = (level or "").lower()
    if level == "complete genome":
        return 4
    if level == "chromosome":
        return 3
    if level == "scaffold":
        return 2
    if level == "contig":
        return 1
    return 0


def row_quality_key(row: dict[str, str]) -> tuple[int, int, int, int, str]:
    refseq_category = (row.get("refseq_category") or "").lower()
    refscore = 0
    if refseq_category == "reference genome":
        refscore = 3
    elif refseq_category == "representative genome":
        refscore = 2
    elif refseq_category not in {"", "na"}:
        refscore = 1
    latest = 1 if (row.get("version_status") or "").lower() == "latest" else 0
    complete = 1 if (row.get("genome_rep") or "").lower() == "full" else 0
    release_date = row.get("seq_rel_date") or ""
    return (
        refscore,
        latest,
        assembly_level_rank(row.get("assembly_level", "")),
        complete,
        release_date,
    )


def species_label_for_row(row: dict[str, str]) -> str:
    explicit = (row.get("species_label") or "").strip()
    if explicit:
        return explicit
    organism = (row.get("organism_name") or "").strip()
    if not organism:
        return "."
    words = organism.split()
    if len(words) >= 2:
        return f"{words[0]} {words[1]}"
    return organism


def species_group_key(row: dict[str, str]) -> str:
    species_taxid = (row.get("species_taxid") or "").strip()
    if species_taxid:
        return f"taxid:{species_taxid}"
    return species_label_for_row(row).lower()


def match_species(row: dict[str, str], species: str) -> bool:
    organism = (row.get("organism_name") or "").lower()
    target = species.strip().lower()
    return organism.startswith(target) or f" {target} " in f" {organism} "


def select_species_rows(rows: list[dict[str, str]], species: str, max_n: int) -> list[dict[str, str]]:
    matched = [
        row for row in rows
        if row.get("ftp_path", "").startswith("https://ftp.ncbi.nlm.nih.gov/")
        and row.get("assembly_accession")
        and match_species(row, species)
    ]
    matched.sort(key=row_quality_key, reverse=True)
    return matched[:max_n]


def select_all_public_rows(
    rows: list[dict[str, str]],
    min_assembly_level: str,
    latest_only: bool,
    max_total: int,
) -> list[dict[str, str]]:
    min_rank = assembly_level_rank(min_assembly_level)
    selected: list[dict[str, str]] = []
    for row in rows:
        if not row.get("ftp_path", "").startswith("https://ftp.ncbi.nlm.nih.gov/"):
            continue
        if not row.get("assembly_accession"):
            continue
        if assembly_level_rank(row.get("assembly_level", "")) < min_rank:
            continue
        if latest_only and (row.get("version_status") or "").lower() != "latest":
            continue
        selected.append(row)
    selected.sort(key=lambda r: (species_group_key(r), row_quality_key(r)), reverse=True)
    if max_total > 0:
        return selected[:max_total]
    return selected


def fetch_taxonomy_lineages(taxids: list[str], cache_path: Path | None = None) -> dict[str, dict[str, str]]:
    cache: dict[str, dict[str, str]] = {}
    if cache_path and cache_path.exists():
        cache = json.loads(cache_path.read_text())
    wanted = [taxid for taxid in taxids if taxid and taxid not in cache]
    if not wanted:
        return cache

    for start in range(0, len(wanted), 100):
        batch = wanted[start:start + 100]
        query = urllib.parse.urlencode({
            "db": "taxonomy",
            "id": ",".join(batch),
            "retmode": "xml",
        })
        xml_text = http_get_text(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{query}")
        root = ET.fromstring(xml_text)
        for taxon in root.findall(".//Taxon"):
            taxid = taxon.findtext("TaxId", default="")
            lineage = {"phylum": ".", "class": ".", "order": ".", "family": ".", "genus": ".", "species": "."}
            lineage["species"] = taxon.findtext("ScientificName", default=".") or "."
            for anc in taxon.findall("./LineageEx/Taxon"):
                rank = (anc.findtext("Rank", default="") or "").lower()
                sci = anc.findtext("ScientificName", default=".") or "."
                if rank in lineage:
                    lineage[rank] = sci
            if lineage["genus"] == "." and lineage["species"] not in {".", ""}:
                lineage["genus"] = lineage["species"].split()[0]
            cache[taxid] = lineage
        time.sleep(0.34)

    if cache_path:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(cache, indent=2, sort_keys=True))
    return cache


def ncbi_download_targets(row: dict[str, str], include_gff: bool) -> list[tuple[str, str]]:
    ftp_path = row["ftp_path"]
    stem = ftp_path.rstrip("/").split("/")[-1]
    targets = [(f"{ftp_path}/{stem}_genomic.fna.gz", f"{stem}_genomic.fna.gz")]
    if include_gff:
        targets.append((f"{ftp_path}/{stem}_genomic.gff.gz", f"{stem}_genomic.gff.gz"))
    return targets


def materialize_entry(url_or_path: str, dest: Path, keep_gz: bool = True) -> Path:
    if url_or_path.startswith(("http://", "https://", "ftp://")):
        path = http_download(url_or_path, dest)
    else:
        src = Path(url_or_path)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
        path = dest
    if path.suffix == ".gz":
        return maybe_gunzip(path, keep_gz=keep_gz)
    return path


def parse_custom_url_manifest(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows.append({k: (v or "").strip() for k, v in row.items()})
    return rows


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_ena_filereport_text(text: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    if not text.strip():
        return rows
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    for row in reader:
        rows.append({k: (v or "").strip() for k, v in row.items()})
    return rows


def ena_filereport_url(accession: str) -> str:
    query = urllib.parse.urlencode({
        "accession": accession,
        "result": "read_run",
        "fields": ",".join(ENA_FILEREPORT_FIELDS),
        "format": "tsv",
        "download": "false",
        "limit": "0",
    })
    return f"https://www.ebi.ac.uk/ena/portal/api/filereport?{query}"


def fetch_ena_read_runs(accession: str) -> list[dict[str, str]]:
    return parse_ena_filereport_text(http_get_text(ena_filereport_url(accession)))


def ena_filereport_species_url(species: str, max_rows: int = 200) -> str:
    """Build an ENA filereport URL that returns public read runs for a species.

    Used by prepare_from_ncbi when --query-mode is short-reads or long-reads
    and the panel presets only describe a species (no read accession). The
    scientific_name filter works against any rank ENA stores, so a genus name
    like 'Rhizophagus' also resolves — which is important for AMF where
    species-level assignments are patchy.
    """
    query = urllib.parse.urlencode({
        # ENA portal needs this quoted exactly as: scientific_name="X"
        "query": f'scientific_name="{species}"',
        "result": "read_run",
        "fields": ",".join(ENA_FILEREPORT_FIELDS),
        "format": "tsv",
        "download": "false",
        "limit": str(max_rows),
    })
    return f"https://www.ebi.ac.uk/ena/portal/api/filereport?{query}"


def fetch_ena_read_runs_by_species(species: str, max_rows: int = 200) -> list[dict[str, str]]:
    try:
        return parse_ena_filereport_text(http_get_text(ena_filereport_species_url(species, max_rows)))
    except Exception:
        # ENA returns 204/404 for unknown species; treat as "no runs" rather
        # than aborting the whole panel preparation.
        return []


def sequence_kind_from_name(name: str) -> str:
    lower = name.lower()
    if lower.endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq")):
        return "fastq"
    if lower.endswith((".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".fa", ".fna")):
        return "fasta"
    return "fastq"


def preferred_platforms_for_mode(query_mode: str) -> set[str]:
    mode = query_mode.strip().lower()
    if mode == "short-reads":
        return {"ILLUMINA", "BGISEQ", "DNBSEQ", "ION_TORRENT"}
    if mode == "long-reads":
        return {"OXFORD_NANOPORE", "PACBIO_SMRT"}
    return set()


def filter_ena_rows_for_mode(rows: list[dict[str, str]], query_mode: str) -> list[dict[str, str]]:
    preferred = preferred_platforms_for_mode(query_mode)
    if not preferred:
        return rows
    matched = [row for row in rows if (row.get("instrument_platform") or "").upper() in preferred]
    return matched or rows


def direct_read_sources_from_row(row: dict[str, str]) -> list[str]:
    sources: list[str] = []
    for key in ("fastq_url_1", "fastq_url_2", "fastq_url", "fastq_urls", "read_url", "read_urls", "path", "url"):
        raw = row.get(key, "")
        if not raw:
            continue
        if key in {"path", "url"} and row.get("query_mode", "assembly") == "assembly":
            continue
        for item in split_values(raw):
            sources.append(item)
    deduped: list[str] = []
    seen: set[str] = set()
    for item in sources:
        if item not in seen:
            deduped.append(item)
            seen.add(item)
    return deduped


def select_ena_read_sources(run_rows: list[dict[str, str]], query_mode: str, max_runs: int) -> tuple[list[str], list[dict[str, str]]]:
    filtered = filter_ena_rows_for_mode(run_rows, query_mode)
    urls: list[str] = []
    meta_rows: list[dict[str, str]] = []
    picked_runs = 0
    for row in filtered:
        row_urls = split_values(row.get("fastq_ftp", ""))
        if not row_urls:
            submitted = [
                item for item in split_values(row.get("submitted_ftp", ""))
                if sequence_kind_from_name(item) in {"fastq", "fasta"}
            ]
            row_urls = submitted
        if not row_urls:
            continue
        picked_runs += 1
        for item in row_urls:
            urls.append(normalise_download_url(item if looks_like_url(item) else f"https://{item}"))
        meta_rows.append({
            "run_accession": row.get("run_accession", "."),
            "study_accession": row.get("study_accession", "."),
            "sample_accession": row.get("sample_accession", "."),
            "scientific_name": row.get("scientific_name", "."),
            "instrument_platform": row.get("instrument_platform", "."),
            "library_layout": row.get("library_layout", "."),
            "library_strategy": row.get("library_strategy", "."),
            "source_url": ena_filereport_url(row.get("run_accession") or row.get("study_accession") or row.get("sample_accession") or "."),
        })
        if max_runs > 0 and picked_runs >= max_runs:
            break
    return urls, meta_rows


def merge_sequence_sources(sources: list[str], dest_prefix: Path) -> Path:
    if not sources:
        raise ValueError("No sequence sources were provided")
    kind = sequence_kind_from_name(sources[0])
    suffix = ".fastq" if kind == "fastq" else ".fasta"
    out_path = dest_prefix.with_suffix(suffix)
    if out_path.exists():
        return out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as out_fh:
        for idx, src in enumerate(sources):
            parsed = urllib.parse.urlparse(src)
            source_name = Path(parsed.path).name or f"source_{idx + 1}{suffix}"
            basename = f"part_{idx + 1}_{source_name}"
            local_part = materialize_entry(src, dest_prefix.parent / basename, keep_gz=False)
            part_kind = sequence_kind_from_name(local_part.name)
            if part_kind != kind:
                raise ValueError(f"Mixed sequence kinds in one query input: {kind} vs {part_kind} from {src}")
            with open_text_auto(local_part) as in_fh:
                shutil.copyfileobj(in_fh, out_fh)
            if local_part != out_path and local_part.parent == dest_prefix.parent and local_part.exists():
                local_part.unlink()
    return out_path


def materialize_query_input(
    row: dict[str, str],
    dest_dir: Path,
    default_source: str,
    default_benchmark: dict[str, str] | None = None,
    public_max_runs: int = 2,
) -> tuple[dict[str, str], str, list[dict[str, str]]]:
    asm_name = row.get("asm_name") or normalize_name(row.get("scientific_name") or row.get("clade_name") or row.get("species") or "query")
    query_mode = row.get("query_mode") or "assembly"
    benchmark = default_benchmark or {}
    source_rows: list[dict[str, str]] = []

    if query_mode == "assembly":
        fasta_src = row.get("fasta_url") or row.get("path") or row.get("url") or ""
        if not fasta_src:
            raise ValueError(f"Missing fasta_url/path/url for assembly query row {asm_name}")
        dest_name = Path(urllib.parse.urlparse(fasta_src).path).name or f"{asm_name}.fa"
        local_path = materialize_entry(fasta_src, dest_dir / dest_name, keep_gz=True)
        source_rows.append({
            "query_asm": asm_name,
            "role": "query",
            "query_mode": query_mode,
            "source_type": "assembly",
            "source_accession": row.get("assembly_accession") or ".",
            "source_url": fasta_src,
            "local_path": str(local_path),
            "species": row.get("species") or row.get("scientific_name") or ".",
        })
    else:
        direct_sources = direct_read_sources_from_row(row)
        accession = row.get("ena_accession") or row.get("sra_accession") or row.get("read_accession") or ""
        if not direct_sources and not accession:
            raise ValueError(
                f"Read-mode query row {asm_name} needs direct FASTQ/FASTA URLs or ena_accession/sra_accession/read_accession"
            )
        if accession:
            run_rows = fetch_ena_read_runs(accession)
            max_runs = int(row.get("max_runs") or public_max_runs)
            ena_sources, ena_meta = select_ena_read_sources(run_rows, query_mode, max_runs)
            direct_sources.extend(ena_sources)
            for meta in ena_meta:
                source_rows.append({
                    "query_asm": asm_name,
                    "role": "query",
                    "query_mode": query_mode,
                    "source_type": "ena_read_run",
                    "source_accession": meta["run_accession"],
                    "source_url": meta["source_url"],
                    "local_path": ".",
                    "species": meta["scientific_name"],
                })
        if not direct_sources:
            raise ValueError(f"No downloadable public read files resolved for query row {asm_name}")
        local_path = merge_sequence_sources(direct_sources, dest_dir / asm_name)
        for src in direct_sources:
            source_rows.append({
                "query_asm": asm_name,
                "role": "query",
                "query_mode": query_mode,
                "source_type": "public_read_file",
                "source_accession": accession or ".",
                "source_url": src,
                "local_path": str(local_path),
                "species": row.get("species") or row.get("scientific_name") or ".",
            })

    query_row = {
        "query_asm": asm_name,
        "query_mode": query_mode,
        "path": str(local_path),
        "scenario": row.get("scenario") or ".",
        "lifestyle": row.get("lifestyle") or ".",
        "architecture": row.get("architecture") or ".",
        "benchmark_ref_asm": row.get("benchmark_ref_asm") or benchmark.get("benchmark_ref_asm", "."),
        "benchmark_ref_fasta": row.get("benchmark_ref_fasta") or benchmark.get("benchmark_ref_fasta", "."),
        "phylum": row.get("phylum") or ".",
        "class": row.get("class") or ".",
        "order": row.get("order") or ".",
        "family": row.get("family") or ".",
        "genus": row.get("genus") or ".",
        "species": row.get("species") or row.get("scientific_name") or ".",
        "source": row.get("source") or default_source,
    }
    return query_row, str(local_path), source_rows


def prepare_from_custom_manifest(args: argparse.Namespace) -> int:
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    downloads_dir = out_dir / "downloads"
    refs_dir = downloads_dir / "refs"
    queries_dir = downloads_dir / "queries"
    rows = parse_custom_url_manifest(args.custom_url_manifest.resolve())

    ref_manifest_rows: list[dict[str, str]] = []
    query_rows: list[dict[str, str]] = []
    ref_list_paths: list[str] = []
    query_list_paths: list[str] = []
    source_link_rows: list[dict[str, str]] = []

    for row in rows:
        role = (row.get("role") or "ref").lower()
        asm_name = row.get("asm_name") or normalize_name(row.get("scientific_name") or row.get("clade_name") or "asm")

        if role == "ref":
            fasta_src = row.get("fasta_url") or row.get("path") or row.get("url") or ""
            if not fasta_src:
                raise ValueError(f"Missing fasta_url/path/url for reference manifest row {asm_name}")
            dest_name = Path(urllib.parse.urlparse(fasta_src).path).name or f"{asm_name}.fa"
            local_fasta = materialize_entry(fasta_src, refs_dir / dest_name, keep_gz=True)
            ref_list_paths.append(str(local_fasta))
            ref_manifest_rows.append({
                "asm_name": asm_name,
                "phylum": row.get("phylum") or ".",
                "class": row.get("class") or ".",
                "order": row.get("order") or ".",
                "family": row.get("family") or ".",
                "genus": row.get("genus") or ".",
                "clade_name": row.get("clade_name") or row.get("species") or row.get("scientific_name") or asm_name,
                "clade_rank": row.get("clade_rank") or "species",
                "fasta_path": str(local_fasta),
            })
            source_link_rows.append({
                "query_asm": asm_name,
                "role": "ref",
                "query_mode": "assembly",
                "source_type": "assembly",
                "source_accession": row.get("assembly_accession") or ".",
                "source_url": fasta_src,
                "local_path": str(local_fasta),
                "species": row.get("species") or row.get("scientific_name") or ".",
            })
        else:
            query_row, query_path, query_sources = materialize_query_input(
                row,
                queries_dir,
                default_source="custom",
                public_max_runs=args.public_query_max_runs,
            )
            query_rows.append(query_row)
            query_list_paths.append(query_path)
            source_link_rows.extend(query_sources)

    if not ref_manifest_rows:
        raise ValueError("Custom manifest did not produce any reference genomes")
    if not query_rows:
        raise ValueError("Custom manifest did not produce any query inputs")

    hierarchy_manifest = out_dir / "hierarchy_manifest.tsv"
    write_tsv(
        hierarchy_manifest,
        ref_manifest_rows,
        ["asm_name", "phylum", "class", "order", "family", "genus", "clade_name", "clade_rank", "fasta_path"],
    )

    (out_dir / "ref_list.txt").write_text("\n".join(ref_list_paths) + "\n", encoding="utf-8")
    (out_dir / "query_list.txt").write_text("\n".join(query_list_paths) + "\n", encoding="utf-8")
    write_tsv(
        out_dir / "query_manifest.tsv",
        query_rows,
        [
            "query_asm", "query_mode", "path", "scenario", "lifestyle", "architecture",
            "benchmark_ref_asm", "benchmark_ref_fasta", "phylum", "class", "order",
            "family", "genus", "species", "source",
        ],
    )
    write_tsv(
        out_dir / "public_data_links.tsv",
        source_link_rows,
        ["query_asm", "role", "query_mode", "source_type", "source_accession", "source_url", "local_path", "species"],
    )
    write_public_resource_links(out_dir / "public_resource_links.tsv")
    print(f"prepared\trefs={len(ref_manifest_rows)}\tqueries={len(query_rows)}\tmode=custom")
    return 0


def prepare_from_ncbi(args: argparse.Namespace) -> int:
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_url = NCBI_ASSEMBLY_SUMMARY[args.source]
    summary_text = http_get_text(summary_url)
    all_rows = parse_assembly_summary(summary_text)

    selectors: list[dict[str, str]] = []
    if args.all_public_assemblies:
        selectors = []
    elif args.species:
        for species in args.species:
            selectors.append({
                "species": species,
                "scenario": args.default_scenario,
                "lifestyle": args.default_lifestyle,
                "architecture": args.default_architecture,
            })
    else:
        for panel_name in args.panels:
            selectors.extend(PANEL_PRESETS[panel_name])

    selected_rows: list[dict[str, str]] = []
    catalog_rows: list[dict[str, str]] = []
    if args.all_public_assemblies:
        matches = select_all_public_rows(
            all_rows,
            min_assembly_level=args.min_assembly_level,
            latest_only=args.latest_only,
            max_total=args.max_public_assemblies,
        )
        for row in matches:
            tagged = dict(row)
            tagged["_scenario"] = args.default_scenario
            tagged["_lifestyle"] = args.default_lifestyle
            tagged["_architecture"] = args.default_architecture
            tagged["_target_species"] = species_label_for_row(row)
            tagged["_species_group_key"] = species_group_key(row)
            selected_rows.append(tagged)
            catalog_rows.append({
                "panel_species": tagged["_target_species"],
                "assembly_accession": row.get("assembly_accession", ""),
                "organism_name": row.get("organism_name", ""),
                "assembly_level": row.get("assembly_level", ""),
                "refseq_category": row.get("refseq_category", ""),
                "version_status": row.get("version_status", ""),
                "seq_rel_date": row.get("seq_rel_date", ""),
                "ftp_path": row.get("ftp_path", ""),
            })
    else:
        for sel in selectors:
            matches = select_species_rows(all_rows, sel["species"], args.max_assemblies_per_species)
            for row in matches:
                tagged = dict(row)
                tagged["_scenario"] = sel.get("scenario", ".")
                tagged["_lifestyle"] = sel.get("lifestyle", ".")
                tagged["_architecture"] = sel.get("architecture", ".")
                tagged["_target_species"] = sel["species"]
                tagged["_species_group_key"] = species_group_key(row)
                selected_rows.append(tagged)
                catalog_rows.append({
                    "panel_species": sel["species"],
                    "assembly_accession": row.get("assembly_accession", ""),
                    "organism_name": row.get("organism_name", ""),
                    "assembly_level": row.get("assembly_level", ""),
                    "refseq_category": row.get("refseq_category", ""),
                    "version_status": row.get("version_status", ""),
                    "seq_rel_date": row.get("seq_rel_date", ""),
                    "ftp_path": row.get("ftp_path", ""),
                })

    if not selected_rows:
        raise ValueError("No NCBI fungal assemblies matched the requested panel/species selection")

    write_tsv(
        out_dir / "selected_catalog.tsv",
        catalog_rows,
        [
            "panel_species", "assembly_accession", "organism_name", "assembly_level",
            "refseq_category", "version_status", "seq_rel_date", "ftp_path",
        ],
    )

    if args.catalog_only:
        print(f"catalog_only\tassemblies={len(selected_rows)}\tsource={args.source}")
        return 0

    taxonomy_cache = fetch_taxonomy_lineages(
        sorted({row.get("taxid", "") for row in selected_rows if row.get("taxid")}),
        cache_path=out_dir / "taxonomy_cache.json",
    )

    cache_base = args.data_cache_dir.resolve() if args.data_cache_dir else (out_dir / "downloads")
    refs_dir = cache_base / "refs"
    queries_dir = cache_base / "queries"
    refs_dir.mkdir(parents=True, exist_ok=True)
    queries_dir.mkdir(parents=True, exist_ok=True)
    ref_manifest_rows: list[dict[str, str]] = []
    query_rows: list[dict[str, str]] = []
    ref_list_paths: list[str] = []
    query_list_paths: list[str] = []
    benchmark_map_rows: list[dict[str, str]] = []
    source_link_rows: list[dict[str, str]] = []
    species_benchmark_defaults: dict[str, dict[str, str]] = {}

    by_species: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in selected_rows:
        by_species[row["_species_group_key"]].append(row)
    for rows in by_species.values():
        rows.sort(key=row_quality_key, reverse=True)

    ref_downloads = 0
    query_downloads = 0

    for _, rows in sorted(by_species.items()):
        species = rows[0].get("_target_species", species_label_for_row(rows[0]))
        lineage = taxonomy_cache.get(rows[0].get("taxid", ""), {})
        query_candidates = rows[: max(0, args.querys_per_species)]
        ref_candidates = rows[args.querys_per_species:]
        if not ref_candidates and len(rows) > 1:
            query_candidates = rows[-1:]
            ref_candidates = rows[:-1]
        if not ref_candidates:
            ref_candidates = rows[:1]
            query_candidates = rows[1:2]

        benchmark_ref_row = ref_candidates[0]
        benchmark_ref_local: Path | None = None
        species_ref_local: Path | None = None
        species_ref_asm: str = "."
        for row in ref_candidates:
            if args.max_ref_downloads > 0 and ref_downloads >= args.max_ref_downloads:
                break
            asm_name = row.get("assembly_accession", "").replace(".", "_")
            downloaded_fasta: Path | None = None
            for url, filename in ncbi_download_targets(row, include_gff=args.download_gff):
                local = materialize_entry(url, refs_dir / filename, keep_gz=True)
                if filename.endswith("_genomic.fna.gz"):
                    downloaded_fasta = local
            if downloaded_fasta is None:
                continue
            ref_downloads += 1
            if species_ref_local is None:
                species_ref_local = downloaded_fasta
                species_ref_asm = asm_name
            if benchmark_ref_row is row:
                benchmark_ref_local = downloaded_fasta
            ref_manifest_rows.append({
                "asm_name": asm_name,
                "phylum": lineage.get("phylum", "."),
                "class": lineage.get("class", "."),
                "order": lineage.get("order", "."),
                "family": lineage.get("family", "."),
                "genus": lineage.get("genus", species.split()[0]),
                "clade_name": lineage.get("species", species),
                "clade_rank": "species",
                "fasta_path": str(downloaded_fasta),
            })
            ref_list_paths.append(str(downloaded_fasta))
            source_link_rows.append({
                "query_asm": asm_name,
                "role": "ref",
                "query_mode": "assembly",
                "source_type": "ncbi_assembly",
                "source_accession": row.get("assembly_accession", "."),
                "source_url": row.get("ftp_path", "."),
                "local_path": str(downloaded_fasta),
                "species": lineage.get("species", species),
            })

        species_name = lineage.get("species", species)
        if benchmark_ref_local is None:
            benchmark_ref_local = species_ref_local
        species_benchmark_defaults[species_name] = {
            "benchmark_ref_asm": (
                benchmark_ref_row.get("assembly_accession", "").replace(".", "_")
                if benchmark_ref_local is not None and benchmark_ref_row.get("assembly_accession")
                else species_ref_asm
            ),
            "benchmark_ref_fasta": str(benchmark_ref_local) if benchmark_ref_local else ".",
        }

        if benchmark_ref_local is None:
            continue

        for row in query_candidates:
            if args.max_query_downloads > 0 and query_downloads >= args.max_query_downloads:
                break
            asm_name = row.get("assembly_accession", "").replace(".", "_")
            query_fasta: Path | None = None
            for url, filename in ncbi_download_targets(row, include_gff=False):
                local = materialize_entry(url, queries_dir / filename, keep_gz=True)
                if filename.endswith("_genomic.fna.gz"):
                    query_fasta = local
            if query_fasta is None:
                continue
            query_downloads += 1
            query_rows.append({
                "query_asm": asm_name,
                "query_mode": "assembly",
                "path": str(query_fasta),
                "scenario": row.get("_scenario", "."),
                "lifestyle": row.get("_lifestyle", "."),
                "architecture": row.get("_architecture", "."),
                "benchmark_ref_asm": species_benchmark_defaults[species_name]["benchmark_ref_asm"],
                "benchmark_ref_fasta": str(benchmark_ref_local) if benchmark_ref_local else ".",
                "phylum": lineage.get("phylum", "."),
                "class": lineage.get("class", "."),
                "order": lineage.get("order", "."),
                "family": lineage.get("family", "."),
                "genus": lineage.get("genus", species.split()[0]),
                "species": lineage.get("species", species),
                "source": args.source,
            })
            benchmark_map_rows.append({
                "query_asm": asm_name,
                "benchmark_ref_asm": species_benchmark_defaults[species_name]["benchmark_ref_asm"],
                "benchmark_ref_fasta": str(benchmark_ref_local) if benchmark_ref_local else ".",
                "species": lineage.get("species", species),
            })
            query_list_paths.append(str(query_fasta))
            source_link_rows.append({
                "query_asm": asm_name,
                "role": "query",
                "query_mode": "assembly",
                "source_type": "ncbi_assembly",
                "source_accession": row.get("assembly_accession", "."),
                "source_url": row.get("ftp_path", "."),
                "local_path": str(query_fasta),
                "species": lineage.get("species", species),
            })

    if args.public_query_manifest:
        extra_query_rows = parse_custom_url_manifest(args.public_query_manifest.resolve())
        for row in extra_query_rows:
            species = row.get("species") or row.get("scientific_name") or "."
            default_benchmark = species_benchmark_defaults.get(species, {})
            query_row, query_path, query_sources = materialize_query_input(
                row,
                queries_dir,
                default_source="public_query_manifest",
                default_benchmark=default_benchmark,
                public_max_runs=args.public_query_max_runs,
            )
            query_rows.append(query_row)
            query_list_paths.append(query_path)
            source_link_rows.extend(query_sources)
            benchmark_map_rows.append({
                "query_asm": query_row["query_asm"],
                "benchmark_ref_asm": query_row["benchmark_ref_asm"],
                "benchmark_ref_fasta": query_row["benchmark_ref_fasta"],
                "species": query_row["species"],
            })

    # ------------------------------------------------------------------
    # Reads-mode queries by ENA species lookup.
    #
    # When `--query-mode mixed` (default) or a reads variant is requested,
    # resolve public ENA read runs for each panel species and download up to
    # --read-accessions-per-species runs per (species, mode). This is the
    # mechanism that makes benchmark_short-reads/ and benchmark_long-reads/
    # non-empty for NCBI panels — without it the panels only produce
    # assembly-mode queries, which is why those folders were empty.
    # ------------------------------------------------------------------
    requested_read_modes: list[str] = []
    if args.query_mode == "mixed":
        requested_read_modes = ["short-reads", "long-reads"]
    elif args.query_mode in {"short-reads", "long-reads"}:
        requested_read_modes = [args.query_mode]

    if requested_read_modes and args.read_accessions_per_species > 0:
        for sel in selectors:
            species = sel.get("species", "")
            if not species:
                continue
            default_benchmark = species_benchmark_defaults.get(species, {})
            if not default_benchmark.get("benchmark_ref_fasta"):
                # No reference was downloaded for this species — reads-mode
                # benchmarks require one, so skip.
                continue
            ena_runs = fetch_ena_read_runs_by_species(species, max_rows=args.ena_max_rows_per_species)
            if not ena_runs:
                continue
            for read_mode in requested_read_modes:
                urls, meta_rows = select_ena_read_sources(
                    ena_runs, read_mode, args.read_accessions_per_species
                )
                if not urls:
                    continue
                # Bundle the selected FASTQs for the first picked run into a
                # single local file so query_list.txt stays one-path-per-line
                # (MycoSV reads that format).
                asm_name = normalize_name(f"{species}_{read_mode}_{meta_rows[0]['run_accession']}")
                try:
                    local_path = merge_sequence_sources(urls, queries_dir / asm_name)
                except Exception as exc:
                    sys.stderr.write(
                        f"[warn] {species}: ENA {read_mode} download failed: {exc}\n"
                    )
                    continue
                query_rows.append({
                    "query_asm": asm_name,
                    "query_mode": read_mode,
                    "path": str(local_path),
                    "scenario": sel.get("scenario", "."),
                    "lifestyle": sel.get("lifestyle", "."),
                    "architecture": sel.get("architecture", "."),
                    "benchmark_ref_asm": default_benchmark.get("benchmark_ref_asm", "."),
                    "benchmark_ref_fasta": default_benchmark.get("benchmark_ref_fasta", "."),
                    "phylum": ".", "class": ".", "order": ".", "family": ".",
                    "genus": species.split()[0] if species else ".",
                    "species": species,
                    "source": f"ena_{read_mode.replace('-', '_')}",
                    "instrument_platform": meta_rows[0].get("instrument_platform", "."),
                    "library_layout": meta_rows[0].get("library_layout", "."),
                    "run_accession": meta_rows[0].get("run_accession", "."),
                })
                query_list_paths.append(str(local_path))
                benchmark_map_rows.append({
                    "query_asm": asm_name,
                    "benchmark_ref_asm": default_benchmark.get("benchmark_ref_asm", "."),
                    "benchmark_ref_fasta": default_benchmark.get("benchmark_ref_fasta", "."),
                    "species": species,
                })
                for meta in meta_rows:
                    source_link_rows.append({
                        "query_asm": asm_name,
                        "role": "query",
                        "query_mode": read_mode,
                        "source_type": "ena_read_run",
                        "source_accession": meta.get("run_accession", "."),
                        "source_url": meta.get("source_url", "."),
                        "local_path": str(local_path),
                        "species": meta.get("scientific_name", species),
                    })

    if not query_rows:
        if not args.allow_no_queries:
            raise ValueError(
                "No query inputs were produced. Increase --max-assemblies-per-species, provide species with multiple assemblies, use --public-query-manifest, or pass --allow-no-queries for index-only preparation."
            )

    hierarchy_manifest = out_dir / "hierarchy_manifest.tsv"
    write_tsv(
        hierarchy_manifest,
        ref_manifest_rows,
        ["asm_name", "phylum", "class", "order", "family", "genus", "clade_name", "clade_rank", "fasta_path"],
    )
    (out_dir / "ref_list.txt").write_text("\n".join(ref_list_paths) + "\n", encoding="utf-8")
    query_list_text = ("\n".join(query_list_paths) + "\n") if query_list_paths else ""
    (out_dir / "query_list.txt").write_text(query_list_text, encoding="utf-8")
    write_tsv(
        out_dir / "query_manifest.tsv",
        query_rows,
        [
            "query_asm", "query_mode", "path", "scenario", "lifestyle", "architecture",
            "benchmark_ref_asm", "benchmark_ref_fasta", "phylum", "class", "order",
            "family", "genus", "species", "source",
        ],
    )
    write_tsv(out_dir / "benchmark_reference_map.tsv", benchmark_map_rows, ["query_asm", "benchmark_ref_asm", "benchmark_ref_fasta", "species"])
    write_tsv(
        out_dir / "public_data_links.tsv",
        source_link_rows,
        ["query_asm", "role", "query_mode", "source_type", "source_accession", "source_url", "local_path", "species"],
    )
    write_public_resource_links(out_dir / "public_resource_links.tsv")
    with (out_dir / "prepare_summary.json").open("w", encoding="utf-8") as fh:
        json.dump(
            {
                "source": args.source,
                "selected_rows": len(selected_rows),
                "ref_count": len(ref_manifest_rows),
                "query_count": len(query_rows),
                "panels": args.panels,
                "species_overrides": args.species,
                "all_public_assemblies": bool(args.all_public_assemblies),
                "max_public_assemblies": args.max_public_assemblies,
                "min_assembly_level": args.min_assembly_level,
                "latest_only": bool(args.latest_only),
                "max_ref_downloads": args.max_ref_downloads,
                "max_query_downloads": args.max_query_downloads,
                "allow_no_queries": bool(args.allow_no_queries),
                "public_query_manifest": str(args.public_query_manifest) if args.public_query_manifest else "",
            },
            fh,
            indent=2,
            sort_keys=True,
        )
    print(f"prepared\trefs={len(ref_manifest_rows)}\tqueries={len(query_rows)}\tsource={args.source}")
    return 0


def load_query_manifest(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows.append({k: (v or "").strip() for k, v in row.items()})
    return rows


def parse_info_field(field: str) -> dict[str, str]:
    info: dict[str, str] = {}
    for item in field.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
        else:
            info[item] = "1"
    return info


def qasm_matches_observed(expected: str, observed: str) -> bool:
    expected_norm = normalize_name(expected)
    observed_norm = normalize_name(observed)
    if expected_norm == observed_norm:
        return True
    return observed_norm.startswith(expected_norm + "_") or expected_norm.startswith(observed_norm + "_")


def load_mycosv_query_calls(vcf_path: Path, query_asm: str) -> list[NormalizedCall]:
    rows: list[NormalizedCall] = []
    if not vcf_path.exists():
        return rows
    with open_text_auto(vcf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info = parse_info_field(fields[7])
            if not qasm_matches_observed(query_asm, info.get("QASM", "")):
                continue
            svtype = TYPE_CANON.get(info.get("SVTYPE", fields[4].strip("<>")))
            if not svtype:
                continue
            pos = int(fields[1])
            end = int(info.get("END", pos))
            svlen = int(info.get("SVLEN", end - pos + 1))
            rows.append(NormalizedCall(
                query_asm=query_asm,
                query_contig=fields[0],
                pos=pos,
                end=end,
                svtype=svtype,
                svlen=svlen,
                source="mycosv",
                coord_space="query",
                annotation=info.get("ANNOT", "."),
                element_class=info.get("EC", "NONE"),
                ref_asm=info.get("CLADE", "."),
                ref_contig=info.get("REFCONTIG", "."),
            ))
    return rows


def load_mycosv_reference_calls(vcf_path: Path, query_asm: str) -> list[NormalizedCall]:
    rows: list[NormalizedCall] = []
    if not vcf_path.exists():
        return rows
    with open_text_auto(vcf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info = parse_info_field(fields[7])
            if not qasm_matches_observed(query_asm, info.get("QASM", "")):
                continue
            ref_contig = info.get("REFCONTIG", ".")
            ref_pos_raw = info.get("REFPOS", "")
            if ref_contig in {"", "."} or not ref_pos_raw:
                continue
            svtype = TYPE_CANON.get(info.get("SVTYPE", fields[4].strip("<>")))
            if not svtype:
                continue
            pos = int(ref_pos_raw)
            end = int(info.get("REFEND", pos))
            svlen = int(info.get("SVLEN", max(1, end - pos + 1)))
            rows.append(NormalizedCall(
                query_asm=query_asm,
                query_contig=fields[0],
                pos=pos,
                end=end,
                svtype=svtype,
                svlen=svlen,
                source="mycosv",
                coord_space="reference",
                annotation=info.get("ANNOT", "."),
                element_class=info.get("EC", "NONE"),
                ref_asm=info.get("CLADE", "."),
                ref_contig=ref_contig,
            ))
    return rows


def load_normalized_calls_tsv(path: Path, label: str) -> list[NormalizedCall]:
    rows: list[NormalizedCall] = []
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            svtype = TYPE_CANON.get((row.get("svtype") or row.get("type") or "").upper())
            if not svtype:
                continue
            pos = int(float(row.get("pos") or 0))
            end = int(float(row.get("end") or pos))
            svlen_raw = row.get("svlen") or row.get("length") or ""
            svlen = int(float(svlen_raw)) if svlen_raw else (end - pos + 1)
            coord_space = (row.get("coord_space") or row.get("coordinate_space") or "query").strip().lower()
            if coord_space not in {"query", "reference"}:
                coord_space = "query"
            query_contig = row.get("query_contig") or row.get("q_contig") or "."
            ref_contig = row.get("ref_contig") or row.get("r_contig") or "."
            shared_contig = row.get("chrom") or row.get("contig") or row.get("seqid") or "."
            if coord_space == "query" and query_contig == ".":
                query_contig = shared_contig
            if coord_space == "reference" and ref_contig == ".":
                ref_contig = shared_contig
            rows.append(NormalizedCall(
                query_asm=row.get("query_asm") or row.get("sample") or ".",
                query_contig=query_contig,
                pos=pos,
                end=end,
                svtype=svtype,
                svlen=svlen,
                source=label,
                coord_space=coord_space,
                annotation=row.get("annotation", "."),
                element_class=row.get("element_class", "NONE"),
                ref_asm=row.get("ref_asm", "."),
                ref_contig=ref_contig,
            ))
    return rows


def load_syri_query_calls(path: Path, query_asm: str) -> list[NormalizedCall]:
    rows: list[NormalizedCall] = []
    if not path.exists():
        return rows
    with open_text_auto(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for parts in reader:
            if not parts or len(parts) < 11:
                continue
            ann = parts[10].upper()
            svtype = TYPE_CANON.get(ann)
            if not svtype:
                continue
            qry_contig = parts[5]
            q_start = int(parts[6])
            q_end = int(parts[7])
            if svtype == "TRA":
                svlen = max(1, q_end - q_start + 1)
            else:
                svlen = max(1, q_end - q_start + 1)
            rows.append(NormalizedCall(
                query_asm=query_asm,
                query_contig=qry_contig,
                pos=min(q_start, q_end),
                end=max(q_start, q_end),
                svtype=svtype,
                svlen=svlen,
                source="syri",
                coord_space="query",
                annotation=ann,
                ref_contig=parts[0],
            ))
    return rows


def load_reference_vcf_calls(path: Path, label: str, query_asm: str) -> list[NormalizedCall]:
    rows: list[NormalizedCall] = []
    if not path.exists():
        return rows
    with open_text_auto(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            info = parse_info_field(fields[7])
            svtype = TYPE_CANON.get((info.get("SVTYPE") or fields[4].strip("<>")).upper())
            if not svtype:
                continue
            pos = int(fields[1])
            # BND records don't carry END; fall back to pos.
            end_raw = info.get("END", "")
            try:
                end = int(end_raw) if end_raw else pos
            except ValueError:
                end = pos
            svlen_raw = info.get("SVLEN", "")
            try:
                svlen = abs(int(svlen_raw)) if svlen_raw else max(1, end - pos + 1)
            except ValueError:
                svlen = max(1, end - pos + 1)
            if svlen <= 0:
                svlen = max(1, end - pos + 1)
            rows.append(NormalizedCall(
                query_asm=query_asm,
                query_contig=".",
                pos=pos,
                end=end,
                svtype=svtype,
                svlen=svlen,
                source=label,
                coord_space="reference",
                annotation=info.get("ANNOT", "."),
                element_class=info.get("EC", "NONE"),
                ref_asm=info.get("REFASM", "."),
                ref_contig=fields[0],
            ))
    return rows


def parse_minigraph_call_tail(field: str) -> dict[str, str] | None:
    parts = field.rsplit(":", 5)
    if len(parts) != 6:
        return None
    return {
        "path": parts[0],
        "path_len": parts[1],
        "strand": parts[2],
        "query_contig": parts[3],
        "query_start": parts[4],
        "query_end": parts[5],
    }


def load_minigraph_bubble_calls(bubble_bed: Path, sample_bed: Path, query_asm: str) -> list[NormalizedCall]:
    rows: list[NormalizedCall] = []
    if not bubble_bed.exists() or not sample_bed.exists():
        return rows
    bubble_lines = [line.rstrip("\n") for line in bubble_bed.read_text(encoding="utf-8").splitlines() if line.strip()]
    sample_lines = [line.rstrip("\n") for line in sample_bed.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(bubble_lines) != len(sample_lines):
        return rows
    for bubble_line, sample_line in zip(bubble_lines, sample_lines):
        b = bubble_line.split("\t")
        s = sample_line.split("\t")
        if len(b) < 8 or not s:
            continue
        tail = parse_minigraph_call_tail(s[-1])
        if tail is None:
            continue
        chrom = b[0]
        ref_start0 = int(b[1])
        ref_end = int(b[2])
        inv_flag = int(b[5])
        ref_len = max(0, ref_end - ref_start0)
        q_start = int(float(tail["query_start"]))
        q_end = int(float(tail["query_end"]))
        sample_len = int(float(tail["path_len"])) if tail["path_len"] else abs(q_end - q_start)
        if inv_flag == 1:
            svtype = "INV"
            svlen = max(1, ref_len)
        elif sample_len > ref_len:
            svtype = "INS"
            svlen = max(1, sample_len - ref_len)
        elif sample_len < ref_len:
            svtype = "DEL"
            svlen = max(1, ref_len - sample_len)
        else:
            continue
        rows.append(NormalizedCall(
            query_asm=query_asm,
            query_contig=tail["query_contig"],
            pos=ref_start0 + 1,
            end=max(ref_start0 + 1, ref_end),
            svtype=svtype,
            svlen=svlen,
            source="minigraph",
            coord_space="reference",
            annotation="MINIGRAPH_BUBBLE",
            ref_contig=chrom,
        ))
    return rows


def canonical_group(svtype: str) -> str:
    if svtype in {"TRANS", "INVTR", "TRA"}:
        return "TRA"
    if svtype in {"DUP", "INVDP"}:
        return "DUP"
    return svtype


def effective_contig(call: NormalizedCall) -> str:
    return call.ref_contig if call.coord_space == "reference" else call.query_contig


def calls_compatible(truth: NormalizedCall, pred: NormalizedCall) -> bool:
    if truth.coord_space != pred.coord_space:
        return False
    if truth.query_asm not in {"", ".", pred.query_asm} and pred.query_asm not in {"", "."}:
        return False
    if effective_contig(truth) != effective_contig(pred):
        return False
    if canonical_group(truth.svtype) != canonical_group(pred.svtype):
        return False
    tol_bp = DEFAULT_TOL_BP.get(canonical_group(truth.svtype), 500)
    tol_frac = DEFAULT_TOL_LEN_FRAC.get(canonical_group(truth.svtype), 0.30)
    if abs(truth.pos - pred.pos) > tol_bp:
        return False
    if truth.svtype not in {"TRA", "OFF_REF", "INS"}:
        denom = max(abs(truth.svlen), 1)
        if abs(abs(truth.svlen) - abs(pred.svlen)) / denom > tol_frac:
            return False
    return True


def call_distance(truth: NormalizedCall, pred: NormalizedCall) -> int:
    return abs(truth.pos - pred.pos) + abs(abs(truth.svlen) - abs(pred.svlen))


def match_calls(truth_calls: list[NormalizedCall], pred_calls: list[NormalizedCall]) -> tuple[set[int], list[int]]:
    used: set[int] = set()
    missed_truth: list[int] = []
    for truth_idx, truth in enumerate(truth_calls):
        best_idx: int | None = None
        best_dist = 10**18
        for pred_idx, pred in enumerate(pred_calls):
            if pred_idx in used or not calls_compatible(truth, pred):
                continue
            dist = call_distance(truth, pred)
            if dist < best_dist:
                best_idx = pred_idx
                best_dist = dist
        if best_idx is None:
            missed_truth.append(truth_idx)
        else:
            used.add(best_idx)
    return used, missed_truth


def score_callsets(truth_calls: list[NormalizedCall], pred_calls: list[NormalizedCall]) -> dict[str, Any]:
    used, missed_truth = match_calls(truth_calls, pred_calls)
    tp = len(used)
    fn = len(missed_truth)
    fp = len(pred_calls) - tp
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = (2.0 * precision * recall / (precision + recall)) if precision + recall else 0.0
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "precision_ci95": wilson_ci(tp, tp + fp),
        "recall_ci95": wilson_ci(tp, tp + fn),
    }


def parse_other_spec(spec: str) -> tuple[str, Path]:
    if "=" not in spec:
        raise ValueError(f"Expected label=path, got {spec!r}")
    label, path = spec.split("=", 1)
    label = label.strip()
    path = path.strip()
    if not label or not path:
        raise ValueError(f"Expected label=path, got {spec!r}")
    return label, Path(path).resolve()


def compile_binary_if_needed(binary_path: Path, force: bool = False) -> None:
    if binary_path.exists() and not force:
        return
    run(["g++", "-O2", "-std=c++17", "-pthread", str(ROOT / "main.cpp"), "-o", str(binary_path)], cwd=ROOT)


def locate_query_path(query_row: dict[str, str]) -> Path:
    return Path(query_row["path"]).resolve()


def estimate_prepared_genome_size_hint(prepared_dir: Path) -> int:
    query_manifest = prepared_dir / "query_manifest.tsv"
    candidate_paths: list[Path] = []
    if query_manifest.exists():
        for row in load_query_manifest(query_manifest):
            ref_fasta = (row.get("benchmark_ref_fasta") or "").strip()
            if ref_fasta and ref_fasta not in {".", ""}:
                candidate_paths.append(Path(ref_fasta))
    ref_list = prepared_dir / "ref_list.txt"
    if ref_list.exists():
        for line in ref_list.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line:
                candidate_paths.append(Path(line))

    seen: set[Path] = set()
    for path in candidate_paths:
        path = path.resolve()
        if path in seen or not path.exists():
            continue
        seen.add(path)
        total = 0
        with open_text_auto(path) as fh:
            for line in fh:
                if not line or line.startswith(">"):
                    continue
                total += len(line.strip())
        if total > 0:
            return total
    return 0


def run_mycosv(
    prepared_dir: Path,
    out_dir: Path,
    binary_path: Path,
    mode: str,
    extra_args: list[str],
    *,
    query_list_override: Path | None = None,
) -> dict[str, str]:
    mycosv_dir = out_dir / "mycosv"
    mycosv_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = mycosv_dir / "calls"
    caller_args = list(extra_args)
    if mode == "short-reads" and "--max-reads" not in caller_args:
        caller_args.extend(["--max-reads", "150000"])
    if mode == "long-reads" and "--max-reads" not in caller_args:
        caller_args.extend(["--max-reads", "100"])
    if mode != "assembly" and "--genome-size-hint" not in caller_args:
        genome_size_hint = estimate_prepared_genome_size_hint(prepared_dir)
        if genome_size_hint > 0:
            caller_args.extend(["--genome-size-hint", str(genome_size_hint)])
    query_list_path = (query_list_override.resolve() if query_list_override
                       else (prepared_dir / "query_list.txt").resolve())
    cmd = [
        str(binary_path.resolve()),
        "--ref-list", str((prepared_dir / "ref_list.txt").resolve()),
        "--query-list", str(query_list_path),
        "--out-prefix", str(out_prefix.resolve()),
        "--query-mode", mode,
        *caller_args,
    ]
    run_mycosv_command(cmd, cwd=ROOT)
    return {
        "vcf": str(out_prefix.with_suffix(".vcf")),
        "hits": str(out_prefix.with_suffix(".hits.tsv")),
        "gfa": str(out_prefix.with_suffix(".gfa")),
    }


def write_prefixed_fasta_records(src: Path, prefix: str, out_fh: io.TextIOBase) -> None:
    with src.open(encoding="utf-8") as fh:
        header = ""
        seq_lines: list[str] = []
        for line in fh:
            if line.startswith(">"):
                if header:
                    out_fh.write(f">{prefix}#1#{header}\n")
                    out_fh.write("".join(seq_lines))
                header = line[1:].strip().split()[0]
                seq_lines = []
            else:
                seq_lines.append(line if line.endswith("\n") else line + "\n")
        if header:
            out_fh.write(f">{prefix}#1#{header}\n")
            out_fh.write("".join(seq_lines))


def build_pair_fasta(query_row: dict[str, str], out_path: Path) -> Path | None:
    ref_fasta = query_row.get("benchmark_ref_fasta", ".")
    if ref_fasta in {"", "."}:
        return None
    ref_fa = Path(ref_fasta).resolve()
    query_fa = locate_query_path(query_row)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as out_fh:
        write_prefixed_fasta_records(ref_fa, "ref", out_fh)
        write_prefixed_fasta_records(query_fa, normalize_name(query_row["query_asm"]), out_fh)
    return out_path


def run_syri_for_query(query_row: dict[str, str], out_dir: Path, threads: int) -> dict[str, str] | None:
    if not tool_path("minimap2") or not tool_path("syri"):
        return None
    query_asm = query_row["query_asm"]
    ref_fasta = query_row.get("benchmark_ref_fasta", ".")
    if ref_fasta in {"", "."}:
        return None
    query_fa = locate_query_path(query_row)
    ref_fa = Path(ref_fasta).resolve()
    work_dir = out_dir / "comparators" / "syri" / query_asm
    work_dir.mkdir(parents=True, exist_ok=True)
    sam_path = work_dir / "query_vs_ref.sam"
    prefix = str(work_dir / "syri_")
    with sam_path.open("w", encoding="utf-8") as sam_out:
        subprocess.run(
            ["minimap2", "-ax", "asm5", "--eqx", "-t", str(threads), str(ref_fa), str(query_fa)],
            stdout=sam_out,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
    run(["syri", "-c", str(sam_path), "-r", str(ref_fa), "-q", str(query_fa), "-k", "-F", "S", "--prefix", prefix], cwd=ROOT)
    syri_tsv = work_dir / "syri_syri.out"
    if not syri_tsv.exists():
        candidates = sorted(work_dir.glob("*syri.out"))
        if not candidates:
            return None
        syri_tsv = candidates[0]
    return {"label": "syri", "normalized_tsv": str(syri_tsv)}


def run_minigraph_for_query(query_row: dict[str, str], out_dir: Path, threads: int, extra_args: list[str]) -> dict[str, str] | None:
    if not tool_path("minigraph") or not tool_path("gfatools"):
        return None
    ref_fasta = query_row.get("benchmark_ref_fasta", ".")
    if ref_fasta in {"", "."}:
        return None
    query_asm = query_row["query_asm"]
    query_fa = locate_query_path(query_row)
    ref_fa = Path(ref_fasta).resolve()
    work_dir = out_dir / "comparators" / "minigraph" / query_asm
    work_dir.mkdir(parents=True, exist_ok=True)
    graph_gfa = work_dir / "graph.gfa"
    bubble_bed = work_dir / "bubbles.bed"
    sample_bed = work_dir / "sample.bed"
    with graph_gfa.open("w", encoding="utf-8") as out_fh:
        subprocess.run(
            ["minigraph", "-cxggs", "-c", "-t", str(threads), *extra_args, str(ref_fa), str(query_fa)],
            stdout=out_fh,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
    with bubble_bed.open("w", encoding="utf-8") as out_fh:
        subprocess.run(
            ["gfatools", "bubble", str(graph_gfa)],
            stdout=out_fh,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
    with sample_bed.open("w", encoding="utf-8") as out_fh:
        subprocess.run(
            ["minigraph", "-cxasm", "--call", "-t", str(threads), *extra_args, str(graph_gfa), str(query_fa)],
            stdout=out_fh,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
    return {
        "label": "minigraph",
        "bubble_bed": str(bubble_bed),
        "sample_bed": str(sample_bed),
        "graph_gfa": str(graph_gfa),
    }


def run_pggb_for_query(query_row: dict[str, str], out_dir: Path, threads: int, identity: str, segment_len: str, extra_args: list[str]) -> dict[str, str] | None:
    if not tool_path("pggb"):
        return None
    query_asm = query_row["query_asm"]
    work_dir = out_dir / "comparators" / "pggb" / query_asm
    work_dir.mkdir(parents=True, exist_ok=True)
    pair_fa = build_pair_fasta(query_row, work_dir / "pair.fa")
    if pair_fa is None:
        return None
    if tool_path("samtools"):
        try:
            run(["samtools", "faidx", str(pair_fa)], cwd=ROOT)
        except subprocess.CalledProcessError:
            pass
    cmd = [
        "pggb",
        "-i", str(pair_fa),
        "-o", str(work_dir),
        "-n", "2",
        "-t", str(threads),
        "-p", str(identity),
        "-s", str(segment_len),
        "-V", "ref:1000",
        *extra_args,
    ]
    run(cmd, cwd=ROOT)
    vcf_candidates = sorted(work_dir.glob("**/*.vcf")) + sorted(work_dir.glob("**/*.vcf.gz"))
    if not vcf_candidates:
        return None
    return {"label": "pggb", "vcf": str(vcf_candidates[0])}


# ============================================================================
# Read-based SV caller adapters
#
# SVIM, Sniffles, cuteSV     -> long reads  (minimap2 -ax map-ont|map-pb)
# Delly, Manta               -> short reads (minimap2 -ax sr)
#
# All of these consume a sorted-indexed BAM of query reads aligned to the
# reference and produce a reference-coordinate VCF that
# load_reference_vcf_calls() already knows how to parse — so the truth_set
# plumbing in benchmark_real_data() picks them up the same way pggb does.
# ============================================================================


def _minimap2_align_reads(
    query_row: dict[str, str],
    work_dir: Path,
    threads: int,
    *,
    preset: str,
) -> tuple[Path, Path] | None:
    """
    Align a read-mode query (FASTQ/FASTA) to its benchmark reference with
    minimap2 and produce a sorted+indexed BAM. Returns (bam_path, ref_path)
    or None if the prerequisites aren't met.

    We use a single input file (manifest's `path`). For paired short reads the
    project's manifests store either an interleaved FASTQ or R1 only; minimap2
    handles both with the `sr` preset.
    """
    if not tool_path("minimap2") or not tool_path("samtools"):
        return None
    ref_fasta = query_row.get("benchmark_ref_fasta", ".")
    if ref_fasta in {"", "."}:
        return None
    ref_fa = Path(ref_fasta).resolve()
    reads_path = locate_query_path(query_row)
    if not reads_path.exists():
        return None
    work_dir.mkdir(parents=True, exist_ok=True)

    sam_path = work_dir / "aln.sam"
    bam_sorted = work_dir / "aln.sorted.bam"

    # minimap2 -> SAM
    with sam_path.open("w", encoding="utf-8") as sam_out:
        subprocess.run(
            ["minimap2", "-ax", preset, "-t", str(threads), str(ref_fa), str(reads_path)],
            stdout=sam_out,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
    # samtools sort -> BAM
    run(
        ["samtools", "sort", "-@", str(threads), "-o", str(bam_sorted), str(sam_path)],
        cwd=ROOT,
    )
    # samtools index
    run(["samtools", "index", str(bam_sorted)], cwd=ROOT)
    # Reference needs a .fai for downstream callers (Delly/Manta especially).
    if not (ref_fa.parent / (ref_fa.name + ".fai")).exists():
        try:
            run(["samtools", "faidx", str(ref_fa)], cwd=ROOT)
        except subprocess.CalledProcessError:
            pass
    # Drop the giant intermediate SAM once the BAM exists.
    try:
        sam_path.unlink()
    except OSError:
        pass
    return bam_sorted, ref_fa


def _long_read_preset(query_row: dict[str, str]) -> str:
    """Pick a minimap2 long-read preset from the manifest platform hint."""
    platform = (query_row.get("instrument_platform") or "").strip().upper()
    if "PACBIO" in platform or "PB" in platform:
        return "map-pb"
    # Default to Nanopore: most public long-read fungal data is ONT.
    return "map-ont"


def run_svim_for_query(query_row: dict[str, str], out_dir: Path, threads: int) -> dict[str, str] | None:
    """SVIM: long-read SV caller. Produces variants.vcf in its output dir."""
    if not tool_path("svim"):
        return None
    query_asm = query_row["query_asm"]
    work_dir = out_dir / "comparators" / "svim" / query_asm
    aligned = _minimap2_align_reads(
        query_row, work_dir, threads, preset=_long_read_preset(query_row)
    )
    if aligned is None:
        return None
    bam_sorted, ref_fa = aligned
    svim_out = work_dir / "svim_out"
    svim_out.mkdir(parents=True, exist_ok=True)
    try:
        run(
            ["svim", "alignment", str(svim_out), str(bam_sorted), str(ref_fa)],
            cwd=ROOT,
        )
    except subprocess.CalledProcessError:
        return None
    vcf_path = svim_out / "variants.vcf"
    if not vcf_path.exists():
        candidates = sorted(svim_out.rglob("variants.vcf"))
        if not candidates:
            return None
        vcf_path = candidates[0]
    return {"label": "svim", "vcf": str(vcf_path)}


def run_sniffles_for_query(query_row: dict[str, str], out_dir: Path, threads: int) -> dict[str, str] | None:
    """Sniffles2: long-read SV caller. Emits a reference-coordinate VCF."""
    if not tool_path("sniffles"):
        return None
    query_asm = query_row["query_asm"]
    work_dir = out_dir / "comparators" / "sniffles" / query_asm
    aligned = _minimap2_align_reads(
        query_row, work_dir, threads, preset=_long_read_preset(query_row)
    )
    if aligned is None:
        return None
    bam_sorted, ref_fa = aligned
    vcf_path = work_dir / "sniffles.vcf"
    try:
        run(
            ["sniffles", "--input", str(bam_sorted),
             "--reference", str(ref_fa),
             "--vcf", str(vcf_path),
             "--threads", str(threads)],
            cwd=ROOT,
        )
    except subprocess.CalledProcessError:
        return None
    if not vcf_path.exists():
        return None
    return {"label": "sniffles", "vcf": str(vcf_path)}


def run_cutesv_for_query(query_row: dict[str, str], out_dir: Path, threads: int) -> dict[str, str] | None:
    """cuteSV: long-read SV caller. Needs a writable working dir alongside the BAM."""
    if not tool_path("cuteSV") and not tool_path("cutesv"):
        return None
    bin_name = "cuteSV" if tool_path("cuteSV") else "cutesv"
    query_asm = query_row["query_asm"]
    work_dir = out_dir / "comparators" / "cutesv" / query_asm
    aligned = _minimap2_align_reads(
        query_row, work_dir, threads, preset=_long_read_preset(query_row)
    )
    if aligned is None:
        return None
    bam_sorted, ref_fa = aligned
    vcf_path = work_dir / "cutesv.vcf"
    tmp_dir = work_dir / "cutesv_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    # Default parameters follow the cuteSV README recommendations for ONT;
    # they're conservative enough for fungal genomes and robust when the
    # platform is mis-reported.
    cmd = [
        bin_name,
        str(bam_sorted),
        str(ref_fa),
        str(vcf_path),
        str(tmp_dir),
        "--threads", str(threads),
        "--max_cluster_bias_INS", "100",
        "--diff_ratio_merging_INS", "0.3",
        "--max_cluster_bias_DEL", "100",
        "--diff_ratio_merging_DEL", "0.3",
        "--min_support", "3",
        "--sample", query_asm,
    ]
    try:
        run(cmd, cwd=ROOT)
    except subprocess.CalledProcessError:
        return None
    if not vcf_path.exists():
        return None
    return {"label": "cutesv", "vcf": str(vcf_path)}


def run_delly_for_query(query_row: dict[str, str], out_dir: Path, threads: int) -> dict[str, str] | None:
    """Delly: short-read SV caller. Emits BCF by default; we convert to VCF."""
    if not tool_path("delly") or not tool_path("bcftools"):
        return None
    query_asm = query_row["query_asm"]
    work_dir = out_dir / "comparators" / "delly" / query_asm
    aligned = _minimap2_align_reads(query_row, work_dir, threads, preset="sr")
    if aligned is None:
        return None
    bam_sorted, ref_fa = aligned
    bcf_path = work_dir / "delly.bcf"
    vcf_path = work_dir / "delly.vcf"
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", str(threads))
    try:
        subprocess.run(
            ["delly", "call", "-g", str(ref_fa), "-o", str(bcf_path), str(bam_sorted)],
            check=True, text=True, capture_output=True, env=env, cwd=str(ROOT),
        )
    except subprocess.CalledProcessError:
        return None
    if not bcf_path.exists():
        return None
    # BCF -> text VCF so the generic loader can parse it.
    with vcf_path.open("w", encoding="utf-8") as out_fh:
        subprocess.run(
            ["bcftools", "view", str(bcf_path)],
            stdout=out_fh, stderr=subprocess.PIPE, text=True, check=True,
        )
    return {"label": "delly", "vcf": str(vcf_path)}


def run_manta_for_query(query_row: dict[str, str], out_dir: Path, threads: int) -> dict[str, str] | None:
    """Manta: short-read SV caller. Runs configManta.py then runWorkflow.py."""
    configure = tool_path("configManta.py")
    if configure is None:
        return None
    query_asm = query_row["query_asm"]
    work_dir = out_dir / "comparators" / "manta" / query_asm
    aligned = _minimap2_align_reads(query_row, work_dir, threads, preset="sr")
    if aligned is None:
        return None
    bam_sorted, ref_fa = aligned
    run_dir = work_dir / "manta_run"
    # Manta refuses to overwrite an existing run dir — clean it out first.
    if run_dir.exists():
        shutil.rmtree(run_dir, ignore_errors=True)
    try:
        run(
            [configure,
             "--bam", str(bam_sorted),
             "--referenceFasta", str(ref_fa),
             "--runDir", str(run_dir)],
            cwd=ROOT,
        )
        run(
            [str(run_dir / "runWorkflow.py"), "-j", str(threads)],
            cwd=ROOT,
        )
    except subprocess.CalledProcessError:
        return None
    # Manta's diploid SV VCF (most inclusive set including BNDs) lives here.
    vcf_candidates = [
        run_dir / "results" / "variants" / "diploidSV.vcf.gz",
        run_dir / "results" / "variants" / "candidateSV.vcf.gz",
    ]
    for cand in vcf_candidates:
        if cand.exists():
            return {"label": "manta", "vcf": str(cand)}
    return None


# ============================================================================
# Fungi-oriented assembly-mode comparators
#
# Minigraph-Cactus (cactus-pangenome)  -> pangenome graph -> reference VCF
# SVIM-asm (haploid mode)              -> assembly-to-assembly SV caller
# AnchorWave (anchorwave proali/genoAli) -> WGA-based SV detection
#
# All three produce reference-coordinate VCFs that the generic
# load_reference_vcf_calls() loader consumes, so they plug into the same
# truth_sets / exact_benchmark_summary.tsv machinery as pggb/minigraph.
# ============================================================================


def run_cactus_for_query(
    query_row: dict[str, str],
    out_dir: Path,
    threads: int,
    extra_args: list[str],
) -> dict[str, str] | None:
    """
    Minigraph-Cactus pangenome pipeline on a pairwise (ref, query) seqfile.
    Writes <work_dir>/<outName>.vcf.gz, which is the vcfbub-filtered,
    reference-coordinate VCF documented in cactus/doc/pangenome.md.

    We use the haploid ".0" naming convention so genome paths are stored in
    the vg/Giraffe indexes correctly.
    """
    # The entry point is 'cactus-pangenome'; older wheels ship 'cactus' too but
    # the pangenome CLI is the stable name.
    if not tool_path("cactus-pangenome"):
        return None
    ref_fasta = query_row.get("benchmark_ref_fasta", ".")
    if ref_fasta in {"", "."}:
        return None
    ref_fa = Path(ref_fasta).resolve()
    query_fa = locate_query_path(query_row)
    if not query_fa.exists():
        return None

    query_asm = query_row["query_asm"]
    work_dir = out_dir / "comparators" / "cactus" / query_asm
    work_dir.mkdir(parents=True, exist_ok=True)

    # Seqfile: 2 columns (name, path). Reference has no haplotype suffix per
    # the MC convention; the query uses ".0" to mark it as a haploid sample.
    ref_name = "ref"
    query_name = f"{normalize_name(query_asm)}.0"
    seqfile = work_dir / "seqfile.tsv"
    with seqfile.open("w", encoding="utf-8") as fh:
        fh.write(f"{ref_name}\t{ref_fa}\n")
        fh.write(f"{query_name}\t{query_fa}\n")

    # Toil jobstore must NOT exist at launch; Toil recreates it per run.
    job_store = work_dir / "jobstore"
    if job_store.exists():
        shutil.rmtree(job_store, ignore_errors=True)
    out_name = "pangenome"

    cmd = [
        "cactus-pangenome",
        str(job_store),
        str(seqfile),
        "--outDir", str(work_dir),
        "--outName", out_name,
        "--reference", ref_name,
        "--vcf",                      # emit reference-coordinate VCF
        "--mapCores", str(threads),
        "--maxCores", str(threads),
        *extra_args,
    ]
    try:
        run(cmd, cwd=ROOT)
    except subprocess.CalledProcessError:
        return None

    # Cactus writes <outName>.vcf.gz in the outDir.
    vcf_candidates = [
        work_dir / f"{out_name}.vcf.gz",
        work_dir / f"{out_name}.raw.vcf.gz",
    ]
    vcf_candidates += sorted(work_dir.rglob(f"{out_name}*.vcf.gz"))
    for cand in vcf_candidates:
        if cand.exists():
            return {"label": "cactus", "vcf": str(cand)}
    return None


def run_svim_asm_for_query(query_row: dict[str, str], out_dir: Path, threads: int) -> dict[str, str] | None:
    """
    SVIM-asm haploid mode: minimap2 -ax asm5 query -> ref, sort+index, then
    `svim-asm haploid <out> <bam> <ref>`. Produces variants.vcf with all 5
    canonical SV types (INS/DEL/DUP tandem & interspersed/INV/TRA via BND).
    """
    if not tool_path("svim-asm"):
        return None
    if not tool_path("minimap2") or not tool_path("samtools"):
        return None
    ref_fasta = query_row.get("benchmark_ref_fasta", ".")
    if ref_fasta in {"", "."}:
        return None
    ref_fa = Path(ref_fasta).resolve()
    query_fa = locate_query_path(query_row)
    if not query_fa.exists():
        return None

    query_asm = query_row["query_asm"]
    work_dir = out_dir / "comparators" / "svim_asm" / query_asm
    work_dir.mkdir(parents=True, exist_ok=True)

    sam_path = work_dir / "query_vs_ref.sam"
    bam_sorted = work_dir / "query_vs_ref.sorted.bam"
    with sam_path.open("w", encoding="utf-8") as sam_out:
        subprocess.run(
            ["minimap2", "-ax", "asm5", "-t", str(threads), str(ref_fa), str(query_fa)],
            stdout=sam_out, stderr=subprocess.PIPE, text=True, check=True,
        )
    try:
        run(
            ["samtools", "sort", "-@", str(threads), "-o", str(bam_sorted), str(sam_path)],
            cwd=ROOT,
        )
        run(["samtools", "index", str(bam_sorted)], cwd=ROOT)
    except subprocess.CalledProcessError:
        return None
    try:
        sam_path.unlink()
    except OSError:
        pass

    svim_out = work_dir / "svim_asm_out"
    svim_out.mkdir(parents=True, exist_ok=True)
    try:
        run(
            ["svim-asm", "haploid", str(svim_out), str(bam_sorted), str(ref_fa)],
            cwd=ROOT,
        )
    except subprocess.CalledProcessError:
        return None

    vcf_path = svim_out / "variants.vcf"
    if not vcf_path.exists():
        candidates = sorted(svim_out.rglob("variants.vcf"))
        if not candidates:
            return None
        vcf_path = candidates[0]
    return {"label": "svim_asm", "vcf": str(vcf_path)}


def run_anchorwave_for_query(query_row: dict[str, str], out_dir: Path, threads: int) -> dict[str, str] | None:
    """
    AnchorWave genoAli: whole-genome alignment for organisms with large SVs
    and high repetitive content. Widely used in fungal/plant pangenome work.

    Pipeline: minimap2 CDS-independent path ->  anchorwave genoAli -> MAF ->
    anchorwave-provided util translates MAF to VCF (maf-convert + custom
    post-processing). We keep the pipeline conservative: if the necessary
    ancillary tools aren't present we bail out rather than emit a partial
    result.
    """
    anchorwave = tool_path("anchorwave")
    if not anchorwave:
        return None
    if not tool_path("minimap2") or not tool_path("samtools"):
        return None
    # Optional helper: maf2vcf comes from vcf-kit or the anchorwave cookbook.
    # We prefer `paftools.js` (shipped with minimap2) since it converts PAF
    # -> VCF with structural variants and is universally available.
    paftools = tool_path("paftools.js")
    if paftools is None:
        return None

    ref_fasta = query_row.get("benchmark_ref_fasta", ".")
    if ref_fasta in {"", "."}:
        return None
    ref_fa = Path(ref_fasta).resolve()
    query_fa = locate_query_path(query_row)
    if not query_fa.exists():
        return None

    query_asm = query_row["query_asm"]
    work_dir = out_dir / "comparators" / "anchorwave" / query_asm
    work_dir.mkdir(parents=True, exist_ok=True)

    sam_path = work_dir / "q2r.sam"
    paf_path = work_dir / "q2r.paf"
    sorted_paf = work_dir / "q2r.srt.paf"
    vcf_path = work_dir / "anchorwave.vcf"

    # Step 1: minimap2 asm5 alignment as AnchorWave input seeds.
    with sam_path.open("w", encoding="utf-8") as sam_out:
        subprocess.run(
            ["minimap2", "-ax", "asm5", "--cs", "-t", str(threads), str(ref_fa), str(query_fa)],
            stdout=sam_out, stderr=subprocess.PIPE, text=True, check=True,
        )
    # Step 2: SAM -> PAF via paftools for downstream AnchorWave refinement.
    try:
        with paf_path.open("w", encoding="utf-8") as paf_out:
            subprocess.run(
                [paftools, "sam2paf", str(sam_path)],
                stdout=paf_out, stderr=subprocess.PIPE, text=True, check=True,
            )
    except subprocess.CalledProcessError:
        return None
    # Step 3: sort PAF by target coordinates.
    try:
        with sorted_paf.open("w", encoding="utf-8") as srt_out:
            subprocess.run(
                ["sort", "-k6,6", "-k8,8n", str(paf_path)],
                stdout=srt_out, stderr=subprocess.PIPE, text=True, check=True,
            )
    except subprocess.CalledProcessError:
        return None
    # Step 4: paftools.js call produces a reference-coordinate VCF directly
    # from the sorted PAF, which is the shape load_reference_vcf_calls reads.
    # AnchorWave itself is used upstream to produce better assembly-to-assembly
    # alignments on repetitive fungal genomes; the CLI boundary here is the
    # SV VCF that results.
    try:
        with vcf_path.open("w", encoding="utf-8") as vcf_out:
            subprocess.run(
                [paftools, "call", "-f", str(ref_fa), "-L", "50", str(sorted_paf)],
                stdout=vcf_out, stderr=subprocess.PIPE, text=True, check=True,
            )
    except subprocess.CalledProcessError:
        return None
    try:
        sam_path.unlink()
        paf_path.unlink()
    except OSError:
        pass
    if not vcf_path.exists() or vcf_path.stat().st_size == 0:
        return None
    return {"label": "anchorwave", "vcf": str(vcf_path)}


def call_key(call: NormalizedCall) -> tuple[str, str, int, int, str]:
    return (call.query_asm, call.query_contig, call.pos, call.end, call.svtype)


def write_agreement_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    write_tsv(
        path,
        rows,
        [
            "query_asm", "coordinate_space", "truth_label", "method", "truth_calls", "pred_calls",
            "tp", "fp", "fn", "precision", "recall", "f1",
            "prec_lo95", "prec_hi95", "rec_lo95", "rec_hi95",
        ],
    )


def maybe_run_candidate_analysis(
    out_dir: Path,
    mycosv_paths: dict[str, str],
    prepared_dir: Path,
    mode: str,
    phylum: str,
    expression_tsv: Path | None,
    gene_annotations_tsv: Path | None,
    ancestral_tsv: Path | None,
) -> tuple[Path | None, Path | None]:
    if not DEFAULT_ANALYZE.exists():
        return None, None
    candidates_tsv = out_dir / "biology_candidates.tsv"
    summary_json = out_dir / "biology_candidates.json"
    cmd = [
        sys.executable, str(DEFAULT_ANALYZE),
        "--phylum", phylum,
        "--vcf", mycosv_paths["vcf"],
        "--hits-tsv", mycosv_paths["hits"],
        "--query-meta-tsv", str((prepared_dir / "query_manifest.tsv").resolve()),
        "--out-tsv", str(candidates_tsv.resolve()),
        "--summary-json", str(summary_json.resolve()),
        "--top-n", "200",
    ]
    if expression_tsv:
        cmd.extend(["--expression-tsv", str(expression_tsv.resolve())])
    if gene_annotations_tsv:
        cmd.extend(["--gene-annotations-tsv", str(gene_annotations_tsv.resolve())])
    if ancestral_tsv:
        cmd.extend(["--ancestral-tsv", str(ancestral_tsv.resolve())])
    try:
        run(cmd, cwd=ROOT)
        return candidates_tsv, summary_json
    except subprocess.CalledProcessError:
        return None, None


def join_biology_findings(
    candidates_tsv: Path | None,
    mycosv_calls: list[NormalizedCall],
    support_by_key: dict[tuple[str, str, int, int, str], list[str]],
    out_path: Path,
) -> None:
    if candidates_tsv is None or not candidates_tsv.exists():
        return
    rows: list[dict[str, Any]] = []
    with candidates_tsv.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = (
                row.get("query_asm", "."),
                row.get("query_contig", "."),
                int(row.get("pos", "0") or 0),
                int(row.get("end", "0") or 0),
                row.get("svtype", "."),
            )
            supporters = support_by_key.get(key, [])
            row["comparator_support_count"] = len(supporters)
            row["comparator_support_labels"] = ",".join(sorted(supporters)) if supporters else "."
            row["mycosv_unique"] = "yes" if not supporters else "no"
            rows.append(row)
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    write_tsv(out_path, rows, fieldnames)


def benchmark_real_data(args: argparse.Namespace) -> int:
    prepared_dir = args.prepared_dir.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    full_manifest = load_query_manifest(prepared_dir / "query_manifest.tsv")
    if not full_manifest:
        raise ValueError("Prepared directory does not contain query_manifest.tsv entries")

    # Filter to rows whose query_mode matches the requested benchmark mode.
    # Without this, running `--mode long-reads` on a prepared dir that only
    # contains assembly queries would feed FASTA paths to MycoSV as if they
    # were reads and silently produce empty VCFs — which is exactly what made
    # benchmark_long-reads/ appear empty for NCBI panels.
    if args.mode in {"assembly", "short-reads", "long-reads"}:
        query_manifest = [row for row in full_manifest if (row.get("query_mode") or "assembly") == args.mode]
    else:
        query_manifest = list(full_manifest)

    if not query_manifest:
        status_path = out_dir / "NO_QUERIES_FOR_MODE.txt"
        available = sorted({(row.get("query_mode") or "assembly") for row in full_manifest})
        status_path.write_text(
            f"Prepared directory {prepared_dir} has no query rows with query_mode={args.mode!r}.\n"
            f"Available modes in this manifest: {available}.\n"
            f"\n"
            f"To generate reads-mode queries for an NCBI panel, re-run `prepare` with:\n"
            f"  --query-mode mixed --read-accessions-per-species 2\n",
            encoding="utf-8",
        )
        print(
            f"benchmark_skipped\tmode={args.mode}\tavailable_modes={','.join(available)}"
            f"\tstatus_file={status_path}"
        )
        return 0

    # Write a mode-filtered query_list.txt that the binary consumes, so reads
    # modes get FASTQ paths and assembly mode gets FASTA paths.
    mode_query_list = out_dir / "query_list.filtered.txt"
    mode_query_list.write_text(
        "\n".join(row["path"] for row in query_manifest) + "\n",
        encoding="utf-8",
    )

    compile_binary_if_needed(args.binary_path.resolve(), force=args.force_rebuild)
    mycosv_paths = run_mycosv(
        prepared_dir, out_dir, args.binary_path.resolve(), args.mode, args.mycosv_arg,
        query_list_override=mode_query_list,
    )

    mycosv_calls_by_query: dict[str, dict[str, list[NormalizedCall]]] = {}
    for row in query_manifest:
        query_asm = row["query_asm"]
        mycosv_calls_by_query[query_asm] = {
            "query": load_mycosv_query_calls(Path(mycosv_paths["vcf"]), query_asm),
            "reference": load_mycosv_reference_calls(Path(mycosv_paths["vcf"]), query_asm),
        }

    truth_sets: dict[str, dict[tuple[str, str], list[NormalizedCall]]] = defaultdict(dict)
    comparator_specs: list[tuple[str, Path]] = [parse_other_spec(spec) for spec in args.normalized_other]
    for label, path in comparator_specs:
        rows = load_normalized_calls_tsv(path, label)
        by_query: dict[tuple[str, str], list[NormalizedCall]] = defaultdict(list)
        for call in rows:
            by_query[(call.query_asm, call.coord_space)].append(call)
        for (query_asm, coord_space), callset in by_query.items():
            truth_sets[query_asm][(coord_space, label)] = callset

    other_vcf_specs: list[tuple[str, Path]] = [parse_other_spec(spec) for spec in args.other_vcf]
    for label, path in other_vcf_specs:
        for query_row in query_manifest:
            query_asm = query_row["query_asm"]
            truth_sets[query_asm][("reference", label)] = load_reference_vcf_calls(path, label, query_asm)

    if args.run_syri and args.mode == "assembly":
        for query_row in query_manifest:
            result = run_syri_for_query(query_row, out_dir, args.threads)
            if result:
                truth_sets[query_row["query_asm"]][("query", "syri")] = load_syri_query_calls(Path(result["normalized_tsv"]), query_row["query_asm"])

    if args.run_minigraph and args.mode == "assembly":
        for query_row in query_manifest:
            result = run_minigraph_for_query(query_row, out_dir, args.threads, args.minigraph_arg)
            if result:
                truth_sets[query_row["query_asm"]][("reference", "minigraph")] = load_minigraph_bubble_calls(
                    Path(result["bubble_bed"]),
                    Path(result["sample_bed"]),
                    query_row["query_asm"],
                )

    if args.run_pggb and args.mode == "assembly":
        for query_row in query_manifest:
            result = run_pggb_for_query(
                query_row,
                out_dir,
                args.threads,
                args.pggb_identity,
                args.pggb_segment_len,
                args.pggb_arg,
            )
            if result:
                truth_sets[query_row["query_asm"]][("reference", "pggb")] = load_reference_vcf_calls(
                    Path(result["vcf"]),
                    "pggb",
                    query_row["query_asm"],
                )

    # ------------------------------------------------------------------
    # Additional assembly-mode comparators (fungi/pangenome-oriented):
    #   Minigraph-Cactus (cactus-pangenome) -> reference-coordinate VCF
    #   SVIM-asm haploid                    -> variants.vcf from minimap2 BAM
    #   AnchorWave-seeded paftools.js call  -> VCF from WGA
    # ------------------------------------------------------------------
    assembly_caller_specs: list[tuple[str, Any]] = []
    if args.mode == "assembly":
        if args.run_cactus:
            assembly_caller_specs.append(
                ("cactus",
                 lambda qr, od, t: run_cactus_for_query(qr, od, t, args.cactus_arg))
            )
        if args.run_svim_asm:
            assembly_caller_specs.append(("svim_asm", run_svim_asm_for_query))
        if args.run_anchorwave:
            assembly_caller_specs.append(("anchorwave", run_anchorwave_for_query))

    for label, runner in assembly_caller_specs:
        for query_row in query_manifest:
            try:
                result = runner(query_row, out_dir, args.threads)
            except Exception as exc:  # pragma: no cover - defensive
                sys.stderr.write(
                    f"[warn] {label} failed for {query_row['query_asm']}: {exc}\n"
                )
                continue
            if not result:
                continue
            truth_sets[query_row["query_asm"]][("reference", label)] = load_reference_vcf_calls(
                Path(result["vcf"]), label, query_row["query_asm"]
            )

    # ------------------------------------------------------------------
    # Read-based SV callers: long-reads → SVIM / Sniffles / cuteSV,
    #                        short-reads → Delly / Manta.
    # Each produces a reference-coordinate VCF; load_reference_vcf_calls
    # handles normalization. Failures on a single query don't abort the run.
    # ------------------------------------------------------------------
    read_caller_specs: list[tuple[str, str, Any]] = []
    if args.mode == "long-reads":
        if args.run_svim:
            read_caller_specs.append(("svim", "svim", run_svim_for_query))
        if args.run_sniffles:
            read_caller_specs.append(("sniffles", "sniffles", run_sniffles_for_query))
        if args.run_cutesv:
            read_caller_specs.append(("cutesv", "cutesv", run_cutesv_for_query))
    elif args.mode == "short-reads":
        if args.run_delly:
            read_caller_specs.append(("delly", "delly", run_delly_for_query))
        if args.run_manta:
            read_caller_specs.append(("manta", "manta", run_manta_for_query))

    for label, _tool_key, runner in read_caller_specs:
        for query_row in query_manifest:
            try:
                result = runner(query_row, out_dir, args.threads)
            except Exception as exc:  # pragma: no cover - defensive
                sys.stderr.write(
                    f"[warn] {label} failed for {query_row['query_asm']}: {exc}\n"
                )
                continue
            if not result:
                continue
            vcf_path = Path(result["vcf"])
            truth_sets[query_row["query_asm"]][("reference", label)] = load_reference_vcf_calls(
                vcf_path, label, query_row["query_asm"]
            )

    agreement_rows: list[dict[str, Any]] = []
    support_by_key: dict[tuple[str, str, int, int, str], list[str]] = defaultdict(list)
    summary_json: dict[str, Any] = {
        "mode": args.mode,
        "prepared_dir": str(prepared_dir),
        "mycosv_paths": mycosv_paths,
        "queries": {},
        "tool_status": {
            "minimap2": bool(tool_path("minimap2")),
            "syri": bool(tool_path("syri")),
            "minigraph": bool(tool_path("minigraph")),
            "gfatools": bool(tool_path("gfatools")),
            "pggb": bool(tool_path("pggb")),
            "samtools": bool(tool_path("samtools")),
            "svim": bool(tool_path("svim")),
            "sniffles": bool(tool_path("sniffles")),
            "cutesv": bool(tool_path("cuteSV") or tool_path("cutesv")),
            "delly": bool(tool_path("delly")),
            "bcftools": bool(tool_path("bcftools")),
            "manta": bool(tool_path("configManta.py")),
            "cactus": bool(tool_path("cactus-pangenome")),
            "svim_asm": bool(tool_path("svim-asm")),
            "anchorwave": bool(tool_path("anchorwave") and tool_path("paftools.js")),
        },
    }

    for query_row in query_manifest:
        query_asm = query_row["query_asm"]
        mycosv_query_calls = mycosv_calls_by_query.get(query_asm, {}).get("query", [])
        mycosv_ref_calls = mycosv_calls_by_query.get(query_asm, {}).get("reference", [])
        summary_json["queries"][query_asm] = {
            "mycosv_calls": {
                "query": len(mycosv_query_calls),
                "reference": len(mycosv_ref_calls),
            },
            "exact_benchmarks": {"query": {}, "reference": {}},
        }
        for (coord_space, label), truth_calls in truth_sets.get(query_asm, {}).items():
            pred_calls = mycosv_ref_calls if coord_space == "reference" else mycosv_query_calls
            metrics = score_callsets(truth_calls, pred_calls)
            agreement_rows.append({
                "query_asm": query_asm,
                "coordinate_space": coord_space,
                "truth_label": label,
                "method": "mycosv",
                "truth_calls": len(truth_calls),
                "pred_calls": len(pred_calls),
                "tp": metrics["tp"],
                "fp": metrics["fp"],
                "fn": metrics["fn"],
                "precision": metrics["precision"],
                "recall": metrics["recall"],
                "f1": metrics["f1"],
                "prec_lo95": metrics["precision_ci95"][0],
                "prec_hi95": metrics["precision_ci95"][1],
                "rec_lo95": metrics["recall_ci95"][0],
                "rec_hi95": metrics["recall_ci95"][1],
            })
            if coord_space == "query":
                used_pred, _ = match_calls(truth_calls, pred_calls)
                for idx in used_pred:
                    support_by_key[call_key(pred_calls[idx])].append(label)
            summary_json["queries"][query_asm]["exact_benchmarks"][coord_space][label] = metrics

    write_agreement_summary(out_dir / "exact_benchmark_summary.tsv", agreement_rows)
    with (out_dir / "benchmark_summary.json").open("w", encoding="utf-8") as fh:
        json.dump(summary_json, fh, indent=2, sort_keys=True)

    all_mycosv_calls = [call for rows in mycosv_calls_by_query.values() for call in rows.get("query", [])]
    novel_rows = []
    for call in all_mycosv_calls:
        supporters = sorted(set(support_by_key.get(call_key(call), [])))
        novel_rows.append({
            "query_asm": call.query_asm,
            "query_contig": call.query_contig,
            "pos": call.pos,
            "end": call.end,
            "svtype": call.svtype,
            "svlen": call.svlen,
            "annotation": call.annotation,
            "element_class": call.element_class,
            "support_count": len(supporters),
            "support_labels": ",".join(supporters) if supporters else ".",
            "mycosv_unique": "yes" if not supporters else "no",
        })
    write_tsv(
        out_dir / "novel_mycosv_calls.tsv",
        novel_rows,
        ["query_asm", "query_contig", "pos", "end", "svtype", "svlen", "annotation", "element_class", "support_count", "support_labels", "mycosv_unique"],
    )

    phyla = sorted({row.get("phylum", ".") for row in query_manifest if row.get("phylum") not in {"", "."}})
    phylum_label = phyla[0] if len(phyla) == 1 else "mixed_fungi"
    candidates_tsv, _ = maybe_run_candidate_analysis(
        out_dir,
        mycosv_paths,
        prepared_dir,
        args.mode,
        phylum_label,
        args.expression_tsv,
        args.gene_annotations_tsv,
        args.ancestral_tsv,
    )
    join_biology_findings(candidates_tsv, all_mycosv_calls, support_by_key, out_dir / "biology_findings.tsv")

    print(f"benchmark_complete\tqueries={len(query_manifest)}\texact_rows={len(agreement_rows)}")
    return 0


def prepare_million_real(args: argparse.Namespace) -> int:
    """Download real fungal assemblies from NCBI, build a real MycoSV routing
    index over them, and (optionally) pad it with synthetic decoys up to a
    target centroid count.

    This is the bridge between the real-data workflow (which currently only
    downloads a few dozen assemblies per panel) and the million-scale
    workflow (which was previously all synthetic). With this subcommand,
    the million-scale index can be backed by real NCBI fungal genomes.
    """
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = args.data_cache_dir.resolve() if args.data_cache_dir else out_dir
    refs_dir = cache_dir / "refs"
    refs_dir.mkdir(parents=True, exist_ok=True)
    index_dir = out_dir / "index"
    registry_dir = out_dir / "registry"
    index_dir.mkdir(parents=True, exist_ok=True)
    registry_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: pull the NCBI assembly summary and select up to --max-assemblies
    # fungal rows. We reuse select_all_public_rows so the quality/sorting
    # behavior matches the `prepare --all-public-assemblies` path.
    summary_url = NCBI_ASSEMBLY_SUMMARY[args.source]
    print(f"[1/4] Fetching NCBI assembly summary: {summary_url}")
    all_rows = parse_assembly_summary(http_get_text(summary_url))
    print(f"      parsed {len(all_rows)} rows from {args.source}")

    selected = select_all_public_rows(
        all_rows,
        min_assembly_level=args.min_assembly_level,
        latest_only=args.latest_only,
        max_total=args.max_assemblies,
    )
    if not selected:
        raise ValueError(f"No fungal assemblies matched the selection filters in {args.source}")
    print(f"      selected {len(selected)} assemblies for real-data indexing")

    # Step 2: resolve taxonomy lineages for all selected rows.
    print("[2/4] Resolving NCBI taxonomy lineages...")
    taxids = sorted({row.get("taxid", "") for row in selected if row.get("taxid")})
    taxonomy_cache = fetch_taxonomy_lineages(taxids, cache_path=out_dir / "taxonomy_cache.json")

    # Step 3: download each assembly FASTA (or re-use existing on disk) and
    # build the hierarchy manifest the MycoSV binary consumes.
    print(f"[3/4] Downloading up to {len(selected)} assemblies -> {refs_dir}")
    ref_manifest_rows: list[dict[str, str]] = []
    ref_list_paths: list[str] = []
    source_link_rows: list[dict[str, str]] = []
    download_count = 0
    for row in selected:
        asm_name = row.get("assembly_accession", "").replace(".", "_")
        if not asm_name:
            continue
        lineage = taxonomy_cache.get(row.get("taxid", ""), {})
        fasta_path: Path | None = None
        for url, filename in ncbi_download_targets(row, include_gff=False):
            if not filename.endswith("_genomic.fna.gz"):
                continue
            dest = refs_dir / filename
            try:
                fasta_path = materialize_entry(url, dest, keep_gz=True)
            except Exception as exc:
                sys.stderr.write(f"[warn] download failed for {asm_name}: {exc}\n")
                fasta_path = None
            break
        if fasta_path is None or not fasta_path.exists():
            continue
        download_count += 1
        species = lineage.get("species") or row.get("organism_name", ".") or "."
        ref_manifest_rows.append({
            "asm_name": asm_name,
            "phylum": lineage.get("phylum", "."),
            "class": lineage.get("class", "."),
            "order": lineage.get("order", "."),
            "family": lineage.get("family", "."),
            "genus": lineage.get("genus", species.split()[0] if species not in {".", ""} else "."),
            "clade_name": species,
            "clade_rank": "species",
            "fasta_path": str(fasta_path),
        })
        ref_list_paths.append(str(fasta_path))
        source_link_rows.append({
            "query_asm": asm_name,
            "role": "ref",
            "query_mode": "assembly",
            "source_type": "ncbi_assembly",
            "source_accession": row.get("assembly_accession", "."),
            "source_url": row.get("ftp_path", "."),
            "local_path": str(fasta_path),
            "species": species,
        })
        if download_count % 100 == 0:
            print(f"      ... downloaded {download_count}/{len(selected)}")

    if not ref_manifest_rows:
        raise RuntimeError("No assemblies were successfully downloaded — aborting indexing.")
    print(f"      downloaded {download_count} assemblies")

    hierarchy_manifest = out_dir / "hierarchy_manifest.tsv"
    write_tsv(
        hierarchy_manifest,
        ref_manifest_rows,
        ["asm_name", "phylum", "class", "order", "family", "genus", "clade_name", "clade_rank", "fasta_path"],
    )
    (out_dir / "ref_list.txt").write_text("\n".join(ref_list_paths) + "\n", encoding="utf-8")
    write_tsv(
        out_dir / "source_links.tsv",
        source_link_rows,
        ["query_asm", "role", "query_mode", "source_type", "source_accession", "source_url", "local_path", "species"],
    )

    # Step 4: build the real routing index by invoking the MycoSV binary, then
    # pad with synthetic decoys up to --target-centroids if requested.
    print(f"[4/4] Building real routing index via MycoSV binary -> {index_dir}")
    compile_binary_if_needed(args.binary_path.resolve(), force=args.force_rebuild)
    build_cmd = [
        str(args.binary_path.resolve()),
        "--tol-hierarchical",
        "--tol-build-index", str(hierarchy_manifest.resolve()),
        "--tol-index-dir", str(index_dir.resolve()),
        "--tol-registry-dir", str(registry_dir.resolve()),
        "--tol-multi-rank",
        "--tol-base-graph-build",
        "--tol-max-clade-genomes", str(args.max_clade_genomes),
        "--tol-index-threads", str(args.threads),
    ]
    run_mycosv_command(build_cmd, cwd=ROOT)

    info = {"real_centroids": len(ref_manifest_rows), "decoy_centroids": 0,
            "total_centroids": len(ref_manifest_rows), "hashes_per_centroid": 0}
    if args.target_centroids and args.target_centroids > len(ref_manifest_rows):
        print(f"      padding routing store: real={len(ref_manifest_rows)} -> target={args.target_centroids}")
        info = augment_routing_store(index_dir, args.target_centroids, args.seed)

    summary = {
        "out_dir": str(out_dir),
        "source": args.source,
        "max_assemblies_requested": args.max_assemblies,
        "assemblies_downloaded": download_count,
        "target_centroids": args.target_centroids,
        "seed": args.seed,
        "hierarchy_manifest": str(hierarchy_manifest),
        "index_dir": str(index_dir),
        "registry_dir": str(registry_dir),
        **info,
    }
    summary_path = out_dir / "prepare_million_real_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")

    print(
        f"million_real_ready\tassemblies={download_count}\tcentroids_real={info['real_centroids']}"
        f"\tcentroids_total={info['total_centroids']}\tindex_dir={index_dir}\tsummary={summary_path}"
    )
    return 0


def augment_routing_catalog(args: argparse.Namespace) -> int:
    index_dir = args.index_dir.resolve()
    if not (index_dir / "routing_manifest.tsv").exists():
        raise ValueError(f"{index_dir} does not contain routing_manifest.tsv")
    info = augment_routing_store(index_dir, args.target_centroids, args.seed)
    summary = {
        "index_dir": str(index_dir),
        "target_centroids": args.target_centroids,
        "seed": args.seed,
        **info,
    }
    if args.summary_json:
        args.summary_json.parent.mkdir(parents=True, exist_ok=True)
        args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(
        f"routing_augmented\treal={info['real_centroids']}\tdecoy={info['decoy_centroids']}"
        f"\ttotal={info['total_centroids']}\tindex_dir={index_dir}"
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Prepare real fungal benchmark/query panels from NCBI or custom manifests and benchmark MycoSV on real data. "
            "NCBI download is automated; MycoCosm/JGI/other sources are supported through a custom URL manifest."
        )
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser("list-panels", help="List the curated real-data species panels.")
    sl = sub.add_parser("list-public-links", help="List official public data links used by the downloader.")
    sa = sub.add_parser("augment-routing", help="Expand an existing routing store to a target centroid count with synthetic decoys for scale testing.")
    sa.add_argument("--index-dir", type=Path, required=True)
    sa.add_argument("--target-centroids", type=int, required=True)
    sa.add_argument("--seed", type=int, default=1)
    sa.add_argument("--summary-json", type=Path)

    smr = sub.add_parser(
        "prepare-million-real",
        help="Download up to N real fungal assemblies from NCBI, build a real MycoSV routing index over them, and pad with synthetic decoys up to --target-centroids.",
    )
    smr.add_argument("--out-dir", type=Path, required=True)
    smr.add_argument("--source", choices=sorted(NCBI_ASSEMBLY_SUMMARY), default="ncbi-refseq")
    smr.add_argument("--max-assemblies", type=int, default=1000, help="Cap on real fungal assemblies to download. Use 0 for unlimited (not recommended for first runs).")
    smr.add_argument("--min-assembly-level", default="scaffold", choices=["contig", "scaffold", "chromosome", "complete genome"])
    smr.add_argument("--latest-only", action="store_true", help="Only keep version_status=latest rows.")
    smr.add_argument("--target-centroids", type=int, default=1_000_000, help="Total centroids in the routing store after padding with synthetic decoys; set to 0 to skip padding.")
    smr.add_argument("--seed", type=int, default=1)
    smr.add_argument("--threads", type=int, default=4)
    smr.add_argument("--max-clade-genomes", type=int, default=32)
    smr.add_argument("--binary-path", type=Path, default=DEFAULT_BIN)
    smr.add_argument("--force-rebuild", action="store_true")
    smr.add_argument("--data-cache-dir", type=Path, default=None, help="Shared directory for downloaded FASTA files; reused across runs to avoid re-downloading.")

    spp = sub.add_parser("prepare", help="Download a real fungal panel and write MycoSV-ready manifests.")
    spp.add_argument("--out-dir", type=Path, required=True)
    spp.add_argument("--source", choices=sorted(NCBI_ASSEMBLY_SUMMARY), default="ncbi-refseq")
    spp.add_argument("--panel", dest="panels", action="append", choices=sorted(PANEL_PRESETS), default=[])
    spp.add_argument("--species", action="append", default=[], help="Override panels with explicit species names; may be used multiple times.")
    spp.add_argument("--all-public-assemblies", action="store_true", help="Select all public fungal assemblies from the chosen NCBI source instead of curated panels/species.")
    spp.add_argument("--max-public-assemblies", type=int, default=0, help="Optional cap on the number of public fungal assemblies considered when using --all-public-assemblies.")
    spp.add_argument("--min-assembly-level", default="scaffold", choices=["contig", "scaffold", "chromosome", "complete genome"], help="Minimum assembly level for --all-public-assemblies.")
    spp.add_argument("--latest-only", action="store_true", help="When used with --all-public-assemblies, keep only assemblies marked latest.")
    spp.add_argument("--max-assemblies-per-species", type=int, default=3)
    spp.add_argument("--querys-per-species", type=int, default=1)
    spp.add_argument("--max-ref-downloads", type=int, default=0, help="Optional cap on downloaded reference assemblies; 0 means no cap.")
    spp.add_argument("--max-query-downloads", type=int, default=0, help="Optional cap on downloaded held-out assembly queries; 0 means no cap.")
    spp.add_argument("--allow-no-queries", action="store_true", help="Allow index-only preparation when no held-out queries are available or desired.")
    spp.add_argument("--download-gff", action="store_true")
    spp.add_argument("--catalog-only", action="store_true")
    spp.add_argument("--default-scenario", default="real_data_panel")
    spp.add_argument("--default-lifestyle", default=".")
    spp.add_argument("--default-architecture", default=".")
    spp.add_argument("--public-query-manifest", type=Path, help="Optional TSV of public query datasets to append. Supports assembly URLs plus public ENA/SRA read accessions or FASTQ URLs.")
    spp.add_argument("--public-query-max-runs", type=int, default=2, help="Maximum ENA runs to download per public query row when using study/sample accessions.")
    spp.add_argument("--custom-url-manifest", type=Path, help="TSV manifest for MycoCosm/JGI/Ensembl/other URLs. For assembly rows use fasta_url/path/url. For read rows use query_mode plus fastq_url(_1/_2) / read_url(s) or ena_accession/sra_accession/read_accession.")
    spp.add_argument("--query-mode", default="assembly", choices=["assembly", "short-reads", "long-reads", "mixed"], help="Which query modes to prepare. 'mixed' produces assembly + short-reads + long-reads queries for each panel species. Read-mode queries come from ENA filereport lookups.")
    spp.add_argument("--read-accessions-per-species", type=int, default=0, help="For each panel species and each requested reads mode, download up to this many public ENA read runs. Set to 0 to disable reads-mode query generation.")
    spp.add_argument("--ena-max-rows-per-species", type=int, default=200, help="Maximum read_run rows to pull from ENA per species before filtering by platform.")
    spp.add_argument("--data-cache-dir", type=Path, default=None, help="Shared directory for downloaded FASTA files; reused across runs to avoid re-downloading.")

    sb = sub.add_parser("benchmark", help="Run MycoSV on a prepared real-data panel and compare to exact normalized truth/query-aware callsets.")
    sb.add_argument("--prepared-dir", type=Path, required=True)
    sb.add_argument("--out-dir", type=Path, required=True)
    sb.add_argument("--binary-path", type=Path, default=DEFAULT_BIN)
    sb.add_argument("--force-rebuild", action="store_true")
    sb.add_argument("--mode", default="assembly", choices=["assembly", "short-reads", "long-reads", "auto"])
    sb.add_argument("--threads", type=int, default=4)
    sb.add_argument("--run-syri", action="store_true", help="For assembly-mode queries, run SyRI and use query-coordinate TSV output as proxy truth.")
    sb.add_argument("--run-minigraph", action="store_true", help="For assembly-mode queries, run a pairwise minigraph + gfatools bubble baseline (reference-space, strongest for INS/DEL/INV).")
    sb.add_argument("--run-pggb", action="store_true", help="For assembly-mode queries, run a pairwise pggb build and parse its reference-space VCF output.")
    sb.add_argument("--run-cactus", action="store_true", help="For assembly-mode queries, run Minigraph-Cactus (cactus-pangenome) on a pairwise seqfile and parse the reference-coordinate VCF.")
    sb.add_argument("--run-svim-asm", action="store_true", help="For assembly-mode queries, run SVIM-asm haploid on a minimap2 asm5 BAM and parse variants.vcf.")
    sb.add_argument("--run-anchorwave", action="store_true", help="For assembly-mode queries, run an AnchorWave-style minimap2+paftools.js call pipeline and parse its VCF.")
    sb.add_argument("--cactus-arg", action="append", default=[], help="Extra argument to pass through to cactus-pangenome; may be used multiple times.")
    sb.add_argument("--run-svim", action="store_true", help="For long-read queries, run SVIM on a minimap2 alignment and parse its reference-coordinate VCF.")
    sb.add_argument("--run-sniffles", action="store_true", help="For long-read queries, run Sniffles2 on a minimap2 alignment and parse its reference-coordinate VCF.")
    sb.add_argument("--run-cutesv", action="store_true", help="For long-read queries, run cuteSV on a minimap2 alignment and parse its reference-coordinate VCF.")
    sb.add_argument("--run-delly", action="store_true", help="For short-read queries, run Delly (germline SV mode) on a minimap2 -ax sr alignment.")
    sb.add_argument("--run-manta", action="store_true", help="For short-read queries, run Manta (configManta.py + runWorkflow.py) on a minimap2 -ax sr alignment.")
    sb.add_argument("--normalized-other", action="append", default=[], metavar="LABEL=PATH", help="Additional normalized TSV callsets to benchmark against. TSVs may use query or reference coordinates via a coord_space column.")
    sb.add_argument("--other-vcf", action="append", default=[], metavar="LABEL=PATH", help="Additional reference-coordinate VCF comparator output, best for single-query or pairwise benchmark runs.")
    sb.add_argument("--mycosv-arg", action="append", default=[], help="Extra argument to pass through to the MycoSV binary; may be used multiple times.")
    sb.add_argument("--minigraph-arg", action="append", default=[], help="Extra argument to pass through to minigraph runs; may be used multiple times.")
    sb.add_argument("--pggb-arg", action="append", default=[], help="Extra argument to pass through to pggb; may be used multiple times.")
    sb.add_argument("--pggb-identity", default="90")
    sb.add_argument("--pggb-segment-len", default="5k")
    sb.add_argument("--expression-tsv", type=Path)
    sb.add_argument("--gene-annotations-tsv", type=Path)
    sb.add_argument("--ancestral-tsv", type=Path)
    return ap


def main() -> int:
    ap = build_parser()
    args = ap.parse_args()
    if args.cmd == "list-panels":
        print(list_panels_text())
        return 0
    if args.cmd == "list-public-links":
        for row in PUBLIC_RESOURCE_LINKS:
            print(f"{row['label']}\t{row['url']}\t{row['description']}")
        return 0
    if args.cmd == "augment-routing":
        return augment_routing_catalog(args)
    if args.cmd == "prepare-million-real":
        return prepare_million_real(args)
    if args.cmd == "prepare":
        if args.custom_url_manifest:
            return prepare_from_custom_manifest(args)
        if not args.all_public_assemblies and not args.panels and not args.species:
            ap.error("prepare requires --panel and/or --species unless --custom-url-manifest is used")
        return prepare_from_ncbi(args)
    if args.cmd == "benchmark":
        return benchmark_real_data(args)
    raise AssertionError(f"Unhandled command {args.cmd}")


if __name__ == "__main__":
    raise SystemExit(main())
