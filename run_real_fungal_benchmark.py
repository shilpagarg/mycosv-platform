#!/usr/bin/env python3
# Designed for Linux

from __future__ import annotations

import argparse
import bisect
import csv
import ctypes
import gzip
import io
import json
import math
import os
import random
import re
import shutil
import statistics
import subprocess
import sys
import threading
import time
import uuid
import hashlib
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

from run_million_mode_query_benchmark import augment_routing_store
from sv_pr_utils import DEFAULT_TOL_BP, DEFAULT_TOL_LEN_FRAC, expand_to_multisample_vcf, wilson_ci


ROOT = Path(__file__).resolve().parent
DEFAULT_BIN = ROOT / "fungi_graphsv_tol_bin"
DEFAULT_ANALYZE = ROOT / "analyze_new_biology_candidates.py"
MYCOSV_BRIDGE_CPP = ROOT / "mycosv_cli_bridge.cpp"
DEFAULT_DATA_CACHE = ROOT / "data_cache"

NCBI_ASSEMBLY_SUMMARY = {
    "ncbi-refseq": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt",
    "ncbi-genbank": "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/assembly_summary.txt",
}
NCBI_BEST_SOURCE = "ncbi-best"
NCBI_SOURCE_CHOICES = sorted([*NCBI_ASSEMBLY_SUMMARY, NCBI_BEST_SOURCE])

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
        "label": "ensembl_fungi_ftp",
        "url": "https://fungi.ensembl.org/info/data/ftp/index.html",
        "description": "Ensembl Fungi public FTP downloads for FASTA, GTF/GFF3, GenBank and other annotation files",
    },
    {
        "label": "ncbi_datasets",
        "url": "https://www.ncbi.nlm.nih.gov/datasets/",
        "description": "NCBI Datasets genome packages and programmatic assembly download API",
    },
    {
        "label": "mycocosm",
        "url": "https://mycocosm.jgi.doe.gov/",
        "description": "JGI MycoCosm fungal genome portal",
    },
    # Gene/expression data sources used by analyze_new_biology_candidates.py
    # when populating expression_supported / expression_log2_fc / expression_padj.
    # Operators populate expression.tsv from these manually for now (no
    # canonical species->experiment mapping for fungi), and prepare-step picks
    # it up automatically as prepared_dir/expression.tsv on the next benchmark.
    {
        "label": "ensembl_fungi_rest",
        "url": "https://rest.ensembl.org/documentation/info/lookup",
        "description": "Ensembl Fungi REST: gene + ortholog metadata for annotated fungal species",
    },
    {
        "label": "ensembl_fungi_ftp",
        "url": "https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/current/gff3/",
        "description": "Ensembl Fungi current GFF3 dumps (fallback when NCBI has no GFF for an assembly)",
    },
    {
        "label": "expression_atlas_baseline",
        "url": "https://www.ebi.ac.uk/gxa/experiments?experimentType=baseline&kingdom=fungi",
        "description": "EBI Expression Atlas baseline experiments (per-tissue/condition gene expression for fungi)",
    },
    {
        "label": "expression_atlas_differential",
        "url": "https://www.ebi.ac.uk/gxa/experiments?experimentType=differential&kingdom=fungi",
        "description": "EBI Expression Atlas differential experiments (log2FC/padj across conditions; primary source for expression_log2_fc / expression_padj)",
    },
    {
        "label": "ncbi_geo",
        "url": "https://www.ncbi.nlm.nih.gov/geo/browse/?view=series&search=fungi",
        "description": "NCBI GEO public RNA-seq / microarray series (use SRA/ENA for raw, recount3/Atlas for processed)",
    },
    {
        "label": "fungidb",
        "url": "https://fungidb.org/fungidb/app",
        "description": "FungiDB / VEuPathDB: integrated genotype + phenotype + expression for pathogenic fungi",
    },
    {
        "label": "phi_base",
        "url": "http://www.phi-base.org/",
        "description": "PHI-base: pathogen-host interaction phenotypes; gene-level pathogenicity calls",
    },
    {
        "label": "fungaltraits",
        "url": "https://github.com/traitecoevo/fungaltraits",
        "description": "FungalTraits: curated genus-/species-level lifestyle and trait database (Polõme 2020)",
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
    "read_count",
    "base_count",
    "submitted_ftp",
    "submitted_md5",
    "submitted_bytes",
]

# Minimum read count for an ENA run to be considered for SV calling. Below this
# the run is almost certainly a sentinel/test upload (e.g. SRR33624766 reported
# 1 read) and would fail every comparator with "no alignments". Picked at 1 000
# so a single MiSeq lane (~10⁶ reads) trivially passes while obvious junk like
# "1 read" / "1000 reads" deposits don't get selected and waste a download.
_ENA_MIN_READS = 1000

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

BIOLOGY_FINDINGS_EXTRA_FIELDS = [
    "comparator_support_count",
    "comparator_support_labels",
    "single_reference_equivalent",
    "mycosv_unique",
    "evidence_tier",
]

READ_VALIDATION_FIELDS = [
    "query_asm",
    "ref_contig",
    "pos",
    "end",
    "svtype",
    "source",
    "coord_space",
    "read_support",
    "validation_support",
    "support_source",
    "read_validated",
    "status",
]

# ── Long-read platform detection ───────────────────────────────────────────
# PacBio HiFi (CCS): Revio, Sequel IIe/II CCS — ≥99 % accuracy, 10–25 kb.
#   minimap2 preset:  map-hifi
#   SV callers:       sniffles2, cuteSV (HiFi-tuned cluster params), SVIM
# PacBio CLR (Sequel I, RS II): lower per-read accuracy.
#   minimap2 preset:  map-pb
# ONT R10.4.1 simplex: ~Q20 accuracy on PromethION / GridION / MinION Mk1C.
#   minimap2 preset:  map-ont  (same as R9.4.1; Sniffles2 --long-read-model
#                               ont_r10_q20 optional for v2.2+)
#   WhatsHap phasing: applicable for dikaryotic / diploid fungi (e.g. Puccinia,
#                     Leptosphaeria, Zymoseptoria) once SNP calls are available.
#
# Short reads (Illumina NovaSeq / HiSeq):
#   bwa-mem2 (or minimap2 -ax sr) → samtools sort+index → Delly / Manta.
#   bwa-mem2 is the gold standard for Delly/Manta BAM inputs; minimap2 sr is
#   used here for a single-tool dependency.

# Instrument model keywords indicating HiFi CCS output.
_HIFI_MODEL_KW: frozenset[str] = frozenset({"revio", "sequel iie", "sequel 2e"})
# ENA library_strategy values that directly assert CCS / HiFi processing.
_HIFI_STRATEGY_KW: frozenset[str] = frozenset({"hifi", "ccs", "hi-fi"})


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
    read_support: int | None = None
    mate_contig: str = "."
    mate_pos: int = 0
    mate_end: int = 0


_TOOL_TIMEOUT = int(os.environ.get("MYCOSV_TOOL_TIMEOUT", "14400"))
# Soft cap applied to subprocess calls inside the per-query comparator runners
# (minigraph / svim_asm / svim / sniffles / cuteSV / delly / manta / anchorwave
# / cactus and their shared _minimap2_align_reads helper). Smaller than
# _TOOL_TIMEOUT so a single hung tool can't eat the whole mode budget when
# comparators are scheduled in parallel across queries; 45 min was empirically
# the wall time above which the bench almost always missed its 4 h cap.
# Override with MYCOSV_COMPARATOR_TIMEOUT (seconds) for long-walltime runs.
_COMPARATOR_TIMEOUT = int(os.environ.get("MYCOSV_COMPARATOR_TIMEOUT", "2700"))
_BIOLOGY_TIMEOUT = int(os.environ.get("MYCOSV_BIOLOGY_TIMEOUT", "120"))
_GENE_ANNOTATION_MAX_BYTES = int(os.environ.get(
    "MYCOSV_GENE_ANNOTATION_MAX_BYTES",
    str(512 * 1024 * 1024),
))


def _stderr_tail(exc: BaseException, max_lines: int = 8) -> str:
    """Last `max_lines` of captured stderr from a CalledProcessError, joined
    into a single newline-prefixed string for inclusion in `[warn]` lines.
    Returns "" if no stderr is available."""
    raw = getattr(exc, "stderr", None)
    if not raw:
        return ""
    if isinstance(raw, bytes):
        try:
            raw = raw.decode("utf-8", errors="replace")
        except Exception:
            raw = repr(raw)
    lines = [ln for ln in raw.splitlines() if ln.strip()]
    if not lines:
        return ""
    tail = lines[-max_lines:]
    return "\n  | " + "\n  | ".join(tail)


def _persist_stderr(work_dir: Path, label: str, exc: BaseException) -> Path | None:
    """Dump full stderr (and stdout, if non-empty) to a per-failure log so the
    exit-1/2 reason is recoverable after the run finishes. Returns the log
    path on success, None if nothing to write."""
    raw_err = getattr(exc, "stderr", None) or ""
    raw_out = getattr(exc, "stdout", None) or ""
    if isinstance(raw_err, bytes):
        raw_err = raw_err.decode("utf-8", errors="replace")
    if isinstance(raw_out, bytes):
        raw_out = raw_out.decode("utf-8", errors="replace")
    if not raw_err and not raw_out:
        return None
    try:
        work_dir.mkdir(parents=True, exist_ok=True)
        path = work_dir / f"{label}.stderr.log"
        with path.open("w", encoding="utf-8") as fh:
            if raw_err:
                fh.write("=== stderr ===\n")
                fh.write(raw_err)
                if not raw_err.endswith("\n"):
                    fh.write("\n")
            if raw_out:
                fh.write("=== stdout ===\n")
                fh.write(raw_out)
                if not raw_out.endswith("\n"):
                    fh.write("\n")
        return path
    except OSError:
        return None


def _log_comparator_failure(out_dir: Path, label: str, query_asm: str, reason: str) -> None:
    """Append a structured failure row so the visualization (and any operator
    grep) can see exactly which (query, comparator) pairs lost coverage and
    why, without parsing free-form stderr.
    """
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "comparator_failures.tsv"
        new_file = not path.exists()
        with path.open("a", encoding="utf-8") as fh:
            if new_file:
                fh.write("query_asm\tcomparator\treason\n")
            fh.write(f"{query_asm}\t{label}\t{reason.replace(chr(9), ' ').replace(chr(10), ' ')}\n")
    except OSError:
        pass


def run(
    cmd: list[str],
    cwd: Path | None = None,
    timeout: int = _TOOL_TIMEOUT,
    memory_limit_bytes: int | None = None,
) -> subprocess.CompletedProcess[str]:
    preexec_fn = None
    if memory_limit_bytes and memory_limit_bytes > 0 and os.name != "nt":
        def _limit_child() -> None:
            try:
                import resource

                resource.setrlimit(
                    resource.RLIMIT_AS,
                    (memory_limit_bytes, memory_limit_bytes),
                )
            except Exception:
                pass

        preexec_fn = _limit_child
    # Capture as bytes and decode manually with errors="replace": some tools
    # (svim/sniffles/cutesv on long-read FASTQs with non-ASCII headers) emit
    # non-UTF-8 bytes on stderr, which crashes text=True with UnicodeDecodeError
    # before the wrapper can inspect rc/stderr_tail.
    completed = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        capture_output=True,
        check=False,
        timeout=timeout,
        preexec_fn=preexec_fn,
    )
    stdout_text = completed.stdout.decode("utf-8", errors="replace") if completed.stdout else ""
    stderr_text = completed.stderr.decode("utf-8", errors="replace") if completed.stderr else ""
    if completed.returncode != 0:
        raise subprocess.CalledProcessError(
            completed.returncode, completed.args, stdout_text, stderr_text,
        )
    return subprocess.CompletedProcess(
        completed.args, completed.returncode, stdout_text, stderr_text,
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


def _detect_cgroup_memory_max_bytes() -> int | None:
    """Return this process's cgroup v2 memory.max in bytes, or None if unset.

    A SIGKILL with no stderr from the binary on a 754 GiB host is almost
    always a cgroup OOM kill — the user-slice on shared HPC login nodes is
    typically capped at 12 GiB. Surfacing this up front turns a confusing
    silent kill into an actionable error.
    """
    try:
        with open("/proc/self/cgroup", encoding="utf-8") as fh:
            line = fh.readline().strip()
    except OSError:
        return None
    # cgroup v2 lines look like "0::/user.slice/user-1234.slice/session-...scope"
    if "::" not in line:
        return None
    rel = line.split("::", 1)[1].lstrip("/")
    cur = Path("/sys/fs/cgroup")
    parts = [p for p in rel.split("/") if p]
    # Walk from the leaf upward; the first ancestor that has memory.max != "max"
    # is the binding limit.
    for depth in range(len(parts), -1, -1):
        candidate = cur.joinpath(*parts[:depth]) / "memory.max"
        try:
            value = candidate.read_text(encoding="utf-8").strip()
        except OSError:
            continue
        if value and value != "max":
            try:
                return int(value)
            except ValueError:
                continue
    return None


def _preflight_memory_check(cmd: list[str]) -> None:
    limit = _detect_cgroup_memory_max_bytes()
    if limit is None or limit <= 0:
        return
    gib = limit / (1024 ** 3)
    # 12 GiB is the typical user-slice cap on shared HPC nodes; the AMF
    # assembly stage routinely peaks at >20 GiB. Warn loudly so the operator
    # knows to switch to a compute node or raise the cgroup limit instead of
    # chasing a phantom binary bug.
    if gib < 24:
        sys.stderr.write(
            f"[mycosv preflight] WARNING: cgroup memory limit is {gib:.1f} GiB. "
            f"The MycoSV binary loads multiple references into RAM and may be "
            f"SIGKILLed by the cgroup OOM-killer. Recommend running on a node "
            f"with >=24 GiB available (slurm/srun --mem=32G), or raise the "
            f"user-slice memory.max.\n"
        )
        sys.stderr.write(f"[mycosv preflight] cmd: {' '.join(cmd)}\n")


def _mycosv_child_memory_limit_bytes() -> int | None:
    limit = _detect_cgroup_memory_max_bytes()
    if limit is None or limit <= 0:
        return None
    # Keep headroom for Python, tee/logging, shared libraries, and Slurm's own
    # accounting so a runaway binary fails inside the benchmark wrapper instead
    # of taking down the whole batch step via cgroup OOM.
    headroom = max(2 * 1024 ** 3, limit // 10)
    child_limit = limit - headroom
    if child_limit < 4 * 1024 ** 3:
        return None
    return child_limit


def _stream_subprocess_to_files(
    cmd: list[str],
    *,
    cwd: Path | None,
    timeout: int,
    stdout_path: Path,
    stderr_path: Path,
    memory_limit_bytes: int | None,
) -> subprocess.CompletedProcess[str]:
    """Run `cmd` with stdout/stderr piped DIRECTLY to disk files (not
    buffered in Python memory), so `tail -f stderr_path` shows live
    progress and a slurm time-out leaves a real log instead of an empty
    one.  Returns a CompletedProcess whose stdout/stderr fields contain the
    final on-disk tail (so existing callers that scan `result.stderr` for
    error patterns continue to work). On timeout, kills the process group
    and re-raises subprocess.TimeoutExpired with the same tail attached.
    """
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    stderr_path.parent.mkdir(parents=True, exist_ok=True)
    preexec_fn = None
    if memory_limit_bytes and memory_limit_bytes > 0 and os.name != "nt":
        def _limit_child() -> None:
            try:
                import resource
                resource.setrlimit(
                    resource.RLIMIT_AS,
                    (memory_limit_bytes, memory_limit_bytes),
                )
            except Exception:
                pass
            try:
                os.setsid()
            except Exception:
                pass
        preexec_fn = _limit_child
    elif os.name != "nt":
        def _new_session() -> None:
            try:
                os.setsid()
            except Exception:
                pass
        preexec_fn = _new_session

    def _tail(path: Path, n: int = 400) -> str:
        try:
            with path.open("rb") as fh:
                fh.seek(0, os.SEEK_END)
                size = fh.tell()
                read = min(size, 64 * 1024)
                fh.seek(size - read)
                data = fh.read(read).decode("utf-8", errors="replace")
        except OSError:
            return ""
        lines = data.splitlines()
        return "\n".join(lines[-n:])

    with stdout_path.open("wb") as out_fh, stderr_path.open("wb") as err_fh:
        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd) if cwd else None,
            stdout=out_fh,
            stderr=err_fh,
            preexec_fn=preexec_fn,
        )
        try:
            rc = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            # Try graceful SIGTERM to the process group, then SIGKILL.
            try:
                if os.name != "nt":
                    os.killpg(os.getpgid(proc.pid), 15)
                else:
                    proc.terminate()
            except Exception:
                pass
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                try:
                    if os.name != "nt":
                        os.killpg(os.getpgid(proc.pid), 9)
                    else:
                        proc.kill()
                except Exception:
                    pass
                proc.wait()
            tail_err = _tail(stderr_path)
            tail_out = _tail(stdout_path)
            sys.stderr.write(
                f"[mycosv] streamed subprocess timed out after {timeout}s — "
                f"killed via SIGKILL. Stderr tail at {stderr_path}:\n"
            )
            for line in tail_err.splitlines()[-40:]:
                sys.stderr.write(f"  {line}\n")
            raise subprocess.TimeoutExpired(
                cmd=cmd, timeout=timeout, output=tail_out, stderr=tail_err,
            )

    stdout_text = _tail(stdout_path, n=2000)
    stderr_text = _tail(stderr_path, n=2000)
    if rc != 0:
        raise subprocess.CalledProcessError(rc, cmd, stdout_text, stderr_text)
    return subprocess.CompletedProcess(cmd, rc, stdout_text, stderr_text)


def run_mycosv_command(
    cmd: list[str],
    cwd: Path | None = None,
    *,
    stream_stdout_path: Path | None = None,
    stream_stderr_path: Path | None = None,
) -> subprocess.CompletedProcess[str]:
    _preflight_memory_check(cmd)
    try:
        if stream_stdout_path is not None and stream_stderr_path is not None:
            # Stream mode: stdout/stderr go straight to disk so the operator
            # can `tail -f` the log and see whether MycoSV is alive or stuck
            # on a specific query. This is the fix for "the benchmark looks
            # hung but is actually working" — subprocess.run with capture_output
            # buffers all bytes until exit, which made any long per-query
            # routing/loading silently invisible.
            sys.stderr.write(
                f"[mycosv] streaming stderr -> {stream_stderr_path}\n"
                f"[mycosv]   `tail -f {stream_stderr_path}` to watch live\n"
            )
            return _stream_subprocess_to_files(
                cmd,
                cwd=cwd,
                timeout=_TOOL_TIMEOUT,
                stdout_path=stream_stdout_path,
                stderr_path=stream_stderr_path,
                memory_limit_bytes=_mycosv_child_memory_limit_bytes(),
            )
        return run(cmd, cwd=cwd, memory_limit_bytes=_mycosv_child_memory_limit_bytes())
    except subprocess.CalledProcessError as exc:
        # Surface SIGKILL as an OOM signal rather than a generic exit code.
        if exc.returncode in (-9, 137):
            limit = _detect_cgroup_memory_max_bytes()
            limit_str = (
                f"{limit / (1024 ** 3):.1f} GiB" if limit else "unbounded host RAM"
            )
            sys.stderr.write(
                f"[mycosv] binary was SIGKILLed (rc={exc.returncode}). This is "
                f"almost certainly a cgroup OOM kill (cgroup memory.max = "
                f"{limit_str}). Re-run with more memory or fewer refs.\n"
            )
        # The binary's stderr/stdout are otherwise lost because we capture
        # them into the CompletedProcess.  Replay the tail to the script's
        # stderr so the operator sees the actual failure mode (parse error,
        # missing index, segfault stack, etc.) instead of a bare exit code.
        for stream_name, stream in (("stderr", exc.stderr), ("stdout", exc.stdout)):
            if not stream:
                continue
            tail = stream.splitlines()[-40:]
            sys.stderr.write(f"[mycosv] binary {stream_name} (last {len(tail)} lines):\n")
            for line in tail:
                sys.stderr.write(f"  {line}\n")
        raise
    except subprocess.TimeoutExpired as exc:
        sys.stderr.write(
            f"[mycosv] binary timed out after {exc.timeout}s. Re-run with a "
            f"longer --tool-timeout, fewer queries, or smaller --tol-max-clade-genomes.\n"
        )
        raise


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


_FALLBACK_ENV_BINS: tuple[Path, ...] = (
    # Same default the install_tools.sh / Apptainer wrappers use. Looking here
    # makes the comparator pre-flight see installed tools even when the user
    # invoked python3 from a non-activated shell (e.g. via run_all_experiments.sh
    # which does not source conda.sh). Honors $CONDA_PREFIX / $MYCOSV_ENV_PATH.
    Path(os.environ.get("MYCOSV_ENV_PATH", os.environ.get("CONDA_PREFIX", "/dev/null"))) / "bin",
    Path("/mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/tools/envs/envs/fungi_graph_sv/bin"),
)


def _prepend_env_bins_to_path() -> None:
    """Prepend the project's known conda env bin to PATH (idempotent).

    Comparator subprocess calls (`minimap2 ...`, `syri ...`, `delly ...`) use
    bare tool names, so they require the binaries to be on PATH at exec time.
    When run_all_experiments.sh launches python3 without activating conda, PATH
    only contains system dirs and the comparators silently fail with FileNotFoundError.
    Prepending here lets the operator skip `conda activate` and still get a full
    comparator run.
    """
    current = os.environ.get("PATH", "").split(os.pathsep)
    current_set = set(current)
    additions: list[str] = []
    for env_bin in _FALLBACK_ENV_BINS:
        if not env_bin.is_dir():
            continue
        s = str(env_bin)
        if s and s not in current_set:
            additions.append(s)
            current_set.add(s)
    if additions:
        os.environ["PATH"] = os.pathsep.join(additions + current)


_prepend_env_bins_to_path()


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


_HTTP_TIMEOUT = 300  # seconds; prevents hanging on slow/unresponsive NCBI

# NCBI / EBI assembly mirrors return 503 in bursts when the server pool is
# saturated, and 429 when we cross their rate limit.  These statuses are
# transient: a short wait + retry recovers the vast majority of failures
# the bulk download loop sees.  4xx other than 429 are permanent (e.g. 404
# for an assembly that doesn't expose a GFF) and are not retried.
_HTTP_RETRY_STATUS = {429, 500, 502, 503, 504}
# 8 attempts at 2,4,8,16,32,60,60,60 = ~242 s end-to-end. The previous 5-attempt
# schedule (max 30s total) gave up while NCBI's ftp.ncbi.nlm.nih.gov mirror was
# still in a 503 burst — saw 4 hard failures (GCA_000507425_3, GCA_003184365_1,
# GCA_020503465_1, GCA_051529335_1) in the prep run. Going to 8 attempts plus a
# longer backoff cap keeps the failure rate near zero on a typical server hiccup
# without blowing wall time on a permanent outage.
_HTTP_MAX_ATTEMPTS = 8
_HTTP_BACKOFF_BASE = 2.0
_HTTP_BACKOFF_CAP  = 90.0
# Max workers for ftp.ncbi.nlm.nih.gov. 8 was the prior default but the
# 2026-05-15 prep run (slurm-14936460) still saw a sustained 503 storm —
# hundreds of [retry] lines, multiple URLs climbing to attempt 4/8. Dropping
# to 6 + the per-host semaphore/cooldown below keeps GFF throughput close to
# the 8-worker run while staying inside NCBI's well-behaved client window.
_HTTP_MAX_PARALLEL_FTP_NCBI = 6

# Hosts we treat as a single overloadable resource. When a 503/429 fires on
# any of these, every other in-flight request to the same host pauses for a
# short cooldown so workers stop hammering the same throttled pool in
# lockstep. ftp.ncbi.nlm.nih.gov is the only one that bursts in practice; the
# tuple shape keeps it cheap to add ENA/Ensembl if they ever start throttling.
_THROTTLED_HOSTS: tuple[str, ...] = ("ftp.ncbi.nlm.nih.gov",)
_NCBI_HOST_SEM = threading.BoundedSemaphore(_HTTP_MAX_PARALLEL_FTP_NCBI)
_NCBI_COOLDOWN_LOCK = threading.Lock()
_NCBI_COOLDOWN_UNTIL = 0.0  # monotonic deadline; 0 means no cooldown active
# How long a 503/429 blocks every NCBI request. Floor at 3 s (NCBI's
# throttle window is short) and let server-supplied Retry-After lengthen it
# up to the backoff cap.
_NCBI_COOLDOWN_MIN_SECONDS = 3.0


def _request_throttle_host(req: urllib.request.Request) -> str | None:
    """Return the throttled host for ``req``, or None if no throttle applies."""
    try:
        host = (req.host or "").lower()
    except AttributeError:
        host = ""
    if not host:
        return None
    for throttled in _THROTTLED_HOSTS:
        if host == throttled or host.endswith("." + throttled):
            return throttled
    return None


def _wait_for_ncbi_cooldown() -> None:
    deadline = _NCBI_COOLDOWN_UNTIL
    if deadline <= 0.0:
        return
    remaining = deadline - time.monotonic()
    if remaining > 0:
        time.sleep(remaining)


def _trigger_ncbi_cooldown(seconds: float) -> None:
    if seconds <= 0:
        return
    global _NCBI_COOLDOWN_UNTIL
    with _NCBI_COOLDOWN_LOCK:
        new_deadline = time.monotonic() + seconds
        if new_deadline > _NCBI_COOLDOWN_UNTIL:
            _NCBI_COOLDOWN_UNTIL = new_deadline


def _http_retry_sleep(attempt: int, retry_after: str | None) -> float:
    # Honor server-supplied Retry-After (seconds form; HTTP-date form is rare
    # for these mirrors and we'd just fall back to the exponential schedule).
    if retry_after:
        try:
            wait = float(retry_after.strip())
            if wait > 0:
                return min(wait, _HTTP_BACKOFF_CAP)
        except ValueError:
            pass
    # Add small jitter so a thread-pool of N workers does not retry in lockstep
    # and re-trigger the same 503 burst.
    base = min(_HTTP_BACKOFF_BASE * (2 ** attempt), _HTTP_BACKOFF_CAP)
    jitter = random.uniform(0.0, min(base * 0.25, 5.0))
    return base + jitter


def _http_open_with_retry(req: urllib.request.Request):
    """Open `req` with retry-on-transient-error.  Returns the urlopen response;
    caller is responsible for closing it (use as a context manager)."""
    throttled_host = _request_throttle_host(req)
    last_exc: Exception | None = None
    for attempt in range(_HTTP_MAX_ATTEMPTS):
        # Honor the shared cooldown before every attempt, including the first:
        # a worker that arrives after another worker already tripped 503 should
        # not immediately retry into the same overload window.
        if throttled_host is not None:
            _wait_for_ncbi_cooldown()
            _NCBI_HOST_SEM.acquire()
        # pending_sleep is the backoff we should sleep AFTER releasing the
        # semaphore — keeping the sleep outside the host-limited region means
        # the semaphore caps in-flight requests, not idle-waiting workers.
        pending_sleep: float = 0.0
        try:
            try:
                return urllib.request.urlopen(req, timeout=_HTTP_TIMEOUT)
            except urllib.error.HTTPError as exc:
                if exc.code not in _HTTP_RETRY_STATUS or attempt == _HTTP_MAX_ATTEMPTS - 1:
                    raise
                retry_after = exc.headers.get("Retry-After") if exc.headers else None
                wait = _http_retry_sleep(attempt, retry_after)
                # 503/429 means NCBI is throttling our IP — pause every other
                # worker for the longer of the per-request backoff and a 3 s
                # floor so the server-side window can clear.
                if throttled_host is not None and exc.code in {429, 503}:
                    _trigger_ncbi_cooldown(max(wait, _NCBI_COOLDOWN_MIN_SECONDS))
                sys.stderr.write(
                    f"[retry] HTTP {exc.code} for {req.full_url}; sleeping {wait:.1f}s "
                    f"(attempt {attempt + 1}/{_HTTP_MAX_ATTEMPTS})\n"
                )
                sys.stderr.flush()
                pending_sleep = wait
                last_exc = exc
            except (urllib.error.URLError, TimeoutError) as exc:
                if attempt == _HTTP_MAX_ATTEMPTS - 1:
                    raise
                wait = _http_retry_sleep(attempt, None)
                sys.stderr.write(
                    f"[retry] network error for {req.full_url}: {exc}; "
                    f"sleeping {wait:.1f}s (attempt {attempt + 1}/{_HTTP_MAX_ATTEMPTS})\n"
                )
                sys.stderr.flush()
                pending_sleep = wait
                last_exc = exc
        finally:
            if throttled_host is not None:
                _NCBI_HOST_SEM.release()
        if pending_sleep > 0:
            time.sleep(pending_sleep)
    # Unreachable: the final attempt either returns or raises above.
    assert last_exc is not None
    raise last_exc


def http_get_text(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "MycoSV-real-benchmark/1.0"})
    with _http_open_with_retry(req) as resp:
        return resp.read().decode("utf-8")


def http_get_text_cached(url: str, cache_path: Path) -> str:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    if cache_path.exists() and cache_path.stat().st_size > 0:
        return cache_path.read_text(encoding="utf-8")
    text = http_get_text(url)
    # See http_download for the rationale on the per-process tmp suffix.
    tmp = cache_path.with_suffix(
        cache_path.suffix + f".part.{os.getpid()}.{uuid.uuid4().hex[:8]}"
    )
    try:
        tmp.write_text(text, encoding="utf-8")
        os.replace(tmp, cache_path)
    except Exception:
        if tmp.exists():
            tmp.unlink()
        raise
    return text


def http_download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    req = urllib.request.Request(url, headers={"User-Agent": "MycoSV-real-benchmark/1.0"})
    # Per-process tmp suffix: when two concurrent runs share a data_cache and
    # both download the same Ensembl GFF (one Ensembl record can back many NCBI
    # accessions), a static `.part` collides. The first rename wins, the second
    # then sees its `.part` gone and crashes with FileNotFoundError mid-rename,
    # killing the whole prepare step. Salting with PID+rand makes each writer's
    # tmp file private; the cache-hit guard above keeps races from re-downloading.
    tmp = dest.with_suffix(
        dest.suffix + f".part.{os.getpid()}.{uuid.uuid4().hex[:8]}"
    )
    try:
        with _http_open_with_retry(req) as resp, tmp.open("wb") as out:
            content_length = resp.getheader("Content-Length")
            shutil.copyfileobj(resp, out)
        # Validate Content-Length when the server advertised one. ENA / NCBI
        # mirrors occasionally close the connection after writing a partial
        # stream without raising IncompleteRead, which previously cached a
        # truncated .gz that later blew up at decompress time as
        # "Compressed file ended before the end-of-stream marker was reached".
        if content_length is not None:
            try:
                expected = int(content_length)
            except ValueError:
                expected = None
            if expected is not None:
                got = tmp.stat().st_size
                if got != expected:
                    tmp.unlink(missing_ok=True)
                    raise IOError(
                        f"Truncated download for {url}: got {got} bytes, expected {expected}"
                    )
        # os.replace is atomic and overwrites if a concurrent writer already
        # finished and published the same dest — both writers end up with a
        # valid file at dest, and neither raises.
        os.replace(tmp, dest)
    except Exception:
        if tmp.exists():
            tmp.unlink()
        raise
    return dest


def data_cache_base(args: argparse.Namespace, out_dir: Path) -> Path:
    """Return the persistent download cache for prepare-style commands."""
    configured = getattr(args, "data_cache_dir", None)
    return configured.resolve() if configured else (out_dir / "downloads")


def cached_filename_for_source(source: str, fallback_name: str) -> str:
    """Stable cache filename for arbitrary URLs and local manifest paths.

    NCBI assembly downloads already carry accession-specific basenames. Custom
    public manifests are less predictable, so URL-backed inputs get a short hash
    prefix to avoid collisions such as multiple providers exposing genome.fa.gz.
    """
    parsed = urllib.parse.urlparse(source)
    name = Path(parsed.path).name if parsed.path else ""
    if not name:
        name = fallback_name
    if looks_like_url(source):
        digest = hashlib.sha1(source.encode("utf-8")).hexdigest()[:12]
        return f"{digest}_{name}"
    return name


def maybe_gunzip(path: Path, keep_gz: bool = True) -> Path:
    # Keep .gz in cache; the C++ binary reads gzip natively via popen.
    # Never write a decompressed copy alongside the archive.
    return path


def open_text_auto(path: Path):
    # Use errors="replace" so malformed bytes inside ENA-fetched FASTQs (mixed
    # encodings, partial gzip streams, or rare non-ASCII headers) don't kill
    # the bounded-FASTQ writer or the comparator pipeline. Replacement chars
    # only appear in headers/qual lines that downstream tools tolerate.
    try:
        with path.open("rb") as fh:
            magic = fh.read(2)
    except OSError:
        magic = b""
    if path.suffix == ".gz" or magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _fastq_record_iter(path: Path):
    with open_text_auto(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not qual:
                return
            yield header, seq, plus, qual


def subset_fastq_records(src: Path, dest: Path, max_records: int) -> tuple[Path, int, bool]:
    """Write at most max_records FASTQ records to dest.

    Public ENA runs can be many gigabytes. MycoSV already caps reads internally,
    but external comparators map the FASTQ path directly; feeding them a bounded
    subset keeps benchmark wall time predictable without changing the prepared
    manifest on disk.
    """
    if max_records <= 0:
        return src, 0, False
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest, max_records, True
    count = 0
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        with tmp.open("w", encoding="utf-8") as out_fh:
            for header, seq, plus, qual in _fastq_record_iter(src):
                if count >= max_records:
                    break
                out_fh.write(header)
                out_fh.write(seq)
                out_fh.write(plus)
                out_fh.write(qual)
                count += 1
        tmp.rename(dest)
    except Exception:
        if tmp.exists():
            tmp.unlink()
        raise
    return dest, count, count > 0


def cap_read_query_inputs(
    query_manifest: list[dict[str, str]],
    out_dir: Path,
    mode: str,
    max_short_reads: int,
    max_long_reads: int,
    *,
    mycosv_use_full_reads: bool = False,
) -> list[dict[str, str]]:
    if mode not in {"auto", "short-reads", "long-reads"}:
        return query_manifest

    capped_rows: list[dict[str, str]] = []
    subset_dir = out_dir / "read_subsets"
    for row in query_manifest:
        row_mode = mode if mode in {"short-reads", "long-reads"} else (row.get("query_mode") or "assembly")
        if row_mode not in {"short-reads", "long-reads"}:
            capped_rows.append(row)
            continue
        max_records = max_long_reads if row_mode == "long-reads" else max_short_reads
        if max_records <= 0:
            capped_rows.append(row)
            continue
        original = locate_query_path(row)
        if sequence_kind_from_name(original.name) != "fastq":
            capped_rows.append(row)
            continue
        # Keep the bounded copy name stable and sample-like. Both MycoSV and
        # the external comparators use this bounded FASTQ by default; otherwise
        # million-real short-read runs can hand MycoSV multi-GB public FASTQs
        # while comparators see tiny subsets, causing hour-scale stalls,
        # cgroup kills, and empty calls.vcf files.
        dest = subset_dir / f"{normalize_name(row.get('query_asm', original.stem))}.fastq"
        try:
            capped_path, kept, capped = subset_fastq_records(original, dest, max_records)
        except Exception as exc:
            sys.stderr.write(
                f"[reads-mode] could not create bounded FASTQ for "
                f"{row.get('query_asm', original.name)}: {exc}; using original\n"
            )
            capped_rows.append(row)
            continue
        new_row = dict(row)
        new_row["path"] = str(capped_path)
        if mycosv_use_full_reads:
            new_row["mycosv_path"] = str(original)
        capped_rows.append(new_row)
        if capped:
            sys.stderr.write(
                f"[reads-mode] benchmark input capped for "
                f"{row.get('query_asm', original.name)}: {kept} reads -> {capped_path}\n"
            )
    return capped_rows


def filter_assembly_query_inputs(
    query_manifest: list[dict[str, str]],
    out_dir: Path,
    max_contigs: int,
    max_bp: int,
) -> list[dict[str, str]]:
    """Drop assembly queries that are too fragmented/large for matrix runs.

    Some public fungal assemblies are MAG-style or highly fragmented drafts.
    They are valid biological inputs, but they can make the hierarchical
    assembly caller spend hours before producing a first VCF row. Filtering
    them here keeps the panel matrix moving while recording exactly what was
    skipped so the full query can be rerun separately on a long-walltime node.
    """
    if max_contigs <= 0 and max_bp <= 0:
        return query_manifest

    kept: list[dict[str, str]] = []
    skipped: list[dict[str, str]] = []
    for row in query_manifest:
        query_path = locate_query_path(row)
        n_contigs, total_bp, longest_bp = _fasta_stats(query_path)
        reasons: list[str] = []
        if max_contigs > 0 and n_contigs > max_contigs:
            reasons.append(f"contigs>{max_contigs}")
        if max_bp > 0 and total_bp > max_bp:
            reasons.append(f"bp>{max_bp}")
        if reasons:
            skipped.append({
                "query_asm": row.get("query_asm", query_path.name),
                "path": str(query_path),
                "contigs": str(n_contigs),
                "total_bp": str(total_bp),
                "longest_bp": str(longest_bp),
                "reason": ",".join(reasons),
            })
            continue
        kept.append(row)

    if skipped:
        skipped_path = out_dir / "SKIPPED_ASSEMBLY_QUERIES.tsv"
        write_tsv(
            skipped_path,
            skipped,
            ["query_asm", "path", "contigs", "total_bp", "longest_bp", "reason"],
        )
        sys.stderr.write(
            f"[assembly-filter] skipped {len(skipped)} oversized/fragmented "
            f"assembly query(s); kept {len(kept)}. Details: {skipped_path}\n"
        )
    return kept


_FASTA_CONTIG_CACHE: dict[str, frozenset[str]] = {}


def fasta_contig_names(fasta_path: Path) -> frozenset[str]:
    """Return the set of sequence IDs in a (possibly gzipped) FASTA.

    Used to scope MycoSV's reference-coordinate calls to a specific
    benchmark reference per query, so precision/recall against pairwise
    comparators (which only see one reference) is measured fairly.
    """
    key = str(fasta_path)
    cached = _FASTA_CONTIG_CACHE.get(key)
    if cached is not None:
        return cached
    contigs: set[str] = set()
    if not fasta_path.exists():
        _FASTA_CONTIG_CACHE[key] = frozenset()
        return frozenset()
    try:
        with open_text_auto(fasta_path) as fh:
            for line in fh:
                if line.startswith(">"):
                    contigs.add(line[1:].split()[0] if line[1:].strip() else "")
    except OSError:
        _FASTA_CONTIG_CACHE[key] = frozenset()
        return frozenset()
    contigs.discard("")
    frozen = frozenset(contigs)
    _FASTA_CONTIG_CACHE[key] = frozen
    return frozen


INPUT_PREFLIGHT_FIELDS = [
    "role", "query_asm", "path", "expected_format", "status", "reason",
]


def _preflight_sequence_file(path: Path, expected_format: str) -> tuple[bool, str]:
    """Cheaply verify that a FASTA/FASTQ path exists, is readable, and starts
    like the format the benchmark will hand to MycoSV/comparators.

    This intentionally reads only the first record. It catches corrupt cached
    gzip payloads and manifest/path mixups before the expensive MycoSV and
    comparator launches allocate large reference/query indexes.
    """
    if not path.exists():
        return False, "missing"
    try:
        if path.stat().st_size == 0:
            return False, "empty"
    except OSError as exc:
        return False, f"stat_failed:{type(exc).__name__}"

    try:
        with open_text_auto(path) as fh:
            if expected_format in {"fastq", "fastq_or_fasta"}:
                header = fh.readline()
                if expected_format == "fastq_or_fasta" and header.startswith(">"):
                    saw_sequence = False
                    for i, line in enumerate(fh):
                        if i > 10000:
                            break
                        if line.startswith(">"):
                            if saw_sequence:
                                break
                            continue
                        if line.strip():
                            saw_sequence = True
                            break
                    return (True, "ok") if saw_sequence else (False, "fasta_sequence_missing")
                seq = fh.readline()
                plus = fh.readline()
                qual = fh.readline()
                if not header:
                    return False, "no_fastq_record"
                if not header.startswith("@"):
                    return False, "fastq_header_not_at"
                if not seq.strip():
                    return False, "fastq_empty_sequence"
                if not plus.startswith("+"):
                    return False, "fastq_plus_line_missing"
                if not qual:
                    return False, "fastq_quality_missing"
                return True, "ok"

            # FASTA: allow comments/blank preamble, but require one header and
            # at least one sequence line so gzip errors surface immediately.
            saw_header = False
            saw_sequence = False
            for i, line in enumerate(fh):
                if i > 10000:
                    break
                stripped = line.strip()
                if not stripped or stripped.startswith(";"):
                    continue
                if stripped.startswith(">"):
                    saw_header = True
                    continue
                if saw_header:
                    saw_sequence = True
                    break
            if not saw_header:
                return False, "fasta_header_missing"
            if not saw_sequence:
                return False, "fasta_sequence_missing"
            return True, "ok"
    except (OSError, EOFError, gzip.BadGzipFile, UnicodeError) as exc:
        return False, f"read_failed:{type(exc).__name__}"


def preflight_benchmark_inputs(
    query_manifest: list[dict[str, str]],
    out_dir: Path,
    mode: str,
) -> list[dict[str, str]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    seen_refs: set[str] = set()
    query_expected = "fastq_or_fasta" if mode in {"auto", "short-reads", "long-reads"} else "fasta"

    def add(role: str, query_asm: str, raw_path: str, expected_format: str) -> None:
        path = Path(raw_path)
        ok, reason = _preflight_sequence_file(path, expected_format)
        rows.append({
            "role": role,
            "query_asm": query_asm,
            "path": str(path),
            "expected_format": expected_format,
            "status": "ok" if ok else "fail",
            "reason": reason,
        })

    for row in query_manifest:
        query_asm = row.get("query_asm", ".")
        query_path = (row.get("mycosv_path") or row.get("path") or "").strip()
        if query_path:
            add("query", query_asm, query_path, query_expected)
        else:
            rows.append({
                "role": "query",
                "query_asm": query_asm,
                "path": "",
                "expected_format": query_expected,
                "status": "fail",
                "reason": "missing_manifest_path",
            })
        bench_ref = (row.get("benchmark_ref_fasta") or "").strip()
        if bench_ref and bench_ref != "." and bench_ref not in seen_refs:
            seen_refs.add(bench_ref)
            add("benchmark_ref", query_asm, bench_ref, "fasta")

    write_tsv(out_dir / "INPUT_PREFLIGHT.tsv", rows, INPUT_PREFLIGHT_FIELDS)
    failed = [r for r in rows if r["status"] != "ok"]
    if failed:
        marker = out_dir / "INPUT_PREFLIGHT_FAILED.txt"
        marker.write_text(
            "Benchmark input preflight failed before launching expensive "
            "MycoSV/comparator work. See INPUT_PREFLIGHT.tsv for paths and "
            "reasons. Remove corrupt cached files or fix query_manifest.tsv, "
            "then rerun using --reuse-index-dir when applicable.\n",
            encoding="utf-8",
        )
        examples = "; ".join(
            f"{r['role']}:{r['query_asm']}:{r['reason']}"
            for r in failed[:5]
        )
        raise ValueError(
            f"Benchmark input preflight failed for {len(failed)} file(s): {examples}. "
            f"Details: {out_dir / 'INPUT_PREFLIGHT.tsv'}"
        )
    return rows


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


def ncbi_source_components(source: str) -> list[str]:
    if source == NCBI_BEST_SOURCE:
        return ["ncbi-refseq", "ncbi-genbank"]
    if source not in NCBI_ASSEMBLY_SUMMARY:
        raise ValueError(f"Unknown NCBI source {source!r}")
    return [source]


def assembly_summary_cache_path(cache_base: Path, source: str) -> Path:
    return cache_base / "assembly_summaries" / f"{source}_fungi_assembly_summary.txt"


def source_priority(row: dict[str, str]) -> int:
    source = row.get("_catalog_source") or row.get("source_catalog") or ""
    if source == "ncbi-refseq":
        return 2
    if source == "ncbi-genbank":
        return 1
    return 0


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


def row_quality_key(row: dict[str, str]) -> tuple[int, int, int, int, int, str]:
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
        latest,
        assembly_level_rank(row.get("assembly_level", "")),
        complete,
        refscore,
        source_priority(row),
        release_date,
    )


def assembly_pair_key(row: dict[str, str]) -> str:
    accession = (row.get("assembly_accession") or "").strip()
    paired = (row.get("gbrs_paired_asm") or "").strip()
    if paired and paired.lower() not in {"na", ".", "none"}:
        return "|".join(sorted([accession, paired]))
    return accession


def deduplicate_best_assembly_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    best_by_pair: dict[str, dict[str, str]] = {}
    for row in rows:
        key = assembly_pair_key(row)
        if not key:
            continue
        current = best_by_pair.get(key)
        if current is None or row_quality_key(row) > row_quality_key(current):
            best_by_pair[key] = row
    return list(best_by_pair.values())


def fetch_ncbi_assembly_rows(
    source: str,
    cache_base: Path,
    *,
    progress: bool = False,
) -> tuple[list[dict[str, str]], list[Path]]:
    rows: list[dict[str, str]] = []
    cache_paths: list[Path] = []
    for component in ncbi_source_components(source):
        summary_url = NCBI_ASSEMBLY_SUMMARY[component]
        cache_path = assembly_summary_cache_path(cache_base, component)
        cache_paths.append(cache_path)
        if progress:
            print(f"[1/4] Fetching NCBI assembly summary: {summary_url}", flush=True)
        summary_text = http_get_text_cached(summary_url, cache_path)
        parsed = parse_assembly_summary(summary_text)
        for row in parsed:
            row["_catalog_source"] = component
        rows.extend(parsed)
        if progress:
            print(f"      parsed {len(parsed)} rows from {component}", flush=True)
    if len(ncbi_source_components(source)) > 1:
        before = len(rows)
        rows = deduplicate_best_assembly_rows(rows)
        if progress:
            print(
                f"      retained {len(rows)} best assemblies after RefSeq/GenBank pairing "
                f"deduplication ({before} raw rows)",
                flush=True,
            )
    return rows, cache_paths


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


# NCBI Taxonomy renames that have left panel presets out of sync. Each entry
# maps a panel-preset species name to additional accepted names that NCBI's
# assembly_summary may use today. Add bidirectionally so legacy and current
# names both match. Sourced from the 2019–2023 Saccharomycetes reclassification
# (Takashima & Sugita) and the 2024 ICTF list updates.
SPECIES_ALIASES: dict[str, list[str]] = {
    "candida glabrata": ["nakaseomyces glabratus", "nakaseomyces glabrata", "[candida] glabrata"],
    "candida krusei": ["pichia kudriavzevii", "issatchenkia orientalis"],
    "candida lusitaniae": ["clavispora lusitaniae"],
    "candida guilliermondii": ["meyerozyma guilliermondii"],
    "candida tropicalis": ["[candida] tropicalis"],
    # Cryptococcus species complex split — keep both names recognised.
    "cryptococcus neoformans": ["cryptococcus neoformans var. grubii", "cryptococcus deneoformans"],
    # Leptosphaeria maculans — split into species complex (LepmaJN3, LepmaPHW1
    # etc.). NCBI assemblies still use `Leptosphaeria maculans` in
    # organism_name but several MAGs / re-annotations use `Plenodomus lingam`.
    # Adding both keeps the panel selection robust against the 2017 reclass.
    "leptosphaeria maculans": [
        "leptosphaeria maculans 'brassicae'",
        "leptosphaeria maculans 'lepidii'",
        "plenodomus lingam",
    ],
    # Ustilago maydis — modern teleomorph name is Mycosarcoma maydis (2018
    # ICTF revision). NCBI assemblies still mostly carry `Ustilago maydis`
    # but newer Mexican/U.S. plant-pathology submissions use Mycosarcoma.
    "ustilago maydis": [
        "mycosarcoma maydis",
        "[ustilago] maydis",
    ],
}


def match_species(row: dict[str, str], species: str) -> bool:
    organism = (row.get("organism_name") or "").lower()
    target = species.strip().lower()
    candidates = [target] + SPECIES_ALIASES.get(target, [])
    for cand in candidates:
        if organism.startswith(cand) or f" {cand} " in f" {organism} ":
            return True
    return False


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
    selected.sort(key=lambda r: (row_quality_key(r), species_group_key(r)), reverse=True)
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


# NCBI BioSample attribute names that carry ecologically relevant phenotypic info.
_PHENOTYPE_ATTRS = {
    "isolation_source", "host", "disease", "geographic_location",
    "collection_date", "env_biome", "env_feature", "env_material",
    "pathogenicity", "trophic_level", "lifestyle", "tissue",
    "culture_collection", "strain", "substrain",
}


_VALID_BIOSAMPLE_RE = re.compile(r"^SAM[NED][A-Z]?\d+$")


def _is_valid_biosample_id(bid: str) -> bool:
    # NCBI's assembly_summary populates missing fields with the literal "na".
    # Sending those (or other placeholders) to efetch returns HTTP 400 and
    # spams the operator log. Real BioSample accessions match SAMN/SAMEA/SAMD.
    return bool(bid) and _VALID_BIOSAMPLE_RE.match(bid.strip()) is not None


def fetch_ncbi_biosample_phenotypes(
    biosample_ids: list[str],
    cache_path: Path | None = None,
) -> dict[str, dict[str, str]]:
    """Download BioSample phenotypic attributes for a list of BioSample accessions.

    Results are stored in *cache_path* (JSON) so subsequent runs skip the
    network round-trip.  Pass cache_path=data_cache_dir/'phenotypic_metadata.json'
    to persist across experiments.
    """
    cache: dict[str, dict[str, str]] = {}
    if cache_path and cache_path.exists():
        try:
            cache = json.loads(cache_path.read_text(encoding="utf-8"))
        except Exception:
            cache = {}

    wanted = [bid for bid in biosample_ids if _is_valid_biosample_id(bid) and bid not in cache]
    skipped = [bid for bid in biosample_ids if bid and not _is_valid_biosample_id(bid)]
    if skipped:
        sample = ", ".join(sorted(set(skipped))[:5])
        sys.stderr.write(
            f"[phenotype] skipping {len(skipped)} non-BioSample id(s) "
            f"(e.g. {sample}); these are NCBI 'na' placeholders or non-conforming.\n"
        )
    if not wanted:
        return cache

    # Batch size: NCBI eutils tolerates up to ~500 IDs via GET, but BioSample
    # accessions in this pipeline are sometimes mixed with sample-set IDs that
    # individually expand to multi-record XML; large GET batches can hit the
    # endpoint's per-request size cap and return HTTP 400. 100/batch is the
    # documented safe value for efetch and stays well under URL limits.
    for start in range(0, len(wanted), 100):
        batch = wanted[start : start + 100]
        query = urllib.parse.urlencode({
            "db": "biosample",
            "id": ",".join(batch),
            "retmode": "xml",
            # NCBI usage policy: identify the client. Without these, eutils
            # silently throttles and occasionally responds with 400 instead
            # of 429 when the upstream rejects the batch.
            "tool": "mycosv-benchmark",
            "email": "shilpa.garg2k7@gmail.com",
        })
        try:
            xml_text = http_get_text(
                f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{query}"
            )
        except Exception as exc:
            sys.stderr.write(
                f"[warn] BioSample phenotype fetch failed for batch "
                f"({len(batch)} ids, first={batch[0]!r}): "
                f"{type(exc).__name__}: {exc}\n"
            )
            continue
        try:
            root = ET.fromstring(xml_text)
        except Exception:
            continue
        for bs in root.findall(".//BioSample"):
            accession = bs.get("accession", "")
            if not accession:
                for attr in bs.findall(".//Id[@db='BioSample']"):
                    accession = (attr.text or "").strip()
                    break
            if not accession:
                continue
            attrs: dict[str, str] = {}
            for attr_el in bs.findall(".//Attribute"):
                name = (attr_el.get("attribute_name") or attr_el.get("harmonized_name") or "").lower()
                val = (attr_el.text or "").strip()
                if name in _PHENOTYPE_ATTRS and val:
                    attrs[name] = val
            cache[accession] = attrs
        time.sleep(0.34)

    if cache_path:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(cache, indent=2, sort_keys=True), encoding="utf-8")
    return cache


def ncbi_download_targets(row: dict[str, str], include_gff: bool) -> list[tuple[str, str]]:
    ftp_path = row["ftp_path"].rstrip("/")
    stem = ftp_path.split("/")[-1]
    targets = [(f"{ftp_path}/{stem}_genomic.fna.gz", f"{stem}_genomic.fna.gz")]
    if include_gff:
        targets.append((f"{ftp_path}/{stem}_genomic.gff.gz", f"{stem}_genomic.gff.gz"))
    return targets


def ncbi_genbank_target(row: dict[str, str]) -> tuple[str, str]:
    ftp_path = row["ftp_path"].rstrip("/")
    stem = ftp_path.split("/")[-1]
    return f"{ftp_path}/{stem}_genomic.gbff.gz", f"{stem}_genomic.gbff.gz"


# Ensembl Fungi FTP layout: /pub/fungi/<release>/gff3/<species_dir>/<File>.gff3.gz
# `current` symlinks to the latest release. The release number is also baked
# into the GFF filename (e.g. ...R64-1-1.62.gff3.gz for release 62), so we
# discover it once from the directory listing of pub/fungi/ and cache it.
#
# Per-species index: species_EnsemblFungi.txt (TSV, ~300 KB) lists every
# species directory along with its NCBI assembly_accession and taxonomy_id —
# small enough to fetch once and re-use across panels. We deliberately do not
# pull species_metadata_EnsemblFungi.json (376 MB) since it is overkill for
# accession-keyed lookup.
_ENSEMBL_FUNGI_FTP_BASE = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/current"
_ENSEMBL_FUNGI_SPECIES_TSV_URL = (
    f"{_ENSEMBL_FUNGI_FTP_BASE}/species_EnsemblFungi.txt"
)
_ENSEMBL_FUNGI_PARENT_URL = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/"


def _ensembl_fungi_release(cache_dir: Path) -> str | None:
    """Discover the current Ensembl Fungi release number (e.g. "62").

    Cached under data_cache to avoid an HTTP roundtrip per prepare. Returns
    None when the FTP listing is unreachable; the caller skips the Ensembl
    fallback in that case.
    """
    cache_path = cache_dir / "ensembl_fungi_release.txt"
    if cache_path.exists() and cache_path.stat().st_size > 0:
        value = cache_path.read_text(encoding="utf-8").strip()
        if value.isdigit():
            return value
    try:
        listing = http_get_text(_ENSEMBL_FUNGI_PARENT_URL)
    except Exception as exc:  # noqa: BLE001 - best-effort
        sys.stderr.write(f"[ensembl-fungi] release discovery failed: {exc}\n")
        return None
    releases = sorted(
        {int(m) for m in re.findall(r'href="release-(\d+)/"', listing)},
        reverse=True,
    )
    if not releases:
        return None
    release = str(releases[0])
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(release + "\n", encoding="utf-8")
    return release


def _ensembl_fungi_species_index(cache_dir: Path) -> dict[str, dict[str, str]]:
    """Return Ensembl Fungi species TSV indexed by NCBI accession + taxid.

    Empty dict on any fetch / parse failure: the Ensembl fallback is best-effort
    and must never block the primary NCBI GFF + GBFF path.
    """
    cache_path = cache_dir / "ensembl_fungi_species.tsv"
    cache_dir.mkdir(parents=True, exist_ok=True)
    if not cache_path.exists() or cache_path.stat().st_size == 0:
        try:
            http_download(_ENSEMBL_FUNGI_SPECIES_TSV_URL, cache_path)
        except Exception as exc:  # noqa: BLE001 - best-effort
            sys.stderr.write(
                f"[ensembl-fungi] species TSV fetch failed: {exc}\n"
            )
            return {}
    index: dict[str, dict[str, str]] = {}
    try:
        with cache_path.open(encoding="utf-8") as fh:
            # Header line begins with '#'; csv.DictReader can read it after
            # stripping the leading '#'.
            first = fh.readline().lstrip("#").rstrip("\n")
            field_names = first.split("\t")
            reader = csv.DictReader(fh, fieldnames=field_names, delimiter="\t")
            for rec in reader:
                if not isinstance(rec, dict):
                    continue
                # Keys we want to look up by: NCBI assembly_accession (with
                # and without version), NCBI taxonomy_id.
                acc = (rec.get("assembly_accession") or "").strip()
                taxid = (rec.get("taxonomy_id") or "").strip()
                if acc:
                    index.setdefault(acc, rec)
                    base = acc.split(".", 1)[0]
                    if base and base != acc:
                        index.setdefault(base, rec)
                if taxid:
                    index.setdefault(taxid, rec)
    except (OSError, csv.Error) as exc:
        sys.stderr.write(f"[ensembl-fungi] species TSV parse failed: {exc}\n")
        return {}
    return index


_ENSEMBL_FUNGI_COLLECTION_RE = re.compile(r"^(fungi_[a-z0-9]+_collection)_core_")


def _ensembl_fungi_collection(rec: dict[str, str]) -> str | None:
    """Extract Ensembl's collection bucket from species_EnsemblFungi.txt's
    core_db field. Top-level model species (S. cerevisiae, N. crassa, ...)
    have no collection; non-model species are grouped under
    fungi_<phylum><N>_collection (e.g. fungi_mucoromycota1_collection).
    Without this hint, the GFF URL silently 404s for ~80 % of fungal asms.
    """
    core_db = (rec.get("core_db") or "").strip()
    if not core_db:
        return None
    m = _ENSEMBL_FUNGI_COLLECTION_RE.match(core_db)
    return m.group(1) if m else None


def _ensembl_fungi_gff_url(
    row: dict[str, str],
    cache_dir: Path,
) -> tuple[str, str] | None:
    """Resolve an Ensembl Fungi GFF3 URL for one NCBI assembly row.

    Returns (url, filename) or None if the assembly has no Ensembl Fungi
    counterpart, the release cannot be discovered, or the species index is
    unreachable. The URL points at the whole-genome GFF
    (`<Species>.<assembly>.<release>.gff3.gz`); per-chromosome files are
    skipped because the downstream parser expects a single GFF per asm.
    """
    release = _ensembl_fungi_release(cache_dir)
    if release is None:
        return None
    index = _ensembl_fungi_species_index(cache_dir)
    if not index:
        return None
    acc = (row.get("assembly_accession") or "").strip()
    candidates: list[str] = []
    if acc:
        candidates.append(acc)
        candidates.append(acc.split(".", 1)[0])
    for key in ("taxid", "species_taxid"):
        v = (row.get(key) or "").strip()
        if v:
            candidates.append(v)
    rec: dict[str, str] | None = None
    for key in candidates:
        if key and key in index:
            rec = index[key]
            break
    if rec is None:
        return None
    species_dir = (rec.get("species") or "").strip().lower()
    assembly = (rec.get("assembly") or "").strip()
    if not species_dir or not assembly:
        return None
    # The species TSV occasionally records `assembly` values with embedded
    # whitespace ("version 1", "CBS 141442 assembly", "Neocallimastix sp. G1
    # v1.0"). Those map to filenames Ensembl publishes with spaces collapsed to
    # underscores; passing the raw string would also produce a urllib request
    # that fails with "URL can't contain control characters" *before* the HTTP
    # call is even attempted. Sanitise both URL path components defensively.
    def _safe(part: str) -> str:
        # Collapse any run of whitespace to a single underscore — matches the
        # actual filenames Ensembl serves and dodges the urllib pre-flight check.
        return re.sub(r"\s+", "_", part)

    species_dir = _safe(species_dir)
    assembly = _safe(assembly)
    # The whole-genome GFF filename is "<Species_Cap>.<assembly>.<release>.gff3.gz".
    # Capitalize only the first character (Ensembl convention: do NOT title-case
    # each word — strain suffixes like "_cen_pk113_7d_gca_000269885" must stay
    # lowercase).
    species_cap = species_dir[:1].upper() + species_dir[1:]
    filename = f"{species_cap}.{assembly}.{release}.gff3.gz"
    collection = _ensembl_fungi_collection(rec)
    if collection:
        url = f"{_ENSEMBL_FUNGI_FTP_BASE}/gff3/{collection}/{species_dir}/{filename}"
    else:
        url = f"{_ENSEMBL_FUNGI_FTP_BASE}/gff3/{species_dir}/{filename}"
    return url, filename


# Tally of per-asm annotation outcomes. download_ncbi_gene_annotation_source
# used to emit one stderr line per assembly without a NCBI GFF (1400+ lines
# in a 2000-assembly run, drowning real errors). The function now records
# (kind, asm_name) tuples here; the caller prints a single summary line.
ANNOTATION_SOURCE_TALLY: dict[str, list[str]] = {
    "ncbi_gff": [],
    "ensembl_fungi_gff": [],
    "ncbi_gbff": [],
    "none": [],
}


def reset_annotation_source_tally() -> None:
    for k in ANNOTATION_SOURCE_TALLY:
        ANNOTATION_SOURCE_TALLY[k].clear()


def annotation_source_summary() -> str:
    counts = {k: len(v) for k, v in ANNOTATION_SOURCE_TALLY.items()}
    total = sum(counts.values())
    if total == 0:
        return ""
    parts = [
        f"ncbi_gff={counts['ncbi_gff']}",
        f"ensembl_fungi_gff={counts['ensembl_fungi_gff']}",
        f"ncbi_gbff={counts['ncbi_gbff']}",
        f"none={counts['none']}",
    ]
    return f"gene_annotation_sources: total={total}  " + "  ".join(parts)


def download_ncbi_gene_annotation_source(
    row: dict[str, str],
    asm_name: str,
    dest_dir: Path,
    role_label: str,
    *,
    ensembl_cache_dir: Path | None = None,
) -> Path | None:
    """Download the best public gene-annotation source for one NCBI assembly.

    Preference order:
      1. NCBI GFF3 (genomic.gff.gz)               — when NCBI auto-annotated it
      2. Ensembl Fungi GFF3 (species_metadata)    — for GenBank-only assemblies
                                                    Ensembl re-annotates
      3. NCBI GenBank flatfile (genomic.gbff.gz)  — gene/CDS records embedded
                                                    in the deposited record

    A missing/empty annotation source is not a fatal download error: many
    GenBank assemblies are sequence-only. Per-assembly outcomes are recorded
    in ANNOTATION_SOURCE_TALLY so the caller can emit a single end-of-stage
    summary instead of 1000+ informational stderr lines.
    """
    gff_url, gff_filename = ncbi_download_targets(row, include_gff=True)[1]
    try:
        path = materialize_entry(gff_url, dest_dir / gff_filename, keep_gz=True)
        ANNOTATION_SOURCE_TALLY["ncbi_gff"].append(asm_name)
        return path
    except urllib.error.HTTPError as exc:
        if getattr(exc, "code", None) != 404:
            raise

    if ensembl_cache_dir is not None:
        ens = _ensembl_fungi_gff_url(row, ensembl_cache_dir)
        if ens is not None:
            ens_url, ens_filename = ens
            try:
                path = materialize_entry(ens_url, dest_dir / ens_filename, keep_gz=True)
                ANNOTATION_SOURCE_TALLY["ensembl_fungi_gff"].append(asm_name)
                return path
            except urllib.error.HTTPError as exc:
                # 404 here is normal (Ensembl may list the species but not
                # ship a GFF in `current`); any other HTTP error falls through
                # to the NCBI GBFF path rather than blocking the assembly.
                if getattr(exc, "code", None) not in {404, 403}:
                    sys.stderr.write(
                        f"[ensembl-fungi] {asm_name}: {exc} ({ens_url})\n"
                    )
            except Exception as exc:  # noqa: BLE001 - best-effort fallback
                sys.stderr.write(
                    f"[ensembl-fungi] {asm_name}: {exc}\n"
                )

    gbff_url, gbff_filename = ncbi_genbank_target(row)
    try:
        gbff = materialize_entry(gbff_url, dest_dir / gbff_filename, keep_gz=True)
    except urllib.error.HTTPError as exc:
        if getattr(exc, "code", None) == 404:
            ANNOTATION_SOURCE_TALLY["none"].append(asm_name)
            sys.stderr.write(
                f"[gene-annot] no public gene annotation source for {role_label}{asm_name} "
                "(NCBI has no GFF/GBFF; Ensembl Fungi did not match)\n"
            )
            return None
        raise
    ANNOTATION_SOURCE_TALLY["ncbi_gbff"].append(asm_name)
    return gbff


# ============================================================================
# One-shot caches for downstream multi-omics signal
#
# Both gene-expression (Expression Atlas differential) and ecological-trait
# (FungalTraits) data are species-keyed and slow to refetch — the EBI Atlas
# JSON endpoint is paginated and FungalTraits is hosted as a single CSV. We
# cache them under data_cache/ and reuse across panels; the prepare-step
# stitches per-panel `expression.tsv` from the cache so the existing
# auto-pickup in benchmark_real_data() needs no plumbing change.
# ============================================================================


_EXPRESSION_ATLAS_BASE = "https://www.ebi.ac.uk/gxa/json/experiments"


def _atlas_species_slug(species: str) -> str:
    return species.strip().lower().replace(" ", "_")


def fetch_expression_atlas_for_species(
    species_list: list[str],
    cache_dir: Path,
) -> dict[str, Path]:
    """Pull EBI Expression Atlas differential experiment summaries per species.

    For each species we hit /gxa/json/experiments?species=...&experimentType=differential
    once and persist the JSON list to cache_dir/<slug>.json. The JSON carries
    per-experiment accessions; the analyzer's expression.tsv only needs
    gene-level log2_fc / padj, which is sourced separately when an experiment
    accession is followed (TSV at /gxa/experiments/<acc>/download/all-analytics).
    The on-disk cache is the unit of reuse; we never refetch when the file
    exists, matching the phenotype-cache convention.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    out: dict[str, Path] = {}
    for species in species_list:
        if not species or species == ".":
            continue
        slug = _atlas_species_slug(species)
        cached = cache_dir / f"{slug}.json"
        if cached.exists() and cached.stat().st_size > 0:
            out[species] = cached
            continue
        params = urllib.parse.urlencode({
            "species": species,
            "experimentType": "differential",
        })
        url = f"{_EXPRESSION_ATLAS_BASE}?{params}"
        try:
            text = http_get_text(url)
        except Exception as exc:
            sys.stderr.write(
                f"[expression-atlas] {species!r}: lookup failed "
                f"({type(exc).__name__}: {exc}); skipping\n"
            )
            continue
        cached.write_text(text, encoding="utf-8")
        out[species] = cached
        time.sleep(0.5)
    return out


def _looks_like_atlas_analytics_tsv(text: str) -> bool:
    stripped = text.lstrip()
    if not stripped:
        return False
    if stripped.startswith("<") or "<!doctype html" in stripped[:200].lower():
        return False
    first = stripped.splitlines()[0] if stripped.splitlines() else ""
    return "\t" in first and "Gene ID" in first


def fetch_atlas_experiment_analytics(
    experiment_accession: str,
    cache_dir: Path,
) -> Path | None:
    """Download the per-gene analytics TSV (log2_fc + padj per contrast) for
    a single Expression Atlas experiment. Cached at
    cache_dir/<accession>.analytics.tsv; reused across runs.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    cached = cache_dir / f"{experiment_accession}.analytics.tsv"
    unavailable = cached.with_suffix(cached.suffix + ".unavailable")
    if unavailable.exists():
        return None
    if cached.exists() and cached.stat().st_size > 0:
        try:
            cached_text = cached.read_text(encoding="utf-8", errors="replace")
        except OSError:
            cached_text = ""
        if _looks_like_atlas_analytics_tsv(cached_text):
            return cached
        bad = cached.with_suffix(cached.suffix + ".invalid.html")
        try:
            cached.rename(bad)
        except OSError:
            try:
                cached.unlink()
            except OSError:
                pass
        sys.stderr.write(
            f"[expression-atlas] {experiment_accession}: cached analytics "
            f"was not TSV; quarantined as {bad.name}\n"
        )
        try:
            unavailable.write_text("cached analytics was not a TSV\n", encoding="utf-8")
        except OSError:
            pass
        return None
    url = (
        f"https://www.ebi.ac.uk/gxa/experiments/{experiment_accession}"
        f"/download/all-analytics?accessKey="
    )
    try:
        text = http_get_text(url)
    except Exception as exc:
        sys.stderr.write(
            f"[expression-atlas] {experiment_accession}: download failed "
            f"({type(exc).__name__}: {exc})\n"
        )
        return None
    if not _looks_like_atlas_analytics_tsv(text):
        bad = cached.with_suffix(cached.suffix + ".invalid.html")
        try:
            bad.write_text(text, encoding="utf-8")
        except OSError:
            pass
        try:
            unavailable.write_text(
                "analytics endpoint did not return a gene analytics TSV\n",
                encoding="utf-8",
            )
        except OSError:
            pass
        sys.stderr.write(
            f"[expression-atlas] {experiment_accession}: analytics endpoint "
            f"did not return a gene analytics TSV; cached as unavailable\n"
        )
        return None
    cached.write_text(text, encoding="utf-8")
    return cached


def assemble_expression_tsv(
    species_to_listing: dict[str, Path],
    species_to_query_asms: dict[str, list[str]],
    analytics_cache_dir: Path,
    out_path: Path,
    *,
    max_experiments_per_species: int = 1,
) -> int:
    """Materialise prepared_dir/expression.tsv from the cached Atlas data.

    Per panel-species we expand the cached experiment listing to at most
    max_experiments_per_species analytics tables, pivot each (gene_id, log2fc,
    padj) record into the analyzer's expected schema (query_asm, query_contig,
    gene_id, gene_name, log2_fc, padj, condition), and join across all
    query_asms that map to that species. Returns the number of rows written.

    The analyzer keys on query_asm, so emitting one block per query_asm
    avoids forcing the caller to learn the species → assembly mapping.
    """
    rows: list[dict[str, Any]] = []
    for species, listing_path in species_to_listing.items():
        try:
            listing = json.loads(listing_path.read_text(encoding="utf-8"))
        except Exception:
            continue
        experiments = listing.get("experiments") or []
        if not experiments:
            continue
        accessions = [e.get("experimentAccession") for e in experiments if e.get("experimentAccession")]
        accessions = accessions[:max_experiments_per_species]
        for acc in accessions:
            analytics = fetch_atlas_experiment_analytics(acc, analytics_cache_dir)
            if analytics is None:
                continue
            try:
                lines = analytics.read_text(encoding="utf-8").splitlines()
            except Exception:
                continue
            if not lines:
                continue
            header = lines[0].split("\t")
            try:
                gid_idx = header.index("Gene ID")
            except ValueError:
                continue
            gname_idx = header.index("Gene Name") if "Gene Name" in header else -1
            # Atlas analytics carry one (log2foldchange, p-value) pair per
            # contrast; the column header encodes the contrast name as
            # "<Contrast>.log2foldchange" / "<Contrast>.p-value". We pivot
            # the first contrast only — additional contrasts are ignored
            # to keep expression.tsv flat, matching the analyzer schema.
            contrast: str | None = None
            l2_idx: int | None = None
            p_idx: int | None = None
            for i, col in enumerate(header):
                if col.endswith(".log2foldchange"):
                    contrast = col[: -len(".log2foldchange")]
                    l2_idx = i
                    try:
                        p_idx = header.index(f"{contrast}.p-value")
                    except ValueError:
                        p_idx = None
                    break
            if contrast is None or l2_idx is None or p_idx is None:
                continue
            for raw in lines[1:]:
                parts = raw.split("\t")
                if len(parts) <= max(l2_idx, p_idx, gid_idx):
                    continue
                gene_id = parts[gid_idx].strip()
                if not gene_id:
                    continue
                gene_name = parts[gname_idx].strip() if gname_idx >= 0 else gene_id
                log2_fc = parts[l2_idx].strip()
                padj = parts[p_idx].strip()
                if not log2_fc or not padj:
                    continue
                for qasm in species_to_query_asms.get(species, []):
                    rows.append({
                        "query_asm": qasm,
                        "query_contig": ".",   # gene-level, not contig-resolved
                        "gene_id": gene_id,
                        "gene_name": gene_name,
                        "distance_bp": ".",
                        "log2_fc": log2_fc,
                        "padj": padj,
                        "condition": contrast,
                    })
    if not rows:
        return 0
    write_tsv(
        out_path,
        rows,
        ["query_asm", "query_contig", "gene_id", "gene_name",
         "distance_bp", "log2_fc", "padj", "condition"],
    )
    return len(rows)


_FUNGALTRAITS_URL = (
    "https://raw.githubusercontent.com/traitecoevo/fungaltraits/master/"
    "funtothefun.csv"
)


def fetch_fungaltraits_table(cache_dir: Path) -> Path | None:
    """Download the FungalTraits genus-level lifestyle/trait CSV once and
    cache. Returns the local path, or None on failure. Reused across panels
    via the persistent data_cache directory.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    cached = cache_dir / "fungaltraits.csv"
    if cached.exists() and cached.stat().st_size > 0:
        return cached
    try:
        text = http_get_text(_FUNGALTRAITS_URL)
    except Exception as exc:
        sys.stderr.write(
            f"[fungaltraits] download failed ({type(exc).__name__}: {exc}); "
            f"continuing without trait enrichment\n"
        )
        return None
    cached.write_text(text, encoding="utf-8")
    return cached


def write_ecological_summary_tsv(
    fungaltraits_csv: Path | None,
    species_to_query_asms: dict[str, list[str]],
    out_path: Path,
) -> int:
    """Project the cached FungalTraits CSV onto every panel species, writing
    a flat TSV the visualization report and biology analyzer can join on.
    Returns the number of rows written, 0 if no traits matched.
    """
    if fungaltraits_csv is None or not fungaltraits_csv.exists():
        return 0
    try:
        text = fungaltraits_csv.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return 0
    lines = text.splitlines()
    if len(lines) < 2:
        return 0
    # FungalTraits CSV is comma-delimited with a single header row.
    header = next(csv.reader([lines[0]]))
    by_genus: dict[str, dict[str, str]] = {}
    by_species: dict[str, dict[str, str]] = {}
    for raw in lines[1:]:
        parts = next(csv.reader([raw]))
        if len(parts) < len(header):
            parts.extend([""] * (len(header) - len(parts)))
        record = dict(zip(header, parts))
        genus = (record.get("GENUS") or record.get("Genus") or "").strip().lower()
        species = (record.get("SPECIES") or record.get("Species") or "").strip().lower()
        if not species:
            species = (record.get("speciesMatched") or record.get("species") or "").strip().replace("_", " ").lower()
        if not genus and species:
            genus = species.split()[0]
        trait_name = (record.get("trait_name") or "").strip().lower()
        trait_value = (record.get("value") or "").strip()
        if trait_name and trait_value:
            # The current traitecoevo/fungaltraits export is long-form
            # (speciesMatched, trait_name, value). Pivot the ecology fields we
            # use so the downstream join can consume it like the older wide
            # FungalTraits v1.2 CSV.
            pivoted = {
                "GENUS": genus,
                "SPECIES": species,
            }
            if trait_name == "trophic_mode_fg":
                pivoted["Trophic_mode"] = trait_value
                pivoted["primary_lifestyle"] = trait_value
            elif trait_name == "substrate":
                pivoted["Substrate"] = trait_value
            else:
                pivoted[trait_name] = trait_value
            record = pivoted
        if genus:
            by_genus.setdefault(genus, {}).update({k: v for k, v in record.items() if v})
        if species:
            by_species.setdefault(species, {}).update({k: v for k, v in record.items() if v})
    rows: list[dict[str, Any]] = []
    for species, qasms in species_to_query_asms.items():
        species_low = species.strip().lower()
        genus_low = species_low.split(" ")[0] if species_low else ""
        record = by_species.get(species_low) or by_genus.get(genus_low) or {}
        if not record:
            continue
        primary_lifestyle = (
            record.get("primary_lifestyle")
            or record.get("Primary_lifestyle")
            or record.get("PRIMARY_LIFESTYLE")
            or "."
        )
        secondary_lifestyle = (
            record.get("Secondary_lifestyle")
            or record.get("secondary_lifestyle")
            or "."
        )
        substrate = (
            record.get("Plant_pathogenic_capacity_template")
            or record.get("Substrate")
            or record.get("substrate")
            or "."
        )
        trophic_mode = (
            record.get("Trophic_mode")
            or record.get("trophic_mode")
            or "."
        )
        for qasm in qasms:
            rows.append({
                "query_asm": qasm,
                "species": species,
                "primary_lifestyle": primary_lifestyle,
                "secondary_lifestyle": secondary_lifestyle,
                "trophic_mode": trophic_mode,
                "substrate_or_host": substrate,
            })
    if not rows:
        return 0
    write_tsv(
        out_path,
        rows,
        ["query_asm", "species", "primary_lifestyle", "secondary_lifestyle",
         "trophic_mode", "substrate_or_host"],
    )
    return len(rows)


def materialize_entry(url_or_path: str, dest: Path, keep_gz: bool = True) -> Path:
    if url_or_path.startswith(("http://", "https://", "ftp://")):
        path = http_download(url_or_path, dest)
    else:
        src = Path(url_or_path)
        dest.parent.mkdir(parents=True, exist_ok=True)
        try:
            same_file = src.resolve() == dest.resolve()
        except OSError:
            same_file = False
        if not same_file:
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


MYCOSV_HITS_FIELDS = [
    "query_asm", "query_contig", "type", "ref_asm", "ref_contig",
    "ref_pos", "ref_end", "pos", "end", "svlen", "block_score", "anchors",
    "genotype", "gq", "annotation", "alignment_mode", "query_mode",
    "fused_posterior_alt", "fused_logodds_alt", "fused_effective_depth",
    "fused_layers", "read_support",
]


def write_mycosv_failure_outputs(out_prefix: Path, reason: str) -> dict[str, str]:
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    vcf_path = out_prefix.with_suffix(".vcf")
    hits_path = out_prefix.with_suffix(".hits.tsv")
    gfa_path = out_prefix.with_suffix(".gfa")
    safe_reason = reason.replace("\n", " ").replace("\t", " ")
    vcf_path.write_text(
        "##fileformat=VCFv4.3\n"
        "##source=fungi_graphsv_tol_v3\n"
        f"##mycosv_status=failed:{safe_reason}\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
        encoding="utf-8",
    )
    write_tsv(hits_path, [{
        "query_asm": ".",
        "query_contig": ".",
        "type": "NO_CALLS",
        "ref_asm": ".",
        "ref_contig": ".",
        "ref_pos": 0,
        "ref_end": 0,
        "pos": 0,
        "end": 0,
        "svlen": 0,
        "block_score": 0,
        "anchors": 0,
        "genotype": "./.",
        "gq": 0,
        "annotation": f"MYCOSV_FAILED:{safe_reason}",
        "alignment_mode": "diagnostic",
        "query_mode": ".",
        "fused_posterior_alt": 0,
        "fused_logodds_alt": 0,
        "fused_effective_depth": 0,
        "fused_layers": 0,
        "read_support": -1,
    }], MYCOSV_HITS_FIELDS)
    gfa_path.write_text(f"H\tVN:Z:1.0\tST:Z:MYCOSV_FAILED\tRS:Z:{safe_reason}\n", encoding="utf-8")
    return {"vcf": str(vcf_path), "hits": str(hits_path), "gfa": str(gfa_path)}


def snapshot_mycosv_outputs(out_prefix: Path) -> dict[Path, bytes]:
    """Keep prior MycoSV artifacts in memory before a rerun truncates them.

    The C++ binary streams directly to calls.vcf/calls.hits.tsv/calls.gfa. If a
    rerun is killed mid-write, the old successful callset used to be lost and
    replaced by a partial or failure-header file. These artifacts are small
    enough for the real-data benchmark panels, and preserving them makes failed
    reruns diagnosable instead of destructive.
    """
    snapshots: dict[Path, bytes] = {}
    for path in (
        out_prefix.with_suffix(".vcf"),
        out_prefix.with_suffix(".hits.tsv"),
        out_prefix.with_suffix(".gfa"),
        out_prefix.parent / (out_prefix.name + ".multisample.vcf"),
    ):
        if path.exists() and path.is_file() and path.stat().st_size > 0:
            snapshots[path] = path.read_bytes()
    return snapshots


def restore_mycosv_outputs(snapshots: dict[Path, bytes]) -> None:
    for path, payload in snapshots.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)


def vcf_data_record_count(path: Path) -> int:
    if not path.exists() or not path.is_file() or path.stat().st_size == 0:
        return 0
    n = 0
    try:
        with open_text_auto(path) as fh:
            for line in fh:
                if line.strip() and not line.startswith("#"):
                    n += 1
    except OSError:
        return 0
    return n


def promote_hierarchical_checkpoint(out_prefix: Path) -> dict[str, str] | None:
    """Use completed per-contig checkpoint output as the canonical callset.

    A long panel can be killed after hierarchical_call_assembly has flushed
    thousands of per-contig calls but before the final calls.vcf is complete.
    In that case the checkpoint is the best available pangenome caller output
    and should feed the benchmark/reporting layers instead of a stale or empty
    canonical VCF.
    """
    hier_vcf = out_prefix.parent / (out_prefix.name + ".hierarchical.vcf")
    hier_hits = out_prefix.parent / (out_prefix.name + ".hierarchical.hits.tsv")
    if vcf_data_record_count(hier_vcf) == 0:
        return None
    canonical_vcf = out_prefix.with_suffix(".vcf")
    canonical_hits = out_prefix.with_suffix(".hits.tsv")
    canonical_vcf.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(hier_vcf, canonical_vcf)
    if hier_hits.exists() and hier_hits.stat().st_size > 0:
        shutil.copy2(hier_hits, canonical_hits)
    return {
        "vcf": str(canonical_vcf),
        "hits": str(canonical_hits),
        "gfa": str(out_prefix.with_suffix(".gfa")),
    }


_GFF_GENE_TYPES: frozenset[str] = frozenset({
    # NCBI / Ensembl Fungi gene-level feature types we want to surface to the
    # SV biology analyzer. ncRNA / tRNA / rRNA are deliberately excluded — the
    # candidate scoring already discriminates protein-coding loci.
    "gene", "protein_coding_gene", "pseudogene",
})


def _parse_gff_attributes(attr_field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for entry in attr_field.strip().split(";"):
        entry = entry.strip()
        if not entry or "=" not in entry:
            continue
        k, _, v = entry.partition("=")
        out[k.strip().lower()] = urllib.parse.unquote(v.strip())
    return out


_GBFF_GENE_TYPES = frozenset({
    "gene",
    "protein_coding_gene",
    "pseudogene",
    "cds",
})


def _annotation_source_kind(path: Path) -> str:
    name = path.name.lower()
    if name.endswith(".gbff") or name.endswith(".gbff.gz"):
        return "gbff"
    return "gff"


def _parse_genbank_location(location: str) -> tuple[int, int, str] | None:
    nums = [int(x) for x in re.findall(r"\d+", location)]
    if not nums:
        return None
    start = min(nums)
    end = max(nums)
    strand = "-" if "complement" in location.lower() else "+"
    return start, end, strand


def _clean_genbank_qualifier_value(value: str) -> str:
    value = value.strip()
    if value.startswith('"') and value.endswith('"') and len(value) >= 2:
        value = value[1:-1]
    return value.replace('""', '"')


def _flush_genbank_feature(
    rows: list[dict[str, Any]],
    seen: set[tuple[str, str, str]],
    asm_name: str,
    contig: str,
    feature: dict[str, Any] | None,
) -> None:
    if not feature or not contig:
        return
    ftype = str(feature.get("type", "")).lower()
    if ftype not in _GBFF_GENE_TYPES:
        return
    loc = _parse_genbank_location(str(feature.get("location", "")))
    if loc is None:
        return
    start, end, strand = loc
    quals: dict[str, str] = feature.get("qualifiers", {})
    gene_id = (
        quals.get("locus_tag")
        or quals.get("gene")
        or quals.get("old_locus_tag")
        or quals.get("protein_id")
        or quals.get("db_xref")
        or ""
    )
    if not gene_id:
        gene_id = f"{ftype}:{contig}:{start}-{end}:{strand}"
    key = (asm_name, contig, gene_id)
    if key in seen:
        return
    seen.add(key)
    rows.append({
        "query_asm": asm_name,
        "query_contig": contig,
        "gene_id": gene_id,
        "gene_name": quals.get("gene") or quals.get("locus_tag") or gene_id,
        "start": start,
        "end": end,
        "strand": strand,
        "biotype": quals.get("gene_biotype") or quals.get("gbkey") or ftype,
        "product": quals.get("product", ""),
    })


def gbff_to_gene_annotations(gbff_paths: list[tuple[str, Path]]) -> list[dict[str, Any]]:
    """Parse NCBI GenBank flatfiles into the same row schema as GFF.

    Many GenBank fungal assemblies expose a GBFF feature table but no public
    GFF. This lightweight parser extracts gene/CDS spans without pulling in a
    BioPython dependency for the benchmark environment.
    """
    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str]] = set()
    feature_re = re.compile(r"^     (\S+)\s+(.+)$")
    qualifier_re = re.compile(r"^                     /([^=]+)(?:=(.*))?$")
    continuation_re = re.compile(r"^                     ([^/].*)$")
    for asm_name, gbff_path in gbff_paths:
        if not gbff_path.exists():
            continue
        opener = gzip.open if gbff_path.suffix == ".gz" else open
        contig = ""
        in_features = False
        current: dict[str, Any] | None = None
        current_qualifier: str | None = None
        try:
            with opener(gbff_path, "rt", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    line = line.rstrip("\n")
                    if line.startswith("LOCUS"):
                        _flush_genbank_feature(rows, seen, asm_name, contig, current)
                        current = None
                        current_qualifier = None
                        in_features = False
                        parts = line.split()
                        contig = parts[1] if len(parts) > 1 else contig
                        continue
                    if line.startswith("VERSION"):
                        parts = line.split()
                        if len(parts) > 1:
                            contig = parts[1]
                        continue
                    if line.startswith("FEATURES"):
                        in_features = True
                        continue
                    if line.startswith("ORIGIN") or line == "//":
                        _flush_genbank_feature(rows, seen, asm_name, contig, current)
                        current = None
                        current_qualifier = None
                        in_features = False
                        continue
                    if not in_features:
                        continue
                    match = feature_re.match(line)
                    if match:
                        _flush_genbank_feature(rows, seen, asm_name, contig, current)
                        current = {
                            "type": match.group(1),
                            "location": match.group(2).strip(),
                            "qualifiers": {},
                        }
                        current_qualifier = None
                        continue
                    if current is None:
                        continue
                    qmatch = qualifier_re.match(line)
                    if qmatch:
                        qkey = qmatch.group(1).strip().lower()
                        qval = _clean_genbank_qualifier_value(qmatch.group(2) or "")
                        current["qualifiers"][qkey] = qval
                        current_qualifier = qkey
                        continue
                    cmatch = continuation_re.match(line)
                    if cmatch and current_qualifier:
                        prev = current["qualifiers"].get(current_qualifier, "")
                        cont = _clean_genbank_qualifier_value(cmatch.group(1))
                        current["qualifiers"][current_qualifier] = f"{prev} {cont}".strip()
        except OSError as exc:
            sys.stderr.write(f"[gene-annot] skip {gbff_path}: {type(exc).__name__}: {exc}\n")
            continue
    return rows


def gff_to_gene_annotations(gff_paths: list[tuple[str, Path]]) -> list[dict[str, Any]]:
    """Convert a list of (asm_name, annotation path) tuples into the row schema
    expected by analyze_new_biology_candidates.load_gene_annotations.

    Output columns: query_asm, query_contig, gene_id, gene_name, start, end.
    The asm_name is stored verbatim so the analyzer's per-(asm, contig) lookup
    finds genes for either ref-coordinate or query-coordinate breakpoints
    (the analyzer falls back to '.' when no asm-keyed row matches).
    """
    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str]] = set()
    gbff_paths: list[tuple[str, Path]] = []
    for asm_name, gff_path in gff_paths:
        if _annotation_source_kind(gff_path) == "gbff":
            gbff_paths.append((asm_name, gff_path))
            continue
        if not gff_path.exists():
            continue
        opener = gzip.open if gff_path.suffix == ".gz" else open
        try:
            with opener(gff_path, "rt", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if not line or line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 9:
                        continue
                    contig, _src, ftype, start_s, end_s, _score, strand, _phase, attrs = parts[:9]
                    if ftype.lower() not in _GFF_GENE_TYPES:
                        continue
                    try:
                        start = int(start_s)
                        end = int(end_s)
                    except ValueError:
                        continue
                    if end < start:
                        start, end = end, start
                    parsed = _parse_gff_attributes(attrs)
                    gene_id = (
                        parsed.get("id")
                        or parsed.get("locus_tag")
                        or parsed.get("gene_id")
                        or parsed.get("name")
                        or ""
                    )
                    if not gene_id:
                        continue
                    key = (asm_name, contig, gene_id)
                    if key in seen:
                        continue
                    seen.add(key)
                    rows.append({
                        "query_asm": asm_name,
                        "query_contig": contig,
                        "gene_id": gene_id,
                        "gene_name": parsed.get("name") or parsed.get("gene") or gene_id,
                        "start": start,
                        "end": end,
                        "strand": strand if strand in {"+", "-"} else ".",
                        "biotype": parsed.get("biotype") or parsed.get("gene_biotype") or ftype,
                        "product": parsed.get("product", ""),
                    })
        except OSError as exc:
            sys.stderr.write(f"[gene-annot] skip {gff_path}: {type(exc).__name__}: {exc}\n")
            continue
    if gbff_paths:
        rows.extend(gbff_to_gene_annotations(gbff_paths))
    return rows


GENE_ANNOTATION_COLUMNS = [
    "query_asm", "query_contig", "gene_id", "gene_name",
    "start", "end", "strand", "biotype", "product",
]


def write_gene_annotations_tsv(out_path: Path, gff_paths: list[tuple[str, Path]]) -> Path | None:
    rows = gff_to_gene_annotations(gff_paths)
    if not rows:
        sys.stderr.write(
            f"[gene-annot] no gene records parsed from {len(gff_paths)} annotation source(s); "
            f"skipping {out_path.name}\n"
        )
        return None
    write_tsv(out_path, rows, GENE_ANNOTATION_COLUMNS)
    sys.stderr.write(
        f"[gene-annot] wrote {len(rows)} gene records to {out_path}\n"
    )
    return out_path


def stream_gene_annotations_to_tsv(
    out_path: Path,
    gff_pairs: list[tuple[str, Path]],
    asm_aliases: dict[str, set[str]],
    ref_to_queries: dict[str, list[str]],
    progress_every: int = 100,
) -> int:
    """Stream gene rows source-by-source straight to disk with alias expansion.

    The non-streaming path used to collect every parsed gene from every
    GFF/GBFF source plus its alias duplicates into a single in-memory list
    before calling write_tsv. With 2000 GBFF sources (~10K genes each) and
    per-row alias expansion to ~3-5 owners, that materialised 60M+ dicts and
    OOM-killed prepare_million_real under the 12 GiB cgroup cap before any
    TSV bytes hit disk. Streaming bounds memory to one source (~10K rows)
    at a time.

    Returns the number of rows written (header excluded).
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    sources_total = len(gff_pairs)
    sources_processed = 0
    sources_with_records = 0
    total_rows = 0
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=GENE_ANNOTATION_COLUMNS, delimiter="\t")
        writer.writeheader()
        for asm_name, annot_path in gff_pairs:
            source_rows = gff_to_gene_annotations([(asm_name, annot_path)])
            sources_processed += 1
            if source_rows:
                sources_with_records += 1
                owners = set(asm_aliases.get(asm_name, {asm_name}))
                for q_asm in ref_to_queries.get(asm_name, []):
                    owners.update(asm_aliases.get(q_asm, {q_asm}))
                # Per-source (owner, contig, gene_id) dedup is sufficient: each
                # source has a unique asm_name so cross-source duplicates are
                # impossible by construction.
                seen: set[tuple[str, str, str]] = set()
                for gene_row in source_rows:
                    contig = gene_row["query_contig"]
                    gene_id = gene_row["gene_id"]
                    for owner in owners:
                        key = (owner, contig, gene_id)
                        if key in seen:
                            continue
                        seen.add(key)
                        out_row = dict(gene_row)
                        out_row["query_asm"] = owner
                        writer.writerow(out_row)
                        total_rows += 1
            if progress_every and sources_processed % progress_every == 0:
                print(
                    f"      gene_annotations: parsed {sources_processed}/{sources_total} sources "
                    f"({sources_with_records} non-empty), wrote {total_rows} rows so far",
                    flush=True,
                )
    return total_rows


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
    """Build an ENA portal URL that returns public read runs for a species.

    Used by prepare_from_ncbi when --query-mode is short-reads or long-reads
    and the panel presets only describe a species (no read accession).

    Note: the ENA portal API splits filtering and accession-based lookup
    across two endpoints. /filereport accepts a single ``accession=`` and
    rejects ``query=`` with HTTP 400. Filter-by-name lives on /search, so
    species lookups go through that endpoint.

    The scientific_name filter works against any rank ENA stores, so a genus
    name like 'Rhizophagus' also resolves — important for AMF where
    species-level assignments are patchy.
    """
    query = urllib.parse.urlencode({
        # ENA portal expects this quoted exactly as: scientific_name="X"
        "query": f'scientific_name="{species}"',
        "result": "read_run",
        "fields": ",".join(ENA_FILEREPORT_FIELDS),
        "format": "tsv",
        "download": "false",
        "limit": str(max_rows),
    })
    return f"https://www.ebi.ac.uk/ena/portal/api/search?{query}"


def fetch_ena_read_runs_by_species(species: str, max_rows: int = 200) -> list[dict[str, str]]:
    # Try the preset name first, then SPECIES_ALIASES (NCBI/ENA renames such as
    # Candida glabrata -> Nakaseomyces glabratus). ENA's filereport keys on the
    # current scientific_name, so a panel still using the legacy name otherwise
    # silently returns 0 runs even though reads exist under the new genus.
    candidates: list[str] = [species]
    candidates.extend(SPECIES_ALIASES.get(species.strip().lower(), []))
    seen: set[str] = set()
    for cand in candidates:
        key = cand.strip().lower()
        if not key or key in seen:
            continue
        seen.add(key)
        try:
            rows = parse_ena_filereport_text(http_get_text(ena_filereport_species_url(cand, max_rows)))
        except Exception as exc:
            sys.stderr.write(
                f"[reads-mode] ENA species lookup error for {cand!r}: "
                f"{type(exc).__name__}: {exc}\n"
            )
            continue
        if rows:
            if cand != species:
                sys.stderr.write(
                    f"[reads-mode] resolved {species!r} -> {cand!r} via SPECIES_ALIASES\n"
                )
            return rows
    return []


def sequence_kind_from_name(name: str) -> str:
    lower = name.lower()
    if lower.endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq")):
        return "fastq"
    if lower.endswith((".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".fa", ".fna")):
        return "fasta"
    # Older ENA submissions for PacBio RSII/Sequel I expose `.bas.h5`,
    # `.bax.h5`, and `.metadata.xml` payloads via submitted_ftp. Earlier this
    # function defaulted unknown extensions to "fastq", so the downloader
    # happily concatenated the HDF5/XML bytes into a `.fastq` file and every
    # downstream comparator died with `'utf-8' codec can't decode byte ...`.
    # Return "unknown" so callers can filter explicitly.
    return "unknown"


def preferred_platforms_for_mode(query_mode: str) -> set[str]:
    mode = query_mode.strip().lower()
    if mode == "short-reads":
        return {"ILLUMINA", "BGISEQ", "DNBSEQ", "ION_TORRENT"}
    if mode == "long-reads":
        return {"OXFORD_NANOPORE", "PACBIO_SMRT"}
    return set()


def _is_pacbio_hifi(row: dict[str, str]) -> bool:
    """Return True when an ENA run row represents PacBio HiFi (CCS) reads.

    Detection order:
    1. library_strategy contains a known HiFi keyword ("CCS", "HiFi", "Hi-Fi").
    2. instrument_model matches Revio or Sequel IIe — both are HiFi-only.
    3. Sequel II with WGS strategy defaults to HiFi; Sequel II CLR runs
       typically carry "CLR" in library_strategy.
    """
    platform = (row.get("instrument_platform") or "").upper()
    if "PACBIO" not in platform and "SMRT" not in platform:
        return False
    strategy = (row.get("library_strategy") or "").lower()
    if any(kw in strategy for kw in _HIFI_STRATEGY_KW):
        return True
    model = (row.get("instrument_model") or "").lower()
    if any(kw in model for kw in _HIFI_MODEL_KW):
        return True
    # Sequel II without an explicit CLR flag → assume modern HiFi workflow.
    if ("sequel ii" in model or "sequel 2" in model) and "clr" not in strategy:
        return True
    return False


def _long_read_platform_score(row: dict[str, str]) -> int:
    """Rank ENA long-read runs for selection priority (higher = preferred).

    3 — PacBio HiFi (Revio / Sequel IIe / Sequel II CCS)
          Highest per-read accuracy; minimap2 map-hifi + sniffles2 / cuteSV.
    2 — ONT PromethION or GridION
          High-depth, likely R10.4.1 simplex in recent submissions (~Q20).
    1 — ONT MinION / Mk1C
          Valid long-read data; recent kits may carry R10.4.1 chemistry.
    0 — PacBio CLR (RS II, Sequel I)
          Lower base accuracy; minimap2 map-pb still produces usable BAMs.
   -1 — Any other platform that passed the long-reads filter.
    """
    if _is_pacbio_hifi(row):
        return 3
    platform = (row.get("instrument_platform") or "").upper()
    model = (row.get("instrument_model") or "").lower()
    if "OXFORD_NANOPORE" in platform:
        if "promethion" in model or "gridion" in model:
            return 2
        return 1  # MinION, Mk1C, or unspecified ONT model
    if "PACBIO" in platform or "SMRT" in platform:
        return 0  # CLR fallback
    return -1


def _ena_read_count(row: dict[str, str]) -> int:
    """Return the run's reported read_count, or 0 if missing/non-numeric."""
    raw = (row.get("read_count") or "").strip()
    if not raw:
        return 0
    try:
        return int(raw)
    except ValueError:
        return 0


def filter_ena_rows_for_mode(rows: list[dict[str, str]], query_mode: str) -> list[dict[str, str]]:
    preferred = preferred_platforms_for_mode(query_mode)
    if not preferred:
        return rows
    matched = [row for row in rows if (row.get("instrument_platform") or "").upper() in preferred]
    # Hard filter: an Illumina run must NOT be returned for long-reads mode (and
    # vice versa). The previous `matched or rows` fallback caused species with
    # no long-read submissions (e.g. L. kluyveri) to receive the same Illumina
    # accessions for both short-reads and long-reads modes — the resulting
    # files were byte-identical, mislabeled, and broke downstream platform
    # detection in the SV callers.
    # Drop runs with read_count below the noise floor — saw `SRR33624766: 1
    # reads` in production, which is a sentinel ENA upload and aligns to nothing.
    # When read_count is absent (older submissions), keep the row so we don't
    # over-filter species with unreliable metadata.
    matched = [row for row in matched
               if _ena_read_count(row) == 0 or _ena_read_count(row) >= _ENA_MIN_READS]
    result = matched
    # For long reads, rank so PacBio HiFi > ONT PromethION > ONT MinION > PacBio CLR.
    if query_mode == "long-reads":
        result = sorted(result, key=_long_read_platform_score, reverse=True)
    return result


def _ena_fastq_bytes(row: dict[str, str]) -> int:
    total = 0
    for raw in split_values(row.get("fastq_bytes", "")):
        try:
            total += int(raw)
        except ValueError:
            continue
    return total


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
    if query_mode == "long-reads":
        filtered = sorted(
            filtered,
            key=lambda row: (-_long_read_platform_score(row), _ena_fastq_bytes(row) or 10**18),
        )
    elif query_mode == "short-reads":
        filtered = sorted(filtered, key=lambda row: _ena_fastq_bytes(row) or 10**18)
    urls: list[str] = []
    meta_rows: list[dict[str, str]] = []
    picked_runs = 0
    for row in filtered:
        # Only accept ENA's curated FASTQ mirrors. submitted_ftp can carry
        # `.bas.h5`, `.bax.h5`, and `.metadata.xml` for older PacBio runs (e.g.
        # ERR3500124), which the previous fallback happily concatenated into a
        # `.fastq` file — every downstream caller then died with a UTF-8
        # decode error. If fastq_ftp is empty we have no usable reads.
        row_urls = [
            item for item in split_values(row.get("fastq_ftp", ""))
            if sequence_kind_from_name(item) == "fastq"
        ]
        if not row_urls:
            continue
        picked_runs += 1
        normalized_urls = [
            normalise_download_url(item if looks_like_url(item) else f"https://{item}")
            for item in row_urls
        ]
        urls.extend(normalized_urls)
        meta_rows.append({
            "run_accession": row.get("run_accession", "."),
            "study_accession": row.get("study_accession", "."),
            "sample_accession": row.get("sample_accession", "."),
            "scientific_name": row.get("scientific_name", "."),
            "instrument_platform": row.get("instrument_platform", "."),
            "library_layout": row.get("library_layout", "."),
            "library_strategy": row.get("library_strategy", "."),
            "source_url": ena_filereport_url(row.get("run_accession") or row.get("study_accession") or row.get("sample_accession") or "."),
            # Keep the exact file group for this run. Callers that retry
            # run-by-run must preserve paired-end mates instead of zipping one
            # metadata row against a flattened URL list.
            "selected_urls": ";".join(normalized_urls),
        })
        if max_runs > 0 and picked_runs >= max_runs:
            break
    return urls, meta_rows


def selected_urls_from_ena_meta(meta: dict[str, str]) -> list[str]:
    return split_values(meta.get("selected_urls", ""))


def merge_sequence_sources(sources: list[str], dest_prefix: Path) -> Path:
    # Byte-mode concatenation. ENA mirrors occasionally serve compressed FASTQs
    # whose extension does not match their magic bytes (bz2/zstd/partial), and
    # the previous text-mode merge crashed at the first non-UTF-8 byte (seen in
    # production as "'utf-8' codec can't decode byte 0x89" on Rhizophagus
    # long-reads). Treating the payload as opaque bytes is safe because every
    # downstream consumer (mycosv, samtools, minimap2, comparators) reads
    # gzip-/text- formats via their own auto-detection.
    if not sources:
        raise ValueError("No sequence sources were provided")
    kind = sequence_kind_from_name(sources[0])
    if kind == "unknown":
        raise ValueError(
            f"Cannot identify sequence kind from filename: {sources[0]}. "
            "Expected .fastq[.gz]/.fq[.gz]/.fasta[.gz]/.fa[.gz]/.fna[.gz]."
        )
    suffix = ".fastq" if kind == "fastq" else ".fasta"
    out_path = dest_prefix.with_suffix(suffix)
    if out_path.exists() and out_path.stat().st_size > 0:
        if kind == "fastq":
            try:
                _validate_fastq_payload(out_path)
            except ValueError as exc:
                # Cached file from a previous run is corrupt (e.g. PacBio
                # bas.h5 / metadata.xml saved as `.fastq`). Drop it so the
                # caller can retry against a different ENA accession instead
                # of silently re-using a 10 GB blob of garbage forever.
                sys.stderr.write(f"[reads-mode] dropping corrupt cached payload: {exc}\n")
                # _validate_fastq_payload already unlinked the file before
                # raising — fall through to the re-download path.
        if out_path.exists() and out_path.stat().st_size > 0:
            return out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("wb") as out_fh:
        for idx, src in enumerate(sources):
            parsed = urllib.parse.urlparse(src)
            source_name = Path(parsed.path).name or f"source_{idx + 1}{suffix}"
            basename = f"part_{idx + 1}_{source_name}"
            local_part = materialize_entry(src, dest_prefix.parent / basename, keep_gz=False)
            part_kind = sequence_kind_from_name(local_part.name)
            if part_kind != kind:
                raise ValueError(f"Mixed sequence kinds in one query input: {kind} vs {part_kind} from {src}")
            part_size = local_part.stat().st_size if local_part.exists() else 0
            if part_size < 20:
                # Bail on empty / HTML-error-page payloads (gzip header alone
                # is 10 bytes; a one-record gzipped FASTQ is ~40-100 bytes,
                # so 20 catches "obviously broken" without rejecting tiny
                # fixtures). The previous silent skip let mycosv read 0
                # sequences and still report success.
                if local_part != out_path and local_part.parent == dest_prefix.parent and local_part.exists():
                    local_part.unlink()
                raise ValueError(
                    f"Downloaded part looks truncated ({part_size} bytes) for {src}"
                )
            try:
                with local_part.open("rb") as in_fh:
                    # Strip a leading gzip header by streaming through gzip if the
                    # part is gz-magic but the merged output is plain text. This
                    # preserves the "one plain FASTQ" output contract for the C++
                    # binary while still tolerating gz parts (materialize_entry
                    # with keep_gz=False already gunzips, so this is belt-and-
                    # suspenders for mixed mirrors).
                    if in_fh.read(2) == b"\x1f\x8b":
                        in_fh.seek(0)
                        with gzip.open(in_fh, "rb") as gz_fh:
                            shutil.copyfileobj(gz_fh, out_fh)
                    else:
                        in_fh.seek(0)
                        shutil.copyfileobj(in_fh, out_fh)
            except (EOFError, OSError, gzip.BadGzipFile):
                # The cached part is a corrupt/truncated gzip — drop it so the
                # next pool URL (or a future re-run) re-downloads instead of
                # repeatedly failing on the same broken bytes.
                if local_part != out_path and local_part.parent == dest_prefix.parent and local_part.exists():
                    local_part.unlink()
                raise
            if local_part != out_path and local_part.parent == dest_prefix.parent and local_part.exists():
                local_part.unlink()
    if kind == "fastq":
        _validate_fastq_payload(out_path)
    return out_path


def _validate_fastq_payload(path: Path) -> None:
    """Sanity-check that a freshly merged FASTQ actually starts with `@`.

    ENA `submitted_ftp` and a few legacy mirrors occasionally serve PacBio
    `.bas.h5` / `.metadata.xml` payloads (XML preamble + HDF5 binary) that
    look superficially like a generic blob. The previous code wrote those
    bytes to a `.fastq` file and every downstream comparator then died with
    `'utf-8' codec can't decode byte 0x89 in position N: invalid start byte`.
    Rather than emit a corrupt file, fail fast so the caller can drop the
    run and pick a different ENA accession.
    """
    try:
        with path.open("rb") as fh:
            head = fh.read(4)
    except OSError:
        return
    if not head:
        return
    if head[:2] == b"\x1f\x8b":  # gzipped FASTQ — content checked at consumption.
        return
    if head[:1] != b"@":
        path.unlink(missing_ok=True)
        snippet = head.decode("ascii", errors="replace") if all(0x20 <= b < 0x7F or b in {0x09, 0x0A, 0x0D} for b in head) else head.hex()
        raise ValueError(
            f"Downloaded FASTQ at {path} does not start with '@' (got {snippet!r}); "
            "ENA likely returned a non-FASTQ payload (e.g. PacBio bas.h5 or metadata.xml)."
        )


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
        dest_name = cached_filename_for_source(fasta_src, f"{asm_name}.fa")
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

    # When the query came from an ENA read run, propagate the platform /
    # layout / accession from the first selected meta row so the manifest
    # schema lines up with the prepare_from_ncbi / prepare_million_real
    # outputs. The benchmark step uses instrument_platform to gate reads-
    # mode aligner choice (5008/5086) and emitting the column unconditionally
    # avoids a schema mismatch when callers concat manifests across paths.
    first_ena_meta: dict[str, str] = {}
    for sr in source_rows:
        if sr.get("source_type") == "ena_read_run":
            first_ena_meta = sr
            break
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
        "instrument_platform": row.get("instrument_platform") or first_ena_meta.get("instrument_platform") or ".",
        "library_layout": row.get("library_layout") or first_ena_meta.get("library_layout") or ".",
        "run_accession": row.get("run_accession") or first_ena_meta.get("source_accession") or ".",
    }
    return query_row, str(local_path), source_rows


def prepare_from_custom_manifest(args: argparse.Namespace) -> int:
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_base = data_cache_base(args, out_dir)
    cache_base.mkdir(parents=True, exist_ok=True)
    refs_dir = cache_base / "refs"
    queries_dir = cache_base / "queries"
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
            dest_name = cached_filename_for_source(fasta_src, f"{asm_name}.fa")
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
            "instrument_platform", "library_layout", "run_accession",
        ],
    )
    write_tsv(
        out_dir / "public_data_links.tsv",
        source_link_rows,
        ["query_asm", "role", "query_mode", "source_type", "source_accession", "source_url", "local_path", "species"],
    )
    write_public_resource_links(out_dir / "public_resource_links.tsv")
    (out_dir / "prepare_summary.json").write_text(
        json.dumps(
            {
                "mode": "custom",
                "ref_count": len(ref_manifest_rows),
                "query_count": len(query_rows),
                "data_cache_dir": str(cache_base),
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    print(f"prepared\trefs={len(ref_manifest_rows)}\tqueries={len(query_rows)}\tmode=custom")
    return 0


def prepare_from_ncbi(args: argparse.Namespace) -> int:
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_base = data_cache_base(args, out_dir)
    cache_base.mkdir(parents=True, exist_ok=True)
    reset_annotation_source_tally()
    all_rows, assembly_summary_caches = fetch_ncbi_assembly_rows(args.source, cache_base)

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
                "source_catalog": row.get("_catalog_source", args.source),
            })
    else:
        for sel in selectors:
            matches = select_species_rows(all_rows, sel["species"], args.max_assemblies_per_species)
            if not matches:
                # Surface this loudly so the operator does not silently lose a
                # panel reference (e.g. Leptosphaeria maculans / Ustilago
                # maydis disappeared from a previous run because no entry in
                # assembly_summary started with those exact names — the
                # species are present in NCBI but as `Plenodomus` synonyms or
                # `Mycosarcoma maydis` / `[U.] maydis`. Add to SPECIES_ALIASES
                # to recover them.
                sys.stderr.write(
                    f"[panel-select] WARNING: no NCBI assembly matched species "
                    f"{sel['species']!r} for panel scenario={sel.get('scenario', '.')!r}. "
                    f"Add a SPECIES_ALIASES entry if NCBI uses a different name today.\n"
                )
                sys.stderr.flush()
                continue
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
                    "source_catalog": row.get("_catalog_source", args.source),
                })

    if not selected_rows:
        raise ValueError("No NCBI fungal assemblies matched the requested panel/species selection")

    write_tsv(
        out_dir / "selected_catalog.tsv",
        catalog_rows,
        [
            "panel_species", "assembly_accession", "organism_name", "assembly_level",
            "refseq_category", "version_status", "seq_rel_date", "ftp_path", "source_catalog",
        ],
    )

    if args.catalog_only:
        print(f"catalog_only\tassemblies={len(selected_rows)}\tsource={args.source}")
        return 0

    taxonomy_cache = fetch_taxonomy_lineages(
        sorted({row.get("taxid", "") for row in selected_rows if row.get("taxid")}),
        cache_path=cache_base / "taxonomy_cache.json",
    )

    # Download BioSample phenotypic metadata once into data_cache (reused across runs).
    phenotype_cache_path = cache_base / "phenotypic_metadata.json"
    biosample_ids = sorted({row.get("biosample", "") for row in selected_rows if row.get("biosample")})
    if biosample_ids:
        phenotype_meta = fetch_ncbi_biosample_phenotypes(biosample_ids, cache_path=phenotype_cache_path)
        print(f"[phenotype] cached {len(phenotype_meta)} BioSample records -> {phenotype_cache_path}")
    else:
        phenotype_meta: dict[str, dict[str, str]] = {}

    # One-shot caches for downstream multi-omics signal. Both live under
    # data_cache/ so subsequent panel preparations reuse the JSON / CSV
    # without refetching from EBI (Atlas) or GitHub (FungalTraits). The
    # per-panel `expression.tsv` and `ecological_traits.tsv` are stitched
    # from these caches at the end of prepare, after query_rows are known.
    expression_atlas_cache = cache_base / "expression_atlas"
    expression_analytics_cache = cache_base / "expression_atlas_analytics"
    panel_species = sorted({sel.get("species", "") for sel in selectors if sel.get("species")})
    expression_listings: dict[str, Path] = {}
    if panel_species:
        expression_listings = fetch_expression_atlas_for_species(
            panel_species, expression_atlas_cache,
        )
        if expression_listings:
            print(
                f"[expression-atlas] cached listings for "
                f"{len(expression_listings)}/{len(panel_species)} species "
                f"-> {expression_atlas_cache}"
            )
    fungaltraits_csv = fetch_fungaltraits_table(cache_base)
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
    # Per-ref (asm_name, gff.gz path) pairs for the GFF -> gene_annotations.tsv
    # converter. Populated when either a GFF or GBFF annotation source succeeds
    # for a given assembly.
    gff_pairs: list[tuple[str, Path]] = []

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
            for url, filename in ncbi_download_targets(row, include_gff=False):
                local = materialize_entry(url, refs_dir / filename, keep_gz=True)
                if filename.endswith("_genomic.fna.gz"):
                    downloaded_fasta = local
            if downloaded_fasta is None:
                continue
            ref_downloads += 1
            if args.download_gff:
                annotation_source = download_ncbi_gene_annotation_source(
                    row, asm_name, refs_dir, "",
                    ensembl_cache_dir=cache_base,
                )
                if annotation_source is not None:
                    gff_pairs.append((asm_name, annotation_source))
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
        defaults_entry = {
            "benchmark_ref_asm": (
                benchmark_ref_row.get("assembly_accession", "").replace(".", "_")
                if benchmark_ref_local is not None and benchmark_ref_row.get("assembly_accession")
                else species_ref_asm
            ),
            "benchmark_ref_fasta": str(benchmark_ref_local) if benchmark_ref_local else ".",
        }
        # Index under both the NCBI taxonomy species name and the panel-preset
        # species name. The reads-mode loop below looks up by the panel-preset
        # name; without the alias, NCBI strain suffixes (e.g. "Lachancea
        # kluyveri NRRL Y-12651") silently miss and read queries are dropped.
        species_benchmark_defaults[species_name] = defaults_entry
        if species and species != species_name:
            species_benchmark_defaults.setdefault(species, defaults_entry)

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
            if args.download_gff:
                annotation_source = download_ncbi_gene_annotation_source(
                    row, asm_name, queries_dir, "query ",
                    ensembl_cache_dir=cache_base,
                )
                if annotation_source is not None:
                    gff_pairs.append((asm_name, annotation_source))
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
                "source": row.get("_catalog_source", args.source),
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
        sys.stderr.write(
            f"[reads-mode] resolving ENA runs for {len(selectors)} panel species "
            f"(modes={requested_read_modes}, max_runs/species={args.read_accessions_per_species})\n"
        )
        for sel in selectors:
            species = sel.get("species", "")
            if not species:
                continue
            default_benchmark = species_benchmark_defaults.get(species, {})
            if not default_benchmark.get("benchmark_ref_fasta"):
                # No reference was downloaded for this species — reads-mode
                # benchmarks require one, so skip.
                sys.stderr.write(
                    f"[reads-mode] skip {species!r}: no reference assembly available "
                    f"(known species keys: {sorted(species_benchmark_defaults.keys())[:5]}...)\n"
                )
                continue
            ena_runs = fetch_ena_read_runs_by_species(species, max_rows=args.ena_max_rows_per_species)
            if not ena_runs:
                sys.stderr.write(
                    f"[reads-mode] skip {species!r}: ENA filereport returned 0 runs "
                    f"(check network or species name)\n"
                )
                continue
            sys.stderr.write(f"[reads-mode] {species!r}: ENA returned {len(ena_runs)} candidate runs\n")
            for read_mode in requested_read_modes:
                # Pull up to 4x the requested count so we can retry past runs
                # whose payloads fail FASTQ validation (e.g. PacBio bas.h5
                # accessions whose fastq_ftp is empty drop out in
                # select_ena_read_sources, but the fastq_ftp value can still
                # point to mirrors that occasionally serve corrupt content).
                pool_size = max(args.read_accessions_per_species * 4, args.read_accessions_per_species)
                pool_urls, pool_meta = select_ena_read_sources(
                    ena_runs, read_mode, pool_size
                )
                if not pool_urls:
                    sys.stderr.write(
                        f"[reads-mode] {species!r} mode={read_mode}: no eligible runs after platform filter\n"
                    )
                    continue
                local_path = None
                meta_rows: list[dict[str, str]] = []
                # Each meta row carries the complete URL group for one run
                # (including paired-end mates); walk run-by-run until one
                # downloads and validates as FASTQ.
                attempts = 0
                for meta in pool_meta:
                    if attempts >= args.read_accessions_per_species:
                        break
                    run_acc = meta["run_accession"]
                    urls = selected_urls_from_ena_meta(meta)
                    if not urls:
                        continue
                    asm_name = normalize_name(f"{species}_{read_mode}_{run_acc}")
                    try:
                        candidate_path = merge_sequence_sources(urls, queries_dir / asm_name)
                    except Exception as exc:
                        sys.stderr.write(
                            f"[warn] {species}: ENA {read_mode} download failed for {run_acc}: {exc}\n"
                        )
                        continue
                    local_path = candidate_path
                    meta_rows = [meta]
                    attempts += 1
                    break  # one validated run is enough for this (species, mode)
                if local_path is None:
                    sys.stderr.write(
                        f"[warn] {species}: no usable ENA {read_mode} run after validation; skipping\n"
                    )
                    continue
                sys.stderr.write(
                    f"[reads-mode] {species!r} mode={read_mode}: picked {len(meta_rows)} run(s)\n"
                )
                asm_name = normalize_name(f"{species}_{read_mode}_{meta_rows[0]['run_accession']}")
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
            "instrument_platform", "library_layout", "run_accession",
        ],
    )
    write_tsv(out_dir / "benchmark_reference_map.tsv", benchmark_map_rows, ["query_asm", "benchmark_ref_asm", "benchmark_ref_fasta", "species"])
    write_tsv(
        out_dir / "public_data_links.tsv",
        source_link_rows,
        ["query_asm", "role", "query_mode", "source_type", "source_accession", "source_url", "local_path", "species"],
    )
    write_public_resource_links(out_dir / "public_resource_links.tsv")

    # Write a symlink / copy of the phenotypic metadata into the prepared dir for
    # downstream analysis scripts that expect it alongside other manifests.
    if phenotype_meta:
        phenotype_out = out_dir / "phenotypic_metadata.json"
        if not phenotype_out.exists():
            phenotype_out.write_text(
                json.dumps(phenotype_meta, indent=2, sort_keys=True), encoding="utf-8"
            )

    # Auto-build prepared_dir/gene_annotations.tsv from any GFF/GBFF annotation
    # sources we downloaded alongside FASTA. The benchmark step will pick this
    # up automatically (no need for the caller to pass --gene-annotations-tsv).
    #
    # Indexing trick: the SV biology analyzer keys gene lookups by
    # (query_asm, contig), but our GFFs come from refs. Ref-coordinate calls
    # in calls.hits.tsv carry query_asm=<query asm name> and contig=<ref
    # contig>, so we duplicate every gene row under every query_asm that
    # uses that ref (per benchmark_map_rows). This makes (query_asm, ref_contig)
    # lookups hit without forcing the analyzer to learn ref↔query joins.
    gene_annotations_count = 0
    if gff_pairs:
        gene_annotations_path = out_dir / "gene_annotations.tsv"
        # Build (asm_name -> [aliases...]) so each gene row is emitted under
        # every form the analyzer might lookup against. Aliases include:
        #   1. the prepared manifest's asm_name itself (GCA_000146045_2)
        #   2. the FASTA basename without .gz (GCA_000146045.2_R64_genomic.fna)
        #      because calls.hits.tsv writes query_asm as the FASTA filename
        #   3. every query_asm that uses this ref as benchmark_ref_asm
        # Without all three, ref-coord candidates from query X won't find
        # genes annotated against ref Y, and the nearest_gene fallback is silent.
        asm_aliases: dict[str, set[str]] = defaultdict(set)
        for ref_row in ref_manifest_rows:
            asm = ref_row.get("asm_name", "")
            if not asm:
                continue
            asm_aliases[asm].add(asm)
            fasta_path = ref_row.get("fasta_path", "")
            if fasta_path:
                basename = Path(fasta_path).name
                if basename.endswith(".gz"):
                    basename = basename[:-3]
                asm_aliases[asm].add(basename)
        for q_row in query_rows:
            q_asm = q_row.get("query_asm", "")
            q_path = q_row.get("path", "")
            if q_asm and q_path:
                basename = Path(q_path).name
                if basename.endswith(".gz"):
                    basename = basename[:-3]
                asm_aliases[q_asm].add(q_asm)
                asm_aliases[q_asm].add(basename)
        ref_to_queries: dict[str, list[str]] = defaultdict(list)
        for bm in benchmark_map_rows:
            ref_asm = bm.get("benchmark_ref_asm", "")
            q_asm = bm.get("query_asm", "")
            if ref_asm and q_asm:
                ref_to_queries[ref_asm].append(q_asm)
        gene_annotations_count = stream_gene_annotations_to_tsv(
            gene_annotations_path, gff_pairs, asm_aliases, ref_to_queries,
        )
        if gene_annotations_count:
            sys.stderr.write(
                f"[gene-annot] wrote {gene_annotations_count} gene records "
                f"(ref-keyed + per-query / FASTA-basename aliases) to {gene_annotations_path}\n"
            )
        else:
            try:
                gene_annotations_path.unlink()
            except FileNotFoundError:
                pass

    # Stitch the cached Expression Atlas / FungalTraits data into the
    # per-panel TSVs the analyzer auto-picks-up at benchmark time. Building
    # them here (rather than re-downloading on every benchmark) means each
    # benchmark run reuses the cache without any network round-trip.
    species_to_query_asms: dict[str, list[str]] = defaultdict(list)
    for q_row in query_rows:
        species = q_row.get("species") or "."
        species_to_query_asms[species].append(q_row["query_asm"])

    expression_rows_written = 0
    if expression_listings and species_to_query_asms:
        expression_rows_written = assemble_expression_tsv(
            expression_listings,
            species_to_query_asms,
            expression_analytics_cache,
            out_dir / "expression.tsv",
        )
        if expression_rows_written:
            sys.stderr.write(
                f"[expression-atlas] stitched {expression_rows_written} "
                f"gene-level rows into {out_dir / 'expression.tsv'}\n"
            )

    ecological_rows_written = 0
    if species_to_query_asms:
        ecological_rows_written = write_ecological_summary_tsv(
            fungaltraits_csv,
            species_to_query_asms,
            out_dir / "ecological_traits.tsv",
        )
        if ecological_rows_written:
            sys.stderr.write(
                f"[fungaltraits] joined {ecological_rows_written} "
                f"trait records into {out_dir / 'ecological_traits.tsv'}\n"
            )

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
                "data_cache_dir": str(cache_base),
                "assembly_summary_cache": ";".join(str(p) for p in assembly_summary_caches),
                "assembly_summary_caches": [str(p) for p in assembly_summary_caches],
                "phenotypic_records_cached": len(phenotype_meta),
                "expression_rows_written": expression_rows_written,
                "ecological_rows_written": ecological_rows_written,
                "expression_atlas_species_cached": len(expression_listings),
                "fungaltraits_csv_cached": bool(fungaltraits_csv and fungaltraits_csv.exists()),
            },
            fh,
            indent=2,
            sort_keys=True,
        )
    print(f"prepared\trefs={len(ref_manifest_rows)}\tqueries={len(query_rows)}\tsource={args.source}")
    return 0


def load_query_manifest(path: Path) -> list[dict[str, str]]:
    # Surface a clear, actionable error when the manifest is missing — the raw
    # FileNotFoundError used to bubble up from `path.open()` deep in benchmark()
    # with no hint that the upstream `prepare` step never finished (e.g. another
    # job partially cleaned the prepared dir, or the user pointed --prepared-dir
    # at the wrong place).
    if not path.exists():
        raise FileNotFoundError(
            f"query_manifest.tsv not found at {path}. The prepared directory is "
            f"incomplete or was never populated by the `prepare` sub-command. "
            f"Re-run prepare against {path.parent} before benchmarking."
        )
    rows: list[dict[str, str]] = []
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows.append({k: (v or "").strip() for k, v in row.items()})
    return rows


def _normalise_query_group_token(token: str) -> str:
    t = re.sub(r"[^a-z0-9]+", "", (token or "").lower())
    aliases = {
        "penciliium": "penicillium",
        "penicilium": "penicillium",
        "mychrrzia": "mycorrhiza",
        "mychrrhiza": "mycorrhiza",
        "mycorrhzia": "mycorrhiza",
        "trichderma": "trichoderma",
        "trichdermia": "trichoderma",
    }
    return aliases.get(t, t)


def _query_row_taxon_blob(row: dict[str, str]) -> str:
    fields = (
        "phylum", "class", "order", "family", "genus", "species",
        "query_asm", "path", "benchmark_ref_asm", "benchmark_ref_fasta",
        "lifestyle", "architecture",
    )
    return " ".join((row.get(f) or "") for f in fields).lower()


def select_one_query_per_group(
    query_manifest: list[dict[str, str]],
    group_spec: str,
    out_dir: Path,
    hierarchy_manifest: Path | None = None,
) -> list[dict[str, str]]:
    requested = [
        _normalise_query_group_token(tok)
        for tok in re.split(r"[,;\s]+", group_spec or "")
        if tok.strip()
    ]
    if not requested:
        return query_manifest

    def matches(row: dict[str, str], target: str) -> bool:
        blob = _query_row_taxon_blob(row)
        compact = re.sub(r"[^a-z0-9]+", "", blob)
        if target == "mycorrhiza":
            return any(term in compact for term in (
                "mycorrhiza", "mycorrhizal", "rhizophagus", "glomus",
                "glomeromycotina", "glomeromycetes", "arbuscular"
            ))
        return target in compact

    def hierarchy_rows() -> list[dict[str, str]]:
        if hierarchy_manifest is None or not hierarchy_manifest.exists():
            return []
        rows: list[dict[str, str]] = []
        with hierarchy_manifest.open() as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                cols = line.split("\t")
                if len(cols) < 9:
                    continue
                rows.append({
                    "asm": cols[0], "phylum": cols[1], "class": cols[2],
                    "order": cols[3], "family": cols[4], "genus": cols[5],
                    "species": cols[6], "rank": cols[7], "path": cols[8],
                })
        return rows

    hrows = hierarchy_rows()

    def synthesize_from_hierarchy(target: str) -> dict[str, str] | None:
        candidates = [r for r in hrows if r.get("rank") == "species" and matches(r, target)]
        if not candidates:
            return None
        q = candidates[0]
        ref = next(
            (r for r in candidates[1:] if r.get("asm") != q.get("asm")),
            None,
        )
        if ref is None:
            ref = next(
                (r for r in hrows
                 if r.get("rank") == "species"
                 and r.get("asm") != q.get("asm")
                 and r.get("genus") == q.get("genus")),
                None,
            )
        if ref is None:
            ref = next(
                (r for r in hrows
                 if r.get("rank") == "species"
                 and r.get("asm") != q.get("asm")
                 and r.get("family") == q.get("family")),
                None,
            )
        ref = ref or q
        return {
            "query_asm": q.get("asm", ""),
            "query_mode": "assembly",
            "path": q.get("path", ""),
            "scenario": "million_real_group_target",
            "lifestyle": ".",
            "architecture": ".",
            "benchmark_ref_asm": ref.get("asm", "."),
            "benchmark_ref_fasta": ref.get("path", "."),
            "phylum": q.get("phylum", "."),
            "class": q.get("class", "."),
            "order": q.get("order", "."),
            "family": q.get("family", "."),
            "genus": q.get("genus", "."),
            "species": q.get("species", "."),
            "source": "hierarchy_manifest_group_target",
            "instrument_platform": ".",
            "library_layout": ".",
            "run_accession": ".",
        }

    selected: list[dict[str, str]] = []
    selected_ids: set[str] = set()
    report_rows: list[dict[str, str]] = []
    for target in requested:
        row = next((r for r in query_manifest if matches(r, target)), None)
        status = "found_query_manifest"
        if row is None:
            row = synthesize_from_hierarchy(target)
            status = "synthesized_from_hierarchy"
        if row is None:
            report_rows.append({
                "requested_group": target, "status": "missing",
                "query_asm": ".", "genus": ".", "species": ".", "path": ".",
            })
            continue
        qid = row.get("query_asm") or row.get("path") or repr(row)
        if qid not in selected_ids:
            selected_ids.add(qid)
            selected.append(row)
        report_rows.append({
            "requested_group": target, "status": status,
            "query_asm": row.get("query_asm", "."),
            "genus": row.get("genus", "."),
            "species": row.get("species", "."),
            "path": row.get("path", "."),
        })

    report = out_dir / "REQUESTED_QUERY_GROUPS.tsv"
    write_tsv(report, report_rows, [
        "requested_group", "status", "query_asm", "genus", "species", "path",
    ])
    missing = [r["requested_group"] for r in report_rows if r["status"] == "missing"]
    if missing:
        sys.stderr.write(
            "[benchmark] requested query group(s) absent from query_manifest.tsv "
            f"and hierarchy_manifest.tsv: {', '.join(missing)}\n"
        )
    if selected:
        sys.stderr.write(
            "[benchmark] selected one query per requested group where available "
            "(query_manifest first, hierarchy_manifest fallback): "
            f"{len(selected)} of {len(requested)} group(s); report={report}\n"
        )
        return selected
    sys.stderr.write(
        "[benchmark] none of the requested query groups were present; keeping "
        "the mode-matched query manifest unchanged.\n"
    )
    return query_manifest


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


_BND_ALT_RE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


def _parse_bnd_alt_mate(alt: str) -> tuple[str, int] | None:
    m = _BND_ALT_RE.search(alt or "")
    if not m:
        return None
    try:
        return m.group(1), int(m.group(2))
    except ValueError:
        return None


def _mate_fields(info: dict[str, str], alt: str = "") -> tuple[str, int, int]:
    mate_contig = info.get("CHR2") or info.get("MATE_CONTIG") or "."
    mate_pos = _parse_int_field(info.get("POS2") or info.get("END2"))
    if (not mate_contig or mate_contig == "." or mate_pos is None) and ("[" in alt or "]" in alt):
        parsed = _parse_bnd_alt_mate(alt)
        if parsed is not None:
            mate_contig, mate_pos = parsed
    mate_end = _parse_int_field(info.get("END2"))
    if mate_pos is None:
        mate_pos = 0
    if mate_end is None:
        mate_end = mate_pos
    return mate_contig or ".", mate_pos, mate_end


def qasm_matches_observed(expected: str, observed: str) -> bool:
    expected_norm = normalize_name(expected)
    observed_norm = normalize_name(observed)
    if expected_norm == observed_norm:
        return True
    return observed_norm.startswith(expected_norm + "_") or expected_norm.startswith(observed_norm + "_")


def _parse_int_field(value: str | None) -> int | None:
    if value is None or value == "":
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def _parse_pseudocontig_support_name(name: str) -> int | None:
    for tag in ("_n", "_mf"):
        pos = name.rfind(tag)
        if pos < 0:
            continue
        start = pos + len(tag)
        end = start
        while end < len(name) and name[end].isdigit():
            end += 1
        if end == start:
            continue
        try:
            return int(name[start:end])
        except ValueError:
            continue
    return None


def _mycosv_intrinsic_read_support(info: dict[str, str], query_contig: str) -> int | None:
    support = _parse_int_field(info.get("SUPPORT"))
    if support is not None:
        return support
    return _parse_pseudocontig_support_name(query_contig)


def load_mycosv_query_calls(vcf_path: Path, query_asm: str) -> list[NormalizedCall]:
    rows: list[NormalizedCall] = []
    if not vcf_path.exists():
        return rows
    with open_text_auto(vcf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            info = parse_info_field(fields[7])
            if not qasm_matches_observed(query_asm, info.get("QASM", "")):
                continue
            svtype = TYPE_CANON.get(info.get("SVTYPE", fields[4].strip("<>")))
            if not svtype:
                continue
            pos = int(fields[1])
            end = int(info.get("END", pos))
            svlen = int(info.get("SVLEN", end - pos + 1))
            mate_contig, mate_pos, mate_end = _mate_fields(info, fields[4])
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
                read_support=_mycosv_intrinsic_read_support(info, fields[0]),
                mate_contig=mate_contig,
                mate_pos=mate_pos,
                mate_end=mate_end,
            ))
    return rows


def fasta_bp_and_contigs(path: Path | None) -> tuple[int, int]:
    if path is None or not path.exists():
        return 0, 0
    total_bp = 0
    contigs = 0
    with open_text_auto(path) as fh:
        for line in fh:
            if line.startswith(">"):
                contigs += 1
            else:
                total_bp += len(line.strip())
    return total_bp, contigs


def sequence_bp_and_records(path: Path | None) -> tuple[int, int]:
    if path is None or not path.exists():
        return 0, 0
    if sequence_kind_from_name(path.name) == "fastq":
        total_bp = 0
        records = 0
        with open_text_auto(path) as fh:
            while True:
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline()
                plus = fh.readline()
                qual = fh.readline()
                if not qual:
                    break
                total_bp += len(seq.strip())
                records += 1
        return total_bp, records
    return fasta_bp_and_contigs(path)


def mycosv_qasm_counts(vcf_path: Path) -> tuple[dict[str, int], int, int]:
    counts: dict[str, int] = defaultdict(int)
    missing_qasm = 0
    malformed = 0
    if not vcf_path.exists():
        return counts, missing_qasm, malformed
    with open_text_auto(vcf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                malformed += 1
                continue
            info = parse_info_field(fields[7])
            qasm = info.get("QASM", "")
            if not qasm:
                missing_qasm += 1
                qasm = "."
            counts[qasm] += 1
    return counts, missing_qasm, malformed


def write_sv_volume_audit(
    out_path: Path,
    query_manifest: list[dict[str, str]],
    mycosv_calls_by_query: dict[str, dict[str, Any]],
    truth_sets: dict[str, dict[tuple[str, str], list[NormalizedCall]]],
    mycosv_vcf: Path,
    *,
    mode: str,
    mycosv_failed: bool,
) -> list[dict[str, object]]:
    qasm_counts, missing_qasm, malformed = mycosv_qasm_counts(mycosv_vcf)
    rows: list[dict[str, object]] = []
    for row in query_manifest:
        query_asm = row["query_asm"]
        query_path = Path(row.get("mycosv_path") or row["path"])
        query_bp, query_contigs = sequence_bp_and_records(query_path)
        bench_ref = (row.get("benchmark_ref_fasta") or "").strip()
        ref_bp, ref_contigs = fasta_bp_and_contigs(
            Path(bench_ref) if bench_ref and bench_ref != "." else None
        )
        estimated_cov = (float(query_bp) / float(ref_bp)) if ref_bp > 0 else 0.0
        mycosv_count = len(mycosv_calls_by_query.get(query_asm, {}).get("query", []))
        truth_counts = {
            f"{coord}:{label}": len(calls)
            for (coord, label), calls in truth_sets.get(query_asm, {}).items()
        }
        max_truth = max(truth_counts.values(), default=0)
        mycosv_query_calls = list(mycosv_calls_by_query.get(query_asm, {}).get("query", []))
        n_off_ref = sum(1 for c in mycosv_query_calls if getattr(c, "svtype", "") == "OFF_REF")
        n_anchored = len(mycosv_query_calls) - n_off_ref
        qmode = (row.get("query_mode") or mode or "").lower()
        scenario = (row.get("scenario") or "").lower()
        species = (row.get("species") or "").lower()
        is_read_capped = (
            qmode in {"short-reads", "long-reads"}
            and "read_subsets" in Path(row.get("path", "")).parts
        )
        if mycosv_failed:
            status = "fail"
            diagnosis = "mycosv_failed"
        elif malformed or missing_qasm:
            status = "fail"
            diagnosis = "malformed_or_missing_qasm_records"
        elif mycosv_count == 0:
            status = "fail"
            diagnosis = "no_calls_for_query"
        elif max_truth > 0 and n_anchored == 0 and n_off_ref > 0:
            status = "fail"
            diagnosis = "off_ref_only_against_anchored_truth"
        elif qmode == "long-reads" and 0.0 < estimated_cov < 5.0:
            status = "fail"
            diagnosis = "long_read_input_too_low_coverage_for_sv_volume"
        elif qmode == "short-reads" and 0.0 < estimated_cov < 10.0:
            status = "fail"
            diagnosis = "short_read_input_too_low_coverage_for_sv_volume"
        elif qmode in {"short-reads", "long-reads"} and is_read_capped:
            status = "ok"
            diagnosis = "bounded_benchmark_read_input"
        elif max_truth > 0 and mycosv_count < max_truth:
            status = "low"
            diagnosis = "below_comparator_burden"
        elif max_truth == 0:
            status = "needs_model"
            diagnosis = "no_truth_or_comparator_volume_model"
        else:
            status = "ok"
            diagnosis = "consistent_with_available_comparator_burden"
        rows.append({
            "query_asm": query_asm,
            "query_mode": qmode or ".",
            "scenario": row.get("scenario", "."),
            "species": row.get("species", "."),
            "query_bp": query_bp,
            "query_contigs": query_contigs,
            "benchmark_ref_bp": ref_bp,
            "benchmark_ref_contigs": ref_contigs,
            "estimated_query_coverage": f"{estimated_cov:.2f}" if estimated_cov > 0 else "0",
            "mycosv_query_calls": mycosv_count,
            "max_comparator_or_truth_calls": max_truth,
            "observed_qasm_count_in_vcf": sum(
                count for qasm, count in qasm_counts.items()
                if qasm_matches_observed(query_asm, qasm)
            ),
            "missing_qasm_records_in_vcf": missing_qasm,
            "malformed_records_in_vcf": malformed,
            "status": status,
            "diagnosis": diagnosis,
            "truth_counts": json.dumps(truth_counts, sort_keys=True),
        })
    write_tsv(
        out_path,
        rows,
        [
            "query_asm", "query_mode", "scenario", "species",
            "query_bp", "query_contigs", "benchmark_ref_bp", "benchmark_ref_contigs",
            "estimated_query_coverage",
            "mycosv_query_calls", "max_comparator_or_truth_calls",
            "observed_qasm_count_in_vcf", "missing_qasm_records_in_vcf",
            "malformed_records_in_vcf", "status", "diagnosis", "truth_counts",
        ],
    )
    return rows


def load_mycosv_paired_calls(
    vcf_path: Path, query_asm: str,
) -> tuple[list[NormalizedCall], list[NormalizedCall], list[tuple[str, str, int, int, str]]]:
    """Single-pass loader that returns (query_calls, ref_calls, ref_to_query_keys).

    `ref_to_query_keys[i]` is the call_key tuple of the QUERY-coord
    NormalizedCall produced from the same VCF row that yielded `ref_calls[i]`.
    Lets the benchmark loop, after matching reference-coord truth vs
    reference-coord mycosv predictions, attribute the comparator-support label
    back to the QUERY-coord call_key that `support_by_key` is indexed by — the
    bridge across coord spaces that was previously missing and caused every
    mycosv call in biology_findings.tsv to be flagged `mycosv_unique=yes`.
    """
    query_calls: list[NormalizedCall] = []
    ref_calls: list[NormalizedCall] = []
    ref_to_query_keys: list[tuple[str, str, int, int, str]] = []
    if not vcf_path.exists():
        return query_calls, ref_calls, ref_to_query_keys
    with open_text_auto(vcf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            info = parse_info_field(fields[7])
            if not qasm_matches_observed(query_asm, info.get("QASM", "")):
                continue
            svtype = TYPE_CANON.get(info.get("SVTYPE", fields[4].strip("<>")))
            if not svtype:
                continue
            try:
                query_pos = int(fields[1])
            except ValueError:
                continue
            query_end = int(info.get("END", query_pos))
            svlen = int(info.get("SVLEN", query_end - query_pos + 1))
            mate_contig, mate_pos, mate_end = _mate_fields(info, fields[4])
            query_contig = fields[0]
            ref_contig = info.get("REFCONTIG", ".")
            ref_pos_raw = info.get("REFPOS", "")
            q_call = NormalizedCall(
                query_asm=query_asm,
                query_contig=query_contig,
                pos=query_pos,
                end=query_end,
                svtype=svtype,
                svlen=svlen,
                source="mycosv",
                coord_space="query",
                annotation=info.get("ANNOT", "."),
                element_class=info.get("EC", "NONE"),
                ref_asm=info.get("CLADE", "."),
                ref_contig=ref_contig,
                read_support=_mycosv_intrinsic_read_support(info, fields[0]),
                mate_contig=mate_contig,
                mate_pos=mate_pos,
                mate_end=mate_end,
            )
            query_calls.append(q_call)
            if ref_contig not in {"", "."} and ref_pos_raw:
                try:
                    ref_pos = int(ref_pos_raw)
                except ValueError:
                    continue
                ref_end = int(info.get("REFEND", ref_pos))
                ref_calls.append(NormalizedCall(
                    query_asm=query_asm,
                    query_contig=query_contig,
                    pos=ref_pos,
                    end=ref_end,
                    svtype=svtype,
                    svlen=svlen,
                    source="mycosv",
                    coord_space="reference",
                    annotation=info.get("ANNOT", "."),
                    element_class=info.get("EC", "NONE"),
                    ref_asm=info.get("CLADE", "."),
                    ref_contig=ref_contig,
                    read_support=_mycosv_intrinsic_read_support(info, fields[0]),
                    mate_contig=mate_contig,
                    mate_pos=mate_pos,
                    mate_end=mate_end,
                ))
                ref_to_query_keys.append(
                    (query_asm, query_contig, query_pos, query_end, svtype)
                )
    return query_calls, ref_calls, ref_to_query_keys


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
            mate_contig, mate_pos, mate_end = _mate_fields(info, fields[4])
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
                read_support=_mycosv_intrinsic_read_support(info, fields[0]),
                mate_contig=mate_contig,
                mate_pos=mate_pos,
                mate_end=mate_end,
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
                mate_contig=row.get("mate_contig") or row.get("chr2") or ".",
                mate_pos=int(float(row.get("mate_pos") or row.get("pos2") or 0)),
                mate_end=int(float(
                    row.get("mate_end") or row.get("end2")
                    or row.get("mate_pos") or row.get("pos2") or 0
                )),
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


_REF_VCF_MIN_SV_BP = 30  # paftools/anchorwave emit single-bp variants too; skip <30 bp.


def _infer_svtype_from_alleles(ref_allele: str, alt_allele: str) -> tuple[str | None, int]:
    """Derive (svtype, svlen) from REF/ALT when SVTYPE INFO is absent.

    Used for VCFs emitted by paftools.js (AnchorWave pipeline) and svim-asm
    in non-symbolic mode, which represent INS/DEL with explicit allele
    sequences rather than `<INS>` / `<DEL>` symbolic ALT.
    """
    if not ref_allele or not alt_allele or alt_allele in {".", "*"}:
        return None, 0
    # BND notation: ALT contains breakend brackets.
    if "[" in alt_allele or "]" in alt_allele:
        return "TRA", 0
    # Symbolic ALT (e.g. <INS>, <DEL>) — already handled by caller; bail.
    if alt_allele.startswith("<"):
        return None, 0
    diff = len(alt_allele) - len(ref_allele)
    if diff > 0:
        return "INS", diff
    if diff < 0:
        return "DEL", -diff
    return None, 0  # SNV or MNV — not an SV.


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
            ref_allele = fields[3]
            alt_allele = fields[4]
            svtype = TYPE_CANON.get((info.get("SVTYPE") or alt_allele.strip("<>")).upper())
            inferred_svlen: int | None = None
            if not svtype:
                # paftools.js / svim-asm non-symbolic VCFs do not set SVTYPE
                # and use explicit allele sequences. Derive from REF/ALT so
                # those callers contribute to the comparator baseline set.
                svtype, inferred_svlen = _infer_svtype_from_alleles(ref_allele, alt_allele)
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
                svlen = abs(int(svlen_raw)) if svlen_raw else 0
            except ValueError:
                svlen = 0
            if svlen == 0 and inferred_svlen is not None:
                svlen = inferred_svlen
            if svlen == 0:
                svlen = max(1, end - pos + 1)
            # For DEL inferred from alleles, set END = pos + svlen so that
            # downstream reference-coord matching has a meaningful interval.
            if inferred_svlen is not None and end_raw == "":
                if svtype == "DEL":
                    end = pos + svlen
                else:
                    end = pos
            # Skip sub-SV-size events (paftools emits SNV/MNV-like rows too).
            if svtype in {"INS", "DEL", "DUP", "INV"} and svlen < _REF_VCF_MIN_SV_BP:
                continue
            mate_contig, mate_pos, mate_end = _mate_fields(info, alt_allele)
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
                mate_contig=mate_contig,
                mate_pos=mate_pos,
                mate_end=mate_end,
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
    """Load minigraph bubble calls as a normalized comparator callset.

    minigraph emits one bubble per local divergence, including microsatellite
    expansions and 30 bp indels — typical yeast assemblies produce ~1 000 bubbles
    per sample, of which 60–80 % are sub-50 bp events that mycosv's chain caller
    intentionally collapses. Comparing 1 000 minigraph bubbles against
    ~20 mycosv chain-level events caps recall at ~2 % independent of correctness.

    To bring the truth-set granularity in line with mycosv (and with how human
    SV benchmarks usually score):
      1. optionally drop bubbles below MINIGRAPH_MIN_SV_BP (default 0 keeps
         minigraph's raw calls; set 50 to match the usual SV-size convention);
      2. coalesce adjacent same-type bubbles within MINIGRAPH_MERGE_GAP_BP
         (default 1 000 bp) on the same contig into a single representative event
         whose pos = first bubble start, end = last bubble end, svlen = sum.

    Both thresholds are overridable via env vars so a panel that wants the raw
    fine-grained bubbles back can opt out (set MINIGRAPH_MIN_SV_BP=0 and
    MINIGRAPH_MERGE_GAP_BP=0).
    """
    rows: list[NormalizedCall] = []
    if not bubble_bed.exists() or not sample_bed.exists():
        return rows
    bubble_lines = [line.rstrip("\n") for line in bubble_bed.read_text(encoding="utf-8").splitlines() if line.strip()]
    sample_lines = [line.rstrip("\n") for line in sample_bed.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(bubble_lines) != len(sample_lines):
        return rows
    try:
        min_sv_bp = int(os.environ.get("MINIGRAPH_MIN_SV_BP", "0"))
    except ValueError:
        min_sv_bp = 0
    try:
        merge_gap_bp = int(os.environ.get("MINIGRAPH_MERGE_GAP_BP", "1000"))
    except ValueError:
        merge_gap_bp = 1000

    raw: list[NormalizedCall] = []
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
        if min_sv_bp > 0 and svlen < min_sv_bp:
            continue
        raw.append(NormalizedCall(
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
    if merge_gap_bp <= 0 or not raw:
        return raw
    # Coalesce adjacent same-type bubbles on the same contig within merge_gap_bp.
    raw.sort(key=lambda c: (c.ref_contig, c.svtype, c.pos))
    for call in raw:
        if not rows:
            rows.append(call)
            continue
        last = rows[-1]
        same_block = (
            last.ref_contig == call.ref_contig
            and last.svtype == call.svtype
            and (call.pos - last.end) <= merge_gap_bp
        )
        if not same_block:
            rows.append(call)
            continue
        rows[-1] = replace(
            last,
            end=max(last.end, call.end),
            svlen=last.svlen + call.svlen,
        )
    return rows


def canonical_group(svtype: str) -> str:
    if svtype in {"TRANS", "INVTR", "TRA"}:
        return "TRA"
    if svtype in {"DUP", "INVDP"}:
        return "DUP"
    return svtype


def effective_contig(call: NormalizedCall) -> str:
    return call.ref_contig if call.coord_space == "reference" else call.query_contig


def _call_span_end(call: NormalizedCall) -> int:
    return max(call.end, call.pos + abs(call.svlen))


def _call_span_contains(span_call: NormalizedCall, pos: int, *, pad: int = 0) -> bool:
    return (span_call.pos - pad) <= pos <= (_call_span_end(span_call) + pad)


def _has_mate(call: NormalizedCall) -> bool:
    return call.mate_contig not in {"", "."} and call.mate_pos > 0


_SPAN_CONTAIN_TYPES = {"INV", "TRA", "DEL", "DUP"}


def _span_contain_applies(truth: NormalizedCall, pred: NormalizedCall) -> bool:
    """Return True when the predicted call's span legitimately covers the truth
    breakpoint AND the type group supports span-based matching.

    MycoSV emits chain-level DEL/DUP blocks that span the genomic interval
    between consecutive MEM anchors; read-level callers (svim, sniffles,
    cutesv, delly, manta) emit per-event fine-grained breakpoints. Allowing
    span-containment lets one chain-level pred match a single fine-grained
    truth event nested inside it. The truth-loop in match_calls() still
    consumes the pred greedily, so a single coarse pred is credited for at
    most one TP — it cannot inflate TP by claiming many fine-grained truths
    at once.
    """
    group = canonical_group(truth.svtype)
    if group not in _SPAN_CONTAIN_TYPES and canonical_group(pred.svtype) not in _SPAN_CONTAIN_TYPES:
        return False
    if truth.svtype != "INV" and pred.svtype != "INV" \
            and canonical_group(truth.svtype) != "TRA" \
            and canonical_group(pred.svtype) != "TRA":
        # DEL/DUP: only allow span-contain when pred is at least 2x the truth
        # length. Otherwise the regular pos+length tolerance is the right gate
        # and span-contain would silently relax length agreement on co-located
        # but length-disagreeing same-scale calls.
        if abs(pred.svlen) < 2 * max(1, abs(truth.svlen)):
            return False
    tol_bp = DEFAULT_TOL_BP.get(group, 500)
    return _call_span_contains(pred, truth.pos, pad=tol_bp)


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
    span_match = _span_contain_applies(truth, pred)
    pos_within_tol = abs(truth.pos - pred.pos) <= tol_bp
    if not (pos_within_tol or span_match):
        return False
    if truth.svtype not in {"INV", "TRA", "OFF_REF", "INS"} and not span_match:
        # Skip strict length agreement when the pred is a coarse chain-level
        # block that span-contains the truth: by construction |pred.svlen| is
        # much larger than |truth.svlen| in that case, and length disagreement
        # is the expected signal rather than a mismatch.
        denom = max(abs(truth.svlen), 1)
        if abs(abs(truth.svlen) - abs(pred.svlen)) / denom > tol_frac:
            return False
    if canonical_group(truth.svtype) == "TRA" or canonical_group(pred.svtype) == "TRA":
        truth_has_mate = _has_mate(truth)
        pred_has_mate = _has_mate(pred)
        if truth_has_mate or pred_has_mate:
            if not (truth_has_mate and pred_has_mate):
                return False
            if truth.mate_contig != pred.mate_contig:
                return False
            if abs(truth.mate_pos - pred.mate_pos) > tol_bp:
                return False
    return True


def call_distance(truth: NormalizedCall, pred: NormalizedCall) -> int:
    if _span_contain_applies(truth, pred):
        pos_d = 0
    else:
        inv_or_tra = canonical_group(truth.svtype) == "TRA" or canonical_group(pred.svtype) == "TRA" or truth.svtype == "INV" or pred.svtype == "INV"
        if inv_or_tra:
            pos_d = abs(truth.pos - (pred.pos + abs(pred.svlen) // 2))
        else:
            pos_d = abs(truth.pos - pred.pos)
    len_d = 0 if canonical_group(truth.svtype) == "TRA" else abs(abs(truth.svlen) - abs(pred.svlen))
    mate_d = 0
    if (
        (canonical_group(truth.svtype) == "TRA" or canonical_group(pred.svtype) == "TRA")
        and _has_mate(truth)
        and _has_mate(pred)
    ):
        mate_d = abs(truth.mate_pos - pred.mate_pos)
    return pos_d + len_d + mate_d


def build_consensus_truth(
    callsets: list[list[NormalizedCall]],
    *,
    min_support: int = 2,
) -> list[NormalizedCall]:
    """Return SV calls supported by at least min_support of the input callsets.

    Two calls "support" each other iff calls_compatible() — same coord space,
    same canonical SV type, position within DEFAULT_TOL_BP, length within
    DEFAULT_TOL_LEN_FRAC. The returned set has one representative per cluster
    (the earliest-encountered call), so |consensus| <= |smallest input|.
    """
    if min_support < 1 or not callsets:
        return []
    flat: list[tuple[int, NormalizedCall]] = []
    for src_idx, calls in enumerate(callsets):
        for c in calls:
            flat.append((src_idx, c))
    if not flat:
        return []
    n = len(flat)
    parent = list(range(n))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        a, b = find(i), find(j)
        if a != b:
            parent[a] = b

    for i in range(n):
        for j in range(i + 1, n):
            # Consensus means support from independent callsets. Do not let
            # duplicate/nearby calls from one comparator bridge two otherwise
            # separate clusters: A1~B and A2~B should not make A1~A2 support
            # each other as if they were independent evidence.
            if flat[i][0] == flat[j][0]:
                continue
            if calls_compatible(flat[i][1], flat[j][1]):
                union(i, j)

    clusters: dict[int, list[tuple[int, NormalizedCall]]] = defaultdict(list)
    for i, item in enumerate(flat):
        clusters[find(i)].append(item)

    out: list[NormalizedCall] = []
    for members in clusters.values():
        sources = {src for src, _ in members}
        if len(sources) >= min_support:
            out.append(members[0][1])
    return out


def match_calls(truth_calls: list[NormalizedCall], pred_calls: list[NormalizedCall]) -> tuple[set[int], list[int]]:
    pairs: list[tuple[int, int, int]] = []
    for truth_idx, truth in enumerate(truth_calls):
        for pred_idx, pred in enumerate(pred_calls):
            if calls_compatible(truth, pred):
                pairs.append((call_distance(truth, pred), truth_idx, pred_idx))
    pairs.sort()

    used: set[int] = set()
    used_truth: set[int] = set()
    for _dist, truth_idx, pred_idx in pairs:
        if truth_idx in used_truth or pred_idx in used:
            continue
        used_truth.add(truth_idx)
        used.add(pred_idx)

    missed_truth = [idx for idx in range(len(truth_calls)) if idx not in used_truth]
    return used, missed_truth


def score_callsets(truth_calls: list[NormalizedCall], pred_calls: list[NormalizedCall]) -> dict[str, Any]:
    if not truth_calls:
        return {
            "tp": float("nan"),
            "fp": float("nan"),
            "fn": float("nan"),
            "precision": float("nan"),
            "recall": float("nan"),
            "f1": float("nan"),
            "precision_ci95": (float("nan"), float("nan")),
            "recall_ci95": (float("nan"), float("nan")),
            "status": "no_truth",
        }
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
        "status": "ok",
    }


# All canonical SV types we stratify per-(query, comparator) by. ALL is the
# aggregate row (already produced by score_callsets). The biology_findings /
# visualization layers expect every per-svtype row to be present even when
# the (truth, pred) intersection is empty for that type, so the per-type
# loop below always emits a row per known canonical type.
_CANONICAL_SV_TYPES: tuple[str, ...] = ("INS", "DEL", "INV", "DUP", "TRA", "OFF_REF")


def score_callsets_by_svtype(
    truth_calls: list[NormalizedCall],
    pred_calls: list[NormalizedCall],
) -> dict[str, dict[str, Any]]:
    """Return {svtype -> metrics} stratified per canonical SV type, plus an
    "ALL" key with aggregate metrics. Filters both truth and predictions to
    each type before scoring so precision is computed against same-type
    predictions only — a comparator that emits zero INS calls but many DEL
    calls scores 0/0 on INS rather than diluting the DEL row.
    """
    out: dict[str, dict[str, Any]] = {"ALL": score_callsets(truth_calls, pred_calls)}
    for svtype in _CANONICAL_SV_TYPES:
        t_sub = [c for c in truth_calls if c.svtype == svtype]
        p_sub = [c for c in pred_calls if c.svtype == svtype]
        out[svtype] = score_callsets(t_sub, p_sub)
    return out


# ── Fungal-specific leave-one-out comparator-variance benchmark ──────────
#
# Single-number F1 against comparator consensus hides how much of the metric is
# driven by *which* comparators happened to be in the pool. LOO replays the
# score K times, each time excluding one comparator; the F1 dispersion is
# the "comparator-induced" component. Folds are also stratified by fungal-
# specific axes — length bin (sub-TE / TE element / TE cluster / arm) and
# element class (TE_LTR / TE_TIR / STARSHIP / HGT / RIP / ...) — so the
# reader can tell whether the variance is in boring small INDELs or in the
# biologically interesting HGT / STARSHIP / >50 kb arm rearrangements.
#
# Length-bin boundaries reflect fungal SV biology:
#   <500 bp        : sub-TE fragments / small accessory
#   500 bp–5 kb    : full TE element (LTR retrotransposon, helitron, TIR)
#   5–50 kb        : TE cluster, Starship cargo, accessory chromosome chunk
#   >50 kb         : whole-arm rearrangement / accessory chromosome
_FUNGAL_LEN_BINS: tuple[tuple[int, int, str], ...] = (
    (50,        500,         "lt_500bp_subTE"),
    (500,       5_000,       "500bp_5kb_TE_element"),
    (5_000,     50_000,      "5kb_50kb_TE_cluster_starship"),
    (50_000,    10**12,      "gt_50kb_arm_or_accessory"),
)


def _fungal_length_bin(svlen: int) -> str:
    a = abs(int(svlen or 0))
    for lo, hi, label in _FUNGAL_LEN_BINS:
        if lo <= a < hi:
            return label
    return "lt_50bp_below_threshold"


_FUNGAL_ELEMENT_CLASSES: tuple[str, ...] = (
    "TE_LTR", "TE_TIR", "TE_LINE", "TE_SINE",
    "STARSHIP", "HGT", "RIP", "REPEAT", "NONE",
)


def _stratify_calls_fungal(
    calls: list[NormalizedCall],
    query_phylum: str,
) -> dict[str, dict[str, list[NormalizedCall]]]:
    out: dict[str, dict[str, list[NormalizedCall]]] = {
        "LENGTH_BIN":    defaultdict(list),
        "ELEMENT_CLASS": defaultdict(list),
        "PHYLUM":        defaultdict(list),
    }
    phylum = (query_phylum or ".").strip() or "."
    for c in calls:
        # TRA calls have svlen=1 (breakpoint) per the VCF parser at
        # line ~4503, so the numeric bin would mis-route every TRA into
        # the "<50 bp sub-TE" bucket. TRAs get their own stratum so they
        # are scored as the structural events they are, independent of
        # the breakpoint coordinate distance.
        if (c.svtype or "").upper() == "TRA":
            out["LENGTH_BIN"]["TRA_breakpoint"].append(c)
        else:
            out["LENGTH_BIN"][_fungal_length_bin(c.svlen)].append(c)
        ec = (c.element_class or "NONE").upper()
        if ec not in _FUNGAL_ELEMENT_CLASSES:
            ec = "NONE"
        out["ELEMENT_CLASS"][ec].append(c)
        out["PHYLUM"][phylum].append(c)
    return out


def _score_strata_fungal(
    truth_calls: list[NormalizedCall],
    pred_calls: list[NormalizedCall],
    query_phylum: str,
) -> dict[str, dict[str, dict[str, Any]]]:
    """Per-fungal-stratum scoring: each stratum filters BOTH truth and
    predictions to the same key before scoring, so e.g. STARSHIP precision
    is computed against STARSHIP predictions only.
    """
    t = _stratify_calls_fungal(truth_calls, query_phylum)
    p = _stratify_calls_fungal(pred_calls,  query_phylum)
    out: dict[str, dict[str, dict[str, Any]]] = {}
    for sk in ("LENGTH_BIN", "ELEMENT_CLASS", "PHYLUM"):
        out[sk] = {}
        for k in sorted(set(t[sk].keys()) | set(p[sk].keys())):
            out[sk][k] = score_callsets(t[sk].get(k, []), p[sk].get(k, []))
    return out


def score_loo_consensus(
    pred_calls: list[NormalizedCall],
    comparator_callsets: dict[str, list[NormalizedCall]],
    *,
    min_support: int = 2,
    query_phylum: str = ".",
) -> dict[str, Any]:
    """Leave-one-out comparator-variance benchmark.

    For K comparator callsets:
      • Baseline: K-comparator consensus_(min_support) of K.
      • For each comparator i, exclude i and compute consensus on the
        remaining K-1 → leave-one-out fold; score pred against it.
    Returns per-fold metrics + F1 mean/stdev/range + the most influential
    comparator (the one whose exclusion shifts F1 most). All folds are
    additionally broken down by fungal length bin / element class / phylum.

    Skipped when K < min_support + 1 — without a margin, LOO folds would
    drop below the support threshold and emit empty truth.
    """
    K = len(comparator_callsets)
    if K < min_support + 1:
        return {
            "status":      "skipped",
            "reason":      f"need >= {min_support + 1} comparators for LOO at min_support={min_support}, got {K}",
            "comparators": sorted(comparator_callsets.keys()),
        }

    labels = sorted(comparator_callsets.keys())

    full_truth  = build_consensus_truth(
        [comparator_callsets[lbl] for lbl in labels], min_support=min_support,
    )
    full_metrics = score_callsets(full_truth, pred_calls)
    full_strata  = _score_strata_fungal(full_truth, pred_calls, query_phylum)

    folds: dict[str, dict[str, Any]] = {}
    f1_vals: list[float] = []
    for excluded in labels:
        loo_truth = build_consensus_truth(
            [comparator_callsets[lbl] for lbl in labels if lbl != excluded],
            min_support=min_support,
        )
        loo_metrics = score_callsets(loo_truth, pred_calls)
        folds[excluded] = {
            "truth_n":  len(loo_truth),
            "metrics":  loo_metrics,
            "strata":   _score_strata_fungal(loo_truth, pred_calls, query_phylum),
        }
        f1 = loo_metrics.get("f1")
        if isinstance(f1, (int, float)) and not math.isnan(f1):
            f1_vals.append(float(f1))

    if f1_vals:
        f1_mean  = statistics.fmean(f1_vals)
        f1_stdev = statistics.pstdev(f1_vals) if len(f1_vals) > 1 else 0.0
        f1_min, f1_max = min(f1_vals), max(f1_vals)
        f1_swing = f1_max - f1_min
    else:
        f1_mean = f1_stdev = f1_min = f1_max = f1_swing = float("nan")

    base_f1 = full_metrics.get("f1")
    most_influential: dict[str, Any] | None = None
    if isinstance(base_f1, (int, float)) and not math.isnan(base_f1):
        best_abs = 0.0
        for lbl, info in folds.items():
            f1 = info["metrics"].get("f1")
            if not isinstance(f1, (int, float)) or math.isnan(f1):
                continue
            delta = float(f1) - float(base_f1)
            if abs(delta) > best_abs:
                best_abs = abs(delta)
                most_influential = {"label": lbl, "delta_f1": delta}

    # Verdict requires ≥2 valid fold F1 values, otherwise the "robust" label
    # is unearned (every fold's consensus could have been empty). >5pp swing
    # over ≥2 folds → comparator-driven; reader should not quote a single F1
    # number without the swing alongside it.
    if len(f1_vals) < 2:
        verdict = "insufficient_folds"
    elif isinstance(f1_swing, float) and not math.isnan(f1_swing) and f1_swing > 0.05:
        verdict = "high_variance_comparator_driven"
    else:
        verdict = "low_variance_robust"

    return {
        "status":                       "ok",
        "comparators":                  labels,
        "comparators_n":                K,
        "min_support":                  min_support,
        "phylum":                       query_phylum or ".",
        "baseline_full_consensus": {
            "truth_n": len(full_truth),
            "metrics": full_metrics,
            "strata":  full_strata,
        },
        "loo_folds":                    folds,
        "f1_full_baseline":             base_f1,
        "f1_mean_loo":                  f1_mean,
        "f1_stdev_loo":                 f1_stdev,
        "f1_range_loo":                 [f1_min, f1_max],
        "f1_swing_loo":                 f1_swing,
        "most_influential_comparator":  most_influential,
        "verdict":                      verdict,
    }


def _emit_loo_summary_rows(
    *,
    query_asm:     str,
    query_phylum:  str,
    coord_space:   str,
    loo:           dict[str, Any],
) -> list[dict[str, Any]]:
    """Flatten a score_loo_consensus() result into per-(fold, stratum) TSV
    rows for loo_consensus_summary.tsv. Always emits the baseline as
    excluded_comparator='NONE' so downstream plots have a reference point.
    """
    rows: list[dict[str, Any]] = []
    if loo.get("status") != "ok":
        return rows

    def push(excluded: str, stratum_type: str, stratum_value: str,
             metrics: dict[str, Any]) -> None:
        rows.append({
            "query_asm":           query_asm,
            "phylum":              query_phylum or ".",
            "coordinate_space":    coord_space,
            "excluded_comparator": excluded,
            "stratum_type":        stratum_type,
            "stratum_value":       stratum_value,
            "truth_n":             metrics.get("tp", 0) + metrics.get("fn", 0)
                                    if isinstance(metrics.get("tp"), int) else 0,
            "tp":                  metrics.get("tp", float("nan")),
            "fp":                  metrics.get("fp", float("nan")),
            "fn":                  metrics.get("fn", float("nan")),
            "precision":           metrics.get("precision", float("nan")),
            "recall":              metrics.get("recall", float("nan")),
            "f1":                  metrics.get("f1", float("nan")),
            "status":              metrics.get("status", "unknown"),
        })

    base = loo["baseline_full_consensus"]
    push("NONE", "ALL", "ALL", base["metrics"])
    for sk, by_key in base["strata"].items():
        for k, m in by_key.items():
            push("NONE", sk, k, m)

    for excluded, info in loo["loo_folds"].items():
        push(excluded, "ALL", "ALL", info["metrics"])
        for sk, by_key in info["strata"].items():
            for k, m in by_key.items():
                push(excluded, sk, k, m)

    return rows


def _summarize_loo_variance_global(
    by_query: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    """Aggregate per-query LOO variance into one panel-wide summary.

    Reports separately for `query` and `reference` coord spaces. Per coord
    space: mean baseline F1, mean LOO stdev, mean LOO swing, and the most
    common "most-influential comparator" across queries (a vote).
    """
    out: dict[str, Any] = {}
    coord_spaces: set[str] = set()
    for per_cs in by_query.values():
        coord_spaces.update(per_cs.keys())

    for cs in sorted(coord_spaces):
        base_f1s: list[float] = []
        stdevs:   list[float] = []
        swings:   list[float] = []
        influencer_votes: dict[str, int] = defaultdict(int)
        n_queries = 0
        n_high_variance = 0
        for per_cs in by_query.values():
            entry = per_cs.get(cs)
            if not entry or entry.get("status") != "ok":
                continue
            n_queries += 1
            bf = entry.get("f1_full_baseline")
            st = entry.get("f1_stdev_loo")
            sw = entry.get("f1_swing_loo")
            if isinstance(bf, (int, float)) and not math.isnan(bf): base_f1s.append(float(bf))
            if isinstance(st, (int, float)) and not math.isnan(st): stdevs.append(float(st))
            if isinstance(sw, (int, float)) and not math.isnan(sw): swings.append(float(sw))
            if entry.get("verdict") == "high_variance_comparator_driven":
                n_high_variance += 1
            mi = entry.get("most_influential_comparator")
            if isinstance(mi, dict) and mi.get("label"):
                influencer_votes[str(mi["label"])] += 1

        top_influencer = None
        if influencer_votes:
            best = max(influencer_votes.items(), key=lambda kv: kv[1])
            top_influencer = {"label": best[0], "vote_count": best[1]}

        def _mean(xs: list[float]) -> float:
            return statistics.fmean(xs) if xs else float("nan")

        out[cs] = {
            "queries_with_loo":          n_queries,
            "queries_high_variance":     n_high_variance,
            "mean_baseline_f1":          _mean(base_f1s),
            "mean_loo_stdev":            _mean(stdevs),
            "mean_loo_swing":            _mean(swings),
            "most_influential_overall":  top_influencer,
            "influencer_vote_counts":    dict(influencer_votes),
        }
    return out


def _sanitize_for_json(obj: Any) -> Any:
    """Recursively convert NaN / inf floats to None and tuples to lists so
    the LOO doc dumps to *strict* JSON (parseable by jq, JS JSON.parse, and
    schema validators). json.dump's default `allow_nan=True` writes literal
    `NaN` tokens which Python parses back but every other JSON parser
    rejects — so we walk the doc explicitly instead of relying on the
    `default=` callback (which is never invoked for floats since floats
    are technically serializable)."""
    if isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            return None
        return obj
    if isinstance(obj, dict):
        return {k: _sanitize_for_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_sanitize_for_json(v) for v in obj]
    return obj


def diagnose_match_failures(
    truth_calls: list[NormalizedCall],
    pred_calls: list[NormalizedCall],
) -> list[dict[str, Any]]:
    """For every predicted call that did NOT match any truth call, return a
    row describing the closest truth candidate (same contig + same svtype
    group) and the specific reason matching failed:

      contig_mismatch   — no truth call shares the predicted contig
      type_mismatch     — same contig but no truth call shares the SV type
      pos_out_of_tol    — same contig+type but breakpoint > DEFAULT_TOL_BP
      svlen_out_of_tol  — close enough on pos but svlen ratio > tol_frac
      mate_mismatch     — TRA mate contig/pos disagrees
      no_truth_for_type — empty truth subset for the SV type

    These rows feed a match_failures.tsv that lets the operator see at a
    glance why precision/recall is zero without re-running the benchmark.
    """
    used, _missed = match_calls(truth_calls, pred_calls)
    truth_by_contig_type: dict[tuple[str, str], list[tuple[int, NormalizedCall]]] = defaultdict(list)
    for ti, t in enumerate(truth_calls):
        truth_by_contig_type[(effective_contig(t), canonical_group(t.svtype))].append((ti, t))
    rows: list[dict[str, Any]] = []
    for pi, pred in enumerate(pred_calls):
        if pi in used:
            continue
        contig = effective_contig(pred)
        group = canonical_group(pred.svtype)
        cohort = truth_by_contig_type.get((contig, group), [])
        if not cohort:
            same_contig_any_type = any(effective_contig(t) == contig for t in truth_calls)
            reason = "type_mismatch" if same_contig_any_type else "contig_mismatch"
            closest_pos_delta: str = "."
            closest_svlen_delta: str = "."
            closest_truth_idx: str = "."
        else:
            best = min(cohort, key=lambda r: abs(r[1].pos - pred.pos))
            t_idx, t = best
            tol_bp = DEFAULT_TOL_BP.get(group, 500)
            tol_frac = DEFAULT_TOL_LEN_FRAC.get(group, 0.30)
            pos_delta = abs(t.pos - pred.pos)
            len_delta = abs(abs(t.svlen) - abs(pred.svlen))
            denom = max(abs(t.svlen), 1)
            span_match = _span_contain_applies(t, pred)
            if pos_delta > tol_bp and not span_match:
                reason = "pos_out_of_tol"
            elif group not in {"INV", "TRA", "OFF_REF", "INS"} \
                    and not span_match \
                    and len_delta / denom > tol_frac:
                reason = "svlen_out_of_tol"
            elif group == "TRA" and (_has_mate(t) or _has_mate(pred)) and (
                t.mate_contig != pred.mate_contig or abs(t.mate_pos - pred.mate_pos) > tol_bp
            ):
                reason = "mate_mismatch"
            else:
                # Should be a TP if we got here — usually means pred was bumped
                # by another pred grabbing the same truth in the greedy match.
                reason = "claimed_by_other_pred"
            closest_pos_delta = str(pos_delta)
            closest_svlen_delta = str(len_delta)
            closest_truth_idx = str(t_idx)
        rows.append({
            "pred_contig": contig,
            "pred_pos": pred.pos,
            "pred_end": pred.end,
            "pred_svtype": pred.svtype,
            "pred_svlen": pred.svlen,
            "reason": reason,
            "closest_truth_idx": closest_truth_idx,
            "closest_pos_delta": closest_pos_delta,
            "closest_svlen_delta": closest_svlen_delta,
        })
    return rows


def validation_basis_for_label(label: str) -> str:
    """Classify what a benchmark row is actually anchored to.

    The TSV keeps the historical `truth_label` column for compatibility, but
    for real fungal data comparator calls are baseline agreement, not ground
    truth.  Consumers should use this column to decide whether a row represents
    independent raw-read validation or only agreement with another algorithm.
    """
    if label == "read_level_union":
        return "raw_read_validated"
    if label.endswith("_read_supported"):
        return "comparator_agreement_read_supported"
    if label.startswith("consensus_"):
        return "comparator_agreement"
    if label == "no_comparator":
        return "no_independent_validation"
    return "comparator_baseline"


def _emit_per_svtype_rows(
    *,
    query_asm: str,
    coord_space: str,
    truth_label: str,
    method: str,
    truth_calls: list[NormalizedCall],
    pred_calls: list[NormalizedCall],
) -> list[dict[str, Any]]:
    """Return one agreement row per (svtype) plus the aggregate ALL row.

    Lifts the row-shape from a single score_callsets() call into an
    agreement-table-ready list, with each row carrying the `svtype` column
    that the visualization layer uses to score "MycoSV vs comparator per SV
    type" — the headline view requested for the real-data benchmark.
    """
    metrics_by_type = score_callsets_by_svtype(truth_calls, pred_calls)
    rows: list[dict[str, Any]] = []
    for svtype, m in metrics_by_type.items():
        if svtype == "ALL":
            t_count = len(truth_calls)
            p_count = len(pred_calls)
        else:
            t_count = sum(1 for c in truth_calls if c.svtype == svtype)
            p_count = sum(1 for c in pred_calls if c.svtype == svtype)
        rows.append({
            "query_asm": query_asm,
            "coordinate_space": coord_space,
            "truth_label": truth_label,
            "validation_basis": validation_basis_for_label(truth_label),
            "svtype": svtype,
            "method": method,
            "truth_calls": t_count,
            "pred_calls": p_count,
            "tp": m["tp"],
            "fp": m["fp"],
            "fn": m["fn"],
            "precision": m["precision"],
            "recall": m["recall"],
            "f1": m["f1"],
            "prec_lo95": m["precision_ci95"][0],
            "prec_hi95": m["precision_ci95"][1],
            "rec_lo95": m["recall_ci95"][0],
            "rec_hi95": m["recall_ci95"][1],
            "status": m.get("status", "ok"),
        })
    return rows


def parse_other_spec(spec: str) -> tuple[str, Path]:
    if "=" not in spec:
        raise ValueError(f"Expected label=path, got {spec!r}")
    label, path = spec.split("=", 1)
    label = label.strip()
    path = path.strip()
    if not label or not path:
        raise ValueError(f"Expected label=path, got {spec!r}")
    return label, Path(path).resolve()


# ============================================================================
# Read-level (raw FASTQ) independent validation of truth-set SVs
#
# Algorithm-based truth (minigraph, cactus, svim_asm, anchorwave, pggb, syri)
# all share the same blind spot: they call SVs from assemblies, so any event
# that is an assembly artefact (mis-join, hetero collapse, contig-edge gap)
# inherits "truth" status. Read-level validation re-anchors each candidate
# in the raw FASTQ by counting split / supplementary alignments that span
# the breakpoint coordinates. SVs without raw-read evidence are dropped from
# the candidate set before scoring, removing the largest source of caller-bias.
# ============================================================================


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def _cigar_ref_len(cigar: str) -> int:
    return sum(int(n) for n, op in CIGAR_RE.findall(cigar) if op in {"M", "D", "N", "=", "X"})


def _within_bp(a: int, b: int, flank_bp: int) -> bool:
    return abs(a - b) <= flank_bp


def _len_compatible(observed_len: int, expected_len: int | None, svtype: str | None) -> bool:
    if expected_len is None or expected_len <= 0 or svtype in {None, "INS", "TRA", "OFF_REF"}:
        return True
    tol_frac = DEFAULT_TOL_LEN_FRAC.get(canonical_group(svtype), 0.30)
    return abs(observed_len - expected_len) / max(expected_len, 1) <= max(tol_frac, 0.50)


def _cigar_indel_supports_call(
    cigar: str,
    ref_start: int,
    pos: int,
    end: int,
    *,
    svtype: str | None,
    svlen: int | None,
    flank_bp: int,
) -> bool:
    ref_cursor = ref_start
    expected_len = abs(svlen) if svlen is not None else None
    for raw_n, op in CIGAR_RE.findall(cigar):
        n = int(raw_n)
        if op in {"M", "=", "X"}:
            ref_cursor += n
            continue
        if op in {"D", "N"}:
            event_start = ref_cursor
            event_end = ref_cursor + n - 1
            ref_cursor += n
            if svtype not in {None, "DEL"}:
                continue
            same_locus = (
                _within_bp(event_start, pos, flank_bp)
                or _within_bp(event_end, end, flank_bp)
                or (pos - flank_bp) <= event_start <= (end + flank_bp)
                or (pos - flank_bp) <= event_end <= (end + flank_bp)
            )
            if same_locus and _len_compatible(n, expected_len, svtype):
                return True
            continue
        if op == "I":
            event_pos = max(ref_start, ref_cursor - 1)
            if svtype not in {None, "INS", "DUP", "OFF_REF"}:
                continue
            if (_within_bp(event_pos, pos, flank_bp) or _within_bp(event_pos, end, flank_bp)) and _len_compatible(n, expected_len, svtype):
                return True
            continue
        if op in {"S", "H", "P"}:
            continue
    return False


def _split_or_clip_supports_call(
    cigar: str,
    ref_start: int,
    ref_end: int,
    fields: list[str],
    pos: int,
    end: int,
    *,
    svtype: str | None,
    flank_bp: int,
    min_clip: int,
) -> bool:
    has_sa = any(f.startswith("SA:Z:") for f in fields[11:])
    left_clip = 0
    right_clip = 0
    m = re.match(r"^(\d+)[SH]", cigar)
    if m:
        left_clip = int(m.group(1))
    m = re.search(r"(\d+)[SH]$", cigar)
    if m:
        right_clip = int(m.group(1))
    clipped = left_clip >= min_clip or right_clip >= min_clip
    if not (has_sa or clipped):
        return False
    bp_set = {pos, end}
    if left_clip >= min_clip and any(_within_bp(ref_start, bp, flank_bp) for bp in bp_set):
        return True
    if right_clip >= min_clip and any(_within_bp(ref_end, bp, flank_bp) for bp in bp_set):
        return True
    if has_sa:
        if any(_within_bp(ref_start, bp, flank_bp) or _within_bp(ref_end, bp, flank_bp) for bp in bp_set):
            return True
        # Whole-block INV/TRA alignments may carry supplementary evidence while
        # the embedded breakpoint lies inside the aligned interval.
        if svtype in {"INV", "TRA"} and any((ref_start - flank_bp) <= bp <= (ref_end + flank_bp) for bp in bp_set):
            return True
    return False


def _samtools_count_breakpoint_support(
    bam_path: Path,
    contig: str,
    pos: int,
    end: int,
    *,
    svtype: str | None = None,
    svlen: int | None = None,
    flank_bp: int = 250,
    min_clip: int = 30,
) -> int:
    """Return the number of reads with split/clipped alignments spanning the
    candidate breakpoint window. Uses `samtools view` and counts reads whose
    CIGAR carries a soft/hard clip ≥ min_clip on either side, or that have an
    SA tag (supplementary alignment), within ±flank_bp of the breakpoint.

    A breakpoint with ≥ K supporting reads (K=3 default elsewhere) is treated
    as raw-data confirmed; below K it is dropped from truth.
    """
    if not tool_path("samtools"):
        return 0
    region = f"{contig}:{max(1, pos - flank_bp)}-{end + flank_bp}"
    try:
        proc = subprocess.run(
            ["samtools", "view", "-F", "0x900", str(bam_path), region],
            text=True, capture_output=True, check=True,
            timeout=_TOOL_TIMEOUT,
        )
    except subprocess.CalledProcessError:
        return 0
    except subprocess.TimeoutExpired:
        return 0
    support = 0
    for line in proc.stdout.splitlines():
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            continue
        try:
            read_pos = int(fields[3])
        except ValueError:
            continue
        cigar = fields[5]
        ref_end = read_pos + max(0, _cigar_ref_len(cigar) - 1)
        if _cigar_indel_supports_call(
            cigar, read_pos, pos, end,
            svtype=svtype, svlen=svlen, flank_bp=flank_bp,
        ) or _split_or_clip_supports_call(
            cigar, read_pos, ref_end, fields, pos, end,
            svtype=svtype, flank_bp=flank_bp, min_clip=min_clip,
        ):
            support += 1
    return support


def _build_validation_bam(
    query_row: dict[str, str],
    work_dir: Path,
    threads: int,
) -> tuple[Path, Path] | None:
    """Build a sorted+indexed BAM of the query's raw reads vs benchmark ref.

    For assembly-mode queries this uses the assembly contigs themselves as
    "reads" (asm20 alignment — the benchmark ref is a diverged sibling
    clade), which still surfaces split-alignment evidence at every
    assembly-supported breakpoint. For reads-mode queries the appropriate
    long/short minimap2 preset is selected.
    """
    mode = (query_row.get("query_mode") or "assembly").lower()
    if mode == "assembly":
        if not tool_path("minimap2") or not tool_path("samtools"):
            return None
        ref_fasta = query_row.get("benchmark_ref_fasta", ".")
        if ref_fasta in {"", "."}:
            return None
        ref_fa = Path(ref_fasta).resolve()
        query_fa = locate_query_path(query_row)
        if not query_fa.exists():
            return None
        work_dir.mkdir(parents=True, exist_ok=True)
        ref_fa_plain = _ensure_plain_fasta(ref_fa, work_dir)
        if ref_fa_plain is None:
            return None
        sam_path = work_dir / "validation.sam"
        bam_sorted = work_dir / "validation.sorted.bam"
        if bam_sorted.exists() and (work_dir / "validation.sorted.bam.bai").exists():
            return bam_sorted, ref_fa_plain
        # asm20, not asm5: the benchmark reference is a held-out sibling-clade
        # genome. asm5 (<=5% divergence) fragments the contig-vs-ref alignment
        # for cross-species fungal pairs into short MAPQ-0 blocks, so the
        # split/clip breakpoint signal never lands where the call expects it
        # — every assembly-supported SV then fails read-validation and is
        # dropped from read_validated_truth.tsv. asm20 keeps the alignment
        # coherent across the genus-level divergence these pairs actually
        # have; matches the svim_asm / syri / clade-lift presets.
        with sam_path.open("wb") as sam_out:
            subprocess.run(
                ["minimap2", "-ax", "asm20", "--cs", "-t", str(threads),
                 str(ref_fa_plain), str(query_fa)],
                stdout=sam_out, stderr=subprocess.PIPE, check=True,
                timeout=_TOOL_TIMEOUT,
            )
        try:
            run(["samtools", "sort", "-@", str(threads), "-o", str(bam_sorted), str(sam_path)], cwd=ROOT)
            run(["samtools", "index", str(bam_sorted)], cwd=ROOT)
        except subprocess.CalledProcessError:
            return None
        try:
            sam_path.unlink()
        except OSError:
            pass
        return bam_sorted, ref_fa_plain

    preset = _long_read_preset(query_row) if mode == "long-reads" else "sr"
    aligned = _minimap2_align_reads(query_row, work_dir, threads, preset=preset)
    return aligned


def validate_calls_with_reads(
    truth_calls: list[NormalizedCall],
    query_row: dict[str, str],
    work_dir: Path,
    *,
    threads: int,
    min_support: int,
    flank_bp: int,
    force_external: bool = False,
) -> tuple[list[NormalizedCall], list[dict[str, Any]]]:
    """Re-anchor algorithm-derived candidate calls in the raw query data.

    force_external=True disables the "trust MycoSV's intrinsic support"
    shortcut: EVERY candidate — comparator or MycoSV — must clear the
    external split/clipped-read threshold. This is required when the input
    is a tool-agnostic union (read_level_union): letting MycoSV calls
    self-validate via their own anchor/cluster counts would make precision
    against that truth circular.

    Returns (filtered_truth, per_call_support_rows). filtered_truth contains
    only calls with >= min_support split/clipped read support at the
    breakpoint. per_call_support_rows is the per-SV evidence record for the
    on-disk read_validated_truth.tsv (always written even for dropped calls).
    """
    def internal_support(call: NormalizedCall) -> int | None:
        return call.read_support if call.source == "mycosv" and call.read_support is not None else None

    def support_source(call: NormalizedCall, intrinsic: int | None) -> str:
        if intrinsic is None:
            return "external_validation"
        mode = (query_row.get("query_mode") or "assembly").lower()
        if mode == "assembly":
            return "mycosv_assembly_anchors"
        if mode == "long-reads":
            return "mycosv_long_read_cluster"
        if mode == "short-reads":
            return "mycosv_short_read_kmer"
        return "mycosv_internal"

    def support_row(
        call: NormalizedCall,
        *,
        intrinsic: int | None,
        validation_support: int | None,
        validated: bool | None,
        status: str,
    ) -> dict[str, Any]:
        return {
            "query_asm": query_row.get("query_asm", "."),
            "ref_contig": call.ref_contig if call.coord_space == "reference" else call.query_contig,
            "pos": call.pos,
            "end": call.end,
            "svtype": call.svtype,
            "source": call.source,
            "coord_space": call.coord_space,
            "read_support": intrinsic if intrinsic is not None else -1,
            "validation_support": validation_support if validation_support is not None else -1,
            "support_source": support_source(call, intrinsic),
            "read_validated": "yes" if validated else ("unknown" if validated is None else "no"),
            "status": status,
        }

    def effective_flank(call: NormalizedCall) -> int:
        # The matcher (calls_compatible) uses DEFAULT_TOL_BP per SV type
        # (DEL/DUP 2500, INV/TRA 10000); a flat 250 bp validation window
        # rejected real DEL/DUP/INV/TRA calls whose split-read evidence sat
        # 300-2000 bp from the comparator's reported breakpoint, causing
        # asymmetric truth shrinkage and ~15-25 F1 hits on the larger SV
        # classes. Use 1/5 of the matcher tolerance (DEL/DUP -> 500 bp,
        # INV/TRA -> 2000 bp), floored at the user-set flank_bp so callers
        # who explicitly pass a wider window keep that wider behavior.
        per_type_tol = DEFAULT_TOL_BP.get(canonical_group(call.svtype), 0)
        return max(flank_bp, per_type_tol // 5) if per_type_tol > 0 else flank_bp

    aligned = _build_validation_bam(query_row, work_dir, threads)
    rows: list[dict[str, Any]] = []
    if aligned is None:
        has_intrinsic_support = any(internal_support(call) is not None for call in truth_calls)
        for call in truth_calls:
            intrinsic = internal_support(call)
            validated = (
                False if force_external
                else (intrinsic is not None and intrinsic >= min_support)
            )
            rows.append(support_row(
                call,
                intrinsic=intrinsic,
                validation_support=None,
                validated=None if force_external else (
                    validated if intrinsic is not None else None),
                status="validation_unavailable",
            ))
        # Without an alignment there is no external evidence — under
        # force_external nothing can be confirmed, so the validated set is empty
        # rather than falling back to (tool-biased) intrinsic counts.
        if force_external:
            return [], rows
        kept = [
            call for call in truth_calls
            if internal_support(call) is not None
            and internal_support(call) >= min_support
        ]
        if has_intrinsic_support:
            return kept, rows
        return list(truth_calls), rows
    bam_sorted, _ref_plain = aligned
    # Build the set of contigs actually present in the validation BAM header
    # once, up front. samtools view on an absent contig returns 0 reads
    # silently — without this guard we mis-attribute "no reads at breakpoint"
    # to the call, when in reality the contig is on a sibling clade we never
    # aligned to. Used below to skip the breakpoint scan for those calls.
    bam_contigs: frozenset[str] = frozenset()
    if tool_path("samtools"):
        try:
            hdr = subprocess.run(
                ["samtools", "view", "-H", str(bam_sorted)],
                text=True, capture_output=True, check=True, timeout=_TOOL_TIMEOUT,
            )
            bam_contigs = frozenset(
                ln.split("\t")[1][3:]
                for ln in hdr.stdout.splitlines()
                if ln.startswith("@SQ\t")
                and len(ln.split("\t")) >= 2
                and ln.split("\t")[1].startswith("SN:")
            )
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
            bam_contigs = frozenset()
    kept: list[NormalizedCall] = []
    for call in truth_calls:
        intrinsic = internal_support(call)
        mode = (query_row.get("query_mode") or "assembly").lower()
        # MycoSV's C++ pipeline already clusters at SUPPORT>=2; any call with
        # a usable intrinsic count is "raw-data confirmed" by construction.
        # Re-aligning short-read kmer/long-read cluster calls to a single
        # benchmark_ref BAM and re-counting split reads loses ~50% of them
        # purely because:
        #   (a) the call lives on a sibling clade the BAM never indexed, or
        #   (b) short-read assembly anchors land on a synthetic sr_unitig
        #       contig that the BAM does not contain.
        # Trust the intrinsic count, but still record the external split-read
        # count when available (for diagnostic value in the TSV).
        # force_external skips this shortcut entirely — see docstring.
        if (not force_external
                and call.source == "mycosv"
                and intrinsic is not None and intrinsic >= min_support):
            contig = call.ref_contig if call.coord_space == "reference" and call.ref_contig not in {"", "."} else call.query_contig
            ext_support: int | None
            if mode == "assembly" and call.coord_space == "query":
                ext_support = None
            elif bam_contigs and contig not in bam_contigs:
                ext_support = None  # sibling clade contig — would always read 0
            else:
                ext_support = _samtools_count_breakpoint_support(
                    bam_sorted, contig, call.pos, call.end,
                    svtype=call.svtype, svlen=call.svlen,
                    flank_bp=effective_flank(call), min_clip=30,
                )
            status_label = (
                "query_space_not_reference_validated"
                if mode == "assembly" and call.coord_space == "query"
                else ("validated_intrinsic" if call.coord_space == "query" else "validated")
            )
            rows.append(support_row(
                call,
                intrinsic=intrinsic,
                validation_support=ext_support,
                validated=True,
                status=status_label,
            ))
            kept.append(call)
            continue
        contig = call.ref_contig if call.coord_space == "reference" and call.ref_contig not in {"", "."} else call.query_contig
        if bam_contigs and contig not in bam_contigs:
            # Distinguish "external validation impossible" from
            # "external validation failed" so the TSV stays interpretable.
            # Under force_external a call we cannot externally check is NOT
            # admitted to the validated set (no intrinsic fallback).
            validated_flag = (
                False if force_external
                else (intrinsic is not None and intrinsic >= min_support)
            )
            rows.append(support_row(
                call,
                intrinsic=intrinsic,
                validation_support=None,
                validated=None if force_external else (
                    validated_flag if intrinsic is not None else None),
                status="contig_absent_from_validation_bam",
            ))
            if validated_flag:
                kept.append(call)
            continue
        validation_support = _samtools_count_breakpoint_support(
            bam_sorted, contig, call.pos, call.end,
            svtype=call.svtype, svlen=call.svlen,
            flank_bp=effective_flank(call), min_clip=30,
        )
        if call.svtype == "TRA" and _has_mate(call):
            mate_support = _samtools_count_breakpoint_support(
                bam_sorted,
                call.mate_contig,
                call.mate_pos,
                call.mate_end or call.mate_pos,
                svtype=call.svtype, svlen=call.svlen,
                flank_bp=effective_flank(call), min_clip=30,
            )
            validation_support = max(validation_support, mate_support)
        # force_external: external split-read count is the ONLY evidence that
        # counts — intrinsic (tool-derived) support is ignored for truth
        # membership so the union truth stays tool-agnostic.
        effective_support = (
            validation_support if force_external
            else max(validation_support, intrinsic or 0)
        )
        validated = effective_support >= min_support
        rows.append(support_row(
            call,
            intrinsic=intrinsic,
            validation_support=validation_support,
            validated=validated,
            status="validated" if validated else "not_validated",
        ))
        if validated:
            kept.append(call)
    return kept, rows


# ============================================================================
# MycoSV-to-benchmark-reference liftover: translate MycoSV's sibling-clade
# REFCONTIG/REFPOS into benchmark_ref_fasta coordinates so per-comparator PR
# scoring is apples to apples.
#
# MyCoSV's pangenomic routing picks the closest clade per region and reports
# REFPOS against THAT clade. The benchmark callers see exactly one
# benchmark_ref_fasta per query, so a call reported on NC_012864.1
# (S. cerevisiae) at pos 10567 cannot be matched against a comparator call on
# NC_089928.1 (the user-chosen ref) without coordinate translation. We
# minimap2-asm20 align the sibling clade FASTA to benchmark ref, parse the PAF,
# and translate REFCONTIG/REFPOS/REFEND for each call. PAFs are cached in
# <out_dir>/lift_cache/ so a re-run of the same panel does not pay the
# alignment cost twice.
# ============================================================================


def _clade_fasta_path(clade_name: str, data_cache_dir: Path) -> Path | None:
    """Resolve a mycosv INFO CLADE= field (basename like
    'GCF_000026945.1_ASM2694v1_genomic.fna') to its cached FASTA path.
    """
    if not clade_name or clade_name in {".", "OFF_REFERENCE"}:
        return None
    refs_dir = data_cache_dir / "refs"
    candidates = [
        refs_dir / clade_name,
        refs_dir / f"{clade_name}.gz",
    ]
    for cand in candidates:
        if cand.exists():
            return cand
    return None


def _ensure_clade_to_bench_paf(
    clade_fasta: Path,
    benchmark_ref_fasta: Path,
    cache_dir: Path,
    threads: int,
) -> Path | None:
    if not tool_path("minimap2"):
        return None
    if not clade_fasta.exists() or not benchmark_ref_fasta.exists():
        return None
    cache_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{clade_fasta.name.replace('.gz', '').replace('.fna', '').replace('.fasta', '')}"
    bench_stem = f"{benchmark_ref_fasta.name.replace('.gz', '').replace('.fna', '').replace('.fasta', '')}"
    paf_path = cache_dir / f"{stem}__to__{bench_stem}.paf"
    if paf_path.exists() and paf_path.stat().st_size > 0:
        return paf_path
    tmp_path = paf_path.with_suffix(".paf.part")
    # asm20 (≥80 % identity, divergence ≤20 %) captures cross-species fungal
    # synteny — Saccharomyces sensu stricto, Candida sister lineages, Puccinia
    # rusts — that asm5 (intra-species) misses. Without it the lift drops most
    # sibling-clade calls because mycosv routes broadly across the family tree.
    try:
        with tmp_path.open("wb") as out:
            subprocess.run(
                ["minimap2", "-cx", "asm20", "--secondary=no", "-t", str(threads),
                 str(benchmark_ref_fasta), str(clade_fasta)],
                stdout=out, stderr=subprocess.PIPE, check=True,
                timeout=_TOOL_TIMEOUT,
            )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        tmp_path.unlink(missing_ok=True)
        return None
    tmp_path.rename(paf_path)
    return paf_path


def _load_lift_table(paf_path: Path) -> dict[str, list[tuple[int, int, str, int, int, int]]]:
    """Parse a PAF into {query_contig: [(qs, qe, tname, ts, te, strand), …]}
    sorted by qs. Strand is +1 / -1.
    """
    table: dict[str, list[tuple[int, int, str, int, int, int]]] = defaultdict(list)
    with paf_path.open() as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue
            try:
                qs = int(cols[2])
                qe = int(cols[3])
                ts = int(cols[7])
                te = int(cols[8])
            except ValueError:
                continue
            qname = cols[0]
            tname = cols[5]
            strand = 1 if cols[4] == "+" else -1
            table[qname].append((qs, qe, tname, ts, te, strand))
    for q in table:
        table[q].sort(key=lambda r: r[0])
    return dict(table)


def _lift_pos(
    lift_table: dict[str, list[tuple[int, int, str, int, int, int]]],
    qname: str,
    qpos: int,
) -> tuple[str, int] | None:
    intervals = lift_table.get(qname)
    if not intervals:
        return None
    # Walk linearly — n_intervals per contig is small in practice (asm5 emits
    # one block per contiguous syntenic block). Binary search is only a win
    # for >>100 blocks, which would itself signal an unstable alignment.
    for qs, qe, tname, ts, te, strand in intervals:
        if qs <= qpos <= qe:
            offset = qpos - qs
            if strand > 0:
                return tname, ts + offset
            return tname, te - offset
    return None


def _lift_calls_to_benchmark_ref(
    calls: list[NormalizedCall],
    benchmark_ref_fasta: Path,
    data_cache_dir: Path,
    cache_dir: Path,
    threads: int,
) -> list[NormalizedCall]:
    """Return calls with REFCONTIG/REFPOS/REFEND translated to
    benchmark_ref_fasta coords where possible.

    Calls already on a benchmark-ref contig pass through unchanged. Calls that
    cannot be lifted are also retained on their native clade so diagnostic
    any-clade rows still see the full MycoSV prediction set; the strict
    benchmark-ref filter later excludes them from single-reference PR scoring.
    """
    if not calls:
        return []
    bench_contigs = fasta_contig_names(benchmark_ref_fasta)
    lift_tables_by_clade: dict[str, dict | None] = {}
    out: list[NormalizedCall] = []
    for call in calls:
        if call.coord_space != "reference":
            out.append(call)
            continue
        if call.ref_contig in bench_contigs:
            out.append(call)
            continue
        clade = call.ref_asm or "."
        if clade in {".", "", "OFF_REFERENCE"}:
            # No clade info — keep the call on its native ref_contig so it
            # still scores for `reference_any_clade` rows; the bench_contigs
            # filter at the caller will drop it from the strict reference
            # row. Previous behaviour silently discarded these.
            out.append(call)
            continue
        if clade not in lift_tables_by_clade:
            clade_fa = _clade_fasta_path(clade, data_cache_dir)
            if clade_fa is None:
                lift_tables_by_clade[clade] = None
            else:
                paf = _ensure_clade_to_bench_paf(
                    clade_fa, benchmark_ref_fasta, cache_dir, threads
                )
                lift_tables_by_clade[clade] = _load_lift_table(paf) if paf else None
        table = lift_tables_by_clade.get(clade)
        if not table:
            # PAF could not be built (minimap2 missing, clade FASTA absent,
            # or alignment too diverged for asm20). Keep the call on its
            # native ref_contig so the `reference_any_clade` row still sees
            # it; the bench_contigs filter at the caller drops it from the
            # strict reference row. Previously: silent drop, which made
            # pred_calls in exact_benchmark_summary.tsv much smaller than
            # the actual mycosv prediction set.
            out.append(call)
            continue
        lifted_start = _lift_pos(table, call.ref_contig, call.pos)
        if lifted_start is None:
            # PAF exists but does not cover this breakpoint. Keep the call
            # on its original clade contig — see comment above.
            out.append(call)
            continue
        lifted_end = _lift_pos(table, call.ref_contig, call.end) or lifted_start
        new_contig = lifted_start[0]
        new_pos = lifted_start[1]
        new_end = lifted_end[1]
        if new_end < new_pos:
            new_pos, new_end = new_end, new_pos
        if new_contig not in bench_contigs:
            # Lift landed on a non-benchmark contig; keep the original call
            # for any-clade scoring rather than dropping outright.
            out.append(call)
            continue
        mate_contig = call.mate_contig
        mate_pos = call.mate_pos
        mate_end = call.mate_end
        if _has_mate(call):
            lifted_mate = _lift_pos(table, call.mate_contig, call.mate_pos)
            if lifted_mate is None:
                out.append(call)
                continue
            lifted_mate_end = (
                _lift_pos(table, call.mate_contig, call.mate_end or call.mate_pos)
                or lifted_mate
            )
            mate_contig = lifted_mate[0]
            mate_pos = lifted_mate[1]
            mate_end = lifted_mate_end[1]
            if mate_end < mate_pos:
                mate_pos, mate_end = mate_end, mate_pos
            if mate_contig not in bench_contigs:
                out.append(call)
                continue
        out.append(replace(
            call,
            pos=new_pos,
            end=new_end,
            ref_contig=new_contig,
            mate_contig=mate_contig,
            mate_pos=mate_pos,
            mate_end=mate_end,
        ))
    return out


def reference_projection_locus_key(
    call: NormalizedCall,
    bucket_bp: int = 100,
) -> tuple[str, str, str, int, int, int, str, int]:
    """Approximate identity key for projected benchmark-reference calls.

    MycoSV may emit several raw pairwise observations for one biological SV
    because the same query locus is compared to multiple close pangenome refs.
    After liftover to the single benchmark reference, those observations should
    contribute one prediction to the single-reference PR comparison. The raw
    observations remain in calls.vcf and the pangenome layers; this key only
    defines the projected single-reference view.
    """
    bucket = max(1, bucket_bp)
    mate_contig = "."
    mate_bucket = -1
    if _has_mate(call):
        mate_contig = call.mate_contig or "."
        mate_bucket = max(0, call.mate_pos) // bucket
    return (
        call.query_asm,
        call.ref_contig or ".",
        call.svtype,
        max(0, call.pos) // bucket,
        max(0, call.end) // bucket,
        abs(call.svlen) // bucket,
        mate_contig,
        mate_bucket,
    )


def deduplicate_projected_reference_calls(
    paired_calls: list[tuple[NormalizedCall, tuple[str, str, int, int, str] | None]],
) -> tuple[
    list[NormalizedCall],
    list[tuple[str, str, int, int, str] | None],
    set[tuple[str, str, int, int, str]],
]:
    """Collapse projected MycoSV reference calls to one call per benchmark locus.

    Returns (deduped_ref_calls, representative_query_keys, all_query_keys).
    The representative keys stay parallel to deduped_ref_calls for comparator
    support propagation; all_query_keys marks every raw pangenome observation
    that is single-reference-equivalent for biology/novelty labelling.
    """
    by_locus: dict[
        tuple[str, str, str, int, int, int, str, int],
        tuple[NormalizedCall, tuple[str, str, int, int, str] | None],
    ] = {}
    all_query_keys: set[tuple[str, str, int, int, str]] = set()
    for call, qkey in paired_calls:
        if qkey is not None:
            all_query_keys.add(qkey)
        key = reference_projection_locus_key(call)
        if key not in by_locus:
            by_locus[key] = (call, qkey)
            continue
        prev_call, _prev_qkey = by_locus[key]
        prev_support = prev_call.read_support if prev_call.read_support is not None else -1
        new_support = call.read_support if call.read_support is not None else -1
        if new_support > prev_support:
            by_locus[key] = (call, qkey)
    deduped = list(by_locus.values())
    return (
        [call for call, _qkey in deduped],
        [qkey for _call, qkey in deduped],
        all_query_keys,
    )


def compile_binary_if_needed(binary_path: Path, force: bool = False) -> None:
    sources = [ROOT / "main.cpp", *ROOT.glob("*.hpp")]
    needs_build = force or not binary_path.exists()
    if not needs_build and binary_path.exists():
        bin_mtime = binary_path.stat().st_mtime
        needs_build = any(src.exists() and src.stat().st_mtime > bin_mtime for src in sources)
    if not needs_build:
        return
    run(["g++", "-O2", "-std=c++17", "-pthread", str(ROOT / "main.cpp"), "-o", str(binary_path)], cwd=ROOT)


def locate_query_path(query_row: dict[str, str]) -> Path:
    return Path(query_row["path"]).resolve()


def estimate_prepared_genome_size_hint(prepared_dir: Path) -> int:
    """Return the median genome size (bp) across every query's benchmark ref.

    The hint feeds the C++ binary's reads-mode auto-tuner (coverage estimate,
    minimum cluster support). Returning the *first* ref's size was wrong for
    mixed-phyla panels — a 12 Mbp yeast hint applied to a 30 Mbp Aspergillus
    query (or, worse, a 700 Mbp Gigaspora) makes coverage estimates orders
    of magnitude off. Median across all query benchmark refs is the
    representative aggregate; falls through to the corpus ref_list median
    when no query manifest exists.
    """
    query_manifest = prepared_dir / "query_manifest.tsv"
    candidate_paths: list[Path] = []
    if query_manifest.exists():
        for row in load_query_manifest(query_manifest):
            ref_fasta = (row.get("benchmark_ref_fasta") or "").strip()
            if ref_fasta and ref_fasta not in {".", ""}:
                candidate_paths.append(Path(ref_fasta))
    if not candidate_paths:
        ref_list = prepared_dir / "ref_list.txt"
        if ref_list.exists():
            for line in ref_list.read_text(encoding="utf-8").splitlines():
                line = line.strip()
                if line:
                    candidate_paths.append(Path(line))

    sizes: list[int] = []
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
            sizes.append(total)
    if not sizes:
        return 0
    sizes.sort()
    return sizes[len(sizes) // 2]


def run_mycosv(
    prepared_dir: Path,
    out_dir: Path,
    binary_path: Path,
    mode: str,
    extra_args: list[str],
    *,
    query_list_override: Path | None = None,
    threads: int = 8,
    max_clade_genomes: int = 8,
    reuse_index_dir: Path | None = None,
    reuse_registry_dir: Path | None = None,
    ref_list_override: Path | None = None,
) -> dict[str, str]:
    mycosv_dir = out_dir / "mycosv"
    mycosv_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = mycosv_dir / "calls"
    caller_args = list(extra_args)
    if ref_list_override is not None:
        # A benchmark-scoped ref list can still be large in million-real runs:
        # 256 fungal FASTAs are enough to spend tens of GiB on raw sequences +
        # suffix-array cache and can spin for a long time on the first query.
        # Keep flat fallback only for small override lists unless explicitly
        # requested. Hierarchical-only mode still produces a VCF/report, while
        # avoiding the memory-heavy all-ref MEM-chain fallback.
        # The million-real benchmark builds a bounded ref subset for downstream
        # comparison, but 256 fungal assemblies can still expand to tens of
        # thousands of contigs and several GiB of raw sequence. Keep flat
        # fallback to genuinely small debug/benchmark subsets by default; raise
        # MYCOSV_FLAT_REF_FALLBACK_MAX_REFS or set MYCOSV_FORCE_FLAT_REF_FALLBACK=1
        # when running on a node sized for all-vs-all MEM-chain rescue.
        flat_ref_limit = int(os.environ.get("MYCOSV_FLAT_REF_FALLBACK_MAX_REFS", "64"))
        try:
            with ref_list_override.open(encoding="utf-8") as fh:
                override_ref_count = sum(1 for line in fh if line.strip())
        except OSError:
            override_ref_count = flat_ref_limit + 1
        force_flat = os.environ.get("MYCOSV_FORCE_FLAT_REF_FALLBACK", "0") == "1"
        if force_flat or override_ref_count <= flat_ref_limit:
            caller_args = [a for a in caller_args if a != "--no-flat-ref-fallback"]
            sys.stderr.write(
                f"[mycosv] flat-ref fallback allowed for {override_ref_count} "
                f"benchmark refs (limit={flat_ref_limit})\n"
            )
        elif "--no-flat-ref-fallback" not in caller_args:
            caller_args.append("--no-flat-ref-fallback")
            sys.stderr.write(
                f"[mycosv] keeping flat-ref fallback disabled for "
                f"{override_ref_count} benchmark refs (limit={flat_ref_limit}); "
                "set MYCOSV_FORCE_FLAT_REF_FALLBACK=1 to override\n"
            )
    # Give read-mode auto-tuning a genome-size denominator even when the
    # query-list points at bounded benchmark FASTQs.
    if mode != "assembly" and "--genome-size-hint" not in caller_args:
        genome_size_hint = estimate_prepared_genome_size_hint(prepared_dir)
        if genome_size_hint > 0:
            caller_args.extend(["--genome-size-hint", str(genome_size_hint)])
    query_list_path = (query_list_override.resolve() if query_list_override
                       else (prepared_dir / "query_list.txt").resolve())

    # Use hierarchical routing when a hierarchy_manifest.tsv exists (real-data panels
    # with multi-species refs benefit greatly from routing; without it the binary would
    # try to load all refs into a flat graph and produce 0 calls across phyla).
    hierarchy_manifest = prepared_dir / "hierarchy_manifest.tsv"
    if hierarchy_manifest.exists() and (prepared_dir / "ref_list.txt").stat().st_size > 0:
        # If the caller passed a pre-built index (e.g. from
        # prepare-million-real), point the binary at it directly instead of
        # rebuilding. Saves the multi-hour rebuild on the million-real flow,
        # where the index has already been written next to the prepared dir.
        if reuse_index_dir is not None:
            idx_dir = reuse_index_dir.resolve()
            if not (idx_dir / "routing_manifest.tsv").exists():
                raise FileNotFoundError(
                    f"--reuse-index-dir {idx_dir} does not contain routing_manifest.tsv"
                )
            reg_dir = (reuse_registry_dir.resolve() if reuse_registry_dir
                       else (idx_dir.parent / "registry").resolve())
            sys.stderr.write(
                f"[mycosv] reusing pre-built routing index at {idx_dir} "
                f"(registry={reg_dir})\n"
            )
            # For small benchmark ref subsets, flat fallback gives anchored
            # reference-coordinate calls. For large subsets, the guard above
            # keeps --no-flat-ref-fallback so the million-real benchmark does
            # not load hundreds of FASTAs into the MEM-chain path.
            if ref_list_override is None and "--no-flat-ref-fallback" not in caller_args:
                caller_args.append("--no-flat-ref-fallback")
                sys.stderr.write(
                    "[mycosv] disabling flat reference fallback for reused "
                    "hierarchical index to keep memory bounded\n"
                )
            elif ref_list_override is not None:
                if "--no-flat-ref-fallback" in caller_args:
                    sys.stderr.write(
                        f"[mycosv] flat-ref fallback DISABLED for {ref_list_override}; "
                        "using hierarchical calls only for this benchmark subset\n"
                    )
                else:
                    sys.stderr.write(
                        f"[mycosv] flat-ref fallback ENABLED against {ref_list_override} "
                        "(benchmark-only ref subset; safe for memory)\n"
                    )
        else:
            idx_dir = mycosv_dir / "idx"
            reg_dir = mycosv_dir / "reg"
            idx_dir.mkdir(parents=True, exist_ok=True)
            reg_dir.mkdir(parents=True, exist_ok=True)
            # Build index only if not already present.
            if not (idx_dir / "routing_manifest.tsv").exists():
                build_cmd = [
                    str(binary_path.resolve()),
                    "--tol-hierarchical",
                    "--tol-build-index", str(hierarchy_manifest.resolve()),
                    "--tol-index-dir", str(idx_dir.resolve()),
                    "--tol-registry-dir", str(reg_dir.resolve()),
                    "--tol-multi-rank",
                    "--tol-base-graph-build",
                    "--tol-max-clade-genomes", str(max_clade_genomes),
                    "--tol-index-threads", str(threads),
                ]
                run_mycosv_command(build_cmd, cwd=ROOT)
        ref_list_path = (ref_list_override.resolve() if ref_list_override is not None
                         else (prepared_dir / "ref_list.txt").resolve())
        cmd = [
            str(binary_path.resolve()),
            "--tol-hierarchical",
            "--tol-index-dir", str(idx_dir.resolve()),
            "--tol-registry-dir", str(reg_dir.resolve()),
            "--ref-list", str(ref_list_path),
            "--query-list", str(query_list_path),
            "--out-prefix", str(out_prefix.resolve()),
            "--query-mode", mode,
            "--threads", str(threads),
            "--tol-index-threads", str(threads),
            *caller_args,
        ]
    else:
        ref_list_path = (ref_list_override.resolve() if ref_list_override is not None
                         else (prepared_dir / "ref_list.txt").resolve())
        cmd = [
            str(binary_path.resolve()),
            "--ref-list", str(ref_list_path),
            "--query-list", str(query_list_path),
            "--out-prefix", str(out_prefix.resolve()),
            "--query-mode", mode,
            *caller_args,
        ]
    # Stream stdout/stderr to disk in real time so a hung query is visible
    # via `tail -f calls.stderr.log` instead of looking dead until the
    # subprocess exits (or SLURM kills it).
    calls_stdout_log = mycosv_dir / "calls.stdout.log"
    calls_stderr_log = mycosv_dir / "calls.stderr.log"
    result = run_mycosv_command(
        cmd,
        cwd=ROOT,
        stream_stdout_path=calls_stdout_log,
        stream_stderr_path=calls_stderr_log,
    )
    # The streaming path already wrote to disk; the in-memory `result.stdout`
    # / `result.stderr` carry only the tail. The non-streaming path (older
    # callers) returns full bytes — write them too so existing log layout
    # is preserved either way.
    if result.stdout and not calls_stdout_log.exists():
        calls_stdout_log.write_text(result.stdout, encoding="utf-8")
    if result.stderr and not calls_stderr_log.exists():
        calls_stderr_log.write_text(result.stderr, encoding="utf-8")
    vcf_path = out_prefix.with_suffix(".vcf")
    hits_path = out_prefix.with_suffix(".hits.tsv")
    if not vcf_path.exists() or vcf_path.stat().st_size == 0:
        raise subprocess.CalledProcessError(
            90,
            cmd,
            stderr=f"MycoSV produced an empty or missing VCF: {vcf_path}",
        )
    data_records = 0
    try:
        with vcf_path.open(encoding="utf-8") as fh:
            has_vcf_header = False
            for line in fh:
                if line.startswith("#CHROM\t"):
                    has_vcf_header = True
                elif line and not line.startswith("#"):
                    data_records += 1
    except OSError as exc:
        raise subprocess.CalledProcessError(
            90,
            cmd,
            stderr=f"MycoSV VCF could not be read: {vcf_path}: {exc}",
        ) from exc
    if not has_vcf_header:
        raise subprocess.CalledProcessError(
            90,
            cmd,
            stderr=f"MycoSV VCF is malformed; missing #CHROM header: {vcf_path}",
        )
    stderr_text = getattr(result, "stderr", "") or ""
    fatal_empty_markers = (
        "std::bad_alloc",
        "no sequences after preprocessing",
        "cannot decompress",
        "query preprocessing/calling failure",
    )
    if data_records == 0 and any(marker in stderr_text for marker in fatal_empty_markers):
        raise subprocess.CalledProcessError(
            91,
            cmd,
            stderr=(
                f"MycoSV produced a header-only VCF after query failures: "
                f"{vcf_path}. See {mycosv_dir / 'calls.stderr.log'}"
            ),
        )
    if not hits_path.exists():
        raise subprocess.CalledProcessError(
            90,
            cmd,
            stderr=f"MycoSV produced no hits TSV: {hits_path}",
        )
    return {
        "vcf": str(vcf_path),
        "hits": str(hits_path),
        "gfa": str(out_prefix.with_suffix(".gfa")),
    }


def write_prefixed_fasta_records(src: Path, prefix: str, out_fh: io.TextIOBase) -> None:
    if str(src).endswith(".gz"):
        import gzip as _gzip
        fh_ctx = _gzip.open(str(src), "rt", encoding="utf-8")
    else:
        fh_ctx = src.open(encoding="utf-8")
    with fh_ctx as fh:
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
    # SyRI calls into pysam.FastaFile and refuses gzipped refs; same for the
    # query when it's piped from .gz. Decompress both into the work dir.
    ref_fa_plain = _ensure_plain_fasta(ref_fa, work_dir)
    if ref_fa_plain is None:
        return None
    query_fa_plain = _ensure_plain_fasta(query_fa, work_dir)
    if query_fa_plain is None:
        return None

    # SyRI's chrmatch heuristic (synsearchFunctions.pyx:438-465) auto-renames
    # query contigs to the ref chromosome they best align to, then sys.exit(1)s
    # with "Homologous chromosomes were not identified correctly" when more
    # than one query contig competes for the same ref chromosome.  That's the
    # default failure mode for fragmented MAG/bin queries (hundreds of small
    # contigs) against a chromosome-level ref (~16 chromosomes for a yeast).
    # Skip cleanly with a structured failure row rather than letting SyRI
    # deterministically fail and emit a 30-line cmdline warning.  The 4x
    # threshold is empirical: SyRI handles modest fragmentation but breaks
    # once the query is dramatically more fragmented than the reference.
    ref_n, _, _ = _fasta_stats(ref_fa_plain)
    qry_n, _, _ = _fasta_stats(query_fa_plain)
    if ref_n > 0 and qry_n > max(20, ref_n * 4):
        sys.stderr.write(
            f"[skip] syri {query_asm}: query is too fragmented "
            f"(query={qry_n} contigs vs ref={ref_n}); SyRI's chrmatch heuristic "
            "requires comparable contig counts.  Use minigraph / pggb / svim_asm instead.\n"
        )
        _log_comparator_failure(
            out_dir, "syri", query_asm,
            f"skipped_incompatible_chromosome_count qry_contigs={qry_n} ref_contigs={ref_n}",
        )
        return None

    sam_path = work_dir / "query_vs_ref.sam"
    # SyRI joins --dir and --prefix as plain strings before opening its log
    # file, so an absolute --prefix combined with the default --dir (cwd) gets
    # the doubled "<cwd>/<abs-prefix>" path that crashes dictConfig with
    # "Unable to configure handler 'log_file'". Pass --dir explicitly and keep
    # --prefix as a basename so the join always lands inside work_dir.
    # asm20, not asm5: the benchmark reference is a held-out sibling-clade
    # genome (cross-species fungal pair). asm5 fragments such alignments and
    # SyRI then rejects the pair as non-syntenic; asm20 (<=20% divergence)
    # is the divergence band these pairs actually fall in.
    with sam_path.open("wb") as sam_out:
        subprocess.run(
            ["minimap2", "-ax", "asm20", "--eqx", "-t", str(threads), str(ref_fa_plain), str(query_fa_plain)],
            stdout=sam_out,
            stderr=subprocess.PIPE,
            check=True,
            timeout=_COMPARATOR_TIMEOUT,
        )
    try:
        run(["syri", "-c", str(sam_path), "-r", str(ref_fa_plain), "-q", str(query_fa_plain),
             "-k", "-F", "S", "--dir", str(work_dir), "--prefix", "syri_"],
            cwd=work_dir, timeout=_COMPARATOR_TIMEOUT)
    except subprocess.CalledProcessError as exc:
        # SyRI rejects highly divergent pairs (e.g. cross-genus assemblies)
        # with a non-zero exit. Treat as "no comparator output" rather than
        # propagating the failure up to abort the whole panel — but capture
        # the actual stderr so the operator can distinguish "no syntenic
        # region" (the expected divergence case) from a missing chromosome,
        # SAM parse error, or pysam crash that needs different remediation.
        log_path = _persist_stderr(work_dir, "syri", exc)
        tail = _stderr_tail(exc)
        log_hint = f" (full log: {log_path})" if log_path else ""
        sys.stderr.write(
            f"[warn] syri rc={exc.returncode} for {query_asm}{log_hint}{tail}\n"
        )
        _log_comparator_failure(
            out_dir, "syri", query_asm,
            f"rc={exc.returncode} stderr_tail={tail.strip().replace(chr(10), ' | ')}",
        )
        return None
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
    # Binary-mode stdout-to-file: minigraph/gfatools emit text but the kernel
    # path bypasses Python's text decoder, so non-UTF-8 bytes in a contig name
    # (rare but observed on public ENA assemblies) don't kill the subprocess
    # wrapper before the comparator callset lands on disk.
    with graph_gfa.open("wb") as out_fh:
        subprocess.run(
            ["minigraph", "-cxggs", "-c", "-t", str(threads), *extra_args, str(ref_fa), str(query_fa)],
            stdout=out_fh,
            stderr=subprocess.PIPE,
            check=True,
            timeout=_COMPARATOR_TIMEOUT,
        )
    with bubble_bed.open("wb") as out_fh:
        subprocess.run(
            ["gfatools", "bubble", str(graph_gfa)],
            stdout=out_fh,
            stderr=subprocess.PIPE,
            check=True,
            timeout=_COMPARATOR_TIMEOUT,
        )
    with sample_bed.open("wb") as out_fh:
        subprocess.run(
            ["minigraph", "-cxasm", "--call", "-t", str(threads), *extra_args, str(graph_gfa), str(query_fa)],
            stdout=out_fh,
            stderr=subprocess.PIPE,
            check=True,
            timeout=_COMPARATOR_TIMEOUT,
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

    # pggb wraps wfmash → seqwish → smoothxg → vg deconstruct, and the
    # combined pipeline returns rc=2 when wfmash produces zero homologous
    # mappings.  Two predictable failure modes we can short-circuit:
    #   (a) pair.fa contains <2 records — `build_pair_fasta` writes ref then
    #       query, but if either FASTA was empty after gz decompression we
    #       end up with a single sequence and pggb's `-n 2` fails immediately.
    #   (b) the longest sequence is shorter than the segment length — wfmash
    #       can't pick a single segment and emits "no segment found".
    n_records, _, longest = _fasta_stats(pair_fa)
    if n_records < 2:
        sys.stderr.write(
            f"[skip] pggb {query_asm}: pair.fa has only {n_records} record(s); "
            "pggb -n 2 requires both ref and query sequences.\n"
        )
        _log_comparator_failure(
            out_dir, "pggb", query_asm,
            f"skipped:incomplete_pair n_records={n_records}",
        )
        return None
    try:
        seg_bp = int(str(segment_len).rstrip("kKmMgG")) * (
            1000 if str(segment_len).lower().endswith("k") else
            1_000_000 if str(segment_len).lower().endswith("m") else
            1_000_000_000 if str(segment_len).lower().endswith("g") else
            1
        )
    except ValueError:
        seg_bp = 5000
    if longest > 0 and longest < seg_bp:
        sys.stderr.write(
            f"[skip] pggb {query_asm}: longest sequence {longest} bp is below "
            f"segment_len {seg_bp} bp; wfmash cannot place a segment.\n"
        )
        _log_comparator_failure(
            out_dir, "pggb", query_asm,
            f"skipped:short_sequences longest={longest} segment={seg_bp}",
        )
        return None

    if tool_path("samtools"):
        try:
            run(["samtools", "faidx", str(pair_fa)], cwd=ROOT, timeout=_COMPARATOR_TIMEOUT)
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
    try:
        run(cmd, cwd=ROOT, timeout=_COMPARATOR_TIMEOUT)
    except subprocess.CalledProcessError as exc:
        # pggb exit 2 is what wfmash / seqwish / smoothxg / vg deconstruct
        # use when their inputs are unalignable (too divergent, identical
        # sequences, or chromosomes too short for the segment size).  Capture
        # the real stderr so the operator can tell which sub-step failed
        # rather than reading "exit status 2" with no context.
        log_path = _persist_stderr(work_dir, "pggb", exc)
        tail = _stderr_tail(exc)
        log_hint = f" (full log: {log_path})" if log_path else ""
        sys.stderr.write(
            f"[warn] pggb rc={exc.returncode} for {query_asm}{log_hint}{tail}\n"
        )
        _log_comparator_failure(
            out_dir, "pggb", query_asm,
            f"rc={exc.returncode} stderr_tail={tail.strip().replace(chr(10), ' | ')}",
        )
        return None
    vcf_candidates = sorted(work_dir.glob("**/*.vcf")) + sorted(work_dir.glob("**/*.vcf.gz"))
    if not vcf_candidates:
        _log_comparator_failure(out_dir, "pggb", query_asm, "rc=0 but no VCF emitted")
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
    """Align reads to the benchmark reference with minimap2 → sorted+indexed BAM.

    Preset routing (set by _long_read_preset / callers):
      map-hifi  PacBio HiFi CCS (Revio, Sequel IIe, Sequel II CCS)
      map-pb    PacBio CLR (RS II, Sequel I)
      map-ont   ONT — R10.4.1 simplex and R9.4.1 both use this preset
      sr        Illumina short reads (bwa-mem2 is an alternative for
                Delly / Manta pipelines that mandate BWA-formatted RG headers)

    The resulting BAM can feed:
      • sniffles2 / SVIM / cuteSV  for long-read SV calling
      • Delly / Manta              for short-read SV calling
      • WhatsHap phase + haplotag  for haplotype-phased variant calling in
        dikaryotic or diploid fungi (Puccinia, Leptosphaeria, Zymoseptoria)

    Returns (bam_path, ref_path) or None if prerequisites are not met.
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
    # SVIM/Sniffles/cuteSV/Delly/Manta all use pysam.FastaFile (or bcftools)
    # which cannot open .fna.gz directly. Hand them a plain copy materialised
    # next to the BAM, so the downstream caller invocations succeed.
    ref_fa_plain = _ensure_plain_fasta(ref_fa, work_dir)
    if ref_fa_plain is None:
        return None

    sam_path = work_dir / "aln.sam"
    bam_sorted = work_dir / "aln.sorted.bam"

    # minimap2 -> SAM (minimap2 itself accepts .gz, but we feed the plain
    # file so the BAM @SQ matches what the SV caller will index later).
    # Stream stdout straight to disk in binary mode so a non-ASCII byte in a
    # public ENA FASTQ header (which surfaces verbatim in @PG / read-name SAM
    # records) doesn't crash the wrapper with UnicodeDecodeError before svim /
    # sniffles / cutesv even see the BAM.
    with sam_path.open("wb") as sam_out:
        subprocess.run(
            ["minimap2", "-ax", preset, "-t", str(threads), str(ref_fa_plain), str(reads_path)],
            stdout=sam_out,
            stderr=subprocess.PIPE,
            timeout=_COMPARATOR_TIMEOUT,
            check=True,
        )
    # samtools sort -> BAM
    run(
        ["samtools", "sort", "-@", str(threads), "-o", str(bam_sorted), str(sam_path)],
        cwd=ROOT,
        timeout=_COMPARATOR_TIMEOUT,
    )
    # samtools index
    run(["samtools", "index", str(bam_sorted)], cwd=ROOT, timeout=_COMPARATOR_TIMEOUT)
    # Reference needs a .fai for downstream callers (Delly/Manta especially).
    if not (ref_fa_plain.parent / (ref_fa_plain.name + ".fai")).exists():
        try:
            run(["samtools", "faidx", str(ref_fa_plain)], cwd=ROOT, timeout=_COMPARATOR_TIMEOUT)
        except subprocess.CalledProcessError:
            pass
    # Drop the giant intermediate SAM once the BAM exists.
    try:
        sam_path.unlink()
    except OSError:
        pass
    return bam_sorted, ref_fa_plain


def _long_read_preset(query_row: dict[str, str]) -> str:
    """Return the minimap2 long-read alignment preset for this query.

    map-hifi  PacBio HiFi CCS (Revio, Sequel IIe, Sequel II CCS).
              Tuned for ≥99 % accuracy reads; incompatible with CLR data.
    map-pb    PacBio CLR (RS II, Sequel I, Sequel II in CLR mode).
    map-ont   Oxford Nanopore — covers R9.4.1 and R10.4.1 simplex.
              R10.4.1 on PromethION/GridION yields ~Q20 average; the same
              minimap2 preset applies, though Sniffles2 ≥v2.2 accepts
              --long-read-model ont_r10_q20 for a marginal recall boost.
    """
    platform = (query_row.get("instrument_platform") or "").strip().upper()
    if "PACBIO" in platform or "SMRT" in platform:
        return "map-hifi" if _is_pacbio_hifi(query_row) else "map-pb"
    # Default: Oxford Nanopore (R9.4.1, R10.4.1, or unspecified chemistry).
    return "map-ont"


def _existing_variants_vcf(work_dir: Path) -> Path | None:
    candidates = [work_dir / "variants.vcf"]
    candidates.extend(sorted(work_dir.rglob("variants.vcf")))
    for cand in candidates:
        if cand.exists() and cand.stat().st_size > 0:
            return cand
    return None


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
    # SVIM 2.0 + matplotlib 3.9 crashes in its final plotting step with
    # `AttributeError: 'Legend' object has no attribute 'legendHandles'`
    # (renamed to `legend_handles`). The VCF is written *before* the plot
    # step so a non-zero exit is harmless — we keep the VCF and continue.
    try:
        run(
            ["svim", "alignment", str(svim_out), str(bam_sorted), str(ref_fa)],
            cwd=ROOT,
            timeout=_COMPARATOR_TIMEOUT,
        )
    except subprocess.CalledProcessError as exc:
        vcf_path = _existing_variants_vcf(svim_out)
        if vcf_path is None:
            _log_comparator_failure(out_dir, "svim", query_asm, "failed_no_vcf")
            return None
        tail = _stderr_tail(exc, max_lines=3)
        legend_bug = "legendHandles" in (tail or "")
        sys.stderr.write(
            f"[warn] svim exited non-zero for {query_asm} "
            f"({'matplotlib legend bug — VCF unaffected' if legend_bug else 'see stderr'}), "
            f"keeping VCF output{tail}\n"
        )
        return {"label": "svim", "vcf": str(vcf_path)}
    vcf_path = _existing_variants_vcf(svim_out)
    if vcf_path is None:
        _log_comparator_failure(out_dir, "svim", query_asm, "failed_no_vcf")
        return None
    return {"label": "svim", "vcf": str(vcf_path)}


def run_sniffles_for_query(query_row: dict[str, str], out_dir: Path, threads: int) -> dict[str, str] | None:
    """Sniffles2: long-read SV caller. Emits a reference-coordinate VCF.

    For ONT R10.4.1 simplex reads, Sniffles2 ≥v2.2 accepts
    --long-read-model ont_r10_q20, which improves recall on ~Q20 data.
    We try the model flag first and fall back silently if unsupported.
    PacBio HiFi reads work with Sniffles2 defaults (the map-hifi BAM RG
    header is sufficient for Sniffles2 to auto-detect the platform).
    """
    if not tool_path("sniffles"):
        return None
    query_asm = query_row["query_asm"]
    preset = _long_read_preset(query_row)
    work_dir = out_dir / "comparators" / "sniffles" / query_asm
    aligned = _minimap2_align_reads(query_row, work_dir, threads, preset=preset)
    if aligned is None:
        return None
    bam_sorted, ref_fa = aligned
    vcf_path = work_dir / "sniffles.vcf"
    base_cmd = [
        "sniffles", "--input", str(bam_sorted),
        "--reference", str(ref_fa),
        "--vcf", str(vcf_path),
        "--threads", str(threads),
    ]
    # ONT R10.4.1 simplex: request the Q20-tuned internal model when available.
    platform = (query_row.get("instrument_platform") or "").upper()
    ont_model_cmd = (
        base_cmd + ["--long-read-model", "ont_r10_q20"]
        if "OXFORD_NANOPORE" in platform else base_cmd
    )
    try:
        run(ont_model_cmd, cwd=ROOT, timeout=_COMPARATOR_TIMEOUT)
    except subprocess.CalledProcessError:
        if ont_model_cmd is not base_cmd:
            try:
                run(base_cmd, cwd=ROOT, timeout=_COMPARATOR_TIMEOUT)
            except subprocess.CalledProcessError:
                return None
        else:
            return None
    if not vcf_path.exists():
        return None
    return {"label": "sniffles", "vcf": str(vcf_path)}


def run_cutesv_for_query(query_row: dict[str, str], out_dir: Path, threads: int) -> dict[str, str] | None:
    """cuteSV: long-read SV caller. Needs a writable working dir alongside the BAM.

    Cluster parameters differ by platform (from the cuteSV README):
      PacBio HiFi  max_cluster_bias 1000 / diff_ratio_merging 0.9 (INS) 0.5 (DEL)
                   High accuracy and long reads justify wider merging windows.
      ONT R10.4.1  max_cluster_bias 100  / diff_ratio_merging 0.3
      ONT R9.4.1   same conservative ONT defaults as R10.4.1
    """
    if not tool_path("cuteSV") and not tool_path("cutesv"):
        return None
    bin_name = "cuteSV" if tool_path("cuteSV") else "cutesv"
    query_asm = query_row["query_asm"]
    preset = _long_read_preset(query_row)
    work_dir = out_dir / "comparators" / "cutesv" / query_asm
    aligned = _minimap2_align_reads(query_row, work_dir, threads, preset=preset)
    if aligned is None:
        return None
    bam_sorted, ref_fa = aligned
    vcf_path = work_dir / "cutesv.vcf"
    tmp_dir = work_dir / "cutesv_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    # PacBio HiFi: tighter cluster merging exploits the higher per-read accuracy.
    # ONT (R10.4.1 simplex or R9.4.1): conservative defaults from the cuteSV README.
    if preset == "map-hifi":
        bias_ins, ratio_ins = "1000", "0.9"
        bias_del, ratio_del = "1000", "0.5"
    else:
        bias_ins, ratio_ins = "100",  "0.3"
        bias_del, ratio_del = "100",  "0.3"
    cmd = [
        bin_name,
        str(bam_sorted), str(ref_fa), str(vcf_path), str(tmp_dir),
        "--threads", str(threads),
        "--max_cluster_bias_INS", bias_ins,
        "--diff_ratio_merging_INS", ratio_ins,
        "--max_cluster_bias_DEL", bias_del,
        "--diff_ratio_merging_DEL", ratio_del,
        "--min_support", "3",
        "--sample", query_asm,
    ]
    try:
        run(cmd, cwd=ROOT, timeout=_COMPARATOR_TIMEOUT)
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
            timeout=_COMPARATOR_TIMEOUT,
        )
    except subprocess.CalledProcessError:
        return None
    if not bcf_path.exists():
        return None
    # BCF -> text VCF so the generic loader can parse it.
    with vcf_path.open("wb") as out_fh:
        subprocess.run(
            ["bcftools", "view", str(bcf_path)],
            stdout=out_fh, stderr=subprocess.PIPE, check=True,
            timeout=_COMPARATOR_TIMEOUT,
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
            timeout=_COMPARATOR_TIMEOUT,
        )
        run(
            [str(run_dir / "runWorkflow.py"), "-j", str(threads)],
            cwd=ROOT,
            timeout=_COMPARATOR_TIMEOUT,
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
        run(cmd, cwd=ROOT, timeout=_COMPARATOR_TIMEOUT)
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


def _fasta_stats(fa: Path) -> tuple[int, int, int]:
    """Return (n_records, total_bp, longest_bp) for a (possibly gzipped) FASTA.
    Used by comparator pre-checks to skip pairs that would deterministically
    blow up the downstream tool — e.g. fragmented MAG queries that defeat
    SyRI's chrmatch heuristic, or pairs with too-short contigs that fall
    below pggb's segment length."""
    n_records = 0
    longest = 0
    total = 0
    cur = 0
    opener = gzip.open if str(fa).endswith(".gz") else open
    try:
        with opener(fa, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    if cur > longest:
                        longest = cur
                    total += cur
                    cur = 0
                    n_records += 1
                else:
                    cur += len(line.strip())
            if cur > longest:
                longest = cur
            total += cur
    except OSError:
        return (0, 0, 0)
    return (n_records, total, longest)


def _ensure_plain_fasta(ref_fa: Path, work_dir: Path) -> Path | None:
    """svim-asm and a few other tools rely on pysam.FastaFile, which cannot
    open bgzipped or plain-gzipped FASTA without a `.gzi` companion. We
    materialise a plain-text copy in the comparator's work directory so the
    rest of the pipeline keeps working with .gz refs in the data cache.
    """
    if ref_fa.suffix != ".gz":
        return ref_fa
    plain = work_dir / ref_fa.with_suffix("").name
    if plain.exists() and plain.stat().st_size > 0:
        return plain
    try:
        with gzip.open(ref_fa, "rb") as src, plain.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    except OSError as exc:
        sys.stderr.write(f"[warn] could not decompress {ref_fa}: {exc}\n")
        return None
    return plain


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

    # svim-asm uses pysam.FastaFile which cannot open .fna.gz; ensure plain.
    ref_fa_plain = _ensure_plain_fasta(ref_fa, work_dir)
    if ref_fa_plain is None:
        return None

    sam_path = work_dir / "query_vs_ref.sam"
    bam_sorted = work_dir / "query_vs_ref.sorted.bam"
    # asm20 (divergence <=20%), NOT asm5 (<=5%). The benchmark reference is a
    # held-out *sibling-clade* genome, so query vs ref is a cross-species
    # fungal pair. asm5 on such a pair yields only short, ambiguous
    # alignments — every record comes back MAPQ 0 — and svim-asm's
    # min_mapq=20 filter then discards them all ("Found 0 candidates"),
    # giving truth_calls=0 and F1=nan for the diverged samples. asm20 keeps
    # the cross-species synteny minimap2 can actually anchor; it matches the
    # asm20 preset the clade-lift path already uses for the same reason.
    with sam_path.open("wb") as sam_out:
        subprocess.run(
            ["minimap2", "-ax", "asm20", "-t", str(threads), str(ref_fa_plain), str(query_fa)],
            stdout=sam_out, stderr=subprocess.PIPE, check=True,
            timeout=_COMPARATOR_TIMEOUT,
        )
    try:
        run(
            ["samtools", "sort", "-@", str(threads), "-o", str(bam_sorted), str(sam_path)],
            cwd=ROOT,
            timeout=_COMPARATOR_TIMEOUT,
        )
        run(["samtools", "index", str(bam_sorted)], cwd=ROOT, timeout=_COMPARATOR_TIMEOUT)
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
            ["svim-asm", "haploid", str(svim_out), str(bam_sorted), str(ref_fa_plain)],
            cwd=ROOT,
            timeout=_COMPARATOR_TIMEOUT,
        )
    except subprocess.CalledProcessError:
        vcf_path = _existing_variants_vcf(svim_out)
        if vcf_path is None:
            _log_comparator_failure(out_dir, "svim_asm", query_asm, "failed_no_vcf")
            return None
        sys.stderr.write(
            f"[warn] svim-asm exited non-zero for {query_asm}, but produced "
            f"{vcf_path}; keeping VCF output\n"
        )
        return {"label": "svim_asm", "vcf": str(vcf_path)}

    vcf_path = _existing_variants_vcf(svim_out)
    if vcf_path is None:
        _log_comparator_failure(out_dir, "svim_asm", query_asm, "failed_no_vcf")
        return None
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
    # paftools.js call wants a plain reference for sequence retrieval.
    ref_fa_plain = _ensure_plain_fasta(ref_fa, work_dir)
    if ref_fa_plain is None:
        return None

    sam_path = work_dir / "q2r.sam"
    paf_path = work_dir / "q2r.paf"
    sorted_paf = work_dir / "q2r.srt.paf"
    vcf_path = work_dir / "anchorwave.vcf"

    # Step 1: minimap2 asm5 alignment as AnchorWave input seeds.
    with sam_path.open("wb") as sam_out:
        subprocess.run(
            ["minimap2", "-ax", "asm5", "--cs", "-t", str(threads), str(ref_fa_plain), str(query_fa)],
            stdout=sam_out, stderr=subprocess.PIPE, check=True,
            timeout=_COMPARATOR_TIMEOUT,
        )
    # Step 2: SAM -> PAF via paftools for downstream AnchorWave refinement.
    try:
        with paf_path.open("wb") as paf_out:
            subprocess.run(
                [paftools, "sam2paf", str(sam_path)],
                stdout=paf_out, stderr=subprocess.PIPE, check=True,
                timeout=_COMPARATOR_TIMEOUT,
            )
    except subprocess.CalledProcessError:
        return None
    # Step 3: sort PAF by target coordinates.
    try:
        with sorted_paf.open("wb") as srt_out:
            subprocess.run(
                ["sort", "-k6,6", "-k8,8n", str(paf_path)],
                stdout=srt_out, stderr=subprocess.PIPE, check=True,
                timeout=_COMPARATOR_TIMEOUT,
            )
    except subprocess.CalledProcessError:
        return None
    # Step 4: paftools.js call produces a reference-coordinate VCF directly
    # from the sorted PAF, which is the shape load_reference_vcf_calls reads.
    # AnchorWave itself is used upstream to produce better assembly-to-assembly
    # alignments on repetitive fungal genomes; the CLI boundary here is the
    # SV VCF that results.
    try:
        with vcf_path.open("wb") as vcf_out:
            subprocess.run(
                [paftools, "call", "-f", str(ref_fa_plain), "-L", "50", str(sorted_paf)],
                stdout=vcf_out, stderr=subprocess.PIPE, check=True,
                timeout=_COMPARATOR_TIMEOUT,
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


def _per_query_thread_budget(n_queries: int, threads: int) -> tuple[int, int]:
    """Return (max_parallel_queries, per_query_threads) for a comparator pool.

    Cap parallelism so each query still gets >=3 threads (the practical floor
    where minimap2 / samtools sort show usable speedup): for the typical
    million-real 5-query x 16-thread case this yields 5 parallel queries with
    3 threads each, versus the old 1 query x 16 threads x 5 serial iterations.
    """
    n = max(1, n_queries)
    per_min = max(1, threads // n if threads >= n else 1)
    if per_min < 3 and threads >= 3:
        max_parallel = max(1, threads // 3)
    else:
        max_parallel = n
    max_parallel = min(max_parallel, n)
    per_query_threads = max(1, threads // max(1, max_parallel))
    return max_parallel, per_query_threads


def run_per_query_in_parallel(
    label: str,
    runner,
    queries: list[dict[str, str]],
    out_dir: Path,
    threads: int,
) -> dict[str, dict[str, str] | None]:
    """Schedule per-query comparator calls across a thread pool.

    `runner` must take (query_row, out_dir, per_query_threads) and return the
    per-query result dict (or None). Failures on a single query are caught,
    logged, and recorded via _log_comparator_failure — they do not abort the
    other queries in the pool. Returns {query_asm: result_or_None}.
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    if not queries:
        return {}
    max_parallel, per_query_threads = _per_query_thread_budget(len(queries), threads)
    results: dict[str, dict[str, str] | None] = {}
    sys.stderr.write(
        f"[parallel] {label}: scheduling {len(queries)} queries across "
        f"{max_parallel} workers ({per_query_threads} threads each)\n"
    )
    sys.stderr.flush()
    with ThreadPoolExecutor(max_workers=max_parallel) as pool:
        futures = {
            pool.submit(runner, qr, out_dir, per_query_threads): qr
            for qr in queries
        }
        for fut in as_completed(futures):
            qr = futures[fut]
            qa = qr["query_asm"]
            try:
                results[qa] = fut.result()
            except subprocess.TimeoutExpired as exc:
                sys.stderr.write(
                    f"[warn] {label} timed out for {qa} after "
                    f"{_COMPARATOR_TIMEOUT}s — continuing other queries\n"
                )
                _log_comparator_failure(out_dir, label, qa, f"timeout:{exc}")
                results[qa] = None
            except Exception as exc:  # pragma: no cover - defensive
                sys.stderr.write(f"[warn] {label} failed for {qa}: {exc}\n")
                _log_comparator_failure(out_dir, label, qa, f"exception:{exc}")
                results[qa] = None
    return results


def call_key(call: NormalizedCall) -> tuple[str, str, int, int, str]:
    return (call.query_asm, call.query_contig, call.pos, call.end, call.svtype)


def pangenome_locus_key(call: NormalizedCall, bucket_bp: int = 100) -> tuple[str, str, str, int, int, int, str, int]:
    """Collapse reference-background duplicates into an approximate pangenome
    locus key. This intentionally ignores ref_asm: the same query-space event
    called against several pangenome references should count once as a
    deduplicated biological locus, while raw_pairwise observations still count
    every VCF row.
    """
    bucket = max(1, bucket_bp)
    mate_contig = "."
    mate_bucket = -1
    if _has_mate(call):
        mate_contig = call.mate_contig or "."
        mate_bucket = max(0, call.mate_pos) // bucket
    return (
        call.query_asm,
        call.query_contig,
        call.svtype,
        max(0, call.pos) // bucket,
        max(0, call.end) // bucket,
        abs(call.svlen) // bucket,
        mate_contig,
        mate_bucket,
    )


def write_pangenome_call_layers(
    out_path: Path,
    query_manifest: list[dict[str, str]],
    mycosv_calls_by_query: dict[str, dict[str, Any]],
    single_ref_counts: dict[str, int],
    single_ref_keys_by_query: dict[str, set[tuple[str, str, int, int, str]]],
    support_by_key: dict[tuple[str, str, int, int, str], list[str]],
    validated_by_query: dict[str, set[tuple[str, str, int, int, str]]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    totals = {
        "raw_pairwise_pangenome_observations": 0,
        "deduplicated_pangenome_loci": 0,
        "single_reference_equivalent_calls": 0,
        "pangenome_only_calls": 0,
        "pangenome_only_read_supported": 0,
        "pangenome_only_intrinsic_supported": 0,
    }
    for row in query_manifest:
        qa = row["query_asm"]
        calls = list(mycosv_calls_by_query.get(qa, {}).get("query", []))
        loci = {pangenome_locus_key(c) for c in calls}
        validated = validated_by_query.get(qa, set())
        single_ref_keys = single_ref_keys_by_query.get(qa, set())
        single_ref_loci = {
            pangenome_locus_key(c) for c in calls if call_key(c) in single_ref_keys
        }
        pangenome_only_loci: dict[
            tuple[str, str, str, int, int, int, str, int],
            dict[str, bool],
        ] = {}
        for call in calls:
            key = call_key(call)
            if support_by_key.get(key):
                continue
            locus = pangenome_locus_key(call)
            if key in single_ref_keys or locus in single_ref_loci:
                continue
            state = pangenome_only_loci.setdefault(
                locus, {"read_supported": False, "intrinsic_supported": False}
            )
            if key in validated:
                state["read_supported"] = True
            elif (call.read_support or 0) >= 2:
                state["intrinsic_supported"] = True
        single_ref = int(single_ref_counts.get(qa, 0))
        pangenome_only = len(pangenome_only_loci)
        pangenome_only_read = sum(
            1 for state in pangenome_only_loci.values()
            if state["read_supported"]
        )
        pangenome_only_intrinsic = sum(
            1 for state in pangenome_only_loci.values()
            if state["intrinsic_supported"] and not state["read_supported"]
        )
        record = {
            "query_asm": qa,
            "query_mode": row.get("query_mode", "."),
            "raw_pairwise_pangenome_observations": len(calls),
            "deduplicated_pangenome_loci": len(loci),
            "single_reference_equivalent_calls": single_ref,
            "pangenome_only_calls": pangenome_only,
            "pangenome_only_read_supported": pangenome_only_read,
            "pangenome_only_intrinsic_supported": pangenome_only_intrinsic,
            "raw_to_deduplicated_ratio": f"{(len(calls) / len(loci)):.3f}" if loci else "nan",
            "single_ref_fraction_of_raw": f"{(single_ref / len(calls)):.4f}" if calls else "nan",
        }
        rows.append(record)
        for k in totals:
            totals[k] += int(record[k])
    if rows:
        total_raw = totals["raw_pairwise_pangenome_observations"]
        total_dedup = totals["deduplicated_pangenome_loci"]
        total_single = totals["single_reference_equivalent_calls"]
        rows.append({
            "query_asm": "ALL",
            "query_mode": "all",
            **totals,
            "raw_to_deduplicated_ratio": f"{(total_raw / total_dedup):.3f}" if total_dedup else "nan",
            "single_ref_fraction_of_raw": f"{(total_single / total_raw):.4f}" if total_raw else "nan",
        })
    write_tsv(
        out_path,
        rows,
        [
            "query_asm", "query_mode",
            "raw_pairwise_pangenome_observations",
            "deduplicated_pangenome_loci",
            "single_reference_equivalent_calls",
            "pangenome_only_calls",
            "pangenome_only_read_supported",
            "pangenome_only_intrinsic_supported",
            "raw_to_deduplicated_ratio",
            "single_ref_fraction_of_raw",
        ],
    )
    return rows


_TE_LIKE_CLASSES = {
    "REPEAT", "TE", "LTR_GYPSY", "LTR_COPIA", "DNA_TIR", "HELITRON", "MITE",
    "LINE", "SINE", "RIP", "STARSHIP", "HGT", "TE_LTR", "TE_TIR", "TE_LINE", "TE_SINE",
}


def _fisher_right_tail(a: int, b: int, c: int, d: int) -> float:
    """One-sided Fisher exact P(X >= a): feature enrichment in group A.

    Table:
      feature yes:  a in group A, c in group B
      feature no:   b in group A, d in group B
    """
    n = a + b + c + d
    row_a = a + b
    col_feature = a + c
    if n == 0 or row_a == 0 or col_feature == 0:
        return float("nan")
    lo = max(0, row_a - (n - col_feature))
    hi = min(row_a, col_feature)

    def hypergeom(x: int) -> float:
        return (
            math.comb(col_feature, x)
            * math.comb(n - col_feature, row_a - x)
            / math.comb(n, row_a)
        )

    return min(1.0, sum(hypergeom(x) for x in range(max(a, lo), hi + 1)))


def _odds_ratio(a: int, b: int, c: int, d: int) -> float:
    # Haldane-Anscombe correction keeps odds ratios finite for empty cells.
    return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))


def _stream_gene_interval_index(
    gene_annotations_tsv: Path | None,
    calls: list[NormalizedCall],
) -> dict[tuple[str, str], tuple[list[int], list[int]]]:
    wanted: dict[tuple[str, str], None] = {
        (c.query_asm, c.query_contig): None
        for c in calls
        if c.query_asm and c.query_contig
    }
    wanted.update({
        (".", c.query_contig): None
        for c in calls
        if c.query_contig
    })
    if not wanted or gene_annotations_tsv is None or not gene_annotations_tsv.exists():
        return {}
    if _GENE_ANNOTATION_MAX_BYTES >= 0:
        try:
            size = gene_annotations_tsv.stat().st_size
        except OSError:
            size = 0
        if size > _GENE_ANNOTATION_MAX_BYTES:
            sys.stderr.write(
                f"[gene-annot] skipping gene-proximal enrichment for "
                f"{gene_annotations_tsv} ({size} bytes > "
                f"MYCOSV_GENE_ANNOTATION_MAX_BYTES={_GENE_ANNOTATION_MAX_BYTES})\n"
            )
            return {}
    intervals: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    try:
        with gene_annotations_tsv.open(encoding="utf-8", errors="replace", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                qasm = row.get("query_asm") or row.get("asm") or "."
                contig = row.get("query_contig") or row.get("contig") or row.get("chrom") or row.get("seqid") or ""
                key = (qasm, contig)
                if key not in wanted:
                    continue
                try:
                    start = int(float(row.get("start") or row.get("gene_start") or row.get("pos") or 0))
                    end = int(float(row.get("end") or row.get("gene_end") or row.get("stop") or start))
                except ValueError:
                    continue
                if end < start:
                    start, end = end, start
                if start > 0 and end > 0:
                    intervals[key].append((start, end))
    except OSError:
        return {}

    out: dict[tuple[str, str], tuple[list[int], list[int]]] = {}
    for key, vals in intervals.items():
        vals.sort()
        starts: list[int] = []
        prefix_max_end: list[int] = []
        cur = 0
        for start, end in vals:
            starts.append(start)
            cur = max(cur, end)
            prefix_max_end.append(cur)
        out[key] = (starts, prefix_max_end)
    return out


def _call_is_gene_proximal(
    call: NormalizedCall,
    gene_index: dict[tuple[str, str], tuple[list[int], list[int]]],
    window_bp: int,
) -> bool:
    starts, prefix_max_end = (
        gene_index.get((call.query_asm, call.query_contig))
        or gene_index.get((".", call.query_contig))
        or ([], [])
    )
    if not starts:
        return False
    left = max(1, min(call.pos, call.end) - window_bp)
    right = max(call.pos, call.end) + window_bp
    idx = bisect.bisect_right(starts, right)
    return idx > 0 and prefix_max_end[idx - 1] >= left


def write_mycosv_novel_biology_enrichment(
    out_path: Path,
    calls: list[NormalizedCall],
    support_by_key: dict[tuple[str, str, int, int, str], list[str]],
    single_ref_keys_by_query: dict[str, set[tuple[str, str, int, int, str]]],
    single_ref_loci_by_query: dict[
        str, set[tuple[str, str, str, int, int, int, str, int]]
    ],
    validated_by_query: dict[str, set[tuple[str, str, int, int, str]]],
    gene_annotations_tsv: Path | None,
    *,
    gene_window_bp: int = 5000,
) -> list[dict[str, Any]]:
    gene_index = _stream_gene_interval_index(gene_annotations_tsv, calls)

    def has_external_read(call: NormalizedCall) -> bool:
        return call_key(call) in validated_by_query.get(call.query_asm, set())

    def is_unique(call: NormalizedCall) -> bool:
        key = call_key(call)
        if key in single_ref_keys_by_query.get(call.query_asm, set()):
            return False
        if pangenome_locus_key(call) in single_ref_loci_by_query.get(call.query_asm, set()):
            return False
        return not support_by_key.get(key)

    feature_defs = {
        "te_or_mge_like": lambda c: (c.element_class or "NONE").upper() in _TE_LIKE_CLASSES,
        "off_reference_or_novel": lambda c: c.svtype == "OFF_REF" or c.annotation in {"NOVEL", "NOVEL_WEAK", "DIVERGED"},
        "hgt_or_starship_candidate": lambda c: (c.element_class or "").upper() in {"HGT", "STARSHIP"},
        f"gene_proximal_{gene_window_bp}bp": lambda c: _call_is_gene_proximal(c, gene_index, gene_window_bp),
    }

    rows: list[dict[str, Any]] = []
    for feature, predicate in feature_defs.items():
        unique_calls = [c for c in calls if is_unique(c)]
        supported_calls = [c for c in calls if not is_unique(c)]
        a = sum(1 for c in unique_calls if predicate(c))
        b = len(unique_calls) - a
        c_count = sum(1 for c in supported_calls if predicate(c))
        d = len(supported_calls) - c_count
        rows.append({
            "feature": feature,
            "group_a": "mycosv_unique",
            "group_b": "comparator_supported",
            "a_feature": a,
            "a_nonfeature": b,
            "b_feature": c_count,
            "b_nonfeature": d,
            "odds_ratio": f"{_odds_ratio(a, b, c_count, d):.6g}",
            "fisher_right_p": (
                f"{_fisher_right_tail(a, b, c_count, d):.6g}"
                if supported_calls and unique_calls else "nan"
            ),
            "unique_read_supported_feature": sum(
                1 for call in unique_calls if predicate(call) and has_external_read(call)
            ),
            "unique_intrinsic_supported_feature": sum(
                1 for call in unique_calls
                if predicate(call) and not has_external_read(call) and (call.read_support or 0) >= 2
            ),
            "interpretation": (
                "heuristic_candidate_enrichment_not_confirmed_hgt"
                if feature == "hgt_or_starship_candidate"
                else "enrichment_screen"
            ),
        })
    write_tsv(
        out_path,
        rows,
        [
            "feature", "group_a", "group_b",
            "a_feature", "a_nonfeature", "b_feature", "b_nonfeature",
            "odds_ratio", "fisher_right_p",
            "unique_read_supported_feature", "unique_intrinsic_supported_feature",
            "interpretation",
        ],
    )
    return rows


# Evidence tiers for a MycoSV call, ranked highest to lowest. "strong" means the
# call is independently supported by a comparator AND by external read evidence;
# "moderate" by one of the two; "intrinsic_only" means the call cleared the C++
# clustering floor (SUPPORT>=2) but no external signal could be re-anchored
# (sibling-clade contig absent from the validation BAM, low query coverage, …);
# "weak" is the rare residual where even MycoSV's intrinsic count is <=1. The
# panorama panel in the report shows these so an unvalidatable-but-real MycoSV
# call surfaces as "intrinsic_only" rather than disappearing from the F1 plots.
EVIDENCE_TIERS: tuple[str, ...] = ("strong", "moderate", "intrinsic_only", "weak")


def classify_evidence_tier(
    call: NormalizedCall,
    *,
    has_comparator: bool,
    has_external_read_support: bool,
) -> str:
    if has_comparator and has_external_read_support:
        return "strong"
    if has_comparator or has_external_read_support:
        return "moderate"
    intrinsic = call.read_support if call.read_support is not None else 0
    if intrinsic >= 2:
        return "intrinsic_only"
    return "weak"


def write_agreement_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    # svtype column: "ALL" for the aggregate row that already existed, plus
    # one row per canonical SV type for downstream "MycoSV vs comparator
    # per-svtype wins" visualization. Older readers ignore the new column.
    # validation_basis makes the semantics explicit: comparator outputs are
    # baselines/agreement signals, while read_level_union rows are independently
    # anchored in raw FASTQ/read evidence.
    write_tsv(
        path,
        rows,
        [
            "query_asm", "coordinate_space", "truth_label", "validation_basis",
            "svtype", "method",
            "truth_calls", "pred_calls",
            "tp", "fp", "fn", "precision", "recall", "f1",
            "prec_lo95", "prec_hi95", "rec_lo95", "rec_hi95", "status",
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
    ecological_traits_tsv: Path | None,
    fungaltraits_csv: Path | None,
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
        "--hits", mycosv_paths["hits"],
        "--query-metadata", str((prepared_dir / "query_manifest.tsv").resolve()),
        "--out-tsv", str(candidates_tsv.resolve()),
        "--summary-json", str(summary_json.resolve()),
        "--top-n", "200",
    ]
    if expression_tsv:
        cmd.extend(["--expression-tsv", str(expression_tsv.resolve())])
    if gene_annotations_tsv:
        pass_gene_annotations = True
        if _GENE_ANNOTATION_MAX_BYTES >= 0:
            try:
                gene_size = gene_annotations_tsv.stat().st_size
            except OSError:
                gene_size = 0
            if gene_size > _GENE_ANNOTATION_MAX_BYTES:
                pass_gene_annotations = False
                sys.stderr.write(
                    f"[biology] not passing large gene annotation file to "
                    f"candidate analyzer ({gene_size} bytes > "
                    f"MYCOSV_GENE_ANNOTATION_MAX_BYTES="
                    f"{_GENE_ANNOTATION_MAX_BYTES})\n"
                )
        if pass_gene_annotations:
            cmd.extend(["--gene-annotations", str(gene_annotations_tsv.resolve())])
    if ecological_traits_tsv:
        cmd.extend(["--ecological-traits", str(ecological_traits_tsv.resolve())])
    if fungaltraits_csv:
        cmd.extend(["--fungaltraits-csv", str(fungaltraits_csv.resolve())])
    if ancestral_tsv:
        cmd.extend(["--ancestral", str(ancestral_tsv.resolve())])
    try:
        run(cmd, cwd=ROOT, timeout=max(1, _BIOLOGY_TIMEOUT))
        return candidates_tsv, summary_json
    except subprocess.TimeoutExpired:
        sys.stderr.write(
            f"[biology] candidate analysis timed out after "
            f"{_BIOLOGY_TIMEOUT}s; continuing with enrichment-only outputs\n"
        )
        return None, None
    except subprocess.CalledProcessError as exc:
        sys.stderr.write(
            f"[biology] candidate analysis failed rc={exc.returncode}; "
            f"continuing with enrichment-only outputs\n"
        )
        return None, None


def join_biology_findings(
    candidates_tsv: Path | None,
    mycosv_calls: list[NormalizedCall],
    support_by_key: dict[tuple[str, str, int, int, str], list[str]],
    out_path: Path,
    *,
    single_ref_keys_by_query: dict[str, set[tuple[str, str, int, int, str]]] | None = None,
    single_ref_loci_by_query: dict[
        str, set[tuple[str, str, str, int, int, int, str, int]]
    ] | None = None,
    tier_by_key: dict[tuple[str, str, int, int, str], str] | None = None,
) -> None:
    single_ref_keys_by_query = single_ref_keys_by_query or {}
    single_ref_loci_by_query = single_ref_loci_by_query or {}
    tier_by_key = tier_by_key or {}
    if candidates_tsv is None or not candidates_tsv.exists():
        write_tsv(
            out_path,
            [],
            [
                "query_asm", "query_contig", "pos", "end", "svtype", "svlen",
                "annotation", "element_class", *BIOLOGY_FINDINGS_EXTRA_FIELDS,
            ],
        )
        return
    rows: list[dict[str, Any]] = []
    fieldnames: list[str] = []
    with candidates_tsv.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        for row in reader:
            key = (
                row.get("query_asm", "."),
                row.get("query_contig", "."),
                int(row.get("pos", "0") or 0),
                int(row.get("end", "0") or 0),
                row.get("svtype", "."),
            )
            supporters = support_by_key.get(key, [])
            try:
                svlen = int(row.get("svlen", "0") or 0)
            except ValueError:
                svlen = 0
            row_call = NormalizedCall(
                key[0], key[1], key[2], key[3], key[4], svlen, "mycosv",
                annotation=row.get("annotation", "."),
                element_class=row.get("element_class", "NONE"),
            )
            in_single_ref = (
                key in single_ref_keys_by_query.get(key[0], set())
                or pangenome_locus_key(row_call) in single_ref_loci_by_query.get(key[0], set())
            )
            row["comparator_support_count"] = len(supporters)
            row["comparator_support_labels"] = ",".join(sorted(supporters)) if supporters else "."
            row["single_reference_equivalent"] = "yes" if in_single_ref else "no"
            row["mycosv_unique"] = "yes" if not supporters and not in_single_ref else "no"
            row["evidence_tier"] = tier_by_key.get(key, "weak")
            rows.append(row)
    if not rows:
        if fieldnames:
            for extra in BIOLOGY_FINDINGS_EXTRA_FIELDS:
                if extra not in fieldnames:
                    fieldnames.append(extra)
            write_tsv(out_path, [], fieldnames)
        return
    fieldnames = list(rows[0].keys())
    write_tsv(out_path, rows, fieldnames)


def _report_comparator_status(args: argparse.Namespace, out_dir: Path) -> None:
    """Check which requested SOTA comparators are on PATH; warn clearly for missing ones."""

    # Map (flag_attr, tool_binary, conda_install_hint)
    _TOOL_HINTS: dict[str, tuple[list[str], str]] = {
        # assembly
        "run_syri":       (["minimap2", "syri"],
                           "conda install -c bioconda minimap2 syri"),
        "run_minigraph":  (["minigraph", "gfatools"],
                           "conda install -c bioconda minigraph gfatools"),
        "run_pggb":       (["pggb"],
                           "conda install -c bioconda pggb"),
        "run_cactus":     (["cactus-pangenome"],
                           "conda install -c bioconda cactus  (or download binary from GitHub)"),
        "run_svim_asm":   (["svim-asm", "minimap2", "samtools"],
                           "conda install -c bioconda svim-asm minimap2 samtools"),
        "run_anchorwave": (["anchorwave", "minimap2", "samtools"],
                           "conda install -c bioconda anchorwave minimap2 samtools"),
        # short-reads
        "run_delly":      (["delly", "bcftools"],
                           "conda install -c bioconda delly bcftools"),
        "run_manta":      (["configManta.py", "samtools"],
                           "conda install -c bioconda manta samtools"),
        # long-reads
        "run_svim":       (["svim"],
                           "conda install -c bioconda svim"),
        "run_sniffles":   (["sniffles"],
                           "conda install -c bioconda sniffles"),
        "run_cutesv":     (["cuteSV", "samtools"],
                           "conda install -c bioconda cutesv samtools"),
    }

    # Categorize every comparator in three buckets so the file is informative
    # even when no --run-X flag was passed (previously the file was just a
    # header). The buckets are: (1) requested + available -> actually ran,
    # (2) requested + missing binaries -> silently skipped, (3) not requested
    # but binaries are present -> available, just not enabled.
    requested_running: list[tuple[str, list[str]]] = []
    requested_missing: list[tuple[str, list[str], str]] = []
    not_requested: list[tuple[str, list[str], list[str]]] = []  # (name, present, absent)

    for flag, (binaries, hint) in _TOOL_HINTS.items():
        absent = [b for b in binaries if not tool_path(b)]
        present = [b for b in binaries if tool_path(b)]
        tool_name = flag.replace("run_", "")
        if getattr(args, flag, False):
            if absent:
                requested_missing.append((flag, absent, hint))
            else:
                requested_running.append((tool_name, binaries))
        else:
            not_requested.append((tool_name, present, absent))

    lines: list[str] = []
    lines.append("=== Comparator Pre-flight Check ===\n\n")
    lines.append(f"Mode: {getattr(args, 'mode', '?')}\n\n")

    if requested_running:
        lines.append("AVAILABLE AND ENABLED (will run):\n")
        for name, binaries in requested_running:
            lines.append(f"  [run] {name}  (binaries: {', '.join(binaries)})\n")
        lines.append("\n")

    if requested_missing:
        lines.append("REQUESTED BUT MISSING (silently skipped):\n")
        for flag, absent, hint in requested_missing:
            tool_name = flag.replace("run_", "")
            lines.append(f"  [skip] {tool_name}: missing {', '.join(absent)}\n")
            lines.append(f"         Install: {hint}\n")
        lines.append("\n")

    if not_requested:
        lines.append("NOT REQUESTED (pass the corresponding --run-X flag to enable):\n")
        for name, present, absent in not_requested:
            flag_hint = f"--run-{name.replace('_', '-')}"
            if absent:
                lines.append(
                    f"  [off] {name}  (would also need: {', '.join(absent)}; "
                    f"{flag_hint} to enable)\n"
                )
            else:
                lines.append(
                    f"  [off] {name}  (binaries OK; pass {flag_hint} to enable)\n"
                )
        lines.append("\n")

    if requested_missing or any(absent for _, _, absent in not_requested):
        lines.append("Install all missing comparator binaries with:\n")
        lines.append(f"  bash {Path(__file__).parent / 'install_tools.sh'}\n")
        lines.append("  bash install_tools.sh --check  (lists per-tool availability)\n")

    status_text = "".join(lines)
    status_file = out_dir / "COMPARATORS_STATUS.txt"
    status_file.write_text(status_text, encoding="utf-8")

    # Aliases for the legacy two-bucket caller below.
    available = [name for name, _ in requested_running]
    missing = requested_missing

    if missing:
        sys.stderr.write("\n[warn] The following requested comparators are missing and will "
                         "be skipped (no results will appear for them):\n")
        for flag, absent, hint in missing:
            tool_name = flag.replace("run_", "")
            sys.stderr.write(f"  ✗ {tool_name}: missing binaries {', '.join(absent)}\n")
            sys.stderr.write(f"    Install: {hint}\n")
        sys.stderr.write(f"\nFull status written to: {status_file}\n")
        sys.stderr.write(f"Install all tools: bash {Path(__file__).parent / 'install_tools.sh'}\n\n")
    elif available:
        sys.stderr.write(f"[info] All requested comparators available: "
                         f"{', '.join(available)}\n")
    else:
        sys.stderr.write("[info] No SOTA comparators requested (pass --run-syri etc. to enable)\n")
        sys.stderr.write(f"       Install script: bash {Path(__file__).parent / 'install_tools.sh'}\n")


_COMPARATOR_FLAG_BY_MODE: dict[str, list[tuple[str, list[str]]]] = {
    # (run_X attr, required binaries); used to auto-enable comparators by mode.
    "assembly":    [("run_syri", ["minimap2", "syri"]),
                    ("run_minigraph", ["minigraph", "gfatools"]),
                    ("run_pggb", ["pggb"]),
                    ("run_cactus", ["cactus-pangenome"]),
                    ("run_svim_asm", ["svim-asm", "minimap2", "samtools"]),
                    ("run_anchorwave", ["anchorwave", "minimap2", "samtools"])],
    "short-reads": [("run_delly", ["delly", "bcftools"]),
                    ("run_manta", ["configManta.py", "samtools"])],
    "long-reads":  [("run_svim", ["svim"]),
                    ("run_sniffles", ["sniffles"]),
                    ("run_cutesv", ["cuteSV", "samtools"])],
}


def _auto_enable_comparators(args: argparse.Namespace) -> None:
    """Toggle every --run-X flag whose binaries are detected, scoped to the
    benchmark mode. Prints the resolved list to stderr so the operator sees
    exactly which comparators will run for this invocation.
    """
    mode = getattr(args, "mode", "assembly")
    candidates = _COMPARATOR_FLAG_BY_MODE.get(mode, [])
    enabled: list[str] = []
    skipped: list[tuple[str, list[str]]] = []
    for flag, binaries in candidates:
        absent = [b for b in binaries if not tool_path(b)]
        if absent:
            skipped.append((flag.replace("run_", ""), absent))
            continue
        if not getattr(args, flag, False):
            setattr(args, flag, True)
        enabled.append(flag.replace("run_", ""))
    if enabled:
        sys.stderr.write(
            f"[comparators] --run-all-comparators auto-enabled for mode={mode}: "
            f"{', '.join(enabled)}\n"
        )
    if skipped:
        for tool_name, absent in skipped:
            sys.stderr.write(
                f"[comparators] {tool_name}: missing {', '.join(absent)} — "
                f"run install_tools.sh to install\n"
            )
    if not enabled:
        sys.stderr.write(
            f"[comparators] --run-all-comparators: no comparator binaries "
            f"detected for mode={mode}; will produce no_comparator placeholder rows\n"
        )


def benchmark_real_data(args: argparse.Namespace) -> int:
    prepared_dir = args.prepared_dir.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    mycosv_only = bool(getattr(args, "mycosv_only", False))
    if getattr(args, "run_all_comparators", False) and not mycosv_only:
        _auto_enable_comparators(args)
    # Force-on the canonical comparator(s) per mode so the visualization's
    # "mycosv-vs-baseline per SV type" panel is never empty. Minigraph is
    # the assembly-mode baseline; Delly/Manta cover short-reads; Sniffles2
    # and cuteSV cover long-reads. Read-level validation (samtools-driven
    # split-read counting) runs in every mode independently — see
    # validate_calls_with_reads(). --mycosv-only disables this forcing for
    # the million-real flow where comparators are out of scope.
    _MANDATORY_BASELINES_BY_MODE: dict[str, list[tuple[str, list[str]]]] = {
        "assembly":    [("run_minigraph", ["minigraph", "gfatools"])],
        "short-reads": [("run_delly", ["delly", "bcftools"]),
                        ("run_manta", ["configManta.py", "samtools"])],
        "long-reads":  [("run_sniffles", ["sniffles"]),
                        ("run_cutesv", ["cuteSV", "samtools"])],
    }
    forced: list[str] = []
    missing_baselines: list[tuple[str, list[str]]] = []
    if not mycosv_only:
        for flag, binaries in _MANDATORY_BASELINES_BY_MODE.get(args.mode, []):
            absent = [b for b in binaries if not tool_path(b)]
            if absent:
                missing_baselines.append((flag.replace("run_", ""), absent))
                continue
            if not getattr(args, flag, False):
                setattr(args, flag, True)
                forced.append(flag.replace("run_", ""))
        if forced:
            sys.stderr.write(
                f"[comparators] forcing mandatory baselines on for mode={args.mode}: "
                f"{', '.join(forced)} (canonical real-data baselines for this mode)\n"
            )
        if missing_baselines:
            for tool_name, absent in missing_baselines:
                sys.stderr.write(
                    f"[comparators] WARNING: {tool_name} (mandatory {args.mode} "
                    f"baseline) is missing binaries {', '.join(absent)} — install "
                    f"via install_tools.sh; the per-SV-type wins panel will be "
                    f"thin without it.\n"
                )
    else:
        # In mycosv-only mode, hard-disable every --run-X flag in case the
        # caller mixed flags. This makes the no-comparator path explicit.
        for flag in (
            "run_syri", "run_minigraph", "run_pggb", "run_cactus",
            "run_svim_asm", "run_anchorwave",
            "run_svim", "run_sniffles", "run_cutesv",
            "run_delly", "run_manta",
        ):
            if hasattr(args, flag):
                setattr(args, flag, False)
        sys.stderr.write(
            "[comparators] --mycosv-only: skipping every algorithmic comparator. "
            "exact_benchmark_summary.tsv will use no_comparator placeholder rows; "
            "biology_findings.tsv / novel_mycosv_calls.tsv / TE classification "
            "still flow through.\n"
        )
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

    query_groups = getattr(args, "benchmark_query_genera", "") or ""
    if query_groups:
        query_manifest = select_one_query_per_group(
            query_manifest,
            query_groups,
            out_dir,
            prepared_dir / "hierarchy_manifest.tsv",
        )

    max_benchmark_queries = int(getattr(args, "max_benchmark_queries", 0) or 0)
    if max_benchmark_queries > 0 and len(query_manifest) > max_benchmark_queries:
        original_n = len(query_manifest)
        query_manifest = query_manifest[:max_benchmark_queries]
        sys.stderr.write(
            f"[benchmark] limiting query manifest to first {len(query_manifest)} "
            f"of {original_n} query row(s) via --max-benchmark-queries\n"
        )

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

    query_manifest = cap_read_query_inputs(
        query_manifest,
        out_dir,
        args.mode,
        getattr(args, "max_comparator_short_reads", 500000),
        getattr(args, "max_comparator_long_reads", 200000),
        mycosv_use_full_reads=getattr(args, "mycosv_use_full_reads", False),
    )
    if args.mode == "assembly":
        query_manifest = filter_assembly_query_inputs(
            query_manifest,
            out_dir,
            getattr(args, "max_assembly_query_contigs", 0),
            getattr(args, "max_assembly_query_bp", 0),
        )
        if not query_manifest:
            status_path = out_dir / "NO_QUERIES_AFTER_ASSEMBLY_FILTER.txt"
            status_path.write_text(
                f"All assembly queries were filtered before benchmarking.\n"
                f"max_assembly_query_contigs={getattr(args, 'max_assembly_query_contigs', 0)}\n"
                f"max_assembly_query_bp={getattr(args, 'max_assembly_query_bp', 0)}\n"
                f"See SKIPPED_ASSEMBLY_QUERIES.tsv for the skipped query list.\n",
                encoding="utf-8",
            )
            print(
                f"benchmark_skipped\tmode={args.mode}\treason=assembly_query_filter"
                f"\tstatus_file={status_path}"
            )
            return 0

    if not getattr(args, "skip_input_preflight", False):
        preflight_benchmark_inputs(query_manifest, out_dir, args.mode)

    # Pre-flight: report which comparators are available / missing and write a
    # COMPARATORS_STATUS.txt so the user does not need to grep through logs.
    _report_comparator_status(args, out_dir)

    # Write a mode-filtered query_list.txt that the binary consumes, so reads
    # modes get FASTQ paths and assembly mode gets FASTA paths.
    mode_query_list = out_dir / "query_list.filtered.txt"
    mode_query_list.write_text(
        "\n".join(row.get("mycosv_path") or row["path"] for row in query_manifest) + "\n",
        encoding="utf-8",
    )

    # Build a benchmark-scoped ref_list that contains the references each
    # mode-filtered query is supposed to be compared against, *plus* a bounded
    # neighborhood of phylogenetically related refs from prepared_dir.
    #
    # History: limiting this list to one benchmark_ref_fasta per query (5 refs
    # for the million_real held-out queries) capped MycoSV's per-query call
    # volume an order of magnitude below the comparators. The C++ binary
    # honors --max-ref-memory-mb so an over-broad list is safe — refs are
    # loaded sequentially until the cap (now 8 GB; ≈200 fungal refs) — but a
    # list that misses every related ref forces the per-query MEM-chain
    # search to fall back on cross-genus refs (~5 % overlap) instead of the
    # genus-mate the user actually deposited. Two queries (Lodderomyces,
    # Dactylellina) ended up with 0 calls against their own benchmark refs.
    #
    # Strategy: start with the per-query benchmark refs, then give each query
    # a fair share of NEIGHBOR_REF_CAP — adding genus → family → order → class
    # neighbors per query before moving on. A per-query budget is essential
    # because one query's genus can dominate the corpus (~1100 Saccharomyces
    # refs in the million_real manifest); a global walk would let it swallow
    # every slot and leave the other four queries with no neighbors.
    # The C++ side then drops anything below the k=16 prefilter overlap
    # floor (0.02), so distant refs don't poison the chain shortlist.
    NEIGHBOR_REF_CAP = max(1, int(getattr(args, "benchmark_ref_cap", 512)))
    bench_refs: list[str] = []
    seen_bench_refs: set[str] = set()
    # Pre-seed seen_bench_refs with every query's OWN fasta path so the
    # genus/family/order/class neighbor walk below cannot re-add the query
    # itself as a "neighbor." The query is registered in hierarchy_manifest.tsv
    # under its own genus, so without this exclusion the chain-search anchors
    # snap to the query's own contigs (perfect self-identity), suppressing
    # every DUP/TRA call and forcing Path C into the OFF_REF tile sweep —
    # observed on the F. falciforme vs F. oxysporum Fo47 run as 0 DUP, 0 TRA,
    # and 3,584 phantom NOVEL_WEAK windows.
    for row in query_manifest:
        qp = (row.get("path") or "").strip()
        if qp and qp != ".":
            seen_bench_refs.add(qp)
    per_query_taxa: list[dict[str, str]] = []
    for row in query_manifest:
        bench_ref = (row.get("benchmark_ref_fasta") or "").strip()
        if bench_ref and bench_ref != "." and bench_ref not in seen_bench_refs \
                and Path(bench_ref).exists():
            seen_bench_refs.add(bench_ref)
            bench_refs.append(bench_ref)
        per_query_taxa.append({
            rank: (row.get(rank) or "").strip()
            for rank in ("genus", "family", "order", "class")
        })

    hierarchy_path = prepared_dir / "hierarchy_manifest.tsv"
    if hierarchy_path.exists() and per_query_taxa:
        try:
            with hierarchy_path.open(encoding="utf-8") as fh:
                manifest_rows = list(csv.DictReader(fh, delimiter="\t"))
        except (OSError, csv.Error) as exc:
            sys.stderr.write(
                f"[bench-ref] could not read {hierarchy_path}: {exc}; "
                f"falling back to per-query benchmark refs only\n"
            )
            manifest_rows = []
        n_queries = max(1, len(per_query_taxa))
        # Reserve at least 4 slots per query so even the smallest genera get
        # neighbors; the global cap is the hard ceiling.
        per_query_budget = max(4, NEIGHBOR_REF_CAP // n_queries)
        for taxa in per_query_taxa:
            if len(bench_refs) >= NEIGHBOR_REF_CAP:
                break
            added_for_query = 0
            for rank in ("genus", "family", "order", "class"):
                if added_for_query >= per_query_budget:
                    break
                if len(bench_refs) >= NEIGHBOR_REF_CAP:
                    break
                wanted = taxa.get(rank, "")
                if not wanted or wanted == ".":
                    continue
                for hrow in manifest_rows:
                    if added_for_query >= per_query_budget:
                        break
                    if len(bench_refs) >= NEIGHBOR_REF_CAP:
                        break
                    v = (hrow.get(rank) or "").strip()
                    if v != wanted:
                        continue
                    fasta = (hrow.get("fasta_path") or "").strip()
                    if fasta and fasta not in seen_bench_refs and Path(fasta).exists():
                        seen_bench_refs.add(fasta)
                        bench_refs.append(fasta)
                        added_for_query += 1

    bench_ref_list_path: Path | None = None
    if bench_refs:
        bench_ref_list_path = out_dir / "ref_list.benchmark.txt"
        bench_ref_list_path.write_text("\n".join(bench_refs) + "\n", encoding="utf-8")
        sys.stderr.write(
            f"[bench-ref] wrote {len(bench_refs)} refs to {bench_ref_list_path} "
            f"(per-query benchmark refs + genus/family/order/class neighbors, "
            f"capped at {NEIGHBOR_REF_CAP})\n"
        )

    compile_binary_if_needed(args.binary_path.resolve(), force=args.force_rebuild)
    mycosv_failed = False
    out_prefix = out_dir / "mycosv" / "calls"
    prior_mycosv_outputs = snapshot_mycosv_outputs(out_prefix)
    try:
        mycosv_paths = run_mycosv(
            prepared_dir, out_dir, args.binary_path.resolve(), args.mode, args.mycosv_arg,
            query_list_override=mode_query_list,
            threads=args.threads,
            max_clade_genomes=args.max_clade_genomes,
            reuse_index_dir=getattr(args, "reuse_index_dir", None),
            reuse_registry_dir=getattr(args, "reuse_registry_dir", None),
            ref_list_override=bench_ref_list_path,
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        # Don't abort the whole panel/mode pipeline when the binary crashes
        # (e.g. cgroup OOM kill mid-write produces a 0-byte calls.vcf, then
        # the bare exception bubbles to main and skips downstream report
        # writes). Mark the failure on disk, then prefer the per-contig
        # hierarchical checkpoint if it contains calls; that checkpoint is the
        # completed pangenome caller output for all flushed contigs and is more
        # informative than a stale canonical VCF or an empty failure callset.
        mycosv_failed = True
        rc = getattr(exc, "returncode", "timeout")
        checkpoint_paths = promote_hierarchical_checkpoint(out_prefix)
        if checkpoint_paths is not None:
            mycosv_paths = checkpoint_paths
            checkpoint_n = vcf_data_record_count(Path(checkpoint_paths["vcf"]))
            failure_action = (
                f"promoted the hierarchical per-contig checkpoint "
                f"({checkpoint_n} calls) to the canonical MycoSV outputs"
            )
        elif prior_mycosv_outputs:
            restore_mycosv_outputs(prior_mycosv_outputs)
            mycosv_paths = {
                "vcf": str(out_prefix.with_suffix(".vcf")),
                "hits": str(out_prefix.with_suffix(".hits.tsv")),
                "gfa": str(out_prefix.with_suffix(".gfa")),
            }
            failure_action = (
                "restored the previous MycoSV outputs instead of replacing "
                "them with an empty failure callset"
            )
        else:
            mycosv_paths = write_mycosv_failure_outputs(out_prefix, f"rc={rc}")
            failure_action = "wrote an empty failure callset"
        marker = out_dir / "MYCOSV_FAILED.txt"
        marker.write_text(
            f"MycoSV binary failed (rc={rc}). The benchmark {failure_action} "
            f"so panel-level reports still render.\n"
            f"Common causes: cgroup OOM kill (see [mycosv] line above), "
            f"missing index files, or unreadable FASTA.\n",
            encoding="utf-8",
        )
        sys.stderr.write(
            f"[benchmark] mycosv failed for mode={args.mode!r}; continuing "
            f"after {failure_action}. Marker: {marker}\n"
        )

    # Materialize a multi-sample sibling of calls.vcf so spot-checks see one
    # column per query asm (the binary writes a single SAMPLE column with
    # provenance only in the QASM info field).
    try:
        vcf_path = Path(mycosv_paths["vcf"])
        multisample_vcf = expand_to_multisample_vcf(
            vcf_path,
            vcf_path.with_suffix(".multisample.vcf"),
            [row["query_asm"] for row in query_manifest],
        )
    except Exception as exc:
        sys.stderr.write(f"[multisample-vcf] expand failed: {exc}\n")

    mycosv_calls_by_query: dict[str, dict[str, Any]] = {}
    for row in query_manifest:
        query_asm = row["query_asm"]
        q_calls, r_calls, r_to_q_keys = load_mycosv_paired_calls(
            Path(mycosv_paths["vcf"]), query_asm,
        )
        # `ref_to_query_keys` is a parallel list (same length, same order as
        # r_calls). The phase-1 lift fallback preserves order and length —
        # every input ref-call yields exactly one output, with at most pos/end
        # rewritten — so this index mapping survives the lift step at 6755
        # below. After bench_contigs filtering we propagate the original index
        # alongside each kept ref-call so the support-tracking step at 6905
        # can still reach back to the QUERY-coord call_key for the same VCF
        # row, fixing the always-`mycosv_unique=yes` bug in biology_findings.
        mycosv_calls_by_query[query_asm] = {
            "query": q_calls,
            "reference": r_calls,
            "ref_to_query_keys": r_to_q_keys,
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

    # SyRI / minigraph / pggb produce comparator callsets when they succeed; any tool
    # failure on a single query (SyRI's CalledProcessError on weak alignments,
    # minigraph crashes on huge gaps, pggb timeouts) is captured here so it
    # does not abort the entire benchmark — surviving comparators still
    # contribute to exact_benchmark_summary.tsv.
    # Comparators are scheduled across queries via run_per_query_in_parallel so
    # the 5-query x N-tool bench fan-out lands inside the per-mode wall budget;
    # the old serial loop spent ~N x single_query_time and routinely blew the
    # 4 h cap with one slow svim_asm / minigraph / sniffles per mode.
    if args.run_syri and args.mode == "assembly":
        results = run_per_query_in_parallel(
            "syri", run_syri_for_query, query_manifest, out_dir, args.threads,
        )
        for query_row in query_manifest:
            result = results.get(query_row["query_asm"])
            if result:
                truth_sets[query_row["query_asm"]][("query", "syri")] = load_syri_query_calls(
                    Path(result["normalized_tsv"]), query_row["query_asm"]
                )

    if args.run_minigraph and args.mode == "assembly":
        minigraph_args = args.minigraph_arg
        def _minigraph_runner(qr, od, t):
            return run_minigraph_for_query(qr, od, t, minigraph_args)
        results = run_per_query_in_parallel(
            "minigraph", _minigraph_runner, query_manifest, out_dir, args.threads,
        )
        for query_row in query_manifest:
            result = results.get(query_row["query_asm"])
            if result:
                truth_sets[query_row["query_asm"]][("reference", "minigraph")] = load_minigraph_bubble_calls(
                    Path(result["bubble_bed"]),
                    Path(result["sample_bed"]),
                    query_row["query_asm"],
                )

    if args.run_pggb and args.mode == "assembly":
        pggb_identity = args.pggb_identity
        pggb_segment_len = args.pggb_segment_len
        pggb_args = args.pggb_arg
        def _pggb_runner(qr, od, t):
            return run_pggb_for_query(qr, od, t, pggb_identity, pggb_segment_len, pggb_args)
        results = run_per_query_in_parallel(
            "pggb", _pggb_runner, query_manifest, out_dir, args.threads,
        )
        for query_row in query_manifest:
            result = results.get(query_row["query_asm"])
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
            cactus_args = args.cactus_arg
            def _cactus_runner(qr, od, t):
                return run_cactus_for_query(qr, od, t, cactus_args)
            assembly_caller_specs.append(("cactus", _cactus_runner))
        if args.run_svim_asm:
            assembly_caller_specs.append(("svim_asm", run_svim_asm_for_query))
        if args.run_anchorwave:
            assembly_caller_specs.append(("anchorwave", run_anchorwave_for_query))

    for label, runner in assembly_caller_specs:
        results = run_per_query_in_parallel(
            label, runner, query_manifest, out_dir, args.threads,
        )
        for query_row in query_manifest:
            result = results.get(query_row["query_asm"])
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
        results = run_per_query_in_parallel(
            label, runner, query_manifest, out_dir, args.threads,
        )
        for query_row in query_manifest:
            result = results.get(query_row["query_asm"])
            if not result:
                continue
            vcf_path = Path(result["vcf"])
            truth_sets[query_row["query_asm"]][("reference", label)] = load_reference_vcf_calls(
                vcf_path, label, query_row["query_asm"]
            )

    agreement_rows: list[dict[str, Any]] = []
    read_validated_truth_rows: list[dict[str, Any]] = []
    # Leave-one-out comparator-variance benchmark (fungal-specific).
    # loo_summary_rows is the per-(fold, stratum) flattened TSV;
    # loo_variance_by_query is the per-query aggregate (mean / stdev / swing /
    # most-influential comparator) consumed by loo_consensus_variance.json.
    loo_summary_rows: list[dict[str, Any]] = []
    loo_variance_by_query: dict[str, dict[str, Any]] = {}
    # Per-call diagnostic: for every predicted mycosv call that failed to
    # TP against a truth call, record why (contig_mismatch / pos_out_of_tol
    # / svlen_out_of_tol / type_mismatch / mate_mismatch / claimed_by_other_pred)
    # and the closest truth candidate. Written to match_failures.tsv next to
    # exact_benchmark_summary.tsv so an operator can see why recall is 0
    # without re-running.
    match_failure_rows: list[dict[str, Any]] = []
    # Write a header-only placeholder up front so a SLURM time-out (or any
    # mid-run kill) still leaves the visualization a parseable file. The final
    # write_tsv at the end of benchmark() overwrites this with the populated
    # rows when we reach it.
    write_tsv(out_dir / "read_validated_truth.tsv", [], READ_VALIDATION_FIELDS)
    support_by_key: dict[tuple[str, str, int, int, str], list[str]] = defaultdict(list)
    # call_keys of MycoSV preds that passed external read validation, keyed by
    # query_asm. Populated inside the per-query loop and consumed after it to
    # tier every MycoSV call into strong/moderate/intrinsic_only/weak.
    validated_mycosv_keys_by_query: dict[str, set[tuple[str, str, int, int, str]]] = defaultdict(set)
    summary_json: dict[str, Any] = {
        "mode": args.mode,
        "prepared_dir": str(prepared_dir),
        "mycosv_paths": mycosv_paths,
        "mycosv_status": "failed" if mycosv_failed else "ok",
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
    sv_volume_audit_rows = write_sv_volume_audit(
        out_dir / "sv_volume_audit.tsv",
        query_manifest,
        mycosv_calls_by_query,
        truth_sets,
        Path(mycosv_paths["vcf"]),
        mode=args.mode,
        mycosv_failed=mycosv_failed,
    )
    summary_json["sv_volume_audit"] = sv_volume_audit_rows

    data_cache_dir = (
        getattr(args, "data_cache_dir", None) or DEFAULT_DATA_CACHE
    ).resolve()
    lift_cache_dir = out_dir / "lift_cache"
    single_ref_equivalent_counts: dict[str, int] = {}
    single_ref_equivalent_keys_by_query: dict[
        str, set[tuple[str, str, int, int, str]]
    ] = defaultdict(set)
    single_ref_equivalent_loci_by_query: dict[
        str, set[tuple[str, str, str, int, int, int, str, int]]
    ] = defaultdict(set)
    for query_row in query_manifest:
        query_asm = query_row["query_asm"]
        mycosv_query_calls = mycosv_calls_by_query.get(query_asm, {}).get("query", [])
        mycosv_ref_calls_all = mycosv_calls_by_query.get(query_asm, {}).get("reference", [])
        # Parallel to mycosv_ref_calls_all (same length, same order) — gives us
        # the QUERY-coord call_key for each ref-coord call, even after the lift
        # rewrites pos/end. Used below to propagate per-comparator support back
        # to the right entry in support_by_key.
        mycosv_ref_to_query_keys = list(
            mycosv_calls_by_query.get(query_asm, {}).get("ref_to_query_keys", [])
        )
        # MycoSV is pangenomic in every mode: it may anchor a query locus on a
        # close sibling ref rather than the one single-reference comparators
        # used. Project all reference-coordinate MycoSV calls to the benchmark
        # ref first, then filter to benchmark contigs. Without this, assembly
        # mode silently discarded sibling-ref calls instead of comparing the
        # single-reference-equivalent projection.
        bench_ref = query_row.get("benchmark_ref_fasta") or "."
        bench_contigs: frozenset[str] = frozenset()
        if bench_ref not in {"", "."}:
            bench_contigs = fasta_contig_names(Path(bench_ref))
        if bench_contigs and tool_path("minimap2") is not None:
            mycosv_ref_calls_all = _lift_calls_to_benchmark_ref(
                mycosv_ref_calls_all,
                Path(bench_ref),
                data_cache_dir,
                lift_cache_dir / query_asm,
                threads=args.threads,
            )
        if bench_contigs:
            paired_keep = [
                (c, mycosv_ref_to_query_keys[i] if i < len(mycosv_ref_to_query_keys) else None)
                for i, c in enumerate(mycosv_ref_calls_all)
                if c.ref_contig in bench_contigs
            ]
            mycosv_ref_projected_raw = len(paired_keep)
            (
                mycosv_ref_calls,
                mycosv_ref_query_keys_kept,
                single_ref_all_query_keys,
            ) = deduplicate_projected_reference_calls(paired_keep)
        else:
            paired_keep = [
                (c, mycosv_ref_to_query_keys[i] if i < len(mycosv_ref_to_query_keys) else None)
                for i, c in enumerate(mycosv_ref_calls_all)
            ]
            mycosv_ref_projected_raw = len(paired_keep)
            (
                mycosv_ref_calls,
                mycosv_ref_query_keys_kept,
                single_ref_all_query_keys,
            ) = deduplicate_projected_reference_calls(paired_keep)
        single_ref_equivalent_counts[query_asm] = len(mycosv_ref_calls)
        single_ref_equivalent_keys_by_query[query_asm] = single_ref_all_query_keys
        single_ref_equivalent_loci_by_query[query_asm] = {
            pangenome_locus_key(c)
            for c in mycosv_query_calls
            if call_key(c) in single_ref_equivalent_keys_by_query[query_asm]
        }
        # Count OFF_REF events that exist in the query-coord set but are
        # un-matchable (no REFPOS, so they were dropped from the ref-coord set).
        # These are real predictions that the single-reference truth callers
        # can't represent, so reporting them separately stops them from being
        # silently invisible in exact_benchmark_summary.tsv.
        mycosv_off_ref_dropped = sum(
            1 for c in mycosv_query_calls if c.svtype == "OFF_REF"
        )
        # Misrouted: ref-coord calls whose REFCONTIG is from a sibling clade,
        # not the benchmark target. They're discarded by the bench_contigs
        # filter but counted here so the operator sees the routing penalty.
        mycosv_misrouted = max(0, len(mycosv_ref_calls_all) - mycosv_ref_projected_raw)
        summary_json["queries"][query_asm] = {
            "mycosv_calls": {
                "query": len(mycosv_query_calls),
                "reference": len(mycosv_ref_calls),
                "reference_projected_raw": mycosv_ref_projected_raw,
                "reference_total": len(mycosv_ref_calls_all),
                "benchmark_ref_contigs": len(bench_contigs),
                "off_ref_dropped": mycosv_off_ref_dropped,
                "misrouted_to_sibling_clade": mycosv_misrouted,
                "benchmark_ref_filter_kept_fraction": (
                    round(mycosv_ref_projected_raw / len(mycosv_ref_calls_all), 4)
                    if mycosv_ref_calls_all else 1.0
                ),
            },
            "exact_benchmarks": {"query": {}, "reference": {}, "reference_any_clade": {}},
        }
        truth_for_query = truth_sets.get(query_asm, {})
        for (coord_space, label), truth_calls in truth_for_query.items():
            pred_calls = mycosv_ref_calls if coord_space == "reference" else mycosv_query_calls
            agreement_rows.extend(_emit_per_svtype_rows(
                query_asm=query_asm,
                coord_space=coord_space,
                truth_label=label,
                method="mycosv",
                truth_calls=truth_calls,
                pred_calls=pred_calls,
            ))
            for mf in diagnose_match_failures(truth_calls, pred_calls):
                match_failure_rows.append({
                    "query_asm": query_asm,
                    "coordinate_space": coord_space,
                    "truth_label": label,
                    "method": "mycosv",
                    **mf,
                })
            metrics = score_callsets(truth_calls, pred_calls)
            if coord_space == "query":
                used_pred, _ = match_calls(truth_calls, pred_calls)
                for idx in used_pred:
                    support_by_key[call_key(pred_calls[idx])].append(label)
            elif coord_space == "reference":
                # Bridge ref-coord matches back to the query-coord call_key
                # via the parallel mycosv_ref_query_keys_kept list. Without
                # this, every mycosv call in biology_findings.tsv was flagged
                # `mycosv_unique=yes` because comparator-vs-mycosv matches
                # only happen in reference space for most real-data tools
                # (svim/sniffles/cutesv/delly/manta/minigraph all emit
                # ref-coord truth) and that side of the matching used to
                # write nothing to support_by_key.
                used_pred, _ = match_calls(truth_calls, pred_calls)
                for idx in used_pred:
                    qkey = (
                        mycosv_ref_query_keys_kept[idx]
                        if idx < len(mycosv_ref_query_keys_kept) else None
                    )
                    if qkey is not None:
                        support_by_key[qkey].append(label)
            summary_json["queries"][query_asm]["exact_benchmarks"][coord_space][label] = metrics

            # Per-comparator read-level validation: anchor THIS comparator's
            # truth in the raw query data and re-score. Emits a `<label>
            # _read_supported` row alongside the algorithmic-truth row, so
            # the visualization can ask "does mycosv beat minigraph after we
            # filter minigraph's calls to those with raw read support?". On
            # by default; off via --no-validate-with-reads.
            if (
                getattr(args, "validate_with_reads", False)
                and coord_space == "reference"
                and tool_path("samtools") is not None
                and tool_path("minimap2") is not None
            ):
                val_dir = out_dir / "read_validation" / query_asm / label
                kept_truth, support_rows = validate_calls_with_reads(
                    truth_calls,
                    query_row,
                    val_dir,
                    threads=args.threads,
                    min_support=args.read_validation_min_support,
                    flank_bp=args.read_validation_flank_bp,
                )
                read_validated_truth_rows.extend(support_rows)
                rs_label = f"{label}_read_supported"
                agreement_rows.extend(_emit_per_svtype_rows(
                    query_asm=query_asm,
                    coord_space=coord_space,
                    truth_label=rs_label,
                    method="mycosv",
                    truth_calls=kept_truth,
                    pred_calls=pred_calls,
                ))
                summary_json["queries"][query_asm]["exact_benchmarks"][coord_space][rs_label] = {
                    **score_callsets(kept_truth, pred_calls),
                    "comparator": label,
                    "comparator_calls": len(truth_calls),
                    "read_validated": len(kept_truth),
                    "min_split_reads": args.read_validation_min_support,
                }

        # Fix B: "any-clade" rows — score every reference-space comparator
        # against the *unfiltered* mycosv ref-coord set so the operator can
        # tell apart "mycosv called the wrong thing" from "mycosv called the
        # right thing on the wrong reference clade". Same truth, different
        # predictions; only emitted in reference space.
        for (coord_space, label), truth_calls in truth_for_query.items():
            if coord_space != "reference":
                continue
            agreement_rows.extend(_emit_per_svtype_rows(
                query_asm=query_asm,
                coord_space="reference_any_clade",
                truth_label=label,
                method="mycosv",
                truth_calls=truth_calls,
                pred_calls=mycosv_ref_calls_all,
            ))
            summary_json["queries"][query_asm]["exact_benchmarks"]["reference_any_clade"][label] = score_callsets(truth_calls, mycosv_ref_calls_all)

        # Comparator consensus per coord_space — a candidate call is in the
        # consensus iff it is supported by >=2 comparators (calls_compatible
        # across position+length+type tolerance). This dilutes single-tool
        # bias (e.g. minigraph's bubble-extraction conservatism, svim_asm's
        # alt-allele preference) by keeping only events that ≥2 independent
        # callers agree on, then scores mycosv against that consensus.
        for coord_space in ("query", "reference"):
            ref_labels = [
                lbl for (cs, lbl) in truth_for_query.keys()
                if cs == coord_space and lbl != "consensus_2of_n"
            ]
            if len(ref_labels) < 2:
                continue
            consensus_calls = build_consensus_truth(
                [truth_for_query[(coord_space, lbl)] for lbl in ref_labels],
                min_support=2,
            )
            pred_calls = mycosv_ref_calls if coord_space == "reference" else mycosv_query_calls
            consensus_label = f"consensus_2of_{len(ref_labels)}"
            agreement_rows.extend(_emit_per_svtype_rows(
                query_asm=query_asm,
                coord_space=coord_space,
                truth_label=consensus_label,
                method="mycosv",
                truth_calls=consensus_calls,
                pred_calls=pred_calls,
            ))
            summary_json["queries"][query_asm]["exact_benchmarks"][coord_space][
                consensus_label
            ] = score_callsets(consensus_calls, pred_calls)

            # ── Leave-one-out comparator-variance benchmark ──────────────
            # Replays the consensus K times, each time excluding one
            # comparator, and reports the F1 dispersion + the most
            # influential comparator. Fungal-specific strata (length bin /
            # element class / phylum) let the reader see whether the
            # variance lives in boring small INDELs or in the biologically
            # interesting HGT / STARSHIP / arm-level events. Needs ≥3
            # comparators so each LOO fold still has ≥2 left to vote.
            phylum_label = str(query_row.get("phylum", ".") or ".")
            if len(ref_labels) >= 3:
                loo_inputs = {
                    lbl: truth_for_query[(coord_space, lbl)] for lbl in ref_labels
                }
                loo = score_loo_consensus(
                    pred_calls,
                    loo_inputs,
                    min_support=2,
                    query_phylum=phylum_label,
                )
                summary_json["queries"][query_asm]["exact_benchmarks"][coord_space][
                    f"{consensus_label}_loo"
                ] = loo
                loo_variance_by_query.setdefault(query_asm, {})[coord_space] = {
                    k: v for k, v in loo.items() if k != "loo_folds"
                }
                loo_summary_rows.extend(_emit_loo_summary_rows(
                    query_asm=query_asm,
                    query_phylum=phylum_label,
                    coord_space=coord_space,
                    loo=loo,
                ))
            else:
                # LOO needs ≥3 comparators (each fold has K-1 ≥ 2 left so
                # consensus_2of_(K-1) is well-defined). Without it we'd
                # silently emit an empty loo_consensus_summary.tsv and the
                # operator would have no idea why. Emit an explicit skipped
                # row so `grep SKIPPED loo_consensus_summary.tsv` reveals
                # the reason at a glance, and record the same reason in the
                # JSON so visualization can show "LOO not run: 2 callers".
                skip_reason = (
                    f"only {len(ref_labels)} comparator(s) available "
                    f"({','.join(sorted(ref_labels))}); need >=3 for LOO. "
                    f"Enable a 3rd comparator (e.g. MILLION_REAL_RUN_CACTUS=1 "
                    f"or --run-anchorwave for assembly; gridss/lumpy for short-reads)."
                )
                summary_json["queries"][query_asm]["exact_benchmarks"][coord_space][
                    f"{consensus_label}_loo"
                ] = {
                    "status":      "skipped",
                    "reason":      skip_reason,
                    "comparators": sorted(ref_labels),
                    "phylum":      phylum_label,
                }
                loo_variance_by_query.setdefault(query_asm, {})[coord_space] = {
                    "status":      "skipped",
                    "reason":      skip_reason,
                    "comparators": sorted(ref_labels),
                    "phylum":      phylum_label,
                }
                loo_summary_rows.append({
                    "query_asm":           query_asm,
                    "phylum":              phylum_label,
                    "coordinate_space":    coord_space,
                    "excluded_comparator": "NONE",
                    "stratum_type":        "SKIPPED",
                    "stratum_value":       skip_reason,
                    "truth_n":             0,
                    "tp":                  0,
                    "fp":                  0,
                    "fn":                  0,
                    "precision":           float("nan"),
                    "recall":              float("nan"),
                    "f1":                  float("nan"),
                    "status":              "skipped_low_comparator_count",
                })

            # Independent read-level validation: re-anchor the consensus
            # candidates in the raw query data (FASTQ for reads-mode, contigs for
            # assembly-mode) by counting split / clipped reads spanning the
            # breakpoint. SVs without raw-data support are dropped from the
            # candidate set before scoring, removing the assembly-only artefacts
            # that all algorithm comparators inherit. Only emitted when
            # `--validate-with-reads` is on (default), and only for the
            # reference coord space because that's where read-level support
            # has unambiguous coordinates.
            if (
                getattr(args, "validate_with_reads", False)
                and coord_space == "reference"
                and tool_path("samtools") is not None
                and tool_path("minimap2") is not None
            ):
                val_dir = out_dir / "read_validation" / query_asm / "consensus"
                kept_calls, support_rows = validate_calls_with_reads(
                    consensus_calls,
                    query_row,
                    val_dir,
                    threads=args.threads,
                    min_support=args.read_validation_min_support,
                    flank_bp=args.read_validation_flank_bp,
                )
                read_validated_truth_rows.extend(support_rows)
                rv_label = f"{consensus_label}_read_supported"
                agreement_rows.extend(_emit_per_svtype_rows(
                    query_asm=query_asm,
                    coord_space=coord_space,
                    truth_label=rv_label,
                    method="mycosv",
                    truth_calls=kept_calls,
                    pred_calls=pred_calls,
                ))
                summary_json["queries"][query_asm]["exact_benchmarks"][coord_space][rv_label] = {
                    **score_callsets(kept_calls, pred_calls),
                    "consensus_input": len(consensus_calls),
                    "read_validated": len(kept_calls),
                    "min_split_reads": args.read_validation_min_support,
                }

                # Apples-to-apples comparator scoring: each algorithmic
                # comparator is now ALSO scored as a predictor against the
                # same read-validated truth, per SV type. This populates the
                # "MycoSV vs comparator per SV type — wins matrix" panel in
                # the visualization (method=<comparator_label> rows alongside
                # the existing method=mycosv row, all under the shared
                # truth_label=consensus_..._read_supported).
                for cmp_label in ref_labels:
                    cmp_calls = truth_for_query.get((coord_space, cmp_label), [])
                    if not cmp_calls:
                        continue
                    agreement_rows.extend(_emit_per_svtype_rows(
                        query_asm=query_asm,
                        coord_space=coord_space,
                        truth_label=rv_label,
                        method=cmp_label,
                        truth_calls=kept_calls,
                        pred_calls=cmp_calls,
                    ))

        # Always validate MycoSV's own predictions against the held-out
        # reads/assembly. Validate reference and query coordinate calls
        # separately: a query may have a few anchored REF calls plus many
        # OFF_REF query-space calls, and the old `ref_calls or query_calls`
        # fallback silently skipped the latter.
        if getattr(args, "validate_with_reads", False):
            read_validation_summary: dict[str, Any] = {
                "min_split_reads": args.read_validation_min_support,
            }
            supported_preds_by_coord: dict[str, list[NormalizedCall]] = {}
            for validation_source, calls_to_validate in (
                ("mycosv_reference", mycosv_ref_calls),
                ("mycosv_query", mycosv_query_calls),
            ):
                if not calls_to_validate:
                    continue
                val_dir = out_dir / "read_validation" / query_asm / validation_source
                kept_calls, support_rows = validate_calls_with_reads(
                    calls_to_validate,
                    query_row,
                    val_dir,
                    threads=args.threads,
                    min_support=args.read_validation_min_support,
                    flank_bp=args.read_validation_flank_bp,
                )
                read_validated_truth_rows.extend(support_rows)
                coord_key = "reference" if validation_source.endswith("reference") else "query"
                supported_preds_by_coord[coord_key] = kept_calls
                for c in kept_calls:
                    validated_mycosv_keys_by_query[query_asm].add(call_key(c))
                read_validation_summary[validation_source] = {
                    "input_calls": len(calls_to_validate),
                    "read_validated": len(kept_calls),
                }

            if read_validation_summary.keys() - {"min_split_reads"}:
                summary_json["queries"][query_asm]["read_validation"] = read_validation_summary

            # ── Tool-agnostic raw-read validation (read_level_union) ───────
            #   An SV is independently validated iff the RAW READS support it — regardless of
            #   which caller (if any) proposed it. We pool every caller's
            #   candidate loci (all comparators + MycoSV) into a per-coord
            #   union, then keep only those clearing the external split/
            #   clipped-read threshold (force_external=True, so MycoSV's own
                #   anchor/cluster counts cannot self-validate). MycoSV and every
                #   comparator are then scored against this same raw-read
                #   validated set per SV type, removing the comparator-universe bias of the
            #   consensus / per-tool rows.
            for coord_space in ("reference", "query"):
                union_by_key: dict[tuple[str, str, int, int, str], NormalizedCall] = {}
                for (cs, _lbl), cset in truth_for_query.items():
                    if cs != coord_space:
                        continue
                    for c in cset:
                        union_by_key.setdefault(call_key(c), c)
                mycosv_self = (mycosv_ref_calls if coord_space == "reference"
                               else mycosv_query_calls)
                for c in mycosv_self:
                    union_by_key.setdefault(call_key(c), c)
                union_calls = list(union_by_key.values())
                if not union_calls:
                    continue
                val_dir = (out_dir / "read_validation" / query_asm
                           / f"read_level_union_{coord_space}")
                read_level_truth, support_rows = validate_calls_with_reads(
                    union_calls,
                    query_row,
                    val_dir,
                    threads=args.threads,
                    min_support=args.read_validation_min_support,
                    flank_bp=args.read_validation_flank_bp,
                    force_external=True,
                )
                read_validated_truth_rows.extend(support_rows)
                if not read_level_truth:
                    continue
                rl_label = "read_level_union"
                agreement_rows.extend(_emit_per_svtype_rows(
                    query_asm=query_asm,
                    coord_space=coord_space,
                    truth_label=rl_label,
                    method="mycosv",
                    truth_calls=read_level_truth,
                    pred_calls=mycosv_self,
                ))
                summary_json["queries"][query_asm]["exact_benchmarks"][coord_space][rl_label] = {
                    **score_callsets(read_level_truth, mycosv_self),
                    "union_input": len(union_calls),
                    "read_validated": len(read_level_truth),
                    "min_split_reads": args.read_validation_min_support,
                }
                # Same tool-agnostic truth, every comparator as predictor —
                # populates the per-SV-type wins matrix on equal footing.
                for (cs, cmp_label), cmp_calls in truth_for_query.items():
                    if cs != coord_space:
                        continue
                    agreement_rows.extend(_emit_per_svtype_rows(
                        query_asm=query_asm,
                        coord_space=coord_space,
                        truth_label=rl_label,
                        method=cmp_label,
                        truth_calls=read_level_truth,
                        pred_calls=cmp_calls,
                    ))

            # Diagnostic scoring row: same truth, but MycoSV predictions
            # filtered to read-supported calls. This tells apart low F1 caused
            # by unsupported predictions from low F1 despite read support.
            for (coord_space, label), truth_calls in truth_for_query.items():
                supported_preds = supported_preds_by_coord.get(coord_space)
                if supported_preds is None:
                    continue
                agreement_rows.extend(_emit_per_svtype_rows(
                    query_asm=query_asm,
                    coord_space=coord_space,
                    truth_label=label,
                    method="mycosv_read_supported",
                    truth_calls=truth_calls,
                    pred_calls=supported_preds,
                ))

        # Flush per-query so a mid-loop SLURM timeout / OOM leaves the TSV
        # populated with whatever was validated so far, instead of just the
        # header placeholder written at line ~6354.
        write_tsv(
            out_dir / "read_validated_truth.tsv",
            read_validated_truth_rows,
            READ_VALIDATION_FIELDS,
        )

    # If no comparator was available (e.g. all of pggb/minigraph/syri/svim_asm
    # missing on this host) the per-query exact benchmarks above produce zero
    # rows, leaving exact_benchmark_summary.tsv header-only and the merged
    # real_merged.tsv empty. Emit a MycoSV-only placeholder per query so the
    # visualization report has something to render and the operator gets a
    # clear "no comparator was run" signal rather than silent emptiness.
    if not agreement_rows:
        any_tool_present = any(summary_json["tool_status"].values())
        sys.stderr.write(
            f"[benchmark] no comparator produced a baseline callset "
            f"(tools_available={any_tool_present}). Emitting MycoSV-only "
            f"placeholder rows so downstream reports stay populated.\n"
        )
        for query_row in query_manifest:
            query_asm = query_row["query_asm"]
            for coord_space, calls_key in (("query", "query"), ("reference", "reference")):
                preds = mycosv_calls_by_query.get(query_asm, {}).get(calls_key, [])
                # Without comparator baselines or raw-read validation rows,
                # tp / fp / fn are undefined, so use NaN rather than inventing
                # zero false positives for unscored predictions.
                agreement_rows.append({
                    "query_asm": query_asm,
                    "coordinate_space": coord_space,
                    "truth_label": "no_comparator",
                    "validation_basis": validation_basis_for_label("no_comparator"),
                    "svtype": "ALL",
                    "method": "mycosv",
                    "truth_calls": float("nan"),
                    "pred_calls": len(preds),
                    "tp": float("nan"),
                    "fp": float("nan"),
                    "fn": float("nan"),
                    "precision": float("nan"),
                    "recall": float("nan"),
                    "f1": float("nan"),
                    "prec_lo95": float("nan"),
                    "prec_hi95": float("nan"),
                    "rec_lo95": float("nan"),
                    "rec_hi95": float("nan"),
                    "status": "no_truth",
                })

    write_tsv(
        out_dir / "read_validated_truth.tsv",
        read_validated_truth_rows,
        READ_VALIDATION_FIELDS,
    )

    write_tsv(
        out_dir / "match_failures.tsv",
        match_failure_rows,
        [
            "query_asm", "coordinate_space", "truth_label", "method",
            "pred_contig", "pred_pos", "pred_end", "pred_svtype", "pred_svlen",
            "reason", "closest_truth_idx", "closest_pos_delta", "closest_svlen_delta",
        ],
    )

    phyla = sorted({row.get("phylum", ".") for row in query_manifest if row.get("phylum") not in {"", "."}})
    phylum_label = phyla[0] if len(phyla) == 1 else "mixed_fungi"

    # Auto-pickup of gene annotations and expression: if the operator did not
    # pass --gene-annotations-tsv / --expression-tsv explicitly, look for the
    # files prepare wrote next to the manifests. Lets the biology analyzer
    # populate expression_gene / expression_distance_bp without manual flags.
    gene_annotations_tsv = args.gene_annotations_tsv
    if gene_annotations_tsv is None:
        candidate = prepared_dir / "gene_annotations.tsv"
        if candidate.exists():
            gene_annotations_tsv = candidate
            sys.stderr.write(f"[gene-annot] using {candidate}\n")
    expression_tsv = args.expression_tsv
    if expression_tsv is None:
        candidate = prepared_dir / "expression.tsv"
        if candidate.exists():
            expression_tsv = candidate
            sys.stderr.write(f"[expression] using {candidate}\n")
    ecological_traits_tsv = prepared_dir / "ecological_traits.tsv"
    if not ecological_traits_tsv.exists():
        ecological_traits_tsv = None
    fungaltraits_csv = data_cache_dir / "fungaltraits.csv"
    if not fungaltraits_csv.exists():
        fungaltraits_csv = None

    # Evidence-tier panorama. Tier every MycoSV call (query + reference coord)
    # into strong / moderate / intrinsic_only / weak so the visualization can
    # show real-but-unvalidatable calls explicitly instead of letting them
    # vanish into FP territory in the F1 plots. The tier map is keyed by the
    # query-coord call_key so it joins cleanly with novel_mycosv_calls.tsv /
    # biology_findings.tsv (both indexed in query coordinates).
    all_mycosv_calls = [call for rows in mycosv_calls_by_query.values() for call in rows.get("query", [])]
    tier_by_key: dict[tuple[str, str, int, int, str], str] = {}
    for qa, sets in mycosv_calls_by_query.items():
        validated = validated_mycosv_keys_by_query.get(qa, set())
        # Tier each query-coord call. Reference-coord rows from the same SV
        # are siblings of the query-coord row; tiering once on query keys
        # avoids double counting in the panorama panel.
        for call in sets.get("query", []):
            supporters = support_by_key.get(call_key(call), [])
            tier_by_key[call_key(call)] = classify_evidence_tier(
                call,
                has_comparator=bool(supporters),
                has_external_read_support=call_key(call) in validated,
            )

    tier_count_rows: list[dict[str, Any]] = []
    counts: dict[tuple[str, str, str], int] = defaultdict(int)
    for call in all_mycosv_calls:
        tier = tier_by_key.get(call_key(call), "weak")
        counts[(call.query_asm, call.svtype, tier)] += 1
    for (qa, svt, tier), n in sorted(counts.items()):
        tier_count_rows.append({
            "query_asm": qa, "svtype": svt, "tier": tier, "n_calls": n,
        })
    write_tsv(
        out_dir / "mycosv_evidence_tiers.tsv",
        tier_count_rows,
        ["query_asm", "svtype", "tier", "n_calls"],
    )

    pangenome_layer_rows = write_pangenome_call_layers(
        out_dir / "pangenome_call_layers.tsv",
        query_manifest,
        mycosv_calls_by_query,
        single_ref_equivalent_counts,
        single_ref_equivalent_keys_by_query,
        support_by_key,
        validated_mycosv_keys_by_query,
    )
    summary_json["pangenome_call_layers"] = pangenome_layer_rows

    novel_rows = []
    for call in all_mycosv_calls:
        key = call_key(call)
        supporters = sorted(set(support_by_key.get(key, [])))
        read_supported = key in validated_mycosv_keys_by_query.get(call.query_asm, set())
        intrinsic_supported = (call.read_support or 0) >= 2
        in_single_ref = (
            key in single_ref_equivalent_keys_by_query.get(call.query_asm, set())
            or pangenome_locus_key(call) in single_ref_equivalent_loci_by_query.get(call.query_asm, set())
        )
        mycosv_unique = not supporters and not in_single_ref
        if in_single_ref:
            discovery_bucket = "single_reference_equivalent"
        elif supporters:
            discovery_bucket = "comparator_supported"
        elif read_supported:
            discovery_bucket = "pangenome_only_read_supported"
        elif intrinsic_supported:
            discovery_bucket = "pangenome_only_intrinsic_supported"
        else:
            discovery_bucket = "pangenome_only_weak"
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
            "single_reference_equivalent": "yes" if in_single_ref else "no",
            "mycosv_unique": "yes" if mycosv_unique else "no",
            "read_supported": "yes" if read_supported else "no",
            "intrinsic_supported": "yes" if intrinsic_supported else "no",
            "evidence_tier": tier_by_key.get(key, "weak"),
            "discovery_bucket": discovery_bucket,
        })
    write_tsv(
        out_dir / "novel_mycosv_calls.tsv",
        novel_rows,
        ["query_asm", "query_contig", "pos", "end", "svtype", "svlen",
         "annotation", "element_class", "support_count", "support_labels",
         "single_reference_equivalent",
         "mycosv_unique", "read_supported", "intrinsic_supported",
         "evidence_tier", "discovery_bucket"],
    )

    enrichment_rows = write_mycosv_novel_biology_enrichment(
        out_dir / "mycosv_novel_biology_enrichment.tsv",
        all_mycosv_calls,
        support_by_key,
        single_ref_equivalent_keys_by_query,
        single_ref_equivalent_loci_by_query,
        validated_mycosv_keys_by_query,
        gene_annotations_tsv,
    )
    summary_json["mycosv_novel_biology_enrichment"] = enrichment_rows

    candidates_tsv, _ = maybe_run_candidate_analysis(
        out_dir,
        mycosv_paths,
        prepared_dir,
        args.mode,
        phylum_label,
        expression_tsv,
        gene_annotations_tsv,
        ecological_traits_tsv,
        fungaltraits_csv,
        args.ancestral_tsv,
    )
    join_biology_findings(
        candidates_tsv, all_mycosv_calls, support_by_key,
        out_dir / "biology_findings.tsv",
        single_ref_keys_by_query=single_ref_equivalent_keys_by_query,
        single_ref_loci_by_query=single_ref_equivalent_loci_by_query,
        tier_by_key=tier_by_key,
    )

    write_agreement_summary(out_dir / "exact_benchmark_summary.tsv", agreement_rows)
    with (out_dir / "benchmark_summary.json").open("w", encoding="utf-8") as fh:
        json.dump(summary_json, fh, indent=2, sort_keys=True)

    # Leave-one-out comparator-variance outputs (fungal-specific). Always
    # written so consumers can detect "no LOO ran" (header-only file) vs
    # "ran but every query had < 3 comparators" (empty body).
    loo_tsv_fields = [
        "query_asm", "phylum", "coordinate_space", "excluded_comparator",
        "stratum_type", "stratum_value", "truth_n", "tp", "fp", "fn",
        "precision", "recall", "f1", "status",
    ]
    write_tsv(out_dir / "loo_consensus_summary.tsv", loo_summary_rows, loo_tsv_fields)
    loo_variance_doc: dict[str, Any] = {
        "queries":  loo_variance_by_query,
        "global":   _summarize_loo_variance_global(loo_variance_by_query),
    }
    with (out_dir / "loo_consensus_variance.json").open("w", encoding="utf-8") as fh:
        json.dump(_sanitize_for_json(loo_variance_doc), fh, indent=2, sort_keys=True, allow_nan=False)

    print(f"benchmark_complete\tqueries={len(query_manifest)}\texact_rows={len(agreement_rows)}"
          f"\tloo_rows={len(loo_summary_rows)}\tloo_queries={len(loo_variance_by_query)}")
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
    cache_dir = data_cache_base(args, out_dir)
    refs_dir = cache_dir / "refs"
    refs_dir.mkdir(parents=True, exist_ok=True)
    index_dir = out_dir / "index"
    registry_dir = out_dir / "registry"
    index_dir.mkdir(parents=True, exist_ok=True)
    registry_dir.mkdir(parents=True, exist_ok=True)
    reset_annotation_source_tally()

    # Step 1: pull the NCBI assembly summary and select up to --max-assemblies
    # fungal rows. We reuse select_all_public_rows so the quality/sorting
    # behavior matches the `prepare --all-public-assemblies` path.
    # flush=True on every progress print: when stdout is piped through `tee`
    # (run_all_experiments.sh does this), it switches from line- to block-
    # buffered, so per-100 download progress was hidden for hours and made
    # the step look frozen even when it was making real progress.
    all_rows, assembly_summary_caches = fetch_ncbi_assembly_rows(
        args.source,
        cache_dir,
        progress=True,
    )

    selected = select_all_public_rows(
        all_rows,
        min_assembly_level=args.min_assembly_level,
        latest_only=args.latest_only,
        max_total=args.max_assemblies,
    )
    if not selected:
        raise ValueError(f"No fungal assemblies matched the selection filters in {args.source}")
    print(f"      selected {len(selected)} assemblies for real-data indexing", flush=True)

    # Step 2: resolve taxonomy lineages for all selected rows.
    print("[2/4] Resolving NCBI taxonomy lineages...", flush=True)
    taxids = sorted({row.get("taxid", "") for row in selected if row.get("taxid")})
    taxonomy_cache = fetch_taxonomy_lineages(taxids, cache_path=cache_dir / "taxonomy_cache.json")

    # Step 3: download each assembly FASTA (or re-use existing on disk) and
    # build the hierarchy manifest the MycoSV binary consumes.
    # Tighter progress reporting + per-row time budget: with ~10000 assemblies
    # and a shared cache, most rows are no-op cache hits but a single hung
    # download (NCBI 5xx, slow mirror) used to silently stall the whole loop.
    # Emit one progress line every 200 examined rows + one every 60 s of wall
    # clock so the operator can tell the difference between "downloading" and
    # "stuck", and bail individual rows that exceed _HTTP_TIMEOUT.
    # Parallel download. The serial loop spent ~3 h on ~3300/10000 rows in
    # production because every NCBI fetch was synchronous; threading is safe
    # since materialize_entry is per-URL and disk-cached. Workers are tunable
    # via MILLION_REAL_DOWNLOAD_WORKERS. Lowered from 8 → 6 after the
    # 2026-05-15 prep run (slurm-14936460) still produced hundreds of 503
    # retries with 8 workers; defense-in-depth comes from the
    # _NCBI_HOST_SEM/_NCBI_COOLDOWN pair above, which throttles even when
    # this env var pushes the worker count back up.
    download_workers = max(
        1, int(os.environ.get("MILLION_REAL_DOWNLOAD_WORKERS",
                              str(_HTTP_MAX_PARALLEL_FTP_NCBI)))
    )
    print(
        f"[3/4] Downloading up to {len(selected)} assemblies -> {refs_dir} "
        f"(workers={download_workers})",
        flush=True,
    )
    ref_manifest_rows: list[dict[str, str]] = []
    ref_list_paths: list[str] = []
    source_link_rows: list[dict[str, str]] = []
    # Per-asm gene-annotation source path (GFF preferred, GBFF fallback) so
    # the prepared dir can grow a gene_annotations.tsv the benchmark step
    # auto-picks up via load_gene_annotations(). Only populated when
    # --million-real-download-gff is on.
    gff_pairs: list[tuple[str, Path]] = []
    download_gff = bool(getattr(args, "million_real_download_gff", False))

    def _download_one(
        row: dict[str, str],
    ) -> tuple[dict[str, str] | None, dict[str, str] | None, tuple[str, Path] | None]:
        asm_name = row.get("assembly_accession", "").replace(".", "_")
        if not asm_name:
            return None, None, None
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
                sys.stderr.flush()
                fasta_path = None
            break
        if fasta_path is None or not fasta_path.exists():
            return None, None, None
        species = lineage.get("species") or row.get("organism_name", ".") or "."
        manifest_row = {
            "asm_name": asm_name,
            "phylum": lineage.get("phylum", "."),
            "class": lineage.get("class", "."),
            "order": lineage.get("order", "."),
            "family": lineage.get("family", "."),
            "genus": lineage.get("genus", species.split()[0] if species not in {".", ""} else "."),
            "clade_name": species,
            "clade_rank": "species",
            "fasta_path": str(fasta_path),
        }
        source_row = {
            "query_asm": asm_name,
            "role": "ref",
            "query_mode": "assembly",
            "source_type": "ncbi_assembly",
            "source_accession": row.get("assembly_accession", "."),
            "source_url": row.get("ftp_path", "."),
            "local_path": str(fasta_path),
            "species": species,
        }
        gff_pair: tuple[str, Path] | None = None
        if download_gff:
            try:
                annotation_source = download_ncbi_gene_annotation_source(
                    row, asm_name, refs_dir, "",
                    ensembl_cache_dir=cache_dir,
                )
            except Exception as exc:
                sys.stderr.write(
                    f"[warn] gene-annot fetch failed for {asm_name}: {exc}\n"
                )
                sys.stderr.flush()
                annotation_source = None
            if annotation_source is not None:
                gff_pair = (asm_name, annotation_source)
        return manifest_row, source_row, gff_pair

    examined = 0
    download_count = 0
    gff_count = 0
    last_progress_t = time.monotonic()
    # Pre-warm the Ensembl Fungi caches BEFORE the ThreadPoolExecutor starts:
    # the release-discovery + species-TSV downloads write shared files under
    # cache_dir, and N workers racing on first-fetch would all call
    # http_download() against the same .part file. Warming serially upfront
    # populates the on-disk cache, after which every worker hits the fast
    # path. Errors here are non-fatal — _ensembl_fungi_gff_url tolerates a
    # missing index and falls through to NCBI GBFF.
    if download_gff:
        try:
            _ensembl_fungi_release(cache_dir)
            _ensembl_fungi_species_index(cache_dir)
        except Exception as exc:  # noqa: BLE001 - best-effort warm-up
            sys.stderr.write(f"[ensembl-fungi] cache warm-up failed: {exc}\n")
    from concurrent.futures import ThreadPoolExecutor, as_completed

    with ThreadPoolExecutor(max_workers=download_workers) as pool:
        futures = {pool.submit(_download_one, row): row for row in selected}
        for fut in as_completed(futures):
            examined += 1
            try:
                manifest_row, source_row, gff_pair = fut.result()
            except Exception as exc:
                sys.stderr.write(f"[warn] download worker raised: {exc}\n")
                sys.stderr.flush()
                manifest_row = source_row = None
                gff_pair = None
            if manifest_row is not None and source_row is not None:
                ref_manifest_rows.append(manifest_row)
                ref_list_paths.append(manifest_row["fasta_path"])
                source_link_rows.append(source_row)
                download_count += 1
                if gff_pair is not None:
                    gff_pairs.append(gff_pair)
                    gff_count += 1
            now = time.monotonic()
            if examined % 200 == 0 or (now - last_progress_t) >= 60.0:
                print(
                    f"      ... examined {examined}/{len(selected)} "
                    f"available={download_count} gff={gff_count}",
                    flush=True,
                )
                last_progress_t = now

    if not ref_manifest_rows:
        raise RuntimeError("No assemblies were successfully downloaded — aborting indexing.")
    print(f"      downloaded/cached {download_count} assemblies (examined {examined})", flush=True)
    annot_summary = annotation_source_summary()
    if annot_summary:
        print(f"      {annot_summary}", flush=True)

    # ── Hold out a small subset as MycoSV-only benchmark queries ─────────────
    # The downstream `benchmark` sub-command is fed by query_manifest.tsv +
    # ref_list.txt; without these, the million-real artifact is just an index
    # with no end-to-end SV-call/biology/visualization payload. Reserve K
    # assemblies as queries (sampled stride-uniformly across phyla so we get
    # taxonomic diversity even when --max-assemblies is small), and exclude
    # them from the index manifest. Per-query benchmark_ref_fasta is the
    # closest sibling in the same genus → family → phylum, falling back to
    # the first ref so read-level validation has SOMETHING to align against.
    n_queries = max(0, int(getattr(args, "million_real_queries", 0) or 0))
    n_queries = min(n_queries, max(0, len(ref_manifest_rows) - 1))
    query_manifest_rows: list[dict[str, str]] = []
    query_list_paths: list[str] = []
    if n_queries > 0:
        # Exclude metagenomic / unidentified placeholders from the held-out
        # pool. ENA's read-runs lookup returns 0 hits for taxa starting with
        # "uncultured " / "unidentified " / "unclassified ", and the panel
        # then prints `[reads-mode] skip 'uncultured Pseudogymnoascus': ENA
        # filereport returned 0 runs` for every reads-mode pass. Filtering
        # them up front avoids that noise + ensures every held-out query has
        # a chance of matching public reads for read-level validation.
        def _is_metagenomic_placeholder(species: str) -> bool:
            sl = (species or "").strip().lower()
            return any(sl.startswith(prefix) for prefix in (
                "uncultured ", "unidentified ", "unclassified ", "environmental ",
            ))

        eligible_indices = [
            i for i, r in enumerate(ref_manifest_rows)
            if not _is_metagenomic_placeholder(r.get("clade_name", ""))
        ]
        if not eligible_indices:
            eligible_indices = list(range(len(ref_manifest_rows)))
        n_queries = min(n_queries, max(0, len(eligible_indices) - 1))
        if n_queries == 0:
            print(
                "      [skip] no non-metagenomic assemblies left after filter; skipping holdout selection",
                flush=True,
            )
            query_indices: list[int] = []
        else:
            # Seed-aware holdout. The old stride scheme picked the same 5
            # assemblies for every --seed value, so re-running with a
            # different seed couldn't surface selection bias. random.sample
            # with a Random(args.seed) instance gives reproducible-but-
            # varied selection: same seed -> same set, different seed ->
            # disjoint coverage.
            rng = random.Random(int(getattr(args, "seed", 0) or 0))
            query_indices = sorted(rng.sample(eligible_indices, n_queries))
        # Lookup helpers indexed by lineage so we can pick the closest sibling.
        by_genus: dict[str, list[int]] = defaultdict(list)
        by_family: dict[str, list[int]] = defaultdict(list)
        by_phylum: dict[str, list[int]] = defaultdict(list)
        for idx, r in enumerate(ref_manifest_rows):
            by_genus[r.get("genus", ".") or "."].append(idx)
            by_family[r.get("family", ".") or "."].append(idx)
            by_phylum[r.get("phylum", ".") or "."].append(idx)
        query_set = set(query_indices)

        def pick_benchmark_ref(qi: int) -> tuple[str, str]:
            qrow = ref_manifest_rows[qi]
            # Prefer a non-metagenomic candidate at every taxonomic level
            # before falling back to any candidate. A `uncultured ...`
            # benchmark ref pairs the held-out query against a MAG of unknown
            # assembly quality, which produces noisy SV truth — keep them as a
            # last resort rather than the first match.
            for bucket, key in (
                (by_genus, qrow.get("genus", ".") or "."),
                (by_family, qrow.get("family", ".") or "."),
                (by_phylum, qrow.get("phylum", ".") or "."),
            ):
                cands = [c for c in bucket.get(key, [])
                         if c != qi and c not in query_set]
                if not cands:
                    continue
                clean = [c for c in cands
                         if not _is_metagenomic_placeholder(
                             ref_manifest_rows[c].get("clade_name", ""))]
                if clean:
                    cand = ref_manifest_rows[clean[0]]
                    return cand.get("asm_name", "."), cand["fasta_path"]
                cand = ref_manifest_rows[cands[0]]
                return cand.get("asm_name", "."), cand["fasta_path"]
            # No genus/family/phylum sibling in the corpus. The previous
            # behavior fell through to ref_manifest_rows[0] — an alphabetically
            # first ascomycete against a phylum-isolated basidiomycete /
            # cryptomycotan query — producing zero homologous alignments and
            # silent NaN F1. Returning (".", ".") propagates honestly: the
            # query lands in query_manifest.tsv with an empty bench ref, and
            # the comparator runners + read-validation paths all skip cleanly
            # on benchmark_ref_fasta in {"", "."} (e.g. line 6498).
            sys.stderr.write(
                f"[holdout] {qrow.get('asm_name', '?')}: no genus/family/phylum "
                f"sibling in corpus (phylum={qrow.get('phylum', '.')}); "
                f"skipping benchmark ref selection. Comparators and "
                f"read-level validation will skip this query.\n"
            )
            return ".", "."

        for qi in query_indices:
            qrow = ref_manifest_rows[qi]
            bench_ref_asm, bench_ref = pick_benchmark_ref(qi)
            query_manifest_rows.append({
                "query_asm": qrow["asm_name"],
                "query_mode": "assembly",
                "path": qrow["fasta_path"],
                "scenario": "million_real",
                "lifestyle": ".",
                "architecture": ".",
                "benchmark_ref_asm": bench_ref_asm,
                "benchmark_ref_fasta": bench_ref,
                "phylum": qrow.get("phylum", "."),
                "class": qrow.get("class", "."),
                "order": qrow.get("order", "."),
                "family": qrow.get("family", "."),
                "genus": qrow.get("genus", "."),
                "species": qrow.get("clade_name", "."),
                "source": args.source,
                "instrument_platform": ".",
                "library_layout": ".",
                "run_accession": ".",
            })
            query_list_paths.append(qrow["fasta_path"])
        # Drop queries from the ref/hierarchy manifest so the index doesn't
        # see its own truth (would silently boost recall).
        keep_mask = [i not in query_set for i in range(len(ref_manifest_rows))]
        # Bug-fix: source_link_rows for the held-out asms were appended with
        # role="ref" inside _download_one (we don't know the holdout decision
        # at that point). Re-tag them to role="query" so public_data_links
        # / source_links honestly reflect the manifest. Without this, every
        # downstream join against role="query" misses the assembly query
        # rows entirely.
        held_out_asm_names = {ref_manifest_rows[qi]["asm_name"] for qi in query_indices}
        for sl in source_link_rows:
            if sl.get("role") == "ref" and sl.get("query_asm") in held_out_asm_names:
                sl["role"] = "query"
        ref_manifest_rows = [r for r, k in zip(ref_manifest_rows, keep_mask) if k]
        ref_list_paths = [p for p, k in zip(ref_list_paths, keep_mask) if k]
        print(
            f"      held out {len(query_manifest_rows)} assemblies as MycoSV-only "
            f"benchmark queries; {len(ref_manifest_rows)} remain in the index",
            flush=True,
        )

    # ── Optional: pull public ENA reads for each held-out query species so
    # ────── the million-real bench step can also exercise short-reads /
    # ────── long-reads modes (matching the per-panel real-data flow).
    # ────── Without this, only `benchmark_assembly/` is ever populated.
    # ────── Each (query, mode) pair appends a new row to the query manifest;
    # ────── benchmark filters by mode at run time so the assembly-mode
    # ────── invocation continues to see only the assembly rows.
    include_reads = bool(getattr(args, "million_real_include_reads", False))
    if include_reads and query_manifest_rows:
        modes_arg = getattr(args, "million_real_read_modes", "both") or "both"
        if modes_arg == "both":
            requested_read_modes = ["short-reads", "long-reads"]
        elif modes_arg in {"short-reads", "long-reads"}:
            requested_read_modes = [modes_arg]
        else:
            requested_read_modes = []
        runs_per_query = max(1, int(getattr(args, "million_real_read_runs_per_query", 1) or 1))
        ena_max_rows = max(1, int(getattr(args, "ena_max_rows_per_species", 200) or 200))
        queries_dir = cache_dir / "queries"
        queries_dir.mkdir(parents=True, exist_ok=True)
        # Snapshot the assembly rows; we mutate query_manifest_rows below.
        assembly_query_rows = list(query_manifest_rows)
        sys.stderr.write(
            f"[reads-mode] resolving ENA runs for {len(assembly_query_rows)} held-out "
            f"queries (modes={requested_read_modes}, runs/query={runs_per_query})\n"
        )
        for arow in assembly_query_rows:
            species = arow.get("species", ".") or "."
            if species in {".", ""}:
                continue
            try:
                ena_runs = fetch_ena_read_runs_by_species(species, max_rows=ena_max_rows)
            except Exception as exc:
                sys.stderr.write(
                    f"[reads-mode] ENA species lookup failed for {species!r}: {exc}\n"
                )
                continue
            if not ena_runs:
                sys.stderr.write(
                    f"[reads-mode] skip {species!r}: ENA filereport returned 0 runs\n"
                )
                continue
            sys.stderr.write(
                f"[reads-mode] {species!r}: ENA returned {len(ena_runs)} candidate runs\n"
            )
            for read_mode in requested_read_modes:
                pool_size = max(runs_per_query * 4, runs_per_query)
                pool_urls, pool_meta = select_ena_read_sources(ena_runs, read_mode, pool_size)
                if not pool_urls:
                    sys.stderr.write(
                        f"[reads-mode] {species!r} mode={read_mode}: "
                        f"no eligible runs after platform filter\n"
                    )
                    continue
                attempts = 0
                for meta in pool_meta:
                    if attempts >= runs_per_query:
                        break
                    run_acc = meta.get("run_accession", "na")
                    urls = selected_urls_from_ena_meta(meta)
                    if not urls:
                        continue
                    candidate_name = normalize_name(
                        f"{arow['query_asm']}_{read_mode}_{run_acc}"
                    )
                    try:
                        local_path = merge_sequence_sources(urls, queries_dir / candidate_name)
                    except Exception as exc:
                        sys.stderr.write(
                            f"[warn] reads download failed for {species} {read_mode} "
                            f"{run_acc}: {exc}\n"
                        )
                        continue
                    read_row = dict(arow)
                    read_row.update({
                        "query_asm": candidate_name,
                        "query_mode": read_mode,
                        "path": str(local_path),
                        "source": f"ena_{read_mode.replace('-', '_')}",
                        "instrument_platform": meta.get("instrument_platform", "."),
                        "library_layout": meta.get("library_layout", "."),
                        "run_accession": run_acc,
                    })
                    query_manifest_rows.append(read_row)
                    query_list_paths.append(str(local_path))
                    source_link_rows.append({
                        "query_asm": candidate_name,
                        "role": "query",
                        "query_mode": read_mode,
                        "source_type": "ena_read_run",
                        "source_accession": run_acc,
                        "source_url": meta.get("source_url", "."),
                        "local_path": str(local_path),
                        "species": meta.get("scientific_name", species),
                    })
                    attempts += 1
                    sys.stderr.write(
                        f"[reads-mode] {species!r} mode={read_mode}: "
                        f"added query {candidate_name}\n"
                    )
                if attempts == 0:
                    sys.stderr.write(
                        f"[warn] {species}: no usable ENA {read_mode} run after validation; skipping\n"
                    )

    hierarchy_manifest = out_dir / "hierarchy_manifest.tsv"
    write_tsv(
        hierarchy_manifest,
        ref_manifest_rows,
        ["asm_name", "phylum", "class", "order", "family", "genus", "clade_name", "clade_rank", "fasta_path"],
    )
    (out_dir / "ref_list.txt").write_text("\n".join(ref_list_paths) + "\n", encoding="utf-8")
    if query_manifest_rows:
        write_tsv(
            out_dir / "query_manifest.tsv",
            query_manifest_rows,
            ["query_asm", "query_mode", "path", "scenario", "lifestyle", "architecture",
             "benchmark_ref_asm", "benchmark_ref_fasta", "phylum", "class", "order",
             "family", "genus", "species", "source",
             "instrument_platform", "library_layout", "run_accession"],
        )
        (out_dir / "query_list.txt").write_text(
            "\n".join(query_list_paths) + "\n", encoding="utf-8"
        )
    source_link_columns = [
        "query_asm", "role", "query_mode", "source_type",
        "source_accession", "source_url", "local_path", "species",
    ]
    write_tsv(out_dir / "source_links.tsv", source_link_rows, source_link_columns)
    # Mirror under the per-panel name so downstream tools that look up
    # `prepared/public_data_links.tsv` (the per-panel convention) succeed
    # without an extra rename step.
    write_tsv(out_dir / "public_data_links.tsv", source_link_rows, source_link_columns)
    write_public_resource_links(out_dir / "public_resource_links.tsv")

    # selected_catalog.tsv mirrors the per-panel prepare output: a flat row
    # per NCBI assembly considered for the index, including the ones held
    # out as queries. Lets operators reproduce/audit the selection without
    # re-running fetch_ncbi_assembly_rows.
    write_tsv(
        out_dir / "selected_catalog.tsv",
        [
            {
                "panel_species": species_label_for_row(row),
                "assembly_accession": row.get("assembly_accession", ""),
                "organism_name": row.get("organism_name", ""),
                "assembly_level": row.get("assembly_level", ""),
                "refseq_category": row.get("refseq_category", ""),
                "version_status": row.get("version_status", ""),
                "seq_rel_date": row.get("seq_rel_date", ""),
                "ftp_path": row.get("ftp_path", ""),
                "source_catalog": row.get("_catalog_source", args.source),
            }
            for row in selected
        ],
        [
            "panel_species", "assembly_accession", "organism_name", "assembly_level",
            "refseq_category", "version_status", "seq_rel_date", "ftp_path", "source_catalog",
        ],
    )

    # benchmark_reference_map.tsv: one row per query asm, pointing at the
    # benchmark-ref fasta that read-level validation will align against. We
    # already computed `benchmark_ref_fasta` per held-out query above, plus
    # any reads-mode children inherit the same ref via dict(arow), so this
    # is a direct projection of query_manifest_rows.
    benchmark_map_rows = [
        {
            "query_asm": qrow.get("query_asm", "."),
            "benchmark_ref_asm": qrow.get("benchmark_ref_asm", "."),
            "benchmark_ref_fasta": qrow.get("benchmark_ref_fasta", "."),
            "species": qrow.get("species", "."),
        }
        for qrow in query_manifest_rows
    ]
    if benchmark_map_rows:
        write_tsv(
            out_dir / "benchmark_reference_map.tsv",
            benchmark_map_rows,
            ["query_asm", "benchmark_ref_asm", "benchmark_ref_fasta", "species"],
        )

    # phenotypic_metadata.json + ecological_traits.tsv mirror the per-panel
    # prepare path: pull BioSample phenotype attributes for the selected
    # rows and join FungalTraits onto the held-out species set so the
    # biology analyzer / visualization report pick up substrate /
    # primary_lifestyle / etc. without manual flags. Caches live under
    # data_cache so subsequent runs are no-ops; first run pays the network.
    phenotype_meta: dict[str, dict[str, str]] = {}
    ecological_rows_written = 0
    if getattr(args, "million_real_phenotypes", False):
        biosample_ids = sorted({row.get("biosample", "") for row in selected if row.get("biosample")})
        phenotype_cache_path = cache_dir / "phenotypic_metadata.json"
        if biosample_ids:
            try:
                phenotype_meta = fetch_ncbi_biosample_phenotypes(
                    biosample_ids, cache_path=phenotype_cache_path,
                )
                print(
                    f"      phenotype: cached {len(phenotype_meta)} BioSample records "
                    f"-> {phenotype_cache_path}",
                    flush=True,
                )
            except Exception as exc:
                sys.stderr.write(
                    f"[warn] phenotype lookup failed: {type(exc).__name__}: {exc}; "
                    f"continuing without phenotypic_metadata.json\n"
                )
                phenotype_meta = {}
        if phenotype_meta:
            (out_dir / "phenotypic_metadata.json").write_text(
                json.dumps(phenotype_meta, indent=2, sort_keys=True), encoding="utf-8"
            )
        try:
            fungaltraits_csv = fetch_fungaltraits_table(cache_dir)
        except Exception as exc:
            sys.stderr.write(
                f"[warn] fungaltraits fetch failed: {type(exc).__name__}: {exc}; "
                f"continuing without ecological_traits.tsv\n"
            )
            fungaltraits_csv = None
        species_to_query_asms: dict[str, list[str]] = defaultdict(list)
        for qrow in query_manifest_rows:
            sp = qrow.get("species") or "."
            species_to_query_asms[sp].append(qrow.get("query_asm", "."))
        if species_to_query_asms:
            ecological_rows_written = write_ecological_summary_tsv(
                fungaltraits_csv,
                species_to_query_asms,
                out_dir / "ecological_traits.tsv",
            )
            if ecological_rows_written:
                print(
                    f"      ecological_traits: joined {ecological_rows_written} "
                    f"trait records -> {out_dir / 'ecological_traits.tsv'}",
                    flush=True,
                )

    # gene_annotations.tsv: aggregate every GFF/GBFF source we materialised
    # during the parallel download. Same alias-expansion trick as the
    # per-panel prepare path so analyze_new_biology_candidates' lookup
    # against (query_asm, ref_contig) always hits a row when the gene
    # exists in the ref annotation, regardless of whether the call was
    # written under the manifest asm_name or the FASTA basename.
    gene_annotations_count = 0
    if gff_pairs:
        gene_annotations_path = out_dir / "gene_annotations.tsv"
        asm_aliases: dict[str, set[str]] = defaultdict(set)
        for ref_row in ref_manifest_rows:
            asm = ref_row.get("asm_name", "")
            if not asm:
                continue
            asm_aliases[asm].add(asm)
            fasta_path = ref_row.get("fasta_path", "")
            if fasta_path:
                basename = Path(fasta_path).name
                if basename.endswith(".gz"):
                    basename = basename[:-3]
                asm_aliases[asm].add(basename)
        for q_row in query_manifest_rows:
            q_asm = q_row.get("query_asm", "")
            q_path = q_row.get("path", "")
            if q_asm and q_path:
                basename = Path(q_path).name
                if basename.endswith(".gz"):
                    basename = basename[:-3]
                asm_aliases[q_asm].add(q_asm)
                asm_aliases[q_asm].add(basename)
        ref_to_queries: dict[str, list[str]] = defaultdict(list)
        for bm in benchmark_map_rows:
            ref_asm = bm.get("benchmark_ref_asm", "")
            q_asm = bm.get("query_asm", "")
            if ref_asm and q_asm:
                ref_to_queries[ref_asm].append(q_asm)
        print(
            f"      gene_annotations: streaming {len(gff_pairs)} GFF/GBFF sources "
            f"-> {gene_annotations_path}",
            flush=True,
        )
        gene_annotations_count = stream_gene_annotations_to_tsv(
            gene_annotations_path, gff_pairs, asm_aliases, ref_to_queries,
        )
        if gene_annotations_count:
            print(
                f"      gene_annotations: wrote {gene_annotations_count} gene records "
                f"(from {len(gff_pairs)} GFF/GBFF sources) -> {gene_annotations_path}",
                flush=True,
            )
        else:
            sys.stderr.write(
                f"[gene-annot] no gene records parsed from {len(gff_pairs)} annotation source(s); "
                f"skipping gene_annotations.tsv\n"
            )
            try:
                gene_annotations_path.unlink()
            except FileNotFoundError:
                pass

    # Step 4: build the real routing index by invoking the MycoSV binary, then
    # pad with synthetic decoys up to --target-centroids if requested.
    print(f"[4/4] Building real routing index via MycoSV binary -> {index_dir}", flush=True)
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

    # `queries_held_out` historically counted every row in the query manifest
    # (assembly + per-mode reads children), which conflated "species held
    # out" with "query rows produced". Split them so the summary is honest:
    # query_assemblies_held_out = unique assemblies, query_rows_total = the
    # rows the benchmark step iterates over.
    held_out_assembly_rows = sum(1 for q in query_manifest_rows if q.get("query_mode") == "assembly")
    summary = {
        "out_dir": str(out_dir),
        "data_cache_dir": str(cache_dir),
        "assembly_summary_cache": ";".join(str(p) for p in assembly_summary_caches),
        "assembly_summary_caches": [str(p) for p in assembly_summary_caches],
        "source": args.source,
        "max_assemblies_requested": args.max_assemblies,
        "assemblies_downloaded": download_count,
        "gff_sources_downloaded": len(gff_pairs),
        "gene_annotations_rows_written": gene_annotations_count,
        "phenotypic_records_cached": len(phenotype_meta),
        "ecological_rows_written": ecological_rows_written,
        "query_assemblies_held_out": held_out_assembly_rows,
        "query_rows_total": len(query_manifest_rows),
        # Back-compat alias for older log scrapers that key on this field.
        "queries_held_out": len(query_manifest_rows),
        "refs_in_index": len(ref_manifest_rows),
        "target_centroids": args.target_centroids,
        "seed": args.seed,
        "hierarchy_manifest": str(hierarchy_manifest),
        "index_dir": str(index_dir),
        "registry_dir": str(registry_dir),
        **info,
    }
    summary_path = out_dir / "prepare_million_real_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    # Mirror under prepare_summary.json so tooling that expects the
    # per-panel filename also finds the run summary in the million-real
    # output dir.
    (out_dir / "prepare_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )

    print(
        f"million_real_ready\tassemblies={download_count}"
        f"\tqueries_held_out={len(query_manifest_rows)}"
        f"\tcentroids_real={info['real_centroids']}"
        f"\tcentroids_total={info['total_centroids']}\tindex_dir={index_dir}\tsummary={summary_path}",
        flush=True,
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
    smr.add_argument(
        "--source",
        choices=NCBI_SOURCE_CHOICES,
        default=NCBI_BEST_SOURCE,
        help=(
            "NCBI assembly source. ncbi-best combines RefSeq + GenBank, "
            "deduplicates paired GCF/GCA assemblies, and keeps the best "
            "latest/highest-level/full/curated record."
        ),
    )
    smr.add_argument("--max-assemblies", type=int, default=1000, help="Cap on real fungal assemblies to download. Use 0 for unlimited (not recommended for first runs).")
    smr.add_argument("--min-assembly-level", default="scaffold", choices=["contig", "scaffold", "chromosome", "complete genome"])
    smr.add_argument("--latest-only", action="store_true", help="Only keep version_status=latest rows.")
    smr.add_argument("--target-centroids", type=int, default=1_000_000, help="Total centroids in the routing store after padding with synthetic decoys; set to 0 to skip padding.")
    smr.add_argument("--seed", type=int, default=1)
    smr.add_argument("--threads", type=int, default=32)
    smr.add_argument("--max-clade-genomes", type=int, default=32)
    smr.add_argument("--binary-path", type=Path, default=DEFAULT_BIN)
    smr.add_argument("--force-rebuild", action="store_true")
    smr.add_argument("--data-cache-dir", type=Path, default=DEFAULT_DATA_CACHE, help="Shared directory for downloaded FASTA/GFF/FASTQ files and metadata caches; reused across runs to avoid re-downloading. Defaults to data_cache/ next to this script.")
    smr.add_argument(
        "--million-real-queries", type=int, default=0,
        help=(
            "Hold out N assemblies as MycoSV-only benchmark queries (sampled "
            "stride-uniformly across phyla). Writes query_manifest.tsv + "
            "query_list.txt next to the index so the standard `benchmark` "
            "subcommand can run end-to-end (SV calls, TE classification, "
            "biology candidates, visualization) without algorithmic "
            "comparators. Default 0 = index-only, no held-out queries."
        ),
    )
    # Reads-mode queries for the held-out species: pull public ENA runs so
    # the million-real bench can also exercise short-reads (Illumina) and
    # long-reads (PacBio HiFi / ONT R10.4 etc.) modes — without these,
    # only benchmark_assembly/ ever gets populated.  The per-mode queries
    # share the held-out assembly's benchmark_ref_fasta (closest sibling),
    # so per-query truth alignment still works.  Off by default to keep
    # prep wall-time bounded; toggle via --million-real-include-reads.
    smr.add_argument(
        "--million-real-include-reads", action="store_true",
        help="For each held-out query species, also resolve and download "
             "ENA reads to seed reads-mode benchmark queries.",
    )
    smr.add_argument(
        "--million-real-read-modes", default="both",
        choices=["both", "short-reads", "long-reads"],
        help="Which reads modes to materialise when --million-real-include-reads "
             "is set (default both).",
    )
    smr.add_argument(
        "--million-real-read-runs-per-query", type=int, default=1,
        help="Number of ENA runs to fetch per (held-out query, reads mode); "
             "the FASTQs are merged into a single per-query input (default 1).",
    )
    smr.add_argument(
        "--ena-max-rows-per-species", type=int, default=200,
        help="Maximum read_run rows to pull from ENA per species before "
             "filtering by platform.",
    )
    # GFF/GBFF -> gene_annotations.tsv next to the manifest, mirroring the
    # per-panel `prepare` default. NCBI hosts genomic.gff.gz alongside the
    # FASTA so the marginal per-assembly cost is small; for ~2k assemblies
    # this is on the order of minutes once warmed. Pass --no-million-real-
    # download-gff if bandwidth- or wall-time-constrained.
    smr.add_argument(
        "--million-real-download-gff", dest="million_real_download_gff",
        action="store_true", default=True,
        help="Download GFF (or GBFF fallback) annotations alongside reference "
             "FASTA so the prepared dir gets a gene_annotations.tsv (default: on).",
    )
    smr.add_argument(
        "--no-million-real-download-gff", dest="million_real_download_gff",
        action="store_false",
        help="Skip GFF/GBFF download (disables auto gene_annotations.tsv).",
    )
    # FungalTraits + BioSample phenotypes are species-level lookups used by
    # the biology analyzer and the visualization report. The CSVs are cached
    # under data_cache/ so the first run pays the network cost and later
    # runs are no-ops; default ON keeps million-real on par with the per-
    # panel prepare path.
    smr.add_argument(
        "--million-real-phenotypes", dest="million_real_phenotypes",
        action="store_true", default=True,
        help="Fetch BioSample phenotype attributes + FungalTraits and emit "
             "phenotypic_metadata.json / ecological_traits.tsv (default: on).",
    )
    smr.add_argument(
        "--no-million-real-phenotypes", dest="million_real_phenotypes",
        action="store_false",
        help="Skip BioSample/FungalTraits enrichment (disables phenotypic_"
             "metadata.json + ecological_traits.tsv).",
    )

    spp = sub.add_parser("prepare", help="Download a real fungal panel and write MycoSV-ready manifests.")
    spp.add_argument("--out-dir", type=Path, required=True)
    spp.add_argument(
        "--source",
        choices=NCBI_SOURCE_CHOICES,
        default=NCBI_BEST_SOURCE,
        help=(
            "NCBI assembly source. ncbi-best combines RefSeq + GenBank, "
            "deduplicates paired GCF/GCA assemblies, and keeps the best "
            "latest/highest-level/full/curated record."
        ),
    )
    spp.add_argument("--panel", "--panels", dest="panels", action="append", choices=sorted(PANEL_PRESETS), default=[])
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
    # GFF defaults to ON so that downstream gene-near-breakpoint analysis has
    # something to join to without a second `prepare` round-trip. NCBI hosts
    # GFF.gz next to FASTA.gz at the same FTP path, so the marginal cost is
    # small. Pass --no-download-gff to skip when bandwidth-constrained.
    spp.add_argument("--download-gff", dest="download_gff", action="store_true", default=True,
                     help="Download GFF annotations alongside reference FASTA "
                          "(default: on; used to build gene_annotations.tsv).")
    spp.add_argument("--no-download-gff", dest="download_gff", action="store_false",
                     help="Skip GFF download (disables auto gene_annotations.tsv).")
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
    spp.add_argument("--data-cache-dir", type=Path, default=DEFAULT_DATA_CACHE, help="Shared directory for downloaded FASTA/GFF/FASTQ files and metadata caches; reused across runs to avoid re-downloading. Defaults to data_cache/ next to this script.")

    sb = sub.add_parser("benchmark", help="Run MycoSV on a prepared real-data panel and compare to exact normalized truth/query-aware callsets.")
    sb.add_argument("--prepared-dir", type=Path, required=True)
    sb.add_argument("--out-dir", type=Path, required=True)
    sb.add_argument("--binary-path", type=Path, default=DEFAULT_BIN)
    sb.add_argument("--force-rebuild", action="store_true")
    sb.add_argument("--mode", default="assembly", choices=["assembly", "short-reads", "long-reads", "auto"])
    sb.add_argument("--threads", type=int, default=32)
    # Lower than the millon-real default of 32: hot clades with many genomes
    # can spike the binary's per-clade graph build past a 12 GiB cgroup. 8 is
    # safe for 12 GiB; raise on larger nodes via --max-clade-genomes.
    sb.add_argument("--max-clade-genomes", type=int, default=8,
                    help="Per-clade cap on genomes loaded into RAM by the MycoSV binary "
                         "during the hierarchical index build (default: 8, safe for ~12 GiB).")
    sb.add_argument("--run-all-comparators", action="store_true",
                    help="Auto-enable every --run-X flag whose tool binaries are detected on PATH "
                         "(or in the project conda env). The most robust way to get truth-set rows "
                         "in exact_benchmark_summary.tsv without naming each comparator individually.")
    sb.add_argument("--mycosv-only", action="store_true",
                    help="Run MycoSV only — no algorithmic comparators (minigraph/syri/cactus/"
                         "svim-asm/anchorwave/sniffles/cuteSV/svim/Delly/Manta) are launched, "
                         "and the per-mode mandatory-baseline forcing in benchmark_real_data is "
                         "disabled. Used by the million-real flow where comparators are out of "
                         "scope. exact_benchmark_summary.tsv falls back to the no_comparator "
                         "placeholder rows; biology candidates / TE classification / read-level "
                         "validation still run.")
    sb.add_argument("--reuse-index-dir", type=Path,
                    help="Reuse a pre-built MycoSV routing index instead of rebuilding from "
                         "the prepared dir's hierarchy_manifest.tsv. Path must contain "
                         "routing_manifest.tsv (e.g. ${MILLION_REAL_DIR}/index from "
                         "prepare-million-real).")
    sb.add_argument("--reuse-registry-dir", type=Path,
                    help="Reuse a pre-built clade registry alongside --reuse-index-dir.")
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
    sb.add_argument("--max-comparator-short-reads", type=int, default=500000,
                    help="Cap short-read FASTQ records used by external comparators "
                         "and read validation. 0 disables subsetting (default: 500000).")
    sb.add_argument("--max-comparator-long-reads", type=int, default=200000,
                    help="Cap long-read FASTQ records used by external comparators "
                         "and read validation. 0 disables subsetting (default: 200000). "
                         "Bumped from 20000 because at the lower cap svim / sniffles / "
                         "cuteSV produce zero SV calls on ~60–80 Mbp fungal genomes "
                         "(0.02–2× coverage) — which is why the headline TSVs were "
                         "dominated by status=no_truth.")
    sb.add_argument("--mycosv-use-full-reads", action="store_true",
                    help="In read modes, pass the original full FASTQ to MycoSV "
                         "while comparators use the capped subset. Off by default "
                         "because million-real public FASTQs can be multi-GB and "
                         "previously caused empty VCFs after timeout/OOM kills.")
    sb.add_argument("--max-assembly-query-contigs", type=int, default=0,
                    help="Skip assembly-mode query FASTAs with more than this many "
                         "records before launching MycoSV/comparators. 0 disables "
                         "the filter (default: 0). Skipped rows are written to "
                         "SKIPPED_ASSEMBLY_QUERIES.tsv.")
    sb.add_argument("--max-assembly-query-bp", type=int, default=0,
                    help="Skip assembly-mode query FASTAs above this total bp before "
                         "launching MycoSV/comparators. 0 disables the filter "
                         "(default: 0).")
    sb.add_argument("--skip-input-preflight", action="store_true",
                    help="Skip the cheap FASTA/FASTQ readability preflight. "
                         "Normally keep this on: it catches corrupt cached gzip "
                         "files and manifest path mixups before expensive MycoSV/"
                         "comparator work starts.")
    sb.add_argument("--benchmark-ref-cap", type=int, default=512,
                    help="Maximum benchmark reference FASTAs to pass to MycoSV "
                         "after adding per-query refs plus genus/family/order/"
                         "class neighbors (default: 512). Million-real runs can "
                         "raise this to recover expected call volume.")
    sb.add_argument("--max-benchmark-queries", type=int, default=0,
                    help="Run only the first N mode-matched query rows from "
                         "query_manifest.tsv. 0 keeps all queries. Useful for "
                         "short smoke jobs that only need MycoSV SV counts.")
    sb.add_argument("--benchmark-query-genera", default="",
                    help="Comma/space separated target fungal groups. Benchmark "
                         "selects one mode-matched query row per requested group "
                         "when present in query_manifest.tsv; missing assembly "
                         "groups are synthesized from hierarchy_manifest.tsv "
                         "when possible. Writes REQUESTED_QUERY_GROUPS.tsv and "
                         "reports missing groups. "
                         "Also accepts common misspellings and mycorrhiza as a "
                         "Rhizophagus/Glomus/Glomeromycetes-style biology group.")
    sb.add_argument("--normalized-other", action="append", default=[], metavar="LABEL=PATH", help="Additional normalized TSV callsets to benchmark against. TSVs may use query or reference coordinates via a coord_space column.")
    sb.add_argument("--other-vcf", action="append", default=[], metavar="LABEL=PATH", help="Additional reference-coordinate VCF comparator output, best for single-query or pairwise benchmark runs.")
    sb.add_argument("--mycosv-arg", action="append", default=[], help="Extra argument to pass through to the MycoSV binary; may be used multiple times.")
    sb.add_argument("--minigraph-arg", action="append", default=[], help="Extra argument to pass through to minigraph runs; may be used multiple times.")
    sb.add_argument("--pggb-arg", action="append", default=[], help="Extra argument to pass through to pggb; may be used multiple times.")
    # 90% identity is wfmash's tutorial default for human haplotype panels;
    # cross-strain fungal pairs routinely sit at 85-95% nucleotide identity,
    # so the stricter default produced rc=2 (zero homologous mappings) on
    # most yeast and basidiomycete panels.  80 is the value PanSN tutorials
    # recommend for "diverged but related" inputs and recovers those pairs.
    sb.add_argument("--pggb-identity", default="80")
    sb.add_argument("--pggb-segment-len", default="5k")
    sb.add_argument("--expression-tsv", type=Path)
    sb.add_argument("--gene-annotations-tsv", type=Path)
    sb.add_argument("--ancestral-tsv", type=Path)
    # Read-level (FASTQ-anchored) independent validation of candidate calls.
    # On by default — algorithm comparators inherit assembly artefacts, so
    # raw-read validated rows are the preferred fungal validation basis.
    sb.add_argument("--validate-with-reads", dest="validate_with_reads",
                    action="store_true", default=True,
                    help="Re-anchor candidate SV calls in raw query reads / "
                         "contigs via samtools split-read counting "
                         "(default: on; needs minimap2 + samtools).")
    sb.add_argument("--no-validate-with-reads", dest="validate_with_reads",
                    action="store_false")
    sb.add_argument("--read-validation-min-support", type=int, default=2,
                    help="Minimum split/clipped reads required at the "
                         "breakpoint to keep an SV in the read-validated "
                         "validated callset (default: 2). The C++ MycoSV pipeline "
                         "clusters at SUPPORT=2; defaulting validation to 3 "
                         "auto-failed ~50–90 %% of mycosv short-read calls "
                         "even though each call already had ≥2 cluster reads.")
    sb.add_argument("--read-validation-flank-bp", type=int, default=250,
                    help="Window around each breakpoint where supporting "
                         "split/clipped reads are counted (default: 250 bp).")
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
