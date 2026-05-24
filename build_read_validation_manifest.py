#!/usr/bin/env python3
"""Build a RAW_READ_VALIDATION_TSV for the fungal benchmark.

For every assembly query in query_manifest.tsv, resolve a public ENA read run
(preferring long reads), download the FASTQ to the local reads cache, cap to a
byte budget so a single PacBio/ONT run does not consume terabytes, and emit a
TSV with the columns expected by load_raw_read_validation_manifest():
  query_asm, path, query_mode, instrument_platform, library_layout, run_accession

Queries with no public reads are recorded in an availability log alongside the
manifest (and omitted from the manifest, which causes the downstream pipeline
to leave them as `validation_unavailable` rather than failing).
"""
from __future__ import annotations

import argparse
import csv
import os
import shutil
import sys
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_real_fungal_benchmark import (  # noqa: E402  (sibling import)
    ena_filereport_url,
    fetch_ena_read_runs,
    fetch_ena_read_runs_by_species,
    normalise_download_url,
    select_ena_read_sources,
    sequence_kind_from_name,
)

MANIFEST_FIELDS = [
    "query_asm",
    "path",
    "query_mode",
    "instrument_platform",
    "library_layout",
    "run_accession",
    "source_url",
    "status",
]


def _pick_runs_for_query(row: dict[str, str], *, max_runs: int) -> tuple[list[str], list[dict[str, str]], str]:
    """Return (urls, ena_meta_rows, mode) for a query row, preferring long reads."""
    run_acc = (row.get("run_accession") or "").strip()
    candidates: list[dict[str, str]] = []
    if run_acc and run_acc != ".":
        try:
            candidates = fetch_ena_read_runs(run_acc)
        except Exception as exc:  # network or parse error
            sys.stderr.write(f"[reads] {row['query_asm']}: run lookup {run_acc} failed: {exc}\n")
    if not candidates:
        species = (row.get("species") or row.get("scientific_name") or "").strip()
        if species and species != ".":
            try:
                candidates = fetch_ena_read_runs_by_species(species, max_rows=200)
            except Exception as exc:
                sys.stderr.write(f"[reads] {row['query_asm']}: species lookup '{species}' failed: {exc}\n")
    if not candidates:
        return [], [], "long-reads"

    for mode in ("long-reads", "short-reads"):
        urls, meta = select_ena_read_sources(candidates, mode, max_runs)
        if urls:
            return urls, meta, mode
    return [], [], "long-reads"


def _download_capped(url: str, dest: Path, *, max_bytes: int) -> tuple[bool, int]:
    """Stream `url` into `dest`, stopping after max_bytes. Idempotent.

    Returns (success, bytes_written). If the file already exists with non-zero
    size we treat it as already-downloaded and skip — re-runs are cheap.
    """
    if dest.exists() and dest.stat().st_size > 0:
        return True, dest.stat().st_size
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    written = 0
    req = urllib.request.Request(url, headers={"User-Agent": "mycosv-bench/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=120) as resp, tmp.open("wb") as out_fh:
            while True:
                chunk = resp.read(1 << 20)
                if not chunk:
                    break
                if max_bytes > 0 and written + len(chunk) > max_bytes:
                    out_fh.write(chunk[: max_bytes - written])
                    written = max_bytes
                    break
                out_fh.write(chunk)
                written += len(chunk)
        tmp.rename(dest)
        return True, written
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        sys.stderr.write(f"[reads] download failed for {url}: {exc}\n")
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass
        return False, written


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--query-manifest", type=Path, required=True,
                    help="prepared/query_manifest.tsv")
    ap.add_argument("--reads-cache", type=Path, required=True,
                    help="directory where FASTQs are stored (e.g. data_cache/raw_reads)")
    ap.add_argument("--out-tsv", type=Path, required=True,
                    help="output read-validation manifest TSV")
    ap.add_argument("--max-bases", type=int, default=300_000_000,
                    help="approximate cap per query, treated as bytes; default 300MB "
                         "(~1-2M PacBio HiFi or ~600k ONT reads)")
    ap.add_argument("--max-runs-per-query", type=int, default=1,
                    help="cap number of ENA runs fetched per query (default 1)")
    ap.add_argument("--force", action="store_true",
                    help="re-emit manifest even if it exists")
    args = ap.parse_args()

    if args.out_tsv.exists() and not args.force:
        sys.stderr.write(f"[reads] manifest already exists at {args.out_tsv}; pass --force to rebuild\n")
        return 0

    if not args.query_manifest.exists():
        sys.stderr.write(f"[reads] query manifest not found: {args.query_manifest}\n")
        return 2

    args.reads_cache.mkdir(parents=True, exist_ok=True)
    rows_out: list[dict[str, str]] = []
    availability_log: list[dict[str, str]] = []

    with args.query_manifest.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        queries = [r for r in reader if (r.get("query_mode") or "").strip() == "assembly"]

    sys.stderr.write(f"[reads] {len(queries)} assembly queries to resolve\n")
    for q in queries:
        qasm = (q.get("query_asm") or "").strip()
        if not qasm:
            continue
        urls, meta, mode = _pick_runs_for_query(q, max_runs=args.max_runs_per_query)
        if not urls:
            sys.stderr.write(f"[reads] {qasm}: no public ENA reads found\n")
            availability_log.append({
                "query_asm": qasm,
                "status": "reads_unavailable",
                "reason": "no_ena_run_for_species_or_accession",
            })
            continue
        # Download the first FASTQ URL only (multi-file pair handled by mycosv).
        url = urls[0]
        run_acc = meta[0].get("run_accession", ".") if meta else "."
        fname = Path(urllib.request.url2pathname(url.split("/")[-1])).name or f"{qasm}.fastq.gz"
        dest = args.reads_cache / qasm / fname
        ok, nbytes = _download_capped(normalise_download_url(url), dest, max_bytes=args.max_bases)
        if not ok:
            availability_log.append({
                "query_asm": qasm,
                "status": "reads_unavailable",
                "reason": f"download_failed:{url}",
            })
            continue
        rows_out.append({
            "query_asm": qasm,
            "path": str(dest),
            "query_mode": mode,
            "instrument_platform": meta[0].get("instrument_platform", "."),
            "library_layout": meta[0].get("library_layout", "."),
            "run_accession": run_acc,
            "source_url": url,
            "status": f"downloaded:{nbytes}B",
        })
        availability_log.append({
            "query_asm": qasm,
            "status": "available",
            "reason": f"{run_acc}:{nbytes}B",
        })
        sys.stderr.write(f"[reads] {qasm}: {run_acc} {mode} {nbytes}B -> {dest}\n")

    # Write manifest (only rows with usable paths — load_raw_read_validation_manifest
    # skips entries without a path).
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_tsv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=MANIFEST_FIELDS, delimiter="\t")
        writer.writeheader()
        for row in rows_out:
            writer.writerow(row)

    # Write availability log next to the manifest.
    log_path = args.out_tsv.with_suffix(".availability.tsv")
    with log_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["query_asm", "status", "reason"], delimiter="\t")
        writer.writeheader()
        for row in availability_log:
            writer.writerow(row)

    avail = sum(1 for r in availability_log if r["status"] == "available")
    unavail = len(availability_log) - avail
    sys.stderr.write(f"[reads] manifest written to {args.out_tsv}: {avail} available, {unavail} unavailable\n")
    sys.stderr.write(f"[reads] availability log: {log_path}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
