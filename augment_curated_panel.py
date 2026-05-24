#!/usr/bin/env python3
"""Augment a prepare-million-real query_manifest with curated + supplemental
panel accessions.

Inputs:
  --prepared-dir       directory produced by `prepare-million-real`
  --curated-tsv        user-curated Pezizomycotina accessions (Aspergillus/
                       Trichoderma/Penicillium/Talaromyces; pinned as queries)
  --supplemental-tsv   one rep per AMF/EMF/yeast/mushroom/truffle/endophyte/
                       other-filamentous species (auto-picked from catalog)
  --catalog            ncbi-genbank assembly_summary file used to fill in
                       taxonomy fields and local FASTA paths for the injected
                       queries.

Effect:
  - REWRITES prepared-dir/query_manifest.tsv so the queries are exactly
    `curated + supplemental` (deduplicated by accession).
  - For each injected query, picks the best benchmark_ref_asm using:
      conspecific (curated/supplemental pool, then catalog), else
      congeneric, else
      whatever prepare-million-real already chose.
  - REWRITES prepared-dir/benchmark_reference_map.tsv to match.
  - Writes injected_panel.tsv (audit log of which species got which ref).
  - Returns the count of queries written (printed on the last line of stdout
    so the bash wrapper can use it to size the SLURM array).

Anything else under prepared-dir (index/, registry/, gene_annotations.tsv,
source_links.tsv, …) is left untouched so the routing index and gene tables
that prepare-million-real built can be reused.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

REFS_CACHE = Path("/mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/AMF/scale/data_cache/refs")


def parse_catalog(path: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 20:
                continue
            acc = cols[0]
            rows[acc] = {
                "accession": acc,
                "level": cols[11],
                "organism": cols[7],
                "taxid": cols[5],
                "species_taxid": cols[6],
                "asm_name": cols[15],
                "ftp_path": cols[19],
            }
    return rows


def normalize_query_asm(acc: str) -> str:
    """Convert GCA_xxx.y → GCA_xxx_y to match prepare-million-real convention."""
    return acc.replace(".", "_")


def local_path_for(acc: str, asm_name: str) -> Path:
    """Map an NCBI accession + asm_name to the local data_cache/refs path."""
    return REFS_CACHE / f"{acc}_{asm_name}_genomic.fna.gz"


def load_user_panel(path: Path) -> list[dict[str, str]]:
    rows = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            acc = (r.get("accession") or r.get("input_accession") or "").strip()
            if not acc:
                continue
            rows.append({
                "accession": acc,
                "organism": (r.get("organism") or "").strip(),
                "asm_name": (r.get("assembly_name") or r.get("asm_name") or "").strip(),
                "guild": (r.get("guild") or "").strip(),
            })
    return rows


def species_from_organism(org: str) -> str:
    parts = org.split()
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
    return org


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prepared-dir", type=Path, required=True)
    ap.add_argument("--curated-tsv", type=Path, required=True)
    ap.add_argument("--supplemental-tsv", type=Path, required=True)
    ap.add_argument("--catalog", type=Path,
                    default=Path("data_cache/assembly_summaries/ncbi-genbank_fungi_assembly_summary.txt"))
    ap.add_argument("--catalog-refseq", type=Path,
                    default=Path("data_cache/assembly_summaries/ncbi-refseq_fungi_assembly_summary.txt"))
    args = ap.parse_args()

    prep_qm = args.prepared_dir / "query_manifest.tsv"
    prep_brm = args.prepared_dir / "benchmark_reference_map.tsv"
    if not prep_qm.exists():
        sys.stderr.write(f"[augment] prepared query_manifest not found at {prep_qm}\n")
        return 2

    # Load prepared manifest so we can reuse its taxonomy fields and existing
    # benchmark refs as fallbacks.
    prep_rows: list[dict[str, str]] = []
    with prep_qm.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        prep_rows = list(reader)
    prep_by_species: dict[str, dict[str, str]] = {}
    prep_by_genus: dict[str, list[dict[str, str]]] = defaultdict(list)
    for r in prep_rows:
        sp = species_from_organism(r.get("species") or "")
        if sp:
            prep_by_species.setdefault(sp, r)
            g = sp.split()[0]
            prep_by_genus[g].append(r)

    # Load catalogs for taxonomy fill-in.
    catalog = parse_catalog(args.catalog)
    catalog.update(parse_catalog(args.catalog_refseq))

    # Compose user panel = curated + supplemental, dedup by accession (keep first occurrence).
    seen_acc: set[str] = set()
    user_panel: list[dict[str, str]] = []
    for src in (args.curated_tsv, args.supplemental_tsv):
        for row in load_user_panel(src):
            if row["accession"] in seen_acc:
                continue
            seen_acc.add(row["accession"])
            user_panel.append(row)
    sys.stderr.write(f"[augment] user-panel size after dedup: {len(user_panel)}\n")

    # Build injected query rows.
    out_rows: list[dict[str, str]] = []
    audit_rows: list[dict[str, str]] = []
    by_species_user: dict[str, dict[str, str]] = {}
    for u in user_panel:
        sp = species_from_organism(u["organism"])
        by_species_user.setdefault(sp, u)

    def pick_benchmark_ref(query_sp: str, query_acc: str) -> tuple[str, str, str]:
        """Return (ref_asm_id_underscore, ref_fasta_path, decision_note)."""
        # 1. conspecific from user panel (not self)
        for sp2, u in by_species_user.items():
            if sp2 == query_sp and u["accession"] != query_acc:
                meta = catalog.get(u["accession"])
                if meta:
                    return (normalize_query_asm(u["accession"]),
                            str(local_path_for(u["accession"], meta["asm_name"])),
                            f"conspecific_user:{u['accession']}")
        # 2. conspecific from catalog
        candidates = [m for m in catalog.values() if species_from_organism(m["organism"]) == query_sp and m["accession"] != query_acc]
        if candidates:
            # Prefer best level
            rank = {"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}
            candidates.sort(key=lambda r: -rank.get(r["level"], 0))
            m = candidates[0]
            return (normalize_query_asm(m["accession"]),
                    str(local_path_for(m["accession"], m["asm_name"])),
                    f"conspecific_catalog:{m['accession']}")
        # 3. fall back to prepare-million-real's pick for this species
        prep = prep_by_species.get(query_sp)
        if prep:
            return (prep.get("benchmark_ref_asm", ""),
                    prep.get("benchmark_ref_fasta", ""),
                    "prepare_species_match")
        # 4. fall back to any prepare pick in the same genus
        genus = query_sp.split()[0]
        prep_g = prep_by_genus.get(genus)
        if prep_g:
            prep = prep_g[0]
            return (prep.get("benchmark_ref_asm", ""),
                    prep.get("benchmark_ref_fasta", ""),
                    f"prepare_genus_match:{prep.get('species', '?')}")
        # 5. last resort: use the catalog's first congeneric
        congen = [m for m in catalog.values() if m["organism"].startswith(genus + " ") and m["accession"] != query_acc]
        if congen:
            congen.sort(key=lambda r: -{"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}.get(r["level"], 0))
            m = congen[0]
            return (normalize_query_asm(m["accession"]),
                    str(local_path_for(m["accession"], m["asm_name"])),
                    f"congeneric_catalog:{m['accession']}")
        return ("", "", "no_ref_found")

    for u in user_panel:
        sp = species_from_organism(u["organism"])
        meta = catalog.get(u["accession"])
        if not meta:
            sys.stderr.write(f"[augment] WARNING: {u['accession']} not in catalog; skipping\n")
            continue
        ref_id, ref_fa, note = pick_benchmark_ref(sp, u["accession"])
        # Pull taxonomy fields from any prepare row for this species if available.
        prep = prep_by_species.get(sp, {})
        row = {fn: "" for fn in fieldnames}
        row.update({
            "query_asm": normalize_query_asm(u["accession"]),
            "query_mode": "assembly",
            "path": str(local_path_for(u["accession"], meta["asm_name"])),
            "scenario": "million_real",
            "lifestyle": prep.get("lifestyle", ".") or ".",
            "architecture": prep.get("architecture", ".") or ".",
            "benchmark_ref_asm": ref_id,
            "benchmark_ref_fasta": ref_fa,
            "phylum": prep.get("phylum", "") or ".",
            "class": prep.get("class", "") or ".",
            "order": prep.get("order", "") or ".",
            "family": prep.get("family", "") or ".",
            "genus": prep.get("genus", "") or sp.split()[0],
            "species": prep.get("species", "") or sp,
            "source": "ncbi-best",
            "instrument_platform": ".",
            "library_layout": ".",
            "run_accession": ".",
        })
        out_rows.append(row)
        audit_rows.append({
            "query_acc": u["accession"],
            "species": sp,
            "guild": u.get("guild", "."),
            "benchmark_ref_asm": ref_id,
            "ref_decision": note,
        })

    # Sanity: ensure FASTA paths actually exist on disk; flag missing.
    missing = [r["path"] for r in out_rows if not Path(r["path"]).exists()]
    if missing:
        sys.stderr.write(f"[augment] {len(missing)} query FASTA(s) not yet in {REFS_CACHE}; they'll be downloaded by the bootstrap reads phase or first benchmark task.\n")
        for m in missing[:10]:
            sys.stderr.write(f"  missing: {m}\n")
    missing_ref = [r["benchmark_ref_fasta"] for r in out_rows if r["benchmark_ref_fasta"] and not Path(r["benchmark_ref_fasta"]).exists()]
    if missing_ref:
        sys.stderr.write(f"[augment] {len(missing_ref)} benchmark-ref FASTA(s) not on disk yet.\n")

    # Backup, then overwrite query_manifest.tsv and benchmark_reference_map.tsv.
    backup_qm = prep_qm.with_suffix(".tsv.bak.preaugment")
    if not backup_qm.exists():
        backup_qm.write_bytes(prep_qm.read_bytes())
    backup_brm = prep_brm.with_suffix(".tsv.bak.preaugment")
    if prep_brm.exists() and not backup_brm.exists():
        backup_brm.write_bytes(prep_brm.read_bytes())

    with prep_qm.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    brm_fields = ["query_asm", "benchmark_ref_asm", "benchmark_ref_fasta", "species"]
    with prep_brm.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=brm_fields, delimiter="\t")
        writer.writeheader()
        for row in out_rows:
            writer.writerow({
                "query_asm": row["query_asm"],
                "benchmark_ref_asm": row["benchmark_ref_asm"],
                "benchmark_ref_fasta": row["benchmark_ref_fasta"],
                "species": row["species"],
            })

    audit_path = args.prepared_dir / "injected_panel_audit.tsv"
    with audit_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["query_acc", "species", "guild", "benchmark_ref_asm", "ref_decision"], delimiter="\t")
        writer.writeheader()
        writer.writerows(audit_rows)

    sys.stderr.write(f"[augment] wrote {len(out_rows)} queries to {prep_qm}\n")
    sys.stderr.write(f"[augment] wrote {len(out_rows)} ref rows to {prep_brm}\n")
    sys.stderr.write(f"[augment] audit: {audit_path}\n")
    # Last line of stdout = the count, parseable by bash wrappers.
    print(len(out_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
