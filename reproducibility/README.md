# Reproducibility

Exact inputs, commands, tool versions and checksums for reproducing the MycoSV
main figures, tables and benchmarks. Addresses the reviewers' code-availability
requests (accession lists, exact commands, expected outputs, checksums).

## 1. Accession lists (source data)
- **Five-genome discovery panel + long-read-validated genomes**: accessions in
  the manuscript (Table 1, Table 2, Table 2b).
- **200-genome panel**: `manuscript/supplementary/Supplementary_Data_1_200_query_accessions.tsv`
  (query accession, NCBI-derived ranks, benchmark reference, taxonomy source).
- **24,217 indexed assemblies** (scalability): `manuscript/supplementary/Supplementary_Data_2_24217_indexed_accessions.tsv`.
- Taxonomy source for all: NCBI RefSeq/GenBank `assembly_summary` + NCBI Taxonomy
  (assembly taxid); no sequence-derived or algorithmic taxonomy was used.

## 2. Core caller
MycoSV is a single self-contained C++ binary (no third-party tool in the calling
path). Build from the repository root:
```
g++ -O2 -std=c++17 -pthread main.cpp -o fungi_graphsv_tol_bin
```
Run (assembly mode) against a prepared clade index/registry:
```
./fungi_graphsv_tol_bin --tol-hierarchical --tol-index-dir <index> \
  --tol-registry-dir <registry> --ref-list <refs.txt> --query-list <queries.txt> \
  --out-prefix <out>/calls --query-mode assembly --threads 8 --graph-native-mode
```

## 3. Element-annotation benchmark (Reviewer comment 8)
Third-party tools are used only as **comparators**, never in MycoSV's calls.
See `experiments/biology_benchmark/README.md` for the full workflow. Scripts are
mirrored in `scripts/`:
- `reannotate_elements.cpp` — re-classify element_class with the fixed classifier
  (structural-first; Pezizomycotina Starship gate; filamentous-only RIP gate;
  adaptive HGT threshold). Build: `g++ -O2 -std=c++17 -pthread -I<repo_root> reannotate_elements.cpp -o reannotate_elements`.
- `score_biology_comparators.py` → `biology_comparator_pr.tsv` (RIP/TE/Starship
  precision-recall vs RIPCAL, OcculterCut, TEsorter, MycoMobilome, starfish, starbase).
- `score_ltr_benchmark.py` → `ltr_benchmark_summary.tsv` (SV-discovery enrichment
  in LTR regions vs de-novo LTRharvest and the MycoMobilome LTR subset).
- `reannotate_panel200.py` — 200-genome element reannotation.

## 4. Population analysis (Reviewer comment 3)
`scripts/build_population_sfs.py <by_query_dir> <out_prefix>` → occupancy matrix +
site-frequency spectrum (Fig. 4). Figure: `manuscript/build_figure4_population.py`.

## 5. Main figures
Run figure builders from the **repository root**:
```
python3 manuscript/build_figure1_assembly.py    # Figure 2 (five-genome panel)
python3 manuscript/build_figure2_panel.py        # Figure 3 (200-genome panel)
python3 manuscript/build_figure4_population.py    # Figure 4 (population SFS)
```

## 6. Versions and checksums
- `tool_versions.txt` — comparator and core-tool versions.
- `checksums.md5` — md5 of the key result tables feeding the manuscript.
  Verify with `md5sum -c reproducibility/checksums.md5` from the repository root.

## 7. Container (scaffold)
`Dockerfile` builds the MycoSV caller and the Python analysis environment. The
comparator tools (RIPCAL, TEsorter, starfish, OcculterCut, LTRharvest) and the
MycoMobilome/starbase databases are large external resources fetched separately
(see `experiments/biology_benchmark/README.md` and the Zenodo DOIs above).
