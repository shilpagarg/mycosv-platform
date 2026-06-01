# MycoSV-platform

MycoSV-platform is a hierarchical, taxonomy-aware pangenome framework for
structural-variant discovery in fungal genomes. It routes assemblies or reads
through clade-aware pangenome graph layers to recover canonical variants and
off-reference sequence that single-reference pipelines can miss.

## Key Features

- Three-layer Tree-of-Life architecture: phylum-sharded routing, clade graph
  registry/cache, and graph-native SV calling.
- Supports assembly, long-read, short-read, and auto-detected query modes.
- Detects insertions, deletions, duplications, inversions, translocations, and
  off-reference novel sequence.
- Annotates repeat, TE, RIP, HGT, and Starship-like candidate biology.
- Includes fungal benchmark, million-scale simulation, read-validation, and
  manuscript-figure/report utilities.

## Repository Contents

| Path | Purpose |
| --- | --- |
| `main.cpp` | Main C++ CLI and orchestration entry point |
| `layer1_clade_graph.hpp` | Per-clade pangenome graph and SV calling logic |
| `layer2_registry.hpp` | Clade graph registry and LRU cache |
| `layer3_routing_index.hpp` | Hierarchical routing index |
| `query_input_handler.hpp` | Assembly/read mode preprocessing |
| `fungi_tol_bridge.hpp` | Tree-of-Life integration layer |
| `run_real_fungal_benchmark.py` | Real fungal panel preparation and benchmark runs |
| `run_million_mode_query_benchmark.py` | Million-catalog scaling benchmark |
| `plot_mycosv_pangenome_calls.py` | Pangenome-call summaries, tables, and plots |
| `sv_visualization_report.py` | HTML/PNG reporting utilities |
| `manuscript/` | Manuscript tables and figure-generation scripts |

## Requirements

- Linux or HPC environment
- `g++` with C++17 support
- Python 3.9+
- `pytest` for tests
- Optional external genomics tools for full benchmark workflows; see
  `install_tools.sh`

## Build

From the repository root:

```bash
g++ -O2 -std=c++17 -pthread -I. main.cpp -o fungi_graphsv_tol_bin
```

To inspect available CLI options:

```bash
./fungi_graphsv_tol_bin --help
```

## Minimal Usage

Run MycoSV on assembly FASTA query lists:

```bash
./fungi_graphsv_tol_bin \
  --ref-list refs.txt \
  --query-list queries.txt \
  --query-mode assembly \
  --out-prefix results/mycosv
```

Run with long reads or short reads:

```bash
./fungi_graphsv_tol_bin \
  --ref-list refs.txt \
  --query-list long_reads.txt \
  --query-mode long-reads \
  --out-prefix results/mycosv_long_reads

./fungi_graphsv_tol_bin \
  --ref-list refs.txt \
  --query-list short_reads.txt \
  --query-mode short-reads \
  --out-prefix results/mycosv_short_reads
```

Build a hierarchical Tree-of-Life index from a manifest:

```bash
./fungi_graphsv_tol_bin \
  --tol-hierarchical \
  --tol-build-index manifest.tsv \
  --tol-index-dir tol_index \
  --tol-registry-dir tol_registry \
  --tol-index-threads 16
```

Run using the hierarchical index:

```bash
./fungi_graphsv_tol_bin \
  --tol-hierarchical \
  --tol-index-dir tol_index \
  --tol-registry-dir tol_registry \
  --query-list queries.txt \
  --query-mode auto \
  --out-prefix results/mycosv_tol
```

Primary outputs are written with the chosen output prefix, including VCF and TSV
files such as `PREFIX.vcf` and `PREFIX.hits.tsv`.

## Benchmarks And Tests

Run the quick smoke tests:

```bash
pytest -q test_golden_smoke.py test_sv_report_smoke.py
```

Run the main small validation suite:

```bash
python3 -m pytest test_pipeline_features.py test_amf.py test_all_use_cases.py -v
python3 -m pytest test_real_fungal_benchmark.py test_new_biology_candidates.py -v
```

Run the project experiment wrapper:

```bash
bash run_all_experiments.sh --small
bash run_all_experiments.sh --large
bash run_all_experiments.sh --real
```

For full workflow details, see `EXPERIMENTS_GUIDE.md` and
`QUICK_COMMAND_REFERENCE.md`.

## Documentation

- `MYCOSV_ALGORITHM.md`: full algorithm and implementation specification
- `MYCOSV_QUICK_REFERENCE.md`: concise architecture, parameters, and outputs
- `DOCUMENTATION_INDEX.md`: map of project scripts, tests, and docs
- `SV_TYPE_COVERAGE.md`: structural-variant type coverage notes
- `BENCHMARKING_STRATEGY.md`: benchmark design and interpretation

## Citation And Credit

If you use MycoSV-platform, please credit this repository:

```text
MycoSV-platform: hierarchical pangenome structural-variant discovery for fungal genomes.
https://github.com/shilpagarg/mycosv-platform
```

Citation metadata is also provided in `CITATION.cff`.

## License

MycoSV-platform is released under the MIT License. See `LICENSE`.

## Contact

For issues, questions, feature requests, or benchmark contributions, email:

```text
shilpa.garg2k7@gmail.com
```
