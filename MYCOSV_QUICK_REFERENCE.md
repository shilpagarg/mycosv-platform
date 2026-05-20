# MycoSV Quick Reference

## Algorithm Overview

**Three-Layer Hierarchical SV Caller** for fungal pangenomes:

```
Query (Assembly/Reads)
    ↓
[Query Input Handler] — Auto-detect mode; convert reads → consensus
    ├─ Assembly: use directly
    ├─ Long-reads: k-mer consensus clustering (k=12)
    └─ Short-reads: de-Bruijn unitigs (k=21)
    ↓
[Layer 3] Phylum-Sharded Routing — VP-tree O(log N) nearest neighbors
    ↓
[Layer 2] Clade Registry — LRU cache O(1) per-clade graph loading
    ↓
[Layer 1] SV Calling — Syncmer seeding, chain-refine, bubble classification
    ↓
VCF + TSV output (5 SV types: INS/DEL/DUP/INV/TRA)
```

---

## Key Algorithms

### 1. Syncmer Seeding (Layer 1)
- **Definition**: Order-preserving k-mers (k=21, s=11).
- **Time**: O(N) linear scan.
- **Purpose**: Exact-match anchors between query and reference.

```python
def syncmers(seq, k=21, s=11):
    """Generate syncmers from sequence."""
    anchors = []
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        prefix = kmer[:s]
        suffix = kmer[k-s:k]
        if prefix < suffix and prefix < kmer[s:]:
            anchors.append((i, hash(kmer)))
    return anchors
```

### 2. Chain-and-Refine Alignment (Layer 1)
- **Chaining**: Connect syncmers into blocks (score = anchor_count - gap_penalty).
- **Refinement**: Local Smith-Waterman alignment (~1000 bp window).
- **Time**: O(S log S) chaining + O(W²) refinement per block.

### 3. VP-Tree Routing (Layer 3)
- **Build**: O(N log N) hierarchical partitioning.
- **Query**: O(log N) nearest neighbors with distance pruning.
- **Purpose**: Route query to candidate clades without scanning all references.

```
Query Sketch → VP-Tree Root
              ↙         ↘
         Closer        Farther
        Clades         Clades
        (prune by distance)
              ↓
    Candidate Clades (top-N)
```

### 4. Bubble Detection & Classification (Layer 1)
- **Bubble**: Pair of vertex-disjoint paths in pangenome graph.
- **Types**:
  - **Trivial**: Single edge vs edge sequence (INS/DEL).
  - **Triallelic**: 3+ paths (complex rearrangements).
  - **Non-ref**: Paths absent from reference (neo-alleles).

**Classification**:
```
Bubble → Is collinear? → YES → Multi-allele (polymorph)
           ↓ NO
         Has inversion? → YES → Inversion
           ↓ NO
         Spans contigs? → YES → Translocation
           ↓ NO
         Unclassified (complex)
```

### 5. LRU Cache (Layer 2)
- **O(1) insert/evict**: Hashtable + linked list.
- **Capacity**: 16 GB default (tunable).
- **Purpose**: Avoid re-loading per-clade pangenome graphs from disk.

```cpp
// Pseudocode
struct LRUCache {
    unordered_map<string, CacheEntry> cache;  // O(1) lookup
    list<string> lru_order;                    // Eviction order (LRU at tail)
    
    void insert(const string& clade, const Graph& graph) {
        bytes_used += size(graph);
        while (bytes_used > max_bytes) {
            evict();  // Remove tail
        }
        cache[clade] = {graph, now()};
        lru_order.push_front(clade);  // MRU at head
    }
};
```

---

## SV Types (5 Classes)

| Type | Definition | Detection Method |
|------|-----------|------------------|
| **INS** | Insertion | alt_path longer; novel sequence |
| **DEL** | Deletion | ref_path longer |
| **DUP** | Duplication | Collinear multi-copy |
| **INV** | Inversion | alt_path = reverse-complement(ref) |
| **TRA** | Translocation | Spans contigs or reverses order |

---

## Query Input Modes

| Mode | Input Format | Processing | Output |
|------|--------------|-----------|--------|
| **Assembly** | FASTA contigs | Use directly | Query directly |
| **Long-reads** | FASTA/FASTQ (ONT/PacBio) | k-mer consensus clusters (k=12) | Pseudo-contigs |
| **Short-reads** | FASTA/FASTQ (Illumina) | de-Bruijn unitigs (k=21) | Pseudo-contigs |

### Mode-Specific Auto-Tuning

| Parameter | Assembly | Long-Reads | Short-Reads |
|-----------|----------|-----------|------------|
| k | 21 | 15 | 21 |
| minAnchors | 2 | 1 | 3 |
| minBlockScore | 6.0 | 3.0 | 6.0 |
| chainGapBand | 5000 | 15000 | 5000 |
| Use secondary seeds | No | Yes | Yes |

---

## Repeat & TE Annotation (8 Detectors)

| Detector | Rule | Example |
|----------|------|---------|
| **Tandem Repeat** | Period 2–12 bp, ≥5 copies, ≥50 bp | (ATTGC)×10 |
| **LTR** | Direct terminals ≥50 bp + high-GC core | LTR (Gypsy, Copia) |
| **TIR** | Inverted terminals ≥30 bp | DNA transposons |
| **LINE/Helitron** | AT-rich + poly-A tails | Autonomous elements |
| **SINE** | Short (~300 bp) + high-GC + terminal | Derived from RNA |
| **STARSHIP** | AT-rich hull + genic cargo ≥1 kb | Ascomycete-specific |
| **HGT** | GC deviation > ±0.08 / ≥500 bp | Horizontal gene transfer |
| **RIP** | C/G ratio > 2.5 in 500 bp window | Repeat-induced mutation |

---

## Precision & Recall

### Wilson 95% Confidence Intervals

$$\text{Precision} = \frac{TP}{TP + FP} \quad \text{Recall} = \frac{TP}{TP + FN}$$

**Where**:
- TP = call matches truth within 10 bp (endpoint tolerance)
- FP = call with no matching truth
- FN = truth with no matching call

### Expected Performance

| Mode | Precision | Recall | Notes |
|------|-----------|--------|-------|
| Assembly | >97% | >97% | Reference-quality inputs |
| Long-reads | >99% | 90–95% | High TP detection; recall limited by coverage |
| Short-reads | 85–90% | 80–88% | Unitigs lost in repetitive regions |

### Mode-Specific Challenges

**Assembly**:
- ✗ Complex nested SVs (alignment ambiguity)
- ✓ Triallelic topology classification (DS-11) mitigation

**Long-reads**:
- ✗ Chimeric reads (spurious SVs)
- ✓ `lrMaxReadLen=300kb` guard
- ✗ Low coverage <10× (poor clustering)
- ✓ Auto-tune `lrMinCluster=1` at low coverage

**Short-reads**:
- ✗ de-Bruijn assembly loses unitigs (repetitive regions)
- ✗ Large SVs >5 kb (fewer spanning reads)
- ✓ `srMinUnitigLen=200` tunable

---

## Data Structures (DS-1 through DS-19)

| DS | Name | Time | Space | Purpose |
|----|------|------|-------|---------|
| 1 | O(1) LRU Cache | insert: O(1), evict: O(1) | O(G·C) | Per-clade graph caching |
| 4 | VP-Tree Routing | query: O(log N) | O(N) | Phylum-sharded routing |
| 5 | Bloom Filter Prefilter | test: O(k) | 64 KB | Quick-reject before distance |
| 7 | PathPositionIndex | insert: O(log N) | O(M) | TraIntra detection |
| 10 | Sparse-Table LCA | query: O(1) | O(N log N) | Lowest common ancestor |
| 11 | Triallelic Topology Classifier | classify: O(1) | O(1) | TRA_INTRA vs TRA_INTER |
| 13 | Suffix Array + LCP | build: O(N log N), query: O(log N) | O(N log Σ) | Maximal exact matches |
| 15 | VEB Tree | op: O(log log U) | O(U) | Predecessor/successor |
| 16 | Merge-Sort Tree | query: O(log² N) | O(N log N) | Interval stabbing |
| 17 | Fenwick Tree (BIT) | op: O(log N) | O(N) | Prefix-sum queries |
| 18 | Chain Treap | build: O(N log N), query: O(1) avg | O(N) | Seed chaining |

---

## Complexity Analysis

### Time Complexity per Query

```
Layer 3 Routing:     O(log N)       [VP-tree query + Bloom prefilter]
Layer 2 Cache:       O(1)           [Hashtable LRU lookup]
Layer 1 Seeding:     O(L)           [Linear scan, L=query length]
Layer 1 Chaining:    O(S log S)     [S=syncmer count; band constraint]
Layer 1 Refinement:  O(K·W²)        [K=block count, W=window size]
─────────────────────────────────────
Total:              O(log N + L + S log S + K·W²)
```

### Space Complexity

```
Layer 3 Index:      O(N)           [VP-tree + sketches + Bloom filters]
Layer 2 Cache:      O(min(G·C, cache_bytes))  [LRU-bounded]
Layer 1 Working:    O(L + B)       [Query + blocks]
─────────────────────────────────────
Total:              O(N + cache + L)
```

### Scaling

- **Queries**: Linear in query count (batch-processed).
- **References**: Logarithmic (VP-tree O(log N) routing).
- **Clades accessed**: O(1)–O(10) (depends on query specificity).
- **Cache**: Bounded by max_cache_gb (LRU enforcement).

---

## Key Parameters

### Essential

```bash
--ref-list refs.txt              # One assembly per line
--query-list queries.txt         # One assembly/reads per line
--out-prefix results/callset     # Output prefix
--tol-index-dir indexes/tol      # Layer 3 index directory
```

### Query Mode

```bash
--query-mode assembly       # Pre-assembled (default)
--query-mode long-reads     # ONT/PacBio (auto-detectable)
--query-mode short-reads    # Illumina (auto-detectable)
```

### Syncmer Seeding

```bash
--k 21              # Primary k-mer length
--s 11              # Smer window (must be < k)
--seed-stride 1     # Every n-th seed (1 = all)
```

### Calling Thresholds

```bash
--min-sv-len 40             # Minimum SV length (bp)
--max-sv-len 1000000        # Maximum SV length (bp)
--min-block-score 6.0       # Minimum log-odds score
--min-anchors 2             # Minimum syncmer anchors per call
```

### Caching (Layer 2)

```bash
--tol-cache-gb 16          # Max cache size (GB)
--tol-cache-entries 128    # Max clades in memory
```

### Parallelism

```bash
--threads 16                # CPU cores to use
```

---

## Performance Benchmarks

### Simulated Data
- **Test suite**: 13 ecological scenarios (23 domain corrections).
- **Query modes**: Assembly, long-reads, short-reads.
- **SV types**: All 5 (INS/DEL/DUP/INV/TRA).
- **Results**: 24/24 PASS (100%).

### Real Fungal Data
- **Panels**: 5 (compact yeast, AMF, cross-phylum HGT, TE-rich, two-speed).
- **Integration**: NCBI RefSeq + ENA samples + minigraph truth.
- **Validation**: Precision/recall vs truth VCF.
- **Results**: 13/13 PASS (100%).

### Timing (1M references)

```
Single query vs 1M refs (routing):  ~2–5 sec
Per-clade SV calling:               ~1–10 sec (per clade)
Total (10 clades avg):              ~5–15 sec per query
1000 queries:                        ~4 hours (16-core)
```

### Memory

```
Layer 1 per-clade:                  50–500 MB (compressed)
Layer 2 cache (16 GB):              ~30–60 clades
Layer 3 index (1M refs):            ~1–2 GB
```

---

## Common Issues

| Issue | Cause | Solution |
|-------|-------|----------|
| Low recall (short-reads) | Unitigs lost in repetitive regions | Increase `--sr-min-unitig-len` |
| High FP (long-reads) | Chimeric reads create spurious calls | Decrease `--lr-max-read-len` |
| Slow query | Many candidate clades | Check VP-tree balance; rebuild if needed |
| Out of memory | Cache too large | Reduce `--tol-cache-gb` |
| Missing assembly | Query mode mis-detected | Explicitly set `--query-mode assembly` |

---

## Output Files

For `--out-prefix PREFIX`, one unified call set is written:

```
PREFIX.hits.tsv         # All SV calls (one per row); alignment_mode column tags origin
PREFIX.vcf              # Same calls in VCF 4.3 (CLADE, CLADE_RANK, OFFREF, OFF_REF_TIER)
PREFIX.gfa              # Pangenome graph fragment (graph-native mode only)
PREFIX.ancestral.tsv    # ToL ancestral-state annotations (optional)
PREFIX.te_predictions.tsv  # TE classifier output
```

**Pangenome- vs single-reference calls share the same files.** Split via the `alignment_mode` column:

| Origin | `alignment_mode` values |
|---|---|
| Pangenome (multi-ref / graph-native) | `mem_chain_cached_single_ref_multi`, `..._multi;secondary_seed_rescue`, `graph_native_offref_window` |
| Single-reference (pairwise) | `mem_chain_cached_single_ref`, `...;secondary_seed_rescue`, `simple_length_fallback`, `simple_offref_fallback`, `reads_mode_kmer_fallback` |

Pangenome rows are additionally identifiable by a non-species `CLADE_RANK` (phylum/class/order/family/genus) emitted by `hierarchical_call_assembly_multirank`.

---

## TE Classification Quick Reference

### Method
k-mer nearest-centroid via VPTree (same infrastructure as Layer 3 routing).
- **Sketch**: canonical k-mer FracMin (k=21, p=0.05 → ~5% of k-mers kept)
- **Distance**: Jaccard (1 − |A∩B|/|A∪B|)
- **Taxonomy**: class / order / superfamily (PanTEon/RepBase label format)

### Label format
```
>ID#Class/Order/Superfamily     e.g.  >TE001#LTR/Gypsy/Chromovirus
>ID#Class/Superfamily           e.g.  >TE002#LTR/Copia
```

### Commands
```bash
# Build index from labeled training FASTA
echo train.fasta > train.lst
./fungi_graphsv_tol --te-train --query-list train.lst \
    --te-index-prefix models/te_clf

# Classify unknown sequences
echo unknown.fasta > test.lst
./fungi_graphsv_tol --te-classify --query-list test.lst \
    --te-index-prefix models/te_clf --out-prefix results/te
# → results/te.te_predictions.tsv

# Benchmark vs PanTEon SOTA tools
python3 run_te_benchmark.py \
    --train-fasta repbase_fungi_train.fa \
    --test-fasta  repbase_fungi_test.fa \
    --out-dir     te_benchmark/

# Demo (no data download needed)
python3 run_te_benchmark.py --download-fungi-demo --out-dir te_demo/

# Install SOTA TE tools (NeuralTE, DeepTE, TERL, Terrier, ClassifyTE, CREATE, TEClass2)
bash install_tools.sh --te-only
```

### Parameters
| Flag | Default | Notes |
|------|---------|-------|
| `--te-k` | 21 | k-mer length |
| `--te-fracmin-p` | 0.05 | Sketch density (0–1) |
| `--te-max-hashes` | 4096 | Max hashes per centroid |
| `--te-index-prefix` | — | Path stem for .vptree/.meta files |

### PanTEon benchmark reference (fungi, NeuralTE best)
| Level | F1 (best paper) |
|-------|----------------|
| Class | 0.88 |
| Order | 0.79 |
| Superfamily | 0.72 |

---

## Tool Installation

Install all SV callers and TE classifiers into the `mycosv` conda environment:

```bash
# Full install (SV tools + TE tools + build MycoSV binary)
bash install_tools.sh

# Check what is already available
bash install_tools.sh --check

# SV tools only (SyRI, minigraph, PGGB, Delly, Manta, SVIM, Sniffles, cuteSV)
bash install_tools.sh --sv-only

# TE tools only (NeuralTE, DeepTE, TERL, Terrier, ClassifyTE, CREATE, TEClass2)
bash install_tools.sh --te-only

# Build MycoSV binary only (assumes conda env already active)
bash install_tools.sh --mycosv-only
```

---

## References

- Hong, C., et al. (2016). "Optimizing seed size yields improved sensitivity in read mapping"
- Yianilos, P. (1993). "Data structures and algorithms for nearest neighbor search in general metric spaces"
- Slot, J. & Rokas, A. (2011). "Horizontal transfer of a large and highly toxic secondary metabolite gene cluster between fungi"
- Urquhart, A., et al. (2023). "Giant transposons with structured cargo of metabolic genes"

---

**MycoSV v1.0** | Last Updated: 20 April 2026
