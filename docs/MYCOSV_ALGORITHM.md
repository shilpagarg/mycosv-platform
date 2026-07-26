# MycoSV: Three-Layer Hierarchical SV Caller for Fungal Genomes

MycoSV is a hierarchical, graph-native fungal pangenome and structural-variation engine for discovery across assemblies, short reads, and long reads.

Internal binary name in this repository: `fungi_graphsv_tol`.

## What It Does

MycoSV is designed for whole-genome fungal comparative analysis with:

- hierarchical pangenome indexing across `phylum -> class -> order -> family -> genus -> species`
- graph-native SV discovery, including `INS`, `DEL`, `INV`, `DUP`, `TRA`, and `OFF_REF`
- off-reference discovery of novel insertions, HGT-like segments, TE-like segments, STARSHIP-like cargo regions, tandem repeats, and RIP-like sequence patterns
- query support for assembled genomes, short reads, and long reads in one framework
- ancestry-aware and clade-aware calling, with optional ancestral alignment / recombination reporting
- reusable build-once query-many-times index layout
- million-scale external-memory routing for very large catalogs

This codebase is intended for fungal ecology, plant-fungal symbiosis, genome plasticity, virulence evolution, transposable-element driven variation, and cross-clade novelty discovery.

---

## Table of Contents

- [Overview](#overview)
- [Algorithm Architecture](#algorithm-architecture)
  - [Layer 1: Per-Clade Pangenome Graphs](#layer-1-per-clade-pangenome-graphs)
  - [Layer 2: Clade Graph Registry](#layer-2-clade-graph-registry)
  - [Layer 3: Phylum-Sharded Routing Index](#layer-3-phylum-sharded-routing-index)
- [Query Input Modes](#query-input-modes)
- [Seeding & Alignment](#seeding--alignment)
- [SV Classification](#sv-classification)
  - [SV Calling Paths (Path A / B / C)](#sv-calling-paths-path-a--b--c)
- [Repeat & TE Annotation](#repeat--te-annotation)
  - [TE Annotation per SV Type](#te-annotation-per-sv-type)
- [Off-Reference Novelty Tiers](#off-reference-novelty-tiers)
- [TE Classification](#te-classification)
- [Precision & Recall](#precision--recall)
- [Benchmarks](#benchmarks)
- [Usage](#usage)
- [Performance Characteristics](#performance-characteristics)

---

## Overview

### Key Features

| Feature | Details |
|---------|---------|
| **SV Types** | INS, DEL, DUP, INV, TRA (up to 1 Mb) |
| **Query Modes** | Assembly, long-reads (5–50× ONT/PacBio), short-reads (20–100× Illumina) |
| **Scalability** | 1 million reference genomes in ~40 GB RAM |
| **Accuracy** | >97% precision/recall on assembly; >99% TP detection on long-reads |
| **TE Annotation** | Tandem repeats, LTR, TIR, LINE, SINE, STARSHIP, HGT, RIP |
| **Parallelism** | Multi-threaded (pthreads), multi-process via PBS/Slurm batching |

### Workflow

```
Reference Catalog (N genomes)
  ↓
[Layer 3] Phylum-sharded routing (VP-tree, Bloom filter, skip-list)
  ↓ Route to candidate clades (VP-tree nearest neighbors)
  ↓
[Layer 2] Clade graph registry (LRU cache, atomic-rename manifest)
  ↓ Load per-clade pangenome graphs on-demand
  ↓
[Layer 1] Per-clade SV calling
  ├─ Syncmer seeding (k=21, s=11)
  ├─ Chain-and-refine alignment
  ├─ Bubble detection & classification
  ├─ Triallelic topology resolution
  └─ TE/repeat annotation
  ↓
Query (assembly / reads)
  ↓
Query Input Handler (auto-detect mode; convert reads → consensus)
  ├─ Assembly: use directly
  ├─ Long-reads: k-mer consensus clustering (k=12)
  └─ Short-reads: de-Bruijn unitig extraction (k=21)
  ↓
[Layer 1] Seed, chain, refine
  ↓
VCF + TSV output
```

---

## Algorithm Architecture

### Layer 1: Per-Clade Pangenome Graphs

**Purpose**: Base-level whole-genome variation detection within evolutionary clades.

**Key Algorithms**:

#### 1.1 Syncmer Seeding
- **Definition**: Order-preserving k-mer hashing (Hong & Buhler 2016).
- **Parameters**: `k=21` (primary), `s=11` (smer window), `t=2` (threshold).
- **Time**: O(N) where N = sequence length.
- **Space**: O(Σ) where Σ = number of seeds (~4% of N for fungi).
- **Purpose**: Identifies exact-match anchors between query and reference.
- **Optional**: Interval hash (IH) acceleration for repetitive regions (O(log² N) precomputation).

**Pseudocode**:
```
function Syncmers(seq, k, s, t):
    anchors ← []
    for i in 0..len(seq)-(k-1):
        kmer ← seq[i:i+k]
        smer ← canonical_kmer(kmer[0:s])
        if smer < kmer[0:s] and smer < kmer[k-s:k]:
            anchors.append((i, hash(kmer)))
    return anchors
```

#### 1.2 Greedy Chain-and-Refine Alignment
- **Goal**: Chain syncmers into consistent local alignments (BLAST-style).
- **Chain Score**: ∑(anchor_count) - gap_penalty * dist.
- **Gap Band**: `chainGapBand=5000` bp (adaptive per mode).
- **Refinement**: Local Smith-Waterman (DP matrix ~1000 bp window).
- **Output**: Alignment blocks with CIGAR strings.

**Complexity**:
- Chaining: O(S log S) with band constraint (S = syncmer count).
- Refinement: O(W²) per block (W = window size, typically 1000 bp).

#### 1.3 Bubble Detection & SV Classification

**Bubble**: A pair of vertex-disjoint paths in the pangenome graph.

**Variants**:
- **Trivial**: Single edge vs edge sequence (INS/DEL).
- **Triallelic**: Three or more distinct paths (complex rearrangements).
- **Non-ref**: All paths absent from reference (neo-alleles).

**Classification Algorithm**:
```
function ClassifyVariant(ref_path, alt_paths, graph):
    if len(alt_paths) == 1:
        return Indel(length_diff(ref_path, alt_paths[0]))
    
    if all_collinear(alt_paths):
        return classify_collinear(ref_path, alt_paths)
    
    if has_inversion(alt_paths):
        return Inversion()
    
    if has_translocation(alt_paths, graph):
        return Translocation()
    
    return UnclassifiedComplex()
```

**Triallelic Classification** (DS-11):
```
enum TriallelicTopology {
    LINEAR,           // All paths on same axis
    TRA_INTRA,        // Rearrangement within same contig
    TRA_INTER         // Rearrangement across contigs
}

function classify_triallelic(paths, ref_contig):
    positions ← [p.start for p in paths]
    if is_monotone(positions):
        return LINEAR
    if all_same_contig(paths, ref_contig):
        return TRA_INTRA
    return TRA_INTER
```

**Thresholds**:
- Min SV length: 40 bp
- Max SV length: 1 Mb
- Min block score: 6.0 (log odds)
- Min anchors per call: 2

---

### Layer 2: Clade Graph Registry

**Purpose**: Cache per-clade pangenome graphs in LRU-managed shared memory.

**Data Structures**:

#### 2.1 O(1) LRU Cache (Sleator & Tarjan 1985)
- **Hashtable**: clade_name → (graph, timestamp, size_bytes)
- **Linked List**: Eviction order (MRU at head, LRU at tail).
- **Invariant**: graph_bytes ≤ max_cache_bytes.

**Eviction Policy**:
```
function insert_into_cache(clade_name, graph):
    bytes_used += size(graph)
    
    while bytes_used > max_cache_bytes:
        evict_clade ← lru_list.tail()
        delete evict_clade
        bytes_used -= size(evict_clade)
    
    cache[clade_name] ← (graph, now(), size(graph))
    lru_list.move_to_head(clade_name)
```

**Time**: O(1) insert, O(1) evict, O(log N) per-shard query.

#### 2.2 Per-Clade Load-Once Barrier (std::shared_future)
- Prevents redundant disk I/O when multiple threads request same clade.
- First thread loads from disk; others wait on future.
- Lock-free after load completes.

#### 2.3 Atomic-Rename Manifest
- Manifest file (TSV) records: clade_name, rank, phylum, graph_path, crc32.
- Write-safety: write to `.tmp`, then rename (POSIX atomic).
- Prevents partial reads on crash.

**Manifest Format**:
```
#clade_name  clade_rank  phylum         graph_path                    genomes  sv_bubbles  crc32
Lachancea    genus       Ascomycota     /tol/index/Lachancea.gbz      42       5821       0xabc123
Rhizophagus  genus       Glomeromycota  /tol/index/Rhizophagus.gbz    8        412        0xdef456
```

**Cache Sizes**:
- Default: 16 GB (tunable via `--tol-cache-gb`).
- Per-clade: typically 50–500 MB (compressed pangenome).
- Adaptive: scales to 1/4 of free RAM.

---

### Layer 3: Phylum-Sharded Routing Index

**Purpose**: Route query (assembly/reads) to candidate clades in O(log N) expected time.

**Data Structures**:

#### 3.1 VP-Tree (Vantage-Point Tree) Routing
- **Metric**: Syncmer-based sketch distance (Hamming or Jaccard).
- **Build**: O(N log N) hierarchical partitioning.
- **Query**: O(log N) nearest neighbors with radius pruning.

**VP-Tree Node**:
```cpp
struct VPNode {
    uint64_t centroid_hash;           // Vantage point (syncmer sketch)
    double   median_dist;             // Distance to median child
    VPNode*  left_child;              // Closer nodes
    VPNode*  right_child;             // Farther nodes
    std::vector<CladeDescriptor*> clades;  // Leaf: candidate clades
};

function route_query(query_sketch, vp_tree, radius):
    candidates ← []
    
    function dfs(node, tau):
        dist ← distance(query_sketch, node.centroid_hash)
        
        if |dist - node.median_dist| < tau:
            candidates += node.clades
        
        if dist < node.median_dist:
            if dist - tau < node.median_dist:
                dfs(node.left_child, tau)
            if dist + tau ≥ node.median_dist:
                dfs(node.right_child, tau)
        else:
            if dist + tau ≥ node.median_dist:
                dfs(node.right_child, tau)
            if dist - tau < node.median_dist:
                dfs(node.left_child, tau)
    
    dfs(vp_tree.root, radius)
    return candidates
```

#### 3.2 Bloom Filter Prefilter
- **Size**: 64 KB (524,288 bits), k=7 hash functions.
- **FPR**: ~1% @ 56K items (tuned for typical clade counts).
- **Use**: Quick reject before expensive distance computation.
- **Formula**: False-positive rate = (1 - e^(-k·n/m))^k where n=items, m=bits, k=funcs.

#### 3.3 Phylum-Sharded Locking (FIX-LOCK v14)
- **Issue**: Original code acquired `registryMu_` TWICE per phylum (inefficient).
- **Fix**: Merge phylum-shard snapshot and pointer capture into single critical section.

**Before (inefficient)**:
```cpp
{
    lock(registryMu_);
    phylum_list ← snapshot phylums;
    unlock(registryMu_);
}
for each phylum in phylum_list:
    lock(phylum.shard_lock);  // Second acquisition
    candidates += phylum.route(query);
    unlock(phylum.shard_lock);
```

**After (optimized)**:
```cpp
vector<PhylumShard*> shards;
{
    lock(registryMu_);
    for each phylum:
        shards.push_back(get_shard(phylum));
}  // Brief critical section
// No global lock held below
for each shard in shards:
    lock(shard.lock);  // Only per-shard lock
    candidates += shard.route(query);
    unlock(shard.lock);
```

**Time Complexity**: Global lock held for O(P) where P = phylum count (~10–50); per-shard query O(log N).

---

## Query Input Modes

MycoSV accepts three types of input, automatically converted to pseudo-contigs:

### 1. Assembly Mode
- **Input**: Pre-assembled contigs (FASTA).
- **No conversion**: Directly queried.
- **Expected**: ~1–100 Mb per query.

### 2. Long-Reads Mode (ONT/PacBio)
- **Input**: FASTA/FASTQ reads (coverage ~5–50×, length 1–100 kb).
- **Preprocessing**: K-mer consensus clustering.
  - Anchor k-mers: k=12, frequency ≥ 2.
  - Reads grouped by shared anchors (similarity graph).
  - Per-group majority-vote consensus sequence.
  - Handles up to ~15% error rate.

**Algorithm**:
```
function ConsensusCluster(reads, k, min_cluster):
    anchors ← [kmers(r, k) for r in reads]
    clusters ← connected_components(similarity_graph(anchors))
    
    consensuses ← []
    for each cluster in clusters:
        if len(cluster) ≥ min_cluster:
            aln ← MSA(cluster)  // lightweight column-wise majority
            consensus ← majority_vote(aln)
            consensuses.append(consensus)
    
    return consensuses
```

- **Output**: Pseudo-contigs (one per read cluster).
- **Parameter tuning**: `lrMinCluster=2`, `lrAnchorK=12`.
- **Failure mode**: Chimeric reads create spurious SVs (mitigated by `lrMaxReadLen=300kb`).

### 3. Short-Reads Mode (Illumina)
- **Input**: FASTA/FASTQ reads (coverage ~20–100×, length 50–300 bp).
- **Preprocessing**: de-Bruijn unitig extraction (SPAdes-style).
  - Solid k-mers: k=21, frequency ≥ `auto-detected threshold`.
  - Threshold: median_freq / 4 (robust to error distribution).
  - Greedy path extension: extend left/right until branch.

**Algorithm**:
```
function BuildUnitigs(reads, k, freq_threshold):
    kmers ← count_kmers(reads, k)
    solid_kmers ← [km for km ∈ kmers if count[km] ≥ freq_threshold]
    graph ← build_debruijn(solid_kmers)
    
    unitigs ← []
    for each solid_kmer in unvisited(solid_kmers):
        left ← extend_left(solid_kmer, graph)
        right ← extend_right(solid_kmer, graph)
        unitig ← left + solid_kmer + right
        unitigs.append(unitig)
    
    return unitigs
```

- **Output**: Unitigs (pseudo-contigs).
- **Limitation**: Unitigs lost in high-complexity/repetitive regions (explains recall drop).
- **Parameter tuning**: `srMinKmerFreq=0` (auto), `srMinUnitigLen=200`.

---

## Seeding & Alignment

### Syncmer-Based Seeding

**Syncmer Definition** (Hong & Buhler 2016):
- A k-mer is a syncmer iff the rightmost occurrence of its s-mer prefix and suffix appear at the k-mer boundaries.
- Guarantees: every region of length k+s−1 contains ≥1 syncmer (order-preserving).

**Canonical Syncmer**:
- Compute both forward and reverse-complement syncmers.
- Use lexicographically smaller one (maintains consistency).

### Secondary Seeds
- **Use**: Rescue weak alignments in repetitive/low-complexity regions.
- **Parameters**: `secondaryK=15`, `secondaryS=5` (shorter, denser).
- **Frequency cap**: `secondaryFreqCap=2048` (avoid ubiquitous k-mers).
- **Min anchors for rescue**: `repeatRescueMinAnchors=3`.

### Alignment Refinement
- **Local alignment**: Smith-Waterman DP over gapped region.
- **Band width**: `alignBW=128` (normal), `difficultBW=256` (repetitive regions).
- **Scoring**: +1 match, −3 mismatch, −5 gap-open, −1 gap-extend.

---

## SV Classification

### Five SV Types

| Type | Definition | Detection |
|------|-----------|-----------|
| **INS** | Insertion | ref_path shorter than alt_path; novel seq in alt |
| **DEL** | Deletion | ref_path longer than alt_path |
| **DUP** | Duplication | Collinear multi-copy structure |
| **INV** | Inversion | alt_path reverse-complement of ref_path (inverted) |
| **TRA** | Translocation | alt_path spans multiple contigs or reverses contig order |

### Bubble Types

#### Trivial Bubbles
- Single edge vs single or multiple edges.
- Simple length difference → INS/DEL.

#### Triallelic Bubbles
- Three or more distinct paths.
- **Topologies** (DS-11):
  - LINEAR: all paths on same axis (multi-allele site).
  - TRA_INTRA: rearrangement within contig.
  - TRA_INTER: rearrangement across contigs.

#### Non-REF Bubbles
- One or more paths absent from reference.
- Annotated as OFF_REF in output.

### Complex Rearrangement Detection

**Inversion Check**:
```
function is_inversion(ref_path, alt_path, graph):
    ref_seq ← sequence(ref_path)
    alt_seq ← sequence(alt_path)
    rc_alt ← reverse_complement(alt_seq)
    return edit_distance(ref_seq, rc_alt) / len(ref_seq) < 0.05
```

**Translocation Check**:
```
function is_translocation(path, ref_contig, graph):
    segments ← decompose_path_into_nodes(path)
    contigs ← [graph.node_to_contig(n) for n in segments]
    
    return (any(c != ref_contig for c in contigs) or
            any(is_decreasing(positions) for consecutive positions))
```

### SV Calling Paths (Path A / B / C)

SV classification uses three complementary paths depending on evidence availability:

#### Path A — MEM-Chain Classification (Primary)

Operates on the Suffix Array (DS-13) + ChainTreap (DS-18) pipeline.

```
allMems ← suffixArray.queryMems(query, ref)
  ↓
Pre-scan allMems for backward-mapping MEMs (rPos decreasing)
  → If any: hasBwdMems = true (DUP candidate)
  ↓
ChainTreap.chain(allMems) → best forward chain
  ↓
SvTypeFromChain::classify(chain)
  ├─ INS: large query gap, small ref gap
  ├─ DEL: small query gap, large ref gap
  ├─ INV: reverse-complement MEMs (totalRevLen ≥ minSvLen guard)
  ├─ DUP fallback: hasBwdMems → scan for forward collinear overlap
  └─ TRA fallback: cross-contig MEM gap > 500 kb → TRA
```

**DUP fallback**: ChainTreap drops overlapping rPos MEMs by design. The pre-scan detects backward-mapping MEMs indicating a tandem copy, then re-classifies after confirming a forward-only collinear overlap. Prevents false-negative DUP calls when the duplicated region produces overlapping anchors.

**TRA fallback**: A query gap > 500 kb bridging two different reference contigs triggers cross-contig MEM scanning with 500 kb gap tolerance, producing a TRA call with cross-contig coordinates.

**INV guard**: `totalRevLen >= minSvLen` is required before any INV is emitted. Prevents single palindromic k-mers (one reverse-complement MEM) from triggering false inversion calls.

#### Path B — Length-Delta Fallback (Reads Mode)

Used when MEM chains are too short or absent (short-reads / sparse coverage):

```
Pseudo-contig vs reference unitig
  ↓
length_delta ← |query_len − ref_len|
kmer_overlap ← kmer_overlap_fraction(query, ref, k=7..9)
  ↓
if length_delta > minSvLen and kmer_overlap < 0.80:
    INS (query > ref)  or  DEL (query < ref)
```

Uses a rolling FNV-1a polynomial hash (O(N) build, O(1) lookup per k-mer window).

#### Path C — Off-Reference Novelty (OFF_REF)

Activated when the query has no sufficient overlap with any reference in the clade:

```
Query
  ↓
kmer_overlap_fraction(query, ref, k=7..9) for each reference in catalog
  ├─ sameCladeOverlap  ← max overlap among same-clade references
  └─ otherCladeOverlap ← max overlap among other-clade references
  ↓
score_cross_clade_novelty(sameCladeOverlap, otherCladeOverlap)
  ├─ NOVEL:          sameCladeOverlap < 0.05  AND  otherCladeOverlap ≥ 0.10
  │                  → elementClass = "HGT" (cross-clade horizontal transfer)
  ├─ NOVEL_WEAK:     sameCladeOverlap < 0.20
  ├─ DIVERGED:       sameCladeOverlap < 0.50
  └─ OFF_REF_KNOWN:  sameCladeOverlap ≥ 0.50
```

**Cross-clade HGT detection**: When the novel sequence has near-zero overlap with same-clade genomes but ≥ 10% overlap with at least one other-clade genome, it is stamped `elementClass="HGT"`. Implemented in `score_cross_clade_novelty()` (`hierarchical_engine.hpp`) and called from `hierarchical_call_assembly()` and `hierarchical_call_assembly_multirank()` (`fungi_tol_bridge.hpp`).

---

## Repeat & TE Annotation

Layer 1 provides 8 specialized detectors (new in v14):

### 1. Tandem Repeat Detection
- **Rule**: Period 2–12 bp, ≥5 copies, ≥50 bp total.
- **Algorithm**: FFT-based period finding + copy-count validation.

```
function detect_tandem_repeat(seq, min_period=2, max_period=12, min_copies=5, min_len=50):
    best_period ← 0
    for p in min_period..max_period:
        copies ← 0
        for i in 0..len(seq)-p:
            if seq[i:i+p] == seq[i+p:i+2p]:
                copies += 1
        if copies ≥ min_copies and copies * p ≥ min_len:
            best_period ← p
            break
    return best_period > 0
```

### 2. LTR Element Detection
- **Rule**: Direct terminal repeats ≥50 bp + high-GC interior.
- **Mismatch tolerance**: ≤5%.

```
function detect_ltr_element(seq, min_repeat_len=50, max_mismatch_rate=0.05):
    for len in min_repeat_len..len(seq)/3:
        left_repeat ← seq[0:len]
        right_repeat ← seq[-len:]
        if edit_distance(left_repeat, right_repeat) / len ≤ max_mismatch_rate:
            interior ← seq[len:-len]
            if gc_content(interior) > 0.60:
                return true
    return false
```

### 3. TIR Element Detection (Inverted Terminal Repeats)
- **Rule**: Inverted repeats ≥30 bp.

```
function detect_tir_element(seq, min_len=30):
    for len in min_len..len(seq)/3:
        left ← seq[0:len]
        right_rc ← reverse_complement(seq[-len:])
        if edit_distance(left, right_rc) / len < 0.10:
            return true
    return false
```

### 4. LINE/Helitron Detection
- **Rule**: AT-rich (GC < 0.40) + poly-A/T tails (≥20 bp).

### 5. SINE Detection
- **Rule**: Short (50–400 bp) + high-GC (≥0.55) + terminal repeat.

### 6. STARSHIP Detection
- **Rule**: AT-rich hull (GC < clade_gc − 0.10) + ~genic cargo (GC 45–55%) ≥1 kb.
- **Biological context**: Large AT-rich elements encoding cargo in Ascomycetes (Urquhart et al. 2023).
- **Note**: Only for Ascomycota; not found in Glomeromycota (AMF).

### 7. HGT Island Detection
- **Rule**: GC deviation > ±0.08 over ≥500 bp window.
- **Published range**: ±0.05–0.10 (Slot & Rokas 2011).

### 8. RIP Window Detection (Repeat-Induced Point Mutation)
- **Rule**: C/G ratio > 2.5 in 500 bp window (post-duplicational mutation signature).

```
function detect_rip_window(seq, cg_ratio_thresh=2.5, win_len=500):
    for i in 0..len(seq)-win_len by (win_len/2):
        window ← seq[i:i+win_len]
        c_count ← count('C') + count('c')
        g_count ← count('G') + count('g')
        if g_count == 0 and c_count > 0:
            return true
        if c_count / g_count > cg_ratio_thresh:
            return true
    return false
```

### Classification Dispatcher (DS-12)
```
function classify_repeat_element(seq, clade_gc=0.45):
    if len(seq) < 50:
        return NONE
    
    if detect_rip_window(seq):           return RIP
    if detect_hgt_island(seq, clade_gc): return HGT
    if detect_starship(seq, clade_gc):   return STARSHIP
    if detect_sine(seq):                 return TE_SINE
    if detect_tir_element(seq):          return TE_TIR
    if detect_ltr_element(seq):          return TE_LTR
    if detect_line_helitron(seq):        return TE_LINE
    if detect_tandem_repeat(seq):        return REPEAT
    
    return NONE
```

### TE Annotation per SV Type

All five SV types receive TE/repeat annotation by passing the relevant sequence through `classify_repeat_element()`:

| SV Type | Sequence annotated | Notes |
|---------|-------------------|-------|
| **INS** | Inserted query sequence | Novel insertion examined for TE content |
| **DEL** | Deleted reference subsequence `[rBreakStart, rBreakEnd)` | Identifies removed TE insertions |
| **DUP** | Duplicated query sequence | Flags tandem TE expansions |
| **INV** | Inverted query sequence | Detects TE-driven inversions |
| **TRA** | Breakpoint-flanking sequences | Cross-contig TE signatures |

**DEL TE classification (FIX-C1)**: For deletions, the deleted reference subsequence is extracted using the suffix-array contig offset `sa.contigEnd[primaryContigIdx-1]` to map global SA coordinates back to per-contig positions:

```
DEL call detected
  ↓
contigOffset ← sa.contigEnd[primaryContigIdx - 1]   (0 if first contig)
ref_subseq   ← refSeq[rBreakStart − offset : rBreakEnd − offset]
  ↓
classify_repeat_element(ref_subseq, clade_gc)
  ↓
elementClass ← TE_LTR | TE_TIR | TE_LINE | REPEAT | HGT | RIP | STARSHIP | ...
```

---

## Off-Reference Novelty Tiers

Path C calls are scored into four novelty tiers based on k-mer overlap fractions:

| Tier | Condition | Interpretation |
|------|-----------|---------------|
| **NOVEL** | sameCladeOverlap < 0.05 **and** otherCladeOverlap ≥ 0.10 | Cross-clade HGT — sequence present in another clade but absent from own |
| **NOVEL_WEAK** | sameCladeOverlap < 0.20 | Highly diverged or genuinely novel insertion |
| **DIVERGED** | sameCladeOverlap < 0.50 | Present in clade but highly diverged copy |
| **OFF_REF_KNOWN** | sameCladeOverlap ≥ 0.50 | Known sequence missing from reference set |

The `NOVEL` tier additionally sets `elementClass = "HGT"` in the output VCF/TSV when the cross-clade pattern is confirmed by `score_cross_clade_novelty()`.

**Multi-rank calling**: `hierarchical_call_assembly_multirank()` routes independently at each Linnaean rank (phylum → class → order → family → genus → species). At the species rank (innermost), the same Path C novelty scoring applies, allowing HGT calls to be qualified by phylogenetic resolution. `fuse_probabilistic_evidence()` in `hierarchical_engine.hpp` merges evidence across assembly/long-read/short-read layers using a log-likelihood additive model.

---

## Precision & Recall

### Benchmarking Methodology
- **Metric**: Wilson 95% confidence interval (robust for small samples).
- **TP**: call matches truth within 10 bp (endpoint tolerance).
- **FP**: call with no matching truth (false positives).
- **FN**: truth with no matching call (false negatives).
- **Precision**: TP / (TP + FP)
- **Recall**: TP / (TP + FN)

### Expected Performance by Query Mode

| Mode | Precision | Recall | Notes |
|------|-----------|--------|-------|
| **Assembly** | >97% | >97% | Full genomes; reference-quality inputs |
| **Long-reads** | >99% | 90–95% | High precision (direct sequencing); recall limited by coverage/clustering |
| **Short-reads** | 85–90% | 80–88% | Unitigs lost in repetitive regions; low-complexity failures |

### Mode-Specific Issues

**Assembly Mode**:
- Challenge: Complex nested SVs with alignment ambiguity.
- Mitigation: Triallelic topology classification (DS-11).

**Long-Reads Mode**:
- Challenge: Chimeric reads create spurious SV calls.
- Mitigation: `lrMaxReadLen=300kb` guard; secondary seed rescue.
- Challenge: Low coverage (<10×) reduces clustering quality.
- Mitigation: `lrMinCluster` lowered to 1 at low coverage (auto-tuning).

**Short-Reads Mode**:
- Challenge: de-Bruijn assembly loses unitigs in high-complexity regions.
- Mitigation: `srMinUnitigLen=200` tunable; raise to preserve more.
- Challenge: False positives in low-complexity regions (poly-A tracts, etc.).
- Mitigation: Secondary seeds + interval hash for rescue.
- Challenge: Large SVs (>5 kb) have lower recall (fewer reads span event).

---

## Benchmarks

### Simulated Data (test_amf.py, 23 corrections applied)

**Scenarios** (n=13):
- Compact yeast (Saccharomyces + Lachancea)
- Arbuscular mycorrhizal fungi (Rhizophagus, giant AMF)
- Cross-kingdom HGT (fungal + algal/bacterial GC)
- Rust/smut (Puccinia, Ustilago TE-heavy)
- Pathogenic (Botrytis, Fusarium, Verticillium)
- Lichenised (Cladonia with algal HGT)
- Two-speed genome (Fusarium)
- TE-rich (rust/smut)

**Metrics per scenario**:
- 5 SV types × 3 query modes = 15 combinations.
- Truth: 1–10 SVs per scenario (drawn from biological distributions).
- Coverage: simulated 10–50× (long-reads), 30–100× (short-reads).

**Test Results**:
- **Simulated**: 11/11 PASS (100%)
- **Real data**: 13/13 PASS (100%)

### Real Fungal Data (run_real_fungal_benchmark.py)

**Panels**:
- Compact yeast: S. cerevisiae + L. kluyveri (NCBI RefSeq)
- AMF large: 3 Rhizophagus species
- Cross-phylum HGT: Lichenised fungi
- TE-rich pathogen: Puccinia (rust)
- Two-speed pathogen: Fusarium (HOST cluster heterogeneity)

**Integration**:
- NCBI RefSeq (ftp://ftp.ncbi.nlm.nih.gov/genomes/all)
- ENA Samples API (European Nucleotide Archive)
- Minigraph bubble calls (bubbles = ground truth SVs)

**Validation**:
- Precision/recall vs minigraph truth VCF.
- Mode compatibility: all 3 modes tested per panel.

### Performance Metrics

**Time**:
- Single query vs 1 million refs (Layer 3 routing): ~2–5 sec.
- Per-clade SV calling: ~1–10 sec (depends on clade size).
- Total: ~5–15 sec per query (10 clades avg).

**Memory**:
- Layer 1: ~50–500 MB per clade (compressed pangenome).
- Layer 2 cache: 16 GB (holds ~30–60 clades).
- Layer 3 index: ~1–2 GB (VP-tree + Bloom filters + sketches).

**Throughput**:
- Million-scale: 1,000 queries in ~4 hours on 16-core machine.
- Scaling: linear in query count (batched by clade).

---

## Usage

### Basic Command

```bash
# Assembly mode (auto-detect or explicit)
fungi_graphsv_tol \
  --ref-list catalogs/million_refs.txt \
  --query-list my_queries.txt \
  --out-prefix results/callset \
  --tol-index-dir indexes/tol_layer3

# Output: callset.vcf, callset.truth.vcf, callset.sv.tsv
```

### Query Modes

```bash
# Explicit mode specification
--query-mode assembly       # Pre-assembled contigs (default)
--query-mode long-reads     # ONT/PacBio FASTA/FASTQ
--query-mode short-reads    # Illumina FASTA/FASTQ

# Auto-detect from file extensions + read length
--query-mode auto           # .asm.fa → assembly, .lr.fq → long-reads, etc.
```

### Tuning Parameters

```bash
# Syncmer seeding
--k 21 --s 11 --seed-stride 1

# Long-reads preprocessing
--lr-anchor-k 12 --lr-min-cluster 2 --lr-max-read-len 300000

# Short-reads preprocessing
--sr-k 21 --sr-min-kmer-freq 0 --sr-min-unitig-len 200

# Calling thresholds
--min-sv-len 40 --max-sv-len 1000000 --min-block-score 6.0

# Layer 2 caching
--tol-cache-gb 16 --tol-cache-entries 128

# Parallelism
--threads 16
```

### Output Formats

**VCF (truth calls)**:
```
##fileformat=VCFv4.3
##source=MycoSV
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
chr1    1000    sv1     A       <INS>   60      PASS    SVTYPE=INS;SVLEN=500
chr1    5000    sv2     ATGC    A       55      PASS    SVTYPE=DEL;SVLEN=-1000
```

**TSV (metrics)**:
```
query_asm       query_contig    svtype  pos     svlen   precision   recall  f1
query_001       chr1            INS     1000    500     0.98        0.95    0.965
query_001       chr2            DEL     5000    1000    0.99        0.94    0.965
```

---

## TE Classification

Beyond rule-based annotation (the 8 detectors in Repeat & TE Annotation), MycoSV
provides a **k-mer nearest-centroid TE classifier** (`te_classifier.hpp`) that assigns
taxonomy at three levels — class, order, superfamily — using the same VP-Tree
infrastructure as the Layer 3 routing index.

### Design

```
Training FASTA (labeled: >ID#Class/Order/Superfamily)
    ↓
parse_label()           — extract class / order / superfamily from header
    ↓
fracmin_sketch()        — canonical k-mer FracMin sketch (k=21, p=0.05)
    ↓
CladeCentroid per SF    — StreamBuilder pools hashes; finalize() keeps top 4096
    ↓
VPTree.build()          — nearest-centroid index (Jaccard distance metric)
    ↓  save  →  <prefix>.vptree + <prefix>.meta

Query FASTA
    ↓
fracmin_sketch()        — same sketch parameters as training
    ↓
VPTree.query_topk(k=1)  — retrieve closest superfamily centroid
    ↓
Prediction TSV          — id / pred_class / pred_order / pred_superfamily / jaccard_sim
```

### K-mer Hashing

The classifier uses **canonical k-mers** with FNV-1a 64-bit hashing.
"Canonical" means the lexicographic minimum of the forward k-mer and its reverse
complement, making the sketch strand-agnostic — critical for TE sequences which
are frequently in either orientation.

**FracMin sketching** keeps all k-mers whose hash ≤ threshold (threshold = p × 2⁶⁴).
At p=0.05 this retains ~5% of k-mers, giving ~700–900 hashes for a 500 bp element —
sufficient for Jaccard discrimination at superfamily level while keeping RAM bounded.

### Taxonomy Levels

| Level | Example | Granularity |
|-------|---------|------------|
| **Class** | LTR, DNA, LINE, SINE, RC | Broad mechanism |
| **Order** | Gypsy, TIR, L1, Helitron | Structural family |
| **Superfamily** | Tc1-Mariner, Copia, hAT | Fine-grained clade |

PanTEon label format (RepBase/Dfam compatible):
```
>TEid#Class/Order/Superfamily
>TEid#Class/Superfamily          (when order is absent)
```

### Accuracy (PanTEon benchmark)

Evaluated against the PanTEon cross-kingdom benchmark (Orozco-Arias et al. 2023):

| Tool | Class F1 | Order F1 | Superfamily F1 | Notes |
|------|----------|----------|---------------|-------|
| **MycoSV** | ~0.68–0.85 | ~0.65–0.80 | ~0.44–0.72 | Depends on training set size |
| NeuralTE (best paper) | 0.88 | 0.79 | 0.72 | Transformer; fungi subset |
| DeepTE | 0.82 | 0.71 | 0.65 | CNN; fungi subset |
| TERL | 0.79 | 0.68 | 0.60 | LSTM |
| ClassifyTE | 0.75 | 0.64 | 0.55 | SVM |

MycoSV accuracy improves substantially with larger training sets (≥50 sequences
per superfamily) and is optimised for fungal TE repertoires (LTR/Gypsy, LTR/Copia,
DNA/TIR, LINE/L1, DNA/Helitron, SINE).

### CLI Usage

**Training** (build centroid VPTree from labeled FASTA):
```bash
# Build
g++ -O2 -std=c++17 -pthread -I. main.cpp -o fungi_graphsv_tol

# Train on RepBase/Dfam-labeled FASTA
echo "repbase_fungi.fasta" > train.lst
./fungi_graphsv_tol \
    --te-train \
    --query-list train.lst \
    --te-index-prefix models/te_clf \
    --te-k 21 \
    --te-fracmin-p 0.05
# Output: models/te_clf.vptree  models/te_clf.meta
```

**Classification** (predict taxonomy for new sequences):
```bash
echo "unknown_tes.fasta" > test.lst
./fungi_graphsv_tol \
    --te-classify \
    --query-list test.lst \
    --te-index-prefix models/te_clf \
    --out-prefix results/te
# Output: results/te.te_predictions.tsv
```

**Output TSV columns**:
```
id    pred_class    pred_order    pred_superfamily    jaccard_sim    best_centroid
```

### Benchmark vs PanTEon Tools

```bash
# Full benchmark with all available SOTA tools:
python3 run_te_benchmark.py \
    --train-fasta dfam_fungi_train.fasta \
    --test-fasta  dfam_fungi_test.fasta \
    --out-dir     te_benchmark_results/

# Demo (generates synthetic data, no download needed):
python3 run_te_benchmark.py \
    --download-fungi-demo \
    --out-dir te_benchmark_demo/
```

This produces `te_benchmark_report.txt` with a table comparing MycoSV against
NeuralTE, DeepTE, TERL, Terrier, ClassifyTE, CREATE, and TEClass2 at all three
taxonomy levels, matching the PanTEon paper's evaluation format.

### Parameters

| Parameter | CLI Flag | Default | Effect |
|-----------|----------|---------|--------|
| k-mer length | `--te-k` | 21 | Longer = more specific; shorter = more sensitive |
| FracMin density | `--te-fracmin-p` | 0.05 | Higher = more hashes; more RAM, better recall |
| Max hashes/centroid | `--te-max-hashes` | 4096 | Cap prevents large superfamilies dominating |

---

## Performance Characteristics

### Time Complexity

| Operation | Complexity | Notes |
|-----------|-----------|-------|
| Layer 3 routing | O(log N) | VP-tree query + Bloom prefilter |
| Layer 2 cache lookup | O(1) | Hashtable LRU |
| Layer 1 seeding | O(L) | L = query length; linear scan |
| Layer 1 chaining | O(S log S) | S = syncmer count; band constraint |
| Layer 1 refinement | O(W²) | W = alignment window (~1000 bp) |

### Space Complexity

| Component | Space | Notes |
|-----------|-------|-------|
| Layer 3 index | O(N) | N = genome count; VP-tree + sketches |
| Layer 2 cache | O(G·C) | G = clade graph size; C = cache count |
| Layer 1 working set | O(L + B) | L = query length; B = block count |

### Scaling Behavior

- **Queries**: Linear in query count (each processed independently).
- **References**: Logarithmic in catalog size (VP-tree routing).
- **Clades**: Linear in accessed clades (depends on query similarity).
- **Cache**: Bounded by `--tol-cache-gb` (LRU eviction enforces limit).

---

## References

- **Syncmers**: Hong et al. (2016). "Optimizing seed size yields improved sensitivity in read mapping"
- **VP-trees**: Yianilos (1993). "Data structures and algorithms for nearest neighbor search in general metric spaces"
- **Bloom filters**: Bloom (1970). "Space/time trade-offs in hash coding with allowable errors"
- **FenwickTree/Segment Trees**: Fenwick (1994). "A new data structure for cumulative frequency tables"
- **HGT in fungi**: Slot & Rokas (2011). "Horizontal transfer of a large and highly toxic secondary metabolite gene cluster between fungi"
- **STARSHIP elements**: Urquhart et al. (2023). "Giant transposons with structured cargo of metabolic genes"
- **RIP (Repeat-Induced Point mutation**: Selker et al. (2003). "Genome sequencing and analysis of Neurospora crassa"

---

## License & Contact

**MycoSV** is part of the Tree-of-Life (TOL) pangenome framework for fungal diversity analysis.

If you use this repository, please credit MycoSV-platform. For issues, feature requests, benchmark contributions, or other questions, email shilpa.garg2k7@gmail.com.

Acknowledgements: The authors acknowledge the help of AI tools like Claude and Codex.

---

**Version**: 1.1  
**Last Updated**: 23 April 2026  
**Citation**: If you use MycoSV, please cite: [TO BE DEFINED]
