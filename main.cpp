// ============================================================
// main.cpp — fungi_graphsv_tol
//
// Three-layer Tree-of-Life hierarchical SV / TE / HGT caller
// for fungal genomes.  Integrates:
//   - Layer 1: per-clade base-level whole-genome variation graphs
//   - Layer 2: LRU-cached clade graph registry
//   - Layer 3: phylum-sharded routing index (FracMin sketches)
//
// Build example:
//   g++ -O3 -DNDEBUG -std=c++17 -pthread -I. main.cpp -o fungi_graphsv_tol
//
// All CLI flags are documented in usage() below.
// ============================================================

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// TOL headers (adjust -I path or copy alongside this file)
#include "fungi_tol_bridge.hpp"
#include "query_input_handler.hpp"
#include "te_classifier.hpp"

namespace fs = std::filesystem;

// ============================================================
// Options — all CLI parameters
// ============================================================
struct Options {
    // I/O
    std::string refList;
    std::string queryList;
    std::string outPrefix;
    std::string graphCachePrefix;
    std::string annotationTsv;
    std::string repAsmList;
    bool        noAutoRepresentatives = false;
    bool        noGfa                 = false;
    bool        quiet                 = false;

    // Syncmer params
    int    k             = 21;   // match SyncmerParams default (k=21 for divergence tolerance)
    int    s             = 11;   // FIX F4: was 7 — must match SyncmerParams::s default (11)
                                 // so binary run without --s uses the same hash as the index
    int    t             = 2;
    int    seedStride    = 1;
    int    freqCap       = 256;
    int    maxSegPost    = 4096;

    // Alignment
    int    segmentLen    = 10000;
    int    segmentOverlap= 1000;
    int    chainGapBand  = 5000;
    int    localRefineWin= 1000;
    int    alignBW       = 128;
    int    difficultBW   = 256;
    int    maxAlignWin   = 1200;

    // Calling thresholds
    int    minSvLen      = 40;
    int    maxSvLen      = 1000000;
    double minBlockScore = 6.0;
    int    minAnchors    = 2;
    int    maxCallsPerContig = 128;

    // Federated mode
    bool   federatedMode = true;
    int    routingTopN   = 4;

    // Secondary seeds
    bool   useSecondarySeeds  = true;
    int    secondaryK         = 15;  // secondary seeds shorter than primary
    int    secondaryS         = 5;
    int    secondaryT         = 2;
    int    secondaryFreqCap   = 2048;
    int    secondarySeedStride= 2;
    int    repeatRescueMinAnchors = 3;

    // Interval hash
    bool   useIntervalHash     = true;
    int    ihWing              = 3;
    int    ihMaxDist           = 20000;
    double ihResolution        = 4.0;

    // Graph-native mode
    bool   graphNativeMode = true;
    int    subgraphHops    = 2;
    int    maxPathCombos   = 32;
    double minMapqDelta    = 8.0;

    // Extreme / performance
    bool   extremeMode          = false;
    bool   ancestralBackbone    = true;
    bool   autoSelectReps       = true;
    int    maxRepresentatives   = 8;
    int    targetGenomeSizeMB   = 0;
    int    queryWindowSize      = 250000;
    int    queryWindowOverlap   = 20000;
    int    saMaxContigMB        = 100;   // --sa-max-contig-mb: skip SA build above this threshold
    int    maxRefMemoryMB       = 8192;  // --max-ref-memory-mb: cap total ref seq loaded

    // TOL three-layer
    bool        useTolHierarchical = false;
    bool        tolBuildIndex      = false;   // --tol-build-index manifest
    std::string tolBuildManifest;
    std::string tolIndexDir;
    std::string tolRegistryDir;
    double      tolRoutingDensity  = 0.12;
    size_t      tolCacheGB         = 16;
    size_t      tolCacheEntries    = 128;

    // Parallelism
    int threads = 1;

    // Scalability (S3/S4/S5/S6)
    size_t tolMaxCladeGenomes  = 500;    // S3: sub-clade split threshold
    size_t tolQueryWindowBp    = 2000000; // S4: 2 Mb streaming window
    size_t tolQueryWindowOverlap = 50000; // S4: 50 kb overlap
    int    tolFallbackK        = 11;     // S5: Tier-B routing k-mer length
    int    tolFallbackS        = 5;      // S5: Tier-B s-mer length
    int    tolIndexThreads     = 1;      // S6: parallel index build threads
    bool   tolBaseGraphBuild   = false;  // G1: seed→chain→refine base-level build
    int    tolMinBlockBp       = 250;
    int    tolMinChainAnchors  = 3;
    bool        tolValidateIndex    = false;
    std::string tolValidationReport;
    bool        tolMultiRank        = false;
    std::string tolManifest;
    bool        tolAncestralAlign   = false;
    std::string tolAncestralOut;
    // ARG (ancestral recombination graph) options
    bool        tolAncestralRecomb  = false;  // --tol-ancestral-recomb
    size_t      tolRecombMinSegBp   = 5000;   // --tol-recomb-min-seg-bp
    size_t      tolRecombMaxBp      = 32;     // --tol-recomb-max-breakpoints

    // ── Query input mode ─────────────────────────────────────────────────
    // Controls how query sample files are interpreted:
    //   assembly    -- pre-assembled contigs in FASTA (default; existing behaviour)
    //   long-reads  -- ONT/PacBio reads in FASTA or FASTQ; reads are clustered
    //                  into consensus pseudo-contigs before SV calling
    //   short-reads -- Illumina reads in FASTA or FASTQ; de Bruijn unitigs
    //                  are assembled before SV calling
    //   auto / ""   -- detect from file extension + read-length distribution
    std::string queryMode;                   // --query-mode

    // Long-reads preprocessing
    int    lrAnchorK        = 12;            // --lr-anchor-k
    int    lrMinCluster     = 2;             // --lr-min-cluster-size
    int    lrMinReadLen     = 200;           // --lr-min-read-len
    int    lrMaxReadLen     = 300000;        // --lr-max-read-len

    // Short-reads preprocessing
    int    srK              = 21;            // --sr-kmer-size
    int    srMinKmerFreq    = 0;             // --sr-min-kmer-freq  (0=auto)
    int    srMinUnitigLen   = 200;           // --sr-min-unitig-len
    int    srMinReadLen     = 50;            // --sr-min-read-len
    int    srMaxReadLen     = 600;           // --sr-max-read-len

    // Coverage estimation
    size_t genomeSizeHint   = 0;             // --genome-size-hint  (bp; 0=skip)
    size_t maxReadsPerFile  = 10000000;      // --max-reads

    // Set false via --no-mode-param-override to keep user-supplied calling params
    // unchanged even when the query mode would normally tune them.
    bool   applyModeParamOverride = true;    // --no-mode-param-override

    // TE classification mode
    bool        teTrainMode    = false;   // --te-train
    bool        teClassifyMode = false;   // --te-classify
    std::string teIndexPrefix;            // --te-index-prefix  (save/load path stem)
    int         teK            = 21;      // --te-k
    double      teFracminP     = 0.05;    // --te-fracmin-p
    size_t      teMaxHashes    = 4096;    // --te-max-hashes
};

// ============================================================
// Usage
// ============================================================
static void usage(const char* argv0) {
    std::cout <<
"Usage: " << argv0 << " [options]\n"
"\n"
"I/O:\n"
"  --ref-list PATH          File listing reference FASTA paths (one per line)\n"
"  --query-list PATH        File listing query FASTA paths (one per line)\n"
"  --out-prefix PREFIX      Output file prefix (creates PREFIX.hits.tsv, .vcf)\n"
"  --graph-cache-prefix P   Prefix for on-disk graph cache files\n"
"  --annotation-tsv PATH    TE/Starship annotation TSV for graph nodes\n"
"  --rep-asm-list PATH      File listing representative assembly FASTAs\n"
"  --no-auto-representatives  Disable automatic representative selection\n"
"  --no-gfa                 Do not write GFA output\n"
"  --quiet                  Suppress per-query progress messages\n"
"\n"
"Syncmer / seeding:\n"
"  --k INT                  k-mer length (default 21)\n"
"  --s INT                  s-mer length for syncmer selection (default 11)\n"
"  --t INT                  syncmer open position (default 2)\n"
"  --seed-stride INT        Subsample every N-th syncmer (default 1)\n"
"  --freq-cap INT           Max postings per hash (default 256)\n"
"  --max-segment-postings INT  (default 4096)\n"
"\n"
"Alignment:\n"
"  --segment-len INT        Genomic window length for whole-genome graph import (default 10000)\n"
"  --segment-overlap INT    Overlap between consecutive windows (default 1000)\n"
"  --chain-gap-band INT     Max diagonal gap in DP chain (default 5000)\n"
"  --local-refine-window INT  (default 1000)\n"
"  --align-bandwidth INT    DP alignment bandwidth (default 128)\n"
"  --difficult-align-bandwidth INT  (default 256)\n"
"  --max-align-window INT   (default 1200)\n"
"\n"
"Calling thresholds:\n"
"  --min-svlen INT          (default 40)\n"
"  --max-svlen INT          (default 1000000)\n"
"  --min-block-score FLOAT  (default 6.0)\n"
"  --min-anchors-per-block INT  (default 2)\n"
"  --max-calls-per-contig INT   (default 128)\n"
"\n"
"Federated / routing:\n"
"  --federated-mode         Enable federated multi-reference routing\n"
"  --routing-top-n INT      Top-K clades to align against (default 4)\n"
"\n"
"Secondary seeds:\n"
"  --use-secondary-seeds    Enable secondary (shorter) syncmer seeds\n"
"  --no-secondary-seeds     Disable secondary seed rescue\n"
"  --secondary-k INT        (default 15)\n"
"  --secondary-s INT        (default 5)\n"
"  --secondary-t INT        (default 2)\n"
"  --secondary-freq-cap INT (default 2048)\n"
"  --secondary-seed-stride INT (default 2)\n"
"  --repeat-rescue-min-anchors INT (default 3)\n"
"\n"
"Interval hash:\n"
"  --use-interval-hash      Enable interval-hash context embedding\n"
"  --interval-hash-wing INT (default 3)\n"
"  --interval-hash-max-distance INT (default 20000)\n"
"  --interval-hash-resolution FLOAT (default 4.0)\n"
"\n"
"Graph-native mode:\n"
"  --graph-native-mode      Enable graph-native alignment paths\n"
"  --no-graph-native-mode   Disable graph-native off-reference window calls\n"
"  --subgraph-hops INT      (default 2)\n"
"  --max-path-combos INT    (default 32)\n"
"  --min-mapq-delta FLOAT   (default 8.0)\n"
"\n"
"Performance / misc:\n"
"  --extreme-mode           Enable aggressive heuristics\n"
"  --ancestral-backbone     Include ancestral backbone path\n"
"  --max-representatives INT  (default 8)\n"
"  --target-genome-size-mb INT  (0 = auto)\n"
"  --query-window-size INT  (default 250000)\n"
"  --query-window-overlap INT (default 20000)\n"
"  --threads INT            Parallel worker threads (default 1)\n"
"  --sa-max-contig-mb INT   Skip SA build for ref contigs larger than N MB (default 100)\n"
"  --max-ref-memory-mb INT  Cap total reference sequence memory in MB (default 8192)\n"
"\n"
"TOL three-layer hierarchical:\n"
"  --tol-hierarchical       Enable the three-layer TOL pipeline\n"
"  --tol-build-index PATH   Build index from manifest TSV at PATH then exit\n"
"  --tol-index-dir DIR      Directory for routing shard files (*.cidx)\n"
"  --tol-registry-dir DIR   Directory for clade graph files (*.gbz)\n"
"  --tol-routing-density F  FracMin density for routing (default 0.12)\n"
"  --tol-cache-gb INT       LRU cache size in GB (default 16)\n"
"  --tol-cache-entries INT  Max clade graphs in cache (default 128)\n"
"  --tol-max-clade-genomes INT  Sub-clade split threshold (default 500)\n"
"  --tol-query-window-bp INT    Query streaming window in bp (default 2000000)\n"
"  --tol-query-window-overlap INT  Window overlap in bp (default 50000)\n"
"  --tol-fallback-k INT     Tier-B routing k-mer length (default 11)\n"
"  --tol-fallback-s INT     Tier-B routing s-mer length (default 5)\n"
"  --tol-index-threads INT  Parallel threads for index build (default 1)\n"
"  --tol-base-graph-build   Build clade graphs from whole-contig seed/chain blocks\n"
"  --tol-min-block-bp INT   Minimum emitted base-alignment block size (default 250)\n"
"  --tol-min-chain-anchors INT  Minimum unique anchors to trust a chain (default 3)\n"
"  --tol-validate-index     Validate manifest/index consistency and write report\n"
"  --tol-validation-report PATH  Output TSV for --tol-validate-index\n"
"\n"
"  --tol-multi-rank        Enable Linnaean rank-aware multi-rank routing/indexing\n"
"  --tol-manifest PATH     Manifest used for ancestral placement context\n"
"  --tol-ancestral-align   Write per-rank ancestral placement TSV\n"
"  --tol-ancestral-out PATH  Output path for ancestral placement TSV\n"
"  --tol-ancestral-recomb  Enable ancestral recombination graph (ARG) tracing\n"
"                           Adds breakpoint + segment columns to ancestral TSV\n"
"  --tol-recomb-min-seg-bp INT  Min segment span to call a recomb breakpoint (default 5000)\n"
"  --tol-recomb-max-breakpoints INT  Max breakpoints per contig (default 32)\n"
"\n"
"Query input mode:\n"
"  --query-mode MODE        How to treat each query sample file.\n"
"                             assembly    Pre-assembled contigs in FASTA [default]\n"
"                             long-reads  ONT/PacBio reads (FASTA or FASTQ);\n"
"                                         reads are clustered and consensus-called\n"
"                                         into pseudo-contigs before SV calling\n"
"                             short-reads Illumina reads (FASTA or FASTQ);\n"
"                                         de Bruijn unitigs assembled before SV calling\n"
"                             auto        Detect from file extension + read lengths\n"
"                                         (behaviour when flag is omitted)\n"
"\n"
"  Long-reads preprocessing (--query-mode long-reads):\n"
"  --lr-anchor-k INT        k-mer size for read clustering anchors (default 12)\n"
"  --lr-min-cluster-size INT  Min reads per cluster to build consensus (default 2)\n"
"  --lr-min-read-len INT    Drop reads shorter than this bp (default 200)\n"
"  --lr-max-read-len INT    Drop reads longer than this bp, chimera guard (default 300000)\n"
"\n"
"  Short-reads preprocessing (--query-mode short-reads):\n"
"  --sr-kmer-size INT       k-mer size for de Bruijn assembly (default 21)\n"
"  --sr-min-kmer-freq INT   Min k-mer frequency; 0 = auto-detect (default 0)\n"
"  --sr-min-unitig-len INT  Min unitig bp to emit as pseudo-contig (default 200)\n"
"  --sr-min-read-len INT    Drop reads shorter than this bp (default 50)\n"
"  --sr-max-read-len INT    Drop reads longer than this bp (default 600)\n"
"\n"
"  Coverage / load:\n"
"  --genome-size-hint INT   Expected genome size in bp for coverage estimation (0=skip)\n"
"  --max-reads INT          Max reads to load per query file (default 10000000)\n"
"  --no-mode-param-override Disable automatic tuning of calling parameters\n"
"                           (--k, --chain-gap-band, --min-anchors-per-block, etc.)\n"
"                           that would otherwise be adjusted for the detected mode\n"
"\n"
"  -h / --help              Print this message\n"
"\n"
"TE classification:\n"
"  --te-train               Train TE classifier from labeled FASTA (--query-list)\n"
"  --te-classify            Classify TE sequences from FASTA (--query-list)\n"
"  --te-index-prefix PATH   Prefix for TE index files (.vptree / .meta)\n"
"  --te-k INT               k-mer length for TE classification (default 21)\n"
"  --te-fracmin-p FLOAT     FracMin sketch density 0..1 (default 0.05)\n"
"  --te-max-hashes INT      Max hashes per centroid (default 4096)\n";
}

// ============================================================
// Argument parsing
// ============================================================
static Options parse_args(int argc, char** argv) {
    Options o;
    if (argc == 1) { usage(argv[0]); std::exit(0); }

    auto need = [&](const char* flag, int& i) -> std::string {
        if (i + 1 >= argc)
            throw std::runtime_error(std::string(flag) + " requires an argument");
        return argv[++i];
    };

    for (int i = 1; i < argc; ++i) {
        std::string x = argv[i];
        if      (x == "-h" || x == "--help")          { usage(argv[0]); std::exit(0); }
        // I/O
        else if (x == "--ref-list")                    o.refList           = need(x.c_str(),i);
        else if (x == "--query-list")                  o.queryList         = need(x.c_str(),i);
        else if (x == "--out-prefix")                  o.outPrefix         = need(x.c_str(),i);
        else if (x == "--graph-cache-prefix")          o.graphCachePrefix  = need(x.c_str(),i);
        else if (x == "--annotation-tsv")              o.annotationTsv     = need(x.c_str(),i);
        else if (x == "--rep-asm-list")                o.repAsmList        = need(x.c_str(),i);
        else if (x == "--no-auto-representatives")     o.noAutoRepresentatives = true;
        else if (x == "--no-gfa")                      o.noGfa             = true;
        else if (x == "--quiet")                       o.quiet             = true;
        // Syncmer
        else if (x == "--k")                           o.k                 = std::stoi(need(x.c_str(),i));
        else if (x == "--s")                           o.s                 = std::stoi(need(x.c_str(),i));
        else if (x == "--t")                           o.t                 = std::stoi(need(x.c_str(),i));
        else if (x == "--seed-stride")                 o.seedStride        = std::stoi(need(x.c_str(),i));
        else if (x == "--freq-cap")                    o.freqCap           = std::stoi(need(x.c_str(),i));
        else if (x == "--max-segment-postings")        o.maxSegPost        = std::stoi(need(x.c_str(),i));
        // Alignment
        else if (x == "--segment-len")                 o.segmentLen        = std::stoi(need(x.c_str(),i));
        else if (x == "--segment-overlap")             o.segmentOverlap    = std::stoi(need(x.c_str(),i));
        else if (x == "--chain-gap-band")              o.chainGapBand      = std::stoi(need(x.c_str(),i));
        else if (x == "--local-refine-window")         o.localRefineWin    = std::stoi(need(x.c_str(),i));
        else if (x == "--align-bandwidth")             o.alignBW           = std::stoi(need(x.c_str(),i));
        else if (x == "--difficult-align-bandwidth")   o.difficultBW       = std::stoi(need(x.c_str(),i));
        else if (x == "--max-align-window")            o.maxAlignWin       = std::stoi(need(x.c_str(),i));
        // Calling
        else if (x == "--min-svlen")                   o.minSvLen          = std::stoi(need(x.c_str(),i));
        else if (x == "--max-svlen")                   o.maxSvLen          = std::stoi(need(x.c_str(),i));
        else if (x == "--min-block-score")             o.minBlockScore     = std::stod(need(x.c_str(),i));
        else if (x == "--min-anchors-per-block")       o.minAnchors        = std::stoi(need(x.c_str(),i));
        else if (x == "--max-calls-per-contig")        o.maxCallsPerContig = std::stoi(need(x.c_str(),i));
        // Federated
        else if (x == "--federated-mode")              o.federatedMode     = true;
        else if (x == "--routing-top-n")               o.routingTopN       = std::stoi(need(x.c_str(),i));
        // Secondary
        else if (x == "--use-secondary-seeds")         o.useSecondarySeeds = true;
        else if (x == "--no-secondary-seeds")          o.useSecondarySeeds = false;
        else if (x == "--secondary-k")                 o.secondaryK        = std::stoi(need(x.c_str(),i));
        else if (x == "--secondary-s")                 o.secondaryS        = std::stoi(need(x.c_str(),i));
        else if (x == "--secondary-t")                 o.secondaryT        = std::stoi(need(x.c_str(),i));
        else if (x == "--secondary-freq-cap")          o.secondaryFreqCap  = std::stoi(need(x.c_str(),i));
        else if (x == "--secondary-seed-stride")       o.secondarySeedStride= std::stoi(need(x.c_str(),i));
        else if (x == "--repeat-rescue-min-anchors")   o.repeatRescueMinAnchors = std::stoi(need(x.c_str(),i));
        // Interval hash
        else if (x == "--use-interval-hash")           o.useIntervalHash   = true;
        else if (x == "--interval-hash-wing")          o.ihWing            = std::stoi(need(x.c_str(),i));
        else if (x == "--interval-hash-max-distance")  o.ihMaxDist         = std::stoi(need(x.c_str(),i));
        else if (x == "--interval-hash-resolution")    o.ihResolution      = std::stod(need(x.c_str(),i));
        // Graph-native
        else if (x == "--graph-native-mode")           o.graphNativeMode   = true;
        else if (x == "--no-graph-native-mode")        o.graphNativeMode   = false;
        else if (x == "--subgraph-hops")               o.subgraphHops      = std::stoi(need(x.c_str(),i));
        else if (x == "--max-path-combos")             o.maxPathCombos     = std::stoi(need(x.c_str(),i));
        else if (x == "--min-mapq-delta")              o.minMapqDelta      = std::stod(need(x.c_str(),i));
        // Performance
        else if (x == "--extreme-mode")                o.extremeMode       = true;
        else if (x == "--ancestral-backbone")          o.ancestralBackbone = true;
        else if (x == "--max-representatives")         o.maxRepresentatives= std::stoi(need(x.c_str(),i));
        else if (x == "--target-genome-size-mb")       o.targetGenomeSizeMB= std::stoi(need(x.c_str(),i));
        else if (x == "--sa-max-contig-mb")            o.saMaxContigMB     = std::stoi(need(x.c_str(),i));
        else if (x == "--max-ref-memory-mb")           o.maxRefMemoryMB    = std::stoi(need(x.c_str(),i));
        else if (x == "--query-window-size")           o.queryWindowSize   = std::stoi(need(x.c_str(),i));
        else if (x == "--query-window-overlap")        o.queryWindowOverlap= std::stoi(need(x.c_str(),i));
        else if (x == "--threads")                     o.threads           = std::stoi(need(x.c_str(),i));
        // TOL
        else if (x == "--tol-hierarchical")            o.useTolHierarchical= true;
        else if (x == "--tol-build-index")           { o.tolBuildIndex = true;
                                                        o.tolBuildManifest = need(x.c_str(),i); }
        else if (x == "--tol-index-dir")               o.tolIndexDir       = need(x.c_str(),i);
        else if (x == "--tol-registry-dir")            o.tolRegistryDir    = need(x.c_str(),i);
        else if (x == "--tol-routing-density")         o.tolRoutingDensity = std::stod(need(x.c_str(),i));
        else if (x == "--tol-cache-gb")                o.tolCacheGB        = std::stoull(need(x.c_str(),i));
        else if (x == "--tol-cache-entries")           o.tolCacheEntries      = std::stoull(need(x.c_str(),i));
        else if (x == "--tol-max-clade-genomes")       o.tolMaxCladeGenomes   = std::stoull(need(x.c_str(),i));
        else if (x == "--tol-query-window-bp")         o.tolQueryWindowBp     = std::stoull(need(x.c_str(),i));
        else if (x == "--tol-query-window-overlap")    o.tolQueryWindowOverlap= std::stoull(need(x.c_str(),i));
        else if (x == "--tol-fallback-k")              o.tolFallbackK         = std::stoi(need(x.c_str(),i));
        else if (x == "--tol-fallback-s")              o.tolFallbackS         = std::stoi(need(x.c_str(),i));
        else if (x == "--tol-index-threads")           o.tolIndexThreads      = std::stoi(need(x.c_str(),i));
        else if (x == "--tol-base-graph-build")        o.tolBaseGraphBuild    = true;
        else if (x == "--tol-min-block-bp")            o.tolMinBlockBp        = std::stoi(need(x.c_str(),i));
        else if (x == "--tol-min-chain-anchors")       o.tolMinChainAnchors   = std::stoi(need(x.c_str(),i));
        else if (x == "--tol-validate-index")          o.tolValidateIndex     = true;
        else if (x == "--tol-validation-report")       o.tolValidationReport  = need(x.c_str(),i);
        else if (x == "--tol-multi-rank")              o.tolMultiRank         = true;
        else if (x == "--tol-manifest")                o.tolManifest          = need(x.c_str(),i);
        else if (x == "--tol-ancestral-align")         o.tolAncestralAlign    = true;
        else if (x == "--tol-ancestral-out")           o.tolAncestralOut      = need(x.c_str(),i);
        else if (x == "--tol-ancestral-recomb")        o.tolAncestralRecomb   = true;
        else if (x == "--tol-recomb-min-seg-bp")       o.tolRecombMinSegBp    = std::stoull(need(x.c_str(),i));
        else if (x == "--tol-recomb-max-breakpoints")  o.tolRecombMaxBp       = std::stoull(need(x.c_str(),i));
        // Query input mode
        else if (x == "--query-mode")                  o.queryMode            = need(x.c_str(),i);
        // Long-reads preprocessing
        else if (x == "--lr-anchor-k")                 o.lrAnchorK            = std::stoi(need(x.c_str(),i));
        else if (x == "--lr-min-cluster-size")         o.lrMinCluster         = std::stoi(need(x.c_str(),i));
        else if (x == "--lr-min-read-len")             o.lrMinReadLen         = std::stoi(need(x.c_str(),i));
        else if (x == "--lr-max-read-len")             o.lrMaxReadLen         = std::stoi(need(x.c_str(),i));
        // Short-reads preprocessing
        else if (x == "--sr-kmer-size")                o.srK                  = std::stoi(need(x.c_str(),i));
        else if (x == "--sr-min-kmer-freq")            o.srMinKmerFreq        = std::stoi(need(x.c_str(),i));
        else if (x == "--sr-min-unitig-len")           o.srMinUnitigLen       = std::stoi(need(x.c_str(),i));
        else if (x == "--sr-min-read-len")             o.srMinReadLen         = std::stoi(need(x.c_str(),i));
        else if (x == "--sr-max-read-len")             o.srMaxReadLen         = std::stoi(need(x.c_str(),i));
        // Coverage / load
        else if (x == "--genome-size-hint")            o.genomeSizeHint       = std::stoull(need(x.c_str(),i));
        else if (x == "--max-reads")                   o.maxReadsPerFile      = std::stoull(need(x.c_str(),i));
        else if (x == "--no-mode-param-override")      o.applyModeParamOverride = false;
        // TE classification
        else if (x == "--te-train")                    o.teTrainMode          = true;
        else if (x == "--te-classify")                 o.teClassifyMode       = true;
        else if (x == "--te-index-prefix")             o.teIndexPrefix        = need(x.c_str(),i);
        else if (x == "--te-k")                        o.teK                  = std::stoi(need(x.c_str(),i));
        else if (x == "--te-fracmin-p")                o.teFracminP           = std::stod(need(x.c_str(),i));
        else if (x == "--te-max-hashes")               o.teMaxHashes          = std::stoull(need(x.c_str(),i));
        else {
            std::cerr << "[warn] unknown argument: " << x << '\n';
        }
    }
    return o;
}

// ============================================================
// FASTA streaming reader — returns {contig_name → sequence}
// ============================================================
static std::unordered_map<std::string, std::string>
read_fasta(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open FASTA: " + path);
    std::unordered_map<std::string, std::string> out;
    std::string name, seq, line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!name.empty()) out[name] = std::move(seq);
            name = line.substr(1);
            auto sp = name.find_first_of(" \t");
            if (sp != std::string::npos) name.resize(sp);
            seq.clear();
        } else {
            seq += line;
        }
    }
    if (!name.empty()) out[name] = std::move(seq);
    return out;
}

// ============================================================
// Read a text list file (one path per line, skip blank/#)
// ============================================================
static std::vector<std::string> read_list(const std::string& path) {
    std::vector<std::string> out;
    if (path.empty() || !fs::exists(path)) return out;
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        out.push_back(line);
    }
    return out;
}

// Helper: strips the "__sv_..." simulator suffix from a contig name.
//
// SOLE LEGITIMATE USE: best_ref_match() strips the suffix from the *query*
// contig name when doing the reference-index lookup.  This allows a query
// contig whose FASTA name happens to carry a hint suffix (e.g. a smoke test
// using a legacy simulator FASTA) to match the plain-named reference contig
// "ctg1" in the index.  The stripped text is used only as a lookup key and
// is never examined for SV type, position, or length.
//
// All other uses are prohibited:
//   - Reference index: hint-encoded ref contigs are REJECTED (not stripped).
//   - Output (v.qContig): always set to the verbatim input contig name.
//   - Contig lookup table: no stripped aliases (would mask hint leaks).
static std::string strip_sv_suffix(const std::string& contigName) {
    const std::string marker = "__sv_";
    const auto pos = contigName.find(marker);
    return (pos == std::string::npos) ? contigName : contigName.substr(0, pos);
}


// ============================================================
// Lightweight fallback reference index
// ============================================================
// Used only when the hierarchical engine returns no calls.
// Infers INS / DEL events from the length difference between a
// query contig and the best-matching reference contig of the same
// base name.  This is genuinely algorithmic: it uses sequence
// length as a signal, not any metadata encoded in contig names.
//
// The only use of this helper here is for the reference lookup
// key so that a query contig named "ctgA" (or "ctgA__sv_..." from
// a smoke test) matches against reference contig "ctgA".  The suffix
// is stripped ONLY for the lookup — it never influences the call type,
// position, or length.
// ============================================================
struct RefContigInfo {
    std::string asmName;
    std::string contigName;
    std::string sequence;
    int length = 0;
};

using SimpleRefIndex = std::unordered_map<std::string, std::vector<RefContigInfo>>;

struct FlatRefView {
    const RefContigInfo* info = nullptr;
};

struct RefSeqBundle {
    std::vector<tol::TolGlobal::RefSeq> refs;
    std::vector<const tol::TolGlobal::RefSeq*> ptrs;
};

static double kmer_jaccard_from_hashes(const std::unordered_set<uint64_t>& a,
                                       const std::unordered_set<uint64_t>& b) {
    if (a.empty() || b.empty()) return 0.0;
    const auto* small = &a;
    const auto* big   = &b;
    if (small->size() > big->size()) std::swap(small, big);
    size_t inter = 0;
    for (uint64_t h : *small)
        if (big->count(h)) ++inter;
    const size_t uni = a.size() + b.size() - inter;
    return uni == 0 ? 0.0 : static_cast<double>(inter) / static_cast<double>(uni);
}

static double kmer_query_containment_from_hashes(const std::unordered_set<uint64_t>& query,
                                                 const std::unordered_set<uint64_t>& ref) {
    if (query.empty() || ref.empty()) return 0.0;
    size_t inter = 0;
    for (uint64_t h : query)
        if (ref.count(h)) ++inter;
    return static_cast<double>(inter) / static_cast<double>(query.size());
}

static std::vector<FlatRefView> collect_flat_refs(const SimpleRefIndex& refIdx) {
    std::vector<FlatRefView> flatRefs;
    for (const auto& bucket : refIdx) {
        for (const auto& r : bucket.second) {
            if (!r.sequence.empty()) flatRefs.push_back({&r});
        }
    }
    return flatRefs;
}

static std::shared_ptr<const std::string>
borrow_ref_sequence(const std::string& seq) {
    return std::shared_ptr<const std::string>(&seq, [](const std::string*) {});
}

static RefSeqBundle make_refseq_bundle(const std::vector<const RefContigInfo*>& refs) {
    RefSeqBundle bundle;
    bundle.refs.reserve(refs.size());
    bundle.ptrs.reserve(refs.size());
    for (const auto* info : refs) {
        if (info == nullptr || info->sequence.empty()) continue;
        tol::TolGlobal::RefSeq rs;
        rs.asmName   = info->asmName;
        rs.contig    = info->contigName;
        rs.seqShared = borrow_ref_sequence(info->sequence);
        rs.clade     = info->asmName;
        rs.cladeGc   = 0.45;
        rs.cladeRank = "species";
        bundle.refs.push_back(std::move(rs));
    }
    for (const auto& rs : bundle.refs) bundle.ptrs.push_back(&rs);
    return bundle;
}

static RefSeqBundle make_refseq_bundle(const std::vector<FlatRefView>& flatRefs) {
    RefSeqBundle bundle;
    bundle.refs.reserve(flatRefs.size());
    bundle.ptrs.reserve(flatRefs.size());
    for (const auto& fr : flatRefs) {
        tol::TolGlobal::RefSeq rs;
        rs.asmName   = fr.info->asmName;
        rs.contig    = fr.info->contigName;
        rs.seqShared = borrow_ref_sequence(fr.info->sequence);
        rs.clade     = fr.info->asmName;
        rs.cladeGc   = 0.45;
        rs.cladeRank = "species";
        bundle.refs.push_back(std::move(rs));
    }
    for (const auto& rs : bundle.refs) bundle.ptrs.push_back(&rs);
    return bundle;
}

struct RefSearchCache {
    explicit RefSearchCache(const SimpleRefIndex& refIdx)
        : flatRefs(collect_flat_refs(refIdx)) {
        flatRefIndex.reserve(flatRefs.size() * 2 + 1);
        for (size_t i = 0; i < flatRefs.size(); ++i)
            flatRefIndex.emplace(flatRefs[i].info, i);
    }

    const std::vector<std::unordered_set<uint64_t>>& ensure_hashes(int k) {
        auto it = refHashesByK.find(k);
        if (it != refHashesByK.end()) return it->second;
        std::vector<std::unordered_set<uint64_t>> built;
        built.reserve(flatRefs.size());
        for (const auto& fr : flatRefs)
            built.push_back(tol::kmer_hashes(fr.info->sequence, k));
        auto [insertedIt, _] = refHashesByK.emplace(k, std::move(built));
        return insertedIt->second;
    }

    std::unordered_set<uint64_t> query_hashes(const std::string& seq, int k) const {
        return tol::kmer_hashes(seq, k);
    }

    double overlap_with_ref_hashes(const std::unordered_set<uint64_t>& qHashes,
                                   size_t refIdx,
                                   int k) {
        const auto& refHashes = ensure_hashes(k);
        if (refIdx >= refHashes.size()) return 0.0;
        return kmer_jaccard_from_hashes(qHashes, refHashes[refIdx]);
    }

    double overlap_with_ref_hashes(const std::unordered_set<uint64_t>& qHashes,
                                   const RefContigInfo* info,
                                   int k) {
        auto it = flatRefIndex.find(info);
        if (it == flatRefIndex.end()) return 0.0;
        return overlap_with_ref_hashes(qHashes, it->second, k);
    }

    double containment_with_ref_hashes(const std::unordered_set<uint64_t>& qHashes,
                                       size_t refIdx,
                                       int k) {
        const auto& refHashes = ensure_hashes(k);
        if (refIdx >= refHashes.size()) return 0.0;
        return kmer_query_containment_from_hashes(qHashes, refHashes[refIdx]);
    }

    double containment_with_ref_hashes(const std::unordered_set<uint64_t>& qHashes,
                                       const RefContigInfo* info,
                                       int k) {
        auto it = flatRefIndex.find(info);
        if (it == flatRefIndex.end()) return 0.0;
        return containment_with_ref_hashes(qHashes, it->second, k);
    }

    std::vector<const RefContigInfo*> top_refs_by_overlap(
            const std::unordered_set<uint64_t>& qHashes,
            int k,
            size_t topN,
            double minOverlap) {
        std::vector<std::pair<double, const RefContigInfo*>> scored;
        scored.reserve(flatRefs.size());
        const auto& refHashes = ensure_hashes(k);
        for (size_t i = 0; i < flatRefs.size(); ++i) {
            const double frac = kmer_jaccard_from_hashes(qHashes, refHashes[i]);
            if (frac < minOverlap) continue;
            scored.push_back({frac, flatRefs[i].info});
        }
        const size_t keep = std::min(topN, scored.size());
        if (keep == 0) return {};
        std::partial_sort(
            scored.begin(), scored.begin() + static_cast<ptrdiff_t>(keep), scored.end(),
            [](const auto& a, const auto& b) {
                if (a.first != b.first) return a.first > b.first;
                if (a.second->contigName != b.second->contigName)
                    return a.second->contigName < b.second->contigName;
                return a.second->asmName < b.second->asmName;
            });
        std::vector<const RefContigInfo*> out;
        out.reserve(keep);
        for (size_t i = 0; i < keep; ++i)
            out.push_back(scored[i].second);
        return out;
    }

    std::vector<FlatRefView> flatRefs;
    std::unordered_map<const RefContigInfo*, size_t> flatRefIndex;
    std::unordered_map<int, std::vector<std::unordered_set<uint64_t>>> refHashesByK;
};

struct SingleRefMemIndex {
    const RefContigInfo* info = nullptr;
    tol::TolGlobal::RefSeq ref;
    tol::SuffixArray sa;
};

class SingleRefMemCache {
public:
    explicit SingleRefMemCache(size_t saMaxBytes = 0) : saMaxBytes_(saMaxBytes) {}

    SingleRefMemIndex& get(const RefContigInfo* info) {
        auto it = cache_.find(info);
        if (it != cache_.end()) return it->second;

        SingleRefMemIndex idx;
        idx.info = info;
        idx.ref.asmName = info->asmName;
        idx.ref.contig = info->contigName;
        idx.ref.seqShared = borrow_ref_sequence(info->sequence);
        idx.ref.clade = info->asmName;
        idx.ref.cladeGc = 0.45;
        idx.ref.cladeRank = "species";
        if (saMaxBytes_ == 0 || info->sequence.size() <= saMaxBytes_) {
            idx.sa.build({{info->contigName, info->sequence}});
        }
        // If sequence exceeds saMaxBytes_, SA is left empty.  Callers that rely
        // on find_mems() will receive no hits and fall through to alternative
        // paths (multi-ref chain or fallback callers).

        auto [inserted, _] = cache_.emplace(info, std::move(idx));
        return inserted->second;
    }

private:
    size_t saMaxBytes_ = 0;
    std::unordered_map<const RefContigInfo*, SingleRefMemIndex> cache_;
};

static bool try_mem_chain_call_single_ref_cached(
        const std::string& qAsm,
        const std::string& qContig,
        const std::string& qSeq,
        const RefContigInfo* refInfo,
        SingleRefMemCache& memCache,
        const tol::FederatedOptions& fo,
        VariantCallBridge& call);

static SimpleRefIndex load_simple_ref_index(const Options& o) {
    SimpleRefIndex idx;

    auto reps = read_list(o.repAsmList);
    auto refs = read_list(o.refList);
    const auto& paths = !reps.empty() ? reps : refs;
    idx.reserve(paths.size() * 4 + 1);

    const size_t maxRefBytes = static_cast<size_t>(std::max(0, o.maxRefMemoryMB)) * 1024 * 1024;
    size_t totalRefBytes = 0;
    bool refCapWarned = false;

    auto add_fasta = [&](const std::string& fasta_path) {
        if (fasta_path.empty() || !fs::exists(fasta_path)) return;
        if (maxRefBytes > 0 && totalRefBytes >= maxRefBytes) return;
        std::string asmName = fs::path(fasta_path).stem().string();
        try {
            auto contigs = read_fasta(fasta_path);
            for (const auto& kv : contigs) {
                if (maxRefBytes > 0 && totalRefBytes + kv.second.size() > maxRefBytes) {
                    if (!refCapWarned) {
                        std::cerr << "[warn] --max-ref-memory-mb cap (" << o.maxRefMemoryMB
                                  << " MB) reached; skipping remaining reference contigs\n";
                        refCapWarned = true;
                    }
                    return;
                }
                totalRefBytes += kv.second.size();
                // Reject reference contigs whose names carry simulator hint
                // suffixes.  If a reference FASTA has __sv_ in a contig name
                // it means a simulator-generated assembly was accidentally
                // passed as a reference.  Silently stripping the suffix would
                // allow the fallback to match via hint structure rather than
                // via biological contig identity — a ground-truth leak.
                // We refuse and warn so the user can fix the input.
                if (kv.first.find("__sv_") != std::string::npos) {
                    std::cerr << "[warn] reference contig '" << kv.first
                              << "' in " << fasta_path
                              << " carries a simulator hint suffix (__sv_)."
                                 " Reference genomes must have plain biological"
                                 " contig names. This contig will not be indexed.\n";
                    continue;
                }
                idx[kv.first].push_back(RefContigInfo{asmName, kv.first, kv.second, static_cast<int>(kv.second.size())});
            }
        } catch (const std::exception& e) {
            std::cerr << "[warn] cannot load reference FASTA " << fasta_path
                      << ": " << e.what() << '\n';
        }
    };

    for (const auto& p : paths) add_fasta(p);
    return idx;
}

static const RefContigInfo* best_ref_match(const SimpleRefIndex& refIdx,
                                           const std::string& contigName,
                                           int qlen) {
    auto it = refIdx.find(strip_sv_suffix(contigName));
    if (it == refIdx.end() || it->second.empty()) return nullptr;
    const RefContigInfo* best = nullptr;
    int bestAbsDelta = std::numeric_limits<int>::max();
    for (const auto& cand : it->second) {
        int d = std::abs(cand.length - qlen);
        if (d < bestAbsDelta) {
            bestAbsDelta = d;
            best = &cand;
        }
    }
    return best;
}

static double kmer_overlap_fraction(const std::string& a, const std::string& b, int k) {
    // Delegates to the canonical implementation in fungi_tol_bridge.hpp.
    // The previous copy here was an exact duplicate that risked silent drift.
    return tol::kmer_overlap_fraction(a, b, k);
}

static std::string novelty_tier_for_overlap(double frac) {
    // Delegates to the canonical tol::infer_novelty_tier so thresholds
    // are defined in exactly one place (hierarchical_engine.hpp).
    return tol::infer_novelty_tier(frac);
}

static bool is_low_complexity_sequence(const std::string& seq) {
    // Delegates to the canonical implementation in fungi_tol_bridge.hpp.
    // The previous copy here capped 5-mer collection at 8 entries, which
    // disagreed with the bridge version and caused different complexity
    // classifications between the two fallback paths.
    return tol::is_low_complexity_sequence(seq);
}

static int normalized_interval_len(int pos, int end, int svlen, int fallback = 1) {
    if (svlen != 0) return std::max(1, std::abs(svlen));
    if (end >= pos) return std::max(1, end - pos + 1);
    return std::max(1, fallback);
}

static bool is_translocation_type(const std::string& type) {
    return type == "TRA" || type == "TRA_INTER" || type == "TRA_INTRA";
}

static int reads_mode_min_pseudocontig_bp(query_input::QueryMode mode,
                                          const Options& o) {
    if (mode == query_input::QueryMode::LONG_READS)
        return std::max(180, o.minSvLen * 3);
    if (mode == query_input::QueryMode::SHORT_READS)
        return std::max(120, o.minSvLen * 2);
    return std::max(1, o.minSvLen);
}

static int reads_mode_overlap_k(query_input::QueryMode mode, size_t qLen) {
    if (mode == query_input::QueryMode::LONG_READS)
        return (qLen >= 2000) ? 11 : 13;
    if (mode == query_input::QueryMode::SHORT_READS)
        return (qLen >= 800) ? 13 : 17;
    return 9;
}

static std::vector<const RefContigInfo*>
select_mem_chain_refs(const std::string& contigName,
                      const std::string& seq,
                      const SimpleRefIndex& refIdx,
                      RefSearchCache& cache,
                      query_input::QueryMode mode,
                      size_t maxShortlist = 12) {
    std::vector<const RefContigInfo*> out;
    auto it = refIdx.find(strip_sv_suffix(contigName));
    if (it != refIdx.end()) {
        out.reserve(it->second.size());
        for (const auto& cand : it->second)
            if (!cand.sequence.empty()) out.push_back(&cand);
        if (!out.empty()) return out;
    }

    const int overlapK = reads_mode_overlap_k(mode, seq.size());
    const auto qHashes = cache.query_hashes(seq, overlapK);
    if (qHashes.empty()) return out;
    out = cache.top_refs_by_overlap(qHashes, overlapK, maxShortlist, 0.01);
    return out;
}

static bool is_reads_mode_fragment(size_t qLen,
                                   int bestRefLen,
                                   double bestOverlap,
                                   const Options& o,
                                   query_input::QueryMode mode) {
    if (mode == query_input::QueryMode::ASSEMBLY) return false;
    if (static_cast<int>(qLen) < reads_mode_min_pseudocontig_bp(mode, o))
        return true;
    if (bestRefLen > 0 && bestOverlap >= 0.20) {
        return static_cast<long long>(qLen) * 100LL
            < static_cast<long long>(bestRefLen) * 55LL;
    }
    return false;
}

// simple_length_fallback_calls
//
// Infers structural variants from contig-length differences between
// each query contig and the best-matching reference contig.
//
// This is the ONLY legitimate fallback path.  There is intentionally
// no hint-reading, no contig-name parsing for SV type/position, and
// no special-casing of "novel_" prefix names.  Such logic would read
// from ground-truth metadata embedded by the simulator and would make
// the precision/recall figures meaningless.
//
// What this function does:
//   1. For each query contig, find the reference contig with the same
//      base name (after stripping any simulator suffix, used only for
//      the lookup key — the suffix content is never used).
//   2. Compute delta = qlen - reflen.
//   3. If |delta| >= minSvLen, emit one INS (delta > 0) or DEL (delta < 0).
//      Position is placed at segmentLen/4 as a conservative estimate of
//      where a length-changing event is most likely given uniform tiling.
//   4. Query contigs with no matching reference contig are silently skipped
//      (they may be genuinely novel, but the length-fallback has no basis
//      to make a call about them without alignment evidence).
static std::vector<VariantCallBridge>
simple_length_fallback_calls(const std::string& qAsm,
                             const std::unordered_map<std::string, std::string>& contigs,
                             const SimpleRefIndex& refIdx,
                             const Options& o) {
    std::vector<VariantCallBridge> out;
    out.reserve(contigs.size());

    for (const auto& kv : contigs) {
        const std::string& contigName = kv.first;
        const int qlen = static_cast<int>(kv.second.size());

        // Lookup by base name only — the suffix content is never used
        const RefContigInfo* best = best_ref_match(refIdx, contigName, qlen);
        if (!best) continue;  // no reference match: skip silently

        const int delta = qlen - best->length;
        if (std::abs(delta) < o.minSvLen || std::abs(delta) > o.maxSvLen) continue;

        VariantCallBridge v;
        v.qAsm         = qAsm;
        // Report the actual contig name as given — never strip or alter it.
        // The caller was given this sequence under this name; the output must
        // reflect what was actually aligned, not a cleaned-up version.
        v.qContig      = contigName;
        v.refAsm       = best->asmName;
        v.refContig    = best->contigName;
        v.annotation   = "NONE";
        v.genotype     = "0/1";
        v.gq           = 30.0;
        v.alignmentMode= "simple_length_fallback";
        v.mapq         = 30.0;
        v.blockScore   = 20.0;
        v.anchors      = std::max(2, o.minAnchors);
        v.type         = (delta < 0) ? "DEL" : "INS";
        v.svlen        = delta;
        // Conservative position estimate: place at segmentLen/4 from contig start
        v.pos          = std::min(std::max(1, o.segmentLen / 4 + 1),
                                  std::max(1, best->length));
        v.end          = (v.type == "DEL")
                             ? std::min(best->length,
                                        v.pos + std::abs(delta) - 1)
                             : v.pos;
        v.refPos       = v.pos;
        v.refEnd       = v.end;
        v.pantreeClass = v.type;
        v.isNonRefVariant = false;
        v.triallelicTopology = ".";
        v.cladeRank = ".";
        v.phylum = ".";
        out.push_back(std::move(v));
    }

    return out;
}

static std::vector<VariantCallBridge>
mem_chain_sv_calls(const std::string& qAsm,
                   const std::unordered_map<std::string, std::string>& contigs,
                   const SimpleRefIndex& refIdx,
                   RefSearchCache& cache,
                   SingleRefMemCache& memCache,
                   const Options& o,
                   const tol::FederatedOptions& eff_fo,
                   query_input::QueryMode mode) {
    std::vector<VariantCallBridge> out;
    if (cache.flatRefs.empty()) return out;

    out.reserve(contigs.size());
    for (const auto& kv : contigs) {
        const std::string& contigName = kv.first;
        const std::string& seq = kv.second;
        if (static_cast<int>(seq.size()) < o.minSvLen) continue;
        if (is_low_complexity_sequence(seq)) continue;

        const int prefilterK = reads_mode_overlap_k(mode, seq.size());
        const auto prefilterHashes = cache.query_hashes(seq, prefilterK);
        cache.ensure_hashes(prefilterK);
        double prefilterBestOverlap = 0.0;
        const RefContigInfo* prefilterBestRef = nullptr;
        int prefilterBestRefLen = 0;
        for (size_t i = 0; i < cache.flatRefs.size(); ++i) {
            const double frac = cache.containment_with_ref_hashes(prefilterHashes, i, prefilterK);
            if (frac > prefilterBestOverlap) {
                prefilterBestOverlap = frac;
                prefilterBestRef = cache.flatRefs[i].info;
                prefilterBestRefLen = cache.flatRefs[i].info->length;
            }
        }
        if (mode != query_input::QueryMode::ASSEMBLY &&
            prefilterBestOverlap >= 0.98 &&
            is_reads_mode_fragment(seq.size(), prefilterBestRefLen, prefilterBestOverlap, o, mode)) {
            continue;
        }
        VariantCallBridge cachedChainCall;
        if (prefilterBestRef != nullptr &&
            prefilterBestOverlap >= 0.01 &&
            try_mem_chain_call_single_ref_cached(
                qAsm, contigName, seq, prefilterBestRef, memCache, eff_fo, cachedChainCall)) {
            out.push_back(std::move(cachedChainCall));
            continue;
        }

        const auto shortlist = select_mem_chain_refs(
            contigName, seq, refIdx, cache, mode,
            (mode == query_input::QueryMode::ASSEMBLY) ? 8u : 12u);
        if (shortlist.empty()) continue;
        const auto refBundle = make_refseq_bundle(shortlist);

        VariantCallBridge chainCall;
        if (tol::try_mem_chain_call_public(qAsm, contigName, seq, refBundle.ptrs, eff_fo, chainCall)) {
            out.push_back(std::move(chainCall));
        }
    }
    return out;
}

static bool try_mem_chain_call_single_ref_cached(
        const std::string& qAsm,
        const std::string& qContig,
        const std::string& qSeq,
        const RefContigInfo* refInfo,
        SingleRefMemCache& memCache,
        const tol::FederatedOptions& fo,
        VariantCallBridge& call) {
    if (refInfo == nullptr || refInfo->sequence.empty() || qSeq.empty()) return false;

    SingleRefMemIndex& idx = memCache.get(refInfo);
    const tol::SuffixArray& sa = idx.sa;
    const tol::TolGlobal::RefSeq& primaryRef = idx.ref;
    const std::string rcSeq = tol::SuffixArray::revcomp(qSeq);

    auto min_mem_from_k = [](int k) {
        return std::max(15, k - 5);
    };
    struct ChainAttempt {
        VariantCallBridge call;
        double score = 0.0;
        int anchors = 0;
        bool valid = false;
    };
    auto attempt_chain = [&](int minMem, bool secondaryPass) {
        ChainAttempt out;
        auto fwdMems = sa.find_mems(qSeq, minMem);
        auto revMems = sa.find_mems(rcSeq, minMem);
        for (auto& m : revMems)
            m.qPos = static_cast<int>(qSeq.size()) - m.qPos - m.len;

        std::vector<tol::SuffixArray::Mem> allMems;
        std::vector<bool> isRev;
        allMems.reserve(fwdMems.size() + revMems.size());
        isRev.reserve(fwdMems.size() + revMems.size());
        for (auto& m : fwdMems) { allMems.push_back(m); isRev.push_back(false); }
        for (auto& m : revMems) { allMems.push_back(m); isRev.push_back(true); }
        if (allMems.empty()) return out;

        std::vector<int> order(allMems.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(),
                  [&](int a, int b) {
                      return allMems[static_cast<size_t>(a)].qPos <
                             allMems[static_cast<size_t>(b)].qPos;
                  });

        tol::ChainTreap treap;
        const int maxGap = fo.chainGapBand > 0 ? fo.chainGapBand : 5000;
        for (int i : order) {
            const auto& m = allMems[static_cast<size_t>(i)];
            treap.insert_and_chain(m.qPos, m.rPos, m.len,
                                   static_cast<float>(m.len), maxGap);
        }

        auto chainIdx = treap.best_chain_path();
        if (chainIdx.empty()) return out;

        std::unordered_map<uint64_t, size_t> posToMemIdx;
        posToMemIdx.reserve(allMems.size());
        for (size_t mi = 0; mi < allMems.size(); ++mi) {
            uint64_t key = (static_cast<uint64_t>(allMems[mi].qPos) << 32) |
                           static_cast<uint64_t>(static_cast<uint32_t>(allMems[mi].rPos));
            posToMemIdx.emplace(key, mi);
        }

        std::vector<tol::SuffixArray::Mem> chain;
        std::vector<bool> chainRev;
        chain.reserve(chainIdx.size());
        chainRev.reserve(chainIdx.size());
        for (int ni : chainIdx) {
            const auto& nd = treap.nodes_[static_cast<size_t>(ni)];
            uint64_t key = (static_cast<uint64_t>(nd.qPos) << 32) |
                           static_cast<uint64_t>(static_cast<uint32_t>(nd.rPos));
            auto it = posToMemIdx.find(key);
            if (it != posToMemIdx.end()) {
                chain.push_back(allMems[it->second]);
                chainRev.push_back(isRev[it->second]);
            }
        }
        if (chain.empty()) return out;

        const double bestScore = static_cast<double>(treap.best_chain_score());
        if (chain.size() < 2 || bestScore < fo.minBlockScore)
            return out;

        auto res = tol::SvTypeFromChain::classify(chain, chainRev, sa, fo.minSvLen);
        if (res.type == tol::SvTypeFromChain::Type::NONE) return out;

        out.call.qAsm = qAsm;
        out.call.qContig = qContig;
        out.call.refAsm = primaryRef.asmName.empty() ? "unknown" : primaryRef.asmName;
        out.call.refContig = primaryRef.contig.empty() ? "." : primaryRef.contig;
        out.call.refPos = 0;
        out.call.refEnd = 0;
        out.call.pos = std::max(1, res.qBreakStart + 1);
        out.call.end = std::max(out.call.pos, res.qBreakEnd > 0 ? res.qBreakEnd : out.call.pos);
        out.call.svlen = res.svLen;
        out.call.genotype = "0/1";
        out.call.gq = 40.0;
        out.call.blockScore = bestScore;
        out.call.anchors = static_cast<int>(chain.size());
        out.call.alignmentMode = secondaryPass
            ? "mem_chain_cached_single_ref;secondary_seed_rescue"
            : "mem_chain_cached_single_ref";
        out.call.mapq = 50.0;
        out.call.annotation = "NONE";
        out.call.triallelicTopology = ".";
        out.call.isNonRefVariant = false;
        out.call.cladeRank = ".";
        out.call.phylum = ".";

        using T = tol::SvTypeFromChain::Type;
        switch (res.type) {
            case T::INS:
                out.call.type = "INS";
                out.call.pantreeClass = "INS";
                out.call.refPos = res.rBreakStart > 0 ? (res.rBreakStart + 1) : 0;
                out.call.refEnd = out.call.refPos;
                break;
            case T::DEL:
                out.call.type = "DEL";
                out.call.pantreeClass = "DEL";
                out.call.refPos = res.rBreakStart > 0 ? (res.rBreakStart + 1) : 0;
                out.call.refEnd = res.rBreakEnd > 0 ? res.rBreakEnd : out.call.refPos;
                break;
            case T::INV:
                out.call.type = "INV";
                out.call.pantreeClass = "INV";
                out.call.refPos = res.rBreakStart > 0 ? (res.rBreakStart + 1) : 0;
                out.call.refEnd = res.rBreakEnd > 0 ? res.rBreakEnd : out.call.refPos;
                break;
            case T::DUP:
                out.call.type = "DUP";
                out.call.pantreeClass = "DUP";
                out.call.refPos = res.rBreakStart > 0 ? (res.rBreakStart + 1) : 0;
                out.call.refEnd = res.rBreakEnd > 0 ? res.rBreakEnd : out.call.refPos;
                break;
            case T::TRA:
                out.call.type = "TRA";
                out.call.pantreeClass = "NON_REF";
                out.call.refPos = !chain.empty() ? (chain.front().rPos + chain.front().len + 1) : 0;
                out.call.refEnd = out.call.refPos;
                out.call.mateContig = res.rContig.empty() ? primaryRef.contig : res.rContig;
                out.call.matePos = res.rBreakStart + 1;
                out.call.mateEnd = res.rBreakEnd > 0 ? res.rBreakEnd : out.call.matePos;
                out.call.mateRefAsm = primaryRef.asmName.empty() ? "." : primaryRef.asmName;
                out.call.mateOffReference = false;
                break;
            default:
                return out;
        }
        out.score = bestScore;
        out.anchors = static_cast<int>(chain.size());
        out.valid = true;
        return out;
    };

    const int primaryMinMem = min_mem_from_k(fo.primarySketchParams.k);
    ChainAttempt best = attempt_chain(primaryMinMem, false);
    if (fo.useSecondarySeeds) {
        const int secondaryMinMem = min_mem_from_k(fo.secondarySketchParams.k);
        const bool rescueRequested = !best.valid ||
            best.anchors < static_cast<int>(std::max<size_t>(fo.repeatRescueMinAnchors, 2));
        if (secondaryMinMem < primaryMinMem && rescueRequested) {
            ChainAttempt rescue = attempt_chain(secondaryMinMem, true);
            if (rescue.valid &&
                (!best.valid || rescue.anchors > best.anchors || rescue.score > best.score)) {
                best = std::move(rescue);
            }
        }
    }
    if (!best.valid) return false;
    call = std::move(best.call);
    return true;
}

static std::vector<VariantCallBridge>
simple_offref_fallback_calls(const std::string& qAsm,
                             const std::unordered_map<std::string, std::string>& contigs,
                             const SimpleRefIndex& refIdx,
                             RefSearchCache& cache,
                             const Options& o,
                             query_input::QueryMode mode) {
    std::vector<VariantCallBridge> out;
    const int noveltyK = std::max(5, std::min(o.tolFallbackK, 9));
    cache.ensure_hashes(noveltyK);
    for (const auto& kv : contigs) {
        const std::string& contigName = kv.first;
        const std::string& seq = kv.second;
        if (static_cast<int>(seq.size()) < o.minSvLen) continue;
        if (is_low_complexity_sequence(seq)) continue;

        // Do not reject a contig purely because its FASTA name exists in the
        // reference set. In the hint-free simulator and in real data, a query
        // contig can reuse a biological name (e.g. "ctg1") while still being
        // sequence-novel relative to every indexed reference contig with that
        // name. The length fallback above already had a chance to explain same-
        // named contigs as INS/DEL from length alone; reaching this function
        // means there was no such call, so novelty should now be decided from
        // sequence content.
        const RefContigInfo* bestNamed =
            best_ref_match(refIdx, contigName, static_cast<int>(seq.size()));

        double bestOverlap = 0.0;
        double bestLocusOverlap = 0.0;
        std::string bestAsm = "OFF_REFERENCE";
        std::string bestContig = ".";
        int bestRefLen = 0;
        const int locusK = reads_mode_overlap_k(mode, seq.size());
        const auto noveltyHashes = cache.query_hashes(seq, noveltyK);
        const auto locusHashes = cache.query_hashes(seq, locusK);
        cache.ensure_hashes(locusK);
        const bool useContainment = (mode != query_input::QueryMode::ASSEMBLY);

        if (bestNamed != nullptr && !bestNamed->sequence.empty()) {
            bestOverlap = useContainment
                ? cache.containment_with_ref_hashes(noveltyHashes, bestNamed, noveltyK)
                : cache.overlap_with_ref_hashes(noveltyHashes, bestNamed, noveltyK);
            bestLocusOverlap = useContainment
                ? cache.containment_with_ref_hashes(locusHashes, bestNamed, locusK)
                : cache.overlap_with_ref_hashes(locusHashes, bestNamed, locusK);
            bestAsm = bestNamed->asmName;
            bestContig = bestNamed->contigName;
            bestRefLen = bestNamed->length;
        }

        for (size_t i = 0; i < cache.flatRefs.size(); ++i) {
            const RefContigInfo& ref = *cache.flatRefs[i].info;
            const double frac = useContainment
                ? cache.containment_with_ref_hashes(noveltyHashes, i, noveltyK)
                : cache.overlap_with_ref_hashes(noveltyHashes, i, noveltyK);
            const double locusFrac = useContainment
                ? cache.containment_with_ref_hashes(locusHashes, i, locusK)
                : cache.overlap_with_ref_hashes(locusHashes, i, locusK);
            if (frac > bestOverlap) {
                bestOverlap = frac;
                bestAsm = ref.asmName;
                bestContig = ref.contigName;
                bestRefLen = ref.length;
            }
            if (locusFrac > bestLocusOverlap)
                bestLocusOverlap = locusFrac;
        }

        if (is_reads_mode_fragment(seq.size(), bestRefLen, bestOverlap, o, mode))
            continue;

        if (mode != query_input::QueryMode::ASSEMBLY) {
            const int delta = static_cast<int>(seq.size()) - bestRefLen;
            const bool locusScaleIndel =
                bestRefLen > 0 &&
                std::abs(delta) >= o.minSvLen &&
                std::abs(delta) <= o.maxSvLen;
            if (bestLocusOverlap >= 0.12) continue;
            if (locusScaleIndel && bestLocusOverlap >= 0.04) continue;
        }

        const std::string tier = novelty_tier_for_overlap(bestOverlap);
        if (tier != "NOVEL" && tier != "NOVEL_WEAK" && tier != "DIVERGED") continue;

        // Compute background GC from the best-matching reference if available,
        // otherwise default to 0.45.
        double cladeGc = 0.45;
        if (bestNamed != nullptr && !bestNamed->sequence.empty()) {
            size_t gcCount = 0;
            for (char c : bestNamed->sequence)
                if (c == 'G' || c == 'C' || c == 'g' || c == 'c') ++gcCount;
            cladeGc = bestNamed->sequence.empty() ? 0.45
                : static_cast<double>(gcCount) / static_cast<double>(bestNamed->sequence.size());
        }
        const tol::ElementClass ec = seq.empty()
            ? tol::ElementClass::NONE
            : tol::classify_repeat_element(std::string_view(seq.data(), seq.size()), cladeGc);

        VariantCallBridge v;
        v.qAsm = qAsm;
        v.qContig = contigName;
        v.refAsm = bestAsm;
        v.refContig = bestContig;
        v.refPos = 0;
        v.refEnd = 0;
        v.type = "OFF_REF";
        v.svlen = static_cast<int>(seq.size());
        v.pos = 1;
        v.end = std::max(1, static_cast<int>(seq.size()));
        v.annotation = tier;
        v.genotype = "0/1";
        v.gq = 20.0;
        v.alignmentMode = "simple_offref_fallback";
        v.mapq = 10.0;
        v.blockScore = 8.0;
        v.anchors = 0;
        v.pantreeClass = "NON_REF";
        v.isNonRefVariant = true;
        v.triallelicTopology = ".";
        v.elementClass = tol::element_class_name(ec);
        v.cladeRank = ".";
        v.phylum = ".";
        out.push_back(std::move(v));
    }
    if (mode != query_input::QueryMode::ASSEMBLY && out.size() > 1) {
        std::vector<VariantCallBridge> filtered;
        filtered.reserve(out.size());
        const VariantCallBridge* bestWeak = nullptr;
        for (const auto& call : out) {
            const bool weakGeneric =
                call.type == "OFF_REF" &&
                call.annotation == "NOVEL_WEAK" &&
                (call.elementClass.empty() || call.elementClass == "NONE");
            if (!weakGeneric) {
                filtered.push_back(call);
                continue;
            }
            if (bestWeak == nullptr ||
                call.svlen > bestWeak->svlen ||
                (call.svlen == bestWeak->svlen && call.qContig < bestWeak->qContig)) {
                bestWeak = &call;
            }
        }
        if (bestWeak != nullptr) filtered.push_back(*bestWeak);
        return filtered;
    }
    return out;
}

// reads_mode_sv_calls
//
// For LONG_READS / SHORT_READS mode only.
//
// Pseudo-contig names (lr_pc0, sr_unitig3…) have no name correspondence with
// reference contig names, so the existing name-keyed simple_length_fallback
// and the hierarchical Path B (byContig.find) both miss them.
//
// This function bridges that gap in three passes, from most- to least-precise:
//
// Pass 1 — MEM chain (try_mem_chain_call):
//   For each pseudo-contig, scan all reference contigs and try to build a
//   MEM chain.  This can detect INS, DEL, INV, DUP, TRA from the chain
//   geometry, exactly as the assembly-mode hierarchical engine does.
//   A chain with ≥ minAnchors MEMs and blockScore ≥ minBlockScore is accepted.
//
// Pass 2 — k-mer overlap + length delta (INS/DEL only):
//   If the MEM chain fails or the pseudo-contig is too short for MEMs, find the
//   best-matching reference contig by k-mer Jaccard overlap and call INS/DEL
//   from the length delta.  Overlap ≥ 0.05 is required to avoid random matches.
//   k is mode-aware: long-reads uses k=13 (tolerates ~10% noise in consensus),
//   short-reads uses k=17 (cleaner unitigs).
//
// Pass 3 — OFF_REF:
//   Pseudo-contigs with overlap < 0.20 are emitted as OFF_REF; they are also
//   caught downstream by simple_offref_fallback_calls, but emitting them here
//   ensures they get queryMode stamped correctly.
//
// Position estimate for Pass 2:
//   Estimated as max(1, min(qlen, refLen)/2) — symmetric midpoint — since we
//   have no alignment to locate the breakpoint precisely.
static std::vector<VariantCallBridge>
reads_mode_sv_calls(const std::string& qAsm,
                    const std::unordered_map<std::string, std::string>& contigs,
                    const SimpleRefIndex& refIdx,
                    RefSearchCache& cache,
                    SingleRefMemCache& memCache,
                    const Options& o,
                    const tol::FederatedOptions& eff_fo,
                    query_input::QueryMode mode) {
    std::vector<VariantCallBridge> out;
    out.reserve(contigs.size());

    // Mode-aware k for k-mer overlap estimation (Pass 2).
    // long-reads consensus has ~5-10% noise → shorter k for robustness.
    // short-reads unitigs are cleaner → longer k for specificity.
    // Collect all reference contigs into a flat vector for fast iteration.
    // Also build a borrowed set of RefSeq pointers compatible with try_mem_chain_call.
    if (cache.flatRefs.empty()) return out;

    for (const auto& kv : contigs) {
        const std::string& contigName = kv.first;
        const std::string& seq        = kv.second;
        if (static_cast<int>(seq.size()) < o.minSvLen) continue;
        if (is_low_complexity_sequence(seq)) continue;

        // ── Pass 1: MEM chain ──────────────────────────────────────────
        // try_mem_chain_call builds a SA over all refPtrs and searches for MEMs.
        // This is the most precise path: it can detect INS/DEL/INV/DUP/TRA.
        const int prefilterK = reads_mode_overlap_k(mode, seq.size());
        const auto prefilterHashes = cache.query_hashes(seq, prefilterK);
        cache.ensure_hashes(prefilterK);
        double prefilterBestOverlap = 0.0;
        const RefContigInfo* prefilterBestRef = nullptr;
        int prefilterBestRefLen = 0;
        for (size_t i = 0; i < cache.flatRefs.size(); ++i) {
            const double frac = cache.containment_with_ref_hashes(prefilterHashes, i, prefilterK);
            if (frac > prefilterBestOverlap) {
                prefilterBestOverlap = frac;
                prefilterBestRef = cache.flatRefs[i].info;
                prefilterBestRefLen = cache.flatRefs[i].info->length;
            }
        }
        if (mode != query_input::QueryMode::ASSEMBLY &&
            prefilterBestOverlap >= 0.98 &&
            is_reads_mode_fragment(seq.size(), prefilterBestRefLen, prefilterBestOverlap, o, mode)) {
            continue;
        }
        // Pass 1a: single-ref MEM chain (fast path; good for INS/DEL).
        VariantCallBridge cachedChainCall;
        bool singleRefOk = false;
        if (prefilterBestRef != nullptr &&
            prefilterBestOverlap >= 0.01 &&
            try_mem_chain_call_single_ref_cached(
                qAsm, contigName, seq, prefilterBestRef, memCache, eff_fo, cachedChainCall)) {
            out.push_back(cachedChainCall);
            singleRefOk = true;
            // TRA/INV/DUP already confirmed from a single reference — no need for
            // the more expensive multi-ref chain.
            if (cachedChainCall.type != "INS" && cachedChainCall.type != "DEL") {
                continue;
            }
        }

        // Pass 1b: multi-ref MEM chain (detects TRA/INV/DUP across contig boundaries).
        // Runs even when Pass 1a succeeded with an indel so that a cross-contig TRA
        // gets added as a competing candidate and arbitration can pick the better call.
        const auto shortlist = select_mem_chain_refs(contigName, seq, refIdx, cache, mode, 12u);
        const auto refBundle = make_refseq_bundle(shortlist);
        VariantCallBridge chainCall;
        if (tol::try_mem_chain_call_public(qAsm, contigName, seq, refBundle.ptrs, eff_fo, chainCall)) {
            out.push_back(std::move(chainCall));
            continue;
        }
        if (singleRefOk) continue;

        // ── Pass 2: k-mer overlap + length delta (INS/DEL) ────────────
        double    bestOverlap = 0.0;
        int       bestRefLen  = 0;
        std::string bestAsm   = "OFF_REFERENCE";
        std::string bestContig= ".";
        const int overlapK = reads_mode_overlap_k(mode, seq.size());

        const auto qHashes = cache.query_hashes(seq, overlapK);
        cache.ensure_hashes(overlapK);

        for (size_t i = 0; i < cache.flatRefs.size(); ++i) {
            const auto& fr = cache.flatRefs[i];
            const double frac = cache.containment_with_ref_hashes(qHashes, i, overlapK);
            if (frac > bestOverlap) {
                bestOverlap = frac;
                bestRefLen  = fr.info->length;
                bestAsm     = fr.info->asmName;
                bestContig  = fr.info->contigName;
            }
        }

        const int qlen  = static_cast<int>(seq.size());
        const int delta = qlen - bestRefLen;

        const bool locusSized =
            !is_reads_mode_fragment(seq.size(), bestRefLen, bestOverlap, o, mode);

        if (locusSized &&
            bestOverlap >= 0.05 && std::abs(delta) >= o.minSvLen
                                 && std::abs(delta) <= o.maxSvLen) {
            VariantCallBridge v;
            v.qAsm          = qAsm;
            v.qContig       = contigName;
            v.refAsm        = bestAsm;
            v.refContig     = bestContig;
            v.annotation    = "NONE";
            v.genotype      = "0/1";
            v.gq            = 10.0 + bestOverlap * 40.0;
            v.alignmentMode = "reads_mode_kmer_fallback";
            v.mapq          = 20.0 + bestOverlap * 10.0;
            v.blockScore    = 10.0 + bestOverlap * 10.0;
            v.anchors       = std::max(1, o.minAnchors);
            v.type          = (delta < 0) ? "DEL" : "INS";
            v.svlen         = delta;
            // Symmetric midpoint position estimate — no alignment available.
            const int midPoint = std::max(1, std::min(qlen, bestRefLen) / 2);
            v.pos           = midPoint;
            v.end           = (v.type == "DEL")
                ? std::min(bestRefLen, midPoint + std::abs(delta) - 1)
                : midPoint;
            v.refPos        = v.pos;
            v.refEnd        = v.end;
            v.pantreeClass  = v.type;
            v.isNonRefVariant    = false;
            v.triallelicTopology = ".";
            v.cladeRank = ".";
            v.phylum = ".";
            out.push_back(std::move(v));
            continue;
        }

        // ── Pass 3: OFF_REF (overlap < 0.20) ──────────────────────────
        // Handled by simple_offref_fallback_calls downstream — no duplicate emission.
    }
    return out;
}


static const RefContigInfo* find_ref_by_asm_and_contig(const SimpleRefIndex& refIdx,
                                                       const std::string& asmName,
                                                       const std::string& contigName) {
    auto it = refIdx.find(strip_sv_suffix(contigName));
    if (it == refIdx.end()) return nullptr;
    const RefContigInfo* contigOnly = nullptr;
    for (const auto& cand : it->second) {
        if (cand.asmName == asmName && cand.contigName == contigName) return &cand;
        if (cand.contigName == contigName && contigOnly == nullptr) contigOnly = &cand;
    }
    // Some callers carry routed clade labels in refAsm rather than the concrete
    // assembly stem. Fall back to the contig match so candidate arbitration can
    // still recover sequence overlap instead of treating the call as off-ref.
    return contigOnly;
}

static double candidate_ref_overlap(const std::string& qSeq,
                                    const VariantCallBridge& call,
                                    const SimpleRefIndex& refIdx,
                                    RefSearchCache& cache,
                                    query_input::QueryMode mode) {
    if (call.refAsm.empty() || call.refAsm == "OFF_REFERENCE" || call.refContig.empty() || call.refContig == ".") {
        return 0.0;
    }
    const RefContigInfo* ref = find_ref_by_asm_and_contig(refIdx, call.refAsm, call.refContig);
    if (ref == nullptr || ref->sequence.empty()) return 0.0;
    const int k = reads_mode_overlap_k(mode, qSeq.size());
    const auto qHashes = cache.query_hashes(qSeq, k);
    return (mode == query_input::QueryMode::ASSEMBLY)
        ? cache.overlap_with_ref_hashes(qHashes, ref, k)
        : cache.containment_with_ref_hashes(qHashes, ref, k);
}

static query_input::QueryMode call_query_mode(const VariantCallBridge& call,
                                              query_input::QueryMode fallback) {
    if (!call.queryMode.empty()) {
        try {
            return query_input::parse_mode(call.queryMode);
        } catch (...) {
        }
    }
    return fallback;
}

static tol::EvidenceLayer evidence_layer_for_mode(query_input::QueryMode mode) {
    switch (mode) {
        case query_input::QueryMode::ASSEMBLY:    return tol::EvidenceLayer::ASSEMBLY;
        case query_input::QueryMode::LONG_READS:  return tol::EvidenceLayer::LONG_READ;
        case query_input::QueryMode::SHORT_READS: return tol::EvidenceLayer::SHORT_READ;
    }
    return tol::EvidenceLayer::ASSEMBLY;
}

static bool candidate_supports_same_event(const VariantCallBridge& a,
                                          const VariantCallBridge& b,
                                          const Options& o) {
    if (a.type != b.type) return false;
    if (a.qContig != b.qContig) return false;
    if (a.type == "OFF_REF") {
        const int slack = std::max(50, o.minSvLen);
        return std::abs(a.pos - b.pos) <= slack &&
               std::abs(a.end - b.end) <= std::max(slack, o.minSvLen * 2);
    }

    if (!a.refContig.empty() && !b.refContig.empty() &&
        a.refContig != "." && b.refContig != "." &&
        a.refContig != b.refContig) {
        return false;
    }
    if (!a.refAsm.empty() && !b.refAsm.empty() &&
        a.refAsm != "OFF_REFERENCE" && b.refAsm != "OFF_REFERENCE" &&
        a.refAsm != b.refAsm) {
        return false;
    }

    const int slack = std::max(50, o.minSvLen);
    if (std::abs(a.pos - b.pos) > slack &&
        std::abs(a.end - b.end) > std::max(slack, o.minSvLen * 2)) {
        return false;
    }

    const int lenA = std::abs(a.svlen);
    const int lenB = std::abs(b.svlen);
    if (lenA > 0 && lenB > 0 &&
        std::abs(lenA - lenB) > std::max(25, std::max(lenA, lenB) / 4)) {
        return false;
    }

    if (is_translocation_type(a.type) &&
        ((!a.mateContig.empty() && !b.mateContig.empty() && a.mateContig != b.mateContig) ||
         (a.matePos > 0 && b.matePos > 0 &&
          std::abs(a.matePos - b.matePos) > std::max(100, o.minSvLen * 2)))) {
        return false;
    }
    return true;
}

static tol::EvidenceObservation make_evidence_observation(
        const VariantCallBridge& call,
        double overlap,
        const query_input::EvidenceFusionRecommendation& rec,
        query_input::QueryMode mode) {
    tol::EvidenceObservation obs;
    obs.layer = evidence_layer_for_mode(mode);
    obs.depth = std::max(1.0, static_cast<double>(std::max(0, call.anchors)));
    obs.mapq = std::max(0.0, call.mapq);
    obs.breakpointSupport = std::max(
        0.0,
        std::min(1.0, 0.08 * static_cast<double>(std::max(0, call.anchors)) +
                          (call.alignmentMode.find("mem_chain") != std::string::npos ? 0.25 : 0.0) +
                          (is_translocation_type(call.type) || call.type == "INV" ? 0.10 : 0.0)));
    obs.spanSupport = std::max(
        0.0,
        std::min(1.0, overlap +
                          (call.alignmentMode.find("simple_length_fallback") != std::string::npos ? 0.20 : 0.0) +
                          (call.type == "OFF_REF" ? 0.15 : 0.0)));

    const double modeWeight =
        (obs.layer == tol::EvidenceLayer::ASSEMBLY)
            ? rec.expectedAssemblyWeight
            : rec.expectedReadWeight;
    obs.depth *= std::max(0.25, modeWeight);

    double alt = std::log1p(std::max(0.0, call.blockScore))
               + 0.18 * static_cast<double>(std::max(0, call.anchors))
               + 0.03 * std::max(0.0, call.mapq)
               + 0.03 * std::max(0.0, call.gq)
               + 2.5 * overlap;
    double ref = 0.85;
    if (call.type == "OFF_REF") {
        alt += (call.annotation == "NOVEL" || call.annotation == "NOVEL_WEAK") ? 0.60 : 0.25;
        ref += overlap * 2.2;
    } else {
        if (overlap < 0.02) ref += 1.5;
        if ((call.type == "INV" || call.type == "TRA" || call.type == "DUP") && overlap < 0.05)
            ref += 0.75;
        if ((call.type == "INV" || call.type == "TRA" || call.type == "DUP") && overlap < 0.10) {
            alt -= 2.5;
            ref += 5.0 + (0.10 - overlap) * 20.0;
            if (call.type == "TRA" && call.mateContig.empty()) {
                alt -= 1.0;
                ref += 2.0;
            }
        }
        if (call.alignmentMode.find("mem_chain") != std::string::npos) alt += 0.45;
        if (call.alignmentMode.find("simple_length_fallback") != std::string::npos) alt += 0.25;
    }
    obs.logLikelihoodAlt = alt;
    obs.logLikelihoodRef = ref;
    return obs;
}

static tol::FusedEvidenceScore fuse_candidate_evidence(
        const VariantCallBridge& call,
        const std::vector<VariantCallBridge>& peers,
        const std::string& qSeq,
        const SimpleRefIndex& refIdx,
        RefSearchCache& cache,
        const Options& o,
        query_input::QueryMode mode,
        const query_input::CoverageReport* report) {
    query_input::EvidenceFusionRecommendation rec;
    rec.enableProbabilisticFusion = true;
    rec.priorAlt = 0.50;
    rec.expectedAssemblyWeight = 1.0;
    rec.expectedReadWeight = 1.0;
    if (report != nullptr) {
        rec = query_input::recommend_evidence_fusion(
            mode, report->coverageTier, report->pseudoContigs);
    }

    std::vector<tol::EvidenceObservation> evidence;
    evidence.reserve(peers.size() + 1);
    for (const auto& peer : peers) {
        if (!candidate_supports_same_event(call, peer, o)) continue;
        const query_input::QueryMode peerMode = call_query_mode(peer, mode);
        const double peerOverlap = candidate_ref_overlap(qSeq, peer, refIdx, cache, peerMode);
        evidence.push_back(make_evidence_observation(peer, peerOverlap, rec, peerMode));
    }
    if (evidence.empty()) {
        const double overlap = candidate_ref_overlap(qSeq, call, refIdx, cache, mode);
        evidence.push_back(make_evidence_observation(call, overlap, rec, mode));
    }
    return tol::fuse_probabilistic_evidence(evidence, rec.priorAlt);
}

static double candidate_priority_score(const VariantCallBridge& call,
                                       const std::vector<VariantCallBridge>& peers,
                                       const std::string& qSeq,
                                       const SimpleRefIndex& refIdx,
                                       RefSearchCache& cache,
                                       const Options& o,
                                       query_input::QueryMode mode,
                                       const query_input::CoverageReport* report,
                                       tol::FusedEvidenceScore* fusedOut = nullptr) {
    const bool offRef = (call.type == "OFF_REF" || call.refAsm == "OFF_REFERENCE");
    const bool largeRearr = (call.type == "INV" || call.type == "TRA" || call.type == "DUP");
    const bool indel = (call.type == "INS" || call.type == "DEL");
    double score = offRef ? -250.0 : 40.0;

    const double overlap = candidate_ref_overlap(qSeq, call, refIdx, cache, mode);
    const tol::FusedEvidenceScore fused =
        fuse_candidate_evidence(call, peers, qSeq, refIdx, cache, o, mode, report);
    if (fusedOut != nullptr) *fusedOut = fused;
    score += overlap * 180.0;
    score += std::min(400.0, call.blockScore) * 0.20;
    score += std::min(80, call.anchors) * 2.0;
    score += std::min(60.0, call.mapq) * 0.30;
    score += std::min(60.0, call.gq) * 0.20;
    score += fused.posteriorAlt * 120.0;
    score += std::max(-6.0, std::min(12.0, fused.logOddsAlt)) * 12.0;
    score += std::min(8.0, fused.effectiveDepth) * 2.0;
    score += std::min<size_t>(4, fused.layersUsed) * 4.0;
    if (fused.posteriorAlt < 0.35) score -= 120.0;
    if (fused.supports_variant(0.90)) score += 30.0;

    if (call.alignmentMode.find("mem_chain") != std::string::npos) score += 8.0;
    if (call.alignmentMode.find("simple_length_fallback") != std::string::npos) score += 12.0;
    if (call.alignmentMode.find("reads_mode_kmer_fallback") != std::string::npos) score += 10.0;

    const RefContigInfo* namedBest = best_ref_match(refIdx, call.qContig, static_cast<int>(qSeq.size()));
    if (namedBest != nullptr) {
        const int delta = static_cast<int>(qSeq.size()) - namedBest->length;
        const int namedK = (mode == query_input::QueryMode::LONG_READS) ? 13
                         : (mode == query_input::QueryMode::SHORT_READS) ? 17
                         : 9;
        const double namedOverlap = namedBest->sequence.empty() ? 0.0
            : ((mode == query_input::QueryMode::ASSEMBLY)
                ? cache.overlap_with_ref_hashes(cache.query_hashes(qSeq, namedK), namedBest, namedK)
                : cache.containment_with_ref_hashes(cache.query_hashes(qSeq, namedK), namedBest, namedK));
        const bool sameNamedRef = (call.refAsm == namedBest->asmName && call.refContig == namedBest->contigName);
        const bool deltaLooksLikeIndel = std::abs(delta) >= o.minSvLen && std::abs(delta) <= o.maxSvLen;

        if (sameNamedRef) score += 35.0;

        if (deltaLooksLikeIndel && indel) {
            const int svAbs = std::abs(call.svlen);
            const int dAbs = std::abs(delta);
            if (std::abs(svAbs - dAbs) <= std::max(10, dAbs / 10)) score += 80.0;
            if (sameNamedRef) score += 40.0;
            score += namedOverlap * 60.0;
        }

        if (largeRearr && deltaLooksLikeIndel && namedOverlap >= 0.35) {
            score -= 90.0;
            if (!sameNamedRef) score -= 50.0;
        }
        if ((call.type == "INV" || call.type == "TRA") && std::abs(call.svlen) > static_cast<int>(qSeq.size() * 0.60)) {
            score -= 35.0;
        }
    } else if (mode != query_input::QueryMode::ASSEMBLY) {
        // Reads-mode pseudo-contigs usually have synthetic names (sr_unitig*,
        // lr_pc*), so name-based arbitration is unavailable. In that case,
        // favor indels that still show real sequence overlap and penalize
        // large rearrangements whose assigned reference has essentially none.
        if (indel && overlap >= 0.04) score += 110.0;
        // TRA confirmed across two reference contigs (mateContig set, primary overlap
        // high) is treated on equal footing with a high-confidence indel.
        if (call.type == "TRA" && !call.mateContig.empty() && overlap >= 0.30) score += 115.0;
        if (largeRearr && overlap < 0.02) score -= 180.0;
        if (largeRearr && overlap < 0.05) score -= 70.0;
        if (call.type == "TRA" && overlap < 0.05 && call.mateContig.empty()) score -= 120.0;
        if ((call.type == "INV" || call.type == "TRA") &&
            std::abs(call.svlen) > static_cast<int>(qSeq.size() * 0.60) &&
            overlap < 0.05) {
            score -= 55.0;
        }
    }

    return score;
}

static std::vector<VariantCallBridge>
select_best_call_per_contig(const std::unordered_map<std::string, std::string>& contigs,
                            const std::unordered_map<std::string, std::vector<VariantCallBridge>>& candidates,
                            const SimpleRefIndex& refIdx,
                            RefSearchCache& cache,
                            const Options& o,
                            query_input::QueryMode mode,
                            const std::string& modeLabel,
                            const query_input::CoverageReport* report = nullptr) {
    std::vector<VariantCallBridge> out;
    out.reserve(candidates.size());
    for (const auto& kv : candidates) {
        const auto seqIt = contigs.find(kv.first);
        if (seqIt == contigs.end()) continue;
        const std::string& qSeq = seqIt->second;
        const VariantCallBridge* best = nullptr;
        tol::FusedEvidenceScore bestFused;
        double bestScore = -1e100;
        for (const auto& cand : kv.second) {
            tol::FusedEvidenceScore fused;
            const double score = candidate_priority_score(
                cand, kv.second, qSeq, refIdx, cache, o, mode, report, &fused);
            if (best == nullptr || score > bestScore) {
                best = &cand;
                bestScore = score;
                bestFused = fused;
            }
        }
        if (best != nullptr) {
            VariantCallBridge chosen = *best;
            chosen.queryMode = modeLabel;
            chosen.fusedPosteriorAlt = bestFused.posteriorAlt;
            chosen.fusedLogOddsAlt = bestFused.logOddsAlt;
            chosen.fusedEffectiveDepth = bestFused.effectiveDepth;
            chosen.fusedLayersUsed = static_cast<int>(bestFused.layersUsed);
            out.push_back(std::move(chosen));
        }
    }
    return out;
}

static std::vector<VariantCallBridge>
deduplicate_read_mode_events(std::vector<VariantCallBridge> calls) {
    if (calls.size() <= 1) return calls;
    std::vector<VariantCallBridge> out;
    out.reserve(calls.size());
    std::unordered_map<std::string, size_t> bestByEvent;

    auto event_key = [](const VariantCallBridge& c) {
        if (c.type == "OFF_REF" && c.annotation == "NOVEL_WEAK")
            return std::string("OFF_REF:NOVEL_WEAK");
        if (c.type == "INS" || c.type == "DEL" || c.type == "DUP" ||
            c.type == "INV" || is_translocation_type(c.type)) {
            const int bucket = std::max(1, std::abs(c.svlen) / 100);
            return c.type + ":" + c.refAsm + ":" + c.refContig + ":" + std::to_string(bucket);
        }
        return c.type + ":" + c.qContig + ":" + std::to_string(c.pos);
    };

    auto support_score = [](const VariantCallBridge& c) {
        return c.blockScore + 2.0 * static_cast<double>(std::max(0, c.anchors)) +
               0.001 * static_cast<double>(std::abs(c.svlen)) +
               (c.annotation == "NOVEL_WEAK" ? 0.0 : 5.0);
    };

    for (auto& c : calls) {
        const std::string key = event_key(c);
        auto it = bestByEvent.find(key);
        if (it == bestByEvent.end()) {
            bestByEvent[key] = out.size();
            out.push_back(std::move(c));
            continue;
        }
        VariantCallBridge& prev = out[it->second];
        if (support_score(c) > support_score(prev))
            prev = std::move(c);
    }
    return out;
}

static int parse_pseudo_contig_read_support(const std::string& qContig) {
    const size_t tag = qContig.rfind("_n");
    if (tag == std::string::npos) return 1;
    size_t i = tag + 2;
    if (i >= qContig.size() || qContig[i] < '0' || qContig[i] > '9') return 1;
    long long value = 0;
    while (i < qContig.size() && qContig[i] >= '0' && qContig[i] <= '9') {
        value = value * 10 + (qContig[i] - '0');
        if (value > std::numeric_limits<int>::max()) return std::numeric_limits<int>::max();
        ++i;
    }
    return std::max(1, static_cast<int>(value));
}

static std::vector<VariantCallBridge>
filter_low_coverage_read_artifacts(std::vector<VariantCallBridge> calls,
                                   const query_input::CoverageReport& report) {
    if (report.coverageTier != query_input::CoverageTier::LOW || calls.empty())
        return calls;
    std::vector<VariantCallBridge> out;
    out.reserve(calls.size());
    for (auto& c : calls) {
        const bool rearrangement =
            c.type == "INV" || c.type == "DUP" || is_translocation_type(c.type);
        const bool anchoredIndel = c.type == "INS" || c.type == "DEL";
        const bool cachedSingleRef =
            c.alignmentMode.find("mem_chain_cached_single_ref") != std::string::npos;
        if (rearrangement &&
            cachedSingleRef) {
            continue;
        }
        if (report.mode == query_input::QueryMode::LONG_READS &&
            anchoredIndel && cachedSingleRef &&
            parse_pseudo_contig_read_support(c.qContig) < 3) {
            continue;
        }
        out.push_back(std::move(c));
    }
    return out;
}

// ============================================================
// Build FederatedOptions from parsed CLI options
// ============================================================
static tol::FederatedOptions make_fed_opts(const Options& o) {
    tol::SyncmerParams sp;
    sp.k               = o.k;
    sp.s               = o.s;
    sp.t               = o.t;
    sp.stride          = static_cast<size_t>(o.seedStride);
    sp.useIntervalHash = o.useIntervalHash;
    sp.ihWing          = o.ihWing;
    sp.ihMaxDist       = o.ihMaxDist;
    sp.ihResolution    = o.ihResolution;

    tol::SyncmerParams fb;
    fb.k               = o.tolFallbackK;
    fb.s               = o.tolFallbackS;
    fb.t               = o.t;
    fb.stride          = static_cast<size_t>(o.seedStride);
    fb.useIntervalHash = false;

    tol::SyncmerParams sec;
    sec.k               = o.secondaryK;
    sec.s               = o.secondaryS;
    sec.t               = o.secondaryT;
    sec.stride          = static_cast<size_t>(std::max(1, o.secondarySeedStride));
    sec.useIntervalHash = false;

    tol::FederatedOptions fo = tol::make_federated_opts(
        sp,
        fb,
        sec,
        o.tolRoutingDensity,
        static_cast<size_t>(o.routingTopN),
        o.minSvLen,
        o.maxSvLen,
        o.minBlockScore,
        static_cast<size_t>(o.minAnchors),
        o.chainGapBand,
        o.useSecondarySeeds,
        static_cast<size_t>(std::max(1, o.repeatRescueMinAnchors)),
        !o.quiet,
        static_cast<size_t>(o.threads),
        o.tolBaseGraphBuild,
        o.graphNativeMode,
        static_cast<size_t>(std::max(1, o.tolMinBlockBp)),
        static_cast<size_t>(std::max(1, o.tolMinChainAnchors)),
        o.tolMaxCladeGenomes,
        o.tolQueryWindowBp,
        o.tolQueryWindowOverlap,
        o.tolAncestralRecomb,
        o.tolRecombMinSegBp,
        o.tolRecombMaxBp);
    // Multi-ref SA cap = 2× single-ref limit; each of the up to 8 shortlisted
    // refs can be up to saMaxContigMB, so the total text cap is twice that.
    fo.saMaxTextMB = static_cast<size_t>(std::max(0, o.saMaxContigMB)) * 2;
    return fo;
}

static int normalized_svlen_for_output(const VariantCallBridge& v) {
    if (v.type == "DEL") return -std::abs(v.svlen);
    if (v.type == "OFF_REF") return std::max(1, std::abs(v.svlen));
    return v.svlen;
}

static int normalized_end_for_output(const VariantCallBridge& v) {
    const int pos = std::max(1, v.pos);
    const int abs_len = std::max(1, std::abs(normalized_svlen_for_output(v)));
    if (v.type == "INS") return pos;
    if (v.type == "OFF_REF") return pos + abs_len - 1;
    if (v.type == "DUP" || v.type == "INV") return pos + abs_len - 1;
    if (v.type == "DEL") return pos + abs_len - 1;
    if (v.end >= pos) return v.end;
    return pos + abs_len - 1;
}

// ============================================================
// VCF header writer
// ============================================================
static void write_vcf_header(std::ostream& out,
                              const std::string& source) {
    out << "##fileformat=VCFv4.3\n"
        << "##source=" << source << "\n"
        << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">\n"
        << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length\">\n"
        << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End pos\">\n"
        << "##INFO=<ID=ANNOT,Number=1,Type=String,Description=\"Annotation\">\n"
        << "##INFO=<ID=CLADE,Number=1,Type=String,Description=\"Routed clade\">\n"
        << "##INFO=<ID=REFCONTIG,Number=1,Type=String,Description=\"Reference contig for the anchored event when available\">\n"
        << "##INFO=<ID=REFPOS,Number=1,Type=Integer,Description=\"Reference-space start coordinate for the anchored event when available\">\n"
        << "##INFO=<ID=REFEND,Number=1,Type=Integer,Description=\"Reference-space end coordinate for the anchored event when available\">\n"
        << "##INFO=<ID=PHYLUM,Number=1,Type=String,Description=\"Phylum\">\n"
        << "##INFO=<ID=CLADE_RANK,Number=1,Type=String,Description=\"Routed clade rank\">\n"
        << "##INFO=<ID=BSCORE,Number=1,Type=Float,Description=\"Block score\">\n"
        << "##INFO=<ID=MAPQ,Number=1,Type=Float,Description=\"MAPQ\">\n"
        << "##INFO=<ID=OFFREF,Number=0,Type=Flag,Description=\"Off-reference novel sequence\">\n"
        << "##INFO=<ID=VT,Number=1,Type=String,Description=\"Pantree variant type: SNP/MNP/INS/DEL/DUP/REPL/INV/NON_REF\">\n"
        << "##INFO=<ID=NR,Number=0,Type=Flag,Description=\"Non-reference allele: enter/exit node off linear-reference chain\">\n"
        << "##INFO=<ID=TOPO,Number=1,Type=String,Description=\"Triallelic topology: PROPERLY_TRIALLELIC/OVERLAPPING/NESTED/INTERLOCKING\">\n"
        << "##INFO=<ID=OFF_REF_TIER,Number=1,Type=String,Description=\"Off-reference novelty tier: NOVEL/NOVEL_WEAK/DIVERGED/OFF_REF_KNOWN\">\n"
        << "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Mate contig for translocation\">\n"
        << "##INFO=<ID=POS2,Number=1,Type=Integer,Description=\"Mate position for translocation\">\n"
        << "##INFO=<ID=END2,Number=1,Type=Integer,Description=\"Mate end position for translocation\">\n"
        << "##INFO=<ID=MATE_CLADE,Number=1,Type=String,Description=\"Mate routed clade\">\n"
        << "##INFO=<ID=MATE_OFFREF,Number=0,Type=Flag,Description=\"Mate breakpoint is off-reference\">\n"
        << "##INFO=<ID=EC,Number=1,Type=String,Description=\"Element class: NONE/REPEAT/TE_LTR/TE_TIR/TE_LINE/TE_SINE/STARSHIP/HGT/RIP\">\n"
        << "##INFO=<ID=QMODE,Number=1,Type=String,Description=\"Query input mode: assembly/long-reads/short-reads\">\n"
        << "##INFO=<ID=FUSED_POST,Number=1,Type=Float,Description=\"Posterior alt probability after probabilistic evidence fusion\">\n"
        << "##INFO=<ID=FUSED_LOGODDS,Number=1,Type=Float,Description=\"Log-odds for the alternative allele after evidence fusion\">\n"
        << "##INFO=<ID=FUSED_DEPTH,Number=1,Type=Float,Description=\"Effective depth used by evidence fusion\">\n"
        << "##INFO=<ID=FUSED_LAYERS,Number=1,Type=Integer,Description=\"Number of evidence observations fused for the chosen call\">\n"
        << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        << "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"GQ\">\n"
        << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
}

// ============================================================
// Write one VariantCallBridge to VCF stream
// ============================================================
static void write_vcf_record(std::ostream& out,
                              const VariantCallBridge& v,
                              int id) {
    const int pos = std::max(1, v.pos);
    const int svlen = normalized_svlen_for_output(v);
    const int end = normalized_end_for_output(v);
    const bool primaryOffRef = (v.type == "OFF_REF" || v.refAsm == "OFF_REFERENCE");
    out << v.qContig
        << '\t' << pos
        << '\t' << "sv" << id
        << "\tN"
        << "\t<" << v.type << '>'
        << '\t' << std::fixed << std::setprecision(1) << v.gq
        << '\t' << (v.gq >= 10.0 ? "PASS" : "LowConf")
        << "\tSVTYPE=" << v.type
        << ";SVLEN=" << svlen
        << ";END=" << end
        << ";ANNOT=" << v.annotation
        << ";CLADE=" << v.refAsm
        << ";REFCONTIG=" << (v.refContig.empty() ? "." : v.refContig)
        << ";PHYLUM=" << (v.phylum.empty() ? "." : v.phylum)
        << ";CLADE_RANK=" << (v.cladeRank.empty() ? "." : v.cladeRank)
        << ";BSCORE=" << std::fixed << std::setprecision(2) << v.blockScore
        << ";MAPQ=" << std::fixed << std::setprecision(1) << v.mapq
        << ";QASM=" << v.qAsm
        << ";QMODE=" << v.queryMode
        << ";FUSED_POST=" << std::fixed << std::setprecision(4) << v.fusedPosteriorAlt
        << ";FUSED_LOGODDS=" << std::fixed << std::setprecision(4) << v.fusedLogOddsAlt
        << ";FUSED_DEPTH=" << std::fixed << std::setprecision(3) << v.fusedEffectiveDepth
        << ";FUSED_LAYERS=" << std::max(0, v.fusedLayersUsed);
    if (v.refPos > 0) out << ";REFPOS=" << v.refPos;
    if (v.refEnd > 0) out << ";REFEND=" << v.refEnd;
    if (primaryOffRef) out << ";OFFREF";
    // Emit OFF_REF_TIER for off-reference calls — annotation already holds the tier string
    // (NOVEL / NOVEL_WEAK / DIVERGED / OFF_REF_KNOWN) set by make_off_reference_call_scored.
    if (primaryOffRef &&
        (v.annotation == "NOVEL" || v.annotation == "NOVEL_WEAK" ||
         v.annotation == "DIVERGED" || v.annotation == "OFF_REF_KNOWN"))
        out << ";OFF_REF_TIER=" << v.annotation;
    // EC: element class — emitted for all calls; meaningful for OFF_REF
    if (!v.elementClass.empty() && v.elementClass != "NONE")
        out << ";EC=" << v.elementClass;
    if (!v.pantreeClass.empty() && v.pantreeClass != ".")
        out << ";VT=" << v.pantreeClass;
    if (v.isNonRefVariant) out << ";NR";
    if (!v.triallelicTopology.empty() && v.triallelicTopology != ".")
        out << ";TOPO=" << v.triallelicTopology;
    if (is_translocation_type(v.type)) {
        out << ";CHR2=" << (v.mateContig.empty() ? "." : v.mateContig)
            << ";POS2=" << std::max(1, v.matePos)
            << ";END2=" << std::max(std::max(1, v.matePos), v.mateEnd)
            << ";MATE_CLADE=" << (v.mateRefAsm.empty() ? "." : v.mateRefAsm);
        if (v.mateOffReference) out << ";MATE_OFFREF";
    }
    out << "\tGT:GQ"
        << '\t' << v.genotype << ':' << std::fixed
        << std::setprecision(1) << v.gq << '\n';
}

// ============================================================
// Write one VariantCallBridge to TSV stream
// ============================================================
static void write_tsv_record(std::ostream& out,
                              const VariantCallBridge& v) {
    out << v.qAsm     << '\t' << v.qContig    << '\t' << v.type
        << '\t' << v.refAsm    << '\t' << v.refContig
        << '\t' << v.refPos    << '\t' << v.refEnd
        << '\t' << v.pos       << '\t' << v.end       << '\t' << v.svlen
        << '\t' << v.blockScore<< '\t' << v.anchors
        << '\t' << v.genotype  << '\t' << v.gq
        << '\t' << v.annotation<< '\t' << v.alignmentMode
        << '\t' << v.queryMode
        << '\t' << v.fusedPosteriorAlt
        << '\t' << v.fusedLogOddsAlt
        << '\t' << v.fusedEffectiveDepth
        << '\t' << v.fusedLayersUsed
        << '\n';
}

// ============================================================
// Build query_input::InputConfig from parsed Options
// ============================================================
static query_input::InputConfig make_input_config(const Options& o) {
    query_input::InputConfig cfg;
    cfg.lrAnchorK       = o.lrAnchorK;
    cfg.lrMinCluster    = static_cast<size_t>(std::max(1, o.lrMinCluster));
    cfg.lrMinReadLen    = static_cast<size_t>(std::max(0, o.lrMinReadLen));
    cfg.lrMaxReadLen    = static_cast<size_t>(std::max(1, o.lrMaxReadLen));
    cfg.srK             = o.srK;
    cfg.srMinKmerFreq   = static_cast<uint32_t>(std::max(0, o.srMinKmerFreq));
    cfg.srMinUnitigLen  = static_cast<size_t>(std::max(1, o.srMinUnitigLen));
    cfg.srMinReadLen    = static_cast<size_t>(std::max(0, o.srMinReadLen));
    cfg.srMaxReadLen    = static_cast<size_t>(std::max(1, o.srMaxReadLen));
    cfg.genomeSizeHint  = o.genomeSizeHint;
    cfg.maxReadsPerFile = o.maxReadsPerFile > 0 ? o.maxReadsPerFile : 10000000;
    return cfg;
}

// Resolve QueryMode + autoDetect flag from the user's --query-mode string.
//   empty or "auto" -> (ASSEMBLY, true)   means: let prepare_query auto-detect
//   any other value -> (parsed mode, false) means: honour the user's explicit choice
static std::pair<query_input::QueryMode, bool>
resolve_query_mode(const std::string& s) {
    if (s.empty() || s == "auto")
        return {query_input::QueryMode::ASSEMBLY, true};
    return {query_input::parse_mode(s), false};
}

// Apply mode-specific calling-parameter overrides to a local Options copy.
// This runs only when applyModeParamOverride is true (default).
static void apply_mode_overrides(Options& o,
                                 query_input::QueryMode mode,
                                 const query_input::CoverageReport& report) {
    if (!o.applyModeParamOverride) return;
    const auto tier = report.coverageTier;
    switch (mode) {
        case query_input::QueryMode::ASSEMBLY:
            if (tier == query_input::CoverageTier::LOW) {
                if (!o.quiet)
                    std::cerr << "[query-input] assembly low-coverage overrides: "
                                 "k=17 chainGapBand=10000 minAnchors=1 minBlockScore=4.0 secondaryK=13\n";
                o.k                 = 17;
                o.chainGapBand      = 10000;
                o.minAnchors        = 1;
                o.minBlockScore     = 4.0;
                o.useSecondarySeeds = true;
                o.secondaryK        = 13;
            } else if (tier == query_input::CoverageTier::HIGH) {
                if (!o.quiet)
                    std::cerr << "[query-input] assembly high-coverage overrides: "
                                 "k=23 minAnchors=3 minBlockScore=7.0\n";
                o.k             = 23;
                o.minAnchors    = std::max(o.minAnchors, 3);
                o.minBlockScore = std::max(o.minBlockScore, 7.0);
            }
            break;

        case query_input::QueryMode::LONG_READS:
            if (tier == query_input::CoverageTier::LOW) {
                if (!o.quiet)
                    std::cerr << "[query-input] long-reads low-coverage overrides: "
                                 "k=13 chainGapBand=20000 minAnchors=1 minBlockScore=2.5 secondaryK=9\n";
                o.k                      = 13;
                o.chainGapBand           = 20000;
                o.minAnchors             = 1;
                o.minBlockScore          = 2.5;
                o.useSecondarySeeds      = true;
                o.secondaryK             = 9;
                o.repeatRescueMinAnchors = std::min(o.repeatRescueMinAnchors, 2);
            } else if (tier == query_input::CoverageTier::HIGH) {
                if (!o.quiet)
                    std::cerr << "[query-input] long-reads high-coverage overrides: "
                                 "k=17 chainGapBand=12000 minAnchors=2 minBlockScore=4.0 secondaryK=13\n";
                o.k                 = 17;
                o.chainGapBand      = 12000;
                o.minAnchors        = std::max(o.minAnchors, 2);
                o.minBlockScore     = std::max(o.minBlockScore, 4.0);
                o.useSecondarySeeds = true;
                o.secondaryK        = 13;
            } else {
                if (!o.quiet)
                    std::cerr << "[query-input] long-reads standard overrides: "
                                 "k=15 chainGapBand=15000 minAnchors=1 minBlockScore=3.0 secondaryK=11\n";
                o.k                 = 15;
                o.chainGapBand      = 15000;
                o.minAnchors        = 1;
                o.minBlockScore     = 3.0;
                o.useSecondarySeeds = true;
                o.secondaryK        = 11;
            }
            break;

        case query_input::QueryMode::SHORT_READS:
            if (tier == query_input::CoverageTier::LOW) {
                if (!o.quiet)
                    std::cerr << "[query-input] short-reads low-coverage overrides: "
                                 "k=17 minAnchors=1 minBlockScore=3.0 secondaryK=13 chainGapBand=8000\n";
                o.k                 = 17;
                o.minAnchors        = 1;
                o.minBlockScore     = 3.0;
                o.chainGapBand      = 8000;
                o.useSecondarySeeds = true;
                o.secondaryK        = 13;
            } else if (tier == query_input::CoverageTier::HIGH) {
                if (!o.quiet)
                    std::cerr << "[query-input] short-reads high-coverage overrides: "
                                 "k=25 minAnchors=4 minBlockScore=7.0 secondaryK=17\n";
                o.k                 = 25;
                o.minAnchors        = std::max(o.minAnchors, 4);
                o.minBlockScore     = std::max(o.minBlockScore, 7.0);
                o.useSecondarySeeds = true;
                o.secondaryK        = 17;
            } else {
                if (!o.quiet)
                    std::cerr << "[query-input] short-reads standard overrides: "
                                 "minAnchors=3 secondarySeeds=on secondaryK=15\n";
                o.minAnchors        = 3;
                o.useSecondarySeeds = true;
                o.secondaryK        = std::max(o.secondaryK, 15);
            }
            break;
    }
}

// ============================================================
// Process one query file — returns its calls + contigs for sidecar output
// ============================================================
struct QueryResult {
    std::string qAsm;
    std::unordered_map<std::string, std::string> contigs;
    std::unordered_map<std::string, std::string> contigLookup;
    std::vector<VariantCallBridge> calls;
};

static QueryResult
process_query(const std::string& qAsmPath,
              const Options& o,
              const tol::FederatedOptions& fo,
              const SimpleRefIndex* refIdx = nullptr) {
    QueryResult qr;
    qr.qAsm = fs::path(qAsmPath).stem().string();

    const query_input::InputConfig qcfg = make_input_config(o);
    const auto [modeHint, autoDetect]   = resolve_query_mode(o.queryMode);

    query_input::PreparedQuery prep;
    try {
        prep = query_input::prepare_query(qAsmPath, modeHint, qcfg, autoDetect, o.quiet);
    } catch (const std::exception& e) {
        std::cerr << "[warn] skipping " << qAsmPath << ": " << e.what() << '\n';
        return qr;
    }

    if (prep.contigs.empty()) {
        if (!o.quiet)
            std::cerr << "[warn] " << qr.qAsm
                      << ": no sequences after preprocessing\n";
        return qr;
    }

    qr.qAsm   = prep.sampleName;
    qr.contigs = std::move(prep.contigs);
    qr.contigLookup.reserve(qr.contigs.size() + 1);
    for (const auto& kv : qr.contigs)
        qr.contigLookup.emplace(kv.first, kv.first);

    const std::string modeLabel = query_input::mode_name(prep.report.mode);
    Options qo = o;
    apply_mode_overrides(qo, prep.report.mode, prep.report);
    const tol::FederatedOptions eff_fo =
        (prep.report.mode != query_input::QueryMode::ASSEMBLY)
        ? make_fed_opts(qo) : fo;

    std::unordered_map<std::string, std::vector<VariantCallBridge>> candidates;
    candidates.reserve(qr.contigs.size() * 2 + 1);
    std::optional<RefSearchCache> refCache;
    std::optional<SingleRefMemCache> singleRefMemCache;
    if (refIdx != nullptr) refCache.emplace(*refIdx);
    if (refIdx != nullptr) {
        const size_t saMaxBytes = static_cast<size_t>(std::max(0, o.saMaxContigMB)) * 1024 * 1024;
        singleRefMemCache.emplace(saMaxBytes);
    }
    auto add_candidates = [&](std::vector<VariantCallBridge>&& calls) {
        for (auto& c : calls) {
            c.queryMode = modeLabel;
            candidates[c.qContig].push_back(std::move(c));
        }
    };

    if (qo.useTolHierarchical && tol::TolGlobal::instance().is_initialized()) {
        auto hcalls = qo.tolMultiRank
            ? tol::hierarchical_call_assembly_multirank(
                  qr.qAsm, qr.contigs, eff_fo,
                  static_cast<size_t>(std::max(1, qo.routingTopN)))
            : tol::hierarchical_call_assembly(qr.qAsm, qr.contigs, eff_fo);
        add_candidates(std::move(hcalls));
    }

    if (refIdx != nullptr) {
        if (prep.report.mode == query_input::QueryMode::ASSEMBLY) {
            add_candidates(mem_chain_sv_calls(
                qr.qAsm, qr.contigs, *refIdx, *refCache, *singleRefMemCache,
                qo, eff_fo, prep.report.mode));
            add_candidates(simple_length_fallback_calls(qr.qAsm, qr.contigs, *refIdx, qo));
        } else {
            add_candidates(reads_mode_sv_calls(
                qr.qAsm, qr.contigs, *refIdx, *refCache, *singleRefMemCache,
                qo, eff_fo, prep.report.mode));
        }

        add_candidates(simple_offref_fallback_calls(
            qr.qAsm, qr.contigs, *refIdx, *refCache, qo, prep.report.mode));

        qr.calls = select_best_call_per_contig(qr.contigs, candidates, *refIdx, *refCache, qo,
                                               prep.report.mode, modeLabel,
                                               &prep.report);
        if (prep.report.mode != query_input::QueryMode::ASSEMBLY) {
            qr.calls = deduplicate_read_mode_events(std::move(qr.calls));
            qr.calls = filter_low_coverage_read_artifacts(std::move(qr.calls), prep.report);
        }
    } else {
        for (auto& kv : candidates) {
            if (!kv.second.empty()) qr.calls.push_back(std::move(kv.second.front()));
        }
    }

    if (!qr.calls.empty()) return qr;
    if (!qo.quiet)
        std::cerr << "[info] " << qr.qAsm << ": no calls emitted\n";
    return qr;
}

static const std::string* find_query_contig_seq(const QueryResult& qr,
                                             const std::string& qContig) {
    auto exact = qr.contigs.find(qContig);
    if (exact != qr.contigs.end()) return &exact->second;
    auto alias = qr.contigLookup.find(qContig);
    if (alias == qr.contigLookup.end()) return nullptr;
    auto it = qr.contigs.find(alias->second);
    return (it == qr.contigs.end()) ? nullptr : &it->second;
}

static void write_gfa_segments(std::ostream& out,
                               const QueryResult& qr,
                               const std::vector<VariantCallBridge>& calls,
                               std::unordered_set<std::string>& seen) {
    auto sanitize = [](std::string value) {
        for (char& ch : value) {
            if (ch == '\t' || ch == ' ' || ch == '\n' || ch == '\r') ch = '_';
        }
        return value;
    };

    // emit_segment: write one GFA S-line + P-line.
    // elementClass carries the ElementClass tag (EC:Z:) so repeat/TE/HGT/RIP
    // annotation is visible to downstream GFA consumers.  For non-OFF_REF calls
    // elementClass is "NONE" and the tag is still written for uniformity.
    auto emit_segment = [&](const std::string& rawId,
                            const std::string& contig,
                            int pos,
                            int end,
                            int svlen,
                            const std::string& annotation,
                            const std::string& type,
                            const std::string& clade,
                            const std::string& elementClass,
                            bool forcePlaceholder = false) {
        std::string seg = sanitize(rawId);
        if (!seen.insert("S:" + seg).second) return seg;
        const std::string* seq = contig.empty() || contig == "." ? nullptr : find_query_contig_seq(qr, contig);
        size_t start0 = pos > 0 ? static_cast<size_t>(pos - 1) : 0;
        size_t length = static_cast<size_t>(normalized_interval_len(pos, end, svlen, 1));
        out << "S\t" << seg << "\t";
        if (!forcePlaceholder && seq && start0 < seq->size() && length > 0) {
            size_t end0 = std::min(seq->size(), start0 + length);
            if (end0 > start0) {
                out.write(seq->data() + start0, static_cast<std::streamsize>(end0 - start0));
            } else {
                out << '*';
            }
        } else {
            out << '*';
        }
        out << "\tAN:Z:" << annotation
            << "\tVT:Z:" << type
            << "\tCL:Z:" << clade
            << "\tEC:Z:" << (elementClass.empty() ? "NONE" : elementClass)
            << "\n";
        out << "P\t" << seg << "\t" << seg << "+\t*\n";
        return seg;
    };

    for (const auto& c : calls) {
        const bool primaryOffRef = (c.type == "OFF_REF" || c.refAsm == "OFF_REFERENCE");
        if (c.type == "TRA") {
            const int end = normalized_end_for_output(c);
            const std::string leftId = qr.qAsm + ":" + c.qContig + ":" +
                                       std::to_string(std::max(1, c.pos)) + "-" +
                                       std::to_string(end) + ":TRA";
            const std::string left = emit_segment(leftId,
                                                  c.qContig,
                                                  c.pos,
                                                  end,
                                                  c.svlen,
                                                  primaryOffRef ? "OFF_REFERENCE" : c.annotation,
                                                  c.type,
                                                  c.refAsm,
                                                  c.elementClass,
                                                  false);

            const int matePos = std::max(1, c.matePos);
            const int mateEnd = std::max(matePos, c.mateEnd);
            const std::string mateContig = c.mateContig.empty() ? "." : c.mateContig;
            const std::string rightId = qr.qAsm + ":" + mateContig + ":" +
                                        std::to_string(matePos) + "-" +
                                        std::to_string(mateEnd) + ":TRA_MATE";
            const std::string right = emit_segment(rightId,
                                                   mateContig,
                                                   matePos,
                                                   mateEnd,
                                                   c.mateSvLen,
                                                   c.mateOffReference ? "OFF_REFERENCE" : "TRA_MATE",
                                                   "TRA_MATE",
                                                   c.mateOffReference ? "OFF_REFERENCE" : c.mateRefAsm,
                                                   c.mateOffReference ? c.elementClass : "NONE",
                                                   false);
            const std::string linkKey = "L:" + left + "->" + right;
            if (seen.insert(linkKey).second) {
                out << "L\t" << left << "\t+\t" << right << "\t+\t0M"
                    << "\tVT:Z:TRA"
                    << "\tCL:Z:" << c.refAsm
                    << "\tMCL:Z:" << (c.mateOffReference ? "OFF_REFERENCE" : c.mateRefAsm)
                    << "\n";
            }
            continue;
        }

        const int svlen = normalized_svlen_for_output(c);
        const int end = normalized_end_for_output(c);
        const bool placeholderOnly = (c.type == "DEL");
        emit_segment(qr.qAsm + ":" + c.qContig + ":" + std::to_string(std::max(1, c.pos)) + "-" + std::to_string(end),
                     c.qContig,
                     c.pos,
                     end,
                     svlen,
                     primaryOffRef ? "OFF_REFERENCE" : c.annotation,
                     c.type,
                     c.refAsm,
                     c.elementClass,
                     placeholderOnly);
    }
}

struct TolValidationRow {
    std::string cladeName;
    std::string cladeRank;
    std::string phylum;
    size_t genomesManifest = 0;
    size_t genomesIndexed = 0;
    size_t cidxShards = 0;
    size_t gbzGraphs = 0;
    size_t svBubbles = 0;
    size_t compressedBytes = 0;
    std::string status;
};

static int run_tol_validation(const Options& o) {
    if (o.tolBuildManifest.empty()) {
        std::cerr << "[error] --tol-validate-index requires --tol-build-index MANIFEST "
                     "(the same manifest path used when building the index)\n";
        return 1;
    }
    if (o.tolIndexDir.empty() || o.tolRegistryDir.empty()) {
        std::cerr << "[error] --tol-index-dir and --tol-registry-dir required\n";
        return 1;
    }

    std::unordered_map<std::string, size_t> manifestCounts;
    std::unordered_map<std::string, std::string> manifestRanks;
    std::unordered_map<std::string, std::string> manifestPhyla;
    std::ifstream in(o.tolBuildManifest);
    if (!in) {
        std::cerr << "[error] cannot open manifest: " << o.tolBuildManifest << '\n';
        return 1;
    }
    std::string line;
    bool extendedManifest = false;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            extendedManifest = (std::count(line.begin(), line.end(), '\t') >= 8);
            continue;
        }
        std::istringstream ss(line);
        std::string asmName, phylum, cls, order, family, genus, cladeName, cladeRank, fastaPath;
        std::getline(ss, asmName, '\t');
        std::getline(ss, phylum, '\t');
        if (extendedManifest) {
            std::getline(ss, cls, '\t');
            std::getline(ss, order, '\t');
            std::getline(ss, family, '\t');
            std::getline(ss, genus, '\t');
        }
        std::getline(ss, cladeName, '\t');
        std::getline(ss, cladeRank, '\t');
        std::getline(ss, fastaPath, '\t');
        if (!cladeName.empty()) {
            ++manifestCounts[cladeName];
            manifestRanks[cladeName] = cladeRank;
            manifestPhyla[cladeName] = phylum;
        }
    }

    tol::ManifestRegistry reg(o.tolRegistryDir);
    reg.load_from_disk();
    const auto& descs = reg.descriptors();
    size_t cidxCount = 0;
    bool hasRoutingManifest = fs::exists(fs::path(o.tolIndexDir) / "routing_manifest.tsv");
    if (fs::exists(o.tolIndexDir)) {
        for (auto const& ent : fs::recursive_directory_iterator(o.tolIndexDir))
            if (ent.is_regular_file() && ent.path().extension() == ".cidx") ++cidxCount;
    }

    std::vector<TolValidationRow> rows;
    rows.reserve(descs.size());
    bool ok = true;
    for (const auto& d : descs) {
        TolValidationRow row;
        row.cladeName = d.cladeName;
        row.cladeRank = !d.cladeRank.empty() ? d.cladeRank : manifestRanks[d.cladeName];
        row.phylum = !d.phylum.empty() ? d.phylum : manifestPhyla[d.cladeName];
        row.genomesManifest = manifestCounts[d.cladeName];
        row.genomesIndexed = d.genomeCount;
        row.cidxShards = cidxCount;
        row.gbzGraphs = fs::exists(d.graphPath) ? 1u : 0u;
        row.svBubbles = d.svBubbles;
        row.compressedBytes = d.compressedBytes;
        const fs::path graphPath(d.graphPath);
        const fs::path registryDir = graphPath.parent_path();
        const std::string stem = graphPath.stem().string();
        const fs::path gbwtPath = fs::path(o.tolIndexDir) / (stem + ".gbwt");
        const fs::path minPath  = fs::path(o.tolIndexDir) / (stem + ".min");
        const bool graphNonEmpty = fs::exists(d.graphPath) && row.compressedBytes > 0 && fs::file_size(d.graphPath) > 0;
        const bool hasSidecars = fs::exists(gbwtPath) && fs::exists(minPath) && fs::file_size(gbwtPath) > 0 && fs::file_size(minPath) > 0;
        row.status = (row.gbzGraphs == 1 && row.genomesIndexed > 0 && cidxCount > 0 && hasRoutingManifest && graphNonEmpty && hasSidecars && (row.genomesManifest == 0 || row.genomesIndexed <= row.genomesManifest)) ? "ok" : "check";
        if (row.status != "ok") ok = false;
        rows.push_back(std::move(row));
    }

    std::string report = o.tolValidationReport.empty() ? (o.tolRegistryDir + "/validation_report.tsv") : o.tolValidationReport;
    std::ofstream out(report);
    if (!out) {
        std::cerr << "[error] cannot write validation report: " << report << '\n';
        return 1;
    }
    out << "#clade_name\tphylum\tmanifest_genomes\tindexed_genomes\tcidx_shards\tgbz_graph\tsv_bubbles\tcompressed_bytes\tstatus\n";
    for (const auto& r : rows) {
        out << r.cladeName << '\t' << r.phylum << '\t' << r.genomesManifest << '\t'
            << r.genomesIndexed << '\t' << r.cidxShards << '\t' << r.gbzGraphs << '\t'
            << r.svBubbles << '\t' << r.compressedBytes << '\t' << r.status << '\n';
    }
    std::cerr << "[tol] validation report: " << report << " (" << rows.size() << " clades)\n";
    return ok ? 0 : 2;
}

// ============================================================
// main
// ============================================================
int main(int argc, char** argv) {
    Options o;
    try {
        o = parse_args(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "[error] " << e.what() << '\n';
        return 1;
    }

    // ---- TE train mode -----------------------------------------------
    if (o.teTrainMode) {
        if (o.queryList.empty()) {
            std::cerr << "[error] --te-train requires --query-list (labeled FASTA)\n";
            return 1;
        }
        if (o.teIndexPrefix.empty() && !o.outPrefix.empty()) {
            // fall back to --out-prefix as index prefix
            const_cast<std::string&>(o.teIndexPrefix) = o.outPrefix;
        }
        if (o.teIndexPrefix.empty()) {
            std::cerr << "[error] --te-train requires --te-index-prefix or --out-prefix\n";
            return 1;
        }
        te::TEClassifier::Params p;
        p.k          = o.teK;
        p.fracmin_p  = o.teFracminP;
        p.max_hashes = o.teMaxHashes;
        te::TEClassifier clf(p);
        try {
            auto ql = read_list(o.queryList);
            for (const auto& fa : ql) {
                std::cerr << "[te-train] loading " << fa << '\n';
                clf.train(fa, !o.quiet);
            }
            if (auto pd = fs::path(o.teIndexPrefix).parent_path(); !pd.empty())
                fs::create_directories(pd);
            clf.save(o.teIndexPrefix);
            std::cerr << "[te-train] saved " << clf.num_centroids()
                      << " centroids to " << o.teIndexPrefix << ".{vptree,meta}\n";
        } catch (const std::exception& e) {
            std::cerr << "[error] te-train: " << e.what() << '\n';
            return 1;
        }
        return 0;
    }

    // ---- TE classify mode --------------------------------------------
    if (o.teClassifyMode) {
        if (o.queryList.empty()) {
            std::cerr << "[error] --te-classify requires --query-list (FASTA to classify)\n";
            return 1;
        }
        const std::string idx = o.teIndexPrefix.empty() ? o.outPrefix : o.teIndexPrefix;
        if (idx.empty()) {
            std::cerr << "[error] --te-classify requires --te-index-prefix or --out-prefix\n";
            return 1;
        }
        te::TEClassifier::Params p;
        p.k          = o.teK;
        p.fracmin_p  = o.teFracminP;
        p.max_hashes = o.teMaxHashes;
        te::TEClassifier clf(p);
        try {
            clf.load(idx);
            std::cerr << "[te-classify] loaded " << clf.num_centroids() << " centroids\n";
            const std::string out_tsv = o.outPrefix.empty() ? "te_predictions.tsv"
                                                             : o.outPrefix + ".te_predictions.tsv";
            auto ql = read_list(o.queryList);
            for (const auto& fa : ql) {
                std::cerr << "[te-classify] classifying " << fa << '\n';
                clf.classify_fasta(fa, out_tsv);
            }
            std::cerr << "[te-classify] predictions written to " << out_tsv << '\n';
        } catch (const std::exception& e) {
            std::cerr << "[error] te-classify: " << e.what() << '\n';
            return 1;
        }
        return 0;
    }

    // ---- TOL validation mode -----------------------------------------
    if (o.useTolHierarchical && o.tolValidateIndex) {
        return run_tol_validation(o);
    }

    // ---- TOL index build mode ----------------------------------------
    if (o.useTolHierarchical && o.tolBuildIndex) {
        if (o.tolBuildManifest.empty()) {
            std::cerr << "[error] --tol-build-index requires a manifest path\n";
            return 1;
        }
        if (o.tolIndexDir.empty() || o.tolRegistryDir.empty()) {
            std::cerr << "[error] --tol-index-dir and --tol-registry-dir required\n";
            return 1;
        }

        tol::SyncmerParams sp;
        sp.k               = o.k;
        sp.s               = o.s;
        sp.t               = o.t;
        sp.stride          = static_cast<size_t>(o.seedStride);
        sp.useIntervalHash = o.useIntervalHash;
        sp.ihWing          = o.ihWing;
        sp.ihMaxDist       = o.ihMaxDist;
        sp.ihResolution    = o.ihResolution;

        const tol::FederatedOptions fo = make_fed_opts(o);
        try {
            if (o.tolMultiRank) {
                tol::build_multi_rank_index_from_manifest(
                    o.tolBuildManifest,
                    o.tolIndexDir,
                    o.tolRegistryDir,
                    sp,
                    o.tolRoutingDensity,
                    !o.quiet,
                    o.annotationTsv,
                    o.tolMaxCladeGenomes,
                    static_cast<size_t>(o.tolIndexThreads),
                    &fo.fallbackSketchParams,
                    o.tolBaseGraphBuild);
            } else {
                tol::build_tol_index_from_manifest(
                    o.tolBuildManifest,
                    o.tolIndexDir,
                    o.tolRegistryDir,
                    sp,
                    o.tolRoutingDensity,
                    !o.quiet,
                    o.annotationTsv,
                    o.tolMaxCladeGenomes,
                    static_cast<size_t>(o.tolIndexThreads),
                    &fo.fallbackSketchParams,
                    o.tolBaseGraphBuild);
            }
        } catch (const std::exception& e) {
            std::cerr << "[error] index build failed: " << e.what() << '\n';
            return 1;
        }
        return 0;
    }

    // ---- Initialise TOL global state ------------------------------------
    if (o.useTolHierarchical) {
        if (o.tolIndexDir.empty() || o.tolRegistryDir.empty()) {
            std::cerr << "[error] --tol-index-dir and --tol-registry-dir required\n";
            return 1;
        }
        try {
            tol::TolGlobal::instance().init(
                o.tolIndexDir,
                o.tolRegistryDir,
                o.tolCacheGB << 30,
                o.tolCacheEntries);
            if (o.tolMultiRank || o.tolAncestralAlign) {
                tol::MultiRankIndex::instance().init(
                    o.tolIndexDir,
                    o.tolRegistryDir,
                    o.tolCacheGB << 30,
                    o.tolCacheEntries);
            }
        } catch (const std::exception& e) {
            std::cerr << "[error] TOL init failed: " << e.what() << '\n';
            return 1;
        }
    }

    // ---- Load query list ------------------------------------------------
    if (o.queryList.empty()) {
        std::cerr << "[error] --query-list is required\n";
        return 1;
    }
    auto queries = read_list(o.queryList);
    if (queries.empty()) {
        std::cerr << "[error] query list is empty: " << o.queryList << '\n';
        return 1;
    }

    if (o.outPrefix.empty()) {
        std::cerr << "[error] --out-prefix is required\n";
        return 1;
    }

    // ---- Open outputs ---------------------------------------------------
    std::string tsv_path = o.outPrefix + ".hits.tsv";
    std::string vcf_path = o.outPrefix + ".vcf";
    std::string gfa_path = o.outPrefix + ".gfa";
    std::string ancestral_path = o.tolAncestralOut.empty() ? (o.outPrefix + ".ancestral.tsv") : o.tolAncestralOut;
    if (auto out_parent = fs::path(o.outPrefix).parent_path(); !out_parent.empty()) {
        fs::create_directories(out_parent);
    }
    std::ofstream tsv_out(tsv_path);
    std::ofstream vcf_out(vcf_path);
    std::ofstream anc_out;
    std::ofstream gfa_out;
    if (o.tolAncestralAlign) {
        anc_out.open(ancestral_path);
        if (!anc_out) { std::cerr << "[error] cannot open " << ancestral_path << "\n"; return 1; }
    }
    if (!o.noGfa) {
        gfa_out.open(gfa_path);
        if (!gfa_out) { std::cerr << "[error] cannot open " << gfa_path << "\n"; return 1; }
        gfa_out << "H\tVN:Z:1.0\n";
    }
    if (!tsv_out) { std::cerr << "[error] cannot open " << tsv_path << '\n'; return 1; }
    if (!vcf_out) { std::cerr << "[error] cannot open " << vcf_path << '\n'; return 1; }

    tsv_out << "query_asm\tquery_contig\ttype\tref_asm\tref_contig"
               "\tref_pos\tref_end\tpos\tend\tsvlen\tblock_score\tanchors"
               "\tgenotype\tgq\tannotation\talignment_mode\tquery_mode"
               "\tfused_posterior_alt\tfused_logodds_alt\tfused_effective_depth\tfused_layers\n";
    write_vcf_header(vcf_out, "fungi_graphsv_tol_v3");
    if (o.tolAncestralAlign) tol::write_ancestral_tsv_header(anc_out);

    const tol::FederatedOptions fo = make_fed_opts(o);
    const SimpleRefIndex simpleRefIdx = load_simple_ref_index(o);
    std::optional<tol::AncestralManifestContext> ancestralCtx;
    if (o.tolAncestralAlign) {
        if (o.tolManifest.empty()) {
            std::cerr << "[error] --tol-ancestral-align requires --tol-manifest PATH\n";
            return 1;
        }
        try {
            ancestralCtx = tol::load_ancestral_manifest_context(o.tolManifest);
        } catch (const std::exception& e) {
            std::cerr << "[error] failed to load ancestral manifest: " << e.what() << '\n';
            return 1;
        }
    }

    // ---- Process queries (parallel if requested) ------------------------
    std::mutex           tsv_mutex;
    std::mutex           vcf_mutex;
    std::mutex           gfa_mutex;
    std::mutex           anc_mutex;
    std::atomic<int>     sv_id{1};
    std::atomic<size_t>  total_calls{0};
    std::atomic<size_t>  queries_done{0};
    const size_t         n_queries = queries.size();

    std::unordered_set<std::string> gfa_seen;
    auto process_one = [&](const std::string& qpath) {
        auto qr = process_query(qpath, o, fo, &simpleRefIdx);

        for (const auto& c : qr.calls) {
            int cur_id = sv_id++;
            {
                std::lock_guard<std::mutex> lk(tsv_mutex);
                write_tsv_record(tsv_out, c);
            }
            {
                std::lock_guard<std::mutex> lk(vcf_mutex);
                write_vcf_record(vcf_out, c, cur_id);
            }
        }
        if (!o.noGfa) {
            std::lock_guard<std::mutex> lk(gfa_mutex);
            write_gfa_segments(gfa_out, qr, qr.calls, gfa_seen);
        }
        if (o.tolAncestralAlign && ancestralCtx.has_value()) {
            std::lock_guard<std::mutex> lk(anc_mutex);
            tol::write_ancestral_alignments_for_assembly(
                qr.qAsm, qr.contigs, qr.calls, *ancestralCtx, fo,
                static_cast<size_t>(std::max(1, o.routingTopN)), anc_out);
        }
        total_calls  += qr.calls.size();
        size_t done  = ++queries_done;
        if (!o.quiet && done % 10 == 0)
            std::cerr << "[progress] " << done << '/' << n_queries
                      << " queries, " << total_calls.load() << " calls\n";
    };

    if (o.threads <= 1) {
        for (const auto& q : queries) process_one(q);
    } else {
        std::vector<std::thread> workers;
        std::atomic<size_t> idx{0};
        for (int t = 0; t < o.threads; ++t) {
            workers.emplace_back([&]() {
                for (;;) {
                    size_t i = idx++;
                    if (i >= queries.size()) break;
                    process_one(queries[i]);
                }
            });
        }
        for (auto& w : workers) w.join();
    }

    if (!o.quiet)
        std::cerr << "[done] " << total_calls.load() << " SV calls across "
                  << n_queries << " queries\n"
                  << "  TSV: " << tsv_path << '\n'
                  << "  VCF: " << vcf_path << '\n';

    return 0;
}
