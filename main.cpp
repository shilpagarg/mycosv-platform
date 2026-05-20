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
#include <cmath>
#include <csignal>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
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
#include <unistd.h>
#include <utility>
#include <vector>

// TOL headers (adjust -I path or copy alongside this file)
#include "fungi_tol_bridge.hpp"
#include "query_input_handler.hpp"
#include "te_classifier.hpp"

namespace fs = std::filesystem;

// ============================================================
// Graceful shutdown on SIGTERM/SIGINT.
//
// The wrapper run_real_fungal_benchmark.py sends SIGTERM (then SIGKILL after
// a 30 s grace period) when --mycosv-tool-timeout fires. Without a handler,
// the SIGKILL discards any in-memory candidate calls that haven't reached
// the per-query flush in main(). The handler below sets a flag the worker
// loop polls, and best-effort flushes the open ofstreams so partial output
// is preserved up to the last per-query flush.
// ============================================================
static std::atomic<bool> g_shutdown_requested{false};
static std::ofstream*    g_signal_tsv_out      = nullptr;
static std::ofstream*    g_signal_vcf_out      = nullptr;
static std::ofstream*    g_signal_hier_tsv_out = nullptr;
static std::ofstream*    g_signal_hier_vcf_out = nullptr;
static std::ofstream*    g_signal_gfa_out      = nullptr;
static std::ofstream*    g_signal_anc_out      = nullptr;

static void on_terminate_signal(int sig) {
    // async-signal-safe portion: set the atomic flag and write a fixed message.
    g_shutdown_requested.store(true, std::memory_order_relaxed);
    static const char msg[] =
        "[mycosv] caught signal, flushing partial outputs before exit\n";
    ssize_t r = ::write(2, msg, sizeof(msg) - 1);
    (void)r;
    // ofstream::flush() is not formally async-signal-safe, but in glibc/libstdc++
    // is in practice and is the only thing standing between us and lost data
    // before the wrapper escalates to SIGKILL.
    if (g_signal_tsv_out      && g_signal_tsv_out->good())      g_signal_tsv_out->flush();
    if (g_signal_vcf_out      && g_signal_vcf_out->good())      g_signal_vcf_out->flush();
    if (g_signal_hier_tsv_out && g_signal_hier_tsv_out->good()) g_signal_hier_tsv_out->flush();
    if (g_signal_hier_vcf_out && g_signal_hier_vcf_out->good()) g_signal_hier_vcf_out->flush();
    if (g_signal_gfa_out      && g_signal_gfa_out->good())      g_signal_gfa_out->flush();
    if (g_signal_anc_out      && g_signal_anc_out->good())      g_signal_anc_out->flush();
    // Re-raise with default handler so the parent sees the correct exit code.
    std::signal(sig, SIG_DFL);
    std::raise(sig);
}

static void install_signal_handlers() {
    std::signal(SIGTERM, on_terminate_signal);
    std::signal(SIGINT,  on_terminate_signal);
}

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
    // 10 Mb upper bound: fungal accessory chromosome PAV and STARSHIP cargo
    // reach into the megabase range (Aspergillus accessory chroms, Fusarium
    // dispensables). The previous 1 Mb cap silently dropped these.
    int    maxSvLen      = 10000000;
    double minBlockScore = 6.0;
    int    minAnchors    = 2;
    // Default 10 000 (was 128). Real fungal isolates can carry many
    // genome-wide SVs, so 128 capped recall at a fraction of the observable
    // burden. The cap is now a safety net rather than the dominant gate;
    // selection happens via the minBlockScore floor.
    int    maxCallsPerContig = 10000;

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
    // --max-ref-memory-mb caps the *raw* reference sequence pool stored in
    // the SimpleRefIndex; per-thread suffix-array caches are bounded
    // independently by --single-ref-cache-mb (LRU).  Hitting this cap
    // silently drops ref contigs from the search space, so we want it large
    // enough to cover a full routed sub-graph on a typical fungal panel:
    // 32 GiB fits ~30 multi-Gbp basidiomycete refs or ~2000 yeast refs in
    // a single sub-graph load and is still well within a 128 GiB cgroup.
    int    maxRefMemoryMB       = 32768; // --max-ref-memory-mb: cap total ref seq loaded
    // --single-ref-cache-mb: per-query-thread cap on the SingleRefMemCache
    // that holds suffix arrays for refs touched while calling SVs against a
    // query.  Each cached SA is ~13x the raw seq bytes, so an unbounded cache
    // crossed with `--threads` parallel queries blows past job memory on
    // multi-Gbp ref panels.  When this cap is exceeded the LRU entry is
    // evicted; 0 disables the cap (legacy behavior).
    int    singleRefCacheMB     = 1024;
    bool   noFlatRefFallback    = false; // --no-flat-ref-fallback: skip memory-heavy flat ref fallback
    // --skip-flat-if-hier-calls N: if the hierarchical phase produced at least N
    // calls for a query, skip the flat-MEM-chain fallback for that query.
    // Lets a user keep the flat fallback as a safety net for queries that
    // hierarchical can't route, without paying its 6 000+ contig cost when
    // hierarchical already produced a usable callset. 0 = disabled.
    int    skipFlatIfHierCalls  = 0;
    // --max-flat-ref-contigs N: hard cap on the total unique contig names
    // retained in the SimpleRefIndex used by the flat-MEM-chain fallback.
    // Fungal reference FASTAs commonly contain hundreds of unplaced
    // scaffolds, so an "8-ref" cap can still load 6 000+ contigs.
    // 0 = unbounded (legacy behaviour).
    int    maxFlatRefContigs    = 0;

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
    int    srMinReadLen     = 35;            // --sr-min-read-len
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

    // Diagnostics: fast offline checks that exit before any SV calling.
    //   --diagnose registry   : load --tol-registry-dir / --tol-index-dir, report
    //                            descriptors, FASTAs found, FASTAs on disk, CR-in-path,
    //                            and the resulting TolGlobal allRefs_ size. Catches
    //                            the CRLF manifest / missing-FASTA class of bug in
    //                            seconds without going through SV calling.
    std::string diagnoseMode;
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
"  --max-ref-memory-mb INT  Cap total reference sequence memory in MB (default 32768)\n"
"  --single-ref-cache-mb INT Per-thread SA cache cap in MB (default 1024; 0 = unbounded)\n"
"  --no-flat-ref-fallback Skip flat ref loading/fallback when hierarchical routing is enabled\n"
"  --skip-flat-if-hier-calls INT  Skip flat-MEM-chain fallback for a query when the\n"
"                            hierarchical phase already produced >= INT calls (0=off)\n"
"  --max-flat-ref-contigs INT  Cap unique ref contigs retained in the flat fallback\n"
"                            SimpleRefIndex (0=unbounded). Useful for fungal refs that\n"
"                            ship hundreds of unplaced scaffolds per FASTA.\n"
"\n"
"Diagnostics:\n"
"  --diagnose MODE          Run a fast offline health check and exit.\n"
"                            MODE=registry: replay TolGlobal::init and report\n"
"                            descriptor count, FASTA existence, CR-in-path, and the\n"
"                            resulting allRefs_ size. Exit 0 OK, 1 broken, 2 usage.\n"
"                            Catches CRLF-mangled manifests and missing-FASTA bugs\n"
"                            in seconds, before any SV calling.\n"
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
"  --sr-min-read-len INT    Drop reads shorter than this bp (default 35)\n"
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
        else if (x == "--single-ref-cache-mb")         o.singleRefCacheMB  = std::stoi(need(x.c_str(),i));
        else if (x == "--no-flat-ref-fallback")        o.noFlatRefFallback = true;
        else if (x == "--skip-flat-if-hier-calls")     o.skipFlatIfHierCalls = std::stoi(need(x.c_str(),i));
        else if (x == "--max-flat-ref-contigs")        o.maxFlatRefContigs = std::stoi(need(x.c_str(),i));
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
        else if (x == "--diagnose")                    o.diagnoseMode         = need(x.c_str(),i);
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
    FastaStream fs(path);
    std::istream& in = fs.get();
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

static double kmer_best_containment_from_hashes(const std::unordered_set<uint64_t>& a,
                                                const std::unordered_set<uint64_t>& b) {
    if (a.empty() || b.empty()) return 0.0;
    const auto* small = &a;
    const auto* big   = &b;
    if (small->size() > big->size()) std::swap(small, big);
    size_t inter = 0;
    for (uint64_t h : *small)
        if (big->count(h)) ++inter;
    const double aContain = static_cast<double>(inter) / static_cast<double>(a.size());
    const double bContain = static_cast<double>(inter) / static_cast<double>(b.size());
    return std::max(aContain, bContain);
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
            const double jaccard = kmer_jaccard_from_hashes(qHashes, refHashes[i]);
            const double contain = kmer_best_containment_from_hashes(qHashes, refHashes[i]);
            const double frac = std::max(jaccard, contain);
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

// LRU cache of per-ref suffix arrays.  Bounded by `cacheMaxBytes` (sum of
// estimated SA + LCP + ISA + text bytes across cached entries) to keep the
// per-thread footprint predictable when --threads parallel queries each
// touch a large routed sub-graph.  Eviction order is least-recently-used.
class SingleRefMemCache {
public:
    explicit SingleRefMemCache(size_t saMaxBytes = 0, size_t cacheMaxBytes = 0)
        : saMaxBytes_(saMaxBytes), cacheMaxBytes_(cacheMaxBytes) {}

    SingleRefMemIndex& get(const RefContigInfo* info) {
        auto it = cache_.find(info);
        if (it != cache_.end()) {
            // Hit: promote to LRU front.
            lru_.splice(lru_.begin(), lru_, it->second.lruIt);
            return it->second.idx;
        }

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

        const size_t entryBytes = estimate_entry_bytes(idx);
        if (cacheMaxBytes_ > 0) {
            evict_until_fits(entryBytes);
        }

        lru_.push_front(info);
        Entry e{std::move(idx), lru_.begin(), entryBytes};
        currentBytes_ += entryBytes;
        auto [inserted, _] = cache_.emplace(info, std::move(e));
        return inserted->second.idx;
    }

private:
    struct Entry {
        SingleRefMemIndex                                idx;
        std::list<const RefContigInfo*>::iterator        lruIt;
        size_t                                           bytes = 0;
    };

    static size_t estimate_entry_bytes(const SingleRefMemIndex& idx) {
        // SA + LCP + ISA are vector<int> (4 bytes each), text is the
        // concatenated reference (1 byte/char).  Constant overhead for the
        // RefSeq metadata is negligible compared to multi-megabyte SAs.
        const size_t intArrays = (idx.sa.sa.size() + idx.sa.lcp.size() + idx.sa.isa.size())
                                 * sizeof(int);
        return idx.sa.text.size() + intArrays;
    }

    void evict_until_fits(size_t incoming) {
        // Always make room for the incoming entry; never evict everything if
        // the single incoming entry is itself larger than the budget — accept
        // the overrun rather than thrashing.
        while (!lru_.empty()
               && currentBytes_ + incoming > cacheMaxBytes_) {
            const RefContigInfo* victim = lru_.back();
            auto vit = cache_.find(victim);
            if (vit == cache_.end()) {
                // Defensive: keep LRU and map in sync even if we somehow get
                // out of step (shouldn't happen but cheaper than crashing).
                lru_.pop_back();
                continue;
            }
            currentBytes_ -= vit->second.bytes;
            lru_.pop_back();
            cache_.erase(vit);
        }
    }

    size_t saMaxBytes_      = 0;
    size_t cacheMaxBytes_   = 0;
    size_t currentBytes_    = 0;
    std::list<const RefContigInfo*>                          lru_;
    std::unordered_map<const RefContigInfo*, Entry>          cache_;
};

static bool try_mem_chain_call_single_ref_cached(
        const std::string& qAsm,
        const std::string& qContig,
        const std::string& qSeq,
        const RefContigInfo* refInfo,
        SingleRefMemCache& memCache,
        const tol::FederatedOptions& fo,
        VariantCallBridge& call);

// Multi-emit variant: returns ALL per-pair INS/DEL events from the MEM chain
// (not just the dominant one). TRA/INV/DUP still emit a single representative
// event, matching the legacy semantics. This is the primary recall fix for
// real diverged fungal genomes where each query chain contains tens of small
// SVs that the single-emit path collapsed into one call.
static bool try_mem_chain_call_single_ref_cached_multi(
        const std::string& qAsm,
        const std::string& qContig,
        const std::string& qSeq,
        const RefContigInfo* refInfo,
        SingleRefMemCache& memCache,
        const tol::FederatedOptions& fo,
        std::vector<VariantCallBridge>& calls);

static SimpleRefIndex load_simple_ref_index(const Options& o) {
    SimpleRefIndex idx;

    auto reps = read_list(o.repAsmList);
    auto refs = read_list(o.refList);
    const auto& paths = !reps.empty() ? reps : refs;
    idx.reserve(paths.size() * 4 + 1);

    const size_t maxRefBytes = static_cast<size_t>(std::max(0, o.maxRefMemoryMB)) * 1024 * 1024;

    // ── Visibility: announce the load up front. ──────────────────────────
    // Without this, the operator saw nothing between Python's "[mycosv]
    // flat-ref fallback ENABLED" line and the first per-query "[progress]"
    // line — 5–15 minutes of silence on million_real (256 fungal ref
    // FASTAs × ~30 Mb on RDS shared storage). That looked like a hang
    // even when the binary was working correctly.
    // Threads sized below; emit the announcement after sizing so the
    // message reflects what will actually happen (serial vs N threads).
    const auto t0 = std::chrono::steady_clock::now();

    // ── Parallel FASTA load. ─────────────────────────────────────────────
    // The previous serial loop was the dominant cost in the silent gap
    // before the first per-query progress line (one syscall + ~30 Mb
    // string allocation per ref × 256 refs ≈ 5–15 min on shared storage).
    // Each worker builds a thread-local SimpleRefIndex chunk, then we
    // merge serially. The memory cap is re-checked at merge time — a
    // worker may have produced more bytes than the cap allows, but we
    // accept the soft overshoot rather than blocking workers on a shared
    // atomic counter. The cap is a budget, not a hard guard.
    const unsigned hw = std::max(1u, std::thread::hardware_concurrency());
    unsigned nThreads = std::min<unsigned>(8u, std::min<unsigned>(hw, static_cast<unsigned>(paths.size())));
    if (const char* env = std::getenv("MYCOSV_REFLOAD_THREADS")) {
        try { nThreads = std::max<unsigned>(1u, std::min<unsigned>(64u, static_cast<unsigned>(std::stoi(env)))); }
        catch (...) {}
    }
    if (paths.size() <= 1) nThreads = 1;

    if (!o.quiet) {
        std::cerr << "[refload] reading " << paths.size() << " reference FASTAs ("
                  << "max " << o.maxRefMemoryMB << " MB)";
        if (nThreads > 1) std::cerr << " with " << nThreads << " threads";
        else              std::cerr << " serially";
        std::cerr << "\n";
    }

    struct ChunkResult {
        SimpleRefIndex chunk;
        size_t bytes = 0;
    };
    auto worker = [&](size_t lo, size_t hi, ChunkResult* out, std::atomic<size_t>* doneCounter,
                      std::atomic<size_t>* lastReportSec, std::chrono::steady_clock::time_point t0) {
        for (size_t i = lo; i < hi; ++i) {
            const std::string& fasta_path = paths[i];
            if (fasta_path.empty() || !fs::exists(fasta_path)) {
                ++(*doneCounter);
                continue;
            }
            std::string asmName = fs::path(fasta_path).stem().string();
            try {
                auto contigs = read_fasta(fasta_path);
                for (auto& kv : contigs) {
                    if (kv.first.find("__sv_") != std::string::npos) {
                        std::cerr << "[warn] reference contig '" << kv.first
                                  << "' in " << fasta_path
                                  << " carries a simulator hint suffix (__sv_)."
                                     " Reference genomes must have plain biological"
                                     " contig names. This contig will not be indexed.\n";
                        continue;
                    }
                    const int len = static_cast<int>(kv.second.size());
                    out->bytes += kv.second.size();
                    out->chunk[kv.first].push_back(RefContigInfo{asmName, kv.first, std::move(kv.second), len});
                }
            } catch (const std::exception& e) {
                std::cerr << "[warn] cannot load reference FASTA " << fasta_path
                          << ": " << e.what() << '\n';
            }
            const size_t done = ++(*doneCounter);
            if (!o.quiet) {
                const auto now = std::chrono::steady_clock::now();
                const size_t elapsed = static_cast<size_t>(
                    std::chrono::duration_cast<std::chrono::seconds>(now - t0).count());
                // Per-25-ref OR per-30-second progress, whichever comes first.
                bool emit = (done % 25 == 0) || (done == paths.size());
                if (!emit) {
                    size_t prev = lastReportSec->load(std::memory_order_relaxed);
                    if (elapsed >= prev + 30) {
                        if (lastReportSec->compare_exchange_strong(prev, elapsed)) emit = true;
                    }
                }
                if (emit) {
                    std::cerr << "[refload] " << done << "/" << paths.size()
                              << " refs loaded (" << elapsed << "s elapsed)\n";
                }
            }
        }
    };

    std::vector<ChunkResult> chunks(nThreads);
    std::vector<std::thread> threads;
    std::atomic<size_t> doneCounter{0};
    std::atomic<size_t> lastReportSec{0};
    const size_t per = (paths.size() + nThreads - 1) / nThreads;
    for (unsigned t = 0; t < nThreads; ++t) {
        const size_t lo = std::min<size_t>(t * per, paths.size());
        const size_t hi = std::min<size_t>(lo + per, paths.size());
        if (lo >= hi) break;
        threads.emplace_back(worker, lo, hi, &chunks[t], &doneCounter, &lastReportSec, t0);
    }
    for (auto& th : threads) th.join();

    // ── Serial merge with memory-cap enforcement. ────────────────────────
    size_t totalRefBytes = 0;
    bool refCapWarned = false;
    // Bug 4 fix: cap unique contig names retained in the flat-fallback
    // SimpleRefIndex. The per-FASTA cap upstream doesn't constrain the
    // load — a single fungal ref can ship hundreds of unplaced scaffolds
    // and a 8-FASTA cap can still expand to 6 000+ contigs at the flat
    // fallback's iteration unit.
    const size_t maxRefContigs = static_cast<size_t>(std::max(0, o.maxFlatRefContigs));
    bool refContigCapWarned = false;
    for (auto& c : chunks) {
        for (auto& kv : c.chunk) {
            for (auto& contig : kv.second) {
                if (maxRefBytes > 0 && totalRefBytes + contig.sequence.size() > maxRefBytes) {
                    if (!refCapWarned) {
                        std::cerr << "[warn] --max-ref-memory-mb cap (" << o.maxRefMemoryMB
                                  << " MB) reached; skipping remaining reference contigs\n";
                        refCapWarned = true;
                    }
                    continue;
                }
                if (maxRefContigs > 0 && idx.size() >= maxRefContigs &&
                    idx.find(kv.first) == idx.end()) {
                    if (!refContigCapWarned) {
                        std::cerr << "[warn] --max-flat-ref-contigs cap ("
                                  << o.maxFlatRefContigs
                                  << ") reached; skipping remaining unique reference contigs\n";
                        refContigCapWarned = true;
                    }
                    continue;
                }
                totalRefBytes += contig.sequence.size();
                idx[kv.first].push_back(std::move(contig));
            }
        }
    }

    if (!o.quiet) {
        const auto t1 = std::chrono::steady_clock::now();
        const auto sec = std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
        std::cerr << "[refload] done: " << idx.size() << " unique contig names, "
                  << (totalRefBytes / (1024 * 1024)) << " MB, "
                  << sec << "s\n";
    }
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
    // Keep overlap thresholds and edge-case handling single-sourced.
    return tol::kmer_overlap_fraction(a, b, k);
}

static std::string novelty_tier_for_overlap(double frac) {
    // Delegates to the canonical tol::infer_novelty_tier so thresholds
    // are defined in exactly one place (hierarchical_engine.hpp).
    return tol::infer_novelty_tier(frac);
}

static bool is_low_complexity_sequence(const std::string& seq) {
    // Keep read-mode and assembly fallback complexity filters aligned.
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
    // ASSEMBLY: k=9 saturated the 4^9=262 k k-mer space for any multi-Mb
    // fungal contig, so containment_with_ref_hashes returned ~1.0 for EVERY
    // reference contig — the prefilter could not tell a true syntenic
    // homolog from an unrelated genome, and ref selection (perRefRanked /
    // select_mem_chain_refs) degenerated to an arbitrary partial_sort tie
    // break. Query chromosomes got chained against the wrong genomes and
    // the per-sample call volume collapsed. Assembly contigs are error-free,
    // so a large k is both safe and discriminative: at k=16 a genus-level
    // homolog still scores ~0.05–0.4 containment while non-homologs fall to
    // ~0, restoring meaningful ranking. Short assembly fragments keep enough
    // 16-mers to rank against; sub-16 bp contigs yield no hashes either way.
    return 16;
}

static int dynamic_calls_per_contig_cap(const Options& o,
                                        query_input::QueryMode mode,
                                        size_t contigLen) {
    const int base = std::max(1, o.maxCallsPerContig);
    if (mode != query_input::QueryMode::ASSEMBLY) return base;
    // Density-aware safety cap. Real fungal SV densities (Fg, Zt, Fo
    // pangenomes) reach ~1 SV per 1–2 kb in TE-rich accessory compartments
    // and one per 5–10 kb in core compartments. Scale at the high-density
    // end (1 / 1500 bp) so chromosome-sized contigs aren't artificially
    // throttled.
    const size_t byLength = std::max<size_t>(static_cast<size_t>(base),
                                              contigLen / 1500);
    return static_cast<int>(std::min<size_t>(50000, std::max<size_t>(base, byLength)));
}

static std::vector<const RefContigInfo*>
select_mem_chain_refs(const std::string& contigName,
                      const std::string& seq,
                      const SimpleRefIndex& refIdx,
                      RefSearchCache& cache,
                      query_input::QueryMode mode,
                      size_t maxShortlist = 12) {
    std::vector<const RefContigInfo*> out;
    const std::string baseName = strip_sv_suffix(contigName);
    std::unordered_set<std::string> matchingAsms;
    auto it = refIdx.find(baseName);
    if (it != refIdx.end()) {
        out.reserve(it->second.size());
        for (const auto& cand : it->second)
            if (!cand.sequence.empty()) {
                out.push_back(&cand);
                matchingAsms.insert(cand.asmName);
            }
    }

    // Augment shortlist with the other contigs of each assembly that contains
    // a same-named contig. Required for inter-contig TRA detection: with only
    // same-named contigs in the SA, MEM chains never see the cross-contig
    // anchor needed by SvTypeFromChain to emit a TRA.
    if (!matchingAsms.empty()) {
        for (const auto& bucket : refIdx) {
            if (bucket.first == baseName) continue;
            for (const auto& cand : bucket.second) {
                if (cand.sequence.empty()) continue;
                if (matchingAsms.count(cand.asmName))
                    out.push_back(&cand);
            }
        }
        return out;
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
    // Bug-fix (2026-05-12): the old guard only fired at bestOverlap ≥ 0.20,
    // which left a window 0.05 ≤ overlap < 0.20 where Pass 2 would happily
    // emit DEL = (qLen − bestRefLen) as a single giant deletion the size of
    // the entire ref contig. We were seeing 596 kb / 929 kb "DEL" calls
    // produced this way — they are not deletions, they are pseudo-contigs
    // covering only a sliver of a much longer reference contig. At low
    // overlap, also treat the pseudocontig as a fragment whenever it is
    // far shorter than the matched reference. The 0.55 ratio matches the
    // high-overlap branch above for consistency.
    if (bestRefLen > 0 && bestOverlap > 0.0 && bestOverlap < 0.20) {
        if (static_cast<long long>(qLen) * 100LL
                < static_cast<long long>(bestRefLen) * 55LL) {
            return true;
        }
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
        // Same rationale as reads_mode_kmer_fallback above: the length-only
        // path picked a ref contig from a name match without running an
        // alignment, so the position is a fixed placeholder
        // (segmentLen/4 + 1). Setting refPos / refEnd to 0 lets the VCF
        // emitter skip those INFO fields and the benchmark loader exclude
        // the call from reference-coord scoring rather than report a
        // breakpoint at a non-real position. The call is still emitted in
        // query-coord; downstream callers can promote it back if a real
        // alignment becomes available.
        v.refPos       = 0;
        v.refEnd       = 0;
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
    // Only pathologically fragmented draft assemblies (>1000 contigs) trip
    // the per-fragment guard, and even then we skip just the truly tiny
    // (<10 kb) contigs. The original >100-contig / <1 Mb rule discarded
    // EVERY sub-Mb contig of a hybrid assembly — e.g. GCA_030512215 has 26
    // contigs between 10 kb and 1 Mb carrying real genomic sequence that
    // were all silently dropped. With the k=16 prefilter each fragment now
    // matches only its true homolog, so the per-fragment SA cost the old
    // guard worried about is naturally bounded; the guard need only exclude
    // sub-10 kb scaffolds that cannot carry meaningful chain breakpoints.
    const bool fragmentedAssembly =
        mode == query_input::QueryMode::ASSEMBLY && contigs.size() > 1000;
    for (const auto& kv : contigs) {
        const std::string& contigName = kv.first;
        const std::string& seq = kv.second;
        if (static_cast<int>(seq.size()) < o.minSvLen) continue;
        if (is_low_complexity_sequence(seq)) continue;
        if (fragmentedAssembly && seq.size() < 10000) {
            // Sub-10 kb fragments of a >1000-contig draft assembly: the
            // length / off-reference fallbacks below still report coarse
            // structural signal, but an exact-MEM chain on a fragment this
            // short produces noise, not breakpoint evidence.
            continue;
        }

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
        // ALGORITHMIC FIX (per-ref pairwise + non-waterfall):
        //   Previously two short-circuits dropped most of the SV burden:
        //     (a) the multi-ref MEM chain ran ONCE over the top-K refs
        //         concatenated into a single suffix-array text; the chain
        //         mosaic across refs is NOT the same as per-ref pairwise
        //         calls, and one query-vs-best-mosaic loses the SVs that
        //         only exist relative to OTHER refs.
        //     (b) if the multi-ref path emitted ANY call, the single-ref
        //         cached path was skipped (`continue`).
        //   The literature standard (svim_asm, minigraph, SyRI) is per-ref
        //   pairwise. We now run BOTH the multi-ref pangenome chain AND a
        //   per-ref cached chain against every reference whose prefilter
        //   k-mer containment exceeds a modest floor. Calls union and are
        //   later type+svlen-deduplicated by select_best_call_per_contig.
        const auto shortlist = select_mem_chain_refs(
            contigName, seq, refIdx, cache, mode,
            (mode == query_input::QueryMode::ASSEMBLY) ? 24u : 32u);
        if (!shortlist.empty()) {
            const auto refBundle = make_refseq_bundle(shortlist);
            std::vector<VariantCallBridge> chainCalls;
            if (tol::try_mem_chain_call_multi_public(qAsm, contigName, seq,
                                                     refBundle.ptrs, eff_fo, chainCalls)) {
                for (auto& chainCall : chainCalls)
                    out.push_back(std::move(chainCall));
            }
        }

        // Per-ref cached chain over the top-K prefilter-positive references.
        // The SingleRefMemCache amortises SA cost across queries; the
        // top-K cap is needed for million_real where N_refs can exceed
        // 1000 and an unbounded loop blows up runtime even with caching.
        //
        // Per-ref selection floor + cap.
        //
        // History: the 2026-05-14 "performance fix" raised the floor to 0.10
        // and the cap to 8 to bound runtime. But with the old k=9 prefilter
        // every ref scored ~1.0, so 0.10 filtered nothing and the cap of 8
        // just kept 8 ARBITRARY refs — the genuine homolog was routinely not
        // among them, which is what collapsed per-sample call volume.
        //
        // Now that reads_mode_overlap_k uses k=16 for assembly, containment
        // is discriminative: a true genus-level homolog scores ~0.05–0.4
        // while unrelated genomes fall to ~0. A small floor (0.02) therefore
        // excludes non-homologs on its own — the cap is only a safety net.
        //
        // Per-assembly diversity: with the broadened bench_ref_list (hundreds
        // of ref ASSEMBLIES × several contigs each), a global
        // top-K is heavily skewed toward whichever ref has the most homologous
        // chromosomes. We first reserve N slots per ref ASSEMBLY (best contig
        // each) so every loaded ref gets a chance, then fill remaining slots
        // with the global top by overlap. Without this, sister-species refs
        // can shadow more-distant refs entirely and queries from sparser
        // genera (Lodderomyces, Dactylellina — 6-7 corpus refs each) end up
        // with 0 contigs from their own genus in the per-ref pairwise loop.
        constexpr double kPerRefMinOverlap = 0.02;
        constexpr size_t kMaxPerRefRefs   = 256;
        constexpr size_t kPerAsmReserved  = 3;
        struct PerRefCand { size_t idx; double frac; };
        std::vector<PerRefCand> perRefRanked;
        perRefRanked.reserve(cache.flatRefs.size());
        for (size_t i = 0; i < cache.flatRefs.size(); ++i) {
            const double frac = cache.containment_with_ref_hashes(
                prefilterHashes, i, prefilterK);
            if (frac < kPerRefMinOverlap) continue;
            perRefRanked.push_back({i, frac});
        }
        if (perRefRanked.size() > kMaxPerRefRefs) {
            std::sort(perRefRanked.begin(), perRefRanked.end(),
                      [](const PerRefCand& a, const PerRefCand& b) {
                          return a.frac > b.frac;
                      });
            // Reserve slots per ref ASSEMBLY in overlap order.
            std::unordered_map<std::string, size_t> kept_per_asm;
            std::vector<PerRefCand> diverse;
            std::vector<PerRefCand> remainder;
            diverse.reserve(kMaxPerRefRefs);
            remainder.reserve(perRefRanked.size());
            for (const auto& pr : perRefRanked) {
                const RefContigInfo* info = cache.flatRefs[pr.idx].info;
                if (info == nullptr) { remainder.push_back(pr); continue; }
                if (kept_per_asm[info->asmName] < kPerAsmReserved) {
                    diverse.push_back(pr);
                    ++kept_per_asm[info->asmName];
                } else {
                    remainder.push_back(pr);
                }
            }
            for (const auto& pr : remainder) {
                if (diverse.size() >= kMaxPerRefRefs) break;
                diverse.push_back(pr);
            }
            if (diverse.size() > kMaxPerRefRefs) diverse.resize(kMaxPerRefRefs);
            perRefRanked = std::move(diverse);
        }
        for (const auto& pr : perRefRanked) {
            const RefContigInfo* refInfo = cache.flatRefs[pr.idx].info;
            if (refInfo == nullptr) continue;
            std::vector<VariantCallBridge> perRefCalls;
            if (try_mem_chain_call_single_ref_cached_multi(
                    qAsm, contigName, seq, refInfo, memCache,
                    eff_fo, perRefCalls)) {
                for (auto& c : perRefCalls) out.push_back(std::move(c));
            }
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
        // DUP fallback: ChainTreap requires strictly-increasing rPos and silently
        // drops backward-mapping MEMs that are the signature of tandem duplications.
        // With diverged genomes the primary chain classifies DUPs as INS.
        if (res.type == tol::SvTypeFromChain::Type::NONE ||
            res.type == tol::SvTypeFromChain::Type::INS) {
            std::vector<tol::SuffixArray::Mem> fwdSorted;
            fwdSorted.reserve(fwdMems.size());
            for (int i : order)
                if (!isRev[static_cast<size_t>(i)])
                    fwdSorted.push_back(allMems[static_cast<size_t>(i)]);
            if (fwdSorted.size() >= 2) {
                auto dupRes = tol::SvTypeFromChain::classify(
                    fwdSorted, std::vector<bool>(fwdSorted.size(), false),
                    sa, fo.minSvLen);
                if (dupRes.type == tol::SvTypeFromChain::Type::DUP) {
                    double dupScore = 0.0;
                    for (const auto& m : fwdSorted) dupScore += static_cast<double>(m.len);
                    if (dupScore >= fo.minBlockScore) {
                        res      = dupRes;
                        chain    = fwdSorted;
                        chainRev.assign(fwdSorted.size(), false);
                    }
                }
            }
        }
        if (res.type == tol::SvTypeFromChain::Type::NONE) return out;

        out.call.qAsm = qAsm;
        out.call.qContig = qContig;
        out.call.refAsm = primaryRef.asmName.empty() ? "unknown" : primaryRef.asmName;
        out.call.refContig = primaryRef.contig.empty() ? "." : primaryRef.contig;
        out.call.refPos = 0;
        out.call.refEnd = 0;
        out.call.pos = std::max(1, res.qBreakStart + 1);
        out.call.end = std::max(out.call.pos, res.qBreakEnd >= 0 ? res.qBreakEnd + 1 : out.call.pos);
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
                out.call.refPos = res.rBreakStart >= 0 ? (res.rBreakStart + 1) : 0;
                out.call.refEnd = out.call.refPos;
                break;
            case T::DEL:
                out.call.type = "DEL";
                out.call.pantreeClass = "DEL";
                out.call.refPos = res.rBreakStart >= 0 ? (res.rBreakStart + 1) : 0;
                out.call.refEnd = res.rBreakEnd > 0 ? res.rBreakEnd : out.call.refPos;
                break;
            case T::INV:
                out.call.type = "INV";
                out.call.pantreeClass = "INV";
                out.call.refPos = res.rBreakStart >= 0 ? (res.rBreakStart + 1) : 0;
                out.call.refEnd = res.rBreakEnd > 0 ? res.rBreakEnd : out.call.refPos;
                break;
            case T::DUP:
                out.call.type = "DUP";
                out.call.pantreeClass = "DUP";
                out.call.refPos = res.rBreakStart >= 0 ? (res.rBreakStart + 1) : 0;
                out.call.refEnd = res.rBreakEnd > 0 ? res.rBreakEnd : out.call.refPos;
                break;
            case T::TRA:
                out.call.type = "TRA";
                out.call.pantreeClass = "NON_REF";
                {
                    const int srcRPos = !chain.empty()
                        ? (chain.front().rPos + chain.front().len) : 0;
                    int srcCi = -1;
                    for (int ci = 0; ci < static_cast<int>(sa.contigEnd.size()); ++ci) {
                        if (srcRPos < sa.contigEnd[static_cast<size_t>(ci)]) {
                            srcCi = ci;
                            break;
                        }
                    }
                    if (srcCi < 0 && !sa.contigEnd.empty())
                        srcCi = static_cast<int>(sa.contigEnd.size()) - 1;
                    const int srcOff = (srcCi > 0)
                        ? sa.contigEnd[static_cast<size_t>(srcCi) - 1] : 0;
                    out.call.refPos = srcRPos > 0 ? (srcRPos - srcOff + 1) : 0;
                    out.call.refEnd = out.call.refPos;
                }
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

// Multi-emit chain caller (assembly + reads-mode INS/DEL recall fix).
// Mirrors try_mem_chain_call_single_ref_cached but uses SvTypeFromChain::classify_all
// to yield every per-pair INS/DEL gap instead of only the chain's dominant event.
// On real diverged fungal genomes a single chain typically contains 10–100 small SVs
// that the single-emit path collapsed into one call, driving recall to ~5%.
static bool try_mem_chain_call_single_ref_cached_multi(
        const std::string& qAsm,
        const std::string& qContig,
        const std::string& qSeq,
        const RefContigInfo* refInfo,
        SingleRefMemCache& memCache,
        const tol::FederatedOptions& fo,
        std::vector<VariantCallBridge>& calls) {
    calls.clear();
    if (refInfo == nullptr || refInfo->sequence.empty() || qSeq.empty()) return false;

    SingleRefMemIndex& idx = memCache.get(refInfo);
    const tol::SuffixArray& sa = idx.sa;
    const tol::TolGlobal::RefSeq& primaryRef = idx.ref;
    const std::string rcSeq = tol::SuffixArray::revcomp(qSeq);

    auto min_mem_from_k = [](int k) {
        return std::max(15, k - 5);
    };

    struct MultiChainAttempt {
        std::vector<VariantCallBridge> calls;
        double score   = 0.0;
        int    anchors = 0;
        bool   valid   = false;
    };

    auto attempt_chain = [&](int minMem, bool secondaryPass) {
        MultiChainAttempt out;
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

        std::unordered_map<uint64_t, size_t> posToMemIdx;
        posToMemIdx.reserve(allMems.size());
        for (size_t mi = 0; mi < allMems.size(); ++mi) {
            uint64_t key = (static_cast<uint64_t>(allMems[mi].qPos) << 32) |
                           static_cast<uint64_t>(static_cast<uint32_t>(allMems[mi].rPos));
            posToMemIdx.emplace(key, mi);
        }

        // Fix 4: short pseudocontigs (typical sr_unitig in short-reads mode,
        // lr_pcN clusters in long-reads mode) often carry one strong MEM
        // anchor plus weak/divergent flanks — the chain.size() >= 2 +
        // minBlockScore = 6 gate forced 98 % of them to OFF_REF. Soften both
        // gates when qSeq is short so the chain caller can fire on a single
        // good anchor.
        const bool shortPseudocontig = qSeq.size() < 1000;
        const double scoreFloor = shortPseudocontig
            ? std::min(fo.minBlockScore, 3.0)
            : fo.minBlockScore;
        const size_t minChainAnchors = shortPseudocontig ? 1u : 2u;

        // Fix 1: iterative multi-chain extraction. A chromosome-length contig
        // typically has 5–50 independent syntenic blocks separated by SVs
        // (each block is its own MEM chain because the gap between blocks
        // exceeds chainGapBand=5 kbp). The old code took only
        // treap.best_chain_path() — i.e. one block per (qcontig, ref) pair —
        // so a sample with 250 SVs scattered across the genome got reduced
        // to ~17 emitted calls. Now we extract chains best-first, mark their
        // MEMs as consumed, rebuild the treap from the rest, and repeat.
        //
        // Cap scales with query length at ~1 chain per 50 kb (matches the
        // typical fungal syntenic-block spacing), capped at 1024 for safety.
        const int maxChainsPerCall = std::max(
            64,
            std::min(1024, static_cast<int>(qSeq.size() / 50000)));
        const int maxGap = fo.chainGapBand > 0 ? fo.chainGapBand : 5000;
        std::vector<bool> memConsumed(allMems.size(), false);
        std::vector<tol::SvTypeFromChain::Result> aggregatedEvents;
        double topChainScore = 0.0;
        int    topChainAnchors = 0;
        std::vector<tol::SuffixArray::Mem> topChain;
        std::vector<bool> topChainRev;

        for (int chainIter = 0; chainIter < maxChainsPerCall; ++chainIter) {
            tol::ChainTreap treap;
            int memsLeft = 0;
            for (int i : order) {
                if (memConsumed[static_cast<size_t>(i)]) continue;
                ++memsLeft;
                const auto& m = allMems[static_cast<size_t>(i)];
                treap.insert_and_chain(m.qPos, m.rPos, m.len,
                                       static_cast<float>(m.len), maxGap);
            }
            if (memsLeft == 0) break;

            auto chainIdx = treap.best_chain_path();
            if (chainIdx.empty()) break;

            const double bestScore = static_cast<double>(treap.best_chain_score());
            if (bestScore < scoreFloor) break;

            std::vector<tol::SuffixArray::Mem> chain;
            std::vector<bool> chainRev;
            std::vector<size_t> chainMemIdx;
            chain.reserve(chainIdx.size());
            chainRev.reserve(chainIdx.size());
            chainMemIdx.reserve(chainIdx.size());
            for (int ni : chainIdx) {
                const auto& nd = treap.nodes_[static_cast<size_t>(ni)];
                uint64_t key = (static_cast<uint64_t>(nd.qPos) << 32) |
                               static_cast<uint64_t>(static_cast<uint32_t>(nd.rPos));
                auto it = posToMemIdx.find(key);
                if (it != posToMemIdx.end()) {
                    chain.push_back(allMems[it->second]);
                    chainRev.push_back(isRev[it->second]);
                    chainMemIdx.push_back(it->second);
                }
            }

            // ALGORITHMIC FIX (recall): treap.best_chain_path() returns the
            // highest-SCORING chain, not the longest. A single long MEM
            // (score = len) routinely outranks a multi-anchor chain of
            // shorter MEMs. The old code did `break` here, so as soon as
            // the top-scoring residual chain had < minChainAnchors anchors
            // we abandoned ALL further iterations — including any genuine
            // multi-anchor chains hiding behind the long singleton. Mark
            // the offending MEMs consumed and `continue` instead: the next
            // iteration's treap rebuild will then surface the multi-anchor
            // chains we actually want. The `chainIdx.empty()` and
            // `bestScore < scoreFloor` breaks above already guarantee
            // forward progress (the loop bound also fences against runaway).
            if (chain.size() < minChainAnchors) {
                for (size_t mi : chainMemIdx) memConsumed[mi] = true;
                continue;
            }

            // Mark this chain's MEMs as consumed so the next iteration
            // searches the remaining genome instead of re-extracting the
            // same chain.
            for (size_t mi : chainMemIdx) memConsumed[mi] = true;

            auto events = tol::SvTypeFromChain::classify_all(
                chain, chainRev, sa, fo.minSvLen);

            // DUP fallback applies only on the first chain — the rescue is
            // about catching a tandem-dup pattern that ChainTreap silently
            // drops because it requires strictly-increasing rPos. After we
            // start iterating, the previously consumed MEMs would distort
            // this signal.
            if (chainIter == 0) {
                const bool needsDupRescue = events.empty() ||
                    (events.size() == 1 &&
                     (events.front().type == tol::SvTypeFromChain::Type::INS));
                if (needsDupRescue) {
                    std::vector<tol::SuffixArray::Mem> fwdSorted;
                    fwdSorted.reserve(fwdMems.size());
                    for (int i : order)
                        if (!isRev[static_cast<size_t>(i)])
                            fwdSorted.push_back(allMems[static_cast<size_t>(i)]);
                    if (fwdSorted.size() >= 2) {
                        auto dupRes = tol::SvTypeFromChain::classify(
                            fwdSorted, std::vector<bool>(fwdSorted.size(), false),
                            sa, fo.minSvLen);
                        if (dupRes.type == tol::SvTypeFromChain::Type::DUP) {
                            double dupScore = 0.0;
                            for (const auto& m : fwdSorted) dupScore += static_cast<double>(m.len);
                            if (dupScore >= fo.minBlockScore) {
                                events.assign(1, dupRes);
                                chain    = fwdSorted;
                                chainRev.assign(fwdSorted.size(), false);
                            }
                        }
                    }
                }
            }
            if (events.empty()) continue;

            // ALGORITHMIC FIX (recall): capture topChain on the first chain
            // that PRODUCES EVENTS, not strictly chainIter==0. The old check
            // would leave topChain empty if iteration 0 yielded no events
            // (now also possible because the <minChainAnchors continue we
            // added above can skip iteration 0). With topChain empty,
            // fill_common would emit every aggregatedEvent with
            // blockScore=0 and anchors=0 — and select_best_call_per_contig
            // would drop them all via its blockScore < minBlockScore (=6)
            // gate, silently throwing away every event from this query
            // contig × ref pair.
            if (topChain.empty()) {
                topChainScore   = bestScore;
                topChainAnchors = static_cast<int>(chain.size());
                topChain        = chain;
                topChainRev     = chainRev;
            }
            for (auto& ev : events) aggregatedEvents.push_back(std::move(ev));
        }

        if (aggregatedEvents.empty()) return out;

        // Sub-event aggregation: when classify_all emits many tiny INS/DEL within
        // a few hundred bp of each other (typical of microsatellite mis-anchoring),
        // merging them into a single representative call reduces FP rate without
        // hurting recall — the comparator truth set also collapses these.
        // Decoupled from minSvLen — see fungi_tol_bridge.hpp:mergeWindow note.
        const int mergeWindow = 80;
        std::sort(aggregatedEvents.begin(), aggregatedEvents.end(),
                  [](const tol::SvTypeFromChain::Result& a,
                     const tol::SvTypeFromChain::Result& b) {
                      if (a.qBreakStart != b.qBreakStart) return a.qBreakStart < b.qBreakStart;
                      return a.svLen > b.svLen;
                  });
        std::vector<tol::SvTypeFromChain::Result> merged;
        merged.reserve(aggregatedEvents.size());
        for (auto& ev : aggregatedEvents) {
            if (!merged.empty()) {
                auto& last = merged.back();
                if (last.type == ev.type &&
                    last.rContig == ev.rContig &&
                    std::abs(ev.qBreakStart - last.qBreakEnd) <= mergeWindow) {
                    last.qBreakEnd = std::max(last.qBreakEnd, ev.qBreakEnd);
                    last.rBreakEnd = std::max(last.rBreakEnd, ev.rBreakEnd);
                    last.svLen     = std::max(last.svLen, ev.svLen);
                    continue;
                }
            }
            merged.push_back(std::move(ev));
        }

        // Use the top chain's metadata for the call header (alignment mode,
        // score, anchor count). Per-event chain identity isn't currently
        // recorded; downstream metrics only use the per-VCF aggregate.
        const double bestScore = topChainScore;
        std::vector<tol::SuffixArray::Mem>& chain = topChain;
        std::vector<bool>& chainRev = topChainRev;
        (void)chainRev;  // silence unused warning when no TRA emitted

        auto fill_common = [&](VariantCallBridge& v) {
            v.qAsm = qAsm;
            v.qContig = qContig;
            v.refAsm = primaryRef.asmName.empty() ? "unknown" : primaryRef.asmName;
            v.refContig = primaryRef.contig.empty() ? "." : primaryRef.contig;
            v.refPos = 0;
            v.refEnd = 0;
            v.genotype = "0/1";
            v.gq = 40.0;
            v.blockScore = bestScore;
            v.anchors = static_cast<int>(chain.size());
            v.alignmentMode = secondaryPass
                ? "mem_chain_cached_single_ref_multi;secondary_seed_rescue"
                : "mem_chain_cached_single_ref_multi";
            v.mapq = 50.0;
            v.annotation = "NONE";
            v.triallelicTopology = ".";
            v.isNonRefVariant = false;
            v.cladeRank = ".";
            v.phylum = ".";
        };

        using T = tol::SvTypeFromChain::Type;
        for (const auto& res : merged) {
            if (res.type == T::NONE) continue;
            VariantCallBridge v;
            fill_common(v);
            v.pos = std::max(1, res.qBreakStart + 1);
            v.end = std::max(v.pos, res.qBreakEnd >= 0 ? res.qBreakEnd + 1 : v.pos);
            v.svlen = res.svLen;
            switch (res.type) {
                case T::INS:
                    v.type = "INS"; v.pantreeClass = "INS";
                    v.refPos = res.rBreakStart >= 0 ? (res.rBreakStart + 1) : 0;
                    v.refEnd = v.refPos;
                    break;
                case T::DEL:
                    v.type = "DEL"; v.pantreeClass = "DEL";
                    v.refPos = res.rBreakStart >= 0 ? (res.rBreakStart + 1) : 0;
                    v.refEnd = res.rBreakEnd > 0 ? res.rBreakEnd : v.refPos;
                    break;
                case T::INV:
                    v.type = "INV"; v.pantreeClass = "INV";
                    v.refPos = res.rBreakStart >= 0 ? (res.rBreakStart + 1) : 0;
                    v.refEnd = res.rBreakEnd > 0 ? res.rBreakEnd : v.refPos;
                    break;
                case T::DUP:
                    v.type = "DUP"; v.pantreeClass = "DUP";
                    v.refPos = res.rBreakStart >= 0 ? (res.rBreakStart + 1) : 0;
                    v.refEnd = res.rBreakEnd > 0 ? res.rBreakEnd : v.refPos;
                    break;
                case T::TRA:
                    v.type = "TRA"; v.pantreeClass = "NON_REF";
                    {
                        const int srcRPos = !chain.empty()
                            ? (chain.front().rPos + chain.front().len) : 0;
                        int srcCi = -1;
                        for (int ci = 0; ci < static_cast<int>(sa.contigEnd.size()); ++ci) {
                            if (srcRPos < sa.contigEnd[static_cast<size_t>(ci)]) {
                                srcCi = ci; break;
                            }
                        }
                        if (srcCi < 0 && !sa.contigEnd.empty())
                            srcCi = static_cast<int>(sa.contigEnd.size()) - 1;
                        const int srcOff = (srcCi > 0)
                            ? sa.contigEnd[static_cast<size_t>(srcCi) - 1] : 0;
                        v.refPos = srcRPos > 0 ? (srcRPos - srcOff + 1) : 0;
                        v.refEnd = v.refPos;
                    }
                    v.mateContig = res.rContig.empty() ? primaryRef.contig : res.rContig;
                    v.matePos = res.rBreakStart + 1;
                    v.mateEnd = res.rBreakEnd > 0 ? res.rBreakEnd : v.matePos;
                    v.mateRefAsm = primaryRef.asmName.empty() ? "." : primaryRef.asmName;
                    v.mateOffReference = false;
                    break;
                default:
                    continue;
            }
            out.calls.push_back(std::move(v));
        }

        if (out.calls.empty()) return out;
        out.score   = bestScore;
        out.anchors = static_cast<int>(chain.size());
        out.valid   = true;
        return out;
    };

    const int primaryMinMem = min_mem_from_k(fo.primarySketchParams.k);
    MultiChainAttempt best = attempt_chain(primaryMinMem, false);
    if (fo.useSecondarySeeds) {
        const int secondaryMinMem = min_mem_from_k(fo.secondarySketchParams.k);
        const bool rescueRequested = !best.valid ||
            best.anchors < static_cast<int>(std::max<size_t>(fo.repeatRescueMinAnchors, 2));
        if (secondaryMinMem < primaryMinMem && rescueRequested) {
            MultiChainAttempt rescue = attempt_chain(secondaryMinMem, true);
            if (rescue.valid &&
                (!best.valid ||
                 rescue.calls.size() > best.calls.size() ||
                 rescue.anchors > best.anchors ||
                 rescue.score   > best.score)) {
                best = std::move(rescue);
            }
        }
    }
    if (!best.valid) return false;
    calls = std::move(best.calls);
    return !calls.empty();
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
    const bool fragmentedAssembly =
        mode == query_input::QueryMode::ASSEMBLY && contigs.size() > 100;
    for (const auto& kv : contigs) {
        const std::string& contigName = kv.first;
        const std::string& seq = kv.second;
        if (static_cast<int>(seq.size()) < o.minSvLen) continue;
        if (is_low_complexity_sequence(seq)) continue;
        if (fragmentedAssembly && seq.size() < 1000000) {
            // Whole-contig novelty scans are not meaningful for highly
            // fragmented draft assemblies: hundreds of short contigs would
            // each be compared to every benchmark reference contig, dominating
            // runtime while producing coarse OFF_REF calls. Length fallback
            // still captures obvious PAV/indel-scale differences.
            continue;
        }

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
        // Emit all four tiers including OFF_REF_KNOWN.
        if (tier != "NOVEL" && tier != "NOVEL_WEAK" && tier != "DIVERGED" &&
            tier != "OFF_REF_KNOWN") continue;

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
        // ALGORITHMIC FIX (per-ref pairwise + non-waterfall) — mirrors the
        // assembly path. Run BOTH the multi-ref pangenome MEM chain AND a
        // per-ref cached MEM chain against every prefilter-positive
        // reference, AND the length-delta Pass 2 below. No `continue`
        // short-circuits between passes.
        const auto shortlist = select_mem_chain_refs(contigName, seq, refIdx, cache, mode, 32u);
        if (!shortlist.empty()) {
            const auto refBundle = make_refseq_bundle(shortlist);
            std::vector<VariantCallBridge> chainCalls;
            if (tol::try_mem_chain_call_multi_public(qAsm, contigName, seq,
                                                     refBundle.ptrs, eff_fo, chainCalls)) {
                for (auto& chainCall : chainCalls) out.push_back(std::move(chainCall));
            }
        }
        // Per-assembly diversity (see assembly-mode comment above).
        constexpr double kPerRefMinOverlap = 0.02;
        constexpr size_t kMaxPerRefRefs   = 160;
        constexpr size_t kPerAsmReserved  = 3;
        struct PerRefCand { size_t idx; double frac; };
        std::vector<PerRefCand> perRefRanked;
        perRefRanked.reserve(cache.flatRefs.size());
        for (size_t i = 0; i < cache.flatRefs.size(); ++i) {
            const double frac = cache.containment_with_ref_hashes(
                prefilterHashes, i, prefilterK);
            if (frac < kPerRefMinOverlap) continue;
            perRefRanked.push_back({i, frac});
        }
        if (perRefRanked.size() > kMaxPerRefRefs) {
            std::sort(perRefRanked.begin(), perRefRanked.end(),
                      [](const PerRefCand& a, const PerRefCand& b) {
                          return a.frac > b.frac;
                      });
            std::unordered_map<std::string, size_t> kept_per_asm;
            std::vector<PerRefCand> diverse;
            std::vector<PerRefCand> remainder;
            diverse.reserve(kMaxPerRefRefs);
            remainder.reserve(perRefRanked.size());
            for (const auto& pr : perRefRanked) {
                const RefContigInfo* info = cache.flatRefs[pr.idx].info;
                if (info == nullptr) { remainder.push_back(pr); continue; }
                if (kept_per_asm[info->asmName] < kPerAsmReserved) {
                    diverse.push_back(pr);
                    ++kept_per_asm[info->asmName];
                } else {
                    remainder.push_back(pr);
                }
            }
            for (const auto& pr : remainder) {
                if (diverse.size() >= kMaxPerRefRefs) break;
                diverse.push_back(pr);
            }
            if (diverse.size() > kMaxPerRefRefs) diverse.resize(kMaxPerRefRefs);
            perRefRanked = std::move(diverse);
        }
        for (const auto& pr : perRefRanked) {
            const RefContigInfo* refInfo = cache.flatRefs[pr.idx].info;
            if (refInfo == nullptr) continue;
            std::vector<VariantCallBridge> perRefCalls;
            if (try_mem_chain_call_single_ref_cached_multi(
                    qAsm, contigName, seq, refInfo, memCache,
                    eff_fo, perRefCalls)) {
                for (auto& c : perRefCalls) out.push_back(std::move(c));
            }
        }

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

        // Plausibility guard for Pass 2: even when the pair survives
        // is_reads_mode_fragment, |delta| should fit within the part of the
        // ref the pseudocontig could actually shadow. Cap by min(qlen,
        // bestRefLen): a 1.4 kb pseudocontig cannot legitimately anchor a
        // 950 kb DEL — that's a partial overlap, not a deletion. Without
        // this, the kmer-fallback path emits one giant DEL per fragmentary
        // pseudocontig and dwarfs the real call set.
        const int spanLimit =
            (mode == query_input::QueryMode::ASSEMBLY)
                ? std::numeric_limits<int>::max()
                : std::max(o.minSvLen, std::min(qlen, bestRefLen));
        if (locusSized &&
            bestOverlap >= 0.05 && std::abs(delta) >= o.minSvLen
                                 && std::abs(delta) <= o.maxSvLen
                                 && std::abs(delta) <= spanLimit) {
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
            // Reference coordinates are NOT meaningful in the kmer-fallback
            // path — we picked the best-containment ref contig by k-mer
            // sketch but ran no alignment, so the midpoint is a placeholder
            // and would always miss a comparator's actual breakpoint. Leave
            // refPos/refEnd at 0 so write_vcf_record omits the REFPOS /
            // REFEND fields and load_mycosv_reference_calls excludes the
            // call from the reference-coord truth/pred set. The call still
            // appears in query-coord scoring and the multisample VCF.
            v.refPos        = 0;
            v.refEnd        = 0;
            v.pantreeClass  = v.type;
            v.isNonRefVariant    = false;
            v.triallelicTopology = ".";
            v.cladeRank = ".";
            v.phylum = ".";
            out.push_back(std::move(v));
            // Algorithmic fix: removed the `continue` that short-circuited
            // Pass 3 below. Pass 3 OFF_REF emission is handled by
            // simple_offref_fallback_calls downstream, but we no longer
            // skip the rest of the per-contig accumulator either way.
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
        // TRA with strong chain evidence (high blockScore + anchors) inherently has
        // low reference overlap because the breakpoint spans two different reference
        // contigs — penalising on overlap alone would suppress all inter-contig TRAs.
        const bool strongChainTra = (call.type == "TRA" &&
                                     call.blockScore >= 6.0 &&
                                     call.anchors >= 2);
        if ((call.type == "INV" || call.type == "DUP" ||
             (call.type == "TRA" && !strongChainTra)) && overlap < 0.05)
            ref += 0.75;
        if ((call.type == "INV" || call.type == "DUP" ||
             (call.type == "TRA" && !strongChainTra)) && overlap < 0.10) {
            alt -= 2.5;
            ref += 5.0 + (0.10 - overlap) * 20.0;
        }
        if (call.type == "TRA" && !strongChainTra &&
            call.mateContig.empty() && overlap < 0.10) {
            alt -= 1.0;
            ref += 2.0;
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
    // Strong-evidence TRA calls have inherently low posterior because reference
    // overlap is near zero for inter-contig events; skip the penalty for them.
    const bool strongChainTra2 = (call.type == "TRA" &&
                                   call.blockScore >= 6.0 &&
                                   call.anchors >= 2);
    if (fused.posteriorAlt < 0.35 && !strongChainTra2) score -= 120.0;
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

        // Algorithmic fix: dropped the DUP "looks-like-indel" -90/-140
        // penalty. In fungi, tandem duplications of effector clusters or
        // TE arrays *do* produce a length delta matching the duplicated
        // copy, and the previous code systematically downvoted them to
        // below the secondary floor, eliminating real DUP calls in
        // two-speed pathogens (Fusarium, Zymoseptoria).
    } else if (mode != query_input::QueryMode::ASSEMBLY) {
        // Reads-mode pseudo-contigs (sr_unitig*, lr_pc*) lack stable contig
        // names so name-based arbitration is unavailable. The penalties
        // below are appropriate for reads mode because near-zero overlap
        // with the assigned reference signals a chimeric assembly artifact
        // (genuine reads-mode SVs always anchor at one end). For ASSEMBLY
        // mode these don't fire (the if-branch above runs instead), so
        // they don't suppress real low-overlap fungal pangenome SVs.
        if (indel && overlap >= 0.04) score += 110.0;
        if (call.type == "TRA" && !call.mateContig.empty() && overlap >= 0.30) score += 115.0;
        if (largeRearr && overlap < 0.02) score -= 180.0;
        if (largeRearr && overlap < 0.05) score -= 70.0;
        if (call.type == "TRA" && overlap < 0.05 &&
            call.mateContig.empty() && !strongChainTra2) score -= 120.0;
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

        // ALGORITHMIC FIX (recall):
        //   The previous design picked ONE top winner per contig and then
        //   added up to (contigCap-1) "secondaries" gated by a hard
        //   secondaryFloor = bestSc.score - 200. On real fungal pangenomes
        //   (e.g. F. graminearum 52 420 SVs / panel, ≈500 per pairwise) most
        //   genuine SVs score well below the top pick and were silently
        //   dropped. The pipeline became architecturally bounded at
        //   O(maxCallsPerContig * n_contigs) regardless of true SV burden.
        //
        //   New behaviour: emit EVERY candidate whose blockScore meets the
        //   user-set minBlockScore floor and whose svlen/anchor counts are
        //   above the configured minima, with type-aware non-overlap so a
        //   co-located INS+DEL+DUP+OFF_REF set (typical at TE breakpoints)
        //   all survive. The single contigCap is retained only as a runaway
        //   safety net at 10× the previous default.
        struct ScoredCand {
            const VariantCallBridge* call = nullptr;
            tol::FusedEvidenceScore fused;
            double score = 0.0;
        };
        std::vector<ScoredCand> scored;
        scored.reserve(kv.second.size());

        for (const auto& cand : kv.second) {
            tol::FusedEvidenceScore fused;
            const double score = candidate_priority_score(
                cand, kv.second, qSeq, refIdx, cache, o, mode, report, &fused);
            scored.push_back({ &cand, fused, score });
        }
        if (scored.empty()) continue;

        std::sort(scored.begin(), scored.end(),
                  [](const ScoredCand& a, const ScoredCand& b) {
                      return a.score > b.score;
                  });

        // Safety cap retained but very large; the meaningful floor is now
        // minBlockScore + minAnchors on the candidate itself.
        const int contigCap = dynamic_calls_per_contig_cap(o, mode, qSeq.size());
        // mergeSlop kept very small — only collapse two truly identical
        // SVs at the same coordinate. SAME-TYPE adjacency is the only
        // legitimate fusion case; cross-type calls always coexist.
        const int mergeSlop = std::max(20, o.minSvLen / 2);

        // ── ALGORITHMIC FIX (precision): per-contig reference consolidation ──
        //   The hierarchical multi-ref search emits candidate calls against
        //   EVERY reference assembly a query contig touched. A length/off-ref
        //   fallback that a real chain alignment has already superseded, or a
        //   2-anchor `secondary_seed_rescue` stub on the SAME ref the chain
        //   already explains, is a search artifact. We drop those when a true
        //   chain has out-scored them by a wide margin.
        //
        //   We deliberately DO NOT drop genuine cross-reference chain calls
        //   here, even when they are far weaker than the dominant ref's chain.
        //   In a fungal multi-reference benchmark each query contig may carry
        //   real, distinct SVs against several related refs (sister-species
        //   pangenome; the same query position can be DEL vs ref_A and INS
        //   vs ref_B). The earlier 4× cross-ref drop collapsed per-query call
        //   counts ~10× on million_real because the strongest sister-species
        //   chain shadowed every other ref's signal.
        double dominantBlock = 0.0;
        for (const auto& sc : scored) {
            if (sc.call && sc.call->blockScore > dominantBlock) {
                dominantBlock = sc.call->blockScore;
            }
        }
        const double dominanceFactor = 4.0;
        auto is_simple_fallback = [](const std::string& m) {
            return m.rfind("simple_length_fallback", 0) == 0 ||
                   m.rfind("simple_offref_fallback", 0) == 0;
        };
        // `secondary_seed_rescue` is a recall-oriented second pass that emits
        // marginal 2–4 anchor stubs. Such a stub is legitimate only when no
        // primary chain explains the contig — once a real chain dominates
        // (4× block score), a far-weaker rescue stub is redundant noise even
        // when it lands on the dominant reference.
        auto is_secondary_rescue = [](const std::string& m) {
            return m.find("secondary_seed_rescue") != std::string::npos;
        };

        // Per-type kept intervals: an INS at pos X does not block a DEL
        // or OFF_REF at the same position (literature: TE breakpoints
        // routinely carry co-located insertions, deletions and novelty).
        struct KeptInterval { int pos; int end; int svlen; };
        std::unordered_map<std::string, std::vector<KeptInterval>> keptByType;

        // Score floor must be below the OFF_REF starting score (-250). The
        // real per-candidate gate is the minBlockScore check below, which
        // is exempt for OFF_REF (those have no chain-derived blockScore).
        const double scoreFloor = -400.0;

        int emittedThisContig = 0;
        for (const auto& sc : scored) {
            if (sc.call == nullptr) continue;
            if (sc.score < scoreFloor) continue;
            // Apply the user-set minimum-block-score floor as the real gate,
            // not a runtime-relative cutoff. Off-ref/novelty calls do not
            // have a chain-derived blockScore so are exempt.
            if (sc.call->type != "OFF_REF" &&
                sc.call->blockScore < o.minBlockScore) continue;

            // Reference-consolidation gate (see comment above): drop only the
            // search artifacts — simple_* fallbacks and secondary_seed_rescue
            // stubs — when they are far weaker than the dominant chain. Genuine
            // cross-reference chain calls are kept; a sister-species pangenome
            // legitimately surfaces distinct SVs against multiple refs.
            const bool farWeaker =
                sc.call->blockScore * dominanceFactor < dominantBlock;
            bool lengthDeltaSupportedIndel = false;
            if (is_simple_fallback(sc.call->alignmentMode) &&
                (sc.call->type == "INS" || sc.call->type == "DEL")) {
                const RefContigInfo* namedBest =
                    best_ref_match(refIdx, sc.call->qContig, static_cast<int>(qSeq.size()));
                if (namedBest != nullptr &&
                    sc.call->refAsm == namedBest->asmName &&
                    sc.call->refContig == namedBest->contigName) {
                    const int delta = static_cast<int>(qSeq.size()) - namedBest->length;
                    const int deltaAbs = std::abs(delta);
                    const int svAbs = std::abs(sc.call->svlen);
                    lengthDeltaSupportedIndel =
                        deltaAbs >= o.minSvLen &&
                        deltaAbs <= o.maxSvLen &&
                        std::abs(svAbs - deltaAbs) <= std::max(10, deltaAbs / 10);
                }
            }
            if (farWeaker &&
                (is_simple_fallback(sc.call->alignmentMode) ||
                 is_secondary_rescue(sc.call->alignmentMode)) &&
                !lengthDeltaSupportedIndel) {
                continue;
            }

            // ALGORITHMIC FIX (recall): bucket by (type, refAsm), not type
            // alone. A DEL of the query at position X relative to sister-
            // species ref_A and a DEL at the same X relative to a more
            // distant ref_B are GENUINELY DIFFERENT calls — distinct
            // (qAsm, refAsm) pairs in the multi-sample VCF. The earlier
            // type-only key collapsed all per-ref pairwise emissions into a
            // single call per (type, pos), undoing the cross-ref dominance
            // gate fix above: 11 of every 12 cross-ref signals were lost in
            // dedup, even though we explicitly kept them through the gate.
            auto& kept = keptByType[sc.call->type + "\x1f" + sc.call->refAsm];

            const int cs = sc.call->pos, ce = sc.call->end;
            const int svl = std::abs(sc.call->svlen);
            bool overlaps = false;
            for (const auto& kp : kept) {
                const int ks = kp.pos, ke = kp.end;
                // Same-type / same-ref proximity window: two same-type calls
                // on one query contig against the SAME ref whose positions
                // fall within ~half the event length AND whose lengths agree
                // (within 2×) are the SAME SV reported by two alignment
                // modes — e.g. mem_chain_ds13 and mem_chain_cached emitting
                // the identical 500 bp DEL a few dozen bp apart. Collapsing
                // them is the dominant same-type FP source. Genuinely
                // distinct co-located SVs (different size or >window apart)
                // still coexist via the bare mergeSlop.
                const int sameSvl = std::abs(kp.svlen);
                // Only collapse two same-type calls that are genuinely the
                // SAME event reported twice: near-identical size (within
                // 1.5x + 50 bp) AND physically close. The previous rule
                // (2x + 100 bp size tolerance, window up to 1000 bp) fused
                // distinct co-located SVs in TE-rich fungal compartments —
                // exactly the high-density regions that carry most of the
                // real SV burden — and was a major drag on call volume.
                const bool lenAgree =
                    2 * std::max(svl, sameSvl) <= 3 * std::min(svl, sameSvl) + 100;
                const int win = lenAgree
                    ? std::max(mergeSlop, std::min(std::max(svl, sameSvl) / 6, 100))
                    : mergeSlop;
                if (cs <= ke + win && ce + win >= ks) {
                    overlaps = true;
                    break;
                }
            }
            if (overlaps) continue;

            VariantCallBridge chosen = *sc.call;
            chosen.queryMode = modeLabel;
            chosen.fusedPosteriorAlt = sc.fused.posteriorAlt;
            chosen.fusedLogOddsAlt = sc.fused.logOddsAlt;
            chosen.fusedEffectiveDepth = sc.fused.effectiveDepth;
            chosen.fusedLayersUsed = static_cast<int>(sc.fused.layersUsed);
            out.push_back(std::move(chosen));
            kept.push_back({cs, ce, sc.call->svlen});
            if (++emittedThisContig >= contigCap) break;
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
        // Tightened position bucket from 250 → 25 bp. The wider bucket fused
        // genuine adjacent SVs in TE-rich loci where the literature finds the
        // densest signal (Fusarium two-speed compartments, Zymoseptoria
        // accessory chromosomes). svlen bucket is 50 bp — together they
        // require both close position and similar size to collapse.
        if (c.type == "OFF_REF") {
            const int posBucket = std::max(0, c.pos) / 25;
            return c.type + ":" + c.annotation + ":" + c.qContig + ":" +
                   std::to_string(posBucket);
        }
        if (c.type == "INS" || c.type == "DEL" || c.type == "DUP" ||
            c.type == "INV" || is_translocation_type(c.type)) {
            const int bucket = std::max(1, std::abs(c.svlen) / 50);
            const int posBucket = std::max(0, c.refPos > 0 ? c.refPos : c.pos) / 25;
            return c.type + ":" + c.refAsm + ":" + c.refContig + ":" +
                   std::to_string(posBucket) + ":" + std::to_string(bucket);
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
        // Previously: every cachedSingleRef INV/DUP/TRA was wiped at low
        // coverage, which killed real low-coverage ONT rearrangements (a
        // primary use case for fungal pangenome work). Keep them when chain
        // evidence is reasonable; only drop the most fragile.
        // Algorithmic fix: anchors<3 AND blockScore<8 at low coverage is
        // exactly what real fungal 5–10× ONT data looks like — the chain
        // is short by sequencing depth, not by artifact. Use a much
        // gentler gate that only drops the weakest combination, and
        // recognise type-specific evidence (mated TRA, mem_chain INV
        // with cross-strand evidence) as adequate even at minimal anchors.
        if (rearrangement && cachedSingleRef) {
            const bool stronglyEvidenced =
                (c.type == "TRA" && !c.mateContig.empty()) ||
                c.blockScore >= 6.0 ||
                c.anchors >= 2 ||
                parse_pseudo_contig_read_support(c.qContig) >= 2;
            if (!stronglyEvidenced) continue;
        }
        if (report.mode == query_input::QueryMode::LONG_READS &&
            anchoredIndel && cachedSingleRef &&
            parse_pseudo_contig_read_support(c.qContig) < 1) {
            // Single-read INS/DEL at 5–10× ONT *can* be real, especially in
            // TE-rich loci where coverage collapses. Only drop if pseudo
            // contig carries explicit read support of 0.
            continue;
        }
        out.push_back(std::move(c));
    }
    return out;
}

// ----------------------------------------------------------------------------
// Low-coverage cosine-similarity genotyper.
//
// Even at 1–2× depth, the *pattern* of evidence on the candidate haplotype
// graph is informative — a HET site has reads supporting both REF and ALT
// nodes in roughly equal proportion, a HOM_ALT site has reads almost
// exclusively on ALT nodes, and a REF site has reads almost exclusively on
// REF nodes. Cosine similarity is invariant to the L2 norm of the vector,
// which is exactly the depth, so it gives stable genotype calls where a
// raw allele-balance threshold would flap with one extra read.
//
// We do not currently expose per-graph-node coverage to this layer, so
// the per-call evidence vector is reconstructed from quantities the
// federated caller already produces:
//   • fusedPosteriorAlt          → expected ALT mass
//   • 1 − fusedPosteriorAlt      → expected REF mass
//   • fusedLayersUsed            → number of supporting observations
//                                   (treated as how many "graph nodes"
//                                    contributed evidence, the proxy for
//                                    the read-vs-haplotype overlap
//                                    cardinality)
//
// Template haplotype patterns at the same dimensionality (REF-mass,
// ALT-mass, support-cardinality):
//   REF      = (1, 0, k)       reads land on REF nodes only
//   HET      = (1, 1, k)       reads land on both
//   HOM_ALT  = (0, 1, k)       reads land on ALT nodes only
// where k normalises the cardinality channel; the cosine then collapses
// to the angle in the (REF, ALT) plane while the support channel stays
// constant — exactly the "depth-insensitive" behaviour the algorithm
// note from 2026-05-12 calls for. GQ becomes the cosine margin in dB.
static void assign_low_coverage_genotype_cosine(
        std::vector<VariantCallBridge>& calls,
        const query_input::CoverageReport& report) {
    if (report.coverageTier != query_input::CoverageTier::LOW || calls.empty())
        return;
    auto cos_sim = [](double a0, double a1, double a2,
                      double b0, double b1, double b2) -> double {
        const double dot = a0 * b0 + a1 * b1 + a2 * b2;
        const double na  = std::sqrt(a0 * a0 + a1 * a1 + a2 * a2);
        const double nb  = std::sqrt(b0 * b0 + b1 * b1 + b2 * b2);
        if (na <= 0.0 || nb <= 0.0) return 0.0;
        return dot / (na * nb);
    };
    for (auto& c : calls) {
        const double pAlt = std::min(1.0, std::max(0.0, c.fusedPosteriorAlt));
        const double pRef = 1.0 - pAlt;
        // Anchor the support channel at 1.0 so a single-observation call
        // still contributes a non-zero magnitude on that axis. Without
        // this floor, the (1,0,0) and (0,1,0) templates fully dominate
        // and HET is unreachable for layersUsed=0 calls.
        const double k    = 1.0 + std::min(8.0, static_cast<double>(c.fusedLayersUsed));
        const double cosRef = cos_sim(pRef, pAlt, k, 1.0, 0.0, k);
        const double cosHet = cos_sim(pRef, pAlt, k, 1.0, 1.0, k);
        const double cosHom = cos_sim(pRef, pAlt, k, 0.0, 1.0, k);
        const char* gt = "0/1";
        double best = cosHet;
        double runnerUp = std::max(cosRef, cosHom);
        if (cosHom >= cosHet && cosHom >= cosRef) {
            gt = "1/1";
            best = cosHom;
            runnerUp = std::max(cosRef, cosHet);
        } else if (cosRef >= cosHet && cosRef >= cosHom) {
            gt = "0/0";
            best = cosRef;
            runnerUp = std::max(cosHet, cosHom);
        }
        c.genotype = gt;
        // GQ as cosine-margin in dB-style units: clamp to [1, 60] so it
        // composes with the existing GQ histogram. A flat tie (margin≈0)
        // becomes GQ≈1; a clean separation (margin≈0.3) becomes ~30.
        const double margin = std::max(0.0, best - runnerUp);
        c.gq = std::min(60.0, std::max(1.0, 1.0 + 100.0 * margin));
    }
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
    // Multi-ref SA is a pangenome-rescue path, not the main pairwise caller.
    // It still needs enough text to hold ordinary fungal chromosomes/scaffolds;
    // a 1 MB cap silently skipped nearly every real contig and disabled
    // cross-contig TRA evidence.  Keep it bounded for the old compact_yeast
    // hang, but allow a modest chromosome-scale bundle.
    const size_t userMultiRefCap = static_cast<size_t>(std::max(0, o.saMaxContigMB)) * 2;
    fo.saMaxTextMB = userMultiRefCap == 0 ? 64 : std::min<size_t>(userMultiRefCap, 64);
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
        << "##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"Call support: assembly anchor count, long-read cluster size, or short-read k-mer/unitig support\">\n"
        << "##INFO=<ID=FUSED_POST,Number=1,Type=Float,Description=\"Posterior alt probability after probabilistic evidence fusion\">\n"
        << "##INFO=<ID=FUSED_LOGODDS,Number=1,Type=Float,Description=\"Log-odds for the alternative allele after evidence fusion\">\n"
        << "##INFO=<ID=FUSED_DEPTH,Number=1,Type=Float,Description=\"Effective depth used by evidence fusion\">\n"
        << "##INFO=<ID=FUSED_LAYERS,Number=1,Type=Integer,Description=\"Number of evidence observations fused for the chosen call\">\n"
        << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        << "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"GQ\">\n"
        << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
}

static int effective_call_support(const VariantCallBridge& v) {
    if (v.readSupport >= 0) return v.readSupport;
    if (v.queryMode == "assembly") {
        if (v.anchors > 0) return v.anchors;
        if (v.fusedLayersUsed > 0) return v.fusedLayersUsed;
        return 0;
    }
    return v.readSupport;
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
    const int support = effective_call_support(v);
    if (support >= 0) out << ";SUPPORT=" << support;
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
        << '\t' << effective_call_support(v)
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

// CheckpointWriters: optional sink for per-phase partial outputs. When provided,
// process_query writes the raw hierarchical-phase calls to hier_vcf / hier_tsv
// before the flat-MEM-chain fallback starts, so a SIGKILL during the fallback
// does not discard work the hierarchical phase already completed.
struct CheckpointWriters {
    std::ofstream*    hier_tsv = nullptr;
    std::ofstream*    hier_vcf = nullptr;
    std::mutex*       mu       = nullptr;
    std::atomic<int>* sv_id    = nullptr;
};

static QueryResult
process_query(const std::string& qAsmPath,
              const Options& o,
              const tol::FederatedOptions& fo,
              const SimpleRefIndex* refIdx = nullptr,
              const CheckpointWriters* ckpt = nullptr) {
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
        const size_t saMaxBytes    = static_cast<size_t>(std::max(0, o.saMaxContigMB)) * 1024 * 1024;
        const size_t cacheMaxBytes = static_cast<size_t>(std::max(0, o.singleRefCacheMB)) * 1024 * 1024;
        singleRefMemCache.emplace(saMaxBytes, cacheMaxBytes);
    }
    auto add_candidates = [&](std::vector<VariantCallBridge>&& calls) {
        for (auto& c : calls) {
            c.queryMode = modeLabel;
            candidates[c.qContig].push_back(std::move(c));
        }
    };

    size_t hierarchical_call_count = 0;
    if (qo.useTolHierarchical && tol::TolGlobal::instance().is_initialized()) {
        if (!qo.quiet)
            std::cerr << "[query-progress] " << qr.qAsm
                      << ": hierarchical routing/calling start (contigs="
                      << qr.contigs.size() << ")\n";
        const bool per_contig_checkpoint =
            ckpt && ckpt->hier_tsv && ckpt->hier_vcf && ckpt->mu && ckpt->sv_id;
        size_t per_contig_checkpoint_calls = 0;
        std::unique_ptr<tol::ScopedPerContigFlushHook> scoped_flush_hook;
        if (per_contig_checkpoint) {
            scoped_flush_hook = std::make_unique<tol::ScopedPerContigFlushHook>(
                [ckpt, &modeLabel, &per_contig_checkpoint_calls](
                    const std::string&,
                    const std::string&,
                    const std::vector<VariantCallBridge>& contigCalls) {
                    if (contigCalls.empty()) return;
                    std::lock_guard<std::mutex> lk(*ckpt->mu);
                    for (const auto& c : contigCalls) {
                        VariantCallBridge cc = c;
                        if (cc.queryMode.empty()) cc.queryMode = modeLabel;
                        write_tsv_record(*ckpt->hier_tsv, cc);
                        write_vcf_record(*ckpt->hier_vcf, cc, ckpt->sv_id->fetch_add(1));
                        ++per_contig_checkpoint_calls;
                    }
                    ckpt->hier_tsv->flush();
                    ckpt->hier_vcf->flush();
                });
        }
        auto hcalls = qo.tolMultiRank
            ? tol::hierarchical_call_assembly_multirank(
                  qr.qAsm, qr.contigs, eff_fo,
                  static_cast<size_t>(std::max(1, qo.routingTopN)))
            : tol::hierarchical_call_assembly(qr.qAsm, qr.contigs, eff_fo);
        scoped_flush_hook.reset();
        hierarchical_call_count = hcalls.size();
        if (!qo.quiet)
            std::cerr << "[query-progress] " << qr.qAsm
                      << ": hierarchical routing/calling done (calls="
                      << hierarchical_call_count << ")\n";
        // Bug 1 fix: write raw hierarchical calls to the checkpoint outputs
        // BEFORE the flat-MEM-chain fallback runs. The flat fallback can
        // iterate thousands of contigs and exceed the wrapper's
        // mycosv_tool_timeout; without this write, SIGKILL discards the
        // hierarchical work entirely. The checkpoint files are append-only
        // and represent a "guaranteed at least this much" view.
        if (per_contig_checkpoint) {
            if (!qo.quiet)
                std::cerr << "[query-progress] " << qr.qAsm
                          << ": hierarchical per-contig checkpoint flushed ("
                          << per_contig_checkpoint_calls << " calls)\n";
        } else if (ckpt && ckpt->hier_tsv && ckpt->hier_vcf && ckpt->mu && ckpt->sv_id) {
            std::lock_guard<std::mutex> lk(*ckpt->mu);
            for (auto& c : hcalls) {
                VariantCallBridge cc = c;
                if (cc.queryMode.empty()) cc.queryMode = modeLabel;
                write_tsv_record(*ckpt->hier_tsv, cc);
                write_vcf_record(*ckpt->hier_vcf, cc, ckpt->sv_id->fetch_add(1));
            }
            ckpt->hier_tsv->flush();
            ckpt->hier_vcf->flush();
            if (!qo.quiet)
                std::cerr << "[query-progress] " << qr.qAsm
                          << ": hierarchical checkpoint flushed ("
                          << hierarchical_call_count << " calls)\n";
        }
        add_candidates(std::move(hcalls));
    }

    // Bug 2 fix: optional gate that skips the flat-MEM-chain fallback when
    // the hierarchical phase already produced enough calls. Set via
    // --skip-flat-if-hier-calls N. The flat fallback iterates every contig
    // in the SimpleRefIndex (thousands on a fungal ref panel) and is the
    // dominant per-query cost; this gate lets the operator keep the
    // fallback as a safety net for queries that didn't route, without
    // paying its cost on queries hierarchical handled.
    const bool skip_flat_due_to_hier =
        (qo.skipFlatIfHierCalls > 0 &&
         hierarchical_call_count >= static_cast<size_t>(qo.skipFlatIfHierCalls));
    if (skip_flat_due_to_hier && refIdx != nullptr && !qo.quiet) {
        std::cerr << "[query-progress] " << qr.qAsm
                  << ": flat MEM-chain fallback SKIPPED (hierarchical calls="
                  << hierarchical_call_count
                  << " >= --skip-flat-if-hier-calls "
                  << qo.skipFlatIfHierCalls << ")\n";
    }

    if (refIdx != nullptr && !skip_flat_due_to_hier) {
        if (prep.report.mode == query_input::QueryMode::ASSEMBLY) {
            if (!qo.quiet)
                std::cerr << "[query-progress] " << qr.qAsm
                          << ": flat MEM-chain fallback start (refs="
                          << refIdx->size() << " contig names)\n";
            try {
                add_candidates(mem_chain_sv_calls(
                    qr.qAsm, qr.contigs, *refIdx, *refCache, *singleRefMemCache,
                    qo, eff_fo, prep.report.mode));
            } catch (const std::bad_alloc&) {
                std::cerr << "[warn] " << qr.qAsm
                          << ": flat MEM-chain fallback skipped after std::bad_alloc; "
                             "continuing with hierarchical/simple candidates\n";
            } catch (const std::exception& e) {
                std::cerr << "[warn] " << qr.qAsm
                          << ": flat MEM-chain fallback skipped: "
                          << e.what() << '\n';
            }
            if (!qo.quiet)
                std::cerr << "[query-progress] " << qr.qAsm
                          << ": flat MEM-chain fallback done\n";
            try {
                add_candidates(simple_length_fallback_calls(qr.qAsm, qr.contigs, *refIdx, qo));
            } catch (const std::exception& e) {
                std::cerr << "[warn] " << qr.qAsm
                          << ": simple length fallback skipped: "
                          << e.what() << '\n';
            }
        } else {
            if (!qo.quiet)
                std::cerr << "[query-progress] " << qr.qAsm
                          << ": reads-mode flat fallback start (refs="
                          << refIdx->size() << " contig names)\n";
            try {
                add_candidates(reads_mode_sv_calls(
                    qr.qAsm, qr.contigs, *refIdx, *refCache, *singleRefMemCache,
                    qo, eff_fo, prep.report.mode));
            } catch (const std::bad_alloc&) {
                std::cerr << "[warn] " << qr.qAsm
                          << ": reads-mode flat fallback skipped after std::bad_alloc; "
                             "continuing with hierarchical/simple candidates\n";
            } catch (const std::exception& e) {
                std::cerr << "[warn] " << qr.qAsm
                          << ": reads-mode flat fallback skipped: "
                          << e.what() << '\n';
            }
            if (!qo.quiet)
                std::cerr << "[query-progress] " << qr.qAsm
                          << ": reads-mode flat fallback done\n";
        }

        if (!qo.quiet)
            std::cerr << "[query-progress] " << qr.qAsm
                      << ": off-reference fallback/select start\n";
        try {
            add_candidates(simple_offref_fallback_calls(
                qr.qAsm, qr.contigs, *refIdx, *refCache, qo, prep.report.mode));
        } catch (const std::bad_alloc&) {
            std::cerr << "[warn] " << qr.qAsm
                      << ": off-reference fallback skipped after std::bad_alloc; "
                         "selecting existing candidates\n";
        } catch (const std::exception& e) {
            std::cerr << "[warn] " << qr.qAsm
                      << ": off-reference fallback skipped: "
                      << e.what() << '\n';
        }

        qr.calls = select_best_call_per_contig(qr.contigs, candidates, *refIdx, *refCache, qo,
                                               prep.report.mode, modeLabel,
                                               &prep.report);
        if (!qo.quiet)
            std::cerr << "[query-progress] " << qr.qAsm
                      << ": off-reference fallback/select done (selected="
                      << qr.calls.size() << ")\n";
        if (prep.report.mode != query_input::QueryMode::ASSEMBLY) {
            qr.calls = deduplicate_read_mode_events(std::move(qr.calls));
            qr.calls = filter_low_coverage_read_artifacts(std::move(qr.calls), prep.report);
            // Re-genotype low-coverage reads-mode calls using cosine
            // similarity over the (REF mass, ALT mass, support cardinality)
            // template patterns; no-op when coverageTier != LOW.
            assign_low_coverage_genotype_cosine(qr.calls, prep.report);
        }
    } else {
        // No flat-ref fallback active (either --no-flat-ref-fallback + tol-hierarchical,
        // or --skip-flat-if-hier-calls fired on this query because the hierarchical
        // phase already produced enough calls).
        // candidates here only carry hierarchical calls; without refIdx we
        // can't rerun candidate_priority_score, but we must still pick the
        // BEST hierarchical candidate per query contig — taking front()
        // discarded higher-scoring calls and was the headline driver of the
        // ~9-call-per-query recall collapse on heavy panels (amf_large /
        // te_rich_pathogen / two_speed_pathogen).
        auto fallback_score = [](const VariantCallBridge& c) {
            // Block score is the primary chain-support metric MycoSV reports;
            // tie-break on |svlen| so larger structural events outrank tiny
            // boundary noise when block scores are equal.
            return c.blockScore +
                   0.001 * static_cast<double>(std::abs(c.svlen)) +
                   (c.type == "OFF_REF" && c.annotation == "NOVEL" ? 1.0 : 0.0);
        };
        const int mergeSlop = std::max(20, qo.minSvLen / 2);
        for (auto& kv : candidates) {
            if (kv.second.empty()) continue;
            const auto seqIt = qr.contigs.find(kv.first);
            const size_t contigLen = (seqIt == qr.contigs.end()) ? 0 : seqIt->second.size();
            const int contigCap = dynamic_calls_per_contig_cap(qo, prep.report.mode, contigLen);
            std::vector<std::pair<const VariantCallBridge*, double>> scored;
            scored.reserve(kv.second.size());
            for (const auto& cand : kv.second) {
                scored.emplace_back(&cand, fallback_score(cand));
            }
            std::sort(scored.begin(), scored.end(),
                      [](const auto& a, const auto& b) { return a.second > b.second; });
            // Top-K non-overlapping per-contig: same idea as the flat-ref path.
            // Single-pick collapsed multi-emit chain output (and the per-chain
            // INS/DEL events surfaced by classify_all) down to one call,
            // throttling recall on chromosome-sized contigs that carry many
            // distinct SVs.
            struct KeptInterval { int pos; int end; int svlen; };
            std::unordered_map<std::string, std::vector<KeptInterval>> keptByType;
            int keptCount = 0;
            for (const auto& sc : scored) {
                if (keptCount >= contigCap) break;
                if (sc.first == nullptr) continue;
                const int cs = sc.first->pos, ce = sc.first->end;
                const int svl = std::abs(sc.first->svlen);
                auto& kept = keptByType[sc.first->type + "\x1f" + sc.first->refAsm];
                bool overlaps = false;
                for (const auto& kp : kept) {
                    const int ks = kp.pos, ke = kp.end;
                    const int sameSvl = std::abs(kp.svlen);
                    const bool lenAgree =
                        2 * std::max(svl, sameSvl) <= 3 * std::min(svl, sameSvl) + 100;
                    const int win = lenAgree
                        ? std::max(mergeSlop, std::min(std::max(svl, sameSvl) / 4, 250))
                        : mergeSlop;
                    if (cs <= ke + win && ce + win >= ks) {
                        overlaps = true;
                        break;
                    }
                }
                if (overlaps) continue;
                VariantCallBridge chosen = *sc.first;
                chosen.queryMode = modeLabel;
                qr.calls.push_back(std::move(chosen));
                kept.push_back({cs, ce, sc.first->svlen});
                ++keptCount;
            }
        }
        if (prep.report.mode != query_input::QueryMode::ASSEMBLY) {
            qr.calls = deduplicate_read_mode_events(std::move(qr.calls));
            qr.calls = filter_low_coverage_read_artifacts(std::move(qr.calls), prep.report);
            // Re-genotype low-coverage reads-mode calls using cosine
            // similarity over the (REF mass, ALT mass, support cardinality)
            // template patterns; no-op when coverageTier != LOW.
            assign_low_coverage_genotype_cosine(qr.calls, prep.report);
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

    // emit_ref_anchor: write a placeholder S-line (sequence "*") for the
    // reference contig that a variant is anchored to, plus a path containing
    // just that anchor.  The S-line carries CL:Z:<refAsm> so downstream tools
    // can group anchors by reference assembly.  Returns the segment id which
    // is used by the L-edge that links the variant to its reference position.
    //
    // Using a per-(refAsm,refContig) segment keeps the augmented GFA compact
    // even when many variants share the same reference contig — the segment
    // is emitted once and every variant that maps to it adds an L-line.
    auto emit_ref_anchor = [&](const std::string& refAsm,
                                const std::string& refContig) -> std::string {
        const std::string id = sanitize(refAsm + ":" + refContig + ":REF");
        if (seen.insert("S:" + id).second) {
            out << "S\t" << id << "\t*"
                << "\tAN:Z:REFERENCE"
                << "\tVT:Z:REF"
                << "\tCL:Z:" << (refAsm.empty() ? std::string(".") : refAsm)
                << "\tEC:Z:NONE"
                << "\n";
            out << "P\t" << id << "\t" << id << "+\t*\n";
        }
        return id;
    };
    auto emit_anchor_link = [&](const std::string& refSeg,
                                 const std::string& varSeg,
                                 const std::string& type,
                                 const std::string& clade) {
        if (refSeg.empty() || varSeg.empty()) return;
        const std::string fwd = "L:" + refSeg + "->" + varSeg;
        if (seen.insert(fwd).second) {
            out << "L\t" << refSeg << "\t+\t" << varSeg << "\t+\t0M"
                << "\tVT:Z:" << type
                << "\tCL:Z:" << (clade.empty() ? std::string(".") : clade)
                << "\tAN:Z:LEFT_FLANK\n";
        }
        const std::string rev = "L:" + varSeg + "->" + refSeg;
        if (seen.insert(rev).second) {
            out << "L\t" << varSeg << "\t+\t" << refSeg << "\t+\t0M"
                << "\tVT:Z:" << type
                << "\tCL:Z:" << (clade.empty() ? std::string(".") : clade)
                << "\tAN:Z:RIGHT_FLANK\n";
        }
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
            // Anchor each TRA breakend to its reference contig so the GFA is
            // an augmentation of the reference backbone, not a free-floating
            // catalog of variant segments.  A breakend without a refContig
            // (e.g. OFF_REF) is skipped — there is no anchor to attach to.
            if (!primaryOffRef && !c.refAsm.empty() && !c.refContig.empty()) {
                const std::string refSeg = emit_ref_anchor(c.refAsm, c.refContig);
                emit_anchor_link(refSeg, left, c.type, c.refAsm);
            }
            if (!c.mateOffReference && !c.mateRefAsm.empty() && !c.mateContig.empty()
                && c.mateContig != ".") {
                const std::string mateRefSeg = emit_ref_anchor(c.mateRefAsm, c.mateContig);
                emit_anchor_link(mateRefSeg, right, "TRA_MATE", c.mateRefAsm);
            }
            continue;
        }

        const int svlen = normalized_svlen_for_output(c);
        const int end = normalized_end_for_output(c);
        const bool placeholderOnly = (c.type == "DEL");
        const std::string varSeg = emit_segment(
            qr.qAsm + ":" + c.qContig + ":" + std::to_string(std::max(1, c.pos)) + "-" + std::to_string(end),
            c.qContig,
            c.pos,
            end,
            svlen,
            primaryOffRef ? "OFF_REFERENCE" : c.annotation,
            c.type,
            c.refAsm,
            c.elementClass,
            placeholderOnly);
        // Augmentation: every reference-anchored variant gets two L-edges
        // (ref→var, var→ref) to its parent reference contig.  OFF_REF variants
        // have no ref anchor, so they remain disconnected — that is the correct
        // graph topology for novel sequence in a pangenome augmentation.
        if (!primaryOffRef && !c.refAsm.empty() && !c.refContig.empty()) {
            const std::string refSeg = emit_ref_anchor(c.refAsm, c.refContig);
            emit_anchor_link(refSeg, varSeg, c.type, c.refAsm);
        }
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
    install_signal_handlers();
    Options o;
    try {
        o = parse_args(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "[error] " << e.what() << '\n';
        return 1;
    }

    // ---- Diagnose mode -----------------------------------------------
    // Login-node-fast health checks; exits before any SV calling. Use this
    // when the multisample VCF is empty/all-OFF_REF or when changing
    // routing-index inputs — catches CRLF-mangled paths, missing FASTAs,
    // and unloaded registry sequences in seconds.
    if (!o.diagnoseMode.empty()) {
        if (o.diagnoseMode == "registry") {
            if (o.tolRegistryDir.empty()) {
                std::cerr << "[diagnose] --tol-registry-dir is required for --diagnose registry\n";
                return 2;
            }
            const fs::path manifestPath = fs::path(o.tolRegistryDir) / "clade_manifest.tsv";
            std::cout << "[diagnose] registry: " << o.tolRegistryDir << "\n";
            std::cout << "[diagnose] manifest: " << manifestPath.string()
                      << "  exists=" << (fs::exists(manifestPath) ? "yes" : "NO") << "\n";

            // 1) Raw manifest CR check — counts ANY \r in the file, not just
            // trailing CR. The original bug ([Bug 5]) had \r mid-line, embedded
            // inside column 8 (fasta_paths), so a trailing-only check missed it.
            size_t cr_anywhere = 0, total_lines = 0;
            if (fs::exists(manifestPath)) {
                std::ifstream in(manifestPath);
                std::string line;
                while (std::getline(in, line)) {
                    ++total_lines;
                    for (char c : line) if (c == '\r') ++cr_anywhere;
                }
            }
            std::cout << "[diagnose] manifest: " << total_lines
                      << " lines, " << cr_anywhere
                      << " stray CR character(s) in body";
            if (cr_anywhere > 0) {
                std::cout << "  (CRLF/embedded-CR detected — readers MUST strip \\r; "
                             "see split_tab/split_csv in fungi_tol_bridge.hpp; "
                             "ignoring this corrupts fs::exists on the final path "
                             "of every multi-FASTA cell)";
            }
            std::cout << "\n";

            // 2) Manifest-only sanity (cheap, no FASTA reads). Always runs.
            tol::ManifestRegistry quickReg(o.tolRegistryDir);
            try { quickReg.load_from_disk(); }
            catch (const std::exception& e) {
                std::cerr << "[diagnose] manifest load threw: " << e.what() << "\n";
                return 1;
            }
            const auto& descs = quickReg.descriptors();

            size_t descs_with_fastas = 0, total_fastapaths = 0, paths_exist = 0,
                   paths_with_cr = 0, paths_with_ws = 0;
            std::unordered_set<std::string> unique_paths;
            for (const auto& d : descs) {
                if (!d.fastaPaths.empty()) ++descs_with_fastas;
                for (const auto& fp : d.fastaPaths) {
                    ++total_fastapaths;
                    unique_paths.insert(fp);
                    if (fs::exists(fp)) ++paths_exist;
                    if (fp.find('\r') != std::string::npos) ++paths_with_cr;
                    if (!fp.empty() && (std::isspace(static_cast<unsigned char>(fp.front())) ||
                                         std::isspace(static_cast<unsigned char>(fp.back()))))
                        ++paths_with_ws;
                }
            }

            std::cout << "[diagnose] descriptors=" << descs.size()
                      << " with_fastas=" << descs_with_fastas
                      << "\n";
            std::cout << "[diagnose] fasta_paths total=" << total_fastapaths
                      << " unique=" << unique_paths.size()
                      << " exist_on_disk=" << paths_exist
                      << " with_embedded_CR=" << paths_with_cr
                      << " with_leading_or_trailing_WS=" << paths_with_ws
                      << "\n";

            // Sample 3 paths for eyeballing.
            std::cout << "[diagnose] sample registry paths:\n";
            size_t shown = 0;
            for (const auto& d : descs) {
                if (shown >= 3) break;
                if (d.fastaPaths.empty()) continue;
                const std::string& p = d.fastaPaths.front();
                std::cout << "  clade=" << d.cladeName
                          << " rank=" << d.cladeRank
                          << " fasta=" << p
                          << " exists=" << (fs::exists(p) ? "yes" : "NO")
                          << "\n";
                ++shown;
            }

            // 3) Optional full-load check — only if --ref-list provides a
            // constraint. Without it, TolGlobal::init would try to load all
            // 9 000+ FASTAs (~300 GB) which OOMs on a login node. With a
            // constraint, we load only those refs and report allRefs_ size.
            size_t allRefsAfter = 0, byContigAfter = 0;
            bool ranInit = false;
            if (!o.refList.empty()) {
                std::unordered_set<std::string> allowedTolFastas;
                for (const auto& p : read_list(o.refList))
                    if (!p.empty()) allowedTolFastas.insert(p);
                if (!allowedTolFastas.empty()) {
                    std::cout << "[diagnose] running TolGlobal::init with --ref-list filter ("
                              << allowedTolFastas.size() << " allowed FASTAs)...\n";
                    try {
                        tol::TolGlobal::instance().init(
                            o.tolIndexDir.empty() ? o.tolRegistryDir : o.tolIndexDir,
                            o.tolRegistryDir,
                            o.tolCacheGB << 30,
                            o.tolCacheEntries,
                            &allowedTolFastas);
                        ranInit = true;
                        allRefsAfter   = tol::TolGlobal::instance().all_refs().size();
                        byContigAfter  = tol::TolGlobal::instance().refs_by_contig().size();
                        std::cout << "[diagnose] TolGlobal.allRefs_=" << allRefsAfter
                                  << " refs_by_contig=" << byContigAfter << "\n";
                    } catch (const std::exception& e) {
                        std::cerr << "[diagnose] TolGlobal::init threw: " << e.what() << "\n";
                        return 1;
                    }
                }
            } else {
                std::cout << "[diagnose] skipping TolGlobal::init (no --ref-list "
                             "supplied; pass --ref-list to verify allRefs_ load)\n";
            }

            // Verdict — non-zero exit on any clear-cut breakage so callers
            // (smoke scripts, pre-flight checks in submit.sh) can branch on rc.
            int rc = 0;
            std::cout << "[diagnose] verdict: ";
            if (descs.empty()) {
                std::cout << "BROKEN — no descriptors parsed from manifest";
                rc = 1;
            } else if (paths_exist == 0 && total_fastapaths > 0) {
                std::cout << "BROKEN — every manifest FASTA missing on disk "
                             "(CR-mangled?=" << paths_with_cr << ", total="
                          << total_fastapaths << ")";
                rc = 1;
            } else if (ranInit && allRefsAfter == 0) {
                std::cout << "BROKEN — TolGlobal::init produced allRefs_=0 even with --ref-list";
                rc = 1;
            } else if (cr_anywhere > 0) {
                std::cout << "DEGRADED — manifest contains " << cr_anywhere
                          << " stray CR char(s); readers must strip "
                             "(split_tab/split_csv in fungi_tol_bridge.hpp do this since [Bug 5] fix)";
                rc = 0;
            } else if (paths_with_cr > 0) {
                std::cout << "DEGRADED — " << paths_with_cr
                          << " path(s) contain \\r; readers must strip (split_tab/split_csv do)";
                rc = 0;
            } else {
                std::cout << "OK";
            }
            std::cout << "\n";
            return rc;
        }
        std::cerr << "[diagnose] unknown mode '" << o.diagnoseMode
                  << "'.  Supported: registry\n";
        return 2;
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
            std::unordered_set<std::string> allowedTolFastas;
            if (!o.refList.empty()) {
                for (const auto& p : read_list(o.refList)) {
                    if (!p.empty()) allowedTolFastas.insert(p);
                }
            }
            const std::unordered_set<std::string>* allowedTolFastasPtr =
                allowedTolFastas.empty() ? nullptr : &allowedTolFastas;
            tol::TolGlobal::instance().init(
                o.tolIndexDir,
                o.tolRegistryDir,
                o.tolCacheGB << 30,
                o.tolCacheEntries,
                allowedTolFastasPtr);
            if (o.tolMultiRank || o.tolAncestralAlign) {
                tol::MultiRankIndex::instance().init(
                    o.tolIndexDir,
                    o.tolRegistryDir,
                    o.tolCacheGB << 30,
                    o.tolCacheEntries,
                    allowedTolFastasPtr);
            }
        } catch (const std::exception& e) {
            std::cerr << "[error] TOL init failed: " << e.what() << '\n';
            return 1;
        }

        // Pre-flight: refuse to run SV calling when no reference sequences
        // loaded. Without this guard the hierarchical caller would silently
        // fall through to the whole-contig OFF_REF safety net (one
        // SVTYPE=OFF_REF per query contig, NOVEL tier, no real SVs) and
        // produce a misleading "successful" run. Common causes: CRLF-encoded
        // manifest paths, --ref-list filter that excludes every registry
        // FASTA, missing FASTAs on disk after a partial download.
        const auto& tolAllRefs = tol::TolGlobal::instance().all_refs();
        if (tolAllRefs.empty()) {
            std::cerr << "[error] --tol-hierarchical: zero reference sequences loaded "
                         "into TolGlobal after init. The hierarchical caller would emit "
                         "only OFF_REF whole-contig calls.\n"
                         "        Run `" << argv[0] << " --diagnose registry "
                         "--tol-registry-dir " << o.tolRegistryDir;
            if (!o.tolIndexDir.empty()) std::cerr << " --tol-index-dir " << o.tolIndexDir;
            if (!o.refList.empty())     std::cerr << " --ref-list " << o.refList;
            std::cerr << "` for a per-cause breakdown.\n";
            return 3;
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
    // Bug 1 fix: hierarchical-only checkpoint outputs written before the
    // flat-MEM-chain fallback starts. These are durable across SIGTERM/SIGKILL
    // so the operator always has the hierarchical-phase callset available.
    std::string hier_tsv_path = o.outPrefix + ".hierarchical.hits.tsv";
    std::string hier_vcf_path = o.outPrefix + ".hierarchical.vcf";
    if (auto out_parent = fs::path(o.outPrefix).parent_path(); !out_parent.empty()) {
        fs::create_directories(out_parent);
    }
    std::ofstream tsv_out(tsv_path);
    std::ofstream vcf_out(vcf_path);
    std::ofstream hier_tsv_out(hier_tsv_path);
    std::ofstream hier_vcf_out(hier_vcf_path);
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
    if (!hier_tsv_out) { std::cerr << "[error] cannot open " << hier_tsv_path << '\n'; return 1; }
    if (!hier_vcf_out) { std::cerr << "[error] cannot open " << hier_vcf_path << '\n'; return 1; }

    tsv_out << "query_asm\tquery_contig\ttype\tref_asm\tref_contig"
               "\tref_pos\tref_end\tpos\tend\tsvlen\tblock_score\tanchors"
               "\tgenotype\tgq\tannotation\talignment_mode\tquery_mode"
               "\tfused_posterior_alt\tfused_logodds_alt\tfused_effective_depth\tfused_layers"
               "\tread_support\n";
    hier_tsv_out << "query_asm\tquery_contig\ttype\tref_asm\tref_contig"
                    "\tref_pos\tref_end\tpos\tend\tsvlen\tblock_score\tanchors"
                    "\tgenotype\tgq\tannotation\talignment_mode\tquery_mode"
                    "\tfused_posterior_alt\tfused_logodds_alt\tfused_effective_depth\tfused_layers"
                    "\tread_support\n";
    write_vcf_header(vcf_out, "fungi_graphsv_tol_v3");
    write_vcf_header(hier_vcf_out, "fungi_graphsv_tol_v3_hierarchical_checkpoint");
    tsv_out.flush();
    vcf_out.flush();
    hier_tsv_out.flush();
    hier_vcf_out.flush();
    if (o.tolAncestralAlign) tol::write_ancestral_tsv_header(anc_out);

    // Bug 3 fix: register open streams with the signal handler so SIGTERM
    // can best-effort flush partial output before the wrapper's SIGKILL.
    g_signal_tsv_out      = &tsv_out;
    g_signal_vcf_out      = &vcf_out;
    g_signal_hier_tsv_out = &hier_tsv_out;
    g_signal_hier_vcf_out = &hier_vcf_out;
    if (!o.noGfa)            g_signal_gfa_out = &gfa_out;
    if (o.tolAncestralAlign) g_signal_anc_out = &anc_out;

    const tol::FederatedOptions fo = make_fed_opts(o);
    std::optional<SimpleRefIndex> simpleRefIdx;
    const bool skipFlatRefFallback = o.noFlatRefFallback && o.useTolHierarchical;
    if (skipFlatRefFallback) {
        if (!o.quiet)
            std::cerr << "[info] --no-flat-ref-fallback: skipping flat reference "
                         "index load; using hierarchical calls only\n";
    } else {
        if (o.noFlatRefFallback && !o.quiet)
            std::cerr << "[warn] --no-flat-ref-fallback ignored because "
                         "--tol-hierarchical is not enabled\n";
        simpleRefIdx.emplace(load_simple_ref_index(o));
    }
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
    // Bug 1 fix: dedicated mutex + atomic id for the hierarchical checkpoint
    // stream. Independent of the main vcf_mutex so the checkpoint flush
    // inside process_query never contends with the per-query final flush.
    std::mutex           hier_mutex;
    std::atomic<int>     hier_sv_id{1};
    std::atomic<int>     sv_id{1};
    std::atomic<size_t>  total_calls{0};
    std::atomic<size_t>  queries_done{0};
    std::atomic<size_t>  query_failures{0};
    const size_t         n_queries = queries.size();
    CheckpointWriters    ckpt{&hier_tsv_out, &hier_vcf_out, &hier_mutex, &hier_sv_id};

    std::unordered_set<std::string> gfa_seen;
    auto process_one = [&](const std::string& qpath) {
        // Bug 3 fix: respect shutdown requests so a SIGTERM mid-run stops
        // launching new queries instead of trying to start fresh ones that
        // will just be killed by the impending SIGKILL.
        if (g_shutdown_requested.load(std::memory_order_relaxed)) {
            if (!o.quiet)
                std::cerr << "[mycosv] shutdown requested, skipping remaining query "
                          << qpath << "\n";
            ++queries_done;
            return;
        }
        QueryResult qr;
        try {
            qr = process_query(qpath, o, fo,
                               simpleRefIdx ? &(*simpleRefIdx) : nullptr,
                               &ckpt);
        } catch (const std::bad_alloc&) {
            ++query_failures;
            if (!o.quiet)
                std::cerr << "[warn] skipping " << qpath
                          << ": std::bad_alloc while processing query\n";
            size_t done = ++queries_done;
            if (!o.quiet)
                std::cerr << "[progress] " << done << '/' << n_queries
                          << " queries, " << total_calls.load() << " calls\n";
            return;
        } catch (const std::exception& e) {
            ++query_failures;
            if (!o.quiet)
                std::cerr << "[warn] skipping " << qpath
                          << ": " << e.what() << '\n';
            size_t done = ++queries_done;
            if (!o.quiet)
                std::cerr << "[progress] " << done << '/' << n_queries
                          << " queries, " << total_calls.load() << " calls\n";
            return;
        }

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
        {
            std::lock_guard<std::mutex> lk(tsv_mutex);
            tsv_out.flush();
        }
        {
            std::lock_guard<std::mutex> lk(vcf_mutex);
            vcf_out.flush();
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
        if (!o.quiet)
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
                  << "  VCF: " << vcf_path << '\n'
                  << "  Hierarchical checkpoint TSV: " << hier_tsv_path << '\n'
                  << "  Hierarchical checkpoint VCF: " << hier_vcf_path << '\n';

    // Bug 3 fix: clear signal-handler pointers before the streams go out of
    // scope, so a late SIGTERM doesn't reach a destroyed ofstream.
    g_signal_tsv_out      = nullptr;
    g_signal_vcf_out      = nullptr;
    g_signal_hier_tsv_out = nullptr;
    g_signal_hier_vcf_out = nullptr;
    g_signal_gfa_out      = nullptr;
    g_signal_anc_out      = nullptr;

    if (total_calls.load() == 0 && query_failures.load() > 0) {
        if (!o.quiet)
            std::cerr << "[error] all emitted callsets are empty after "
                      << query_failures.load()
                      << " query preprocessing/calling failure(s)\n";
        return 2;
    }

    return 0;
}
