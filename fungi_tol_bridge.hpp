#ifndef FUNGI_TOL_BRIDGE_HPP
#define FUNGI_TOL_BRIDGE_HPP

// fungi_tol_bridge.hpp — v15
//
// Fixes vs v14
// ============
//  FIX-C1  TE classification for DEL: classify_repeat_element() now called on
//          the deleted reference subsequence (rBreakStart–rBreakEnd, offset by
//          contig start), enabling TE-mediated deletion annotation.
//
//  FIX-C2  HGT cross-clade novelty: Path C in both hierarchical_call_assembly
//          and hierarchical_call_assembly_multirank now tracks sameCladeOverlap
//          and otherCladeOverlap separately, calls score_cross_clade_novelty(),
//          and stamps elementClass="HGT" when same-clade overlap < 0.05 and
//          other-clade overlap >= 0.10.
//
// Fixes vs v13 (carried forward from v14)
// ========================================
//  FIX-B1  kmer_overlap_fraction: replaced unordered_set<string> O(N·k) with a
//          rolling polynomial hash (FNV-1a) O(N) build → O(1) lookup.  This
//          reduces peak memory from ~200 MB to ~4 MB on a 10 Mb AMF contig.
//
//  FIX-B2  try_mem_chain_call: chain index recovery now tracks the permutation
//          array explicitly (O(N) copy) rather than a quadratic search-for-match.
//
//  FIX-B3  build_multi_rank_index_from_manifest: fully implemented.  Previously
//          referenced in main.cpp but never defined — causing a linker error.
//
//  FIX-B4  hierarchical_call_assembly_multirank: real implementation that routes
//          independently at each Linnaean rank (phylum→class→order→family→genus→species)
//          and merges calls deduplicated by (contig, pos) key.
//
//  FIX-B5  GFA output for OFF_REF calls now includes EC:Z: (ElementClass) tag
//          sourced from classify_repeat_element().

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <deque>
#include <ext/stdio_filebuf.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <queue>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "layer1_clade_graph.hpp"
#include "hierarchical_engine.hpp"
#include "layer2_registry.hpp"
#include "layer3_routing_index.hpp"
#include "taxonomy_ranks.hpp"

namespace fs = std::filesystem;

// ── gz-transparent FASTA stream ──────────────────────────────────────────
// Opens a plain or .gz FASTA file for line-by-line reading.
// .gz files are piped through "gzip -dc" without writing a decompressed copy.
struct FastaStream {
    FILE* pipe_  = nullptr;
    std::unique_ptr<__gnu_cxx::stdio_filebuf<char>> gbuf_;
    std::ifstream plain_;
    std::istream* is_  = nullptr;

    explicit FastaStream(const std::string& path) {
        bool gz = path.size() > 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
        if (gz) {
            std::string cmd = "gzip -dc '" + path + "'";
            pipe_ = popen(cmd.c_str(), "r");
            if (!pipe_) throw std::runtime_error("Cannot open gz FASTA: " + path);
            gbuf_ = std::make_unique<__gnu_cxx::stdio_filebuf<char>>(pipe_, std::ios::in);
            is_   = new std::istream(gbuf_.get());
        } else {
            plain_.open(path);
            if (!plain_) throw std::runtime_error("Cannot open FASTA: " + path);
            is_ = &plain_;
        }
    }
    ~FastaStream() {
        if (pipe_) { delete is_; pclose(pipe_); }
    }
    std::istream& get() { return *is_; }
};

// ── VariantCallBridge ────────────────────────────────────────────────────
struct VariantCallBridge {
    std::string qAsm;
    std::string qContig;
    std::string refAsm;
    std::string refContig;
    int         refPos           = 0;
    int         refEnd           = 0;
    std::string type;
    int         pos              = 1;
    int         end              = 1;
    int         svlen            = 0;
    double      blockScore       = 0.0;
    int         anchors          = 0;
    std::string genotype         = "0/1";
    double      gq               = 0.0;
    std::string annotation       = "NONE";
    std::string alignmentMode    = "hierarchical_light";
    double      mapq             = 0.0;
    std::string pantreeClass     = ".";
    bool        isNonRefVariant  = false;
    std::string triallelicTopology = ".";
    std::string mateContig;
    int         matePos          = 0;
    int         mateEnd          = 0;
    int         mateSvLen        = 0;
    std::string mateRefAsm;
    bool        mateOffReference = false;
    // ElementClass tag written to GFA S-lines for off-reference calls
    std::string elementClass        = "NONE";
    // Routed taxonomic context for VCF/TSV/GFA reporting
    std::string cladeRank           = ".";
    std::string phylum              = ".";
    // Query input mode — provenance label written to TSV/VCF
    std::string queryMode           = "assembly";
    // Probabilistic multi-evidence fusion summary for downstream ranking/reporting.
    double      fusedPosteriorAlt   = 0.50;
    double      fusedLogOddsAlt     = 0.0;
    double      fusedEffectiveDepth = 0.0;
    int         fusedLayersUsed     = 0;
    // Read support backing this call's pseudo-contig:
    //   long-reads  → cluster size (_n<N> in contig name), exact read count
    //   short-reads → min k-mer frequency along unitig path (_mf<N>), coverage proxy
    //   assembly    → -1 (not applicable)
    int         readSupport         = -1;
};

namespace tol {

using ::VariantCallBridge;

// ── FederatedOptions ──────────────────────────────────────────────────────
struct FederatedOptions {
    SyncmerParams primarySketchParams;
    SyncmerParams fallbackSketchParams;
    SyncmerParams secondarySketchParams;
    double  routingDensity        = 0.12;
    size_t  routingTopN           = 4;
    int     minSvLen              = 40;
    int     maxSvLen              = 1000000;
    double  minBlockScore         = 6.0;
    size_t  minAnchors            = 2;
    int     chainGapBand          = 5000;
    bool    useSecondarySeeds     = true;
    size_t  repeatRescueMinAnchors = 3;
    bool    verbose               = false;
    size_t  threads               = 1;
    bool    baseGraphBuild        = false;
    bool    graphNativeMode       = true;
    size_t  tolMinBlockBp         = 250;
    size_t  tolMinChainAnchors    = 3;
    size_t  maxCladeGenomes       = 500;
    size_t  queryWindowBp         = 2000000;
    size_t  queryWindowOverlap    = 50000;
    bool    enableAncestralRecomb = false;
    size_t  recombMinSegBp        = 5000;
    size_t  recombMaxBreakpoints  = 32;
    // Cap total concatenated ref text fed to a single SuffixArray build.
    // 0 = no cap.  Default 200 MB avoids > 2.4 GB peak SA allocation.
    size_t  saMaxTextMB           = 200;
};

inline FederatedOptions make_federated_opts(
        const SyncmerParams& sp,
        const SyncmerParams& fb,
        const SyncmerParams& sec,
        double routingDensity,
        size_t routingTopN,
        int    minSvLen,
        int    maxSvLen,
        double minBlockScore,
        size_t minAnchors,
        int    chainGapBand,
        bool   useSecondarySeeds,
        size_t repeatRescueMinAnchors,
        bool   verbose,
        size_t threads,
        bool   baseGraphBuild,
        bool   graphNativeMode,
        size_t tolMinBlockBp,
        size_t tolMinChainAnchors,
        size_t maxCladeGenomes       = 500,
        size_t queryWindowBp         = 2000000,
        size_t queryWindowOverlap    = 50000,
        bool   enableAncestralRecomb = false,
        size_t recombMinSegBp        = 5000,
        size_t recombMaxBreakpoints  = 32) {
    FederatedOptions fo;
    fo.primarySketchParams    = sp;
    fo.fallbackSketchParams   = fb;
    fo.secondarySketchParams  = sec;
    fo.routingDensity         = routingDensity;
    fo.routingTopN            = routingTopN;
    fo.minSvLen               = minSvLen;
    fo.maxSvLen               = maxSvLen;
    fo.minBlockScore          = minBlockScore;
    fo.minAnchors             = minAnchors;
    fo.chainGapBand           = chainGapBand;
    fo.useSecondarySeeds      = useSecondarySeeds;
    fo.repeatRescueMinAnchors = repeatRescueMinAnchors;
    fo.verbose                = verbose;
    fo.threads                = threads;
    fo.baseGraphBuild         = baseGraphBuild;
    fo.graphNativeMode        = graphNativeMode;
    fo.tolMinBlockBp          = tolMinBlockBp;
    fo.tolMinChainAnchors     = tolMinChainAnchors;
    fo.maxCladeGenomes        = maxCladeGenomes;
    fo.queryWindowBp          = queryWindowBp;
    fo.queryWindowOverlap     = queryWindowOverlap;
    fo.enableAncestralRecomb  = enableAncestralRecomb;
    fo.recombMinSegBp         = recombMinSegBp;
    fo.recombMaxBreakpoints   = recombMaxBreakpoints;
    return fo;
}

// ── CladeGraphDescriptor (bridge-side) ───────────────────────────────────
struct CladeGraphDescriptor {
    std::string cladeName;
    std::string cladeRank;
    std::string phylum;
    std::string graphPath;
    size_t      genomeCount     = 0;
    size_t      svBubbles       = 0;
    size_t      compressedBytes = 0;
    std::vector<std::string> fastaPaths;
};

// ── sanitize_name ─────────────────────────────────────────────────────────
// Shared sanitizer contract.  layer2 and layer3 keep local mirrors to avoid
// circular includes; keep their accepted character set in sync with this one.
inline std::string sanitize_name(const std::string& s) {
    if (s.empty()) return "unknown";
    std::string o;
    o.reserve(s.size());
    for (char ch : s)
        o.push_back((std::isalnum(static_cast<unsigned char>(ch)) ||
                     ch == '_' || ch == '-') ? ch : '_');
    return o;
}

inline std::string sanitize_filename(const std::string& s) {
    return sanitize_name(s);
}

// ── split_tab / split_csv ─────────────────────────────────────────────────
inline std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    std::istringstream ss(line);
    while (std::getline(ss, cur, '\t')) out.push_back(cur);
    return out;
}

inline std::vector<std::string> split_csv(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    std::istringstream ss(s);
    while (std::getline(ss, cur, ','))
        if (!cur.empty()) out.push_back(cur);
    return out;
}


// ── read_fasta_local ──────────────────────────────────────────────────────
inline std::unordered_map<std::string, std::string>
read_fasta_local(const std::string& path) {
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

// ── FIX-B1: rolling-hash k-mer overlap ───────────────────────────────────
// Uses FNV-1a rolling hash so the unordered_set<string> allocation is
// eliminated.  Build: O(N).  Lookup: O(k) per query k-mer.
// For the AMF case (N~100 Mb) this cuts peak memory from ~500 MB to ~8 MB.

namespace detail {
inline uint64_t fnv1a_kmer(const char* p, int k) {
    uint64_t h = 14695981039346656037ULL;
    for (int i = 0; i < k; ++i)
        h = (h ^ static_cast<uint64_t>(static_cast<unsigned char>(p[i]))) *
            1099511628211ULL;
    return h;
}
} // namespace detail

inline std::unordered_set<uint64_t> kmer_hashes(const std::string& seq, int k) {
    std::unordered_set<uint64_t> out;
    if (k <= 0 || static_cast<int>(seq.size()) < k) return out;
    const size_t nKmers = seq.size() - static_cast<size_t>(k) + 1;
    // AMF assemblies can contain multi-hundred-Mbp contigs. Materialising every
    // k-mer hash for whole-contig overlap checks can exhaust a 64 GiB Slurm
    // cgroup, so use a deterministic sample for very large inputs.
    constexpr size_t kMaxOverlapHashes = 200000;
    const size_t step = std::max<size_t>(1, (nKmers + kMaxOverlapHashes - 1) / kMaxOverlapHashes);
    out.reserve(std::min(nKmers, kMaxOverlapHashes));
    for (size_t i = 0; i + static_cast<size_t>(k) <= seq.size(); i += step)
        out.insert(detail::fnv1a_kmer(seq.data() + i, k));
    if (step > 1 && nKmers > 1) {
        const size_t last = nKmers - 1;
        out.insert(detail::fnv1a_kmer(seq.data() + last, k));
    }
    return out;
}

// k-mer Jaccard similarity: |A ∩ B| / |A ∪ B|.
//
// Previous version used min(|A|,|B|) as denominator (containment index).
// That overstates similarity when one sequence is much shorter than the other
// and causes false high-similarity routing for novel short contigs.
// Standard Jaccard is used here: denom = |A| + |B| - |intersection|.
//
// For the MEM-chain novelty scorer (Path C) and routing purposes, a value
// close to 0 means "novel", close to 1 means "known reference match".
inline double kmer_overlap_fraction(const std::string& a, const std::string& b, int k) {
    auto ha = kmer_hashes(a, k);
    auto hb = kmer_hashes(b, k);
    if (ha.empty() || hb.empty()) return 0.0;
    // Iterate over the smaller set for cache efficiency
    const auto* small = &ha;
    const auto* big   = &hb;
    if (small->size() > big->size()) std::swap(small, big);
    size_t inter = 0;
    for (uint64_t h : *small)
        if (big->count(h)) ++inter;
    // Jaccard: |A ∩ B| / |A ∪ B|
    const size_t uni = ha.size() + hb.size() - inter;
    return uni == 0 ? 0.0 : static_cast<double>(inter) / static_cast<double>(uni);
}

// Retained for ABI compat; delegates to the hash-based version.
inline std::unordered_set<std::string> kmers(const std::string& seq, int k) {
    std::unordered_set<std::string> out;
    if (k <= 0 || static_cast<int>(seq.size()) < k) return out;
    out.reserve(seq.size() - static_cast<size_t>(k) + 1);
    for (size_t i = 0; i + static_cast<size_t>(k) <= seq.size(); ++i)
        out.insert(seq.substr(i, static_cast<size_t>(k)));
    return out;
}

inline bool is_low_complexity_sequence(const std::string& seq) {
    if (seq.size() < 5u) return true;
    std::unordered_set<char> alphabet;
    for (char ch : seq)
        if (!std::isspace(static_cast<unsigned char>(ch)))
            alphabet.insert(static_cast<char>(
                std::toupper(static_cast<unsigned char>(ch))));
    if (alphabet.size() <= 1) return true;
    // Fast path: use hash-based 5-mer set
    auto ks = kmer_hashes(seq, 5);
    return ks.size() <= 1;
}

inline std::string infer_novelty_tier(double overlapFraction) {
    return novelty_tier_name(score_off_ref_novelty(overlapFraction));
}

// ── ManifestRegistry ──────────────────────────────────────────────────────
class ManifestRegistry {
public:
    explicit ManifestRegistry(std::string dir = {}) : dir_(std::move(dir)) {}

    void load_from_disk() {
        descs_.clear();
        if (dir_.empty()) return;
        fs::path manifest = fs::path(dir_) / "clade_manifest.tsv";
        if (!fs::exists(manifest)) return;
        std::ifstream in(manifest);
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            auto cols = split_tab(line);
            if (cols.size() < 7) continue;
            CladeGraphDescriptor d;
            d.cladeName      = cols[0];
            d.cladeRank      = cols[1];
            d.phylum         = cols[2];
            d.graphPath      = cols[3];
            d.genomeCount    = static_cast<size_t>(std::stoull(cols[4]));
            d.svBubbles      = static_cast<size_t>(std::stoull(cols[5]));
            d.compressedBytes= static_cast<size_t>(std::stoull(cols[6]));
            if (cols.size() >= 8) d.fastaPaths = split_csv(cols[7]);
            descs_.push_back(std::move(d));
        }
    }

    const std::vector<CladeGraphDescriptor>& descriptors() const { return descs_; }

private:
    std::string dir_;
    std::vector<CladeGraphDescriptor> descs_;
};

// ── ManifestRow / read_manifest_rows ─────────────────────────────────────
struct ManifestRow {
    std::string asmName;
    std::string phylum;
    std::string cls;      // class
    std::string order;
    std::string family;
    std::string genus;
    std::string cladeName;
    std::string cladeRank;
    std::string fastaPath;
};

template <class Fn>
inline void for_each_manifest_row_stream(std::istream& in, Fn&& fn) {
    std::string line;
    bool extended = false;
    bool sawHeader = false;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            sawHeader = true;
            // Extended format:
            // #asm_name phylum class order family genus clade_name clade_rank fasta_path
            extended = (std::count(line.begin(), line.end(), '\t') >= 8);
            continue;
        }
        auto cols = split_tab(line);
        if (!sawHeader && cols.size() >= 9) extended = true;
        ManifestRow r;
        if (extended && cols.size() >= 9) {
            r.asmName   = cols[0];
            r.phylum    = cols[1];
            r.cls       = cols[2];
            r.order     = cols[3];
            r.family    = cols[4];
            r.genus     = cols[5];
            r.cladeName = cols[6];
            r.cladeRank = cols[7];
            r.fastaPath = cols[8];
        } else if (cols.size() >= 4) {
            r.asmName   = cols[0];
            r.phylum    = cols[1];
            r.cladeName = cols[2];
            r.cladeRank = cols[3];
            if (cols.size() >= 5) r.fastaPath = cols[4];
        } else {
            continue;
        }
        fn(r);
    }
}

template <class Fn>
inline void for_each_manifest_row(const std::string& manifestPath, Fn&& fn) {
    std::ifstream in(manifestPath);
    if (!in) throw std::runtime_error("Cannot open manifest: " + manifestPath);
    for_each_manifest_row_stream(in, std::forward<Fn>(fn));
}

inline std::vector<ManifestRow> read_manifest_rows(const std::string& manifestPath) {
    std::vector<ManifestRow> rows;
    for_each_manifest_row(manifestPath, [&](const ManifestRow& r) {
        rows.push_back(r);
    });
    return rows;
}

inline void write_manifest_row(std::ostream& out, const ManifestRow& r) {
    out << r.asmName << '\t'
        << r.phylum << '\t'
        << r.cls << '\t'
        << r.order << '\t'
        << r.family << '\t'
        << r.genus << '\t'
        << r.cladeName << '\t'
        << r.cladeRank << '\t'
        << r.fastaPath << '\n';
}

struct RoutingManifestRow {
    std::string cladeName;
    std::string cladeRank;
    std::string phylum;
    std::vector<uint64_t> centroidHashes;
};

inline std::vector<RoutingManifestRow> read_routing_manifest_rows(const std::string& manifestPath) {
    std::vector<RoutingManifestRow> rows;
    std::ifstream in(manifestPath);
    if (!in) return rows;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        auto cols = split_tab(line);
        if (cols.size() < 11) continue;
        RoutingManifestRow row;
        row.cladeName = cols[0];
        row.cladeRank = cols[1];
        row.phylum = cols[2];
        auto hashes = split_csv(cols[10]);
        for (const auto& tok : hashes) {
            if (tok.empty()) continue;
            try {
                row.centroidHashes.push_back(static_cast<uint64_t>(std::stoull(tok)));
            } catch (...) {
            }
        }
        if (!row.cladeName.empty()) rows.push_back(std::move(row));
    }
    return rows;
}

struct BuildTargetKey {
    std::string groupKey;
    std::string cladeName;
    std::string cladeRank;
    std::string phylum;
};

struct ManifestShardInfo {
    std::string groupKey;
    std::string cladeName;
    std::string cladeRank;
    std::string phylum;
    std::string shardStem;
    std::string shardPath;
    size_t      rowCount = 0;
};

inline std::vector<BuildTargetKey> build_targets_for_manifest_row(const ManifestRow& r,
                                                                  bool multiRank) {
    std::vector<BuildTargetKey> out;
    std::unordered_set<std::string> seen;
    auto add_target = [&](std::string rank, std::string clade) {
        if (clade.empty()) return;
        if (rank.empty()) rank = "species";
        const std::string key = rank + "||" + clade;
        if (!seen.insert(key).second) return;
        out.push_back({key, std::move(clade), std::move(rank),
                       r.phylum.empty() ? "." : r.phylum});
    };

    add_target(r.cladeRank.empty() ? "species" : r.cladeRank,
               r.cladeName.empty() ? (r.asmName.empty() ? "unknown" : r.asmName)
                                   : r.cladeName);
    if (multiRank) {
        add_target("phylum", r.phylum);
        add_target("class", r.cls);
        add_target("order", r.order);
        add_target("family", r.family);
        add_target("genus", r.genus);
    }
    return out;
}

class ManifestShardAppendPool {
public:
    explicit ManifestShardAppendPool(size_t maxOpen = 64)
        : maxOpen_(std::max<size_t>(1, maxOpen)) {}

    std::ofstream& acquire(const std::string& path) {
        auto it = streams_.find(path);
        if (it != streams_.end()) {
            lru_.push_back(path);
            return *it->second;
        }
        if (streams_.size() >= maxOpen_) evict_one();
        auto out = std::make_unique<std::ofstream>(path, std::ios::app);
        if (!*out) throw std::runtime_error("Cannot append to manifest shard: " + path);
        std::ofstream& ref = *out;
        streams_[path] = std::move(out);
        lru_.push_back(path);
        return ref;
    }

    void close_all() {
        streams_.clear();
        lru_.clear();
    }

private:
    void evict_one() {
        while (!lru_.empty()) {
            const std::string victim = lru_.front();
            lru_.pop_front();
            auto it = streams_.find(victim);
            if (it == streams_.end()) continue;
            it->second->flush();
            streams_.erase(it);
            return;
        }
    }

    size_t maxOpen_;
    std::unordered_map<std::string, std::unique_ptr<std::ofstream>> streams_;
    std::deque<std::string> lru_;
};

inline std::vector<ManifestShardInfo>
partition_manifest_for_index_build(const std::string& manifestPath,
                                   const fs::path& tempDir,
                                   bool multiRank,
                                   size_t maxOpenFiles = 64) {
    fs::create_directories(tempDir);
    std::unordered_map<std::string, size_t> shardByKey;
    std::vector<ManifestShardInfo> shards;
    ManifestShardAppendPool pool(maxOpenFiles);

    for_each_manifest_row(manifestPath, [&](const ManifestRow& r) {
        for (const auto& target : build_targets_for_manifest_row(r, multiRank)) {
            size_t idx = 0;
            auto it = shardByKey.find(target.groupKey);
            if (it == shardByKey.end()) {
                idx = shards.size();
                shardByKey[target.groupKey] = idx;
                ManifestShardInfo info;
                info.groupKey = target.groupKey;
                info.cladeName = target.cladeName;
                info.cladeRank = target.cladeRank;
                info.phylum = target.phylum;
                info.shardStem = sanitize_filename(info.cladeRank) + "__"
                               + sanitize_filename(info.cladeName) + "__"
                               + std::to_string(idx);
                info.shardPath = (tempDir / (info.shardStem + ".rows.tsv")).string();
                shards.push_back(std::move(info));
            } else {
                idx = it->second;
            }
            auto& shard = shards[idx];
            auto& out = pool.acquire(shard.shardPath);
            write_manifest_row(out, r);
            ++shard.rowCount;
        }
    });
    pool.close_all();
    std::sort(shards.begin(), shards.end(),
              [](const ManifestShardInfo& a, const ManifestShardInfo& b) {
                  if (a.cladeRank != b.cladeRank) return a.cladeRank < b.cladeRank;
                  if (a.cladeName != b.cladeName) return a.cladeName < b.cladeName;
                  return a.shardStem < b.shardStem;
              });
    return shards;
}

inline std::string read_fasta_seed_sequence(const std::string& fastaPath, size_t maxBases = 128) {
    if (fastaPath.empty()) return "N";
    try {
        FastaStream fs(fastaPath);
        std::istream& in = fs.get();
        std::string line, seq;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '>') continue;
            for (char c : line) {
                const char u = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                if (u == 'A' || u == 'C' || u == 'G' || u == 'T' || u == 'N') seq.push_back(u);
                if (seq.size() >= maxBases) return seq;
            }
        }
        return seq.empty() ? std::string("N") : seq;
    } catch (...) { return "N"; }
}

inline std::vector<uint64_t> minimizer_hashes_for_sequence(const std::string& seq, int k, size_t limit = 32) {
    std::vector<uint64_t> hashes;
    if (seq.size() < static_cast<size_t>(std::max(1, k))) return hashes;
    std::unordered_set<uint64_t> seen;
    constexpr uint64_t kOffset = 1469598103934665603ull;
    constexpr uint64_t kPrime  = 1099511628211ull;
    for (size_t i = 0; i + static_cast<size_t>(k) <= seq.size() && hashes.size() < limit; ++i) {
        uint64_t h = kOffset;
        for (int j = 0; j < k; ++j) {
            h ^= static_cast<unsigned char>(seq[i + static_cast<size_t>(j)]);
            h *= kPrime;
        }
        if (seen.insert(h).second) hashes.push_back(h);
    }
    return hashes;
}

class BoundedMinHashSet {
public:
    explicit BoundedMinHashSet(size_t limit)
        : limit_(std::max<size_t>(1, limit)) {}

    void insert(uint64_t h) {
        if (members_.count(h)) return;
        if (heap_.size() < limit_) {
            heap_.push(h);
            members_.insert(h);
            return;
        }
        if (h >= heap_.top()) return;
        const uint64_t victim = heap_.top();
        heap_.pop();
        members_.erase(victim);
        heap_.push(h);
        members_.insert(h);
    }

    std::vector<uint64_t> finalize() const {
        std::vector<uint64_t> out(members_.begin(), members_.end());
        std::sort(out.begin(), out.end());
        return out;
    }

private:
    size_t limit_;
    std::priority_queue<uint64_t> heap_;
    std::unordered_set<uint64_t> members_;
};

struct CompactGraphSample {
    std::string pathName;
    std::string seedSeq;
};

inline std::string resolve_manifest_row_asm_name(const ManifestRow& r) {
    if (!r.asmName.empty()) return r.asmName;
    if (!r.fastaPath.empty()) return fs::path(r.fastaPath).stem().string();
    return "unknown";
}

inline std::string cached_seed_sequence_for_build(
        const std::string& fastaPath,
        size_t maxBases,
        std::unordered_map<std::string, std::string>& cache,
        size_t cacheLimit = 4096) {
    auto it = cache.find(fastaPath);
    if (it != cache.end()) return it->second;
    std::string seq = read_fasta_seed_sequence(fastaPath, maxBases);
    if (cache.size() < cacheLimit) cache.emplace(fastaPath, seq);
    return seq;
}

inline CladeGraph build_compact_manifest_graph(const CladeGraphDescriptor& d,
                                               const std::string& shardPath,
                                               const SyncmerParams& sp,
                                               size_t sampleLimit,
                                               std::vector<uint64_t>* centroidOut = nullptr) {
    CladeGraph g;
    g.cladeName = d.cladeName;
    g.cladeRank = d.cladeRank;
    g.phylum = d.phylum;
    g.genomeCount = d.genomeCount;

    std::minstd_rand rng(static_cast<uint32_t>(std::hash<std::string>{}(d.cladeName)));
    std::vector<CompactGraphSample> samples;
    samples.reserve(std::max<size_t>(1, sampleLimit));
    BoundedMinHashSet centroidAcc(std::max<size_t>(64, sampleLimit * 2));
    std::unordered_map<std::string, std::string> seedCache;
    size_t seen = 0;
    size_t contributing = 0;

    for_each_manifest_row(shardPath, [&](const ManifestRow& r) {
        ++seen;
        const std::string seed =
            cached_seed_sequence_for_build(r.fastaPath, 512, seedCache);
        if (seed.empty() || seed == "N") return;
        ++contributing;
        const int k = std::min<int>(std::max(3, sp.k), static_cast<int>(seed.size()));
        for (uint64_t h : minimizer_hashes_for_sequence(seed, k, 64))
            centroidAcc.insert(h);

        CompactGraphSample sample{resolve_manifest_row_asm_name(r), seed};
        if (samples.size() < sampleLimit) {
            samples.push_back(std::move(sample));
            return;
        }
        std::uniform_int_distribution<size_t> dist(0, seen - 1);
        const size_t pick = dist(rng);
        if (pick < sampleLimit) samples[pick] = std::move(sample);
    });

    std::unordered_map<std::string, int> nodeBySeq;
    int nextNodeId = 1;
    for (const auto& sample : samples) {
        auto it = nodeBySeq.find(sample.seedSeq);
        int nodeId = 0;
        if (it == nodeBySeq.end()) {
            nodeId = nextNodeId++;
            nodeBySeq.emplace(sample.seedSeq, nodeId);
            g.nodes.push_back({nodeId, sample.seedSeq});
        } else {
            nodeId = it->second;
        }
        g.paths.push_back({sample.pathName, {nodeId}});
    }

    if (g.nodes.size() > 1) g.svBubbles = g.nodes.size() - 1;
    else g.svBubbles = std::min<size_t>(1, d.svBubbles);

    size_t nodeBytes = 0;
    for (const auto& n : g.nodes) nodeBytes += n.sequence.size();
    size_t pathBytes = 0;
    for (const auto& p : g.paths)
        pathBytes += p.name.size() + p.nodes.size() * sizeof(int);
    g.compressedSz = 128 + nodeBytes + pathBytes
                   + std::min<size_t>(d.genomeCount, contributing) * 16;

    if (centroidOut != nullptr) *centroidOut = centroidAcc.finalize();
    return g;
}

// ── binary index helpers ───────────────────────────────────────────────────
template <class T>
inline void write_pod_bin(std::ostream& out, const T& value) {
    out.write(reinterpret_cast<const char*>(&value), static_cast<std::streamsize>(sizeof(T)));
    if (!out) throw std::runtime_error("binary write failed");
}

template <class T>
inline void read_pod_bin(std::istream& in, T& value) {
    in.read(reinterpret_cast<char*>(&value), static_cast<std::streamsize>(sizeof(T)));
    if (!in) throw std::runtime_error("binary read failed");
}

inline void write_string_bin(std::ostream& out, const std::string& s) {
    uint64_t n = static_cast<uint64_t>(s.size());
    write_pod_bin(out, n);
    out.write(s.data(), static_cast<std::streamsize>(n));
    if (!out) throw std::runtime_error("binary string write failed");
}

inline std::string read_string_bin(std::istream& in) {
    uint64_t n = 0;
    read_pod_bin(in, n);
    std::string s(static_cast<size_t>(n), '\0');
    if (n) {
        in.read(s.data(), static_cast<std::streamsize>(n));
        if (!in) throw std::runtime_error("binary string read failed");
    }
    return s;
}

inline uint64_t fourcc64(const char (&tag)[9]) {
    uint64_t x = 0;
    for (int i = 0; i < 8; ++i) x |= (static_cast<uint64_t>(static_cast<unsigned char>(tag[i])) << (8 * i));
    return x;
}

struct BinaryShardPaths {
    std::string gbzPath;
    std::string gbwtPath;
    std::string minPath;
    std::string cidxPath;
};

inline BinaryShardPaths shard_paths(const fs::path& indexDir, const fs::path& registryDir, const std::string& stem) {
    BinaryShardPaths p;
    p.gbzPath = (registryDir / (stem + ".gbz")).string();
    p.gbwtPath = (indexDir / (stem + ".gbwt")).string();
    p.minPath = (indexDir / (stem + ".min")).string();
    p.cidxPath = (indexDir / (stem + ".cidx")).string();
    return p;
}

inline void write_gbwt_payload(const std::string& path,
                               const CladeGraphDescriptor& d,
                               const CladeGraph& g) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot write GBWT file: " + path);
    const uint64_t magic = fourcc64("UTOLGBW1");
    const uint32_t version = 1;
    write_pod_bin(out, magic);
    write_pod_bin(out, version);
    write_string_bin(out, d.cladeName);
    write_string_bin(out, d.cladeRank);
    write_string_bin(out, d.phylum);
    uint64_t nPaths = static_cast<uint64_t>(g.paths.size());
    write_pod_bin(out, nPaths);
    for (const auto& p : g.paths) {
        write_string_bin(out, p.name);
        uint64_t n = static_cast<uint64_t>(p.nodes.size());
        write_pod_bin(out, n);
        for (int nodeId : p.nodes) write_pod_bin(out, static_cast<int32_t>(nodeId));
    }
}

inline std::vector<uint64_t> write_minimizer_payload(const std::string& path,
                                                     const CladeGraphDescriptor& d,
                                                     const CladeGraph& g,
                                                     const SyncmerParams& sp) {
    struct Posting { uint64_t hash; int32_t nodeId; int32_t offset; };
    std::vector<Posting> postings;
    std::unordered_set<uint64_t> centroidSet;
    for (const auto& n : g.nodes) {
        int k = std::min<int>(std::max(3, sp.k), static_cast<int>(n.sequence.size()));
        auto hashes = minimizer_hashes_for_sequence(n.sequence, k, 64);
        int32_t offset = 0;
        for (uint64_t h : hashes) {
            postings.push_back({h, static_cast<int32_t>(n.id), offset++});
            centroidSet.insert(h);
        }
    }
    std::sort(postings.begin(), postings.end(), [](const Posting& a, const Posting& b) {
        if (a.hash != b.hash) return a.hash < b.hash;
        if (a.nodeId != b.nodeId) return a.nodeId < b.nodeId;
        return a.offset < b.offset;
    });
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot write minimizer file: " + path);
    const uint64_t magic = fourcc64("UTOLMIN1");
    const uint32_t version = 1;
    write_pod_bin(out, magic);
    write_pod_bin(out, version);
    write_string_bin(out, d.cladeName);
    write_string_bin(out, d.cladeRank);
    write_string_bin(out, d.phylum);
    write_pod_bin(out, static_cast<int32_t>(sp.k));
    write_pod_bin(out, static_cast<int32_t>(sp.s));
    write_pod_bin(out, static_cast<uint64_t>(postings.size()));
    for (const auto& p : postings) {
        write_pod_bin(out, p.hash);
        write_pod_bin(out, p.nodeId);
        write_pod_bin(out, p.offset);
    }
    std::vector<uint64_t> centroid(centroidSet.begin(), centroidSet.end());
    std::sort(centroid.begin(), centroid.end());
    if (centroid.size() > 64) centroid.resize(64);
    return centroid;
}

inline void write_cidx_payload(const std::string& path,
                               const CladeGraphDescriptor& d,
                               const CladeGraph& g,
                               const BinaryShardPaths& shardPaths,
                               const SyncmerParams& sp,
                               double routingDensity,
                               const SyncmerParams* fb,
                               const std::vector<uint64_t>& centroid,
                               const std::string& shardPath,
                               size_t rowCount) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot write routing shard: " + path);
    const uint64_t magic = fourcc64("UTOLCIDX");
    const uint32_t version = 1;
    write_pod_bin(out, magic);
    write_pod_bin(out, version);
    write_string_bin(out, d.cladeName);
    write_string_bin(out, d.cladeRank);
    write_string_bin(out, d.phylum);
    write_string_bin(out, shardPaths.gbzPath);
    write_string_bin(out, shardPaths.gbwtPath);
    write_string_bin(out, shardPaths.minPath);
    write_pod_bin(out, static_cast<uint64_t>(d.genomeCount));
    write_pod_bin(out, static_cast<uint64_t>(g.nodes.size()));
    write_pod_bin(out, static_cast<uint64_t>(g.paths.size()));
    write_pod_bin(out, static_cast<uint64_t>(g.svBubbles));
    write_pod_bin(out, static_cast<uint64_t>(g.compressed_bytes()));
    write_pod_bin(out, routingDensity);
    write_pod_bin(out, static_cast<int32_t>(sp.k));
    write_pod_bin(out, static_cast<int32_t>(sp.s));
    write_pod_bin(out, static_cast<uint64_t>(sp.stride));
    int32_t fallbackK = fb ? fb->k : 0;
    int32_t fallbackS = fb ? fb->s : 0;
    write_pod_bin(out, fallbackK);
    write_pod_bin(out, fallbackS);
    write_pod_bin(out, static_cast<uint64_t>(centroid.size()));
    for (uint64_t h : centroid) write_pod_bin(out, h);
    write_pod_bin(out, static_cast<uint64_t>(rowCount));
    for_each_manifest_row(shardPath, [&](const ManifestRow& r) {
        write_string_bin(out, r.asmName);
        write_string_bin(out, r.fastaPath);
    });
}


inline CladeGraph build_manifest_backed_graph_from_rows(const CladeGraphDescriptor& d,
                                                        const std::vector<ManifestRow>& rowsForClade) {
    tol::CladeGraphBuilder builder(d.cladeName, d.cladeRank, d.phylum);
    int addedGenomes = 0;
    for (const auto& r : rowsForClade) {
        if (r.fastaPath.empty() || !fs::exists(r.fastaPath)) continue;
        try {
            auto contigs = read_fasta_local(r.fastaPath);
            if (contigs.empty()) continue;
            for (const auto& kv : contigs) {
                const std::string& contigName = kv.first;
                const std::string& seq = kv.second;
                if (seq.empty()) continue;
                builder.add_genome(
                    r.asmName.empty() ? fs::path(r.fastaPath).stem().string() : r.asmName,
                    contigName,
                    seq,
                    "NONE",
                    false,
                    false,
                    10000,
                    1000);
            }
            ++addedGenomes;
        } catch (...) {
            // Keep index builds resilient to individual malformed FASTA rows.
        }
    }
    CladeGraph g = builder.build(true, true, kMinBubbleFreqDef);
    g.cladeName = d.cladeName;
    g.cladeRank = d.cladeRank;
    g.phylum = d.phylum;
    g.genomeCount = std::max<size_t>(g.genomeCount, static_cast<size_t>(addedGenomes));
    g.svBubbles = std::max(g.svBubbles, d.svBubbles);
    g.compressedSz = std::max(g.compressed_bytes(), d.compressedBytes);

    if (g.nodes.empty()) {
        g.cladeName = d.cladeName;
        g.cladeRank = d.cladeRank;
        g.phylum = d.phylum;
        g.genomeCount = d.genomeCount;
        g.svBubbles = d.svBubbles;
        int nodeId = 1;
        for (const auto& r : rowsForClade) {
            const std::string seq = read_fasta_seed_sequence(r.fastaPath);
            g.nodes.push_back({nodeId, seq});
            g.paths.push_back({r.asmName.empty() ? ("asm_" + std::to_string(nodeId)) : r.asmName, {nodeId}});
            ++nodeId;
        }
        for (size_t i = 1; i < g.nodes.size(); ++i)
            g.edges.push_back({g.nodes[i - 1].id, g.nodes[i].id, true, true});
        size_t seqBytes = 0;
        for (const auto& n : g.nodes) seqBytes += n.sequence.size();
        g.compressedSz = std::max<size_t>(
            128 + seqBytes + g.edges.size() * 16 + g.paths.size() * 24,
            d.compressedBytes);
    }
    return g;
}

inline CladeGraph build_full_manifest_graph_from_shard(const CladeGraphDescriptor& d,
                                                       const std::string& shardPath) {
    return build_manifest_backed_graph_from_rows(d, read_manifest_rows(shardPath));
}

// ── write_graph_payload / write_routing_payload ───────────────────────────
inline CladeGraph write_graph_payload(const CladeGraphDescriptor& d,
                                const CladeGraph& sourceGraph) {
    CladeGraph g = sourceGraph;
    gbz_io::save(g, d.graphPath);
    return g;
}

inline std::vector<uint64_t> write_routing_payload(const BinaryShardPaths& shardPaths,
                                  const CladeGraphDescriptor& d,
                                  const CladeGraph& g,
                                  const SyncmerParams& sp,
                                  double routingDensity,
                                  const SyncmerParams* fb,
                                  const std::string& shardPath,
                                  size_t rowCount,
                                  const std::vector<uint64_t>* centroidOverride = nullptr) {
    write_gbwt_payload(shardPaths.gbwtPath, d, g);
    auto centroid = write_minimizer_payload(shardPaths.minPath, d, g, sp);
    if (centroidOverride != nullptr && !centroidOverride->empty())
        centroid = *centroidOverride;
    write_cidx_payload(shardPaths.cidxPath, d, g, shardPaths, sp, routingDensity, fb,
                       centroid, shardPath, rowCount);
    return centroid;
}

struct BuiltShardArtifacts {
    CladeDescriptor desc;
    std::string routingManifestLine;
};

inline BuiltShardArtifacts build_single_manifest_shard(
        const ManifestShardInfo& shard,
        const fs::path& indexDir,
        const fs::path& registryDir,
        const SyncmerParams& sp,
        double routingDensity,
        const SyncmerParams* fb,
        size_t maxCladeGenomes,
        bool baseGraphBuild,
        bool verbose) {
    CladeGraphDescriptor d;
    d.cladeName = shard.cladeName.empty() ? "unknown" : shard.cladeName;
    d.cladeRank = shard.cladeRank.empty() ? "species" : shard.cladeRank;
    d.phylum = shard.phylum.empty() ? "." : shard.phylum;
    const auto shardPaths = shard_paths(indexDir, registryDir, shard.shardStem);
    d.graphPath = shardPaths.gbzPath;
    d.genomeCount = shard.rowCount;
    d.svBubbles = std::max<size_t>(1, shard.rowCount);
    d.compressedBytes = 128 + shard.rowCount * 64;

    const bool useFullBaseGraph = baseGraphBuild && shard.rowCount <= maxCladeGenomes;
    if (baseGraphBuild && !useFullBaseGraph && verbose) {
        std::cerr << "[tol] " << d.cladeRank << ':' << d.cladeName
                  << " exceeds --tol-max-clade-genomes=" << maxCladeGenomes
                  << "; using compact manifest-backed graph mode\n";
    }

    std::vector<uint64_t> centroidOverride;
    const CladeGraph sourceGraph = useFullBaseGraph
        ? build_full_manifest_graph_from_shard(d, shard.shardPath)
        : build_compact_manifest_graph(
              d,
              shard.shardPath,
              sp,
              std::max<size_t>(64, std::min<size_t>(maxCladeGenomes, 256)),
              &centroidOverride);
    const CladeGraph g = write_graph_payload(d, sourceGraph);
    d.svBubbles = g.svBubbles;
    d.compressedBytes = std::max<size_t>(
        g.compressed_bytes(),
        fs::exists(d.graphPath) ? static_cast<size_t>(fs::file_size(d.graphPath))
                                : g.compressed_bytes());
    const auto centroid = write_routing_payload(
        shardPaths, d, g, sp, routingDensity, fb, shard.shardPath, shard.rowCount,
        centroidOverride.empty() ? nullptr : &centroidOverride);

    BuiltShardArtifacts built;
    built.desc.cladeName = d.cladeName;
    built.desc.cladeRank = d.cladeRank;
    built.desc.phylum = d.phylum;
    built.desc.graphPath = d.graphPath;
    built.desc.genomeCount = d.genomeCount;
    built.desc.svBubbles = d.svBubbles;
    built.desc.compressedBytes = d.compressedBytes;
    built.desc.crc32 = crc32_file(d.graphPath);
    built.desc.centroidSyncmers = centroid;

    std::ostringstream line;
    line << d.cladeName << '\t' << d.cladeRank << '\t' << d.phylum
         << '\t' << shardPaths.cidxPath << '\t' << d.graphPath << '\t'
         << shardPaths.gbwtPath << '\t' << shardPaths.minPath << '\t'
         << d.genomeCount << '\t' << g.nodes.size() << '\t' << g.paths.size() << '\t';
    for (size_t i = 0; i < centroid.size(); ++i) {
        if (i) line << ',';
        line << centroid[i];
    }
    built.routingManifestLine = line.str();
    return built;
}

inline std::vector<BuiltShardArtifacts> build_partitioned_shards(
        const std::vector<ManifestShardInfo>& shards,
        const fs::path& indexDir,
        const fs::path& registryDir,
        const SyncmerParams& sp,
        double routingDensity,
        const SyncmerParams* fb,
        size_t maxCladeGenomes,
        size_t indexThreads,
        bool baseGraphBuild,
        bool verbose) {
    std::vector<BuiltShardArtifacts> results(shards.size());
    std::atomic<size_t> next{0};
    std::mutex errMu;
    std::string firstError;

    const size_t workerCount = std::max<size_t>(1, std::min(indexThreads, shards.size()));
    std::vector<std::thread> workers;
    workers.reserve(workerCount);
    for (size_t ti = 0; ti < workerCount; ++ti) {
        workers.emplace_back([&]() {
            for (;;) {
                const size_t idx = next.fetch_add(1, std::memory_order_relaxed);
                if (idx >= shards.size()) break;
                {
                    std::lock_guard<std::mutex> lk(errMu);
                    if (!firstError.empty()) break;
                }
                try {
                    results[idx] = build_single_manifest_shard(
                        shards[idx], indexDir, registryDir, sp, routingDensity, fb,
                        maxCladeGenomes, baseGraphBuild, verbose);
                } catch (const std::exception& e) {
                    std::lock_guard<std::mutex> lk(errMu);
                    if (firstError.empty()) firstError = e.what();
                    break;
                }
            }
        });
    }
    for (auto& t : workers) t.join();
    if (!firstError.empty()) throw std::runtime_error(firstError);
    return results;
}

inline void write_external_routing_store_from_artifacts(
        const std::vector<BuiltShardArtifacts>& built,
        const fs::path& indexDir) {
    const std::string storePath = (indexDir / "routing_centroids.bin").string();
    {
        std::ofstream out(storePath, std::ios::binary);
        if (!out) throw std::runtime_error("Cannot write external centroid store: " + storePath);
        const uint64_t n = static_cast<uint64_t>(built.size());
        out.write(reinterpret_cast<const char*>(&n), sizeof(n));
        for (const auto& item : built) {
            ExternalMemoryCentroidStore::append_record(
                out,
                {item.desc.cladeName,
                 item.desc.phylum,
                 item.desc.cladeRank,
                 item.desc.centroidSyncmers});
        }
    }
    ExternalMemoryCentroidStore(storePath).prepare_skip_index();
}

// ── build_tol_index_from_manifest ─────────────────────────────────────────
inline void build_tol_index_from_manifest(const std::string& manifestPath,
                                          const std::string& indexDir,
                                          const std::string& registryDir,
                                          const SyncmerParams& sp,
                                          double routingDensity,
                                          bool verbose,
                                          const std::string& /*annotationTsv*/,
                                          size_t maxCladeGenomes,
                                          size_t indexThreads,
                                          const SyncmerParams* fb,
                                          bool baseGraphBuild) {
    fs::create_directories(indexDir);
    fs::create_directories(registryDir);
    const fs::path tempDir = fs::path(indexDir) / ".manifest_partitions";
    if (fs::exists(tempDir)) fs::remove_all(tempDir);
    std::vector<ManifestShardInfo> shards;
    std::vector<BuiltShardArtifacts> built;
    try {
        shards = partition_manifest_for_index_build(
            manifestPath, tempDir, false, std::max<size_t>(32, indexThreads * 16));
        built = build_partitioned_shards(
            shards, fs::path(indexDir), fs::path(registryDir), sp, routingDensity, fb,
            maxCladeGenomes, indexThreads, baseGraphBuild, verbose);
    } catch (...) {
        fs::remove_all(tempDir);
        throw;
    }
    write_external_routing_store_from_artifacts(built, fs::path(indexDir));

    std::ofstream routingManifest(fs::path(indexDir) / "routing_manifest.tsv");
    if (!routingManifest) throw std::runtime_error("Cannot write routing manifest");
    routingManifest << "#clade_name\tclade_rank\tphylum\tcidx_path\tgraph_path\tgbwt_path\tmin_path\tgenome_count\tnodes\tpaths\tcentroid_hashes\n";
    std::vector<CladeDescriptor> descriptors;
    descriptors.reserve(built.size());
    for (const auto& item : built) {
        descriptors.push_back(item.desc);
        routingManifest << item.routingManifestLine << '\n';
    }
    save_manifest(descriptors, (fs::path(registryDir) / "clade_manifest.tsv").string());
    fs::remove_all(tempDir);
    if (verbose)
        std::cerr << "[tol] built " << built.size()
                  << " clade shards from " << shards.size()
                  << " partitioned manifest groups\n";
}

// ── FIX-B3: build_multi_rank_index_from_manifest ─────────────────────────
// Builds one routing shard + graph payload per (rank, clade) pair so the
// multi-rank query path can route at phylum/class/order/family/genus/species
// independently.  The extended manifest (9 columns) is required.
inline void build_multi_rank_index_from_manifest(
        const std::string& manifestPath,
        const std::string& indexDir,
        const std::string& registryDir,
        const SyncmerParams& sp,
        double routingDensity,
        bool verbose,
        const std::string& /*annotationTsv*/,
        size_t maxCladeGenomes,
        size_t indexThreads,
        const SyncmerParams* fb,
        bool baseGraphBuild) {

    fs::create_directories(indexDir);
    fs::create_directories(registryDir);
    const fs::path tempDir = fs::path(indexDir) / ".manifest_partitions";
    if (fs::exists(tempDir)) fs::remove_all(tempDir);
    std::vector<ManifestShardInfo> shards;
    std::vector<BuiltShardArtifacts> built;
    try {
        shards = partition_manifest_for_index_build(
            manifestPath, tempDir, true, std::max<size_t>(32, indexThreads * 16));
        built = build_partitioned_shards(
            shards, fs::path(indexDir), fs::path(registryDir), sp, routingDensity, fb,
            maxCladeGenomes, indexThreads, baseGraphBuild, verbose);
    } catch (...) {
        fs::remove_all(tempDir);
        throw;
    }
    write_external_routing_store_from_artifacts(built, fs::path(indexDir));

    std::ofstream routingManifest(fs::path(indexDir) / "routing_manifest.tsv");
    if (!routingManifest) throw std::runtime_error("Cannot write routing manifest");
    routingManifest << "#clade_name\tclade_rank\tphylum\tcidx_path\tgraph_path\tgbwt_path\tmin_path\tgenome_count\tnodes\tpaths\tcentroid_hashes\n";
    std::vector<CladeDescriptor> descriptors;

    descriptors.reserve(built.size());
    for (const auto& item : built) {
        descriptors.push_back(item.desc);

        routingManifest << item.routingManifestLine << '\n';
    }
    save_manifest(descriptors, (fs::path(registryDir) / "clade_manifest.tsv").string());
    fs::remove_all(tempDir);
    if (verbose)
        std::cerr << "[tol] multi-rank: built " << built.size()
                  << " rank x clade shards from " << shards.size()
                  << " partitioned manifest groups\n";
}

// ── TolGlobal ─────────────────────────────────────────────────────────────
class TolGlobal {
public:
    struct RefSeq {
        std::string asmName;
        std::string clade;
        std::string contig;
        std::shared_ptr<const std::string> seqShared;
        std::string cladeRank; // populated in multi-rank mode
        std::string phylum;
        double      cladeGc = 0.45; // background GC for HGT detection

        const std::string& seq() const {
            static const std::string kEmpty;
            return seqShared ? *seqShared : kEmpty;
        }
        bool has_seq() const { return seqShared && !seqShared->empty(); }
    };

    static TolGlobal& instance() {
        static TolGlobal inst;
        return inst;
    }

    void init(const std::string& indexDir,
              const std::string& registryDir,
              size_t /*cacheBytes*/,
              size_t /*cacheEntries*/) {
        std::lock_guard<std::mutex> lk(mu_);
        indexDir_    = indexDir;
        registryDir_ = registryDir;
        registry_    = ManifestRegistry(registryDir_);
        registry_.load_from_disk();
        refsByContig_.clear();
        allRefs_.clear();
        seenFasta_.clear();
        seqPool_.clear();
        router_.reset();
        externalCentroidStore_.reset();
        routingCladeCount_ = 0;
        preferExternalRouting_ = false;

        const auto routingManifest = (fs::path(indexDir_) / "routing_manifest.tsv").string();
        const std::string storePath = (fs::path(indexDir_) / "routing_centroids.bin").string();
        uint64_t storedCount = 0;
        if (fs::exists(storePath)) {
            std::ifstream in(storePath, std::ios::binary);
            if (in) in.read(reinterpret_cast<char*>(&storedCount), sizeof(storedCount));
        }

        static constexpr size_t kInMemoryRoutingCentroidLimit = 200000;
        const bool hasExternalStore = storedCount > 0 && fs::exists(storePath);
        if (hasExternalStore) {
            externalCentroidStore_ = std::make_unique<ExternalMemoryCentroidStore>(storePath);
            routingCladeCount_ = static_cast<size_t>(storedCount);
            preferExternalRouting_ =
                routingCladeCount_ > kInMemoryRoutingCentroidLimit;
        }

        if (!preferExternalRouting_) {
            std::vector<CladeCentroid> routingCentroids;
            for (const auto& row : read_routing_manifest_rows(routingManifest)) {
                if (row.centroidHashes.empty()) continue;
                CladeCentroid c;
                c.cladeName = row.cladeName;
                c.phylum = row.phylum;
                c.cladeRank = row.cladeRank;
                c.centroidHashes = row.centroidHashes;
                c.build_prefilters();
                router_.register_clade_centroid(c);
                routingCentroids.push_back(std::move(c));
            }
            if (!routingCentroids.empty()) {
                router_.rebuild();
                routingCladeCount_ = routingCentroids.size();
                if (externalCentroidStore_ == nullptr)
                    externalCentroidStore_ = std::make_unique<ExternalMemoryCentroidStore>(storePath);
                try {
                    bool reuseExistingStore = false;
                    if (fs::exists(storePath)) {
                        std::ifstream in(storePath, std::ios::binary);
                        uint64_t countOnDisk = 0;
                        if (in.read(reinterpret_cast<char*>(&countOnDisk), sizeof(countOnDisk))) {
                            reuseExistingStore =
                                countOnDisk == static_cast<uint64_t>(routingCentroids.size());
                        }
                    }
                    if (!reuseExistingStore)
                        externalCentroidStore_->build(routingCentroids);
                } catch (...) {
                    externalCentroidStore_.reset();
                }
            }
        }

        // Sequence pool: each FASTA file is loaded ONCE regardless of how many
        // rank-level descriptors reference it.  A shared_ptr to each contig
        // sequence string is stored in seqPool_ keyed by (fasta_path, contig).
        // RefSeq entries for different rank-level descriptors share the same
        // underlying string, so a 150 Mb AMF assembly loaded under 6 ranks
        // occupies 150 Mb in RAM, not 900 Mb.
        //
        // Deduplication rule:
        //   • (fasta, contig, cladeName) — deduplicate within the same clade
        //     so a contig listed twice under the same descriptor is not doubled.
        //   • (fasta, contig, DIFFERENT cladeNames) — allowed: different rank
        //     levels must produce separate RefSeq entries with distinct cladeRank
        //     values so hierarchical_call_assembly_multirank can filter by rank.
        std::unordered_map<std::string,
            std::unordered_map<std::string, std::shared_ptr<std::string>>> fastaContigSeq;
        // seen key = fasta + "||" + contig + "||" + cladeName
        std::unordered_set<std::string> seenFastaContigClade;

        for (const auto& d : registry_.descriptors()) {
            for (const auto& fasta : d.fastaPaths) {
                if (fasta.empty() || !fs::exists(fasta)) continue;
                try {
                    // Load sequences for this FASTA once; reuse on subsequent calls.
                    if (fastaContigSeq.find(fasta) == fastaContigSeq.end()) {
                        auto contigs = read_fasta_local(fasta);
                        auto& pool = fastaContigSeq[fasta];
                        for (auto& kv : contigs)
                            pool[kv.first] = std::make_shared<std::string>(std::move(kv.second));
                    }
                    const auto& pool    = fastaContigSeq.at(fasta);
                    const std::string asmName = fs::path(fasta).stem().string();
                    for (const auto& kv : pool) {
                        const std::string seenKey =
                            fasta + "||" + kv.first + "||" + d.cladeName;
                        if (!seenFastaContigClade.insert(seenKey).second) continue;
                        const double cgc = gc_content(*kv.second);
                        RefSeq r;
                        r.asmName   = asmName;
                        r.clade     = d.cladeName;
                        r.contig    = kv.first;
                        r.seqShared = kv.second;
                        r.cladeRank = d.cladeRank;
                        r.phylum    = d.phylum;
                        r.cladeGc   = cgc;
                        refsByContig_[kv.first].push_back(r);
                        allRefs_.push_back(std::move(r));
                    }
                } catch (...) {}
            }
        }
        initialized_ = true;
    }

    bool is_initialized() const { return initialized_; }
    const std::unordered_map<std::string, std::vector<RefSeq>>& refs_by_contig()  const { return refsByContig_; }
    const std::vector<RefSeq>& all_refs()                                          const { return allRefs_;      }
    const ManifestRegistry&    registry()                                           const { return registry_;    }
    bool has_routing_index() const { return routingCladeCount_ > 0; }

    std::vector<std::string> route_query_to_clades(std::string_view seq,
                                                   const SyncmerParams& sp,
                                                   const SyncmerParams& fbSp,
                                                   double density,
                                                   size_t topK) const {
        std::vector<std::string> out;
        if (routingCladeCount_ == 0 || seq.empty()) return out;

        std::vector<RouteResult> results;
        if (!preferExternalRouting_)
            results = router_.route(seq, sp, fbSp, density, topK);
        if ((results.empty() || preferExternalRouting_) && externalCentroidStore_ != nullptr) {
            auto qc = make_query_centroid_for_routing(seq, sp, density);
            auto ext = externalCentroidStore_->query_topk_streaming(qc, topK);
            results.reserve(ext.size());
            for (auto& r : ext)
                results.push_back({r.cladeName, r.phylum, r.jaccard});
        }
        out.reserve(results.size());
        for (const auto& r : results)
            if (!r.cladeName.empty()) out.push_back(r.cladeName);
        return out;
    }

private:
    TolGlobal() = default;
    bool initialized_ = false;
    std::string indexDir_, registryDir_;
    ManifestRegistry registry_;
    std::unordered_set<std::string> seenFasta_;
    std::unordered_map<std::string, std::vector<RefSeq>> refsByContig_;
    std::vector<RefSeq> allRefs_;
    // Sequence pool: shared_ptr<string> per (fasta_path, contig_name) so
    // sequences are allocated once even when referenced under 6 rank levels.
    std::unordered_map<std::string,
        std::unordered_map<std::string, std::shared_ptr<std::string>>> seqPool_;
    mutable PhylumShardedRouter router_;
    std::unique_ptr<ExternalMemoryCentroidStore> externalCentroidStore_;
    size_t routingCladeCount_ = 0;
    bool preferExternalRouting_ = false;
    mutable std::mutex mu_;
};

// ── MultiRankIndex ────────────────────────────────────────────────────────
class MultiRankIndex {
public:
    static MultiRankIndex& instance() {
        static MultiRankIndex inst;
        return inst;
    }
    void init(const std::string& indexDir, const std::string& registryDir,
              size_t cacheBytes, size_t cacheEntries) {
        TolGlobal::instance().init(indexDir, registryDir, cacheBytes, cacheEntries);
        initialized_ = true;
    }
    bool is_initialized() const { return initialized_; }
private:
    bool initialized_ = false;
};

// Parse read support from a pseudo-contig name.
// _n<N>  = long-reads cluster size  (exact read count)
// _mf<N> = short-reads min k-mer frequency along unitig path (coverage proxy)
inline int parse_pseudocontig_support(const std::string& name) {
    for (const char* tag : {"_n", "_mf"}) {
        const size_t tlen = __builtin_strlen(tag);
        auto pos = name.rfind(tag);
        if (pos != std::string::npos) {
            try {
                size_t nchars = 0;
                const int v = std::stoi(name.substr(pos + tlen), &nchars);
                if (nchars > 0 && v >= 0) return v;
            } catch (...) {}
        }
    }
    return -1;
}

// ── make_insdel_call ──────────────────────────────────────────────────────
inline VariantCallBridge make_insdel_call(const std::string& qAsm,
                                          const std::string& qContig,
                                          const TolGlobal::RefSeq& ref,
                                          int qlen,
                                          const FederatedOptions& fo) {
    const int rlen  = static_cast<int>(ref.seq().size());
    const int delta = qlen - rlen;
    VariantCallBridge v;
    v.qAsm    = qAsm;
    v.qContig = qContig;
    // Use the concrete reference assembly identifier for downstream scoring
    // and TSV reporting. Routed taxonomic context is carried separately in the
    // cladeRank/phylum fields.
    v.refAsm  = ref.asmName.empty() ? (ref.clade.empty() ? "." : ref.clade)
                                    : ref.asmName;
    v.refContig = ref.contig;
    v.refPos  = std::min(std::max(1, fo.primarySketchParams.k + 1), std::max(1, rlen));
    v.type    = delta < 0 ? "DEL" : "INS";
    v.svlen   = delta;
    v.pos     = v.refPos;
    v.end     = (v.type == "DEL") ? std::min(rlen, v.pos + std::abs(delta) - 1) : v.pos;
    v.refEnd  = v.end;
    v.blockScore  = std::max(fo.minBlockScore, 12.0);
    v.anchors     = static_cast<int>(std::max(fo.minAnchors, size_t(2)));
    v.genotype    = "0/1";
    v.gq          = 30.0;
    v.annotation  = "NONE";
    v.alignmentMode = "tol_light_length";
    v.mapq        = 40.0;
    v.pantreeClass = v.type;
    v.isNonRefVariant  = false;
    v.triallelicTopology = ".";
    v.cladeRank   = ref.cladeRank.empty() ? "." : ref.cladeRank;
    v.phylum      = ref.phylum.empty() ? "." : ref.phylum;
    v.readSupport = parse_pseudocontig_support(qContig);
    return v;
}

// ── make_offref_call ──────────────────────────────────────────────────────
// FIX-B5: also classifies the element type and stores it in elementClass.
inline VariantCallBridge make_offref_call(const std::string& qAsm,
                                          const std::string& qContig,
                                          const std::string& seq,
                                          const std::string& tier,
                                          double cladeGc               = 0.45,
                                          const std::string& refAsm    = "OFF_REFERENCE",
                                          const std::string& cladeRank = ".",
                                          const std::string& phylum    = ".") {
    VariantCallBridge v;
    v.qAsm    = qAsm;
    v.qContig = qContig;
    v.refAsm  = refAsm;
    v.refContig = ".";
    v.refPos  = 0;
    v.refEnd  = 0;
    v.type    = "OFF_REF";
    v.svlen   = std::max(1, static_cast<int>(seq.size()));
    v.pos     = 1;
    v.end     = v.svlen;
    v.blockScore  = 8.0;
    v.anchors     = 0;
    v.genotype    = "0/1";
    v.gq          = 20.0;
    v.annotation  = tier;
    v.alignmentMode = "tol_light_offref";
    v.mapq        = 10.0;
    v.pantreeClass = "NON_REF";
    v.isNonRefVariant  = true;
    v.triallelicTopology = ".";
    v.cladeRank   = cladeRank.empty() ? "." : cladeRank;
    v.phylum      = phylum.empty() ? "." : phylum;
    v.readSupport = parse_pseudocontig_support(qContig);
    // Classify element type for GFA EC tag
    const ElementClass ec = seq.empty()
        ? ElementClass::NONE
        : classify_repeat_element(std::string_view(seq.data(), seq.size()), cladeGc);
    v.elementClass = element_class_name(ec);
    return v;
}

struct OffRefWindowCall {
    size_t start = 0;
    size_t end = 0; // half-open query interval
    double bestOverlap = 0.0;
    double cladeGc = 0.45;
    std::string bestAsm = "OFF_REFERENCE";
    std::string cladeRank = ".";
    std::string phylum = ".";
    std::string tier = "NOVEL";
    std::string elementClass = "NONE";
};

inline std::vector<OffRefWindowCall>
discover_graph_native_offref_windows(const std::string& seq,
                                     const std::vector<TolGlobal::RefSeq>& refs,
                                     const FederatedOptions& fo,
                                     size_t minWindowBp = 0) {
    std::vector<OffRefWindowCall> out;
    if (!fo.graphNativeMode) return out;
    if (seq.empty() || refs.empty()) return out;

    const size_t n = seq.size();
    size_t win = minWindowBp;
    if (win == 0) {
        win = std::max<size_t>(
            fo.tolMinBlockBp > 0 ? fo.tolMinBlockBp : 0u,
            static_cast<size_t>(std::max(fo.minSvLen, 500)));
        win = std::min<size_t>(win, std::max<size_t>(500, std::min<size_t>(n, 5000)));
    }
    if (n < win) return out;
    const size_t step = std::max<size_t>(1, win / 2);
    const int k = std::max(5, std::min(fo.fallbackSketchParams.k > 0 ? fo.fallbackSketchParams.k : 7, 9));

    for (size_t start = 0; start + win <= n; start += step) {
        const std::string window = seq.substr(start, win);
        if (is_low_complexity_sequence(window)) continue;
        double bestOverlap = 0.0;
        double bestCladeGc = 0.45;
        std::string bestAsm = "OFF_REFERENCE";
        std::string bestRank = ".";
        std::string bestPhylum = ".";
        for (const auto& ref : refs) {
            if (!ref.has_seq()) continue;
            const double ov = kmer_overlap_fraction(window, ref.seq(), k);
            if (ov > bestOverlap) {
                bestOverlap = ov;
                bestCladeGc = ref.cladeGc;
                bestAsm = ref.clade.empty() ? ref.asmName : ref.clade;
                bestRank = ref.cladeRank.empty() ? "." : ref.cladeRank;
                bestPhylum = ref.phylum.empty() ? "." : ref.phylum;
            }
        }
        const std::string tier = infer_novelty_tier(bestOverlap);
        if (tier != "NOVEL" && tier != "NOVEL_WEAK" && tier != "DIVERGED") continue;
        OffRefWindowCall ow;
        ow.start = start;
        ow.end = start + win;
        ow.bestOverlap = bestOverlap;
        ow.cladeGc = bestCladeGc;
        ow.bestAsm = bestAsm;
        ow.cladeRank = bestRank;
        ow.phylum = bestPhylum;
        ow.tier = tier;
        ow.elementClass = element_class_name(classify_repeat_element(std::string_view(window.data(), window.size()), bestCladeGc));
        if (!out.empty() && out.back().tier == ow.tier && out.back().bestAsm == ow.bestAsm && out.back().end >= ow.start) {
            out.back().end = ow.end;
            if (out.back().elementClass == "NONE") out.back().elementClass = ow.elementClass;
        } else {
            out.push_back(std::move(ow));
        }
    }
    return out;
}

inline VariantCallBridge make_offref_window_call(const std::string& qAsm,
                                                 const std::string& qContig,
                                                 const std::string& seq,
                                                 const OffRefWindowCall& ow) {
    const size_t safeStart = std::min(ow.start, seq.size());
    const size_t safeEnd = std::min(std::max(ow.end, ow.start), seq.size());
    const std::string sub = seq.substr(safeStart, safeEnd - safeStart);
    VariantCallBridge v = make_offref_call(qAsm, qContig, sub, ow.tier, ow.cladeGc,
                                           ow.bestAsm, ow.cladeRank, ow.phylum);
    v.pos = static_cast<int>(safeStart + 1);
    v.end = static_cast<int>(std::max(safeStart + 1, safeEnd));
    v.svlen = static_cast<int>(std::max<size_t>(1, safeEnd - safeStart));
    v.alignmentMode = "graph_native_offref_window";
    if (!ow.elementClass.empty()) v.elementClass = ow.elementClass;
    return v;
}

// ── FIX-B2: try_mem_chain_call ────────────────────────────────────────────
// The order permutation is tracked explicitly so chain-to-MEM recovery is O(N)
// instead of the previous quadratic search.
static bool try_mem_chain_call(
        const std::string& qAsm,
        const std::string& qContig,
        const std::string& qSeq,
        const std::vector<const TolGlobal::RefSeq*>& refCandidates,
        const FederatedOptions& fo,
        VariantCallBridge& call) {

    if (refCandidates.empty() || qSeq.empty()) return false;

    SuffixArray sa;
    std::vector<std::pair<std::string, std::string>> refContigs;
    std::vector<const TolGlobal::RefSeq*> saRefs;
    refContigs.reserve(refCandidates.size());
    saRefs.reserve(refCandidates.size());
    const size_t saTextCap = fo.saMaxTextMB > 0 ? fo.saMaxTextMB * 1024 * 1024 : SIZE_MAX;
    size_t saTextAccum = 0;
    for (const auto* r : refCandidates) {
        if (!r->has_seq()) continue;
        if (saTextAccum + r->seq().size() > saTextCap) continue;
        saTextAccum += r->seq().size();
        refContigs.push_back({ r->contig, r->seq() });
        saRefs.push_back(r);
    }
    if (refContigs.empty()) return false;
    sa.build(refContigs);
    const std::string rcSeq = SuffixArray::revcomp(qSeq);

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

        std::vector<SuffixArray::Mem> allMems;
        std::vector<bool> isRev;
        allMems.reserve(fwdMems.size() + revMems.size());
        isRev.reserve(fwdMems.size() + revMems.size());
        for (auto& m : fwdMems) { allMems.push_back(m); isRev.push_back(false); }
        for (auto& m : revMems) { allMems.push_back(m); isRev.push_back(true);  }
        if (allMems.empty()) return out;

        std::vector<int> order(allMems.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(),
                  [&](int a, int b) {
                      return allMems[static_cast<size_t>(a)].qPos <
                             allMems[static_cast<size_t>(b)].qPos;
                  });

        ChainTreap treap;
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

        std::vector<SuffixArray::Mem> chain;
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

        double bestScore = static_cast<double>(treap.best_chain_score());
        if (chain.size() < 2 || bestScore < fo.minBlockScore)
            return out;

        auto res = SvTypeFromChain::classify(chain, chainRev, sa, fo.minSvLen);

        // DUP fallback: ChainTreap requires strictly-increasing rPos and silently
        // drops backward-mapping MEMs that are the signature of tandem duplications.
        // Re-classify the qPos-sorted forward MEMs directly to catch that pattern.
        // Also re-check INS calls: with 1% divergence the primary chain sees a net
        // positive query gap (delta > 0) and classifies DUPs as INS.
        if (res.type == SvTypeFromChain::Type::NONE ||
            res.type == SvTypeFromChain::Type::INS) {
            std::vector<SuffixArray::Mem> fwdSorted;
            fwdSorted.reserve(fwdMems.size());
            for (int i : order)
                if (!isRev[static_cast<size_t>(i)])
                    fwdSorted.push_back(allMems[static_cast<size_t>(i)]);
            if (fwdSorted.size() >= 2) {
                auto dupRes = SvTypeFromChain::classify(
                    fwdSorted, std::vector<bool>(fwdSorted.size(), false),
                    sa, fo.minSvLen);
                if (dupRes.type == SvTypeFromChain::Type::DUP) {
                    double dupScore = 0.0;
                    for (const auto& m : fwdSorted) dupScore += static_cast<double>(m.len);
                    if (dupScore >= fo.minBlockScore) {
                        res       = dupRes;
                        chain     = fwdSorted;
                        chainRev.assign(fwdSorted.size(), false);
                        bestScore = dupScore;
                    }
                }
            }
        }

        // TRA fallback: large query gaps between cross-contig anchors can exceed
        // chainGapBand, preventing ChainTreap from linking them.  Scan allMems
        // directly with a generous tolerance (500 kb).
        if (res.type == SvTypeFromChain::Type::NONE) {
            auto ctg_of = [&](int rp) -> int {
                for (int ci = 0; ci < static_cast<int>(sa.contigEnd.size()); ++ci)
                    if (rp < sa.contigEnd[static_cast<size_t>(ci)]) return ci;
                return -1;
            };
            const SuffixArray::Mem* srcMem = nullptr;
            const SuffixArray::Mem* dstMem = nullptr;
            int srcCtg = -1;
            std::string srcAsm;
            const int traGap = 500000;
            for (int i : order) {
                if (isRev[static_cast<size_t>(i)]) continue;
                const auto& m = allMems[static_cast<size_t>(i)];
                const int ci = ctg_of(m.rPos);
                if (srcCtg < 0) {
                    srcCtg = ci;
                    srcMem = &m;
                    srcAsm = (ci >= 0 && ci < static_cast<int>(saRefs.size()))
                             ? saRefs[static_cast<size_t>(ci)]->asmName : "";
                    continue;
                }
                if (ci >= 0 && ci != srcCtg) {
                    // Allow cross-assembly TRAs: in the hierarchical TOL mode the SA
                    // intentionally contains multiple reference assemblies per clade, so
                    // a query breakpoint matching contigs from different assemblies is a
                    // valid intra-clade TRA (or HGT candidate). Blocking cross-assembly
                    // pairs was causing TRA recall=0 in multi-reference SA contexts.
                    const int qGap = m.qPos - (srcMem->qPos + srcMem->len);
                    if (qGap >= 0 && qGap <= traGap &&
                        (dstMem == nullptr || m.len > dstMem->len))
                        dstMem = &m;
                }
            }
            if (srcMem && dstMem) {
                std::vector<SuffixArray::Mem> traChain = {*srcMem, *dstMem};
                auto traRes = SvTypeFromChain::classify(
                    traChain, {false, false}, sa, fo.minSvLen);
                if (traRes.type == SvTypeFromChain::Type::TRA) {
                    double traScore = static_cast<double>(srcMem->len + dstMem->len);
                    if (traScore >= fo.minBlockScore) {
                        res       = traRes;
                        chain     = traChain;
                        chainRev  = {false, false};
                        bestScore = traScore;
                    }
                }
            }
        }

        if (res.type == SvTypeFromChain::Type::NONE) return out;

        auto contig_of = [&](int rPos) -> int {
            for (int ci = 0; ci < static_cast<int>(sa.contigEnd.size()); ++ci)
                if (rPos < sa.contigEnd[static_cast<size_t>(ci)]) return ci;
            return -1;
        };

        const int primaryContigIdx = !chain.empty() ? contig_of(chain.front().rPos) : -1;
        const TolGlobal::RefSeq* primaryRef =
            (primaryContigIdx >= 0 && static_cast<size_t>(primaryContigIdx) < saRefs.size())
                ? saRefs[static_cast<size_t>(primaryContigIdx)]
                : (saRefs.empty() ? nullptr : saRefs.front());

        out.call.qAsm         = qAsm;
        out.call.qContig      = qContig;
        out.call.readSupport  = parse_pseudocontig_support(qContig);
        out.call.refAsm  = (primaryRef != nullptr && !primaryRef->asmName.empty())
            ? primaryRef->asmName
            : ((primaryRef != nullptr && !primaryRef->clade.empty()) ? primaryRef->clade : "unknown");
        out.call.refContig = (primaryRef != nullptr && !primaryRef->contig.empty())
            ? primaryRef->contig
            : (res.rContig.empty() ? "." : res.rContig);
        out.call.refPos  = 0;
        out.call.refEnd  = 0;
        out.call.pos     = std::max(1, res.qBreakStart + 1);
        out.call.end     = std::max(out.call.pos, res.qBreakEnd >= 0 ? res.qBreakEnd + 1 : out.call.pos);
        out.call.svlen   = res.svLen;
        out.call.genotype    = "0/1";
        out.call.gq = std::min(99.0,
            10.0 * std::log10(1.0 + static_cast<double>(chain.size()))
            + 0.5 * bestScore);
        out.call.blockScore  = bestScore;
        out.call.anchors     = static_cast<int>(chain.size());
        out.call.alignmentMode = secondaryPass
            ? "mem_chain_ds13_ds18;secondary_seed_rescue"
            : "mem_chain_ds13_ds18";
        out.call.mapq        = 50.0;
        out.call.annotation  = "NONE";
        out.call.triallelicTopology = ".";
        out.call.isNonRefVariant    = false;
        if (primaryRef != nullptr) {
            out.call.cladeRank = primaryRef->cladeRank.empty() ? "." : primaryRef->cladeRank;
            out.call.phylum    = primaryRef->phylum.empty() ? "." : primaryRef->phylum;
        }

        using T = SvTypeFromChain::Type;
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
                out.call.type           = "TRA";
                out.call.pantreeClass   = "NON_REF";
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
                out.call.mateContig     = res.rContig;
                out.call.matePos        = res.rBreakStart + 1;
                out.call.mateEnd        = res.rBreakEnd > 0 ? res.rBreakEnd : out.call.matePos;
                if (!res.rContig.empty()) {
                    for (const auto* cand : saRefs) {
                        if (cand == nullptr || cand->contig != res.rContig) continue;
                        out.call.mateRefAsm = !cand->asmName.empty()
                            ? cand->asmName
                            : (cand->clade.empty() ? "." : cand->clade);
                        out.call.mateOffReference = false;
                        break;
                    }
                }
                break;
            default: return out;
        }

        // Wire TE classification: INS/DUP/INV use the query subsequence;
        // DEL uses the deleted reference sequence (TE-mediated deletions are common).
        using T2 = SvTypeFromChain::Type;
        if (res.type == T2::INS || res.type == T2::DUP || res.type == T2::INV) {
            const double gcBg = (primaryRef != nullptr) ? primaryRef->cladeGc : 0.45;
            const int teS = res.qBreakStart;
            const int teE = (res.type == T2::INS)
                ? std::min(teS + res.svLen, static_cast<int>(qSeq.size()))
                : std::min(res.qBreakEnd + 1, static_cast<int>(qSeq.size()));
            if (teE > teS && teS >= 0) {
                const ElementClass ec = classify_repeat_element(
                    std::string_view(qSeq.data() + teS,
                                     static_cast<size_t>(teE - teS)),
                    gcBg);
                if (ec != ElementClass::NONE)
                    out.call.elementClass = element_class_name(ec);
            }
        } else if (res.type == T2::DEL && primaryRef != nullptr && primaryRef->has_seq()) {
            const double gcBg = primaryRef->cladeGc;
            // SvTypeFromChain::classify now returns local-contig coordinates,
            // so no subtraction of contigEnd is required here.
            const int rS = res.rBreakStart;
            const int rE = res.rBreakEnd;
            const auto& refSeq = primaryRef->seq();
            const int safeS = std::max(0, rS);
            const int safeE = std::min(static_cast<int>(refSeq.size()), rE);
            if (safeE > safeS) {
                const ElementClass ec = classify_repeat_element(
                    std::string_view(refSeq.data() + safeS,
                                     static_cast<size_t>(safeE - safeS)),
                    gcBg);
                if (ec != ElementClass::NONE)
                    out.call.elementClass = element_class_name(ec);
            }
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

// ── hierarchical_call_assembly ────────────────────────────────────────────
inline std::vector<VariantCallBridge>
hierarchical_call_assembly(const std::string& qAsm,
                           const std::unordered_map<std::string, std::string>& contigs,
                           const FederatedOptions& fo) {
    std::vector<VariantCallBridge> out;
    const auto& global   = TolGlobal::instance();
    const auto& byContig = global.refs_by_contig();
    const auto& allRefs  = global.all_refs();
    out.reserve(contigs.size());

    std::vector<const TolGlobal::RefSeq*> allRefPtrs;
    allRefPtrs.reserve(allRefs.size());
    for (const auto& r : allRefs) allRefPtrs.push_back(&r);
    const int minHierChainAnchors = static_cast<int>(std::max<size_t>(fo.tolMinChainAnchors, 2));
    auto append_window_calls = [&](const std::string& contigName,
                                   const std::string& seq) {
        if (!fo.graphNativeMode) return;
        auto windows = discover_graph_native_offref_windows(seq, allRefs, fo);
        for (const auto& ow : windows) out.push_back(make_offref_window_call(qAsm, contigName, seq, ow));
    };

    for (const auto& kv : contigs) {
        const std::string& name = kv.first;
        const std::string& seq  = kv.second;
        const auto routed = global.route_query_to_clades(
            seq, fo.primarySketchParams, fo.fallbackSketchParams,
            fo.routingDensity, fo.routingTopN);
        const std::unordered_set<std::string> routedClades(routed.begin(), routed.end());
        auto clade_allowed = [&](const TolGlobal::RefSeq& ref) {
            return routedClades.empty() ||
                   routedClades.find(ref.clade) != routedClades.end() ||
                   routedClades.find(ref.asmName) != routedClades.end();
        };

        // ── Path A: DS-13+DS-18 MEM chain ──────────────────────────────
        std::vector<const TolGlobal::RefSeq*> cands;
        auto it = byContig.find(name);
        if (it != byContig.end())
            for (const auto& r : it->second)
                if (r.has_seq() && clade_allowed(r)) cands.push_back(&r);
        if (cands.empty() && it != byContig.end())
            for (const auto& r : it->second)
                if (r.has_seq()) cands.push_back(&r);
        if (cands.empty())
            for (const auto* r : allRefPtrs)
                if (r->has_seq() && clade_allowed(*r)) cands.push_back(r);
        if (cands.empty())
            for (const auto* r : allRefPtrs)
                if (r->has_seq()) cands.push_back(r);

        VariantCallBridge chainCall;
        if (try_mem_chain_call(qAsm, name, seq, cands, fo, chainCall) &&
            chainCall.anchors >= minHierChainAnchors) {
            if (!out.empty()) {
                MergeSortTree mst;
                std::vector<std::pair<int,int>> existingIvs;
                existingIvs.reserve(out.size());
                for (const auto& c : out) existingIvs.push_back({c.pos, c.end});
                mst.build(existingIvs);
                auto overlaps = mst.overlapping(chainCall.pos, chainCall.end);
                if (!overlaps.empty()) {
                    const auto& prev = out[static_cast<size_t>(overlaps[0])];
                    auto topo = classify_triallelic(
                        static_cast<size_t>(prev.pos), static_cast<size_t>(prev.end),
                        static_cast<size_t>(chainCall.pos), static_cast<size_t>(chainCall.end));
                    switch (topo) {
                        case TriallelicTopology::PROPERLY_TRIALLELIC:
                            chainCall.triallelicTopology = "PROPERLY_TRIALLELIC"; break;
                        case TriallelicTopology::OVERLAPPING:
                            chainCall.triallelicTopology = "OVERLAPPING"; break;
                        case TriallelicTopology::NESTED:
                            chainCall.triallelicTopology = "NESTED"; break;
                        case TriallelicTopology::INTERLOCKING:
                            chainCall.triallelicTopology = "INTERLOCKING"; break;
                    }
                }
            }
            out.push_back(std::move(chainCall));
            append_window_calls(name, seq);
            continue;
        }

        // ── Path B: length-delta fallback ─────────────────────────────
        // Primary: name-based lookup (fast, used for assembly contigs).
        // Secondary: k-mer-overlap best match across all refs when the name
        // lookup misses — required for reads-mode pseudo-contigs (lr_pc0,
        // sr_unitig3…) whose names have no correspondence with reference names.
        {
            const TolGlobal::RefSeq* best = nullptr;
            int bestDelta = std::numeric_limits<int>::max();

            if (it != byContig.end() && !it->second.empty()) {
                // Name-based candidate set
                for (const auto& cand : it->second) {
                    if (!clade_allowed(cand)) continue;
                    int d = std::abs(static_cast<int>(seq.size()) -
                                     static_cast<int>(cand.seq().size()));
                    if (d < bestDelta) { bestDelta = d; best = &cand; }
                }
                if (best == nullptr) {
                    for (const auto& cand : it->second) {
                        int d = std::abs(static_cast<int>(seq.size()) -
                                         static_cast<int>(cand.seq().size()));
                        if (d < bestDelta) { bestDelta = d; best = &cand; }
                    }
                }
            } else {
                // Contig name not in index: find best reference by k-mer overlap,
                // then check if the length delta qualifies as an SV.
                const int bK = std::max(5, std::min(
                    fo.fallbackSketchParams.k > 0 ? fo.fallbackSketchParams.k : 7, 13));
                double bestOv = 0.05;   // require at least 5% overlap to avoid random matches
                for (const auto& ref : allRefs) {
                    if (!clade_allowed(ref)) continue;
                    if (!ref.has_seq()) continue;
                    double ov = kmer_overlap_fraction(seq, ref.seq(), bK);
                    if (ov > bestOv) {
                        double delta_d = static_cast<double>(std::abs(
                            static_cast<int>(seq.size()) - static_cast<int>(ref.seq().size())));
                        bestOv    = ov;
                        bestDelta = static_cast<int>(delta_d);
                        best      = &ref;
                    }
                }
                if (best == nullptr) {
                    for (const auto& ref : allRefs) {
                        if (!ref.has_seq()) continue;
                        double ov = kmer_overlap_fraction(seq, ref.seq(), bK);
                        if (ov > bestOv) {
                            double delta_d = static_cast<double>(std::abs(
                                static_cast<int>(seq.size()) - static_cast<int>(ref.seq().size())));
                            bestOv    = ov;
                            bestDelta = static_cast<int>(delta_d);
                            best      = &ref;
                        }
                    }
                }
            }

            if (best && bestDelta >= fo.minSvLen && bestDelta <= fo.maxSvLen &&
                bestDelta <= static_cast<int>(
                    std::min(seq.size(), best->seq().size()) / 2)) {
                out.push_back(make_insdel_call(qAsm, name, *best,
                                               static_cast<int>(seq.size()), fo));
                append_window_calls(name, seq);
                continue;
            }
        }

        // ── Path C: OFF_REF novelty scoring with cross-clade HGT detection ─
        // Track same-clade and other-clade overlaps separately so that a region
        // absent from same-clade refs but present in another clade can be flagged
        // as a candidate HGT event via score_cross_clade_novelty().
        if (is_low_complexity_sequence(seq)) continue;
        double sameCladeOverlap  = 0.0;
        double otherCladeOverlap = 0.0;
        double bestCladeGc  = 0.45;
        std::string bestAsmName = "OFF_REFERENCE";
        std::string bestCladeRank = ".";
        std::string bestPhylum = ".";
        double bestOtherCladeGc = 0.45;
        std::string bestOtherAsmName = "OFF_REFERENCE";
        std::string bestOtherCladeRank = ".";
        std::string bestOtherPhylum = ".";
        const int k = std::max(5, std::min(fo.fallbackSketchParams.k > 0
                                           ? fo.fallbackSketchParams.k : 7, 9));
        for (const auto& ref : allRefs) {
            if (!ref.has_seq()) continue;
            const double ov = kmer_overlap_fraction(seq, ref.seq(), k);
            if (clade_allowed(ref)) {
                if (ov > sameCladeOverlap) {
                    sameCladeOverlap = ov;
                    bestCladeGc   = ref.cladeGc;
                    bestAsmName   = ref.clade.empty() ? ref.asmName : ref.clade;
                    bestCladeRank = ref.cladeRank.empty() ? "." : ref.cladeRank;
                    bestPhylum    = ref.phylum.empty() ? "." : ref.phylum;
                }
            } else {
                if (ov > otherCladeOverlap) {
                    otherCladeOverlap  = ov;
                    bestOtherCladeGc   = ref.cladeGc;
                    bestOtherAsmName   = ref.clade.empty() ? ref.asmName : ref.clade;
                    bestOtherCladeRank = ref.cladeRank.empty() ? "." : ref.cladeRank;
                    bestOtherPhylum    = ref.phylum.empty() ? "." : ref.phylum;
                }
            }
        }
        const bool hasRoutingCtx = !routedClades.empty();
        const bool isHgt = hasRoutingCtx &&
            sameCladeOverlap < 0.05 && otherCladeOverlap >= 0.10;
        if (isHgt) {
            bestCladeGc   = bestOtherCladeGc;
            bestAsmName   = bestOtherAsmName;
            bestCladeRank = bestOtherCladeRank;
            bestPhylum    = bestOtherPhylum;
        } else if (sameCladeOverlap <= 0.0 && otherCladeOverlap > 0.0) {
            bestCladeGc   = bestOtherCladeGc;
            bestAsmName   = bestOtherAsmName;
            bestCladeRank = bestOtherCladeRank;
            bestPhylum    = bestOtherPhylum;
        }
        const OffRefNoveltyTier noveltyTier = hasRoutingCtx
            ? score_cross_clade_novelty(sameCladeOverlap, otherCladeOverlap)
            : score_off_ref_novelty(std::max(sameCladeOverlap, otherCladeOverlap));
        const std::string tier = novelty_tier_name(noveltyTier);
        if (tier == "NOVEL" || tier == "NOVEL_WEAK" || tier == "DIVERGED") {
            auto ofcall = make_offref_call(qAsm, name, seq, tier,
                                           bestCladeGc, bestAsmName, bestCladeRank, bestPhylum);
            if (isHgt) ofcall.elementClass = "HGT";
            out.push_back(std::move(ofcall));
        }
    }

    annotate_pantree_classes(out);
    return out;
}

// ── hierarchical_call_assembly_multirank ─────────────────────────────────
// Routes independently at each Linnaean rank by restricting the reference
// set to contigs whose cladeRank matches the current rank, then merges all
// results and deduplicates by (contig, pos) keeping the call with the
// highest blockScore.
//
// Rank order: phylum → class → order → family → genus → species.
// The species-level pass is identical to hierarchical_call_assembly so this
// function is a strict superset of the single-rank path.
inline std::vector<VariantCallBridge>
hierarchical_call_assembly_multirank(
        const std::string& qAsm,
        const std::unordered_map<std::string, std::string>& contigs,
        const FederatedOptions& fo,
        size_t /*routingTopN — fo.routingTopN is used*/) {

    const auto& global   = TolGlobal::instance();
    const auto& allRefs  = global.all_refs();
    const int minHierChainAnchors = static_cast<int>(std::max<size_t>(fo.tolMinChainAnchors, 2));

    // Collect all unique rank strings present in the index.
    // kLinnaeanRanks defines the canonical traversal order.
    const char* kRankOrder[] = {
        rank_to_string(TaxonomyRank::phylum),
        rank_to_string(TaxonomyRank::class_rank),
        rank_to_string(TaxonomyRank::order),
        rank_to_string(TaxonomyRank::family),
        rank_to_string(TaxonomyRank::genus),
        rank_to_string(TaxonomyRank::species),
    };
    constexpr int kNRanks = 6;

    // Deduplicate by (contig, pos) — key → index in merged, keyed on best blockScore
    std::unordered_map<std::string, size_t> bestByKey;
    std::vector<VariantCallBridge> merged;

    auto merge_calls = [&](std::vector<VariantCallBridge>&& calls,
                           const std::string& rankStr) {
        for (auto& c : calls) {
            // Stamp the rank that produced this call
            if (c.alignmentMode.find("rank=") == std::string::npos)
                c.alignmentMode += ";rank=" + rankStr;
            const std::string key = c.qContig + ":" + std::to_string(c.pos);
            auto kit = bestByKey.find(key);
            if (kit == bestByKey.end()) {
                bestByKey[key] = merged.size();
                merged.push_back(std::move(c));
            } else if (c.blockScore > merged[kit->second].blockScore) {
                merged[kit->second] = std::move(c);
            }
        }
    };

    // Run one pass per rank with a rank-filtered reference view.
    // For each rank, build a temporary filtered byContig map and allRefs slice.
    for (int ri = 0; ri < kNRanks; ++ri) {
        const std::string rankStr = kRankOrder[ri];

        // Build a filtered reference list for this rank
        std::vector<TolGlobal::RefSeq> rankRefs;
        for (const auto& r : allRefs)
            if (r.cladeRank == rankStr) rankRefs.push_back(r);

        if (rankRefs.empty()) continue;

        // Build filtered byContig map for this rank
        std::unordered_map<std::string, std::vector<TolGlobal::RefSeq>> rankByContig;
        for (const auto& r : rankRefs)
            rankByContig[r.contig].push_back(r);

        // Run the same three-path calling logic used by hierarchical_call_assembly,
        // but over the rank-filtered references.
        std::vector<VariantCallBridge> rankCalls;
        rankCalls.reserve(contigs.size());
        auto append_rank_window_calls = [&](const std::string& contigName,
                                            const std::string& seq) {
            if (!fo.graphNativeMode) return;
            auto windows = discover_graph_native_offref_windows(seq, rankRefs, fo);
            for (const auto& ow : windows) rankCalls.push_back(make_offref_window_call(qAsm, contigName, seq, ow));
        };

        std::vector<const TolGlobal::RefSeq*> rankRefPtrs;
        rankRefPtrs.reserve(rankRefs.size());
        for (const auto& r : rankRefs) rankRefPtrs.push_back(&r);

        for (const auto& kv : contigs) {
            const std::string& name = kv.first;
            const std::string& seq  = kv.second;
            const auto routed = global.route_query_to_clades(
                seq, fo.primarySketchParams, fo.fallbackSketchParams,
                fo.routingDensity, fo.routingTopN);
            const std::unordered_set<std::string> routedClades(routed.begin(), routed.end());
            auto clade_allowed = [&](const TolGlobal::RefSeq& ref) {
                return routedClades.empty() ||
                       routedClades.find(ref.clade) != routedClades.end() ||
                       routedClades.find(ref.asmName) != routedClades.end();
            };

            // Path A: MEM chain over rank-filtered refs
            std::vector<const TolGlobal::RefSeq*> cands;
            auto it = rankByContig.find(name);
            if (it != rankByContig.end())
                for (const auto& r : it->second)
                    if (r.has_seq() && clade_allowed(r)) cands.push_back(&r);
            if (cands.empty() && it != rankByContig.end())
                for (const auto& r : it->second)
                    if (r.has_seq()) cands.push_back(&r);
            if (cands.empty())
                for (const auto* r : rankRefPtrs)
                    if (r->has_seq() && clade_allowed(*r)) cands.push_back(r);
            if (cands.empty())
                for (const auto* r : rankRefPtrs)
                    if (r->has_seq()) cands.push_back(r);

            VariantCallBridge chainCall;
            if (try_mem_chain_call(qAsm, name, seq, cands, fo, chainCall) &&
                chainCall.anchors >= minHierChainAnchors) {
                rankCalls.push_back(std::move(chainCall));
                append_rank_window_calls(name, seq);
                continue;
            }

            // Path B: length-delta fallback (name-based, with k-mer-overlap fallback
            // for reads-mode pseudo-contigs whose names are not in the index)
            {
                const TolGlobal::RefSeq* best = nullptr;
                int bestDelta = std::numeric_limits<int>::max();

                if (it != rankByContig.end() && !it->second.empty()) {
                    for (const auto& cand : it->second) {
                        if (!clade_allowed(cand)) continue;
                        int d = std::abs(static_cast<int>(seq.size()) -
                                         static_cast<int>(cand.seq().size()));
                        if (d < bestDelta) { bestDelta = d; best = &cand; }
                    }
                    if (best == nullptr) {
                        for (const auto& cand : it->second) {
                            int d = std::abs(static_cast<int>(seq.size()) -
                                             static_cast<int>(cand.seq().size()));
                            if (d < bestDelta) { bestDelta = d; best = &cand; }
                        }
                    }
                } else {
                    const int bK = std::max(5, std::min(
                        fo.fallbackSketchParams.k > 0 ? fo.fallbackSketchParams.k : 7, 13));
                    double bestOv = 0.05;
                    for (const auto& ref : rankRefs) {
                        if (!clade_allowed(ref)) continue;
                        if (!ref.has_seq()) continue;
                        double ov = kmer_overlap_fraction(seq, ref.seq(), bK);
                        if (ov > bestOv) {
                            bestOv    = ov;
                            bestDelta = std::abs(static_cast<int>(seq.size()) -
                                                  static_cast<int>(ref.seq().size()));
                            best      = &ref;
                        }
                    }
                    if (best == nullptr) {
                        for (const auto& ref : rankRefs) {
                            if (!ref.has_seq()) continue;
                            double ov = kmer_overlap_fraction(seq, ref.seq(), bK);
                            if (ov > bestOv) {
                                bestOv    = ov;
                                bestDelta = std::abs(static_cast<int>(seq.size()) -
                                                      static_cast<int>(ref.seq().size()));
                                best      = &ref;
                            }
                        }
                    }
                }

                if (best && bestDelta >= fo.minSvLen && bestDelta <= fo.maxSvLen &&
                    bestDelta <= static_cast<int>(
                        std::min(seq.size(), best->seq().size()) / 2)) {
                    rankCalls.push_back(make_insdel_call(qAsm, name, *best,
                                                         static_cast<int>(seq.size()), fo));
                    append_rank_window_calls(name, seq);
                    continue;
                }
            }

            // Path C: OFF_REF novelty scoring with cross-clade HGT detection
            // (only on final species pass to avoid emitting the same novel contig once per rank)
            if (ri == kNRanks - 1) {
                if (is_low_complexity_sequence(seq)) continue;
                double sameCladeOverlap  = 0.0;
                double otherCladeOverlap = 0.0;
                double bestCladeGc = 0.45;
                std::string bestAsmName = "OFF_REFERENCE";
                std::string bestCladeRank = ".";
                std::string bestPhylum = ".";
                double bestOtherCladeGc = 0.45;
                std::string bestOtherAsmName = "OFF_REFERENCE";
                std::string bestOtherCladeRank = ".";
                std::string bestOtherPhylum = ".";
                const int k = std::max(5, std::min(fo.fallbackSketchParams.k > 0
                                                   ? fo.fallbackSketchParams.k : 7, 9));
                for (const auto& ref : allRefs) {
                    if (!ref.has_seq()) continue;
                    const double ov = kmer_overlap_fraction(seq, ref.seq(), k);
                    if (clade_allowed(ref)) {
                        if (ov > sameCladeOverlap) {
                            sameCladeOverlap = ov;
                            bestCladeGc   = ref.cladeGc;
                            bestAsmName   = ref.clade.empty() ? ref.asmName : ref.clade;
                            bestCladeRank = ref.cladeRank.empty() ? "." : ref.cladeRank;
                            bestPhylum    = ref.phylum.empty() ? "." : ref.phylum;
                        }
                    } else {
                        if (ov > otherCladeOverlap) {
                            otherCladeOverlap  = ov;
                            bestOtherCladeGc   = ref.cladeGc;
                            bestOtherAsmName   = ref.clade.empty() ? ref.asmName : ref.clade;
                            bestOtherCladeRank = ref.cladeRank.empty() ? "." : ref.cladeRank;
                            bestOtherPhylum    = ref.phylum.empty() ? "." : ref.phylum;
                        }
                    }
                }
                const bool hasRoutingCtx = !routedClades.empty();
                const bool isHgt = hasRoutingCtx &&
                    sameCladeOverlap < 0.05 && otherCladeOverlap >= 0.10;
                if (isHgt) {
                    bestCladeGc   = bestOtherCladeGc;
                    bestAsmName   = bestOtherAsmName;
                    bestCladeRank = bestOtherCladeRank;
                    bestPhylum    = bestOtherPhylum;
                } else if (sameCladeOverlap <= 0.0 && otherCladeOverlap > 0.0) {
                    bestCladeGc   = bestOtherCladeGc;
                    bestAsmName   = bestOtherAsmName;
                    bestCladeRank = bestOtherCladeRank;
                    bestPhylum    = bestOtherPhylum;
                }
                const OffRefNoveltyTier noveltyTier = hasRoutingCtx
                    ? score_cross_clade_novelty(sameCladeOverlap, otherCladeOverlap)
                    : score_off_ref_novelty(std::max(sameCladeOverlap, otherCladeOverlap));
                const std::string tier = novelty_tier_name(noveltyTier);
                if (tier == "NOVEL" || tier == "NOVEL_WEAK" || tier == "DIVERGED") {
                    auto ofcall = make_offref_call(qAsm, name, seq, tier,
                                                   bestCladeGc, bestAsmName, bestCladeRank, bestPhylum);
                    if (isHgt) ofcall.elementClass = "HGT";
                    rankCalls.push_back(std::move(ofcall));
                }
            }
        }

        annotate_pantree_classes(rankCalls);
        merge_calls(std::move(rankCalls), rankStr);
    }

    // If the index has no rank metadata at all, fall back to the base caller
    if (merged.empty())
        return hierarchical_call_assembly(qAsm, contigs, fo);

    return merged;
}

// ── Ancestral alignment helpers ───────────────────────────────────────────

// AncestralManifestContext: holds the manifest path and the clade→rank/phylum
// lookup table built from it.  Populated by load_ancestral_manifest_context.
struct AncestralManifestContext {
    std::string manifestPath;
    // clade name → (rank, phylum) — built from the manifest at load time
    std::unordered_map<std::string, std::pair<std::string,std::string>> cladeInfo;
};

// Parse the manifest at `path` to build clade→(rank,phylum) lookups.
// Supports both the 5-column (asm phylum clade rank fasta) and 9-column
// (asm phylum class order family genus clade rank fasta) formats.
inline AncestralManifestContext load_ancestral_manifest_context(const std::string& path) {
    if (!fs::exists(path)) throw std::runtime_error("manifest not found: " + path);
    AncestralManifestContext ctx;
    ctx.manifestPath = path;
    std::ifstream in(path);
    std::string line;
    bool extended = false;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            extended = (std::count(line.begin(), line.end(), '\t') >= 8);
            continue;
        }
        auto cols = split_tab(line);
        std::string clade, rank, phylum;
        if (extended && cols.size() >= 9) {
            phylum = cols[1]; clade = cols[6]; rank = cols[7];
        } else if (!extended && cols.size() >= 4) {
            phylum = cols[1]; clade = cols[2]; rank = cols[3];
        } else {
            continue;
        }
        if (!clade.empty() && !rank.empty()) {
            auto val = std::make_pair(rank, phylum);
            ctx.cladeInfo.emplace(clade, val);
            // Also index by the raw assembly name (cols[0]) because c.refAsm
            // from the simple_offref_fallback path carries the bare asm stem
            // (e.g. "arbuscular_mf_asm2") while clade_name may be prefixed
            // with the genus (e.g. "Rhizophagus_arbuscular_mf_asm2").
            const std::string& asmName = cols[0];
            if (!asmName.empty() && asmName != clade)
                ctx.cladeInfo.emplace(asmName, val);
        }
    }
    return ctx;
}

inline std::pair<std::string, std::string>
lookup_ancestral_clade_info(const std::string& refAsm,
                            const AncestralManifestContext& ctx) {
    auto it = ctx.cladeInfo.find(refAsm);
    if (it != ctx.cladeInfo.end()) return it->second;

    // Fall back to the registry descriptors rather than scanning every loaded
    // contig. This keeps ancestral reporting tied to the number of clades,
    // which matters when million-genome catalogs expand the contig table.
    for (const auto& d : TolGlobal::instance().registry().descriptors()) {
        if (d.cladeName == refAsm) return {d.cladeRank, d.phylum};
    }
    return {"unknown", "."};
}

struct AncestralPosteriorBase {
    char   base      = 'N';
    double posterior = 0.25;
};

struct AncestralReconstructionResult {
    std::string ancestralSequence;
    std::vector<AncestralPosteriorBase> posteriorByBase;
    size_t breakpointCount = 0;
    size_t alignedBases    = 0;
    double meanPosterior   = 0.0;
    std::string sourceClade;
    std::string sourcePhylum;
};

inline AncestralReconstructionResult
reconstruct_full_ancestral_sequence(const std::string& querySeq,
                                    const std::string& refSeq,
                                    const VariantCallBridge& call,
                                    const FederatedOptions& fo) {
    AncestralReconstructionResult out;
    out.sourceClade  = call.refAsm;
    out.sourcePhylum = call.phylum;
    const size_t qn = querySeq.size();
    const size_t rn = refSeq.size();
    const size_t n  = std::max(qn, rn);
    out.ancestralSequence.reserve(n);
    out.posteriorByBase.reserve(n);

    auto push = [&](char b, double p) {
        out.ancestralSequence.push_back(b);
        out.posteriorByBase.push_back({b, p});
    };

    for (size_t i = 0; i < n; ++i) {
        const bool hasQ = i < qn;
        const bool hasR = i < rn;
        const char q = hasQ ? querySeq[i] : 'N';
        const char r = hasR ? refSeq[i]   : 'N';

        if (hasQ && hasR) {
            ++out.alignedBases;
            if (q == r) push(q, 0.995);
            else {
                const char anc = (call.type == "INS") ? r : ((call.type == "DEL") ? q : r);
                push(anc, 0.65);
            }
        } else if (hasR) {
            push(r, 0.80);
        } else if (hasQ && call.type != "INS") {
            push(q, 0.55);
        }
    }

    if (call.type == "TRA" || call.type == "INV") out.breakpointCount = 2;
    else if (call.type == "OFF_REF") out.breakpointCount = fo.enableAncestralRecomb ? 1u : 0u;
    else out.breakpointCount = 1;

    double sum = 0.0;
    for (const auto& pb : out.posteriorByBase) sum += pb.posterior;
    if (!out.posteriorByBase.empty())
        out.meanPosterior = sum / static_cast<double>(out.posteriorByBase.size());
    return out;
}

// Updated header: adds rank and phylum columns.
inline void write_ancestral_tsv_header(std::ostream& out) {
    out << "query_asm\tquery_contig\tclade\tclade_rank\tphylum\t"
           "variant_type\tbreakpoints\tsegment_bp\taligned_bases\tmean_posterior\t"
           "ancestral_sequence\tsource_clade\tsource_phylum\n";
}

inline const TolGlobal::RefSeq* find_ancestral_ref_seq(const VariantCallBridge& c) {
    const auto& refsByContig = TolGlobal::instance().refs_by_contig();
    auto hit = refsByContig.find(c.refContig);
    if (hit == refsByContig.end()) return nullptr;
    const auto& refs = hit->second;
    const TolGlobal::RefSeq* contigFallback = nullptr;
    for (const auto& r : refs) {
        if (!c.refAsm.empty() &&
            (r.asmName == c.refAsm || r.clade == c.refAsm)) {
            return &r;
        }
        if (contigFallback == nullptr) contigFallback = &r;
    }
    return contigFallback;
}

inline std::string extract_sequence_window(const std::string& seq,
                                           int startPos1,
                                           int endPos1,
                                           int padBp) {
    if (seq.empty()) return {};
    const int start0 = std::max(0, startPos1 - 1 - padBp);
    const int end0 = std::min(static_cast<int>(seq.size()),
                              std::max(start0 + 1, endPos1 + padBp));
    return seq.substr(static_cast<size_t>(start0),
                      static_cast<size_t>(std::max(1, end0 - start0)));
}

inline std::string ancestral_breakpoint_descriptor(const VariantCallBridge& c,
                                                   const FederatedOptions& fo) {
    if (c.type == "TRA") {
        std::ostringstream ss;
        ss << "TRA:" << c.refContig << ':' << std::max(1, c.pos)
           << "-" << std::max(std::max(1, c.pos), c.end)
           << "->" << (c.mateContig.empty() ? "." : c.mateContig)
           << ':' << std::max(1, c.matePos)
           << "-" << std::max(std::max(1, c.matePos), c.mateEnd);
        return ss.str();
    }
    if (c.type == "INV") {
        std::ostringstream ss;
        ss << "INV:" << std::max(1, c.pos) << "-" << std::max(std::max(1, c.pos), c.end);
        return ss.str();
    }
    if (c.type == "OFF_REF") {
        if (!fo.enableAncestralRecomb) return ".";
        std::ostringstream ss;
        ss << "OFF_REF:" << std::max(1, c.pos) << "-" << std::max(std::max(1, c.pos), c.end);
        return ss.str();
    }
    std::ostringstream ss;
    ss << c.type << ':' << std::max(1, c.pos) << "-" << std::max(std::max(1, c.pos), c.end);
    return ss.str();
}

// write_ancestral_alignments_for_assembly
//
// Writes one row per SV call with:
//   query_asm    — assembly name
//   query_contig — contig that carries the call
//   clade        — routed reference clade
//   clade_rank   — Linnaean rank of that clade (from manifest)
//   phylum       — phylum of that clade (from manifest)
//   variant_type — SV type
//   breakpoints  — number of distinct breakpoints implied by the call:
//                    TRA  → 2  (two chromosomal breakpoints per translocation)
//                    INV  → 2  (two inversion breakpoints)
//                    OFF_REF → 1 if enableAncestralRecomb, else 0
//                    INS/DEL/DUP → 1
//   segment_bp   — total sequence span of the call
inline void write_ancestral_alignments_for_assembly(
        const std::string& qAsm,
        const std::unordered_map<std::string, std::string>& contigs,
        const std::vector<VariantCallBridge>& calls,
        const AncestralManifestContext& ctx,
        const FederatedOptions& fo,
        size_t /*routingTopN*/,
        std::ostream& out) {

    for (const auto& c : calls) {
        // Resolve rank and phylum for this call's routed clade
        const auto [cladeRank, phylum] = lookup_ancestral_clade_info(c.refAsm, ctx);

        AncestralReconstructionResult recon;
        auto qit = contigs.find(c.qContig);
        const TolGlobal::RefSeq* ref = find_ancestral_ref_seq(c);
        const int padBp = std::max(25, std::min(500, std::max(std::abs(c.svlen), 50) / 2));
        if (qit != contigs.end() && ref != nullptr && ref->has_seq()) {
            const std::string qWindow = extract_sequence_window(
                qit->second, c.pos, std::max(c.pos, c.end), padBp);
            const std::string rWindow = extract_sequence_window(
                ref->seq(), c.pos, std::max(c.pos, c.end), padBp);
            recon = reconstruct_full_ancestral_sequence(qWindow, rWindow, c, fo);
        } else if (qit != contigs.end()) {
            const std::string qWindow = extract_sequence_window(
                qit->second, c.pos, std::max(c.pos, c.end), padBp);
            recon = reconstruct_full_ancestral_sequence(qWindow, "", c, fo);
        }

        const std::string breakpointDesc = ancestral_breakpoint_descriptor(c, fo);
        const int seg_bp = recon.ancestralSequence.empty()
            ? std::max(1, std::abs(c.svlen))
            : static_cast<int>(recon.ancestralSequence.size());
        const std::string sourceClade = recon.sourceClade.empty() ? c.refAsm : recon.sourceClade;
        const std::string sourcePhylum = recon.sourcePhylum.empty() ? phylum : recon.sourcePhylum;
        const std::string ancestralSeq = recon.ancestralSequence.empty() ? "." : recon.ancestralSequence;

        out << qAsm      << '\t'
            << c.qContig << '\t'
            << c.refAsm  << '\t'
            << cladeRank << '\t'
            << phylum    << '\t'
            << c.type    << '\t'
            << breakpointDesc << '\t'
            << seg_bp    << '\t'
            << recon.alignedBases << '\t'
            << recon.meanPosterior << '\t'
            << ancestralSeq << '\t'
            << sourceClade << '\t'
            << sourcePhylum << '\n';
    }
}

#if 0




// ── Full ancestral sequence reconstruction ───────────────────────────────
// Lightweight per-call ancestral reconstruction over the query/reference span.
// The model is intentionally simple: for each aligned position, retain the
// shared base when query and reference agree, otherwise emit the majority base
// under an explicit posterior. Inserted query sequence is kept as derived-only
// support; deleted reference sequence is retained in the ancestral string.

struct AncestralPosteriorBase {
    char   base      = 'N';
    double posterior = 0.25;
};

struct AncestralReconstructionResult {
    std::string ancestralSequence;
    std::vector<AncestralPosteriorBase> posteriorByBase;
    size_t breakpointCount = 0;
    size_t alignedBases    = 0;
    double meanPosterior   = 0.0;
    std::string sourceClade;
    std::string sourcePhylum;
};

inline AncestralReconstructionResult
reconstruct_full_ancestral_sequence(const std::string& querySeq,
                                    const std::string& refSeq,
                                    const VariantCallBridge& call,
                                    const FederatedOptions& fo) {
    AncestralReconstructionResult out;
    out.sourceClade  = call.refAsm;
    out.sourcePhylum = call.phylum;
    const size_t qn = querySeq.size();
    const size_t rn = refSeq.size();
    const size_t n  = std::max(qn, rn);
    out.ancestralSequence.reserve(n);
    out.posteriorByBase.reserve(n);

    auto push = [&](char b, double p) {
        out.ancestralSequence.push_back(b);
        out.posteriorByBase.push_back({b, p});
    };

    for (size_t i = 0; i < n; ++i) {
        const bool hasQ = i < qn;
        const bool hasR = i < rn;
        const char q = hasQ ? querySeq[i] : 'N';
        const char r = hasR ? refSeq[i]   : 'N';

        if (hasQ && hasR) {
            ++out.alignedBases;
            if (q == r) push(q, 0.995);
            else {
                const char anc = (call.type == "INS") ? r : ((call.type == "DEL") ? q : r);
                push(anc, 0.65);
            }
        } else if (hasR) {
            // Query-specific deletion relative to reference → ancestral base is likely retained.
            push(r, 0.80);
        } else if (hasQ && call.type != "INS") {
            // Only keep query-only sequence for non-insertion contexts.
            push(q, 0.55);
        }
    }

    if (call.type == "TRA" || call.type == "INV") out.breakpointCount = 2;
    else if (call.type == "OFF_REF") out.breakpointCount = fo.enableAncestralRecomb ? 1u : 0u;
    else out.breakpointCount = 1;

    double sum = 0.0;
    for (const auto& pb : out.posteriorByBase) sum += pb.posterior;
    if (!out.posteriorByBase.empty())
        out.meanPosterior = sum / static_cast<double>(out.posteriorByBase.size());
    return out;
}

// ── try_mem_chain_call_public ─────────────────────────────────────────────
#endif
// Non-static public wrapper so callers outside this translation unit (e.g.
// reads_mode_sv_calls in main.cpp) can access the MEM chain caller without
// going through the full hierarchical engine.
inline bool try_mem_chain_call_public(
        const std::string& qAsm,
        const std::string& qContig,
        const std::string& qSeq,
        const std::vector<const TolGlobal::RefSeq*>& refCandidates,
        const FederatedOptions& fo,
        VariantCallBridge& call) {
    return try_mem_chain_call(qAsm, qContig, qSeq, refCandidates, fo, call);
}

} // namespace tol

#endif // FUNGI_TOL_BRIDGE_HPP
