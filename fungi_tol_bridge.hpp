#ifndef FUNGI_TOL_BRIDGE_HPP
#define FUNGI_TOL_BRIDGE_HPP

// fungi_tol_bridge.hpp - v15
//
// Integration layer for hierarchical fungal SV calling. Key responsibilities:
//  - annotate DEL and OFF_REF sequence with repeat/TE/HGT element classes;
//  - score cross-clade novelty using separate same-clade and other-clade overlap;
//  - build and query multi-rank indices across Linnaean ranks;
//  - recover MEM-chain provenance without quadratic chain-to-MEM lookup;
//  - emit GFA records with element-class tags for downstream reports.

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <deque>
#include <ext/stdio_filebuf.h>
#include <filesystem>
#include <fstream>
#include <functional>
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
#include <string_view>
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

// gz-transparent FASTA stream
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

// VariantCallBridge
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
    // Query input mode - provenance label written to TSV/VCF
    std::string queryMode           = "assembly";
    // Probabilistic multi-evidence fusion summary for downstream ranking/reporting.
    double      fusedPosteriorAlt   = 0.50;
    double      fusedLogOddsAlt     = 0.0;
    double      fusedEffectiveDepth = 0.0;
    int         fusedLayersUsed     = 0;
    // Evidence support backing this call:
    //   assembly    -> anchor count / fused evidence layers, filled before output
    //   long-reads  -> cluster size (_n<N> in contig name), exact read count
    //   short-reads -> min k-mer frequency along unitig path (_mf<N>), coverage proxy
    int         readSupport         = -1;
    // Sequence carried by the variant allele/affected segment when compact
    // enough for VCF INFO output. For INS/DUP/INV/OFF_REF this is query
    // sequence; for DEL it is the deleted reference sequence.
    std::string variantSeq;
};

inline void set_variant_sequence_excerpt(VariantCallBridge& v,
                                         std::string_view seq,
                                         size_t maxLen = 5000) {
    if (seq.empty()) return;
    const size_t n = std::min(seq.size(), maxLen);
    v.variantSeq.assign(seq.substr(0, n));
    for (char& ch : v.variantSeq) {
        ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        if (ch != 'A' && ch != 'C' && ch != 'G' && ch != 'T' && ch != 'N')
            ch = 'N';
    }
}

namespace tol {

using ::VariantCallBridge;

// FederatedOptions
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

// CladeGraphDescriptor (bridge-side)
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

// sanitize_name
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

// split_tab / split_csv
// Registry manifests written with CRLF line endings put a stray '\r' on the
// last cell of every line. Strip it so downstream path lookups see clean cells.
inline std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    std::istringstream ss(line);
    while (std::getline(ss, cur, '\t')) {
        if (!cur.empty() && cur.back() == '\r') cur.pop_back();
        out.push_back(std::move(cur));
    }
    return out;
}

inline std::vector<std::string> split_csv(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    std::istringstream ss(s);
    while (std::getline(ss, cur, ',')) {
        if (!cur.empty() && cur.back() == '\r') cur.pop_back();
        if (!cur.empty()) out.push_back(std::move(cur));
    }
    return out;
}


// read_fasta_local
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

// for_each_fasta_record
// Streaming variant of read_fasta_local: hands each (name, seq) pair to
// `fn` one at a time, reusing the same buffers between records. Callers
// that do not need every contig held in memory at once (e.g. the full-base
// graph builder, which slices each contig into segments and never reuses
// it) should prefer this over read_fasta_local to halve transient RSS.
// Returns true if at least one record was emitted.
template <class Fn>
inline bool for_each_fasta_record(const std::string& path, Fn&& fn) {
    FastaStream fs(path);
    std::istream& in = fs.get();
    std::string name, seq, line;
    bool any = false;
    auto emit = [&]() {
        if (name.empty() || seq.empty()) return;
        fn(static_cast<const std::string&>(name),
           static_cast<const std::string&>(seq));
        any = true;
    };
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '>') {
            emit();
            name.assign(line, 1, std::string::npos);
            auto sp = name.find_first_of(" \t");
            if (sp != std::string::npos) name.resize(sp);
            seq.clear();
        } else {
            seq += line;
        }
    }
    emit();
    return any;
}

// Rolling-hash k-mer overlap
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

inline std::unordered_set<uint64_t> kmer_hashes(std::string_view seq, int k) {
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

inline std::unordered_set<uint64_t> kmer_hashes(const std::string& seq, int k) {
    return kmer_hashes(std::string_view(seq.data(), seq.size()), k);
}

// k-mer Jaccard similarity: |A intersect B| / |A union B|. Standard Jaccard
// avoids containment-style overstatement when one sequence is much shorter
// than the other.
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
    // Jaccard: |A  intersect  B| / |A  union  B|
    const size_t uni = ha.size() + hb.size() - inter;
    return uni == 0 ? 0.0 : static_cast<double>(inter) / static_cast<double>(uni);
}

inline double kmer_best_containment_fraction(const std::string& a, const std::string& b, int k) {
    auto ha = kmer_hashes(a, k);
    auto hb = kmer_hashes(b, k);
    if (ha.empty() || hb.empty()) return 0.0;
    const auto* small = &ha;
    const auto* big   = &hb;
    if (small->size() > big->size()) std::swap(small, big);
    size_t inter = 0;
    for (uint64_t h : *small)
        if (big->count(h)) ++inter;
    const double aContain = static_cast<double>(inter) / static_cast<double>(ha.size());
    const double bContain = static_cast<double>(inter) / static_cast<double>(hb.size());
    return std::max(aContain, bContain);
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

inline bool is_low_complexity_sequence(std::string_view seq) {
    if (seq.size() < 5u) return true;
    std::unordered_set<char> alphabet;
    for (char ch : seq)
        if (!std::isspace(static_cast<unsigned char>(ch)))
            alphabet.insert(static_cast<char>(
                std::toupper(static_cast<unsigned char>(ch))));
    if (alphabet.size() <= 1) return true;
    // Reject pure homopolymer and very-low-diversity 5-mer profiles, while
    // keeping the floor low enough that AT-rich STARSHIP hulls (typically
    // 6-30 distinct 5-mers per 500 bp window) and TE termini pass.
    auto ks = kmer_hashes(seq, 5);
    return ks.size() < 4;
}

inline bool is_low_complexity_sequence(const std::string& seq) {
    return is_low_complexity_sequence(std::string_view(seq.data(), seq.size()));
}

inline std::string infer_novelty_tier(double overlapFraction) {
    return novelty_tier_name(score_off_ref_novelty(overlapFraction));
}

// ManifestRegistry
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
        std::unordered_map<std::string, size_t> col;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            if (line[0] == '#') {
                auto header = split_tab(line.substr(1));
                for (size_t i = 0; i < header.size(); ++i) col[header[i]] = i;
                continue;
            }
            auto cols = split_tab(line);
            if (cols.size() < 7) continue;
            CladeGraphDescriptor d;
            auto get = [&](const std::string& name, size_t fallback) -> std::string {
                auto it = col.find(name);
                const size_t idx = (it == col.end()) ? fallback : it->second;
                return idx < cols.size() ? cols[idx] : std::string();
            };
            d.cladeName       = get("clade_name", 0);
            d.cladeRank       = get("clade_rank", 1);
            d.phylum          = get("phylum", 2);
            d.graphPath       = get("graph_path", 3);
            d.genomeCount     = static_cast<size_t>(std::stoull(get("genome_count", 4)));
            d.svBubbles       = static_cast<size_t>(std::stoull(get("sv_bubbles", 5)));
            d.compressedBytes = static_cast<size_t>(std::stoull(get("compressed_bytes", 6)));
            const std::string fastaCol = get("fasta_paths", SIZE_MAX);
            if (!fastaCol.empty()) d.fastaPaths = split_csv(fastaCol);
            // Compatibility: manifests without a fasta_paths column store
            // crc32 in cols[7]. Treat that as metadata rather than a path.
            if (d.fastaPaths.empty() && cols.size() >= 8 && col.find("fasta_paths") == col.end()
                && cols[7].find('/') != std::string::npos) {
                d.fastaPaths = split_csv(cols[7]);
            }
            descs_.push_back(std::move(d));
        }
    }

    const std::vector<CladeGraphDescriptor>& descriptors() const { return descs_; }

private:
    std::string dir_;
    std::vector<CladeGraphDescriptor> descs_;
};

// ManifestRow / read_manifest_rows
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

inline std::vector<std::string> collect_unique_fasta_paths_from_shard(const std::string& shardPath) {
    std::vector<std::string> paths;
    std::unordered_set<std::string> seen;
    for_each_manifest_row(shardPath, [&](const ManifestRow& r) {
        if (r.fastaPath.empty()) return;
        if (seen.insert(r.fastaPath).second) paths.push_back(r.fastaPath);
    });
    return paths;
}

inline void recover_descriptor_fasta_paths_from_hierarchy_manifest(
        std::vector<CladeGraphDescriptor>& descs,
        const fs::path& hierarchyManifest) {
    if (descs.empty() || !fs::exists(hierarchyManifest)) return;
    bool needsRecovery = false;
    for (const auto& d : descs) {
        if (d.fastaPaths.empty()) {
            needsRecovery = true;
            break;
        }
    }
    if (!needsRecovery) return;

    std::unordered_map<std::string, std::vector<std::string>> byRankClade;
    std::unordered_map<std::string, std::unordered_set<std::string>> seen;
    for_each_manifest_row(hierarchyManifest.string(), [&](const ManifestRow& r) {
        for (const auto& target : build_targets_for_manifest_row(r, true)) {
            if (r.fastaPath.empty()) continue;
            const std::string key = target.cladeRank + "||" + target.cladeName;
            if (seen[key].insert(r.fastaPath).second)
                byRankClade[key].push_back(r.fastaPath);
        }
    });
    for (auto& d : descs) {
        if (!d.fastaPaths.empty()) continue;
        const std::string key = d.cladeRank + "||" + d.cladeName;
        auto it = byRankClade.find(key);
        if (it != byRankClade.end()) d.fastaPaths = it->second;
    }
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

// BuilderFastaHashCache
// Process-global, thread-safe map FASTA-path -> minimizer hashes of the
// FASTA's leading 512 bp seed. Each FASTA is opened (and gzip-decoded)
// exactly once across the whole multi-rank index build. Upper-rank shards
// (genus / family / order / class / phylum) that exceed the full-base cap
// reuse the species-level hashes instead of re-opening every constituent
// FASTA. Save: 5x I/O on average (1 read per FASTA vs ~6 across ranks),
// which is the single largest wall-clock contributor for runs with deep
// taxonomic redundancy like panel200 (24k FASTAs).
struct BuilderFastaHashCache {
    static BuilderFastaHashCache& instance() {
        static BuilderFastaHashCache inst;
        return inst;
    }

    // Returns <=`limit` minimizer hashes computed over the FASTA's first
    // 512 bp seed using FNV-1a k-mer hashing (matches the inline hashing
    // in build_compact_manifest_graph). Recomputes on every call only if
    // the per-(path, k, limit) tuple was never cached; otherwise O(1) +
    // a vector copy of <=limit uint64.
    std::vector<uint64_t> hashes_for(const std::string& path, int k, size_t limit) {
        const Key key{path, static_cast<int32_t>(std::max(3, k)),
                      static_cast<uint32_t>(std::max<size_t>(1, limit))};
        {
            std::lock_guard<std::mutex> lk(mu_);
            auto it = cache_.find(key);
            if (it != cache_.end()) return it->second;
        }
        const std::string seed = read_fasta_seed_sequence(path, 512);
        std::vector<uint64_t> h;
        if (!seed.empty() && seed != "N") {
            const int kk = std::min<int>(key.k, static_cast<int>(seed.size()));
            h = minimizer_hashes_for_sequence(seed, kk, limit);
        }
        std::lock_guard<std::mutex> lk(mu_);
        // bounded cache: 200k entries x ~120 B/entry + <=64x8B hashes ~ 100 MB max
        if (cache_.size() < 200000) cache_.emplace(key, h);
        return h;
    }

    void clear() {
        std::lock_guard<std::mutex> lk(mu_);
        cache_.clear();
    }

private:
    struct Key {
        std::string path;
        int32_t k;
        uint32_t limit;
        bool operator==(const Key& o) const {
            return k == o.k && limit == o.limit && path == o.path;
        }
    };
    struct KeyHash {
        size_t operator()(const Key& key) const noexcept {
            size_t h = std::hash<std::string>{}(key.path);
            h ^= static_cast<size_t>(key.k) * 0x9e3779b97f4a7c15ULL;
            h ^= static_cast<size_t>(key.limit) * 0xbf58476d1ce4e5b9ULL + (h << 6) + (h >> 2);
            return h;
        }
    };
    std::unordered_map<Key, std::vector<uint64_t>, KeyHash> cache_;
    std::mutex mu_;
};

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

    // Centroid roll-up: reuse minimizer hashes already computed for this
    // FASTA at any lower rank (BuilderFastaHashCache is process-global).
    // Seed sequences are only fetched for the reservoir-sampled `samples`
    // that need a node payload (<=sampleLimit per shard, typically 64).
    const int kBase = std::max(3, sp.k);
    auto& hashCache = BuilderFastaHashCache::instance();

    for_each_manifest_row(shardPath, [&](const ManifestRow& r) {
        ++seen;
        if (r.fastaPath.empty()) return;
        auto cachedHashes = hashCache.hashes_for(r.fastaPath, kBase, 64);
        if (cachedHashes.empty()) return;  // FASTA missing / empty
        ++contributing;
        for (uint64_t h : cachedHashes) centroidAcc.insert(h);

        // Materialize the seed sequence only if this row is going to be
        // retained in the reservoir sample (most rows in a large shard
        // are evicted by the reservoir without ever needing the seed).
        bool willKeep = false;
        size_t replaceIdx = sampleLimit;
        if (samples.size() < sampleLimit) {
            willKeep = true;
        } else {
            std::uniform_int_distribution<size_t> dist(0, seen - 1);
            const size_t pick = dist(rng);
            if (pick < sampleLimit) {
                willKeep = true;
                replaceIdx = pick;
            }
        }
        if (!willKeep) return;

        const std::string seed =
            cached_seed_sequence_for_build(r.fastaPath, 512, seedCache);
        if (seed.empty() || seed == "N") return;

        CompactGraphSample sample{resolve_manifest_row_asm_name(r), seed};
        if (samples.size() < sampleLimit) {
            samples.push_back(std::move(sample));
        } else {
            samples[replaceIdx] = std::move(sample);
        }
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

// binary index helpers
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
        const std::string asmLabel =
            r.asmName.empty() ? fs::path(r.fastaPath).stem().string() : r.asmName;
        bool any = false;
        try {
            // Stream contigs one at a time instead of materializing the
            // whole genome as an unordered_map<string,string>. For a 30 MB
            // fungal genome this cuts transient RSS by ~30 MB per worker
            // and avoids the unordered_map's ~1.5x overhead - important
            // when up to 32 genomes pass through per full-base shard and
            // up to 4 such shards run concurrently (tiered worker pool).
            any = for_each_fasta_record(
                r.fastaPath,
                [&](const std::string& contigName, const std::string& seq) {
                    if (seq.empty()) return;
                    builder.add_genome(asmLabel, contigName, seq,
                                       "NONE", false, false, 10000, 1000);
                });
        } catch (...) {
            // Keep index builds resilient to individual malformed FASTA rows.
        }
        if (any) ++addedGenomes;
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

// write_graph_payload / write_routing_payload
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
    built.desc.fastaPaths = collect_unique_fasta_paths_from_shard(shard.shardPath);
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

// rank ordering helper
// Sort shards so lower ranks (species first) build before upper ranks
// (phylum last). This populates BuilderFastaHashCache during species/genus
// processing, so upper-rank compact shards reuse cached hashes instead of
// re-reading FASTAs - the biggest wall-clock win for runs with deep
// taxonomic redundancy.
inline int shard_rank_order(const std::string& rank) {
    if (rank == "species") return 0;
    if (rank == "genus") return 1;
    if (rank == "family") return 2;
    if (rank == "order") return 3;
    if (rank == "class") return 4;
    if (rank == "phylum") return 5;
    return 6;  // unknown / root - process last
}

// True if this shard will use full-base-graph mode (rowCount <= cap AND
// baseGraphBuild requested). Full-base shards are the memory-heavy ones
// (~1-2 GB peak per shard with 30 MB fungal genomes); they get scheduled
// onto a separate worker lane with limited concurrency.
inline bool shard_is_heavy(const ManifestShardInfo& s,
                           size_t maxCladeGenomes,
                           bool baseGraphBuild) {
    return baseGraphBuild && s.rowCount <= maxCladeGenomes && s.rowCount >= 2;
}

// build_partitioned_shards (tiered worker pool)
// Two-lane scheduler that bounds peak RSS:
//   - heavy lane: at most max(1, indexThreads/4) workers; runs the
//     full-base-graph shards (<=cap genomes, reads each FASTA in full,
//     ~1 GB peak per shard).
//   - light lane: the remaining workers; runs compact-mode shards
//     (>cap genomes; only the 512 bp seed per FASTA is read, <=50 MB
//     peak per shard).
// Shards are also rank-sorted so species/genus process before phylum
// (warms the BuilderFastaHashCache for centroid roll-up at upper ranks).
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

    // Build heavy / light index lists in rank order (species first).
    std::vector<size_t> heavyIdx, lightIdx;
    heavyIdx.reserve(shards.size());
    lightIdx.reserve(shards.size());
    for (size_t i = 0; i < shards.size(); ++i) {
        if (shard_is_heavy(shards[i], maxCladeGenomes, baseGraphBuild))
            heavyIdx.push_back(i);
        else
            lightIdx.push_back(i);
    }
    auto byRank = [&](size_t a, size_t b) {
        const int ra = shard_rank_order(shards[a].cladeRank);
        const int rb = shard_rank_order(shards[b].cladeRank);
        if (ra != rb) return ra < rb;
        // Within a rank, smaller shards first -> keeps memory steady and
        // gets centroids into the cache quickly.
        if (shards[a].rowCount != shards[b].rowCount)
            return shards[a].rowCount < shards[b].rowCount;
        return a < b;
    };
    std::sort(heavyIdx.begin(), heavyIdx.end(), byRank);
    std::sort(lightIdx.begin(), lightIdx.end(), byRank);

    std::mutex errMu;
    std::string firstError;
    auto record_error = [&](const std::string& msg) {
        std::lock_guard<std::mutex> lk(errMu);
        if (firstError.empty()) firstError = msg;
    };
    auto have_error = [&]() {
        std::lock_guard<std::mutex> lk(errMu);
        return !firstError.empty();
    };

    auto run_lane = [&](const std::vector<size_t>& laneIdx,
                        std::atomic<size_t>& cursor,
                        const char* /*laneName*/) {
        for (;;) {
            const size_t pos = cursor.fetch_add(1, std::memory_order_relaxed);
            if (pos >= laneIdx.size()) break;
            if (have_error()) break;
            const size_t shardIdx = laneIdx[pos];
            try {
                results[shardIdx] = build_single_manifest_shard(
                    shards[shardIdx], indexDir, registryDir, sp, routingDensity, fb,
                    maxCladeGenomes, baseGraphBuild, verbose);
            } catch (const std::exception& e) {
                record_error(e.what());
                break;
            }
        }
    };

    // Lane sizing: full-base shards can hold ~1 GB per worker, compact
    // shards ~50 MB. With 16 threads + 256 GB cap, allow at most 4 heavy
    // workers (~4 GB peak) so the other 12 can rip through compact shards.
    const size_t totalThreads = std::max<size_t>(1, indexThreads);
    size_t heavyWorkers = std::max<size_t>(1, totalThreads / 4);
    if (heavyIdx.empty()) heavyWorkers = 0;
    if (heavyWorkers > heavyIdx.size()) heavyWorkers = heavyIdx.size();
    size_t lightWorkers = totalThreads > heavyWorkers ? totalThreads - heavyWorkers : 1;
    if (lightIdx.empty()) lightWorkers = 0;
    if (lightWorkers > lightIdx.size()) lightWorkers = lightIdx.size();
    if (heavyWorkers + lightWorkers == 0) heavyWorkers = std::max<size_t>(1, totalThreads);

    if (verbose) {
        std::cerr << "[tol] build_partitioned_shards: " << shards.size()
                  << " shards (" << heavyIdx.size() << " heavy / "
                  << lightIdx.size() << " light), workers="
                  << heavyWorkers << " heavy + " << lightWorkers << " light\n";
    }

    std::atomic<size_t> heavyCursor{0};
    std::atomic<size_t> lightCursor{0};
    std::vector<std::thread> workers;
    workers.reserve(heavyWorkers + lightWorkers);
    for (size_t i = 0; i < heavyWorkers; ++i)
        workers.emplace_back([&]() { run_lane(heavyIdx, heavyCursor, "heavy"); });
    for (size_t i = 0; i < lightWorkers; ++i)
        workers.emplace_back([&]() { run_lane(lightIdx, lightCursor, "light"); });
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

// build_tol_index_from_manifest
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

// build_multi_rank_index_from_manifest
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

// TolGlobal
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
              size_t /*cacheEntries*/,
              const std::unordered_set<std::string>* allowedFastaPaths = nullptr) {
        std::lock_guard<std::mutex> lk(mu_);
        indexDir_    = indexDir;
        registryDir_ = registryDir;
        registry_    = ManifestRegistry(registryDir_);
        registry_.load_from_disk();
        auto registryDescriptors = registry_.descriptors();
        recover_descriptor_fasta_paths_from_hierarchy_manifest(
            registryDescriptors,
            fs::path(registryDir_).parent_path() / "hierarchy_manifest.tsv");
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
        //   - (fasta, contig, cladeName) - deduplicate within the same clade
        //     so a contig listed twice under the same descriptor is not doubled.
        //   - (fasta, contig, DIFFERENT cladeNames) - allowed: different rank
        //     levels must produce separate RefSeq entries with distinct cladeRank
        //     values so hierarchical_call_assembly_multirank can filter by rank.
        std::unordered_map<std::string,
            std::unordered_map<std::string, std::shared_ptr<std::string>>> fastaContigSeq;
        // seen key = fasta + "||" + contig + "||" + cladeName
        std::unordered_set<std::string> seenFastaContigClade;

        const bool debugTolInit = std::getenv("MYCOSV_DEBUG_HIER") != nullptr;
        size_t dbg_descs = 0, dbg_fastas_seen = 0, dbg_skipped_missing = 0,
               dbg_skipped_allowlist = 0, dbg_loaded = 0, dbg_load_threw = 0;
        if (debugTolInit) {
            std::cerr << "[tol-init-dbg] descriptors=" << registryDescriptors.size()
                      << " allowed=" << (allowedFastaPaths ? allowedFastaPaths->size() : 0)
                      << "\n";
        }

        for (const auto& d : registryDescriptors) {
            ++dbg_descs;
            for (const auto& fasta : d.fastaPaths) {
                ++dbg_fastas_seen;
                if (fasta.empty() || !fs::exists(fasta)) { ++dbg_skipped_missing; continue; }
                if (allowedFastaPaths != nullptr &&
                    allowedFastaPaths->find(fasta) == allowedFastaPaths->end()) {
                    ++dbg_skipped_allowlist;
                    continue;
                }
                try {
                    // Load sequences for this FASTA once; reuse on subsequent calls.
                    if (fastaContigSeq.find(fasta) == fastaContigSeq.end()) {
                        auto contigs = read_fasta_local(fasta);
                        auto& pool = fastaContigSeq[fasta];
                        for (auto& kv : contigs)
                            pool[kv.first] = std::make_shared<std::string>(std::move(kv.second));
                        ++dbg_loaded;
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
                } catch (...) { ++dbg_load_threw; }
            }
        }
        if (debugTolInit) {
            std::cerr << "[tol-init-dbg] descs_visited=" << dbg_descs
                      << " fastas_seen=" << dbg_fastas_seen
                      << " skipped_missing=" << dbg_skipped_missing
                      << " skipped_allowlist=" << dbg_skipped_allowlist
                      << " fastas_loaded=" << dbg_loaded
                      << " load_threw=" << dbg_load_threw
                      << " allRefs_=" << allRefs_.size()
                      << " refsByContig_=" << refsByContig_.size()
                      << "\n";
            // Show first 5 descriptors with at least 1 fastaPath and the path itself
            size_t shown = 0;
            for (const auto& d : registryDescriptors) {
                if (shown >= 5) break;
                if (d.fastaPaths.empty()) continue;
                std::cerr << "[tol-init-dbg]   desc clade=" << d.cladeName
                          << " rank=" << d.cladeRank
                          << " nFastas=" << d.fastaPaths.size()
                          << " first=" << d.fastaPaths.front() << "\n";
                ++shown;
            }
            if (allowedFastaPaths && !allowedFastaPaths->empty()) {
                size_t shown2 = 0;
                for (const auto& p : *allowedFastaPaths) {
                    if (shown2 >= 3) break;
                    std::cerr << "[tol-init-dbg]   allowed[" << shown2 << "]=" << p << "\n";
                    ++shown2;
                }
            }
        }
        initialized_ = true;
    }

    bool is_initialized() const { return initialized_; }
    const std::unordered_map<std::string, std::vector<RefSeq>>& refs_by_contig()  const { return refsByContig_; }
    const std::vector<RefSeq>& all_refs()                                          const { return allRefs_;      }
    const ManifestRegistry&    registry()                                           const { return registry_;    }
    bool has_routing_index() const { return routingCladeCount_ > 0; }

    // Synthetic decoy centroids written by augment_routing_store() in the
    // million-real flow carry these prefixes. They have random 64-bit hashes
    // and no backing FASTA - if a decoy wins a routing top-K slot by a
    // coincidental hash collision (small Jaccard denominator inflates ratio),
    // that slot is wasted and effective real-clade coverage shrinks. We strip
    // them here AND request extra slots upstream so the post-filter still
    // returns topK real candidates when decoys polluted the raw top-K.
    static bool is_decoy_clade_name(const std::string& n) {
        return n.rfind("decoy_clade_", 0) == 0;
    }
    static bool is_decoy_phylum_name(const std::string& p) {
        return p.rfind("DecoyPhylum_", 0) == 0;
    }

    std::vector<std::string> route_query_to_clades(std::string_view seq,
                                                   const SyncmerParams& sp,
                                                   const SyncmerParams& fbSp,
                                                   double density,
                                                   size_t topK) const {
        std::vector<std::string> out;
        if (routingCladeCount_ == 0 || seq.empty()) return out;
        auto bounded_route_view = [](std::string_view s) -> std::string {
            constexpr size_t kMaxRouteBp = 1500000;
            constexpr size_t kRouteWindowBp = 500000;
            if (s.size() <= kMaxRouteBp) return std::string(s);
            std::string bounded;
            bounded.reserve(kMaxRouteBp);
            bounded.append(s.substr(0, kRouteWindowBp));
            const size_t mid = s.size() / 2;
            const size_t midStart = mid > kRouteWindowBp / 2 ? mid - kRouteWindowBp / 2 : 0;
            bounded.append(s.substr(midStart, std::min(kRouteWindowBp, s.size() - midStart)));
            const size_t tailStart = s.size() > kRouteWindowBp ? s.size() - kRouteWindowBp : 0;
            bounded.append(s.substr(tailStart));
            return bounded;
        };
        const std::string routeSeq = bounded_route_view(seq);
        const std::string_view routeView(routeSeq.data(), routeSeq.size());

        // Over-fetch when the external store is in play so the real-clade
        // floor still meets topK after decoys are filtered. 4x headroom
        // matches the worst-case decoy contamination we've seen in routing
        // probes on the 1M-centroid store.
        const size_t fetchK = (externalCentroidStore_ != nullptr) ? topK * 4 : topK;

        std::vector<RouteResult> results;
        if (!preferExternalRouting_)
            results = router_.route(routeView, sp, fbSp, density, fetchK);
        if ((results.empty() || preferExternalRouting_) && externalCentroidStore_ != nullptr) {
            auto qc = make_query_centroid_for_routing(routeView, sp, density);
            auto ext = externalCentroidStore_->query_topk_streaming(qc, fetchK);
            results.reserve(ext.size());
            for (auto& r : ext)
                results.push_back({r.cladeName, r.phylum, r.jaccard});
        }
        out.reserve(std::min(topK, results.size()));
        for (const auto& r : results) {
            if (out.size() >= topK) break;
            if (r.cladeName.empty()) continue;
            if (is_decoy_clade_name(r.cladeName)) continue;
            if (is_decoy_phylum_name(r.phylum)) continue;
            out.push_back(r.cladeName);
        }
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

// MultiRankIndex
class MultiRankIndex {
public:
    static MultiRankIndex& instance() {
        static MultiRankIndex inst;
        return inst;
    }
    void init(const std::string& indexDir, const std::string& registryDir,
              size_t cacheBytes, size_t cacheEntries,
              const std::unordered_set<std::string>* allowedFastaPaths = nullptr) {
        TolGlobal::instance().init(indexDir, registryDir, cacheBytes, cacheEntries,
                                   allowedFastaPaths);
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

inline bool is_reads_pseudocontig_name(const std::string& name) {
    return name.find("_mf") != std::string::npos ||
           name.find("_n")  != std::string::npos ||
           name.rfind("sr_unitig", 0) == 0 ||
           name.rfind("lr_pc", 0) == 0;
}

inline std::string refseq_sequence_cache_key(const TolGlobal::RefSeq* r,
                                             int k,
                                             size_t sampleLimit) {
    if (r == nullptr) return {};
    std::string key;
    key.reserve(r->asmName.size() + r->contig.size() + 48);
    key += r->asmName;
    key.push_back('\x1f');
    key += r->contig;
    key.push_back('\x1f');
    key += std::to_string(k);
    key.push_back('\x1f');
    key += std::to_string(sampleLimit);
    return key;
}

inline std::shared_ptr<const std::vector<uint64_t>>
sampled_ref_hashes_cached(const TolGlobal::RefSeq* r,
                          int k,
                          size_t sampleLimit) {
    static const auto empty = std::make_shared<const std::vector<uint64_t>>();
    if (r == nullptr || !r->has_seq() || k <= 0 || sampleLimit == 0) return empty;
    const std::string key = refseq_sequence_cache_key(r, k, sampleLimit);
    static std::mutex mu;
    static std::unordered_map<std::string, std::shared_ptr<const std::vector<uint64_t>>> cache;
    static std::deque<std::string> insertionOrder;
    constexpr size_t kMaxCachedRefSketches = 8192;
    {
        std::lock_guard<std::mutex> lk(mu);
        auto it = cache.find(key);
        if (it != cache.end()) return it->second;
    }

    std::vector<uint64_t> hashes;
    const std::string& refSeq = r->seq();
    if (refSeq.size() >= static_cast<size_t>(k)) {
        const size_t nKmers = refSeq.size() - static_cast<size_t>(k) + 1;
        const size_t step = std::max<size_t>(1, (nKmers + sampleLimit - 1) / sampleLimit);
        hashes.reserve(std::min(nKmers, sampleLimit + 1));
        for (size_t i = 0; i + static_cast<size_t>(k) <= refSeq.size(); i += step)
            hashes.push_back(detail::fnv1a_kmer(refSeq.data() + i, k));
        if (step > 1 && nKmers > 1) {
            const size_t last = nKmers - 1;
            hashes.push_back(detail::fnv1a_kmer(refSeq.data() + last, k));
        }
        std::sort(hashes.begin(), hashes.end());
        hashes.erase(std::unique(hashes.begin(), hashes.end()), hashes.end());
    }

    std::lock_guard<std::mutex> lk(mu);
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;
    auto ptr = std::make_shared<const std::vector<uint64_t>>(std::move(hashes));
    cache.emplace(key, ptr);
    insertionOrder.push_back(key);
    while (cache.size() > kMaxCachedRefSketches && !insertionOrder.empty()) {
        cache.erase(insertionOrder.front());
        insertionOrder.pop_front();
    }
    return ptr;
}

inline std::vector<const TolGlobal::RefSeq*>
top_ref_candidates_for_chaining(const std::string& qSeq,
                                const std::vector<const TolGlobal::RefSeq*>& refs,
                                const FederatedOptions& fo) {
    struct ScoredRef {
        const TolGlobal::RefSeq* ref = nullptr;
        double score = 0.0;
        size_t ordinal = 0;
    };
    auto ref_dedupe_key = [](const TolGlobal::RefSeq* r) {
        std::string key;
        if (r == nullptr) return key;
        key.reserve(r->asmName.size() + r->contig.size() + 2);
        key += r->asmName;
        key.push_back('\x1f');
        key += r->contig;
        return key;
    };
    std::vector<const TolGlobal::RefSeq*> uniqueRefs;
    uniqueRefs.reserve(refs.size());
    std::unordered_set<std::string> seenRefs;
    seenRefs.reserve(refs.size());
    for (const auto* r : refs) {
        if (r == nullptr || !r->has_seq()) continue;
        const std::string key = ref_dedupe_key(r);
        if (seenRefs.insert(key).second) uniqueRefs.push_back(r);
    }

    std::vector<ScoredRef> scored;
    scored.reserve(uniqueRefs.size());
    const int k = std::max(5, std::min(fo.fallbackSketchParams.k > 0
                                       ? fo.fallbackSketchParams.k : 7, 13));
    const auto qHashes = kmer_hashes(qSeq, k);
    if (qHashes.empty()) return {};
    auto sparse_ref_score = [&](const TolGlobal::RefSeq* r) {
        constexpr size_t kMaxSparseRefKmers = 2048;
        const auto refHashesPtr = sampled_ref_hashes_cached(r, k, kMaxSparseRefKmers);
        const auto& refHashes = *refHashesPtr;
        if (refHashes.empty()) return 0.0;
        size_t hits = 0;
        for (uint64_t h : refHashes)
            if (qHashes.find(h) != qHashes.end()) ++hits;
        const double hitFrac = static_cast<double>(hits) / static_cast<double>(refHashes.size());
        const double qLen = static_cast<double>(std::max<size_t>(1, qSeq.size()));
        const double rLen = static_cast<double>(std::max<size_t>(1, r ? r->seq().size() : 0));
        const double lenRatio = std::min(qLen, rLen) / std::max(qLen, rLen);
        return hitFrac + 0.05 * lenRatio;
    };
    size_t ordinal = 0;
    for (const auto* r : uniqueRefs) {
        const double ov = sparse_ref_score(r);
        scored.push_back({r, ov, ordinal++});
    }
    std::sort(scored.begin(), scored.end(),
              [](const ScoredRef& a, const ScoredRef& b) {
                  if (a.score != b.score) return a.score > b.score;
                  return a.ordinal < b.ordinal;
              });
    const size_t keep = std::min(scored.size(),
        std::max<size_t>(8, std::min<size_t>(32, fo.routingTopN * 4)));
    std::vector<const TolGlobal::RefSeq*> out;
    out.reserve(keep);
    for (size_t i = 0; i < keep; ++i) out.push_back(scored[i].ref);
    return out;
}

inline std::string stable_refseq_key(const TolGlobal::RefSeq* r) {
    if (r == nullptr) return {};
    std::string key;
    key.reserve(r->asmName.size() + r->contig.size() + r->clade.size() +
                r->cladeRank.size() + 8);
    key += r->asmName;
    key.push_back('\x1f');
    key += r->contig;
    key.push_back('\x1f');
    key += r->clade;
    key.push_back('\x1f');
    key += r->cladeRank;
    return key;
}

inline void refine_local_chain_breakpoint(SvTypeFromChain::Result& r,
                                          const std::string& qSeq,
                                          const std::string& refSeq) {
    using T = SvTypeFromChain::Type;
    if (qSeq.empty() || refSeq.empty()) return;
    if (r.type != T::INS && r.type != T::DEL) return;
    if (r.qBreakStart < 0 || r.rBreakStart < 0) return;

    int q = std::min(r.qBreakStart, static_cast<int>(qSeq.size()));
    int rp = std::min(r.rBreakStart, static_cast<int>(refSeq.size()));

    // A deliberately small local left-normalization around exact-MEM chain gaps.
    // It improves repeat-adjacent breakpoint consistency without paying for a
    // full DP alignment for every candidate in large AMF panels.
    int shifts = 0;
    constexpr int kMaxRefineShift = 128;
    while (q > 0 && rp > 0 && shifts < kMaxRefineShift &&
           qSeq[static_cast<size_t>(q - 1)] == refSeq[static_cast<size_t>(rp - 1)]) {
        --q;
        --rp;
        ++shifts;
    }
    if (shifts == 0) return;

    r.qBreakStart = q;
    r.rBreakStart = rp;
    if (r.type == T::INS) {
        r.qBreakEnd = std::max(r.qBreakStart, r.qBreakEnd - shifts);
        r.rBreakEnd = r.rBreakStart;
    } else {
        r.qBreakEnd = r.qBreakStart;
        r.rBreakEnd = std::max(r.rBreakStart, r.rBreakEnd - shifts);
    }
}

// make_insdel_call
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

// make_offref_call
// Classifies the element type and stores it in elementClass.
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
        : classify_repeat_element(std::string_view(seq.data(), seq.size()), cladeGc, v.phylum);
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
        // Allow windows down to max(minSvLen, 100): large enough that random
        // k-mer overlap is informative, small enough to surface real microSVs.
        win = std::max<size_t>(
            fo.tolMinBlockBp > 0 ? fo.tolMinBlockBp : 0u,
            static_cast<size_t>(std::max(fo.minSvLen, 100)));
        win = std::min<size_t>(win, std::max<size_t>(100, std::min<size_t>(n, 5000)));
    }
    if (n < win) return out;
    size_t step = std::max<size_t>(1, win / 2);
    constexpr size_t kMaxOffRefWindows = 256;
    if (n > win) {
        const size_t rawWindows = ((n - win) / step) + 1;
        if (rawWindows > kMaxOffRefWindows) {
            step = std::max<size_t>(step, (n - win + kMaxOffRefWindows - 1) / kMaxOffRefWindows);
        }
    }
    const int k = std::max(5, std::min(fo.fallbackSketchParams.k > 0 ? fo.fallbackSketchParams.k : 7, 9));

    std::vector<const TolGlobal::RefSeq*> scanRefs;
    constexpr size_t kMaxOffRefScanRefs = 256;
    if (refs.size() <= kMaxOffRefScanRefs) {
        scanRefs.reserve(refs.size());
        for (const auto& ref : refs)
            if (ref.has_seq()) scanRefs.push_back(&ref);
    } else {
        scanRefs.reserve(kMaxOffRefScanRefs);
        const size_t stride = std::max<size_t>(1, refs.size() / kMaxOffRefScanRefs);
        for (size_t i = 0; i < refs.size() && scanRefs.size() < kMaxOffRefScanRefs; i += stride)
            if (refs[i].has_seq()) scanRefs.push_back(&refs[i]);
        if (scanRefs.empty()) {
            for (const auto& ref : refs) {
                if (!ref.has_seq()) continue;
                scanRefs.push_back(&ref);
                if (scanRefs.size() >= kMaxOffRefScanRefs) break;
            }
        }
    }

    auto sparse_window_ref_overlap = [&](const std::unordered_set<uint64_t>& windowHashes,
                                         const TolGlobal::RefSeq* ref) {
        if (windowHashes.empty() || ref == nullptr || !ref->has_seq()) return 0.0;
        constexpr size_t kMaxSparseRefKmers = 1024;
        const auto refHashesPtr = sampled_ref_hashes_cached(ref, k, kMaxSparseRefKmers);
        const auto& refHashes = *refHashesPtr;
        if (refHashes.empty()) return 0.0;
        size_t hits = 0;
        for (uint64_t h : refHashes)
            if (windowHashes.find(h) != windowHashes.end()) ++hits;
        return static_cast<double>(hits) / static_cast<double>(refHashes.size());
    };

    for (size_t start = 0; start + win <= n; start += step) {
        const std::string_view window(seq.data() + start, win);
        if (is_low_complexity_sequence(window)) continue;
        const auto windowHashes = kmer_hashes(window, k);
        double bestOverlap = 0.0;
        double bestCladeGc = 0.45;
        std::string bestAsm = "OFF_REFERENCE";
        std::string bestRank = ".";
        std::string bestPhylum = ".";
        for (const auto* ref : scanRefs) {
            if (ref == nullptr || !ref->has_seq()) continue;
            const double ov = sparse_window_ref_overlap(windowHashes, ref);
            if (ov > bestOverlap) {
                bestOverlap = ov;
                bestCladeGc = ref->cladeGc;
                bestAsm = ref->clade.empty() ? ref->asmName : ref->clade;
                bestRank = ref->cladeRank.empty() ? "." : ref->cladeRank;
                bestPhylum = ref->phylum.empty() ? "." : ref->phylum;
            }
        }
        const std::string tier = infer_novelty_tier(bestOverlap);
        // OFF_REF_KNOWN (overlap >= 0.20) windows are also SVs in the
        // pangenome sense: reference segments missing from the query reference
        // set or present in a divergent location.
        if (tier != "NOVEL" && tier != "NOVEL_WEAK" && tier != "DIVERGED" &&
            tier != "OFF_REF_KNOWN") continue;
        OffRefWindowCall ow;
        ow.start = start;
        ow.end = start + win;
        ow.bestOverlap = bestOverlap;
        ow.cladeGc = bestCladeGc;
        ow.bestAsm = bestAsm;
        ow.cladeRank = bestRank;
        ow.phylum = bestPhylum;
        ow.tier = tier;
        ow.elementClass = element_class_name(classify_repeat_element(
            window, bestCladeGc, bestPhylum));
        if (!out.empty() && out.back().tier == ow.tier && out.back().bestAsm == ow.bestAsm && out.back().end >= ow.start) {
            out.back().end = ow.end;
            if (out.back().elementClass == "NONE") out.back().elementClass = ow.elementClass;
        } else {
            out.push_back(std::move(ow));
        }
    }
    // NOVEL_WEAK at sparse k-mer sampling (k<=9, <=1024 ref hashes) is a noise
    // band: ANY shared low-complexity / homopolymer / common-motif k-mer lands
    // a random fungal window in 0 < overlap < 0.05. Drop *isolated*
    // single-window NOVEL_WEAK tiles - they uniformly fill the per-contig
    // kMaxOffRefWindows=256 cap with SUPPORT=0 BSCORE=8 stubs (3,584 phantom
    // calls on the F. falciforme vs F. oxysporum benchmark). Genuine
    // NOVEL_WEAK regions span >1 window and survive merging, so keep those.
    // NOVEL/DIVERGED/OFF_REF_KNOWN single-window calls are still informative
    // (overlap==0 or >=0.05) and pass through unchanged.
    if (!out.empty()) {
        std::vector<OffRefWindowCall> kept;
        kept.reserve(out.size());
        for (auto& ow : out) {
            const bool isSingleWindow = (ow.end - ow.start) <= win;
            if (ow.tier == "NOVEL_WEAK" && isSingleWindow) continue;
            kept.push_back(std::move(ow));
        }
        out.swap(kept);
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
    set_variant_sequence_excerpt(v, sub);
    return v;
}

// try_mem_chain_call
// The order permutation is tracked explicitly so chain-to-MEM recovery is O(N)
// instead of quadratic.
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
    const auto chainRefs = top_ref_candidates_for_chaining(qSeq, refCandidates, fo);
    refContigs.reserve(chainRefs.size());
    saRefs.reserve(chainRefs.size());
    const size_t saTextCap = fo.saMaxTextMB > 0 ? fo.saMaxTextMB * 1024 * 1024 : SIZE_MAX;
    size_t saTextAccum = 0;
    for (const auto* r : chainRefs) {
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
        const bool readsPseudo = is_reads_pseudocontig_name(qContig);
        const bool shortPseudo = readsPseudo && qSeq.size() <= 2000;
        const double effectiveMinBlockScore = shortPseudo
            ? std::min(fo.minBlockScore, qSeq.size() <= 300 ? 2.0 : 3.0)
            : fo.minBlockScore;

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
        if (chain.size() < 2 || bestScore < effectiveMinBlockScore)
            return out;

        auto res = SvTypeFromChain::classify(chain, chainRev, sa, fo.minSvLen);

        // DUP fallback: ChainTreap requires strictly increasing rPos and drops
        // backward-mapping MEMs that are the signature of tandem duplications.
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
                    // Allow cross-assembly TRAs: hierarchical TOL SAs contain
                    // multiple reference assemblies per clade, so a query
                    // breakpoint matching different assemblies is a valid
                    // intra-clade TRA or HGT candidate.
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
                set_variant_sequence_excerpt(
                    out.call,
                    std::string_view(qSeq.data() + teS,
                                     static_cast<size_t>(teE - teS)));
                const ElementClass ec = classify_repeat_element(
                    std::string_view(qSeq.data() + teS,
                                     static_cast<size_t>(teE - teS)),
                    gcBg,
                    (primaryRef != nullptr) ? primaryRef->phylum : ".");
                if (ec != ElementClass::NONE)
                    out.call.elementClass = element_class_name(ec);
            }
        } else if (res.type == T2::DEL && primaryRef != nullptr && primaryRef->has_seq()) {
            const double gcBg = primaryRef->cladeGc;
            // SvTypeFromChain::classify returns local-contig coordinates here.
            const int rS = res.rBreakStart;
            const int rE = res.rBreakEnd;
            const auto& refSeq = primaryRef->seq();
            const int safeS = std::max(0, rS);
            const int safeE = std::min(static_cast<int>(refSeq.size()), rE);
            if (safeE > safeS) {
                set_variant_sequence_excerpt(
                    out.call,
                    std::string_view(refSeq.data() + safeS,
                                     static_cast<size_t>(safeE - safeS)));
                const ElementClass ec = classify_repeat_element(
                    std::string_view(refSeq.data() + safeS,
                                     static_cast<size_t>(safeE - safeS)),
                    gcBg,
                    primaryRef->phylum);
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

// multi-ref suffix-array cache
// try_mem_chain_call_multi builds a SuffixArray over the concatenated top-K
// reference contigs. The ref set is selected from the query contig's
// k-mer-similar shortlist, so consecutive query contigs of the same genome
// route to the same or near-identical ref set. SuffixArray::build is
// O(N log^2 N), so cache built SAs across contigs.
//
// This LRU cache keys a built SA by the exact ordered set of RefSeq pointers
// that went into it (pointers are stable for the program lifetime). It is a
// single process-wide cache guarded by a mutex, not thread_local, so the cache
// is bounded once for the whole process. Returned values are shared_ptr, so a
// thread keeps its SA
// alive even if another thread evicts that entry concurrently - only the
// lookup/insert is serialised, the (expensive) SA *use* is lock-free.
struct MultiRefSaCache {
    struct Entry {
        std::string                   key;
        std::shared_ptr<SuffixArray>  sa;
        size_t                        bytes = 0;
    };
    std::mutex         mu_;
    std::vector<Entry> entries_;            // front = most-recently-used
    size_t             totalBytes_ = 0;
    // Process-wide budget. SA footprint is ~13 B/char of concatenated ref
    // text; 2 GiB holds several distinct chromosome-scale ref sets, and the
    // benchmark access pattern (consecutive query contigs of one genome route
    // to the same ref set) means even a handful of entries gives a high hit
    // rate.
    //
    // Runtime-overridable via MYCOSV_SA_CACHE_MB (default 2048). Read-mode
    // pangenome passes iterate up to N benchmark refs per pseudo-contig; with
    // thousands of pseudo-contigs and large ref SAs, small budgets thrash.
    // Sizing the cache to hold all benchmark-ref SAs makes per-contig cost a
    // cache hit. Pure cache sizing: identical calls, fewer rebuilds.
    static size_t max_bytes() {
        static const size_t v = [] {
            size_t mb = 2048;
            if (const char* e = std::getenv("MYCOSV_SA_CACHE_MB"); e && *e) {
                try {
                    long long x = std::stoll(e);
                    if (x > 0) mb = static_cast<size_t>(x);
                } catch (...) {}
            }
            return mb * static_cast<size_t>(1024) * 1024;
        }();
        return v;
    }

    static size_t footprint(const SuffixArray& sa) {
        // text (1B) + sa/lcp/isa (4B each) per character.
        return sa.text.size() * (1 + 4 + 4 + 4);
    }

    std::shared_ptr<SuffixArray>
    get_or_build(const std::string& key,
                 const std::vector<std::pair<std::string, std::string>>& refContigs) {
        {
            std::lock_guard<std::mutex> lk(mu_);
            for (size_t i = 0; i < entries_.size(); ++i) {
                if (entries_[i].key == key) {
                    if (i != 0) std::rotate(entries_.begin(),
                                            entries_.begin() + static_cast<ptrdiff_t>(i),
                                            entries_.begin() + static_cast<ptrdiff_t>(i) + 1);
                    return entries_.front().sa;
                }
            }
        }
        // Build outside the lock so concurrent threads building *different*
        // ref sets do not serialise on each other. A duplicate concurrent
        // build of the same key is harmless - last writer wins, both callers
        // get a valid SA.
        auto sa = std::make_shared<SuffixArray>();
        sa->build(refContigs);
        const size_t bytes = footprint(*sa);
        {
            std::lock_guard<std::mutex> lk(mu_);
            for (size_t i = 0; i < entries_.size(); ++i) {
                if (entries_[i].key == key) {  // another thread won the race
                    if (i != 0) std::rotate(entries_.begin(),
                                            entries_.begin() + static_cast<ptrdiff_t>(i),
                                            entries_.begin() + static_cast<ptrdiff_t>(i) + 1);
                    return entries_.front().sa;
                }
            }
            while (!entries_.empty() && totalBytes_ + bytes > max_bytes()) {
                totalBytes_ -= entries_.back().bytes;
                entries_.pop_back();
            }
            entries_.insert(entries_.begin(), Entry{key, sa, bytes});
            totalBytes_ += bytes;
        }
        return sa;
    }
};

// try_mem_chain_call_multi
// Multi-emit variant of try_mem_chain_call. Uses SvTypeFromChain::classify_all
// to extract every per-pair INS/DEL gap inside a MEM chain instead of only the
// dominant one. On diverged real fungal genomes this is the difference between
// ~5% recall (single-emit) and capturing the 50-500 small SVs per chain that
// minigraph/svim_asm/cactus surface from the alignment ops.
static bool try_mem_chain_call_multi(
        const std::string& qAsm,
        const std::string& qContig,
        const std::string& qSeq,
        const std::vector<const TolGlobal::RefSeq*>& refCandidates,
        const FederatedOptions& fo,
        std::vector<VariantCallBridge>& outCalls) {
    outCalls.clear();
    if (refCandidates.empty() || qSeq.empty()) return false;

    std::vector<const TolGlobal::RefSeq*> saRefs;
    const auto chainRefs = top_ref_candidates_for_chaining(qSeq, refCandidates, fo);
    saRefs.reserve(chainRefs.size());
    const size_t saTextCap = fo.saMaxTextMB > 0 ? fo.saMaxTextMB * 1024 * 1024 : SIZE_MAX;
    size_t saTextAccum = 0;
    for (const auto* r : chainRefs) {
        if (!r->has_seq()) continue;
        if (saTextAccum + r->seq().size() > saTextCap) continue;
        saTextAccum += r->seq().size();
        saRefs.push_back(r);
    }
    if (saRefs.empty()) return false;
    // Canonicalise the ref order before keying/concatenating the SA.  Use a
    // stable biological identity, not RefSeq* addresses: multirank calling
    // builds temporary rankRefs vectors, so pointer addresses can be reused
    // for different references and would otherwise return a stale cached SA.
    std::sort(saRefs.begin(), saRefs.end(),
              [](const TolGlobal::RefSeq* a, const TolGlobal::RefSeq* b) {
                  return stable_refseq_key(a) < stable_refseq_key(b);
              });
    std::vector<std::pair<std::string, std::string>> refContigs;
    refContigs.reserve(saRefs.size());
    std::string saKey;
    for (const auto* r : saRefs) {
        refContigs.push_back({ r->contig, r->seq() });
        saKey += stable_refseq_key(r);
        saKey += ':';
    }
    static MultiRefSaCache saCache;  // process-wide, mutex-guarded
    std::shared_ptr<SuffixArray> saPtr = saCache.get_or_build(saKey, refContigs);
    const SuffixArray& sa = *saPtr;
    const std::string rcSeq = SuffixArray::revcomp(qSeq);

    auto min_mem_from_k = [](int k) {
        return std::max(15, k - 5);
    };
    struct MultiAttempt {
        std::vector<VariantCallBridge> calls;
        double score   = 0.0;
        int    anchors = 0;
        bool   valid   = false;
    };

    auto attempt_chain = [&](int minMem, bool secondaryPass) {
        MultiAttempt out;
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
        const bool readsPseudo = is_reads_pseudocontig_name(qContig);
        const bool shortPseudo = readsPseudo && qSeq.size() <= 2000;
        const double effectiveMinBlockScore = shortPseudo
            ? std::min(fo.minBlockScore, qSeq.size() <= 300 ? 2.0 : 3.0)
            : fo.minBlockScore;

        std::vector<int> order(allMems.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(),
                  [&](int a, int b) {
                      return allMems[static_cast<size_t>(a)].qPos <
                             allMems[static_cast<size_t>(b)].qPos;
                  });

        const int maxGap = fo.chainGapBand > 0 ? fo.chainGapBand : 5000;
        std::unordered_map<uint64_t, size_t> posToMemIdx;
        posToMemIdx.reserve(allMems.size());
        for (size_t mi = 0; mi < allMems.size(); ++mi) {
            uint64_t key = (static_cast<uint64_t>(allMems[mi].qPos) << 32) |
                           static_cast<uint64_t>(static_cast<uint32_t>(allMems[mi].rPos));
            posToMemIdx.emplace(key, mi);
        }

        // Iterative primary-chain extraction. treap.best_chain_path() returns
        // the highest-scoring chain, so a single long MEM can outrank a
        // multi-anchor chain of shorter MEMs. Skip singleton primaries by
        // marking their MEMs consumed and rebuilding until a chain has >=2
        // anchors or the score floor is exhausted.
        std::vector<bool> primaryConsumed(allMems.size(), false);
        std::vector<SuffixArray::Mem> chain;
        std::vector<bool> chainRev;
        double bestScore = 0.0;
        bool primaryFound = false;
        constexpr int kMaxPrimaryAttempts = 32;
        for (int attempt = 0; attempt < kMaxPrimaryAttempts; ++attempt) {
            ChainTreap treap;
            int memsLeft = 0;
            for (int i : order) {
                if (primaryConsumed[static_cast<size_t>(i)]) continue;
                ++memsLeft;
                const auto& m = allMems[static_cast<size_t>(i)];
                treap.insert_and_chain(m.qPos, m.rPos, m.len,
                                       static_cast<float>(m.len), maxGap);
            }
            if (memsLeft == 0) break;

            auto chainIdx = treap.best_chain_path();
            if (chainIdx.empty()) break;

            const double curScore = static_cast<double>(treap.best_chain_score());
            if (curScore < effectiveMinBlockScore) break;

            chain.clear();
            chainRev.clear();
            chain.reserve(chainIdx.size());
            chainRev.reserve(chainIdx.size());
            std::vector<size_t> chainMemIdx;
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
            if (chain.empty()) break;

            // Found a multi-anchor chain (or any chain at all if shortPseudo
            // accepts singletons). Take it.
            if (chain.size() >= 2 || shortPseudo) {
                bestScore = curScore;
                primaryFound = true;
                break;
            }
            // Singleton chain on a non-shortPseudo query: consume it and
            // try the next-best chain.
            for (size_t mi : chainMemIdx) primaryConsumed[mi] = true;
        }
        if (!primaryFound) return out;

        std::vector<SvTypeFromChain::Result> events;
        if (chain.size() < 2) {
            if (!shortPseudo) return out;
            const auto& m = chain.front();
            const int leftFlank = m.qPos;
            const int rightFlank = static_cast<int>(qSeq.size()) - (m.qPos + m.len);
            const int flank = std::max(leftFlank, rightFlank);
            if (flank < fo.minSvLen) return out;
            auto contig_of_single = [&](int rPos) -> int {
                for (int ci = 0; ci < static_cast<int>(sa.contigEnd.size()); ++ci)
                    if (rPos < sa.contigEnd[static_cast<size_t>(ci)]) return ci;
                return -1;
            };
            const int ci = contig_of_single(m.rPos);
            const int ctgOff = (ci > 0) ? sa.contigEnd[static_cast<size_t>(ci) - 1] : 0;
            SvTypeFromChain::Result ev;
            ev.type = SvTypeFromChain::Type::INS;
            ev.svLen = flank;
            if (leftFlank >= rightFlank) {
                ev.qBreakStart = 0;
                ev.qBreakEnd = leftFlank - 1;
                ev.rBreakStart = std::max(0, m.rPos - ctgOff);
            } else {
                ev.qBreakStart = m.qPos + m.len;
                ev.qBreakEnd = static_cast<int>(qSeq.size()) - 1;
                ev.rBreakStart = std::max(0, m.rPos + m.len - ctgOff);
            }
            ev.rBreakEnd = ev.rBreakStart;
            if (ci >= 0) ev.rContig = sa.contigName[static_cast<size_t>(ci)];
            events.push_back(std::move(ev));
        } else {
            events = SvTypeFromChain::classify_all(chain, chainRev, sa, fo.minSvLen);
        }

        // DUP fallback: ChainTreap drops backward-mapping MEMs that signal
        // tandem dups; re-classify the forward MEMs to catch that pattern.
        const bool needsDupRescue = events.empty() ||
            (events.size() == 1 && events.front().type == SvTypeFromChain::Type::INS);
        if (needsDupRescue) {
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
                    if (dupScore >= effectiveMinBlockScore) {
                        events.assign(1, dupRes);
                        chain     = fwdSorted;
                        chainRev.assign(fwdSorted.size(), false);
                        bestScore = dupScore;
                    }
                }
            }
        }

        // TRA fallback: cross-contig anchors whose qGap exceeds chainGapBand never
        // make it into the treap chain; scan allMems directly.
        if (events.empty()) {
            auto ctg_of = [&](int rp) -> int {
                for (int ci = 0; ci < static_cast<int>(sa.contigEnd.size()); ++ci)
                    if (rp < sa.contigEnd[static_cast<size_t>(ci)]) return ci;
                return -1;
            };
            const SuffixArray::Mem* srcMem = nullptr;
            const SuffixArray::Mem* dstMem = nullptr;
            int srcCtg = -1;
            const int traGap = 500000;
            for (int i : order) {
                if (isRev[static_cast<size_t>(i)]) continue;
                const auto& m = allMems[static_cast<size_t>(i)];
                const int ci = ctg_of(m.rPos);
                if (srcCtg < 0) {
                    srcCtg = ci;
                    srcMem = &m;
                    continue;
                }
                if (ci >= 0 && ci != srcCtg) {
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
                    if (traScore >= effectiveMinBlockScore) {
                        events.assign(1, traRes);
                        chain     = traChain;
                        chainRev  = {false, false};
                        bestScore = traScore;
                    }
                }
            }
        }

        if (events.empty()) return out;

        // Merge adjacent small INS/DEL events of the same type within mergeWindow bp
        // to keep the call set comparable to comparator truth (svim_asm collapses
        // microsatellite-adjacent indels into a single record).
        // The window must not scale with minSvLen: TE-burst regions can carry
        // multiple real events within a few hundred bp.
        const int mergeWindow = 80;
        std::sort(events.begin(), events.end(),
                  [](const SvTypeFromChain::Result& a,
                     const SvTypeFromChain::Result& b) {
                      if (a.qBreakStart != b.qBreakStart) return a.qBreakStart < b.qBreakStart;
                      return a.svLen > b.svLen;
                  });
        std::vector<SvTypeFromChain::Result> merged;
        merged.reserve(events.size());
        for (auto& ev : events) {
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

        // Iterative multi-chain extraction: after the primary chain has been
        // converted to events, remove its query footprint and mine additional
        // non-overlapping chains for the same (query contig, reference) pair.
        // This rescues fragmented/read pseudo-contigs where one strong local
        // path can otherwise suppress weaker but valid SV-bearing paths.
        std::vector<std::pair<int,int>> usedQ;
        usedQ.reserve(chain.size());
        for (const auto& m : chain)
            usedQ.push_back({m.qPos, m.qPos + m.len});
        auto overlaps_used_q = [&](const SuffixArray::Mem& m) {
            const int s = m.qPos;
            const int e = m.qPos + m.len;
            for (const auto& iv : usedQ)
                if (s < iv.second && e > iv.first) return true;
            return false;
        };
        // Keep extra local-chain mining bounded. The primary chain plus
        // per-reference passes carry the benchmarkable SV signal; mining one
        // extra chain per 50 kb across a multi-reference bundle caused
        // chromosome-sized fungal contigs to spend most of their runtime in
        // repetitive local chains and inflated raw pangenome observations.
        // One extra chain per 250 kb preserves accessory/TE-burst rescue while
        // keeping the multi-ref pangenome rescue path from dominating runtime.
        const int kMaxExtraChains = std::max(
            12,
            std::min(128, static_cast<int>(qSeq.size() / 250000)));
        for (int extra = 0; extra < kMaxExtraChains; ++extra) {
            std::vector<int> extraOrder;
            extraOrder.reserve(order.size());
            for (int i : order)
                if (!overlaps_used_q(allMems[static_cast<size_t>(i)]))
                    extraOrder.push_back(i);
            if (extraOrder.empty()) break;

            ChainTreap extraTreap;
            std::vector<int> inserted;
            inserted.reserve(extraOrder.size());
            for (int i : extraOrder) {
                const auto& m = allMems[static_cast<size_t>(i)];
                extraTreap.insert_and_chain(m.qPos, m.rPos, m.len,
                                            static_cast<float>(m.len), maxGap);
                inserted.push_back(i);
            }
            const double extraScore = static_cast<double>(extraTreap.best_chain_score());
            if (extraScore < effectiveMinBlockScore) break;
            const auto extraIdx = extraTreap.best_chain_path();
            if (extraIdx.empty()) break;

            std::vector<SuffixArray::Mem> extraChain;
            std::vector<bool> extraRev;
            extraChain.reserve(extraIdx.size());
            extraRev.reserve(extraIdx.size());
            for (int ni : extraIdx) {
                if (ni < 0 || static_cast<size_t>(ni) >= inserted.size()) continue;
                const int mi = inserted[static_cast<size_t>(ni)];
                extraChain.push_back(allMems[static_cast<size_t>(mi)]);
                extraRev.push_back(isRev[static_cast<size_t>(mi)]);
            }
            // treap.best_chain_path() returns the highest-scoring chain. A
            // single long MEM can outrank a multi-anchor chain of shorter MEMs;
            // mark this chain's q-footprint consumed and continue so the next
            // treap rebuild can surface multi-anchor chains. The empty-order,
            // low-score, and loop-bound exits still guarantee termination.
            if (extraChain.size() < 2) {
                for (const auto& m : extraChain)
                    usedQ.push_back({m.qPos, m.qPos + m.len});
                continue;
            }
            auto extraEvents = SvTypeFromChain::classify_all(
                extraChain, extraRev, sa, fo.minSvLen);
            if (extraEvents.empty()) {
                for (const auto& m : extraChain)
                    usedQ.push_back({m.qPos, m.qPos + m.len});
                continue;
            }
            for (auto& ev : extraEvents)
                if (ev.type != SvTypeFromChain::Type::NONE)
                    merged.push_back(std::move(ev));
            for (const auto& m : extraChain)
                usedQ.push_back({m.qPos, m.qPos + m.len});
        }

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

        // TRA recall fix: cross-contig translocation as an ADDITIONAL event
        //   A translocation-bearing query contig still produces a strong
        //   primary chain (the bulk of the contig aligns cleanly), so the
        //   `events.empty()` TRA fallback above never fires for real TRAs - it
        //   only ever caught contigs with NO alignment at all (i.e. noise).
        //   Here we scan the forward MEMs that landed on a reference contig
        //   OTHER than the primary one: a coherent off-contig cluster spanning
        //   >= minSvLen of the query, backed by MEMs over >= half that span,
        //   is a translocated segment. We APPEND it to `merged` rather than
        //   replacing the primary chain.
        if (primaryContigIdx >= 0) {
            struct OffCluster {
                int qLo = std::numeric_limits<int>::max();
                int qHi = std::numeric_limits<int>::min();
                int rLo = std::numeric_limits<int>::max();
                int rHi = std::numeric_limits<int>::min();
                long long cov = 0;
            };
            std::unordered_map<int, OffCluster> offByContig;
            for (int i : order) {
                if (isRev[static_cast<size_t>(i)]) continue;
                const auto& m = allMems[static_cast<size_t>(i)];
                const int ci = contig_of(m.rPos);
                if (ci < 0 || ci == primaryContigIdx) continue;
                auto& oc = offByContig[ci];
                oc.qLo = std::min(oc.qLo, m.qPos);
                oc.qHi = std::max(oc.qHi, m.qPos + m.len);
                oc.rLo = std::min(oc.rLo, m.rPos);
                oc.rHi = std::max(oc.rHi, m.rPos + m.len);
                oc.cov += m.len;
            }
            int bestCi = -1;
            long long bestCov = 0;
            for (const auto& kv : offByContig)
                if (kv.second.cov > bestCov) { bestCov = kv.second.cov; bestCi = kv.first; }
            if (bestCi >= 0) {
                const OffCluster& oc = offByContig[bestCi];
                const int qSpan = oc.qHi - oc.qLo;
                if (qSpan >= fo.minSvLen &&
                    bestCov * 2 >= static_cast<long long>(qSpan)) {
                    const int ctgStart = (bestCi > 0)
                        ? sa.contigEnd[static_cast<size_t>(bestCi) - 1] : 0;
                    SvTypeFromChain::Result tra;
                    tra.type        = SvTypeFromChain::Type::TRA;
                    tra.qBreakStart = oc.qLo;
                    tra.qBreakEnd   = oc.qHi;
                    tra.rBreakStart = oc.rLo - ctgStart;
                    tra.rBreakEnd   = oc.rHi - ctgStart;
                    tra.svLen       = qSpan;
                    tra.rContig     = sa.contigName[static_cast<size_t>(bestCi)];
                    bool dup = false;
                    for (const auto& ev : merged)
                        if (ev.type == SvTypeFromChain::Type::TRA &&
                            ev.rContig == tra.rContig) { dup = true; break; }
                    if (!dup) merged.push_back(tra);
                }
            }
        }

        auto fill_common = [&](VariantCallBridge& v) {
            v.qAsm        = qAsm;
            v.qContig     = qContig;
            v.readSupport = parse_pseudocontig_support(qContig);
            v.refAsm  = (primaryRef != nullptr && !primaryRef->asmName.empty())
                ? primaryRef->asmName
                : ((primaryRef != nullptr && !primaryRef->clade.empty()) ? primaryRef->clade : "unknown");
            v.refContig = (primaryRef != nullptr && !primaryRef->contig.empty())
                ? primaryRef->contig : ".";
            v.refPos = 0;
            v.refEnd = 0;
            v.genotype = "0/1";
            v.gq = std::min(99.0,
                10.0 * std::log10(1.0 + static_cast<double>(chain.size()))
                + 0.5 * bestScore);
            v.blockScore = bestScore;
            v.anchors    = static_cast<int>(chain.size());
            v.alignmentMode = secondaryPass
                ? "mem_chain_ds13_ds18_multi;secondary_seed_rescue"
                : "mem_chain_ds13_ds18_multi";
            v.mapq = 50.0;
            v.annotation = "NONE";
            v.triallelicTopology = ".";
            v.isNonRefVariant = false;
            if (primaryRef != nullptr) {
                v.cladeRank = primaryRef->cladeRank.empty() ? "." : primaryRef->cladeRank;
                v.phylum    = primaryRef->phylum.empty() ? "." : primaryRef->phylum;
            }
        };

        using T = SvTypeFromChain::Type;
        for (auto res : merged) {
            if (res.type == T::NONE) continue;
            const TolGlobal::RefSeq* eventRef = primaryRef;
            if (!res.rContig.empty()) {
                for (const auto* cand : saRefs) {
                    if (cand != nullptr && cand->contig == res.rContig) {
                        eventRef = cand;
                        break;
                    }
                }
            }
            if (eventRef != nullptr && eventRef->has_seq())
                refine_local_chain_breakpoint(res, qSeq, eventRef->seq());
            VariantCallBridge v;
            fill_common(v);
            v.pos   = std::max(1, res.qBreakStart + 1);
            v.end   = std::max(v.pos, res.qBreakEnd >= 0 ? res.qBreakEnd + 1 : v.pos);
            v.svlen = res.svLen;
            if (res.type != T::TRA)
                v.refContig = res.rContig.empty() ? v.refContig : res.rContig;
            // When a multi-ref chain spans contigs from different reference
            // assemblies, use the event contig's refAsm/cladeRank/phylum so
            // dedup and multisample provenance stay keyed to the true anchor.
            if (res.type != T::TRA && eventRef != nullptr && eventRef != primaryRef) {
                v.refAsm = !eventRef->asmName.empty()
                    ? eventRef->asmName
                    : (eventRef->clade.empty() ? v.refAsm : eventRef->clade);
                if (!eventRef->cladeRank.empty()) v.cladeRank = eventRef->cladeRank;
                if (!eventRef->phylum.empty())    v.phylum    = eventRef->phylum;
            }
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
                case T::TRA: {
                    v.type = "TRA"; v.pantreeClass = "NON_REF";
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
                    if (srcCi >= 0 && static_cast<size_t>(srcCi) < saRefs.size()) {
                        const auto* srcRef = saRefs[static_cast<size_t>(srcCi)];
                        if (srcRef != nullptr) {
                            v.refAsm = !srcRef->asmName.empty()
                                ? srcRef->asmName
                                : (srcRef->clade.empty() ? v.refAsm : srcRef->clade);
                            v.refContig = srcRef->contig.empty() ? v.refContig : srcRef->contig;
                            if (!srcRef->cladeRank.empty()) v.cladeRank = srcRef->cladeRank;
                            if (!srcRef->phylum.empty())    v.phylum    = srcRef->phylum;
                        }
                    }
                    v.refPos = srcRPos > 0 ? (srcRPos - srcOff + 1) : 0;
                    v.refEnd = v.refPos;
                    v.mateContig = res.rContig;
                    v.matePos = res.rBreakStart + 1;
                    v.mateEnd = res.rBreakEnd > 0 ? res.rBreakEnd : v.matePos;
                    if (!res.rContig.empty()) {
                        for (const auto* cand : saRefs) {
                            if (cand == nullptr || cand->contig != res.rContig) continue;
                            v.mateRefAsm = !cand->asmName.empty()
                                ? cand->asmName
                                : (cand->clade.empty() ? "." : cand->clade);
                            v.mateOffReference = false;
                            break;
                        }
                    }
                    break;
                }
                default: continue;
            }

            // TE classification (mirror single-emit path).
            using T2 = SvTypeFromChain::Type;
            if (res.type == T2::INS || res.type == T2::DUP || res.type == T2::INV) {
                const double gcBg = (eventRef != nullptr) ? eventRef->cladeGc : 0.45;
                const int teS = res.qBreakStart;
                const int teE = (res.type == T2::INS)
                    ? std::min(teS + res.svLen, static_cast<int>(qSeq.size()))
                    : std::min(res.qBreakEnd + 1, static_cast<int>(qSeq.size()));
                if (teE > teS && teS >= 0) {
                    set_variant_sequence_excerpt(
                        v,
                        std::string_view(qSeq.data() + teS,
                                         static_cast<size_t>(teE - teS)));
                    const ElementClass ec = classify_repeat_element(
                        std::string_view(qSeq.data() + teS,
                                         static_cast<size_t>(teE - teS)),
                        gcBg,
                        (eventRef != nullptr) ? eventRef->phylum : ".");
                    if (ec != ElementClass::NONE)
                        v.elementClass = element_class_name(ec);
                }
            } else if (res.type == T2::DEL && eventRef != nullptr && eventRef->has_seq()) {
                const double gcBg = eventRef->cladeGc;
                const int rS = res.rBreakStart;
                const int rE = res.rBreakEnd;
                const auto& refSeq = eventRef->seq();
                const int safeS = std::max(0, rS);
                const int safeE = std::min(static_cast<int>(refSeq.size()), rE);
                if (safeE > safeS) {
                    set_variant_sequence_excerpt(
                        v,
                        std::string_view(refSeq.data() + safeS,
                                         static_cast<size_t>(safeE - safeS)));
                    const ElementClass ec = classify_repeat_element(
                        std::string_view(refSeq.data() + safeS,
                                         static_cast<size_t>(safeE - safeS)),
                        gcBg,
                        eventRef->phylum);
                    if (ec != ElementClass::NONE)
                        v.elementClass = element_class_name(ec);
                }
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
    MultiAttempt best = attempt_chain(primaryMinMem, false);
    if (fo.useSecondarySeeds) {
        const int secondaryMinMem = min_mem_from_k(fo.secondarySketchParams.k);
        const bool rescueRequested = !best.valid ||
            best.anchors < static_cast<int>(std::max<size_t>(fo.repeatRescueMinAnchors, 2));
        if (secondaryMinMem < primaryMinMem && rescueRequested) {
            MultiAttempt rescue = attempt_chain(secondaryMinMem, true);
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
    outCalls = std::move(best.calls);
    return !outCalls.empty();
}

// per-contig checkpoint hook
// Thread-local: lets the caller (main.cpp process_query) inject a flush
// callback that runs the moment each contig is fully processed inside
// hierarchical_call_assembly, so checkpointing preserves completed contigs
// even if a job ends mid-query.
using PerContigFlushHook = std::function<void(const std::string& qAsm,
                                              const std::string& contigName,
                                              const std::vector<VariantCallBridge>& contigCalls)>;
inline PerContigFlushHook& per_contig_flush_hook() {
    static thread_local PerContigFlushHook hook;
    return hook;
}

struct ScopedPerContigFlushHook {
    PerContigFlushHook previous;

    explicit ScopedPerContigFlushHook(PerContigFlushHook next)
        : previous(std::move(per_contig_flush_hook())) {
        per_contig_flush_hook() = std::move(next);
    }

    ~ScopedPerContigFlushHook() {
        per_contig_flush_hook() = std::move(previous);
    }
};

// RAII guard: runs the hook at the end of each contig iteration regardless
// of which `continue` (Path A / B / C / off-ref) ended the iteration.
struct ContigFlushGuard {
    std::vector<VariantCallBridge>& out;
    size_t startIdx;
    const std::string& qAsm;
    const std::string& contig;
    const std::string* rankStr = nullptr;
    void flush() {
        auto& hook = per_contig_flush_hook();
        if (!hook) return;
        if (out.size() <= startIdx) return;
        std::vector<VariantCallBridge> sub(out.begin() + startIdx, out.end());
        if (rankStr) {
            for (auto& c : sub)
                if (c.alignmentMode.find("rank=") == std::string::npos)
                    c.alignmentMode += ";rank=" + *rankStr;
        }
        try { hook(qAsm, contig, sub); } catch (...) {}
        startIdx = out.size();
    }
    ~ContigFlushGuard() {
        flush();
    }
};

// hierarchical_call_assembly
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
    const int minHierChainAnchors = static_cast<int>(std::max<size_t>(fo.tolMinChainAnchors, 1));
    auto append_window_calls = [&](const std::string& contigName,
                                   const std::string& seq) {
        if (!fo.graphNativeMode) return;
        auto windows = discover_graph_native_offref_windows(seq, allRefs, fo);
        for (const auto& ow : windows) out.push_back(make_offref_window_call(qAsm, contigName, seq, ow));
    };

    // Debug visibility: env var MYCOSV_DEBUG_HIER=1 turns on per-contig
    // breadcrumbs through the Path A/B/C decision tree so we can see why
    // the multisample VCF is dominated by OFF_REF whole-contig calls.
    const bool debugHier = std::getenv("MYCOSV_DEBUG_HIER") != nullptr;
    if (debugHier) {
        std::cerr << "[hier-dbg] qAsm=" << qAsm
                  << " allRefs=" << allRefPtrs.size()
                  << " withSeq=" << std::count_if(allRefPtrs.begin(), allRefPtrs.end(),
                       [](const TolGlobal::RefSeq* r){ return r && r->has_seq(); })
                  << "\n";
    }
    const bool skipPerContigRouting = allRefPtrs.size() <= 50000;
    for (const auto& kv : contigs) {
        const std::string& name = kv.first;
        const std::string& seq  = kv.second;
        // Per-contig flush guard: completed-contig calls are flushed to the
        // checkpoint sink when this iteration ends, regardless of which
        // `continue` path (A/B/C) terminates it.
        ContigFlushGuard _flushGuard{out, out.size(), qAsm, name};
        const auto routed = skipPerContigRouting
            ? std::vector<std::string>{}
            : global.route_query_to_clades(
                seq, fo.primarySketchParams, fo.fallbackSketchParams,
                fo.routingDensity, fo.routingTopN);
        const std::unordered_set<std::string> routedClades(routed.begin(), routed.end());
        auto clade_allowed = [&](const TolGlobal::RefSeq& ref) {
            return routedClades.empty() ||
                   routedClades.find(ref.clade) != routedClades.end() ||
                   routedClades.find(ref.asmName) != routedClades.end();
        };

        // Path A: DS-13+DS-18 MEM chain
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

        std::vector<VariantCallBridge> chainCalls;
        if (debugHier) {
            std::cerr << "[hier-dbg]   contig=" << name
                      << " qLen=" << seq.size()
                      << " routedClades=" << routedClades.size()
                      << " cands=" << cands.size()
                      << " prefilter=start\n";
        }
        const auto chainPrefilter = top_ref_candidates_for_chaining(seq, cands, fo);
        if (debugHier) {
            std::cerr << "[hier-dbg]   contig=" << name
                      << " qLen=" << seq.size()
                      << " routedClades=" << routedClades.size()
                      << " cands=" << cands.size()
                      << " chainPrefilter=" << chainPrefilter.size()
                      << "\n";
        }
        std::vector<VariantCallBridge> multiCalls;
        size_t multiOk = 0, multiKept = 0, singleOk = 0, singleKept = 0;
        if (try_mem_chain_call_multi(qAsm, name, seq, chainPrefilter, fo, multiCalls)) {
            multiOk = multiCalls.size();
            for (auto& c : multiCalls)
                if (c.anchors >= minHierChainAnchors ||
                    (is_reads_pseudocontig_name(name) && c.anchors >= 1)) {
                    chainCalls.push_back(std::move(c));
                    ++multiKept;
                }
        }
        for (const auto* cand : chainPrefilter) {
            std::vector<const TolGlobal::RefSeq*> oneCand{cand};
            std::vector<VariantCallBridge> oneCalls;
            if (try_mem_chain_call_multi(qAsm, name, seq, oneCand, fo, oneCalls)) {
                singleOk += oneCalls.size();
                for (auto& c : oneCalls)
                    if (c.anchors >= minHierChainAnchors ||
                        (is_reads_pseudocontig_name(name) && c.anchors >= 1)) {
                        chainCalls.push_back(std::move(c));
                        ++singleKept;
                    }
            }
        }
        if (debugHier) {
            std::cerr << "[hier-dbg]   PathA multi: "
                      << multiOk << " raw / " << multiKept << " kept (>=anchors "
                      << minHierChainAnchors << "); single-ref: "
                      << singleOk << " raw / " << singleKept
                      << " kept; chainCalls.size=" << chainCalls.size() << "\n";
        }
        if (!chainCalls.empty()) {
            MergeSortTree mst;
            std::vector<std::pair<int,int>> existingIvs;
            existingIvs.reserve(out.size());
            for (const auto& c : out) existingIvs.push_back({c.pos, c.end});
            const bool haveExisting = !existingIvs.empty();
            if (haveExisting) mst.build(existingIvs);
            for (auto& chainCall : chainCalls) {
                if (haveExisting) {
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
            }
            _flushGuard.flush();
            append_window_calls(name, seq);
            continue;
        }

        // Path B: length-delta fallback
        // Primary: name-based lookup (fast, used for assembly contigs).
        // Secondary: k-mer-overlap best match across all refs when the name
        // lookup misses - required for reads-mode pseudo-contigs (lr_pc0,
        // sr_unitig3...) whose names have no correspondence with reference names.
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

            // Accessory-chromosome PAV can exceed half the shorter contig by
            // definition, so length fallback only requires minSvLen/maxSvLen.
            if (best && bestDelta >= fo.minSvLen && bestDelta <= fo.maxSvLen) {
                out.push_back(make_insdel_call(qAsm, name, *best,
                                               static_cast<int>(seq.size()), fo));
                _flushGuard.flush();
                append_window_calls(name, seq);
                continue;
            }
            // Best ref matched but the contig-wide length delta is too small
            // to call as one SV. Emit per-window OFF_REF/INDEL calls so
            // sub-contig variants are not absorbed by Path C's whole-contig
            // fallback.
            if (best) {
                size_t windowsBefore = out.size();
                append_window_calls(name, seq);
                if (out.size() != windowsBefore) continue;
            }
        }

        // Path C: OFF_REF novelty scoring with cross-clade HGT detection
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
        // Path C novelty k: was clamped to 5-9, which made k-mer Jaccard
        // saturate for any two fungal sequences (4^7 = 16384 keys vs 30+ Mb
        // of reference) so the NOVEL tier almost never fired. Raise the
        // ceiling to 17 so the sketch can actually discriminate, and use
        // best-containment in addition to Jaccard so an embedded HGT /
        // STARSHIP block does not get masked by long host flanks.
        const int k = std::max(11, std::min(fo.fallbackSketchParams.k > 0
                                            ? fo.fallbackSketchParams.k : 15, 17));
        for (const auto& ref : allRefs) {
            if (!ref.has_seq()) continue;
            const double ovJ = kmer_overlap_fraction(seq, ref.seq(), k);
            const double ovC = kmer_best_containment_fraction(seq, ref.seq(), k);
            const double ov  = std::max(ovJ, ovC);
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
        // Match score_cross_clade_novelty thresholds: sameClade < 0.10,
        // otherClade >= 0.08, other - same >= 0.05.
        const bool isHgt = hasRoutingCtx &&
            sameCladeOverlap < 0.10 &&
            otherCladeOverlap >= 0.08 &&
            (otherCladeOverlap - sameCladeOverlap) >= 0.05;
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
        // Emit all four tiers including OFF_REF_KNOWN - see note in
        // discover_graph_native_offref_windows.
        if (tier == "NOVEL" || tier == "NOVEL_WEAK" || tier == "DIVERGED" ||
            tier == "OFF_REF_KNOWN") {
            // Prefer per-window OFF_REF calls so a divergent contig produces
            // many small events that can match a per-position truth set.
            // The contig-wide call is kept only as a safety net when graph-
            // native windowing is disabled or yields no windows.
            size_t windowsBefore = out.size();
            append_window_calls(name, seq);
            if (out.size() == windowsBefore) {
                auto ofcall = make_offref_call(qAsm, name, seq, tier,
                                               bestCladeGc, bestAsmName, bestCladeRank, bestPhylum);
                if (isHgt) ofcall.elementClass = "HGT";
                out.push_back(std::move(ofcall));
            } else if (isHgt) {
                for (size_t i = windowsBefore; i < out.size(); ++i)
                    out[i].elementClass = "HGT";
            }
        }
    }

    annotate_pantree_classes(out);
    return out;
}

// hierarchical_call_assembly_multirank
// Routes independently at each Linnaean rank by restricting the reference
// set to contigs whose cladeRank matches the current rank, then merges all
// results and deduplicates by (contig, pos) keeping the call with the
// highest blockScore.
//
// Rank order: phylum -> class -> order -> family -> genus -> species.
// The species-level pass is identical to hierarchical_call_assembly so this
// function is a strict superset of the single-rank path.
inline std::vector<VariantCallBridge>
hierarchical_call_assembly_multirank(
        const std::string& qAsm,
        const std::unordered_map<std::string, std::string>& contigs,
        const FederatedOptions& fo,
        size_t /*routingTopN - fo.routingTopN is used*/) {

    const auto& global   = TolGlobal::instance();
    const auto& allRefs  = global.all_refs();
    const int minHierChainAnchors = static_cast<int>(std::max<size_t>(fo.tolMinChainAnchors, 1));

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

    // Deduplicate by (contig, pos) - key -> index in merged, keyed on best blockScore
    std::unordered_map<std::string, size_t> bestByKey;
    std::vector<VariantCallBridge> merged;

    auto merge_calls = [&](std::vector<VariantCallBridge>&& calls,
                           const std::string& rankStr) {
        for (auto& c : calls) {
            // Stamp the rank that produced this call
            if (c.alignmentMode.find("rank=") == std::string::npos)
                c.alignmentMode += ";rank=" + rankStr;
            // Dedup key must include SV type, svlen bucket, AND refAsm.
            // Earlier history: the (qContig:pos) key collapsed a co-located
            // INS+DEL - common at TE breakpoints - into a single record,
            // dropping every secondary co-located event across ranks. Adding
            // type+svlen fixed that. But omitting refAsm collapses the
            // SAME-type same-position SV against different refs into one
            // record (whichever rank's call wins blockScore), discarding
            // genuine per-(refAsm) signal that the multi-sample VCF needs to
            // surface as separate rows. Cross-ref calls at the same query
            // coordinate are biologically distinct events and must survive.
            const int svlenBucket = std::max(1, std::abs(c.svlen) / 100);
            const std::string key = c.qContig + ":" + std::to_string(c.pos) +
                                    ":" + c.type +
                                    ":" + std::to_string(svlenBucket) +
                                    ":" + c.refAsm;
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
        const bool skipPerContigRouting = rankRefPtrs.size() <= 50000;

        for (const auto& kv : contigs) {
            const std::string& name = kv.first;
            const std::string& seq  = kv.second;
            // Same checkpoint contract as the single-rank caller, but scoped
            // to rankCalls so each completed rank/contig payload is flushed
            // before the rank-wide merge/dedup pass can be interrupted.
            ContigFlushGuard _flushGuard{rankCalls, rankCalls.size(), qAsm, name, &rankStr};
            const auto routed = skipPerContigRouting
                ? std::vector<std::string>{}
                : global.route_query_to_clades(
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

            std::vector<VariantCallBridge> chainCalls;
            const auto chainPrefilter = top_ref_candidates_for_chaining(seq, cands, fo);
            std::vector<VariantCallBridge> multiCalls;
            if (try_mem_chain_call_multi(qAsm, name, seq, chainPrefilter, fo, multiCalls)) {
                for (auto& c : multiCalls)
                    if (c.anchors >= minHierChainAnchors ||
                        (is_reads_pseudocontig_name(name) && c.anchors >= 1))
                        chainCalls.push_back(std::move(c));
            }
            for (const auto* cand : chainPrefilter) {
                std::vector<const TolGlobal::RefSeq*> oneCand{cand};
                std::vector<VariantCallBridge> oneCalls;
                if (try_mem_chain_call_multi(qAsm, name, seq, oneCand, fo, oneCalls)) {
                    for (auto& c : oneCalls)
                        if (c.anchors >= minHierChainAnchors ||
                            (is_reads_pseudocontig_name(name) && c.anchors >= 1))
                            chainCalls.push_back(std::move(c));
                }
            }
            if (!chainCalls.empty()) {
                for (auto& c : chainCalls) rankCalls.push_back(std::move(c));
                _flushGuard.flush();
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

                if (best && bestDelta >= fo.minSvLen && bestDelta <= fo.maxSvLen) {
                    rankCalls.push_back(make_insdel_call(qAsm, name, *best,
                                                         static_cast<int>(seq.size()), fo));
                    _flushGuard.flush();
                    append_rank_window_calls(name, seq);
                    continue;
                }
                if (best) {
                    size_t windowsBefore = rankCalls.size();
                    append_rank_window_calls(name, seq);
                    if (rankCalls.size() != windowsBefore) continue;
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
                // Same Path C novelty-k change as in hierarchical_call_assembly:
                // raise the clamp from 5-9 to 11-17 and combine Jaccard with
                // best-containment so HGT/STARSHIP blocks bordered by host
                // sequence are not masked by long flanks.
                const int k = std::max(11, std::min(fo.fallbackSketchParams.k > 0
                                                    ? fo.fallbackSketchParams.k : 15, 17));
                for (const auto& ref : allRefs) {
                    if (!ref.has_seq()) continue;
                    const double ovJ = kmer_overlap_fraction(seq, ref.seq(), k);
                    const double ovC = kmer_best_containment_fraction(seq, ref.seq(), k);
                    const double ov  = std::max(ovJ, ovC);
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
                // Same thresholds as score_cross_clade_novelty.
                const bool isHgt = hasRoutingCtx &&
                    sameCladeOverlap < 0.10 &&
                    otherCladeOverlap >= 0.08 &&
                    (otherCladeOverlap - sameCladeOverlap) >= 0.05;
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
                if (tier == "NOVEL" || tier == "NOVEL_WEAK" ||
                    tier == "DIVERGED" || tier == "OFF_REF_KNOWN") {
                    size_t windowsBefore = rankCalls.size();
                    append_rank_window_calls(name, seq);
                    if (rankCalls.size() == windowsBefore) {
                        auto ofcall = make_offref_call(qAsm, name, seq, tier,
                                                       bestCladeGc, bestAsmName, bestCladeRank, bestPhylum);
                        if (isHgt) ofcall.elementClass = "HGT";
                        rankCalls.push_back(std::move(ofcall));
                    } else if (isHgt) {
                        for (size_t i = windowsBefore; i < rankCalls.size(); ++i)
                            rankCalls[i].elementClass = "HGT";
                    }
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

// Ancestral alignment helpers

// AncestralManifestContext: holds the manifest path and the clade->rank/phylum
// lookup table built from it.  Populated by load_ancestral_manifest_context.
struct AncestralManifestContext {
    std::string manifestPath;
    // clade name -> (rank, phylum) - built from the manifest at load time
    std::unordered_map<std::string, std::pair<std::string,std::string>> cladeInfo;
};

// Parse the manifest at `path` to build clade->(rank,phylum) lookups.
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
//   query_asm    - assembly name
//   query_contig - contig that carries the call
//   clade        - routed reference clade
//   clade_rank   - Linnaean rank of that clade (from manifest)
//   phylum       - phylum of that clade (from manifest)
//   variant_type - SV type
//   breakpoints  - number of distinct breakpoints implied by the call:
//                    TRA  -> 2  (two chromosomal breakpoints per translocation)
//                    INV  -> 2  (two inversion breakpoints)
//                    OFF_REF -> 1 if enableAncestralRecomb, else 0
//                    INS/DEL/DUP -> 1
//   segment_bp   - total sequence span of the call
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




// Full ancestral sequence reconstruction
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
            // Query-specific deletion relative to reference -> ancestral base is likely retained.
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

// try_mem_chain_call_public
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

inline bool try_mem_chain_call_multi_public(
        const std::string& qAsm,
        const std::string& qContig,
        const std::string& qSeq,
        const std::vector<const TolGlobal::RefSeq*>& refCandidates,
        const FederatedOptions& fo,
        std::vector<VariantCallBridge>& calls) {
    return try_mem_chain_call_multi(qAsm, qContig, qSeq, refCandidates, fo, calls);
}

} // namespace tol

#endif // FUNGI_TOL_BRIDGE_HPP
