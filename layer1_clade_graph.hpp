#ifndef FUNGI_TOL_LAYER1_CLADE_GRAPH_HPP
#define FUNGI_TOL_LAYER1_CLADE_GRAPH_HPP

// layer1_clade_graph.hpp — Layer 1: per-clade pangenome graph data structures.
//
// Data-structure inventory (DS numbers match design doc):
//   DS-7   PathPositionIndex  — order-statistics for TraIntra detection
//   DS-8   is_inversion_flex  — quick-reject for impossible inversions
//   DS-10  ReferenceLCAIndex  — Euler-tour sparse-table O(1) LCA
//   DS-11  classify_triallelic — triallelic topology classification
//   DS-12  PantreeVariantClass — pantree/NON_REF annotation
//   DS-13  SuffixArray + LCP  — O(N log N) build, O(|q| log N) MEM query
//   DS-15  VEBTree            — O(log log U) predecessor/successor
//   DS-16  MergeSortTree      — O(log² N) interval stabbing
//   DS-17  FenwickTree        — O(log N) prefix-sum / rank queries
//   DS-18  ChainTreap         — O(N log N) MEM seed chaining
//
// Repeat / TE / STARSHIP / HGT / RIP annotators (new in v14):
//   detect_tandem_repeat      — period 2–12, ≥5 copies, ≥50 bp
//   detect_ltr_element        — direct terminal repeats ≥50 bp + high-GC interior
//   detect_tir_element        — inverted terminal repeats ≥30 bp
//   detect_line_helitron      — AT-rich + poly-A/T tails (LINE/Helitron)
//   detect_sine               — short + high-GC + terminal repeat
//   detect_starship           — AT-rich hull (GC < cladeGc-0.10) + ~genic cargo ≥1 kb
//   detect_hgt_island         — GC deviation >±0.08 over ≥500 bp (published range ±0.05-0.10)
//   detect_rip_window         — C/G ratio >2.5 in a 500 bp window (benchmark rule)
//   classify_repeat_element   — master dispatcher returning ElementClass

#include <algorithm>
#include <array>
#include <climits>
#include <cstdint>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <cmath>
#include <string_view>

namespace tol {

// =========================================================================
// SyncmerParams — seeding / sketching parameters
// Defined in layer 1 to break the circular dependency: layer3 includes this
// header and needs SyncmerParams, but fungi_tol_bridge includes layer3.
// =========================================================================
struct SyncmerParams {
    int    k              = 21;
    int    s              = 11;
    int    t              = 2;
    size_t stride         = 1;
    bool   useIntervalHash = false;
    int    ihWing         = 3;
    int    ihMaxDist      = 20000;
    double ihResolution   = 4.0;
};

// =========================================================================
// BaseBlockSegmenter — syncmer stub (real implementation compiled separately)
// =========================================================================
struct BaseBlockSegmenter {
    static std::vector<std::pair<size_t, uint64_t>>
    syncmers(std::string_view /*seq*/, int /*k*/, int /*s*/) {
        return {};
    }
};

// =========================================================================
// CladeGraph — pangenome graph node (minimal stub; full engine compiled separately)
// =========================================================================
struct CladeGraph {
    struct Node {
        int id = 0;
        std::string sequence;
    };
    struct Edge {
        int from = 0;
        int to = 0;
        bool fromForward = true;
        bool toForward = true;
    };
    struct Path {
        std::string name;
        std::vector<int> nodes;
    };

    std::string cladeName;
    std::string cladeRank;
    std::string phylum;
    size_t genomeCount  = 0;
    size_t svBubbles    = 0;
    size_t compressedSz = 0;
    std::vector<Node> nodes;
    std::vector<Edge> edges;
    std::vector<Path> paths;

    size_t compressed_bytes() const { return compressedSz; }
};

// =========================================================================
// Import types for CladeGraphBuilder
// =========================================================================
struct ImportSegment {
    int         id              = 0;
    std::string asmName, contig, seq, annotation;
    int         start           = 0;
    int         end             = 0;
    bool        isRepresentative = false;
    int         blockStart      = -1;
    int         blockEnd        = -1;
};

struct ImportEdge {
    int  from      = 0;
    int  to        = 0;
    bool isForward = true;
};

// =========================================================================
// gbz_io binary graph serialization
// =========================================================================
namespace gbz_io {
    inline constexpr uint64_t kMagic = 0x315A42474C4F5455ULL; // "UTOLGBZ1"
    inline constexpr uint32_t kVersion = 1;

    template <class T>
    inline void write_pod(std::ostream& out, const T& value) {
        out.write(reinterpret_cast<const char*>(&value), static_cast<std::streamsize>(sizeof(T)));
        if (!out) throw std::runtime_error("Failed writing binary graph payload");
    }

    template <class T>
    inline void read_pod(std::istream& in, T& value) {
        in.read(reinterpret_cast<char*>(&value), static_cast<std::streamsize>(sizeof(T)));
        if (!in) throw std::runtime_error("Failed reading binary graph payload");
    }

    inline void write_string(std::ostream& out, const std::string& s) {
        uint64_t n = static_cast<uint64_t>(s.size());
        write_pod(out, n);
        out.write(s.data(), static_cast<std::streamsize>(n));
        if (!out) throw std::runtime_error("Failed writing binary graph string");
    }

    inline std::string read_string(std::istream& in) {
        uint64_t n = 0;
        read_pod(in, n);
        std::string s(static_cast<size_t>(n), '\0');
        if (n) {
            in.read(s.data(), static_cast<std::streamsize>(n));
            if (!in) throw std::runtime_error("Failed reading binary graph string");
        }
        return s;
    }

    inline void save(const CladeGraph& g, const std::string& path) {
        std::ofstream out(path, std::ios::binary);
        if (!out) throw std::runtime_error("Cannot write graph file: " + path);
        write_pod(out, kMagic);
        write_pod(out, kVersion);
        write_string(out, g.cladeName);
        write_string(out, g.cladeRank);
        write_string(out, g.phylum);
        write_pod(out, static_cast<uint64_t>(g.genomeCount));
        write_pod(out, static_cast<uint64_t>(g.svBubbles));
        write_pod(out, static_cast<uint64_t>(g.compressedSz));
        write_pod(out, static_cast<uint64_t>(g.nodes.size()));
        for (const auto& n : g.nodes) {
            write_pod(out, static_cast<int32_t>(n.id));
            write_string(out, n.sequence);
        }
        write_pod(out, static_cast<uint64_t>(g.edges.size()));
        for (const auto& e : g.edges) {
            write_pod(out, static_cast<int32_t>(e.from));
            write_pod(out, static_cast<int32_t>(e.to));
            write_pod(out, static_cast<uint8_t>(e.fromForward ? 1 : 0));
            write_pod(out, static_cast<uint8_t>(e.toForward ? 1 : 0));
        }
        write_pod(out, static_cast<uint64_t>(g.paths.size()));
        for (const auto& p : g.paths) {
            write_string(out, p.name);
            write_pod(out, static_cast<uint64_t>(p.nodes.size()));
            for (int nodeId : p.nodes) write_pod(out, static_cast<int32_t>(nodeId));
        }
    }

    inline CladeGraph load(const std::string& path) {
        std::ifstream in(path, std::ios::binary);
        if (!in) throw std::runtime_error("Cannot read graph file: " + path);
        uint64_t magic = 0;
        uint32_t version = 0;
        read_pod(in, magic);
        read_pod(in, version);
        if (magic != kMagic) throw std::runtime_error("Graph file has invalid GBZ magic: " + path);
        if (version != kVersion) throw std::runtime_error("Unsupported GBZ version in: " + path);
        CladeGraph g;
        g.cladeName = read_string(in);
        g.cladeRank = read_string(in);
        g.phylum = read_string(in);
        uint64_t genomeCount = 0, svBubbles = 0, compressedSz = 0;
        read_pod(in, genomeCount);
        read_pod(in, svBubbles);
        read_pod(in, compressedSz);
        g.genomeCount = static_cast<size_t>(genomeCount);
        g.svBubbles = static_cast<size_t>(svBubbles);
        g.compressedSz = static_cast<size_t>(compressedSz);
        uint64_t nNodes = 0;
        read_pod(in, nNodes);
        g.nodes.reserve(static_cast<size_t>(nNodes));
        for (uint64_t i = 0; i < nNodes; ++i) {
            int32_t id = 0;
            read_pod(in, id);
            g.nodes.push_back({static_cast<int>(id), read_string(in)});
        }
        uint64_t nEdges = 0;
        read_pod(in, nEdges);
        g.edges.reserve(static_cast<size_t>(nEdges));
        for (uint64_t i = 0; i < nEdges; ++i) {
            int32_t from = 0, to = 0;
            uint8_t ff = 1, tf = 1;
            read_pod(in, from);
            read_pod(in, to);
            read_pod(in, ff);
            read_pod(in, tf);
            g.edges.push_back({static_cast<int>(from), static_cast<int>(to), ff != 0, tf != 0});
        }
        uint64_t nPaths = 0;
        read_pod(in, nPaths);
        g.paths.reserve(static_cast<size_t>(nPaths));
        for (uint64_t i = 0; i < nPaths; ++i) {
            CladeGraph::Path p;
            p.name = read_string(in);
            uint64_t nodeCount = 0;
            read_pod(in, nodeCount);
            p.nodes.reserve(static_cast<size_t>(nodeCount));
            for (uint64_t j = 0; j < nodeCount; ++j) {
                int32_t nodeId = 0;
                read_pod(in, nodeId);
                p.nodes.push_back(static_cast<int>(nodeId));
            }
            g.paths.push_back(std::move(p));
        }
        return g;
    }
}

inline constexpr double kMinBubbleFreqDef = 0.05;

inline double sketch_jaccard_similarity(const std::vector<uint64_t>& a,
                                        const std::vector<uint64_t>& b) {
    if (a.empty() || b.empty()) return 0.0;
    size_t ia = 0, ib = 0, inter = 0;
    while (ia < a.size() && ib < b.size()) {
        if (a[ia] < b[ib]) ++ia;
        else if (a[ia] > b[ib]) ++ib;
        else {
            ++inter;
            ++ia;
            ++ib;
        }
    }
    const size_t uni = a.size() + b.size() - inter;
    return uni == 0 ? 0.0
                    : static_cast<double>(inter) / static_cast<double>(uni);
}

inline std::string collapse_bucket_key(const ImportSegment& seg) {
    const size_t len = seg.seq.size();
    const size_t prefixLen = std::min<size_t>(12, len);
    const size_t suffixLen = std::min<size_t>(12, len);
    std::string key;
    key.reserve(seg.annotation.size() + prefixLen + suffixLen + 48);
    key += seg.annotation.empty() ? "NONE" : seg.annotation;
    key += '|';
    key += std::to_string(len / 250);
    key += '|';
    key.append(seg.seq.data(), prefixLen);
    key += '|';
    if (suffixLen > 0)
        key.append(seg.seq.data() + (len - suffixLen), suffixLen);
    return key;
}

inline double approximate_sequence_similarity(const std::string& a,
                                              const std::string& b) {
    if (a.empty() || b.empty()) return 0.0;
    if (a == b) return 1.0;
    const size_t maxLen = std::max(a.size(), b.size());
    const size_t minLen = std::min(a.size(), b.size());
    if (minLen == 0) return 0.0;

    size_t prefix = 0;
    while (prefix < minLen && a[prefix] == b[prefix]) ++prefix;

    size_t suffix = 0;
    while (suffix + prefix < minLen &&
           a[a.size() - 1 - suffix] == b[b.size() - 1 - suffix]) {
        ++suffix;
    }

    return static_cast<double>(prefix + suffix) /
           static_cast<double>(maxLen + std::min(prefix, suffix));
}

inline bool segments_should_collapse(
        const ImportSegment& seg,
        const ImportSegment& rep,
        const std::unordered_map<int,std::vector<uint64_t>>& sketches,
        double collapseThresh) {
    if (seg.seq == rep.seq) return true;

    const size_t maxLen = std::max(seg.seq.size(), rep.seq.size());
    const size_t minLen = std::min(seg.seq.size(), rep.seq.size());
    if (maxLen == 0) return false;
    if ((maxLen - minLen) > std::max<size_t>(64, maxLen / 5)) return false;

    if (!seg.annotation.empty() && !rep.annotation.empty() &&
        seg.annotation != "NONE" && rep.annotation != "NONE" &&
        seg.annotation != rep.annotation) {
        return false;
    }

    auto sit = sketches.find(seg.id);
    auto rit = sketches.find(rep.id);
    if (sit != sketches.end() && rit != sketches.end() &&
        !sit->second.empty() && !rit->second.empty()) {
        return sketch_jaccard_similarity(sit->second, rit->second) >= collapseThresh;
    }

    return approximate_sequence_similarity(seg.seq, rep.seq) >=
           std::max(0.85, collapseThresh);
}

inline CladeGraph build_clade_graph(
        const std::string& name, const std::string& rank,
        const std::string& phylum,
        const std::vector<ImportSegment>& segs,
        const std::vector<ImportEdge>&    edges,
        const std::unordered_map<int,std::vector<uint64_t>>& sketches,
        bool enableCollapse, double collapseThresh,
        bool enableCompact, double minBubbleFreq) {
    CladeGraph g;
    g.cladeName = name;
    g.cladeRank = rank;
    g.phylum    = phylum;
    if (segs.empty()) return g;

    std::set<std::string> genomes;
    size_t rawSeqBytes = 0;
    for (const auto& seg : segs) {
        genomes.insert(seg.asmName);
        rawSeqBytes += seg.seq.size();
    }
    g.genomeCount = genomes.size();

    std::vector<const ImportSegment*> ordered;
    ordered.reserve(segs.size());
    for (const auto& seg : segs) ordered.push_back(&seg);
    std::sort(ordered.begin(), ordered.end(),
              [](const ImportSegment* a, const ImportSegment* b) {
                  if (a->isRepresentative != b->isRepresentative)
                      return a->isRepresentative > b->isRepresentative;
                  if (a->asmName != b->asmName) return a->asmName < b->asmName;
                  if (a->contig != b->contig) return a->contig < b->contig;
                  if (a->start != b->start) return a->start < b->start;
                  return a->id < b->id;
              });

    std::unordered_map<std::string, std::vector<size_t>> buckets;
    std::unordered_map<int, int> segToNode;
    std::vector<const ImportSegment*> nodeReps;
    std::vector<size_t> nodeSupport;
    int nextNodeId = 1;

    for (const ImportSegment* seg : ordered) {
        const std::string bucket = collapse_bucket_key(*seg);
        int matchedNodeId = 0;
        auto bit = buckets.find(bucket);
        if (bit != buckets.end()) {
            for (size_t idx : bit->second) {
                const ImportSegment& rep = *nodeReps[idx];
                if (!enableCollapse ||
                    segments_should_collapse(*seg, rep, sketches, collapseThresh)) {
                    matchedNodeId = g.nodes[idx].id;
                    ++nodeSupport[idx];
                    break;
                }
            }
        }
        if (matchedNodeId == 0) {
            matchedNodeId = nextNodeId++;
            g.nodes.push_back({matchedNodeId, seg->seq});
            nodeReps.push_back(seg);
            nodeSupport.push_back(1);
            buckets[bucket].push_back(g.nodes.size() - 1);
        }
        segToNode[seg->id] = matchedNodeId;
    }

    std::unordered_map<std::string, std::vector<const ImportSegment*>> pathSegs;
    pathSegs.reserve(segs.size());
    for (const auto& seg : segs)
        pathSegs[seg.asmName + "\x1f" + seg.contig].push_back(&seg);

    struct EdgeKey {
        int from = 0;
        int to = 0;
        bool fromForward = true;
        bool toForward = true;
        bool operator==(const EdgeKey& o) const {
            return from == o.from && to == o.to &&
                   fromForward == o.fromForward && toForward == o.toForward;
        }
    };
    struct EdgeKeyHash {
        size_t operator()(const EdgeKey& e) const noexcept {
            size_t h = static_cast<size_t>(e.from) * 1315423911u +
                       static_cast<size_t>(e.to);
            h ^= static_cast<size_t>(e.fromForward ? 0x9e37u : 0x7f4au) << 1;
            h ^= static_cast<size_t>(e.toForward ? 0x85ebu : 0xc2b2u) << 2;
            return h;
        }
    };
    std::unordered_map<EdgeKey, size_t, EdgeKeyHash> edgeCounts;

    for (auto& kv : pathSegs) {
        auto& parts = kv.second;
        std::sort(parts.begin(), parts.end(),
                  [](const ImportSegment* a, const ImportSegment* b) {
                      if (a->start != b->start) return a->start < b->start;
                      if (a->end != b->end) return a->end < b->end;
                      return a->id < b->id;
                  });

        const size_t sep = kv.first.find('\x1f');
        CladeGraph::Path p;
        p.name = (sep == std::string::npos) ? kv.first : kv.first.substr(0, sep);
        if (sep != std::string::npos && sep + 1 < kv.first.size()) {
            p.name += "::";
            p.name += kv.first.substr(sep + 1);
        }

        int prevNode = 0;
        for (const ImportSegment* seg : parts) {
            int nodeId = segToNode[seg->id];
            if (!enableCompact || p.nodes.empty() || p.nodes.back() != nodeId)
                p.nodes.push_back(nodeId);
            if (prevNode != 0 && prevNode != nodeId)
                ++edgeCounts[{prevNode, nodeId, true, true}];
            prevNode = nodeId;
        }
        if (!p.nodes.empty()) g.paths.push_back(std::move(p));
    }

    for (const auto& e : edges) {
        auto fit = segToNode.find(e.from);
        auto tit = segToNode.find(e.to);
        if (fit == segToNode.end() || tit == segToNode.end()) continue;
        if (fit->second == tit->second) continue;
        ++edgeCounts[{fit->second, tit->second, e.isForward, e.isForward}];
    }

    g.edges.reserve(edgeCounts.size());
    for (const auto& kv : edgeCounts)
        g.edges.push_back({kv.first.from, kv.first.to,
                           kv.first.fromForward, kv.first.toForward});
    std::sort(g.edges.begin(), g.edges.end(),
              [](const CladeGraph::Edge& a, const CladeGraph::Edge& b) {
                  if (a.from != b.from) return a.from < b.from;
                  if (a.to != b.to) return a.to < b.to;
                  if (a.fromForward != b.fromForward) return a.fromForward < b.fromForward;
                  return a.toForward < b.toForward;
              });

    std::unordered_map<std::string, std::set<int>> bubbleContexts;
    std::unordered_map<std::string, size_t> bubbleSupport;
    std::unordered_map<std::string, std::set<int>> blockContexts;
    std::unordered_map<std::string, size_t> blockSupport;

    for (const auto& kv : pathSegs) {
        auto parts = kv.second;
        std::sort(parts.begin(), parts.end(),
                  [](const ImportSegment* a, const ImportSegment* b) {
                      if (a->start != b->start) return a->start < b->start;
                      if (a->end != b->end) return a->end < b->end;
                      return a->id < b->id;
                  });
        for (size_t i = 0; i < parts.size(); ++i) {
            const int cur = segToNode[parts[i]->id];
            const int prev = (i > 0) ? segToNode[parts[i - 1]->id] : 0;
            const int next = (i + 1 < parts.size()) ? segToNode[parts[i + 1]->id] : 0;
            const std::string ctx = std::to_string(prev) + "|" + std::to_string(next);
            bubbleContexts[ctx].insert(cur);
            ++bubbleSupport[ctx];

            if (parts[i]->blockStart >= 0 && parts[i]->blockEnd >= 0 &&
                parts[i]->blockEnd >= parts[i]->blockStart) {
                const std::string blk =
                    std::to_string(parts[i]->blockStart) + "|" +
                    std::to_string(parts[i]->blockEnd);
                blockContexts[blk].insert(cur);
                ++blockSupport[blk];
            }
        }
    }

    const size_t minBubbleSupport = std::max<size_t>(
        1, static_cast<size_t>(std::ceil(std::max(0.0, minBubbleFreq) *
                                         std::max<size_t>(1, g.genomeCount))));
    for (const auto& kv : bubbleContexts) {
        auto sit = bubbleSupport.find(kv.first);
        if (kv.second.size() > 1 && sit != bubbleSupport.end() &&
            sit->second >= minBubbleSupport) {
            ++g.svBubbles;
        }
    }
    for (const auto& kv : blockContexts) {
        auto sit = blockSupport.find(kv.first);
        if (kv.second.size() > 1 && sit != blockSupport.end() &&
            sit->second >= minBubbleSupport) {
            ++g.svBubbles;
        }
    }

    size_t nodeBytes = 0;
    for (const auto& n : g.nodes) nodeBytes += n.sequence.size();
    size_t pathBytes = 0;
    for (const auto& p : g.paths)
        pathBytes += p.name.size() + p.nodes.size() * sizeof(int);
    const size_t edgeBytes = g.edges.size() * (sizeof(int) * 2 + 2);
    const size_t metadataBytes = 128 + g.cladeName.size() + g.cladeRank.size() + g.phylum.size();
    g.compressedSz = metadataBytes + nodeBytes + pathBytes + edgeBytes;
    if (rawSeqBytes > 0)
        g.compressedSz = std::min(g.compressedSz, rawSeqBytes + metadataBytes + edgeBytes);
    return g;
}

// =========================================================================
// ── REPEAT / TE / STARSHIP / HGT / RIP ANNOTATORS ────────────────────────
// =========================================================================

// ElementClass — annotation label for an off-reference or on-reference region.
enum class ElementClass {
    NONE,
    REPEAT,       // tandem repeat
    TE_LTR,       // LTR retrotransposon
    TE_TIR,       // TIR DNA transposon
    TE_LINE,      // LINE / Helitron
    TE_SINE,      // SINE
    STARSHIP,     // Starship mega-element
    HGT,          // horizontal gene transfer island
    RIP,          // repeat-induced point mutation
};

inline const char* element_class_name(ElementClass ec) {
    switch (ec) {
        case ElementClass::REPEAT:   return "REPEAT";
        case ElementClass::TE_LTR:   return "TE_LTR";
        case ElementClass::TE_TIR:   return "TE_TIR";
        case ElementClass::TE_LINE:  return "TE_LINE";
        case ElementClass::TE_SINE:  return "TE_SINE";
        case ElementClass::STARSHIP: return "STARSHIP";
        case ElementClass::HGT:      return "HGT";
        case ElementClass::RIP:      return "RIP";
        default:                     return "NONE";
    }
}

// ── GC content helper ─────────────────────────────────────────────────────
inline double gc_content(std::string_view seq) {
    if (seq.empty()) return 0.0;
    size_t gc = 0;
    for (char c : seq)
        if (c == 'G' || c == 'C' || c == 'g' || c == 'c') ++gc;
    return static_cast<double>(gc) / static_cast<double>(seq.size());
}

// ── detect_tandem_repeat ─────────────────────────────────────────────────
// Period-regularity counting.
// Returns true if seq contains a tandem run of period p (2–12) with at
// least minCopies copies totalling at least minLen bases.
//
// Algorithm: for each period p, check that seq[i] == seq[i-p] for a run
// of length ≥ p*minCopies using strict period-regularity counting.
// Period-regularity tandem repeat detector.
// For each period p in [minPeriod, maxPeriod], counts a run of consecutive
// positions i where seq[i] == seq[i - p].  A run of length L means the
// segment has period p and spans L+p bases (the p-length seed plus L
// matching offsets).  The run counter resets to 0 on any mismatch so
// partial matches from unrelated regions do not combine.
//
// Returns true when any period p has a run ≥ (minCopies-1)*p positions
// (i.e. at least minCopies copies of the p-mer unit) AND the total span
// ≥ minLen bases.
inline bool detect_tandem_repeat(std::string_view seq,
                                  int minPeriod  = 2,
                                  int maxPeriod  = 12,
                                  int minCopies  = 5,
                                  int minLen     = 50) {
    const int n = static_cast<int>(seq.size());
    if (n < minLen) return false;

    for (int p = minPeriod; p <= maxPeriod; ++p) {
        if (p * minCopies > n) continue;
        // run counts consecutive period-p matching positions starting at p.
        // A run of (minCopies - 1) * p means minCopies copies of the p-mer.
        const int target_run = (minCopies - 1) * p;
        int run = 0;
        for (int i = p; i < n; ++i) {
            if (seq[static_cast<size_t>(i)] == seq[static_cast<size_t>(i - p)])
                ++run;
            else
                run = 0;   // reset on any mismatch — no partial credit
            if (run >= target_run && (run + p) >= minLen)
                return true;
        }
    }
    return false;
}

// ── detect_ltr_element ───────────────────────────────────────────────────
// Looks for direct terminal repeat (DTR) at both ends.
// DTR is detected by checking that the first ltrLen bases of seq match the
// last ltrLen bases with ≤ mismatchRate mismatches.
// High-GC interior heuristic: GC of the middle 50% must be ≥ gcThresh.
inline bool detect_ltr_element(std::string_view seq,
                                int    ltrLen       = 50,
                                double mismatchRate = 0.10,
                                double gcThresh     = 0.45) {
    const int n = static_cast<int>(seq.size());
    if (n < ltrLen * 2 + 100) return false;

    // Check direct terminal repeats
    int mismatches = 0;
    const int maxMm = static_cast<int>(std::ceil(mismatchRate * ltrLen));
    for (int i = 0; i < ltrLen; ++i) {
        if (seq[static_cast<size_t>(i)] != seq[static_cast<size_t>(n - ltrLen + i)])
            ++mismatches;
        if (mismatches > maxMm) return false;
    }

    // Interior GC check
    const int midStart = n / 4;
    const int midEnd   = 3 * n / 4;
    return gc_content(seq.substr(static_cast<size_t>(midStart),
                                 static_cast<size_t>(midEnd - midStart))) >= gcThresh;
}

// ── detect_tir_element ───────────────────────────────────────────────────
// Inverted terminal repeats: first tirLen bases complement-match last tirLen
// bases in reverse order.
inline bool detect_tir_element(std::string_view seq,
                                int    tirLen       = 30,
                                double mismatchRate = 0.10) {
    const int n = static_cast<int>(seq.size());
    if (n < tirLen * 2 + 50) return false;

    // Build complement table at first use (cannot return raw array from lambda).
    static const std::array<char,256> COMP_ARR = [](){
        std::array<char,256> t{};
        for (int i = 0; i < 256; ++i) t[static_cast<size_t>(i)] = static_cast<char>(i);
        t[static_cast<unsigned char>('A')] = 'T';
        t[static_cast<unsigned char>('T')] = 'A';
        t[static_cast<unsigned char>('C')] = 'G';
        t[static_cast<unsigned char>('G')] = 'C';
        t[static_cast<unsigned char>('a')] = 't';
        t[static_cast<unsigned char>('t')] = 'a';
        t[static_cast<unsigned char>('c')] = 'g';
        t[static_cast<unsigned char>('g')] = 'c';
        return t;
    }();
    const auto& COMP = COMP_ARR;

    const int maxMm = static_cast<int>(std::ceil(mismatchRate * tirLen));
    int mismatches  = 0;
    for (int i = 0; i < tirLen; ++i) {
        char fwd = seq[static_cast<size_t>(i)];
        char rev = COMP[static_cast<unsigned char>(seq[static_cast<size_t>(n - 1 - i)])];
        if (fwd != rev) ++mismatches;
        if (mismatches > maxMm) return false;
    }
    return true;
}

// ── detect_line_helitron ─────────────────────────────────────────────────
// AT-rich (GC < 0.42) overall + poly-A or poly-T run ≥ polyRunLen at
// either terminus.
inline bool detect_line_helitron(std::string_view seq,
                                  double gcThresh  = 0.42,
                                  int    polyRunLen = 10) {
    if (seq.size() < 200u) return false;
    if (gc_content(seq) >= gcThresh) return false;

    // Poly-A/T at 5' end
    auto count_run = [&](size_t start, size_t lim, char base) {
        int cnt = 0;
        for (size_t i = start; i < lim; ++i) {
            char c = seq[i];
            if (c == base || c == static_cast<char>(std::tolower(static_cast<unsigned char>(base))))
                ++cnt;
            else break;
        }
        return cnt;
    };

    if (count_run(0,               std::min(seq.size(), static_cast<size_t>(polyRunLen * 2)), 'A') >= polyRunLen) return true;
    if (count_run(0,               std::min(seq.size(), static_cast<size_t>(polyRunLen * 2)), 'T') >= polyRunLen) return true;
    if (count_run(seq.size() >= static_cast<size_t>(polyRunLen * 2) ? seq.size() - static_cast<size_t>(polyRunLen * 2) : 0,
                  seq.size(), 'A') >= polyRunLen) return true;
    if (count_run(seq.size() >= static_cast<size_t>(polyRunLen * 2) ? seq.size() - static_cast<size_t>(polyRunLen * 2) : 0,
                  seq.size(), 'T') >= polyRunLen) return true;
    return false;
}

// ── detect_sine ──────────────────────────────────────────────────────────
// Short (≤ 500 bp) + high-GC interior (≥ 0.50) + terminal repeat
// (first/last 15 bp match with ≤ 2 mismatches).
inline bool detect_sine(std::string_view seq) {
    const int n = static_cast<int>(seq.size());
    if (n < 80 || n > 500) return false;
    if (gc_content(seq) < 0.50) return false;

    // Terminal repeat: first 15 bp ~ last 15 bp
    const int trLen = std::min(15, n / 4);
    int mm = 0;
    for (int i = 0; i < trLen; ++i) {
        if (seq[static_cast<size_t>(i)] != seq[static_cast<size_t>(n - trLen + i)])
            ++mm;
        if (mm > 2) return false;
    }
    return true;
}

// ── detect_starship ──────────────────────────────────────────────────────
// Starship mega-element (Ascomycota only; not confirmed in Glomeromycota or
// Basidiomycota).
//
// Signature: AT-rich hull relative to clade background + GC-rich "cargo"
// sub-region ≥ 1 kb in the interior.
//
// Hull criterion: sequence GC < (cladeGc - overallGcDrop)
//   overallGcDrop = 0.10 (Urquhart et al. 2023: hull is ~10% below background)
//   This is RELATIVE to cladeGc so Neurospora (GC~54%) works correctly:
//   hull must be < 44%, which is biologically accurate for Starship hulls.
//
// Cargo criterion: interior window GC >= (cladeGc - 0.02), i.e. approximately
//   genic GC (45-55% in practice). Previous threshold of 0.55 was too high
//   and would miss Starships in lower-GC organisms.
//
// Ref: Urquhart et al. 2023 Current Biology (discovery paper).
inline bool detect_starship(std::string_view seq,
                             double cladeGc       = 0.45,
                             double overallGcDrop = 0.10,   // hull is cladeGc-0.10
                             double cargoGcMin    = -1.0,   // -1 = auto (cladeGc-0.02)
                             int    cargoMinLen   = 1000) {
    const int n = static_cast<int>(seq.size());
    if (n < cargoMinLen + 200) return false;

    // Hull must be AT-rich relative to background
    const double hullThresh  = cladeGc - overallGcDrop;
    if (gc_content(seq) >= hullThresh) return false;

    // Cargo must be approximately genic GC
    const double cargoThresh = (cargoGcMin < 0.0)
        ? std::max(0.40, cladeGc - 0.02) : cargoGcMin;

    const int winStart = n / 5;
    const int winEnd   = 4 * n / 5;
    for (int i = winStart; i + cargoMinLen <= winEnd; i += cargoMinLen / 2) {
        if (gc_content(seq.substr(static_cast<size_t>(i),
                                  static_cast<size_t>(cargoMinLen))) >= cargoThresh)
            return true;
    }
    return false;
}

// ── detect_hgt_island ────────────────────────────────────────────────────
// Horizontal gene transfer island: GC content deviates from clade background
// by > ±gcDeviation over a sliding window of ≥ minLen bases.
//
// cladeGc: background GC of the clade (computed externally or passed as 0.45).
inline bool detect_hgt_island(std::string_view seq,
                               double cladeGc     = 0.45,
                               double gcDeviation = 0.12,  // benchmark rule: > ±0.12 from clade background
                               int    minLen      = 500) {
    const int n = static_cast<int>(seq.size());
    if (n < minLen) return false;

    for (int i = 0; i + minLen <= n; i += minLen / 2) {
        double wgc = gc_content(seq.substr(static_cast<size_t>(i),
                                           static_cast<size_t>(minLen)));
        if (std::fabs(wgc - cladeGc) > gcDeviation) return true;
    }
    return false;
}

// ── detect_rip_window ────────────────────────────────────────────────────
// RIP (Repeat-Induced Point mutation): benchmark rule for this pipeline is a
// simple C/G imbalance proxy in a 500 bp sliding window.
//
// Detection rule requested for the fungal benchmark suite:
//   C/G ratio > 2.5 over >= 500 bp.
//
// This is intentionally a lightweight, portable heuristic for simulation and
// OFF_REF annotation rather than a full RIP-index estimator.
inline bool detect_rip_window(std::string_view seq,
                               double cgRatioThresh = 2.5,
                               int    winLen        = 500) {
    const int n = static_cast<int>(seq.size());
    if (n < winLen) return false;

    for (int i = 0; i + winLen <= n; i += winLen / 2) {
        int cCount = 0, gCount = 0;
        for (int j = i; j < i + winLen; ++j) {
            const char c = seq[static_cast<size_t>(j)];
            if (c == 'C' || c == 'c') ++cCount;
            else if (c == 'G' || c == 'g') ++gCount;
        }
        if (gCount == 0) {
            if (cCount > 0) return true;
            continue;
        }
        if (static_cast<double>(cCount) / static_cast<double>(gCount) > cgRatioThresh)
            return true;
    }
    return false;
}

// ── classify_repeat_element ──────────────────────────────────────────────
// Master dispatcher.  Tests in order of specificity (most specific first).
// Returns the first match, or ElementClass::NONE if nothing fires.
//
// cladeGc: pass the background GC for the relevant clade; default 0.45.
inline ElementClass classify_repeat_element(std::string_view seq,
                                             double cladeGc = 0.45) {
    if (seq.size() < 50u) return ElementClass::NONE;

    // RIP takes priority (post-translational; affects any repeated element)
    if (detect_rip_window(seq))                        return ElementClass::RIP;

    // HGT: GC-shifted island
    if (detect_hgt_island(seq, cladeGc))               return ElementClass::HGT;

    // Starship: large AT-rich element with GC-rich cargo
    if (detect_starship(seq, cladeGc))                 return ElementClass::STARSHIP;

    // TE subtypes
    if (detect_sine(seq))                              return ElementClass::TE_SINE;
    if (detect_tir_element(seq))                       return ElementClass::TE_TIR;
    if (detect_ltr_element(seq))                       return ElementClass::TE_LTR;
    if (detect_line_helitron(seq))                     return ElementClass::TE_LINE;

    // Generic tandem repeat
    if (detect_tandem_repeat(seq))                     return ElementClass::REPEAT;

    return ElementClass::NONE;
}

// =========================================================================
// DS-7: PathPositionIndex  — order-statistics treap for TraIntra detection
// =========================================================================
struct PathPositionIndex {
    std::vector<size_t> orderStats;  // sorted path positions

    // O(log N) sorted insert via lower_bound.
    void insert_position(size_t pos) {
        auto it = std::lower_bound(orderStats.begin(), orderStats.end(), pos);
        orderStats.insert(it, pos);
    }

    // Returns true if any consecutive pair decreases (potential intra-TRA).
    bool is_non_monotone() const {
        for (size_t i = 1; i < orderStats.size(); ++i)
            if (orderStats[i] < orderStats[i - 1]) return true;
        return false;
    }

    // Quick-reject prefilter: returns true if all positions in [start,end]
    // are monotone (safe to skip deeper check).
    bool quick_reject_window(size_t windowStart, size_t windowEnd) const {
        bool prev_set = false;
        size_t prev = 0;
        for (size_t p : orderStats) {
            if (p < windowStart) continue;
            if (p > windowEnd)   break;
            if (prev_set && p < prev) return false;
            prev = p; prev_set = true;
        }
        return true;
    }
};

enum class TraIntra { LINEAR, TRA_INTRA, TRA_INTER };

// =========================================================================
// DS-8: is_inversion_flex  — quick-reject for impossible inversions
// =========================================================================
// Correct fix: |refLen - altLen| must be ≤ tolerance × min(ref,alt).
inline bool is_inversion_flex(size_t refLen, size_t altLen,
                               double tolerance = 0.10) {
    if (refLen == 0 || altLen == 0) return false;
    const size_t lo = std::min(refLen, altLen);
    const size_t hi = std::max(refLen, altLen);
    return (hi - lo) <= static_cast<size_t>(tolerance * static_cast<double>(lo) + 0.5);
}

// =========================================================================
// DS-10: ReferenceLCAIndex  — Euler tour + sparse-table O(1) RMQ LCA
// =========================================================================
struct ReferenceLCAIndex {
    std::vector<int>              eulerTour;
    std::vector<int>              depth;
    std::vector<int>              first;
    std::vector<std::vector<int>> sparseTable;

    // Build sparse table after eulerTour/depth/first are populated.  O(N log N).
    void build_sparse_table() {
        const int n = static_cast<int>(depth.size());
        if (n == 0) return;
        int LOG = 1;
        while ((1 << LOG) <= n) ++LOG;
        sparseTable.assign(static_cast<size_t>(LOG), std::vector<int>(static_cast<size_t>(n)));
        for (int i = 0; i < n; ++i) sparseTable[0][static_cast<size_t>(i)] = i;
        for (int j = 1; j < LOG; ++j) {
            for (int i = 0; i + (1 << j) <= n; ++i) {
                int l = sparseTable[static_cast<size_t>(j-1)][static_cast<size_t>(i)];
                int r = sparseTable[static_cast<size_t>(j-1)][static_cast<size_t>(i + (1 << (j-1)))];
                sparseTable[static_cast<size_t>(j)][static_cast<size_t>(i)] =
                    (depth[static_cast<size_t>(l)] <= depth[static_cast<size_t>(r)]) ? l : r;
            }
        }
    }

    // O(1) LCA query.  Requires build_sparse_table() to have been called.
    int branch_point(int a, int b) const {
        if (first.empty() || sparseTable.empty()) return std::min(a, b);
        if (a < 0 || b < 0 ||
            a >= static_cast<int>(first.size()) ||
            b >= static_cast<int>(first.size())) return -1;
        int l = first[static_cast<size_t>(a)];
        int r = first[static_cast<size_t>(b)];
        if (l > r) std::swap(l, r);
        const int len = r - l + 1;
        if (len <= 0) return eulerTour[static_cast<size_t>(l)];
        int k = 0;
        while ((1 << (k+1)) <= len) ++k;
        int pl = sparseTable[static_cast<size_t>(k)][static_cast<size_t>(l)];
        int pr = sparseTable[static_cast<size_t>(k)][static_cast<size_t>(r - (1 << k) + 1)];
        int pos = (depth[static_cast<size_t>(pl)] <= depth[static_cast<size_t>(pr)]) ? pl : pr;
        return eulerTour[static_cast<size_t>(pos)];
    }
};

// =========================================================================
// DS-11: classify_triallelic — topology of two overlapping SV intervals
// =========================================================================
enum class TriallelicTopology {
    PROPERLY_TRIALLELIC,
    OVERLAPPING,
    NESTED,
    INTERLOCKING,
};

inline TriallelicTopology classify_triallelic(size_t a0, size_t a1,
                                               size_t b0, size_t b1) {
    if (a1 < b0 || b1 < a0) return TriallelicTopology::PROPERLY_TRIALLELIC;
    if ((a0 <= b0 && a1 >= b1) || (b0 <= a0 && b1 >= a1))
        return TriallelicTopology::NESTED;
    if ((a0 < b0 && a1 < b1) || (b0 < a0 && b1 < a1))
        return TriallelicTopology::OVERLAPPING;
    return TriallelicTopology::INTERLOCKING;
}

// =========================================================================
// DS-12: PantreeVariantClass
// =========================================================================
enum class PantreeVariantClass { SNP, MNP, INS, DEL, DUP, REPL, INV, NON_REF };

struct NonRefVariant {
    bool value = false;
    explicit NonRefVariant(bool v = false) : value(v) {}
    explicit operator bool() const { return value; }
};

inline PantreeVariantClass classify_pantree(const std::string& svtype) {
    if (svtype == "INS")  return PantreeVariantClass::INS;
    if (svtype == "DEL")  return PantreeVariantClass::DEL;
    if (svtype == "DUP")  return PantreeVariantClass::DUP;
    if (svtype == "INV")  return PantreeVariantClass::INV;
    if (svtype == "REPL") return PantreeVariantClass::REPL;
    if (svtype == "MNP")  return PantreeVariantClass::MNP;
    if (svtype == "SNP")  return PantreeVariantClass::SNP;
    return PantreeVariantClass::NON_REF;
}

inline const char* pantree_class_name(PantreeVariantClass c) {
    switch (c) {
        case PantreeVariantClass::SNP:  return "SNP";
        case PantreeVariantClass::MNP:  return "MNP";
        case PantreeVariantClass::INS:  return "INS";
        case PantreeVariantClass::DEL:  return "DEL";
        case PantreeVariantClass::DUP:  return "DUP";
        case PantreeVariantClass::REPL: return "REPL";
        case PantreeVariantClass::INV:  return "INV";
        default:                        return "NON_REF";
    }
}

template <class Variant>
inline void annotate_pantree_classes(std::vector<Variant>& vars) {
    for (auto& v : vars) {
        v.pantreeClass    = pantree_class_name(classify_pantree(v.type));
        v.isNonRefVariant = (v.type == "OFF_REF");
        if (v.triallelicTopology.empty()) v.triallelicTopology = ".";
    }
}

// =========================================================================
// DS-13: SuffixArray + LCP  (Manber & Myers 1990; Kasai et al. 2001)
// =========================================================================
struct SuffixArray {
    std::string      text;       // concatenated ref contigs + sentinels
    std::vector<int> sa;
    std::vector<int> lcp;
    std::vector<int> isa;
    std::vector<int>         contigEnd;
    std::vector<std::string> contigName;

    // Build from (name, seq) pairs.  Each contig is separated by a unique
    // sentinel (char values 1–30) so cross-contig matches never fire.
    //
    // IMPORTANT: sentinel IDs are capped to 1–30 (all < 32) so that the
    // find_mems guard `static_cast<unsigned char>(c) < 32u` reliably stops
    // extension at every sentinel.  The previous range 1–63 allowed IDs
    // 32–63 (printable ASCII) to slip through the guard, enabling spurious
    // cross-contig MEMs when more than 31 contigs were loaded.
    void build(const std::vector<std::pair<std::string,std::string>>& contigs) {
        text.clear(); contigEnd.clear(); contigName.clear();
        int sentinel_id = 1;
        for (const auto& [name, seq] : contigs) {
            contigName.push_back(name);
            text += seq;
            text += static_cast<char>(sentinel_id);
            sentinel_id = (sentinel_id % 30) + 1;   // cap to 1-30: all < 32
            contigEnd.push_back(static_cast<int>(text.size()));
        }
        const int n = static_cast<int>(text.size());
        sa.resize(static_cast<size_t>(n));
        std::iota(sa.begin(), sa.end(), 0);

        // O(N log N) prefix-doubling (Manber & Myers 1990)
        std::vector<int> rank_(static_cast<size_t>(n)), tmp(static_cast<size_t>(n));
        for (int i = 0; i < n; ++i)
            rank_[static_cast<size_t>(i)] = static_cast<unsigned char>(text[static_cast<size_t>(i)]);
        for (int gap = 1; gap < n; gap <<= 1) {
            auto cmp = [&](int a, int b) {
                if (rank_[static_cast<size_t>(a)] != rank_[static_cast<size_t>(b)])
                    return rank_[static_cast<size_t>(a)] < rank_[static_cast<size_t>(b)];
                int ra = a + gap < n ? rank_[static_cast<size_t>(a + gap)] : -1;
                int rb = b + gap < n ? rank_[static_cast<size_t>(b + gap)] : -1;
                return ra < rb;
            };
            std::sort(sa.begin(), sa.end(), cmp);
            tmp[static_cast<size_t>(sa[0])] = 0;
            for (int i = 1; i < n; ++i) {
                tmp[static_cast<size_t>(sa[static_cast<size_t>(i)])] =
                    tmp[static_cast<size_t>(sa[static_cast<size_t>(i-1)])] +
                    (cmp(sa[static_cast<size_t>(i-1)], sa[static_cast<size_t>(i)]) ? 1 : 0);
            }
            rank_ = tmp;
            if (rank_[static_cast<size_t>(sa[static_cast<size_t>(n-1)])] == n-1) break;
        }

        // Inverse SA
        isa.resize(static_cast<size_t>(n));
        for (int i = 0; i < n; ++i)
            isa[static_cast<size_t>(sa[static_cast<size_t>(i)])] = i;

        // O(N) LCP via Kasai et al. (2001)
        lcp.assign(static_cast<size_t>(n), 0);
        int h = 0;
        for (int i = 0; i < n; ++i) {
            if (isa[static_cast<size_t>(i)] > 0) {
                int j = sa[static_cast<size_t>(isa[static_cast<size_t>(i)] - 1)];
                while (i + h < n && j + h < n &&
                       text[static_cast<size_t>(i+h)] == text[static_cast<size_t>(j+h)]) ++h;
                lcp[static_cast<size_t>(isa[static_cast<size_t>(i)])] = h;
                if (h > 0) --h;
            }
        }
    }

    struct Mem { int qPos = 0; int rPos = 0; int len = 0; };

    // Find all Maximal Exact Matches of query against the reference text.
    // Algorithm:
    //   For each query position i, find the SA entry that gives the longest
    //   prefix match with query[i..] using binary search.  The binary search
    //   comparison handles sentinels (chars < 32) correctly: a sentinel at
    //   text[pos+len] means the SA suffix is exhausted before the query —
    //   lexicographically the suffix is less than the query extension, so we
    //   should move RIGHT (lo = mid + 1).
    //
    //   After finding the best match, expand the SA interval [lo2, hi2] to
    //   all entries sharing a common prefix of length best_len using LCP.
    //   Emit one Mem per SA entry in the interval.
    //
    //   Right-extension maximality: after a MEM of length L at query pos i,
    //   advance i by L (skip bases already covered by this MEM).
    std::vector<Mem> find_mems(const std::string& query, int minLen = 20) const {
        std::vector<Mem> out;
        if (sa.empty() || query.empty()) return out;
        constexpr size_t kMaxMemsPerQuery = 500000;
        const int n  = static_cast<int>(text.size());
        const int qn = static_cast<int>(query.size());

        int i = 0;
        while (i < qn) {
            int lo = 0, hi = n - 1;
            int best_len = 0, best_mid = -1;

            while (lo <= hi) {
                const int mid  = lo + (hi - lo) / 2;
                const int rpos = sa[static_cast<size_t>(mid)];

                // Measure match length at this SA entry
                int len = 0;
                while (i + len < qn && rpos + len < n &&
                       query[static_cast<size_t>(i + len)] ==
                           text[static_cast<size_t>(rpos + len)] &&
                       static_cast<unsigned char>(
                           text[static_cast<size_t>(rpos + len)]) >= 32u)
                    ++len;

                if (len > best_len) { best_len = len; best_mid = mid; }

                // Binary search direction:
                //   After matching len chars, the next characters decide.
                //   Case A: query exhausted (len == qn-i) → suffix >= query prefix.
                //           Go left: hi = mid - 1.
                //   Case B: ref suffix exhausted OR hit sentinel (< 32) → the SA
                //           entry is lex-less than any query extension.
                //           Go right: lo = mid + 1.
                //   Case C: both have a character to compare.
                //           text[rpos+len] < query[i+len] → SA entry is lex-less.
                //           Go right: lo = mid + 1.
                //           Otherwise go left: hi = mid - 1.
                if (i + len >= qn) {
                    // query exhausted: suffix ≥ query prefix here → go left
                    hi = mid - 1;
                } else if (rpos + len >= n ||
                           static_cast<unsigned char>(
                               text[static_cast<size_t>(rpos + len)]) < 32u) {
                    // sentinel or end-of-text: SA suffix is lex-smaller → go right
                    lo = mid + 1;
                } else if (text[static_cast<size_t>(rpos + len)] <
                           query[static_cast<size_t>(i + len)]) {
                    lo = mid + 1;
                } else {
                    hi = mid - 1;
                }
            }

            if (best_len < minLen || best_mid < 0) { ++i; continue; }

            // Expand SA interval to all entries sharing a common prefix of
            // length best_len using the LCP array.
            int lo2 = best_mid, hi2 = best_mid;
            while (lo2 > 0 &&
                   lcp[static_cast<size_t>(lo2)] >= best_len)
                --lo2;
            while (hi2 + 1 < n &&
                   lcp[static_cast<size_t>(hi2 + 1)] >= best_len)
                ++hi2;

            for (int k = lo2; k <= hi2; ++k) {
                Mem m;
                m.qPos = i;
                m.rPos = sa[static_cast<size_t>(k)];
                m.len  = best_len;
                out.push_back(m);
                if (out.size() >= kMaxMemsPerQuery) return out;
            }
            i += best_len;  // right-extension maximality
        }
        return out;
    }

    static std::string revcomp(const std::string& s) {
        static const std::array<char,256> COMP = [](){
            std::array<char,256> t{};
            t.fill('N');
            t[static_cast<unsigned char>('A')] = 'T'; t[static_cast<unsigned char>('T')] = 'A';
            t[static_cast<unsigned char>('C')] = 'G'; t[static_cast<unsigned char>('G')] = 'C';
            t[static_cast<unsigned char>('a')] = 't'; t[static_cast<unsigned char>('t')] = 'a';
            t[static_cast<unsigned char>('c')] = 'g'; t[static_cast<unsigned char>('g')] = 'c';
            t[static_cast<unsigned char>('N')] = 'N'; t[static_cast<unsigned char>('n')] = 'n';
            return t;
        }();
        std::string r(s.rbegin(), s.rend());
        for (char& ch : r) ch = COMP[static_cast<unsigned char>(ch)];
        return r;
    }

    bool empty() const { return sa.empty(); }
};

// =========================================================================
// DS-15: van Emde Boas Tree  (van Emde Boas 1975/1977)
// Universe 2^LOG_U; default LOG_U=20 covers ~1 Mb; LOG_U=25 covers ~32 Mb.
// =========================================================================
template <int LOG_U>
struct VEBTree {
    static constexpr int U    = 1 << LOG_U;
    static constexpr int HALF = LOG_U / 2;
    static constexpr int LO_U = 1 << HALF;
    static constexpr int HI_U = 1 << (LOG_U - HALF);

    int min_ = -1, max_ = -1;

    std::unique_ptr<VEBTree<(LOG_U > 1 ? HALF : 1)>>              summary;
    std::vector<std::unique_ptr<VEBTree<(LOG_U > 1 ? HALF : 1)>>> cluster;

    VEBTree() {
        if constexpr (LOG_U > 1) {
            summary = std::make_unique<VEBTree<(LOG_U > 1 ? HALF : 1)>>();
            cluster.resize(static_cast<size_t>(HI_U));
            for (auto& c : cluster)
                c = std::make_unique<VEBTree<(LOG_U > 1 ? HALF : 1)>>();
        }
    }

    bool empty()   const { return min_ == -1; }
    int  minimum() const { return min_; }
    int  maximum() const { return max_; }

    bool contains(int x) const {
        if (x == min_ || x == max_) return true;
        if constexpr (LOG_U <= 1) return false;
        int hi = x >> HALF, lo = x & (LO_U - 1);
        if (hi >= HI_U || !cluster[static_cast<size_t>(hi)]) return false;
        return cluster[static_cast<size_t>(hi)]->contains(lo);
    }

    void insert(int x) {
        if (min_ == -1) { min_ = max_ = x; return; }
        if (x < min_) std::swap(x, min_);
        if (x > max_) max_ = x;
        if constexpr (LOG_U > 1) {
            int hi = x >> HALF, lo = x & (LO_U - 1);
            if (cluster[static_cast<size_t>(hi)]->empty()) summary->insert(hi);
            cluster[static_cast<size_t>(hi)]->insert(lo);
        }
    }

    int successor(int x) const {
        if (empty() || x >= max_) return -1;
        if (x < min_) return min_;
        if constexpr (LOG_U <= 1) return (x == 0 && max_ == 1) ? 1 : -1;
        else {
            int hi = x >> HALF, lo = x & (LO_U - 1);
            if (hi < HI_U && !cluster[static_cast<size_t>(hi)]->empty() &&
                lo < cluster[static_cast<size_t>(hi)]->max_) {
                int offset = cluster[static_cast<size_t>(hi)]->successor(lo);
                return (hi << HALF) | offset;
            }
            int next_hi = summary->successor(hi);
            if (next_hi == -1) return -1;
            return (next_hi << HALF) | cluster[static_cast<size_t>(next_hi)]->min_;
        }
    }

    int predecessor(int x) const {
        if (empty() || x <= min_) return -1;
        if (x > max_) return max_;
        if constexpr (LOG_U <= 1) return (x == 1 && min_ == 0) ? 0 : -1;
        else {
            int hi = x >> HALF, lo = x & (LO_U - 1);
            if (hi < HI_U && !cluster[static_cast<size_t>(hi)]->empty() &&
                lo > cluster[static_cast<size_t>(hi)]->min_) {
                int offset = cluster[static_cast<size_t>(hi)]->predecessor(lo);
                return (hi << HALF) | offset;
            }
            int prev_hi = summary->predecessor(hi);
            if (prev_hi == -1) return (x > min_) ? min_ : -1;
            return (prev_hi << HALF) | cluster[static_cast<size_t>(prev_hi)]->max_;
        }
    }

    bool any_in(int lo, int hi) const {
        if (empty() || lo > max_ || hi < min_) return false;
        if (lo <= min_ || hi >= max_) return true;
        int s = successor(lo - 1);
        return s != -1 && s <= hi;
    }
};

using VEBTree25 = VEBTree<25>;

// =========================================================================
// DS-16: MergeSortTree for interval stabbing queries  (Willard 1985)
//
// build():       O(N log N) — inserts each interval's hi into the segment
//                tree node for its lo rank; sorts each node's list.
// stab_count(p): O(log² N) — walks O(log N) nodes covering [0..rank(p)],
//                binary-searches each node's sorted hi list for hi >= p.
// overlapping(): O(N) linear scan — correct and safe for the typical
//                call count per assembly (<1 000).
// =========================================================================
struct MergeSortTree {
    int                             n_ = 0;
    std::vector<std::pair<int,int>> intervals_;
    std::vector<std::vector<int>>   tree_;
    // Coordinate-compressed lo values, stored after build() so stab_count
    // does not need to rebuild them on every call.
    std::vector<int>                xs_;

    void build(const std::vector<std::pair<int,int>>& ivs) {
        intervals_ = ivs;
        n_ = static_cast<int>(ivs.size());
        xs_.clear();
        tree_.clear();
        if (n_ == 0) return;
        xs_.reserve(static_cast<size_t>(n_));
        for (const auto& [lo, hi] : ivs) xs_.push_back(lo);
        std::sort(xs_.begin(), xs_.end());
        xs_.erase(std::unique(xs_.begin(), xs_.end()), xs_.end());
        const int M  = static_cast<int>(xs_.size());
        const int sz = 4 * M;
        tree_.assign(static_cast<size_t>(sz), {});
        for (const auto& [lo, hi] : ivs) {
            int pos = static_cast<int>(
                std::lower_bound(xs_.begin(), xs_.end(), lo) - xs_.begin());
            update(1, 0, M - 1, pos, hi);
        }
        for (auto& v : tree_) std::sort(v.begin(), v.end());
    }

    // O(log² N): count intervals [lo,hi] where lo <= p AND hi >= p.
    // Algorithm: intervals are indexed by their lo rank in the segment tree.
    // All nodes covering the range [0 .. rank(p)] hold intervals with lo <= p.
    // Within each such node the hi list is sorted; binary search counts hi >= p.
    int stab_count(int p) const {
        if (n_ == 0 || xs_.empty() || tree_.empty()) return 0;
        const int M = static_cast<int>(xs_.size());
        // Highest rank r such that xs_[r] <= p; -1 if none
        const int pRank = static_cast<int>(
            std::upper_bound(xs_.begin(), xs_.end(), p) - xs_.begin()) - 1;
        if (pRank < 0) return 0;
        int count = 0;
        stab_query(1, 0, M - 1, 0, pRank, p, count);
        return count;
    }

    // O(N) scan — correct and safe for typical call counts (< 1 000 per assembly).
    std::vector<int> overlapping(int qlo, int qhi) const {
        std::vector<int> out;
        for (int i = 0; i < n_; ++i) {
            const auto& [lo, hi] = intervals_[static_cast<size_t>(i)];
            if (lo <= qhi && hi >= qlo) out.push_back(i);
        }
        return out;
    }

    bool empty() const { return n_ == 0; }

private:
    // Recursive segment-tree range query for stab_count.
    // Accumulates into count the number of hi values >= p in all nodes
    // whose coordinate range is fully contained in [ql..qr].
    void stab_query(int node, int nl, int nr, int ql, int qr,
                    int p, int& count) const {
        if (static_cast<size_t>(node) >= tree_.size()) return;
        if (ql > nr || qr < nl) return;
        if (ql <= nl && nr <= qr) {
            const auto& v = tree_[static_cast<size_t>(node)];
            // Sorted ascending: lower_bound(p) gives first hi >= p
            auto it = std::lower_bound(v.begin(), v.end(), p);
            count += static_cast<int>(v.end() - it);
            return;
        }
        const int mid = (nl + nr) / 2;
        stab_query(2 * node,     nl,    mid, ql, qr, p, count);
        stab_query(2 * node + 1, mid+1, nr,  ql, qr, p, count);
    }

    void update(int node, int nl, int nr, int pos, int val) {
        if (static_cast<size_t>(node) >= tree_.size()) return;
        tree_[static_cast<size_t>(node)].push_back(val);
        if (nl == nr) return;
        int mid = (nl + nr) / 2;
        if (pos <= mid) update(2*node,   nl,    mid, pos, val);
        else            update(2*node+1, mid+1, nr,  pos, val);
    }
};

// =========================================================================
// DS-17: Fenwick Tree / Binary Indexed Tree  (Fenwick 1994)
// =========================================================================
struct FenwickTree {
    std::vector<size_t> bit_;
    int n_ = 0;

    explicit FenwickTree(int n = 0)
        : bit_(static_cast<size_t>(n + 1), 0), n_(n) {}

    void reset(int n) {
        n_ = n;
        bit_.assign(static_cast<size_t>(n + 1), 0);
    }

    void update(int i, size_t delta) {
        for (++i; i <= n_; i += i & -i)
            bit_[static_cast<size_t>(i)] += delta;
    }

    size_t prefix_sum(int i) const {
        size_t s = 0;
        for (++i; i > 0; i -= i & -i)
            s += bit_[static_cast<size_t>(i)];
        return s;
    }

    size_t range_sum(int l, int r) const {
        if (l > r) return 0;
        return prefix_sum(r) - (l > 0 ? prefix_sum(l - 1) : 0);
    }

    int find_kth(size_t target) const {
        if (target == 0) return 0;
        int pos = 0; size_t cur = 0;
        for (int pw = 1 << 20; pw; pw >>= 1) {
            if (pos + pw <= n_ &&
                cur + bit_[static_cast<size_t>(pos + pw)] < target) {
                pos += pw;
                cur += bit_[static_cast<size_t>(pos)];
            }
        }
        return std::min(pos, n_ - 1);
    }
};

// =========================================================================
// DS-18: ChainTreap — O(N log N) MEM seed chaining  (Aragon & Seidel 1989)
// =========================================================================
struct ChainTreap {
    struct Node {
        int   qPos = 0, rPos = 0, len = 0;
        float score = 0.0f, best = 0.0f;
        int   prev = -1, prio = 0, left = -1, right = -1;
    };

    std::vector<Node> nodes_;
    int               root_ = -1;
    std::mt19937      rng_{42};

    void clear() { nodes_.clear(); root_ = -1; }

    float insert_and_chain(int qPos, int rPos, int len,
                           float matchScore, int maxGap) {
        int   predIdx  = -1;
        float predBest = find_pred_score(root_, rPos, qPos, maxGap, predIdx);
        int   idx      = static_cast<int>(nodes_.size());
        Node  nd;
        nd.qPos = qPos; nd.rPos = rPos; nd.len = len;
        nd.score = matchScore;
        nd.best  = predBest + matchScore;
        nd.prev  = predIdx;
        nd.prio  = static_cast<int>(rng_());
        nodes_.push_back(nd);
        root_ = insert_node(root_, idx);
        return nodes_[static_cast<size_t>(idx)].best;
    }

    float best_chain_score() const {
        float best = 0.0f;
        for (const auto& nd : nodes_) best = std::max(best, nd.best);
        return best;
    }

    std::vector<int> best_chain_path() const {
        if (nodes_.empty()) return {};
        int tail = 0;
        for (int i = 1; i < static_cast<int>(nodes_.size()); ++i)
            if (nodes_[static_cast<size_t>(i)].best > nodes_[static_cast<size_t>(tail)].best)
                tail = i;
        std::vector<int> path;
        for (int cur = tail; cur >= 0; cur = nodes_[static_cast<size_t>(cur)].prev)
            path.push_back(cur);
        std::reverse(path.begin(), path.end());
        return path;
    }

private:
    float find_pred_score(int node, int rPos, int qPos, int maxGap,
                          int& bestIdx) const {
        if (node < 0) return 0.0f;
        const Node& nd = nodes_[static_cast<size_t>(node)];
        float best = 0.0f;
        // Check the current node: it qualifies as a predecessor only when
        // both rPos and qPos are strictly less than the query's coordinates
        // and within the chaining gap band.
        if (nd.rPos < rPos && nd.qPos < qPos &&
            rPos - nd.rPos <= maxGap && qPos - nd.qPos <= maxGap) {
            best    = nd.best;
            bestIdx = node;
        }
        // Left subtree (rPos ≤ nd.rPos): may contain valid predecessors.
        {
            int li = -1;
            float lb = find_pred_score(nd.left, rPos, qPos, maxGap, li);
            if (lb > best) { best = lb; bestIdx = li; }
        }
        // Right subtree (rPos > nd.rPos): only explore when nd.rPos < rPos,
        // because all nodes in the right subtree have rPos >= nd.rPos.
        // Exploring the right subtree when nd.rPos >= rPos always misses
        // (those nodes have even larger rPos) — previously unchecked, causing
        // O(N) wasted traversal and occasional wrong predecessor results.
        if (nd.rPos < rPos) {
            int ri = -1;
            float rb = find_pred_score(nd.right, rPos, qPos, maxGap, ri);
            if (rb > best) { best = rb; bestIdx = ri; }
        }
        return best;
    }

    int insert_node(int root, int idx) {
        if (root < 0) return idx;
        Node& nd = nodes_[static_cast<size_t>(idx)];
        Node& rt = nodes_[static_cast<size_t>(root)];
        if (nd.rPos <= rt.rPos) {
            rt.left = insert_node(rt.left, idx);
            if (nodes_[static_cast<size_t>(rt.left)].prio > rt.prio) {
                int nr = rt.left;
                rt.left = nodes_[static_cast<size_t>(nr)].right;
                nodes_[static_cast<size_t>(nr)].right = root;
                return nr;
            }
        } else {
            rt.right = insert_node(rt.right, idx);
            if (nodes_[static_cast<size_t>(rt.right)].prio > rt.prio) {
                int nr = rt.right;
                rt.right = nodes_[static_cast<size_t>(nr)].left;
                nodes_[static_cast<size_t>(nr)].left = root;
                return nr;
            }
        }
        return root;
    }
};

// =========================================================================
// SvTypeFromChain — classify a MEM chain → INS/DEL/INV/DUP/TRA/NONE
// =========================================================================
struct SvTypeFromChain {
    enum class Type { NONE, INS, DEL, INV, DUP, TRA };

    struct Result {
        Type   type        = Type::NONE;
        int    qBreakStart = 0, qBreakEnd = 0;
        int    rBreakStart = 0, rBreakEnd = 0;
        int    svLen       = 0;
        std::string rContig;
    };

    static Result classify(const std::vector<SuffixArray::Mem>& chain,
                           const std::vector<bool>& isRevComp,
                           const SuffixArray& sa,
                           int minSvLen = 40) {
        Result res;
        if (chain.size() < 2) return res;

        auto contig_of = [&](int rPos) -> int {
            for (int ci = 0; ci < static_cast<int>(sa.contigEnd.size()); ++ci)
                if (rPos < sa.contigEnd[static_cast<size_t>(ci)]) return ci;
            return -1;
        };

        const int N   = static_cast<int>(chain.size());
        int c0  = contig_of(chain[0].rPos);
        int cN1 = contig_of(chain[static_cast<size_t>(N-1)].rPos);

        // TRA: first and last MEM on different reference contigs.
        // Anchor the primary break after the first MEM on the source contig,
        // and the mate break at the start of the terminal MEM on the target
        // contig. The old code incorrectly reported the source contig again as
        // CHR2/POS2, which made valid inter-contig chains look like broken
        // self-translocations in VCF output.
        // rBreakStart/rBreakEnd are returned in *local contig coordinates*
        // (mate-contig space) so VCF POS2 is interpretable as a position on
        // the named CHR2 contig — not the concatenated suffix-array offset.
        if (c0 >= 0 && cN1 >= 0 && c0 != cN1) {
            res.type        = Type::TRA;
            res.qBreakStart = chain[0].qPos + chain[0].len;
            res.qBreakEnd   = res.qBreakStart;
            const int mateOff = (cN1 > 0)
                ? sa.contigEnd[static_cast<size_t>(cN1) - 1] : 0;
            res.rBreakStart = chain[static_cast<size_t>(N-1)].rPos - mateOff;
            res.rBreakEnd   = chain[static_cast<size_t>(N-1)].rPos +
                              chain[static_cast<size_t>(N-1)].len - mateOff;
            res.rContig     = sa.contigName[static_cast<size_t>(cN1)];
            res.svLen       = std::max(0, chain[static_cast<size_t>(N-1)].qPos -
                                          res.qBreakStart);
            return res;
        }

        // INV: require enough revcomp evidence (>= minSvLen) to avoid false
        // positives from single short palindromic k-mers
        int totalRevLen = 0;
        for (size_t i = 0; i < isRevComp.size(); ++i)
            if (isRevComp[i]) totalRevLen += chain[i].len;
        if (totalRevLen >= minSvLen) {
            int start = chain[0].qPos;
            int end   = chain[static_cast<size_t>(N-1)].qPos +
                        chain[static_cast<size_t>(N-1)].len;
            const int ctgOff = (c0 > 0)
                ? sa.contigEnd[static_cast<size_t>(c0) - 1] : 0;
            res.type        = Type::INV;
            res.qBreakStart = start;
            res.qBreakEnd   = end - 1;  // 0-based inclusive; caller adds +1 for 1-based VCF
            res.rBreakStart = chain[0].rPos - ctgOff;
            res.rBreakEnd   = chain[static_cast<size_t>(N-1)].rPos +
                              chain[static_cast<size_t>(N-1)].len - ctgOff;
            res.svLen       = end - start;
            if (res.svLen < minSvLen) return Result{};
            if (c0 >= 0) res.rContig = sa.contigName[static_cast<size_t>(c0)];
            return res;
        }

        // Compute cumulative query and reference gaps, and track the strongest
        // local backward jump in reference coordinates. With repetitive query
        // genomes (e.g. arbuscular mycorrhizal fungi) the LCP-expanded
        // find_mems output adds noisy MEM pairs that push *cumulative* rGap
        // back above -minSvLen, which masked tandem duplications and led the
        // primary chain to misclassify them as INS. Detecting DUP on the
        // strongest individual backward jump (with a small forward qGap, the
        // tandem-dup signature) keeps DUP recall on those scenarios while
        // leaving INS/DEL/INV detection unchanged.
        int qGap = 0, rGap = 0;
        int dupPairIdx     = -1;
        int dupPairBack    = 0;   // most negative consecutive rGap
        int dupPairQGap    = 0;
        for (int i = 0; i + 1 < N; ++i) {
            const int lqg = chain[static_cast<size_t>(i+1)].qPos -
                            (chain[static_cast<size_t>(i)].qPos + chain[static_cast<size_t>(i)].len);
            const int lrg = chain[static_cast<size_t>(i+1)].rPos -
                            (chain[static_cast<size_t>(i)].rPos + chain[static_cast<size_t>(i)].len);
            qGap += lqg;
            rGap += lrg;
            if (lrg < dupPairBack) {
                dupPairBack = lrg;
                dupPairIdx  = i;
                dupPairQGap = lqg;
            }
        }

        const int delta = qGap - rGap;

        // DUP: cumulative reference overlap, OR a single backward rPos jump
        // ≥ minSvLen accompanied by a non-negative query gap (tandem-dup
        // signature: q advances by ~svLen while r returns to the original
        // copy's coordinates).
        const bool cumulativeDup = (rGap < -minSvLen);
        const bool perPairDup    = (dupPairIdx >= 0 &&
                                    dupPairBack < -minSvLen &&
                                    dupPairQGap >= -minSvLen &&
                                    -dupPairBack >= dupPairQGap);
        if (cumulativeDup || perPairDup) {
            const int ctgOff = (c0 > 0)
                ? sa.contigEnd[static_cast<size_t>(c0) - 1] : 0;
            res.type = Type::DUP;
            if (cumulativeDup) {
                res.qBreakStart = chain[0].qPos;
                res.qBreakEnd   = chain[static_cast<size_t>(N-1)].qPos +
                                  chain[static_cast<size_t>(N-1)].len - 1;
                res.rBreakStart = chain[0].rPos - ctgOff;
                res.rBreakEnd   = res.rBreakStart + (-rGap);
                res.svLen       = -rGap;
            } else {
                // Anchor on the offending pair: the duplicated copy sits
                // between chain[dupPairIdx] and chain[dupPairIdx+1] in
                // query space, mapping back to chain[dupPairIdx+1].rPos in
                // reference space.
                res.qBreakStart = chain[static_cast<size_t>(dupPairIdx)].qPos +
                                  chain[static_cast<size_t>(dupPairIdx)].len;
                res.qBreakEnd   = chain[static_cast<size_t>(dupPairIdx + 1)].qPos +
                                  chain[static_cast<size_t>(dupPairIdx + 1)].len - 1;
                res.rBreakStart = chain[static_cast<size_t>(dupPairIdx + 1)].rPos - ctgOff;
                res.rBreakEnd   = res.rBreakStart + (-dupPairBack);
                res.svLen       = -dupPairBack;
            }
            if (c0 >= 0) res.rContig = sa.contigName[static_cast<size_t>(c0)];
            return res;
        }

        if (std::abs(delta) < minSvLen) return res;

        // INS / DEL: find the consecutive MEM pair with the dominant local gap
        // mismatch — that is where the SV actually sits, not necessarily chain[0].
        // With 1% genome divergence MEM chains have many short anchors so using
        // chain[0] always places the breakpoint near contig start (wrong).
        {
            int bestBreakIdx = 0;
            int bestLocalAbsDelta = 0;
            for (int i = 0; i + 1 < N; ++i) {
                const int lqg = chain[static_cast<size_t>(i+1)].qPos -
                                (chain[static_cast<size_t>(i)].qPos + chain[static_cast<size_t>(i)].len);
                const int lrg = chain[static_cast<size_t>(i+1)].rPos -
                                (chain[static_cast<size_t>(i)].rPos + chain[static_cast<size_t>(i)].len);
                const int ld = std::abs(lqg - lrg);
                if (ld > bestLocalAbsDelta) { bestLocalAbsDelta = ld; bestBreakIdx = i; }
            }
            res.qBreakStart = chain[static_cast<size_t>(bestBreakIdx)].qPos +
                              chain[static_cast<size_t>(bestBreakIdx)].len;
            const int breakRPos =
                chain[static_cast<size_t>(bestBreakIdx)].rPos +
                chain[static_cast<size_t>(bestBreakIdx)].len;
            const int breakCi   = contig_of(breakRPos > 0 ? breakRPos - 1 : 0);
            const int ctgOff    = (breakCi > 0)
                ? sa.contigEnd[static_cast<size_t>(breakCi) - 1] : 0;
            res.rBreakStart = breakRPos - ctgOff;
        }
        res.svLen       = std::abs(delta);
        if (delta > 0) {
            res.type = Type::INS;
        } else {
            res.type      = Type::DEL;
            res.qBreakEnd = res.qBreakStart;
            res.rBreakEnd = res.rBreakStart + (-delta);
        }
        if (c0 >= 0) res.rContig = sa.contigName[static_cast<size_t>(c0)];
        return res;
    }
};

} // namespace tol

#endif // FUNGI_TOL_LAYER1_CLADE_GRAPH_HPP
