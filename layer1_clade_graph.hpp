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
//   detect_hgt_island         — GC deviation >±0.10 over ≥500 bp
//   detect_rip_window         — lightweight RIP-product index in a 500 bp window
//   classify_repeat_element   — master dispatcher returning ElementClass

#include <algorithm>
#include <array>
#include <climits>
#include <cstdint>
#include <cctype>
#include <deque>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
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
// BaseBlockSegmenter — Hong & Buhler (2016) open syncmer seeding.
//
// A k-mer is a syncmer iff the minimum-hash s-mer among its (k - s + 1)
// constituent s-mers occurs at offset `t` inside the k-mer. Hashing is
// canonical (lex-min of forward and reverse-complement) so seeds are
// strand-agnostic. Implementation: monotonic-deque sliding minimum,
// O(N) over the input sequence. Defaults match SyncmerParams (k=21, s=11,
// t=2) so existing call sites passing only (seq, k, s) keep working.
// k-mers containing any non-ACGT base are skipped.
// =========================================================================
struct BaseBlockSegmenter {
    static std::vector<std::pair<size_t, uint64_t>>
    syncmers(std::string_view seq, int k, int s, int t = 2) {
        std::vector<std::pair<size_t, uint64_t>> out;
        const size_t N = seq.size();
        if (k <= 0 || s <= 0 || s > k || N < static_cast<size_t>(k)) return out;
        const int W = k - s + 1;                  // s-mers per k-mer
        if (t < 0 || t > k - s) t = (k - s) / 2;

        auto complement = [](char c) -> char {
            switch (c) {
                case 'A': case 'a': return 'T';
                case 'C': case 'c': return 'G';
                case 'G': case 'g': return 'C';
                case 'T': case 't': return 'A';
                default: return 'N';
            }
        };
        auto fnv1a = [](const char* p, size_t len) -> uint64_t {
            uint64_t h = 14695981039346656037ULL;
            for (size_t i = 0; i < len; ++i) {
                h ^= static_cast<uint64_t>(static_cast<unsigned char>(p[i]));
                h *= 1099511628211ULL;
            }
            return h;
        };
        auto canonical_hash = [&](size_t pos, int len) -> uint64_t {
            std::string fwd(seq.data() + pos, static_cast<size_t>(len));
            std::string rev(static_cast<size_t>(len), 'N');
            for (int i = 0; i < len; ++i)
                rev[static_cast<size_t>(len - 1 - i)] =
                    complement(fwd[static_cast<size_t>(i)]);
            const std::string& canon = (fwd <= rev) ? fwd : rev;
            return fnv1a(canon.data(), canon.size());
        };

        // Cumulative count of non-ACGT bases for O(1) window validity checks.
        std::vector<size_t> nPref(N + 1, 0);
        for (size_t i = 0; i < N; ++i) {
            const char c = seq[i];
            const bool isACGT = (c == 'A' || c == 'C' || c == 'G' || c == 'T'
                              || c == 'a' || c == 'c' || c == 'g' || c == 't');
            nPref[i + 1] = nPref[i] + (isACGT ? 0 : 1);
        }

        const size_t M = N - static_cast<size_t>(s) + 1;   // # s-mer positions
        std::vector<uint64_t> sHash(M, 0);
        std::vector<unsigned char> sOk(M, 0);
        for (size_t i = 0; i < M; ++i) {
            if (nPref[i + static_cast<size_t>(s)] - nPref[i] == 0) {
                sOk[i] = 1;
                sHash[i] = canonical_hash(i, s);
            }
        }

        // Monotonic deque: front = argmin s-mer index in the current k-mer
        // window. Any N inside an s-mer breaks the chain (clear deque).
        std::deque<size_t> dq;
        out.reserve(M / static_cast<size_t>(W) + 8);
        for (size_t i = 0; i < M; ++i) {
            if (!sOk[i]) { dq.clear(); continue; }
            while (!dq.empty() && sHash[dq.back()] >= sHash[i]) dq.pop_back();
            dq.push_back(i);
            if (i + 1 < static_cast<size_t>(W)) continue;
            const size_t kStart = i + 1 - static_cast<size_t>(W);
            while (!dq.empty() && dq.front() < kStart) dq.pop_front();
            if (dq.empty()) continue;
            if (nPref[kStart + static_cast<size_t>(k)] - nPref[kStart] != 0)
                continue;
            if (dq.front() == kStart + static_cast<size_t>(t)) {
                out.emplace_back(kStart, canonical_hash(kStart, k));
            }
        }
        return out;
    }
};

// =========================================================================
// CladeGraph — per-clade pangenome graph (nodes, oriented edges, paths).
// Built by build_clade_graph() below (collapses ImportSegments into nodes via
// sketch-Jaccard buckets, reconstructs per-asm/contig paths, aggregates
// oriented edge counts, and tallies bubble / block contexts that downstream
// callers use as SV-bubble candidates). CladeGraphBuilder in
// layer2_registry.hpp wraps it for streaming construction.
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
// Strict period-regularity detector. For each period p, a run counts
// consecutive positions where seq[i] == seq[i - p]. A run of
// (minCopies - 1) * p positions spans minCopies copies of the p-mer unit;
// the run resets on mismatch so unrelated partial matches never combine.
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
                run = 0;
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
//   genic GC (45-55% in practice), including lower-GC organisms.
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

inline bool starship_supported_phylum(std::string_view phylum) {
    if (phylum.empty() || phylum == "." || phylum == "unknown" || phylum == "UNKNOWN")
        return true;  // no taxonomic context available; preserve sequence-only classifier.
    return phylum == "Ascomycota";
}

// ── detect_hgt_island ────────────────────────────────────────────────────
// Horizontal gene transfer island: GC content deviates from clade background
// by > ±gcDeviation over a sliding window of ≥ minLen bases.
//
// cladeGc: background GC of the clade (computed externally or passed as 0.45).
inline bool detect_hgt_island(std::string_view seq,
                               double cladeGc     = 0.45,
                               double gcDeviation = 0.10,  // candidate screen, not proof of HGT
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
// RIP (Repeat-Induced Point mutation): lightweight RIP-product proxy in a
// sliding window.  RIP leaves a C->T / G->A signature in repeated DNA, commonly
// summarized with dinucleotide indices such as TpA/CpA or composite RIP
// indices.  This detector intentionally stays portable and alignment-free, but
// it uses the expected CpA depletion / TpA enrichment signal rather than
// generic C/G imbalance, which can mark merely C-rich sequence as RIP.
inline bool detect_rip_window(std::string_view seq,
                               double productIndexThresh = 1.5,
                               int    winLen        = 500) {
    const int n = static_cast<int>(seq.size());
    if (n < winLen) return false;

    for (int i = 0; i + winLen <= n; i += winLen / 2) {
        int cpa = 0, tpa = 0, tpg = 0, apc = 0, gpt = 0, apt = 0;
        for (int j = i; j + 1 < i + winLen; ++j) {
            const char a = static_cast<char>(std::toupper(static_cast<unsigned char>(seq[static_cast<size_t>(j)])));
            const char b = static_cast<char>(std::toupper(static_cast<unsigned char>(seq[static_cast<size_t>(j + 1)])));
            if (a == 'C' && b == 'A') ++cpa;
            else if (a == 'T' && b == 'A') ++tpa;
            else if (a == 'T' && b == 'G') ++tpg;
            else if (a == 'A' && b == 'C') ++apc;
            else if (a == 'G' && b == 'T') ++gpt;
            else if (a == 'A' && b == 'T') ++apt;
        }

        const double productIndex =
            static_cast<double>(tpa + tpg + 1) / static_cast<double>(cpa + 1);
        const double compositeIndex =
            static_cast<double>(tpa + tpg + 1) / static_cast<double>(apc + gpt + apt + 1);
        if (productIndex >= productIndexThresh && compositeIndex >= 0.8 && tpa >= 20)
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
                                             double cladeGc = 0.45,
                                             std::string_view phylum = ".") {
    if (seq.size() < 50u) return ElementClass::NONE;

    // RIP takes priority (post-translational; affects any repeated element)
    if (detect_rip_window(seq))                        return ElementClass::RIP;

    // Starship: large AT-rich element with GC-rich cargo.  Treat this as a
    // taxon-aware label when phylum is known; without context use the
    // sequence-only classifier for standalone detector tests and callers that
    // do not pass taxonomy.
    if (starship_supported_phylum(phylum) &&
        detect_starship(seq, cladeGc))                 return ElementClass::STARSHIP;

    // HGT: GC-shifted island
    if (detect_hgt_island(seq, cladeGc))               return ElementClass::HGT;

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
// Accept inversions only when allele lengths differ within the relative tolerance.
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
    // Sentinel IDs stay in 1-30 (all < 32) so find_mems reliably stops
    // extension at every contig boundary and cannot form cross-contig MEMs.
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

            // Repetitive fungal sequence can yield thousands of equally long
            // SA hits for one query position. Emitting the whole interval
            // creates a MEM cloud that downstream chaining repeatedly mines
            // without adding breakpoint information. Keep a deterministic,
            // evenly spaced slice per query position; the global cap remains
            // as a final safety net across the whole query.
            constexpr int kMaxHitsPerQueryPos = 64;
            const int intervalN = hi2 - lo2 + 1;
            const int emitN = std::min(intervalN, kMaxHitsPerQueryPos);
            const int stride = std::max(1, intervalN / emitN);
            int emitted = 0;
            for (int k = lo2; k <= hi2 && emitted < emitN; k += stride, ++emitted) {
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
        // subtreeBest = max(best) over this node's whole subtree. Maintained
        // on insert/rotation so find_pred_score can prune entire subtrees
        // that cannot improve the running best — turning the predecessor
        // search from O(n) (full left-subtree walk) into ~O(log n) amortised
        // and the whole chain build from O(n^2) into ~O(n log n). On
        // chromosome-scale fungal query contigs the O(n^2) walk was the
        // dominant cost and made assembly-mode runs exceed the 4 h budget.
        float subtreeBest = 0.0f;
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
        nd.subtreeBest = nd.best;
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
    // Recompute subtreeBest for `idx` from its own best plus the children's
    // subtreeBest. Cheap O(1); call bottom-up after any structural change.
    void pull(int idx) {
        Node& n = nodes_[static_cast<size_t>(idx)];
        float b = n.best;
        if (n.left  >= 0) b = std::max(b, nodes_[static_cast<size_t>(n.left)].subtreeBest);
        if (n.right >= 0) b = std::max(b, nodes_[static_cast<size_t>(n.right)].subtreeBest);
        n.subtreeBest = b;
    }

    float find_pred_score(int node, int rPos, int qPos, int maxGap,
                          int& bestIdx) const {
        float best = 0.0f;
        bestIdx = -1;
        find_pred_rec(node, rPos, qPos, maxGap, best, bestIdx);
        return best;
    }

    // Exact predecessor search prunes:
    //   1. subtreeBest <= best: no entry in this subtree can beat the running best.
    //   2. nd.rPos < rPos - maxGap: every left-subtree rPos is out of band.
    void find_pred_rec(int node, int rPos, int qPos, int maxGap,
                       float& best, int& bestIdx) const {
        if (node < 0) return;
        const Node& nd = nodes_[static_cast<size_t>(node)];
        if (nd.subtreeBest <= best) return;
        // Current node qualifies as a predecessor only when both coordinates
        // are strictly smaller and within the chaining gap band.
        if (nd.rPos < rPos && nd.qPos < qPos &&
            rPos - nd.rPos <= maxGap && qPos - nd.qPos <= maxGap) {
            if (nd.best > best) { best = nd.best; bestIdx = node; }
        }
        // Left subtree (all rPos ≤ nd.rPos): worth visiting only if nd.rPos
        // itself is still within the gap band's lower bound — otherwise every
        // left node is too far back.
        if (nd.rPos >= rPos - maxGap)
            find_pred_rec(nd.left, rPos, qPos, maxGap, best, bestIdx);
        // Right subtree (all rPos ≥ nd.rPos): only when nd.rPos < rPos, since
        // nodes with rPos ≥ rPos can never be predecessors.
        if (nd.rPos < rPos)
            find_pred_rec(nd.right, rPos, qPos, maxGap, best, bestIdx);
    }

    int insert_node(int root, int idx) {
        if (root < 0) { pull(idx); return idx; }
        if (nodes_[static_cast<size_t>(idx)].rPos <=
            nodes_[static_cast<size_t>(root)].rPos) {
            int newLeft = insert_node(nodes_[static_cast<size_t>(root)].left, idx);
            nodes_[static_cast<size_t>(root)].left = newLeft;
            if (nodes_[static_cast<size_t>(newLeft)].prio >
                nodes_[static_cast<size_t>(root)].prio) {
                int nr = newLeft;
                nodes_[static_cast<size_t>(root)].left =
                    nodes_[static_cast<size_t>(nr)].right;
                nodes_[static_cast<size_t>(nr)].right = root;
                pull(root);
                pull(nr);
                return nr;
            }
            pull(root);
        } else {
            int newRight = insert_node(nodes_[static_cast<size_t>(root)].right, idx);
            nodes_[static_cast<size_t>(root)].right = newRight;
            if (nodes_[static_cast<size_t>(newRight)].prio >
                nodes_[static_cast<size_t>(root)].prio) {
                int nr = newRight;
                nodes_[static_cast<size_t>(root)].right =
                    nodes_[static_cast<size_t>(nr)].left;
                nodes_[static_cast<size_t>(nr)].left = root;
                pull(root);
                pull(nr);
                return nr;
            }
            pull(root);
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
        // contig. rBreakStart/rBreakEnd are returned in local contig coordinates
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

        // The query span the chain physically covers. A genuine tandem
        // duplication — its extra copy lives inside the query contig — cannot
        // be larger than this. Cumulative rGap, by contrast, sums EVERY noisy
        // backward MEM pair across the chain; on repeat-rich query genomes
        // (arbuscular mycorrhizal fungi) that aggregate explodes to many× the
        // contig length and produced 100–300 kb "DUP" calls on 20 kb contigs.
        const int qChainSpan = chain[static_cast<size_t>(N-1)].qPos +
                               chain[static_cast<size_t>(N-1)].len - chain[0].qPos;

        // Coherence guard for the per-pair DUP signal: a genuine tandem
        // duplication RE-TRAVERSES the duplicated window — after the backward
        // rPos jump the chain continues forward THROUGH the same reference
        // region. A lone repeat MEM, by contrast, jumps back once and is
        // immediately followed by a large compensating FORWARD leap back onto
        // the main diagonal. Require the successor MEM to stay inside the
        // duplicated window (small forward step) rather than leap past it;
        // an unverifiable backward jump at the very end of the chain is
        // treated as non-coherent (precision-favouring).
        bool perPairDupCoherent = false;
        if (dupPairIdx >= 0 && dupPairIdx + 2 < N) {
            const int afterRGap = chain[static_cast<size_t>(dupPairIdx + 2)].rPos -
                (chain[static_cast<size_t>(dupPairIdx + 1)].rPos +
                 chain[static_cast<size_t>(dupPairIdx + 1)].len);
            perPairDupCoherent = (afterRGap * 2 < -dupPairBack);
        }

        // DUP: cumulative reference overlap, OR a single backward rPos jump
        // ≥ minSvLen accompanied by a non-negative query gap (tandem-dup
        // signature: q advances by ~svLen while r returns to the original
        // copy's coordinates). Both forms are bounded by qChainSpan so a
        // noise-inflated rGap can never masquerade as a chromosome-scale DUP.
        const bool perPairDup    = (dupPairIdx >= 0 &&
                                    dupPairBack < -minSvLen &&
                                    dupPairQGap >= -minSvLen &&
                                    -dupPairBack >= dupPairQGap &&
                                    -dupPairBack <= qChainSpan &&
                                    perPairDupCoherent);
        const bool cumulativeDup = (rGap < -minSvLen) && (-rGap <= qChainSpan);
        if (cumulativeDup || perPairDup) {
            const int ctgOff = (c0 > 0)
                ? sa.contigEnd[static_cast<size_t>(c0) - 1] : 0;
            res.type = Type::DUP;
            // Prefer the dominant single backward jump: the per-pair anchor is
            // the real tandem-dup signature (one copy, one size), whereas
            // cumulative rGap is a noise-prone aggregate. Fall back to the
            // cumulative view only when no single dominant backward pair
            // exists (and only then within the qChainSpan bound above).
            if (perPairDup) {
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
            } else {
                res.qBreakStart = chain[0].qPos;
                res.qBreakEnd   = chain[static_cast<size_t>(N-1)].qPos +
                                  chain[static_cast<size_t>(N-1)].len - 1;
                res.rBreakStart = chain[0].rPos - ctgOff;
                res.rBreakEnd   = res.rBreakStart + (-rGap);
                res.svLen       = -rGap;
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
        // An insertion's novel sequence is physically present in the query
        // contig, so |svLen| cannot exceed the query span the chain covers.
        // A larger value means `delta` is an artifact of many small gaps
        // summed across a noisy chain, not a single event — reject it rather
        // than emit a chromosome-scale phantom INS.
        if (delta > 0 && res.svLen > qChainSpan) return Result{};
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

    // classify_all: emit ALL local INS/DEL gaps in the chain rather than only the
    // dominant one. The original classify() returns a single Result, which on real
    // diverged fungal genomes collapses 50–500 small SVs per query down to ~1–17
    // calls (assembly-mode F1 0.007–0.04 on compact_yeast). This walker keeps the
    // global TRA/INV/DUP semantics — those still return a single event — but for
    // INS/DEL it scans every consecutive MEM pair and emits an event whenever the
    // per-pair (qGap - rGap) difference is ≥ minSvLen, ignoring pairs that cross
    // a contig boundary (TRA territory), flip strand (INV territory), or jump
    // backward in r (DUP territory). All event coordinates are returned in local
    // contig space, same convention as classify().
    static std::vector<Result> classify_all(const std::vector<SuffixArray::Mem>& chain,
                                            const std::vector<bool>& isRevComp,
                                            const SuffixArray& sa,
                                            int minSvLen = 40) {
        std::vector<Result> events;
        if (chain.size() < 2) return events;

        Result global = classify(chain, isRevComp, sa, minSvLen);
        if (global.type == Type::TRA ||
            global.type == Type::INV ||
            global.type == Type::DUP) {
            events.push_back(global);
            return events;
        }

        auto contig_of = [&](int rPos) -> int {
            for (int ci = 0; ci < static_cast<int>(sa.contigEnd.size()); ++ci)
                if (rPos < sa.contigEnd[static_cast<size_t>(ci)]) return ci;
            return -1;
        };

        const int N = static_cast<int>(chain.size());
        for (int i = 0; i + 1 < N; ++i) {
            // Strand transitions belong to the global INV path.
            if (i < static_cast<int>(isRevComp.size()) &&
                i + 1 < static_cast<int>(isRevComp.size()) &&
                isRevComp[static_cast<size_t>(i)] != isRevComp[static_cast<size_t>(i + 1)])
                continue;

            const int rEndA = chain[static_cast<size_t>(i)].rPos +
                              chain[static_cast<size_t>(i)].len;
            const int qEndA = chain[static_cast<size_t>(i)].qPos +
                              chain[static_cast<size_t>(i)].len;
            const int lqg   = chain[static_cast<size_t>(i + 1)].qPos - qEndA;
            const int lrg   = chain[static_cast<size_t>(i + 1)].rPos - rEndA;

            // Cross-contig pairs are TRA territory.
            const int ciA = contig_of(chain[static_cast<size_t>(i)].rPos);
            const int ciB = contig_of(chain[static_cast<size_t>(i + 1)].rPos);
            if (ciA < 0 || ciB < 0 || ciA != ciB) continue;

            // Backward rPos jumps are DUP territory.
            if (lrg < 0) continue;
            // Forward overlap in query (negative qGap) means adjacent MEM
            // footprints overlap, likely from a near-tandem repeat; skip.
            if (lqg < 0) continue;

            const int ld = lqg - lrg;
            if (std::abs(ld) < minSvLen) continue;

            Result r;
            const int breakRPos = rEndA;
            const int breakCi   = contig_of(breakRPos > 0 ? breakRPos - 1 : 0);
            const int ctgOff    = (breakCi > 0)
                ? sa.contigEnd[static_cast<size_t>(breakCi) - 1] : 0;
            r.qBreakStart = qEndA;
            r.rBreakStart = breakRPos - ctgOff;
            r.svLen       = std::abs(ld);
            if (ld > 0) {
                r.type        = Type::INS;
                r.qBreakEnd   = r.qBreakStart + ld;
                r.rBreakEnd   = r.rBreakStart;
            } else {
                r.type        = Type::DEL;
                r.qBreakEnd   = r.qBreakStart;
                r.rBreakEnd   = r.rBreakStart + (-ld);
            }
            if (breakCi >= 0) r.rContig = sa.contigName[static_cast<size_t>(breakCi)];
            events.push_back(std::move(r));
        }

        // If per-pair scan found nothing but the global dominant signal is INS/DEL,
        // fall back to the single-event view rather than dropping the chain.
        if (events.empty() &&
            (global.type == Type::INS || global.type == Type::DEL)) {
            events.push_back(global);
        }
        return events;
    }
};

// =========================================================================
// PangenomeBubbleSvCaller — bubble-walking SV caller over a CladeGraph.
//
// Walks every per-(asm,contig) path against a REF path (the longest path by
// total node-sequence length, alphabetical tiebreak). Each ALT node is
// either a REF-anchor (its id appears in the REF node sequence) or a
// divergent node. A bubble closes whenever:
//   (a) a REF-anchor is hit with a non-empty divergent buffer, OR
//   (b) the REF-anchor's position is not exactly lastRefIdx + 1
//       (skipped REF span → DEL signature; backward jump → DUP signature).
//
// Classification (per bubble):
//   INV  : refLen == altLen > 0 and revcomp(altSeq) == refSeq
//   DUP  : altSeq begins with two or more consecutive copies of refSeq
//   INS  : altLen > refLen (or refLen == 0)
//   DEL  : altLen < refLen (or altLen == 0)
//   else : isComplex = true (same length, different sequence)
//
// Identical bubbles seen on multiple ALT paths are collapsed into one record
// with `supportingPaths` listing each contributing genome path.
//
// TRA: not emitted — CladeGraph paths here are per (asm, contig), so a
// genuine cross-contig translocation lives in a different path and isn't a
// single bubble. The MEM/chain-based path (SvTypeFromChain) covers TRA.
// =========================================================================
struct PangenomeBubbleSV {
    SvTypeFromChain::Type    type            = SvTypeFromChain::Type::NONE;
    std::string              cladeName;
    std::string              refPathName;
    int                      refNodeStart    = 0;     // REF anchor index (exclusive lower)
    int                      refNodeEnd      = 0;     // REF anchor index (exclusive upper)
    int                      refLenBp        = 0;
    int                      altLenBp        = 0;
    int                      svLen           = 0;     // signed altLen - refLen
    std::string              refSeq;
    std::string              altSeq;
    std::vector<std::string> supportingPaths;
    bool                     isComplex       = false;
    // TRA: at least one divergent node is also used on a contig other than
    // the ALT path's own contig — material has moved across reference
    // contigs. `traPartnerContigs` lists the foreign contigs implicated.
    std::vector<std::string> traPartnerContigs;
    // OFF_REF: ALT k-mer overlap with the union of REF-path sequences is
    // below 5% — bubble represents sequence largely absent from the clade
    // reference set (Path-C novelty candidate). The chain-based caller's
    // score_cross_clade_novelty() can subsequently qualify this as HGT /
    // NOVEL_WEAK / DIVERGED vs other clades.
    bool                     isOffRef        = false;
    double                   refKmerOverlap  = 1.0;
};

struct PangenomeBubbleSvCaller {
    int    minSvLen          = 40;
    int    maxSvLen          = 1'000'000;
    int    offRefKmerK       = 15;
    double offRefMaxOverlap  = 0.05;
    int    offRefMinAltBp    = 100;

    // CladeGraph path naming convention is "asm::contig" (see build_clade_graph).
    static std::pair<std::string, std::string>
    split_path_name(const std::string& name) {
        const auto pos = name.find("::");
        if (pos == std::string::npos) return {name, ""};
        return {name.substr(0, pos), name.substr(pos + 2)};
    }

    // Inline FNV-1a k-mer hashing for the OFF_REF overlap check. Keeping it
    // local avoids pulling fungi_tol_bridge.hpp into layer1 (which would
    // create a header cycle: bridge already includes layer3 → layer1).
    static void kmer_hashes_into(std::unordered_set<uint64_t>& out,
                                 const std::string& s, int k) {
        if (k <= 0 || static_cast<int>(s.size()) < k) return;
        const uint64_t basis = 14695981039346656037ULL;
        const uint64_t prime = 1099511628211ULL;
        for (size_t i = 0; i + static_cast<size_t>(k) <= s.size(); ++i) {
            uint64_t h = basis;
            for (int j = 0; j < k; ++j) {
                h ^= static_cast<uint64_t>(static_cast<unsigned char>(s[i + j]));
                h *= prime;
            }
            out.insert(h);
        }
    }

    static std::string revcomp(const std::string& s) {
        std::string out(s.size(), 'N');
        for (size_t i = 0; i < s.size(); ++i) {
            const char c = s[s.size() - 1 - i];
            switch (c) {
                case 'A': case 'a': out[i] = 'T'; break;
                case 'C': case 'c': out[i] = 'G'; break;
                case 'G': case 'g': out[i] = 'C'; break;
                case 'T': case 't': out[i] = 'A'; break;
                default:            out[i] = 'N'; break;
            }
        }
        return out;
    }

    // True iff `alt` begins with two or more consecutive copies of `ref`.
    static bool is_tandem_dup(const std::string& alt, const std::string& ref) {
        if (ref.empty() || alt.size() < 2 * ref.size()) return false;
        if (alt.compare(0, ref.size(), ref) != 0) return false;
        if (alt.compare(ref.size(), ref.size(), ref) != 0) return false;
        return true;
    }

    std::vector<PangenomeBubbleSV> call(const CladeGraph& g) const {
        std::vector<PangenomeBubbleSV> out;
        if (g.paths.empty() || g.nodes.empty()) return out;

        std::unordered_map<int, const std::string*> nodeSeq;
        nodeSeq.reserve(g.nodes.size() * 2);
        for (const auto& n : g.nodes) nodeSeq.emplace(n.id, &n.sequence);

        auto path_bp = [&](const CladeGraph::Path& p) -> size_t {
            size_t s = 0;
            for (int nid : p.nodes) {
                auto it = nodeSeq.find(nid);
                if (it != nodeSeq.end()) s += it->second->size();
            }
            return s;
        };
        const CladeGraph::Path* refPath = &g.paths.front();
        size_t refBp = path_bp(*refPath);
        for (const auto& p : g.paths) {
            const size_t b = path_bp(p);
            if (b > refBp || (b == refBp && p.name < refPath->name)) {
                refPath = &p;
                refBp   = b;
            }
        }

        std::unordered_map<int, std::vector<int>> refPos;
        refPos.reserve(refPath->nodes.size() * 2);
        for (size_t i = 0; i < refPath->nodes.size(); ++i) {
            refPos[refPath->nodes[i]].push_back(static_cast<int>(i));
        }

        // Per-node contig set (across all paths). A node whose set exceeds
        // {altContig} signals cross-contig material movement → TRA.
        std::unordered_map<int, std::unordered_set<std::string>> nodeContigs;
        nodeContigs.reserve(g.nodes.size() * 2);
        for (const auto& p : g.paths) {
            const std::string contig = split_path_name(p.name).second;
            for (int nid : p.nodes) nodeContigs[nid].insert(contig);
        }
        const std::string refContig = split_path_name(refPath->name).second;

        // REF k-mer union for the OFF_REF overlap check. Use REF path nodes
        // as the in-clade reference proxy (CladeGraph collapses identical
        // segments across genomes, so REF path nodes already approximate
        // the conserved core).
        std::unordered_set<uint64_t> refKmers;
        for (int nid : refPath->nodes) {
            auto it = nodeSeq.find(nid);
            if (it != nodeSeq.end())
                kmer_hashes_into(refKmers, *it->second, offRefKmerK);
        }

        auto concat_span = [&](int firstIdxExcl, int lastIdxExcl) -> std::string {
            std::string s;
            const int lo = firstIdxExcl + 1;
            const int hi = std::min(lastIdxExcl,
                                    static_cast<int>(refPath->nodes.size()));
            for (int i = std::max(lo, 0); i < hi; ++i) {
                auto it = nodeSeq.find(refPath->nodes[static_cast<size_t>(i)]);
                if (it != nodeSeq.end()) s += *it->second;
            }
            return s;
        };
        auto concat_nodes = [&](const std::vector<int>& ids) -> std::string {
            std::string s;
            for (int nid : ids) {
                auto it = nodeSeq.find(nid);
                if (it != nodeSeq.end()) s += *it->second;
            }
            return s;
        };

        struct BubbleKey {
            int         refStart;
            int         refEnd;
            std::string altSeq;
            bool operator==(const BubbleKey& o) const {
                return refStart == o.refStart && refEnd == o.refEnd &&
                       altSeq == o.altSeq;
            }
        };
        struct BubbleKeyHash {
            size_t operator()(const BubbleKey& k) const noexcept {
                size_t h = static_cast<size_t>(k.refStart) * 0x9E3779B185EBCA87ULL;
                h ^= static_cast<size_t>(k.refEnd) + 0x165667B19E3779F9ULL +
                     (h << 6) + (h >> 2);
                h ^= std::hash<std::string>{}(k.altSeq);
                return h;
            }
        };
        std::unordered_map<BubbleKey, size_t, BubbleKeyHash> seen;

        auto emit_bubble = [&](int refStartIdx, int refEndIdx,
                               const std::vector<int>& altDivergent,
                               const std::string& altPathName) {
            std::string refSeq = concat_span(refStartIdx, refEndIdx);
            std::string altSeq = concat_nodes(altDivergent);
            const int refLen = static_cast<int>(refSeq.size());
            const int altLen = static_cast<int>(altSeq.size());
            if (refLen == 0 && altLen == 0) return;
            const int absLen = std::max(refLen, altLen);
            if (absLen < minSvLen || absLen > maxSvLen) return;

            using T = SvTypeFromChain::Type;
            PangenomeBubbleSV sv;
            sv.cladeName    = g.cladeName;
            sv.refPathName  = refPath->name;
            sv.refNodeStart = refStartIdx;
            sv.refNodeEnd   = refEndIdx;
            sv.refLenBp     = refLen;
            sv.altLenBp     = altLen;
            sv.svLen        = altLen - refLen;
            sv.refSeq       = refSeq;
            sv.altSeq       = altSeq;
            sv.supportingPaths.push_back(altPathName);

            // Cross-contig TRA detection: union of contigs that host any
            // divergent node, minus the alt path's own contig. Any non-empty
            // remainder means the bubble material was moved across REF
            // contigs in this genome.
            const std::string altContig = split_path_name(altPathName).second;
            std::unordered_set<std::string> partnerSet;
            for (int nid : altDivergent) {
                auto it = nodeContigs.find(nid);
                if (it == nodeContigs.end()) continue;
                for (const auto& c : it->second)
                    if (c != altContig) partnerSet.insert(c);
            }
            sv.traPartnerContigs.assign(partnerSet.begin(), partnerSet.end());
            std::sort(sv.traPartnerContigs.begin(), sv.traPartnerContigs.end());
            const bool isCrossContig = !sv.traPartnerContigs.empty();

            // OFF_REF: ALT k-mer overlap with REF path's k-mer union.
            sv.refKmerOverlap = 1.0;
            if (altLen >= offRefMinAltBp) {
                std::unordered_set<uint64_t> altKmers;
                kmer_hashes_into(altKmers, altSeq, offRefKmerK);
                if (!altKmers.empty()) {
                    size_t inter = 0;
                    for (uint64_t h : altKmers) if (refKmers.count(h)) ++inter;
                    sv.refKmerOverlap =
                        static_cast<double>(inter) /
                        static_cast<double>(altKmers.size());
                    if (sv.refKmerOverlap < offRefMaxOverlap) sv.isOffRef = true;
                }
            }

            // Classification priority:
            //   TRA  (cross-contig signal — overrides others)
            //   INV  (length-matched reverse-complement allele)
            //   DUP  (tandem expansion of REF span)
            //   INS / DEL  (length-delta sign)
            //   COMPLEX  (equal length, non-INV substitution)
            if (isCrossContig) {
                sv.type = T::TRA;
            } else if (refLen > 0 && altLen > 0 && refLen == altLen &&
                       revcomp(altSeq) == refSeq) {
                sv.type = T::INV;
            } else if (is_tandem_dup(altSeq, refSeq)) {
                sv.type = T::DUP;
            } else if (refLen == 0 && altLen > 0) {
                sv.type = T::INS;
            } else if (altLen == 0 && refLen > 0) {
                sv.type = T::DEL;
            } else if (altLen > refLen) {
                sv.type = T::INS;
            } else if (altLen < refLen) {
                sv.type = T::DEL;
            } else {
                sv.type      = T::NONE;
                sv.isComplex = true;
            }

            BubbleKey key{refStartIdx, refEndIdx, altSeq};
            auto it = seen.find(key);
            if (it != seen.end()) {
                out[it->second].supportingPaths.push_back(altPathName);
            } else {
                seen.emplace(std::move(key), out.size());
                out.push_back(std::move(sv));
            }
        };

        for (const auto& p : g.paths) {
            if (&p == refPath) continue;
            int lastRefIdx = -1;
            std::vector<int> bubble;
            for (int nid : p.nodes) {
                auto it = refPos.find(nid);
                if (it == refPos.end()) {
                    bubble.push_back(nid);
                    continue;
                }
                // Snap to the REF occurrence closest to the monotone-next slot.
                int chosen  = it->second.front();
                int bestGap = std::abs(chosen - (lastRefIdx + 1));
                for (int candidate : it->second) {
                    const int gap = std::abs(candidate - (lastRefIdx + 1));
                    if (gap < bestGap) { bestGap = gap; chosen = candidate; }
                }
                const int  newRefIdx    = chosen;
                const bool isSkipForward = (newRefIdx > lastRefIdx + 1);
                const bool isBackward    = (newRefIdx <= lastRefIdx);
                const bool haveBubble    = !bubble.empty();

                if (haveBubble || isSkipForward || isBackward) {
                    const int rStart = lastRefIdx;
                    const int rEnd   = isBackward ? lastRefIdx + 1 : newRefIdx;
                    emit_bubble(rStart, rEnd, bubble, p.name);
                    bubble.clear();
                }
                lastRefIdx = newRefIdx;
            }
            if (!bubble.empty()) {
                emit_bubble(lastRefIdx,
                            static_cast<int>(refPath->nodes.size()),
                            bubble, p.name);
            }
        }

        return out;
    }
};

} // namespace tol

#endif // FUNGI_TOL_LAYER1_CLADE_GRAPH_HPP
