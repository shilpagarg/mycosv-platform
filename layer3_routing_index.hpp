#pragma once
// layer3_routing_index.hpp — v14
// Phylum-Sharded Routing Index  (Layer 3)
//
// DS-4   VP-tree nearest-clade routing  (Uhlmann 1991 / Yianilos 1993)
// DS-5   Bloom filter hash-suppression prefilter  (Bloom 1970)
// DS-6   DSU sub-clade splitter  (Tarjan 1975)
// DS-19  Skip-list style sparse disk directory  (Pugh 1989) for million-scale
//        centroid shortlist routing without materializing the full catalog
//        in RAM.
// DS-14  WaveletTree over BWT for O(k log σ) alphabet prefilter
//
// FIX-LOCK v14: route() previously acquired registryMu_ TWICE per phylum
// (once to snapshot the phylum list, once for the per-shard pointer).
// The double-acquire is now eliminated: the phylum-shard snapshot and
// pointer capture are merged into a single brief critical section that
// produces a vector<PhylumShard*>; all per-shard queries then run under
// only the per-shard shared_lock, with no global lock held.

#include "layer1_clade_graph.hpp"

#include <algorithm>
#include <atomic>
#include <bitset>
#include <cmath>
#include <deque>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <memory>
#include <mutex>
#include <random>
#include <shared_mutex>
#include <string_view>
#include <unordered_map>
#include <utility>

namespace fs = std::filesystem;

namespace tol {

// =========================================================================
// DS-5: Bloom filter  (64 KB, k=7 hash functions, ~1% FPR @ 56 K items)
// =========================================================================
struct BloomFilter {
    static constexpr size_t kBits  = 524288; // 64 KB
    static constexpr int    kFuncs = 7;

    std::vector<uint64_t> words;
    BloomFilter() : words(kBits / 64, 0) {}

    void insert(uint64_t h) {
        for (int i = 0; i < kFuncs; ++i) {
            uint64_t x = hash_i(h, i) % kBits;
            words[x / 64] |= (1ULL << (x % 64));
        }
    }

    bool probably_contains(uint64_t h) const {
        for (int i = 0; i < kFuncs; ++i) {
            uint64_t x = hash_i(h, i) % kBits;
            if (!(words[x / 64] & (1ULL << (x % 64)))) return false;
        }
        return true;
    }

    bool empty() const {
        for (auto w : words) if (w) return false;
        return true;
    }

    void write(std::ostream& o) const {
        uint32_t n = static_cast<uint32_t>(words.size());
        o.write(reinterpret_cast<const char*>(&n), 4);
        o.write(reinterpret_cast<const char*>(words.data()),
                static_cast<std::streamsize>(n * 8));
    }

    static BloomFilter read(std::istream& in) {
        BloomFilter bf;
        uint32_t n = 0;
        in.read(reinterpret_cast<char*>(&n), 4);
        if (n != kBits / 64) throw std::runtime_error("Bloom filter size mismatch");
        in.read(reinterpret_cast<char*>(bf.words.data()),
                static_cast<std::streamsize>(n * 8));
        return bf;
    }

private:
    static uint64_t hash_i(uint64_t h, int i) {
        // Kirsch & Mitzenmacher 2006 double-hashing
        uint64_t h1 = h ^ (h >> 33);
        h1 *= 0xff51afd7ed558ccdULL; h1 ^= h1 >> 33;
        h1 *= 0xc4ceb9fe1a85ec53ULL; h1 ^= h1 >> 33;
        uint64_t h2 = h ^ (h >> 17);
        h2 *= 0xbf58476d1ce4e5b9ULL; h2 ^= h2 >> 31;
        return h1 + static_cast<uint64_t>(i) * h2;
    }
};

struct CompressedBucketBitmap {
    static constexpr uint16_t kBucketBits = 12;
    static constexpr uint16_t kBucketMask = (1u << kBucketBits) - 1u;

    std::vector<uint16_t> buckets;

    void build_from_hashes(const std::vector<uint64_t>& hashes) {
        buckets.clear();
        buckets.reserve(hashes.size());
        for (uint64_t h : hashes) {
            const uint16_t bucket = static_cast<uint16_t>((h >> (64 - kBucketBits)) & kBucketMask);
            buckets.push_back(bucket);
        }
        std::sort(buckets.begin(), buckets.end());
        buckets.erase(std::unique(buckets.begin(), buckets.end()), buckets.end());
    }

    bool intersects(const CompressedBucketBitmap& o) const {
        size_t ia = 0, ib = 0;
        while (ia < buckets.size() && ib < o.buckets.size()) {
            if (buckets[ia] < o.buckets[ib]) ++ia;
            else if (buckets[ia] > o.buckets[ib]) ++ib;
            else return true;
        }
        return false;
    }

    bool empty() const { return buckets.empty(); }
};

// =========================================================================
// DS-14: WaveletTree alphabet pre-filter over DNA text
// σ=4 (A=0,C=1,G=2,T=3); rank(c,i) in O(1) via prefix popcount table.
// =========================================================================
struct WaveletTree {
    static constexpr int SIGMA     = 4;
    static constexpr int LOG_SIGMA = 2;

    static int encode(char c) {
        switch (c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default:             return 0;
        }
    }

    int n_ = 0;
    std::array<std::vector<uint64_t>, SIGMA> bits_;
    std::array<std::vector<int>,      SIGMA> cnt_;

    void build(const std::string& text) {
        n_ = static_cast<int>(text.size());
        if (n_ == 0) return;
        const int words = (n_ + 63) / 64;
        for (int c = 0; c < SIGMA; ++c) {
            bits_[static_cast<size_t>(c)].assign(static_cast<size_t>(words), 0ULL);
            cnt_[static_cast<size_t>(c)].assign(static_cast<size_t>(words + 1), 0);
        }
        for (int i = 0; i < n_; ++i) {
            int code = encode(text[static_cast<size_t>(i)]);
            bits_[static_cast<size_t>(code)][static_cast<size_t>(i / 64)]
                |= (1ULL << (i % 64));
        }
        for (int c = 0; c < SIGMA; ++c) {
            cnt_[static_cast<size_t>(c)][0] = 0;
            for (int w = 0; w < words; ++w) {
                cnt_[static_cast<size_t>(c)][static_cast<size_t>(w + 1)] =
                    cnt_[static_cast<size_t>(c)][static_cast<size_t>(w)] +
                    __builtin_popcountll(bits_[static_cast<size_t>(c)][static_cast<size_t>(w)]);
            }
        }
    }

    // O(1) rank via prefix popcount + partial word popcount.
    int rank(int code, int i) const {
        if (n_ == 0 || i <= 0 || code < 0 || code >= SIGMA) return 0;
        i = std::min(i, n_);
        const int fw  = i / 64;
        const int rem = i % 64;
        int result = cnt_[static_cast<size_t>(code)][static_cast<size_t>(fw)];
        if (rem > 0) {
            uint64_t mask = (1ULL << rem) - 1ULL;
            result += __builtin_popcountll(
                bits_[static_cast<size_t>(code)][static_cast<size_t>(fw)] & mask);
        }
        return result;
    }

    // Zero-cost alphabet prefilter: return 0 if any character in pattern
    // never appears in text (impossible k-mer match).
    int count_kmer(const std::string& pattern) const {
        if (n_ == 0 || pattern.empty()) return 0;
        for (char ch : pattern)
            if (rank(encode(ch), n_) == 0) return 0;
        return 1;
    }

    bool has_char(char c) const { return n_ > 0 && rank(encode(c), n_) > 0; }
    bool empty()          const { return n_ == 0; }
};

// =========================================================================
// CladeCentroid — FracMin sketch + Bloom filter + WaveletTree prefilter
// =========================================================================
struct CladeCentroid {
    std::string            cladeName;
    std::string            phylum;
    std::string            cladeRank;
    std::vector<uint64_t>  centroidHashes;
    size_t                 genomeSyncmers = 0;
    BloomFilter            suppressFilter;
    WaveletTree            seqFilter;
    CompressedBucketBitmap membershipSketch;

    void build_prefilters() {
        suppressFilter = BloomFilter();
        for (uint64_t h : centroidHashes) suppressFilter.insert(h);
        membershipSketch.build_from_hashes(centroidHashes);
    }

    bool may_share_hash_with(const CladeCentroid& o) const {
        if (!membershipSketch.empty() && !o.membershipSketch.empty() &&
            !membershipSketch.intersects(o.membershipSketch))
            return false;

        const auto* probe = &centroidHashes;
        const auto* filt  = &o.suppressFilter;
        if (probe->size() > o.centroidHashes.size()) {
            probe = &o.centroidHashes;
            filt = &suppressFilter;
        }
        if (probe->empty() || filt->empty()) return true;
        for (uint64_t h : *probe)
            if (filt->probably_contains(h))
                return true;
        return false;
    }

    double jaccard_distance(const CladeCentroid& o) const {
        if (centroidHashes.empty() && o.centroidHashes.empty()) return 0.0;
        if (centroidHashes.empty() || o.centroidHashes.empty()) return 1.0;
        if (!may_share_hash_with(o)) return 1.0;
        size_t inter = 0, ia = 0, ib = 0;
        const auto& a = centroidHashes;
        const auto& b = o.centroidHashes;
        while (ia < a.size() && ib < b.size()) {
            if      (a[ia] < b[ib]) ++ia;
            else if (a[ia] > b[ib]) ++ib;
            else { ++inter; ++ia; ++ib; }
        }
        const size_t uni = a.size() + b.size() - inter;
        return uni > 0 ? 1.0 - static_cast<double>(inter) / static_cast<double>(uni) : 0.0;
    }

    struct StreamBuilder {
        std::string name, phylum, cladeRank;
        size_t maxH = 4096;
        std::vector<uint64_t> pool;
        size_t rawSyncmerCount = 0;

        void accumulate(const std::vector<uint64_t>& hashes) {
            pool.insert(pool.end(), hashes.begin(), hashes.end());
            rawSyncmerCount += hashes.size();
        }

        CladeCentroid finalize() {
            std::sort(pool.begin(), pool.end());
            pool.erase(std::unique(pool.begin(), pool.end()), pool.end());
            if (pool.size() > maxH) pool.resize(maxH);
            CladeCentroid c;
            c.cladeName       = name;
            c.phylum          = phylum;
            c.cladeRank       = cladeRank;
            c.centroidHashes  = std::move(pool);
            c.genomeSyncmers  = rawSyncmerCount;
            c.build_prefilters();
            return c;
        }
    };
};

inline CladeCentroid make_query_centroid_for_routing(std::string_view seq,
                                                     const SyncmerParams& sp,
                                                     double densityThresh) {
    const uint64_t thresh =
        static_cast<uint64_t>(densityThresh * std::ldexp(1.0, 64) - 1.0);
    CladeCentroid qc;
    qc.cladeName = "query";
    auto smers = BaseBlockSegmenter::syncmers(seq, sp.k, sp.s);
    for (const auto& [pos, h] : smers)
        if (h && h <= thresh) qc.centroidHashes.push_back(h);
    std::sort(qc.centroidHashes.begin(), qc.centroidHashes.end());
    qc.centroidHashes.erase(
        std::unique(qc.centroidHashes.begin(), qc.centroidHashes.end()),
        qc.centroidHashes.end());
    qc.build_prefilters();
    return qc;
}

// =========================================================================
// FracMinSketch — per-genome sketch for bootstrap clustering
// =========================================================================
struct FracMinSketch {
    std::string            asmName, cladeName;
    std::vector<uint64_t>  hashes;
    size_t                 genomeSyncmers = 0;
};

// =========================================================================
// DS-6: DSU  (Tarjan 1975) — path compression + union by rank
// =========================================================================
struct DSU {
    // FIX: use int throughout and cast only at the vector subscript boundary
    // to silence -Wsign-conversion without changing the algorithm or ABI.
    // The public API remains int-typed (callers pass genome indices as int).
    std::vector<int> parent, rank_;

    explicit DSU(int n)
        : parent(static_cast<size_t>(n)), rank_(static_cast<size_t>(n), 0) {
        std::iota(parent.begin(), parent.end(), 0);
    }

    int find(int x) {
        while (parent[static_cast<size_t>(x)] != x) {
            parent[static_cast<size_t>(x)] =
                parent[static_cast<size_t>(parent[static_cast<size_t>(x)])]; // path halving
            x = parent[static_cast<size_t>(x)];
        }
        return x;
    }

    void unite(int a, int b) {
        a = find(a); b = find(b);
        if (a == b) return;
        if (rank_[static_cast<size_t>(a)] < rank_[static_cast<size_t>(b)])
            std::swap(a, b);
        parent[static_cast<size_t>(b)] = a;
        if (rank_[static_cast<size_t>(a)] == rank_[static_cast<size_t>(b)])
            ++rank_[static_cast<size_t>(a)];
    }

    bool same(int a, int b) { return find(a) == find(b); }
};

// DS-6: O(N log N) sub-clade splitter via random projection + DSU.
inline std::vector<std::vector<size_t>>
split_into_subclades(const std::vector<FracMinSketch>& sketches,
                      size_t maxPerGroup,
                      uint64_t seed = 0xdeadbeefcafeULL) {
    const int N = static_cast<int>(sketches.size());
    if (N == 0) return {};
    if (N <= static_cast<int>(maxPerGroup)) {
        std::vector<size_t> g;
        g.reserve(static_cast<size_t>(N));
        for (int i = 0; i < N; ++i) g.push_back(static_cast<size_t>(i));
        return { std::move(g) };
    }

    static constexpr int kDims = 32;
    std::mt19937_64 rng(seed);

    // Union of all hashes — used to sample pivot positions
    std::vector<uint64_t> allH;
    allH.reserve(sketches.size() * 64);
    for (const auto& s : sketches)
        allH.insert(allH.end(), s.hashes.begin(), s.hashes.end());
    std::sort(allH.begin(), allH.end());
    allH.erase(std::unique(allH.begin(), allH.end()), allH.end());

    if (allH.empty()) {
        // No hashes — distribute round-robin
        std::vector<std::vector<size_t>> groups;
        for (int i = 0; i < N; ++i) {
            if (groups.empty() || groups.back().size() >= maxPerGroup)
                groups.emplace_back();
            groups.back().push_back(static_cast<size_t>(i));
        }
        return groups;
    }

    std::uniform_int_distribution<size_t> pick(0, allH.size() - 1);
    std::vector<uint64_t> pivots(kDims);
    for (int d = 0; d < kDims; ++d) pivots[static_cast<size_t>(d)] = allH[pick(rng)];
    std::sort(pivots.begin(), pivots.end());

    // Projection key per genome
    std::vector<uint32_t> keys(static_cast<size_t>(N), 0);
    for (int i = 0; i < N; ++i)
        for (int d = 0; d < kDims; ++d)
            if (std::binary_search(sketches[static_cast<size_t>(i)].hashes.begin(),
                                   sketches[static_cast<size_t>(i)].hashes.end(),
                                   pivots[static_cast<size_t>(d)]))
                keys[static_cast<size_t>(i)] |= (1u << d);

    std::vector<int> order(static_cast<size_t>(N));
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return keys[static_cast<size_t>(a)] < keys[static_cast<size_t>(b)]; });

    DSU dsu(N);
    static constexpr int kHammingThr = 4;
    for (int i = 0; i + 1 < N; ++i) {
        const int oi  = order[static_cast<size_t>(i)];
        const int oi1 = order[static_cast<size_t>(i + 1)];
        int diff = static_cast<int>(__builtin_popcount(
            keys[static_cast<size_t>(oi)] ^
            keys[static_cast<size_t>(oi1)]));
        if (diff <= kHammingThr)
            dsu.unite(oi, oi1);
    }

    // Collect groups from DSU roots
    std::unordered_map<int, std::vector<size_t>> byRoot;
    for (int i = 0; i < N; ++i)
        byRoot[dsu.find(i)].push_back(static_cast<size_t>(i));

    // Split oversized groups recursively (one level is enough in practice)
    std::vector<std::vector<size_t>> result;
    for (auto& [root, grp] : byRoot) {
        if (grp.size() <= maxPerGroup) {
            result.push_back(std::move(grp));
        } else {
            // Deterministic round-robin split
            size_t chunks = (grp.size() + maxPerGroup - 1) / maxPerGroup;
            for (size_t c = 0; c < chunks; ++c) {
                std::vector<size_t> sub;
                for (size_t j = c; j < grp.size(); j += chunks)
                    sub.push_back(grp[j]);
                result.push_back(std::move(sub));
            }
        }
    }
    return result;
}

// =========================================================================
// DS-4: VP-Tree  (Uhlmann 1991 / Yianilos 1993)
// =========================================================================
struct VPTree {
    struct Node { int idx = -1; double mu = 0.0; int left = -1, right = -1; };

    std::vector<CladeCentroid> centroids_;
    std::vector<Node>          nodes_;

    void build(const std::vector<CladeCentroid>& cs) {
        centroids_ = cs;
        for (auto& c : centroids_) c.build_prefilters();
        nodes_.clear();
        if (centroids_.empty()) return;
        std::vector<int> inds(centroids_.size());
        std::iota(inds.begin(), inds.end(), 0);
        build_node(inds, 0, static_cast<int>(inds.size()));
    }

    struct RouteResult {
        std::string cladeName;
        std::string phylum;
        double      jaccard = 0.0;
    };

    std::vector<RouteResult> query_topk(const CladeCentroid& q, size_t k) const {
        if (nodes_.empty()) return {};
        std::vector<std::pair<double, int>> heap;
        heap.reserve(k + 1);
        double tau = std::numeric_limits<double>::max();
        search(0, q, k, heap, tau);
        std::sort_heap(heap.begin(), heap.end());
        std::vector<RouteResult> res;
        res.reserve(heap.size());
        for (const auto& [dist, idx] : heap) {
            RouteResult r;
            r.cladeName = centroids_[static_cast<size_t>(idx)].cladeName;
            r.phylum    = centroids_[static_cast<size_t>(idx)].phylum;
            r.jaccard   = 1.0 - dist; // Jaccard similarity = 1 - Jaccard distance
            res.push_back(std::move(r));
        }
        return res;
    }

    bool   empty() const { return centroids_.empty(); }
    size_t size()  const { return centroids_.size(); }

    void save(const std::string& path) const {
        std::ofstream out(path, std::ios::binary);
        if (!out) return;
        uint32_t nc = static_cast<uint32_t>(centroids_.size());
        out.write(reinterpret_cast<const char*>(&nc), 4);
        for (const auto& c : centroids_) {
            uint32_t nl = static_cast<uint32_t>(c.cladeName.size());
            out.write(reinterpret_cast<const char*>(&nl), 4);
            out.write(c.cladeName.data(), static_cast<std::streamsize>(nl));
            uint32_t np = static_cast<uint32_t>(c.phylum.size());
            out.write(reinterpret_cast<const char*>(&np), 4);
            out.write(c.phylum.data(), static_cast<std::streamsize>(np));
            uint32_t nh = static_cast<uint32_t>(c.centroidHashes.size());
            out.write(reinterpret_cast<const char*>(&nh), 4);
            out.write(reinterpret_cast<const char*>(c.centroidHashes.data()),
                      static_cast<std::streamsize>(nh * 8));
        }
        uint32_t nn = static_cast<uint32_t>(nodes_.size());
        out.write(reinterpret_cast<const char*>(&nn), 4);
        for (const auto& n : nodes_) {
            out.write(reinterpret_cast<const char*>(&n.idx),   4);
            out.write(reinterpret_cast<const char*>(&n.mu),    8);
            out.write(reinterpret_cast<const char*>(&n.left),  4);
            out.write(reinterpret_cast<const char*>(&n.right), 4);
        }
    }

    static VPTree load(const std::string& path) {
        VPTree t;
        std::ifstream in(path, std::ios::binary);
        if (!in) return t;
        uint32_t nc = 0;
        in.read(reinterpret_cast<char*>(&nc), 4);
        t.centroids_.resize(nc);
        for (auto& c : t.centroids_) {
            uint32_t nl = 0; in.read(reinterpret_cast<char*>(&nl), 4);
            c.cladeName.resize(nl); in.read(c.cladeName.data(), static_cast<std::streamsize>(nl));
            uint32_t np = 0; in.read(reinterpret_cast<char*>(&np), 4);
            c.phylum.resize(np); in.read(c.phylum.data(), static_cast<std::streamsize>(np));
            uint32_t nh = 0; in.read(reinterpret_cast<char*>(&nh), 4);
            c.centroidHashes.resize(nh);
            in.read(reinterpret_cast<char*>(c.centroidHashes.data()),
                    static_cast<std::streamsize>(nh * 8));
            c.build_prefilters();
        }
        uint32_t nn = 0; in.read(reinterpret_cast<char*>(&nn), 4);
        t.nodes_.resize(nn);
        for (auto& n : t.nodes_) {
            in.read(reinterpret_cast<char*>(&n.idx),   4);
            in.read(reinterpret_cast<char*>(&n.mu),    8);
            in.read(reinterpret_cast<char*>(&n.left),  4);
            in.read(reinterpret_cast<char*>(&n.right), 4);
        }
        return t;
    }

private:
    int build_node(std::vector<int>& inds, int lo, int hi) {
        if (lo >= hi) return -1;
        const int nodeIdx = static_cast<int>(nodes_.size());
        nodes_.push_back({});
        Node& node = nodes_.back();

        if (hi - lo == 1) {
            node.idx = inds[static_cast<size_t>(lo)];
            return nodeIdx;
        }

        std::swap(inds[static_cast<size_t>(lo)],
                  inds[static_cast<size_t>(lo + (hi - lo) / 2)]);
        node.idx = inds[static_cast<size_t>(lo)];
        const CladeCentroid& vp = centroids_[static_cast<size_t>(node.idx)];

        std::vector<std::pair<double, int>> dists;
        dists.reserve(static_cast<size_t>(hi - lo - 1));
        for (int i = lo + 1; i < hi; ++i) {
            const int idx_i = inds[static_cast<size_t>(i)];
            dists.push_back({ vp.jaccard_distance(centroids_[static_cast<size_t>(idx_i)]), idx_i });
        }

        const size_t mid = dists.size() / 2;
        std::nth_element(dists.begin(), dists.begin() + static_cast<ptrdiff_t>(mid), dists.end());
        node.mu = dists[mid].first;

        std::vector<int> left_inds, right_inds;
        for (size_t i = 0; i < dists.size(); ++i)
            (i <= mid ? left_inds : right_inds).push_back(dists[i].second);

        const int lLo = lo + 1;
        const int lHi = lLo + static_cast<int>(left_inds.size());
        const int rLo = lHi;
        const int rHi = rLo + static_cast<int>(right_inds.size());
        for (int i = 0; i < static_cast<int>(left_inds.size());  ++i)
            inds[static_cast<size_t>(lLo + i)] = left_inds[static_cast<size_t>(i)];
        for (int i = 0; i < static_cast<int>(right_inds.size()); ++i)
            inds[static_cast<size_t>(rLo + i)] = right_inds[static_cast<size_t>(i)];

        // node ref may be invalidated by push_back in recursion — use saved index
        int lChild = build_node(inds, lLo, lHi);
        int rChild = build_node(inds, rLo, rHi);
        nodes_[static_cast<size_t>(nodeIdx)].left  = lChild;
        nodes_[static_cast<size_t>(nodeIdx)].right = rChild;
        return nodeIdx;
    }

    void search(int nIdx, const CladeCentroid& q, size_t k,
                std::vector<std::pair<double, int>>& heap, double& tau) const {
        if (nIdx < 0 || nIdx >= static_cast<int>(nodes_.size())) return;
        const Node& node = nodes_[static_cast<size_t>(nIdx)];
        const double d   = centroids_[static_cast<size_t>(node.idx)].jaccard_distance(q);

        if (heap.size() < k || d < heap.front().first) {
            heap.push_back({ d, node.idx });
            std::push_heap(heap.begin(), heap.end());
            if (heap.size() > k) {
                std::pop_heap(heap.begin(), heap.end());
                heap.pop_back();
            }
            if (heap.size() == k) tau = heap.front().first;
        }

        if (node.left < 0 && node.right < 0) return;
        const bool goLeft  = (d - tau <= node.mu);
        const bool goRight = (d + tau >= node.mu);
        if (d <= node.mu) {
            if (goLeft)  search(node.left,  q, k, heap, tau);
            if (goRight) search(node.right, q, k, heap, tau);
        } else {
            if (goRight) search(node.right, q, k, heap, tau);
            if (goLeft)  search(node.left,  q, k, heap, tau);
        }
    }
};

// =========================================================================
// PhylumShardedRouter  — per-phylum VP-tree + Tier-B fallback
//
// FIX-LOCK v14:
//   route() now captures (PhylumShard*, bool) pairs under a single brief
//   registryMu_ hold, then releases registryMu_ before doing any per-shard
//   work.  This eliminates the double-acquire that existed in v12/v13.
// =========================================================================
class PhylumShardedRouter {
public:
    struct RouteResult {
        std::string cladeName;
        std::string phylum;
        double      jaccard = 0.0;
    };

    struct QueryWindowRouteSummary {
        size_t windowStart = 0;
        size_t windowEnd = 0;
        std::vector<RouteResult> routes;
    };

    void set_route_cache_capacity(size_t n) { routeCacheCapacity_ = std::max<size_t>(1, n); }

    void reset() {
        std::lock_guard<std::mutex> lk(registryMu_);
        dirty_.store(false, std::memory_order_release);
        pendingByPhylum_.clear();
        pendingFbByPhylum_.clear();
        shards_.clear();
        fbShards_.clear();
        {
            std::lock_guard<std::mutex> cacheLk(routeCacheMu_);
            routeCache_.clear();
            routeCacheOrder_.clear();
        }
    }

    std::vector<QueryWindowRouteSummary> route_windows(std::string_view seq,
                                                       const SyncmerParams& sp,
                                                       const SyncmerParams& fbSp,
                                                       double densityThresh,
                                                       size_t topK,
                                                       size_t windowBp,
                                                       size_t windowOverlap) {
        std::vector<QueryWindowRouteSummary> out;
        if (seq.empty()) return out;
        windowBp = std::max<size_t>(1, windowBp);
        windowOverlap = std::min(windowOverlap, windowBp - 1);
        const size_t step = std::max<size_t>(1, windowBp - windowOverlap);
        if (seq.size() <= windowBp) {
            QueryWindowRouteSummary q;
            q.windowStart = 0;
            q.windowEnd = seq.size();
            q.routes = route(seq, sp, fbSp, densityThresh, topK);
            out.push_back(std::move(q));
            return out;
        }
        for (size_t start = 0; start < seq.size(); start += step) {
            const size_t end = std::min(seq.size(), start + windowBp);
            QueryWindowRouteSummary q;
            q.windowStart = start;
            q.windowEnd = end;
            q.routes = route(seq.substr(start, end - start), sp, fbSp, densityThresh, topK);
            out.push_back(std::move(q));
            if (end == seq.size()) break;
        }
        return out;
    }

    void register_clade_centroid(const CladeCentroid& c) {
        std::lock_guard<std::mutex> lk(registryMu_);
        pendingByPhylum_[c.phylum].push_back(c);
        dirty_.store(true, std::memory_order_release);
    }

    void register_clade_centroids(const CladeCentroid& primary,
                                   const CladeCentroid& fallback) {
        std::lock_guard<std::mutex> lk(registryMu_);
        pendingByPhylum_[primary.phylum].push_back(primary);
        pendingFbByPhylum_[fallback.phylum].push_back(fallback);
        dirty_.store(true, std::memory_order_release);
    }

    void rebuild() {
        std::lock_guard<std::mutex> lk(registryMu_);
        rebuild_locked();
    }

    // FIX-LOCK v14: single registryMu_ acquisition for snapshot.
    std::vector<RouteResult> route(std::string_view seq,
                                    const SyncmerParams& sp,
                                    const SyncmerParams& fbSp,
                                    double densityThresh,
                                    size_t topK) {
        if (dirty_.load(std::memory_order_acquire)) {
            std::lock_guard<std::mutex> lk(registryMu_);
            if (dirty_.load(std::memory_order_relaxed)) rebuild_locked();
        }

        const CladeCentroid qc = make_query_centroid(seq, sp, densityThresh);
        bool haveFallbackQc = false;
        CladeCentroid fallbackQc;
        const uint64_t cacheKey = centroid_signature(qc, topK);
        {
            std::lock_guard<std::mutex> lk(routeCacheMu_);
            auto it = routeCache_.find(cacheKey);
            if (it != routeCache_.end()) return it->second;
        }

        // ── Snapshot shard pointers under ONE brief lock ─────────────────
        struct ShardSnapshot { PhylumShard* shard; bool isFallback; };
        std::vector<ShardSnapshot> snap;
        {
            std::lock_guard<std::mutex> lk(registryMu_);
            snap.reserve(shards_.size() + fbShards_.size());
            for (auto& [p, s] : shards_)   snap.push_back({ s.get(), false });
            for (auto& [p, s] : fbShards_)  snap.push_back({ s.get(), true  });
        }
        // ── Per-shard queries with only per-shard shared_lock ────────────
        std::vector<std::pair<double, RouteResult>> all, fbAll;
        for (auto& [shard, isFb] : snap) {
            if (!shard) continue;
            std::shared_lock<std::shared_mutex> sl(shard->mu);
            if (shard->tree.empty()) continue;
            if (isFb && !haveFallbackQc) {
                fallbackQc = make_query_centroid(seq, fbSp, densityThresh);
                haveFallbackQc = true;
            }
            auto results = shard->tree.query_topk(isFb ? fallbackQc : qc, topK);
            auto& dest = isFb ? fbAll : all;
            for (auto& r : results) {
                RouteResult rr;
                rr.cladeName = std::move(r.cladeName);
                rr.phylum    = std::move(r.phylum);
                rr.jaccard   = r.jaccard;
                dest.emplace_back(rr.jaccard, std::move(rr));
            }
        }

        // Sort and pick top-K from primary results
        std::sort(all.begin(), all.end(),
                  [](const auto& a, const auto& b){ return a.first > b.first; });
        std::vector<RouteResult> results;
        results.reserve(std::min(topK, all.size()));
        for (size_t i = 0; i < std::min(topK, all.size()); ++i)
            results.push_back(std::move(all[i].second));

        // Tier-B fallback if primary returned nothing
        if (results.empty() && !fbAll.empty()) {
            std::sort(fbAll.begin(), fbAll.end(),
                      [](const auto& a, const auto& b){ return a.first > b.first; });
            for (size_t i = 0; i < std::min(topK, fbAll.size()); ++i)
                results.push_back(std::move(fbAll[i].second));
        }
        {
            std::lock_guard<std::mutex> lk(routeCacheMu_);
            if (routeCache_.size() >= routeCacheCapacity_ && !routeCacheOrder_.empty()) {
                routeCache_.erase(routeCacheOrder_.front());
                routeCacheOrder_.pop_front();
            }
            routeCache_[cacheKey] = results;
            routeCacheOrder_.push_back(cacheKey);
        }
        return results;
    }

    size_t clade_count() const {
        std::lock_guard<std::mutex> lk(registryMu_);
        size_t n = 0;
        for (const auto& [p, v] : pendingByPhylum_) n += v.size();
        return n;
    }

    void save(const std::string& indexDir) const {
        fs::create_directories(indexDir);
        std::lock_guard<std::mutex> lk(registryMu_);
        for (const auto& [phy, shard] : shards_) {
            std::shared_lock<std::shared_mutex> sl(shard->mu);
            shard->tree.save(indexDir + "/vptree_" + sanitize(phy) + ".bin");
        }
        for (const auto& [phy, shard] : fbShards_) {
            std::shared_lock<std::shared_mutex> sl(shard->mu);
            shard->tree.save(indexDir + "/vptree_fb_" + sanitize(phy) + ".bin");
        }
    }

    void load(const std::string& indexDir) {
        if (!fs::exists(indexDir)) return;
        std::lock_guard<std::mutex> lk(registryMu_);
        for (const auto& entry : fs::directory_iterator(indexDir)) {
            const auto fn = entry.path().filename().string();
            if (fn.rfind("vptree_fb_", 0) == 0 && fn.size() > 14) {
                const std::string phy = fn.substr(10, fn.size() - 14);
                auto& shard = fbShards_[phy];
                if (!shard) shard = std::make_unique<PhylumShard>();
                std::unique_lock<std::shared_mutex> ul(shard->mu);
                shard->tree = VPTree::load(entry.path().string());
            } else if (fn.rfind("vptree_", 0) == 0 && fn.size() > 11) {
                const std::string phy = fn.substr(7, fn.size() - 11);
                auto& shard = shards_[phy];
                if (!shard) shard = std::make_unique<PhylumShard>();
                std::unique_lock<std::shared_mutex> ul(shard->mu);
                shard->tree = VPTree::load(entry.path().string());
            }
        }
        dirty_.store(false, std::memory_order_release);
    }

private:
    struct PhylumShard {
        mutable std::shared_mutex mu;
        VPTree                    tree;
    };

    static std::string sanitize(const std::string& s) {
        std::string o; o.reserve(s.size());
        for (char c : s)
            o.push_back(std::isalnum(static_cast<unsigned char>(c)) ||
                        c == '_' || c == '-' ? c : '_');
        return o;
    }

    static CladeCentroid make_query_centroid(std::string_view seq,
                                              const SyncmerParams& sp,
                                              double densityThresh) {
        return make_query_centroid_for_routing(seq, sp, densityThresh);
    }

    void rebuild_locked() {
        for (auto& [phy, cents] : pendingByPhylum_) {
            auto& shard = shards_[phy];
            if (!shard) shard = std::make_unique<PhylumShard>();
            std::unique_lock<std::shared_mutex> ul(shard->mu);
            shard->tree = VPTree{};
            shard->tree.build(cents);
        }
        for (auto& [phy, cents] : pendingFbByPhylum_) {
            auto& shard = fbShards_[phy];
            if (!shard) shard = std::make_unique<PhylumShard>();
            std::unique_lock<std::shared_mutex> ul(shard->mu);
            shard->tree = VPTree{};
            shard->tree.build(cents);
        }
        dirty_.store(false, std::memory_order_release);
    }

    static uint64_t centroid_signature(const CladeCentroid& qc, size_t topK) {
        uint64_t sig = 1469598103934665603ULL ^ static_cast<uint64_t>(topK);
        const size_t take = std::min<size_t>(qc.centroidHashes.size(), 16);
        for (size_t i = 0; i < take; ++i) {
            sig ^= qc.centroidHashes[i] + 0x9e3779b97f4a7c15ULL + (sig << 6) + (sig >> 2);
        }
        sig ^= static_cast<uint64_t>(qc.centroidHashes.size());
        return sig;
    }

    mutable std::mutex registryMu_;
    std::atomic<bool>  dirty_{false};
    std::mutex routeCacheMu_;
    size_t routeCacheCapacity_ = 1024;
    std::unordered_map<uint64_t, std::vector<RouteResult>> routeCache_;
    std::deque<uint64_t> routeCacheOrder_;

    std::unordered_map<std::string, std::vector<CladeCentroid>> pendingByPhylum_;
    std::unordered_map<std::string, std::vector<CladeCentroid>> pendingFbByPhylum_;
    std::unordered_map<std::string, std::unique_ptr<PhylumShard>> shards_;
    std::unordered_map<std::string, std::unique_ptr<PhylumShard>> fbShards_;
};



// =========================================================================
// ExternalMemoryCentroidStore — streaming centroid operations for catalogs
// too large to keep resident in RAM. The store writes one compact centroid
// record at a time and supports chunked sequential scans, which keeps peak
// memory proportional to the chunk size rather than the full catalog size.
// =========================================================================
struct ExternalMemoryCentroidStore {
    struct DiskRecord {
        std::string           cladeName;
        std::string           phylum;
        std::string           cladeRank;
        std::vector<uint64_t> hashes;
    };

    struct SkipEntry {
        uint64_t key = 0;
        uint64_t offset = 0;
        uint64_t ordinal = 0;
    };

    struct SkipIndex {
        uint32_t stride = 16;
        std::vector<std::vector<SkipEntry>> levels;

        bool empty() const {
            return levels.empty() || levels.front().empty();
        }

        size_t size() const {
            return empty() ? 0u : levels.front().size();
        }
    };

    std::string path;
    mutable std::shared_ptr<SkipIndex> skipIndexCache;
    mutable std::mutex skipIndexMu;

    explicit ExternalMemoryCentroidStore(std::string p = {}) : path(std::move(p)) {}

    static void append_record(std::ostream& out, const DiskRecord& r) {
        auto write_string = [&](const std::string& s) {
            uint32_t n = static_cast<uint32_t>(s.size());
            out.write(reinterpret_cast<const char*>(&n), 4);
            out.write(s.data(), static_cast<std::streamsize>(n));
        };
        write_string(r.cladeName);
        write_string(r.phylum);
        write_string(r.cladeRank);
        uint32_t nh = static_cast<uint32_t>(r.hashes.size());
        out.write(reinterpret_cast<const char*>(&nh), 4);
        if (nh)
            out.write(reinterpret_cast<const char*>(r.hashes.data()),
                      static_cast<std::streamsize>(nh * sizeof(uint64_t)));
    }

    static bool read_record(std::istream& in, DiskRecord& r) {
        auto read_string = [&](std::string& s) -> bool {
            uint32_t n = 0;
            if (!in.read(reinterpret_cast<char*>(&n), 4)) return false;
            s.resize(n);
            return static_cast<bool>(in.read(s.data(), static_cast<std::streamsize>(n)));
        };
        r = {};
        if (!read_string(r.cladeName)) return false;
        if (!read_string(r.phylum))    return false;
        if (!read_string(r.cladeRank)) return false;
        uint32_t nh = 0;
        if (!in.read(reinterpret_cast<char*>(&nh), 4)) return false;
        r.hashes.resize(nh);
        if (nh && !in.read(reinterpret_cast<char*>(r.hashes.data()),
                           static_cast<std::streamsize>(nh * sizeof(uint64_t))))
            return false;
        return true;
    }

    static bool read_record_at(std::istream& in, uint64_t offset, DiskRecord& r) {
        in.clear();
        in.seekg(static_cast<std::streamoff>(offset), std::ios::beg);
        if (!in) return false;
        return read_record(in, r);
    }

    static uint64_t signature_key_from_hashes(const std::vector<uint64_t>& hashes) {
        uint64_t sig = 0x9e3779b97f4a7c15ULL;
        const size_t take = std::min<size_t>(hashes.size(), 8);
        for (size_t i = 0; i < take; ++i)
            sig ^= hashes[i] + 0x9e3779b97f4a7c15ULL + (sig << 6) + (sig >> 2);
        sig ^= static_cast<uint64_t>(hashes.size()) * 0xbf58476d1ce4e5b9ULL;
        return sig;
    }

    static uint64_t signature_key(const DiskRecord& r) {
        return signature_key_from_hashes(r.hashes);
    }

    std::string skip_path() const {
        return path + ".skip";
    }

    static void write_skip_index_file(const std::string& skipPath,
                                      const SkipIndex& index) {
        std::ofstream out(skipPath, std::ios::binary);
        if (!out) throw std::runtime_error("Cannot write centroid skip index: " + skipPath);
        const uint64_t magic = 0x31504b534c4f5455ULL; // "UTOLSKP1"
        const uint32_t version = 1;
        const uint32_t levelCount = static_cast<uint32_t>(index.levels.size());
        out.write(reinterpret_cast<const char*>(&magic), sizeof(magic));
        out.write(reinterpret_cast<const char*>(&version), sizeof(version));
        out.write(reinterpret_cast<const char*>(&index.stride), sizeof(index.stride));
        out.write(reinterpret_cast<const char*>(&levelCount), sizeof(levelCount));
        for (const auto& level : index.levels) {
            uint64_t n = static_cast<uint64_t>(level.size());
            out.write(reinterpret_cast<const char*>(&n), sizeof(n));
            for (const auto& e : level) {
                out.write(reinterpret_cast<const char*>(&e.key), sizeof(e.key));
                out.write(reinterpret_cast<const char*>(&e.offset), sizeof(e.offset));
                out.write(reinterpret_cast<const char*>(&e.ordinal), sizeof(e.ordinal));
            }
        }
    }

    static SkipIndex load_skip_index_file(const std::string& skipPath) {
        std::ifstream in(skipPath, std::ios::binary);
        if (!in) throw std::runtime_error("Cannot open centroid skip index: " + skipPath);
        uint64_t magic = 0;
        uint32_t version = 0;
        uint32_t stride = 0;
        uint32_t levelCount = 0;
        in.read(reinterpret_cast<char*>(&magic), sizeof(magic));
        in.read(reinterpret_cast<char*>(&version), sizeof(version));
        in.read(reinterpret_cast<char*>(&stride), sizeof(stride));
        in.read(reinterpret_cast<char*>(&levelCount), sizeof(levelCount));
        if (!in || magic != 0x31504b534c4f5455ULL || version != 1)
            throw std::runtime_error("Invalid centroid skip index: " + skipPath);
        SkipIndex index;
        index.stride = std::max<uint32_t>(2u, stride);
        index.levels.resize(levelCount);
        for (uint32_t li = 0; li < levelCount; ++li) {
            uint64_t n = 0;
            in.read(reinterpret_cast<char*>(&n), sizeof(n));
            auto& level = index.levels[static_cast<size_t>(li)];
            level.resize(static_cast<size_t>(n));
            for (auto& e : level) {
                in.read(reinterpret_cast<char*>(&e.key), sizeof(e.key));
                in.read(reinterpret_cast<char*>(&e.offset), sizeof(e.offset));
                in.read(reinterpret_cast<char*>(&e.ordinal), sizeof(e.ordinal));
            }
        }
        if (!in) throw std::runtime_error("Corrupt centroid skip index: " + skipPath);
        return index;
    }

    void build_skip_index(uint32_t stride = 16) const {
        std::ifstream in(path, std::ios::binary);
        if (!in) throw std::runtime_error("Cannot open external centroid store: " + path);
        uint64_t n = 0;
        in.read(reinterpret_cast<char*>(&n), sizeof(n));
        std::vector<SkipEntry> base;
        base.reserve(static_cast<size_t>(n));
        DiskRecord r;
        while (in) {
            const std::streamoff pos = in.tellg();
            if (pos < 0) break;
            if (!read_record(in, r)) break;
            base.push_back({signature_key(r), static_cast<uint64_t>(pos), 0});
        }
        std::sort(base.begin(), base.end(),
                  [](const SkipEntry& a, const SkipEntry& b) {
                      if (a.key != b.key) return a.key < b.key;
                      return a.offset < b.offset;
                  });
        for (size_t i = 0; i < base.size(); ++i)
            base[i].ordinal = static_cast<uint64_t>(i);

        SkipIndex index;
        index.stride = std::max<uint32_t>(2u, stride);
        index.levels.push_back(base);
        while (index.levels.back().size() > 1) {
            const auto& prev = index.levels.back();
            std::vector<SkipEntry> next;
            next.reserve((prev.size() + index.stride - 1) / index.stride);
            for (size_t i = 0; i < prev.size(); i += index.stride)
                next.push_back(prev[i]);
            if (!prev.empty() && next.back().ordinal != prev.back().ordinal)
                next.push_back(prev.back());
            if (next.size() >= prev.size()) {
                next.clear();
                next.push_back(prev.front());
            }
            index.levels.push_back(std::move(next));
        }

        write_skip_index_file(skip_path(), index);
        std::lock_guard<std::mutex> lk(skipIndexMu);
        skipIndexCache = std::make_shared<SkipIndex>(std::move(index));
    }

    std::shared_ptr<SkipIndex> ensure_skip_index() const {
        {
            std::lock_guard<std::mutex> lk(skipIndexMu);
            if (skipIndexCache) return skipIndexCache;
        }
        const std::string skipPath = skip_path();
        try {
            auto loaded = std::make_shared<SkipIndex>(load_skip_index_file(skipPath));
            std::lock_guard<std::mutex> lk(skipIndexMu);
            if (!skipIndexCache) skipIndexCache = loaded;
            return skipIndexCache;
        } catch (...) {
        }
        try {
            build_skip_index();
            std::lock_guard<std::mutex> lk(skipIndexMu);
            return skipIndexCache;
        } catch (...) {
        }
        return {};
    }

    void prepare_skip_index(uint32_t stride = 16) const {
        build_skip_index(stride);
    }

    void build(const std::vector<CladeCentroid>& centroids) const {
        {
            std::ofstream out(path, std::ios::binary);
            if (!out) throw std::runtime_error("Cannot write external centroid store: " + path);
            uint64_t n = static_cast<uint64_t>(centroids.size());
            out.write(reinterpret_cast<const char*>(&n), sizeof(n));
            for (const auto& c : centroids)
                append_record(out, {c.cladeName, c.phylum, c.cladeRank, c.centroidHashes});
        }
        build_skip_index();
    }

    template <class Fn>
    void for_each_record(Fn&& fn) const {
        std::ifstream in(path, std::ios::binary);
        if (!in) throw std::runtime_error("Cannot open external centroid store: " + path);
        uint64_t n = 0;
        in.read(reinterpret_cast<char*>(&n), sizeof(n));
        DiskRecord r;
        for (uint64_t i = 0; i < n && read_record(in, r); ++i) fn(r);
    }

    std::vector<VPTree::RouteResult>
    query_topk_streaming(const CladeCentroid& q, size_t k, size_t chunkRecords = 4096) const {
        auto push_best = [&](std::vector<std::pair<double, VPTree::RouteResult>>& best,
                             const DiskRecord& r) {
            CladeCentroid c;
            c.cladeName = r.cladeName;
            c.phylum = r.phylum;
            c.cladeRank = r.cladeRank;
            c.centroidHashes = r.hashes;
            c.build_prefilters();
            const double d = c.jaccard_distance(q);
            VPTree::RouteResult rr;
            rr.cladeName = r.cladeName;
            rr.phylum    = r.phylum;
            rr.jaccard   = 1.0 - d;
            best.push_back({d, rr});
            std::push_heap(best.begin(), best.end(),
                           [](const auto& a, const auto& b) { return a.first < b.first; });
            if (best.size() > k) {
                std::pop_heap(best.begin(), best.end(),
                              [](const auto& a, const auto& b) { return a.first < b.first; });
                best.pop_back();
            }
        };

        auto finalize_best = [](std::vector<std::pair<double, VPTree::RouteResult>> best) {
            std::sort(best.begin(), best.end(),
                      [](const auto& a, const auto& b) { return a.first < b.first; });
            std::vector<VPTree::RouteResult> out;
            out.reserve(best.size());
            for (auto& kv : best) out.push_back(std::move(kv.second));
            return out;
        };

        auto indexed = ensure_skip_index();
        if (indexed && !indexed->empty()) {
            const auto& base = indexed->levels.front();
            const uint64_t qKey = signature_key_from_hashes(q.centroidHashes);
            auto locate_floor_ordinal = [&](uint64_t key) {
                size_t ord = 0;
                for (size_t li = indexed->levels.size(); li-- > 0;) {
                    const auto& level = indexed->levels[li];
                    auto it = std::upper_bound(
                        level.begin(), level.end(), key,
                        [](uint64_t lhs, const SkipEntry& rhs) { return lhs < rhs.key; });
                    if (it != level.begin()) {
                        --it;
                        ord = std::max<size_t>(ord, static_cast<size_t>(it->ordinal));
                    }
                }
                while (ord + 1 < base.size() && base[ord + 1].key <= key) ++ord;
                return ord;
            };

            const size_t center = locate_floor_ordinal(qKey);
            size_t window = std::min(base.size(), std::max<size_t>(128, chunkRecords));
            const size_t maxWindow =
                std::min(base.size(), std::max<size_t>(4096, chunkRecords * 8));
            for (;;) {
                const size_t half = window / 2;
                const size_t lo = center > half ? center - half : 0;
                const size_t hi = std::min(base.size(), lo + window);
                std::vector<uint64_t> offsets;
                offsets.reserve(hi - lo);
                for (size_t i = lo; i < hi; ++i) offsets.push_back(base[i].offset);
                std::sort(offsets.begin(), offsets.end());

                std::ifstream in(path, std::ios::binary);
                std::vector<std::pair<double, VPTree::RouteResult>> best;
                best.reserve(k + 1);
                DiskRecord r;
                for (uint64_t offset : offsets) {
                    if (!read_record_at(in, offset, r)) continue;
                    push_best(best, r);
                }
                auto out = finalize_best(std::move(best));
                if (!out.empty() && out.front().jaccard >= 0.999999999999) return out;
                if (window >= maxWindow || window >= base.size()) break;
                window = std::min(base.size(), window * 2);
            }
        }

        std::vector<std::pair<double, VPTree::RouteResult>> best;
        best.reserve(k + 1);

        std::ifstream in(path, std::ios::binary);
        if (!in) return {};
        uint64_t n = 0;
        in.read(reinterpret_cast<char*>(&n), sizeof(n));
        DiskRecord r;
        size_t buffered = 0;
        for (uint64_t i = 0; i < n && read_record(in, r); ++i) {
            push_best(best, r);
            if (++buffered >= std::max<size_t>(1, chunkRecords)) buffered = 0;
        }
        return finalize_best(std::move(best));
    }
};

using ScalableRoutingResult = PhylumShardedRouter::RouteResult;
using RouteResult           = PhylumShardedRouter::RouteResult;
using ScalableRoutingIndex  = struct _ScalableRoutingIndexTag {};

inline std::vector<ScalableRoutingResult>
route_query(PhylumShardedRouter& router,
            std::string_view seq,
            const SyncmerParams& sp,
            const SyncmerParams& fbSp,
            double density,
            size_t topK) {
    return router.route(seq, sp, fbSp, density, topK);
}

} // namespace tol
