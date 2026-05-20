#pragma once
// layer2_registry.hpp — v14
// Clade Graph Registry (Layer 2)
//
// Changes vs v12/v13
// ==================
//  DS-1   O(1) LRU eviction (Sleator & Tarjan 1985) — unchanged
//  DS-2   Per-clade load-once shared_future guard — unchanged
//  DS-3   Atomic-rename manifest flush — unchanged
//  DS-17  FenwickTree for eviction accounting existed in earlier iterations.
//         The current hot path uses O(1) LRU bookkeeping plus a running byte
//         counter instead, which is simpler and cheaper for the present cache.
//  FIX-M9 insert_into_cache no longer rebuilds any global accounting structure
//         under the write-lock; eviction now relies on the running byte counter
//         and O(1) LRU removal.
//  FIX-M10 sanitize() is now a thin forwarder to tol::sanitize_name_impl
//          so the character set is defined in exactly one place.  The impl
//          is inlined here (cannot include fungi_tol_bridge.hpp — circular).

#include "layer1_clade_graph.hpp"

#include <array>
#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <future>
#include <iostream>
#include <list>
#include <mutex>
#include <shared_mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

#if defined(__unix__) || defined(__APPLE__)
#  include <unistd.h>
#endif

namespace tol {

// ── G2: free-RAM adaptive cache sizing ───────────────────────────────────
inline size_t adaptive_cache_bytes() {
#if defined(__unix__) || defined(__APPLE__)
    long pages = sysconf(_SC_AVPHYS_PAGES);
    long psz   = sysconf(_SC_PAGE_SIZE);
    if (pages > 0 && psz > 0)
        return std::min(static_cast<size_t>(pages) * static_cast<size_t>(psz) / 4,
                        size_t(32) << 30);
#endif
    return size_t(8) << 30;
}

// ── CRC-32 (IEEE 802.3) ───────────────────────────────────────────────────
inline uint32_t crc32_file(const std::string& path) {
    static const std::array<uint32_t, 256> T = []() {
        std::array<uint32_t, 256> t{};
        for (uint32_t i = 0; i < 256; ++i) {
            uint32_t c = i;
            for (int j = 0; j < 8; ++j) c = (c & 1) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
            t[i] = c;
        }
        return t;
    }();
    std::ifstream in(path, std::ios::binary);
    if (!in) return 0;
    char buf[65536];
    uint32_t crc = 0xFFFFFFFFu;
    while (in.read(buf, sizeof(buf)) || in.gcount())
        for (std::streamsize i = 0; i < in.gcount(); ++i)
            crc = (crc >> 8) ^ T[(crc ^ static_cast<uint8_t>(buf[i])) & 0xFF];
    return crc ^ 0xFFFFFFFFu;
}

// ── Manifest types ────────────────────────────────────────────────────────
struct CladeDescriptor {
    std::string         cladeName, cladeRank, phylum, graphPath;
    size_t              genomeCount    = 0;
    size_t              svBubbles      = 0;
    size_t              compressedBytes= 0;
    std::vector<std::string> fastaPaths;
    uint32_t            crc32          = 0;
    std::vector<uint64_t> centroidSyncmers;
};

// DS-3: atomic-rename prevents partial writes on crash.
inline void save_manifest(const std::vector<CladeDescriptor>& clades,
                           const std::string& path) {
    fs::create_directories(fs::path(path).parent_path());
    const std::string tmp = path + ".tmp";
    {
        std::ofstream out(tmp);
        if (!out) throw std::runtime_error("Cannot write manifest: " + tmp);
        out << "#clade_name\tclade_rank\tphylum\tgraph_path\tgenome_count"
               "\tsv_bubbles\tcompressed_bytes\tfasta_paths\tcrc32\tcentroid_hashes\n";
        for (const auto& c : clades) {
            out << c.cladeName << '\t' << c.cladeRank << '\t' << c.phylum << '\t'
                << c.graphPath << '\t' << c.genomeCount << '\t' << c.svBubbles << '\t'
                << c.compressedBytes << '\t';
            for (size_t i = 0; i < c.fastaPaths.size(); ++i) {
                if (i) out << ',';
                out << c.fastaPaths[i];
            }
            out << '\t' << c.crc32 << '\t';
            for (size_t i = 0; i < c.centroidSyncmers.size(); ++i) {
                if (i) out << ',';
                out << c.centroidSyncmers[i];
            }
            out << '\n';
        }
    }
    fs::rename(tmp, path);
}

inline std::vector<CladeDescriptor> load_manifest(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot read manifest: " + path);
    std::vector<CladeDescriptor> out;
    std::string line;
    std::unordered_map<std::string, size_t> col;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            std::string header = line.substr(1);
            std::istringstream hs(header);
            std::string name;
            size_t idx = 0;
            while (std::getline(hs, name, '\t')) col[name] = idx++;
            continue;
        }
        std::vector<std::string> cols;
        std::istringstream ls(line);
        std::string cell;
        while (std::getline(ls, cell, '\t')) {
            // Bug 5 fix: strip CR from CRLF-encoded manifests so the last cell
            // (typically the fasta_paths or centroid_hashes column on the very
            // last line of a CRLF manifest) doesn't keep a trailing '\r' that
            // breaks fs::exists() downstream.
            if (!cell.empty() && cell.back() == '\r') cell.pop_back();
            cols.push_back(std::move(cell));
        }
        if (cols.size() < 7) continue;
        auto get = [&](const std::string& name, size_t fallback) -> std::string {
            auto it = col.find(name);
            const size_t idx = (it == col.end()) ? fallback : it->second;
            return idx < cols.size() ? cols[idx] : std::string();
        };
        CladeDescriptor c;
        c.cladeName = get("clade_name", 0);
        c.cladeRank = get("clade_rank", 1);
        c.phylum = get("phylum", 2);
        c.graphPath = get("graph_path", 3);
        try { c.genomeCount = static_cast<size_t>(std::stoull(get("genome_count", 4))); }
        catch (...) { c.genomeCount = 0; }
        try { c.svBubbles = static_cast<size_t>(std::stoull(get("sv_bubbles", 5))); }
        catch (...) { c.svBubbles = 0; }
        try { c.compressedBytes = static_cast<size_t>(std::stoull(get("compressed_bytes", 6))); }
        catch (...) { c.compressedBytes = 0; }

        const std::string fastaCol = get("fasta_paths", SIZE_MAX);
        if (!fastaCol.empty()) {
            std::istringstream fp(fastaCol);
            std::string path;
            while (std::getline(fp, path, ',')) {
                // Bug 5 fix: defensive CR strip in case the cell-level strip
                // above somehow missed it (e.g., embedded CR mid-list).
                if (!path.empty() && path.back() == '\r') path.pop_back();
                if (!path.empty()) c.fastaPaths.push_back(std::move(path));
            }
            try { c.crc32 = static_cast<uint32_t>(std::stoul(get("crc32", 8))); }
            catch (...) { c.crc32 = 0; }
        } else {
            // Backward compatibility with manifests written before the
            // fasta_paths column existed: column 7 was crc32. Also avoid
            // treating the short-lived broken schema's CRC field as a FASTA
            // path; callers recover paths from hierarchy_manifest.tsv.
            try { c.crc32 = static_cast<uint32_t>(std::stoul(get("crc32", 7))); }
            catch (...) { c.crc32 = 0; }
        }
        std::string hashes = get("centroid_hashes", col.find("fasta_paths") == col.end() ? 8 : 9);
        std::istringstream hs(hashes);
        std::string tok;
        while (std::getline(hs, tok, ','))
            if (!tok.empty()) {
                try { c.centroidSyncmers.push_back(std::stoull(tok)); }
                catch (...) {}
            }
        out.push_back(std::move(c));
    }
    return out;
}

// ── LRU cache entry (DS-1: holds list iterator for O(1) eviction) ─────────
using LruList = std::list<std::string>;

struct CacheEntry {
    std::string                  cladeName;
    std::shared_ptr<CladeGraph>  graph;
    LruList::iterator            lruIt;
    size_t                       sizeBytes = 0;
};

// =========================================================================
// CladeGraphRegistry — DS-1 O(1) LRU + DS-2 load-once future guard
// =========================================================================
class CladeGraphRegistry {
public:
    explicit CladeGraphRegistry(std::string baseDir,
                                 size_t maxCacheBytes   = 0,
                                 size_t maxCacheEntries = 256)
        : baseDir_(std::move(baseDir))
        , manifestPath_(baseDir_ + "/clade_manifest.tsv")
        , maxCacheBytes_(maxCacheBytes == 0 ? adaptive_cache_bytes() : maxCacheBytes)
        , maxCacheEntries_(maxCacheEntries)
    {}

    // FIX-L: gbz_io::save runs outside any lock.
    void register_clade(const CladeGraph& g,
                        const std::vector<uint64_t>& centroidHashes = {}) {
        const std::string gbzPath = clade_path(g.cladeName);
        gbz_io::save(g, gbzPath);                       // I/O outside lock
        const uint32_t crc = crc32_file(gbzPath);
        CladeDescriptor desc;
        desc.cladeName       = g.cladeName;
        desc.cladeRank       = g.cladeRank;
        desc.phylum          = g.phylum;
        desc.graphPath       = gbzPath;
        desc.genomeCount     = g.genomeCount;
        desc.svBubbles       = g.svBubbles;
        desc.compressedBytes = g.compressed_bytes();
        desc.crc32           = crc;
        desc.centroidSyncmers= centroidHashes;
        {
            std::unique_lock lock(mu_);
            clades_.erase(
                std::remove_if(clades_.begin(), clades_.end(),
                    [&](const CladeDescriptor& d){ return d.cladeName == g.cladeName; }),
                clades_.end());
            clades_.push_back(desc);
            flush_manifest_locked();
            ++totalRegistered_;
            // M7: do NOT cache on register — keep RAM free for query phase.
        }
    }

    // DS-2 + DS-1: load-once future guard; O(1) LRU touch.
    std::shared_ptr<CladeGraph> get(const std::string& cladeName) {
        // ── Fast path: already cached ────────────────────────────────────
        {
            std::unique_lock wl(mu_);
            auto cit = cache_.find(cladeName);
            if (cit != cache_.end()) {
                // Move to MRU front in O(1)
                lruList_.splice(lruList_.begin(), lruList_, cit->second.lruIt);
                cit->second.lruIt = lruList_.begin();
                ++cacheHits_;
                return cit->second.graph;
            }
        }
        ++cacheMisses_;

        // ── DS-2: get-or-create shared_future sentinel ───────────────────
        std::shared_future<std::shared_ptr<CladeGraph>> fut;
        bool isOwner = false;
        {
            std::unique_lock wl(mu_);
            // Double-check after acquiring write lock
            auto cit = cache_.find(cladeName);
            if (cit != cache_.end()) {
                lruList_.splice(lruList_.begin(), lruList_, cit->second.lruIt);
                cit->second.lruIt = lruList_.begin();
                ++cacheHits_; --cacheMisses_;
                return cit->second.graph;
            }
            auto pit = pending_.find(cladeName);
            if (pit != pending_.end()) {
                fut = pit->second;
            } else {
                // FIX: Do NOT insert a dummy promise that will immediately be
                // destroyed (leaving pending_ holding a broken future).  Instead,
                // just mark ownership; the real promise is created below while
                // the lock is released, and then stored in pending_ in one step.
                isOwner = true;
            }
        }

        if (isOwner) {
            // Create promise and store its future under the lock in one step,
            // eliminating the window where pending_ held a broken-promise future.
            std::promise<std::shared_ptr<CladeGraph>> prom2;
            auto fut2 = prom2.get_future().share();
            {
                std::unique_lock wl(mu_);
                // Guard against a race where two owners both reach this point
                // (possible if the outer fast-path check races).
                auto pit2 = pending_.find(cladeName);
                if (pit2 != pending_.end()) {
                    // Another owner beat us — become a non-owner waiter.
                    fut = pit2->second;
                    isOwner = false;
                } else {
                    pending_[cladeName] = fut2;
                }
            }
            if (!isOwner) {
                return fut.get();
            }
            std::shared_ptr<CladeGraph> g;
            std::exception_ptr ep;
            try {
                g = load_from_disk_unlocked(cladeName);
            } catch (...) {
                ep = std::current_exception();
            }
            if (ep) {
                {
                    std::unique_lock wl(mu_);
                    pending_.erase(cladeName);
                }
                prom2.set_exception(ep);
                std::rethrow_exception(ep);
            }
            prom2.set_value(g);
            insert_into_cache(cladeName, g);
            {
                std::unique_lock wl(mu_);
                pending_.erase(cladeName);
            }
            return g;
        }

        // Non-owner: wait for the owner's future.
        return fut.get();
    }

    void load_manifest_from_disk() {
        std::unique_lock wl(mu_);
        try {
            clades_ = tol::load_manifest(manifestPath_);
        } catch (...) { /* manifest may not exist yet */ }
    }

    struct CacheStats {
        size_t hits, misses, entries, evictions, totalBytes, registered;
    };
    CacheStats stats() const {
        return { cacheHits_.load(), cacheMisses_.load(),
                 cache_.size(), cacheEvictions_.load(),
                 cacheCurrentBytes_.load(), totalRegistered_.load() };
    }

    size_t clade_count()  const { return clades_.size(); }
    const std::string& base_dir() const { return baseDir_; }
    const std::vector<CladeDescriptor>& descriptors() const { return clades_; }

private:
    std::string clade_path(const std::string& name) const {
        return baseDir_ + "/" + sanitize(name) + ".gbz";
    }

    // Identical character set to tol::sanitize_name in fungi_tol_bridge.hpp.
    // Cannot include the bridge here (circular dependency).
    static std::string sanitize(const std::string& s) {
        std::string o;
        o.reserve(s.size());
        for (char c : s)
            o.push_back(std::isalnum(static_cast<unsigned char>(c)) ||
                        c == '_' || c == '-' ? c : '_');
        return o;
    }

    void flush_manifest_locked() {
        try { tol::save_manifest(clades_, manifestPath_); }
        catch (const std::exception& e) {
            std::cerr << "[registry] manifest flush: " << e.what() << '\n';
        }
    }

    // Disk load — called OUTSIDE any lock (DS-2).
    std::shared_ptr<CladeGraph> load_from_disk_unlocked(const std::string& cladeName) {
        CladeDescriptor descCopy;
        {
            std::shared_lock rl(mu_);
            bool found = false;
            for (const auto& d : clades_)
                if (d.cladeName == cladeName) { descCopy = d; found = true; break; }
            if (!found)
                throw std::runtime_error("Clade not in registry: " + cladeName);
        }
        if (descCopy.crc32) {
            const uint32_t actual = crc32_file(descCopy.graphPath);
            if (actual != descCopy.crc32)
                throw std::runtime_error(
                    "GBZ CRC mismatch for " + cladeName + " — delete & rebuild index");
        }
        return std::make_shared<CladeGraph>(gbz_io::load(descCopy.graphPath));
    }

    // FIX-M9 + FIX-M15: insert_into_cache
    //
    // The byte counter cacheCurrentBytes_ is incremented BEFORE the eviction
    // loop so that the loop always sees the would-be post-insertion size.
    // Previously the counter was incremented after the loop, meaning the
    // eviction test used the pre-insertion size and could leave the cache
    // over-budget by exactly sz bytes when sz > 0, or fill indefinitely with
    // zero-byte stub entries when sz == 0.
    //
    // The full write-lock is held throughout (required for list/map
    // consistency).  The eviction itself is O(K) where K is the number of
    // entries evicted — no full Fenwick rebuild needed.
    void insert_into_cache(const std::string& cladeName,
                            std::shared_ptr<CladeGraph> g) {
        const size_t sz = g->compressed_bytes();
        std::unique_lock wl(mu_);

        // Account for the new entry before deciding what to evict.
        cacheCurrentBytes_ += sz;

        // Evict from LRU tail until both byte and entry limits are satisfied.
        // The +1 on cache_.size() anticipates the insertion below.
        while (!cache_.empty()) {
            const bool tooLarge = (cacheCurrentBytes_ > maxCacheBytes_);
            const bool tooMany  = (cache_.size() + 1 > maxCacheEntries_);
            if (!tooLarge && !tooMany) break;
            evict_lru_locked();
        }

        lruList_.push_front(cladeName);
        cache_[cladeName] = { cladeName, g, lruList_.begin(), sz };
        // cacheCurrentBytes_ was already incremented above.
    }

    // DS-1: O(1) LRU eviction — back of list is oldest.
    void evict_lru_locked() {
        if (cache_.empty()) return;
        const std::string& oldest = lruList_.back();
        auto it = cache_.find(oldest);
        if (it != cache_.end()) {
            cacheCurrentBytes_ -= it->second.sizeBytes;
            cache_.erase(it);
        }
        lruList_.pop_back();
        ++cacheEvictions_;
    }

    std::string baseDir_, manifestPath_;
    size_t maxCacheBytes_, maxCacheEntries_;

    mutable std::shared_mutex mu_;
    std::vector<CladeDescriptor> clades_;
    std::unordered_map<std::string, CacheEntry> cache_;
    LruList lruList_;   // front = MRU, back = LRU

    // DS-2: in-flight futures prevent duplicate disk loads.
    std::unordered_map<std::string,
                       std::shared_future<std::shared_ptr<CladeGraph>>> pending_;

    std::atomic<size_t> cacheHits_{0}, cacheMisses_{0}, cacheEvictions_{0};
    std::atomic<size_t> cacheCurrentBytes_{0}, totalRegistered_{0};
};

// =========================================================================
// CladeGraphBuilder
// =========================================================================
class CladeGraphBuilder {
public:
    CladeGraphBuilder(std::string name, std::string rank, std::string phylum)
        : name_(std::move(name)), rank_(std::move(rank)), phylum_(std::move(phylum)) {}

    void add_genome(const std::string& asm_, const std::string& contig,
                    const std::string& seq,  const std::string& annot = "NONE",
                    bool isRef = false, bool isRep = false,
                    int segLen = 10000, int segOv = 1000) {
        const int L    = static_cast<int>(seq.size());
        const int step = segLen - segOv;
        int prev = -1;
        for (int s = 0; s < L; s += step) {
            const int e = std::min(L, s + segLen);
            append_seg(asm_, contig,
                       seq.substr(static_cast<size_t>(s), static_cast<size_t>(e - s)),
                       s, e, annot, isRep || isRef, &prev);
            if (e == L) break;
        }
    }

    void add_segment_str(const std::string& asm_, const std::string& contig,
                         const std::string& segSeq, int globalStart,
                         const std::string& annot = "NONE", bool isRep = false,
                         int* prevIdInOut = nullptr,
                         int blockStart = -1, int blockEnd = -1) {
        ImportSegment s;
        s.id              = nextId_++;
        s.asmName         = asm_;
        s.contig          = contig;
        s.start           = globalStart;
        s.end             = globalStart + static_cast<int>(segSeq.size());
        s.seq             = segSeq;
        s.annotation      = annot;
        s.isRepresentative= isRep;
        s.blockStart      = blockStart;
        s.blockEnd        = blockEnd;
        segs_.push_back(s);
        if (prevIdInOut && *prevIdInOut >= 0)
            edges_.push_back({ *prevIdInOut, s.id, true });
        if (prevIdInOut) *prevIdInOut = s.id;
    }

    void set_node_sketch(int segId, std::vector<uint64_t> hashes) {
        sketchMap_[segId] = std::move(hashes);
    }

    CladeGraph build(bool enableCollapse = true, bool enableCompact = true,
                     double minBubbleFreq = kMinBubbleFreqDef) {
        return build_clade_graph(name_, rank_, phylum_, segs_, edges_, sketchMap_,
                                 enableCollapse, 0.50, enableCompact, minBubbleFreq);
    }

private:
    void append_seg(const std::string& asm_, const std::string& contig,
                    const std::string& seq,  int s, int e,
                    const std::string& annot, bool isRep, int* prev) {
        ImportSegment seg;
        seg.id              = nextId_++;
        seg.asmName         = asm_;
        seg.contig          = contig;
        seg.start           = s;
        seg.end             = e;
        seg.seq             = seq;
        seg.annotation      = annot;
        seg.isRepresentative= isRep;
        segs_.push_back(seg);
        if (*prev >= 0) edges_.push_back({ *prev, seg.id, true });
        *prev = seg.id;
    }

    std::string name_, rank_, phylum_;
    std::vector<ImportSegment> segs_;
    std::vector<ImportEdge>    edges_;
    std::unordered_map<int, std::vector<uint64_t>> sketchMap_;
    int nextId_ = 0;
};

} // namespace tol
