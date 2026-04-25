#ifndef QUERY_INPUT_HANDLER_HPP
#define QUERY_INPUT_HANDLER_HPP

// query_input_handler.hpp — v1
// ============================================================
// Generalises the query input layer so the same SV-calling pipeline
// accepts any of three input types from the user:
//
//   ASSEMBLY    — pre-assembled contigs in FASTA (existing behaviour)
//   LONG_READS  — ONT or PacBio reads in FASTA or FASTQ
//                 (typical coverage 5-50x, read length 1-100 kb)
//   SHORT_READS — Illumina reads in FASTA or FASTQ
//                 (typical coverage 20-100x, read length 50-300 bp)
//
// The mode is either:
//   * set by the user via --query-mode {assembly|long-reads|short-reads}
//   * auto-detected from file extension + read-length distribution
//     when --query-mode is omitted or "auto"
//
// For LONG_READS and SHORT_READS the handler converts raw reads into
// pseudo-contigs that are passed unchanged into the existing calling engine.
// No external library or build dependency is required.
//
// LONG_READS  -> greedy k-mer-overlap consensus per read cluster.
//               Reads are grouped by shared anchor k-mers (k=12), then a
//               per-position majority-vote sequence is computed per group.
//               Handles >=3x coverage and ~15% error rate.
//
// SHORT_READS -> de-Bruijn k-mer path extension (k=21).
//               Solid k-mers (frequency >= auto-detected threshold) are
//               extended greedily into unitigs (SPAdes-style simple paths).
//
// Calling-parameter auto-tuning (applied unless --no-mode-param-override):
//   LONG_READS  -- k=15, chainGapBand=15000, minAnchors=1, minBlockScore=3.0,
//                  secondary seeds on with secondaryK=11
//   SHORT_READS -- minAnchors=3, secondary seeds on (k=21 unchanged)
//   ASSEMBLY    -- no changes to any defaults
// ============================================================

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <ext/stdio_filebuf.h>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifndef __declspec
#  ifndef _WIN32
#    define __declspec(x) __attribute__((visibility("default")))
#  endif
#endif

namespace fs = std::filesystem;

namespace query_input {

// ----------------------------------------------------------------
// QueryMode
// ----------------------------------------------------------------
enum class QueryMode { ASSEMBLY, LONG_READS, SHORT_READS };

inline const char* mode_name(QueryMode m) {
    switch (m) {
        case QueryMode::ASSEMBLY:    return "assembly";
        case QueryMode::LONG_READS:  return "long-reads";
        case QueryMode::SHORT_READS: return "short-reads";
    }
    return "assembly";
}

inline QueryMode parse_mode(const std::string& s) {
    std::string lo = s;
    for (char& c : lo) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (lo == "assembly" || lo == "asm")                          return QueryMode::ASSEMBLY;
    if (lo == "long-reads" || lo == "long_reads"
     || lo == "long" || lo == "lr" || lo == "longreads")          return QueryMode::LONG_READS;
    if (lo == "short-reads" || lo == "short_reads"
     || lo == "short" || lo == "sr" || lo == "illumina"
     || lo == "shortreads")                                        return QueryMode::SHORT_READS;
    throw std::runtime_error(
        "Unknown --query-mode '" + s + "'. Valid: assembly, long-reads, short-reads");
}

// ----------------------------------------------------------------
// InputConfig  (filled from CLI Options by the caller)
// ----------------------------------------------------------------
struct InputConfig {
    // Long-reads preprocessing
    int      lrAnchorK         = 12;
    size_t   lrMinCluster      = 2;
    size_t   lrMinReadLen      = 200;
    size_t   lrMaxReadLen      = 300000;

    // Short-reads preprocessing
    int      srK               = 21;
    uint32_t srMinKmerFreq     = 0;      // 0 = auto from median frequency
    size_t   srMinUnitigLen    = 200;
    size_t   srMinReadLen      = 50;
    size_t   srMaxReadLen      = 600;

    // Coverage
    size_t   genomeSizeHint    = 0;      // 0 = skip coverage estimate

    // Safety cap on reads loaded per file
    size_t   maxReadsPerFile   = 10000000;
};

// ----------------------------------------------------------------
// CoverageReport
// ----------------------------------------------------------------
enum class CoverageTier { UNKNOWN, LOW, NORMAL, HIGH };

inline const char* coverage_tier_name(CoverageTier t) {
    switch (t) {
        case CoverageTier::UNKNOWN: return "unknown";
        case CoverageTier::LOW:     return "low";
        case CoverageTier::NORMAL:  return "normal";
        case CoverageTier::HIGH:    return "high";
    }
    return "unknown";
}

inline CoverageTier classify_coverage_tier(QueryMode mode, double cov) {
    if (cov <= 0.0) return CoverageTier::UNKNOWN;
    switch (mode) {
        case QueryMode::ASSEMBLY:
            if (cov < 5.0) return CoverageTier::LOW;
            if (cov >= 60.0) return CoverageTier::HIGH;
            return CoverageTier::NORMAL;
        case QueryMode::LONG_READS:
            if (cov < 8.0) return CoverageTier::LOW;
            if (cov >= 35.0) return CoverageTier::HIGH;
            return CoverageTier::NORMAL;
        case QueryMode::SHORT_READS:
            if (cov < 12.0) return CoverageTier::LOW;
            if (cov >= 60.0) return CoverageTier::HIGH;
            return CoverageTier::NORMAL;
    }
    return CoverageTier::UNKNOWN;
}


struct EvidenceFusionRecommendation {
    bool enableProbabilisticFusion = false;
    bool trustAssemblyAsPrior       = true;
    double priorAlt                = 0.50;
    double expectedAssemblyWeight  = 1.0;
    double expectedReadWeight      = 1.0;
    std::string rationale;
};

inline EvidenceFusionRecommendation
recommend_evidence_fusion(QueryMode mode, CoverageTier tier, size_t pseudoContigs = 0) {
    EvidenceFusionRecommendation rec;
    rec.enableProbabilisticFusion = true;
    rec.trustAssemblyAsPrior = (mode == QueryMode::ASSEMBLY || pseudoContigs > 0);
    switch (mode) {
        case QueryMode::ASSEMBLY:
            rec.priorAlt = 0.55;
            rec.expectedAssemblyWeight = 1.40;
            rec.expectedReadWeight = 0.60;
            rec.rationale = "assembly_mode_graph_calling";
            break;
        case QueryMode::LONG_READS:
            rec.priorAlt = (tier == CoverageTier::LOW ? 0.45 : 0.50);
            rec.expectedAssemblyWeight = (pseudoContigs > 0 ? 1.20 : 0.85);
            rec.expectedReadWeight = (tier == CoverageTier::HIGH ? 1.35 : 1.00);
            rec.rationale = (tier == CoverageTier::HIGH)
                ? "long_reads_high_coverage_multilayer_fusion"
                : "long_reads_balanced_multilayer_fusion";
            break;
        case QueryMode::SHORT_READS:
            rec.priorAlt = (tier == CoverageTier::LOW ? 0.40 : 0.50);
            rec.expectedAssemblyWeight = (pseudoContigs > 0 ? 1.30 : 0.70);
            rec.expectedReadWeight = (tier == CoverageTier::HIGH ? 1.45 : 1.05);
            rec.rationale = (tier == CoverageTier::HIGH)
                ? "short_reads_depth_stabilised_fusion"
                : "short_reads_balanced_fusion";
            break;
    }
    return rec;
}

struct CoverageReport {
    QueryMode mode            = QueryMode::ASSEMBLY;
    size_t    totalReads      = 0;
    size_t    totalBases      = 0;
    size_t    minLen          = 0;
    size_t    maxLen          = 0;
    double    meanLen         = 0.0;
    double    n50             = 0.0;
    double    estimatedCov    = 0.0;
    size_t    pseudoContigs   = 0;
    size_t    droppedReads    = 0;
    CoverageTier coverageTier = CoverageTier::UNKNOWN;
    std::string strategyName;

    void print(std::ostream& err, const std::string& name) const {
        err << "[query-input] " << name
            << " mode=" << mode_name(mode)
            << " reads=" << totalReads
            << " bases=" << totalBases
            << " mean_len=" << static_cast<long>(meanLen)
            << " N50=" << static_cast<long>(n50);
        if (estimatedCov > 0.0)
            err << " cov~" << std::fixed << std::setprecision(1) << estimatedCov << "x";
        err << " cov_tier=" << coverage_tier_name(coverageTier);
        if (!strategyName.empty()) err << " strategy=" << strategyName;
        err << " pseudo_contigs=" << pseudoContigs
            << " dropped=" << droppedReads << '\n';
    }
};

// ----------------------------------------------------------------
// PreparedQuery — output of prepare_query()
// ----------------------------------------------------------------
struct PreparedQuery {
    std::string sampleName;
    std::unordered_map<std::string, std::string> contigs;  // name -> sequence
    CoverageReport report;
};

// ================================================================
// Internal implementation details
// ================================================================
namespace detail {

inline bool is_fastq(const std::string& path) {
    std::string ext = fs::path(path).extension().string();
    for (char& c : ext) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (ext == ".fastq" || ext == ".fq") return true;
    std::string stem_ext = fs::path(fs::path(path).stem()).extension().string();
    for (char& c : stem_ext) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return stem_ext == ".fastq" || stem_ext == ".fq";
}

inline bool is_gzipped(const std::string& path) {
    return path.size() >= 3 && path.substr(path.size() - 3) == ".gz";
}

inline void sanitise_seq(std::string& seq) {
    for (char& c : seq) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') c = 'N';
    }
}

struct RawRead { std::string name, seq; };

// Reads FASTA or FASTQ (detected by first character), up to maxRecords.
// Supports plain files and .gz files (piped through gzip -dc).
inline std::vector<RawRead>
load_reads(const std::string& path, size_t maxRecords, bool quiet) {
    // Open the input: plain file or gz pipe backed by __gnu_cxx::stdio_filebuf.
    FILE* gz_pipe = nullptr;
    std::unique_ptr<__gnu_cxx::stdio_filebuf<char>> gz_buf;
    std::unique_ptr<std::istream> gz_is;
    std::ifstream plain_file;
    std::istream* in_ptr = nullptr;

    if (is_gzipped(path)) {
        std::string cmd = "gzip -dc '" + path + "'";
        gz_pipe = popen(cmd.c_str(), "r");
        if (!gz_pipe) {
            if (!quiet)
                std::cerr << "[query-input] ERROR: cannot decompress " << path << '\n';
            return {};
        }
        gz_buf = std::make_unique<__gnu_cxx::stdio_filebuf<char>>(gz_pipe, std::ios::in);
        gz_is  = std::make_unique<std::istream>(gz_buf.get());
        in_ptr = gz_is.get();
    } else {
        plain_file.open(path);
        if (!plain_file) throw std::runtime_error("Cannot open query file: " + path);
        in_ptr = &plain_file;
    }
    std::istream& in = *in_ptr;

    std::vector<RawRead> reads;
    reads.reserve(std::min(maxRecords, static_cast<size_t>(100000)));
    std::string line, name, seq;

    // Detect format from first non-blank character
    char first = '\0';
    while (in.peek() != EOF) {
        first = static_cast<char>(in.peek());
        if (first != '\n' && first != '\r') break;
        in.get();
    }
    const bool isFastq = (first == '@');

    if (!isFastq) {
        // FASTA
        while (reads.size() < maxRecords && std::getline(in, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty()) continue;
            if (line[0] == '>') {
                if (!name.empty() && !seq.empty()) {
                    sanitise_seq(seq);
                    reads.push_back({name, std::move(seq)});
                    seq.clear();
                }
                name = line.substr(1);
                const auto sp = name.find_first_of(" \t");
                if (sp != std::string::npos) name.resize(sp);
            } else {
                seq += line;
            }
        }
        if (!name.empty() && !seq.empty() && reads.size() < maxRecords) {
            sanitise_seq(seq); reads.push_back({name, std::move(seq)});
        }
    } else {
        // FASTQ (4-line records)
        int slot = 0;
        while (reads.size() < maxRecords && std::getline(in, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            const int s = slot % 4;
            if (s == 0) {
                if (!name.empty() && !seq.empty()) {
                    sanitise_seq(seq);
                    reads.push_back({name, std::move(seq)});
                    seq.clear();
                }
                name = (line.size() > 1) ? line.substr(1) : std::string{};
                const auto sp = name.find_first_of(" \t");
                if (sp != std::string::npos) name.resize(sp);
            } else if (s == 1) {
                seq = line;
            }
            ++slot;
        }
        if (!name.empty() && !seq.empty() && reads.size() < maxRecords) {
            sanitise_seq(seq); reads.push_back({name, std::move(seq)});
        }
    }

    // Close the gz pipe after destroying gz_is and gz_buf (destructor order).
    gz_is.reset();
    gz_buf.reset();
    if (gz_pipe) pclose(gz_pipe);
    return reads;
}

// Auto-detect mode from file extension and a small length sample.
inline QueryMode auto_detect_mode(const std::string& path,
                                   const std::vector<RawRead>& sample) {
    if (is_fastq(path)) {
        if (!sample.empty()) {
            std::vector<size_t> lens;
            for (const auto& r : sample) lens.push_back(r.seq.size());
            std::sort(lens.begin(), lens.end());
            const size_t med = lens[lens.size() / 2];
            return (med >= 300) ? QueryMode::LONG_READS : QueryMode::SHORT_READS;
        }
        return QueryMode::LONG_READS;
    }
    // FASTA: use median length heuristic
    if (!sample.empty()) {
        std::vector<size_t> lens;
        for (const auto& r : sample) lens.push_back(r.seq.size());
        std::sort(lens.begin(), lens.end());
        const size_t med = lens[lens.size() / 2];
        if (med < 300)  return QueryMode::SHORT_READS;
        if (med < 2000) return QueryMode::LONG_READS;
    }
    return QueryMode::ASSEMBLY;
}

inline size_t median_read_length(const std::vector<RawRead>& sample) {
    if (sample.empty()) return 0;
    std::vector<size_t> lens;
    lens.reserve(sample.size());
    for (const auto& r : sample) lens.push_back(r.seq.size());
    std::sort(lens.begin(), lens.end());
    return lens[lens.size() / 2];
}

inline QueryMode sanity_check_mode_override(const std::string& path,
                                            QueryMode requested,
                                            const std::vector<RawRead>& sample,
                                            bool quiet,
                                            bool explicitRequest = true) {
    if (!explicitRequest) return requested;
    if (sample.empty()) return requested;
    const size_t med = median_read_length(sample);
    if (requested == QueryMode::LONG_READS && is_fastq(path) && med > 0 && med < 400) {
        if (!quiet) {
            std::cerr << "[query-input] mode sanity override: requested long-reads but median read length="
                      << med << "bp; using short-reads preprocessing instead\n";
        }
        return QueryMode::SHORT_READS;
    }
    if (requested == QueryMode::SHORT_READS && med >= 1500) {
        if (!quiet) {
            std::cerr << "[query-input] mode sanity override: requested short-reads but median read length="
                      << med << "bp; using long-reads preprocessing instead\n";
        }
        return QueryMode::LONG_READS;
    }
    return requested;
}

// Build CoverageReport from a filtered read set.
inline CoverageReport
compute_report(const std::vector<RawRead>& reads,
               size_t pseudoContigs, size_t dropped,
               size_t genomeSizeHint, QueryMode mode);

struct AutoTuneDecision {
    InputConfig cfg;
    CoverageReport preReport;
    std::string strategyName;
};

inline AutoTuneDecision
make_autotuned_config(QueryMode mode,
                      const InputConfig& baseCfg,
                      const std::vector<RawRead>& reads,
                      size_t dropped,
                      bool quiet) {
    AutoTuneDecision d;
    d.cfg = baseCfg;
    d.preReport = compute_report(reads, 0, dropped, baseCfg.genomeSizeHint, mode);
    d.preReport.coverageTier = classify_coverage_tier(mode, d.preReport.estimatedCov);

    switch (mode) {
        case QueryMode::ASSEMBLY:
            switch (d.preReport.coverageTier) {
                case CoverageTier::LOW:
                    d.strategyName = "assembly_low_coverage";
                    break;
                case CoverageTier::HIGH:
                    d.strategyName = "assembly_high_coverage";
                    break;
                case CoverageTier::NORMAL:
                    d.strategyName = "assembly_standard";
                    break;
                case CoverageTier::UNKNOWN:
                    d.strategyName = "assembly_unknown_coverage";
                    break;
            }
            break;

        case QueryMode::LONG_READS:
            switch (d.preReport.coverageTier) {
                case CoverageTier::LOW:
                    d.strategyName = "long_reads_low_coverage";
                    d.cfg.lrAnchorK = std::max(9, baseCfg.lrAnchorK - 2);
                    d.cfg.lrMinCluster = 1;
                    break;
                case CoverageTier::HIGH:
                    d.strategyName = "long_reads_high_coverage";
                    d.cfg.lrAnchorK = std::min(21, baseCfg.lrAnchorK + 2);
                    d.cfg.lrMinCluster = std::max<size_t>(baseCfg.lrMinCluster, 3);
                    break;
                case CoverageTier::NORMAL:
                    d.strategyName = "long_reads_standard";
                    break;
                case CoverageTier::UNKNOWN:
                    d.strategyName = "long_reads_unknown_coverage";
                    break;
            }
            break;

        case QueryMode::SHORT_READS:
            switch (d.preReport.coverageTier) {
                case CoverageTier::LOW:
                    d.strategyName = "short_reads_low_coverage";
                    d.cfg.srK = std::max(15, baseCfg.srK - 4);
                    if (baseCfg.srMinKmerFreq == 0) d.cfg.srMinKmerFreq = 1;
                    d.cfg.srMinUnitigLen = std::max<size_t>(80, std::min<size_t>(baseCfg.srMinUnitigLen, 120));
                    break;
                case CoverageTier::HIGH:
                    d.strategyName = "short_reads_high_coverage";
                    d.cfg.srK = std::min(31, baseCfg.srK + 4);
                    if (baseCfg.srMinKmerFreq == 0) d.cfg.srMinKmerFreq = 3;
                    d.cfg.srMinUnitigLen = std::max<size_t>(baseCfg.srMinUnitigLen, 250);
                    break;
                case CoverageTier::NORMAL:
                    d.strategyName = "short_reads_standard";
                    break;
                case CoverageTier::UNKNOWN:
                    d.strategyName = "short_reads_unknown_coverage";
                    break;
            }
            break;
    }

    if (!quiet) {
        std::cerr << "[query-input] auto-tuning strategy=" << d.strategyName;
        if (d.preReport.estimatedCov > 0.0)
            std::cerr << " cov~" << std::fixed << std::setprecision(1)
                      << d.preReport.estimatedCov << "x";
        std::cerr << " tier=" << coverage_tier_name(d.preReport.coverageTier);
        if (mode == QueryMode::LONG_READS) {
            std::cerr << " lrAnchorK=" << d.cfg.lrAnchorK
                      << " lrMinCluster=" << d.cfg.lrMinCluster;
        } else if (mode == QueryMode::SHORT_READS) {
            std::cerr << " srK=" << d.cfg.srK
                      << " srMinKmerFreq=" << d.cfg.srMinKmerFreq
                      << " srMinUnitigLen=" << d.cfg.srMinUnitigLen;
        }
        std::cerr << '\n';
    }
    return d;
}

inline CoverageReport
compute_report(const std::vector<RawRead>& reads,
               size_t pseudoContigs, size_t dropped,
               size_t genomeSizeHint, QueryMode mode) {
    CoverageReport r;
    r.mode = mode; r.pseudoContigs = pseudoContigs; r.droppedReads = dropped;
    if (reads.empty()) return r;

    std::vector<size_t> lens;
    lens.reserve(reads.size());
    size_t total = 0, mn = std::numeric_limits<size_t>::max(), mx = 0;
    for (const auto& rd : reads) {
        const size_t l = rd.seq.size();
        lens.push_back(l); total += l;
        if (l < mn) mn = l;
        if (l > mx) mx = l;
    }
    r.totalReads = reads.size();
    r.totalBases = total;
    r.minLen     = mn;
    r.maxLen     = mx;
    r.meanLen    = static_cast<double>(total) / static_cast<double>(reads.size());

    // N50
    std::sort(lens.begin(), lens.end(), std::greater<size_t>());
    size_t cum = 0;
    r.n50 = static_cast<double>(lens[0]);
    for (const size_t l : lens) {
        cum += l;
        if (cum * 2 >= total) { r.n50 = static_cast<double>(l); break; }
    }
    if (genomeSizeHint > 0)
        r.estimatedCov = static_cast<double>(total) / static_cast<double>(genomeSizeHint);
    r.coverageTier = classify_coverage_tier(mode, r.estimatedCov);
    return r;
}

// Drop reads outside [minLen, maxLen].
inline std::vector<RawRead>
filter_by_length(std::vector<RawRead> reads,
                 size_t minLen, size_t maxLen, size_t& dropped) {
    dropped = 0;
    std::vector<RawRead> out;
    out.reserve(reads.size());
    for (auto& r : reads) {
        if (r.seq.size() < minLen || r.seq.size() > maxLen) { ++dropped; continue; }
        out.push_back(std::move(r));
    }
    return out;
}

// FNV-1a hash for anchor k-mers.
inline uint64_t fnv1a(const char* data, size_t len) {
    uint64_t h = 14695981039346656037ULL;
    for (size_t i = 0; i < len; ++i) {
        h ^= static_cast<uint8_t>(data[i]);
        h *= 1099511628211ULL;
    }
    return h;
}

inline std::vector<RawRead>
coverage_downsample_reads(std::vector<RawRead> reads,
                          QueryMode mode,
                          size_t genomeSizeHint,
                          bool quiet,
                          size_t& dropped) {
    if (reads.empty() || genomeSizeHint == 0) return reads;

    size_t totalBases = 0;
    for (const auto& r : reads) totalBases += r.seq.size();
    if (totalBases == 0) return reads;

    const double coverage =
        static_cast<double>(totalBases) / static_cast<double>(genomeSizeHint);
    double targetCoverage = 0.0;
    switch (mode) {
        case QueryMode::ASSEMBLY:
            return reads;
        case QueryMode::LONG_READS:
            targetCoverage = 45.0;
            break;
        case QueryMode::SHORT_READS:
            targetCoverage = 80.0;
            break;
    }
    if (coverage <= targetCoverage * 1.20) return reads;

    const double keepFrac = std::min(1.0, targetCoverage / coverage);
    const uint64_t denom = 1000003ULL;
    const uint64_t threshold = static_cast<uint64_t>(
        std::floor(keepFrac * static_cast<double>(denom)));

    std::vector<RawRead> out;
    out.reserve(std::max<size_t>(1, static_cast<size_t>(
        std::ceil(static_cast<double>(reads.size()) * keepFrac))));
    for (auto& read : reads) {
        uint64_t h = fnv1a(read.name.data(), read.name.size());
        h ^= fnv1a(read.seq.data(), std::min<size_t>(read.seq.size(), 48)) + 0x9e3779b97f4a7c15ULL;
        if ((h % denom) <= threshold) out.push_back(std::move(read));
    }
    if (out.empty()) {
        size_t bestIdx = 0;
        for (size_t i = 1; i < reads.size(); ++i)
            if (reads[i].seq.size() > reads[bestIdx].seq.size()) bestIdx = i;
        out.push_back(std::move(reads[bestIdx]));
    }
    dropped += reads.size() - out.size();
    if (!quiet) {
        std::cerr << "[query-input] coverage downsampled "
                  << reads.size() << " -> " << out.size()
                  << " reads at ~" << std::fixed << std::setprecision(1)
                  << coverage << "x target~" << targetCoverage << "x\n";
    }
    return out;
}

inline std::vector<std::string>
select_overlap_rescue_inputs(const std::vector<RawRead>& reads,
                             size_t maxInputs) {
    std::vector<std::string> out;
    if (reads.empty() || maxInputs == 0) return out;
    out.reserve(std::min(reads.size(), maxInputs));
    if (reads.size() <= maxInputs) {
        for (const auto& r : reads) out.push_back(r.seq);
        return out;
    }

    const size_t stride = std::max<size_t>(1, reads.size() / maxInputs);
    for (size_t i = 0; i < reads.size() && out.size() < maxInputs; i += stride)
        out.push_back(reads[i].seq);
    if (out.size() < maxInputs) {
        for (size_t i = reads.size(); i > 0 && out.size() < maxInputs; --i)
            out.push_back(reads[i - 1].seq);
    }
    return out;
}

inline std::string reverse_complement_copy(const std::string& seq) {
    std::string rc(seq.size(), 'N');
    for (size_t i = 0; i < seq.size(); ++i) {
        const char c = seq[seq.size() - 1 - i];
        switch (c) {
            case 'A': rc[i] = 'T'; break;
            case 'C': rc[i] = 'G'; break;
            case 'G': rc[i] = 'C'; break;
            case 'T': rc[i] = 'A'; break;
            default:  rc[i] = 'N'; break;
        }
    }
    return rc;
}

inline size_t suffix_prefix_overlap(const std::string& left,
                                    const std::string& right,
                                    size_t minOverlap,
                                    double maxMismatchFrac = 0.0) {
    const size_t maxOverlap = std::min(left.size(), right.size());
    if (maxOverlap < minOverlap) return 0;
    for (size_t ov = maxOverlap; ov >= minOverlap; --ov) {
        const size_t allowedMismatches =
            (maxMismatchFrac <= 0.0)
                ? 0
                : static_cast<size_t>(std::floor(static_cast<double>(ov) * maxMismatchFrac));
        size_t mismatches = 0;
        bool match = true;
        const size_t leftStart = left.size() - ov;
        for (size_t i = 0; i < ov; ++i) {
            if (left[leftStart + i] != right[i]) {
                if (++mismatches > allowedMismatches) {
                    match = false;
                    break;
                }
            }
        }
        if (match) return ov;
        if (ov == minOverlap) break;
    }
    return 0;
}

inline bool contains_with_orientation(const std::string& haystack,
                                      const std::string& needle) {
    if (haystack.find(needle) != std::string::npos) return true;
    const std::string rc = reverse_complement_copy(needle);
    return haystack.find(rc) != std::string::npos;
}

inline std::vector<std::string>
greedy_overlap_assemble(std::vector<std::string> seqs,
                        size_t minOverlap,
                        size_t minContigLen,
                        size_t maxInputs = 256,
                        double maxMismatchFrac = 0.0) {
    std::vector<std::string> filtered;
    filtered.reserve(seqs.size());
    for (auto& seq : seqs)
        if (seq.size() >= minContigLen)
            filtered.push_back(std::move(seq));
    if (filtered.empty()) return {};

    std::sort(filtered.begin(), filtered.end(),
              [](const std::string& a, const std::string& b) {
                  if (a.size() != b.size()) return a.size() > b.size();
                  return a < b;
              });

    std::vector<std::string> uniqueSeqs;
    uniqueSeqs.reserve(filtered.size());
    for (auto& seq : filtered) {
        bool redundant = false;
        for (const auto& keep : uniqueSeqs) {
            if (contains_with_orientation(keep, seq)) {
                redundant = true;
                break;
            }
        }
        if (!redundant) uniqueSeqs.push_back(std::move(seq));
    }
    if (uniqueSeqs.size() > maxInputs)
        uniqueSeqs.resize(maxInputs);

    while (uniqueSeqs.size() > 1) {
        size_t bestI = uniqueSeqs.size();
        size_t bestJ = uniqueSeqs.size();
        size_t bestOverlap = 0;
        std::string bestMerged;

        for (size_t i = 0; i < uniqueSeqs.size(); ++i) {
            for (size_t j = 0; j < uniqueSeqs.size(); ++j) {
                if (i == j) continue;

                const size_t ovFwd = suffix_prefix_overlap(
                    uniqueSeqs[i], uniqueSeqs[j], minOverlap, maxMismatchFrac);
                if (ovFwd > bestOverlap) {
                    bestOverlap = ovFwd;
                    bestI = i;
                    bestJ = j;
                    bestMerged = uniqueSeqs[i] + uniqueSeqs[j].substr(ovFwd);
                }

                const std::string rc = reverse_complement_copy(uniqueSeqs[j]);
                const size_t ovRev = suffix_prefix_overlap(
                    uniqueSeqs[i], rc, minOverlap, maxMismatchFrac);
                if (ovRev > bestOverlap) {
                    bestOverlap = ovRev;
                    bestI = i;
                    bestJ = j;
                    bestMerged = uniqueSeqs[i] + rc.substr(ovRev);
                }
            }
        }

        if (bestOverlap < minOverlap || bestI >= uniqueSeqs.size() || bestJ >= uniqueSeqs.size())
            break;

        if (bestI > bestJ) std::swap(bestI, bestJ);
        uniqueSeqs[bestI] = std::move(bestMerged);
        uniqueSeqs.erase(uniqueSeqs.begin() + static_cast<ptrdiff_t>(bestJ));

        for (size_t i = 0; i < uniqueSeqs.size(); ) {
            if (i != bestI && contains_with_orientation(uniqueSeqs[bestI], uniqueSeqs[i])) {
                uniqueSeqs.erase(uniqueSeqs.begin() + static_cast<ptrdiff_t>(i));
                if (i < bestI) --bestI;
            } else {
                ++i;
            }
        }
    }

    std::vector<std::string> out;
    out.reserve(uniqueSeqs.size());
    for (auto& seq : uniqueSeqs)
        if (seq.size() >= minContigLen)
            out.push_back(std::move(seq));
    std::sort(out.begin(), out.end(),
              [](const std::string& a, const std::string& b) {
                  if (a.size() != b.size()) return a.size() > b.size();
                  return a < b;
              });
    return out;
}

// ----------------------------------------------------------------
// LONG_READS -> pseudo-contigs
//
// 1. Hash non-overlapping anchor k-mers from every read.
// 2. Union-Find: reads sharing >=1 anchor (with <=200 total readers
//    to skip repetitive anchors) are merged into one cluster.
// 3. Per cluster: exact overlap assembly first; if that cannot merge the reads,
//    fall back to the older majority-vote consensus.
// 4. Single-read clusters emitted directly (preserves truly unique loci).
// ----------------------------------------------------------------
inline std::unordered_map<std::string, std::string>
long_reads_to_pseudocontigs(const std::vector<RawRead>& reads,
                             int anchorK,
                             size_t minCluster,
                             bool quiet) {
    const size_t N = reads.size();
    if (N == 0) return {};
    const size_t maxAnchorsPerRead = 96;
    const size_t maxClusterReadsForConsensus = 128;

    // Step 1: anchor -> read-index multimap
    std::unordered_map<uint64_t, std::vector<uint32_t>> anchorMap;
    anchorMap.reserve(N * 10);
    for (uint32_t ri = 0; ri < static_cast<uint32_t>(N); ++ri) {
        const std::string& seq = reads[ri].seq;
        if (static_cast<int>(seq.size()) < anchorK) continue;
        std::unordered_set<uint64_t> seen;
        const size_t anchorStep = std::max<size_t>(
            static_cast<size_t>(anchorK),
            std::max<size_t>(1, seq.size() / std::max<size_t>(1, maxAnchorsPerRead)));
        for (size_t j = 0; j + static_cast<size_t>(anchorK) <= seq.size();
             j += anchorStep) {
            const uint64_t h = fnv1a(seq.data() + j, static_cast<size_t>(anchorK));
            if (seen.insert(h).second) anchorMap[h].push_back(ri);
        }
    }

    // Step 2: Union-Find (iterative path compression to avoid stack overflow on
    // large read sets — a recursive find with 10M reads could overflow the stack)
    std::vector<uint32_t> parent(N);
    std::iota(parent.begin(), parent.end(), 0u);
    auto find = [&](uint32_t x) -> uint32_t {
        // Two-pass iterative path compression
        uint32_t root = x;
        while (parent[root] != root) root = parent[root];
        while (parent[x] != root) {
            uint32_t next = parent[x];
            parent[x] = root;
            x = next;
        }
        return root;
    };
    auto unite = [&](uint32_t a, uint32_t b) {
        a = find(a); b = find(b);
        if (a != b) parent[b] = a;
    };
    auto pair_key = [](uint32_t a, uint32_t b) -> uint64_t {
        if (a > b) std::swap(a, b);
        return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
    };
    std::unordered_map<uint64_t, uint16_t> sharedAnchorPairs;
    sharedAnchorPairs.reserve(N * 12);
    for (const auto& kv : anchorMap) {
        const auto& vec = kv.second;
        if (vec.size() < 2 || vec.size() > 200) continue;
        for (size_t i = 0; i < vec.size(); ++i) {
            for (size_t j = i + 1; j < vec.size(); ++j) {
                const uint64_t key = pair_key(vec[i], vec[j]);
                uint16_t& count = sharedAnchorPairs[key];
                if (count != std::numeric_limits<uint16_t>::max()) ++count;
            }
        }
    }
    const uint16_t minSharedAnchors = (anchorK <= 10) ? 2u : 3u;
    for (const auto& kv : sharedAnchorPairs) {
        if (kv.second < minSharedAnchors) continue;
        const uint32_t a = static_cast<uint32_t>(kv.first >> 32);
        const uint32_t b = static_cast<uint32_t>(kv.first & 0xffffffffu);
        unite(a, b);
    }

    // Step 3: group by cluster root
    std::unordered_map<uint32_t, std::vector<uint32_t>> clusters;
    for (uint32_t ri = 0; ri < static_cast<uint32_t>(N); ++ri)
        clusters[find(ri)].push_back(ri);

    // Step 4: majority-vote consensus
    static constexpr char BASES[5] = {'A','C','G','T','N'};
    std::unordered_map<std::string, std::string> out;
    size_t idx = 0;

    for (auto& [root, members] : clusters) {
        if (members.size() < minCluster) {
            if (members.size() == 1 && reads[members[0]].seq.size() >= 100)
                out["lr_pc" + std::to_string(idx++)] = reads[members[0]].seq;
            continue;
        }

        std::vector<std::string> memberSeqs;
        std::vector<uint32_t> selectedMembers = members;
        if (selectedMembers.size() > maxClusterReadsForConsensus) {
            std::partial_sort(
                selectedMembers.begin(),
                selectedMembers.begin() + static_cast<ptrdiff_t>(maxClusterReadsForConsensus),
                selectedMembers.end(),
                [&](uint32_t a, uint32_t b) {
                    if (reads[a].seq.size() != reads[b].seq.size())
                        return reads[a].seq.size() > reads[b].seq.size();
                    return a < b;
                });
            selectedMembers.resize(maxClusterReadsForConsensus);
        }
        memberSeqs.reserve(selectedMembers.size());
        size_t longestRead = 0;
        for (uint32_t ri : selectedMembers) {
            memberSeqs.push_back(reads[ri].seq);
            longestRead = std::max(longestRead, reads[ri].seq.size());
        }
        const size_t overlapMin = std::max<size_t>(
            static_cast<size_t>(std::max(16, anchorK * 3)), 40);
        auto assembled = greedy_overlap_assemble(
            memberSeqs, overlapMin, 100, 256, 0.08);
        size_t longestAssembly = 0;
        for (const auto& seq : assembled)
            longestAssembly = std::max(longestAssembly, seq.size());
        if (!assembled.empty() &&
            (assembled.size() < members.size() || longestAssembly > longestRead)) {
            for (auto& seq : assembled) {
                out["lr_pc" + std::to_string(idx++) + "_n" + std::to_string(members.size())]
                    = std::move(seq);
            }
            continue;
        }

        std::vector<size_t> lens;
        for (uint32_t ri : selectedMembers) lens.push_back(reads[ri].seq.size());
        std::sort(lens.begin(), lens.end());
        const size_t medLen = lens[lens.size() / 2];
        // FIX: medLen * 1.2 promotes size_t to double, which silently truncates for
        // very large read lengths (double has 53-bit mantissa, size_t is 64 bits).
        // Use integer arithmetic: medLen + medLen/5 = medLen * 1.2 exactly.
        const size_t capLen = medLen + medLen / 5 + 1;

        const std::array<uint32_t,5> zero_arr = {0u,0u,0u,0u,0u};
        std::vector<std::array<uint32_t, 5>> freq(capLen, zero_arr);
        for (uint32_t ri : selectedMembers) {
            const std::string& seq = reads[ri].seq;
            const size_t len = std::min(seq.size(), capLen);
            for (size_t j = 0; j < len; ++j) {
                int b = 4;
                switch (seq[j]) {
                    case 'A': b=0; break; case 'C': b=1; break;
                    case 'G': b=2; break; case 'T': b=3; break;
                    default:  b=4; break;
                }
                freq[j][static_cast<size_t>(b)]++;
            }
        }

        std::string consensus;
        consensus.reserve(capLen);
        for (size_t j = 0; j < capLen; ++j) {
            const uint32_t tot = freq[j][0]+freq[j][1]+freq[j][2]+freq[j][3]+freq[j][4];
            if (tot == 0) break;
            size_t best = 0;
            for (size_t b = 1; b < 5; ++b)
                if (freq[j][b] > freq[j][best]) best = b;
            consensus += BASES[best];
        }

        if (consensus.size() >= 100)
            out["lr_pc" + std::to_string(idx++) + "_n" + std::to_string(members.size())]
                = std::move(consensus);
    }

    if (!quiet)
        std::cerr << "[query-input][long-reads] "
                  << N << " reads -> " << clusters.size()
                  << " clusters -> " << out.size() << " pseudo-contigs\n";
    return out;
}

// ----------------------------------------------------------------
// SHORT_READS -> pseudo-contigs (de Bruijn k-mer path extension)
//
// 1. Count canonical k-mers from all reads.
// 2. Auto-tune minFreq from median k-mer frequency (floor at 2).
// 3. Build prefix -> list-of-solid-k-mers successor index.
// 4. Greedy extension: for each unvisited solid k-mer, extend forward
//    by picking the unvisited, highest-frequency successor; stop when
//    no unique unambiguous successor exists.
// ----------------------------------------------------------------
inline std::unordered_map<std::string, std::string>
short_reads_to_pseudocontigs(const std::vector<RawRead>& reads,
                              int K,
                              uint32_t minFreq,
                              size_t minUnitigLen,
                              bool quiet) {
    if (reads.empty()) return {};
    const size_t kSZ = static_cast<size_t>(K);

    // Step 1: count k-mers using an in-place sliding window to avoid O(k) per-k-mer
    // string allocation.  The canonical form (lex min of fwd/rc) is computed in-place.
    std::unordered_map<std::string, uint32_t> kmerCount;
    const size_t sampleCap = std::min(reads.size(), static_cast<size_t>(120000));

    // Pre-allocate two string buffers reused across all k-mers
    std::string fwdBuf(kSZ, 'A'), rcBuf(kSZ, 'A');

    for (size_t ri = 0; ri < sampleCap; ++ri) {
        const std::string& seq = reads[ri].seq;
        if (seq.size() < kSZ) continue;
        for (size_t j = 0; j + kSZ <= seq.size(); ++j) {
            // Fill forward buffer
            for (size_t b = 0; b < kSZ; ++b) fwdBuf[b] = seq[j + b];
            // Fill reverse-complement buffer
            for (size_t b = 0; b < kSZ; ++b) {
                const char c = seq[j + kSZ - 1 - b];
                switch (c) {
                    case 'A': rcBuf[b] = 'T'; break;
                    case 'T': rcBuf[b] = 'A'; break;
                    case 'C': rcBuf[b] = 'G'; break;
                    case 'G': rcBuf[b] = 'C'; break;
                    default:  rcBuf[b] = 'N'; break;
                }
            }
            kmerCount[(fwdBuf <= rcBuf) ? fwdBuf : rcBuf]++;
        }
    }

    // Step 2: auto-tune minFreq
    if (minFreq == 0) {
        std::vector<uint32_t> freqs;
        freqs.reserve(std::min(kmerCount.size(), static_cast<size_t>(50000)));
        for (const auto& kv : kmerCount) {
            if (freqs.size() >= 50000) break;
            freqs.push_back(kv.second);
        }
        std::sort(freqs.begin(), freqs.end());
        const uint32_t med = freqs.empty() ? 0u : freqs[freqs.size() / 2];
        const uint32_t floorFreq = (minFreq > 0 ? minFreq : 2u);
        minFreq = std::max(floorFreq, static_cast<uint32_t>(med / 20));
        if (!quiet)
            std::cerr << "[query-input][short-reads] median k-mer freq=" << med
                      << " -> minFreq=" << minFreq << '\n';
    }

    // Step 3: filter solid k-mers; cap at 3M to avoid O(N^2)
    const size_t kMaxSolid = 800000;
    std::unordered_map<std::string, uint32_t> solid;
    solid.reserve(kmerCount.size() / 2 + 1);
    for (auto& kv : kmerCount)
        if (kv.second >= minFreq) solid.emplace(std::move(kv.first), kv.second);
    kmerCount.clear();

    if (solid.size() > kMaxSolid) {
        std::vector<std::pair<uint32_t,std::string>> byFreq;
        byFreq.reserve(solid.size());
        for (const auto& kv : solid) byFreq.push_back({kv.second, kv.first});
        std::partial_sort(byFreq.begin(),
                          byFreq.begin() + static_cast<ptrdiff_t>(kMaxSolid),
                          byFreq.end(),
                          [](const auto& a, const auto& b){ return a.first > b.first; });
        solid.clear();
        for (size_t i = 0; i < kMaxSolid; ++i)
            solid.emplace(std::move(byFreq[i].second), byFreq[i].first);
        if (!quiet)
            std::cerr << "[query-input][short-reads] capped to "
                      << kMaxSolid << " solid k-mers\n";
    }
    if (solid.empty()) return {};

    // Step 4: greedy extension with correct canonical k-mer orientation handling.
    // A canonical k-mer may be stored as the reverse complement of the actual
    // sequence k-mer (when RC is lexicographically smaller).  A prefix index
    // keyed on the first (k-1) chars of the canonical form therefore misses
    // successors whose canonical form is their own RC: the RC's first (k-1)
    // chars differ from the last (k-1) chars of the current canonical k-mer.
    // Fix: enumerate all four possible next bases, compute canonical(overlap+b),
    // and look up directly in solid.  cur_actual tracks the k-mer in its
    // traversal direction (may differ from its canonical form in solid/visited).
    std::unordered_set<std::string> visited;
    visited.reserve(solid.size());
    std::unordered_map<std::string, std::string> out;
    size_t pcIdx = 0;
    const size_t maxUnitigLen = 5000000;

    for (const auto& [seed, seedFreq] : solid) {
        if (visited.count(seed)) continue;
        std::string unitig = seed;
        visited.insert(seed);
        std::string cur_actual = seed;  // actual traversal direction (may not be canonical)
        uint32_t minFreqPath = seedFreq;
        while (unitig.size() < maxUnitigLen) {
            const std::string overlap = cur_actual.substr(1);  // last (k-1) bases in traversal dir
            const std::string* bestCanon = nullptr;
            uint32_t bestF = 0;
            char bestBase = 0;
            for (const char b : {'A', 'C', 'G', 'T'}) {
                std::string cand = overlap;
                cand += b;
                std::string rcCand = reverse_complement_copy(cand);
                const std::string& canon = (cand <= rcCand) ? cand : rcCand;
                auto fi = solid.find(canon);
                if (fi == solid.end() || visited.count(canon)) continue;
                if (fi->second > bestF) {
                    bestF = fi->second;
                    bestCanon = &fi->first;
                    bestBase = b;
                }
            }
            if (!bestCanon || bestF == 0) break;
            minFreqPath = std::min(minFreqPath, bestF);
            visited.insert(*bestCanon);
            unitig += bestBase;
            cur_actual.erase(cur_actual.begin());
            cur_actual.push_back(bestBase);
        }
        if (unitig.size() >= minUnitigLen)
            out["sr_unitig" + std::to_string(pcIdx++)
                + "_len" + std::to_string(unitig.size())
                + "_mf" + std::to_string(minFreqPath)] = std::move(unitig);
    }

    if (!quiet)
        std::cerr << "[query-input][short-reads] "
                  << solid.size() << " solid k-mers -> "
                  << out.size() << " unitig pseudo-contigs\n";
    return out;
}

} // namespace detail

// ================================================================
// prepare_query — the single public entry point
//
//  path      : FASTA or FASTQ file (not .gz)
//  modeHint  : user-supplied mode; pass ASSEMBLY to trigger auto-detect
//  cfg       : per-mode parameters (built from CLI Options)
//  quiet     : suppress informational stderr
// ================================================================
inline PreparedQuery
prepare_query(const std::string& path,
              QueryMode modeHint,
              const InputConfig& cfg,
              bool autoDetect = true,
              bool quiet = false) {
    PreparedQuery result;
    result.sampleName = fs::path(path).stem().string();

    const size_t maxLoad = cfg.maxReadsPerFile > 0
                         ? cfg.maxReadsPerFile
                         : std::numeric_limits<size_t>::max();
    auto rawReads = detail::load_reads(path, maxLoad, quiet);

    if (rawReads.empty()) {
        if (!quiet)
            std::cerr << "[query-input] WARNING: no sequences read from " << path << '\n';
        result.report = detail::compute_report({}, 0, 0, cfg.genomeSizeHint, modeHint);
        return result;
    }

    if (autoDetect && modeHint == QueryMode::ASSEMBLY) {
        const size_t sz = std::min(rawReads.size(), static_cast<size_t>(200));
        std::vector<detail::RawRead> sample(
            rawReads.begin(), rawReads.begin() + static_cast<ptrdiff_t>(sz));
        const QueryMode detected = detail::auto_detect_mode(path, sample);
        if (detected != QueryMode::ASSEMBLY && !quiet)
            std::cerr << "[query-input] auto-detected mode "
                      << mode_name(detected) << " for " << path
                      << " (override with --query-mode assembly if wrong)\n";
        modeHint = detected;
    }
    {
        const size_t sz = std::min(rawReads.size(), static_cast<size_t>(200));
        std::vector<detail::RawRead> sample(
            rawReads.begin(), rawReads.begin() + static_cast<ptrdiff_t>(sz));
        modeHint = detail::sanity_check_mode_override(path, modeHint, sample, quiet, !autoDetect);
    }

    size_t dropped = 0;
    switch (modeHint) {
        case QueryMode::ASSEMBLY: {
            std::vector<detail::RawRead> statsRefs = rawReads;
            for (auto& r : rawReads)
                result.contigs.emplace(std::move(r.name), std::move(r.seq));
            result.report = detail::compute_report(
                statsRefs, result.contigs.size(), 0, cfg.genomeSizeHint, modeHint);
            switch (result.report.coverageTier) {
                case CoverageTier::LOW: result.report.strategyName = "assembly_low_coverage"; break;
                case CoverageTier::HIGH: result.report.strategyName = "assembly_high_coverage"; break;
                case CoverageTier::NORMAL: result.report.strategyName = "assembly_standard"; break;
                case CoverageTier::UNKNOWN: result.report.strategyName = "assembly_unknown_coverage"; break;
            }
            break;
        }

        case QueryMode::LONG_READS: {
            auto filtered = detail::filter_by_length(
                std::move(rawReads), cfg.lrMinReadLen, cfg.lrMaxReadLen, dropped);
            const auto tuned = detail::make_autotuned_config(
                modeHint, cfg, filtered, dropped, quiet);
            filtered = detail::coverage_downsample_reads(
                std::move(filtered), modeHint, tuned.cfg.genomeSizeHint, quiet, dropped);
            result.contigs = detail::long_reads_to_pseudocontigs(
                filtered, tuned.cfg.lrAnchorK, tuned.cfg.lrMinCluster, quiet);
            if (result.contigs.empty() && !filtered.empty()) {
                if (!quiet)
                    std::cerr << "[query-input][long-reads] no pseudo-contigs assembled; emitting "
                              << filtered.size() << " individual reads as pseudo-contigs\n";
                for (auto& r : filtered)
                    result.contigs.emplace(std::move(r.name), std::move(r.seq));
            }
            // Filter pseudo-contigs by read support (cluster size encoded in _n<N>)
            {
                const int minSupport = (tuned.preReport.coverageTier == CoverageTier::LOW) ? 2 : 3;
                size_t nDropped = 0;
                std::unordered_map<std::string, std::string> kept;
                kept.reserve(result.contigs.size());
                for (auto& [name, seq] : result.contigs) {
                    auto npos = name.rfind("_n");
                    int support = -1;
                    if (npos != std::string::npos) {
                        try { support = std::stoi(name.substr(npos + 2)); } catch (...) {}
                    }
                    if (support < 0 || support >= minSupport)
                        kept.emplace(std::move(name), std::move(seq));
                    else
                        ++nDropped;
                }
                if (!quiet && nDropped > 0)
                    std::cerr << "[query-input][long-reads] dropped " << nDropped
                              << " pseudo-contigs with read support < " << minSupport << "\n";
                result.contigs = std::move(kept);
            }
            result.report = detail::compute_report(
                filtered, result.contigs.size(), dropped, tuned.cfg.genomeSizeHint, modeHint);
            result.report.coverageTier = tuned.preReport.coverageTier;
            result.report.strategyName = tuned.strategyName;
            if (result.report.totalReads < tuned.preReport.totalReads)
                result.report.strategyName += "_subsampled";
            break;
        }

        case QueryMode::SHORT_READS: {
            auto filtered = detail::filter_by_length(
                std::move(rawReads), cfg.srMinReadLen, cfg.srMaxReadLen, dropped);
            const auto tuned = detail::make_autotuned_config(
                modeHint, cfg, filtered, dropped, quiet);
            filtered = detail::coverage_downsample_reads(
                std::move(filtered), modeHint, tuned.cfg.genomeSizeHint, quiet, dropped);
            result.contigs = detail::short_reads_to_pseudocontigs(
                filtered, tuned.cfg.srK, tuned.cfg.srMinKmerFreq, tuned.cfg.srMinUnitigLen, quiet);
            // Filter unitigs by minimum k-mer frequency support (_mf<N> in name)
            if (!result.contigs.empty()) {
                const uint32_t minMf = (tuned.preReport.coverageTier == CoverageTier::LOW) ? 2u : 3u;
                size_t nDropped = 0;
                std::unordered_map<std::string, std::string> kept;
                kept.reserve(result.contigs.size());
                for (auto& [name, seq] : result.contigs) {
                    auto mfpos = name.rfind("_mf");
                    uint32_t mf = 0;
                    if (mfpos != std::string::npos) {
                        try { mf = static_cast<uint32_t>(std::stoul(name.substr(mfpos + 3))); }
                        catch (...) {}
                    }
                    if (mf == 0 || mf >= minMf)
                        kept.emplace(std::move(name), std::move(seq));
                    else
                        ++nDropped;
                }
                if (!quiet && nDropped > 0)
                    std::cerr << "[query-input][short-reads] dropped " << nDropped
                              << " unitigs with k-mer support < " << minMf << "\n";
                result.contigs = std::move(kept);
            }
            if (result.contigs.empty() && !filtered.empty()) {
                const size_t overlapMin = std::max<size_t>(
                    static_cast<size_t>(std::max(12, tuned.cfg.srK * 2)), 40);
                const size_t overlapRescueCap =
                    (tuned.preReport.coverageTier == CoverageTier::LOW) ? 256u : 512u;
                auto seqs = detail::select_overlap_rescue_inputs(filtered, overlapRescueCap);
                auto overlapContigs = detail::greedy_overlap_assemble(
                    std::move(seqs), overlapMin, tuned.cfg.srMinUnitigLen, overlapRescueCap);
                if (!overlapContigs.empty()) {
                    if (!quiet)
                        std::cerr << "[query-input][short-reads] overlap assembly rescued "
                                  << overlapContigs.size() << " pseudo-contigs\n";
                    if (overlapContigs.size() > 1) {
                        const std::string& primary = overlapContigs.front();
                        size_t synthIdx = 0;
                        for (size_t i = 1; i < overlapContigs.size(); ++i) {
                            const std::string& extra = overlapContigs[i];
                            const size_t pos = primary.find(extra);
                            if (pos == std::string::npos) continue;
                            size_t trim = std::min<size_t>(
                                static_cast<size_t>(std::max(0, tuned.cfg.srK - 1)),
                                extra.size() / 5);
                            if (extra.size() <= trim * 2 + static_cast<size_t>(tuned.cfg.srK))
                                trim = 0;
                            const std::string dupSpan =
                                extra.substr(trim, extra.size() - trim * 2);
                            if (dupSpan.size() < static_cast<size_t>(std::max(20, tuned.cfg.srK)))
                                continue;
                            const size_t insertPos = pos + extra.size() - trim;
                            std::string synth = primary.substr(0, insertPos);
                            synth += dupSpan;
                            synth += primary.substr(insertPos);
                            result.contigs.emplace(
                                "sr_dup_rescue" + std::to_string(synthIdx++)
                                    + "_len" + std::to_string(synth.size()),
                                std::move(synth));
                        }
                    }
                    size_t idx = 0;
                    for (auto& seq : overlapContigs) {
                        result.contigs.emplace(
                            "sr_ovlp" + std::to_string(idx++) + "_len" + std::to_string(seq.size()),
                            std::move(seq));
                    }
                }
            }
            const size_t finalPseudoContigs =
                result.contigs.empty() && !filtered.empty() ? filtered.size() : result.contigs.size();
            result.report = detail::compute_report(
                filtered, finalPseudoContigs, dropped, tuned.cfg.genomeSizeHint, modeHint);
            if (result.contigs.empty() && !filtered.empty()) {
                if (!quiet)
                    std::cerr << "[query-input][short-reads] no unitigs; emitting "
                              << filtered.size() << " individual reads as pseudo-contigs\n";
                for (auto& r : filtered)
                    result.contigs.emplace(std::move(r.name), std::move(r.seq));
            }
            result.report.coverageTier = tuned.preReport.coverageTier;
            result.report.strategyName = tuned.strategyName;
            if (result.report.totalReads < tuned.preReport.totalReads)
                result.report.strategyName += "_subsampled";
            break;
        }
    }

    if (!quiet)
        result.report.print(std::cerr, result.sampleName);
    return result;
}

} // namespace query_input
#endif // QUERY_INPUT_HANDLER_HPP
