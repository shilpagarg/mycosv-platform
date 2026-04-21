#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "layer3_routing_index.hpp"

namespace fs = std::filesystem;

struct Args {
    size_t nCentroids = 1000000;
    size_t hashesPerCentroid = 32;
    size_t queries = 32;
    size_t topK = 4;
    size_t phylumCount = 8;
    size_t chunkRecords = 4096;
    uint64_t seed = 1;
    std::string storePath = "million_scale_routing_store.bin";
    std::string reportTsv = "million_scale_routing_report.tsv";
};

static void usage(const char* argv0) {
    std::cerr
        << "Usage: " << argv0 << " [options]\n"
        << "  --n-centroids INT\n"
        << "  --hashes-per-centroid INT\n"
        << "  --queries INT\n"
        << "  --top-k INT\n"
        << "  --phylum-count INT\n"
        << "  --chunk-records INT\n"
        << "  --seed INT\n"
        << "  --store PATH\n"
        << "  --report-tsv PATH\n";
}

static std::string need_arg(const char* flag, int argc, char** argv, int& i) {
    if (i + 1 >= argc) throw std::runtime_error(std::string("Missing value for ") + flag);
    return argv[++i];
}

static Args parse_args(int argc, char** argv) {
    Args args;
    for (int i = 1; i < argc; ++i) {
        const std::string flag = argv[i];
        if (flag == "-h" || flag == "--help") {
            usage(argv[0]);
            std::exit(0);
        } else if (flag == "--n-centroids") {
            args.nCentroids = static_cast<size_t>(std::stoull(need_arg(flag.c_str(), argc, argv, i)));
        } else if (flag == "--hashes-per-centroid") {
            args.hashesPerCentroid = static_cast<size_t>(std::stoull(need_arg(flag.c_str(), argc, argv, i)));
        } else if (flag == "--queries") {
            args.queries = static_cast<size_t>(std::stoull(need_arg(flag.c_str(), argc, argv, i)));
        } else if (flag == "--top-k") {
            args.topK = static_cast<size_t>(std::stoull(need_arg(flag.c_str(), argc, argv, i)));
        } else if (flag == "--phylum-count") {
            args.phylumCount = static_cast<size_t>(std::stoull(need_arg(flag.c_str(), argc, argv, i)));
        } else if (flag == "--chunk-records") {
            args.chunkRecords = static_cast<size_t>(std::stoull(need_arg(flag.c_str(), argc, argv, i)));
        } else if (flag == "--seed") {
            args.seed = static_cast<uint64_t>(std::stoull(need_arg(flag.c_str(), argc, argv, i)));
        } else if (flag == "--store") {
            args.storePath = need_arg(flag.c_str(), argc, argv, i);
        } else if (flag == "--report-tsv") {
            args.reportTsv = need_arg(flag.c_str(), argc, argv, i);
        } else {
            throw std::runtime_error("Unknown argument: " + flag);
        }
    }
    args.nCentroids = std::max<size_t>(1, args.nCentroids);
    args.hashesPerCentroid = std::max<size_t>(1, args.hashesPerCentroid);
    args.queries = std::max<size_t>(1, args.queries);
    args.topK = std::max<size_t>(1, args.topK);
    args.phylumCount = std::max<size_t>(1, args.phylumCount);
    args.chunkRecords = std::max<size_t>(1, args.chunkRecords);
    return args;
}

static uint64_t splitmix64(uint64_t& x) {
    x += 0x9e3779b97f4a7c15ULL;
    uint64_t z = x;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static std::vector<uint64_t> make_hashes(size_t n, uint64_t seed) {
    std::vector<uint64_t> hashes;
    hashes.reserve(n);
    while (hashes.size() < n) {
        hashes.push_back(splitmix64(seed));
    }
    std::sort(hashes.begin(), hashes.end());
    hashes.erase(std::unique(hashes.begin(), hashes.end()), hashes.end());
    while (hashes.size() < n) {
        hashes.push_back(splitmix64(seed));
        std::sort(hashes.begin(), hashes.end());
        hashes.erase(std::unique(hashes.begin(), hashes.end()), hashes.end());
    }
    if (hashes.size() > n) hashes.resize(n);
    return hashes;
}

static tol::ExternalMemoryCentroidStore::DiskRecord
make_record(size_t idx, size_t hashesPerCentroid, size_t phylumCount, uint64_t seed) {
    static const char* kRanks[] = {"phylum", "class", "order", "family", "genus", "species"};
    tol::ExternalMemoryCentroidStore::DiskRecord r;
    r.cladeName = "clade_" + std::to_string(idx);
    r.phylum = "phylum_" + std::to_string(idx % phylumCount);
    r.cladeRank = kRanks[idx % 6];
    uint64_t localSeed = seed ^ (static_cast<uint64_t>(idx) * 0x9e3779b97f4a7c15ULL);
    r.hashes = make_hashes(hashesPerCentroid, localSeed);
    return r;
}

int main(int argc, char** argv) {
    const Args args = parse_args(argc, argv);
    const fs::path storePath(args.storePath);
    const fs::path reportPath(args.reportTsv);
    if (!storePath.parent_path().empty()) fs::create_directories(storePath.parent_path());
    if (!reportPath.parent_path().empty()) fs::create_directories(reportPath.parent_path());

    const auto writeStart = std::chrono::steady_clock::now();
    {
        std::ofstream out(args.storePath, std::ios::binary);
        if (!out) throw std::runtime_error("Cannot open store path: " + args.storePath);
        const uint64_t n = static_cast<uint64_t>(args.nCentroids);
        out.write(reinterpret_cast<const char*>(&n), sizeof(n));
        for (size_t i = 0; i < args.nCentroids; ++i) {
            const auto rec = make_record(i, args.hashesPerCentroid, args.phylumCount, args.seed);
            tol::ExternalMemoryCentroidStore::append_record(out, rec);
        }
    }
    const double writeSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - writeStart).count();

    tol::ExternalMemoryCentroidStore store(args.storePath);
    const auto indexStart = std::chrono::steady_clock::now();
    store.prepare_skip_index();
    const double indexSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - indexStart).count();
    size_t topHitMatches = 0;
    const auto queryStart = std::chrono::steady_clock::now();
    const size_t stride = std::max<size_t>(1, args.nCentroids / args.queries);
    for (size_t qi = 0; qi < args.queries; ++qi) {
        const size_t idx = std::min(args.nCentroids - 1, qi * stride);
        const auto rec = make_record(idx, args.hashesPerCentroid, args.phylumCount, args.seed);
        tol::CladeCentroid q;
        q.cladeName = "query_" + std::to_string(qi);
        q.phylum = rec.phylum;
        q.cladeRank = rec.cladeRank;
        q.centroidHashes = rec.hashes;
        const auto best = store.query_topk_streaming(q, args.topK, args.chunkRecords);
        if (!best.empty() && best.front().cladeName == rec.cladeName) ++topHitMatches;
    }
    const double querySeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - queryStart).count();

    const uintmax_t storeBytes = fs::exists(args.storePath) ? fs::file_size(args.storePath) : 0;
    const uintmax_t skipIndexBytes =
        fs::exists(args.storePath + ".skip") ? fs::file_size(args.storePath + ".skip") : 0;
    const double qps = querySeconds > 0.0 ? static_cast<double>(args.queries) / querySeconds : 0.0;
    const double topHitRecall = static_cast<double>(topHitMatches) / static_cast<double>(args.queries);

    std::ofstream report(args.reportTsv);
    if (!report) throw std::runtime_error("Cannot open report path: " + args.reportTsv);
    report << "n_centroids\thashes_per_centroid\tqueries\ttop_k\tphylum_count\tchunk_records\tstore_bytes\tskip_index_bytes\twrite_seconds\tindex_seconds\tquery_seconds\tqueries_per_second\ttop_hit_recall\n";
    report << args.nCentroids << '\t'
           << args.hashesPerCentroid << '\t'
           << args.queries << '\t'
           << args.topK << '\t'
           << args.phylumCount << '\t'
           << args.chunkRecords << '\t'
           << storeBytes << '\t'
           << skipIndexBytes << '\t'
           << writeSeconds << '\t'
           << indexSeconds << '\t'
           << querySeconds << '\t'
           << qps << '\t'
           << topHitRecall << '\n';

    std::cout << "n_centroids\t" << args.nCentroids << "\n"
              << "store_bytes\t" << storeBytes << "\n"
              << "skip_index_bytes\t" << skipIndexBytes << "\n"
              << "write_seconds\t" << writeSeconds << "\n"
              << "index_seconds\t" << indexSeconds << "\n"
              << "query_seconds\t" << querySeconds << "\n"
              << "queries_per_second\t" << qps << "\n"
              << "top_hit_recall\t" << topHitRecall << "\n";
    return 0;
}
