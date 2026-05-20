#pragma once
// te_classifier.hpp — k-mer nearest-centroid TE classifier
//
// Parses PanTEon/RepBase-format FASTA headers:
//   >ID#Class/Order/Superfamily  (e.g. >TE1#DNA/TIR/Tc1-Mariner)
//   >ID#Class/Superfamily        (e.g. >TE2#LTR/Copia)
//   >ID  (unlabeled — skipped during training)
//
// Builds per-superfamily CladeCentroid objects using FracMin sketching,
// then organises them in a VPTree for O(log N) nearest-centroid lookup.
// Predictions are returned at three levels: class, order, superfamily.
//
// Build the classifier:
//   TEClassifier clf;
//   clf.train("pangeon_train.fasta", "te_index/");
//
// Classify sequences:
//   clf.load("te_index/");
//   auto res = clf.classify_fasta("test.fasta", "predictions.tsv");
//
// CLI (wired into main.cpp via --te-train / --te-classify):
//   mycosv --te-train  --query-list train.fa --out-prefix te_index/clf
//   mycosv --te-classify --query-list test.fa --te-index-prefix te_index/clf
//          --out-prefix te_index/clf
//
// Dependencies: layer3_routing_index.hpp (CladeCentroid, VPTree)

#include "layer3_routing_index.hpp"
#include "fungi_tol_bridge.hpp"  // FastaStream for gz-transparent FASTA reading

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace te {

// =========================================================================
// Label parsing
// =========================================================================
struct TELabel {
    std::string id;
    std::string te_class;    // e.g. "DNA", "LTR", "LINE", "SINE", "RC"
    std::string te_order;    // e.g. "TIR", "Gypsy/Ty3", "L1", ""
    std::string superfamily; // e.g. "Tc1-Mariner", "Copia"
    bool labeled = false;
};

// Normalise common RepBase/PanTEon spelling variants.
inline std::string norm_class(const std::string& s) {
    if (s == "Class_I" || s == "Retrotransposon") return "LTR";
    if (s == "Class_II" || s == "TIR") return "DNA";
    if (s == "MITE") return "DNA";
    return s;
}

inline TELabel parse_label(const std::string& header) {
    TELabel lbl;
    // header may start with '>'
    const std::string h = (header[0] == '>') ? header.substr(1) : header;

    // ID is everything up to the first '#' or space
    const size_t hash_pos  = h.find('#');
    const size_t space_pos = h.find(' ');
    const size_t id_end    = std::min(hash_pos, space_pos);
    lbl.id = h.substr(0, id_end);

    if (hash_pos == std::string::npos || hash_pos >= h.size() - 1) {
        lbl.labeled = false;
        return lbl;
    }

    // Everything after '#'
    std::string tax = h.substr(hash_pos + 1);
    // Strip any trailing whitespace / description after space
    const size_t sp = tax.find(' ');
    if (sp != std::string::npos) tax = tax.substr(0, sp);

    // Split on '/'
    std::vector<std::string> parts;
    std::stringstream ss(tax);
    std::string tok;
    while (std::getline(ss, tok, '/')) {
        if (!tok.empty()) parts.push_back(tok);
    }

    if (parts.empty()) { lbl.labeled = false; return lbl; }

    lbl.te_class    = norm_class(parts[0]);
    lbl.te_order    = (parts.size() >= 3) ? parts[1] : "";
    lbl.superfamily = parts.back();
    lbl.labeled     = true;
    return lbl;
}

// Combine class/order/superfamily into a single centroid name (for VPTree lookup)
inline std::string centroid_key(const TELabel& lbl) {
    if (!lbl.te_order.empty())
        return lbl.te_class + "/" + lbl.te_order + "/" + lbl.superfamily;
    return lbl.te_class + "/" + lbl.superfamily;
}

// =========================================================================
// K-mer hashing (FNV-1a 64-bit, canonical k-mer)
// TE classification uses FracMin sketching rather than syncmer seeding:
// FracMin keeps all k-mer hashes ≤ p·2^64 which gives a denser, density-
// controlled sketch suited to short TE elements (~500 bp), whereas the
// Hong-Buhler syncmer in BaseBlockSegmenter::syncmers is optimised for
// long-genome routing centroids.
// =========================================================================
inline uint64_t fnv1a_64(const char* data, size_t len) {
    uint64_t h = 14695981039346656037ULL;
    for (size_t i = 0; i < len; ++i) {
        h ^= static_cast<uint64_t>(static_cast<unsigned char>(data[i]));
        h *= 1099511628211ULL;
    }
    return h;
}

inline char complement(char c) {
    switch (c) {
        case 'A': case 'a': return 'T';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        case 'T': case 't': return 'A';
        default: return 'N';
    }
}

// Canonical k-mer: lexicographic min of forward and reverse-complement
inline uint64_t canonical_kmer_hash(const std::string& seq, size_t pos, int k) {
    std::string fwd(seq.begin() + static_cast<long>(pos),
                    seq.begin() + static_cast<long>(pos + static_cast<size_t>(k)));
    std::string rev(k, ' ');
    for (int i = 0; i < k; ++i)
        rev[static_cast<size_t>(k - 1 - i)] = complement(fwd[static_cast<size_t>(i)]);

    const std::string& canonical = (fwd < rev) ? fwd : rev;
    return fnv1a_64(canonical.data(), static_cast<size_t>(k));
}

// FracMin sketch: keep all k-mer hashes <= threshold (bottom fraction p)
inline std::vector<uint64_t> fracmin_sketch(const std::string& seq, int k, double p = 0.05) {
    const uint64_t thresh =
        static_cast<uint64_t>(p * static_cast<double>(std::numeric_limits<uint64_t>::max()));
    std::vector<uint64_t> hashes;
    const size_t n = seq.size();
    if (n < static_cast<size_t>(k)) return hashes;

    // Deduplicate with unordered_set before returning
    std::unordered_set<uint64_t> seen;
    for (size_t i = 0; i + static_cast<size_t>(k) <= n; ++i) {
        const uint64_t h = canonical_kmer_hash(seq, i, k);
        if (h <= thresh && seen.insert(h).second)
            hashes.push_back(h);
    }
    std::sort(hashes.begin(), hashes.end());
    return hashes;
}

// =========================================================================
// Build a CladeCentroid from a collection of sequences for one TE label
// =========================================================================
inline tol::CladeCentroid build_centroid(const std::string& label,
                                          const std::vector<std::string>& seqs,
                                          int k = 21,
                                          double fracmin_p = 0.05,
                                          size_t max_hashes = 4096) {
    tol::CladeCentroid::StreamBuilder sb;
    sb.name      = label;
    sb.phylum    = "TE";
    sb.cladeRank = "superfamily";
    sb.maxH      = max_hashes;

    for (const auto& seq : seqs) {
        auto hashes = fracmin_sketch(seq, k, fracmin_p);
        sb.accumulate(hashes);
    }
    return sb.finalize();
}

// =========================================================================
// FASTA reader — yields (header, sequence) pairs
// =========================================================================
inline void read_fasta(const std::string& path,
                        const std::function<void(const std::string&, const std::string&)>& cb) {
    FastaStream fs(path);
    std::istream& in = fs.get();
    std::string line, hdr, seq;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!hdr.empty() && !seq.empty()) cb(hdr, seq);
            hdr = line;
            seq.clear();
        } else {
            seq += line;
        }
    }
    if (!hdr.empty() && !seq.empty()) cb(hdr, seq);
}

// =========================================================================
// TEClassifier
// =========================================================================
class TEClassifier {
public:
    struct Params {
        int    k          = 21;
        double fracmin_p  = 0.05;
        size_t max_hashes = 4096;
    };

    struct Prediction {
        std::string id;
        std::string pred_class;
        std::string pred_order;
        std::string pred_superfamily;
        double      jaccard_sim  = 0.0;  // 1 - Jaccard distance to best centroid
        std::string best_centroid_key;   // full "Class/Order/Superfamily" key
    };

    TEClassifier() : params_(Params{}) {}
    explicit TEClassifier(Params p) : params_(p) {}

    // ---- Training -------------------------------------------------------
    // Parse labeled FASTA, accumulate sequences per superfamily, build VPTree.
    void train(const std::string& fasta_path, bool verbose = true) {
        std::unordered_map<std::string, std::vector<std::string>> bucket; // key → seqs

        size_t total = 0, labeled = 0;
        read_fasta(fasta_path, [&](const std::string& hdr, const std::string& seq) {
            ++total;
            const TELabel lbl = parse_label(hdr);
            if (!lbl.labeled) return;
            ++labeled;
            bucket[centroid_key(lbl)].push_back(seq);
            // Remember decomposition for later
            key_to_label_[centroid_key(lbl)] = lbl;
        });

        if (verbose) {
            std::cerr << "[te-train] sequences: " << total
                      << "  labeled: " << labeled
                      << "  superfamilies: " << bucket.size() << '\n';
        }

        std::vector<tol::CladeCentroid> centroids;
        centroids.reserve(bucket.size());
        for (const auto& [key, seqs] : bucket) {
            centroids.push_back(
                build_centroid(key, seqs, params_.k, params_.fracmin_p, params_.max_hashes));
            if (verbose) {
                std::cerr << "  centroid " << key
                          << "  seqs=" << seqs.size()
                          << "  hashes=" << centroids.back().centroidHashes.size() << '\n';
            }
        }
        tree_.build(centroids);
        trained_ = true;
    }

    // ---- Save / load ----------------------------------------------------
    void save(const std::string& index_prefix) const {
        if (!trained_) throw std::runtime_error("TEClassifier: not trained");
        tree_.save(index_prefix + ".vptree");

        std::ofstream meta(index_prefix + ".meta");
        meta << "k=" << params_.k << '\n';
        meta << "fracmin_p=" << params_.fracmin_p << '\n';
        meta << "max_hashes=" << params_.max_hashes << '\n';
        // Write centroid key → label decomposition
        for (const auto& [key, lbl] : key_to_label_) {
            meta << "L\t" << key << '\t' << lbl.te_class << '\t'
                 << lbl.te_order << '\t' << lbl.superfamily << '\n';
        }
    }

    void load(const std::string& index_prefix) {
        tree_ = tol::VPTree::load(index_prefix + ".vptree");
        if (tree_.empty()) throw std::runtime_error("TEClassifier: failed to load VPTree");

        std::ifstream meta(index_prefix + ".meta");
        if (!meta) throw std::runtime_error("TEClassifier: cannot open " + index_prefix + ".meta");
        std::string line;
        while (std::getline(meta, line)) {
            if (line.substr(0, 2) == "k=")
                params_.k = std::stoi(line.substr(2));
            else if (line.substr(0, 10) == "fracmin_p=")
                params_.fracmin_p = std::stod(line.substr(10));
            else if (line.substr(0, 11) == "max_hashes=")
                params_.max_hashes = std::stoul(line.substr(11));
            else if (!line.empty() && line[0] == 'L') {
                // L\tkey\tclass\torder\tsuperfamily
                std::istringstream ss(line.substr(2));
                std::string key, cls, ord, sf;
                std::getline(ss, key, '\t');
                std::getline(ss, cls, '\t');
                std::getline(ss, ord, '\t');
                std::getline(ss, sf,  '\t');
                TELabel lbl;
                lbl.te_class = cls; lbl.te_order = ord; lbl.superfamily = sf; lbl.labeled = true;
                key_to_label_[key] = lbl;
            }
        }
        trained_ = true;
    }

    // ---- Classify one sequence ------------------------------------------
    Prediction classify(const std::string& id, const std::string& seq) const {
        if (!trained_) throw std::runtime_error("TEClassifier: not trained or loaded");
        Prediction p;
        p.id = id;
        if (seq.size() < static_cast<size_t>(params_.k)) return p;

        // Build query centroid
        tol::CladeCentroid qc;
        qc.cladeName = id;
        qc.centroidHashes = fracmin_sketch(seq, params_.k, params_.fracmin_p);
        qc.build_prefilters();

        auto hits = tree_.query_topk(qc, 1);
        if (hits.empty()) return p;

        p.best_centroid_key = hits[0].cladeName;
        p.jaccard_sim       = hits[0].jaccard;

        auto it = key_to_label_.find(p.best_centroid_key);
        if (it != key_to_label_.end()) {
            p.pred_class       = it->second.te_class;
            p.pred_order       = it->second.te_order;
            p.pred_superfamily = it->second.superfamily;
        } else {
            // Key is the centroid name itself — parse it
            const TELabel lbl = parse_label("#" + p.best_centroid_key);
            p.pred_class       = lbl.te_class;
            p.pred_order       = lbl.te_order;
            p.pred_superfamily = lbl.superfamily;
        }
        return p;
    }

    // ---- Classify entire FASTA, write TSV -------------------------------
    std::vector<Prediction> classify_fasta(const std::string& fasta_path,
                                            const std::string& out_tsv) const {
        std::vector<Prediction> results;
        std::ofstream out(out_tsv);
        if (!out) throw std::runtime_error("Cannot open output TSV: " + out_tsv);

        out << "id\tpred_class\tpred_order\tpred_superfamily\tjaccard_sim\tbest_centroid\n";

        read_fasta(fasta_path, [&](const std::string& hdr, const std::string& seq) {
            std::string id_full = (hdr[0] == '>') ? hdr.substr(1) : hdr;
            const size_t h_pos = id_full.find('#');
            const size_t s_pos = id_full.find(' ');
            const size_t id_end = std::min(
                h_pos != std::string::npos ? h_pos : id_full.size(),
                s_pos != std::string::npos ? s_pos : id_full.size());
            const std::string id = id_full.substr(0, id_end);
            auto pred = classify(id, seq);
            results.push_back(pred);
            out << pred.id << '\t'
                << pred.pred_class << '\t'
                << pred.pred_order << '\t'
                << pred.pred_superfamily << '\t'
                << std::fixed << std::setprecision(4) << pred.jaccard_sim << '\t'
                << pred.best_centroid_key << '\n';
        });
        return results;
    }

    bool trained() const { return trained_; }
    size_t num_centroids() const { return tree_.size(); }

private:
    Params     params_;
    tol::VPTree tree_;
    bool        trained_ = false;
    std::unordered_map<std::string, TELabel> key_to_label_;
};

} // namespace te
