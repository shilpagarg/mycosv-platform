#ifndef FUNGI_TOL_TAXONOMY_RANKS_HPP
#define FUNGI_TOL_TAXONOMY_RANKS_HPP

// taxonomy_ranks.hpp — Single authoritative vocabulary for Linnaean ranks.
//
// ALL rank strings used by manifest readers, routing shards, and hier_rows
// iteration pass through rank_from_string / rank_to_string so a typo is
// caught in one place at compile-time or parse-time.
//
// C++ mirror of the Python RANKS list in test_amf.py; both must be kept in sync.

#include <string>

namespace tol {

// ── TaxonomyRank ──────────────────────────────────────────────────────────
enum class TaxonomyRank {
    phylum,
    class_rank,   // stored on disk as "class" (reserved C++ word)
    order,
    family,
    genus,
    species,
    unknown,      // catch-all for unrecognised strings
};

// ── rank_to_string ────────────────────────────────────────────────────────
inline const char* rank_to_string(TaxonomyRank r) {
    switch (r) {
        case TaxonomyRank::phylum:     return "phylum";
        case TaxonomyRank::class_rank: return "class";
        case TaxonomyRank::order:      return "order";
        case TaxonomyRank::family:     return "family";
        case TaxonomyRank::genus:      return "genus";
        case TaxonomyRank::species:    return "species";
        default:                       return "unknown";
    }
}

// ── rank_from_string ──────────────────────────────────────────────────────
// Returns TaxonomyRank::unknown for unrecognised strings so callers can warn.
inline TaxonomyRank rank_from_string(const std::string& s) {
    if (s == "phylum")  return TaxonomyRank::phylum;
    if (s == "class")   return TaxonomyRank::class_rank;
    if (s == "order")   return TaxonomyRank::order;
    if (s == "family")  return TaxonomyRank::family;
    if (s == "genus")   return TaxonomyRank::genus;
    if (s == "species") return TaxonomyRank::species;
    return TaxonomyRank::unknown;
}

// Ordered broadest-to-narrowest.  Defined once; used by hier_rows iteration
// and the multi-rank index builder.
inline constexpr TaxonomyRank kLinnaeanRanks[] = {
    TaxonomyRank::phylum,
    TaxonomyRank::class_rank,
    TaxonomyRank::order,
    TaxonomyRank::family,
    TaxonomyRank::genus,
    TaxonomyRank::species,
};
inline constexpr int kLinnaeanRankCount = 6;

} // namespace tol

#endif // FUNGI_TOL_TAXONOMY_RANKS_HPP
