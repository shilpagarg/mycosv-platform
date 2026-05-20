#ifndef FUNGI_TOL_HIERARCHICAL_ENGINE_HPP
#define FUNGI_TOL_HIERARCHICAL_ENGINE_HPP

// hierarchical_engine.hpp — Off-reference novelty scoring and variant stamping.
//
// Responsibilities:
//   • OffRefNoveltyTier enum + helpers (score_off_ref_novelty, novelty_tier_name)
//   • make_off_reference_call_scored<V>() — stamps off-ref fields on any call record
//   • UncoveredWindow — half-open [start,end) interval of query bases with no ref match

#include <cmath>
#include <string>
#include <vector>

namespace tol {

// ── OffRefNoveltyTier ─────────────────────────────────────────────────────
// Thresholds based on k-mer Jaccard overlap fraction with the best matching
// reference contig.  Boundaries chosen empirically on fungal genome benchmarks:
//   NOVEL      : overlap == 0   → truly sequence-novel (new clade, HGT island …)
//   NOVEL_WEAK : 0 < overlap < 0.05  → very diverged, likely new species
//   DIVERGED   : 0.05–0.20          → diverged within same genus/family
//   OFF_REF_KNOWN: ≥ 0.20           → known off-reference (Starship, TE, RIP …)
enum class OffRefNoveltyTier {
    NOVEL,
    NOVEL_WEAK,
    DIVERGED,
    OFF_REF_KNOWN,
};

inline const char* novelty_tier_name(OffRefNoveltyTier t) {
    switch (t) {
        case OffRefNoveltyTier::NOVEL:         return "NOVEL";
        case OffRefNoveltyTier::NOVEL_WEAK:    return "NOVEL_WEAK";
        case OffRefNoveltyTier::DIVERGED:      return "DIVERGED";
        default:                               return "OFF_REF_KNOWN";
    }
}

inline OffRefNoveltyTier score_off_ref_novelty(double overlapFraction) {
    if (overlapFraction <= 0.0)  return OffRefNoveltyTier::NOVEL;
    if (overlapFraction < 0.05)  return OffRefNoveltyTier::NOVEL_WEAK;
    if (overlapFraction < 0.20)  return OffRefNoveltyTier::DIVERGED;
    return OffRefNoveltyTier::OFF_REF_KNOWN;
}

// Cross-clade novelty: a region absent from same-clade references but present
// in a different clade is a candidate HGT event.  Returns NOVEL so the calling
// site treats the locus as highly interesting; the elementClass field should be
// set to "HGT" separately via classify_repeat_element().
//
// Thresholds widened: same-clade < 0.10 (was 0.05) and other-clade ≥ 0.08
// (was 0.10), AND the other-clade signal must exceed same-clade by ≥ 0.05.
// The previous strict AND-gate combined with the k=7 Jaccard inflation meant
// real HGT islands bordered by host sequence almost never reached NOVEL —
// the calling sites now use higher k and containment, so the thresholds
// can be loosened without inflating FPs.
inline OffRefNoveltyTier score_cross_clade_novelty(
        double sameCladeOverlap,
        double highestOtherCladeOverlap) {
    const bool absentInSameClade  = sameCladeOverlap < 0.10;
    const bool presentInOther     = highestOtherCladeOverlap >= 0.08;
    const bool otherExceedsSelf   =
        (highestOtherCladeOverlap - sameCladeOverlap) >= 0.05;
    if (absentInSameClade && presentInOther && otherExceedsSelf)
        return OffRefNoveltyTier::NOVEL;
    return score_off_ref_novelty(sameCladeOverlap);
}

// ── UncoveredWindow ───────────────────────────────────────────────────────
// Half-open interval [start, end) in query coordinates with no reference MEMs.
// Collected during MEM-chain construction; later used to classify
// REPEAT / TE / STARSHIP / HGT / RIP annotation classes.
struct UncoveredWindow {
    size_t start = 0;
    size_t end   = 0;
    size_t length() const { return end > start ? end - start : 0; }
};

// ── make_off_reference_call_scored ───────────────────────────────────────
// Stamps off-reference fields onto any call record V that has the same
// fields as VariantCallBridge.  Template so it works with both the bridge
// struct and future specialisations.
//
// Sets:
//   annotation        ← novelty tier name
//   type              ← "OFF_REF"
//   pantreeClass      ← "NON_REF"
//   isNonRefVariant   ← true
//   triallelicTopology← "." if currently empty
template <class Variant>
inline Variant make_off_reference_call_scored(Variant v, OffRefNoveltyTier tier) {
    v.annotation       = novelty_tier_name(tier);
    v.type             = "OFF_REF";
    v.pantreeClass     = "NON_REF";
    v.isNonRefVariant  = true;
    if (v.triallelicTopology.empty()) v.triallelicTopology = ".";
    return v;
}



// ── Probabilistic evidence fusion ────────────────────────────────────────
// A depth-aware fusion model that combines read- and assembly-level evidence
// on the log-likelihood scale. The implementation is deliberately lightweight
// and header-only so it can be used in harnesses without any external math
// dependency.

enum class EvidenceLayer {
    ASSEMBLY,
    LONG_READ,
    SHORT_READ,
};

inline const char* evidence_layer_name(EvidenceLayer l) {
    switch (l) {
        case EvidenceLayer::ASSEMBLY:   return "assembly";
        case EvidenceLayer::LONG_READ:  return "long_read";
        case EvidenceLayer::SHORT_READ: return "short_read";
    }
    return "assembly";
}

struct EvidenceObservation {
    EvidenceLayer layer            = EvidenceLayer::ASSEMBLY;
    double        logLikelihoodRef = 0.0;
    double        logLikelihoodAlt = 0.0;
    double        depth            = 1.0;
    double        mapq             = 60.0;
    double        breakpointSupport= 0.0;  // [0,1]
    double        spanSupport      = 0.0;  // [0,1]
};

struct FusedEvidenceScore {
    double posteriorAlt   = 0.5;
    double posteriorRef   = 0.5;
    double logOddsAlt     = 0.0;
    double effectiveDepth = 0.0;
    size_t layersUsed     = 0;

    bool supports_variant(double minPosterior = 0.90) const {
        return posteriorAlt >= minPosterior;
    }
};

inline double clamp_unit(double x) {
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    return x;
}

inline double stable_sigmoid(double x) {
    if (x >= 0.0) {
        const double z = std::exp(-x);
        return 1.0 / (1.0 + z);
    }
    const double z = std::exp(x);
    return z / (1.0 + z);
}

inline double evidence_layer_reliability(const EvidenceObservation& obs) {
    const double mapqW = clamp_unit(obs.mapq / 60.0);
    const double support = 0.5 * clamp_unit(obs.breakpointSupport)
                         + 0.5 * clamp_unit(obs.spanSupport);
    switch (obs.layer) {
        case EvidenceLayer::ASSEMBLY:
            return 0.70 + 0.30 * std::max(mapqW, support);
        case EvidenceLayer::LONG_READ:
            return 0.55 + 0.45 * (0.60 * mapqW + 0.40 * support);
        case EvidenceLayer::SHORT_READ:
            return 0.45 + 0.55 * (0.60 * mapqW + 0.40 * support);
    }
    return 0.5;
}

inline FusedEvidenceScore
fuse_probabilistic_evidence(const std::vector<EvidenceObservation>& evidence,
                            double priorAlt = 0.5) {
    FusedEvidenceScore out;
    if (evidence.empty()) return out;

    priorAlt = std::max(1e-9, std::min(1.0 - 1e-9, priorAlt));
    double logOdds = std::log(priorAlt / (1.0 - priorAlt));
    double effDepth = 0.0;

    for (const auto& obs : evidence) {
        const double reliability = evidence_layer_reliability(obs);
        // Saturating depth term: additional depth helps at all depths but with
        // diminishing returns, which prevents ultra-deep short-read piles from
        // swamping higher-level assembly evidence.
        const double depthWeight = std::log1p(std::max(0.0, obs.depth));
        const double llr = (obs.logLikelihoodAlt - obs.logLikelihoodRef);
        logOdds += reliability * depthWeight * llr;
        effDepth += reliability * std::max(0.0, obs.depth);
        ++out.layersUsed;
    }

    out.logOddsAlt     = logOdds;
    out.posteriorAlt   = stable_sigmoid(logOdds);
    out.posteriorRef   = 1.0 - out.posteriorAlt;
    out.effectiveDepth = effDepth;
    return out;
}
} // namespace tol

#endif // FUNGI_TOL_HIERARCHICAL_ENGINE_HPP
