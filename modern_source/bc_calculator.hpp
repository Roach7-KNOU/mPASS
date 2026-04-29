#pragma once
// bc_calculator.hpp – Buriedness Count (BC) calculation
//
// The BC score of a probe (or a protein atom) is the number of protein atoms
// within a 14 Å sphere centred on that probe.  Cells whose farthest corner is
// ≤ 14 Å away are counted wholesale; cells that straddle the boundary are
// checked atom by atom.
//
// This replaces the BC logic scattered across bumpcheckAtoms2, bumpcheckAtoms3,
// checkBC, and the inline BC block inside execDbscanClusteringFirst.

#include "data_types.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>

namespace mpass {

// ─── Per-probe BC computation ────────────────────────────────────────────────

/// Result returned by calculateBuriedness.
struct BuriednessResult {
    int    buriedness{};        // protein atoms within BC_DISTANCE_CUTOFF (14 Å)
    double closestAtomDist{9999.0};
    double averageContactDist{};
    int    averageContactCount{};
};

/// Compute the BC score for a given 3-D position.
/// Uses the same split strategy as the original: cells fully inside the 14 Å
/// sphere are counted en-masse; cells that may straddle the boundary are
/// checked atom by atom.
[[nodiscard]] inline BuriednessResult calculateBuriedness(
        const Vector3&            probePos,
        const SimulationState&    state,
        double                    cellSize = 2.0) {

    BuriednessResult result;
    int cellId = cellIndexFor(probePos,
                              state.minX, state.minY, state.minZ,
                              state.gridMaxX, state.gridMaxY, state.gridMaxZ,
                              cellSize);
    if (cellId < 0) return result;

    int wholesaleCount = 0;
    int exactCount     = 0;
    double totalContactDist = 0.0;
    int    contactCount     = 0;

    for (const auto& gp : state.gridProperties) {
        int nid = neighborCellId(cellId, gp,
                                 state.gridMaxX, state.gridMaxY, state.gridMaxZ);
        if (nid < 0) continue;
        const GridCell& cell = state.cells[nid];

        if (gp.maxDist <= BC_DISTANCE_CUTOFF) {
            // All atoms in this cell are guaranteed within 14 Å
            wholesaleCount += static_cast<int>(cell.atomIndices.size());
        } else {
            // Check each atom individually
            for (int atomIdx : cell.atomIndices) {
                double d = distance(probePos, state.atoms[atomIdx].position);
                if (d <= BC_DISTANCE_CUTOFF) {
                    ++exactCount;
                }
                if (d < result.closestAtomDist)
                    result.closestAtomDist = d;
                // Average dist uses non-metal atoms within 3.8 Å
                if (d < 3.8001 &&
                    state.atoms[atomIdx].polarity != Polarity::Metal) {
                    totalContactDist += d;
                    ++contactCount;
                }
            }
        }
    }

    result.buriedness = wholesaleCount + exactCount;
    if (contactCount > 0) {
        result.averageContactDist  = totalContactDist / contactCount;
        result.averageContactCount = contactCount;
    }
    return result;
}

// ─── Probe bump-check with atom layer ────────────────────────────────────────

/// Check whether placing a probe at probePos clashes with any protein atom or
/// existing probe (within their respective minimum distance thresholds).
/// If accepted, fills in probe geometry fields.
///
/// @param newProbe   Probe to test; position must already be set.
/// @param state      Current simulation state.
/// @param probeCutoff  Minimum probe–probe distance (1.0 Å for layers 1-3,
///                     1.2 Å for outer layers).
/// @param bcCutoffActive  When true, reject probes below g_bcCutoff.
/// @return true if the probe is accepted.
[[nodiscard]] inline bool bumpCheckAndScore(
        Probe&                 newProbe,
        const SimulationState& state,
        double                 probeCutoff     = 1.0,
        bool                   bcCutoffActive  = false,
        double                 cellSize        = 2.0) {

    int cellId = cellIndexFor(newProbe.position,
                              state.minX, state.minY, state.minZ,
                              state.gridMaxX, state.gridMaxY, state.gridMaxZ,
                              cellSize);
    if (cellId < 0) return false;

    double minDist       = 9999.0;
    double totalDist     = 0.0;
    int    contactCount  = 0;
    int    wholesaleBC   = 0;
    int    exactBC       = 0;

    for (const auto& gp : state.gridProperties) {
        int nid = neighborCellId(cellId, gp,
                                 state.gridMaxX, state.gridMaxY, state.gridMaxZ);
        if (nid < 0) continue;
        const GridCell& cell = state.cells[nid];

        // ── Probe–probe bump check ──────────────────────────────────────────
        if (gp.minDist <= 2.5 && !cell.probeIndices.empty()) {
            for (int pid : cell.probeIndices) {
                double d = distance(newProbe.position, state.probes[pid].position);
                if (d < probeCutoff - 1e-4)
                    return false;  // too close to existing probe
            }
        }

        // ── Atom checks ────────────────────────────────────────────────────
        for (int atomIdx : cell.atomIndices) {
            double d = distance(newProbe.position, state.atoms[atomIdx].position);

            // Clash with protein atom
            double minAllowed = minProbeAtomDist(state.atoms[atomIdx].polarity);
            if (d < minAllowed - 1e-4)
                return false;

            if (d < minDist) minDist = d;

            // Soft contact radius (non-metal atoms within 3.8 Å)
            if (d < 3.8001 && state.atoms[atomIdx].polarity != Polarity::Metal) {
                totalDist += d;
                ++contactCount;
            }

            // BC counting
            if (gp.maxDist <= BC_DISTANCE_CUTOFF) {
                // Covered wholesale below
            } else {
                if (d <= BC_DISTANCE_CUTOFF) ++exactBC;
            }
        }

        // Wholesale BC count
        if (gp.maxDist <= BC_DISTANCE_CUTOFF) {
            wholesaleBC += static_cast<int>(cell.atomIndices.size());
        }
    }

    newProbe.closestAtomDist   = minDist;
    newProbe.averageAtomDist   = (contactCount > 0) ? totalDist / contactCount : 0.0;
    newProbe.averageAtomCount  = contactCount;
    newProbe.buriedness        = wholesaleBC + exactBC;

    if (bcCutoffActive && newProbe.buriedness < state.bcCutoff)
        return false;

    return true;
}

// ─── Protein-atom BC ("proteinBc") ───────────────────────────────────────────

/// Compute the BC value for every protein atom and store it in
/// ProteinAtom::buriedness.  Used to weight surface accessibility.
inline void computeAtomBuriedness(SimulationState& state, double cellSize = 2.0) {
    for (auto& atom : state.atoms) {
        auto res = calculateBuriedness(atom.position, state, cellSize);
        atom.buriedness = res.buriedness;
    }
}

// ─── BC threshold derivation ─────────────────────────────────────────────────

/// After the first-layer probes have been scored, derive bcCutoff from
/// the probe BC distribution exactly as in execDbscanClusteringFirst.
inline void deriveBcCutoff(SimulationState& state) {
    if (state.probeBcList.empty()) return;

    auto& bcList = state.probeBcList;
    std::sort(bcList.begin(), bcList.end(), std::greater<int>());

    // Top-10 mean
    int sampleSize = std::min(10, static_cast<int>(bcList.size()));
    double sum = std::accumulate(bcList.begin(), bcList.begin() + sampleSize, 0.0);
    state.top10Mean = sum / sampleSize;

    // 25th-percentile quartile
    state.quartile = bcList[static_cast<std::size_t>(bcList.size() * 0.25)];

    // Linear interpolation of the BC cutoff ratio from protein size
    double upper, lower;
    const double upper_ratio = 0.5, lower_ratio = 0.8;
    std::size_t n = state.atoms.size();

    upper = 637.0 + (-895.18674 * std::exp(-0.00231 * n));
    lower = (n <= 999) ? 314.0 : 132.86 * std::log(static_cast<double>(n)) - 601.0;

    double alpha = (lower_ratio - upper_ratio) / (lower - upper);
    double beta  = upper_ratio - alpha * upper;

    double ratio = alpha * state.top10Mean + beta;
    ratio = std::clamp(ratio, upper_ratio, lower_ratio);

    if (state.bcCutoffRatioOverride > 0.0)
        ratio = state.bcCutoffRatioOverride;

    state.bcCutoffRatio = ratio;
    state.bcCutoff      = state.top10Mean * ratio;
}

} // namespace mpass
