#pragma once
// probe_placer.hpp – Probe placement: triple-sphere intersection + bump check
//
// Replaces GetTripleVertex() and the probe-generation loop logic in cal.cpp.
//
// Algorithm (GetTripleVertex):
//   Given three spheres (centres Ri, Rj, Rk with combined radii sgmi, sgmj,
//   sgmk) and a probe radius sgmp, find the two points equidistant from all
//   three centres at the sum distance (sgm + sgmp).  The derivation is
//   geometric: build a local coordinate frame (Xu = Rj-Ri normalised, Zu =
//   normal to the plane), find the projection point Rb in the plane, then
//   lift ±h out of the plane to get the two solutions.

#include "data_types.hpp"
#include "geometry.hpp"
#include <optional>
#include <array>

namespace mpass {

// ─── Triple-vertex result ────────────────────────────────────────────────────

/// Both candidate probe positions arising from three touching spheres.
struct TripleVertexResult {
    Vector3 point1;
    Vector3 point2;
};

// ─── GetTripleVertex (geometry engine) ───────────────────────────────────────

/// Compute the two positions where a probe of radius `probeRadius` can touch
/// all three atoms simultaneously.
///
/// Returns std::nullopt when:
///   - the three atoms are coplanar with no real lift height (discriminant < 0)
///   - the result fails the ≤ 0.01 Å tolerance check (degenerate geometry)
///
/// @param centerI  Centre of atom i
/// @param centerJ  Centre of atom j
/// @param centerK  Centre of atom k
/// @param radiusI  VDW radius of atom i
/// @param radiusJ  VDW radius of atom j
/// @param radiusK  VDW radius of atom k
/// @param probeRadius  Radius of the probe sphere
[[nodiscard]] inline std::optional<TripleVertexResult>
getTripleVertex(const Vector3& centerI, const Vector3& centerJ, const Vector3& centerK,
                double radiusI, double radiusJ, double radiusK,
                double probeRadius) noexcept {

    // ── Local coordinate frame ─────────────────────────────────────────────
    Vector3 ijVec = centerJ - centerI;
    if (ijVec.lengthSquared() < 1e-15) return std::nullopt; // coincident atoms
    Vector3 xu = ijVec.normalized();
    Vector3 temp = centerK - centerI;

    // Zu = normal to the plane spanned by (I, J, K)
    Vector3 zu = xu.cross(temp);
    if (zu.lengthSquared() < 1e-15) return std::nullopt; // degenerate (collinear)
    zu = zu.normalized();

    Vector3 yu = zu.cross(xu);

    // ── Projection points (midpoint + correction along each edge) ─────────
    auto projectionMidpoint = [&](const Vector3& a, const Vector3& b,
                                  double ra, double rb) -> Vector3 {
        double distSq  = distanceSquared(a, b);
        double t       = (sq(ra + probeRadius) - sq(rb + probeRadius)) / (2.0 * distSq);
        return (a + b) * 0.5 + (b - a) * t;
    };

    Vector3 tij = projectionMidpoint(centerI, centerJ, radiusI, radiusJ);
    Vector3 tik = projectionMidpoint(centerI, centerK, radiusI, radiusK);

    // ── Find base point Rb in the plane ────────────────────────────────────
    // Solve: Rb = Tij + t * Yu  where Rb lies on the perpendicular bisector
    // plane from tik.
    Vector3 tikMtij = tik - tij;
    Vector3 tikMri  = tik - centerI;
    double  denom   = tikMri.dot(yu);
    if (std::abs(denom) < 1e-15) return std::nullopt;

    double  t  = tikMtij.dot(tikMri) / denom;
    Vector3 rb = tij + yu * t;

    // ── Lift height h out of the plane ─────────────────────────────────────
    double hSq = sq(radiusI + probeRadius) - distanceSquared(rb, centerI);
    if (hSq < 0.0) return std::nullopt;
    double h = std::sqrt(hSq);

    TripleVertexResult result;
    result.point1 = rb + zu * h;
    result.point2 = rb - zu * h;

    // ── Validation: all three contact distances must be within 0.01 Å ─────
    auto check = [&](const Vector3& p) -> bool {
        return std::abs(distance(centerI, p) - radiusI - probeRadius) <= 0.01 &&
               std::abs(distance(centerJ, p) - radiusJ - probeRadius) <= 0.01 &&
               std::abs(distance(centerK, p) - radiusK - probeRadius) <= 0.01;
    };
    if (!check(result.point1)) return std::nullopt;

    return result;
}

// ─── Convenience wrapper for atom indices ────────────────────────────────────

/// Same as getTripleVertex but takes atom indices from the simulation state.
[[nodiscard]] inline std::optional<TripleVertexResult>
probePositionsFromThreeAtoms(int i, int j, int k,
                             double probeRadius,
                             const SimulationState& state) noexcept {
    return getTripleVertex(
        state.atoms[i].position, state.atoms[j].position, state.atoms[k].position,
        state.atoms[i].vdwRadius, state.atoms[j].vdwRadius, state.atoms[k].vdwRadius,
        probeRadius);
}

// ─── Probe neighbour assignment ───────────────────────────────────────────────

/// Populate Probe::neighborProbes for a batch of probe indices.
/// Probes within `cutoffSq` squared distance are recorded as neighbours.
///
/// @param probeIndices  Indices of probes to process (empty = all).
/// @param cutoffSq      Squared distance threshold (e.g. 1.5² = 2.25 for DBSCAN).
/// @param cellSize      Grid cell size used for lookup.
inline void assignProbeNeighbors(SimulationState& state,
                                 const std::vector<int>& probeIndices,
                                 double cutoffSq,
                                 int    maxGridOffset, // number of GridProperty entries
                                 double cellSize = 2.0) {
    for (int pi : probeIndices) {
        Probe& probe = state.probes[pi];
        probe.neighborProbes.clear();

        int cellId = cellIndexFor(probe.position,
                                  state.minX, state.minY, state.minZ,
                                  state.gridMaxX, state.gridMaxY, state.gridMaxZ,
                                  cellSize);
        if (cellId < 0) continue;

        int limit = std::min(maxGridOffset,
                             static_cast<int>(state.gridProperties.size()));
        for (int k = 0; k < limit; ++k) {
            int nid = neighborCellId(cellId, state.gridProperties[k],
                                     state.gridMaxX, state.gridMaxY, state.gridMaxZ);
            if (nid < 0) continue;
            for (int npi : state.cells[nid].probeIndices) {
                if (npi == pi) continue;
                double dsq = distanceSquared(probe.position,
                                             state.probes[npi].position);
                if (dsq > 0 && dsq < cutoffSq)
                    probe.neighborProbes.push_back(npi);
            }
        }
    }
}

} // namespace mpass
