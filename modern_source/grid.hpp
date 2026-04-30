#pragma once
// grid.hpp – Spatial grid / cell-list infrastructure
//
// Replaces the raw index arithmetic spread throughout cal.cpp:
//   getCellId / getCellIdXyz / initializeCellList
//
// A uniform 3-D grid is built once from the protein bounding box.
// Each grid cell stores the indices of protein atoms and probes that fall
// inside it.  Neighbor-cell traversal uses the pre-loaded GridProperty
// offsets, exactly as in the original code, but the indexing logic is now
// encapsulated in one place.

#include "data_types.hpp"
#include <cmath>
#include <vector>
#include <algorithm>

namespace mpass {

// ─── Cell-index helpers ──────────────────────────────────────────────────────

/// Return the 1-D grid cell index for a point given the grid origin and
/// per-dimension size.  Returns -1 when the point is outside the grid.
[[nodiscard]] inline int cellIndexFor(const Vector3& point,
                                      double minX, double minY, double minZ,
                                      int gridMaxX, int gridMaxY, int gridMaxZ,
                                      double cellSize) noexcept {
    auto ix = static_cast<int>(std::ceil((point.x - minX) / cellSize));
    auto iy = static_cast<int>(std::ceil((point.y - minY) / cellSize));
    auto iz = static_cast<int>(std::ceil((point.z - minZ) / cellSize));
    if (ix < 0 || iy < 0 || iz < 0 ||
        ix >= gridMaxX || iy >= gridMaxY || iz >= gridMaxZ)
        return -1;
    return ix + gridMaxX * iy + gridMaxX * gridMaxY * iz;
}

[[nodiscard]] inline int cellIndexForCoords(double x, double y, double z,
                                            double minX, double minY, double minZ,
                                            int gridMaxX, int gridMaxY, int gridMaxZ,
                                            double cellSize) noexcept {
    return cellIndexFor({x, y, z}, minX, minY, minZ,
                        gridMaxX, gridMaxY, gridMaxZ, cellSize);
}

/// Compute neighbor cell id given a base cell id and a relative offset stored
/// in a GridProperty entry.  Returns -1 when the resulting id is out of range.
[[nodiscard]] inline int neighborCellId(int baseCellId,
                                        const GridProperty& gp,
                                        int gridMaxX, int gridMaxY, int gridMaxZ) noexcept {
    int newId = baseCellId
              + gp.gridIndex[0]
              + gridMaxX * gp.gridIndex[1]
              + gridMaxX * gridMaxY * gp.gridIndex[2];
    int totalCells = gridMaxX * gridMaxY * gridMaxZ;
    return (newId >= 0 && newId < totalCells) ? newId : -1;
}

// ─── Grid initialisation ────────────────────────────────────────────────────

/// Build the uniform 3-D grid from the protein bounding box and insert all
/// protein atoms into the appropriate cells.  Mirrors initializeCellList()
/// and sets SimulationState::gridMaxX/Y/Z plus SimulationState::cells.
///
/// @param state   Simulation state with atoms and bounding-box fields already
///                populated by the PDB reader.
/// @param cellSize Size of one grid cell in Angstroms (default 2.0 Å).
/// @return Total number of near-atom pairs found (atoms within 7.7 Å of each
///         other) – matches the original return value from initializeCellList.
inline int buildCellList(SimulationState& state, double cellSize = 2.0) {
    // Grid dimensions from bounding box
    state.gridMaxX = static_cast<int>(std::ceil((state.maxX - state.minX) / cellSize)) + 1;
    state.gridMaxY = static_cast<int>(std::ceil((state.maxY - state.minY) / cellSize)) + 1;
    state.gridMaxZ = static_cast<int>(std::ceil((state.maxZ - state.minZ) / cellSize)) + 1;

    int totalCells = state.gridMaxX * state.gridMaxY * state.gridMaxZ;
    state.cells.assign(totalCells, GridCell{});
    state.atomCellList.clear();

    // Insert atoms into cells
    for (int i = 0; i < static_cast<int>(state.atoms.size()); ++i) {
        int cid = cellIndexFor(state.atoms[i].position,
                               state.minX, state.minY, state.minZ,
                               state.gridMaxX, state.gridMaxY, state.gridMaxZ,
                               cellSize);
        if (cid < 0) continue;
        state.cells[cid].atomIndices.push_back(i);
        state.atomCellList.push_back(cid);
    }

    // Deduplicate atomCellList
    std::sort(state.atomCellList.begin(), state.atomCellList.end());
    state.atomCellList.erase(
        std::unique(state.atomCellList.begin(), state.atomCellList.end()),
        state.atomCellList.end());

    // Count near-atom pairs (d ≤ 7.7 Å) – original return value
    const int nearPairRangeLimit = 275; // original used 275 grid-property entries
    int pairCount = 0;
    const double pairCutoff = 7.7;

    for (int i = 0; i < static_cast<int>(state.atoms.size()); ++i) {
        int cid = cellIndexFor(state.atoms[i].position,
                               state.minX, state.minY, state.minZ,
                               state.gridMaxX, state.gridMaxY, state.gridMaxZ,
                               cellSize);
        if (cid < 0) continue;

        int limit = std::min(nearPairRangeLimit,
                             static_cast<int>(state.gridProperties.size()));
        for (int k = 0; k < limit; ++k) {
            int nid = neighborCellId(cid, state.gridProperties[k],
                                     state.gridMaxX, state.gridMaxY, state.gridMaxZ);
            if (nid < 0) continue;
            for (int j : state.cells[nid].atomIndices) {
                if (j <= i) continue; // count each pair once
                if (distance(state.atoms[i].position, state.atoms[j].position)
                        <= pairCutoff)
                    ++pairCount;
            }
        }
    }
    return pairCount;
}

/// Insert a newly created probe into the grid.
inline void insertProbeIntoGrid(SimulationState& state, int probeIndex,
                                double cellSize = 2.0) {
    int cid = cellIndexFor(state.probes[probeIndex].position,
                           state.minX, state.minY, state.minZ,
                           state.gridMaxX, state.gridMaxY, state.gridMaxZ,
                           cellSize);
    if (cid < 0) return;
    state.cells[cid].probeIndices.push_back(probeIndex);
    state.probeCellList.push_back(cid);
}

/// Deduplicate probeCellList (call after a batch of insertions).
inline void deduplicateProbeCellList(SimulationState& state) {
    std::sort(state.probeCellList.begin(), state.probeCellList.end());
    state.probeCellList.erase(
        std::unique(state.probeCellList.begin(), state.probeCellList.end()),
        state.probeCellList.end());
}

/// Rebuild probeCellList and grid probe-index lists from state.probes.
/// Call this after probe compaction (weedOverlappingProbes) to ensure
/// grid cells only reference valid, re-indexed probe indices.
inline void rebuildProbeCellList(SimulationState& state, double cellSize = 2.0) {
    for (auto& cell : state.cells)
        cell.probeIndices.clear();
    state.probeCellList.clear();
    for (int i = 0; i < static_cast<int>(state.probes.size()); ++i)
        insertProbeIntoGrid(state, i, cellSize);
    deduplicateProbeCellList(state);
}

/// Populate ProteinAtom::neighborAtoms for all atoms (within 7.7 Å of each
/// other).  Must be called after buildCellList().
inline void buildAtomNeighborList(SimulationState& state, double cellSize = 2.0) {
    const double cutoffSq = 7.7 * 7.7;
    for (int i = 0; i < static_cast<int>(state.atoms.size()); ++i) {
        state.atoms[i].neighborAtoms.clear();
        int cid = cellIndexFor(state.atoms[i].position,
                               state.minX, state.minY, state.minZ,
                               state.gridMaxX, state.gridMaxY, state.gridMaxZ,
                               cellSize);
        if (cid < 0) continue;
        for (const auto& gp : state.gridProperties) {
            int nid = neighborCellId(cid, gp,
                                     state.gridMaxX, state.gridMaxY, state.gridMaxZ);
            if (nid < 0) continue;
            for (int j : state.cells[nid].atomIndices) {
                if (j == i) continue;
                if (distanceSquared(state.atoms[i].position,
                                    state.atoms[j].position) <= cutoffSq)
                    state.atoms[i].neighborAtoms.push_back(j);
            }
        }
        auto& na = state.atoms[i].neighborAtoms;
        std::sort(na.begin(), na.end());
        na.erase(std::unique(na.begin(), na.end()), na.end());
    }
}

/// Populate ProteinAtom::neighborProbes for all atoms.
/// Finds probes within `cutoff` Å of each atom (scans up to `maxGridOffset`
/// grid-property entries, matching the original assignAtomProbe behaviour).
/// Must be called after probes are in the grid (insertProbeIntoGrid / rebuild).
inline void buildAtomProbeNeighborList(SimulationState& state,
                                       double cutoff         = 5.2,
                                       int    maxGridOffset  = 275,
                                       double cellSize       = 2.0) {
    const double cutoffSq = cutoff * cutoff;
    for (auto& atom : state.atoms)
        atom.neighborProbes.clear();

    for (int i = 0; i < static_cast<int>(state.atoms.size()); ++i) {
        ProteinAtom& atom = state.atoms[i];
        int cid = cellIndexFor(atom.position,
                               state.minX, state.minY, state.minZ,
                               state.gridMaxX, state.gridMaxY, state.gridMaxZ,
                               cellSize);
        if (cid < 0) continue;
        int limit = std::min(maxGridOffset,
                             static_cast<int>(state.gridProperties.size()));
        for (int k = 0; k < limit; ++k) {
            int nid = neighborCellId(cid, state.gridProperties[k],
                                     state.gridMaxX, state.gridMaxY, state.gridMaxZ);
            if (nid < 0) continue;
            for (int probeIdx : state.cells[nid].probeIndices) {
                if (distanceSquared(atom.position,
                                    state.probes[probeIdx].position) <= cutoffSq)
                    atom.neighborProbes.push_back(probeIdx);
            }
        }
        auto& np = atom.neighborProbes;
        std::sort(np.begin(), np.end());
        np.erase(std::unique(np.begin(), np.end()), np.end());
    }
}

} // namespace mpass
