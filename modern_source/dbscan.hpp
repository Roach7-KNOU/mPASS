#pragma once
// dbscan.hpp – DBSCAN clustering for probe sets
//
// Replaces dbscanClustering() and dbscanClusteringForSurface() in cal.cpp.
//
// Algorithm (standard DBSCAN):
//   1. For each probe compute density = |{j : dist(i,j) < eps, j in probeIds}|
//      using the pre-built nearProbes lists.
//   2. Visit each unvisited probe.  If density >= minPts start a new cluster
//      and expand it by BFS through all reachable core points.
//   3. Probes that never reach minPts neighbours remain with clusterId == 0
//      (noise / unassigned).
//
// Two variants are provided:
//   - clusterProbes()        – writes Probe::clusterId (matches dbscanClustering)
//   - clusterProbesSurface() – writes Probe::clusterIdFirst (matches
//                              dbscanClusteringForSurface)

#include "data_types.hpp"
#include "geometry.hpp"
#include <vector>
#include <algorithm>

namespace mpass {

// ─── Internal BFS expander ───────────────────────────────────────────────────

namespace detail {

/// Expand a DBSCAN cluster starting from probe `seed`.
/// @param seed       Index of the core probe that starts the cluster.
/// @param clusterId  The cluster label to assign.
/// @param minPts     Minimum-density threshold.
/// @param state      Simulation state.
/// @param getClusterId   Lambda: (const Probe&) -> int  to read cluster-id field.
/// @param setClusterId   Lambda: (Probe&, int) -> void to write cluster-id field.
template <typename GetId, typename SetId>
inline void expandCluster(int seed, int clusterId, int minPts,
                          SimulationState& state,
                          GetId getClusterId, SetId setClusterId) {
    // Working queue = mergedProbes already computed for `seed`
    // (includes all probes within eps distance)
    std::vector<int> queue = state.probes[seed].mergedProbes;

    for (std::size_t qi = 0; qi < queue.size(); ++qi) {
        int k = queue[qi];
        if (state.probes[k].isVisited == 0) {
            state.probes[k].isVisited = 1;
            if (state.probes[k].density >= minPts) {
                // Add k's neighbours into the queue for further expansion
                for (int nb : state.probes[k].mergedProbes)
                    queue.push_back(nb);
            }
        }
        // Assign cluster if still unassigned
        if (getClusterId(state.probes[k]) == 0)
            setClusterId(state.probes[k], clusterId);
    }
}

} // namespace detail

// ─── Main clustering function ────────────────────────────────────────────────

/// Cluster `probeIds` using DBSCAN(eps, minPts).
/// Writes results to Probe::clusterId.
///
/// @param probeIds   Subset of probe indices to cluster (may be all probes).
/// @param eps        Epsilon neighbourhood distance.
/// @param minPts     Minimum points to form a core.
/// @param startId    First cluster-ID value (incremented before each new cluster).
/// @return           Last cluster-ID assigned (== startId when nothing clustered).
[[nodiscard]] inline int clusterProbes(SimulationState& state,
                                       const std::vector<int>& probeIds,
                                       double eps, int minPts,
                                       int    startId = 0) {
    // ── Step 1: Density calculation ────────────────────────────────────────
    for (int pi : probeIds) {
        Probe& probe = state.probes[pi];
        probe.isVisited = 0;
        probe.clusterId = 0;
        probe.density   = 0;
        probe.mergedProbes.clear();

        for (int nb : probe.neighborProbes) {
            if (nb == pi) continue;
            // Restrict to probeIds membership (fast path: nearProbes already
            // filtered to the right distance range by assignProbeNeighbors)
            if (std::find(probeIds.begin(), probeIds.end(), nb) == probeIds.end())
                continue;
            double d = distance(probe.position, state.probes[nb].position);
            if (d < eps) {
                ++probe.density;
                probe.mergedProbes.push_back(nb);
            }
        }
    }

    // ── Step 2: Cluster expansion ──────────────────────────────────────────
    int clusterId = startId;
    for (int pi : probeIds) {
        Probe& probe = state.probes[pi];
        if (probe.isVisited != 0) continue;
        probe.isVisited = 1;

        if (probe.density < minPts) continue; // noise – leave clusterId = 0

        ++clusterId;
        probe.clusterId = clusterId;
        detail::expandCluster(
            pi, clusterId, minPts, state,
            [](const Probe& p) { return p.clusterId; },
            [](Probe& p, int id) { p.clusterId = id; });
    }
    return clusterId;
}

/// Variant that writes to Probe::clusterIdFirst (for surface / opening detection).
[[nodiscard]] inline int clusterProbesSurface(SimulationState& state,
                                              const std::vector<int>& probeIds,
                                              double eps, int minPts) {
    // Density using ALL pairs within probeIds (O(n²) – probeIds is small here)
    for (int pi : probeIds) {
        Probe& probe = state.probes[pi];
        probe.isVisited     = 0;
        probe.clusterIdFirst = 0;
        probe.density       = 0;
        probe.mergedProbes.clear();

        for (int pj : probeIds) {
            if (pj == pi) continue;
            double d = distance(probe.position, state.probes[pj].position);
            if (d < eps) {
                ++probe.density;
                probe.mergedProbes.push_back(pj);
            }
        }
    }

    int clusterId = 0;
    for (int pi : probeIds) {
        Probe& probe = state.probes[pi];
        if (probe.isVisited != 0) continue;
        probe.isVisited = 1;
        if (probe.density < minPts) continue;

        ++clusterId;
        probe.clusterIdFirst = clusterId;
        detail::expandCluster(
            pi, clusterId, minPts, state,
            [](const Probe& p) { return p.clusterIdFirst; },
            [](Probe& p, int id) { p.clusterIdFirst = id; });
    }
    return clusterId;
}

// ─── Cluster-set builder ─────────────────────────────────────────────────────

/// After clustering, group probe indices by their Probe::clusterId and
/// populate state.clusters.
inline void buildClusterSet(SimulationState& state) {
    std::map<int, Cluster> tmpMap;

    for (int i = 0; i < static_cast<int>(state.probes.size()); ++i) {
        int cid = state.probes[i].clusterId;
        if (cid <= 0) continue;

        auto [it, inserted] = tmpMap.emplace(cid, Cluster{});
        Cluster& cl = it->second;
        if (inserted) {
            cl.clusterId     = cid;
            cl.minBuriedness = 9999;
        }

        cl.probeIds.push_back(i);
        cl.sumBuriedness += state.probes[i].buriedness;
        ++cl.probeCount;
        cl.maxLayer = std::max(cl.maxLayer, state.probes[i].numLayer);

        if (state.probes[i].buriedness > cl.maxBuriedness) {
            cl.maxBuriedness     = state.probes[i].buriedness;
            cl.maxBuriedProbeId  = i;
        }
        if (state.probes[i].buriedness < cl.minBuriedness)
            cl.minBuriedness = state.probes[i].buriedness;

        switch (state.probes[i].numLayer) {
            case 1: ++cl.layer1; break;
            case 2: ++cl.layer2; break;
            case 3: ++cl.layer3; break;
            default: ++cl.layer4; break;
        }
    }

    state.clusters.clear();
    state.clusters.reserve(tmpMap.size());
    for (auto& [_, cl] : tmpMap) {
        if (cl.probeCount > 0)
            cl.averageBuriedness = static_cast<double>(cl.sumBuriedness)
                                   / cl.probeCount;
        state.clusters.push_back(std::move(cl));
    }
}

// ─── First-layer weed-by-density ─────────────────────────────────────────────

/// Remove redundant probes that are too close to each other (< 0.7 Å),
/// keeping the one with higher density (or smaller average distance on ties).
/// Matches weedFirstLayer() / checkProbeBump() behaviour.
inline void weedOverlappingProbes(SimulationState& state, double weedRadius = 0.7) {
    // For each surviving probe, mark nearby probes as removed
    std::sort(state.probes.begin(), state.probes.end(),
              [](const Probe& a, const Probe& b) { return a.density > b.density; });

    // Re-index after sort – probes' own `index` fields need updating
    for (int i = 0; i < static_cast<int>(state.probes.size()); ++i)
        state.probes[i].index = i;

    for (int i = 0; i < static_cast<int>(state.probes.size()); ++i) {
        if (!state.probes[i].isSurvived) continue;

        for (int nb : state.probes[i].neighborProbes) {
            if (nb == i || nb >= static_cast<int>(state.probes.size())) continue;
            Probe& other = state.probes[nb];
            if (!other.isSurvived) continue;

            double d = distance(state.probes[i].position, other.position);
            if (d < weedRadius - 1e-4) {
                // Keep the one with higher density; break ties by average dist
                if (state.probes[i].density > other.density) {
                    other.isSurvived = false;
                    other.clusterId  = -4;
                } else if (state.probes[i].density < other.density) {
                    state.probes[i].isSurvived = false;
                    state.probes[i].clusterId  = -4;
                    break; // current probe is dead
                } else {
                    if (state.probes[i].averageAtomDist <= other.averageAtomDist) {
                        other.isSurvived = false;
                        other.clusterId  = -4;
                    } else {
                        state.probes[i].isSurvived = false;
                        state.probes[i].clusterId  = -4;
                        break;
                    }
                }
            }
        }
    }

    // Compact: remove non-surviving probes
    state.probes.erase(
        std::remove_if(state.probes.begin(), state.probes.end(),
                       [](const Probe& p) { return !p.isSurvived; }),
        state.probes.end());

    // Re-index
    for (int i = 0; i < static_cast<int>(state.probes.size()); ++i)
        state.probes[i].index = i;
}

} // namespace mpass
