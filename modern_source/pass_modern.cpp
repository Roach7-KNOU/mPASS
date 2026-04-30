// pass_modern.cpp – Modern C++17 entry point for mPASS
//
// Orchestrates the same pipeline as the original pass.cpp:
//   1. Read configuration and PDB files
//   2. Build cell list (atom neighbours and BC)
//   3. Generate probe layers 1-4+
//   4. DBSCAN clustering
//   5. Score and output

#include "file_reader.hpp"
#include "grid.hpp"
#include "bc_calculator.hpp"
#include "probe_placer.hpp"
#include "dbscan.hpp"
#include "residue_properties.hpp"

#include <fstream>
#include <iostream>
#include <ctime>
#include <cstdio>
#include <string>
#include <algorithm>
#include <vector>

using namespace mpass;

// ─── Timing helper ───────────────────────────────────────────────────────────

static double elapsedSeconds(clock_t begin, clock_t end) {
    return static_cast<double>(end - begin) / CLOCKS_PER_SEC;
}

// ─── Usage ───────────────────────────────────────────────────────────────────

static void printUsage(const char* prog) {
    std::cerr << "Usage: " << prog
              << " <pdb_file> [bc_cutoff_ratio]\n"
              << "  bc_cutoff_ratio  optional override (0.0 – 1.0, default auto)\n";
}

// ─── Layer 1 surface count ────────────────────────────────────────────────────
// Counts all protein atoms within 7 Å of a probe using up to 667 grid-property
// entries.  Matches execDbscanClusteringFirst (no proteinBc filter).

static int surfaceCountLayer1(const Vector3& pos, const SimulationState& state) {
    int count = 0;
    int cellId = cellIndexFor(pos,
                              state.minX, state.minY, state.minZ,
                              state.gridMaxX, state.gridMaxY, state.gridMaxZ, 2.0);
    if (cellId < 0) return 0;
    int limit = std::min(667, static_cast<int>(state.gridProperties.size()));
    for (int k = 0; k < limit; ++k) {
        int nid = neighborCellId(cellId, state.gridProperties[k],
                                 state.gridMaxX, state.gridMaxY, state.gridMaxZ);
        if (nid < 0) continue;
        for (int atomIdx : state.cells[nid].atomIndices) {
            if (distance(pos, state.atoms[atomIdx].position) <= 7.0)
                ++count;
        }
    }
    return count;
}

// ─── Layer 1 bump check (atoms only) ─────────────────────────────────────────
// Checks atom clashes with no probe-probe or BC filter.
// Fills probe.neighborAtoms (all atoms within 5.2 Å) and geometry fields.
// Returns false if any atom is too close.

static bool bumpCheckLayer1(Probe& newProbe, const SimulationState& state) {
    int cellId = cellIndexFor(newProbe.position,
                              state.minX, state.minY, state.minZ,
                              state.gridMaxX, state.gridMaxY, state.gridMaxZ, 2.0);
    if (cellId < 0) return false;

    double minDist      = 9999.0;
    double totalDist    = 0.0;
    int    contactCount = 0;

    int limit = std::min(275, static_cast<int>(state.gridProperties.size()));
    for (int k = 0; k < limit; ++k) {
        int nid = neighborCellId(cellId, state.gridProperties[k],
                                 state.gridMaxX, state.gridMaxY, state.gridMaxZ);
        if (nid < 0) continue;
        for (int atomIdx : state.cells[nid].atomIndices) {
            double d = distance(newProbe.position, state.atoms[atomIdx].position);

            double cutoff = minProbeAtomDist(state.atoms[atomIdx].polarity);
            if (d < cutoff - 1e-4) return false;

            if (d <= 5.20)
                newProbe.neighborAtoms.push_back(atomIdx);

            if (d < 3.8001 && state.atoms[atomIdx].polarity != Polarity::Metal) {
                totalDist += d;
                ++contactCount;
            }
            if (d < minDist) minDist = d;
        }
    }

    // Deduplicate (the same atom may appear in multiple neighbour cells)
    auto& na = newProbe.neighborAtoms;
    std::sort(na.begin(), na.end());
    na.erase(std::unique(na.begin(), na.end()), na.end());

    newProbe.closestAtomDist  = minDist;
    newProbe.averageAtomDist  = (contactCount > 0) ? totalDist / contactCount : 0.0;
    newProbe.averageAtomCount = contactCount;
    return true;
}

// ─── Layer 1: atom-atom-atom probe generation ─────────────────────────────────
// Mirrors generateFirstLayer() in the original.  For every ordered triplet
// (i < j, i < k, j < k) drawn from each atom's nearAtom list, it attempts to
// place a probe tangent to all three atoms.

static void generateLayer1(SimulationState& state) {
    state.currentLayer = 1;
    state.probeCellList.clear();

    int atomCount = static_cast<int>(state.atoms.size());
    for (int i = 0; i < atomCount; ++i) {
        // Copy neighbour list to avoid invalidation if state.probes grows
        const std::vector<int> nbrs = state.atoms[i].neighborAtoms;
        int polarI  = static_cast<int>(state.atoms[i].polarity);

        for (int a = 0; a < static_cast<int>(nbrs.size()); ++a) {
            int j = nbrs[a];
            if (j <= i) continue;

            for (int b = a + 1; b < static_cast<int>(nbrs.size()); ++b) {
                int k = nbrs[b];
                if (k <= j) continue;

                int polarJ = static_cast<int>(state.atoms[j].polarity);
                int polarK = static_cast<int>(state.atoms[k].polarity);

                int  isPolarSum = polarI + polarJ + polarK;
                bool isMetal    = (polarI == 2 || polarJ == 2 || polarK == 2);

                // Distances between the three atom centres (+ epsilon)
                double distIJ = distance(state.atoms[i].position,
                                         state.atoms[j].position) + 0.0001;
                double distIK = distance(state.atoms[i].position,
                                         state.atoms[k].position) + 0.0001;
                double distJK = distance(state.atoms[j].position,
                                         state.atoms[k].position) + 0.0001;

                double iR = 0.0, jR = 0.0, kR = 0.0, probeRadius = 0.0;
                int    check_ok = 0, n = 0, n_a = 0;

                // ── Case: at least one polar/metal atom ───────────────────
                if (isPolarSum > 0) {
                    for (int step = 0; step < 8 && !check_ok; ++step) {
                        n           = step;
                        probeRadius = PROBE_RADIUS_POLAR + RADIUS_INCREMENT * n;

                        iR = (polarI == 1) ? PROBE_RADIUS_POLAR + RADIUS_INCREMENT * n
                           : (polarI == 2) ? PROBE_RADIUS_METAL + RADIUS_INCREMENT * n
                           :                 state.atoms[i].vdwRadius + 0.35
                                             - RADIUS_INCREMENT * n;
                        jR = (polarJ == 1) ? PROBE_RADIUS_POLAR + RADIUS_INCREMENT * n
                           : (polarJ == 2) ? PROBE_RADIUS_METAL + RADIUS_INCREMENT * n
                           :                 state.atoms[j].vdwRadius + 0.35
                                             - RADIUS_INCREMENT * n;
                        kR = (polarK == 1) ? PROBE_RADIUS_POLAR + RADIUS_INCREMENT * n
                           : (polarK == 2) ? PROBE_RADIUS_METAL + RADIUS_INCREMENT * n
                           :                 state.atoms[k].vdwRadius + 0.35
                                             - RADIUS_INCREMENT * n;

                        if (iR + jR + 2*probeRadius > distIJ &&
                            iR + kR + 2*probeRadius > distIK &&
                            jR + kR + 2*probeRadius > distJK)
                            check_ok = 1;
                    }
                    // Secondary try: vary only apolar radii outward
                    if (!check_ok) {
                        for (int step = 0; step < 7 && !check_ok; ++step) {
                            n_a = step;
                            if (polarI == 0)
                                iR = state.atoms[i].vdwRadius + RADIUS_INCREMENT * n_a;
                            if (polarJ == 0)
                                jR = state.atoms[j].vdwRadius + RADIUS_INCREMENT * n_a;
                            if (polarK == 0)
                                kR = state.atoms[k].vdwRadius + RADIUS_INCREMENT * n_a;

                            if (iR + jR + 2*probeRadius > distIJ &&
                                iR + kR + 2*probeRadius > distIK &&
                                jR + kR + 2*probeRadius > distJK)
                                check_ok = 1;
                        }
                    }
                }

                // ── Case: all apolar ──────────────────────────────────────
                if (isPolarSum == 0) {
                    probeRadius = PROBE_RADIUS_APOLAR;  // 1.6
                    for (int step = 0; step < 7 && !check_ok; ++step) {
                        n_a = step;
                        iR  = state.atoms[i].vdwRadius + RADIUS_INCREMENT * n_a;
                        jR  = state.atoms[j].vdwRadius + RADIUS_INCREMENT * n_a;
                        kR  = state.atoms[k].vdwRadius + RADIUS_INCREMENT * n_a;

                        if (iR + jR + 2*probeRadius > distIJ &&
                            iR + kR + 2*probeRadius > distIK &&
                            jR + kR + 2*probeRadius > distJK)
                            check_ok = 1;
                    }
                }

                if (!check_ok) continue;

                auto tv = getTripleVertex(state.atoms[i].position,
                                          state.atoms[j].position,
                                          state.atoms[k].position,
                                          iR, jR, kR, probeRadius);
                if (!tv) continue;

                for (int pos = 0; pos < 2; ++pos) {
                    Vector3 pt = (pos == 0) ? tv->point1 : tv->point2;

                    Probe probe;
                    probe.position           = pt;
                    probe.radius             = probeRadius;
                    probe.isPolar            = isPolarSum;
                    probe.numLayer           = 1;
                    probe.clusterId1         = 1;
                    probe.isSurvived         = true;
                    probe.contactAtomIndices = {i, j, k};
                    if (isMetal) probe.isPolar = 7;

                    if (bumpCheckLayer1(probe, state)) {
                        probe.index = static_cast<int>(state.probes.size());
                        state.probes.push_back(std::move(probe));
                        insertProbeIntoGrid(state,
                            static_cast<int>(state.probes.size()) - 1);
                    }
                }
            }
        }
    }

    deduplicateProbeCellList(state);
}

// ─── Layer 2: probe-atom-atom probe generation ────────────────────────────────
// Mirrors generateSubFirstLayer().  Uses each layer-1 probe's atom neighbour
// list to form (probe, atom_j, atom_k) triplets.

static void generateLayer2(SimulationState& state) {
    state.currentLayer = 2;

    // Snapshot probe count – only process the probes that already exist
    int probeCount = static_cast<int>(state.probes.size());

    for (int i = 0; i < probeCount; ++i) {
        // Copy the fields needed to avoid invalidation when state.probes grows
        const Vector3          probePos     = state.probes[i].position;
        const std::vector<int> atomNeighbors = state.probes[i].neighborAtoms;
        int atomsSize = static_cast<int>(atomNeighbors.size());

        for (int a = 0; a < atomsSize; ++a) {
            int j = atomNeighbors[a];
            for (int b = a + 1; b < atomsSize; ++b) {
                int k = atomNeighbors[b];

                int polarJ = static_cast<int>(state.atoms[j].polarity);
                int polarK = static_cast<int>(state.atoms[k].polarity);

                int  isPolarSum = polarJ + polarK;
                bool isMetal    = (polarJ == 2 || polarK == 2);

                // Distances (probe i is treated as the third "atom")
                double distJK = distance(state.atoms[j].position,
                                         state.atoms[k].position) + 0.0001;
                double distIJ = distance(probePos,
                                         state.atoms[j].position) + 0.0001;
                double distIK = distance(probePos,
                                         state.atoms[k].position) + 0.0001;

                double iR = 0.0, jR = 0.0, kR = 0.0, probeRadius = 0.0;
                int    check_ok = 0, n = 0, n_a = 0;

                if (isPolarSum > 0) {
                    for (int step = 0; step < 8 && !check_ok; ++step) {
                        n           = step;
                        probeRadius = PROBE_RADIUS + RADIUS_INCREMENT * n;
                        iR          = PROBE_RADIUS - RADIUS_INCREMENT * n;

                        jR = (polarJ == 1) ? 1.8 + RADIUS_INCREMENT * n
                           : (polarJ == 2) ? 1.5 + RADIUS_INCREMENT * n
                           :                 state.atoms[j].vdwRadius + 0.9
                                             - RADIUS_INCREMENT * n;
                        kR = (polarK == 1) ? 1.8 + RADIUS_INCREMENT * n
                           : (polarK == 2) ? 1.5 + RADIUS_INCREMENT * n
                           :                 state.atoms[k].vdwRadius + 0.9
                                             - RADIUS_INCREMENT * n;

                        if (jR + kR + 2*probeRadius > distJK &&
                            iR + jR + 2*probeRadius > distIJ &&
                            iR + kR + 2*probeRadius > distIK)
                            check_ok = 1;
                    }
                }
                if (!check_ok) {
                    probeRadius = PROBE_RADIUS;
                    for (int step = 0; step < 7 && !check_ok; ++step) {
                        n_a = step;
                        iR  = PROBE_RADIUS;
                        if (polarJ == 0)
                            jR = state.atoms[j].vdwRadius + 0.9 + RADIUS_INCREMENT * n_a;
                        if (polarK == 0)
                            kR = state.atoms[k].vdwRadius + 0.9 + RADIUS_INCREMENT * n_a;

                        if (jR + kR + 2*probeRadius > distJK &&
                            iR + jR + 2*probeRadius > distIJ &&
                            iR + kR + 2*probeRadius > distIK)
                            check_ok = 1;
                    }
                }

                if (!check_ok) continue;

                // Order: atom j, atom k, probe i  (matches original call)
                auto tv = getTripleVertex(state.atoms[j].position,
                                          state.atoms[k].position,
                                          probePos,
                                          jR, kR, iR, probeRadius);
                if (!tv) continue;

                for (int pos = 0; pos < 2; ++pos) {
                    Vector3 pt = (pos == 0) ? tv->point1 : tv->point2;

                    Probe newProbe;
                    newProbe.position           = pt;
                    newProbe.radius             = probeRadius;
                    newProbe.isPolar            = isPolarSum;
                    newProbe.numLayer           = 2;
                    newProbe.clusterId1         = 1;
                    newProbe.isSurvived         = true;
                    newProbe.contactAtomIndices = {j, k, i};
                    if (isMetal) newProbe.isPolar = 7;

                    if (bumpCheckAndScore(newProbe, state, 1.0, true)) {
                        newProbe.index = static_cast<int>(state.probes.size());
                        state.probes.push_back(std::move(newProbe));
                        insertProbeIntoGrid(state,
                            static_cast<int>(state.probes.size()) - 1);
                    }
                }
            }
        }
    }

    deduplicateProbeCellList(state);
    // Populate atom→probe neighbour lists needed by layer 3
    buildAtomProbeNeighborList(state);
}

// ─── Layer 3: atom-probe-probe probe generation ───────────────────────────────
// Mirrors generateSubSecondLayer().  Uses each atom's probe neighbour list
// (populated by buildAtomProbeNeighborList) to form (atom_i, probe_j, probe_k).

static void generateLayer3(SimulationState& state) {
    state.currentLayer = 3;

    int atomCount = static_cast<int>(state.atoms.size());
    for (int i = 0; i < atomCount; ++i) {
        // Copy the probe neighbour list; atom data itself won't move
        const std::vector<int> probeNbrs = state.atoms[i].neighborProbes;
        int  polarI     = static_cast<int>(state.atoms[i].polarity);
        int  probesSize = static_cast<int>(probeNbrs.size());

        for (int p = 0; p < probesSize; ++p) {
            int j = probeNbrs[p];
            // Copy probe positions before any potential vector reallocation
            const Vector3 posJ = state.probes[j].position;

            for (int q = p + 1; q < probesSize; ++q) {
                int k = probeNbrs[q];
                if (j == k) continue;

                const Vector3 posK = state.probes[k].position;

                double distJK = distance(posJ, posK) + 0.0001;
                double distIJ = distance(state.atoms[i].position, posJ) + 0.0001;
                double distIK = distance(state.atoms[i].position, posK) + 0.0001;

                int isPolarSum = polarI;
                int stepIdx    = (isPolarSum > 0) ? 8 : 7;

                double iR = 0.0, jR = 0.0, kR = 0.0, probeRadius = 0.0;
                int    check_ok = 0, n = 0;

                for (int step = 0; step < stepIdx && !check_ok; ++step) {
                    n           = step;
                    probeRadius = PROBE_RADIUS + RADIUS_INCREMENT * n;

                    iR = (polarI == 1) ? 1.8 + RADIUS_INCREMENT * n
                       : (polarI == 2) ? 1.5 + RADIUS_INCREMENT * n
                       :                 state.atoms[i].vdwRadius + 0.9;

                    jR = PROBE_RADIUS - RADIUS_INCREMENT * n;
                    kR = PROBE_RADIUS - RADIUS_INCREMENT * n;

                    if (jR + kR + 2*probeRadius > distJK &&
                        iR + jR + 2*probeRadius > distIJ &&
                        iR + kR + 2*probeRadius > distIK)
                        check_ok = 1;
                }

                if (!check_ok) continue;

                // Order: atom i, probe j, probe k  (matches original call)
                auto tv = getTripleVertex(state.atoms[i].position,
                                          posJ, posK,
                                          iR, jR, kR, probeRadius);
                if (!tv) continue;

                for (int pos = 0; pos < 2; ++pos) {
                    Vector3 pt = (pos == 0) ? tv->point1 : tv->point2;

                    Probe newProbe;
                    newProbe.position           = pt;
                    newProbe.radius             = probeRadius;
                    newProbe.isPolar            = isPolarSum;
                    newProbe.numLayer           = 3;
                    newProbe.clusterId          = 1;
                    newProbe.isSurvived         = true;
                    newProbe.contactAtomIndices = {i, j, k};
                    if (polarI == 2) newProbe.isPolar = 7;

                    if (bumpCheckAndScore(newProbe, state, 1.0, true)) {
                        newProbe.index = static_cast<int>(state.probes.size());
                        state.probes.push_back(std::move(newProbe));
                        insertProbeIntoGrid(state,
                            static_cast<int>(state.probes.size()) - 1);
                    }
                }
            }
        }
    }

    deduplicateProbeCellList(state);
}

// ─── Layer 4+: probe-probe-probe probe generation (iterating) ─────────────────
// Mirrors generateNextLayer().  Iterates until no new probes are added or the
// layer limit is reached.  Uses probe.neighborProbes (2.8 Å) built by
// assignProbeNeighbors before the first call.

static void generateLayer4(SimulationState& state) {
    state.currentLayer = 4;

    // Build initial probe-probe neighbour lists (2.8 Å, 117 grid properties)
    {
        std::vector<int> allIds;
        allIds.reserve(state.probes.size());
        for (int i = 0; i < static_cast<int>(state.probes.size()); ++i)
            allIds.push_back(i);
        assignProbeNeighbors(state, allIds, 2.8 * 2.8, 117);
    }

    int  totalLayer = 0;
    int  start      = 0;
    int  nobump     = 0;

    do {
        std::vector<int> newProbeIds;
        int end = static_cast<int>(state.probes.size());

        if (nobump != 0)
            start = end - nobump + 1;
        nobump = 0;

        for (int i = start; i < end; ++i) {
            // Copy fields before any potential reallocation
            const Vector3          posI   = state.probes[i].position;
            const std::vector<int> nbrsI  = state.probes[i].neighborProbes;

            for (int oi = 0; oi < static_cast<int>(nbrsI.size()); ++oi) {
                int j = nbrsI[oi];
                const Vector3 posJ = state.probes[j].position;

                for (int pi = oi + 1; pi < static_cast<int>(nbrsI.size()); ++pi) {
                    int k = nbrsI[pi];
                    if (i == j || i == k || j == k) continue;

                    const Vector3 posK = state.probes[k].position;

                    double distIJ = distance(posI, posJ) + 0.0001;
                    double distIK = distance(posI, posK) + 0.0001;
                    double distJK = distance(posJ, posK) + 0.0001;

                    if (distIJ > 2.8 + 0.0001 ||
                        distIK > 2.8 + 0.0001 ||
                        distJK > 2.8 + 0.0001) continue;

                    double iR = PROBE_RADIUS, jR = PROBE_RADIUS, kR = PROBE_RADIUS;
                    double probeRadius = PROBE_RADIUS;

                    if (iR + jR + 2*probeRadius <= distIJ ||
                        iR + kR + 2*probeRadius <= distIK ||
                        jR + kR + 2*probeRadius <= distJK) continue;

                    auto tv = getTripleVertex(posI, posJ, posK,
                                              iR, jR, kR, probeRadius);
                    if (!tv) continue;

                    for (int pos = 0; pos < 2; ++pos) {
                        Vector3 pt = (pos == 0) ? tv->point1 : tv->point2;

                        Probe newProbe;
                        newProbe.position           = pt;
                        newProbe.radius             = probeRadius;
                        newProbe.isPolar            = 0;
                        newProbe.numLayer           = totalLayer + 4;
                        newProbe.clusterId          = 0;
                        newProbe.isSurvived         = true;
                        newProbe.contactAtomIndices = {i, j, k};

                        // Layer 4+: use combined buriedness+surfaceCount filter,
                        // surface multiplier = 3 and probe-probe cutoff = 1.2 Å.
                        if (bumpCheckAndScore(newProbe, state,
                                              1.2, true, true, 3)) {
                            newProbe.index = static_cast<int>(state.probes.size());
                            state.probes.push_back(std::move(newProbe));
                            int newIdx = static_cast<int>(state.probes.size()) - 1;
                            insertProbeIntoGrid(state, newIdx);
                            newProbeIds.push_back(newIdx);
                            ++nobump;
                        }
                    }
                }
            }
        }

        deduplicateProbeCellList(state);
        // Assign 2.8 Å neighbours to newly added probes only
        if (!newProbeIds.empty())
            assignProbeNeighbors(state, newProbeIds, 2.8 * 2.8, 117);

        std::printf(
            "REMARK %3d layer (%3d - %3d), Survived: %3d, Total probes: %5zu\n",
            totalLayer + 4, start, end, nobump,
            state.probes.size());

        ++totalLayer;
    } while (nobump > 0 && totalLayer < MAX_LAYERS - 3);

    // Final pass: rebuild all probe-probe neighbour lists at 2.8 Å
    {
        std::vector<int> allIds;
        allIds.reserve(state.probes.size());
        for (int i = 0; i < static_cast<int>(state.probes.size()); ++i)
            allIds.push_back(i);
        assignProbeNeighbors(state, allIds, 2.8 * 2.8, 117);
    }
}

// ─── main ────────────────────────────────────────────────────────────────────

int main(int argc, char* argv[]) {
    if (argc < 2) { printUsage(argv[0]); return 1; }

    const std::string pdbFilename = argv[1];
    SimulationState state;

    // Optional BC cutoff ratio override
    if (argc > 2) state.bcCutoffRatioOverride = std::stod(argv[2]);

    clock_t startTotal = clock();

    // ── 1. Load configuration ─────────────────────────────────────────────
    {
        std::ifstream ifs("atom_property");
        if (!ifs) { std::cerr << "Cannot open atom_property\n"; return 1; }
        readAtomPropertyFile(ifs, state);
    }
    {
        std::ifstream ifs("grid_property");
        if (!ifs) { std::cerr << "Cannot open grid_property\n"; return 1; }
        readGridPropertyFile(ifs, state);
    }

    // ── 2. Read PDB ───────────────────────────────────────────────────────
    {
        std::ifstream ifs(pdbFilename);
        if (!ifs) { std::cerr << "Cannot open " << pdbFilename << "\n"; return 1; }
        readPdbFile(ifs, state);
    }

    std::printf("REMARK Read PDB [%6zu atoms]: "
                "Min(x,y,z): %8.3f %8.3f %8.3f  "
                "Max(x,y,z): %8.3f %8.3f %8.3f\n",
                state.atoms.size(),
                state.minX, state.minY, state.minZ,
                state.maxX, state.maxY, state.maxZ);

    // ── 3. Build cell list ────────────────────────────────────────────────
    clock_t t0 = clock();
    int pairCount = buildCellList(state);
    buildAtomNeighborList(state);    // populate atom.neighborAtoms (7.7 Å)
    computeAtomBuriedness(state);    // populate atom.buriedness (for surface counts)
    clock_t t1 = clock();
    std::printf("REMARK [%5.2fs] Near-atom pairs (<=7.7A): %d\n",
                elapsedSeconds(t0, t1), pairCount);

    // ── 4. Residue properties ─────────────────────────────────────────────
    initResidueProperties(state);

    // ── 5. Layer 1 probe generation ───────────────────────────────────────
    t0 = clock();
    generateLayer1(state);
    // Assign probe-probe neighbours at 1.5 Å for the initial DBSCAN
    {
        std::vector<int> allProbeIds;
        allProbeIds.reserve(state.probes.size());
        for (int i = 0; i < static_cast<int>(state.probes.size()); ++i)
            allProbeIds.push_back(i);
        assignProbeNeighbors(state, allProbeIds, 1.5 * 1.5, 275);
    }
    t1 = clock();
    std::printf("REMARK [%5.2fs] 1st Layer + probe-pair assignment: %zu probes\n",
                elapsedSeconds(t0, t1), state.probes.size());

    // ── 6. Initial DBSCAN + weed + BC filter ─────────────────────────────
    t0 = clock();
    std::size_t beforeWeed = state.probes.size();
    {
        std::vector<int> ids;
        ids.reserve(state.probes.size());
        for (int i = 0; i < static_cast<int>(state.probes.size()); ++i)
            ids.push_back(i);
        clusterProbes(state, ids, 0.7001, 1, 0);
    }
    weedOverlappingProbes(state, 0.7);
    rebuildProbeCellList(state);

    // Compute buriedness (numBc1) for surviving probes and derive threshold
    for (auto& p : state.probes) {
        auto bcRes = calculateBuriedness(p.position, state);
        p.buriedness = bcRes.buriedness;
        state.probeBcList.push_back(p.buriedness);
    }
    deriveBcCutoff(state);

    // Compute surface count (numBc2) and filter: keep if numBc1 + numBc2 >= bcCutoff
    {
        std::vector<Probe> survived;
        survived.reserve(state.probes.size());
        for (auto& p : state.probes) {
            int sc = surfaceCountLayer1(p.position, state);
            p.surfaceCount = static_cast<double>(sc);
            if (p.buriedness + sc >= state.bcCutoff)
                survived.push_back(std::move(p));
        }
        state.probes = std::move(survived);
    }
    rebuildProbeCellList(state);

    t1 = clock();
    std::printf("REMARK [%5.2fs] Weed/BC (< 0.7A) %zu -> %zu\n",
                elapsedSeconds(t0, t1), beforeWeed, state.probes.size());

    // ── 7. Layer 2: probe-atom-atom ───────────────────────────────────────
    t0 = clock();
    std::size_t before2 = state.probes.size();
    generateLayer2(state);
    t1 = clock();
    std::printf("REMARK [%5.2fs] 2nd Layer: %zu -> %zu probes\n",
                elapsedSeconds(t0, t1), before2, state.probes.size());

    // ── 8. Layer 3: atom-probe-probe ──────────────────────────────────────
    t0 = clock();
    std::size_t before3 = state.probes.size();
    generateLayer3(state);
    t1 = clock();
    std::printf("REMARK [%5.2fs] 3rd Layer: %zu -> %zu probes\n",
                elapsedSeconds(t0, t1), before3, state.probes.size());

    // ── 9. Layer 4+: probe-probe-probe ────────────────────────────────────
    t0 = clock();
    std::size_t before4 = state.probes.size();
    generateLayer4(state);
    t1 = clock();
    std::printf("REMARK [%5.2fs] 4th and further layers: %zu -> %zu probes\n",
                elapsedSeconds(t0, t1), before4, state.probes.size());

    // ── 10. Final DBSCAN clustering ───────────────────────────────────────
    {
        std::vector<int> allProbeIds;
        allProbeIds.reserve(state.probes.size());
        for (int i = 0; i < static_cast<int>(state.probes.size()); ++i)
            allProbeIds.push_back(i);
        assignProbeNeighbors(state, allProbeIds, 2.8 * 2.8, 117);
        clusterProbes(state, allProbeIds, 2.5, 16, 0);
        buildClusterSet(state);
    }

    // ── 11. PLB scores ────────────────────────────────────────────────────
    calculatePlbScores(state);

    clock_t endTotal = clock();
    std::printf("REMARK [%5.2fs] Total execution (%zu protein atoms)\n",
                elapsedSeconds(startTotal, endTotal), state.atoms.size());

    return 0;
}
