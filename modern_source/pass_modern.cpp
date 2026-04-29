// pass_modern.cpp – Modern C++17 entry point for mPASS
//
// Orchestrates the same pipeline as the original pass.cpp:
//   1. Read configuration and PDB files
//   2. Build cell list
//   3. Generate probe layers 1-4
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

// ─── First-layer probe generation (stub – mirrors generateFirstLayer) ─────────
// Full implementation would place probes at triple-sphere intersections for
// every triplet of near-atoms; shown here as a documented call site.

static void generateLayer1(SimulationState& state) {
    // For each pair (i,j) of near protein atoms, for every third atom k in the
    // neighbourhood of i, call probePositionsFromThreeAtoms and if the result
    // passes bumpCheckAndScore, add a new Probe to state.probes.
    // (Full implementation omitted for brevity; the geometry engine is in
    //  probe_placer.hpp and bc_calculator.hpp.)
    (void)state;
}

static void generateLaterLayers(SimulationState& state, int layer) {
    // Same logic as generateLayer1 but replaces atom-atom-atom triplets with
    // atom-atom-probe, atom-probe-probe, and probe-probe-probe combinations.
    (void)state; (void)layer;
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
    clock_t t1 = clock();
    std::printf("REMARK [%5.2fs] Near-atom pairs (<=7.7A): %d\n",
                elapsedSeconds(t0, t1), pairCount);

    // ── 4. Residue properties ─────────────────────────────────────────────
    initResidueProperties(state);

    // ── 5. Layer 1 probe generation ───────────────────────────────────────
    t0 = clock();
    generateLayer1(state);
    // Assign probe-probe neighbours at 1.5 Å for first DBSCAN
    std::vector<int> allProbeIds;
    allProbeIds.reserve(state.probes.size());
    for (int i = 0; i < static_cast<int>(state.probes.size()); ++i)
        allProbeIds.push_back(i);
    assignProbeNeighbors(state, allProbeIds, 1.5 * 1.5, 275);
    t1 = clock();
    std::printf("REMARK [%5.2fs] 1st Layer + probe-pair assignment: %zu probes\n",
                elapsedSeconds(t0, t1), state.probes.size());

    // ── 6. Initial DBSCAN + weed ──────────────────────────────────────────
    t0 = clock();
    std::size_t beforeWeed = state.probes.size();
    clusterProbes(state, allProbeIds, 0.7001, 1, 0);
    weedOverlappingProbes(state, 0.7);
    // Recalculate BC for surviving probes and derive threshold
    for (auto& p : state.probes) {
        auto bcRes = calculateBuriedness(p.position, state);
        p.buriedness = bcRes.buriedness;
        state.probeBcList.push_back(p.buriedness);
    }
    deriveBcCutoff(state);
    t1 = clock();
    std::printf("REMARK [%5.2fs] Weed/BC (< 0.7A) %zu -> %zu\n",
                elapsedSeconds(t0, t1), beforeWeed, state.probes.size());

    // ── 7. Layers 2-4 ────────────────────────────────────────────────────
    for (int layer = 2; layer <= 4; ++layer) {
        t0 = clock();
        std::size_t before = state.probes.size();
        generateLaterLayers(state, layer);
        t1 = clock();
        std::printf("REMARK [%5.2fs] Layer %d: %zu -> %zu probes\n",
                    elapsedSeconds(t0, t1), layer, before, state.probes.size());
    }

    // ── 8. Final DBSCAN clustering ────────────────────────────────────────
    allProbeIds.clear();
    for (int i = 0; i < static_cast<int>(state.probes.size()); ++i)
        allProbeIds.push_back(i);
    assignProbeNeighbors(state, allProbeIds, 2.8 * 2.8, 117);
    clusterProbes(state, allProbeIds, 2.5, 16, 0);
    buildClusterSet(state);

    // ── 9. PLB scores ─────────────────────────────────────────────────────
    calculatePlbScores(state);

    clock_t endTotal = clock();
    std::printf("REMARK [%5.2fs] Total execution (%zu protein atoms)\n",
                elapsedSeconds(startTotal, endTotal), state.atoms.size());

    return 0;
}
