// tests_modern.cpp – Unit tests for the modern mPASS C++17 rewrite
//
// Build:
//   g++ -std=c++17 -O2 -I. tests_modern.cpp -o run_tests && ./run_tests

#include "geometry.hpp"
#include "data_types.hpp"
#include "grid.hpp"
#include "bc_calculator.hpp"
#include "probe_placer.hpp"
#include "dbscan.hpp"
#include "residue_properties.hpp"
#include "file_reader.hpp"

#include <sstream>
#include <cmath>
#include <iostream>
#include <functional>
#include <stdexcept>
#include <string>
#include <algorithm>

using namespace mpass;

// ─── Minimal test runner ─────────────────────────────────────────────────────

static int g_passed = 0, g_failed = 0;

static void runTest(const char* name, std::function<void()> fn) {
    std::cout << "[ RUN ] " << name << "\n";
    try {
        fn();
        std::cout << "[ OK  ] " << name << "\n";
        ++g_passed;
    } catch (const std::exception& e) {
        std::cout << "[FAIL ] " << name << " -- " << e.what() << "\n";
        ++g_failed;
    } catch (...) {
        std::cout << "[FAIL ] " << name << " -- unknown exception\n";
        ++g_failed;
    }
}

// Use regular functions for test cases to avoid comma-in-body macro issues
#define ASSERT_NEAR(a, b, tol) \
    do { if (std::abs((double)(a)-(double)(b)) > (tol)) \
        throw std::runtime_error(std::string("ASSERT_NEAR failed at line ") \
            + std::to_string(__LINE__) + ": " \
            + std::to_string((double)(a)) + " vs " + std::to_string((double)(b))); } while(0)

#define ASSERT_TRUE(expr) \
    do { if (!(expr)) throw std::runtime_error( \
        std::string("ASSERT_TRUE failed at line ") + std::to_string(__LINE__) \
        + ": " #expr); } while(0)

#define ASSERT_FALSE(expr) \
    do { if (expr) throw std::runtime_error( \
        std::string("ASSERT_FALSE failed at line ") + std::to_string(__LINE__) \
        + ": " #expr); } while(0)

#define ASSERT_EQ(a, b) \
    do { if ((a) != (b)) throw std::runtime_error( \
        std::string("ASSERT_EQ failed at line ") + std::to_string(__LINE__) \
        + ": " + std::to_string(a) + " != " + std::to_string(b)); } while(0)

// ─── Shared test data helpers ─────────────────────────────────────────────────

static std::string minimalAtomPropertyContent() {
    // Columns: 0-2=resName, 4-7=atomName, 9-10=atomType, 11-16=charge,
    //          20-24=vdwRadius, 27=isPolar, 29-31=stdResName
    // Verified layout: "ALA CA  CT  0.0000   1.800 0 ALA" (indices 0-31)
    return "ALA CA  CT  0.0000   1.800 0 ALA\n"
           "ALA CB  CT  0.0000   1.800 0 ALA\n"
           "GLY CA  CT  0.0000   1.800 0 GLY\n";
}

static std::string minimalPdbContent() {
    return
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00\n"
        "ATOM      2  CB  ALA A   1       2.000   3.000   4.000  1.00  0.00\n"
        "ATOM      3  CA  GLY A   2       5.000   6.000   7.000  1.00  0.00\n";
}

static SimulationState makeLinearProbes(int n, double spacing) {
    SimulationState state;
    state.minX = -10; state.maxX = n * spacing + 10;
    state.minY = -10; state.maxY = 10;
    state.minZ = -10; state.maxZ = 10;
    state.gridMaxX = 20; state.gridMaxY = 5; state.gridMaxZ = 5;
    for (int i = 0; i < n; ++i) {
        Probe p;
        p.index      = i;
        p.position   = {i * spacing, 0.0, 0.0};
        p.isSurvived = true;
        p.buriedness = 100;
        state.probes.push_back(p);
    }
    return state;
}

static void assignLinearNeighbors(SimulationState& state, double eps) {
    int n = static_cast<int>(state.probes.size());
    for (int i = 0; i < n; ++i) {
        state.probes[i].neighborProbes.clear();
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            if (distance(state.probes[i].position, state.probes[j].position) < eps)
                state.probes[i].neighborProbes.push_back(j);
        }
    }
}

// ═════════════════════════════════════════════════════════════════════════════
// Section 1 – geometry.hpp
// ═════════════════════════════════════════════════════════════════════════════

static void test_vector3_arithmetic() {
    Vector3 a{1, 2, 3}; Vector3 b{4, 5, 6};
    Vector3 s = a + b;
    ASSERT_NEAR(s.x, 5.0, 1e-9);
    ASSERT_NEAR(s.y, 7.0, 1e-9);
    ASSERT_NEAR(s.z, 9.0, 1e-9);
    ASSERT_NEAR((b - a).x, 3.0, 1e-9);
    ASSERT_NEAR(a.dot(b), 32.0, 1e-9);
    Vector3 cross = a.cross(b);
    ASSERT_NEAR(cross.x, -3.0, 1e-9);
    ASSERT_NEAR(cross.y,  6.0, 1e-9);
    ASSERT_NEAR(cross.z, -3.0, 1e-9);
}

static void test_distance() {
    Vector3 a{0,0,0}; Vector3 b{3,4,0};
    ASSERT_NEAR(distance(a, b), 5.0, 1e-9);
    ASSERT_NEAR(distanceSquared(a, b), 25.0, 1e-9);
}

static void test_normalize() {
    Vector3 v{0, 0, 5};
    Vector3 n = v.normalized();
    ASSERT_NEAR(n.z, 1.0, 1e-9);
    ASSERT_NEAR(n.length(), 1.0, 1e-9);
}

static void test_normalize_zero_throws() {
    Vector3 z{0, 0, 0};
    bool threw = false;
    try { z.normalized(); } catch (const std::domain_error&) { threw = true; }
    ASSERT_TRUE(threw);
}

static void test_sq_helper() {
    ASSERT_NEAR(sq(3.0), 9.0, 1e-9);
    ASSERT_NEAR(sq(-4.0), 16.0, 1e-9);
}

static void test_angle_orthogonal() {
    Vector3 origin{0,0,0}; Vector3 px{1,0,0}; Vector3 py{0,1,0};
    ASSERT_NEAR(angleDegrees(origin, px, py), 90.0, 1e-6);
}

// ═════════════════════════════════════════════════════════════════════════════
// Section 2 – file_reader.hpp
// ═════════════════════════════════════════════════════════════════════════════

static void test_read_atom_property_file() {
    SimulationState state;
    std::istringstream iss(minimalAtomPropertyContent());
    readAtomPropertyFile(iss, state);
    ASSERT_EQ(static_cast<int>(state.atomPropertyMap.size()), 3);
    ASSERT_TRUE(state.atomPropertyMap.count("ALA-CA") == 1);
    ASSERT_NEAR(state.atomPropertyMap.at("ALA-CA").vdwRadius, 1.80, 0.01);
}

static void test_read_pdb_atom_and_residue_count() {
    SimulationState state;
    std::istringstream apf(minimalAtomPropertyContent());
    readAtomPropertyFile(apf, state);
    std::istringstream pdb(minimalPdbContent());
    readPdbFile(pdb, state);
    ASSERT_EQ(static_cast<int>(state.atoms.size()), 3);
    ASSERT_EQ(static_cast<int>(state.residues.size()), 2);
}

static void test_read_pdb_bounding_box() {
    SimulationState state;
    std::istringstream apf(minimalAtomPropertyContent());
    readAtomPropertyFile(apf, state);
    std::istringstream pdb(minimalPdbContent());
    readPdbFile(pdb, state);
    ASSERT_NEAR(state.minX, 1.0 - 5.0, 1e-6);
    ASSERT_NEAR(state.maxX, 5.0 + 5.0, 1e-6);
}

static void test_read_pdb_unknown_atom_skipped() {
    SimulationState state;
    std::istringstream apf("ALA CA  CT  0.0000  0 1.80 0 ALA\n");
    readAtomPropertyFile(apf, state);
    std::istringstream pdb(minimalPdbContent());
    readPdbFile(pdb, state);
    // CB and GLY-CA are not in the map; only ALA-CA (res 1) loads
    ASSERT_TRUE(state.atoms.size() <= 2);
}

static void test_read_pdb_empty_stream() {
    SimulationState state;
    std::istringstream apf(minimalAtomPropertyContent());
    readAtomPropertyFile(apf, state);
    std::istringstream pdb("");
    readPdbFile(pdb, state);
    ASSERT_EQ(static_cast<int>(state.atoms.size()), 0);
    ASSERT_EQ(static_cast<int>(state.residues.size()), 0);
}

// ═════════════════════════════════════════════════════════════════════════════
// Section 3 – probe_placer.hpp
// ═════════════════════════════════════════════════════════════════════════════

static void test_triple_vertex_equilateral() {
    // Three identical atoms (r=1.7 Å) with pairwise touching distances (d=3.4 Å).
    // Circumradius of equilateral triangle = d/sqrt(3) ≈ 1.964 Å.
    // Probe (rp=0.7) must sit at distance r+rp=2.4 Å from each centre.
    // Since 1.964 < 2.4, a real solution (h > 0) exists above the plane.
    double r    = 1.7;
    double rp   = 0.7;
    double side = 2.0 * r;              // touching atoms: d = r_i + r_j = 3.4
    double h    = side * std::sqrt(3.0) / 2.0;
    Vector3 A{0.0,      0.0, 0.0};
    Vector3 B{side,     0.0, 0.0};
    Vector3 C{side/2.0, h,   0.0};
    auto result = getTripleVertex(A, B, C, r, r, r, rp);
    ASSERT_TRUE(result.has_value());
    ASSERT_NEAR(distance(A, result->point1), r + rp, 0.02);
    ASSERT_NEAR(distance(B, result->point1), r + rp, 0.02);
    ASSERT_NEAR(distance(C, result->point1), r + rp, 0.02);
}

static void test_triple_vertex_collinear() {
    Vector3 A{0,0,0}; Vector3 B{2,0,0}; Vector3 C{4,0,0};
    auto result = getTripleVertex(A, B, C, 1.7, 1.7, 1.7, 0.7);
    ASSERT_FALSE(result.has_value());
}

static void test_triple_vertex_coincident() {
    Vector3 A{0,0,0}; Vector3 B{0,0,0}; Vector3 C{4,0,0};
    auto result = getTripleVertex(A, B, C, 1.7, 1.7, 1.7, 0.7);
    ASSERT_FALSE(result.has_value());
}

// ═════════════════════════════════════════════════════════════════════════════
// Section 4 – dbscan.hpp
// ═════════════════════════════════════════════════════════════════════════════

static void test_dbscan_5_tight_probes_1_cluster() {
    auto state = makeLinearProbes(5, 0.5);
    assignLinearNeighbors(state, 1.5);
    std::vector<int> ids{0,1,2,3,4};
    int lastId = clusterProbes(state, ids, 1.5, 1, 0);
    ASSERT_EQ(lastId, 1);
    for (int i = 0; i < 5; ++i)
        ASSERT_EQ(state.probes[i].clusterId, 1);
}

static void test_dbscan_two_groups_2_clusters() {
    auto state = makeLinearProbes(6, 0.5);
    state.probes[3].position.x = 10.0;
    state.probes[4].position.x = 10.5;
    state.probes[5].position.x = 11.0;
    assignLinearNeighbors(state, 1.5);
    std::vector<int> ids{0,1,2,3,4,5};
    int lastId = clusterProbes(state, ids, 1.5, 1, 0);
    ASSERT_EQ(lastId, 2);
}

static void test_dbscan_single_probe_is_noise() {
    SimulationState state;
    Probe p; p.index = 0; p.position = {0,0,0}; p.isSurvived = true;
    state.probes.push_back(p);
    std::vector<int> ids{0};
    int lastId = clusterProbes(state, ids, 1.5, 2, 0);
    ASSERT_EQ(lastId, 0);
    ASSERT_EQ(state.probes[0].clusterId, 0);
}

static void test_dbscan_empty_set_returns_startid() {
    SimulationState state;
    std::vector<int> ids;
    int lastId = clusterProbes(state, ids, 1.5, 1, 42);
    ASSERT_EQ(lastId, 42);
}

// ═════════════════════════════════════════════════════════════════════════════
// Section 5 – bc_calculator.hpp
// ═════════════════════════════════════════════════════════════════════════════

static void test_bc_cutoff_clamped() {
    SimulationState state;
    for (int i = 0; i < 500; ++i) {
        ProteinAtom a; a.position = {static_cast<double>(i), 0, 0};
        state.atoms.push_back(a);
    }
    for (int i = 0; i < 10; ++i)  state.probeBcList.push_back(200);
    for (int i = 0; i < 90; ++i)  state.probeBcList.push_back(50);
    deriveBcCutoff(state);
    ASSERT_NEAR(state.top10Mean, 200.0, 1e-6);
    ASSERT_TRUE(state.bcCutoff >= 200.0 * 0.5 - 1e-6);
    ASSERT_TRUE(state.bcCutoff <= 200.0 * 0.8 + 1e-6);
}

static void test_bc_cutoff_override() {
    SimulationState state;
    state.bcCutoffRatioOverride = 0.65;
    for (int i = 0; i < 500; ++i) { ProteinAtom a; state.atoms.push_back(a); }
    for (int i = 0; i < 100; ++i)  state.probeBcList.push_back(100);
    deriveBcCutoff(state);
    ASSERT_NEAR(state.bcCutoff, 100.0 * 0.65, 1e-6);
}

// ═════════════════════════════════════════════════════════════════════════════
// Section 6 – residue_properties.hpp
// ═════════════════════════════════════════════════════════════════════════════

static void test_residue_props_20_entries() {
    SimulationState state;
    initResidueProperties(state);
    ASSERT_EQ(static_cast<int>(state.sogaIndexMap.size()), 20);
    ASSERT_NEAR(state.sogaIndexMap.at("TRP"), 2.518, 1e-4);
    ASSERT_NEAR(state.hydrophobicityMap.at("ILE"), 4.5, 1e-4);
}

static void test_plb_ala_ile() {
    SimulationState state;
    initResidueProperties(state);
    Residue rAla; rAla.residueName = "ALA"; state.residues.push_back(rAla);
    Residue rIle; rIle.residueName = "ILE"; state.residues.push_back(rIle);
    Cluster cl; cl.residuesContact = {0, 1};
    state.clusters.push_back(cl);
    calculatePlbScores(state);
    ASSERT_NEAR(state.clusters[0].plb,             0.701 + 1.006, 1e-4);
    ASSERT_NEAR(state.clusters[0].hydrophobicity,   1.8   + 4.5,  1e-4);
}

static void test_plb_unknown_residue_skipped() {
    SimulationState state;
    initResidueProperties(state);
    Residue rUnk; rUnk.residueName = "UNK"; state.residues.push_back(rUnk);
    Cluster cl; cl.residuesContact = {0};
    state.clusters.push_back(cl);
    calculatePlbScores(state);
    ASSERT_NEAR(state.clusters[0].plb, 0.0, 1e-9);
}

// ─── main ────────────────────────────────────────────────────────────────────

int main() {
    std::cout << "\n===== mPASS Modern C++ Tests =====\n\n";

    // geometry
    runTest("Vector3 arithmetic",                test_vector3_arithmetic);
    runTest("distance",                          test_distance);
    runTest("normalize",                         test_normalize);
    runTest("normalize zero vector throws",      test_normalize_zero_throws);
    runTest("sq helper",                         test_sq_helper);
    runTest("angle orthogonal = 90 deg",         test_angle_orthogonal);

    // file_reader
    runTest("readAtomPropertyFile basic",        test_read_atom_property_file);
    runTest("readPdbFile atom+residue count",    test_read_pdb_atom_and_residue_count);
    runTest("readPdbFile bounding box",          test_read_pdb_bounding_box);
    runTest("readPdbFile unknown atom skipped",  test_read_pdb_unknown_atom_skipped);
    runTest("readPdbFile empty stream",          test_read_pdb_empty_stream);

    // probe_placer
    runTest("getTripleVertex equilateral",       test_triple_vertex_equilateral);
    runTest("getTripleVertex collinear",         test_triple_vertex_collinear);
    runTest("getTripleVertex coincident",        test_triple_vertex_coincident);

    // dbscan
    runTest("DBSCAN 5 tight probes 1 cluster",  test_dbscan_5_tight_probes_1_cluster);
    runTest("DBSCAN two groups 2 clusters",      test_dbscan_two_groups_2_clusters);
    runTest("DBSCAN single probe is noise",      test_dbscan_single_probe_is_noise);
    runTest("DBSCAN empty set returns startId",  test_dbscan_empty_set_returns_startid);

    // bc_calculator
    runTest("BC cutoff clamped [0.5, 0.8]",     test_bc_cutoff_clamped);
    runTest("BC cutoff override",                test_bc_cutoff_override);

    // residue_properties
    runTest("20 residues populated",             test_residue_props_20_entries);
    runTest("PLB ALA+ILE accumulation",          test_plb_ala_ile);
    runTest("PLB unknown residue skipped",       test_plb_unknown_residue_skipped);

    std::cout << "\nResults: " << g_passed << " passed, " << g_failed << " failed.\n";
    return g_failed == 0 ? 0 : 1;
}
