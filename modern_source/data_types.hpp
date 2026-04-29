#pragma once
// data_types.hpp – Clean, modern data structures for mPASS
// Replaces macro-heavy pass.h; uses value types and owning smart pointers.

#include "geometry.hpp"
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <memory>
#include <optional>

namespace mpass {

// ─── Physical constants / algorithm parameters ──────────────────────────────

inline constexpr double DIST_POLAR_APOLAR_MIN  = 3.3;
inline constexpr double DIST_POLAR_POLAR_MIN   = 2.5;
inline constexpr double DIST_POLAR_METAL_MIN   = 2.2;
inline constexpr double PROBE_RADIUS           = 0.7;
inline constexpr double PROBE_RADIUS_APOLAR    = 1.6;
inline constexpr double PROBE_RADIUS_POLAR     = DIST_POLAR_POLAR_MIN / 2.0;
inline constexpr double BC_DISTANCE_CUTOFF     = 14.0;  // sphere radius for BC counting
inline constexpr double MAX_LAYERS             = 30;

// Polarity codes stored in AtomProperty::isPolar
enum class Polarity : int { Apolar = 0, Polar = 1, Metal = 2 };

// Minimum probe–atom distance per polarity
[[nodiscard]] inline double minProbeAtomDist(Polarity polarity) noexcept {
    switch (polarity) {
        case Polarity::Polar:  return DIST_POLAR_POLAR_MIN;
        case Polarity::Metal:  return DIST_POLAR_METAL_MIN;
        default:               return DIST_POLAR_APOLAR_MIN;
    }
}

// ─── AtomProperty ───────────────────────────────────────────────────────────

struct AtomProperty {
    std::string residueName;
    std::string atomName;
    std::string atomType;
    std::string standardResidueName;
    double      charge{};
    double      vdwRadius{};
    Polarity    polarity{Polarity::Apolar};
};

// Key: "RESIDUE-ATOMNAME"
using AtomPropertyMap = std::unordered_map<std::string, AtomProperty>;

// ─── GridProperty ───────────────────────────────────────────────────────────

struct GridProperty {
    int    gridIndex[3]{};   // relative cell offsets (dx, dy, dz)
    int    gridFlags[3]{};
    double maxDist{};
    double minDist{};
};

// ─── Residue ────────────────────────────────────────────────────────────────

struct Residue {
    int         residueIndex{};       // sequential index within the protein
    int         realResidueNumber{};  // number as in PDB
    std::string residueName;
    Vector3     centroid;
    double      maxRadius{};
    int         maxBuriedness{};
    bool        isBinding{false};
    std::vector<int> atomIndices;
};

// ─── ProteinAtom ────────────────────────────────────────────────────────────

struct ProteinAtom {
    int         serialNumber{};
    std::string atomName;
    std::string residueName;
    std::string standardResidueName;
    std::string chain;
    int         residueNumber{};      // sequential (0-based index)
    int         realResidueNumber{};
    Vector3     position;
    double      vdwRadius{};
    double      charge{};
    Polarity    polarity{Polarity::Apolar};
    int         buriedness{};         // protein atom BC value
    std::vector<int> neighborAtoms;
    std::vector<int> neighborProbes;
};

// ─── Probe ──────────────────────────────────────────────────────────────────

struct Probe {
    int     index{};
    Vector3 position;
    std::string type;          // "C", "N", "O" for output
    double  radius{};
    bool    isSurvived{false};
    bool    isPolar{false};
    double  charge{};
    int     numLayer{};        // generation layer (1 = first, 2 = second, …)
    int     clusterId{};
    int     clusterIdFirst{};
    int     clusterId1{};
    int     isVisited{};

    // Buriedness metrics
    int    buriedness{};          // numBc1 (protein atoms within 14 Å)
    double surfaceCount{};        // numBc2 (surface-accessible contribution)
    int    buriednessFlag{};      // numBc3 (-1 = opening marker)

    // Geometry
    double closestAtomDist{9999.0};
    double closestProbeDist{9999.0};
    double averageAtomDist{};
    int    averageAtomCount{};
    int    density{};

    std::vector<int> contactAtomIndices;   // up to 3 atoms that generated this probe
    std::vector<int> neighborAtoms;
    std::vector<int> neighborProbes;
    std::vector<int> mergedProbes;         // probes within eps for DBSCAN
};

// ─── Grid Cell ──────────────────────────────────────────────────────────────

struct GridCell {
    std::vector<int> atomIndices;
    std::vector<int> probeIndices;
};

// ─── Cluster ────────────────────────────────────────────────────────────────

struct Cluster {
    int    clusterId{};
    int    maxBuriedness{};
    int    maxBuriedProbeId{};
    int    minBuriedness{9999};
    int    sumBuriedness{};
    int    probeCount{};
    int    maxLayer{};
    int    layer1{}, layer2{}, layer3{}, layer4{};
    int    opening{};
    double averageBuriedness{};
    double stdBuriedness{};
    double depth{};   // min dist from opening to deepest probe
    double width{};   // max pairwise probe distance
    double plb{};     // SOGA-based score
    double hydrophobicity{};
    std::vector<int> probeIds;
    std::vector<int> residues;         // side-chain contacts
    std::vector<int> residuesContact;  // all contacts

    // BC bins (10 equal decile bands between bcCutoff and top10Mean)
    std::vector<int> bcBins[10];
    std::vector<int> bcBinTop;   // probes with buriednessFlag == -1
};

// ─── Global-state bundle (replaces file-scope globals) ──────────────────────

struct SimulationState {
    AtomPropertyMap               atomPropertyMap;
    std::vector<GridProperty>     gridProperties;
    std::unordered_map<std::string, double> sogaIndexMap;
    std::unordered_map<std::string, double> hydrophobicityMap;

    std::vector<ProteinAtom>      atoms;
    std::vector<Residue>          residues;
    std::vector<Probe>            probes;
    std::vector<GridCell>         cells;
    std::vector<Cluster>          clusters;

    // Bounding box (expanded by ±5 Å)
    double minX{9999}, minY{9999}, minZ{9999};
    double maxX{-9999}, maxY{-9999}, maxZ{-9999};

    // Grid dimensions
    int gridMaxX{}, gridMaxY{}, gridMaxZ{};
    int currentLayer{};

    // BC thresholds derived from probe statistics
    double top10Mean{};
    double quartile{};
    double bcCutoff{-1.0};
    double bcCutoffRatio{0.7};
    double bcCutoffRatioOverride{-1.0};

    std::vector<int> atomCellList;
    std::vector<int> probeCellList;
    std::vector<int> probeBcList;
    std::vector<std::string> chainList;
    std::unordered_map<std::string, int> chainSizeMap;
};

} // namespace mpass
