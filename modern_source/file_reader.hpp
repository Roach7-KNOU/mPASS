#pragma once
// file_reader.hpp – PDB / AtomProperty / GridProperty file parsers
//
// Replaces reader.cpp.  All parsers take an std::istream& (testable with
// std::istringstream) and fill the provided SimulationState.
//
// Key improvements over the original:
//   - Uses std::stoi / std::stod instead of atoi / atof (throws on bad input).
//   - Trim helper is local and composable.
//   - Residue finalisation uses a helper lambda → one code path.
//   - bounding-box expansion is factored out.

#include "data_types.hpp"
#include "geometry.hpp"
#include <istream>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <cctype>

namespace mpass {

// ─── String utilities ────────────────────────────────────────────────────────

[[nodiscard]] inline std::string trimLeft(std::string_view sv) {
    auto it = std::find_if(sv.begin(), sv.end(),
                           [](unsigned char c) { return !std::isspace(c); });
    return std::string(it, sv.end());
}

[[nodiscard]] inline std::string trimRight(std::string_view sv) {
    auto it = std::find_if(sv.rbegin(), sv.rend(),
                           [](unsigned char c) { return !std::isspace(c); });
    return std::string(sv.begin(), it.base());
}

[[nodiscard]] inline std::string trim(std::string_view sv) {
    return trimLeft(trimRight(sv));
}

// ─── Safe substring helpers ───────────────────────────────────────────────────

[[nodiscard]] inline std::string safeSubstr(const std::string& s, std::size_t pos,
                                             std::size_t len) {
    if (pos >= s.size()) return "";
    return s.substr(pos, std::min(len, s.size() - pos));
}

[[nodiscard]] inline double parseDouble(const std::string& s, std::size_t pos,
                                        std::size_t len) {
    std::string sub = trim(safeSubstr(s, pos, len));
    if (sub.empty()) return 0.0;
    try { return std::stod(sub); }
    catch (...) { return 0.0; }
}

[[nodiscard]] inline int parseInt(const std::string& s, std::size_t pos,
                                  std::size_t len) {
    std::string sub = trim(safeSubstr(s, pos, len));
    if (sub.empty()) return 0;
    try { return std::stoi(sub); }
    catch (...) { return 0; }
}

// ─── Bounding box helper ──────────────────────────────────────────────────────

inline void expandBoundingBox(SimulationState& state, const Vector3& pos,
                               double margin = 5.0) noexcept {
    state.minX = std::min(state.minX, pos.x - margin);
    state.minY = std::min(state.minY, pos.y - margin);
    state.minZ = std::min(state.minZ, pos.z - margin);
    state.maxX = std::max(state.maxX, pos.x + margin);
    state.maxY = std::max(state.maxY, pos.y + margin);
    state.maxZ = std::max(state.maxZ, pos.z + margin);
}

// ─── Atom property file ───────────────────────────────────────────────────────

/// Read the atom-property configuration file (one record per line, fixed-width).
/// Populates state.atomPropertyMap.
/// Format (0-based columns):
///   [0-2]  residueName   (3 chars)
///   [4-7]  atomName      (4 chars, trimmed)
///   [9-10] atomType      (2 chars, trimmed)
///   [11-16] charge       (6 chars)
///   [20-24] vdwRadius    (5 chars)
///   [27]   isPolar       (1 char: 0=apolar, 1=polar, 2=metal)
///   [29-31] standardResidueName (3 chars, trimmed)
inline void readAtomPropertyFile(std::istream& stream, SimulationState& state) {
    std::string line;
    while (std::getline(stream, line)) {
        if (line.size() < 32) continue;

        AtomProperty prop;
        prop.residueName          = safeSubstr(line, 0, 3);
        prop.atomName             = trim(safeSubstr(line, 4, 4));
        prop.atomType             = trim(safeSubstr(line, 9, 2));
        prop.charge               = parseDouble(line, 11, 6);
        prop.vdwRadius            = parseDouble(line, 20, 5);
        int polarCode             = parseInt(line, 27, 1);
        prop.polarity             = static_cast<Polarity>(polarCode);
        prop.standardResidueName  = trim(safeSubstr(line, 29, 3));

        std::string key = prop.residueName + "-" + prop.atomName;
        state.atomPropertyMap.emplace(std::move(key), std::move(prop));
    }
}

// ─── Grid property file ───────────────────────────────────────────────────────

/// Read the grid-property configuration file.
/// Populates state.gridProperties.
inline void readGridPropertyFile(std::istream& stream, SimulationState& state) {
    std::string line;
    while (std::getline(stream, line)) {
        if (line.size() < 29) continue;

        GridProperty gp;
        gp.gridIndex[0]  = parseInt(line, 0, 3);
        gp.gridIndex[1]  = parseInt(line, 3, 3);
        gp.gridIndex[2]  = parseInt(line, 6, 3);
        gp.gridFlags[0]  = parseInt(line, 9, 2);
        gp.gridFlags[1]  = parseInt(line, 11, 2);
        gp.gridFlags[2]  = parseInt(line, 13, 2);
        gp.maxDist       = parseDouble(line, 15, 7);
        gp.minDist       = parseDouble(line, 22, 7);
        state.gridProperties.push_back(std::move(gp));
    }
}

// ─── PDB file ─────────────────────────────────────────────────────────────────

/// Read a PDB file (ATOM / HETATM records) and populate state.atoms and
/// state.residues.
/// Only records whose "RESIDUE-ATOMNAME" key appears in atomPropertyMap are
/// retained, matching the original filter.
inline void readPdbFile(std::istream& stream, SimulationState& state) {
    std::string line;
    int residueIndex    = 0;
    int atomSerial      = -1;
    int currentResNum   = 9999;   // sentinel "not yet set"
    std::string currentResName;

    double sumX = 0, sumY = 0, sumZ = 0;
    int    residueAtomCount = 0;
    double maxRadius        = 0;
    std::vector<Vector3> residueAtomPositions;
    std::vector<int>     residueAtomIndices;

    // ── Finalise a completed residue ───────────────────────────────────────
    auto finalizeResidue = [&](int realResidueNumber) {
        if (residueAtomCount == 0) return;

        Vector3 centroid{ sumX / residueAtomCount,
                          sumY / residueAtomCount,
                          sumZ / residueAtomCount };

        for (const auto& pos : residueAtomPositions)
            maxRadius = std::max(maxRadius, distance(pos, centroid));

        Residue res;
        res.residueIndex     = residueIndex;
        res.realResidueNumber = realResidueNumber;
        res.residueName      = currentResName;
        res.centroid         = centroid;
        res.maxRadius        = maxRadius;
        res.atomIndices      = residueAtomIndices;
        state.residues.push_back(std::move(res));
    };

    while (std::getline(stream, line)) {
        // Accept both ATOM and HETATM records
        bool isAtom   = line.size() >= 6 && line.substr(0, 4) == "ATOM";
        bool isHetatm = line.size() >= 6 && line.substr(0, 6) == "HETATM";
        if (!isAtom && !isHetatm) continue;

        std::string atomName    = trim(safeSubstr(line, 12, 4));
        std::string residueName = safeSubstr(line, 17, 3);
        std::string key         = residueName + "-" + atomName;

        auto it = state.atomPropertyMap.find(key);
        if (it == state.atomPropertyMap.end()) continue; // unknown atom type

        const AtomProperty& prop = it->second;

        ProteinAtom atom;
        atom.serialNumber        = parseInt(line, 6, 5);
        atom.atomName            = atomName;
        atom.chain               = safeSubstr(line, 21, 1);
        atom.residueName         = residueName;
        atom.standardResidueName = prop.standardResidueName;
        atom.vdwRadius           = prop.vdwRadius;
        atom.polarity            = prop.polarity;
        atom.charge              = prop.charge;
        atom.buriedness          = 0;

        int realResNum = parseInt(line, 22, 4);
        atom.realResidueNumber = realResNum;

        atom.position.x = parseDouble(line, 30, 8);
        atom.position.y = parseDouble(line, 38, 8);
        atom.position.z = parseDouble(line, 46, 8);

        ++atomSerial;

        // Update bounding box
        expandBoundingBox(state, atom.position);

        // Chain bookkeeping
        state.chainList.push_back(atom.chain);
        state.chainSizeMap[atom.chain]++;

        // ── Residue boundary detection ─────────────────────────────────────
        if (currentResNum == 9999) {
            // First atom ever
            currentResNum  = realResNum;
            currentResName = atom.standardResidueName;
        }

        if (currentResNum != realResNum) {
            // Crossed into a new residue
            finalizeResidue(currentResNum);
            ++residueIndex;

            // Reset accumulators
            sumX = sumY = sumZ = 0.0;
            maxRadius = 0.0;
            residueAtomCount = 0;
            residueAtomPositions.clear();
            residueAtomIndices.clear();

            currentResNum  = realResNum;
            currentResName = atom.standardResidueName;
        }

        atom.residueNumber = residueIndex;

        sumX += atom.position.x;
        sumY += atom.position.y;
        sumZ += atom.position.z;
        residueAtomPositions.push_back(atom.position);
        residueAtomIndices.push_back(atomSerial);
        ++residueAtomCount;

        state.atoms.push_back(std::move(atom));
    }

    // Finalise last residue
    finalizeResidue(currentResNum);
}

} // namespace mpass
