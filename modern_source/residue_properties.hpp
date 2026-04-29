#pragma once
// residue_properties.hpp – SOGA index and hydrophobicity tables
//
// Replaces the hard-coded map initialisation in soga.cpp and the
// calculatePlb() function.

#include "data_types.hpp"
#include <unordered_map>
#include <string>
#include <stdexcept>

namespace mpass {

// ─── Static tables ───────────────────────────────────────────────────────────

/// Populate the SOGA (Soga et al.) index and Kyte-Doolittle hydrophobicity
/// tables in the SimulationState.  Replaces residuePropertySetInit().
inline void initResidueProperties(SimulationState& state) {
    // SOGA binding-site discrimination index (Soga et al., 2007)
    state.sogaIndexMap = {
        {"ALA", 0.701}, {"CYS", 1.650}, {"ASP", 1.015}, {"GLU", 0.956},
        {"PHE", 1.952}, {"GLY", 0.788}, {"HIS", 2.286}, {"ILE", 1.006},
        {"LYS", 0.468}, {"LEU", 1.045}, {"MET", 1.894}, {"ASN", 0.811},
        {"PRO", 0.212}, {"GLN", 0.669}, {"ARG", 0.751}, {"SER", 0.880},
        {"THR", 0.701}, {"VAL", 0.585}, {"TRP", 2.518}, {"TYR", 2.039},
    };

    // Kyte-Doolittle hydrophobicity scale
    state.hydrophobicityMap = {
        {"ALA",  1.8}, {"CYS",  2.5}, {"ASP", -3.5}, {"GLU", -3.5},
        {"PHE",  2.8}, {"GLY", -0.4}, {"HIS", -3.2}, {"ILE",  4.5},
        {"LYS", -3.9}, {"LEU",  3.8}, {"MET",  1.9}, {"ASN", -3.5},
        {"PRO", -1.6}, {"GLN", -3.5}, {"ARG", -4.5}, {"SER", -0.8},
        {"THR", -0.7}, {"VAL",  4.2}, {"TRP", -0.9}, {"TYR", -1.3},
    };
}

// ─── PLB (Protein Ligand-ability Binding) score ──────────────────────────────

/// Compute the SOGA-based PLB score and cumulative hydrophobicity for every
/// cluster.  Replaces calculatePlb().
///
/// For each cluster, the function sums the SOGA index and hydrophobicity of
/// all residues in contact (residuesContact list).  Unknown residue names are
/// silently skipped, matching the original `itRa != end()` guard.
inline void calculatePlbScores(SimulationState& state) {
    for (auto& cluster : state.clusters) {
        double sumSoga  = 0.0;
        double sumHydro = 0.0;

        for (int residueIdx : cluster.residuesContact) {
            const std::string& resName = state.residues[residueIdx].residueName;

            auto itSoga  = state.sogaIndexMap.find(resName);
            auto itHydro = state.hydrophobicityMap.find(resName);

            if (itSoga != state.sogaIndexMap.end()) {
                sumSoga  += itSoga->second;
                sumHydro += itHydro->second;  // same residues share both maps
            }
        }

        cluster.plb          = sumSoga;
        cluster.hydrophobicity = sumHydro;
    }
}

} // namespace mpass
