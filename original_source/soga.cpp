#include "pass.h"
#include "cal.h"

void residuePropertySetInit() {
    g_sogaIndexMap["ALA"] = 0.701;
    g_sogaIndexMap["CYS"] = 1.650;   
    g_sogaIndexMap["ASP"] = 1.015;   
    g_sogaIndexMap["GLU"] = 0.956;   
    g_sogaIndexMap["PHE"] = 1.952;   
    g_sogaIndexMap["GLY"] = 0.788;   
    g_sogaIndexMap["HIS"] = 2.286;   
    g_sogaIndexMap["ILE"] = 1.006;   
    g_sogaIndexMap["LYS"] = 0.468;   
    g_sogaIndexMap["LEU"] = 1.045;   
    g_sogaIndexMap["MET"] = 1.894;   
    g_sogaIndexMap["ASN"] = 0.811;   
    g_sogaIndexMap["PRO"] = 0.212;   
    g_sogaIndexMap["GLN"] = 0.669;   
    g_sogaIndexMap["ARG"] = 0.751;   
    g_sogaIndexMap["SER"] = 0.880;   
    g_sogaIndexMap["THR"] = 0.701;   
    g_sogaIndexMap["VAL"] = 0.585;   
    g_sogaIndexMap["TRP"] = 2.518;   
    g_sogaIndexMap["TYR"] = 2.039;   

    g_hydrophobicityMap["ALA"] = 1.8;
    g_hydrophobicityMap["CYS"] = 2.5;   
    g_hydrophobicityMap["ASP"] = -3.5;   
    g_hydrophobicityMap["GLU"] = -3.5;   
    g_hydrophobicityMap["PHE"] = 2.8;   
    g_hydrophobicityMap["GLY"] = -0.4;   
    g_hydrophobicityMap["HIS"] = -3.2;   
    g_hydrophobicityMap["ILE"] = 4.5;   
    g_hydrophobicityMap["LYS"] = -3.9;   
    g_hydrophobicityMap["LEU"] = 3.8;   
    g_hydrophobicityMap["MET"] = 1.9;   
    g_hydrophobicityMap["ASN"] = -3.5;   
    g_hydrophobicityMap["PRO"] = -1.6;   
    g_hydrophobicityMap["GLN"] = -3.5;   
    g_hydrophobicityMap["ARG"] = -4.5;   
    g_hydrophobicityMap["SER"] = -0.8;   
    g_hydrophobicityMap["THR"] = -0.7;   
    g_hydrophobicityMap["VAL"] = 4.2;   
    g_hydrophobicityMap["TRP"] = -0.9;   
    g_hydrophobicityMap["TYR"] = -1.3;
}

double calculatePlb(double ratio) {
    int count = 0;
    double sumRa = 0, sumHyd = 0;
    ResiduePropertyMap::iterator itRa, itHyd;
    string residueName;

    for (int i = 0; i < (int)g_clusters.size(); i++) {
        sumRa = 0;
        sumHyd = 0;
        for (int j = 0; j < (int)g_clusters[i]->residuesContact.size(); j++) {
            residueName = g_residues[g_clusters[i]->residuesContact[j]]->residueName;
            itRa = g_sogaIndexMap.find(residueName);
            itHyd = g_hydrophobicityMap.find(residueName);

            if (itRa != g_sogaIndexMap.end()) {
                sumRa += itRa->second;
                sumHyd += itHyd->second;
            }
        }
        g_clusters[i]->plb = sumRa;
        g_clusters[i]->hydrophobicity = sumHyd;
    }
    return 0;
}
