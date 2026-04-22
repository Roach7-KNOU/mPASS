#include "pass.h"
#include "reader.h"

using namespace std;

/* Read protein atom from PDB file */
void readPdbFile(ifstream& ifs) {
    string buffer, tempBuffer, strKey, tempResidueName;
    int residueIdx = 0, tempResidueNumber = 9999, residueSize = 0, atomSerial = -1;
    double sumX = 0, sumY = 0, sumZ = 0, tempDist = 0, maxDist = 0;
    vector<Vector3> residueAtoms;
    vector<int> atomNumbers;
    Vector3 tempCoord;

    map<string, int>::iterator itChain;
    AtomPropertyMap::iterator itProp;

    auto finalizeResidue = [&](int realResidueNumber) {
        if (residueSize == 0) return;

        tempCoord.x = sumX / residueSize;
        tempCoord.y = sumY / residueSize;
        tempCoord.z = sumZ / residueSize;
        for (int i = 0; i < (int)residueAtoms.size(); i++) {
            tempDist = calculateDistance(residueAtoms[i], tempCoord);
            if (tempDist > maxDist) maxDist = tempDist;
        }

        Residue* residueProp = new Residue;
        residueProp->residueNumber = residueIdx;
        residueProp->realResidueNumber = realResidueNumber;
        residueProp->residueName = tempResidueName;
        residueProp->point = tempCoord;
        residueProp->maxRadius = maxDist;
        residueProp->atomNumbers = atomNumbers;
        residueProp->maxBc = 0;
        residueProp->isOk = 0;
        g_residues.push_back(residueProp);
    };

    while (getline(ifs, buffer)) {
        tempBuffer = buffer;
        if (tempBuffer.substr(0, 4).compare("ATOM") != 0 && tempBuffer.substr(0, 6).compare("HETATM") != 0) {
            continue;
        }

        const string atomName = trim(tempBuffer.substr(12, 4));
        const string residueName = tempBuffer.substr(17, 3);
        strKey = residueName + "-" + atomName;
        itProp = g_atomPropertyMap.find(strKey);
        if (itProp == g_atomPropertyMap.end()) {
            continue;
        }

        Atom* tempAtom = new Atom;
        tempAtom->serialNumber = atoi(tempBuffer.substr(6, 5).c_str());
        tempAtom->atomName = atomName;
        tempAtom->chain = tempBuffer.substr(21, 1);
        tempAtom->residueName = residueName;

        g_chainList.push_back(tempAtom->chain);
        itChain = g_chainSizeMap.find(tempAtom->chain);
        if (itChain == g_chainSizeMap.end()) {
            g_chainSizeMap[tempAtom->chain] = 1;
        } else {
            itChain->second++;
        }

        const int realResidueNumber = atoi(tempBuffer.substr(22, 4).c_str());
        tempAtom->residueIdx = realResidueNumber;
        tempAtom->point.x = atof(tempBuffer.substr(30, 8).c_str());
        tempAtom->point.y = atof(tempBuffer.substr(38, 8).c_str());
        tempAtom->point.z = atof(tempBuffer.substr(46, 8).c_str());

        tempAtom->vdwRadius = (itProp->second)->vdwRadius;
        tempAtom->isPolar = (itProp->second)->isPolar;
        tempAtom->standardResidueName = (itProp->second)->standardResidueName;
        tempAtom->proteinBc = 0;
        atomSerial++;

        g_minX = (g_minX > tempAtom->point.x - 5) ? tempAtom->point.x - 5 : g_minX;
        g_minY = (g_minY > tempAtom->point.y - 5) ? tempAtom->point.y - 5 : g_minY;
        g_minZ = (g_minZ > tempAtom->point.z - 5) ? tempAtom->point.z - 5 : g_minZ;

        g_maxX = (g_maxX < tempAtom->point.x + 5) ? tempAtom->point.x + 5 : g_maxX;
        g_maxY = (g_maxY < tempAtom->point.y + 5) ? tempAtom->point.y + 5 : g_maxY;
        g_maxZ = (g_maxZ < tempAtom->point.z + 5) ? tempAtom->point.z + 5 : g_maxZ;

        if (tempResidueNumber == 9999 || tempResidueNumber == realResidueNumber) {
            if (residueIdx == 0) {
                tempResidueNumber = realResidueNumber;
            }

            tempResidueName = tempAtom->standardResidueName;
            sumX += tempAtom->point.x;
            sumY += tempAtom->point.y;
            sumZ += tempAtom->point.z;

            tempCoord.x = tempAtom->point.x;
            tempCoord.y = tempAtom->point.y;
            tempCoord.z = tempAtom->point.z;

            residueAtoms.push_back(tempCoord);
            atomNumbers.push_back(atomSerial);
            residueSize++;
        } else {
            finalizeResidue(tempResidueNumber);
            residueIdx++;

            maxDist = 0;
            residueAtoms.clear();
            atomNumbers.clear();

            tempResidueNumber = realResidueNumber;
            tempResidueName = tempAtom->standardResidueName;
            sumX = tempAtom->point.x;
            sumY = tempAtom->point.y;
            sumZ = tempAtom->point.z;

            tempCoord.x = tempAtom->point.x;
            tempCoord.y = tempAtom->point.y;
            tempCoord.z = tempAtom->point.z;
            residueAtoms.push_back(tempCoord);
            atomNumbers.push_back(atomSerial);
            residueSize = 1;
        }

        tempAtom->residueNumber = residueIdx;
        g_proteinAtoms.push_back(tempAtom);
    }

    finalizeResidue(tempResidueNumber);

    residueAtoms.clear();
    atomNumbers.clear();
}

/* Read atom property */
void readAtomPropertyFile(ifstream& ifs) {
    string buffer, tempBuffer, strKey;
    AtomProperty* tempProp;
    while (getline(ifs, buffer)) {
        tempBuffer = buffer;
        tempProp = new AtomProperty;
        tempProp->residueName = tempBuffer.substr(0, 3);
        tempProp->atomName = trim(tempBuffer.substr(4, 4));
        tempProp->atomType = trim(tempBuffer.substr(9, 2));
        tempProp->charge = atof(tempBuffer.substr(11, 6).c_str());
        tempProp->vdwRadius = atof(tempBuffer.substr(20, 5).c_str());
        tempProp->isPolar = atoi(tempBuffer.substr(27, 1).c_str());
        tempProp->standardResidueName = trim(tempBuffer.substr(29, 3).c_str());

        strKey = tempProp->residueName + "-" + tempProp->atomName;
        if (!g_atomPropertyMap.insert(AtomPropertyMap::value_type(strKey, tempProp)).second) {
            delete tempProp;
        }
    }
}

/* Read grid property */
void readGridPropertyFile(ifstream& ifs) {
    string buffer, tempBuffer;
    GridProperty* tempProp;
    while (getline(ifs, buffer)) {
        tempBuffer = buffer;
        tempProp = new GridProperty;
        tempProp->gridIndex[0] = atoi(tempBuffer.substr(0, 3).c_str());
        tempProp->gridIndex[1] = atoi(tempBuffer.substr(3, 3).c_str());
        tempProp->gridIndex[2] = atoi(tempBuffer.substr(6, 3).c_str());
        tempProp->gridProperty[0] = atoi(tempBuffer.substr(9, 2).c_str());
        tempProp->gridProperty[1] = atoi(tempBuffer.substr(11, 2).c_str());
        tempProp->gridProperty[2] = atoi(tempBuffer.substr(13, 2).c_str());
        tempProp->max = atof(tempBuffer.substr(15, 7).c_str());
        tempProp->min = atof(tempBuffer.substr(22, 7).c_str());
        g_gridProperties.push_back(tempProp);
    }
}
