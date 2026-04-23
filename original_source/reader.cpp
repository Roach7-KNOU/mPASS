#include "pass.h"
#include "reader.h"
#include <memory>

using namespace std;

/* Read protein atom from PDB file */
void readPdbFile(ifstream& ifs) {
    const size_t kInitialResidueCapacity = 64;
    string buffer, strKey, tempResidueName;
    int residueIdx = 0, tempResidueNumber = 9999, realResidueNumber, residueSize = 0, atomSerial = -1;
    double sumX = 0, sumY = 0, sumZ = 0, tempDist = 0, maxDist = 0;
    vector<Vector3> residueAtoms;
    vector<int> atomNumbers;
    Vector3 tempCoord;
    Residue* residueProp;

    map<string, int>::iterator itChain;
    AtomPropertyMap::iterator itProp;
    residueAtoms.reserve(kInitialResidueCapacity);
    atomNumbers.reserve(kInitialResidueCapacity);

    while (getline(ifs, buffer)) {
        if (buffer.compare(0, 4, "ATOM") == 0 || buffer.compare(0, 6, "HETATM") == 0) {
            const string atomName = trim(buffer.substr(12, 4));
            const string chain = buffer.substr(21, 1);
            const string residueName = buffer.substr(17, 3);

            g_chainList.push_back(chain);
            itChain = g_chainSizeMap.find(chain);
            if (itChain == g_chainSizeMap.end()) {
                g_chainSizeMap[chain] = 1;
            } else {
                itChain->second++;
            }

            strKey = residueName + "-" + atomName;

            itProp = g_atomPropertyMap.find(strKey);

            if (itProp != g_atomPropertyMap.end()) {
                std::unique_ptr<Atom> tempAtom = std::make_unique<Atom>();
                tempAtom->serialNumber = atoi(buffer.substr(6, 5).c_str());
                tempAtom->atomName = atomName;
                tempAtom->chain = chain;
                tempAtom->residueName = residueName;
                tempAtom->residueIdx = atoi(buffer.substr(22, 4).c_str());
                realResidueNumber = tempAtom->residueIdx;
                tempAtom->point.x = atof(buffer.substr(30, 8).c_str());
                tempAtom->point.y = atof(buffer.substr(38, 8).c_str());
                tempAtom->point.z = atof(buffer.substr(46, 8).c_str());
            
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
                    if (residueIdx == 0) 
                        tempResidueNumber = realResidueNumber;

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
                    tempCoord.x = sumX / residueSize;
                    tempCoord.y = sumY / residueSize;
                    tempCoord.z = sumZ / residueSize;
                    for (int i = 0; i < (int)residueAtoms.size(); i++) {
                        tempDist = calculateDistance(residueAtoms[i], tempCoord);
                        if (tempDist > maxDist) maxDist = tempDist;
                    }

                    residueProp = new Residue;
                    residueProp->residueNumber = residueIdx;
                    residueProp->realResidueNumber = realResidueNumber;
                    residueProp->residueName = tempResidueName;
                    residueProp->point = tempCoord;
                    residueProp->maxRadius = maxDist;
                    residueProp->atomNumbers = atomNumbers;
                    residueProp->maxBc = 0;
                    residueProp->isOk = 0;
                    g_residues.push_back(residueProp);

                    residueIdx++;
                    maxDist = 0;
                    residueAtoms.clear();
                    atomNumbers.clear();

                    tempResidueNumber = realResidueNumber;
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
                g_proteinAtoms.push_back(tempAtom.release());

            } else {
                //cout << "REMARK  ERROR : " << strKey << endl;
            }
        }
    }

    tempCoord.x = sumX / residueSize;
    tempCoord.y = sumY / residueSize;
    tempCoord.z = sumZ / residueSize;
    for (int i = 0; i < (int)residueAtoms.size(); i++) {
        tempDist = calculateDistance(residueAtoms[i], tempCoord);
        if (tempDist > maxDist) maxDist = tempDist;
    }
    residueProp = new Residue;
    residueProp->residueNumber = residueIdx;
    residueProp->point = tempCoord;
    residueProp->residueName = tempResidueName;
    residueProp->maxRadius = maxDist;
    residueProp->atomNumbers = atomNumbers;
    residueProp->realResidueNumber = realResidueNumber;
    residueProp->maxBc = 0;
    residueProp->isOk = 0;

    g_residues.push_back(residueProp);

    residueAtoms.clear();
    atomNumbers.clear();
}

/* Read atom property */
void readAtomPropertyFile(ifstream& ifs) {
    string buffer, strKey;
    AtomProperty* tempProp;
    while (getline(ifs, buffer)) {
        tempProp = new AtomProperty;
        tempProp->residueName = buffer.substr(0, 3);
        tempProp->atomName = trim(buffer.substr(4, 4));
        tempProp->atomType = trim(buffer.substr(9, 2));
        tempProp->charge = atof(buffer.substr(11, 6).c_str());
        tempProp->vdwRadius = atof(buffer.substr(20, 5).c_str());
        tempProp->isPolar = atoi(buffer.substr(27, 1).c_str());
        tempProp->standardResidueName = trim(buffer.substr(29, 3).c_str());

        strKey = tempProp->residueName + "-" + tempProp->atomName;
        if (!g_atomPropertyMap.insert(AtomPropertyMap::value_type(strKey, tempProp)).second) {
            delete tempProp;
        }
    }
}

/* Read grid property */
void readGridPropertyFile(ifstream& ifs) {
    string buffer;
    GridProperty* tempProp;
    while (getline(ifs, buffer)) {
        tempProp = new GridProperty;
        tempProp->gridIndex[0] = atoi(buffer.substr(0, 3).c_str());
        tempProp->gridIndex[1] = atoi(buffer.substr(3, 3).c_str());
        tempProp->gridIndex[2] = atoi(buffer.substr(6, 3).c_str());
        tempProp->gridProperty[0] = atoi(buffer.substr(9, 2).c_str());
        tempProp->gridProperty[1] = atoi(buffer.substr(11, 2).c_str());
        tempProp->gridProperty[2] = atoi(buffer.substr(13, 2).c_str());
        tempProp->max = atof(buffer.substr(15, 7).c_str());
        tempProp->min = atof(buffer.substr(22, 7).c_str());
        g_gridProperties.push_back(tempProp);
    }
}
