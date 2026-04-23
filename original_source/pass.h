#ifndef PASS_H
#define PASS_H

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <map>
#include <cmath>
#include <algorithm>

#define TRUE		1
#define FALSE		0
#define SQUARE(x)	(x)*(x)
#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)<(b)?(b):(a))

#define prAtom(i)	g_proteinAtoms[i]->point
#define prRadius(i)	g_proteinAtoms[i]->vdwRadius
#define prX(i)		g_proteinAtoms[i]->point.x
#define prY(i)		g_proteinAtoms[i]->point.y
#define prZ(i)		g_proteinAtoms[i]->point.z
#define prType(i)   g_proteinAtoms[i]->atomName
#define prPolar(i)	g_proteinAtoms[i]->isPolar
#define prChain(i)	g_proteinAtoms[i]->chain

#define ResiNum(i)	g_proteinAtoms[i]->residueNumber
#define ResiIdx(i)	g_proteinAtoms[i]->residueNumber
#define ResiName(i)	g_proteinAtoms[i]->standardResidueName
#define prCharge(i)	g_proteinAtoms[i]->charge
#define nrAtom(i,j)	g_proteinAtoms[i]->nearAtoms[j]
#define prisOK(i)	g_proteinAtoms[i]->isOk

#define residueAtoms(i)	    g_residues[i]->atomNumbers
#define residueAtom(i,j)    g_residues[i]->atomNumbers[j]
#define residueName(i)	    g_residues[i]->residueName
#define residueNumber(i)	g_residues[i]->residueNumber
#define residueType(i)		g_residues[i]->isOk
#define residueProbeGen(i)  g_residues[selectedResidue]->probeGenResidue

#define dist_check(i,j,r) vsDist(prAtom(i), prAtom(j)) <= (prRadius(i)+prRadius(j)+2*r)
#define dist_check_probe(i,j,r) vsDist(prAtom(i), pAtom(j)) <= (prRadius(i)+pRadius(j)+2*r)
#define dist_check_probe_2(i,j,r) vsDist(pAtom(i), pAtom(j)) <= (4*r)

#define pAtom(i)	    g_probes[i]->point
#define pASP(i)		    g_probes[i]->checkAsp
#define pX(i)		    g_probes[i]->point.x
#define pY(i)		    g_probes[i]->point.y
#define pZ(i)		    g_probes[i]->point.z
#define pCAtoms(i,j)	g_probes[i]->contactAtoms[j]
#define pType(i)	    g_probes[i]->type
#define pVDW_r(i)	    g_probes[i]->radius
#define pRadius(i)	    g_probes[i]->radius
#define pCharge(i)	    g_probes[i]->charge
#define pPolar(i)	    g_probes[i]->isPolar
#define pDist(i)	    g_probes[i]->closestDist
#define pPDist(i)	    g_probes[i]->closestProbeDist
#define pADist(i)	    g_probes[i]->averageDist
#define pANum(i)	    g_probes[i]->averageNum
#define pIsSurvived(i)	g_probes[i]->isSurvived
#define pClustID(i)	    g_probes[i]->clusterId
#define pNLayer(i)	    g_probes[i]->numLayer
#define LOWKEY(x)	    (x & 0xffff)
#define HIGHKEY(x)	    (((x) >> 16) & 0xffff)

#define DISTANCE_POLAR_APOLAR_MIN	3.3  
#define DISTANCE_POLAR_POLAR_MIN	2.5  
#define DISTANCE_POLAR_METAL_MIN	2.2  
#define RADIUS_INCREMENT		0.05 
#define RADIUS_INCREMENT_STEP	10

#define RADIUS_APOLAR_SP2		2.05
#define RADIUS_APOLAR_SP3		2.25

#define PROBE_RADIUS_PROBE		0.7	
#define PROBE_RADIUS_APOLAR		1.6	
#define PROBE_RADIUS_POLAR		DISTANCE_POLAR_POLAR_MIN/2
#define PROBE_RADIUS_METAL		DISTANCE_POLAR_POLAR_MIN/2

#define APOLAR_RANGE1			DISTANCE_POLAR_APOLAR_MIN - PROBE_RADIUS_POLAR 
#define APOLAR_RANGE2			DISTANCE_POLAR_APOLAR_MIN - PROBE_RADIUS_POLAR-0.2

#define BC_MIN_SCALE			0.3 
#define	BC_MIN				160 

#define	MAX_LAYER			30
#define MAX_GRID 			125000

#define PI					3.141

using namespace std;

struct Vector3 {
    double x, y, z;
}; 

struct AtomProperty {
    string residueName;    
    string atomName;	    
    string atomType;	    
    double charge;	    
    double vdwRadius;	    
    int isPolar;		    
    string standardResidueName; 
};

struct ResidueProperty {
    double ca, sa, ra, hydrophobicity; 
};

struct GridProperty {
	int gridIndex[3];
	int gridProperty[3];
	double max;
	double min;
};

struct Atom {
    int serialNumber;	    
    string atomName;	    
    string residueName;    
    int residueNumber;	    
    Vector3 point;	    
    double vdwRadius;	    
    double charge;	    
    int isPolar;		    
    string atomType;
    string standardResidueName;  
	string chain;
    int residueIdx;	    
    vector<int> nearAtoms;   	
	vector<int> nearProbes;
	int proteinBc;
};

struct Residue {
    int residueNumber;		    
    int realResidueNumber;
	string residueName;	    
    vector<int> atomNumbers;	       
    Vector3 point;		    
    double maxRadius;		    
    int numBc1, numBc2, numBc3; 
	int maxBc;
    int isOk;
};

struct Probe {
    int serialNumber;		    
    Vector3 point;		    
    string type;		    
    double radius;		    
    int position;		    
    int isPolar;			    
    double charge;		    
    int clusterIdFirst;	    
    int clusterId;		    
	int clusterId1;
	int clusterId2;
	int clusterId3;
    int isSurvived;		    
    vector<int> contactAtoms;  
    double closestDist;		    
    double closestProbeDist;	    
    double averageDist;		    
    int averageNum;		    
    int closestAtom;		    
    vector<int> mergedAtoms;   
    int checkAsp;
    int numBc1, numBc3, numBc4, numBc5;
	double numBc2;
    double numBcProbe, numBcProbeTotal, bcPercent;
	double numWtBcNum;
    vector<int> nearProbes;
	vector<int> nearProbesTemp;	    
    int numLayer;		    
    int isVisited;
    vector<int> nearAtoms;	    
    int density;		    
	int density1;
	int density2;
	int step;
};

struct Cell {
	vector<int> atoms;
	vector<int> probes;
	vector<int> bumpAtomList;
	vector<int> bcAtomList;
};

struct GridPoint {
	Vector3 point;
	vector<int> probes;
	vector<double> dist;
	int flag;
};

struct Cluster {
    int clusterId;
    int probeSize;
    vector<int> volCount1, volCount2, volCount3, volCount4;
	int binStart;
	int binNumber;
	vector<int> bcBin1, bcBin2, bcBin3, bcBin4, bcBin5;
	vector<int> bcBin6, bcBin7, bcBin8, bcBin9, bcBin10;
	vector<int> bcBinTop, bcBinBottom;
    double plb; 
    double zPlb;
    double hydrophobicity;
	double depth;
	double width;
    vector<int> probeIds;
    vector<int> residues;
	vector<int> residuesContact;
	int maxBcNum;
	int maxBcId;
    int minBcNum;
    int sumBcNum;
	int sumBcNum4;
	int sumProbeNum;
	int maxLayer;
	int minLayer;
	int layer1, layer2, layer3, layer4;
	int layer1Polar, layer2Polar, layer3Polar;
	int opening;
    double averageBcNum;
	double averageBcNum4;
	double averageBcNumS;
    double averageWtBcNum;
	double sumWtBcNum;
	double sumSd;
	double zBcmM, zBcAv, zBcS, zBcWt, zSd, zVol, zwtSd;
    double stdBcNum, stdBcNum4, stdBcNumS, stdWtBcNum;
    double zBcNum;
	double totalScore;
	double it;
    int type;
    int flag;
	int isOccluded;
};

struct ProbePropertySnapshot {
	int index;
	int serialNumber;
	string type;
	double x;
	double y;
	double z;
	double radius;
	double charge;
	int isPolar;
	int isSurvived;
	int clusterId;
	int clusterId1;
	int clusterId2;
	int clusterId3;
	int numLayer;
	int density;
	int density1;
	int density2;
	int step;
	int numBc1;
	double numBc2;
	int numBc3;
	int numBc4;
	int numBc5;
	double numBcProbe;
	double numBcProbeTotal;
	double bcPercent;
	double numWtBcNum;
	double closestDist;
	double closestProbeDist;
	double averageDist;
	int averageNum;
	int closestAtom;
	vector<int> contactAtoms;
	vector<int> nearProbes;
	vector<int> nearAtoms;
};

typedef map<string, AtomProperty*> AtomPropertyMap;
extern AtomPropertyMap g_atomPropertyMap;

typedef map<string, double> ResiduePropertyMap;
extern ResiduePropertyMap g_sogaIndexMap;
extern ResiduePropertyMap g_hydrophobicityMap;

typedef map<string, string> ResidueNameMap;
extern ResidueNameMap g_residueNames;

typedef vector<GridProperty*> GridProperties;
extern GridProperties g_gridProperties;

typedef vector<vector<int>> GridCells;
extern GridCells g_gridCells;

typedef vector<Atom*> ProteinAtoms;
extern ProteinAtoms g_proteinAtoms;

typedef vector<Cell*> Cells;
extern Cells g_cells;

typedef vector<Residue*> Residues;
extern Residues g_residues;

typedef vector<Probe*> Probes;
extern Probes g_probes;

typedef vector<Cluster*> Clusters;
extern Clusters g_clusters;

typedef map<int, Cluster*> ClusterMap;
extern ClusterMap g_clusterMap;
extern ClusterMap g_subClusterMap;

typedef vector<GridPoint*> GridVectors;
extern GridVectors g_gridVectors;

extern double g_minX, g_minY, g_minZ, g_maxX, g_maxY, g_maxZ;
extern int g_gridMaxX, g_gridMaxY, g_gridMaxZ;
extern int g_layer;
extern vector<string> g_chainList;
extern map<string, int> g_chainSizeMap;
extern double g_bcCutoffRatio;
extern double g_bcCutoffRatioG;
extern int g_maxBcNum, g_minBcNum, g_bcCutoffProbe;
extern double g_maxResidueX, g_maxResidueY, g_maxResidueZ, g_minResidueX, g_minResidueY, g_minResidueZ;
extern double g_top10Mean, g_quartile, g_bcCutoff;

typedef unsigned int PairKey;
typedef map<PairKey, float> PairDistMap;
typedef map<float, PairKey, greater<float>> DistPairMap;
typedef map<PairKey, vector<int>> Grid;

extern PairDistMap g_residuePairMap;
extern PairDistMap g_probePairMap;
extern DistPairMap g_distPairMap;

extern vector<int> g_atomCellList;
extern vector<int> g_probeCellList;
extern vector<int> g_probeBcList;

inline double calculateDistance(Vector3& a, Vector3& b) {
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

void readPdbFile(ifstream& ifs);
void readAtomPropertyFile(ifstream& ifs);
void readGridPropertyFile(ifstream& ifs);

int assignAtomPair();
void calculateProbePair();
void insertGridProteinNew(double gridSpacing, Grid& proteinGrid);
void insertGridProbeNew(int probeIndex, double gridSpacing, Grid& probeGrid);
void insertGridProtein();
void generateFirstLayer();
void generateFirstLayerNew();

void assignResiduePair();
void assignProbePair();
void displayGridProbe();
void displayGrid();

void calProbePairDist();
void generateProbe();

void generateSubFirstLayer();
void generateSubFirstLayerNew();
void generateSubSecondLayer();

void weedCutoffBc();
void weedFirstLayer();
void checkProbeNearAtoms();
void generateNextLayer();

void calculateBc();
void calculateVolume();
void checkResidueGrid();
void checkAtomGrid();
void bumpProbeWeed(vector<Probe*>& probes);

void calculateBcResidue();
void insertClusterSet();
void insertClusterSetFinal(int step);

void displayResidueCenter();
void displayProbeAtom();

int dbscanClustering(vector<int> probeIds, double eps, int minPts, int step);
int dbscanClusteringForSurface(vector<int> probeIds, double eps, int minPts);

void execDbscanClusteringFirst();
void execDbscanClusteringSec();
int execDbscanClustering(double eps, int minPts, int start); 
int execDbscanClusteringSub(double eps, int minPts, int start); 

void execClusterWeed(vector<Probe*>& probes, vector<vector<int>>& clusterSet);
void execClusterWeedRe(vector<Probe*>& probes, vector<vector<int>>& clusterSet);
void weedCluster(vector<Probe*>& probes, vector<vector<int>>& clusterSet, int cutoff);

void combineCluster();
void appendCluster();

AtomProperty* getCharge(vector<AtomProperty*>& atomProps, string residueName, string atomName);

void calculateProbeBc();
void extractProbeProperty(const Probe* probe, int index, ProbePropertySnapshot& snapshot);
bool saveProbePropertiesToJson(const string& outputPath);
void usage();
void help();

void residuePropertySetInit();
double calculatePlb(double ratio);

#endif
