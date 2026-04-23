#include "pass.h"
#include "reader.h"
#include "cal.h"
#include "soga.h"

/* global variables */
AtomPropertyMap g_atomPropertyMap;						// Atom property list
GridProperties g_gridProperties;						// Grid Cell information
ResiduePropertyMap g_sogaIndexMap,g_hydrophobicityMap ;		// Reisidue Property ( Soga index, hydrophibicity)
Residues g_residues;  				// Residue list

ResidueNameMap g_residueNames;

ProteinAtoms g_proteinAtoms;					// Protein atom list
Probes g_probes;							// Probe list
Cells g_cells;					// Grid Cell list

PairDistMap g_probePairMap;						// Probe - Probe distance (key : g_probes index)
DistPairMap g_distPairMap;						// Probe - Probe distance (key : distance)

Clusters g_clusters;						// Probe Cluster list
ClusterMap g_clusterMap;					// Probe Clusters Set (key : g_clusters index)
ClusterMap g_subClusterMap;	

Grid probe_grid_15;

GridVectors g_gridVectors;

vector <int> g_atomCellList;
vector <int> g_probeCellList;
vector <int> g_probeBcList;
vector <string> g_chainList;
map <string, int> g_chainSizeMap;


double g_minX=9999, g_minY=9999, g_minZ=9999 , g_maxX=-9999, g_maxY=-9999, g_maxZ=-9999; // Max, Min ( x,y,z) coordinate of protein atom
int g_maxBcNum = 0, g_minBcNum = 9999, g_bcCutoffProbe=0;
int g_gridMaxX, g_gridMaxY, g_gridMaxZ;
int g_layer;
double g_top10Mean, g_quartile,g_bcCutoff;
double g_bcCutoffRatioG;

int main(int argc, char *argv[]) {
	int count,probe_num;
	ifstream ifs, ifs_atom, ifs_grid; 
	string filename;
	clock_t begin,end,start;

	if( argc < 2) { 
		usage();
		return 0;
	}

	filename = argv[1]; // Input file name

	if( argc > 2) {
		g_bcCutoffRatioG = atof(argv[2]);
	} else {
		g_bcCutoffRatioG = -1;
	}
	g_bcCutoff = -1;

	start = clock();
	ifs_atom.open("atom_property"); // Open atom property file

	if(!ifs_atom) {
		cerr << "Pass" << ": could not find config file AtomProperty" <<endl;
		return 0;
	}

	readAtomPropertyFile(ifs_atom); // Read protein atom property(atom name, residue name, VDW radius, charge, polarity) from file (AtomProperty)

	ifs_grid.open("grid_property"); // Open atom property file

	if(!ifs_grid) {
		cerr << "Pass" << ": could not find config file GridProperty" <<endl;
		return 0;
	}
	readGridPropertyFile(ifs_grid); // Read protein atom property(atom name, residue name, VDW radius, charge, polarity) from file (AtomProperty)

	ifs.open(filename.c_str());

	if(!ifs) {
		cerr << "Pass" << ": could not find file '" << filename << "'." <<endl;
		return 0;
	}
	readPdbFile(ifs); // read PDB file 

	printf("REMARK Read PDB (protein (& cofactor if present) atoms) [%6d]: Min(x,y,z) : %8.3f %8.3f %8.3f Max(x,y,z) : %8.3f %8.3f %8.3f\n", g_proteinAtoms.size(),g_minX,g_minY,g_minZ,g_maxX,g_maxY,g_maxZ);

	begin = clock();
	count = initializeCellList(); 
	end = clock();
	printf("REMARK [%5.2fs] Find Near-Atom Pairs (<=7.7A) : %d pairs\n",((double)(end-begin)) / CLOCKS_PER_SEC, count);  

	residuePropertySetInit();
	begin = clock();
	generateFirstLayer();
	assignProbePair();
	end = clock();
	printf("REMARK [%5.2fs] 1st Layer Generation & Finding Near-Probe Pairs: [%d]\n",((double)(end-begin)) / CLOCKS_PER_SEC, g_probes.size());  
	
	begin = clock();
	probe_num = g_probes.size();
	execDbscanClusteringFirst();
	end = clock();
	printf("REMARK [%5.2fs] Weeding Probes/BC Calculation ( < 0.7A) %5d -> %5d\n", ((double)(end-begin)) / CLOCKS_PER_SEC, probe_num, g_probes.size());
	
	begin = clock();
	probe_num = g_probes.size();
	generateSubFirstLayer();
	end = clock();
	printf("REMARK [%5.2fs] 2nd Layer Generation/BC Calculation [%6d] : %6d -> %6d\n", ((double)(end-begin)) / CLOCKS_PER_SEC,g_probes.size()-probe_num,probe_num,g_probes.size());  
	
	begin = clock();
	probe_num = g_probes.size();
	generateSubSecondLayer();
	end = clock();
	printf("REMARK [%5.2fs] 3rd Layer Generation/BC Calculation [%6d] : %6d -> %6d\n", ((double)(end-begin)) / CLOCKS_PER_SEC,g_probes.size()-probe_num,probe_num,g_probes.size());  

	probe_num = g_probes.size();
	generateNextLayer();
	end = clock();
	printf("REMARK [%5.2fs] 4th and further Layer Generation [%6d] : %6d -> %6d \n", ((double)(end-begin)) / CLOCKS_PER_SEC,g_probes.size()-probe_num,probe_num,g_probes.size());  


	
//	execDbscanClustering(1.5,1,0); 

//	insertClusterSet();


	/*
		dist  : max density
		1.4   : 9
		1.6   : 13
		1.8   : 14
		2.0	  : 18
		2.2	  : 22
		2.4   : 27
		2.5   : 30
	*/


//	execDbscanClusteringSub(2.5,16,1);

//	calculateVolume();
//	insertClusterSetFinal(2);

//	calculatePlb(0.65); 
	end = clock();
	printf("REMARK [%5.2fs] Total Execution time for %5d Protein Atoms\n",((double)(end-start)) / CLOCKS_PER_SEC,g_proteinAtoms.size() );
	displayProbeAtom();
	if (!saveProbePropertiesToJson("probe_properties.json")) {
		cerr << "REMARK probe properties json export failed." << endl;
	}
}

void usage() {
	cout << "Usage : " << endl;
}
