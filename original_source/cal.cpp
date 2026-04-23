#include "pass.h"
#include "cal.h"
#include <cerrno>
#include <cstring>

namespace {
	string escapeJsonString(const string& input) {
		string escaped;
		escaped.reserve(input.size());
		for (size_t i = 0; i < input.size(); ++i) {
			const char ch = input[i];
			switch (ch) {
				case '\"': escaped += "\\\""; break;
				case '\\': escaped += "\\\\"; break;
				case '\b': escaped += "\\b"; break;
				case '\f': escaped += "\\f"; break;
				case '\n': escaped += "\\n"; break;
				case '\r': escaped += "\\r"; break;
				case '\t': escaped += "\\t"; break;
				default: escaped += ch; break;
			}
		}
		return escaped;
	}

	void writeIntArrayJson(ofstream& ofs, const vector<int>& values) {
		ofs << "[";
		for (size_t i = 0; i < values.size(); ++i) {
			if (i > 0) {
				ofs << ",";
			}
			ofs << values[i];
		}
		ofs << "]";
	}
}




bool compareProbeSurviving(Probe* a, Probe* b) {
	return a->isSurvived >  b->isSurvived;
}

bool compareProbeClusterId(Probe* a, Probe* b) {
	return a->clusterId1*100+a->clusterId <  b->clusterId1*100+b->clusterId;
}

bool compareProbeBc(Probe* a, Probe* b) {
	return a->numBc1 >  b->numBc1;
}

bool compareClusterSize(Cluster* a, Cluster* b) {
	return a->probeIds.size() <  b->probeIds.size();
}


bool compareNumber(int i, int j) {
	return (i<j);
}

bool compareGrid(GridPoint* a, GridPoint* b) {
	return (a->flag > b->flag);
}


double getCellId( Vector3 point, double size) {
	return ceil((point.x-g_minX)/size)+(g_gridMaxX*ceil((point.y-g_minY)/size))+(g_gridMaxX*g_gridMaxY*ceil((point.z-g_minZ)/size));
}

double getCellIdXyz( double x, double y, double z, double size) {
	return ceil((x-g_minX)/size)+(g_gridMaxX*ceil((y-g_minY)/size))+(g_gridMaxX*g_gridMaxY*ceil((z-g_minZ)/size));
}

int getCellIndex_shift ( Vector3 point,int size, int x, int y, int z) {
	return (ceil((point.z-g_minZ+10)/size)+z)+(g_gridMaxY*(ceil((point.y-g_minY+10)/size)+y))+(g_gridMaxX*(ceil((point.x-g_minX+10)/size))+x);
}

double checkProbeDist(int i , int j) {



	return vsDist(pAtom(i),pAtom(j));

	/*

	   PairDistMap::iterator pairIter;
	   pairIter = g_probePairMap.find(getPairKey(i,j));
	   if(g_probePairMap.end() == pairIter)
	   return vsDist(pAtom(i),pAtom(j));
	   else
	   return (*pairIter).second;
	 */
}

void firstClusterIdAssign() {
	for(int i=0; i<g_probes.size(); i++){
		g_probes[i]->clusterId1 = g_probes[i]->clusterId;
	}
}


void assignClusterSet() {
	ClusterMap::iterator clusterIter;
	GridVectors::iterator gridIter;
	int a,b,c,m,proteinKey,new_proteinKey;
	double dist;
	g_clusterMap.clear();

	for(int i=0; i < g_probes.size() ; i++) {
		if(g_probes[i]->clusterId >0) {
			clusterIter = g_clusterMap.find(g_probes[i]->clusterId);
			if(g_clusterMap.end() == clusterIter) {
				Cluster* temp = new Cluster;
				clusterIter = g_clusterMap.insert(ClusterMap::value_type(g_probes[i]->clusterId,temp)).first;
			} 
			clusterIter->second->clusterId = g_probes[i]->clusterId;
			clusterIter->second->probeIds.push_back(i);
			//	clusterIter->second->probeSize = clusterIter->second->probeIds.size();
		}
		//	g_probes[i]->clusterId1 = g_probes[i]->clusterId;
	}
}
void insertClusterSetSub(vector<int> probeList) {

	ClusterMap::iterator clusterIter;
	GridVectors::iterator gridIter;
	int a,b,c,i,clusterId;

	double x1=0,y1=0,z1=0;

	double dist,bcCutoff50,alpha,beta,weight,tempBcNum;
	g_subClusterMap.clear();

	for(int a=0; a < probeList.size() ; a++) {
		i = probeList[a];	
		clusterId = g_probes[i]->clusterId;
		clusterIter = g_subClusterMap.find(clusterId);
		if(g_subClusterMap.end() == clusterIter) {
			Cluster* temp = new Cluster;
			clusterIter = g_subClusterMap.insert(ClusterMap::value_type(clusterId,temp)).first;
			clusterIter->second->type = 0;
			clusterIter->second->clusterId = clusterId;
			clusterIter->second->probeIds.clear();
		} 

		if(g_probes[i]->density >= 18) {
			clusterIter->second->type++;
		}
		clusterIter->second->probeIds.push_back(i);

	}
}

void insertClusterSet() {
	ClusterMap::iterator clusterIter,clusterSubIter;
	GridVectors::iterator gridIter;
	int i,clusterId;

	double dist,bcCutoff50,alpha,beta,weight,tempBcNum;
	g_subClusterMap.clear();
	g_clusterMap.clear();

	for(i=0; i < g_probes.size() ; i++) {
		clusterId = g_probes[i]->clusterId1;

		clusterSubIter = g_subClusterMap.find(clusterId);

		if(g_subClusterMap.end() == clusterSubIter) {
			Cluster* temp = new Cluster;
			clusterSubIter = g_subClusterMap.insert(ClusterMap::value_type(clusterId,temp)).first;
		} 
		clusterSubIter->second->probeIds.push_back(i);
	}

	for(clusterSubIter = g_subClusterMap.begin(); clusterSubIter != g_subClusterMap.end() ; clusterSubIter++){

		if(clusterSubIter->second->probeIds.size() > 12 && clusterSubIter->first > 0) {
			g_clusterMap.insert(ClusterMap::value_type(clusterSubIter->first,clusterSubIter->second));
			for(i=0 ; i < clusterSubIter->second->probeIds.size() ; i++) {
				g_probes[clusterSubIter->second->probeIds[i]]->clusterId = 0;
			}
			//			cout << "REMARK cID : " << clusterSubIter->first << " : " << clusterSubIter->second->probeIds.size() << endl;
		} else {
			for(i=0 ; i < clusterSubIter->second->probeIds.size() ; i++) {
				g_probes[clusterSubIter->second->probeIds[i]]->clusterId = -1;
			}
		}
	}
}

void insertClusterSetFinal( int step) {
	ClusterMap::iterator clusterIter,clusterSubIter;

	GridVectors::iterator gridIter;
	int a,b,c,m,cellId,newCellId,clusterId,maxClusterId,maxBc=0;

	int newCid1,newCid2,tempCid;
	double dist,bcCutoff50,alpha,beta,weight,tempBcNum;
	g_clusterMap.clear();

	for(int i=0; i < g_probes.size() ; i++) {
		if(g_probes[i]->clusterId >-1) 
			clusterId = g_probes[i]->clusterId1*100 + g_probes[i]->clusterId;
		else 
			clusterId = g_probes[i]->clusterId1*100;

		clusterIter = g_clusterMap.find(clusterId);

		if(g_clusterMap.end() == clusterIter) {
			Cluster* temp = new Cluster;
			clusterIter = g_clusterMap.insert(ClusterMap::value_type(clusterId,temp)).first;
			clusterIter->second->maxBcNum = 0;
			clusterIter->second->maxLayer=0;
			clusterIter->second->minBcNum = 9999;
			clusterIter->second->sumBcNum = 0;
			clusterIter->second->sumWtBcNum =0;
			clusterIter->second->sumProbeNum =0;
			clusterIter->second->hydrophobicity = 0;
			clusterIter->second->plb = 0;
			clusterIter->second->flag = 0;
			clusterIter->second->layer1 = 0;
			clusterIter->second->layer2 = 0;
			clusterIter->second->layer3 = 0;
			clusterIter->second->layer4 = 0;
			clusterIter->second->sumBcNum4=0;
			clusterIter->second->layer1Polar = 0;
			clusterIter->second->layer2Polar = 0;
			clusterIter->second->layer3Polar = 0;
		} 


		if(g_probes[i]->numBc3 == -1) {
			clusterIter->second->bcBinTop.push_back(i);
		}

		if(g_probes[i]->numBc1 < g_bcCutoff+(g_top10Mean-g_bcCutoff)*0.1)
			clusterIter->second->bcBin1.push_back(i);
		if(g_probes[i]->numBc1 > g_bcCutoff+(g_top10Mean-g_bcCutoff)*0.1 && g_probes[i]->numBc1 <= g_bcCutoff+(g_top10Mean-g_bcCutoff)*0.2)
			clusterIter->second->bcBin2.push_back(i);
		if(g_probes[i]->numBc1 > g_bcCutoff+(g_top10Mean-g_bcCutoff)*0.2 && g_probes[i]->numBc1 <= g_bcCutoff+(g_top10Mean-g_bcCutoff)*0.3)
			clusterIter->second->bcBin3.push_back(i);
		if(g_probes[i]->numBc1 > g_bcCutoff+(g_top10Mean-g_bcCutoff)*0.3 && g_probes[i]->numBc1 <= g_bcCutoff+(g_top10Mean-g_bcCutoff)*0.4)
			clusterIter->second->bcBin4.push_back(i);
		if(g_probes[i]->numBc1 > g_bcCutoff+(g_top10Mean-g_bcCutoff)*0.4 && g_probes[i]->numBc1 <= g_bcCutoff+(g_top10Mean-g_bcCutoff)*0.5)
			clusterIter->second->bcBin5.push_back(i);
		if(g_probes[i]->numBc1 > g_bcCutoff+(g_top10Mean-g_bcCutoff)*0.5 && g_probes[i]->numBc1 <= g_bcCutoff+(g_top10Mean-g_bcCutoff)*0.6)
			clusterIter->second->bcBin6.push_back(i);
		if(g_probes[i]->numBc1 > g_bcCutoff+(g_top10Mean-g_bcCutoff)*0.6 && g_probes[i]->numBc1 <= g_bcCutoff+(g_top10Mean-g_bcCutoff)*0.7)
			clusterIter->second->bcBin7.push_back(i);
		if(g_probes[i]->numBc1 > g_bcCutoff+(g_top10Mean-g_bcCutoff)*0.7 && g_probes[i]->numBc1 <= g_bcCutoff+(g_top10Mean-g_bcCutoff)*0.8)
			clusterIter->second->bcBin8.push_back(i);
		if(g_probes[i]->numBc1 > g_bcCutoff+(g_top10Mean-g_bcCutoff)*0.8 && g_probes[i]->numBc1 <= g_bcCutoff+(g_top10Mean-g_bcCutoff)*0.9)
			clusterIter->second->bcBin9.push_back(i);
		if(g_probes[i]->numBc1 > g_bcCutoff+(g_top10Mean-g_bcCutoff)*0.9)
			clusterIter->second->bcBin10.push_back(i);



		clusterIter->second->clusterId = clusterId;
		clusterIter->second->probeIds.push_back(i);
		clusterIter->second->maxLayer = MAX(pNLayer(i),clusterIter->second->maxLayer);
		clusterIter->second->probeSize = clusterIter->second->probeIds.size();
		//		clusterIter->second->maxBcNum = ( clusterIter->second->maxBcNum < g_probes[i]->numBc1  ) ? g_probes[i]->numBc1 : clusterIter->second->maxBcNum;
		clusterIter->second->minBcNum = ( clusterIter->second->minBcNum > g_probes[i]->numBc1 ) ? g_probes[i]->numBc1 : clusterIter->second->minBcNum;
		clusterIter->second->sumBcNum += g_probes[i]->numBc1;
		clusterIter->second->sumProbeNum++;

		if(clusterIter->second->maxBcNum < g_probes[i]->numBc1 ) {
			clusterIter->second->maxBcNum = g_probes[i]->numBc1;
			clusterIter->second->maxBcId = i;
		}


		switch (pNLayer(i)) {
			case 1 : clusterIter->second->layer1++; 
					 if(g_probes[i]->isPolar >0 || g_probes[i]->closestDist < 3.2001) 
						 clusterIter->second->layer1Polar++;
					 break;
			case 2 : clusterIter->second->layer2++; 
					 if(g_probes[i]->isPolar >0 || g_probes[i]->closestDist < 3.2001)
						 clusterIter->second->layer2Polar++;
					 break;
			case 3 : clusterIter->second->layer3++; 
					 if(g_probes[i]->isPolar >0 || g_probes[i]->closestDist < 3.2001)
						 clusterIter->second->layer3Polar++;
					 break;
			default : clusterIter->second->layer4++; 
					  clusterIter->second->sumBcNum4 += g_probes[i]->numBc1;
					  break;

		}


		cellId = getCellId(g_probes[i]->point,2);

		for(int k=0; k< 125; k++) {
			newCellId = cellId + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);
			if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ) {
				if(g_cells[newCellId]->atoms.size() >0) {
					for(int l=0; l < g_cells[newCellId]->atoms.size(); l++) {
						dist = vsDist(g_probes[i]->point,prAtom(g_cells[newCellId]->atoms[l]));
						if(dist <= 3.8001) {
							if( (ResiName(g_cells[newCellId]->atoms[l]) == "GLY" && prType(g_cells[newCellId]->atoms[l])  == "CA") ||
									(ResiName(g_cells[newCellId]->atoms[l]) == "ALA" && (prType(g_cells[newCellId]->atoms[l]) == "CA" || prType(g_cells[newCellId]->atoms[l]) == "CB")) ||
									(ResiName(g_cells[newCellId]->atoms[l]) == "PRO" && (prType(g_cells[newCellId]->atoms[l]) == "CA" || prType(g_cells[newCellId]->atoms[l]) == "CB" || prType(g_cells[newCellId]->atoms[l]) == "CG" || prType(g_cells[newCellId]->atoms[l]) == "CD")) ||
									(prType(g_cells[newCellId]->atoms[l]) != "C" && prType(g_cells[newCellId]->atoms[l]) != "N" && prType(g_cells[newCellId]->atoms[l]) != "O" && prType(g_cells[newCellId]->atoms[l]) != "CA" && prType(g_cells[newCellId]->atoms[l]) != "CB" )
							  ) {
								clusterIter->second->residues.push_back(ResiIdx(g_cells[newCellId]->atoms[l]));
							}
							clusterIter->second->residuesContact.push_back(ResiIdx(g_cells[newCellId]->atoms[l]));

						}
					}
				}
			}
		}
	}


	for(int i=0 ; i< g_gridVectors.size(); i++ ) {
			for(int j=0 ; j< g_gridVectors[i]->probes.size(); j++) {
				clusterIter = g_clusterMap.find(g_probes[g_gridVectors[i]->probes[j]]->clusterId1*100+g_probes[g_gridVectors[i]->probes[j]]->clusterId);
				if(clusterIter != g_clusterMap.end()) {
					if(g_gridVectors[i]->flag == 5) 
						clusterIter->second->volCount2.push_back(i);

					if(g_gridVectors[i]->flag > 10)
						clusterIter->second->volCount3.push_back(i);

					if(g_gridVectors[i]->flag == 4) 
						clusterIter->second->volCount4.push_back(i);

					clusterIter->second->volCount1.push_back(i);
				}
			}
	}

		/*
		   for(int i=0 ; i< g_gridVectors.size(); i++ ) {

		   if(g_gridVectors[i]->flag == 1 || g_gridVectors[i]->flag > 3) {
		   for(int j=0 ; j <g_gridVectors[i]->probes.size() ; j++) {
		   clusterIter = g_clusterMap.find(g_probes[g_gridVectors[i]->probes[j]]->clusterId1*100+g_probes[g_gridVectors[i]->probes[j]]->clusterId);
		   if(g_gridVectors[i]->dist[j] <= 1.4 && clusterIter != g_clusterMap.end()) {
		   clusterIter->second->volCount1.push_back(i);

		   }
		   }
		   }

		   if(g_gridVectors[i]->flag % 10 == 5 ) {
		   for(int j=0 ; j <g_gridVectors[i]->probes.size() ; j++) {
		   clusterIter = g_clusterMap.find(g_probes[g_gridVectors[i]->probes[j]]->clusterId1*100+g_probes[g_gridVectors[i]->probes[j]]->clusterId);
		   if( g_gridVectors[i]->dist[j] <= 1.4 && clusterIter != g_clusterMap.end()) {
		   clusterIter->second->volCount3.push_back(i);

		   }
		   }
		   }

		   if(g_gridVectors[i]->flag > 10) {
		   for(int j=0 ; j <g_gridVectors[i]->probes.size() ; j++) {
		   clusterIter = g_clusterMap.find(g_probes[g_gridVectors[i]->probes[j]]->clusterId1*100+g_probes[g_gridVectors[i]->probes[j]]->clusterId);
		   if(g_gridVectors[i]->dist[j] <= 1.4 && clusterIter != g_clusterMap.end()) {
		   clusterIter->second->volCount2.push_back(i);
		   }
		   }

		   }
		   }
		 */


		for(clusterIter = g_clusterMap.begin(); clusterIter != g_clusterMap.end() ; clusterIter++){
			vector<int>::iterator it;
			sort(clusterIter->second->residues.begin(),clusterIter->second->residues.end(),compareNumber);
			it = unique (clusterIter->second->residues.begin(), clusterIter->second->residues.end()); 
			clusterIter->second->residues.resize( it - clusterIter->second->residues.begin() );    

			sort(clusterIter->second->residuesContact.begin(),clusterIter->second->residuesContact.end(),compareNumber);
			it = unique (clusterIter->second->residuesContact.begin(), clusterIter->second->residuesContact.end()); 
			clusterIter->second->residuesContact.resize( it - clusterIter->second->residuesContact.begin() );    

			clusterIter->second->averageBcNum = (double)clusterIter->second->sumBcNum /(double) clusterIter->second->sumProbeNum ;
			if(clusterIter->second->layer4>0) {
				clusterIter->second->averageBcNum4  = (double)clusterIter->second->sumBcNum4 / clusterIter->second->layer4;
			}




			double tempSum=0,tempSum4=0,tempWtSum=0;

			vector <int> topProbes, bottomProbes, midProbes;
			for(int i=0; i< clusterIter->second->probeIds.size() ; i++) {
				tempSum += SQUARE(g_probes[clusterIter->second->probeIds[i]]->numBc1 - clusterIter->second->averageBcNum);
				if(pNLayer(clusterIter->second->probeIds[i]) >3) {
					tempSum4 += SQUARE(g_probes[clusterIter->second->probeIds[i]]->numBc1 - clusterIter->second->averageBcNum4);
				}


				if(g_probes[clusterIter->second->probeIds[i]]->numBc1 >  clusterIter->second->maxBcNum - (clusterIter->second->maxBcNum - clusterIter->second->minBcNum ) * 0.05) 
					bottomProbes.push_back(clusterIter->second->probeIds[i]);


			}


			double openingNumber=0;


			clusterIter->second->binNumber = 0;
			if(clusterIter->second->bcBin10.size() > 0) {
				clusterIter->second->type = 10; 
				topProbes.assign( clusterIter->second->bcBin10.begin(), clusterIter->second->bcBin10.end() ); 
				clusterIter->second->binStart = 10;
				clusterIter->second->binNumber++;
				
			}
			if(clusterIter->second->bcBin9.size() > 0) {
				clusterIter->second->type = 9; 
				topProbes.assign( clusterIter->second->bcBin9.begin(), clusterIter->second->bcBin9.end() );
				clusterIter->second->binStart = 9;
				clusterIter->second->binNumber++;

			}
			if(clusterIter->second->bcBin8.size() > 0) {
				clusterIter->second->type = 8;
				topProbes.assign( clusterIter->second->bcBin8.begin(), clusterIter->second->bcBin8.end() );
				clusterIter->second->binStart = 8;
				clusterIter->second->binNumber++;

			}
			if(clusterIter->second->bcBin7.size() > 0){
				clusterIter->second->type = 7;
				topProbes.assign( clusterIter->second->bcBin7.begin(), clusterIter->second->bcBin7.end() );
				clusterIter->second->binStart = 7;
				clusterIter->second->binNumber++;

			}
			if(clusterIter->second->bcBin6.size() > 0) {
				clusterIter->second->type = 6;
				topProbes.assign( clusterIter->second->bcBin6.begin(), clusterIter->second->bcBin6.end() );
				clusterIter->second->binStart = 6;
				clusterIter->second->binNumber++;

			}
			if(clusterIter->second->bcBin5.size() > 0) {
				clusterIter->second->type = 5;
				topProbes.assign( clusterIter->second->bcBin5.begin(), clusterIter->second->bcBin5.end() );
				clusterIter->second->binStart = 5;
				clusterIter->second->binNumber++;

			}
			if(clusterIter->second->bcBin4.size() > 0) {
				clusterIter->second->type = 4;
				topProbes.assign( clusterIter->second->bcBin4.begin(), clusterIter->second->bcBin4.end() );
				clusterIter->second->binStart = 4;
				clusterIter->second->binNumber++;

			}
			if(clusterIter->second->bcBin3.size() > 0){
				clusterIter->second->type = 3;
				topProbes.assign( clusterIter->second->bcBin3.begin(), clusterIter->second->bcBin3.end() );
				clusterIter->second->binStart = 3;
				clusterIter->second->binNumber++;

			}
			if(clusterIter->second->bcBin2.size() > 0){
				clusterIter->second->type = 2;
				topProbes.assign( clusterIter->second->bcBin2.begin(), clusterIter->second->bcBin2.end() );
				clusterIter->second->binStart = 2;
				clusterIter->second->binNumber++;

			}
			if(clusterIter->second->bcBin1.size() > 0) {
				clusterIter->second->type = 1;
				clusterIter->second->binStart = 1;
				clusterIter->second->binNumber++;

				topProbes.assign( clusterIter->second->bcBin1.begin(), clusterIter->second->bcBin1.end() );
			}

			openingNumber = dbscanClusteringForSurface(clusterIter->second->bcBinTop,4,1);


			double minDist=9999;

			if(openingNumber >0) {
				for(int i=0; i<clusterIter->second->bcBinTop.size();i++) {
					dist = vsDist(pAtom(clusterIter->second->bcBinTop[i]),pAtom(clusterIter->second->maxBcId));
					minDist = (minDist < dist) ? minDist : dist;
				}
			} else {
				for(int i=0; i<topProbes.size();i++) {
					dist = vsDist(pAtom(topProbes[i]),pAtom(clusterIter->second->maxBcId));
					minDist = (minDist < dist) ? minDist : dist;
				}

			}



			double maxDist=-9999;
			for(int i=0 ; i<  clusterIter->second->probeIds.size(); i++) {
				for(int j=i+1 ; j< clusterIter->second->probeIds.size(); j++) {
					dist = vsDist(pAtom( clusterIter->second->probeIds[i]),pAtom( clusterIter->second->probeIds[j]));
					maxDist = (maxDist > dist) ? maxDist : dist;
				}
			}


			clusterIter->second->depth=minDist;
			clusterIter->second->width=maxDist;
			clusterIter->second->opening = openingNumber;




			clusterIter->second->stdBcNum = sqrt(tempSum/clusterIter->second->probeIds.size());
			if(clusterIter->second->layer4>0) {
				clusterIter->second->stdBcNum4= sqrt(tempSum4/clusterIter->second->layer4);
			}

		}

		newCid1 = 0;
		newCid2 = 1;
		tempCid = -1;
		g_subClusterMap.clear();
		for(clusterIter = g_clusterMap.begin(); clusterIter != g_clusterMap.end() ; clusterIter++){
			if( clusterIter->second->probeIds.size() > 12 && clusterIter->first !=0 ) {
				if(tempCid != clusterIter->first/100) { 
					newCid1++;
					tempCid = clusterIter->first/100;

				}
				clusterSubIter = g_subClusterMap.insert(ClusterMap::value_type(newCid1*100 + clusterIter->first%100,clusterIter->second)).first;
				for(int i=0; i<clusterIter->second->probeIds.size() ; i++) {
					g_probes[clusterIter->second->probeIds[i]]->clusterId1 = newCid1;
				}

			} else {

				//			clusterSubIter = g_subClusterMap.insert(ClusterMap::value_type(newCid2,clusterIter->second)).first;
				for(int i=0; i<clusterIter->second->probeIds.size() ; i++) {
					g_probes[clusterIter->second->probeIds[i]]->clusterId1 = newCid2;
					g_probes[clusterIter->second->probeIds[i]]->clusterId = -1;
				}

				newCid2++;

			}
			//		cout << "REMARK CID : " << clusterIter->first << " : " <<clusterIter->second->probeIds.size() <<" : "<< clusterIter->second->flag << endl;
		}

		g_clusterMap.clear();
		for(clusterSubIter = g_subClusterMap.begin(); clusterSubIter != g_subClusterMap.end() ; clusterSubIter++){

			g_clusterMap.insert(ClusterMap::value_type(clusterSubIter->first,clusterSubIter->second));
			clusterSubIter->second->clusterId  = clusterSubIter->first;
			//		cout << "REMARK New CID : " << clusterSubIter->first << " : " <<clusterSubIter->second->probeIds.size() <<" : "<< clusterSubIter->second->flag << endl;
		}



	}
	int checkBC(double x,double y, double z) {
		int cellId,newCellId;
		double dist,cutoff,cutoff_probe;
		double totalDist=0,minDist=9999;

		int bc,range,countTotal=0,calTotal=0,averageNum=0;

		vector <int> tempAtoms,tempCell,tempProbes;

		cellId = getCellIdXyz(x,y,z,2);

		range = g_gridProperties.size();

		for(int k=0; k< range; k++) {
			newCellId = cellId + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);
			if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ) {
				if(g_gridProperties[k]->max <= 14) {
					countTotal += g_cells[newCellId]->atoms.size();
				} else {
					if(g_cells[newCellId]->atoms.size() >0) { 
						for(int l=0; l < g_cells[newCellId]->atoms.size(); l++) {
							dist = distPoint(x,y,z,prAtom(g_cells[newCellId]->atoms[l]));
							if(dist <= 14) {
								calTotal++;
							}
						}
					}
				}
			}
		}

		return countTotal+calTotal;
	}

	double getAngle(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3) {
		double v1, v2, v1_mag, v2_mag, v1_norm_x,v1_norm_y,v1_norm_z, v2_norm_x,v2_norm_y,v2_norm_z, res;
		double xv1,yv1,zv1,xv2,yv2,zv2;

		xv1 = x2-x1;
		yv1 = y2-y1;
		zv1 = z2-z1;

		xv2 = x3-x1;
		yv2 = y3-y1;
		zv2 = z3-z1;

		v1_mag = sqrt(xv1*xv1 + yv1*yv1 + zv1*zv1);
		v2_mag = sqrt(xv2*xv2 + yv2*yv2 + zv2*zv2);

		v1_norm_x = xv1 / v1_mag;
		v1_norm_y = yv1 / v1_mag;
		v1_norm_z = zv1 / v1_mag;
		v2_norm_x = xv2 / v2_mag;
		v2_norm_y = yv2 / v2_mag;
		v2_norm_z = zv2 / v2_mag;



		res = v1_norm_x * v2_norm_x + v1_norm_y * v2_norm_y + v1_norm_z * v2_norm_z;

		return acos(res)*180/PI;







	}



	void calculateVolume() {
		double x,y,z,dist,distMinProbe,distMinAtom,cutoff,checkDist;
		x = (int)g_minX;
		int idx = 0,type=0,id,resi,probeNum,atomNum,probeSurNum,cellCount,total_max=0;
		int cellId, newCellId;
		int targetCid,targetCid1;
		int bc,flag,flagProbe,flagProbeCount,flagAtom,flagAtomCount,flagVacant,flagVacantProtein,flagVacantCount,flagCount;


		int minProbe,minAtom;
		vector <int> probeIds;
		vector <double> probeDist;

		ClusterMap::iterator clusterIter;

		GridPoint* tempGrid;
		g_gridVectors.clear();

		while(x < g_maxX+1) {
			y=(int)g_minY;
			while(y < g_maxY+1) {
				z = (int)g_minZ;
				while(z < g_maxZ+1) {
					id = idx%9999;
					resi = idx/9999;
					probeNum = 0;
					probeSurNum=0;
					atomNum = 0;
					distMinProbe = 9999;
					distMinAtom = 9999;
					cellId = getCellIdXyz(x,y,z,2);
					cellCount=0;
					type = -1;
					flag = 0;

					probeIds.clear();
					probeDist.clear();

					tempGrid = new GridPoint;


					for(int k=0; k<27 ; k++) {
						newCellId = cellId + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);

						if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ) {
							if(g_cells[newCellId]->probes.size() >0) { 
								for(int l=0; l < g_cells[newCellId]->probes.size(); l++) {
									dist = distPoint(x,y,z,pAtom(g_cells[newCellId]->probes[l]));

									/*								if(distMinProbe > dist) {
																	distMinProbe = dist;
																	minProbe = g_cells[newCellId]->probes[l];
																	} */

									if(dist <= g_probes[g_cells[newCellId]->probes[l]]->closestDist/2 || dist < 1.2) {

										tempGrid->probes.push_back(g_cells[newCellId]->probes[l]);
										probeNum++;

									}


								}
							}
						}
					}

					/*
					   if(probeNum >0) {
					   for(int k=0; k<125 ; k++) {
					   newCellId = cellId + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);
					   if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ) {
					   if(g_cells[newCellId]->atoms.size() >0) {
					   for(int l=0; l < g_cells[newCellId]->atoms.size(); l++) {
					   dist = distPoint(x,y,z,prAtom(g_cells[newCellId]->atoms[l]));
					   if(dist < prRadius(g_cells[newCellId]->atoms[l]))  
					   atomNum++;

					   if(distMinAtom > dist) {
					   distMinAtom = dist;
					   minAtom = g_cells[newCellId]->atoms[l];
					   }
					   }
					   }
					   }
					   }
					   }
					 */





					if(probeNum >0) { // && atomNum==0) {
						//double angle = getAngle(x,y,z,pAtom(minProbe).x,pAtom(minProbe).y,pAtom(minProbe).z,prAtom(minAtom).x,prAtom(minAtom).y,prAtom(minAtom).z);
						//	if(angle >90 ) {

						//	dist = vsDist(pAtom(minProbe),prAtom(minAtom));

						//	cout << "REMARK min dist : " << idx << " :: " << distMinProbe << " ( " << minProbe << " ) : " << distMinAtom << " ( " << minAtom << " ) " << " : " << dist << " : " << angle << endl;
						tempGrid->point.x = x;
						tempGrid->point.y = y;
						tempGrid->point.z = z;
						tempGrid->flag = 1;

						g_gridVectors.push_back(tempGrid);


						idx++;
						//	}
					}

					z += 1;
					}
					y += 1;
				}
				x += 1;
			}

			idx = 0;

			//	sort (g_gridVectors.begin(),g_gridVectors.end(),compareGrid);
			cout << "REMARK Total Grid Point : " << g_gridVectors.size() << endl;
			for(int i=0 ; i< g_gridVectors.size(); i++ ) {

				if(g_gridVectors[i]->flag == 3) {
					for(int j=0 ; j< g_gridVectors[i]->probes.size(); j++) {
						if( g_gridVectors[i]->dist[j] <= 2.1) {
							g_probes[g_gridVectors[i]->probes[j]]->numBc3 = -1;
						}
					}
				}

				if(g_gridVectors[i]->flag == 1) {

					cellId = getCellIdXyz(g_gridVectors[i]->point.x+1,g_gridVectors[i]->point.y,g_gridVectors[i]->point.z,2);

					flagProbe = 0;
					flagAtom = 0;
					flagProbeCount = 0;
					flagAtomCount = 0;
					flagVacantCount = 0;
					flagVacant = 0;
					flagVacantProtein = 0;
					flagCount =0;

					for(int k=0; k<125 ; k++) {
						newCellId = cellId + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);

						if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ) {
							if(g_cells[newCellId]->probes.size() >0) { 
								for(int l=0; l < g_cells[newCellId]->probes.size(); l++) {


									dist = distPoint(g_gridVectors[i]->point.x+1,g_gridVectors[i]->point.y,g_gridVectors[i]->point.z,pAtom(g_cells[newCellId]->probes[l]));

									if(dist <= g_probes[g_cells[newCellId]->probes[l]]->closestDist/2 || dist < 1.2) 
										flagProbe=1;



								}
							}
							if(g_cells[newCellId]->atoms.size() >0) {
								for(int l=0; l < g_cells[newCellId]->atoms.size(); l++) {
									dist = distPoint(g_gridVectors[i]->point.x+1,g_gridVectors[i]->point.y,g_gridVectors[i]->point.z,prAtom(g_cells[newCellId]->atoms[l]));

									if(dist < prRadius(g_cells[newCellId]->atoms[l])+0.7)  {
										flagAtom = 1;
									}
									if( dist < prRadius(g_cells[newCellId]->atoms[l])+1.2 )  {
										flagVacant =1;
									}
								}
							}
						}
					}

					if( flagProbe == 1) { 
						flagProbeCount++;
						flagVacant = 0;
					}

					if( flagAtom ==1)
						flagAtomCount++;

					if( flagProbe == 0 && flagVacant ==0)
						flagVacantCount++;

					flagProbe = 0;
					flagAtom = 0;
					flagVacant =0;
					flagVacantProtein = 0;

					cellId = getCellIdXyz(g_gridVectors[i]->point.x-1,g_gridVectors[i]->point.y,g_gridVectors[i]->point.z,2);

					for(int k=0; k<125 ; k++) {
						newCellId = cellId + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);

						if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ) {
							if(g_cells[newCellId]->probes.size() >0) { 
								for(int l=0; l < g_cells[newCellId]->probes.size(); l++) {
									dist = distPoint(g_gridVectors[i]->point.x-1,g_gridVectors[i]->point.y,g_gridVectors[i]->point.z,pAtom(g_cells[newCellId]->probes[l]));
									if(dist <= g_probes[g_cells[newCellId]->probes[l]]->closestDist/2 || dist < 1.2) 
										flagProbe=1;


								}
							}
							if(g_cells[newCellId]->atoms.size() >0) {
								for(int l=0; l < g_cells[newCellId]->atoms.size(); l++) {
									dist = distPoint(g_gridVectors[i]->point.x-1,g_gridVectors[i]->point.y,g_gridVectors[i]->point.z,prAtom(g_cells[newCellId]->atoms[l]));

									if(dist < prRadius(g_cells[newCellId]->atoms[l])+0.7)  {
										flagAtom = 1;
									}
									if( dist < prRadius(g_cells[newCellId]->atoms[l])+1.2)  {
										flagVacant =1;
									}
								}
							}
						}
					}

					if( flagProbe == 1) { 
						flagProbeCount++;
						flagVacant = 0;
					}

					if(flagAtom ==1)
						flagAtomCount++;

					if( flagProbe == 0 && flagVacant ==0)
						flagVacantCount++;

					flagProbe = 0;
					flagAtom = 0;
					flagVacant =0;
					flagVacantProtein=0;

					cellId = getCellIdXyz(g_gridVectors[i]->point.x,g_gridVectors[i]->point.y+1,g_gridVectors[i]->point.z,2);

					for(int k=0; k<125 ; k++) {
						newCellId = cellId + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);

						if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ) {
							if(g_cells[newCellId]->probes.size() >0) { 
								for(int l=0; l < g_cells[newCellId]->probes.size(); l++) {
									dist = distPoint(g_gridVectors[i]->point.x,g_gridVectors[i]->point.y+1,g_gridVectors[i]->point.z,pAtom(g_cells[newCellId]->probes[l]));
									if(dist <= g_probes[g_cells[newCellId]->probes[l]]->closestDist/2 || dist < 1.2) 
										flagProbe=1;					
								}
							}
							if(g_cells[newCellId]->atoms.size() >0) {
								for(int l=0; l < g_cells[newCellId]->atoms.size(); l++) {
									dist = distPoint(g_gridVectors[i]->point.x,g_gridVectors[i]->point.y+1,g_gridVectors[i]->point.z,prAtom(g_cells[newCellId]->atoms[l]));


									if(dist < prRadius(g_cells[newCellId]->atoms[l])+0.7)  {
										flagAtom = 1;
									}
									if( dist < prRadius(g_cells[newCellId]->atoms[l])+1.2)  {
										flagVacant =1;
									}
								}
							}
						}
					}

					if( flagProbe == 1) { 
						flagProbeCount++;
						flagVacant = 0;
					}

					if( flagAtom ==1)
						flagAtomCount++;

					if( flagProbe == 0 && flagVacant ==0)
						flagVacantCount++;

					flagProbe = 0;
					flagAtom = 0;
					flagVacant =0;
					flagVacantProtein = 0;

					cellId = getCellIdXyz(g_gridVectors[i]->point.x,g_gridVectors[i]->point.y-1,g_gridVectors[i]->point.z,2);

					for(int k=0; k<125 ; k++) {
						newCellId = cellId + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);

						if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ) {
							if(g_cells[newCellId]->probes.size() >0) { 
								for(int l=0; l < g_cells[newCellId]->probes.size(); l++) {
									dist = distPoint(g_gridVectors[i]->point.x,g_gridVectors[i]->point.y-1,g_gridVectors[i]->point.z,pAtom(g_cells[newCellId]->probes[l]));
									if(dist <= g_probes[g_cells[newCellId]->probes[l]]->closestDist/2 || dist < 1.2) 
										flagProbe=1;			
								}
							}
							if(g_cells[newCellId]->atoms.size() >0) {
								for(int l=0; l < g_cells[newCellId]->atoms.size(); l++) {
									dist = distPoint(g_gridVectors[i]->point.x,g_gridVectors[i]->point.y-1,g_gridVectors[i]->point.z,prAtom(g_cells[newCellId]->atoms[l]));


									if(dist < prRadius(g_cells[newCellId]->atoms[l])+0.7)  {
										flagAtom = 1;
									}
									if( dist < prRadius(g_cells[newCellId]->atoms[l])+1.2)  {
										flagVacant =1;
									}

								}
							}
						}
					}

					if( flagProbe == 1) { 
						flagProbeCount++;
						flagVacant = 0;
					}
					if(flagAtom ==1)
						flagAtomCount++;

					if( flagProbe == 0 && flagVacant ==0)
						flagVacantCount++;

					flagProbe = 0;
					flagAtom = 0;
					flagVacant =0;
					flagVacantProtein =0;

					cellId = getCellIdXyz(g_gridVectors[i]->point.x,g_gridVectors[i]->point.y,g_gridVectors[i]->point.z+1,2);



					for(int k=0; k<125 ; k++) {
						newCellId = cellId + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);

						if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ) {
							if(g_cells[newCellId]->probes.size() >0) { 
								for(int l=0; l < g_cells[newCellId]->probes.size(); l++) {
									dist = distPoint(g_gridVectors[i]->point.x,g_gridVectors[i]->point.y,g_gridVectors[i]->point.z+1,pAtom(g_cells[newCellId]->probes[l]));
									if(dist <= g_probes[g_cells[newCellId]->probes[l]]->closestDist/2 || dist < 1.2) 
										flagProbe=1;
								}
							}
							if(g_cells[newCellId]->atoms.size() >0) {
								for(int l=0; l < g_cells[newCellId]->atoms.size(); l++) {
									dist = distPoint(g_gridVectors[i]->point.x,g_gridVectors[i]->point.y,g_gridVectors[i]->point.z+1,prAtom(g_cells[newCellId]->atoms[l]));

									if(dist < prRadius(g_cells[newCellId]->atoms[l])+0.7)  {
										flagAtom = 1;
									}
									if( dist < prRadius(g_cells[newCellId]->atoms[l])+1.2)  {
										flagVacant =1;
									}

								}
							}
						}
					}


					if( flagProbe == 1) { 
						flagProbeCount++;
						flagVacant = 0;
					}

					if(flagAtom ==1)
						flagAtomCount++;

					if( flagProbe == 0 && flagVacant ==0)
						flagVacantCount++;

					flagProbe = 0;
					flagAtom = 0;
					flagVacant =0;
					flagVacantProtein=0;

					cellId = getCellIdXyz(g_gridVectors[i]->point.x,g_gridVectors[i]->point.y,g_gridVectors[i]->point.z-1,2);


					for(int k=0; k<125 ; k++) {
						newCellId = cellId + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);

						if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ) {
							if(g_cells[newCellId]->probes.size() >0) { 
								for(int l=0; l < g_cells[newCellId]->probes.size(); l++) {
									dist = distPoint(g_gridVectors[i]->point.x,g_gridVectors[i]->point.y,g_gridVectors[i]->point.z-1,pAtom(g_cells[newCellId]->probes[l]));
									if(dist <= g_probes[g_cells[newCellId]->probes[l]]->closestDist/2 || dist < 1.2) 
										flagProbe=1;
								}
							}
							if(g_cells[newCellId]->atoms.size() >0) {
								for(int l=0; l < g_cells[newCellId]->atoms.size(); l++) {
									dist = distPoint(g_gridVectors[i]->point.x+1,g_gridVectors[i]->point.y,g_gridVectors[i]->point.z-1,prAtom(g_cells[newCellId]->atoms[l]));

									if(dist < prRadius(g_cells[newCellId]->atoms[l])+0.7)  {
										flagAtom = 1;
									}
									if( dist < prRadius(g_cells[newCellId]->atoms[l])+1.2)  {
										flagVacant =1;
									}

								}
							}
						}
					}

					if( flagProbe == 1) { 
						flagProbeCount++;
						flagVacant = 0;
					}

					if( flagAtom ==1)
						flagAtomCount++;

					if( flagProbe == 0 && flagVacant ==0)
						flagVacantCount++;




					if(flagProbeCount == 6) {
						g_gridVectors[i]->flag = 4;
					} 
					else {
						if(flagAtomCount >0) {
							g_gridVectors[i]->flag = 5;
						}
						if(flagVacantCount > 0 || g_gridVectors[i]->flag ==1) {
							g_gridVectors[i]->flag += 10;
						}

					}
				}

				flagVacant = 0;
				if(g_gridVectors[i]->flag > 10) {
					for(int j=0 ; j< g_gridVectors[i]->probes.size(); j++) {
							if(g_probes[g_gridVectors[i]->probes[j]]->numBc1 > g_bcCutoff + g_top10Mean * 0.1) {
								flagVacant = 1;
							} 
					}
					if(flagVacant == 1) {
						g_gridVectors[i]->flag -= 10;
					} else {
						for(int j=0 ; j< g_gridVectors[i]->probes.size(); j++) {
								g_probes[g_gridVectors[i]->probes[j]]->numBc3 = -1;
						}
					}
			}
				id = idx % 9999;
				if(g_gridVectors[i]->flag ) {

					/*
					   for(int j=0 ; j< g_gridVectors[i]->probes.size(); j++) {
					   clusterIter = g_clusterMap.find(g_probes[g_gridVectors[i]->probes[j]]->clusterId1*100+g_probes[g_gridVectors[i]->probes[j]]->clusterId);
					   if(clusterIter != g_clusterMap.end()) {
					   if(g_gridVectors[i]->flag == 5) 
					   clusterIter->second->volCount2.push_back(i);

					   if(g_gridVectors[i]->flag > 10)
					   clusterIter->second->volCount3.push_back(i);

					   clusterIter->second->volCount1.push_back(i);
					   }
					   }
					 */


					printf("ATOM  %5d  H   VOL  %4d    %8.3f%8.3f%8.3f                                  %4d %4d %4d %5d\n",
							id,
							g_gridVectors[i]->flag,
							g_gridVectors[i]->point.x,
							g_gridVectors[i]->point.y,
							g_gridVectors[i]->point.z,

							flagProbeCount,
							flagAtomCount,
							flagVacantCount,
							g_gridVectors[i]->probes.size()



						  );

					idx++;


					if(g_gridVectors[i]->flag == 100 ) {
						printf("ATOM  %5d  H   VOL  %4d    %8.3f%8.3f%8.3f      %4d\n",
								id,
								100,
								g_gridVectors[i]->point.x+1,
								g_gridVectors[i]->point.y,
								g_gridVectors[i]->point.z,
								flagVacantCount
							  );

						idx++;
						printf("ATOM  %5d  H   VOL  %4d    %8.3f%8.3f%8.3f      %4d\n",
								id,
								100,
								g_gridVectors[i]->point.x-1,
								g_gridVectors[i]->point.y,
								g_gridVectors[i]->point.z,
								flagVacantCount
							  );

						idx++;
						printf("ATOM  %5d  H   VOL  %4d    %8.3f%8.3f%8.3f      %4d\n",
								id,
								100,
								g_gridVectors[i]->point.x,
								g_gridVectors[i]->point.y+1,
								g_gridVectors[i]->point.z,
								flagVacantCount
							  );

						idx++;
						printf("ATOM  %5d  H   VOL  %4d    %8.3f%8.3f%8.3f      %4d\n",
								id,
								100,
								g_gridVectors[i]->point.x,
								g_gridVectors[i]->point.y-1,
								g_gridVectors[i]->point.z,
								flagVacantCount
							  );

						idx++;
						printf("ATOM  %5d  H   VOL  %4d    %8.3f%8.3f%8.3f      %4d\n",
								id,
								100,
								g_gridVectors[i]->point.x,
								g_gridVectors[i]->point.y,
								g_gridVectors[i]->point.z+1,
								flagVacantCount
							  );

						idx++;
						printf("ATOM  %5d  H   VOL  %4d    %8.3f%8.3f%8.3f      %4d\n",
								id,
								100,
								g_gridVectors[i]->point.x,
								g_gridVectors[i]->point.y,
								g_gridVectors[i]->point.z-1,
								flagVacantCount
							  );

						idx++;

					}
				}
			}
			cout << "END " << endl;
		}



		void checkProbeBump(int probe_id) {
			int key, newKey,size=1;
			double dist,cutoff=0.7;
			int cellId, newCellId;
			cellId = getCellId(pAtom(probe_id),2);

			for(int k=0; k< 27; k++) {
				newCellId = cellId + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);
				if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ)
					if(g_cells[newCellId]->probes.size() >0) { 
						for(int l=0; l < g_cells[newCellId]->probes.size(); l++) {
							dist = vsDist(pAtom(probe_id),pAtom(g_cells[newCellId]->probes[l]));
							if(dist < cutoff-0.0001 && probe_id != g_cells[newCellId]->probes[l])  {
								g_probes[g_cells[newCellId]->probes[l]]->clusterId = -4;
								pIsSurvived(g_cells[newCellId]->probes[l]) = -1;
							}
						}
					}
			}
		}

		void weedFirstLayer() {

			int a,b;
			DistPairMap::iterator distIter;

			for( distIter = g_distPairMap.begin(); distIter != g_distPairMap.end() ; ++distIter) {
				a = LOWKEY((*distIter).second);
				b = HIGHKEY((*distIter).second);

				if((*distIter).first > 0.7) {
					if(pClustID(a) >0 && pClustID(b)> 0 && pClustID(a) == pClustID(b) && pClustID(a) != -4 && pClustID(b) != -4) {
						checkProbeBump(a);
						checkProbeBump(b);

					}
				} else {
					if(pClustID(a) != -4 && pClustID(b) != -4) {
						if(g_probes[a]->density == g_probes[b]->density) {
							if(pADist(a) < pADist(b)) {
								g_probes[b]->clusterId = -4;
								pIsSurvived(b) = -1;
							} else {
								g_probes[a]->clusterId = -4;
								pIsSurvived(a) = -1;
							}
						} else{
							if(g_probes[a]->density > g_probes[a]->density) {
								g_probes[b]->clusterId = -4;
								pIsSurvived(b) = -1;
							}
							else {
								g_probes[a]->clusterId = -4;
								pIsSurvived(a) = -1;
							}
						}
					}
				}
			}
		}

		/*
		   store distance (g_probes-g_probes) ( < 5A) for selection
		   insert probeIds(1.5A) into g_probes->nearProbes id for clustering;

		 */
		void assignProbePair() { 
			int cellId,newCellId,ori,count=0;
			Cell* temp;

			vector<int>::iterator it;

			int probePairKey;
			PairDistMap::iterator pairIter;

			clock_t begin,end;
			double dist;

			g_probePairMap.clear();
			for ( int i=0; i < g_probeCellList.size() ; i++) {
				for(int j=0; j < g_cells[g_probeCellList[i]]->probes.size(); j++) {
					g_probes[g_cells[g_probeCellList[i]]->probes[j]]->nearProbes.clear();
					for(int k=0; k< 275; k++) { // minimum distance < 5.66
						newCellId = g_probeCellList[i] + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);
						if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ)
							if(g_cells[newCellId]->probes.size() >0) { 

								for(int l=0; l < g_cells[newCellId]->probes.size(); l++) {
									dist = vsDist(pAtom(g_cells[g_probeCellList[i]]->probes[j]),pAtom(g_cells[newCellId]->probes[l]));
									if(g_cells[g_probeCellList[i]]->probes[j] != g_cells[newCellId]->probes[l] && dist > 0 && dist < 1.5) // for clustering (dist 1.5A and 0.7)
										g_probes[g_cells[g_probeCellList[i]]->probes[j]]->nearProbes.push_back(g_cells[newCellId]->probes[l]);
								}
							}
					}
				}
			}
		}


		void assignAtomProbe() {

			int cellId,newCellId,ori,count=0;
			Cell* temp;

			vector<int>::iterator it;

			double dist;
			for ( int i=0; i < g_atomCellList.size() ; i++) {
				for(int j=0; j < g_cells[g_atomCellList[i]]->atoms.size(); j++) {
					for(int k=0; k< 275; k++) { // minimun distance < 5.66
						newCellId = g_atomCellList[i] + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);
						if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ)
							if(g_cells[newCellId]->probes.size() >0) { 
								for(int l=0; l < g_cells[newCellId]->probes.size(); l++) {
									dist = vsDist(prAtom(g_cells[g_atomCellList[i]]->atoms[j]),pAtom(g_cells[newCellId]->probes[l]));
									if(dist <= 5.2){ // g_probes-protein 3.8A + g_probes-g_probes 1.4A
										g_proteinAtoms[g_cells[g_atomCellList[i]]->atoms[j]]->nearProbes.push_back(g_cells[newCellId]->probes[l]);
										//								cout << "OK" << endl;
									}
								}
							}
					}
				}
			}
		}

		/*
		   void assignProbePair_re() {
		   int cellId,newCellId,ori,count=0;
		   Cell* temp;

		   vector<int>::iterator it;

		   int probePairKey;
		   PairDistMap::iterator pairIter;

		   clock_t begin,end;
		   double dist;

		   g_probePairMap.clear();
		   g_distPairMap.clear();
		   for ( int i=0; i < g_probeCellList.size() ; i++) {
		   for(int j=0; j < g_cells[g_probeCellList[i]]->probes.size(); j++) {
		   g_probes[g_cells[g_probeCellList[i]]->probes[j]]->nearProbes.clear();
		   for(int k=0; k< 275; k++) {
		   newCellId = g_probeCellList[i] + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);
		   if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ) {
		   if(g_cells[newCellId]->probes.size() >0) { 
		   for(int l=0; l < g_cells[newCellId]->probes.size(); l++) {
		   probePairKey =getPairKey(g_cells[g_probeCellList[i]]->probes[j],g_cells[newCellId]->probes[l]);
		   pairIter = g_probePairMap.find(probePairKey);
		   if(g_probePairMap.end() == pairIter) {
		   dist = vsDist(pAtom(g_cells[g_probeCellList[i]]->probes[j]),pAtom(g_cells[newCellId]->probes[l]));
		   if(dist > 0 && dist <=2.8) {
		   g_probePairMap.insert(PairDistMap::value_type(probePairKey,dist));
		   g_distPairMap.insert(DistPairMap::value_type(dist,probePairKey));
		   }
		   } else {
		   dist = (*pairIter).second;
		   }
		   if( g_cells[g_probeCellList[i]]->probes[j]!= g_cells[newCellId]->probes[l] && dist > 0 && dist <= 2.8)  { 
		   g_probes[g_cells[g_probeCellList[i]]->probes[j]]->nearProbes.push_back(g_cells[newCellId]->probes[l]);
		   }
		   }
		   }
		   if(g_cells[newCellId]->atoms.size() >0) { 
		   for(int l=0; l < g_cells[newCellId]->atoms.size(); l++) {
		   dist = vsDist(pAtom(g_cells[g_probeCellList[i]]->probes[j]),prAtom(g_cells[newCellId]->atoms[l]));
		   if(dist <= 5.2)  // g_probes - new g_probes : 1.4A  + protein - new g_probes : 3.8A for generation
		   g_proteinAtoms[g_cells[newCellId]->atoms[l]]->nearProbes.push_back(g_cells[g_probeCellList[i]]->probes[j]);
		   }
		   }
		   }
		   }
		   }
		   }
		   }
		 */
		void assignProbePair_next(vector<int> probeIds) { // for 4>= g_layer probes
			int cellId,newCellId,ori,count=0;
			Cell* temp;

			vector<int>::iterator it;

			int probePairKey;
			PairDistMap::iterator pairIter;

			clock_t begin,end;
			double dist;

			for ( int i=0; i < probeIds.size() ; i++) {
				g_probes[probeIds[i]]->nearProbes.clear();
				cellId = getCellId(pAtom(probeIds[i]),2);
				for(int k=0; k< 117; k++) {
					newCellId = cellId + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);
					if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ) {
						if(g_cells[newCellId]->probes.size() >0) { 
							for(int l=0; l < g_cells[newCellId]->probes.size(); l++) {
								dist = vsDist(pAtom(probeIds[i]),pAtom(g_cells[newCellId]->probes[l]));
								if(dist <= 2.80)
									g_probes[probeIds[i]]->nearProbes.push_back(g_cells[newCellId]->probes[l]);
							}
						}
					}
				}
			}
		}
		void assignProbePair_sub(vector<int> probeList) {
			int cellId,newCellId,nProbeID,count=0;
			Cell* temp;

			vector<int>::iterator it;

			int probePairKey;
			PairDistMap::iterator pairIter;

			double dist;

			for(int i=0; i< probeList.size(); i++) {
				nProbeID = probeList[i];
				g_probes[nProbeID]->nearProbes.clear();
				cellId = getCellId(g_probes[nProbeID]->point,2);
				for(int k=0; k< 117; k++) { // minimun < 2.83
					newCellId = cellId + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);
					if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ) {
						if(g_cells[newCellId]->probes.size() >0) { 
							for(int l=0; l < g_cells[newCellId]->probes.size(); l++) {
								dist = vsDist(pAtom(nProbeID),pAtom(g_cells[newCellId]->probes[l]));
								if( nProbeID != g_cells[newCellId]->probes[l] && dist > 0 && dist <= 2.8)  {
									g_probes[nProbeID]->nearProbes.push_back(g_cells[newCellId]->probes[l]);
								}
							}
						}
					}
				}

			}
		}


		void assignProbePair_final() {
			int cellId,newCellId,ori;
			Cell* temp;

			vector<int>::iterator it;

			int probePairKey,count,all;
			PairDistMap::iterator pairIter;

			double dist,checkDist;

			g_probePairMap.clear();

			for ( int i=0; i < g_probeCellList.size() ; i++) {

				for(int j=0; j < g_cells[g_probeCellList[i]]->probes.size(); j++) {

					count=all=0;
					g_probes[g_cells[g_probeCellList[i]]->probes[j]]->nearProbes.clear();

					for(int k=0; k< 117; k++) { // minimun < 2.83
						newCellId = g_probeCellList[i] + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);
						if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ) {
							if(g_cells[newCellId]->probes.size() >0) { 
								for(int l=0; l < g_cells[newCellId]->probes.size(); l++) {
									dist = vsDist(pAtom(g_cells[g_probeCellList[i]]->probes[j]),pAtom(g_cells[newCellId]->probes[l]));
									if( g_cells[g_probeCellList[i]]->probes[j] != g_cells[newCellId]->probes[l] && dist > 0 && dist <= 2.8)  {
										g_probes[g_cells[g_probeCellList[i]]->probes[j]]->nearProbes.push_back(g_cells[newCellId]->probes[l]);
										if(dist < 1.41 && pNLayer(g_cells[newCellId]->probes[l]) == 1) 
											count++;

										if(dist < 1.41 && pNLayer(g_cells[newCellId]->probes[l]) > 1) // step 4
											all++;

									}
								}
							}
						}

						g_probes[g_cells[g_probeCellList[i]]->probes[j]]->density1 = all;
						g_probes[g_cells[g_probeCellList[i]]->probes[j]]->density2 = count;
					}
				}
			}
		}


		int bumpcheckAtoms(Probe* new_probe) { // bump checking for first g_layer g_probes with protein atom
			int cellId,newCellId;
			double dist,cutoff,totalDist=0,minDist=9999;
			int bc,range,countTotal=0,calTotal=0,averageNum=0;

			vector <int> tempAtoms,tempCell;
			cellId = getCellId(new_probe->point,2);
			for(int k=0; k< 275; k++) { // minimun distance < 5.66
				newCellId = cellId + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);
				if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ) {
					if(g_cells[newCellId]->atoms.size() >0) { 
						for(int l=0; l < g_cells[newCellId]->atoms.size(); l++) {
							dist = vsDist(new_probe->point,prAtom(g_cells[newCellId]->atoms[l]));

							if(g_proteinAtoms[g_cells[newCellId]->atoms[l]]->isPolar == 1)
								cutoff = 2.5;
							else if(g_proteinAtoms[g_cells[newCellId]->atoms[l]]->isPolar == 2)
								cutoff = 2.2;
							else
								cutoff = 3.3;

							if(dist < cutoff-0.0001){
								return FALSE;
							}

							if(dist <= 5.20) { // g_probes - new g_probes : 1.4A  + protein - new g_probes : 3.8A for generation
								new_probe->nearAtoms.push_back(g_cells[newCellId]->atoms[l]);
							}
							if(dist  < 3.8001 && g_proteinAtoms[g_cells[newCellId]->atoms[l]]->isPolar < 2) { // need thinking
								totalDist += dist;
								averageNum++;

							}
							minDist = MIN( minDist , dist );
						}
					}
				}
			}
			new_probe->averageDist = totalDist/averageNum;
			new_probe->averageNum = averageNum;
			new_probe->closestDist = minDist;

			return TRUE;
		}


		int bumpcheckAtoms2(Probe* new_probe) {
			int cellId,newCellId;
			double dist,cutoff,cutoff_probe;
			double totalDist=0,minDist=9999;
			double count_surface=0;

			int bc,range,countTotal=0,calTotal=0,averageNum=0;

			vector <int> tempAtoms,tempCell,tempProbes;
			cellId = getCellId(new_probe->point,2);

			range = g_gridProperties.size();
			//	cutoff_probe = (PROBE_RADIUS_PROBE+new_probe->radius)*0.95 ;


			//	if(new_probe->numLayer < 4)
			//		cutoff_probe *= 0.9;

			cutoff_probe = 1.00; //for 1,2,3 g_layer g_probes bump check*

			for(int k=0; k< range; k++) {
				newCellId = cellId + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);
				if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ) {

					if(g_cells[newCellId]->probes.size() >0 ) { 
						for(int l=0; l < g_cells[newCellId]->probes.size(); l++) {
							dist = vsDist(new_probe->point,pAtom(g_cells[newCellId]->probes[l]));
							if(dist < cutoff_probe-0.0001)
								return FALSE;
						}
					}

					if(g_cells[newCellId]->atoms.size() >0) { 
						for(int l=0; l < g_cells[newCellId]->atoms.size(); l++) {
							dist = vsDist(new_probe->point,prAtom(g_cells[newCellId]->atoms[l]));

							if(g_proteinAtoms[g_cells[newCellId]->atoms[l]]->isPolar == 1)
								cutoff = 2.5;
							else if(g_proteinAtoms[g_cells[newCellId]->atoms[l]]->isPolar == 2)
								cutoff = 2.2;
							else
								cutoff = 3.3;

							if(dist < cutoff-0.0001 ){
								return FALSE;
							}

							if(dist  < 3.8001 && g_proteinAtoms[g_cells[newCellId]->atoms[l]]->isPolar < 2) { // need thinking
								totalDist += dist;
								averageNum++;

							}
							minDist = MIN( minDist , dist );

							/*
							   if(dist <= 4.9001) {
							   new_probe->nearAtoms.push_back(g_cells[newCellId]->atoms[l]);
							   tempAtoms.push_back(g_cells[newCellId]->atoms[l]);
							   }
							 */

						}
					}

					if(g_gridProperties[k]->max <= 14) {
						countTotal += g_cells[newCellId]->atoms.size();
						if(g_gridProperties[k]->min<= 10) {
							if(g_cells[newCellId]->atoms.size() >0) {
								for(int l=0; l < g_cells[newCellId]->atoms.size(); l++) {
									dist = vsDist(new_probe->point,prAtom(g_cells[newCellId]->atoms[l]));
									if(dist <= 7) {
										//								count_surface +=  5 - (2-g_proteinAtoms[g_cells[newCellId]->atoms[l]]->proteinBc)*0.5;
										if(g_proteinAtoms[g_cells[newCellId]->atoms[l]]->proteinBc < g_bcCutoff) {
											count_surface += 1;
										}
									}
								}
							}
						}

					} else {
						if(g_cells[newCellId]->atoms.size() >0) { 
							for(int l=0; l < g_cells[newCellId]->atoms.size(); l++) {
								dist = vsDist(new_probe->point,prAtom(g_cells[newCellId]->atoms[l]));
								if(dist <= 14) {
									calTotal++;
								}
							}
						}
					}
				}
			}



			/*
			   if(g_layer == 2)
			   for(int i=0; i< tempAtoms.size(); i++) {
			   g_proteinAtoms[tempAtoms[i]]->nearProbes.push_back(g_probes.size());
			   }
			 */
			new_probe->averageDist = totalDist/averageNum;
			new_probe->averageNum = averageNum;
			new_probe->closestDist = minDist;


			new_probe->numBc1 = countTotal+calTotal;
			new_probe->numBc2 = count_surface;
			if(countTotal+calTotal < g_bcCutoff)
				return FALSE;
			g_probeBcList.push_back(new_probe->numBc1);
			return TRUE;
		}

		int bumpcheckProbe(Probe* new_probe,int i) {
			int cellId,newCellId;
			double dist,cutoff;
			int bc,range,countTotal=0,calTotal=0;

			vector <int> tempAtoms,tempCell,tempProbes;


			cutoff = 1.20;

			for(int j=0 ; j < g_probes[i]->nearProbes.size(); j++) {
				dist = vsDist(new_probe->point,pAtom(g_probes[i]->nearProbes[j]));
				if(dist < cutoff - 0.0001)
					return FALSE;
			}

			return TRUE;
		}



		int bumpcheckAtoms3(Probe* new_probe,int p1, int p2, int p3) {
			int cellId,newCellId;
			double dist,cutoff,cutoff_probe;
			int bc,range,countTotal=0,calTotal=0;
			double count_surface=0;

			vector <int> tempAtoms,tempCell,tempProbes;
			cellId = getCellId(new_probe->point,2);

			range = g_gridProperties.size();

			if( !bumpcheckProbe ( new_probe, p1) )
				return FALSE;
			if( !bumpcheckProbe ( new_probe, p2) )
				return FALSE;
			if( !bumpcheckProbe ( new_probe, p3) )
				return FALSE;

			cutoff_probe = 1.20;


			for(int k=0; k< range; k++) {
				newCellId = cellId + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);
				if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ) {

					if( g_gridProperties[k]->min <= 2.5  && g_cells[newCellId]->probes.size() >0 ) { 
						for(int l=0; l < g_cells[newCellId]->probes.size(); l++) {
							dist = vsDist(new_probe->point,pAtom(g_cells[newCellId]->probes[l]));
							if(dist < cutoff_probe-0.0001)
								return FALSE;
						}
					}

					if(g_cells[newCellId]->atoms.size() >0) { 
						for(int l=0; l < g_cells[newCellId]->atoms.size(); l++) {
							dist = vsDist(new_probe->point,prAtom(g_cells[newCellId]->atoms[l]));

							if(g_proteinAtoms[g_cells[newCellId]->atoms[l]]->isPolar == 1)
								cutoff = 2.5;
							else if(g_proteinAtoms[g_cells[newCellId]->atoms[l]]->isPolar == 2)
								cutoff = 2.2;
							else
								cutoff = 3.3;

							if(dist < cutoff-0.0001){
								return FALSE;
							}
						}
					}

					if(g_gridProperties[k]->max <= 14) {
						countTotal += g_cells[newCellId]->atoms.size();

						if(g_gridProperties[k]->min <= 10) {
							if(g_cells[newCellId]->atoms.size() >0) {
								for(int l=0; l < g_cells[newCellId]->atoms.size(); l++) {
									dist = vsDist(new_probe->point,prAtom(g_cells[newCellId]->atoms[l]));
									if(dist <= 7) {
										if(g_proteinAtoms[g_cells[newCellId]->atoms[l]]->proteinBc < g_bcCutoff)  
											count_surface += 3;
									}
								}
							}

						}


					} else {
						if(g_cells[newCellId]->atoms.size() >0) { 
							for(int l=0; l < g_cells[newCellId]->atoms.size(); l++) {
								dist = vsDist(new_probe->point,prAtom(g_cells[newCellId]->atoms[l]));
								if(dist <= 14) {
									calTotal++;
								}
							}
						}
					}
				}
			}
			new_probe->numBc1 = countTotal+calTotal;
			new_probe->numBc2 = count_surface;

			if(new_probe->numBc1 + new_probe->numBc2  < g_bcCutoff)
				return FALSE; 

			return TRUE;
		}


		int dbscanClustering(std::vector<int> probeIds, double eps, int minPts, int step) {
			int a,b,c,i,j,k,clusterId= step ;
			double dist;
			std::vector <int> g_clusters;

			//printf("REMARK DBSCAN Clustering : g_probes size = %5d, eps = %6.2f, minPts = %3d\t Cluster_id = %3d\n",probeIds.size(), eps, minPts,step);
			//Calcluate density

			for(a=0; a < probeIds.size(); a++) {
				i = probeIds[a];
				g_probes[i]->isVisited = 0;
				g_probes[i]->clusterId = 0;
				g_probes[i]->density =0 ;
				g_probes[i]->mergedAtoms.clear();
				for(b=0; b<g_probes[i]->nearProbes.size(); b++) {
					j = g_probes[i]->nearProbes[b];
					if(i != j && find(probeIds.begin(), probeIds.end(), j) != probeIds.end() ) {

						dist = checkProbeDist(i,j);

						if(dist != -1 && dist < eps ) {
							g_probes[i]->density++;
							g_probes[i]->mergedAtoms.push_back(j);
						}
						if(dist == -1) 
							cout << "DIST error : " << i << "," << j << endl;
					}
				} //loop end


			}


			for(a=0; a < probeIds.size(); a++) {
				i = probeIds[a];
				if(g_probes[i]->isVisited == 0) {
					g_probes[i]->isVisited = 1;

					if(g_probes[i]->density < minPts) {
					} else {
						clusterId++;
						g_probes[i]->clusterId = clusterId;
						for(b=0; b < g_probes[i]->mergedAtoms.size(); b++) {
							k = g_probes[i]->mergedAtoms[b];
							if(g_probes[k]->isVisited == 0) {
								g_probes[k]->isVisited = 1;
								if(g_probes[k]->density >= minPts) {
									for(int l=0; l < g_probes[k]->mergedAtoms.size(); l++) {
										g_probes[i]->mergedAtoms.push_back(g_probes[k]->mergedAtoms[l]);
									}
								}
							}
							if(g_probes[k]->clusterId == 0 ) {
								g_probes[k]->clusterId = clusterId;
							}
						}
					}
				}
			}
			return clusterId;
		}

		int dbscanClusteringForSurface(std::vector<int> probeIds, double eps, int minPts) {
			int a,b,c,i,j,k,clusterId= 0 ;
			double dist;
			std::vector <int> g_clusters;

			for(a=0; a < probeIds.size(); a++) {
				i = probeIds[a];
				g_probes[i]->isVisited = 0;
				g_probes[i]->clusterIdFirst = 0;
				g_probes[i]->density =0 ;
				g_probes[i]->mergedAtoms.clear();
				for(b=0; b<probeIds.size(); b++) {
					j = probeIds[b];
					if(i != j ) {
						dist = checkProbeDist(i,j);
						if(dist != -1 && dist < eps ) {
							g_probes[i]->density++;
							g_probes[i]->mergedAtoms.push_back(j);
						}
						if(dist == -1) 
							cout << "DIST error : " << i << "," << j << endl;
					}
				} //loop end

			}

			for(a=0; a < probeIds.size(); a++) {
				i = probeIds[a];
				if(g_probes[i]->isVisited == 0) {
					g_probes[i]->isVisited = 1;

					if(g_probes[i]->density < minPts) {
					} else {
						clusterId++;
						g_probes[i]->clusterIdFirst = clusterId;
						for(b=0; b < g_probes[i]->mergedAtoms.size(); b++) {
							k = g_probes[i]->mergedAtoms[b];
							if(g_probes[k]->isVisited == 0) {
								g_probes[k]->isVisited = 1;
								if(g_probes[k]->density >= minPts) {
									for(int l=0; l < g_probes[k]->mergedAtoms.size(); l++) {
										g_probes[i]->mergedAtoms.push_back(g_probes[k]->mergedAtoms[l]);
									}
								}
							}
							if(g_probes[k]->clusterIdFirst == 0 ) {
								g_probes[k]->clusterIdFirst = clusterId;
							}
						}
					}
				}
			}
			return clusterId;
		}



		void execDbscanClusteringFirst() {
			vector<int> probeIds,firstProbeIds;
			int clusterId = 0,probe_size=0,i,j;
			int cellId,newCellId,countTotal,calTotal;
			double dist;
			int start, end;
			double upper, lower,upper_ratio,lower_ratio,alpha,beta,g_bcCutoffRatio;

			int probePairKey;

			int count_surface; 
			ClusterMap::iterator clusterIter;
			vector<int>::iterator it;


			for( i=0; i<g_probes.size(); i++){
				firstProbeIds.push_back(i);
			}

			clusterId = dbscanClustering( firstProbeIds,0.7001,1,0);

			assignClusterSet();

			probeIds.clear();

			for(clusterIter = g_clusterMap.begin(); clusterIter != g_clusterMap.end() ; clusterIter++){
				probeIds.insert(probeIds.end(),clusterIter->second->probeIds.begin(),clusterIter->second->probeIds.end());

				for(i=0 ; i < probeIds.size(); i++) {
					for(j=i+1; j < probeIds.size(); j++) {
						if(i!=j) {
							probePairKey = getPairKey(probeIds[i],probeIds[j]);
							dist = vsDist(pAtom(probeIds[i]),pAtom(probeIds[j]));
							g_distPairMap.insert(DistPairMap::value_type(dist,probePairKey));
						}

					}

				}

				weedFirstLayer();
				g_distPairMap.clear();
				probeIds.clear();
			}


			for ( int i=0; i < g_probeCellList.size() ; i++) {
				g_cells[g_probeCellList[i]]->probes.clear(); {
				}
			}

			i=0;
			sort(g_probes.begin(),g_probes.end(),compareProbeSurviving);
			g_probeCellList.clear();
			//Selection surviving g_probes and calculation BC
			while(g_probes[i]->isSurvived==1){
				cellId = getCellId(pAtom(i),2);
				g_probes[i]->numBc1 = g_probes[i]->numBc2 =0;
				countTotal = 0;
				calTotal = 0;

				count_surface=0;
				for(int k=0; k< g_gridProperties.size(); k++) {
					newCellId = cellId + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);
					if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ) {
						if(g_gridProperties[k]->max <= 14) {
							countTotal += g_cells[newCellId]->atoms.size();
						} else {
							if(g_cells[newCellId]->atoms.size() >0) { 
								for(int l=0; l < g_cells[newCellId]->atoms.size(); l++) {
									dist = vsDist(pAtom(i),prAtom(g_cells[newCellId]->atoms[l]));
									if(dist <= 14) {
										calTotal++;
									}
								}
							}
						}
					}
				}

				g_probes[i]->numBc1 = countTotal+calTotal;
				g_probeBcList.push_back(g_probes[i]->numBc1);
				i++;

			}


			//	g_probes.erase(g_probes.begin()+i,g_probes.end());

			sort(g_probeBcList.begin(),g_probeBcList.end(),std::greater<int>());

			double tempSum=0;
			for(i=0; i<10; i++) {
				tempSum += g_probeBcList[i];
			}

			g_top10Mean = tempSum/10.0;

			g_quartile = g_probeBcList[g_probeBcList.size()*0.25];


			if(g_proteinAtoms.size() <= 999) {
				upper = 637 + (-895.18674*exp(-0.00231*g_proteinAtoms.size()));
				lower = 314;
				upper_ratio = 0.5000;
				lower_ratio = 0.8000;
			}
			if(g_proteinAtoms.size() >999) {
				upper = 637 + (-895.18674*exp(-0.00231*g_proteinAtoms.size()));
				lower = 132.86*log(g_proteinAtoms.size())-601;
				upper_ratio = 0.5000;
				lower_ratio = 0.8000;
			}
			alpha = (lower_ratio-upper_ratio)/(lower-upper);
			beta = upper_ratio - alpha*upper;

			g_bcCutoffRatio = alpha*g_top10Mean + beta;

			if(g_bcCutoffRatio < upper_ratio) 
				g_bcCutoffRatio = upper_ratio;
			if(g_bcCutoffRatio > lower_ratio) 
				g_bcCutoffRatio = lower_ratio;


			//	g_bcCutoffRatio = 0.70;

			if( g_bcCutoffRatioG > 0) {
				g_bcCutoffRatio = g_bcCutoffRatioG;
			}
			g_bcCutoff = g_top10Mean * g_bcCutoffRatio;



			printf("REMARK Total 1 g_layer select g_probes number : %5d ,protein atom : %5d, upper : %5.2f(%5.2f), lower %5.2f(%5.2f)\nREMARK top 10 mean : %5.2f,  g_quartile : %4.1f ,  bc cutoff ratio : %5.4f  , g_bcCutoff , %5.2f\n",g_probes.size(),g_proteinAtoms.size(),upper,upper_ratio,lower,lower_ratio,g_top10Mean,g_quartile , g_bcCutoffRatio, g_bcCutoff);

			sort(g_probes.begin(),g_probes.end(),compareProbeBc);
			g_probeCellList.clear();

			Probes temp_probe;

			for(i=0 ; i<g_probes.size(); i++) {

				cellId = getCellId(pAtom(i),2);
				count_surface=0;
				for(int k=0; k< 667 ; k++) {
					newCellId = cellId + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);
					if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ) {
						if(g_cells[newCellId]->atoms.size() >0) { 
							for(int l=0; l < g_cells[newCellId]->atoms.size(); l++) {
								dist = vsDist(pAtom(i),prAtom(g_cells[newCellId]->atoms[l]));
								if(dist <= 7) {
									count_surface += 1;
								}
							}
						}
					}
				}

				g_probes[i]->numBc2 = count_surface;


				if(g_probes[i]->numBc1 + g_probes[i]->numBc2 >= g_bcCutoff){
					temp_probe.push_back(g_probes[i]);
				}
			}

			g_probes.clear();
			g_probes.assign(temp_probe.begin(),temp_probe.end());


			for(i=0 ; i< g_probes.size(); i++) {
				cellId = getCellId(pAtom(i),2);
				g_cells[cellId]->probes.push_back(i);
				g_probeCellList.push_back(cellId);
			}



			printf("REMARK Total 1 g_layer bc cutoff g_probes number : %5d\n",g_probes.size());

			sort(g_probeCellList.begin(),g_probeCellList.end());
			it = unique(g_probeCellList.begin(),g_probeCellList.end());
			g_probeCellList.resize(it - g_probeCellList.begin());

		}

		int execDbscanClustering(double eps, int minPts ,int start) {
			vector<int> probeIds,firstProbeIds;
			int clusterId = start,probe_number,a,b,i,j,min_cid;
			int append_num;

			double dist,minDist;
			ClusterMap::iterator clusterIter;
			probeIds.clear();

			for(i=0; i<g_probes.size(); i++){
				pClustID(i)=0; 
				g_probes[i]->clusterIdFirst = -1;
				g_probes[i]->clusterId1 = -1;
				g_probes[i]->clusterId2 = -1;
				g_probes[i]->clusterId3 = -1;
				g_probes[i]->density = 0;
				probeIds.push_back(i);
			}

			printf("REMARK Final DBSCAN Clustering : g_probes size = %5d, eps = %6.2f, minPts = %3d\n",probeIds.size(), eps, minPts);
			clusterId = dbscanClustering( probeIds,eps,minPts,0);

			for(int i=0; i<g_probes.size(); i++){
				g_probes[i]->clusterId1 = pClustID(i);
				g_probes[i]->clusterId=-1;
			}

			return clusterId;

		}

		int execDbscanClusteringSub(double eps, int minPts ,int start) {
			vector<int> probeIds,probe_ids_sub,contact_probe,list;
			int clusterId = start,probe_number,a,b,c,i,j,k,near_probes,cid,cid_sub;
			int append_num,probePairKey,count,count_all,flag;
			ClusterMap::iterator clusterIter,clusterSubIter;
			double dist,minDist,dist_sub,ratio,count_contact;
			int targetCid, targetCid1;
			PairDistMap::iterator pairIter;
			map<int, double> contact;
			map<int, vector <int> > contact_target; 
			map<int, double>::iterator it;
			map<int, vector <int> >::iterator it_target;


			printf("REMARK Final Sub DBSCAN Clustering : g_probes size = %5d, eps = %6.2f, minPts = %3d\n",g_probes.size(), eps, minPts);

			for(clusterIter = g_clusterMap.begin(); clusterIter != g_clusterMap.end() ; clusterIter++){
				probeIds.clear();
				probeIds.insert(probeIds.end(),clusterIter->second->probeIds.begin(),clusterIter->second->probeIds.end());
				clusterId = dbscanClustering(probeIds, eps,minPts,0);

				insertClusterSetSub(probeIds); // g_clusters result 2.5, 18


				g_clusters.clear();

				probeIds.clear();
				count =0;
				for(clusterSubIter = g_subClusterMap.begin(); clusterSubIter != g_subClusterMap.end(); clusterSubIter++) {
					ratio = (double) clusterSubIter->second->type / (double) clusterSubIter->second->probeIds.size();
					if( clusterSubIter->second->type < 2 || clusterSubIter->second->probeIds.size() < 51 || clusterSubIter->first == 0) {
						//					cout << "REMARK Cluster List : " << clusterIter->first << " - " <<clusterSubIter->first << " : " << clusterSubIter->second->probeIds.size() << " : " << clusterSubIter->second->type << " * " << endl;
						probeIds.insert(probeIds.end(),clusterSubIter->second->probeIds.begin(),clusterSubIter->second->probeIds.end());
					} else {
						//				cout << "REMARK Cluster List : " << clusterIter->first << " - " << clusterSubIter->first << " : " << clusterSubIter->second->probeIds.size() << " : " << clusterSubIter->second->type << endl;
						clusterSubIter->second->flag = 1;
						count++;
					}


				}

				//		cout << "REMARK Noise Probes : " << probeIds.size() << endl;


				printf("REMARK Final Sub 2 DBSCAN Clustering : g_probes size = %5d, eps = %6.2f, minPts = %3d\n",probeIds.size(), 1.8, 3);
				if(count == 0) {
					clusterId = 0;
				}

				dbscanClustering(probeIds,1.8,3,clusterId);

				insertClusterSetSub(probeIds); // noise g_probes g_clusters 1.8,3
				for(clusterSubIter = g_subClusterMap.begin(); clusterSubIter != g_subClusterMap.end(); clusterSubIter++) {
					//			cout << "REMARK Cluster Sub 2 List : " << clusterIter->first << " - " << clusterSubIter->first << " : " << clusterSubIter->second->probeIds.size() << " : " << clusterSubIter->second->type << endl;
				}
				// Noise Probe to Cluster
				clusterSubIter = g_subClusterMap.find(0);
				if( clusterSubIter != g_subClusterMap.end()) {
					//			cout << "REMARK Sub 2 g_clusters probes: " << clusterSubIter->second->clusterId << " : " <<clusterSubIter->second->probeIds.size() << endl;
					for(j=0 ; j < clusterSubIter->second->probeIds.size() ; j++) {
						a = clusterSubIter->second->probeIds[j];
						//						cout << "REMARK Noise : " << a << " : "  << g_probes[a]->nearProbes.size() << " : " <<  g_probes[a]->clusterId <<" -> " ;
						minDist = 9999;
						for ( k = 0; k < g_probes[a]->nearProbes.size(); k++) {
							b = g_probes[a]->nearProbes[k];
							if( g_probes[b]->clusterId > 0 &&  g_probes[a]->clusterId != g_probes[b]->clusterId) {
								dist = vsDist(pAtom(a),pAtom(b));
								//														cout << "  " << dist << " , " ;
								if(minDist > dist) {
									minDist = dist;
									g_probes[a]->clusterId = g_probes[b]->clusterId;
								}
							}
						}
						//				cout << g_probes[a]->clusterId << " : " << minDist << endl;
						if(g_subClusterMap.find(g_probes[a]->clusterId) != g_subClusterMap.end())
							g_subClusterMap.find(g_probes[a]->clusterId)->second->probeIds.push_back(a);
					}
				}

				insertClusterSetSub(clusterIter->second->probeIds); // noise g_probes g_clusters 1.8,3



				for(clusterSubIter = g_subClusterMap.begin(); clusterSubIter != g_subClusterMap.end(); clusterSubIter++) {

					//			if(clusterSubIter->second->probeIds.size() < 51)
					//				cout << "REMARK Cluster List : " << clusterIter->first << " - " << clusterSubIter->first << " : " << clusterSubIter->second->probeIds.size() << " : " << clusterSubIter->second->type << " * " << endl;
					//			else 
					//				cout << "REMARK Cluster List : " << clusterIter->first << " - " << clusterSubIter->first << " : " << clusterSubIter->second->probeIds.size() << " : " << clusterSubIter->second->type << endl;

					if(clusterSubIter->first > clusterId) {

						g_clusters.push_back(clusterSubIter->second);
					}
				}


				for(i=0; i<g_clusters.size(); i++) {
					//			cout << "REMARK Sub g_clusters : " << g_clusters[i]->clusterId << " : " << g_clusters[i]->probeIds.size() << endl;
					count=0;
					contact.clear();
					contact_target.clear();
					for(j=0 ; j < g_clusters[i]->probeIds.size() ; j++) {
						a = g_clusters[i]->probeIds[j];
						//				cout << "REMARK Noise Cluster: " << a << " :: "  << g_probes[a]->nearProbes.size() << " : " <<  g_probes[a]->clusterId <<" -> " ;
						minDist = 9999;
						targetCid = -1;
						for ( k = 0; k < g_probes[a]->nearProbes.size(); k++) {
							b = g_probes[a]->nearProbes[k];
							if( g_probes[a]->clusterId != g_probes[b]->clusterId && g_probes[b]->clusterId <= clusterId &&  g_probes[b]->clusterId >0) {
								dist = vsDist(pAtom(a),pAtom(b));
								if(minDist > dist) {
									minDist = dist;
									targetCid = g_probes[b]->clusterId;
									//							cout << "  " << dist << "(" << targetCid  <<") , ";
								}
								if(dist <=1.9) {
									if(contact_target.find(g_probes[b]->clusterId) == contact_target.end()) {
										list.clear();
										list.push_back(b);
										contact_target.insert( map < int, vector <int> >::value_type(g_probes[b]->clusterId,list));
									}
									else {
										contact_target.find(g_probes[b]->clusterId)->second.push_back(b);
									}
								}
							}
						}

						if(minDist <= 1.9) {
							if(contact.find(targetCid) == contact.end()) {
								contact[targetCid] = 10 + (1.9-minDist);
							} else {
								dist_sub = contact[targetCid]- ((int)(contact[targetCid]/10) * 10) + 1.9;
								dist_sub = MIN(dist_sub,minDist);
								contact[targetCid] = ((int)(contact[targetCid]/10) * 10) + 10 + (1.9 - dist_sub);
							}
						}
						//				cout << "->" << targetCid  << " : " << minDist << endl;
					}

					for( it_target = contact_target.begin() ; it_target != contact_target.end(); it_target++) {
						//				cout << "REMARK Contact target: "  << it_target->first << " : ";

						sort(it_target->second.begin(),it_target->second.end());
						it_target->second.resize(unique(it_target->second.begin(),it_target->second.end()) - it_target->second.begin());

						for(int id=0; id < it_target->second.size() ;id++) {
							//					cout << it_target->second[id] << ",";
						}
						//			cout << endl;
					}

					targetCid = -1;
					count_contact = -1;
					for( it = contact.begin() ; it != contact.end(); it++) {
						//				cout << "REMARK Contact : "  << it->first << " : " << it->second << endl;

						if(g_subClusterMap.find(it->first) != g_subClusterMap.end()){


							int count_cutoff = 20;

							if(g_clusters[i]->probeIds.size() > 10) 
								count_cutoff = 30;

							if(g_clusters[i]->probeIds.size() > 20) 
								count_cutoff = 40;

							if(g_clusters[i]->probeIds.size() > 30) 
								count_cutoff = 50;

							if(g_clusters[i]->probeIds.size() > 40) 
								count_cutoff = 60;





							if(it->second > count_contact && ( it->second + contact_target.find(it->first)->second.size()*10) > count_cutoff) {
								//						if(contact_target.find(it->first)->second.size() >= 4) {
								targetCid = it->first;
								count_contact = it->second;
								//						}
							}

						}
					}

					if(targetCid > -1) {
						//				cout << "REMARK combine : "  << targetCid << " : " << count_contact << endl;
						for(j=0 ; j < g_clusters[i]->probeIds.size() ; j++) {
							//					cout << "REMARK " << g_clusters[i]->probeIds[j] << " : " << g_probes[g_clusters[i]->probeIds[j]]->clusterId << " ->" << targetCid << endl;
							g_probes[g_clusters[i]->probeIds[j]]->clusterId = targetCid;
							g_subClusterMap.find(targetCid)->second->probeIds.push_back(g_clusters[i]->probeIds[j]);
						}
						g_subClusterMap.find(g_clusters[i]->clusterId)->second->flag = 0;
					}



				}


			}

		}

		int GetTripleVertex(Vector3 Ri, Vector3 Rj, Vector3 Rk, Vector3* Rp1, Vector3* Rp2, double sgmi, double sgmj, double sgmk, double sgmp) {

			Vector3  Tij, Tjk, Tik, Rb, Ru, Xu, Yu, Zu, Temp, Temp1, Temp2;
			double  h, temp;

			Xu=normalizeVector(vectorDiff(Rj, Ri));
			Zu=normalizeVector(vectorCross(Xu, vectorDiff(Rk, Ri)));
			Yu=vectorCross(Zu, Xu);

			temp=(SQUARE(sgmi+sgmp) - SQUARE(sgmj+sgmp)) / (2 * vsDistSq(Rj, Ri));
			Temp1=vectorScale(vectorSum(Ri, Rj), 0.5);
			Temp2=vectorScale(vectorDiff(Rj,Ri), temp);
			Tij=vectorSum(Temp1, Temp2);

			temp=(SQUARE(sgmj+sgmp) - SQUARE(sgmk+sgmp)) / (2 * vsDistSq(Rk, Rj));
			Temp1=vectorScale(vectorSum(Rj, Rk), 0.5);
			Temp2=vectorScale(vectorDiff(Rk,Rj), temp);
			Tjk=vectorSum(Temp1, Temp2);

			temp=(SQUARE(sgmi+sgmp) - SQUARE(sgmk+sgmp)) / (2 * vsDistSq(Rk, Ri));
			Temp1=vectorScale(vectorSum(Ri, Rk), 0.5);
			Temp2=vectorScale(vectorDiff(Rk,Ri), temp);
			Tik=vectorSum(Temp1, Temp2);

			temp=vectorDot(vectorDiff(Tik,Tij),vectorDiff(Tik,Ri))/ vectorDot(vectorDiff(Tik,Ri), Yu);

			Ru=vectorScale(Yu, temp);

			Rb=vectorSum(Tij, Ru);
			temp=SQUARE(sgmi+sgmp)-vsDistSq(Rb,Ri);
			if(temp<0)
				return FALSE;
			h=sqrt(temp);
			(*Rp1)=vectorSum(Rb, vectorScale(Zu, h));
			(*Rp2)=vectorDiff(Rb, vectorScale(Zu, h));

			if(fabs(vsDist(Ri, (*Rp1))-sgmi-sgmp)>0.01)
				return -1;
			if(fabs(vsDist(Rj, (*Rp1))-sgmj-sgmp)>0.01)
				return -2;
			if(fabs(vsDist(Rk, (*Rp1))-sgmk-sgmp)>0.01)
				return -3;


			//if(g_layer==2)
			//	printf("GETRIPLE : %8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f\n",Ri.x,Ri.y,Ri.z,Rj.x,Rj.y,Rj.z,Rk.x,Rk.y,Rk.z,sgmi,sgmj,sgmk,sgmp);
			return TRUE;

		}

		void extractProbeProperty(const Probe* probe, int index, ProbePropertySnapshot& snapshot) {
			snapshot.index = index;
			snapshot.serialNumber = probe->serialNumber;
			snapshot.type = probe->type;
			snapshot.x = probe->point.x;
			snapshot.y = probe->point.y;
			snapshot.z = probe->point.z;
			snapshot.radius = probe->radius;
			snapshot.charge = probe->charge;
			snapshot.isPolar = probe->isPolar;
			snapshot.isSurvived = probe->isSurvived;
			snapshot.clusterId = probe->clusterId;
			snapshot.clusterId1 = probe->clusterId1;
			snapshot.clusterId2 = probe->clusterId2;
			snapshot.clusterId3 = probe->clusterId3;
			snapshot.numLayer = probe->numLayer;
			snapshot.density = probe->density;
			snapshot.density1 = probe->density1;
			snapshot.density2 = probe->density2;
			snapshot.step = probe->step;
			snapshot.numBc1 = probe->numBc1;
			snapshot.numBc2 = probe->numBc2;
			snapshot.numBc3 = probe->numBc3;
			snapshot.numBc4 = probe->numBc4;
			snapshot.numBc5 = probe->numBc5;
			snapshot.numBcProbe = probe->numBcProbe;
			snapshot.numBcProbeTotal = probe->numBcProbeTotal;
			snapshot.bcPercent = probe->bcPercent;
			snapshot.numWtBcNum = probe->numWtBcNum;
			snapshot.closestDist = probe->closestDist;
			snapshot.closestProbeDist = probe->closestProbeDist;
			snapshot.averageDist = probe->averageDist;
			snapshot.averageNum = probe->averageNum;
			snapshot.closestAtom = probe->closestAtom;
			snapshot.contactAtoms = probe->contactAtoms;
			snapshot.nearProbes = probe->nearProbes;
			snapshot.nearAtoms = probe->nearAtoms;
		}

		bool saveProbePropertiesToJson(const string& outputPath) {
			ofstream ofs(outputPath.c_str());
			if (!ofs) {
				cerr << "Pass: failed to write probe properties json '" << outputPath << "': " << strerror(errno) << endl;
				return false;
			}

			ofs << "{\"probeCount\":" << g_probes.size() << ",\"probes\":[";
			ProbePropertySnapshot snapshot;
			for (size_t i = 0; i < g_probes.size(); ++i) {
				extractProbeProperty(g_probes[i], static_cast<int>(i), snapshot);
				if (i > 0) {
					ofs << ",";
				}
				ofs << "{";
				ofs << "\"index\":" << snapshot.index;
				ofs << ",\"serialNumber\":" << snapshot.serialNumber;
				ofs << ",\"type\":\"" << escapeJsonString(snapshot.type) << "\"";
				ofs << ",\"point\":{\"x\":" << snapshot.x << ",\"y\":" << snapshot.y << ",\"z\":" << snapshot.z << "}";
				ofs << ",\"radius\":" << snapshot.radius;
				ofs << ",\"charge\":" << snapshot.charge;
				ofs << ",\"isPolar\":" << snapshot.isPolar;
				ofs << ",\"isSurvived\":" << snapshot.isSurvived;
				ofs << ",\"clusterId\":" << snapshot.clusterId;
				ofs << ",\"clusterId1\":" << snapshot.clusterId1;
				ofs << ",\"clusterId2\":" << snapshot.clusterId2;
				ofs << ",\"clusterId3\":" << snapshot.clusterId3;
				ofs << ",\"numLayer\":" << snapshot.numLayer;
				ofs << ",\"density\":" << snapshot.density;
				ofs << ",\"density1\":" << snapshot.density1;
				ofs << ",\"density2\":" << snapshot.density2;
				ofs << ",\"step\":" << snapshot.step;
				ofs << ",\"numBc1\":" << snapshot.numBc1;
				ofs << ",\"numBc2\":" << snapshot.numBc2;
				ofs << ",\"numBc3\":" << snapshot.numBc3;
				ofs << ",\"numBc4\":" << snapshot.numBc4;
				ofs << ",\"numBc5\":" << snapshot.numBc5;
				ofs << ",\"numBcProbe\":" << snapshot.numBcProbe;
				ofs << ",\"numBcProbeTotal\":" << snapshot.numBcProbeTotal;
				ofs << ",\"bcPercent\":" << snapshot.bcPercent;
				ofs << ",\"numWtBcNum\":" << snapshot.numWtBcNum;
				ofs << ",\"closestDist\":" << snapshot.closestDist;
				ofs << ",\"closestProbeDist\":" << snapshot.closestProbeDist;
				ofs << ",\"averageDist\":" << snapshot.averageDist;
				ofs << ",\"averageNum\":" << snapshot.averageNum;
				ofs << ",\"closestAtom\":" << snapshot.closestAtom;
				ofs << ",\"contactAtoms\":";
				writeIntArrayJson(ofs, snapshot.contactAtoms);
				ofs << ",\"nearProbes\":";
				writeIntArrayJson(ofs, snapshot.nearProbes);
				ofs << ",\"nearAtoms\":";
				writeIntArrayJson(ofs, snapshot.nearAtoms);
				ofs << "}";
			}
			ofs << "]}";
			if (!ofs) {
				cerr << "Pass: failed while writing probe properties json '" << outputPath << "'." << endl;
				return false;
			}
			return true;
		}

		void displayProbeAtom() {
			string Probe,contact_name_1,contact_name_2,contact_name_3,contact_1,contact_2, contact_3,resn;
			int contact_resn_1,contact_resn_2,contact_resn_3;
			int idx=1,ter_idx=0,clusterId;
			string chain; 

			ResidueNameMap::iterator residue_name_iter;
			//	sort(g_probes.begin(),g_probes.end(),compareProbeClusterId);
			for(int i=0 ; i< g_probes.size() ; i++) {
				if(g_probes[i]->numLayer == 1) {
					contact_name_1 = prType(pCAtoms(i,0));
					contact_name_2 = prType(pCAtoms(i,1));
					contact_name_3 = prType(pCAtoms(i,2));

					residue_name_iter = g_residueNames.find(ResiName(pCAtoms(i,0)));
					contact_1 = residue_name_iter->second;
					residue_name_iter = g_residueNames.find(ResiName(pCAtoms(i,1)));
					contact_2 = residue_name_iter->second;
					residue_name_iter = g_residueNames.find(ResiName(pCAtoms(i,2)));
					contact_3 = residue_name_iter->second;

					contact_resn_1 = ResiIdx(pCAtoms(i,0))+1;
					contact_resn_2 = ResiIdx(pCAtoms(i,1))+1;
					contact_resn_3 = ResiIdx(pCAtoms(i,2))+1;


				}

				if(g_probes[i]->numLayer == 2) {
					contact_name_1 = prType(pCAtoms(i,0));
					contact_name_2 = prType(pCAtoms(i,1));
					contact_name_3 = "PRB";
					residue_name_iter = g_residueNames.find(ResiName(pCAtoms(i,0)));
					contact_1 = residue_name_iter->second;
					residue_name_iter = g_residueNames.find(ResiName(pCAtoms(i,1)));
					contact_2 = residue_name_iter->second;
					contact_3 = "p"; //pCAtoms(i,2);
					contact_resn_1 = ResiIdx(pCAtoms(i,0))+1;
					contact_resn_2 = ResiIdx(pCAtoms(i,1))+1;
					contact_resn_3 = pCAtoms(i,2)+1;


				} 
				if(g_probes[i]->numLayer == 3) {
					contact_name_1 = prType(pCAtoms(i,0));
					contact_name_2 = "PRB";
					contact_name_3 = "PRB";
					residue_name_iter = g_residueNames.find(ResiName(pCAtoms(i,0)));
					contact_1 = residue_name_iter->second;

					contact_2 = "p"; //pCAtoms(i,1);
					contact_3 = "p"; //pCAtoms(i,2);
					contact_resn_1 = ResiIdx(pCAtoms(i,0))+1;
					contact_resn_2 = pCAtoms(i,1)+1;
					contact_resn_3 = pCAtoms(i,2)+1;

				} 

				if(g_probes[i]->numLayer >3) {
					contact_name_1 = "PRB";
					contact_name_2 = "PRB";
					contact_name_3 = "PRB";
					contact_1 = "p"; //pCAtoms(i,0);
					contact_2 = "p"; //pCAtoms(i,1);
					contact_3 = "p"; //pCAtoms(i,2);
					contact_resn_1 = 1;//pCAtoms(i,0)+1;
					contact_resn_2 = 1;//pCAtoms(i,1)+1;
					contact_resn_3 = 1;//pCAtoms(i,2)+1;
				} 
				if(g_probes[i]->isSurvived > 0  ) {
					if(g_probes[i]->isPolar == 0) {
						pType(i) = "C";
					} else if (g_probes[i]->isPolar ==1 ) {
						pType(i) = "N";
					} else {
						pType(i) = "O";
					}
				}

				if(g_probes[i]->numLayer < 4 && g_probes[i]->isPolar == 0 && g_probes[i]->closestDist < 3.2) {
					pType(i) = "F";
				}
				if(g_probes[i]->numLayer >3)
					pType(i) = "C";

				resn = "PRB A";

				//g_bcCutoff = g_minBcNum + (g_maxBcNum-g_minBcNum)*(1-ratio);



				/*if(g_probes[i]->clusterId1 < 9000) {
				  if(g_probes[i]->clusterId3 == -1)
				  clusterId = g_probes[i]->clusterId1*100+g_probes[i]->clusterId2;
				  else
				  clusterId = g_probes[i]->clusterId1*100+g_probes[i]->clusterId2*10+g_probes[i]->clusterId3;
				  } else*/
				if(g_probes[i]->clusterId <0 ) {
					resn = "PRB Z";

					clusterId = g_probes[i]->clusterId1;
				} else {
					if(g_probes[i]->clusterId >0) {
						clusterId = g_probes[i]->clusterId1*100+ g_probes[i]->clusterId;
					} else {
						clusterId = g_probes[i]->clusterId1*100;
					}
				}
				//			clusterId = g_probes[i]->clusterId1;


				//			clusterId = g_probes[i]->clusterId1;



				if(g_probes[i]->numBc3 != -1) {
					g_probes[i]->numBc3 = 0;
				}
				printf("ATOM  %5d  %-3s %s%4d    %8.3f%8.3f%8.3f      %d : %5d %5d : %3d    %4d :: %5.1f %5.1f",
						idx,
						pType(i).c_str(),
						resn.c_str(),
						//						clusterId,
						clusterId,
						pX(i),
						pY(i),
						pZ(i),
						g_probes[i]->numBc1,
						g_probes[i]->clusterId1,
						g_probes[i]->clusterId,
						g_probes[i]->numLayer,
						g_probes[i]->nearProbes.size(),
						g_probes[i]->averageDist,
						g_probes[i]->closestDist
					  );


				/*
				   printf("ATOM  %5d  %-3s %s%4d    %8.3f%8.3f%8.3f      %4d                %4.1f %s%-4d(%3s) %s%-4d(%3s) %s%-4d(%3s) : %2d  :%3d %3d %3d : %3d %3d: %3d :: %5.2f %5.2f%2d %2d %3d",
				   idx,
				   pType(i).c_str(),
				   resn.c_str(),
				//						clusterId,
				clusterId,
				pX(i),
				pY(i),
				pZ(i),
				g_probes[i]->density,
				//g_probes[i]->numBc1, //g_probes[i]->numLayer,//g_probes[i]->numBc1,
				g_probes[i]->numBc2,
				//	g_probes[i]->numLayer,
				//						g_probes[i]->numBc2 + (g_probes[i]->numBcProbe/2),
				contact_1.c_str(),
				contact_resn_1,
				contact_name_1.c_str(),
				contact_2.c_str(),
				contact_resn_2,
				contact_name_2.c_str(),
				contact_3.c_str(),
				contact_resn_3,
				contact_name_3.c_str(),
				g_probes[i]->numLayer,
				g_probes[i]->nearProbes.size(),
				g_probes[i]->nearAtoms.size(),
				g_probes[i]->density,
				g_probes[i]->clusterId1,
				g_probes[i]->clusterId,
				g_probes[i]->isSurvived,
				g_probes[i]->averageDist,
				g_probes[i]->closestDist,
				pANum(i),
				g_probes[i]->isPolar,i
				//					g_probes[i]->step
				);

				 */
				/*
				   for(int j=0 ; j< g_probes[i]->nearProbes.size(); j++)
				   printf("%4d,",g_probes[i]->nearProbes[j]);
				 */

				cout << endl;
				idx++;

			}
			cout << "END" << endl;
			for(int i=0; i < g_proteinAtoms.size() ; i++) {
				printf("ATOM  %5d  %-3s %s %s%4d    %8.3f%8.3f%8.3f      %4d",
						i+1,
						prType(i).c_str(),
						ResiName(i).c_str(),
						prChain(i).c_str(),
						ResiIdx(i) + 1,
						prX(i),
						prY(i),
						prZ(i),
						g_proteinAtoms[i]->proteinBc
						//					g_residues[ResiIdx(i)]->maxBc
						//ResiIdx(i)
					  );

				//		for(int j=0 ; j< g_proteinAtoms[i]->nearProbes.size(); j++)
				//			printf("%4d,",g_proteinAtoms[i]->nearProbes[j]);

				cout << endl;
			}

			printf("END\n");
		}



		int initializeCellList() {

			int cellId,newCellId,count=0;
			Cell* temp;

			vector<int>::iterator it;

			double dist;
			g_gridMaxX = ceil((g_maxX-g_minX)/2);
			g_gridMaxY = ceil((g_maxY-g_minY)/2);
			g_gridMaxZ = ceil((g_maxZ-g_minZ)/2);


			/* Initialize cell list */
			for(int i=0; i < g_gridMaxX*g_gridMaxY*g_gridMaxZ; i++) {
				temp = new Cell;
				g_cells.push_back(temp);
			}

			// Insert cell id which have atom into cell list and iPassnsert atom id into cell
			g_atomCellList.clear();
			for (int i=0; i < g_proteinAtoms.size() ; i++) {
				cellId = getCellId(prAtom(i),2);
				g_cells[cellId]->atoms.push_back(i);
				g_atomCellList.push_back(cellId);
			}

			sort(g_atomCellList.begin(),g_atomCellList.end());
			it = unique(g_atomCellList.begin(),g_atomCellList.end());
			g_atomCellList.resize(it - g_atomCellList.begin());


			// Calculate atom-atom distance and assign near atom within 7.7A
			for ( int i=0; i < g_atomCellList.size() ; i++) {
				for(int j=0; j < g_cells[g_atomCellList[i]]->atoms.size(); j++) {
					for(int k=0; k< 614; k++) {
						newCellId = g_atomCellList[i] + ( g_gridProperties[k]->gridIndex[0] + g_gridMaxX*g_gridProperties[k]->gridIndex[1] + g_gridMaxX*g_gridMaxY*g_gridProperties[k]->gridIndex[2]);
						if(newCellId >-1 && newCellId < g_gridMaxX*g_gridMaxY*g_gridMaxZ)
							if(g_cells[newCellId]->atoms.size() >0) { 
								for(int l=0; l < g_cells[newCellId]->atoms.size(); l++) {
									dist = vsDist(prAtom(g_cells[g_atomCellList[i]]->atoms[j]),prAtom(g_cells[newCellId]->atoms[l]));
									if(dist <= 7.7){ // g_probes-protein 3.8A * 2 and for checking bump
										g_proteinAtoms[g_cells[g_atomCellList[i]]->atoms[j]]->nearAtoms.push_back(g_cells[newCellId]->atoms[l]);
										count++;
									}
								}
							}
					}
					g_proteinAtoms[g_cells[g_atomCellList[i]]->atoms[j]]->proteinBc = checkBC(prX(g_cells[g_atomCellList[i]]->atoms[j]),prY(g_cells[g_atomCellList[i]]->atoms[j]),prZ(g_cells[g_atomCellList[i]]->atoms[j]));
					g_residues[ResiIdx(g_cells[g_atomCellList[i]]->atoms[j])]->maxBc += checkBC(prX(g_cells[g_atomCellList[i]]->atoms[j]),prY(g_cells[g_atomCellList[i]]->atoms[j]),prZ(g_cells[g_atomCellList[i]]->atoms[j]));
					g_residues[ResiIdx(g_cells[g_atomCellList[i]]->atoms[j])]->isOk +=1;

				}
			}

			for (int i=0; i < g_proteinAtoms.size() ; i++) {
				sort(g_proteinAtoms[i]->nearAtoms.begin(),g_proteinAtoms[i]->nearAtoms.end());
				it = unique(g_proteinAtoms[i]->nearAtoms.begin(),g_proteinAtoms[i]->nearAtoms.end());
				g_proteinAtoms[i]->nearAtoms.resize(it - g_proteinAtoms[i]->nearAtoms.begin());

			}

			printf("REMARK Cell ID (x,y,z) : %4d%4d%4d , Total Number of 2A Grid Cell : %7d\n",g_gridMaxX,g_gridMaxY,g_gridMaxZ,g_gridMaxX*g_gridMaxY*g_gridMaxZ);
			return count;
		}

		void generateFirstLayer() {
			int a,b,c,i,j,k,n, isPolar,isMetal,check_ok=0,nTriple,nobump,probe_flag_1=0,probe_flag_2=0;
			int proteinKey, new_proteinKey,temp,n_a;
			int cellId,step;

			vector<int>::iterator it;


			double i_R,j_R,k_R,nProbe_radius,distIj,distJk,distIk;
			GridVectors::iterator grid_iter_protein,grid_iter_sub;
			Vector3 Rp1,Rp2;
			Probe* nProbe;

			g_probeCellList.clear();
			g_layer=1;
			for(i=0; i < g_proteinAtoms.size(); i++) {
				for(a=0; a < g_proteinAtoms[i]->nearAtoms.size(); a++) {
					j = g_proteinAtoms[i]->nearAtoms[a];
					if( j > i) {
						for(b=a+1; b < g_proteinAtoms[i]->nearAtoms.size(); b++) {
							k = g_proteinAtoms[i]->nearAtoms[b];
							if( k > j) {
								distIj = calculateDistance(prAtom(i),prAtom(j))+0.0001;
								distIk = calculateDistance(prAtom(i),prAtom(k))+0.0001;
								distJk = calculateDistance(prAtom(j),prAtom(k))+0.0001;

								isMetal =0;
								isPolar = prPolar(i)+prPolar(j)+prPolar(k);
								if(prPolar(i)==2 || prPolar(j)==2 || prPolar(k)==2) {
									isMetal = 1;

								}
								// to check the possibility of accretion for three atoms;  
								check_ok=0;
								n=0;
								if(isPolar > 0) {  // for all cases except three apolar atom triplet;
									for( step=0 ; step < 8; step++) {	
										n= step;
										nProbe_radius = PROBE_RADIUS_POLAR + RADIUS_INCREMENT*n; // (temporary)g_probes radius [1.25 - 1.6 A]

										if(prPolar(i)==1) 
											i_R = PROBE_RADIUS_POLAR + RADIUS_INCREMENT*n; // (temporary)isPolar atom radius [1.25 - 1.60 A] --> to make the distance 2.50 - 3.2 A
										if(prPolar(i)==2) 
											i_R = PROBE_RADIUS_METAL + RADIUS_INCREMENT*n; // (temporary)metal atom radius [0.95 - 1.30 A] --> to make the distance 2.2 - 2.9A
										if(prPolar(i)==0)  
											i_R = prRadius(i) + 0.35 - RADIUS_INCREMENT*n; // apolar atom radius [1.70(Csp2)+0.35 or 1.90(Csp3)+0.35] --> to make the distance 3.30 or 3.50 A

										if(prPolar(j)==1)
											j_R = PROBE_RADIUS_POLAR + RADIUS_INCREMENT*n;
										if(prPolar(j)==2)
											j_R = PROBE_RADIUS_METAL + RADIUS_INCREMENT*n;
										if(prPolar(j)==0)
											j_R = prRadius(j) + 0.35 - RADIUS_INCREMENT*n;

										if(prPolar(k)==1) 
											k_R = PROBE_RADIUS_POLAR + RADIUS_INCREMENT*n;
										if(prPolar(k)==2)
											k_R = PROBE_RADIUS_METAL + RADIUS_INCREMENT*n;
										if(prPolar(k)==0)
											k_R = prRadius(k) + 0.35 - RADIUS_INCREMENT*n;

										if(j_R + k_R + 2*nProbe_radius > distJk && i_R + j_R + 2*nProbe_radius > distIj && i_R + k_R + 2*nProbe_radius > distIk){
											check_ok =1;
											break;
										}
									}

									n_a=0;
									if(!check_ok) {
										for( step=0 ; step < 7 ; step++) {
											n_a = step;
											if(prPolar(i)==0) 	
												i_R = prRadius(i) + RADIUS_INCREMENT*n_a; // 1.7 - 2.0A, 1.9 - 2.2A  to make the distance 3.30 - 3.60 or 3.50 - 3.80 A for Csp2 or Csp3, respectively
											if(prPolar(j)==0)
												j_R = prRadius(j) + RADIUS_INCREMENT*n_a;
											if(prPolar(k)==0)
												k_R = prRadius(k) + RADIUS_INCREMENT*n_a;

											if(j_R + k_R + 2*nProbe_radius > distJk && i_R + j_R + 2*nProbe_radius > distIj && i_R + k_R + 2*nProbe_radius > distIk) {
												check_ok =1;
												break;
											}
										}
									}
								}
								//	if(check_ok != 1)   either unable to accrete in the previous step (nProbe_radius = 1.6/i_R=1.6 or 1.3 for isPolar or metal) or when three atoms are all apolar )
								if(isPolar == 0){  // for three apolar atom cases
									nProbe_radius=1.6;
									for( step=0 ; step < 7 ; step++) {
										n_a = step;
										i_R = prRadius(i) + RADIUS_INCREMENT*n_a; // 1.7 - 2.0A, 1.9 - 2.2A  to make the distance 3.30 - 3.60 or 3.50 - 3.80 A for Csp2 or Csp3, respectively
										j_R = prRadius(j) + RADIUS_INCREMENT*n_a;
										k_R = prRadius(k) + RADIUS_INCREMENT*n_a;

										if(j_R + k_R + 2*nProbe_radius > distJk && i_R + j_R + 2*nProbe_radius > distIj && i_R + k_R + 2*nProbe_radius > distIk) {
											check_ok =1;
											break;
										}
									}
								}

								if(check_ok == 1) {
									if(TRUE==GetTripleVertex(prAtom(i), prAtom(j), prAtom(k),&Rp1, &Rp2, i_R,j_R,k_R, nProbe_radius)) {
										nTriple++;
										Probe* nProbe1 = new Probe;
										Probe* nProbe2 = new Probe;

										nProbe1->point = Rp1; // Assign generated coordinate to new g_probes(position 1)
										nProbe2->point = Rp2; // Assign generated coordinate to new g_probes(position 2)
										nProbe1->radius = nProbe_radius;
										nProbe2->radius = nProbe_radius;
										nProbe1->isPolar = isPolar;
										nProbe2->isPolar = isPolar;
										nProbe1->numLayer = 1;
										nProbe2->numLayer = 1;
										nProbe1->clusterId1 = 1;
										nProbe2->clusterId1 = 1;
										nProbe1->numBc1=0;
										nProbe2->numBc1=0;

										if(prPolar(i)==2 || prPolar(j)==2 || prPolar(k)==2) {
											nProbe1->isPolar = 7;
											nProbe2->isPolar = 7;
										}

										if(bumpcheckAtoms(nProbe1) ) {
											probe_flag_1=1;
											nobump++;
											nProbe1->position = 1;
											nProbe1->isSurvived = 1;
											nProbe1->step = n*100+n_a;
											nProbe1->contactAtoms.push_back(i);
											nProbe1->contactAtoms.push_back(j);
											nProbe1->contactAtoms.push_back(k);
											g_probes.push_back(nProbe1);
											cellId = getCellId(pAtom(g_probes.size()-1),2);
											g_cells[cellId]->probes.push_back(g_probes.size()-1);
											g_probeCellList.push_back(cellId);
										}

										if(bumpcheckAtoms(nProbe2) ) {
											probe_flag_2 = 1;
											nobump++;
											nProbe2->position = 2;
											nProbe2->isSurvived = 1;
											nProbe2->step = n*100+n_a;
											nProbe2->contactAtoms.push_back(i);
											nProbe2->contactAtoms.push_back(j);
											nProbe2->contactAtoms.push_back(k);
											g_probes.push_back(nProbe2);
											cellId = getCellId(pAtom(g_probes.size()-1),2);
											g_cells[cellId]->probes.push_back(g_probes.size()-1);
											g_probeCellList.push_back(cellId);

										}
									}
								} 
							} 
						}
					}
				}
			}

			sort(g_probeCellList.begin(),g_probeCellList.end());
			it = unique(g_probeCellList.begin(),g_probeCellList.end());
			g_probeCellList.resize(it - g_probeCellList.begin());

		}
		void generateSubFirstLayer() {
			int a,b,c,i,j,k,m,n,o,p,q,r, isPolar,isMetal,check_ok=0,nTriple,nobump,probe_flag_1=0,probe_flag_2=0;
			int proteinKey, new_proteinKey,temp;
			double i_R,j_R,k_R,nProbe_radius=1.9,distIj,distJk,distIk;
			int cellId;
			int polar_atom, apolar_atom,step;

			GridVectors::iterator grid_iter_protein,grid_iter_sub;
			Vector3 Rp1,Rp2;
			Probe* nProbe;

			int n_a=0;

			std::vector <int> merged_list;
			vector<int>::iterator it;

			g_layer = 2;
			temp = g_probes.size();
			for(i=0 ; i< temp  ; i++) {
				for(a=0; a < g_probes[i]->nearAtoms.size() ; a++) {   //i; g_probes, j,k ; protein atom
					for(b=a+1; b < g_probes[i]->nearAtoms.size() ; b++) {
						j = g_probes[i]->nearAtoms[a];
						k = g_probes[i]->nearAtoms[b];
						distJk = calculateDistance(prAtom(j),prAtom(k))+0.0001;
						distIj = calculateDistance(prAtom(j),pAtom(i))+0.0001;
						distIk = calculateDistance(prAtom(k),pAtom(i))+0.0001;

						isMetal =0;
						isPolar = prPolar(j)+prPolar(k);
						if(prPolar(j)==2 || prPolar(k)==2) {
							isMetal = 1;
						}

						check_ok=0;
						n=0;

						if(isPolar >0) {
							for( step=0 ; step < 8; step++) {	
								n= step;
								nProbe_radius = PROBE_RADIUS_PROBE + RADIUS_INCREMENT*n;  //PROBE_RADIUS_PROBE = 0.7A

								i_R = PROBE_RADIUS_PROBE - RADIUS_INCREMENT*n;

								if(prPolar(j)==1)
									j_R = 1.8 + RADIUS_INCREMENT*n;
								if(prPolar(j)==2)
									j_R = 1.5 + RADIUS_INCREMENT*n;
								if(prPolar(j)==0)
									j_R = prRadius(j) + 0.9 - RADIUS_INCREMENT*n;

								if(prPolar(k)==1) 
									k_R = 1.8 + RADIUS_INCREMENT*n;
								if(prPolar(k)==2)
									k_R = 1.5 + RADIUS_INCREMENT*n;
								if(prPolar(k)==0)
									k_R = prRadius(k) + 0.9 - RADIUS_INCREMENT*n;

								if(j_R + k_R + 2*nProbe_radius > distJk && i_R + j_R + 2*nProbe_radius > distIj && i_R + k_R + 2*nProbe_radius > distIk) {
									check_ok =1;
									break;
								}

							}
						}
						n_a=0;
						if(check_ok != 1) {
							nProbe_radius=PROBE_RADIUS_PROBE;

							for( step=0 ; step < 7 ; step++) {
								n_a = step;

								i_R = PROBE_RADIUS_PROBE;

								if(prPolar(j)==0)
									j_R = prRadius(j) +0.9 + RADIUS_INCREMENT*n_a;
								if(prPolar(k)==0)
									k_R = prRadius(k) +0.9 + RADIUS_INCREMENT*n_a;

								if(j_R + k_R + 2*nProbe_radius > distJk && i_R + j_R + 2*nProbe_radius > distIj && i_R + k_R + 2*nProbe_radius > distIk) {
									check_ok =1;
									break;
								}
							}

						}


						if(check_ok==1 && TRUE==GetTripleVertex(prAtom(j), prAtom(k), pAtom(i),&Rp1, &Rp2, j_R,k_R,i_R, nProbe_radius)){


							nTriple++;
							Probe* nProbe1 = new Probe;
							Probe* nProbe2 = new Probe;

							nProbe1->point = Rp1; // Assign generated coordinate to new g_probes(position 1)
							nProbe2->point = Rp2; // Assign generated coordinate to new g_probes(position 2)
							nProbe1->radius = nProbe_radius;
							nProbe2->radius = nProbe_radius;
							nProbe1->isPolar = isPolar;
							nProbe2->isPolar = isPolar;
							nProbe1->numLayer = 2;
							nProbe2->numLayer = 2;
							nProbe1->clusterId1 = 1;
							nProbe2->clusterId1 = 1;
							if(prPolar(j)==2 || prPolar(k)==2) {
								nProbe1->isPolar = 7;
								nProbe2->isPolar = 7;
							}

							if( bumpcheckAtoms2(nProbe1)) {

								nobump++;
								nProbe1->position = 1;
								nProbe1->isSurvived = 1;
								nProbe1->step = n*10+n_a;
								nProbe1->contactAtoms.push_back(j);
								nProbe1->contactAtoms.push_back(k);
								nProbe1->contactAtoms.push_back(i);
								g_probes.push_back(nProbe1);
								cellId = getCellId(pAtom(g_probes.size()-1),2);
								g_cells[cellId]->probes.push_back(g_probes.size()-1);
								g_probeCellList.push_back(cellId);
								//							printf("Combi ER : %5d,%5d,%5d\t%5.2f,%5.2f,%5.2f,%5.2f\t%5.2f,%5.2f,%5.2f\n",i,j,k,i_R,j_R,k_R,nProbe_radius,distJk,distIj,distIk);
								//							printf("Coor : %20.18f,%20.18f,%20.18f\n",nProbe1->point.x,nProbe1->point.y,nProbe1->point.z);

							}
							if(bumpcheckAtoms2(nProbe2)) {
								nobump++;
								nProbe2->position = 2;
								nProbe2->isSurvived = 1;
								nProbe2->step = n*10+n_a;
								nProbe2->contactAtoms.push_back(j);
								nProbe2->contactAtoms.push_back(k);
								nProbe2->contactAtoms.push_back(i);
								g_probes.push_back(nProbe2);
								cellId = getCellId(pAtom(g_probes.size()-1),2);
								g_cells[cellId]->probes.push_back(g_probes.size()-1);
								g_probeCellList.push_back(cellId);
								//							printf("Combi ER : %5d,%5d,%5d\t%5.2f,%5.2f,%5.2f,%5.2f\t%5.2f,%5.2f,%5.2f\n",i,j,k,i_R,j_R,k_R,nProbe_radius,distJk,distIj,distIk);
								//                          printf("Coor : %20.18f,%20.18f,%20.18f\n",nProbe2->point.x,nProbe2->point.y,nProbe2->point.z);

							}
						}
					}

				}
			}

			sort(g_probeCellList.begin(),g_probeCellList.end());
			it = unique(g_probeCellList.begin(),g_probeCellList.end());
			g_probeCellList.resize(it - g_probeCellList.begin());


			assignAtomProbe();
		}

		void generateSubSecondLayer() {
			int a,b,c,i,j,k,m,n,o,p,q,r, isPolar,isMetal,check_ok,nTriple,nobump,probe_flag_1,probe_flag_2;
			int proteinKey, new_proteinKey, new_probeKey,key_15;
			double i_R,j_R,k_R,nProbe_radius=1.9,distIj,distJk,distIk;
			GridVectors::iterator grid_iter_protein,grid_iter_sub;
			GridVectors::iterator grid_iter_probe,grid_iter_15;
			Vector3 Rp1,Rp2;
			Probe* nProbe;
			int cellId,step,step_idx;
			std::vector <int> merged_list;
			vector<int>::iterator it;
			g_layer=3;

			for(i=0; i < g_proteinAtoms.size() ; i++) {

				a=g_proteinAtoms[i]->nearProbes.size();

				for(p= 0 ; p < a ; p++) {
					for(q= p+1 ; q < a ; q++) {

						j = g_proteinAtoms[i]->nearProbes[p];
						k = g_proteinAtoms[i]->nearProbes[q];

						if( j != k) {

							distJk = calculateDistance(pAtom(j),pAtom(k))+0.0001;
							distIj = calculateDistance(prAtom(i),pAtom(j))+0.0001;
							distIk = calculateDistance(prAtom(i),pAtom(k))+0.0001;

							isPolar = prPolar(i);
							check_ok=0;



							if( isPolar )
								step_idx=8;
							else
								step_idx=7;
							for( step=0 ; step < step_idx; step++) 	{
								n = step;
								nProbe_radius = PROBE_RADIUS_PROBE + RADIUS_INCREMENT*n;

								if(prPolar(i)==1 )
									i_R = 1.8 + RADIUS_INCREMENT*n;
								if(prPolar(i)==2)
									i_R = 1.5 + RADIUS_INCREMENT*n;
								if(prPolar(i)==0)
									i_R = prRadius(i) + 0.9;

								j_R = PROBE_RADIUS_PROBE - RADIUS_INCREMENT*n;
								k_R = PROBE_RADIUS_PROBE - RADIUS_INCREMENT*n;


								if(j_R + k_R + 2*nProbe_radius > distJk && i_R + j_R + 2*nProbe_radius > distIj && i_R + k_R + 2*nProbe_radius > distIk) {
									check_ok = 1;
									break;
								}

							} 


							if(check_ok==1 && TRUE==GetTripleVertex(prAtom(i), pAtom(j), pAtom(k),&Rp1, &Rp2, i_R,j_R,k_R, nProbe_radius)){

								nTriple++;
								Probe* nProbe1 = new Probe;
								Probe* nProbe2 = new Probe;

								nProbe1->point = Rp1; // Assign generated coordinate to new g_probes(position 1)
								nProbe2->point = Rp2; // Assign generated coordinate to new g_probes(position 2)
								nProbe1->radius = nProbe_radius;
								nProbe2->radius = nProbe_radius;
								nProbe1->isPolar = isPolar;
								nProbe2->isPolar = isPolar;
								nProbe1->numLayer = 3;
								nProbe2->numLayer = 3;
								nProbe1->clusterId = 1;
								nProbe2->clusterId = 1;
								if(prPolar(i)==2) {
									nProbe1->isPolar = 7;
									nProbe2->isPolar = 7;
								}

								if( bumpcheckAtoms2(nProbe1)) {
									probe_flag_1 =1;
									nobump++;
									nProbe1->position = 1;
									nProbe1->isSurvived = 1;
									nProbe1->step=n;
									nProbe1->contactAtoms.push_back(i);
									nProbe1->contactAtoms.push_back(j);
									nProbe1->contactAtoms.push_back(k);
									g_probes.push_back(nProbe1);
									cellId = getCellId(pAtom(g_probes.size()-1),2);
									g_cells[cellId]->probes.push_back(g_probes.size()-1);
									g_probeCellList.push_back(cellId);

								}

								if(bumpcheckAtoms2(nProbe2)) {
									probe_flag_2=1;
									nobump++;
									nProbe2->position = 2;
									nProbe2->isSurvived = 1;
									nProbe2->step=n;
									nProbe2->contactAtoms.push_back(i);
									nProbe2->contactAtoms.push_back(j);
									nProbe2->contactAtoms.push_back(k);
									g_probes.push_back(nProbe2);

									cellId = getCellId(pAtom(g_probes.size()-1),2);
									g_cells[cellId]->probes.push_back(g_probes.size()-1);
									g_probeCellList.push_back(cellId);


								}
							}

						}

					}
				}
			}

			sort(g_probeCellList.begin(),g_probeCellList.end());
			it = unique(g_probeCellList.begin(),g_probeCellList.end());
			g_probeCellList.resize(it - g_probeCellList.begin());

		}

		void generateNextLayer() {

			int a,b,c,i,j,k,m,n=0,o,p,q,r,total_layer=0, isPolar,isMetal,check_ok=0,nTriple,nobump,probe_flag_1,probe_flag_2;
			int probeKey, new_probeKey, key_35;
			int start=0, end = g_probes.size();
			int cellId;
			double i_R,j_R,k_R,nProbe_radius=1.9,distIj,distJk,distIk;
			GridVectors::iterator grid_iter_protein,grid_iter_sub,grid_iter_35;
			GridVectors::iterator grid_iter_probe,grid_iter_probe_sub;
			Vector3 Rp1,Rp2;
			Probe* nProbe;

			std::vector <int> probeIds,nprobe_ids;
			vector<int>::iterator it;

			probeIds.clear();

			for(i=0 ; i< g_probes.size(); i++) {
				probeIds.push_back(i);
			}

			assignProbePair_sub(probeIds);
			nobump = 0;

			do {
				nprobe_ids.clear();
				end = g_probes.size();

				if(nobump != 0) 
					start = end - nobump+1;
				nobump = 0;

				for( i=start ; i < end; i++) {
					for( o = 0; o < g_probes[i]->nearProbes.size() ; o++) {
						j = g_probes[i]->nearProbes[o];
						for(p = o+1 ; p < g_probes[i]->nearProbes.size() ; p++) {
							k = g_probes[i]->nearProbes[p];
							if( j!=k && i!=j && i!=k) {

								distJk = calculateDistance(pAtom(j),pAtom(k))+0.0001;
								distIj = calculateDistance(pAtom(i),pAtom(j))+0.0001;
								distIk = calculateDistance(pAtom(i),pAtom(k))+0.0001;

								if(distJk <= 2.8 && distIj <= 2.8 && distIk <= 2.8) {

									nProbe_radius = i_R = j_R = k_R = PROBE_RADIUS_PROBE;

									probe_flag_1 = probe_flag_2 = check_ok = n = 0; 

									if(j_R + k_R + 2*nProbe_radius > distJk && i_R + j_R + 2*nProbe_radius > distIj && i_R + k_R + 2*nProbe_radius > distIk)
										check_ok=1;

									if( check_ok ==1 & TRUE==GetTripleVertex(pAtom(i), pAtom(j), pAtom(k),&Rp1, &Rp2, i_R,j_R,k_R, nProbe_radius)){
										nTriple++;
										Probe* nProbe1 = new Probe;
										Probe* nProbe2 = new Probe;

										nProbe1->point = Rp1; // Assign generated coordinate to new g_probes(position 1)
										nProbe2->point = Rp2; // Assign generated coordinate to new g_probes(position 2)
										nProbe1->radius = nProbe_radius;
										nProbe2->radius = nProbe_radius;
										nProbe1->isPolar = 0;
										nProbe2->isPolar = 0;
										nProbe1->numLayer = total_layer+4;
										nProbe2->numLayer = total_layer+4;
										nProbe1->clusterId = 0;
										nProbe2->clusterId = 0;


										if( bumpcheckAtoms3(nProbe1,i,j,k)) {
											probe_flag_1=1;
											nobump++;
											nProbe1->position = 1;
											nProbe1->isSurvived = 1;
											nProbe1->contactAtoms.push_back(i);
											nProbe1->contactAtoms.push_back(j);
											nProbe1->contactAtoms.push_back(k);
											g_probes.push_back(nProbe1);
											nprobe_ids.push_back(g_probes.size()-1);
											cellId = getCellId(pAtom(g_probes.size()-1),2);
											g_cells[cellId]->probes.push_back(g_probes.size()-1);
											g_probeCellList.push_back(cellId);

										}

										if( bumpcheckAtoms3(nProbe2,i,j,k)) {
											probe_flag_2=1;
											nobump++;
											nProbe1->position = 2;
											nProbe1->isSurvived = 1;
											nProbe1->contactAtoms.push_back(i);
											nProbe1->contactAtoms.push_back(j);
											nProbe1->contactAtoms.push_back(k);
											g_probes.push_back(nProbe2);
											nprobe_ids.push_back(g_probes.size()-1);
											cellId = getCellId(pAtom(g_probes.size()-1),2);
											g_cells[cellId]->probes.push_back(g_probes.size()-1);
											g_probeCellList.push_back(cellId);

										}
									}
								}
							}
						}

					}
				}
				sort(g_probeCellList.begin(),g_probeCellList.end());
				it = unique(g_probeCellList.begin(),g_probeCellList.end());
				g_probeCellList.resize(it - g_probeCellList.begin());

				assignProbePair_sub(nprobe_ids);

				printf ("REMARK %3d g_layer (%3d - %3d) , Survived : %3d, Total g_probes : %5d\n",total_layer+4,start,end,nobump,g_probes.size());
				total_layer++;
			} while(nobump > 0 && total_layer < MAX_LAYER-3);
			assignProbePair_final();
		}
