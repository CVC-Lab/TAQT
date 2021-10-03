#include <math.h>

#include "dual.h"
#include "disjointset.h"

#include <fstream>
#include <iostream>

using namespace std;

float DualNode::all_min = 0;
float DualNode::all_max = 1;

DualGraph::DualGraph(float min, float max, int nr, float t_vol) {
	nrang = nr;
	if(nr > 0) {
		nodes = new dynamic_array<DualNode*>[nrang];
		cut_vals = new float[nrang+1];
	} else {
		nodes = NULL;
		cut_vals = NULL;
	}
	nnodes = 0;
	m_min = min;
	m_max = max;
	total_vol = t_vol;
	data = NULL;
}

void DualGraph::build(VolumeReg3Critical* p_data, AugmentedContourTree* actree)
{
	const static double MIN_THRE = 1e-6; 
	data = p_data;
	int i, j;
	CutVertex **new_verts = (CutVertex **)malloc(sizeof(CutVertex *)*(nrang+1));

	int *nverts = new int[nrang+1];
	int dim[3];
	data->getDimension(dim);
	int max_id = MAX(dim[0]*dim[1]*dim[2], actree->getMaxID());
	int cur_id = max_id;
	float factor = nrang/(m_max - m_min);
	DualNode::setGlobalMinMax(m_min, m_max);
	map<int, CutVertex*> *cut_maps = new map<int, CutVertex*>[nrang+1];
	cut_vals[0] = m_min; cut_vals[nrang] = m_max;
	for(i = 1; i < nrang; i++) {
		cut_vals[i] = m_min + i*(m_max-m_min)/nrang;
	}
	for (i = 0; i <= nrang; i++) {
		new_verts[i] = actree->cut(data, cut_vals[i], nverts[i]);
#ifdef _DEBUG
		printf("# of cuts at value %f: %d\n", cut_vals[i], nverts[i]);
#endif
		for (j = 0; j < nverts[i]; j++) {
			new_verts[i][j].id = cur_id ++;
			(cut_maps[i])[new_verts[i][j].super_arc->object.id] = &(new_verts[i][j]);
		}
	}
	//assert(nverts[nrang] == 0);

	int *levels = new int[max_id];
	dynamic_array<int> *cpnt_sets = new dynamic_array<int>[nrang];
	// find the level of all critical points
	LL_Node<CriticalPoint>* ptr = actree->node_ids.head();
	while (ptr != NULL) {
		int cid = ptr->object.id;
		//printf("vid: %d val = %f\n", cid, ptr->object.val);
		ptr->object.val = data->getValue(cid);
		
        levels[cid] = getCutLevel(ptr->object.val);
        if(levels[cid] >= 0 && levels[cid] < nrang) {
			(cpnt_sets[levels[cid]]).insert(cid);
		}
		ptr = ptr->next;
	}
	
	for (i = 0; i < nrang; i++) {
		buildRange(i, cpnt_sets[i], new_verts[i], nverts[i], new_verts[i+1], nverts[i+1],
				   cut_maps, levels, p_data, actree);
	}

	delete[] levels;
	// connect the neighboring dual nodes
	for (i = 1; i < nrang; i++) {
		for (j = 0; j < nverts[i]; j++) {
			DualNode::connect(new_verts[i][j].ngb_l, new_verts[i][j].ngb_u);
		}
	}
	//printf("***********total # of nodes************** = %d\n", nnodes);
	delete[] cut_maps;
	delete[] nverts;
	delete[] cpnt_sets;
	for(i = 0; i <= nrang; i++) {
		delete[] new_verts[i];
	}
	free(new_verts);
}

void DualGraph::buildRange(int ir, dynamic_array<int>& crit_pnts, CutVertex* low_verts, int nl,
						   CutVertex* up_verts, int nu, map<int, CutVertex*> *cut_maps, int levels[], 
						   VolumeReg3Critical* p_data, AugmentedContourTree* actree)
{
	int j, k;
	int total = crit_pnts.size() + nl + nu;
	DisjointSet *djs = new DisjointSet(total);
	map<int, int> refs;
	for (j = 0; j < nl; j++) {
		refs[low_verts[j].id] = j;
		djs->makeSet(j);
	}
	for (k = 0; k < crit_pnts.size(); k++, j++) {
		refs[crit_pnts[k]] = j;
		djs->makeSet(j);
	}
	for (k = 0; k < nu; k++, j++) {
		refs[up_verts[k].id] = j;
		djs->makeSet(j);
	}

	for (j = 0; j < crit_pnts.size(); j++) {
		int cid = crit_pnts[j];
		SuperNode *snode = (actree->nodes)[cid];
		int ngb = snode->degree();
		for (k = 0; k < ngb; k++) {
			int nid = snode->getNeighbor(k);
			if (levels[nid] == ir) {
				if (refs[nid] < refs[cid]) {
					djs->unionSet(refs[cid], refs[nid]);
				}
			} else if (levels[nid] < ir) {			// lower neighbor crossing the sub-range
				LL_Node<SuperArc>* sarc = snode->getArcPtr(nid);
				nid = ((cut_maps[ir])[sarc->object.id])->id;
				assert(refs[nid] < refs[cid]);
				djs->unionSet(refs[cid], refs[nid]);
			}
		}
	}
	for (j = 0; j < nu; j++) {
		LL_Node<SuperArc>* sarc = up_verts[j].super_arc;
		int nid;
		if(data->getValue(sarc->object.v2) <= data->getValue(sarc->object.v1)){
			nid = sarc->object.v2;
		} else {
			nid = sarc->object.v1;
		}
		if (levels[nid] == ir) {		// >= because the cricital point nid can be on the cutting edge
				djs->unionSet(refs[up_verts[j].id], refs[nid]);
		} else { // would not consider if nid is above current segment 
				
				if(levels[nid] > ir) {
									printf("ir = %d, levels[%d] = %d\n", ir, nid, levels[nid]);
									printf("critical val = %f\n", data->getValue(nid));
									printf("n1 = %d val1 = %f, n2 = %d val2 = %f\n", 
											sarc->object.v1, data->getValue(sarc->object.v1),
											sarc->object.v2, data->getValue(sarc->object.v2));
								}
				assert(levels[nid] < ir);
				nid = ((cut_maps[ir])[sarc->object.id])->id;
				djs->unionSet(refs[up_verts[j].id], refs[nid]);
			}
		}

	map<int, DualNode*> dualrefs;
	int dual_id;
	DualNode *dualnode;
	for (j = 0; j < nl; j++) {
		dual_id = djs->find(j);
		if (dualrefs.find(dual_id) == dualrefs.end()) {
			dualnode = new DualNode(ir);
			dualnode->id = nnodes ++;
			dualnode->bl = low_verts[j].super_arc->object.getBettiNumbers();
			dualnode->setFuncMax(cut_vals[ir]);
			dualnode->setFuncMin(cut_vals[ir]);
			dualrefs[dual_id] = dualnode;
			nodes[ir].insert(dualnode);
			low_verts[j].ngb_u = dualnode;
		} else {
			dualnode = dualrefs[dual_id];
			dualnode->bl += low_verts[j].super_arc->object.getBettiNumbers();
			low_verts[j].ngb_u = dualnode;
		}
	}
	for (j = 0; j < nu; j++) {
		dual_id = djs->find(j + nl + crit_pnts.size());
		if (dualrefs.find(dual_id) == dualrefs.end()) {
			dualnode = new DualNode(ir);
			dualnode->id = nnodes ++;
			dualnode->bu = up_verts[j].super_arc->object.getBettiNumbers();
			///dualnode->setFuncMax((ir+1)/(float)nrang);
			///dualnode->setFuncMin((ir+1)/(float)nrang);
			dualnode->setFuncMax(cut_vals[ir+1]);
			dualnode->setFuncMin(cut_vals[ir+1]);
			dualrefs[dual_id] = dualnode;
			nodes[ir].insert(dualnode);
			up_verts[j].ngb_l = dualnode;
		} else {
			dualnode = dualrefs[dual_id];
			///dualnode->setFuncMax((ir+1)/(float)nrang);
			dualnode->setFuncMax(cut_vals[ir+1]);
			dualnode->bu += up_verts[j].super_arc->object.getBettiNumbers();
			up_verts[j].ngb_l = dualnode;
		}
	}

	// set the function min and max of dual nodes
	for (k = 0, j = nl; k < crit_pnts.size(); k++, j++) {
			dual_id = djs->find(j);
			dualnode = dualrefs[dual_id];
			CriticalPoint* cpnt = &((*(actree->p_cpmap))[crit_pnts[k]]->object);
			///dualnode->setFuncMax((cpnt->val-m_min) / (m_max-m_min));
			///dualnode->setFuncMin((cpnt->val-m_min) / (m_max-m_min));
			dualnode->setFuncMax(cpnt->val);
			dualnode->setFuncMin(cpnt->val);
			dualnode->ncrit ++;
	}

	// set the volumes of dual nodes
	for (j = 0; j < nl; j++) {
		dual_id = djs->find(j);
		dualnode = dualrefs[dual_id];
		int tid;

		float cut_val;
		if (!dualnode->isVolumeComputed()) {
			LL_Node<SuperArc> *arc = low_verts[j].super_arc;
			int v1, v2;
			int idx[3], nidx[3];
			float range[2];
			///range[0] = m_min + dualnode->fmin*(m_max-m_min);
			///range[1] = m_min + dualnode->fmax*(m_max-m_min);
			///cut_val = m_min + ir*(m_max-m_min)/nrang;
			range[0] = dualnode->fmin;
			range[1] = dualnode->fmax;
			cut_val = cut_vals[ir];
			if (arc->object.intra_node < 0) {
				if(data->getValue(arc->object.v2) < data->getValue(arc->object.v1)) {
					v1 = arc->object.v2;
					v2 = arc->object.v1;
				} else {
					v1 = arc->object.v1;
					v2 = arc->object.v2;
				}
				
				if (!data->areNeighbors(v1, v2)) {
					int xid;
#ifdef _DEBUG
					data->id2Index(v1, idx);
					data->id2Index(v2, nidx);

					printf("v1 : (%d %d %d) = %f, ", idx[0], idx[1], idx[2], data->getValue(v1));
					printf("v2 : (%d %d %d) = %f, ", nidx[0], nidx[1], nidx[2], data->getValue(v2));
#endif
					SuperNode* snode1 = actree->getNode(v1);
					int x1 = snode1->getNeighborTag(v2);
					SuperNode* snode2 = actree->getNode(v2);
					int x2 = snode2->getNeighborTag(v1);
					if (data->getValue(x1) > data->getValue(v1) && data->areNeighbors(v1, x1)) {
						xid = x1;
						v2 = xid;
					} else if(data->getValue(x2) < data->getValue(v2) && data->areNeighbors(v2, x2)) {
						xid = x2;
						v1 = xid;
					}
				}
				assert(data->areNeighbors(v1, v2));
				//printf("v1 degree = %d, v2 degree = %d\n", (actree->nodes)[v1]->degree(),
				//   (actree->nodes[v2])->degree());
			} else {
				v1 = arc->object.intra_node;
				v2 = data->getCutEdge(cut_val, range, v1);
				if (data->getValue(v1) > data->getValue(v2)) {
					int tmp = v1;
					v1 = v2;
					v2 = tmp;
				}
				
#ifdef _DEBUG
			printf("v1: %d %f, v2: %d %f, cut_val: %f, range: (%f %f)\n", v1, data->getValue(v1),
					v2, data->getValue(v2), cut_val, range[0], range[1]);
#endif
			}
			data->id2Index(v1, idx);
			data->id2Index(v2, nidx);
			//printf("idx : (%d %d %d) %f, nidx : (%d %d %d) %f\n", idx[0], idx[1], idx[2], data->getValue(v1),
			//		nidx[0], nidx[1], nidx[2], data->getValue(v2));
			tid = data->getTetraIndex(idx, nidx);

			
#ifdef _MOMENTS
			dualnode->moments = data->calcMoments(range, tid);
			dualnode->volume = dualnode->moments.vol / data->getVolume();
			dualnode->fint = dualnode->moments.Fint / data->getVolume();
#else			
			float fint;
			dualnode->volume = data->calcVolume(range, tid, fint);
			dualnode->fint = fint;		
			//
			if(fint >= 100000) {
				printf("dualnode %d : (%f %f) vol = %f, fint = %f\n", dualnode->id, range[0], range[1], 
					dualnode->volume, fint);
			}
#endif
			//assert(dualnode->volume >= 0);
		}
	}
		
	for(j = 0; j < nu; j++) {
		dual_id = djs->find(j + nl + crit_pnts.size());
		dualnode = dualrefs[dual_id];
		int tid;
		if(!dualnode->isVolumeComputed()) {
			LL_Node<SuperArc>* arc = up_verts[j].super_arc;
			int v1, v2;
			int idx[3], nidx[3];
			float range[2];
			///range[0] = m_min + dualnode->fmin*(m_max-m_min);
			///range[1] = m_min + dualnode->fmax*(m_max-m_min);
			range[0] = dualnode->fmin;
			range[1] = dualnode->fmax;
			if(arc->object.intra_node < 0) {
				if(data->getValue(arc->object.v2) < data->getValue(arc->object.v1)) {
					v1 = arc->object.v2;
					v2 = arc->object.v1;
				} else {
					v1 = arc->object.v1;
					v2 = arc->object.v2;
				}

				if (!data->areNeighbors(v1, v2)) {
					int xid;
#ifdef _DEBUG
					data->id2Index(v1, idx);
					data->id2Index(v2, nidx);
					printf("*******v1 : (%d %d %d) = %f, ", idx[0], idx[1], idx[2], data->getValue(v1));
					printf("v2 : (%d %d %d) = %f, ", nidx[0], nidx[1], nidx[2], data->getValue(v2));
#endif
					SuperNode* snode1 = actree->getNode(v1);
					int x1 = snode1->getNeighborTag(v2);
					SuperNode* snode2 = actree->getNode(v2);
					int x2 = snode2->getNeighborTag(v1);
					if (data->getValue(x1) > data->getValue(v1) && data->areNeighbors(v1, x1)) {
						xid = x1;
						v2 = xid;
					} else if(data->getValue(x2) < data->getValue(v2) && data->areNeighbors(v2, x2)) {
						xid = x2;
						v1 = xid;
					}
				}
			} else {
				v1 = arc->object.intra_node;
				float cut_val = cut_vals[ir+1];
				v2 = data->getCutEdge(cut_val, range, v1);
				if (data->getValue(v1) > data->getValue(v2)) {
					int tmp = v1;
					v1 = v2;
					v2 = tmp;
				}
#ifdef _DEBUG
				if(cut_val == 0) printf("up v1: %d %f, v2: %d %f, cut_val: %f, range: (%f %f)\n", v1, data->getValue(v1),
						v2, data->getValue(v2), cut_val, range[0], range[1]);
#endif
			}
			data->id2Index(v1, idx);
			data->id2Index(v2, nidx);

			//printf("idx : (%d %d %d) %f, nidx : (%d %d %d) %f\n", idx[0], idx[1], idx[2], data->getValue(v1),
			//		nidx[0], nidx[1], nidx[2], data->getValue(v2));
			tid = data->getTetraIndex(idx, nidx);
			
#ifdef _MOMENTS
			dualnode->moments = data->calcMoments(range, tid);
			dualnode->volume = dualnode->moments.vol / data->getVolume();
			dualnode->fint = dualnode->moments.Fint / data->getVolume();
#else
			float fint;
			dualnode->volume = data->calcVolume(range, tid, fint);
			dualnode->fint = fint;
#endif
			//printf("dualnode %d : (%f %f) vol = %f\n", dualnode->id, range[0], range[1], dualnode->volume);
			//assert(dualnode->volume >= 0);
		}
	}

	delete djs;
}

DualGraph* DualGraph::mergeRanges()
{
	// assume nrange is power of 2
	DualGraph* dual = new DualGraph(m_min, m_max, nrang/2);
	dual->m_min = m_min;
	dual->m_max = m_max;
	
	int i, j, k, n;

	map<int, int> refs;
	map<int, DualNode *> dualrefs;
	int cur_id = 0;
	for (i = 0; i < nrang/2; i++) {
		DisjointSet *djs = new DisjointSet(nodes[2*i].size() + nodes[2*i+1].size());
		for (j = 0; j < nodes[2*i].size(); j++) {
			refs[(nodes[2*i])[j]->id] = j;
			djs->makeSet(j);
		}
		for (k = 0; k < nodes[2*i+1].size(); k++, j++) {
			refs[(nodes[2*i+1])[k]->id] = j;
			djs->makeSet(j);
			int ngb = (nodes[2*i+1])[k]->lowerDegree();
			for (n = 0; n < ngb; n++) {
				int nid = (nodes[2*i+1])[k]->lowerNeighbor(n)->id;
				djs->unionSet(j, refs[nid]);
			}
		}
		for (j = 0; j < (nodes[2*i]).size(); j++) {
			DualNode *dualnode;
			int dual_id = djs->find(j);
			if (dualrefs.find(dual_id) == dualrefs.end()) {
				dualnode = new DualNode(i);
				dualnode->id = cur_id++;
				dualrefs[dual_id] = dualnode;
				dual->insert(i, dualnode);
			} else {
				dualnode = dualrefs[dual_id];
			}
			dualnode->merge((nodes[2*i])[j], 0);
		}

		for (j = 0; j < (nodes[2*i+1]).size(); j++) {
			DualNode *dualnode;
			int dual_id = djs->find(j+nodes[2*i].size());
			if (dualrefs.find(dual_id) == dualrefs.end()) {
				dualnode = new DualNode(i);
				dualnode->id = cur_id++;
				dualrefs[dual_id] = dualnode;
				dual->insert(i, dualnode);
			} else {
				dualnode = dualrefs[dual_id];
			}
			dualnode->merge((nodes[2*i+1])[j], 1);
		}
		refs.clear();
		dualrefs.clear();
		delete djs;
	}

	// connect the coarsers dual graph
	for (i = 0; i < nrang/2-1; i++) {
		for (j = 0; j < nodes[2*i+1].size(); j++) {
			int ngb = (nodes[2*i+1])[j]->upperDegree();
			for (n = 0; n < ngb; n++) {
				DualNode::connect((nodes[2*i+1])[j]->parent, (nodes[2*i+1])[j]->upperNeighbor(n)->parent);
				//(nodes[2*i+1])[j]->parent->addUpperNeighbor((nodes[2*i+1])[j]->ngb_us[n]->parent);
				//(nodes[2*i+1])[j]->ngb_us[n]->parent->addLowerNeighbor((nodes[2*i+1])[j]->parent);
			}
		}
	}
	return dual;
}

void DualGraph::print()
{
	int i, j;
	int count = 0;
	for (i = 0; i < nrang; i++) {
		printf("Range %d .........................................\n", i);
		for (j = 0; j < nodes[i].size(); j++) {
			(nodes[i])[j]->print();
			count ++;
		}
	}
	printf("The number of nodes: %d\n", count);
}

void DualGraph::dump(const char* fname)
{
	ofstream out(fname);
	out << 3 << endl;
	
	int i, j;
	int nr = nRanges();
	int nn = totalNodes();

	fprintf(stderr, "\t\t dump(): %f\n\n", total_vol);
	out << total_vol << endl;
	
	out << nr << " " << nn << endl;
	out << m_min << " " << m_max << endl;
	for(i = 0; i < nr; i++) {
		out << nodes[i].length() << " ";
	}
	out << endl;

	for(i = 0; i < nr; i++) {
		for(j = 0; j < nNodes(i); j++) {
			out << *((nodes[i])[j]);
		}
	}
	out.close();		
}

void DualGraph::prune(float threshold) 
{
	int i, nd = nnodes;
	
	// The coarsest level
	if(nrang == 1) {
		dynamic_array<DualNode *> tnode;
		tnode.insert(nodes[0]);
		nodes[0].clear();

		for(i = 0; i < nd; i++) {
			if(tnode[i]->volume >= threshold) {
				//printf("accepted volume = %f\n", tnode[i]->volume);
				nodes[0].insert(tnode[i]);
			} else {
				//printf("pruned volume = %f\n", tnode[i]->volume);
				pruned.insert(tnode[i]);
				tnode[i]->pruned = true;
				nnodes --;
			}
		}
	} else {
		// For other levels, a node is pruned if its parent is pruned
		int nr;
		for(nr = 0; nr < nrang ; nr++) {
			dynamic_array<DualNode *> tnode;
			tnode.insert(nodes[nr]);
			nd = nodes[nr].length();
			nodes[nr].clear();
			for(i = 0; i < nd; i++) {
				if(tnode[i]->parent->pruned) {
					// prune
					pruned.insert(tnode[i]);
					tnode[i]->pruned = true;
					nnodes --;
				} else {
					nodes[nr].insert(tnode[i]);
				}
			}
		}
	}
}


bool DualGraph::loadDualFile(const char* fname)
{
	ifstream in(fname);
	if(!in.is_open()) {				
		//printf("\n loadDualFile is_open segmentation fault\n");
		cerr << "cannot open file " << fname << endl;		
		return false;
	}
	
	int head;
	in >> head;
	if (head != 3) {
		cerr << "Dual contour tree version mismatch, aborted\n";
		cerr << "Verison 3 is expected\n";
		in.close();
		return false;
	}

	in >> total_vol;
	fprintf(stderr, "\t\t loadDualFile(): %f\n\n", total_vol);
	in >> nrang >> nnodes;
	in >> m_min >> m_max;

	clean();
	int i, j, k;
	nodes = new dynamic_array<DualNode*>[nrang];
	cut_vals = new float[nrang+1];
	DualNode **arry = (DualNode **) malloc(sizeof(DualNode*)*nnodes);
	for(i = 0; i < nnodes; i++) {
		arry[i] = new DualNode;
	}

	for(i = 0, k = 0; i < nrang; i++) {
		int nn; 
		in >> nn;
		for(j = 0; j < nn; j++, k++) {
			nodes[i].insert(arry[k]);
		}
	}

	for(i = 0; i < nrang; i++) {
		for(j = 0; j < nodes[i].length(); j++) {
			DualNode* dualnode = (nodes[i])[j];
			//in >> dualnode->id >> dualnode->level;
			//in >> dualnode->fmin >> dualnode->fmax >> dualnode->volume >> dualnode->ncrit >> dualnode->fint;
			//in >> dualnode->bl >> dualnode->bu;
			//in >> dualnode->moments;
			in >> *dualnode;
			int nl, nu, nid;
			in >> nl >> nu;
			for(k = 0; k < nl; k++) {
				in >> nid;
				dualnode->addLowerNeighbor(arry[nid]);
			}
			for(k = 0; k < nu; k++) {
				in >> nid;
				dualnode->addUpperNeighbor(arry[nid]);
			}
		}
	}
	in.close();
	free(arry);
	return true;
}

void DualGraph::clean()
{
	int i, j;
	for(i = 0; i < nrang && nodes != NULL; i++) {
		for(j = 0; j < nodes[i].length(); j++) {
			delete (nodes[i])[j];
		}
	}
	delete[] nodes;
	delete[] cut_vals;
	for(i = 0; i < pruned.size(); i++) {
		delete pruned[i];
	}
}

DualGraph::~DualGraph()
{
	clean();	
}
