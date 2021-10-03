#include "contourtree.h"
#include "VolCritical.h"
#include "CriticalPoint.h"
#include "disjointset.h"

const int VolCritical::MAX_NEI_NUM = 14;

VolCritical:: ~VolCritical(void) {
	if (p_vcp) delete p_vcp;
}
 
 ContourTree* VolCritical::computeContourTree() {

	ContourTree* p_jtree = joinTree();
	ContourTree* p_stree = splitTree();

	ContourTree* p_ctree = ContourTree::mergeTree(p_jtree, p_stree);
	delete p_jtree;
	delete p_stree;
	return p_ctree;
}

 /**
  *
  */
 void VolCritical::sortCriticalPoints(vector<CriticalPoint>*& p_v)
{
	int i;
	if (p_v != NULL) {
		delete p_v;
		p_v = NULL;
	}

	p_v = new vector<CriticalPoint>;

	int nv = getNVerts();
	for (i = 0; i < nv; i++) {
		float val = getValue(i);
		CriticalPoint cp(i, val);
		p_v->push_back(cp);
	}
	sort(p_v->begin(), p_v->end());
	}


ContourTree* VolCritical::joinTree(void)
{
	ContourTree* p_tree = new ContourTree();
	int i, nv = getNVerts();

	if (p_vcp == NULL) {
		//sortCriticalPoints(p_vcp);
		p_vcp = calcLUStars();
	}
	int *p_map = new int[nv];
	for (i = 0; i < nv; i++) {
		p_map[(*p_vcp)[i].id] = i;
	}

	DisjointSet* djs = new DisjointSet(nv);
	vector<CriticalPoint>::iterator it;
	for (it = p_vcp->begin(); it != p_vcp->end(); ++it) {
		int id = (*it).id;
		p_tree->addNode(*it);
		if (isMinima(*it)) {
			djs->makeSet(p_map[id]);
			continue;
		}
		// (*it) is not a minimum point. Each vertex has 14 neighboring vertices
		int nNeighbors, neighbors[MAX_NEI_NUM];
		nNeighbors = findNeighbors(id, neighbors); 
		for (int j = 0; j < nNeighbors; j++) {
			int nid = neighbors[j];
			if (p_map[nid] < p_map[id]) {
				int mid, pj;
				int pi = djs->find(p_map[id]);
				mid = pj = djs->find(p_map[nid]);
				if (pi == pj) continue;
				if (pi == -1) djs->addElement(p_map[id], pj);
				else mid = djs->link(pi, pj);
				int k = (*p_vcp)[djs->getMaxNode(pj)].id;
#ifdef _DEBUG
				//printf("add arc (%d %d)\n", id, k);
#endif			
				p_tree->addArc(id, k, nid);
				djs->setMaxNode(p_map[id], mid);
			}
		}
	}

	delete djs;
	delete[] p_map;

	p_tree->setMaxID(getNVerts());
	return p_tree;
}

ContourTree* VolCritical::splitTree(void)
{
	ContourTree* p_tree = new ContourTree();
	int i, nv = getNVerts();

	if (p_vcp == NULL) {
		p_vcp = calcLUStars();
	}

	int *p_map = new int[nv];
	for (i = 0; i < nv; i++) {
		p_map[(*p_vcp)[i].id] = nv - 1 - i;
	}

	DisjointSet djs(nv);
	vector<CriticalPoint>::reverse_iterator rit;
	for (rit = p_vcp->rbegin(); rit != p_vcp->rend(); ++rit) {
		int id = (*rit).id;
		p_tree->addNode(*rit);
		if (isMaxima(*rit)) {
			djs.makeSet(p_map[id]);
			continue;
		}
		
		int nNeighbors, neighbors[MAX_NEI_NUM];
		nNeighbors = findNeighbors(id, neighbors); 
		for (int j = 0; j < nNeighbors; j++) {
			int nid = neighbors[j];
				if (p_map[nid] < p_map[id]) {
					int mid, pj;
					int pi = djs.find(p_map[id]);
					mid = pj = djs.find(p_map[nid]);
					if (pi == pj) continue;
					if (pi == -1) djs.addElement(p_map[id], pj);
					else mid = djs.link(pi, pj);
					int k = (*p_vcp)[nv-1-djs.getMaxNode(pj)].id;
#ifdef _DEBUG
					//printf("add arc (%d %d)\n", id, k);
#endif
					p_tree->addArc(id, k, nid);
					djs.setMaxNode(p_map[id], mid);
			}
		}
	}
	delete[] p_map;

	p_tree->setMaxID(getNVerts());
	return p_tree;
}