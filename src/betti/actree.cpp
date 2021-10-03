#include "cut.h"
#include "actree.h"
#include <fstream>
#include <iostream>

using namespace std;

AugmentedContourTree::AugmentedContourTree()
: p_intree(NULL) {

}

AugmentedContourTree::~AugmentedContourTree()
{
	if(p_intree != NULL) delete p_intree;
}

void AugmentedContourTree::calcBettiNumbers()
{
	LL_Node<SuperArc>* ptr = arcs.head();
	while(ptr != NULL) {
		ptr->object.calcBettiNumber();
		ptr = ptr->next;
	}
}

void AugmentedContourTree::reduce()
{
	LL_Node<CriticalPoint>* ptr = node_ids.head();
	while(ptr != NULL) {
		int id = ptr->object.id;
		SuperNode* p_node = (*(nodes.find(id))).second;
		int LS = ptr->object.LS;
		int US = ptr->object.US;
		int dbe = ptr->object.dbe;
		ptr = ptr->next;
		if(p_node->degree() == 2 && LS == 0 && US == 0) {
			delNode(id);
		}
	}
}

void AugmentedContourTree::delNode(int n)
{
	supernode_map::iterator it_sup = nodes.find(n);
	critpoint_map::iterator it_lln = p_cpmap->find(n);

	assert(it_sup != nodes.end() && it_lln != p_cpmap->end());
	SuperNode* p_supnode = (*it_sup).second;
	LL_Node<CriticalPoint>* p_llnode = (*it_lln).second;
	assert(p_supnode->degree() <= 2);
	if(p_supnode->degree() == 1) {
		int nid = p_supnode->getNeighbor(0);
		supernode_map::iterator it = nodes.find(nid);
		SuperNode* p_snode = (*it).second;
		p_snode->removeNeighbor(n);
		p_snode->removeArc(n);				// new
		LL_Node<SuperArc>* p_arc = p_supnode->getArcPtr(nid);
		arcs.remove(p_arc);
	} else if(p_supnode->degree() == 2) {
		int j = p_supnode->getNeighbor(0);
		int k = p_supnode->getNeighbor(1);
		int xe, be;
		LL_Node<SuperArc>* p_arc1 = p_supnode->getArcPtr(j);
		LL_Node<SuperArc>* p_arc2 = p_supnode->getArcPtr(k);
		if(p_arc1->object.betti1() != p_arc2->object.betti1() || p_arc1->object.betti2() != p_arc2->object.betti2()) {
			// cannot delete the node
			return;
		}
		xe = p_arc1->object.xe;
		be = p_arc1->object.be;

		supernode_map::iterator it = nodes.find(j);
		SuperNode* p_node1 = (*it).second;
		it = nodes.find(k);
		SuperNode* p_node2 = (*it).second;

		arcs.remove(p_arc1);
		arcs.remove(p_arc2);

		p_node1->removeNeighbor(n);
		p_node1->removeArc(n);
		p_node2->removeNeighbor(n);
		p_node2->removeArc(n);
		addArc(j, k, xe, be, n);
	}
	delete p_supnode;
	nodes.erase(it_sup);
	p_cpmap->erase(it_lln);
	node_ids.remove(p_llnode);
}

void AugmentedContourTree::addArc(int v1, int v2, int xe, int be, int nid, int x1, int x2)
{
	supernode_map::iterator in, im;

	in = nodes.find(v1);
	im = nodes.find(v2);
	assert(in != nodes.end() && im != nodes.end());
	((*in).second)->addNeighbor(v2, x1);
	((*im).second)->addNeighbor(v1, x2);

	SuperArc sa(v1, v2, xe, be);
	if(nid >= 0) sa.addIntraNode(nid);
	LL_Node<SuperArc>* p_arcnode = new LL_Node<SuperArc>(sa);

	arcs.insert(p_arcnode);
	((*in).second)->addArc(v2, p_arcnode);
	((*im).second)->addArc(v1, p_arcnode);
}

void AugmentedContourTree::dump(const char* fname)
{
	ofstream out(fname);
	out << 1 << endl;
	int nc = (int)nodes.size(), na = (int)arcs.size();
	out << nc << " " << na << endl;
	LL_Node<CriticalPoint>* ptr = node_ids.head();
	while(ptr != NULL) {
		out << ptr->object.id << " " << ptr->object.val << " ";
		SuperNode* p_node = (*(nodes.find(ptr->object.id))).second;

		out << p_node->degree() << " ";

		for(int i = 0; i < p_node->degree(); i++) {

			out << p_node->getNeighbor(i) << " " << p_node->getTag(i) << " ";

		}

		out << endl;

		ptr = ptr->next;

	}
	// output arcs
	LL_Node<SuperArc>* arcptr = arcs.head();
	while(arcptr != NULL) {
		out << arcptr->object.v1 << " " << arcptr->object.v2 << " ";
		out << arcptr->object.xe << " " << arcptr->object.be << "  ";
		out << 1 << " " << arcptr->object.intra_node << endl;
		arcptr = arcptr->next;
	}
	out.close();
	
}

void AugmentedContourTree::load(const char* fname)
{
	int nc, na, i;
		
	ifstream in(fname);
	int head;
	in >> head;
	if (head != 1) {
		cerr << "not an augmented contour tree, aborted\n";
		in.close();
		return;
	}
	in >> nc; in >> na;
	for(i = 0; i < nc; i++) {
		int id, nid, xid, degree;
		float val;
		in >> id;
		in >> val;
		in >> degree;
		addNode(id, val);
		for(int j = 0; j < degree; j++) {
			in >> nid;

			in >> xid;

			SuperNode* p_node = (*(nodes.find(id))).second;
			p_node->addNeighbor(nid, xid);
		}
	}

	for(i = 0; i < na; i++) {
		int v1, v2, xe, be, count, intra_node;
		in >> v1;
		in >> v2; 
		in >> xe;
		in >> be;
		in >> count;
		in >> intra_node;
		addArc(v1, v2, xe, be, intra_node);
	}
	in.close();
	done();
}

void AugmentedContourTree::done()
{
	int i, j, nv = (int)nodes.size();
	float *vals = new float[nv];

	// traverse the critical point list
	LL_Node<CriticalPoint>* ptr = node_ids.head();
	i = 0;
	while(ptr != NULL) {
		vals[i] = (ptr->object).val;
		ptr = ptr->next;
		i++;
	}
	// compacting critical vals
	for(i = 1, j = 0; i < nv; i++) {
		if(vals[i] != vals[j]) {
			vals[++j] = vals[i];
		} 
	}

	if(p_intree != NULL) delete p_intree;
	p_intree = new IntervalTree<LL_Node<SuperArc>* >(j+1, vals);
	delete[] vals;
	
	// add arcs to the interval tree
	LL_Node<SuperArc>* arc_ptr = arcs.head();
	while(arc_ptr != NULL) {
		int n1 = arc_ptr->object.v1;
		int n2 = arc_ptr->object.v2;
		float x1 = (*(p_cpmap->find(n1))).second->object.val;
		float x2 = (*(p_cpmap->find(n2))).second->object.val;
		p_intree->insert(arc_ptr, f_min(x1, x2), f_max(x1, x2));
		arc_ptr = arc_ptr->next;
	}
	p_intree->done();
		
	// set IDs the super arcs 
	i = 0;
	arc_ptr = arcs.head();
	while (arc_ptr != NULL) {
		arc_ptr->object.id = i;
		i++;
		arc_ptr = arc_ptr->next;
	} 	
}

BettiNumber AugmentedContourTree::getBettiNumbers(float x)
{
	BettiNumber betti;

	if(p_intree == NULL) {
		cerr << "construct interval tree before querying Betti Numbers\n";
		return betti;
	}
	vector<LL_Node<SuperArc> *>* p_results = p_intree->query(x);
	vector<LL_Node<SuperArc> *>::iterator it;
	for(it = p_results->begin(); it != p_results->end(); ++it) {
		LL_Node<SuperArc>* ptr = (*it);
		betti.b0 += 1;
		betti.b1 += ptr->object.betti1();
		betti.b2 += ptr->object.betti2();
	}
	delete p_results;
	return betti;
}

CutVertex* AugmentedContourTree::cut(VolumeReg3Critical* p_data, float x, int& nv)
{
	const double minimal = 1e-8;
	if(p_intree == NULL) {
		cerr << "need to construct the interval tree first\n";
		return NULL;
	}
	int i = 0;

	vector<LL_Node<SuperArc> *>* p_results = p_intree->query(x);
	nv = p_results->size();
	CutVertex* cutverts = new CutVertex[nv];
	vector<LL_Node<SuperArc> *>::iterator it;
	int count = 0;
	for(it = p_results->begin(), i = 0; it != p_results->end(); ++it, i++) {
		LL_Node<SuperArc>* ptr = (*it);
				
		int n1 = ptr->object.v1;
		int n2 = ptr->object.v2;
		float x1 = p_data->getValue(n1);
		float x2 = p_data->getValue(n2);
		if(x1 == x2) {
			assert(x1 == x && x2 == x);
			//cutverts[count++].super_arc = ptr;
		} else if (x1 == x) {
			if (x2 > x) cutverts[count++].super_arc = ptr;
		} else if (x2 == x) {
			if (x1 > x) cutverts[count++].super_arc = ptr;
		} else {
			cutverts[count++].super_arc = ptr;
		}
		//printf("v1: %d %f, v2: %d %f, cut_val: %f\n", n1, x1, n2, x2, x);
	}
	nv = count;
	delete p_results;
	
	return cutverts;
}
