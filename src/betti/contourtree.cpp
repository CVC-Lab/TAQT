#include <stdio.h>
#include <assert.h>
#include "contourtree.h"
#include "actree.h"

#include <deque>
#include <iostream>
#include <fstream>
using namespace std;

ContourTree::ContourTree()
{
	p_cpmap = new critpoint_map;
}

ContourTree::~ContourTree()
{
	supernode_map::iterator it = nodes.begin();

	while(it != nodes.end()) {
		SuperNode* p_node = (*it).second;
		nodes.erase(it++);
		delete p_node;
	}
	if(p_cpmap != NULL) {
		delete p_cpmap;
	}
}

void ContourTree::addNode(int n, float val) 
{
	SuperNode *p_node = new SuperNode(n);
	nodes[n] = p_node;
	CriticalPoint cp(n, val);
	LL_Node<CriticalPoint>* p_ll = new LL_Node<CriticalPoint>(cp);
	(*p_cpmap)[n] = p_ll;
	node_ids.insert(p_ll);
}

void ContourTree::addNode(const CriticalPoint& cp)
{
	int n = cp.id;
	SuperNode *p_node = new SuperNode(n);
	nodes[n] = p_node;
	LL_Node<CriticalPoint>* p_ll = new LL_Node<CriticalPoint>(cp);
	(*p_cpmap)[n] = p_ll;
	node_ids.insert(p_ll);
}

void ContourTree::addArc(int n, int m, int x1, int x2)
{
	supernode_map::iterator in, im;

	in = nodes.find(n);
	im = nodes.find(m);
#ifdef _DEBUG
	assert(in != nodes.end() && im != nodes.end());
#endif
	((*in).second)->addNeighbor(m, x1);
	((*im).second)->addNeighbor(n, x2);
}

SuperNode* ContourTree::getNode(int n)
{
	supernode_map::iterator it = nodes.find(n);

	if(it == nodes.end()) return NULL;
	return (*it).second;
}

void ContourTree::reduce(void)
{
	LL_Node<CriticalPoint>* ptr = node_ids.head();
	while(ptr != NULL) {
		int id = ptr->object.id;
		SuperNode* p_node = (*(nodes.find(id))).second;
		ptr = ptr->next;
		if(p_node->degree() == 2) {
			delNode(id);
		}
	}
	
}

void ContourTree::print(void)
{
	int nd = (int)nodes.size();
	supernode_map::iterator it;

	LL_Node<CriticalPoint>* p_cur = node_ids.head();
	for(int i = 0; i < nd; i++) {
		int id = (p_cur->object).id;
		p_cur = p_cur->next;
		it = nodes.find(id);
		SuperNode* cnode = (*it).second;
		printf("%d: ", (*it).first);
		for(int j = 0; j < cnode->degree(); j++) {
			printf("%d ", cnode->getNeighbor(j));
		}
		printf("\n");
	}
}

void ContourTree::delNode(int n)
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
	} else if(p_supnode->degree() == 2) {
		int j = p_supnode->getNeighbor(0);
		int k = p_supnode->getNeighbor(1);
		supernode_map::iterator it = nodes.find(j);

		SuperNode* p_node1 = (*it).second;
		it = nodes.find(k);
		SuperNode* p_node2 = (*it).second;
		int x1 = p_node1->getNeighborTag(n);
		int x2 = p_node2->getNeighborTag(n);
		p_node1->removeNeighbor(n);
		p_node2->removeNeighbor(n);
		addArc(j, k, x1, x2);
	}
	delete p_supnode;
	nodes.erase(it_sup);
	p_cpmap->erase(it_lln);
	node_ids.remove(p_llnode);
}

bool ContourTree::isLeaf(int n)
{
	supernode_map::iterator it_sup = nodes.find(n);
	critpoint_map::iterator it_lln = p_cpmap->find(n);

	assert(it_sup != nodes.end() && it_lln != p_cpmap->end());
	SuperNode* p_supnode = (*it_sup).second;
	LL_Node<CriticalPoint>* p_llnode = (*it_lln).second;

	return (p_supnode->degree() == 1 && !node_ids.isTail(p_llnode));
}

ContourTree* ContourTree::mergeTree(ContourTree* p_jtree, ContourTree* p_stree)
{
	printf("start merging join and split trees...\n");
	ContourTree* p_tree = new ContourTree();

	deque<LeafNode> leafs;
	// iterate through the join tree
	LL_Node<CriticalPoint>* ptr = p_jtree->node_ids.head();

	while(ptr != NULL) {
		int id = ptr->object.id;
		p_tree->addNode(ptr->object);
		if(p_jtree->isLeaf(id) && p_stree->getNode(id)->degree() <= 2) {
			LeafNode leaf(id, 0);
			leafs.push_back(leaf);
		} else if(p_stree->isLeaf(id) && p_jtree->getNode(id)->degree() <= 2) {
			LeafNode leaf(id, 1);
			leafs.push_back(leaf);
		}
		ptr = ptr->next;
	}
	
	// remove nodes from join and split trees, and add superarcs to the merged tree.
	while(!leafs.empty()) {
		int id, nid;
		SuperNode* p_snode;
		LeafNode leaf = leafs.front();
		leafs.pop_front();
		id = leaf.id;
		switch(leaf.type) {
		case 0:
			p_snode = p_jtree->getNode(id);			
			break;
		case 1:
			p_snode = p_stree->getNode(id);
			break;
		}

		if(p_snode == NULL) continue;

		int degree = p_snode->degree();
		if(degree > 0) { // not the last node in the tree
			nid = p_snode->getNeighbor(0);
			int x1, x2;
			SuperNode *tmp;
			switch(leaf.type) {
			case 0:
				tmp = p_jtree->getNode(nid);
				x2 = tmp->getNeighborTag(id);
				break;
			case 1:
				tmp = p_stree->getNode(nid);
				x2 = tmp->getNeighborTag(id);
				break;
			}
			x1 = p_snode->getTag(0);
#ifdef _DEBUG
			//printf("add arc (%d %d) ", id, nid);
			//printf("%d, %d \n", leaf.type, p_stree->getNode(26552));
			//fflush(stdout);
#endif
			p_tree->addArc(id, nid, x1, x2);
		}
		p_jtree->delNode(id);
		p_stree->delNode(id);

		if(degree > 0) {
			switch(leaf.type) {
			case 0:
				if(p_jtree->isLeaf(nid) && p_stree->getNode(nid)->degree() <= 2) {
					LeafNode leaf(nid, 0);
					leafs.push_back(leaf);
				} else if(p_stree->isLeaf(nid) && p_jtree->getNode(nid)->degree() <=2) {
					LeafNode leaf(nid, 1);
					leafs.push_back(leaf);
				}
				break;
			case 1:
				if(p_stree->isLeaf(nid) && p_jtree->getNode(nid)->degree() <= 2) {
					LeafNode leaf(nid, 1);
					leafs.push_back(leaf);
				} else if(p_jtree->isLeaf(nid) && p_stree->getNode(nid)->degree() <=2) {
					LeafNode leaf(nid, 0);
					leafs.push_back(leaf);
				}
				break;
			}
		}
	}
		
	p_tree->setMaxID(MAX(p_jtree->getMaxID(), p_stree->getMaxID()));
	return p_tree;
}

// Remove SuperArc (n, m)
void ContourTree::removeArc(int n, int m)
{
	//todo
}

AugmentedContourTree* ContourTree::augment()
{
	AugmentedContourTree* p_actree = new AugmentedContourTree();

	int i, j, nv = (int)nodes.size();
	int *x = new int[nv];
	int *b = new int[nv];
	int *p_map = new int[nv];
	for(i = 0; i < nv; i++) {
		x[i] = 0;
		b[i] = 0;
	}

	deque<int> leafs;
	// iterate through the current tree
	i = 0;
	LL_Node<CriticalPoint>* ptr = node_ids.head();
	while(ptr != NULL) {
		int id = ptr->object.id;
		p_map[id] = i;
		p_actree->addNode(ptr->object);
		if(isLeaf(id)) {
			leafs.push_back(id);
		}
		ptr = ptr->next;
		i++;
	}

	while(!leafs.empty()) {
		int nid, id = leafs.front();
		leafs.pop_front();
		SuperNode* p_snode = getNode(id);
		int degree = p_snode->degree();
		if(degree > 0) { // not the last node
			nid = p_snode->getNeighbor(0);
			int x1 = p_snode->getTag(0);
			SuperNode *tmp = getNode(nid);
			int x2 = tmp->getNeighborTag(id);
			i = p_map[id]; 
			j = p_map[nid];
			//LL_Node<CriticalPoint>* p_llnode = (*(p_cpmap->find(nid))).second; 
			//int USj = p_llnode->object.US;
			//int LSj = p_llnode->object.LS;
			//int dbej = p_llnode->object.dbe;
			LL_Node<CriticalPoint>* p_llnode = (*(p_cpmap->find(id))).second;
			int LSi = p_llnode->object.LS;
			int USi = p_llnode->object.US;
			int dbei = p_llnode->object.dbe;
			
			int delta = (i < j)? 1:-1;
			int xe = delta * (x[i] - USi + LSi);
			int be = delta * (b[i] + dbei);
			x[j] += delta * xe;
			b[j] += delta * be;
#ifdef _DEBUG
            //printf("add acr (%d %d) with xe = %d, be = %d\n", id, nid, xe, be);
#endif
			p_actree->addArc(id, nid, xe, be, -1, x1, x2);
		}
		delNode(id);
		if(degree > 0) {
			if(isLeaf(nid)) leafs.push_back(nid);
		}
	}
	delete[] x;
	delete[] b;
	delete[] p_map;

	p_actree->setMaxID(getMaxID());
	return p_actree;
}

void ContourTree::dump(const char* fname)
{
	ofstream out(fname);
	out << 0 << endl;
	LL_Node<CriticalPoint>* ptr = node_ids.head();
	while(ptr != NULL) {
		out << ptr->object.id << " " << ptr->object.val << " ";
		SuperNode* p_node = (*(nodes.find(ptr->object.id))).second;
		out << p_node->degree() << " ";
		for(int i = 0; i < p_node->degree(); i++) {
			out << p_node->getNeighbor(i) << " ";
		}
		out << endl;
		ptr = ptr->next;
	}
	out.close();
}

void ContourTree::load(const char* fname)
{
	ifstream in(fname);
	int head;
	in >> head;
	assert(head == 0);
	while(!in.eof()) {
		int id, nid, degree;
		float val;
		in >> id;
		if(in.eof()) break;
		in >> val;
		in >> degree;
		
		addNode(id, val);
		for(int i = 0; i < degree; i++) {
			in >> nid;
			SuperNode* p_node = (*(nodes.find(id))).second;
			p_node->addNeighbor(nid);
		}
	}
	in.close();
}

void ContourTree::truncateLower(float x)
{
	// traverse every node in contour tree incremently, 
	// delete the node if its value and the values of all its neighbors 
	// are below x
	LL_Node<CriticalPoint>* ptr = node_ids.head();
	while(ptr != NULL) {
		int id = ptr->object.id;
		SuperNode* p_node = (*(nodes.find(id))).second;
		
		ptr = ptr->next;		
		critpoint_map::iterator it = p_cpmap->find(p_node->getID());
		float val = (*it).second->object.val;
		if(val >= x) break;
		bool truncFlag = true;
		for(int i = 0; i < p_node->degree(); i++) {
			SuperNode* p_nei = (*(nodes.find(p_node->getNeighbor(i)))).second;
			it = p_cpmap->find(p_nei->getID());
			if ((*it).second->object.val > x) {
				truncFlag = false;
				break;
			}
		}
		if(truncFlag && p_node->degree() <= 2) {
			delNode(id);
		}
		

	}
}
