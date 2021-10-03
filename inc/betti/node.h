#ifndef CONTOUR_NODE_H
#define CONTOUR_NODE_H

#include <assert.h>

#include <map>
#include <vector>
#include "arc.h"
#include "linkedlist.h"

using namespace std;

/**
 * SuperNode -- A class representing a node in a contour tree
 */
class SuperNode {
public:
	/**
	 * Construct a node with id i.
	 */
	SuperNode(int i);

	///
	~SuperNode() {}

	/**
	 * Get the number of neighboring nodes.
	 */
	int degree() const { return (int)(neighbors.size()); }

	/**
	 * Add a neighor of id n.
	 * @param x: edge (n, x) causes the add of the neighbor
	 */
	void addNeighbor(int n, int x = -1);

	/**
	 * Remove the neighbor with id n
	 */
	void removeNeighbor(int n);

	/**
	 * Get the nth neighbor.
	 * @note n < node's degree
	 */
	int getNeighbor(int n);

	int getTag(int n);
	
	int getNeighborTag(int nid);

	/**
	 * Get the ID of a supernode.
	 */
	int getID() const { 
		return id;
	}

	/**
	 */
	void addArc(int n, LL_Node<SuperArc>* p_arc) {
		arcs[n] = p_arc;
	}

	/**
	 * Return the pointer to the superarc between this node and node n
	 */
	LL_Node<SuperArc>* getArcPtr(int n) {
		map<int, LL_Node<SuperArc>*>::iterator it = arcs.find(n);
		assert(it != arcs.end());
		return (*it).second;
	}

	/**
	 */
	void removeArc(int n) {
		map<int, LL_Node<SuperArc>*>::iterator it = arcs.find(n);
		assert(it != arcs.end());
		arcs.erase(it);
	}

private:
	struct sn_pair {
		sn_pair(int n, int x) : nei(n), tag(x) {}
		
		int nei;
		int tag;
	};
	int id;
	vector<sn_pair> neighbors;
	map<int, LL_Node<SuperArc>*> arcs;	
};

#endif


