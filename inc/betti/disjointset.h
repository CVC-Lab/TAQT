#ifndef DISJOINT_SET_H
#define DISJOINT_SET_H

//#include <dynarray.h>

/*
#ifdef NO_HASHMAP
#include <map>
#else
#include <hash_map>
#endif

using namespace std;*/


class DisjointSet {
public:
	/**
	 * Create an empty disjoint set
	 */
	DisjointSet(int _size);

	~DisjointSet(void);

	void unionSet(int x , int y) {
		link(find(x) , find(y));
	}

	/**
	 * Link two sets into a new set.
	 * @return the canonical element of the new set
	 */
	int link(int x, int y) {
		if (nodes[x].rank > nodes[y].rank) {
			int t = x; x = y; y = t;
		} else if (nodes[x].rank == nodes[y].rank) {
			(nodes[y]).rank++;
		}
		nodes[x].p = y;
		return y;
	}
	/**
	 * Add an elment x to the set with parent p
	 */
	int addElement(int x, int p) {
		nodes[x].p = p;
		nodes[x].rank = (nodes[p].rank == 0)? 1:nodes[p].rank;

		return p;
	}
	/**
	 * Create a new set with element t
	 */
	void makeSet(int t){
		DJSNode node(t);
		nodes[t] = node;
	}

	/**
	 * @return -1 if t is not in the sets
	 */
	int find(int t) {
		//if (t >= (int)nodes.size() || t < 0) return -1;
		if (nodes[t].p < 0)	return -1;
		if (t != nodes[t].p) nodes[t].p = find(nodes[t].p);
		return nodes[t].p;
	}

	int getMaxNode(int p) {
		return nodes[p].max;
	}

	void setMaxNode(int x, int p) {
		nodes[p].max = x;
	}

private:
	struct DJSNode {
		DJSNode(int t = -1)
		: rank(0), p(t), max(t){
		}

		int rank;		// union heurist parameter
		int p;			// conanical member of the set	
		int max;		// the maximum node in the set
	};
	//dynamic_array<DJSNode> nodes;
	DJSNode *nodes;
	int size;
};

#endif



