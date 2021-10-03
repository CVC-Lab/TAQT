#ifndef CONTOUR_TREE_H
#define CONTOUR_TREE_H
#include "node.h"

#include <vector>

#define NO_HASHMAP

#ifdef NO_HASHMAP
#include <map>
#else
#include <hash_map>
#endif

#include "basic.h"
#include "CriticalPoint.h"
#include "linkedlist.h"

using namespace std;

class AugmentedContourTree;

/**
 * ContourTree -- A class representing the contour tree structure.
 * @note Nodes are always added in increasing order
 */
class ContourTree {
public:
	/**
	 * Create an empty contour tree
	 */
	ContourTree();

	///
	virtual ~ContourTree();

	/**
	 * Add a node of id n into the contour tree.
     * @param n Node id
     * @param val function of the Node
	 */
	void addNode(int n, float val);

	void addNode(const CriticalPoint& cp);

	/**
	 * Remove a node with id n in the contour tree.
	 * @param n Node ID
	 */
	virtual void delNode(int n);

	/**
	 * Add an arc (n, m) into the contour tree.
	 * @note Node n and m should already be in the tree.

	 * @param x: (n, x) is the edge that cause the addition of this arc.
	 */
	void addArc(int n, int m, int x1 = -1, int x2 = -1);

	virtual void reduce(void);
	
	void print(void);
	
	/**
	 * Construct an Augment Contourtree from current contour tree.
	 * @note all nodes in the current tree are removed after augmentation.
	 */
	AugmentedContourTree* augment();

	static ContourTree* mergeTree(ContourTree* p_jtree, ContourTree* p_stree);
	
	/**
	 * Remove SuperArc (n, m)
	 */
	void removeArc(int n, int m);

	/**
	 * Check if a node in the contour tree is leaf.
	 * A node is a leaf if it isn't the root and has only one neighbor.
	 */
	bool isLeaf(int n);

	/**
	 * Get number of nodes.
	 */
	int getNumOfNodes() const { return (int)(nodes.size()); }

	/**
	 * Get the Node with id n in the contour tree.
	 * @return NULL if no node with id n in the tree
	 */
	SuperNode* getNode(int n);

	/**
	 * Get the sorted vector of all critical points in the contour tree
	 */
	vector<CriticalPoint>* getSortedPoints() {
		vector<CriticalPoint>* p_vcp = new vector<CriticalPoint>;
		if(!node_ids.empty()) {
			LL_Node<CriticalPoint>* ptr = node_ids.head();
			while(ptr != NULL) {
				p_vcp->push_back(ptr->object);
				ptr = ptr->next;
			}
		}
		return p_vcp;
	}

	/**
	 * Truncate all vertices lower than x;
	 */
	void truncateLower(float x);

	/**
	 * Write the contour tree to a file
	 */
	virtual void dump(const char* fname);

	/**
	 * Read contents of a contour tree from a file
	 */
	virtual void load(const char* fname);
	
	/**
	 * Get the range of function values
	 */
	void getMinMax(float& min, float &max) {
		LL_Node<CriticalPoint>* head = node_ids.head();
		LL_Node<CriticalPoint>* tail = node_ids.tail();

		min = head->object.val;
		max = tail->object.val;
	}
	
	/*
	 *	Set the maximum ID of the critical points
	 */
	void setMaxID(int _id) {
		max_id = _id;
	}

	int getMaxID() const {
		return max_id;
	}
	
//protected:
	typedef struct LeafNode {
		LeafNode(int _id, unsigned char _type) {
			id = _id;
			type = _type;
		}

		int id;
		unsigned char type;		// 0: join tree, 1: split tree
	} LeafNode;
#ifdef NO_HASHMAP
	typedef map<int, SuperNode*> supernode_map;
	typedef map<int, LL_Node<CriticalPoint>* > critpoint_map;
	typedef map<int, LL_Node<CriticalPoint>* >* critpoint_map_p;
#else
	typedef hash_map<int, SuperNode*> supernode_map;
	typedef hash_map<int, LL_Node<CriticalPoint>* > critpoint_map;
	typedef hash_map<int, LL_Node<CriticalPoint>* >* critpoint_map_p;
#endif

	supernode_map nodes;
	critpoint_map_p p_cpmap;
	LinkedList<CriticalPoint> node_ids;
	int max_id;
};

#endif


