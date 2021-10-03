#ifndef MULTI_CONTREE_H
#define MULTI_CONTREE_H

#include "dual.h"

/*
 *	Multi-resolution contour tree
 */
class MultiConTree
{
public:
	MultiConTree(DualGraph *dtree, float total_vol = 1.0f);
	
	~MultiConTree();
	
	/**
	 *	Compute the similarity between two multi-resolution contour trees.
	 */
	float matching(MultiConTree* multree);
	
	/**
	 *	Maximum level of the MCT
	 */
	int maxLevel() const {
		return max_level;
	}

	/**
	 * get the total volume of the tree
	 */
	float totalVolume() const {
		float vol = 0;
		for (int i = 0; i < dual_trees[0]->nNodes(0); i++) {
			vol += dual_trees[0]->getNode(0, i)->volume;
		}
		return vol;
	}

	/**
	 * Normalize the volume such that the total is one
	 */
	 void normalize();
	 
	/**
	 * Prune small nodes of volume < threshold
	 */
	 void prune(float threshold); 
	 
	/**
	 * Get the center of mass and principal axes of nodes at level l and range k
	 */
	void getOrientation(int l, int k, float *ctr, float **axes);
	
protected:
	int max_range, max_level;
	DualGraph** dual_trees;
	float total_volume;
private:
	int log2(int);

	/**
	 *	Search for the best match for a dual node
	 */
	DualNode* matchSearch(DualNode* dualnode, DualGraph* graph);

	/**
	 *	The matching score between two nodes
	 */
	float matchScore(DualNode *node1, DualNode *node2);
};
#endif
