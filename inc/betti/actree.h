#ifndef AUGMENTED_CONTOURTREE_H
#define AUGMENTED_CONTOURTREE_H

#include "arc.h"
#include "contourtree.h"
#include "linkedlist.h"
#include "interval_tree.h"

#include <dynarray.h>
#include "volume3dcritical.h"

class CutVertex;

/**
 *	Augmented Contour Tree with Betti Numbers
 */
class AugmentedContourTree : public ContourTree {
public:
	/**
	 * Construct an Augmented Contourtree 
	 */
	AugmentedContourTree();

	virtual ~AugmentedContourTree();

	void calcBettiNumbers();

	virtual void reduce();

	virtual void delNode(int n);

	/*
	 *	
	 */
	void addArc(int v1, int v2, int xe, int be, int intra = -1, int x1 = -1, int x2 = -1);

	/**
	 * Finailize the tree and construct search structures
	 */
	void done();

	/**
	 * Write the contour tree to a file
	 */
	virtual void dump(const char* fname);

	/**
	 * Read contents of a contour tree from a file
	 */
	virtual void load(const char* fname);
	
	/**
	 * Compute the BettiNumbers for the level set of isovalue x
	 */
	BettiNumber getBettiNumbers(float x);
	
	/**
	 * Find the intersection points with a range cut
	 */
	CutVertex* cut(VolumeReg3Critical* p_data, float x, int& nv);
	

	friend class DualGraph;

    LinkedList<SuperArc> arcs;
	IntervalTree<LL_Node<SuperArc>* >* p_intree;
protected:
	float f_min(float x, float y) {
		return (x < y)? x:y;
	}
	float f_max(float x, float y) {
		return (x > y)? x:y;
	}
};
#endif


