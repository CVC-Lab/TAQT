#ifndef RANGE_CUT_H
#define RANGE_CUT_H

#include "actree.h"
#include "linkedlist.h"

class DualNode;

/**
 * An edge in a subrange
 */
class CutEdge {
public:
	CutEdge(int _vl, int _vu) {
		vl = _vl;
		vu = _vu;
		flag_l = flag_u = false;
	}
	
	int 	vl, vu;			// lower and upper vertices of the edge.
	bool    flag_l, flag_u;	// if vl or vu is a CutVertex.	
};

/**
 * A vertex of the contour tree cut by the functional ranges
 */
class CutVertex {
public:
	
	DualNode *ngb_l, *ngb_u;				// lower and upper dual node of the vertex
	int level, id;			
	LL_Node<SuperArc>* super_arc;			// the super arc that contains the cut vertex
};


#endif

