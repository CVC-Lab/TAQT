#ifndef _NN_KDTREE_H
#define _NN_KDTREE_H

#include "nnstruct.h"

#include <ANN/ANNx.h>

class NNKDTree : public NNStruct {
public:
	NNKDTree(int np, Point3D *points);

	~NNKDTree();

	virtual int NNSearch(const Point3D &pnt, int k, int nn_idx[], float dist[]);

	/*
	 *	Search for neighbors within a sphere.
	 *  Up to k neastes neighbors will be returned.
	 */
	virtual int NNSearch(const Point3D& pnt, int k, float r2, int nn_idx[], float dist[]);

private:
	ANNkd_tree *p_tree;
	ANNpointArray m_pntarray;
};

#endif

