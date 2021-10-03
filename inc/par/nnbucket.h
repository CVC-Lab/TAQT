#ifndef _NN_BUCKET_H
#define _NN_BUCKET_H

#include "nnstruct.h"
#include "dynarray.h"

/**
 * Use 3D buckets to help nearest neighbor search.
 *
 */
class NNBucket : public NNStruct
{
public:
	NNBucket(int np, Point3D *points, int buck_dim[3] = 0);
		
	virtual ~NNBucket();

	virtual int NNSearch(const Point3D& pnt, int k, int nn_idx[], float dist[]);

	virtual int NNSearch(const Point3D& pnt, float r2, dynamic_array<int>& nn_idx);

	virtual int NNSearch(const Point3D& pnt, int n, float r2, int nn_idx[], float dist[]);

	/*
	 *	Insert the point with pid in the point array to the buckets and replace
	 *  the existing array with the new point array.
	 */
	void insertAndReplace(int np, Point3D* _points, int pid);
	
protected:
	void init();

	int index2buck(int idx[3]) {
		return (idx[0] + idx[1]*(m_dim[0]-1) + idx[2]*(m_dim[0]-1)*(m_dim[1]-1)); 
	}

	int index2buck(int i, int j, int k) {
		return (i + j*(m_dim[0]-1) + k*(m_dim[0]-1)*(m_dim[1]-1)); 
	}

	bool valid(int idx[3]) {
		return (idx[0] >= 0 && idx[0] < (m_dim[0]-1) &&
				idx[1] >= 0 && idx[1] < (m_dim[1]-1) &&
				idx[2] >= 0 && idx[2] < (m_dim[2]-1));
	}

	int m_dim[3];
	float m_min[3], m_max[3];
	dynamic_array<int> *p_bkidx;
};
#endif

