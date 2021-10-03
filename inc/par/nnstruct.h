#ifndef _NEAREST_NEIGHBOR_H
#define _NEAREST_NEIGHBOR_H

#include <assert.h>
#include "geom.h"
#include "point3d.h"
#include "dynarray.h"

class NNeighbor {
public:

	bool operator < (const NNeighbor& pnt) {
		return dist < pnt.dist;
	}
	int id;
	float dist;
};

int compareNNeighbors(const void *p1, const void *p2);

/**
 * A base class for finding nearest neighbors of a point in 3D. 
 * Given a point cloud in 3D, it provides functions to find k neart neighboring points
 * of any given point.
 * @note This is the query interface. Concrete search algorithms are implemented in 
 *  subclasses.
 */
class NNStruct
{
public:
	NNStruct(int np, Point3D *points) : p_pnts(points), m_np(np) {
		assert(np >= 1);
	}
	
	virtual ~NNStruct() {}

	/**
	 * Search for k nearst neighbors.
	 * @return actual number of neighbors found.
	 * @note The neighboring points may include the query point itself. 
	 */
	virtual int NNSearch(const Point3D& pnt, int k, int nn_idx[], float dist[]) {
		return 0;
	}

	/*
	 * Search for neighbors within a sphere of a point.
	 * 
	 * @param r2: the square of radius r^2 
	 */
	virtual int NNSearch(const Point3D& pnt, float r2, dynamic_array<int>& nn_idx) {
		return 0;
	}

	/*
	 *	Search for neighbors within a sphere.
	 *  Up to k neastes neighbors will be returned.
	 */
	virtual int NNSearch(const Point3D& pnt, int k, float r2, int nn_idx[], float dist[]) {
		return 0;
	}

protected:
	int m_np;
	/**
	 * Array of all points in the 3D cloud.
	 *
	 * @note NNStruct doesn't have the ownership of this array.
	 *		 Its memory should not be released by NNStruct.
	 */
	Point3D *p_pnts;
};

#endif

