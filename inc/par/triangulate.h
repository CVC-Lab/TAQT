#ifndef _PAR_TRIANGULATE_H
#define _PAR_TRIANGULATE_H

#include "dynarray.h"
#include "parsurface.h"

void sortPolarCC(int n, Polar2Coord* polars);

/**
 * 
 */
class PointConnect {
public:
	PointConnect(int nn = 6) {
		p_nbrs = new dynamic_array<int> (nn);
	}

	~PointConnect() {
		delete p_nbrs;
	}

	void insert(int n) {
		p_nbrs->insert(n);
	}
	 
	int numOfNeighbors() {
		return p_nbrs->length();
	}
	
	// counter-clock wise list of connected points
	dynamic_array<int> *p_nbrs;
};

/**
 * Triangulate a set of points sampling a surface.
 * Position and normal information is given for each point.
 */
class Triangulate {
public:
	Triangulate(PointSurface* psurf);
	
	~Triangulate();

	/**
	 * Triangulate the points based on the nearest neighbor information.
	 */
	void NNTriangulate();

	/**
	 *	
	 */
	void NNTriangulate(float r);

	void NNTriangulate2(float r);
	
	static void setNumOfNeighbors(int n) {
		nn = n;
	}

	int* getNeighbors(int i, int &nb) {
		nb = (*p_cnct)[i]->numOfNeighbors();
		return (*p_cnct)[i]->p_nbrs->data();
	}

	dynamic_array<PointConnect*>* p_cnct;
private:
	static int nn;	// number of numbering points 

	/**
	 * Point sample surface
	 *
	 * @note Triangulate doesn't have the ownership of this surface.
	 *		 Its memory should not be released by Triangulate.
	 */
	PointSurface* p_pntsrf;
};
#endif

