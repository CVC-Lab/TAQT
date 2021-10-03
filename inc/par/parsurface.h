#ifndef _PARTICLE_SURFACE_H
#define _PARTICLE_SURFACE_H

#include "function.h"
#include "dynarray.h"
#include "point3d.h"
#include "mtxlib.h"

typedef enum {
	SPHERE_FUN = 0,
	REG3D_FUN,
	PDB_FUN
} VolType;

/**
 * Point based isosurface of a scalar function. 	
 */

class PointSurface
{
public:
	/*
	 * Construct a 3D function for different vol types
	 */ 
	PointSurface(VolType vt);						// analytical functions
	PointSurface(const char* fname, VolType vt = REG3D_FUN);	// function from file

	~PointSurface();

	void makeSurface(float isoval, int nsample = -1);

	/*
	 *	Trace a point to the surface
	 */
	void tracePoint(Point3D &pnt);
		
	void upSampling(int vd1, int vd2, int vd3, int n, Point3D* samples);
	
	void repulse();

	int getNumOfPoints() const {
		return p_pnts->length();
	}

	Point3D* getPoints(int &np) {
		np = p_pnts->length();
		return p_pnts->data();
	}

	void getBoundingBox(float min[3], float max[3]) {
		p_fun->getBoundingBox(min, max);
	}
	
	void flipNormals() {
		int n = p_pnts->length();
		for(int i = 0; i < n; i++) {
			(*p_pnts)[i].normal[0] = - (*p_pnts)[i].normal[0];
			(*p_pnts)[i].normal[1] = - (*p_pnts)[i].normal[1];
			(*p_pnts)[i].normal[2] = - (*p_pnts)[i].normal[2];			
		}
	}

	/*
	 *	Add a sampling point.
	 *
	 */
	int insertPoint(float pos[3]);

	/*
	 *	Distance between point i and j
	 */
	float pointDistance(int i, int j) {
		return (float)sqrt(((*p_pnts)[i].pos[0] -(*p_pnts)[j].pos[0])*((*p_pnts)[i].pos[0] -(*p_pnts)[j].pos[0])
						  +((*p_pnts)[i].pos[1] -(*p_pnts)[j].pos[1])*((*p_pnts)[i].pos[1] -(*p_pnts)[j].pos[1])
						  +((*p_pnts)[i].pos[2] -(*p_pnts)[j].pos[2])*((*p_pnts)[i].pos[2] -(*p_pnts)[j].pos[2]));
	}

	void getFuncMinMax(float &min, float &max) {
		p_fun->getFuncMinMax(min, max);
	}

	float isoValue() const {
		return m_isoval;
	}

	bool hasPoints() const {
		return (p_pnts && p_pnts->length() > 0);
	}
	
protected:
	void init();
	
	vector3 repulsionForce(int p1, int p2, float rr1, float rr2);

	float	m_isoval;
	int		m_nsamp;
	dynamic_array<Point3D>* p_pnts;	
	Function * p_fun;
	float m_threshold;		// error threshold
};
#endif

