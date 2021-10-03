#ifndef _SURFACE3D_H
#define _SURFACE3D_H

#include "point3d.h"
#include "dynarray.h"

class TriIndex {
public:
	TriIndex(int v1 = 0, int v2 = 0, int v3 = 0) {
		vid[0] = v1; vid[1] = v2; vid[2] = v3;
	}
	
	int vid[3];
};

/**
 * A polygonal surface consisting a collection of triangles. 
 * @note It is represented as a face index set
 * @author Xiaoyu Zhang
 * @version 1.1
 */
class Surface3D {
public:
    Surface3D(float _lower[3] = NULL, float _upper[3] = NULL);
	
    virtual ~Surface3D();
	
    int getNumOfVerts() const { return verts.length(); }
	
    int getNumOfTris() const {return tris.length(); }
	
    inline void getBoundingBox(float min[3], float max[3]) const;
	
    inline void setBoundingBox(float min[3], float max[3]);
	
    /**
	 * Insert a new vertex with given position and normal to the surface.
	 * @param c the value used to color the vertex. 
	 */
    int addVert(float pos[3], float norm[3], int c = 0);
	
    /**
	 * Insert a new triangle (v1, v2, v3) to the surface.
	 */
    int addTri(int v1, int v2, int v3);
	
	const Point3D* getVerts() { return verts.data(); }
	const TriIndex* getTris() { return tris.data(); }

	void getTriIndex(int nt, TriIndex& tid) {
		tid = tris[nt];
	}

	void getVertex(int nv, Point3D& pnt) {
		pnt = verts[nv];
	}

	void flipNormals() {
		int nv = getNumOfVerts();
		for(int i = 0; i < nv; i++) {
			verts[i].normal[0] = -verts[i].normal[0];
			verts[i].normal[1] = -verts[i].normal[1];
			verts[i].normal[2] = -verts[i].normal[2];
		}
	}

    /**
	* Reset (delete all vertices and triangles). 
	*/
    void reset();
	
protected:    
    // lower and upper bound of the surface geometry
    float   m_minext[3], m_maxext[3];
    bool    m_colored;            // If the surface colored using function

	dynamic_array<Point3D> verts;
	dynamic_array<TriIndex> tris;
};

void Surface3D::getBoundingBox(float min[3], float max[3]) const
{
	for(int i = 0; i < 3; i++) {
		min[i] = m_minext[i];
		max[i] = m_maxext[i];
	}
}

void Surface3D::setBoundingBox(float min[3], float max[3])
{
	for(int i = 0; i < 3; i++) {
		m_minext[i] = min[i];
		m_maxext[i] = max[i];
	}
}
#endif


