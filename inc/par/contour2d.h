//  ______________________________________________________________________
//
//    FILE
//      contour2d.h - Class for a 2d isocontour polyline
//
//      Copyright (c) 1998 Emilio Camahort, Dan Schikore
//
//    DESCRIPTION
//      contour2d is a class for representing a 2d isocontour polyline.
//  ______________________________________________________________________

// $Id: contour2d.h,v 1.1 2006/03/29 22:30:41 xyzhang Exp $

#ifndef _CONTOUR_2D_H
#define _CONTOUR_2D_H

#include <string.h>
#include <sys/types.h>

#include "dynarray.h"


class Contour2D {
public:
	typedef struct Edge {
		Edge(int _v1, int _v2) {
			v1 = _v1; v2 = _v2;
		}
		int v1, v2;
	} Edge;

	typedef struct Vertex {
		Vertex(float _x, float _y, float _z) {
			x = _x; y = _y; z = _z;
		}
		float x, y, z;
	} Vertex;

	// constructor
	Contour2D();

	// destructor
	~Contour2D();

	// reset (delete all vertices and triangles)
	void Reset(void);
	void Done(void);
	int  isDone(void) { return(done); }

	// add a vertex with the given position and normal
	int addVert(float p[2])
	{ return(addVert(p[0], p[1])); }
	int addVert(float, float);

	// add an edge indexed by the given 2 vertices
	int addEdge(int v[2])   { return(addEdge(v[0], v[1])); }
	int addEdge(int, int);

	// get the number of vertices or edges
	int getNVert(void)        { return(m_NVert); }
	int getNEdge(void)        { return(m_NEdge);  }

	// write vertices and triangles to a file
	int write(const char *filename);

	void setExtent(float min[3], float max[3])
	{
		memcpy(minext, min, sizeof(float[3]));
		memcpy(maxext, max, sizeof(float[3]));
	}

protected :

	int done;				// done with isocontour ??


	// the number of vertices and edges
	int m_NVert, m_NEdge;

	float minext[3], maxext[3];

	// arrays of vertices, and edges
	dynamic_array<Vertex>* p_Verts;
	dynamic_array<Edge>* p_Edges;
};

#endif
