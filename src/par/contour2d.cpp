//  ______________________________________________________________________
//
//    NAME
//      Contour2D - Class for a 2d contour curve
//
//      Copyright (c) 1998 Emilio Camahort, Dan Schikore
//
//    SYNOPSIS
//      #include <contour2d.h>
//  ______________________________________________________________________

// $Id: contour2d.cpp,v 1.1 2006/03/29 22:30:41 xyzhang Exp $

#include <stdio.h>
#ifndef MACOS_X
#include <malloc.h>
#endif
#include <stdlib.h>

#include "contour2d.h"

//------------------------------------------------------------------------
//
// Contour2D() - basic constructor
//
//------------------------------------------------------------------------
Contour2D::Contour2D()
{
	done = 0;
	m_NVert = 0;
	m_NEdge = 0;
	p_Verts = new dynamic_array<Vertex>;
	p_Edges = new dynamic_array<Edge>;
}

//------------------------------------------------------------------------
//
// ~Contour2D() - free allocated memory
//
//------------------------------------------------------------------------
Contour2D::~Contour2D()
{
	delete p_Verts;
	delete p_Edges;
}

//------------------------------------------------------------------------
//
// addVert() - add a vertex with the given (unit) normal
//
//------------------------------------------------------------------------
int Contour2D::addVert(float x, float y)
{
	int n = m_NVert++;

	Vertex vert(x, y, maxext[2]);
	p_Verts->insert(vert);

	return(n);
}

//------------------------------------------------------------------------
//
// AddEdge() - add an edge indexed by its 2 vertices
//
//------------------------------------------------------------------------
int Contour2D::addEdge(int v1, int v2)
{
	int n = m_NEdge++;

	Edge edge(v1, v2);
	p_Edges->insert(edge);

	return(n);
}

//------------------------------------------------------------------------
//
// Reset() - clear vertex and edge info
//
//------------------------------------------------------------------------
void
Contour2D::Reset(void)
{
	m_NVert = 0;
	m_NEdge = 0;
	done = 0;
	p_Verts->clear();
	p_Edges->clear();
}

void Contour2D::Done(void)
{
	done = 1;
}

//------------------------------------------------------------------------
//
// write() - write vertex and triangles to a file
//
//------------------------------------------------------------------------

int	Contour2D::write(const char *filename)
{
	FILE *fp;
	int v, t;

	fp = fopen(filename, "w");

	// silent failure --> changed by Emilio: return 1 = ERROR
	if (fp == NULL)
		return 1;

	//fprintf(fp, "%d %d 0 0 0 0 0\n0 0 0\n", m_NVert, m_NEdge);
	fprintf(fp, "%d %d \n", m_NVert, m_NEdge);
	for (v = 0; v < m_NVert; v++) {
		Vertex vert = (*p_Verts)[v];
		fprintf(fp, "%g %g %g\n", vert.x, vert.y, vert.z);
	}
	//fprintf(fp, "0 0\n");

	for (t = 0; t < m_NEdge; t++) {
		Edge edge = (*p_Edges)[t];
		fprintf(fp, "%d %d\n", edge.v1, edge.v2);
	}
	fclose(fp);

	return 0;
}

