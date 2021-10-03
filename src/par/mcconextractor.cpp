#include <stdlib.h>
#include <math.h>

#include "cubes.h"
#include "vtkMarchingCubesCases.h"
#include "vecmath.h"
#include "mcconextractor.h"

/**
* Information to specify an edge.
*/
typedef struct {
	int dir;
	int di,dj,dk;
	int d1,d2;
} EdgeInfo;

static EdgeInfo edgeinfo[12] = {
	{ 0, 0, 0, 0, 0, 1 },
	{ 2, 1, 0, 0, 1, 2 },
	{ 0, 0, 0, 1, 3, 2 },
	{ 2, 0, 0, 0, 0, 3 },
	{ 0, 0, 1, 0, 4, 5 },
	{ 2, 1, 1, 0, 5, 6 },
	{ 0, 0, 1, 1, 7, 6 },
	{ 2, 0, 1, 0, 4, 7 },
	{ 1, 0, 0, 0, 0, 4 },
	{ 1, 1, 0, 0, 1, 5 },
	{ 1, 0, 0, 1, 3, 7 },
	{ 1, 1, 0, 1, 2, 6 }
};


MCContourExtractor::MCContourExtractor()
{
}

int MCContourExtractor::interpRegEdges(float val[8], float grad[8][3], float isoval, 
									   int i, int j, int k, int edge, Reg3Data* p_reg3data)
{
	float pt[3];
	float norm[3];
	int   clr = 0;
	EdgeInfo *ei = &edgeinfo[edge];
	int v = 0;
	float ival;			// intersection point
	float orig[3], span[3];

    EdgeIndex eid(i+ei->di, j+ei->dj, k+ei->dk, ei->dir);
    map<EdgeIndex, int, LTEdgeIndex>::iterator it = edgeset->find(eid);
    if(it != edgeset->end()) {
        return (*it).second;
    } else {
        pair<EdgeIndex, int> newedge(eid, nvert);
		edgeset->insert(newedge);
        nvert ++;
    }

	p_reg3data->getOrig(orig);
	p_reg3data->getSpan(span);
	// do the interpolation
	if(val[ei->d1] == val[ei->d2]) ival = 0.5f;
	else ival = (isoval - val[ei->d1])/(val[ei->d2] - val[ei->d1]);
	pt[0] = orig[0] + (i + ei->di)*span[0];
	pt[1] = orig[1] + (j + ei->dj)*span[1];
	pt[2] = orig[2] + (k + ei->dk)*span[2];
	norm[0] = grad[ei->d1][0]*(1.0f-ival) + grad[ei->d2][0]*ival;
	norm[1] = grad[ei->d1][1]*(1.0f-ival) + grad[ei->d2][1]*ival;
	norm[2] = grad[ei->d1][2]*(1.0f-ival) + grad[ei->d2][2]*ival;

	switch (ei->dir) {
	case 0:
		pt[0] += ival * span[0];
		break;
	case 1:
		pt[1] += ival * span[1];
		break;
	case 2:
		pt[2] += ival * span[2];
		break;
	}

	// add to the coordinate list
	normalize3f(norm);
	v = surf->addVert(pt, norm, clr);
	return(v-1);
}


Surface3D * MCContourExtractor::contourReg3Data(Reg3Data* p_reg3data, float isoval)
{
	int dim[3];
	int edge_v[12];
	float val[8], grad[8][3];
	int v[3];
	int i, j, k;
	int e, t;
	int code;
	int edge;

	surf = new Surface3D();
    edgeset = new map<EdgeIndex, int, LTEdgeIndex>();
    nvert = 0;
	int dx, dy, dz;
	p_reg3data->getDim(dim);
	dx = dim[0]-1;
	dy = dim[1]-1;
	dz = dim[2]-1;
	for(k = 0; k < dz; k++) {
		for(j = 0; j < dy; j++) {
			for(i = 0; i < dx; i++) {
				p_reg3data->getCellValues(i, j, k, val);

				code = 0;
				if (val[0] < isoval) code |= 0x01;
				if (val[1] < isoval) code |= 0x02;
				if (val[2] < isoval) code |= 0x04;
				if (val[3] < isoval) code |= 0x08;
				if (val[4] < isoval) code |= 0x10;
				if (val[5] < isoval) code |= 0x20;
				if (val[6] < isoval) code |= 0x40;
				if (val[7] < isoval) code |= 0x80;

				if (cubeedges[code][0] == 0) continue;
				p_reg3data->getCellGrads(i, j, k, grad);

				for (e=0; e<cubeedges[code][0]; e++) {
					edge = cubeedges[code][1+e];
					edge_v[edge] = interpRegEdges(val, grad, isoval, i, j, k, edge, p_reg3data);
				}

				for (t=0; triCases[code].edges[t] != -1; ) {
					v[0] = edge_v[triCases[code].edges[t++]];
					v[1] = edge_v[triCases[code].edges[t++]];
					v[2] = edge_v[triCases[code].edges[t++]];
					surf->addTri(v[0], v[1], v[2]);
				}
			}
		}
	}
    delete edgeset;
	return surf;
}


