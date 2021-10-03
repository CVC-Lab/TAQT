#include <math.h>
#include <cstdio>
#include <assert.h>
#include <cstdlib>
#include <cstring>

#include "volume3dcritical.h"
#include "CriticalPoint.h"
#include "disjointset.h"
#include <nrutil.h>

#include <algorithm>
#include <deque>

using namespace std;

const int VolumeReg3Critical::vert_neighboring[15][15] = {
	{7, 1, 3, 4, 8, 10, 12, 14},
	{7, 0, 2, 5, 9, 10, 12, 14},
	{7, 1, 3, 6, 9, 10, 13, 14},
	{7, 0, 2, 7, 8, 10, 13, 14},
	{7, 0, 5, 7, 8, 11, 12, 14},
	{7, 1, 4, 6, 9, 11, 12, 14},
	{7, 2, 5, 7, 9, 11, 13, 14},
	{7, 3, 4, 6, 8, 11, 13, 14},
	{5, 0, 3, 4, 7, 14},
	{5, 1, 2, 5, 6, 14},
	{5, 0, 1, 2, 3, 14},
	{5, 4, 5, 6, 7, 14},
	{5, 0, 1, 4, 5, 14},
	{5, 2, 3, 6, 7, 14},
	{14, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13}
};

// six simplex decomposition of a regular cell
const int VolumeReg3Critical::six_simplex_neighbors[14][3] = {
	{0, 0, 1}, {0, 1, 0}, {1, 0, 0},
	{0, 0, -1}, {0, -1, 0}, {-1, 0, 0},
	{0, 1, 1}, {1, 0, 1}, {1, 1, 0},
	{0, -1, -1}, {-1, 0, -1}, {-1, -1, 0},
	{1, 1, 1}, {-1, -1, -1}
};

const int VolumeReg3Critical::six_tetra_index[6][4][3] = {
	{{0, 0, 0}, {1, 0, 0}, {1, 0, 1}, {1, 1, 1}},
	{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {1, 1, 1}},
	{{0, 0, 0}, {0, 1, 0}, {1, 1, 0}, {1, 1, 1}},
	{{0, 0, 0}, {0, 1, 0}, {0, 1, 1}, {1, 1, 1}},
	{{0, 0, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}},
	{{0, 0, 0}, {0, 0, 1}, {0, 1, 1}, {1, 1, 1}}
};

const int VolumeReg3Critical::six_tetra_verts[6][4] = {
	{0, 1, 5, 7},
	{0, 1, 3, 7},
	{0, 2, 3, 7},
	{0, 2, 6, 7},
	{0, 4, 5, 7},
	{0, 4, 6, 7}
};

const int VolumeReg3Critical::six_tetra_neighbors[6][4][4] = {
	{{1, 0, 0, 5}, {0, 0, 0, 4}, {0, 0, 0, 1}, {0, -1, 0, 2}},
	{{1, 0, 0, 3}, {0, 0, 0, 2}, {0, 0, 0, 0}, {0, 0, -1, 4}},
	{{0, 1, 0, 0}, {0, 0, 0, 1}, {0, 0, 0, 3}, {0, 0, -1, 5}},
	{{0, 1, 0, 4}, {0, 0, 0, 5}, {0, 0, 0, 2}, {-1, 0, 0, 1}},
	{{0, 0, 1, 1}, {0, 0, 0, 0}, {0, 0, 0, 5}, {0, -1, 0, 3}},
	{{0, 0, 1, 2}, {0, 0, 0, 3}, {0, 0, 0, 4}, {-1, 0, 0, 0}}
};

int VolumeReg3Critical::g_dim[3];

static int compareInt(const void* p1, const void* p2)
{
	int x1 = *((int *)p1);
	int x2 = *((int *)p2);

	if (x1 < x2) return -1;
	if (x1 > x2) return 1;
	return 0;
}

bool VolumeReg3Critical::isMin(int id, float vals[8]) {
	if (id >= 8) return false;
	CriticalPoint pnt;
	pnt.id = id;
	pnt.val = vals[id];
	for (int i = 1; i <= 3; i++) {
		CriticalPoint p2(vert_neighboring[id][i], vals[vert_neighboring[id][i]]);
		if (p2 < pnt) return false;
	}
	return true;
}

VolumeReg3Critical::VolumeReg3Critical(Reg3Data * _fun, Reg3Data * _fun2) 
: p_fun(_fun), p_pot(_fun2)
{
	orig[0] = orig[1] = orig[2] = 0;
	p_fun->getDim(dim);
	p_fun->getOrig(m_orig);
	p_fun->getSpan(m_span);
	for (int i = 0; i < 3; i++) {
		g_dim[i] = dim[i];
	}
	p_vcp = NULL;
}

VolumeReg3Critical::VolumeReg3Critical(Reg3Data* _fun, int _dim[3], int _orig[3], Reg3Data *_fun2)
: p_fun(_fun), p_pot(_fun2)
{
	for (int i = 0; i < 3; i++) {
		dim[i] = _dim[i];
		orig[i] = _orig[1];
	}
	p_fun->getOrig(m_orig);
	p_fun->getSpan(m_span);
	p_vcp = NULL;
}

VolumeReg3Critical::~VolumeReg3Critical(void) 
{
}


void VolumeReg3Critical::split(VolumeReg3Critical * & p_mesh1, VolumeReg3Critical * & p_mesh2) {
	if (isCell()) {
		assert(0);		// never split a cell
	}
	char axis = 'x';
	// split along the longest dimension
	if (dim[0] >= dim[1]) {
		if (dim[0] >= dim[2]) axis = 'x';
		else axis = 'z';
	} else {
		if (dim[1] >= dim[2]) axis = 'y';
		else axis = 'z';
	}

	int d1[3], d2[3], o2[3];
	switch (axis) {
	case 'x':
		d1[0] = (dim[0] >> 1) + 1;
		d2[0] = ((dim[0]-1) >> 1) + 1;
		d1[1] = d2[1] = dim[1];
		d1[2] = d2[2] = dim[2];
		o2[0] = orig[0] + (d1[0]-1);
		o2[1] = orig[1];
		o2[2] = orig[2];
		break;
	case 'y':
		d1[0] = d2[0] = dim[0];
		d1[1] = (dim[1] >> 1) + 1;
		d2[1] = ((dim[1]-1) >> 1) + 1;
		d1[2] = d2[2] = dim[2];
		o2[0] = orig[0];
		o2[1] = orig[1] + (d1[1]-1);
		o2[2] = orig[2];
		break;
	case 'z':
		d1[0] = d2[0] = dim[0];
		d1[1] = d2[1] = dim[1];
		d1[2] = (dim[2] >> 1) + 1;
		d2[2] = ((dim[2]-1) >> 1) + 1;
		o2[0] = orig[0];
		o2[1] = orig[1];
		o2[2] = orig[2] + (d1[2]-1);
		break;
	}
	p_mesh1 = new VolumeReg3Critical(p_fun, d1, orig);
	p_mesh2 = new VolumeReg3Critical(p_fun, d2, o2);
}


vector<CriticalPoint>* VolumeReg3Critical::calcLUStars() {
	//if (p_vcp == NULL) {
	//	p_vcp = sortCriticalPoints();
	//}
	if (p_vcp == NULL) {
		sortCriticalPoints(p_vcp);
	}
	int i, j, k, nv = dim[0]*dim[1]*dim[2];
	int *p_map = new int[nv];
	for (i = 0; i < nv; i++) {
		p_map[(*p_vcp)[i].id] = i;
	}

	// for each edge
	for (i = 0; i < nv; i++) {
		int id = (*p_vcp)[i].id;
		int idx[3], nidx[3];
		id2Index(id, idx);
		for (j = 0; j < 14; j++) {
			for (k = 0; k < 3; k++) {
				nidx[k] = idx[k] + six_simplex_neighbors[j][k];
			}
			if (isValidIndex(nidx)) {
				int nid = index2ID(nidx);
				if (i < p_map[nid]) {
					int nj = p_map[nid];
					(*p_vcp)[nj].LS --;
					(*p_vcp)[i].US --;
				}
			}
		}
		int t1 = (*p_vcp)[i].US;
		int t2 = (*p_vcp)[i].LS;
	}

	// for each triangle
	// x planes
	for (i = 0; i < dim[0]; i++) {
		for (j = 0; j < dim[1]-1; j++) {
			for (k = 0; k < dim[2]-1; k++) {
				// each rectangle has two triangles
				bool flag = (i == 0 || i == dim[0]-1);
				LUTriangle(p_vcp, p_map, index2ID(i, j, k), index2ID(i, j+1, k),
						   index2ID(i, j+1, k+1), flag);
				LUTriangle(p_vcp, p_map, index2ID(i, j, k), index2ID(i, j, k+1),
						   index2ID(i, j+1, k+1), flag);
			}
		}
	}
	// y planes
	for (j = 0; j < dim[1]; j++) {
		for (i = 0; i < dim[0]-1; i++) {
			for (k = 0; k < dim[2]-1; k++) {
				bool flag = (j == 0 || j == dim[1]-1);
				LUTriangle(p_vcp, p_map, index2ID(i, j, k), index2ID(i+1, j, k),
						   index2ID(i+1, j, k+1), flag);
				LUTriangle(p_vcp, p_map, index2ID(i, j, k), index2ID(i, j, k+1),
						   index2ID(i+1, j, k+1), flag);
			}
		}
	}
	// z planes
	for (k = 0; k < dim[2]; k++) {
		for (i = 0; i < dim[0]-1; i++) {
			for (j = 0; j < dim[1]-1; j++) {
				bool flag = (k == 0 || k == dim[2]-1);
				LUTriangle(p_vcp, p_map, index2ID(i, j, k), index2ID(i+1, j, k),
						   index2ID(i+1, j+1, k), flag);
				LUTriangle(p_vcp, p_map, index2ID(i, j, k), index2ID(i, j+1, k),
						   index2ID(i+1, j+1, k), flag);
			}
		}
	}
	// inside each cell
	for (i = 0; i < dim[0]-1; i++) {
		for (j = 0; j < dim[1]-1; j++) {
			for (k = 0; k < dim[2]-1; k++) {
				LUTriangle(p_vcp, p_map, index2ID(i, j, k), index2ID(i, j, k+1),
						   index2ID(i+1, j+1, k+1), false);
				LUTriangle(p_vcp, p_map, index2ID(i, j, k), index2ID(i+1, j+1, k),
						   index2ID(i+1, j+1, k+1), false);
				LUTriangle(p_vcp, p_map, index2ID(i, j, k), index2ID(i, j+1, k),
						   index2ID(i+1, j+1, k+1), false);
				LUTriangle(p_vcp, p_map, index2ID(i, j, k), index2ID(i+1, j, k+1),
						   index2ID(i+1, j+1, k+1), false);
				LUTriangle(p_vcp, p_map, index2ID(i, j, k), index2ID(i+1, j, k),
						   index2ID(i+1, j+1, k+1), false);
				LUTriangle(p_vcp, p_map, index2ID(i, j, k), index2ID(i, j+1, k+1),
						   index2ID(i+1, j+1, k+1), false);
			}
		}
	}
	// for each tetrahedron
	for (i = 0; i < dim[0]-1; i++) {
		for (j = 0; j < dim[1]-1; j++) {
			for (k = 0; k < dim[2]-1; k++) {
				LUTetrahedron(p_vcp, p_map, index2ID(i, j, k), index2ID(i+1, j, k),
							  index2ID(i+1, j, k+1), index2ID(i+1, j+1, k+1));
				LUTetrahedron(p_vcp, p_map, index2ID(i, j, k), index2ID(i+1, j, k),
							  index2ID(i+1, j+1, k), index2ID(i+1, j+1, k+1));
				LUTetrahedron(p_vcp, p_map, index2ID(i, j, k), index2ID(i, j+1, k),
							  index2ID(i+1, j+1, k), index2ID(i+1, j+1, k+1));
				LUTetrahedron(p_vcp, p_map, index2ID(i, j, k), index2ID(i, j+1, k),
							  index2ID(i, j+1, k+1), index2ID(i+1, j+1, k+1));
				LUTetrahedron(p_vcp, p_map, index2ID(i, j, k), index2ID(i, j, k+1),
							  index2ID(i+1, j, k+1), index2ID(i+1, j+1, k+1));
				LUTetrahedron(p_vcp, p_map, index2ID(i, j, k), index2ID(i, j, k+1),
							  index2ID(i, j+1, k+1), index2ID(i+1, j+1, k+1));
			}
		}
	}
	delete[] p_map;

	return p_vcp;
}

bool VolumeReg3Critical::isMinima(const CriticalPoint& cp) {
	int idx[3], nidx[3];
	id2Index(cp.id, idx);
	for (int j = 0; j < 14; j++) {
		for (int k = 0; k < 3; k++)	nidx[k] = idx[k] + six_simplex_neighbors[j][k];
		if (isValidIndex(nidx)) {
			int nid = index2ID(nidx);
			CriticalPoint ncp(nid, (*p_fun)[nid]);
			if (ncp < cp) return false;
		}
	}
	return true;
}

bool VolumeReg3Critical::isMaxima(const CriticalPoint& cp) {
	int idx[3], nidx[3];
	id2Index(cp.id, idx);
	for (int j = 0; j < 14; j++) {
		for (int k = 0; k < 3; k++)	nidx[k] = idx[k] + six_simplex_neighbors[j][k];
		if (isValidIndex(nidx)) {
			int nid = index2ID(nidx);
			CriticalPoint ncp(nid, (*p_fun)[nid]);
			if (cp < ncp) return false;
		}
	}
	return true;
}

int VolumeReg3Critical::findNeighbors(int vid, int* neighbors)
{
	int nNeighbors = 0;
	// At most 14 neighbors for a vertex
	int idx3[3], nidx3[3];	

	id2Index(vid, idx3);
	for (int j = 0; j < 14; j++) {
		for (int k = 0; k < 3; k++)	
			nidx3[k] = idx3[k] + six_simplex_neighbors[j][k];
		if (isValidIndex(nidx3)) {
			int nid = index2ID(nidx3);
			neighbors[nNeighbors++] = nid;
		}
	}
	return nNeighbors;
}


void VolumeReg3Critical::orderTriangleVerts(int* p_map, int v1, int v2, int v3, int ordered[3])
{
	ordered[0] = p_map[v1];
	ordered[1] = p_map[v2];
	ordered[2] = p_map[v3];
	qsort(ordered, 3, sizeof(int), compareInt);
}

void VolumeReg3Critical::orderTetrahedronVerts(int *p_map, int v1, int v2, int v3, int v4, int ordered[4])
{
	ordered[0] = p_map[v1];
	ordered[1] = p_map[v2];
	ordered[2] = p_map[v3];
	ordered[3] = p_map[v4];
	qsort(ordered, 4, sizeof(int), compareInt);
}

void VolumeReg3Critical::LUTriangle(vector<CriticalPoint>* p_v, int *p_map, int v1, int v2, int v3, bool bounflag)
{
	int tri[3];
	orderTriangleVerts(p_map, v1, v2, v3, tri);
	(*p_v)[tri[2]].LS ++;
	(*p_v)[tri[0]].US ++;
	// boundary triangle
	if (bounflag) {
		(*p_v)[tri[2]].dbe --;
		(*p_v)[tri[0]].dbe ++;
	}
}

void VolumeReg3Critical::LUTetrahedron(vector<CriticalPoint>* p_v, int *p_map, int v1, int v2, int v3, int v4)
{
	int tetra[4];
	orderTetrahedronVerts(p_map, v1, v2, v3, v4, tetra);
	(*p_v)[tetra[3]].LS --;
	(*p_v)[tetra[0]].US --;
}

bool VolumeReg3Critical::isValidIndex(int idx[3])
{
	return(idx[0] >= 0 && idx[1] >= 0 && idx[2] >= 0 &&
		   idx[0] < g_dim[0] && idx[1] < g_dim[1] && idx[2] < g_dim[2]);
}

void VolumeReg3Critical::id2Index(int id, int idx3[3])
{
	idx3[0] = id % g_dim[0];
	idx3[1] = (id / g_dim[0]) % g_dim[1];
	idx3[2] = id / (g_dim[0]*g_dim[1]);
}

int VolumeReg3Critical::tetNeighbor(int tid, int side)
{
	int cid = tid / 6;
	int tet = tid % 6;
	int idx[3];

	cell2Index(cid, idx);

	idx[0] += six_tetra_neighbors[tet][side][0];
	idx[1] += six_tetra_neighbors[tet][side][1];
	idx[2] += six_tetra_neighbors[tet][side][2];

	if (!isValidCell(idx)) return -1;
	return(6*index2cell(idx[0], idx[1], idx[2]) + six_tetra_neighbors[tet][side][3]);
}


int VolumeReg3Critical::getMaxValNeighbor(int vid)
{
	int i, idx[3], nidx[3], v2;
	float lcl_max;
	id2Index(vid, idx);

	v2 = -1;

	for (i = 0; i < 14; i++) {

		nidx[0] = idx[0] + six_simplex_neighbors[i][0];

		nidx[1] = idx[1] + six_simplex_neighbors[i][1];

		nidx[2] = idx[2] + six_simplex_neighbors[i][2];

		if (isValidIndex(nidx)) {

			int nid = index2ID(nidx);

			if (v2 < 0 || getValue(nid) > lcl_max ) {

				v2 = nid;

				lcl_max = getValue(nid);

			}

		}

	}

	return v2;

}



int VolumeReg3Critical::getMinValNeighbor(int vid)
{
	int i, idx[3], nidx[3], v2;
	float lcl_min;

	id2Index(vid, idx);
	v2 = -1;
	for (i = 0; i < 14; i++) {
		nidx[0] = idx[0] + six_simplex_neighbors[i][0];

		nidx[1] = idx[1] + six_simplex_neighbors[i][1];

		nidx[2] = idx[2] + six_simplex_neighbors[i][2];

		if (isValidIndex(nidx)) {

			int nid = index2ID(nidx);

			if (v2 < 0 || getValue(nid) < lcl_min ) {

				v2 = nid;

				lcl_min = getValue(nid);

			}

		}

	}

	return v2;

}


int VolumeReg3Critical::getCutEdge(float cut_val, float range[2], int &v1) {
	int v2 = -1;


	if (getValue(v1) == cut_val) {

		assert(getValue(v1) == range[0] || getValue(v1) == range[1]);

		if (getValue(v1) <= range[0]) {

			v2 = getMaxValNeighbor(v1);

		} else if (getValue(v1) >= range[1]) {

			v2 = getMinValNeighbor(v1);

		}

		return v2;

	}


	if (getValue(v1) < cut_val) {
		if (getValue(v1) >= range[0]) {

			v2 = getMaxValNeighbor(v1);

			return v2;

		}

		v2 = v1;
		while (getValue(v2) <= cut_val) {
			//search for v1's neighbors that has the max value > f(v1)
			v1 = v2;

			v2 = getMaxValNeighbor(v1);
		}
	} else if (getValue(v1) > cut_val) {
		if (getValue(v1) <= range[1]) {

			v2 = getMinValNeighbor(v1);

			return v2;
		}

		v2 = v1;

		while (getValue(v2) > cut_val) {

			v1 = v2;

			v2 = getMinValNeighbor(v1);

		}
	}

	return v2;
}

int VolumeReg3Critical::getTetraIndex(int idx1[3], int idx2[3])
{
	if (!areNeighbors(idx1, idx2)) {
		printf("non-neighboring verts, special case....\n");
		int k, dx[3], nidx[3];
		for (k = 0; k < 3; k++) {
			dx[k] = idx2[k] - idx1[k];
		}
		bool flag = false;
		int angle = 0;
		float val = getValue(index2ID(idx1));
		for (k = 0; k < 14; k++) {
			nidx[0] = idx1[0] + six_simplex_neighbors[k][0];
			nidx[1] = idx1[1] + six_simplex_neighbors[k][1];
			nidx[2] = idx1[2] + six_simplex_neighbors[k][2];
			if (getValue(index2ID(nidx)) > val) {
				if ((dx[0]*six_simplex_neighbors[k][0]+dx[1]*six_simplex_neighbors[k][1]
					 +dx[2]*six_simplex_neighbors[k][2]) < angle || (!flag)) {
					idx2[0] = nidx[0];
					idx2[1] = nidx[1];
					idx2[2] = nidx[2];
					angle = dx[0]*six_simplex_neighbors[k][0]+dx[1]*six_simplex_neighbors[k][1]+dx[2]*six_simplex_neighbors[k][2];
					flag = true;
				}
			}
		}
		assert(flag);
	}

	int i, j, cidx[3];
	cidx[0] = MIN(idx1[0], idx2[0]);
	cidx[1] = MIN(idx1[1], idx2[1]);
	cidx[2] = MIN(idx1[2], idx2[2]);

	for (i = 0; i < 3; i++)
		cidx[i] = MIN(cidx[i], g_dim[i]-2);
	int lv1 = (idx1[0]-cidx[0]) + ((idx1[1]-cidx[1]) << 1) + ((idx1[2]-cidx[2]) << 2);
	int lv2 = (idx2[0]-cidx[0]) + ((idx2[1]-cidx[1]) << 1) + ((idx2[2]-cidx[2]) << 2);

	for (i = 0; i < 6; i++) {
		bool flag1, flag2;
		flag1 = flag2 = false;
		for (j = 0; j < 4; j++) {
			if (lv1 == six_tetra_verts[i][j]) flag1 = true;
			if (lv2 == six_tetra_verts[i][j]) flag2 = true;
		}
		if (flag1 && flag2)	break;
	}

	return 6*index2cell(cidx[0], cidx[1], cidx[2])+i;
}

float VolumeReg3Critical::calcVolume(float range[2], int start_tid, float &fint)
{
	float volum = 0;
	int i, idx[3];
	float face_rang[2];

	char *touched = new char[6*(g_dim[0]-1)*(g_dim[1]-1)*(g_dim[2]-1)];
	memset(touched, 0, 6*(g_dim[0]-1)*(g_dim[1]-1)*(g_dim[2]-1));
	deque<int> tetque;
	touched[start_tid] = 1;
	tetque.push_back(start_tid);

	fint = 0;
	while (!tetque.empty()) {
		int tid = tetque.front();
		tetque.pop_front();
		int cid = tid / 6;
		int tet = tid %6;
		cell2Index(cid, idx);
		float vals[4];
		vals[0] = getValue(index2ID(idx[0]+six_tetra_index[tet][0][0], idx[1]+six_tetra_index[tet][0][1], idx[2]+six_tetra_index[tet][0][2]));
		vals[1] = getValue(index2ID(idx[0]+six_tetra_index[tet][1][0], idx[1]+six_tetra_index[tet][1][1], idx[2]+six_tetra_index[tet][1][2]));
		vals[2] = getValue(index2ID(idx[0]+six_tetra_index[tet][2][0], idx[1]+six_tetra_index[tet][2][1], idx[2]+six_tetra_index[tet][2][2]));
		vals[3] = getValue(index2ID(idx[0]+six_tetra_index[tet][3][0], idx[1]+six_tetra_index[tet][3][1], idx[2]+six_tetra_index[tet][3][2]));
		float x0[3], x1[3], x2[3], x3[3];
		for (i = 0; i < 3; i++) {
			x0[i] = m_orig[i]+(idx[i] + six_tetra_index[tet][0][i])*m_span[i];
			x1[i] = m_orig[i]+(idx[i] + six_tetra_index[tet][1][i])*m_span[i];
			x2[i] = m_orig[i]+(idx[i] + six_tetra_index[tet][2][i])*m_span[i];
			x3[i] = m_orig[i]+(idx[i] + six_tetra_index[tet][3][i])*m_span[i];
		}
		// adjust values if they are equal
		float tmax = MAX(MAX(vals[0], vals[1]), MAX(vals[2], vals[3]));
		float tmin = MIN(MIN(vals[0], vals[1]), MIN(vals[2], vals[3]));
		//epsilon = MIN((p_fun->getFuncMax()-p_fun->getFuncMin())/(p_fun->getVertCount()*100), (tmax-tmin)/400);
		//printf("fmax: %f %f, fmin: %f %f, esp: %f\n", p_fun->getFuncMax(), tmax, p_fun->getFuncMin(), tmin, epsilon);
		float delta = subVolume(x0, x1, x2, x3, vals, range);
		//float delta = tetVolume(x0, x1, x2, x3, vals[0], vals[1], vals[2], vals[3], range[1]) -
		//			  tetVolume(x0, x1, x2, x3, vals[0], vals[1], vals[2], vals[3], range[0]);
		volum = volum + delta;
		if (p_pot) {
			float f[4];
			f[0] = p_pot->getValue(index2ID(idx[0]+six_tetra_index[tet][0][0], idx[1]+six_tetra_index[tet][0][1], idx[2]+six_tetra_index[tet][0][2]));
			f[1] = p_pot->getValue(index2ID(idx[0]+six_tetra_index[tet][1][0], idx[1]+six_tetra_index[tet][1][1], idx[2]+six_tetra_index[tet][1][2]));
			f[2] = p_pot->getValue(index2ID(idx[0]+six_tetra_index[tet][2][0], idx[1]+six_tetra_index[tet][2][1], idx[2]+six_tetra_index[tet][2][2]));
			f[3] = p_pot->getValue(index2ID(idx[0]+six_tetra_index[tet][3][0], idx[1]+six_tetra_index[tet][3][1], idx[2]+six_tetra_index[tet][3][2]));    
			float tint = subFuncInteg(x0, x1, x2, x3, vals, f, range);
			fint += tint;
			//printf("dv: %f, df: %f; vol = %f, fint = %f\n", delta, tint, volum, fint);
		}
#ifdef DEBUG_VOLVOL
		//if (delta < 0) {
		printf("tetra %d, range (%f, %f), vals = (%f, %f, %f, %f), vol = %f, delta = %f\n", 
			   tid, range[0], range[1], vals[0], vals[1], vals[2], vals[3], volum, delta);
		/*printf("x0 = (%f, %f, %f)\n", x0[0], x0[1], x0[2]);
		printf("x1 = (%f, %f, %f)\n", x1[0], x1[1], x1[2]);
		printf("x2 = (%f, %f, %f)\n", x2[0], x2[1], x2[2]);
		printf("x3 = (%f, %f, %f)\n", x3[0], x3[1], x3[2]);
		}*/
#endif
		for (i = 0; i < 4; i++) {
			face_rang[0] = MIN(vals[(i+1)%4], MIN(vals[(i+2)%4], vals[(i+3)%4]));
			face_rang[1] = MAX(vals[(i+1)%4], MAX(vals[(i+2)%4], vals[(i+3)%4]));
			if (!disjointRange(range, face_rang)) {
				int ngb_tid = tetNeighbor(tid, i);
				if (ngb_tid >= 0 && !touched[ngb_tid]) {
					touched[ngb_tid] = 1;
					tetque.push_back(ngb_tid);
				}
			}
		}
	}

	delete[] touched;

	double vol_norm = 1.0 / getVolume();
	fint = float (vol_norm * fint);
	return float (volum * vol_norm);
	//return volum;
}

#define N_MOM 20

VolMoments VolumeReg3Critical::calcMoments(float range[2], int start_tid) 
{
	VolMoments mom;

	double volum = 0;
	int i, idx[3];
	float face_rang[2];

	char *touched = new char[6*(g_dim[0]-1)*(g_dim[1]-1)*(g_dim[2]-1)];
	memset(touched, 0, 6*(g_dim[0]-1)*(g_dim[1]-1)*(g_dim[2]-1));
	deque<int> tetque;
	touched[start_tid] = 1;
	tetque.push_back(start_tid);

	while (!tetque.empty()) {
		int tid = tetque.front();
		tetque.pop_front();
		int cid = tid / 6;
		int tet = tid %6;
		cell2Index(cid, idx);
		float vals[4], pot[4] = {0, 0, 0, 0};
		vals[0] = getValue(index2ID(idx[0]+six_tetra_index[tet][0][0], idx[1]+six_tetra_index[tet][0][1], idx[2]+six_tetra_index[tet][0][2]));
		vals[1] = getValue(index2ID(idx[0]+six_tetra_index[tet][1][0], idx[1]+six_tetra_index[tet][1][1], idx[2]+six_tetra_index[tet][1][2]));
		vals[2] = getValue(index2ID(idx[0]+six_tetra_index[tet][2][0], idx[1]+six_tetra_index[tet][2][1], idx[2]+six_tetra_index[tet][2][2]));
		vals[3] = getValue(index2ID(idx[0]+six_tetra_index[tet][3][0], idx[1]+six_tetra_index[tet][3][1], idx[2]+six_tetra_index[tet][3][2]));
		if(p_pot) {
			pot[0] = p_pot->getValue(index2ID(idx[0]+six_tetra_index[tet][0][0], idx[1]+six_tetra_index[tet][0][1], idx[2]+six_tetra_index[tet][0][2]));
			pot[1] = p_pot->getValue(index2ID(idx[0]+six_tetra_index[tet][1][0], idx[1]+six_tetra_index[tet][1][1], idx[2]+six_tetra_index[tet][1][2]));
			pot[2] = p_pot->getValue(index2ID(idx[0]+six_tetra_index[tet][2][0], idx[1]+six_tetra_index[tet][2][1], idx[2]+six_tetra_index[tet][2][2]));
			pot[3] = p_pot->getValue(index2ID(idx[0]+six_tetra_index[tet][3][0], idx[1]+six_tetra_index[tet][3][1], idx[2]+six_tetra_index[tet][3][2]));
		}
		float x0[3], x1[3], x2[3], x3[3];
		for (i = 0; i < 3; i++) {
			x0[i] = m_orig[i]+(idx[i] + six_tetra_index[tet][0][i])*m_span[i];
			x1[i] = m_orig[i]+(idx[i] + six_tetra_index[tet][1][i])*m_span[i];
			x2[i] = m_orig[i]+(idx[i] + six_tetra_index[tet][2][i])*m_span[i];
			x3[i] = m_orig[i]+(idx[i] + six_tetra_index[tet][3][i])*m_span[i];
		}
		
		float f1[N_MOM], f2[N_MOM], f3[N_MOM], f4[N_MOM], out[N_MOM];
		f1[0] = f2[0] = f3[0] = f4[0] = 1;
        f1[1] = x0[0]; f2[1] = x1[0]; f3[1] = x2[0]; f4[1] = x3[0];
        f1[2] = x0[1]; f2[2] = x1[1]; f3[2] = x2[1]; f4[2] = x3[1];
        f1[3] = x0[2]; f2[3] = x1[2]; f3[3] = x2[2]; f4[3] = x3[2];
        f1[4] = x0[0]*x0[0]; f2[4] = x1[0]*x1[0]; f3[4] = x2[0]*x2[0]; f4[4] = x3[0]*x3[0];
        f1[5] = x0[1]*x0[1]; f2[5] = x1[1]*x1[1]; f3[5] = x2[1]*x2[1]; f4[5] = x3[1]*x3[1];
        f1[6] = x0[2]*x0[2]; f2[6] = x1[2]*x1[2]; f3[6] = x2[2]*x2[2]; f4[6] = x3[2]*x3[2];
        f1[7] = x0[0]*x0[1]; f2[7] = x1[0]*x1[1]; f3[7] = x2[0]*x2[1]; f4[7] = x3[0]*x3[1];
        f1[8] = x0[0]*x0[2]; f2[8] = x1[0]*x1[2]; f3[8] = x2[0]*x2[2]; f4[8] = x3[0]*x3[2];
        f1[9] = x0[1]*x0[2]; f2[9] = x1[1]*x1[2]; f3[9] = x2[1]*x2[2]; f4[9] = x3[1]*x3[2];                   
		f1[10] = pot[0]*x0[0]; f2[10] = pot[1]*x1[0]; f3[10] = pot[2]*x2[0]; f4[10] = pot[3]*x3[0];
        f1[11] = pot[0]*x0[1]; f2[11] = pot[1]*x1[1]; f3[11] = pot[2]*x2[1]; f4[11] = pot[3]*x3[1];
        f1[12] = pot[0]*x0[2]; f2[12] = pot[1]*x1[2]; f3[12] = pot[2]*x2[2]; f4[12] = pot[3]*x3[2];
        f1[13] = pot[0]*x0[0]*x0[0]; f2[13] = pot[1]*x1[0]*x1[0]; f3[13] = pot[2]*x2[0]*x2[0]; f4[13] = pot[3]*x3[0]*x3[0];
        f1[14] = pot[0]*x0[1]*x0[1]; f2[14] = pot[1]*x1[1]*x1[1]; f3[14] = pot[2]*x2[1]*x2[1]; f4[14] = pot[3]*x3[1]*x3[1];
        f1[15] = pot[0]*x0[2]*x0[2]; f2[15] = pot[1]*x1[2]*x1[2]; f3[15] = pot[2]*x2[2]*x2[2]; f4[15] = pot[3]*x3[2]*x3[2];
        f1[16] = pot[0]*x0[0]*x0[1]; f2[16] = pot[1]*x1[0]*x1[1]; f3[16] = pot[2]*x2[0]*x2[1]; f4[16] = pot[3]*x3[0]*x3[1];
        f1[17] = pot[0]*x0[0]*x0[2]; f2[17] = pot[1]*x1[0]*x1[2]; f3[17] = pot[2]*x2[0]*x2[2]; f4[17] = pot[3]*x3[0]*x3[2];
        f1[18] = pot[0]*x0[1]*x0[2]; f2[18] = pot[1]*x1[1]*x1[2]; f3[18] = pot[2]*x2[1]*x2[2]; f4[18] = pot[3]*x3[1]*x3[2];                       
		f1[19] = pot[0]; f2[19] = pot[1]; f3[19] = pot[2]; f4[19] = pot[3];
		
		subFuncIntegMul(x0, x1, x2, x3, vals, f1, f2, f3, f4, out, N_MOM, range);
		mom.vol += out[0];
		mom.Px += out[1];
		mom.Py += out[2];
		mom.Pz += out[3];
		mom.Pxx += out[4];
		mom.Pyy += out[5];
		mom.Pzz += out[6];
		mom.Pxy += out[7];
		mom.Pxz += out[8];
		mom.Pyz += out[9];
		mom.Fx += out[10];
		mom.Fy += out[11];
		mom.Fz += out[12];
		mom.Fxx += out[13];
		mom.Fyy += out[14];
		mom.Fzz += out[15];
		mom.Fxy += out[16];
		mom.Fxz += out[17];
		mom.Fyz += out[18];
		mom.Fint += out[19];
		
		if(mom.Fmin > MIN(MIN(pot[0], pot[1]), MIN(pot[2], pot[3]))) mom.Fmin = MIN(MIN(pot[0], pot[1]), MIN(pot[2], pot[3]));
		if(mom.Fmax < MAX(MAX(pot[0], pot[1]), MAX(pot[2], pot[3]))) mom.Fmax = MAX(MAX(pot[0], pot[1]), MAX(pot[2], pot[3]));
		
		for (i = 0; i < 4; i++) {
			face_rang[0] = MIN(vals[(i+1)%4], MIN(vals[(i+2)%4], vals[(i+3)%4]));
			face_rang[1] = MAX(vals[(i+1)%4], MAX(vals[(i+2)%4], vals[(i+3)%4]));
			if (!disjointRange(range, face_rang)) {
				int ngb_tid = tetNeighbor(tid, i);
				if (ngb_tid >= 0 && !touched[ngb_tid]) {
					touched[ngb_tid] = 1;
					tetque.push_back(ngb_tid);
				}
			}
		}
	}

	delete[] touched;
	//mom.print();
	return mom;
}

void VolumeReg3Critical::sortVerts(float* &x0, float* &x1, float* &x2, float* &x3, float v[4], float f[4])
{
	float *_t, t;
	if (v[3] < v[2]) {
		_t = x3; x3 = x2; x2 = _t;
		t = v[3]; v[3] = v[2]; v[2] = t;
		t = f[3]; f[3] = f[2]; f[2] = t;
	}
	if (v[2] < v[1]) {
		_t = x2; x2 = x1; x1 = _t;
		t = v[2]; v[2] = v[1]; v[1] = t;
		t = f[2]; f[2] = f[1]; f[1] = t;
	}
	if (v[1] < v[0]) {
		_t = x1; x1 = x0; x0 = _t;
		t = v[1]; v[1] = v[0]; v[0] = t;
		t = f[1]; f[1] = f[0]; f[0] = t;
	}
	if (v[3] < v[2]) {
		_t = x3; x3 = x2; x2 = _t;
		t = v[3]; v[3] = v[2]; v[2] = t;
		t = f[3]; f[3] = f[2]; f[2] = t;
	}
	if (v[2] < v[1]) {
		_t = x2; x2 = x1; x1 = _t;
		t = v[2]; v[2] = v[1]; v[1] = t;
		t = f[2]; f[2] = f[1]; f[1] = t;
	}
	if (v[3] < v[2]) {
		_t = x3; x3 = x2; x2 = _t;
		t = v[3]; v[3] = v[2]; v[2] = t;
		t = f[3]; f[3] = f[2]; f[2] = t;
	}
	
	//float epsilon = MIN((p_fun->getFuncMax()-p_fun->getFuncMin())/(p_fun->getVertCount()*100), (v[3]-v[0])/400);
	//epsilon = 0.1*MAX(epsilon, (v[3]-v[0])/4000);
	float epsilon = MAX(1e-5f, (v[3]-v[1])/4000);
	if (v[1] <= v[0]+epsilon) v[1] += epsilon;
	if (v[2] <= v[1]+epsilon) v[2] += 2 * epsilon;
	if (v[3] <= v[2]+epsilon) v[3] += 4 * epsilon;
}

/**
 * Calculate the volume of the portion of a tetrahedron with f < fx
 * @note Assume v1 < v2 < v3 < v4
 */ 
float tetVolume(float x1[3], float x2[3], float x3[3], float x4[3], float v1, float v2,
				float v3, float v4, float fx)
{
	float *_t;
	float mid[3], mid2[3];
	float vec1[3], vec2[3], vec3[3];
	float ival, val;
	float s, s2, s3;
	float area1, area2, midarea, volume;
	float cum1;
	float t;
	float cp[3];
	unsigned int first;

	//assert(v1 < v2 && v2 < v3 && v3 < v4);
	//normalize the values to [0, 1];
	/* v2 = (v2-v1)/(v4-v1);
	v3 = (v3-v1)/(v4-v1);
	fx = (fx-v1)/(v4-v1);
	v1 = 0; v4 = 1; */
	
	if (fx <= v1) return 0;

	vec1[0] = x2[0]-x1[0];
	vec1[1] = x2[1]-x1[1];
	vec1[2] = x2[2]-x1[2];
	vec2[0] = x3[0]-x1[0];
	vec2[1] = x3[1]-x1[1];
	vec2[2] = x3[2]-x1[2];
	vec3[0] = x4[0]-x1[0];
	vec3[1] = x4[1]-x1[1];
	vec3[2] = x4[2]-x1[2];
	cp[0] = vec3[0]*(vec1[1]*vec2[2]-vec1[2]*vec2[1]);
	cp[1] = vec3[1]*(vec1[2]*vec2[0]-vec1[0]*vec2[2]);
	cp[2] = vec3[2]*(vec1[0]*vec2[1]-vec1[1]*vec2[0]);
	//volume = sqrt(cp[0]*cp[0] + cp[1]*cp[1] + cp[2]*cp[2])/6.0;
	volume = fabs(cp[0] + cp[1] + cp[2])/6.0;

	if (fx >= v4) return volume;

	// compute the first area
	if (v1 != v3)
		ival = (v3-v2)/(v3-v1);
	else
		ival = 0.0;
	mid[0] = (1.0-ival)*x3[0] + (ival)*x1[0];
	mid[1] = (1.0-ival)*x3[1] + (ival)*x1[1];
	mid[2] = (1.0-ival)*x3[2] + (ival)*x1[2];
	if (v1 != v4)
		ival = (v4-v2)/(v4-v1);
	else
		ival = 0.0;
	mid2[0] = (1.0-ival)*x4[0] + (ival)*x1[0];
	mid2[1] = (1.0-ival)*x4[1] + (ival)*x1[1];
	mid2[2] = (1.0-ival)*x4[2] + (ival)*x1[2];

	vec1[0] = mid[0]-x2[0];
	vec1[1] = mid[1]-x2[1];
	vec1[2] = mid[2]-x2[2];
	vec2[0] = mid2[0]-x2[0];
	vec2[1] = mid2[1]-x2[1];
	vec2[2] = mid2[2]-x2[2];
	cp[0] = vec1[1]*vec2[2]-vec1[2]*vec2[1];
	cp[1] = vec1[2]*vec2[0]-vec1[0]*vec2[2];
	cp[2] = vec1[0]*vec2[1]-vec1[1]*vec2[0];
	area1 = 0.5 * fabs(sqrt(cp[0]*cp[0] + cp[1]*cp[1] + cp[2]*cp[2]));


	// compute the second area
	if (v2 != v4)
		ival = (v4-v3)/(v4-v2);
	else
		ival = 0.0;
	mid[0] = (1.0-ival)*x4[0] + (ival)*x2[0];
	mid[1] = (1.0-ival)*x4[1] + (ival)*x2[1];
	mid[2] = (1.0-ival)*x4[2] + (ival)*x2[2];
	if (v4 != v1)
		ival = (v4-v3)/(v4-v1);
	else
		ival = 0.0;
	mid2[0] = (1.0-ival)*x4[0] + (ival)*x1[0];
	mid2[1] = (1.0-ival)*x4[1] + (ival)*x1[1];
	mid2[2] = (1.0-ival)*x4[2] + (ival)*x1[2];

	vec1[0] = mid[0]-x3[0];
	vec1[1] = mid[1]-x3[1];
	vec1[2] = mid[2]-x3[2];
	vec2[0] = mid2[0]-x3[0];
	vec2[1] = mid2[1]-x3[1];
	vec2[2] = mid2[2]-x3[2];
	cp[0] = vec1[1]*vec2[2]-vec1[2]*vec2[1];
	cp[1] = vec1[2]*vec2[0]-vec1[0]*vec2[2];
	cp[2] = vec1[0]*vec2[1]-vec1[1]*vec2[0];
	area2 = 0.5 * fabs(sqrt(cp[0]*cp[0] + cp[1]*cp[1] + cp[2]*cp[2]));

	if (v2-v1 >= v4-v3)
		midarea = area1 * (1.0+(v3-v2)/(v2-v1));
	else if (v4-v3 > v2-v1)
		midarea = area2 * (1.0+(v3-v2)/(v4-v3));
	else {
		assert(0);
		// This formula is wrong -- xiaoyu
		// have to compute the midarea
		vec1[0] = (x2[0]-x1[0])/2;
		vec1[1] = (x2[1]-x1[1])/2;
		vec1[2] = (x2[2]-x1[2])/2;
		vec2[0] = (x4[0]-x3[0])/2;
		vec2[1] = (x4[1]-x3[1])/2;
		vec2[2] = (x4[2]-x3[2])/2;
		cp[0] = (vec1[1]*vec2[2]-vec1[2]*vec2[1]);
		cp[1] = (vec1[2]*vec2[0]-vec1[0]*vec2[2]);
		cp[2] = (vec1[0]*vec2[1]-vec1[1]*vec2[0]);
		midarea = 2*(sqrt(cp[0]*cp[0] + cp[1]*cp[1] + cp[2]*cp[2])) - (area1+area2)/2;
#ifdef DEBUG_VOLVOL
		printf("1: %f 2: %f mid: %f\n", area1, area2, midarea);
		printf("v1: %f v2: %f v3: %f v4: %f\n", v1, v2, v3, v4);
#endif
	}

#ifdef DEBUG_VOLVOL
	printf("1: %f 2: %f mid: %f\n", area1, area2, midarea);
	printf("v1: %f v2: %f v3: %f v4: %f\n", v1, v2, v3, v4);
#endif
	float vol2 = ((v3-v1)*area1 + (v3-v2)*midarea + (v4-v2)*area2)/3;
	float factor = volume / vol2;

	// compute the
	if (fx < v2) {
		if (v1 == v2)
			return 0;
		else {
			s = (fx-v1)/(v2-v1);
			return factor*s*s*s*area1*(v2-v1)/3.0;
		}
	}

	cum1 = area1*(v2-v1)/3.0;

	if (fx < v3) {
		s = (fx - v2)/(float)(v3-v2);
		s2 = s*s;
		s3 = s2*s;
		val = cum1 + (area1*(s - s2 + s3/3.0)
					  + 2*midarea*(s2/2.0 - s3/3.0)
					  + area2*(s3/3.0)) * (v3-v2);
		return factor*val;
	}

	cum1 += (area1/3.0 + midarea/3.0 + area2/3.0) * (v3-v2);

	if (fx < v4) {
		if (v4==v2)
			val = area2;
		else {
			s = (fx-v3)/(float)(v4-v3);
			s2 = s*s;
			s3 = s2*s;
			val = cum1 + (area2*(s - s2 + s3/3.0)) * (v4-v3);
#ifdef DEBUG_VOLVOL
			printf("val (%f) = %f\n", fx, (1.0-s)*(1.0-s) * area2);
#endif
		}
		return factor * val;
	}

	cum1 += area2/3.0 * (v4-v3);

	return factor * cum1;
} 				
		
/****************************************************************/
/* Carbo Index 													*/
/****************************************************************/

float CarboIndex(const VolumeReg3Critical& reg1, float **R1, float *C1,
				 const VolumeReg3Critical& reg2, float **R2, float *C2)
{
	int i, j, k;
	// A = R2*R1^T
	float **A = matrix(0, 2, 0, 2);
	for(i = 0; i < 3; i++) {
		for(j = 0; j < 3; j++) {
			A[i][j] = 0;
			for(k = 0; k < 3; k++) {
				A[i][j] += R2[i][k]*R1[j][k];
			}
		}
	}
	
	// compute |f1| and |f2|
	float cv1 = reg1.m_span[0]*reg1.m_span[1]*reg1.m_span[2];
	float cv2 = reg2.m_span[0]*reg2.m_span[1]*reg2.m_span[2];
	double x = 0;
	for(k = 0; k < reg1.dim[2]; k++) {
		for(j = 0; j < reg1.dim[1]; j++) {
			for(i = 0; i < reg1.dim[0]; i++) {
				x += reg1.getValue(i, j, k, 1) * reg1.getValue(i, j, k, 1);
				//printf("%d %d %d: f = %f, x = %f\n", i, j, k, reg1.getValue(i, j, k, 1), x);
			}
		}
	}
	float norm_f1 = sqrt(x*cv1);
	x = 0;
	for(k = 0; k < reg2.dim[2]; k++) {
		for(j = 0; j < reg2.dim[1]; j++) {
			for(i = 0; i < reg2.dim[0]; i++) {
				x += reg2.getValue(i, j, k, 1) * reg2.getValue(i, j, k, 1);
				//printf("%d %d %d: f = %f, x = %f\n", i, j, k, reg2.getValue(i, j, k, 1), x);
			}
		}
	}
	float norm_f2 = sqrt(x*cv2);
	printf("f1 norm = %.8f, f2 norm = %.8f\n", norm_f1, norm_f2);
	
	// compute dot product <f1, f2>
	x = 0;
	float coord[3], pnt[3];
	for(k = 0; k < reg1.dim[2]; k++) {
		for(j = 0; j < reg1.dim[1]; j++) {
			for(i = 0; i < reg1.dim[0]; i++) {
				coord[0] = reg1.m_orig[0] + i*reg1.m_span[0] - C1[0];
				coord[1] = reg1.m_orig[1] + j*reg1.m_span[1] - C1[1];
				coord[2] = reg1.m_orig[2] + k*reg1.m_span[2] - C1[2];
				float f1 = reg1.getValue(i, j, k, 1);
				for(int l = 0; l < 3; l++) {
					pnt[l] = A[l][0]*coord[0]+A[l][1]*coord[1]+A[l][2]*coord[2] + C2[l];
				}
				float f2 = reg2.getValue(pnt, 1);
				//printf("p1 = (%f %f %f) f1 = %f -> p2 = (%f %f %f) f2 = %f\n", coord[0]+C1[0], coord[1]+C1[1],
				//		coord[2]+C1[2], f1, pnt[0], pnt[1], pnt[2], f2);
				x += f1*f2;
			}
		}
	}
	x *= cv1;
	float carbo = x / (norm_f1*norm_f2);
	printf("Carbo index: %.8f\n", carbo);
	free_matrix(A, 0, 2, 0, 2);
	return carbo;
}
