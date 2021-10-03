#ifndef VOLUME_REG3_H
#define VOLUME_REG3_H

#include <iostream>

#include <reg3data.h>
#include "contourtree.h"
#include "CriticalPoint.h"
#include "moment.h"

#include "VolCritical.h"

using namespace std;

struct MinMax {
	float min, max;
};

float tetVolume(float x1[3], float x2[3], float x3[3], float x4[3], float v1, float v2,
				float v3, float v4, float fx);

float tetraFuncIntegral(float x1[3], float x2[3], float x3[3], float x4[3], 
						float v1, float v2, float v3, float v4,
						float f1, float f2, float f3, float f4, float fx);	

void tetIntegMultiFuncs(float x1[3], float x2[3], float x3[3], float x4[3],
						float v1, float v2, float v3, float v4,
						float *f1, float *f2, float *f3, float *f4, float *out, 
						int n, float fx);
						
MinMax tetFunMinMax(float x1[3], float x2[3], float x3[3], float x4[3], 
					float v1, float v2, float v3, float v4,
					float f1, float f2, float f3, float f4, float range[2]);

/**
* VolumeReg3Critical -- A Class to build contour tree for regular 3D volume
* @author Xiaoyu Zhang
*/

class VolumeReg3Critical : public VolCritical 
{
public:
	VolumeReg3Critical(Reg3Data * _fun, Reg3Data *_fun2 = NULL);
	
	VolumeReg3Critical(Reg3Data* _fun, int _dim[3], int _orig[3], Reg3Data *_fun2 = NULL);
	
	~VolumeReg3Critical(void);
	
	void split(VolumeReg3Critical * & p_mesh1, VolumeReg3Critical * & p_mesh2);
	
	bool isCell(void) {
		return (dim[0] == 2 && dim[1] == 2 && dim[2] == 2);
	}
	
	//virtual ContourTree* joinTree(void);
	//virtual ContourTree* splitTree(void);
	
	void getDimension(int _dim[3]) {
		_dim[0] = dim[0];
		_dim[1] = dim[1];
		_dim[2] = dim[2];
	}
	
	/********************************************************************/
	/* Overridden virtual functions                                     */
	/********************************************************************/
	virtual float getValue(int id) const {
		return (*p_fun)[id];
	}
	
	virtual int getNVerts() const {
		return dim[0]*dim[1]*dim[2];
	}

	virtual vector<CriticalPoint>* calcLUStars();

	float getValue(int i, int j, int k, int nf = 0) const {
		switch(nf) {
			case 0:
			return p_fun->getValue(i, j, k);
			case 1:
			if(p_pot == NULL) return 0;
			return p_pot->getValue(i, j, k);
		}
		return 0;
	}
	
	/**
	* Get the function value at a given coordinates
	* @note return 0 if the coordinates are outside of the bounding box
	*/
	float getValue(float coord[3], int nf = 0) const {
		switch(nf) {
			case 0:
			return p_fun->getValue(coord);
			case 1:
			if(p_pot == NULL) return 0;
			return p_pot->getValue(coord);
		}
		return 0;
	}
	
	void setValue(int id, float x) {
		p_fun->setValue(id, x);
	}
	/**
	*	Calculate the volume of a connected component with given range
	*  @param fint Integral of a second function on the volume
	*/
	float calcVolume(float range[2], int start_tid, float &fint);
	
	/**
	* Calculate volumetric moments for a shell of given range
	*/
	VolMoments calcMoments(float range[2], int start_id);
	
	/*
	*	Return a tetrahedron that contains the given edge
	*	@parameter 
	*/
	int getTetraIndex(int idx1[3], int idx2[3]);
	
	/*
	*	Find an edge that intersect a cut value
	*/
	int getCutEdge(float cut_val, float range[2], int &v1);
	
	/*
	*	Find the neigbhoring vertex of the maximum value
	*/
	int getMaxValNeighbor(int vid);
	
	/*
	*	Find the neigbhoring vertex of the minimum value
	*/
	int getMinValNeighbor(int vid);
	
	/*
	*	The volume of the reg3D data
	*/
	
	float getVolume() const {
		// assumed cube cells
		return (g_dim[0]-1)*(g_dim[1]-1)*(g_dim[2]-1)*m_span[0]*m_span[1]*m_span[2];
	}
	
	/**
	* Calculate Carbo index for a given orientation
	*/
	friend float CarboIndex(const VolumeReg3Critical& reg1, float **R1, float *C1,
	const VolumeReg3Critical& reg2, float **R2, float *C2);
	
	static const int vert_neighboring[15][15]; 
	static const int six_simplex_neighbors[14][3];
	static const int six_tetra_index[6][4][3];
	static const int six_tetra_neighbors[6][4][4];
	static const int six_tetra_verts[6][4];
	
	inline void cell2Index(int id, int idx[3]);
	inline bool isValidCell(int idx[3]);
	inline int index2cell(int i, int j, int k);
	
	bool isValidIndex(int idx[3]);
	inline int index2ID(int idx3[3]);
	inline int index2ID(int i, int j, int k);
	
	void id2Index(int id, int idx3[3]);
	
	virtual bool isMinima(const CriticalPoint& cp);
    virtual bool isMaxima(const CriticalPoint& cp);
	
	virtual int findNeighbors(int vid, int* neighbors);

	inline bool areNeighbors(int v1, int v2);
	inline bool areNeighbors(int idx1[3], int idx2[3]);
	protected:
	int dim[3];
	int orig[3];
	float m_orig[3], m_span[3];
	
	static int g_dim[3];			// global dimension of the whole volume
	Reg3Data*  p_fun;
	Reg3Data * p_pot;				// second data volume (e.g. potential)
	
	bool isMin(int id, float vals[8]);
	
	bool isValid(float x) {
		return (x > 0 && x < 1);
	}

	/**
	* Sort tetra vertices by their values and remove any degeneracy 
	*/
	void sortVerts(float* &x0, float* &x1, float* &x2, float* &x3, float v[4], float f[4]);


	float subFuncInteg(float x1[3], float x2[3], float x3[3], float x4[3], float v[4], float f[4], float range[2])
	{
		float *p1 = x1, *p2 = x2, *p3 = x3, *p4 = x4;
		float vals[4], func[4];
		for(int i = 0; i < 4; i++) { // make the copies so we don't mess with the original 
			vals[i] = v[i];
			func[i] = f[i];
		}
		sortVerts(p1, p2, p3, p4, vals, func);
		return tetraFuncIntegral(p1, p2, p3, p4, vals[0], vals[1], vals[2], vals[3], func[0], func[1], func[2], func[3], range[1]) - 
			   tetraFuncIntegral(p1, p2, p3, p4, vals[0], vals[1], vals[2], vals[3], func[0], func[1], func[2], func[3], range[0]);
	}
	
	void subFuncIntegMul(float x1[3], float x2[3], float x3[3], float x4[3], float v[4],
						  float *f1, float *f2, float *f3, float *f4, float *out, int n, float range[2])
	{
		float *p1 = x1, *p2 = x2, *p3 = x3, *p4 = x4;
		float vals[4], func[4];
		for(int i = 0; i < 4; i++) { // make the copies so we don't mess with the original 
			vals[i] = v[i];
		}
		sortVerts(p1, p2, p3, p4, vals, func);
		float *g1 = f1, *g2 = f2, *g3 = f3, *g4 = f4;
		float vals2[4], func2[4];
		for(int i = 0; i < 4; i++) { // make the copies so we don't mess with the original 
			vals2[i] = v[i];
		}
		sortVerts(g1, g2, g3, g4, vals2, func2);
		float *ou = new float[n];
		float *ol = new float[n];
		tetIntegMultiFuncs(p1, p2, p3, p4, vals[0], vals[1], vals[2], vals[3], g1, g2, g3, g4, ou, n, range[1]);
		tetIntegMultiFuncs(p1, p2, p3, p4, vals[0], vals[1], vals[2], vals[3], g1, g2, g3, g4, ol, n, range[0]);
		for(int i = 0; i < n; i++) {
			out[i] = ou[i] -ol[i];
			//printf("%d %f\n", i, out[i]);
		}
		//assert(out[0] >= 0);
		delete[] ou;
		delete[] ol;
		return;
	}	
		
	float subVolume(float x1[3], float x2[3], float x3[3], float x4[3], float v[4], float range[2])
	{
		float *p1 = x1, *p2 = x2, *p3 = x3, *p4 = x4;
		float vals[4], func[4];
		for(int i = 0; i < 4; i++) { // make the copies so we don't mess with the original 
			vals[i] = v[i];
		}
		sortVerts(p1, p2, p3, p4, vals, func);
		return tetVolume(p1, p2, p3, p4, vals[0], vals[1], vals[2], vals[3], range[1]) - 
			   tetVolume(p1, p2, p3, p4, vals[0], vals[1], vals[2], vals[3], range[0]);
	}
	

	/**
	* The globally unique ID that represents a critical point.
	* @note the critical points can be either a vertex or a saddle point.
	*  It only applies to a cell.
	* @param n the local index of a critical point
	*     
	*                              ----------
	*							  /         /
	*                            /         /	
	*                            ---------
	* @note we use the last three digits to represent the type of a critical point.
	*  vertex : 0, face saddle point: 1, 2, 3 depending on the direction of the face,
	*  body saddle point: 4
	*/
	inline int globalID(int n);
	
	/**
	* Find indices of v1, v2, and v3 in the p_map and output them
	* in increasing order.
	*/
	void orderTriangleVerts(int* p_map, int v1, int v2, int v3, int ordered[3]);
	
	/**
	* Find indices of v1, v2, v3, and v4 in the p_map and output them
	* in increasing order.
	*/
	void orderTetrahedronVerts(int* p_map, int v1, int v2, int v3, int v4, int ordered[4]);
	
	void LUTriangle(vector<CriticalPoint>* p_v, int *p_map, int v1, int v2, int v3, bool bounflag);
	
	void LUTetrahedron(vector<CriticalPoint>* p_v, int *p_map, int v1, int v2, int v3, int v4);
	
	/**
	* get the neighbors of a tetrahedron
	* @param side: 
	* @return the neighboring tetrahedron id 
	*/
	int tetNeighbor(int tid, int side);
	
	/*
	* Check if two ranges intersect each other.	
	*/
	inline bool disjointRange(float rang1[2], float rang2[2]);
};

int VolumeReg3Critical::globalID(int n)
{
	int base = orig[0] + orig[1]*g_dim[0] + orig[2]*g_dim[0]*g_dim[1];
	switch(n) {
		case 0:
		return (base << 3);
		case 1:
		return (base+1) << 3;
		case 2:
		return (base+g_dim[0]+1) << 3;
		case 3:
		return (base+g_dim[0]) << 3;
		case 4:
		return (base+g_dim[0]*g_dim[1]) << 3;
		case 5:
		return (base+g_dim[0]*g_dim[1]+1) << 3;
		case 6:
		return (base+g_dim[0]*g_dim[1]+g_dim[0]+1) << 3;
		case 7:
		return (base+g_dim[0]*g_dim[1]+g_dim[0]) << 3;
		case 8:			// saddle point starts here
		return (base << 3) + 1;				// 1 means orthogonal to x axis
		case 9:
		return ((base+1) << 3) + 1;
		case 10:
		return (base << 3) + 2;				// 2 means orthogonal to y axis
		case 11:
		return ((base+g_dim[0]) << 3) + 2;
		case 12:
		return (base << 3) + 3;				// 3 means orthogonal to z axis
		case 13:
		return((base+g_dim[0]*g_dim[1])) + 3;
		case 14:		// body saddle point
		return (base << 3) + 4;
		default:		// should never come here
		assert(0);
	}
	return 0;
}

int VolumeReg3Critical::index2ID(int idx3[3])
{
	return (idx3[0] + idx3[1]*g_dim[0] + idx3[2]*g_dim[0]*g_dim[1]);
}

int VolumeReg3Critical::index2ID(int i, int j, int k)
{
	return (i + j*g_dim[0] + k*g_dim[0]*g_dim[1]);
}

int VolumeReg3Critical::index2cell(int i, int j, int k)
{
	return (i + j*(g_dim[0]-1) + k*(g_dim[0]-1)*(g_dim[1]-1));
}

void VolumeReg3Critical::cell2Index(int id, int idx3[3])
{
	idx3[0] = id % (g_dim[0]-1);
	idx3[1] = (id / (g_dim[0]-1)) % (g_dim[1]-1);
	idx3[2] = id / ((g_dim[0]-1)*(g_dim[1]-1));
}

bool VolumeReg3Critical::isValidCell(int idx[3])
{
	return (idx[0] >= 0 && idx[1] >= 0 && idx[2] >= 0 &&
	idx[0] < (g_dim[0]-1) && idx[1] < (g_dim[1]-1) && idx[2] < (g_dim[2]-1));
}

bool VolumeReg3Critical::disjointRange(float rang1[2], float rang2[2])
{
	if(rang2[0] > rang1[1]) return true;
	if(rang2[1] < rang1[0]) return true;
	if(rang1[0] > rang2[1]) return true;
	if(rang1[1] < rang2[0]) return true;
	return false;
}

bool VolumeReg3Critical::areNeighbors(int v1, int v2)
{
	int idx1[3], idx2[3];
	id2Index(v1, idx1);
	id2Index(v2, idx2);
	
	for(int i = 0; i < 14; i++) {
		if((idx1[0]+six_simplex_neighbors[i][0] == idx2[0]) &&
			(idx1[1]+six_simplex_neighbors[i][1] == idx2[1]) &&
		(idx1[2]+six_simplex_neighbors[i][2] == idx2[2]))
		return true;
	}
	return false;
}

bool VolumeReg3Critical::areNeighbors(int idx1[3], int idx2[3])
{
	for(int i = 0; i < 14; i++) {
		if((idx1[0]+six_simplex_neighbors[i][0] == idx2[0]) &&
			(idx1[1]+six_simplex_neighbors[i][1] == idx2[1]) &&
		(idx1[2]+six_simplex_neighbors[i][2] == idx2[2]))
		return true;
	}
	return false;
}
#endif











