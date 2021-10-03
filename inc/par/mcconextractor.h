#ifndef _MC_CONTOUREXTRACTOR_H
#define _MC_CONTOUREXTRACTOR_H

#include "conextractor.h"

#include <map>
using namespace std;

/**
* Each edge in regular mesh is indexed by its starting point and direction.
*/
class EdgeIndex {
public:
    EdgeIndex(int i, int j, int k, int _dir) {
        vid[0] = i;
        vid[1] = j;
        vid[2] = k;
        dir = _dir;
    }
	
    int vid[3];
    int dir;
};

/**
* Less than operator for edgeindex.
*/
struct LTEdgeIndex {
    bool operator ()(const EdgeIndex& e1, const EdgeIndex& e2) const
    {
        return ((e1.vid[2] < e2.vid[2]) ||
			((e1.vid[2] == e2.vid[2]) && (e1.vid[1] < e2.vid[1])) ||
			((e1.vid[2] == e2.vid[2]) && (e1.vid[1] == e2.vid[1]) && (e1.vid[0] < e2.vid[0])) ||
			((e1.vid[2] == e2.vid[2]) && (e1.vid[1] == e2.vid[1]) && (e1.vid[0] == e2.vid[0]) && (e1.dir < e2.dir)));
    }
};

/**
 * Marching Cube based contour extraction.
 */
class MCContourExtractor : public ContourExtractor {
public:
	
    MCContourExtractor();
	
	virtual ~MCContourExtractor() {}

    virtual Surface3D * contourReg3Data(Reg3Data* p_reg3data, float isoval);

protected:
	
    int interpRegEdges(float val[8], float grad[8][3], float isoval, 
		int i, int j, int k, int edge, Reg3Data* p_reg3data);
	
private:    
    Surface3D* surf;
    int nvert;
    map<EdgeIndex, int, LTEdgeIndex>* edgeset;
};

#endif

