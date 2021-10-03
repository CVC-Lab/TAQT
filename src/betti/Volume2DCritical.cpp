#include "Volume2DCritical.h"
#include "disjointset.h"

// simplex decomposition of a 2D regular cell
const int Volume2DCritical::simplex_neighbors[6][2] = {
	{0, 1}, {1, 0}, {0, -1}, {-1, 0},
	{1, 1}, {-1, -1}
};

Volume2DCritical::Volume2DCritical(Reg2Data* _fun)
{
	p_fun =_fun;
	p_fun->getDim(m_dim);
	p_fun->getOrig(m_orig);
	p_fun->getSpan(m_span);
}

Volume2DCritical::~Volume2DCritical(void)
{
}

/********************************************************************/
/* Overridden virtual functions                                     */
/********************************************************************/
int Volume2DCritical::getNVerts() const {
	return m_dim[0]*m_dim[1];
}
	
float Volume2DCritical::getValue(int i) const{
	return p_fun->getValue(i);
}

bool Volume2DCritical::isMinima(const CriticalPoint& cp) {
	int nNei, neighbors[6];
	int vid = cp.id;
	
	nNei = findNeighbors(vid, neighbors);
	for(int i = 0; i < nNei; i++) {
		CriticalPoint ncp(neighbors[i], p_fun->getValue(neighbors[i]));
		if(ncp < cp) return false;
	}
	return true;
}
 
bool Volume2DCritical::isMaxima(const CriticalPoint& cp) {
	int nNei, neighbors[6];
	int vid = cp.id;
	
	nNei = findNeighbors(vid, neighbors);
	for(int i = 0; i < nNei; i++) {
		CriticalPoint ncp(neighbors[i], p_fun->getValue(neighbors[i]));
		if(cp < ncp) return false;
	}
	return true;
}

int Volume2DCritical::findNeighbors(int vid, int* neighbors) {
	int nNei = 0;

	// At most 6 neighbors for a 2Dvertex
	int idx2[2], nidx2[2];	

	id2Index(vid, idx2);
	for (int j = 0; j < 6; j++) {
		for (int k = 0; k < 2; k++)	
			nidx2[k] = idx2[k] + simplex_neighbors[j][k];
		if (isValidIndex(nidx2)) {
			int nid = index2ID(nidx2);
			neighbors[nNei++] = nid;
		}
	}
	return nNei;
}
