#pragma once

#include <vector>
#include <algorithm>

#include "CriticalPoint.h"

class ContourTree;

using namespace std;

/**
 * Volumetric Data structure for contour tree construction
 */
class VolCritical
{
public:

	VolCritical(void){
		p_vcp = 0;
	}

	virtual ~VolCritical(void); 

	virtual ContourTree* joinTree(void);
	virtual ContourTree* splitTree(void);
	virtual ContourTree* computeContourTree();

protected:
	vector<CriticalPoint>* p_vcp;

	virtual int getNVerts() const = 0;
	virtual float getValue(int i) const = 0;

	virtual bool isMinima(const CriticalPoint& cp) = 0;
    virtual bool isMaxima(const CriticalPoint& cp) = 0;
	virtual int findNeighbors(int vid, int* neighbors) = 0;

	// Sort the critical points of volume
	void sortCriticalPoints(vector<CriticalPoint>*& p_v);
	virtual vector<CriticalPoint>* calcLUStars() = 0;

	static const int MAX_NEI_NUM;
};
