#pragma once
#include "VolCritical.h"
#include <reg2data.h>

class Volume2DCritical :
	public VolCritical
{
public:
	Volume2DCritical(Reg2Data* _fun);
	virtual ~Volume2DCritical(void);

	//virtual ContourTree* joinTree();
	//virtual ContourTree* splitTree();

protected:
	/********************************************************************/
	/* Overridden virtual functions                                     */
	/********************************************************************/
	virtual int getNVerts() const;
	virtual float getValue(int i) const;

	virtual bool isMinima(const CriticalPoint& cp);
    virtual bool isMaxima(const CriticalPoint& cp);
	virtual int findNeighbors(int vid, int* neighbors);

	// Sort the critical points of volume
	virtual vector<CriticalPoint>* calcLUStars() {
		if (p_vcp == NULL) {
			sortCriticalPoints(p_vcp);
		}
		return p_vcp;
	}

	static const int simplex_neighbors[6][2];

	Reg2Data* p_fun;
	unsigned int m_dim[2];
	float m_orig[2], m_span[2];

private:
	/**
	 * Convert vertex id to (ix, iy) index
	 */
	void id2Index(int vid, int idx[2]) {
		idx[0] = vid % m_dim[0];
		idx[1] = vid / m_dim[0];
	}
	/**
	 * Convert from (id[0], id[1]) vertex index to a ID
	 */
	int index2ID(int idx[2]) {
		return (idx[0] + idx[1]*m_dim[0]);
	}
	/**
	 * Determine whether a vertex index (ix, iy) is valid
	 */
	bool isValidIndex(int idx[2]) {
		return ( idx[0] >= 0 && idx[0] < m_dim[0] &&
				 idx[1] >= 0 && idx[1] < m_dim[1]);
	}


};
