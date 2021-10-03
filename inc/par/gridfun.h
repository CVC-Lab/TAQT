#ifndef _GRID_FUNCTION_H
#define _GRID_FUNCTION_H

#include "function.h"
#include "reg3data.h"

/**
 *
 */
class GridFun : public Function {
public:
	GridFun(const char* fname);

	virtual ~GridFun();

	virtual void evalFunction(const coord3d& pnt, float &val);

	virtual void evalGradient(const coord3d& pnt, float grad[3]);

	virtual void evalFuncGrad(const coord3d& pnt, float& val, float grad[3]);

	virtual bool valid(const coord3d& pnt);

	virtual void getBoundingBox(float min[3], float max[3]) {
		for(int i = 0; i < 3; i++) {
			min[i] = m_min[3];
			max[i] = m_max[3];
		}
	}

	virtual void getFuncMinMax(float& min, float& max);

protected:
	Reg3Data* p_reg3;
	float m_min[3], m_max[3];
	int m_dim[3];
	float m_orig[3], m_span[3];
};
#endif

