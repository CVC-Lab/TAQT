#ifndef _SCALAR_FUNC_H
#define _SCALAR_FUNC_H

#include "geom.h"

/**
 * 
 */
class Function {
public:
	Function() {}

	virtual ~Function() {}

	virtual void evalFunction(const coord3d& pnt, float& val) = 0;

	virtual void evalGradient(const coord3d& pnt, float grad[3]) = 0;

	// evaluate function value and gradient in one function call
	virtual void evalFuncGrad(const coord3d& pnt, float& val, float grad[3]) = 0;

	virtual bool valid(const coord3d& pnt) { return true; }

	virtual void getBoundingBox(float min[3], float max[3]) {
		min[0] = min[1] = min[2] = -1;
		max[0] = max[1] = max[2] = 1;
	}

	virtual void getFuncMinMax(float &min, float &max) {
		min = 0; max = 1;
	}
};

#endif

