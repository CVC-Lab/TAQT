#ifndef _SPHERE_FUNC_H
#define _SPHERE_FUNC_H

#include "function.h"

/**
 *
 */
class SphereFun : public Function {
public:
	SphereFun() {}

	virtual ~SphereFun() {}

	virtual void evalFunction(const coord3d& pnt, float& val);

	virtual void evalGradient(const coord3d& pnt, float grad[3]);

	virtual void evalFuncGrad(const coord3d& pnt, float& val, float grad[3]);

	virtual bool valid(const coord3d& pnt);

	virtual void getFuncMinMax(float &min, float &max);
};

#endif

