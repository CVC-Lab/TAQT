#include "spherefun.h"

void SphereFun::evalFunction(const coord3d& pnt, float &val)
{
	//val = pnt[0]*pnt[0] + pnt[1]*pnt[1] + pnt[2]*pnt[2];
	val = pnt[0]*pnt[0] + pnt[1]*pnt[1] + pnt[2]*pnt[2];
}

void SphereFun::evalGradient(const coord3d& pnt, float grad[3])
{
	grad[0] = 2 * pnt[0];
	grad[1] = 2 * pnt[1];
	grad[2] = 2 * pnt[2];
}

void SphereFun::evalFuncGrad(const coord3d& pnt, float& val, float grad[3])
{
	val = pnt[0]*pnt[0] + pnt[1]*pnt[1] + pnt[2]*pnt[2];;
	grad[0] = 2 * pnt[0];
	grad[1] = 2 * pnt[1];
	grad[2] = 2 * pnt[2];
}

bool SphereFun::valid(const coord3d& pnt)
{
	return (pnt[0] >= -1 && pnt[0] <= 1 && pnt[1] >= -1 && pnt[1] <= 1 && pnt[2] >= -1 && pnt[2] <= 1);
}

void SphereFun::getFuncMinMax(float &min, float &max)
{
	min = 0;
	max = 3;
}

