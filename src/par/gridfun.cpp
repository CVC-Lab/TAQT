#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <bio.h>
#include "gridfun.h"

GridFun::GridFun(const char* fname)
{
	p_reg3 = new Reg3Data(fname);

	p_reg3->getBoundingBox(m_min, m_max);
	p_reg3->getDim(m_dim);
	p_reg3->getOrig(m_orig);
	p_reg3->getSpan(m_span);
	//p_reg3->calcGradient();
}

GridFun::~GridFun()
{
	delete p_reg3;
}

void GridFun::evalFunction(const coord3d& pnt, float &val)
{
	// trilinear interpolation
	int i, j, k;
	float u[2], v[2], w[2];
	float data[8];

	if(!valid(pnt)) {
		val = -1;
		return;
	}
	
	i = (int)((pnt[0] - m_orig[0]) / m_span[0]);
	j = (int)((pnt[1] - m_orig[1]) / m_span[1]);
	k = (int)((pnt[2] - m_orig[2]) / m_span[2]);
	p_reg3->getCellValues(i, j, k, data);

	u[1] = (pnt[0] - m_orig[0])/m_span[0] - i; u[0] = 1 - u[1];
	v[1] = (pnt[1] - m_orig[1])/m_span[1] - j; v[0] = 1 - v[1];
	w[1] = (pnt[2] - m_orig[2])/m_span[2] - k; w[0] = 1 - w[1];

	val = 0;
	for(int kk = 0; kk < 2; kk++)
		for(int jj = 0; jj < 2; jj++)
			for(int ii = 0; ii < 2; ii++) {
				val += 	data[vertinfo[kk][jj][ii]]*u[ii]*v[jj]*w[kk];
			}
}

void GridFun::evalGradient(const coord3d& pnt, float grad[3])
{
	// trilinear interpolation
	int i, j, k;
	float u[2], v[2], w[2];
	float grads[8][3];

	if(!valid(pnt)) {
		grad[0] = grad[1] = grad[2] = 0;
		return;
	}
	i = (pnt[0] - m_orig[0]) / m_span[0];
	j = (pnt[1] - m_orig[1]) / m_span[1];
	k = (pnt[2] - m_orig[2]) / m_span[2];
	p_reg3->getCellGrads(i, j, k, grads);

	u[1] = (pnt[0] - m_orig[0])/m_span[0] - i; u[0] = 1 - u[1];
	v[1] = (pnt[1] - m_orig[1])/m_span[1] - j; v[0] = 1 - v[1];
	w[1] = (pnt[2] - m_orig[2])/m_span[2] - k; w[0] = 1 - w[1];

	grad[0] = grad[1] = grad[2] = 0;
	for(int kk = 0; kk < 2; kk++)
		for(int jj = 0; jj < 2; jj++)
			for(int ii = 0; ii < 2; ii++) {
				grad[0] += grads[vertinfo[kk][jj][ii]][0]*u[ii]*v[jj]*w[kk]; 
				grad[1] += grads[vertinfo[kk][jj][ii]][1]*u[ii]*v[jj]*w[kk];
				grad[2] += grads[vertinfo[kk][jj][ii]][2]*u[ii]*v[jj]*w[kk];
			}
}

void GridFun::evalFuncGrad(const coord3d& pnt, float& val, float grad[3])
{
	int i, j, k;
	float u[2], v[2], w[2];
	float data[8], grads[8][3];

	if(!valid(pnt)) {
		val = -1;
		grad[0] = grad[1] = grad[2] = 0;
		return;
	}
	
	i = (pnt[0] - m_orig[0]) / m_span[0];
	j = (pnt[1] - m_orig[1]) / m_span[1];
	k = (pnt[2] - m_orig[2]) / m_span[2];
	p_reg3->getCellValues(i, j, k, data);
	p_reg3->getCellGrads(i, j, k, grads);

	u[1] = (pnt[0] - m_orig[0])/m_span[0] - i; u[0] = 1 - u[1];
	v[1] = (pnt[1] - m_orig[1])/m_span[1] - j; v[0] = 1 - v[1];
	w[1] = (pnt[2] - m_orig[2])/m_span[2] - k; w[0] = 1 - w[1];
	
	val = 0; grad[0] = grad[1] = grad[2] = 0;
	for(int kk = 0; kk < 2; kk++)
		for(int jj = 0; jj < 2; jj++)
			for(int ii = 0; ii < 2; ii++) {
				val += 	data[vertinfo[kk][jj][ii]]*u[ii]*v[jj]*w[kk];
				grad[0] += grads[vertinfo[kk][jj][ii]][0]*u[ii]*v[jj]*w[kk]; 
				grad[1] += grads[vertinfo[kk][jj][ii]][1]*u[ii]*v[jj]*w[kk];
				grad[2] += grads[vertinfo[kk][jj][ii]][2]*u[ii]*v[jj]*w[kk];
			}	
}

bool GridFun::valid(const coord3d& pnt)
{
	return (pnt[0] >= m_min[0] && pnt[0] < m_max[0] && 
			pnt[1] >= m_min[1] && pnt[1] < m_max[1] && 
			pnt[2] >= m_min[2] && pnt[2] < m_max[2]);
}

void GridFun::getFuncMinMax(float &min, float& max)
{
	p_reg3->getFuncMinMax(min, max);
}

