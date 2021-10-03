#include <string.h>
#include <math.h>

#include "particlefun.h"
#include "elements.h"

ParticleFun::ParticleFun(const char* pdb_file)
{
	p_mol = new Protein(pdb_file);
	p_mol->getBoundingBox(bbox_min, bbox_max);
}

ParticleFun::~ParticleFun()
{
	delete p_mol;
}

void ParticleFun::evalFunction(const coord3d& pnt, float& val)
{
	val = p_mol->evalDensity(pnt);
}

void ParticleFun::evalGradient(const coord3d& pnt, float grad[3])
{
	p_mol->evalGradient(pnt, grad);
}

void ParticleFun::evalFuncGrad(const coord3d& pnt, float& val, float grad[3])
{
	p_mol->evalDenAndGrad(pnt, val, grad);
}

bool ParticleFun::valid(const coord3d& pnt)
{
	return (pnt[0] >= bbox_min[0] && pnt[0] < bbox_max[0] && 
			pnt[1] >= bbox_min[1] && pnt[1] < bbox_max[1] && 
			pnt[2] >= bbox_min[2] && pnt[2] < bbox_max[2]);
}


