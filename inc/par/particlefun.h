#ifndef _PARTICLE_FUNC_H
#define _PARTICLE_FUNC_H

#include "function.h"
#include "dynarray.h"
#include "protein.h"

/**
 *
 */
class ParticleFun : public Function {
public:
	/**
	 * Read atoms from a PDB file.	
	 */
	ParticleFun(const char* pdb_file);

	virtual ~ParticleFun();

	virtual void evalFunction(const coord3d& pnt, float& val);

	virtual void evalGradient(const coord3d& pnt, float grad[3]);

	virtual void evalFuncGrad(const coord3d& pnt, float& val, float grad[3]);

	virtual bool valid(const coord3d& pnt);

	virtual void getBoundingBox(float min[3], float max[3]) {
		for(int i = 0; i < 3; i++) {
			min[i] = bbox_min[i];
			max[i] = bbox_max[i];
		}
	}
protected:
	Protein *p_mol;
	float bbox_min[3], bbox_max[3]; 
};

#endif

