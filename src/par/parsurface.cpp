#include <time.h>
#include <stdlib.h>

#include "gridfun.h"
#include "particlefun.h"
#include "spherefun.h"
#include "parsurface.h"
#include "vecmath.h"
#include "nnbucket.h"
#include "nnkdtree.h"
#include "rand.h"
#include "barycentric.h"

PointSurface::PointSurface(VolType vt)
{
	switch(vt) {
	case SPHERE_FUN:
		p_fun = new SphereFun();
		break;
	default:
		p_fun = NULL;
	}
	init();
}
PointSurface::PointSurface(const char* fname, VolType vt)
{
	VolType type;

	size_t len = strlen(fname);
	if (len > 4 && strcmp(fname+len-4, ".r3d") == 0) {
		type = REG3D_FUN;
	} else if (len > 6 && strcmp(fname+len-6, ".rawiv") == 0) {
		type = REG3D_FUN;
	} else if (len > 4 && strcmp(fname+len-4, ".pdb") == 0) {
		type = PDB_FUN;
	} else {
		type = vt;
	}

	switch (type) {
	case SPHERE_FUN:
		p_fun = new SphereFun();
		break;
	case REG3D_FUN:
		p_fun = new GridFun(fname);
		break;
	case PDB_FUN:
		p_fun = new ParticleFun(fname);
		break;
	default:
		p_fun = NULL;
	}
	init();
}

void PointSurface::init()
{
	p_pnts = NULL;
	m_threshold = 1e-4f;
	m_nsamp = 10000;
}

PointSurface::~PointSurface()
{
	delete p_fun;
	if (p_pnts)	delete p_pnts;
}

void PointSurface::makeSurface(float isoval, int nsample)
{
	if(p_pnts && isoval == m_isoval && nsample == m_nsamp) return;

	if (p_pnts != NULL)	delete p_pnts;
	m_isoval = isoval;
	if(nsample > 0) m_nsamp = nsample;
	p_pnts = new dynamic_array<Point3D>(m_nsamp);

	int i, j, count = 0;
	float min[3], max[3];
	bool *alive, *done;

	p_fun->getBoundingBox(min, max);
	Point3D *points = p_pnts->data();
	alive = new bool[m_nsamp];
	done = new bool[m_nsamp];

	int t1 = time(NULL);
	long seed = -t1;
	ran2(&seed);
	for (i = 0; i < m_nsamp; i++) {
		points[i].pos[0] = min[0] + ran2(&seed) * (max[0]-min[0]);
		points[i].pos[1] = min[1] + ran2(&seed) * (max[1]-min[1]);
		points[i].pos[2] = min[2] + ran2(&seed) * (max[2]-min[2]);
		alive[i] = true;
		done[i] = false;
	}

	for (j = 0; j < 100; j++) {
		for (i = 0; i < m_nsamp; i++) {
			float val, *grad;
			if (done[i] || !alive[i]) continue;
			if (!p_fun->valid(points[i].pos)) {
				//alive[i] = false;
				//continue;
				points[i].pos[0] = min[0] + ran2(&seed) * (max[0]-min[0]);
				points[i].pos[1] = min[1] + ran2(&seed) * (max[1]-min[1]);
				points[i].pos[2] = min[2] + ran2(&seed) * (max[2]-min[2]);
			}
			grad = points[i].normal;
			p_fun->evalFuncGrad(points[i].pos, val, grad);
			//printf("point %d, (%f %f %f) = %f\n", i, points[i].x, points[i].y, points[i].z, val);
			if (fabs(val-isoval) < m_threshold) {
				done[i] = true;
				count++;
				continue;
			}
			float len = (grad[0]*grad[0]+grad[1]*grad[1]+grad[2]*grad[2]);
			float dif = val - isoval;
			points[i].pos[0] += -grad[0]*dif/len;
			points[i].pos[1] += -grad[1]*dif/len;
			points[i].pos[2] += -grad[2]*dif/len;
		}
		if (count == m_nsamp) break;
	}
	unsigned int t2 = time(NULL);
	printf("%d particles converged in %d secs\n", count, t2-t1);

	for (i = 0; i < m_nsamp; i++) {
		if (done[i]) {
			normalize3f(points[i].normal);
			p_pnts->insert(points[i]);
		}
	}

	delete[] done;
	delete[] alive; 
}

void PointSurface::tracePoint(Point3D &pnt)
{
	Point3D tmp = pnt;
	int i;
	float val, grad[3];

	for(i = 0; i < 100; i++) {
		p_fun->evalFuncGrad(pnt.pos, val, grad);
		if(fabs(val-m_isoval) < m_threshold) {
			return;			
		} else if(!p_fun->valid(pnt.pos)) {
			break;
		} else {
			float len = (grad[0]*grad[0]+grad[1]*grad[1]+grad[2]*grad[2]);
			float dif = val - m_isoval;
			pnt.pos[0] += -grad[0]*dif/len;
			pnt.pos[1] += -grad[1]*dif/len;
			pnt.pos[2] += -grad[2]*dif/len;
		}
	}
	pnt = tmp;
}

void PointSurface::upSampling(int vd1, int vd2, int vd3, int n, Point3D* samples)
{
	int i, j;
	float (*bary)[3];
	Point3D* points = p_pnts->data();

	bary = (float (*)[3]) malloc(sizeof(float[3]) * n *n);
	uniTriSample(n, bary);

	for(i = 0; i < n*n; i++) {
		for(j = 0; j < 3; j++) {
			samples[i].pos[j] = bary[i][0]*points[vd1].pos[j] + bary[i][1]*points[vd2].pos[j] + bary[i][2]*points[vd3].pos[j];
			samples[i].normal[j] = bary[i][0]*points[vd1].normal[j] + bary[i][1]*points[vd2].normal[j] + bary[i][2]*points[vd3].normal[j];
		}
		tracePoint(samples[i]);		
	}
	free(bary);
}

void PointSurface::repulse()
{
	int i, j, k, np = p_pnts->length();
	Point3D *points = p_pnts->data();
	int *nn_idx;
	float *dist;
	NNStruct* nnstrt = new NNBucket(np, points);

	float *rr2 = new float[np];
	vector3 *forces = new vector3[np];
	float min[3], max[3];

	p_fun->getBoundingBox(min, max);
	float avgr2 = (2*((max[2]-min[2])*(max[1]-min[1]) + (max[2]-min[2])*(max[0]-min[0]) + 
					  (max[1]-min[1])*(max[0]-min[0])) / np);
	for (i = 0; i < np; i++) {
		rr2[i] = avgr2;
		forces[i].x = forces[i].y = forces[i].z = 0;
	}

	int nn = 10;
	nn_idx = new int[nn];
	dist = new float[nn];
	for (i = 0; i < np; i++) {
		k = nnstrt->NNSearch(points[i], nn, rr2[i], nn_idx, dist);
		for (j = 0; j < k; j++) {
			if (nn_idx[j] != i && dotProduct3f(points[i].normal, points[nn_idx[j]].normal) > 0)
				forces[i] += repulsionForce(i, nn_idx[j], rr2[i], rr2[nn_idx[j]]);
		}
	}

	Point3D pnt;
	for (i = 0; i < np; i++) {
		float val, *grad;
		forces[i] = points[i].project(forces[i]);
		pnt.pos[0] = points[i].pos[0] + forces[i].x * (float)sqrt(rr2[i]) / 6;
		pnt.pos[1] = points[i].pos[1] + forces[i].y * (float)sqrt(rr2[i]) / 6;
		pnt.pos[2] = points[i].pos[2] + forces[i].z * (float)sqrt(rr2[i]) / 6;

		for (j = 0; j < 100; j++) {
			if (!p_fun->valid(pnt.pos)) {
				break;
			}
			grad = pnt.normal;
			p_fun->evalFuncGrad(pnt.pos, val, grad);
			if (fabs(val-m_isoval) < m_threshold) {
				points[i] = pnt;
				break;
			}
			float len = (grad[0]*grad[0]+grad[1]*grad[1]+grad[2]*grad[2]);
			float dif = val - m_isoval;
			pnt.pos[0] += -grad[0]*dif/len;
			pnt.pos[1] += -grad[1]*dif/len;
			pnt.pos[2] += -grad[2]*dif/len;
		}
	}

	delete[] nn_idx;
	delete[] dist;
	delete nnstrt;
	delete[] forces;
	delete[] rr2;
}

int PointSurface::insertPoint(float pos[3])
{
	float val, *grad;
	Point3D pnt(pos);

	grad = pnt.normal;
	for (int i = 0; i < 100; i++) {
		p_fun->evalFuncGrad(pnt.pos, val, grad);
		if (fabs(val-m_isoval) < m_threshold) {
			break;
		}
		float len = (grad[0]*grad[0]+grad[1]*grad[1]+grad[2]*grad[2]);
		float dif = val - m_isoval;
		pnt.pos[0] += -grad[0]*dif/len;
		pnt.pos[1] += -grad[1]*dif/len;
		pnt.pos[2] += -grad[2]*dif/len;
	}
	p_pnts->insert(pnt);

	return p_pnts->length()-1;
}

vector3 PointSurface::repulsionForce(int p1, int p2, float rr1, float rr2)
{
	vector3 force(0, 0, 0);

	if (p1 == p2) return force;
	float v12[3];

	for (int i = 0; i < 3; i++) {
		v12[i] =  (*p_pnts)[p1].pos[i] - (*p_pnts)[p2].pos[i];
	}
	float len = normalize3f(v12);
	float fmag = (float) exp(-3*len*len/rr1);

	force.x = v12[0] * fmag;
	force.y = v12[1] * fmag;
	force.z = v12[2] * fmag;

	return force;
}

