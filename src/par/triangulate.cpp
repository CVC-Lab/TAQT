#include "triangulate.h"
#include "mtxlib.h"
#include "vecmath.h"
#include "nnbucket.h"
#include "nnkdtree.h"

int Triangulate::nn = 6;

int comparePolar(const void *p1, const void *p2)
{
	Polar2Coord polar1 = *((Polar2Coord *) p1);
	Polar2Coord polar2 = *((Polar2Coord *) p2);

	if(polar1 < polar2) return -1;
	else if(polar2 < polar1) return 1;
	return 0;
}

void sortPolarCC(int n, Polar2Coord* polars) 
{
	qsort(polars, n, sizeof(Polar2Coord), comparePolar);
}

/*
 * Check if (r1, theta1) is in the triangle formed by (r0, 0) and (r2, theta2)
 * @note (r, theta) here is the polar coordinates
 */
bool inTriangle(float r0, float r1, float theta1, float r2, float theta2)
{
	return ((r2*sin(theta2)*r1*cos(theta1) + (r0-r2)*cos(theta2)*r1*sin(theta1)
			- r0*r2*sin(theta2)) < 0); 
}

Triangulate::Triangulate(PointSurface* psurf)
: p_pntsrf(psurf) 
{
	int np = psurf->getNumOfPoints();
	p_cnct = new dynamic_array<PointConnect *>((int)(np*1.1));
	for(int i = 0; i < np; i++) {
		PointConnect *pcnt = new PointConnect;
		p_cnct->insert(pcnt);
	}
}

Triangulate::~Triangulate()
{
	int np = p_cnct->length();
	for(int i = 0; i < np; i++) {
		delete (*p_cnct)[i];
	}
	delete p_cnct;
}

void Triangulate::NNTriangulate()
{
	int i, j, l, np;
	Point3D* points = p_pntsrf->getPoints(np);
	NNBucket *nnbuck = new NNBucket(np, points);
	int *nnidx = new int[nn+1];
	float *dist = new float[nn+1];
	vector3 *vecs = new vector3[nn+1];

	Polar2Coord *polars = new Polar2Coord[nn];
	for(i = 0; i < p_pntsrf->getNumOfPoints(); i++) {
		int k = nnbuck->NNSearch(points[i], nn+1, nnidx, dist);
		for(j = 0, l = 0; j < k; j++ ) {
			if(nnidx[j] != i && dotProduct3f(points[i].normal, points[nnidx[j]].normal) > 0) {
				polars[l] = points[i].project2dPolar(points[nnidx[j]]);
				polars[l].id = nnidx[j];
				l++;
			}
		}
		sortPolarCC(l, polars);
		for(j = 0; j < l; j++) 
			(*p_cnct)[i]->insert(polars[j].id);
	}

	delete[] polars;
	delete[] nnidx;
	delete[] dist;
	delete nnbuck;
}

void Triangulate::NNTriangulate(float r)
{
	int i, j, k, m, nn, np;
	Point3D* points = p_pntsrf->getPoints(np);
	NNStruct *nnstruct = new NNKDTree(np, points);
	float r2 = r*r;

	nn = 15;
	int *nn_idx = new int[nn];
	float *dist = new float[nn];
	Polar2Coord* polars = new Polar2Coord[nn+1];
	for(i = 0; i < np; i++) {
		m = nnstruct->NNSearch(points[i], nn, r2, nn_idx, dist);
		for(j = 0, k = 0; j < m; j++) {
			if(nn_idx[j] != i && dotProduct3f(points[i].normal, points[nn_idx[j]].normal) > 0) {
				polars[k] = points[i].project2dPolar(points[nn_idx[j]]);
				polars[k].id = nn_idx[j];
				k++;
			}
		}
		sortPolarCC(k, polars);
		// find the neighbor that is closest to the point
		int		t = 0; 
		float rds = polars[t].r;
		for(j = 1; j < k; j++) {
			if(polars[j].r < rds) {
				t = j;
				rds = polars[j].r;
			}
		}

		(*p_cnct)[i]->insert(polars[t].id);
		float oldtheta = polars[t].theta;
		float oldr = polars[t].r;
		for(j = 0; j <= t; j++) {
			polars[j].theta += 2*PI;
		}
		for(j = 1; j < k; j++) {
			int n1 = (j+t) % k;
			int n2 = (j+t+1) % k;
			if(polars[n2].theta - oldtheta > 2) {
				(*p_cnct)[i]->insert(polars[n1].id);
				oldtheta = polars[n1].theta;
				oldr = polars[n1].r;
			} else if(inTriangle(oldr, polars[n1].r, polars[n1].theta-oldtheta, 
								 polars[n2].r, polars[n2].theta-oldtheta)) {
				(*p_cnct)[i]->insert(polars[n1].id);
				oldtheta = polars[n1].theta;
				oldr = polars[n1].r;
			}
		}
	}

	delete[] polars;
	delete[] nn_idx;
	delete[] dist;
	delete nnstruct;
}

void Triangulate::NNTriangulate2(float r)
{
	int i, j, k, np;

	Point3D* points = p_pntsrf->getPoints(np);
	NNBucket *nnbuck = new NNBucket(np, points);
	dynamic_array<int> nn_idx;
	Polar2Coord *polars;

	float r2 = r*r;
	for(i = 0; i < p_pntsrf->getNumOfPoints(); i++) {
		int m = nnbuck->NNSearch(points[i], r2, nn_idx);
		polars = new Polar2Coord[m+1];
		for(j = 0, k = 0; j < m; j++) {
			if(nn_idx[j] != i && dotProduct3f(points[i].normal, points[nn_idx[j]].normal) > 0){
				polars[k] = points[i].project2dPolar(points[nn_idx[j]]);
				polars[k].id = nn_idx[j];
				k++;
			}
		}

		sortPolarCC(k, polars);
		(*p_cnct)[i]->insert(polars[0].id);
		polars[k].id = polars[0].id;
		polars[k].theta = 2*PI + polars[0].theta;
		polars[k].r	= polars[0].r;
		float oldtheta = polars[0].theta;
		for(j = 1; j < k; j++) {
			if(polars[j].theta - oldtheta > 1) {
				while(polars[j].theta - oldtheta > 1.57) {
					// add new points to the surface
					float dx = (float)(polars[j].r*(cos(oldtheta+1.03) - cos(polars[j].theta)));
					float dy = (float)(polars[j].r*(sin(oldtheta+1.03) - sin(polars[j].theta)));
					float xax[3], yax[3], pos[3];
					points[i].localAxes(xax, yax);
					pos[0] = points[polars[j].id].pos[0] + dx*xax[0] + dy*yax[0];
					pos[1] = points[polars[j].id].pos[1] + dx*xax[1] + dy*yax[1];
					pos[2] = points[polars[j].id].pos[2] + dx*xax[2] + dy*yax[2];
					int pid = p_pntsrf->insertPoint(pos);
					PointConnect *pc = new PointConnect;
					p_cnct->insert(pc);
					(*p_cnct)[i]->insert(pid);
					points = p_pntsrf->getPoints(np);
					nnbuck->insertAndReplace(np, points, pid);
					Polar2Coord newpnt = points[i].project2dPolar(points[pid]);
					newpnt.id = pid;
					oldtheta = newpnt.theta;
				}
				(*p_cnct)[i]->insert(polars[j].id);
				oldtheta = polars[j].theta;
			} else if(polars[j+1].theta - oldtheta > 1.57) {
				(*p_cnct)[i]->insert(polars[j].id);
				oldtheta = polars[j].theta;
			}
		}

		delete[] polars;
	}			

	delete nnbuck;
}



