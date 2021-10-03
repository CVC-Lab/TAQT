#include <math.h>
#include <stdlib.h>

#include "vecorder.h"

static vector3 xaxis, yaxis;
static const float CTAN_MAX = 1e10;
static const float CTAN_MIN = -1e10;

static int angleCompare(const void* p1, const void* p2)
{
	vector3 *v[2];
	v[0] = (vector3 *)p1;
	v[1] = (vector3 *)p2;

	int sgn[2], i;
	float ctan[2];
	
	for(i = 0; i < 2; i++) {
		float vy = DotProduct(*v[i], yaxis);
		if(vy == 0) {
			sgn[i] = 1;
			ctan[i] = (DotProduct(*v[i], xaxis) > 0)? CTAN_MAX:CTAN_MIN;
		} else if(vy < 0) {
			sgn[i] = -1;
			ctan[i] = DotProduct(*v[i], xaxis) / DotProduct(*v[i], yaxis);
		} else {
			sgn[i] = 1;
			ctan[i] = DotProduct(*v[i], xaxis) / DotProduct(*v[i], yaxis);
		}
	}
	
	if(sgn[0] > sgn[1]) return -1;
	else if(sgn[0] < sgn[1]) return 1;
	else if(ctan[0] > ctan[1]) return -1;
	else if(ctan[0] < ctan[1]) return 1;

	return 0;
}

bool operator < (const PlaneAngle& a1, const PlaneAngle& a2)
{
	if(a1.m_sgny < a2.m_sgny) return false;
	else if(a1.m_sgny > a2.m_sgny) return true;
	return (a1.m_ctan > a2.m_ctan);
}

void sortVectorCC(vector3 vecs[], int nv, const vector3 normal)
{
	// use vecs[0] as the x axix;
	yaxis = CrossProduct(normal, vecs[0]);
	xaxis = vecs[0];
	xaxis.normalize();
	yaxis.normalize();
	
	qsort(&(vecs[1]), nv-1, sizeof(vector3), angleCompare);
}

