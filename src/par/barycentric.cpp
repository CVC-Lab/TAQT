#include "barycentric.h"

void uniTriSample(int n, float (*bary)[3])
{
	int i, j, k;
	float denom = 3.0f * n;

	for(i = 0; i < n; i++) {
		k = i*i;
		for(j = 0; j < 2*i+1; j++, k++) {
			bary[k][0] = 1 - (3*i+2)/denom;
			bary[k][1] = (((j+3)>>1) + ((j >> 1) << 1)) / denom;
			bary[k][2] = 1 - bary[k][0] - bary[k][1];
		}
	}
}

