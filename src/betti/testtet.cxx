#include <stdio.h>
#include <stdlib.h>
#include "volume.h"

int main(int argc, char* argv[])
{
	float x1[3] = {0, 0, 0};
	float x2[3] = {0, 1, 0};
	float x3[3] = {0, 1, 1};
	float x4[3] = {1, 1, 1};
	float v1 = 1, v2 = 2, v3 = 3, v4 = 4;
	float f1 = -1, f2 = -1, f3 = -1, f4 = -2;
	
	float x = atof(argv[1]);
	
	float vol = tetVolume(x1, x2, x3, x4, v1, v2, v3, v4, x);
	float fint = tetraFuncIntegral(x1, x2, x3, x4, v1, v2, v3, v4, f1, f2, f3, f4, x);
	
	printf("vol = %f, fint = %f\n", vol, fint);
	
	return 0;
}

