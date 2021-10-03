#include "moment.h"
#include <nr.h>
#include <math.h>

void MomtAttrib::print()
{
	printf("Inertia: %f %f %f, ", I1, I2, I3);
	printf("Function integral = %f, min = %f, max = %f, ", Fint, Fmin, Fmax);
	printf("Dipole moment len = %f, direction = %f, ", Dlen, Dang);
	printf("Quarduple moment = %f %f %f\n", Q1, Q2, Q3);
}

void VolMoments::principalAxes(float *eigval, float **eigvec)
{
	int i, j, nrot;
	float **mtx = matrix(1, 3, 1, 3);
	getCovMatrix(mtx);
	//double len = DMAX(1, (fabs(mtx[1][1])+fabs(mtx[2][2])+fabs(mtx[3][3]))/3);

	jacobi(mtx, 3, eigval, eigvec, &nrot);
	eigsrt(eigval, eigvec, 3);
	free_matrix(mtx, 1, 3, 1, 3);
	
	for (i = 1; i <= 3; i++) {
		printf("eigenval %d = %f\n", i, eigval[i]);
		printf("eigenvec = (%f %f %f) \n", eigvec[1][i], eigvec[2][i], eigvec[3][i]);
	}
}

MomtAttrib VolMoments::toAttributes()
{
	MomtAttrib attrib;
	
	if(vol <= 0) return attrib;
	
	float ctr[3];
	getCOM(ctr);
	float *Fm1 = fvector(1, 3);
	float **Fm2 = matrix(1, 3, 1, 3);
	Fm1[1] = Fx - Fint*ctr[0];
	Fm1[2] = Fy - Fint*ctr[1];
	Fm1[3] = Fz - Fint*ctr[2];

	float **mtx = matrix(1, 3, 1, 3);
	getCovMatrix(mtx);
	float *eigval = fvector(1, 3);
	float **eigvec = matrix(1, 3, 1, 3);
	int nrot;
	jacobi(mtx, 3, eigval, eigvec, &nrot);
	eigsrt(eigval, eigvec, 3);
	free_matrix(mtx, 1, 3, 1, 3);
	
	//printf("Function integral = %f, min = %f, max = %f\n", Fint, Fmin, Fmax);
	float len = sqrt(Fm1[1]*Fm1[1] + Fm1[2]*Fm1[2] + Fm1[3]*Fm1[3]);
	float angle = (Fm1[1]*eigvec[1][1] + Fm1[2]*eigvec[2][1] + Fm1[3]*eigvec[3][1]) / len;
	//printf("Func 1st moment len = %f, direction = %f\n", len, angle); 
	Fm2[1][1] = Fxx - 2*Fx*ctr[0] + ctr[0]*ctr[0]*Fint;
	Fm2[2][2] = Fyy - 2*Fy*ctr[1] + ctr[1]*ctr[1]*Fint;
	Fm2[3][3] = Fzz - 2*Fz*ctr[2] + ctr[2]*ctr[2]*Fint;
	Fm2[1][2] = Fm2[2][1] = Fxy - Fx*ctr[1] - Fy*ctr[0] + ctr[0]*ctr[1]*Fint;
	Fm2[1][3] = Fm2[3][1] = Fxz - Fx*ctr[2] - Fz*ctr[0] + ctr[0]*ctr[2]*Fint;
	Fm2[2][3] = Fm2[3][2] = Fyz - Fy*ctr[2] - Fz*ctr[1] + ctr[1]*ctr[2]*Fint;
	float *peigval = fvector(1, 3);
	float **peigvec = matrix(1, 3, 1, 3);
	jacobi(Fm2, 3, peigval, peigvec, &nrot);
	eigsrt(peigval, peigvec, 3);
	//printf("potential eigenvals = %f %f %f\n", peigval[1], peigval[2], peigval[3]);
	
	attrib.I1 = eigval[1] / vol;
	attrib.I2 = eigval[2] / vol;
	attrib.I3 = eigval[3] / vol;
	attrib.Fint = Fint / vol;
	attrib.Fmin = Fmin;
	attrib.Fmax = Fmax;
	attrib.Dlen = len / vol;
	attrib.Dang = angle;
	attrib.Q1 = peigval[1] / vol;
	attrib.Q2 = peigval[2] / vol;
	attrib.Q3 = peigval[3] / vol;
	
	free_vector(eigval, 1, 3);
	free_matrix(eigvec, 1, 3, 1, 3);
	free_matrix(peigvec, 1, 3, 1, 3);
	free_vector(peigval, 1, 3);	
	free_matrix(Fm2, 1, 3, 1, 3);
	free_vector(Fm1, 1, 3);

	return attrib;
}
