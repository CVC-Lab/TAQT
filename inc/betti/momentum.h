#ifndef VOL_MOMENTUM_H
#define VOL_MOMENTUM_H

#include <stdio.h>

/************************************************************************/
/* The first and second order momentums of a volume                     */
/************************************************************************/

class VolMomentum {
public:
	VolMomentum() {
		vol = 0;
		Px = Py = Pz = 0;
		Pxx = Pyy = Pzz = Pxy = Pxz = Pyz = 0;
	}

	/**
	 * Get the center of mass of the volume
	 */
	inline void getCOM(float com[3]);

	/**
	 * Get the covariance matrix 
	 */
	inline void getCovMatrix(float *mtx);

	/**
	 * Print out in a nice format
	 */
	inline void print();

	inline void printRaw();

	VolMomentum& operator = (const VolMomentum & mom) {
		if (this == &mom) {
			return *this;
		}
		vol = mom.vol;
		Px = mom.Px; Py = mom.Py; Pz = mom.Pz;
		Pxx = mom.Pxx; Pyy = mom.Pyy; Pzz = mom.Pzz;
		Pxy = mom.Pxy; Pxz = mom.Pxz; Pyz = mom.Pyz;
	}

	friend VolMomentum operator + (const VolMomentum & m1, const VolMomentum & m2) {
		VolMomentum mom;
		mom.vol = m1.vol + m2.vol;
		mom.Px = m1.Px + m2.Px;
		mom.Py = m1.Py + m2.Py;
		mom.Pz = m1.Pz + m2.Pz;
		mom.Pxx = m1.Pxx + m2.Pxx;
		mom.Pyy = m1.Pyy + m2.Pyy;
		mom.Pzz = m1.Pzz + m2.Pzz;
		mom.Pxy = m1.Pxy + m2.Pxy;
		mom.Pxz = m1.Pxz + m2.Pxz;
		mom.Pyz = m1.Pyz + m2.Pyz;
		return mom;
	}
	/************************************************************************/
	/* Member variables                                                     */
	/************************************************************************/
	float vol;
	float Px, Py, Pz;
	float Pxx, Pyy, Pzz, Pxy, Pxz, Pyz;
}; 

void VolMomentum::getCOM(float com[3])
{
	if (vol <= 0) {
		com[0] = com[1] = com[2];
		return;
	}
	com[0] = Px / vol;
	com[1] = Py / vol;
	com[2] = Pz / vol;
}

void VolMomentum::getCovMatrix(float *mtx)
{
	int i;
	if (vol <= 0) {
		for (i = 0; i < 9; i++) {
			mtx[i] = 0;
		}
		return;
	}
	float x[3];
	getCOM(x);
	mtx[0] = Pxx - vol*x[0]*x[0];
	mtx[4] = Pyy - vol*x[1]*x[1];
	mtx[8] = Pzz - vol*x[2]*x[2];
	mtx[1] = mtx[3] = Pxy - vol*x[0]*x[1];
	mtx[2] = mtx[6] = Pxz - vol*x[0]*x[2];
	mtx[5] = mtx[7] = Pyz - vol*x[1]*x[2];
}

void VolMomentum::print()
{
	float ctr[3], mtx[9];
	getCOM(ctr);
	getCovMatrix(mtx);
	printf("Volume = %f\n", vol);
	printf("Center: (%f %f %f)\n", ctr[0], ctr[1], ctr[2]);
	for (int i = 0; i < 3; i++) {
		printf("%f\t %f\t %f\n", mtx[3*i], mtx[3*i+1], mtx[3*i+2]);
	}
}

void VolMomentum::printRaw()
{
	printf("Volume = %f\n", vol);
	printf("Px = %f, Py = %f, Pz = %f\n", Px, Py, Pz);

	printf("%f\t%f\t%f\n", Pxx, Pxy, Pxz);
	printf("%f\t%f\t%f\n", Pxy, Pyy, Pyz);
	printf("%f\t%f\t%f\n", Pxz, Pyz, Pzz);
}
#endif 

		  
