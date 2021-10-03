#ifndef VOL_MOMENTUM_H
#define VOL_MOMENTUM_H

#include <math.h>
#include <stdio.h>
#include <nrutil.h>

#include <iostream>

using namespace std;

/**
 * Attributes related to the volume moments
 */
class MomtAttrib{
public:
	float I1, I2, I3;	// inertia
	float Fint;
	float Fmin, Fmax;	// integral of f over the volume
	float Dlen;			// the length of the dipole moment
	float Dang;			// Cos of the angle between the dipole moment and the principal axis
	float Q1, Q2, Q3;	// major components of the quardruple moment
	
	void print();
};

/************************************************************************/
/* The first and second order moments of a volume                       */
/************************************************************************/

class VolMoments {
public:
	VolMoments() {
		vol = 0;
		Px = Py = Pz = 0;
		Pxx = Pyy = Pzz = Pxy = Pxz = Pyz = 0;
		Fint = 0;
		Fmin = 0; Fmax = 0;
		Fx = 0; Fy = 0; Fz = 0;
		Fxx = Fyy = Fzz = Fxy = Fxz = Fyz = 0;
	}

	/**
	 * Get the center of mass of the volume
	 */
	inline void getCOM(float com[3]);

	/**
	 * Get the 3x3 covariance matrix 
	 */
	inline void getCovMatrix(float **mtx);

	/**
	 * Print out in a nice format
	 */
	inline void print();

	inline void printRaw();

	/**
	 * Convert to matchable attributes
	 */
	MomtAttrib toAttributes();
	
	/**
	 * Compute the eigenvalues and eigenvectors of the covariance matrix
	 * The eigenvalues are sorted by their values in descending order
	 */
	void principalAxes(float *eigval, float **eigvec);
	
	VolMoments& operator = (const VolMoments & mom) {
		if (this == &mom) {
			return *this;
		}
		vol = mom.vol;
		Px = mom.Px; Py = mom.Py; Pz = mom.Pz;
		Pxx = mom.Pxx; Pyy = mom.Pyy; Pzz = mom.Pzz;
		Pxy = mom.Pxy; Pxz = mom.Pxz; Pyz = mom.Pyz;
		Fint = mom.Fint; 
		Fmin = mom.Fmin; Fmax = mom.Fmax;
		Fx = mom.Fx; Fy = mom.Fy; Fz = mom.Fz;
		Fxx = mom.Fxx; Fyy = mom.Fyy; Fzz = mom.Fzz;
		Fxy = mom.Fxy; Fxz = mom.Fxz; Fyz = mom.Fyz;
		return *this;
	}

	friend VolMoments operator + (const VolMoments & m1, const VolMoments & m2) {
		VolMoments mom;
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
		mom.Fint = m1.Fint + m2.Fint;
		mom.Fx = m1.Fx + m2.Fx;
		mom.Fy = m1.Fy + m2.Fy;
		mom.Fz = m1.Fz + m2.Fz;
		mom.Fmin = FMIN(m1.Fmin, m2.Fmin);
		mom.Fmax = FMAX(m1.Fmax, m2.Fmax);
		mom.Fxx = m1.Fxx + m2.Fxx;
		mom.Fyy = m1.Fyy + m2.Fyy;
		mom.Fzz = m1.Fzz + m2.Fzz;
		mom.Fxy = m1.Fxy + m2.Fxy;
		mom.Fxz = m1.Fxz + m2.Fxz;
		mom.Fyz = m1.Fyz + m2.Fyz;
		return mom;
	}
	
	friend ostream& operator << (ostream& out, const VolMoments& mom) {
		out << mom.vol << " " << mom.Px << " " << mom.Py << " " << mom.Pz << " ";
		out << mom.Pxx << " " << mom.Pyy << " " << mom.Pzz << " ";
		out << mom.Pxy << " " << mom.Pxz << " " << mom.Pyz << " ";
		out << mom.Fint << " " << mom.Fmin << " " << mom.Fmax << " ";
		out << mom.Fx << " " << mom.Fy << " " << mom.Fz << " ";
		out << mom.Fxx << " " << mom.Fyy << " " << mom.Fzz << " ";
		out << mom.Fxy << " " << mom.Fxz << " " << mom.Fyz << endl;
		return out;
	}
	
	friend istream& operator >> (istream& in, VolMoments& mom) {
		in >> mom.vol >> mom.Px >> mom.Py >> mom.Pz;
		in >> mom.Pxx >> mom.Pyy >> mom.Pzz;
		in >> mom.Pxy >> mom.Pxz >> mom.Pyz;
		in >> mom.Fint >> mom.Fmin >> mom.Fmax;
		in >> mom.Fx >> mom.Fy >> mom.Fz;
		in >> mom.Fxx >> mom.Fyy >> mom.Fzz;
		in >> mom.Fxy >> mom.Fxz >> mom.Fyz;
		return in;
	}
		
	/************************************************************************/
	/* Member variables                                                     */
	/************************************************************************/
	float vol;
	float Px, Py, Pz;
	float Pxx, Pyy, Pzz, Pxy, Pxz, Pyz;
	/* Moments of the potential f*/
	float Fint;			// integral of f over the volume
	float Fmin, Fmax;
	float Fx, Fy, Fz;
	float Fxx, Fyy, Fzz, Fxy, Fxz, Fyz;
}; 

void VolMoments::getCOM(float com[3])
{
	if (vol <= 0) {
		com[0] = com[1] = com[2];
		return;
	}
	com[0] = Px / vol;
	com[1] = Py / vol;
	com[2] = Pz / vol;
}

void VolMoments::getCovMatrix(float **mtx)
{
	int i, j;
	if (vol <= 0) {
		for (i = 0; i < 3; i++) {
			for(j = 0; j < 3; j++) {
				mtx[i][j] = 0;
			}
		}
		return;
	}
	float x[3];
	getCOM(x);
	mtx[1][1] = Pxx - vol*x[0]*x[0];
	mtx[2][2] = Pyy - vol*x[1]*x[1];
	mtx[3][3] = Pzz - vol*x[2]*x[2];
	mtx[1][2] = mtx[2][1] = Pxy - vol*x[0]*x[1];
	mtx[1][3] = mtx[3][1] = Pxz - vol*x[0]*x[2];
	mtx[2][3] = mtx[3][2] = Pyz - vol*x[1]*x[2];
}

void VolMoments::print()
{
	float ctr[3], **mtx;
	mtx = matrix(1, 3, 1, 3);
	getCOM(ctr);
	getCovMatrix(mtx);
	printf("Volume = %f\n", vol);
	printf("Center: (%f %f %f)\n", ctr[0], ctr[1], ctr[2]);
	for (int i = 1; i <= 3; i++) {
		printf("%f\t %f\t %f\n", mtx[i][1], mtx[i][2], mtx[i][3]);
	}
	
	free_matrix(mtx, 1, 3, 1, 3);
}

void VolMoments::printRaw()
{
	printf("Volume = %f\n", vol);
	printf("Px = %f, Py = %f, Pz = %f\n", Px, Py, Pz);

	printf("%f\t%f\t%f\n", Pxx, Pxy, Pxz);
	printf("%f\t%f\t%f\n", Pxy, Pyy, Pyz);
	printf("%f\t%f\t%f\n", Pxz, Pyz, Pzz);
}
#endif 

		  
