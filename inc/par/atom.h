#ifndef _XYZ_ATOM_H
#define _XYZ_ATOM_H 

/**
 * Atom -- A Class representing an atom
 */
class Atom {
public:
	Atom(int _id, float _center[3]);

	Atom(int _id = 0, float x = 0, float y = 0, float z = 0);

	void setBlobby(float _blob) {
		blobby = _blob;
	}
	
	float getBlobby() {
		return blobby;
	}

	/**
	 * Evaluate the electron density at pnt using Gaussian distribution 	
	 */
	double elecDenGaussian(const float pnt[3]);

	void elecGradGaussian(const float pnt[3], float grad[3]);

	void elecDenGradGaussian(const float pnt[3], float& den, float grad[3]);

	/*
	 *	Evaluate the electron density using metaball approximation
	 */
	double elecDenMeta(const float pnt[3]);

	void elecGradMeta(const float pnt[3], float grad[3]);

	void elecDenGradMeta(const float pnt[3], float& den, float grad[3]);

	/*
	 *	Add contributions from this atom to the input grid function
	 */
	void accumDensity(float *data, int dim[3], float orig[3], float span[3]);

	void accumGradient(float *grads, int dim[3], float orig[3], float span[3]);
	
	void accumDenAndGrad(float *data, float *grads, int dim[3], float orig[3], float span[3]);
	
	int id;					// id number in the periodical table
	float center[3];		// coordinates of the atom center

	/*
	 *	Atom shell distribution
	 *  f = 0: r > Ro
	 *  f > 0: Ri < r < Ro
	 *  f < 0: r < Ri
	 */
	float evalShell(const float pnt[3]);

	void accumShell(float *data, int dim[3], float orig[3], float span[3]);

	void setThickness(float _thick) {
		thickness = _thick;
	}

	float getThickness() const {
		return thickness;
	}

protected:
	void init();
	
	float blobby;			// blobby factor

	static float a, b;		// parameters used in metaball approximation			

	float thickness;
};

#endif

