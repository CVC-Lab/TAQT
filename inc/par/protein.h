#ifndef PDB_PROTEIN_H
#define PDB_PROTEIN_H

#include "atom.h"
#include "dynarray.h"

/**
 * Protein -- A Class representing a protein molecule.
 */
class Protein {
public:
	// Construct from a PDB file
	Protein(const char* fname);

	///
	virtual ~Protein();

	/**
	 * Get the number of atoms in the molecule
	 */
	int numOfAtoms() const { return atoms.length(); }

	/**
	 * Get the nth atom in the molecule
	 */
	Atom getAtom(int n) const;

	/**
	 * Add an atom to the molecule
	 */
	void addAtom(const Atom& atom);

	/**
	 * Get bounding box of the molecule
	 */
	void getBoundingBox(float min[3], float max[3]);

	/**
	 * Compute the ElectronDensity 
	 * @return in Function3D format for topological compuation
	 */
	void getElectronDensity(int dim[3], float *dens);

	/**
	 * Compute the gridient of Electron Density.
	 */
	void getGradient(int dim[3], float *grads);

	/*
	 *	Compute the electron density and gradient.
	 */
	void getElecDenAndGradient(int dim[3], float *dens, float *grads);

	/**
	 * Evaluate density and gradient at a given point in space
	 */
	void evalDenAndGrad(const float pnt[3], float& den, float grad[3]);

	float evalDensity(const float pnt[3]);

	void  evalGradient(const float pnt[3], float grad[3]);
	
	float evalLaplace(const float pnt[3]);


	/*
	 *	Compute the molecular surface shell and interior representation
	 */
	void getMolShell(int dim[3], float *data);

	void getMolInter(int dim[3], float *data);

	/*
	 * Set the blobby factor of all atoms.
	 */
	 void setBlobby(float _blobby) {
		 blobby = _blobby;
		 for(int i = 0; i < atoms.length(); i++) {
			 atoms[i].setBlobby(_blobby);
		 }
	 }
	 
protected:
	//
	void readPDB(const char* fname);
	
	/**
	 * Compute Electron Potential value at mesh points
	 * @param min Return the minimum of the potential
	 * @param max Return the maximum of the potential
	 */
	void calcElectronDensity(float* data, int dim[3], float orig[3], 
							 float span[3], float& min, float& max);

	void calcGradient(float* data, int dim[3], float orig[3], 
							 float span[3], float min[3], float max[3]);

	void calcLaplace(float* data, int dim[3], float orig[3], 
							float span[3], float& min, float& max);


	static float blobby;
	dynamic_array<Atom> atoms;			// array of atoms of the molecule	
	float m_min[3], m_max[3];			// bounding box
};

#endif


