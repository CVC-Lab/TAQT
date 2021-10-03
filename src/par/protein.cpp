#include <string.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include <bio/bio.h>

#include "protein.h"
#include "elements.h"
#include "reg3data.h"

float Protein::blobby = -4;

Protein::Protein(const char* fname)
{
	int len = strlen(fname);
	if (!(len > 4 && strcmp(fname+len-4, ".pdb") == 0)) {
		error("UnKown Molecule File");
	}
	readPDB(fname);
}

Protein::~Protein()
{
}

Atom Protein::getAtom(int n) const
{
	assert(n < numOfAtoms()); 
	return atoms[n];
}

void Protein::addAtom(const Atom& atom)
{
	atoms.insert(atom);
}

void Protein::getBoundingBox(float min[3], float max[3])
{
	for (int i = 0; i < 3; i++) {
		min[i] = m_min[i];
		max[i] = m_max[i];
	}
}

void Protein::getElectronDensity(int dim[3], float *vals)
{
	float min[3], max[3];
	getBoundingBox(min, max);

	float fmin, fmax;
	float orig[3], span[3];

	orig[0] = min[0]; orig[1] = min[1]; orig[2] = min[2];
	span[0] = (max[0]-min[0])/(dim[0]-1); 
	span[1] = (max[1]-min[1])/(dim[1]-1); 
	span[2] = (max[2]-min[2])/(dim[2]-1);

	// electron density
	calcElectronDensity(vals, dim, orig, span, fmin, fmax);
	printf("min value = %f, max value = %f\n", fmin, fmax); 
}

void Protein::getGradient(int dim[3], float *grads)
{
	float min[3], max[3];
	getBoundingBox(min, max);

	float fmin[3], fmax[3];
	float orig[3], span[3];

	orig[0] = min[0]; orig[1] = min[1]; orig[2] = min[2];
	span[0] = (max[0]-min[0])/(dim[0]-1); 
	span[1] = (max[1]-min[1])/(dim[1]-1); 
	span[2] = (max[2]-min[2])/(dim[2]-1);

	calcGradient(grads, dim, orig, span, fmin, fmax);
}

void Protein::getElecDenAndGradient(int dim[3], float *dens, float *grads)
{
	float min[3], max[3];
	getBoundingBox(min, max);
	
	float orig[3], span[3];
	
	orig[0] = min[0]; orig[1] = min[1]; orig[2] = min[2];
	span[0] = (max[0]-min[0])/(dim[0]-1); 
	span[1] = (max[1]-min[1])/(dim[1]-1); 
	span[2] = (max[2]-min[2])/(dim[2]-1);

	int i, nv;
	nv = dim[0] * dim[1] * dim[2];
	for (i = 0; i < nv; i++) {
		dens[i] = 0;
		grads[3*i] = 0;
		grads[3*i+1] = 0;
		grads[3*i+2] = 0;
	}
	for (i = 0; i < numOfAtoms(); i++) {
		atoms[i].accumDenAndGrad(dens, grads, dim, orig, span);
	}
}

void Protein::calcElectronDensity(float *data, int dim[3], float orig[3], 
								  float span[3], float &min, float &max)
{
	int i, nv;

	nv = dim[0]*dim[1]*dim[2];

	for (i = 0; i < nv; i++) data[i] = 0;
	for (i = 0; i < numOfAtoms(); i++) {
		atoms[i].accumDensity(data, dim, orig, span);
	}

/*
	int j, k;
	for (k = 0; k < dim[2]; k++)
		for (j = 0; j < dim[1]; j++)
			for (i = 0; i < dim[0]; i++) {
				int l = i + j*dim[0] + k*dim[0]*dim[1];
				float pnt[3];
				pnt[0] = orig[0] + i*span[0];
				pnt[1] = orig[1] + j*span[1];
				pnt[2] = orig[2] + k*span[2];
				data[l] = evalDensity(pnt);
			}*/

	min = max = data[0];
	for (i = 1; i < nv; i++) {
		min = (data[i] < min)? data[i]:min;
		max = (data[i] > max)? data[i]:max;
	}
}

void Protein::readPDB(const char* fname)
{
	FILE* fp;

	if ((fp = fopen(fname, "r")) == NULL) {
		error("Cannot Open the PDB File");
	}

	char    line[82], s[5], id[4];
	int     i, j;
	double  x, y, z;

	while (fgets( line, sizeof(line), fp ) != NULL ) {
		if ( strncmp( line, "ATOM", 4 ) == 0 ) {
			sscanf(&line[6],"%5d", &i );
			sscanf(&line[13],"%3s", id );
			sscanf(&line[17],"%4s", s );
			sscanf(&line[30],"%8lf", &x );
			sscanf(&line[38],"%8lf", &y );
			sscanf(&line[46],"%8lf", &z );
			id[0] = toupper( id[0] );
			id[1] = tolower( id[1] );
			if ( id[1] == 0 ) {
				id[1] = ' ';
			}
			if (id[1] >= '0' && id[1] <= '9') {
				id[1] = ' ';
			}
			if ( s[0] != 0 ) {
				id[1] = ' ';
			}
			for ( i=0; i < MAXELEMENTS; i++ ) {
				if ((id[0] == element[i].symbol[0]) && 
					(id[1] == element[i].symbol[1]) ) {
					break;
				}
			}
			if ( i == MAXELEMENTS ) {
				error("Element in PDF file not identified");
			} else {
				Atom atom(i, (float)x, (float)y, (float)z);
				addAtom(atom);
			}
		}
	}       
	fclose(fp);

	if (numOfAtoms() == 0) {
		m_min[0] = m_min[1] = m_min[2] = 0;
		m_max[0] = m_max[1] = m_max[2] = 1;
		return;
	}
	for (i = 0; i < 3; i++) {
		m_min[i] = atoms[0].center[i];
		m_max[i] = atoms[0].center[i];
	}
	for (j = 1; j < atoms.length(); j++) {
		for (i = 0; i < 3; i++) {
			if (atoms[j].center[i] < m_min[i]) m_min[i] = atoms[j].center[i];
			if (atoms[j].center[i] > m_max[i]) m_max[i] = atoms[j].center[i];
		}
	}
	for (i = 0; i < 3; i++) {
		m_min[i] -= 4;
		m_max[i] += 4;
	}
}

void Protein::calcGradient(float* data, int dim[3], float orig[3], float span[3], float min[3], float max[3])
{
	int i, nv;

	nv = dim[0]*dim[1]*dim[2];
	for (i = 0; i < 3*nv; i++) {
		data[i] = 0;
	}

	for (i = 0; i < numOfAtoms(); i++) {
		atoms[i].accumGradient(data, dim, orig, span);
	}
/*
	int j, k;
  	float grad[3];
	for (k = 0; k < dim[2]; k++)
		for (j = 0; j < dim[1]; j++)
			for (i = 0; i < dim[0]; i++) {
				int l = i + j*dim[0] + k*dim[0]*dim[1];
				float pnt[3];
				pnt[0] = orig[0] + i*span[0];
				pnt[1] = orig[1] + j*span[1];
				pnt[2] = orig[2] + k*span[2];
				evalGradient(pnt, grad);
				data[0+3*l] = grad[0];
				data[1+3*l] = grad[1];
				data[2+3*l] = grad[2];
			}
	for (i = 0; i < 3; i++) {
		min[i] = data[i];
		max[i] = data[i];
	}
	for (i = 1; i < nv; i++) {
		for (j = 0; j < 3; j++) {
			min[j] = (data[j+3*i] < min[j])? data[j+3*i]:min[j];
			max[j] = (data[j+3*i] > max[j])? data[j+3*i]:max[j];
		}
	}*/

}

void Protein::calcLaplace(float* data, int dim[3], float orig[3], float span[3], float& min, float& max)
{
	int i, j, k, nv;

	nv = dim[0]*dim[1]*dim[2];

	for (i = 0; i < nv; i++) data[i] = 0;
	for (k = 0; k < dim[2]; k++)
		for (j = 0; j < dim[1]; j++)
			for (i = 0; i < dim[0]; i++) {
				int l = i + j*dim[0] + k*dim[0]*dim[1];
				float pnt[3];
				pnt[0] = orig[0] + i*span[0];
				pnt[1] = orig[1] + j*span[1];
				pnt[2] = orig[2] + k*span[2];
				data[l] = evalLaplace(pnt);
			}
	min = max = data[0];
	for (i = 1; i < nv; i++) {
		min = (data[i] < min)? data[i]:min;
		max = (data[i] > max)? data[i]:max;
	}
}

float Protein::evalDensity(const float pnt[3])
{
	int na = numOfAtoms();
	float val = 0;
	for (int i = 0; i < na; i++) {
		val += (float)atoms[i].elecDenGaussian(pnt);
	}
	return val;
}

void Protein::evalGradient(const float pnt[3], float grad[3])
{
	float atom_grad[3];
	grad[0] = grad[1] = grad[2] = 0;
	int na = numOfAtoms();
	for (int i = 0; i < na; i++) {
		atoms[i].elecGradGaussian(pnt, atom_grad);
		grad[0] += atom_grad[0];
		grad[1] += atom_grad[1];
		grad[2] += atom_grad[2];

		//double blobby = -1 * element[atom.id].ionenergy;
		/*
		double blobby = -4;
				double r =  (atom.center[0] - pnt[0]) * (atom.center[0] - pnt[0]) +
							(atom.center[1] - pnt[1]) * (atom.center[1] - pnt[1]) +
							(atom.center[2] - pnt[2]) * (atom.center[2] - pnt[2]);
				double r0 = element[atom.id].VdWradius * element[atom.id].VdWradius;
				grad[0] += (float)(2*blobby*(pnt[0]-atom.center[0])*exp(blobby*r/r0 - blobby)/r0);
				grad[1] += (float)(2*blobby*(pnt[1]-atom.center[1])*exp(blobby*r/r0 - blobby)/r0);
				grad[2] += (float)(2*blobby*(pnt[2]-atom.center[2])*exp(blobby*r/r0 - blobby)/r0);*/
		
	}
}

void Protein::evalDenAndGrad(const float pnt[3], float& den, float grad[3])
{
	den = 0;
	grad[0] = grad[1] = grad[2] = 0;
	Atom atom;
	//double blobby = -1 * element[atom.id].ionenergy;
	double blobby = -4;
	int na = numOfAtoms();
	for (int i = 0; i < na; i++) {
		atom = atoms[i];
		double r =  (atom.center[0] - pnt[0]) * (atom.center[0] - pnt[0]) +
			(atom.center[1] - pnt[1]) * (atom.center[1] - pnt[1]) +
			(atom.center[2] - pnt[2]) * (atom.center[2] - pnt[2]);
		double r0 = element[atom.id].VdWradius * element[atom.id].VdWradius;
		grad[0] += (float)(2*blobby*(pnt[0]-atom.center[0])*exp(blobby*r/r0 - blobby)/r0);
		grad[1] += (float)(2*blobby*(pnt[1]-atom.center[1])*exp(blobby*r/r0 - blobby)/r0);
		grad[2] += (float)(2*blobby*(pnt[2]-atom.center[2])*exp(blobby*r/r0 - blobby)/r0);
		den += (float)exp(blobby*r/r0 - blobby);
	}
}

float Protein::evalLaplace(const float pnt[3])
{
	float lap = 0;

	int na = numOfAtoms();
	for (int i = 0; i < na; i++) {
		Atom atom = atoms[i];
		//double blobby = -1 * element[atom.id].ionenergy;
		double blobby = -4;
		double r =  (atom.center[0] - pnt[0]) * (atom.center[0] - pnt[0]) +
					(atom.center[1] - pnt[1]) * (atom.center[1] - pnt[1]) +
					(atom.center[2] - pnt[2]) * (atom.center[2] - pnt[2]);
		double r0 = element[atom.id].VdWradius * element[atom.id].VdWradius;
		lap += (float)((6+4*blobby*r/r0) * blobby/r0 * exp(blobby*r/r0 - blobby));
	}
	return lap;
}

void Protein::getMolShell(int dim[3], float *data)
{
	float min[3], max[3];
	getBoundingBox(min, max);
	
	float orig[3], span[3];
	
	orig[0] = min[0]; orig[1] = min[1]; orig[2] = min[2];
	span[0] = (max[0]-min[0])/(dim[0]-1); 
	span[1] = (max[1]-min[1])/(dim[1]-1); 
	span[2] = (max[2]-min[2])/(dim[2]-1);

	int i, nv = dim[0]*dim[1]*dim[2];
	
	for (i = 0; i < nv; i++) data[i] = 0;
	for (i = 0; i < numOfAtoms(); i++) {
		atoms[i].accumShell(data, dim, orig, span);
	}
}



