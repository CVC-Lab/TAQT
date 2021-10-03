#include <cmath>

#include <par/atom.h>
#include <par/geom.h>
#include <par/elements.h>

float Atom::a = 3;
float Atom::b = 2;

Atom::Atom(int _id, float x, float y, float z)
{
	id = _id;
	center[0] = x;
	center[1] = y;
	center[2] = z;
	init();
}

Atom::Atom(int _id, float _center[3])
{
	id = _id;
	for(int i = 0; i < 3; i++) {
		center[i] = _center[i];
	}
	init();
}

void Atom::init()
{
	blobby = -4;
	thickness = 0.7;
}

double Atom::elecDenGaussian(const float pnt[3])
{
	double r2 = (center[0] - pnt[0]) * (center[0] - pnt[0]) +
				(center[1] - pnt[1]) * (center[1] - pnt[1]) +
				(center[2] - pnt[2]) * (center[2] - pnt[2]);
	double r02 = element[id].VdWradius * element[id].VdWradius;
	return exp(blobby*r2/r02 - blobby);
}

void Atom::elecGradGaussian(const float pnt[3], float grad[3])
{
	//double blobby = -1 * element[atom.id].ionenergy;
	double r2 = (center[0] - pnt[0]) * (center[0] - pnt[0]) +
				(center[1] - pnt[1]) * (center[1] - pnt[1]) +
				(center[2] - pnt[2]) * (center[2] - pnt[2]);
	double r02 = element[id].VdWradius * element[id].VdWradius;
	grad[0] = (float)(2*blobby*(pnt[0]-center[0])*exp(blobby*r2/r02 - blobby)/r02);
	grad[1] = (float)(2*blobby*(pnt[1]-center[1])*exp(blobby*r2/r02 - blobby)/r02);
	grad[2] = (float)(2*blobby*(pnt[2]-center[2])*exp(blobby*r2/r02 - blobby)/r02);
}

void Atom::elecDenGradGaussian(const float pnt[3], float& den, float grad[3])
{
	double r2 = (center[0] - pnt[0]) * (center[0] - pnt[0]) +
		(center[1] - pnt[1]) * (center[1] - pnt[1]) +
		(center[2] - pnt[2]) * (center[2] - pnt[2]);
	double r02 = element[id].VdWradius * element[id].VdWradius;
	den = (float)exp(blobby*r2/r02 - blobby);
	grad[0] = (float)(2*blobby*(pnt[0]-center[0])*exp(blobby*r2/r02 - blobby)/r02);
	grad[1] = (float)(2*blobby*(pnt[1]-center[1])*exp(blobby*r2/r02 - blobby)/r02);
	grad[2] = (float)(2*blobby*(pnt[2]-center[2])*exp(blobby*r2/r02 - blobby)/r02);
}

double Atom::elecDenMeta(const float pnt[3])
{
	double r =  sqrt((center[0] - pnt[0]) * (center[0] - pnt[0]) +
				(center[1] - pnt[1]) * (center[1] - pnt[1]) +
				(center[2] - pnt[2]) * (center[2] - pnt[2]));
	r = r / (b * element[id].VdWradius);

	if(r <= 0.333333333333) {
		return a * (1 - r*r);
	} else if (r <= 1) {
		return a * (1-r) * (1-r);
	}
	return 0;
}

void Atom::elecGradMeta(const float pnt[3], float grad[3])
{
	/*
	 *	Todo
	 */
}
void Atom::elecDenGradMeta(const float pnt[3], float& den, float grad[3])
{
	/*
	 *	Todo
	 */
}

void Atom::accumDensity(float *data, int dim[3], float orig[3], float span[3])
{
	// The influence radius of the atom is b*VdWradius
	int min_idx[3], max_idx[3];
	int i, j, k, n;
	float pnt[3];

	for(i = 0; i < 3; i++) {
		min_idx[i] = (int)MAX(0, (center[i]-b*element[id].VdWradius-orig[i])/span[i]);
		max_idx[i] = (int)MIN(dim[i]-1, (center[i]+b*element[id].VdWradius-orig[i])/span[i]+0.5);
	}

	for(k = min_idx[2]; k <= max_idx[2]; k++) {
		for(j = min_idx[1]; j <= max_idx[1]; j++) {
			for(i = min_idx[0]; i <= max_idx[0]; i++) {
				pnt[0] = orig[0] + i*span[0];
				pnt[1] = orig[1] + j*span[1];
				pnt[2] = orig[2] + k*span[2];
				n = i + j*dim[0] + k*dim[0]*dim[1];
				data[n] += (float) elecDenGaussian(pnt);
			}
		}
	}
}

void Atom::accumGradient(float *grads, int dim[3], float orig[3], float span[3])
{
	int min_idx[3], max_idx[3];
	int i, j, k, n;
	float pnt[3], grad[3];
	
	for(i = 0; i < 3; i++) {
		min_idx[i] = (int)MAX(0, (center[i]-b*element[id].VdWradius-orig[i])/span[i]);
		max_idx[i] = (int)MIN(dim[i]-1, (center[i]+b*element[id].VdWradius-orig[i])/span[i]+0.5);
	}
	
	for(k = min_idx[2]; k <= max_idx[2]; k++) {
		for(j = min_idx[1]; j <= max_idx[1]; j++) {
			for(i = min_idx[0]; i <= max_idx[0]; i++) {
				pnt[0] = orig[0] + i*span[0];
				pnt[1] = orig[1] + j*span[1];
				pnt[2] = orig[2] + k*span[2];
				n = i + j*dim[0] + k*dim[0]*dim[1];
				elecGradGaussian(pnt, grad);
				grads[3*n] += grad[0];
				grads[3*n+1] += grad[1];
				grads[3*n+2] += grad[2];
			}
		}
	}
}

void Atom::accumDenAndGrad(float *data, float *grads, int dim[3], float orig[3], float span[3])
{
	int min_idx[3], max_idx[3];
	int i, j, k, n;
	float pnt[3], grad[3], den;
	
	for(i = 0; i < 3; i++) {
		min_idx[i] = (int)MAX(0, (center[i]-b*element[id].VdWradius-orig[i])/span[i]);
		max_idx[i] = (int)MIN(dim[i]-1, (center[i]+b*element[id].VdWradius-orig[i])/span[i]+0.5);
	}
	
	for(k = min_idx[2]; k <= max_idx[2]; k++) {
		for(j = min_idx[1]; j <= max_idx[1]; j++) {
			for(i = min_idx[0]; i <= max_idx[0]; i++) {
				pnt[0] = orig[0] + i*span[0];
				pnt[1] = orig[1] + j*span[1];
				pnt[2] = orig[2] + k*span[2];
				n = i + j*dim[0] + k*dim[0]*dim[1];
				elecDenGradGaussian(pnt, den, grad);
				data[n] += den;
				grads[3*n] += grad[0];
				grads[3*n+1] += grad[1];
				grads[3*n+2] += grad[2];
			}
		}
	}
}

float Atom::evalShell(const float pnt[3])
{
	double r = sqrt((center[0] - pnt[0]) * (center[0] - pnt[0]) +
					(center[1] - pnt[1]) * (center[1] - pnt[1]) +
					(center[2] - pnt[2]) * (center[2] - pnt[2]));

	if(r > element[id].VdWradius+thickness) return 0;
	else if( r >= element[id].VdWradius) return 1;
	return -10;
}

void Atom::accumShell(float *data, int dim[3], float orig[3], float span[3])
{
	int i, j, k, n, min_idx[3], max_idx[3];
	float pnt[3];

	for(i = 0; i < 3; i++) {
		min_idx[i] = (int)MAX(0, (center[i]-2*element[id].VdWradius-orig[i])/span[i]);
		max_idx[i] = (int)MIN(dim[i]-1, (center[i]+2*element[id].VdWradius-orig[i])/span[i]+0.5);
	}
	for(k = min_idx[2]; k <= max_idx[2]; k++) {
		for(j = min_idx[1]; j <= max_idx[1]; j++) {
			for(i = min_idx[0]; i <= max_idx[0]; i++) {
				pnt[0] = orig[0] + i*span[0];
				pnt[1] = orig[1] + j*span[1];
				pnt[2] = orig[2] + k*span[2];
				n = i + j*dim[0] + k*dim[0]*dim[1];
				data[n] += evalShell(pnt);
			}
		}
	}
}

