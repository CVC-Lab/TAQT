#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>

//#include "mesh.h"
//#include "free.h"
#include "function.h"
#include "protein.h"
#include "pcio.h"

Function3D::Function3D(const char* fname)
: p_data(NULL)
{
	dim[0] = dim[1] = dim[2] = 0;

	// vtm file
	size_t len = strlen(fname);
	if (len > 4 && strcmp(fname+len-4, ".vtm") == 0) {
		readVTM(fname);
	} else if (len > 6 && strcmp(fname+len-6, ".rawiv") == 0) {
		readRawiv(fname);
	} else if (len > 4 && strcmp(fname+len-4, ".pdb") == 0) {
		readPDB(fname);
	}
	else {
		fprintf(stderr, "Unknown file type\n");
	}
}

Function3D::Function3D(int _dim[3], float* vals) 
{
	for (int i = 0; i < 3; i++)	dim[i] = _dim[i];
	p_data = new float[dim[0]*dim[1]*dim[2]];
	memcpy(p_data, vals, sizeof(float)*dim[0]*dim[1]*dim[2]);
}

Function3D::~Function3D(void)
{
	if (p_data != NULL)	delete[] p_data;
}

void Function3D::readVTM(const char* fname)
{
	// assume single timestep and single variable
	char* file = strdup(fname);
	VTmesh* p_mesh = ::VTreadMesh(file);
	free(file);
	if (p_mesh == NULL)	return;
	// check if data is rectilinear 3D
	if (VTmeshMeshType(p_mesh) == VT_MESH_RECT3D) {
		VTrect3Dmesh* p_rect3d = p_mesh->mesh.rect3d;
		dim[0] = p_rect3d->xdim;
		dim[1] = p_rect3d->ydim;
		dim[2] = p_rect3d->zdim;
		int i, nv = dim[0]*dim[1]*dim[2];
		p_data = new float[nv];
		switch (VTmeshVarType(p_mesh, 0)) {
		case VT_UCHAR:
			for (i = 0; i < nv; i++) {
				p_data[i] = p_rect3d->data[0]->data[0]->data.ucdata[i];
			}
			break;
		case VT_CHAR:
			for (i = 0; i < nv; i++) {
				p_data[i] = p_rect3d->data[0]->data[0]->data.cdata[i];
			}
			break;
		case VT_SHORT:
			for (i = 0; i < nv; i++) {
				p_data[i] = p_rect3d->data[0]->data[0]->data.sdata[i];
			}
			break;
		case VT_LONG:
			for (i = 0; i < nv; i++) {
				p_data[i] = (float)p_rect3d->data[0]->data[0]->data.ldata[i];
			}
			break;
		case VT_FLOAT:
			for (i = 0; i < nv; i++) {
				p_data[i] = p_rect3d->data[0]->data[0]->data.fdata[i];
			}
			break;
		}
	}
	VTfreeMesh(p_mesh);
}

static float orig[3], span[3];

void Function3D::readRawiv(const char* fname)
{
	int i, nverts, ncells;
	float minext[3], maxext[3];
	// float orig[3], span[3];

	// determine file data type
	struct stat filestat;
	if (stat(fname, &filestat) < 0) {
		fprintf(stderr, "cannot find data file %s\n", fname);
		return;
	}
	int sz = filestat.st_size;  

	FILE *fp = fopen(fname, "rb");
	if (fp == NULL) {
		fprintf(stderr, "cannot open %s\n", fname);
		return;
	}
#ifdef _LITTLE_ENDIAN
	getFloat(minext, 3, fp);
	getFloat(maxext, 3, fp);
	getInt(&nverts, 1, fp);
	getInt(&ncells, 1, fp);
	getInt(dim, 3, fp);
	getFloat(orig, 3, fp);
	getFloat(span, 3, fp);
#else
	fread(minext, sizeof(float[3]), 1, fp);
	fread(maxext, sizeof(float[3]), 1, fp); 
	fread(&nverts, sizeof(int), 1, fp);
	fread(&ncells, sizeof(int), 1, fp);
	fread(dim, sizeof(int), 3, fp);
	fread(orig, sizeof(float), 3, fp);
	fread(span, sizeof(float), 3, fp);
#endif  
	nverts = dim[0]*dim[1]*dim[2]; 
	int dtype = sz / (2 * nverts);
	printf("dim: %d %d %d, data type: %d\n", dim[0], dim[1], dim[2], dtype);
	printf("orig: %f %f %f\n", orig[0], orig[1], orig[2]);
	printf("span: %f %f %f\n", span[0], span[1], span[2]);

	p_data = new float[nverts];
	switch (dtype) {
	case 0:	{	// unsigned char 
		unsigned char* ucdata = new unsigned char[nverts];
		fread(ucdata, sizeof(unsigned char), nverts, fp);
		for(i = 0; i < nverts; i++) {
			p_data[i] = ucdata[i];
		}
		delete[] ucdata;
		break;
	}
	case 1:	{	// unsigned short
		unsigned short* usdata = new unsigned short[nverts];
#ifdef _LITTLE_ENDIAN
		getShort((short *)usdata, nverts, fp);
#else
		fread(usdata, sizeof(short), nverts, fp);
#endif
		for(i = 0; i < nverts; i++) {
			p_data[i] = usdata[i];
		}
		delete[] usdata;
		break;
	}
	case 2:		// float
#ifdef _LITTLE_ENDIAN
		getFloat(p_data, nverts, fp);
#else
		fread(p_data, sizeof(float), nverts, fp);
#endif
		break;
	}
}

void Function3D::readPDB(const char* fname)
{
	Protein* mol = ProteinFactory::makeProtein(fname);
	assert(mol);

	dim[0] = 64;
	dim[1] = 64;
	dim[2] = 64;

	printf("protein has %d atoms\n", mol->getNumOfAtoms());
	
	p_data = mol->getElectronDensity(dim);
	delete mol;
}
	
void Function3D::calcGradient(const char* fname)
{
	float *dx, *dy, *dz;

	dx = new float[dim[0]*dim[1]*dim[2]];
	dy = new float[dim[0]*dim[1]*dim[2]];
	dz = new float[dim[0]*dim[1]*dim[2]];

	finiteDiff(dx, dy, dz);	

	FILE* fp = fopen(fname, "wb");
	if(fp != NULL) {
		// need fix for little endian
		fwrite(dx, sizeof(float), dim[0]*dim[1]*dim[2], fp);
		fwrite(dy, sizeof(float), dim[0]*dim[1]*dim[2], fp);
		fwrite(dz, sizeof(float), dim[0]*dim[1]*dim[2], fp);
		fclose(fp);
	}

	delete[] dx;
	delete[] dy;
	delete[] dz;
}

void Function3D::calcGradientLength(const char* fname)
{
	float *dx, *dy, *dz, *len;

	dx = new float[dim[0]*dim[1]*dim[2]];
	dy = new float[dim[0]*dim[1]*dim[2]];
	dz = new float[dim[0]*dim[1]*dim[2]];
	len = new float[dim[0]*dim[1]*dim[2]];

	finiteDiff(dx, dy, dz);	

	for(int i = 0; i < dim[0]*dim[1]*dim[2]; i++) {
		len[i] = (float) sqrt(dx[i]*dx[i] + dy[i]*dy[i] + dz[i]*dz[i]);
	}
	FILE* fp = fopen(fname, "wb");
	if(fp != NULL) {
		// need fix for little endian
		fwrite(len, sizeof(float), dim[0]*dim[1]*dim[2], fp);
		fclose(fp);
	}

	delete[] len;
	delete[] dx;
	delete[] dy;
	delete[] dz;
}

void Function3D::calcLaplacian(const char* fname)
{
	int i, j, k;
	float *dx, *dy, *dz, *laplace, sum;

	dx = new float[dim[0]*dim[1]*dim[2]];
	dy = new float[dim[0]*dim[1]*dim[2]];
	dz = new float[dim[0]*dim[1]*dim[2]];
	laplace = new float[dim[0]*dim[1]*dim[2]];

	finiteDiff(dx, dy, dz);
	int n = 0;
	for(k = 0; k < dim[2]; k++) {
		for(j = 0; j < dim[1]; j++) {
			for(i = 0; i < dim[0]; i++) {
				sum = 0;
				if(i == 0) {
					sum += (dx[index2vert(i+1, j, k)] - dx[index2vert(i, j, k)]) / span[0];
				} else if(i == (int)dim[0]-1) {
					sum += (dx[index2vert(i, j, k)] - dx[index2vert(i-1, j, k)]) / span[0];
				} else {
					sum += 0.5f * (dx[index2vert(i+1, j, k)] - dx[index2vert(i-1, j, k)]) / span[0];
				}

				if(j == 0) {
					sum += (dy[index2vert(i, j+1, k)] - dy[index2vert(i, j, k)]) / span[1];
				} else if(j == (int)dim[1]-1) {
					sum += (dy[index2vert(i, j, k)] - dy[index2vert(i, j-1, k)]) / span[1];
				} else {
					sum += 0.5f * (dy[index2vert(i, j+1, k)] - dy[index2vert(i, j-1, k)]) / span[1];
				}

				if(k == 0) {
					sum += (dz[index2vert(i, j, k+1)] - dz[index2vert(i, j, k)]) / span[2];
				} else if(k == (int)dim[2]-1) {
					sum += (dz[index2vert(i, j, k)] - dz[index2vert(i, j, k-1)]) / span[2];
				} else {
					sum += 0.5f * (dz[index2vert(i, j, k+1)] - dz[index2vert(i, j, k-1)]) / span[2];
				}

				laplace[n] = sum;
				n++;
			}
		}
	}

	FILE* fp = fopen(fname, "wb");
	if(fp != NULL) {
		fwrite(laplace, sizeof(float), dim[0]*dim[1]*dim[2], fp);
		fclose(fp);
	}

	delete[] dx;
	delete[] dy;
	delete[] dz;
	delete[] laplace;
}

void Function3D::finiteDiff(float* dx, float* dy, float* dz)
{
	int i, j, k;

	int n = 0;
	for(k = 0; k < dim[2]; k++) {
		for(j = 0; j < dim[1]; j++) {
			for(i = 0; i < dim[0]; i++) {
				if (i==0) {
					/* use right difference */
					dx[n] = getValue(i+1, j, k) - getValue(i, j, k) ;
				} else if (i == int(dim[0]-1)) {
					/* use left difference */
					dx[n] = getValue(i, j, k) - getValue(i-1, j, k);
				} else {
					/* use central difference */
					dx[n] = (getValue(i+1, j, k) - getValue(i-1, j, k)) * 0.5f;
				}

				if (j==0) {
					dy[n] = getValue(i, j+1, k) - getValue(i, j, k);
				} else if (j == int(dim[1]-1)) {
					dy[n] = getValue(i, j, k) - getValue(i, j-1, k);
				} else {
					dy[n] = (getValue(i, j+1, k) - getValue(i, j-1, k)) * 0.5f;
				}

				if (k==0) {
					dz[n] = getValue(i, j, k+1) - getValue(i, j, k);
				} else if (k == int(dim[2]-1)) {
					dz[n] = getValue(i, j, k) - getValue(i, j, k-1);
				} else {
					dz[n] = (getValue(i, j, k+1) - getValue(i, j, k-1)) * 0.5f;
				}

				dx[n] = dx[n] / span[0];
				dy[n] = dy[n] / span[1];
				dz[n] = dz[n] / span[2];
			
				n++;
			}
		}
	}
}


