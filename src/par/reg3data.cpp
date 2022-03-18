#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "reg3data.h"
#include "protein.h"
#include <bio.h>
#include "cubes.h"
#include "endians.h"

static int bitsize(int i)
{
    int b = 1, size=0;
	
    while (b <= i) {
        b<<=1;
        size++;
    }
    return(size);
}

Reg3Data::Reg3Data(const char* fname, int dim[3])
{
	m_dim[0] = m_dim[1] = m_dim[2] = 2;
	// m_orig is declared in header.
	m_orig[0] = m_orig[1] = m_orig[2] = 0;
	m_span[0] = m_span[1] = m_span[2] = 1;
	m_type = FLOAT_TYPE;
	p_grad = NULL;
	if(fname == NULL) {
		p_data.float_ptr = new float[8];
	} else {
		if(dim != NULL) {
			m_dim[0] = dim[0];
			m_dim[1] = dim[1];
			m_dim[2] = dim[2];
		}
		readFile(fname);
	}

	init();
}

Reg3Data::Reg3Data(int dim[3], void* array, DataType type)
{
	int i;

	for(i = 0; i < 3; i++) {
		m_dim[i] = dim[i];
		m_orig[i] = 0;
		m_span[i] = 1;
	}
	m_type = type;

	switch(m_type) {
	case UCHAR_TYPE:
		p_data.uchar_ptr = new unsigned char[m_dim[0]*m_dim[1]*m_dim[2]];
		memcpy(p_data.uchar_ptr, array, sizeof(char)*m_dim[0]*m_dim[1]*m_dim[2]);
		break;
	case USHORT_TYPE:
		p_data.ushort_ptr = new unsigned short[m_dim[0]*m_dim[1]*m_dim[2]];
		memcpy(p_data.ushort_ptr, array, sizeof(short)*m_dim[0]*m_dim[1]*m_dim[2]);
		break;
	case INT_TYPE:
		p_data.int_ptr = new int[m_dim[0]*m_dim[1]*m_dim[2]];
		memcpy(p_data.int_ptr, array, sizeof(int)*m_dim[0]*m_dim[1]*m_dim[2]);
		break;
	case FLOAT_TYPE:
		p_data.float_ptr = new float[m_dim[0]*m_dim[1]*m_dim[2]];
		memcpy(p_data.float_ptr, array, sizeof(float)*m_dim[0]*m_dim[1]*m_dim[2]);
		break;
	case DOUBLE_TYPE:
		p_data.double_ptr = new double[m_dim[0]*m_dim[1]*m_dim[2]];
		memcpy(p_data.double_ptr, array, sizeof(double)*m_dim[0]*m_dim[1]*m_dim[2]);
		break;
	default:
		error("Unknown Data Type");
	}
	p_grad = NULL;
	init();
}

Reg3Data::Reg3Data(int dim[3], Protein* _prtn)
{
	float min[3], max[3];

	m_type = FLOAT_TYPE;
	_prtn->getBoundingBox(min, max);
	for(int i = 0; i < 3; i++) {
		m_dim[i] = dim[i];
		m_orig[i] = min[i];
		m_span[i] = (max[i] - min[i]) / (m_dim[i]-1);
	}
	p_data.float_ptr = new float[m_dim[0]*m_dim[1]*m_dim[2]];
	//_prtn->getMolShell(m_dim, p_data.float_ptr);
	//p_grad = NULL;
	p_grad = new float[3*m_dim[0]*m_dim[1]*m_dim[2]];
	_prtn->getElecDenAndGradient(m_dim, p_data.float_ptr, p_grad);

	init();
}

Reg3Data::~Reg3Data()
{
	switch(m_type) {
	case UCHAR_TYPE:
		delete[] p_data.uchar_ptr;
		break;
	case USHORT_TYPE:
		delete[] p_data.ushort_ptr;
		break;
	case FLOAT_TYPE:
		delete[] p_data.float_ptr;
		break;
	case INT_TYPE:
		delete[] p_data.int_ptr;
		break;
	case DOUBLE_TYPE:
		delete[] p_data.double_ptr;
	}
	if(p_grad) delete[] p_grad;
}

float Reg3Data::getValue(float coord[3])
{
	int idx[3];
	float u[2], v[2], w[2];
	float data[8];
	//if(coord[0] < m_orig[0] || coord[0] >= (m_orig[0]+(dim[0]-1)*m_span[0]) return 0;
	//if(coord[1] < m_orig[1] || coord[1] >= (m_orig[1]+(dim[1]-1)*m_span[1]) return 0;
	//if(coord[2] < m_orig[2] || coord[2] >= (m_orig[2]+(dim[0]-1)*m_span[2]) return 0;
	for(int i = 0; i < 3; i++) {
		idx[i] = (coord[i] - m_orig[i])/m_span[i];
	}
	if(idx[0] < 0 || idx[0] >= (m_dim[0]-1)) return 0;
	if(idx[1] < 0 || idx[1] >= (m_dim[1]-1)) return 0;
	if(idx[2] < 0 || idx[2] >= (m_dim[2]-1)) return 0;
	
	getCellValues(idx[0], idx[1], idx[2], data);
	
	u[1] = (coord[0] - m_orig[0])/m_span[0] - idx[0]; u[0] = 1 - u[1];
	v[1] = (coord[1] - m_orig[1])/m_span[1] - idx[1]; v[0] = 1 - v[1];
	w[1] = (coord[2] - m_orig[2])/m_span[2] - idx[2]; w[0] = 1 - w[1];
	
	//printf("i, j, k = %d %d %d, u v w = %f %f %f\n", idx[0], idx[1], idx[2], u[1], v[1], w[1]);
	float val = 0;
	for(int k = 0; k < 2; k++) {
		for(int j = 0; j < 2; j++) {
			for(int i = 0; i < 2; i++) {
				val += 	data[vertinfo[k][j][i]]*u[i]*v[j]*w[k];
			}
		}
	}
	return val;
}

void Reg3Data::getCellValues(int i, int j, int k, float vals[8]) const
{
	int v = index2vert(i, j, k);
	vals[0] = Data::getValue(v);
	vals[1] = Data::getValue(v + 1);
	vals[2] = Data::getValue(v + m_xy + 1);
	vals[3] = Data::getValue(v + m_xy);
	vals[4] = Data::getValue(v + m_dim[0]);
	vals[5] = Data::getValue(v + m_dim[0] + 1);
	vals[6] = Data::getValue(v + m_xy + m_dim[0] + 1);
	vals[7] = Data::getValue(v + m_xy + m_dim[0]);
}

void Reg3Data::getCellGrads(int i, int j, int k, float grads[8][3]) 
{
    getVertGrad(i,     j,  k  , grads[0]);
    getVertGrad(i+1,   j,  k  , grads[1]);
    getVertGrad(i+1,   j,  k+1, grads[2]);
    getVertGrad(i,     j,  k+1, grads[3]);
    getVertGrad(i,   j+1,  k  , grads[4]);
    getVertGrad(i+1, j+1,  k  , grads[5]);
    getVertGrad(i+1, j+1,  k+1, grads[6]);
    getVertGrad(i,   j+1,  k+1, grads[7]);
}

void Reg3Data::calcGradient()
{
	if(hasGradient()) return;

	// Use Bspline approximation to calc gradients
    float v[27];
    int ix[3], iy[3], iz[3];
    int i, j, k, l, m, n, t;

	p_grad = new float[3*m_dim[0]*m_dim[1]*m_dim[2]];
	memset(p_grad, 0, sizeof(float)*3*m_dim[0]*m_dim[1]*m_dim[2]);
	for(k = 0; k < m_dim[2]; k++)
		for(j = 0; j < m_dim[1]; j++)
			for(i = 0; i < m_dim[0]; i++)
			{
				ix[0] = (i-1 >= 0)? i-1:0;
				ix[1] = i;
				ix[2] = (i+1 < m_dim[0])? i+1:i;
				iy[0] = (j-1 >= 0)? j-1:0;
				iy[1] = j;
				iy[2] = (j+1 < m_dim[1])? j+1:j;
				iz[0] = (k-1 >= 0)? k-1:0;
				iz[1] = k;
				iz[2] = (k+1 < m_dim[2])? k+1:k;

				t = 0;
				for (n = 0; n < 3; n++) {
					for (m = 0; m < 3; m++) {
						for (l = 0; l < 3; l++) {
							v[t] = getValue(ix[l], iy[m], iz[n]);
							t++;
						}
					}
				}

				float *ptr = p_grad + 3*index2vert(i, j, k); 
				for (l = 0; l < 27; l++) {
					*ptr += x_grad_mask[l]*v[l];
					*(ptr+1) += y_grad_mask[l]*v[l];
					*(ptr+2) += z_grad_mask[l]*v[l];
				}
				*ptr /= m_span[0];
				*(ptr+1) /= m_span[1];
				*(ptr+2) /= m_span[2];
			}
}

/************************************************************************/
/* Private Functions                                                    */
/************************************************************************/
void Reg3Data::init()
{
	m_xy = m_dim[1]*m_dim[0];
    m_xbits = bitsize(m_dim[0] - 2);
    m_ybits = bitsize(m_dim[1] - 2);
    m_zbits = bitsize(m_dim[2] - 2);
    if(m_xbits == 0)
        m_xbits = 1;
    if (m_ybits == 0)
        m_ybits = 1;
    if (m_zbits == 0)
        m_zbits = 1;
	
    m_yshift = m_xbits;
    m_zshift = m_xbits + m_ybits;
	
    m_xmask = (1<<m_xbits) - 1;
    m_ymask = (1<<m_ybits) - 1;
    m_zmask = (1<<m_zbits) - 1;

	calcGradient();

	int i, n = m_dim[0]*m_dim[1]*m_dim[2];
	switch(m_type) {
	case UCHAR_TYPE:
		m_fmin = m_fmax = p_data.uchar_ptr[0];
		for(i = 1; i < n; i++) {
			if(p_data.uchar_ptr[i] < m_fmin) m_fmin = p_data.uchar_ptr[i];
			else if(p_data.uchar_ptr[i] > m_fmax) m_fmax = p_data.uchar_ptr[i];
		}
		break;
	case USHORT_TYPE:
		m_fmin = m_fmax = p_data.ushort_ptr[0];
		for(i = 1; i < n; i++) {
			if(p_data.ushort_ptr[i] < m_fmin) m_fmin = p_data.ushort_ptr[i];
			else if(p_data.ushort_ptr[i] > m_fmax) m_fmax = p_data.ushort_ptr[i];
		}
		break;
	case FLOAT_TYPE:
		m_fmin = m_fmax = p_data.float_ptr[0];
		for(i = 1; i < n; i++) {
			if(p_data.float_ptr[i] < m_fmin) m_fmin = p_data.float_ptr[i];
			else if(p_data.float_ptr[i] > m_fmax) m_fmax = p_data.float_ptr[i];
		}
		break;
	case INT_TYPE:
		m_fmin = m_fmax = (float)p_data.int_ptr[0];
		for(i = 1; i < n; i++) {
			if(p_data.int_ptr[i] < m_fmin) m_fmin = (float)p_data.int_ptr[i];
			else if(p_data.int_ptr[i] > m_fmax) m_fmax = (float)p_data.int_ptr[i];
		}
		break;
	case DOUBLE_TYPE:
		m_fmin = m_fmax = (float)p_data.double_ptr[0];
		for(i = 1; i < n; i++) {
			if(p_data.double_ptr[i] < m_fmin) m_fmin = (float)p_data.double_ptr[i];
			else if(p_data.double_ptr[i] > m_fmax) m_fmax = (float)p_data.double_ptr[i];
		}
		break;
	}
}

void Reg3Data::readFile(const char* fname)
{
	// Does the system support multiple file formats?
	size_t len = strlen(fname);
	if (len > 4 && strcmp(fname+len-4, ".r3d") == 0) {
		readR3D(fname);
	} else if (len > 6 && strcmp(fname+len-6, ".rawiv") == 0) {
		readRawiv(fname);
	} else if (len > 4 && strcmp(fname+len-4, ".pdb") == 0) {
		readPDB(fname);
	} else {
		error("Unknown File Type");
	}
}

/**
 * Read regular 3D scalar function from .r3d file.
 * R3D file is a simple file format for regular 3D data. It supports multiple
 * time steps and multiple variables. The header of .r3d files has the following
 * structure:
 * { 
 *	int dim[3]; float orig[3]; float span[3];
 *	int	ntime, int nvar;
 *	iterate for each variable { char var_name[8]; int var_type;}
 * }
 * The Raw data follows the header in the order of variable-wise first and time-wise second.
 */
void Reg3Data::readR3D(const char* fname)
{
	int		i, nvar, ntime;
	char	(*var_names)[9];
	int		*var_types;
	
	DiskIO *pio = new BufferedIO(fname);
	if(!pio->open()) {
		error("Data File Open Failed");
	
	}
	pio->get(m_dim, 3);
	pio->get(m_orig, 3);
	pio->get(m_span, 3);
	pio->get(&ntime, 1);
	pio->get(&nvar, 1);

	var_names = (char (*)[9])malloc(sizeof(char[9])*nvar);
	memset(var_names, 0, 9*nvar);
	var_types = (int *)malloc(sizeof(int)*nvar);
	for(i = 0; i < nvar; i++) {
		pio->get(var_names[i], 8);
		pio->get(&(var_types[i]), 1);
	}
	
	int nverts = m_dim[0]*m_dim[1]*m_dim[2];
	switch(var_types[0]) {
	case 0:
		m_type = UCHAR_TYPE;
		p_data.uchar_ptr = new unsigned char[nverts];
		pio->get(p_data.uchar_ptr, nverts);
		break;
	case 1:
		m_type = USHORT_TYPE;
		p_data.ushort_ptr = new unsigned short[nverts];
		pio->get(p_data.ushort_ptr, nverts);
		break;
	case 2:
		m_type = FLOAT_TYPE;
		p_data.float_ptr = new float[nverts];
		pio->get(p_data.float_ptr, nverts);
		break;
	case 3:
		m_type = INT_TYPE;
		p_data.int_ptr = new int[nverts];
		pio->get(p_data.int_ptr, nverts);
		break;
	case 4:
		m_type = DOUBLE_TYPE;
		p_data.double_ptr = new double[nverts];
		pio->get(p_data.double_ptr, nverts);
		break;
	default:
		error("Data Type Error in R3D File");
	}
	// r3d file contains gradient
	if(nvar >= 4 && var_types[1] == 2 && var_types[2] == 2 && var_types[3] == 2 &&
	   (strcmp(var_names[1], "DX") == 0) && 
	   (strcmp(var_names[2], "DY") == 0) &&
	   (strcmp(var_names[3], "DZ") == 0) ) 
	{
		p_grad = new float[3*nverts];
		float *p0 = new float[nverts];
		float *p1 = new float[nverts];
		float *p2 = new float[nverts];
		pio->get(p0, nverts);
		pio->get(p1, nverts);
		pio->get(p2, nverts);
		for(i = 0; i < nverts; i++) {
			p_grad[3*i] = p0[i];
			p_grad[3*i+1] = p1[i];
			p_grad[3*i+2] = p2[i];
		}
		delete[] p0;
		delete[] p1;
		delete[] p2;
	}
	free(var_types);
	free(var_names);
	pio->close();
	delete pio;
}

/**
 * Read regular 3D function from the aged rawiv format.
 */
void Reg3Data::readRawiv(const char* fname)
{
	int nverts, ncells;
	float minext[3], maxext[3];
	
	// determine file data type
	struct stat filestat;
	if (stat(fname, &filestat) < 0) {
		fprintf( stderr, "cannot find data file: %s\n", fname );
	}
	int sz = filestat.st_size; 
	
	DiskIO *pio = new BufferedIO(fname);
	if(!pio->open()) {
		fprintf( stderr, "data file open fail: %s\n", fname );
	}
	pio->get(minext, 3);
	pio->get(maxext, 3);
	pio->get(&nverts, 1);
	pio->get(&ncells, 1);
	pio->get(m_dim, 3);
	pio->get(m_orig, 3);
	pio->get(m_span, 3);

	if(!big_endian())
        {
   	  for(int i=0; i<3; i++) SWAP_32(&(minext[i]));
	  for(int i=0; i<3; i++) SWAP_32(&(maxext[i]));
	  SWAP_32(&(nverts));
	  SWAP_32(&(ncells));
	  for(int i=0; i<3; i++) SWAP_32(&(m_dim[i]));
	  for(int i=0; i<3; i++) SWAP_32(&(m_orig[i]));
	  for(int i=0; i<3; i++) SWAP_32(&(m_span[i]));
        }

	nverts = m_dim[0]*m_dim[1]*m_dim[2]; 
	unsigned long int voxelTypeSize = (unsigned long int)( (unsigned long int)(sz-68)/((unsigned long int)(m_dim[0])*(unsigned long int)(m_dim[1])*(unsigned long int)(m_dim[2])));

	switch (voxelTypeSize) {
	case 1:
		m_type = UCHAR_TYPE;
		p_data.uchar_ptr = new unsigned char[nverts];
		pio->get(p_data.uchar_ptr, nverts);
		break;
	case 2:
		m_type = USHORT_TYPE;
		p_data.ushort_ptr = new unsigned short[nverts];
		pio->get(p_data.ushort_ptr, nverts);
		break;
	case 4:
		m_type = FLOAT_TYPE;
		p_data.float_ptr = new float[nverts];
		pio->get(p_data.float_ptr, nverts);
		break;
	default:
		error("Type Error in Rawiv File");
	}

	/* swap the volume data if on little endian machine */
	if(!big_endian())
	{
	   size_t len = m_dim[0] * m_dim[1] * m_dim[2];
	   switch(m_type)
	   {
	     case USHORT_TYPE: for(int i=0;i<len;i++) SWAP_16(p_data.ushort_ptr + i); break;
	     case FLOAT_TYPE:  for(int i=0;i<len;i++) SWAP_32(p_data.float_ptr + i); break;
	     default: break; /* no swapping needed for unsigned char data, and unsigned int is not defined for rawiv */
	   }
	}

#ifdef _DEBUG
	printf("ext: (%f, %f, %f) to (%f, %f, %f)\n", minext[0], minext[1], minext[2],
		maxext[0], maxext[1], maxext[2]);
	printf("dim: %d %d %d, data type: %d\n", m_dim[0], m_dim[1], m_dim[2], m_type);
	printf("orig: %f %f %f\n", m_orig[0], m_orig[1], m_orig[2]);
	printf("span: %f %f %f\n", m_span[0], m_span[1], m_span[2]);
#endif

	pio->close();
	delete pio;
}

/**
 * Construct regular 3D sampling from a PDB (protein data bank) file.
 */
void Reg3Data::readPDB(const char* fname)
{
	float blobby = -1;
	Protein mol(fname);
	float min[3], max[3];
	mol.getBoundingBox(min, max);
	for(int i = 0; i < 3; i++) {
		m_orig[i] = min[i];
		m_span[i] = (max[i] - min[i]) / (m_dim[i]-1);
	}
	p_data.float_ptr = new float[m_dim[0]*m_dim[1]*m_dim[2]];
	mol.setBlobby(blobby);
	//mol.getMolShell(m_dim, p_data.float_ptr);
	//p_grad = NULL;
	mol.getElectronDensity(m_dim, p_data.float_ptr);
	p_grad = new float[3*m_dim[0]*m_dim[1]*m_dim[2]];
	mol.getGradient(m_dim, p_grad);
}

void Reg3Data::writeRawiv(const char* fname)
{
	DiskIO *fio = new BufferedIO(fname, DiskIO::WRITE);
	if(!fio->open()) {
		error("cannot open the rawiv file");
	}

	float minext[3], maxext[3];
	int nverts, ncells;

	for(int i = 0; i < 3; i++) {
		minext[i] = m_orig[i];
		maxext[i] = m_orig[i] + (m_dim[i]-1)*m_span[i];
	}
	fio->put(minext, 3);
	fio->put(maxext, 3);
	nverts = m_dim[0]*m_dim[1]*m_dim[2];
	ncells = (m_dim[0]-1)*(m_dim[1]-1)*(m_dim[2]-1);
	fio->put(&nverts, 1);
	fio->put(&ncells, 1);
	fio->put(m_dim, 3);
	fio->put(m_orig, 3);
	fio->put(m_span, 3);

	switch(m_type) {
	case UCHAR_TYPE:
		fio->put(p_data.uchar_ptr, nverts);
		break;
	case USHORT_TYPE:
		fio->put(p_data.ushort_ptr, nverts);
		break;
	case FLOAT_TYPE:
		fio->put(p_data.float_ptr, nverts);
		break;
	default:
		break;
	}
	fio->close(false);
	delete fio;
}

