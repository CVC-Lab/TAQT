#ifndef _REG3_DATA_H
#define _REG3_DATA_H

#include <stdio.h>
#include <stdlib.h>

#include "data.h"

class Protein;
/**
 * Reg3Data: A class representing a scalar field on a regular 3D grid.
 *
 * @note This class is based on many previous versions of regular 3D data class.
 */

class Reg3Data : public Data

{
public:
	// Contruct the scalar function from a file
	Reg3Data(const char* fname = NULL, int dim[3] = NULL);

	// Construct from a data array
	Reg3Data(int dim[3], void* array, DataType type = Data::FLOAT_TYPE);

	// Construct from a protein molecule
	Reg3Data(int dim[3], Protein* _prtn);

	// Destructor
	~Reg3Data();

	void getDim(int dim[3]) const {
		for(int i = 0; i < 3; i++) dim[i] = m_dim[i];
	}

	void getOrig(float orig[3]) const {
		for(int i = 0; i < 3; i++) orig[i] = m_orig[i];
	}

	void setOrig(const float orig[3]) {
		for(int i = 0; i < 3; i++) m_orig[i] = orig[i];
	}

	void getSpan(float span[3]) const {
		for(int i = 0; i < 3; i++) span[i] = m_span[i];
	}

	void setSpan(const float span[3]) {
		for(int i = 0; i < 3; i++) m_span[i] = span[i];
	}

	// get the number of vertices 
	int getNVerts() const {
		return m_dim[0]*m_dim[1]*m_dim[2];
	}

	int getNCells() const {
		return (m_dim[0]-1)*(m_dim[1]-1)*(m_dim[2]-1);
	}
	float getFuncMin() const {
		return m_fmin;
	}

	float getFuncMax() const {
		return m_fmax;
	}

	void getFuncMinMax(float& min, float& max) const {
		min = m_fmin; max = m_fmax;
	}

	void getBoundingBox(float min[3], float max[3]) const {
		for(int i = 0; i < 3; i++) {
			min[i] = m_orig[i];
			max[i] = m_orig[i] + (m_dim[i]-1)*m_span[i];
		}
	}
	
	int getMaxCellIndex() const {
		return index2cell(m_dim[0]-1, m_dim[1]-1, m_dim[2]-1);
	}

	int getVertCount() const {
		return m_dim[0] * m_dim[1] * m_dim[2];
	}
	
	void getCoord(int id, float coord[3]) const {
		int i, j, k;
		vert2index(id, i, j, k);
		coord[0] = m_orig[0] + m_span[0]*i;
		coord[1] = m_orig[1] + m_span[1]*j;
		coord[2] = m_orig[2] + m_span[2]*k;
	}

	void setValue(int nv, float x) {
		switch (m_type) {
		case UCHAR_TYPE:
			p_data.uchar_ptr[nv] = (unsigned char)x;
			break;
		case USHORT_TYPE:
			p_data.ushort_ptr[nv] = (unsigned short)x;
			break;
		case FLOAT_TYPE:
			p_data.float_ptr[nv] = x;
			break;
		case INT_TYPE:
			p_data.int_ptr[nv] = (int)x;
			break;
		case DOUBLE_TYPE:
			p_data.double_ptr[nv] = x;
			break;
		}
	}

	float operator [] (int nv) const {
		return Data::getValue(nv);
	}

	float getValue(int nv) const {
		return Data::getValue(nv);
	}

	float getValue(int i, int j, int k) const {
		return Data::getValue(index2vert(i, j, k));
	}
	
	/**
	 * Get the function value at a given coordinates
	 * @note return 0 if the coordinates are outside of the bounding box
	 */
	float getValue(float coord[3]);
	
	inline void getVertGrad(int i, int j, int k, float grad[3]); 

	void getCellValues(int i, int j, int k, float vals[8]) const;

	void getCellValues(int c, float vals[8]) {
		int i, j, k;
		cell2index(c, i, j, k);
		getCellValues(i, j, k, vals);
	}

	void getCellGrads(int i, int j, int k, float grads[8][3]);

	void getCellGrads(int c, float grads[8][3]) {
		int i, j, k;
		cell2index(c, i, j, k);
		getCellGrads(i, j, k, grads);
	}

	/**
	 * Compute the gradient vector at all grid points.
	 */
	void calcGradient();

	/**
	 * Save the regular 3D data to a rawiv file
	 * @param dim Dimension of the volume
	 * @param write the ith variable
	 */
	void writeRawiv(const char* fname);

protected:
	// min max values of the scalar function
	float		m_fmin, m_fmax;	
	int			m_dim[3];
	float		m_orig[3], m_span[3];

private:
	// m_xy = dim[1]*dim[0];
	int m_xy;                       
    // Those variables are to facilitate index2cell and cell2index funcs
    int m_xbits, m_ybits, m_zbits;
    int m_xmask, m_ymask, m_zmask;
    int m_yshift, m_zshift;

	void init();

	void readFile(const char* fname);

	void readR3D(const char* fname);

	void readRawiv(const char* fname);

	void readPDB(const char* fname);
	
	bool hasGradient() const {
		return (p_grad != NULL);
	}

	int index2vert(int i, int j, int k) const {
		return (k*m_xy+j*m_dim[0]+i);
	}

	void vert2index(int v, int &i, int &j, int &k) const {
		i = v % m_dim[0];
        j = (v / m_dim[0]) % m_dim[1];
        k = v / (m_dim[0]*m_dim[1]);
	}

	int index2cell(int i, int j, int k) const {
		return((k << m_zshift) | (j << m_yshift) | i);		
	}

	void cell2index(int c, int &i, int &j, int &k) const {
	    int _left;
        i = c&m_xmask;
        _left = c>>m_xbits;
        j = _left&m_ymask;
        _left = _left>>m_ybits;
        k = _left&m_zmask;
    }
};

void Reg3Data::getVertGrad(int i, int j, int k, float grad[3])
{
	int v = index2vert(i, j, k); 
	if(hasGradient()) {
		// read from saved gradients
		grad[0] = p_grad[3*v];
		grad[1] = p_grad[3*v+1];
		grad[2] = p_grad[3*v+2];
	} else {
		// central differences
		if (i==0) {
			// use right difference
			grad[0] = Data::getValue(v+1) - Data::getValue(v);
        } else if (i == m_dim[0]-1) {
			// use left difference
			grad[0] = Data::getValue(v) - Data::getValue(v-1);
        } else {
			// use central difference
			grad[0] = (Data::getValue(v+1) - Data::getValue(v-1)) / 2;
        }

		if (j == 0) {
			grad[1] = Data::getValue(v+m_dim[0]) - Data::getValue(v);
		} else if (j  ==  m_dim[1]-1) {
			grad[1] = Data::getValue(v) - Data::getValue(v-m_dim[0]);
		} else {
			grad[1] = (Data::getValue(v+m_dim[0]) - Data::getValue(v-m_dim[0])) / 2;
        }

		if (k == 0) {
			grad[2] = Data::getValue(v+m_xy) - Data::getValue(v);
        } else if (k == m_dim[2]-1) {
			grad[2] = Data::getValue(v) - Data::getValue(v-m_xy);
		} else {
			grad[2] = (Data::getValue(v+m_xy) - Data::getValue(v-m_xy)) / 2;
		}
		grad[0] /= m_span[0];
		grad[1] /= m_span[1];
		grad[2] /= m_span[2];
	}
}

#endif


