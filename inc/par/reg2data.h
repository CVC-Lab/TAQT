//--------------------------------------------------------------------
//
// Reg2Datad - class for a 2d regular grid of scalar data
//
// Copyright (c) 1997 Dan Schikore - updated by Emilio Camahort, 1999
//
//--------------------------------------------------------------------

// $Id: reg2data.h,v 1.2 2006/03/29 22:30:41 xyzhang Exp $

#ifndef Reg2DataD_H
#define Reg2DataD_H

#include <string.h>

#include "data.h"

// added by Emilio: these are forward declarations necessary to avoid
//		    complaints by the new SGI C++ compiler

typedef unsigned int u_int;
typedef unsigned char u_char;
static const int BAD_INDEX = -1;

//--------------------------------------------------------------------
//
// Reg2Data - a volume of scalar data.
//
//--------------------------------------------------------------------
class Reg2Data : public Data
{
protected:
	u_int m_dim[2];			// data members
	float m_orig[2];
	float m_span[2];

	int xbits, ybits;
    int xmask, ymask;
    int yshift;

public:				
	// constructors and destructors
	Reg2Data(const char *rawfile, DataType type = Data::FLOAT_TYPE);
	Reg2Data(int *dim, void *data, DataType type = Data::FLOAT_TYPE );
	~Reg2Data() {}

	// member access methods

	void	getDim(u_int *v)  const { memcpy(v, m_dim,  2 * sizeof(u_int)); }
	void	getOrig(float *v) const { memcpy(v, m_orig, 2 * sizeof(float)); }
	void	getSpan(float *v) const { memcpy(v, m_span, 2 * sizeof(float)); }

	void	setDim(int *v) {m_dim[0] = v[0]; m_dim[1] = v[1];}
	void	setOrig(float *v) {m_orig[0] = v[0]; m_orig[1] = v[1];}
	void	setSpan(float *v) {m_span[0] = v[0]; m_span[1] = v[1];}

	int		getNVerts() const { return m_dim[0]*m_dim[1];}

	void getFuncMinMax(float& min, float& max) const {
		min = max = getValue(0);
		for(int i = 1; i < m_dim[0]*m_dim[1]; i++) {
			if(max < getValue(i)) {max = getValue(i);}
			else if(min > getValue(i)) {min = getValue(i);}
		}
	}

	float xCoord(int i) { return(m_orig[0] + i*m_span[0]); }
	float yCoord(int j) { return(m_orig[1] + j*m_span[1]); }

	int maxCellIndex(void) { return(index2cell(m_dim[0]-2, m_dim[1]-2)); }

	float maxext[3], minext[3];
	
public :	// get data or gradient approximations (by differencing)

	void getCellValues(int c, float *val)
	{ 
		int i,j;
		cell2index(c,i,j);
		getCellValues(i,j,val);
	}

	void getCellValues(int i, int j, float *val)
	{
		val[0] = getValue(index2vert(i,j));
		val[1] = getValue(index2vert(i+1,j));
		val[2] = getValue(index2vert(i+1,j+1));
		val[3] = getValue(index2vert(i,j+1));
	}

	float getVertValue(int i, int j) {
		return getValue(index2vert(i, j));
	}

	u_int   getCellVert(int c, int v)
	{ 
		int i, j;
		cell2index(c,i,j);
		switch (v) {
		case 0:
			return(index2vert(i,j));
		case 1:
			return(index2vert(i+1,j));
		case 2:
			return(index2vert(i+1,j+1));
		case 3:
			return(index2vert(i,j+1));
		}
		error("BAD_INDEX");
		return -1;
	}

	u_int getNCellVerts(void) { return(4); }
	u_int getNCellFaces(void) { return(4); }
	int getCellAdj(int c, int f)
	{ 
		int i, j;
		cell2index(c,i,j);
		switch (f) {
		case 0:
			return(j==0? -1 : index2cell(i,j-1));
		case 1:
			return(i==(signed int)m_dim[0]-2? -1 : index2cell(i+1,j));
		case 2:
			return(j==(signed int)m_dim[1]-2? -1 : index2cell(i,j+1));
		case 3:
			return(i==0? -1 : index2cell(i-1,j));
		}
		return(-1);
	}

	void getCellRange(int c, float &min, float &max)
	{
		float t;
		u_int i;
		max = min = getValue(getCellVert(c,0));
		for (i=1; i<getNCellVerts(); i++)
			if ((t=getValue(getCellVert(c,i))) < min)
				min = t;
			else if (t > max)
				max = t;
	}

	void getFaceRange(u_int c, u_int f, float &min, float &max)
	{
		float t;
		min = max = getValue(getCellVert(c,f));
		if ((t=getValue(getCellVert(c,f<3?f+1:0))) < min)
			min = t;
		else if (t > max)
			max = t;
	}


	void cell2index(int c, int &i, int &j)
	{ 
		int _left;
		i = c&xmask;
		_left = c>>xbits;
		j = _left&ymask;
	}

	int index2cell(int i, int j)
	{ return((j << yshift) | i); }

	int index2vert(int i, int j)
	{ return(i + j*m_dim[0]); }

	void vert2index(int v, int &i, int &j) const {
		i = v % m_dim[0];
        j = v / m_dim[0];
	}

	bool writeRawiv(const char* fname);

protected:
	void init();

	// read data from a file
	void readFile(const char* fname);

	void readRawiv(const char* fname);
};

#endif
