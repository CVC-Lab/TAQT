#ifndef PAR_DATA_H
#define PAR_DATA_H

#include <stdio.h>
#include <stdlib.h>

#define error(x) {fprintf(stderr, "%s\n", x); exit(1);}
typedef unsigned long size_n;

static const int vertinfo[2][2][2] = {
	{{0, 1}, {4, 5}},
	{{3, 2}, {7, 6}}
};

class Data {
public:
	/**
	* supported data types for scalar field.
	*/
	typedef enum {
		UCHAR_TYPE = 0,
		USHORT_TYPE,
		FLOAT_TYPE,
		INT_TYPE,
		DOUBLE_TYPE
	} DataType;

	/**
	* union of data pointers.
	*/
	union DataPtr {
		unsigned char*	uchar_ptr;
		unsigned short* ushort_ptr;
		float*			float_ptr;
		int*			int_ptr;
		double*			double_ptr;
	};

	float getValue(int nv) const {
		switch (m_type) {
		case UCHAR_TYPE:
			return p_data.uchar_ptr[nv];
		case USHORT_TYPE:
			return p_data.ushort_ptr[nv];
		case FLOAT_TYPE:
			return p_data.float_ptr[nv];
		case INT_TYPE:
			return (float)(p_data.int_ptr[nv]);
		case DOUBLE_TYPE:
			return (float)(p_data.double_ptr[nv]);
		}
		// should not come here
		error("Unknown Type");
	}

protected:
	DataType	m_type;
	DataPtr		p_data;
	float*		p_grad;

};

#endif