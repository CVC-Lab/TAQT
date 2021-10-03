//------------------------------------------------------------------------
//
// Reg2Data.C - class for a regular 2d grid of scalar data
//
// Copyright (c) 1997 Dan Schikore - modified by Emilio Camahort, 1999
//
//------------------------------------------------------------------------

// $Id: reg2data.cpp,v 1.2 2006/03/29 22:30:41 xyzhang Exp $

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef MACOS_X
#include <malloc.h>
#endif
#include <string.h>

#ifndef WIN32
#include <unistd.h>
#endif

#include <bio.h>
#include "reg2data.h"

#define SQR(x) ((x)*(x))

#define FSAMPLES 256

#define VERBOSE 1

//------------------------------------------------------------------------
//
// bitsize - return the number of bits in an int (easier way?)
//
//------------------------------------------------------------------------
static int bitsize(unsigned int i)
{
   u_int b = 1, size=0;

   while (b <= i) {
      b<<=1;
      size++;
   }
   return(size);
}

//------------------------------------------------------------------------
//
// Reg2Data() - constructor to initialize the data
//
//------------------------------------------------------------------------

void Reg2Data::init()
{
	xbits = bitsize(m_dim[0] - 2);
	ybits = bitsize(m_dim[1] - 2);
	if (xbits == 0)
		xbits = 1;
	if (ybits == 0)
		ybits = 1;

	yshift = xbits;

	xmask = (1<<xbits) - 1;
	ymask = (1<<ybits) - 1;
}

Reg2Data::Reg2Data(const char *fn, DataType t)
{
	m_type = t;
#ifdef VERBOSE
	printf("reading dimensions\n");
#endif
	readFile(fn);
	p_grad = NULL;
	init();
}

void Reg2Data::readFile(const char* fname)
{
	size_t len = strlen(fname);
	if (len > 7 && strcmp(fname+len-7, ".rawiv2") == 0) {
		readRawiv(fname);
	} else {
		error("Unknown File Type");
	}
}

void Reg2Data::readRawiv(const char* fname)
{
	int nverts, ncells;

	DiskIO *pio = new BufferedIO(fname);
	if(!pio->open()) {
		error("Data File Open Failed");
	
	}
	pio->get(minext, 3);
	pio->get(maxext, 3);
	pio->get(&nverts, 1);
	pio->get(&ncells, 1);
	pio->get(m_dim, 2);
	pio->get(m_orig, 2);
	pio->get(m_span, 2);

#ifdef _DEBUG
	printf("dim: %d %d\n", m_dim[0], m_dim[1]);
	printf("orig: %f %f\n", m_orig[0], m_orig[1]);
	printf("span: %f %f\n", m_span[0], m_span[1]);
#endif
	
	nverts = m_dim[0]*m_dim[1];
	switch (m_type) {
	case UCHAR_TYPE:
		p_data.uchar_ptr = new unsigned char[nverts];
		pio->get(p_data.uchar_ptr, nverts);
		break;
	case USHORT_TYPE:
		p_data.ushort_ptr = new unsigned short[nverts];
		pio->get(p_data.ushort_ptr, nverts);
		break;
	case FLOAT_TYPE:
		p_data.float_ptr = new float[nverts];
		pio->get(p_data.float_ptr, nverts);
		break;
	default:
		error("Type Error in Rawiv File");
	}
	pio->close();
	delete pio;
}

bool Reg2Data::writeRawiv(const char* fname)
{
	DiskIO *fio = new BufferedIO(fname, DiskIO::WRITE);
	if(!fio->open()) {
		error("cannot open the rawiv file");
		return false;
	}

	int nverts, ncells;

	for(int i = 0; i < 2; i++) {
		minext[i] = m_orig[i];
		maxext[i] = m_orig[i] + (m_dim[i]-1)*m_span[i];
	}

	fio->put(minext, 3);
	fio->put(maxext, 3);
	nverts = m_dim[0]*m_dim[1];
	ncells = (m_dim[0]-1)*(m_dim[1]-1);
	fio->put(&nverts, 1);
	fio->put(&ncells, 1);
	fio->put(m_dim, 2);
	fio->put(m_orig, 2);
	fio->put(m_span, 2);

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
	return true;
}

//------------------------------------------------------------------------
//
// Reg2Data() - alternative constructor for the libcontour library
//
//------------------------------------------------------------------------

Reg2Data::Reg2Data(int *dim, void *array, Data::DataType t)
{
    int nverts = dim[0] * dim[1];
    int ncells = (dim[0]-1) * (dim[1]-1);

    //fread(&nverts, sizeof(int), 1, fp);
    //fread(&ncells, sizeof(int), 1, fp);

#ifdef _DEBUG
    printf("%d verts, %d cells\n", nverts, ncells);
	printf("reading dimensions\n");
#endif

   m_dim[0] = dim[0];
   m_dim[1] = dim[1];

   m_orig[0] = m_orig[1] = 0.0f;
   m_span[0] = m_span[1] = 1.0f;

#ifdef _DEBUG
	printf("dim: %d %d\n", m_dim[0], m_dim[1]);
	printf("orig: %f %f\n", m_orig[0], m_orig[1]);
	printf("span: %f %f\n", m_span[0], m_span[1]);
#endif
	m_type = t;

	switch(m_type) {
	case UCHAR_TYPE:
		p_data.uchar_ptr = new unsigned char[nverts];
		memcpy(p_data.uchar_ptr, array, sizeof(char)*nverts);
		break;
	case USHORT_TYPE:
		p_data.ushort_ptr = new unsigned short[nverts];
		memcpy(p_data.ushort_ptr, array, sizeof(short)*nverts);
		break;
	case INT_TYPE:
		p_data.int_ptr = new int[nverts];
		memcpy(p_data.int_ptr, array, sizeof(int)*nverts);
		break;
	case FLOAT_TYPE:
		p_data.float_ptr = new float[nverts];
		memcpy(p_data.float_ptr, array, sizeof(float)*nverts);
		break;
	case DOUBLE_TYPE:
		p_data.double_ptr = new double[nverts];
		memcpy(p_data.double_ptr, array, sizeof(double)*nverts);
		break;
	default:
		error("Unknown Data Type");
	}
	p_grad = NULL;

	init();

#ifdef _DEBUG
	printf("xbits %d, ybits %d\n", xbits, ybits);
	printf("yshift %d\n", yshift);
	printf("xmask %d\n", xmask);
	printf("ymask %d\n", ymask);
#endif

}

