#include <iostream>
#include <stdlib.h>
#include <bio.h>
#include <reg3data.h>

/******
* this program uses the diskIO class to decimate a 3D regular data, which
* may be larger than main memory
*/

using namespace std;

typedef float dtype;

int main(int argc, char* argv[])
{
	int  i, j, k, fac[3];
	
	Reg3Data* p_data = new Reg3Data(argv[1]);

	// read decimation factor from command line
	for(i = 0; i < 3; i++) {
		fac[i] = atoi(argv[3+i]);
	}
	cout << "decimation factor:" << fac[0] << "," << fac[1] << "," << fac[2] << endl;

	// read header info from input
	float minext[3], maxext[3], orig[3], span[3];
	int nverts, ncells, dim[3];
	p_data->getOrig(orig);
	p_data->getSpan(span);
	p_data->getDim(dim);

	printf("dim: %d %d %d\n", dim[0], dim[1], dim[2]);
	printf("orig: %f %f %f\n", orig[0], orig[1], orig[2]);
	printf("span: %f %f %f\n", span[0], span[1], span[2]);

	int ndim[3];
	float nspan[3];
	for(i = 0; i < 3; i++) {
		ndim[i] = dim[i]/fac[i];
		nspan[i] = span[i]*fac[i];
	}

	dtype* dout = new dtype[ndim[0]*ndim[1]*ndim[2]];
	Reg3Data* p_dec = new Reg3Data(ndim, dout);
	p_dec->setSpan(nspan);
	p_dec->setOrig(orig);
	
	// write header info to output
	nverts = ndim[0]*ndim[1]*ndim[2];
	ncells = (ndim[0]-1)*(ndim[1]-1)*(ndim[2]-1);

	// write out data
	for(k = 0; k < ndim[2]; k++) {  
		for(j = 0; j < ndim[1]; j++) {
			for(i = 0; i < ndim[0]; i++) {
				int vid = i + j*ndim[0] + k*ndim[0]*ndim[1];
				p_dec->setValue(vid, p_data->getValue(i*fac[0], j*fac[1], k*fac[2]));
			}
		}
	}
	
	p_dec->writeRawiv(argv[2]);

	delete p_data;
	delete p_dec;

	return 0;
}
