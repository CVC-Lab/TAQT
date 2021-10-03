#include <betti.h>
#include <iostream>
#include <Reg2DConExtractor.h>
#include <contour2d.h>

using namespace std;

#define DIM 64

int main(int argc, char* argv[])
{
	if(argc < 4) {
		cerr << "Usage: " << argv[0] << " <rawiv file>"  << " <output file>" << " <isoval>" << endl;
		exit(1);
	}
	int mode = 2;
	if(argc > 4) {
		mode = atoi(argv[4]);
	}
	char ctree_fname[256];
	AugmentedContourTree *p_actree = NULL;
	sprintf(ctree_fname, "%s.ct\0", argv[1]);

	// create a simple 2D function and test 2D contour tree
	if(mode == 2) {	// 2D function
		Reg2Data* p_d2d = new Reg2Data(argv[1]);
		/*int dim[2] = {DIM, DIM};
		float *vals = new float[DIM*DIM];
		for(int j = 0; j < DIM; j++) {
		for(int i = 0; i < DIM; i++) {
		vals[i+j*DIM] = (i*i + j*j) / 10.0f;
		}
		}
		Reg2Data* p_d2d = new Reg2Data(dim, vals);
		p_d2d->maxext[2] = 0;
		p_d2d->writeRawiv("circle.rawiv2");*/

		Volume2DCritical* p_vol = new Volume2DCritical(p_d2d);

		ContourTree *p_tree = p_vol->computeContourTree();

		printf("reducing contour tree ...\n");

		p_tree->reduce();


		p_tree->dump(ctree_fname);

		delete p_tree;

		// Extracting 2D contours
		Reg2DConExtractor reg2Con(p_d2d);
		float isoval = 0.1f;
		if(argc > 3) {
			isoval = atof(argv[3]);
		}
		Contour2D* p_con = reg2Con.extractContour2D(isoval);
		p_con->write(argv[2]);
		//p_con->write("circle.c2d");
		delete p_d2d;
		delete p_con;
	}
	else if(mode == 3) {
		Reg3Data* p_fun = NULL;
		Reg3Data* p_pot = NULL;
		p_fun = new Reg3Data(argv[1]);

		VolumeReg3Critical* p_vol = new VolumeReg3Critical(p_fun, p_pot);
		ContourTree* p_ctree = p_vol->computeContourTree();
		printf("reducing contour tree ...\n");

		p_ctree->reduce();
		p_ctree->dump(ctree_fname);

		delete p_ctree;
		delete p_vol;
	}
}

