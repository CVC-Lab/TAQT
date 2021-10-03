#define STANDALONE

#include <XmlRpc.h>

#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>

#include <time.h>
#ifdef _WIN32
         #include <sys/types.h>
         #include <sys/timeb.h>
#else
         #include <sys/time.h>
#endif

#include <betti.h>

using namespace std;

// The server
XmlRpc::XmlRpcServer s;


float sim_w[6];

//Get the current time in seconds as a double value 
double getTime()
{
	#ifdef _WIN32
		 time_t ltime;
		 _timeb tstruct;
		 time( &ltime );
		 _ftime( &tstruct );
		 return (double) (ltime + 1e-3*(tstruct.millitm));
	#else
		 struct timeval t;
		 gettimeofday( &t, NULL );
		 return (double)(t.tv_sec + 1e-6*t.tv_usec);
	#endif
}

MultiConTree* buildMultiCon(const char* fname, int level, int dim[3], const char* potfile, float gmin, 
							float gmax, VolumeReg3Critical* &p_vol);

MultiConTree* buildMultiCon_XmlRpc(const char* fname, int level, float gmin, float gmax);


float matchRawiv(const char* fname1, const char* fname2)
{
	double start, end;
	start = getTime();
	
	sim_w[0] = 0.5;
	sim_w[1] = 0.5;
	sim_w[2] = 0;
	sim_w[3] = 0;
	sim_w[4] = 0;
	sim_w[5] = 0;
	

	int dim[3] = {8, 8, 8};
	int level = 16;
	printf("matching %s %s\n", fname1, fname2);
	if(fname1 == NULL || fname1 == NULL) {
		fprintf(stderr, "input file is NULL\n");
		exit(1);
	}
//	if(strcmp (fname1, fname2) == 0) return 1.0;

	char *pfile1 = NULL, *pfile2 = NULL;
	float gmin = 1, gmax = -1;

	VolumeReg3Critical *p_vol1, *p_vol2;
	MultiConTree *multicon1 = buildMultiCon(fname1, level, dim, pfile1, gmin, gmax, p_vol1);
	MultiConTree *multicon2 = buildMultiCon(fname2, level, dim, pfile2, gmin, gmax, p_vol2);

	float score = multicon1->matching(multicon2);
	printf("matching score = %f\n", score);
	FILE *fp;
	fp = fopen("score", "w");
	
	fprintf(fp, "%f", score);
	fclose(fp);

	end = getTime();
	//float t = end.tv_sec-start.tv_sec + (end.tv_usec-start.tv_usec)/1000000.0f;
	fprintf(stderr, "\nelapsed time = %f seconds\n", end-start);
	
	delete multicon1;
	delete multicon2;
	delete p_vol1;
	delete p_vol2;

	return score;	
}

float matchRawiv(const char* fname1, const char* fname2, float fmin, float fmax, char* potFile1, char* potFile2)
{
/*
	if(argc < 3) {
		fprintf( stderr, "Wrong arguments.\n" );
		fprintf( stderr, "Usage: %s <data file 1> <data file 2> [<fmin> <fmax> <pot file 1> <pot file 2>] \n", argv[0]);
		return( 1 );                                                                                                    
		                                                                                                                
                                                                                      
	}
*/	
	double start, end;
	start = getTime();
	
	sim_w[0] = 0.5;
	sim_w[1] = 0.5;
	sim_w[2] = 0.0;
	sim_w[3] = 0.0;
	sim_w[4] = 0.0;
	sim_w[5] = 0.0;
	

	int dim[3] = {8, 8, 8};
	int level = 16;
	printf("matching %s %s\n", fname1, fname2);

	char *pfile1 = NULL, *pfile2 = NULL;
	float gmin = 1, gmax = -1;
	gmin = fmin;
	gmax = fmax;
	

	pfile1 = potFile1;
	pfile2 = potFile2;
	if(pfile1 != NULL && pfile2 != NULL) {
		sim_w[0] = 0.03;
		sim_w[1] = 0.08;
		sim_w[2] = 0.44;
		sim_w[3] = 0.21;
		sim_w[4] = 0.1;
		sim_w[5] = 0.14;
	}

	/*
	if(argc > 12) {
		for(int i = 0; i < 6; i++) {
			sim_w[i] = atof(argv[i+7]);
		}
	}
	*/
	VolumeReg3Critical *p_vol1, *p_vol2;
	MultiConTree *multicon1 = buildMultiCon(fname1, level, dim, pfile1, gmin, gmax, p_vol1);
	MultiConTree *multicon2 = buildMultiCon(fname2, level, dim, pfile2, gmin, gmax, p_vol2);

	float score = multicon1->matching(multicon2);
	printf("matching score = %f\n", score);
	FILE *fp;
	fp = fopen("score", "w");
	
	fprintf(fp, "%f", score);
	fclose(fp);

	end = getTime();
	//float t = end.tv_sec-start.tv_sec + (end.tv_usec-start.tv_usec)/1000000.0f;
	fprintf(stderr, "\nelapsed time = %f seconds\n", end-start);
	
	delete multicon1;
	delete multicon2;
	delete p_vol1;
	delete p_vol2;

	return score;	
}

double matchDct(const char* fname1, const char* fname2)
{
	int port = 8500;

/*	if(argc < 3) {
		fprintf( stderr, "Wrong arguments.\n" );
		fprintf( stderr, "Usage: %s <dct file 1> <dct file 2> [<fmin> <fmax>] \n", argv[0]);
		return( 1 );
	}
*/	
	double start, end;
	start = getTime();
	
	sim_w[0] = 0.5;
	sim_w[1] = 0.5;
	sim_w[2] = 0;
	sim_w[3] = 0;
	sim_w[4] = 0;
	sim_w[5] = 0;
	

	int dim[3] = {8, 8, 8};
	int level = 16;
	printf("matching %s %s\n", fname1, fname2);

	float gmin = 0, gmax = 1;
	float score = 1.0;
	MultiConTree *multicon1 = buildMultiCon_XmlRpc(fname1, level, gmin, gmax);
	if(multicon1==NULL) { return score=0.0; printf("inputPDB doesn't exist\n");}
	MultiConTree *multicon2 = buildMultiCon_XmlRpc(fname2, level, gmin, gmax);
	if(multicon2==NULL) { return score=0.0; printf("inputPDB doesn't exist\n");}
	
	score = multicon1->matching(multicon2);
	printf("matching score = %f\n", score);
	
	FILE *fp;
	fp = fopen("score", "w");
	
	fprintf(fp, "%f", score);
	fclose(fp);

	end = getTime();
	//float t = end.tv_sec-start.tv_sec + (end.tv_usec-start.tv_usec)/1000000.0f;
	fprintf(stderr, "\nelapsed time = %f seconds\n", end-start);

	delete multicon1;
	delete multicon2;

	return score;	
}


using namespace XmlRpc;

// execute match param: argc, *argv[]
class executeDCTMatch : public XmlRpcServerMethod
{
public:
	executeDCTMatch(XmlRpcServer* s) : XmlRpcServerMethod("executeDCTMatch", s) {}
	
	void execute(XmlRpcValue& params, XmlRpcValue& result)
	{
		fprintf(stderr, "executeDCTMatch is called\n");
		
		string fname1 = params[0];
		string fname2 = params[1];
		
		double score = double(matchDct(fname1.c_str(), fname2.c_str()));
		fprintf(stderr, "executeDCTMatch score is %f\n", score);
		result = double(score);
	}
} executeDCTMatch(&s);

// execute match param: argc, *argv[]
class executeRAWIVMatch : public XmlRpcServerMethod
{
public:
	executeRAWIVMatch(XmlRpcServer* s) : XmlRpcServerMethod("executeRAWIVMatch", s) {}
	
	void execute(XmlRpcValue& params, XmlRpcValue& result)
	{
		fprintf(stderr, "executeRAWIVMatch is called\n");
		
		string fname1 = std::string(params[0]);
		string fname2 = std::string(params[1]);
		
		double score = double(matchRawiv(fname1.c_str(), fname2.c_str()));
		fprintf(stderr, "executeDCTMatch score is %f\n", score);
		result = double(score);
	}
} executeRAWIVMatch(&s);

// execute match param: argc, *argv[]
class executeClusterDCTMatch : public XmlRpcServerMethod
{
public:
	executeClusterDCTMatch(XmlRpcServer* s) : XmlRpcServerMethod("executeClusterDCTMatch", s) {}
	
	void execute(XmlRpcValue& params, XmlRpcValue& result)
	{
		fprintf(stderr, "executeClusterDCTMatch is called\n");
		int i=0,j=0;
		double score=0.0;
		int numOfResult = int(params[0]);
		string *fname = new string[numOfResult];
		for(i=0; i<numOfResult; i++) {
			fname[i] = std::string(params[i+1]);
		}
		
		int m = numOfResult+1;
		for(j=0; j<numOfResult; j++) {
			for(i=0; i<numOfResult; i++) {
				if(i==j) {
					result[j*m+i] = 1.0;
				}
				else if(i>j) {
					result[j*m+i]=result[i*m+j]=double(matchDct(fname[j].c_str(), fname[i].c_str()));
				}
			}
		}
		
	/*	
		for(j=0; j<numOfResult; j++) {
			for(i=0; i<numOfResult; i++) {
				printf("%f ", double(result[j*m+i]));
			}
			printf("\n");
		}
	*/	
	}
} executeClusterDCTMatch(&s);



int main(int argc, char *argv[])
{
#ifdef STANDALONE
	fprintf( stderr, "Running standalone build\n");
	if(argc < 3) {
		fprintf( stderr, "Usage: %s <data file 1> <data file 2> <data type: potential or density> [<fmin> <fmax> <weights: 6 floats> <out file>]\n", argv[0]);		
		return( 1 );
	}
        double start, end;
	start = getTime();
	
	sim_w[0] = 0.03;
	sim_w[1] = 0.08;
	sim_w[2] = 0.44;
	sim_w[3] = 0.21;
	sim_w[4] = 0.1;
	sim_w[5] = 0.14;
	
	int dim[3] = {8, 8, 8};
	int level = 16;
	fprintf( stderr, "*************************\n");
	fprintf( stderr, "matching: %s %s\n", argv[1], argv[2]);

	bool potential_data = false;
	if( !strncmp( argv[3], "-potential", 9 ) )
		potential_data = true;

	char *pfile1 = NULL, *pfile2 = NULL;
	float gmin = 1, gmax = -1;
	if (argc > 4) {
		gmin = atof(argv[4]);
		gmax = atof(argv[5]);
	        fprintf( stderr, "min, max: %f %f\n", gmin, gmax);
	}	
	
	if( !potential_data )
	{
		sim_w[0] = 0.5;
		sim_w[1] = 0.5;
		sim_w[2] = 0;
		sim_w[3] = 0;
		sim_w[4] = 0;
		sim_w[5] = 0;
	}

	// read weight
	if(argc > 11) {
		for(int i = 0; i < 6; i++) {
			sim_w[i] = atof(argv[i+6]);
		}
	}
	
	VolumeReg3Critical *p_vol1, *p_vol2;
	MultiConTree *multicon1 = buildMultiCon(argv[1], level, dim, NULL, gmin, gmax, p_vol1);
	MultiConTree *multicon2 = buildMultiCon(argv[2], level, dim, NULL, gmin, gmax, p_vol2);

	float score = multicon1->matching(multicon2);

	printf("weights: %f %f %f %f %f %f\n", sim_w[0], sim_w[1], sim_w[2], sim_w[3], sim_w[4], sim_w[5] );	
	printf("matching score = %f\n", score);
	FILE *fp;
	if(argc > 12) {
		fp = fopen(argv[13], "w");
	} else {
		fp = fopen("score", "w");
	}
	fprintf(fp, "%f", score);
	fclose(fp);

	end = getTime();
	//float t = end.tv_sec-start.tv_sec + (end.tv_usec-start.tv_usec)/1000000.0f;
	fprintf(stderr, "\nelapsed time = %f seconds\n", end-start);
	
	/* float *c1 = fvector(0, 2);
	float *c2 = fvector(0, 2);
	float **ax1 = matrix(0, 2, 0, 2);
	float **ax2 = matrix(0, 2, 0, 2);
	multicon1->getOrientation(0, 0, c1, ax1);
	multicon2->getOrientation(0, 0, c2, ax2);

	printf("c1 = (%f, %f, %f), c2 = (%f %f %f)\n", c1[0], c1[1], c1[2], c2[0], c2[1], c2[2]);
	for (int i = 0; i < 3; i++) {
		printf("R1: 	%f\t%f\t%f\n", ax1[i][0], ax1[i][1], ax1[i][2]);
	}
	for (int i = 0; i < 3; i++) {
		printf("R2: 	%f\t%f\t%f\n", ax2[i][0], ax2[i][1], ax2[i][2]);
	} */	

	/*double PI = 3.1415926;
	float cx = 0;
	float **mtx = matrix(0, 2, 0, 2);
	int RES = 40;
	for (int i = 0; i <= RES; i++) {
		for (int j = 0; j <= RES; j++) {
			for (int k = 0; k <= RES; k++) {
				float s = i / ((float)RES);
				float t1 = 2*PI*j/RES;
				float t2 = 2*PI*k/RES;
				float d1 = sqrt(1-s);
				float d2 = sqrt(s);
				float w = cos(t2)*d2;
				float x = sin(t1)*d1;
				float y = cos(t1)*d1;
				float z = sin(t2)*d2;
				float wx, wy, wz, xx, yy, yz, xy, xz, zz, x2, y2, z2; 	

				// calculate coefficients
				x2 = x + x; y2 = y + y; 
				z2 = z + z;
				xx = x * x2; xy = x * y2; xz = x * z2;
				yy = y * y2; yz = y * z2; zz = z * z2;
				wx = w * x2; wy = w * y2; wz = w * z2;

				ax2[0][0] = 1.0 - (yy + zz); 
				ax2[0][1] = xy - wz;
				ax2[0][2] = xz + wy; 
				ax2[1][0] = xy + wz; 
				ax2[1][1] = 1.0 - (xx + zz);
				ax2[1][2] = yz - wx; 
				ax2[2][0] = xz - wy; 
				ax2[2][1] = yz + wx;
				ax2[2][2] = 1.0 - (xx + yy); 

				float carbo = CarboIndex(*p_vol1, ax1, c1, *p_vol2, ax2, c2);
				if( carbo > cx) {
					cx = carbo;
					for(int ii = 0; ii < 3; ii++) {
						for(int jj = 0; jj < 3; jj++) {
							mtx[ii][jj] = ax2[ii][jj];
						}
					}
				}
			}
		}
	}
	printf("largest carbo index = %.8f\n", cx);
	for(int i = 0; i < 3; i++) {
		printf("%f\t%f\t%f\n", mtx[i][0], mtx[i][1], mtx[i][2]);
	}
	free_matrix(mtx, 0, 2, 0, 2);*/	

	delete multicon1;
	delete multicon2;
	delete p_vol1;
	delete p_vol2;
#else	
	fprintf( stderr, "Running server build\n");
	int port = 8500;
	if(argc < 2) {
		fprintf( stderr, "Wrong arguments.\n" );
		fprintf( stderr, "Usage: %s <port number(8500)>\n", argv[0]);
		return( 1 );
	}
	port = atoi(argv[1]);
	
	// XML RPC Server
	// Create the server socket on the specified port
	s.bindAndListen(port);
		
	// Enable introspection
	s.enableIntrospection(true);
			
	// Wait for requests indefinitely
	int it = 0;
	fprintf(stderr, "Wait = %d\n", it++);
	s.work(-1.0);
#endif
	return 0;
}

MultiConTree* buildMultiCon(const char* fname, int level, int dim[3],  
			const char* potfile, float gmin, float gmax, 
			VolumeReg3Critical* &p_vol) {
	// read the first data
	Reg3Data* p_fun = NULL;
	Reg3Data* p_pot = NULL;
	p_vol = NULL;

	char ctree_fname[256], dual_fname[256];
	AugmentedContourTree *p_actree = NULL;
#ifdef STANDALONE
	sprintf(ctree_fname, "%s.ctr\0", fname);
	sprintf(dual_fname, "%s.dct\0", fname);

	printf("ctree: %s\n", ctree_fname);
	printf("dual: %s\n", dual_fname);
#else
	printf("length = %d\n", strlen(fname));
	
	// June added /////////////////////////////////////////////////////
	char* pch;
	int index = 0;
	pch = strrchr((char *)fname, '/');
	index = pch-fname+1;
//	printf("found at %d\n", index);
	
	
	char fnametmp[256];
	int i;
	for(i = 0; i < strlen(fname)-index; i++) {
//		printf("copy %c\n", fname[index+i]);
		fnametmp[i] = fname[index+i];
	}
	fnametmp[i] ='\0';
	printf("%s\n", fnametmp);
	///////////////////////////////////////////////////////////////////
	if (potfile) {
//		sprintf(ctree_fname, "/home/work/junenim/MOBIOS/VirtualScreen/ChemCompoundMactData/CtrFiles/DX_%s.ctr\0", fnametmp);
//		sprintf(dual_fname, "/home/work/junenim/MOBIOS/VirtualScreen/ChemCompoundMactData/DctFiles/DX_%s.dct\0", fnametmp);
		//sprintf(ctree_fname, "/home/work/junenim/MOBIOS/GeoMactData/CtrFiles/DX_%s.ctr\0", fnametmp);
		//sprintf(dual_fname, "/home/work/MOBIOS/GeoMactData/DctFiles/DX_%s.dct\0", fnametmp);
		sprintf(ctree_fname, "/share/work/junenim/MOBIOS/MactData/CtrFiles/DX_%s.ctr\0", fnametmp);
		sprintf(dual_fname, "/share/work/junenim/MOBIOS/MactData/DctFiles/DX_%s.dct\0", fnametmp);
	}
	else {
//		sprintf(ctree_fname, "/home/work/junenim/MOBIOS/VirtualScreen/ChemCompoundMactData/CtrFiles/%s.ctr\0", fnametmp);
//		sprintf(dual_fname, "/home/work/junenim/MOBIOS/VirtualScreen/ChemCompoundMactData/DctFiles/%s.dct\0", fnametmp);
		//sprintf(ctree_fname, "/home/work/junenim/MOBIOS/GeoMactData/CtrFiles/%s.ctr\0", fnametmp);
		//sprintf(dual_fname, "/home/work/junenim/MOBIOS/GeoMactData/DctFiles/%s.dct\0", fnametmp);
		sprintf(ctree_fname, "/share/work/junenim/MOBIOS/MactData/CtrFiles/%s.ctr\0", fnametmp);
		sprintf(dual_fname, "/share/work/junenim/MOBIOS/MactData/DctFiles/%s.dct\0", fnametmp);
	}
	
//	sprintf(ctree_fname, "/home/work/junenim/MOBIOS/TAQT/%s.ctr\0", fnametmp);
//	sprintf(dual_fname, "/home/work/junenim/MOBIOS/TAQT/%s.dct\0", fnametmp);
#endif
	FILE *fp, *dfp;
	DualGraph *dual;
	float global_min, global_max;

	printf("\nbuildMultiCon: Dual contour tree of %s\n", fname);

	p_fun = new Reg3Data(fname, dim);

	if (gmin < gmax) {
		global_min = gmin;
		global_max = gmax;
	} else {
		float min, max;
		p_fun->getFuncMinMax(min, max);
		global_min = min + 0.000001;
		global_max = max + 0.000001;
	}
	printf("fmin = %f, fmax = %f\n", global_min, global_max);

	DualNode::setGlobalMinMax(global_min, global_max);

	p_fun = new Reg3Data(fname, dim);
	if (potfile) {
			p_pot = new Reg3Data(potfile, dim);
	} else {
		p_pot = new Reg3Data(fname, dim);
		for(int i = 0; i < p_pot->getNVerts(); i++) {
			p_pot->setValue(i, 1.0f);
		} 
	}
	p_vol = new VolumeReg3Critical(p_fun, p_pot);
	
//	fprintf(stderr, "\n\n dual Min,Max: (%f, %f)\n",global_min, global_max);
//	fprintf(stderr, "p_vol->getVolume(): %f\n\n", p_vol->getVolume());
					 
	dual = new DualGraph(global_min, global_max, level, p_vol->getVolume());  
	dual->setVolume(p_vol->getVolume());
	
	if (!dual->loadDualFile(dual_fname)) {
		if ((fp = fopen(ctree_fname, "r")) != NULL) {
			fclose(fp);
			p_actree = new AugmentedContourTree();
			p_actree->load(ctree_fname);
		} else {
			float min, max;
			p_fun->getFuncMinMax(min, max);
			fprintf(stderr, "function min = %f, max = %f\n", min, max);

			ContourTree *p_stree = p_vol->splitTree();
			ContourTree *p_jtree = p_vol->joinTree();
			ContourTree *p_tree = ContourTree::mergeTree(p_jtree, p_stree);
			delete p_stree;
			delete p_jtree;

			p_actree = p_tree->augment();
			delete p_tree;

			printf("reducing contour tree ...\n");
			p_actree->reduce();
			//p_actree->print();
			p_actree->done();
			p_actree->dump(ctree_fname);
		}
		//float min, max;
		//p_actree->getMinMax(min, max);
		//fprintf(stderr, "func min = %f, max = %f\n", min, max);
		//perturb min and max
		//min -= (max-min)*0.000013;
		//max += (max-min)*0.000013;
		//min = global_min;
		//max = global_max;
		DualNode::setGlobalMinMax(global_min, global_max);
		dual->build(p_vol, p_actree);
		dual->dump(dual_fname);
		delete p_actree;
		//delete p_vol;
		//delete p_pot;
	}
	
	fprintf(stderr, "\t\tdual->getVolume(): %f\n\n", dual->getTotalVolume());
	MultiConTree *multicon = new MultiConTree(dual, dual->getTotalVolume());
//	delete p_vol;

	return multicon;
}

MultiConTree* buildMultiCon_XmlRpc(const char* fname, int level, float gmin, float gmax) {

	FILE *fp, *dfp;
	DualGraph *dual;
	float global_min, global_max;

	printf("\nbuildMultiCon_XmlRpc: Dual contour tree of %s\n", fname);

	global_min = gmin;
	global_max = gmax;
	//printf("fmin = %f, fmax = %f\n", global_min, global_max);

	DualNode::setGlobalMinMax(global_min, global_max);

	dual = new DualGraph(global_min, global_max, level);  
		
	if(!dual->loadDualFile(fname)) {
		fprintf(stderr, "%s dual contour tree does NOT exist\n", fname);
		return 0;
		//exit(1);
	}
	
	MultiConTree *multicon = new MultiConTree(dual, dual->getVolume());

	return multicon;
}


