#include <reg3data.h>
#include <stdio.h>

int main(int argc, char* argv[]) 
{
	if(argc < 3) {
		fprintf(stderr, "Usage: %s <Rawiv File> <NewVolume Config File>\n", argv[0]);
		exit(1);
	}
	Reg3Data *pData = new Reg3Data(argv[1]);
	printf("func min = %f, max = %f\n", pData->getFuncMin(), pData->getFuncMax());
	float min = pData->getFuncMin();
	float max = pData->getFuncMax();
	if(!(min < 0 && max > 0)) {
		fprintf(stderr, "ERROR: potential is not in an appropriate range\n");
		exit(1);
	}
	float zero = -min / (max - min);	
	// Write out the configure file
	FILE *fp = fopen(argv[2], "w");
	fprintf(fp, "Anthony and Vinay are Great.\n");
	fprintf(fp, "Alphamap\n");
	fprintf(fp, "Number of nodes\n%d\n", 9);
	fprintf(fp, "Position and opacity\n");
	fprintf(fp, "0 0.4\n");
	fprintf(fp, "1 0.4\n");
	fprintf(fp, "%f 0.3\n", zero-0.02);
	fprintf(fp, "%f 0.15\n", zero-0.01);
	fprintf(fp, "%f 0.001\n", zero-0.001);
	fprintf(fp, "%f 0\n", zero);
	fprintf(fp, "%f 0.001\n", zero+0.001);
	fprintf(fp, "%f 0.15\n", zero+0.01);
	fprintf(fp, "%f 0.3\n", zero+0.02);
	fprintf(fp, "ColorMap\n");
	fprintf(fp, "Number of nodes\n%d\n", 7);
	fprintf(fp, "Position and RGB\n");
	fprintf(fp, "0 %d %d %d\n", 0, 0, 1);
	fprintf(fp, "1 %d %d %d\n", 1, 0, 0);
	fprintf(fp, "%f %d %d %d\n", zero-0.02, 0, 0, 1);
	fprintf(fp, "%f 0.2 0.2 1\n", zero-0.01);
	fprintf(fp, "%f 1 1 1\n", zero, 1, 1, 1);
	fprintf(fp, "%f 1 0.2 0.2\n", zero+0.01, 1, 0.2, 0.2);
	fprintf(fp, "%f %d %d %d\n", zero+0.02, 1, 0, 0);
	fprintf(fp, "IsocontourMap\nNumber of nodes\n%d\n", 0);
	delete pData;
	return 0;
}

