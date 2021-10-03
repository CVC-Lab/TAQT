#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "token.h"
#include "lex.yy.c"

#include <bio/bio.h>

int nd, nc, dim[3];
float span[3], orig[3];
FILE *ifp;
DiskIO *out;

static int parse_str(char* str);
static int read_data();

int main(int argc, char* argv[])
{
	if(argc < 3) {
		fprintf(stderr, "%s <input dx file> <output rawiv file>\n", argv[0]);
		exit(1);
	}

	if((ifp = freopen(argv[1], "r", stdin)) == NULL) {
		fprintf(stderr, "cannot open input dx file %s\n", argv[1]);
		exit(1);
	}

	out = new BufferedIO(argv[2], DiskIO::WRITE);	
	if(!out->open()) {
		fprintf(stderr, "cannot open output file %s\n", argv[2]);
		exit(1);
	}

	char cbuf[512];

	nd = nc = 0;
	span[0] = span[1] = span[2] = 0;

	while(nd == 0) {
		int res = yylex();
		if(res == 0) {
			fclose(ifp);
			return 0;
		} else if(res == ID) {
			parse_str(yylval);
		}

	}

	fclose(ifp);
	return 0;
}

int parse_str(char* str)
{
	if(strcmp(str, "counts") == 0) {
		if(nc == 0) {
			for(int i = 0; i < 3; i++) {
				int res = yylex();
				dim[i] = atoi(yylval);
				nc = 1;
			}
		}
		printf("dimension: %d %d %d\n", dim[0], dim[1], dim[2]);
	} else if (strcmp(str, "delta") == 0) {
		for(int i = 0; i < 3; i++) {
			int res = yylex();
			span[i] = (span[i] < atof(yylval))? atof(yylval):span[i];
		}
		printf("span: %f %f %f\n", span[0], span[1], span[2]);
	} else if(strcmp(str, "follows") == 0) {
		read_data();
		nd = 1;
	} else if(strcmp(str, "origin") == 0) {
		for(int i = 0; i < 3; i++) {
			int res = yylex();
			orig[i] = atof(yylval);
		}
		printf("origin: %f %f %f\n", orig[0], orig[1], orig[2]);
	} else {
		/* ingore */
	}

	return 1;
}

int read_data()
{
	float *data;
	int i, j, k;

	data = (float *) malloc(sizeof(float)*dim[0]*dim[1]*dim[2]);
	for(i = 0; i < dim[0]; i++) {
		for(j = 0; j < dim[1]; j++) {
			for(k = 0; k < dim[2]; k++) {
				int n = i + j*dim[0] + k*dim[0]*dim[1];
				int res = yylex();
				data[n] = (float)atof(yylval);
			}
		}
	}

	float minext[3], maxext[3];
	for(i = 0; i < 3; i++) {
		minext[i] = orig[i];
		maxext[i] = orig[i] + (dim[i]-1)*span[i];
	}
	out->put(minext, 3);
	out->put(maxext, 3);
	
	int nverts = dim[0]*dim[1]*dim[2];
	int ncells = (dim[0]-1)*(dim[1]-1)*(dim[2]-1);
	out->put(&nverts, 1);
	out->put(&ncells, 1);
	out->put(dim, 3);
	out->put(orig, 3);
	out->put(span, 3);
	
	
	out->put(data, dim[0]*dim[1]*dim[2]);
	out->close(0);
	return 1;
}
