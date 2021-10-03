#include "Reg2DConExtractor.h"
#include "contour2d.h"

//------------------------------------------------------------------------
//
// rect2d - a local static variable
//
//------------------------------------------------------------------------

static int rect2d[16][5] = {
   {0},
   {1, 3, 0},
   {1, 0, 1},
   {1, 3, 1},
   {1, 1, 2},
   {2, 3, 2, 1, 0},
   {1, 0, 2},
   {1, 3, 2},
   {1, 2, 3},
   {1, 2, 0},
   {2, 0, 3, 2, 1},
   {1, 2, 1},
   {1, 1, 3},
   {1, 1, 0},
   {1, 0, 3},
   {0}
};

Reg2DConExtractor::Reg2DConExtractor(Reg2Data* _data)
{
	p_data = _data;
	p_EdgeMap = 0;
}

Reg2DConExtractor::~Reg2DConExtractor(void)
{
}

Contour2D* Reg2DConExtractor::extractContour2D(float isoval)
{
	if(p_data == NULL) return NULL;

	Contour2D* p_con = new Contour2D();
	p_con->setExtent(p_data->minext, p_data->maxext);

	unsigned int dim[2];
	int nCellX, nCellY;
	p_data->getDim(dim);
	nCellX = dim[0] - 1;
	nCellY = dim[1] - 1;

	float vals[4];
	p_EdgeMap = new map<int, int>();
	for(int j = 0; j < nCellY; j++) {
		for(int i = 0; i < nCellX; i++) {
			p_data->getCellValues(i, j, vals);
			int code = cellCode(vals, isoval);
			// iterate the num of conour edges in the cell
			for(int e=0; e<rect2d[code][0]; e++) { 
				int v1, v2;
				if(!isInEdgeCache(rect2d[code][2*e+1], i, j, v1)) { 
					v1 = interpEdge(rect2d[code][2*e+1], vals, isoval, i, j, p_con);
					addToCache(rect2d[code][2*e+1], i, j, v1);
				}
				if(!isInEdgeCache(rect2d[code][2*e+2], i, j, v2)) {
					v2 = interpEdge(rect2d[code][2*e+2], vals, isoval, i, j, p_con);
					addToCache(rect2d[code][2*e+2], i, j, v2);
				}
				p_con->addEdge(v1, v2);
			}
		}
	}

	delete p_EdgeMap;
	return p_con;
}

int Reg2DConExtractor::cellCode(float vals[4], float isoval)
{
	int code = 0;
	if (vals[0] < isoval) code += 0x01;
	if (vals[1] < isoval) code += 0x02;
	if (vals[2] < isoval) code += 0x04;
	if (vals[3] < isoval) code += 0x08;

	return code;
}

int Reg2DConExtractor::interpEdge(int edge, float *val, float isovalue, int i, int j, Contour2D* p_Con)
{
   float ival = 0.5f;
   float pt[2];

   switch (edge) {
      case 0:
         if(val[0] != val[1]) ival = (isovalue-val[1])/(val[0]-val[1]);
         pt[0] = (1.0f-ival) * p_data->xCoord(i+1) +
                       ival  * p_data->xCoord(i+0);
         pt[1] = p_data->yCoord(j);
         break;
      case 1:
         if(val[1] != val[2]) ival = (isovalue-val[2])/(val[1]-val[2]);
         pt[0] = p_data->xCoord(i+1);
         pt[1] = (1.0f-ival) * p_data->yCoord(j+1) +
                       ival  * p_data->yCoord(j+0);
         break;
      case 2:
         if(val[2] != val[3]) ival = (isovalue-val[3])/(val[2]-val[3]);
         pt[0] = (1.0f-ival) * p_data->xCoord(i+0) +
                       ival  * p_data->xCoord(i+1);
         pt[1] = p_data->yCoord(j+1);
         break;
      case 3:
         if(val[3] != val[0]) ival = (isovalue-val[0])/(val[3]-val[0]);
         pt[0] = p_data->xCoord(i);
         pt[1] = (1.0f-ival) * p_data->yCoord(j+0) +
                       ival  * p_data->yCoord(j+1);
         break;
   }
   return(p_Con->addVert(pt));
}

bool Reg2DConExtractor::isInEdgeCache(int edge, int i, int j, int &v)
{
	//return false;
	int edge_id = edgeID(edge, i, j);
	map<int, int>::iterator it = p_EdgeMap->find(edge_id);
	if(it !=  p_EdgeMap->end()) {
		v = (*it).second;
		return true;
	}
	return false;
}

void Reg2DConExtractor::addToCache(int edge, int i, int j, int v)
{
	int edge_id = edgeID(edge, i, j);
	(*p_EdgeMap)[edge_id] = v;
}

int Reg2DConExtractor::edgeID(int edge, int i, int j) {
	int ix, iy, dir;
	unsigned int dim[2];
	p_data->getDim(dim);
	int nEdgesX = (dim[0]-1)*dim[1];
	switch(edge) {
		case 0:
			ix = i;
			iy = j;
			dir = 0;
			break;
		case 1:
			ix = i+1;
			iy = j;
			dir = 1;
			break;
		case 2:
			ix = i;
			iy = j+1;
			dir = 0;
			break;
		case 3:
			ix = i;
			iy = j;
			dir = 1;
			break;
	}
	int eid;
	switch(dir) {
		case 0:
			eid = ix + iy *(dim[0] -1);
			break;
		case 1:
			eid = nEdgesX + iy + ix*(dim[1]-1);
			break;
	}
	return eid;
}