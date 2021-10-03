#pragma once
#include "Contour2DExtractor.h"
#include "reg2data.h"

#include <map>
using namespace std;

class Reg2DConExtractor :
	public Contour2DExtractor
{
public:
	Reg2DConExtractor(Reg2Data* _data);
	virtual ~Reg2DConExtractor(void);

	virtual Contour2D* extractContour2D(float isoval);

protected:
	Reg2Data* p_data;

	int cellCode(float vals[4], float isoval);

	///\r the the interploated vertex id
	int interpEdge(int edge, float *vals, float isovalue, int i, int j, Contour2D* p_Con);

	bool isInEdgeCache(int edge, int i, int j, int& v);

	void addToCache(int edge, int i, int j, int v);

private:
	map<int, int>* p_EdgeMap;

	int edgeID(int edge, int i, int j);
};
