#include "disjointset.h"

DisjointSet::DisjointSet(int _size)
{
	size = _size;
	nodes = new DJSNode[size];
}

DisjointSet::~DisjointSet()
{	
	delete[] nodes;
}











