#ifndef UNIONFIND_H
#define UNIONFIND_H

#include<stdio.h>

class DisjointSet
{
public:
    /**
     * Create an empty disjoint set
     */
	DisjointSet(int _size) {
		size = _size;
		p = new int[size];
		rank = new int[size];
		max = new int[size];
		
		int i;
		for (i = 0 ; i<size ; i++) {
			p[i] = -1; rank[i]=0; LowestVtxArray[i]=-1;
		}
	}

	~DisjointSet(void) {
		delete[] p;
		delete[] rank;
		delete[] max;
	}

	/**
	 * Link two sets into a new set.
	 * @return the canonical element of the new set
	 */
    void link_0(int x, int y) {
		if (rank[x] > rank[y]) p[y] = x;
		else p[x] = y;
		if (rank[x] == rank[y]) 
			rank[y]++;
	}

	void link(int x , int y) {
		link_0(find(x) , find(y));
	}

	/**
	 * Add an elment x to the set with parent p
	 */
	int addElement(int x, int p) {

	}

	/**
	 * Create a new set with element t
	 */
	void makeSet(int t);

	/**
	 * @return -1 if t is not in the sets
	 */
    int find(int t);

	int getMaxNode(int x);
	void setMaxNode(int x, int p);

private:
	int* p;
	int* rank;
	int size;
	int* max;

	struct DJSNode {
		DJSNode() {
			rank = 0;
			p = -1;
			max = -1;
		}

		int rank;		// union heurist parameter
		int p;			// conanical member of the set	
		int max;		// the maximum node in the set
	};
	dynamic_array<DJSNode> nodes;

/*
#ifdef NO_HASHMAP
	map<int, DJSNode> nodes;
#else
	hash_map<int, DJSNode> nodes;
#endif*/

};

class UnionFind
{
	int* p;
	int* rank;
	//float* val;
	int size;
	int* LowestVtxArray;

public :
	UnionFind(int Size) 
	{ 
		size = Size;
		p = new int[Size];
		rank = new int[Size];
		LowestVtxArray = new int[Size];
	}

	void Clean()
	{
		int i;
		for (i = 0 ; i<size ; i++) {
			p[i]=0; rank[i]=0; LowestVtxArray[i]=0;
		}
	}


	~UnionFind() 
	{
		delete p;
		delete rank;
		delete LowestVtxArray;
	}

	void MakeSet(int x) 
	{
		p[x] = x ;
		rank[x] = 0;
	}
	
	void Union(int x , int y)
	{
		Link(FindSet(x) , FindSet(y));
	}
	
	void Link(int x , int y)
	{
		if (rank[x] > rank[y]) p[y] = x;
		else p[x] = y;
		if (rank[x] == rank[y]) 
			rank[y]++;
	}

	int FindSet(int x)
	{
		if (x != p[x])
			p[x] = FindSet(p[x]);
		return p[x];
	}

	void LowestVertex(int v , int vid)
	{	
		LowestVtxArray[FindSet(v)] = vid; 
	}

	int getLowestVertex(int v)
	{
		return LowestVtxArray[FindSet(v)];
	}

	void HighestVertex(int v , int vid)
	{
		LowestVtxArray[FindSet(v)] = vid; 
	}

	int getHighestVertex(int v)
	{
		return LowestVtxArray[FindSet(v)];
	}
};
#endif
