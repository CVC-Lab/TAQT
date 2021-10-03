/***************************************************************************
                          interval_tree.h  -  description
                             -------------------
    begin                : Sun Oct 19 2003
    copyright            : (C) 2003 by Xiaoyu Zhang
    email                : xiaoyu@csusm.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 
#ifndef _INTERVAL_TREE_H
#define _INTERVAL_TREE_H

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <vector>
#include <algorithm>

using namespace std;

/**
 *
 */
template <class T>
class IntervalNode {
public:
	IntervalNode(const T& x, float _min, float _max) {
		min = _min;
		max = _max;
		obj = x;
	}

	float min, max;
	T obj;
};

template <class T>
class IntervalBucket {
public:
	///
	void insert(unsigned int n) {
		vec.push_back(n);
	}

	int size() {
		return (int)vec.size();
	}

	unsigned int get(int n) {
		return vec[n];
	}

	/**
	 * Sort the bucket in the increasing order of min 
	 * in the node vector. 
	 */
	void sortMin();

	/**
	 * Sort the bucket in the decreasing order of max 
	 * in the node vector. 
	 */
	void sortMax();

private:
	vector<unsigned int> vec;
};

template <class T>
class MinLess {
public:
	static void setNodePointer(vector<IntervalNode<T> >* _ptr) {
		p_nodes = _ptr;
	}
	bool operator() (const unsigned int& n1, const unsigned int& n2) {
		float x1 = (*p_nodes)[n1].min;
		float x2 = (*p_nodes)[n2].min;

		return x1 < x2;
	}

	static vector<IntervalNode<T> >* p_nodes;
 private:
};

template <class T>
class MaxGreater {
public:
	static void setNodePointer(vector<IntervalNode<T> >* _ptr) {
		p_nodes = _ptr;
	}
	bool operator() (const unsigned int& n1, const unsigned int& n2) {
		float x1 = (*p_nodes)[n1].max;
		float x2 = (*p_nodes)[n2].max;

		return x1 > x2;
	}

private:
	static vector<IntervalNode<T> >* p_nodes;
};

template <class T>
class IntervalTree {
public:
	IntervalTree(unsigned int n = 0, float *_vals = NULL);

	virtual ~IntervalTree(void);

	void init(unsigned int n, float *vals);

	void insert(const T& x, float min, float max);

	void done();

	/**
	 * query() -- traverse the tree, storing the cell id's of all
	 *            segments containing the given value in a list
	 * @return Vector of satisfying intervals
	 */
	vector<T>* query(float val);

protected:
	unsigned int addNode(const T& x, float min, float max) {
		unsigned int n = (unsigned int)nodes.size();
		IntervalNode<T> node(x, min, max);
		nodes.push_back(node);
		return n;
	}

	static int cmp_float(const void* p1, const void* p2);

	float getMin(int n) {
		return nodes[n].min;
	}

	float getMax(int n) {
		return nodes[n].max;
	}

private:
	unsigned int nleaf;	// number of splitting values
	float *vals;		// splitting values

	vector<IntervalNode<T> >nodes;		// Array of interval nodes
	IntervalBucket<T> *minlist;			// minimum list of each splitting node
	IntervalBucket<T> *maxlist;			// maximum list of each splitting node
};

template <class T>
void IntervalBucket<T>::sortMin()
{
	sort(vec.begin(), vec.end(), MinLess<T>());
}

template <class T>
void IntervalBucket<T>::sortMax()
{
	sort(vec.begin(), vec.end(), MaxGreater<T>());
}

template <class T>
IntervalTree<T>::IntervalTree(unsigned int n, float *_vals)
{
	if (n == 0) {
		nleaf = 0;
		vals = NULL;
		minlist = NULL;
		maxlist = NULL;
	} else {
		init(n, _vals);
	}
}

template <class T>
IntervalTree<T>::~IntervalTree(void)
{
	if (vals != NULL) delete[] vals;
	if (minlist != NULL) delete[] minlist;
	if (maxlist != NULL) delete[] maxlist;
}

template <class T>
int IntervalTree<T>::cmp_float(const void* p1, const void* p2)
{
	float x1 = *((float *)p1);
	float x2 = *((float *)p2);
	if (x1 < x2) return -1;
	if (x1 > x2) return 1;
	return 0;
}

template <class T>
void IntervalTree<T>::init(unsigned int n, float* _vals)
{
	nleaf = n;
	vals = new float[nleaf];
	memcpy(vals, _vals, sizeof(float)*nleaf);
	qsort(vals, nleaf, sizeof(float), IntervalTree<T>::cmp_float);
	minlist = new IntervalBucket<T>[nleaf];
	maxlist = new IntervalBucket<T>[nleaf];
}

template <class T>
void IntervalTree<T>::insert(const T& x, float min, float max)
{
	unsigned int left, right, root;
	unsigned int n;

	n = addNode(x, min, max);

	left = 0;
	right = nleaf-1;

#ifdef DEBUG_TREE
	printf("inserting node %d (%f %f)\n", n, min, max);
#endif
	while (left < right) {
		root = (left + right) >> 1;

#ifdef DEBUG_TREE
		printf("comparing with split value %f (node %d)\n", vals[root], root);
#endif

		if (min <= vals[root] && vals[root] <= max) {
			minlist[root].insert(n);
			maxlist[root].insert(n);
			return;
		}

		if (min > vals[root]) {
			left=root+1;
#ifdef DEBUG_TREE
			printf("left->root+1\n");
#endif
		} else { /* max < vals[root] */
			right=root-1;
#ifdef DEBUG_TREE
			printf("right->root-1\n");
#endif
		}
	}

	// left == right
	minlist[left].insert(n);
	maxlist[left].insert(n);
}

template <class T>
void IntervalTree<T>::done()
{
	MinLess<T>::setNodePointer(&nodes);
	MaxGreater<T>::setNodePointer(&nodes);

	for (unsigned int i = 0; i < nleaf; i++) {
		minlist[i].sortMin();
		maxlist[i].sortMax();
	}
}

template <class T>
vector<T>* IntervalTree<T>::query(float val)
{
	vector<T>* p_vec = new vector<T>;

	int left, right, root;
	int i;

	left = 0;
	right = nleaf-1;

	while (left < right) {
		root = (left + right) >> 1;

		if (vals[root] >= val) {
			// for all cells in minlist, we know max > val
			// search the minlist for all cells with min < val
			for (i = 0; i < minlist[root].size(); i++) {
				//if (getMin(minlist[root].get(i)) <= val && getMax(minlist[root].get(i)) > val) {
				if (getMin(minlist[root].get(i)) <= val) {
					p_vec->push_back(nodes[minlist[root].get(i)].obj);
				} else
					break;
			}
			right=root-1;
		} else {
			// for all cells in maxlist, we know min < val
			// search the maxlist for all cells with max > val
			for (i = 0; i < maxlist[root].size(); i++)
				if (getMax(maxlist[root].get(i)) > val)
					p_vec->push_back(nodes[maxlist[root].get(i)].obj);
				else
					break;
			left=root+1;
		}
	}

    return p_vec;    
}

#endif

