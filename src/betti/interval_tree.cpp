#include "interval_tree.h"
#include "actree.h"

#include <cstdio>
#include <cstring>

using namespace std;
#ifdef _WIN32
//template <class T> vector<IntervalNode<T> >* MinLess<T>::p_nodes; 
//template <class T> vector<IntervalNode<T> >* MaxGreater<T>::p_nodes; 
vector<IntervalNode<int> >* MinLess<int>::p_nodes;
vector<IntervalNode<int> >* MaxGreater<int>::p_nodes;

vector<IntervalNode<LL_Node<SuperArc>* > >* MinLess<LL_Node<SuperArc>* >::p_nodes;
vector<IntervalNode<LL_Node<SuperArc>* > >* MaxGreater<LL_Node<SuperArc>* >::p_nodes;

#else

// arand: hacking this on 3/15/2012

//template<int> vector<IntervalNode<int> >* MinLess<int>::p_nodes;
//template<int> vector<IntervalNode<int> >* MaxGreater<int>::p_nodes;

template<class T> vector<IntervalNode<T> >* MinLess<T>::p_nodes;
template<class T> vector<IntervalNode<T> >* MaxGreater<T>::p_nodes;

//template vector<IntervalNode<int> >* MinLess<int>::p_nodes;
//template vector<IntervalNode<int> >* MaxGreater<int>::p_nodes;

//template<> vector<IntervalNode<int> >* MinLess<int>::p_nodes;
//template<> vector<IntervalNode<int> >* MaxGreater<int>::p_nodes;

template vector<IntervalNode<LL_Node<SuperArc>* > >* MinLess<LL_Node<SuperArc>* >::p_nodes;
template vector<IntervalNode<LL_Node<SuperArc>* > >* MaxGreater<LL_Node<SuperArc>* >::p_nodes;

#endif

