/*****************************************************************************\
 *
 * pqueue.h -- 
 *
 *
 * Author:      Fausto Bernardini (fxb@cs.purdue.edu)
 *
 * Created - June 15, 1993
 * Ported to C++ by Raymund Merkert - June 1995
 * Changes by Fausto Bernardini - Sept 1995 
 *
\*****************************************************************************/

#ifndef __PQUEUE_H
#define __PQUEUE_H

//#include "basic.h"
#include "bin.h"


template<class T> class PriorityQueue;
template<class T> class SortedPriorityQueue;

template<class T> class Pqueuerec {
	friend class PriorityQueue<T>;
	friend class SortedPriorityQueue<T>;
private:
	// DATA
	double key;
	T item;
};


/** A templated priority queue.
  A priority queue based on a partially ordered tree.
  The number of #Pqueuerecs# that can be put into this queue is given as
  a parameter when the object is created ( or a default is used when no
  parameter is given).  This number does not change. */
template<class T>
class PriorityQueue {
public:
	/// Constructor.
	PriorityQueue( int bocksize=0 );

	/// Reinitialize.
	void cleanUp();

	/// Insert #item# into the queue with priority #priority#.
	void insert( T& item, double priority );

	/// Removes the element with highest priority from the queue.
	double extract(T& item);

	/// Return the element with highest priority from the queue.
	double max(T& item);

	/// Return #true# if the queus is empty, false otherwise.
	bool isEmpty();

	void remove(int);

	int numItems(void) { return(_q.numItems()); }

	T *nthEntry(int n) { return(&_q[n].item); }

protected:

	Pqueuerec<T> *nthItem(int n) { return(&_q[n]); }

	// Used by extract to keep the tree partially ordered.
	void sink(int i);

private:
	// DATA
	Bin< Pqueuerec<T> > _q;		// Bin that holds the Pqueuerecs
};


/*****************************************************************************\
 *
\*****************************************************************************/

template<class T>
inline PriorityQueue<T>::PriorityQueue( int blocksize )
{
	_q.setBlocksize(blocksize);
}


template<class T>
inline void PriorityQueue<T>::cleanUp()
{
	_q.cleanUp();
}


template<class T>
inline void PriorityQueue<T>::insert(T& new_item, double new_priority )
{
	int i, p;

	i = _q.numItems();
	_q.add();
	p = (i-1)/2;
	while( i > 0 && _q[p].key < new_priority ) {
		_q[i] = _q[p];
		i = p;
		p = (i-1)/2;
	}
	_q[i].item = new_item;
	_q[i].key = new_priority;
}


template<class T>
inline void PriorityQueue<T>::sink(int i)
{
	int l, r, max;

	while( 1 ) {
		l = 2*i+1; r = 2*i+2;
		max = ( (l < _q.numItems()) && (_q[l].key > _q[i].key) ? l : i );
		max = ( (r < _q.numItems()) && (_q[r].key > _q[max].key) ? r : max );
		if( max == i ) {
			break;
		} else {
			swap(_q[max], _q[i]);
			i = max;
		}
	}
}


template<class T>
inline double PriorityQueue<T>::extract(T& item)
{
	double p = _q[0].key;
	item = _q[0].item;

	_q[0] = _q[_q.numItems()-1];
	_q.removeLast();
	sink(0);
	return p;
}


template<class T>
inline void PriorityQueue<T>::remove(int n)
{
	_q[n] = _q[_q.numItems()-1];
	_q.removeLast();
	sink(n);
}

template<class T>
inline double PriorityQueue<T>::max(T& item)
{
	item = _q[0].item;
	return _q[0].key;
}


template<class T>
inline bool PriorityQueue<T>::isEmpty()
{
	return _q.numItems() == 0;
}


#endif

