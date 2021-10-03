#ifndef _SUPER_ARC_H
#define _SUPER_ARC_H

#include <iostream>

#include "basic.h"

using namespace std;

/**
 *	
 */
class BettiNumber {
public:
	BettiNumber() {
		b0 = b1 = b2 = 0;
	}

	const BettiNumber& operator += (const BettiNumber& bn) {
		b0 += bn.b0;
		b1 += bn.b1;
		b2 += bn.b2;
		return *this;
	}

	/*
	 *	Similarity between two set of Betti numbers 
	 */
	float simil(const BettiNumber& betti) {
		float s0, s1, s2;

		s0 = (MAX(b0, betti.b0) == 0)? 1.0f : MIN(b0, betti.b0) / (float)MAX(b0, betti.b0);
		s1 = (MAX(b1, betti.b1) == 0)? 1.0f : MIN(b1, betti.b1) / (float)MAX(b1, betti.b1);
		s2 = (MAX(b2, betti.b2) == 0)? 1.0f : MIN(b2, betti.b2) / (float)MAX(b2, betti.b2);

		return 0.5f*s0 + 0.3f*s1 + 0.2f*s2;
	}

	friend ostream& operator << (ostream& os, const BettiNumber& bn) {
		os << bn.b0 << " " << bn.b1 << " " << bn.b2 << " ";
		return os;
	}

	friend istream& operator >> (istream& is, BettiNumber& bn) {
		is >> bn.b0;
		is >> bn.b1;
		is >> bn.b2;
		return is;
	}
	
	int b0, b1, b2;
};

/*
 *	
 */
class SuperArc
{
public:
	SuperArc();

	SuperArc(int _v1, int _v2, int _xe, int _be);
	
	~SuperArc(void) {}

	int betti1() {
		return (1-xe) + ((be == 0)? 1:0);
	}

	int betti2() {
		return ((be == 0)? 1:0);
	}
	
	void calcBettiNumber() {
		//b2 = (be == 0)? 1:0;
		//b1 = 1 + b2 - xe;
	}

	BettiNumber getBettiNumbers();
	
	void addIntraNode(int nid) {
		intra_node = nid;
	}

	bool hasIntraNode() {
		return (intra_node >= 0);
	}

	int v1, v2;			// SuperNodes of the Arc
	int xe;				// Euler Characteristic
	int be;				// Number of boundary edges
	//int b1;			// Betti Number 1
	//int b2;			// Betti Number 2;
	int id;				// super arc ID
	int intra_node;		// Internal critical point (if there is one);
};
#endif


