#include "arc.h"

SuperArc::SuperArc()
: intra_node(-1)
{
	v1 = v2 = 0;
	xe = be = 0;
	//b2 = 1;
	//b1 = 1 + b2 - xe;
}

SuperArc::SuperArc(int _v1, int _v2, int _xe, int _be)
: intra_node(-1)
{
	v1 = _v1;
	v2 = _v2;
	xe = _xe;
	be = _be;
	//b2 = (be == 0)? 1:0;
	//b1 = 1 + b2 - xe;
}

BettiNumber SuperArc::getBettiNumbers()
{
	BettiNumber bn;
	bn.b0 = 1;
	bn.b1 = betti1();
	bn.b2 = betti2();
	return bn;
}


