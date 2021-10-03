#include "nnstruct.h"

int compareNNeighbors(const void *p1, const void *p2)
{
	NNeighbor ne1 = *((NNeighbor*) p1);
	NNeighbor ne2 = *((NNeighbor*) p2);
	if(ne1 < ne2) return -1;
	else if(ne2 < ne1) return 1;
	return 0;
}

