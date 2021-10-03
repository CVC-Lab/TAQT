#include <assert.h>

#include "node.h"

SuperNode::SuperNode(int i) {
	id = i;
}

void SuperNode::addNeighbor(int n, int x) {
	neighbors.push_back(sn_pair(n, x));
}

int SuperNode::getNeighbor(int n) {
#ifdef _DEBUG
	assert(n < degree());
#endif
	return neighbors[n].nei;
}


int SuperNode::getTag(int n) {

	return neighbors[n].tag;

}



int SuperNode::getNeighborTag(int nid) {

	vector<sn_pair>::iterator it;

	for(it = neighbors.begin(); it != neighbors.end(); ++it) {

		if((*it).nei == nid) return (*it).tag;

	}

	assert(0);

	return -1;

}


void SuperNode::removeNeighbor(int n)
{
	vector<sn_pair>::iterator it;
	for(it = neighbors.begin(); it != neighbors.end(); ++it) {
		if((*it).nei == n) {
			neighbors.erase(it);

			return;
		}
	}
	// should never be here
	assert(0);
}
