#include "nnkdtree.h"

NNKDTree::NNKDTree(int np, Point3D *points)
: NNStruct(np, points)
{
	m_pntarray = annAllocPts(np, 3);
	for(int i = 0; i < np; i++) {
		m_pntarray[i][0] = points[i].pos[0];
		m_pntarray[i][1] = points[i].pos[1];
		m_pntarray[i][2] = points[i].pos[2];
	}

	p_tree = new ANNkd_tree(m_pntarray, np, 3);
}

NNKDTree::~NNKDTree()
{
	delete p_tree;
	annDeallocPts(m_pntarray);
}

int NNKDTree::NNSearch(const Point3D &pnt, int k, int nn_idx[], float dist[]) 
{
	ANNidxArray 	p_idx = new ANNidx[k];
	ANNdistArray  	p_dist = new ANNdist[k];
	ANNpoint		point = annAllocPt(3);

	point[0] = pnt.pos[0];
	point[1] = pnt.pos[1];
	point[2] = pnt.pos[2];

	p_tree->annkSearch(point, k, p_idx, p_dist);

	for(int i = 0; i < k; i++) {
		nn_idx[i] = p_idx[i];
		dist[i] = p_dist[i];
	}
	delete[] p_dist;
	delete[] p_idx;
	annDeallocPt(point);

	return k;
}

int NNKDTree::NNSearch(const Point3D& pnt, int k, float r2, int nn_idx[], float dist[])
{
	ANNidxArray 	p_idx = new ANNidx[k];
	ANNdistArray  	p_dist = new ANNdist[k];
	ANNpoint		point = annAllocPt(3);
	
	point[0] = pnt.pos[0];
	point[1] = pnt.pos[1];
	point[2] = pnt.pos[2];
	
	p_tree->annkSearch(point, k, p_idx, p_dist);
	
	float r = (float)sqrt(r2);
	int i;
	for(i = 0; i < k; i++) {
		if(p_dist[i] < r) {
			nn_idx[i] = p_idx[i];
			dist[i] = p_dist[i];
		} else 
			break;
	}
	delete[] p_dist;
	delete[] p_idx;
	annDeallocPt(point);
	
	return i;
}

