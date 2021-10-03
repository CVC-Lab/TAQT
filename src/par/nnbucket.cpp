#include <stdlib.h>
#include <math.h>

#include "nnbucket.h"

NNBucket::NNBucket(int np, Point3D *points, int buck_dim[3] /* = 0 */)
: NNStruct(np, points) 
{
	if(buck_dim == 0) {		// default bucket dimensions
		int n = (int)MIN(33, exp(log((double)np)/3.0)+1);
		m_dim[0] = m_dim[1] = m_dim[2] = n;
	} else {
		m_dim[0] = buck_dim[0];
		m_dim[1] = buck_dim[1];
		m_dim[2] = buck_dim[2];
	}
	init();
}

NNBucket::~NNBucket()
{
	delete[] p_bkidx;
}

void NNBucket::insertAndReplace(int np, Point3D* _points, int pid)
{
	m_np = np;
	p_pnts = _points;
	int bk[3];
	for(int i = 0; i < 3; i++) {
		bk[i] = (int)MIN(m_dim[i]-2, (p_pnts[pid].pos[i] - m_min[i])*(m_dim[i]-1)/(m_max[i]-m_min[i]));
	}
	p_bkidx[index2buck(bk)].insert(pid);
}

int NNBucket::NNSearch(const Point3D& pnt, int k, int nn_idx[], float dist[])
{
	int i, ii, jj, kk, bk[3], m = 0, l = 0;
	
	for(i = 0; i < 3; i++) {
		bk[i] = (int)MIN(m_dim[i]-2, (pnt.pos[i] - m_min[i])*(m_dim[i]-1) / (m_max[i]-m_min[i]));
	}
	
	while (m <= k) {
		m = 0; l++;
		if(l >= m_dim[0]-1) break;
		for(kk = -l; kk <=l; kk++) {
			for(jj = -l; jj <= l; jj++) {
				for(ii = -l; ii <= l; ii++) {
					int bid[3];
					bid[0] = bk[0] + ii;
					bid[1] = bk[1] + jj;
					bid[2] = bk[2] + kk;
					if(valid(bid)) {
						m += p_bkidx[index2buck(bid)].length();
					}
				}
			}
		}
	}
	NNeighbor *neighbors = new NNeighbor[m];

	int n = 0;
	for(kk = -l; kk <=l; kk++) {
		for(jj = -l; jj <= l; jj++) {
			for(ii = -l; ii <= l; ii++) {
				int bid[3];
				bid[0] = bk[0] + ii;
				bid[1] = bk[1] + jj;
				bid[2] = bk[2] + kk;
				if(valid(bid)) {
					int cnt = p_bkidx[index2buck(bid)].length();
					int *idx = p_bkidx[index2buck(bid)].data(); 
					for(i = 0; i < cnt; i++) {
						neighbors[n].id = idx[i];
						neighbors[n].dist = pnt.distance(p_pnts[idx[i]]);
						n++;
					}
				}
			}
		}
	}
	qsort(neighbors, m, sizeof(NNeighbor), compareNNeighbors);

	for(i = 0; i < MIN(m, k); i++) {
		nn_idx[i] = neighbors[i].id;
		dist[i] = neighbors[i].dist;
	}

	delete[] neighbors;
	return MIN(m, k);
}


int NNBucket::NNSearch(const Point3D& pnt, float r2, dynamic_array<int>& nn_idx)
{
	int i, j, k, bk[3], bk_min[3], bk_max[3];
	dynamic_array<int> pid;

	float r = (float) sqrt(r2);
	for(i = 0; i < 3; i++) {
		bk[i] = (int)MIN(m_dim[i]-2, (pnt.pos[i] - m_min[i])*(m_dim[i]-1) / (m_max[i]-m_min[i]));
		bk_min[i] = (int)MAX(0, (pnt.pos[i] - r - m_min[i]) * (m_dim[i]-1) / (m_max[i]-m_min[i]));
		bk_max[i] = (int)MIN(m_dim[i]-2, (pnt.pos[i] + r - m_min[i]) * (m_dim[i]-1) / (m_max[i]-m_min[i]));
	}

	for(k = bk_min[2]; k <= bk_max[2]; k++) {
		for(j = bk_min[1]; j <= bk_max[1]; j++) {
			for(i = bk_min[0]; i <= bk_max[0]; i++) {
				pid.insert(p_bkidx[index2buck(i, j, k)]);
			}
		}
	}

	float dist;
	int m = pid.length();
	NNeighbor *neighbors = new NNeighbor[m];
	nn_idx.clear();

	for(i = 0, j = 0; i < m; i++) {
		dist = pnt.distance2(p_pnts[pid[i]]);
		if(dist < r2) {
			neighbors[j].id = pid[i];
			neighbors[j].dist = dist;
			j++;
		}
	}
	qsort(neighbors, j, sizeof(NNeighbor), compareNNeighbors);
	
	for(i = 0; i < j; i++) {
		nn_idx.insert(neighbors[i].id);
	}
	
	delete[] neighbors;

	return nn_idx.length();
}


int NNBucket::NNSearch(const Point3D& pnt, int n, float r2, int nn_idx[], float dist[]) 
{
	int i, j, k, bk[3], bk_min[3], bk_max[3];
	dynamic_array<int> pid;

	float r = (float)sqrt(r2);
	for(i = 0; i < 3; i++) {
		bk[i] = (int)MIN(m_dim[i]-2, (pnt.pos[i] - m_min[i])*(m_dim[i]-1) / (m_max[i]-m_min[i]));
		bk_min[i] = (int)MAX(0, (pnt.pos[i] - r - m_min[i]) * (m_dim[i]-1) / (m_max[i]-m_min[i]));
		bk_max[i] = (int)MIN(m_dim[i]-2, (pnt.pos[i] + r - m_min[i]) * (m_dim[i]-1) / (m_max[i]-m_min[i]));
	}

	for(k = bk_min[2]; k <= bk_max[2]; k++) {
		for(j = bk_min[1]; j <= bk_max[1]; j++) {
			for(i = bk_min[0]; i <= bk_max[0]; i++) {
				pid.insert(p_bkidx[index2buck(i, j, k)]);
			}
		}
	}
	int m = pid.length();
	NNeighbor *neighbors = new NNeighbor[m];

	for(i = 0, j = 0; i < m; i++) {
		float d2 = pnt.distance2(p_pnts[pid[i]]);
		if(d2 < r2) {
			neighbors[j].id = pid[i];
			neighbors[j].dist = (float)sqrt(d2);
			j++;
		}
	}
	qsort(neighbors, j, sizeof(NNeighbor), compareNNeighbors);

	for(i = 0; i < MIN(j, n); i++) {
		nn_idx[i] = neighbors[i].id;
		dist[i] = neighbors[i].dist;
	}
	delete[] neighbors;

	return MIN(j, n);
}

/************************************************************************/
/*  Private Functions                                                   */
/************************************************************************/
void NNBucket::init()
{
	int i, j;

	for(i = 0; i < 3; i++) {
		m_min[i] = p_pnts[0].pos[i];
		m_max[i] = p_pnts[0].pos[i];
	}
	for(j = 1; j < m_np; j++) {
		for(i = 0; i < 3; i++) {
			if(m_min[i] > p_pnts[j].pos[i]) m_min[i] = p_pnts[j].pos[i];
			if(m_max[i] < p_pnts[j].pos[i]) m_max[i] = p_pnts[j].pos[i];
		}
	}
	for(i = 0; i < 3; i++) {
		m_min[i] -= (m_max[i]-m_min[i]) * 0.01f;
		m_max[i] += (m_max[i]-m_min[i]) * 0.01f;
	}
	p_bkidx = new dynamic_array<int>[(m_dim[0]-1)*(m_dim[1]-1)*(m_dim[2]-1)];
	
	int bk[3];
	for(j = 0; j < m_np; j++) {
		for(i = 0; i < 3; i++) {
			bk[i] = (int)MIN(m_dim[i]-2, (p_pnts[j].pos[i] - m_min[i])*(m_dim[i]-1)/(m_max[i]-m_min[i]));
		}
		p_bkidx[index2buck(bk)].insert(j);
	}
}

