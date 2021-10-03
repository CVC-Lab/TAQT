#include "surface3d.h"

Surface3D::Surface3D(float _lower[3] /* = NULL */, float _upper[3] /* = NULL */)
: m_colored(false)
{
	int i;
	for(i = 0; i < 3; i++) {
		m_minext[i] = 0;
		m_maxext[i] = 1;
	}

	if(_lower && _upper) {
		for(i = 0; i < 3; i++) {
			m_minext[i] = _lower[i];
			m_maxext[i] = _upper[i];
		}
	}
}

Surface3D::~Surface3D()
{
}

int Surface3D::addVert(float pos[3], float norm[3], int c)
{
	for(int i = 0; i < 3; i++) {
		if(pos[i] < m_minext[i]) m_minext[i] = pos[i];
		if(pos[i] > m_maxext[i]) m_maxext[i] = pos[i];
	} 
	Point3D point(pos, norm, c);
	return	verts.insert(point);
}

int Surface3D::addTri(int v1, int v2, int v3)
{
	TriIndex tri(v1, v2, v3);
	return tris.insert(tri);
}

void Surface3D::reset()
{
	verts.clear();
	tris.clear();
}

