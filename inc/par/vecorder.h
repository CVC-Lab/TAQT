#ifndef _VECTOR_ORDER_H
#define _VECTOR_ORDER_H

#include "mtxlib.h"

/**
 * Sort vectors in a 2-D plane couter-clockwise.
 * @param normal is the normal vector of the plane.
 */
void sortVectorCC(vector3 vecs[], int nv, const vector3 normal);

/**
 * The plane angle between two vectors.
 */
class PlaneAngle
{
public:
	// The angle is pointed from v0 to v1.
	PlaneAngle(float ctan, int sgny) : m_ctan(ctan), m_sgny(sgny) {
	};

	friend bool operator < (const PlaneAngle& a1, const PlaneAngle& a2);

protected:
	float	m_ctan;			// value of ctan(theta)
	int		m_sgny;			// > 0 if theta < PI and < 0 if theta > PI
};

#endif


