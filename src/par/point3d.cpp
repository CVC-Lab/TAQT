#include "point3d.h"
#include "mtxlib.h"
#include "vecmath.h"

vector3 Point3D::project(const Point3D& pnt)
{
	vector3 vec;
	float pvec[3];
	for(int i = 0; i < 3; i++) {
		pvec[i] = pnt.pos[i] - pos[i];
	}

	float dot = dotProduct3f(normal, pvec);
	vec.x = pvec[0] - dot * normal[0];
	vec.y = pvec[1] - dot * normal[1];
	vec.z = pvec[2] - dot * normal[2];

	return vec;
}

vector3 Point3D::project(const vector3& vec)
{
	vector3 vtan;
	float dot = normal[0]*vec.x + normal[1]*vec.y + normal[2]*vec.z;
	vtan.x = vec.x - dot * normal[0];
	vtan.y = vec.y - dot * normal[1];
	vtan.z = vec.z - dot * normal[2];

	return vtan;
}

vector2 Point3D::project2d(const Point3D& pnt)
{
	vector2 vec;
	float xaxis[3], yaxis[3], pvec[3];

	localAxes(xaxis, yaxis);

	for(int i = 0; i < 3; i++) {
		pvec[i] = pnt.pos[i] - pos[i];
	}

	vec.x = dotProduct3f(pvec, xaxis);
	vec.y = dotProduct3f(pvec, yaxis);

	return vec;
}

Polar2Coord Point3D::project2dPolar(const Point3D& pnt)
{
	vector2 vec;
	float xaxis[3], yaxis[3], pvec[3];
	
	localAxes(xaxis, yaxis);

	for(int i = 0; i < 3; i++) {
		pvec[i] = pnt.pos[i] - pos[i];
	}
	
	vec.x = dotProduct3f(pvec, xaxis);
	vec.y = dotProduct3f(pvec, yaxis);

	Polar2Coord polar(vec);

	return polar;
}

void Point3D::localAxes(float xaxis[3], float yaxis[3])
{
	if(normal[2] < 0.1) {
		xaxis[0] = -normal[1];
		xaxis[1] = normal[0];
		xaxis[2] = 0;
	} else {
		xaxis[0] = 0;
		xaxis[1] = -normal[2];
		xaxis[2] = normal[1];
	}
	
	crossProduct3f(normal, xaxis, yaxis);
	normalize3f(xaxis);
	normalize3f(yaxis);
}

