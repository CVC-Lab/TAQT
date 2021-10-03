#ifndef _POINT3D_H
#define _POINT3D_H

#include <math.h>

#include "mtxlib.h"
#include "vecmath.h"

/*
 *	Polar 2D coordinates
 */

class Polar2Coord
{
public:
	Polar2Coord(int _id = 0) : id(_id) {
		r = 1;
		theta = 0;
	}

	// Construct polar coordinates from Eculiean coordinates
	Polar2Coord(float x[2], int _id = 0) : id(_id) {
		r = (float) sqrt(x[0]*x[0] + x[1]*x[1]);
		if(r == 0) {
			theta = 0;
		} else {
			theta = (float) acos(x[0]/r);
			if(x[1] < 0) theta = 2*PI-theta;
		}
	}

	Polar2Coord(const vector2 vec) {
		r = (float)sqrt(vec.x*vec.x + vec.y*vec.y);
		if(r == 0) {
			theta = 0;
		} else {
			theta = (float) acos(vec.x/r);
			if(vec.y < 0) theta = 2*PI-theta;
		}
	}

	bool operator < (const Polar2Coord& p1) const {
		if(theta < p1.theta) return true;
		return false;
	}

	int id;
	float r;
	float theta;		// 0 <= theta < 2*PI
};

/**
 * A point in 3D space.
 */
class Point3D {
public:
	Point3D() {
		pos[0] = pos[1] = pos[2] = 0;
		normal[0] = normal[1] = 0; normal[2] = 1;
		color = 0;
	}

	Point3D(float _pos[3]) {
		for(int i = 0; i < 3; i++) {
			pos[i] = _pos[i];
		}
		normal[0] = normal[1] = 0; normal[2] = 1;
		color = 0;
	}

	Point3D(float _pos[3], float _norm[3], int  _clr = 0) {
		for(int i = 0; i < 3; i++) {
			pos[i] = _pos[i];
			normal[i] = _norm[i];
			color = _clr;
		}
	}

	~Point3D() {}

	void setPos(float _pos[3]) {
		pos[0] = _pos[0];
		pos[1] = _pos[1];
		pos[2] = _pos[2];
	}

	void setNormal(float _norm[3]) {
		normal[0] = _norm[0];
		normal[1] = _norm[1];
		normal[2] = _norm[2];
	}

	void flipNormal() {
		for(int i = 0; i < 3; i++) 
			normal[i] = -normal[i];
	}

	// r g b values have 8 bytes
	void setColor(int r, int g, int b) {
		color = (r << 16) | (g << 8) | b; 
	}

	void getColor(int &r, int &g, int &b) {
		b = color & 0x11111111;
		g = (color >> 8) & 0x11111111;
		r = (color >> 16) & 0x11111111;
	}

	/**
	 * project a point to the tangent plane passing through a point.
	 * @return A 3D vector from this point to the projected point
	 * @note The normal at this point should have been normalized.
	 */
	vector3 project(const Point3D& pnt);

	/**
	 * Project a vector to the tangent plane of a point.
	 */
	vector3 project(const vector3& vec);

	/**
	 * Project a point into a 2D vector in the tangent plane.
	 */
	vector2 project2d(const Point3D& pnt);

	/*
	 *	Project a point into 2D with polar coordinates
	 */
	Polar2Coord project2dPolar(const Point3D& pnt);

	/*
	 *	X and Y axes of local tangent plane
	 */
	void localAxes(float xaxis[3], float yaxis[3]);

	/*
	 *	Compute the distance between two points.
	 */
	float distance(const Point3D& pnt) const {
		return (float) sqrt((pos[0]-pnt.pos[0])*(pos[0]-pnt.pos[0]) +
							(pos[1]-pnt.pos[1])*(pos[1]-pnt.pos[1]) +
							(pos[2]-pnt.pos[2])*(pos[2]-pnt.pos[2]));
	}

	/*
	 *	Compute the square distance between two points.
	 */
	float distance2(const Point3D& pnt) const {
		return	(pos[0]-pnt.pos[0])*(pos[0]-pnt.pos[0]) +
				(pos[1]-pnt.pos[1])*(pos[1]-pnt.pos[1]) +
				(pos[2]-pnt.pos[2])*(pos[2]-pnt.pos[2]);
	}
	
	float	pos[3];			// position
	float	normal[3];		// normal
	int		color;			// rgb color or index to color map	
};
#endif

