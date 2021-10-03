
#include <stdio.h>
#include <math.h>
#include "vecmath.h"

float normalize3f(float n[3])
{
	float  result;
	result = (float) sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
	if (result == 0.0) {
		/* default vector */
		n[0] = n[1] = 0; n[2] = 1;
		return 0;
	}
	n[0] = n[0]/result;
	n[1] = n[1]/result;
	n[2] = n[2]/result;

	return result;
}

double normalize3d(double vector[3])
{
	double  length;

	length = sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
	if (length == 0.0) {
		vector[0] = vector[1] = 0; vector[2] = 1;
		return length;
	}
	vector[0] = vector[0]/length;
	vector[1] = vector[1]/length;
	vector[2] = vector[2]/length;

	return length;
}

/**
 * dotProduct3f -- compute dot product of two 3d floating-point vectors.
 */
float dotProduct3f(const float p1[3], const float p2[3])
{
	return (p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2]);
}


/**
 * dotProduct3d -- compute dot product of two 3d double precision vectors.                                                  
 */
double dotProduct3d(const double u[3], const double v[3])
{
	double w;
	w = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
	return(w);
}

/**
 * innerAngle3f -- compute cos of the inner angle between two vectors.
 */
float innerAngle3f(float p1[3], float p2[3])
{
	float   result;

	result = p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2];
	if(result == 0) return 0;

    /* normalization */
	result = (float) (result/sqrt((p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2])) ); 
	result = (float) (result/sqrt((p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2])) ); 

	return(result);
}
/**
 * normalProduct3d -- compute normalized vector W as the cross product of U and V.
 * @note The direction of W is chose so that W.tan > 0.
 */
void normalProduct3d(double U[3], double V[3], double tan[3], double W[3]) 
{
	double   len, InnerPro;
	/* cross product               */
	W[0] = U[1]*V[2] - U[2]*V[1];
	W[1] = U[2]*V[0] - U[0]*V[2];
	W[2] = U[0]*V[1] - U[1]*V[0];
	/* normal                      */
	len = sqrt(W[0]* W[0] + W[1]* W[1] + W[2]* W[2]);
	if(len == 0) return;

	/* inner product of W and tan  */
	InnerPro =  W[0]* tan[0] + W[1]* tan[1] + W[2]* tan[2];
	if (InnerPro < 0.0) len = - len;

	W[0] = W[0]/len;
	W[1] = W[1]/len;
	W[2] = W[2]/len;
}


/**
 * crossProductNormalized3d -- compute normalized cross product W to two vectors U and V.                                                    
 */
void crossProductNormalized3d(double U[3], double V[3], double W[3]) 
{
	double   len;
	/* cross product               */
	W[0] = U[1]*V[2] - U[2]*V[1];
	W[1] = U[2]*V[0] - U[0]*V[2];
	W[2] = U[0]*V[1] - U[1]*V[0];
	/* normal                      */
	len = sqrt(W[0]* W[0] + W[1]* W[1] + W[2]* W[2]);
	if(len == 0) return;

	W[0] = W[0]/len;
	W[1] = W[1]/len;
	W[2] = W[2]/len;
}


void crossProduct3f(const float u[3], const float v[3], float w[3])
{
	w[0] = u[1]*v[2] - u[2]*v[1];
	w[1] = u[2]*v[0] - v[2]*u[0];
	w[2] = u[0]*v[1] - v[0]*u[1];
}

/** 
 * projectVectorToPlane3f -- Project a vector to palne defined by 
 *                (p-point)^T normal = 0                                     
 * @param vector: the given vector to be projected
 * @parm  normal: the plane normal             
 * @param provector: the result of the projection                                  
 */
void projectVectorToPlane3f(float vector[3], float normal[3], float point[3], float provector[3])
{
	double matrix[3][3],nx,ny,nz;
	float  b,c;
	int    i;

	nx = normal[0]; 
	ny = normal[1];
	nz = normal[2];
	RotateMatrix_z(nx,ny,nz,matrix);
	b = 0.0f;  
	c = 0.0f;
	for (i = 0; i < 3; i++) {
		b = float (b + vector[i]*matrix[i][0]);
		c = float (c + vector[i]*matrix[i][1]);
	}

	for (i = 0; i < 3; i++) {
		provector[i] = float (b*matrix[i][0] + c*matrix[i][1]);
	}
}

/**
 * determinant3f -- compute the determinant of three vectors [v1,v2,v3] in R^3
 */
float determinant3f(float v1[3], float v2[3], float v3[3])
{
	float result;

	result = v1[0]*v2[1]*v3[2] + v1[1]*v2[2]*v3[0] + v1[2]*v2[0]*v3[1] - 
		v1[2]*v2[1]*v3[0] - v1[1]*v2[0]*v3[2] - v1[0]*v2[2]*v3[1];
	return(result);
} 

/**
 * determinant3d-- compute the determinant of three vectors [v1,v2,v3] in R^3 
 */
double  determinant3d(double v1[3], double v2[3], double v3[3])
{
	double result;

	result = v1[0]*v2[1]*v3[2] + v1[1]*v2[2]*v3[0] + v1[2]*v2[0]*v3[1] -
		v1[2]*v2[1]*v3[0] - v1[1]*v2[0]*v3[2] - v1[0]*v2[2]*v3[1];
	return(result);
}

/**
 * DDeterminant_Matrix -- compute the determinant of  a 3x3 matrix M 
 */
double  determinantMatrix3d(double M[3])
{
	double result;

	result = M[0]*M[4]*M[8] + M[3]*M[7]*M[2] + M[6]*M[1]*M[5] -
		M[6]*M[4]*M[2] - M[3]*M[1]*M[8] - M[0]*M[7]*M[5];
	return(result);
}

/**
 * RotateMatrix -- Compute a rotate matrix M that rotate  vector n1 to n2
 *                 along the vector n3, that is, M * n1 = n2             
 *                 where n1, n2, n3  are normalized vectors and n3 is    
 *                 perpendicular to n1 and n2.                           
 */
void RotateMatrix(double n1[3], double n2[3], double n3[3], double matrix[3][3])
{
	int     i;
	double  matrix1[3][3], matrix2[3][3],n13[3],n23[3];

	crossProductNormalized3d(n1,n3,n13);
	crossProductNormalized3d(n2,n3,n23);

	matrix1[0][0] = n2[0];
	matrix1[1][0] = n2[1];
	matrix1[2][0] = n2[2];

	matrix1[0][1] = n23[0];
	matrix1[1][1] = n23[1];
	matrix1[2][1] = n23[2];

	matrix1[0][2] = n3[0];
	matrix1[1][2] = n3[1];
	matrix1[2][2] = n3[2];

	matrix2[0][0] = n1[0];
	matrix2[0][1] = n1[1];
	matrix2[0][2] = n1[2];

	matrix2[1][0] = n13[0];
	matrix2[1][1] = n13[1];
	matrix2[1][2] = n13[2];

	matrix2[2][0] = n3[0];
	matrix2[2][1] = n3[1];
	matrix2[2][2] = n3[2];

	for (i = 0; i < 3; i++) {
		matrix[i][0] = matrix1[i][0]*matrix2[0][0]
		+ matrix1[i][1]*matrix2[1][0]
		+ matrix1[i][2]*matrix2[2][0];

		matrix[i][1] = matrix1[i][0]*matrix2[0][1]
		+ matrix1[i][1]*matrix2[1][1]
		+ matrix1[i][2]*matrix2[2][1];

		matrix[i][2] = matrix1[i][0]*matrix2[0][2]
		+ matrix1[i][1]*matrix2[1][2]
		+ matrix1[i][2]*matrix2[2][2];
	}
}

/**
 * RotateMatrix_z -- Rotate the vector (nx,ny,nz)^T to be z-axis      
 *                 That is, matrix^T*(nx,ny,nz)^T = [0,0, |n|]          
 */
void RotateMatrix_z(double nx, double ny, double nz, double matrix[3][3])
{
	double  c1,c2,s1,s2, normal,normalz;

	normal = sqrt(nx*nx + ny*ny + nz*nz);
	normalz = sqrt(nx*nx + ny*ny);
	c1 = nz/normal;       c2 = -ny/normalz;
	s1 = normalz/normal;  s2 = nx/normalz;

	if (normalz < 0.001) {
		matrix[0][0] = 1.0;    matrix[0][1] = 0.0;    matrix[0][2] = 0.0;
		matrix[1][0] = 0.0;    matrix[1][1] = 1.0;    matrix[1][2] = 0.0;
		matrix[2][0] = 0.0;    matrix[2][1] = 0.0;    matrix[2][2] = 1.0;
	}
	if (normalz >= 0.001) {
		matrix[0][0] = c2;     matrix[0][1] = -c1*s2; matrix[0][2] = s1*s2;
		matrix[1][0] = s2;     matrix[1][1] = c1*c2;  matrix[1][2] = -s1*c2;
		matrix[2][0] = 0.0;    matrix[2][1] = s1;     matrix[2][2] = c1;
	}
}

/**
 * TriangleNormal3f -- compute the normal of a triangle by the right-handed rule.
 * @param pi: the vertices of the triangle
 */
void TriangleNormal(float p1[3], float p2[3], float p3[3], float normal[3])
{
	float       x1,y1,z1,x2,y2,z2;

	x1 = p1[0] - p3[0];
	y1 = p1[1] - p3[1];
	z1 = p1[2] - p3[2];

	x2 = p2[0] - p3[0];
	y2 = p2[1] - p3[1];
	z2 = p2[2] - p3[2];

	normal[0] = y1*z2 - y2*z1;
	normal[1] = x2*z1 - x1*z2;
	normal[2] = x1*y2 - x2*y1;
	normalize3f(normal);
}

/**
 * Barycentric_Coord_3D_Triangle -- compute the barycentric coordinate of a 
 * point in term a 3D triangle in space. 
 */
void Barycentric_Coord_3D_Triangle(float p[3], float p1[3], float p2[3], float p3[3], float b123[3])
{
/*	int   i,info; 
	float pp3[3], pp13[3], pp23[3];
	float c1,c2,c3,c4,y1,y2; 

	for (i = 0; i < 3; i++) {
		pp3[i]  = p[i]  - p3[i]; 
		pp13[i] = p1[i] - p3[i]; 
		pp23[i] = p2[i] - p3[i]; 
	}
	c1 = innerAngle3f(pp13,pp13); 
	c2 = innerAngle3f(pp23,pp13); 
	c3 = c2; 
	c4 = innerAngle3f(pp23,pp23);
	y1 = innerAngle3f(pp3,pp13);
	y2 = innerAngle3f(pp3,pp23);
	if (c1 > c4) {
		c2 = c2/c1; 
		c3 = c3/c1;
		c4 = c4/c1;  
		y1 = y1/c1; 
		y2 = y2/c1; 
		c1 = 1.0; 
	} else {
		c1 = c1/c4;  
		c2 = c2/c4; 
		c3 = c3/c4;
		y1 = y1/c4; 
		y2 = y2/c4; 
		c4 = 1.0; 
	}

	info = LinearSystemOrder2(c1,c2,c3,c4,y1,y2,b123,b123+1); 
	if (info == 0) {
		printf("The triangle is degenerate, in Barycentric_Coord_3D_Triangle\n"); 
		b123[0] = 0.3333333; 
		b123[1] = 0.3333333; 
		b123[2] = 0.3333333; 
		printf("c1,c2,c3,c4,y1,y2 = %f,%f, %f, %f,%f, %f\n", c1,c2,c3,c4,y1,y2); 
		printf("b123 = %f,%f, %f\n", b123[0],b123[1],b123[2]); 
		printf("p = %f,%f,%f\n", p[0],p[1],p[2]); 
		return; 
	}
	b123[2] = 1.0 - b123[0] - b123[1]; */
}


/**
 * pointDistance2f -- Distance of two 2D points
 */
float pointDistance2f(float p1[2], float p2[2])
{
	float result; 

	result = (float) sqrt((p1[0] - p2[0])*(p1[0] - p2[0]) + 
		(p1[1] - p2[1])*(p1[1] - p2[1])); 

	return(result); 
}

/**
 * pointDistance3f -- Distance of two 3D points
 */
float pointDistance3f(float p1[3], float p2[3])
{
	float result;

	result = (float) sqrt((p1[0] - p2[0])*(p1[0] - p2[0]) +
		(p1[1] - p2[1])*(p1[1] - p2[1]) +
		(p1[2] - p2[2])*(p1[2] - p2[2]) );

	return(result);
}





