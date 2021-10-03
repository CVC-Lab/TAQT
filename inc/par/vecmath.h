#ifndef VECMATH_H
#define VECMATH_H

#ifndef PI
#define PI (3.14159265f)
#endif
/**
 * @return length of the vector before normalization	
 */
float normalize3f(float n[3]);
double normalize3d(double vector[3]);
/**
 * dotProduct3f -- compute dot product of two 3d floating-point vectors.
 */
float dotProduct3f(const float p1[3], const float p2[3]);
/**
 * dotProduct3d -- compute dot product of two 3d double precision vectors.                                                  
 */
 double dotProduct3d(const double u[3], const double v[3]);

/**
 * normalProduct3d -- compute normalized vector W as the cross product of U and V.
 * @note The direction of W is chose so that W.tan > 0.
 */
void normalProduct3d(double U[3], double V[3], double tan[3], double W[3]);

/**
 * crossProductNormalized3d -- compute normalized cross product W of two vectors U and V.                                                    
 */
void crossProductNormalized3d(double U[3], double V[3], double W[3]);

void crossProduct3f(const float u[3], const float v[3], float w[3]);

/** 
 * projectVectorToPlane3f -- Project a vector to palne defined by 
 *                (p-point)^T normal = 0                                     
 * @param vector: the given vector to be projected
 * @parm  normal: the plane normal             
 * @param provector: the result of the projection                                  
 */
void projectVectorToPlane3f(float vector[3], float normal[3], float point[3], float provector[3]);

/**
 * determinant3f -- compute the determinant of three vectors [v1,v2,v3] in R^3
 */
float determinant3f(float v1[3], float v2[3], float v3[3]);

 /**
 * determinant3d-- compute the determinant of three vectors [v1,v2,v3] in R^3 
 */
double  determinant3d(double v1[3], double v2[3], double v3[3]);

/**
 * DDeterminant_Matrix -- compute the determinant of  a 3x3 matrix M 
 */
double  determinantMatrix3d(double M[3]);

/**
 * RotateMatrix -- Compute a rotate matrix M that rotate  vector n1 to n2
 *                 along the vector n3, that is, M * n1 = n2             
 *                 where n1, n2, n3  are normalized vectors and n3 is    
 *                 perpendicular to n1 and n2.                           
 */
void RotateMatrix(double n1[3], double n2[3], double n3[3], double matrix[3][3]);

/**
 * RotateMatrix_z -- Rotate the vector (nx,ny,nz)^T to be z-axis      
 *                 That is, matrix^T*(nx,ny,nz)^T = [0,0, |n|]          
 */
void RotateMatrix_z(double nx, double ny, double nz, double matrix[3][3]);

/**
 * TriangleNormal3f -- compute the normal of a triangle by the right-handed rule.
 * @param pi: the vertices of the triangle
 */
void TriangleNormal(float p1[3], float p2[3], float p3[3], float normal[3]);

/**
 * pointDistance2f -- Distance of two 2D points
 */
float pointDistance2f(float p1[2], float p2[2]);

/**
 * pointDistance3f -- Distance of two 3D points
 */
float pointDistance3f(float p1[3], float p2[3]);

#endif

