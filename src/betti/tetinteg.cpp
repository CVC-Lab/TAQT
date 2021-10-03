#include <assert.h>
#include "volume3dcritical.h"

//#define DEBUG_VOLVOL 1
float tetraFuncIntegral(float x1[3], float x2[3], float x3[3], float x4[3], 
						float v1, float v2, float v3, float v4,
						float f1, float f2, float f3, float f4,
						float fx)
{
	float *_t, t;
	
	 float f[4], v[4];
	/* v[0] = v1; v[1] = v2; v[2] = v3; v[3] = v4;
	f[0] = f1; f[1] = f2; f[2] = f3; f[3] = f4;
	degenerate(x1, x2, x3, x4, v, f);
	v1 = v[0]; v2 = v[1]; v3 = v[2]; v4 = v[3];
	f1 = f[0]; f2 = f[1]; f3 = f[2]; f4 = f[3]; 
	assert(v1 < v2 && v2 < v3 && v3 < v4); */
	if (fx <= v1) return 0;

	float vec1[3], vec2[3], vec3[3], cp[3], mid[3], mid1[3], mid2[3];
	float volume, area1, area2, midarea, ival;

	vec1[0] = x2[0]-x1[0];
	vec1[1] = x2[1]-x1[1];
	vec1[2] = x2[2]-x1[2];
	vec2[0] = x3[0]-x1[0];
	vec2[1] = x3[1]-x1[1];
	vec2[2] = x3[2]-x1[2];
	vec3[0] = x4[0]-x1[0];
	vec3[1] = x4[1]-x1[1];
	vec3[2] = x4[2]-x1[2];
	cp[0] = vec3[0]*(vec1[1]*vec2[2]-vec1[2]*vec2[1]);
	cp[1] = vec3[1]*(vec1[2]*vec2[0]-vec1[0]*vec2[2]);
	cp[2] = vec3[2]*(vec1[0]*vec2[1]-vec1[1]*vec2[0]);
	//volume = sqrt(cp[0]*cp[0] + cp[1]*cp[1] + cp[2]*cp[2])/6.0;
	volume = fabs(cp[0] + cp[1] + cp[2])/6.0;

	if (fx >= v4) return volume*(f1+f2+f3+f4)/4;

	// a cell with constant value will not contribute to the volume
#ifdef DEBUG_VOLVOL
	printf("v1 = %.10f, v2 = %.10f, v3 = %.10f, v4 = %.10f, epsilon = %f\n", v1, v2, v3, v4, epsilon);
	printf(">>> v1 = %.10f, v2 = %.10f, v3 = %.10f, v4 = %.10f\n", v1, v2, v3, v4);
#endif

	// compute the first area
	if (v1 != v3)
		ival = (v3-v2)/(v3-v1);
	else
		ival = 0.0;
	mid[0] = (1.0-ival)*x3[0] + (ival)*x1[0];
	mid[1] = (1.0-ival)*x3[1] + (ival)*x1[1];
	mid[2] = (1.0-ival)*x3[2] + (ival)*x1[2];
	if (v1 != v4)
		ival = (v4-v2)/(v4-v1);
	else
		ival = 0.0;
	mid2[0] = (1.0-ival)*x4[0] + (ival)*x1[0];
	mid2[1] = (1.0-ival)*x4[1] + (ival)*x1[1];
	mid2[2] = (1.0-ival)*x4[2] + (ival)*x1[2];

	vec1[0] = mid[0]-x2[0];
	vec1[1] = mid[1]-x2[1];
	vec1[2] = mid[2]-x2[2];
	vec2[0] = mid2[0]-x2[0];
	vec2[1] = mid2[1]-x2[1];
	vec2[2] = mid2[2]-x2[2];
	cp[0] = vec1[1]*vec2[2]-vec1[2]*vec2[1];
	cp[1] = vec1[2]*vec2[0]-vec1[0]*vec2[2];
	cp[2] = vec1[0]*vec2[1]-vec1[1]*vec2[0];
	area1 = 0.5 * fabs(sqrt(cp[0]*cp[0] + cp[1]*cp[1] + cp[2]*cp[2]));


	// compute the second area
	if (v2 != v4)
		ival = (v4-v3)/(v4-v2);
	else
		ival = 0.0;
	mid[0] = (1.0-ival)*x4[0] + (ival)*x2[0];
	mid[1] = (1.0-ival)*x4[1] + (ival)*x2[1];
	mid[2] = (1.0-ival)*x4[2] + (ival)*x2[2];
	if (v4 != v1)
		ival = (v4-v3)/(v4-v1);
	else
		ival = 0.0;
	mid2[0] = (1.0-ival)*x4[0] + (ival)*x1[0];
	mid2[1] = (1.0-ival)*x4[1] + (ival)*x1[1];
	mid2[2] = (1.0-ival)*x4[2] + (ival)*x1[2];

	vec1[0] = mid[0]-x3[0];
	vec1[1] = mid[1]-x3[1];
	vec1[2] = mid[2]-x3[2];
	vec2[0] = mid2[0]-x3[0];
	vec2[1] = mid2[1]-x3[1];
	vec2[2] = mid2[2]-x3[2];
	cp[0] = vec1[1]*vec2[2]-vec1[2]*vec2[1];
	cp[1] = vec1[2]*vec2[0]-vec1[0]*vec2[2];
	cp[2] = vec1[0]*vec2[1]-vec1[1]*vec2[0];
	area2 = 0.5 * fabs(sqrt(cp[0]*cp[0] + cp[1]*cp[1] + cp[2]*cp[2]));

	if (v2-v1 >= v4-v3)
		midarea = area1 * (1.0+(v3-v2)/(v2-v1));
	else if (v4-v3 > v2-v1)
		midarea = area2 * (1.0+(v3-v2)/(v4-v3));
	else {
		assert(0);
		// have to compute the midarea
		// This formula is wrong -- xiaoyu
		vec1[0] = (x2[0]-x1[0])/2;
		vec1[1] = (x2[1]-x1[1])/2;
		vec1[2] = (x2[2]-x1[2])/2;
		vec2[0] = (x4[0]-x3[0])/2;
		vec2[1] = (x4[1]-x3[1])/2;
		vec2[2] = (x4[2]-x3[2])/2;
		cp[0] = (vec1[1]*vec2[2]-vec1[2]*vec2[1]);
		cp[1] = (vec1[2]*vec2[0]-vec1[0]*vec2[2]);
		cp[2] = (vec1[0]*vec2[1]-vec1[1]*vec2[0]);
		midarea = 2*(sqrt(cp[0]*cp[0] + cp[1]*cp[1] + cp[2]*cp[2])) - (area1+area2)/2;
	}
#ifdef DEBUG_VOLVOL
	printf("area1: %f area2: %f mid: %f\n", area1, area2, midarea);
	printf("v1: %f v2: %f v3: %f v4: %f\n", v1, v2, v3, v4);
#endif

	float vol2 = ((v3-v1)*area1 + (v3-v2)*midarea + (v4-v2)*area2)/3;
	float factor = volume / vol2;

	float s, vol, fint, func[4];
	func[0] = f1;
	func[1] = f2;
	func[2] = (v2-v1)*f4/(v4-v1) + (v4-v2)*f1/(v4-v1);
	func[3] = (v2-v1)*f3/(v3-v1) + (v3-v2)*f1/(v3-v1);
	//printf("v: %f %f %f %f, f: %f %f %f %f\n", v1, v2, v3, v4, f1, f2, f3, f4);
	if (fx <= v2) {
		s = (fx - v1) / (v2 - v1);
		vol = factor*area1*(v2-v1)*s*s*s / 3;
		//printf("vol: %f, func: %f %f %f %f\n", vol, func[0], func[1], func[2], func[3]);
		return vol * ((4-3*s)*func[0] + s*func[1] + s*func[2] + s*func[3])/4;
	}
	float cum1 = factor*area1*(v2-v1)*(func[0]+func[1]+func[2]+func[3])/12;

	float A = f3-f2 + (v3-v2)*(f4-f1)/(v4-v1) + (v3-v2)*(f3-f1)/(v3-v1);
	float B = f2 + (v2-v1)*f4/(v4-v1) + (v4-v2)*f1/(v4-v1) + (v2-v1)*f3/(v3-v1) + (v3-v2)*f1/(v3-v1);
	float r = (v3-v2)/(v2-v1);
	float C = f3-f2 + (v3-v2)*(f4-f1)/(v4-v1) + (v3-v2)*(f4-f2)/(v4-v2);
	float D = f2 + (v2-v1)*f4/(v4-v1) + (v4-v2)*f1/(v4-v1) + f2;
	t = (v3-v2)/(v4-v3);        
	if (fx <= v3) {
		s = (fx-v2)/(v3-v2);
		float s2 = s*s;
		float s3 = s2*s;
		float s4 = s3*s;
		//float fint1 = factor*area1*(v3-v2)*(B*s+ (A+B*(r-1))*s2/2 + ((r-1)*A-r*B)*s3/3 - A*r*s4/4)/3;
		float fint1 = factor*(v3-v2)*s*(A*s*(area1*(6-8*s+3*s2)+(4-3*s)*s*midarea)+2*B*(2*area1*(3-3*s+s2)+(3-2*s)*s*midarea))/36;
		//float fint2 = factor*area2*(v3-v2)*(D*(1+t)*s2/2+(C*(1+t)-D*t)*s3/3-C*t*s4/4)/3;
		float fint2 = factor*(v3-v2)*s2*(C*s*(3*s*(area2-midarea)+4*midarea)+D*(4*s*(area2-midarea)+6*midarea))/36;
#ifdef DEBUG_VOLVOL
		printf("A = %f, B = %f, C = %f, D = %f, r = %f, t = %f, s = %f\n", A, B, C, D, r, t,s );
		printf("fint1 = %.10f, fint2 = %.10f, cum1 = %f\n", fint1, fint2, cum1+fint1+fint2);
#endif
		return cum1+fint1+fint2;
	}
	//cum1 += ((B/2 + B*r/6 + A/6 + A*r/12)*area1 + (D/2 + D*t/6 + C/3 + C*t/12)*area2)*factor*(v3-v2)/3;
	cum1 += factor*(v3-v2)*((A*(area1+midarea)+2*B*(2*area1+midarea))/36 + (2*D*(2*area2+midarea)+C*(3*area2+midarea))/36);

	// v3 < fx < v4
	func[0] = f4;
	func[1] = f3;
	func[2] = (v4-v3)*f1/(v4-v1) + (v3-v1)*f4/(v4-v1);
	func[3] = (v4-v3)*f2/(v4-v2) + (v3-v2)*f4/(v4-v2);
	s = (fx-v3)/(v4-v3);
	fint = factor*area2*(v4-v3)*((func[0]+func[1]+func[2]+func[3]) - 
								 (1-s)*(1-s)*(1-s)*((1+3*s)*func[0]+(1-s)*func[1]+(1-s)*func[2]+(1-s)*func[3]))/12;

#ifdef DEBUG_VOLVOL
	printf("fx = %f, cum1 = %f\n", fx, cum1+fint);
#endif
	return cum1+fint;
}

void tetIntegMultiFuncs(float x1[3], float x2[3], float x3[3], float x4[3],
						float v1, float v2, float v3, float v4,
						float *f1, float *f2, float *f3, float *f4, float *out, int n, float fx)
{
	int i;
	
#ifdef _DEBUG
	assert(v1 < v2 && v2 < v3 && v3 < v4);
#endif

	for(i = 0; i < n; i++) out[i] = 0;
	
	if (fx <= v1) return;

	float vec1[3], vec2[3], vec3[3], cp[3], mid[3], mid1[3], mid2[3];
	float volume, area1, area2, midarea, ival;

	vec1[0] = x2[0]-x1[0];
	vec1[1] = x2[1]-x1[1];
	vec1[2] = x2[2]-x1[2];
	vec2[0] = x3[0]-x1[0];
	vec2[1] = x3[1]-x1[1];
	vec2[2] = x3[2]-x1[2];
	vec3[0] = x4[0]-x1[0];
	vec3[1] = x4[1]-x1[1];
	vec3[2] = x4[2]-x1[2];
	cp[0] = vec3[0]*(vec1[1]*vec2[2]-vec1[2]*vec2[1]);
	cp[1] = vec3[1]*(vec1[2]*vec2[0]-vec1[0]*vec2[2]);
	cp[2] = vec3[2]*(vec1[0]*vec2[1]-vec1[1]*vec2[0]);
	volume = fabs(cp[0] + cp[1] + cp[2])/6.0;

	if (fx >= v4) {
		for(i = 0; i < n; i++) {
			out[i] = volume * (f1[i]+f2[i]+f3[i]+f4[i])/4;
		}
		return;
	}

	// a cell with constant value will not contribute to the volume
#ifdef DEBUG_VOLVOL
	printf("v1 = %.10f, v2 = %.10f, v3 = %.10f, v4 = %.10f, epsilon = %f\n", v1, v2, v3, v4, epsilon);
	printf(">>> v1 = %.10f, v2 = %.10f, v3 = %.10f, v4 = %.10f\n", v1, v2, v3, v4);
#endif

	// compute the first area
	if (v1 != v3)
		ival = (v3-v2)/(v3-v1);
	else
		ival = 0.0;
	mid[0] = (1.0-ival)*x3[0] + (ival)*x1[0];
	mid[1] = (1.0-ival)*x3[1] + (ival)*x1[1];
	mid[2] = (1.0-ival)*x3[2] + (ival)*x1[2];
	if (v1 != v4)
		ival = (v4-v2)/(v4-v1);
	else
		ival = 0.0;
	mid2[0] = (1.0-ival)*x4[0] + (ival)*x1[0];
	mid2[1] = (1.0-ival)*x4[1] + (ival)*x1[1];
	mid2[2] = (1.0-ival)*x4[2] + (ival)*x1[2];

	vec1[0] = mid[0]-x2[0];
	vec1[1] = mid[1]-x2[1];
	vec1[2] = mid[2]-x2[2];
	vec2[0] = mid2[0]-x2[0];
	vec2[1] = mid2[1]-x2[1];
	vec2[2] = mid2[2]-x2[2];
	cp[0] = vec1[1]*vec2[2]-vec1[2]*vec2[1];
	cp[1] = vec1[2]*vec2[0]-vec1[0]*vec2[2];
	cp[2] = vec1[0]*vec2[1]-vec1[1]*vec2[0];
	area1 = 0.5 * fabs(sqrt(cp[0]*cp[0] + cp[1]*cp[1] + cp[2]*cp[2]));


	// compute the second area
	if (v2 != v4) ival = (v4-v3)/(v4-v2);
	else ival = 0.0;
	mid[0] = (1.0-ival)*x4[0] + (ival)*x2[0];
	mid[1] = (1.0-ival)*x4[1] + (ival)*x2[1];
	mid[2] = (1.0-ival)*x4[2] + (ival)*x2[2];
	if (v4 != v1) ival = (v4-v3)/(v4-v1);
	else ival = 0.0;
	mid2[0] = (1.0-ival)*x4[0] + (ival)*x1[0];
	mid2[1] = (1.0-ival)*x4[1] + (ival)*x1[1];
	mid2[2] = (1.0-ival)*x4[2] + (ival)*x1[2];

	vec1[0] = mid[0]-x3[0];
	vec1[1] = mid[1]-x3[1];
	vec1[2] = mid[2]-x3[2];
	vec2[0] = mid2[0]-x3[0];
	vec2[1] = mid2[1]-x3[1];
	vec2[2] = mid2[2]-x3[2];
	cp[0] = vec1[1]*vec2[2]-vec1[2]*vec2[1];
	cp[1] = vec1[2]*vec2[0]-vec1[0]*vec2[2];
	cp[2] = vec1[0]*vec2[1]-vec1[1]*vec2[0];
	area2 = 0.5 * fabs(sqrt(cp[0]*cp[0] + cp[1]*cp[1] + cp[2]*cp[2]));

	if (v2-v1 >= v4-v3)
		midarea = area1 * (1.0+(v3-v2)/(v2-v1));
	else 
		midarea = area2 * (1.0+(v3-v2)/(v4-v3));
	
#ifdef DEBUG_VOLVOL
	printf("area1: %f area2: %f mid: %f\n", area1, area2, midarea);
	printf("v1: %f v2: %f v3: %f v4: %f\n", v1, v2, v3, v4);
#endif

	float vol2 = ((v3-v1)*area1 + (v3-v2)*midarea + (v4-v2)*area2)/3;
	float factor = volume / vol2;

	float s, vol, fint, func[4];
	float *cum1 = new float[n];
	
	//printf("v: %f %f %f %f, f: %f %f %f %f\n", v1, v2, v3, v4, f1, f2, f3, f4);
	if (fx <= v2) {
		s = (fx - v1) / (v2 - v1);
		vol = factor*area1*(v2-v1)*s*s*s / 3;
		
		for(i = 0; i < n; i++) {
			func[0] = f1[i];
			func[1] = f2[i];
			func[2] = (v2-v1)*f4[i]/(v4-v1) + (v4-v2)*f1[i]/(v4-v1);
			func[3] = (v2-v1)*f3[i]/(v3-v1) + (v3-v2)*f1[i]/(v3-v1);
			out[i] = vol * ((4-3*s)*func[0] + s*func[1] + s*func[2] + s*func[3])/4;
		}
		delete[] cum1;
		return;
	}
	for(i = 0; i < n; i++) {
		func[0] = f1[i];
		func[1] = f2[i];
		func[2] = (v2-v1)*f4[i]/(v4-v1) + (v4-v2)*f1[i]/(v4-v1);
		func[3] = (v2-v1)*f3[i]/(v3-v1) + (v3-v2)*f1[i]/(v3-v1);
		cum1[i] = factor*area1*(v2-v1)*(func[0]+func[1]+func[2]+func[3])/12;
	}
	
	float A, B, C, D, r;      
	if (fx <= v3) {
		s = (fx-v2)/(v3-v2);
		float s2 = s*s;
		float s3 = s2*s;
		float s4 = s3*s;
		
		for(i = 0; i < n; i++) {
			A = f3[i]-f2[i] + (v3-v2)*(f4[i]-f1[i])/(v4-v1) + (v3-v2)*(f3[i]-f1[i])/(v3-v1);
			B = f2[i] + (v2-v1)*f4[i]/(v4-v1) + (v4-v2)*f1[i]/(v4-v1) + (v2-v1)*f3[i]/(v3-v1) + (v3-v2)*f1[i]/(v3-v1);
			r = (v3-v2)/(v2-v1);
			C = f3[i]-f2[i] + (v3-v2)*(f4[i]-f1[i])/(v4-v1) + (v3-v2)*(f4[i]-f2[i])/(v4-v2);
			D = f2[i] + (v2-v1)*f4[i]/(v4-v1) + (v4-v2)*f1[i]/(v4-v1) + f2[i];
	
			float fint1 = factor*(v3-v2)*s*(A*s*(area1*(6-8*s+3*s2)+(4-3*s)*s*midarea)+2*B*(2*area1*(3-3*s+s2)+(3-2*s)*s*midarea))/36;
			float fint2 = factor*(v3-v2)*s2*(C*s*(3*s*(area2-midarea)+4*midarea)+D*(4*s*(area2-midarea)+6*midarea))/36;
#ifdef DEBUG_VOLVOL
		printf("A = %f, B = %f, C = %f, D = %f, r = %f, t = %f, s = %f\n", A, B, C, D, r, t,s );
		printf("fint1 = %.10f, fint2 = %.10f, cum1 = %f\n", fint1, fint2, cum1+fint1+fint2);
#endif
			out[i] = cum1[i]+fint1+fint2;
		}
		delete[] cum1;
		return;
	}
	for(i = 0; i < n; i++) {
		A = f3[i]-f2[i] + (v3-v2)*(f4[i]-f1[i])/(v4-v1) + (v3-v2)*(f3[i]-f1[i])/(v3-v1);
		B = f2[i] + (v2-v1)*f4[i]/(v4-v1) + (v4-v2)*f1[i]/(v4-v1) + (v2-v1)*f3[i]/(v3-v1) + (v3-v2)*f1[i]/(v3-v1);
		r = (v3-v2)/(v2-v1);
		C = f3[i]-f2[i] + (v3-v2)*(f4[i]-f1[i])/(v4-v1) + (v3-v2)*(f4[i]-f2[i])/(v4-v2);
		D = f2[i] + (v2-v1)*f4[i]/(v4-v1) + (v4-v2)*f1[i]/(v4-v1) + f2[i];
		cum1[i] += factor*(v3-v2)*((A*(area1+midarea)+2*B*(2*area1+midarea))/36 + (2*D*(2*area2+midarea)+C*(3*area2+midarea))/36);
	}
	
	// v3 < fx < v4
	s = (fx-v3)/(v4-v3);
	for(i = 0; i < n; i++) {
		func[0] = f4[i];
		func[1] = f3[i];
		func[2] = (v4-v3)*f1[i]/(v4-v1) + (v3-v1)*f4[i]/(v4-v1);
		func[3] = (v4-v3)*f2[i]/(v4-v2) + (v3-v2)*f4[i]/(v4-v2);
		fint = factor*area2*(v4-v3)*((func[0]+func[1]+func[2]+func[3]) - 
								 (1-s)*(1-s)*(1-s)*((1+3*s)*func[0]+(1-s)*func[1]+(1-s)*func[2]+(1-s)*func[3]))/12;
		out[i] = cum1[i] + fint;
#ifdef DEBUG_VOLVOL
		printf("fx = %f, cum1 = %f\n", fx, cum1[i]+fint);
#endif
	}
	delete[] cum1;
	return;
}

