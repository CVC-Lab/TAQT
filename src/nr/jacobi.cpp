#include <math.h>
#include "nr.h"
#include "nrutil.h"
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

/**
Computes all eigenvalues and eigenvectors of a real symmetric matrix a[1..n][1..n]. On
output, elements of a above the diagonal are destroyed. d[1..n] returns the eigenvalues of a.
v[1..n][1..n] is a matrix whose columns contain, on output, the normalized eigenvectors of
a. nrot returns the number of Jacobi rotations that were required. 
*/
void jacobi(float **a, int n, float d[], float **v, int *nrot)
{
	int j,iq,ip,i;
	float tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
	b=fvector(1,n);
	z=fvector(1,n);
	for (ip=1;ip<=n;ip++) { //Initialize to the identity matrix.
		for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=1;ip<=n;ip++) { //Initialize b and d to the diagonal of a. 
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0; //This vector will accumulate terms of the form tapq as in equation (11.1.14).
	}
	*nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		
		for (ip=1;ip<=n-1;ip++) { //Sum off-diagonal elements.
			for (iq=ip+1;iq<=n;iq++){
				sm += fabs(a[ip][iq]);
			}
		}
						
		if (sm == 0.0) { //The normal return, which relies on quadratic convergence to machine underflow.
			free_vector(z,1,n);
			free_vector(b,1,n);
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n); // ...on the first three sweeps.
		else
			tresh=0.0; //...thereafter.
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				//After four sweeps, skip the rotation if the off-diagonal element is small.
				if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])
					&& (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
					a[ip][iq]=0.0;

				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((float)(fabs(h)+g) == (float)fabs(h))
						t=(a[ip][iq])/h; // t = 1/(2?)
					else {
						theta=0.5*h/(a[ip][iq]); // Equation (11.1.10).
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=1;j<=ip-1;j++) { //Case of rotations 1 ? j < p.
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) { //Case of rotations p < j < q.
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++) { //Case of rotations q < j ? n.
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
					//printf("%f %f %f\n", d[1], a[1][2], a[1][3]);
					//printf("%f %f %f\n", a[1][2], d[2], a[2][3]);
					//printf("%f %f %f\n", a[1][3], a[2][3], d[3]);
				}
			}
		}
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip]; //Update d with the sum of tapq,
			z[ip]=0.0; //and reinitialize z.
		}
		
		
	}
	printf("sm = %.10f\n", sm);
	nrerror("Too many iterations in routine jacobi");
}

/**
 * Given the eigenvalues d[1..n] and eigenvectors v[1..n][1..n] as output from jacobi
 * (§11.1) or tqli (§11.3), this routine sorts the eigenvalues into descending order, and rearranges
 * the columns of v correspondingly. The method is straight insertion.
 */
 void eigsrt(float d[], float **v, int n)
 {
	 int k,j,i;
	 float p;
	 for (i=1;i<n;i++) {
		 p=d[k=i];
		 for (j=i+1;j<=n;j++)
			 if (d[j] >= p) p=d[k=j];
		 if (k != i) {
			 d[k]=d[i];
			 d[i]=p;
			 for (j=1;j<=n;j++) {
				 p=v[j][i];
				 v[j][i]=v[j][k];
				 v[j][k]=p;
			 }
		 }
	 }
 }
