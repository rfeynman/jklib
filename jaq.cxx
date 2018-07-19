#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>




#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

void jacobi(double ** a,int n, double * d, double ** v,int & nrot)
{
	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c,*b,*z,*vector();

	b= ( new double[n]) -1;
	z= ( new double[n]) -1;
	for (ip=1;ip<=n;ip++) {
		for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=1;ip<=n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			delete [] (z+1);
			delete [] (b+1);
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && fabs(d[ip])+g == fabs(d[ip])
					&& fabs(d[iq])+g == fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if (fabs(h)+g == fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
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
					for (j=1;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++nrot;
				}
			}
		}
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	fprintf(stderr,"Too many iterations in routine JACOBI");
}

#undef ROTATE
main(int argc, char ** argv)
{
	double sigma1[4]={	   2.1628231e-05  ,  1.6919944e-06 , -5.8866938e-07 , -3.5349456e-06 };
	double sigma2[4]={	   1.6919944e-06  ,  6.6981154e-07 ,  3.4054769e-06 , -1.8863767e-08 };
	double sigma3[4]={	  -5.8866938e-07  ,  3.4054769e-06 ,  2.2354493e-05 ,  1.7614396e-06 };
	double sigma4[4]={	  -3.5349456e-06  , -1.8863767e-08 ,  1.7614396e-06 ,  7.0839148e-07 };

// 	double sigma1[4]={	  1.  ,  8. ,  2. ,  3. };
// 	double sigma2[4]={	  8.  ,  7. ,  0. ,  4. };
// 	double sigma3[4]={	  2.  ,  0. ,  7. ,  4. };
// 	double sigma4[4]={	  3.  ,  4. ,  4. ,  5. };


	double ** a = (new double*[4])-1;
	double ** b = (new double*[4])-1;
	double ** c = (new double*[4])-1;
	double ** v = (new double*[4])-1;
	int nrot;
	double *d = (new double[4])-1;
	for(int i=1; i<= 4; i++) a[i]=  (new double[4])-1;
	for(int i=1; i<= 4; i++) b[i]=  (new double[4])-1;
	for(int i=1; i<= 4; i++) c[i]=  (new double[4])-1;
	for(int i=1; i<= 4; i++) v[i]=  (new double[4])-1;
	for(int i=1; i<= 4; i++)
	{
		a[1][i]=sigma1[i-1];
		a[2][i]=sigma2[i-1];
		a[3][i]=sigma3[i-1];
		a[4][i]=sigma4[i-1];
		b[1][i]=sigma1[i-1];
		b[2][i]=sigma2[i-1];
		b[3][i]=sigma3[i-1];
		b[4][i]=sigma4[i-1];
	}

	printf("b=\n");
	for(int i=1; i<= 4; i++)
	{
		for(int j=1; j<= 4; j++) printf("%e  ", b[i][j]);
		printf("\n");
	}
	printf("\n\n");

	jacobi(a, 4, d, v, nrot);

	printf("a=\n");
	for(int i=1; i<= 4; i++)
	{
		for(int j=1; j<= 4; j++) printf("%e  ", a[i][j]);
		printf("\n");
	}
	printf("\n\n");

	printf("v=\n");
	for(int i=1; i<= 4; i++)
	{
		for(int j=1; j<= 4; j++) printf("%e  ", v[i][j]);
		printf("\n");
	}
	printf("\n\n");

	printf("d=\n");
	for(int j=1; j<= 4; j++) printf("%e  ", d[j]);
	printf("\n");
	printf("\n\n");

	for(int i=1; i<= 4; i++)
	for(int j=1; j<= 4; j++) 
	{
		a[i][j]=0.;
		for(int k=1; k<= 4; k++) a[i][j] += b[i][k]*v[k][j];
	}

	printf("a=\n");
	for(int i=1; i<= 4; i++)
	{
		for(int j=1; j<= 4; j++) printf("%e  ", a[i][j]);
		printf("\n");
	}
	printf("\n\n");

	printf("vd=\n");
	for(int i=1; i<= 4; i++)
	{
		for(int j=1; j<= 4; j++) printf("%e  ", v[i][j]*d[j]);
		printf("\n");
	}
	printf("\n\n");

	for(int i=1; i<= 4; i++)
	for(int j=1; j<= 4; j++) 
	{
		c[i][j]=0.;
		for(int k=1; k<= 4; k++) c[i][j] += v[k][i]*a[k][j];
	}

	printf("c=\n");
	for(int i=1; i<= 4; i++)
	{
		for(int j=1; j<= 4; j++) printf("%e  ", c[i][j]);
		printf("\n");
	}
	printf("\n\n");

}


#undef ROTATE

