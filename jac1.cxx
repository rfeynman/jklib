#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>


main(int argc, char ** argv)
{
double sigma1[4]={	   2.1628231e-05  ,  1.6919944e-06 , -5.8866938e-07 , -3.5349456e-06 };
double sigma2[4]={	   1.6919944e-06  ,  6.6981154e-07 ,  3.4054769e-06 , -1.8863767e-08 };
double sigma3[4]={	  -5.8866938e-07  ,  3.4054769e-06 ,  2.2354493e-05 ,  1.7614396e-06 };
double sigma4[4]={	  -3.5349456e-06  , -1.8863767e-08 ,  1.7614396e-06 ,  7.0839148e-07 };


double * a[5];
for(int i=0; i< 5; i++) a[i]=new double[5];
}

void jacobi(double **a,int n,double d[],double **v,int & nrot)
{
	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c,*b,*z,*vector();

	double * b = new double[n+1]; 
	double * z = new double[n+1]; 
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
			delete [] z;
			delete [] b;
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
						double g=a[i][j];
						double h=a[k][l];
						a[i][j]=g-s*(h+g*tau);
						a[k][l]=h+s*(g-h*tau);
					}
					for (j=ip+1;j<=iq-1;j++) {
						double g=a[i][j];
						double h=a[k][l];
						a[i][j]=g-s*(h+g*tau);
						a[k][l]=h+s*(g-h*tau);
					}
					for (j=iq+1;j<=n;j++) {
						double g=a[i][j];
						double h=a[k][l];
						a[i][j]=g-s*(h+g*tau);
						a[k][l]=h+s*(g-h*tau);
					}
					for (j=1;j<=n;j++) {
						double g=a[i][j];
						double h=a[k][l];
						a[i][j]=g-s*(h+g*tau);
						a[k][l]=h+s*(g-h*tau);
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

