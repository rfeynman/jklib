#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>



inline double square(double x) { return x*x;}

inline void ROTATE(double **a,int i,int j,int k,int l, double s, double tau)
{
	double g=a[i][j];
	double h=a[k][l];
	a[i][j]=g-s*(h+g*tau);
	a[k][l]=h+s*(g-h*tau);
}

void jacobi(double ** a,int n, double * d, double ** v,int & nrot)
{
	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c,*b,*z,*vector();

	b= new double[n];
	z= new double[n];
	for (ip=0;ip<n;ip++) {
		for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=0;ip<n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	nrot=0;
	for (i=0;i<50;i++) {
		sm=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++)
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
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && fabs(d[ip])+g == fabs(d[ip]) && fabs(d[iq])+g == fabs(d[iq]))
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
	
					for (int j=0;   j<ip;   j++)  ROTATE(a,j,ip,j,iq,s,tau);
					for (int j=ip+1;j<iq;   j++)  ROTATE(a,ip,j,j,iq,s,tau);
					for (int j=iq+1;j<n;    j++)  ROTATE(a,ip,j,iq,j,s,tau);
					for (int j=0;   j<n;    j++)  ROTATE(v,j,ip,j,iq,s,tau);
					++nrot;
				}
			}
		}
		for (ip=0;ip<n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	fprintf(stderr,"Too many iterations in routine JACOBI");
}



// a = b
void copy(double ** a, double ** b)
{
	for(int i=0; i<4; i++)
	for(int j=0; j<4; j++) 
	{
		a[i][j]=b[i][j];
	}
}



// a = b~
void transpose(double ** a, double ** b)
{
	for(int i=0; i<4; i++)
	for(int j=0; j<4; j++) 
	{
		a[i][j]=b[j][i];
	}
}



// a = 1
void unit(double ** a)
{
	for(int i=0; i<4; i++)
	for(int j=0; j<4; j++) 
	{
		a[i][j]=0.;
	}
	for(int i=0; i<4; i++)
		a[i][i]=1.;
}


// a = b*c
double mabs(double ** a)
{
	double r=0.;
	for(int i=0; i<4; i++)
	for(int j=0; j<4; j++) 
	{
		 r += square(a[i][j]);
	}
	return r;
}


// a = b*c
void mult_mm(double ** a, double ** b, double **c)
{
	for(int i=0; i<4; i++)
	for(int j=0; j<4; j++) 
	{
		a[i][j]=0.;
		for(int k=0; k<4; k++) a[i][j] += b[i][k]*c[k][j];
	}
}


// a = b~*c
void mult_mtm(double ** a, double ** b, double **c)
{
	for(int i=0; i<4; i++)
	for(int j=0; j<4; j++) 
	{
		a[i][j]=0.;
		for(int k=0; k<4; k++) a[i][j] += b[k][i]*c[k][j];
	}
}


void subtract(double ** a, double ** b, double **c)
{
	for(int i=0; i<4; i++)
	for(int j=0; j<4; j++) 
	{
		a[i][j] = b[i][j]-c[i][j];
	}
}


void print_m(char * s,  double ** a, int n)
{
	printf("%s=\n", s);
	for(int i=0; i<n; i++)
	{
		for(int j=0; j<n; j++) printf("%e  ", a[i][j]);
		printf("\n");
	}
	printf("\n\n");
}



double ** matrix(int n)
{
	double ** a = new double*[n];
	for(int i=0; i<n; i++) a[i]=  new double[n];
	return a;
}


double det(double ** a, int n)
{
	if(n == 1) return a[0][0];
	double d=0.;
	double sign=1.;
	double **b =matrix(n-1);
	for(int i=0; i<n; i++)
	{
		for(int j=0; j<i; j++)
		for(int k=1; k<n; k++) b[j][k-1]=a[j][k];
		for(int j=i+1; j<n; j++)
		for(int k=1; k<n; k++) b[j-1][k-1]=a[j][k];
		d += a[i][0]*det(b,n-1)*sign;
		sign= -sign;
	}
// 	printf("det=====\n");
// 	for(int i=0; i<n; i++)
// 	{
// 		for(int j=0; j<n; j++) printf("%e  ", a[i][j]);
// 		printf("\n");
// 	}
// 	printf("det=>%e\n========\n\n", d);
	return d;
}

main(int argc, char ** argv)
{

	double sigma1[4]={	   2.1628231e-05  ,  1.6919944e-06 , -5.8866938e-07 , -3.5349456e-06 };
	double sigma2[4]={	   1.6919944e-06  ,  6.6981154e-07 ,  3.4054769e-06 , -1.8863767e-08 };
	double sigma3[4]={	  -5.8866938e-07  ,  3.4054769e-06 ,  2.2354493e-05 ,  1.7614396e-06 };
	double sigma4[4]={	  -3.5349456e-06  , -1.8863767e-08 ,  1.7614396e-06 ,  7.0839148e-07 };

// double sigma1[4]={ 5.455409e-06,   3.071146e-06,   5.001294e-07,   2.861242e-06 };
// double sigma2[4]={ 3.071146e-06,   1.206978e-05,  -1.875257e-05,   5.514145e-06 };
// double sigma3[4]={ 5.001294e-07,  -1.875257e-05,   3.508139e-05,  -6.922559e-06 };
// double sigma4[4]={ 2.861242e-06,   5.514145e-06,  -6.922559e-06,   2.974086e-06 };

// double sigma1[4]={   5.110123e-07,  -4.747946e-11,   1.811352e-08,   1.895015e-10};
// double sigma2[4]={  -4.747946e-11,   4.383729e-11,  -2.088349e-10,   6.518844e-13};
// double sigma3[4]={   1.811352e-08,  -2.088349e-10,   4.912626e-07,   4.644163e-11};
// double sigma4[4]={   1.895015e-10,   6.518844e-13,   4.644163e-11,   4.286299e-11};


	double ** sigma = matrix(4);
	double ** s = matrix(4);
	for(int i=0; i<4; i++)
	{
		for(int j=0; j<4; j++) s[i][j]=0.;
		sigma[0][i]=sigma1[i];
		sigma[1][i]=sigma2[i];
		sigma[2][i]=sigma3[i];
		sigma[3][i]=sigma4[i];
	}
	s[0][1]= -1.;
	s[2][3]= -1.;
	s[1][0]=  1.;
	s[3][2]=  1.;


	print_m("sigma", sigma, 4);

	double ** a = matrix(4);
	copy(a, sigma);
	double ** eigen = matrix(4);
	int nrot;
	double *d = new double[4];
	jacobi(a, 4, d, eigen, nrot);

	print_m("eigen", eigen, 4);
	double ** eigenT = matrix(4);
	transpose(eigenT, eigen);

	double ** c = matrix(4);

	printf("eigenvalues=\n");
	for(int j=0; j<4; j++) printf("%e  ", d[j]);
	printf("\n");
	printf("\n\n");



	printf("vd=\n");
	for(int i=0; i<4; i++)
	{
		for(int j=0; j<4; j++) printf("%e  ", eigen[i][j]*d[j]);
		printf("\n");
	}
	printf("\n\n");
	mult_mm(a,sigma,eigen);
	print_m("sigma*eigen", a, 4);

	mult_mm(a,sigma,eigen);
	mult_mm(c,eigenT,a);
	
	print_m("c= eigenT*sigma*eigen", c, 4);

	double bx=sqrt(c[0][0]/1e-6);
	double by=sqrt(c[2][2]/1e-6);


	unit(c);
	c[0][0]=1/bx;
	c[1][1]=bx;
	c[2][2]=1/by;
	c[3][3]=by;

	double ** z = matrix(4);
	mult_mm(z, c, eigenT);
	// 	double sol=16.;
	// 	unit(z);
	// 	z[1][2]=  sol;
	// 	z[3][0]= -sol;
	double ** zT = matrix(4);
	transpose(zT,z);


	mult_mm(a,sigma,zT);
	mult_mm(c,z,a);
	

	print_m("c= z*sigma*zT", c, 4);

	for(double sol= -20.; sol<= 20.; sol += 0.1)
	// for(double sol= -2.258; sol<= -2.256; sol += 0.00001)
	// double sol=16.;
	{

		unit(c);
		c[1][2]=  sol;
		c[3][0]= -sol;
	// 	print_m("c", c, 4);
		double ** zr = matrix(4);
		mult_mm(zr, c, z);
	// 	copy(zr, c);

		mult_mm(c,s,zr);
		mult_mm(a,zr,c);
		subtract(c,a,s);
	// 	print_m("sT*zr*s -s", c, 4);
	// 	printf(" sol %f det c=%e\n", sol, det(c,4));
		printf("  %e %e\n", sol, mabs(c));
	}
}


// make symplectisch
//
// m:= (3*m + m*s*mt*s*m)/2 until m=symplectic
//
//
//
// mt*s*m== s;
