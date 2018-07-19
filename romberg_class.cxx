#include "romberg_class.hxx"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


Romberg::Romberg(double _EPS, int _JMAX, int _K)
{
	EPS= _EPS;
	JMAX = _JMAX;
	JMAXP = JMAX +1;
	K = _K;
}



// Returns the integral of the function func from a to b. Integration is performed by Romberg’s
// method of order 2K, where, e.g., K=2 is Simpson’s rule.
double Romberg::qromb(double a, double b)
{
	double ss,dss;
	double s[JMAXP],h[JMAXP+1]; 	// These store the successive trapezoidal approximations and their relative stepsizes.

	h[0]=1.0;
	for (int j=0;j<JMAX;j++)
	{
		s[j]=trapzd(a,b,j);
// 		printf("j=%d  s=%e\n", j, s[j]);
		if (j+1 >= K)
		{
			polint(&h[j+1-K],&s[j+1-K],K,0.0,ss,dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=0.25*h[j];
	}
	fprintf(stderr, "Too many steps in routine qromb\n");
	return 0.0;
}


// Given arrays xa[1..n] and ya[1..n], and given a value x, this routine returns a value y, and
// an error estimate dy. If P(x) is the polynomial of degree N − 1 such that P(xa[i] ) = ya[i] , 
// i = 1, . . . , n, then the returned value y = P(x)
void Romberg::polint(double xa[], double ya[], int n, double x, double &y, double &dy)
{
	double * c= new double[n];
	double * d= new double[n];

// 	for(int i=0; i<n; i++) printf("%e  ", xa[i]); printf("\n");
// 	for(int i=0; i<n; i++) printf("%e  ", ya[i]); printf("\n");
	int ns=0;
	double dif=fabs(x-xa[0]);
	for (int i=0;i<n;i++)
	{
		double dift=fabs(x-xa[i]);
		if (dift < dif)
		{
			ns=i;
			dif=dift;
		}
		c[i] = d[i] = ya[i];
	}
	y=ya[ns];
// 	printf("ns=%d y = %e\n", ns, y);
	for (int m=1; m<n; m++)
	{
		for (int i=0;i<=n-m-1;i++)
		{
// 			printf("i=%d m=%d i+m = %d\n", i, m, i+m);
			double ho=xa[i]-x;
			double hp=xa[i+m]-x;
			double w=c[i+1]-d[i];
			double den=ho-hp;
			if ( den == 0.0)
			{
				fprintf(stderr, "Error in routine polint\n");
				exit(-1);
			}
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
// 			printf("c[%d] %e %e \n", i, c[i], d[i]);
		}
		dy=(2*ns < (n-m) ? c[ns] : d[ns-1]);
		ns--;
		y += dy;
// 		printf("y %e %e %d\n", y, dy, ns);
	}
	delete [] c;
	delete [] d;
}

// This routine computes the nth stage of refinement of an extended trapezoidal rule. func is 
// the function to be integrated between limits a and b. When called with 
// n=1, the routine returns the crudest estimate of integral_a^b f (x)dx.
// Subsequent calls with n=2,3,...
// (in that sequential order) will improve the accuracy by adding 2^(n-2) additional interior points.
double Romberg::trapzd(double a, double b, int n)
{

	if (n == 0)
	{
		trap_s=0.5*(b-a)*(func(a)+func(b));
	}
	else
	{
		int it=1;
		for (int j=1; j<n; j++) it <<= 1;
// 		printf("it = %d\n", it);
		double tnm=it;
		double del=(b-a)/tnm;
		double x=a+0.5*del;
		double sum =0.0;
		for (int j=0; j<it; j++,x+=del) sum += func(x);
		trap_s=0.5*(trap_s+(b-a)*sum/tnm);
	}
	return trap_s;
}
