#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "rungeKutta_class.hxx"

inline double square(double x) { return x*x; }
inline double cube(double x) { return x*x*x; }


double nup[100], zw[100];
int nw;

class MyRu : public RungeKutta
{
	double M;
	double en;
	double beta;
	double gamma;
	double ks;
	double kr;
public:
	MyRu() : RungeKutta(2)
	{
		M =1.500;
		en=0.e-1;
		gamma=10.;
		beta=sqrt(1.-1./square(gamma));
		ks=14000./(cube(beta*gamma)*17000.*2.);
		ks=1;
		kr=0.00;
	}

               virtual void derivs(double x,           double * y,         double * dydx       )
	       {
		       dydx[0]=y[1];
		       dydx[1]= -kr*y[0] + square((M+en)/(beta*gamma))/cube(y[0]) + ks/y[0];
		       if(y[1]<= 0.) zw[nw]=x;
// 		       printf("%e %e %e %e %e\n", x, y[0], y[1], dydx[0],dydx[1]);
	       }
} myRu;

// 72:el=6 slice=0 sigr=1.341681e-04 sig_rp=2.455441e-04 gamma=1.033100e+01 current=4.100000e-01
// 75:el=6 slice=1 sigr=2.218199e-04 sig_rp=3.984151e-04 gamma=1.033449e+01 current=1.150000e+00
// 78:el=6 slice=2 sigr=2.372778e-04 sig_rp=4.187865e-04 gamma=1.030908e+01 current=1.330000e+00
// 81:el=6 slice=3 sigr=2.242985e-04 sig_rp=4.003873e-04 gamma=1.025903e+01 current=1.290000e+00
// 84:el=6 slice=4 sigr=2.089419e-04 sig_rp=3.908010e-04 gamma=1.016934e+01 current=8.200000e-01
//

main(int argc, char ** argv)
{
	int nok, nbad;
	double ystart[2] ={1., -0.1};
	double dx=0.005;

	nw=0;
	for(int i=1; i<20; i++)
	{
		ystart[0]=1;
		ystart[1]= -1.-0.05*i;
		nup[nw]=ystart[1];
		for(double x1=0., x2=dx; x1<2000*dx; x1 = x2, x2 += dx)
		{
// 	                   a.odeint( ystart,double x1,double x2,double eps,double h1,double hmin,int & nok,int & nbad);
			myRu.odeint( ystart,       x1,     x2,     1.e-9  , 1.e-6    ,   1.e-18   ,      nok,      nbad);
			printf("%e %e %e\n", x2, ystart[0], ystart[1]);
		}
		printf("&\n");
		fprintf(stderr,"%f %f\n", nup[nw], zw[nw]);
		nw++;
	}
}


