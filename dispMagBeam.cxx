#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "rungeKutta_class.hxx"

inline double square(double x) { return x*x; }
inline double cube(double x) { return x*x*x; }
const int nmag=2;


double nup[100], zw[100];
int nw;
//
double slice0[4]={ 1.158309e-02, 1.566725e-02, 1.033100e+01, 4.100000e-01};
double slice1[4]={ 1.489362e-02, 1.995890e-02, 1.033449e+01, 1.150000e+00};
double slice2[4]={ 1.540382e-02, 2.046198e-02, 1.030908e+01, 1.330000e+00};
double slice3[4]={ 1.497660e-02, 2.000561e-02, 1.025903e+01, 1.290000e+00};
double slice4[4]={ 1.445482e-02, 1.976402e-02, 1.016934e+01, 8.200000e-01};


double sigr[5];
double sigp[5];
double gamm[5];
double curr[5];


void cpy()
{
sigr[0]= slice0[0];
sigr[1]= slice1[0];
sigr[2]= slice2[0];
sigr[3]= slice3[0];
sigr[4]= slice4[0];


sigp[0]= slice0[1];
sigp[1]= slice1[1];
sigp[2]= slice2[1];
sigp[3]= slice3[1];
sigp[4]= slice4[1];


gamm[0]= slice0[2];
gamm[1]= slice1[2];
gamm[2]= slice2[2];
gamm[3]= slice3[2];
gamm[4]= slice4[2];


curr[0]= slice0[3];
curr[1]= slice1[3];
curr[2]= slice2[3];
curr[3]= slice3[3];
curr[4]= slice4[3];
}


class MyRu : public RungeKutta
{
	double M;
	double en;
	double beta;
	double ks;
	double kr;
	double km;
	double current;
public:
	double gamma;
	double curr;
	double mag[nmag][4];

	MyRu() : RungeKutta(2)
	{
		M=0;
		M =1.500e-3;
		en=0.e-1;
		gamma=10.;
		current=14.02;
		kr=0.00;
	}

	void setup()
	{
		beta=sqrt(1.-1./square(gamma));
		ks=current*curr/(cube(beta*gamma)*17000.*2.);
		km=square((M+en)/(beta*gamma));
// 		for(int i=0; i<nmag; i++) mag[i][3] = elchar/(elmass*gamma)*mag[i][2];
	}

        virtual void derivs(double x,           double * y,         double * dydx       )
	{
		kr =0.;
		for(int i=0; i<nmag; i++)
			if(x >= mag[i][0] && x <= mag[i][1]) kr = mag[i][3];

	       	dydx[0]=y[1];
	       	dydx[1]= -kr*y[0] + km/cube(y[0]) + ks/y[0];
	       	if(y[1]<= 0.) zw[nw]=x;
// 	        printf("%e %e %e %e %e\n", x, y[0], y[1], dydx[0],dydx[1]);
	}
} myRu;

double funk(double p[])
{
	int nok, nbad;
	double ystart[2];
	double dx=0.001;
	cpy();

	myRu.mag[0][0] = 0.1;
	myRu.mag[0][1] = 0.2;
	myRu.mag[0][2] = 1000.;
	myRu.mag[0][3] = 45.2+p[2];

	myRu.mag[1][0] = 1.1+p[0];
	myRu.mag[1][1] = 1.2+p[0];
	myRu.mag[1][2] = 1000.;
	myRu.mag[1][3] = 45.2+p[3];

	nw=0;
	for(int i=1; i<5; i++)
	{
		ystart[0]=sigr[i];
		ystart[1]=sigp[i];
		myRu.gamma=gamm[i];
		myRu.curr=curr[i];
		myRu.setup();

		nup[nw]=ystart[1];
		printf("0 %e\n",  ystart[0]);
		for(double x1=0., x2=dx; x1<2000*dx; x1 = x2, x2 += dx)
		{
// 	                   a.odeint( ystart,double x1,double x2,double eps,double h1,double hmin,int & nok,int & nbad);
			myRu.odeint( ystart,       x1,     x2,     1.e-9  , 1.e-6    ,   1.e-18   ,      nok,      nbad);
			printf("%e %e\n", x2, ystart[0]);
		}
		printf("&\n");
		fprintf(stderr,"%f %f\n", nup[nw], zw[nw]);
		nw++;
	}
	return 0.;
}


main(int argc, char ** argv)
{
	double p[4] = {0., 0., 0., 0.};
	fprintf(stderr,"f=%e\n", funk(p));
}
