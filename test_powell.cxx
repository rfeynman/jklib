#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "powell_class.hxx"



class mypowell : public powell
{
public:
	virtual double func(double p[]);
	mypowell() : powell(3) {}
};

double mypowell::func(double p[])
{
	double x1 = p[0] -5.;
	double x2 = p[1] -4.;
	double x3 = p[2] -7.;
	return (  x1*x1 + x2 * x2 + x3*x3);
}


main(int argc, char ** argv)
{
	mypowell mp;
	double p[3] = { 6., 5., 8.};
	double ftol= 0.0001; int iter; double fret;

	mp.match(p, ftol, iter, fret);
	printf("%f %f %f\n", p[0], p[1], p[2]);
	printf("%f %f %d\n", ftol , fret, iter);
}
