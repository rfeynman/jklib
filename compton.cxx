#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "complex_class.hxx"


double cross(Complex * e1, Complex * e2)
{
	double ww = 2.;
	Complex e1e2 = e1[0]*e2[0] + e1[1]*e2[1];
	double e1e2abs=e1e2.r*e1e2.r+e1e2.i*e1e2.i;
	Complex e1e2star = e1[0]*conjugate(e2[0]) + e1[1]*conjugate(e2[1]);
	double e1e2starabs=e1e2star.r*e1e2star.r+e1e2star.i*e1e2star.i;
	double cr = (1. + e1e2starabs -e1e2abs) * ww + 2. * ( e1e2abs + e1e2starabs -1.);
	return cr;
}

main(int argc, char ** argv)
{
	Complex e1[2], e2[2];
	e1[0]= Complex(1.,0.);
	e1[1]= Complex(0.,1.);
	e2[0]= Complex(1.,0.);
	e2[1]= Complex(0.,-1.);
	printf("%e\n", cross(e1,e2));

}
