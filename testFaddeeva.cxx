#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "Faddeeva.hxx"

main(int argc, char ** argv)
{
// 	FILE * fb = fopen(argv[1],"r");
// 	if(!fb)
// 	{
// 		perror(argv[1]);
// 		exit(-1);
// 	}


	std::complex<double> z(0.001, -0.002);
	std::complex<double> f = Faddeeva::w(z, 1e-10);
	double r = f.real();
	double i = f.imag();
	printf("%e %e\n", r,i);
}
