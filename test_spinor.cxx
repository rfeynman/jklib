#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "spinor_class.hxx"

main(int argc, char ** argv)
{

	double nu=20.333*M_PI/180.;
	double ux=0.35;
	double uy=0.85;
	double uz=0.55;
	double ul=sqrt(ux*ux+uy*uy+uz*uz);
	ux /= ul;
	uy /= ul;
	uz /= ul;
	printf("u= %f %f %f\n", ux, uy, uz);

	Spinor s;
	s.unit();
	
	s.rotate(nu, ux, uy, uz);
// 	s.vrotate(M_PI);
// 	s.hrotate(nu);
// 	s.vrotate(M_PI);
// 	s.hrotate(nu);
// 	s.vrotate(M_PI);
// 	s.hrotate(nu);
// 	s.srotate(M_PI);
	s.naxis();
	s.mat();
	printf("angle = %f\n", s.angle()*180/M_PI);
	s.rotmat(nu, ux, uy, uz);




	
	
}
