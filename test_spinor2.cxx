#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "spinor2_class.hxx"
// #include "dreh.hxx"

main(int argc, char ** argv)
{

	double nu=0.456;
	double ux=0.35;
	double uy=0.25;
	double uz= sqrt(1. - ux*ux - uy*uy);
	printf("u= %f %f %f\n", ux, uy, uz);
	printf("cos nu = %f sin nu = %f \n ", cos(nu), sin(nu));

	Spinor s;
	Rotation r,sr;
	s.unit();
	r.unit();
	s.print("init");
	r.print("init");
	
	s.srotate(  nu    );
	r.srotate(  nu    );
	s.print("s v");
	r.print("r v");
	sr = s;
	sr.print("sr v");
	exit(0);
	r=s;
	r.print();
	s.hrotate(  M_PI  );         // snake 1
	s.print("h");
	r=s;
	r.print();
	s.vrotate(  nu    );
	s.print("v");
	r=s;
	r.print();
	s.srotate(  M_PI  );         // snake 2
	s.print("s");
	r=s;
	r.print();


	double nn[3];
	s.naxis(nn);
	printf("n-axis %e %e %e\n", nn[0], nn[1], nn[2]);
	printf("angle = %f\n", s.angle()*180/M_PI);


	double v[3]= {0., 1., 0.};
	r.rotateVector(v);
	double a=s.angle();
	printf("v= %e %e %e\n", v[0], v[1], v[2]);
     

	r.unit();
	
	r.vrotate(  nu    );
	r.print();
	r.hrotate(  M_PI  );         // snake 1
	r.print();
	r.vrotate(  nu    );
	r.print();
	r.srotate(  M_PI  );         // snake 2
	r.print();

	r.unit();
	s.unit();
	for(int i=0; i<50; i++)
	{
		double ang = random()*M_PI/RAND_MAX;
// 		r.vrotate(ang);
// 		s.vrotate(ang);
		ang = random()*M_PI/RAND_MAX;
		r.hrotate(ang);
		s.hrotate(ang);
// 		ang = random()*M_PI/RAND_MAX;
// 		r.srotate(ang);
// 		s.srotate(ang);
	}
	r.print("3x3");
	r=s;
	r.print("2x2");

#define TEST1
#ifdef TEST1
	double ang= M_PI/4.;

	s.c0=0.5; s.c1=0.2; s.c2=0.3; s.c3=0.7;
	double wu = sqrt(s.c0*s.c0 + s.c1*s.c1 + s.c2*s.c2 + s.c3*s.c3);
	s.c0 /= wu; s.c1 /= wu; s.c2 /= wu; s.c3 /= wu;
	r=s;
	r.print("hin");
	r.hrotate(ang);
	r.print("hrot");
	r=s;
	r.rotate(ang, 1., 0., 0.);
	r.print("roth");
	s.hrotate(ang);
	r=s;
	r.print("hss");

	s.c0=0.5; s.c1=0.2; s.c2=0.3; s.c3=0.7;
	wu = sqrt(s.c0*s.c0 + s.c1*s.c1 + s.c2*s.c2 + s.c3*s.c3);
	s.c0 /= wu; s.c1 /= wu; s.c2 /= wu; s.c3 /= wu;
	r=s;
	r.print("vin");
	r.vrotate(ang);
	r.print("vrot");
	r=s;
	r.rotate(ang, 1., 0., 0.);
	r.print("rotv");
	s.vrotate(ang);
	r=s;
	r.print("vss");

	s.c0=0.5; s.c1=0.2; s.c2=0.3; s.c3=0.7;
	wu = sqrt(s.c0*s.c0 + s.c1*s.c1 + s.c2*s.c2 + s.c3*s.c3);
	s.c0 /= wu; s.c1 /= wu; s.c2 /= wu; s.c3 /= wu;
	r=s;
	r.print("sin");
	r.srotate(ang);
	r.print("srot");
	r=s;
	r.rotate(ang, 1., 0., 0.);
	r.print("rots");
	s.srotate(ang);
	r=s;
	r.print("sss");


#endif
	double w[3]= {0., 1., 0.};
	r.rotateVector(w);
	printf("w= %e %e %e\n", w[0], w[1], w[2]);
     


	
	
}
