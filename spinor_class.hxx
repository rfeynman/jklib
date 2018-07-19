#ifndef __spinor_class__
#define __spinor_class__
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

class Spinor
{
	double c0, c1, c2, c3;
	double mat3[3][3];
	double n[3];

public:
	void unit()
	{
		c0 =1;
		c1 = c2 = c3 = 0.;
	}
	
	void rotmat(double phi, double ux, double uy, double uz)
	{
		// matrix of a rotation by the angle phi around the axis given by the unit vector (ux,uy,uz)
		// from http://en.wikipedia.orh/wiki/Rotation_matrix
		// this does not use the spinor formalism but can be used to compare the two methods
		double co = cos(phi);
		double si = sin(phi);
		double co1 = 1. - co;

		mat3[0][0] = co + ux*ux*co1;
		mat3[0][1] = ux*uy*co1 - uz*si;
		mat3[0][2] = ux*uz*co1 + uy*si;
		mat3[1][0] = ux*uy*co1 + uz*si;
		mat3[1][1] = co + uy*uy*co1;
		mat3[1][2] = uy*uz*co1 - ux*si;
		mat3[2][0] = ux*uz*co1 - uy*si;
		mat3[2][1] = uy*uz*co1 + ux*si;
		mat3[2][2] = co + uz*uz*co1;
		printf("R= %f %f %f\n   %f %f %f\n    %f %f %f\n ", mat3[0][0], mat3[0][1], mat3[0][2],
								    mat3[1][0], mat3[1][1], mat3[1][2],
								    mat3[2][0], mat3[2][1], mat3[2][2]);
	}
	
	void rotate(double phi, double ux, double uy, double uz)
	{
		//  add a rotation given by the angle phi around the axis given by the unit vector (ux,uy,uz) to this spinor
		double co = cos(phi/2.);
		double si = sin(phi/2.);

		double c0 = co;
		double c1 = si * ux;
		double c2 = si * uy;
		double c3 = si * uz;
		crotate(c0, c1, c2, c3);
	}
	void hrotate(double phi)
	{
		//  add a rotation given by the angle phi around the x axis  to this spinor
		double co = cos(phi/2.);
		double si = sin(phi/2.);
		crotate(co, 0., si, 0.);
	}
	void vrotate(double phi)
	{
		//  add a rotation given by the angle phi around the y axis  to this spinor
		double co = cos(phi/2.);
		double si = sin(phi/2.);
		crotate(co, si, 0., 0.);
	}
	void srotate(double phi)
	{
		//  add a rotation given by the angle phi around the s (z) axis  to this spinor
		double co = cos(phi/2.);
		double si = sin(phi/2.);
		crotate(co, 0., 0., si);
	}
	void crotate(double k0, double k1, double k2, double k3)
	{
		// multiply this spinor with a spinor (k0, k1, k2, k3)
		double x0 = k0*c0 - k3*c3 - k2*c2 - k1*c1;
		double x1 = k1*c0 + k2*c3 + k3*c2 + k0*c1;
		double x2 = k2*c0 - k1*c3 + k0*c2 - k3*c1;
		double x3 = k3*c0 + k0*c3 + k1*c2 - k2*c1;

		c0 = x0; c1 = x1; c2 = x2; c3 = x3;
	}

	void naxis()
	{
		double rn = sqrt( c1*c1 + c2*c2 + c3*c3 );
		n[0] = c1/rn;
		n[1] = c2/rn;
		n[2] = c3/rn;
		printf("n= %f %f %f\n", n[0], n[1], n[2]);
	}

	double angle()
	{
		double s = sqrt(c1*c1+c2*c2+c3*c3);
		return 2.*atan2(s, c0);
// 		return 2.*acos(c0);
	}

	double mat()
	{
		mat3[0][0] = c0*c0 + c1*c1 - c2*c2 - c3*c3;
		mat3[0][1] = 2.*(c1*c2-c0*c3);
		mat3[0][2] = 2.*(c1*c3+c0*c2);
		mat3[1][0] = 2.*(c1*c2+c0*c3);
		mat3[1][1] = c0*c0 - c1*c1 + c2*c2 - c3*c3;
		mat3[1][2] = 2.*(c2*c3-c0*c1);
		mat3[2][0] = 2.*(c1*c3-c0*c2);
		mat3[2][1] = 2.*(c2*c3+c0*c1);
		mat3[2][2] = c0*c0 - c1*c1 - c2*c2 + c3*c3;
		printf("m= %f %f %f\n   %f %f %f\n    %f %f %f\n ", mat3[0][0], mat3[0][1], mat3[0][2],
								    mat3[1][0], mat3[1][1], mat3[1][2],
								    mat3[2][0], mat3[2][1], mat3[2][2]);
	}

	void rvect(double * v)
	{
		double v1 = mat3[0][0]*v[0]+mat3[0][1]*v[1]+mat3[0][2]*v[2];
		double v2 = mat3[1][0]*v[0]+mat3[1][1]*v[1]+mat3[1][2]*v[2];
		double v3 = mat3[2][0]*v[0]+mat3[2][1]*v[1]+mat3[2][2]*v[2];

		v[0] = v1; v[1] = v2; v[2] = v3;
		printf("v= %f %f %f\n", v[0], v[1], v[2]);
	}
};

#endif
