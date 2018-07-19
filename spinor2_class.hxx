#ifndef __spinor_class__
#define __spinor_class__
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

class Rotation;

class Spinor
{
public:
	double c0, c1, c2, c3;
	inline void unit();
	inline void rotate(double phi, double ux, double uy, double uz);
	inline void hrotate(double phi);
	inline void vrotate(double phi);
	inline void srotate(double phi);
	inline double angleNaxis(double * n);
	inline double det();
	inline double print(const char * text = "spinor");
	inline double operator=(Rotation & S);
};



class Rotation
{
public:
	double mat[3][3];
	inline double det();
	inline double angleNaxis(double * n);
	inline void rotateVector(double * v);
	inline double print(const char * text = "mat");
	inline double operator=(Spinor & S);
	inline void srotate(double phi);
	inline void vrotate(double phi);
	inline void hrotate(double phi);
	inline void rotate(double phi, double ux, double uy, double uz);
	inline void multiply(double m[3][3]);
	inline void unit();
};

inline void Rotation::unit()
{
	mat[0][0] = 1.;
	mat[0][1] = 0.;
	mat[0][2] = 0.;
	mat[1][0] = 0.;
	mat[1][1] = 1.;
	mat[1][2] = 0.;
	mat[2][0] = 0.;
	mat[2][1] = 0.;
	mat[2][2] = 1.;
}

inline void Rotation::multiply(double m[3][3])
{
	double tmp00 = m[0][0]*mat[0][0] + m[0][1]*mat[1][0] + m[0][2]*mat[2][0];
	double tmp01 = m[0][0]*mat[0][1] + m[0][1]*mat[1][1] + m[0][2]*mat[2][1];
	double tmp02 = m[0][0]*mat[0][2] + m[0][1]*mat[1][2] + m[0][2]*mat[2][2];
	double tmp10 = m[1][0]*mat[0][0] + m[1][1]*mat[1][0] + m[1][2]*mat[2][0];
	double tmp11 = m[1][0]*mat[0][1] + m[1][1]*mat[1][1] + m[1][2]*mat[2][1];
	double tmp12 = m[1][0]*mat[0][2] + m[1][1]*mat[1][2] + m[1][2]*mat[2][2];
	double tmp20 = m[2][0]*mat[0][0] + m[2][1]*mat[1][0] + m[2][2]*mat[2][0];
	double tmp21 = m[2][0]*mat[0][1] + m[2][1]*mat[1][1] + m[2][2]*mat[2][1];
	double tmp22 = m[2][0]*mat[0][2] + m[2][1]*mat[1][2] + m[2][2]*mat[2][2];

	mat[0][0] = tmp00;
	mat[0][1] = tmp01;
	mat[0][2] = tmp02;
	mat[1][0] = tmp10;
	mat[1][1] = tmp11;
	mat[1][2] = tmp12;
	mat[2][0] = tmp20;
	mat[2][1] = tmp21;
	mat[2][2] = tmp22;
}

inline void Rotation::rotate(double phi, double ux, double uy, double uz)
{
	// matrix of a rotation by the angle phi around the axis given by the unit vector (ux,uy,uz)
	// from http://en.wikipedia.orh/wiki/Rotation_matrix
	double co = cos(phi);
	double si = sin(phi);
	double co1 = 1. - co;

	double mat00 = co + ux*ux*co1;
	double mat01 = ux*uy*co1 - uz*si;
	double mat02 = ux*uz*co1 + uy*si;
	double mat10 = ux*uy*co1 + uz*si;
	double mat11 = co + uy*uy*co1;
	double mat12 = uy*uz*co1 - ux*si;
	double mat20 = ux*uz*co1 - uy*si;
	double mat21 = uy*uz*co1 + ux*si;
	double mat22 = co + uz*uz*co1;

	double tmp00 = mat00*mat[0][0] + mat01*mat[1][0] + mat02*mat[2][0];
	double tmp01 = mat00*mat[0][1] + mat01*mat[1][1] + mat02*mat[2][1];
	double tmp02 = mat00*mat[0][2] + mat01*mat[1][2] + mat02*mat[2][2];
	double tmp10 = mat10*mat[0][0] + mat11*mat[1][0] + mat12*mat[2][0];
	double tmp11 = mat10*mat[0][1] + mat11*mat[1][1] + mat12*mat[2][1];
	double tmp12 = mat10*mat[0][2] + mat11*mat[1][2] + mat12*mat[2][2];
	double tmp20 = mat20*mat[0][0] + mat21*mat[1][0] + mat22*mat[2][0];
	double tmp21 = mat20*mat[0][1] + mat21*mat[1][1] + mat22*mat[2][1];
	double tmp22 = mat20*mat[0][2] + mat21*mat[1][2] + mat22*mat[2][2];

	mat[0][0] = tmp00;
	mat[0][1] = tmp01;
	mat[0][2] = tmp02;
	mat[1][0] = tmp10;
	mat[1][1] = tmp11;
	mat[1][2] = tmp12;
	mat[2][0] = tmp20;
	mat[2][1] = tmp21;
	mat[2][2] = tmp22;
}

inline void Rotation::hrotate(double phi)
{
	// rotate mat around the horizontal axis 
	//
	//	[   1   0   0 ]   [  m00  m01  m02 ]   [      m00           m01           m02      ]  
	//	[   0   c  -s ] * [  m10  m11  m12 ] = [  c*m10-s*m20   c*m11-s*m21   c*m12-s*m22  ]
	//	[   0   s   c ]   [  m20  m21  m22 ]   [  s*m10+c*m20   s*m11+c*m21   s*m12+c*m22  ] 
	//
	double co = cos(phi);
	double si = sin(phi);
	double t = mat[1][0];
	mat[1][0] = co*t - si*mat[2][0];
	mat[2][0] = si*t + co*mat[2][0];
	t = mat[1][1];
	mat[1][1] = co*t - si*mat[2][1];
	mat[2][1] = si*t + co*mat[2][1];
	t = mat[1][2];
	mat[1][2] = co*t - si*mat[2][2];
	mat[2][2] = si*t + co*mat[2][2];
}
	
inline void Rotation::vrotate(double phi)
{
	// rotate mat around the horizontal axis 
	//
	//	[   c   0   s ]   [  m00  m01  m02 ]   [  c*m00+s*m20   c*m01+s*m21   c*m02+s*m22  ]  
	//	[   0   1   0 ] * [  m10  m11  m12 ] = [      m10           m11           m12      ]
	//	[  -s   0   c ]   [  m20  m21  m22 ]   [ -s*m00+c*m20  -s*m01+c*m21  -s*m02+c*m22  ] 
	//
	double co = cos(phi);
	double si = sin(phi);
	double t = mat[0][0];
	mat[0][0] =  co*t + si*mat[2][0];
	mat[2][0] = -si*t + co*mat[2][0];
	t = mat[0][1];
	mat[0][1] =  co*t + si*mat[2][1];
	mat[2][1] = -si*t + co*mat[2][1];
	t = mat[0][2];
	mat[0][2] =  co*t + si*mat[2][2];
	mat[2][2] = -si*t + co*mat[2][2];
}

inline void Rotation::srotate(double phi)
{
	// rotate mat around the horizontal axis 
	//
	//	[   c  -s   0 ]   [  m00  m01  m02 ]   [  c*m00-s*m10   c*m01-s*m11   c*m02-s*m12  ]  
	//	[   s   c   0 ] * [  m10  m11  m12 ] = [  s*m00+c*m10   s*m01+c*m11   s*m02+s*m12  ]
	//	[   0   0   1 ]   [  m20  m21  m22 ]   [      m20           m21           m22      ] 
	//
	double co = cos(phi);
	double si = sin(phi);
	double t = mat[0][0];
	mat[0][0] = co*t - si*mat[1][0];
	mat[1][0] = si*t + co*mat[1][0];
	t = mat[0][1];
	mat[0][1] = co*t - si*mat[1][1];
	mat[1][1] = si*t + co*mat[1][1];
	t = mat[0][2];
	mat[0][2] = co*t - si*mat[1][2];
	mat[1][2] = si*t + co*mat[1][2];
}


inline double Rotation::operator=(Spinor & S)
{
	double c0 = S.c0;
	double c1 = S.c1;
	double c2 = S.c2;
	double c3 = S.c3;
	mat[0][0] = c0*c0 + c1*c1 - c2*c2 - c3*c3;
	mat[0][1] = 2.*(c1*c2-c0*c3);
	mat[0][2] = 2.*(c1*c3+c0*c2);
	mat[1][0] = 2.*(c1*c2+c0*c3);
	mat[1][1] = c0*c0 - c1*c1 + c2*c2 - c3*c3;
	mat[1][2] = 2.*(c2*c3-c0*c1);
	mat[2][0] = 2.*(c1*c3-c0*c2);
	mat[2][1] = 2.*(c2*c3+c0*c1);
	mat[2][2] = c0*c0 - c1*c1 - c2*c2 + c3*c3;
}

inline double Rotation::print(const char * text)
{
	char * blanks = strdup(text);
	for(char * p=blanks; *p; p++) *p=' ';
	printf("%s= %10f %10f %10f\n", text,   mat[0][0], mat[0][1], mat[0][2]);
	printf("%s= %10f %10f %10f\n", blanks, mat[1][0], mat[1][1], mat[1][2]);
	printf("%s= %10f %10f %10f\n", blanks, mat[2][0], mat[2][1], mat[2][2]);
	free(blanks);
}

inline void Rotation::rotateVector(double * v)
{
	double v1 = mat[0][0]*v[0]+mat[0][1]*v[1]+mat[0][2]*v[2];
	double v2 = mat[1][0]*v[0]+mat[1][1]*v[1]+mat[1][2]*v[2];
	double v3 = mat[2][0]*v[0]+mat[2][1]*v[1]+mat[2][2]*v[2];

	v[0] = v1; v[1] = v2; v[2] = v3;
// 	printf("v= %f %f %f\n", v[0], v[1], v[2]);
}
	
inline double Rotation::det()
{
	// calculate the determinante of the matrix
	double d = mat[0][0]  * (mat[1][1]*mat[2][2] - mat[1][2]*mat[2][1] )
	         + mat[0][1]  * (mat[1][2]*mat[2][0] - mat[1][0]*mat[2][2] )
	         + mat[0][2]  * (mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0] );
	return d;
}
	
inline double Rotation::angleNaxis(double * n)
{
	// calculate the rotation axis and angle
	// there are two solutionsi: (angle and vector) and (-angle and - vector).
	// we pick the one where sin(angle) is positive
	double n0 = mat[0][1]-mat[1][0];
	double n1 = mat[0][2]-mat[2][0];
	double n2 = mat[2][1]-mat[1][2];
	double si = sqrt(n0*n0 + n1*n1 + n2*n2 );
	n[0] = n0/si;
	n[1] = n1/si;
	n[2] = n2/si;
	double co = (mat[0][0]+mat[1][1]+mat[2][2] - 1.)/2.;
	double angle = atan2(si,co);
	return angle;
}
	
//--------------------------------------------------------------------------------------------------------------------

inline void Spinor::unit()
{
	c0 =1;
	c1 = c2 = c3 = 0.;
}

inline void Spinor::rotate(double phi, double ux, double uy, double uz)
{

	//  add a rotation given by the angle phi around the axis given by the unit vector (ux,uy,uz) to this spinor
	double co = cos(phi/2.);
	double si = sin(phi/2.);

	double k0 = co;
	double k1 = si * ux;
	double k2 = si * uy;
	double k3 = si * uz;

	double x0 = k0*c0 - k3*c3 - k2*c2 - k1*c1;
	double x1 = k1*c0 + k2*c3 - k3*c2 + k0*c1;
	double x2 = k2*c0 - k1*c3 + k0*c2 + k3*c1;
	double x3 = k3*c0 + k0*c3 + k1*c2 - k2*c1;

	c0 = x0; c1 = x1; c2 = x2; c3 = x3;
}

inline void Spinor::hrotate(double phi)
{
	//  add a rotation given by the angle phi around the x axis  to this spinor
	double co = cos(phi/2.);
	double si = sin(phi/2.);
// 	printf("sin cos %f %f\n", si, co);
	double t = c0;
	c0 = co*c0 - si*c1;
	c1 = si*t  + co*c1;
	t = c2;
	c2 = co*t  - si*c3;
	c3 = co*c3 + si*t;
}

inline void Spinor::vrotate(double phi)
{
	//  add a rotation given by the angle phi around the y axis  to this spinor
	double co = cos(phi/2.);
	double si = sin(phi/2.);
// 	printf("sin cos %f %f\n", si, co);
	double t = c0;
	c0 = co*t  - si*c2;
	c2 = si*t  + co*c2;
	t = c1;
	c1 = si*c3 + co*t;
	c3 = co*c3 - si*t;
}

inline void Spinor::srotate(double phi)
{
	//  add a rotation given by the angle phi around the s (z) axis  to this spinor
	double co = cos(phi/2.);
	double si = sin(phi/2.);
// 	printf("sin cos %f %f\n", si, co);
	double t = c0;
	c0 = co*t  - si*c3;
	c3 = si*t  + co*c3;
	t = c1;
	c1 = co*t  - si*c2;
	c2 = co*c2 + si*t ;
}



inline double Spinor::angleNaxis(double * n)
{
	// calculate the rotation axis and angle
	// there are two solutionsi: (angle and vector) and (-angle and - vector).
	// we pick the one where sin(angle) is positive
	double si = sqrt( c1*c1 + c2*c2 + c3*c3 );
	n[0] = c1/si;
	n[1] = c2/si;
	n[2] = c3/si;
// 	printf("n= %f %f %f\n", n[0], n[1], n[2]);
	return 2.*atan2(si, c0);
}

inline double Spinor::det()
{
	// calculate the determinante of the matrix
	double d =  c0*c0 + c1*c1 + c2*c2 + c3*c3;
	return d;
}
	

inline double Spinor::print(const char * text)
{
	printf("%s= %f %f %f %f\n", text, c0, c1, c2, c3);
}

inline double Spinor::operator=(Rotation & S)
{
	// calculate the rotation axis
	double c022 = S.mat[0][0] + S.mat[1][1] + S.mat[2][2] - 1.;
	c0 = sqrt(c022)/2.;
	c0 = (S.mat[1][0]-S.mat[0][1])/c022;
	c1 = (S.mat[2][0]-S.mat[0][2])/c022;
	c2 = (S.mat[1][2]-S.mat[2][1])/c022;
}

#endif
