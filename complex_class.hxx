#ifndef _complex_class_
#define _complex_class_

#include <math.h>

class Complex
{
public:
	double r;
	double i;

	Complex(double re=0., double im=0.) { r=re; i=im; }
	Complex operator-() { Complex x(-r, -i); return x;}

};

double abs(const Complex &x);

Complex conjugate(const Complex &x);

Complex sqrt(const Complex &x);

Complex operator+(const Complex &x, double y);

Complex operator+(double y, const Complex &x);

Complex operator-(const Complex &x, double y);

Complex operator-(double y, const Complex &x);

Complex operator*(const Complex &x, double y);

Complex operator*(double y, const Complex &x);

Complex operator/(const Complex &x, double y);

Complex operator+(const Complex &x, const Complex &y);

Complex operator-(const Complex &x, const Complex &y);

Complex operator*(const Complex &x, const Complex &y);

#endif
