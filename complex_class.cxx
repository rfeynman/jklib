#include "complex_class.hxx"


Complex conjugate(const Complex &x)
{
	Complex a(x.r,-x.i);
	return a;
}

double abs(const Complex &x)
{
	return sqrt( x.r*x.r + x.i*x.i);
}

Complex sqrt(const Complex &x)
{
	double phi = atan2(x.i,x.r)/2.;
	double a = pow(x.i*x.i + x.r*x.r, 0.25);
	Complex r(a*cos(phi), a*sin(phi));
	return r;
}

Complex operator+(const Complex &x, double y)
{
	Complex r(x.r + y, x.i);
	return r;
}

Complex operator+(double y, const Complex &x)
{
	Complex r(x.r + y, x.i);
	return r;
}

Complex operator-(const Complex &x, double y)
{
	Complex r(x.r - y, x.i);
	return r;
}

Complex operator-(double y, const Complex &x)
{
	Complex r(x.r - y, x.i);
	return r;
}

Complex operator*(const Complex &x, double y)
{
	Complex r(x.r * y, x.i * y);
	return r;
}

Complex operator*(double y, const Complex &x)
{
	Complex r(x.r * y, x.i * y);
	return r;
}

Complex operator/(const Complex &x, double y)
{
	Complex r(x.r / y, x.i / y);
	return r;
}


Complex operator+(const Complex &x, const Complex &y)
{
	Complex r(x.r + y.r, x.i + y.i);
	return r;
}

Complex operator-(const Complex &x, const Complex &y)
{
	Complex r(x.r - y.r, x.i - y.i);
	return r;
}

Complex operator*(const Complex &x, const Complex &y)
{
	Complex r(x.r * y.r - x.i * y.i , x.i * y.r + x.r * y.i);
	return r;
}

