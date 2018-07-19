#include <stdio.h>
#include "spline_class.hxx"


main()
{
	const int n = 10;
	double x[n], y[n];

	for(int i=0; i<n; i++)
	{
		x[i] = i*0.63;
		y[i] = x[i]*x[i]*x[i];
		printf("%lf \t %lf \n",x[i],y[i]);
	}

	printf("\n\n");

// Spline::Spline(double *x,double *y,int n,double yp1, double ypn,int storage)

// 	Spline y2(x,y,n,1.e33,1.e33,0);
	Spline y2(x,y,n,0.,1.e33,0);

	for(int i=0; i<100; i++)
	{
		double xx = 0.063*i;
		double yy = y2 * xx;
		double zz = xx*xx;
		double ww = (y2 % xx)/(xx*xx);
		printf("%lf \t %lf \t %lf \t %lf\n",xx,yy, zz,ww);
	}
}
