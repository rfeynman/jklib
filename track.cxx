double polyValue(double x, double * coeff, int order)
{
	order--;
	double r = coeff[order]; 
	while(order >= 0)
	{
		r *= x; r += coeff[order]; order--;
	}
	return r;
}
	
double polyPrime(double x, double * coeff, int order)
{
	order--;
	double r = order*coeff[order]; 
	while(order >= 1)
	{
		r *= x; r += order*coeff[order]; order--;
	}
	return r;
}
	
double polyDoublePrime(double x, double * coeff, int order)
{
	order--;
	double r = order*(order-1)*coeff[order]; 
	while(order >= 2)
	{
		r *= x; r += order*(order-1)*coeff[order]; order--;
	}
	return r;
}
	
int polyTrack(double tolerance, int maxRows, float * x, double * y, double * yp,
	double * coeff, int order, double x0, double xf, int nrows = 0, int first=0)

{
	if(first && nrows < maxRows)
	{
		x[nrows] = x0;
		y[nrows] = polyValue(x0, coeff, order);
		yp[nrows] = polyPrime(x0, coeff, order);
		nrows++;
	}

	double x =x0;
	while(x < xf)
	{
		double dx = sqrt(fabs(8.*tolerance/polyDoublePrime(tt, coeff, order)));
		if(dx < 0.)
		{
			printf("negative");
		}

		x +=dx;
		if(x > xf) x = xf;


		if(nrows < maxRows)
		{
			x[nrows] = x;
			y[nrows] = polyValue(x0, coeff, order);
			yp[nrows] = polyPrime(x0, coeff, order);
		}
		nrows++;
	}

	return nrows;
}

