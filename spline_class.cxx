#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "spline_class.hxx"

// this version has been modified to let arrays start with the index [0]

typedef long long int64;
// typedef long int64;  // !!!! changed to make tcl work, must be changed back !!!!

Spline::Spline()
{
	printf("spline empty constructor...dont use that\n");
}

Spline::~Spline()
{
	if(storage == SPLINE_ARRAY_TAKEOVER || storage == SPLINE_ARRAY_COPY)
	{
		delete [] x;
		delete [] y;
		delete [] yp;
	}
	delete [] cfa;
	delete [] cfb;
	delete [] integral;
	if(!yp) delete [] cfc;
}

Spline::Spline(const Spline & s)
{
	printf("spline copy constructor...dont use that\n");
	exit(-1);
}
Spline& Spline::operator=(const Spline & s)
{
	printf("spline operator=...dont use that\n");
	exit(-1);
}


Spline::Spline(double *x,double *y,double * yp, int n, int storage)
{
// cubic interpolation, known y0,y1,yp0,y'1

	alloc(x, y, yp, n, storage);
	cfc = yp;
	cfd = y;
	for(int i=0; i< n-1; i++)
	{
		double h = x[i+1]-x[i];
		double h2 = h*h;
		double h3 = h*h2;

		cfa[i] = -2.*(y[i+1]-y[i])/h3 +  (yp[i+1]+yp[i])/h2;
		cfb[i] =  3.*(y[i+1]-y[i])/h2 - (yp[i+1]+2.*yp[i])/h;
		integral[i] = (((cfa[i+1]/4.*h + cfb[i+1]/3.)*h + cfc[i+1]/2.)*h+cfd[i+1])*h;
	}
	cfa[n-1] =0.;
	cfb[n-1] =0.;
}


Spline::Spline(double *x,double *y,int n,double yp1, double ypn,int storage, int special)
{
	alloc(x, y, NULL, n, storage);
	this->yp1 = yp1;
	this->ypn = ypn;
	double * y2 = new double[n];
	if(!y2) { printf("out of memory\n"); exit(-1); }
	double *u  = new double[n];
	if(!u) { printf("out of memory\n"); exit(-1); }
	if(special)
	{
		int from = 0;
		double yp11 = yp1;
		for(int i=from+1; i<n-1; i++)
		{
			if(    (y[i] == y[i-1] || y[i] == y[i+1])
			    || (y[i] >= y[i-1] && y[i] >= y[i+1])
			    || (y[i] >= y[i-1] && y[i] >= y[i+1]) )
			{
				//local minimum found
				makeCoeff(u, y2, from, i, yp11, 0.);
				yp11 = 0.;
				from = i;
			}
		}
		makeCoeff(u, y2, from, n-1, 0., ypn);
	}
	else
	{
		makeCoeff(u, y2, 0, n-1, yp1, ypn);
	}
	delete [] u;
	delete [] y2;
}


void Spline::setIntegral()
{

	integral[0] = 0.;
	for(int i=1; i< n; i++)
	{
		double h = x[i]-x[i-1];
		integral[i] = integral[i-1] + (((cfa[i]/4.*h + cfb[i]/3.)*h + cfc[i]/2.)*h+cfd[i])*h;
	}
}


void Spline::makeCoeff(double * u, double * y2, int from, int to, double yp1, double ypn)
{
	for (int i=from;i<to;i++) {
		if(x[i+1] <= x[i])
		{
			fprintf(stderr, "spline class: x is not increasing, x[%d]=%e, x[%d]=%e\n", i,  x[i], i+1, x[i+1]);
			exit(-1);
		}
	}

	if (yp1 > 0.99e30)
		y2[from]=u[from]=0.0;
	else {
		y2[from] = -0.5;
		u[from]=(3.0/(x[from+1]-x[from]))*((y[from+1]-y[from])/(x[from+1]-x[from])-yp1);
	}
	for (int i=from+1;i<to;i++) {
		double sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		double p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}

	double qn,un;
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[to]-x[to-1]))*(ypn-(y[to]-y[to-1])/(x[to]-x[to-1]));
	}
	y2[to]=(un-qn*u[to-1])/(qn*y2[to-1]+1.0);
	for (int k=to-1;k>=from;k--) y2[k]=y2[k]*y2[k+1]+u[k];
// end of finding y2. now calc the coefficients

        for (int i=from;i<to;i++) {
                double dx = x[i+1] - x[i]; 
                double dy = y[i+1] - y[i]; 
                double d2y = y2[i+1] - y2[i];
                cfd[i] = y[i]; 
                cfc[i] = dy/dx - dx*y2[i]/2. - dx*d2y/6.;
                cfb[i] = y2[i]/2.;
                cfa[i] = d2y/(dx*6);
//              printf("coeff a..d %d %15e %15e %15e %15e\n",
//                      i, cfa[i], cfb[i], cfc[i], cfd[i]);
        }
	cfa[to] = 0.;
	cfb[to] = 0.;
	cfc[to] = ypn;
	cfd[to] = y[to];
}


void Spline::normalize(double value, int & mantissa, int & shift)
{
	double frac = frexp(value, &shift)*0x3fffffff; shift -= 30;
        if( frac < 0.) frac -= 0.5; else frac += 0.5;
        mantissa = int(frac);
        if(!mantissa) shift=0;
        printf("normalize: value = %e mantissa = %d = %x shift = %d\n",
                value, mantissa ,mantissa, shift);
}



void Spline::alloc(double *x,double *y,double * yp, int n, int storage)
{
	if(storage == SPLINE_ARRAY_USE || storage == SPLINE_ARRAY_TAKEOVER)
	{
		this->x = x;
		this->y = y;
		this->yp = yp;
	}
	else if(storage == SPLINE_ARRAY_COPY)
	{
		this->x = new double[n];
		if(!this->x) { printf("out of memory\n"); exit(-1); }
		memcpy(this->x,x,n*sizeof(double));

		this->y = new double[n];
		if(!this->y) { printf("out of memory\n"); exit(-1); }
		memcpy(this->y,y,n*sizeof(double));

		if(yp)
		{
		this->yp = new double[n];
		if(!this->yp) { printf("out of memory\n"); exit(-1); }
		memcpy(this->yp,yp,n*sizeof(double));
		}
		else this->yp = NULL;
		yp = this->yp;
	}
	else
	{
		printf(" \"storage\" argument is out of range\n");
		exit(-1);
	}
	this->n = n;
	this->storage = storage;
	cfa = new double[n];
	cfb = new double[n];
	integral = new double[n];
	if(yp) cfc = yp; else cfc = new double[n];
	cfd = y;

	if(x[0] > x[1])
	{
		xAscend = -1;
		for(int i=2; i<n; i++)
		{
			if(x[i-1] < x[i])
			{
				xAscend = 0;
				break;
			}
		}
	}
	else
	if(x[0] < x[1])
	{
		xAscend = 1;
		for(int i=2; i<n; i++)
		{
			if(x[i-1] > x[i])
			{
				xAscend = 0;
				break;
			}
		}
	}
	else
	xAscend = 0;


	if(y[0] > y[1])
	{
		yAscend = -1;
		for(int i=2; i<n; i++)
		{
			if(y[i-1] < y[i])
			{
				yAscend = 0;
				break;
			}
		}
	}
	else
	if(y[0] < y[1])
	{
		yAscend = 1;
		for(int i=2; i<n; i++)
		{
			if(y[i-1] > y[i])
			{
				yAscend = 0;
				break;
			}
		}
	}
	else
	yAscend = 0;

}




void Spline::wfgTable(int * table)
{

// insert new step stones, if necessary. inserted if:
//
// if there is a local extreme:
//	a*x*x*x + b*x*x + c*x + d = ext
// 	3*a*x*x + 2*b*x + c = 0
// 	x = -b/3*a +- sqrt( b*b/(9*a*a) - c/3*a)
// or
// if there is a local extreme in the inetrmediate result:
//	a*x*x + b*x = ext
//	2*a*x + b  = 0
//	x = -b/2*a (hey, we have done that already)


  const int maxn = 42;
  double polyCff[maxn][4];
  for (int i=0;i<n;i++)
    {
      polyCff[i][0]=cfa[i];
      polyCff[i][1]=cfb[i];
      polyCff[i][2]=cfc[i];
      polyCff[i][3]=cfd[i];
    }
  
  // wnat split?
  for (int i=0;i<n-1;i++)
    {
      double dx = x[i+1] - x[i];
      double dxx = 1.;
      for(int j=3; j>=0; j--)
	{
	  if( fabs( polyCff[i][j] * dxx) < 1.) polyCff[i][j] = 0.;
	  dxx *= dx;
	}
      
      double xmin = dx;
      int split = 0;
      if(polyCff[i][0] != 0. && polyCff[i][1] != 0.)
	{
	  double x3 = - polyCff[i][1]/(2.*polyCff[i][0]);
	  if(x3 >= 1. && x3 <= (xmin - 1. ) )
	    {
	      xmin  = x3;
	      split = 3;
	    }
	}
      
      
      double r1,r2,r3,r4,r5;
      r1 = - polyCff[i][1]/(3.*polyCff[i][0]);
      r2 = r1*r1 - polyCff[i][2]/(3.*polyCff[i][0]);
      if(r2 > 0.)
	{
	  r3 = sqrt(r2);
	  r4 = r1 + r3;
	  if(r4 >= 1. && r4 <= (xmin-1.) ) { xmin =r4; split=4;}
	  r5 = r1 - r3;
	  if(r5 >= 1. && r5 <= (xmin-1.) ) { xmin =r5; split=5;}
	}
      else
	{
	  r4 = 0.;
	  r5 = 0.;
	}
      
      if(split)
	{
	  if(n >= maxn) 
	    {
	      printf("array overflow\n");
	      exit(1);
	    }
	  
	  for(int ii = n; ii > i+1; ii--)
	    {
	      int imi = ii-1;
	      x[ii] = x[imi];
	      y[ii] = y[imi];
	      for(int k=0; k<4; k++) polyCff[ii][k]= polyCff[imi][k];
	    }
	  n++;
	  
	  x[i+1] = xmin + x[i];
	  y[i+1] = (( polyCff[i][0]*xmin+polyCff[i][1])*xmin+polyCff[i][2])*xmin+polyCff[i][3];
	  polyCff[i+1][0] = polyCff[i][0];
	  polyCff[i+1][1] = 3.*polyCff[i][0]*xmin + polyCff[i][1];
	  polyCff[i+1][2] = 3.*polyCff[i][0]*xmin*xmin + 2.*polyCff[i][1]* xmin + polyCff[i][2];
	  polyCff[i+1][3] = y[i+1];
	  printf("split %d new step stone: %f %f \n", split,x[i+1],y[i+1]);
	  printf("coeff %d %15e %15e %15e %15e\n",
		 i, polyCff[i][3], polyCff[i][2], polyCff[i][1], polyCff[i][0]);
	  printf("coeff %d %15e %15e %15e %15e\n",
		 i+1, polyCff[i+1][3], polyCff[i+1][2], polyCff[i+1][1], polyCff[i+1][0]);
	  
	}
      
    }
  
  
  // now calculate the integer coefficients
  
  
  
  for (int i=0;i<n-1;i++)
    {
      int j, k;
      double dx[4]; dx[0]=dx[1]=dx[2] = x[i+1] - x[i]; dx[3]=1.;
      int    mantix;		//  mantissa of delta x
      int    shiftx;  	//  shift of delta x
      int    mantiz[4]; 	//  mantissa of the coefficient/result mantissa
      int    shiftz[4];	//  result shift
      int    shiftc[4]; 	//  shift of the coefficient
      double * polyCoff = polyCff[i];
      printf("\n\n\ncoeff %d %15e %15e %15e %15e x=%e\n",
	     i, polyCoff[3], polyCoff[2], polyCoff[1], polyCoff[0], dx[0]);

      
      // positive shift is a right shift
      
      for(j=0; j<4; j++)
	{
	  printf("Coeff %d: ",j);
	  normalize(polyCoff[j], mantiz[j], shiftc[j]);
	  shiftz[j]=0;
	}
      printf("dx: ");
      normalize(dx[0], mantix, shiftx);
      
      
      // j is the index of the coefficient, k is the index of the shift.
      for(j=0; j<4; j++) if(shiftc[j]) break;
      if(j == 4) continue; // everything zero in this interval, may happen!
      
      double yatend = polyCoff[j] * dx[j]; // at the end of the interval
      int isshift = shiftc[j];
      
      // now look at the next coefficient, but use the old shift
      k = j;
      for(j=j+1; j<4; j++)
	{
	  
// 	printf("\nnorm resul %d:\n",j);
			int  tmpmantissa, tmpshift;
			double val1 = fabs(yatend);
			double valmax = val1;
			double val2 = fabs(polyCoff[j]);
			if( val2 > valmax) valmax = val2;
			double val3 = fabs(yatend + polyCoff[j]);
			if( val3 > valmax) valmax = val3;
			normalize(valmax, tmpmantissa, tmpshift);
			printf("valmax,1,2,3 = %e %e %e %e\n", valmax, val1, val2, val3);
// tmtshift is the right shift for the result
			if(tmpshift > 0) tmpshift=0;  // dont shift more that ints
			shiftz[k] = tmpshift - isshift;
			printf("tmpshift,isshift,shiftz[k]= %d %d %d\n", tmpshift,isshift,shiftz[k]);
// shiftz[k] is the minimum shift to keep the result within an int.


			if(mantiz[j] == 0) shiftc[j] = tmpshift;
			while(tmpshift > shiftc[j])	// create new b with enuf shift
			{
				printf("tmpshift,shiftc[j],mantiz[j] = %d %d %d\n", tmpshift,shiftc[j],mantiz[j]);
				mantiz[j] >>= 1;
				shiftc[j]++;
			}
			isshift = tmpshift;


			yatend = (yatend +polyCoff[j])*dx[j];
			k = j;
		}

		shiftz[3]= - isshift; 
		printf("%12d %12d %12d %12d %12d %12d %12d %12d %12d\n",
			int(x[i]), mantiz[0], mantiz[1], mantiz[2], mantiz[3],
			shiftz[0], shiftz[1], shiftz[2], shiftz[3]);

	int idx = int(dx[0]);
	idx = 0;
	int r = mulsh(mantiz[0],idx,shiftz[0]);
	r = addl(r,mantiz[1]);
	r = mulsh(r,idx,shiftz[1]);
	r = addl(r,mantiz[2]);
	r = mulsh(r,idx,shiftz[2]);
	r = addl(r,mantiz[3]);
	r = mulsh(r,1,shiftz[3]);
	printf("a %e\n", polyCoff[0]*dx[0]*dx[0]*dx[0]);
	printf("b %e\n", polyCoff[1]*dx[0]*dx[0]);
	printf("c %e\n", polyCoff[2]*dx[0]);
	printf("d %e\n", polyCoff[3]);




// make the table...

		table[6*i+0] = x[i] < 0 ? int(x[i]-0.5) : int(x[i]+0.5);
// 		table[6*i+0] = int(x[i]);
		table[6*i+1] = mantiz[3];  // was mantiz[0]
		table[6*i+2] = mantiz[2];
		table[6*i+3] = mantiz[1];
		table[6*i+4] = mantiz[0];
		int kk = shiftz[0];  // was shiftz[3]
		kk <<=8; kk |= shiftz[1];
		kk <<=8; kk |= shiftz[2];
		kk <<=8; kk |= shiftz[3];
		table[6*i+5] = kk;
   	}
  int ix = 6*(n-1);
  table[ix +0] = x[n-1] < 0 ? int(x[n-1]-0.5) : int(x[n-1]+0.5);
  table[ix +1] = int(polyCff[n-1][3]);
  table[ix +2] = 0;
  table[ix +3] = 0;
  table[ix +4] = 0;
  table[ix +5] = 0;
}



int Spline::mulsh(int a, int b, int s)
{
	int  sig = a >> 31;
	if( b < 0) sig = ! sig;
	int64 aa = abs(a);
	int64 bb = abs(b);
	int64 cc = aa*bb;
	if(s >= 0)   cc >>= s; else  cc <<= (-s);
	int64 dd = cc & 0xffffffff00000000LL;
	if(dd) printf("mulsh cc = %llx\n", dd);
	if(sig) cc =  - cc;
	int c = cc;
// 	printf("mulsh %d * %d >> %d = %d\n", a,b,s,c);
	return c;
}

int Spline::addl(int a , int b)
{
	int64 aa = a;
	int64 bb = b;
	int64 cc = aa+bb;
	int64 dd = labs(cc) & 0xffffffff00000000LL;
	if(dd) printf("mulsh cc = %llx\n", dd);
	int c = cc;
// 	printf("addl %d + %d = %d\n", a,b,c);
	return c;
}
 

int Spline::wfgLookup(int * table, int xa)
{
	int klo=0;
	int khi=n-1;
	while (khi-klo > 1) {
		int k=(khi+klo) >> 1;
		if (table[6*k] > xa) khi=k;
		else klo=k;
	}
	int xv = table[6*klo+0];
	int fa = table[6*klo+4];  // was table[6*klo+1]
	int fb = table[6*klo+3];
	int fc = table[6*klo+2];
	int fd = table[6*klo+1];
	int fs = table[6*klo+5];
	int ssd = char(fs & 0xff); fs >>= 8;
	int ssc = char(fs & 0xff); fs >>= 8;
	int ssb = char(fs & 0xff); fs >>= 8;
	int ssa = char(fs & 0xff); fs >>= 8;
	int dx = xa-xv;

// 	printf("%d %12d %12d %12d %12d %12d %12d %12d %12d %12d\n",klo,dx,fa,fb,fc,fd,ssa,ssb,ssc,ssd);
	int r = mulsh(fa,dx,ssa);
	r = addl(r,fb);
	r = mulsh(r,dx,ssb);
	r = addl(r,fc);
	r = mulsh(r,dx,ssc);
	r = addl(r,fd);
	r = mulsh(r,1,ssd);
	return r;
}



double Spline::integralLookup(double xa) // integral from x[0] to x
{
// 	printf("spline n=%d  xa = %f asc = %d\n",n, xa, xAscend);
	int klo=0;
	int khi=n-1;
	while (khi-klo > 1)
	{
		int k=(khi+klo) >> 1;
		if ( (x[k] - xa)*xAscend > 0) khi=k; else klo=k;
	}
// 	printf("k %d %d x %f %f y %f %f\n", klo,khi, x[klo],x[khi],y[klo],y[khi]);
	double b=xa-x[klo];
	double r= (((cfa[klo]/4.*b + cfb[klo]/3.)*b + cfc[klo]/2.)*b+cfd[klo])*b;
// 	printf("a,b,r,h %f %f %f %f\n", a,b,r,h);
	return r+integral[klo];
}

double Spline::lookup(double xa) // spline lookup
{
	int klo=0;
	int khi=n-1;
	while (khi-klo > 1)
	{
		int k=(khi+klo) >> 1;
		if ( (x[k] - xa)*xAscend > 0) khi=k; else klo=k;
	}
	double b=xa-x[klo];
// 	printf("%d x=%e c=%e %e %e %e\n",klo,b,cfa[klo],cfb[klo],cfc[klo],cfd[klo]); 
	double r= ((cfa[klo]*b + cfb[klo])*b + cfc[klo])*b+cfd[klo];
	return r;
}


double Spline::inverseLookup(double ya) // inverse lookup
{
	if(! yAscend) return 0;

	int klo=0;
	int khi=n-1;
	while (khi-klo > 1)
	{
		int k=(khi+klo) >> 1;
		if ((y[k] - ya)*yAscend > 0) khi=k; else klo=k;
	}
	double h=x[khi]-x[klo];
	double v=y[khi]-y[klo];
	double xa = (ya - y[klo])*h/v + x[klo];

	for(int i=0; i<3; i++)
	{
	double b=xa-x[klo];
	double y1=  ((cfa[klo]*b + cfb[klo])*b + cfc[klo])*b+cfd[klo];
	double yprime = ( 3.*cfa[klo]*b + 2.*cfb[klo])*b + cfc[klo];
	xa += ( ya -y1 )/yprime;
	}
	return xa;
}



double Spline::inverseLinearLookup(double ya) // inverse linear lookup
{
	if(! yAscend) return 0;

	int klo=0;
	int khi=n-1;
	while (khi-klo > 1)
	{
		int k=(khi+klo) >> 1;
		if ((y[k] - ya)*yAscend > 0) khi=k; else klo=k;
	}
	double h=x[khi]-x[klo];
	double v=y[khi]-y[klo];
	double xa = (ya - y[klo])*h/v + x[klo];

	return xa;
}


double Spline::doublePrimeLookup(double xa) // derivative lookup
{
	int klo=0;
	int khi=n-1;
	while (khi-klo > 1)
	{
		int k=(khi+klo) >> 1;
		if ( (x[k] - xa)*xAscend > 0) khi=k; else klo=k;
	}
// 	printf("k %d %d x %f %f y %f %f\n", klo,khi, x[klo],x[khi],y[klo],y[khi]);
	double b=xa-x[klo];
	double ydprime =  6.*cfa[klo]*b + 2.*cfb[klo];
	return ydprime;
}

double Spline::primeLookup(double xa) // derivative lookup
{
// 	printf("spline n=%d  xa = %f asc = %d\n",n, xa, xAscend);
	int klo=0;
	int khi=n-1;
	while (khi-klo > 1)
	{
		int k=(khi+klo) >> 1;
		if ( (x[k] - xa)*xAscend > 0) khi=k; else klo=k;
	}
// 	printf("k %d %d x %f %f y %f %f\n", klo,khi, x[klo],x[khi],y[klo],y[khi]);
	double b=xa-x[klo];
	double yprime = ( 3.*cfa[klo]*b + 2.*cfb[klo])*b + cfc[klo];
// 	printf("a,b,r,h %f %f %f %f\n", a,b,r,h);
	return yprime;
}
double Spline::linearLookup(double xa) // linear lookup
{
// 	printf("spline n=%d  xa = %f asc = %d\n",n, xa, xAscend);
	int klo=0;
	int khi=n-1;
	while (khi-klo > 1)
	{
		int k=(khi+klo) >> 1;
		if ( (x[k] - xa)*xAscend > 0) khi=k; else klo=k;
	}
// 	printf("k %d %d x %f %f y %f %f\n", klo,khi, x[klo],x[khi],y[klo],y[khi]);
	double h=x[khi]-x[klo];
	double b=xa-x[klo];
	double r = y[klo] + (y[khi]-y[klo])*b/h;
// 	printf("a,b,r,h %f %f %f %f\n", a,b,r,h);
	return r;
}





	
int Spline::trackFill(double ymax, int maxRows, double * xr, double * yr, double * yp,
	int n, double x0, double xf, int nrows, int first)
{
  double tolerance = 2.e-4*ymax;
  int nrows1 = nrows;
  for(int j=0; j<10; j++)
    {
      nrows = track(tolerance, maxRows, xr, yr, yp, n, x0, xf, nrows1, first);
      if(n < 3)  break;// all linear
      printf(" tol = %f, # of rows = %d\n", tolerance, nrows);
      if(nrows <= maxRows && nrows > maxRows*0.9) break;
      double f= double(maxRows*0.95)/double(nrows);
      tolerance /= f*f;
//       if(tolerance < 1.e-5*ymax) break;
    }
    return nrows;
}




int Spline::track(double tolerance, int maxRows, double * xr, double * yr, double * yp,
	int n, double x0, double xf, int nrows, int first)
{
	if(first && nrows < maxRows)
	{
		xr[nrows] = x0;
		yr[nrows] = lookup(x0);
		if(yp) yp[nrows] = primeLookup(x0);
		nrows++;
	}

	if(n < 3) // all linear
	{
		xr[nrows] = xf;
		yr[nrows] = lookup(xf);
		if(yp) yp[nrows] = primeLookup(xf);
		nrows++;
		return nrows;
	}

	double xa =x0;
	while(xa < xf)
	{
		double dp = doublePrimeLookup(xa);
		double dx;
		if(dp == 0.) dx = (xf - x0)/maxRows;
		else dx = sqrt(fabs(8.*tolerance/dp));
		printf("dp = %e dx = %e tol = %e\n", dp,dx,tolerance);
		xa +=dx;
		if(xa > xf) xa = xf;


		if(nrows < maxRows)
		{
			xr[nrows]  = xa;
			yr[nrows] = lookup(xa);
			if(yp) yp[nrows] = primeLookup(xa);
		}
		nrows++;
	}

	return nrows;
}

