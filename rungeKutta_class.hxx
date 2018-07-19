#include <math.h>
#include <stdlib.h>
#include <stdio.h>



class RungeKutta
{
	const int nvar;
        double * dysav;
        double * ysav;
        double * ytemp;
	double * dym;
	double * dyt;
	double * yt;
	double * yscal;
	double * dydx;
	double * y;

	const double dxsav;	// step size for result storage
	const int kmax;	// size of result storage




	void rkalloc(int n)
	{
        	dysav= new double[n];
        	ysav=  new double[n];
        	ytemp= new double[n];
		dym=   new double[n];
		dyt=   new double[n];
		yt=    new double[n];
		yscal= new double[n];
		dydx=  new double[n];
		y=     new double[n];
		xp=0;yp=0;
	}

public:
	int kount;	// after odint, kount containes the number of values in xp
	double *xp;	// contains the x values for saved results
	double **yp;	// contains the y[_n] values for saved results


	RungeKutta()  : nvar(0), kmax(0), dxsav(0.)
	{
		fprintf(stderr, "dont use the RungeKutta default constructor!!");
		exit(-1);
	}


	// _n = number of differentail equations
	// _dxsav = save result every _dssav step
	// _kmax = max size of the save array
	RungeKutta(const int _n, const int _kmax, const double _dxsav) : nvar(_n), kmax(_kmax), dxsav(_dxsav)
	{
		rkalloc(_n);
		xp=     new double[kmax];
		yp=     new double*[_n];
		for(int i=0; i< _n; i++) yp[i]=     new double[kmax];
	}

	RungeKutta(const int _n) : nvar(_n), kmax(0), dxsav(0.)
	{
		rkalloc(_n);
		xp=0;
		yp=0;
	}

	~RungeKutta()
	{
		delete [] dysav;
		delete [] ysav;
		delete [] ytemp;
		delete [] yt;
		delete [] dyt;
		delete [] dym;
		delete [] yscal;
		delete [] dydx;
		delete [] y;
		delete [] xp;
		if(yp) for(int i=0; i< nvar; i++) delete [] yp[i];
		delete [] yp;
	}

	// runge-kutta driver with adaptive step size control. Integtrate starting values ystart[0..n-1] from x1 to
	// x2 with accuracy eps, storing immediate results in ...
	// h1 is the guessed initial step size, hmin is the minimum step size (can be zero).
	// nok and nbad are the number of good and bad (but fixed) steps
	void odeint(double * ystart,double x1,double x2,double eps,double h1,double hmin,int & nok,int & nbad);


	int rkqc(	double & x,	// initial independent variable ( e.g. time)
			double htry,	// initial stepsize in time
			double eps,
			double & hdid,
			double & hnext
         	);
	void rk4( 	double * y,	// initial function values
			double * dydx,  // initial derivatives calculated by derivs()
			double  x,	// initial independent variable ( e.g. time)
			double  h,	// initial stepsize in time
			double * yout);	// new y after this step
	// user function:
	// calculate the derivatives of a first order coupled  system of differential equations
	virtual void derivs(double x,		// independent variable ( e.g. time)
			    double * y,  	// function values at this time
			    double * dydx	// result dy/dx = f(y, x) 
			   )=0;
};	

