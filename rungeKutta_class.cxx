#include "rungeKutta_class.hxx"

int RungeKutta::rkqc(double & x, double htry, double eps, double & hdid, double & hnext)
{
	const double PGROW = -0.20;
	const double PSHRNK = -0.25;
	const double FCOR  = 1.0/15.0 ;
	const double SAFETY = 0.9;
	const double ERRCON = 6.0e-4;

	double xsav=x;
	for (int i=0;i<nvar;i++)
	{
		ysav[i]=y[i];
		dysav[i]=dydx[i];
	}
	double h=htry;
	for (;;) {
		double hh=0.5*h;
		rk4(ysav,dysav,xsav,hh,ytemp);
		x=xsav+hh;
		derivs(x,ytemp,dydx);
		rk4(ytemp,dydx,x,hh,y);
		x=xsav+h;
		if (x == xsav)
		{
// 			fprintf(stderr,"Step size too small in routine RKQC\n");
			return -1;
		}
		rk4(ysav,dysav,xsav,h,ytemp);
		double errmax=0.0;
		for (int i=0;i<nvar;i++)
		{
			ytemp[i]=y[i]-ytemp[i];
			double temp=fabs(ytemp[i]/yscal[i]);
			if (errmax < temp) errmax=temp;
		}
		errmax /= eps;
		if (errmax <= 1.0)
		{
			hdid=h;
			hnext=(errmax > ERRCON ?  SAFETY*h*exp(PGROW*log(errmax)) : 4.0*h);
			break;
		}
		h=SAFETY*h*exp(PSHRNK*log(errmax));
	}
	for (int i=0;i<nvar;i++) y[i] += ytemp[i]*FCOR;
	return 0;
}




void RungeKutta::rk4(double * y, double * dydx, double  x, double  h, double * yout)
{

	double hh=h*0.5;
	double h6=h/6.0;
	double xh=x+hh;
	for (int i=0;i<nvar;i++) yt[i]=y[i]+hh*dydx[i];
	derivs(xh,yt,dyt);
	for (int i=0;i<nvar;i++) yt[i]=y[i]+hh*dyt[i];
	derivs(xh,yt,dym);
	for (int i=0;i<nvar;i++)
	{
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	derivs(x+h,yt,dyt);
	for (int i=0;i<nvar;i++)
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
}









void RungeKutta::odeint(double * ystart,double x1,double x2,double eps,double h1,double hmin,int & nok,int & nbad)
{
	const int MAXSTP=10000;
	const double  TINY=1.0e-30;
	double xsav,x,hnext,hdid,h;

	x=x1;
	h=(x2 > x1) ? fabs(h1) : -fabs(h1);
	nok = nbad = kount = 0;
	for (int i=0;i<nvar;i++) y[i]=ystart[i];
	if (kmax > 0) xsav=x-dxsav*2.0;
	for (int nstp=0;nstp<MAXSTP;nstp++)
	{
		derivs(x,y,dydx);
		for (int i=0;i<nvar;i++)
			yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		if (kmax > 0)
		{
			if (fabs(x-xsav) > fabs(dxsav))
			{
				if (kount < kmax-1)
				{
					xp[++kount]=x;
					for (int i=0;i<nvar;i++) yp[i][kount]=y[i];
					xsav=x;
				}
			}
		}
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
		rkqc(x,h,eps,hdid,hnext);
		if (hdid == h) ++nok; else ++nbad;
		if ((x-x2)*(x2-x1) >= 0.0) {
			for (int i=0;i<nvar;i++) ystart[i]=y[i];
			if (kmax) {
				xp[++kount]=x;
				for (int i=0;i<nvar;i++) yp[i][kount]=y[i];
			}
			return;
		}
// 		if (fabs(hnext) <= hmin)
// 		{
// 			fprintf(stderr,"Step size too small in ODEINT\n");
// 		}
		h=hnext;
	}
	fprintf(stderr,"Too many steps in routine ODEINT\n");
}
