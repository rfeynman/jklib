#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "powell_class.hxx"


static const int ITMAX=200;
static const double TOL=2.0e-4;
static const double GOLD=1.618034;
static const double GLIMIT=100.0;
static const double TINY=1.0e-20;
static const double CGOLD=0.3819660;
static const double ZEPS=1.0e-10;


int powell::match(double p[], double ftol, int &iter, double &fret)
{

	for (int i=0;i<n;i++)
		for (int j=0;j<n;j++) xi[j*n+i]= i== j ? 1: 0;
	fret=func(p);
	for (int j=0;j<n;j++) pt[j]=p[j];
	for (iter=0;iter < ITMAX;iter++) {
		double fp=fret;
		int ibig=0;
		double del=0.0;
		for (int i=0;i<n;i++) {
			for (int j=0;j<n;j++) xit[j]=xi[j*n+i];
			double fptt=fret;
			linmin(p,xit,fret);
			if (fabs(fptt-fret) > del) {
				del=fabs(fptt-fret);
				ibig=i;
			}
		}
		if (2.0*fabs(fp-fret) <= ftol*(fabs(fp)+fabs(fret))) return 0; 

		for (int j=0;j<n;j++) {
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}
		double fptt=func(ptt);
		if (fptt < fp) {
			double fft = fp-fptt;
			double ffd = fp-fret-del;
			double t=2.0*(fp-2.0*fret+fptt)*ffd*ffd - del*fft*fft;
			if (t < 0.0) {
				linmin(p,xit,fret);
				for (int j=0;j<n;j++) xi[j*n+ibig]=xit[j];
			}
		}
	}
	fprintf(stderr, "Too many iterations in routine POWELL");
	return 1;
}



void powell::linmin(double p[],double xi[],double& fret)
{
	double xx,fx,fb,fa,bx,ax,xmin;

	for(int j=0; j<n; j++) {
		pcom[j] = p[j];
		xicom[j] = xi[j];
	}

	ax=0.0; xx=1.0; bx=2.0;
	mnbrak(ax,xx,bx,fa,fx,fb);
	fret=brent(ax,xx,bx,xmin);
	for(int j=0; j<n; j++) {
		p[j] += ( xi[j] *= xmin );
	}
}




void powell::mnbrak(double & ax,double & bx,double & cx,double & fa,double & fb,double & fc)
{
	double fu,dum;

	fa=f1dim(ax);
	fb=f1dim(bx);
	if (fb > fa) {
		EXCHG(ax,bx);
		EXCHG(fb,fa);
	}
	cx=(bx)+GOLD*(bx-ax);
	fc=f1dim(cx);
	while (fb > fc) {
		double r=(bx-ax)*(fb-fc);
		double q=(bx-cx)*(fb-fa);
		double u=bx-((bx-cx)*q-(bx-ax)*r)/ (2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
		double ulim=bx+GLIMIT*(cx-bx);
		if ((bx-u)*(u-cx) > 0.0) {
			fu=f1dim(u);
			if (fu < fc) {
				ax=(bx);
				bx=u;
				fa=(fb);
				fb=fu;
				return;
			} else if (fu > fb) {
				cx=u;
				fc=fu;
				return;
			}
			u=(cx)+GOLD*(cx-bx);
			fu=f1dim(u);
		} else if ((cx-u)*(u-ulim) > 0.0) {
			fu=f1dim(u);
			if (fu < fc) {
				SHFT(bx,cx,u,cx+GOLD*(cx-bx));
				SHFT(fb,fc,fu,f1dim(u));
			}
		} else if ((u-ulim)*(ulim-cx) >= 0.0) {
			u=ulim;
			fu=f1dim(u);
		} else {
			u=(cx)+GOLD*(cx-bx);
			fu=f1dim(u);
		}
		SHFT(ax,bx,cx,u);
		SHFT(fa,fb,fc,fu);
	}
}



double powell::brent(double ax,double bx,double cx, double & xmin)
{
	int iter;
	double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;

	a=((ax < cx) ? ax : cx);
	b=((ax > cx) ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=f1dim(x);
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=TOL*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			xmin = x;
			return fx;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=f1dim(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u);
			SHFT(fv,fw,fx,fu);
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	fprintf(stderr,"Too many iterations in BRENT");
	return fx;
}




double powell::f1dim(double x)
{
	for (int j=0;j<n;j++) xt[j]=pcom[j]+x*xicom[j];
	double f=func(xt);
	return f;
}
