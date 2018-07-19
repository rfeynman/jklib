class powell
{
	const int n;
	double * pt;
	double * ptt;
	double * xit;
	double * xi;
	double * xt;
	double *pcom,*xicom;
	virtual double func(double *) = 0;
	inline double MAX(double a,double b) {return a > b ? a : b;}
	inline double SIGN(double a,double b) {return b > 0.0 ? fabs(a) : -fabs(a);}
	inline void SHFT(double &a,double &b,double &c,double d) {a=b;b=c;c=d;}
	inline void EXCHG(double &a,double &b) {double t=a;a=b;b=t;}
	powell() : n(0) {};
public:
	powell(int n) : n(n)
	{
		pt= new double[n];
		ptt= new double[n];
		xit= new double[n];
		xi= new double[n*n];
		xt= new double[n];
		pcom= new double[n];
		xicom= new double[n];
	}
	~powell()
	{
		delete [] pt;
		delete [] ptt;
		delete [] xit;
		delete [] xi;
		delete [] xt;
		delete [] pcom;
		delete [] xicom;
	}
	void mnbrak(double & ax,double & bx,double & cx,double & fa,double & fb,double & fc);
	double f1dim(double x);
	int match(double p[], double ftol, int &iter, double &fret);
	double brent(double ax,double bx,double cx, double & xmin);
	void linmin(double p[],double xi[],double& fret);

};
