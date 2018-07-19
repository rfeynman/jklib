#define SPLINE_ARRAY_USE      0
#define SPLINE_ARRAY_TAKEOVER 1
#define SPLINE_ARRAY_COPY     2

// spline function from numerical receipes with additional functionality
//
// J kewisch, 2000
//
//
// the Wave Form generators at RHIC use mirco-processors without floating point instructions
// Instead, 64-bit integer operations are used to simulate floating point operations.
// This package contains routines that run on high level computers to calculate integer parameters
// for the WFG.
//
//
class Spline
{
public:
	double * x;
	double * y;
	double * yp;
	double * cfa;
	double * cfb;
	double * cfc;
	double * cfd;
	double * integral;
	double yp1;
	double ypn;
	int xAscend;
	int yAscend;
	int n;
	int storage;
	int addl(int a , int b);
	int mulsh(int a, int b, int s);
	void normalize(double d, int & ic, int & shift);
	void alloc(double *x,double *y,double * yp, int n, int storage);
	void makeCoeff(double * u, double * y2, int from, int to, double yp1, double ypn);
	Spline();
public:
	Spline(double *x,double *y,int n,double yp1,		// constructor.
		double ypn,int storage, int special=0);    	// when special is true extra points are inserted in local minima
								// to make the integer arithmetic work
								// yp1 and ypn are the slopes at the ends. if yp1 and/or ypn
								// >= 1e30 then the y"  at the end is assumed

	Spline(double *x,double *y,double * yp,			// constructor, when slopes are defined at the data points.
			int n, int storage);
	~Spline();
	Spline(const Spline & s);
	Spline& operator=(const Spline & s);
	void setIntegral();

// 	the first version of this package used operator overloading, which was stupid
	double lookup(double x); // spline lookup, was spline*double
	double inverseLookup(double x); // inverse lookup, was spline/double
	double primeLookup(double x); // derivative lookup, was spline%double
	double doublePrimeLookup(double x); // second derivative lookup, new
	double linearLookup(double x); // linear lookup, was spline^double
	double inverseLinearLookup(double y); // inverse linear lookup, was spline|double
	double integralLookup(double x); // integral from the x[0] to x over the function.
	int wfgLookup(int * table, int xa);
	void wfgTable(int * table);
	int trackFill(double ymax, int maxRows, double * xr, double * yr, double * yp,
	         int n, double x0, double xf, int nrows = 0, int first=0);

	int track(double tolerance, int maxRows, double * xr, double * yr, double * yp,
	         int n, double x0, double xf, int nrows = 0, int first=0);

};

