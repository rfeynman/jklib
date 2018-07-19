// simplex optimization according to numerical recipes.
// modified: arrays start with index zero.
// amoeba is now an abstact class. The user derives his own class
// in order to oevrload the quality function "funk" with his own
// problems. member data of this derived class can also be used
// to pass information into funk without using global variables.

// The amoeba constuctor takes ndim as an argument and allocates 
// all necessary arrays. The function init_p can be used to initialize the
// array p and y. base[ndim] is the start point, delta[ndim] is the variation. 
// The function minimize is the original amoeba function.
// minimize  returns error code:
//			err = 0  ok
//			err > 0  no convergenve in "err" number of calls
// minimize has no arguments since all arguments of the old amoeba function
// are now public members of the class.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

class Amoeba
{
	int ndim;
	int mpts;
	int NMAX;
	double ALPHA;
	double BETA;
	double GAMMA;
	double * psum;
	double * ptry;
	double amotry(int ihi, double fac);
public:
	double **p;
	double *y;
	double ftol;
	virtual double funk(double *)=0;
	int nfunk;
	Amoeba(int ndim);
	virtual ~Amoeba();
	void init_p(double * base, double * delta, double tolerance);
	int minimize();
};

