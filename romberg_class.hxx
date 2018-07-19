class Romberg
{
	double  EPS;
	int JMAX;
	int JMAXP;
	int K;
	double trap_s;
public:
	double qromb(double a, double b);
	void polint(double xa[], double ya[], int n, double x, double &y, double &dy);
	double trapzd(double a, double b, int n);
	Romberg(double _EPS, int _JMAX = 20, int _K = 5);

	virtual double func(double x) = 0;

};


