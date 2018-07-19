#ifndef __RampFit__
#define __RampFit__

class RampFit
{
  void stromt();
  
  RampFit(); 
  void strom(double up, double i_in, double ip_in,
	     double t1, double& i_out, double& ip_out);
public:
  double l;
  double r;
  double o;
  double upmax;
  double up[3];
  double current[4];
  double current_p[4];
  double dt[4];
  double goal;
  double goal_p;
  double umax;
  double imax;
  double imin;
  RampFit(double l, double r, double umax, double upmax,double imax, double imin);
  RampFit(double * vec, double c0, double cp0);
  int minimize(double t);
  int minimize3(double, double, double);
  int shortRamp( int debug = 0 );
  void print_result(const char * what);
  int track(double tolerance, int maxRows, float * x, double * y, int nrows = 0, int first=0);
  int trackp(double tolerance, int maxRows, float * x, double * y, double * yp, int nrows = 0, int first=0);
  int trackpp(double tolerance, int maxRows, float * x, double * y, double * yp, double * ypp, int nrows = 0, int first=0);
  double cur2dt(double);
  void dt2cur(double tt, double &i1, double &ip1);
  double dt2cur(double tt);
  void set_ramp(double i0, double ip0, double i1, double ip1);
  double * dT();
};

#endif
