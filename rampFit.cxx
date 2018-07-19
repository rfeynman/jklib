#include "rampFit.hxx"
#include <math.h>
#include <stdio.h>


// l = ?, r = ?, 
RampFit::RampFit(double l, double r, double umax, double upmax, double imax, double imin) 
{
  this->l = l;
  this->r = r;
  o = r/l;
  this->umax = umax * 0.7; // allow only 90% for safty margin
  imax = fabs(imax);
  imin = fabs(imin);
  if(imin > imax) imax = imin;
  if(imax < 1) imax = 1.;
  this->imax = imax;
  this->imin = imin;
  this->upmax = fabs(upmax);
};

double * RampFit::dT()
{
  double * ddt = new double[4];
  for(int i=0; i< 4; i++) ddt[i]=dt[i];
  return ddt;
}

// vec is result[] array
// it is called like this: RampFit f(result[step], pseudo[0], pseudoP[0]);
RampFit::RampFit(double * vec, double c0, double cp0)
{
  // from SDDS file: T1, T2, T3, Ts = times for 3 regions + their sum
  for(int i=0; i<4; i++) dt[i] = vec[i];
  // from SDDS file: Resistance
  r = vec[4];
  // from SDDS file: Inductance
  l = vec[5]; o = r/l;
  // from SDDS file: Uprime and Uprime2
  up[0] = vec[6]; up[1] = 0.; up[2] = vec[9];
  upmax = fabs(vec[6]);
  // from SDDS file: Imax
  imax = vec[7];
  current[0] = c0;
  current_p[0] = cp0;
  
  for(int j=0; j<3; j++)
    {
      int j1 = j+1;
      double ur = up[j]/r;
      double uro = ur/o;
      double ex0 = exp(-o*dt[j]);
      double ex1 = 1. - ex0;
      current[j1] = current[j] + current_p[j]/o * ex1
	+ uro * (o*dt[j]-ex1);
      current_p[j1] = current_p[j]*ex0 + ur*ex1;
    }
  goal = current[3];
  goal_p = current_p[3];
}

void RampFit::set_ramp(double i0, double ip0, double i1, double ip1)
{
  current[0] = i0;
  current_p[0] = ip0;
  goal = i1;
  goal_p = ip1;
}

void RampFit::print_result(const char * what)
{
  printf("RampFit results: %s\n",what);
  double voltage[4];
  for(int i=0; i<4; i++)
    {
      voltage[i] = r * current[i] + l * current_p[i];
      printf("t%d,i%d,ip%di,u%d = %f\t%f\t%f\t%f\n",
	     i,i,i,i, dt[i], current[i], current_p[i], voltage[i]);
    }
  double vgoal = r * goal + l * goal_p;
  printf("t%d,i%d,ip%di,u%d = %f\t%f\t%f\t%f\n",
	 5,5,5,5,0.,goal,goal_p,vgoal);
}

void RampFit::strom(double up, double i_in, double ip_in,
		double t1, double& i_out, double& ip_out)
{
  double ur = up/r;
  double uro = ur/o;
  double ex1 = exp(-o*t1);
  // printf(" in = %lf %lf\n", i_in, ip_in);
  i_out = i_in + ip_in/o * (1.-ex1) + uro * (o * t1 - 1. + ex1);
  ip_out = ip_in * ex1 + ur * (1. - ex1);
  // printf("out = %lf %lf\n", i_out, ip_out);
}

double RampFit::dt2cur(double tt)
{
  double i1;
  for(int i=0; i<3; i++)
    {
      if(tt <= dt[i])
	{
	  i1 = current[i] + current_p[i]*tt
	    + (up[i]/r -current_p[i]) * o * tt *tt /2.;
	  return i1;
	}
      else tt -= dt[i];
    }
  return current[3];
}

void RampFit::dt2cur(double tt, double &i1, double &ip1)
{
  i1 = current[3];
  ip1= current_p[3];
  for(int i=0; i<3; i++)
    {
      if(tt <= dt[i])
	{
	  strom( up[i], current[i],   current_p[i], tt, i1, ip1);
	  break;
	}
      else tt -= dt[i];
    }
}

double RampFit::cur2dt(double cu)
{
  double tt=0.;
  for(int i=0; i<3; i++)
    {
      if(  ( (current[i+1] > current[i]) &&
	     (cu >= current[i]) &&
	     (cu < current[i+1])
	     ) || 
	   ( (current[i+1] < current[i]) &&
	     (cu >= current[i+1]) &&
	     (cu < current[i])
	     )
	   )
	{
	  double i1=current[i];
	  double ip1=current_p[i];
	  double A = o * (up[i]/r -ip1)/2.;
	  double Q = -ip1/(2.*A) ;
	  double t3= sqrt(fabs(Q*Q- (i1-cu)/A));
	  double t2= Q + t3;
	  double t1= Q - t3;
	  // printf("t1/t2 = %lf %lf\n",t1,t2);
	  if(t1 >= 0. && t1 < dt[i]) t3 = t1;
	  else
	    if(t2 >= 0. && t2 < dt[i]) t3 = t2;
	    else
	      {
		printf(" no solution in cur2dt\n");
		tt -= 100000.;
		break;
	      }
	  for(int iter=0; iter <3; iter++)
	    {
	      double i2, ip2;
	      strom(up[i],i1,ip1,t3,i2,ip2);
	      double diff = cu - i2;
	      double t4 = diff/ip2;
	      // printf("correction cu=%lf, i2=%lf, diff=%lf, t4=%lf, t3=%lf\n", cu,i2,diff,t4,t3);
	      t3 += t4;
	    }
	  
	  tt += t3;
	  break;
	}
      else
	{
	  tt += dt[i];
	}
    }
  return tt;
}


void RampFit::stromt()
{
  dt[3] = 0.;
  for(int j=0; j<3; j++)
    {
      dt[3] += dt[j];
      int j1 = j+1;
      double ur = up[j]/r;
      double uro = ur/o;
      double ex0 = exp(-o*dt[j]);
      double ex1 = 1. - ex0;
      //double ex1 = -expm1(-o*dt[j]);
      current[j1] = current[j] + current_p[j]/o * ex1 + uro *(o*dt[j]-ex1);
      current_p[j1] = current_p[j]*ex0 + ur*ex1;
    }
}


int RampFit::shortRamp( int debug )
{
#if 0
  up[0] = upmax;
  up[1] = 0.;
  up[2] = -upmax;
  dt[0] = 1./720.;
  dt[0] = 0.;
  dt[0] = 1./720.;
  stromt();
  printf("current(0,1,2) = %f %f %f\n",current[0],current[1],current[2]);
#endif
  // --- no change ----
  // make limit magnet dependent
  if( fabs(current[0] - goal)+fabs(current_p[0] - goal_p) < 3.e-3) // = 3 mA
    {
      up[0] = upmax;
      up[1] = 0.;
      up[2] = -upmax;
      dt[0] = 0.;
      for(int i=1; i<4; i++)
	{
	  dt[i] = 0.;
	  current[i] = current[0];
	  current_p[i] = current_p[0];
	}
      return 0;
    }
  double u0 = r * current[0] + l * current_p[0];
  double ug = r * goal + l * goal_p;
  double a = upmax/(r*o);
  double ei = exp( (o*current[0]+current_p[0]) * r/upmax );
  double eg = exp( (o*goal + goal_p) * r/upmax );
  double eu = exp( umax/(a*r) );
  if( debug ) printf("u0 = %e, ug = %e, a = %e, ei = %e, eg = %e, eu = %e\n",
	 u0, ug, a, ei, eg, eu);
  
  // --- uds ---- up-down short
  //                             Eg (ip0 Ei + A w Eu - A Ei w)
  //                       Et := -----------------------------
  //                             Eu (A w Eu + ipg Eu - A Eg w)
  
  double et = eg*(current_p[0]*ei + a*o*eu - a*ei*o)/
    (eu*(a*o*eu + goal_p*eu - a*eg*o));
  if( debug ) printf("uds:   et = %e\n", et);
  if(et >= 1.)   // time has to be positive
    {
      double t2 = log(et)/o;
      dt[0] = (umax - u0)/upmax;
      dt[2] = (umax - ug)/upmax;
      dt[1] = t2;
      dt[3] = dt[0]+dt[1]+dt[2];
      if( debug ) printf("uds:   dt = %e %e %e\n", dt[0], dt[1], dt[2]);
      if(dt[0] >= 0. && dt[1] >= 0. && dt[2] >= 0.)
	{
	  up[0] = upmax;
	  up[1] = 0.;
	  up[2] = -upmax;
	  stromt();
	  if( fabs(current[3] - goal) < 1.e-3)  return 0;
	}
    }
  // --- dus ---- down-up short
  //                                 - ip0 + A w Eu Ei - A w
  //                  Et := - -------------------------------------
  //                          Eu Ei (- A w Eu Eg + ipg Eu Eg + A w)
  
  
  et =  - ( -current_p[0] + a*o*eu*ei - a*o )
    / (ei*eu*( - a*o*eu*eg + goal_p*eu*eg + a*o));
  if( debug ) printf("dus:   et = %e\n", et);
  if(et >= 1.)   // time has to be positive
    {
      double t2 = log(et)/o;
      dt[0] = (u0 + umax)/upmax;
      dt[2] = (ug + umax)/upmax;
      dt[1] = t2;
      dt[3] = dt[0]+dt[1]+dt[2];
      if( debug ) printf("dus:   dt = %e %e %e\n", dt[0], dt[1], dt[2]);
      if(dt[0] >= 0. && dt[1] >= 0. && dt[2] >= 0.)
	{
	  up[0] = -upmax;
	  up[1] = 0.;
	  up[2] = upmax;
	  stromt();
	  if( fabs(current[3] - goal) < 1.e-3)  return 0;
	}
    }
  
  // --- udl ---- up-down long
  // 	 D= A*Eg*w/(A*w + ipg )
  // 	 F= (A*w*Ei*Eg - ip0*Ei*Eg)/(A*w + ipg)
  //                  2     1/2        2     1/2
  //       z := D - (D  - F)   , D + (D  - F)
  
  double d = a*eg*o/(a*o + goal_p);
  double f = eg*ei*(a*o - current_p[0])/(a*o + goal_p);
  double w = d*d - f;
  if( debug ) printf("udl:   d=%e, f=%e, w=%e\n", d,f,w);
  //if(fabs(w) < 1.e-10) w=0.;
  if(w >= 0.)
    {
      w = sqrt(w);
      
      double eu1 = d+w;
      if( debug ) printf("udl:   eu1 = %e\n", eu1);
      if(eu1 > 0.0001)
	{
	  double u1 = log(eu1)*upmax/o;
	  dt[0] = (u1 - u0)/upmax;
	  dt[2] = (u1 - ug)/upmax;
	  dt[1] = 0.;
	  dt[3] = dt[0]+dt[1]+dt[2];
	  if( debug ) printf("udl:   dt = %e %e %e\n", dt[0], dt[1], dt[2]);
	  if(dt[0] >= 0. && dt[1] >= 0. && dt[2] >= 0.)
	    {
	      up[0] = upmax;
	      up[1] = 0.;
	      up[2] = -upmax;
	      stromt();
	      if( fabs(current[3] - goal) < 1.e-3)  return 0;
	    }
	}
      
      eu1 = d-w;
      if( debug ) printf("udl:   eu1 = %e\n", eu1);
      if(eu1 > 0.0001)
	{
	  double u1 = log(eu1)*upmax/o;
	  dt[0] = (u1 - u0)/upmax;
	  dt[2] = (u1 - ug)/upmax;
	  dt[1] = 0.;
	  dt[3] = dt[0]+dt[1]+dt[2];
	  if( debug ) printf("udl:   dt = %e %e %e\n", dt[0], dt[1], dt[2]);
	  if(dt[0] >= 0. && dt[1] >= 0. && dt[2] >= 0.)
	    {
	      up[0] = upmax;
	      up[1] = 0.;
	      up[2] = -upmax;
	      stromt();
	      if( fabs(current[3] - goal) < 1.e-3)  return 0;
	    }
	}
    }
  
  // --- dul ---- down-up long
  // 	 D =  ( A*w*Eg*Ei)/(Eg*(ip0+A*w))
  // 	 F =  -( ipg/w - A)* w*Ei*Eg/(ip0+A*w)
  //                  2     1/2        2     1/2
  //       z := D - (D  - F)   , D + (D  - F)
  
  d = ( a * o * ei )/(current_p[0] + a * o);
  f = -( goal_p - a * o ) * ei * eg / (current_p[0] + a * o);
  w = d*d - f;
  
  if( debug ) printf("dul:   d=%e, f=%e, w=%e\n", d,f,w);
  //if(fabs(w) < 1.e-10) w=0.;
  if(w >= 0.)
    {
      w = sqrt(w);
      
      double eu1 = d+w;
      if( debug ) printf("dul:   eu1 = %e\n", eu1);
      if(eu1 > 0.0001)
	{
	  double u1 = log(eu1)*upmax/o;
	  dt[0] = (u0 - u1)/upmax;
	  dt[2] = (ug - u1)/upmax;
	  dt[1] = 0.;
	  dt[3] = dt[0]+dt[1]+dt[2];
	  if( debug ) printf("dul:   dt = %e %e %e\n", dt[0], dt[1], dt[2]);
	  if(dt[0] >= 0. && dt[1] >= 0. && dt[2] >= 0.)
	    {
	      up[0] = -upmax;
	      up[1] = 0.;
	      up[2] = upmax;
	      stromt();
	      if( fabs(current[3] - goal) < 1.e-3)  return 0;
	    }
	}
      
      eu1 = d-w;
      if( debug ) printf("dul:   eu1 = %e\n", eu1);
      if(eu1 > 0.0001)
	{
	  double u1 = log(eu1)*upmax/o;
	  dt[0] = (u0 - u1)/upmax;
	  dt[2] = (ug - u1)/upmax;
	  dt[1] = 0.;
	  dt[3] = dt[0]+dt[1]+dt[2];
	  if( debug ) printf("dul:   dt = %e %e %e\n", dt[0], dt[1], dt[2]);
	  if(dt[0] >= 0. && dt[1] >= 0. && dt[2] >= 0.)
	    {
	      up[0] = -upmax;
	      up[1] = 0.;
	      up[2] = upmax;
	      stromt();
	      if( fabs(current[3] - goal) < 1.e-3)  return 0;
	    }
	}
      if( debug ) printf("w=%f,eu1=%f,c3=%f\n",w,eu1,current[3]);
    }
  
  printf("short:i0=%f,ip0=%f,ig=%f,ipg=%f,r=%f,l=%f,umax=%f,imax=%f,up=%f\n",
	 current[0], current_p[0], goal, goal_p, r, l, umax, imax, upmax);
  
  return 1;
}



int RampFit::minimize(double t)
{
  double u0 = r * current[0] + l * current_p[0];
  double ug = r * goal + l * goal_p;
  double a  = upmax/(r*o);
  double e0 = exp(o*t);
  double ee = e0*e0;
  double ei = exp( o*t +(o*current[0]+current_p[0]) * r/upmax);
  double eg = exp( (o*goal +goal_p) *r/upmax);

  // --- dd ----
  //                                                   2
  //                                    A w (Ei - Eg E0 )
  //                     Eu := ----------------------------------
  //                           E0 (- A w E0 - ipg E0 + ip0 + A w)

  double eu1 = a*o*(ei-eg*ee)/
    (e0*( - a*o*e0 - goal_p*e0 + current_p[0] + a*o));
  // printf("dd = %f\n", eu1);
  if(eu1 > 0.0001)
    {
      double u1 = log(eu1)*upmax/o;
      // printf("u1 = %f\n", u1);
      dt[0] = (u0 - u1)/upmax;
      dt[2] = (u1 - ug)/upmax;
      dt[1] = t - dt[0] - dt[2];
      if(dt[0] >= 0. && dt[1] >= 0. && dt[2] >= 0.)
	{
	  up[0] = -upmax;
	  up[1] = 0.;
	  up[2] = -upmax;
	  stromt();
	  if( fabs(current[3] - goal) < 1.e-3)  return 0;
	}
    }
  // --- uu ----
  //                          (A E0 w - ipg E0 + ip0 - A w) Ei Eg
  //                     z := -----------------------------------
  //                                   E0 w A (- Eg + Ei)
  
  eu1 = ei*eg*(a*e0*o - goal_p*e0 + current_p[0] -a*o)
    / (o*e0*a*(ei-eg));
  // printf("uu = %f\n", eu1);
  if(eu1 > 0.0001)
    {
      double u1 = log(eu1)*upmax/o;
      // printf("u1 = %f\n", u1);
      dt[0] = (u1 - u0)/upmax;
      dt[2] = (ug - u1)/upmax;
      dt[1] = t - dt[0] - dt[2];
      if(dt[0] >= 0. && dt[1] >= 0. && dt[2] >= 0.)
	{
	  up[0] = upmax;
	  up[1] = 0.;
	  up[2] = upmax;
	  stromt();
	  if( fabs(current[3] - goal) < 1.e-3)  return 0;
	}
    }
  
  // --- ud ----
  //	D = - (- ipg*E0*Ei - A*w*E0*Ei + ip0*Ei - A*w*Ei )/(A*w*E0*2)
  //	F = Eg*Ei;
  //                     2         1/2        2         1/2
  //          z := D - (D  - Eg Ei)   , D + (D  - Eg Ei)
  
  
  double d = ei*( goal_p*e0 +a*o*e0 -current_p[0] + a*o)
    / (a*o*e0*2);
  double f = eg*ei;
  double w = d*d - f;
  // printf("w = %f\n", w);
  if(fabs(w) < 1.e-10) w=0.;
  if(w >= 0.)
    {
      w = sqrt(w);
      
      double eu1 = d+w;
      // printf("ud1 = %f\n", eu1);
      if(eu1 > 0.0001)
	{
	  double u1 = log(eu1)*upmax/o;
	  // printf("u1 = %f\n", u1);
	  dt[0] = (u1 - u0)/upmax;
	  dt[2] = (u1 - ug)/upmax;
	  dt[1] = t - dt[0] - dt[2];
	  if(dt[0] >= 0. && dt[1] >= 0. && dt[2] >= 0.)
	    {
	      up[0] = upmax;
	      up[1] = 0.;
	      up[2] = -upmax;
	      stromt();
	      if( fabs(current[3] - goal) < 1.e-3)  return 0;
	    }
	}
      
      eu1 = d-w;
      // printf("ud2 = %f\n", eu1);
      if(eu1 > 0.0001)
	{
	  double u1 = log(eu1)*upmax/o;
	  // printf("u1 = %f\n", u1);
	  dt[0] = (u1 - u0)/upmax;
	  dt[2] = (u1 - ug)/upmax;
	  dt[1] = t - dt[0] - dt[2];
	  if(dt[0] >= 0. && dt[1] >= 0. && dt[2] >= 0.)
	    {
	      up[0] = upmax;
	      up[1] = 0.;
	      up[2] = -upmax;
	      stromt();
	      if( fabs(current[3] - goal) < 1.e-3)  return 0;
	    }
	}
    }
  
  // --- du ----
  //	D=( A*w*E0*E0*Eg - ipg*E0*E0*Eg+ip0*E0*Eg+A*w*E0*Eg)/ (2*A*w*E0*E0)
  //	F=Ei*Eg/E0/E0
  //                                 2     1/2        2     1/2
  //                      z := D - (D  - F)   , D + (D  - F)
  //  
  
  d = eg*(a*o*e0 - goal_p*e0 +current_p[0] + a*o)/
    (2.*a*o*e0);;
  f = ei*eg/ee;
  w = d*d - f;
  // printf("w = %f\n", w);
  if(fabs(w) < 1.e-10) w=0.;
  if(w >= 0.)
    {
      w = sqrt(w);
      
      double du1 = d+w;
      // printf("du1 = %f\n", du1);
      if(du1 > 0.0001)
	{
	  double u1 = log(du1)*upmax/o;
	  // printf("u1 = %f\n", u1);
	  dt[0] = (u0 - u1)/upmax;
	  dt[2] = (ug - u1)/upmax;
	  dt[1] = t - dt[0] - dt[2];
	  if(dt[0] >= 0. && dt[1] >= 0. && dt[2] >= 0.)
	    {
	      up[0] = -upmax;
	      up[1] = 0.;
	      up[2] = upmax;
	      stromt();
	      if( fabs(current[3] - goal) < 1.e-3)  return 0;
	    }
	}
      
      double du2 = d-w;
      // printf("du2 = %f\n", du2);
      if(du2 > 0.0001)
	{
	  double u1 = log(du2)*upmax/o;
	  // printf("u1 = %f\n", u1);
	  dt[0] = (u0 - u1)/upmax;
	  dt[2] = (ug - u1)/upmax;
	  dt[1] = t - dt[0] - dt[2];
	  if(dt[0] >= 0. && dt[1] >= 0. && dt[2] >= 0.)
	    {
	      up[0] = -upmax;
	      up[1] = 0.;
	      up[2] = upmax;
	      stromt();
	      if( fabs(current[3] - goal) < 1.e-3)  return 0;
	    }
	}
    }
  return 1;
}


// #   TTT                                
// #                                
// u0 := i0 *r + ip0 * l;
// ug := ig *r + ipg * l;
// l  := r/w;
// 
// u1 := u0 + up1 * t1;
// up3 := expand((ug - u1)/ t3);
// 
// x1 := expand(-w * t1);
// x2 := expand(-w * t2);
// x3 := expand(-w * t3);
// 
// # e1 := exp(x1);
// # e2 := exp(x2);
// # e3 := exp(x3);
// #   up1 := - r (w ig - w ig e3 + ipg - ipg e3 - w i0 + w i0 e3 - ip0 + ip0 e3
// # 
// #        + w t3 e2 ip0 e1 e3 - w t3 ipg)/
// # 
// #       (w (- t1 + t1 e3 + e2 e3 t3 - e2 e1 e3 t3))
// #
// #
// #
// #                         ig r   ipg r   i0 r   ip0 r   up1 t1
// #                  up3 := ---- + ----- - ---- - ----- - ------
// #                          t3     t3 w    t3     t3 w     t3
// #
// #
// #
// 




int RampFit::minimize3(double t1, double t2, double t3)
{
  double u0 = r * current[0] + l * current_p[0];
  
  double e1 = exp(-o*t1);
  double e2 = exp(-o*t2);
  double e3 = exp(-o*t3);
  
  up[0] = - r * ( (o*goal + goal_p - o*current[0] - current_p[0]) * (1-e3)
		  + o*t3*(current_p[0]*e1*e2*e3 -goal_p) ) /
    (o*( -t1*(1-e3) + t3*e2*e3*(1-e1)));
  
  up[1] = 0.;
  
  up[2] = ( r*(goal*o + goal_p - current[0]*o -current_p[0]) - up[0]*o*t1) / (t3*o); 
  
  double u1 = u0 + up[0]*t1;
  dt[0] = t1; dt[1] = t2; dt[2] = t3; dt[3] = t1+t2+t3;
  
  printf("up  = %f %f %f upmax = %f u1 = %f umax = %f\n",
	 up[0], up[1], up[2], upmax, u1, umax);
  
  stromt();
  // print_result("minimize3");
  return 0;
}

int RampFit::track(double tolerance, int maxRows, float * x, double * y, int nrows, int first)
{
  double i0 = current[0];
  double ip0 =current_p[0];
  double u0 = r * current[0] + l * current_p[0];
  if(first && nrows < maxRows)
    {
      x[nrows] = 1./720.;
      y[nrows] = i0;
      nrows++;
    }
  double ipp0 = 0.;
  double ttt = 0.;
  
  int cnt[3] = {0,0,0};
  for(int i=0; i<3; i++)
    {
      double tt =0;
      double ur = up[i]/r;
      double uro = ur/o;
      while(tt < dt[i])
	{
	  double t1 = sqrt(fabs(4.*tolerance*imax/(o*(ip0-up[i]/r))));
	  if(t1 < 1./720.) t1 = 1./720.;
	  if(tt+t1 > dt[i]) t1 = dt[i]-tt;
	  if(t1 < 0.)
	    {
	      printf("negative");
	    }
	  
	  tt +=t1;
	  ttt += t1;
	  
	  double ex1 = exp(-o*tt);
	  i0 = current[i] + current_p[i]/o * (1.-ex1) + uro *(o*tt-1.+ex1);
	  ip0 = current_p[i]*ex1+ur*(1.-ex1);
	  ipp0 = o*(ur -  current_p[i]) * ex1;
	  u0 = r * i0 + l * ip0;
	  
	  if(nrows < maxRows)
	    {
	      x[nrows] = t1;
	      y[nrows] = i0;
	    }
	  nrows++;
	  cnt[i]++;
	}
    }
  
  double ii =0;
  for(int i=0; i<3; i++)
    {
      double ur = up[i]/r;
      double iii = sqrt( fabs((ur-current_p[i]))) * (exp(-o*dt[i]/2.)-1.);
      ii += iii;
    }
  
  return nrows;
}

// used like:
// wfgRowsN = f.trackp(tol, 1000, wfgRowsX, wfgRowsY, wfgRowsYP, 1);
//#include <stream.h>
int RampFit::trackp(double tolerance, int maxRows, float * x, double * y, double * yp, int nrows, int first)
{
  double i0  = current[0];
  double ip0 = current_p[0];
  double u0 = r * current[0] + l * current_p[0];
  if(first && nrows < maxRows)
    {
      x[nrows] = 1./720.;
      y[nrows] = i0;
      yp[nrows] = ip0;
      nrows++;
    }
  double ipp0 = 0.;
  double ttt = 0.;
  
  int cnt[3] = {0,0,0};
  for(int i=0; i<3; i++)
    {
      double tt = 0;
      double ur = up[i]/r;
      double uro = ur/o;
      while(tt < dt[i])
	{
	  double t1 = sqrt(fabs(4.*tolerance*imax/(o*(ip0-up[i]/r))));
	  if(t1 < 1./720.) t1 = 1./720.;
	  if(tt+t1 > dt[i]) t1 = dt[i]-tt;
	  if(t1 < 0.)
	    {
	      printf("negative");
	    }
	  
	  tt +=t1;
	  ttt += t1;
	  
	  double ex1 = exp(-o*tt);
	  i0 = current[i] + current_p[i]/o * (1.-ex1) + uro *(o*tt-1.+ex1);
	  ip0 = current_p[i]*ex1+ur*(1.-ex1);
	  ipp0 = o*(ur -  current_p[i]) * ex1;
	  u0 = r * i0 + l * ip0;
	  
	  if(nrows < maxRows)
	    {
	      x[nrows] = t1;
	      y[nrows] = i0;
	      yp[nrows] = ip0;
	    }
	  // cout << "nrows = " << nrows << "   t1 = " << t1 << endl;
	  nrows++;
	  cnt[i]++;
	}
    }
  
  double ii =0;
  for(int i=0; i<3; i++)
    {
      double ur = up[i]/r;
      double iii = sqrt( fabs((ur-current_p[i]))) * (exp(-o*dt[i]/2.)-1.);
      ii += iii;
    }
  
  return nrows;
}

// used like:
// wfgRowsN = f.trackp(tol, 1000, wfgRowsX, wfgRowsY, wfgRowsYP, 1);
//#include <stream.h>
int RampFit::trackpp(double tolerance, int maxRows, float * x, double * y, double * yp, double *ypp, int nrows, int first)
{
  double i0  = current[0];
  double ip0 = current_p[0];
  double ipp0 = 0.;
  double u0 = r * current[0] + l * current_p[0];
  if(first && nrows < maxRows)
    {
      x[nrows] = 1./720.;
      y[nrows] = i0;
      yp[nrows] = ip0;
      ypp[nrows] = ipp0;
      nrows++;
    }
  double ttt = 0.;
  
  int cnt[3] = {0,0,0};
  for(int i=0; i<3; i++)
    {
      double tt = 0;
      double ur = up[i]/r;
      double uro = ur/o;
      while(tt < dt[i])
	{
	  double t1 = sqrt(fabs(4.*tolerance*imax/(o*(ip0-up[i]/r))));
	  if(t1 < 1./720.) t1 = 1./720.;
	  if(tt+t1 > dt[i]) t1 = dt[i]-tt;
	  if(t1 < 0.)
	    {
	      printf("negative");
	    }
	  
	  tt +=t1;
	  ttt += t1;
	  
	  double ex1 = exp(-o*tt);
	  i0 = current[i] + current_p[i]/o * (1.-ex1) + uro *(o*tt-1.+ex1);
	  ip0 = current_p[i]*ex1+ur*(1.-ex1);
	  ipp0 = o*(ur -  current_p[i]) * ex1;
	  u0 = r * i0 + l * ip0;
	  
	  if(nrows < maxRows)
	    {
	      x[nrows] = t1;
	      y[nrows] = i0;
	      yp[nrows] = ip0;
	      ypp[nrows] = ipp0;
	    }
	  // cout << "nrows = " << nrows << "   t1 = " << t1 << endl;
	  nrows++;
	  cnt[i]++;
	}
    }
  
  double ii =0;
  for(int i=0; i<3; i++)
    {
      double ur = up[i]/r;
      double iii = sqrt( fabs((ur-current_p[i]))) * (exp(-o*dt[i]/2.)-1.);
      ii += iii;
    }
  
  return nrows;
}
