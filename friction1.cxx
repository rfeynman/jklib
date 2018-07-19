#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "rungeKutta_class.hxx"


//constants
const double elmass =9.10938188e-31;	// kg
const double prmass =1.6726e-27;	// kg
const double elcharge=1.6022e-19;	// C
const double clight=2.997924562e8;	//  m/s
const double epsilon_0=8.8542e-12;	//  C^2/N*m^2
const double elradius=elcharge*elcharge/(4.*M_PI*epsilon_0*elmass*clight*clight);
const double prradius=elcharge*elcharge/(4.*M_PI*epsilon_0*prmass*clight*clight);
const double sq2=sqrt(2.);
double beta(double gamma){ return sqrt(1.-1./(gamma*gamma));}


// Ion paramters
const double nions=1e9;
const double Z=79;
const double A=197;
const double gamma_i=107.2;
const double frev=78000;
const double emit95=15e-6;   // m
const double ionSigmaZ=0.3;
const double ionEnergySpread=0.001;
const double numIons=1.e9;
const double ionBetaFunc=60.; //m
const double B=5;			// Tesla
const double Lcool=13.;		// length of the cooler [m]

// derived ion parameters
double emitn=emit95/6.;
double emit=emitn/gamma_i;
double ionSigmaX=sqrt(emit*ionBetaFunc);
double ionSigmaY=ionSigmaX;
double Vit=beta(gamma_i)*clight*gamma_i*sqrt(emit/(2.*ionBetaFunc));
double Vit22=2.*Vit*Vit;
double Vil=ionEnergySpread*beta(gamma_i)*clight;
double Vil22=2.*Vil*Vil;
double ionDensityFactor=8./(pow(2*M_PI,3./2.)*Vit*Vit*Vil);
double flightTime=Lcool/(beta(gamma_i)*gamma_i*clight);



// Wiggler and amplifier parameter
const double B0=0.125
const double lambda_u=0.0175
const double M=3.;   
const double K=elcharge*B0*lambda_u*clight/(2.*M_PI*elmass*clight*clight);
const double xsi=K*K/(2.+K*K);
const double lambda=lambda_u*(1.+K*K/2.)/(2.*gamma_i*gamma_i);
const double k_lambda=2.*M_PI/lambda;
const double B0=0.125





// functions

class MyRuInner : public RungeKutta
{
	double gg(double vx, double vy, double vz, double veff);
public:
	double veff, RhoLar, elAx,omegaPlasma; 
	double z;
	MyRuInner() : RungeKutta(1) {}

               virtual void derivs(double x,           double * y,         double * dydx       )
	       {
		       dydx[0]=M_PI/2.*gg(x, 0., z, veff)*x;
	       }
} myRuInner;

class MyRuOuter : public RungeKutta
{
public:
	MyRuOuter() : RungeKutta(1) {}

               virtual void derivs(double x,           double * y,         double * dydx       )
	       {
			int nok, nbad;
			double ystart[1] ={0.};
			myRuInner.z=x;
			myRuInner.odeint( ystart,  0.,   100.*Vit,   1.e-5  , 10,   1.e-10,  nok,   nbad);
		       	dydx[0]=ystart[0];
	       }
} myRuOuter;

double MyRuInner::gg(double vx, double vy, double vz, double veff)
{
	double v2 = vx*vx+vy*vy+vz*vz;
	double v = sqrt(v2);
	double rhi = ionDensityFactor*exp(-vx*vx/Vit22 - vy*vy/Vit22 - vz*vz/Vil22);
	double rmin = elradius*clight*clight + RhoLar*v2;
	double r1 = v/omegaPlasma;
	double r2 = v*flightTime;
	double r3 = 2.* elAx;
	double r=r1;
	if(r2 < r) r=r2;
	if(r3 < r) r=r3;
	double a= (r*v2 + rmin) /  rmin;
	double lcx = log(a);
	double vv = sqrt(v2+veff*veff);
	return  v2*rhi*lcx/(vv*vv*vv);
}



main(int argc, char ** argv)
{

// electron parameters
double numElectrons=1.2e11;
double overlap=1.;
double elEpsXn= 40.e-6;	// rad*m
double elEpsYn= 40.e-6;
double elSigmaZ=0.075;//m
double elEnergySpread=1e-4;

// derived electron parameters
// double elAx=ionSigmaX*overlap;
// double elAy=ionSigmaY*overlap;
double elAx=0.001;
double elAy=0.001;
double elAz=2*elSigmaZ;
// double elSigmaX=ionSigmaX*overlap;
// double elSigmaY=ionSigmaY*overlap;
double elSigmaX=elAx/2.;
double elSigmaY=elAx/2.;
double elDensity=numElectrons/(M_PI*elAx*elAy*elAz*gamma_i);
double Vet=clight/sqrt(2.)*elEpsXn/elSigmaX;
double Vel=clight*elEnergySpread;
double vd=2.*M_PI*elcharge*elDensity*elAx/(B*epsilon_0);
double veff=sqrt(vd*vd+Vel*Vel);
double RhoLar=Vet*elmass/(elcharge*B);
double omegaPlasma=clight*sqrt(4.*M_PI*elDensity*elradius);
double vel=elEnergySpread*clight;



double C1 = 4*Z*Z*M_PI*elradius*prradius*prmass*elDensity*clight*clight*clight*clight*flightTime;
double C1x = erf(elAx/(sq2*ionSigmaX));
double C1y = erf(elAy/(sq2*ionSigmaY));
double C1z = erf(elAz/(sq2*ionSigmaZ));
printf("C1 = %e\n", C1);
printf("C1x = %e\n", C1x);
printf("C1y = %e\n", C1y);
printf("C1z = %e\n", C1z);
printf("emitn = %e\n", emitn);
printf("emit = %e\n", emit);
printf("ionSigmaX = %e\n", ionSigmaX);
printf("ionSigmaY = %e\n", ionSigmaY);
printf("Vit = %e\n", Vit);
printf("Vit22 = %e\n", Vit22);
printf("Vil = %e\n", Vil);
printf("Vil22 = %e\n", Vil22);
printf("ionDensityFactor = %e\n", ionDensityFactor);

printf("elAy = %e\n", elAy);
printf("elSigmaX = %e\n", elSigmaX);
printf("elSigmaY = %e\n", elSigmaY);
printf("Vet = %e\n", Vet);
printf("RhoLar = %e\n", RhoLar);
printf("veff = %e\n", veff);
printf("elDensity = %e\n", elDensity);
printf("omegaPlasma = %e\n", omegaPlasma);
printf("vel = %e\n", vel);
printf("flightTime = %e\n", flightTime);


double D0=1e-6;
double nop=0.6;
double beta_e=200.;

double Id=2.*D0*nop+sqrt(1.-4.*k_lambda*k_lambda*D0*D0*elepsx/beta_e)/(elEnergySpread*k_lambda);


	myRuInner.veff=veff;
	myRuInner.RhoLar=RhoLar;
	myRuInner.elAx=elAx;
	myRuInner.omegaPlasma=omegaPlasma;
	int nok, nbad;
	double ystart[1] ={0.};
// 	myRuOuter.odeint( ystart,double x1,double x2,double eps,double h1,double hmin,int & nok,int & nbad);
	myRuOuter.odeint( ystart,       0.,   100.*Vil,   1.e-5  , 1    ,   1.e-10   ,      nok,      nbad);
	printf("II  %d %d %e\n", nok, nbad, ystart[0]);
	double deltaE =  C1*C1x*C1y*C1z*ystart[0];
	printf("deltaE  %e\n",  deltaE);
	double tauDamping=A*prmass*(2*Vit*Vit+Vil*Vil)/(deltaE*frev);
	printf("tauDamping = %e\n",  tauDamping);


}
