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
double beta(double gamma){ return sqrt(1.-1./(gamma*gamma));}

// Ion paramters
const double Z=79;
const double A=197;
const double gamma_i=108;
const double emit95=15e-6;   // m
const double ionBunchLength=0.3;
const double ionEnergySpread=0.001;
const double numIons=1.e9;
const double ionBetaFunc=60.; //m
const double B=5;			// Tesla
const double Lcool=13.;		// length of the cooler [m]

// derived ion parameters
const double emitn=emit95/6.;
const double emit=emitn/gamma_i;
const double ionSigmaX=sqrt(emit*ionBetaFunc);
const double ionSigmaY=ionSigmaX;
const double Vit=beta(gamma_i)*clight*gamma_i*sqrt(emit/(2.*ionBetaFunc));
const double Vit22=2.*Vit*Vit;
const double Vil=ionEnergySpread*beta(gamma_i)*clight;
const double Vil22=2.*Vil*Vil;
const double ionDensityFactor=8./(pow(2*M_PI,3./2.)*Vit*Vit*Vil);
const double flightTime=Lcool/(beta(gamma_i)*gamma_i*clight);

// electron parameters
double numElectrons=1.2e11;
double overlap=1.;
double elEpsXn= 40.e-6;	// rad*m
double elEpsYn= 40.e-6;
double elBunchLength=0.15;	//m
double elEnergySpread=1e-4;

double elSigmaX=ionSigmaX*overlap;
double elSigmaY=ionSigmaY*overlap;
double Vet=clight/sqrt(2.)*elEpsXn/elSigmaX;
double RhoLar=Vet*elmass/(elcharge*B);
double veff=0.;
double elDensity=numElectrons/(M_PI*elSigmaX*elSigmaY*elBunchLength*2.5*gamma_i);
double omegaPlasma=clight*sqrt(4.*M_PI*elDensity*elradius);
double vel=elEnergySpread*clight;

double Rhomax(double vx, double vy, double vz);
double ionDensity(double vx,double vy,double vz);
double Lcx(double vx, double vy, double vz);
double V2(double vx, double vy, double vz);
double gg(double vx, double vy, double vz, double veff);

class MyRuInner : public RungeKutta
{
	double z;
public:
	MyRuInner(double z) : RungeKutta(1) {this->z=z;}

               virtual void derivs(double x,           double * y,         double * dydx       )
	       {
		       dydx[0]=gg(x, 0., z, veff)*x;
	       }
} myRuInner(1);

class MyRuOuter : public RungeKutta
{
public:
	MyRuOuter() : RungeKutta(1) {}

               virtual void derivs(double x,           double * y,         double * dydx       )
	       {
			int nok, nbad;
			double ystart[1] ={0.};
			myRuInner.odeint( ystart,       0.,   100.*Vit,   1.e-5  , 10    ,   1.e-10   ,      nok,      nbad);
		       	dydx[0]=ystart[0];
	       }
} myRuOuter;


main(int argc, char ** argv)
{
	int nok, nbad;
	double ystart[1] ={0.};
// 	a.odeint( ystart,double x1,double x2,double eps,double h1,double hmin,int & nok,int & nbad);
	myRuOuter.odeint( ystart,       0.,   100.*Vit,   1.e-5  , 10    ,   1.e-10   ,      nok,      nbad);
	printf("%d %d %e\n", nok, nbad, ystart[0]);

}



double gg(double vx, double vy, double vz, double veff)
{
	double v2 = V2(vx,vy,vz);
	double rhi = ionDensity(vx,vy,vz);
	double lcx = Lcx(vx,vy,vz);
	return  v2*rhi*lcx/pow(v2+veff*veff,3./2.);
}

double V2(double vx, double vy, double vz)
{
	return vx*vx+vy*vy+vz*vz;
}
	
double Lcx(double vx, double vy, double vz)
{
	double v2 = V2(vx,vy,vz);
	return log(Rhomax(vx,vy,vz)*v2+elradius*clight*clight+RhoLar*v2/(RhoLar*v2 +elradius*clight*clight));
}

double Rhomax(double vx, double vy, double vz)
{
	double v2 = V2(vx,vy,vz);
	double r1 = v2/omegaPlasma;
	double r2 = v2*flightTime;
	double r3 = 2.* elSigmaX;
	double r=r1;
	if(r2 < r) r=r2;
	if(r3 < r) r=r3;
	return r;
}



double ionDensity(double vx,double vy,double vz)
{
	return ionDensityFactor*exp(-vx*vx/Vit22 - vy*vy/Vit22 - vz*vz/Vil22);
}



