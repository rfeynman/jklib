#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include  "fourier_class.hxx"


void Fourier::addpoint(double f)
{
	if(!data)
	{
		fprintf(stderr," fourier data array not allocated\n");
		exit(-1);
	}
	if(sdata == ndata)
	{
		fprintf(stderr," fourier data array overflow\n");
		exit(-1);
	}
	data[ndata++]=f;
}

void Fourier::hanningFilter()
{
	for(int i=0; i< ndata; i++)
	{
		double s= sin(M_PI*i/ndata);
		data[i] *= 2.*s*s;
	}
}
// this function does the reverse fft the hard way, just for testing, and prints the points to the file in xmgrace format
void  Fourier::makeFunction(FILE * fo, double dt, int maxCoeff)
{

	if(maxCoeff < 0) maxCoeff = ndata/2;
        for(int i=0; i<ndata; i++)
        {   
		double fu=data[0]/2.;
                for(int k=1; k<ndata/2; k++)
                {   
			if(k > maxCoeff) break;
                        double omega = 2.*M_PI*k;
                        fu += data[2*k]*cos(omega*i/ndata) + data[2*k+1]*sin(omega*i/ndata);
                }   
		if( maxCoeff == ndata/2) fu += data[1]*cos(2.*M_PI*i*ndata/2);
                fprintf(fo, "%e %e\n", i*dt, fu*2/ndata);
        }   
        fprintf(fo, "&\n");

}

double  Fourier::calcAmplitudes(double dt)
{
	// dt i the time step between data points.
         amplitudes[0]=data[0];
         amplitudes[ndata/2]=data[1];
	 phases[0]=0.;
	 phases[ndata/2]=0.;
         ipeak=-1;
         peak=0.;
	 int nmax=ndata/2;
         for(int i=1; i<nmax; i++)
         {
              int j=2*i; int k= j+1; 
                amplitudes[i]=sqrt(data[j]*data[j]+data[k]*data[k]);
		phases[i]=atan2(data[k], data[j]);
                if(amplitudes[i] > peak)
                {
                        ipeak=i;
                        peak=amplitudes[i];
			cosPeak=data[j];
			sinPeak=data[k];
                }
         }
         period = (ndata-1)*dt;		// period is the duration of the lowest frequency
	 // f[0] = zero (constant offset)
	 // f[1] =  1/( (ndata-1)*dt )
	 // f[i] =  i/( (ndata-1)*dt )
	 // Q[i] =  f[i] / (1/dt) = f[i]*dt = i/(ndata-1)   /// tune if data are turn-by-turn data




	 // find frequency assuming that the function at the peak is a parabula
	 double jpeak=ipeak;
	 if(ipeak > 0 && ipeak < nmax-1)
	 {
		 jpeak -= ( amplitudes[ipeak+1] - amplitudes[ipeak-1] )/
			  ( 2. *( amplitudes[ipeak+1] + amplitudes[ipeak-1] ) );
	 }


         freq= jpeak/period;
         omega=2.*M_PI*freq;
	 tune = jpeak/(ndata-1);
	 dtune = 1.0/(ndata-1);
         dfreq= 1./period;
         domega=2.*M_PI*dfreq;
	 return freq;
}






// four1 replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1; or replaces
// data[1..2*nn] by nn times its inverse discrete Fourier transform, if isign is input as −1.
// data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn MUST
// be an integer power of 2 (this is not checked for!).

void Fourier::four1(int nn, int isign)
{
	int n=nn << 1;
	int j=1;
	for (int i=1;i<n;i+=2) {
		if (j > i) {
			double temp=data[j];data[j]=data[i]; data[i]=temp;
			temp=data[j-1];data[j-1]=data[i-1]; data[i-1]=temp;
		}
		int m=n/2;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}



	// Here begins the Danielson-Lanczos section of the routine.
	int mmax=2;
	while (n > mmax) {
		int istep=2*mmax;
		double theta=2.*M_PI/(isign*mmax);
		double wtemp=sin(0.5*theta);
		double wpr = -2.0*wtemp*wtemp;
		double wpi=sin(theta);
		double wr=1.0;
		double wi=0.0;
		for (int m=1;m<mmax;m+=2) {
			for (int i=m;i<=n;i+=istep) {
// 				printf("a %5d %5d %5d %e %e %e %e\n", m, i, j, data[i-1], data[i], data[j-1], data[j]);
				j=i+mmax;
				double tempr=wr*data[j-1]-wi*data[j];
				double tempi=wr*data[j]+wi*data[j-1];
				data[j-1]=data[i-1]-tempr;
				data[j]=data[i]-tempi;
				data[i-1] += tempr;
				data[i] += tempi;
// 				printf("b %5d %5d %5d %e %e %e %e\n", m, i, j, data[i-1], data[i], data[j-1], data[j]);
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
// 	for(i=0; i<n; i++) printf("%e\n", data[i]); 
// 	printf("\n");
// 	exit(0);
}


// Calculates the Fourier transform of a set of n real-valued data
// points. Replaces this data (which is stored in array data[1..n])
// by the positive frequency half of its complex Fourier transform.
// The real-valued first and last components of the complex transform are
// returned as elements data[1] and data[2], respectively. n must be a power
// of 2. This routine also calculates the inverse transform of a complex
// data array if it is the transform of real data. (Result in this case
// must be multiplied by 2/n.)

void Fourier::realft(int n, int isign)
{
	double theta=M_PI/double(n);
	double c1=0.5;
	double c2= - 0.5*isign;
	if (isign == 1) {
		four1(n,1);
	} else {
		theta = -theta;
	}
	double wtemp=sin(0.5*theta);
	double wpr = -2.0*wtemp*wtemp;
	double wpi=sin(theta);
	double wr=1.0+wpr;
	double wi=wpi;
	int n2p3=2*n+1;
	for (int i1=2;i1<n;i1+=2) {
		int i2=i1+1;
		int i3=n2p3-i2;
		int i4=i3+1;
		double h1r=   c1*(data[i1]+data[i3]);
		double h1i=   c1*(data[i2]-data[i4]);
		double h2r = -c2*(data[i2]+data[i4]);
		double h2i=   c2*(data[i1]-data[i3]);
		data[i1]=   h1r+wr*h2r-wi*h2i;
		data[i2]=   h1i+wr*h2i+wi*h2r;
		data[i3]=   h1r-wr*h2r+wi*h2i;
		data[i4] = -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		double h1r=data[0];
		data[0] += data[1];
		data[1] = h1r-data[1];
	} else {
		double h1r=data[0];
		data[0]=c1*(data[0]+data[1]);
		data[1]=c1*(h1r-data[1]);
		four1(n,-1);
	}
}
