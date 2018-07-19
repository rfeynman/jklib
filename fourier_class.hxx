#include <stdio.h>
#include <stdlib.h>
#include <math.h>
class Fourier
{
public:
	int ndata; // number of data points filled
	int sdata; // size of the buffer
	double * data;
	double * amplitudes;
	double * phases;
        double period;		// duration of the lowest freqency period
        double freq;		// frequency of the highest peak in the power specturm
        double dfreq;		// frequency difference of two adjacent points
        double omega;
        double domega;
	double tune;		// tune if data is turn-by-turn data
	double dtune;		// tune if data is turn-by-turn data
	double peak;		// peak value in the power spectrum
	int ipeak;		// peak index in the power spectrum
	double cosPeak, sinPeak;	// 
	Fourier()
	{
		data=0;
		ndata=0;
		sdata=0;
		amplitudes = 0;
	}

	Fourier(int n)
	{
		data= new double[n];
		amplitudes = new double[n];
		phases = new double[n];
		ndata=0;
		sdata=n;
	}
	virtual ~Fourier()
	{
		delete [] data;
		delete [] amplitudes;
		delete [] phases;
	}
// realft() calculates the Fourier transform of a set of n real-valued data
// points. Replaces this data (which is stored in array data[1..n])
// by the positive frequency half of its complex Fourier transform.
// The real-valued first and last components of the complex transform are
// returned as elements data[0] and data[1], respectively. n must be a power
// of 2. This routine also calculates the inverse transform of a complex
// data array if it is the transform of real data. (Result in this case
// must be multiplied by 2/n.)

	void realft(int n, int isign);

// four1() replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1; or replaces
// data[1..2*nn] by nn times its inverse discrete Fourier transform, if isign is input as −1.
// data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn MUST
// be an integer power of 2 (this is not checked for!).

	void four1(int nn, int isign);

// addpoint appends the data point f to the internal array data.
	void addpoint(double f);

// calcAmplitudes() calculates the amplitudes of the data array and store them in the amplitudes array. 
// dt is the time step between data points. calcAmplitudes calculates from that and the peak value the frequency.
	double  calcAmplitudes(double dt);
// makeFunction function does the reverse fft the hard way, just for testing, and prints the points to the file in xmgrace format
	void  makeFunction(FILE * fo, double dt=1., int maxCoeff=-1);
// hanningFilter applies the Hanning filter 2*sin^2(2 pi i/n) to the data
	void hanningFilter();

};


