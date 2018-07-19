#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include  "fourier_class.hxx"

main(int argc, char ** argv)
{

        FILE * fo = fopen("_f.agr", "w");
        if(!fo)
        {   
                perror("_f.agr");
                exit(-1);
        }   


        FILE * fp = fopen("_p.agr", "w");
        if(!fp)
        {   
                perror("_p.agr");
                exit(-1);
        }   


	const int n=128;
	double data[n];

        for(int i=0; i<n; i++)
        {
                data[i] = 5.*sin(2.*M_PI*i/n) + 7.*cos(2.*M_PI*9.*i/n);
		fprintf(fo, "%d %e\n", i, data[i]);
        }
	fprintf(fo, "&\n");

	Fourier g(n);
        for(int i=0; i<n; i++) g.addpoint( data[i]);

        g.realft(n/2,1);
	double freq=g.calcAmplitudes(1.);
	fprintf(stderr, "freq=%e\n", freq);
        for(int i=0; i<n/2; i++) fprintf(fp, "%d %e\n", i, g.amplitudes[i]);

        for(int i=0; i<n; i++)
        {
		printf("%d %e\n", i, g.data[i]);
	}

	g.makeFunction( fo);

//          for(int i=1; i<n/2; i++)
//          {   
//               int j=2*i; int k= j+1; 
// 	      double t=g.data[j]; g.data[j]=g.data[k]; g.data[k]=t;
// 	 }
        g.realft(n/2,-1);
	fprintf(fo, "&\n");
	g.makeFunction( fo);

	fprintf(fp, "&\n");
        for(int i=0; i<n/2; i++) fprintf(fp, "%d %e\n", i, g.amplitudes[i]);


}

