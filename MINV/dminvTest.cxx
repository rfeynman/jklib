#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

double  dminv(double *, int);

main(int argc, char ** argv)
{
	double m1[5][5];
	double m2[5][5];
	double m3[5][5];
// 	memset(m1, 0, 25*sizeof(double));
// 	for(int i=0; i<5; i++) m1[i][i]=i;


        for(int i = 0; i < 5; i++)
        {
                for(int j = 0; j < 5; j++)
                {
			m1[i][j] = (j+1)/3.-(i+1)/5.-(i+1)*(j+1)/2.;
			if(i == j) m1[i][j] +=1.;
		}
	}



	printf("M1 original matrix\n");
        for(int i=0; i<5; i++)
        {
                for(int j=0; j<5; j++) printf(" %16.8f ", m1[i][j]);
                printf("\n");
        }


	printf("M2 copy of M1\n");
	memcpy(m2, m1, 25*sizeof(double));

// 	for(int i=0; i<5; i++)
// 		for(int j=i+1; j<5; j++)
// 		{
// 			double tenp = m1[i][j];
// 			m1[i][j] = m1[j][i];
// 			m1[j][i] = tenp;
// 		}
	printf("M2 inverse of M1\n");
	double det = dminv( (double *) m2, 5);
	printf("det = %e\n", det);
// 	for(int i=0; i<5; i++)
// 		for(int j=i+1; j<5; j++)
// 		{
// 			double tenp = m1[i][j];
// 			m1[i][j] = m1[j][i];
// 			m1[j][i] = tenp;
// 		}

	printf("M3 =m2* M1\n");
        for(int i = 0; i < 5; i++)
        {
                for(int j = 0; j < 5; j++)
                {
                        m3[i][j]=0.;
                        for(int k = 0; k < 5; k++)
                        {
                                m3[i][j] += m2[i][k] * m1[k][j];
                        }
                }
        }

        for(int i=0; i<5; i++)
        {
                for(int j=0; j<5; j++) printf(" %16.8e ", m3[i][j]);
                printf("\n");
        }

}
