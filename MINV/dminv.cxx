//   14/11/86 611282107  member name  dminv    (petlib)      fortran77
//
//     .................................................................
//
//        subroutine dminv
//
//        purpose
//           invert a matrix
//
//        usage
//           double d = dminv(a,n)
//
//        description of parameters
//           a - input matrix, destroyed in computation and replaced by
//               resultant inverse.
//           n - order of matrix a
//           d - resultant determinant
//           l - work vector of length n
//           m - work vector of length n
//
//        remarks
//           matrix a must be a general matrix
//
//        subroutines and function subprograms required
//           none
//
//        method
//           the standard gauss-jordan method is used. the determinant
//           is also calculated. a determinant of zero indicates that
//           the matrix is singular.
//
//        translated into fortran77 11.86 jkl
//        translated into c++       05.98 jkl
//     .................................................................
//

#include <math.h>
#include <stdio.h>


double dminv(double *a, int n)
{
      int * l = new int[n];
      int * m = new int[n];
      double d=1.0;
	printf("matrix of %d\n",n);
	for(int i=0; i<n*n; i++) printf("%3.1lf  ", a[i]);
	printf("\n");

      int nk=-n;
      for(int k=0; k<n; k++)
      {
	 printf("k  %8d\n", k +1);
         nk=nk+n;
         l[k]=k;
         m[k]=k;
         int kk=nk+k;
	 printf("kk %8d\n", kk+1);
         double biga=a[kk];
// 	 printf("%15.8e %d\n", biga, kk);
	 for(int j=k; j<n; j++)
	 {
	 printf("j  %8d\n", j +1);
            int iz=n*j;
	 printf("iz %8d\n", iz);
	    for(int i=k; i<n; i++)
	    {
               int ij=iz+i;
	 printf("ij %8d\n", ij+1);
               if(fabs(biga) < fabs(a[ij]))
	       {
                  biga=a[ij];
                  l[k]=i;
                  m[k]=j;
               }
	    }
	 }
	 for(int s=0; s<5; s++) printf("%8d  ", l[s]+1); printf("\n");
	 for(int s=0; s<5; s++) printf("%8d  ", m[s]+1); printf("\n");
         int j=l[k];
         if(j > k) {
            int ki=k-n;
	    for(int i=0; i<n; i++)
	    {
               ki=ki+n;
	 printf("ki %8d\n", ki+1);
               double hold = -a[ki];
// 		  printf("h1 %15.8e\n", hold);
               int ji=ki-k+j;
	 printf("ji %8d\n", ji+1);
               a[ki]=a[ji];
               a[ji] =hold;
	    }
         }
         int i=m[k];
         if(i > k) {
            int jp=n*i;
	    for(int j=0; j<n; j++)
	    {
               int jk=nk+j;
	 printf("jk %8d\n", jk+1);
               int ji=jp+j;
	 printf("ji %8d\n", ji+1);
               double hold= - a[jk];
               a[jk]=a[ji];
               a[ji] =hold;
	    }
         }
         if(biga == 0.0) {
	    delete [] l;
	    delete [] m;
            return 0.0;
         }
	 for(int i=0; i<n; i++)
	 {
            if(i != k) {
               int ik=nk+i;
	 printf("ik %8d\n", ik+1);
               a[ik]=a[ik]/(-biga);
            }
	 }

	 for(i=0; i<n; i++)
	 {
            int ik=nk+i;
	 printf("ik %8d\n", ik+1);
            double hold=a[ik];
            int ij=i-n;
	 printf("ij %8d\n", ij+1);
	    for(int j=0; j<n; j++)
	    {
               ij=ij+n;
	 printf("ij %8d\n", ij+1);
               if(i != k && j != k) {
                  int kj=ij-i+k;
	 printf("kj %8d\n", kj+1);
                  a[ij]=hold*a[kj]+a[ij];
               }
	    }
	 }
         int kj=k-n;
	 printf("kj %8d\n", kj+1);
	 for(j=0; j<n; j++)
	 {
            kj=kj+n;
	 printf("kj %8d\n", kj+1);
            if(j != k) a[kj]=a[kj]/biga;
	 }
         d=d*biga;
	 printf("kk %8d\n", kk+1);
         a[kk]=1.0/biga;
      }

      for(int k=n-2; k>=0; k--)
      {
         int i=l[k];
         if(i > k) {
            int jq=n*k;
            int jr=n*i;
	    for(int j=0; j<n; j++)
	    {
               int jk=jq+j;
	 printf("jk %8d\n", jk+1);
               double hold=a[jk];
               int ji=jr+j;
	 printf("ji %8d\n", ji+1);
               a[jk]=-a[ji];
               a[ji] =hold;
	    }
         }
         int j=m[k];
         if(j > k) {
            int ki=k-n;
	    for(int i=0; i<n; i++)
	    {
               ki=ki+n;
	 printf("ki %8d\n", ki+1);
               double hold=a[ki];
               int ji=ki-k+j;
	 printf("ji %8d\n", ji+1);
               a[ki]=-a[ji];
               a[ji] =hold;
	    }
         }
      }
      delete [] l;
      delete [] m;
      return d;
}
