
int mxequ(double *a, double *b, int i, int j, int idim)
{

/*
c cern proglib# f109    m ==	      .version kernfor  1.0   680220
c orig. 01/01/64 rkb et al
c   solution of a linear == ation system
c   call form:
c   ierr=mxequ(a,b,i,j,idim)
c   a    idim by idim matrix, the i by i part of that matrix is used.
c   b    i by j matrix
c   a    is transformed, so that the product of the diagonal elements
c	  gives the determinat, the result is stored in b.
c	  in case of error the column no with a(l,l)=0 is stored in i. NOT!
c
*/
	int n,ll,lm,l,m,l1m,l1,l1l,lplus1,mn,ml,ln,nn;
	double diag;

	if(i == 0) return(-1);

/* 	          transformation of a-matrix; */
	ll=0;
	for(l=0; l < i; l++)
	{
		diag=a[ll];
		if(diag == 0.) return(ll+1);
		diag=-1./diag;
		lm=l;

		for(m=0; m < i; m++)
		{
			if(l != m) a[lm] *= diag;
			lm=lm+idim;
		}

		lplus1=l+1;
		l1l=ll+1;

	   	for(l1=lplus1; l1 < i; l1++)
	   	{
			l1m=l1;
			lm=l;
	
	   		for(m=0;m <i; m++)
	   		{
				if(l != m) a[l1m] += a[l1l]*a[lm];
   				l1m=l1m+idim;
   				lm=lm+idim;
			}
   			l1l=l1l+1;
		}
   		ll=ll+idim+1;
	}


/* 	          transformation of b-matrix; */
	if(i*j == 0) return(-1);

	ml=0;
	for(l=0; l < i; l++)
	for(m=0; m < j; m++)
	{
		mn=m;
		ln=l;

	   	for(n=0; n < i; n++)
	   	{
			if(l != n) b[mn] += b[ml]*a[ln];
			mn=mn+j;
			ln=ln+idim;
		}
		ml=ml+1;
	}


/* 	          puts final result in b-matrix; */
	nn=0;
	mn=0;
	for(n=0; n <i; n++)
	{
		diag=1./a[nn];
		for(m=0; m < j; m++)
		{
			b[mn]=b[mn]*diag;
			mn=mn+1;
		}
		nn=nn+idim+1;
	}
	return(0);
}
