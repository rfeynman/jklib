#include "amoeba_class.hxx"

Amoeba::Amoeba(int ndim)
{
	Amoeba::ndim = ndim;
	mpts = ndim+1;
	p = new double*[mpts];
	for(int i=0; i<mpts; i++)
		p[i]=new double[ndim];
	y = new double[mpts];
	psum = new double[ndim];
	ptry = new double[ndim];
	NMAX=5000;
	ALPHA=1.0;
	BETA=0.5;
	GAMMA=2.0;
}


Amoeba::~Amoeba()
{
	delete psum;
	delete ptry;
	delete y;
	for(int i=0; i<mpts; i++)
		delete p[i];
	delete p;
}



double Amoeba::amotry(int ihi, double fac)
{
	double fac1=(1.0-fac)/ndim;
	double fac2=fac1-fac;
	for(int j=0;j<ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
	double ytry=funk(ptry);
	nfunk++;
	if(ytry < y[ihi]) {
		y[ihi]=ytry;
		for(int j=0;j<ndim;j++) {
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j]=ptry[j];
		}
	}
	return ytry;
}


int Amoeba::minimize()
{
	int i,ilo,ihi,inhi;

	nfunk=0;

	for(int j=0;j<ndim;j++)
	{
		double sum=0.0;
		for(i=0;i<mpts;i++) sum += p[i][j];
		psum[j]=sum;
	}
	while(1)
	{
		ilo=0;
		ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
		for(i=0;i<mpts;i++)
		{
			if(y[i] < y[ilo]) ilo=i;
			if(y[i] > y[ihi])
			{
				inhi=ihi;
				ihi=i;
			}
			else
			if(y[i] > y[inhi])
				if(i != ihi) inhi=i;
		}
		double rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
		if(rtol < ftol) break;
		if(nfunk >= NMAX) return nfunk;
		double ytry=amotry(ihi,-ALPHA);
		if(ytry <= y[ilo])
			ytry=amotry(ihi,GAMMA);
		else if(ytry >= y[inhi])
		{
			double ysave=y[ihi];
			ytry=amotry(ihi,BETA);
			if(ytry >= ysave)
			{
				for(i=0;i<mpts;i++)
				{
					if(i != ilo)
					{
						for(int j=0;j<ndim;j++)
						{
							psum[j]=0.5*(p[i][j]+p[ilo][j]);
							p[i][j]=psum[j];
						}
						y[i]=funk(psum);
					}
				}
				nfunk += ndim;
				for(int j=0;j<ndim;j++)
				{
					double sum=0.0;
					for(i=0;i<mpts;i++) sum += p[i][j];
					psum[j]=sum;
				}
			}
		}
	}
	return 0;
}



void Amoeba::init_p(double *base, double *delta, double tolerance)
{
	ftol = tolerance;
	for(int i=0; i<=ndim; i++)
		for(int j=0; j<ndim; j++) p[i][j]=base[j];
	for(int j=0; j<ndim; j++)
		p[j+1][j] += delta[j];
	for(int i=0; i<=ndim; i++)
		y[i] = funk(p[i]);
}

