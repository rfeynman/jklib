#include "rfn3.hxx"
#include "wireup.hxx"


Mat::Mat()
{
	unitA();
}


void  Mat::unitA()
{
	memset(a,0, 25*sizeof(double));
	a[0][0] = 1.;
	a[1][1] = 1.;
	a[2][2] = 1.;
	a[3][3] = 1.;
	a[4][4] = 1.;
};
void  Mat::inver1(double v[][4])
{

//
//    "inver1" calculates the matrix v=(1-a)**-1 from the 4x4 part of the matrix a.
//
//
//    j. kewisch 29.11.79
//    changes: ft77,real*8    5.85
      int ll[3],kk[3];
      double m[4][4];
      for(int i=0; i<4; i++) for(int j=0; j<4; j++) m[i][j]= -a[i][j];
      for(int i=0; i<4; i++) m[i][i]=m[i][i]+1.0;


      double det=(m[0][0]*m[1][1]-m[0][1]*m[1][0])*(m[2][2]*m[3][3]-m[2][3]*m[3][2])
          	-(m[0][0]*m[1][2]-m[0][2]*m[1][0])*(m[2][1]*m[3][3]-m[2][3]*m[3][1])
          	+(m[0][0]*m[1][3]-m[0][3]*m[1][0])*(m[2][1]*m[3][2]-m[2][2]*m[3][1])
          	+(m[0][1]*m[1][2]-m[0][2]*m[1][1])*(m[2][0]*m[3][3]-m[2][3]*m[3][0])
          	-(m[0][1]*m[1][3]-m[0][3]*m[1][1])*(m[2][0]*m[3][2]-m[2][2]*m[3][0])
          	+(m[0][2]*m[1][3]-m[0][3]*m[1][2])*(m[2][0]*m[3][1]-m[2][1]*m[3][0]);


      for(int l=0; l<4; l++)
      for(int k=0; k<4; k++)
      {
      	int il= 0;
      	int ik= 0;
      	for(int i=0; i<4; i++)
      	{
      		if(i != l) ll[il++]=i;
      		if(i != k) kk[ik++]=i;
      	}
      	double t =(m[kk[0]][ll[0]]*m[kk[1]][ll[1]]*m[kk[2]][ll[2]]
                  -m[kk[0]][ll[0]]*m[kk[2]][ll[1]]*m[kk[1]][ll[2]]
                  -m[kk[0]][ll[1]]*m[kk[1]][ll[0]]*m[kk[2]][ll[2]]
                  +m[kk[0]][ll[1]]*m[kk[2]][ll[0]]*m[kk[1]][ll[2]]
                  +m[kk[0]][ll[2]]*m[kk[1]][ll[0]]*m[kk[2]][ll[1]]
                  -m[kk[0]][ll[2]]*m[kk[1]][ll[1]]*m[kk[2]][ll[0]])/det;
      	v[l][k]=  (l+k) & 1 ? -t : t;
      }


}

int Mat::eigenwerte()
{



//                ------------calculation of eigen values of matrix a
      spur1=a[0][0]+a[1][1];
      spur2=a[2][2]+a[3][3];
      double det=(a[0][2]+a[3][1])*(a[1][3]+a[2][0])
         -(a[3][0]-a[1][2])*(a[2][1]-a[0][3]);
      double f0=spur1-spur2;
      double f1=spur1+spur2;
      double f2=f0*f0+4.0*det;


//                 -----------if f2 negative, the  abs(eigenvalue) is not 1
      if(f2 < 0.0) return 1;


      f2=sqrt(f2);
      if(f0 >= 0.0)
      {
         lamg1=(f1+f2)*0.5;
         lamg2=(f1-f2)*0.5;
      }
      else
      {
         lamg1=(f1-f2)*0.5;
         lamg2=(f1+f2)*0.5;
      }


      f1=lamg1*lamg1-4.0;
      f2=lamg2*lamg2-4.0;
      Complex ca1(f1);
      Complex ca2(f2);
      ca1=sqrt(ca1);
      ca2=sqrt(ca2);
      clam1=(lamg1+ca1)*0.5;
      clam2=(lamg2+ca2)*0.5;


      if(lamg1*lamg1 >= 4.0) return 2;
      if(lamg2*lamg2 >= 4.0) return 2;


//      printf("eigewerte (%f, %f) (%f,%f)\n", clam1.r, clam1.i, clam2.r, clam2.i);
//      double qx= atan2(clam1.i,(1.0+clam1.r))/M_PI;
//     	double qz= atan2(clam2.i,(1.0+clam2.r))/M_PI;
//	printf(" qx = %f , qz = %f\n", qx,qz);



      return 0;
}
int Mat::eigenvector()
{
//                ---------calculation of eigenvectors
      	double detb=a[2][0]*a[3][1]-a[3][0]*a[2][1];
      	double detc=a[0][2]*a[1][3]-a[1][2]*a[0][3];
      	double fak1=spur2-lamg1;
      	double fak2=spur1-lamg2;


	Complex cw1[4],cw2[4];

      	cw1[0]=-(a[2][0]*a[1][2]+a[3][0]*a[1][3])/fak1+a[1][0];
      	cw1[1]=(a[2][0]*a[0][2]+a[3][0]*a[0][3]+detb)/fak1-(a[0][0]-clam1);
      	cw1[2]=-((a[0][2]+a[3][1])*cw1[0]+(a[1][2]-a[3][0])*cw1[1])/fak1;
      	cw1[3]=-((a[0][3]-a[2][1])*cw1[0]+(a[1][3]+a[2][0])*cw1[1])/fak1;

      	cw2[2]=-(a[0][2]*a[3][0]+a[1][2]*a[3][1])/fak2+a[3][2];
      	cw2[3]=(a[0][2]*a[2][0]+a[1][2]*a[2][1]+detc)/fak2-(a[2][2]-clam2);
      	cw2[0]=-((a[2][0]+a[3][1])*cw2[2]+(a[0][3]-a[2][1])*cw2[3])/fak2;
      	cw2[1]=-((a[2][1]-a[0][3])*cw2[2]+(a[3][1]+a[0][2])*cw2[3])/fak2;

//   prove that a * cw = lambda * cw
	Complex ce1[4],ce2[4];
	ce1[0] = a[0][0]* cw1[0] + a[1][0]*cw1[1] + a[2][0]*cw1[2] + a[3][0]*cw1[3];
	ce1[1] = a[0][1]* cw1[0] + a[1][1]*cw1[1] + a[2][1]*cw1[2] + a[3][1]*cw1[3];
	ce1[2] = a[0][2]* cw1[0] + a[1][2]*cw1[1] + a[2][2]*cw1[2] + a[3][2]*cw1[3];
	ce1[3] = a[0][3]* cw1[0] + a[1][3]*cw1[1] + a[2][3]*cw1[2] + a[3][3]*cw1[3];

	ce2[0] = a[0][0]* cw2[0] + a[1][0]*cw2[1] + a[2][0]*cw2[2] + a[3][0]*cw2[3];
	ce2[1] = a[0][1]* cw2[0] + a[1][1]*cw2[1] + a[2][1]*cw2[2] + a[3][1]*cw2[3];
	ce2[2] = a[0][2]* cw2[0] + a[1][2]*cw2[1] + a[2][2]*cw2[2] + a[3][2]*cw2[3];
	ce2[3] = a[0][3]* cw2[0] + a[1][3]*cw2[1] + a[2][3]*cw2[2] + a[3][3]*cw2[3];

	Complex cl1[4],cl2[4];
	cl1[0] =  cw1[0] *clam1;
	cl1[1] =  cw1[1] *clam1;
	cl1[2] =  cw1[2] *clam1;
	cl1[3] =  cw1[3] *clam1;

	cl2[0] =  cw2[0] *clam2;
	cl2[1] =  cw2[1] *clam2;
	cl2[2] =  cw2[2] *clam2;
	cl2[3] =  cw2[3] *clam2;

// 	for(int i=0; i<4; i++)
//	printf("%f %f %f %f %f %f \n", cw1[i].r, cw1[i].i, ce1[i].r, ce1[i].i, cl1[i].r, cl1[i].i);

// 	for(int i=0; i<4; i++)
// 	printf("%f %f %f %f %f %f \n", cw2[i].r, cw2[i].i, ce2[i].r, ce2[i].i, cl2[i].r, cl2[i].i);

// 	printf("abs clam = %f %f\n", abs(clam1), abs(clam2));
//                 --------change from complex to real values
      	for(int i=0; i<4; i++)
      	{
      		eigvec[i][4]=0.0;
      		eigvec[0][i]=cw1[i].r;
      		eigvec[1][i]=cw1[i].i;
      		eigvec[2][i]=cw2[i].r;
      		eigvec[3][i]=cw2[i].i;
      	}
      	eigvec[4][4]=1.0;
      	double ainv[4][4];
      	inver1(ainv);
      	jam441(eigvec[4],ainv[0],a[4]);
// 	printMat("eigen", eigvec);

//                ---------normalization of the eigenvektoren
      	double rn1= eigvec[0][0]*eigvec[1][1]-eigvec[0][1]*eigvec[1][0]
                   +eigvec[0][2]*eigvec[1][3]-eigvec[0][3]*eigvec[1][2];


// 	printf("rn1 = %f\n", rn1);


      	if(rn1 == 0.)  return 1;
      	if(rn1 < 0. )
      	{
		for(int i=0; i<4; i++) eigvec[1][i]=  -eigvec[1][i];
      		clam1.i = -clam1.i;
      	}


      	double sqrn1=sqrt(fabs(rn1));
// 	printf("sqrn1 = %f\n", sqrn1);
      	for(int i=0; i<4; i++)
      	{
      		eigvec[0][i]=eigvec[0][i]/sqrn1;
      		eigvec[1][i]=eigvec[1][i]/sqrn1;
      	}






      	double rn2= eigvec[2][0]*eigvec[3][1]-eigvec[2][1]*eigvec[3][0]
                   +eigvec[2][2]*eigvec[3][3]-eigvec[2][3]*eigvec[3][2];

// 	printf("rn2 = %f\n", rn2);

      	if(rn2 == 0.)  return 2;
      	if(rn2 < 0. )
      	{
		for(int i=0; i<4; i++) eigvec[3][i]=  -eigvec[3][i];
      		clam2.i = -clam2.i;
      	}


      	double sqrn2=sqrt(fabs(rn2));
// 	printf("sqrn2 = %f\n", sqrn2);
      	for(int i=0; i<4; i++)
      	{
      		eigvec[2][i]=eigvec[2][i]/sqrn2;
      		eigvec[3][i]=eigvec[3][i]/sqrn2;
      	}


// 	printMat("eigen2", eigvec);
// 	printMat("matrix", a);
// 	double betax = eigvec[0][0]*eigvec[0][0] + eigvec[1][0]* eigvec[1][0];
// 	double alfax = - (eigvec[0][0]*eigvec[0][1] + eigvec[1][0]* eigvec[1][1]);
// 	double betaz = eigvec[2][2]*eigvec[2][2] + eigvec[3][2]* eigvec[3][2];
// 	double alfaz = - (eigvec[2][2]*eigvec[2][3] + eigvec[3][2]* eigvec[3][3]);
// 	printf("betax,alfax  %f %f\n", betax, alfax);
// 	printf("betaz,alfaz  %f %f\n", betaz, alfaz);

// 	double b[5][5];
// 	jam555(b[0],a[0],eigvec[0]);
// 	double betax2 = b[0][0]*b[0][0] + b[1][0]* b[1][0];
// 	double alfax2 = - (b[0][0]*b[0][1] + b[1][0]* b[1][1]);
// 	double betaz2 = b[2][2]*b[2][2] + b[3][2]* b[3][2];
// 	double alfaz2 = - (b[2][2]*b[2][3] + b[3][2]* b[3][3]);
// 	printf("betax2,alfax  %f %f\n", betax2, alfax2);
// 	printf("betaz2,alfaz  %f %f\n", betaz2, alfaz2);
// 	printMat("eigen3", b);

      	return 0;
}


int Mat::eigenvector(double betax0, double alfax0, double betaz0, double alfaz0, double * dspini)
{
	if(dspini) for(int i=0; i<4; i++) eigvec[4][i]=dspini[i];
	else for(int i=0; i<4; i++) eigvec[4][i]=0.;
	eigvec[4][4] = 1.;


        double wbx=sqrt(betax0);
        eigvec [0][0]=wbx;
        eigvec [1][0]=0.0;
        eigvec [0][1]=-alfax0/wbx;
        eigvec [1][1]=1.0/wbx;
        
        eigvec [0][3]=0.0;
        eigvec [1][3]=0.0;
        eigvec [0][4]=0.0;
        eigvec [1][4]=0.0;
        

        double wbz=sqrt(betaz0);
        eigvec [2][2]=wbz;
        eigvec [3][2]=0.0;
        eigvec [2][3]=-alfaz0/wbz;
        eigvec [3][3]=1.0/wbz;
        
        eigvec [3][0]=0.0;
        eigvec [4][0]=0.0;
        eigvec [3][1]=0.0;
        eigvec [4][1]=0.0;
// 	printMat("eigeni", eigvec);
// 	double betax2 = eigvec[0][0]*eigvec[0][0] + eigvec[1][0]* eigvec[1][0];
// 	double alfax2 = - (eigvec[0][0]*eigvec[0][1] + eigvec[1][0]* eigvec[1][1]);
// 	double betaz2 = eigvec[2][2]*eigvec[2][2] + eigvec[3][2]* eigvec[3][2];
// 	double alfaz2 = - (eigvec[2][2]*eigvec[2][3] + eigvec[3][2]* eigvec[3][3]);
// 	printf("betax,alfax  %f %f\n", betax0, alfax0);
// 	printf("betaz,alfaz  %f %f\n", betaz0, alfaz0);
// 	printf("betax2,alfax  %f %f\n", betax2, alfax2);
// 	printf("betaz2,alfaz  %f %f\n", betaz2, alfaz2);
	return 0;
}

void Mat::printTwiss()
{
	double betax2 = eigvec[0][0]*eigvec[0][0] + eigvec[1][0]* eigvec[1][0];
	double alfax2 = - (eigvec[0][0]*eigvec[0][1] + eigvec[1][0]* eigvec[1][1]);
	double betaz2 = eigvec[2][2]*eigvec[2][2] + eigvec[3][2]* eigvec[3][2];
	double alfaz2 = - (eigvec[2][2]*eigvec[2][3] + eigvec[3][2]* eigvec[3][3]);
	printf("%12.4f  %12.4f  %12.4f  %12.4f\n", betax2, alfax2, betaz2, alfaz2);
}

void Mat::printMat(const char * what,double b[5][5])
{
	printf("%s:\n", what);
	for(int i=0; i<5; i++)
	{
		for(int j=0; j<5; j++) printf(" %f ", b[j][i]);
		printf("\n");
	}
	double dex = (b[0][0]*b[1][1]-b[0][1]*b[1][0]);
	double dez = (b[2][2]*b[3][3]-b[2][3]*b[3][2]);
	printf("det = %f %f\n", dex,dez);
}
void Mat::printA()
{
	for(int i=0; i<5; i++)
	{
		for(int j=0; j<5; j++) printf(" %f ", a[j][i]);
		printf("\n");
	}
}













Element::Element(const char * name, double tl, int ty, double dl, double ml, double st, double b2, int ix)
{
	swn = strdup(name);
	type = ty;
	totl = tl;
	dlength = dl;
	mlength = ml;
	strength = st;
	p2 = p2 = 0.;
	next = NULL;
	qk=qk45=hx=hz=0.;
	this->b2 = b2;
	this->ix = ix;
	if(type == 8)
	{
		hx = -strength;
		lintr3(mlength, 0.,  hx, 0.);
	}
	if(type == 9)
	{
		qk = -strength;
		lintr3(mlength,  qk, 0., 0.);
	}
}



void Element::setStrength(WireUp * w, int step)
{
	if(ix < 0) return;
	double str = type == 12 ? strength : strength*mlength;
	w->dllParms[ix]->set_want_strength(step, str);
}



void Element::getStrength(WireUp * w, int step)
{
	if(ix < 0) return;
	double str = w->dllParms[ix]->want_strength[step];
// 	if(type == 8)
	str /= mlength;
// 	if(strength != 0.)
// 	{
// 		if( fabs(strength - str)  > 1.e-4)
// 		{
// 			double diff = 100. * (strength - str)/strength;
// 			printf("%s: change %f%% strength from %f to %f, l = %f\n", swn, diff,strength, str,mlength);
// 		}
// 	}
// 	else
// 	{
// 		if( fabs(str)  > 1.e-7)
// 		{
// 			printf("%s: invert %f\n", swn, str);
// 		}
// 	}
	if(fabs(strength-str) < 1.e-5) return;
	strength = str;
	if(type == 8) lintr3(mlength, 0.,  -strength, 0.);
	if(type == 9) lintr3(mlength,  -strength, 0., 0.);
}



void Element::doMat(double a[][5])
{
	drift(a, dlength);
	if(type == 8)
		doLintr3(a);
	else
	if(type == 9)
		doLintr3(a);
	else
	if(type == 12)
		thinquad(a,  - strength);
	else
		drift(a, mlength);
}



void Element::doMatInt(Mat & a, double & mocomp, double & chromx, double & chromz, double & qx, double & qz)
{
	double m110 = a.eigvec[0][0];
	double m120 = a.eigvec[1][0];
	double m330 = a.eigvec[2][2];
	double m340 = a.eigvec[3][2];

	drift(a.eigvec, dlength);

	double qsx = a.eigvec[1][0]*m110 - a.eigvec[0][0]*m120 ;
	double qcx = a.eigvec[0][0]*m110 + a.eigvec[1][0]*m120 ;
	double qsz = a.eigvec[3][2]*m330 - a.eigvec[2][2]*m340 ;
	double qcz = a.eigvec[2][2]*m330 + a.eigvec[3][2]*m340 ;
	qx += atan2(qsx,qcx);
	qz += atan2(qsz,qcz);



	betax0 = a.eigvec[0][0]*a.eigvec[0][0] + a.eigvec[1][0]* a.eigvec[1][0];
	alfax0 = - (a.eigvec[0][0]*a.eigvec[0][1] + a.eigvec[1][0]* a.eigvec[1][1]);
	betaz0 = a.eigvec[2][2]*a.eigvec[2][2] + a.eigvec[3][2]* a.eigvec[3][2];
	alfaz0 = - (a.eigvec[2][2]*a.eigvec[2][3] + a.eigvec[3][2]* a.eigvec[3][3]);
	dispx0 = a.eigvec[4][0];
	dispxp0 = a.eigvec[4][1];
	m110 = a.eigvec[0][0];
	m120 = a.eigvec[1][0];
	m330 = a.eigvec[2][2];
	m340 = a.eigvec[3][2];
// 	printf("%-12s  %8.3f  ", "drift", "drift");
// 	printf("%12.4f  %12.4f  %12.4f  %12.4f  %12.4f  %12.4f \n",
// 		betax0, alfax0, betaz0, alfaz0, dispx0, dispxp0);
	if(type == 8)
		doLintr3(a.eigvec);
	else
	if(type == 9)
		doLintr3(a.eigvec);
	else
	if(type == 12)
		thinquad(a.eigvec,  - strength);
	else
		drift(a.eigvec, mlength);
	betax = a.eigvec[0][0]*a.eigvec[0][0] + a.eigvec[1][0]* a.eigvec[1][0];
	alfax = - (a.eigvec[0][0]*a.eigvec[0][1] + a.eigvec[1][0]* a.eigvec[1][1]);
	betaz = a.eigvec[2][2]*a.eigvec[2][2] + a.eigvec[3][2]* a.eigvec[3][2];
	alfaz = - (a.eigvec[2][2]*a.eigvec[2][3] + a.eigvec[3][2]* a.eigvec[3][3]);
	dispx = a.eigvec[4][0];
	dispxp = a.eigvec[4][1];
	betaxm = ( betax + betax0)/2. + (alfax - alfax0)*mlength / 6.;
	betazm = ( betaz + betaz0)/2. + (alfaz - alfaz0)*mlength / 6.;
	dispm = ( a.eigvec[4][0] +dispx0)/2. + (dispxp0 - a.eigvec[4][1])*mlength/12.;
	if(type == 8)
	{
		chromx += b2*dispm*betaxm;
		chromz -= b2*dispm*betazm;
	}
		
	if(type == 10)
	{
		chromx += strength*dispm*betaxm * mlength;
		chromz -= strength*dispm*betazm * mlength;
	}
	qsx = a.eigvec[1][0]*m110 - a.eigvec[0][0]*m120 ;
	qcx = a.eigvec[0][0]*m110 + a.eigvec[1][0]*m120 ;
	qsz = a.eigvec[3][2]*m330 - a.eigvec[2][2]*m340 ;
	qcz = a.eigvec[2][2]*m330 + a.eigvec[3][2]*m340 ;
	qx += atan2(qsx,qcx);
	qz += atan2(qsz,qcz);


	mocomp += hx * dispm *mlength;
	chromx += (qk-2.*hx*hx) * betaxm * mlength;
	chromz -=  qk           * betazm * mlength;

// 	printf("%-12s  %8.3f  ", swn, totl);
// 	printf("%12.4f  %12.4f  %12.4f  %12.4f  %12.4f  %12.4f\n",
// 		betax, alfax, betaz, alfaz, dispx, dispxp);
}



void Element::lintr3(double el, double qk, double hx, double hz)
{

//     x komponente


      	double kx=qk-hx*hx;
      	double sqk=sqrt(fabs(kx));
      	if(sqk > 1.e-10)
      	{
      		double phix=el*sqk;
		double s,c;
      		if(kx < 0. )
    		{	
   			s=sin(phix);
      			c=cos(phix);
      		}
		else
		{
   			s=sinh(phix);
      			c=cosh(phix);
		}
   		mx00=c;
      		mx01=s/sqk;
      		mx02= -hx*(1-c)/(sqk*sqk);
      		mx10= -s*sqk;
      		mx11=c;
      		mx12= -hx*s/sqk;
      		if(kx > 0. )
		{
   			mx02= -mx02;
      			mx10= -mx10;
      		}
	}
	else
	{
//     limes (k-1/rho**2)  gegen 0
      		mx00=1.0;
      		mx01=el;
      		mx02= -hx*el*el/2;
      		mx10=0.0;
      		mx11=1.0;
      		mx12= -hx*el;
	}




//     z komponente


      	double kz= -qk-hz*hz;
      	sqk=sqrt(fabs(kz));
      	if(sqk > 1.e-10)
	{
      		double phiz=el*sqk;
		double s,c;
      		if(kz < 0. )
		{
      			s=sin(phiz);
      			c=cos(phiz);
		}
		else
		{
      			s=sinh(phiz);
      			c=cosh(phiz);
		}
      		mz00=c;
      		mz01=s/sqk;
      		mz02= -hz*(1-c)/(sqk*sqk);
      		mz10= -s*sqk;
      		mz11=c;
      		mz12= -hz*s/sqk;
      		if(kz > 0. )
		{
      			mz02= -mz02;
      			mz10= -mz10;
   		} 
	}
	else
	{
//     limes (k-1/rho**2)  gegen 0

    		mz00=1.0;
      		mz01=el;
      		mz02= -hz*el*el/2;
      		mz10=0.0;
      		mz11=1.0;
      		mz12= -hz*el;
	}
}

void Element::doLintr3(double a[][5])
{

	 double b00 = mx00 * a[0][0] + mx01 * a[0][1];
	 double b10 = mx00 * a[1][0] + mx01 * a[1][1];
	 double b40 = mx00 * a[4][0] + mx01 * a[4][1] + mx02;
	 double b01 = mx10 * a[0][0] + mx11 * a[0][1];
	 double b11 = mx10 * a[1][0] + mx11 * a[1][1];
	 double b41 = mx10 * a[4][0] + mx11 * a[4][1] + mx12;
	 a[0][0] = b00;
	 a[1][0] = b10;
	 a[4][0] = b40;
	 a[0][1] = b01;
	 a[1][1] = b11;
	 a[4][1] = b41;


	 double b22 = mz00 * a[2][2] + mz01 * a[2][3];
	 double b32 = mz00 * a[3][2] + mz01 * a[3][3];
	 double b42 = mz00 * a[4][2] + mz01 * a[4][3] + mz02;
	 double b23 = mz10 * a[2][2] + mz11 * a[2][3];
	 double b33 = mz10 * a[3][2] + mz11 * a[3][3];
	 double b43 = mz10 * a[4][2] + mz11 * a[4][3] + mz12;
	 a[2][2] = b22;
	 a[3][2] = b32;
	 a[4][2] = b42;
	 a[2][3] = b23;
	 a[3][3] = b33;
	 a[4][3] = b43;

}


void Element::thinquad(double a[][5], double k1l)
{

//     x komponente


	 a[0][1] += k1l * a[0][0];
	 a[1][1] += k1l * a[1][0];
	 a[4][1] += k1l * a[4][0];





//     z komponente


	 a[2][3] -= k1l * a[2][2];
	 a[3][3] -= k1l * a[3][2];
	 a[4][3] -= k1l * a[4][2];
}





void Element::drift(double a[][5], double l)
{
	for(int i=0; i<5; i++)
	{
		a[i][0] += l * a[i][1];
		a[i][2] += l * a[i][3];
	}
}




Opticalc::Opticalc(const char * file, WireUp * w)
{
	FILE * fp = fopen(file,"r");
	if(!fp)
	{
		fprintf(stderr,"Cant open file %s\n", file);
		ok = errno;
		return;
	}

	last = chain = new Element("start", tl, 12, dl,0.,0.,0.,  -1);
	tl=0;
	dl = 0.;
	double r;
	nelem = 1;

	while(fscanf(fp,"%d%s%s%s%lf%lf%lf", &d,latnam,swn,nam, &s, &l, &r) == 7)
	{
		if( !strcmp(swn,"none") ) strcpy(swn,latnam);
		tl += l;
		int ix = w->name2Index(swn);
		if(d == 8) // sbend
		{
			double b2 = 0.;
			if( !strcmp(latnam, "d") ) b2 = -0.0709974248207;
			last->next = new Element(swn, tl, d, dl, l,s/l, b2,  ix);
			nelem++;
			last = last->next;
			dl = 0.;
		}
		else
		if(d == 9 || d == 10) // quadrupole or sextupole
		{
			last->next = new Element(swn, tl, d, dl, l,s, 0., ix);
			nelem++;
			last = last->next;
			dl = 0.;
		}
		else
		if(d == 12 && ! strcmp(nam,"k1l") ) // quadrupole
		{
			last->next = new Element(swn, tl, d, dl, l,s, 0., ix);
			nelem++;
			last = last->next;
			dl = 0.;
		}
		else
		{
			dl += l;
		}

	}
	last->next = new Element("ende", tl, 12, dl,0.,0.,0.,  -1);
	nelem++;
	fclose(fp);

	line = new Element*[nelem];
	int ielem = 0;
	last=chain;
	while(last)
	{
		line[ielem] = last;
		last = last->next;
	}



	ok = 0;
}

void Opticalc::setStrength(WireUp * w, int step)
{
	if(!w) return;
	last=chain;
	while(last)
	{
		last->setStrength(w,step);
		last = last->next;
	}
}

void Opticalc::calc(WireUp * w, int step)
{
	last=chain;
	a.unitA();
	while(last)
	{
		if(w) last->getStrength(w,step);
		last->doMat(a.a);
		last = last->next;
	}

// 	a.printA();
	a.eigenwerte();
	a.eigenvector();
// 	 printf("%f\n", tl);
// 	double betax= 3.00;
// 	double betaz= 3.00;
// 	double alfax = 0.00;
// 	double alfaz = 0.00;
// 	a.eigenvector(betax,alfax,betaz,alfaz);

	mocomp = 0.;
	chromx = 0.;
	chromz = 0.;
	qx = 0.;
	qz = 0.;
	last=chain;
	while(last)
	{
		last->doMatInt(a,mocomp, chromx, chromz, qx, qz);
		last = last->next;
	}
	chromx /= (4.*M_PI);
	chromz /= (4.*M_PI);
	qx /= (2.*M_PI);
	qz /= (2.*M_PI);
	mocomp /= tl;
	double gammat = 1./sqrt(fabs(mocomp));
	printf(" chrom = %f %f moc %f gammat %f qx %f qz %f\n", chromx,chromz,mocomp,gammat, qx, qz);

}
