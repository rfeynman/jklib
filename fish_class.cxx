#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include "fish_class.hxx"

// 	The doubles come in the fortran format with line 1D1 insteadt of 1e1.
double FISH::readDbl()
{
	char line[100];
	double d;

	if( fscanf(fp,"%s", line) != 1)
	{
		fprintf(stderr, "Could not read the next number\n");
		exit(-1);
	}
	for(char * p = line; *p; p++) if(*p == 'D') *p = 'e';
	if( sscanf(line,"%le", &d) != 1)
	{
		fprintf(stderr, "Could not convert the next number\n");
		exit(-1);
	}
	return d;
}

FISH::FISH(const char cavName[])
{
	nfiles=0;
	char line[301];

// the T& files are numbered 1,2,3.. or 01,02,03.. 
//  We first have to find out which one it is.
	char testfilename1[100];
	char testfilename2[100];
	sprintf(testfilename1,"%s%d.T7", cavName, 1);
	fp=fopen(testfilename1, "r");
	fprintf(stderr, "test open file %s\n", testfilename1);
	if(fp) oldnames=0; 
	else 
	{
		sprintf(testfilename2,"%s%02d.T7", cavName, 1);
		fprintf(stderr, "test open file %s\n", testfilename2);
		fp=fopen(testfilename2, "r");
		if(!fp)
		{
			fprintf(stderr, "Cant open fish files\n");
			exit(-1);
		}
		oldnames=1;
	}
	fprintf(stderr,"oldnames = %d\n", oldnames);
	fclose(fp);


	for(int n=0; n<tausend; n++)
	{
		if(oldnames) sprintf(files[n],"%s%2d.T7", cavName, n+1);
		else        sprintf(files[n],"%s%d.T7", cavName, n+1);
		for(int j=0; j<strlen(files[n]); j++) if(files[n][j] == ' ') files[n][j]='0';
		fprintf(stderr,"open %s\n", files[n]);

		fp=fopen(files[n], "r");
		if(!(fp))
		{
			if(errno != ENOENT)
			{
				fprintf(stderr, " cant open file %s\n", files[n]);
				perror(">");
			}
			break;
		}
		fclose(fp);
		nfiles++;
	}

	if(!nfiles)
	{
		fprintf(stderr, " could not open any files\n");
		exit(-1);
	}
	else
	{
		fprintf(stderr, "%d filesfound\n", nfiles);
	}



	rmina = new double[nfiles];
 	rmaxa = new double[nfiles];
 	zmina = new double[nfiles];
 	zmaxa = new double[nfiles];
 	freqa = new double[nfiles];
	rstepa= new int[nfiles];
	zstepa= new int[nfiles];

	// pass 1: read the headers to compute the array sizes
	for(int n=0; n<nfiles; n++)
	{
		fp=fopen(files[n], "r");
		if(!fp)
		{
			fprintf(stderr, " cant open file %s\n", files[n]);
			break;
		}
		zmina[n] = readDbl()/100.;	// convert from cm to m
		zmaxa[n] = readDbl()/100.;
		fscanf(fp, "%d", zstepa+n);
		fprintf(stderr,"%f %lf %d\n", zmina[n], zmaxa[n], zstepa[n]);
		freqa[n] = readDbl();
		fprintf(stderr,"%f\n", freqa[n]);
		rmina[n] = readDbl()/100.;
		rmaxa[n] = readDbl()/100.;
		fscanf(fp, "%d", rstepa+n);
		fprintf(stderr,"%f %lf %d\n", rmina[n], rmaxa[n], rstepa[n]);

		if(n)
		{
			if( zmax != zmina[n] )
			{
				fprintf(stderr, "zmin does not connect %d: %f %f\n", n, zmax, zmina[n]);
				exit(-1);
			}
			if(  rmax != rmaxa[n])
			{
				fprintf(stderr, "rmax does not connect %d\n", n);
				exit(-1);
			}
			if(  rmin != rmina[n])
			{
				fprintf(stderr, "rmin does not connect %d\n", n);
				exit(-1);
			}
			if(  rstep != rstepa[n])
			{
				fprintf(stderr, "rstep does not connect %d\n", n);
				exit(-1);
			}
			if(  freq != freqa[n])
			{
				fprintf(stderr, "freq does not connect %d: %f %f\n", n, freq, freqa[n]);
				exit(-1);
			}
			zmax = zmaxa[n];
			zstep += zstepa[n]-1;
		}
		else
		{
			zmax = zmaxa[n];
			zmin = zmina[n];
			rmax = rmaxa[n];
			rmin = rmina[n];
			freq = freqa[n];
			zstep = zstepa[n];
			rstep = rstepa[n];
			dr=(rmax-rmin)/(rstep-1);
		}
		fclose(fp);
	}

	fprintf(stderr, "zstep %d  zmax %f rstep %d\n", zstep, zmax, rstep);


	zloc=new double[zstep+1];
	rloc=new double[rstep+1];
	ezero=new double[zstep+1];

	ez = new double*[rstep+1];
	er = new double*[rstep+1];
	ee = new double*[rstep+1];
	bt = new double*[rstep+1];

	for(int r=0; r<=rstep; r++)
	{
		ez[r] = new double[zstep+2];
		er[r] = new double[zstep+2];
		ee[r] = new double[zstep+2];
		bt[r] = new double[zstep+2];
		rloc[r]= rmin+ r*dr; 
	}


	int k=0;
	int zz=0;
	Emax=0;
	Ezmax=0;
	Ermax=0;
	Bmax=0;
	// pass 2: read the data
	for(int n=0; n<nfiles; n++)
	{
		fprintf(stderr, "open file %s\n", files[n]);
		fp=fopen(files[n], "r");
		if(!fp)
		{
			fprintf(stderr, " cant open file %s\n", files[n]);
			break;
		}

		fgets(line, 300, fp);	// skip the headers
		fgets(line, 300, fp);
		fgets(line, 300, fp);
		for(int r=0; r<= rstepa[n]; r++)
		{
			for(int z=zz; z<= zstepa[n]+zz; z++)
			{
// 				fprintf(stderr, "r %d z %d\n", r, z);

				ez[r][z] = readDbl();
				er[r][z] = readDbl();
				ee[r][z] = readDbl();
				bt[r][z] = readDbl();
				if(fabs(ee[r][z]) > Emax)   Emax  = fabs(ee[r][z]);
				if(fabs(er[r][z]) > Ermax)  Ermax = fabs(er[r][z]);
				if(fabs(ez[r][z]) > Ezmax)  Ezmax = fabs(ez[r][z]);
				if(fabs(bt[r][z]) > Bmax)   Bmax  = fabs(bt[r][z]);
			}
		}
		for(int z=zz; z<= zstepa[n]+zz; z++)
		{
			double dz=(zmaxa[n]-zmina[n])/(zstepa[n]-1);
			zloc[z]= zmina[n]+ (z-zz)*dz; 
			ezero[z]=ez[0][z];
			if(z >= zstep)
			printf("%d %e %e\n", z,  zloc[z], ezero[z]); 
		}

		zz += zstepa[n]-1;
		fclose(fp);
	}
	dz=(zmax-zmin)/(zstep-1);

}

FISH::~FISH()
{

	delete [] rmina ;
 	delete [] rmaxa ;
 	delete [] zmina ;
 	delete [] zmaxa ;
 	delete [] freqa ;
	delete [] rstepa;
	delete [] zstepa;
	delete [] zloc;
	delete [] ezero;

	for(int r=0; r<=rstep; r++)
	{
		delete [] ez[r];
		delete [] er[r];
		delete [] ee[r];
		delete [] bt[r];
	}
	delete [] ez;
	delete [] er;
	delete [] ee;
	delete [] bt;

}
