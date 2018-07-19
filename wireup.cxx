#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "jkLib/wireup.hxx"


class Magnet  **WireUp::mag = NULL;
int WireUp::nmag =0;
class Power   **WireUp::pps = NULL;
int WireUp::npps =0;
class Wires **WireUp::kabel = NULL;
int WireUp::nwires =0;
class Fields  **WireUp::ff = NULL;
int WireUp::nff =0;
char WireUp::current_species[15];
int WireUp::q_ags;
int WireUp::atom_A;
int WireUp::atom_Z;
double WireUp::mass_per_u;
int WireUp::nspecies = 0;
Speclist * WireUp::specs = NULL;

const double WireUp::clight = 2.997924562e8;

class RfParm ** WireUp::rfparm;
int WireUp::nrfparm;
class Cavity ** WireUp::rfcavity;
int WireUp::nrfcavity;
class CWire1 ** WireUp::rfcwire1;
int WireUp::nrfcwire1;
class CWire2 ** WireUp::rfcwire2;
int WireUp::nrfcwire2;
class CRadius** WireUp::rfcradius;
int WireUp::nrfcradius;
class CFreq** WireUp::rfcfreq;
int WireUp::nrfcfreq;

class Parameter ** WireUp::allParms;
int WireUp::nparms;


/*------- utilities  --------------------------------------------*/
 
char * strDup(const char * s)
{
	int l=strlen(s)+1;
	char * r = new char[l];
	memcpy(r,s,l);
	return r;
}

double rdd(FILE * f)
{
	double i;
	char x[21];
	if( fscanf(f,"%lf",&i) != 1) 
	{
		fscanf(f,"%20s",x);
		printf("no double: %s\n",x);
		exit(1);
	}
	return(i);
}
 
int rdi(FILE * f)
{
	int i;
	char x[21];
	if( fscanf(f,"%d",&i) != 1) 
	{
		fscanf(f,"%20s",x);
		printf("no int: %s\n", x);
		exit(1);
	}
	return(i);
}
 
/*------- wireup class  --------------------------------------------*/


WireUp::WireUp(const char * file)
{
	if(mag) 
	{
		fprintf(stderr,"There can only be one wireUp master!");
		exit(-1);
	}

	nstep=0;
	nwires = 0;
	jumpGamma = 24.;

	FILE * fp = fopen(file,"r");
	if(!fp)
	{
	  fprintf(stderr,"error: could not open file: %s \n",file);
	  exit(1);
	}
 
	double version;
	char section[100];
	char s1[100];
	fgets(section,100,fp);
	printf("%s\n",section);
	sscanf(section,"%s%s%lf",s1,s1,&version);
	if( int(version) != 8)
	{
		fprintf(stderr,"error: wrong version of the wireup file!\n");
		exit(-1);
	}

	fgets(section,100,fp);
	printf("%s\n",section);
 
 
        fscanf(fp,"%s",section);
        printf("section: %s\n",section);
        nspecies = rdi(fp);
        specs= new Speclist[nspecies];
        fscanf(fp,"%s",current_species);

        for(int i=0; i<nspecies; i++)
        {
                fscanf(fp,"%s",specs[i].species);
                specs[i].mass_per_u = rdd(fp);
                specs[i].atom_A = rdi(fp);
                specs[i].atom_Z = rdi(fp);
                specs[i].q_ags  = rdi(fp);
                if( ! strcmp(current_species, specs[i].species))
                {
                        mass_per_u = specs[i].mass_per_u;
                        atom_A = specs[i].atom_A;
                        atom_Z = specs[i].atom_Z;
                        q_ags  = specs[i].q_ags;
                }
        }


 

 
 
/*------- set up foils ------------------------------------------*/
	fscanf(fp,"%s",section);
	printf("section: %s\n",section);
	int nfoils = rdi(fp);
	for(int n=0; n<nfoils; n++)
	{
		char inst[40];
		fscanf(fp,"%s",inst);
	}
 
 
/*------- read the Field information ----------------------------*/
 
	fscanf(fp,"%s",section);
	printf("section: %s\n",section);
	nff = rdi(fp);
	ff = new Fields*[nff];
	for(int n=0; n< nff; n++)
	{
		int k=rdi(fp);

		ff[n] = new Fields(k);
		for(int j=0; j<k; j++)
		{
			ff[n]->current[j]  =rdd(fp);
			ff[n]->int_field[j]=rdd(fp);
		}
		ff[n]->make_spline();
	}
 
/*------- set up magnets ----------------------------------------*/
	fscanf(fp,"%s",section);
	printf("section: %s\n",section);
	nmag = rdi(fp);
	mag = new Magnet*[nmag];
	dag = NULL;
	dps = NULL;
	for(int n=0; n<nmag; n++)
	{
		char inst[40];
		fscanf(fp,"%s",inst);
		int lattice_index = rdi(fp);
		int field_index   = rdi(fp);
		int foil_index    = rdi(fp);
		int type_index    = rdi(fp);
		int pebble_index  = rdi(fp);
		mag[n]= new Magnet(inst, this, lattice_index,
			field_index, foil_index,
			type_index, pebble_index, n);
	}
 
/*------- set up power supplies----------------------------------*/
	fscanf(fp,"%s",section);
	printf("section: %s\n",section);
	npps = rdi(fp);
	pps = new Power*[npps];
	for(int n=0; n<npps; n++)
	{
		char inst[40], model[40], polar[40];
		fscanf(fp,"%s",inst);
		fscanf(fp,"%s",model);
		int christie = ! strcmp(model,"Christie");
		fscanf(fp,"%s",polar);
                int sub1               = rdi(fp);
                int sub2               = rdi(fp);
                int mag1               = rdi(fp);
                int mag2               = rdi(fp);
                int pol1               = rdi(fp);

		int init               = rdi(fp);
		double max_current     = rdd(fp);
		double min_current     = rdd(fp);
		double max_advice      = rdd(fp);
		double min_advice      = rdd(fp);
		double i_rating        = rdd(fp);
		double v_rating        = rdd(fp);
		double bit_per_amp_rdb = rdd(fp);
		double mvolts          = rdd(fp);
		double bit_per_amp     = rdd(fp);
		double max_ramp_speed  = rdd(fp);
		double resistance      = rdd(fp);
		double inductance      = rdd(fp);
                // big cludge------------------------------------------
                if(resistance > inductance) inductance = 2.* resistance;
                // big cludge------------------------------------------
		pps[n]= new Power(inst, init, max_current, min_current,
			max_advice, min_advice, i_rating, v_rating,
			bit_per_amp_rdb, mvolts, bit_per_amp,
			max_ramp_speed, resistance, inductance, christie);
	}
 
 
 
/*------- set up families ---------------------------------------*/
	fscanf(fp,"%s",section);
	printf("section: %s\n",section);
	int nfamily = rdi(fp);
	for(int j=0; j<nfamily; j++)
	{
		int count = rdi(fp);
		int namx  = rdi(fp);
		mag[namx]->fam = new int[count+1];
		mag[namx]->nfam=count+1;
		for(int i=0; i<count; i++)
		{
			int k = rdi(fp);
			int p = rdi(fp);
			mag[namx]->fam[i] = k;
			mag[k]->family_head=namx;
			mag[k]->family_polar=p;
		}
		mag[namx]->fam[count] = namx;
	}
	for(int n=0; n< nmag; n++)
	{
		if(mag[n]->family_head <0 && !(mag[n]->fam) )
		{
			mag[n]->fam = new int[1];
			mag[n]->fam[0] = n;
			mag[n]->nfam = 1;
		}
	}

 
 
/*------- set up power Wires ------------------------------------*/
 
	fscanf(fp,"%s",section);
	printf("section: %s\n",section);
	kabel = new Wires*[npps];  // wild guess how many
	for(int n=rdi(fp); n > 0; n=rdi(fp) )
	{
		Wires *a = new Wires(n);
		kabel[nwires++] = a;
		for(int i=0; i<n; i++) a->mags[i] = rdi(fp);
		for(int i=0; i<n; i++) a->ppss[i] = rdi(fp);
		a->mag_con = new Connection*[n];
		for(int i=0; i<n; i++)
		{
			int count = rdi(fp);
			a->mag_con[i] = new Connection(count);
			for(int j=0; j<count; j++)
			{
				a->mag_con[i]->index[j]    = rdi(fp);
				a->mag_con[i]->polarity[j] = rdi(fp);
			}
		}

		a->pps_con = new Connection*[n];
		for(int i=0; i<n; i++)
		{
			int count = rdi(fp);
			a->pps_con[i] = new Connection(count);
			for(int j=0; j<count; j++)
			{
				a->pps_con[i]->index[j]    = rdi(fp);
				a->pps_con[i]->polarity[j] = rdi(fp);
			}
		}


 
		for(int i=0; i<n; i++) mag[a->mags[i]]->mat=a;
		for(int i=0; i<n; i++) pps[a->ppss[i]]->mat=a;
	}
//      for(int n=0; n<nwires; n++)
//      {
//              Wires *a = kabel[n];
//              int count =a->n;
//              printf("kabel[%d] count %d    ", n, count);
//              printf("mags= ");
//              for(int i=0; i<count;i++) printf(" %d ", a->mags[i]);
//              printf("   ppss= ");
//              for(int i=0; i<count;i++) printf(" %d ", a->ppss[i]);
//              printf("\n");
//      }
 

/*------- set up RF Cavity ---------------------------------------*/
	fscanf(fp,"%s",section);
	printf("section: %s\n",section);
	nrfcavity = rdi(fp);
	rfcavity = new Cavity*[nrfcavity];
	for(int j=0; j<nrfcavity; j++)
	{
		char inst[40];
		fscanf(fp,"%s",inst);
		int n1 = rdi(fp);
		int n2 = rdi(fp);
		int n3 = rdi(fp);
		rfcavity[j] = new Cavity(inst);
	}



/*------- end of input -------------------------------------------*/
	fscanf(fp,"%s",section);
	printf("section: %s\n",section);

	fclose(fp);

/*------- make namelookup m---------------------------------------*/

        nparms = nmag + nrfparm;
        allParms = new Parameter*[nparms];
        for(int i=0; i<nmag; i++)
		allParms[i] = (Parameter *) mag[i];
        for(int i=0; i<nrfparm; i++)
		allParms[i+nmag] = (Parameter *) rfparm[i];
      
	// sort 
        for(int i=0; i<nparms-1; i++)
        for(int j=i+1; j<nparms; j++)
        {
                int g = strcmp( allParms[i]->swn, allParms[j]->swn );
                if( g > 0)
                {
                        Parameter * tmp = allParms[i];
                        allParms[i]= allParms[j];
                        allParms[j]=tmp;
                }
	}      
}





WireUp::WireUp()
{
	if(!mag)
	{
		fprintf(stderr,"no master wireUp\n");
		exit(-1);
	}
	this->nstep=0;
	gamma_e = NULL;
	gamma_e_p = NULL;
	konstante = NULL;
	konstante_p = NULL;

	dag = NULL;
	dps = NULL;
	dfparm = NULL;
	dfcavity = NULL;
	dllParms = NULL;
	jumpGamma = 24.;
}


int WireUp::wInit(int nstep)
{
	if(!mag)
	{
		fprintf(stderr,"no master wireUp\n");
		exit(-1);
	}
	this->nstep=nstep;
	gamma_e = new double[nstep];
	gamma_e_p = new double[nstep];
	konstante = new double[nstep];
	konstante_p = new double[nstep];

	dag = new MagnetD*[nmag];
	for(int n=0; n<nmag; n++) dag[n] = new MagnetD(mag[n], this);

	dps = new PowerD*[npps];
	for(int n=0; n<npps; n++) dps[n] = new PowerD(pps[n], this);

	dfparm = new RfParmD*[nrfparm];
	for(int n=0; n<nrfparm; n++) dfparm[n] = new RfParmD(rfparm[n], this);

	dfcavity = new CavityD*[nrfcavity];
	for(int n=0; n<nrfcavity; n++) dfcavity[n] = new CavityD(rfcavity[n], this);

	dllParms = new ParameterD*[nparms];
	for(int n=0; n<nparms; n++)
	{
              	int j = allParms[n]->myindex;
       		switch(allParms[n]->mytype)
       		{
               		case 0:
				dllParms[n] = dag[j];
                       		break;
               		case 1:
				dllParms[n] = dfparm[j];
                       		break;
               		default:
				dllParms[n] = NULL;
                       		break;
       		}
	}
	return 0;
}


WireUp::~WireUp()
{
	if(dag)
	{
		delete [] gamma_e;
		delete [] gamma_e_p;
		delete [] konstante;
		delete [] konstante_p;
		for(int n=0; n<nmag; n++) delete dag[n];
		delete [] dag;
		for(int n=0; n<npps; n++) delete dps[n];
		delete [] dps;
		for(int n=0; n<nrfparm; n++) delete dfparm[n];
		delete [] dfparm;
		for(int n=0; n<nrfcavity; n++) delete dfcavity[n];
		delete [] dfcavity;
		delete [] dllParms;
	}
	else
	{
		for(int n = 0; n<nwires; n++) delete kabel[n];
		delete [] kabel;
	
		for(int n=0; n< nff; n++) delete ff[n];
		delete [] ff;

		for(int n=0; n<npps; n++) delete pps[n];
		delete [] pps;
		for(int n=0; n<nmag; n++) delete mag[n];
		delete [] mag;

		for(int n=0; n<nrfparm; n++) delete rfparm[n];
		delete [] rfparm;
		for(int n=0; n<nrfcavity; n++) delete rfcavity[n];
		delete [] rfcavity;
		for(int n=0; n<nrfcwire1; n++) delete rfcwire1[n];
		delete [] rfcwire1;
		for(int n=0; n<nrfcwire2; n++) delete rfcwire2[n];
		delete [] rfcwire2;
	}
}


void WireUp::set_gamma(int step, double g)
{
	gamma_e[step] = g;
	set_konstante(step);
}

void WireUp::set_jumpGamma(double g)
{
	jumpGamma = g;
}


void WireUp::set_konstante(int st)
{
	konstante[st] = atom_A * mass_per_u
		* sqrt(gamma_e[st]*gamma_e[st]-1.) / (clight*1.e-9);
	// 10e-9 because integrated field is given in tesla, not gauss
}
 
double WireUp::Konstante(int st)
{
	return konstante[st];
}
 
double WireUp::Konstmax()
{
	return konstmax;
}


void WireUp::clear_all()
{
	for(int step=0; step<nstep; step++)
		clear_all(step);
}



void WireUp::clear_all(int step)
{
	for(int n=0; n<nmag; n++) dag[n]->clear_all(step);
	for(int n=0; n<nstep; n++)dps[n]->clear_all(step);
	for(int n=0; n<nrfparm; n++) dfparm[n]->clear_all(step);
	for(int n=0; n<nrfcavity; n++) dfcavity[n]->clear_all(step);
}


void MagnetD::clear_all(int step)
{
	want_strength[step] = 0.;
	trim_strength[step] = 0.;
	summ_strength[step] = 0.;
	want_current[step] = 0.;
	trim_current[step] = 0.;
	summ_current[step] = 0.;
}

void PowerD::clear_all(int step)
{
	current[step]=0.;
	currentp[step]=0.;
}

void ParameterD::clear_all(int step)
{
	want_strength[step]=0.;
	trim_strength[step]=0.;
	summ_strength[step]=0.;
}

void CavityD::clear_all(int step)
{
	voltage[step]=0.;
}

int WireUp::name2Index(const char * name)
{
	int high =nparms, low=0;

	while(high > low)
	{
		int mid=low +(high - low)/2;
		int r = strcmp(name, allParms[mid]->swn);
		if(r == 0) return mid;
		if(r < 0)
		{
			high = mid;
		}
		else
		{
			low = mid + 1;
		}
	}
	return -1;
}


int WireUp::psName2Index(const char * name)
{
	int high =npps, low=0;

	while(high > low)
	{
		int mid=low +(high - low)/2;
		int r = strcmp(name, pps[mid]->swn);
		if(r == 0) return mid;
		if(r < 0)
		{
			high = mid;
		}
		else
		{
			low = mid + 1;
		}
	}
	return -1;
}



int WireUp::set_want(int step, const char * swn, double want)
{
	int i = name2Index(swn);
//      printf("index = %d\n",i);
	if(i >= 0)
	{
		dllParms[i]->set_want_strength(step, want);
		return 0;
	}
	return -1;
}


int WireUp::set_trim(int step, const char * swn, double trim)
{
	int i = name2Index(swn);
	if(i >= 0)
	{
		dllParms[i]->set_trim_strength(step, trim);
		return 0;
	}
	return -1;
}


int WireUp::set_wantNtrim(int step, const char * swn, double want, double trim)
{
	int i = name2Index(swn);

	if(i >= 0)
	{
		dllParms[i]->set_want_strength(step, want);
		dllParms[i]->set_trim_strength(step, trim);
		return 0;
	}
	return -1;
}


int WireUp::set_psCurrent(int step, const char * swn, double c)
{
	int i = psName2Index(swn);

	if(i >= 0)
	{
		dps[i]->set_current(step, c);
		return 0;
	}
	return -1;
}


void WireUp::make_splines(int n, double * x)
{
	for(int i=0; i<nmag; i++) dag[i]->make_spline(n,x);
	for(int i=0; i<nrfparm; i++) dfparm[i]->make_spline(n,x);
}

void ParameterD::make_spline(int n, double * x)
{
	delete trim_spline;
	trim_spline = new Spline(x, trim_strength, n,
		0., 0.,SPLINE_ARRAY_USE);
	delete want_spline;
	want_spline = new Spline(x, want_strength, n,
		0., 0.,SPLINE_ARRAY_USE);
}



void WireUp::fill_step(int step, double p)
{
	for(int i=0; i<nmag; i++) dag[i]->fill_step(step,p);
	for(int i=0; i<nrfparm; i++) dfparm[i]->fill_step(step,p);
}


void MagnetD::fill_step(int step, double p)
{
	double want = want_spline->lookup(p);
	double trim = trim_spline->lookup(p);
	set_want_strength(step, want);
	set_trim_strength(step, trim);
}


void ParameterD::fill_step(int step, double p)
{
	want_strength[step] = want_spline->lookup(p);
	trim_strength[step] = trim_spline->lookup(p);
}







void WireUp::copy_pebble(int t, int f, int pebble)
{
	for(int m=0; m<nmag; m++)
	{
		if(mag[m]->family_head < 0 && mag[m]->pebble_index == pebble)
		{
			dag[m]->trim_strength[t] = dag[m]->trim_strength[f];
			dag[m]->set_want_strength(t, dag[m]->want_strength[f]);
		}
	}
	for(int m=0; m<nrfparm; m++)
	{
		if(rfparm[m]->pebble_index == pebble)
		{
			dfparm[m]->want_strength[t] = dfparm[m]->want_strength[f];
			dfparm[m]->trim_strength[t] = dfparm[m]->trim_strength[f];
			dfparm[m]->summ_strength[t] = dfparm[m]->summ_strength[f];
		}
	}

}

void WireUp::copy_step(int t, int f)
{
	for(int m=0; m<nmag; m++) dag[m]->copy_step(t, f);
	for(int m=0; m<npps; m++) dps[m]->copy_step(t, f);
	for(int m=0; m<nrfparm; m++) dfparm[m]->copy_step(t, f);
	for(int m=0; m<nrfcavity; m++) dfcavity[m]->copy_step(t, f);

	gamma_e[t] =  gamma_e[f];
	gamma_e_p[t] =  gamma_e_p[f];
	konstante[t] =  konstante[f];
	konstante_p[t] =  konstante_p[f];
}

void MagnetD::copy_step(int t, int f)
{
	want_strength[t] =  want_strength[f];
	trim_strength[t] =  trim_strength[f];
	summ_strength[t] =  summ_strength[f];
	want_current[t] =  want_current[f];
	trim_current[t] =  trim_current[f];
	summ_current[t] =  summ_current[f];
}


void PowerD::copy_step(int t, int f)
{
	current[t] =  current[f];
	currentp[t] =  currentp[f];
}


void ParameterD::copy_step(int t, int f)
{
		want_strength[t] = want_strength[f];
		trim_strength[t] = trim_strength[f];
		summ_strength[t] = summ_strength[f];
}


void CavityD::copy_step(int t, int f)
{
		voltage[t] = voltage[f];
}





void WireUp::insert_step(int i, int n)
{
// insert empty step i shifting steps i+1 to n-1
	for(int m=0; m<nmag; m++) dag[m]->insert_step(i, n);
	for(int m=0; m<npps; m++) dps[m]->insert_step(i, n);
	for(int m=0; m<nrfparm; m++) dfparm[m]->insert_step(i,n);
	for(int m=0; m<nrfcavity; m++) dfcavity[m]->insert_step(i, n);

	for(int j= n;j>i; j--)
	{
		int m=j-1;      
		gamma_e[j] =  gamma_e[m];
		gamma_e_p[j] =  gamma_e_p[m];
		konstante[j] =  konstante[m];
		konstante_p[j] =  konstante_p[m];
	}
}


void MagnetD::insert_step(int i, int n)
{
	for(int j= n;j>i; j--)
	{
		int m=j-1;      
		want_strength[j] =  want_strength[m];
		trim_strength[j] =  trim_strength[m];
		summ_strength[j] =  summ_strength[m];
		want_current[j] =  want_current[m];
		trim_current[j] =  trim_current[m];
		summ_current[j] =  summ_current[m];
	}
}


void PowerD::insert_step(int i, int n)
{
	for(int j= n;j>i; j--)
	{
		int m=j-1;      
		current[j] =  current[m];
		currentp[j] =  currentp[m];
	}
}



void ParameterD::insert_step(int i, int n)
{
	for(int j= n;j>i; j--)
	{
		int m=j-1;      
		want_strength[j] = want_strength[m];
		trim_strength[j] = trim_strength[m];
		summ_strength[j] = summ_strength[m];
	}
}


void CavityD::insert_step(int i, int n)
{
	for(int j= n;j>i; j--)
	{
		int m=j-1;      
		voltage[j] = voltage[m];
	}
}








void WireUp::delete_step(int i, int n)
{
	for(int m=0; m<nmag; m++) dag[m]->delete_step(i, n);
	for(int m=0; m<npps; m++) dps[m]->delete_step(i, n);
	for(int m=0; m<nrfparm; m++) dfparm[m]->delete_step(i,n);
	for(int m=0; m<nrfcavity; m++) dfcavity[m]->delete_step(i, n);

	for(int j= i;j<n-1; j++)
	{
		int m=j+1;      
		gamma_e[j] =  gamma_e[m];
		gamma_e_p[j] =  gamma_e_p[m];
		konstante[j] =  konstante[m];
		konstante_p[j] =  konstante_p[m];
	}

}


void MagnetD::delete_step(int i, int n)
{
	for(int j= i;j<n-1; j++)
	{
		int m=j+1;      
		want_strength[j] =  want_strength[m];
		trim_strength[j] =  trim_strength[m];
		summ_strength[j] =  summ_strength[m];
		want_current[j] =  want_current[m];
		trim_current[j] =  trim_current[m];
		summ_current[j] =  summ_current[m];
	}
}


void PowerD::delete_step(int i, int n)
{
	for(int j= i;j<n-1; j++)
	{
		int m=j+1;      
		current[j] =  current[m];
		currentp[j] =  currentp[m];
	}
}


void ParameterD::delete_step(int i, int n)
{
	for(int j= i;j<n-1; j++)
	{
		int m=j+1;      
		want_strength[j] = want_strength[m];
		trim_strength[j] = trim_strength[m];
		summ_strength[j] = summ_strength[m];
	}
}


void CavityD::delete_step(int i, int n)
{
	for(int j= i;j<n-1; j++)
	{
		int m=j+1;      
		voltage[j] = voltage[m];
	}
}












void WireUp::get_indexArray(int all, int pebble, IndexArray& a)
{
	int i=0;
	for(int n=0; n<nparms; n++) i += allParms[n]->is_pebble(all, pebble);

	a.a = new int[i];
	a.t = new int[i];
	a.names = new char*[i+1];
	i=0;
	for(int n=0; n<nparms; n++)
	{
		if(allParms[n]->is_pebble(all, pebble))
	        {
			a.a[i]=allParms[n]->myindex;
			a.t[i]=allParms[n]->mytype;
			a.names[i++]=allParms[n]->swn;
	      	}
	}
	a.names[i] = NULL;
	a.n=i;
}




void WireUp::readPebble(int file, int pebble, int step)
{
	int nn;
	int bytes = read(file, &nn, sizeof(int));
	if(bytes != sizeof(int))
	{
		printf("readerror\n");
	}

	StrengthItem *s = new StrengthItem[nn];
	bytes = read(file, s, nn*sizeof(StrengthItem));
	if(bytes != nn*(int)sizeof(StrengthItem))
	{
		printf("readerror\n");
	}

        int i=0, n=0;
        int nmis=0; int * missing = new int[nmag];
        while(n<nmag && i < nn)
	{
// 		printf("readp %d %s %d %s\n", i, s[i].swn, n, allParms[n]->swn);

		if(allParms[n]->is_pebble(0, pebble) )
		{
			int cmp = strcmp(s[i].swn, allParms[n]->swn);
			if(! cmp)
			{
                        	dllParms[n]->set_want_strength(step, s[i].want);
                        	dllParms[n]->set_trim_strength(step, s[i].trim);
				i++; n++;
			}
			else if( cmp < 0)
			{
				printf("no such magnet: %s %f %f\n",
					s[i].swn, s[i].want, s[i].trim);
				i++;
			}
			else
			{
				if(pebble < 6) missing[nmis++]=n;
				printf("magnet %s set to zero\n", mag[n]->swn);
                        	dllParms[n]->set_want_strength(step, 0.0);
                        	dllParms[n]->set_trim_strength(step, 0.0);
				n++;
			}
		}
		else n++;
	}
        while(n < nmag)
        {
                if(mag[n]->family_head <0 && pebble == mag[n]->pebble_index)
                {
                        missing[nmis++]=n;
                        dllParms[n]->set_want_strength(step, 0.0);
                        dllParms[n]->set_trim_strength(step, 0.0);
                        printf("magnet %s set to zero\n", mag[n]->swn);
                }
                n++;   
        }
        while(i < nn)
        {
                printf("no such magnet: %s %f %f\n", s[i].swn, s[i].want, s[i].trim);
                i++;
        }
	//        for(int i=0; i<nmis; i++)
        //{
        //        int n = missing[i];
	//                 mag[n]->fixMissing(step);
        //}
        delete [] missing;
	delete [] s;
}



int Magnet::is_pebble(int all, int pebble)
{
	return (all || family_head <0)
		&& (pebble == 0 || pebble == pebble_index);
}

int Parameter::is_pebble(int /* all */, int pebble)
{
	return (pebble == 0 || pebble == pebble_index);
}


void WireUp::writePebble(int file, int pebble, int step)
{
	StrengthItem *s = new StrengthItem[npps];

	int i=0;
	for(int n=0; n<nparms; n++)
	{
		if( allParms[n]->is_pebble(0, pebble) )
		{
			strcpy(s[i].swn, allParms[n]->swn);
			s[i].want = dllParms[n]->want_strength[step];
			s[i].trim = dllParms[n]->trim_strength[step];
			i++;
		}
	}
// 	for(int ii=0; ii<i; ii++)
// 		printf("new: %s %f %f\n",s[ii].swn, s[ii].want, s[ii].trim );
	write(file, &i, sizeof(int));
	write(file, s, i*sizeof(StrengthItem));
	delete [] s;
}

char * WireUp::magnet_names(int all, int pebble)
{
	int l=0;
	for(int n=0; n<nmag; n++)
	{
		if(all || mag[n]->family_head <0)
		if(pebble == 0 || pebble == mag[n]->pebble_index)
		{
			l += strlen(mag[n]->swn)+1;
		}
	}
	for(int n=0; n<nrfparm; n++)
	{
		if(pebble == 0 || pebble == rfparm[n]->pebble_index)
		{
			l += strlen(rfparm[n]->swn)+1;
		}
	}
	if(l == 0) return strDup("-");

	char * s = new char[l];
	l=0;
	for(int n=0; n<nmag; n++)
	{
		if(all || mag[n]->family_head <0)
		if(pebble == 0 || pebble == mag[n]->pebble_index)
		{
			strcpy(s+l,  mag[n]->swn);
			l += strlen(mag[n]->swn);
			s[l++]=',';
		}
	}
	s[l-1] = 0;
	return s;
}


char * WireUp::magnet_names(IndexArray& a)
{
	int l=0;
	for(int i=0; i<a.n; i++)
	{
		l += strlen(mag[a.a[i]]->swn)+1;
	}
	if(l == 0) return strDup("-");

	char * s = new char[l];
	l=0;
	for(int i=0; i<a.n; i++)
	{
		strcpy(s+l,  mag[a.a[i]]->swn);
		l += strlen(mag[a.a[i]]->swn);
		s[l++]=',';
	}
	s[l-1] = 0;
	return s;
}




void WireUp::get_trimArray(IndexArray& a, int step, double * values)
{
	for(int i=0; i<a.n; i++)
		values[i] = dag[a.a[i]]->trim_strength[step];
}
void WireUp::get_summArray(IndexArray& a, int step, double * values)
{
	for(int i=0; i<a.n; i++)
		values[i] = dag[a.a[i]]->summ_strength[step];
}
void WireUp::get_wantArray(IndexArray& a, int step, double * values)
{
	for(int i=0; i<a.n; i++)
		values[i] = dag[a.a[i]]->want_strength[step];
}
void WireUp::set_trimArray(IndexArray& a, int step, double * values)
{
	for(int i=0; i<a.n; i++)
		dag[a.a[i]]->set_trim_strength(step, values[i]);
}
void WireUp::set_wantArray(IndexArray& a, int step, double * values)
{
	for(int i=0; i<a.n; i++)
		dag[a.a[i]]->set_want_strength(step, values[i]);
}
void WireUp::add_trimArray(IndexArray& a, int step, double * values)
{
	for(int i=0; i<a.n; i++)
	{
		double d = values[i] + dag[a.a[i]]->trim_strength[step];
		dag[a.a[i]]->set_trim_strength(step, d);
	}
}
void WireUp::add_wantArray(IndexArray& a, int step, double * values)
{
	for(int i=0; i<a.n; i++)
	{
		double d = values[i] + dag[a.a[i]]->want_strength[step];
		dag[a.a[i]]->set_want_strength(step, d);
	}
}


void WireUp::max(int step)
{
	for(int i=0; i<nwires; i++) kabel[i]->max(this, step);
}

void WireUp::min(int step)
{
	for(int i=0; i<nwires; i++) kabel[i]->min(this, step);
}



void WireUp::mangle(int step)
{
	printf("** mangle step %d\n", step);
	for(int n=0; n<nmag; n++) dag[n]->mangle(step);
}



void MagnetD::fix(int step, int what)
{
	if(m->family_head < 0)
	if( fabs(want_current[step]+trim_current[step]-summ_current[step])  > 1.e-5)
	{
	   double diff = want_current[step]+trim_current[step]-summ_current[step];
	   printf("%-12s: fix currents %9.3f %9.3f -> %9.3f diff = %9.3f\n",
		m->swn, want_current[step],
		trim_current[step],summ_current[step],diff);
	   if(what)
	   {
		set_want_current(step, summ_current[step]);
		set_trim_current(step, 0.);
	   }
	   else
	   {
		set_trim_current(step, summ_current[step]-want_current[step]);
	   }
	}
}

void WireUp::fix(int step, int what)
{
	// if what == 1 the want current is set and trim current is zero
	// if what == 0 the want current is unchanged and the trim current is set.
	for(int n=0; n<nmag; n++)
		dag[n]->fix(step, what);
}



int WireUp::check(int step)
{
	int trunc = 0;
	for(int i=0; i<nwires; i++)
	{
// 		printf("family %d\n", i);
		trunc += kabel[i]->check(this, step);
	}
	return trunc;
}

int WireUp::check1(int step)
{
	int trunc = 0;
	for(int i=0; i<nwires; i++)
	{
// 		printf("family %d\n", i);
		trunc += kabel[i]->check1(this, step);
	}
	return trunc;
}

int WireUp::check1w(int step)
{
	int trunc = 0;
	for(int i=0; i<nwires; i++)
	{
// 		printf("family %d\n", i);
		trunc += kabel[i]->check1w(this, step);
	}
	return trunc;
}

void WireUp::check2(int step)
{
	for(int i=0; i<nwires; i++)
	{
// 		printf("family %d\n", i);
		kabel[i]->check2(this, step);
	}
}


void WireUp::check2w(int step)
{
	for(int i=0; i<nwires; i++)
	{
// 		printf("family %d\n", i);
		kabel[i]->check2w(this, step);
	}
}



void WireUp::checkp(int step)
{
	for(int i=0; i<nwires; i++)
		kabel[i]->checkp(this, step);
}



void WireUp::checkpp(int step)
{
        for(int i=0; i<nwires; i++)
                kabel[i]->checkpp(this, step);
}







//--------------------------------------------------------/

Parameter::Parameter(char * swn, WireUp * w, 
		int pebble_index, int myindex, int mytype)
{
	this->w = w;
	this->myindex = myindex;
	this->mytype = mytype;
	this->swn = strDup(swn);
	this->pebble_index = pebble_index;
}


RfParm::RfParm(char * swn, WireUp * w, int pebble_index, int myindex)
	:Parameter(swn, w, pebble_index, myindex, 1)
{
}

Magnet::Magnet(char * swn, WireUp * w, int lattice_index,
		int field_index, int foil_index, 
		int type_index, int pebble_index, int myindex)
	: Parameter(swn, w, pebble_index, myindex, 0)
{
	this->lattice_index = lattice_index;
	this->field_index = field_index;
	this->foil_index = foil_index;
	this->type_index = type_index;
	family_head= -1;
	family_polar=  1;

	nfam=0;
	fam = NULL;
	mat=NULL;
}



Magnet::~Magnet()
{
	delete [] fam;
}

MagnetD::~MagnetD()
{
	delete [] want_current;
	delete [] trim_current;
	delete [] summ_current;
	delete [] summ_currentp;
}


Parameter::~Parameter()
{
	delete [] swn;
}





const int * Magnet::famIndices(int & n)
{
	if(family_head >= 0) return w->mag[family_head]->famIndices(n);
	else
	{
		n=nfam;
		return fam;
	}
}
 
 
double Magnet::s2i(WireUp * ww, int step, double s)
{
	if(field_index < 0)
	{
//              printf("no field table for %s\n", swn);
		return(0.);
	}
	int q = w->atom_Z;
	double ko = ww->Konstante(step);
	double integ_field = s * ko/q;
	double c = w->ff[field_index]->strom(integ_field);
	if(type_index == 21)
		if(ww->gamma_e[step] > ww->jumpGamma)
			c = -c;
//      printf("integ %lf -> strom %lf\n",integ_field,c);
// 	if( !strcmp(swn,"yi7-qd10") )
// 	{
// 		printf("now s=%f konst = %f q= %d \n",s,ko,q);
// 	}
	return(c);
}

double Magnet::s2didb(WireUp * ww, int step, double s)
{
	if(field_index < 0)
	{
//              printf("no field table for %s\n", swn);
		return(0.);
	}
	double integ_field = s * ww->Konstante(step) / w->atom_Z;
	double c = w->ff[field_index]->dIdb(integ_field);
//      printf("integ %lf -> strom %lf\n",integ_field,c);
	return(c);
}


double Magnet::i2s(WireUp * ww, int step, double c)
{
	if(field_index < 0)
	{
//              printf("no field table for %s\n", swn);
		return(0.);
	}
	if(type_index == 21)
		if(ww->gamma_e[step] > ww->jumpGamma)
			c = -c;
	double intg_field = w->ff[field_index]->integ_field(c);
	double x = intg_field * w->atom_Z / ww->Konstante(step);
	return(x);
}



//--------------------------------------------------------/
ParameterD::ParameterD(Parameter * pa, WireUp * w)
{
	this->pa = pa;
	this->w = w;
	int nstep = w->nstep;
	want_strength = new double[nstep];
	trim_strength = new double[nstep];
	summ_strength = new double[nstep];


	for(int i=0; i<nstep; i++)
	{
		want_strength[i] = 0.;
		trim_strength[i] = 0.;
		summ_strength[i] = 0.;
	}

	want_spline=NULL;
	trim_spline=NULL;
}

MagnetD::MagnetD(Magnet * m, WireUp * w)
	: ParameterD( (Parameter *) m, w)
{
	this->m = m;
	this->w = w;
	int nstep = w->nstep;
	want_current = new double[nstep];
	trim_current = new double[nstep];
	summ_current = new double[nstep];
	summ_currentp = new double[nstep];


	for(int i=0; i<nstep; i++)
	{
		want_current[i] = 0.;
		trim_current[i] = 0.;
		summ_current[i] = 0.;
	}
}


RfParmD::RfParmD(RfParm * r, WireUp * w)
	: ParameterD( (Parameter *) r, w)
{
	this->r = r;
}

ParameterD::~ParameterD()
{
	delete [] want_strength;
	delete [] trim_strength;
	delete [] summ_strength;
	delete want_spline;
	delete trim_spline;
}


void ParameterD::set_want_strengthX(int step, double s)
{
	want_strength[step] = s;
}
void ParameterD::set_trim_strengthX(int step, double s)
{
	trim_strength[step] = s;
}

void ParameterD::set_want_strength(int step, double s)
{
	want_strength[step] = s;
}
void ParameterD::set_trim_strength(int step, double s)
{
	trim_strength[step] = s;
}

void MagnetD::set_want_strength(int step, double s)
{
	want_strength[step] = s;
	double c=m->s2i(w, step, s);
	set_want_current(step, c);
	set_trim_strength(step,trim_strength[step]);
}


void MagnetD::mangle(int step)
{
	if(m->family_head < 0)
	{
		double c=0.;
		for(int n=0; n< m->nfam; n++)
		{
			double cc = w->dag[m->fam[n]]->get_want_current(step);
			printf("mangle %s I = %f\n", w->mag[m->fam[n]]->swn, cc);
			c += cc;
		}
		c /= m->nfam;
		printf("mangle %s Iavg = %f\n", m->swn, c);
		set_want_current(step, c);

		c=0.;
		for(int n=0; n< m->nfam; n++)
		{
			double cc = w->dag[m->fam[n]]->get_trim_current(step);
			c += cc;
		}
		c /= m->nfam;
		set_trim_current(step, c);
	}
}


double MagnetD::get_want_current(int step)
{
	double c = m->s2i(w, step, want_strength[step]);
	c *= m->family_polar;
	return c;
}

double MagnetD::get_trim_current(int step)
{
	double c = m->s2i(w, step, trim_strength[step] + want_strength[step]) - want_current[step];
	c *= m->family_polar;
	return c;
}


void MagnetD::set_want_current(int step, double c)
{
	if(m->family_head >= 0)
	{
		w->dag[m->family_head]->set_want_current(step,
			c*m->family_polar);
	}
	else
	{
		for(int n=0; n< m->nfam; n++)
			 w->dag[m->fam[n]]->xset_want_current(step, c);
	}
}



void MagnetD::xset_want_current(int step, double c)
{
	c *= m->family_polar;
	want_current[step]=c;
	want_strength[step] = m->i2s(w, step, c);
}







void MagnetD::set_trim_strength(int step, double s)
{
	trim_strength[step] = s;
	double c=m->s2i(w, step, s + want_strength[step]) - want_current[step];
	set_trim_current(step, c);
}


void MagnetD::set_trim_current(int step, double c)
{
	if(m->family_head >= 0)
	{
		w->dag[m->family_head]->set_trim_current(step,
			c*m->family_polar);
	}
	else
	{
		for(int n=0; n< m->nfam; n++)
			 w->dag[m->fam[n]]->xset_trim_current(step, c);
	}
}



void MagnetD::xset_trim_current(int step, double c)
{
	c *= m->family_polar;
	trim_current[step]=c;
	trim_strength[step] = m->i2s(w, step, c + want_current[step])
		-want_strength[step];
}



void MagnetD::set_summ_current(int step, double c)
{
	if(m->family_head >= 0)
	{
		w->dag[m->family_head]->set_summ_current(step,
			c*m->family_polar);
	}
	else
	{
		for(int n=0; n< m->nfam; n++)
			 w->dag[m->fam[n]]->xset_summ_current(step, c);
	}
}



void MagnetD::xset_summ_current(int step, double c)
{
	c *= m->family_polar;
	summ_current[step]=c;
	summ_strength[step] = m->i2s(w, step, c);
// 	printf("%s: summ_current = %f (%f %f)\n",
// 		m->swn, summ_current[step], want_current[step], trim_current[step]);
}



void WireUp::make_fit( int step1, int step0, double /* stretch */, double p0, double p1, double * result)
{
	double dtmax[4];
	int who[4];
	for(int i=0; i<4; i++)
	{
		dtmax[i]=0.;
		who[i]= 0;
	}

	for(int n=0; n< npps; n++)
	{
		double * dt =dps[n]->make_fit(step1, step0);
		dt[3] = dt[0]+dt[1]+dt[2];

		for(int i=0; i<4; i++)
		{
			if(dt[i] > dtmax[i])
			{
// 				printf("%-15s: %9.3f %9.3f %9.3f %9.3f", pps[n]->swn,
// 					dps[n]->current[step0], dps[n]->currentp[step0],
// 					dps[n]->current[step1], dps[n]->currentp[step1]);
// 				printf(" time: %9.3f %9.3f %9.3f\n", dt[0], dt[1], dt[2]);
				dtmax[i]=dt[i];
				who[i]=n;
			}
		}
	}
	for(int i=0; i<4; i++)
	{
		printf("step %d %d time %d = %s %f\n", step1,step0, i, pps[who[i]]->swn, dtmax[i]);
	}


// 	printf("step %d %d  maxtime: %d %f %d %f %d %f\n",
// 		step1,step0, who[0], dtmax[0], who[1], dtmax[1], who[2], dtmax[2]);

	int j=who[3];
	RampFit * fit = dps[j]->fit;
	printf("%s\n", pps[j]->swn);
	fit->print_result("w::makefit best ");

	int err = fit->minimize3(dtmax[0], dtmax[1], dtmax[2]);
	printf("minimize err = %d\n", err);		
	printf("%s\n", pps[j]->swn);
	fit->print_result("w::makefit minimize3 ");
	for(int i=0; i<4; i++) result[i]=fit->dt[i];
	result[4] = fit->r;
	result[5] = fit->l;
	double g;
	double di = (fit->current[3] - fit->current[0]);
	if(fabs(di) < 1.e-5) g = 1;
	else g = (p1 - p0)/ (fit->current[3] - fit->current[0]);
	result[6] = fit->up[0] * g;
	result[9] = fit->up[2] * g;
	result[7] = fit->imax * g;
}

double * PowerD::make_fit(int step1, int step0)
{
	double * dt = fit->dt;
	fit->set_ramp(current[step0], currentp[step0], current[step1], currentp[step1]);
// 	if(current[step0] !=current[step1])
// 		printf("%-15s: %9.3f %9.3f %9.3f %9.3f", p->swn,
// 			current[step0], currentp[step0], current[step1], currentp[step1]);
	int err = fit->shortRamp(0);
	if(err) printf("shortramp err %s %d\n", p->swn, err);
// 	if(current[step0] !=current[step1])
// 		printf(" err %2d time: %9.3f %9.3f %9.3f\n", err, dt[0], dt[1], dt[2]);
	return dt;
}

void PowerD::testpp(int step, double dpdt)
{
	double diff = current[29] - current[28];
	double dIdp = diff /0.02;
	double dIdt = dIdp * dpdt;
	if(fabs(currentp[step]) > 0.001)
	{
		printf("%s didp = %f %f +: %f -: %f =%f \n", p->swn, dIdt, currentp[step],
			current[28],current[29], current[step]);
	}


}

void WireUp::set_currentp( int step, double dGammadt, double pseudo, double pseudoP)
{
	double konstq = konstante[step] / atom_Z;
	double dKonstdt = dGammadt * konstq *  gamma_e[step] / (gamma_e[step]*gamma_e[step]-1.);
	for(int n=0; n<nmag; n++)
	{
		dag[n]->set_currentp(step, konstq, dKonstdt, pseudo, pseudoP);
	}
	checkp(step);
}
void MagnetD::set_currentp(int step, double konstq, double dKonstdt, double pseudo, double pseudoP)
{
	if(m->family_head >= 0) return;

	double dkdp = want_spline->primeLookup(pseudo) + trim_spline->primeLookup(pseudo);
	double dkdt = dkdp * pseudoP;
	summ_strength[step] = want_strength[step]+ trim_strength[step];
	double dBdt = dKonstdt*summ_strength[step] + konstq*dkdt;

	double dIdB = m->s2didb(w, step,summ_strength[step]); // !!

// 	double im = m->s2i(w, step,summ_strength[step]-0.01);
// 	double i0 = m->s2i(w, step,summ_strength[step]);
// 	double ip = m->s2i(w, step,summ_strength[step]+0.01);
// 	double ddIdk = (ip-im)/0.02;
// 	double dIdk= dIdB*konstq;
// 	printf("%s: didk = %f ddidk = %f\n", m->swn, dIdk, ddIdk);

	double dIdt = dIdB * dBdt;
	summ_currentp[step] = dIdt;
	
// 	if(fabs(dIdt) > 1.e-6)
// 	{
//         	printf("%s: I = %f, dIdt= %f",m->swn, summ_current[step], dIdt);
// 		if(step > 0) printf(", c- = %f", summ_current[step-1]);
// 		printf(", c+ = %f\n", summ_current[step+1]);
// 	}
}


void MagnetD::testpp(int step, double dpdt, double ddk)
{
	if(step != 1) return;
	double diff = summ_current[29] - summ_current[28];
	double dIdp = diff /0.02;
	double dIdt = dIdp * dpdt;

	summ_strength[28] = want_strength[28]+ trim_strength[28];
	summ_strength[29] = want_strength[29]+ trim_strength[29];
	diff = summ_strength[29] - summ_strength[28];
	double dkdp = diff /0.02;
	double dkdt = dkdp * dpdt;

	if(fabs(summ_currentp[step]) > 0.001)
	{
		printf("%s mag didp = di:%f sup:%f psu=%f 28:%f 29:%f 1:%f \n",
			m->swn, dIdt, summ_currentp[step], dpdt,
			summ_current[28],summ_current[29], summ_current[step]);
		printf("%s mag dkdp = dk:%f ddk:%f psu=%f 28:%f 29:%f 1:%f \n",
			m->swn, dkdt, ddk, dpdt,
			summ_strength[28],summ_strength[29], summ_strength[step]);
	}

}







Power::Power(char * swn, int init, double max_current, double min_current,
	double max_advice, double min_advice, double i_rating, double v_rating,
	double bit_per_amp_rdb, double mvolts, double bit_per_amp,
	double max_ramp_speed, double resistance, double inductance,
	int christie)
{
	this->swn = strDup(swn);
	mat=NULL;
	this->init            = init;
	this->max_current     = max_current;
	this->min_current     = min_current;
	this->max_advice      = max_advice;
	this->min_advice      = min_advice;
	this->i_rating        = i_rating;
	this->v_rating        = v_rating;
	this->bit_per_amp_rdb = bit_per_amp_rdb;
	this->mvolts          = mvolts;
	this->bit_per_amp     = bit_per_amp;
	this->max_ramp_speed  = max_ramp_speed;
	this->resistance      = resistance;
	this->inductance      = inductance;
	this->christie        = christie;

}

Power::~Power()
{
	delete [] swn;
}






int PowerD::set_current(int step, double c)
{
// !!!!!! this code is disabled until the data base is correct!!! 
//         printf("setcurrent: %s = %f\n", p->swn, c);
	current[step]=c;
// 	if(fabs(c) < 1.e-5) return 0;	// ps is off, thats ok.
	if(c > p->max_current+1.e-5)
	{
		printf("%s: truncated from %f to %f\n",
			p->swn, c, p->max_current);
// 		if(c > 0. && p->max_current < 0.)
// 		{
// 			current[step] = 0.;
// 		}
// 		else
// 		{
			current[step]=p->max_current;
// 		}
		return 1;
	}
	if(c < p->min_current-1.e-5)
	{
		printf("%s: truncated from %f to %f\n",
			p->swn, c, p->min_current);
// 		if(c < 0. && p->min_current > 0. ) 
// 		{
// 			current[step] = 0.;
// 		}
// 		else
// 		{
			current[step]=p->min_current;
// 		}
		return 1;
	}
	return 0;
}



PowerD::PowerD(Power * p, WireUp * w)
{
	this->p = p;
	this->w = w;
	int nstep = w->nstep;
	current = new double[nstep];
	currentp = new double[nstep];

	for(int i=0; i<nstep; i++)
	{
		current[i] =  0.;
		currentp[i] =  0.;
	}

	fit = new RampFit(p->inductance, p->resistance, p->v_rating, p->max_ramp_speed,
			p->max_current, p->min_current);
	fit2 = new RampFit(p->inductance, p->resistance, p->v_rating, p->max_ramp_speed,
			p->max_current, p->min_current);
}



PowerD::~PowerD()
{
	delete [] current;
	delete [] currentp;
	delete fit;
	delete fit2;
}





/*------- Fields class methods -----------------------------------*/
 Fields::Fields(int k)
 {
	count = k;
	current =    new double[k];
	int_field=   new double[k];
	sp = NULL;
 }

 Fields::~Fields()
 {
	delete [] current;
	delete [] int_field;
	delete sp;
 }
  
void Fields::make_spline()
{
	delete sp;
	sp = new Spline(int_field, current, count,
		1.e33, 1.e33, SPLINE_ARRAY_USE);
}
 
double Fields::integ_field(double c)
{
	return sp->inverseLookup(c);
}
 
double Fields::strom(double int_f)
{
	return sp->lookup(int_f);
}
 
 
double Fields::dIdb(double int_f)
{
	return sp->primeLookup(int_f);
}
 
 
/*------- Wires class methods -----------------------------------*/
Wires::Wires(int n)
{
	Wires::n=n;
	mags = new int[n];
	ppss = new int[n];
	mag_con = NULL;
	pps_con = NULL;
}
 
Wires::~Wires()
{
	delete [] mags;
	delete [] ppss;
	for(int i=0; i<n; i++) delete mag_con[i];
	for(int i=0; i<n; i++) delete pps_con[i];
	delete [] mag_con;
	delete [] pps_con;
}


void Wires::checkpp( WireUp * w, int step)
{
        for(int i=0; i<n; i++)
        {
                double m=0.;
		int pr = !strcmp(w->mag[ mags[i] ]->swn, "bi5-qf3");
                for(int j=0; j<mag_con[i]->count; j++)
                {
                        int k=mag_con[i]->index[j];
                        int v=mag_con[i]->polarity[j];
                        double c = w->dps[k]->currentp[step];
			if(pr) printf("bi5-qf3 += %s %d %f \n", w->pps[k]->swn, v, c);
                        if(mag_con[i]->polarity[j] > 0) m += c; else m -= c; 
                }
                w->dag[ mags[i] ]->summ_currentp[step]=m;
        }
}

 
 
void Wires::checkp( WireUp * w, int step)
{
	for(int i=0; i<n; i++)
	{
		double p=0.;
// 		int pr = !strcmp(w->pps[ ppss[i] ]->swn, "y6-q7-ps");
		for(int j=0; j<pps_con[i]->count; j++)
		{
                        int k=pps_con[i]->index[j];
                        int v=pps_con[i]->polarity[j];
		     	double s =  w->dag[k]->summ_currentp[step];
// 			if(pr) printf("y6-q7-ps += %s %d %f \n", w->mag[k]->swn, v, s);
		     	if(v > 0) p += s; else p -= s;
		}
		w->dps[ ppss[i] ]->currentp[step]= p;
	}
}


 
int Wires::check( WireUp * w, int step)
{
	int trunc = check1(w, step);
	check2(w, step);
	return trunc;

}

 
int Wires::check1( WireUp * w, int step)
{
	int trunc = 0;
	for(int i=0; i<n; i++)
	{
// 		printf("pps %s = ", w->pps[ppss[i]]->swn);
		double p=0.;
		for(int j=0; j<pps_con[i]->count; j++)
		{
			int k=pps_con[i]->index[j];
			double c = w->dag[k]->want_current[step] +
				   w->dag[k]->trim_current[step];
// 			int si;
// 			if(pps_con[i]->polarity[j] == 1) si = '+';
// 			else
// 			if(pps_con[i]->polarity[j] == -1) si = '-';
// 			else si = '%';
// 			printf("  %c%s(%f) ", si, w->mag[k]->swn, w->dag[k]->want_current[step]);
			if(pps_con[i]->polarity[j] > 0) p += c; else p -= c;
		}
//         	printf(" = %f\n", p);
		trunc += w->dps[ ppss[i] ] ->set_current(step, p);
	}
	return trunc;
}

 
int Wires::check1w( WireUp * w, int step)
{
	for(int i=0; i<n; i++)
	{
		double p=0.;
		for(int j=0; j<pps_con[i]->count; j++)
		{
			int k=pps_con[i]->index[j];
			double c = w->dag[k]->want_current[step] +
				   w->dag[k]->trim_current[step];
			if(pps_con[i]->polarity[j] > 0) p += c; else p -= c;
		}
		w->dps[ ppss[i] ] -> current[step]=p;
	}
	return 0;  // no limit checking
}

 
void Wires::check2( WireUp * w, int step)
{

	for(int i=0; i<n; i++)
	{
// 		printf("mag %s = ", w->mag[mags[i]]->swn);
		double m=0.;
		for(int j=0; j<mag_con[i]->count; j++)
		{
			int k=mag_con[i]->index[j];
			double c = w->dps[k]->current[step];
// 			int si;
// 			if(mag_con[i]->polarity[j] == 1) si = '+';
// 			else
// 			if(mag_con[i]->polarity[j] == -1) si = '-';
// 			else si = '%';
// 			printf("  %c%s(%f) ", si, w->pps[k]->swn, w->dps[k]->current[step]);
			if(mag_con[i]->polarity[j] > 0) m += c; else m -= c;
		}
		int k = mags[i];
		w->dag[k]->set_summ_current(step, m);
//         	printf(" = %f\n", m);
	}

// 	printf("%5d\n",n);
// 	for(int i=0; i<n; i++)
// 	{
// 		int m = mags[i];
// 		int p = ppss[i];
// 		printf("%-12s %15f %-12s %15f %15f\n",
// 			w->mag[m]->swn,
// 			w->dag[m]->want_current[step],
// 			w->pps[p]->swn,
// 			w->dps[p]->current[step],
// 			w->dag[m]->summ_current[step]);
// 	}

}
 
 
 
void Wires::check2w( WireUp * w, int step)
{

	for(int i=0; i<n; i++)
	{
// 		printf("mag %s = ", w->mag[mags[i]]->swn);
		double m=0.;
		for(int j=0; j<mag_con[i]->count; j++)
		{
			int k=mag_con[i]->index[j];
			double c = w->dps[k]->current[step];
			if(mag_con[i]->polarity[j] > 0) m += c; else m -= c;
		}
		int k = mags[i];
		w->dag[k]->set_want_current(step, m);
	}

}
 
 
void Wires::max( WireUp * w, int step)
{
	for(int i=0; i<n; i++)
	{
		double m=0.;
		for(int j=0; j<mag_con[i]->count; j++)
		{
			int k=mag_con[i]->index[j];
			double c = w->pps[k]->max_current;
			if(mag_con[i]->polarity[j] > 0) m += c; else m -= c;
		}
		w->dag[ mags[i] ]->set_summ_current(step, m);
	}
}
 
 
 
void Wires::min( WireUp * w, int step)
{
	for(int i=0; i<n; i++)
	{
		double m=0.;
		for(int j=0; j<mag_con[i]->count; j++)
		{
			int k=mag_con[i]->index[j];
			double c = w->pps[k]->min_current;
			if(mag_con[i]->polarity[j] > 0) m += c; else m -= c;
		}
		w->dag[ mags[i] ]->set_summ_current(step, m);
	}
}
 
 
 
 
RfParmD::~RfParmD()
{
}


Cavity::Cavity(char * swn)
{
	this->swn = strDup(swn);
}
Cavity::~Cavity()
{
	delete [] swn;
}


CavityD::CavityD(Cavity * c,WireUp *w)
{
	this->c = c;
	this->w = w;
	int nstep = w->nstep;
	voltage = new double[nstep];
}

CavityD::~CavityD()
{
	delete [] voltage;
}


CWire1::CWire1(int nc1, int nc2, int ncom)
{
        this->nc1 = nc1;
	this->nc2 = nc2;
	this->ncom = ncom;
        this->cav1 = new int[nc1];
        this->cav2 = new int[nc2];
        this->common = new int[ncom];

	
}
CWire1::~CWire1()
{
	delete [] cav1;
	delete [] cav2;
	delete [] common;
}
CWire2::CWire2(int nc1)
{
        this->nc1 = nc1;
        this->cav1 = new int[nc1];
}
CWire2::~CWire2()
{
	delete [] cav1;
}


void WireUp::recalc(int nStones, int ramptype, double stretch, double * /* gamma */,
	double * pseudo, double * pseudoP, int * pseudoPflag, double result[][MAXRESULT])
{
// 	set_gamma(28, gamma[1]);
// 	set_gamma(29, gamma[1]);

	make_splines(nStones, pseudo);
// 	fill_step(28,pseudo[1]-0.01); check(28);
// 	fill_step(29,pseudo[1]+0.01); check(29);
// 	for(int n=0; n<nmag; n++)
// 	{
// 		double yl = dag[n]->want_strength[28];
// 		double ym = dag[n]->want_strength[ 1];
// 		double yh = dag[n]->want_strength[29];
// 		double d1 = fabs(yl - ym);
// 		double d2 = fabs(yh - ym);
// 		if( (d1 > 1.e-8) || (d2 > 1.e-8 )  )
// 			printf("%s fill %15f < %15f < %15f\n", mag[n]->swn, yl,ym,yh);
// 	}


	for(int step=0; step<nStones; step++)
	{
		printf("******** check step %d *******\n", step);
		check(step);
		for(int i=0; i< MAXRESULT; i++) result[step][i] = 0.;
	}
	int last = nStones-1;

	{
	make_fit(last, 0, stretch, pseudo[0], pseudo[last], result[0]);
	RampFit f(result[0], pseudo[0], pseudoP[0]);
	f.print_result("find global pseudop");
	for(int i=0; i< MAXRESULT; i++) result[1][i] = result[0][i];
	}
	int step0 =0;
	for(int step=1; step<nStones; step++)
	{                                               
		if(step == last || pseudoPflag[step] )
		{
			if( (step - step0) > 2)
			{
				make_fit(step, step0, stretch, pseudo[step0], pseudo[step], result[step]);
				RampFit f(result[step], pseudo[step0], pseudoP[step0]);
// 				f.print_result("find pseudop");
				for(int s = step0+1; s<step; s++)
				{
					double tt = f.cur2dt(pseudo[s]);
					double p,pp;
					f.dt2cur(tt, p, pp );
					printf("step %d time %f pseudo %f p %f pp %f\n", s, tt,pseudo[s],p,pp);
					pseudoP[s]=pp;
					result[s][3] = result[s][1] = tt;
					result[s][0] = result[s][2] = 0.;
					
				}
			}
			step0 = step +1;
		}

	}
	for(int step=0; step<nStones; step++)
	{                                               
		printf("%f %f ", pseudo[step], pseudoP[step]);
		for(int i=0; i< MAXRESULT; i++) printf("%f ", result[step][i]);
		printf("\n");
	}


	make_splines(nStones, pseudo);
	for(int step=0; step<nStones; step++)
	{
		printf("***recalc set_currentp: step %d\n",step);
		double dGammadt = ramptype ? pseudoP[step] : 0;
		set_currentp( step, dGammadt, pseudo[step],pseudoP[step]);
// 		if(step == 1)
// 		{
// 			fill_step(28,pseudo[step]-0.01);
// 			fill_step(29,pseudo[step]+0.01);
// 			check(28);
// 			check(29);
// 			for(int n=0; n<npps; n++) dps[n]->testpp(step, pseudoP[step]);
// 		}
	}




	wfgRowsX[0]=0.;
	wfgRowsY[0]= pseudo[0];
	int nr= 0;
	double ttot = 0.;
	double tolerance = 3.e-4;
// 	FILE * fp = fopen("_plot","w");
#define KLUGE_W
#ifdef KLUGE_W
	RampFit f(result[0], pseudo[0], pseudoP[0]);
	f.print_result("pseudo track");

        for(int j=0; j<10; j++)
        {
		wfgRowsN = 1;
		wfgRowsN = f.track(tolerance, 1000, wfgRowsX, wfgRowsY, wfgRowsN);
                printf(" tol = %f, # of rows = %d\n", tolerance, wfgRowsN);
                if(wfgRowsN <= 128 && wfgRowsN > 90) break;
                if(wfgRowsN > 1000) 
                {
                        printf("max rows > 1000\n");
                        return;
                }
                double f= double(110)/double(wfgRowsN);
                tolerance /= f*f;
        }
        if(wfgRowsN > 128) 
        {
                printf("max rows > 128\n");
                return;
        }




	result[0][8] = tolerance;
	for(int i=nr; i<wfgRowsN; i++)
	{
		printf("wfg i%d %f %f\n", i, wfgRowsX[i], wfgRowsY[i]);
		ttot += wfgRowsX[i];
// 		fprintf(fp,"%f %f\n",  ttot, wfgRowsY[i]);
	}
	nr = wfgRowsN;
#else
	for(int step=1; step<nStones; step++)
	{                                               
		printf("***recalc make_fit: step %d\n",step);
		make_fit(step, step-1, stretch, pseudo[step-1], pseudo[step], result[step]);
	}

	for(int step=1; step<nStones; step++)
	{                                               
		printf("***recalc plot: step %d\n",step);
		RampFit f(result[step], pseudo[step-1], pseudoP[step-1]);

        	f.set_ramp(pseudo[step-1], pseudoP[step-1], pseudo[step], pseudoP[step]);
		int err = f.minimize(f.dt[3]);
		if(err) printf("pseudo track err\n");


		f.print_result("pseudo track");
		wfgRowsN = f.track(1.e-3, 5000, wfgRowsX, wfgRowsY, wfgRowsN);
		printf("\nstep=%d, result = ",step);
		for(int i=0; i<MAXRESULT; i++) printf(" %f ", result[step-1][i]);
		printf("\n");
		for(int i=nr; i<wfgRowsN; i++)
		{
			printf("wfg i%d %f %f\n", i, wfgRowsX[i], wfgRowsY[i]);
			ttot += wfgRowsX[i];
// 			fprintf(fp,"%f %f\n",  ttot, wfgRowsY[i]);
		}
		nr = wfgRowsN;
	}
#endif
// 	fclose(fp);
}




void WireUp::recalc1( int nStones, int ramptype, double stretch, double * gamma,
	double * pseudo, double * pseudoP, int * pseudoPflag, double result[][MAXRESULT])
{
	for(int step=0; step<nStones; step++) check(step);
	make_splines(nStones, pseudo);

	int last = nStones-1;
	pseudoP[0] = pseudoP[last]= 0.;
	pseudoPflag[0] = pseudoPflag[last]= 0;
	for(int step=1; step<last; step++)
	{                                               
		if(!pseudoPflag[step])
		{
			for(int n=0; n<npps; n++) dps[n]->recalc1(step);
			checkpp(step);

			double konstq = konstante[step] / atom_Z;
			double gaga = ramptype ? gamma[step]/(gamma[step]*gamma[step]-1.) : 0.;
			double pp = 1.e20;
			double pp1;
			int np = -1;
			for(int n=0; n<nmag; n++)
			{
				int f = dag[n]->recalc1(step, pseudo[step], gaga, konstq, pp1);
				if(f)
				{
					if(fabs(pp1) < pp)
					{
						pp = pp1;
						np = n;
					}
				}
			}
			pseudoP[step] = pp;
			printf("\n*** step %d best pp %f from %s\n\n", step, pp, mag[np]->swn);
		}
		else
		{
			pseudoP[step] = 0;
		}
	}

	for(int step=0; step<nStones; step++)
	{
		printf("***recalc set_currentp: step %d\n",step);
		double dGammadt = ramptype ? pseudoP[step] : 0;
		set_currentp( step, dGammadt, pseudo[step],pseudoP[step]);
	}




	wfgRowsN = 1;
	wfgRowsX[0]=0.;
	wfgRowsY[0]= pseudo[0];
	int nr= 0;
	double ttot = 0.;
	FILE * fp = fopen("_plot","w");
	for(int step=1; step<nStones; step++)
	{                                               
		printf("***recalc make_fit: step %d\n",step);
		make_fit(step, step-1, stretch, pseudo[step-1], pseudo[step], result[step]);
		printf("***recalc plot: step %d\n",step);
		RampFit f(result[step], pseudo[step-1], pseudoP[step-1]);

        	f.set_ramp(pseudo[step-1], pseudoP[step-1], pseudo[step], pseudoP[step]);
		int err = f.minimize(f.dt[3]);
		if(err) printf("pseudo track err\n");


		f.print_result("pseudo track");
		wfgRowsN = f.track(1.e-3, 5000, wfgRowsX, wfgRowsY, wfgRowsN);
		printf("\nstep=%d, result = ",step);
		for(int i=0; i<MAXRESULT; i++) printf(" %f ", result[step-1][i]);
		printf("\n");
		for(int i=nr; i<wfgRowsN; i++)
		{
			printf("wfg i%d %f %f\n", i, wfgRowsX[i], wfgRowsY[i]);
			ttot += wfgRowsX[i];
			fprintf(fp,"%f %f\n",  ttot, wfgRowsY[i]);
		}
		nr = wfgRowsN;
	}
	fclose(fp);
}





int MagnetD::recalc1(int step, double p, double gaga, double konstq, double & dPdt)
{
	if(m->family_head >= 0) return 0;
	double dkdp = want_spline->primeLookup(p) + trim_spline->primeLookup(p);
	double y = dkdp + gaga;
	if(fabs(y) < 1.e-10) return 0;
	summ_strength[step] = want_strength[step]+ trim_strength[step];
	double dIdB = m->s2didb(w, step,summ_strength[step]);
	dPdt = summ_currentp[step]/(dIdB * konstq * y);
	printf("Step %d move %s dp %f\n", step, m->swn, dPdt); 
	return 1;
}

void PowerD::recalc1(int step)
{
	int stepp = step + 1;
	int stepm = step - 1;
	fit->set_ramp(current[stepm], 0., current[step], 0.);
	int err = fit->shortRamp(0);
	double t1 = fit->dt[3];
	if(err) printf("shortramp err %s %d\n", p->swn, err);
	currentp[step] = fit->current_p[2];
	fit2->set_ramp(current[step], 0., current[stepp], 0.);
	err = fit2->shortRamp(0);
	double t2 = fit2->dt[3];
	if(err) printf("shortramp err %s %d\n", p->swn, err);
	currentp[step] += fit2->current_p[1];
	currentp[step] /= 2.;
	double tt = t1+t2;
	if(tt > 0.01)
	{
		double sl = (current[stepp] - current[stepm])/tt;
		printf("Step %3d slope %-15s dp %12.6f sl %12.6f tt %12.6f\n",
			step, p->swn, currentp[step], sl, tt); 
// 		fit->print_result("first");
// 		fit2->print_result("second");
	}

}
