#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>


#include "rampFit.hxx"
#include "rampcontainer.hxx"

int RampContainer::newPebbleCounter = 0;


RampContainer:: RampContainer(const char *rampfile, int what=0)
{

// what = 0 : old  format
// what = 1 : sdds format

	w = new WireUp();
	w->wInit(MAXSTONE+1);
	for(int i=0; i<MAXSTONE; i++)
	for(int l=0; l<MAXPEBBLE; l++)
	{
		oldPebbleHistory[i][l]=NULL;
		pebbleHistory[i][l]=NULL;
	}
	this->rampfile=strDup(rampfile);
	for(int i=0; i<MAXPEBBLE; i++)
		w->get_indexArray(!i, i, index_array[i]);
	ok = 0;


	if(!what)
	{
        	char *appstore = getenv("APP_STORE");
        	if (!appstore)
		{
			printf("APP_STORE is not defined\n");
               		exit(-1);
        	}
        	char s[1024];

        	strcpy(s, appstore);
        	strcat(s, "/ramps/oldform/");
		int err = chdir(s);
		if(err) { ok = err; return;}

		FILE * rf = fopen(rampfile, "r");
		if(!rf) return;
		char s1[100], s2[100];

		fscanf(rf, "%s%s", s1, s2);
		version = strDup(s2);
		printf("version: %s\n", version); 

		fscanf(rf, "%s%s", s1, s2);
		nickname = strDup(s2);
		printf("nickname: %s\n", nickname); 

		fscanf(rf, "%s%d", s1, &ramptype);
		printf("ramptype: %d\n", ramptype); 

		fscanf(rf, "%s%d", s1, &histLen);
		old_history = new char[histLen+1];
        	while( getc(rf) != '\n');
        	for(int i=0; i<histLen; i++) old_history[i]=getc(rf);
        	old_history[histLen]=0;
		histLen = 14 + strlen(rampfile);
		history = new char[histLen+1];
		strcpy(history, "Created from ");
		strcat(history, rampfile);
		strcat(history, "\n");
	
		fscanf(rf, "%s%lf", s1, &stretch);
		fscanf(rf, "%s%lf", s1, &injectionVoltageBlue);
		fscanf(rf, "%s%lf", s1, &injectionFrequencyBlue);
		fscanf(rf, "%s%lf", s1, &injectionRadiusBlue);
		fscanf(rf, "%s%lf", s1, &jumpGammaBlue);


		fscanf(rf, "%s%d", s1, &nStones);
// 	printf("stones: %s %d\n", s1,nStones); 
		w->clear_all();
		for(int i=0; i<nStones; i++)
		{
			char s[60];
			for(int r=0; r<MAXRESULT; r++)
                		if(fscanf(rf,"%lf",result[i]+r) != 1)
				{
					ok=1;
					return;
				}
			if( fscanf(rf,"%lf%lf%lf%d%s",
				gamma+i,pseudo+i,pseudoP+i,pseudoPflag+i,s) != 5)
			{
				ok=1;
				return;
			}
			if(gamma[i] <= 1.)
			{
				printf("gamma lessequal 1, idiot!!\n");
				ok=1;
				return;
			}
// 		printf("stone %s\n",s);
			w->set_gamma(i, gamma[i]);
			readStone(i, s);
		}
		fclose(rf);

	}
	else
	{


        	char *appstore = getenv("APP_STORE");
        	if (!appstore)
		{
			printf("APP_STORE is not defined\n");
                	exit(-1);
        	}
        	char s[1024];

        	strcpy(s, appstore);
        	strcat(s, "/ramps/data/");
        	strcat(s, rampfile);

		SDDS_TABLE Input;
		int err = SDDS_InitializeInput(&Input, s);
		if(!err)
		{
			printf("cant open file %s\n", s);
			ok=1;
			return;
		}

		if( SDDS_ReadTable(&Input) <= 0 )
		{
			SDDS_PrintErrors(stderr, 0);
			SDDS_Terminate(&Input);
			ok=1;
			return;
		}

        	SDDS_SetColumnFlags(&Input, 1);
        	SDDS_SetRowFlags(&Input, 1);
        	nStones = SDDS_CountRowsOfInterest(&Input);
        	if (nStones == -1)
		{
                	SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
			ok=1;
                	return;
        	}
        	// We need the type of the ramp to be there
        	if (!SDDS_GetParameter(&Input, "Type", (void *) &ramptype))
		{
                	printf("No type for ramp\n");
                	SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
			ok=1;
                	return;
        	}
        	w->clear_all();

		// pull the stretch out, if available
		stretch = 1.;
		injectionVoltageBlue= 1.;
		injectionFrequencyBlue = 1.;
		injectionRadiusBlue = 1.;
		injectionVoltageYellow= 1.;
		injectionFrequencyYellow = 1.;
		injectionRadiusYellow = 1.;
		jumpGammaBlue = 0.;
		jumpGammaYellow = 0.;
		rfRampTime = 10;

        	SDDS_GetParameter(&Input, "RfRampTime", (void *) &rfRampTime);
        	SDDS_GetParameter(&Input, "Stretch", (void *) &stretch);
        	SDDS_GetParameter(&Input, "InjVoltsBlue", (void *) &injectionVoltageBlue);
        	SDDS_GetParameter(&Input, "InjFreqBlue", (void *) &injectionFrequencyBlue);
        	SDDS_GetParameter(&Input, "InjRadiusBlue", (void *) &injectionRadiusBlue);
        	SDDS_GetParameter(&Input, "InjVoltsYel", (void *) &injectionVoltageYellow);
        	SDDS_GetParameter(&Input, "InjFreqYel", (void *) &injectionFrequencyYellow);
        	SDDS_GetParameter(&Input, "InjRadiusYel", (void *) &injectionRadiusYellow);
        	SDDS_GetParameter(&Input, "JumpGammaBlue", (void *) &jumpGammaBlue);
        	SDDS_GetParameter(&Input, "JumpGammaYellow", (void *) &jumpGammaYellow);


        	double * PseudoP = SDDS_GetColumnInDoubles(&Input, "PseudoP");
        	if (PseudoP)
		{
			for(int i=0; i< nStones; i++) pseudoP[i]=PseudoP[i];
        	}
		else
		{
			for(int i=0; i< nStones; i++) pseudoP[i]= 0.;
        	}
		free(PseudoP);

        	if (ramptype == 0)
		{ // A correction/squeeze ramp
                	printf("Correction Ramp\n");
			double Gamma;
                	if (!SDDS_GetParameter(&Input, "Gamma", (void *) &Gamma)) {
                		printf("No gamma for ramp\n");
                        	return;
                	}
			for(int i=0; i< nStones; i++)
			{
				gamma[i]= Gamma;
				w->set_gamma(i, gamma[i]); 
			}
                	double * Pseudo = SDDS_GetColumnInDoubles(&Input, "Pseudo");
                	if (!Pseudo) {
                        	printf("Not a proper Acceleration Ramp file\n");
				ok=1;
                        	return;
                	}
			for(int i=0; i< nStones; i++) pseudo[i]=Pseudo[i];
			free(Pseudo);
        	}
		else
		{ // An acceleration ramp
                	double * Pseudo = SDDS_GetColumnInDoubles(&Input, "Gamma");
 
                	if (!Pseudo)
			{
                        	printf("Not a proper Acceleration Ramp file\n");
				ok=1;
                        	return;
                	}
			for(int i=0; i< nStones; i++)
			{
				pseudo[i] = gamma[i]= Pseudo[i];
				w->set_gamma(i, gamma[i]); 
			}
        	}
        	double * T1 = SDDS_GetColumnInDoubles(&Input, "T1");
        	double * T2 = SDDS_GetColumnInDoubles(&Input, "T2");
        	double * T3 = SDDS_GetColumnInDoubles(&Input, "T3");
        	double * Ts = SDDS_GetColumnInDoubles(&Input, "Tsum");
       		double * Re = SDDS_GetColumnInDoubles(&Input, "Resistance");
        	double * In = SDDS_GetColumnInDoubles(&Input, "Inductance");
        	double * Up = SDDS_GetColumnInDoubles(&Input, "Uprime");
        	double * Up2= SDDS_GetColumnInDoubles(&Input, "Uprime2");
        	double * Im = SDDS_GetColumnInDoubles(&Input, "Imax");
        	double * To = SDDS_GetColumnInDoubles(&Input, "Tolerance");
        	long * Fl = (long *) SDDS_GetColumn(&Input, "Flag");
		for(int i=0; i< nStones; i++)
		{
			if(T1)   result[i][0] = T1[i];
			if(T2)   result[i][1] = T2[i];
			if(T3)   result[i][2] = T3[i];
			if(Ts)   result[i][3] = Ts[i];
			if(Re)   result[i][4] = Re[i];
			if(In)   result[i][5] = In[i];
			if(Up)   result[i][6] = Up[i];
			if(Up2)  result[i][9] = Up2[i];
			else     result[i][9] = -Up[i];
			if(Im)   result[i][7] = Im[i];
			if(To)   result[i][8] = To[i];
			if(Fl)   pseudoPflag[i] = Fl[i];
			printf("resu %d ", i);
			for(int j=0; j<MAXRESULT; j++) printf(" %f ", result[i][j]);
			printf("\n");
		}
		free(T1); free(T2); free(T3); free(Ts); free(Re);
		free(In); free(Up); free(Im); free(To); free(Fl);




        	old_history = strdup("johannes has not told me yet");

	
        	char ** Stones = (char**) SDDS_GetColumn(&Input, "StepStone");
		char** BlueGlobal = (char**) SDDS_GetColumn(&Input, "BlueGlobal");
		char** YellowGlobal = (char**) SDDS_GetColumn(&Input, "YellowGlobal");
		char** GreenGlobal = (char**) SDDS_GetColumn(&Input, "GreenGlobal");
		char** BlueCorrector = (char**) SDDS_GetColumn(&Input, "BlueCorrector");
		char** YellowCorrector = (char**) SDDS_GetColumn(&Input, "YellowCorrector");
		char** RF = (char**) SDDS_GetColumn(&Input, "RF");

        	if (!Stones || !BlueGlobal || !YellowGlobal || !GreenGlobal || !BlueCorrector || !YellowCorrector || !RF)
		{
                	printf("Not a proper Acceleration Ramp file\n");
			ok=1;
                	return;
        	}
		for(int i=0; i< nStones; i++)
		{
			strcpy(stoneName[i][0], Stones[i]);
			strcpy(stoneName[i][1], BlueGlobal[i]);
			strcpy(stoneName[i][2], YellowGlobal[i]);
			strcpy(stoneName[i][3], GreenGlobal[i]);
			strcpy(stoneName[i][4], BlueCorrector[i]);
			strcpy(stoneName[i][5], YellowCorrector[i]);
			strcpy(stoneName[i][6], RF[i]);

			printf("%s ip,pp,g %f %f %f\n",
				stoneName[i][0], pseudo[i], pseudoP[i], gamma[i]);
	
			Pebble bg("BlueGlobal",BlueGlobal[i]);
			sddsMagPebble(i,bg, 1);
			Pebble bc("BlueCorrector",BlueCorrector[i]);
			sddsMagPebble(i,bc, 3);
			Pebble gg("GreenGlobal",GreenGlobal[i]);
			sddsMagPebble(i,gg, 5);
			Pebble yg("YellowGlobal",YellowGlobal[i]);
			sddsMagPebble(i,yg, 2);
			Pebble yc("YellowCorrector",YellowCorrector[i]);
			sddsMagPebble(i,yc, 4);
			Pebble rf("RF",RF[i]);
			sddsRFPebble(i,rf, 6);

			for(int pebble=0; pebble<MAXPEBBLE; pebble++)
			printf("%s ", stoneName[i][pebble]);
			printf("\n");
		}
		SDDS_FreeStringArray(Stones, nStones);

	}
	for(int pebble=0; pebble<MAXPEBBLE; pebble++)
	{
		for(int i=0; i<nStones; i++)
		{
			stoneNameArray[pebble][i]=stoneName[i][pebble];
		}
		stoneNameArray[pebble][nStones] = NULL;
	}
}

void RampContainer::sddsRFPebble(int step, Pebble & p, int pebble)
{
 
}
 
void RampContainer::sddsMagPebble(int step, Pebble & p, int pebble)
{
        for(int n=0; n<p.Nrows; n++)
        {
		int element_index = w->name2Index(p.Swns[n]);
		if(element_index >= 0)
		{
			w->dllParms[element_index]->set_want_strength(step,p.Want[n]);
			w->dllParms[element_index]->set_trim_strength(step,p.Trim[n]);
		}
        }
}


int RampContainer::set_stepStone(int step)
{
	if(step < 0 || step >= nStones) return -1;
	return step;
}
int RampContainer::set_pebble(int pebble)
{
	if(pebble < 1 || pebble > 5) return -1;
	return pebble;
}





void RampContainer::set_modified(int step, int pebble)
{
	if( modified[step][pebble]) return;
	newPebbleCounter++;
	char pn[40];
	sprintf(pn,"new%04d",newPebbleCounter);
	strcpy(stoneName[step][pebble], pn);
	modified[step][pebble] = 1;
	set_modified(step, 0);
}



void RampContainer::readStone(int step, char * newName)
{
// 		printf("readstone %d %s\n", step, newName);
		FILE * db = fopen("stonelist","r");
		if(!db)
		{
			printf("cant open stonelist\n");
			ok = 2;
			return;
		}

		char s1[MAXPEBBLE][40];
		while( fscanf(db,"%s%s%s%s%s%s%s",
			s1[0],s1[1],s1[2],s1[3],s1[4],s1[5],s1[6]) == MAXPEBBLE)
		{
// 			printf("stone %s stonelist %s -> %s %s %s %s %s %s\n",
// 				newName, s1[0],s1[1],s1[2],s1[3],s1[4],s1[5],s1[6]);
			if( !strcmp(newName,s1[0]) )
			{
				for(int j=1; j<MAXPEBBLE; j++)
				{
					strcpy(stoneName[step][j], s1[j]);
					modified[step][j] = 0;
				}
				ok=0;
// 				printf("found stone %s -> %s %s %s %s %s %s\n",
// 					s1[0],s1[1],s1[2],s1[3],s1[4],s1[5],s1[6]);
				goto found;
			}
		}
		printf("stone %s not found in stonelist\n", newName);
		ok = 3;
		fclose(db);
		return;

found:		fclose(db);
		strcpy(stoneName[step][0], newName);
		modified[step][0] = 0;
//	 	w->clear_all(step);
		for(int pebble=1; pebble<MAXPEBBLE; pebble++)
		{
			readPebble(step, pebble, stoneName[step][pebble]);
// 			if(stoneName[step][pebble][0] == 'A') break;
		}
}


RampContainer::RampContainer (RampContainer * r , double t)
{
	w = new WireUp();
	w->wInit(MAXSTONE+1);
	printf("%p %f %p\n", r,t,w);
}

void RampContainer::sddsSave(const char * file)
{
       	char *appstore = getenv("APP_STORE");
       	if (!appstore)
	{
		printf("APP_STORE is not defined\n");
               	exit(-1);
       	}
       	char s[1024];

       	strcpy(s, appstore);
       	strcat(s, "/ramps/data/");
       	strcat(s, file);

	SDDS_TABLE Output;
        int err = SDDS_InitializeOutput(&Output, SDDS_BINARY, 1, NULL, "Ramp version 1.0", s);
	if(!err)
	{
		printf("cant make sdds output ramp file: %s\n",s);
		return;
	}
        SDDS_DefineParameter(&Output, "Type",            NULL, NULL, NULL, NULL, SDDS_LONG,   NULL);
        SDDS_DefineParameter(&Output, "Stretch",         NULL, NULL, NULL, NULL, SDDS_DOUBLE, NULL);
        SDDS_DefineParameter(&Output, "InjVoltsBlue",    NULL, NULL, NULL, NULL, SDDS_DOUBLE, NULL);
        SDDS_DefineParameter(&Output, "InjFreqBlue",     NULL, NULL, NULL, NULL, SDDS_DOUBLE, NULL);
        SDDS_DefineParameter(&Output, "InjRadiusBlue",   NULL, NULL, NULL, NULL, SDDS_DOUBLE, NULL);
        SDDS_DefineParameter(&Output, "JumpGammaBlue",   NULL, NULL, NULL, NULL, SDDS_DOUBLE, NULL);
        SDDS_DefineParameter(&Output, "InjVoltsYel",     NULL, NULL, NULL, NULL, SDDS_DOUBLE, NULL);
        SDDS_DefineParameter(&Output, "InjFreqYel",      NULL, NULL, NULL, NULL, SDDS_DOUBLE, NULL);
        SDDS_DefineParameter(&Output, "InjRadiusYel",    NULL, NULL, NULL, NULL, SDDS_DOUBLE, NULL);
        SDDS_DefineParameter(&Output, "JumpGammaYellow", NULL, NULL, NULL, NULL, SDDS_DOUBLE, NULL);
        if (ramptype == 0)
	{
                SDDS_DefineColumn(&Output, "Pseudo", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0);
        	SDDS_DefineParameter(&Output, "Gamma", NULL, NULL, NULL, NULL, SDDS_DOUBLE, NULL);
        }
	else
	{
                SDDS_DefineColumn(&Output, "Gamma", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0);
        }

//      double *T1, *T2, *T3, *Tsum, *Resistance, *Inductance, *Uprime, *Imax, *Tolerance; long *Flag
        SDDS_DefineColumn(&Output, "T1",         NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0);       
        SDDS_DefineColumn(&Output, "T2",         NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0);       
        SDDS_DefineColumn(&Output, "T3",         NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0);       
        SDDS_DefineColumn(&Output, "Tsum",       NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0);     
        SDDS_DefineColumn(&Output, "Resistance", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0);       
        SDDS_DefineColumn(&Output, "Inductance", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0);       
        SDDS_DefineColumn(&Output, "Uprime",     NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0);   
        SDDS_DefineColumn(&Output, "Imax",       NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0);     
        SDDS_DefineColumn(&Output, "Tolerance",  NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0);        
        SDDS_DefineColumn(&Output, "PseudoP",    NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0);  
        SDDS_DefineColumn(&Output, "Flag",       NULL, NULL, NULL, NULL, SDDS_LONG,   0);

        SDDS_DefineColumn(&Output, "StepStone",  NULL, NULL, NULL, NULL, SDDS_STRING, 0);
        SDDS_DefineColumn(&Output, "BlueGlobal",  NULL, NULL, NULL, NULL, SDDS_STRING, 0);
        SDDS_DefineColumn(&Output, "YellowGlobal",  NULL, NULL, NULL, NULL, SDDS_STRING, 0);
        SDDS_DefineColumn(&Output, "GreenGlobal",  NULL, NULL, NULL, NULL, SDDS_STRING, 0);
        SDDS_DefineColumn(&Output, "BlueCorrector",  NULL, NULL, NULL, NULL, SDDS_STRING, 0);
        SDDS_DefineColumn(&Output, "YellowCorrector",  NULL, NULL, NULL, NULL, SDDS_STRING, 0);
        SDDS_DefineColumn(&Output, "RF",  NULL, NULL, NULL, NULL, SDDS_STRING, 0);

        SDDS_WriteLayout(&Output);
        SDDS_StartTable(&Output, nStones);

	SDDS_SetParameters(&Output, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, 
		"Type",      ramptype,
		"Stretch",   stretch,
		"RfRampTime",   rfRampTime,
		"InjVoltsBlue",  injectionVoltageBlue,
		"InjFreqBlue",   injectionFrequencyBlue,
		"InjRadiusBlue", injectionRadiusBlue,
		"JumpGammaBlue", jumpGammaBlue,
		"InjVoltsYel",  injectionVoltageYellow,
		"InjFreqYel",   injectionFrequencyYellow,
		"InjRadiusYel", injectionRadiusYellow,
		"JumpGammaYellow", jumpGammaYellow,
		NULL);

        for (int i = 0; i < nStones; ++i)
	{
		if(stoneName[i][0][0] == 'n')
		{
			char newName[60];
			sprintf(newName,"step%05d",time(NULL));
			strcpy(stoneName[i][0], newName);
			// make new stone and save it
			for(int p=1; p<MAXPEBBLE; p++)
			{
				if(stoneName[i][p][0] == 'n')
				{
					char oldName[60];
					strcpy(oldName, stoneName[i][p]);
					for(int ii=0; ii<nStones; ii++)
						if(!strcmp(stoneName[ii][p], oldName) )
						{
							strcpy(stoneName[ii][p], newName);
							modified[ii][p] = 0;
						}
					sddsWritePebble(i,p,stoneName[i][p]);
				}
			}
		}
	}

        for (int i = 0; i < nStones; ++i)
	{
                if (ramptype == 0)
		{
			SDDS_SetParameters(&Output, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, 
				"Gamma",      gamma[0], NULL);
                        SDDS_SetRowValues(&Output, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, i,
				"Pseudo",     pseudo[i],
				"T1",         result[i][0],
				"T2",         result[i][1],
				"T3",         result[i][2], 
                                "Tsum",       result[i][3],
				"Resistance", result[i][4],
				"Inductance", result[i][5],
				"Uprime",     result[i][6],
				"Uprime2",    result[i][9],
				"Imax",       result[i][7],
				"Tolerance",  result[i][8],
				"PseudoP",    pseudoP[i],
				"Flag",       pseudoPflag[i],
				"StepStone",  stoneName[i][0],
				NULL);
                }
		else
		{
                        SDDS_SetRowValues(&Output, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, i,
				"Gamma",            pseudo[i],
				"T1",               result[i][0],
				"T2",               result[i][1],
				"T3",               result[i][2], 
                                "Tsum",             result[i][3],
				"Resistance",       result[i][4],
				"Inductance",       result[i][5],
				"Uprime",           result[i][6],
				"Uprime2",          result[i][9],
				"Imax",             result[i][7],
				"Tolerance",        result[i][8],
				"PseudoP",          pseudoP[i],
				"Flag",             pseudoPflag[i],
				"StepStone",        stoneName[i][0],
				"BlueGlobal",       stoneName[i][0],
				"YellowGlobal",     stoneName[i][0],
				"GreenGlobal",      stoneName[i][0],
				"BlueCorrector",    stoneName[i][0],
				"YellowCorrector",  stoneName[i][0],
				"RF",               stoneName[i][0],
				NULL);
				printf("sdds stonename %s\n",stoneName[i][0]);
                }
        }
        if (!SDDS_WriteTable(&Output))
	{
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        }
        SDDS_Terminate(&Output);

}

int RampContainer::sddsWritePebble(int step, int pebble, const char *fname)
{
	char *pe[6] = {"BlueGlobal", "BlueCorrector", "GreenGlobal", "YellowGlobal", "YellowCorrector", "RF"};
        char *appstore = getenv("APP_STORE");
        if (!appstore)
	{
		printf("APP_STORE is not defined\n");
               	exit(-1);
        }
        char s[1024];
	sprintf(s,"%s/ramps/data/%s/%s",appstore,pe[pebble-1],fname);


	SDDS_TABLE Output;
        if(!SDDS_InitializeOutput(&Output, SDDS_BINARY, 1, NULL, "Pebble version 1.0", s))
	{
		printf("cant open pebble %s for write\n",s);
		return -1;
	}

        SDDS_DefineColumn(&Output, "swn", NULL, NULL, NULL, NULL, SDDS_STRING, 0);
        SDDS_DefineColumn(&Output, "want", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0);
        SDDS_DefineColumn(&Output, "trim", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0);

        SDDS_WriteLayout(&Output);
	int Nrows = index_array[pebble].n;
        SDDS_StartTable(&Output, Nrows);
        for (int i = 0; i < Nrows; ++i)
	{
		int index = index_array[pebble].a[i];
                SDDS_SetRowValues(&Output, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, i,
			"swn",  w->allParms[index]->swn,
			"want", w->dllParms[index]->want_strength[step],
			"trim", w->dllParms[index]->trim_strength[step],
			NULL);
        }
        if (!SDDS_WriteTable(&Output)) {
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        }
        SDDS_Terminate(&Output);
        return 0;
}



int RampContainer::save(const char * rampfile)
{
        char *appstore = getenv("APP_STORE");
        if (!appstore)
	{
		printf("APP_STORE is not defined\n");
               	exit(-1);
        }
        char s[1024];

        strcpy(s, appstore);
        strcat(s, "/ramps/oldform/");
	int err = chdir(s);
	if(err) return err;

        FILE * rf = fopen(rampfile, "w");
        if(!rf) return errno;
	fprintf(rf,"version: %s\n", version);
	fprintf(rf,"nickname: %s\n", nickname);
	fprintf(rf,"ramptype: %d\n", ramptype);
	fprintf(rf,"history: %d\n", histLen);
        fprintf(rf,"%s\n", history);
	fprintf(rf,"stretch: %f\n", stretch);
	fprintf(rf,"injectionVoltage: %f\n", injectionVoltageBlue);
	fprintf(rf,"injectionFrequency: %f\n", injectionFrequencyBlue);
	fprintf(rf,"injectionRadius: %f\n", injectionRadiusBlue);
	fprintf(rf,"jumpGamma: %f\n", jumpGammaBlue);

	fprintf(rf,"stones: %d\n", nStones);
        for(int step=0; step<nStones; step++)
        {
		int err = writeStepStone(step);
		if(err) return err;
        }
        for(int step=0; step<nStones; step++)
        {

		for(int r=0; r<MAXRESULT; r++)
                	fprintf(rf,"%f ",result[step][r]);
                fprintf(rf,"%f %f %f %d %s\n",
			gamma[step],
			pseudo[step],
			pseudoP[step],
			pseudoPflag[step],
			stoneName[step][0]);
        }



	fclose(rf);
	return 0;
}

int RampContainer::writeStepStone(int i)
{
	if(stoneName[i][0][0] == 'n')
	{
		// get next name
                FILE * db = fopen("stonelist","r");
                if(!db)
                {
			printf("cant open stonelist\n");
                        return -1;
                }
 
                char s1[MAXPEBBLE][40];
                while( fscanf(db,"%s%s%s%s%s%s%s",
                        s1[0],s1[1],s1[2],s1[3],s1[4],s1[5],s1[6]) == MAXPEBBLE)
			; // keep reading;
                fclose(db);
		int newstone = atoi(s1[0]+5) +1;
		char tmp[20];
		sprintf(tmp,"stepS%05d",newstone);
		strcpy(stoneName[i][0], tmp);
                for(int j=1; j<MAXPEBBLE; j++)
                {
                        int err = writePebble(i, j, newstone);
			if(err) return err;
                }
		db = fopen("stonelist","a");
		fprintf(db,"%s %s %s %s %s %s %s\n",
			stoneName[i][0],
			stoneName[i][1],
			stoneName[i][2],
			stoneName[i][3],
			stoneName[i][4],
			stoneName[i][5],
			stoneName[i][6]);
		fclose(db);
	}
	return 0;
}


int RampContainer::writePebble(int i, int j, int newstone)
{
	char * pebName[MAXPEBBLE] =
		{"stepS","bGlob","yGlob","bCorr","yCorr","gGlob","rfGbl"};
	if(stoneName[i][j][0] == 'n')
	{
		char oldName[40], newName[40];
		strcpy(oldName, stoneName[i][j]);
		sprintf(newName,"%s%05d",pebName[j],newstone);
		for(int ii=0; ii<nStones; ii++)
			if(!strcmp(stoneName[ii][j], oldName) )
			{
				strcpy(stoneName[ii][j], newName);
				modified[ii][j] = 0;
			}
                int db = open(stoneName[i][j], O_WRONLY | O_CREAT, 0444);
                if(db == -1)
                {
			printf("cant open file %s, error %d\n", stoneName[i][j], errno);
                        return -1;
                }
        	char date_str[30];
        	time_t *time_ptr = new time_t;
        	time(time_ptr);
        	strftime(date_str, 30, "%y%m%d_%H%M%S", localtime(time_ptr));
		strcat(date_str, " p 1.01");
        	delete time_ptr;
		write(db,date_str,30);

		if(pebbleHistory[i][j]) 
		{
			int l = strlen(pebbleHistory[i][j])+1;
			write(db,&l,sizeof(int));
			write(db,pebbleHistory[i][j],l);
		}
		else
		{
			int l=0;
			write(db,&l,sizeof(int));
		}

		w->writePebble(db, j, i);
		close(db);
	}
	return 0;
}

void RampContainer::readPebble(int step, int pebble, const char * file)
{
// 	printf("read pebble %d %d %s\n", step, pebble, file);
	if( file[0] == 'A')
	{
		if( ! strcmp(file, "Anull") ) return;
		FILE * fp = fopen(file,"r");
		if(!fp)
        	{
			printf("cant open file %s, error %d\n", file, errno);
                	return;
        	}
		char swn[20];
		double want, trim;
		while( fscanf(fp, "%s %lf %lf", swn,&want,&trim) == 3)
		{
			int element_index = w->name2Index(swn);
			if(element_index >= 0)
			{
				w->dllParms[element_index]->set_want_strengthX(step,want);
				w->dllParms[element_index]->set_trim_strengthX(step,trim);
				int pebble = w->allParms[element_index]->pebble_index;
			}
		}
		fclose(fp);
		w->mangle(step);
		return;
	}

	int rf = open(file, O_RDONLY);
        if(rf == -1)
        {
		printf("cant open file %s, error %d\n", file, errno);
                return;
        }

	char swn[40];
	read(rf,swn,30);
// 	printf("file %s open: %s\n", file, swn);
	int l;
	read(rf,&l, sizeof(int));
	delete [] oldPebbleHistory[step][pebble];
	delete [] pebbleHistory[step][pebble];

	if(l)
	{
		oldPebbleHistory[step][pebble] = new char[l];
		read(rf,oldPebbleHistory[step][pebble], l);
	}
	else 
	{
		oldPebbleHistory[step][pebble] = strDup("empty\n");
	}
	l  = 14 + strlen(file);
	pHistLen[step][pebble] = l;
	char * h = new char[l+1];
	strcpy(h, "Created from ");
	strcat(h, file);
	strcat(h, "\n");
	pebbleHistory[step][pebble] = h;


	w->readPebble(rf, pebble, step);
	close(rf);
}

int RampContainer::check()
{
	int trunc;
	for(int step=0; step<nStones; step++)
		trunc += w->check(step);
	return trunc;
}

int RampContainer::check(int step)
{
	return w->check(step);
}

const char * RampContainer::get_nickname()
{
	return nickname;
}

void RampContainer::set_nickname(const char * name)
{
	delete [] nickname;
	nickname = strDup(name);
}

 

RampContainer::~RampContainer()
{
	delete w;
// 	for(int i=0; i<nStones; i++)
// 	{
// 		printf("stones:%2d:  ",i);
// 		for(int l=0; l<MAXPEBBLE; l++)
// 			printf("   %s  ",stoneName[i][l]);
// 		printf("\n");
// 	}
	delete [] rampfile;
	delete [] version;
	delete [] nickname;
	delete [] history;
	delete [] old_history;
	for(int i=0; i<nStones; i++)
	for(int j=0; j<MAXPEBBLE; j++)
	{
		delete [] pebbleHistory[i][j];
		delete [] oldPebbleHistory[i][j];
	}

}


void RampContainer::get_times(double * t)
{
	for(int i=0; i<nStones; i++)
	{
		t[i] = result[i][3];
	}
}
 
void RampContainer::get_pseudos(double * p)
{
	for(int i=0; i<nStones; i++)
	{
		p[i] = pseudo[i];
	}
}
void RampContainer::get_pseudosP(double * p)
{
	for(int i=0; i<nStones; i++)
	{
		p[i] = pseudoP[i];
	}
}
void RampContainer::get_gammas(double * g)
{
	for(int i=0; i<nStones; i++)
	{
		g[i] = gamma[i];
	}
}

double RampContainer::get_time(int step)
{
		return result[step][3];
}
double RampContainer::get_pseudo(int step)
{
		return pseudo[step];
}
double RampContainer::get_pseudoP(int step)
{
		return pseudoP[step];
}
double RampContainer::get_gamma(int step)
{
		return gamma[step];
}

const char ** RampContainer::get_stoneNames(int pebble)
{
	return (const char **) stoneNameArray[pebble];
}

const char * RampContainer::get_stoneName(int step, int pebble)
{
	return stoneName[step][pebble];
}


const char ** RampContainer::get_siteWideNames()
{
	return (const char **) index_array[0].names;
}


char ** RampContainer::get_swnArray(int pebble)
{
	return index_array[pebble].names;
}


const char * RampContainer::get_oldRampHistory()
{
	return old_history;
}
const char * RampContainer::get_rampHistory()
{
	return history;
}
const char * RampContainer::get_oldPebbleHistory(int step, int pebble)
{
	return (const char * ) oldPebbleHistory[step][pebble];
}
const char * RampContainer::get_pebbleHistory(int step, int pebble)
{
	return (const char * ) pebbleHistory[step][pebble];
}
int RampContainer::get_numMags(int pebble)
{
	return index_array[pebble].n;
}
void RampContainer::get_trimArray(int step, int pebble, double * values)
{
	w->get_trimArray(index_array[pebble], step, values);
}
void RampContainer::get_wantArray(int step, int pebble, double * values)
{
	w->get_wantArray(index_array[pebble], step, values);
}
void RampContainer::get_maxStrength(int step, int pebble, double * values)
{
	w->max(step);
	w->get_summArray(index_array[pebble], step, values);
}
void RampContainer::get_minStrength(int step, int pebble, double * values)
{
	w->min(step);
	w->get_summArray(index_array[pebble], step, values);
}


void RampContainer::set_pebbleHistory(int step, int pebble, char * h)
{
	int l= strlen(h);
	if(l > pHistLen[step][pebble])
	{
		delete [] pebbleHistory[step][pebble];
		pebbleHistory[step][pebble] = strDup(h);
	}
	else
	{
		strcpy(pebbleHistory[step][pebble], h);
	}
	pHistLen[step][pebble] = l;
}
void RampContainer::add_pebbleHistory(int step, int pebble, char * h)
{
	int  l= strlen(h);
	char * r = new char[histLen+l+2];
	strcpy(r, pebbleHistory[step][pebble]);
	delete [] pebbleHistory[step][pebble];
	strcat(r, h);
	histLen += l;
	if(r[histLen] != '\n') r[histLen++]= '\n';
	pebbleHistory[step][pebble]=r;
}
void RampContainer::set_rampHistory(const char * h)
{
	int l= strlen(h);
	if(l > histLen)
	{
		delete [] history;
		history = strDup(h);
	}
	else
	{
		strcpy(history, h);
	}
	histLen = l;
}
void RampContainer::add_rampHistory(const char * h)
{
	int  l= strlen(h);
	char * t = history;
	history = new char[histLen+l+2];
	strcpy(history,t);
	delete [] t;
	strcat(history, h);
	histLen += l;
	if(history[histLen-1] != '\n')
		history[histLen++]= '\n';
}
void RampContainer::set_trimArray(int step, int pebble, double * a)
{
	w->set_trimArray(index_array[pebble], step, a);
	set_modified(step, pebble);
}
void RampContainer::set_wantArray(int step, int pebble, double * a)
{
	w->set_wantArray(index_array[pebble], step, a);
	set_modified(step, pebble);
}
void RampContainer::add_trimArray(int step, int pebble, double * a)
{
	w->add_trimArray(index_array[pebble], step, a);
	set_modified(step, pebble);
}
void RampContainer::add_wantArray(int step, int pebble, double * a)
{
	w->add_wantArray(index_array[pebble], step, a);
	set_modified(step, pebble);
}
double RampContainer::get_trimValue(const char * name, int step)
{
	int element_index = w->name2Index(name);
	if(element_index >= 0) return 0.;
	return w->dllParms[element_index]->trim_strength[step];
}
double RampContainer::get_wantValue(const char * name, int step)
{
	int element_index = w->name2Index(name);
	if(element_index >= 0) return 0.;
	return w->dllParms[element_index]->want_strength[step];
}
int RampContainer::set_trimValue(const char * name, int step, double s)
{
	int element_index = w->name2Index(name);
	if(element_index < 0) return -1;
	w->dllParms[element_index]->set_trim_strength(step,s);
	int pebble = w->allParms[element_index]->pebble_index;
	set_modified(step, pebble);
	return pebble;
}
int RampContainer::set_wantValue(const char * name, int step, double s)
{
	int element_index = w->name2Index(name);
	if(element_index < 0) return -1;
	w->dllParms[element_index]->set_want_strength(step,s);
	int pebble = w->allParms[element_index]->pebble_index;
	set_modified(step, pebble);
	return pebble;
}
int RampContainer::add_trimValue(const char * name, int step, double s)
{
	int element_index = w->name2Index(name);
	if(element_index < 0) return -1;
	s += w->dag[element_index]->trim_strength[step];
	w->dllParms[element_index]->set_trim_strength(step,s);
	int pebble = w->allParms[element_index]->pebble_index;
	set_modified(step, pebble);
	return pebble;
}
int RampContainer::add_wantValue(const char * name, int step, double s)
{
	int element_index = w->name2Index(name);
	if(element_index < 0) return -1;
	s += w->dag[element_index]->want_strength[step];
	w->dllParms[element_index]->set_want_strength(step,s);
	int pebble = w->allParms[element_index]->pebble_index;
	set_modified(step, pebble);
	return pebble;
}



void RampContainer::set_ramptype(int type)
{
	if(ramptype == type) return;
	ramptype = type;

	if(ramptype) // it is an acceleration ramp
	{
		for(int n=0; n<nStones; n++)
		{
			pseudo[n] = gamma[n];
			pseudoP[n] = 1.e33;
		}
		pseudoP[0] = pseudoP[nStones-1] = 0.;
	}
	else
	{
		for(int n=0; n<nStones; n++)
		{
			gamma[n] = gamma[0];
		}
	}
	recalc();
}






int RampContainer::set_atStone(const char * p)
{
	for(int i=0; i<nStones; i++)
	{
		if( !strcmp(stoneName[i][0], p)) return i;
	}
	return -1;
}
void RampContainer::replace_stone(int step, int pebble, char * newName)
{
	for(int n=0; n<nStones; n++)
	if( !strcmp(stoneName[n][pebble], newName) )
	{
		if(pebble == 0)
		{
		   for(int i=0; i<MAXPEBBLE; i++)
		   {
			strcpy(stoneName[step][i], stoneName[n][i]);
			modified[step][i] =0;
			modified[n][i] =0;
			w->copy_pebble(step,n,i);
		   }
		}
		else
		{
			strcpy(stoneName[step][pebble], stoneName[n][pebble]);
			modified[step][pebble] =0;
			modified[n][pebble] =0;
			w->copy_pebble(step,n,pebble);
		}	
		recalc();
		return;
	}
	if(pebble == 0) readStone(step, newName);
	else            readPebble(step, pebble, newName);
	recalc();
}


int RampContainer::delete_stone(int step)
{
	if(nStones <= 2) return -1;

	w->delete_step(step, nStones);
	nStones--;
	for(int j=step; j< nStones; j++)
	{
		int m=j+1;      
		pseudo[j] =  pseudo[m];
		pseudoP[j] =  pseudoP[m];
		gamma[j] =  gamma[m];
		for(int k=0; k<MAXRESULT; k++) result[j][k]=result[m][k];
		for(int pebble=0; pebble<MAXPEBBLE; pebble++)
		{
			strcpy(stoneName[j][pebble], stoneName[m][pebble]);
			modified[j][pebble] = modified[m][pebble];
		}
	}
	for(int pebble=0; pebble<MAXPEBBLE; pebble++)
	{
		stoneName[nStones][pebble][0]=0;
		stoneNameArray[pebble][nStones] = NULL;
	}
	if(step >= nStones) step = nStones-1;
	return step;
}


int RampContainer::set_atPseudo(double p)
{
	int step;
	for(step=0; step<nStones; step++)
	{
		if(fabs(p-pseudo[step]) < 1.e-3) return step;
		if(p < pseudo[step]) break;
	}

	w->set_gamma(MAXSTONE, ramptype ? gamma[0] : p);
	w->fill_step(MAXSTONE,p);	
	w->insert_step(step, nStones);
	w->copy_step(step, MAXSTONE);

	insert_step(step);
	pseudo[step] = p;
	pseudoP[step] = 1.e33;
	for(int k=0; k<MAXRESULT; k++) result[step][k]=0.;
	for(int pebble=1; pebble<MAXPEBBLE; pebble++)
		set_modified(step, pebble);
	nStones ++;
	for(int pebble=0; pebble<MAXPEBBLE; pebble++)
	{
		stoneName[nStones][pebble][0]=0;
		stoneNameArray[pebble][nStones] = NULL;
		stoneNameArray[pebble][nStones-1] = stoneName[nStones-1][pebble];
	}
	w->make_splines(nStones, pseudo);
	return step;
}


void RampContainer::insert_step(int i)
{
	for(int j= nStones;j>i; j--)
	{
		int m=j-1;      
		pseudo[j] =  pseudo[m];
		pseudoP[j] =  pseudoP[m];
		gamma[j] =  gamma[m];
		for(int k=0; k<MAXRESULT; k++) result[j][k]=result[m][k];
		for(int pebble=0; pebble<MAXPEBBLE; pebble++)
		{
			strcpy(stoneName[j][pebble], stoneName[m][pebble]);
			modified[j][pebble] = modified[m][pebble];
		}
	}
	for(int pebble=0; pebble<MAXPEBBLE; pebble++)
	{
		stoneName[i][pebble][0] = 0;
		modified[i][pebble] = 0;
	}
}




void RampContainer::set_gamma(int step, double value)
{
	gamma[step] = value;
	if(ramptype)
	{
		pseudo[step]=value;
		w->set_gamma(step,value);
	}
	else
	for(int st=0; st < nStones; st++)
	{
		pseudo[st]=value;
		w->set_gamma(st,value);
	}
		
}
void RampContainer::set_pseudo(int step, double value)
{
	pseudo[step] = value;
}
void RampContainer::set_pseudoP(int step, double value)
{
	pseudoP[step] = value;
}
void RampContainer::set_stretch(double t)
{
	stretch = t;
	recalc();
}
int RampContainer::set_nearTime(double t)
{
	printf("neartime: %f\n", t);
	return 0; // returns step
}
int RampContainer::set_nearPseudo(double p)
{
	double diff = fabs(p-pseudo[0]);
	int step = 0;
	for(int i=1; i<nStones; i++)
	{
		double diff1 = fabs(p-pseudo[i]);
		if(diff1 < diff)
		{
			diff = diff1;
			step = i;
		}
	}
	return step;
}
 


int RampContainer::set_atTime(double time)
{
	printf("attime: %f\n", time);
	return 0; // returns step
}




void RampContainer::recalc()
{
	w->recalc(nStones, ramptype, stretch, gamma, pseudo, pseudoP, pseudoPflag, result);
}







