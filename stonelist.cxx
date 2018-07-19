#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "jkLib/stonelist.hxx"


// Load the stonelist file
StepStoneList::StepStoneList()
{
        char *appstore = getenv("APP_STORE");
        if (!appstore) {
                printf("Need to define $APP_STORE\n");
                exit(-1);
        }
        char s[1024];
        strcpy(s, appstore);
        strcat(s, "/ramps/data/stonelist");
        int err = SDDS_InitializeInput(&Input, (char*) s);
	if(!err)
	{
		MaxRows = Nrows = 0;
		return;
	}
        long ncols;
        char **colnames = SDDS_GetColumnNames(&Input, &ncols);
        if (!colnames)
	{
		printf("stonelist %s hosed\n", s);
		MaxRows = Nrows = 0;
        	SDDS_Terminate(&Input);
                return;
        }
	for(int i=0; i<ncols; i++)
	{
		printf("stonelist cols %s\n", colnames[i]);
	}
	if(ncols != 7)
	{
		printf("stonelist %s hosed\n", s);
		MaxRows = Nrows = 0;
        	SDDS_Terminate(&Input);
                return;
        }
	SDDS_FreeStringArray(colnames, ncols);
 
        SDDS_ReadTable(&Input);
        SDDS_SetColumnFlags(&Input, 1);
        SDDS_SetRowFlags(&Input, 1);
        MaxRows = Nrows = SDDS_CountRowsOfInterest(&Input);
        Stones          = (char**) SDDS_GetColumn(&Input, "StepStone");
        BlueGlobal      = (char**) SDDS_GetColumn(&Input, "BlueGlobal");
        YellowGlobal    = (char**) SDDS_GetColumn(&Input, "YellowGlobal");
        GreenGlobal     = (char**) SDDS_GetColumn(&Input, "GreenGlobal");
        BlueCorrector   = (char**) SDDS_GetColumn(&Input, "BlueCorrector");
        YellowCorrector = (char**) SDDS_GetColumn(&Input, "YellowCorrector");
        RF              = (char**) SDDS_GetColumn(&Input, "RF");
        SDDS_Terminate(&Input);
}
 
StepStoneList::~StepStoneList()
{
	freelist();
}
 
void StepStoneList::write()
{
        char *appstore = getenv("APP_STORE");
        if (!appstore) {
                printf("Need to define $APP_STORE\n");
                exit(-1);
        }
        char s[1024];
        strcpy(s, appstore);
        strcat(s, "/ramps/data/stonelist");
        SDDS_TABLE Output;
        SDDS_InitializeOutput(&Output, SDDS_ASCII, 1, NULL, "StepStoneList version 1.0", s);
        SDDS_DefineColumn(&Output, "StepStone",       NULL, NULL, NULL, NULL, SDDS_STRING, 0);        
        SDDS_DefineColumn(&Output, "BlueGlobal",      NULL, NULL, NULL, NULL, SDDS_STRING, 0);       
        SDDS_DefineColumn(&Output, "YellowGlobal",    NULL, NULL, NULL, NULL, SDDS_STRING, 0);     
        SDDS_DefineColumn(&Output, "BlueCorrector",   NULL, NULL, NULL, NULL, SDDS_STRING, 0);    
        SDDS_DefineColumn(&Output, "YellowCorrector", NULL, NULL, NULL, NULL, SDDS_STRING, 0);  
        SDDS_DefineColumn(&Output, "GreenGlobal",     NULL, NULL, NULL, NULL, SDDS_STRING, 0);      
        SDDS_DefineColumn(&Output, "RF",              NULL, NULL, NULL, NULL, SDDS_STRING, 0);
        SDDS_WriteLayout(&Output);
        SDDS_StartTable(&Output, Nrows);

        for (int i = 0; i<Nrows; i++) 
	{
                SDDS_SetRowValues(&Output, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, i,
			 "StepStone",       Stones[i], NULL);
                SDDS_SetRowValues(&Output, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, i,
			 "BlueGlobal",      BlueGlobal[i], NULL); 
		SDDS_SetRowValues(&Output, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, i,
			 "YellowGlobal",    YellowGlobal[i], NULL); 
                SDDS_SetRowValues(&Output, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, i,
			 "BlueCorrector",   GreenGlobal[i], NULL); 
                SDDS_SetRowValues(&Output, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, i,
			 "YellowCorrector", BlueCorrector[i], NULL); 
                SDDS_SetRowValues(&Output, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, i,
			 "GreenGlobal",     YellowCorrector[i], NULL); 
                SDDS_SetRowValues(&Output, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, i,
			 "RF",              RF[i], NULL); 
        }
        SDDS_WriteTable(&Output);
        SDDS_Terminate(&Output);        

}
 
void StepStoneList::freelist()
{
	if(Nrows == 0) return;
	SDDS_FreeStringArray(Stones,          Nrows);
	SDDS_FreeStringArray(BlueGlobal,      Nrows);
	SDDS_FreeStringArray(YellowGlobal,    Nrows);
	SDDS_FreeStringArray(GreenGlobal,     Nrows);
	SDDS_FreeStringArray(BlueCorrector,   Nrows);
	SDDS_FreeStringArray(YellowCorrector, Nrows);
	SDDS_FreeStringArray(RF,              Nrows);
}
 
int StepStoneList::addStone(const char* stone, const char * bg, const char * yg,
	const char * gg, const char * bc, const char * yc, const char * rf)
{
	if(Nrows == MaxRows)
	{
		MaxRows += 20;
		char ** sto = new char*[MaxRows];
		char ** bgl = new char*[MaxRows];
		char ** ygl = new char*[MaxRows];
		char ** ggl = new char*[MaxRows];
		char ** bcl = new char*[MaxRows];
		char ** ycl = new char*[MaxRows];
		char ** rfl = new char*[MaxRows];
		for(int i=0; i<Nrows; i++)
		{
			sto[i] = Stones[i];
			bgl[i] = BlueGlobal[i];
			ygl[i] = YellowGlobal[i];
			ggl[i] = GreenGlobal[i];
			bcl[i] = BlueCorrector[i];
			ycl[i] = YellowCorrector[i];
			rfl[i] = RF[i];
		}
		freelist();
		Stones = sto;
		BlueGlobal = bgl;
		YellowGlobal = ygl;
		GreenGlobal = ggl;
		BlueCorrector = bcl;
		YellowCorrector = ycl;
		RF = rfl;
	}

	Stones[Nrows] = strdup(stone);
	BlueGlobal[Nrows] = strdup(bg);
	YellowGlobal[Nrows] = strdup(yg);
	GreenGlobal[Nrows] = strdup(gg);
	BlueCorrector[Nrows] = strdup(bc);
	YellowCorrector[Nrows] = strdup(yc);
	RF[Nrows] = strdup(rf);
	Nrows++;
	return Nrows;
}
 
int StepStoneList::findIndex(const char* stone)
{
        for (int i = 0; i < Nrows; ++i) {
                if (strcmp(stone, Stones[i]) == 0) return i;
        }
        return -1;
}
 
Pebble::Pebble()
{
        Swns = NULL;
        Want = NULL;
        Trim = NULL;
}
 
Pebble::~Pebble()
{
        SDDS_FreeStringArray(Swns, Nrows);
        free(Want);
        free(Trim);
}
 
Pebble::Pebble(const char* type, const char* name)
{
        char *appstore = getenv("APP_STORE");
        if (!appstore) {
                printf("Need to define $APP_STORE\n");
                exit(-1);
        }
        char s[1024];
        sprintf(s, "%s/ramps/data/%s/%s", appstore, type, name);
 
        SDDS_InitializeInput(&Input, s);
        printf("%s is open\n", s);
 
 
        SDDS_ReadTable(&Input);
        SDDS_SetColumnFlags(&Input, 1);
        SDDS_SetRowFlags(&Input, 1);
 
        Want = SDDS_GetColumnInDoubles(&Input, "want");
        Trim = SDDS_GetColumnInDoubles(&Input, "trim");
        Swns = (char**) SDDS_GetColumn(&Input, "swn");
 
        if (!Swns || !Want || !Trim) {
                printf("peblle file %s hosed \n",s );
        }
        Nrows = SDDS_CountRowsOfInterest(&Input);
        SDDS_Terminate(&Input);
}
 
 

