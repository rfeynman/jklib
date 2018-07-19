#include <SDDS/SDDS.h>




// Class to hold the 'stonelist' info
// This will come from the DB later
class StepStoneList 
{
public:
        StepStoneList();
        virtual ~StepStoneList();
/// return the row# if found or Nrows if the stone name is not found
        int findIndex(const char*);
        int Nrows;
	int MaxRows;
	void write();
	void freelist();
	int addStone(const char* stone, const char * bg, const char * yg,
		const char * gg, const char * bc, const char * yc, const char * rf);


        char** Stones;
        char** BlueGlobal;
        char** YellowGlobal;
        char** GreenGlobal;
        char** BlueCorrector;
        char** YellowCorrector;
        char** RF;
protected:
        SDDS_TABLE Input;
};
 
 
class Pebble
{
public:
        Pebble(const char* type, const char *name);
        Pebble();
        virtual ~Pebble();
        long Nrows;
        char **Swns;
        double *Want, *Trim;
protected:
        SDDS_TABLE Input, Output;
};
 

