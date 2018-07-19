#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

const int maxKnobs = 20;
const int maxKnobListeners = 100;

class KnobContainer
{
public:
	KnobContainer();
	~KnobContainer();
 

        int save(const char * filename);
        int load(const char * filename);

        int addKnob(const char * name);
        int findKnob(const char * name);
        void deleteKnob(int kIndex);
        void deleteAllKnobs();
        void clearKnob(int kIndex);
 
        void set_knobMin(int kIndex, double min);
        double get_knobMin(int kIndex);
 
        void set_knobMax(int kIndex, double max);
        double get_knobMax(int kIndex);
 
        void set_knobUnit(int kIndex, const char * unit);
        const char * get_knobUnit(int kIndex);
 
        char ** knobNames();
        double * knobMins();
        double * knobMaxs();
 
        int knobAddMagnet(int kIndex, char * magnet);
        int knobFindMagnet(int kIndex, char * magnet);
        void knobDeleteMagnet(int mIndex);
 
        void set_knobConst(int mIndex, double cons);
        double get_knobConst(int kIndex, int mIndex);
 
        int nKnobMags(int kIndex);
        char ** knobMagNames(int kIndex);
        double * knobMagConst(int kIndex);
 
// private stuff
// -------------
// private:


 
	int nKnobs;
        char label[maxKnobs][40];
        double min[maxKnobs];
        double max[maxKnobs];
	char unit[maxKnobs][40];
 
        int nKnobListener;
        char swn[maxKnobListeners][20];
        int knobIndex[maxKnobListeners];
        double constant[maxKnobListeners];
};
