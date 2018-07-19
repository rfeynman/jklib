#include "wireup.hxx" 
#include "stonelist.hxx"
#include <fcntl.h>
#include <errno.h>

const int MAXSTONE = 32;
const int MAXPEBBLE = 7; // 0 = stone name, 1-6 = pebble name


class RampContainer
{
public:

	RampContainer (const char * rampfile, int what=0);
	RampContainer (RampContainer *, double);
	virtual ~RampContainer ();
	int ok;
	static int newPebbleCounter;

	int save(const char * filename);
	int check(); // check all
	int check(int step);
	void recalc();


// Ramp type
// ----------

	int get_ramptype() { return ramptype; }
	void set_ramptype(int type);  // type 1 = acceleration, 0=other

// History and comments
// --------------------

	const char * get_nickname(); // do not free
	void set_nickname(const char * name);

	const char * get_oldRampHistory(); // do not free
	const char * get_rampHistory(); // do not free
	void set_rampHistory(const char * name);
	void add_rampHistory(const char * name);

	const char ** get_siteWideNames(); // do not free,
					   // last pointer is NULL

// stone and pebble number
// -----------------------
	
	int set_stepStone(int step); // returns step or -1
	int set_pebble(int pebble);  // returns pebble or -1

// Steping stones
// --------------

	int get_numStones() { return nStones; }

	void get_times(double * t); // must supply t[nStones]

	void get_pseudos(double * p); // must supply p[nStones]

	void get_pseudosP(double * p); // must supply p[nStones]

	void get_gammas(double * g); // must supply g[nStones]

	const char ** get_stoneNames(int pebble); // dont delete


// Select a step stone:
// ----------------------


	double get_time(int step);
	int set_atTime(double); // returns step
	int set_nearTime(double); // returns step

	double get_pseudo(int step);
	int set_atPseudo(double); // returns step
	int set_nearPseudo(double); // returns step

	int set_atStone(const char *stonename); // returns step

// Modify stepstones:
// ------------------

	int delete_stone(int step); // returns step

	double get_stretch() { return stretch; }
	void set_stretch(double);

	void set_pseudo(int step, double value);

	double get_pseudoP(int step);
	void set_pseudoP(int step, double value);

	double get_gamma(int step);
	void set_gamma(int step, double value);

	const char * get_stoneName(int step, int pebble);  //dont delete
	void replace_stone(int step, int pebble, char * newName);

// Pebble operations:
// ----------------------

	const char * get_oldPebbleHistory(int step, int pebble);
	const char * get_pebbleHistory(int step, int pebble);
	void set_pebbleHistory(int step, int pebble, char *);
	void add_pebbleHistory(int step, int pebble, char *);

	int get_numMags(int pebble);

	char ** get_swnArray(int pebble);  // do not free,
					   // last pointer is NULL


					// values is existing array of
					// length 'get_numMags(pebble)'
	void get_trimArray(int step, int pebble, double * values);
	void set_trimArray(int step, int pebble, double * values);
	void add_trimArray(int step, int pebble, double * values);

	void get_wantArray(int step, int pebble, double * values);
	void set_wantArray(int step, int pebble, double * values);
	void add_wantArray(int step, int pebble, double * values);

	void get_maxStrength(int step, int pebble, double * values);
	void get_minStrength(int step, int pebble, double * values);

	double get_trimValue(const char * name, int step);
	int set_trimValue(const char * name, int step, double s); // returns pebble
	int add_trimValue(const char * name, int step, double s); // returns pebble

	double get_wantValue(const char * name, int step);
	int set_wantValue(const char * name, int step, double s); // returns pebble
	int add_wantValue(const char * name, int step, double s); // returns pebble



// private stuff
// -------------
// private:
	WireUp * w;
	char * rampfile;
	char * version;
	char * nickname;
	int ramptype;
	char * history;
	char * old_history;
	int histLen;



	int nStones;
	int nMags;
	IndexArray index_array[MAXPEBBLE];
        double stretch;
        double injectionVoltageBlue;
        double injectionFrequencyBlue;
        double injectionRadiusBlue;
        double injectionVoltageYellow;
        double injectionFrequencyYellow;
        double injectionRadiusYellow;
        double jumpGammaBlue;
        double jumpGammaYellow;
        double rfRampTime;


	double pseudo[MAXSTONE];
	double pseudoP[MAXSTONE];
	int pseudoPflag[MAXSTONE];
	double gamma[MAXSTONE];
	double result[MAXSTONE][MAXRESULT];
	char stoneName[MAXSTONE][MAXPEBBLE][50];
	char * stoneNameArray[MAXPEBBLE][MAXSTONE+1];
	char * oldPebbleHistory[MAXSTONE][MAXPEBBLE];
	char * pebbleHistory[MAXSTONE][MAXPEBBLE];
	int pHistLen[MAXSTONE][MAXPEBBLE];
	int modified[MAXSTONE][MAXPEBBLE];



	void readPebble(int step, int pebble, const char * file);
	void readStone(int step, char * newName);
	void sddsRFPebble(int step, Pebble & p, int pebble);
	void sddsMagPebble(int step, Pebble & p, int pebble);
	void sddsSave(const char * file);
	int sddsWritePebble(int step, int pebble, const char *fname);

	int writePebble(int i, int j, int newstone);
	int writeStepStone(int i);
	void insert_step(int i);
	void set_modified(int step, int pebble);

};

