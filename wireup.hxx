#include <math.h>
#include "spline_class.hxx"
#include "rampFit.hxx"

const int MAXRESULT = 10;

double dminv(double *, int);
double rdd(FILE * f);
int rdi(FILE * f);
char * strDup(const char *);

struct StrengthItem
{
	char swn[16];
	double want;
	double trim;
};

class IndexArray
{
public:
	int n;
	int * a;
	int * t;
	char ** names;
	IndexArray() { a=NULL; t=NULL; names=NULL;}
	~IndexArray() { delete [] a; delete [] t; delete [] names;}
};

struct Speclist
{
        char species[15];
        double mass_per_u;
        int atom_A;
        int atom_Z;
        int q_ags;
};


class WireUp
{
public:
	static class Magnet  **mag;
	static int nmag;
	class MagnetD **dag;
	static class Power   **pps;
	static int npps;
	class PowerD ** dps;
	static class Wires **kabel;
	static int nwires;
	static class Fields  **ff;
	static int nff;
	static char current_species[15];
	static int q_ags;
	static int atom_A;
	static int atom_Z;
	static double mass_per_u;
	static int nspecies;
	static Speclist * specs;

	static const double clight;

	static class RfParm ** rfparm;
	static int nrfparm;
	class RfParmD ** dfparm;
	static class Cavity ** rfcavity;
	static int nrfcavity;
	class CavityD ** dfcavity;
	static class CWire1 ** rfcwire1;
	static int nrfcwire1;
	static class CWire2 ** rfcwire2;
	static int nrfcwire2;
	static class CRadius** rfcradius;
	static int nrfcradius;
	static class CFreq** rfcfreq;
	static int nrfcfreq;

 
	static class Parameter ** allParms;
	static int nparms;
	class ParameterD ** dllParms;

	int nstep;
	double * gamma_e;
	double * gamma_e_p;
	double * konstante;
	double * konstante_p;
	double konstmax;
	double jumpGamma;
        float wfgRowsX[5000];
	double wfgRowsY[5000];
	int wfgRowsN;


 
	void set_gamma(int step, double g);
	void set_jumpGamma(double g);
	void set_konstante(int step);
	double Konstante(int step);
	double Konstmax();
	void max(int step);
	void min(int step);
	int check(int step);
	int check1(int step);
	int check1w(int step);
	void check2(int step);
	void check2w(int step);
	void checkp(int step);
	void checkpp(int step);

	void make_splines(int n, double * x);
	void insert_step(int i, int n);
	void delete_step(int i, int n);
	void copy_step(int t, int f);
	void copy_pebble(int t, int f, int pebble);
	void fill_step(int i, double p);
	void clear_all();
	void clear_all(int step);
	void fix(int step, int what);
	void mangle(int step);
	int name2Index(const char * name);
	int rfName2Index(const char * name);
	int psName2Index(const char * name);

	int set_psCurrent(int step, const char * swn, double c); // return 0 = ok
	int set_wantNtrim(int step, const char * swn, double want, double trim); // return 0 = ok
	int set_want(int step, const char * swn, double want); // return 0 = ok
	int set_trim(int step, const char * swn, double trim); // return 0 = ok
	void get_indexArray(int all, int pebble, IndexArray& a);
	void get_trimArray(IndexArray& a, int step, double * values);
	void get_wantArray(IndexArray& a, int step, double * values);
	void get_summArray(IndexArray& a, int step, double * values);
	void set_trimArray(IndexArray& a, int step, double * values);
	void set_wantArray(IndexArray& a, int step, double * values);
	void add_trimArray(IndexArray& a, int step, double * values);
	void add_wantArray(IndexArray& a, int step, double * values);
	void set_currentp( int step, double dGammadt, double pseudo, double pseudoP);
	void make_fit( int step1, int step0, double stretch, double p0, double p1, double * result);
	void recalc1(int nStones, int ramptype, double stretch, double * gamma,
		double * pseudo, double * pseudoP, int * pseudoPflag, double result[][MAXRESULT]);
	void recalc(int nStones, int ramptype, double stretch, double * gamma,
		double * pseudo, double * pseudoP, int * pseudoPflag, double result[][MAXRESULT]);

	char * magnet_names(IndexArray& a); // caller must delete array
	char * magnet_names(int all, int pebble); // caller must delete array
	void writePebble(int file, int pebble, int step);
	void  readPebble(int file, int pebble, int step);

 
	WireUp(const char * file);
	WireUp();
	int wInit(int nstep);
	virtual ~WireUp();
};
 
struct  Connection
{
	int count;
	int * index;
	int * polarity;
	Connection(int n)
	{
		count = n;
		index = new int[n];
		polarity = new int[n];
	}
	virtual ~Connection()
	{
		delete [] index;
		delete [] polarity;
	}
};

class Wires
{
	Wires() {};
public:
	Wires(int n);
	virtual ~Wires();
	int n;
	int * ppss;
	int * mags;
	Connection ** mag_con;
	Connection ** pps_con;
	void max(WireUp * w, int step);
	void min(WireUp * w, int step);
	int  check(WireUp * w, int step);
	int  check1(WireUp * w, int step);
	int  check1w(WireUp * w, int step);
	void check2(WireUp * w, int step);
	void check2w(WireUp * w, int step);
	void checkp(WireUp * w, int step);
	void checkpp( WireUp * w, int step);
};

 
class Fields
{
public:
	Fields(int n);
	virtual ~Fields();
	int count;
	double *current;
	double *int_field;
	Spline *sp;
	double integ_field(double c);
	double strom(double int_f);
	double dIdb(double int_f);
	void make_spline();
};





class Parameter
{
public:
	Parameter(char * swn, WireUp * w, int pebble, int myindex, int mytype);
	virtual ~Parameter();
	class WireUp * w;
	char * swn;
	int pebble_index;
	int myindex;
	int mytype;
	virtual int is_pebble(int all, int pebble);
};

class Magnet : public Parameter
{

public:
	Magnet(char * swn, WireUp * w, int lattice_index,
		int field_index, int foil_index, 
		int type_index, int pebble_index, int myindex);

	virtual ~Magnet();

	int *fam;
	int nfam;
	Wires *mat;
	int family_head;
	int family_polar;
	int lattice_index;
	int field_index;
	int foil_index;
	int type_index;

	double s2i(WireUp * ww, int step, double s);
	double i2s(WireUp * ww, int step, double c);
	double s2didb(WireUp * ww, int step, double s);
	const int * famIndices(int & n); // caller must *not* delete array
	int is_pebble(int all, int pebble);


};


class RfParm : public Parameter
{
public:
	RfParm(char * swn, WireUp * w, int pebble_index, int myindex);
	class CWire * mat;
};










class ParameterD
{
public:
	ParameterD(Parameter * pa, WireUp * w);
	virtual ~ParameterD();

	Parameter * pa;
	WireUp * w;

	double * want_strength;
	double * trim_strength;
	double * summ_strength;
	Spline * want_spline;
	Spline * trim_spline;

	virtual void make_spline(int n, double * x);
	virtual void clear_all(int step);
	virtual void insert_step(int i, int n);
	virtual void delete_step(int i, int n);
	virtual void copy_step(int t, int f);
	virtual void fill_step(int i, double p);

	virtual void set_want_strength(int step, double s);
	virtual void set_trim_strength(int step, double s);

	virtual void set_want_strengthX(int step, double s);
	virtual void set_trim_strengthX(int step, double s);
};

class RfParmD : public ParameterD
{
public:
	RfParmD(RfParm * r, WireUp * w);
	virtual ~RfParmD();
	RfParm * r;
};



class MagnetD : public ParameterD
{
public:
	MagnetD(Magnet * m, WireUp * w);
	virtual ~MagnetD();
	Magnet * m;

	virtual void clear_all(int step);
	virtual void insert_step(int i, int n);
	virtual void delete_step(int i, int n);
	virtual void copy_step(int t, int f);
	virtual void fill_step(int i, double p);

	virtual void set_want_strength(int step, double s);
	virtual void set_trim_strength(int step, double s);


	double * want_current;
	double * trim_current;
	double * summ_current;
	double * summ_currentp;


	void set_want_current(int step, double c);
	void set_trim_current(int step, double c);
	void set_summ_current(int step, double c);
	void xset_want_current(int step, double c);
	void xset_trim_current(int step, double c);
	void xset_summ_current(int step, double c);
	void mangle(int step);
	void fix(int step, int what);
	double get_want_current(int step);
	double get_trim_current(int step);
	void set_currentp(int step, double konstq, double dKonstdt, double pseudo, double pseudoP);
	int recalc1(int step, double p, double gaga, double konstq, double & dPdt);
	void testpp(int step, double dpdt, double ddk);

};
 












class PowerD
{
public:
	PowerD(Power * p, WireUp * w);
	virtual ~PowerD();


	WireUp * w;
	Power * p;

	double * current;
	double * currentp;
        RampFit *fit;
        RampFit *fit2;


	void clear_all(int step);
	void insert_step(int i, int n);
	void delete_step(int i, int n);
	void copy_step(int t, int f);
	int set_current(int step, double c);
	double * make_fit(int step1, int step0); // caller must *not* delete array
	void recalc1(int step);
	void testpp(int step, double dpdt);
};

class Power
{
public:
	Power(char * swn, int init, double max_current, double min_current,
	double max_advice, double min_advice, double i_rating, double v_rating,
	double bit_per_amp_rdb, double mvolts, double bit_per_amp,
	double max_ramp_speed, double resistance, double inductance,
	int christie);

	virtual ~Power();
	char * swn;
	char model[20];
	char polar[20];
	int init;
	double max_current;
	double min_current;
	double max_advice;
	double min_advice;
	double i_rating;
	double v_rating;
	double bit_per_amp;
	double bit_per_amp_rdb; // same value, maybe different sign
	double mvolts;
        double max_ramp_speed;                                             
        double resistance;
        double inductance;
        int christie;

        Wires *mat;
};      



class Cavity
{
public:
	Cavity(char * swn);
	virtual ~Cavity();
	char * swn;
	WireUp *w;
	class CWire * mat;
	double min_volts;
	double max_volts;
};

class CavityD
{
public:
	CavityD(Cavity * c,WireUp *w);
	virtual ~CavityD();
	Cavity * c;
	WireUp *w;
	double * voltage;
	void clear_all(int step);
	void copy_step(int t, int f);
	void insert_step(int, int);
	void delete_step(int, int);
};


class CWire
{
public:
	WireUp *w;
};


class CWire1 : public CWire
{
public:
	CWire1(int nc1, int nc2, int ncom);
	virtual ~CWire1();
	int parm1, parm2;
	int nc1, nc2, ncom;
	int * cav1;
	int * cav2;
	int * common;
};


class CWire2 : public CWire
{
public:
	CWire2(int nc1);
	virtual ~CWire2();
	int parm1;
	int nc1;
	int *cav1;
};

class CRadius : public CWire
{
public:
	int parm1;
};

class CFreq : public CWire
{
public:
	int parm1;
};
