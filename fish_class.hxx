//	Class to read one or many files from the superfish sf7 post processor.
//	The data are assembled into one array. All length are converted to meters.
//
//	The class identifies automatically the file format and substitutes
//      the D in the fortran exponent to an E when reading the data.
//
//	some computers do not handle thewindows cr+lf and attach instead a zero at the end of the line.
//      The fix for that is not included so far.




class FISH
{
	static const int tausend = 300;
	int nfiles;
	char  files[tausend][150];
	FILE * fp;
	int oldnames;
// 	void readln(char * line, FILE * fp);
	double readDbl();
public:
	double Emax;
	double Ezmax;
	double Ermax;
	double Bmax;
	double rmin, rmax, zmin, zmax, freq;
	int rstep, zstep;
	double * rmina, * rmaxa, * zmina, * zmaxa, * freqa; // values for each file
	int * rstepa, * zstepa; // values for each file
	double dr;
	double dz;

        double * zloc;
        double * rloc;
        double * ezero;
        double **ez ;
        double **er ;
        double **ee ;
        double **bt ;

	FISH(const char cavName[]);
	virtual ~FISH();
	
};
