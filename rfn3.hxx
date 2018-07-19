#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include "complex_class.hxx"
#include "jams.hxx"

class WireUp;



class Mat
{
public:
        double a[5][5];
        double eigvec[5][5];
        Complex clam1;
        Complex clam2;
        double spur1,spur2,lamg1,lamg2;
        void unitA();
        void printA();
        void printMat(const char * what, double b[5][5]);
        int eigenwerte();
        int eigenvector();
        void inver1(double v[][4]);
        int eigenvector(double betax0, double alfax0, double betaz0, double alfaz0, double * dspini=NULL);
        void printTwiss();

        Mat();


};

class Element
{
        double mx00;
        double mx01;
        double mx02;
        double mx10;
        double mx11;
        double mx12;
        double mz00;
        double mz01;
        double mz02;
        double mz10;
        double mz11;
        double mz12;
public:
        char * swn;
        double totl;
        double dlength;
        double mlength;
        double strength;
	double qk, qk45, hx, hz;
        double p1,p2;
	double b2;
        int type;
	int ix;

	double betax0;
	double alfax0;
	double betaz0;
	double alfaz0;
	double dispx0;
	double dispxp0;
	double betax;
	double alfax;
	double betaz;
	double alfaz;
	double dispx;
	double dispxp;
	double betaxm;
	double betazm;
	double dispm;

        void drift(double a[][5], double l);
        void thinquad(double a[][5], double k1l);
        void lintr3(double l, double qk, double hx, double hz);
        void doLintr3(double a[][5]);
        void doMat(double a[][5]);
	void doMatInt(Mat & a, double & mocomp, double & chromx, double & chromz, double & qx, double & qz);
	void getStrength(class WireUp * w, int step);
	void setStrength(class WireUp * w, int step);
        Element * next;
        Element(const char * name, double tl, int ty, double dl, double ml, double st, double b2, int ix);
};


class Opticalc
{
public:
	int ok;
        int d;
        char latnam[30];
        char swn[30];
        char nam[30];
        double s,l;
        Mat a;
        double tl;
        double dl;
	int nelem;
        Element ** line;
        Element * chain;
        Element * last;
	double mocomp;
	double chromx;
	double chromz ;
	double qx;
	double qz;

        Opticalc(const char * file, WireUp * w);
        void calc(WireUp * w, int step);
	void setStrength(WireUp * w, int step);
};

