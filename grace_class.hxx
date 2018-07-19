#ifndef _GRACE_INCLUDED_ 
#define _GRACE_INCLUDED_ 

#include <stdio.h>
#include <stdlib.h>
#include <error.h>
#include <vector>
#include <algorithm>
#include <string>

// public:
// 	GRACE(const char * filename, const char * title, const char * xaxis, const char * yaxis)
// 	void subTitle(const char * s)
// 	void zoom(double xmin, double xmax, double ymin, double ymax)
// 	int addCurve(const char * setName)
// 	void addPoint(int g, double x, double y)
// 	void addPoints(int g, double * x, double * y, int n, int skip=1)
// 	double maxX()
// 	double minX()
// 	double maxY()
// 	double minY()
// 	void scale(int g, double s)
// 	void sortX(int g)
// 	void sortX()

class GRACE
{
	int dozoom;
	double xmin, xmax, ymin, ymax;

	class Curve;
	class Point
	{
		friend class Curve;
		double x, y;
		Point(double _x, double _y) { x = _x; y = _y;}
	};
	typedef std::vector<Point> PointList;


	class Curve
	{
		friend class GRACE;
		PointList data;
		std::string name;
		int dots;

		static bool greaterX(Point an, Point  am);
		static bool greaterY(Point an, Point  am);
		double maxX();
		double minX();
		double maxY();
		double minY();
		void sortX();
		void addPoint(double x, double y);
		void scale(double s);
		void write(FILE * fp);
public:
		Curve();
		Curve(const char * name);
		~Curve();
	};

	typedef std::vector<Curve> CurveList;
	CurveList curveList;



public:
	std::string title;
	std::string subtitle;
	std::string filename;
	std::string xaxis;
	std::string yaxis;
	GRACE(const char * filename, const char * title, const char * xaxis, const char * yaxis);
	virtual ~GRACE();

	int addCurve(const char * setName=0);
	void mkdots(int g);
	void addPoint(int g, double x, double y);
	void addPoints(int g, double * x, double * y, int n, int skip=1);
	double maxX();
	double minX();
	double maxY();
	double minY();
	void scale(int g, double s);
	void sortX(int g);
	void sortX();
	void subTitle(const char * s);
	void zoom(double xmin, double xmax, double ymin, double ymax);
};
#endif
