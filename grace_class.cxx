#include "grace_class.hxx"
#ifdef FOOTER
@    string on
@    string loctype view
@    string 0.08, 0.060
@    string color 1
@    string rot 0
@    string font 0
@    string just 0
@    string char size 1.000000
@    string def "01234567890123456789012345678901234567890123456789012345678901234567890123456789"
#endif


GRACE::Curve::Curve()
{
	name = "";
	dots = 0;
}

GRACE::Curve::Curve(const char * name) 
{
	if(name) this->name = name;
	dots = 0;
}

GRACE::Curve::~Curve()
{
}

bool GRACE::Curve::greaterX(Point an, Point  am)
{
      		return( an.x < am.x) ;
}

bool GRACE::Curve::greaterY(Point an, Point  am)
{
      		return( an.y < am.y) ;
}

double GRACE::Curve::maxX()
{
	if( data.empty() ) return -1e33;
	PointList::iterator  j = std::max_element(data.begin(), data.end(), greaterX);
	return j->x;;
}

double GRACE::Curve::minX()
{
	if( data.empty() ) return 1e33;
	PointList::iterator  j = std::min_element(data.begin(), data.end(), greaterX);
	return j->x;;
}

double GRACE::Curve::maxY()
{
	if( data.empty() ) return -1e33;
	PointList::iterator  j = std::max_element(data.begin(), data.end(), greaterY);
	return j->y;;
}

double GRACE::Curve::minY()
{
	if( data.empty() ) return 1e33;
	PointList::iterator  j = std::min_element(data.begin(), data.end(), greaterY);
	return j->y;;
}

void GRACE::Curve::sortX()
{
	std::sort(data.begin(), data.end(), greaterX);
}

void GRACE::Curve::addPoint(double x, double y)
{
	Point d(x,y);
	data.push_back(d);
}

void GRACE::Curve::scale(double s)
{
	for(PointList::iterator j=data.begin(); j != data.end(); j++)
	{
		j->y *= s;;
	}
}

void GRACE::Curve::write(FILE * fp)
{
	for(PointList::iterator j=data.begin(); j != data.end(); j++)
	{
		Point jj = *j;
		fprintf(fp, "%e %e\n", jj.x, jj.y);
	}
	fprintf(fp, "&\n");
}




GRACE::GRACE(const char * filename, const char * title, const char * xaxis, const char * yaxis)
{
	this->filename = filename;
	this->title = title;
	this->xaxis = xaxis;
	this->yaxis = yaxis;
	dozoom = 0;

}

GRACE::~GRACE()
{
	FILE * fp = fopen(filename.c_str(), "w");
	if(!fp)
	{
		perror(filename.c_str());
		exit(-1);
	}
	fprintf(fp,"@title \"%s\"\n", title.c_str());
	fprintf(fp,"@frame linewidth 2\n");
	fprintf(fp,"@xaxis label \"%s\"\n", xaxis.c_str());
	fprintf(fp,"@yaxis label \"%s\"\n", yaxis.c_str());
	fprintf(fp,"@map color 15 to (255, 255, 0), \"yellow\"\n");
	fprintf(fp,"@map color 5 to (0, 139, 0), \"green4\"\n");
	if(subtitle.size())
	{
		fprintf(fp, "@with g0\n");
		fprintf(fp, "@subtitle \"%s\"\n", subtitle.c_str());
	}
	if(dozoom)
	{
		fprintf(fp,"@g0 on\n");
		fprintf(fp,"@with g0\n");
		fprintf(fp,"@    world xmin %f\n", xmin);
		fprintf(fp,"@    world xmax %f\n", xmax);
		fprintf(fp,"@    world ymin %f\n", ymin);
		fprintf(fp,"@    world ymax %f\n", ymax);
	}

	for(int i=0; i< curveList.size(); i++)
	{
		if(curveList[i].name.size())
		{
       			fprintf(fp, "@    s%d legend \"%s\"\n",i, curveList[i].name.c_str() );
		}

		if(curveList[i].dots)
		{
       			fprintf(fp, "@with g0\n");
       			fprintf(fp, "@    s%d symbol 9\n",i);
       			fprintf(fp, "@    s%d symbol size 0.060000\n",i);
       			fprintf(fp, "@    s%d line type 0\n",i);
			
		}
		curveList[i].write(fp);
	}

	fprintf(fp, "@hardcopy device \"PNG\"\n");
	fprintf(fp, "@device \"PNG\" DPI 1200\n");

	fclose(fp);
}

int GRACE::addCurve(const char * setName)
{
	Curve g(setName);
	curveList.push_back(g);
	return(curveList.size()-1);
}

void GRACE::mkdots(int g)
{
	curveList[g].dots=1;
}

void GRACE::addPoint(int g, double x, double y)
{
	curveList[g].addPoint(x, y);
}

void GRACE::addPoints(int g, double * x, double * y, int n, int skip)
{
	for(int i=0; i<n; i += skip)
		curveList[g].addPoint(x[i], y[i]);
}

double GRACE::maxX()
{
	double maxx=-1e33;
	for(CurveList::iterator j=curveList.begin(); j != curveList.end(); j++)
	{
		double max = j-> maxX();
		if(maxx < max) maxx = max;
	}
	return maxx;
}

double GRACE::minX()
{
	double minx= 1e33;
	for(CurveList::iterator j=curveList.begin(); j != curveList.end(); j++)
	{
		double min = j-> minX();
		if(minx > min) minx = min;
	}
	return minx;
}

double GRACE::maxY()
{
	double maxy=-1e33;
	for(CurveList::iterator j=curveList.begin(); j != curveList.end(); j++)
	{
		double max = j-> maxY();
		if(maxy < max) maxy = max;
	}
	return maxy;
}

double GRACE::minY()
{
	double miny= 1e33;
	for(CurveList::iterator j=curveList.begin(); j != curveList.end(); j++)
	{
		double min = j-> minY();
		if(miny > min) miny = min;
	}
	return miny;
}

void GRACE::scale(int g, double s)
{
	curveList[g].scale(s);
}

void GRACE::sortX(int g)
{
	curveList[g].sortX();
}

void GRACE::sortX()
{
	for(CurveList::iterator j=curveList.begin(); j != curveList.end(); j++)
	{		
		j->sortX();
	}
}

void GRACE::subTitle(const char * s)
{
	subtitle = s;
}

void GRACE::zoom(double xmin, double xmax, double ymin, double ymax)
{
	dozoom = 1;
	this->xmax = xmax;
	this->xmin = xmin;
	this->ymax = ymax;
	this->ymin = ymin;
}


