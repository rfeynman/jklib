#include "jkLib/knobcontainer.hxx"


KnobContainer::KnobContainer()
{
	nKnobs=nKnobListener = 0;
}


KnobContainer::~KnobContainer()
{
}



int KnobContainer::load(const char * filename)
{
	FILE * fp = fopen(filename, "r");
	if(!fp)
	{
		return errno;
	}
	fscanf(fp,"%d", &nKnobs);
	for(int n=0; n<nKnobs; n++)
		fscanf(fp,"%s%lf%lf%s",
			label[n], min+n, max+n, unit[n]);

	fscanf(fp,"%d", &nKnobListener);
	for(int n=0; n<nKnobListener; n++)
		fscanf(fp,"%s%d%lf",
			swn[n], knobIndex+n, constant+n);
	fclose(fp);
	return 0;
}


int KnobContainer::save(const char * filename)
{
	FILE * fp = fopen(filename, "w");
	if(!fp)
	{
		return errno;
	}
	fprintf(fp,"%d\n", nKnobs);
	for(int n=0; n<nKnobs; n++)
		fprintf(fp,"%s %lf %lf %s\n",
			label[n], min[n], max[n], unit[n]);

	fprintf(fp,"%d\n", nKnobListener);
	for(int n=0; n<nKnobListener; n++)
		fprintf(fp,"%s %d %lf\n",
			swn[n], knobIndex[n], constant[n]);
	fclose(fp);
	return 0;
}



int KnobContainer::addKnob(const char * name)
{
	for(int n=0; n<nKnobs; n++)
	{
		if( ! strcmp(label[n], name) ) return n;
	}
	if(nKnobs >= maxKnobs) return -1;
	strcpy(label[nKnobs], name);
	max[nKnobs] = 1;
	min[nKnobs] = -1;
	strcpy(unit[nKnobs], "things");
	int n=nKnobs++;
	return n;
}


int KnobContainer::findKnob(const char * name)
{
	for(int i=0; i<nKnobs; i++)
	{
		if( ! strcmp(label[i], name) ) return i;
	}
	return -1;
}


void KnobContainer::deleteKnob(int kIndex)
{
	clearKnob(kIndex);
	for(int i = kIndex+1; i<nKnobs; i++)
	{
		strcpy(label[i-1], label[i]);
		max[i-1] = max[i];
		min[i-1] = min[i];
		strcpy(unit[i-1], unit[i]);
	}
	nKnobs--;
}


void KnobContainer::deleteAllKnobs()
{
	nKnobs=nKnobListener = 0;
}


void KnobContainer::clearKnob(int kIndex)
{
	int j=0;
	for(int i=0; i<nKnobListener; i++)
	{
		if(knobIndex[i] != kIndex)
		{
			strcpy(swn[j], swn[i]);
			knobIndex[j] = knobIndex[i];
			constant[j++] = constant[i];
		}
	}
	nKnobListener =j;
}



void KnobContainer::set_knobMin(int kIndex, double minimum)
{
	if(kIndex >= 0 && kIndex < nKnobs) min[kIndex] = minimum;
}


double KnobContainer::get_knobMin(int kIndex)
{
	if(kIndex >= 0 && kIndex < nKnobs) return min[kIndex];
	else return 0.;
}



void KnobContainer::set_knobMax(int kIndex, double maximum)
{
	if(kIndex >= 0 && kIndex < nKnobs) max[kIndex] = maximum;
}


double KnobContainer::get_knobMax(int kIndex)
{
	if(kIndex >= 0 && kIndex < nKnobs) return max[kIndex];
	else return 0.;
}



void KnobContainer::set_knobUnit(int kIndex, const char * un)
{
	if(kIndex >= 0 && kIndex < nKnobs) strcpy(unit[kIndex], un);
}


const char * KnobContainer::get_knobUnit(int kIndex)
{
	if(kIndex >= 0 && kIndex < nKnobs) return unit[kIndex];
	else NULL;
}



char ** KnobContainer::knobNames()
{
	char ** l = new char * [nKnobs+1];
	for(int i=0; i< nKnobs; i++)
		l[i] = label[i];
	l[nKnobs]=NULL;
	return l;
}


double * KnobContainer::knobMins()
{
	return min;
}


double * KnobContainer::knobMaxs()
{
	return max;
}



int KnobContainer::knobAddMagnet(int kIndex, char * magnet)
{
	for(int i=0; i<nKnobListener; i++)
	{
		if(kIndex == knobIndex[i] && !strcmp(magnet, swn[i]) )
			return i;
	}
	if(nKnobListener >= maxKnobListeners) return -1;
	strcpy(swn[nKnobListener], magnet);
	knobIndex[nKnobListener] = kIndex;
	constant[nKnobListener] = 1;
	int n = nKnobListener++;
	return n;
}


int KnobContainer::knobFindMagnet(int kIndex, char * magnet)
{
	for(int i=0; i<nKnobListener; i++)
	{
		if(kIndex == knobIndex[i] && !strcmp(magnet, swn[i]) )
			return i;
	}
	 return -1;
}


void KnobContainer::knobDeleteMagnet(int mIndex)
{
	for(int i=mIndex+1; i<nKnobListener; i++)
	{
		strcpy(swn[i-1], swn[i]);
		knobIndex[i-1] = knobIndex[i];
		constant[i-1] = constant[i];
	}
	nKnobListener--;
}



void KnobContainer::set_knobConst(int mIndex, double cons)
{
	if(mIndex >= 0 && mIndex < nKnobListener) constant[mIndex] = cons;
}


double KnobContainer::get_knobConst(int kIndex, int mIndex)
{
	if(mIndex >= 0 && mIndex < nKnobListener) return constant[mIndex];
}



int KnobContainer::nKnobMags(int kIndex)
{
	int n=0;
	for(int i=0; i<nKnobListener; i++)
	{
		if(knobIndex[i] == kIndex) n++;
	}
	return n;
}


char ** KnobContainer::knobMagNames(int kIndex)
{
	int n=nKnobMags(kIndex);
	char ** s = new char*[n+1];
	int j=0;
	for(int i=0; i<nKnobListener; i++)
	{
		if(knobIndex[i] == kIndex) s[j++] = swn[i];
	}
	s[j] = NULL;
	return s;
}


double * KnobContainer::knobMagConst(int kIndex)
{
	int n=nKnobMags(kIndex);
	double * s = new double[n];
	int j=0;
	for(int i=0; i<nKnobListener; i++)
	{
		if(knobIndex[i] == kIndex) s[j++] = constant[i];
	}
	return s;
}



