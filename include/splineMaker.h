#ifndef SPLINE_MAKER_H
#define SPLINE_MAKER_H

#include "allroot.h"

using namespace ROOT::Math;
using namespace std;

class splineAlignment {
public:
	const static int left = 0;
	const static int center = 1;
	const static int right = 2;
};


class splineMaker
{
public:

	// pass through constructor
	splineMaker( const vector< double > &x, const vector< double > &y, Interpolation::Type type = Interpolation::kCSPLINE );
	
	// from histogram
	splineMaker( TH1D* hist, int place = splineAlignment::left, Interpolation::Type type = Interpolation::kCSPLINE, int firstBin = 1, int lastBin = -1 );

	TGraph* graph( double xmin, double xmax, double step );
	//void draw( TH1D* hist, double xmin, double xmax, double step );

	Interpolator* getSpline() { return spline; }

	~splineMaker();

	

private:
	Interpolator* spline;
	double domainMin, domainMax;
};




#endif