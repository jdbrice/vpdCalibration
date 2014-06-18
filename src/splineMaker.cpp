

#include "splineMaker.h"


splineMaker::splineMaker( const vector< double > &x, const vector< double > &y, Interpolation::Type type ){
	spline = NULL;
	spline = new Interpolator( x, y, type);

	domainMin = x[ 0 ];
	domainMax = x[ x.size() - 1 ];

}

splineMaker::splineMaker( TH1D* hist, int place, Interpolation::Type type , int firstBin , int lastBin ){
	spline = NULL;
	if ( !hist ) 
		return;

	int nBins = hist -> GetNbinsX();
	int fBin = firstBin;
	int lBin = lastBin;

	
	if ( lBin <= -1 || lBin >= nBins )
		lBin = nBins;
	if ( fBin < 0 || fBin >= nBins )
		fBin = 1;

	vector<double> x((lBin - fBin) + 1 + 2);
	vector<double> y((lBin - fBin) + 1 + 2);


	int j = 1;
	x[ 0 ] = hist->GetBinLowEdge( fBin );
	y[ 0 ] = hist->GetBinContent( fBin );
	for ( int i = fBin; i <= lBin; i++){

		double bEdge = hist->GetBinLowEdge( i );
		double bWidth = hist->GetBinWidth( i );

		double _y = hist->GetBinContent( i );
		double _x = bEdge;
		if ( splineAlignment::left == place )
			_x = bEdge;
		else if ( splineAlignment::center == place )
			_x = (bEdge + (bWidth/2.0));
		else if ( splineAlignment::right == place )
			_x = ( bEdge + bWidth );

		x[ j ] = _x;
		y[ j ] = _y;

		j++;
	}
	x[ j ] = hist->GetBinLowEdge( lBin ) + hist->GetBinWidth( lBin );
	y[ j ] = hist->GetBinContent( lBin );

	domainMin = x[ 0 ];
	domainMax = x[ x.size() - 1 ];

	spline = new Interpolator( x, y, type);

}

splineMaker::~splineMaker(){

	if ( spline )
		delete spline;
	spline = NULL;


}

TGraph * splineMaker::graph( double xmin, double xmax, double step ){
	
	if ( xmin < domainMin )
		xmin = domainMin;
	if ( xmax > domainMax )
		xmax = domainMax;

   	const Int_t n = ( (xmax - xmin ) / step)  ;
   	Int_t i = 0;
   	Float_t xcoord[n], ycoord[n];

   	for ( double xi = xmin+step; xi < xmax; xi += step) { 
      	xcoord[i] = xi;
      	ycoord[i] = spline->Eval(xi);
      	i++; 
   	}
   	
   	TGraph *gr = new TGraph( n, xcoord, ycoord );

	return gr;
}
