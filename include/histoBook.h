
#ifndef HISTOBOOK_H
#define HISTOBOOK_H

#include "allroot.h"
#include <map>
#include <string>
#include <sstream>
#include <stdarg.h>

using namespace std;

class legendAlignment {
public:
	static const int topLeft = 1;
	static const int topRight = 2;
	static const int topCenter = 3;
	static const int bottomLeft = 4;
	static const int bottomRight = 5;
	static const int bottomCenter = 6;

};

class histoBook {

private:

	string currentDir;

	std::map<string, TH1*> book;
	
	string filename;
	
	TFile *file;

	string styling;

	TLegend * legend;


public:

	histoBook( string name );
	~histoBook();

	string cd( string dir, bool subd = false );
	void add( string name, TH1 * );
	TH1* get( string name, string sdir = "" ) ;
	TH1* get( stringstream sstr, string sdir = "") { return get( sstr.str(), sdir); }
	void fill( string name, double bin, double weight = 1);
	void make1F( string name, string title, uint nBins, double low, double hi );
	void make1D( string name, string title, uint nBins, double low, double hi );
	void make2D( 	string name, string title, 
					uint nBinsX, double lowX, double hiX, uint nBinsY, double lowY, double hiY );
	void make2D( 	string name, string title, 
					uint nBinsX, const Double_t* xBins, uint nBinsY, double lowY, double hiY );

	TLegend* getLegend() { return legend; }

	histoBook* draw(string name, Option_t* opt= "", bool leg = false );
	histoBook* draw( Option_t* opt= "", bool leg = false );
	histoBook* clearLegend() { legend->Clear(); return this; };
	histoBook* placeLegend( int alignment, double width = -1, double height = -1 );
	

	TDirectory* getDirectory( ) { return gDirectory; }

	void save();

	histoBook* style( string hName );
	//histoBook* set( string param, ... );
	histoBook* set( string param, string val );
	histoBook* set( string param, double p1, double p2 = -1, double p3 = -1, double p4 = -1 );




private:
	void globalStyle();

};



#endif