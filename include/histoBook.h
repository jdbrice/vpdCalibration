
#ifndef HISTOBOOK_H
#define HISTOBOOK_H

#include "allroot.h"
#include <map>
#include <string>
#include <sstream>

using namespace std;

class histoBook {

private:

	string currentDir;

	std::map<string, TH1*> book;
	
	string filename;
	
	TFile *file;


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

	TDirectory* getDirectory( ) { return gDirectory; }

	void save();

};



#endif