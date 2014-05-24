

#include "histoBook.h"

// constructor sets the name of the file for saving
histoBook::histoBook( string name ){

	if (name.find(  ".root") != std::string::npos){
		filename = name;	
	} else
		filename = name + ".root";

	currentDir = "/";

	file = new TFile( filename.c_str(), "recreate" );
	file->cd();


}
// destructor
histoBook::~histoBook(){

	save();
	file->Close();
}
void histoBook::save() {

	file->Write();
}


void histoBook::add( string name, TH1* h ){

	string oName = name;
	if ( name.length() <= 1 || !h )
		return;

	name = currentDir + name;
	
	// dont allow duplicated name overites
	if ( book[ name ] ){
		cout << "[histoBook.add] Duplicate histogram name in this directory " << currentDir << " / " << oName << endl;
		return;
	}

	// save the histo to the map
	book[ name ] = h;

}
/*
*
* TODO:: add support for full subdirectory trees
*/
string histoBook::cd( string sdir = "/", bool subd ){

	string old = currentDir;

	char* csdir = (char*)sdir.c_str();
	file->cd();

	if ( file->GetDirectory( csdir ) ){
		file->cd( csdir );
	} else {
		//cout << "[histoBook.cd] creating directory " << sdir << endl;
		TDirectory* d = file->mkdir( csdir );
		file->cd( csdir );

	}

	currentDir = sdir;

	return old;
}

void histoBook::make1F( string name, string title, uint nBins, double low, double hi  ){

	TH1F* h;
	h = new TH1F( name.c_str(), title.c_str(), nBins, low, hi );

	this->add( name, h );
}


void histoBook::make1D( string name, string title, uint nBins, double low, double hi  ){

	TH1D* h;
	h = new TH1D( name.c_str(), title.c_str(), nBins, low, hi );

	this->add( name, h );
}
void histoBook::make2D( string name, string title, uint nBinsX, double lowX, double hiX, uint nBinsY, double lowY, double hiY ){

	TH2D* h;

	h = new TH2D( name.c_str(), title.c_str(), nBinsX, lowX, hiX, nBinsY, lowY, hiY );

	this->add( name, h );
}
void histoBook::make2D( string name, string title, uint nBinsX, const Double_t* xBins, uint nBinsY, double lowY, double hiY ){

	TH2D* h;
	h = new TH2D( name.c_str(), title.c_str(), nBinsX, xBins, nBinsY, lowY, hiY );

	this->add( name, h );
}


TH1* histoBook::get( string name, string sdir  ){
	if ( sdir.compare("") == 0)
		sdir = currentDir;
	return book[ ( sdir  + name  ) ];
}

void histoBook::fill( string name, double bin, double weight ){ 
	if ( get( name ) != 0)
		get( name )->Fill( bin, weight );
}



