
#include "reporter.h"
#include "utils.h"

int reporter::instances = 0;
using namespace jdbUtils;

reporter::reporter( string filename, int w, int h ){

	this->filename = filename;

	canvas = new TCanvas( ("reporter"+ts(instances)).c_str() , "canvas", w, h);
	canvas->Print( ( filename + "[" ).c_str() );
	instances++;
}

reporter::~reporter() {
	// properly close the report file
	canvas->Print( ( filename + "]" ).c_str() );
	delete canvas;
}

void reporter::newPage( int dx, int dy){

	this->dx = dx;
	this->dy = dy;

	canvas->Clear();
	canvas->Divide( dx, dy );
	cd( 1 );
}

void reporter::cd( int pad ){
	canvas->cd( pad );
	currentPad = pad;
}

void reporter::cd( int x, int y){
	cd( ( (y - 1 ) * dx) + x );
}

void reporter::next(){
	currentPad++;
	if ( currentPad > dx*dy){
		savePage();
		newPage( dx, dy);
	} else {
		cd( currentPad );
	}
}

void reporter::savePage( ) {
	canvas->Print( ( filename ).c_str() );
}