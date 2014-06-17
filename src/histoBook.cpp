

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

	
	// make the legend and draw it once to apply styles etc. 
	// for some reason needed to make styling work on the first draw
	legend = new TLegend( 0.65, 0.65, 0.9, 0.9);
	legend->SetFillColor( kWhite );
	legend->Draw();
	legend->Clear();

	globalStyle();

}
// destructor
histoBook::~histoBook(){

	delete legend;

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
string histoBook::cd( string sdir = "/" ){

	string old = currentDir;

	char* csdir = (char*)sdir.c_str();
	file->cd();

	if ( file->GetDirectory( csdir ) ){
		file->cd( csdir );
	} else {
		//cout << "[histoBook.cd] creating directory " << sdir << endl;
		file->mkdir( csdir );
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

void histoBook::make1D( string name, string title, uint nBins, const Double_t* bins  ){

	TH1D* h;
	h = new TH1D( name.c_str(), title.c_str(), nBins, bins );

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
TH2* histoBook::get2D( string name, string sdir  ){
	if ( sdir.compare("") == 0)
		sdir = currentDir;
	return (TH2*)book[ ( sdir  + name  ) ];
}

void histoBook::fill( string name, double bin, double weight ){ 
	if ( get( name ) != 0)
		get( name )->Fill( bin, weight );
}


void histoBook::globalStyle(){

	gStyle->SetCanvasColor(kWhite);     // background is no longer mouse-dropping white
  	gStyle->SetPalette(1,0);            // blue to red false color palette. Use 9 for b/w
  	gStyle->SetCanvasBorderMode(0);     // turn off canvas borders
  	gStyle->SetPadBorderMode(0);
  	gStyle->SetPaintTextFormat("5.2f");  // What precision to put numbers if plotted with "TEXT"


  	// For publishing:
  	gStyle->SetLineWidth(2.);
  	gStyle->SetTextSize(0.7);
  	gStyle->SetLabelSize(0.05,"xy");
  	gStyle->SetTitleSize(0.05,"xy");
  	gStyle->SetTitleOffset(1.0,"x");
  	gStyle->SetTitleOffset(1.5,"y");
  	gStyle->SetPadTopMargin(0.1);
  	gStyle->SetPadRightMargin(0.1);
  	gStyle->SetPadBottomMargin(0.16);
  	gStyle->SetPadLeftMargin(0.2);

  	gStyle->SetFillColor(-1); 
	gStyle->SetFillStyle(4000); 


	
}


histoBook* histoBook::style( string histName ){
	styling = histName;
	return this;
}

histoBook* histoBook::set( string param, string p1, string p2, string p3, string p4 ){
	
	// force the param name to lowercase
	transform(param.begin(), param.end(), param.begin(), ::tolower);

    TH1* h = get( styling );
    if ( h ){

	    if ( "title" == param ){
	    	h->SetTitle( p1.c_str() );
	    } else if ( "x" == param ){
	    	h->GetXaxis()->SetTitle( p1.c_str() );
	    } else if ( "y" == param ){
	    	h->GetYaxis()->SetTitle( p1.c_str() );
	    } else if ( "legend" == param ){
	    	if ( p2 == "")
	    		p2="lpf";
	    	legend->AddEntry( h, p1.c_str(), p2.c_str() );
			legend->Draw();
	    } else if ( "draw" == param ){
	    	drawOption = p1;
	    }
	}

	return this;
}

histoBook* histoBook::set( string param, double p1, double p2, double p3, double p4  ){

	transform(param.begin(), param.end(), param.begin(), ::tolower);

    TH1* h = get( styling );
    if ( h ){

	    if ( "linecolor" == param ){

	    	h->SetLineColor( (int) p1 );
	    } else if ( "domain" == param ){
	    	double min = p1;
	    	double max = p2;
		    h->GetXaxis()->SetRangeUser( min, max );
	    } else if ( "dynamicdomain" == param ){
	    	double thresh = p1;
	    	int min = (int)p2;
	    	int max = (int)p3;
	    	int axis = (int)p4;		// 1 = x, 2 = y

	    	if ( 1 != axis && 2 != axis )
	    		axis = 1;
	    	
	    	if ( thresh >= 0) {
	    		if ( -1 >= min )
	    			min = h->FindFirstBinAbove( thresh, axis );
	    		if ( -1 >= max )
	    			max = h->FindLastBinAbove( thresh, axis );
	    	}
	    	
	    	if ( 1 == axis )
		  	  h->GetXaxis()->SetRange( min, max );
		  	else if ( 2 == axis )
		  		h->GetYaxis()->SetRange( min, max );

	    }  else if ( "range" == param ){

	    	double min = p1;
	    	double max = p2;
	    	
	    	h->GetYaxis()->SetRangeUser( min, max );
	    } else if ( "markercolor" == param ) {
	    	h->SetMarkerColor( (int)p1 );
	    } else if ( "markerstyle" == param ) {
	    	h->SetMarkerStyle( (int)p1 );
	    } else if ( "legend" == param ){
	    	// p1 - alignmentX
	    	// p2 - alignmentY
	    	// p3 - width
	    	// p4 - height

	    	// make sure option is valid
	    	if ( !(legendAlignment::center == p1 || legendAlignment::left == p1 || legendAlignment::right == p1) )
	    		p1 = legendAlignment::best;
	    	if ( !(legendAlignment::center == p2 || legendAlignment::top == p2 || legendAlignment::bottom == p2) )
	    		p2 = legendAlignment::best;
	    	placeLegend( p1, p2, p3, p4 );
	    }

    }
    
    
    


	return this;
}


histoBook* histoBook::draw(string name, Option_t* opt ){

	// no parameters
	if ( name == ""){
		TH1* h = get( styling );
		if ( h ){
			// use the draw option set in its styling
			h->Draw( drawOption.c_str() );
			drawOption = "";
		}	
	} else {
		TH1* h = get( name );
		if ( h ){
			h->Draw( opt );
		}
	}
	
	return this;
}


histoBook* histoBook::placeLegend( int alignmentX, int alignmentY, double width, double height ){

	double mR = 1 - gPad->GetRightMargin();
	double mL = gPad->GetLeftMargin();
	double mT = 1- gPad->GetTopMargin();
	double mB = gPad->GetBottomMargin();

	double x1, x2, y1, y2;

	if ( width <= 0 || width > 1 )
		width = .2;
	if ( height <= 0 || height > 1 )
		height = .2;

	// alignment best needs a current histo
	if ( !(get( styling )) ){
		if ( legendAlignment::best == alignmentX )
			alignmentX = legendAlignment::right;
		if ( legendAlignment::best == alignmentY )
			alignmentY = legendAlignment::top;
	} else {

		//TODO

	}


	if ( 	legendAlignment::left ==  alignmentX ){
		x1 =  mL ;
		x2 =  mL + width;
	}
	if ( 	legendAlignment::right ==  alignmentX ){
		x1 =  mR - width;
		x2 =  mR ;
	}
	if ( 	legendAlignment::center ==  alignmentX ){
		x1 =  0.55 - width/2.0;
		x2 =  0.55 + width/2.0;
	}
	if ( 	legendAlignment::top ==  alignmentY ){
		y1 =  mT - height;
		y2 = mT ;
	}
	if ( 	legendAlignment::bottom ==  alignmentY ){
		y1 =  mB ;
		y2 =  mB + height;
	}
	if ( 	legendAlignment::center ==  alignmentY ){
		y1 =  0.55 - height/2.0;
		y2 =  0.55 + height/2.0;
	}
	legend->SetX1NDC( x1 );
	legend->SetX2NDC( x2 );
	legend->SetY1NDC( y1 );
	legend->SetY2NDC( y2 );

	return this;
}