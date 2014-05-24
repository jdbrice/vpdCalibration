#include <iostream>
#include <fstream>
using namespace std;

void plotCorrections( TString ifn = "cor.dat", int t = 0) {

	TString out = ifn + ".root";
	TFile * f = new TFile( out, "RECREATE" );
	int channel;
	double nBins[ 39 ];
	double* bins[ 39 ];
	double* vals[ 39 ];

	bool good = true;

	ifstream infile;
    infile.open( ifn );
	if ( infile.is_open() ){

		while ( good && !infile.eof()){
			
			infile >> channel;

			int tBins = 0;
			infile >> tBins;
			cout << "Channel[ " << channel << " ] : " << tBins << endl;
			
			if ( t == 0 )
				tBins ++;

			nBins[ channel ] = tBins;
			bins[ channel ] = new double[ tBins ];
			vals[ channel ] = new double[ tBins ];

			cout << "BINS" << endl;
			double tmp = 0;
			for ( int i=0; i < tBins; i++ ){
				infile >> tmp;
				bins[ channel ][ i ] = tmp;
				cout << " " << tmp;
			}
			cout << endl << endl;
			cout << "VALS" << endl;
			for ( int i=0; i < tBins; i++ ){
				infile >> tmp;
				vals[ channel ][ i ] = tmp;
				cout << " " << tmp;
			}
			cout << endl << endl;
		
			if ( channel == 38 ){
				good = false;
			}
		}

	} else {
		cout << "Could not open file \n";
	}


	TH1D* hcor[ 39 ];

	for ( int j = 1; j < 39; j++ ){
		stringstream sstr;
		sstr << "cor" << j;
		hcor[j] = new TH1D( sstr.str().c_str(), sstr.str().c_str(), nBins[ j ] - 1, bins[ j ] );
		for ( int i = 0; i < nBins[ j ]; i++ ){

			hcor[ j ]->SetBinContent( i + 1, vals[ j ][ i ] );

		}
	}


	f->Write();

}
