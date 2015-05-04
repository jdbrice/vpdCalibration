
#include <fstream>


void nudgeVpd( string vpdParams, string qaFile, double fitLow = -10.0, double fitHigh = 10.0 ){

	Reporter rp( "VpdNudge.pdf" );
	rp.newPage();

	ifstream infile;
	infile.open( vpdParams.c_str() );

	ofstream outfile;
	outfile.open( "vpd_4DB.nudged.dat" );

	TFile * qa = new TFile( qaFile.c_str() );

	TH1D * dz = (TH1D*) qa->Get( "zvertexDelta" );
	dz->GetXaxis()->SetRangeUser( fitLow * 1.5, fitHigh * 1.5 );

	TF1 *g = new TF1( "gaus", "gaus" );
	g->SetRange( fitLow, fitHigh );

	dz->Draw();
	gPad->SetLogy(1);
	dz->Fit( g, "R" );
	rp.savePage();

	double mean = g->GetParameter( 1 );
	cout << "Adjust mean by " << -mean << " [cm]" << endl;
	const Float_t c_light = 29.9792458;
	double dt = 2 * (mean) / c_light;
	cout << "dt = " << dt << " [ns]" << endl;


	/**
	 * Read in and write out new
	 */
	
	for ( int i = 0; i < 38; i++ ){

		int channel = -1, nBins = -1;
		infile >> channel;
		infile >> nBins;

		outfile << channel << endl;
		outfile << nBins << endl;

		//copy the bin edges
		for ( int iBin = 0; iBin <= nBins; iBin ++ ){

			double bv = -1;
			infile >> bv;
			outfile << bv << " ";
		}
		outfile << endl;

		// copy the corrections
		for ( int iBin = 0; iBin <= nBins; iBin ++ ){

			double bv = -1;
			infile >> bv;
			if ( channel < 20 ) // only apply correction to one side
				outfile << (bv + dt) << " ";
			else 
				outfile << bv << " ";
		}
		outfile << endl;

	} // loop on channels




	infile.close();
	outfile.close();


	delete qa;
	delete g;



}