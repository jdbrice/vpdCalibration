
#include <fstream>


void adjustVpdZ( double currentVpdVzMean, string vpdParams  ){

	double dt = currentVpdVzMean / ( 29.9 ) * 2;

	ifstream infile;
	infile.open( vpdParams.c_str() );

	ofstream outfile;
	outfile.open( "vpd_4DB.vz_Off.dat" );

	cout << "Adjust Time by " << dt << " [ns]" << endl;
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
			double tdt = dt;
			infile >> bv;
			if ( 0.0 == bv || i <= 19 ){ // ONLY add to WEST corrections
				tdt = 0.0;
			}
			outfile << (bv + tdt) << " ";
		}
		outfile << endl;

	} // loop on channels




	infile.close();
	outfile.close();

}