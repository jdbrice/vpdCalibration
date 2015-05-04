






void compareBill( string inName = "params.dat" ){

	HistoBook * book = new HistoBook( "data.root" );
	ifstream infile;
	
	infile.open( inName.c_str() );
	
	if ( infile.is_open() ){
		vector<double> bins[ 38 ];
		vector<double> corr[ 38 ];
		for ( int i = 0 ; i < 38; i++ ){

			int channel = -1;
			infile >> channel;
			cout << "Channel : " << channel << endl;

			int tBins = 0;
			infile >> tBins;
			cout << "tBins " << tBins << endl;


			book->make1D( "file"+ts(0)+"channel"+ts(channel-1), "", tBins, bins[ channel-1 ].data() );

			if ( 	channel >= 1 && channel <= 38  ){

				double tmp = 0;
				for ( int j=0; j < tBins; j++ ){
					infile >> tmp;
					bins[ channel - 1 ].push_back( tmp );
					if ( channel == 1 )
						cout << " bin[ " << j << " ] = " << tmp << endl;
				}
				for ( int j=0; j < tBins; j++ ){
					infile >> tmp;
					corr[ channel - 1 ].push_back( tmp );
					//book->get( "file"+ts(fi)+"channel"+ts(channel-1) )->SetBinContent( j+1, tmp );
				}

				/*string title = "Channel "+ts(channel)+" Slewing Corrections; " + xLabel + ";" + yLabel ;
				// create a hist
				book->cd( "params" );
				book->make1D( "file"+ts(fi)+"channel"+ts(channel-1), title, tBins, totBins[ channel-1] );


				
				
				for ( int j=0; j < tBins; j++ ){ 
					book->get( "file"+ts(fi)+"channel"+ts(channel-1) )->SetBinContent( j+1, correction[ channel - 1 ][ j ]  );
				}*/

				
			} else {
				cout << "[calib." << __FUNCTION__ <<  "] " << "Bad file format, cannot load parameters" << endl;
				good = false;

			} 


		}
	} else {
		good = false;
	}

	infile.close();
}