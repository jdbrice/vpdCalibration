
#include <iostream>
#include <algorithm>    // std::sort
#include <vector>       // std::vector

using namespace std;

bool found( vector<int> vec, int val ){

	for ( int i = 0; i < vec.size(); i++){
		if ( val == vec[ i ] )
			return true;
	}
	return false;

}


void channelLife( string filename = "cLife.xml" ){


	gStyle->SetOptStat( 0 );

	// get the config file
	// has data location
	// reporter output file
	// HistoBook output file etc.
	XmlConfig * cfg = new XmlConfig( filename );

	
	/**
	 * Setup the logger, reporter, and HistoBook - make all histograms
	 */
	Logger::setGlobalLogLevel( Logger::llAll );

	Logger logger;

	Reporter rp( cfg, "Reporter." );

	HistoBook book( "hbChannelLife.root", cfg );
	book.makeAll( "histograms" );

	string det=cfg->getString( "Detector" , "West" );
	int detOff = 0;
	if ( "East" == det )
		detOff = 19;


	/*
	 * Setup the chain
	 */
	TChain * chain = new TChain( cfg->getString( "input.dst:treeName" ).c_str() );
	ChainLoader::load( chain, cfg->getString( "input.dst:url" ), cfg->getInt( "input.dst:maxFiles", -1 ) );
	TreeMap tm( chain );

	int nEvents = tm.getEntries();

	logger.info() << "nEvents = " << nEvents << endl;
	TaskTimer tt;
	tt.start();

	/**
	 * Loop through events to make a list of run numbers
	 */
	vector<int> runNumbers;
	vector<int> sortedRunNumbers;
	for ( int i = 0; i < nEvents; i++ ){

		tm.getEntry( i );
		int run = tm[ "run" ];
		if ( 16009013 == run ) // this run was marked as junk but i accidentally included it
			continue;
		
		if ( !found( runNumbers, run ) ){
			runNumbers.push_back( run );
			logger.trace() << "Run# " << ts( run ) << endl;
		}
	}

	while( runNumbers.size() > 0 ){

		int minRun = 100000000;
		int index = -1;
		for ( int i = 0; i < runNumbers.size(); i++ ){
			if ( runNumbers[ i ] < minRun ){
				minRun = runNumbers[ i ];
				index = i;
			}
		}

		if ( index >= 0 ){
			sortedRunNumbers.push_back( minRun );
			runNumbers.erase( runNumbers.begin()+index );	
		}
	}

	vector<string> histos = cfg->childrenOf( "histograms" );
	for ( int i = 0; i < histos.size(); i++ ){
		string hName = cfg->getString( histos[ i ] + ":name" );
		for ( int j = 0; j < sortedRunNumbers.size(); j ++ ){
			string rn = ts( sortedRunNumbers[ j ] );
			book.clone( hName, hName + "_" + rn );
		}
		
	}

	book.make2D( "offline", det+" Offline (HPTDC);Run Number;PMT Channel", sortedRunNumbers.size(), 0, sortedRunNumbers.size(), 19, detOff, 19+detOff );
	book.make2D( "online", det + " Online (Bbq);Run Number;PMT Channel", sortedRunNumbers.size(), 0, sortedRunNumbers.size(), 16, 0, 16 );

	book.get2D( "offline" )->SetBit( TH1::kCanRebin );
	book.get2D( "online" )->SetBit( TH1::kCanRebin );

	/**
	 * Loop through events
	 */
	
	for ( int i = 0; i < nEvents; i++ ){

		tm.getEntry( i );
		int run = tm[ "run" ];
		string rn = "_"+ts(run);

		for ( int j= 0; j < 19; j ++ ){

			if ( j < 16 ){
				book.fill( "bbqAdc", j, tm.get( "vpdBbqAdc"+det, j ) );
				book.fill( "bbqTdc", j, tm.get( "vpdBbqTdc"+det, j ) );
			}

			book.fill( "le", j, tm.get( "vpdLe"+det, j ) );
			book.fill( "tot", j, tm.get( "vpdTot"+det, j ) );

			if ( j < 16 ){
				book.fill( "bbqAdc"+rn, j, tm.get( "vpdBbqAdc"+det, j ) );
				book.fill( "bbqTdc"+rn, j, tm.get( "vpdBbqTdc"+det, j ) );
			}

			book.fill( "le"+rn, j, tm.get( "vpdLe"+det, j ) );
			book.fill( "tot"+rn, j, tm.get( "vpdTot"+det, j ) );


		}
	}

	int threshold = 1;
	for ( int j = 0; j < sortedRunNumbers.size(); j ++ ){
		string rn = "_"+ts( sortedRunNumbers[ j ] );

		//cout << "Run " << rn << endl;
		for ( int i = 0; i < 19; i++ ){
			int hits = book.get2D( "tot" + rn )->ProjectionY( "_py", i+1, i+1 )->Integral( 2, -1);

			if ( 5 == i || 13 == i || 14 == i )
				hits = 0;
			if ( hits > threshold ){
				char * binLabel = ts(sortedRunNumbers[ j ]).c_str();
				book.get2D( "offline" )->Fill( binLabel, (int)i+detOff, 1 );
			}

		}
	}

	book.get2D( "offline" )->LabelsDeflate("X");

	for ( int j = 0; j < sortedRunNumbers.size(); j ++ ){
		string rn = "_"+ts( sortedRunNumbers[ j ] );

		//cout << "Run " << rn << endl;
		for ( int i = 0; i < 16; i++ ){
			int hits = book.get2D( "bbqTdc" + rn )->ProjectionY( "_py", i+1, i+1 )->Integral( 10, -1);

			if ( hits > threshold ){
				char * binLabel = ts(sortedRunNumbers[ j ]).c_str();
				book.get2D( "online" )->Fill( binLabel, i, 1 );
			}

		}
	}



	rp.newPage( 1, 2 );
	book.style( "offline" )->set("draw", "colz")->draw();
	rp.next();
	book.style( "online" )->set("draw", "colz")->draw();
	rp.savePage();
	/**
	 * Print out a little report showing each channel's total activity
	 */
	rp.newPage( 1, 2 );
	book.style( "le" )->set( "draw", "colz" )->set( "logZ", 1 )->draw();
	rp.next();
	book.style( "tot" )->set( "draw", "colz" )->set( "logZ", 1 )->draw();
	rp.savePage();

	rp.newPage( 1, 2 );
	book.style( "bbqAdc" )->set( "draw", "colz" )->set( "logZ", 1 )->draw();
	rp.next();
	book.style( "bbqTdc" )->set( "draw", "colz" )->set( "logZ", 1 )->draw();
	rp.savePage();


	for ( int j = 0; j < sortedRunNumbers.size(); j ++ ){
		string rn = "_"+ts( sortedRunNumbers[ j ] );
		

		rp.newPage( 1, 2 );
		book.style( "le"+rn )->set( "draw", "colz" )->set( "logZ", 1 )->
		set("title", ts( sortedRunNumbers[ j ] ) + " le" )->draw();
		rp.next();
		book.style( "tot"+rn )->set( "draw", "colz" )->set( "logZ", 1 )->
		set("title", ts( sortedRunNumbers[ j ] ) + " tot" )->draw();
		rp.savePage();

		rp.newPage( 1, 2 );
		book.style( "bbqAdc"+rn )->set( "draw", "colz" )->set( "logZ", 1 )->
		set("title", ts( sortedRunNumbers[ j ] ) + " bbq Adc" )->draw();
		rp.next();
		book.style( "bbqTdc"+rn )->set( "draw", "colz" )->set( "logZ", 1 )->
		set("title", ts( sortedRunNumbers[ j ] ) + " bbq Tdc" )->draw();
		rp.savePage();


	}


	logger.info() << "Finished in " << tt.elapsedTime() << endl;



}