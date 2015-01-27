

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
	 * Loop through events
	 */
	for ( int i = 0; i < nEvents; i++ ){

		tm.getEntry( i );


		for ( int j= 0; j < 19; j ++ ){

			if ( j < 16 ){
				book.fill( "bbqAdc", j, tm.get( "vpdBbqAdcWest", j ) );
				book.fill( "bbqTdc", j, tm.get( "vpdBbqTdcWest", j ) );
			}

			book.fill( "le", j, tm.get( "vpdLeWest", j ) );
			book.fill( "tot", j, tm.get( "vpdTotWest", j ) );			
		}
	}



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

	logger.info() << "Finished in " << tt.elapsedTime() << endl;



}