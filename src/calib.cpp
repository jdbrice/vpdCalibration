
#include "constants.h"
#include "calib.h"
#include "histoBook.h"
#include <fstream>
#include <sstream>


// provides my own string shortcuts etc.
using namespace jdbUtils;

/*
*	Constructor
*
* 	Parameters:
*	chain:	The chain object containing all data compatible with the TOFrPicoDST format
*	nIterations:	Max number of iterations to run
*	xmlConfig:		The xml configuration defining key aspects of the calibration
*					such as number of tot bins to use, data location etc. See repo Readme
*					for a sample configuration.		
*/
calib::calib( TChain* chain, uint nIterations, xmlConfig con )  {
	cout << "[calib.calib] " << endl;
	
	config = con;

	// default number of tot bins is defined in constants in case non is given in config file
	numTOTBins = constants::numTOTBins;
	
	// set the number of tot bins from the config if given
	
	numTOTBins = config.getAsInt( "numTOTBins", constants::numTOTBins );
	minTOT = config.getAsDouble( "minTOT", constants::minTOT );
	maxTOT = config.getAsDouble( "maxTOT", constants::maxTOT );
	
	// now build arrays that need numTOTBins
	for ( int j = 0; j < constants::nChannels; j++){
		correction[ j ] 	= new double[ numTOTBins + 1 ];
		totBins[ j ] 		= new double[ numTOTBins + 1 ];
		useTOTBin[ j ] 		= new bool[ numTOTBins + 1 ];

		deadDetector[ j ]	= false;
	}

	// zero the corrections & offsets
	for ( int j = 0; j < constants::nChannels; j++){
		for (int k = 0; k < numTOTBins + 1; k++){
			correction[ j ] [ k ] = 0;
		}
		initialOffsets[ j ] = 0;
		spline[ j ] = NULL;
	}
	
	// set the maximum number of iterations
	maxIterations = nIterations;
	// create the histogram book
	book = new histoBook( ( config.getAsString( "baseName" ) + config.getAsString( "rootOutput" ) ) );

	_chain = chain;
	pico = new TOFrPicoDst( _chain );

	gStyle->SetOptStat( 0 );

	currentIteration = 0;

	gErrorIgnoreLevel=kError;
	// create a canvas for report building 
	report = new reporter( config.getAsString( "baseName" ) + config.getAsString( "reportOutput" ) );


	std::vector<double> tmp = config.getAsDoubleVector( "vzOutlierCut" );
	if ( tmp.size() >= 1 && tmp[ 0 ] != config.getAsString( "vzOutlierCut" ) )
		vzOutlierCut = tmp;
	else {
		vzOutlierCut.push_back( 40 );
		vzOutlierCut.push_back( 20 );
		vzOutlierCut.push_back( 15 );
		vzOutlierCut.push_back( 8 );
		vzOutlierCut.push_back( 5 );
	}
	
	tmp = config.getAsDoubleVector( "avgNTimingCut" );
	if ( tmp.size() >= 1 && tmp[ 0 ] != config.getAsString( "avgNTimingCut" ) )
		avgNTimingCut = tmp;
	else {
		avgNTimingCut.push_back( 2 );
		avgNTimingCut.push_back( 1 );
		avgNTimingCut.push_back( 0.6 );
	}

	avgNBackgroundCut = config.getAsDouble( "avgNBackgroundCut", 10 );

}

/*
*	Destructor
*	Deletes the histoBook ensuring it is saved.
*/
calib::~calib() {
	
	delete book;
	delete report;
	
	for ( int j = 0; j < constants::nChannels; j++){
		delete [] correction[j];
		delete [] totBins[j];
		if ( useTOTBin[ j ])
			delete [] useTOTBin[ j ];
		if ( spline [ j ] )
			delete spline[ j ];
	
	}
	cout << "[calib.~calib] " << endl;
}

/*
*	Offsets
*	Calculates the initial offsets for each channel with respect to channel 1 on the west side.
*	Then performs the offset calculation again after all corrections have been applied
*	
*/
void calib::offsets() {

	startTimer();

	if ( !_chain ){
		cout << "[calib." << __FUNCTION__ << "] ERROR: Invalid chain " << endl;
		return;
	}


	Int_t nevents = (Int_t)_chain->GetEntries();
	cout << "[calib." << __FUNCTION__ << "] Loaded: " << nevents << " events " << endl;

	book->cd( "initialOffset" );
	
	// make all the histos the first round
	book->make2D( "tdc", "tdc relative to channel 0", constants::nChannels, -0.5, constants::nChannels-0.5, 600, -100, 100 );
	book->make1D( "tdcMean", "tdc relative to channel 0", constants::nChannels, -0.5, constants::nChannels-0.5 );
	book->make2D( "correctedOffsets", "corrected Initial Offsets", constants::nChannels, -0.5, constants::nChannels-0.5, 200, -100, 100 );
	book->make2D( "tdcRaw", "All tdc Values ", constants::nChannels, 0, constants::nChannels, 1000, 0, 51200 );


	cout << "[calib." << __FUNCTION__ << "] Made Histograms " << endl;

	// loop over all events
	for(Int_t i=0; i<nevents; i++) {
    	_chain->GetEntry(i);
		
		progressBar( i, nevents, 75 );

		// perform the cuts used in calibration step to ensure the distributions match
		double tpcZ = pico->vertexZ;
		float vx = pico->vertexX;
    	float vy = pico->vertexY;
    	float vxy = TMath::Sqrt( vx*vx + vy*vy );
    	if ( vxy > 1 ) continue;
    	if ( pico->nTofHits <= 1 ) continue;
    	if ( TMath::Abs( tpcZ ) > 100 ) continue;


		// channel 1 on the west side is the reference channel
    	double reference = pico->vpdLeWest[0];


		for( int j = constants::startWest; j < constants::endEast; j++) {

			
			// skip dead detectors
			if ( deadDetector[ j ] ) continue;

			int nHits = pico->numHits( j );
			
			if ( nHits < constants::minHits ) 
				continue;

			double tdc = pico->channelTDC( j );
	    	double tot = pico->channelTOT( j );
	    
	    	book->fill( "tdcRaw", j, tdc );

	    	if(tot <= minTOT || tot >= maxTOT) continue;	    


		    book->fill( "tdc", j, tdc - reference );

		}	
	} // end loop on events

	// calculate the offsets
  	TH2D* tdc = (TH2D*) book->get( "tdc" );

	for ( int i = constants::startWest; i < constants::endEast; i++ ){
		TH1D* tmp = tdc->ProjectionY( "tmp", i+1, i+1 );
		cout << "Channel [ " << i+1 << " ] Offset = " << tmp->GetMean() << endl;

		this->initialOffsets[ i ] = tmp->GetMean();
		book->get( "tdcMean" )->SetBinContent( i+1, this->initialOffsets[ i ] );
		book->get( "tdcMean" )->SetBinError( i+1, tmp->GetMeanError() );

		delete tmp;
	}

	

	// get the east / west offset
	TH1D* west = tdc->ProjectionY( "westOffset", 2, 19 );
	book->add( "westOffset", west );
	TH1D* east = tdc->ProjectionY( "eastOffset", 20, 38 );
	book->add( "eastOffset", east );

	westMinusEast = ( west->GetMean() - east->GetMean() );
	
	cout << "West Mean = " << west->GetMean() << endl;
	cout << "East Mean = " << east->GetMean() << endl;
	cout << "West - East offset: " << westMinusEast << endl; 

	report->newPage( 1, 2);


	book->style( "tdc" )
		->set( "title", "Channel TDC wrt West Channel 1")
		->set( "range", -35.0, 35.0 )->draw();
	book->style( "tdcMean" )
		->set( "title", "Channel TDC wrt West Channel 1")
		->set( "markerStyle", 20)->set( "markerColor", 2 )
		->set( "linecolor", 2)->set( "draw", "same ple" )
		->draw();

	report->next();

	//book->style( "westOffset"+ts(currentIteration) )->set( "legend", "berberjh" );
	//book->clearLegend();

	book->style( "westOffset" )->set( "lineColor", kRed)
		->set( "title", "West (Channels 1-19) vs. East (Channels 20-38)" )->draw()->set( "legend", "West");

	book->style( "eastOffset" )->set( "lineColor", kBlue)
		->set( "title", "West (Channels 1-19) vs. East (Channels 20-38)" )
		->set( "draw", "same")
		->draw( )
		->set( "legend", "East")->set( "legend", legendAlignment::right, legendAlignment::top);
	
	report->savePage();

	cout << "[calib." << __FUNCTION__ << "] completed in " << elapsed() << " seconds " << endl;
}

/*
*	binTOT
*	variableBinning:
*		false	- fixed binning the tot space from minTOT to maxTOT
*		true	- calculates variable binning for tot space such that 
*				the number of events is roughly equal for each tot bin
*
*/
void calib::binTOT( bool variableBinning ) {

	cout << "[calib." << __FUNCTION__ << "] Starting " << endl;

	if ( variableBinning )
		cout << "[calib." << __FUNCTION__ << "] Variable Binning TOT Range :  " << minTOT << " -> " << maxTOT << endl;
	else
		cout << "[calib." << __FUNCTION__ << "] Fixed Binning TOT Range :  " << minTOT << " -> " << maxTOT << endl;

	cout << "[calib." << __FUNCTION__ << "] Using " << numTOTBins << " bins for TOT" << endl;

	startTimer();

	if ( variableBinning == false ){

		for(Int_t i=0; i<constants::nChannels; i++) {
			for ( int s = 0; s <= numTOTBins; s++ ){
				double edge = ((maxTOT - minTOT) / (double) numTOTBins) * s;
				edge += minTOT;
				totBins[ i ][ s ] = edge;
			}
		} // loop channels

	} else { // variableBinning == true ( default )

		Int_t nevents = (int)_chain->GetEntries();
		vector<double> tots[ constants::nChannels];

		cout << "[calib." << __FUNCTION__ << "] Processing " <<  nevents << " events" << endl;

		for(Int_t i=0; i<nevents; i++) {
	    	_chain->GetEntry(i);

			progressBar( i, nevents, 75 );
	    	Int_t numEast = pico->numberOfVpdEast;
	      	Int_t numWest = pico->numberOfVpdWest;
	     
		    if( numWest > constants::minHits){
		        
		    	for(Int_t j = 0; j < constants::endWest; j++) {
		        	Double_t tot = pico->channelTOT( j );
		          
		        	if(tot > minTOT && tot < maxTOT ) 
		          		tots[j].push_back(tot);
		        }

		    }

	  		if( numEast > constants::minHits ){
	    
	    		for(Int_t j = constants::startEast; j < constants::endEast; j++) {
	      			Double_t tot = pico->channelTOT( j );
	      
			        if( tot > minTOT && tot < maxTOT) 
			        	tots[j].push_back(tot);
	    		}

	  		}

		} // lopp events 	

		// get a threshold for a dead detector
		int threshold = 0;
		for(Int_t i=0; i<constants::nChannels; i++) {
			Int_t size = tots[i].size();
			threshold += size;
		}
		threshold /= (double)constants::nChannels; // the average of all detectors
		threshold *= .25;

		// loop through the channels and determine binning
		for(Int_t i=0; i<constants::nChannels; i++) {
	      
	    	Int_t size = tots[i].size();
	      	cout << "[calib.binTOT] Channel[ " << i << " ] : " << size << " hits" << endl;
	      	
	      	if( size < threshold ) { // check for dead channels
	        	
	        	Double_t step = ( maxTOT - minTOT ) / numTOTBins;

	        	for(Int_t j=0; j <= numTOTBins; j++) {

	                totBins[ i ][ j ] = ( step * j ) + minTOT; 
	        	}
	        	cout  << "[calib.binTOT] VPD Channel [ " << i << " ] is dead! " << "( " << size << " hits)" <<endl;
	        	
	        	// set this detector to dead
	        	deadDetector[ i ] = true;

	      	} else { // channel not dead

	      		deadDetector[ i ] = false;
	      		Int_t step = size / (numTOTBins + 1 ); 
	    
	      		// sort into ascending order
	        	std::sort( tots[i].begin(), tots[i].end());
	        	
	        	totBins[ i ][0] = tots[ i ].at(0);
	        	totBins[ i ][ numTOTBins ] = maxTOT;
	        	
	        	for( Int_t j = 1; j < numTOTBins ; j++) {

	        		double d1 = tots[i].at( step * j );
	        		double d2 = tots[ i ].at( step * (j - 1) );
	        	
	        		totBins[ i ][ j ] = ( d1 + d2 ) / 2.0;
	            	
	        	}	// loop over tot bins

	      } // end channle not dead

	  	} // end loop channles
	  	
		for(Int_t i = 0; i< constants::nChannels; i++) {
			tots[i].clear();
		}	
	}
	cout << "[calib." << __FUNCTION__ << "] completed in " << elapsed() << " seconds " << endl;

}

double calib::getCorrection( int vpdChannel, double tot ){
	
	//if ( currentIteration <= 0 )
	//	return 0;

	if ( spline[ vpdChannel ] && spline[ vpdChannel ]->getSpline() ){
		
		return spline[ vpdChannel ]->getSpline()->Eval( tot );
	}
	
	int totBin = binForTOT( vpdChannel, tot ); 
	
	return correction[ vpdChannel ][ totBin ];

}

int calib::binForTOT( int vpdChannel, double tot ){

	stringstream sstr; 
	sstr << "channel" << vpdChannel;
	string old = book->cd( sstr.str() );
	sstr.str("");    	sstr << "it" << (currentIteration - 1) <<  "totcor";	
	TH1D* tmp = (TH1D*)book->get( sstr.str() );
	int bin = 0;
	if ( tmp ){
		bin = tmp->GetXaxis()->FindBin( tot );
	} else {
		//cout << "[calib." << __FUNCTION__ << "] Cant Find Tot Bin for tot = " << tot << " in channel : " << vpdChannel << endl;
	}
	book->cd( old );

	return bin;

}


void calib::outlierRejection( bool reject ) {

	// must be called from inside event loop in the calib step
	string iStr = "it"+ts(currentIteration);

	if ( reject == false ){
		// reset the state
		for ( int j = constants::startWest; j < constants::endEast; j++ ){
			useDetector[ j ] = true;
		}
		westIsGood = true;
		eastIsGood = true;
		return;
	}

	// get the TPC z vertex
	double tpcZ = pico->vertexZ;

	double vzCut = 40;
	if ( currentIteration < vzOutlierCut.size() )
		vzCut = vzOutlierCut[ currentIteration ];	// use the cut for this step
	else 
		vzCut = vzOutlierCut[ vzOutlierCut.size() - 1 ];	// after that use the last cut defined for all other steps

	book->cd( "OutlierRejection" );

	int numValidPairs = 0;

	// reset the state
	for ( int j = constants::startWest; j < constants::endEast; j++ ){
		useDetector[ j ] = false;
	}

	eastIsGood = false;
	westIsGood = false;

	double sumEast = 0;
	double sumWest = 0;
	double countEast = 0;
	double countWest = 0;
	stringstream sstr;

	for ( int j = constants::startWest; j < constants::endWest; j++ ){

		if ( deadDetector[ j ] ) continue;

		double tdcWest = pico->channelTDC( j );
	    double totWest = pico->channelTOT( j );

	    tdcWest -= this->initialOffsets[ j ];

	  
	    if( totWest <= minTOT || totWest > maxTOT) continue;
	    
	    double corWest = getCorrection( j, totWest );
	    tdcWest -= corWest;

	    sumWest += tdcWest;
	    countWest++;

		for ( int k = constants::startEast; k < constants::endEast; k++ ){
			
			if ( deadDetector[ k ] ) continue;

			double tdcEast = pico->channelTDC( k );
	    	double totEast = pico->channelTOT( k );

	    	tdcEast -= this->initialOffsets[ k ];

	    	if( totEast <= minTOT || totEast > maxTOT) continue;
	    	
	    	double corEast = getCorrection( k, totEast );
	    	tdcEast -= corEast;

	    	if ( j == constants::startWest ){
	    		sumEast += tdcEast;
	    		countEast++;
	    	}

	    	// calculate the VPD z Vertex
	    	double vpdZ = constants::c * ( tdcEast - tdcWest) / 2.0;

	    	book->get( iStr+"All" )->Fill( tpcZ - vpdZ );
	    	book->get( iStr +"zTPCzVPD" )->Fill( tpcZ, vpdZ );
	    	

	    	if ( TMath::Abs( tpcZ - vpdZ ) < vzCut  ){

	    		// valid pair
	    		useDetector[ k ] = true;
	    		useDetector[ j ] = true;
	    		eastIsGood = true;
	    		westIsGood = true;
	    		numValidPairs ++ ;

	    	} 
		    	  			
		} // loop channel k
	} // loop channel j

	if ( countEast >= 1 && countWest >= 1){
		double vpdZ = constants::c * ( (sumEast/countEast) - (sumWest/countWest)) / 2.0;	
		book->fill( iStr+"zTPCzVPDAvg", tpcZ, vpdZ );
		book->fill( iStr+"avg", ( tpcZ-vpdZ ));
	}

	book->fill( iStr+"nValidPairs", numValidPairs );

	int nAccepted = 0;
	for ( int j = constants::startWest; j < constants::endWest; j++ ){
		if( useDetector[ j ] )
			nAccepted ++;
	}


	book->fill( iStr+"nAcceptedWest", nAccepted );

	nAccepted = 0;
	for ( int j = constants::startEast; j < constants::endEast; j++ ){
		if( useDetector[ j ] )
			nAccepted ++;
	}

	book->fill( iStr+"nAcceptedEast", nAccepted );

}

void calib::prepareStepHistograms() {
	
	// for names
	string iStr = "it"+ts(currentIteration);
	// for titles
	string step = "Step " + ts( currentIteration+1 ) + " : ";

	/*
	* check that our histos are made for this iteration
	*/
	if ( book->get( "correctedOffsets", "initialOffset"  ) == 0){
		book->cd( "initialOffset" );
		
	}
	
	for ( int ch = constants::startWest; ch < constants::endEast; ch++ ){
		
		book->cd( "channel" + ts(ch) );

		// make channel titles start at 1
		string sCh = "Channel "+ts(ch+1);

		book->make2D( 	iStr + "tdctot", step +sCh+" TDC vs TOT ;tot [ns];tdc [ns]", 
							numTOTBins , totBins[ ch ], 1000, -40, 40 );
		book->make2D( 	iStr + "tdccor", step +sCh+" TDC vs TOT ;tot [ns];tdc [ns]", 
							numTOTBins , totBins[ ch ], 1000, -20, 20 );
		book->make1D( 	iStr + "tdc", step + sCh+" TDC;tot [ns];tdc [ns]", 
							500, -10, 10 );
		book->make2D( 	iStr + "avgN", step + sCh + " : 1 - <N>;# of Detectors;tdc [ns]", 
						constants::nChannels/2, 1, constants::nChannels/2, 1000, -20, 20 );
		book->make2D( 	iStr + "cutAvgN", step + sCh + " : 1 - <N>;# of Detectors;tdc [ns]", 
							constants::nChannels/2, 1, constants::nChannels/2, 1000, -20, 20 );
	}

	/*
	* outlier rejection histos
	*/
	book->cd( "OutlierRejection" );

	int zBins = 600, zRange = 200;

	book->make1D( 	iStr + "All", step + "Outlier Rejection; z_{TPC} - z_{VPD}; [#]", zBins, -zRange, zRange );
	book->make1D( 	iStr + "avg", step + "TPC vs. VPD z Vertex using <East> & <West>; z_{TPC} - z_{VPD} [cm]; [#]", 	zBins, -zRange, zRange );
	
	book->make2D( 	iStr + "zTPCzVPD", step + "TPC vs. VPD z Vertex; z_{TPC};z_{VPD}", zBins/2, -zRange/2, zRange/2, zBins/2, -zRange/2, zRange/2 );
	book->make2D( 	iStr + "zTPCzVPDAvg", step + "TPC vs. VPD z Vertex using <East> & <West>; z_{TPC} [cm];z_{VPD} [cm]", zBins/2, -zRange/2, zRange/2, zBins/2, -zRange/2, zRange/2 );

	book->make1D( 	iStr + "nValidPairs", step + "# of Valid Pairs; # of Pairs; [#]", 500, 0, 500 );
	book->make1D( 	iStr + "nAcceptedWest", step + "# of Accepted Detectors; # of Detectors; [#] ",
							19, -0.5, 18.5 );
	book->make1D( 	iStr + "nAcceptedEast", step + "# of Accepted Detectors; # of Detectors; [#] ",
							19, -0.5, 18.5 );
	/*
	* outlier rejection histos
	*/

	// offsets
	book->cd( "initialOffset" );
	book->make2D( 	iStr + "Offsets", step + " tdc wrt Channel 1; Detector ; [#] ",
							constants::nChannels, -0.5, constants::nChannels-0.5, 600, -100, 100 );

	cout << "[calib." << __FUNCTION__ << "[" << currentIteration << "]] " << " Histograms Booked " << endl;

}

void calib::step( ) {

	cout << "[calib." << __FUNCTION__ << "[" << currentIteration << "]] " << " Start " << endl;
	
	startTimer();

	bool outliers =  config.getAsBool( "outlierRejection" );
	bool removeOffset = config.getAsBool( "removeOffset" );

	double outlierCut = 2;
	if ( currentIteration < avgNTimingCut.size() )
		outlierCut = avgNTimingCut[ currentIteration ];	// use the cut for this step
	else 
		outlierCut = avgNTimingCut[ avgNTimingCut.size() - 1 ];	// after that use the last cut defined for all other steps

	
	// the data we will use over and over 
	double tot[ constants::nChannels ];		// tot value
	double tdc[ constants::nChannels ];		// tdc value
	double off[ constants::nChannels ];		// offset value
	double tAll[ constants::nChannels ];	// tdc with all corrections ( offset and correction)

	// correction based on channel and tot value
	double corr[ constants::nChannels ];	
	// reference tdc time => the 1st channel on the west side
	double reference;

	string iStr = "it"+ts(currentIteration);
	stringstream sstr;

	// make sure the histograms are ready
	prepareStepHistograms();
	

	Int_t nevents = (int)_chain->GetEntries();
	for(Int_t i = 0; i < nevents; i++) {
    	_chain->GetEntry(i);

		progressBar( i, nevents, 75 );

    
    	float vx = pico->vertexX;
    	float vy = pico->vertexY;
    	float vxy = TMath::Sqrt( vx*vx + vy*vy );
    	if ( vxy > 1 ) continue;

    	double tpcZ = pico->vertexZ;
    	if ( pico->nTofHits <= 1 ) continue;
    	if ( TMath::Abs( tpcZ ) > 100 ) continue;

    	// perform outlier rejection for this event
    	outlierRejection( outliers );

    	if ( removeOffset )
  		  	averageN();
 		   	

    	// Alias the values for this event for ease
    	for( int j = constants::startWest; j < constants::endEast; j++) {
    		
    		if ( deadDetector[ j ] ) continue;
			if ( !useDetector[ j ] ) continue;

    		tot[ j ] = pico->channelTOT( j );
    		tdc[ j ] = pico->channelTDC( j );
    		if ( removeOffset )
   	 			off[ j ] = this->initialOffsets[ j ];
   	 		else 
   	 			off[ j ] = 0;
   			tAll[ j ] = tdc[ j ] - off[ j ];

    		if(tot[ j ] <= minTOT || tot[ j ] >= maxTOT) continue;
  			
  			corr[ j ] = getCorrection( j, tot[ j ] );
  			tAll[ j ] -= corr[ j ];
    	}
    	reference = pico->vpdLeWest[0];

		// loop over every channel on the west and then on the east side
		for( int j = constants::startWest; j < constants::endEast; j++) {
			
			// skip dead detectors
			if ( deadDetector[ j ] ) continue;
			if ( !useDetector[ j ] ) continue;
			// require the tot is within range
	    	if(tot[ j ] <= minTOT || tot[ j ] > maxTOT) continue;


	    	double tdcSumWest = 0;
			double tdcSumEast = 0;
	    	double countEast = 0;
	    	double countWest = 0;

	    	for( int k = constants::startWest; k < constants::endEast; k++) {

	    		// skip dead detectors
				if ( deadDetector[ k ] ) continue;
				if ( !useDetector[ k ] ) continue;
	    		
	    		if(tot[ k ] <= minTOT || tot[ k ] > maxTOT) continue;
	    		if ( j == k ) continue;
	    		
	    		if ( k >= constants::startWest && k < constants::endWest ){
	    			tdcSumWest += ( tAll[ k ]);
	    			countWest ++;
	    		} else if ( k >= constants::startEast && k < constants::endEast ){
	    			tdcSumEast += ( tAll[ k ]);
	    			countEast ++;
	    		}

			}	// loop on vpdChannel k


			/*
			*	Now recalculate the average times using the previously calculated average to
			*	apply a cut on the range of variation
			*/
			double cutSumWest = 0;
			double cutSumEast = 0;
	    	double cutCountEast = 0;
	    	double cutCountWest = 0;

	    	for( int k = constants::startWest; k < constants::endEast; k++) {

	    		// skip dead detectors
				if ( deadDetector[ k ] ) continue;
				if ( !useDetector[ k ] ) continue;
	    		if ( j == k ) continue;
	    		if(tot[ k ] <= minTOT || tot[ k ] > maxTOT) continue;

	    		if ( k >= constants::startWest && k < constants::endWest ){
	    			double tAvg = tdcSumWest / countWest;

	    			if ( tAll[ k ] - tAvg < outlierCut && tAll[ k ] - tAvg > -outlierCut ){
	    				cutSumWest += tAll[ k ];
	    				cutCountWest ++;
	    			}
	    		} else if ( k >= constants::startEast && k < constants::endEast ){
	    			double tAvg = tdcSumEast / countEast;

	    			if ( tAll[ k ] - tAvg < outlierCut && tAll[ k ] - tAvg > -outlierCut ){
	    				cutSumEast += tAll[ k ];
	    				cutCountEast ++;
	    			}
	    		}

			}	// loop on vpdChannel k


			book->cd( "initialOffset" );
	 		if ( currentIteration == 0 ){
	 			//Plot the offsets after correction just to be sure it all works
		    	book->fill( "correctedOffsets", j, tdc[ j ] - reference - off[ j ] );
		    }
		    // now fill the offsets to see how it changes with the cuts / outlier rejection
		    book->fill( iStr+"Offsets", j, tdc[ j ] - corr[ j ] - (reference - corr[ 0 ]));

		    // set the avg and count varaibles for this run
		    // if j corresponds to a west channel then use tdcSumWest, countWest
		    // if j corresponds to an east channel then use tdcSumEast, countEast
	    	double avg = (tdcSumWest / countWest );
	    	double cutAvg = ( cutSumWest / cutCountWest );
	    	int count = countWest;
	    	int cutCount = cutCountWest;

	    	if ( !removeOffset ){
	    		cutAvg = avg;
	    		cutCount = count;
	    	}


	 		if ( j >= constants::startEast && j < constants::endEast ){
	    		avg = (tdcSumEast / countEast );
	    		count = countEast;
	    		
	    		if ( removeOffset ){
	    			cutAvg = ( cutSumEast / cutCountEast );
	    			cutCount = cutCountEast;
	    		} else {
	    			cutAvg = avg;
	    			cutCount = count;
	    		}

	    	}

	    	if ( count <= constants::minHits ) continue;

	    	// change into this channels dir for histogram saving
			book->cd( "channel" + ts(j) );	    	
	    	book->fill( iStr+"tdctot", tot[ j ], tdc[ j ] - off[ j ] - cutAvg );
	    	book->fill( iStr+"tdccor", tot[ j ], tAll[ j ] - cutAvg );
	    	book->fill( iStr+"tdc" , tAll[ j ]  - cutAvg );
	
		}	
	}

	cout << "[calib." << __FUNCTION__ << "[" << currentIteration << "]] " << " completed in " << elapsed() << " seconds " << endl;
	
	makeCorrections();
	stepReport();

	currentIteration++;

	
}

/*
*
*	finish
*	called *after* the last iteration to perform final cuts and finish the calibration calculations
*
*/
void calib::finish( ){

	int last = currentIteration - 1;
	string iStr = "it" + ts( last );

	TF1 * g = new TF1( "g", "gaus", -1.0, 1.0 );

	report->newPage( 3, 2);
	
	book->cd ( "final" );

	TH1D*	sigmas = new TH1D( 	"detSigma", "Detector Resolution;Detector; Resolution [ ns ]",
								constants::nChannels, 0.5, constants::nChannels+0.5 ); // add 0.5 to center the bin numbers

	TLatex		*text = new TLatex();
	text->SetNDC();
	text->SetTextSize(0.044);

	for ( int j = constants::startWest; j < constants::endEast; j++ ){

		string iCh = "channel" + ts( j );
		
		// move to the next pad in the report PDF
		if ( j > constants::startWest )
			report->next();

		book->cd( "final/fit" );

		TH2D* tmp = (TH2D*)book->get( iStr + "cutAvgN", iCh );
		tmp->FitSlicesY( g, 0, -1, avgNBackgroundCut );
		TH1D* fsySig = (TH1D*)gDirectory->FindObject( (iStr + "cutAvgN" + "_2").c_str() );
		TH1D* fsyMean = (TH1D*)gDirectory->FindObject( (iStr + "cutAvgN" + "_1").c_str() );
		
		
		book->cd ( "final" );

		TH1D* sigFit = (TH1D*)fsySig->Clone( (iCh + "sigmaFit").c_str() );
		book->add( iCh + "sigmaFit", sigFit );
		TH1D* mean = (TH1D*)fsyMean->Clone( (iCh + "sigmaMean").c_str() );
		book->add( iCh + "sigmaMean", mean );

		TF1 * fr = new TF1( "fr", calib::detectorResolution, 0, 19, 1);

		sigFit->Fit( "fr", "QR" );

		cout << "Channel [ " << j << " ] Resolution: " << fr->GetParameter( 0 ) << " ns " << endl;
		sigmas->SetBinContent( j + 1, fr->GetParameter( 0 ) );
		sigmas->SetBinError( j + 1, fr->GetParError( 0 ) );

		double max = book->get( iCh + "sigmaFit" )->GetMaximum();

		book->style( iCh + "sigmaFit" )	->set( "title", "Detector # " + ts(j) + " Resolution Fit" )
										->set( "y", "Time [ns]")
										->set( "markerStyle", 17)
										->set( "markerColor", 2)
										->set( "range", -0.2, max + .1 )
										->draw();
		book->style( iCh + "sigmaMean" )->set( "markerStyle", kCircle )->set( "draw", "same")
										->draw( );
		double chi		 	= fr->GetChisquare();
		float np			= fr->GetNumberFitPoints();
		double chiDoF = 0;
		if (np>1){
			chiDoF	= chi/((Float_t)(np-1));
		}

		text->DrawLatex(0.25,0.86, ("#sigma = " + ts( fr->GetParameter( 0 ) )).c_str() );
		text->DrawLatex(0.25,0.81, ("#chi^{2}/DOF = " + ts( chiDoF) ).c_str() );

		
	    if ( j == constants::endEast - 1){
	    	report->savePage();
	    }

	}

	report->newPage();

	double max = sigmas->GetMaximum();


	book->cd ( "final" );
	book->add( "detSigma", sigmas );
	book->get( "detSigma" )->Fit( "pol0", "Q" );
	book->style( "detSigma" )
		->set( "markerStyle", 17)
		->set( "markerColor", 2 )
		->set( "range", 0, max + .015 )
		->draw();

	report->savePage();	


	// get the final channel offsets
	book->cd( "initialOffset" );

	// calculate the offsets
  	TH2D* tdc = (TH2D*) book->get( iStr + "Offsets" );
	for ( int i = constants::startWest; i < constants::endEast; i++ ){
		TH1D* tmp = tdc->ProjectionY( "tmp", i+1, i+1 );
		delete tmp;
	}

	// get the east / west offset
	TH1D* west = tdc->ProjectionY( "westOffsetFinal", 2, 19 );
	book->add( "westOffsetFinal", west );
	TH1D* westAll = tdc->ProjectionY( "westOffsetAllFinal", 1, 19 );
	book->add( "westOffsetAllFinal", westAll );

	finalWestOffset = ( west->GetMean() + westAll->GetMean() ) / 2.0;
	cout << "Final West Offset : " << finalWestOffset << endl;


}

void calib::loop( ) {

	for ( unsigned int i = 0; i < maxIterations; i++ ){
		step();
	}
	finish();

}

void calib::makeCorrections( ){

	cout << "[calib." << __FUNCTION__ << "[" << currentIteration << "]] " << " Start " << endl;
	

	string iStr = "it" + ts( currentIteration );
	startTimer();


	/**
	 * First time generate the list of good / bad bins if using fixed binning
	 */
	if ( 0 == currentIteration && false == config.getAsBool( "variableBinning" ) ){
		

		for( int k = constants::startWest; k < constants::endEast; k++) {
			book->cd( "channel" + ts( k ) );

			// slewing curve without correction applied to channel k
	    	TH2D* pre = (TH2D*) book->get( iStr + "tdctot" );
	    	double total = pre->Integral();

	    	// if there are less than N% the number of events that would be in each bin 
	    	// assuming equal distribution then reject the bin
	    	double minEvents = total / ( numTOTBins * 5.0);
	    	for ( int ib = 1; ib <= numTOTBins ; ib ++ ){
	    		useTOTBin[ k ][ ib ] = true;
	    		double binTotal = pre->Integral( ib, ib );
	    		if ( binTotal <= minEvents ){
	    			useTOTBin[ k ][ ib ] = false;
	    			//cout << " Channel [ " << k << " ] Bin [ " << ib << " ] set to false " << endl;
	    		}
	    	}
		}
	}
	


	report->newPage( 4, 5);

	// get the corrections for the next iteration 
	for( int k = constants::startWest; k < constants::endEast; k++) {

		if ( currentIteration <= 1 || currentIteration == maxIterations - 1){
			if ( deadDetector[ k ] ){
				report->next();
				if ( k == constants::endWest - 1){
			    	report->savePage();
			    	report->newPage( 4, 5);
		    	}
		    	continue;
			}
			if ( k != constants::startWest && k != constants::startEast )
				report->next();
		}

		// switch into channel dir
		book->cd( "channel" + ts( k ) );

		// slewing curve without correction applied to channel k
	    TH2D* pre = (TH2D*) book->get( iStr + "tdctot" );

	    // slewing curve with correction applied to channel k
	    TH2D* post = (TH2D*) book->get( iStr + "tdccor" );

		book->cd( "channel" + ts( k ) + "/fit" );

		// do the fit
	    pre->FitSlicesY();

	    TH1D* preMean = (TH1D*) gDirectory->FindObject( 
	    					(iStr + "tdctot_1").c_str() );

	    // do the fit
	    TF1* g = new TF1( "g", "gaus", -1, 1 );
	    post->FitSlicesY( g );
	    delete g;

	    TH1D* postMean = (TH1D*) gDirectory->FindObject( 
	    					(iStr + "tdccor_1").c_str() );

	    book->cd( "channel" + ts( k ) );

	    TH1D* cor = (TH1D*) preMean->Clone( (iStr + "totcor").c_str() );
	    book->add( (iStr + "totcor").c_str(), cor  );

	    TH1D* dif = (TH1D*) postMean->Clone( (iStr + "difcor").c_str() );
	    book->add( ("it" + ts( currentIteration ) + "difcor").c_str(), dif  );

	    double reset = 0;
	    for ( int ib = 1; ib <= numTOTBins ; ib ++ ){

	    	// reject outlier bins
	    	if ( ib >= 2 && ib <= numTOTBins - 1){

	    		double s1 = cor->GetBinContent( ib - 1 );
	    		double s2 = cor->GetBinContent( ib + 1 );
	    		double avgSides = (s1 + s2 ) / 2.0;
	    		double val = cor->GetBinContent( ib );

	    		if ( 	TMath::Abs( val - avgSides ) >= 1 && TMath::Abs( val - s1 ) >= .75 && TMath::Abs( val - s2 ) >= .75 ) /////// TODO
	    			cor->SetBinContent( ib, avgSides );

	    	} else if ( ib == 1 ){	// first bin
	    		double val = cor->GetBinContent( ib );
	    		double s2 = cor->GetBinContent( ib + 1 );
	    		if ( TMath::Abs( val - s2 ) >= 1 )
	    			cor->SetBinContent( ib, s2 );
	    	} else if ( ib == numTOTBins ){	// last bin
	    		double val = cor->GetBinContent( ib );
	    		double s1 = cor->GetBinContent( ib - 1 );
	    		if ( TMath::Abs( val - s1 ) >= 1 )
	    			cor->SetBinContent( ib, s1 );
	    	}
	    	
	    	if ( useTOTBin[ k ][ ib ]  )
	    		reset = cor->GetBinContent( ib );
	    	else{
	    		cor->SetBinContent( ib, reset );
	    		dif->SetBinContent( ib, 0 );
	    	}

	    	correction[ k ][ ib  ] = cor->GetBinContent( ib );
	    	
	    }


	    if ( spline[ k ])
	    	delete spline[ k ];
	   
	   	// set the spline type
	   	// default to none
	   	Interpolation::Type type = ROOT::Math::Interpolation::kAKIMA;
	   	bool useSpline = false;
	   	if ( "akima" == config.getAsString( "splineType", "akima" ) ){
	    	type = ROOT::Math::Interpolation::kAKIMA;	
	    	useSpline = true; 
	    } else if ( "linear" == config.getAsString( "splineType" ) ){
	    	type = ROOT::Math::Interpolation::kLINEAR;	
	    	useSpline = true;
	    } else if ( "cspline" == config.getAsString( "splineType" ) ){
	    	type =  ROOT::Math::Interpolation::kCSPLINE;	
	    	useSpline = true;
	    }
	    
	    if ( useSpline )
		    spline[ k ] = new splineMaker( cor, splineAlignment::center, type );

		splineMaker* vSpline;
		if ( currentIteration == 0 ) 
			vSpline = spline[ k ];
		else
			vSpline = new splineMaker( dif, splineAlignment::center, type );
	    
	    if ( currentIteration <= 1 || currentIteration == maxIterations - 1){
		    

		    book->style( ("it"+ts(currentIteration)+"tdccor") )
		    	->set( "range", -5.0, 5.0);


		    post->Draw( "colz" );
		    
		    if ( useSpline ){
		    	TGraph* g = vSpline->graph( minTOT, maxTOT, 0.2);
		    	g->GetYaxis()->SetRangeUser( -5, 5);
		    	g->SetMarkerStyle(7);
		    	g->SetMarkerColor( kRed );
		    	g->Draw( "same cp" );
			}
			post->GetYaxis()->SetRangeUser( -5, 5);
			if ( currentIteration != 0 ) 
				delete vSpline;

		
		    if ( k == constants::endEast - 1 || k == constants::endWest - 1){
		    	report->savePage();
		    	report->newPage( 4, 5);
	    	}
    	}
	    

	}

	cout << "[calib." << __FUNCTION__ << "[" << currentIteration << "]] " << " completed in " << elapsed() << " seconds " << endl;
}

void calib::writeParameters(  ){

	cout << "[calib." << __FUNCTION__ << "[" << currentIteration << "]] " << " Start " << endl;

	string outName = config.getAsString( "baseName" ) + config.getAsString( "paramsOutput" );

	if ( outName.length() <= 4 )
		return;

	ofstream f;
	f.open( outName.c_str() );

	bool removeOffset = config.getAsBool( "removeOffset" );

	book->cd( "correctionParameters" );

	for ( int j = constants::startWest; j < constants::endEast; j++ ){

		book->make1D( "slewingCor" +ts(j), "Channel "+ts(j+1)+ "Slewing Correction", numTOTBins, minTOT, maxTOT  );

		f << (j + 1) << endl;
		f << numTOTBins << endl;

		for ( int i = 0; i <= numTOTBins; i++ ){
			f << totBins[ j ][ i ] << " ";
		}
		f << endl;
		for ( int i = 0; i <= numTOTBins; i++ ){
			double off = 0;
			if ( removeOffset ){
				off = initialOffsets[ j ] - finalWestOffset;
				if ( j >= constants::startEast && j < constants::endEast )
					off += westMinusEast;
			}
			double fCor = correction[ j ][ i + 1 ] + off;
			book->get( "slewingCor" +ts(j) )->SetBinContent( i+1, fCor );
			f << ( fCor ) << " ";
		}
		f << endl;
	} 

	f.close();

	outName = config.getAsString( "baseName" ) + "_spline_" + config.getAsString( "paramsOutput" );

	if ( outName.length() <= 4 )
		return;

	f.open( outName.c_str() );

	

	for ( int j = constants::startWest; j < constants::endEast; j++ ){

		book->make1D( "splineSlewingCor" +ts(j), "Channel "+ts(j+1)+ " : Slewing Correction", numTOTBins, minTOT, maxTOT );

		f << (j + 1) << endl;
		f << numTOTBins << endl;

		for ( int i = 0; i <= numTOTBins; i++ ){
			f << totBins[ j ][ i ] << " ";
		}
		f << endl;
		for ( double tot = minTOT; tot <= maxTOT; tot += (maxTOT - minTOT)/(double)numTOTBins ){
			double off = 0;
			if ( removeOffset ){
				off = initialOffsets[ j ] - finalWestOffset;
				if ( j >= constants::startEast && j < constants::endEast )
					off += westMinusEast;
			}
			double cor = spline[ j ]->getSpline()->Eval( tot );
			double fCor = cor + off;
			book->get( "splineSlewingCor" +ts(j) )->SetBinContent( book->get( "splineSlewingCor" +ts(j) )->GetXaxis()->FindBin( tot ), 
				cor );

			f << ( fCor ) << " ";
		}
		f << endl;
	} 

	f.close();

	//draw the parameters
	report->newPage( 3, 4 );

	for ( int j = constants::startWest; j < constants::endEast; j++ ){


		if ( j != constants::startEast && j != constants::startWest )
			report->next();

		book->style( "splineSlewingCor" +ts(j) )->set("lineColor", kRed )->set( "range", -5, 5)->draw();

		if ( j == constants::endWest - 1 || j == constants::endEast - 1){
			report->savePage();
			report->newPage( 3, 4 );
		}

	}
	cout << "[calib." << __FUNCTION__ << "] " << " completed in " << elapsed() << " seconds " << endl;

}

void calib::readParameters(  ){


	

	cout << "[calib." << __FUNCTION__ << "] " << " Start " << endl;
	
	startTimer();

	
	
	if ( config.isVector( "paramsInput") ){


		

		bool good = true;
		std::vector<string> files = config.getAsStringVector( "paramsInput" );
		std::vector<string> lNames = config.getAsStringVector( "paramsInput" );
		if ( config.isVector( "paramsLegend" ) ){
			lNames = config.getAsStringVector( "paramsLegend" );
		}

		for ( uint fi = 0; fi < files.size(); fi ++ ){
			string inName =  files[ fi ];

			ifstream infile;
	    	infile.open( inName.c_str() );
	    	cout << "[calib." << __FUNCTION__ << "] " << " Reading slewing corrections from " << inName << endl;
			if ( infile.is_open() ){
				for ( int i = 0 ; i < constants::nChannels; i++ ){

					int channel = -1;
					infile >> channel;

					int tBins = 0;
					infile >> tBins;
				
					if ( 	channel >= 1 && channel <= 38 && 
							tBins == numTOTBins ){

						double tmp = 0;
						for ( int j=0; j <= tBins; j++ ){
							infile >> tmp;
							totBins[ channel - 1 ][ j ] = tmp;
							
						}
						
						// create a hist
						book->cd( "params" );
						book->make1D( "file"+ts(fi)+"channel"+ts(channel-1), "Channel "+ts(channel)+" Slewing Corrections", tBins, totBins[ channel-1] );


						for ( int j=0; j <= tBins; j++ ){
							infile >> tmp;
							correction[ channel - 1 ][ j ] = tmp;
							book->get( "file"+ts(fi)+"channel"+ts(channel-1) )->SetBinContent( j+1, tmp );
						}
						
					} else {
						cout << "[calib." << __FUNCTION__ <<  "] " << "Bad file format, cannot load parameters" << endl;
						good = false;

					}


				}
			} else {
				cout << "[calib." << __FUNCTION__ <<  "] " << " Cannot open file " << endl;  
				good = false;
			}

			infile.close();
		}

		if ( good ){
			cout << "[calib." << __FUNCTION__ <<  "] " << " Read parameters for all channels " << endl;

		    report->newPage( 2, 3);
			for ( int k = 0; k < constants::nChannels; k++ ){

				for ( uint fi = 0; fi < files.size(); fi ++ ){
				
					
					
				    book->style( "file"+ts(fi)+"channel"+ts(k) )
				    	->set( "lineColor", 1+fi );
				    if ( 0 == fi ){
				    	book->clearLegend();
					    book->style( "file"+ts(fi)+"channel"+ts(k) )
							->draw()
							->set( "legend", lNames[ fi ] )->set( "legend", legendAlignment::right, legendAlignment::top, .4);
					} else {
						book->style( "file"+ts(fi)+"channel"+ts(k) )
							->set( "draw", "same" )
							->draw()
							->set( "legend", lNames[ fi ]);
					}

				}

				report->next();

			    if ( k == constants::nChannels-1)
			    	report->savePage();
			}
			currentIteration = 1;
		}
	}

	cout << "[calib." << __FUNCTION__ <<  "] " << " completed in " << elapsed() << " seconds " << endl;
}

void calib::stepReport() {

	string iStr = "it"+ts(currentIteration);


	double vzCut = 40;
	if ( currentIteration < vzOutlierCut.size() )
		vzCut = vzOutlierCut[ currentIteration ];	// use the cut for this step
	else 
		vzCut = vzOutlierCut[ vzOutlierCut.size() - 1 ];	// after that use the last cut defined for all other steps


	// build filled area histo
	book->cd( "OutlierRejection" );
	TH1D* fill = (TH1D*)book->get( iStr+"All" )->Clone( "fill" );
	fill->SetFillColor( kGreen );
	fill->GetXaxis()->SetRangeUser( -vzCut, vzCut );

	// report the outlier rejection
	report->newPage(2, 2);
	
	gPad->SetLogy(1);
	book->style( iStr+"All" )->draw();
	fill->Draw("same" );

	report->cd( 2, 1 );
	gPad->SetLogy(1);
	book->style( iStr+"avg" )->set( "draw", "" )->draw();

	report->cd( 1, 2 );
	book->style( iStr+"zTPCzVPD" )->set( "draw", "colz" )->draw();

	report->cd( 2, 2 );
	book->style( iStr+"zTPCzVPDAvg" )->set( "draw", "colz" )->draw();
	report->savePage();

	
	book->clearLegend();
	
	report->newPage(1, 2);

	gPad->SetLogy(1);
	book->style( iStr+"nAcceptedWest" )->draw( "draw", "")
		->draw()->set( "legend", "West" )
		->set("legend", legendAlignment::right, legendAlignment::top);
	
	book->style( iStr+"nAcceptedEast" )
		->set( "lineColor", kRed )->set( "draw", "same" )->draw()->set( "legend", "East");

	report->next();
	gPad->SetLogy(1);
	book->style( iStr+"nValidPairs" )->set( "dynamicDomain", 0.0f, 1, -1 )->draw();

	report->savePage();
}

void calib::averageN() {

	string iStr = "it"+ts(currentIteration);
	double outlierCut = 2;

	if ( currentIteration < avgNTimingCut.size() )
		outlierCut = avgNTimingCut[ currentIteration ];	// use the cut for this step
	else 
		outlierCut = avgNTimingCut[ avgNTimingCut.size() - 1 ];	// after that use the last cut defined for all other steps
	
	// the data we will use over and over 
	double tot[ constants::nChannels ];		// tot value
	double tdc[ constants::nChannels ];		// tdc value
	double off[ constants::nChannels ];		// offset value
	double tAll[ constants::nChannels ];	// tdc with all corrections ( offset and correction)

	// correction based on channel and tot value
	double corr[ constants::nChannels ];	
	// reference tdc time => the 1st channel on the west side
	double reference;
	stringstream sstr;

	// Alias the values for this event for ease
	for( int j = constants::startWest; j < constants::endEast; j++) {
		
		if ( deadDetector[ j ] ) continue;
		if ( !useDetector[ j ] ) continue;

		tot[ j ] = pico->channelTOT( j );
		tdc[ j ] = pico->channelTDC( j );
		off[ j ] = this->initialOffsets[ j ];	
		tAll[ j ] = tdc[ j ] - off[ j ];
		
		if(tot[ j ] <= minTOT || tot[ j ] > maxTOT) continue;
			
			corr[ j ] = getCorrection( j, tot[ j ] );
			tAll[ j ] -= corr[ j ];
	}

	reference = pico->vpdLeWest[0];
	int start = constants::startWest;
	int end = constants::endWest;

	// loop over West then East and calculate the two sides seperately 
	for ( int sides = 0; sides < 2; sides ++ ){
		if ( 1 == sides){
			start = constants::startEast;
			end = constants::endEast;
		}
		for ( int i = start; i < end; i++ ){

			book->cd( "channel" + ts( i ) );

			double count = 0;
			double avg = 0;
			double c = 0, a = 0; // tmp count and average variables used before cut

			if ( !useDetector[ i ] ) continue;

			if ( westIsGood && eastIsGood ){

				// get the count and average with no cuts
				for ( int j = start; j < end; j++ ){
					if ( useDetector[ j ] && i != j ){
						++c;
						a += tAll[ j ];
					}
				}

				
				if ( c > 0 ){
					// calculate the average
					a /= c;

					// fills the <N> variation within channel
					for ( int j = start; j < end; j++ ){
						if ( useDetector[ j ] && i != j ){
							book->get( iStr + "avgN" )->Fill( c, tAll[ j ] - a );	
						}
					}
				}
				
				// now calculate the count and average with the timing cut 
				for ( int j = start; j < end; j++ ){
					if ( useDetector[ j ] && i != j && c > 0){	
					// now reject events too far from the average time and redetermine count and average					
						if ( 	tAll[ j ] - a > -outlierCut && tAll[ j ] - a < outlierCut ){
							++count;
							avg += tAll[ j ];
						}
					}
				}

				if ( count )
					avg /= count;
				else
					avg = -9999;

				if ( count ){
					book->get( iStr + "cutAvgN" )->Fill( count, tAll[ i ] - avg );
				}

			} // West and East Good
		} // loop det channel
	} // loop sides
}


Double_t calib::detectorResolution(Double_t *x, Double_t *par){
	Double_t resval = 0.0;
	if (x[0]>0){
		resval	= TMath::Sqrt( x[0] / ( 1.0 + x[0] ) );
		return ( par[0] / resval );
	} else {
		return 0.0;
	}
}
