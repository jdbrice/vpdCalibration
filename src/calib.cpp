
#include "constants.h"
#include "calib.h"
#include "histoBook.h"
#include <fstream>
#include <sstream>

// provides my own string shortcuts etc.
using namespace jdbUtils;


/**
 * Constructor - Initializes all of the calibration parameters from the configuration file
 * @param chain       The chain object containing all data compatible with the TOFrPicoDST format
 * @param nIterations Max number of iterations to run
 * @param con         The xml configuration defining key aspects of the calibration
 *					such as number of tot bins to use, data location etc. See repo Readme
 *					for a sample configuration.
 */
calib::calib( TChain* chain, uint nIterations, xmlConfig con )  {
	cout << "[calib.calib] " << endl;
	
	gErrorIgnoreLevel=kError;

	config = con;

	// default number of tot bins is defined in constants in case non is given in config file
	numTOTBins = constants::numTOTBins;
	
	// set the number of tot bins and the range from the config if given 
	numTOTBins = config.getAsInt( "numTOTBins", constants::numTOTBins );
	minTOT = config.getAsDouble( "minTOT", constants::minTOT );
	maxTOT = config.getAsDouble( "maxTOT", constants::maxTOT );

	// set the histogram info verbosity to show nothing
	gStyle->SetOptStat( 0 );
	
	// create the histogram book
	book = new histoBook( ( config.getAsString( "baseName" ) + config.getAsString( "rootOutput" ) ) );
	
	// create a report builder 
	report = new reporter( config.getAsString( "baseName" ) + config.getAsString( "reportOutput" ) );

	

	// now build arrays that need numTOTBins
	for ( int j = 0; j < constants::nChannels; j++){
		correction[ j ] 	= new double[ numTOTBins + 1 ];
		totBins[ j ] 		= new double[ numTOTBins + 1 ];

		deadDetector[ j ]	= false;
	}

	// zero the corrections & offsets
	for ( int j = 0; j < constants::nChannels; j++){
		for (int k = 0; k < numTOTBins + 1; k++){
			correction[ j ] [ k ] = 0;
		}
		initialOffsets[ j ] = 0;
		outlierOffsets[ j ] = 0;
		spline[ j ] = NULL;
	}
	
	// set the maximum number of iterations
	maxIterations = nIterations;


	if ( "paramReport" != config.getAsString( "jobType" ) ){
		// keep the chain variable and make the picoDST var
		_chain = chain;
		pico = new TOFrPicoDst( _chain );
	}

	

	// start at step 0
	currentIteration = 0;

	



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


	// set the spline type
   	// default to akima
   	Interpolation::Type type = ROOT::Math::Interpolation::kAKIMA;
   	useSpline = false;
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

    splineType = type;




    // set the variable types
    xVariable = config.getAsString( "xVariable", "tof-tot" );
    yVariable = config.getAsString( "yVariable", "tof-le" );

    if ( (string)"tof-tot" == xVariable )
    	xLabel = xVariable + " [ns] ";
    else
    	xLabel = xVariable;
    yLabel = yVariable + " [ns] ";


    for ( int i = 0; i < constants::nChannels; i++ ){
    	triggerToTofMap[ i ] = -1;
    	tofToTriggerMap[ i ] = -1;
    }
    if ( 	config.nodeExists( "channelMap" ) ){
    	readTriggerToTofMap();
    }

    mapTriggerToTof = config.getAsBool( "mapTriggerToTof", false );

}

/**
 *	Destructor - Deletes the histoBook ensuring it is saved.
 */
calib::~calib() {
	
	delete book;
	delete report;
	
	for ( int j = 0; j < constants::nChannels; j++){
		delete [] correction[j];
		delete [] totBins[j];
		if ( spline [ j ] )
			delete spline[ j ];
	
	}
	cout << "[calib.~calib] " << endl;
}

/**
 * Generalizes the calibration method to allow any combination of the trigger
 * or time-of-flight electronic's descriminators to be used.
 * @param  channel - the VPD channel for which the value should be retrieved. West = 1-19, East = 20-38
 * @return         returns the channel's timing value
 */
double calib::getX( int channel ) {
	
	if ( mapTriggerToTof ){
		if ( (string)"bbq-adc" == xVariable || (string)"mxq-adc" == xVariable ){
			channel = tofToTriggerMap[ channel ];
		}
	}
	


	if ( (string)"tof-tot" == xVariable )
		return pico->channelTOT( channel );
	else if ( (string)"bbq-adc" == xVariable ){
		return ( (double)pico->bbqADC( channel )  );
	} else if ( (string)"mxq-adc" == xVariable ){
		return ( (double)pico->mxqADC( channel )  );
	}

	// the y varaibles in case you want to do non-standard comparisons
	if ( (string)"tof-tdc" == xVariable )
		return pico->channelTDC( channel );

	// default to the tof tot
	return pico->channelTOT( channel );	
}

/**
 * Gerealizes the calibration method to allow any combination of the trigger
 * or time-of-flight electronic's descriminators to be used
 * @param  channel - the VPD channel for which the value should be retrieved. West = 1-19, East = 20-38
 * @return         returns the channel's timing value
 */
double calib::getY( int channel ){

	if ( mapTriggerToTof ){
		if ( (string)"bbq-tdc" == yVariable || (string)"mxq-tdc" == yVariable ){
			channel = tofToTriggerMap[ channel ];
		}
	}

	if ( (string)"tof-le" == yVariable )
		return pico->channelTDC( channel );
	else if ( (string)"bbq-tdc" == yVariable ){
		return ( (double)pico->bbqTDC( channel ) * constants::tacToNS );
	} else if ( (string)"mxq-tdc" == yVariable ){
		return ( (double)pico->mxqTDC( channel ) * constants::tacToNS );
	}

	// the x varaibles in case you want to do non-standard comparisons
	if ( (string)"tof-tot" == yVariable )
		return pico->channelTOT( channel );

	// default to the tof tot
	return pico->channelTDC( channel );		
}

/**
 *	Offsets
 *	Calculates the initial offsets for each channel with respect to channel 1 on the west side.
 *	Then performs the offset calculation again after all corrections have been applied
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
	if ( ! book->get( "tdc" ) ){
		book->make2D( "tdc", yVariable + " relative to West Channel 1; Detector ; " + yLabel, constants::nChannels, -0.5, constants::nChannels-0.5, 2000, -100, 100 );
		book->make1D( "tdcMean", yVariable + " relative to West Channel 1; Detector ; " + yLabel, constants::nChannels, -0.5, constants::nChannels-0.5 );
		book->make2D( "correctedOffsets", "Corrected Initial Offsets", constants::nChannels, -0.5, constants::nChannels-0.5, 2000, -100, 100 );
		book->make2D( "tdcRaw", "All tdc Values ", constants::nChannels, 0, constants::nChannels, 1000, 0, 51200 );	
		book->make2D( "tdcOffsetRemoved", yVariable + " relative to West Channel 1; Detector ; " + yLabel, constants::nChannels, -0.5, constants::nChannels-0.5, 2000, -100, 100 );
	}
	


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
    	double reference = getY( 0 );
    	

		for( int j = constants::startWest; j < constants::endEast; j++) {

			
			// skip dead detectors
			if ( deadDetector[ j ] ) continue;

			int nHits = pico->numHits( j );
			
			if ( nHits < constants::minHits ) 
				continue;

			double tdc = getY( j );
	    	double tot = getX( j );

	    	book->fill( "tdcRaw", j, tdc );

	    	if(tot <= minTOT || tot >= maxTOT) continue;	    


		    book->fill( "tdc", j, tdc - reference );

		}	
	} // end loop on events

	// calculate the offsets
  	TH2D* tdc = (TH2D*) book->get( "tdc" );

	for ( int i = constants::startWest; i < constants::endEast; i++ ){
		TH1D* tmp = tdc->ProjectionY( "tmp", i+1, i+1 );

		double max = tmp->GetBinCenter( tmp->GetMaximumBin() );
		double rms = 1.2*tmp->GetRMS();
		tmp->GetXaxis()->SetRangeUser( max - rms, max + rms  );

		if ( i == constants::startWest )
			this->initialOffsets[ i ] = 0;
		else
			this->initialOffsets[ i ] = tmp->GetMean();

		cout << "Channel [ " << i+1 << " ] Offset = " << this->initialOffsets[ i ] << " ns " << endl;

		book->get( "tdcMean" )->SetBinContent( i+1, this->initialOffsets[ i ] );
		book->get( "tdcMean" )->SetBinError( i+1, tmp->GetMeanError() );

		delete tmp;
	}

	// loop over all events to draw them with offsets removed
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
    	double reference = getY( 0 );
    	

		for( int j = constants::startWest; j < constants::endEast; j++) {

			
			// skip dead detectors
			if ( deadDetector[ j ] ) continue;

			int nHits = pico->numHits( j );
			
			if ( nHits < constants::minHits ) 
				continue;

			double tdc = getY( j );
	    	double tot = getX( j );

	    	if(tot <= minTOT || tot >= maxTOT) continue;	    

		    book->fill( "tdcOffsetRemoved", j, tdc - reference - initialOffsets[ j ] );

		}	
	} // end loop on events

	// calculate what the final mean of the west side channels would be if no offsets where removed.
	// It is the average using equal weight for each channel 
	double totalWest = 0;
	double numWest = 0;
	for ( int k = constants::startWest; k < constants::endWest; k++){
		if ( deadDetector[ k ] ) continue;
		totalWest += this->initialOffsets[ k ];
		numWest++;
	}
	finalWestOffset = (totalWest / numWest );
	cout << "Channel 0 floating offset : " << finalWestOffset << " ns " << endl;

	double totalEast = 0;
	double numEast = 0;
	for ( int k = constants::startEast; k < constants::endEast; k++){
		if ( deadDetector[ k ] ) continue;
		totalEast += this->initialOffsets[ k ];
		numEast++;
	}
	double eastOffset = (totalEast / numEast );
	eastWestOffset = (finalWestOffset - eastOffset);
	cout << "West - East : " << (finalWestOffset - eastOffset) << " ns " << endl;

	

	// get the east / west offset just for plotting
	double westMinusEast = 0;
	TH1D* west = tdc->ProjectionY( "westOffset", constants::startWest+2, constants::endWest );
	book->add( "westOffset", west );
	TH1D* east = tdc->ProjectionY( "eastOffset", constants::startEast+1, constants::endEast );
	book->add( "eastOffset", east );

	westMinusEast = ( west->GetMean() - east->GetMean() );
	
	cout << "West Mean = " << west->GetMean() << " ns " << endl;
	cout << "East Mean = " << east->GetMean() << " ns " <<endl;
	cout << "West - East offset: " << westMinusEast << " ns " << endl; 

	report->newPage( 1, 2);


	book->style( "tdc" )
		->set( "range", -35.0, 35.0 )->draw();
	book->style( "tdcMean" )
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

void calib::finalOffsets() {

	startTimer();

	if ( !_chain ){
		cout << "[calib." << __FUNCTION__ << "] ERROR: Invalid chain " << endl;
		return;
	}


	Int_t nevents = (Int_t)_chain->GetEntries();
	cout << "[calib." << __FUNCTION__ << "] Loaded: " << nevents << " events " << endl;

	book->cd( "finalOffset" );
	
	// make all the histos the first round
	if ( ! book->get( "tdc" ) ){
		book->make2D( "tdc", yVariable + " relative to West Channel 1; Detector ; " + yLabel, constants::nChannels, -0.5, constants::nChannels-0.5, 2000, -100, 100 );
		book->make1D( "tdcMean", yVariable + " relative to West Channel 1; Detector ; " + yLabel, constants::nChannels, -0.5, constants::nChannels-0.5 );
		book->make2D( "correctedOffsets", "Corrected Initial Offsets", constants::nChannels, -0.5, constants::nChannels-0.5, 2000, -100, 100 );
		book->make2D( "tdcRaw", "All tdc Values ", constants::nChannels, 0, constants::nChannels, 1000, 0, 51200 );	
		book->make2D( "tdcOffsetRemoved", yVariable + " relative to West Channel 1; Detector ; " + yLabel, constants::nChannels, -0.5, constants::nChannels-0.5, 2000, -100, 100 );
	}
	


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
    	double reference = getY( 0 ) - getCorrection( 0, getX( 0 ) );
    	

		for( int j = constants::startWest; j < constants::endEast; j++) {

			
			// skip dead detectors
			if ( deadDetector[ j ] ) continue;

			int nHits = pico->numHits( j );
			
			if ( nHits < constants::minHits ) 
				continue;
			double tot = getX( j );
			double tdc = getY( j ) - getCorrection( j, tot );

	    	if(tot <= minTOT || tot >= maxTOT) continue;	    


		    book->fill( "tdc", j, tdc - reference );

		}	
	} // end loop on events

	// calculate the offsets
  	TH2D* tdc = (TH2D*) book->get( "tdc" );

	for ( int i = constants::startWest; i < constants::endEast; i++ ){
		TH1D* tmp = tdc->ProjectionY( "tmp", i+1, i+1 );

		double max = tmp->GetBinCenter( tmp->GetMaximumBin() );
		double rms = 1.2*tmp->GetRMS();
		tmp->GetXaxis()->SetRangeUser( max - rms, max + rms  );

		if ( i == constants::startWest )
			this->initialOffsets[ i ] = 0;
		else
			this->initialOffsets[ i ] = tmp->GetMean();

		//if ( i >= constants::startEast && i < constants::endEast ){

		//}

		cout << "Channel [ " << i+1 << " ] Offset = " << this->initialOffsets[ i ] << " ns " << endl;

		book->get( "tdcMean" )->SetBinContent( i+1, this->initialOffsets[ i ] );
		book->get( "tdcMean" )->SetBinError( i+1, tmp->GetMeanError() );

		delete tmp;
	}

	cout << "[calib." << __FUNCTION__ << "] completed in " << elapsed() << " seconds " << endl;
}

/**
 * Produces the bins to use when histogramming the values in TOT space
 * @param variableBinning 
 *        True  	Calculates variable binning for tot space such that 
 *        			the number of events is roughly equal for each tot bin
 *        False 	Fixed binning the tot space from minTOT to maxTOT
 */
void calib::binTOT( bool variableBinning ) {

	cout << "[calib." << __FUNCTION__ << "] Starting " << endl;

	if ( variableBinning )
		cout << "[calib." << __FUNCTION__ << "] Variable Binning TOT Range :  " << minTOT << " -> " << maxTOT << endl;
	else
		cout << "[calib." << __FUNCTION__ << "] Fixed Binning TOT Range :  " << minTOT << " -> " << maxTOT << endl;

	cout << "[calib." << __FUNCTION__ << "] Using " << numTOTBins << " bins for TOT" << endl;

	startTimer();



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
	        	Double_t tot = getX( j );
	          
	        	if(tot > minTOT && tot < maxTOT ) 
	          		tots[j].push_back(tot);
	        }

	    }

  		if( numEast > constants::minHits ){
    
    		for(Int_t j = constants::startEast; j < constants::endEast; j++) {
      			Double_t tot = getX( j );
      
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

      		if ( variableBinning ){
	      		
	      		Int_t step = size / (numTOTBins + 1 ); 
	    
	      		// sort into ascending order
	        	std::sort( tots[i].begin(), tots[i].end());
	        	
	        	totBins[ i ][0] = minTOT;
	        	totBins[ i ][ numTOTBins ] = maxTOT;
	        	
	        	for( Int_t j = 1; j < numTOTBins ; j++) {

	        		double d1 = tots[i].at( step * j );
	        		totBins[ i ][ j ] = d1;
	            	
	        	}	// loop over tot bins
        	}	// end variable binning 
        	else { // fixed binning

				for ( int s = 0; s <= numTOTBins; s++ ){
					double edge = ((maxTOT - minTOT) / (double) numTOTBins) * s;
					edge += minTOT;
					totBins[ i ][ s ] = edge;
				}
			
        	}

      } // end channle not dead

  	} // end loop channles
  	
	for(Int_t i = 0; i< constants::nChannels; i++) {
		tots[i].clear();
	}	
	
	cout << "[calib." << __FUNCTION__ << "] completed in " << elapsed() << " seconds " << endl;
}

/**
 * Retrieves the Slewing correction for a given tot value for a given channel
 * @param  vpdChannel VPD Channel for the slewing correction 
 * @param  tot        The TOT value for the correction
 * @return            The TDC correction for the channel at the given tot value
 */
double calib::getCorrection( int vpdChannel, double tot ){
	

	// use splines to get the correction value if set to
	if ( useSpline && spline[ vpdChannel ] && spline[ vpdChannel ]->getSpline() ){
		//return spline[ vpdChannel ]->getSpline()->Eval( tot );
		return spline[ vpdChannel ]->eval( tot );
	}
	
	// if not fall back to doing bin based corrections
	int totBin = binForTOT( vpdChannel, tot ); 
	return correction[ vpdChannel ][ totBin ];

}

/**
 * Determines which bin corresponds to a given TOT value
 * @param  vpdChannel The Channel whose binning should be used
 * @param  tot        The TOT value
 * @return            Returns the bin for the given tot value
 */
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

/**
 * Performs the outlier rejection calculations for each event.
 * Calculates the VPD zVertex for all combinations of east and west detectors
 * and uses the detectors only in pairs which are within a given cut from the TPC
 * zVertex
 * @param reject 
 *        True 		Performs outlier reject
 *        False 	Uses all detectors, all events
 */	
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

		double tdcWest = getY( j );
	    double totWest = getX( j );

	    tdcWest -= (this->initialOffsets[ j ] + this->outlierOffsets[ j ]);

	  
	    if( totWest <= minTOT || totWest > maxTOT) continue;
	    
	    double corWest = getCorrection( j, totWest );
	    tdcWest -= corWest;

	    sumWest += tdcWest;
	    countWest++;

		for ( int k = constants::startEast; k < constants::endEast; k++ ){
			
			if ( deadDetector[ k ] ) continue;

			double tdcEast = getY( k );
	    	double totEast = getX( k );

	    	tdcEast -= (this->initialOffsets[ k ] + this->outlierOffsets[ k ]);

	    	if( totEast <= minTOT || totEast > maxTOT) continue;
	    	
	    	double corEast = getCorrection( k, totEast );
	    	tdcEast -= corEast;

	    	if ( j == constants::startWest ){
	    		sumEast += tdcEast;
	    		countEast++;
	    	}

	    	// calculate the VPD z Vertex
	    	double vpdZ = constants::c * ( tdcEast - tdcWest) / 2.0;
	    	if ( (string)"bbq-tdc" == yVariable || (string)"mxq-tdc" == yVariable )
	    		vpdZ = constants::c * ( tdcWest - tdcEast) / 2.0;

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
		if ( (string)"bbq-tdc" == yVariable || (string)"mxq-tdc" == yVariable )
			vpdZ = constants::c * ( (sumWest/countWest) - (sumEast/countEast)) / 2.0;	
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

/**
 * Prepares the histograms for each step of the calibration
 */
void calib::prepareStepHistograms() {
	
	// for names
	string iStr = "it"+ts(currentIteration);
	// for titles
	string step = "Step " + ts( currentIteration+1 ) + " : ";


	for ( int ch = constants::startWest; ch < constants::endEast; ch++ ){
		
		book->cd( "channel" + ts(ch) );

		// make channel titles start at 1
		string sCh = "Channel "+ts(ch+1);
		string title2D = step + sCh + " " + yVariable + " vs " + xVariable + ";" + xLabel + ";" + yLabel ;
		string title1D = step + sCh + " " + yVariable + ";" + yLabel + "; [ # ] "  ;

		book->make2D( 	iStr + "tdctot", 	title2D, numTOTBins , totBins[ ch ], 1000, -40, 40 );
		book->make2D( 	iStr + "tdccor", 	title2D, numTOTBins , totBins[ ch ], 1000, -20, 20 );
		book->make1D( 	iStr + "tdc", 		title1D, 500, -10, 10 );
		book->make2D( 	iStr + "avgN", 		step + sCh + " : 1 - <N>;# of Detectors;" + yLabel, 
						constants::nChannels/2, 1, constants::nChannels/2, 1000, -20, 20 );
		book->make2D( 	iStr + "cutAvgN", 	step + sCh + " : 1 - <N>;# of Detectors;" + yLabel, 
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
	book->make2D( 	iStr + "Offsets", step + yLabel + " wrt West Channel 1; Detector ; [#] ",
							constants::nChannels, -0.5, constants::nChannels-0.5, 2000, -100, 100 );

	cout << "[calib." << __FUNCTION__ << "[" << currentIteration << "]] " << " Histograms Booked " << endl;

}


/**
 * Performs a single step of the slewing calibration. 
 * Loops over all events, performs the outlier rejection and averageN calulations.
 * Then each channel is looped over. For each channel the average of all times from other
 * channels is calculated. The slewing curve for each channel is then plotted against the 
 * average from other channels.
 */
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

	cout << "[calib." << __FUNCTION__ << "[" << currentIteration << "]] " << " Calibrating " << endl;

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

    		tot[ j ] = getX( j );
    		tdc[ j ] = getY( j );
    		if ( removeOffset )
   	 			off[ j ] = this->initialOffsets[ j ];
   	 		else 
   	 			off[ j ] = 0;
   			tAll[ j ] = tdc[ j ] - off[ j ];

    		if(tot[ j ] <= minTOT || tot[ j ] >= maxTOT) continue;
  			
  			corr[ j ] = getCorrection( j, tot[ j ] );
  			tAll[ j ] -= corr[ j ];
    	}
    	reference = getY( 0 ) - getCorrection( 0, tot[ 0 ] );

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

void calib::checkStep( ) {

	cout << "[calib." << __FUNCTION__ << "[" << currentIteration << "]] " << " Start " << endl;
	
	string iStr = "it"+ts(currentIteration);
	int zBins = 600, zRange = 200;

	book->cd( "Vertex_Z");
	book->make1D( 	iStr+"all", iStr + "z_{TPC} - z_{VPD}; [#]", zBins, -zRange, zRange );
	book->make1D( 	iStr+"avg", "TPC vs. VPD z Vertex using <East> & <West>; z_{TPC} - z_{VPD} [cm]; [#]", 	zBins, -zRange, zRange );
	book->make2D( 	iStr+"zTPCzVPD", "TPC vs. VPD z Vertex; z_{TPC};z_{VPD}", zBins/2, -zRange/2, zRange/2, zBins/2, -zRange/2, zRange/2 );
	book->make2D( 	iStr+"zTPCzVPDAvg", "TPC vs. VPD z Vertex using <East> & <West>; z_{TPC} [cm];z_{VPD} [cm]", zBins/2, -zRange/2, zRange/2, zBins/2, -zRange/2, zRange/2 );

	startTimer();

	for ( int j = constants::startWest; j < constants::endEast; j++ ){
		//cout << "Offset[ " << j << " ] = " << initialOffsets[ j ] << endl;
	}

	cout << "[calib." << __FUNCTION__ << "[" << currentIteration << "]] " << " Calibrating " << endl;

	double avgCountEast = 0;
	double neEast = 0;
	double avgCountWest = 0;
	double neWest = 0;

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

    	double sumEast = 0;
		double sumWest = 0;
		double countEast = 0;
		double countWest = 0;
		stringstream sstr;

		for ( int j = constants::startWest; j < constants::endWest; j++ ){

			double tdcWest = getY( j );
		    double totWest = getX( j );
		  
		    if( totWest <= minTOT || totWest > maxTOT) continue;
		    
		    double corWest = getCorrection( j, totWest ) +initialOffsets[ j ];
		    tdcWest -= corWest;

		    sumWest += tdcWest;
		    countWest++;

			for ( int k = constants::startEast; k < constants::endEast; k++ ){

				double tdcEast = getY( k );
		    	double totEast = getX( k );

		    	if( totEast <= minTOT || totEast > maxTOT) continue;
		    	
		    	double corEast = getCorrection( k, totEast ) +initialOffsets[ k ];
		    	tdcEast -= corEast;

		    	if ( j == constants::startWest ){
		    		sumEast += tdcEast;
		    		countEast++;
		    	}

		    	// calculate the VPD z Vertex
		    	double vpdZ = constants::c * ( tdcEast - tdcWest) / 2.0;
		    	if ( (string)"bbq-tdc" == yVariable || (string)"mxq-tdc" == yVariable )
		    		vpdZ = constants::c * ( tdcWest - tdcEast) / 2.0;

		    	book->get( iStr+"all" )->Fill( tpcZ - vpdZ );
		    	book->get( iStr+"zTPCzVPD" )->Fill( tpcZ, vpdZ );
			    	  			
			} // loop channel k
		} // loop channel j

		double ovc = 1.0 / 29.9792458;
		if ( countEast >= 1 && countWest >= 1){

			//cout << "SumWest : " << sumWest << endl ;
			//cout << "" << (ovc * -1.29 * 2 * countWest) << endl;
			//sumWest = sumWest + (ovc * 3.35 * 2 * countWest);
			double vpdZ = constants::c * ( (sumEast/countEast) - (sumWest/countWest)) / 2.0;	
			//if ( (string)"bbq-tdc" == yVariable || (string)"mxq-tdc" == yVariable )
			//	vpdZ = constants::c * ( (sumWest/countWest) - (sumEast/countEast)) / 2.0;	
			book->fill( iStr+"zTPCzVPDAvg", tpcZ, vpdZ );
			book->fill( iStr+"avg", ( tpcZ-vpdZ ));

			

		}

		if ( countEast >= 1 ){
			neEast++;
    		avgCountEast += countEast;
		}
		if ( countWest >= 1){
			neWest++;
			avgCountWest += countWest;
		}
    	
    	
 		   	

    		
	}


	cout << " Avg Count East " << (avgCountEast / neEast ) << endl;
	cout << " Avg Count West " << (avgCountWest / neWest ) << endl;

	cout << "[calib." << __FUNCTION__ << "[" << currentIteration << "]] " << " completed in " << elapsed() << " seconds " << endl;
	
	// report the outlier rejection
	report->newPage(2, 2);
	
	gPad->SetLogy(1);
	book->style( iStr+"all" )
		->set( "numberOfTicks", 5, 5)
		->draw();
	

	report->cd( 2, 1 );
	gPad->SetLogy(1);
	book->style( iStr+"avg" )
		->set( "numberOfTicks", 5, 5)
		->set( "draw", "" )
		->set( "domain", -100, 100)
		->draw();

	report->cd( 1, 2 );
	book->style( iStr+"zTPCzVPD" )
		->set( "numberOfTicks", 5, 5)
		->set( "draw", "colz" )->draw();

	report->cd( 2, 2 );
	book->style( iStr+"zTPCzVPDAvg" )
		->set( "numberOfTicks", 5, 5 )
		->set( "draw", "colz" )

		->draw();
	report->savePage();

	
	book->clearLegend();


	/** Determine the Offset and contrain it for the outlier rejection*/
	double toff = book->get( iStr+"avg") -> GetMean();
	cout << "TPC - VPD = " << toff << endl;
	double ovc = 1.0 / constants::c;
	double nOff = ovc * toff * 2.0;
	cout << "Channel [ west ] += " <<nOff << endl;
	for ( int i = constants::startWest; i < constants::endWest; i++ ){
		this->initialOffsets[ i ] += nOff;
	}

	currentIteration++;

}

/**
 * Performs the fitting after all calibration steps have been completed
 * and outputs the report plots for each detectors resolution.
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
										->set( "numberOfTicks", 10, 2 )
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

		text->DrawLatex(0.25,0.86, ("#sigma = " + ts( fr->GetParameter( 0 ) ) + " [ns] "  ).c_str() );
		text->DrawLatex(0.25,0.81, ("#chi^{2}/DOF = " + ts( chiDoF) ).c_str() );

		
	    if ( j == constants::endEast - 1){
	    	report->savePage();
	    }

	}

	report->newPage();

	double max = sigmas->GetMaximum();

	TF1 * avgFit = new TF1( "avgFit", "pol0", 0, 39 );


	book->cd ( "final" );
	book->add( "detSigma", sigmas );
	book->get( "detSigma" )->Fit( "avgFit", "QR" );
	book->style( "detSigma" )
		->set( "markerStyle", 17)
		->set( "markerColor", 2 )
		->set( "range", 0, max + .015 )
		->draw();

	text->DrawLatex(0.25,0.25, ("Average #sigma = " + ts( avgFit->GetParameter( 0 ) ) + " [ns] ").c_str() );

	report->savePage();	

}

/**
 * Executes the given number of calibration steps and then
 * calls the finish() function to fit the single detector resolution.
 */
void calib::loop( ) {

	for ( unsigned int i = 0; i < maxIterations; i++ ){
		step();
	}
	finish();

}


/**
 * After each calibration step the corrections are calculated from the slewing curves.
 * If Splines are used the spline is perpared and both the slewing curve and spline are drawn.
 */
void calib::makeCorrections( ){

	cout << "[calib." << __FUNCTION__ << "[" << currentIteration << "]] " << " Start " << endl;
	startTimer();

	bool removeOffset = config.getAsBool( "removeOffset" );

	string iStr = "it" + ts( currentIteration );
	// Bins with greater error on the fitslicesY mean will not be used in final correction
	double maxError = config.getAsDouble( "binMaxError", 0.10);
	
	
	report->newPage( 4, 5);

	// get the corrections for the next iteration 
	for( int k = constants::startWest; k < constants::endEast; k++) {


		TF1* g = new TF1( "g", "gaus", -5, 5 );
		if ( !removeOffset )
			g->SetRange( -100, 100 );

		if ( currentIteration <= 1 || currentIteration == maxIterations - 1){
			if ( deadDetector[ k ] ){
				report->next();
				if ( k == constants::endWest - 1 || constants::endEast - 1 == k){
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
	    pre->FitSlicesY( g );

	    TH1D* preMean = (TH1D*) gDirectory->FindObject( 
	    					(iStr + "tdctot_1").c_str() );

	    // do the fit
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

	    	if ( removeOffset ){
				// reject low statistics bins and instead interpolate between good bins
		    	if ( ib > 1){
			    	if ( preMean->GetBinError( ib ) < maxError && preMean->GetBinError( ib ) != 0 )
			    		reset = cor->GetBinContent( ib );
			    	else {
			    		cor->SetBinContent( ib, reset );
			    		dif->SetBinContent( ib, 0 );
			    	}
		    	} 
	    	}

	    	if ( currentIteration >= 1 && TMath::Abs( dif->GetBinContent( ib ) ) > 1 ){
	    		dif->SetBinContent( ib, 0 );
	    	}

	    	correction[ k ][ ib  ] = cor->GetBinContent( ib );
	    	
	    }


	    if ( spline[ k ])
	    	delete spline[ k ];
	    
	    if ( useSpline )
		    spline[ k ] = new splineMaker( cor, splineAlignment::center, splineType );

		// make a spline for drawing
		splineMaker* vSpline;
		vSpline = new splineMaker( dif, splineAlignment::center, splineType );
	    
	    if ( currentIteration <= 1 || currentIteration == maxIterations - 1){
		    
		    if ( removeOffset )
		    	book->style( ("it"+ts(currentIteration)+"tdccor") )->set( "range", -5.0, 5.0);
		    else
		    	book->style( ("it"+ts(currentIteration)+"tdccor") )->set( "dynamicdomain", 1, -1, -1, 2 );
		    	
		    if ( xVariable.find( "adc" ) != string::npos )
		    	book->style( ("it"+ts(currentIteration)+"tdccor") )->set( "numberOfTicks", 5, 5);

		    post->Draw( "colz" );
		    
		    if ( useSpline ){
		    	TGraph* g = vSpline->graph( minTOT, maxTOT, 0.2);
		    	g->GetYaxis()->SetRangeUser( -5, 5);
		    	g->SetMarkerStyle(7);
		    	g->SetMarkerColor( kRed );
		    	g->Draw( "same cp" );
			}
			
		    if ( k == constants::endEast - 1 || k == constants::endWest - 1){
		    	report->savePage();
		    	report->newPage( 4, 5);
	    	}

	    	// delete the spline for drawing
			delete vSpline;
    	}
	    

	}

	cout << "[calib." << __FUNCTION__ << "[" << currentIteration << "]] " << " completed in " << elapsed() << " seconds " << endl;
}

/**
 * Outputs the slewing curve corrections in DB format
 * Adds the offsets back as well as the channel 1 relative offset
 * so that the corrected time is true to the initial time.
 */
void calib::writeParameters(  ){

	cout << "[calib." << __FUNCTION__ << "[" << currentIteration << "]] " << " Start " << endl;
	cout << " min, max, nBins : " << minTOT << ",  " << maxTOT << ", " << numTOTBins << endl;

	string outName = config.getAsString( "baseName" ) + config.getAsString( "paramsOutput", "params.dat" );

	ofstream f;
	f.open( outName.c_str() );

	bool removeOffset = config.getAsBool( "removeOffset" );

	book->cd( "correctionParameters" );

	for ( int j = constants::startWest; j < constants::endEast; j++ ){
		

		string title = "Channel "+ts(j+1)+ " : Slewing Correction ;" + xLabel + ";" + yLabel;
		string sTitle = "Channel "+ts(j+1)+ " : Spline Slewing Correction ;" + xLabel + ";" + yLabel;

		book->make1D( "slewingCor" +ts(j), title, numTOTBins, totBins[ j ]  );
		if ( useSpline )
			book->make1D( "splineSlewingCor" +ts(j), sTitle, numTOTBins, totBins[ j ]  );
		book->make1D( "slewingNoOffset" +ts(j), title, numTOTBins, totBins[ j ] );

		f << (j + 1) << endl;
		f << numTOTBins << endl;

		for ( int i = 0; i <= numTOTBins; i++ ){
			f << totBins[ j ][ i ] << " ";
		}
		f << endl;
		for ( int i = 0; i <= numTOTBins; i++ ){
			double totStep = (( maxTOT - minTOT ) / (double)(numTOTBins));
			double tot = minTOT + totStep * (double)(i);

			if ( 0 == i )
				tot = minTOT;
			else if ( numTOTBins == i )
				tot = maxTOT;
			else
				tot = (totBins[ j ][ i ] + totBins[ j ][ i + 1 ] ) / 2.0;

			double off = 0;
			if ( removeOffset ){
				off = initialOffsets[ j ] - finalWestOffset;
			} else {
				if ( j >= constants::startEast && j < constants::endEast )
					off = 0 - eastWestOffset;
			}
			
			// bin based corrections
			double bCor = 0;
			
			if ( !deadDetector[ j ] )
				bCor = correction[ j ][ i + 1 ] + off;
			
			book->get( "slewingCor" +ts(j) )->SetBinContent( i+1, bCor );

			

			// spline based corrections
			double sCor = 0;
			if ( !deadDetector[ j ] )
				sCor = spline[ j ]->eval( tot ) + off;

			book->get( "splineSlewingCor" +ts(j) )->SetBinContent( 
				book->get( "splineSlewingCor" +ts(j) )->GetXaxis()->FindBin( tot ), sCor );

			if ( useSpline ){
				f << ( sCor ) << " ";
				book->get( "slewingNoOffset" +ts(j))->SetBinContent( i, sCor - off );
			} else {	// use the bin based corrections only
				f << ( bCor ) << " ";
				book->get( "slewingNoOffset" +ts(j))->SetBinContent( i, bCor - off );
			}
		}
		f << endl;
	} 

	f.close();


	//draw the parameters
	report->newPage( 3, 4 );

	for ( int j = constants::startWest; j < constants::endEast; j++ ){


		if ( j != constants::startEast && j != constants::startWest )
			report->next();

		if ( removeOffset )
			book->style( "slewingNoOffset"+ts(j) )->set("lineColor", kRed )->set( "range", -5, 5)->draw();
		else
			book->style( "slewingNoOffset"+ts(j) )->set("lineColor", kRed )->draw();

		if ( xVariable.find( "adc" ) != string::npos )
			book->style( "slewingNoOffset"+ts(j) )->set( "numberOfTicks", 5, 5);

		if ( j == constants::endWest - 1 || j == constants::endEast - 1){
			report->savePage();
			report->newPage( 3, 4 );
		}

	}
	cout << "[calib." << __FUNCTION__ << "] " << " completed in " << elapsed() << " seconds " << endl;

}


/**
 * Reads parameter files in either to compare or to run
 * and produce calibration plots from the existing parameter files
 */
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
					cout << "tBins " << tBins << endl;
					if ( 	channel >= 1 && channel <= 38 && 
							tBins == numTOTBins ){

						double tmp = 0;
						for ( int j=0; j <= tBins; j++ ){
							infile >> tmp;
							totBins[ channel - 1 ][ j ] = tmp;
							if ( channel == 1 )
								cout << " bin[ " << j << " ] = " << tmp << endl;
						}
						
						string title = "Channel "+ts(channel)+" Slewing Corrections; " + xLabel + ";" + yLabel ;
						// create a hist
						book->cd( "params" );
						book->make1D( "file"+ts(fi)+"channel"+ts(channel-1), title, tBins, totBins[ channel-1] );


						for ( int j=0; j <= tBins; j++ ){
							infile >> tmp;
							correction[ channel - 1 ][ j ] = tmp;
							//cout << "Cor[ " << channel - 1 << " ][ " << j << " ] = " << tmp << endl;
							book->get( "file"+ts(fi)+"channel"+ts(channel-1) )->SetBinContent( j+1, tmp );
						}
						// TODO : tmp workaround for variable binning bug
						//correction[ channel - 1 ][ 0 ] = correction[ channel - 1 ][ 1 ];
						for ( int j=0; j <= tBins; j++ ){ 
							book->get( "file"+ts(fi)+"channel"+ts(channel-1) )->SetBinContent( j+1, correction[ channel - 1 ][ j ] /*- correction[ channel - 1 ][ 5 ]*/ );
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
				
					if ( files.size() == 1 && (string)"checkParams" == config.getAsString( "jobType" ) ){
						spline[ k ] = new splineMaker( (TH1D*)book->get("file"+ts(fi)+"channel"+ts(k)), splineAlignment::left, splineType );
					}
					
					
				    book->style( "file"+ts(fi)+"channel"+ts(k) )
				    	->set( "lineColor", 1+fi );
				    if ( 0 == fi ){
				    	book->clearLegend();
					    book->style( "file"+ts(fi)+"channel"+ts(k) )
					    	->set( "domain", minTOT, maxTOT )
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
			currentIteration = 0;
		}
	}

	cout << "[calib." << __FUNCTION__ <<  "] " << " completed in " << elapsed() << " seconds " << endl;


	cout << " min, max, nBins : " << minTOT << ",  " << maxTOT << ", " << numTOTBins << endl;
	for ( int i = 0; i < 30; i ++ ){
		double totStep = (( maxTOT - minTOT ) / (double)(numTOTBins));
		double tot = minTOT + totStep * (double)(i);
		cout << getCorrection( 0, tot );
	}
}


/**
 * Outputs the outlier rejection and slewing curve summary figures to a pdf report
 */
void calib::stepReport() {

	string iStr = "it"+ts(currentIteration);

	if ( config.getAsBool( "outlierRejection" ) == false )
		return;

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
	book->style( iStr+"All" )
		->set( "numberOfTicks", 5, 5)
		->draw();
	fill->Draw("same" );

	report->cd( 2, 1 );
	gPad->SetLogy(1);
	book->style( iStr+"avg" )
		->set( "numberOfTicks", 5, 5)
		->set( "draw", "" )
		->set( "domain", -vzCut, vzCut)
		->draw();

	report->cd( 1, 2 );
	book->style( iStr+"zTPCzVPD" )
		->set( "numberOfTicks", 5, 5)
		->set( "draw", "colz" )->draw();

	report->cd( 2, 2 );
	book->style( iStr+"zTPCzVPDAvg" )
		->set( "numberOfTicks", 5, 5 )
		->set( "draw", "colz" )

		->draw();
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

	
	//book->style( iStr+"avg" )
		//->set( "domain", -100, 100);
	/** Determine the Offset and contrain it for the outlier rejection*/
	double toff = book->get( iStr+"avg") -> GetMean();
	cout << " VPD - TPC = " << toff << endl;
	double ovc = 1.0 / constants::c;
	double nOff = ovc * toff * 2;
	cout << "Channel [ west ] += " <<nOff << endl;
	for ( int i = constants::startWest; i < constants::endWest; i++ ){
		this->initialOffsets[ i ] += nOff;
	}



}

/**
 * Plots the 1 - < N > distribution for all channels
 */
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

		tot[ j ] = getX( j );
		tdc[ j ] = getY( j );
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

/**
 * The functional form for the difference between 
 * sampling a gaussian and the mean of n samples of the same gaussian.
 * Note that especially in high multiplicity the low N range of the fit 
 * Doesn not follow the fit closesly since detector-to-detector variations play a larger role.
 *
 * Exact form would include weighting each detector by its sigma (unknown). 
 * However this form converges to the same value. For large enough statistics this isn't
 * a poblem
 * 
 * @param  x   The number of gaus distributions being sampled
 * @param  par The single detector sigma
 * @return     measured variance in the 1 - <N>
 */
Double_t calib::detectorResolution(Double_t *x, Double_t *par){
	Double_t resval = 0.0;
	if (x[0]>0){
		resval	= TMath::Sqrt( x[0] / ( 1.0 + x[0] ) );
		return ( par[0] / resval );
	} else {
		return 0.0;
	}
}

/**
 * Reads in the trigger to tof map produced by the DB reader for the vpdTriggerToTofMap Table
 * Used to compare the trigger (bbq and mxq) channels directly to the tof channels
 */
void calib::readTriggerToTofMap(){

	cout << "[calib." << __FUNCTION__ << "] " << " Start " << endl;
	startTimer();

	// maximum possible number of trigger tubes per side of vpd
  	// currently there are only 16 PMTs active on the trigger side
	short 	eastTriggerToTofMap[ constants::nhChannels ],
	    	westTriggerToTofMap[ constants::nhChannels ];

	for(int i = 0; i < constants::nhChannels; i++ ) {
		eastTriggerToTofMap[ i ] = -1;
		westTriggerToTofMap[ i ] = -1;
	}

	ifstream inData;
	inData.open( config.getAsString( "channelMap" ).c_str() );

	string eastWest = "";
	short nPMTs = 0;
	inData >> eastWest;    // "east" or "west" lowercase
	inData >> nPMTs;

	for ( short j = 0; j < nPMTs; j++ ) {
		short triggerIndex = -1;
		short tofIndex = -1;

		inData >> triggerIndex >> tofIndex;

		// ensure that the trigger index is valid
		if ( triggerIndex >= 0 ){
		  
		  if ( "east" == eastWest )
		    eastTriggerToTofMap[ triggerIndex  ] = tofIndex ;
		  else if ( "west" == eastWest )
		    westTriggerToTofMap[ triggerIndex  ] = tofIndex ;
		    

		} else {
		  cout << " Malformed map file " << endl;
		  return;
		}

	}

	// now read the other side
	eastWest = "";

	inData >> eastWest;    // "east" or "west" lowercase
	inData >> nPMTs;
	for ( short j = 0; j < nPMTs; j++ ) {
		short triggerIndex = -1;
		short tofIndex = -1;

		inData >> triggerIndex >> tofIndex;

		// ensure that the trigger index is valid
		if ( triggerIndex >= 0 ){
		  
		  if ( "east" == eastWest )
		    eastTriggerToTofMap[ triggerIndex  ] = tofIndex ;
		  else if ( "west" == eastWest )
		    westTriggerToTofMap[ triggerIndex  ] = tofIndex ;

		} else {
		  cout << " Malformed map file " << endl;
		  return;
		}

	}

	inData.close();


	// now copy the maps into the class array
	for ( int i = 0; i < constants::nhChannels; i++ ){
		triggerToTofMap[ i ] = westTriggerToTofMap[ i ]  ;
		triggerToTofMap[ i + constants::nhChannels ] = eastTriggerToTofMap[ i ] + constants::nhChannels;
		
		if ( 0 <= westTriggerToTofMap[ i ])
			tofToTriggerMap[ westTriggerToTofMap[ i ] ] = i;
		if ( 0 <= eastTriggerToTofMap[ i ] )
			tofToTriggerMap[ eastTriggerToTofMap[ i ] + constants::nhChannels ] = i + constants::nhChannels;

		/*
		if ( westTriggerToTofMap[ i ] != -1 )
			cout << "west " << left << setw( 4 ) << i << " -> " << westTriggerToTofMap[ i ] << endl;
		
		if ( eastTriggerToTofMap[ i ] != -1 )
			cout << "east " << left << setw( 4 ) << i << " -> " << eastTriggerToTofMap[ i ] << endl; 
		*/
	}

	cout << "[calib." << __FUNCTION__ << "] " << " completed in " << elapsed() << " seconds " << endl;


}