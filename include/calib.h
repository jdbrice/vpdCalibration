#ifndef CALIB_H
#define CALIB_H

#include "allroot.h"

#include "histoBook.h"
#include "constants.h"
#include "TOFrPicoDst.h"
#include <vector>

// clock_t, clock, CLOCKS_PER_SEC 
#include <time.h>       

// for testing if stdout is interactive or pipe / file
#include "unistd.h"
#include "xmlConfig.h"


class calib{

private:

	// the total number of iterations to try
	uint maxIterations;

	// the current iteration in the calibration loop
	uint currentIteration;

	// Dont use bad detectors
	// defaults to false, set to true if dead channel is found
	bool deadDetector[ constants::nChannels ];

	// the initial offsets for each channel relative to the 1st channel on the west side
	double initialOffsets[ constants::nChannels ];
	// the west - east offset
	double westMinusEast;

	// number of tot bins to use
	// defaults to constants::numTOTBins if not set in config
	int numTOTBins;

	// the corrections for each channel as a function of ToT
	double * correction[ constants::nChannels ];

	// the main chain object
	TChain * _chain;

	// the histobook that stores all of our calibration histograms
	histoBook *book;

	// the pico dst for simpler chain usage
	TOFrPicoDst * pico;

	// variable bins for tot values -> helps with low statistics
	// calculated in binTOT
	Double_t * totBins[ constants::nChannels ];


	// outlier rejection
	// for each event we will want 
	// # of pairs east vs west

	// list of detectors with prompt hits for this event and usable in calibration
	// calculated for each event in outlierRejection()
	bool useDetector[ constants::nChannels ];


	// use for timing
	clock_t startTime;

	// config file
	xmlConfig config;
	

public:

	// Constructor
	calib( TChain * chain, uint nIterations, xmlConfig config );

	// destructor
	~calib();

	// calculates the inital offsets of each channel
	void offsets( );

	// determines the binning in tot space for each channel
	void binTOT( bool variableBinning = true );

	// calculate the zVtx for each pair of east / west times
	void zVtxPairs();

	// executes the full correction loop
	void loop();

	// executes a single loop of the iterative correction process
	void step( );
	void prepareStepHistograms();

	// get the correction for a given channel, given tot value
	double getCorrection( int vpdChannel, double tot );
	// get the bin for a given tot value ona given channel
	int binForTOT( int vpdChannel, double tot );

	void writeParameters(  );
	void readParameters( );
	
	// reports
	void stepReport();

protected:

	void makeCorrections();

	// performs outlier rejection by selecting detectors on the east and west only when they produce
	// a z vertex that is consistent with a prompt particle ( ie consistent with TPC vertex ).
	void outlierRejection( bool reject = true );

	/*
	*	Utility functions that should be moved soon
	*/ 
	void startTimer( ) { startTime = clock(); }
	double elapsed( ) { return ( (clock() - startTime) / (double)CLOCKS_PER_SEC ); }
	void progressBar( double progress, int max ){
		
		// skip for non interactive output
		if (!isatty(fileno(stdout)) && progress <= 1 )
			return;

		double per = progress  * 100;
		per = TMath::Nint( per );

		cout << "[";
    	for ( int ip = 0; ip < max; ip ++ ){
    		if ( ip < TMath::Nint( (progress * (double)max) ) )
    			cout << "=";
    		else 
    			cout << " ";
    	}
    	if (isatty(fileno(stdout)) ){ 
	 	   	cout << "]" << per << "%" << "\r";
			std::cout.flush();
			if (progress > 1)
				cout << "[" << endl;
		} else {
				cout << "]" << per << "%" << "\n";
		}
		
	}
};



#endif