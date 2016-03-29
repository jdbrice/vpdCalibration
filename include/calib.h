#ifndef CALIB_H
#define CALIB_H

#include "allroot.h"

#include "histoBook.h"
#include "constants.h"
#include "TOFrPicoDst.h"
#include "splineMaker.h"
#include <vector>
#include <map>

// clock_t, clock, CLOCKS_PER_SEC 
#include <time.h>       

// for testing if stdout is interactive or pipe / file
#include "unistd.h"
#include "xmlConfig.h"
#include "utils.h"
#include "reporter.h"


class calib{

private:

	// the canvas used to draw report hidtos
	reporter* report;

	// the total number of iterations to try
	uint maxIterations;

	// the current iteration in the calibration loop
	uint currentIteration;

	// the channel used as the reference for calculating offsets
	uint refChannel;
	// Dont use bad detectors
	// defaults to false, set to true if dead channel is found
	bool deadDetector[ constants::nChannels ];

	// the initial offsets for each channel relative to the 1st channel on the west side
	double tacOffsets[ constants::nChannels ];
	double initialOffsets[ constants::nChannels ];
	double outlierOffsets[ constants::nChannels ];
	// the west - east offset only needed for the case of not removing initial offsets
	double eastWestOffset;
	// the floating channel 1 offset
	double finalWestOffset;

	// number of tot bins to use
	// defaults to constants::numTOTBins if not set in config
	int numTOTBins;

	// the minimum and maximum tot range to calibrate
	double minTOT, maxTOT;

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

	// Interpolation based corrections
	splineMaker * spline[ constants::nChannels ];
	Interpolation::Type splineType;
	bool useSpline;


	// list of detectors with prompt hits for this event and usable in calibration
	// calculated for each event in outlierRejection()
	bool useDetector[ constants::nChannels ];

	bool westIsGood;
	bool eastIsGood;


	// use for timing
	clock_t startTime;

	// config file
	xmlConfig config;

	vector<double> avgNTimingCut;
	vector<double> vzOutlierCut;
	double avgNBackgroundCut;

	// the variables to use for calibration
	// tof-le / tof-tot
	// bbq-tdc / bbq-adc
	// mxq-tdc / mxq-adc
	// or combination
	string xVariable;
	string yVariable;

	// labels built from the x/y variable choice
	string xLabel, yLabel;


	// store the trigger to tof channel map
	int triggerToTofMap[ constants::nChannels ];
	int tofToTriggerMap[ constants::nChannels ];
	bool mapTriggerToTof;
	bool convertTacToNS;
	double TACToNS;

	vector<int> maskedChannels;
	map<int, bool> channelMask;
	int firstRun, lastRun;
	

public:


	// Constructor
	calib( TChain * chain, uint nIterations, xmlConfig config );

	// destructor
	~calib();

	// finish the calibration by conducting fits etc.

	void finish();

	// calculates the inital offsets of each channel
	void offsets( );
	void updateOffsets();
	void finalOffsets( );

	void getInitialOffsets();

	void hardCodeTACOffsets();

	// determines the binning in tot space for each channel
	void binTOT( bool variableBinning = true );

	// executes the full correction loop
	void loop();

	// executes a single loop of the iterative correction process
	void step( );
	void checkStep( );
	void prepareStepHistograms();

	// after everything calculate the reference offset on channel 1 on the west
	void referenceOffset();

	// get the correction for a given channel, given tot value
	double getCorrection( int vpdChannel, double tot );
	// get the bin for a given tot value ona given channel
	int binForTOT( int vpdChannel, double tot );

	void writeParameters(  );
	void writeTriggerParameters( );
	void readParameters( );
	
	// reports
	void stepReport();

	// used for generalizing the slewing correction to the trigger side signals.
	// Use the configuration file to set the variables to use
	// Default are the TOF-side electronics Leading edge time (TDC) and time-over-threshold (TOT)
	double getX( int channel );	// default is TOT
	double getY( int channel );	// default is TDC

	bool doingTrigger() {
		if ( "bbq-tdc" == yVariable || "mxq-tdc" == yVariable )
			return true;
		return false;
	}
	

protected:

	bool runInRange( int run ){
		if ( 0 >= firstRun || 0 >= lastRun )
			return true;
		
		if ( run < firstRun || run > lastRun )
			return false;
		return true;
	}

	void makeCorrections();

	// performs outlier rejection by selecting detectors on the east and west only when they produce
	// a z vertex that is consistent with a prompt particle ( ie consistent with TPC vertex ).
	void outlierRejection( bool reject = true );

	void averageN();

	void readTriggerToTofMap();

	static Double_t detectorResolution(Double_t *x, Double_t *par);

	/*
	*	Utility functions that should be moved soon
	*/ 
	void startTimer( ) { startTime = clock(); }
	double elapsed( ) { return ( (clock() - startTime) / (double)CLOCKS_PER_SEC ); }
};



#endif