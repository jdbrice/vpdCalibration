
#include <iostream>
#include "allroot.h"
#include "constants.h"
#include "histoBook.h"
#include "chainLoader.h"
#include "calib.h"



int main( int argc, char* argv[] ) {

    cout << " VPD Start Time Calibration" << endl;
    
    if ( argc < 2 ){
    	
    	cout << "Call executable with \n "
    			"1) xml config filename"
                << endl;

       	return 0;

    } 

    /* Give a summary of config file */
    xmlConfig config( argv[ 1 ] );
    cout << "[CONFIG] File: " << argv[ 1 ] << endl;
    cout << "[CONFIG]" << " jobType is "        << config.getAsString( "jobType" ) << endl;
    cout << "[CONFIG]" << " rootOutput = "      << config.getAsString( "rootOutput" ) << endl;
    cout << "[CONFIG]" << " dataDir = "         << config.getAsString( "dataDir" ) << endl;
    cout << "[CONFIG]" << " variableBinning = " << config.getAsString( "variableBinning" ) << endl;
    cout << "[CONFIG]" << " maxFiles = "        << config.getAsString( "maxFiles" ) << endl;
    cout << "[CONFIG]" << " numIterations = "   << config.getAsString( "numIterations" ) << endl;
    cout << "[CONFIG]" << " paramsOutput = "    << config.getAsString( "paramsOutput" ) << endl;
    cout << "[CONFIG]" << " paramsInput = "     << config.getAsString( "paramsInput" ) << endl;
    cout << "[CONFIG]" << " removeOffset = "    << config.getAsString( "removeOffset" ) << endl;
    cout << "[CONFIG]" << " outlierRejection = "<< config.getAsString( "outlierRejection" ) << endl;
    cout << "[CONFIG]" << " numTOTBins = "      << config.getAsString( "numTOTBins" ) << endl;

    string jobType = config.getAsString( "jobType" );

    /*   Load the files into the chain */
    TChain * chain = new TChain( "tof" );
    if ( config.getAsInt( "maxFiles" ) ){
        chainLoader::load( chain, (char*)config.getAsString( "dataDir" ).c_str(), config.getAsInt( "maxFiles" ) );
    }
    
    // get the num of iterations
    int numIterations = 5;
    if ( config.getAsInt( "numIterations" ) >= 1)
        numIterations = config.getAsInt( "numIterations" );

    // create a calibration object
    calib vpdCalib( chain, numIterations, config );
 

    if ( jobType == (string)"genReport" ){

        vpdCalib.readParameters( config.getAsString( "paramsInput" ) );
    
        vpdCalib.zVtxPairs();

    }
    if ( jobType == (string)"calibrate" ){

        // determine the variable binning in tot space
        vpdCalib.binTOT( config.getAsBool( "variableBinning")  ); 

        // get the inital offsets
        if ( config.getAsBool( "removeOffset" ) )
            vpdCalib.offsets();

        // look at the VPD vs TPC vertex pairs before calibration
        vpdCalib.zVtxPairs();

        // run the main calibration loop
        vpdCalib.loop();

        // look at VPD vs TPC vertex after calibration
        vpdCalib.zVtxPairs();
        
        // write out the parameters file
        vpdCalib.writeParameters( config.getAsString( "paramsOutput" ) );
    }


	return 0;
}
