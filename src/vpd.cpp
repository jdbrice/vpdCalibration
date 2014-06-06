
#include <iostream>
#include "allroot.h"
#include "constants.h"
#include "histoBook.h"
#include "chainLoader.h"
#include "calib.h"
#include "utils.h"



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

    config.display( "jobType" );
    config.display( "baseName" );
    config.display( "rootOutput" );
    config.display( "reportOutput" );
    config.display( "paramsOutput" );
    config.display( "paramsInput" );
    cout << endl;
    config.display( "dataDir" );
    config.display( "maxFiles" );
    cout << endl;
    config.display( "numIterations" );
    config.display( "variableBinning" );
    config.display( "removeOffset" );
    config.display( "outlierRejection" );
    config.display( "splineType" );
    config.display( "numTOTBins" );
    config.display( "splineType" );
    cout << endl;
    config.display( "vzOutlierCut" );    
    cout << endl;
    config.display( "avgNBackgroundCut" );
    config.display( "avgNTimingCut" );

    return 0;
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

        vpdCalib.readParameters(  );
    
        //vpdCalib.zVtxPairs();

    }
    if ( jobType == (string)"calibrate" ){

        // determine the variable binning in tot space
        vpdCalib.binTOT( config.getAsBool( "variableBinning")  );     
        
        // get the inital offsets
        vpdCalib.offsets();

        // look at the VPD vs TPC vertex pairs before calibration
        //if ( config.getAsBool( "compareVpdTPC") )
        //   vpdCalib.zVtxPairs();

        // run the main calibration loop
        vpdCalib.loop();

        // look at VPD vs TPC vertex after calibration
        //if ( config.getAsBool( "compareVpdTPC") )
        //    vpdCalib.zVtxPairs();
      
        // write out the parameters file
        //vpdCalib.writeParameters(  );
        
    }


	return 0;
}
