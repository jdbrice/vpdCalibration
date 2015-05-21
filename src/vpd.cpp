
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
    cout << endl;
    config.display( "baseName" );
    config.display( "rootOutput" );
    config.display( "reportOutput" );
    config.display( "paramsOutput" );
    config.display( "paramsInput" );
    config.display( "paramsLegend" );
    cout << endl;
    config.display( "channelMap" );
    config.display( "mapTriggerToTof" );
    cout << endl;
    config.display( "xVariable" );
    config.display( "yVariable" );
    cout << endl;
    config.display( "dataDir" );
    config.display( "maxFiles" );
    cout << endl;
    config.display( "numIterations" );
    config.display( "variableBinning" );
    config.display( "binMaxError" );
    config.display( "removeOffset" );
    config.display( "outlierRejection" );
    config.display( "numTOTBins" );
    config.display( "minTOT" );
    config.display( "maxTOT" );
    config.display( "splineType" );
    cout << endl;
    config.display( "vzOutlierCut" );    
    cout << endl;
    config.display( "avgNBackgroundCut" );
    config.display( "avgNTimingCut" );
    /* Give a summary of config file */

    cout << endl << endl << "Beginning Calibration" << endl << endl;

    string jobType = config.getAsString( "jobType", "calibration" );

    // Load the files into the chain 
    TChain * chain = new TChain( "tof" );
   
    chainLoader::load( chain, (char*)config.getAsString( "dataDir" ).c_str(), config.getAsInt( "maxFiles", 10000 ) );
    
    // get the num of iterations
    int numIterations = config.getAsInt( "numIterations", 5 );

    // create a calibration object
    calib vpdCalib( chain, numIterations, config );
 

    if ( (string)"paramReport" == jobType  ){

        vpdCalib.readParameters(  );

    } else if ( (string)"checkParams" == jobType  ){

        vpdCalib.readParameters(  );

        vpdCalib.checkStep();
        vpdCalib.step();

        vpdCalib.writeParameters(  );
        
        

    } else if ( (string)"calibrate" == jobType ){

        // determine the variable binning in tot space
        vpdCalib.binTOT( config.getAsBool( "variableBinning")  );     
        
        // get the inital offsets
        vpdCalib.offsets();

        
        // run the main calibration loop
        vpdCalib.loop();

        //vpdCalib.checkStep();

        // write out the parameters file
        vpdCalib.writeParameters();
        
    }


	return 0;
}
