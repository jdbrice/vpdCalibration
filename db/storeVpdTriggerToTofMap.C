// macro to store vpd trigger to tof channel map table to database
//
// based on
//  http://www.star.bnl.gov/STAR/comp/db/StoreDbTable.cc.html
//
// Daniel Brandenburg   06/19/2014
//

#include <iostream>
#include <string>
#include <fstream>
using namespace std;

void storeVpdTriggerToTofMap( bool test = true, const char* time = "2010-01-01 00:00:00")
{


  //-- load dBase and Table definition libraries
  gSystem->Load("St_base");
  gSystem->Load("StChain");
  gSystem->Load("StUtilities");
  gSystem->Load("St_Tables.so");

  gSystem->Load("StDbLib.so");
  gSystem->Load("libStDb_Tables.so");

  //-- get the singleton manager
  StDbManager* dbManager = StDbManager::Instance();

  //-- connect to the db & get an empty container
  StDbConfigNode* configNode = dbManager->initConfig("Calibrations_tof");

  // set the request time
  string writeTime = time;

  // get the db table
  StDbTable* triggerToTofMap = configNode->addDbTable( "vpdTriggerToTofMap" );



  // maximum possible number of trigger tubes per side of vpd
  // currently there are only 16 PMTs active on the trigger side
  const short maxPMTs = 19;
  short eastTriggerToTofMap[ maxPMTs ],
        westTriggerToTofMap[ maxPMTs ];

  for(int i = 0; i < maxPMTs; i++ ) {
    eastTriggerToTofMap[ i ] = -1;
    westTriggerToTofMap[ i ] = -1;
  }

  ifstream inData;
  inData.open("vpdTriggerToTofMap.dat");

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
      cout << " Malformed map file, exiting! " << endl;
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
      cout << " Malformed map file, exiting! " << endl;
      return;
    }

  }

  inData.close();


  // report the map
  cout << " Map read in is: " << endl << endl;
  cout << "east" << endl;
  cout << "Trigger channel ---> tof channel " << endl;
  for ( int i = 0; i < maxPMTs; i++ ){
    cout << " " << (i ) << " ---> " << eastTriggerToTofMap[ i ] << endl;
  }
  cout << "west" << endl;
  cout << "Trigger channel ---> tof channel " << endl;
  for ( int i = 0; i < maxPMTs; i++ ){
    cout << " " << (i ) << " ---> " << westTriggerToTofMap[ i ] << endl;
  }

  if ( false == test ){

    vpdTriggerToTofMap_st *vpdMap = new vpdTriggerToTofMap_st[ maxPMTs * 2 ];

    for ( int i = 0; i < maxPMTs; i++ ){
      vpdMap[ i ]->eastWest = east;
      vpdMap[ i ]->triggerIndex = i;
      vpdMap[ i ]->tofIndex = eastTriggerToTofMap[ i ];
    }
    for ( int i = 0; i < maxPMTs; i++ ){
      vpdMap[ i + maxPMTs ]->eastWest = west;
      vpdMap[ i + maxPMTs ]->triggerIndex = i;
      vpdMap[ i + maxPMTs ]->tofIndex = westTriggerToTofMap[ i ];
    }

    cout << "Uploading VPD trigger to tof channel map" << endl;
    triggerToTofMap->SetTable( (char*) vpdMap, maxPMTs*2 );

    cout << "Setting store time to : " << writeTime << endl;
    dbManager->setStoreTime( writeTime.c_str() );

    dbManager->storeTable( triggerToTofMap );

    cout << " vpdTriggerToTofMap Uploaded " << endl;

  }



}