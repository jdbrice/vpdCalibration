// macro to read vpd trigger to tof channel map table from database
//
// based on
//  http://www.star.bnl.gov/STAR/comp/db/StoreDbTable.cc.html
//
// Daniel Brandenburg   06/19/2014
//

#include <iostream>
#include <fstream>
using namespace std;

void readVpdTriggerToTofMap(const char* time = "2010-01-01 00:00:00")
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
  string readTime = time;
  dbManager->setRequestTime( readTime.c_str());

  // get the db table
  StDbTable* triggerToTofMap = configNode->addDbTable( "vpdTriggerToTofMap" );
  dbManager->fetchDbTable( triggerToTofMap );

  // output some table details
  cout<<vpdTriggerToTofMap->getVersion()<<endl;
  cout<<vpdTriggerToTofMap->getBeginDateTime()<<endl;
  cout<<vpdTriggerToTofMap->getEndDateTime()<<endl;


  vpdTriggerToTofMap_st* table = static_cast<vpdTriggerToTofMap_st*>(vpdTriggerToTofMap->GetTable());
  
  if( !table ) {
    cout << " Table is invalid, exiting! " << endl;
    return;
  }

  Int_t nRows = vpdTriggerToTofMap->GetNRows();
  
  cout << " Number of rows = " << nRows << endl;

  // maximum possible number of trigger tubes per side of vpd
  // currently there are only 16 PMTs active on the trigger side
  const short maxPMTs = 19;
  short   eastTriggerToTofMap[ maxPMTs ],
          westTriggerToTofMap[ maxPMTs ];

  for(int i = 0; i < maxPMTs; i++ ) {
    eastTriggerToTofMap[ i ] = -1;
    westTriggerToTofMap[ i ] = -1;
  }

  cout<<"<---------------- Read out from DataBase -------------->"<<endl;


  ofstream outData;
  outData.open("vpdTriggerToTofMap_read.dat");
  for ( int i = 0; i < nRows; i++ ){
    

    // ensure the trigger index is valid
    if ( table[ i ].triggerIndex >= 0 ){

      if ( east == table[ i ].eastWest )
        eastTriggerToTofMap[ table[ i ].triggerIndex ] = table[ i ].tofIndex;
      else if ( west == table[ i ].eastWest )
        westTriggerToTofMap[ table[ i ].triggerIndex ] = table[ i ].tofIndex;

    }

  }

  outData << "east" << endl;
  outData << maxPMTs << endl;
  for ( int i = 0; i < maxPMTs; i++ ){
    outData << i << " " << eastTriggerToTofMap[ i ] << endl;
  }

  outData << "west" << endl;
  outData << maxPMTs << endl;
  for ( int i = 0; i < maxPMTs; i++ ){
    outData << i << " " << westTriggerToTofMap[ i ] << endl;
  }
  
  outData.close();

  cout << " VPD Trigger to Tof channel map written to vpdTriggerToTofMap_read.dat" << endl;


}