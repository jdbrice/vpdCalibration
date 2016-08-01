

#include "chainLoader.h"


void chainLoader::load( 
						TChain * chain, 	// the chain object to fill
						char* ntdir, 		// the directory in which to look for ntuples
						uint maxFiles 
						) {
	//cout << " [chainLoader] searching " << ntdir << " for ntuples" << endl;

	if (maxFiles == 0)
		maxFiles = 1000;

	uint nFiles = 0;
	DIR *dir;
	struct dirent *ent;
	bool go = true;
	if ( (dir = opendir ( ntdir ) ) != NULL) {
		
		while ( go && (ent = readdir ( dir) ) != NULL) {

	    	if ( strstr( ent->d_name, "root") ){
	    		
	    		if ( nFiles >= maxFiles ) 
	    			go = false;

	    		char fn[ 1024 ];
	    		sprintf( fn, "%s%s", ntdir, ent->d_name );
	    		cout << "[chainLoader] Adding file " << fn << " to chain" << endl;
	    		chain->Add( fn );

	    		
	    		nFiles++;

	    	}
	  	}
	  	
	  	cout << "[chainLoader] " << (nFiles - 1) << " files loaded into chain" << endl;

	  	closedir (dir);
	}

}


void chainLoader::loadList(  TChain * _chain, string _listFile, int _maxFiles ){
	
	string classname = "ChainLoader";
	cout << "( chain, listFile=" << _listFile << ", maxFiles=" << _maxFiles << " )" << endl;

	int nFiles = 0;

	string line;
	ifstream fListFile( _listFile.c_str());
	if ( fListFile.is_open() ){

		while ( getline( fListFile, line ) ){
			_chain->Add( line.c_str() );
			cout << "Adding File : " << line << endl;
			nFiles++;

			if ( _maxFiles >= 1 && nFiles >= _maxFiles ){
				break;
			}

		}
		fListFile.close();
			

	} else {
		cout << "Could not open " << _listFile  << endl;
	}

	cout << nFiles << " files loaded into chain" << endl;
} // loadList