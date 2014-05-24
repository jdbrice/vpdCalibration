

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
	    		chain->Add( fn );

	    		//cout << "[chainLoader] Adding file " << fn << " to chain" << endl;
	    		nFiles++;

	    	}
	  	}
	  	
	  	cout << "[chainLoader] " << (nFiles - 1) << " files loaded into chain" << endl;

	  	closedir (dir);
	}

}