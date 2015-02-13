/*******************************************************************
*
*	Name:	globalOff.C
*	Description:
*		Loads in the histogram BTofTimeRes_vs_Tray
*		fits the FitSlicesY to a gaussian and 
*		applies the mean as the new offset for each tray
*
*
********************************************************************/


#include <iostream>
#include <fstream>
using namespace std;

static const int TRAYS = 120;
static const int MODULES = 32;
static const int BOARDS = 6;

static const int NUMBER_OF_OFFSETS = 120;

struct t0Offset {

	public:
		double t0;
		int firstTray;	// inclusive
		int lastTray;	// inclusive
		
		int firstModule = 1; // inclusive
		int lastModule = 32; // inclusive

};

t0Offset newOffsets[ NUMBER_OF_OFFSETS ];


double applyOffsets( int tray, int module, double timing ){

	for ( int i = 0; i < NUMBER_OF_OFFSETS; i ++ ){

		if ( tray >= newOffsets[i].firstTray && tray <= newOffsets[i].lastTray ){
			if ( module >= newOffsets[i].firstModule && module <= newOffsets[i].lastModule ){

				return timing + newOffsets[i].t0;
			}
		}

	}

	return timing;
}


void globalOff( string inName = "t0_4DB.dat", string outName = "meanByMean_t0_4DB.dat", double gOffset = 0.0 ) {

	cout << "[INPUT] " << inName << endl;
	cout << "[OUTPUT] " << outName << endl;


	// initialize the offset structures to cover all moduels
	for ( int i = 0; i < NUMBER_OF_OFFSETS; i++ ){
		newOffsets[ i ].firstModule = 1;
		newOffsets[ i ].lastModule = 32;
	}

	// loop through all of the trays and get an offset for each
	for ( int ti = 0; ti < NUMBER_OF_OFFSETS ; ti++){
		int tray = ti+1;
		// build an offset for all of the normal trays
		if ( true ){
			newOffsets[ ti ].t0 = gOffset;
			newOffsets[ ti ].firstTray = tray;
			newOffsets[ ti ].lastTray = tray;
		} else {
			// OFFSETS FOR SPECIAL TRAYS APPLIED IN fitDoubles.C
		}

	}





	double t0[ TRAYS ][ MODULES ][ BOARDS ];

	ifstream infile;
    infile.open(inName.c_str());
	if ( infile.is_open() ){

		while ( !infile.eof()){
			int tray, board, module;
			infile >> tray >> module >> board;
			double timing;
			infile >> timing;
			
			if ( tray <= TRAYS && module <= MODULES && board <= BOARDS ){
				//cout << "[ " << tray << " ] " << "[ " << module << " ] " << "[ " << board << " ]  = " << timing << endl;
				t0[ tray-1 ][ module-1 ][board-1] = timing;
			}

		}

	} else {
		cout << "Could not open file \n";
	}


	// now apply the offsets and write out a new file:
	ofstream out( outName.c_str());
	for ( int it = 0; it < TRAYS; it ++ ){
		for ( int im = 0; im < MODULES; im ++ ){
			for ( int ib = 0; ib < BOARDS; ib ++ ){
				
				double ct = t0[it][im][ib];

				ct = applyOffsets( it+1, im+1, ct );

				out << (it+1) << " " << (im+1) << " " << (ib+1) << endl;
				out << ct << endl;

			}		
		}
	}

	out.close();

	return;
}

