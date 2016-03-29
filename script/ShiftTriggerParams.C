

XmlConfig config;

int eBinEdges[19][10];
int eCorrs[19][10];

int wBinEdges[19][10];
int wCorrs[19][10];
int numBins;

const Int_t startWest = 0;
const Int_t endWest = 19;
const Int_t trgEndWest = 16;

const Int_t startEast = 19;
const Int_t endEast = 38;
const Int_t trgEndEast = 35;

void writeTriggerParameters();
void readParams();

void ShiftTriggerParams( ){

	Logger::setGlobalLogLevel( "all" );
	config.loadFile( "ShiftTriggerParams.xml" );


	readParams();

	int shiftEast = config.getInt( "Shift:east", 0 );
	int shiftWest = config.getInt( "Shift:west", 0 );

	INFO( "", "EAST shift = " << shiftEast );
	INFO( "", "WEST shift = " << shiftWest );

	for ( int j = 0; j < 16; j++){
		for ( int i = 0; i < numBins; i++ ){
			eCorrs[j][i] = eCorrs[j][i] + shiftEast;
			wCorrs[j][i] = wCorrs[j][i] + shiftWest;
		}	
	}
	

	writeTriggerParameters();

}


void readParams(){

	ifstream inf( config.getString( "SlewingIn:url" ).c_str() );
	DEBUG( "readParams", "Reading Slewing Params : " << config.getString( "SlewingIn:url" ) );

	for ( int i = 0; i < 19; i ++ ){
		for ( int j = 0; j < 10; j++){
			eBinEdges[ i ][ j ] = 0;
			wBinEdges[ i ][ j ] = 0;
			eCorrs[ i ][ j ] = 0;
			wCorrs[ i ][ j ] = 0;
		}
	}


	if ( !inf.good() ){
		ERROR( "", "BAD SLEWING" );
		return;
	}

	int iEast = 0;
	int iWest = 0;

	string line;
	while( getline( inf, line ) ){

		DEBUG( "readParams", line );
		stringstream ss( line );

		string boardId;

		ss >> boardId;

		if ( "0x16" == boardId || "0x18" == boardId ){

			int channel;
			int nBins;
			int opt;

			ss >> channel >> nBins >> opt;
			DEBUG( "readParams", "nBins " << nBins );

			numBins = nBins;
			
			if ( 0 == opt ) {// bin edges
				for ( int i = 0; i < nBins; i++ ){
					int be;
					ss >> be;
					DEBUG( "readParams", "Edge [" << i << "] = " << be );
					if ( "0x16" == boardId ) // east
						eBinEdges[iEast][i] = be;
					else if ("0x18" == boardId) // west
						wBinEdges[iWest][i] = be;
				}
			} else if ( 1 == opt ) {// bin corrs
				for ( int i = 0; i < nBins; i++ ){
					int bc;
					ss >> bc;
					DEBUG( "readParams", "Corr [" << i << "] = " << bc );
					if ( "0x16" == boardId ) // east
						eCorrs[iEast][i] = bc;
					else if ("0x18" == boardId) // west
						wCorrs[iWest][i] = bc;
				}
			}
			if ( "0x16" == boardId && 1 == opt ) // east
				iEast ++;
			else if ("0x18" == boardId && 1 == opt ) // west
				iWest ++;	
		}

	}

	inf.close();
}


void writeTriggerParameters(){

	cout << "Writing out Trigger Parameters" << endl;

	string outName = config.getString( "baseName", "" ) + config.getString( "paramsOutput", "params.dat" );

	ofstream f;
	f.open( outName.c_str() );

	/**
	 * Preamble
	 */
	string pa = 
	"# Board	Channel	NumBins	BinOrOffset ADCBinLimit0 ADCBinLimit1 ... \n"
	"# Board	Channel	NumBins	BinOrOffset TACOffsetBin0 TACOffsetBin1 ... \n"
	"#   Board is 0xYY where YY is top two hex digits of board address \n"
	"#   Channel is 0-31; only TAC channels are valid (4,5,6,7,12,13,14,15,20,21,22,23,28,29,30,31) \n"
	"#   BinOrOffset = 0 for a Bin Limit row, 1 for an Offset row \n"
	"#   NumBins must be exactly equal to the number of ADC Bins implemented in the QT VHDL \n"
	"#     All ADC Bin Limits must be defined \n"
	"#     Each ADC Bin Limit must be equal to or greater than the previous bin limit \n"
	"#     An ADC value is in bin 'N' if BinLimit[N-1] < ADC <= BinLimit[N] \n"
	"#     The first bin has an implied lower limit of '0' \n"
	"#     An ADC value of '0' falls into the first bin \n"
	"#     Unused ADC bins should have a limit of 4095 and be the highest ADC bins \n"
	"#     At least one bin limit should be 4095 (full range covered) \n"
	"# Slew Corrections come after the QT LUT (ie after TAC Offset/ADC Pedestal Subtraction) \n"
	"# Slew Corrections come before QT Channel Masks";

	f << pa << endl;
	f << "\n" << "#VP001 East" << endl;

	int qtChannels[] = { 4,5,6,7,12,13,14,15,20,21,22,23,28,29,30,31 };
	int iQT = 0;

	// Write out the East VPD
	for ( int j = 0; j < 16; j++ ){
		// if ( deadDetector[ j ] ) continue;

		int channel = 16;

		int board = qtChannels[ iQT ];

		f << "0x" << channel << "\t" << board << setw( 7 ) << numBins << " 0"; // zero for bin edges

		for ( int i = 0; i < numBins; i++ ){
			f << setw(7) << eBinEdges[ j ][ i ] << " ";
		}
		f << endl;

		f << "0x" << channel << "\t" << board << setw( 7 ) << numBins << " 1"; // 1 for bin corrections
		for ( int i = 0; i < numBins; i++ ){
			int full_cor = eCorrs[ j ][ i ]; 
			f << setw(7) << full_cor << " ";
		}
		f << endl;

		iQT++;
	}


	f << "\n" << "#VP002 West" << endl;

	iQT = 0;
	// Write out the East VPD
	for ( int j = 0; j < 16; j++ ){
		// if ( deadDetector[ j ] ) continue;

		int channel = 18;

		int board = qtChannels[ iQT ];

		f << "0x" << channel << "\t" << board << setw( 7 ) << numBins << " 0"; // zero for bin edges

		for ( int i = 0; i < numBins; i++ ){
			f << setw(7) << wBinEdges[ j ][ i ] << " ";
		}
		f << endl;

		f << "0x" << channel << "\t" << board << setw( 7 ) << numBins << " 1"; // 1 for bin corrections
		for ( int i = 0; i < numBins; i++ ){
			// int full_cor = TMath::Nint( -1 * (correction[ j ][ i ] + this->initialOffsets[ j ] ) / tacToNS);
			int full_cor = wCorrs[ j ][ i ];
			f << setw(7) << full_cor << " ";
		}
		f << endl;

		iQT++;
	}

	f.close();

}