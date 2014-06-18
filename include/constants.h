#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "TROOT.h"
// constants should go in here

class constants{

private:

public:

	static const Int_t nChannels = 38; 
	static const Int_t numTOTBins = 40;
	static const Int_t minTOT = 10;
	static const Int_t maxTOT = 40;
	static const Int_t nhChannels = 19;

	static const Int_t startWest = 0;
	static const Int_t endWest = 19;

	static const Int_t startEast = 19;
	static const Int_t endEast = 38;

	// minimum number of hits 
	static const Int_t minHits = 3;

	static const Double_t c = 29.9792458;	// 29.979 * 10^7 m/s

	static const Double_t tacToNS = 0.01773; // from bill llope


	
};







#endif