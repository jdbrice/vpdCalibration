
#ifndef CHAINLOADER_H
#define CHAINLOADER_H

#include "dirent.h"
#include "allroot.h"

class chainLoader{

public:
	static void load( TChain * chain, char* ntdir, uint maxFiles = 1000 );
};

#endif

