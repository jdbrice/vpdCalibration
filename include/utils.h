#ifndef UTILS_H
#define UTILS_H

#include <string>

using namespace std;


namespace jdbUtils{

	std::string ts( int );
	std::string ts( double );
	std::string ts( uint );
	void progressBar( int i, int nevents, int max );
}




#endif