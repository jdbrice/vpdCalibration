#ifndef UTILS_H
#define UTILS_H

#include <string>

using namespace std;

class utils
{
public:

	static string toString( double );
	static string toString( int );
	static string toString( uint );

	utils();
	~utils();

	/* data */
};


namespace jdbUtils{

	std::string ts( int );
	std::string ts( double );
	std::string ts( uint );
}




#endif