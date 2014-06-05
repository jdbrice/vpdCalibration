
#include "utils.h"


string utils::toString( double d ){
	return to_string( (long double) d);
}
string utils::toString( int d ){
	return to_string( (long long int) d);
}
string utils::toString( uint d ){
	return to_string( (long long uint) d);
}



namespace jdbUtils{

	std::string ts( int i ){
		return utils::toString( i );
	}
	std::string ts( double d ){
		return utils::toString( d );
	}
	std::string ts( uint u ){
		return utils::toString( u );
	}

}