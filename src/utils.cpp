
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