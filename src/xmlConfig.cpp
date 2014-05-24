#include "xmlConfig.h"
#include "rapidXML.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cerrno>
#include <sstream>
#include <stdio.h>
#include <string.h>
using namespace std;
using namespace rapidxml;

/*
*	Constructor
*	Parameters:
*		filename: the config filename
*/
xmlConfig::xmlConfig( const char * filename ) {
	//cout << "[xmlConfig] opening " << filename << endl;

	configFile = getFileContents( filename );
	//cout << "[xmlConfig] : " << endl << configFile << endl;
}



string xmlConfig::getFileContents(const char *filename)
{
  
  ifstream in(filename, ios::in | ios::binary);

  if (in)
  {
    std::string contents;

    in.seekg(0, std::ios::end);
    contents.resize(in.tellg());

    in.seekg(0, std::ios::beg);
    in.read(&contents[0], contents.size());

    in.close();

    return(contents);
  }
  throw(errno);
}

int xmlConfig::getAsInt( char* nName ) {


	string str = getAsString( nName );
	if ( str.length() >= 1 ){
		return atoi( str.c_str() );
	}

	return 0;

}

bool xmlConfig::getAsBool( char* nName ) {


	string str = getAsString( nName );
	if ( str.length() >= 1 ){
		return ( atoi( str.c_str() ) >= 1 );
	}
	return false;
}

string xmlConfig::getAsString( char* nName ) {

	stringstream sstr;

	xml_document<> doc;
	char* cstr = new char[configFile.size() + 1];  	// Create char buffer to store string copy
  	strcpy (cstr, configFile.c_str());             		// Copy string into char buffer

  	doc.parse<0>(cstr); 

  	if ( doc.first_node() ){

	  	//cout << "Name of my first param is: " << doc.first_node()->first_node()->name() << endl;

		xml_node<> *node = doc.first_node()->first_node( nName );

		if ( node ){
			//cout << "Node: " << nName << " = " << node->value() << endl;
			
			sstr << node->value();
			return sstr.str();
		}
		else {
			//cout << "[xmlConfig] Node " << nName << " DNE" << endl;
		}

		return sstr.str();
	}

	return sstr.str();
}

double xmlConfig::getAsDouble( char* nName ) {

	string str = getAsString( nName );
	if ( str.length() >= 1 ){
		return atof( str.c_str() );
	}

	return 0;
}
