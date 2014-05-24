#ifndef XMLCONFIG_H
#define XMLCONFIG_H 


#include "rapidXML.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cerrno>
#include <sstream>

using namespace std;
using namespace rapidxml;

class xmlConfig {

public: 
	/*
	*	Constructor
	*	Parameters:
	*		filename: the config filename
	*/
	xmlConfig( const char * filename = "config.xml" );

	string getFileContents(const char *filename);

	int getAsInt( char* node );
	bool getAsBool( char* node );
	double getAsDouble( char* node );
	string getAsString( char* nName );
	//void printConfig();

private:
	string configFile;

};



#endif