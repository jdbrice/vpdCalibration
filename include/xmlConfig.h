#ifndef XMLCONFIG_H
#define XMLCONFIG_H 


#include "rapidXML.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cerrno>
#include <sstream>
#include <algorithm>
#include <vector>


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

	int getAsInt( char* node, int def = 0 );
	bool getAsBool( char* node, bool def = false );
	double getAsDouble( char* node, double def = 0 );
	string getAsString( char* nName, string def = "" );
	void display( char* nName );
	vector<string> getAsStringVector( char* nName );
	vector<double> getAsDoubleVector( char* nName );

	bool nodeExists( char* nName );
	bool isVector( char* nName );

private:
	string configFile;
	string fname;



};



#endif