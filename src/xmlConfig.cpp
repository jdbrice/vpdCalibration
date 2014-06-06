#include "xmlConfig.h"
#include "rapidXML.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cerrno>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <iomanip>
using namespace std;
using namespace rapidxml;

/*
*	Constructor
*	Parameters:
*		filename: the config filename
*/
xmlConfig::xmlConfig( const char * filename ) {

	fname = filename;
	configFile = getFileContents( filename );

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

int xmlConfig::getAsInt( char* nName, int def ) {

	if ( !nodeExists( nName ) )
		return def;

	string str = getAsString( nName );
	if ( str.length() >= 1 ){
		return atoi( str.c_str() );
	}

	return def;

}

bool xmlConfig::getAsBool( char* nName, bool def ) {

	if ( !nodeExists( nName ) )
		return def;

	string str = getAsString( nName );

	// first check for string literal "true" or "false"
	// first push to lower case
	transform(str.begin(), str.end(), str.begin(), ::tolower);
	if ( str == (string)"false" )
		return false;
	if ( str == (string)"true" )
		return true;

	// if not look for an integer
	// 0 or negative = false
	// any positive = true
	if ( str.length() >= 1 ){
		return ( atoi( str.c_str() ) >= 1 );
	}
	return def;
}

string xmlConfig::getAsString( char* nName, string def ) {

	if ( !nodeExists( nName ) )
		return def;

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

		return def;
	}

	return def;
}

vector<string> xmlConfig::getAsStringVector( char* nName ) {

	stringstream sstr;

	std::vector<string> res;

	xml_document<> doc;
	char* cstr = new char[configFile.size() + 1];  	// Create char buffer to store string copy
  	strcpy (cstr, configFile.c_str());             		// Copy string into char buffer

  	doc.parse<0>(cstr); 

  	if ( doc.first_node() ){

	  	//cout << "Name of my first param is: " << doc.first_node()->first_node()->name() << endl;

		xml_node<> *node = doc.first_node()->first_node( nName );
		if ( !node )
			return res;
		for (xml_node<> *child = node->first_node(  ); child; child = child->next_sibling() )
		{
		    sstr.str("");
		    sstr << child->value();
		    res.push_back( sstr.str() );
		}
		
	}

	return res;
}


std::vector<double> xmlConfig::getAsDoubleVector( char* nName ){


	vector<string> str = getAsStringVector( nName );
	vector<double> res;
	for ( int i = 0; i < str.size(); i++ ){
		res.push_back( atof( str[ i ].c_str() ) );
	}
	return res;

}

double xmlConfig::getAsDouble( char* nName, double def ) {

	if ( !nodeExists( nName ) )
		return def;

	string str = getAsString( nName );
	if ( str.length() >= 1 ){
		return atof( str.c_str() );
	}

	return def;
}

void xmlConfig::display( char* nName ){
	// if it is a vector print it that way
	stringstream sstr;
	std::vector<string> v = getAsStringVector( nName );
	if ( v.size() >= 1 && v[ 0 ] != getAsString( nName ) ){
		
		for ( int i = 0; i < v.size(); i++ ){
			sstr.str("");
			sstr << nName << "[ " << i << " ] ";
			cout << "[" << fname <<  "] " << std::left <<  setw(20) << sstr.str() << v[i] << endl;
		}
		return;
	}

	cout << "[" << fname <<  "] " << std::left << setw(20) <<  nName  << getAsString( nName ) << endl;
}

bool xmlConfig::nodeExists( char* nName ){
	

	xml_document<> doc;
	char* cstr = new char[configFile.size() + 1];  	// Create char buffer to store string copy
  	strcpy (cstr, configFile.c_str());             		// Copy string into char buffer

  	doc.parse<0>(cstr); 

  	if ( doc.first_node() ){

		xml_node<> *node = doc.first_node()->first_node( nName );
		if ( !node )
			return false;
		else 
			return true;
				
	}

	return false;
}