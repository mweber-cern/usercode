#ifndef _Utilities_h
#define _Utilities_h

// C++ / STL includes
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <assert.h>

// ROOT includes
#include "TString.h" // for char *Form(...)
#include "TObject.h"

/** \file Utilities.h
 *
 * \brief Various preprocessor macros, variables, functions, and procedures
 * that misfit elsewhere.
 *
 */

// global Variable to log errors (1), warnings (2), info (3), debug(4,5,...)
extern int gLogLevel;

#ifndef NDEBUG
#define LOG(level, message) { if (gLogLevel >= level) { switch (level) { \
case 1: std::cerr << "ERROR: " << message << std::endl; break; \
case 2: std::cerr << "WARNING: " << message << std::endl; break; \
case 3: std::cout << "INFO: " << message << std::endl; break; \
default: std::cout << "DEBUG: " << message << std::endl; } } }
#else
#define LOG(level, message) ;
#endif

#define ERROR(message) LOG(1, message);
#define WARNING(message) LOG(2, message);
#define INFO(message) LOG(3, message);
#define DEBUG(message) LOG(100, message);

// throw an exception that tells me where the exception happened
#define THROW(errmsg) throw (std::string( __PRETTY_FUNCTION__ )+std::string(" (file: ")+std::string( __FILE__ )+std::string(", line: ")+std::string( Form("%d", __LINE__) )+std::string(") ")+std::string(errmsg));
#define CATCH catch (std::string message) { cerr << "EXCEPTION in " << message << std::endl; }
  
// test if two values are equal within given precision
bool equal(const double val1, const double val2, const double precision);
TObject * get_object(const char * filename, const char * objectname);

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);


#endif
