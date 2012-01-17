#include "Utilities.h"

/** \file Utilities.cc Contains often used functions for vectors,  matrices and other. */

// C++ includes
#include <assert.h>
#include <iostream>

// ROOT includes
#include <TMath.h>

using namespace std;

/** \file Utilities.h
 *
 * \brief Various preprocessor macros, variables, functions, and procedures
 * that misfit elsewhere.
 *
 */

// global logging level
int gLogLevel = 3;

/** Test if two values are equal within a certain precision. */
bool equal(const double val1, const double val2, const double precision)
{
  if (val1 == 0)
    return val2 < precision;
  if (val2 == 0)
    return val1 < precision;
//   cout << "val1= " << val1 << "   val2= " << val2 << "  equality= " << TMath::Abs(2*(val1-val2)/(val1+val2)) << endl;
  return TMath::Abs(2*(val1-val2)/(val1+val2)) < precision;
}

// split a string
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

// split a string
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}
