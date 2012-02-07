#include "Utilities.h"

/** \file Utilities.cc Contains often used functions for vectors,  matrices and other. */

// C++ includes
#include <assert.h>
#include <iostream>

// ROOT includes
#include <TMath.h>
#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TClass.h>

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

TObject * get_object(const char * filename, const char * objectname)
{
  // try to find out if file is already opened
  TFile * f = (TFile *) gROOT->GetListOfFiles()->FindObject(filename);
  if (f != 0) {
//     cout << "File " << filename << " was already open" << endl;
  }
  else {
    f = new TFile(filename, "READ");
    if (f == 0) {
      ERROR("Could not create TFile object for " << filename);
      return 0;
    }
    if (!f->IsOpen()) {
      ERROR("Could not open file " << filename);
      return 0;
    }
  }
  TObject * obj = f->Get(objectname);
  if (obj == 0) {
    ERROR("Could not get object " << objectname << " from file " << filename);
    return 0;
  }
  // put histograms in main directory
  if (obj->IsA()->InheritsFrom("TH1")) {
    TH1 * h = (TH1 *) obj;
    h->SetDirectory(0);
  }
  delete f;
  return obj;
}

