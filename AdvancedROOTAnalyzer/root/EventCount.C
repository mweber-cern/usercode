#include <stdlib.h>

#include <iostream>
#include <iomanip>
#include <fstream>

#include "TROOT.h"
#include "TObject.h"
#include "TFile.h"
#include "TChain.h"
#include "TSystem.h"
#include "TH1F.h"

int gLogLevel = 3;

using namespace std;

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
      cerr << "ERROR: Could not create TFile object for " << filename << endl;
      return 0;
    }
    if (!f->IsOpen()) {
      cerr << "ERROR: Could not open file " << filename << endl;
      return 0;
    }
  }
  TObject * obj = f->Get(objectname);
  if (obj == 0) {
    cerr << "Could not get object " << objectname << " from file " << filename << endl;
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

int main(int argc, char *argv[]) 
{
  if (argc != 3) {
    cout << "Usage: EventCount path counter.txt" << endl;
    cout << "Where path points to a directory containing ROOT files" << endl;
    cout << "and counter.txt summarizes the event count " << endl;
    exit (1);
  }
  // Loop over all files in directory <path> and fill
  // root file names into list <filelist>
  Text_t   lpath[strlen(argv[1])+8];
  sprintf(lpath, "%s", argv[1]);
  Text_t * basepath  = gSystem->ExpandPathName(lpath); // delete it later...
  void   * dirhandle = gSystem->OpenDirectory(basepath);
  const Text_t * basename;
  if (gLogLevel)
    cout << "Reading file names from directory " << basepath << endl;
  TList  * filelist  = new TList;
  while ((basename = gSystem->GetDirEntry(dirhandle))) {
    // Skip non-ROOT files
    if (!strstr(basename, ".root")) 
      continue;
    Text_t * fullname = new Text_t[strlen(basepath)+strlen(basename)+2];
    strcpy(fullname, basepath);
    strcat(fullname, "/");
    strcat(fullname, basename);
    filelist->Add(new TObjString(fullname));
    delete fullname;
  }

  // check if we have read some files
  if (filelist->GetSize() < 1) {
    cerr << "ERR: No *.root files available in " << basepath << " Exiting. " << endl;
    exit(EXIT_FAILURE);
  }

  long int counters[3] = { 0 };

  // Add all root files to a TChain
  if (gLogLevel)
    cout << "Adding files to TChain" << endl;
  
  TChain * chain = new TChain("ACSkimAnalysis/allData");
  TIter next(filelist);
  while (TObjString * obj = (TObjString *) next()) {
    chain->Add(obj->GetString());
    if (gLogLevel)
      cout << "FILE: " << obj->GetString() << " added to TChain" << endl;

    // open file, get histograms, add to counters
    TH1F * h = (TH1F *) get_object(obj->GetString(), "ACSkimAnalysis/h_counters");
    if (h == 0) {
      cerr << "ERR: Could not find histogram h_counters in file" << endl;
    }
    else {
      counters[0] += h->GetBinContent(1);
      counters[1] += h->GetBinContent(2);
      counters[2] += h->GetBinContent(3);
    }
  }

  // Delete all filelist elements and the list itself
  // (They are not deleted automatically when the TList is deleted!!!)
  filelist->Delete();
  delete filelist;
  delete basepath;


  long int events_from_tree = chain->GetEntries();

  long int fcounter[3] = { 0 };
  // open text file with event count
  char buffer[255];
  ifstream infile(argv[2]);
  while (infile.getline(buffer, 255)) {
    long int nev;
    if (sscanf(buffer, " ==>  #events read    :  %ld", &nev ) == 1) {
      fcounter[0] = nev;
    }
    if (sscanf(buffer, " ==>  #events presel  :  %ld", &nev ) == 1) {
      fcounter[1] = nev;
    }
    if (sscanf(buffer, " ==>  #events written :  %ld", &nev ) == 1) {
      fcounter[2] = nev;
    }
  }

  // Output a table
  cout << "Source                      Read         Presel        Written" << endl;
  cout << "counter.txt      "
       << setw(15) << fcounter[0]
       << setw(15) << fcounter[1]
       << setw(15) << fcounter[2] << endl;
  cout << "h_counter        " 
       << setw(15) << counters[0]
       << setw(15) << counters[1]
       << setw(15) << counters[2] << endl;
  cout << "Tree             " 
       << setw(45) << events_from_tree << endl;

  bool OK = true;
  // Analyse numbers from different sources
  if (fcounter[0] != counters[0]) {
    cerr << "ERROR: Read events do not agree" << endl;
    OK = false;
  }
  if (fcounter[1] != counters[1]) {
    cerr << "ERROR: Presel events do not agree" << endl;
    OK = false;
  }
  if (fcounter[2] != counters[2] || fcounter[2] != events_from_tree) {
    cerr << "ERROR: Written events do not agree" << endl;
    OK = false;
  }
  if (OK) {
    cout << endl << "OK: Events match!!" << endl;
  }
  else 
    exit(2);

  return 0;
}
