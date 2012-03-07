#include <stdlib.h>

#include <iostream>
#include <fstream>

#include "TROOT.h"
#include "TChain.h"
#include "TSystem.h"

int gLogLevel = 3;

using namespace std;

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

  // Add all root files to a TChain
  if (gLogLevel)
    cout << "Adding files to TChain" << endl;

  TChain * chain = new TChain("ACSkimAnalysis/allData");
  TIter next(filelist);
  while (TObjString * obj = (TObjString *) next()) {
    chain->Add(obj->GetString());
    if (gLogLevel)
      cout << "FILE: " << obj->GetString() << " added to TChain" << endl;
  }

  // Delete all filelist elements and the list itself
  // (They are not deleted automatically when the TList is deleted!!!)
  filelist->Delete();
  delete filelist;
  delete basepath;

  long int events_from_tree = chain->GetEntries();
  cout << "Found a total of " << events_from_tree << " events in chain" << endl;

  // open text file with event count
  char buffer[255];
  ifstream infile(argv[2]);
  while (infile.getline(buffer, 255)) {
    long int nev;
    if (strstr(buffer, "written"))
      cout << buffer << endl;
    if (sscanf(buffer, " ==>  #events written :  %ld", &nev ) == 1) {
      cout << "Events from file " << argv[2] << " : " << nev << endl;
      if (nev != events_from_tree) {
	cerr << endl << "ERROR: Number of events from tree and file do not match!" << endl;
	exit(2);
      }
      else {
	cout << endl << "OK: Events match!!" << endl;
	break;
      }
    }
  }

  return 0;
}
