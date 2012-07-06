#include <stdlib.h>

#include <iostream>
#include <fstream>

#include "TROOT.h"
#include "TChain.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"

int gLogLevel = 3;

using namespace std;

int main(int argc, char *argv[]) 
{
  if (argc != 3 && argc != 4 && argc != 7) {
    cout << "Usage: quickDraw path drawstring [selection [nbins xmin xmax]]" << endl;
    cout << "       will create a file histogram.pdf from the given files" << endl 
	 << endl;
    cout << "path           points to a directory containing ROOT files" << endl;
    cout << "drawstring     is a string usable in TTree::Draw()" << endl;
    cout << "selection      e.g. (muo_n>0)&&(jet_n>0) or 1. if no selection" << endl;
    cout << "nbins          number of bins in histogram" << endl;
    cout << "xmin           low x edge of histogram" << endl;
    cout << "xmax           high x edge of histogram" << endl;
    cout << endl;
    cout << "Example: quickDraw ~/data \"muo_pt[0]\" \"(jet_n>2)\" 100 0 1000" << endl;
    cout << "          will draw muo_pt of the leading muon (index 0) for all events" << endl;
    cout << "          that contain at least three jets in a histogram with 100 bins," << endl;
    cout << "          ranging from 0 to 1000 GeV" << endl;
    cout << endl;
    cout << "Example: quickDraw ~/data vtx_n global_weight" << endl;
    cout << "         will draw vtx_n for all events in an automatically binned histogram " << endl;
    cout << "         The histogram weight will be taken from the leaf global_weight for each event" << endl;
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
    if (gLogLevel > 3)
      cout << "FILE: " << obj->GetString() << " added to TChain" << endl;
  }

  // Delete all filelist elements and the list itself
  // (They are not deleted automatically when the TList is deleted!!!)
  filelist->Delete();
  delete filelist;
  delete basepath;

  long int events_from_tree = chain->GetEntries();
  cout << "Found a total of " << events_from_tree << " events in chain" << endl;
  
  gROOT->SetBatch();
  TCanvas * c1 = new TCanvas();
  TFile * f  = new TFile("histo.root", "RECREATE");
  Int_t nbins = 200;
  double xmin = 0;
  double xmax = 1;
  if (argc == 7) {
    nbins = atoi(argv[4]);
    xmin = atof(argv[5]);
    xmax = atof(argv[6]);
  }
  TH1D * h1 = new TH1D("h1", argv[2], nbins, xmin, xmax);
  if (argc < 5) {
    h1->SetBit(TH1::kCanRebin);
  }

  const char * selection = "1.";
  if (argc >= 4)
    selection = argv[3];
  
  chain->Draw(Form("%s>>h1", argv[2]), selection);
  h1->Draw();
  c1->Print("histo.pdf");
  f->Write();
  f->Close();
  delete f;
  return 0;
}
