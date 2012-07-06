// local includes
#include "Utilities.h"
#include "Analysis.h"

// ROOT includes
#include "TStopwatch.h"
#include "TSystem.h"
#include "TString.h"
#include "TEnv.h"
#include "TChain.h"
#include "TFile.h"

// C/C++ includes
#include <stdlib.h>
#include <iostream>

using namespace std;
using namespace ROOT;

enum ErrorCodes { E_WRONG_PARAMS, E_OUTPUTFILE, E_NO_TREE, E_NO_ANALYSIS };

int main(int argc, char *argv[])
{
  // process and timing information
  ProcInfo_t info;
  TStopwatch timer;
  timer.Start();

  // get date and time
  time_t rawtime;
  struct tm * timeinfo;
  time(&rawtime);
  timeinfo = localtime ( &rawtime );

  // user greeting
  INFO("");
  INFO("######################################################################");
  INFO("#                                                                    #");
  INFO("# Start executing analyzer at " << asctime(timeinfo)		       );
  INFO("#                                                                    #");
  INFO("######################################################################");
  INFO("");

  // check command line arguments
  if (argc != 2) {
    ERROR("Usage: analyzer configfile ");
    return E_WRONG_PARAMS;
  }

  // read configuration file
  TEnv cfgFile(argv[1]);

  INFO("Using configuration file: " << argv[1]);

  // configure values from config file
  gLogLevel                           = cfgFile.GetValue("LogLevel", 3);
  string inputFiles                   = cfgFile.GetValue("InputFiles", "in.root");
  const char * outputFileName         = cfgFile.GetValue("OutputFileName", "out.root");

  TChain * chain = new TChain("ACSkimAnalysis/allData");
  // try to find out if inputfiles is a directory.
  Text_t * basepath  = gSystem->ExpandPathName(inputFiles.c_str());
  void   * dirhandle = gSystem->OpenDirectory(basepath);
  if (dirhandle != 0) {
    const Text_t * basename;
    while ((basename = gSystem->GetDirEntry(dirhandle))) {
      // Skip non-ROOT files
      if (!strstr(basename, ".root"))
	continue;
      string fullname = basepath + string("/") + basename;
      INFO("Adding file " << fullname);
      chain->Add(fullname.c_str());
    }
  }
  else {
    // OK, this is not a directory, add the given file(s)
    vector<string> files;
    split(inputFiles, ' ', files);
    for (vector<string>::const_iterator it = files.begin(); it != files.end(); it++) {
      INFO("Adding file " << *it);
      chain->Add(it->c_str());
    }
  }
  delete basepath;

  // create output file
  INFO("Creating output file and cloning tree")
  TFile * outFile = new TFile(outputFileName, "RECREATE");
  if (outFile == 0 || !outFile->IsOpen()) {
    ERROR("Could not open output file " << outputFileName);
    return E_OUTPUTFILE;
  }
  outFile->mkdir("ACSkimAnalysis");
  outFile->cd("ACSkimAnalysis");

  // Turn off branches 
  std::vector<string> removeBranches;
  split(cfgFile.GetValue("RemoveBranches", "None"), ',', removeBranches);
  std::vector<string>::const_iterator str;
  for (str = removeBranches.begin(); str != removeBranches.end(); str++) {
    if (*str == "None" || *str == "none")
      continue;
    INFO("Disabling branch(es) " << *str);
    chain->SetBranchStatus(str->c_str(), 0);
  }
  // clone tree structure, do not copy events...
  TTree * outTree = chain->CloneTree(0);
  if (outTree == 0) {
    ERROR("Could not clone tree from input file - maybe no tree in input file?");
    return E_NO_TREE;
  }
  // histograms belong to ROOT directory...
  outFile->cd();
  
  // create analysis object
  Analysis * analysis = 0;
  try {
    analysis = new Analysis(*chain, *outTree, cfgFile);
  }
  CATCH;
  if (analysis == 0) {
    ERROR("Could not create analysis class - out of memory?");
    return E_NO_ANALYSIS;
  }

  // mem info
  gSystem->GetProcInfo(&info);
  INFO("resident mem (MB) : " << info.fMemResident/1000.);
  INFO("virtual  mem (MB) : " << info.fMemVirtual/1000.);

  // start loop
  INFO("starting event loop");
  try {
    analysis->Loop();
  }
  CATCH;
  INFO("end of loop");

  // save data in file
  INFO("Saving data to file and closing...");
  // get currrent file - needed because of output tree spanning different
  // files, via setting TTree::SetMaxTreeSize
  outFile = outTree->GetCurrentFile();
  outFile->Write();
  outFile->Close();
  delete outFile;
  // no, we do not need to delete the objects in the file, this is done on the file delete... 
  
  // delete analysis object
  INFO("Deleting analysis");
  delete analysis;

  // time and memory info
  timer.Stop();
  gSystem->GetProcInfo(&info);
  INFO("");
  INFO("real time (s)     : " << timer.RealTime());
  INFO("CPU time (s)      : " << timer.CpuTime()); 
  INFO("resident mem (MB) : " << info.fMemResident/1000.);
  INFO("virtual  mem (MB) : " << info.fMemVirtual/1000.);
  INFO("");

  // info at end of program
  time(&rawtime);
  timeinfo = localtime ( &rawtime );
  INFO("");
  INFO("######################################################################");;
  INFO("#                                                                    #");;
  INFO("# End executing analyzer at " << asctime(timeinfo)		       );
  INFO("#                                                                    #");;
  INFO("######################################################################");;
  INFO("");
  
  return 0;
}
