#include <sstream>

using namespace std;

#include "TLatex.h"
#include "TColor.h"
#include "TMath.h"
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TIterator.h"
#include "THashList.h"

#include "plot.h"

//////////////////////////////////////////////////////////////////////
// initialize variables from plot.h
Int_t           gLogLevel = 3;
TEnv          * gConfig = 0;
Int_t           gMaxProcess = 0;
TCanvas       * gCanvas = 0;
Int_t           nPad = 0;
Bool_t          gGerman = kFALSE;
Bool_t          gAutoOrder = kTRUE;
Int_t           gMaxPeriod = 0;
TString         gStage = "0";
Int_t           gPadNr = 0;
TH1D         ** gHisto[gMaxPad] = { 0 };   // size gMaxProcess
TH1D         ** gStack[gMaxPad] = { 0 };   // size gMaxProcess
TH1D         ** gShadow[gMaxPad] = { 0 };  // size gMaxProcess
TH2D         ** gHisto2 = 0;           // size gMaxProcess
TH3D         ** gHisto3[gMaxPad] = { 0 };   // size gMaxProcess
TProcess      * gProcess; // size gMaxProcess
TArrow        * ArrowMin[gMaxPad] = { 0 };
TArrow        * ArrowMax[gMaxPad] = { 0 };
Int_t         * gOrder[gMaxPad] = { 0 };
char         ** gPeriod = 0;
TEnv         ** gCuts = 0;
Int_t           gStart = 0;
Int_t           gEnd = 0;
TProcessInfo ** gProcessInfo = 0;
Double_t      * gLumi = 0;
Int_t           gMaxSignal = 0;
const char    * gSubDir = ".";
const char    * gBase = 0;
Bool_t          gIsColor = kTRUE;
Int_t           gGrid = kFALSE;
Bool_t          gMoveOverflow = kFALSE;
Bool_t          gHasOverflow[gMaxPad] = { 0 };
TH1D          * gFitSignal = 0;
TH1D          * gFitBackground = 0;
TH1D          * gFitData = 0;
Int_t           gVersion = 1;


// helper to duplicate existing string with operator new
char * strdup_new(const char * text)
{
  // duplicate string with help of new
  char * temp =  new char[strlen(text)+1];
  strcpy(temp, text);
  return temp;
}

//////////////////////////////////////////////////////////////////////
// options for TStyle, TH1, TH2, TLegend, TGraph, ...

void setopt(TStyle * style)
{
  style->SetPalette(1);
  style->SetOptTitle(0);          // don't show histogram title
  style->SetPadLeftMargin(0.15);
  return;
  // style options
  style->SetStatColor(10);        // white background for statistics box
//   style->SetTitleColor(10);    // title background white
  style->SetFillColor(10);        // fill everything with a white background
  style->SetFillStyle(1001);      // solid fill style
  style->SetCanvasBorderMode(0);  // no borders in canvases please!
  style->SetCanvasBorderSize(0);  // no borders in canvases please!
  style->SetCanvasColor(10);      // white canvases please
  style->SetPadBorderMode(0);     // no borders in Pads please!
  style->SetPadBottomMargin(0.15);
  style->SetPadRightMargin(0.04);
  style->SetPadTopMargin(0.07);
  style->SetPadBorderSize(0);     // no borders in Pads please!
  style->SetPadColor(10);         // white pads please
  style->SetFrameBorderMode(0);
  style->SetFrameBorderSize(0);
  style->SetOptStat(0);           // don't show statistics
  style->SetOptFit(1111);         // show fit parameters, probab., chi^2
  style->SetErrorX(0);            // no horizontal error
  style->SetTitleOffset(3, "Y");  // title offset
  style->SetTitleOffset(3, "X");  // title offset
  style->SetLabelSize(0.03, "Y"); // label size
  style->SetLabelSize(0.03, "X"); // label size
  style->SetTitleXSize(0.03);     // title size
  style->SetTitleYSize(0.03);     // title size
}

void setopt(TCanvas * canvas)
{
  canvas->SetLeftMargin(0.15);
  return;
  // set default options for canvas
  canvas->SetTopMargin(0.03);
  canvas->SetRightMargin(0.05);
  canvas->SetBottomMargin(0.15);
  canvas->SetFillColor(kWhite);
  canvas->SetFillStyle(1001);
  canvas->SetBorderSize(0);
  canvas->SetBorderMode(0);
}

void setopt(TH1 * histo)
{
  // set histo default options
  histo->SetTitleOffset(1.3, "Y"); // title offset
  return;
  histo->SetLabelSize(0.06, "X"); // title offset
  histo->SetLabelSize(0.06, "Y"); // title offset
  histo->SetLabelSize(0.06, "Z"); // title offset
  histo->SetTitleSize(0.06, "X");  // title size
  histo->SetTitleSize(0.06, "Y");  // title size
  histo->SetTitleSize(0.06, "Z");  // title size
  histo->SetTitleOffset(1.1, "X"); // title offset
  histo->SetTitleOffset(1.1, "Z"); // title offset
  histo->SetMarkerStyle(8);
  histo->SetMarkerSize(1.1);
//    histo->GetXaxis()->SetTitleColor(kBlack);
//    histo->GetYaxis()->SetTitleColor(kBlack);
//    histo->GetXaxis()->CenterTitle();
}

void setopt(TLegend * leg)
{
  // set legend default option
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
}

void setopt(TGraph * gr)
{
  return;
  TH1F * histo = gr->GetHistogram();
  if (histo != 0)
    setopt(histo);
  gr->SetMarkerStyle(8);
  gr->SetMarkerSize(1.1);
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

void delete_old_config()
{
  DEBUG("delete_old_config() start");
  if (gConfig) {
    delete gConfig;
    gConfig = 0;
  }
  if (gPeriod) {
    delete[] gPeriod;
    gPeriod = 0;
  }
  if (gCuts) {
    delete[] gCuts;
    gCuts = 0;
  }
  if (gLumi) {
    delete[] gLumi;
    gLumi = 0;
  }
  for (Int_t pad = 0; pad < gMaxPad; pad++) {
    if (gHisto[pad]) {
      delete[] gHisto[pad];
      gHisto[pad] = 0;
    }
    if (gStack[pad]) {
      delete[] gStack[pad];
      gStack[pad] = 0;
    }
    if (gShadow[pad]) {
      delete[] gShadow[pad];
      gShadow[pad] = 0;
    }
    if (gOrder[pad]) {
      delete[] gOrder[pad];
      gOrder[pad] = 0;
    }
    if (gHisto3[pad]) {
      delete[] gHisto3[pad];
      gHisto3[pad] = 0;
    }
  }
  if (gHisto2) {
    delete[] gHisto2;
    gHisto2 = 0;
  }
  // @todo: need to delete everything created with strdup_new here...
  if (gProcess) {
    delete[] gProcess;
    gProcess = 0;
  }
  if (gBase) {
    delete gBase;
    gBase = 0;
  }
  DEBUG("delete_old_config() end");
}

// read all config files
void read_config_files(const char * configFileName = "plot.cfg")
{
  INFO("Opening config file " << configFileName);
  gConfig = new TEnv(configFileName);

  THashList * gTable = (THashList *) gConfig->GetTable();
  TIter next(gTable);
  TEnvRec * obj = 0;
  
  // number of data taking periods, backgrounds, signals
  bool signalToggle = false;
  TString processName = "toast";
  vector<TString> periods;
  vector<TString> backgrounds;
  vector<TString> signals;
 
  

  while ( (obj = (TEnvRec *) next()) ) {
    TString objName = obj->GetName();
    // skip if settings or data
    if (objName.BeginsWith("settings") || objName.BeginsWith("data")) continue;

    // count periods
    if (objName.BeginsWith("period")) {
      periods.push_back(objName);
      continue;
    }

    // signal or background?
    if (objName == "background") {
      signalToggle = false;
      continue;
    }
    else if (objName == "signal") {
      signalToggle = true;
      continue;
    }

    // count signals and periods
    if (!objName.BeginsWith(processName+'.')) {
      
      TObjArray * objStrArr = (TObjArray *) objName.Tokenize(".");
      TObjString * objStrEle = (TObjString *) objStrArr->At(0);
      processName = objStrEle->GetString();

      if (signalToggle) {
	signals.push_back(processName);
      }
      else {
	backgrounds.push_back(processName);
      }
    }
  }

  // get number of periods
  gMaxPeriod = periods.size();
  if (gMaxPeriod < 1) {
    ERROR("No periods found in config file");
    return;
  }
  INFO("Configuring for " << gMaxPeriod << " period(s)");

  // get number of background MCs
  Int_t N_bg = backgrounds.size();
  if (N_bg < 0) {
    ERROR("Number of backgrounds found in config file is < 0");
    return;
  }
  INFO("Number of backgrounds: " << N_bg);

  // get number of signal MCs
  gMaxSignal = signals.size();
  if (gMaxSignal < 0) {
    ERROR("Number of signals N_sg found in config file < 0");
    return;
  }
  INFO("Number of signals: " << gMaxSignal);

  // create and initialize period depending arrays
  gPeriod = new char * [gMaxPeriod];
  gCuts   = new TEnv * [gMaxPeriod];
  gLumi   = new Double_t[gMaxPeriod];
  for (Int_t i = 0; i < gMaxPeriod; i++) {
    DEBUG("Period #" << i << ": " << gConfig->GetValue(Form("%s", periods.at(i).Data()), ""));
    gPeriod[i] = strdup_new(gConfig->GetValue(Form("%s", periods.at(i).Data()), ""));
    gCuts[i] = 0;
    gLumi[i] = 0;
  }

  // plot all periods by default
  gStart = 0; gEnd = gMaxPeriod-1;

  // number of processes
  gMaxProcess = N_bg + gMaxSignal + 1;
  
  // create and initialize arrays depending on gMaxProcess
  for (Int_t pad = 0; pad < gMaxPad; pad++) {
    gHisto[pad]  = new TH1D * [gMaxProcess];
    gStack[pad]  = new TH1D * [gMaxProcess];
    gShadow[pad] = new TH1D * [gMaxProcess];
    gOrder[pad]  = new Int_t[gMaxProcess];
    gHisto3[pad] = new TH3D * [gMaxProcess];
    for (Int_t i = 0; i < gMaxProcess; i++) {
      gHisto[pad][i]  = 0;
      gStack[pad][i]  = 0;
      gShadow[pad][i] = 0;
      gOrder[pad][i] = i;
      gHisto3[pad][i] = 0;
    }
  }
  gHisto2      = new TH2D * [gMaxProcess];
  gProcess     = new TProcess[gMaxProcess];
  for (Int_t i = 0; i < gMaxProcess; i++) {
    gHisto2[i] = 0;
  }

  // automatic color selection
  // read all processes
  // signal and background
  for (Int_t i = 0; i < gMaxProcess-1; i++) {
    const char * selector = 0;
    if (i < gMaxSignal) {
      selector = strdup_new(Form("%s", signals.at(i).Data()));
    }
    else {
      selector = strdup_new(Form("%s", backgrounds.at(i-gMaxSignal).Data()));
    }
    gProcess[i].fname = strdup_new(selector);
    gProcess[i].tname = strdup_new(gConfig->GetValue(Form("%s.label", selector), gProcess[i].fname));
    const char * fcolor = gConfig->GetValue(Form("%s.fcolor", selector), "1");
    gProcess[i].fcolor = gROOT->ProcessLine(fcolor);
    const char * lcolor = gConfig->GetValue(Form("%s.lcolor", selector), fcolor);
    gProcess[i].lcolor = gROOT->ProcessLine(lcolor);
    gProcess[i].lstyle = gConfig->GetValue(Form("%s.lstyle", selector),  1);
    gProcess[i].hcolor = gConfig->GetValue(Form("%s.hcolor", selector),  kWhite);
    gProcess[i].hstyle = gConfig->GetValue(Form("%s.hstyle", selector),  1001);
    gProcess[i].marker = gConfig->GetValue(Form("%s.marker", selector),  8);
    gProcess[i].mcolor = gConfig->GetValue(Form("%s.mcolor", selector),  1);
    gProcess[i].msize  = gConfig->GetValue(Form("%s.msize", selector),  0.7);
    gProcess[i].join   = gConfig->GetValue(Form("%s.join", selector), kFALSE);
    gProcess[i].stack  = gConfig->GetValue(Form("%s.stack", selector), kTRUE);
    if (i == 0 || i == gMaxSignal) {
      if (gProcess[i].join) {
	ERROR("Misconfiguration for process " << selector << ":");
	ERROR("You are not allowed to join the first signal or background histogram");
      }
    }
  }
  // data
  gProcess[gMaxProcess-1].fname     = strdup_new("data_doublemu");
  gProcess[gMaxProcess-1].tname     = strdup_new(gConfig->GetValue("data_doublemu.label", "Data"));
  gProcess[gMaxProcess-1].fcolor    = kWhite;
  gProcess[gMaxProcess-1].lcolor    = kWhite;
  gProcess[gMaxProcess-1].hcolor    = kWhite;
  gProcess[gMaxProcess-1].hstyle    = 1001;
  gProcess[gMaxProcess-1].lstyle    = 1;
  gProcess[gMaxProcess-1].marker    = 8;
  gProcess[gMaxProcess-1].mcolor    = 1;
  gProcess[gMaxProcess-1].msize     = .7;
  gProcess[gMaxProcess-1].join      = kFALSE;
  gProcess[gMaxProcess-1].stack     = kFALSE;

  // misc setup
  gBase = strdup_new(gConfig->GetValue("settings.basedir", "."));

  // read xs and lumi for each process
  gProcessInfo = new TProcessInfo * [gMaxPeriod];
  for (Int_t period = 0; period < gMaxPeriod; period++) {
    gProcessInfo[period] = new TProcessInfo[gMaxProcess-1];
    for (Int_t process = 0; process < gMaxProcess-1; process++) {
      DEBUG("Option " << Form("%s.%s.xs", gProcess[process].fname, gPeriod[period])
	    << " has value " << gConfig->GetValue(Form("%s.%s.xs", gProcess[process].fname, gPeriod[period]), gOptDefault) );
      gProcessInfo[period][process].xs = gConfig->GetValue(Form("%s.%s.xs", gProcess[process].fname, gPeriod[period]), gOptDefault);
      gProcessInfo[period][process].Nev = gConfig->GetValue(Form("%s.%s.Nev", gProcess[process].fname, gPeriod[period]), -1);
      gProcessInfo[period][process].weight = gConfig->GetValue(Form("%s.%s.weight", gProcess[process].fname, gPeriod[period]), 1.);
      
      DEBUG("Nev    = " << gProcessInfo[period][process].Nev);
      DEBUG("xs     = " << gProcessInfo[period][process].xs);
      DEBUG("weight = " << gProcessInfo[period][process].weight);

      if (gProcessInfo[period][process].xs == gOptDefault) {
	ERROR("Unreasonable or no cross-section " << gProcessInfo[period][process].xs
	      << " for process " << gProcess[process].fname << " in file " << configFileName);
      }
      if (gProcessInfo[period][process].Nev < 0) {
	ERROR("Unreasonable event number " << gProcessInfo[period][process].Nev
	      << " for process " << gProcess[process].fname << " in file " << configFileName);
      }
      if (gProcessInfo[period][process].weight < 0) {
	ERROR("Unreasonable weight " << gProcessInfo[period][process].weight
	      << " for process " << gProcess[process].fname << " in file " << configFileName);
      }
    }
    gLumi[period] = gConfig->GetValue(Form("data_doublemu.%s.lumi", gPeriod[period]), gOptDefault);
    if (gLumi[period] == gOptDefault || gLumi[period] < 0) {
      ERROR("Bad lumi " << gLumi[period] << " from file " << configFileName);
    }
  }
}

// read configuration from file
void read_config_file(const char * configFileName = "Overview.cfg")
{
  INFO("Opening config file " << configFileName);
  gConfig = new TEnv(configFileName);

  // number of data taking periods
  gMaxPeriod = gConfig->GetValue("periods", -1);
  if (gMaxPeriod < 1) {
    ERROR("Number of periods not found in config file or periods < 1");
    return;
  }
  INFO("Configuring for " << gMaxPeriod << " periods");

  // create and initialize period depending arrays
  gPeriod = new char * [gMaxPeriod];
  gCuts   = new TEnv * [gMaxPeriod];
  gLumi   = new Double_t[gMaxPeriod];
  for (Int_t i = 0; i < gMaxPeriod; i++) {
    DEBUG("Period #" << i << ": "
	  << gConfig->GetValue(Form("period.%d.name", i), ""));
    gPeriod[i] = strdup_new(gConfig->GetValue(Form("period.%d.name", i), ""));
    gCuts[i] = 0;
    gLumi[i] = 0;
  }
  // plot all periods by default
  gStart = 0; gEnd = gMaxPeriod-1;

  // get number of background MCs
  Int_t N_bg = gConfig->GetValue("N_bg", -1);
  if (N_bg < 0) {
    ERROR("Number of backgrounds N_bg not found in config file or N_bg < 0");
    return;
  }
  INFO("Number of backgrounds: " << N_bg);

  // get number of signal MCs
  gMaxSignal = gConfig->GetValue("N_sg", -1);
  if (gMaxSignal < 0) {
    ERROR("Number of signals N_sg not found in config file or N_sg < 0");
    return;
  }
  INFO("Number of signals: " << gMaxSignal);

  // number of processes
  gMaxProcess = N_bg + gMaxSignal + 1;

  // create and initialize arrays depending on gMaxProcess
  for (Int_t pad = 0; pad < gMaxPad; pad++) {
    gHisto[pad]  = new TH1D * [gMaxProcess];
    gStack[pad]  = new TH1D * [gMaxProcess];
    gShadow[pad] = new TH1D * [gMaxProcess];
    gOrder[pad]  = new Int_t[gMaxProcess];
    gHisto3[pad] = new TH3D * [gMaxProcess];
    for (Int_t i = 0; i < gMaxProcess; i++) {
      gHisto[pad][i]  = 0;
      gStack[pad][i]  = 0;
      gShadow[pad][i] = 0;
      gOrder[pad][i] = i;
      gHisto3[pad][i] = 0;
    }
  }
  gHisto2      = new TH2D * [gMaxProcess];
  gProcess     = new TProcess[gMaxProcess];
  for (Int_t i = 0; i < gMaxProcess; i++) {
    gHisto2[i] = 0;
  }

  // automatic color selection
  // read all processes
  // signal and background
  for (Int_t i = 0; i < gMaxProcess-1; i++) {
    const char * selector = 0;
    if (i < gMaxSignal) {
      selector = strdup_new(Form("sg.%d", i));
    }
    else {
      selector = strdup_new(Form("bg.%d", i-gMaxSignal));
    }
    gProcess[i].fname = strdup_new(gConfig->GetValue(Form("%s.file", selector),
						 strdup_new(selector)));
    gProcess[i].tname = strdup_new(gConfig->GetValue(Form("%s.label", selector),
						 gProcess[i].fname));
    gProcess[i].fcolor = gConfig->GetValue(Form("%s.fcolor", selector),  i);
    gProcess[i].lcolor = gConfig->GetValue(Form("%s.lcolor", selector),  gProcess[i].fcolor);
    gProcess[i].lstyle = gConfig->GetValue(Form("%s.lstyle", selector),  1);
    gProcess[i].hcolor = gConfig->GetValue(Form("%s.hcolor", selector),  kWhite);
    gProcess[i].hstyle = gConfig->GetValue(Form("%s.hstyle", selector),  1001);
    gProcess[i].marker = gConfig->GetValue(Form("%s.marker", selector),  8);
    gProcess[i].mcolor = gConfig->GetValue(Form("%s.mcolor", selector),  1);
    gProcess[i].msize  = gConfig->GetValue(Form("%s.msize", selector),  0.7);
    gProcess[i].join   = gConfig->GetValue(Form("%s.join", selector), kFALSE);
    gProcess[i].stack  = gConfig->GetValue(Form("%s.stack", selector), kTRUE);
    if (i == 0 || i == gMaxSignal) {
      if (gProcess[i].join) {
	ERROR("Misconfiguration for process " << selector << ":");
	ERROR("You are not allowed to join the first signal or background histogram");
      }
    }
  }
  // data
  gProcess[gMaxProcess-1].fname = strdup_new(gConfig->GetValue("file", "data"));
  gProcess[gMaxProcess-1].tname = strdup_new(gConfig->GetValue("label", "Data"));
  gProcess[gMaxProcess-1].fcolor    = kWhite;
  gProcess[gMaxProcess-1].lcolor    = kWhite;
  gProcess[gMaxProcess-1].hcolor    = kWhite;
  gProcess[gMaxProcess-1].hstyle    = 1001;
  gProcess[gMaxProcess-1].lstyle    = 1;
  gProcess[gMaxProcess-1].marker    = 8;
  gProcess[gMaxProcess-1].mcolor    = 1;
  gProcess[gMaxProcess-1].msize     = .7;
  gProcess[gMaxProcess-1].join      = kFALSE;
  gProcess[gMaxProcess-1].stack     = kFALSE;

  // misc setup
  gBase = strdup_new(gConfig->GetValue("basedir", "."));
}

// read all cross-section and luminosity information
void read_xs_and_lumi()
{
  gProcessInfo = new TProcessInfo * [gMaxPeriod];
  for (Int_t period = 0; period < gMaxPeriod; period++) {
    char * fpath = strdup_new(gConfig->GetValue(Form("period.%d.xs", period),
						Form("xs_%s.cfg", gPeriod[period])));
    INFO("Reading cross-section and lumi from file " << fpath);
    char * expath = gSystem->ExpandPathName(fpath);
    TEnv * config = new TEnv(expath);
    delete expath;
    gProcessInfo[period] = new TProcessInfo[gMaxProcess-1];
    for (Int_t process = 0; process < gMaxProcess-1; process++) {
      DEBUG("Option " << Form("%s.xs", gProcess[process].fname)
	    << " has value " << config->GetValue(Form("%s.xs", gProcess[process].fname),
						 gOptDefault) );
      gProcessInfo[period][process].xs = config->GetValue(Form("%s.xs",
							       gProcess[process].fname),
							  gOptDefault);
      gProcessInfo[period][process].Nev = config->GetValue(Form("%s.Nev",
							       gProcess[process].fname),
							   -1);
      gProcessInfo[period][process].weight = config->GetValue(Form("%s.weight",
							       gProcess[process].fname),
							  1.);
      
      DEBUG("Nev    = " << gProcessInfo[period][process].Nev);
      DEBUG("xs     = " << gProcessInfo[period][process].xs);
      DEBUG("weight = " << gProcessInfo[period][process].weight);

      if (gProcessInfo[period][process].xs == gOptDefault) {
	ERROR("Unreasonable or no cross-section " << gProcessInfo[period][process].xs
	      << " for process " << gProcess[process].fname << " in file " << fpath);
      }
      if (gProcessInfo[period][process].Nev < 0) {
	ERROR("Unreasonable event number " << gProcessInfo[period][process].Nev
	      << " for process " << gProcess[process].fname << " in file " << fpath);
      }
      if (gProcessInfo[period][process].weight < 0) {
	ERROR("Unreasonable weight " << gProcessInfo[period][process].weight
	      << " for process " << gProcess[process].fname << " in file " << fpath);
      }
    }
    gLumi[period] = config->GetValue("lumi", gOptDefault);
    if (gLumi[period] == gOptDefault || gLumi[period] < 0) {
      ERROR("Bad lumi " << gLumi[period] << " from file " << fpath);
    }
    delete config;
    delete fpath;
  }
}

//////////////////////////////////////////////////////////////////////
// setup routines

void version(Int_t version)
{
  switch (version) {
    case 0:
    case 1:
      break;
    default:
      ERROR("Version unknown: " << version);
  }
  gVersion = version;
}

// this function should be called e.g. by your rootlogon.C, before any functions of the
// plot.C code are called.
void setup(const char * configFileName)
{
  INFO("Setting up plot macro from config directory");

  // set plot style
  setopt(gStyle);

  // use weighted histograms
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();

  // read configuration file
  const char * expConfigFilePath = gSystem->ExpandPathName(Form("../config/%s", configFileName));
  read_config_files(expConfigFilePath);
  delete expConfigFilePath;
}

// this function should be called e.g. by your rootlogon.C, before any functions of the
// plot.C code are called.
void setup_dumpallplots(const char * configFileName)
{
  INFO("Setting up plot macro");

  // set plot style
  setopt(gStyle);

  // use weighted histograms
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();

  // read configuration file
  const char * expConfigFileName = gSystem->ExpandPathName(Form("../config/%s", configFileName));
  read_config_file(expConfigFileName);
  delete expConfigFileName;

  // read cross-section and lumi
  read_xs_and_lumi();
}


//////////////////////////////////////////////////////////////////////
// routines for drawing histograms

void cd(Int_t canvas)
{
  if (canvas < 1 || canvas > gMaxPad) {
    ERROR("argument \"canvas\" needs to be between 1 and " << gMaxPad);
    return;
  }
  // switch to corresponding canvas
  if (!gCanvas)
    return;
  gPadNr = canvas-1;
  gCanvas->cd(canvas);
}

void MakeCanvas(Int_t dx, Int_t dy)
{
  // create analysis canvas with requested subdivisions
  if (gCanvas) {
    delete gCanvas;
    gCanvas = 0;
    nPad = 0;
  }

  // determine x and y sizes
  Int_t xsize, ysize;
  if (dy > dx) {
    xsize = 600;
    ysize = (Int_t) (600 * TMath::Sqrt(2.));
  }
  else {
    ysize = 600;
    xsize = (Int_t) (600 * TMath::Sqrt(2.));
  }

  // max number of subpads: 3x2
  if (dx * dy > gMaxPad) {
    ERROR("Too many subpads (" << dx*dy << "), max possible:" << gMaxPad);
    return;
  }

  // create canvas
  gCanvas = new TCanvas("gCanvas", "analysis", 0, 0, xsize, ysize);

  // set canvas style
  setopt(gCanvas);

  // make subpads and activate first pad
  gPadNr = 0;
  gCanvas->Divide(dx, dy);
  nPad = dx*dy;
  cd(gPadNr+1);
}

void MakeCanvas2(Int_t dy, Double_t percent = 0.8)
{
  // make another canvas
  if (gCanvas) {
    delete gCanvas;
  }

  // create canvas
  if (dy == 1)
    gCanvas = new TCanvas("gCanvas", "analysis", 0, 0,
			  (Int_t) (600*TMath::Sqrt(2.)*1.25), 600);
  else
    gCanvas = new TCanvas("gCanvas", "analysis", 0, 0,
			  (Int_t) (600*1.25), (Int_t) (600*TMath::Sqrt(2.)));

  // set canvas style

  // make subpads and activate first pad
  gPadNr = 0;
  Double_t ysize = 0.99/dy;
  char name[20];
  for (Int_t i = 0; i < dy; i++) {
    snprintf(name, 20, "gCanvas_%d", i+1);
    TPad * pad = new TPad(name, name, 0.001, 1-(i+1)*ysize, percent-0.01, 1-i*ysize,
			  kWhite);
    pad->SetNumber(i+1);
    pad->Draw();
    snprintf(name, 20, "gCanvas_%d", 10*(i+1));
    pad = new TPad(name, name, percent, 1-(i+1)*ysize, 0.99, 1-i*ysize, kWhite);
    pad->SetNumber(10*(i+1));
    pad->Draw();
  }
  gCanvas->Modified();
  // activate first Pad
  cd(gPadNr+1);
}

void padupdate()
{
  if (gPad == 0) 
    return;
  gPad->Modified();
  gPad->Update();
  gPad->RedrawAxis();
}

Int_t FindFirstHisto()
{
  // find a histogram that has been filled and return its index in gHisto[gPadNr]
  // if no histogram existent, return -1
  for (Int_t i = 0; i < gMaxProcess; i++) {
    Int_t process = gOrder[gPadNr][i];
    if (gHisto[gPadNr][process] != 0) {
      return process;
    }
  }
  WARNING("No histograms found - use plot()");
  return -1;
}

void color(Bool_t on)
{
  // set color on or off
  gIsColor = on;
}

void top()
{
  // order mc's to plot signal on top
  for (Int_t i = 0; i < gMaxProcess; i++) {
    gOrder[gPadNr][i] = i;
  }
}

void bottom()
{
  // plot signal at bottom
  // signal
  for (Int_t i = 0; i < gMaxSignal; i++) {
    gOrder[gPadNr][gMaxProcess-1-gMaxSignal+i] = i;
  }
  // background
  for (Int_t i = gMaxSignal; i < gMaxProcess-1; i++) {
    gOrder[gPadNr][i-gMaxSignal] = i;
  }
  // data remains unchanged
  gOrder[gPadNr][gMaxProcess-1] = gMaxProcess-1;
}

const char * getpath(Int_t period)
{
  static char fpath[400];
  snprintf(fpath, 400, "%s/%s/%s", gBase, gSubDir, gPeriod[period]);
  return fpath;
}

void readcuts()
{
  // read in cut values from option file
  for (Int_t period = gStart; period <= gEnd; period++) {
    // create new TOption object
    if (gCuts[period] != 0)
      delete gCuts[period];
    INFO("Reading option file " << Form("%s/%s/config/analyzer.cfg", gBase, gSubDir));
    gCuts[period] = new TEnv(Form("%s/%s/config/analyzer.cfg", gBase, gSubDir));
  }
}

void selection(const char * subdir)
{
  // save in global variable
  if (gSubDir && strcmp(gSubDir, "."))
    delete gSubDir;
  gSubDir = strdup_new(subdir);

  // read in option file
  readcuts();
}

void period(const char * startperiod, const char * endperiod)
{
  if (endperiod == 0) {
    endperiod = startperiod;
  }

  for (gStart = 0; gStart < gMaxPeriod; gStart++) {
    if (!strcmp(gPeriod[gStart], startperiod))
      break;
  }
  if (gStart == gMaxPeriod) {
    ERROR("Period " << startperiod << " is not a valid period");
    return;
  }

  for (gEnd = 0; gEnd < gMaxPeriod; gEnd++) {
    if (!strcmp(gPeriod[gEnd], endperiod))
      break;
  }
  if (gEnd == gMaxPeriod) {
    ERROR("period " << endperiod << " is not a valid period");
      return;
  }
  if (gEnd < gStart) {
    ERROR("Change order of periods");
    return;
  }

  // read in option file
  readcuts();
}

void stage(TString s)
{
  // set stage
  gStage = s;
}

const char * GetUnit(const char * hname)
{
  DEBUG("GetUnit() start");
  // get unit of histogram name for labelling y-axis automatically
  static char * unit = 0;

  if (unit != 0) {
    DEBUG("deleting old unit: " << unit);
    delete[] unit;
    unit = 0;
  }

  // find position where to separate the two strings
  char * pos = strchr(hname, '@');

  if (pos == 0) {
    // no unit given, return null pointer
    return 0;
  }

  // unit given, go over separator
  pos++;

  // and skip leading blanks
  while (isspace(*pos))
    pos++;

  unit = strdup_new(pos);

  DEBUG("unit found: " << unit);
  return unit;
}

const char * GetXTitle(const TH1D * histo)
{
  // create new string containing only x title of histogram name
  static char * title = 0;

  DEBUG("GetXTitle() start");

  if (title) {
    delete[] title;
    title = 0;
  }

  DEBUG("histo title " <<  histo->GetTitle());

  // find position where to separate the two strings
  const char * pos = strchr(histo->GetTitle(), ':');

  if (pos == 0) {
    pos = histo->GetTitle();
  }
  else {
    // y-title given, go over separator
    pos++;

    // and skip leading blanks
    while (isspace(*pos))
      pos++;
  }

  DEBUG("x title after removing y part: " << pos);

  // find unit
  char * upos = strchr(pos, '@');
  if (upos != 0) {
    const char * unit = GetUnit(pos);
    // unit found
    DEBUG("title len = " << (int) (upos-pos));
    Long_t len = upos-pos;
    title = new char[len+strlen(unit)+4];
    strncpy(title, pos, len);
    title[len] = '\0';
    DEBUG("title now = " << title);
    strcat(title, " [");
    DEBUG("title now = " << title);
    strcat(title, unit);
    DEBUG("title now = " << title);
    strcat(title, "]");
    DEBUG("title now = " << title);
  }
  else {
    // no unit
    title = strdup_new(pos);
  }
  DEBUG("GetXTitle() end");
  return title;
}

const char * GetYTitle(const TH1D * histo)
{
  const Int_t max_size = 4096;
  // create new string containing only y title of histogram name
  static char * title = 0;  // the y-title
  DEBUG("GetYTitle() start");
  if (title) {
    DEBUG("deleting title " << title);
    delete[] title;
    title = 0;
  }
  // find position where to separate the two strings
  char * pos = strchr(histo->GetTitle(), ':');

  if (pos == 0) {
    DEBUG("no y title given");
    // check if bin names are given - in this case do not output bin width
    if (histo->GetXaxis()->GetLabels() != 0) {
      // assume these are just events
      title = new char[max_size];
      snprintf(title, max_size, "Number of events");
      return title;
    }
    // no y-title given
    // get unit
    const char * unit = GetUnit(histo->GetTitle());
    Double_t binsize = histo->GetBinWidth(1);
    // unit given
    if (unit) {
      if (gGerman) {
	title = new char[max_size];
	snprintf(title, max_size, "Anzahl Ereignisse / %4.2f %s", binsize, unit);
      }
      else {
	title = new char[max_size];
	snprintf(title, max_size, "Number of events / %4.2f %s", binsize, unit);
      }
    }
    else {
      if (gGerman) {
	title = new char[max_size];
	snprintf(title, max_size, "Anzahl Ereignisse / %4.2f", binsize);
      }
      else {
	title = new char[max_size];
	snprintf(title, max_size, "Number of events / %4.2f", binsize);
      }
    }
  }
  else {
    // OK, y title given

    // cut title x part away
    Int_t len = strlen(histo->GetTitle())-strlen(pos);
    char * ttitle = new char[len+1];
    strncpy(ttitle, histo->GetTitle(), len);
    ttitle[len] = '\0';
    DEBUG("y title after cutting x-part away: " << ttitle);
    
    // find if it contains a unit
    const char * unit = GetUnit(ttitle);
    if (unit) {
      // unit found
      title = new char[strlen(ttitle)+4];
      strncpy(title, ttitle, strlen(ttitle)-strlen(unit)-1);
      DEBUG("title now = " << title);
      strcat(title, " [");
      DEBUG("title now = " << title);
      strcat(title, unit);
      DEBUG("title now = " << title);
      strcat(title, "]");
      DEBUG("title now = " << title);
    }
    else {
      title = ttitle;
    }
  }

  DEBUG("GetYTitle() end");
  return title;
}

Double_t GetOpt(Int_t period, const char * hname, const char * MinMax)
{
  // get option value with suffix MinMax from histogram hname
  if (!gCuts[period]) {
    ERROR("option file not present, return default " << gOptDefault);
  }
  char optname[256];
  // get value
  snprintf(optname, 256, "%s_%s", hname+6, MinMax);
  Double_t value = gCuts[gStart]->GetValue(optname, gOptDefault);
  if (value == gOptDefault) {
    INFO("Option " << optname << " not found, return default " << gOptDefault);
  }
  return value;
}

Double_t GetOptMin(Int_t period, const char * hname)
{
  // get minimum value from option file for period gPeriod[gStart]
  return GetOpt(period, hname, "min");
}

Double_t GetOptMax(Int_t period, const char * hname)
{
  // get maximum value from option file for period gPeriod[gStart]
  return GetOpt(period, hname, "max");
}


TArrow * arrow(Double_t position, Int_t neighbourbins)
{
  DEBUG("start arrow");
  // print arrow to indicate a cut
  Int_t hstart = FindFirstHisto();
  if (hstart < 0)
    return 0;
  Double_t fmax = gStack[gPadNr][hstart]->GetMaximum();
  Double_t fmin = gStack[gPadNr][hstart]->GetMinimum();
  DEBUG("fmax = " << fmax);
  DEBUG("fmin = " << fmin);
  TAxis * axis = gStack[gPadNr][hstart]->GetXaxis();
  if (axis == 0) {
    ERROR("No axis found");
    return 0;
  }
  Int_t fbin = axis->FindFixBin(position);
  DEBUG("fbin = " << fbin << ", max = " << gStack[gPadNr][hstart]->GetNbinsX());
  Double_t fbottom = 1E-8;
  for (Int_t process = 0; process < gMaxProcess; process++) {
    if (gStack[gPadNr][process]) {
      for (Int_t i = fbin - neighbourbins; i <= fbin + neighbourbins; i++) {
        Double_t temp = gStack[gPadNr][process]->GetBinContent(i);
        if (temp > fbottom)
          fbottom = temp;
      }
    }
  }
  DEBUG("fbottom = " << fbottom);

  Double_t ftop;
  if (gPad->GetLogy() == 0) {
    fbottom = TMath::Min(fmax*0.95, fbottom += fmax * 0.1);
    ftop    = TMath::Min(fmax, fbottom + fmax * 0.25);
  } else {
    fbottom = TMath::Min(fmax/(TMath::Power(10, 0.05)*((TMath::Log10(fmax)
							-TMath::Log10(fmin)))),
			 fbottom*TMath::Power(10, 0.1*(TMath::Log10(fmax)
						       -TMath::Log10(fmin))));
    ftop    = TMath::Min(fmax,
			 fbottom*TMath::Power(10, 0.25*(TMath::Log10(fmax)
							-TMath::Log10(fmin))));
  }
  DEBUG("ftop = " << ftop);
  DEBUG("fbottom = " << fbottom);
  TArrow * a1 = new TArrow(position, ftop, position, fbottom);
  if (a1 == 0)
    return 0;
  a1->SetLineWidth(2);
  a1->SetArrowSize(0.03);
  a1->Draw();
  return a1;
}

void drawcut()
{
  DEBUG("start drawcut");
  // draw arrows representing the cuts in the current histogram

  // delete old arrows
  if (ArrowMin[gPadNr] != 0) {
    delete ArrowMin[gPadNr];
    ArrowMin[gPadNr] = 0;
  }
  if (ArrowMax[gPadNr] != 0) {
    delete ArrowMax[gPadNr];
    ArrowMax[gPadNr] = 0;
  }

  Int_t hstart = FindFirstHisto();
  if (hstart < 0)
    return;
  const char * t = gStack[gPadNr][hstart]->GetName();
  DEBUG("histo name: " << t);
  // only draw arrows for N-1 plots
  if (t[5] != 'n')
    return;
  DEBUG("Getting values from options file");
  // get values from option file
  Double_t xmin   = GetOptMin(gStart, t);
  Double_t xmax   = GetOptMax(gStart, t);
  if (xmin != gOptDefault)
    ArrowMin[gPadNr] = arrow(xmin);
  if (xmax != gOptDefault)
    ArrowMax[gPadNr] = arrow(xmax);
  DEBUG("end drawcut");
}

void drawoverflow()
{
  if (!gHasOverflow[gPadNr] || !gMoveOverflow)
    return;
  
  Double_t lmax = -1E99;
  Double_t bw = 0;
  Double_t xmax = 0, xmin = 0, ymin = 0;
  // find maximum in the overflow bin
  for (Int_t process = 0; process < gMaxProcess; process++) {
    TH1D * h = gStack[gPadNr][process];
    if (h == 0)
      continue;
    Int_t bin  = h->GetNbinsX();
    lmax = TMath::Max(lmax, gStack[gPadNr][process]->GetBinContent(bin));
    bw   = h->GetBinWidth(bin);
    xmax = h->GetXaxis()->GetXmax();
    xmin = h->GetXaxis()->GetXmin();
    ymin = h->GetYaxis()->GetXmin();
  }
  if (lmax < 0)
    return;
  const char * text = "Last bin includes overflow";
  // text above last bin
  // TText * t   = new TText(x+bw/2,1.1*(lmax+sqrt(lmax)), "Overflow bin");
  // t->SetTextAlign(12);
  // text next to overflow bin, right from axis
  TText * t = new TText(xmax+0.01*(xmax-xmin), ymin, text);
  t->SetTextAlign(13);
  t->SetTextSize(0.03);
  t->SetTextAngle(90);
  t->Draw();
}

void draw(Bool_t autotitle = kTRUE)
{
  DEBUG("enter draw()");

  // draw the histograms in the canvas
  Bool_t first = kTRUE;   // needed for drawing first histogram without "same"

  // clear drawing pad and delete histograms
  if (gPad) {
    gPad->cd();
    gPad->Clear();
  }
  else {
    MakeCanvas();
  }

  // enable or disable grid in x and/or y
  gPad->SetGridy((gGrid / 10) % 10);
  gPad->SetGridx(gGrid % 10);

  // now draw all stacked histograms in canvas
  TH1D * histo = 0;
  for (Int_t i = 0; i < gMaxProcess-1; i++) {
    Int_t process = gOrder[gPadNr][i];
    DEBUG("process " << process << " [" << gProcess[process].fname << "]");
    histo = gStack[gPadNr][process];
    if (histo == 0)
      continue;

    histo->SetLineStyle(gProcess[process].lstyle);
    histo->SetLineColor(gProcess[process].lcolor);

    // only draw stacked histograms in first loop
    if (!gProcess[process].stack)
      continue;
    DEBUG("stacked histogram");
    if (gIsColor) {
      histo->SetFillStyle(1001);
      histo->SetFillColor(gProcess[process].fcolor);
    }
    else {
      // delete old shadow histogram
      TH1D * shadow = gShadow[gPadNr][process];
      if (shadow) {
	delete shadow;
      }
      shadow = new TH1D(*histo);
      shadow->SetFillStyle(1001);
      shadow->SetFillColor(kWhite);
      shadow->Draw("same");
      histo->SetFillStyle(gProcess[process].hstyle);
      histo->SetFillColor(gProcess[process].hcolor);
    }
    if (first) {
      if (gProcess[process].join) {
	ERROR("first process to be drawn (" << gProcess[process].fname 
	      << ") has option 'join'");
      }
      DEBUG("draw first");
      // draw first
      if (autotitle) {
	histo->SetXTitle(GetXTitle(histo));
	histo->SetYTitle(GetYTitle(histo));
      }
      setopt(histo);
      histo->Draw("histo");
      histo->GetXaxis()->CenterTitle();
      histo->GetYaxis()->CenterTitle(kFALSE);
      first = kFALSE;
      DEBUG("end draw first");
    }
    else {
      // draw others, but only if not joined
      if (!gProcess[process].join) {
	DEBUG("draw");
	histo->Draw("histosame");
      }
    }
  }

  // draw unstacked histograms as lines, also draw data
  for (Int_t i = 0; i < gMaxProcess; i++) {
    Int_t process = gOrder[gPadNr][i];
    LOG(10, "process " << process << " [" << gProcess[process].fname << "]");
    // does the histogram exist?
    histo = gStack[gPadNr][process];
    if (histo == 0)
      continue;
    // ignore stacked histograms
    if (gProcess[process].stack)
      continue;

    if (process == gMaxProcess - 1) {
      // data
      DEBUG("data");
      histo->SetMarkerStyle(gProcess[process].marker);
      histo->SetMarkerColor(gProcess[process].mcolor);
      histo->SetMarkerSize(gProcess[process].msize);
    }
    else {
      histo->SetLineWidth(2);
    }

    DEBUG("unstacked histogram");
    if (first) {
      if (gProcess[process].join) {
	ERROR("first process to be drawn (" << gProcess[process].fname 
	      << ") has option 'join'");
      }
      // draw first
      DEBUG("draw unstacked first");
      if (autotitle) {
	histo->SetXTitle(GetXTitle(histo));
	histo->SetYTitle(GetYTitle(histo));
      }
      setopt(histo);
      if (process == gMaxProcess-1)
	histo->Draw("e1p");
      else
	histo->Draw("histo");
      histo->GetXaxis()->CenterTitle();
      histo->GetYaxis()->CenterTitle(kFALSE);
      first = kFALSE;
      DEBUG("end draw unstacked first");
    }
    else {
      if (process == gMaxProcess-1) {
	// draw data
	DEBUG("draw data");
	histo->Draw("e1psame");
      } 
      else {
	// draw others, but only if not joined
	if (!gProcess[process].join) {
	  DEBUG("draw");
	  histo->Draw("histosame");
	}
      }
    }
  }  

  // draw cut (if it is there)
  drawcut();
  drawoverflow();
  padupdate();
  DEBUG("leave draw()");
}

void max(Double_t maximum)
{
  // set maximum for all histograms
  FindFirstHisto();
  for (Int_t i = 0; i < gMaxProcess; i++) {
    if (gStack[gPadNr][i] != 0)
      gStack[gPadNr][i]->SetMaximum(maximum);
  }
  padupdate();
}

void min(Double_t minimum)
{
  DEBUG("enter min()");
  // set minimum for all histograms
  if (gPad && gPad->GetLogy() && minimum <= 0) {
    ERROR("Can not set minimum " << minimum << " with logarithmic y-axis");
    return;
  }
  FindFirstHisto();
  for (Int_t i = 0; i < gMaxProcess; i++) {
    if (gStack[gPadNr][i] != 0)
      gStack[gPadNr][i]->SetMinimum(minimum);
  }
  padupdate();
  DEBUG("leave min()");
}

void findbins(Double_t xlow, Double_t xup, Int_t & startbin, Int_t & endbin)
{
  DEBUG("enter findbins");
  // returns the bin numbers to the corresponding x-axis values. If xlow==xup,
  // then the full range is returned

  // find first histo
  Int_t hstart = FindFirstHisto();
  if (hstart < 0) {
    startbin = 1;
    endbin = 1;
    return;
  }

  // full range?
  if (xlow == xup) {
    startbin = 1;
    endbin   = gStack[gPadNr][hstart]->GetNbinsX();
  }
  else {
    startbin = gStack[gPadNr][hstart]->GetXaxis()->FindFixBin(xlow);
    endbin   = gStack[gPadNr][hstart]->GetXaxis()->FindFixBin(xup);
  }
}

void findmax(Double_t low = 0, Double_t up = 0)
{
  DEBUG("enter findmax()");
  // find maximum by iterating over all histograms and counting bin content

  // convert x-axis range to bin range
  Int_t start, end;
  findbins(low, up, start, end);

  Double_t lmax = -1E99;
  for (Int_t process = 0; process < gMaxProcess; process++) {
    if (gStack[gPadNr][process]) {
      for (Int_t i = start; i <= end; i++) {
	Double_t temp = gStack[gPadNr][process]->GetBinContent(i);
	if (temp > lmax)
	  lmax = temp;
      }
    }
  }
  if (lmax > 0)
    max(lmax*(1+1/TMath::Sqrt(lmax))*1.1);
  DEBUG("leave findmax()");
}

void findmin(Double_t low = 0, Double_t up = 0)
{
  DEBUG("enter findmin()");
  // find minimum by iterating over all histograms and counting bin content

  // if normal axis, set minimum to zero
  if (gPad->GetLogy() == 0) {
    min(0);
    return;
  }

  // convert x-axis range to bin range
  Int_t start, end;
  findbins(low, up, start, end);


  Double_t lmin = 1E99;
  for (Int_t process = 0; process < gMaxProcess-1; process++) {
    if (gStack[gPadNr][process]) {
      for (Int_t i = start; i <= end; i++) {
	Double_t temp = gStack[gPadNr][process]->GetBinContent(i);
	if (temp > 0 && temp < lmin)
	  lmin = temp;
      }
    }
  }

  // some space at bottom
  lmin *= 0.5;
  // adjust to normal ranges
  if (lmin < 0.001) {
    lmin = 0.001;
  }

  min(lmin);
  DEBUG("leaved findmin()");
}

void zoom(Double_t low, Double_t up)
{
  // zoom on x-axis
  Double_t overflow = 0.;
  for (Int_t process = 0; process < gMaxProcess-1; process++) {
    TH1D * h = gStack[gPadNr][process];
    if (h) {
      h->SetAxisRange(low, up);
      // count overflows in bins that are not shown
      Int_t maxbin = h->GetXaxis()->FindFixBin(up);
      for (Int_t i = maxbin+1; i <= h->GetNbinsX()+1; i++) {
	overflow += h->GetBinContent(i);
      }
    }
  }
  // find out if now the histogram contains overflows
  gHasOverflow[gPadNr] = overflow > 0.;
  padupdate();
}

void unzoom()
{
  // unzoom x-axis
  Double_t overflow = 0.;
  for (Int_t process = 0; process < gMaxProcess-1; process++) {
    TH1D * h = gStack[gPadNr][process];
    if (h) {
      TAxis * axis = h->GetXaxis();
      if (axis)
	axis->UnZoom();
      // count overflows in bins that are not shown
      overflow += h->GetBinContent(h->GetNbinsX()+1);
    }
  }
  // find out if now the histogram contains overflows
  gHasOverflow[gPadNr] = overflow > 0.;
  padupdate();
}

void rebin(Int_t nbins)
{
  // rebin histograms: combine nbins bins into one bin
  if (FindFirstHisto() < 0)
    return;

  // rebin all histograms - stacked and unstacked
  for (Int_t i = 0; i < gMaxProcess; i++) {
    if (gHisto[gPadNr][i])
      gHisto[gPadNr][i]->Rebin(nbins);
    if (gStack[gPadNr][i])
      gStack[gPadNr][i]->Rebin(nbins);
  }

  // adjust maximum after a rebin and draw
  findmax();
  findmin();
  draw();
}

void shiftbin(Int_t nbins)
{
  // shift whole MC distribution by nbins bins
  if (nbins > 0) {
    for (Int_t process = 0; process < gMaxProcess-1; process++) {
      // unstacked histograms
      TH1D * histo = gHisto[gPadNr][process];
      if (histo == 0)
	continue;
      for (Int_t bin = histo->GetNbinsX(); bin >= 0; bin--) {
	histo->SetBinContent(bin, histo->GetBinContent(bin-nbins));
      }
      // stacked histograms
      histo = gStack[gPadNr][process];
      if (histo == 0)
	continue;
      for (Int_t bin = histo->GetNbinsX(); bin >= 0; bin--) {
	histo->SetBinContent(bin, histo->GetBinContent(bin-nbins));
      }
    }
  }
  else {
    for (Int_t process = 0; process < gMaxProcess-1; process++) {
      // unstacked histograms
      TH1D * histo = gHisto[gPadNr][process];
      if (histo == 0)
	continue;
      for (Int_t bin = 0; bin <= histo->GetNbinsX(); bin++) {
	if (histo)
	  histo->SetBinContent(bin, histo->GetBinContent(bin-nbins));
      }
      // stacked histograms
      histo = gStack[gPadNr][process];
      if (histo == 0)
	continue;
      for (Int_t bin = 0; bin <= histo->GetNbinsX(); bin++) {
	if (histo)
	  histo->SetBinContent(bin, histo->GetBinContent(bin-nbins));
      }
    }
  }
  draw();
}

void mirror()
{
  // mirror all histograms
  for (Int_t process = 0; process < gMaxProcess; process++) {
    // unstacked histograms
    TH1D * histo = gHisto[gPadNr][process];
    if (histo == 0)
      continue;
    TH1D * hnew  = (TH1D *) histo->Clone();
    hnew->Reset();
    Int_t nbinmax = histo->GetNbinsX()+1;
    for (Int_t nbin = 0; nbin <= nbinmax; nbin++) {
      hnew->SetBinContent(nbin, histo->GetBinContent(nbinmax-nbin));
    }
    delete histo;
    gHisto[gPadNr][process] = hnew;
    // stacked histograms
    histo = gStack[gPadNr][process];
    if (histo == 0)
      continue;
    hnew  = (TH1D *) histo->Clone();
    hnew->Reset();
    nbinmax = histo->GetNbinsX()+1;
    for (Int_t nbin = 0; nbin <= nbinmax; nbin++) {
      hnew->SetBinContent(nbin, histo->GetBinContent(nbinmax-nbin));
    }
    delete histo;
    gStack[gPadNr][process] = hnew;
  }
  draw();
}

void print(const char * hname)
{
  // print current canvas in a file "name", substitute current
  // histogram name if no name given
  if (FindFirstHisto() < 0)
    return;

  char fpath[256];
  if (hname == 0) {
    TH1D * h1 = gStack[gPadNr][0];
    if (h1)
      if (gStart != gEnd)
	snprintf(fpath, 256, "%s-%s-%s.pdf",
		h1->GetName(), gPeriod[gStart], gPeriod[gEnd]);
      else
	snprintf(fpath, 256, "%s-%s.pdf", h1->GetName(), gPeriod[gStart]);
    else
      snprintf(fpath, 256, "%s.pdf", "unknown");
  }
  else
    snprintf(fpath, 256, "%s", hname);

  gCanvas->Print(fpath);
}

void pprint(const char * hname)
{
  // print current pad in a file "name", substitute current
  // histogram name if no name given
  if (FindFirstHisto() < 0) {
    if (hname != 0)
      gPad->Print(hname);
    return;
  }

  char fpath[256];
  if (hname == 0) {
    TH1D * h1 = gStack[gPadNr][0];
    if (h1)
      if (gStart != gEnd)
	snprintf(fpath, 256, "%s-%s-%s.pdf",
		h1->GetName(), gPeriod[gStart], gPeriod[gEnd]);
      else
	snprintf(fpath, 256, "%s-%s.pdf", h1->GetName(), gPeriod[gStart]);
    else
      snprintf(fpath, 256, "%s.pdf", "unknown");
  }
  else
    snprintf(fpath, 256, "%s", hname);

  gPad->Print(fpath);
}

/**
 * Set the axis titles for this plot.
 *
 * If no name is given, then the default title is taken
 */
void title(const char * title)
{
  for (Int_t process = 0; process < gMaxProcess-1; process++) {
    TH1D * h = gStack[gPadNr][process];
    if (h) {
      h->SetTitle(title);
    }
  }
  DEBUG("title was set");
  // now use automatic title
  draw(kTRUE);
}

/**
 * Set y-axis to have a logarithmic scale
 */
void logy()
{
  DEBUG("enter logy()");
  if (gCanvas == 0)
    MakeCanvas();
  min(0.001);
  gPad->SetLogy(1);
  draw();
  DEBUG("leave logy()");
}

/**
 * Set y-axis to have a linear scale
 */
void liny()
{
  // set y axis to be linear
  if (gCanvas == 0)
    MakeCanvas();
  gPad->SetLogy(0);
  min(0.);
  padupdate();
}

void drawperiod(Int_t posi = -1)
{
  DEBUG("enter drawperiod()");

  // draw energy in current plot
  static TLatex * t[gMaxPad] = { 0 };
  char etext[50];
  if (gStart != gEnd) {
    snprintf(etext, 50, "CMS %s - %s Work in progress", gPeriod[gStart], gPeriod[gEnd]); //preliminary
  }
  else {
    snprintf(etext, 50, "CMS %s Work in progress", gPeriod[gStart]);
  }

  //Int_t hstart = FindFirstHisto();
  //if (hstart < 0)
  //  return;
  if (t[gPadNr] != 0) {
    delete t[gPadNr];
    t[gPadNr] = 0;
  }
  if (posi > 0) {
    t[gPadNr] = new TLatex(0.955, 0.99, etext);
    t[gPadNr]->SetTextAlign(33);
  }
  else if (posi < 0) {
    t[gPadNr] = new TLatex(0.24, 0.99, etext);
    t[gPadNr]->SetTextAlign(13);
  }
  else {
    t[gPadNr] = new TLatex(0.5, 0.99, etext);
    t[gPadNr]->SetTextAlign(23);
  }
  t[gPadNr]->SetTextSize(0.04);
  t[gPadNr]->SetNDC();
  t[gPadNr]->Draw();
}

void compose()
{
  DEBUG("enter compose");
  // sum up the histograms in the current order, such that they can be
  // drawn on the screen, one before the other

  // loop over processes in reverse order and sum up
  for (Int_t i = gMaxProcess-2; i >= 0; i--) {
    Int_t process = gOrder[gPadNr][i];
    DEBUG("process " << gProcess[process].fname);
    TH1D * histo = gHisto[gPadNr][process];
    // not existing
    if (histo == 0)
      continue;
    // copy histogram to stack
    DEBUG("copy histogram to stack");
    gStack[gPadNr][process] = new TH1D(*histo);
    // not to be stacked - copy only
    if (!gProcess[process].stack) {
      continue;
    }
    // add sum of previous ones
    for (Int_t j = i+1; j < gMaxProcess-1; j++) {
      Int_t backmc = gOrder[gPadNr][j];
      // not stacked - ignore
      if (!gProcess[backmc].stack)
	continue;
      // not existing
      if (gHisto[gPadNr][backmc] == 0)
	continue;
      // Add
      DEBUG("Adding histogram " << gProcess[backmc].fname);
      gStack[gPadNr][process]->Add(gStack[gPadNr][backmc]);
      // break loop, since all previous were added
      break;
    }
  }
  TH1D * hdata = gHisto[gPadNr][gMaxProcess-1];
  if (hdata)
    gStack[gPadNr][gMaxProcess-1] = new TH1D(*hdata);
  DEBUG("leave compose");
}

void legend(Double_t mincontent, Int_t posi, Double_t miny)
{
  DEBUG("enter legend()");
  // draw a legend with currently used colors
  static TLegend * t[gMaxPad] = { 0 };
  Int_t position = FindFirstHisto();
  if (position < 0)
    return;

  if (miny < 0) {
    miny = 0.2;
  }

  // get sum of all processes
  Double_t integral = gStack[gPadNr][position]->Integral();
  DEBUG("integral = " << integral);

  // find processes that contribute with given fraction to the integral
  vector<Int_t> entries;
  for (Int_t i = 0; i < gMaxProcess-1; i++) {
    Int_t process = gOrder[gPadNr][i];
    // not existing
    if (gHisto[gPadNr][process] == 0)
      continue;
    // if joined, it was already considered in previous iteration
    if (gProcess[process].join)
      continue;
    // OK, get integral
    Double_t sum = gHisto[gPadNr][process]->Integral();
    DEBUG(gProcess[process].fname << ": " << sum);
    // check if histograms are joined -> we need to add statistics
    for (Int_t j = i+1; j < gMaxProcess-1; j++) {
      Int_t proc = gOrder[gPadNr][j];
      // only add joined histograms
      if (!gProcess[proc].join)
	break;
      // not existing
      if (gHisto[gPadNr][proc] == 0)
	continue;
      DEBUG(gProcess[proc].fname << ": + " << gHisto[gPadNr][proc]->Integral());
      sum += gHisto[gPadNr][proc]->Integral();
    }
    DEBUG(gProcess[process].fname << ": " << sum);
    if (sum/integral < mincontent) {
      continue;
    }
    entries.push_back(process);
    sum = 0;
  }
  // add data
  if (gHisto[gPadNr][gMaxProcess-1])
    entries.push_back(gMaxProcess-1);

  // determine minimum
  UInt_t num = entries.size();
  miny = TMath::Max(0.92-num*0.06, miny);
  DEBUG("Found " << num << " histograms, miny = " << miny);

  if (t[gPadNr]) {
    delete t[gPadNr];
    t[gPadNr] = 0;
  }
  Double_t minx, maxx;
  if (posi > 0) {
    minx = 0.720; 
    maxx = 0.92;
  }
  else if (posi < 0) {
    minx = 0.173;
    maxx = 0.413;
  }
  else {
    minx = 0.45;
    maxx = 0.69;
  }
  t[gPadNr] = new TLegend(minx, miny, maxx, 0.92);
  setopt(t[gPadNr]);
  for (Int_t i = 0; i < (Int_t) entries.size(); i++) {
    Int_t process = entries[i];
    const char * opt;
    const char * name = gProcess[process].tname;
    if (process == gMaxProcess-1) {
      if (gGerman)
	name = "Daten";
      opt= "p";
    }
    else if (gProcess[process].stack) {
      opt = "f";
    }
    else {
      opt = "l";
    }
    t[gPadNr]->AddEntry(gStack[gPadNr][process], name, opt);
  }
  t[gPadNr]->Draw();

  // insert energy information
  drawperiod(posi != 0 ? -posi : -1);
  DEBUG("exit legend()");
}

// void lumi()
// {
//   const char * name[50];
//   sprintf(name, "#int L dt = %2.1f pb^{-1}", lumi);
//   t[gPadNr]->DrawLatex(0.4,0.92,name);
//}


void subfigure(const char * subfig)
{
  // add entry for subfigure in the current canvas, e.g. "a)"
  TText * t = new TText(0.95, 0.97, subfig);
  t->SetTextAlign(33);
  t->SetNDC();
  t->Draw();
}

void DeleteOld()
{
  // delete existing histograms gStack[gPadNr][] and gHisto[gPadNr][]
  for (Int_t i = 0; i < gMaxProcess; i++) {
    // unstacked
    if (gHisto[gPadNr][i] != 0) {
      delete gHisto[gPadNr][i];
      // mark as deleted
      gHisto[gPadNr][i] = 0;
    }
    // stacked
    if (gStack[gPadNr][i] != 0) {
      delete gStack[gPadNr][i];
      // mark as deleted
      gStack[gPadNr][i] = 0;
    }
  }
}

// TH1 * get_period(Int_t process, const char * hname)
// {
//   // loops over currently selected periods and adds the histogram
//   // "hname" from the corresponding mc, returning the summary histogram
//   char fpath[400];

//   TH1 * htemp = 0;
//   TH1 * hadd  = 0;

//   for (Int_t period = gStart; period <= gEnd; period++) {
//     // open file
//     snprintf(fpath, 400, "%s/%s.root", getpath(period), gProcess[process].fname);
//     DEBUG("Opening file " << fpath);
//     TFile * f = TFile::Open(fpath);
//     if (f == 0 || !f->IsOpen()) {
//       ERROR("Could not open file " << fpath);
//       if (hadd)
// 	delete hadd;
//       return 0;
//     }
//     // get key list from file
//     TList * keylist = f->GetListOfKeys();
//     TIter next(keylist);
//     while (TKey * key = (TKey *) next()) {
//       // look for histogram name
//       TString s(key->GetName());
//       if (s != hname)
// 	continue;
//       DEBUG(key->GetClassName() << ": " << key->GetName());
//       if (strcmp(key->GetClassName(), "TH1D") && 
// 	  strcmp(key->GetClassName(), "TH2D") && 
// 	  strcmp(key->GetClassName(), "TH3D")) {
// 	ERROR("Can not work with class " << key->GetClassName());
// 	if (hadd)
// 	  delete hadd;
// 	delete f;
// 	return 0;
//       }
//       htemp = (TH1 *) key->ReadObj();
//       // no need to look further in key list
//       break;
//     }

//     // scale mc according to luminosity
//     if (process != gMaxProcess-1) { // do not scale data
//       Double_t mcLumi = gProcessInfo[period][process].Nev/gProcessInfo[period][process].xs;
//       Double_t weight = gLumi[period]/mcLumi * gProcessInfo[period][process].weight;
//       DEBUG("mcLumi = " << mcLumi << ", weight = " << weight);
//       htemp->Scale(weight);
//     }

//     // add to summary histogram
//     if (hadd == 0) {
//       hadd = new TH1(*htemp);
//       hadd->SetDirectory(0);
//     }
//     else {
//       hadd->Add(htemp);
//     }
//     // close file
//     delete f;
//   }
//   return hadd;
// }

TH1D * addperiod(Int_t process, const char * hname,
		 const char * selection, Int_t nbins, Double_t min,
		 Double_t max)
{
  // loops over currently selected periods and adds the histogram
  // "hname" from the corresponding mc, returning the summary histogram
  char fpath[400];
  char mname[400];

  TKey * key   = 0; // key from file
  TH1D * htemp = 0;
  TH1D * hadd  = 0;

  for (Int_t period = gStart; period <= gEnd; period++) {
    // open file
    snprintf(fpath, 400, "%s/%s.root", getpath(period), gProcess[process].fname);
    DEBUG("Opening file " << fpath);
    TFile * f = TFile::Open(fpath);
    if (f == 0 || !f->IsOpen()) {
      if (hadd)
	delete hadd;
      ERROR("Could not open file " << fpath);
      return 0;
    }
    // get key from file
    char histname[strlen(hname)+4];
    if (gVersion == 0) {
      sprintf(histname, "h1_%s", hname);
    }
    else if (gVersion == 1) {
      sprintf(histname, "h1_%s_%s", gStage.Data(), hname);
    }
    key = (TKey *) f->GetKey(histname);
    DEBUG("Histogram key = " << key);
    if (key == 0) {
      // key has not been found, generate histogram from tree
      // create specified histogram
      TH1D * hmine = new TH1D("hmine", hname, nbins, min, max);
      // get tree from current file
      const char * treename = gConfig->GetValue("settings.treename", "t");
      TTree * dtree = (TTree *) f->Get(treename);
      if (dtree == 0) {
	// no tree found, try next period
	ERROR("Tree \"" << treename << "\" not found in file " << fpath);
	delete f;
	continue;
      }
      // save tree data in my histogram
      snprintf(mname, 400, "%s>>hmine", hname);
      Int_t result = dtree->Draw(mname, selection, "goff");
      if (result < 0) {
	ERROR("Variable not found in tree");
	delete hmine;
	delete f;
	return 0;
      }
      htemp = (TH1D *) hmine;
      htemp->SetDirectory(0);
    }
    else {
      if (strcmp(key->GetClassName(), "TH1D")) {
	ERROR("Can only plot one-dimensional histograms");
	if (hadd)
	  delete hadd;
	delete f;
	return 0;
      }
      htemp = (TH1D * ) f->Get(histname);
    }
    // check physical minimum
    if (htemp->GetMinimum() < 0) {
      ERROR("Minimum " << htemp->GetMinimum() << " < 0 in file "<< fpath);
    }

    // move overflow bin
    Int_t bin = htemp->GetNbinsX();
    if (htemp->GetBinContent(bin+1) > 0) {
      gHasOverflow[gPadNr] = kTRUE;
      if (gMoveOverflow) {
	htemp->SetBinContent(bin, htemp->GetBinContent(bin)+htemp->GetBinContent(bin+1));
	htemp->SetBinError(bin, TMath::Sqrt(TMath::Power(htemp->GetBinError(bin), 2) +
					    TMath::Power(htemp->GetBinError(bin+1), 2)));
      }
    }

    // scale mc according to luminosity
    if (process != gMaxProcess-1) { // do not scale data
      Double_t mcLumi = gProcessInfo[period][process].Nev/gProcessInfo[period][process].xs;
      Double_t weight = gLumi[period]/mcLumi * gProcessInfo[period][process].weight;
      DEBUG("mcLumi = " << mcLumi << ", weight = " << weight);
      htemp->Scale(weight);
    }

    // add to summary histogram
    if (hadd == 0) {
      hadd = new TH1D(*htemp);
      hadd->SetDirectory(0);
    }
    else {
      hadd->Add(htemp);
    }
    // close file
    delete f;
  }
  return hadd;
}

void print_order()
{
  for (Int_t i = 0; i < gMaxProcess; i++) {
    INFO("gOrder[" << i << "] = " << gOrder[gPadNr][i]);
  }
}

void auto_order()
{
  DEBUG("Enter auto_order()");
  Double_t integral[gMaxProcess];
  // order all but data
  for (Int_t i = gMaxSignal; i < gMaxProcess-1; i++) {
    integral[i] = 0;
    if (!gHisto[gPadNr][i]) {
      continue;
    }
    // is it already joined? 
    if (gProcess[i].join && i > 0) {
      // if it is a joined histogram, keep it below its parent
      integral[i] = integral[i-1]*0.999999;
    }
    else {
      // otherwise order according to its integral.
      Double_t sum = gHisto[gPadNr][i]->Integral();
      // loop over all histograms which are joined with this one
      // and add the integral
      for (Int_t j = i+1; j < gMaxProcess-1; j++) {
	// stop if no more bg's to be added
	if (!gProcess[j].join)
	  break;
	if (gHisto[gPadNr][j] == 0)
	  continue;
	sum += gHisto[gPadNr][j]->Integral();
      }
      // set integral to a small default value keeping a meaningful order if
      // all histograms are zero
      integral[i] = sum+(gMaxProcess-i)*1e-10;
    }
    DEBUG(gProcess[i].fname << ": " << integral[i]);
  }
  TMath::Sort(gMaxProcess-1-gMaxSignal, integral, gOrder[gPadNr]);
  DEBUG("leave auto_order()");
}

void plot(const char * hname, const char * selection,
	  Int_t nbins, Double_t min, Double_t max)
{
  // plot histogram "hname" for the selected selection and period
  // loops over all files.

  // check for canvas
  if (gCanvas == 0) {
    MakeCanvas();
  }

  // delete previously plotted histograms
  DeleteOld();

  // reset overflow display
  gHasOverflow[gPadNr] = kFALSE;

  // loop over all processes
  for (Int_t process = 0; process < gMaxProcess; process++) {
    // get summary histogram from all periods
    gHisto[gPadNr][process] = addperiod(process, hname,
					selection, nbins, min, max);
  }
  // order histograms according to integral
  if (gAutoOrder)
    auto_order();
  // compose histograms for drawing
  compose();
  // find maximum & minimum
  findmax();
  findmin();
  draw();
}

void plotadd(const char * name1, const char * name2)
{
  // add up two different plots
  // check for canvas
  if (gCanvas == 0) {
    MakeCanvas();
  }

  // delete previously plotted histograms
  DeleteOld();

  // loop over all processes
  for (Int_t process = 0; process < gMaxProcess; process++) {
    // get summary histogram from all periods
    gHisto[gPadNr][process] = addperiod(process, name1,
					0, 1, 0, 1);
    gHisto[gPadNr][process]->Add(addperiod(process, name2,
					   0, 1, 0, 1));
  }

  // compose histograms for drawing
  compose();
  // find maximum & minimum
  findmax();
  findmin();
  draw();
}

void plotadd(const char * name1, const char * name2,
	     const char * name3, const char * name4)
{
  // add up two different plots
  // check for canvas
  if (gCanvas == 0) {
    MakeCanvas();
  }

  // delete previously plotted histograms
  DeleteOld();

  // loop over all processes
  for (Int_t process = 0; process < gMaxProcess; process++) {
    // get summary histogram from all periods
    gHisto[gPadNr][process] = addperiod(process, name1,
					0, 1, 0, 1);
    gHisto[gPadNr][process]->Add(addperiod(process, name2,
					   0, 1, 0, 1));
    gHisto[gPadNr][process]->Add(addperiod(process, name3,
					   0, 1, 0, 1));
    gHisto[gPadNr][process]->Add(addperiod(process, name4,
					   0, 1, 0, 1));
  }

  // compose histograms for drawing
  compose();
  // find maximum & minimum
  findmax();
  findmin();
  draw();
}

TH2D * addperiod2(Int_t process, const char * hname,
		  const TH2D * hsample, const char * selection)
{
  // loops over currently selected periods and adds the histogram
  // "hname" from the corresponding mc, returning the summary histogram
  char fpath[400];
  char mname[100];

  TKey * key   = 0; // key from file
  TH2D * htemp = 0;
  TH2D * hadd  = 0;

  for (Int_t period = gStart; period <= gEnd; period++) {
    // open file
    snprintf(fpath, 400, "%s/%s.root", getpath(period), gProcess[process].fname);
    TFile * f = new TFile(fpath, "READ");
    if (!f->IsOpen()) {
      ERROR("File " << fpath << " not found.");
      return 0;
    }
    // get key from file
    char histname[strlen(hname)+4];
    if (gVersion == 0) {
      sprintf(histname, "h2_%s", hname);
    }
    else if (gVersion == 1) {
      sprintf(histname, "h2_%s_%s", gStage.Data(), hname);
    }
    key = (TKey *) f->GetKey(histname);
    DEBUG("Histogram key = " << key);
    if (key == 0) {
      // key has not been found, generate histogram from tree
      // create specified histogram
      // get tree from current file
      const char * treename = gConfig->GetValue("treename", "t");
      TTree * dtree = (TTree *) f->Get(treename);
      if (dtree == 0) {
	// no tree found, try next period
	WARNING("Tree \"" << treename << "\" not found in file" << fpath);
	continue;
      }
      if (!hsample) {
	ERROR("No sample histogram given!");
	delete f;
	return 0;
      }
      // copy structure from sample histogram
      htemp = new TH2D(*hsample);
      htemp->SetName("htemp");
      // save tree data in my histogram
      snprintf(mname, 100, "%s>>htemp", hname);
      dtree->Draw(mname, selection, "goff");
      htemp->SetTitle(hname);
      htemp->SetDirectory(0);
    }
    else {
      if (strcmp(key->GetClassName(), "TH2D")) {
	ERROR("Can only plot two-dimensional histograms");
	delete f;
	return 0;
      }
      htemp = (TH2D * ) f->Get(histname);
    }
    // check physical minimum
    if (htemp->GetMinimum() < 0) {
      ERROR("Minimum " << htemp->GetMinimum() << " < 0 in file " << fpath);
    }

    // scale mc according to luminosity
    if (process != gMaxProcess-1) { // do not scale data
      Double_t mcLumi = gProcessInfo[period][process].Nev/gProcessInfo[period][process].xs;
      Double_t weight = gLumi[period]/mcLumi * gProcessInfo[period][process].weight;
      htemp->Scale(weight);
    }

    // add to summary histogram
    if (hadd == 0) {
      hadd = new TH2D(*htemp);
      hadd->SetDirectory(0);
    }
    else {
      hadd->Add(htemp);
    }
    // close file
    delete f;
  }
  return hadd;
}

void plot2(const char * hname, const TH2D * hsample, const char * selection,
	   const char * drawopt)
{
  // plot a two-dimensional plot in a new Canvas
  static TCanvas * gCanvas2 = 0;

  Int_t nMC = 0;
  // creation loop
  for (Int_t process = 0; process < gMaxProcess; process++) {
    // get summary histogram from all periods
    if (gHisto2[process]) {
      delete gHisto2[process];
      gHisto2[process] = 0;
    }
    gHisto2[process] = addperiod2(process, hname, hsample, selection);
    if (gHisto2[process] == 0) {
      ERROR("Histogram " << hname << " not found for process " << process);
      continue;
    }
    if (gHisto2[process]->Integral() == 0) {
      DEBUG("Empty histogram found for process " <<  process);
      continue;
    }
    nMC++;
    gHisto2[process]->SetXTitle(GetXTitle((TH1D *) gHisto2[process]));
    gHisto2[process]->SetYTitle(GetYTitle((TH1D *) gHisto2[process]));
    gHisto2[process]->SetTitle(gProcess[process].tname);
  }
  if (nMC < 1) {
    ERROR("Only empty histograms found");
    return;
  }

  gStyle->SetOptTitle(1);
  if (gCanvas2 == 0) {
    gCanvas2 = new TCanvas("gCanvas2", "two-dimensional plots",
			   10, 10, 600, (Int_t) (600*TMath::Sqrt(2.0)));
    Int_t dy = (Int_t) TMath::Sqrt((Double_t) nMC);
    INFO("nMC = " << nMC << ", dy = " << dy);
    Int_t dx = nMC / dy;
    if (dx*dy < nMC)
      dy++;
    gCanvas2->Divide(dx, dy);
  }

  // draw loop
  Int_t pad = 1;
  for (Int_t process = 0; process < gMaxProcess; process++) {
    // get summary histogram from all periods
    if (gHisto2[process] == 0)
      continue;
    if (gHisto2[process]->Integral() == 0)
      continue;
    gCanvas2->cd(pad);
    gHisto2[process]->Draw(drawopt);
    pad++;
  }
}

TH3D * addperiod3(Int_t process, const char * hname,
		  const TH3D * hsample, const char * selection)
{
  // loops over currently selected periods and adds the histogram
  // "hname" from the corresponding mc, returning the summary histogram
  char fpath[400];
  char mname[100];

  TKey * key   = 0; // key from file
  TH3D * htemp = 0;
  TH3D * hadd  = 0;

  for (Int_t period = gStart; period <= gEnd; period++) {
    // open file
    snprintf(fpath, 400, "%s/%s.root", getpath(period), gProcess[process].fname);
    TFile * f = new TFile(fpath, "READ");
    if (!f->IsOpen()) {
      ERROR("File " << fpath << " not found.");
      return 0;
    }
    // get key from file
    char histname[strlen(hname)+4];
    if (gVersion == 0) {
      sprintf(histname, "h3_%s", hname);
    }
    else if (gVersion == 1) {
      sprintf(histname, "h3_%s_%s", gStage.Data(), hname);
    }
    key = (TKey *) f->GetKey(histname);
    DEBUG("Histogram key = " << key);
    if (key == 0) {
      // key has not been found, generate histogram from tree
      // create specified histogram
      // get tree from current file
      const char * treename = gConfig->GetValue("settings.treename", "t");
      TTree * dtree = (TTree *) f->Get(treename);
      if (dtree == 0) {
	// no tree found, try next period
	WARNING("Tree \"" << treename << "\" not found in file" << fpath);
	continue;
      }
      if (!hsample) {
	ERROR("No sample histogram given!");
	delete f;
	return 0;
      }
      // copy structure from sample histogram
      htemp = new TH3D(*hsample);
      htemp->SetName("htemp");
      // save tree data in my histogram
      snprintf(mname, 100, "%s>>htemp", hname);
      dtree->Draw(mname, selection, "goff");
      htemp->SetTitle(hname);
      htemp->SetDirectory(0);
    }
    else {
      if (strcmp(key->GetClassName(), "TH3D")) {
	ERROR("Can only plot two-dimensional histograms");
	delete f;
	return 0;
      }
      htemp = (TH3D * ) f->Get(histname);
    }
    // check physical minimum
    if (htemp->GetMinimum() < 0) {
      ERROR("Minimum " << htemp->GetMinimum() << " < 0 in file " << fpath);
    }

    // scale mc according to luminosity
    if (process != gMaxProcess-1) { // do not scale data
      Double_t mcLumi = gProcessInfo[period][process].Nev/gProcessInfo[period][process].xs;
      Double_t weight = gLumi[period]/mcLumi * gProcessInfo[period][process].weight;
      htemp->Scale(weight);
    }

    // add to summary histogram
    if (hadd == 0) {
      hadd = new TH3D(*htemp);
      hadd->SetDirectory(0);
    }
    else {
      hadd->Add(htemp);
    }
    // close file
    delete f;
  }
  return hadd;
}

void plot3(const char * hname, const TH3D * hsample, const char * selection)
{
  // creation loop
  for (Int_t process = 0; process < gMaxProcess; process++) {
    // get summary histogram from all periods
    if (gHisto3[gPadNr][process]) {
      delete gHisto3[gPadNr][process];
      gHisto3[gPadNr][process] = 0;
    }
    gHisto3[gPadNr][process] = addperiod3(process, hname, hsample, selection);
    if (gHisto3[gPadNr][process] == 0) {
      ERROR("Histogram " << hname << " not found for process " << process);
      continue;
    }
    gHisto3[gPadNr][process]->SetTitle(gProcess[process].tname);
  }
}

TH1D * getHisto(const char * name, bool match, int min, int max)
{
  // sum up all histograms that match the given name, starting from min up to
  // (but not including) max.
  // 
  // If match is true all individual signal histograms whose file name
  // contains "name" are summed up. If match is false, all individual signal
  // histograms but those whose file name matches "name" are taken. If name ==
  // 0, all histograms are taken.  
  // 
  // Can be used e.g. for fitting (see xsection.C)

  TH1D * hMC = 0;
  // loop over all processes in given range
  for (Int_t j = min; j < max; j++) {
    TH1D * histo = gHisto[gPadNr][j];
    if (histo == 0)
      continue;
    TString s(gProcess[j].fname);
    if (name != 0 && ((match && !s.Contains(name)) || (!match && s.Contains(name)))) {
      DEBUG("skipping process " << s);
      continue;
    }
    DEBUG("adding process " << s);
    if (hMC == 0) {
      hMC = new TH1D(*histo);
      hMC->SetDirectory(0);
    }
    else {
      hMC->Add(histo, 1);
    }
  }
  return hMC;
}

TH2D * getHisto2(const char * name, bool match, int min, int max)
{
  // sum up all histograms that match the given name, starting from min up to
  // (but not including) max.
  // 
  // If match is true all individual signal histograms whose file name
  // contains "name" are summed up. If match is false, all individual signal
  // histograms but those whose file name matches "name" are taken. If name ==
  // 0, all histograms are taken.  
  // 
  // Can be used e.g. for fitting (see xsection.C)

  TH2D * hMC = 0;
  // loop over all processes in given range
  for (Int_t j = min; j < max; j++) {
    TH2D * histo = gHisto2[j];
    if (histo == 0)
      continue;
    TString s(gProcess[j].fname);
    if (name != 0 && ((match && !s.Contains(name)) || (!match && s.Contains(name)))) {
      DEBUG("skipping process " << s);
      continue;
    }
    DEBUG("adding process " << s);
    if (hMC == 0) {
      hMC = new TH2D(*histo);
      hMC->SetDirectory(0);
    }
    else {
      hMC->Add(histo, 1);
    }
  }
  return hMC;
}

TH1D * signalHisto(const char * name, bool match)
{
  // Return a histogram containing only signal. If match is true all
  // individual signal histograms whose file name contains "name" are summed
  // up. If match is false, all individual signal histograms but those whose
  // file name matches "name" are taken. If name == 0, all histograms are
  // taken.
  //
  // Can be used e.g. for fitting (see xsection.C)

  return getHisto(name, match, 0, gMaxSignal);
}

TH1D * backgroundHisto(const char * name, bool match)
{
  // Return a histogram containing only background. If match is true all
  // individual background histograms whose file name contains "name" are
  // summed up. If match is false, all individual background histograms but
  // those whose file name matches "name" are taken. If name == 0, all
  // histograms are taken.
  //
  // Can be used e.g. for fitting (see xsection.C)

  return getHisto(name, match, gMaxSignal, gMaxProcess-1);
}

TH2D * backgroundHisto2(const char * name, bool match)
{
  // same as above, but for 2D
  return getHisto2(name, match, gMaxSignal, gMaxProcess-1);
}

TH1D * dataHisto()
{
  // creates a histogram containing only data, extracted from global
  // histograms gStack[gPadNr][]. Can be used e.g. for fitting (see xsection.C)
  TH1D * histo = gHisto[gPadNr][gMaxProcess-1];
  if (histo == 0)
    return 0;
  TH1D * hData = new TH1D(*histo);
  hData->SetDirectory(0);
  return hData;
}

TH2D * dataHisto2()
{
  // creates a histogram containing only data, extracted from global
  // histograms gStack[gPadNr][]. Can be used e.g. for fitting (see xsection.C)
  TH2D * histo = gHisto2[gMaxProcess-1];
  if (histo == 0)
    return 0;
  TH2D * hData = new TH2D(*histo);
  hData->SetDirectory(0);
  return hData;
}

void smooth(TH1D * histo, Double_t xlow, Double_t xup, Int_t nbins)
{
  // smooth given histogram in range
  TCanvas * c2 = new TCanvas("c2", "c2", 40, 40, 800, 600);
  c2->cd();
  Int_t startbin;
  Int_t lastbin;
  findbins(xlow, xup, startbin, lastbin);
  // make a copy
  TH1D * htemp = new TH1D(*histo);
  htemp->Reset();
  INFO("startbin = " << startbin << ", lastbin = " << lastbin);
  // build integral distribution
  Double_t sum = 0;
  for (Int_t bin = startbin; bin <= lastbin; bin++) {
    sum += histo->GetBinContent(bin);
    htemp->SetBinContent(bin, (Double_t) sum);
  }
  htemp->Draw();
  DEBUG("smoothing");
  // now loop over integral distribution and smooth it
  TH1D * htemp2 = new TH1D(*htemp);
  htemp2->Reset();
  for (Int_t bin = startbin; bin <= lastbin; bin++) {
    Int_t leftbin = bin - nbins;
    Int_t rightbin = bin + nbins;
//      if (leftbin < startbin) {
//        leftbin  = startbin;
//        rightbin = bin + startbin-leftbin;
//      }
//      if (rightbin > lastbin) {
//        leftbin += rightbin-lastbin;
//        rightbin = lastbin;
//      }
    Double_t newvalue = htemp->Integral(leftbin, rightbin) / (rightbin-leftbin+1);
    INFO("leftbin = " << leftbin << ", rightbin = " << rightbin << ", value = " << newvalue);
    htemp2->SetBinContent(bin, newvalue);
  }
  htemp2->SetLineColor(kRed);
  htemp2->Draw("same");
  DEBUG("Difference");
  // now copy to original histogram
  for (Int_t bin = lastbin-nbins; bin >= startbin+nbins; bin--) {
    histo->SetBinContent(bin, htemp2->GetBinContent(bin)-htemp2->GetBinContent(bin-1));
  }
  DEBUG("delete");
  delete htemp2;
  delete htemp;
}

TH1D * join(Int_t nhistos, TH1D * histo[])
{
  // join the given histograms into one long histogram
  static Int_t joincount = 0;

  joincount++;
  // first loop: count bins
  Int_t maxbins = 0;
  for (Int_t n = 0; n < nhistos; n++) {
    maxbins += histo[n]->GetNbinsX();
  }
  // create new histogram
  TH1D * hnew = new TH1D(Form("hjoin%d", joincount), Form("Composition of %d histograms", nhistos),
			 maxbins, 0, maxbins);
  // fill new histogram
  Int_t binindex = 1;
  for (Int_t n = 0; n < nhistos; n++) {
    for (Int_t j = 0; j < histo[n]->GetNbinsX(); j++) {
      hnew->SetBinContent(binindex+j, histo[n]->GetBinContent(j+1));
    }
    binindex += histo[n]->GetNbinsX();
  }
  return hnew;
}

// execute the same command (given as string) on all plots in current canvas
void all(const char * command)
{
  for (int pad = 0; pad < nPad; pad++) {
    cd(pad+1);
    gROOT->ProcessLine(command);
  }
}

void grid(Int_t xy)
{
  gGrid = xy;
}

void nogrid()
{
  grid(0);
}
