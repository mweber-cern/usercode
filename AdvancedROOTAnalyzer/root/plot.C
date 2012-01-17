#include <sstream>

using namespace std;

#include "TLatex.h"
#include "TColor.h"
#include "TMath.h"
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"

#include "plot.h"

//////////////////////////////////////////////////////////////////////
// initialize variables from plot.h
int             gLogLevel = 3;
TEnv          * gConfig = 0;
Int_t           gMaxProcess = 0;
TCanvas       * gCanvas = 0;
Bool_t          gGerman = kFALSE;
Bool_t          gAutoOrder = kTRUE;
Int_t           gMaxPeriod = 0;
Int_t           gStage = 0;
Int_t           gPadNr = 0;
TH1D         ** gHisto[gMaxPad] = { 0 };   // size gMaxProcess
TH1D         ** gShadow[gMaxPad] = { 0 };  // size gMaxProcess
TH2D         ** gHisto2 = 0;           // size gMaxProcess
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
Int_t           gSignal = 1;
const char    * gSubDir = ".";
const char    * gBase = 0;
Bool_t          gIsColor = kTRUE;
TH1D          * gFitSignal = 0;
TH1D          * gFitBackground = 0;
TH1D          * gFitData = 0;



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
  // style options
  style->SetStatColor(10);        // white background for statistics box
  style->SetPadColor(19);         // pads with grey background
  style->SetOptTitle(0);          // don't show histogram title
//   style->SetTitleColor(10);    // title background white
  style->SetFillColor(10);        // fill everything with a white background
  style->SetFillStyle(1001);      // solid fill style
  style->SetCanvasBorderMode(0);  // no borders in canvases please!
  style->SetCanvasBorderSize(0);  // no borders in canvases please!
  style->SetCanvasColor(10);      // white canvases please
  style->SetPadBorderMode(0);     // no borders in Pads please!
  style->SetPadLeftMargin(0.15);
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
  style->SetLabelSize(0.03, "Y"); // title offset
  style->SetLabelSize(0.03, "X"); // title offset
  style->SetTitleXSize(0.03);     // title size
  style->SetTitleYSize(0.03);     // title size
}

void setopt(TCanvas * canvas)
{
  // set default options for canvas
  canvas->SetTopMargin(0.03);
  canvas->SetRightMargin(0.05);
  canvas->SetBottomMargin(0.15);
  canvas->SetLeftMargin(0.15);
  canvas->SetFillColor(999);
  canvas->SetFillStyle(1001);
  canvas->SetBorderSize(0);
  canvas->SetBorderMode(0);
}

void setopt(TH1 * histo)
{
  // set histo default options
  histo->SetLabelSize(0.06, "X"); // title offset
  histo->SetLabelSize(0.06, "Y"); // title offset
  histo->SetTitleSize(0.06, "X");  // title size
  histo->SetTitleSize(0.06, "Y");  // title size
  histo->SetTitleOffset(1.3, "Y"); // title offset
  histo->SetTitleOffset(1.1, "X"); // title offset
//    histo->GetXaxis()->SetTitleColor(kBlack);
//    histo->GetYaxis()->SetTitleColor(kBlack);
//    histo->GetXaxis()->CenterTitle();
}

// void setopt(TH2D * histo)
// {
//   // set histo default options
//   histo->SetLabelSize(0.06, "X"); // title offset
//   histo->SetLabelSize(0.06, "Y"); // title offset
//   histo->SetTitleSize(0.06, "X");  // title size
//   histo->SetTitleSize(0.06, "Y");  // title size
//   histo->SetTitleOffset(1.1, "Y"); // title offset
//   histo->SetTitleOffset(1.1, "X"); // title offset
// }

void setopt(TLegend * leg)
{
  // set legend default option
  leg->SetBorderSize(0);
  leg->SetFillColor(999);
}

void setopt(TGraph * gr)
{
  TH1F * histo = gr->GetHistogram();
  if (histo != 0)
    setopt(histo);
  gr->SetMarkerStyle(8);
  gr->SetMarkerSize(1.1);
}


// read configuration from file
void read_config_file(const char * configFileName = "Overview.cfg")
{
  if (gConfig) {
    delete gConfig;
    gConfig = 0;
  }
  INFO("Opening config file " << configFileName);
  gConfig = new TEnv(configFileName);

  // number of data taking periods
  gMaxPeriod = gConfig->GetValue("N_period", -1);
  if (gMaxPeriod < 1) {
    ERROR("Number of periods N_period not found in config file or N_period < 1");
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
  if (N_bg < 1) {
    ERROR("Number of backgrounds N_bg not found in config file or N_bg < 1");
    return;
  }
  INFO("Number of backgrounds: " << N_bg);

  // get number of signal MCs
  gMaxSignal = gConfig->GetValue("N_sg", -1);
  if (gMaxSignal < 1) {
    ERROR("Number of signals N_sg not found in config file or N_sg < 1");
    return;
  }
  INFO("Number of signals: " << gMaxSignal);

  // number of processes
  gMaxProcess = N_bg + gMaxSignal + 1;

  // create and initialize arrays depending on gMaxProcess
  for (int pad = 0; pad < gMaxPad; pad++) {
    gHisto[pad]  = new TH1D * [gMaxProcess];
    gShadow[pad] = new TH1D * [gMaxProcess];
    gOrder[pad]  = new Int_t[gMaxProcess];
    for (int i = 0; i < gMaxProcess; i++) {
      gHisto[pad][i]  = 0;
      gShadow[pad][i] = 0;
      gOrder[pad][i] = i;
    }
  }
  gHisto2      = new TH2D * [gMaxProcess];
  gProcess     = new TProcess[gMaxProcess];
  for (int i = 0; i < gMaxProcess; i++) {
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
    gProcess[i].lcolor = gConfig->GetValue(Form("%s.lcolor", selector),  1);
    gProcess[i].lstyle = gConfig->GetValue(Form("%s.lstyle", selector),  1);
    gProcess[i].hcolor = gConfig->GetValue(Form("%s.hcolor", selector),  999);
    gProcess[i].hstyle = gConfig->GetValue(Form("%s.hstyle", selector),  1001);
    gProcess[i].marker = gConfig->GetValue(Form("%s.marker", selector),  8);
    gProcess[i].mcolor = gConfig->GetValue(Form("%s.mcolor", selector),  1);
    gProcess[i].msize  = gConfig->GetValue(Form("%s.msize", selector),  0.7);
    gProcess[i].join   = gConfig->GetValue(Form("%s.join", selector), kFALSE);
  }
  // data
  gProcess[gMaxProcess-1].fname = strdup_new(gConfig->GetValue("file", "data"));
  gProcess[gMaxProcess-1].tname = strdup_new(gConfig->GetValue("label", "Data"));
  gProcess[gMaxProcess-1].fcolor    = 999;
  gProcess[gMaxProcess-1].lcolor    = 999;
  gProcess[gMaxProcess-1].hcolor    = 999;
  gProcess[gMaxProcess-1].hstyle    = 1001;
  gProcess[gMaxProcess-1].lstyle    = 1;
  gProcess[gMaxProcess-1].marker    = 8;
  gProcess[gMaxProcess-1].mcolor    = 1;
  gProcess[gMaxProcess-1].msize     = .7;
  gProcess[gMaxProcess-1].join      = kFALSE;

  // read which histograms to join (old style)
  bool first = true;
  for (Int_t i = 2; i < gMaxProcess; i++) {
    char * dummy = strdup_new(gConfig->GetValue(Form("joinbg.%d", i), "dummy"));
    if (!strcmp(dummy, "dummy")) {
      delete dummy;
      continue;
    }
    // warn user not to use this type of syntax, it is dangerous
    if (first) {
      WARNING("Old-style join syntax \"joinbg." << i << "\" found in config file");
      first = false;
    }
    Int_t id;
    istringstream ss(dummy);
    for (Int_t j = 0; j < i; j++) {
      ss >> id;
      if (id < 0 || id > gMaxProcess || ss.bad()) {
	ERROR("reading parameter #" << j << " of joinbg." << i << " (value " << id);
	delete dummy;
	return;
      }
      // do not join the first one
      if (j > 0)
	gProcess[id].join = kTRUE;
    }
    delete dummy;
  }

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
// setup routine

// this function should be called e.g. by your rootlogon.C, before any functions of the
// plot.C code are called.
void setup(const char * configFileName = "Overview.cfg")
{
  // do initialization only once
  static Bool_t initialized = kFALSE;
  if (initialized)
    return;

  INFO("Setting up plot macro");

  // set plot style
  setopt(gStyle);

  // define a real white color
  TColor mWhite(999, 1., 1., 1., "mWhite");

  // read configuration file
  read_config_file(configFileName);

  // read cross-section and lumi
  read_xs_and_lumi();

  // initialization done
  initialized = kTRUE;
}


//////////////////////////////////////////////////////////////////////
// routines for drawing histograms

void cd(Int_t canvas)
{
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
  char name[14];
  for (Int_t i = 0; i < dy; i++) {
    sprintf(name, "gCanvas_%d", i+1);
    TPad * pad = new TPad(name, name, 0.001, 1-(i+1)*ysize, percent-0.01, 1-i*ysize,
			  999);
    pad->SetNumber(i+1);
    pad->Draw();
    sprintf(name, "gCanvas_%d", 10*(i+1));
    pad = new TPad(name, name, percent, 1-(i+1)*ysize, 0.99, 1-i*ysize, 999);
    pad->SetNumber(10*(i+1));
    pad->Draw();
  }
  gCanvas->Modified();
  // activate first Pad
  cd(gPadNr+1);
}

Int_t FindFirstHisto()
{
  // find a histogram that has been filled and return its index in gHisto[gPadNr]
  // if no histogram existent, return -1
  for (Int_t i = 0; i < gMaxProcess; i++) {
    if (gHisto[gPadNr][gOrder[gPadNr][i]] != 0) {
      return gOrder[gPadNr][i];
    }
  }
  ERROR("No MC found");
  return -1;
}

void color(Bool_t on = kTRUE)
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
  sprintf(fpath, "%s/%s/%s", gBase, gSubDir, gPeriod[period]);
  return fpath;
}

void readcuts()
{
  // read in cut values from option file
  for (Int_t period = gStart; period <= gEnd; period++) {
    // create new TOption object
    if (gCuts[period] != 0)
      delete gCuts[period];
    INFO("Reading option file " << Form("%s/cuts.cfg", getpath(period)));
    gCuts[period] = new TEnv(Form("%s/cuts.cfg", getpath(period)));
  }
}

void selection(const char * subdir)
{
  // set up if not yet done
  setup();

  // save in global variable
  if (gSubDir && strcmp(gSubDir, "."))
    delete gSubDir;
  gSubDir = strdup_new(subdir);

  if (strstr(gSubDir, "final")) {
    // plot signal up for final selection
    top();
  }
  else {
    // plot signal down for preselection
    bottom();
  }

  // read in option file
  readcuts();
}

void periods(char * startperiod, char * endperiod = 0)
{
  // select period range to be displayed
  setup();

  if (endperiod == 0) {
    endperiod = startperiod;
  }

  for (gStart = 0;
       (strcmp(gPeriod[gStart], startperiod) && strcmp(gPeriod[gStart], endperiod));
       gStart++)
    ;
  if (gStart == gMaxPeriod) {
    ERROR("Period " << startperiod << " is not a valid period");
    return;
  }

  for (gEnd = 0; (gPeriod[gEnd] != endperiod) && (gEnd < gMaxPeriod); gEnd++)
    ;
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

const char * GetUnit(const char * hname)
{
  // get unit of histogram name for labelling y-axis automatically
  static char * unit = 0;

  if (unit != 0) {
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

  DEBUG("GetXTitle() start\n");

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

  char * upos = strchr(histo->GetTitle(), '@');
  if (upos != 0) {
    const char * unit = GetUnit(histo->GetTitle());
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
  // create new string containing only y title of histogram name
  static char * title = 0;  // the y-title
  DEBUG("GetYTitle() start");
  if (title) {
    delete[] title;
    title = 0;
  }
  // find position where to separate the two strings
  char * pos = strchr(histo->GetTitle(), ':');

  // get unit
  const char * unit = GetUnit(histo->GetTitle());

  if (pos == 0) {
    // no y-title given
    Double_t binsize = histo->GetBinWidth(1);
    // unit given
    if (unit) {
      if (gGerman) {
	title = new char[26+strlen(unit)];
	sprintf(title, "Anzahl Ereignisse / %4.2f %s", binsize, unit);
      }
      else {
	title = new char[25+strlen(unit)];
	sprintf(title, "Number of events / %4.2f %s", binsize, unit);
      }
    }
    else {
      if (gGerman) {
	title = new char[26];
	sprintf(title, "Anzahl Ereignisse / %4.2f", binsize);
      }
      else {
	title = new char[25];
	sprintf(title, "Number of events / %4.2f", binsize);
      }
    }
  }
  else {
    Int_t len = strlen(histo->GetTitle())-strlen(pos);
    title = new char[len+1];
    strncpy(title, histo->GetTitle(), len);
    title[len] = '\0';
  }

  DEBUG("GetYTitle() end");
  return title;
}

Double_t GetOpt(Int_t ene, const char * hname, const char * MinMax)
{
  // get option value with suffix MinMax from histogram hname
  if (!gCuts[ene]) {
    ERROR("option file not present, return default " << gOptDefault);
  }
  char optname[256];
  // get value
  sprintf(optname, "%s%s", hname+2, MinMax);
  Double_t value = gCuts[gStart]->GetValue(optname, gOptDefault);
  if (value == gOptDefault) {
    ERROR("Option " << optname << " not found, return default " << gOptDefault);
  }
  return value;
}

Double_t GetOptMin(Int_t ene, const char * hname)
{
  // get minimum value from option file for energy gPeriod[gStart]
  return GetOpt(ene, hname, "Min");
}

Double_t GetOptMax(Int_t ene, const char * hname)
{
  // get maximum value from option file for energy gPeriod[gStart]
  return GetOpt(ene, hname, "Max");
}

TArrow * arrow(Double_t position, Int_t neighbourbins = 0)
{
  DEBUG("start arrow");
  // print arrow to indicate a cut
  Int_t hstart = FindFirstHisto();
  if (hstart < 0)
    return 0;
  Double_t fmax = gHisto[gPadNr][hstart]->GetMaximum();
  TAxis * axis = gHisto[gPadNr][hstart]->GetXaxis();
  if (axis == 0) {
    ERROR("No axis found");
    return 0;
  }
  Int_t fbin = axis->FindFixBin(position);
  DEBUG("fbin = " << fbin << ", max = " << gHisto[gPadNr][hstart]->GetNbinsX());
  Double_t fbottom = 1E-8;
  for (Int_t process = 0; process < gMaxProcess; process++) {
    if (gHisto[gPadNr][process]) {
      for (Int_t i = fbin - neighbourbins; i <= fbin + neighbourbins; i++) {
        Double_t temp = gHisto[gPadNr][process]->GetBinContent(i);
        if (temp > fbottom)
          fbottom = temp;
      }
    }
  }
  Double_t ftop;
  if (gPad->GetLogy() == 0) {
    fbottom += fmax * 0.1;
    ftop     = TMath::Min(fbottom + fmax * 0.25, 0.95*fmax);
  } else {
    fbottom  = TMath::Log10(fbottom) + 0.1 * TMath::Log10(fmax);
    ftop     = TMath::Min(fbottom + 0.4 * TMath::Log10(fmax), 0.95*TMath::Log10(fmax));
  }
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
  const char * t = gHisto[gPadNr][hstart]->GetName();
  // only draw arrows for N-1 plots
  if (t[1] != 'n')
    return;
  // get values from option file
  Double_t xmin   = GetOptMin(gStart, t);
  Double_t xmax   = GetOptMax(gStart, t);
  if (xmin != gOptDefault)
    ArrowMin[gPadNr] = arrow(xmin);
  if (xmax != gOptDefault)
    ArrowMax[gPadNr] = arrow(xmax);
  DEBUG("end drawcut\n");
}

void draw(Bool_t autotitle = kTRUE)
{
  DEBUG("start draw");

  // draw the histograms in the canvas
  Bool_t first = kTRUE;   // needed for drawing first histogram without "same"

  // clear drawing pad and delete histograms
  gPad->cd();
  gPad->Clear();

  // now draw all histograms in canvas
  TH1D * histo;
  for (Int_t i = 0; i < gMaxProcess; i++) {
    Int_t process = gOrder[gPadNr][i];
    LOG(10, "process " << process << "[" << gProcess[process].fname << "]");
    histo = gHisto[gPadNr][process];
    if (histo == 0)
      continue;

    DEBUG("draw 1");
    histo->SetTickLength(-0.02);
    histo->SetLineStyle(gProcess[process].lstyle);

    if (process == gMaxProcess - 1) {
      // data
      histo->SetMarkerStyle(gProcess[process].marker);
      histo->SetMarkerColor(gProcess[process].mcolor);
      histo->SetMarkerSize(gProcess[process].msize);
    }
    else {
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
        shadow->SetFillColor(999);
        shadow->Draw("same");
        histo->SetFillStyle(gProcess[process].hstyle);
        histo->SetFillColor(gProcess[process].hcolor);
      }
    }
    if (first) {
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
    else if (process == gMaxProcess-1) {
      // draw data
      histo->Draw("e1psame");
    }
    else {
      // draw others, but only if not joined
      if (!gProcess[process].join)
	histo->Draw("histosame");
    }
  }
  // draw cut (if it is there)
  drawcut();
  gPad->Update();
  gPad->RedrawAxis();
  DEBUG("end draw");
}

void max(Double_t maximum)
{
  // set maximum for all histograms
  for (Int_t i = 0; i < gMaxProcess; i++) {
    if (gHisto[gPadNr][i] != 0)
      gHisto[gPadNr][i]->SetMaximum(maximum);
  }
  draw();
}

void min(Double_t minimum)
{
  // set minimum for all histograms
  for (Int_t i = 0; i < gMaxProcess; i++) {
    if (gHisto[gPadNr][i] != 0)
      gHisto[gPadNr][i]->SetMinimum(minimum);
  }
  draw();
}

void findbins(Double_t xlow, Double_t xup, Int_t & startbin, Int_t & endbin)
{
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
    endbin   = gHisto[gPadNr][hstart]->GetNbinsX();
  }
  else {
    startbin = gHisto[gPadNr][hstart]->GetXaxis()->FindFixBin(xlow);
    endbin   = gHisto[gPadNr][hstart]->GetXaxis()->FindFixBin(xup);
  }
}

void findmax(Double_t low = 0, Double_t up = 0)
{
  // find maximum by iterating over all histograms and counting bin content

  // convert x-axis range to bin range
  Int_t start, end;
  findbins(low, up, start, end);

  Double_t lmax = -1E99;
  for (Int_t process = 0; process < gMaxProcess; process++) {
    if (gHisto[gPadNr][process]) {
      for (Int_t i = start; i <= end; i++) {
	Double_t temp = gHisto[gPadNr][process]->GetBinContent(i);
	if (temp > lmax)
	  lmax = temp;
      }
    }
  }
  if (lmax > 0)
    max(lmax*(1+1/TMath::Sqrt(lmax))*1.1);
}

void findmin(Double_t low = 0, Double_t up = 0)
{
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
    if (gHisto[gPadNr][process]) {
      for (Int_t i = start; i <= end; i++) {
	Double_t temp = gHisto[gPadNr][process]->GetBinContent(i);
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
}

void zoom(Double_t low, Double_t up)
{
  // zoom on x-axis
  for (Int_t process = 0; process < gMaxProcess-1; process++) {
    if (gHisto[gPadNr][process]) {
      gHisto[gPadNr][process]->SetAxisRange(low, up);
    }
  }
  draw();
}

void unzoom()
{
  // unzoom x-axis
  for (Int_t process = 0; process < gMaxProcess-1; process++) {
    if (gHisto[gPadNr][process]) {
      TAxis * axis = gHisto[gPadNr][process]->GetXaxis();
      if (axis)
	axis->UnZoom();
    }
  }
  draw();
}

void rebin(Int_t nbins)
{
  // rebin histograms: combine nbins bins into one bin
  if (FindFirstHisto() < 0)
    return;

  // rebin all histograms
  for (Int_t i = 0; i < gMaxProcess; i++) {
    if (gHisto[gPadNr][i])
      gHisto[gPadNr][i]->Rebin(nbins);
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
      TH1D * histo = gHisto[gPadNr][process];
      if (histo == 0)
	continue;
      for (Int_t bin = histo->GetNbinsX(); bin >= 0; bin--) {
	histo->SetBinContent(bin, histo->GetBinContent(bin-nbins));
      }
    }
  }
  else {
    for (Int_t process = 0; process < gMaxProcess-1; process++) {
      TH1D * histo = gHisto[gPadNr][process];
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
    TH1D * histo = gHisto[gPadNr][process];
    TH1D * hnew  = (TH1D *) histo->Clone();
    hnew->Reset();
    Int_t nbinmax = histo->GetNbinsX()+1;
    for (Int_t nbin = 0; nbin <= nbinmax; nbin++) {
      hnew->SetBinContent(nbin, histo->GetBinContent(nbinmax-nbin));
    }
    delete histo;
    gHisto[gPadNr][process] = hnew;
  }
  draw();
}

void print(const char * hname = 0)
{
  // print current canvas in a file "name", substitute current
  // histogram name if no name given
  if (FindFirstHisto() < 0)
    return;

  char fpath[256];
  if (hname == 0) {
    TH1D * h1 = gHisto[gPadNr][0];
    if (h1)
      if (gStart != gEnd)
	sprintf(fpath, "%s-%s-%s.pdf",
		h1->GetName(), gPeriod[gStart], gPeriod[gEnd]);
      else
	sprintf(fpath, "%s-%s.pdf", h1->GetName(), gPeriod[gStart]);
    else
      sprintf(fpath, "%s.pdf", "unknown");
  }
  else
    sprintf(fpath, "%s", hname);

  gCanvas->Print(fpath);
}

void pprint(const char * hname = 0)
{
  // print current pad in a file "name", substitute current
  // histogram name if no name given
  if (FindFirstHisto() < 0)
    return;

  char fpath[256];
  if (hname == 0) {
    TH1D * h1 = gHisto[gPadNr][0];
    if (h1)
      if (gStart != gEnd)
	sprintf(fpath, "%s-%s-%s.pdf",
		h1->GetName(), gPeriod[gStart], gPeriod[gEnd]);
      else
	sprintf(fpath, "%s-%s.pdf", h1->GetName(), gPeriod[gStart]);
    else
      sprintf(fpath, "%s.pdf", "unknown");
  }
  else
    sprintf(fpath, "%s", hname);

  gPad->Print(fpath);
}

void title(const char * title = 0)
{
  // set the axis titles for this plot.
  // If no name given, then the default title is taken
  Int_t first = FindFirstHisto();
  if (first < 0)
    return;

  gHisto[gPadNr][first]->SetTitle(title);

  // now use automatic title
  draw(kTRUE);
}

void logy()
{
  // set y axis to be logarithmic
  gPad->SetLogy(1);
  // adjust minimum
  findmin();
  draw();
}

void liny()
{
  // set y axis to be linear
  if (gCanvas == 0)
    MakeCanvas();
  gPad->SetLogy(0);
  // adjust minimum
  findmin();
  draw();
}

void ascii(const char * fname = 0)
{
  // print #data, #signal, #back in ascii format to screen and into file
  if (FindFirstHisto() < 0)
    return;

  FILE * outfile = 0;
  if (fname) {
    // assemble path name
    char fpath[256];
    sprintf(fpath, "%s/%s/%s.ascii", gBase, gSubDir, fname);
    // convert to lower case
    for (UInt_t i = 0; i < strlen(fpath); i++)
      fpath[i] = tolower(fpath[i]);

    // open file
    outfile = fopen(fpath, "w");
    if (!outfile) {
      ERROR("Could not open file" << fpath);
      return;
    }
    INFO("Output in file " << fpath);
  }

  // loop over histogram bins and save each bin in one separate line
  Double_t nData, nSignal, nBack, nZZBack, temp;
  for (Int_t i = 1; i <= gHisto[gPadNr][0]->GetNbinsX(); i++) {
    nData = gHisto[gPadNr][gMaxProcess-1]->GetBinContent(i);
    // divide other histograms into signal, ZZ background and other background
    nSignal = 0;
    nBack   = 0;
    nZZBack = 0;
    // loop over MC's
    for (Int_t j = 0; j < gMaxProcess-1; j++) {
      // get number of events for this histogram
      temp = gHisto[gPadNr][gOrder[gPadNr][j]]->GetBinContent(i);
      // subtract following mc
      if (j < gMaxProcess-2)
	temp -= gHisto[gPadNr][gOrder[gPadNr][j+1]]->GetBinContent(i);
      if (gOrder[gPadNr][j] < 3)
	// signal
	nSignal += temp;
      else if (strstr(gProcess[gOrder[gPadNr][j]].fname, "background")) {
	// ZZ background
	nZZBack += temp;
      }
      else {
	// other background
	nBack   += temp;
      }
    }
    printf("%10.7f   %10.7f   %10.7f   %10.7f\n", nData, nSignal, nBack, nZZBack);
    if (fname)
      fprintf(outfile, "%f   %f   %f   %f\n", nData, nSignal, nBack, nZZBack);
  }
  if (outfile)
    fclose(outfile);
}

void drawperiod(Int_t posi = -1)
{
  // draw energy in current plot
  static TLatex * t[gMaxPad] = { 0 };
  char etext[50];
  if (gStart != gEnd) {
    sprintf(etext, "CMS %s - %s", gPeriod[gStart], gPeriod[gEnd]);
  }
  else {
    sprintf(etext, "CMS %s", gPeriod[gStart]);
  }

  Int_t hstart = FindFirstHisto();
  if (hstart < 0)
    return;
  if (t[gPadNr] != 0) {
    delete t[gPadNr];
    t[gPadNr] = 0;
  }
  if (posi > 0) {
    t[gPadNr] = new TLatex(0.955, 0.97, etext);
    t[gPadNr]->SetTextAlign(33);
  }
  else if (posi < 0) {
    t[gPadNr] = new TLatex(0.25, 0.97, etext);
    t[gPadNr]->SetTextAlign(13);
  }
  else {
    t[gPadNr] = new TLatex(0.5, 0.97, etext);
    t[gPadNr]->SetTextAlign(23);
  }
  t[gPadNr]->SetNDC();
  t[gPadNr]->Draw();
}

void legend(Double_t mincontent = 0.01, Int_t posi = 1, Double_t miny = -1)
{
  // draw a legend with currently used colors
  static TLegend * t[gMaxPad] = { 0 };
  Int_t position = FindFirstHisto();
  if (position < 0)
    return;

  if (miny < 0) {
    if (strstr(gSubDir, "findsusyb3"))
      miny = 0.4;
    else
      miny = 0.6;
  }

  // count processes
  Int_t num = 0;
  Double_t integral = gHisto[gPadNr][gOrder[gPadNr][position]]->Integral();
  for (Int_t i = 0; i < gMaxProcess; i++) {
    if (gHisto[gPadNr][gOrder[gPadNr][i]] == 0)
      continue;
    if (i < gMaxProcess-2) {
      if (gHisto[gPadNr][gOrder[gPadNr][i+1]]) {
	if ((gHisto[gPadNr][gOrder[gPadNr][i]]->Integral() -
	     gHisto[gPadNr][gOrder[gPadNr][i+1]]->Integral())/integral < mincontent)
	  continue;
      }
    }
    else if (gHisto[gPadNr][gOrder[gPadNr][i]]->Integral()/integral < mincontent)
      continue;

    num++;
  }
  // determine minimum
  miny = TMath::Max(0.92-num*0.06, miny);
  INFO("Found " << num << " histograms, miny = " << miny);

  if (t[gPadNr]) {
    delete t[gPadNr];
    t[gPadNr] = 0;
  }
  if (posi > 0)
    t[gPadNr] = new TLegend(0.720, miny, 0.92, 0.92);
  else if (posi < 0)
    t[gPadNr] = new TLegend(0.173, miny, 0.413, 0.92);
  else
    t[gPadNr] = new TLegend(0.45, miny, 0.69, 0.92);
  setopt(t[gPadNr]);
  for (Int_t i = 0; i < gMaxProcess; i++) {
    if (gHisto[gPadNr][gOrder[gPadNr][i]] == 0)
      continue;
    if (i < gMaxProcess-2) {
      if (gHisto[gPadNr][gOrder[gPadNr][i+1]]) {
	if ((gHisto[gPadNr][gOrder[gPadNr][i]]->Integral() -
	     gHisto[gPadNr][gOrder[gPadNr][i+1]]->Integral())/integral < mincontent)
	  continue;
      }
    }
    else if (gHisto[gPadNr][gOrder[gPadNr][i]]->Integral()/integral < mincontent)
      continue;

    const char * opt = "f"; // fill area
    if (i == gMaxProcess-1) {
      opt = "p"; // data with marker
      if (gGerman)
	t[gPadNr]->AddEntry(gHisto[gPadNr][gOrder[gPadNr][i]], "Daten", opt);
      else
	t[gPadNr]->AddEntry(gHisto[gPadNr][gOrder[gPadNr][i]], gProcess[gOrder[gPadNr][i]].tname, opt);
    }
    else
      t[gPadNr]->AddEntry(gHisto[gPadNr][gOrder[gPadNr][i]], gProcess[gOrder[gPadNr][i]].tname, opt);
  }
  t[gPadNr]->Draw();
  // insert energy information
  drawperiod(posi != 0 ? -posi : -1);
}

void subfigure(const char * subfig = "a)")
{
  // add entry for subfigure in the current canvas, e.g. "a)"
  TText * t = new TText(0.95, 0.97, subfig);
  t->SetTextAlign(33);
  t->SetNDC();
  t->Draw();
}

void legend2(Double_t mincontent = 0.01)
{
  // draw a legend with currently used colors
  static TLegend * t2[gMaxPad] = { 0 };
  Int_t position = FindFirstHisto();
  if (position < 0)
    return;

  // count processes
  Int_t num = 0;
  Double_t integral = gHisto[gPadNr][gOrder[gPadNr][position]]->Integral();
  for (Int_t i = 0; i < gMaxProcess; i++) {
    if (gHisto[gPadNr][gOrder[gPadNr][i]] == 0)
      continue;
    if (i < gMaxProcess-2) {
      if (gHisto[gPadNr][gOrder[gPadNr][i+1]]) {
	if ((gHisto[gPadNr][gOrder[gPadNr][i]]->Integral() -
	     gHisto[gPadNr][gOrder[gPadNr][i+1]]->Integral())/integral < mincontent)
	  continue;
      }
    }
    else if (gHisto[gPadNr][gOrder[gPadNr][i]]->Integral()/integral < mincontent)
      continue;

    num++;
  }

  // determine minimum
  Double_t miny  = TMath::Max(0.92-num*0.1, 0.15);
  INFO("Found " << num << " histograms, miny = " << miny);

  if (t2[gPadNr]) {
    delete t2[gPadNr];
    t2[gPadNr] = 0;
  }

  t2[gPadNr] = new TLegend(0.01, miny, 0.99, 0.99);
  setopt(t2[gPadNr]);
  t2[gPadNr]->SetTextSize(0.12);
  t2[gPadNr]->SetMargin(0.2);
  for (Int_t i = 0; i < gMaxProcess; i++) {
    if (gHisto[gPadNr][gOrder[gPadNr][i]] == 0)
      continue;
    if (i < gMaxProcess-2) {
      if (gHisto[gPadNr][gOrder[gPadNr][i+1]]) {
	if ((gHisto[gPadNr][gOrder[gPadNr][i]]->Integral() -
	     gHisto[gPadNr][gOrder[gPadNr][i+1]]->Integral())/integral < mincontent)
	  continue;
      }
    }
    else if (gHisto[gPadNr][gOrder[gPadNr][i]]->Integral()/integral < mincontent)
      continue;

    const char * opt = "f"; // fill area
    if (i == gMaxProcess-1) {
      opt = "p"; // data with marker
      if (gGerman)
	t2[gPadNr]->AddEntry(gHisto[gPadNr][gOrder[gPadNr][i]], "Daten", opt);
      else
	t2[gPadNr]->AddEntry(gHisto[gPadNr][gOrder[gPadNr][i]],
			     gProcess[gOrder[gPadNr][i]].tname, opt);
    }
    else
      t2[gPadNr]->AddEntry(gHisto[gPadNr][gOrder[gPadNr][i]],
			   gProcess[gOrder[gPadNr][i]].tname, opt);
  }
  DEBUG("now drawing legend");
  Int_t savenr = gPadNr;
  cd(10*(gPadNr+1));
  t2[savenr]->Draw();
  cd(savenr+1);
  // insert energy information
  drawperiod(1);
}

void DeleteOld()
{
  // delete existing histograms gHisto[gPadNr][]
  for (Int_t i = 0; i < gMaxProcess; i++) {
    if (gHisto[gPadNr][i] != 0) {
      delete gHisto[gPadNr][i];
      // mark as deleted
      gHisto[gPadNr][i] = 0;
    }
  }
}

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
    sprintf(fpath, "%s/%s.root", getpath(period), gProcess[process].fname);
    LOG(5, "Opening file " << fpath);
    TFile * f = TFile::Open(fpath);
    if (f == 0 || !f->IsOpen()) {
      if (hadd)
	delete hadd;
	ERROR("Could not open file " << fpath);
      return 0;
    }
    // get key from file
    key = (TKey *) f->GetKey(hname);
    if (key == 0) {
      // key has not been found, generate histogram from tree
      // create specified histogram
      TH1D * hmine = new TH1D("hmine", hname, nbins, min, max);
      // get tree from current file
      const char * treename = gConfig->GetValue("treename", "t");
      TTree * dtree = (TTree *) f->Get(treename);
      if (dtree == 0) {
	// no tree found, try next period
	WARNING("Tree \"" << treename << "\" not found in file " << fpath);
	continue;
      }
      // save tree data in my histogram
      sprintf(mname, "%s>>hmine", hname);
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
      htemp = (TH1D * ) f->Get(hname);
    }
    // check physical minimum
    if (htemp->GetMinimum() < 0) {
      ERROR("Minimum " << htemp->GetMinimum() << " < 0 in file "<< fpath);
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
  for (int i = 0; i < gMaxProcess; i++) {
    INFO("gOrder[" << i << "] = " << gOrder[gPadNr][i]);
  }
}

void auto_order()
{
  double integral[gMaxProcess];
  // order all but data
  for (Int_t i = 0; i < gMaxProcess-1; i++) {
    if (gHisto[gPadNr][i])
      integral[i] = gHisto[gPadNr][i]->Integral();
    else
      integral[i] = 0;
  }
  TMath::Sort(gMaxProcess-1, integral, gOrder[gPadNr]);
}

void compose()
{
  // sum up the histograms in the current order, such that they can be
  // drawn on the screen, one before the other

  // loop over processes
  for (Int_t i = 0; i < gMaxProcess; i++) {
    Int_t process = gOrder[gPadNr][i];
    // don't sum up data (last one)
    if (process < gMaxProcess - 1) {
      // add to all previos ones
      for (Int_t j = 0; j < i; j++) {
	Int_t backmc = gOrder[gPadNr][j];
	if ((gHisto[gPadNr][backmc] != 0) && (gHisto[gPadNr][process] != 0))
	  gHisto[gPadNr][backmc]->Add(gHisto[gPadNr][process]);
      }
    }
  }
}

void decompose()
{
  // subtract summed histograms that are ready for drawing into histograms
  // corresponding to each process...
  if (FindFirstHisto() < 0)
    return;
  // loop over processes, don't subtract from last mc
  for (Int_t i = 0; i < gMaxProcess-2; i++) {
    if (gHisto[gPadNr][gOrder[gPadNr][i]] == 0)
      continue;
    // find next mc
    Int_t j;
    for (j = i+1; (j < gMaxProcess-1) && (gHisto[gPadNr][gOrder[gPadNr][j]] == 0); j++)
      ;
    if (j < gMaxProcess-1)
      gHisto[gPadNr][gOrder[gPadNr][i]]->Add(gHisto[gPadNr][gOrder[gPadNr][j]], -1);
  }
}

void plot(const char * hname, const char * selection = "global_weight",
	  Int_t nbins = 100, Double_t min = 0, Double_t max = 1)
{
  // plot histogram "hname" for the selected selection and period
  // loops over all files.

  // check for canvas
  if (gCanvas == 0) {
    MakeCanvas();
  }

  // delete previously plotted histograms
  DeleteOld();

  // loop over all processes
  for (Int_t process = 0; process < gMaxProcess; process++) {
    // get summary histogram from all periods
    gHisto[gPadNr][process] = addperiod(process, hname,
					selection, nbins, min, max);
  }

  if (gAutoOrder)
    auto_order();
  // compose histograms for drawing
  compose();
  // find maximum & minimum
  findmax();
  findmin();
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
    sprintf(fpath, "%s/%s.root", getpath(period), gProcess[process].fname);
    TFile * f = new TFile(fpath, "READ");
    if (!f->IsOpen()) {
      ERROR("File " << fpath << " not found.");
      return 0;
    }
    // get key from file
    key = (TKey *) f->GetKey(hname);
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
      sprintf(mname, "%s>>htemp", hname);
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
      htemp = (TH2D * ) f->Get(hname);
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

void plot2(const char * hname,
	   const TH2D * hsample = 0, const char * selection = "global_weight",
	   const char * drawopt = "box")
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
      ERROR("Empty histogram found for process " <<  process);
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

TH1D * signalHisto()
{
  // creates a histogram containing only signal, extracted from global
  // histograms gHisto[gPadNr][]. Can be used e.g. for fitting (see xsection.C)
  decompose();
  // loop over background MC's (not data)
  TH1D * hSignal = 0;
  Bool_t first = kTRUE;
  for (Int_t j = 0; j < 3; j++) {
    if (first) {
      hSignal = new TH1D(*gHisto[gPadNr][j]);
      hSignal->SetDirectory(0);
      first = kFALSE;
    }
    else {
      hSignal->Add((TH1D *) gHisto[gPadNr][j], 1);
    }
  }
  compose();

  return hSignal;
}

TH1D * backgroundHisto()
{
  // creates a histogram containing only background, extracted from global
  // histograms gHisto[gPadNr][]. Can be used e.g. for fitting (see xsection.C)
  decompose();
  // loop over background MC's (not data)
  TH1D * hBackground = 0;
  Bool_t first = kTRUE;
  for (Int_t j = 3; j < gMaxProcess-1; j++) {
    if (first) {
      hBackground = new TH1D(*gHisto[gPadNr][j]);
      hBackground->SetDirectory(0);
      first = kFALSE;
    }
    else {
      hBackground->Add((TH1D *) gHisto[gPadNr][j], 1);
    }
  }
  compose();

  return hBackground;
}

TH1D * dataHisto()
{
  // creates a histogram containing only data, extracted from global
  // histograms gHisto[gPadNr][]. Can be used e.g. for fitting (see xsection.C)
  TH1D * hData = new TH1D(*gHisto[gPadNr][gMaxProcess-1]);
  hData->SetDirectory(0);
  return hData;
}

void showintegral()
{
  // plots integral of all mc's
  TCanvas * cint = new TCanvas("cint", "integral distribution",
			       10,10, (Int_t) (600*TMath::Sqrt(2.)), 600);
  cint->cd();
  static TH1D ** hint;
  hint = new TH1D * [gMaxProcess];
  Double_t maxsum = 0;
  // compute bin-to-bin efficiency and purity
  for (Int_t process = 0; process < gMaxProcess; process++) {
    // delete old histos
    if (hint[process]) {
      delete hint[process];
      hint[process] = 0;
    }
    // skip non-existing mc's
    if (gHisto[gPadNr][gOrder[gPadNr][process]] == 0)
      continue;
    // divide mc's into signal & background
    hint[process] = new TH1D(*gHisto[gPadNr][gOrder[gPadNr][process]]);
    hint[process]->Reset();
    Double_t sum = 0;
    // make integral distribution
    for (Int_t bin = 1; bin <= hint[process]->GetNbinsX(); bin++) {
      sum += gHisto[gPadNr][gOrder[gPadNr][process]]->GetBinContent(bin);
      hint[process]->SetBinContent(bin, sum);
    }
    // find absolute maximum
    if (sum > maxsum)
      maxsum = sum;
    // set line styles
    hint[process]->SetLineColor(hint[process]->GetFillColor() == 10
				 ? 1 : hint[process]->GetFillColor());
    hint[process]->SetFillStyle(0);
  }
  // draw histograms
  for (Int_t process = 0; process < gMaxProcess; process++) {
    hint[process]->SetMaximum(maxsum);
    if (process == 0) // first
      hint[process]->Draw();
    else if (process == gMaxProcess-1) { // data
      hint[process]->SetLineStyle(gProcess[process].lstyle);
      hint[process]->SetLineColor(1);
      hint[process]->SetMarkerStyle(8);
      hint[process]->SetMarkerSize(0.6);
      hint[process]->Draw("e1psame");
    }
    else
      hint[process]->Draw("same");
  }
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
  static int joincount = 0;

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
