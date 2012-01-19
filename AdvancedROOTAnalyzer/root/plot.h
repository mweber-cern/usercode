#ifndef _plot_h__
#define _plot_h__

#include <iostream>

#include "TEnv.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraph.h"

//////////////////////////////////////////////////////////////////////
// preprocessor macros for logging
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

//////////////////////////////////////////////////////////////////////
// global configuration options etc.

// global Variable to log errors (1), warnings (2), info (3), debug(4,5,...)
extern int            gLogLevel;

extern TEnv         * gConfig;

// maximum number of processes (signal+background+data)
extern Int_t          gMaxProcess;

// drawing canvas
extern TCanvas      * gCanvas;

// language for plots
extern Bool_t         gGerman;

// automatically arrange histograms in plot
extern Bool_t         gAutoOrder;

// number of data taking periods
extern Int_t          gMaxPeriod;

// number of analysis stage
extern Int_t          gStage;

// max. number of subpads in canvas
const  Int_t          gMaxPad = 9;
extern Int_t          gPadNr;

// the histograms for each process
extern TH1D        ** gHisto[gMaxPad];   // size gMaxProcess
extern TH1D        ** gShadow[gMaxPad];  // size gMaxProcess
extern TH2D        ** gHisto2;         // size gMaxProcess

// process drawing options and more
struct TProcess {
  char       * fname;     // file name must be the same for each period
  char       * tname;     // title name
  Int_t        fcolor;    // color for fill drawing
  Int_t        lcolor;    // line color
  Style_t      lstyle;    // line style
  Int_t        hcolor;    // hatch color
  Style_t      hstyle;    // hatch style
  Style_t      marker;    // marker style
  Int_t        mcolor;    // marker color
  Size_t       msize;     // marker size
  Bool_t       join;      // join with MC plotted above?
};

// the process names
extern TProcess     * gProcess; // size gMaxProcess

// the arrows 
extern TArrow * ArrowMin[gMaxPad];
extern TArrow * ArrowMax[gMaxPad];

// order in which the MC's are displayed / added
extern Int_t        * gOrder[gMaxPad]; // size gMaxProcess

// LHC periods
extern char  ** gPeriod; // size gMaxPeriod

// option files containing cut values
extern TEnv        ** gCuts; // size gMaxPeriod

// default option value returned if no option found
const  Double_t       gOptDefault = -1234.56;


// start and end in period array
extern Int_t          gStart;
extern Int_t          gEnd;

struct TProcessInfo {
  Double_t     xs; // cross-section
  Long_t       Nev; // Number of events
  Double_t     weight; // weight (fudge factor)
};

extern TProcessInfo ** gProcessInfo; // size [gMaxPeriod][gMaxProcess]

// The luminosity for each period (in 1/pb)
extern Double_t     * gLumi; // size gMaxPeriod

// how many signal MCs exist
extern Int_t          gMaxSignal;

// which signal mc's to regard as signal and which as background
// bit field for signal:
extern Int_t          gSignal;

// subdir to base directory (set by selection())
extern const char   * gSubDir;

// base directory
extern const char   * gBase;

// draw histograms in color or hatch histograms
extern Bool_t         gIsColor;

// variables for poissonian binned likelihood fit
extern TH1D         * gFitSignal;
extern TH1D         * gFitBackground;
extern TH1D         * gFitData;


//////////////////////////////////////////////////////////////////////
// functions used by other macros
void MakeCanvas(Int_t dx = 1, Int_t dy = 2);
void findbins(Double_t xlow, Double_t xup, Int_t & startbin, Int_t & endbin);
Int_t FindFirstHisto();
void compose();
void decompose();
TH1D * addperiod(Int_t process, const char * hname,
		 const char * selection, Int_t nbins, Double_t min,
		 Double_t max);
void cd(Int_t canvas);
void setopt(TStyle * style);
void setopt(TCanvas * canvas);
void setopt(TH1 * histo);
void setopt(TLegend * leg);
void setopt(TGraph * gr);

#endif