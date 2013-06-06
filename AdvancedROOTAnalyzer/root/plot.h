#ifndef _plot_h__
#define _plot_h__

#include <iostream>

#include "TEnv.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TClass.h"

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

// version
extern Int_t          gVersion;

// global Variable to log errors (1), warnings (2), info (3), debug(4,5,...)
extern Int_t          gLogLevel;

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
extern TString        gStage;

// max. number of subpads in canvas
const  Int_t          gMaxPad = 9;
// current active pad
extern Int_t          gPadNr;
// current number of pads in canvas
extern Int_t          nPad;

// the histograms for each process
extern TH1D        ** gHisto[gMaxPad];   // size gMaxProcess, unstacked histograms
extern TH1D        ** gStack[gMaxPad];   // size gMaxProcess, stacked histograms
extern TH1D        ** gShadow[gMaxPad];  // size gMaxProcess, useful when using hatch styles
extern TH2D        ** gHisto2;           // size gMaxProcess
extern TH3D        ** gHisto3[gMaxPad];  // size gMaxProcess

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
  Bool_t       stack;     // stack with other histograms
};

// the process names
extern TProcess     * gProcess; // size gMaxProcess

// the arrows 
extern TArrow       * ArrowMin[gMaxPad];
extern TArrow       * ArrowMax[gMaxPad];

// order in which the MC's are displayed / added
extern Int_t        * gOrder[gMaxPad]; // size gMaxProcess

// LHC periods
extern char        ** gPeriod; // size gMaxPeriod

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

// subdir to base directory (set by selection())
extern const char   * gSubDir;

// base directory
extern const char   * gBase;

// draw histograms in color or hatch histograms
extern Bool_t         gIsColor;

// draw a grid in the plots
extern Int_t          gGrid;

// sum up overflows in last bin of histogram
extern Bool_t         gMoveOverflow;

//////////////////////////////////////////////////////////////////////
// functions used by other macros

// everyday use
void selection(const char * subdir);
void period(const char * startperiod, const char * endperiod = 0);
void stage(TString s);
void plot(const char * hname, const char * selection = "global_weight",
	  Int_t nbins = 100, Double_t min = 0, Double_t max = 1);
void plot2(const char * hname,
	   const TH2D * hsample = 0, const char * selection = "global_weight",
	   const char * drawopt = "box");
void plot3(const char * hname,
	   const TH3D * hsample = 0, const char * selection = "global_weight");
void cd(Int_t canvas);
void zoom(Double_t low, Double_t up);
void unzoom();
void max(Double_t maximum);
void min(Double_t minimum);
void liny();
void logy();
void grid(Int_t xy = 11); // turn on grid in plots (10 => X, 1 => Y)
void nogrid(); // turn off grid in plots
void rebin(Int_t nbins);
void legend(Double_t mincontent = 0.01, Int_t posi = 1, Double_t miny = -1);
TArrow * arrow(Double_t position, Int_t neighbourbins = 0);
void print(const char * hname = 0);
void pprint(const char * hname = 0);

// less often used
void version(Int_t version);
void setup(const char * configFileName = "Overview.cfg");
void MakeCanvas(Int_t dx = 1, Int_t dy = 2);
void title(const char * title = 0);
void shiftbin(Int_t nbins);
void mirror();
void smooth(TH1D * histo, Double_t xlow, Double_t xup, Int_t nbins);
void subfigure(const char * subfig = "a)");
void color(Bool_t on = kTRUE);
void top();
void bottom();
void plotadd(const char * name1, const char * name2);
void plotadd(const char * name1, const char * name2,
	     const char * name3, const char * name4);

// get sums of specific histograms for special purposes, e.g. fitting
TH1D * getHisto(const char * name = 0, bool match = true, int min=0, int max = gMaxProcess-1);
TH1D * signalHisto(const char * name = 0, bool match = true);
TH1D * backgroundHisto(const char * name = 0, bool match = true);
TH2D * backgroundHisto2(const char * name = 0, bool match = true);
TH1D * dataHisto();
TH2D * dataHisto2();

// technical stuff used by other macros, do not use directly
char * strdup_new(const char * text);
void setopt(TStyle * style);
void setopt(TCanvas * canvas);
void setopt(TH1 * histo);
void setopt(TLegend * leg);
void setopt(TGraph * gr);
const char * getpath(Int_t period);
void findbins(Double_t xlow, Double_t xup, Int_t & startbin, Int_t & endbin);
TH1D * join(Int_t nhistos, TH1D * histo[]);
Int_t FindFirstHisto();
TH1D * addperiod(Int_t process, const char * hname,
		 const char * selection, Int_t nbins, Double_t min,
		 Double_t max);

class TObject;
TObject * get_object(const char * filename, const char * objectname);

#endif
