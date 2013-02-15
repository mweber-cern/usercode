#include "TMath.h"
#include "TH1D.h"

#include "plot.h"
#include "rpv.h"
#include "fakerate.h"

#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

struct proc {
  double N;          // number of events
  double staterr;    // MC statistical error
  double systfactor; // systematical uncertainty factor (percent/100) e.g. cross-section
  const char * name; // process name
};

const int nMax = 19;
const int nbins = 6;
const int firststage = 7; // first stage that contains the first bin

proc procs[nbins][nMax];

proc sumprocs[nbins][nMax];

int readprocs(const char * fname)
{
  FILE * infile = fopen(fname, "r");
  if (infile == 0) {
    ERROR("could not open file " << fname);
    return 0;
  }
  char buffer[256];
  int stage = 0;
  int n = 0;
  int i = 0;
  while (fgets(buffer, 256, infile)) {
    if (!strncmp(buffer, " Stage ", 7)) {
      // found a new mass bin, extract stage
      if (sscanf(buffer, " Stage [%d]", & stage) != 1) {
	ERROR("could not get stage: " << buffer);
      }
      INFO("Stage: " << stage);
      // reset counter for processes in this mass bin
      n = 0;
      i = stage-firststage;
      if (i < 0) {
	ERROR("first stage in file must be at least " << firststage);
	return 0;
      }
      else if (i >= nbins) {
	ERROR("last stage must be less than " << firststage + nbins);
	return 0;
      }
    }
    else if (!strncmp(buffer, " +", 2) || !strncmp(buffer, " *", 2)) {
      buffer[1] = '+'; // fix for sscanf below, expects +
      // found a MC or data
      if (stage == 0) {
	ERROR("stage is null!");
	return 0;
      }
      char name[256];
      if (sscanf(buffer, " + %lf +/- %lf %s", 
		 & procs[i][n].N, 
		 & procs[i][n].staterr,
		 name) != 3) {
	ERROR("could not scan line: " << buffer);
	return 0;
      }
      procs[i][n].name = strdup_new(name);
      // assign cross-section error
      if (!strcmp(name, "WW")) {
	procs[i][n].systfactor = 1.5/43.;
      }
      else if (!strcmp(name, "WZ_Q")) {
	procs[i][n].systfactor = 0.7/18.2;
      }
      else if (!strcmp(name, "WZ_Nu")) {
	procs[i][n].systfactor = 0.7/18.2;
      }
      else if (!strcmp(name, "ZZ_nu")) {
	procs[i][n].systfactor = 0.15/5.9;
      }
      else if (!strcmp(name, "ZZ_Q")) {
	procs[i][n].systfactor = 0.15/5.9;
      }
      else if (!strcmp(name, "ZZ_L")) {
	procs[i][n].systfactor = 0.15/5.9;
      }
      else if (!strcmp(name, "fakes")) {
	procs[i][n].systfactor = 0.3;
      }
      else if (!strcmp(name, "data")) {
	procs[i][n].systfactor = 0;
	procs[i][n].staterr = 0; // fix statistical error
      }
      else {
	procs[i][n].systfactor = 0.5;
      }
      INFO(n << " -> " << procs[i][n].name << ": " << procs[i][n].N << " +/- " 
	   << procs[i][n].staterr << " (stat), syst: " << procs[i][n].systfactor);
      n++;
      if (n > nMax) {
	ERROR("Maximum number of processes reached, increase nMax");
      }
    }
    else {
      WARNING("unknown line: " << buffer);
    }
  }
  fclose(infile);
  return n;
}


double fakes[nbins][2];

void fill_fakes(const char * sel = "default29", const char * hname = "btag_jjmm_m")
{
  // six bins with value and error in each bin
  get2dstatistics(fake_estimate_2d(sel, hname), fakes);
}

void summarize_procs(const char * sel = "default29", const char * hname = "btag_jjmm_m")
{
  fill_fakes(sel, hname);
  int maxProc = readprocs("DoublePrompt_MC.txt");
  ofstream out("sum_MC.txt");
  for (int i = 0; i < nbins; i++) {
    out << " Stage [" << i+firststage << "]" << endl;
    double sumN = 0;
    double sumErr2 = 0;
    const char * name = 0;
    for (int j = 0; j < maxProc; j++) {
      name = 0;
      // give a summary at these names
      if (!strcmp(procs[i][j].name, "TTWplus")) {
	name = "VVV";
      }
      else if (!strcmp(procs[i][j].name, "DoublePartonWW")) {
	name = "tt+V";
      }
      else if (!strcmp(procs[i][j].name, "WW")) {
	name = "rare";
      }
      else if (!strcmp(procs[i][j].name, "data")) {
	name = "VV";
      }
      else
	name = 0;
      if (name != 0) {
	out << " + " << sumN << " +/- " << TMath::Sqrt(sumErr2) << " " << name << endl;
	sumN = 0;
	sumErr2 = 0;
      }
      sumN += procs[i][j].N;
      sumErr2 += TMath::Power(procs[i][j].staterr, 2);
    }
    // output fake line
    out << " + " << fakes[i][0] << " +/- " << fakes[i][1] << " fakes" << endl;
    // output data line
    out << " * " << sumN << " +/- " << TMath::Sqrt(sumErr2) << " data" << endl;
    out << "------------------------------------------------------" << endl;
  }
  out.close();
}

void background_systematics()
{
  const int precision = 3;
  cout.precision(precision);
  int maxProc = readprocs("sum_MC.txt");
  // global systematic error, Table 5 of AN
  const double systerrglobal = TMath::Sqrt(1.*1.+2.*2+1.*1.+1.*1.+1.*1.+6.*6.)/100.;
  cout << "global systematic error in %: " << systerrglobal << endl;
  
  double globN = 0;
  double globErrorUncorrelated2 = 0;
  double globErrorCorrelated = 0;
  // walk through table bins and columns and summarize
  for (int bin = 0; bin < nbins; bin++) {
    INFO("bin: " << bin);
    double binN = 0;
    double binErrorUncorrelated2 = 0; 
    double binErrorCorrelated2 = 0; // correlated with other bins
    for (int n = 0; n < maxProc; n++) {
      double locN = 0;
      double locErrorCorrelated2 = 0;
      double locErrorUncorrelated2 = 0;
      if (strcmp(procs[bin][n].name, "data")) {
	locN = procs[bin][n].N;
	// MC cross-section
	locErrorCorrelated2 += TMath::Power(procs[bin][n].N*procs[bin][n].systfactor, 2.);
	// MC statistics
	locErrorUncorrelated2 += TMath::Power(procs[bin][n].staterr, 2.);
	// global systematic error, see above
	locErrorCorrelated2 += TMath::Power(procs[bin][n].N*systerrglobal, 2.);
      }
      INFO(setw(7) << procs[bin][n].name 
	   << ": " << setw(precision+3) << procs[bin][n].N
	   << " +/- " << setw(precision+3) << TMath::Sqrt(locErrorCorrelated2+locErrorUncorrelated2));

      binN += locN;
      binErrorCorrelated2 += locErrorCorrelated2;
      binErrorUncorrelated2 += locErrorUncorrelated2;
    }
    INFO(setw(7) << "sum" 
	 << ": " << setw(precision+3) << binN 
	 << " +/- " << setw(precision+3) 
	 << TMath::Sqrt(binErrorUncorrelated2+binErrorCorrelated2)
	 << " [ +/- " << setw(precision+3) 
	 << TMath::Sqrt(binErrorUncorrelated2) << " (uncorr) "
	 << " +/- " << setw(precision+3) 
	 << TMath::Sqrt(binErrorCorrelated2) << " (corr) ]");
    globN += binN;
    globErrorUncorrelated2 += binErrorUncorrelated2;
    globErrorCorrelated += TMath::Sqrt(binErrorUncorrelated2);
  }
  INFO(setw(7) << "sum" 
       << ": " << setw(precision+3) << globN 
       << " +/- " << setw(precision+3) 
       << TMath::Sqrt(globErrorUncorrelated2+globErrorCorrelated*globErrorCorrelated)
       << " [ +/- " << setw(precision+3) 
       << TMath::Sqrt(globErrorUncorrelated2) << " (uncorr) "
       << " +/- " << setw(precision+3)
       << globErrorCorrelated << " (corr) ]");

  // compute sum over rows with correlated systematics
  globN = 0;
  globErrorUncorrelated2 = 0;
  globErrorCorrelated = 0;
  for (int n = 0; n < maxProc; n++) {
    double procN = 0;
    double procErrorUncorrelated2 = 0; 
    double procErrorCorrelated = 0; // correlated with other bins
    for (int bin = 0; bin < nbins; bin++) {
      double locN = 0;
      double locErrorCorrelated2 = 0;
      double locErrorUncorrelated2 = 0;
      locN += procs[bin][n].N;						       
      if (strcmp(procs[bin][n].name, "data")) {
	locErrorCorrelated2 += TMath::Power(procs[bin][n].N*procs[bin][n].systfactor, 2);  
	locErrorUncorrelated2 += TMath::Power(procs[bin][n].staterr, 2);                 
	locErrorCorrelated2 += TMath::Power(procs[bin][n].N*systerrglobal, 2.);
      }
      procN += locN;
      procErrorCorrelated += TMath::Sqrt(locErrorCorrelated2);
      procErrorUncorrelated2 += locErrorUncorrelated2;
    }
    INFO(setw(7) << procs[0][n].name 
	 << ": " << setw(precision+3) << procN 
	 << " +/- " << setw(precision+3)
	 << TMath::Sqrt(procErrorUncorrelated2+procErrorCorrelated*procErrorCorrelated)
	 << " [ +/- " << setw(precision+3)
	 << TMath::Sqrt(procErrorUncorrelated2) << " (uncorr) "
	 << " +/- " << setw(precision+3)
	 << procErrorCorrelated << " (corr) ]");
    if (strcmp(procs[0][n].name, "data")) {
      globN += procN;
      globErrorCorrelated = TMath::Sqrt(procErrorCorrelated*procErrorCorrelated+procErrorUncorrelated2);
      globErrorUncorrelated2 += procErrorUncorrelated2;
    }
  }
  INFO(setw(7) << "sum" 
       << ": " << setw(precision+3) << globN 
       << " +/- " << setw(precision+3) 
       << TMath::Sqrt(globErrorUncorrelated2+globErrorCorrelated*globErrorCorrelated)
       << " [ +/- " << setw(precision+3) 
       << TMath::Sqrt(globErrorUncorrelated2) << " (uncorr) "
       << " +/- " << setw(precision+3)
       << globErrorCorrelated << " (corr) ]");
}

// /** @todo Need to make this function display MC statistics independent of
//  * other syst. uncertainties. */
// void deltaN()
// {
//   // global systematic error, Table 5 of AN
//   const double systerrglobal = TMath::Sqrt(4.*4.+3.*3.+6.*6.+1.*1.+1.*1.)/100.;
//   cout << "global systematic error in %: " << systerrglobal << endl;
  
  // double totn = 0;
  // double totstaterr2 = 0;
  // double totsysterr2 = 0;
  // double corrsyserr = 0;
  // double corrsyserrtot = 0;
//   for (int n = 0; n < nMax; n++) {
//     double locstaterr2 = 0.;
//     double locsysterr2 = 0.;
//     if (strstr(procs[n].name, "data"))
//       continue;
//     if (!strcmp(procs[n].name, "sf") || !strcmp(procs[n].name, "sf")) {
//       // correlated sys. error from data driven method.  
//       // The error is correlated between channels WJets and QCD and correlated
//       // with other backgrounds due to subtraction procedure
//       corrsyserr = procs[n].N*procs[n].systfactor;
//       corrsyserrtot += corrsyserr;
//     }
//     else {
//       // cross-section uncertainty 
//       corrsyserr = 0;
//       locsysterr2 += TMath::Power(procs[n].N*procs[n].systfactor, 2.);
//       // MC statistics
//       locstaterr2 += TMath::Power(procs[n].staterr, 2.);
//       // global systematic error, see above
//       locsysterr2 += TMath::Power(procs[n].N*systerrglobal, 2.);
//     }
//     totn    += procs[n].N;
//     totsysterr2 += locsysterr2;
//     totstaterr2 += locstaterr2;
//     cout << procs[n].name << ": " << procs[n].N 
// 	 << " +/- " << TMath::Sqrt(locstaterr2) << " (stat)" 
// 	 << " +/- " << TMath::Sqrt(locsysterr2)+corrsyserr << " (syst)" << endl;
//   }
//   cout << "Total number of background events: " << totn 
//        << " +/- " << TMath::Sqrt(totstaterr2) << " (stat)"
//        << " +/- " << TMath::Sqrt(totsysterr2+corrsyserrtot*corrsyserrtot)  << " (syst)" << endl;
//   cout << "Total syst corr err " << corrsyserrtot << endl;
//   cout << "Total syst uncorr err " << TMath::Sqrt(totsysterr2) << endl;
// }

// systematics()
// {
//   deltaN();
// }
