#include "plot.h"
#include "xsection.h"

#include <iomanip>

using namespace std;

#include "TH1D.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TText.h"
#include "TRandom.h"
#include "TFile.h"

// variables for poissonian binned likelihood fit
TH1D * gFitSignal = 0;
TH1D * gFitBackground = 0;
TH1D * gFitData = 0;

TH1D * gSign[3][3] = { { 0 }, { 0 }, { 0 } };
TH1D * gBack[3] = { 0 };
TH1D * gData[3] = { 0 };

TH1D * gFitSig[3] = { 0 };

void BinnedPoissonianLikelihood(Int_t    & npar,  // # of internal parameters
				Double_t * grad,  // array of first derivatives
				Double_t & fval,  // the function value
				Double_t * xvar,  // the input variable values
				Int_t      iflag) // status flag
{
  static Int_t      nBins      = 0;
  static Double_t   nTotSig    = 0; // total number of expect. signal events
  static Double_t   nTotBack   = 0; // total number of background events

  // compute negative logarithm of likelihood for binned poissonian data.
  // Accept only one single variable: The relation of the measured 
  // ZZ cross-section to the SM cross-section.

  fval = 1E9; // Always return some value

  if (iflag == 1) {
    if (npar != 1) {
      ERROR("BinnedPoissonianLikelihood() (init) works with just one parameter, "
	    << npar << " given!");
    }

    // Initialization
    DEBUG("BinnedPoissonianLikelihood() Initialization");
    
    // check inputs
    if (!(gFitSignal || gFitBackground || gFitData)) {
      ERROR("BinnedPoissonianLikelihood() no existing histograms");
      return;
    }

    // initialize local variables
    nBins        = gFitData->GetNbinsX();
    nTotSig      = 0;
    nTotBack     = 0;
    // compute total signal and background values
    for (Int_t i = 1; i <= nBins; i++) {
      nTotSig  += gFitSignal->GetBinContent(i);
      nTotBack += gFitBackground->GetBinContent(i);
    }
  } else if (iflag == 2) {
    // compute gradient of function value fval
    grad = 0;
    ERROR("No gradient computation in BinnedPoissonianLikelihood()");
  } else if (iflag == 3) {
    // Things to do after fit has finished
  } else {
    if (npar != 1) {
      WARNING("BinnedPoissonianLikelihood() (default) works with just one parameter, "
	      << npar << " given!");
    }
    Double_t pval; // partial value
    // compute function value fval
    fval = nTotBack + xvar[0] * nTotSig;
    for (Int_t i = 1; i <= nBins; i++) {
      // evade log(zero)
      pval = gFitBackground->GetBinContent(i)
	+ xvar[0] * gFitSignal->GetBinContent(i);
      if (pval > 0)
	fval -= gFitData->GetBinContent(i) * TMath::Log(pval);
      else if (pval < 0)
	WARNING(" log(<0) in BinnedPoissonianLikelihood()\n");
    }
  }
}

void xsratio(Double_t xs[4], // xs, parabolic error, -error, +error
	     Int_t order, Bool_t draw)
{
  // fit the cross-section ratio (measured/theorie) from global
  // histograms. The input parameter determines the confidence level
  // for the errors, counted in one sigma deviations. I.e. 5 means
  // five sigma deviation, 1 means the usual one sigma error (i.e. CL 67 %)
  if (!(gFitSignal || gFitBackground || gFitData)) {
    ERROR("global fit histograms are not initialized");
    if (gFitSignal == 0)
      ERROR("no signal");
    if (gFitBackground == 0)
      ERROR(" no background");
    if (gFitData == 0)
      ERROR("no data");
    return;
  }

  // create new Minuit with a maximum of five parameters
  if (gMinuit == 0) {
    DEBUG("Creating MINUIT with 5 parameters");
    gMinuit = new TMinuit(5);
  }

  // we need this function to be evaluated
  gMinuit->SetFCN(BinnedPoissonianLikelihood);

  Double_t arglist[10]; // argument list
  Int_t ierflg;         // error return code from Minuit
  
  arglist[0] = gLogLevel < 4 ? -1 : 1;
  DEBUG("Set Printout");
  gMinuit->mnexcm("SET PRI", arglist, 1, ierflg);
  if (ierflg) {
    ERROR("Could not set printout");
    return;
  }
  if (gLogLevel < 2) {
    gMinuit->mnexcm("SET NOWarnings", arglist, 1, ierflg);
  }
  if (ierflg) {
    ERROR("Could not disable warnings");
    return;
  }

  // clear memory
  DEBUG("Clearing Memory");
  gMinuit->mnexcm("CLEAR", arglist, 0, ierflg);
  if (ierflg) {
    ERROR("Could not clear memory");
    return;
  }

  // set error definition for n-sigma deviation (likelihood: order^2/2)
  gMinuit->SetErrorDef(order*order*0.5);
  gMinuit->SetMaxIterations(500);

  // set strategy
  DEBUG("Set Strategy");
  arglist[0] = 2; // very reliable calculation
  gMinuit->mnexcm("SET STR", arglist, 1, ierflg);
  if (ierflg) {
    ERROR("Could not set strategy");
    return;
  }

  // tell Minuit about my parameters
  TString name = "xsratio";
  // expected error
  Double_t experror = gFitData->Integral()-gFitBackground->Integral();
  if (experror > 1)
    experror = 1./TMath::Sqrt(experror);
  else {
    WARNING("Low number of events");
    experror = 1;
  }
  // physical range from 0 to 100
  gMinuit->mnparm(0, name, 1, experror, 0, 100, ierflg);
  if (ierflg) {
    ERROR("Could not set parameter");
    return;
  }

  // initialize data
  DEBUG("Initializing data");
  arglist[0] = 1;
  gMinuit->mnexcm("CALL FCN", arglist, 1, ierflg);
  if (ierflg) {
    ERROR("Could not initialize data");
    return;
  }

  // now minimize
  DEBUG("Minimize");
  gMinuit->mnexcm("MINIMIZE", arglist, 0, ierflg);
  if (ierflg) {
    ERROR("MINIMIZE failed");
//      return;
  }
  DEBUG("Computing errors");
  gMinuit->mnexcm("MINOS", arglist, 0, ierflg);
  if (ierflg) {
    ERROR("MINOS failed");
//      return;
  }

  // get value
  DEBUG("Getting value and error");
  Double_t val, error, end1, end2;
  Int_t ivarbl;
  gMinuit->mnpout(0, name, val, error, end1, end2, ivarbl);
  INFO(name.Data() << ":" 
       << setw(10) << setprecision(5) << val
       << " +- " << setw(10) << setprecision(5) << error);

  // getting minos errors
  DEBUG("Getting minos errors");
  Double_t eplus, eminus, globcc;
  gMinuit->mnerrs(0, eplus, eminus, error, globcc);
  INFO(name.Data() << ":" 
       << setw(10) << setprecision(5) << val
       << " + " << setw(10) << setprecision(5) << eplus
       << " - " << setw(10) << setprecision(5) << eminus
    );
  INFO("global correlation coefficient " 
       << setw(10) << setprecision(5) << globcc);

  // save computed variables
  xs[0] = val;     // fitted value
  xs[1] = error;   // parabolic error
  xs[2] = eminus;  // asymmetric error: negative
  xs[3] = eplus;   // asymmetric error: positive

  if (draw) {
    // now we are ready to plot the curve
    Double_t  likelihood[101];
    Double_t  xval[101];
    Double_t vx, vy;
    Int_t    npar = 1;
    Double_t  xmean = (val+eminus) + (eplus-eminus)/2;
    for (Int_t i = -50; i <= 50; i++) {
      vx = (xmean + ((eplus-eminus)*i*0.51) / 50);
      BinnedPoissonianLikelihood(npar, 0, vy, & vx, 4);
      xval[i+50] = vx;
      likelihood[i+50] = vy;
    }

    // plot -Delta(log(likelihood))
    Double_t minimum = 1E99;
    for (Int_t i = -50; i <= 50; i++) {
      if (likelihood[i+50] < minimum) {
	minimum = likelihood[i+50];
      }
    }
    for (Int_t i = -50; i <= 50; i++) {
      likelihood[i+50] -= minimum;
    }

    TCanvas * cl = new TCanvas("cl", "Likelihood curve", 10, 10, 500, 500);
    setopt(cl);
    TGraph * graph = new TGraph(101, xval, likelihood);
    setopt(graph);
    graph->SetTitle("");
    graph->Draw("AC");
    setopt(graph->GetHistogram());
    graph->GetHistogram()->SetXTitle("f = #sigma_{ZZ}/#sigma^{SM}_{ZZ}");
    graph->GetHistogram()->SetYTitle("-ln L + ln L_{max}");
    graph->GetHistogram()->GetXaxis()->CenterTitle();
    TLatex * text = new TLatex(xmean,
			       graph->GetHistogram()->GetMaximum()*0.99, 
			       "e^{+}e^{-} #rightarrow ZZ #rightarrow l^{+}l^{-}q#bar{q}");
    text->SetTextAlign(23); // top, centered
    text->Draw();
    Text_t etext[256];
    if (gPeriod[gStart] == gPeriod[gEnd]) {
      snprintf(etext, 256, "%s", gPeriod[gStart]);
    }
    else {
      snprintf(etext, 256, "%s - %s", gPeriod[gStart], gPeriod[gEnd]);
    }
    TText * text1 = new TText(xmean, 
			      graph->GetHistogram()->GetMaximum()*0.93,
			      etext);
    text1->SetTextAlign(23); // top, centered
    text1->Draw();
    cl->Update();
    // error
    TLine * l = new TLine(xval[0], 0.5, xval[100], 0.5);
    l->SetLineStyle(3);
    l->SetLineWidth(2);
    l->SetLineColor(kRed);
    l->Draw();
    TArrow * arr = new TArrow(val+eminus, 0.5, val+eminus, 0.);
    arr->SetLineStyle(3);
    arr->SetLineColor(kBlue);
    arr->SetLineWidth(2);
    arr->Draw();
    arr = new TArrow(val+eplus, 0.5, val+eplus, 0.);
    arr->SetLineStyle(3);
    arr->SetLineColor(kBlue);
    arr->SetLineWidth(2);
    arr->Draw();
  }

  // tell FCN to finalize the fit
  arglist[0] = 3;
  gMinuit->mnexcm("CALL FCN", arglist, 1, ierflg);
  if (ierflg) {
    ERROR(" Could not finish the fit\n");
    return;
  }
}

void Likelihood_three(Int_t    & npar,  // # of internal parameters
		      Double_t * grad,  // array of first derivatives
		      Double_t & fval,  // the function value
		      Double_t * xvar,  // the input variable values
		      Int_t      iflag) // status flag
{
  // compute negative logarithm of likelihood for binned poissonian data.
  // Accept three variables for the tree selections: 
  // They are the relation of the measured ZZ cross-section to the SM
  // cross-section.

  fval = 1E9; // Always return some value

  if (npar != 3) {
    ERROR("Likelihood_three works only with three parameters!");
  }

  if (iflag == 1) {
    // Initialization
    DEBUG("BinnedPoissonianLikelihood() Initialization");
  } else if (iflag == 2) {
    // compute gradient of function value fval
    grad = 0;
    ERROR("No gradient computation in BinnedPoissonianLikelihood()");
  } else if (iflag == 3) {
    // Things to do after fit has finished
  } else {
    // initialize local variables
    Double_t pval; // partial value
    // compute function value fval
    fval = 0;
    // loop over histos
    for (Int_t j = 0; j < 3; j++) {
      // loop over bins
      for (Int_t i = 1; i <= gData[j]->GetNbinsX(); i++) {
	// evade log(zero)
	pval = gBack[j]->GetBinContent(i)
	  + xvar[0] * gSign[j][0]->GetBinContent(i)
	  + xvar[1] * gSign[j][1]->GetBinContent(i)
	  + xvar[2] * gSign[j][2]->GetBinContent(i);
	  if (pval > 0) {
	    fval += pval - gData[j]->GetBinContent(i) * TMath::Log(pval);
	  }
	  else if (pval < 0)
	    WARNING(" log(<0) in Likelihood_three()\n");
      }
    }
  }
}

void Likelihood3(Int_t    & npar,  // # of internal parameters
		 Double_t * grad,  // array of first derivatives
		 Double_t & fval,  // the function value
		 Double_t * xvar,  // the input variable values
		 Int_t      iflag) // status flag
{
  // compute negative logarithm of likelihood for binned poissonian data.
  // Accept three variables for the tree selections: 
  // They are the relation of the measured ZZ cross-section to the SM
  // cross-section.

  fval = 1E9; // Always return some value

  if (npar != 3) {
    ERROR("Likelihood_3 works only with three parameters!");
  }

  if (iflag == 1) {
    // Initialization
    DEBUG("BinnedPoissonianLikelihood() Initialization");
  } else if (iflag == 2) {
    // compute gradient of function value fval
    grad = 0;
    ERROR("No gradient computation in BinnedPoissonianLikelihood()");
  } else if (iflag == 3) {
    // Things to do after fit has finished
  } else {
    // initialize local variables
    Double_t pval; // partial value
    // compute function value fval
    fval = 0;
    // loop over bins
    for (Int_t i = 1; i <= gFitData->GetNbinsX(); i++) {
      // evade log(zero)
      pval = gFitBackground->GetBinContent(i)
	+ xvar[0] * gFitSig[0]->GetBinContent(i)
	+ xvar[1] * gFitSig[1]->GetBinContent(i)
	+ xvar[2] * gFitSig[2]->GetBinContent(i);
      if (pval > 0) {
	fval += pval - gFitData->GetBinContent(i) * TMath::Log(pval);
      }
      else if (pval < 0)
	ERROR("log(<0) in Likelihood()");
    }
  }
}

void xsratio3(Double_t xsval[3], Double_t xserrup[3], Double_t xserrlow[3])
{
  // fit all three selection cross sections with covariance matrix
  for (Int_t j = 0; j < 3; j++) {
    if (gFitSig[j] == 0) {
      ERROR("global fit histogram gFitSig[" << j << "] is not initialized");
      return;
    }
  }
  if (gFitBackground == 0) {
    ERROR("gFitBackground == 0");
    return;
  }
  if (gFitData == 0) {
    ERROR("gFitData == 0");
    return;
  }

  // create new Minuit with a maximum of five parameters
  if (gMinuit == 0) {
    INFO("Creating MINUIT with 5 parameters");
    gMinuit = new TMinuit(5);
  }

  // we need this function to be evaluated
  gMinuit->SetFCN(Likelihood3);

  Double_t arglist[10]; // argument list
  for (Int_t i = 0; i < 10; i++) {
    arglist[i] = 0;
  }
  Int_t ierflg;         // error return code from Minuit
  // set printout
  arglist[0] = gLogLevel < 3 ? -1 : 0;
  DEBUG("Set Printout");
  gMinuit->mnexcm("SET PRI", arglist, 1, ierflg);
  if (ierflg) {
    ERROR("Could not set printout");
    return;
  }

  // clear memory
  DEBUG("Clearing Memory");
  gMinuit->mnexcm("CLEAR", arglist, 0, ierflg);
  if (ierflg) {
    ERROR("Could not clear memory");
    return;
  }

  // set error definition for one-sigma deviation in three dimensions
  gMinuit->SetErrorDef(1.763);
  gMinuit->SetMaxIterations(500);

  // set strategy
  DEBUG("Set Strategy");
  arglist[0] = 2; // very reliable calculation
  gMinuit->mnexcm("SET STR", arglist, 1, ierflg);
  if (ierflg) {
    ERROR("Could not set strategy");
    return;
  }

  // tell Minuit about my parameters
  TString name0 = "r1";
  // expected error
  Double_t experror = 1./TMath::Sqrt(gFitData->Integral());
  Double_t expxs = (gFitData->Integral()-gFitBackground->Integral())/gFitSig[0]->Integral();
  if (expxs <= 0) {
    ERROR("zero cross-section expected");
    expxs = 0.1;
  }
  // physical range from 0 to 100
  gMinuit->mnparm(0, name0, expxs, experror, 0, 100, ierflg);
  if (ierflg) {
    ERROR("Could not set parameter");
    return;
  }
  TString name1 = "r2";
  // physical range from 0 to 100
  gMinuit->mnparm(1, name1, expxs, experror, 0, 100, ierflg);
  if (ierflg) {
    ERROR("Could not set parameter");
    return;
  }
  TString name2 = "r3";
  // physical range from 0 to 100
  gMinuit->mnparm(2, name2, expxs, experror, 0, 100, ierflg);
  if (ierflg) {
    ERROR("Could not set parameter");
    return;
  }

  // initialize data
  DEBUG("Initializing data");
  arglist[0] = 1;
  gMinuit->mnexcm("CALL FCN", arglist, 1, ierflg);
  if (ierflg) {
    ERROR("Could not initialize data");
    return;
  }

  // now minimize all three
  DEBUG("Minimize");
  gMinuit->mnexcm("MINIMIZE", arglist, 0, ierflg);
  if (ierflg) {
    ERROR("Minuit MINIMIZE");
  }

//   DEBUG("Computing errors");
//   gMinuit->mnexcm("MINOS", arglist, 0, ierflg);
//   if (ierflg) {
//     ERROR("computing minos errors");
// //      return;
//   }

  // get value
  DEBUG("Getting value and error");
  Double_t val[3], error, end1, end2;
  Int_t ivarbl;
  gMinuit->mnpout(0, name0, val[0], error, end1, end2, ivarbl);
  INFO(name0.Data() << ":" 
       << setw(10) << setprecision(5) << val[0]
       << " +- " << setw(10) << setprecision(5) << error);
  gMinuit->mnpout(1, name1, val[1], error, end1, end2, ivarbl);
  INFO(name1.Data() << ":" 
       << setw(10) << setprecision(5) << val[1]
       << " +- " << setw(10) << setprecision(5) << error);
  gMinuit->mnpout(2, name2, val[2], error, end1, end2, ivarbl);
  INFO(name2.Data() << ":" 
       << setw(10) << setprecision(5) << val[2]
       << " +- " << setw(10) << setprecision(5) << error);

  // getting minos errors
  DEBUG("Getting minos errors");
  Double_t eplus, eminus, globcc;
  gMinuit->mnerrs(0, eplus, eminus, error, globcc);
  xsval[0] = val[0]; xserrup[0] = eplus; xserrlow[0] = eminus;
  INFO(name0.Data() << ":" 
       << setw(10) << setprecision(5) << val[0]
       << " + " << setw(10) << setprecision(5) << eplus
       << " - " << setw(10) << setprecision(5) << eminus
    );
  DEBUG("global correlation coefficient " << globcc);
  gMinuit->mnerrs(1, eplus, eminus, error, globcc);
  xsval[1] = val[1]; xserrup[1] = eplus; xserrlow[1] = eminus;
  INFO(name1.Data() << ":" 
       << setw(10) << setprecision(5) << val[1]
       << " + " << setw(10) << setprecision(5) << eplus
       << " - " << setw(10) << setprecision(5) << eminus
    );
  DEBUG("global correlation coefficient " << globcc);
  gMinuit->mnerrs(2, eplus, eminus, error, globcc);
  xsval[2] = val[2]; xserrup[2] = eplus; xserrlow[2] = eminus;
  INFO(name2.Data() << ":" 
       << setw(10) << setprecision(5) << val[2]
       << " + " << setw(10) << setprecision(5) << eplus
       << " - " << setw(10) << setprecision(5) << eminus
    );
  DEBUG("global correlation coefficient " << globcc);
}

void xsratio_three(Double_t xsval[3], Double_t xserrup[3], Double_t xserrlow[3])
{
  // fit all three selection cross sections with covariance matrix
  for (Int_t j = 0; j < 3; j++) {
    if (!(gSign[j] || gBack[j] || gData[j])) {
      ERROR("global fit histograms " << j << " are not initialized");
      return;
    }
  }

  // create new Minuit with a maximum of five parameters
  if (gMinuit == 0) {
    INFO("Creating MINUIT with 5 parameters");
    gMinuit = new TMinuit(5);
  }

  // we need this function to be evaluated
  gMinuit->SetFCN(Likelihood_three);

  Double_t arglist[10]; // argument list
  for (Int_t i = 0; i < 10; i++) {
    arglist[i] = 0;
  }
  Int_t ierflg;         // error return code from Minuit
  // set printout
  arglist[0] = gLogLevel < 3 ? -1 : 0;
  DEBUG("Set Printout");
  gMinuit->mnexcm("SET PRI", arglist, 1, ierflg);
  if (ierflg) {
    ERROR("Could not set printout");
    return;
  }

  // clear memory
  DEBUG("Clearing Memory");
  gMinuit->mnexcm("CLEAR", arglist, 0, ierflg);
  if (ierflg) {
    ERROR("Could not clear memory");
    return;
  }

  // set error definition for one-sigma deviation in three dimensions
  gMinuit->SetErrorDef(1.763);
  gMinuit->SetMaxIterations(500);

  // set strategy
  DEBUG("Set Strategy");
  arglist[0] = 2; // very reliable calculation
  gMinuit->mnexcm("SET STR", arglist, 1, ierflg);
  if (ierflg) {
    ERROR("Could not set strategy");
    return;
  }

  // tell Minuit about my parameters
  TString name0 = "r_ee";
  // expected error
  Double_t experror = gData[0]->Integral()-gBack[0]->Integral();
  experror = 1./TMath::Sqrt(experror);
  if (experror < 1)
    ERROR("Low number of events");
  Double_t expxs = (gData[0]->Integral()-gBack[0]->Integral())/gSign[0][0]->Integral();
  if (expxs <= 0) {
    ERROR("zero cross-section expected");
    expxs = 0.1;
  }
  // physical range from 0 to 100
  gMinuit->mnparm(0, name0, expxs, experror, 0, 100, ierflg);
  if (ierflg) {
    ERROR("Could not set parameter");
    return;
  }
  TString name1 = "r_mm";
  // expected error
  experror = gData[1]->Integral()-gBack[1]->Integral();
  experror = 1./TMath::Sqrt(experror);
  if (experror < 1)
    ERROR("Low number of events");
  expxs = (gData[1]->Integral()-gBack[1]->Integral())/gSign[1][1]->Integral();
  if (expxs <= 0) {
    ERROR("zero cross-section expected");
    expxs = 0.1;
  }
  // physical range from 0 to 100
  gMinuit->mnparm(1, name1, expxs, experror, 0, 100, ierflg);
  if (ierflg) {
    ERROR("Could not set parameter");
    return;
  }
  TString name2 = "r_tt";
  // expected error
  experror = gData[2]->Integral()-gBack[2]->Integral();
  experror = 1./TMath::Sqrt(experror);
  if (experror < 1)
    ERROR("Low number of events");
  expxs = (gData[2]->Integral()-gBack[2]->Integral())/gSign[2][2]->Integral();
  if (expxs <= 0) {
    ERROR("zero cross-section expected");
    expxs = 0.1;
  }
  // physical range from 0 to 100
  gMinuit->mnparm(2, name2, expxs, experror, 0, 100, ierflg);
  if (ierflg) {
    ERROR("Could not set parameter");
    return;
  }

  // initialize data
  DEBUG("Initializing data");
  arglist[0] = 1;
  gMinuit->mnexcm("CALL FCN", arglist, 1, ierflg);
  if (ierflg) {
    ERROR("Could not initialize data");
    return;
  }

//    // fix parameters...
//    arglist[0] = 1;
//    arglist[1] = 2;
//    arglist[2] = 3;
//    gMinuit->mnexcm("FIX", arglist, 3, ierflg);
//    if (ierflg) {
//      printf("Error fixing parameters\n");
//      Output();
//  //      return;
//    }

//    // now minimize one at a time ...
//    for (Int_t j = 0; j < 3; j++) {
//      // release this parameter
    
//    }

  // release all parameters

  // now minimize all three
  DEBUG("Minimize");
  gMinuit->mnexcm("MINIMIZE", arglist, 0, ierflg);
  if (ierflg) {
    ERROR("Minuit MINIMIZE");
  }
//   DEBUG("Computing errors");
//   gMinuit->mnexcm("MINOS", arglist, 0, ierflg);
//   if (ierflg) {
//     printf("Error computing cross-section\n");
// //      return;
//   }

  // get value
  DEBUG("Getting value and error");
  Double_t val[3], error, end1, end2;
  Int_t ivarbl;
  gMinuit->mnpout(0, name0, val[0], error, end1, end2, ivarbl);
  INFO(name0.Data() << ":" 
       << setw(10) << setprecision(5) << val[0]
       << " +- " << setw(10) << setprecision(5) << error);
  gMinuit->mnpout(1, name1, val[1], error, end1, end2, ivarbl);
  INFO(name1.Data() << ":" 
       << setw(10) << setprecision(5) << val[1]
       << " +- " << setw(10) << setprecision(5) << error);
  gMinuit->mnpout(2, name2, val[2], error, end1, end2, ivarbl);
  INFO(name2.Data() << ":" 
       << setw(10) << setprecision(5) << val[2]
       << " +- " << setw(10) << setprecision(5) << error);

  // getting minos errors
  DEBUG("Getting minos errors");
  Double_t eplus, eminus, globcc;
  gMinuit->mnerrs(0, eplus, eminus, error, globcc);
  xsval[0] = val[0]; xserrup[0] = eplus; xserrlow[0] = eminus;
  INFO(name0.Data() << ":" 
       << setw(10) << setprecision(5) << val[0]
       << " + " << setw(10) << setprecision(5) << eplus
       << " - " << setw(10) << setprecision(5) << eminus
    );
  DEBUG("global correlation coefficient " << globcc);
  gMinuit->mnerrs(1, eplus, eminus, error, globcc);
  xsval[1] = val[1]; xserrup[1] = eplus; xserrlow[1] = eminus;
  INFO(name1.Data() << ":" 
       << setw(10) << setprecision(5) << val[1]
       << " + " << setw(10) << setprecision(5) << eplus
       << " - " << setw(10) << setprecision(5) << eminus
    );
  DEBUG("global correlation coefficient " << globcc);
  gMinuit->mnerrs(2, eplus, eminus, error, globcc);
  xsval[2] = val[2]; xserrup[2] = eplus; xserrlow[2] = eminus;
  INFO(name2.Data() << ":" 
       << setw(10) << setprecision(5) << val[2]
       << " + " << setw(10) << setprecision(5) << eplus
       << " - " << setw(10) << setprecision(5) << eminus
    );
  DEBUG("global correlation coefficient " << globcc);
}

void xsection(Int_t order = 1, const Text_t * histo = "hLm5c")
{
  // compute cross-section for each of the selected energies...
  Int_t start = gStart;
  Int_t end   = gEnd;

  FILE * outfile = fopen(Form("%s.xsection", histo), "w");
  if (!outfile) {
    ERROR("could not create cross-section file " << histo << ".xsection");
    return;
  }
  Double_t xs[4];
  // save cross-section and errors here...
  for (Int_t p = start; p <= end; p++) {
    // only one energy
    period(gPeriod[p]);
    // plot histogram to get cross-section from...
    plot(histo);
    // fill global histos
    gFitSignal     = signalHisto();
    gFitBackground = backgroundHisto();
    gFitData       = dataHisto();
    // compute xs ratio
    xsratio(xs, order, kFALSE);
    // delete histos...
    delete gFitSignal;
    delete gFitBackground;
    delete gFitData;
    // get theoretical cross-section
    Double_t xstheorie = 0;
    for (Int_t i = 0; i < 3; i++) {
      xstheorie += gProcessInfo[p][i].xs;
    }
    printf("%s %10.5f %+10.5f %+10.5f\n", 
	   gPeriod[gStart], 
	   xstheorie * xs[0], -xstheorie * xs[2], xstheorie * xs[3]);
    fprintf(outfile, "%s %10.5f %+10.5f %+10.5f\n", 
	    gPeriod[gStart], 
	    xstheorie * xs[0], -xstheorie * xs[2], xstheorie * xs[3]);
  }
  fclose(outfile);
  // select old energies...
  period(gPeriod[start], gPeriod[end]);
  gFitSignal = 0;
  gFitBackground = 0;
  gFitData = 0;
  INFO("Cross-section results contained in file " << histo << " .xsection");
}

TH1D * experror(Double_t xs[4], const Int_t ntrials = 10000)
{
  // compute expected error from given gFitSignal and gFitBackground histograms.
  // ntrials MC experiments are performed.
  // mean cross-section and error are returned in xs[4]
  // the distribution of MC pseudo-measurments is returned in a histogram


  // check
  if (gFitSignal == 0 || gFitBackground == 0) {
    ERROR("experror(): gFitSignal or gFitBackground undefined");
    return 0;
  }
  if (gFitSignal->GetNbinsX() != gFitBackground->GetNbinsX()) {
    ERROR("experror(): bin number mismatch sig " <<  gFitSignal->GetNbinsX() 
	  << " and back " << gFitBackground->GetNbinsX());
    return 0;
  }

  // initialize random number generator
  gRandom->SetSeed(0);

  // distribution of cross-section
  TH1D * xsdistrib = new TH1D("xsdistrib", "p.d.f. of #sigma_{ZZ}/#sigma_{SM}",
			      300, 0, 3.0);

  // save old data
  TH1D * olddata = gFitData;
  // prepare histograms
  gFitData       = new TH1D(*gFitSignal);
  TH1D * hmc     = new TH1D(*gFitSignal);
  hmc->Add(gFitBackground);
//    printf("Cross check: %8.3f signal %8.3f background %8.3f sum %8.3f mc\n",
//  	 gFitSignal->Integral(), gFitBackground->Integral(), 
//  	 gFitSignal->Integral()+gFitBackground->Integral(),
//  	 hmc->Integral());
  // internal used array for statistics
  Double_t * values = new Double_t[ntrials];
  INFO("Starting loop...");
  // make monte-carlo pseudo-data sample
  Double_t mean = 0;
  Int_t nbins = gFitData->GetNbinsX();
  Int_t i;
  Int_t n;
  for (n = 0; n < ntrials; n++) {
    if ((n % 1000) == 0)
      if (gLogLevel > 1)
	printf("n = %d\n", n);
    // generate poisson around theorie
    for (i = 1; i <= nbins; i++) {
      gFitData->SetBinContent(i, gRandom->Poisson(hmc->GetBinContent(i)));
    }
    // compute xs ratio
    xsratio(xs, 1, kFALSE);
    // fill histogram
    xsdistrib->Fill(xs[0]);
    // also save in array
    values[n] = xs[0];
    // and add to sum to compute mean value
    mean += xs[0];
  }
  INFO("Loop ended");

  // clear memory
  delete gFitData;
  gFitData = olddata;
  delete hmc;

  // do statistics
  Int_t * sindex = new Int_t[ntrials];
  TMath::Sort(ntrials, values, sindex, kFALSE);

  // one-sigma lower limit
  Double_t limlow  = values[sindex[int(ntrials * TMath::Erfc(1./TMath::Sqrt(2.)) / 2.)]];
  // one-sigma upper limit
  Double_t limhigh = values[sindex[int(ntrials*(1-TMath::Erfc(1./TMath::Sqrt(2.)) / 2.))]];
  mean /= ntrials;

  if (gLogLevel > 2)
    printf("Mean value and error with %8d trials: %8.4f +%8.4f -%8.4f\n",
	   ntrials, mean, limhigh-mean, mean-limlow);

  // clear more memory
  delete[] sindex;
  delete[] values;

  // save for later use
  xs[0] = mean;
  xs[1] = (limhigh-limlow) / 2.;
  xs[2] = mean-limlow;
  xs[3] = limhigh-mean;

  // draw histo
  TCanvas * scanvas = new TCanvas("scanvas", "cross-section distribution", 
				   40, 40, 800, 600);
  setopt(xsdistrib);
  xsdistrib->GetXaxis()->SetTitle("#sigma_{MC}/#sigma_{SM}");
  xsdistrib->Draw();

  TLatex * l = new TLatex(3, xsdistrib->GetMaximum(),
			  Form("Mean = %8.4f +%8.4f -%8.4f",
			       mean, limhigh-mean, mean-limlow));
  l->SetTextAlign(33);
  l->Draw();
  scanvas->Update();

  return xsdistrib;
}

void xsExpectedError(Int_t ntrials)
{
  // try ntrials times to measure the cross-section with a fake sample and
  // return statistical error on this computation

  Double_t xs[4];
  gFitSignal = signalHisto();
  gFitBackground = backgroundHisto();
  TH1D * xsdistrib = experror(xs, ntrials);
  xsdistrib->SetName(Form("experror-%s-%s-%d", gPeriod[gStart], gPeriod[gEnd], ntrials));

  // save histo in file
  TFile * f = new TFile("experror.root", "UPDATE");
  xsdistrib->Write();
  f->Close();
  delete f;

  delete gFitSignal;
  delete gFitBackground;
}

void correctfactor(const char * pname, Double_t xlow, Double_t xup, bool setWeight)
{
  // determine correction factor gCorrect from current histogram
  if (gStart != gEnd) {
    ERROR("Working on multiple periods is impossible");
    return;
  }
  // now loop over processes to get "signal" and "background"
  gFitSignal = 0;
  gFitBackground = 0;
  for (Int_t process = 0; process < gMaxProcess-1; process++) {
    if (gHisto[gPadNr][process] == 0)
      continue;
    TString pn(gProcess[process].fname);
    if (pn.Contains(pname)) {
      INFO("Found signal " << pn);
      if (gFitSignal == 0) {
	gFitSignal = new TH1D(*gHisto[gPadNr][process]);
      }
      else {
	gFitSignal->Add(gHisto[gPadNr][process]);
      }
    }
    else {
      // INFO("Found background " << pn);
      if (gFitBackground == 0) {
	gFitBackground = new TH1D(*gHisto[gPadNr][process]);
	gFitBackground->SetDirectory(0);
      }
      else
	gFitBackground->Add(gHisto[gPadNr][process]);
    }
  }
  gFitData = new TH1D(*gHisto[gPadNr][gMaxProcess-1]);
  // only fit in requested range by zeroing out signal values
  Int_t binlow, binup;
  findbins(xlow, xup, binlow, binup);
  for (Int_t i = 1; i < binlow; i++) {
    gFitBackground->SetBinContent(i, 0);
    gFitSignal->SetBinContent(i, 0);
    gFitData->SetBinContent(i, 0);
  }
  for (Int_t i = binup; i <= gFitData->GetNbinsX(); i++) {
    gFitBackground->SetBinContent(i, 0);
    gFitSignal->SetBinContent(i, 0);
    gFitData->SetBinContent(i, 0);
  }
  // now fit cross-section
  Double_t xs[4];
  xsratio(xs);
  for (Int_t process = 0; process < gMaxProcess-1; process++) {
    TString pn(gProcess[process].fname);
    if (!pn.Contains(pname))
      continue;
    Double_t oldweight = gProcessInfo[gStart][process].weight;
    Double_t newweight = xs[0] * oldweight;
    printf("Period %s, process %s: old weight = %8.5f, new weight = %8.5f %+8.5f %+8.5f\n",
	   gPeriod[gStart], gProcess[process].fname,
	   oldweight, 
	   newweight, 
	   xs[3] * oldweight,
	   xs[2] * oldweight);
    if (setWeight) {
      gProcessInfo[gStart][process].weight = newweight;
      cout << "Setting weight done!" << endl;
    }
  }
  
  delete gFitData;
  gFitData = 0;
  delete gFitBackground;
  gFitBackground = 0;
  delete gFitSignal;
  gFitSignal = 0;
}

void allcorrect(Text_t * hname, Text_t * pname, Double_t xlow, Double_t xup)
{
  Int_t oldstart = gStart;
  Int_t oldend   = gEnd;
  for (Int_t p = oldstart; p <= oldend; p++) {
    period(gPeriod[p]);
    plot(hname);
    correctfactor(pname, xlow, xup);
  }
  period(gPeriod[oldstart], gPeriod[oldend]);
}

void getscale()
{
  gFitBackground = 0;

  // loop over processes to get "signal" and "background"
  for (Int_t process = 0; process < gMaxProcess-1; process++) {
    if (gHisto[gPadNr][process] == 0)
      continue;
    TString pn(gProcess[process].fname);
    if (pn.Contains("dyll")) {
      cout << "Found signal 0: " << pn << endl;
      if (gFitSig[0] == 0) {
	gFitSig[0] = new TH1D(*gHisto[gPadNr][process]);
      }
      else {
	gFitSig[0]->Add(gHisto[gPadNr][process]);
      }
    }
    else if (pn.Contains("qcd")) {
      cout << "Found signal 1: " << pn << endl;
      if (gFitSig[1] == 0) {
	gFitSig[1] = new TH1D(*gHisto[gPadNr][process]);
      }
      else {
	gFitSig[1]->Add(gHisto[gPadNr][process]);
      }
    }
    else if (pn.Contains("ttjets")) {
      cout << "Found signal 2: " << pn << endl;
      if (gFitSig[2] == 0) {
	gFitSig[2] = new TH1D(*gHisto[gPadNr][process]);
      }
      else {
	gFitSig[2]->Add(gHisto[gPadNr][process]);
      }
    }
    else {
      // INFO("Found background " << pn);
      if (gFitBackground == 0) {
	gFitBackground = new TH1D(*gHisto[gPadNr][process]);
	gFitBackground->SetDirectory(0);
      }
      else
	gFitBackground->Add(gHisto[gPadNr][process]);
    }
  }
  gFitData = new TH1D(*gHisto[gPadNr][gMaxProcess-1]);

  double xsval[3], xserrup[3], xserrlow[3];
  xsratio3(xsval, xserrup, xserrlow);
  for (Int_t process = 0; process < gMaxProcess-1; process++) {
    TString pn(gProcess[process].fname);
    Int_t index = -1;
    if (pn.Contains("dyll")) {
      index = 0;
    }
    else if (pn.Contains("qcd")) {
      index = 1;
    }
    else if (pn.Contains("ttjets")) {
      index = 2;
    }
    else
      continue;

    Double_t oldweight = gProcessInfo[gStart][process].weight;
    Double_t newweight = xsval[index] * oldweight;
    printf("Period %s, process %s: old weight = %8.5f, new weight = %8.5f %+8.5f %+8.5f\n",
	   gPeriod[gStart], gProcess[process].fname,
	   oldweight, 
	   newweight, 
	   xserrup[index] * oldweight,
	   xserrlow[index] * oldweight);
    if (true) {
      gProcessInfo[gStart][process].weight = newweight;
      cout << "Setting weight done!" << endl;
    }
  }
}
