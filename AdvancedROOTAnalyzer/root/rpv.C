#include <iostream>
#include <iomanip>

#include "rpv.h"

#include "TMath.h"
#include "TROOT.h"
#include "TFile.h"
#include "TLatex.h"
#include "TF1.h"
#include "TH2D.h"
#include "TH1D.h"
#include "plot.h"
#include "fakerate.h"

using namespace std;

/** Compute the value of lambda_prime_211 at any scale, given a value val0 at scale mu0
 *
 * Calculation see Dreiner et al, Phys. Rev. D 75, 035003 (2007)
 *
 * @param double log_mu  The log(scale) at which to compute the value of lambda_prime_211
 * @param double log_mu0 The log(scale) at which the value val0 has been given
 * @param double val0 The value of lambda_prime_211 at the scale mu0
 */
double lambda_prime_211(double log_mu, double log_mu0, double val0)
{
  double C_F = 4./3.;      // color factor
  double alpha_s = 0.1184; // alpha_s(M_Z)
  double pi = TMath::Pi();

  double beta_qcd = 3*C_F; // running with only QCD corrections
  double beta_susy = -C_F; // running with only SUSY corrections
  double beta_tot = beta_qcd+beta_susy; // running with both QCD and SUSY QCD corrections

  return val0/(1+alpha_s/(4*pi)*beta_tot*2*(log_mu-log_mu0));
  // next line is numerically unstable...
  // return val0/(1+alpha_s/(4*pi)*beta_tot*2*(TMath::Log(mu)-TMath::Log(mu0)));
}

void running_lambda_prime_211()
{
  MakeCanvas(1,1);
  cd(1);
  double q_min = 2;
  double q_max = TMath::Log(1.84882e+16)/TMath::Log(10.); // Softsusy 3.17: 1.84882e+16
  cout << "qmax = " << q_max << endl;
  const int nMax = 100;
  double x[nMax];
  double y[nMax];
  TH2F * hframe = new TH2F("hframe", "frame", 1, 1, 17, 1, 0, 0.1);
  setopt(hframe);
  hframe->SetXTitle("log_{10}(#mu/GeV)");
  hframe->SetYTitle("#lambda'_{211}");
  hframe->Draw();
  for (int n = 0; n < nMax; n++) {
    x[n] = (q_max-q_min)*n/(nMax-1)+q_min;
    y[n] = lambda_prime_211(x[n], q_max, 1E-2);
  }
  TGraph * gr = new TGraph(nMax, x, y);
  gr->Draw("lpsame");
  gPad->Print("rpv.pdf");
  cout << "Value of lambda'_211 at electroweak scale: " << y[0] << endl;
  cout << "Value of lambda'_211 at GUT scale: " << y[nMax-1] << endl;
}

Bool_t read_xsfile(const char * fname, int nMax,
		   int p1, int p2, double & lambda, double & sqrts,
		   int & n, double m[], double xs_lo[], double xs_nlo[])
{
  FILE * xsfile = fopen(fname, "r");
  if (xsfile == 0) {
    ERROR("could not open file " << fname);
    return kFALSE;
  }
  char buffer[256];
  double mS, Mlr;
  double LO, NLO, k;
  double LO_c, NLO_c, k_c;
  double LO_tot, NLO_tot;
  n = 0;
  while (fgets(buffer, 256, xsfile)) {
    if (n >= nMax) {
      cerr << "File contains more entries than array (" << nMax << " - enlarge array" << endl;
      break;
    }
    if (buffer[0] == '#')
      continue;
    int i1, i2;
    if (sscanf(buffer, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	       & i1, & i2,
	       & sqrts,
	       & mS, & Mlr, & lambda,
	       & LO, & NLO, & k,
	       & LO_c, & NLO_c, & k_c,
	       & LO_tot, & NLO_tot) != 14) {
      cerr << "Error in format in this line: " << buffer << endl;
      continue;
    }
    if ((TMath::Abs(i1) == TMath::Abs(p1) && TMath::Abs(i2) == TMath::Abs(p2)) ||
	(TMath::Abs(i1) == TMath::Abs(p2) && TMath::Abs(i2) == TMath::Abs(p1))) {
      xs_lo[n] = LO_tot;
      xs_nlo[n] = NLO_tot;
      m[n] = mS;
      n++;
    }
  }
  return kTRUE;
}

void xscurve(const char * fname = "scan_7TeV_ud_smuon.txt")
{
  double lambda, sqrts;
  const int nMax = 1000;
  double m[nMax];
  double xs_lo[nMax];
  double xs_nlo[nMax];
  int n = 0;
  double xs_min = 1E9;
  double xs_max = 0;
  double m_min = 1E9;
  double m_max = 0;
  char quarks[6] = { 'd', 'u', 's', 'c', 'b', 't' };
  int p1, p2;
  p1 = 1; p2 = 2; // ud
  if (!read_xsfile(fname, nMax, p1, p2, lambda, sqrts, n, m, xs_lo, xs_nlo)) {
    cerr << "Could not read file " << fname << endl;
    return;
  }
  cout << "Found " << n << " data points in file " << fname << endl;
  for (int i = 0; i < n; i++) {
    xs_max = TMath::Max(TMath::Max(xs_max, xs_lo[i]), xs_nlo[i]);
    xs_min = TMath::Min(TMath::Min(xs_min, xs_lo[i]), xs_nlo[i]);

    m_min = TMath::Min(m_min, m[i]);
    m_max = TMath::Max(m_max, m[i]);
  }

  TCanvas * c1 = new TCanvas("c1", "Single slepton Cross-section",
			     (Int_t) (400*TMath::Sqrt(2.)), 400);
  setopt(c1);
  gPad->SetLogy();
  TH2F * hframe = new TH2F("hframe", "frame", 1, m_min, m_max, 1, xs_min, xs_max);
  setopt(hframe);
  hframe->SetXTitle("m(#bf{#tilde{#mu}^{#pm}}) [GeV]");
  hframe->SetYTitle(Form("#sigma(%c %c #rightarrow #bf{#tilde{#mu}^{#pm}}) [pb]",
			 quarks[TMath::Abs(p1)-1], quarks[TMath::Abs(p2)-1]));
  hframe->Draw();

  TGraph * gr_lo = new TGraph(n, m, xs_lo);
  // setopt(gr_lo);
  gr_lo->SetLineColor(kRed);
  gr_lo->SetLineStyle(kDashed);
  gr_lo->SetLineWidth(3.);
  gr_lo->Draw("lsame");

  TGraph * gr_nlo = new TGraph(n, m, xs_nlo);
  // setopt(gr_nlo);
  gr_nlo->SetLineColor(kBlue);
  gr_nlo->SetLineStyle(kSolid);
  gr_nlo->SetLineWidth(3.);
  gr_nlo->Draw("lsame");

  TLegend * leg = new TLegend(0.6, 0.76, 0.93, 0.91, Form("#sqrt{s} = %.0f TeV, #lambda'_{211} = %.3f", sqrts/1000., lambda));
  setopt(leg);
  leg->AddEntry(gr_lo, "LO", "l");
  leg->AddEntry(gr_nlo, "NLO", "l");
  leg->Draw();
}

/** Return linear interpolation between points (x0, val0) and (x1, val1) at given x */
double interpolate_linear(double x0, double val0, double x1, double val1, double x)
{
  if (x0 > x1)
    THROW("Integrity violation: x0 > x1");
  if (x > x1 || x < x0)
    THROW("Given x not in interval [x0,x1]");

  double m = (val1-val0)/(x1-x0);
  return val0 + m * (x - x0);
}

void xs_lm1()
{
  double m_smuon = 185.162;
  double m_sneutrino = 167.315;

  double xs_lo  = interpolate_linear(185.0, 1.112930, 190.0, 1.020013, m_smuon);
  double xs_nlo = interpolate_linear(185.0, 1.496843, 190.0, 1.374188, m_smuon);

  xs_lo  += interpolate_linear(165.0, 1.604928, 170.0, 1.457219, m_sneutrino);
  xs_nlo += interpolate_linear(165.0, 2.138423, 170.0, 1.945395, m_sneutrino);

  cout << "xs@lm1 (unscaled):" << endl;
  cout << "xs_lo  = " << xs_lo << " pb" << endl;
  cout << "xs_nlo = " << xs_nlo << " pb" << endl;

  double lp211_gut = 0.01;
  double lp211_ewk = 0.0301796;
  double scale = (lp211_ewk*lp211_ewk)/(lp211_gut*lp211_gut);
  xs_lo *= scale;
  xs_nlo *= scale;

  cout << "xs@lm1: (scaled)" << endl;
  cout << "xs_lo  = " << xs_lo << " pb" << endl;
  cout << "xs_nlo = " << xs_nlo << " pb" << endl;
}

/**
 * Analyze the signal and make some summary plots
 */
void signalplot(Bool_t batchmode = kTRUE)
{
  gROOT->SetBatch(batchmode);

  // plot
  cd(1);
  plot("Sig_nMuon");
  zoom(0, 6);
  legend();
  cd(2);
  plot("Sig_nNeutrino");
  zoom(0, 6);
  print("Sig_nMuNeutrino.pdf");

  // plot
  cd(1);
  plot("Sig_nQuark");
  zoom(0, 6);
  legend();
  cd(2);
  plot("Sig_nHiggs");
  zoom(0, 6);
  print("Sig_nQuarkHiggs.pdf");

  // plot
  cd(1);
  plot("Sig_SleptonMass");
  zoom(980,1100);
  legend();
  cd(2);
  plot("Sig_Resonance");
  zoom(980,1100);
  print("Sig_SleptonMassResonance.pdf");

  cd(1);
  plot("Sig_nBoson");
  legend();
  cd(2);
  plot("Sig_SleptonCharge");
  print("Sig_nBosonSleptonCharge.pdf");

  // plot
  cd(1);
  plot("Sig_EtaMu");
  rebin(2);
  legend();
  cd(2);
  plot("Sig_PhiMu");
  rebin(4);
  print("Sig_MuEtaPhi.pdf");

  // plot
  cd(1);
  plot("Sig_ptMu1");
  rebin(5);
  // zoom(0,80);
  legend(0.01, -1);
  cd(2);
  plot("Sig_ptMu2");
  rebin(5);
  // zoom(0,80);
  print("Sig_MuPt12.pdf");

  // plot
  cd(1);
  plot("Sig_EtaQuark");
  legend();
  cd(2);
  plot("Sig_PhiQuark");
  rebin(4);
  print("Sig_QuarkEtaPhi.pdf");

  // plot
  cd(1);
  plot("Sig_ptQuark1");
  rebin(5);
  // zoom(0,80);
  legend(-1);
  cd(2);
  plot("Sig_ptQuark2");
  rebin(5);
  // zoom(0,80);
  print("Sig_QuarkPt.pdf");

  // plot
  cd(1);
  plot("Sig_MuDEta");
  legend();
  cd(2);
  plot("Sig_MuDPhi");
  rebin(5);
  print("Sig_MuEtaPhi.pdf");

  // plot
  cd(1);
  plot("Sig_MuMass");
  rebin(5);
  legend();
  cd(2);
  plot("Sig_MuDPt");
  rebin(5);
  print("Sig_MuMassDPt.pdf");

  // plot
  cd(1);
  plot("Sig_JetDEta");
  legend();
  cd(2);
  plot("Sig_JetDPhi");
  rebin(5);
  print("Sig_JetEtaPhi.pdf");

  // plot
  cd(1);
  plot("Sig_JetMass");
  rebin(5);
  zoom(0,300);
  legend();
  cd(2);
  plot("Sig_JetDPt");
  rebin(5);
  print("Sig_JetMassDPt.pdf");

  // plot
  cd(1);
  plot("Sig_JetSumPt");
  rebin(5);
  legend();
  cd(2);
  plot("xxx");
  print("Sig_JetSumPt.pdf");

  // plot
  cd(1);
  plot("Sig_Mass3");
  zoom(200, 220);
  legend();
  cd(2);
  plot("Sig_Angle3");
  rebin(5);
  print("Sig_Mass3Angle3.pdf");

  // plot
  cd(1);
  plot("vtx_n");
  legend();
  cd(2);
  plot("vtx_x");
  logy();
  min(0.1);
  print("Vertex_1.pdf");

  // plot
  cd(1);
  plot("vtx_y");
  logy();
  min(0.1);
  legend();
  cd(2);
  plot("vtx_z");
  liny();
  print("Vertex_2.pdf");

  // plot
  cd(1);
  plot("vtx_ntr");
  liny();
  legend();
  cd(2);
  plot("vtx_ndof");
  print("Vertex_3.pdf");
}

void generate_pileup_histogram()
{

// Distribution used for Fall2011 MC.

  Double_t Fall2011[50] = {
    0.003388501,
    0.010357558,
    0.024724258,
    0.042348605,
    0.058279812,
    0.068851751,
    0.072914824,
    0.071579609,
    0.066811668,
    0.060672356,
    0.054528356,
    0.04919354,
    0.044886042,
    0.041341896,
    0.0384679,
    0.035871463,
    0.03341952,
    0.030915649,
    0.028395374,
    0.025798107,
    0.023237445,
    0.020602754,
    0.0180688,
    0.015559693,
    0.013211063,
    0.010964293,
    0.008920993,
    0.007080504,
    0.005499239,
    0.004187022,
    0.003096474,
    0.002237361,
    0.001566428,
    0.001074149,
    0.000721755,
    0.000470838,
    0.00030268,
    0.000184665,
    0.000112883,
    6.74043E-05,
    3.82178E-05,
    2.22847E-05,
    1.20933E-05,
    6.96173E-06,
    3.4689E-06,
    1.96172E-06,
    8.49283E-07,
    5.02393E-07,
    2.15311E-07,
    9.56938E-08
  };

  TFile * file = new TFile("Fall11_MC_pileup_truth.root", "RECREATE");
  TH1F * h1 = new TH1F("pileup", "Pileup MC Fall11 distribution", 50, 0, 50);
  for (int i = 0; i < 50; i++) {
    h1->SetBinContent(i+1, Fall2011[i]);
  }
  file->Write();
  file->Close();
  delete file;
}

/*
 * Compute best value for smearing of energy resolution by calculating a chi2
 * for various scaling parameters and choosing those with lowest chi2 value.
 */

void compute_energy_resolution()
{
  // get histograms from file
  MakeCanvas(2,2);
  plot2("pfmet_scaled");
  TH2D * hData = gHisto2[gMaxProcess-1];
  if (hData == 0) {
    cerr << "ERR: no data found" << endl;
    return;
  }
  // Create empty result histogram, one bin for each factor
  // take binning from existing data histogram
  TH1D * hResult = hData->ProjectionY("hResult", 1, 1);
  hResult->SetTitle("#chi^2");
  hResult->Reset();
  // Create empty histogram with same bining as existing histograms
  TH2D * hBackground = new TH2D(*hData);
  hBackground->Reset();
  // add all background MC contributions together, excluding the not stacked ones
  for (int process = 0; process < gMaxProcess-1; process++) {
    Int_t i = gOrder[gPadNr][process];
    TH2D * histo2 = gHisto2[i];
    if (histo2 == 0) {
      continue;
    }
    if (!gProcess[i].stack) {
      cout << "INF: not adding " << gProcess[i].fname << endl;
      continue;
    }
    cout << "INF: Adding " << gProcess[i].fname << endl;
    hBackground->Add(histo2);
  }
  // for data, always take scale 0 (first bin)
  TH1D * pData = hData->ProjectionX("pData_px", 1, 1);
  INFO("scale factor for data: " << hData->GetYaxis()->GetBinLowEdge(1)+hData->GetYaxis()->GetBinWidth(1)/2.);
  cd(1);
  pData->Draw();
  cd(2);

  Int_t ndof = 0;
  Int_t nlow = 0;
  // loop over all possible scale factor values
  for (int biny = 1; biny < hData->GetNbinsY()+1; biny++) {
    double scale = hData->GetYaxis()->GetBinLowEdge(biny)+hData->GetYaxis()->GetBinWidth(biny)/2.;
    INFO("current MC scale factor is " << scale);
    // get 1d background histogram
    TH1D * pBackground = hBackground->ProjectionX("_px", biny, biny);
    pBackground->Draw();
    gPad->Update();
    // compute chi2 between data and background
    double ChiSquare = 0;
    for (int binx = 1; binx < pData->GetNbinsX()+1 && pData->GetBinLowEdge(binx) < 100; binx++) {
      // sum of all MC's ( theory )
      double temp1 = pBackground->GetBinContent(binx);

      // skip empty bins
      if (temp1 < 1E-3)
	continue;

      // theorie - data
      double temp = temp1 - pData->GetBinContent(binx);
      ChiSquare += temp * temp / temp1;
      ndof++;
      // count bins with less than five entries
      if (temp1 < 5) {
	nlow++;
      }
    }
    // "Fill" Histogram
    hResult->SetBinContent(biny, ChiSquare);
    // in order to be able to fit the histogram, one also needs bin errors
    hResult->SetBinError(biny, 1.);
  }
  cd(3);
  hResult->Draw("hist");

  // obtain minimum from fit
  TF1 * f1 = new TF1("f1", "pol2",
		     hResult->GetBinLowEdge(1),
		     hResult->GetBinLowEdge(hResult->GetNbinsX()+2));
  hResult->Fit(f1, "R0", "goff");
  f1->SetLineColor(kBlue);
  f1->SetLineWidth(2.);
  f1->Draw("same");
  // double a = f1->GetParameter(0);
  double b = f1->GetParameter(1);
  double c = f1->GetParameter(2);
  double minvalue = -b/(2*c);
  double sigma = 1./TMath::Sqrt(c);
  INFO("minimum from fit: f = " << minvalue << " +/- " << sigma
       << ", chi2 = " << f1->Eval(minvalue));

  // get minimum of histogram, take the corresponding scale factor,
  // and plot data and MC for optical comparison
  Int_t biny = hResult->GetMinimumBin();
  Double_t best_scale = hResult->GetBinLowEdge(biny)+hResult->GetBinWidth(biny)/2.;
  Double_t best_chi2 = hResult->GetBinContent(biny);
  INFO("best scale factor: " << best_scale << " with chi2 = " << best_chi2);
  TH1D * pBack = hBackground->ProjectionX("_px", biny, biny);

  cd(4);
  pBack->SetFillColor(kYellow);
  pBack->Draw("hist");
  pData->Draw("epsame");
  TLatex * l = new TLatex(0.5, 0.85, Form("scale = %5.3f", best_scale));
  l->SetNDC();
  l->Draw();
}

// fit width of energy distribution
void compute_energy_resolution_offset()
{
  plot("JER_deltae");
  TH1D * back = backgroundHisto();
//   TF1 * f1 = new TF1("f1", "[0]*TMath::BreitWigner(x,[1],[2])", -50, 50);
  TF1 * f1 = new TF1("f1", "gaus", -50, 50);
  f1->SetParameter(0, back->Integral());
  f1->SetParameter(1, back->GetMean());
  f1->SetParameter(2, back->GetRMS());
  // do a log-likelihood fit (ll), do not plot result (0), and use function range (R)
  back->Fit("f1", "R");
  f1->SetLineColor(kBlue);
  f1->SetLineWidth(2.);
  f1->Draw("same");
}


void plotlinlog(const char * name)
{
  cd(1);
  plot(name);
  liny();
  legend();

  cd(2);
  plot(name);
  logy();
  legend();

  print(Form("%s.pdf", name));
}

void mc_comparison()
{
  setup("../config/plot2.cfg");
  selection("mctest");
  MakeCanvas();

  plotlinlog("check_mumass");
  plotlinlog("check_vtxn");
  plotlinlog("check_nmuon");
  plotlinlog("check_njets");
}

// ttbar and data in one histogram
void signal_xcheck()
{
  gStyle->SetPadRightMargin(0.2);
  plot2("m_smu_chi");
  MakeCanvas(1,1);
  TH2D * hBack = 0;
  // plot background(s)
  for (int i = 0; i < gMaxProcess-1; i++) {
    if (!strncmp("ttjets", gProcess[i].fname, 6)) {
      hBack = gHisto2[i];
      hBack->SetAxisRange(0, 1000);
      hBack->Draw("colz");
    }
  }
  // // plot signal(s)
  // for (int i = 0; i < gMaxProcess-1; i++) {
  //   if (!strncmp("signal", gProcess[i].fname, 6)) {
  // 	TH2D * hSignal = gHisto2[i];
  // 	hSignal->Scale(1./50.);
  // 	hSignal->Draw("boxsame");
  //   }
  // }
  // plot data
  TH2D * hdata = gHisto2[gMaxProcess-1];
  hdata->SetMarkerStyle(5);
  hdata->SetMarkerSize(2.);
  hdata->Draw("same");

  TLegend * leg = new TLegend(0.6, 0.2, 0.75, 0.5);
  leg->AddEntry(hBack, "t#bar{t}", "f");
  // leg->AddEntry(hSignal, "signal", "l");
  leg->AddEntry(hdata, "data", "p");
  leg->Draw();
}

void limitplot(const char * fname = "XscLimitsAndErrorsRooStats.txt")
{
  gROOT->SetBatch(1);
  TH1D * h1 = new TH1D("h1", "sigma", 100, 0., 3.);

  FILE * limitfile = fopen(fname, "r");
  if (limitfile == 0) {
    ERROR("could not open file " << fname);
    return;
  }
  char buffer[256];
  double m0, m12;
  double observed, expected;
  double onesigmaup, onesigmadown;
  double twosigmaup, twosigmadown;
  while (fgets(buffer, 256, limitfile)) {
    if (buffer[0] == '#')
      continue;
    if (sscanf(buffer, "%lf & %lf & %lf & %lf & %lf & %lf & %lf & %lf \\",
	       & m0, & m12,
	       & observed, & expected,
	       & onesigmaup, & onesigmadown,
	       & twosigmaup, & twosigmadown
	  ) != 8) {
      cerr << "Error in format in this line: " << buffer << endl;
      continue;
    }
    h1->Fill((observed-expected)/(onesigmaup-expected));
  }
  h1->Draw();
  gPad->Print("limit.pdf");
}

void twoplots(const char * sel1, const char * sel2, const char * hname)
{
  selection(sel1);
  cd(1);
  plot(hname);
  selection(sel2);
  cd(2);
  plot(hname);
}

void selectionplots()
{
  MakeCanvas();
  cd(1);
  plot("cutflow");
  min(1);
  logy();
  legend();
}

void controlplots(const char * sel)
{
  selection(sel);
  plot("pfmet");
  zoom(0,150);
  legend();
  cd(2);
  plot("CR1_m_mumu");
  zoom(0,200);
  legend();
  print(Form("JER_%s.pdf", sel));
  cd(1);
  plot("CR2_m_mumu");
  zoom(0,200);
  legend();
  cd(2);
  plot("CR3_m_mumu");
  zoom(0,200);
  legend();
  print(Form("JER_%s_btag.pdf", sel));

  // final control region
  cd(1);
  plot("CR5_m_smuon");
  zoom(0,200);
  legend();
  cd(2);
  plot("CR6_m_smuon");
  zoom(0,200);
  legend();
  print(Form("CR_final_%s.pdf", sel));
}

/** From the given histogram, which is required to have 3x3 bins, extract the
 *  lower-diagonal bin entries and errors and store them in the given array of
 *  doubles, where the first entry for each bin is the content and the second
 *  entry the error. */
void get2dstatistics(const TH2D * h2, double values[6][2])
{
  // plot histograms
  // loop over all histograms
  double N = 0;
  double e2 = 0;
  if (h2 == 0) {
    WARNING("Given histogram pointer is zero");
    return;
  }
  if (h2->GetNbinsX() != 3 || h2->GetNbinsY() != 3) {
    ERROR("Wrong number of bins in histogram in get2dstatistics()");
    return;
  }
  int bin = 0;
  for (int binx = 1; binx <= h2->GetNbinsX(); binx++) {
    for (int biny = 1; biny <= binx; biny++) {
      // cout << "binx = " << binx << ", biny = " << biny << endl;
      cout << "x=" << setw(4) << h2->GetXaxis()->GetBinLowEdge(binx)
  	   << ".." << setw(4) << h2->GetXaxis()->GetBinUpEdge(binx) << ", ";
      cout << "y=" << setw(4) << h2->GetYaxis()->GetBinLowEdge(biny)
  	   << ".." << setw(4) << h2->GetYaxis()->GetBinUpEdge(biny) << ", ";
      double n = h2->GetBinContent(binx, biny);
      double e = h2->GetBinError(binx, biny);
      N += n;
      e2 += e*e;
      cout << setw(4) << "N = " << n << " +/- " << e << endl;
      values[bin][0] = n; // value
      values[bin][1] = e; // error
      bin++;
    }
  }
  cout << "Total: " << N << " +/- " << TMath::Sqrt(e2) << endl;
  cout << "From histo: " << h2->Integral() << endl;
}

void paperplots(const char * sel)
{
  // make plots for AN and paper
  selection(sel);

  // T/L ratio selection
  tightlooseplots();

  // T/L ratio plots
  get_subtracted_tight_loose_ratio();
  tight_loose_ratioplot();

  // btag control region
  cd(1);
  plot("CR1_m_mumu");
  zoom(0,300);
  legend();
  pprint();
  cd(2);
  plot("CR2_m_mumu");
  zoom(0,300);
  legend();
  pprint();

  // Fakes (signal region) m_smuon
  TH1D * h_smuon = fake_estimate_1d(sel, "btag_m_smuon");
  cd(1);
  rebin(4);
  legend();
  gPad->Print("m_smuon_singlefake.pdf");
  cd(2);
  rebin(4);
  legend();
  gPad->Print("m_smuon_doublefake.pdf");

  // Fakes (signal region) m_gaugino
  TH1D * h_gaugino = fake_estimate_1d(sel, "btag_m_gaugino");
  cd(1);
  rebin(4);
  legend();
  gPad->Print("m_gaugino_singlefake.pdf");
  cd(2);
  rebin(4);
  legend();
  gPad->Print("m_gaugino_doublefake.pdf");

  // Save fakes to a file
  TFile * f = new TFile("fakes.root", "RECREATE");
  h_smuon->SetName("h1_6_jjmm_m");
  h_smuon->Write();
  h_gaugino->SetName("h1_6_jjm_m");
  h_gaugino->Write();
  f->Close();
  delete f;

  // Fakes (control region)
  fake_estimate_1d(sel, "CR6_m_smuon");
  cd(1);
  rebin(5);
  legend();
  gPad->Print("CR_singlefake.pdf");
  cd(2);
  rebin(5);
  legend();
  gPad->Print("CR_doublefake.pdf");
}
