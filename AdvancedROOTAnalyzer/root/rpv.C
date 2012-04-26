#ifndef __CINT__

#include <iostream>

#include "plot.h"
#include "TMath.h"
#include "TROOT.h"
#include "TFile.h"

#endif // __CINT__

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

Bool_t read_xsfile(const char * name, int nMax, int & n, double m[], double xs_lo[], double xs_nlo[])
{
  FILE * xsfile = fopen(fname, "r");
  if (xsfile == 0) {
    ERROR("could not open file " << fname);
    return kFALSE;
  }
  char buffer[256];
  int p1, p2;
  double mS, Mlr, lambda;
  double LO, NLO, k;
  double LO_c, NLO_c, k_c;
  double LO_tot, NLO_tot;
  while (fgets(buffer, 256, xsfile)) {
    if (n >= nMax) { 
      WARNING("File contains more entries than array - enlarge array");
      break;
    }
    if (buffer[0] == '#')
      continue;
    if (sscanf(buffer, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
	       & p1, & p2,
	       & mS, & Mlr, & lambda,
	       & LO, & NLO, & k,
	       & LO_c, & NLO_c, & k_c,
	       & LO_tot, & NLO_tot) != 13) {
      WARNING("Error in format in this line: " << buffer);
      continue;
    }
    xs_lo[n] = LO_tot;
    xs_nlo[n] = NLO_tot;
    m[n] = mS;
    n++;
  }
}

void xsection(const char * fname = "scan_14TeV_cs.txt")
{
  int nMax = 1000;
  double m[1000];
  double xs_lo[1000];
  double xs_nlo[1000];
  int n = 0;
  double xs_min = 1E9;
  double xs_max = 0;
  double m_min = 1E9;
  double m_max = 0;
  char quarks[6] = { 'd', 'u', 's', 'c', 'b', 't' };
  ERROR("this cannot work / fix code");
  while (false) {
    xs_max = TMath::Max(TMath::Max(xs_max, LO_tot), NLO_tot);
    xs_min = TMath::Min(TMath::Min(xs_min, LO_tot), NLO_tot);
    m_min = TMath::Min(m_min, mS);
    m_max = TMath::Max(m_max, mS);
  }
  INFO("Found " << n << " data points in file " << fname);
  MakeCanvas(1,1);
  gPad->SetLogy();
  TH2F * hframe = new TH2F("hframe", "frame", 1, m_min, m_max, 1, xs_min, xs_max);
  setopt(hframe);
  hframe->SetXTitle("m(#bf{#tilde{#mu}}) [GeV]");
  hframe->SetYTitle(Form("#sigma(%c %c #rightarrow #bf{#tilde{#mu}}) [pb]", quarks[TMath::Abs(p1)-1], quarks[TMath::Abs(p2)-1]));
  hframe->Draw();
  
  TGraph * gr_lo = new TGraph(n, m, xs_lo);
  // setopt(gr_lo);
  gr_lo->SetLineColor(kRed);
  gr_lo->SetLineStyle(kDashed);
  gr_lo->Draw("plsame");

  TGraph * gr_nlo = new TGraph(n, m, xs_nlo);
  // setopt(gr_nlo);
  gr_nlo->SetLineColor(kBlue);
  gr_nlo->SetLineStyle(kSolid);
  gr_nlo->Draw("plsame");

  TLegend * leg = new TLegend(0.6, 0.76, 0.95, 0.91, Form("#sqrt{s} = 14 TeV, #lambda'_{211} = %.3f", lambda));
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

  // plot 1
  cd(1);
  plot("Sig_nMuon");
  zoom(0, 6);
  legend();
  cd(2);
  plot("Sig_nNeutrino");
  zoom(0, 6);
  print("Sig_nMuNeutrino.pdf");

  // plot 2
  cd(1);
  plot("Sig_SleptonMass");
  zoom(120,220);
  legend();
  cd(2);
  plot("Sig_SleptonCharge");
  print("Sig_SleptonMassAndCharge.pdf");

  // plot 3
  cd(1);
  plot("Sig_EtaMu");
  rebin(2);
  legend();
  cd(2);
  plot("Sig_PhiMu");
  rebin(4);
  print("Sig_MuEtaPhi.pdf");

  // plot 4
  cd(1);
  plot("Sig_ptMu1");
  zoom(0,80);
  legend(0.01, -1);
  cd(2); 
  plot("Sig_ptMu2");
  zoom(0,80);
  print("Sig_MuPt12.pdf");

  // plot 5
  cd(1);
  plot("Sig_EtaQuark");
  legend();
  cd(2);
  plot("Sig_PhiQuark");
  rebin(4);
  print("Sig_QuarkEtaPhi.pdf");

  // plot 6
  cd(1);
  plot("Sig_ptQuark1");
  zoom(0,80);
  legend();
  cd(2);
  plot("Sig_ptQuark2");
  zoom(0,80);
  print("Sig_QuarkPt.pdf");

  // plot 7
  cd(1);
  plot("Sig_Mass3");
  zoom(50, 150);
  legend();
  cd(2);
  plot("Sig_MuIsoR");
  print("Sig_MassDeltaR.pdf");

  // plot 8
  cd(1);
  plot("vtx_n");
  legend();
  cd(2);
  plot("vtx_x");
  logy();
  min(0.1);
  print("Vertex_1.pdf");

  // plot 9
  cd(1);
  plot("vtx_y");
  logy();
  min(0.1);
  legend();
  cd(2);
  plot("vtx_z");
  liny();
  print("Vertex_2.pdf");

  // plot 9
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
