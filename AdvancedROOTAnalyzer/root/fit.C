#include "TMath.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"

#include "plot.h"

/**
 * Crystal Ball function
 * 
 * J. E. Gaiser, Charmonium Spectroscopy from Radiative Decays of the J/Psi and Psi-Prime, 
 * Ph.D. Thesis, SLAC-R-255 (1982) Appendix-F 
 * 
 * Used to model lossy processes in high energy physics. First used by 
 * Crystal Ball experiment, hence its name.
 *
 * It consists of a gaussian core portion and a power-law low-end tail 
 * below a certain threshold. 
 *
 * @param x is the value for which the function is evaluated
 *
 * @param alpha is the number of sigmas away from the Gaussian MPV where the matching 
 * is done. 
 * 
 * @param n is the exponent of the power law. 
 *
 * @param mean is the MPV of the Gaussian core.
 *
 * @param sigma is the standard deviation of the Gaussian core.
 *
 * @return value of crystal ball function for given parameters
 *
 */

double crystal_ball(double x, double alpha, double n, double mean, double sigma)
{
  if (x - mean >= - alpha*sigma) {
    // gaussian core
    double val = (x-mean)/sigma;
    return TMath::Exp(-val*val/2.);
  }
  else {
    // power-law tail
    double B = n/alpha - alpha;
    double A = TMath::Power(n/alpha, n) * TMath::Exp(-alpha*alpha/2.);
    return A*TMath::Power(B-(x-mean)/sigma, -n);
  }
}

void test_crystal_ball()
{
  // plotting range
  const int nMax = 1000;
  const double xlow = -10;
  const double xup = 4;

  // parameters for function
  double n = 1;
  double alpha = 1; 
  double mean = 0.;
  double sigma = 1.;
  const char * title = Form("Crystal ball n = %f, #alpha = %f, #bar{x} = %f, #sigma = %f",
			    n, alpha, mean, sigma);
			    
  TH1D * h1 = new TH1D("h1", title, 1000, xlow, xup);
  for (int i = 0; i < nMax; i++) {
    h1->SetBinContent(i+1, crystal_ball(xlow+(xup-xlow)/1000.*i, alpha, n, mean, sigma));
  }

  n = 2;
  title = Form("Crystal ball n = %f, #alpha = %f, #bar{x} = %f, #sigma = %f",
	       n, alpha, mean, sigma);
			    
  TH1D * h2 = new TH1D("h2", title, 1000, xlow, xup);
  for (int i = 0; i < nMax; i++) {
    h2->SetBinContent(i+1, crystal_ball(xlow+(xup-xlow)/1000.*i, alpha, n, mean, sigma));
  }

  n = 1;
  alpha = 2;
  title = Form("Crystal ball n = %f, #alpha = %f, #bar{x} = %f, #sigma = %f",
	       n, alpha, mean, sigma);
			    
  TH1D * h3 = new TH1D("h3", title, 1000, xlow, xup);
  for (int i = 0; i < nMax; i++) {
    h3->SetBinContent(i+1, crystal_ball(xlow+(xup-xlow)/1000.*i, alpha, n, mean, sigma));
  }

  TCanvas * c1 = new TCanvas;
  c1->cd();
  h1->Draw();
  h2->SetLineColor(kBlue);
  h2->Draw("same");
  h3->SetLineColor(kRed);
  h3->Draw("same");
}

void fitmc(TF1 * f1)
{
  TH1D * back = backgroundHisto("LM1");
  // do a log-likelihood fit (ll), do not plot result (0), and use function range (R)
  back->Fit("f1", "R0", "goff");
  f1->SetLineColor(kBlue);
  f1->SetLineWidth(2.);
  f1->Draw("same");
}
