#include "plot.h"
#include "stat.h"

#include <TMath.h>
#include <TRandom.h>
#include <TLatex.h>
#include <TVector.h>
#include <TMatrix.h>

Double_t chi2_int(Int_t start, Int_t end, Int_t & ndof, Int_t & nlow)
{
  // compute chi2 for the global histograms. Used only internally,
  // use chi2 or chi2low instead
  Double_t ChiSquare = 0;
  ndof = 0;
  nlow = 0;
  Double_t temp, temp1;
  TH1D * hmcsum = gStack[gPadNr][gOrder[gPadNr][0]];
  if (hmcsum == 0) 
    return 0;
  TH1D * hdata = gStack[gPadNr][gMaxProcess-1];
  if (hdata == 0)
    return 0;
  for (Int_t i = start; i <= end; i++) {
    // sum of all MC's ( theorie )
    temp1 = hmcsum->GetBinContent(i);

    // skip empty bins
    if (temp1 < 1E-3)
      continue;

    // theorie - data
    temp = temp1 - hdata->GetBinContent(i);
    ChiSquare += temp * temp / temp1;
    ndof++;
    // count bins with less than five entries
    if (temp1 < 5) {
      nlow++;
    }
  }
  return ChiSquare;
}

void stat(Double_t low, Double_t up)
{
  // count number of signal, background and data MC
  if (FindFirstHisto() < 0) {
    return;
  }

  // total counts
  Double_t  nTotSignal    = 0;                // number of signal
  Double_t  nTotal        = 0;                // number of all mc
  Double_t  nTotExpected  = 0;                // number of expected signal
  Double_t  nData         = 0;                // number of data
  Double_t  nTotErrSig    = 0;
  Double_t  nTotErrBack   = 0;

  // for each mc/signal mc
  Double_t  * nEvents   = new Double_t[gMaxProcess-1]; // number of selected events
  Double_t  * nGene     = new Double_t[gMaxProcess-1]; // number of generated events
  Double_t  * nExpected = new Double_t[gMaxProcess-1]; // number of expected events
  Double_t  * weight    = new Double_t[gMaxProcess-1]; // weight of this MC

  // find bins to start and end with
  Int_t start, end;
  findbins(low, up, start, end);
  Int_t hstart = FindFirstHisto();
  if (hstart < 0)
    return;

  // total number of data events
  if (gHisto[gPadNr][gMaxProcess-1]) {
    nData = gHisto[gPadNr][gMaxProcess-1]->Integral(start, end);
  }
  DEBUG("nData = " << nData);

  // compute numbers for each mc separately
  for (Int_t i = 0; i < gMaxProcess-1; i++) {
    // zero out all values
    nEvents[i] = 0;
    nExpected[i] = 0;
    nGene[i] = 0;
    weight[i] = 0;
    // check if mc exists
    if (gHisto[gPadNr][i] == 0) {
      // mc does not exist
      continue;
    }
    // get integral for this mc
    nEvents[i] = gHisto[gPadNr][i]->Integral(start, end);
    nTotal += nEvents[i];
    // get expected number from global process information
    for (Int_t period = gStart; period <= gEnd; period++) {
      // compute expected events
      nExpected[i] += gLumi[period] * gProcessInfo[period][i].xs;
      DEBUG("gLumi[" << period << "]=" << gLumi[period] << ", xs = " 
	    << gProcessInfo[period][i].xs);
      DEBUG("nExpected = " << nExpected[i]);
      // compute generated events
      nGene[i] += gProcessInfo[period][i].Nev;
      double mcLumi = gProcessInfo[period][i].Nev/gProcessInfo[period][i].xs;
      weight[i] += gLumi[period]/mcLumi * gProcessInfo[period][i].weight;
    }
    weight[i] /= (gEnd-gStart+1);
    // something we do for signal only
    if (i < gMaxSignal) {
      // add signal events
      nTotSignal += nEvents[i];
    }
    nTotExpected += nExpected[i];
  }

  Double_t effi; // efficiency
  Double_t puri; // purity
  Double_t nErr; // error on number of selected events

  printf("*******************************************************************************\n");
  Int_t  nfill = 15 - strlen(gStack[gPadNr][hstart]->GetName());
  printf("%15s       N_events +/-    Delta_N      effi      purity  av. weight\n", 
	 gStack[gPadNr][hstart]->GetName());
  printf("*******************************************************************************\n");
  // inform user for each mc (ordered)
  for (Int_t k = 0; k < gMaxProcess-1; k++) {
    Int_t i = gOrder[gPadNr][k];
    if (nExpected[i] < 1E-9)
      effi = 0;
    else
      effi = nEvents[i] / nExpected[i];
    if (nTotal == 0)
      puri = 0;
    else
      puri = nEvents[i] / nTotal;
    if (effi > 0 && nGene[i] > 0) {
      nErr = nEvents[i] * sqrt((1.-effi)/(effi*nGene[i]));
    }
    else
      nErr = 0;
    // INFO("nExpected= " << nExpected[i] << ", nEvents= " << nEvents[i] << ", nGene= " << nGene[i] << ", effi= " << effi << ", radicand= " << (1-effi)/(effi*nGene[i]));
    if (i < gMaxSignal) {
      nTotErrSig += nErr*nErr;
    }
    else {
      nTotErrBack += nErr*nErr;
    }
    printf("%s:", gProcess[i].fname);
    nfill = 16 - strlen(gProcess[i].fname);
    for (Int_t j = 0; j < nfill; j++)
      printf(" ");
    printf(" %12.3f +/- %10.3f  %8.4f %%  %8.4f %%  %8.3f\n", 
	   nEvents[i], nErr, effi*100, puri*100, weight[i]);
  }

  // summaries
  printf("-------------------------------------------------------------------------------\n");
  if (nTotExpected != 0)
    effi = nTotSignal / nTotExpected;
  else 
    effi = 0;
  if (nTotal != 0)
    puri = nTotSignal / nTotal;
  else
    puri = 0;
  printf("total signal:     %12.3f +/- %10.3f  %8.4f %%  %8.4f %%  %8.6f\n", 
	 nTotSignal, TMath::Sqrt(nTotErrSig), effi*100, puri*100, effi*puri);
  printf("total back:       %12.3f +/- %10.3f\n", 
	 nTotal-nTotSignal, TMath::Sqrt(nTotErrBack));
  printf("===============================================================================\n");
  // total, under & overflow mc
  Double_t uflow  = 0;
  Double_t oflow  = 0;
  for (Int_t k = 0; k < gMaxProcess-1; k++) {
    Int_t i = gOrder[gPadNr][k];
    if (gHisto[gPadNr][i] == 0)
      continue;
    // sum up stacked and unstacked histograms
    uflow      += gHisto[gPadNr][i]->GetBinContent(0);
    Int_t obin  = gHisto[gPadNr][i]->GetNbinsX()+1;
    oflow      += gHisto[gPadNr][i]->GetBinContent(obin);
  }
  printf("total mc:         %12.3f +/- %10.3f (u: %7.3f) (o: %7.3f)   Lumi:\n", 
	 nTotal, TMath::Sqrt(nTotErrSig+nTotErrBack), uflow, oflow);
  // total, under & overflow data
  uflow  = 0;
  oflow  = 0;
  if (gStack[gPadNr][gMaxProcess-1]) {
    uflow      = gStack[gPadNr][gMaxProcess-1]->GetBinContent(0);
    Int_t obin = gStack[gPadNr][gMaxProcess-1]->GetNbinsX()+1;
    oflow      = gStack[gPadNr][gMaxProcess-1]->GetBinContent(obin);
  }
  // compute lumi
  Double_t lumi = 0;
  for (Int_t period = gStart; period <= gEnd; period++) {
    lumi += gLumi[period];
  }
  printf("data:             %12.3f +/- %10.3f (u: %7.3f) (o: %7.3f) %6.2f\n", 
	 nData, TMath::Sqrt(Double_t(nData)), uflow, oflow, lumi);
  // deviation and chi2
  Int_t ndof, nlow;
  Double_t ChiSquare = chi2_int(start, end, ndof, nlow);
  Double_t Deviation = 0;
  if (nTotal != 0)
    Deviation = (nData - nTotal) / TMath::Sqrt(nTotal);
  printf("Deviation:        %9.3f sigma   Chi2/ndof: %8.0f/%3d   Prob: %9.6f %%\n", 
	 Deviation, ChiSquare, ndof, 100 * TMath::Prob(ChiSquare, ndof));

  printf("*******************************************************************************\n");

  delete[] nExpected;
  delete[] nGene;
  delete[] nEvents;
  delete[] weight;
}

/** 
 * Output for each process the efficiency and error on the efficiency
 */
void efficiency()
{
  
}

void chi2(Double_t low = 0, Double_t up = 0)
{
  // compute chi2 for this distribution
  if (FindFirstHisto() < 0)
    return;

  // find bins for given range
  Int_t start, end;
  if (low == up) {
    start = 1;
    end   = gStack[gPadNr][0]->GetNbinsX();
  } 
  else {
    start = gStack[gPadNr][0]->GetXaxis()->FindFixBin(low);
    end   = gStack[gPadNr][0]->GetXaxis()->FindFixBin(up);
  }

  // compute chi2
  Int_t ndof, nlow;
  Double_t ChiSquare = chi2_int(start, end, ndof, nlow);
  printf("Chi2/ndof\t%9.3f/%3d\n", ChiSquare, ndof);
  printf("P-Value\t\t%9.3f%%\n", 100 * TMath::Prob(ChiSquare, ndof));

  // check validation of sample
  if (nlow > 0.1 * ndof)
    printf("WARN: You should use chi2low() to get a more reliable result\n");
  
}

void chi2low(Double_t low = 0, Double_t up = 0, Int_t nevents = 100000)
{
  // compute chi2-distribution for a low number of events...
  if (FindFirstHisto() < 0)
    return;

  // inform user that this may take a while
  printf("This may take a while, please be patient\n");

  // find bins for given range
  Int_t start, end;
  if (low == up) {
    start = 1;
    end   = gStack[gPadNr][0]->GetNbinsX();
  } 
  else {
    start = gStack[gPadNr][0]->GetXaxis()->FindFixBin(low);
    end   = gStack[gPadNr][0]->GetXaxis()->FindFixBin(up);
  }

  Double_t ChiSquare;
  Double_t nTheorie, nDelta;
  // now generate new possible data distribution according to theorie and
  // sum up the chi2
  TH1D * chi2distrib = new TH1D("chi2distrib", "Chi2 distribution (mc)",
				nevents/100, 0, 100);
  for (Int_t i = 0; i < nevents; i++) {
    ChiSquare = 0;
    // loop over each bin
    for (Int_t j = start; j <= end; j++) {
      // get theorie value
      nTheorie = gStack[gPadNr][gOrder[gPadNr][0]]->GetBinContent(j);
      if (nTheorie < 1E-9)
	continue;

      // compute poisson distribution around theorie, take it as data
      nDelta = gRandom->Poisson(nTheorie) - nTheorie;
      ChiSquare += (nDelta*nDelta) / nTheorie;
    }
    chi2distrib->Fill(ChiSquare);
  }
  chi2distrib->SetXTitle(chi2distrib->GetTitle());
  chi2distrib->Scale(1.0/nevents);
  chi2distrib->Draw();
  
  // now compute the P-value
  Int_t ndof, nlow;
  ChiSquare = chi2_int(start, end, ndof, nlow);
  printf("Chi2/ndof\t%9.3f/%3d\n", ChiSquare, ndof);
  Double_t mean = chi2distrib->GetMean();
  printf("<chi2>\t\t%9.3f\n", mean);
  printf("P-value\t\t%9.3f%%\n", 100 * (chi2distrib->Integral((Int_t) ChiSquare, 100)));
}

Double_t loglikelihood(TH1D * hmc, TH1D * hdata, Int_t start, Int_t end)
{
  // compute log likelihood with a poisson probability for
  // the global histograms.
  Double_t pval = hmc->Integral(start, end);
  Double_t n, mu;
  for (Int_t i = start; i <= end; i++) {
    n  = hdata->GetBinContent(i);
    mu = hmc->GetBinContent(i);
    if (mu <= 0) {
      if (n != 0) {
	printf("ERR: Data found but no mc, bin %d\n", i);
	return 0;
      }
      // skip empty bins
      continue;
    }
    pval -= n * TMath::Log(mu) + TMath::Log(TMath::Gamma(n+1));
    if (gLogLevel > 0)
      printf("n = %8.4f, mu = %8.4f, pval = %8.4f\n", n, mu, pval);
  }
  return pval;
}

void loglikdistrib(Int_t ntrials = 10000, Bool_t print = kFALSE)
{
  // compute distribution of log likelihood value
  TH1D * hmc   = gStack[gPadNr][gOrder[gPadNr][0]];
  TH1D * hdata = gStack[gPadNr][gMaxProcess-1];
  Int_t nbins = hmc->GetNbinsX();
  Double_t loglik = loglikelihood(hmc, hdata, 1, nbins);
  TH1D * htest = new TH1D(*hdata);
  TH1D * lldistrib = new TH1D("lldistrib", "log(Likelihood) distribution", 
			      1000, loglik-200, loglik+200);
  setopt(lldistrib);
  for (Int_t n = 0; n < ntrials; n++) {
    // generate poisson around theorie
    for (Int_t i = 1; i <= nbins; i++) {
      htest->SetBinContent(i, gRandom->Poisson(hmc->GetBinContent(i)));
    }
    lldistrib->Fill(loglikelihood(hmc, htest, 1, nbins));
  }
  TCanvas * llcanvas = new TCanvas("llcanvas", "Log(Likelihood) distribution", 
				   40, 40, 800, 600);
  setopt(llcanvas);
  lldistrib->SetFillColor(kYellow);
  lldistrib->Draw();
  lldistrib->GetYaxis()->SetTitle("Anzahl Ereignisse");
  lldistrib->GetXaxis()->SetTitle("-ln L");
  // autozoom
  Int_t lowbin = 1;
  while (lldistrib->GetBinContent(lowbin) == 0)
    lowbin++;
  Int_t highbin = lldistrib->GetNbinsX();
  while (lldistrib->GetBinContent(highbin) == 0)
    highbin--;
  lldistrib->SetAxisRange(lldistrib->GetBinLowEdge(lowbin), 
			  lldistrib->GetBinLowEdge(highbin));
  TH1D * hworse = (TH1D *) lldistrib->Clone();
  for (Int_t nbin = 1; nbin < 501; nbin++) {
    hworse->SetBinContent(nbin, 0);
  }
  hworse->SetFillColor(95);
  hworse->Draw("same");
  Double_t pvalue = lldistrib->Integral(501,1000) / lldistrib->Integral();
  TLatex * tex = new TLatex(0.18, 0.96, Form("-ln L_{obs} = %5.2f", loglik));
  tex->SetNDC();
  tex->SetTextAlign(13);
  tex->Draw();
  tex = new TLatex(0.18, 0.86, Form("CL_{obs} = %.3f", pvalue));
  tex->SetNDC();
  tex->SetTextAlign(13);
  tex->Draw();
  TLine * l = new TLine(loglik, 0, loglik, lldistrib->GetMaximum());
  l->SetLineWidth(3);
  l->SetLineColor(kBlue);
  l->Draw();
  llcanvas->Modified();
  llcanvas->Update();
  if (print)
    llcanvas->Print("lldistrib.pdf");
  cd(gPadNr+1);
}

void searchcut(Bool_t low)
{
  // show efficiency versus purity for a specific cut value
  static TH1D ** heffi = 0;
  static TH1D ** hpuri = 0;
  static TH1D ** hep = 0;

  // create histogram store if needed
  if (heffi == 0)
    heffi = new TH1D * [gMaxSignal];
  if (hpuri == 0)
    hpuri = new TH1D * [gMaxSignal];
  if (hep == 0)
    hep   = new TH1D * [gMaxSignal];
  for (Int_t i = 0; i < gMaxSignal; i++) {
    heffi[i] = 0;
    hpuri[i] = 0;
    hep[i] = 0;
  }

  // look for a histogram
  Int_t start = FindFirstHisto();
  if (start < 0) 
    return;

  DEBUG("first histogram found: " << gProcess[start].fname);

  // norm purity to total number of events
  TH1D * htotal = new TH1D(*gStack[gPadNr][start]);
  DEBUG("Total number of MC events: " << htotal->Integral());

  for (Int_t process = 0; process < gMaxSignal; process++) {
    DEBUG("processsing " << gProcess[process].fname);
    // delete old histos
    if (heffi[process]) {
      delete heffi[process];
      heffi[process] = 0;
    }
    if (hpuri[process]) {
      delete hpuri[process];
      hpuri[process] = 0;
    }
    if (hep[process]) {
      delete hep[process];
      hep[process] = 0;
    }
    // skip non-existing mc's
    TH1D * histo = gHisto[gPadNr][process];
    if (histo == 0)
      continue;
    // create corresponding histograms
    heffi[process] = new TH1D(*histo);
    heffi[process]->Reset();
    hpuri[process] = new TH1D(*histo);
    hpuri[process]->Reset();
    hep[process] = new TH1D(*histo);
    hep[process]->Reset();
    Double_t effi;
    Double_t puri;
    Double_t summc    = 0;
    Double_t sumtotal = 0;
    Double_t nMc      = gHisto[gPadNr][process]->Integral();
    DEBUG("This process has " << nMc << " events");
    // make effi and puri distribution
    if (low) {
      for (Int_t bin = heffi[process]->GetNbinsX(); bin > 1; bin--) {
	summc    += gHisto[gPadNr][process]->GetBinContent(bin);
	sumtotal += htotal->GetBinContent(bin);
	effi = nMc ? summc/nMc : 0;
	puri = sumtotal ? summc/sumtotal : 1;
	heffi[process]->SetBinContent(bin, effi);
	hpuri[process]->SetBinContent(bin, puri);
	hep[process]->SetBinContent(bin, effi*puri);
      }
    }
    else {
      for (Int_t bin = 1; bin <= heffi[process]->GetNbinsX(); bin++) {
	summc    += gHisto[gPadNr][process]->GetBinContent(bin);
	sumtotal += htotal->GetBinContent(bin);
	effi = nMc ? summc/nMc : 0;
	puri = sumtotal ? summc/sumtotal : 1;
	heffi[process]->SetBinContent(bin, effi);
	hpuri[process]->SetBinContent(bin, puri);
	hep[process]->SetBinContent(bin, effi*puri);
      }
    }
    // set line styles
    heffi[process]->SetLineColor(kRed);
    heffi[process]->SetLineStyle(process);
    heffi[process]->SetFillStyle(0);
    hpuri[process]->SetLineColor(kBlue);
    hpuri[process]->SetLineStyle(process);
    hpuri[process]->SetFillStyle(0);
    hep[process]->SetLineColor(kGreen);
    hep[process]->SetLineStyle(process);
    hep[process]->SetFillStyle(0);
  }
  delete htotal;
  // draw histograms
  for (Int_t process = 0; process < gMaxSignal; process++) {
    if (gHisto[gPadNr][process] == 0)
      continue;
    TCanvas * cep = new TCanvas(Form("cep_%d", process), 
				Form("efficiency and purity for %s", 
				     gProcess[process].fname), 
				10*process,10*process, 
				(Int_t) (600*TMath::Sqrt(2.)), 
				600);
    cep->cd();
    heffi[process]->SetMaximum(1.1);
    hpuri[process]->SetMaximum(1.1);
    heffi[process]->Draw();
    hpuri[process]->Draw("same");
    hep[process]->Draw("same");
    // Add a legend
    TLegend * l = new TLegend(0.720, 0.5, 0.92, 0.92, gProcess[process].tname);
    setopt(l);
    l->AddEntry(heffi[process], "Efficiency", "l");
    l->AddEntry(hpuri[process], "Purity", "l");
    l->AddEntry(hep[process], "effi*puri", "l");
    l->Draw();
    // output optimal cut values to user
    printf("Optimal cut value for process %s: %8.4f\n", 
	   gProcess[process].fname,
	   hep[process]->GetBinLowEdge(hep[process]->GetMaximumBin()));
  }
  cd(gPadNr+1);
}

void uppercut()
{
  searchcut(kFALSE);
}

void lowercut()
{
  searchcut(kTRUE);
}

void integral(Bool_t up = kTRUE, Bool_t draw=kTRUE)
{
  // make integral distribution plot
  // look for a histogram
  Int_t start = FindFirstHisto();
  if (start < 0)
    return;
  if (gStack[gPadNr][gMaxProcess-1] == 0) {
    printf("ERR: no data histo\n");
    return;
  }

  TH1D * mc = gStack[gPadNr][start];
  TH1D * hintmc = (TH1D *) mc->Clone();
  hintmc->Reset();
  hintmc->SetName("hintmc");
  TH1D * data = gStack[gPadNr][gMaxProcess-1];
  TH1D * hintdata = (TH1D *) data->Clone();
  hintdata->Reset();
  hintdata->SetName("hintdata");
  Double_t summc = 0;
  Double_t sumdata = 0;
  Double_t maximum = 0;
  // make integral
  if (up) {
    // loop over all bins including over- and underflow
    for (Int_t bin = 0; bin <= hintmc->GetNbinsX()+1; bin++) {
      summc   += mc->GetBinContent(bin);
      hintmc->SetBinContent(bin, summc);
      sumdata += data->GetBinContent(bin);
      hintdata->SetBinContent(bin, sumdata);
    }
    maximum = TMath::Max(hintmc->GetBinContent(hintmc->GetNbinsX()+1),
			 hintdata->GetBinContent(hintdata->GetNbinsX()+1));
    maximum += TMath::Sqrt(maximum);
  }
  else {
    // loop over all bins including over- and underflow
    for (Int_t bin = hintmc->GetNbinsX()+1; bin >= 0; bin--) {
      summc   += mc->GetBinContent(bin);
      hintmc->SetBinContent(bin, summc);
      sumdata += data->GetBinContent(bin);
      hintdata->SetBinContent(bin, sumdata);
    }
    maximum = TMath::Max(hintmc->GetBinContent(0),
			 hintdata->GetBinContent(0));
    maximum += TMath::Sqrt(maximum);
  }
  hintmc->SetMaximum(maximum);
  hintdata->SetMaximum(maximum);
  // draw
  if (draw) {
    TCanvas * integ = new TCanvas("integ", "integral distribution", 
				  10,10, (Int_t) (600*TMath::Sqrt(2.)), 600);
    setopt(integ);
    hintmc->Draw();
    hintdata->Draw("epsame");
    integ->Modified();
    integ->Update();
  }
  cd(gPadNr+1);
}

Double_t q(Double_t val)
{
  return val*val;
}

void ComputeCorrelatedError(Int_t npoints, // number of measurements
			    Double_t val[], // the individual measurements
			    Double_t sigma_stat[], // statistical error
			    Double_t sigma_sys_uncorr[], // uncorrelated systematic error
			    Double_t sigma_sys_corr[], // correlated systematic error
			    Double_t & mean, // returned: the mean value
			    Double_t & sigma) // returned: the sigma
{
  // build covariance matrix
  TMatrix V(npoints, npoints);
  Double_t temp;
  for (Int_t i = 0; i < npoints; i++) {
    // diagonal
    V(i,i) = q(sigma_stat[i])+q(sigma_sys_uncorr[i])+q(sigma_sys_corr[i]);
    // off-diagonal
    for (Int_t j = i+1; j < npoints; j++) {
      temp = q(sigma_sys_corr[i]);
      V(i,j) = temp;
      V(j,i) = temp;
    }
  }

  if (gLogLevel > 10)
    V.Print();

  // invert covariance matrix
  TMatrix IV(TMatrix::kInverted, V);

  if (gLogLevel > 10)
    IV.Print();

  // compute weight vector
  TVector w(npoints);
  Double_t ivsum = 0;  // sum of all matrix elements
  Double_t rsum;       // sum of a row
  for (Int_t i = 0; i < npoints; i++) {
    rsum = 0;
    for (Int_t j = 0; j < npoints; j++) {
       rsum += IV(i,j);
    }
    w(i) = rsum;
    ivsum += rsum;
    if (gLogLevel > 10)
      printf("rsum = %8.4f,  ivsum = %8.4f\n", rsum, ivsum);
  }

  // normalize elements of weight vector
  Double_t wsum = 0;
  for (Int_t i = 0; i < npoints; i++) {
    w(i) /= ivsum;
    wsum += w(i);
  }

  if (gLogLevel > 10) {
    w.Print();
    printf("wsum = %8.4f\n", wsum);
  }

  // compute averaged value
  mean = 0;
  for (Int_t i = 0; i < npoints; i++) {
    mean += w(i) * val[i];
  }
  if (gLogLevel > 2)
    printf("mean: %8.4f\n", mean);

  // compute variance
  Double_t variance = 0;
  for (Int_t i = 0; i < npoints; i++) {
    for (Int_t j = 0; j < npoints; j++) {
      variance += w(i) * V(i,j) * w(j);
    }
  }
  sigma = TMath::Sqrt(variance);
  if (gLogLevel > 2)
    printf("variance: %8.4f, sigma = %8.4f\n", variance, sigma);
}

