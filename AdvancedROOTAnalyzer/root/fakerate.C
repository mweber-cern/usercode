#include <iostream>
#include <typeinfo>

using namespace std;

#include "TROOT.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TFile.h"
#include "THashList.h"
#include "TString.h"
#include "TMath.h"
#include "TTree.h"
#include "TKey.h"

#include "plot.h"
#include "stat.h"
#include "xsection.h"

// Project a 3D histogram in 2D by integrating over z bins
// do not project overflow / underflow
TH2D * project3d_xy(TH3D * histo3)
{
  WARNING("project3d_xy(): untested - use at your own risk!");
  if (histo3 == 0) 
    return 0;
  // create 2D histogram with same binning
  TH2D * h2 = new TH2D(Form("%s_yx", histo3->GetName()), 
		       histo3->GetTitle(), 
		       histo3->GetNbinsX(), 
		       histo3->GetXaxis()->GetBinLowEdge(1),
		       histo3->GetXaxis()->GetBinLowEdge(histo3->GetNbinsX()+1),
		       histo3->GetNbinsY(), 
		       histo3->GetYaxis()->GetBinLowEdge(1),
		       histo3->GetYaxis()->GetBinLowEdge(histo3->GetNbinsY()+1));

  // copy histogram contents
  for (Int_t i = 1;  i < histo3->GetNbinsX()+1; i++) {
    for (Int_t j = 1; j < histo3->GetNbinsY()+1; j++) {
      double sum = 0;
      for (Int_t k = 1; k < histo3->GetNbinsZ()+1; k++) {
	sum += histo3->GetBinContent(i, j, k);
      }
      h2->SetBinContent(i, j, sum);
    }
  }  
  return h2;
}

const char * get_bin_edges(const TH1 * histo, Int_t i, Int_t j, Int_t k)
{
  static char buffer[255];
  double x = histo->GetXaxis()->GetBinLowEdge(i);
  double y = histo->GetYaxis()->GetBinLowEdge(j);
  double z = histo->GetZaxis()->GetBinLowEdge(k);
  sprintf(buffer, "%4.2f,%4.2f,%4.2f", x, y, z);
  return buffer;
}

const char * get_bin_edges(const TH1 * histo, Int_t i, Int_t j)
{
  static char buffer[255];
  double x = histo->GetXaxis()->GetBinLowEdge(i);
  double y = histo->GetYaxis()->GetBinLowEdge(j);
  sprintf(buffer, "%4.2f,%4.2f", x, y);
  return buffer;
}

void check_3d_histo(const TH3D * hTight3, const TH3D * hLoose3, bool check_negative)
{
  for (Int_t i = 1;  i < hTight3->GetNbinsX()+1; i++) {
    for (Int_t j = 1; j < hTight3->GetNbinsY()+1; j++) {
      for (Int_t k = 1; k < hTight3->GetNbinsZ()+1; k++) {
	  Double_t nLoose = hLoose3->GetBinContent(i, j, k);
	  Double_t nTight = hTight3->GetBinContent(i, j, k);
	// Double_t nRatio = hRatio2->GetBinContent(i, j, k);

	if (nLoose < nTight) {
	  cerr << "ERR: nLoose < nTight" << endl;
	  cout << "INFO: nLoose(" << i << "," << j << "," << k << ") =" << nLoose << endl;
	  cout << "INFO: nTight(" << i << "," << j << "," << k << ") =" << nTight << endl;
	}
	if (check_negative && nLoose < 0) {
	  cerr << "ERR: nLoose < 0" << endl;
	  cout << "INFO: nLoose(" << i << "," << j << "," << k << ") =" << nLoose << endl;
	  cout << "INFO: nTight(" << i << "," << j << "," << k << ") =" << nTight << endl;
	}
        // if (nLoose == 0) {
        // 	WARNING("nLoose(" << i << "," << j << "," << k << ") = 0");
        // 	WARNING("nRatio = " << nRatio);
        // }
        // if (nLoose != 0 && nRatio != nTight/nLoose) {
        // 	cerr << "ERR: ROOT calculation went wrong" << endl;
        // }
        // if (nRatio < 0 || nRatio >= 1.) {
        // 	cerr << "ERR: Found unreasonable values for T/L ratio" << endl;
        // }
        // if (nRatio != 0) {
        // 	// cout << "INFO: T/L (" << get_bin_edges(hRatio2, i, j) << ") " 
        // 	//      << nRatio << " +- " << error << endl;
        // }
      }
    }
  }  
}

void replace_with_average(TH2D * histo, int i, int j)
{
  int min_i = TMath::Max(1, i-1);
  int max_i = TMath::Min(histo->GetNbinsX(), i+1);
  int min_j = TMath::Max(1, j-1);
  int max_j = TMath::Min(histo->GetNbinsY(), j+1);
  INFO("replacing from x bins = " << min_i << ".." << max_i);
  INFO("replacing from y bins = " << min_j << ".." << max_j);
  double average = 0.;
  int count = 0;
  for (int ni = min_i ; ni <= max_i; ni++) {
    for (int nj = min_j; nj <= max_j; nj++) {
      // do not take into account wrong bin
      if (ni == i || nj == j) 
	continue;
      // build average
      average += histo->GetBinContent(ni, nj);
      count +=1;
    }
  }
  average /= count;
  INFO("Fixing bin " << get_bin_edges(histo, i, j) 
       << " old value = " << histo->GetBinContent(i,j)
       << " to average = " << average);
  histo->SetBinContent(i, j, average);
}

void fix_2d_histo(TH2D * hTight2, TH2D * hLoose2)
{
  for (Int_t i = 1;  i < hTight2->GetNbinsX()+1; i++) {
    for (Int_t j = 1; j < hTight2->GetNbinsY()+1; j++) {
      Double_t nLoose = hLoose2->GetBinContent(i, j);
      Double_t nTight = hTight2->GetBinContent(i, j);

      if (nLoose < 0) {
	ERROR("nLoose < 0");
	INFO("nLoose(" << get_bin_edges(hLoose2, i, j) << " ) =" << nLoose);
	replace_with_average(hLoose2, i, j);
      }
      if (nLoose < nTight || nTight < 0) {
	ERROR("nLoose < nTight || nTight < 0");
	INFO("nLoose(" << get_bin_edges(hLoose2, i, j) << " ) =" << nLoose);
	INFO("nTight(" << get_bin_edges(hTight2, i, j) << " ) =" << nTight);
	replace_with_average(hTight2, i, j);
      }
    }
  }  
}

TH1D * get_subtracted_tight_loose_ratio(bool save=true, bool draw=true)
{
  if (gStart != gEnd) {
    ERROR("tight_loose_ratio() can only be called for one period!");
    return 0;
  }
  if (draw) {
    MakeCanvas(2,2);
  }
  // "plot" histograms -> put into memory
  int oldpad = gPadNr;
  cd(1);
  INFO("Plotting tight muons");
  plot3("TightMuons");
  cd(2);
  INFO("Plotting loose muons");
  plot3("LooseMuons");
  TH3D * hTight3_data = gHisto3[0][gMaxProcess-1];
  if (hTight3_data == 0) {
    ERROR("hTight3_data == 0");
    return 0;
  }
  TH3D * hLoose3_data = gHisto3[1][gMaxProcess-1];
  if (hLoose3_data == 0) {
    ERROR("hLoose3_data == 0");
    return 0;
  }
  cd(oldpad+1);
  double nTightBefore = hTight3_data->Integral();
  double nLooseBefore = hLoose3_data->Integral();
  double nTightAfter  = nTightBefore;
  double nLooseAfter  = nLooseBefore;
  INFO("Statistics before subtraction");
  INFO("Tight muons: " << nTightBefore << ", loose muons: " << nLooseBefore);
  // subtract background contribution from data
  for (Int_t j = 0; j < gMaxProcess-1; j++) {
    TString proc(gProcess[j].fname);
    if (proc.Contains("qcd") || proc.Contains("signal") || proc.Contains("background")) {
      INFO("Skipping " << gProcess[j].fname);
      continue;
    }
    TH3D * hTight3_back = gHisto3[0][j];
    if (hTight3_back == 0) {
      ERROR("gHisto[0]["<<j<<"] == 0");
      continue;
    }
    TH3D * hLoose3_back = gHisto3[1][j];
    if (hLoose3_back == 0) {
      ERROR("gHisto[1]["<<j<<"] == 0");
      continue;
    }
    double nTight = hTight3_back->Integral();
    double nLoose = hLoose3_back->Integral();
    if (nTight != 0 || nLoose != 0) {
      INFO(gProcess[j].fname << ": subtracting " << nTight << " tight and " 
	   << nLoose << " loose events");
      INFO("Testing background histograms");
      check_3d_histo(hTight3_back, hLoose3_back, true);
      hTight3_data->Add(hTight3_back, -1.);
      nTightAfter -= nTight;
      hLoose3_data->Add(hLoose3_back, -1.);
      nLooseAfter -= nLoose;
      INFO("Testing subtracted data histograms");
      check_3d_histo(hTight3_data, hLoose3_data, false);
    }
    else {
      ERROR("Histograms not found");
    }
  }
  // output some statistics
  INFO("After subtraction:");
  INFO("Count: Tight muons: " << nTightAfter << ", loose muons: " << nLooseAfter);
  INFO("Histo: Tight muons: " << hTight3_data->Integral() 
       << ", loose muons: " << hLoose3_data->Integral());

  // compute 3D ratio
  TH3D * hRatio3 = new TH3D(*hTight3_data);
  hRatio3->SetDirectory(0);
  hRatio3->SetName("hRatio3");
  hRatio3->SetTitle("Tight/Loose ratio (T/L) (3D)");

  // divide with binomial errors
  hRatio3->Divide(hTight3_data, hLoose3_data, 1., 1., "B");
  // some sanity checks
  for (Int_t i = 1;  i < hRatio3->GetNbinsX()+1; i++) {
    for (Int_t j = 1; j < hRatio3->GetNbinsY()+1; j++) {
      for (Int_t k = 1; k < hRatio3->GetNbinsZ()+1; k++) {
	Double_t nLoose = hLoose3_data->GetBinContent(i, j, k);
	Double_t nTight = hTight3_data->GetBinContent(i, j, k);
	Double_t nRatio = hRatio3->GetBinContent(i, j, k);
	Double_t error  = hRatio3->GetBinError(i, j, k);

	if (nLoose < nTight || nLoose < 0) {
	  cerr << "ERR: nLoose is wrong in your histogram" << endl;
	  cout << "INFO: nLoose(" << get_bin_edges(hRatio3, i, j, k) << ") " << nLoose << endl;
	  cout << "INFO: nTight(" << get_bin_edges(hRatio3, i, j, k) << ") " << nTight << endl;
	}
	if (nLoose != 0 && nRatio != nTight/nLoose) {
	  cerr << "ERR: ROOT calculation went wrong" << endl;
	}
	if (nRatio < 0 || nRatio >= 1.) {
	  cerr << "ERR: Found unreasonable values for T/L ratio" << endl;
	}
	if (nRatio != 0) {
	  cout << "INFO: T/L (" << get_bin_edges(hRatio3, i, j, k) << ") " 
	       << nRatio << " +- " << error << endl;
	}
      }
    }
  }

  // now work in 2D-Projection (NUF = do not project underflow, NOF = ... overflow)
  TH2D * hTight2 = (TH2D *) hTight3_data->Project3D("yxNUFNOF");
  TH2D * hLoose2 = (TH2D *) hLoose3_data->Project3D("yxNUFNOF");
  if (fabs(hTight2->Integral() - hTight3_data->Integral())/hTight3_data->Integral() > 1e-10) {
    ERROR("Project3D error tight");
    INFO("2D integral: " << hTight2->Integral() 
	 << ", 3d integral = " << hTight3_data->Integral());
    INFO("integral difference: " << hTight3_data->Integral() - hTight2->Integral());
  }
  if (fabs(hLoose2->Integral() - hLoose3_data->Integral())/hLoose3_data->Integral() > 1e-10) {
    ERROR("Project3D error loose");
    INFO("2D integral: " << hLoose2->Integral() 
	 << ", 3d integral = " << hLoose3_data->Integral());
    INFO("integral difference: " << hLoose2->Integral() - hLoose3_data->Integral());
  }
  // fix 2d histogram for insane values
  fix_2d_histo(hTight2, hLoose2);

  // compute 2D ratio
  TH2D * hRatio2 = new TH2D(*hTight2);
  hRatio2->SetDirectory(0);
  hRatio2->SetName("hRatio2");
  hRatio2->SetTitle("Tight/Loose ratio (T/L) (2D)");

  if (draw) {
    cd(3); 
    setopt(hLoose2);
    hLoose2->SetXTitle("muon p_{T} [GeV]");
    hLoose2->SetYTitle("muon #eta");
    hLoose2->SetZTitle("number of loose muons");
    hLoose2->SetTitleOffset(1.5, "X");
    hLoose2->SetTitleOffset(1.5, "Y");
    hLoose2->SetTitleOffset(2.5, "Z");
    hLoose2->Draw("lego2");
    cd(4); 
    setopt(hTight2);
    hTight2->SetXTitle("muon p_{T} [GeV]");
    hTight2->SetYTitle("muon #eta");
    hTight2->SetZTitle("number of tight muons");
    hTight2->SetTitleOffset(1.5, "X");
    hTight2->SetTitleOffset(1.5, "Y");
    hTight2->SetTitleOffset(2.5, "Z");
    hTight2->Draw("lego2");
  }

  // divide with binomial errors
  hRatio2->Divide(hTight2, hLoose2, 1., 1., "B");
  // some sanity checks
  for (Int_t i = 1;  i < hRatio2->GetNbinsX()+1; i++) {
    for (Int_t j = 1; j < hRatio2->GetNbinsY()+1; j++) {
      Double_t nLoose = hLoose2->GetBinContent(i, j);
      Double_t nTight = hTight2->GetBinContent(i, j);
      Double_t nRatio = hRatio2->GetBinContent(i, j);
      Double_t error  = hRatio2->GetBinError(i, j);

      if (nLoose < nTight || nLoose < 0 || nTight < 0) {
	ERROR("Numbers are wrong in your histogram");
	INFO("nLoose(" << get_bin_edges(hRatio2, i, j) << ") " << nLoose);
	INFO("nTight(" << get_bin_edges(hRatio2, i, j) << ") " << nTight);
	if (nLoose < 0) {
	  // hLoose2->SetBinContent(i, j, 0.);
	  // hTight2->SetBinContent(i, j, TMath::Max(nTight, 0.));
	}
	if (nTight < 0) {
	  // hTight2->SetBinContent(i, j, 0.);
	}
      }
      if (nLoose != 0 && nRatio != nTight/nLoose) {
	ERROR("ROOT calculation went wrong");
      }
      if (nRatio < 0 || nRatio >= 1.) {
	ERROR("Found unreasonable values for T/L ratio");
	INFO("T/L (" << get_bin_edges(hRatio2, i, j) << ") " 
	     << nRatio << " +- " << error);
	replace_with_average(hRatio2, i, j);
      }
      if (nRatio != 0) {
	DEBUG("T/L (" << get_bin_edges(hRatio2, i, j) << ") " 
	      << nRatio << " +- " << error);
      }
    }
  }
  if (draw) {
    cd(1);
    setopt(hRatio2);
    hRatio2->SetXTitle("muon p_{T} [GeV]");
    hRatio2->SetYTitle("muon #eta");
    hRatio2->SetZTitle("T/L ratio");
    hRatio2->SetTitleOffset(1.5, "X");
    hRatio2->SetTitleOffset(1.5, "Y");
    hRatio2->SetTitleOffset(2.5, "Z");
    hRatio2->Draw("lego2");
  }

  // compute 1D ratio
  TH1D * hTight1 = (TH1D *) hTight3_data->Project3D("x");
  TH1D * hLoose1 = (TH1D *) hLoose3_data->Project3D("x");
  TH1D * hRatio1 = new TH1D(*hTight1);
  hRatio1->SetDirectory(0);
  hRatio1->SetName("hRatio1");
  hRatio1->SetTitle("Tight/Loose ratio (T/L)");

  // divide with binomial errors
  hRatio1->Divide(hTight1, hLoose1, 1., 1., "B");

  if (save) {
    // overwrite existing file
    const char * fakeRateFileName = Form("../config/FakeRate_%s.root", gSubDir);
    TFile * fakeRateFile = new TFile(fakeRateFileName, "RECREATE");
    if (fakeRateFile == 0 || !fakeRateFile->IsOpen()) {
      ERROR("Could not open fake rate file " << fakeRateFileName);
      return 0;
    }
    hRatio1->Write();
    hRatio2->Write();
    hRatio3->Write();
    fakeRateFile->Close();
    delete fakeRateFile;
  }

  if (draw) {
    cd(2);
    hRatio1->Draw();
  }
  return hRatio1;
}

void check_division(TH1D * hTight1, TH1D * hLoose1, TH1D * hRatio1)
{
  if (hTight1->GetNbinsX() != hLoose1->GetNbinsX() ||
      hTight1->GetNbinsX() != hRatio1->GetNbinsX()) {
    ERROR("histograms have different binning, exiting");
    return;
  }
  for (Int_t i = 1; i < hTight1->GetNbinsX()+1; i++) {
    double t = hTight1->GetBinContent(i);
    double l = hLoose1->GetBinContent(i);
    double r = hRatio1->GetBinContent(i);

    if (l == 0) {
      WARNING("histogram bin " << i << " contains zero events");
      INFO("r = " << r);
    }
    else {
      if (r != t/l)
	ERROR("check_division(): ROOT calculation went wrong");
    }
  }
}

void check_1d_ratio_histogram(TH1D * histo)
{
  if (histo == 0) {
    ERROR("check_1d_histogram(): zero pointer for histogram");
    return;
  }
  for (Int_t i = 1; i < histo->GetNbinsX()+1; i++) {
    double n = histo->GetBinContent(i);
    if (n < 0) {
      ERROR("histogram has < 0 as entry");
    }
    if (n > 1) {
      ERROR("histogram has > 1 as entry");
    }
    if (n == 0) {
      WARNING("histogram has == 0 as entry");
    }
    if (n == 1) {
      WARNING("histogram has == 1 as entry");
    }
  }
}

TH1D * get_1d_ratio(TH3D * hTight3, TH3D * hLoose3) 
{
  // project on x-axis
  TH1D * hTight1 = (TH1D *) hTight3->Project3D("x");
  TH1D * hLoose1 = (TH1D *) hLoose3->Project3D("x");
  TH1D * hRatio1 = new TH1D(*hTight1);
  hRatio1->SetDirectory(0);
  hRatio1->SetName("hRatio1");
  hRatio1->SetTitle("Tight/Loose ratio (T/L)");

  // divide with binomial errors
  hRatio1->Divide(hTight1, hLoose1, 1., 1., "B");

  // check division
  check_division(hTight1, hLoose1, hRatio1);
  return hRatio1;
}

void tight_loose_ratioplot()
{
  DEBUG("tight_loose_ratioplot() start");
  MakeCanvas();

  DEBUG("reading histograms");
  // "plot" histograms, i.e. read into memory
  cd(1);
  top();
  plot3("TightMuons");
  cd(2);
  top();
  plot3("LooseMuons");
  cd(1);

  TLegend * leg = new TLegend(0.53, 0.52, 0.77, 0.84);
  setopt(leg);
  Bool_t first = kTRUE;

  // get summed histograms
  vector<Int_t> entries;
  TH3D * hTightSum = 0;
  TH3D * hLooseSum = 0;
  for (Int_t i = 0; i < gMaxProcess; i++) {
    Int_t process = gOrder[0][i];
    // not existing
    TH3D * hTight3 = gHisto3[0][process];
    TH3D * hLoose3 = gHisto3[1][process];
    if (hTight3 == 0 || hLoose3 == 0) {
      ERROR("could not get histogram # " << i << " from memory");
      continue;
    }
    // if joined, it was already considered in previous iteration
    if (gProcess[process].join)
      continue;
    if (strncmp(gProcess[process].fname, "qcd", 3) &&
	strncmp(gProcess[process].fname, "dyll", 4) &&
	strncmp(gProcess[process].fname, "ttjets", 6) &&
	strncmp(gProcess[process].fname, "wjetstolnu", 6) &&
	strncmp(gProcess[process].fname, "data", 3))
      continue;
    DEBUG("Creating histos for process " << gProcess[process].fname);
    hTightSum = new TH3D(*hTight3);
    hLooseSum = new TH3D(*hLoose3);
    // check if histograms are joined -> we need to add statistics
    for (Int_t j = i+1; j < gMaxProcess; j++) {
      Int_t proc = gOrder[0][j];
      if (gOrder[1][j] != proc) {
	ERROR("wrong order - need to restart!");
	return;
      }
      // only add joined histograms
      if (!gProcess[proc].join)
	break;
      // not existing
      if (gHisto3[0][proc] == 0 || gHisto3[1][proc] == 0)
	continue;
      DEBUG("Adding histos for process " << gProcess[proc].fname);
      hTightSum->Add(gHisto3[0][proc], 1.);
      hLooseSum->Add(gHisto3[1][proc], 1.);
    }
    DEBUG("Getting ratio");
    TH1D * histo = get_1d_ratio(hTightSum, hLooseSum);
    if (histo == 0) {
      ERROR("Division failed");
      return;
    }
    check_1d_ratio_histogram(histo);
    setopt(histo);
    histo->SetLineColor(gProcess[process].lcolor);
    histo->SetLineStyle(gProcess[process].lstyle);
    DEBUG("Adding to legend");
    if (process != gMaxProcess-1) {
      leg->AddEntry(histo, gProcess[process].tname, "l");
    }
    else {
      leg->AddEntry(histo, gProcess[process].tname);
    }
    DEBUG("Draw");
    if (first) {
      histo->SetMaximum(1.);
      histo->SetMinimum(0.);
      histo->SetXTitle("p_{T}(#mu) [GeV]");
      histo->SetYTitle("T/L ratio");
      histo->Draw("ehisto");
      first = kFALSE;
    }
    else {
      if (process == gMaxProcess-1) {
	histo->SetMarkerStyle(gProcess[gMaxProcess-1].marker);
	histo->Draw("epsame");
      }
      else {
	histo->Draw("ehistosame");
      }
    }
  }
  // add subtracted data histogram
  TH1D * hdata_subtracted = get_subtracted_tight_loose_ratio(false, false);
  if (hdata_subtracted == 0)
    return;
  hdata_subtracted->SetMarkerColor(kBlue);
  hdata_subtracted->SetMarkerStyle(8);
  hdata_subtracted->Draw("epsame");
  leg->AddEntry(hdata_subtracted, "data subtr.", "ep");
  leg->Draw();
}

void tightlooseplots(int start = 1, int end = 5)
{
  for (int i = start; i <= end; i++) {
    switch (i) {
      case 1:
	// plot 1
	cd(1);
	plot("nTL_met");
	logy();
	min(0.1);
	arrow(50.);
	legend(1e-3);
	cd(2);
	plot("nTL_ht");
	max(1e6);
	min(0.1);
	logy();
	arrow(50.);
	legend(1e-3);
	print("tightloose_1.pdf");
	break;
      case 2:
	// plot 2
	cd(1);
	plot("nTL_mupt");
	logy();
	min(0.1);
	max(1e6);
	arrow(15.);
	legend(1e-3);

	cd(2);
	plot("nTL_jetpt");
	logy();
	min(0.1);
	max(1e6);
	zoom(0, 700);
	arrow(40.);
	legend(1e-3);
	print("tightloose_2.pdf");
	break;
      case 3:
	// plot 3
	cd(1);
	liny();
	plot("nTL_jetdphi");
	min(0.1);
	arrow(1.);
	legend(1e-3, 0);

	cd(2);
	plot("nTL_mt");
	zoom(0, 100);
	liny();
	arrow(40.);
	legend(1e-3);
	print("tightloose_3.pdf");  
	break;
      case 4:
	// plot 4
	cd(1);
	liny();
	plot("nTL_zmass");
	arrow(71.);
	arrow(111.);
	zoom(0, 200);
	min(0.1);
	legend(1e-3);
	
	cd(2);
	plot("nTL_nloose");
	min(0.1);
	logy();
	arrow(0.5);
	arrow(1.5);
	legend(1e-3);
	print("tightloose_4.pdf");  
	break;

      case 5:
	// plot 5
	plot("TL_metdphi");
	liny();
	min(0.1);
	legend(1e-3);
	print("tightloose_5.pdf");
	break;
    }
  }
}

void doublefakeplots(const char * sele = "doublefake")
{
  selection(sele); 
  cd(1);
  plot("DF_eta"); 
  cd(2); 
  plot("DF_pt");
  print("doublefake_1.pdf");

  
  cd(1);
  plot("cutflow");
  liny();
  max(70);
  cd(2);
  plot("m_mumu");
  rebin(5);
  stat();
  print("doublefake_2.pdf");

  cd(1); 
  plot("fFakeRate", "1.");
  cd(2);
  TF1 * f1 = new TF1("f1", "x/(1.-x)", 0, 0.5);
  f1->Draw();
  print("doublefake_3.pdf");
}

// subtract all background MCs from data but the one specified with "notremove"
double get_fakes(const char * hname, const char * notremove)
{
  // get number of single fakes from data histogram
  plot(hname);
  rebin(5);
  legend();
  TH1D * hData = dataHisto();
  if (hData == 0) {
    ERROR("get_fakes() needs a data histogram");
    return 0;
  }
  double N = hData->Integral();
  INFO("Data events: " << N);

  TH1D * hSub = backgroundHisto(notremove, false);
  if (hSub == 0) {
    ERROR("get_fakes() needs a background histogram");
    return 0;
  }
  N = hSub->Integral();
  INFO("Background events : " << N);
  hData->Add(hSub, -1.);
  N = hData->Integral();
  INFO("Single fakes: " << N);

  return N;
}

/** Estimate number of single fakes
 *
 * Removes all MCs from data but wjets (which is estimated by this)
 *
 * @BUG: this also removes wwjets but contribution is negligible.
 */
double get_singlefakes(const char * hname)
{
  return get_fakes(hname, "wjets");
}

double get_doublefakes(const char * hname)
{
  return get_fakes(hname, "qcd");
}

/* Estimate single- and double-fakes and derive number of QCD and W+jets events. */
void fake_estimate(const char * sel, const char * hname, double & N_W, double & N_QCD)
{
  MakeCanvas();
  selection(Form("%s_singlefake", sel));
  cd(1);
  double N_sf = get_singlefakes(hname);
  selection(Form("%s_doublefake", sel));
  cd(2);
  double N_df = get_doublefakes(hname);
  print(Form("%s-fake_estimate.pdf", sel));
  INFO("N_sf = " << N_sf);
  INFO("N_df = " << N_df);
  // subtract double fakes from single fakes
  N_W   = N_sf - N_df;
  N_QCD = N_df;
  INFO("QCD    events: " << N_QCD);
  INFO("W+jets events: " << N_W);
}

struct syst_struct {
  const char * sel;
  const char * cfg;
};

void fakerate_systematics(int istart = 0, int iend = 999)
{
  // histogram name to be used for plotting
  const char * hname = "m_smuon";

  // reference values from standard analysis
  double N_W_ref;
  double N_QCD_ref;

  // get default values
  fake_estimate("default13", hname, N_W_ref, N_QCD_ref);

  // list of systematics
  const int nMax = 10;
  const syst_struct sel[nMax] = {
    { "reliso_03", "" },
    { "reliso_05", "" },
    { "reliso_06", "" },
    { "reliso_08", "" },

    { "jetptmin_50", "" },
    { "jetptmin_60", "" },
    { "jetptmin_80", "" },

    { "fakeratemethod_zero", "" },
    
    { "triggerbias_singlemu", "singlemu" },
    { "triggerbias_mu8_jet40", "mu8_jet40" }
  };

  // event counters
  double N_W[nMax];
  double N_QCD[nMax];

  // get all numbers
  for (int i = istart; i < TMath::Min(iend, nMax); i++) {
    setup(Form("../config/plot%s.cfg", sel[i].cfg));
    fake_estimate(sel[i].sel, hname, N_W[i], N_QCD[i]);
  }

  // compute differences
  for (int i = istart; i < TMath::Min(iend, nMax); i++) {
    N_QCD[i] -= N_QCD_ref;
    N_W[i] -= N_W_ref;
    // output table line
    cout << sel[i].sel << ": " 
	 << "QCD = " << 100.*N_QCD[i]/N_QCD_ref << "%, " 
	 << "W = " << 100.*N_W[i]/N_W_ref << "%" << endl;
  }
  
}
