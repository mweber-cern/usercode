#include <iostream>
#include <typeinfo>

using namespace std;

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TFile.h"
#include "THashList.h"
#include "TString.h"
#include "TMath.h"

#include "plot.h"
#include "stat.h"


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

TH1D * tight_loose_ratio(const char * process = "data_doublemu", bool save=true)
{
  if (gStart != gEnd) {
    ERROR("tight_loose_ratio() can only be called for one period!");
    return 0;
  }
  TFile f(Form("%s/%s.root", getpath(gStart), process));
  if (!f.IsOpen()) {
    return 0;
  }
  TH3D * hTight3 = (TH3D *) f.Get("h3_TightMuons");
  TH3D * hLoose3 = (TH3D *) f.Get("h3_LooseMuons");
  TH3D * hRatio3 = new TH3D(*hTight3);
  hRatio3->SetDirectory(0);
  hRatio3->SetName("hRatio3");
  hRatio3->SetTitle("Tight/Loose ratio (T/L) (3D)");

  // divide with binomial errors
  hRatio3->Divide(hTight3, hLoose3, 1., 1., "B");
  // some sanity checks
  for (Int_t i = 1;  i < hRatio3->GetNbinsX()+1; i++) {
    for (Int_t j = 1; j < hRatio3->GetNbinsY()+1; j++) {
      for (Int_t k = 1; k < hRatio3->GetNbinsZ()+1; k++) {
	Double_t nLoose = hLoose3->GetBinContent(i, j, k);
	Double_t nTight = hTight3->GetBinContent(i, j, k);
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

  // now do the same in 2D-Projection
  TH2D * hTight2 = (TH2D *) hTight3->Project3D("yx");
  TH2D * hLoose2 = (TH2D *) hLoose3->Project3D("yx");
  TH2D * hRatio2 = new TH2D(*hTight2);
  hRatio2->SetDirectory(0);
  hRatio2->SetName("hRatio2");
  hRatio2->SetTitle("Tight/Loose ratio (T/L) (2D)");

  // divide with binomial errors
  hRatio2->Divide(hTight2, hLoose2, 1., 1., "B");
  // some sanity checks
  for (Int_t i = 1;  i < hRatio2->GetNbinsX()+1; i++) {
    for (Int_t j = 1; j < hRatio2->GetNbinsY()+1; j++) {
      Double_t nLoose = hLoose2->GetBinContent(i, j);
      Double_t nTight = hTight2->GetBinContent(i, j);
      Double_t nRatio = hRatio2->GetBinContent(i, j);
      Double_t error  = hRatio2->GetBinError(i, j);

      if (nLoose < nTight || nLoose < 0) {
	cerr << "ERR: nLoose is wrong in your histogram" << endl;
	cout << "INFO: nLoose(" << get_bin_edges(hRatio2, i, j) << ") " << nLoose << endl;
	cout << "INFO: nTight(" << get_bin_edges(hRatio2, i, j) << ") " << nTight << endl;
      }
      if (nLoose != 0 && nRatio != nTight/nLoose) {
	cerr << "ERR: ROOT calculation went wrong" << endl;
      }
      if (nRatio < 0 || nRatio >= 1.) {
	cerr << "ERR: Found unreasonable values for T/L ratio" << endl;
      }
      if (nRatio != 0) {
	cout << "INFO: T/L (" << get_bin_edges(hRatio2, i, j) << ") " 
	     << nRatio << " +- " << error << endl;
      }
    }
  }

  // now do the same in 1D-Projection
  TH1D * hTight1 = (TH1D *) hTight3->Project3D("x");
  TH1D * hLoose1 = (TH1D *) hLoose3->Project3D("x");
  TH1D * hRatio1 = new TH1D(*hTight1);
  hRatio1->SetDirectory(0);
  hRatio1->SetName("hRatio1");
  hRatio1->SetTitle("Tight/Loose ratio (T/L)");

  // divide with binomial errors
  hRatio1->Divide(hTight1, hLoose1, 1., 1., "B");

  // make plots
  cout << "Period: " << gPeriod[gStart] << " Ratio OK." << endl;
  MakeCanvas();
  cd(1);
  hRatio1->Draw("le");
  cd(2);
  hRatio2->Draw("LEGO");
  setopt(hRatio2);
  hRatio2->GetXaxis()->SetTitle("#mu p_{T} [GeV]");
  hRatio2->GetYaxis()->SetTitle("#mu #eta");
  hRatio2->GetZaxis()->SetTitle("T/L ratio");

  if (save) {
    const char * fakeRateFileName = "FakeRate.root";
    TFile * fakeRateFile = new TFile(fakeRateFileName, "RECREATE");
    if (fakeRateFile == 0 || !fakeRateFile->IsOpen()) {
      cerr << "ERR: Could not open fake rate file " << fakeRateFileName << endl;
      return 0;
    }
    hRatio2->Write();
    hRatio3->Write();
    fakeRateFile->Close();
    delete fakeRateFile;
  }

  // return one-dimensional plot
  return hRatio1;
}

void ratioplot()
{
  selection("tightloose");
  period("2011");

  TCanvas * c1 = new TCanvas("c1", "Tight-to-Loose ratio plots", 10, 10, (int) 600.*TMath::Sqrt(2.), 600);
  setopt(c1);

  // get histograms for each process
  bool first = true;
  TLegend * leg = new TLegend(0.7, 0.5, 0.95, 0.95);
  for (Int_t i = 0; i < gMaxProcess; i++) {
    TH1D * histo = tight_loose_ratio(gProcess[i].fname, false);
    if (histo == 0)
      continue;
    setopt(histo);
    histo->SetLineColor(gProcess[i].lcolor);
    histo->SetLineStyle(gProcess[i].lstyle);
    leg->AddEntry(histo, "l", gProcess[i].fname);
    c1->cd();
    if (first) {
      histo->Draw("histo");
      first = false;
    }
    else {
      histo->Draw("histosame");
    }
  }
}

void tightlooseplots(const char * sele = "tightloose")
{
  selection(sele);

  // plot 1
  cd(1);
  plot("nTL_met");
  logy();
  min(0.1);
  arrow(50.);
  legend(1e-3);
  cd(2);
  plot("nTL_jetpt");
  logy();
  min(0.1);
  arrow(40.);
  legend(1e-3);
  print("tightloose_1.pdf");

  // plot 2
  cd(1);
  plot("nTL_HT");
  min(0.1);
  logy();
  arrow(200);
  legend(1e-3);

  cd(2);
  plot("TL_nloose");
  min(0.1);
  logy();
  arrow(1);
  arrow(2);
  legend(1e-3);
  print("tightloose_2.pdf");
  
  // plot 3
  cd(1);
  plot("TL_jetdphi");
  min(0.1);
  logy();
  arrow(1.);
  legend(1e-3);

  cd(2);
  plot("TL_MT");
  logy();
  min(0.1);
  legend(1e-3);
  print("tightloose_3.pdf");  

  // plot 4
  cd(1);
  plot("TL_metdphi");
  logy();
  min(0.1);
  legend(1e-3);

  cd(2);
  plot("TL_zmass");
  logy();
  min(0.1);
  legend(1e-3);

  print("tightloose_4.pdf");  

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

double get_doublefakes()
{
  MakeCanvas();
  TFile f("/user/mweber/doublefake/2011/data_doublemu.root");
  TH1D * hmu = (TH1D *) f.Get("h1_m_mumu");
  hmu->SetDirectory(0);
  setopt(hmu);
  hmu->Rebin(5);
  hmu->Draw();
  double N = hmu->Integral();
  cout << "Found " << N << " events in final selection" << endl;
  cd(2);
  TH1D * hcut = (TH1D *) f.Get("h1_cutflow");
  hcut->SetDirectory(0);
  setopt(hcut);
  hcut->SetMarkerStyle(8);
  hcut->SetMaximum(70.);
  hcut->Draw("ep");
  N = hcut->GetBinContent(12.);
  cout << "Got " << N << " events from cut flow" << endl;
  return N;
}

const char * stages[] = { 
"weight","pileup rew.","redo skim","trigger","objectID","cleaning","tightloose","jetid","doublefake","singlefake","standard","muonID","met","m_mumu","charge"
};

void check_all_files()
{
  for (int i = 0; i < 167; i++) {
    const char * fname = Form("/user/mweber/doublefake2/2011/data_doublemu_%d.root", i);
    TFile * f = new TFile(fname);
    if (f == 0 || !f->IsOpen()) {
      cout << "Problem opening file " << fname << endl;
      return;
    }
    cout << fname << endl;
    TH1D * h1 = (TH1D *) f->Get("cutflow");
    if (h1 == 0) {
      cout << "Error getting histogram from file" << endl;
      return;
    }
    if (h1->GetNbinsX() != 15) {
      cout << "histogram has " << h1->GetNbinsX() << " bins " << endl;
    }
    // iterate over bin names
    THashList * hl = h1->GetXaxis()->GetLabels();
    int j = 0;
    TIter next(hl);
    TObjString * s;
    while ((s = (TObjString *) next())) {
      if (strcmp(s->GetString().Data(), stages[j])) {
	cout << "string mismatch: expected: " << stages[j] 
	     << " but got instead: " << s->GetString().Data() << endl;
      }
      j++;
    }
    if (h1 ->GetBinContent(7) != 0) {
      cout << "bin 7 contains " << h1->GetBinContent(7) << "events!" << endl;
    }
    delete f;
  }
}

