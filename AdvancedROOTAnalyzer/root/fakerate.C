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

TH1D * get_subtracted_tight_loose_ratio(bool save=true, bool draw=true)
{
  if (gStart != gEnd) {
    ERROR("tight_loose_ratio() can only be called for one period!");
    return 0;
  }
  if (draw) {
    MakeCanvas();
    // MakeCanvas(2,2);
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
    }
    hTight3_data->Add(hTight3_back, -1.);
    nTightAfter -= nTight;
    hLoose3_data->Add(hLoose3_back, -1.);
    nLooseAfter -= nLoose;
  }
  // output some statistics
  INFO("After subtraction:");
  INFO("Count: Tight muons: " << nTightAfter << ", loose muons: " << nLooseAfter);
  INFO("Histo: Tight muons: " << hTight3_data->Integral() 
       << ", loose muons: " << hLoose3_data->Integral());

  // now work in 2D-Projection
  TH2D * hTight2 = (TH2D *) hTight3_data->Project3D("yx");
  TH2D * hLoose2 = (TH2D *) hLoose3_data->Project3D("yx");
  TH2D * hRatio2 = new TH2D(*hTight2);
  hRatio2->SetDirectory(0);
  hRatio2->SetName("hRatio2");
  hRatio2->SetTitle("Tight/Loose ratio (T/L) (2D)");

  if (draw) {
    cd(3); hLoose2->Draw("lego2");
    cd(4); hTight2->Draw("lego2");
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
      }
      if (nRatio != 0) {
	DEBUG("T/L (" << get_bin_edges(hRatio2, i, j) << ") " 
	      << nRatio << " +- " << error);
      }
    }
  }
  if (draw) {
    cd(1);
    hRatio2->Draw("lego2");
  }

  // now do the same in 1D-Projection
  TH1D * hTight1 = (TH1D *) hTight3_data->Project3D("x");
  TH1D * hLoose1 = (TH1D *) hLoose3_data->Project3D("x");
  TH1D * hRatio1 = new TH1D(*hTight1);
  hRatio1->SetDirectory(0);
  hRatio1->SetName("hRatio1");
  hRatio1->SetTitle("Tight/Loose ratio (T/L)");

  // divide with binomial errors
  hRatio1->Divide(hTight1, hLoose1, 1., 1., "B");

  if (save) {
    const char * fakeRateFileName = "FakeRate.root";
    TFile * fakeRateFile = new TFile(fakeRateFileName, "RECREATE");
    if (fakeRateFile == 0 || !fakeRateFile->IsOpen()) {
      cerr << "ERR: Could not open fake rate file " << fakeRateFileName << endl;
      return 0;
    }
    hRatio1->Write();
    hRatio2->Write();
    fakeRateFile->Close();
    delete fakeRateFile;
  }

  if (draw) {
    cd(2);
    hRatio1->Draw();
  }
  return hRatio1;
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
  TH3D * hTight3 = (TH3D *) f.Get("h3_0_TightMuons");
  TH3D * hLoose3 = (TH3D *) f.Get("h3_0_LooseMuons");
  if (hTight3 == 0) {
    cerr << "Could not get hTight3"; 
    return 0;
  }
  if (hLoose3 == 0) {
    cerr << "Could not get hLoose3"; 
    return 0;
  }
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
  hRatio2->Draw("LEGO2");
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
    hRatio1->Write();
    hRatio2->Write();
    hRatio3->Write();
    fakeRateFile->Close();
    delete fakeRateFile;
  }

  // return one-dimensional plot
  return hRatio1;
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
  return hRatio1;
}

void ratioplot(const char * sele = "tightloose5")
{
  DEBUG("ratioplot() start");
  MakeCanvas();
  if (sele)
    selection(sele);

  DEBUG("reading histograms");
  // "plot" histograms, i.e. read into memory
  cd(1);
  plot3("TightMuons");
  cd(2);
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
    if (hTight3 == 0 || hLoose3 == 0)
      continue;
    // if joined, it was already considered in previous iteration
    if (gProcess[process].join)
      continue;
    if (strncmp(gProcess[process].fname, "qcd", 3) &&
	strncmp(gProcess[process].fname, "dyll", 4) &&
	strncmp(gProcess[process].fname, "ttjets", 3) &&
	strncmp(gProcess[process].fname, "data", 3))
      continue;
    DEBUG("Creating histos for process " << gProcess[process].fname);
    hTightSum = new TH3D(*hTight3);
    hLooseSum = new TH3D(*hLoose3);
    // check if histograms are joined -> we need to add statistics
    for (Int_t j = i+1; j < gMaxProcess; j++) {
      Int_t proc = gOrder[0][j];
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
  hdata_subtracted->SetMarkerColor(kBlue);
  hdata_subtracted->SetMarkerStyle(8);
  hdata_subtracted->Draw("epsame");
  leg->AddEntry(hdata_subtracted, "data subtr.", "ep");
  leg->Draw();
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
  plot("nTL_ht");
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
  plot("TL_mt");
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

TH1D * get_doublefakes(const char * hname)
{
  if (gStart != gEnd) {
    ERROR("Cannot work for more than one period");
    return 0;
  }
  MakeCanvas();
  TFile f(Form("%s/%s", getpath(gStart), "data_doublemu.root"));
  TH1D * hmu = (TH1D *) f.Get(hname);
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
  return hcut;
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

void lumicorrection()
{
  plot("nTL_met");
  correctfactor("qcd", 0, 100);
  plot("nTL_met");
  correctfactor("tt", 130, 500);
  plot("nTL_met");
  correctfactor("qcd", 0, 100);
  plot("nTL_met");
  correctfactor("dyll", 0, 100);
  plot("nTL_met");
  correctfactor("qcd", 0, 100);
  plot("nTL_met");
  cd(2);
  plot("nTL_met");
  logy();
}
