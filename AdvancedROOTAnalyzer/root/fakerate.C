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

void tight_loose_ratio(const char * process = "qcd1520mu")
{
  for (Int_t period = 0; period < gMaxPeriod; period++) {
    TFile f(Form("%s/%s.root", getpath(period), process));
    if (!f.IsOpen()) {
      return;
    }
    TH3D * hTight3 = f.Get("h3_TightMuons");
    TH3D * hLoose3 = f.Get("h3_LooseMuons");
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
	    cout << "INFO: nLoose(" << get_bin_edges(hRatio, i, j, k) ") " << nLoose << endl;
	    cout << "INFO: nTight(" << get_bin_edges(hRatio, i, j, k) ") " << nTight << endl;
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
    hTight2 = (TH2D *) hTight3->Project3D("xy");
    hLoose2 = (TH2D *) hLoose3->Project3D("xy");
    TH2D * hRatio2 = new TH2D(*hTight2);
    hRatio2->SetDirectory(0);
    hRatio2->SetName("hRatio2");
    hRatio2->SetTitle("Tight/Loose ratio (T/L) (2D)");

    // divide with binomial errors
    hRatio2->Divide(hTight2, hLoose2, 1., 1., "B");
    // some sanity checks
    for (Int_t i = 1;  i < hRatio2->GetNbinsX()+1; i++) {
      for (Int_t j = 1; j < hRatio2->GetNbinsY()+1; j++) {
	Double_t nLoose = hLoose2->GetBinContent(i, j, k);
	Double_t nTight = hTight2->GetBinContent(i, j, k);
	Double_t nRatio = hRatio2->GetBinContent(i, j, k);
	Double_t error  = hRatio2->GetBinError(i, j, k);

	if (nLoose < nTight || nLoose < 0) {
	  cerr << "ERR: nLoose is wrong in your histogram" << endl;
	  cout << "INFO: nLoose(" << get_bin_edges(hRatio, i, j) ") " << nLoose << endl;
	  cout << "INFO: nTight(" << get_bin_edges(hRatio, i, j) ") " << nTight << endl;
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

    cout << "Period: " << gPeriod[period] << " Ratio OK." << endl;
    MakeCanvas();
    cd(1);
    hRatio3->Draw();
    cd(2);
    hRatio2->Draw();

    const char * fakeRateFileName = "FakeRate.root";
    TFile * fakeRateFile = new TFile(fakeRateFileName, "RECREATE");
    if (fakeRateFile == 0 || !fakeRateFile->IsOpen()) {
      cerr << "ERR: Could not open fake rate file " << fakeRateFileName << endl;
      return;
    }
    hRatio2->Write();
    hRatio3->Write();
    fakeRateFile->Close();
    delete fakeRateFile;
  }
}
