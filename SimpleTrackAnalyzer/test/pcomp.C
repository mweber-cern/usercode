void pcomp()
{
  gStyle->SetOptStat(0);
  TFile * f0T = TFile::Open("histo_0T.root");
  hp0T = (TH1F *) f0T->Get("hpt");
  hp0T->SetLineColor(kBlue);
  hp0T->Draw();
  TFile * f3T = TFile::Open("histo_3T.root");
  hp3T = (TH1F *) f3T->Get("hpt");
  hp3T->SetLineColor(kRed);
  hp3T->Draw("same");
  TLegend * leg = new TLegend(0.5, 0.67, 0.9, 0.88, 0, "brNDC"); // 60, 6, 100, 9);
  leg->AddEntry(hp0T, "ALCARECOTkAlCosmics0T", "l");  
  leg->AddEntry(hp3T, "ALCARECOTkAlCosmics (p > 4 GeV)", "l");
  leg->Draw();
  gPad->Print("0Tvs3T.pdf");
}
