void plot()
{
  TFile f("histo.root");
  TCanvas * c1 = new TCanvas("c1", "c1", 0, 0, 600, 600);
  hfxy->SetMarkerStyle(29);
  hfxy->SetMarkerColor(kBlue);
  hfxy->Draw();
  hlxy->SetMarkerStyle(29);
  hlxy->SetMarkerColor(kRed);
  hlxy->Draw("same");
  c1->Print("c1.pdf");
  gSystem->Exit(0);
}
