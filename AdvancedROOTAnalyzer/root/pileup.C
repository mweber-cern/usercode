void generate_pileup_histograms()
{
  hfall11_fullsim = generate_pileup_fall11_fullsim(true);
  hsummer12_s7 = generate_pileup_summer12_S7_fullsim(true);
  hsummer12_s10 = generate_pileup_summer12_S10_fullsim(true);
  hsummer12_fastsim = generate_pileup_summer12_fastsim(true);
  hstartup12_fastsim = generate_pileup_startup12_fastsim(true);
}

void compare_pileup_histograms()
{
  TH1D * hfall11_fullsim = generate_pileup_fall11_fullsim(false);
  TH1D * hsummer12_s7 = generate_pileup_summer12_S7_fullsim(false);
  TH1D * hsummer12_s10 = generate_pileup_summer12_S10_fullsim(false);
  TH1D * hsummer12_fastsim = generate_pileup_summer12_fastsim(false);
  TH1D * hstartup12_fastsim = generate_pileup_startup12_fastsim(false);

  MakeCanvas(1,1);
  gCanvas->DrawFrame(0, 0, 60, 0.07);
  hsummer12_fastsim->SetLineColor(kBlue);
  hsummer12_fastsim->SetFillStyle(0);
  hsummer12_fastsim->Draw("lsame");
  hstartup12_fastsim->SetLineColor(kGreen);
  hstartup12_fastsim->SetFillStyle(0);
  hstartup12_fastsim->Draw("lsame");
  hsummer12_s10->SetLineColor(kRed);
  hsummer12_s10->SetLineStyle(kDotted);
  hsummer12_s10->SetFillStyle(0);
  hsummer12_s10->Draw("lsame");
  
  TFile * f = TFile::Open("../config/histAllData12.root");
  hdata = (TH1D * ) f->Get("pileup");
  hdata->SetDirectory(0);
  delete f;
  hdata->SetLineColor(kBlack);
  hdata->SetFillStyle(0);
  hdata->DrawNormalized("lsame");

  f = TFile::Open("/user/mweber/fullsim/2012/signal_fullsim_m0_eq_1000_m12_eq_1200.root");
  TH1D * hfullsim = (TH1D * ) f->Get("h1_0_pu_TrueNrInter");
  hfullsim->SetMarkerStyle(22);
  hfullsim->DrawNormalized("epsame");

  f = TFile::Open("/user/mweber/fastsim/2012/signal_fastsim_m0_eq_1000_m12_eq_1200.root");
  TH1D * hfastsim = (TH1D * ) f->Get("h1_0_pu_TrueNrInter");
  hfastsim->SetMarkerStyle(23);
  hfastsim->DrawNormalized("epsame");
  
  TLegend * leg = new TLegend(0.65, 0.7, 0.93, 0.93);
  setopt(leg);
  leg->AddEntry(hsummer12_fastsim, "Summer12 (Fastsim)", "l");
  leg->AddEntry(hstartup12_fastsim, "Startup12 (Fastsim)", "l");
  leg->AddEntry(hsummer12_s10, "Summer12 S10 (Fullsim)", "l");
  leg->AddEntry(hdata, "Data", "l");
  leg->AddEntry(hfullsim, "Used fullsim scenario", "ep");
  leg->AddEntry(hfastsim, "Used fastsim scenario", "ep");
  leg->Draw();
  
}

TH1D * generate_pileup_summer12_S7_fullsim(bool save)
{
  double Summer12_S7[60] = {
    2.344E-05,
    2.344E-05,
    2.344E-05,
    2.344E-05,
    4.687E-04,
    4.687E-04,
    7.032E-04,
    9.414E-04,
    1.234E-03,
    1.603E-03,
    2.464E-03,
    3.250E-03,
    5.021E-03,
    6.644E-03,
    8.502E-03,
    1.121E-02,
    1.518E-02,
    2.033E-02,
    2.608E-02,
    3.171E-02,
    3.667E-02,
    4.060E-02,
    4.338E-02,
    4.520E-02,
    4.641E-02,
    4.735E-02,
    4.816E-02,
    4.881E-02,
    4.917E-02,
    4.909E-02,
    4.842E-02,
    4.707E-02,
    4.501E-02,
    4.228E-02,
    3.896E-02,
    3.521E-02,
    3.118E-02,
    2.702E-02,
    2.287E-02,
    1.885E-02,
    1.508E-02,
    1.166E-02,
    8.673E-03,
    6.190E-03,
    4.222E-03,
    2.746E-03,
    1.698E-03,
    9.971E-04,
    5.549E-04,
    2.924E-04,
    1.457E-04,
    6.864E-05,
    3.054E-05,
    1.282E-05,
    5.081E-06,
    1.898E-06,
    6.688E-07,
    2.221E-07,
    6.947E-08,
    2.047E-08
  };

  TFile * file;
  if (save) { 
    file = new TFile("../config/Pileup_truth_Summer12_MC_fullsim_S7.root", "RECREATE");
  }
  TH1D  * h1 = new TH1D("pileup", "Pileup Summer 12 MC fullsim S7", 60, 0, 60);
  for (int i = 0; i < 60; i++) {
    h1->SetBinContent(i+1, Summer12_S7[i]);
  }
  if (save) {
    h1->Write();
    h1->SetDirectory(0);
    file->Close();
    delete file;
  }
  return h1;
}

TH1D * generate_pileup_summer12_S10_fullsim(bool save)
{
  // Distribution used for S10 Summer2012 MC.
  Double_t Summer12_S10[60] = {
    2.560E-06,
    5.239E-06,
    1.420E-05, 
    5.005E-05,
    1.001E-04,
    2.705E-04,
    1.999E-03,
    6.097E-03,
    1.046E-02,
    1.383E-02,
    1.685E-02,
    2.055E-02,
    2.572E-02,
    3.262E-02,
    4.121E-02,
    4.977E-02,
    5.539E-02,
    5.725E-02,
    5.607E-02,
    5.312E-02,
    5.008E-02,
    4.763E-02,
    4.558E-02,
    4.363E-02,
    4.159E-02,
    3.933E-02,
    3.681E-02,
    3.406E-02,
    3.116E-02,
    2.818E-02,
    2.519E-02,
    2.226E-02,
    1.946E-02,
    1.682E-02,
    1.437E-02,
    1.215E-02,
    1.016E-02,
    8.400E-03,
    6.873E-03,
    5.564E-03,
    4.457E-03,
    3.533E-03,
    2.772E-03,
    2.154E-03,
    1.656E-03,
    1.261E-03,
    9.513E-04,
    7.107E-04,
    5.259E-04,
    3.856E-04,
    2.801E-04,
    2.017E-04,
    1.439E-04,
    1.017E-04,
    7.126E-05,
    4.948E-05,
    3.405E-05,
    2.322E-05,
    1.570E-05,
    5.005E-06
  };

  TFile * file;
  if (save) { 
    file = new TFile("../config/Pileup_truth_Summer12_MC_fullsim_S10.root", "RECREATE");
  }
  TH1D * h1 = new TH1D("pileup", "Pileup Summer 12 MC fullsim S10", 60, 0, 60);
  for (int i = 0; i < 60; i++) {
    h1->SetBinContent(i+1, Summer12_S10[i]);
  }
  if (save) {
    h1->Write();
    h1->SetDirectory(0);
    file->Close();
    delete file;
  }
  return h1;
}

TH1D * generate_pileup_summer12_fastsim(bool save)
{
  Double_t Summer12_fastsim[60] = {
    2.560E-06,
    5.239E-06,
    1.420E-05,
    5.005E-05,
    1.001E-04,
    2.705E-04,
    1.999E-03,
    6.097E-03,
    1.046E-02,
    1.383E-02,
    1.685E-02,
    2.055E-02,
    2.572E-02,
    3.262E-02,
    4.121E-02,
    4.977E-02,
    5.539E-02,
    5.725E-02,
    5.607E-02,
    5.312E-02,
    5.008E-02,
    4.763E-02,
    4.558E-02,
    4.363E-02,
    4.159E-02,
    3.933E-02,
    3.681E-02,
    3.406E-02,
    3.116E-02,
    2.818E-02,
    2.519E-02,
    2.226E-02,
    1.946E-02,
    1.682E-02,
    1.437E-02,
    1.215E-02,
    1.016E-02,
    8.400E-03,
    6.873E-03,
    5.564E-03,
    4.457E-03,
    3.533E-03,
    2.772E-03,
    2.154E-03,
    1.656E-03,
    1.261E-03,
    9.513E-04,
    7.107E-04,
    5.259E-04,
    3.856E-04,
    2.801E-04,
    2.017E-04,
    1.439E-04,
    1.017E-04,
    7.126E-05,
    4.948E-05,
    3.405E-05,
    2.322E-05,
    1.570E-05,
    5.005E-06
  };
  TFile * file;
  if (save) { 
    file = new TFile("../config/Pileup_truth_Summer12_MC_fastsim.root", "RECREATE");
  }
  TH1D * h1 = new TH1D("pileup", "Pileup Summer 12 MC fastsim", 60, 0, 60);
  for (int i = 0; i < 60; i++) {
    h1->SetBinContent(i+1, Summer12_fastsim[i]);
  }
  if (save) {
    h1->Write();
    h1->SetDirectory(0);
    file->Close();
    delete file;
  }
  return h1;
}

TH1D * generate_pileup_startup12_fastsim(bool save)
{
  double Startup12_inTimeOnly[60] = {
    2.344E-05,
    2.344E-05,
    2.344E-05,
    2.344E-05,
    4.687E-04,
    4.687E-04,
    7.032E-04,
    9.414E-04,
    1.234E-03,
    1.603E-03,
    2.464E-03,
    3.250E-03,
    5.021E-03,
    6.644E-03,
    8.502E-03,
    1.121E-02,
    1.518E-02,
    2.033E-02,
    2.608E-02,
    3.171E-02,
    3.667E-02,
    4.060E-02,
    4.338E-02,
    4.520E-02,
    4.641E-02,
    4.735E-02,
    4.816E-02,
    4.881E-02,
    4.917E-02,
    4.909E-02,
    4.842E-02,
    4.707E-02,
    4.501E-02,
    4.228E-02,
    3.896E-02,
    3.521E-02,
    3.118E-02,
    2.702E-02,
    2.287E-02,
    1.885E-02,
    1.508E-02,
    1.166E-02,
    8.673E-03,
    6.190E-03,
    4.222E-03,
    2.746E-03,
    1.698E-03,
    9.971E-04,
    5.549E-04,
    2.924E-04,
    1.457E-04,
    6.864E-05,
    3.054E-05,
    1.282E-05,
    5.081E-06,
    1.898E-06,
    6.688E-07,
    2.221E-07,
    6.947E-08,
    2.047E-08
  };

  TFile * file;
  if (save) { 
    file = new TFile("../config/Pileup_truth_Startup12_MC_fastsim.root", "RECREATE");
  }
  TH1D * h1 = new TH1D("pileup", "Pileup Startup 12 MC fastsim", 60, 0, 60);
  for (int i = 0; i < 60; i++) {
    h1->SetBinContent(i+1, Startup12_inTimeOnly[i]);
  }
  if (save) {
    h1->Write();
    h1->SetDirectory(0);
    file->Close();
    delete file;
  }
  return h1;
}
 

TH1D * generate_pileup_fall11_fullsim(bool save)
{
  // Distribution used for Fall2011 Fullsim MC.
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

  TFile * file;
  if (save) { 
    file = new TFile("../config/Pileup_truth_Fall11_MC_fullsim.root", "RECREATE");
  }
  TH1D * h1 = new TH1D("pileup", "Pileup MC Fall11 distribution", 50, 0, 50);
  for (int i = 0; i < 50; i++) {
    h1->SetBinContent(i+1, Fall2011[i]);
  }
  if (save) {
    h1->Write();
    h1->SetDirectory(0);
    file->Close();
    delete file;
  }
  return h1;
}
