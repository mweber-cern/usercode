void setNiceColorPalette () {
    const Int_t NRGBs = 5;
    const Int_t NCont = 999;


    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

void rootlogon () {
  printf("Loading rootlogon.C...\n");

  // corrections to ugly ROOT standard style
  gROOT->SetStyle("Plain");
  setNiceColorPalette();

  // load and apply CMS TDR style for plots
  gROOT->LoadMacro(Form("%s/root/tdrstyle.C", getenv("ARASYS")));
  setTDRStyle();

  // load main macros
  gROOT->LoadMacro(Form("%s/root/plot.C+", getenv("ARASYS")));
  gROOT->LoadMacro(Form("%s/root/stat.C+", getenv("ARASYS")));
  gROOT->LoadMacro(Form("%s/root/xsection.C+", getenv("ARASYS")));
  gROOT->LoadMacro(Form("%s/root/fakerate.C+", getenv("ARASYS")));
  gROOT->LoadMacro(Form("%s/root/rpv.C+", getenv("ARASYS")));
  gROOT->LoadMacro(Form("%s/root/systematics.C+", getenv("ARASYS")));
  
  // set up everything
  setup("plot.cfg");

  // if you want to debug what is happening, uncomment next line
  // gLogLevel = 100;

  printf("rootlogon.C loaded successfully.\n");
}
