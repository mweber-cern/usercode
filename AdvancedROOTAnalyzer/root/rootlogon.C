{
  printf("Loading rootlogon.C...\n");

  // corrections to ugly ROOT standard style
  // gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);

  // load main macros
  gROOT->LoadMacro(Form("%s/root/plot.C+", getenv("ARASYS")));
  gROOT->LoadMacro(Form("%s/root/stat.C+", getenv("ARASYS")));

  // set up everything
  setup("../config/plot.cfg");

  // if you want to debug what is happening, uncomment next line
  // gLogLevel = 100;

  // load and apply CMS TDR style for plots
  gROOT->LoadMacro(Form("%s/root/tdrstyle.C", getenv("ARASYS")));
  setTDRStyle();

  printf("rootlogon.C loaded successfully.\n");
}
