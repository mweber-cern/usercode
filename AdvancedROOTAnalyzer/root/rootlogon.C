{
  printf("Loading rootlogon.C...\n");

  // load my macros
  gROOT->LoadMacro(Form("%s/root/plot.C+", getenv("ARASYS")));
  gROOT->LoadMacro(Form("%s/root/stat.C+", getenv("ARASYS")));

  // set up everything
  setup("../config/plot.cfg");

  printf("rootlogon.C loaded successfully.\n");
}
