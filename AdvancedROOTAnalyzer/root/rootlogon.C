{
  printf("Loading rootlogon.C...\n");

  // load my macros
  gROOT->LoadMacro(Form("%s/root/plot.C+", getenv("HOME")));
  gROOT->LoadMacro(Form("%s/root/stat.C+", getenv("HOME")));

  // set up everything
  setup("../config/plot.cfg");

  printf("rootlogon.C loaded successfully.\n");
}
