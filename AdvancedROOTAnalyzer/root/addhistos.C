void addhistos(const Int_t number, const char * fname[])
{
  // adds all histograms (1d) from the given files and
  // writes them into one output file

  // input/output files
  TFile * infile[number];
  TFile * outfile = 0;

  if (number < 1) {
    printf("ERR: Give at least one file name\n");
  }

  // open one histogram file (used to read in histogram names)
  printf("Reading histograms from file %s\n", fname[0]);

  TFile * nameFile = TFile::Open(fname[0], "READ");
  if (!nameFile->IsOpen()) {
    printf("ERR: Opening name file %s\n", fname[0]);
    return;
  }

  // get list of keys with histogram names from file
  TList * keyList  = nameFile->GetListOfKeys();
  TList * nameList = new TList();
  TKey  * key;
  TIter   next(keyList);
  Int_t   nhisto = 0;
  while (key = (TKey *) next()) {
    nhisto++;
    // only store histogram names
    char * cname = key->GetName();
    if (cname == 0) {
      printf("Error in cname\n");
      return;
    }
    // only add 1d histograms
    if (!strcmp(key->GetClassName(), "TH1D"))
      nameList->AddLast(new TObjString(key->GetName()));
  }

  nameFile->Close();
  delete nameFile;

  char       * histo;
  TObjString   * hname;
  Int_t          nh;
  TH1D         * sumhisto;
  TH1D         * testhisto;

  TIter          next2(nameList);
  for (Int_t n = 0; n < number; n++) {
    // open all input files
    infile[n] = TFile::Open(fname[n], "READ");
    if (!infile[n]->IsOpen()) {
      printf("Error opening %s\n", fname);
      return;
    }
  }

  // open output file
  outfile = TFile::Open("join.root", "RECREATE");
  if (!outfile->IsOpen()) {
    printf("Error opening existing file %s\n", "join.root");
    return;
  }

  printf("Now copying histograms\n");
  next2.Reset();
  nh = 0;
  // loop over histograms in input file
  while (hname = (TObjString *) next2()) {
    nh++;
    printf("%3d%%...\b\b\b\b\b\b\b", (100 * nh) / nhisto);
    histo = hname->GetName();
    printf("histo %s\n", histo);

    // sum it up
    sumhisto = 0;
    for (Int_t n = 0; n < number; n++) {
      // get histogram from file
      testhisto = (TH1D *) infile[n]->Get(histo);
      if (testhisto == 0) {
	printf("histo %s not found in file %s\n", histo, fname[n]);
	continue;
      }
      if (sumhisto == 0) {
	// create new histogram for first mc
	sumhisto = new TH1D(*testhisto);
	if (sumhisto == 0) {
	  printf("Out of memory\n");
	  return;
	}
	sumhisto->SetDirectory(0);
      }
      else {
	// add to previous histos
	sumhisto->Add(testhisto);
      }
    }

    // save histogram in output file
    outfile->cd();
    sumhisto->Write();
    delete sumhisto;
  }

//    printf("Closing files\n");

  // close output file
  outfile->Close();
  delete outfile;
  outfile = 0;

  // close input files
  for (Int_t n = 0; n < number; n++) {
    delete infile[n];
    infile[n] = 0;
  }
}
