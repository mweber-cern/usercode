#include "Analysis.h"

// ROOT includes
#include <TH2.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>

// C / C++ includes
#include <iostream>
#include <sstream>
#include <cmath>
#include <typeinfo>

// local includes
#include "Combinations.h"
#include "Utilities.h"
#include "TCutList.h"

using namespace std;

// return true if a ROOT file named "fileName" exists
Bool_t check_root_file(const char * fileName)
{
  // get rid of special shell variables like ~ or $ expansion
  const char * fullName = gSystem->ExpandPathName(fileName);
  TFile * rootFile = TFile::Open(fullName, "READ");
  if (rootFile == 0) {
    return false;
  }
  if (!rootFile->IsOpen()) {
    delete rootFile;
    return false;
  }
  delete rootFile;
  return true;
}

Analysis::Analysis(TTree & inputTree, TTree & outputTree, TEnv & cfgFile) 
  : TreeContent(& inputTree) , fInputTree(inputTree), fOutputTree(outputTree), 
    fCfgFile(cfgFile) 
{
  //////////////////////////////////////////////////////////////////////
  // Configuration options

  // fill output tree?
  fFill = fCfgFile.GetValue("FillTree", true);
  // Which type of input (signal, background, data, MC)?
  fType = fCfgFile.GetValue("FileType", "data");
  if (fType != "signal" && fType != "background" && fType != "data" && fType != "mc") 
    THROW("FileType \"" + fType + "\" is not an allowed value");
  // maximum number of events to be processed (-1: process all)
  fMaxEvents = fCfgFile.GetValue("MaxEvents", -1);
  // maximum tree size for output file
  fMaxTreeSize = fCfgFile.GetValue("MaxTreeSize", 100000000);
  INFO("Setting maximum tree size to " << fMaxTreeSize);
  // find duplicates?
  fFindDuplicates = fCfgFile.GetValue("FindDuplicates", true);
  // dump all event information?
  fDumpAll = fCfgFile.GetValue("DumpAll", false);
  // dump MC truth?
  fDumpTruth = fCfgFile.GetValue("DumpTruth", false);
  fSkimActive = fCfgFile.GetValue("SkimActive", true);
  fSkimMuons = fCfgFile.GetValue("SkimMuons", 2.);
  fSkimMuoptfirst = fCfgFile.GetValue("SkimMuoptfirst", 15.);
  fSkimMuoptother = fCfgFile.GetValue("SkimMuoptother", 7.);

  //////////////////////////////////////////////////////////////////////
  // analysis variables
  fIsSignal = false;
  fTrigger = new TString(30);

  //////////////////////////////////////////////////////////////////////
  // initialize pileup reweighting
  fPileupDataFile = gSystem->ExpandPathName(fCfgFile.GetValue("PileupDataFile", "~/config/2011_data_pileup.root"));
  fPileupMCFile = gSystem->ExpandPathName(fCfgFile.GetValue("PileupMCFile", "~/config/Fall11_MC_pileup_truth.root"));
  fPileupWeightFile = gSystem->ExpandPathName(fCfgFile.GetValue("PileupWeightFile", "~/config/Weight3D.root"));
  // see if weight file exits
  if (!check_root_file(fPileupWeightFile)) {
    INFO("Could not open pileup reweighting file " << fPileupWeightFile);
    if (!check_root_file(fPileupMCFile)) {
      ERROR("Could not find reweighting file " << fPileupMCFile);
      THROW("Missing required file: " + string(fPileupMCFile));
    }
    if (!check_root_file(fPileupDataFile)) {
      ERROR("Could not find reweighting file " << fPileupDataFile);
      THROW("Missing required file: " + string(fPileupDataFile));
    }
    LumiWeights_ = reweight::LumiReWeighting(fPileupMCFile, fPileupDataFile, "pileup","pileup");
    LumiWeights_.weight3D_init();
  } 
  else {
    LumiWeights_ = reweight::LumiReWeighting();
    LumiWeights_.weight3D_init(fPileupWeightFile);
  }
}

Analysis::~Analysis()
{
  delete fTrigger;
}

void Analysis::ReconstructFinalState(map<char, vector<int> > & particles, int vertex)
{
  DEBUG("Analysis::ReconstructFinalState(" << vertex << ")");
  // loop over all particles
  for (int j = 0; j < truth_n; j++) {
    // only take particles from this vertex
    if (truth_bvtxid[j] != vertex)
      continue;
    Int_t SusyCode = abs(truth_pdgid[j]) / 1000000;

    // check daughters
    int pdgid = abs(truth_pdgid[j]);
    // take all leptons
    if (pdgid == 11 || pdgid == 13 || pdgid == 15) {
      DEBUG("lepton j = " << j << ", pdgid = " << pdgid);
      particles['l'].push_back(j);
    }
    else if (pdgid == 12 || pdgid == 14 || pdgid == 16) {
      DEBUG("neutrino j = " << j << ", pdgid = " << pdgid);
      particles['n'].push_back(j);
    }
    // take all quarks but at vertex 1 (SUSYAna feature)
    else if ((pdgid >= 1 && pdgid <= 6)) {
      if (vertex != 1) {
	DEBUG("quark j = " << j << ", pdgid = " << pdgid);
	particles['q'].push_back(j);
      }
    }
    // check if there is another SUSY particle
    else if (SusyCode == 1 || SusyCode == 2) {
      DEBUG("sparticle j = " << j << ", pdgid = " << pdgid);
      particles['s'].push_back(j);
      ReconstructFinalState(particles, truth_evtxid[j]);
    }
    // W, Z, Higgs -> decaying in two more particles (leptons/quarks)
    else if (pdgid == 23 || pdgid == 24 || pdgid == 25) {
      DEBUG("boson j = " << j << ", pdgid = " << pdgid);
      particles['b'].push_back(j);
    }
    else if (pdgid == 0) {
      // ignore - Herwig special code
    }
    else {
      ERROR("Unknown particle found: " << j << ", pdgid = " << truth_pdgid[j]);
      TruthDump();
      THROW("Unknown particle found");
    }
  }
}

truth_pt_comparator::truth_pt_comparator(const Analysis & analysis)
  : fAnalysis(analysis)
{
}


bool truth_pt_comparator::operator()(int i, int j)
{
  return fAnalysis.truth_pt[i] > fAnalysis.truth_pt[j];
}

void Analysis::SignalStudy()
{
  // muons, neutrinos and quarks in final state
  int nLeptons = 0;
  int nGauginos = 0;
  int nSleptons = 0;
  int nBosons   = 0;

  // slepton charge
  int charge = 0;

  for (int i = 0; i < truth_n; i++) {
    // find particles produced at first vertex
    if (truth_bvtxid[i] != 1)
      continue;
    switch (abs(truth_pdgid[i])) {
      // muon
      case 13:
	nLeptons++;
	if (truth_pdgid[i] > 0)
	  charge--;
	else
	  charge++;
	break;
	// muon neutrino
      case 14:
	nLeptons++;
	break;
	// ignore quarks and status 0 particles
	// this is because SUSYAna does something strange when reading HERWIG info
      case 0:
      case 1:
      case 2:
	break;
	// chargino
      case 1000024:
	nGauginos++;
	if (truth_pdgid[i] > 0)
	  charge++;
	else
	  charge--;	    
	break;
	// neutralino
      case 1000022:
      case 1000023:
	nGauginos++;
	break;
      case 1000013:
	nSleptons++;
	if (truth_pdgid[i] > 0)
	  charge--;
	else
	  charge++;
	break;
      case 1000014:
	nSleptons++;
	break;
      case 24:
	nBosons++;
	if (truth_pdgid[i] > 0)
	  charge++;
	else
	  charge--;
	break;
      default:
	THROW("Unexpected particle found in Signal sample");
    }
  }

  // do some cross-checks
  if ((nLeptons != 1 || nGauginos != 1) && (nSleptons != 1 || nBosons != 1)) {
    THROW("This does not look like signal - expect one muon/neutrino + one gaugino\n"
	  " or one sneutrino/smuon + one boson");
  }

  // follow the decay chain to find out if there are two muons (no neutrinos)
  map<char, vector<int> > particles;
  ReconstructFinalState(particles, 1);
  TLorentzVector slepton;
  for (map<char, vector<int> >::const_iterator it = particles.begin(); 
       it != particles.end(); ++it) {
    char ptype = it->first;
    DEBUG("Type " << ptype << ":");
    for (vector<int>::const_iterator it2 = it->second.begin();
	 it2 != it->second.end(); ++it2) {
      DEBUG("Final state particle: " << *it2 << ", pdgid = " << truth_pdgid[*it2]);
      // sum up four-momenta of all final state particles
      if (ptype == 'l' || ptype == 'n' || ptype == 'q') {
	TLorentzVector particle(truth_px[*it2], truth_py[*it2], truth_pz[*it2], truth_E[*it2]);
	slepton += particle;
      }
    }
  }
  // is it signal (i.e. two muons and at least two jets) or not?
  if (particles['l'].size() == 2 && particles['q'].size() >= 2) {
    fIsSignal = true;
  }
  else {
    fIsSignal = false;
  }

  // take only signal events or background events if requested
  if ((fType == "signal" && !fIsSignal) || (fType == "background" && fIsSignal))
    return;
 
  Fill("Sig_nMuon", particles['l'].size());
  Fill("Sig_nNeutrino", particles['n'].size());
  Fill("Sig_nQuarksSig", particles['q'].size());
  Fill("Sig_nQuarksBack", particles['b'].size()*2);
  Fill("Sig_nQuarks", particles['q'].size()+particles['b'].size()*2);
  Fill("Sig_SleptonMass", slepton.M());
  Fill("Sig_SleptonCharge", charge);

  // get true muon pt distribution, true quark pt distribution etc.
  double ptMax=-1E9;
  double ptMin=1E9;
  for (vector<int>::const_iterator it = particles['l'].begin();
       it != particles['l'].end(); it++) {
    ptMax = TMath::Max(ptMax, truth_pt[*it]);
    ptMin = TMath::Min(ptMin, truth_pt[*it]);
    Fill("Sig_EtaMu", truth_eta[*it]);
    Fill("Sig_PhiMu", truth_phi[*it]);
  }
  Fill("Sig_ptMu1", ptMax);
  Fill("Sig_ptMu2", ptMin);
  for (vector<int>::const_iterator it = particles['q'].begin();
       it != particles['q'].end(); it++) {
    ptMax = TMath::Max(ptMax, truth_pt[*it]);
    ptMin = TMath::Min(ptMin, truth_pt[*it]);
    Fill("Sig_EtaQuark", truth_eta[*it]);
    Fill("Sig_PhiQuark", truth_phi[*it]);
    Fill("Sig_ptQuark1", ptMax);
    Fill("Sig_ptQuark2", ptMin);
  }

  // everything below this line for signal MC only...
  if (!fIsSignal)
    return;

  // compute lepton kinematics
  truth_pt_comparator comp(*this);
  sort(particles['l'].begin(), particles['l'].end(), comp);
  int id1 = particles['l'][0];
  int id2 = particles['l'][1];
  TLorentzVector mu1(truth_px[id1], truth_py[id1], truth_pz[id1], truth_E[id1]);
  TLorentzVector mu2(truth_px[id2], truth_py[id2], truth_pz[id2], truth_E[id2]);
  Fill("Sig_MuMass", (mu1+mu2).M());
  Fill("Sig_MuDPhi", DeltaPhi(truth_phi[id1], truth_phi[id2]));
  Fill("Sig_MuDEta", TMath::Abs(truth_eta[id1]-truth_eta[id2]));
  Fill("Sig_MuDPt", TMath::Abs(truth_pt[id1]-truth_pt[id2]));
  Fill("Sig_MuSumPt", TMath::Abs(truth_pt[id1]+truth_pt[id2]));
  Fill("Sig_MuPt", truth_pt[id1], truth_pt[id2]);

  // compute quark kinematics
  // if more than two quarks, take the most energetic ones
  sort(particles['q'].begin(), particles['q'].end(), comp);
  int id3 = particles['q'][0];
  int id4 = particles['q'][1];
  TLorentzVector q3(truth_px[id3], truth_py[id3], truth_pz[id3], truth_E[id3]);
  TLorentzVector q4(truth_px[id4], truth_py[id4], truth_pz[id4], truth_E[id4]);
  Fill("Sig_JetMass", (q3+q4).M());
  Fill("Sig_JetDPhi", DeltaPhi(truth_phi[id3], truth_phi[id4]));
  Fill("Sig_JetDEta", TMath::Abs(truth_eta[id3]-truth_eta[id4]));
  Fill("Sig_JetDPt", TMath::Abs(truth_pt[id3]-truth_pt[id4]));
  Fill("Sig_JetSumPt", TMath::Abs(truth_pt[id3]+truth_pt[id4]));
  Fill("Sig_JetPt", truth_pt[id3], truth_pt[id4]);

  // look at jets and leptons
  Fill("Sig_Mass3", (mu2+q3+q4).M());
  double DeltaR = TMath::Min(
    TMath::Min(mu1.DeltaR(q3), mu1.DeltaR(q4)),
    TMath::Min(mu2.DeltaR(q3), mu2.DeltaR(q4))
    );
  Fill("Sig_MuIsoR", DeltaR);

  // what happens to the muons when they are not reconstructed?
  if (muo_n < 2) {
    // match stable generated muon to these muons
    for (int i = 0; i < truthl_n; i++) {
      if (abs(truthl_pdgid[i]) != 13) {
	continue;
      }
      TLorentzVector mu_gen(truthl_px[i], truthl_py[i], truthl_pz[i], truthl_E[i]);
      int matched = -1;
      for (int j = 0; j < muo_n; j++) {
	TLorentzVector mu_rec(muo_px[j], muo_py[j], muo_pz[j], muo_E[j]);
	if (mu_gen.DeltaR(mu_rec) < 0.3)
	  matched = j;
      }
      // We are interested in those we did not reconstruct, i.e. cannot match
      if (matched != -1) {
	// just make sure SUSYana likes me :-)
	if (muo_truth[matched] != i) {
	  WARNING("SUSYAna vs findsusyb3 MC truth muon matching difference");
	}
	continue;
      }
      Fill("Sig_LostMuEta", truthl_eta[i]);
      Fill("Sig_LostMuPt", truthl_pt[i]);
      Fill("Sig_LostMuPtEta", truthl_pt[i], truthl_eta[i]);
    }
  }
}

// main event loop
void Analysis::Loop() 
{
  // process and memory info
  ProcInfo_t info;

  // create histograms
  CreateHistograms();

  // set branch addresses in output tree if necessary
  if (fFill) {
    SetBranchAddresses();
  }

  // main event loop
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  INFO("Analysis: Chain contains " << nentries << " events");
  if (fMaxEvents > 0) 
    INFO("Stopping after fMaxEvents = " << fMaxEvents);
  
  for (Long64_t jentry=0; jentry < nentries && (jentry < fMaxEvents || fMaxEvents < 0) ; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) 
      break;
    nb = fChain->GetEntry(jentry);   
    nbytes += nb;

    // some progress / process indicators
    if (fmod((double)jentry,(double)1000)==0) 
      INFO("Event " << jentry);
    if (fmod((double)jentry,(double)10000)==0) {
      gSystem->GetProcInfo(&info);
      INFO("Memory in MB : " << info.fMemResident/1000. << " (resident) " 
	   << info.fMemVirtual/1000. << " (virtual) ");
    }
    if (fDumpTruth) {
      TruthDump();
    }

    // Dump event/object information, see Dump.cxx
    if (fDumpAll) {
      BasicDump(jentry);
      TriggerDump("*"); // select e.g all muon trigger with "Mu"
      TruthDump();
      VertexDump();
      MuonDump(0);      // [less] detailed information 1 [0]
      TruthJetDump();
      CaloJetDump();
      PFJetDump();
      METDump();
      SCDump();
      EleDump(0);       // [less] detailed information 1 [0]
      PFEleDump(0);       // [less] detailed information 1 [0]
    }

    // find duplicates
    if (fFindDuplicates) {
      if (FindDuplicates(global_run, global_event, lumi_section, global_pthat)) {
	// depending on severity one could stop here
	WARNING("Analysis: found duplicate ");
      }
    }

    //////////////////////////////////////////////////////////////////////
    // Split RPV SUSY in RPV signal and RPV background
    if (fType == "signal" || fType == "background") {
      SignalStudy();
      // take only signal events or background events if requested
      if ((fType == "signal" && !fIsSignal) || (fType == "background" && fIsSignal))
	continue;
    }

    FillNoWeight("log_global_weight", TMath::Log(global_weight)/TMath::Log(10.));
    fWeight = 1; /* @todo after signal bugfix: fWeight = global_weight; */

    //////////////////////////////////////////////////////////////////////
    // pileup reweighting
    if (fType == "mc" || fType == "signal" || fType == "background") {
      fWeight *= LumiWeights_.weight3D(pu_num_int[0], pu_num_int[1], pu_num_int[2]);
    }

    //////////////////////////////////////////////////////////////////////
    // Redo skimmer cuts

    // fill histos before cuts
    Fill("bSkim_muo_n", muo_n);
    if (muo_n > 0)
      Fill("bSkim_muo_pt0", muo_pt[0]);
    if (muo_n > 1) 
      Fill("bSkim_muo_pt1", muo_pt[1]);

    // fill N-1 histograms
    TCutList skimCuts;
    skimCuts.Set("muo_n", muo_n >= fSkimMuons);
    skimCuts.Set("muo_pt0", muo_n > 0 ? muo_pt[0] >= fSkimMuoptfirst : false);
    skimCuts.Set("muo_pt1", muo_n > 1 ? muo_pt[1] >= fSkimMuoptother : false);
    if (skimCuts.PassesAllBut("muo_n")) {
      Fill("nSkim_muo_n", muo_n);
    }
    if (skimCuts.PassesAllBut("muo_pt0")) {
      Fill("nSkim_muo_pt0", muo_pt[0]);
    }
    if (skimCuts.PassesAllBut("muo_pt1")) {
      Fill("nSkim_muo_pt1", muo_pt[1]);
    }

    // apply cuts
    if (fSkimActive) {
      if (!skimCuts.PassesAll()) {
	continue; 
      }
    }

    //////////////////////////////////////////////////////////////////////
    // Trigger selection

    bool rejection = true;
    // loop over all triggers
    for (int i = 0; i < trig_n; i++){
      *fTrigger = unpack(trig_name[i]);
      // isolated muon triggers - not used
      // if (fTrigger->Contains("HLT_IsoMu20_v") || fTrigger->Contains("HLT_IsoMu17_v") || fTrigger->Contains("HLT_IsoMu15_v") || fTrigger->Contains("HLT_IsoMu30_eta2p1_v")){
      // double muon triggers
      if (fTrigger->Contains("HLT_Mu17_TkMu8_v") ||
	  fTrigger->Contains("HLT_Mu17_Mu8_v") ||
	  fTrigger->Contains("HLT_Mu13_Mu8_v") ||
	  fTrigger->Contains("HLT_Mu13_Mu7_v") ||
	  fTrigger->Contains("HLT_DoubleMu7") ||
	  fTrigger->Contains("HLT_DoubleMu6")) {
	if (trig_HLTprescale[i] == 1) {
	  rejection = false;
	  break;
	}
      }
    }
    if (rejection)
      continue;

    //////////////////////////////////////////////////////////////////////
    // Object ID

    // Vertices - all cuts should have already been done by the skimmer
    int nVertex = 0;
    TCutList VertexCuts;
    for (int i = 0; i < vtx_n; i++) {
      // do computation
      double x = vtx_x[i] - bs_x;
      double y = vtx_y[i] - bs_y;
      double r = TMath::Sqrt(x*x+y*y);

      // set cut values
      VertexCuts.Set("vtx_fake", vtx_fake[i] == 0);
      VertexCuts.Set("vtx_ndof", vtx_ndof[i] > 4.);
      VertexCuts.Set("vtx_z", vtx_z[i] < 24.);
      VertexCuts.Set("vtx_bs", r < 2.0);

      // Fill N-1 histograms
      if (VertexCuts.PassesAllBut("vtx_fake"))
	Fill("nObj_vtx_fake", vtx_fake[i]);
      if (VertexCuts.PassesAllBut("vtx_ndof"))
	Fill("nObj_vtx_ndof", vtx_ndof[i]);
      if (VertexCuts.PassesAllBut("vtx_z"))
	Fill("nObj_vtx_z", vtx_z[i]);
      if (VertexCuts.PassesAllBut("vtx_bs"))
	Fill("nObj_vtx_bs", r);

      // apply cuts
      if (VertexCuts.PassesAll()) {
	// check if a cut on radial distance to beamspot would help
	nVertex++;
      }
    }

    // Electrons
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SimpleCutBasedEleID2011
    int nElectron = 0;
    TCutList ElectronCuts;
    for (int i = 0; i <pfele_n; i++) {
      // set cuts
      ElectronCuts.Set("pfele_pt", pfele_pt[i] > 15.);

      // fill N-1
      if (ElectronCuts.PassesAllBut("pfele_pt"))
	Fill("nPfele_pt", pfele_pt[i]);

      // apply cuts
      if (ElectronCuts.PassesAll())
	nElectron++;
    }
    
    // Muons
    // see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId, tight muon selection
    vector <int> muons;
    TCutList MuonCuts;
    for (int i = 0; i < muo_n; i++){
      // set cuts
      MuonCuts.Set("muo_pt", muo_pt[i] > 10.);
      MuonCuts.Set("muo_eta", fabs(muo_eta[i]) <= 2.1);
      MuonCuts.Set("muo_id", muo_ID[i][1] == 1);
      MuonCuts.Set("muo_TrkChiNormCm", muo_TrkChiNormCm[i] < 10.);
      MuonCuts.Set("muo_hitsCm", muo_hitsCm[i] > 0);
      MuonCuts.Set("muo_ChambersMatched", muo_ChambersMatched[i] > 1);
      MuonCuts.Set("muo_d0Tk", fabs(muo_d0Tk[i]) < 0.2);
      // MuonCuts.Set("muo_d0bsCm", fabs(muo_d0bsCm[i]) < 0.2);
      MuonCuts.Set("muo_ValidPixelHitsCm", muo_ValidPixelHitsCm[i] > 0);
      // MuonCuts.Set("muo_TrackerLayersMeasCm", muo_TrackerLayersMeasCm[i] > 8);
      MuonCuts.Set("muo_ValidTrackerHitsCm", muo_ValidTrackerHitsCm[i] > 10);

      // N-1 plots
      if (MuonCuts.PassesAllBut("muo_pt")) 
	Fill("nMuon_pt",  muo_pt[i]);
      if (MuonCuts.PassesAllBut("muo_eta")) 
	Fill("nMuon_eta",  muo_eta[i]);
      if (MuonCuts.PassesAllBut("muo_id")) 
	Fill("nMuon_id",  muo_ID[i][1]);
      if (MuonCuts.PassesAllBut("muo_TrkChiNormCm")) 
	Fill("nMuon_TrkChiNormCm",  muo_TrkChiNormCm[i]);
      if (MuonCuts.PassesAllBut("muo_hitsCm")) 
	Fill("nMuon_hitsCm",  muo_hitsCm[i]);
      if (MuonCuts.PassesAllBut("muo_ChambersMatched")) 
	Fill("nMuon_ChambersMatched",  muo_ChambersMatched[i]);
      if (MuonCuts.PassesAllBut("muo_d0Tk")) 
	Fill("nMuon_d0Tk",  fabs(muo_d0Tk[i]));
      if (MuonCuts.PassesAllBut("muo_ValidPixelHitsCm")) 
	Fill("nMuon_ValidPixelHitsCm",  muo_ValidPixelHitsCm[i]);
      if (MuonCuts.PassesAllBut("muo_ValidTrackerHitsCm")) 
	Fill("nMuon_ValidTrackerHitsCm",  muo_ValidTrackerHitsCm[i]);
      
      // apply cuts
      if (MuonCuts.PassesAll()) {
        muons.push_back(i);
      }
    }
    int nMuon = muons.size();

    // Jets
    // https://twiki.cern.ch/twiki/bin/view/CMS/JetID
    double HT = 0;
    vector <int> jets;
    TCutList JetCuts;
    for (int i = 0; i < pfjet_n; i++) {
      // A VBTF muon might be reconstructed as a PF jet. Therefore take only
      // those PF jets that do not "coincide" with the muon
      double Rmin = 1000;
      TLorentzVector jet(pfjet_px[i], pfjet_py[i], pfjet_pz[i], pfjet_E[i]);
      for (int j = 0; j < nMuon; j++) {
	TLorentzVector mu(muo_px[muons[j]], muo_py[muons[j]], muo_pz[muons[j]], 
			  muo_E[muons[j]]);
	double R = mu.DeltaR(jet);
	if (R < Rmin) 
	  Rmin = R;
      }

      // set cuts
      JetCuts.Set("pfjet_pt", pfjet_pt[i] > 30.);
      JetCuts.Set("pfjet_eta", fabs(pfjet_eta[i]) < 2.4);
      JetCuts.Set("pfjet_const", pfjet_const[i] > 1);
      JetCuts.Set("pfjet_pff_0", pfjet_PFF[i][0] > 0.);
      JetCuts.Set("pfjet_pff_1", pfjet_PFF[i][1] < 0.99);
      JetCuts.Set("pfjet_pff_2", pfjet_PFF[i][2] < 0.99);
      JetCuts.Set("pfjet_pff_3", pfjet_PFF[i][3] < 0.99);
      JetCuts.Set("pfjet_pfn_0", pfjet_PFN[i][0] > 0);
      JetCuts.Set("pfjet_muo_deltaR", Rmin > 0.05);

      // N-1 histograms
      if (JetCuts.PassesAllBut("pfjet_pt"))
	Fill("nPfJet_pt", pfjet_pt[i]);
      if (JetCuts.PassesAllBut("pfjet_eta"))
	Fill("nPfJet_eta", pfjet_eta[i]);
      if (JetCuts.PassesAllBut("pfjet_const"))
	Fill("nPfJet_const", pfjet_const[i]);
      if (JetCuts.PassesAllBut("pfjet_pff_0"))
	Fill("nPfJet_pff_0", pfjet_PFF[i][0]);
      if (JetCuts.PassesAllBut("pfjet_pff_1"))
	Fill("nPfJet_pff_1", pfjet_PFF[i][1]);
      if (JetCuts.PassesAllBut("pfjet_pff_2"))
	Fill("nPfJet_pff_2", pfjet_PFF[i][2]);
      if (JetCuts.PassesAllBut("pfjet_pff_3"))
	Fill("nPfJet_pff_3", pfjet_PFF[i][3]);
      if (JetCuts.PassesAllBut("pfjet_pfn_0"))
	Fill("nPfJet_pfn_0", pfjet_PFN[i][0]);
      if (JetCuts.PassesAllBut("pfjet_muo_deltaR"))
	Fill("nPfJet_muo_deltaR", Rmin);
					  
      // apply cuts
      if (JetCuts.PassesAll()) {
	jets.push_back(i);
	HT += pfjet_Et[i];
      }
    }
    int nJet = jets.size();

    //////////////////////////////////////////////////////////////////////
    // Event cleaning

    // require good vertex, but no electron in event, and no ecal noise
    if (nVertex < 1 || nElectron > 0 || noise_ecal_r9 > 0.9) {
        continue;
    }

    //////////////////////////////////////////////////////////////////////
    // @todo: Trigger matching

    //////////////////////////////////////////////////////////////////////
    // loose muon id
    if (nMuon < 2 || muo_pt[muons[0]] <= 20. || muo_pt[muons[1]] <= 15.) {
      continue;
    }

    //////////////////////////////////////////////////////////////////////
    // loose jet id
    if (nJet < 2 || pfjet_pt[jets[0]] <= 30. || pfjet_pt[jets[1]] <= 30.) {
      continue;
    }

    //////////////////////////////////////////////////////////////////////
    // tight muon id
    vector <int> tight_muons;
    for (int i = 0; i < nMuon; i++) {
      double Rmin = 1000.;
      TLorentzVector muon(muo_px[muons[i]], muo_py[muons[i]], muo_pz[muons[i]], 
			  muo_E[muons[i]]);
      for (int j = 0; j < nJet; j++) {
	TLorentzVector jet(pfjet_px[j], pfjet_py[j], pfjet_pz[j], pfjet_E[j]);
        double R = jet.DeltaR(muon);
        if (R < Rmin) 
	  Rmin = R;
      }
      // muon isolation cuts
      if (muo_TrkIso[muons[i]] <= 2. &&
          muo_ECalIso[muons[i]] <= 2. &&
          muo_HCalIso[muons[i]] <= 2. &&
          Rmin >= 0.4) {
        tight_muons.push_back(muons[i]);
      }
    }

    Int_t nMuonsTight = tight_muons.size();
    if (nMuonsTight != 2 || 
	muo_pt[tight_muons[0]] <= 20. || 
	muo_pt[tight_muons[1]] <= 15. || 
	fabs(muo_dzbsCm[tight_muons[0]]-muo_dzbsCm[tight_muons[1]]) > 0.08) {
      continue;
    }
    TLorentzVector mu0(muo_px[0], muo_py[0], muo_pz[0], muo_E[0]);
    TLorentzVector mu1(muo_px[1], muo_py[1], muo_pz[1], muo_E[1]);
    
    /// Particle flow MET
    if (met_et[3] >= 50.) {
      continue;
    }

    double m_mumu = (mu0+mu1).M();
    /// M(mu, mu)
    if (m_mumu < 50.) {
      continue;
    }

    Fill("m_mumu_cut", m_mumu);
    Fill("muo_n_cut", muo_n);

    Fill("vtx_n", vtx_n);
    for (int i = 0; i < vtx_n; i++) {
      Fill("vtx_ntr", vtx_ntr[i]);
      Fill("vtx_ndof", vtx_ndof[i]);
      Fill("vtx_x", vtx_x[i]);
      Fill("vtx_y", vtx_y[i]);
      Fill("vtx_z", vtx_z[i]);
    }

    /// muon charge
    if ((muo_charge[tight_muons[0]]*muo_charge[tight_muons[1]]) == -1) {
      continue;
    }
    
    Fill("m_mumu", m_mumu);
    Fill("muo_n", muo_n);
    FillNoWeight("weight", fWeight);
    FillNoWeight("log_weight", TMath::Log(fWeight)/TMath::Log(10.));

    //////////////////////////////////////////////////////////////////////
    // fill output tree
    if (fFill)
      fOutputTree.Fill();
  }
}

void Analysis::CreateHistograms()
{
  //////////////////////////////////////////////////////////////////////
  // Histograms
  INFO("Analysis::CreateHistograms()");
  // use weighted histograms
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  // change to the file we want to have the histograms in
  TFile * outFile = fOutputTree.GetCurrentFile();
  outFile->cd();

  // signal histograms
  CreateHisto("Sig_nMuon", "Signal MC truth: number of muons in final state", 10, -0.5, 9.5);
  CreateHisto("Sig_nNeutrino", "Signal MC truth, number of neutrinos in final state", 10, -0.5, 9.5);
  CreateHisto("Sig_nQuarksSig", "Signal MC truth, number of quarks (from sparticles) in final state", 10, -0.5, 9.5);
  CreateHisto("Sig_nQuarksBack", "Signal MC truth, number of quarks from decay chain in final state", 10, -0.5, 9.5);
  CreateHisto("Sig_nQuarks", "Signal MC truth: number of quarks in final states", 10, -0.5, 9.5);
  CreateHisto("Sig_SleptonMass", "Resonance mass of slepton@GeV", 1000, 0, 1000);
  CreateHisto("Sig_SleptonCharge", "Charge of slepton in units of e", 3, -1.5, 1.5);

  CreateHisto("Sig_ptMu1", "p_{T} of leading muon@GeV", 500, 0, 500);
  CreateHisto("Sig_ptMu2", "p_{T} of second leading muon@GeV", 500, 0, 500);
  CreateHisto("Sig_EtaMu", "#eta_{#mu}", 100, -5, 5);
  CreateHisto("Sig_PhiMu", "#phi_{#mu}@rad", 628, -3.14, 3.14);
  CreateHisto("Sig_ptQuark1", "p_{T} of leading quark@GeV", 500, 0, 500);
  CreateHisto("Sig_ptQuark2", "p_{T} of second leading quark@GeV", 500, 0, 500);
  CreateHisto("Sig_EtaQuark", "#eta_{#mu}", 100, -5, 5);
  CreateHisto("Sig_PhiQuark", "#phi_{#mu}@rad", 628, -3.14, 3.14);

  // combine muons from decay chain
  CreateHisto("Sig_MuMass", "m(#mu_{1},#mu_{2})@GeV", 500, 0, 500);
  CreateHisto("Sig_MuDPhi", "#Delta#phi(#mu_{1},#mu_{2})@rad", 314, 0., 3.14);
  CreateHisto("Sig_MuDEta", "#Delta#eta(#mu_{1},#mu_{2})", 100, 0., 10.);
  CreateHisto("Sig_MuDPt",  "#Deltap_{T}(#mu_{1},#mu_{2})@GeV", 500, 0, 500);
  CreateHisto("Sig_MuSumPt", "p_{T}(#mu_{1}) + p_{T}(#mu_{2})@GeV", 500, 0, 500);
  CreateHisto("Sig_MuPt", "p_{T}(#mu_{1}):p_{T}(#mu_{2})@GeV", 50, 0, 100, 50, 0, 100);

  // combine quarks from decay chain
  CreateHisto("Sig_JetMass", "m(j_{1},j_{2})@GeV", 500, 0, 500);
  CreateHisto("Sig_JetDPhi", "#Delta#phi(j_{1},j_{2})@rad", 314, 0., 3.14);
  CreateHisto("Sig_JetDEta", "#Delta#eta(j_{1},j_{2})", 100, 0., 10.);
  CreateHisto("Sig_JetDPt",  "#Deltap_{T}(j_{1},j_{2})@GeV", 500, 0, 500);
  CreateHisto("Sig_JetSumPt", "p_{T}(jet_{1}) + p_{T}(jet_{2})@GeV", 500, 0, 500);
  CreateHisto("Sig_JetPt", "p_{T}(q_{1}):p_{T}(q_{2})@GeV", 50, 0, 100, 50, 0, 100);

  // combine leptons and quarks
  CreateHisto("Sig_Mass3", "m(#mu_{2}, q_{1}, q_{2})@GeV", 500, 0, 500);
  CreateHisto("Sig_MuIsoR", "smallest #DeltaR between any jet and any muon", 100, 0, 5);

  // only for signal
  CreateHisto("Sig_LostMuEta", "Lost muons #eta_{#mu}", 100, -5, 5);
  CreateHisto("Sig_LostMuPt", "Lost muons p_{T}", 500, 0, 500);
  CreateHisto("Sig_LostMuPtEta", "Lost muons p_{T}@GeV:Lost muons #eta_{#mu}", 
	      50, 0, 100, 50, -5, 5);

  // Reskimming
  CreateHisto("bSkim_muo_n", "Number of muons", 10, -0.5, 9.5);
  CreateHisto("bSkim_muo_pt0", "pt of leading muon", 500, 0, 500);
  CreateHisto("bSkim_muo_pt1", "pt of second muon", 500, 0, 500);
  CreateHisto("nSkim_muo_n", "Number of muons", 10, -0.5, 9.5);
  CreateHisto("nSkim_muo_pt0", "pt of leading muon", 500, 0, 500);
  CreateHisto("nSkim_muo_pt1", "pt of second muon (skimmer cuts)", 500, 0, 500);

  // Pileup reweighting, event cleaning
  CreateHisto("vtx_n", "Number of vertices in event", 50, -0.5, 49.5);
  CreateHisto("vtx_ntr", "Number tracks in vertex", 151, -0.5, 150.5);
  CreateHisto("vtx_ndof", "Number of deegrees of freedom of vertex", 151, -0.5, 150.5);
  CreateHisto("vtx_x", "X position of vertex@cm", 100, -2, 2);
  CreateHisto("vtx_y", "Y position of vertex@cm", 100, -2, 2);
  CreateHisto("vtx_z", "Z position of vertex@cm", 120, -30, 30);

  // Object selection: vertex
  CreateHisto("nObj_vtx_fake", "vertex isFake()", 2, 0., 2.);
  CreateHisto("nObj_vtx_ndof", "Number of deegrees of freedom of vertex", 151, -0.5, 150.5);
  CreateHisto("nObj_vtx_z", "Z position of vertex@cm", 120, -30, 30);
  CreateHisto("nObj_vtx_bs", "|r_{bs} - r_{vtx}|@cm", 100, 0, 10);

  // Object selection: electron
  CreateHisto("nPfele_pt", "particle flow electron p_{T}@GeV", 1000, 0, 1000);

  // Object selection: muon
  CreateHisto("nMuon_pt", "#mu p_{T}@GeV", 1000, 0, 1000);
  CreateHisto("nMuon_eta","#eta_{#mu}", 100, -3, 3);
  CreateHisto("nMuon_id", "global muon flag", 2, 0, 2);
  CreateHisto("nMuon_TrkChiNormCm", "#mu #chi^{2}/ndof", 100, 0, 100);
  CreateHisto("nMuon_hitsCm", "combined muon hits", 20, 0, 20);
  CreateHisto("nMuon_ChambersMatched", "matched muon chambers", 20, 0, 20);
  CreateHisto("nMuon_d0Tk", "#mu impact parameter wrt vertex from inner track@cm", 100, 0, 10);
  CreateHisto("nMuon_ValidPixelHitsCm", "#mu pixel hits", 5, 0, 5);
  CreateHisto("nMuon_ValidTrackerHitsCm", "#mu tracker hits", 25, 0, 25);

  // Object selection: jets
  CreateHisto("nPfJet_pt", "jet p_{T}", 1000, 0, 1000);
  CreateHisto("nPfJet_eta", "#eta_{jet}", 100, -3, 3);
  CreateHisto("nPfJet_const", "number of jet constituents", 80, -0.5, 79.5);
  CreateHisto("nPfJet_pff_0", "jet energy fraction: charged hadrons", 100, 0, 1);
  CreateHisto("nPfJet_pff_1", "jet energy fraction: neutral hadrons", 100, 0, 1);
  CreateHisto("nPfJet_pff_2", "jet energy fraction: neutral EM hadrons", 100, 0, 1);
  CreateHisto("nPfJet_pff_3", "jet energy fraction: charged EM hadrons", 100, 0, 1);
  CreateHisto("nPfJet_pfn_0", "jet energy constitutens: charged hadrons", 40, -0.5, 39.5);
  CreateHisto("nPfJet_muo_deltaR", "min #Delta R(pf jet, muon)", 100, 0, 5);

  // create individual histograms
  CreateHisto("noise_ecal_r9", "noise_ecal_r9", 100, 0, 1);
  CreateHisto("m_mumu_cut", "m(#mu^{+}, #mu^{-})@GeV", 500, 0, 500);
  CreateHisto("muo_n_cut", "Number of muons", 30, 0, 30);
  CreateHisto("m_mumu", "m(#mu^{+}, #mu^{-})@GeV", 500, 0, 500);
  CreateHisto("muo_n", "Number of muons", 30, 0, 30);

  // event weights
  CreateHisto("log_global_weight", "log(global_weight)", 100, -5, 5);
  CreateHisto("log_weight", "log(weight)", 100, -5, 5);
  CreateHisto("weight", "weight", 100, 0., 2.);
}

// duplicates finder
bool Analysis::FindDuplicates(int run, int evt, double x1, double x2) {

  pair<int,int> temp1(run,evt);
  pair<double,double> temp2(x1,x2);
  Key key (temp1, temp2);
  
  KeyIter pos = _keys.find (key);
  
  if (pos == _keys.end()) {
    _keys.insert (key);
    return false;
  }
  else {
    ERROR("Analysis: duplicate run " << run << " , evt " << evt 
	  << ", x1 " << x1 << ", x2 " << x2);
    return true;
  }
}

void Analysis::SetBranchAddresses()
{
  //////////////////////////////////////////////////////////////////////
  // initialize output tree branches
  INFO("Analysis::init(): Setting output branch addresses");
  TObjArray * branchlist = fInputTree.GetListOfBranches();
  TIter next(branchlist);
  while (TBranch * branch = (TBranch *) next()) {
    if (fOutputTree.GetBranch(branch->GetName())) {
      fOutputTree.SetBranchAddress(branch->GetName(), branch->GetAddress());
      LOG(4, "branch " << branch->GetName() << " has address " << (void *) branch->GetAddress());
    }
  }
  fOutputTree.SetMaxTreeSize(fMaxTreeSize);
}

void Analysis::CreateHisto(const char * name, const char * title, Int_t nbinsx, Double_t xlow, Double_t xup)
{
  histo[name] = new TH1D(Form("h1_%s", name), title, nbinsx, xlow, xup);
}

void Analysis::CreateHisto(const char * name, const char * title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup)
{
  histo2[name] = new TH2D(Form("h2_%s", name), title, nbinsx, xlow, xup, nbinsy, ylow, yup);
}

void Analysis::Fill(const char * name, double value)
{
  TH1D * h = histo[name];
  if (h != 0)
    h->Fill(value, fWeight);
  else {
    THROW(std::string("Histogram \"") + name + std::string("\" not existing. Did you misspell or forgot to create?"));
  }
}

void Analysis::FillNoWeight(const char * name, double value)
{
  TH1D * h = histo[name];
  if (h != 0)
    h->Fill(value);
  else {
    THROW(std::string("Histogram \"") + name + std::string("\" not existing. Did you misspell or forgot to create?"));
  }
}

void Analysis::Fill(const char * name, double x, double y)
{
  TH2D * h = histo2[name];
  if (h != 0)
    histo2[name]->Fill(x, y, fWeight);
  else {
    THROW(std::string("Histogram \"") + name + std::string("\" not existing. Did you misspell or forgot to create?"));
  }
}

// helper
double Analysis::M2(double E1, double px1, double py1, double pz1, double E2, double px2, double py2, double pz2) {
  double InvMass2 =   (E1 + E2)*(E1 + E2)
    - (px1 + px2)*(px1 + px2)
    - (py1 + py2)*(py1 + py2)
    - (pz1 + pz2)*(pz1 + pz2);
  if (InvMass2 < 0) cout << "Mass Square negative: InvariantMass" << InvMass2 << endl; 
  return std::sqrt(InvMass2);
}

double Analysis::M2v2(double pt1, double eta1, double phi1, double m1, double pt2, double eta2, double phi2, double m2){
  TLorentzVector lorentz1;
  TLorentzVector lorentz2;
  lorentz1.SetPtEtaPhiM(pt1, eta1, phi1, m1);
  lorentz2.SetPtEtaPhiM(pt2, eta2, phi2, m2);
  TLorentzVector lorentz1plus2 = lorentz1 + lorentz2;
  return lorentz1plus2.M();
}

double Analysis::AngleMuons(double pt1, double eta1, double phi1, double m1, double pt2, double eta2, double phi2, double m2){
  TLorentzVector * lorentz1 = new TLorentzVector();
  TLorentzVector * lorentz2 = new TLorentzVector();
  lorentz1->SetPtEtaPhiM(pt1, eta1, phi1, m1);
  lorentz2->SetPtEtaPhiM(pt2, eta2, phi2, m2);
   
  double part1_px = lorentz1->Px();
  double part1_py = lorentz1->Py();
  double part1_pz = lorentz1->Pz();
  double part1_p = lorentz1->P();

  double part2_px = lorentz2->Px();
  double part2_py = lorentz2->Py();
  double part2_pz = lorentz2->Pz();
  double part2_p = lorentz2->P();

  //see exotica muon twiki for documentation
  double angle = acos(((-1.*part1_px * part2_px)+(-1.*part1_py * part2_py)+(-1.*part1_pz * part2_pz))/(part1_p * part2_p));
  return angle;

}

Double_t Analysis::DeltaPhi(double a, double b) {

  double temp = fabs(a-b);
  if (temp <= TMath::Pi())
    return temp;
  else
    return  2.*TMath::Pi() - temp;
}

Double_t Analysis::mT(double et1, double phi1, double et2, double phi2) {

  double mm = 2 * et1 * et2 * ( 1. - cos(phi1 - phi2) );
  return sqrt(mm);

}

Double_t Analysis::MinDEt(const std::vector<TLorentzVector> & objects, 
			 std::vector<UInt_t> * lista, 
			 std::vector<UInt_t> * listb) {
  
  // Find the combination with the lowest DEt
  UInt_t n = objects.size();
  if (n==0) return 0.;
  if (n==1) return objects[0].Et();
  if (n>10)  {
    cout << "MinDEt: n too big : " << n << endl;
    return -1;
  }
  if (lista!=0 && listb!=0) { lista->clear(); listb->clear(); }
  
  double mindiff = 1000000000., diff = 0.;
  
  // Determine the combination that minimises the difference
  std::vector< std::vector<UInt_t> > combinationset1;
  std::vector< std::vector<UInt_t> > combinationset2;
  Combinations::mycombinations(n, combinationset1, combinationset2);
  
  if (combinationset1.size() != combinationset2.size() ) {
    cout << "MinDEt: Combination set sizes to not match - something has gone wrong..." << endl;
  }
  
  for (UInt_t set = 0; set<combinationset1.size(); set++) {
    
    std::vector<UInt_t> la = combinationset1[set]; //!< Temporary list a for calculating best combo
    std::vector<UInt_t> lb = combinationset2[set]; //!< Temporary list b for calculating best combo
    
    Double_t aEt = 0., bEt = 0.;
    for (std::vector<UInt_t>::iterator ia=la.begin();ia!=la.end();++ia) {
      //      cout << (*ia) << " ";
      aEt += objects[ (*ia) ].Et();
    }
    //    cout << ", ";
    for (std::vector<UInt_t>::iterator ib=lb.begin();ib!=lb.end();++ib) {
      bEt += objects[ (*ib) ].Et();
      //      cout << (*ib) << " ";
    }
    //    cout << endl;
    diff = fabs(aEt - bEt);
    //    cout << "Difference in Et is " << diff << endl;
    if (diff < mindiff) {
      mindiff = diff;
      if (lista!=0 && listb!=0) { *lista = la; *listb = lb; }
    }
    la.clear(); lb.clear();
    
  } // end of loop over combination sets
  
  //  cout << "Minimum difference is " << mindiff << endl << endl << endl;
  //  cout << "===========================================" << endl;
  
  return mindiff;  
}

Double_t Analysis::AlphaT(const std::vector<TLorentzVector> & objects) {

  return 0.5*((SumET(objects) - MinDEt(objects))/(MT(objects)));

}

Double_t Analysis::SumET(const std::vector<TLorentzVector> & objects) {

  Double_t sEt = 0;    
  for (std::vector<TLorentzVector>::const_iterator o=objects.begin();o!=objects.end();++o) { 
    sEt+=o->Et(); 
  }
  return sEt;

}

Double_t Analysis::MT(const std::vector<TLorentzVector> & objects) {

  Double_t sEt = 0, sPx = 0, sPy = 0;
  for (std::vector<TLorentzVector>::const_iterator o=objects.begin();o!=objects.end();++o) { 
    sEt+=o->Et();sPx+=o->Px();sPy+=o->Py(); 
  }
  Double_t MTsq = sEt*sEt - sPx*sPx - sPy*sPy;
  return MTsq >= 0. ? sqrt(MTsq) : -sqrt(-MTsq);

}

// JES variation: recalculate Jets and propagate into MET
void Analysis::doJESandRecalculateMET(TString corr) {

  // WARNING : Not propagated into TCMET since this need track corrected Jets
  Double_t caloscale;
  Double_t pfscale;

  Double_t dMEx;
  Double_t dMEy;
  Double_t dSumEt;

  // JME-10-010-pas-v8.pdf
  if (corr == "none")
    return;
  else if (corr == "up") {
    caloscale = 1.04;
    pfscale   = 1.04;
  }  
  else if (corr == "down") {
    caloscale = 0.96;
    pfscale   = 0.96;
  }
  else {
    cout << "Analysis::doJESandRecalculateMET: Option '" << corr << "' not known." << endl;
    cout << "                                 No correction applied." << endl;
    return;
  }
  LOG(4, "Analysis::doJESandRecalculateMET called with option '" << corr << "'");
  
  dMEx   = 0.;
  dMEy   = 0.;
  dSumEt = 0.;
  
  // Calo Jets
  for (int i=0; i<calojet_n; i++) {

    // not final, best guess
    if (calojet_pt_raw[i]>20. && calojet_fem[i]<0.9) {
      dMEx   += (caloscale - 1.) * calojet_px[i];
      dMEy   += (caloscale - 1.) * calojet_py[i];
      dSumEt += (caloscale - 1.) * calojet_Et[i];
    }
    calojet_E[i]  *= caloscale;
    calojet_Et[i] *= caloscale;  
    calojet_p[i]  *= caloscale;   
    calojet_pt[i] *= caloscale;  
    calojet_px[i] *= caloscale;
    calojet_py[i] *= caloscale;
    calojet_pz[i] *= caloscale;
  }

  met_sumet[0] += dSumEt;
  met_ex[0]    -= dMEx;
  met_ey[0]    -= dMEy;

  TVector3 mcalo(met_ex[0], met_ey[0], 0.);

  met_phi[0] = mcalo.Phi();
  met_et[0]  = mcalo.Perp();

  // PF Jets
  dMEx   = 0.;
  dMEy   = 0.;
  dSumEt = 0.;
  
  for (int i=0; i<pfjet_n; i++) {

    // cleaned PF jets, have ptraw > 10 GeV
    dMEx   += (pfscale - 1.) * pfjet_px[i];
    dMEy   += (pfscale - 1.) * pfjet_py[i];
    dSumEt += (pfscale - 1.) * pfjet_Et[i];
    
    pfjet_E[i]  *= pfscale;
    pfjet_Et[i] *= pfscale;  
    pfjet_p[i]  *= pfscale;   
    pfjet_pt[i] *= pfscale;  
    pfjet_px[i] *= pfscale;
    pfjet_py[i] *= pfscale;
    pfjet_pz[i] *= pfscale;
  }
  
  met_sumet[3] += dSumEt;
  met_ex[3]    -= dMEx;
  met_ey[3]    -= dMEy;
  
  TVector3 mpf(met_ex[3], met_ey[3], 0.);

  met_phi[3] = mpf.Phi();
  met_et[3]  = mpf.Perp();

}
