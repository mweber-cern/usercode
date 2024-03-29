#include "Analysis.h"

// ROOT includes
#include <TH2.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TRandom.h>

// C / C++ includes
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <typeinfo>

// local includes
#include "Combinations.h"
#include "Utilities.h"
#include "TCutList.h"
#include "RunLumiRanges.h"
#include "EventFilter.h"

using namespace std;

//////////////////////////////////////////////////////////////////////
// constants

// use type 1 corrected PFMET
const int MET_INDEX = 1;

// externally generated with ../bin/gencutflow.py
extern TH1D * get_cutflow_histogram();

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
    fCfgFile(cfgFile), fBTagging(35535)
{
  //////////////////////////////////////////////////////////////////////
  // Configuration options

  // fill output tree?
  fFill = fCfgFile.GetValue("FillTree", true);

  // Which type of input (signal, background, data, MC)?
  fSample = fCfgFile.GetValue("Sample", "None");
  if (fSample == "None")
    THROW("Value \"Sample\" either \"None\" or not configured in configuration file");
  if (!fSample.compare(0, 4, "data")) {
    fInputType = "data";
  }
  else if (!fSample.compare(0, 6, "signal")) {
    fInputType = "signal";
  }
  else if (!fSample.compare(0, 10, "background")) {
    fInputType = "background";
  }
  else {
    fInputType = "mc";
  }
  INFO("File type: " << fInputType);

  // analysis type
  fAnalysisType = fCfgFile.GetValue("AnalysisType", "default");
  if (fAnalysisType != "default"
      && fAnalysisType != "singlefake" && fAnalysisType != "doublefake") {
    THROW("AnalysisType \"" + fAnalysisType + "\" is not an allowed value");
  }

  split(cfgFile.GetValue("Trigger", "None"), ',', fTrigger);
  fForceUnprescaledTriggers = fCfgFile.GetValue("ForceUnprescaledTriggers", true);
  // maximum number of events to be processed (-1: process all)
  fMaxEvents = fCfgFile.GetValue("MaxEvents", -1);
  // maximum tree size for output file
  fMaxTreeSize = fCfgFile.GetValue("MaxTreeSize", 100000000);
  INFO("Setting maximum tree size to " << fMaxTreeSize);
  // find duplicates?
  fFindDuplicates = fCfgFile.GetValue("FindDuplicates", true);

  // Configure debugging output
  fDumpAll = fCfgFile.GetValue("DumpAll", false);
  fDumpBasic = fCfgFile.GetValue("DumpBasic", false);
  fDumpTruth = fCfgFile.GetValue("DumpTruth", false);
  fDumpTrigger = fCfgFile.GetValue("DumpTrigger", "None");

  fSkimActive = fCfgFile.GetValue("SkimActive", true);
  fSkimMuons = fCfgFile.GetValue("SkimMuons", 2);
  fSkimMuoptfirst = fCfgFile.GetValue("SkimMuoptfirst", 15.);
  fSkimMuoptother = fCfgFile.GetValue("SkimMuoptother", 7.);

  // cuts for T/L ratio
  fLooseMuonRelIso = fCfgFile.GetValue("LooseMuonRelIso", 0.4);
  fTL_met_max = fCfgFile.GetValue("TL_met_max",  50.);
  fTL_ht_min = fCfgFile.GetValue("TL_ht_min",  50.);
  fTL_jetpt_min = fCfgFile.GetValue("TL_jetpt_min",  40.);
  fTL_nloose_min = fCfgFile.GetValue("TL_nloose_min",  0.5);
  fTL_nloose_max = fCfgFile.GetValue("TL_nloose_max",  1.5);
  fTL_mumudz_min = fCfgFile.GetValue("TL_mumudz_min",  0.01);
  fTL_zmass_min = fCfgFile.GetValue("TL_zmass_min",  71.);
  fTL_zmass_max = fCfgFile.GetValue("TL_zmass_max",  111.);
  fTL_mupt_min = fCfgFile.GetValue("TL_mupt_min",  15.);
  fTL_firstmupt_min = fCfgFile.GetValue("TL_firstmupt_min",  20.);
  fTL_jetdphi_min = fCfgFile.GetValue("TL_jetdphi_min",  1.0);
  fTL_mumudphi_max = fCfgFile.GetValue("TL_mumudphi_max",  1.0);
  fTL_mt_max = fCfgFile.GetValue("TL_mt_max",  40.);
  fTL_njets_min = fCfgFile.GetValue("TL_njets_min", 1.);

  INFO("fTL_met_max: " << fTL_met_max);
  INFO("fTL_ht_min: " << fTL_ht_min);
  INFO("fTL_jetpt_min: " << fTL_jetpt_min);
  INFO("fTL_nloose_min: " << fTL_nloose_min);
  INFO("fTL_nloose_max: " << fTL_nloose_max);
  INFO("fTL_mumudz_min: " << fTL_mumudz_min);
  INFO("fTL_zmass_min: " << fTL_zmass_min);
  INFO("fTL_zmass_max: " << fTL_zmass_max);
  INFO("fTL_mupt_min: " << fTL_mupt_min);
  INFO("fTL_firstmupt_min: " << fTL_firstmupt_min);
  INFO("fTL_jetdphi_min: " << fTL_jetdphi_min);
  INFO("fTL_mumudphi_max: " << fTL_mumudphi_max);
  INFO("fTL_mt_max: " << fTL_mt_max);
  INFO("fTL_nJets: " << fTL_njets_min);

  // values for smearing jet energies to obtain correct jet energy resolution
  // (JER)
  fJER_calculation = fCfgFile.GetValue("JER_calculation", false);
  fJER_center = fCfgFile.GetValue("JER_center", 0.);
  fJER_smear = fCfgFile.GetValue("JER_smear", 20.);

  INFO("fJER_calculation = " << fJER_calculation);
  INFO("fJER_center = " << fJER_center);
  INFO("fJER_smear = " << fJER_smear);

  //////////////////////////////////////////////////////////////////////
  // analysis variables
  fIsSignal = false;

  //////////////////////////////////////////////////////////////////////
  // initialize pileup reweighting
  fPileupReweighting = fCfgFile.GetValue("PileupReweighting", true);
  INFO("fPileupReweighting = " << fPileupReweighting);
  const char * PileupDataFile = gSystem->ExpandPathName(fCfgFile.GetValue("PileupDataFile_2012", ""));
  const char * PileupMCFile = gSystem->ExpandPathName(fCfgFile.GetValue("PileupMCFile_2012", ""));
  if (!check_root_file(PileupMCFile)) {
    ERROR("Could not find reweighting file " << PileupMCFile);
    THROW("Missing required file: " + string(PileupMCFile));
  }
  if (!check_root_file(PileupDataFile)) {
    ERROR("Could not find reweighting file " << PileupDataFile);
    THROW("Missing required file: " + string(PileupDataFile));
  }
  LumiWeights_ = reweight::LumiReWeighting(PileupMCFile, PileupDataFile, "pileup","pileup");

  delete PileupDataFile;
  delete PileupMCFile;
  INFO("Initialized pileup reweighting");

  //////////////////////////////////////////////////////////////////////
  // FakeRate
  fFakeRateMethod = fCfgFile.GetValue("FakeRateMethod", "lastbin");
  if (fFakeRateMethod != "lastbin" && fFakeRateMethod != "zero") {
    THROW("Wrong fake rate method \"" + fFakeRateMethod + "\", use \"lastbin\" or \"zero\"");
  }
  fFakeRateDimensions = fCfgFile.GetValue("FakeRateDimensions", 2);
  if (fFakeRateDimensions != 2) { // && fFakeRateDimensions != 3) {
    THROW("Wrong number of fake rate dimensions given: " + std::string(fCfgFile.GetValue("FakeRateDimensions", "")));
  }
  if (fAnalysisType == "singlefake" || fAnalysisType == "doublefake") {
    const char * fakeRateFileName = gSystem->ExpandPathName(fCfgFile.GetValue("FakeRateFile", ""));
    // fFakeRateHisto3D = (TH3D *) get_object(fakeRateFileName, "hRatio3");
    // if (fFakeRateHisto3D == 0) {
    //   THROW("Could not get tight/loose 3D histogram from file");
    // }
    fFakeRateHisto2D = (TH2D *) get_object(fakeRateFileName, "hRatio2");
    if (fFakeRateHisto2D == 0) {
      THROW("Could not get tight/loose 2D histogram from file");
    }
    delete fakeRateFileName;
  }
  fFakeRate[0] = 0.;
  fFakeRate[1] = 0.;
  fSingleFakeWeight = 0.;
  INFO("Initialized fake rate by reading histograms from file");

  // Lumi ranges. Default: empty string, i.e. do not filter
  if (fInputType == "data") {
    fLumiRanges = fCfgFile.GetValue("LumiRanges", "");
    if (fLumiRanges != "") {
      INFO("Loaded golden luminosity ranges from file " << fLumiRanges);
    }
  }
  else {
    fLumiRanges = "";
  }
  
  // Event filter. Default: empty string,i.e. do not filter
  if (fInputType == "data") {
    fEventFilter = fCfgFile.GetValue("EventFilter", "");
    if (fEventFilter != "") {
      INFO("Loaded HCAL laser filter information from file " << fEventFilter);
    }
  }
  else {
    fEventFilter = "";
  }

  fDeltaPhiMin = fCfgFile.GetValue("DeltaPhiMin", 2.0);
  INFO("DeltaPhiMin = " << fDeltaPhiMin);
}

Analysis::~Analysis()
{
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
    if (pdgid == 11) {
      DEBUG("lepton j = " << j << ", pdgid = " << pdgid);
      particles['e'].push_back(j);
      particles['l'].push_back(j);
    }
    else if (pdgid == 13) {
      DEBUG("lepton j = " << j << ", pdgid = " << pdgid);
      particles['m'].push_back(j);
      particles['l'].push_back(j);
    }
    else if (pdgid == 15) {
      DEBUG("lepton j = " << j << ", pdgid = " << pdgid);
      particles['t'].push_back(j);
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
    // + photon
    else if (pdgid >= 21 && pdgid <= 25) {
      DEBUG("boson j = " << j << ", pdgid = " << pdgid);
      particles['b'].push_back(j);
      ReconstructFinalState(particles, truth_evtxid[j]);
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

void Analysis::GetSimplifiedModelParticles(map<string, vector<int> > & particles, int vertex)
{
  //DEBUG("Analysis::CheckForSimplifiedModel(" << vertex << ")");
  //loop over all particles
  for (int i = 0; i < truth_n; i++) {
    // only take particles from this vertex
    if (truth_bvtxid[i] != vertex)
      continue;
    
    int pdgid = abs(truth_pdgid[i]);
    
    if (pdgid == 13) {
      particles["mu"].push_back(i);
    }
    else if (pdgid == 1000024) {
      particles["x+-1"].push_back(i);
      GetSimplifiedModelParticles(particles, truth_evtxid[i]);      
    }
    else if (pdgid == 1000037) {
      particles["x+-2"].push_back(i);
      GetSimplifiedModelParticles(particles, truth_evtxid[i]);      
    }
    else if (pdgid == 1000022) {
      particles["x01"].push_back(i);
      GetSimplifiedModelParticles(particles, truth_evtxid[i]);
    }
    else if (pdgid == 1000023) {
      particles["x02"].push_back(i);
      GetSimplifiedModelParticles(particles, truth_evtxid[i]);
    }
    else if (pdgid == 1000025) {
      particles["x03"].push_back(i);
      GetSimplifiedModelParticles(particles, truth_evtxid[i]);
    }
    else if (pdgid == 1000035) {
      particles["x04"].push_back(i);
      GetSimplifiedModelParticles(particles, truth_evtxid[i]);
    }
    else if (pdgid == 1) {
      if (vertex > 1) {
	particles["q"].push_back(i);
      }
    }
    else if (pdgid == 2) {
      if (vertex > 1) {
	particles["q"].push_back(i);
      }
    }
    else if (pdgid == 24) {
	particles["W"].push_back(i);
	GetSimplifiedModelParticles(particles, truth_evtxid[i]);
    }
    else if (pdgid == 23 || pdgid == 25) {
	particles["ZH"].push_back(i);
	GetSimplifiedModelParticles(particles, truth_evtxid[i]);      
    }
    else if (pdgid == 0) {
      // ignore
    }
    else { //other particle found -- ignore
    }
  }
}

int Analysis::IsSimplifiedModel(map<string, vector<int> > & particles)
{
  GetSimplifiedModelParticles(particles, 1);
  
  if (particles["mu"].size() >= 2 &&
      particles["q"].size() >= 2) {
    
    if (particles["x+-2"].size() > 0) { // Chargino 2
      return 7;
    }

    else if (particles["x+-1"].size() > 0) { // Chargino 1
      return 6;
    }

    else if (particles["x04"].size() > 0) { // Neutralino 4 decay 
      return 5;
    }

    else if (particles["x03"].size() > 0) { // Neutralino 3 decay
      return 4;
    }

    else if (particles["x02"].size() > 0) { // Neutralino 2 decay
      return 3;
    }
    else if (particles["x01"].size() > 0 &&
	     (particles["ZH"].size() > 0 || particles["W"].size() > 0)) { // Neutralino 1 decay
      return 2;
    }
    else if (particles["mu"].size() == 2 &&
	     //truth_pdgid[particles['m'][0]] == truth_pdgid[particles['m'][1]] &&
	     particles["x01"].size() >= 1 &&
	     particles["q"].size() == 2 && 
	     particles["ZH"].size() == 0 &&
	     particles["W"].size() == 0) { // SMS
      return 1;
    }
    else {
      // other case, ignore for now
    }

  }
  return -1;
}

Bool_t Analysis::TriggerMatched(Int_t muonIterator, vector <TString> triggerFilters) {
  // muon trigger matching
  for (unsigned int k = 0; k < (*trig_filter).size(); k++) {
    TString name = (TString) (*trig_filter).at(k);
    for (unsigned int l = 0; l < triggerFilters.size(); l++) {
      // check for matching triggerfilter names and check dR
      if ( name.EqualTo(triggerFilters.at(l)) )
	if ( dRMatched(muonIterator, k) )
	  return true;
    }
  }
  
  return false;  
}

Bool_t Analysis::dRMatched(Int_t muonIterator, Int_t triggerID) {
  // helper method for muon trigger matching
  for(int k = 0; k < trigFilter_n; k++) {
    Double_t mindeltaR = 0.5;
    if( trig_filterid[k] == triggerID ) {
      Double_t dEta = fabs( muo_eta[muonIterator] - trig_eta[k] );
      Double_t dPhi = fabs( DeltaPhi(muo_phi[muonIterator], trig_phi[k]) );
      Double_t dR   = sqrt( pow(dEta, 2) + pow(dPhi, 2) );
         if( dR <= mindeltaR ) {
	   return true;
	 }
    }
  }
  return false;
}

Bool_t Analysis::MuonCuts (Int_t i) {
  // see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId, tight muon selection
  TCutList MuonCuts(histo, global_weight);

  // set cuts
  MuonCuts.Set("nMuon_pt", muo_pt[i] > 7., muo_pt[i]);
  MuonCuts.Set("nMuon_eta", fabs(muo_eta[i]) <= 2.1, fabs(muo_eta[i]));

  MuonCuts.Set("nMuon_id", muo_ID[i][1] == 1, muo_ID[i][1]); // Global Muon
  MuonCuts.Set("nMuon_ispf", muo_isPFMuon[i] == 1, muo_isPFMuon[i]);
  MuonCuts.Set("nMuon_d0Tk", fabs(muo_d0Tk[i]) < 0.2, fabs(muo_d0Tk[i])); //dxy
  MuonCuts.Set("nMuon_ValidPixelHitsCm", muo_ValidPixelHitsCm[i] > 0,
	       muo_ValidPixelHitsCm[i]);
  MuonCuts.Set("nMuon_TrkChiNormCm", muo_TrkChiNormCm[i] < 10., muo_TrkChiNormCm[i]);
  MuonCuts.Set("nMuon_hitsCm", muo_ValidMuonHitsCm[i] > 0, muo_ValidMuonHitsCm[i]);

  MuonCuts.Set("nMuon_StationsMatched", muo_StationsMatched[i] > 1, muo_StationsMatched[i]);
  MuonCuts.Set("nMuon_dzTk", fabs(muo_dzTk[i]) < 0.5, fabs(muo_dzTk[i]));
  MuonCuts.Set("nMuon_TrackerLayersMeasCm", muo_TrackerLayersMeasCm[i] > 5,
	       muo_TrackerLayersMeasCm[i]);


  MuonCuts.FillHistograms();
  return MuonCuts.PassesAll();
}


Bool_t Analysis::JetCuts (Int_t i, vector <int> muons) {
  // https://twiki.cern.ch/twiki/bin/view/CMS/JetID
  TCutList JetCuts(histo, global_weight);

  TLorentzVector jet(pfjet_px[i], pfjet_py[i], pfjet_pz[i], pfjet_E[i]);
    
  // set cuts
  JetCuts.Set("nPfJet_pt", pfjet_pt[i] > 15., pfjet_pt[i]);
  JetCuts.Set("nPfJet_eta", fabs(pfjet_eta[i]) < 2.4, fabs(pfjet_eta[i]));
  JetCuts.Set("nPfJet_const", pfjet_const[i] > 1, pfjet_const[i]);
  JetCuts.Set("nPfJet_pff_0", pfjet_PFF[i][0] > 0., pfjet_PFF[i][0]); // Charged Hadron Energy Fraction
  JetCuts.Set("nPfJet_pff_1", pfjet_PFF[i][1] < 0.99, pfjet_PFF[i][1]); // Neutral Hadron Energy Fraction
  JetCuts.Set("nPfJet_pff_2", pfjet_PFF[i][2] < 0.99, pfjet_PFF[i][2]); // Neutral Em Energy Fraction
  JetCuts.Set("nPfJet_pff_3", pfjet_PFF[i][3] < 0.99, pfjet_PFF[i][3]); // Charged Em Energy Fraction
  JetCuts.Set("nPfJet_pfn_0", pfjet_PFN[i][0] > 0, pfjet_PFN[i][0]); // Charged Hadron Multiplicity

  // A muon might be reconstructed as a PF jet. Therefore take
  // only those PF jets that do not "coincide" with a tight muon
  double Rmin = 1000;
  for (unsigned int mu = 0; mu < muons.size(); mu++) {
    int m = muons[mu];
    TLorentzVector mu(muo_px[m], muo_py[m], muo_pz[m], muo_E[m]);
    Rmin = TMath::Min(Rmin, jet.DeltaR(mu));
  }
  JetCuts.Set("nPfJet_rmin", Rmin > 0.05, Rmin);
  JetCuts.FillHistograms();

  Fill("PfJet_m", jet.M());
  // return true if jet passes all cuts
  return JetCuts.PassesAll();
}


Bool_t Analysis::ElectronCuts (Int_t i) {
  // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaCutBasedIdentification    
  TCutList ElectronCuts(histo, global_weight);
  // cut values for WP80 / medium
  const Double_t cutdEtaSCTrackAtVtxB = 0.004;
  const Double_t cutdEtaSCTrackAtVtxE = 0.007; 
  const Double_t cutdPhiSCTrackAtVtxB = 0.06;
  const Double_t cutdPhiSCTrackAtVtxE = 0.03;
  const Double_t cutSigmaIetaIetaB = 0.01;
  const Double_t cutSigmaIetaIetaE = 0.03;
  const Double_t cutHCalOverEmB = 0.12;
  const Double_t cutHCalOverEmE = 0.10;
  const Double_t cutd0vtx = 0.02;
  const Double_t cutdZvtx = 0.1;
  const Double_t cut1OverEminus1OverP = 0.05;
  const Double_t cutPFisoOverPt = 0.15;
  const Double_t cutMissingHits = 1.; //  
      
  // pfiso
  Double_t rhoPrime = max(global_rho, 0.0);
  Double_t iso_n = 0, AEff = 0;
  if (fabs(ele_SCeta[i]) >= 0.0 && fabs(ele_SCeta[i]) < 1.0 ) AEff = 0.210;
  if (fabs(ele_SCeta[i]) >= 1.0 && fabs(ele_SCeta[i]) < 1.479 ) AEff = 0.210;
  if (fabs(ele_SCeta[i]) >= 1.479 && fabs(ele_SCeta[i]) < 2.0 ) AEff = 0.110;
  if (fabs(ele_SCeta[i]) >= 2.0 && fabs(ele_SCeta[i]) < 2.2 ) AEff = 0.140;
  if (fabs(ele_SCeta[i]) >= 2.2 && fabs(ele_SCeta[i]) < 2.3 ) AEff = 0.180;
  if (fabs(ele_SCeta[i]) >= 2.3 && fabs(ele_SCeta[i]) < 2.4 ) AEff = 0.190;
  if (fabs(ele_SCeta[i]) >= 2.4) AEff = 0.260;
  iso_n = ele_PFiso[i][0] + max( ele_PFiso[i][1] + ele_PFiso[i][2] - rhoPrime*AEff, 0.0 ); 

  const Double_t cutHCalIsoOpt =  0.2;
  const Double_t cutECalIsoOpt =  0.2;
  const Double_t cutTkIsoOpt =  0.2;

  // own cuts
  ElectronCuts.Set("nElectron_pt", ele_pt[i] > 15.0, ele_pt[i]);
      
  // section specific cuts
  if (fabs(ele_SCeta[i])<1.4442) {
    // Barrel
    ElectronCuts.Set("nElectron_DEta", fabs(ele_dEtaSCTrackAtVtx[i]) < cutdEtaSCTrackAtVtxB, fabs(ele_dEtaSCTrackAtVtx[i]));
    ElectronCuts.Set("nElectron_DPhi", fabs(ele_dPhiSCTrackAtVtx[i]) < cutdPhiSCTrackAtVtxB, fabs(ele_dPhiSCTrackAtVtx[i]));
    ElectronCuts.Set("nElectron_SigmaIEtaIEta", ele_SigmaIetaIeta[i] < cutSigmaIetaIetaB, ele_SigmaIetaIeta[i]);
    ElectronCuts.Set("nElectron_HoE", ele_HCalOverEm[i] < cutHCalOverEmB, ele_HCalOverEm[i]);
  }
  else if (1.566 < fabs(ele_SCeta[i]) && fabs(ele_SCeta[i]) < 2.5) {
    // Endcap
    ElectronCuts.Set("nElectron_DEta", fabs(ele_dEtaSCTrackAtVtx[i]) < cutdEtaSCTrackAtVtxE, fabs(ele_dEtaSCTrackAtVtx[i]));
    ElectronCuts.Set("nElectron_DPhi", fabs(ele_dPhiSCTrackAtVtx[i]) < cutdPhiSCTrackAtVtxE, fabs(ele_dPhiSCTrackAtVtx[i]));
    ElectronCuts.Set("nElectron_SigmaIEtaIEta", ele_SigmaIetaIeta[i] < cutSigmaIetaIetaE, ele_SigmaIetaIeta[i]);
    ElectronCuts.Set("nElectron_HoE", ele_HCalOverEm[i] < cutHCalOverEmE, ele_HCalOverEm[i]);
  }
  else {
    // ignore electrons in gap between barrel and endcap
  }
  // common cuts
  ElectronCuts.Set("nElectron_d0vtx", ele_d0vtx[i] <= cutd0vtx, ele_d0vtx[i]);
  ElectronCuts.Set("nElectron_1oE-1oP", fabs(1/ele_E[i] -  ele_EoverP[i]/ele_E[i]) < cut1OverEminus1OverP,
		   fabs(1/ele_E[i] -  ele_EoverP[i]/ele_E[i]));
  ElectronCuts.Set("nElectron_pfIsoOpt", iso_n/ele_Et[i] < cutPFisoOverPt, iso_n/ele_Et[i]);
  ElectronCuts.Set("nElectron_MissingHits", ele_TrkExpHitsInner[i] <= cutMissingHits, ele_TrkExpHitsInner[i]);

  ElectronCuts.Set("nElectron_dZvtx", ele_dzvtx[i] <= cutdZvtx, ele_dzvtx[i]);
  ElectronCuts.Set("nElectron_HasMatchedConversions", ele_HasMatchedConversions[i] == true, ele_HasMatchedConversions[i]);

  // isolation
  ElectronCuts.Set("nElectron_HCalIsoOpt", ele_Dr03HCalSumEt[i]/ele_pt[i] < cutHCalIsoOpt, ele_Dr03HCalSumEt[i]/ele_pt[i]);
  ElectronCuts.Set("nElectron_ECalIsoOpt", ele_Dr03ECalSumEt[i]/ele_pt[i] < cutECalIsoOpt, ele_Dr03ECalSumEt[i]/ele_pt[i]);
  ElectronCuts.Set("nElectron_TkIsoOpt" ,    ele_Dr03TkSumPt[i]/ele_pt[i] < cutTkIsoOpt, ele_Dr03TkSumPt[i]/ele_pt[i]);
      
  ElectronCuts.FillHistograms();
  return ElectronCuts.PassesAll();
}

truth_pt_comparator::truth_pt_comparator(const Analysis & analysis)
  : fAnalysis(analysis)
{
}


bool truth_pt_comparator::operator()(int i, int j)
{
  return fAnalysis.truth_pt[i] > fAnalysis.truth_pt[j];
}

void Analysis::SignalStudy(int & charge)
{
  // muons, neutrinos and quarks in final state
  int nLeptons  = 0;
  int nGauginos = 0;
  int nSleptons = 0;
  int nBosons   = 0;
  int nHiggs    = 0;

  // slepton four-momentum
  TLorentzVector slepton;
  TLorentzVector resonance;
  
  for (int i = 0; i < truth_n; i++) {
    // sum up four-momenta of final state particles
    if (truth_evtxid[i] == -1 && truth_pdgid[i] != 0 &&
	(truth_bvtxid[i] != 1 || abs(truth_pdgid[i]) > 6)) {
      TLorentzVector particle(truth_px[i], truth_py[i], truth_pz[i], truth_E[i]);
      slepton += particle;
    }

    // find particles produced at first vertex
    if (truth_bvtxid[i] != 1)
      continue;

    // sum up particles at vertex 1, cross-check with final state particles
    if (truth_pdgid[i] != 0 && abs(truth_pdgid[i]) > 6) {
      TLorentzVector particle(truth_px[i], truth_py[i], truth_pz[i], truth_E[i]);
      resonance += particle;
    }

    // analyse particles
    switch (abs(truth_pdgid[i])) {
      // muon
      case 13: // \mu^-
	nLeptons++;
	if (truth_pdgid[i] > 0)
	  charge--;
	else
	  charge++;
	break;
	// muon neutrino
      case 14: // \nu_\mu
	nLeptons++;
	break;
	// ignore quarks and status 0 particles
	// this is because SUSYAna does something strange when reading HERWIG info
      case 0:
      case 1:
      case 2:
	break;
	// chargino
      case 1000024: // \tilde{\chi}_1^+
      case 1000037: // \tilde{\chi}_2^+
	nGauginos++;
	if (truth_pdgid[i] > 0)
	  charge++;
	else
	  charge--;
	break;
	// neutralino
      case 1000022: // \tilde{\chi}_1^0
      case 1000023: // \tilde{\chi}_2^0
      case 1000025: // \tilde{\chi}_3^0
      case 1000035: // \tilde{\chi}_4^0
	nGauginos++;
	break;
      case 1000013: // \tilde{\mu}_L^-
      case 2000013: // \tilde{\mu}_R^-
	nSleptons++;
	if (truth_pdgid[i] > 0)
	  charge--;
	else
	  charge++;
	break;
      case 1000014: // \tilde{\nu}_{\mu\,L}
	nSleptons++;
	break;
      case 24: // W^+
	nBosons++;
	if (truth_pdgid[i] > 0)
	  charge++;
	else
	  charge--;
	break;
	// Higgs
      case 25: // h^0
	nHiggs++;
	nBosons++;
	break;
      default:
	TruthDump();
	THROW("Unexpected particle found in Signal sample");
    }
  }

  

  // do some cross-checks
  if ((nLeptons != 1 || nGauginos != 1) && (nSleptons != 1 || nBosons != 1)) {
    TruthDump();
    THROW("This does not look like signal - expect one muon/neutrino + one gaugino\n"
	  " or one sneutrino/smuon + one boson");
  }

  // follow the decay chain to find out if there are two muons (no neutrinos)
  map<char, vector<int> > particles;
  ReconstructFinalState(particles, 1);
  // is it signal (i.e. at least two muons and at least two jets) or not?
  if (particles['m'].size() >= 2 && particles['q'].size() >= 2) {
    fIsSignal = true;
  }
  else {
    fIsSignal = false;
  }

  // take only signal events or background events if requested
  if ((fInputType == "signal" && !fIsSignal) || (fInputType == "background" && fIsSignal))
    return;

  Fill("Sig_nElectron", particles['e'].size());
  Fill("Sig_nMuon", particles['m'].size());
  Fill("Sig_nTau", particles['t'].size());
  Fill("Sig_nLeptons", particles['l'].size());
  Fill("Sig_nNeutrino", particles['n'].size());
  Fill("Sig_nQuark", particles['q'].size());
  Fill("Sig_SleptonMass", slepton.M());
  Fill("Sig_SleptonCharge", charge);
  Fill("Sig_nBoson", nBosons);
  Fill("Sig_nHiggs", nHiggs);
  Fill("Sig_Resonance", resonance.M());

  // get true muon pt distribution, true quark pt distribution etc.
  double ptMax=-1E9;
  double ptMin=1E9;
  for (vector<int>::const_iterator it = particles['m'].begin();
       it != particles['m'].end(); it++) {
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
  
  Fill("Sig_signal_counter", 0.);
  map<string, vector<int> > simpleParticles;
  int modelid = IsSimplifiedModel(simpleParticles);
  for (int i = 1; i<=7; i++) {
    if (modelid == i) {
      if (modelid == 6 && fabs(charge) == 0) {
	Fill("Sig_signal_counter", 8);
	continue;
      }
      if (modelid == 7 && fabs(charge) == 0) {
	Fill("Sig_signal_counter", 9);
	continue;
      }
      Fill("Sig_signal_counter",i);
    } 
    else if (modelid == -1) {
      //TruthDump();
    }
  }

  // compute lepton kinematics
  truth_pt_comparator comp(*this);
  sort(particles['l'].begin(), particles['l'].end(), comp);
  int id1 = particles['l'][0];
  int id2 = particles['l'][1];
  fSigMu0->SetXYZT(truth_px[id1], truth_py[id1], truth_pz[id1], truth_E[id1]);
  fSigMu1->SetXYZT(truth_px[id2], truth_py[id2], truth_pz[id2], truth_E[id2]);
  Fill("Sig_MuMass", (*fSigMu0+*fSigMu1).M());
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
  fSigJet0->SetXYZT(truth_px[id3], truth_py[id3], truth_pz[id3], truth_E[id3]);
  fSigJet1->SetXYZT(truth_px[id4], truth_py[id4], truth_pz[id4], truth_E[id4]);
  Fill("Sig_JetMass", (*fSigJet0+*fSigJet1).M());
  Fill("Sig_JetDPhi", DeltaPhi(truth_phi[id3], truth_phi[id4]));
  Fill("Sig_JetDEta", TMath::Abs(truth_eta[id3]-truth_eta[id4]));
  Fill("Sig_JetDPt", TMath::Abs(truth_pt[id3]-truth_pt[id4]));
  Fill("Sig_JetSumPt", TMath::Abs(truth_pt[id3]+truth_pt[id4]));
  Fill("Sig_JetPt", truth_pt[id3], truth_pt[id4]);

  // look at jets and leptons
  TLorentzVector neutralino = *fSigMu1+*fSigJet0+*fSigJet1;
  Fill("Sig_Mass3", neutralino.M());
  Fill("Sig_Angle3", neutralino.DeltaPhi(*fSigMu0));
  double DeltaR = TMath::Min(
    TMath::Min(fSigMu0->DeltaR(*fSigJet0), fSigMu0->DeltaR(*fSigJet1)),
    TMath::Min(fSigMu1->DeltaR(*fSigJet0), fSigMu1->DeltaR(*fSigJet1))
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
      for (int j = 0; j < muo_n; j++) {
	// does this reco muon have a match to the gen lepton?
	if (muo_truth[j] == i)
	  continue;
      }
      Fill("Sig_LostMuEta", truthl_eta[i]);
      Fill("Sig_LostMuPt", truthl_pt[i]);
      Fill("Sig_LostMuPtEta", truthl_pt[i], truthl_eta[i]);
    }
  }
}

void Analysis::TightLooseRatioCalculation(const vector<int> & loose_muons,
					  const vector<int> & tight_muons,
					  const vector<int> & jets,
					  const double HT)
{
  // compute variables
  int nMuonsLoose = loose_muons.size();
  int nMuonsTight = tight_muons.size();

  // fill histograms before cuts
  Fill("TL_met", met_et[MET_INDEX]);
  Fill("TL_metsig", met_etsignif[MET_INDEX]);
  Fill("TL_jetpt", fJets > 0 ? pfjet_pt[jets[0]] : 0);
  Fill("TL_ht", HT);
  Fill("TL_nloose", nMuonsLoose);

  // global requirements to get more QCD like events
  TCutList QCDCuts(histo, global_weight);
  // mainly against ttbar, and a bit w->l nu
  QCDCuts.Set("nTL_met", met_et[MET_INDEX] < fTL_met_max, met_et[MET_INDEX]);
  // remove Drell-Yan and signal
  QCDCuts.Set("nTL_ht", fTL_ht_min <= HT, HT);

  // trigger
  QCDCuts.Set("nTL_jetpt", fJets > 0 ? fTL_jetpt_min <= pfjet_pt[jets[0]]  : false,
	      fJets ? pfjet_pt[jets[0]] : 0);
  // (for fakes) && (against Drell-Yan)
  QCDCuts.Set("nTL_nloose", (fTL_nloose_min <= nMuonsLoose) && (nMuonsLoose < fTL_nloose_max),
	      nMuonsLoose);
  // require some jets
  QCDCuts.Set("nTL_njets", fTL_njets_min <= fJets, fJets);


  double ST = HT;
  for (Int_t i = 0; i < nMuonsLoose; i++) {
    // trigger specific muon cut
    QCDCuts.Set("nTL_firstmupt", fTL_firstmupt_min <= muo_pt[0] , muo_pt[0]);
    // muon dependent requirements to get more QCD like events
    Fill("TL_metdphi", DeltaPhi(met_phi[MET_INDEX], muo_phi[loose_muons[i]]));
    double dphi = 0.;
    if (fJets > 0)
      dphi = DeltaPhi(pfjet_phi[jets[0]], muo_phi[loose_muons[i]]);
    Fill("TL_jetdphi", dphi);

    TLorentzVector mu0(muo_px[loose_muons[i]], muo_py[loose_muons[i]],
		       muo_pz[loose_muons[i]], muo_E[loose_muons[i]]);
    ST += muo_pt[loose_muons[i]];
    // find another (very very) loose muon and construct "Z mass" hypothesis
    double mindz = 1000;
    double zmass = 0.;
    for (Int_t j = 0; j < muo_n; j++) {
      // skip the same muon
      if (j == loose_muons[i])
	continue;
      TLorentzVector mu1(muo_px[j], muo_py[j], muo_pz[j], muo_E[j]);
      double m = (mu0+mu1).M();
      Fill("TL_zmass", m);
      // find mass closest to Z mass
      zmass = 91.+TMath::Sign(1., m-91.)*TMath::Min(fabs(91.-zmass), fabs(91.-m));


     
      Fill("TL_mumudz", fabs( muo_dzTk[j] - muo_dzTk[loose_muons[i]] ));
      mindz = TMath::Min( mindz, fabs(muo_dzTk[j] - muo_dzTk[loose_muons[i]]) );
    }
    QCDCuts.Set("nTL_zmass", zmass > fTL_zmass_min && fTL_zmass_max > zmass, zmass);
    QCDCuts.Set("nTL_mumudz", mindz > fTL_mumudz_min, mindz);

    // compute transverse mass of MET and mu (W hypothesis)
    double cos_phi12 = TMath::Cos(DeltaPhi(met_phi[MET_INDEX], mu0.Phi()));
    double MT = TMath::Sqrt(2*met_et[MET_INDEX]*mu0.Pt()*(1.-cos_phi12));
    Fill("TL_mt", MT);

    // muon pt requirement - needed for MC and data comparison
    QCDCuts.Set("nTL_mupt", fTL_mupt_min <= muo_pt[loose_muons[i]] , muo_pt[loose_muons[i]]);
    // require muons to be on "opposite" side to jet in phi
    QCDCuts.Set("nTL_jetdphi", fTL_jetdphi_min <= dphi, dphi);
    // reject W->l nu
    QCDCuts.Set("nTL_mt", MT <= fTL_mt_max, MT);
  }

  // Fill N-1 cuts
  QCDCuts.FillHistograms();
  // apply cuts
  if (!QCDCuts.PassesAll())
    return;
  
  // if( (fInputType == "data" || fInputType == "mc") && muo_pt[loose_muons[0]] > 59. ) {
  //     INFO("runid: " << global_run << ":" << lumi_section << ":" << global_event);
  //   }
  
  // double-check that there is only one loose muon
  assert(nMuonsLoose == 1);
  assert(fJets >= 1);

  // some more histograms
  Fill("TL_nmutight", nMuonsTight);
  Fill("TL_st", ST);

  // fill histograms in order to compute tight/loose ratio
  int i = 0;
  Fill("LooseMuons",
       muo_pt[loose_muons[i]],
       fabs(muo_eta[loose_muons[i]]),
       pfjet_pt[jets[0]]);
  // is it also a tight muon?
  for (Int_t j = 0; j < nMuonsTight; j++) {
    if (loose_muons[i] == tight_muons[j]) {
      // is a tight muon!
      Fill("TightMuons",
	   muo_pt[tight_muons[j]],
	   fabs(muo_eta[tight_muons[j]]),
	   pfjet_pt[jets[0]]);
    }
  }
}

// main event loop
void Analysis::Loop()
{
  // process and memory info
  ProcInfo_t info;

  // create variables
  CreateVariables();

  // set branch addresses in output tree if necessary
  if (fFill) {
    SetBranchAddresses();
    CreateBranches();
  }
  
  // create histograms
  CreateHistograms();

  // create lumi filter (based on JSON file from configuration)
  lumi::RunLumiRanges runcfg(fLumiRanges);
  
  // create event filter (based on zipped filter file from configuration)
  EventFilter hcalevtcfg(fEventFilter.c_str());

  // main event loop
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  INFO("Analysis: Chain contains " << nentries << " events");
  if (fMaxEvents > 0)
    INFO("Stopping after fMaxEvents = " << fMaxEvents);

  TString trigger(80);

  bool firstWarnWeight = true;
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
    
    FillNoWeight("runnumber", global_run);
    FillNoWeight("pu_TrueNrInter", pu_TrueNrInter);

    //////////////////////////////////////////////////////////////////////
    // Dump event/object information, see Dump.cxx
    if (fDumpBasic || fDumpAll) {
      BasicDump(jentry);
    }
    if (fDumpTruth || fDumpAll) {
      TruthDump();
    }
    if (fDumpTrigger != "None"  || fDumpAll) {
      if (!fDumpAll)
	TriggerDump(fDumpTrigger); // select e.g all muon trigger with "Mu"
      else
	TriggerDump("*");
    }
    if (fDumpAll) {
      VertexDump();
      MuonDump(0);      // [less] detailed information 1 [0]
      PFJetDump();
      METDump();
      EleDump(0);       // [less] detailed information 1 [0]
    }

    // find duplicates
    if (fFindDuplicates) {
      double tempval = global_qscale;
      if (FindDuplicates(global_run, global_event, lumi_section, tempval)) {
	// depending on severity one could stop here
	WARNING("Analysis: found duplicate ");
      }
    }

    //////////////////////////////////////////////////////////////////////
    // Check weight in MC events
    FillNoWeight("log_global_weight", TMath::Log(global_weight)/TMath::Log(10.));
    if (global_weight != 1 && firstWarnWeight) {
      WARNING("Initial weights are not 1.");
      firstWarnWeight = false;
    }
    global_weight = 1;

    Fill("cutflow", "weight"); 
    DEBUG("cutflow " << "weight"); 

    //////////////////////////////////////////////////////////////////////
    // Split RPV SUSY in RPV signal and RPV background
    fSigMu0->SetXYZT(0, 0, 0, 0);
    fSigMu1->SetXYZT(0, 0, 0, 0);
    fSigJet0->SetXYZT(0, 0, 0, 0);
    fSigJet1->SetXYZT(0, 0, 0, 0);
    
    // slepton charge
    int charge = 0;
  
    if (fInputType == "signal" || fInputType == "background") {
      try {
	SignalStudy(charge);
      }
      catch (std::string message) {
	WARNING(message);
      }
      // take only signal events or background events if requested
      if ((fInputType == "signal" && !fIsSignal) || (fInputType == "background" && fIsSignal))
	continue;
    }
    Fill("cutflow", "signal");

    //////////////////////////////////////////////////////////////////////
    // filter data on luminosity section - typically golden JSON file
    if(!runcfg.check(global_run, lumi_section)) {
      DEBUG("Filtering lumi section" << lumi_section << " in global run " << global_run);
      continue;
    }

    //////////////////////////////////////////////////////////////////////
    // filter data on HCAL laser filter
    if(!hcalevtcfg.filter(global_run, lumi_section, global_event)) {
      DEBUG("Filtering event " << global_event << " in lumi section " <<
	    lumi_section << " in global run " << global_run);
      continue;
    }


    Fill("cutflow", "evt&lumifilter");
    DEBUG("cutflow " << "evt&lumifilter");

    Fill("bReweight_vtx_n", vtx_n);

    //////////////////////////////////////////////////////////////////////
    // pileup reweighting
    if (fPileupReweighting && (fInputType == "mc" || fInputType == "signal"
			       || fInputType == "background")) {
      double weight = LumiWeights_.weight(pu_TrueNrInter);
      // if (weight != 0)
      global_weight *= weight;
      // else {
      // 	WARNING("pileup weight is zero - not using pileup reweighting for this event");
      // }
    }

    Fill("cutflow", "pileup rew.");
    DEBUG("cutflow " << "pileup rew.");

    //////////////////////////////////////////////////////////////////////
    // Redo skimmer cuts

    // fill histos before cuts
    Fill("bSkim_muo_n", muo_n);
    if (muo_n > 0)
      Fill("bSkim_muo_pt0", muo_pt[0]);
    if (muo_n > 1)
      Fill("bSkim_muo_pt1", muo_pt[1]);

    // fill N-1 histograms
    TCutList skimCuts(histo, global_weight);
    skimCuts.Set("nSkim_muo_n", muo_n >= fSkimMuons, muo_n);
    if (fSkimMuons > 0)
      skimCuts.Set("nSkim_muo_pt0", muo_n > 0 ? muo_pt[0] >= fSkimMuoptfirst : false,
		   muo_n > 0 ? muo_pt[0] : -1.);
    else
      skimCuts.Set("nSkim_muo_pt0", true, muo_n > 0 ? muo_pt[0] : -1.);
    if (fSkimMuons > 1)
      skimCuts.Set("nSkim_muo_pt1", muo_n > 1 ? muo_pt[1] >= fSkimMuoptother : false,
		   muo_n > 1 ? muo_pt[1] : -1.);
    else
      skimCuts.Set("nSkim_muo_pt1", true, muo_n > 1 ? muo_pt[1] : -1.);
    skimCuts.FillHistograms();

    // apply cuts
    if (fSkimActive) {
      if (!skimCuts.PassesAll()) {
	continue;
      }
    }
    Fill("cutflow", "redo skim");
    DEBUG("cutflow " << "redo skim");

    //////////////////////////////////////////////////////////////////////
    // Trigger selection

    bool rejection = true;
    // loop over all triggers
    for (int i = 0; i < trig_n && rejection; i++) {
      trigger = (*trig_name)[i].c_str();
      // loop over triggers which are accepted
      for (vector<string>::const_iterator it = fTrigger.begin();
	   it != fTrigger.end(); it++) {
	if (trigger.Contains(it->c_str())) {
	  // OK, trigger found
	  if (!fForceUnprescaledTriggers || trig_HLTprescale[i] == 1)
	    rejection = false;
	  break;
	}
      }
    }
    if (rejection)
      continue;

    Fill("cutflow", "trigger");
    DEBUG("cutflow " << "trigger");

    //////////////////////////////////////////////////////////////////////
    // PF jet smearing
    double old_met = met_et[MET_INDEX];
    PFJetSmearing();
    Fill("cutflow", "jetsmear");
    DEBUG("cutflow " << "jetsmear");

    //////////////////////////////////////////////////////////////////////
    // Object ID

    // Vertices
    int nVertex = 0;
    TCutList VertexCuts(histo, global_weight);
    for (int i = 0; i < vtx_n; i++) {
      // do computation
      double x = vtx_x[i] - bs_x;
      double y = vtx_y[i] - bs_y;
      double r = TMath::Sqrt(x*x+y*y);

      // set cut values
      VertexCuts.Set("nObj_vtx_fake", vtx_fake[i] == 0, vtx_fake[i]);
      VertexCuts.Set("nObj_vtx_ndof", vtx_ndof[i] > 4., vtx_ndof[i]);
      VertexCuts.Set("nObj_vtx_z", fabs(vtx_z[i]) < 24., vtx_z[i]);
      VertexCuts.Set("nObj_vtx_bs", r < 2.0, r);
      VertexCuts.FillHistograms();

      // apply cuts
      if (VertexCuts.PassesAll()) {
	// check if a cut on radial distance to beamspot would help
	nVertex++;
      }
    }




    /////////////////////////////////////////////////////////////////////
    // Muons
    vector <int> muons;
    for (int i = 0; i < muo_n; i++){
      if ( MuonCuts(i) ) {
	muons.push_back(i);
	TLorentzVector mu(muo_px[i], muo_py[i], muo_pz[i], muo_E[i]);
	Fill("Muon_m", mu.M());
      }
    }
    int nMuon = muons.size();

    /////////////////////////////////////////////////////////////////////  
    // Jets
    double HT = 0;
    vector <int> jets;
    
    for (int i = 0; i < pfjet_n; i++) {
      if ( JetCuts(i,muons) ) {
	if ( pfjet_pt[i] > 30. )
	  jets.push_back(i);
	HT += pfjet_pt[i];
      }
    }
    fJets = jets.size();

    //////////////////////////////////////////////////////////////////////
    // Electrons
    int nElectron = 0;

    for (int i = 0; i < ele_n; i++) {
      if ( ElectronCuts(i) )
	nElectron++;   
    }
    Fill("nElectron", nElectron);
    
    //////////////////////////////////////////////////////////////////////
    // loose and tight muon id - add isolation and trigger matching
    vector <int> loose_muons;
    vector <int> tight_muons;
    vector <double> tight_dR;
    vector <double> loose_dR;
    for (int i = 0; i < nMuon; i++) {
      vector <TString> muonTriggerFilters;
      // unique triggerfilters to match to
      muonTriggerFilters.push_back("hltDiMuonGlb17Glb8DzFiltered0p2"); // HLT_Mu17_Mu8_v
      muonTriggerFilters.push_back("hltDiMuonGlb17Trk8DzFiltered0p2"); // HLT_Mu17_TkMu8_v      
      //muonTriggerFilters.push_back("hltDiMuonGlb22Trk8DzFiltered0p2"); // HLT_Mu22_TkMu8_v
      //muonTriggerFilters.push_back("hltDiMuonGlb22Trk22DzFiltered0p2"); // HLT_Mu22_TkMu22_v
      
      //Fill("nMuon_matched", matched);
      if ( !TriggerMatched(i, muonTriggerFilters) )
	continue;

      // muo_PFiso
      // relative isolation

      //I = [ sumChargedHadronPt + max(0.,sumNeutralHadronPt + sumPhotonPt - 0.5 * sumPUPt) ] / pt
      int m = muons[i];
      double relIso = ( muo_PFiso[m][0] + max(0., muo_PFiso[m][2] + muo_PFiso[m][3] - 0.5 * muo_PFiso[m][6]) ) /  muo_pt[m];


      // distance to next jet
      double Rmin = 1000;
      Int_t jet_min = -1;
      TLorentzVector mu(muo_px[m], muo_py[m], muo_pz[m], muo_E[m]);
      for (int j = 0; j < fJets; j++) {
	TLorentzVector jet(pfjet_px[jets[j]], pfjet_py[jets[j]],
			   pfjet_pz[jets[j]], pfjet_E[jets[j]]);
	double R = mu.DeltaR(jet);
	if (R < Rmin) {
	  Rmin = R;
	  jet_min = j;
	}
      }
      Fill("MJ_rmin", Rmin);
      if (jet_min > 0 && Rmin < 0.5) {
	Fill("relIsoInJet", relIso);
	Fill("MJ_EF", muo_E[m]/pfjet_E[jets[jet_min]]);
	Fill("MJ_EE", muo_E[m], pfjet_E[jets[jet_min]]);
	// use MC to help if there is a signal particle close
	if (fIsSignal) {
	  double min_mu_dR = TMath::Min(
	    mu.DeltaR(*fSigMu0), mu.DeltaR(*fSigMu1)
	    );
	  Fill("MJ_mu_dR", min_mu_dR);
	  double min_jet_dR = TMath::Min(
	    mu.DeltaR(*fSigJet0), mu.DeltaR(*fSigMu1)
	    );
	  Fill("MJ_jet_dR", min_jet_dR);
	  Fill("MJ_dR", min_mu_dR, min_jet_dR);
	}
	//continue;
      } else {
	Fill("relIsoOutJet", relIso);
      }

      Fill("nMuon_relIso", relIso);
      if (relIso < fLooseMuonRelIso) {
	loose_muons.push_back(muons[i]);
	loose_dR.push_back(Rmin);
      }
      if (relIso < 0.12 && fabs(muo_d0Tk[i]) < 0.02) {
	tight_muons.push_back(muons[i]);
	tight_dR.push_back(Rmin);
      }
    }
    Int_t nMuonsLoose = loose_muons.size();
    Int_t nMuonsTight = tight_muons.size();


    // TEST VTX_N
    Double_t vtx_mindz = 1000.;
    //Double_t vtx_2ndmindz = 1000.;
    for (int i = 0; i < vtx_n; i++) {
      for (int j = 0; j < vtx_n; j++) {
	if( i != j) {
	  vtx_mindz = min( vtx_mindz, fabs(vtx_z[i]-vtx_z[j]));
	}
      }
    }
    Fill("delta_vtx_z", vtx_mindz);
    

    // Photons
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/Vgamma2011PhotonID
    // Variable 	Method 	recommended cuts
    // Pixel Seed Veto 	hasPixelSeed() 	YES
    // H/E 	hadronicOverEm() 	< 0.05
    // shower shape 	sigmaIetaIeta() 	< 0.011 (EB), < 0.03 (EE)
    // hollow cone track isolation 	trkSumPtHollowConeDR04() 	< 2.0 + 0.001 x Et + 0.0167 x rho25 (EB), < 2.0 + 0.001 x Et + 0.032 x rho25 (EE)
    // Jurrasic ECAL Isolation 	ecalRecHitSumEtConeDR04() 	< 4.2 + 0.006 x Et + 0.183 x rho25 (EB), < 4.2 + 0.006 x Et + 0.090 x rho25 (EE)
    // tower-based HCAL Isolation 	hcalTowerSumEtConeDR04() 	< 2.2 + 0.0025 x Et + 0.062 x rho25 (EB), < 2.2 + 0.0025 x Et + 0.180 x rho25 (EE)

    // Additionally apply spike cleaning: sigmaIEtaIEta > 0.001 and sigmaIPhiIPhi > 0.001 in Barrel region only
    // Et is the photon transverse energy
    // rho25 is obtained by the following piece of code

    // process.load('RecoJets.Configuration.RecoPFJets_cff')
    // process.kt6PFJets25 = process.kt6PFJets.clone( doRhoFastjet = True )
    // process.kt6PFJets25.Rho_EtaMax = cms.double(2.5)
    // process.fjSequence25 = cms.Sequence( process.kt6PFJets25 )
    
    DEBUG("cutflow " << "objectID");

    // Fill vertex plots before event cleaning and object cuts
    Fill("bSkim_vtx_n", vtx_n);
    for (int i = 0; i < vtx_n; i++) {
      Fill("bSkim_vtx_ntr", vtx_ntr[i]);
      Fill("bSkim_vtx_ndof", vtx_ndof[i]);
      Fill("bSkim_vtx_x", vtx_x[i]);
      Fill("bSkim_vtx_y", vtx_y[i]);
      Fill("bSkim_vtx_z", vtx_z[i]);
    }

    //////////////////////////////////////////////////////////////////////
    // Event cleaning    
    TCutList EventFilterList(histo, global_weight);

    for (Int_t i = 0; i < eventfilter_n; i++) {
      TString filterName = eventfilter_names->at(i);
      // ignoring the following filters
      if ( filterName.Contains("scrapingFilter")       ||
	   filterName.Contains("ECALDeadCellFilterBE") ||
	   filterName.Contains("MuonFailureFilter")    ||
	   ( fInputType != "data"  &&
	     filterName.Contains("HCALLaserFilterFromTriggerResult")) ) {
	continue;
      }

      EventFilterList.Set(eventfilter_names->at(i), eventfilter_results[i], eventfilter_results[i]);
    }
    EventFilterList.FillHistograms();
    
    if ( nVertex < 1 || nElectron > 0 || !EventFilterList.PassesAll() ) continue;

    Fill("cutflow", "cleaning");
    DEBUG("cutflow " << "cleaning");

    //////////////////////////////////////////////////////////////////////
    // fake rate study: tight-to-loose ratio calculation
    TightLooseRatioCalculation(loose_muons, tight_muons, jets, HT);

    //////////////////////////////////////////////////////////////////////
    // loose muon id
    if (nMuon < 2 || muo_pt[muons[0]] <= 20. || muo_pt[muons[1]] <= 15.) {
      continue;
    }
    Fill("cutflow", "loose muon");
    DEBUG("cutflow " << "loose muon");

    // for cross-checks: loop over muons, find two muons compatible with z mass
    // and plot number of vertices, mass, number of muons, number of jets, ...
    Fill("check_nmuon", nMuon);
    Fill("check_njets", fJets);
    Fill("check_vtxn", vtx_n);
    for (Int_t i = 0; i < nMuon; i++) {
      for (Int_t j = i+1; j < nMuon; j++) {
	TLorentzVector mu0(muo_px[muons[i]], muo_py[muons[i]], muo_pz[muons[i]],
			   muo_E[muons[i]]);
	TLorentzVector mu1(muo_px[muons[j]], muo_py[muons[j]], muo_pz[muons[j]],
			   muo_E[muons[j]]);
	Fill("check_mumass", (mu0+mu1).M());
      }
    }

    //////////////////////////////////////////////////////////////////////
    // loose jet id
    Fill("Jet_pt0", fJets > 0 ? pfjet_pt[jets[0]] : 0);
    Fill("Jet_pt1", fJets > 1 ? pfjet_pt[jets[1]] : 0);
    Fill("Jet_pt", fJets > 0 ? pfjet_pt[jets[0]] : 0, fJets > 1 ? pfjet_pt[jets[1]] : 0);
    if (fJets < 2 || pfjet_pt[jets[0]] <= 30. || pfjet_pt[jets[1]] <= 30.) {
      continue;
    }
    Fill("cutflow", "loose jet");
    DEBUG("cutflow " << "loose jet");

    //////////////////////////////////////////////////////////////////////
    // two muon requirements
    vector<double> muon_dR;
    Fill("Muon_nloose", nMuonsLoose);
    Fill("Muon_ntight", nMuonsTight);
    Fill("Muon_loosetight", nMuonsLoose, nMuonsTight);
    if (fAnalysisType == "doublefake") {
      // require two loose muons which are not tight
      if (nMuonsLoose != 2 || nMuonsTight != 0)
	continue;

      fMuoId[0] = loose_muons[0];
      muon_dR.push_back(loose_dR[0]);
      fMuoId[1] = loose_muons[1];
      muon_dR.push_back(loose_dR[1]);

      for (int n = 0; n < 2; n++) {
	int m = fMuoId[n];
	fFakeRate[n] = GetFakeRate(muo_pt[m], fabs(muo_eta[m]), pfjet_pt[jets[0]]);
	Fill("DF_pt", muo_pt[m]);
	Fill("DF_eta", fabs(muo_eta[m]));
      }

      // weights
      double w0 = fFakeRate[0]/(1.-fFakeRate[0]);
      double w1 = fFakeRate[1]/(1.-fFakeRate[1]);
      // to be used for event subtraction in single orthogonal sample
      fSingleFakeWeight = w0 + w1;
      // double orthogonal sample
      global_weight *= w0 * w1;
      Fill("cutflow", "doublefake");
      DEBUG("cutflow " << "doublefake");
    }
    else if (fAnalysisType == "singlefake") {
      // require one loose muon which is not tight and one tight muon
      if (nMuonsTight != 1 || nMuonsLoose != 2)
	continue;

      fMuoId[0] = tight_muons[0];
      muon_dR.push_back(tight_dR[0]);
      // find loose muon which is not tight muon
      for (int n = 0; n < nMuonsLoose; n++) {
	if (loose_muons[n] == tight_muons[0])
	  continue;
	const Int_t m = loose_muons[n];
	fMuoId[1] = m;
	muon_dR.push_back(loose_dR[n]);

	// weighting with T/L ratio
	fFakeRate[0] = GetFakeRate(muo_pt[m], fabs(muo_eta[m]), pfjet_pt[jets[0]]);
	fFakeRate[1] = -1.;
	Fill("SF_pt", muo_pt[m]);
	Fill("SF_eta", fabs(muo_eta[m]));
      }
      // assert only one loose muon
      assert(muon_dR.size() == 2);

      // swap IDs if necessary (in order to preserve pt ordering)
      if (fMuoId[0] > fMuoId[1]) {
	int temp = fMuoId[0];
	fMuoId[0] = fMuoId[1];
	fMuoId[1] = temp;
	double temp2 = muon_dR[0];
	muon_dR[0] = muon_dR[1];
	muon_dR[1] = temp2;
      }
      // event weight
      global_weight *= fFakeRate[0]/(1-fFakeRate[0]);
      Fill("cutflow", "singlefake");
      DEBUG("cutflow " << "singlefake");
    }
    else if (fAnalysisType == "default") {
      if (nMuonsTight != 2) {
	continue;
      }
      fMuoId[0] = tight_muons[0];
      muon_dR.push_back(tight_dR[0]);
      fMuoId[1] = tight_muons[1];
      muon_dR.push_back(tight_dR[1]);
      Fill("cutflow", "default");
      DEBUG("cutflow " << "default");
    }
    else {
      THROW("Unknown analysis type encountered: " + fAnalysisType);
    }

    //////////////////////////////////////////////////////////////////////
    // selection

    Fill("Muon_pt0", muo_pt[fMuoId[0]]);
    Fill("Muon_pt1", muo_pt[fMuoId[1]]);
    Fill("Muon_pt", muo_pt[fMuoId[0]], muo_pt[fMuoId[1]]);
    Fill("mumudz", fabs(muo_dzTk[fMuoId[0]]-muo_dzTk[fMuoId[1]]));
    if (muo_pt[fMuoId[0]] <= 20. ||
	muo_pt[fMuoId[1]] <= 15. ||
	muon_dR[0] < 0.4 ||
	muon_dR[1] < 0.4 ||
	fabs(muo_dzTk[fMuoId[0]]-muo_dzTk[fMuoId[1]]) > 0.08
	) {
      continue;
    }
    Fill("cutflow", "muonID");
    DEBUG("cutflow " << "muonID");

    //////////////////////////////////////////////////////////////////////
    // jet smeaing (JER)
    PFJetSmearingCalculation();

    // create Lorentz vectors from selected particles
    fMuon0->SetXYZT(muo_px[fMuoId[0]], muo_py[fMuoId[0]],
		    muo_pz[fMuoId[0]], muo_E[fMuoId[0]]);
    fMuon1->SetXYZT(muo_px[fMuoId[1]], muo_py[fMuoId[1]],
		    muo_pz[fMuoId[1]], muo_E[fMuoId[1]]);
    fJet0->SetXYZT(pfjet_px[jets[0]], pfjet_py[jets[0]], pfjet_pz[jets[0]],
		   pfjet_E[jets[0]]);
    fJet1->SetXYZT(pfjet_px[jets[1]], pfjet_py[jets[1]], pfjet_pz[jets[1]],
		   pfjet_E[jets[1]]);
    for (int n = 0; n < fJets && n < 40; n++) {
      fJetId[n] = jets[n];
    }
    *fGaugino = *fMuon1+*fJet0+*fJet1;
    *fSmuon = *fMuon0+*fGaugino;
    double m_mumu = (*fMuon0+*fMuon1).M();
    double GauginoMass = fGaugino->M();
    double SmuonMass = fSmuon->M();

    // needed because no Drell-Yan MC for m(mu,mu) < 10. GeV
    // M(mu, mu)
    if (m_mumu < 15.) {
      continue;
    }
    Fill("cutflow", "m_mumu");
    DEBUG("cutflow " << "m_mumu");

    //////////////////////////////////////////////////////////////////////
    // B-tagging
    fIsBTagged = false;
    for (int n = 0; n < fJets; n++) {
      int jj = jets[n];
      Int_t truth = pfjet_truth[jj];
      bool tagged = pfjet_btag[jj][6] > 0.679; // CVSM: http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2012_187_v5.pdf
      // fBTagging.isBJet(pfjet_btag[jj][6], 
      // 				   truth < 0 ? 0 : truth_pdgid[truth],
      // 				   pfjet_pt[jj], 
      // 				   pfjet_eta[jj]);
      Fill("isbtag", pfjet_btag[jj][6], tagged ? 1 : 0);
      if (tagged)
	fIsBTagged = true;
    }

    //////////////////////////////////////////////////////////////////////
    // MET cut with control regions
    Fill("pfmet_old", old_met);
    Fill("pfmet", met_et[MET_INDEX]);

    if (met_et[MET_INDEX] >= 50.) {
      // fill some histograms for inverse cut (needed for cross-checks), control region
      Fill("CR1_m_mumu", m_mumu);
      if (fIsBTagged) {
	// CR2: large MET and btag enhanced control region
	Fill("CR2_m_mumu", m_mumu);
      }
      else {
	// CR3: large MET and btag veto control region
	Fill("CR3_m_mumu", m_mumu);
      }
      if ((muo_charge[fMuoId[0]]*muo_charge[fMuoId[1]]) == 1) {
	// CR4: large MET, same charge cut
	Fill("CR4_m_mumu", m_mumu);
	Fill("CR4_m_smuon", SmuonMass);
	Fill("CR4_m_gaugino", GauginoMass);
	if (fIsBTagged) {
	  // CR5: large MET, same charge and btag enhanced control region
	  Fill("CR5_m_mumu", m_mumu);
	  Fill("CR5_m_smuon", SmuonMass);
	  Fill("CR5_m_gaugino", GauginoMass);
	}
	else {
	  // CR6: large MET, same charge and btag veto control region (signal-like)
	  Fill("CR6_m_mumu", m_mumu);
	  Fill("CR6_m_smuon", SmuonMass);
	  Fill("CR6_m_gaugino", GauginoMass);
	}
      }
      continue;
    }
    Fill("cutflow", "met");
    DEBUG("cutflow " << "met");

    // Delta_phi between second mu and Gaugino
    double delta_phi = fabs(fMuon0->DeltaPhi(*fGaugino));
    Fill("DeltaPhi", delta_phi);
    if (delta_phi < fDeltaPhiMin)
      continue;
    Fill("cutflow", "#Delta#phi");
    DEBUG("cutflow " << "#Delta#phi");

    // fill some histograms
    Fill("m_mumu_precut", m_mumu);
    Fill("muo_n_precut", muo_n);
    Fill("vtx_n", vtx_n);
    Fill("ht", HT);
    for (int i = 0; i < vtx_n; i++) {
      Fill("vtx_ntr", vtx_ntr[i]);
      Fill("vtx_ndof", vtx_ndof[i]);
      Fill("vtx_x", vtx_x[i]);
      Fill("vtx_y", vtx_y[i]);
      Fill("vtx_z", vtx_z[i]);
    }

    Fill("m_gaugino_precut", GauginoMass);
    Fill("m_smuon_precut", SmuonMass);

    // b-tag veto
    if (fIsBTagged) {
      continue;
    }
    Fill("cutflow", "btag veto");
    DEBUG("cutflow " << "btag veto");

    Fill("btag_m_mumu", m_mumu);
    Fill("btag_m_gaugino", GauginoMass);
    Fill("btag_m_smuon", SmuonMass);
    Fill("btag_jjmm_m", SmuonMass, GauginoMass);
    Fill("btag_m_smu_chi", SmuonMass, GauginoMass);     

    Fill("last_mumudz", fabs(muo_dzTk[fMuoId[0]]-muo_dzTk[fMuoId[1]]));

    // muon charge
    Fill("Muon_charge", muo_charge[fMuoId[0]], muo_charge[fMuoId[1]]);
    Fill("Muon_ch", muo_charge[fMuoId[0]]*muo_charge[fMuoId[1]]);
    if ((muo_charge[fMuoId[0]]*muo_charge[fMuoId[1]]) == -1) {
      continue;
    }
    Fill("cutflow", "charge");
    DEBUG("cutflow " << "charge");

    Fill("m_mumu", m_mumu);
    Fill("m_gaugino", GauginoMass);
    Fill("m_smuon", SmuonMass);
    Fill("jjmm_m", SmuonMass, GauginoMass); // final binning
    Fill("muo_n", muo_n);
    Fill("m_smu_chi", SmuonMass, GauginoMass);

    FillNoWeight("log_weight", TMath::Log(global_weight)/TMath::Log(10.));


    // Signal counter after cuts
    Fill("Cut_signal_counter", 0.);
    map<string, vector<int> > simpleParticles;
    int modelid = IsSimplifiedModel(simpleParticles);
    for (int i = 1; i<=7; i++) {
      if (modelid == i) {
	if (modelid == 6 && fabs(charge) == 0) {
	  Fill("Cut_signal_counter", 8);
	  continue;
	}
	if (modelid == 7 && fabs(charge) == 0) {
	  Fill("Cut_signal_counter", 9);
	  continue;
	}
	Fill("Cut_signal_counter",i);
      } 
      else if (modelid == 0) {
	//TruthDump();
      }
    }


    //////////////////////////////////////////////////////////////////////
    // fill output tree
    if (fFill) {
      fOutputTree.Fill();
    }
  }
}

void Analysis::CreateVariables()
{
  fSigMu0   = new TLorentzVector;
  fSigMu1   = new TLorentzVector;
  fSigJet0  = new TLorentzVector;
  fSigJet1  = new TLorentzVector;
  fMuon0    = new TLorentzVector;
  fMuon1    = new TLorentzVector;
  fJet0     = new TLorentzVector;
  fJet1     = new TLorentzVector;
  fSmuon    = new TLorentzVector;
  fGaugino  = new TLorentzVector;
}

void Analysis::CreateHistograms()
{
  //////////////////////////////////////////////////////////////////////
  // Histograms
  INFO("Analysis::CreateHistograms()");
  // use weighted histograms
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();

  // change to the file we want to have the histograms in
  TFile * outFile = fOutputTree.GetCurrentFile();
  outFile->cd();

  // cut flow - nasty workaround due to ROOTs problems with histograms filled with strings
  histo["cutflow"] = get_cutflow_histogram();
  const Int_t firstrun = 190000;
  const Int_t lastrun  = 220000;
  CreateHisto("runnumber", "Run number", lastrun-firstrun+1, firstrun, lastrun+1);
  CreateHisto("pu_TrueNrInter", "# true pileup interactions", 60, 0, 60);
  // signal histograms
  CreateHisto("Sig_nElectron", "Signal MC truth: number of electrons in final state", 10, -0.5, 9.5);
  CreateHisto("Sig_nMuon", "Signal MC truth: number of muons in final state", 10, -0.5, 9.5);
  CreateHisto("Sig_nTau", "Signal MC truth: number of #tau in final state", 10, -0.5, 9.5);

  CreateHisto("Sig_nLeptons", "Signal MC truth: number of leptons in final state", 10, -0.5, 9.5);
  CreateHisto("Sig_nNeutrino", "Signal MC truth, number of #nu in final state", 10, -0.5, 9.5);

  CreateHisto("Sig_nQuark", "Signal MC truth: number of quarks in final states", 10, -0.5, 9.5);
  CreateHisto("Sig_nHiggs", "Signal MC truth: Number of intermediate Higgs bosons", 10, -0.5, 9.5);
  CreateHisto("Sig_nBoson", "Signal MC truth: Number of intermediate bosons (including Higgs)", 10, -0.5, 9.5);
  CreateHisto("Sig_SleptonMass", "Resonance mass from final state particles@GeV", 1000, 0, 2000);
  CreateHisto("Sig_SleptonCharge", "Charge of slepton in units of e", 3, -1.5, 1.5);
  CreateHisto("Sig_Resonance", "Resonance mass from first vertex@GeV", 1000, 0, 2000);

  CreateHisto("Sig_ptMu1", "p_{T} of leading muon@GeV", 1000, 0, 1000);
  CreateHisto("Sig_ptMu2", "p_{T} of second leading muon@GeV", 1000, 0, 1000);
  CreateHisto("Sig_EtaMu", "#eta_{#mu}", 100, -5, 5);
  CreateHisto("Sig_PhiMu", "#phi_{#mu}@rad", 628, -3.14, 3.14);
  CreateHisto("Sig_ptQuark1", "p_{T} of leading quark@GeV", 1000, 0, 1000);
  CreateHisto("Sig_ptQuark2", "p_{T} of second leading quark@GeV", 1000, 0, 1000);
  CreateHisto("Sig_EtaQuark", "#eta_{#mu}", 100, -5, 5);
  CreateHisto("Sig_PhiQuark", "#phi_{#mu}@rad", 628, -3.14, 3.14);

  // combine muons from decay chain
  CreateHisto("Sig_MuMass", "m(#mu_{1},#mu_{2})@GeV", 1000, 0, 1000);
  CreateHisto("Sig_MuDPhi", "#Delta#phi(#mu_{1},#mu_{2})@rad", 315, 0., 3.15);
  CreateHisto("Sig_MuDEta", "#Delta#eta(#mu_{1},#mu_{2})", 100, 0., 10.);
  CreateHisto("Sig_MuDPt",  "#Deltap_{T}(#mu_{1},#mu_{2})@GeV", 1000, 0, 1000);
  CreateHisto("Sig_MuSumPt", "p_{T}(#mu_{1}) + p_{T}(#mu_{2})@GeV", 1000, 0, 1000);
  CreateHisto("Sig_MuPt", "p_{T}(#mu_{1}):p_{T}(#mu_{2})@GeV", 50, 0, 100, 50, 0, 100);

  // combine quarks from decay chain
  CreateHisto("Sig_JetMass", "m(j_{1},j_{2})@GeV", 500, 0, 500);
  CreateHisto("Sig_JetDPhi", "#Delta#phi(j_{1},j_{2})@rad", 315, 0., 3.15);
  CreateHisto("Sig_JetDEta", "#Delta#eta(j_{1},j_{2})", 100, 0., 10.);
  CreateHisto("Sig_JetDPt",  "#Deltap_{T}(j_{1},j_{2})@GeV", 500, 0, 500);
  CreateHisto("Sig_JetSumPt", "p_{T}(jet_{1}) + p_{T}(jet_{2})@GeV", 500, 0, 500);
  CreateHisto("Sig_JetPt", "p_{T}(q_{1}):p_{T}(q_{2})@GeV", 50, 0, 100, 50, 0, 100);

  // combine leptons and quarks
  CreateHisto("Sig_Mass3", "m(#tilde{#chi}^{0}_{1})@GeV", 500, 0, 500);
  CreateHisto("Sig_Angle3", "#Delta#phi(#mu, #tilde{#chi}^{0}_{1})@GeV", 315, 0, 3.15);
  CreateHisto("Sig_MuIsoR", "smallest #DeltaR between any jet and any muon", 100, 0, 5);

  // only for signal
  CreateHisto("Sig_LostMuEta", "Lost muons #eta_{#mu}", 100, -5, 5);
  CreateHisto("Sig_LostMuPt", "Lost muons p_{T}", 500, 0, 500);
  CreateHisto("Sig_LostMuPtEta", "Lost muons p_{T}@GeV:Lost muons #eta_{#mu}",
	      50, 0, 100, 50, -5, 5);

  // PF met rescaling for MC
  CreateHisto("pfmet_scaled", "scale factor f:PFMET@GeV", 500, 0, 500, 41, -0.005, 0.405);
  CreateHisto("JER_randomsf", "rescaling factors; random", 600, -2, 4);
  CreateHisto("JER_dptsf", "rescaling factors; from dpt", 600, -2, 4);
  CreateHisto("JER_deltae", "reco jet E - true jet E", 100, -50, 50);
  CreateHisto("JER_deltaepteta", "reco jet E - true jet E:p_{T}:#eta", 40, -100, 100, 200, 0, 1000, 5, 0, 2.5);
  CreateHisto("JER_scale", "(reco jet E - true jet E)/true jet E", 100, -2, 2);

  // Reskimming
  CreateHisto("bReweight_vtx_n", "Number of vertices in event", 100, -0.5, 99.5);
  CreateHisto("bSkim_muo_n", "Number of muons", 10, -0.5, 9.5);
  CreateHisto("bSkim_muo_pt0", "pt of leading muon", 500, 0, 500);
  CreateHisto("bSkim_muo_pt1", "pt of second muon", 500, 0, 500);
  CreateHisto("nSkim_muo_n", "Number of muons", 10, -0.5, 9.5);
  CreateHisto("nSkim_muo_pt0", "pt of leading muon", 500, 0, 500);
  CreateHisto("nSkim_muo_pt1", "pt of second muon (skimmer cuts)", 500, 0, 500);

  // Pileup reweighting BEFORE SKIM
  CreateHisto("bSkim_vtx_n", "Number of vertices in event", 100, -0.5, 99.5);
  CreateHisto("bSkim_vtx_ntr", "Number tracks in vertex", 151, -0.5, 150.5);
  CreateHisto("bSkim_vtx_ndof", "Number of deegrees of freedom of vertex", 151, -0.5, 150.5);
  CreateHisto("bSkim_vtx_x", "X position of vertex@cm", 100, -2, 2);
  CreateHisto("bSkim_vtx_y", "Y position of vertex@cm", 100, -2, 2);
  CreateHisto("bSkim_vtx_z", "Z position of vertex@cm", 120, -30, 30);

  // Pileup reweighting
  CreateHisto("vtx_n", "Number of vertices in event", 100, -0.5, 99.5);
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

  CreateHisto("delta_vtx_z", "z distance of vtx" , 10000, 0., 1.);
  CreateHisto("mumudz", "#Delta dz between 2 muons", 10000, 0, 1.);
  CreateHisto("last_mumudz", "#Delta dz between 2 muons", 10000, 0, 1.);

  // Object selection: electron
  CreateHisto("nElectron", "number of electrons", 10, 0, 10);
  CreateHisto("nElectron_pt", "Electron p_{T}@GeV", 1000, 0, 1000);
  CreateHisto("nElectron_DEta", "Electron #Delta #eta", 50, 0, 3);
  CreateHisto("nElectron_DPhi", "Electron #Delta #phi", 315, 0, 3.15);
  CreateHisto("nElectron_SigmaIEtaIEta", "Electron #sigma I#eta I#eta", 10, 0, 1);
  CreateHisto("nElectron_HoE", "Electron HCal / E", 100, 0, 0.2);
  CreateHisto("nElectron_1oE-1oP", "Electron 1/E - 1/p@GeV", 100, 0, 0.06);
  CreateHisto("nElectron_d0vtx", "Electron d0@cm", 100, 0, 0.025);
  CreateHisto("nElectron_dZvtx", "Electron dz@cm", 100, 0, 0.15);
  CreateHisto("nElectron_HasMatchedConversions", "Electrons with  matched conversions", 2, 0, 1.);
  CreateHisto("nElectron_MissingHits", "Electron missing hits", 2, 0, 2.);
  CreateHisto("nElectron_HCalIsoOpt", "Electron HCal Iso / pt", 100, 0, 0.4);
  CreateHisto("nElectron_ECalIsoOpt", "Electron ECal Iso / pt", 100, 0, 0.4);
  CreateHisto("nElectron_TkIsoOpt", "Electron Tk Iso / pt", 100, 0, 0.4);
  CreateHisto("nElectron_pfIsoOpt", "Electron PF Iso / pt", 100, 0, 0.3);

  // Object selection: muon
  CreateHisto("nMuon_pt", "#mu p_{T}@GeV", 1000, 0, 1000);
  CreateHisto("nMuon_eta","#eta_{#mu}", 100, -3, 3);
  CreateHisto("nMuon_id", "global muon flag", 2, 0, 2);
  CreateHisto("nMuon_ispf", "particle flow muon flag", 2, 0, 2);
  CreateHisto("nMuon_TrkChiNormCm", "#mu #chi^{2}/ndof", 100, 0, 100);
  CreateHisto("nMuon_hitsCm", "combined muon hits", 20, 0, 20);
  CreateHisto("nMuon_StationsMatched", "matched muon stations", 20, 0, 20);
  CreateHisto("nMuon_d0Tk", "#mu xy impact parameter wrt vertex from inner track@cm", 500, 0, 1);
  CreateHisto("nMuon_dzTk", "#mu z impact parameter wrt vertex from inner track@cm", 100, 0, 50);
  CreateHisto("nMuon_ChambersMatched", "matched muon chambers", 20, 0, 20);
  CreateHisto("nMuon_ValidTrackerHitsCm", "#mu tracker hits", 25, 0, 25);
  CreateHisto("nMuon_ValidPixelHitsCm", "#mu pixel hits", 5, -0.5, 4.5);
  CreateHisto("nMuon_TrackerLayersMeasCm", "#mu tracker layers with measurements", 14, -0.5, 13.5);
  CreateHisto("Muon_m", "m(#mu)@GeV", 100, 0, 1.);

  CreateHisto("Muon_nloose", "number of loose muons", 5, -0.5, 4.5);
  CreateHisto("Muon_ntight", "number of tight muons", 5, -0.5, 4.5);
  CreateHisto("Muon_loosetight", "number of tight muons:number of loose muons", 5, -0.5, 4.5, 5, -0.5, 4.5);
  // Object selection: tight/loose muons
  CreateHisto("nMuon_relIso", "muon relative isolation", 200, 0, 10);
  CreateHisto("relIsoOutJet", "muon relative isolation outside of jet", 200, 0, 10);
  CreateHisto("relIsoInJet", "muon relative isolation inside of jet", 200, 0, 10);
  CreateHisto("nMuon_matched", "muon trigger matching", 10, -0.5, 9.5);

  CreateHisto("matched", "muon trigger matching", 5, -0.5, 4.5);

  // Object selection: jets
  CreateHisto("nPfJet_pt", "jet p_{T}", 1000, 0, 1000);
  CreateHisto("nPfJet_eta", "#eta_{jet}", 100, -3, 3);
  CreateHisto("nPfJet_const", "number of jet constituents", 80, -0.5, 79.5);
  CreateHisto("nPfJet_pff_0", "jet energy fraction: charged hadrons", 100, 0, 1);
  CreateHisto("nPfJet_pff_1", "jet energy fraction: neutral hadrons", 100, 0, 1);
  CreateHisto("nPfJet_pff_2", "jet energy fraction: neutral EM hadrons", 100, 0, 1);
  CreateHisto("nPfJet_pff_3", "jet energy fraction: charged EM hadrons", 100, 0, 1);
  CreateHisto("nPfJet_pfn_0", "jet energy constitutens: charged hadrons", 40, -0.5, 39.5);
  CreateHisto("nPfJet_rmin", "min #Delta R(#mu, jet)", 500, 0, 5);
  CreateHisto("PfJet_m", "particle flow jet mass@GeV", 50, 0, 25);

  // Object selection: muon/jet arbitration
  CreateHisto("MJ_rmin", "min #Delta R(#mu, jet)", 100, 0, 5);
  CreateHisto("MJ_EF", "Energy fraction of muon in jet", 110, 0, 1.1);
  CreateHisto("MJ_EE", "E_{#mu}:E_{jet}", 20, 0, 100, 20, 0, 100);
  CreateHisto("MJ_mu_dR", "minimal #Delta R from muon to generated muon", 100, 0, 5);
  CreateHisto("MJ_jet_dR", "minimal #Delta R from muon to generated jet", 100, 0, 5);
  CreateHisto("MJ_dR", "minimal #Delta R from muon to generated jet:minimal #Delta R from muon to generated muon", 20, 0, 4, 20, 0, 4);

  // event cleaning
  CreateHisto("p_CSCHaloFilter", "p_CSCHaloFilter", 2, 0, 2);
  CreateHisto("p_HCALLaserFilter", "p_HCALLaserFilter", 2, 0, 2);
  CreateHisto("p_ecalLaserCorrFilter", "p_ecalLaserCorrFilter", 2, 0, 2);
  CreateHisto("p_HCALLaserFilterFromTriggerResult", "p_HCALLaserFilterFromTriggerResult", 2, 0, 2);
  CreateHisto("p_ECALDeadCellFilterTP", "p_ECALDeadCellFilterTP", 2, 0, 2);
  CreateHisto("p_BadSuperCrystalFilter", "p_BadSuperCrystalFilter", 2, 0, 2);
  CreateHisto("p_TrackingFailureFilter", "p_TrackingFailureFilter", 2, 0, 2);
  CreateHisto("p_HBHENoiseFilter", "p_HBHENoiseFilter", 2, 0, 2);
  CreateHisto("p_XtalsLargeLaser", "p_XtalsLargeLaser", 2, 0, 2);
  CreateHisto("HBHE_Noise", "HBHE noise filter result", 2, 0, 2);
  CreateHisto("noise_hcal_minE2Over10TS", "noise_hcal_minE2Over10TS", 2, 0, 2);
  CreateHisto("noise_hcal_maxE2Over10TS", "noise_hcal_maxE2Over10TS", 2, 0, 2);
  CreateHisto("noise_hcal_maxHPDHits", "noise_hcal_maxHPDHits", 2, 0, 2);
  CreateHisto("noise_hcal_maxRBXHits", "noise_hcal_maxRBXHits", 2, 0, 2);
  CreateHisto("noise_hcal_maxHPDNoOtherHits", "noise_hcal_maxHPDNoOtherHits", 2, 0, 2);
  CreateHisto("noise_hcal_maxZeros", "noise_hcal_maxZeros", 2, 0, 2);
  CreateHisto("noise_hcal_min25GeVHitTime", "noise_hcal_min25GeVHitTime", 2, 0, 2);
  CreateHisto("noise_hcal_max25GeVHitTime", "noise_hcal_max25GeVHitTime", 2, 0, 2);
  CreateHisto("noise_hcal_minRBXEMF", "noise_hcal_minRBXEMF", 2, 0, 2);
  CreateHisto("noise_hcal_numIsolatedNoiseChannels", "noise_hcal_numIsolatedNoiseChannels", 2, 0, 2);
  CreateHisto("noise_hcal_isolatedNoiseSumE", "noise_hcal_isolatedNoiseSumE", 2, 0, 2);
  CreateHisto("noise_hcal_isolatedNoiseSumEt", "noise_hcal_isolatedNoiseSumEt", 2, 0, 2);
  CreateHisto("noise_hcal_HasBadRBXTS4TS5", "noise_hcal_HasBadRBXTS4TS5", 2, 0, 2);

  

  // checks for MC production
  CreateHisto("check_nmuon", "number of muons", 5, -0.5, 4.5);
  CreateHisto("check_njets", "number of jets", 8, -0.5, 7.5);
  CreateHisto("check_mumass", "m(mu0, mu1)@GeV", 100, 0, 200);
  CreateHisto("check_vtxn", "number of vertices", 50, -0.5, 49.5);

  // Tight-to-loose ratio
  CreateHisto("TightMuons", "pt_{#mu} @GeV:#eta_{#mu}:Leading jet pt@GeV",
	      9, 10, 100, 5, 0, 2.5, 6, 40, 100);
  CreateHisto("LooseMuons", "pt_{#mu} @GeV:#eta_{#mu}:Leading jet pt@GeV",
	      9, 10, 100, 5, 0, 2.5, 6, 40, 100);

  CreateHisto("nTL_met", "MET@GeV", 100, 0, 500);
  CreateHisto("nTL_jetpt", "leading jet pt@GeV", 150, 0, 1500);
  CreateHisto("nTL_ht", "event HT@GeV", 150, 0, 1500);
  CreateHisto("nTL_nloose", "number of loose muons", 5, -0.5, 4.5);
  CreateHisto("nTL_firstmupt", "muon p_{T}@GeV", 100, 0, 100.);
  CreateHisto("nTL_mupt", "muon p_{T}@GeV", 100, 0, 100.);
  CreateHisto("nTL_jetdphi", "#Delta#phi between leading jet and loose muon", 315, 0, 3.15);
  CreateHisto("nTL_mt", "m_{T}(#mu, MET)@GeV", 500, 0, 500);
  CreateHisto("nTL_mumudz", "#Delta dz between 2 muons", 1000, 0, 1.);
  CreateHisto("nTL_zmass", "m(#mu^{+}, #mu^{-})@GeV", 500, 0, 500);
  CreateHisto("nTL_njets", "number of jets", 10, -0.5, 9.5);

  CreateHisto("TL_met", "MET@GeV", 100, 0, 500);
  CreateHisto("TL_metsig", "MET significance", 100, 0, 50);
  CreateHisto("TL_metdphi", "#Delta#phi between MET and loose muon", 315, 0, 3.15);
  CreateHisto("TL_jetpt", "leading jet pt@GeV", 150, 0, 1500);
  CreateHisto("TL_firstmupt", "leading muon pt@GeV", 100, 0, 100);
  CreateHisto("TL_jetdphi", "#Delta#phi between leading jet and loose muon", 315, 0, 3.15);
  CreateHisto("TL_ht", "event HT@GeV", 150, 0, 1500);
  CreateHisto("TL_nloose", "number of loose muons", 5, 0, 5);
  CreateHisto("TL_zmass", "m(#mu^{+}, #mu^{-})@GeV", 500, 0, 500);
  CreateHisto("TL_mumudz", "#Delta dz between 2 muons", 1000, 0, 1.);
  CreateHisto("TL_mt", "m_{T}(#mu, MET)@GeV", 500, 0, 500);
  CreateHisto("TL_st", "event ST@GeV", 100, 0, 1000);
  CreateHisto("TL_nmutight", "number of tight muons", 5, -0.5, 4.5);
  CreateHisto("TL_fakes", "fake rate requests:pt_{#mu} @GeV:#eta_{#mu}",
	      37, 15, 200, 6, 0, 3.0);


  // double-fake estimate
  CreateHisto("DF_pt", "p_{T}(#mu) (doublefake)@GeV", 5, 0, 50);
  CreateHisto("DF_eta", "#eta_{#mu} (doublefake)", 5, 0, 2.5);

  // single-fake estimate
  CreateHisto("SF_pt", "p_{T}(#mu) (doublefake)@GeV", 5, 0, 50);
  CreateHisto("SF_eta", "#eta_{#mu} (doublefake)", 5, 0, 2.5);

  // selection
  CreateHisto("Muon_pt0", "#mu p_{T}:#mu p_{T}@GeV", 1000, 0, 1000);
  CreateHisto("Muon_pt1", "#mu p_{T}:#mu p_{T}@GeV", 1000, 0, 1000);
  CreateHisto("Muon_pt", "#mu p_{T}:#mu p_{T}@GeV", 20, 0, 100, 20, 0, 100);
  CreateHisto("Jet_pt0", "jet p_{T}:jet p_{T}@GeV", 1000, 0, 1000);
  CreateHisto("Jet_pt1", "jet p_{T}:jet p_{T}@GeV", 1000, 0, 1000);
  CreateHisto("Jet_pt", "jet p_{T}:jet p_{T}@GeV", 20, 0, 100, 20, 0, 100);
  CreateHisto("pfmet_old", "Particle flow MET (before smearing)@GeV", 1000, 0, 1000);
  CreateHisto("pfmet", "Particle flow MET@GeV", 1000, 0, 1000);
  CreateHisto("met3_met6", "met6@GeV:met3@GeV", 40, 0, 200, 40, 0, 200);
  CreateHisto("met3_met7", "met7@GeV:met3@GeV", 40, 0, 200, 40, 0, 200);
  CreateHisto("met6_met7", "met7@GeV:met6@GeV", 40, 0, 200, 40, 0, 200);
  CreateHisto("Muon_charge", "muon charge:muon charge", 3, -1.5, 1.5, 3, -1.5, 1.5);
  CreateHisto("Muon_ch", "muon charge:muon charge", 2, -1.5, 0.5);

  // create individual histograms
  CreateHisto("DeltaPhi", "#Delta#phi(#mu_{1}, gaugino)", 315, 0., 3.15);
  CreateHisto("m_mumu", "m(#mu^{+}, #mu^{-})@GeV", 500, 0, 1000);
  CreateHisto("m_mumu_precut", "m(#mu^{+}, #mu^{-})@GeV", 500, 0, 1000);
  CreateHisto("m_gaugino", "gaugino mass m(#mu_{1},j_{1},j_{2})", 25, 0, 1000);
  CreateHisto("m_gaugino_precut", "gaugino mass m(#mu_{1},j_{1},j_{2})", 25, 0, 1000);
  CreateHisto("m_smuon", "smuon mass m(#mu_{0},#mu_{1},j_{1},j_{2})", 100, 0, 2000);
  CreateHisto("m_smuon_precut", "smuon mass m(#mu_{0},#mu_{1},j_{1},j_{2})", 100, 0, 2000);
  CreateHisto("muo_n_precut", "Number of muons", 20, -0.5, 19.5);
  CreateHisto("m_mumu_zpeak", "Z-Peak m(#mu^{+}, #mu^{-})@GeV", 500, 0, 1000);
  CreateHisto("CR1_m_mumu", "CR1 m(#mu^{+}, #mu^{-})@GeV", 500, 0, 1000);
  CreateHisto("CR2_m_mumu", "CR2 m(#mu^{+}, #mu^{-})@GeV", 500, 0, 1000);
  CreateHisto("CR3_m_mumu", "CR3 m(#mu^{+}, #mu^{-})@GeV", 500, 0, 1000);
  CreateHisto("CR4_m_mumu", "CR4 m(#mu^{+}, #mu^{-})@GeV", 500, 0, 1000);
  CreateHisto("CR4_m_gaugino", "gaugino mass m(#mu_{1},j_{1},j_{2})", 25, 0, 1000);
  CreateHisto("CR4_m_smuon", "smuon mass m(#mu_{0},#mu_{1},j_{1},j_{2})", 100, 0, 2000);
  CreateHisto("CR5_m_mumu", "CR5 m(#mu^{+}, #mu^{-})@GeV", 500, 0, 1000);
  CreateHisto("CR5_m_gaugino", "gaugino mass m(#mu_{1},j_{1},j_{2})", 25, 0, 1000);
  CreateHisto("CR5_m_smuon", "smuon mass m(#mu_{0},#mu_{1},j_{1},j_{2})", 100, 0, 2000);
  CreateHisto("CR6_m_mumu", "CR6 m(#mu^{+}, #mu^{-})@GeV", 500, 0, 1000);
  CreateHisto("CR6_m_gaugino", "gaugino mass m(#mu_{1},j_{1},j_{2})", 25, 0, 1000);
  CreateHisto("CR6_m_smuon", "smuon mass m(#mu_{0},#mu_{1},j_{1},j_{2})", 100, 0, 2000);

  const double binsb[] = { 0, 3.3, 100 };
  CreateHisto("isbtag", "isbtag:TCHE btag", 2, binsb, 2, 0, 2);
  const double bins[] = { 0, 300, 700, 5000 };
  const int nMax = sizeof(bins)/sizeof(double)-1;
  CreateHisto("jjmm_m", "smuon mass m(#mu_{0},#mu_{1},jets)", nMax, bins, nMax, bins);
  CreateHisto("m_smu_chi", "m(#chi): m(#tilde{#mu})", 110, 0, 2200, 25, 0, 500);
  CreateHisto("muo_n", "Number of muons", 30, 0, 30);
  CreateHisto("ht", "HT@GeV", 100, 0, 500);

  // after btag cut
  CreateHisto("btag_m_mumu", "m(#mu^{+}, #mu^{-})@GeV", 500, 0, 1000);
  CreateHisto("btag_m_gaugino", "gaugino mass m(#mu_{1},j_{1},j_{2})", 25, 0, 1000);
  CreateHisto("btag_m_smuon", "smuon mass m(#mu_{0},#mu_{1},j_{1},j_{2})", 100, 0, 2000);
  CreateHisto("btag_m_smu_chi", "m(#chi): m(#tilde{#mu})", 110, 0, 2200, 25, 0, 500);
  CreateHisto("btag_ht", "HT@GeV", 100, 0, 500);
  CreateHisto("btag_jjmm_m", "smuon mass m(#mu_{0},#mu_{1},jets)", nMax, bins, nMax, bins);

  // event weights
  CreateHisto("weight", "weight", 200, 0., 10.);
  CreateHisto("log_global_weight", "log(global_weight)", 100, -5, 5);
  CreateHisto("log_weight", "log(weight)", 100, -5, 5);
  
  // signalstudy
  CreateHisto("Sig_signal_counter", "Signal Counter", 10, 0., 10.);
  CreateHisto("Cut_signal_counter", "Signal Counter", 10, 0., 10.);
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

void Analysis::CreateBranches()
{
  // signal study, MC truth particles
  fOutputTree.Branch("fSigMu0", & fSigMu0, 32000, 0);
  fOutputTree.Branch("fSigMu1", & fSigMu1, 32000, 0);
  fOutputTree.Branch("fSigJet0", & fSigJet0, 32000, 0);
  fOutputTree.Branch("fSigJet1", & fSigJet1, 32000, 0);

  // MET
  fOutputTree.Branch("fJER_met_old", & fJER_met_old, "fJER_met_old/D");

  // fake rate
  fOutputTree.Branch("fFakeRate", & fFakeRate, "fFakeRate[2]/D");
  fOutputTree.Branch("fSingleFakeWeight", & fSingleFakeWeight, "fSingleFakeWeight/D");

  // selected particles
  fOutputTree.Branch("fMuoId", & fMuoId, "fMuoId[2]/I");
  fOutputTree.Branch("fMuon0", & fMuon0, 32000, 0);
  fOutputTree.Branch("fMuon1", & fMuon1, 32000, 0);
  fOutputTree.Branch("fJets", & fJets, "fJets/I");
  fOutputTree.Branch("fJetId", & fJetId, "fJetId[fJets]/I");
  fOutputTree.Branch("fJet0", & fJet0, 32000, 0);
  fOutputTree.Branch("fJet1", & fJet1, 32000, 0);
  fOutputTree.Branch("fIsBTagged", & fIsBTagged, "fIsBTagged/O");
}

void Analysis::CreateHisto(const char * name, const char * title, Int_t nbinsx, Double_t xlow, Double_t xup)
{
  histo[name] = new TH1D(Form("h1_0_%s", name), title, nbinsx, xlow, xup);
}

void Analysis::CreateHisto(const char * name, const char * title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup)
{
  histo2[name] = new TH2D(Form("h2_0_%s", name), title, nbinsx, xlow, xup, nbinsy, ylow, yup);
}

void Analysis::CreateHisto(const char * name, const char * title, Int_t nbinsx, const Double_t * xbins, Int_t nbinsy, Double_t ylow, Double_t yup)
{
  histo2[name] = new TH2D(Form("h2_0_%s", name), title, nbinsx, xbins, nbinsy, ylow, yup);
}

void Analysis::CreateHisto(const char * name, const char * title, Int_t nbinsx, const Double_t * xbins, Int_t nbinsy, const Double_t * ybins)
{
  histo2[name] = new TH2D(Form("h2_0_%s", name), title, nbinsx, xbins, nbinsy, ybins);
}

void Analysis::CreateHisto(const char * name, const char * title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, Int_t nbinsz, Double_t zlow, Double_t zup)
{
  histo3[name] = new TH3D(Form("h3_0_%s", name), title, nbinsx, xlow, xup, nbinsy, ylow, yup, nbinsz, zlow, zup);
}

void Analysis::Fill(const char * name, double value)
{
  TH1D * h = histo[name];
  if (h != 0)
    h->Fill(value, global_weight);
  else {
    THROW(std::string("Histogram \"") + name + std::string("\" not existing. Did you misspell or forgot to create?"));
  }
}

void Analysis::Fill(const char * name, const char * bin)
{
  TH1D * h = histo[name];
  if (h != 0) {
    h->Fill(bin, global_weight);
  }
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
    histo2[name]->Fill(x, y, global_weight);
  else {
    THROW(std::string("Histogram \"") + name + std::string("\" not existing. Did you misspell or forgot to create?"));
  }
}

void Analysis::Fill(const char * name, double x, double y, double z)
{
  TH3D * h = histo3[name];
  if (h != 0)
    histo3[name]->Fill(x, y, z, global_weight);
  else {
    THROW(std::string("Histogram \"") + name + std::string("\" not existing. Did you misspell or forgot to create?"));
  }
}

double Analysis::GetFakeRate(double muopt, double eta, double jetpt)
{
  double fakeRate = 0.;
  Fill("TL_fakes", muopt, eta);
  // check valid input
  if (muopt < 15) {
    // this can happen, but should not be bad because after the fake rate is
    // applied, muons are required to exceed at least 15 GeV in a later
    // analysis step. So the choice below should not matter
    if (fFakeRateMethod == "lastbin") {
      muopt = 15.001;
    }
    else if (fFakeRateMethod == "zero") {
      return 0;
    }
    else
      THROW("unknown fake rate method: " + fFakeRateMethod);
  }
  if (muopt > 100) {
    if (fFakeRateMethod == "lastbin") {
      muopt = 99.99;
    }
    else if (fFakeRateMethod == "zero") {
      return 0;
    }
    else
      THROW("unknown fake rate method: " + fFakeRateMethod);
  }
  if (fabs(eta) > 2.5) {
    ERROR("muon eta > 2.5 in Analysis::GetFakeRate() - this must not happen");
    if (fFakeRateMethod == "lastbin") {
      eta = 2.499;
    }
    else if (fFakeRateMethod == "zero") {
      return 0;
    }
    else
      THROW("unknown fake rate method: " + fFakeRateMethod);
  }
  if (fFakeRateDimensions == 3) {
    Int_t bin = fFakeRateHisto3D->FindFixBin(muopt, fabs(eta), jetpt);
    fakeRate = fFakeRateHisto3D->GetBinContent(bin);
  }
  else if (fFakeRateDimensions == 2) {
    Int_t bin = fFakeRateHisto2D->FindFixBin(muopt, fabs(eta));
    fakeRate = fFakeRateHisto2D->GetBinContent(bin);
  }
  if (fakeRate <= 0. || fakeRate >= 1.) {
    WARNING("unexpected value " << fakeRate << " for T/L ratio in Analysis::GetFakeRate()");
    WARNING("muo pt= " << muopt << ", eta = " << eta << ", jetpt = " << jetpt);
    if (fakeRate < 0)
      //THROW("T/L ratio equals to or less than 0, cannot proceed");
      return 0.001;
    if (fakeRate >= 1)
      THROW("T/L ratio equals to or larger than 1, cannot proceed");
  }
  return fakeRate;
}

// helper
Double_t Analysis::DeltaPhi(double a, double b) {

  double temp = fabs(a-b);
  if (temp <= TMath::Pi())
    return temp;
  else
    return  2.*TMath::Pi() - temp;
}

double Analysis::GetJERScale(double eta)
{
  eta = TMath::Abs(eta);
  // values taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
  if (eta < 0.5) {
    return 1.052; //+-0.012+0.062-0.061
  }
  else if (eta < 1.1) {
    return 1.057; //+-0.012+0.056-0.055
  }
  else if (eta < 1.7) {
    return 1.096; //+-0.017+0.063-0.062
  }
  else if (eta < 2.3) {
    return 1.134; //+-0.035+0.087-0.085
  }
  else if (eta < 5.0) {
    return 1.288; //+-0.127+0.155-0.153
  }
  else {
    DEBUG("jet eta = " << eta << " out of range 0..5 in Analysis::GetJERScale");
    // return some value
    return 1.288;
  }
}

// Jet smearing: match jets to GenJets and smear them. Propagate into MET
void Analysis::PFJetSmearing()
{
  // save old value for later use
  fJER_met_old = met_et[MET_INDEX];
  // do not smear data
  if (fInputType == "data" || fJER_calculation)
    return;

  // loop over all jets
  for (int j = 0; j < pfjet_n; j++) {
    // correct only jets with pT > 15 GeV
    if (pfjet_pt[j] > 15.) {
      // match to generated jet
      double sum_truth = 0;
      for (int t = 0; t < truthjet_n; t++) {
        double dR = sqrt(pow(pfjet_eta[j]-truthjet_eta[t], 2) +
			 pow(DeltaPhi(pfjet_phi[j],truthjet_phi[t]), 2));
        if (dR < 0.5) {
	  // adding up all the jet components
	  sum_truth += truthjet_pt[t];
        }
      }

      double dpt = 0;
      double scale = 0;
      if (sum_truth == 0) {
	// no truth jets found
	dpt = gRandom->Gaus(fJER_center, TMath::Sqrt( TMath::Power(GetJERScale(pfjet_eta[j]),2) - 1 ) * fJER_smear);
	scale = 1 + dpt/pfjet_pt[j];
	Fill("JER_randomsf", scale);
	DEBUG("Random: scale = " << scale);
      }
      else {
	
	scale = (sum_truth + GetJERScale(pfjet_eta[j]) * (pfjet_pt[j] - sum_truth)) / pfjet_pt[j];
	Fill("JER_dptsf", scale);
	// if (truth_jet != pfjet_truth[j])
	//   WARNING("Calculated true jet: " << truth_jet << ", from SUSYAna: " << pfjet_truth[j]);
	DEBUG("Scaled: dpt = " << dpt << ", scale = " << scale);
      }
      // correct MET
      met_ex[MET_INDEX]   -= pfjet_px[j] * (scale - 1.);
      met_ey[MET_INDEX]   -= pfjet_py[j] * (scale - 1.);

      // correct PF jets
      pfjet_E[j]  *= scale;
      pfjet_Et[j] *= scale;
      pfjet_p[j]  *= scale;
      pfjet_pt[j] *= scale;
      pfjet_px[j] *= scale;
      pfjet_py[j] *= scale;
      pfjet_pz[j] *= scale;
    }
  }
  
  // recalculate MET
  TVector3 mpf(met_ex[MET_INDEX], met_ey[MET_INDEX], 0.);
  met_phi[MET_INDEX] = mpf.Phi();
  met_et[MET_INDEX]  = mpf.Perp();
}

void Analysis::PFJetSmearingCalculation()
{
  TH2D * h = histo2["pfmet_scaled"];
  if (h == 0)
    THROW("histogram pfmet_scaled not found");

  if (fInputType == "data" || !fJER_calculation )
    return;

  // loop over all bins
  for (int biny = 1; biny < h->GetNbinsY()+1; biny++) {
    // get value at bin center
    //double f = h->GetYaxis()->GetBinLowEdge(biny)+h->GetYaxis()->GetBinWidth(biny)/2.;

    // loop over all jets
    for (int j = 0; j < pfjet_n; j++) {
      // correct only jets with pT > 15 GeV
      if (pfjet_pt[j] > 15.) {
	// match to generated jet
	// int truth_jet = -1;
	double sum_truth = 0;
	for (int t = 0; t < truthjet_n; t++) {
	  double dR = sqrt(pow(pfjet_eta[j]-truthjet_eta[t], 2) +
			   pow(DeltaPhi(pfjet_phi[j],truthjet_phi[t]), 2));

	  if (dR < 0.5) {
	    // if (dpT > 0)
	    //   WARNING("Analysis::PFJetSmearing() found more than one matching truth jet");
	    sum_truth += truthjet_pt[t];
	    // truth_jet = t;
	  }
	}

	double dpt = 0;
	double scale = 0;
	if (sum_truth == 0) {
	  // no truth jets found, so no calculation to be done
	}
	else {
	  dpt = pfjet_pt[j] - sum_truth;
	  scale = dpt/pfjet_pt[j];
	  // if (truth_jet != pfjet_truth[j])
	  // 	WARNING("Calculated true jet: " << truth_jet << ", from SUSYAna: " << pfjet_truth[j]);
	  Fill("JER_deltae", dpt);
	  Fill("JER_deltaepteta", dpt, pfjet_pt[j], pfjet_eta[j]);
	  Fill("JER_scale", scale);
	}
      }
    }
  }
}

