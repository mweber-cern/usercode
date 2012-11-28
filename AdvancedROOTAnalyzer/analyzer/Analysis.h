#ifndef __Analyis_h
#define __Analyis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TEnv.h>
#include <TString.h>

#include <algorithm>
#include <iostream>
#include <vector>
#include <set>
#include <map>

#include "BTagging.h"

#if VERSION == 78
#include "TreeContent_v78.h"
#include "LumiReweightingStandAlone_v78.h"
#elif VERSION == 88
#include "TreeContent_v88.h"
#include "LumiReweightingStandAlone_v88.h"
#elif VERSION == 92
#include "TreeContent_v92.h"
#include "LumiReweightingStandAlone_v92.h"
#elif VERSION == 96
#include "TreeContent_v96.h"
#include "LumiReweightingStandAlone_v96.h"
#endif


using namespace std;

class Analysis : public TreeContent {
protected:
  TTree &   fInputTree;
  TTree &   fOutputTree;
  TEnv  &   fCfgFile;

  //////////////////////////////////////////////////////////////////////
  // configuration options
  bool      fFill;      // fill output tree?
  string    fInputType; // type of input data
  string    fSample;    // the sample name
  string    fAnalysisType; // type of analysis to be performed
  Long64_t  fMaxEvents; // max events to process
  Long64_t  fMaxTreeSize; // maximum tree size (output)
  bool      fFindDuplicates; // find duplicate events?
  bool      fDumpAll;   // Dump all event information
  bool      fDumpBasic; // Dump basic information
  bool      fDumpTruth; // Dump MC truth information
  string    fDumpTrigger; // Dump Trigger information
  bool      fPileupReweighting; // use pileup reweighting?
  bool      fSkimActive; // Skimming Active?
  Int_t     fSkimMuons;  // how many muons to require
  Double_t  fSkimMuoptfirst; // cut on muon with highest pt
  Double_t  fSkimMuoptother; // cut on muon with second highest pt
  Double_t  fLooseMuonRelIso; // relative isolation cut for loose muon
  string    fFakeRateMethod; // fake rate method, "lastbin" or "zero"
  Int_t     fFakeRateDimensions; // number of dimensions for fake rate calculation
  vector<string> fTrigger; // trigger selection
  bool      fForceUnprescaledTriggers; // Force unprescaled triggers
  string    fLumiRanges; // lumi ranges (JSON file) for data processing
  Double_t  fDeltaPhiMin; // max delta phi (leading mu, gaugino)

  // cuts for T/L ratio
  Double_t fTL_met_max;
  Double_t fTL_ht_min;
  Double_t fTL_jetpt_min;
  Double_t fTL_nloose_min;
  Double_t fTL_nloose_max;
  Double_t fTL_zmass_min;
  Double_t fTL_zmass_max;
  Double_t fTL_mupt_min;
  Double_t fTL_jetdphi_min;
  Double_t fTL_mt_max;
  Double_t fTL_njets_min;

  // values for smearing jet energies (JER)
  Double_t fJER_scale;
  Double_t fJER_center;
  Double_t fJER_smear;
  Double_t fJER_met_old;

  //////////////////////////////////////////////////////////////////////
  // analysis variables

  // signal study
  bool             fIsSignal;
  TLorentzVector * fSigMu0; // only in case of signal: generated muons
  TLorentzVector * fSigMu1;
  TLorentzVector * fSigJet0; // only in case of signal: generated jets
  TLorentzVector * fSigJet1; 

  // fake rate
  TH3D           * fFakeRateHisto3D;
  TH2D           * fFakeRateHisto2D;
  Double_t         fSingleFakeWeight;  // fake rate weight for one mu fluctuation 
  Double_t         fFakeRate[2];       // fake rate for the two muons

  // save particles for later analysis
  Int_t            fMuoId[2]; // ID's of selected muons
  TLorentzVector * fMuon0;
  TLorentzVector * fMuon1;
  Int_t            fJets;
  Int_t            fJetId[40]; // ID's of selected main jets
  TLorentzVector * fJet0;
  TLorentzVector * fJet1;
  double           fBtag0;
  double           fBtag1;
  Bool_t           fIsBTagged0;
  Bool_t           fIsBTagged1;
  TLorentzVector * fSmuon;
  TLorentzVector * fGaugino;

public:
  Analysis(TTree & inputTree, TTree & outputTree, TEnv & cfgFile);
  virtual ~Analysis();

  void SetBranchAddresses();
  void CreateBranches();
  void CreateVariables();
  void Loop();

protected:
  //////////////////////////////////////////////////////////////////////
  // analysis functions
  void ReconstructFinalState(map<char, vector<int> > & particles, int vertex);
  void SignalStudy();
  void TightLooseRatioCalculation(const vector<int> & loose_muons,
				  const vector<int> & tight_muons,
				  const vector<int> & jets,
				  const double HT);
  double GetFakeRate(double muopt, double eta, double jetpt);
  bool filterHBHENoise();
  void PFJetSmearing();
  void PFJetSmearingCalculation();

  // helper functions
  void CreateHistograms();
  void CreateHisto(const char * name, const char * title, 
		   Int_t nbinsx, Double_t xlow, Double_t xup);
  void CreateHisto(const char * name, const char * title, 
		   Int_t nbinsx, Double_t xlow, Double_t xup, 
		   Int_t nbinsy, Double_t ylow, Double_t yup);
  void CreateHisto(const char * name, const char * title, 
		   Int_t nbinsx, const Double_t * xbins, 
		   Int_t nbinsy, Double_t ylow, Double_t yup);
  void CreateHisto(const char * name, const char * title, 
		   Int_t nbinsx, const Double_t * xbins, 
		   Int_t nbinsy, const Double_t * ybins);
  void CreateHisto(const char * name, const char * title, 
		   Int_t nbinsx, Double_t xlow, Double_t xup, 
		   Int_t nbinsy, Double_t ylow, Double_t yup, 
		   Int_t nbinsz, Double_t zlow, Double_t zup);
  void Fill(const char * name, double value);
  void Fill(const char * name, const char * bin);
  void FillNoWeight(const char * name, double value);
  void Fill(const char * name, double x, double y);
  void Fill(const char * name, double x, double y, double z);

  // taken over from ACSUSYAna
  void BasicDump(int i);
  void TriggerDump(TString sel);
  void MuonDump(bool full=1);
  void CaloJetDump();
  void PFJetDump();
  void TruthJetDump();
  void TruthDump();
  void VertexDump();
  void METDump();
  void SCDump();
  void EleDump(bool full=1);
  void PFEleDump(bool full=1);

  bool FindDuplicates(int run, int evt, double x1, double x2);
  Double_t DeltaPhi(double a, double b);
  
  typedef std::pair< pair<int,int> , pair<double,double> > Key;
  typedef std::set<Key> KeySet;
  typedef KeySet::const_iterator KeyIter;
  KeySet _keys;

  // this is the magic for pileup reweighting
  reweight::LumiReWeighting LumiWeights_;
  /* reweight::PoissonMeanShifter PShiftUp_; */
  /* reweight::PoissonMeanShifter PShiftDown_; */

  // b-tagging correction
  BTagging fBTagging;

  // histogram store
  map<string, TH1D *> histo;
  map<string, TH2D *> histo2;
  map<string, TH3D *> histo3;

  friend class truth_pt_comparator;
};

class truth_pt_comparator 
{
protected:
  const Analysis & fAnalysis;
public:
  truth_pt_comparator(const Analysis & analysis);
  bool operator()(int i, int j);
};

#endif
