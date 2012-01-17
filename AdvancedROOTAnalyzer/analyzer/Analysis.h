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

#include "TreeContent.h"
#include "LumiReweightingStandAlone.h"

using namespace std;

class Analysis : public TreeContent {
protected:
  TTree &   fInputTree;
  TTree &   fOutputTree;
  TEnv  &   fCfgFile;
  double    fWeight; // event weight

  //////////////////////////////////////////////////////////////////////
  // configuration options
  bool      fFill;      // fill output tree?
  string    fType;      // type of input data
  Long64_t  fMaxEvents; // max events to process
  Long64_t  fMaxTreeSize; // maximum tree size (output)
  bool      fFindDuplicates; // find duplicate events?
  bool      fDumpAll;   // Dump all event information
  bool      fDumpTruth; // Dump MC truth information
  bool      fSkimActive; // Skimming Active?
  Int_t     fSkimMuons;  // how many muons to require
  Double_t  fSkimMuoptfirst; // cut on muon with highest pt
  Double_t  fSkimMuoptother; // cut on muon with second highest pt
  const char * fPileupDataFile; 
  const char * fPileupMCFile;
  const char * fPileupWeightFile;

  //////////////////////////////////////////////////////////////////////
  // analysis variables
  bool      fIsSignal;
  TString * fTrigger;
  bool      fPileupInitFromWeights;

public:
  Analysis(TTree & inputTree, TTree & outputTree, TEnv & cfgFile);
  virtual ~Analysis();

  void SetBranchAddresses();
  void Loop();

protected:
  // own functions
  void CreateHistograms();
  void ReconstructFinalState(map<char, vector<int> > & particles, int vertex);
  void CreateHisto(const char * name, const char * title, 
		   Int_t nbinsx, Double_t xlow, Double_t xup);
  void CreateHisto(const char * name, const char * title, 
		   Int_t nbinsx, Double_t xlow, Double_t xup, 
		   Int_t nbinsy, Double_t ylow, Double_t yup);
  void Fill(const char * name, double value);
  void FillNoWeight(const char * name, double value);
  void Fill(const char * name, double x, double y);
  void SignalStudy();

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

  void doJESandRecalculateMET(TString corr);

  bool FindDuplicates(int run, int evt, double x1, double x2);
  
  Double_t DeltaPhi(double a, double b);
  Double_t mT(double et1, double phi1, double et2, double phi2);

  Double_t AlphaT(const std::vector<TLorentzVector> & objects);
  Double_t MinDEt(const std::vector<TLorentzVector> & objects, 
		  std::vector<UInt_t> * lista = NULL, 
		  std::vector<UInt_t> * listb = NULL);
  Double_t SumET(const std::vector<TLorentzVector> & objects);
  Double_t MT(const std::vector<TLorentzVector> & objects);
  double M2(double E1, double px1, double py1, double pz1, double E2, double px2, double py2, double pz2);
  double M2v2(double pt1, double eta1, double phi1, double m1, double pt2, double eta2, double phi2, double m2);
  double AngleMuons(double pt1, double eta1, double phi1, double m1, double pt2, double eta2, double phi2, double m2);

  
  typedef std::pair< pair<int,int> , pair<double,double> > Key;
  typedef std::set<Key> KeySet;
  typedef KeySet::const_iterator KeyIter;
  KeySet _keys;

  //this is the magic for reweighting
  reweight::LumiReWeighting LumiWeights_;
  reweight::PoissonMeanShifter PShiftUp_;
  reweight::PoissonMeanShifter PShiftDown_;

  // histogram store
  map<string, TH1D *> histo;
  map<string, TH2D *> histo2;

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
