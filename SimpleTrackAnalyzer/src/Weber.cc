// -*- C++ -*-
//
// Package:    Weber
// Class:      Weber
// 
/**\class Weber Weber.cc Martin/Weber/src/Weber.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Martin Weber
//         Created:  Thu Mar 20 21:07:21 CET 2008
// $Id: Weber.cc,v 1.1 2008/05/20 07:34:20 mweber Exp $
//
//


#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
// system include files
#include <memory>
#include <iostream>
using namespace std;

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
 
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//
// class decleration
//

class Weber : public edm::EDAnalyzer {
public:
  explicit Weber(const edm::ParameterSet&);
  ~Weber();


private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  TH1F * hfx;
  TH1F * hfy;
  TH1F * hfz;
  TH1F * hlx;
  TH1F * hly;
  TH1F * hlz;
  TH1F * hp;
  TH1F * hpt;
  TH1F * hphi;
  TH1F * heta;
  TH1F * hrho;
  TH1F * h1;
  TH2F * hfxy;
  TH2F * hlxy;
  TH2F * hxy;
  TH2F * hrz;
  TH2F * hxz;
  TH2F * hhits;
  TFile * myFile;

  // ---------configuration----------------------------
  const edm::InputTag inputTag_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Weber::Weber(const edm::ParameterSet& iConfig) :
  inputTag_ ( iConfig.getParameter<edm::InputTag>("src") )
{

  // now do what ever initialization is needed
  myFile = new TFile("histo.root", "RECREATE");
  hfx = new TH1F("hfx", "x", 100, -120, 120);
  hfy = new TH1F("hfy", "y", 100, -120, 120);
  hfz = new TH1F("hfz", "z", 100, -300, 300);
  hlx = new TH1F("hlx", "x", 100, -120, 120);
  hly = new TH1F("hly", "y", 100, -120, 120);
  hlz = new TH1F("hlz", "z", 100, -300, 300);
  hp = new TH1F("hp", "p", 100, 0, 100);
  hpt = new TH1F("hpt", "pt", 100, 0, 100);
  hphi = new TH1F("hphi", "#phi", 100, -1.9, 0);
  heta = new TH1F("heta", "#eta", 100, -3, 3);
  hrho = new TH1F("hrho", "#rho", 100, 0, 1000);
  h1 = new TH1F("h1", "r", 200, -120, 120);
  hfxy = new TH2F("hfxy", "fxy", 100, -120, 120, 100, -120, 120);
  hlxy = new TH2F("hlxy", "lxy", 100, -120, 120, 100, -120, 120);
  hxy = new TH2F("hxy", "xy", 100, -120, 120, 100, -120, 120);
  hrz = new TH2F("hrz", "rz", 100, -300, 300, 100, -120, 120);
  hxz = new TH2F("hxz", "xz", 100, -300, 300, 100, -120, 120);
  hhits = new TH2F("hhits", "hits x vs y", 100, -120, 120, 100, -120, 120);
}


Weber::~Weber()
{
    // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  myFile->cd();
  hfx->Write();
  hfy->Write();
  hfz->Write();
  hlx->Write();
  hly->Write();
  hlz->Write();
  hp->Write();
  hpt->Write();
  hphi->Write();
  heta->Write();
  hrho->Write();
  h1->Write();
  hfxy->Write();
  hlxy->Write();
  hxy->Write();
  hrz->Write();
  hxz->Write();
  myFile->Close();
  delete myFile;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
Weber::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using reco::TrackCollection;

   Handle<TrackCollection> tracks;
   iEvent.getByLabel(inputTag_,tracks);
   for(TrackCollection::const_iterator itTrack = tracks->begin();
       itTrack != tracks->end();
       ++itTrack) {
     const reco::TrackBase::Point firstPoint(itTrack->innerPosition());     
     double fx = firstPoint.x();
     double fy = firstPoint.y();
     double fz = firstPoint.z();
     double r = sqrt(fx*fx+fy*fy);
     const reco::TrackBase::Point lastPoint(itTrack->outerPosition());     
     double lx = lastPoint.x();
     double ly = lastPoint.y();
     double lz = lastPoint.z();
     double px = (*itTrack).momentum().x();
     double py = (*itTrack).momentum().y();
     double pz = (*itTrack).momentum().z();
     double p = sqrt(px*px+py*py+pz*pz);
     double pt = (*itTrack).pt();
     double phi = (*itTrack).phi();
     double eta = (*itTrack).eta();
     double rho = firstPoint.Rho();
     hfx->Fill(fx);
     hfy->Fill(fy);
     hfz->Fill(fz);
     hlx->Fill(lx);
     hly->Fill(ly);
     hlz->Fill(lz);
     hp->Fill(p);
     hpt->Fill(pt);
     hphi->Fill(phi);
     heta->Fill(eta);
     hrho->Fill(rho);
     h1->Fill(r);
     hfxy->Fill(fx, fy);
     hlxy->Fill(lx, ly);
     hxy->Fill(fx, fy);
     hxy->Fill(lx, ly);
     hrz->Fill(fz, r);
     hxz->Fill(fz, fx);
     for (trackingRecHit_iterator rit = (*itTrack).recHitsBegin();
	  rit != (*itTrack).recHitsEnd(); ++rit) {
       TrackingRecHitRef hit = *rit;
     }
   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
Weber::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Weber::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(Weber);
