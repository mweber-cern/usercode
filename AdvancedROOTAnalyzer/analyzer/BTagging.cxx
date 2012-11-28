#include "BTagging.h"

#include "Utilities.h"

// Choose Track Counting High Efficiency with working point medium (TCHEM)
const double btag_cut = 3.3;

BTagging::BTagging(int seed)
{
  // initialize private random number generator with given seed
  rand_ = new TRandom3(seed);
}

BTagging::~BTagging()
{
  // delete private random number generator
  delete rand_;
}

bool BTagging::isBJet(double btag, int pdgIdPart, double pt, double eta)
{
  eta = TMath::Abs(eta);
  pdgIdPart = TMath::Abs(pdgIdPart);
  bool isBTagged = (btag > btag_cut);
  if (pt < 30 || pt > 670 || eta > 2.4) {
    WARNING("jet properties out of allowed range in BTagging::isBJet(): pt = " << pt 
	    << ", eta = " << eta);
    return isBTagged;
  }
  modifyBTagsWithSF(isBTagged, 
		    pdgIdPart, 
		    GetBTagScaleFactorTCHEM(btag),
		    GetBTagEfficiencyMCTCHEM(btag),
		    GetCTagEfficiencyMCTCHEM(btag),
		    GetBMisTagScaleFactorTCHEM(pt, eta),
		    GetBMisTagEfficiencyTCHEM(pt, eta)
    );
  return isBTagged;
}

double BTagging::GetBTagScaleFactorTCHEM(double pt)
{
  return 0.932251*((1.+(0.00335634*pt))/(1.+(0.00305994*pt)));
}

double BTagging::GetBTagScaleFactorErrorTCHEM(double pt)
{
  double ptmin[] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500};
  double ptmax[] = {40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 670};
  double SFb_error[] = {
    0.0311456,
    0.0303825,
    0.0209488,
    0.0216987,
    0.0227149,
    0.0260294,
    0.0205766,
    0.0227065,
    0.0260481,
    0.0278001,
    0.0295361,
    0.0306555,
    0.0367805,
    0.0527368 
  };
  for (unsigned int i = 0; i < sizeof(ptmin)/sizeof(double); i++) {
    if (pt >= ptmin[i] && pt < ptmax[i]) {
      return SFb_error[i];
    }
  }
  WARNING("pt = " << pt << " out of range " << ptmin[0] << "..." << ptmax[13]);
  return 0;
}

double BTagging::GetBMisTagScaleFactorTCHEM(double pt, double eta, bool finebin)
{
  // from https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFlightFuncs.C
  eta = TMath::Abs(eta);
  if (finebin) {
    if (eta < 0.8) {
      return (1.2875*((1+(-0.000356371*pt))+(1.08081e-07*(pt*pt))))+(-6.89998e-11*(pt*(pt*(pt/(1+(-0.0012139*pt))))));
    }
    else if (eta < 1.6) {
      return (1.24986*((1+(-0.00039734*pt))+(5.37486e-07*(pt*pt))))+(-1.74023e-10*(pt*(pt*(pt/(1+(-0.00112954*pt))))));
    }
    else if (eta < 2.4) {
      return (1.10763*((1+(-0.000105805*pt))+(7.11718e-07*(pt*pt))))+(-5.3001e-10*(pt*(pt*(pt/(1+(-0.000821215*pt))))));
    }
    else {
      WARNING("b-tag efficiency requested for jet |eta| > 2.4");
      return 1;
    }
  }
  else {
    return (1.06268*((1+(0.00390509*pt))+(-5.85405e-05*(pt*pt))))+(7.87135e-07*(pt*(pt*(pt/(1+(0.01259*pt))))));
  }
}
    
double BTagging::GetBMisTagEfficiencyTCHEM(double pt, double eta, bool finebin)
{
  eta = TMath::Abs(eta);
  if (finebin) {
    if (eta < 0.8) {
      return (0.000919586+(0.00026266*pt))+(-1.75723e-07*(pt*pt));
    }
    else if (eta < 1.6) {
      return (-0.00364137+(0.000350371*pt))+(-1.89967e-07*(pt*pt));
    }
    else if (eta <= 2.4) {
      return (-0.00483904+(0.000367751*pt))+(-1.36152e-07*(pt*pt));
    }
    else {
      WARNING("b-tag efficiency requested for jet |eta| > 2.4");
      return 1;
    }
  }
  else {
    // covering full range in eta
    return (-0.00256163+(0.000332759*pt))+(-2.39887e-07*(pt*pt));
  }
}

double BTagging::GetBTagEfficiencyDataTCHEM(double btag)
{
  return -3.67153247396e-07*btag*btag*btag*btag +  -2.81599797034e-05*btag*btag*btag
    + 0.00293190163243*btag*btag +  -0.0849600849778*btag +  0.928524440715;
}

double BTagging::GetBTagEfficiencyDataErrorTCHEM(double btag)
{
  return 3.03337430722e-06*btag*btag*btag*btag + -0.000171604835897*btag*btag*btag
    + 0.00474711667943*btag*btag + -0.0929933040514*btag + 0.978347619293
    - GetBTagEfficiencyDataErrorTCHEM(btag);
}

double BTagging::GetBTagEfficiencyMCTCHEM(double btag)
{
  return 3.90732786802e-06*btag*btag*btag*btag +  -0.000239934437355*btag*btag*btag
    +  0.00664986827287*btag*btag +  -0.112578996016*btag +  1.00775721404;
}

double BTagging::GetCTagEfficiencyMCTCHEM(double btag)
{
  return 0.343760640168*exp(-0.00315525164823*btag*btag*btag +
			    0.0805427315196*btag*btag + -0.867625139194*btag
			    + 1.44815935164 );
}

void BTagging::modifyBTagsWithSF(bool & isBTagged, int pdgIdPart,
				 double Btag_SF, double Btag_eff, double Ctag_eff,
				 double Bmistag_SF, double Bmistag_eff)
{
  // assume no modification
  pdgIdPart = TMath::Abs(pdgIdPart);

  if(pdgIdPart == 5) {
    // b quarks
    isBTagged = applySF(isBTagged, Btag_SF, Btag_eff);
  }
  else if (pdgIdPart == 4)  {
    // c quarks
    isBTagged = applySF(isBTagged, Btag_SF, Ctag_eff);
  }
  else if((pdgIdPart >= 1 && pdgIdPart <= 3) || pdgIdPart == 21) {
    // light quarks, gluons
    isBTagged = applySF(isBTagged, Bmistag_SF, Bmistag_eff);
  }
  else if (pdgIdPart == 0) {
    // in data it is 0, do nothing
  }
  else {
    // treat everything else as unaffected by b-tagging
    // INFO("unidentified particle " << pdgIdPart << " given to BTagging::modifyBTagsWithSF()");
  }
}

bool BTagging::applySF(bool isBTagged, double Btag_SF, double Btag_eff)
{
  bool newBTag = isBTagged;

  if (Btag_SF == 1)
    return newBTag; // no correction needed

  // throw dice
  double coin = rand_->Uniform(1.);

  if(Btag_SF > 1){  // use this if SF>1
    if( !isBTagged ) {
      // fraction of jets that need to be upgraded
      double mistagPercent = (1.0 - Btag_SF) / (1.0 - (Btag_SF/Btag_eff) );

      // upgrade to tagged
      if( coin < mistagPercent ) {
	newBTag = true;
      }
    }
  }
  else {  // use this if SF<1
    // downgrade tagged to untagged
    if( isBTagged && coin > Btag_SF ) {
      newBTag = false;
    }
  }

  return newBTag;
}
