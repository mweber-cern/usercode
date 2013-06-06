#ifndef TreeContent_h
#define TreeContent_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <string>

using namespace std;

class TreeContent {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        global_weight;
   Int_t           global_procID;
   Double_t        global_qscale;
   Double_t        global_bfield;
   Int_t           global_store;
   UInt_t          global_run;
   UInt_t          global_event;
   Int_t           global_bx;
   Int_t           global_orbit;
   Int_t           global_exp;
   Int_t           global_isdata;
   Double_t        global_rhoEG;
   Double_t        global_rho;
   Int_t           lumi_section;
   Double_t        lumi_del;
   Double_t        lumi_rec;
   Double_t        lumi_delerr;
   Double_t        lumi_recerr;
   Int_t           pu_n;
   Int_t           pu_vtxn;
   Double_t        pu_TrueNrInter;
   Int_t           pu_num_int[3];   //[pu_n]
   Double_t        pu_inst_Lumi[3][100];   //[pu_n]
   Double_t        pu_zPos[3][100];   //[pu_n]
   Double_t        pu_sumPthi[3][100];   //[pu_n]
   Double_t        pu_sumPtlo[3][100];   //[pu_n]
   Int_t           noise_pLoose;
   Int_t           noise_pTight;
   Int_t           noise_pHigh;
   Double_t        noise_ecal_r9;
   Double_t        noise_ecal_E;
   Double_t        noise_ecal_pt;
   Double_t        noise_ecal_px;
   Double_t        noise_ecal_py;
   Double_t        noise_ecal_pz;
   Double_t        noise_ecal_eta;
   Double_t        noise_ecal_phi;
   Double_t        noise_ecal_time;
   Double_t        noise_ecal_chi;
   Int_t           noise_ecal_flag;
   Int_t           noise_ecal_ieta;
   Int_t           noise_ecal_iphi;
   Int_t           eventfilter_n;
   Int_t          eventfilter_results[100];   //[eventfilter_n]
   vector<string>  *eventfilter_names;
   vector<string>  *trig_name;
   Int_t           trig_n;
   Int_t           trig_L1prescale[7000];   //[trig_n]
   Int_t           trig_HLTprescale[7000];   //[trig_n]
   vector<string>  *trig_filter;
   Int_t           trigFilter_n;
   Int_t           trig_filterid[7000];   //[trigFilter_n]
   Int_t           trig_id[7000];   //[trigFilter_n]
   Double_t        trig_pt[7000];   //[trigFilter_n]
   Double_t        trig_eta[7000];   //[trigFilter_n]
   Double_t        trig_phi[7000];   //[trigFilter_n]
   Int_t           truth_n;
   Int_t           truth_pdgid[100];   //[truth_n]
   Int_t           truth_bvtxid[100];   //[truth_n]
   Int_t           truth_evtxid[100];   //[truth_n]
   Double_t        truth_E[100];   //[truth_n]
   Double_t        truth_Et[100];   //[truth_n]
   Double_t        truth_p[100];   //[truth_n]
   Double_t        truth_pt[100];   //[truth_n]
   Double_t        truth_px[100];   //[truth_n]
   Double_t        truth_py[100];   //[truth_n]
   Double_t        truth_pz[100];   //[truth_n]
   Double_t        truth_eta[100];   //[truth_n]
   Double_t        truth_phi[100];   //[truth_n]
   Double_t        truth_m[100];   //[truth_n]
   Int_t           truth_status[100];   //[truth_n]
   Int_t           truthl_n;
   Int_t           truthl_ori[100];   //[truthl_n]
   Int_t           truthl_pdgid[100];   //[truthl_n]
   Double_t        truthl_E[100];   //[truthl_n]
   Double_t        truthl_Et[100];   //[truthl_n]
   Double_t        truthl_p[100];   //[truthl_n]
   Double_t        truthl_pt[100];   //[truthl_n]
   Double_t        truthl_px[100];   //[truthl_n]
   Double_t        truthl_py[100];   //[truthl_n]
   Double_t        truthl_pz[100];   //[truthl_n]
   Double_t        truthl_eta[100];   //[truthl_n]
   Double_t        truthl_phi[100];   //[truthl_n]
   Int_t           pdf_id1;
   Int_t           pdf_id2;
   Float_t         pdf_x1;
   Float_t         pdf_x2;
   Float_t         pdf_f1;
   Float_t         pdf_f2;
   Float_t         pdf_scale;
   Int_t           vtx_n;
   Int_t           vtx_ntr[100];   //[vtx_n]
   Int_t           vtx_fake[100];   //[vtx_n]
   Double_t        vtx_ndof[100];   //[vtx_n]
   Double_t        vtx_x[100];   //[vtx_n]
   Double_t        vtx_y[100];   //[vtx_n]
   Double_t        vtx_z[100];   //[vtx_n]
   Double_t        vtx_chi[100];   //[vtx_n]
   Double_t        bs_x;
   Double_t        bs_y;
   Double_t        bs_z;
   Int_t           tracks_n;
   Double_t        tracks_hqf;
   Double_t        met_et[3];
   Double_t        met_ex[3];
   Double_t        met_ey[3];
   Double_t        met_phi[3];
   Double_t        met_sumet[3];
   Double_t        met_sumetsig[3];
   Double_t        met_etsignif[3];
   Double_t        met_ChargedEMEtFraction[3];
   Double_t        met_ChargedHadEtFraction[3];
   Double_t        met_MuonEtFraction[3];
   Double_t        met_NeutralEMFraction[3];
   Double_t        met_NeutralHadEtFraction[3];
   Double_t        met_Type6EtFraction[3];
   Double_t        met_Type7EtFraction[3];
   Double_t        Genmet_et[2];
   Double_t        Genmet_ex[2];
   Double_t        Genmet_ey[2];
   Double_t        Genmet_phi[2];
   Double_t        Genmet_sumet[2];
   Double_t        Genmet_sumetsig[2];
   Double_t        Systmet_et[12];
   Double_t        Systmet_phi[12];
   Int_t           SystJet_ResUp_n;
   Double_t        SystJet_ResUp_pt[100];   //[SystJet_ResUp_n]
   Double_t        SystJet_ResUp_eta[100];   //[SystJet_ResUp_n]
   Double_t        SystJet_ResUp_phi[100];   //[SystJet_ResUp_n]
   Int_t           SystJet_ResDown_n;
   Double_t        SystJet_ResDown_pt[100];   //[SystJet_ResDown_n]
   Double_t        SystJet_ResDown_eta[100];   //[SystJet_ResDown_n]
   Double_t        SystJet_ResDown_phi[100];   //[SystJet_ResDown_n]
   Int_t           pfjet_n;
   Double_t        pfjet_E[100];   //[pfjet_n]
   Double_t        pfjet_Et[100];   //[pfjet_n]
   Double_t        pfjet_p[100];   //[pfjet_n]
   Double_t        pfjet_pt[100];   //[pfjet_n]
   Double_t        pfjet_pt_raw[100];   //[pfjet_n]
   Double_t        pfjet_px[100];   //[pfjet_n]
   Double_t        pfjet_py[100];   //[pfjet_n]
   Double_t        pfjet_pz[100];   //[pfjet_n]
   Double_t        pfjet_eta[100];   //[pfjet_n]
   Double_t        pfjet_phi[100];   //[pfjet_n]
   Double_t        pfjet_btag[100][18];   //[pfjet_n]
   Double_t        pfjet_charge[100];   //[pfjet_n]
   Int_t           pfjet_n90[100];   //[pfjet_n]
   Int_t           pfjet_flav[100];   //[pfjet_n]
   Int_t           pfjet_truth[100];   //[pfjet_n]
   Int_t           pfjet_const[100];   //[pfjet_n]
   Int_t           pfjet_PFN[100][7];   //[pfjet_n]
   Double_t        pfjet_PFF[100][7];   //[pfjet_n]
   Int_t           truthjet_n;
   Double_t        truthjet_E[100];   //[truthjet_n]
   Double_t        truthjet_Et[100];   //[truthjet_n]
   Double_t        truthjet_p[100];   //[truthjet_n]
   Double_t        truthjet_pt[100];   //[truthjet_n]
   Double_t        truthjet_px[100];   //[truthjet_n]
   Double_t        truthjet_py[100];   //[truthjet_n]
   Double_t        truthjet_pz[100];   //[truthjet_n]
   Double_t        truthjet_eta[100];   //[truthjet_n]
   Double_t        truthjet_phi[100];   //[truthjet_n]
   Int_t           SC_n;
   Int_t           SC_truth[200];   //[SC_n]
   Double_t        SC_E[200];   //[SC_n]
   Double_t        SC_phi[200];   //[SC_n]
   Double_t        SC_eta[200];   //[SC_n]
   Int_t           pho_n;
   Int_t           pho_truth[100];   //[pho_n]
   Double_t        pho_PFCand_px[100];   //[pho_n]
   Double_t        pho_PFCand_py[100];   //[pho_n]
   Double_t        pho_PFCand_pz[100];   //[pho_n]
   Double_t        pho_PFCand_E[100];   //[pho_n]
   Double_t        pho_PFCand_eta[100];   //[pho_n]
   Double_t        pho_PFCand_phi[100];   //[pho_n]
   Int_t           pho_PFCand_pfid[100];   //[pho_n]
   Double_t        pho_PFCand_DeltaR[100];   //[pho_n]
   Double_t        pho_E[100];   //[pho_n]
   Double_t        pho_RawE[100];   //[pho_n]
   Double_t        pho_Et[100];   //[pho_n]
   Double_t        pho_p[100];   //[pho_n]
   Double_t        pho_pt[100];   //[pho_n]
   Double_t        pho_px[100];   //[pho_n]
   Double_t        pho_py[100];   //[pho_n]
   Double_t        pho_pz[100];   //[pho_n]
   Double_t        pho_eta[100];   //[pho_n]
   Double_t        pho_phi[100];   //[pho_n]
   Double_t        pho_SCeta[100];   //[pho_n]
   Double_t        pho_Dr03TkSumPt[100];   //[pho_n]
   Double_t        pho_Dr04TkSumPt[100];   //[pho_n]
   Double_t        pho_SigmaIetaIeta[100];   //[pho_n]
   Double_t        pho_SwissCross[100];   //[pho_n]
   Double_t        pho_e5x5[100];   //[pho_n]
   Double_t        pho_e1x5[100];   //[pho_n]
   Double_t        pho_e3x3[100];   //[pho_n]
   Double_t        pho_HCalOverEm[100];   //[pho_n]
   Double_t        pho_HTowOverEm[100];   //[pho_n]
   Int_t           pho_isPF[100];   //[pho_n]
   Double_t        pho_EcaloIso[100];   //[pho_n]
   Double_t        pho_HcaloIso[100];   //[pho_n]
   Double_t        pho_TrackIso[100];   //[pho_n]
   Double_t        pho_PFisoEG[100][3];   //[pho_n]
   Double_t        pho_Roundness[100];   //[pho_n]
   Double_t        pho_Angle[100];   //[pho_n]
   Double_t        pho_R9[100];   //[pho_n]
   Double_t        pho_SCEtaWidth[100];   //[pho_n]
   Int_t           pho_HasPixelSeed[100];   //[pho_n]
   Int_t           pho_HasConvTracks[100];   //[pho_n]
   Int_t          pho_HasMatchedPromptElectron[100];   //[pho_n]
   Int_t           ele_n;
   Double_t        ele_E[100];   //[ele_n]
   Double_t        ele_Et[100];   //[ele_n]
   Double_t        ele_p[100];   //[ele_n]
   Double_t        ele_pt[100];   //[ele_n]
   Double_t        ele_TrackptError[100];   //[ele_n]
   Double_t        ele_Trackpt[100];   //[ele_n]
   Double_t        ele_px[100];   //[ele_n]
   Double_t        ele_py[100];   //[ele_n]
   Double_t        ele_pz[100];   //[ele_n]
   Double_t        ele_eta[100];   //[ele_n]
   Double_t        ele_phi[100];   //[ele_n]
   Double_t        ele_charge[100];   //[ele_n]
   Double_t        ele_TrkChiNorm[100];   //[ele_n]
   Double_t        ele_d0vtx[100];   //[ele_n]
   Double_t        ele_dzvtx[100];   //[ele_n]
   Double_t        ele_d0bs[100];   //[ele_n]
   Double_t        ele_sd0[100];   //[ele_n]
   Int_t           ele_hits[100];   //[ele_n]
   Int_t           ele_truth[100];   //[ele_n]
   Int_t           ele_isECal[100];   //[ele_n]
   Int_t           ele_isTracker[100];   //[ele_n]
   Int_t           ele_ValidHitFirstPxlB[100];   //[ele_n]
   Int_t           ele_TrkExpHitsInner[100];   //[ele_n]
   Double_t        ele_HCalOverEm[100];   //[ele_n]
   Double_t        ele_HCalOverEmBc[100];   //[ele_n]
   Double_t        ele_Dr03TkSumPt[100];   //[ele_n]
   Double_t        ele_Dr04HCalSumEt[100];   //[ele_n]
   Double_t        ele_Dr03HCalSumEt[100];   //[ele_n]
   Double_t        ele_Dr04HCalSumEtBc[100];   //[ele_n]
   Double_t        ele_Dr03HCalSumEtBc[100];   //[ele_n]
   Double_t        ele_Dr04ECalSumEt[100];   //[ele_n]
   Double_t        ele_Dr03ECalSumEt[100];   //[ele_n]
   Double_t        ele_SigmaIetaIeta[100];   //[ele_n]
   Double_t        ele_dEtaSCTrackAtVtx[100];   //[ele_n]
   Double_t        ele_dPhiSCTrackAtVtx[100];   //[ele_n]
   Double_t        ele_dr03HcalDepth1[100];   //[ele_n]
   Double_t        ele_dr03HcalDepth2[100];   //[ele_n]
   Double_t        ele_e2x5Max[100];   //[ele_n]
   Double_t        ele_e5x5[100];   //[ele_n]
   Double_t        ele_e1x5[100];   //[ele_n]
   Double_t        ele_caloEt[100];   //[ele_n]
   Double_t        ele_SCeta[100];   //[ele_n]
   Double_t        ele_convdist[100];   //[ele_n]
   Double_t        ele_convdcot[100];   //[ele_n]
   Double_t        ele_convr[100];   //[ele_n]
   Double_t        ele_fbrem[100];   //[ele_n]
   Int_t           ele_SC[100];   //[ele_n]
   Double_t        ele_PFiso[100][3];   //[ele_n]
   Double_t        ele_PFCand_px[100];   //[ele_n]
   Double_t        ele_PFCand_py[100];   //[ele_n]
   Double_t        ele_PFCand_pz[100];   //[ele_n]
   Double_t        ele_PFCand_E[100];   //[ele_n]
   Double_t        ele_PFCand_eta[100];   //[ele_n]
   Double_t        ele_PFCand_phi[100];   //[ele_n]
   Int_t           ele_PFCand_pfid[100];   //[ele_n]
   Double_t        ele_PFCand_DeltaR[100];   //[ele_n]
   Double_t        ele_hcalDepth1TowerSumEt03[100];   //[ele_n]
   Double_t        ele_hcalDepth2TowerSumEt03[100];   //[ele_n]
   Double_t        ele_EcalEnergy[100];   //[ele_n]
   Double_t        ele_TrackMomentumAtVtx[100];   //[ele_n]
   Double_t        ele_SwissCross[100];   //[ele_n]
   Double_t        ele_EoverP[100];   //[ele_n]
   Int_t           ele_Classification[100];   //[ele_n]
   Int_t           ele_HasMatchedConversions[100];   //[ele_n]
   Double_t        ele_SCRawEt[100];   //[ele_n]
   Double_t        ele_SCEt[100];   //[ele_n] 
   Int_t           pfele_n;
   Double_t        pfele_p[100];   //[pfele_n]
   Double_t        pfele_E[100];   //[pfele_n]
   Double_t        pfele_Et[100];   //[pfele_n]
   Double_t        pfele_CaloEt[100];   //[pfele_n]
   Double_t        pfele_pt[100];   //[pfele_n]
   Double_t        pfele_TrackptError[100];   //[pfele_n]
   Double_t        pfele_Trackpt[100];   //[pfele_n]
   Double_t        pfele_px[100];   //[pfele_n]
   Double_t        pfele_py[100];   //[pfele_n]
   Double_t        pfele_pz[100];   //[pfele_n]
   Double_t        pfele_eta[100];   //[pfele_n]
   Double_t        pfele_phi[100];   //[pfele_n]
   Double_t        pfele_charge[100];   //[pfele_n]
   Int_t           pfele_truth[100];   //[pfele_n]
   Int_t           pfele_SC[100];   //[pfele_n]
   Double_t        pfele_SwissCross[100];   //[pfele_n]
   Double_t        pfele_caloEt[100];   //[pfele_n]
   Double_t        pfele_SCeta[100];   //[pfele_n]
   Double_t        pfele_HCalOverEm[100];   //[pfele_n]
   Double_t        pfele_Dr03TkSumPt[100];   //[pfele_n]
   Double_t        pfele_Dr04HCalSumEt[100];   //[pfele_n]
   Double_t        pfele_Dr03HCalSumEt[100];   //[pfele_n]
   Double_t        pfele_Dr04ECalSumEt[100];   //[pfele_n]
   Double_t        pfele_Dr03ECalSumEt[100];   //[pfele_n]
   Double_t        pfele_particleIso[100];   //[pfele_n]
   Double_t        pfele_chadIso[100];   //[pfele_n]
   Double_t        pfele_nhadIso[100];   //[pfele_n]
   Double_t        pfele_gamIso[100];   //[pfele_n]
   Int_t           muo_n;
   Double_t        muo_E[100];   //[muo_n]
   Double_t        muo_Et[100];   //[muo_n]
   Double_t        muo_p[100];   //[muo_n]
   Double_t        muo_pt[100];   //[muo_n]
   Double_t        muo_px[100];   //[muo_n]
   Double_t        muo_py[100];   //[muo_n]
   Double_t        muo_pz[100];   //[muo_n]
   Double_t        muo_eta[100];   //[muo_n]
   Double_t        muo_phi[100];   //[muo_n]
   Double_t        muo_charge[100];   //[muo_n]
   Double_t        muo_RelTrkIso[100];   //[muo_n]
   Double_t        muo_TrkIso[100];   //[muo_n]
   Double_t        muo_ECalIso[100];   //[muo_n]
   Double_t        muo_HCalIso[100];   //[muo_n]
   Double_t        muo_TrkIsoDep[100];   //[muo_n]
   Double_t        muo_ECalIsoDep[100];   //[muo_n]
   Double_t        muo_HCalIsoDep[100];   //[muo_n]
   Double_t        muo_AllIso[100];   //[muo_n]
   Double_t        muo_TrkChiNormTk[100];   //[muo_n]
   Double_t        muo_d0Tk[100];   //[muo_n]
   Double_t        muo_sd0Tk[100];   //[muo_n]
   Double_t        muo_dzTk[100];   //[muo_n]
   Double_t        muo_calocomp[100];   //[muo_n]
   Double_t        muo_calotower_e[100];   //[muo_n]
   Int_t           muo_prompttight[100];   //[muo_n]
   Int_t           muo_hitsTk[100];   //[muo_n]
   Int_t           muo_truth[100];   //[muo_n]
   Int_t           muo_ID[100][24];   //[muo_n]
   Int_t           muo_ChambersMatched[100];   //[muo_n]
   Int_t           muo_StationsMatched[100];   //[muo_n]
   Double_t        muo_Valid_fraction[100];   //[muo_n]
   Double_t        muo_TrkChiNormCm[100];   //[muo_n]
   Int_t           muo_hitsCm[100];   //[muo_n]
   Double_t        muo_d0Cm[100];   //[muo_n]
   Double_t        muo_sd0Cm[100];   //[muo_n]
   Double_t        muo_d0OriginCm[100];   //[muo_n]
   Double_t        muo_vx[100];   //[muo_n]
   Double_t        muo_vy[100];   //[muo_n]
   Double_t        muo_vz[100];   //[muo_n]
   Int_t           muo_ValidMuonHitsCm[100];   //[muo_n]
   Int_t           muo_ValidTrackerHitsCm[100];   //[muo_n]
   Int_t           muo_ValidPixelHitsCm[100];   //[muo_n]
   Int_t           muo_TrackerLayersMeasCm[100];   //[muo_n]
   Int_t           muo_TrackerLayersNotMeasCm[100];   //[muo_n]
   Int_t           muo_ValidMuonHitsTk[100];   //[muo_n]
   Int_t           muo_ValidTrackerHitsTk[100];   //[muo_n]
   Int_t           muo_ValidPixelHitsTk[100];   //[muo_n]
   Int_t           muo_TrackerLayersMeasTk[100];   //[muo_n]
   Int_t           muo_TrackerLayersNotMeasTk[100];   //[muo_n]
   Int_t           muo_LostHits[100];   //[muo_n]
   Int_t           muo_LostHitsTk[100];   //[muo_n]
   Int_t           muo_isPFMuon[100];   //[muo_n]
   Int_t           muo_DiMuonVertexValid[100][100];
   Int_t           muo_DiMuonVertexNdf[100][100];
   Double_t        muo_DiMuonVertexChi2[100][100];
   Double_t        muo_DiMuonVertexMass[100][100];
   Double_t        muo_Cocktail_pt[100];   //[muo_n]
   Double_t        muo_Cocktail_phi[100];   //[muo_n]
   Double_t        muo_Cocktail_eta[100];   //[muo_n]
   Double_t        muo_Tracker_pt[100]; //[muo_n]
   Double_t        muo_Tracker_ptError[100]; //[muo_n]
   Double_t        muo_Tracker_phi[100]; //[muo_n]
   Double_t        muo_Tracker_eta[100]; //[muo_n]
   Double_t        muo_TevReco_pt[100][7];   //[muo_n]
   Double_t        muo_TevReco_ptError[100][7];   //[muo_n]
   Double_t        muo_TevReco_eta[100][7];   //[muo_n]
   Double_t        muo_TevReco_phi[100][7];   //[muo_n]
   Double_t        muo_TevReco_chi2[100][7];   //[muo_n]
   Double_t        muo_TevReco_ndof[100][7];   //[muo_n]
   Double_t        muo_PFiso[100][7];   //[muo_n]
   Double_t        muo_PFCand_px[100];   //[muo_n]
   Double_t        muo_PFCand_py[100];   //[muo_n]
   Double_t        muo_PFCand_pz[100];   //[muo_n]
   Double_t        muo_PFCand_E[100];   //[muo_n]
   Double_t        muo_PFCand_eta[100];   //[muo_n]
   Double_t        muo_PFCand_phi[100];   //[muo_n]
   Int_t           muo_PFCand_pfid[100];   //[muo_n]
   Double_t        muo_PFCand_DeltaR[100];   //[muo_n]
   Int_t           tau_n;
   Double_t        tau_p[100];   //[tau_n]
   Double_t        tau_pt[100];   //[tau_n]
   Double_t        tau_E[100];   //[tau_n]
   Double_t        tau_Et[100];   //[tau_n]
   Double_t        tau_M[100];   //[tau_n]
   Double_t        tau_Mt[100];   //[tau_n]
   Double_t        tau_Px[100];   //[tau_n]
   Double_t        tau_Py[100];   //[tau_n]
   Double_t        tau_Pz[100];   //[tau_n]
   Double_t        tau_Eta[100];   //[tau_n]
   Double_t        tau_Phi[100];   //[tau_n]
   Double_t        tau_AltPt[100];   //[tau_n]
   Int_t           tau_DecayMode[100];   //[tau_n]
   Int_t           tau_Charge[100];   //[tau_n]
   Double_t        tau_ParticleIso[100];   //[tau_n]
   Int_t           tau_PFChargedHadrCands[100];   //[tau_n]
   Int_t           tau_PFGammaCands[100];   //[tau_n]
   Double_t        tau_IsolationPFChargedHadrCandsPtSum[100];   //[tau_n]
   Double_t        tau_EcalStripSumEOverPLead[100];   //[tau_n]
   Double_t        tau_LeadPFChargedHadrCandsignedSipt[100];   //[tau_n]
   Int_t           tau_NSignalTracks[100];   //[tau_n]
   Double_t        tau_PFLeadChargedPT[100];   //[tau_n]
   Double_t        tau_Jet_pt[100];   //[tau_n]
   Double_t        tau_Jet_eta[100];   //[tau_n]
   Double_t        tau_Jet_phi[100];   //[tau_n]
   Double_t        tau_Jet_m[100];   //[tau_n]
   Double_t        tau_vtx_x[100];   //[tau_n]
   Double_t        tau_vtx_y[100];   //[tau_n]
   Double_t        tau_vtx_z[100];   //[tau_n]
   Double_t        tau_id[100][31];   //[tau_n]
   Int_t           tau_GenJet_Match_n;
   Int_t           tau_GenJet_DecayMode[100];   //[tau_GenJet_Match_n]
   Int_t           tau_GenJetMatch_Pos[100];   //[tau_GenJet_Match_n]
   Double_t        tau_GenJet_E[100];   //[tau_GenJet_Match_n]
   Double_t        tau_GenJet_Et[100];   //[tau_GenJet_Match_n]
   Double_t        tau_GenJet_Eta[100];   //[tau_GenJet_Match_n]
   Double_t        tau_GenJet_Phi[100];   //[tau_GenJet_Match_n]
   Double_t        tau_GenJet_M[100];   //[tau_GenJet_Match_n]
   Double_t        tau_GenJet_Mt[100];   //[tau_GenJet_Match_n]
   Double_t        tau_GenJet_P[100];   //[tau_GenJet_Match_n]
   Double_t        tau_GenJet_Pt[100];   //[tau_GenJet_Match_n]
   Double_t        tau_GenJet_Px[100];   //[tau_GenJet_Match_n]
   Double_t        tau_GenJet_Py[100];   //[tau_GenJet_Match_n]
   Double_t        tau_GenJet_Pz[100];   //[tau_GenJet_Match_n]
   Double_t        susyScanM0;
   Double_t        susyScanM12;
   Double_t        susyScanA0;
   Double_t        susyScanCrossSection;
   Double_t        susyScanMu;
   Double_t        susyScanRun;
   Double_t        susyScantanbeta;

   // List of branches
   TBranch        *b_global_weight;   //!
   TBranch        *b_global_procID;   //!
   TBranch        *b_global_qscale;   //!
   TBranch        *b_global_bfield;   //!
   TBranch        *b_global_store;   //!
   TBranch        *b_global_run;   //!
   TBranch        *b_global_event;   //!
   TBranch        *b_global_bx;   //!
   TBranch        *b_global_orbit;   //!
   TBranch        *b_global_exp;   //!
   TBranch        *b_global_isdata;   //!
   TBranch        *b_global_rhoEG;   //!
   TBranch        *b_global_rho;   //!
   TBranch        *b_lumi_section;   //!
   TBranch        *b_lumi_del;   //!
   TBranch        *b_lumi_rec;   //!
   TBranch        *b_lumi_delerr;   //!
   TBranch        *b_lumi_recerr;   //!
   TBranch        *b_pu_n;   //!
   TBranch        *b_pu_vtxn;   //!
   TBranch        *b_pu_TrueNrInter;   //!
   TBranch        *b_pu_num_int;   //!
   TBranch        *b_pu_inst_Lumi;   //!
   TBranch        *b_pu_zPos;   //!
   TBranch        *b_pu_sumPthi;   //!
   TBranch        *b_pu_sumPtlo;   //!
   TBranch        *b_noise_pLoose;   //!
   TBranch        *b_noise_pTight;   //!
   TBranch        *b_noise_pHigh;   //!
   TBranch        *b_noise_ecal_r9;   //!
   TBranch        *b_noise_ecal_E;   //!
   TBranch        *b_noise_ecal_pt;   //!
   TBranch        *b_noise_ecal_px;   //!
   TBranch        *b_noise_ecal_py;   //!
   TBranch        *b_noise_ecal_pz;   //!
   TBranch        *b_noise_ecal_eta;   //!
   TBranch        *b_noise_ecal_phi;   //!
   TBranch        *b_noise_ecal_time;   //!
   TBranch        *b_noise_ecal_chi;   //!
   TBranch        *b_noise_ecal_flag;   //!
   TBranch        *b_noise_ecal_ieta;   //!
   TBranch        *b_noise_ecal_iphi;   //!
   TBranch        *b_eventfilter_n;   //!
   TBranch        *b_eventfilter_results;   //!
   TBranch        *b_eventfilter_names;   //!
   TBranch        *b_trig_name;   //!
   TBranch        *b_trig_n;   //!
   TBranch        *b_trig_L1prescale;   //!
   TBranch        *b_trig_HLTprescale;   //!
   TBranch        *b_trig_filter;   //!
   TBranch        *b_trigFilter_n;   //!
   TBranch        *b_trig_filterid;   //!
   TBranch        *b_trig_id;   //!
   TBranch        *b_trig_pt;   //!
   TBranch        *b_trig_eta;   //!
   TBranch        *b_trig_phi;   //!
   TBranch        *b_truth_n;   //!
   TBranch        *b_truth_pdgid;   //!
   TBranch        *b_truth_bvtxid;   //!
   TBranch        *b_truth_evtxid;   //!
   TBranch        *b_truth_E;   //!
   TBranch        *b_truth_Et;   //!
   TBranch        *b_truth_p;   //!
   TBranch        *b_truth_pt;   //!
   TBranch        *b_truth_px;   //!
   TBranch        *b_truth_py;   //!
   TBranch        *b_truth_pz;   //!
   TBranch        *b_truth_eta;   //!
   TBranch        *b_truth_phi;   //!
   TBranch        *b_truth_m;   //!
   TBranch        *b_truth_status;   //!
   TBranch        *b_truthl_n;   //!
   TBranch        *b_truthl_ori;   //!
   TBranch        *b_truthl_pdgid;   //!
   TBranch        *b_truthl_E;   //!
   TBranch        *b_truthl_Et;   //!
   TBranch        *b_truthl_p;   //!
   TBranch        *b_truthl_pt;   //!
   TBranch        *b_truthl_px;   //!
   TBranch        *b_truthl_py;   //!
   TBranch        *b_truthl_pz;   //!
   TBranch        *b_truthl_eta;   //!
   TBranch        *b_truthl_phi;   //!
   TBranch        *b_pdf_id1;   //!
   TBranch        *b_pdf_id2;   //!
   TBranch        *b_pdf_x1;   //!
   TBranch        *b_pdf_x2;   //!
   TBranch        *b_pdf_f1;   //!
   TBranch        *b_pdf_f2;   //!
   TBranch        *b_pdf_scale;   //!
   TBranch        *b_vtx_n;   //!
   TBranch        *b_vtx_ntr;   //!
   TBranch        *b_vtx_fake;   //!
   TBranch        *b_vtx_ndof;   //!
   TBranch        *b_vtx_x;   //!
   TBranch        *b_vtx_y;   //!
   TBranch        *b_vtx_z;   //!
   TBranch        *b_vtx_chi;   //!
   TBranch        *b_bs_x;   //!
   TBranch        *b_bs_y;   //!
   TBranch        *b_bs_z;   //!
   TBranch        *b_tracks_n;   //!
   TBranch        *b_tracks_hqf;   //!
   TBranch        *b_met_et;   //!
   TBranch        *b_met_ex;   //!
   TBranch        *b_met_ey;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_sumet;   //!
   TBranch        *b_met_sumetsig;   //!
   TBranch        *b_met_etsignif;   //!
   TBranch        *b_met_ChargedEMEtFraction;   //!
   TBranch        *b_met_ChargedHadEtFraction;   //!
   TBranch        *b_met_MuonEtFraction;   //!
   TBranch        *b_met_NeutralEMFraction;   //!
   TBranch        *b_met_NeutralHadEtFraction;   //!
   TBranch        *b_met_Type6EtFraction;   //!
   TBranch        *b_met_Type7EtFraction;   //!
   TBranch        *b_Genmet_et;   //!
   TBranch        *b_Genmet_ex;   //!
   TBranch        *b_Genmet_ey;   //!
   TBranch        *b_Genmet_phi;   //!
   TBranch        *b_Genmet_sumet;   //!
   TBranch        *b_Genmet_sumetsig;   //!
   TBranch        *b_Systmet_et;   //!
   TBranch        *b_Systmet_phi;   //!
   TBranch        *b_SystJet_ResUp_n;   //!
   TBranch        *b_SystJet_ResUp_pt;   //!
   TBranch        *b_SystJet_ResUp_eta;   //!
   TBranch        *b_SystJet_ResUp_phi;   //!
   TBranch        *b_SystJet_ResDown_n;   //!
   TBranch        *b_SystJet_ResDown_pt;   //!
   TBranch        *b_SystJet_ResDown_eta;   //!
   TBranch        *b_SystJet_ResDown_phi;   //!
   TBranch        *b_pfjet_n;   //!
   TBranch        *b_pfjet_E;   //!
   TBranch        *b_pfjet_Et;   //!
   TBranch        *b_pfjet_p;   //!
   TBranch        *b_pfjet_pt;   //!
   TBranch        *b_pfjet_pt_raw;   //!
   TBranch        *b_pfjet_px;   //!
   TBranch        *b_pfjet_py;   //!
   TBranch        *b_pfjet_pz;   //!
   TBranch        *b_pfjet_eta;   //!
   TBranch        *b_pfjet_phi;   //!
   TBranch        *b_pfjet_btag;   //!
   TBranch        *b_pfjet_charge;   //!
   TBranch        *b_pfjet_n90;   //!
   TBranch        *b_pfjet_flav;   //!
   TBranch        *b_pfjet_truth;   //!
   TBranch        *b_pfjet_const;   //!
   TBranch        *b_pfjet_PFN;   //!
   TBranch        *b_pfjet_PFF;   //!
   TBranch        *b_truthjet_n;   //!
   TBranch        *b_truthjet_E;   //!
   TBranch        *b_truthjet_Et;   //!
   TBranch        *b_truthjet_p;   //!
   TBranch        *b_truthjet_pt;   //!
   TBranch        *b_truthjet_px;   //!
   TBranch        *b_truthjet_py;   //!
   TBranch        *b_truthjet_pz;   //!
   TBranch        *b_truthjet_eta;   //!
   TBranch        *b_truthjet_phi;   //!
   TBranch        *b_SC_n;   //!
   TBranch        *b_SC_truth;   //!
   TBranch        *b_SC_E;   //!
   TBranch        *b_SC_phi;   //!
   TBranch        *b_SC_eta;   //!
   TBranch        *b_pho_n;   //!
   TBranch        *b_pho_truth;   //!
   TBranch        *b_pho_PFCand_px;   //!
   TBranch        *b_pho_PFCand_py;   //!
   TBranch        *b_pho_PFCand_pz;   //!
   TBranch        *b_pho_PFCand_E;   //!
   TBranch        *b_pho_PFCand_eta;   //!
   TBranch        *b_pho_PFCand_phi;   //!
   TBranch        *b_pho_PFCand_pfid;   //!
   TBranch        *b_pho_PFCand_DeltaR;   //!
   TBranch        *b_pho_E;   //!
   TBranch        *b_pho_RawE;   //!
   TBranch        *b_pho_Et;   //!
   TBranch        *b_pho_p;   //!
   TBranch        *b_pho_pt;   //!
   TBranch        *b_pho_px;   //!
   TBranch        *b_pho_py;   //!
   TBranch        *b_pho_pz;   //!
   TBranch        *b_pho_eta;   //!
   TBranch        *b_pho_phi;   //!
   TBranch        *b_pho_SCeta;   //!
   TBranch        *b_pho_Dr03TkSumPt;   //!
   TBranch        *b_pho_Dr04TkSumPt;   //!
   TBranch        *b_pho_SigmaIetaIeta;   //!
   TBranch        *b_pho_SwissCross;   //!
   TBranch        *b_pho_e5x5;   //!
   TBranch        *b_pho_e1x5;   //!
   TBranch        *b_pho_e3x3;   //!
   TBranch        *b_pho_HCalOverEm;   //!
   TBranch        *b_pho_HTowOverEm;   //!
   TBranch        *b_pho_isPF;   //!
   TBranch        *b_pho_EcaloIso;   //!
   TBranch        *b_pho_HcaloIso;   //!
   TBranch        *b_pho_TrackIso;   //!
   TBranch        *b_pho_PFisoEG;   //!
   TBranch        *b_pho_Roundness;   //!
   TBranch        *b_pho_Angle;   //!
   TBranch        *b_pho_R9;   //!
   TBranch        *b_pho_SCEtaWidth;   //!
   TBranch        *b_pho_HasPixelSeed;   //!
   TBranch        *b_pho_HasConvTracks;   //!
   TBranch        *b_pho_HasMatchedPromptElectron;   //!
   TBranch        *b_ele_n;   //!
   TBranch        *b_ele_E;   //!
   TBranch        *b_ele_Et;   //!
   TBranch        *b_ele_p;   //!
   TBranch        *b_ele_pt;   //!
   TBranch        *b_ele_TrackptError;   //!
   TBranch        *b_ele_Trackpt;   //!
   TBranch        *b_ele_px;   //!
   TBranch        *b_ele_py;   //!
   TBranch        *b_ele_pz;   //!
   TBranch        *b_ele_eta;   //!
   TBranch        *b_ele_phi;   //!
   TBranch        *b_ele_charge;   //!
   TBranch        *b_ele_TrkChiNorm;   //!
   TBranch        *b_ele_d0vtx;   //!
   TBranch        *b_ele_dzvtx;   //!
   TBranch        *b_ele_d0bs;   //!
   TBranch        *b_ele_sd0;   //!
   TBranch        *b_ele_hits;   //!
   TBranch        *b_ele_truth;   //!
   TBranch        *b_ele_isECal;   //!
   TBranch        *b_ele_isTracker;   //!
   TBranch        *b_ele_ValidHitFirstPxlB;   //!
   TBranch        *b_ele_TrkExpHitsInner;   //!
   TBranch        *b_ele_HCalOverEm;   //!
   TBranch        *b_ele_HCalOverEmBc;   //!
   TBranch        *b_ele_Dr03TkSumPt;   //!
   TBranch        *b_ele_Dr04HCalSumEt;   //!
   TBranch        *b_ele_Dr03HCalSumEt;   //!
   TBranch        *b_ele_Dr04HCalSumEtBc;   //!
   TBranch        *b_ele_Dr03HCalSumEtBc;   //!
   TBranch        *b_ele_Dr04ECalSumEt;   //!
   TBranch        *b_ele_Dr03ECalSumEt;   //!
   TBranch        *b_ele_SigmaIetaIeta;   //!
   TBranch        *b_ele_dEtaSCTrackAtVtx;   //!
   TBranch        *b_ele_dPhiSCTrackAtVtx;   //!
   TBranch        *b_ele_dr03HcalDepth1;   //!
   TBranch        *b_ele_dr03HcalDepth2;   //!
   TBranch        *b_ele_e2x5Max;   //!
   TBranch        *b_ele_e5x5;   //!
   TBranch        *b_ele_e1x5;   //!
   TBranch        *b_ele_caloEt;   //!
   TBranch        *b_ele_SCeta;   //!
   TBranch        *b_ele_convdist;   //!
   TBranch        *b_ele_convdcot;   //!
   TBranch        *b_ele_convr;   //!
   TBranch        *b_ele_fbrem;   //!
   TBranch        *b_ele_SC;   //!
   TBranch        *b_ele_PFiso;   //!
   TBranch        *b_ele_PFCand_px;   //!
   TBranch        *b_ele_PFCand_py;   //!
   TBranch        *b_ele_PFCand_pz;   //!
   TBranch        *b_ele_PFCand_E;   //!
   TBranch        *b_ele_PFCand_eta;   //!
   TBranch        *b_ele_PFCand_phi;   //!
   TBranch        *b_ele_PFCand_pfid;   //!
   TBranch        *b_ele_PFCand_DeltaR;   //!
   TBranch        *b_ele_hcalDepth1TowerSumEt03;   //!
   TBranch        *b_ele_hcalDepth2TowerSumEt03;   //!
   TBranch        *b_ele_EcalEnergy;   //!
   TBranch        *b_ele_TrackMomentumAtVtx;   //!
   TBranch        *b_ele_SwissCross;   //!
   TBranch        *b_ele_EoverP;   //!
   TBranch        *b_ele_Classification;   //!
   TBranch        *b_ele_HasMatchedConversions;   //!
   TBranch        *b_ele_SCRawEt;   //!
   TBranch        *b_ele_SCEt;   //!
   TBranch        *b_pfele_n;   //!
   TBranch        *b_pfele_p;   //!
   TBranch        *b_pfele_E;   //!
   TBranch        *b_pfele_Et;   //!
   TBranch        *b_pfele_CaloEt;   //!
   TBranch        *b_pfele_pt;   //!
   TBranch        *b_pfele_TrackptError;   //!
   TBranch        *b_pfele_Trackpt;   //!
   TBranch        *b_pfele_px;   //!
   TBranch        *b_pfele_py;   //!
   TBranch        *b_pfele_pz;   //!
   TBranch        *b_pfele_eta;   //!
   TBranch        *b_pfele_phi;   //!
   TBranch        *b_pfele_charge;   //!
   TBranch        *b_pfele_truth;   //!
   TBranch        *b_pfele_SC;   //!
   TBranch        *b_pfele_SwissCross;   //!
   TBranch        *b_pfele_caloEt;   //!
   TBranch        *b_pfele_SCeta;   //!
   TBranch        *b_pfele_HCalOverEm;   //!
   TBranch        *b_pfele_Dr03TkSumPt;   //!
   TBranch        *b_pfele_Dr04HCalSumEt;   //!
   TBranch        *b_pfele_Dr03HCalSumEt;   //!
   TBranch        *b_pfele_Dr04ECalSumEt;   //!
   TBranch        *b_pfele_Dr03ECalSumEt;   //!
   TBranch        *b_pfele_particleIso;   //!
   TBranch        *b_pfele_chadIso;   //!
   TBranch        *b_pfele_nhadIso;   //!
   TBranch        *b_pfele_gamIso;   //!
   TBranch        *b_muo_n;   //!
   TBranch        *b_muo_E;   //!
   TBranch        *b_muo_Et;   //!
   TBranch        *b_muo_p;   //!
   TBranch        *b_muo_pt;   //!
   TBranch        *b_muo_px;   //!
   TBranch        *b_muo_py;   //!
   TBranch        *b_muo_pz;   //!
   TBranch        *b_muo_eta;   //!
   TBranch        *b_muo_phi;   //!
   TBranch        *b_muo_charge;   //!
   TBranch        *b_muo_RelTrkIso;   //!
   TBranch        *b_muo_TrkIso;   //!
   TBranch        *b_muo_ECalIso;   //!
   TBranch        *b_muo_HCalIso;   //!
   TBranch        *b_muo_TrkIsoDep;   //!
   TBranch        *b_muo_ECalIsoDep;   //!
   TBranch        *b_muo_HCalIsoDep;   //!
   TBranch        *b_muo_AllIso;   //!
   TBranch        *b_muo_TrkChiNormTk;   //!
   TBranch        *b_muo_d0Tk;   //!
   TBranch        *b_muo_sd0Tk;   //!
   TBranch        *b_muo_dzTk;   //!
   TBranch        *b_muo_calocomp;   //!
   TBranch        *b_muo_calotower_e;   //!
   TBranch        *b_muo_prompttight;   //!
   TBranch        *b_muo_hitsTk;   //!
   TBranch        *b_muo_truth;   //!
   TBranch        *b_muo_ID;   //!
   TBranch        *b_muo_ChambersMatched;   //!
   TBranch        *b_muo_StationsMatched;   //!
   TBranch        *b_muo_Valid_fraction;   //!
   TBranch        *b_muo_TrkChiNormCm;   //!
   TBranch        *b_muo_hitsCm;   //!
   TBranch        *b_muo_d0Cm;   //!
   TBranch        *b_muo_sd0Cm;   //!
   TBranch        *b_muo_d0OriginCm;   //!
   TBranch        *b_muo_vx;   //!
   TBranch        *b_muo_vy;   //!
   TBranch        *b_muo_vz;   //!
   TBranch        *b_muo_ValidMuonHitsCm;   //!
   TBranch        *b_muo_ValidTrackerHitsCm;   //!
   TBranch        *b_muo_ValidPixelHitsCm;   //!
   TBranch        *b_muo_TrackerLayersMeasCm;   //!
   TBranch        *b_muo_TrackerLayersNotMeasCm;   //!
   TBranch        *b_muo_ValidMuonHitsTk;   //!
   TBranch        *b_muo_ValidTrackerHitsTk;   //!
   TBranch        *b_muo_ValidPixelHitsTk;   //!
   TBranch        *b_muo_TrackerLayersMeasTk;   //!
   TBranch        *b_muo_TrackerLayersNotMeasTk;   //!
   TBranch        *b_muo_LostHits;   //!
   TBranch        *b_muo_LostHitsTk;   //!
   TBranch        *b_muo_isPFMuon;   //!
   TBranch        *b_muo_DiMuonVertexValid;   //!
   TBranch        *b_muo_DiMuonVertexNdf;   //!
   TBranch        *b_muo_DiMuonVertexChi2;   //!
   TBranch        *b_muo_DiMuonVertexMass;   //!     
   TBranch        *b_muo_Cocktail_pt;   //!
   TBranch        *b_muo_Cocktail_phi;   //!
   TBranch        *b_muo_Cocktail_eta;   //!
   TBranch        *b_muo_Tracker_pt; //!
   TBranch        *b_muo_Tracker_ptError; //!
   TBranch        *b_muo_Tracker_phi; //!
   TBranch        *b_muo_Tracker_eta; //!
   TBranch        *b_muo_TevReco_pt;   //!
   TBranch        *b_muo_TevReco_ptError;   //!
   TBranch        *b_muo_TevReco_eta;   //!
   TBranch        *b_muo_TevReco_phi;   //!
   TBranch        *b_muo_TevReco_chi2;   //!
   TBranch        *b_muo_TevReco_ndof;   //!
   TBranch        *b_muo_PFiso;   //!
   TBranch        *b_muo_PFCand_px;   //!
   TBranch        *b_muo_PFCand_py;   //!
   TBranch        *b_muo_PFCand_pz;   //!
   TBranch        *b_muo_PFCand_E;   //!
   TBranch        *b_muo_PFCand_eta;   //!
   TBranch        *b_muo_PFCand_phi;   //!
   TBranch        *b_muo_PFCand_pfid;   //!
   TBranch        *b_muo_PFCand_DeltaR;   //!
   TBranch        *b_tau_n;   //!
   TBranch        *b_tau_p;   //!
   TBranch        *b_tau_pt;   //!
   TBranch        *b_tau_E;   //!
   TBranch        *b_tau_Et;   //!
   TBranch        *b_tau_M;   //!
   TBranch        *b_tau_Mt;   //!
   TBranch        *b_tau_Px;   //!
   TBranch        *b_tau_Py;   //!
   TBranch        *b_tau_Pz;   //!
   TBranch        *b_tau_Eta;   //!
   TBranch        *b_tau_Phi;   //!
   TBranch        *b_tau_AltPt;   //!
   TBranch        *b_tau_DecayMode;   //!
   TBranch        *b_tau_Charge;   //!
   TBranch        *b_tau_ParticleIso;   //!
   TBranch        *b_tau_PFChargedHadrCands;   //!
   TBranch        *b_tau_PFGammaCands;   //!
   TBranch        *b_tau_IsolationPFChargedHadrCandsPtSum;   //!
   TBranch        *b_tau_EcalStripSumEOverPLead;   //!
   TBranch        *b_tau_LeadPFChargedHadrCandsignedSipt;   //!
   TBranch        *b_tau_NSignalTracks;   //!
   TBranch        *b_tau_PFLeadChargedPT;   //!
   TBranch        *b_tau_Jet_pt;   //!
   TBranch        *b_tau_Jet_eta;   //!
   TBranch        *b_tau_Jet_phi;   //!
   TBranch        *b_tau_Jet_m;   //!
   TBranch        *b_tau_vtx_x;   //!
   TBranch        *b_tau_vtx_y;   //!
   TBranch        *b_tau_vtx_z;   //!
   TBranch        *b_tau_id;   //!
   TBranch        *b_tau_GenJet_Match_n;   //!
   TBranch        *b_tau_GenJet_DecayMode;   //!
   TBranch        *b_tau_GenJetMatch_Pos;   //!
   TBranch        *b_tau_GenJet_E;   //!
   TBranch        *b_tau_GenJet_Et;   //!
   TBranch        *b_tau_GenJet_Eta;   //!
   TBranch        *b_tau_GenJet_Phi;   //!
   TBranch        *b_tau_GenJet_M;   //!
   TBranch        *b_tau_GenJet_Mt;   //!
   TBranch        *b_tau_GenJet_P;   //!
   TBranch        *b_tau_GenJet_Pt;   //!
   TBranch        *b_tau_GenJet_Px;   //!
   TBranch        *b_tau_GenJet_Py;   //!
   TBranch        *b_tau_GenJet_Pz;   //!
   TBranch        *b_susyScanM0;   //!
   TBranch        *b_susyScanM12;   //!
   TBranch        *b_susyScanA0;   //!
   TBranch        *b_susyScanCrossSection;   //!
   TBranch        *b_susyScanMu;   //!
   TBranch        *b_susyScanRun;   //!
   TBranch        *b_susyScantanbeta;   //!

   TreeContent(TTree *tree=0);
   virtual ~TreeContent();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TreeContent_cxx
TreeContent::TreeContent(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("out.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("out.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("out.root:/ACSkimAnalysis");
      dir->GetObject("allData",tree);

   }
   Init(tree);
}

TreeContent::~TreeContent()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TreeContent::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TreeContent::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TreeContent::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   eventfilter_names = 0;
   trig_name = 0;
   trig_filter = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("global_weight", &global_weight, &b_global_weight);
   fChain->SetBranchAddress("global_procID", &global_procID, &b_global_procID);
   fChain->SetBranchAddress("global_qscale", &global_qscale, &b_global_qscale);
   fChain->SetBranchAddress("global_bfield", &global_bfield, &b_global_bfield);
   fChain->SetBranchAddress("global_store", &global_store, &b_global_store);
   fChain->SetBranchAddress("global_run", &global_run, &b_global_run);
   fChain->SetBranchAddress("global_event", &global_event, &b_global_event);
   fChain->SetBranchAddress("global_bx", &global_bx, &b_global_bx);
   fChain->SetBranchAddress("global_orbit", &global_orbit, &b_global_orbit);
   fChain->SetBranchAddress("global_exp", &global_exp, &b_global_exp);
   fChain->SetBranchAddress("global_isdata", &global_isdata, &b_global_isdata);
   fChain->SetBranchAddress("global_rhoEG", &global_rhoEG, &b_global_rhoEG);
   fChain->SetBranchAddress("global_rho", &global_rho, &b_global_rho);
   fChain->SetBranchAddress("lumi_section", &lumi_section, &b_lumi_section);
   fChain->SetBranchAddress("lumi_del", &lumi_del, &b_lumi_del);
   fChain->SetBranchAddress("lumi_rec", &lumi_rec, &b_lumi_rec);
   fChain->SetBranchAddress("lumi_delerr", &lumi_delerr, &b_lumi_delerr);
   fChain->SetBranchAddress("lumi_recerr", &lumi_recerr, &b_lumi_recerr);
   fChain->SetBranchAddress("pu_n", &pu_n, &b_pu_n);
   fChain->SetBranchAddress("pu_vtxn", &pu_vtxn, &b_pu_vtxn);
   fChain->SetBranchAddress("pu_TrueNrInter", &pu_TrueNrInter, &b_pu_TrueNrInter);
   fChain->SetBranchAddress("pu_num_int", &pu_num_int, &b_pu_num_int);
   fChain->SetBranchAddress("pu_inst_Lumi", &pu_inst_Lumi, &b_pu_inst_Lumi);
   fChain->SetBranchAddress("pu_zPos", &pu_zPos, &b_pu_zPos);
   fChain->SetBranchAddress("pu_sumPthi", &pu_sumPthi, &b_pu_sumPthi);
   fChain->SetBranchAddress("pu_sumPtlo", &pu_sumPtlo, &b_pu_sumPtlo);
   fChain->SetBranchAddress("noise_pLoose", &noise_pLoose, &b_noise_pLoose);
   fChain->SetBranchAddress("noise_pTight", &noise_pTight, &b_noise_pTight);
   fChain->SetBranchAddress("noise_pHigh", &noise_pHigh, &b_noise_pHigh);
   fChain->SetBranchAddress("noise_ecal_r9", &noise_ecal_r9, &b_noise_ecal_r9);
   fChain->SetBranchAddress("noise_ecal_E", &noise_ecal_E, &b_noise_ecal_E);
   fChain->SetBranchAddress("noise_ecal_pt", &noise_ecal_pt, &b_noise_ecal_pt);
   fChain->SetBranchAddress("noise_ecal_px", &noise_ecal_px, &b_noise_ecal_px);
   fChain->SetBranchAddress("noise_ecal_py", &noise_ecal_py, &b_noise_ecal_py);
   fChain->SetBranchAddress("noise_ecal_pz", &noise_ecal_pz, &b_noise_ecal_pz);
   fChain->SetBranchAddress("noise_ecal_eta", &noise_ecal_eta, &b_noise_ecal_eta);
   fChain->SetBranchAddress("noise_ecal_phi", &noise_ecal_phi, &b_noise_ecal_phi);
   fChain->SetBranchAddress("noise_ecal_time", &noise_ecal_time, &b_noise_ecal_time);
   fChain->SetBranchAddress("noise_ecal_chi", &noise_ecal_chi, &b_noise_ecal_chi);
   fChain->SetBranchAddress("noise_ecal_flag", &noise_ecal_flag, &b_noise_ecal_flag);
   fChain->SetBranchAddress("noise_ecal_ieta", &noise_ecal_ieta, &b_noise_ecal_ieta);
   fChain->SetBranchAddress("noise_ecal_iphi", &noise_ecal_iphi, &b_noise_ecal_iphi);
   fChain->SetBranchAddress("eventfilter_n", &eventfilter_n, &b_eventfilter_n);
   fChain->SetBranchAddress("eventfilter_results", eventfilter_results, &b_eventfilter_results);
   fChain->SetBranchAddress("eventfilter_names", &eventfilter_names, &b_eventfilter_names);
   fChain->SetBranchAddress("trig_name", &trig_name, &b_trig_name);
   fChain->SetBranchAddress("trig_n", &trig_n, &b_trig_n);
   fChain->SetBranchAddress("trig_L1prescale", trig_L1prescale, &b_trig_L1prescale);
   fChain->SetBranchAddress("trig_HLTprescale", trig_HLTprescale, &b_trig_HLTprescale);
   fChain->SetBranchAddress("trig_filter", &trig_filter, &b_trig_filter);
   fChain->SetBranchAddress("trigFilter_n", &trigFilter_n, &b_trigFilter_n);
   fChain->SetBranchAddress("trig_filterid", trig_filterid, &b_trig_filterid);
   fChain->SetBranchAddress("trig_id", trig_id, &b_trig_id);
   fChain->SetBranchAddress("trig_pt", trig_pt, &b_trig_pt);
   fChain->SetBranchAddress("trig_eta", trig_eta, &b_trig_eta);
   fChain->SetBranchAddress("trig_phi", trig_phi, &b_trig_phi);
   fChain->SetBranchAddress("truth_n", &truth_n, &b_truth_n);
   fChain->SetBranchAddress("truth_pdgid", &truth_pdgid, &b_truth_pdgid);
   fChain->SetBranchAddress("truth_bvtxid", &truth_bvtxid, &b_truth_bvtxid);
   fChain->SetBranchAddress("truth_evtxid", &truth_evtxid, &b_truth_evtxid);
   fChain->SetBranchAddress("truth_E", &truth_E, &b_truth_E);
   fChain->SetBranchAddress("truth_Et", &truth_Et, &b_truth_Et);
   fChain->SetBranchAddress("truth_p", &truth_p, &b_truth_p);
   fChain->SetBranchAddress("truth_pt", &truth_pt, &b_truth_pt);
   fChain->SetBranchAddress("truth_px", &truth_px, &b_truth_px);
   fChain->SetBranchAddress("truth_py", &truth_py, &b_truth_py);
   fChain->SetBranchAddress("truth_pz", &truth_pz, &b_truth_pz);
   fChain->SetBranchAddress("truth_eta", &truth_eta, &b_truth_eta);
   fChain->SetBranchAddress("truth_phi", &truth_phi, &b_truth_phi);
   fChain->SetBranchAddress("truth_m", &truth_m, &b_truth_m);
   fChain->SetBranchAddress("truth_status", &truth_status, &b_truth_status);
   fChain->SetBranchAddress("truthl_n", &truthl_n, &b_truthl_n);
   fChain->SetBranchAddress("truthl_ori", &truthl_ori, &b_truthl_ori);
   fChain->SetBranchAddress("truthl_pdgid", &truthl_pdgid, &b_truthl_pdgid);
   fChain->SetBranchAddress("truthl_E", &truthl_E, &b_truthl_E);
   fChain->SetBranchAddress("truthl_Et", &truthl_Et, &b_truthl_Et);
   fChain->SetBranchAddress("truthl_p", &truthl_p, &b_truthl_p);
   fChain->SetBranchAddress("truthl_pt", &truthl_pt, &b_truthl_pt);
   fChain->SetBranchAddress("truthl_px", &truthl_px, &b_truthl_px);
   fChain->SetBranchAddress("truthl_py", &truthl_py, &b_truthl_py);
   fChain->SetBranchAddress("truthl_pz", &truthl_pz, &b_truthl_pz);
   fChain->SetBranchAddress("truthl_eta", &truthl_eta, &b_truthl_eta);
   fChain->SetBranchAddress("truthl_phi", &truthl_phi, &b_truthl_phi);
   fChain->SetBranchAddress("pdf_id1", &pdf_id1, &b_pdf_id1);
   fChain->SetBranchAddress("pdf_id2", &pdf_id2, &b_pdf_id2);
   fChain->SetBranchAddress("pdf_x1", &pdf_x1, &b_pdf_x1);
   fChain->SetBranchAddress("pdf_x2", &pdf_x2, &b_pdf_x2);
   fChain->SetBranchAddress("pdf_f1", &pdf_f1, &b_pdf_f1);
   fChain->SetBranchAddress("pdf_f2", &pdf_f2, &b_pdf_f2);
   fChain->SetBranchAddress("pdf_scale", &pdf_scale, &b_pdf_scale);
   fChain->SetBranchAddress("vtx_n", &vtx_n, &b_vtx_n);
   fChain->SetBranchAddress("vtx_ntr", vtx_ntr, &b_vtx_ntr);
   fChain->SetBranchAddress("vtx_fake", vtx_fake, &b_vtx_fake);
   fChain->SetBranchAddress("vtx_ndof", vtx_ndof, &b_vtx_ndof);
   fChain->SetBranchAddress("vtx_x", vtx_x, &b_vtx_x);
   fChain->SetBranchAddress("vtx_y", vtx_y, &b_vtx_y);
   fChain->SetBranchAddress("vtx_z", vtx_z, &b_vtx_z);
   fChain->SetBranchAddress("vtx_chi", vtx_chi, &b_vtx_chi);
   fChain->SetBranchAddress("bs_x", &bs_x, &b_bs_x);
   fChain->SetBranchAddress("bs_y", &bs_y, &b_bs_y);
   fChain->SetBranchAddress("bs_z", &bs_z, &b_bs_z);
   fChain->SetBranchAddress("tracks_n", &tracks_n, &b_tracks_n);
   fChain->SetBranchAddress("tracks_hqf", &tracks_hqf, &b_tracks_hqf);
   fChain->SetBranchAddress("met_et", met_et, &b_met_et);
   fChain->SetBranchAddress("met_ex", met_ex, &b_met_ex);
   fChain->SetBranchAddress("met_ey", met_ey, &b_met_ey);
   fChain->SetBranchAddress("met_phi", met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_sumet", met_sumet, &b_met_sumet);
   fChain->SetBranchAddress("met_sumetsig", met_sumetsig, &b_met_sumetsig);
   fChain->SetBranchAddress("met_etsignif", met_etsignif, &b_met_etsignif);
   fChain->SetBranchAddress("met_ChargedEMEtFraction", met_ChargedEMEtFraction, &b_met_ChargedEMEtFraction);
   fChain->SetBranchAddress("met_ChargedHadEtFraction", met_ChargedHadEtFraction, &b_met_ChargedHadEtFraction);
   fChain->SetBranchAddress("met_MuonEtFraction", met_MuonEtFraction, &b_met_MuonEtFraction);
   fChain->SetBranchAddress("met_NeutralEMFraction", met_NeutralEMFraction, &b_met_NeutralEMFraction);
   fChain->SetBranchAddress("met_NeutralHadEtFraction", met_NeutralHadEtFraction, &b_met_NeutralHadEtFraction);
   fChain->SetBranchAddress("met_Type6EtFraction", met_Type6EtFraction, &b_met_Type6EtFraction);
   fChain->SetBranchAddress("met_Type7EtFraction", met_Type7EtFraction, &b_met_Type7EtFraction);
   fChain->SetBranchAddress("Genmet_et", Genmet_et, &b_Genmet_et);
   fChain->SetBranchAddress("Genmet_ex", Genmet_ex, &b_Genmet_ex);
   fChain->SetBranchAddress("Genmet_ey", Genmet_ey, &b_Genmet_ey);
   fChain->SetBranchAddress("Genmet_phi", Genmet_phi, &b_Genmet_phi);
   fChain->SetBranchAddress("Genmet_sumet", Genmet_sumet, &b_Genmet_sumet);
   fChain->SetBranchAddress("Genmet_sumetsig", Genmet_sumetsig, &b_Genmet_sumetsig);
   fChain->SetBranchAddress("Systmet_et", Systmet_et, &b_Systmet_et);
   fChain->SetBranchAddress("Systmet_phi", Systmet_phi, &b_Systmet_phi);
   fChain->SetBranchAddress("SystJet_ResUp_n", &SystJet_ResUp_n, &b_SystJet_ResUp_n);
   fChain->SetBranchAddress("SystJet_ResUp_pt", &SystJet_ResUp_pt, &b_SystJet_ResUp_pt);
   fChain->SetBranchAddress("SystJet_ResUp_eta", &SystJet_ResUp_eta, &b_SystJet_ResUp_eta);
   fChain->SetBranchAddress("SystJet_ResUp_phi", &SystJet_ResUp_phi, &b_SystJet_ResUp_phi);
   fChain->SetBranchAddress("SystJet_ResDown_n", &SystJet_ResDown_n, &b_SystJet_ResDown_n);
   fChain->SetBranchAddress("SystJet_ResDown_pt", &SystJet_ResDown_pt, &b_SystJet_ResDown_pt);
   fChain->SetBranchAddress("SystJet_ResDown_eta", &SystJet_ResDown_eta, &b_SystJet_ResDown_eta);
   fChain->SetBranchAddress("SystJet_ResDown_phi", &SystJet_ResDown_phi, &b_SystJet_ResDown_phi);
   fChain->SetBranchAddress("pfjet_n", &pfjet_n, &b_pfjet_n);
   fChain->SetBranchAddress("pfjet_E", pfjet_E, &b_pfjet_E);
   fChain->SetBranchAddress("pfjet_Et", pfjet_Et, &b_pfjet_Et);
   fChain->SetBranchAddress("pfjet_p", pfjet_p, &b_pfjet_p);
   fChain->SetBranchAddress("pfjet_pt", pfjet_pt, &b_pfjet_pt);
   fChain->SetBranchAddress("pfjet_pt_raw", pfjet_pt_raw, &b_pfjet_pt_raw);
   fChain->SetBranchAddress("pfjet_px", pfjet_px, &b_pfjet_px);
   fChain->SetBranchAddress("pfjet_py", pfjet_py, &b_pfjet_py);
   fChain->SetBranchAddress("pfjet_pz", pfjet_pz, &b_pfjet_pz);
   fChain->SetBranchAddress("pfjet_eta", pfjet_eta, &b_pfjet_eta);
   fChain->SetBranchAddress("pfjet_phi", pfjet_phi, &b_pfjet_phi);
   fChain->SetBranchAddress("pfjet_btag", pfjet_btag, &b_pfjet_btag);
   fChain->SetBranchAddress("pfjet_charge", pfjet_charge, &b_pfjet_charge);
   fChain->SetBranchAddress("pfjet_n90", pfjet_n90, &b_pfjet_n90);
   fChain->SetBranchAddress("pfjet_flav", pfjet_flav, &b_pfjet_flav);
   fChain->SetBranchAddress("pfjet_truth", pfjet_truth, &b_pfjet_truth);
   fChain->SetBranchAddress("pfjet_const", pfjet_const, &b_pfjet_const);
   fChain->SetBranchAddress("pfjet_PFN", pfjet_PFN, &b_pfjet_PFN);
   fChain->SetBranchAddress("pfjet_PFF", pfjet_PFF, &b_pfjet_PFF);
   fChain->SetBranchAddress("truthjet_n", &truthjet_n, &b_truthjet_n);
   fChain->SetBranchAddress("truthjet_E", &truthjet_E, &b_truthjet_E);
   fChain->SetBranchAddress("truthjet_Et", &truthjet_Et, &b_truthjet_Et);
   fChain->SetBranchAddress("truthjet_p", &truthjet_p, &b_truthjet_p);
   fChain->SetBranchAddress("truthjet_pt", &truthjet_pt, &b_truthjet_pt);
   fChain->SetBranchAddress("truthjet_px", &truthjet_px, &b_truthjet_px);
   fChain->SetBranchAddress("truthjet_py", &truthjet_py, &b_truthjet_py);
   fChain->SetBranchAddress("truthjet_pz", &truthjet_pz, &b_truthjet_pz);
   fChain->SetBranchAddress("truthjet_eta", &truthjet_eta, &b_truthjet_eta);
   fChain->SetBranchAddress("truthjet_phi", &truthjet_phi, &b_truthjet_phi);
   fChain->SetBranchAddress("SC_n", &SC_n, &b_SC_n);
   fChain->SetBranchAddress("SC_truth", SC_truth, &b_SC_truth);
   fChain->SetBranchAddress("SC_E", SC_E, &b_SC_E);
   fChain->SetBranchAddress("SC_phi", SC_phi, &b_SC_phi);
   fChain->SetBranchAddress("SC_eta", SC_eta, &b_SC_eta);
   fChain->SetBranchAddress("pho_n", &pho_n, &b_pho_n);
   fChain->SetBranchAddress("pho_truth", pho_truth, &b_pho_truth);
   fChain->SetBranchAddress("pho_PFCand_px", pho_PFCand_px, &b_pho_PFCand_px);
   fChain->SetBranchAddress("pho_PFCand_py", pho_PFCand_py, &b_pho_PFCand_py);
   fChain->SetBranchAddress("pho_PFCand_pz", pho_PFCand_pz, &b_pho_PFCand_pz);
   fChain->SetBranchAddress("pho_PFCand_E", pho_PFCand_E, &b_pho_PFCand_E);
   fChain->SetBranchAddress("pho_PFCand_eta", pho_PFCand_eta, &b_pho_PFCand_eta);
   fChain->SetBranchAddress("pho_PFCand_phi", pho_PFCand_phi, &b_pho_PFCand_phi);
   fChain->SetBranchAddress("pho_PFCand_pfid", pho_PFCand_pfid, &b_pho_PFCand_pfid);
   fChain->SetBranchAddress("pho_PFCand_DeltaR", pho_PFCand_DeltaR, &b_pho_PFCand_DeltaR);
   fChain->SetBranchAddress("pho_E", pho_E, &b_pho_E);
   fChain->SetBranchAddress("pho_RawE", pho_RawE, &b_pho_RawE);
   fChain->SetBranchAddress("pho_Et", pho_Et, &b_pho_Et);
   fChain->SetBranchAddress("pho_p", pho_p, &b_pho_p);
   fChain->SetBranchAddress("pho_pt", pho_pt, &b_pho_pt);
   fChain->SetBranchAddress("pho_px", pho_px, &b_pho_px);
   fChain->SetBranchAddress("pho_py", pho_py, &b_pho_py);
   fChain->SetBranchAddress("pho_pz", pho_pz, &b_pho_pz);
   fChain->SetBranchAddress("pho_eta", pho_eta, &b_pho_eta);
   fChain->SetBranchAddress("pho_phi", pho_phi, &b_pho_phi);
   fChain->SetBranchAddress("pho_SCeta", pho_SCeta, &b_pho_SCeta);
   fChain->SetBranchAddress("pho_Dr03TkSumPt", pho_Dr03TkSumPt, &b_pho_Dr03TkSumPt);
   fChain->SetBranchAddress("pho_Dr04TkSumPt", pho_Dr04TkSumPt, &b_pho_Dr04TkSumPt);
   fChain->SetBranchAddress("pho_SigmaIetaIeta", pho_SigmaIetaIeta, &b_pho_SigmaIetaIeta);
   fChain->SetBranchAddress("pho_SwissCross", pho_SwissCross, &b_pho_SwissCross);
   fChain->SetBranchAddress("pho_e5x5", pho_e5x5, &b_pho_e5x5);
   fChain->SetBranchAddress("pho_e1x5", pho_e1x5, &b_pho_e1x5);
   fChain->SetBranchAddress("pho_e3x3", pho_e3x3, &b_pho_e3x3);
   fChain->SetBranchAddress("pho_HCalOverEm", pho_HCalOverEm, &b_pho_HCalOverEm);
   fChain->SetBranchAddress("pho_HTowOverEm", pho_HTowOverEm, &b_pho_HTowOverEm);
   fChain->SetBranchAddress("pho_isPF", pho_isPF, &b_pho_isPF);
   fChain->SetBranchAddress("pho_EcaloIso", pho_EcaloIso, &b_pho_EcaloIso);
   fChain->SetBranchAddress("pho_HcaloIso", pho_HcaloIso, &b_pho_HcaloIso);
   fChain->SetBranchAddress("pho_TrackIso", pho_TrackIso, &b_pho_TrackIso);
   fChain->SetBranchAddress("pho_PFisoEG", pho_PFisoEG, &b_pho_PFisoEG);
   fChain->SetBranchAddress("pho_Roundness", pho_Roundness, &b_pho_Roundness);
   fChain->SetBranchAddress("pho_Angle", pho_Angle, &b_pho_Angle);
   fChain->SetBranchAddress("pho_R9", pho_R9, &b_pho_R9);
   fChain->SetBranchAddress("pho_SCEtaWidth", pho_SCEtaWidth, &b_pho_SCEtaWidth);
   fChain->SetBranchAddress("pho_HasPixelSeed", pho_HasPixelSeed, &b_pho_HasPixelSeed);
   fChain->SetBranchAddress("pho_HasConvTracks", pho_HasConvTracks, &b_pho_HasConvTracks);
   fChain->SetBranchAddress("pho_HasMatchedPromptElectron", pho_HasMatchedPromptElectron, &b_pho_HasMatchedPromptElectron);
   fChain->SetBranchAddress("ele_n", &ele_n, &b_ele_n);
   fChain->SetBranchAddress("ele_E", ele_E, &b_ele_E);
   fChain->SetBranchAddress("ele_Et", ele_Et, &b_ele_Et);
   fChain->SetBranchAddress("ele_p", ele_p, &b_ele_p);
   fChain->SetBranchAddress("ele_pt", ele_pt, &b_ele_pt);
   fChain->SetBranchAddress("ele_TrackptError", ele_TrackptError, &b_ele_TrackptError);
   fChain->SetBranchAddress("ele_Trackpt", ele_Trackpt, &b_ele_Trackpt);
   fChain->SetBranchAddress("ele_px", ele_px, &b_ele_px);
   fChain->SetBranchAddress("ele_py", ele_py, &b_ele_py);
   fChain->SetBranchAddress("ele_pz", ele_pz, &b_ele_pz);
   fChain->SetBranchAddress("ele_eta", ele_eta, &b_ele_eta);
   fChain->SetBranchAddress("ele_phi", ele_phi, &b_ele_phi);
   fChain->SetBranchAddress("ele_charge", ele_charge, &b_ele_charge);
   fChain->SetBranchAddress("ele_TrkChiNorm", ele_TrkChiNorm, &b_ele_TrkChiNorm);
   fChain->SetBranchAddress("ele_d0vtx", ele_d0vtx, &b_ele_d0vtx);
   fChain->SetBranchAddress("ele_dzvtx", ele_dzvtx, &b_ele_dzvtx);
   fChain->SetBranchAddress("ele_d0bs", ele_d0bs, &b_ele_d0bs);
   fChain->SetBranchAddress("ele_sd0", ele_sd0, &b_ele_sd0);
   fChain->SetBranchAddress("ele_hits", ele_hits, &b_ele_hits);
   fChain->SetBranchAddress("ele_truth", ele_truth, &b_ele_truth);
   fChain->SetBranchAddress("ele_isECal", ele_isECal, &b_ele_isECal);
   fChain->SetBranchAddress("ele_isTracker", ele_isTracker, &b_ele_isTracker);
   fChain->SetBranchAddress("ele_ValidHitFirstPxlB", ele_ValidHitFirstPxlB, &b_ele_ValidHitFirstPxlB);
   fChain->SetBranchAddress("ele_TrkExpHitsInner", ele_TrkExpHitsInner, &b_ele_TrkExpHitsInner);
   fChain->SetBranchAddress("ele_HCalOverEm", ele_HCalOverEm, &b_ele_HCalOverEm);
   fChain->SetBranchAddress("ele_HCalOverEmBc", ele_HCalOverEmBc, &b_ele_HCalOverEmBc);
   fChain->SetBranchAddress("ele_Dr03TkSumPt", ele_Dr03TkSumPt, &b_ele_Dr03TkSumPt);
   fChain->SetBranchAddress("ele_Dr04HCalSumEt", ele_Dr04HCalSumEt, &b_ele_Dr04HCalSumEt);
   fChain->SetBranchAddress("ele_Dr03HCalSumEt", ele_Dr03HCalSumEt, &b_ele_Dr03HCalSumEt);
   fChain->SetBranchAddress("ele_Dr04HCalSumEtBc", ele_Dr04HCalSumEtBc, &b_ele_Dr04HCalSumEtBc);
   fChain->SetBranchAddress("ele_Dr03HCalSumEtBc", ele_Dr03HCalSumEtBc, &b_ele_Dr03HCalSumEtBc);
   fChain->SetBranchAddress("ele_Dr04ECalSumEt", ele_Dr04ECalSumEt, &b_ele_Dr04ECalSumEt);
   fChain->SetBranchAddress("ele_Dr03ECalSumEt", ele_Dr03ECalSumEt, &b_ele_Dr03ECalSumEt);
   fChain->SetBranchAddress("ele_SigmaIetaIeta", ele_SigmaIetaIeta, &b_ele_SigmaIetaIeta);
   fChain->SetBranchAddress("ele_dEtaSCTrackAtVtx", ele_dEtaSCTrackAtVtx, &b_ele_dEtaSCTrackAtVtx);
   fChain->SetBranchAddress("ele_dPhiSCTrackAtVtx", ele_dPhiSCTrackAtVtx, &b_ele_dPhiSCTrackAtVtx);
   fChain->SetBranchAddress("ele_dr03HcalDepth1", ele_dr03HcalDepth1, &b_ele_dr03HcalDepth1);
   fChain->SetBranchAddress("ele_dr03HcalDepth2", ele_dr03HcalDepth2, &b_ele_dr03HcalDepth2);
   fChain->SetBranchAddress("ele_e2x5Max", ele_e2x5Max, &b_ele_e2x5Max);
   fChain->SetBranchAddress("ele_e5x5", ele_e5x5, &b_ele_e5x5);
   fChain->SetBranchAddress("ele_e1x5", ele_e1x5, &b_ele_e1x5);
   fChain->SetBranchAddress("ele_caloEt", ele_caloEt, &b_ele_caloEt);
   fChain->SetBranchAddress("ele_SCeta", ele_SCeta, &b_ele_SCeta);
   fChain->SetBranchAddress("ele_convdist", ele_convdist, &b_ele_convdist);
   fChain->SetBranchAddress("ele_convdcot", ele_convdcot, &b_ele_convdcot);
   fChain->SetBranchAddress("ele_convr", ele_convr, &b_ele_convr);
   fChain->SetBranchAddress("ele_fbrem", ele_fbrem, &b_ele_fbrem);
   fChain->SetBranchAddress("ele_SC", ele_SC, &b_ele_SC);
   fChain->SetBranchAddress("ele_PFiso", ele_PFiso, &b_ele_PFiso);
   fChain->SetBranchAddress("ele_PFCand_px", ele_PFCand_px, &b_ele_PFCand_px);
   fChain->SetBranchAddress("ele_PFCand_py", ele_PFCand_py, &b_ele_PFCand_py);
   fChain->SetBranchAddress("ele_PFCand_pz", ele_PFCand_pz, &b_ele_PFCand_pz);
   fChain->SetBranchAddress("ele_PFCand_E", ele_PFCand_E, &b_ele_PFCand_E);
   fChain->SetBranchAddress("ele_PFCand_eta", ele_PFCand_eta, &b_ele_PFCand_eta);
   fChain->SetBranchAddress("ele_PFCand_phi", ele_PFCand_phi, &b_ele_PFCand_phi);
   fChain->SetBranchAddress("ele_PFCand_pfid", ele_PFCand_pfid, &b_ele_PFCand_pfid);
   fChain->SetBranchAddress("ele_PFCand_DeltaR", ele_PFCand_DeltaR, &b_ele_PFCand_DeltaR);
   fChain->SetBranchAddress("ele_hcalDepth1TowerSumEt03", ele_hcalDepth1TowerSumEt03, &b_ele_hcalDepth1TowerSumEt03);
   fChain->SetBranchAddress("ele_hcalDepth2TowerSumEt03", ele_hcalDepth2TowerSumEt03, &b_ele_hcalDepth2TowerSumEt03);
   fChain->SetBranchAddress("ele_EcalEnergy", ele_EcalEnergy, &b_ele_EcalEnergy);
   fChain->SetBranchAddress("ele_TrackMomentumAtVtx", ele_TrackMomentumAtVtx, &b_ele_TrackMomentumAtVtx);
   fChain->SetBranchAddress("ele_SwissCross", ele_SwissCross, &b_ele_SwissCross);
   fChain->SetBranchAddress("ele_EoverP", ele_EoverP, &b_ele_EoverP);
   fChain->SetBranchAddress("ele_Classification", ele_Classification, &b_ele_Classification);
   fChain->SetBranchAddress("ele_HasMatchedConversions", ele_HasMatchedConversions, &b_ele_HasMatchedConversions);
   fChain->SetBranchAddress("ele_SCRawEt", ele_SCRawEt, &b_ele_SCRawEt);
   fChain->SetBranchAddress("ele_SCEt", ele_SCEt, &b_ele_SCEt);
   fChain->SetBranchAddress("pfele_n", &pfele_n, &b_pfele_n);
   fChain->SetBranchAddress("pfele_p", &pfele_p, &b_pfele_p);
   fChain->SetBranchAddress("pfele_E", &pfele_E, &b_pfele_E);
   fChain->SetBranchAddress("pfele_Et", &pfele_Et, &b_pfele_Et);
   fChain->SetBranchAddress("pfele_CaloEt", &pfele_CaloEt, &b_pfele_CaloEt);
   fChain->SetBranchAddress("pfele_pt", &pfele_pt, &b_pfele_pt);
   fChain->SetBranchAddress("pfele_TrackptError", &pfele_TrackptError, &b_pfele_TrackptError);
   fChain->SetBranchAddress("pfele_Trackpt", &pfele_Trackpt, &b_pfele_Trackpt);
   fChain->SetBranchAddress("pfele_px", &pfele_px, &b_pfele_px);
   fChain->SetBranchAddress("pfele_py", &pfele_py, &b_pfele_py);
   fChain->SetBranchAddress("pfele_pz", &pfele_pz, &b_pfele_pz);
   fChain->SetBranchAddress("pfele_eta", &pfele_eta, &b_pfele_eta);
   fChain->SetBranchAddress("pfele_phi", &pfele_phi, &b_pfele_phi);
   fChain->SetBranchAddress("pfele_charge", &pfele_charge, &b_pfele_charge);
   fChain->SetBranchAddress("pfele_truth", &pfele_truth, &b_pfele_truth);
   fChain->SetBranchAddress("pfele_SC", &pfele_SC, &b_pfele_SC);
   fChain->SetBranchAddress("pfele_SwissCross", &pfele_SwissCross, &b_pfele_SwissCross);
   fChain->SetBranchAddress("pfele_caloEt", &pfele_caloEt, &b_pfele_caloEt);
   fChain->SetBranchAddress("pfele_SCeta", &pfele_SCeta, &b_pfele_SCeta);
   fChain->SetBranchAddress("pfele_HCalOverEm", &pfele_HCalOverEm, &b_pfele_HCalOverEm);
   fChain->SetBranchAddress("pfele_Dr03TkSumPt", &pfele_Dr03TkSumPt, &b_pfele_Dr03TkSumPt);
   fChain->SetBranchAddress("pfele_Dr04HCalSumEt", &pfele_Dr04HCalSumEt, &b_pfele_Dr04HCalSumEt);
   fChain->SetBranchAddress("pfele_Dr03HCalSumEt", &pfele_Dr03HCalSumEt, &b_pfele_Dr03HCalSumEt);
   fChain->SetBranchAddress("pfele_Dr04ECalSumEt", &pfele_Dr04ECalSumEt, &b_pfele_Dr04ECalSumEt);
   fChain->SetBranchAddress("pfele_Dr03ECalSumEt", &pfele_Dr03ECalSumEt, &b_pfele_Dr03ECalSumEt);
   fChain->SetBranchAddress("pfele_particleIso", &pfele_particleIso, &b_pfele_particleIso);
   fChain->SetBranchAddress("pfele_chadIso", &pfele_chadIso, &b_pfele_chadIso);
   fChain->SetBranchAddress("pfele_nhadIso", &pfele_nhadIso, &b_pfele_nhadIso);
   fChain->SetBranchAddress("pfele_gamIso", &pfele_gamIso, &b_pfele_gamIso);
   fChain->SetBranchAddress("muo_n", &muo_n, &b_muo_n);
   fChain->SetBranchAddress("muo_E", muo_E, &b_muo_E);
   fChain->SetBranchAddress("muo_Et", muo_Et, &b_muo_Et);
   fChain->SetBranchAddress("muo_p", muo_p, &b_muo_p);
   fChain->SetBranchAddress("muo_pt", muo_pt, &b_muo_pt);
   fChain->SetBranchAddress("muo_px", muo_px, &b_muo_px);
   fChain->SetBranchAddress("muo_py", muo_py, &b_muo_py);
   fChain->SetBranchAddress("muo_pz", muo_pz, &b_muo_pz);
   fChain->SetBranchAddress("muo_eta", muo_eta, &b_muo_eta);
   fChain->SetBranchAddress("muo_phi", muo_phi, &b_muo_phi);
   fChain->SetBranchAddress("muo_charge", muo_charge, &b_muo_charge);
   fChain->SetBranchAddress("muo_RelTrkIso", muo_RelTrkIso, &b_muo_RelTrkIso);
   fChain->SetBranchAddress("muo_TrkIso", muo_TrkIso, &b_muo_TrkIso);
   fChain->SetBranchAddress("muo_ECalIso", muo_ECalIso, &b_muo_ECalIso);
   fChain->SetBranchAddress("muo_HCalIso", muo_HCalIso, &b_muo_HCalIso);
   fChain->SetBranchAddress("muo_TrkIsoDep", muo_TrkIsoDep, &b_muo_TrkIsoDep);
   fChain->SetBranchAddress("muo_ECalIsoDep", muo_ECalIsoDep, &b_muo_ECalIsoDep);
   fChain->SetBranchAddress("muo_HCalIsoDep", muo_HCalIsoDep, &b_muo_HCalIsoDep);
   fChain->SetBranchAddress("muo_AllIso", muo_AllIso, &b_muo_AllIso);
   fChain->SetBranchAddress("muo_TrkChiNormTk", muo_TrkChiNormTk, &b_muo_TrkChiNormTk);
   fChain->SetBranchAddress("muo_d0Tk", muo_d0Tk, &b_muo_d0Tk);
   fChain->SetBranchAddress("muo_sd0Tk", muo_sd0Tk, &b_muo_sd0Tk);
   fChain->SetBranchAddress("muo_dzTk", muo_dzTk, &b_muo_dzTk);
   fChain->SetBranchAddress("muo_calocomp", muo_calocomp, &b_muo_calocomp);
   fChain->SetBranchAddress("muo_calotower_e", muo_calotower_e, &b_muo_calotower_e);
   fChain->SetBranchAddress("muo_prompttight", muo_prompttight, &b_muo_prompttight);
   fChain->SetBranchAddress("muo_hitsTk", muo_hitsTk, &b_muo_hitsTk);
   fChain->SetBranchAddress("muo_truth", muo_truth, &b_muo_truth);
   fChain->SetBranchAddress("muo_ID", muo_ID, &b_muo_ID);
   fChain->SetBranchAddress("muo_ChambersMatched", muo_ChambersMatched, &b_muo_ChambersMatched);
   fChain->SetBranchAddress("muo_StationsMatched", muo_StationsMatched, &b_muo_StationsMatched);
   fChain->SetBranchAddress("muo_Valid_fraction", muo_Valid_fraction, &b_muo_Valid_fraction);
   fChain->SetBranchAddress("muo_TrkChiNormCm", muo_TrkChiNormCm, &b_muo_TrkChiNormCm);
   fChain->SetBranchAddress("muo_hitsCm", muo_hitsCm, &b_muo_hitsCm);
   fChain->SetBranchAddress("muo_d0Cm", muo_d0Cm, &b_muo_d0Cm);
   fChain->SetBranchAddress("muo_sd0Cm", muo_sd0Cm, &b_muo_sd0Cm);
   fChain->SetBranchAddress("muo_d0OriginCm", muo_d0OriginCm, &b_muo_d0OriginCm);
   fChain->SetBranchAddress("muo_vx", muo_vx, &b_muo_vx);
   fChain->SetBranchAddress("muo_vy", muo_vy, &b_muo_vy);
   fChain->SetBranchAddress("muo_vz", muo_vz, &b_muo_vz);
   fChain->SetBranchAddress("muo_ValidMuonHitsCm", muo_ValidMuonHitsCm, &b_muo_ValidMuonHitsCm);
   fChain->SetBranchAddress("muo_ValidTrackerHitsCm", muo_ValidTrackerHitsCm, &b_muo_ValidTrackerHitsCm);
   fChain->SetBranchAddress("muo_ValidPixelHitsCm", muo_ValidPixelHitsCm, &b_muo_ValidPixelHitsCm);
   fChain->SetBranchAddress("muo_TrackerLayersMeasCm", muo_TrackerLayersMeasCm, &b_muo_TrackerLayersMeasCm);
   fChain->SetBranchAddress("muo_TrackerLayersNotMeasCm", muo_TrackerLayersNotMeasCm, &b_muo_TrackerLayersNotMeasCm);
   fChain->SetBranchAddress("muo_ValidMuonHitsTk", muo_ValidMuonHitsTk, &b_muo_ValidMuonHitsTk);
   fChain->SetBranchAddress("muo_ValidTrackerHitsTk", muo_ValidTrackerHitsTk, &b_muo_ValidTrackerHitsTk);
   fChain->SetBranchAddress("muo_ValidPixelHitsTk", muo_ValidPixelHitsTk, &b_muo_ValidPixelHitsTk);
   fChain->SetBranchAddress("muo_TrackerLayersMeasTk", muo_TrackerLayersMeasTk, &b_muo_TrackerLayersMeasTk);
   fChain->SetBranchAddress("muo_TrackerLayersNotMeasTk", muo_TrackerLayersNotMeasTk, &b_muo_TrackerLayersNotMeasTk);
   fChain->SetBranchAddress("muo_LostHits", muo_LostHits, &b_muo_LostHits);
   fChain->SetBranchAddress("muo_LostHitsTk", muo_LostHitsTk, &b_muo_LostHitsTk);
   fChain->SetBranchAddress("muo_isPFMuon", muo_isPFMuon, &b_muo_isPFMuon);
   fChain->SetBranchAddress("muo_DiMuonVertexValid", muo_DiMuonVertexValid, &b_muo_DiMuonVertexValid);
   fChain->SetBranchAddress("muo_DiMuonVertexNdf", muo_DiMuonVertexNdf, &b_muo_DiMuonVertexNdf);
   fChain->SetBranchAddress("muo_DiMuonVertexChi2", muo_DiMuonVertexChi2, &b_muo_DiMuonVertexChi2);
   fChain->SetBranchAddress("muo_DiMuonVertexMass", muo_DiMuonVertexMass, &b_muo_DiMuonVertexMass);
   fChain->SetBranchAddress("muo_Cocktail_pt", muo_Cocktail_pt, &b_muo_Cocktail_pt);
   fChain->SetBranchAddress("muo_Cocktail_phi", muo_Cocktail_phi, &b_muo_Cocktail_phi);
   fChain->SetBranchAddress("muo_Cocktail_eta", muo_Cocktail_eta, &b_muo_Cocktail_eta);
   fChain->SetBranchAddress("muo_Tracker_pt", muo_Tracker_pt, &b_muo_Tracker_pt);
   fChain->SetBranchAddress("muo_Tracker_ptError", muo_Tracker_ptError, &b_muo_Tracker_ptError);
   fChain->SetBranchAddress("muo_Tracker_phi", muo_Tracker_phi, &b_muo_Tracker_phi);
   fChain->SetBranchAddress("muo_Tracker_eta", muo_Tracker_eta, &b_muo_Tracker_eta);
   fChain->SetBranchAddress("muo_TevReco_pt", muo_TevReco_pt, &b_muo_TevReco_pt);
   fChain->SetBranchAddress("muo_TevReco_ptError", muo_TevReco_ptError, &b_muo_TevReco_ptError);
   fChain->SetBranchAddress("muo_TevReco_eta", muo_TevReco_eta, &b_muo_TevReco_eta);
   fChain->SetBranchAddress("muo_TevReco_phi", muo_TevReco_phi, &b_muo_TevReco_phi);
   fChain->SetBranchAddress("muo_TevReco_chi2", muo_TevReco_chi2, &b_muo_TevReco_chi2);
   fChain->SetBranchAddress("muo_TevReco_ndof", muo_TevReco_ndof, &b_muo_TevReco_ndof);
   fChain->SetBranchAddress("muo_PFiso", muo_PFiso, &b_muo_PFiso);
   fChain->SetBranchAddress("muo_PFCand_px", muo_PFCand_px, &b_muo_PFCand_px);
   fChain->SetBranchAddress("muo_PFCand_py", muo_PFCand_py, &b_muo_PFCand_py);
   fChain->SetBranchAddress("muo_PFCand_pz", muo_PFCand_pz, &b_muo_PFCand_pz);
   fChain->SetBranchAddress("muo_PFCand_E", muo_PFCand_E, &b_muo_PFCand_E);
   fChain->SetBranchAddress("muo_PFCand_eta", muo_PFCand_eta, &b_muo_PFCand_eta);
   fChain->SetBranchAddress("muo_PFCand_phi", muo_PFCand_phi, &b_muo_PFCand_phi);
   fChain->SetBranchAddress("muo_PFCand_pfid", muo_PFCand_pfid, &b_muo_PFCand_pfid);
   fChain->SetBranchAddress("muo_PFCand_DeltaR", muo_PFCand_DeltaR, &b_muo_PFCand_DeltaR);
   fChain->SetBranchAddress("tau_n", &tau_n, &b_tau_n);
   fChain->SetBranchAddress("tau_p", tau_p, &b_tau_p);
   fChain->SetBranchAddress("tau_pt", tau_pt, &b_tau_pt);
   fChain->SetBranchAddress("tau_E", tau_E, &b_tau_E);
   fChain->SetBranchAddress("tau_Et", tau_Et, &b_tau_Et);
   fChain->SetBranchAddress("tau_M", tau_M, &b_tau_M);
   fChain->SetBranchAddress("tau_Mt", tau_Mt, &b_tau_Mt);
   fChain->SetBranchAddress("tau_Px", tau_Px, &b_tau_Px);
   fChain->SetBranchAddress("tau_Py", tau_Py, &b_tau_Py);
   fChain->SetBranchAddress("tau_Pz", tau_Pz, &b_tau_Pz);
   fChain->SetBranchAddress("tau_Eta", tau_Eta, &b_tau_Eta);
   fChain->SetBranchAddress("tau_Phi", tau_Phi, &b_tau_Phi);
   fChain->SetBranchAddress("tau_AltPt", tau_AltPt, &b_tau_AltPt);
   fChain->SetBranchAddress("tau_DecayMode", tau_DecayMode, &b_tau_DecayMode);
   fChain->SetBranchAddress("tau_Charge", tau_Charge, &b_tau_Charge);
   fChain->SetBranchAddress("tau_ParticleIso", tau_ParticleIso, &b_tau_ParticleIso);
   fChain->SetBranchAddress("tau_PFChargedHadrCands", tau_PFChargedHadrCands, &b_tau_PFChargedHadrCands);
   fChain->SetBranchAddress("tau_PFGammaCands", tau_PFGammaCands, &b_tau_PFGammaCands);
   fChain->SetBranchAddress("tau_IsolationPFChargedHadrCandsPtSum", tau_IsolationPFChargedHadrCandsPtSum, &b_tau_IsolationPFChargedHadrCandsPtSum);
   fChain->SetBranchAddress("tau_EcalStripSumEOverPLead", tau_EcalStripSumEOverPLead, &b_tau_EcalStripSumEOverPLead);
   fChain->SetBranchAddress("tau_LeadPFChargedHadrCandsignedSipt", tau_LeadPFChargedHadrCandsignedSipt, &b_tau_LeadPFChargedHadrCandsignedSipt);
   fChain->SetBranchAddress("tau_NSignalTracks", tau_NSignalTracks, &b_tau_NSignalTracks);
   fChain->SetBranchAddress("tau_PFLeadChargedPT", tau_PFLeadChargedPT, &b_tau_PFLeadChargedPT);
   fChain->SetBranchAddress("tau_Jet_pt", tau_Jet_pt, &b_tau_Jet_pt);
   fChain->SetBranchAddress("tau_Jet_eta", tau_Jet_eta, &b_tau_Jet_eta);
   fChain->SetBranchAddress("tau_Jet_phi", tau_Jet_phi, &b_tau_Jet_phi);
   fChain->SetBranchAddress("tau_Jet_m", tau_Jet_m, &b_tau_Jet_m);
   fChain->SetBranchAddress("tau_vtx_x", tau_vtx_x, &b_tau_vtx_x);
   fChain->SetBranchAddress("tau_vtx_y", tau_vtx_y, &b_tau_vtx_y);
   fChain->SetBranchAddress("tau_vtx_z", tau_vtx_z, &b_tau_vtx_z);
   fChain->SetBranchAddress("tau_id", tau_id, &b_tau_id);
   fChain->SetBranchAddress("tau_GenJet_Match_n", &tau_GenJet_Match_n, &b_tau_GenJet_Match_n);
   fChain->SetBranchAddress("tau_GenJet_DecayMode", &tau_GenJet_DecayMode, &b_tau_GenJet_DecayMode);
   fChain->SetBranchAddress("tau_GenJetMatch_Pos", &tau_GenJetMatch_Pos, &b_tau_GenJetMatch_Pos);
   fChain->SetBranchAddress("tau_GenJet_E", &tau_GenJet_E, &b_tau_GenJet_E);
   fChain->SetBranchAddress("tau_GenJet_Et", &tau_GenJet_Et, &b_tau_GenJet_Et);
   fChain->SetBranchAddress("tau_GenJet_Eta", &tau_GenJet_Eta, &b_tau_GenJet_Eta);
   fChain->SetBranchAddress("tau_GenJet_Phi", &tau_GenJet_Phi, &b_tau_GenJet_Phi);
   fChain->SetBranchAddress("tau_GenJet_M", &tau_GenJet_M, &b_tau_GenJet_M);
   fChain->SetBranchAddress("tau_GenJet_Mt", &tau_GenJet_Mt, &b_tau_GenJet_Mt);
   fChain->SetBranchAddress("tau_GenJet_P", &tau_GenJet_P, &b_tau_GenJet_P);
   fChain->SetBranchAddress("tau_GenJet_Pt", &tau_GenJet_Pt, &b_tau_GenJet_Pt);
   fChain->SetBranchAddress("tau_GenJet_Px", &tau_GenJet_Px, &b_tau_GenJet_Px);
   fChain->SetBranchAddress("tau_GenJet_Py", &tau_GenJet_Py, &b_tau_GenJet_Py);
   fChain->SetBranchAddress("tau_GenJet_Pz", &tau_GenJet_Pz, &b_tau_GenJet_Pz);
   fChain->SetBranchAddress("susyScanM0", &susyScanM0, &b_susyScanM0);
   fChain->SetBranchAddress("susyScanM12", &susyScanM12, &b_susyScanM12);
   fChain->SetBranchAddress("susyScanA0", &susyScanA0, &b_susyScanA0);
   fChain->SetBranchAddress("susyScanCrossSection", &susyScanCrossSection, &b_susyScanCrossSection);
   fChain->SetBranchAddress("susyScanMu", &susyScanMu, &b_susyScanMu);
   fChain->SetBranchAddress("susyScanRun", &susyScanRun, &b_susyScanRun);
   fChain->SetBranchAddress("susyScantanbeta", &susyScantanbeta, &b_susyScantanbeta);
   Notify();
}

Bool_t TreeContent::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TreeContent::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TreeContent::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TreeContent_cxx
