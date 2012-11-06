#ifndef TreeContent_v96_h
#define TreeContent_v96_h

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
   Double_t        pfjet_btag[100];   //[pfjet_n]
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
   TBranch        *b_trig_HLTName;   //!
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
