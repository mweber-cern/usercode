#ifndef TreeContent_h
#define TreeContent_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#define INTSIZE 4

class TreeContent {
public :
   union u64 {
     char c[INTSIZE];
     int i;
   };
   
   int get_size(const int* a) {
     int size=0;
     u64 tmp;
     for(int i=0; ; ++i ) {
       tmp.i = a[i];
       size++;
       for(int j=0; j < INTSIZE ; ++j )
         if(tmp.c[j] == '\x0')
           return size;
     }
     return -1;
   }
   
   std::string unpack(const int* a) {
     int size = get_size(a);
     u64 tmp;
     char* c = new char[size*INTSIZE];
     for(int i=0; i < size ; ++i) {
       tmp.i = a[i];
       for(int j=0; j < INTSIZE; ++j )
         c[i*INTSIZE+j] = tmp.c[j];
     }
     std::string ret( c );
     delete[] c;
     return ret;
   }
 
   
   int* pack(const char* c) {
     u64 tmp;
     int j=0,count=0;
     int size = strlen(c)/INTSIZE+1;
     int* ii=new int[size];
     for(int i=0 ; ; i++) {
       tmp.c[j++]=c[i];
       ii[count]=tmp.i;
       if(j==INTSIZE) {
         j=0;
         count++;
       }
       if(c[i] == '\x0')
         break;
     }
     return ii;
   } 

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        global_weight;
   Int_t           global_procID;
   Double_t        global_pthat;
   Double_t        global_bfield;
   Int_t           global_store;
   UInt_t          global_run;
   UInt_t          global_event;
   Int_t           global_bx;
   Int_t           global_orbit;
   Int_t           global_exp;
   Int_t           global_isdata;
   Double_t        global_rho;
   Int_t           lumi_section;
   Double_t        lumi_del;
   Double_t        lumi_rec;
   Double_t        lumi_delerr;
   Double_t        lumi_recerr;
   Int_t           pu_n;
   Int_t           pu_vtxn;
   Int_t           pu_bunchx[3];   //[pu_n]
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
   Double_t        noise_hcal_eventChargeFraction;
   Double_t        noise_hcal_eventEMEnergy;
   Double_t        noise_hcal_eventEMFraction;
   Double_t        noise_hcal_eventHadEnergy;
   Double_t        noise_hcal_eventTrackEnergy;
   Double_t        noise_hcal_flatNoiseSumE;
   Double_t        noise_hcal_flatNoiseSumEt;
   Int_t           noise_hcal_HasBadRBXTS4TS5;
   Double_t        noise_hcal_isolatedNoiseSumE;
   Double_t        noise_hcal_isolatedNoiseSumEt;
   Double_t        noise_hcal_max10GeVHitTime;
   Double_t        noise_hcal_max25GeVHitTime;
   Double_t        noise_hcal_maxE10TS;
   Double_t        noise_hcal_maxE2Over10TS;
   Double_t        noise_hcal_maxE2TS;
   Int_t           noise_hcal_maxHPDHits;
   Int_t           noise_hcal_maxHPDNoOtherHits;
   Int_t           noise_hcal_maxRBXHits;
   Int_t           noise_hcal_maxZeros;
   Double_t        noise_hcal_min10GeVHitTime;
   Double_t        noise_hcal_min25GeVHitTime;
   Double_t        noise_hcal_minE10TS;
   Double_t        noise_hcal_minE2Over10TS;
   Double_t        noise_hcal_minE2TS;
   Double_t        noise_hcal_minHPDEMF;
   Double_t        noise_hcal_minRBXEMF;
   Int_t           noise_hcal_noiseFilterStatus;
   Int_t           noise_hcal_noiseType;
   Int_t           noise_hcal_num10GeVHits;
   Int_t           noise_hcal_num25GeVHits;
   Int_t           noise_hcal_numFlatNoiseChannels;
   Int_t           noise_hcal_numIsolatedNoiseChannels;
   Int_t           noise_hcal_numProblematicRBXs;
   Int_t           noise_hcal_numSpikeNoiseChannels;
   Int_t           noise_hcal_numTriangleNoiseChannels;
   Int_t           noise_hcal_numTS4TS5NoiseChannels;
   Int_t           noise_hcal_passHighLevelNoiseFilter;
   Int_t           noise_hcal_passLooseNoiseFilter;
   Int_t           noise_hcal_passTightNoiseFilter;
   Double_t        noise_hcal_rms10GeVHitTime;
   Double_t        noise_hcal_rms25GeVHitTime;
   Double_t        noise_hcal_spikeNoiseSumE;
   Double_t        noise_hcal_spikeNoiseSumEt;
   Double_t        noise_hcal_triangleNoiseSumE;
   Double_t        noise_hcal_triangleNoiseSumEt;
   Double_t        noise_hcal_TS4TS5NoiseSumE;
   Double_t        noise_hcal_TS4TS5NoiseSumEt;
   Bool_t          noise_HBHE_filter_result;
   Int_t           trig_HLTName[20];
   Int_t           trig_n;
   Int_t           trig_L1prescale[7000];   //[trig_n]
   Int_t           trig_HLTprescale[7000];   //[trig_n]
   Int_t           trig_name[7000][20];   //[trig_n]
   Int_t           trig_filter[7000][20];   //[trig_n]
   Double_t        trig_pt[7000];   //[trig_n]
   Double_t        trig_eta[7000];   //[trig_n]
   Double_t        trig_phi[7000];   //[trig_n]
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
   Double_t        met_et[8];
   Double_t        met_ex[8];
   Double_t        met_ey[8];
   Double_t        met_phi[8];
   Double_t        met_sumet[8];
   Double_t        met_sumetsig[8];
   Double_t        met_etsignif[8];
   Double_t        met_CaloMETInmHF[8];
   Double_t        met_CaloMETInpHF[8];
   Double_t        met_CaloMETPhiInmHF[8];
   Double_t        met_CaloMETPhiInpHF[8];
   Double_t        met_CaloSETInmHF[8];
   Double_t        met_CaloSETInpHF[8];
   Double_t        met_emEtFraction[8];
   Double_t        met_etFractionHadronic[8];
   Double_t        met_maxEtInEmTowers[8];
   Double_t        met_maxEtInHadTowers[8];
   Double_t        met_emEtInHF[8];
   Double_t        met_emEtInEE[8];
   Double_t        met_emEtInEB[8];
   Double_t        met_hadEtInHF[8];
   Double_t        met_hadEtInHE[8];
   Double_t        met_hadEtInHO[8];
   Double_t        met_hadEtInHB[8];
   Double_t        met_ChargedEMEtFraction[8];
   Double_t        met_ChargedHadEtFraction[8];
   Double_t        met_MuonEtFraction[8];
   Double_t        met_NeutralEMFraction[8];
   Double_t        met_NeutralHadEtFraction[8];
   Double_t        met_Type6EtFraction[8];
   Double_t        met_Type7EtFraction[8];
   Int_t           calojet_n;
   Double_t        calojet_E[100];   //[calojet_n]
   Double_t        calojet_Et[100];   //[calojet_n]
   Double_t        calojet_p[100];   //[calojet_n]
   Double_t        calojet_pt[100];   //[calojet_n]
   Double_t        calojet_pt_raw[100];   //[calojet_n]
   Double_t        calojet_px[100];   //[calojet_n]
   Double_t        calojet_py[100];   //[calojet_n]
   Double_t        calojet_pz[100];   //[calojet_n]
   Double_t        calojet_eta[100];   //[calojet_n]
   Double_t        calojet_phi[100];   //[calojet_n]
   Double_t        calojet_fem[100];   //[calojet_n]
   Double_t        calojet_fhad[100];   //[calojet_n]
   Double_t        calojet_btag[100];   //[calojet_n]
   Double_t        calojet_charge[100];   //[calojet_n]
   Double_t        calojet_fHPD[100];   //[calojet_n]
   Double_t        calojet_fRBX[100];   //[calojet_n]
   Int_t           calojet_n90hits[100];   //[calojet_n]
   Int_t           calojet_n90[100];   //[calojet_n]
   Int_t           calojet_flav[100];   //[calojet_n]
   Int_t           calojet_truth[100];   //[calojet_n]
   Int_t           calojet_const[100];   //[calojet_n]
   Int_t           calojet_ID[100];   //[calojet_n]
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
   Int_t           fatjet_n;
   Int_t           fatjet_nsub[100];   //[fatjet_n]
   Double_t        fatjet_pt[100];   //[fatjet_n]
   Double_t        fatjet_px[100];   //[fatjet_n]
   Double_t        fatjet_py[100];   //[fatjet_n]
   Double_t        fatjet_pz[100];   //[fatjet_n]
   Double_t        fatjet_E[100];   //[fatjet_n]
   Double_t        fatjet_eta[100];   //[fatjet_n]
   Double_t        fatjet_phi[100];   //[fatjet_n]
   Double_t        fatjet_sub_pt[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_px[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_py[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_pz[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_E[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_eta[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_phi[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_fem[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_fhad[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_btag[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_n90[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_fHPD[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_fRBX[100][10];   //[fatjet_n]
   Int_t           SC_n;
   Int_t           SC_truth[200];   //[SC_n]
   Double_t        SC_E[200];   //[SC_n]
   Double_t        SC_phi[200];   //[SC_n]
   Double_t        SC_eta[200];   //[SC_n]
   Int_t           SC_trign[200];  //[SC_n]
   Int_t           SC_trig[200][500]; //[SC_n][500]
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
   Double_t        ele_d0bs[100];   //[ele_n]
   Double_t        ele_sd0[100];   //[ele_n]
   Int_t           ele_hits[100];   //[ele_n]
   Int_t           ele_truth[100];   //[ele_n]
   Int_t           ele_isECal[100];   //[ele_n]
   Int_t           ele_isTracker[100];   //[ele_n]
   Int_t           ele_ValidHitFirstPxlB[100];   //[ele_n]
   Int_t           ele_TrkExpHitsInner[100];   //[ele_n]
   Double_t        ele_HCalOverEm[100];   //[ele_n]
   Double_t        ele_Dr03TkSumPt[100];   //[ele_n]
   Double_t        ele_Dr04HCalSumEt[100];   //[ele_n]
   Double_t        ele_Dr03HCalSumEt[100];   //[ele_n]
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
   Int_t           ele_trign[100];   //[ele_n]
   Int_t           ele_trig[100][500];   //[ele_n]
   Int_t           ele_SC[100];   //[ele_n]
   Double_t        ele_PFiso[100][9];      //[ele_n]
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
   Double_t        ele_SwissCross[100];   //[ele_n]
   Double_t        ele_EoverP[100];   //[ele_n]
   Int_t           pfele_n;
   Double_t        pfele_p[100];   //[pfele_n]
   Double_t        pfele_E[100];   //[pfele_n]
   Double_t        pfele_Et[100];   //[pfele_n]
   Double_t        pfele_CaloEt[100];   //[pfele_n]
   Double_t        pfele_pt[100];   //[pfele_n]
   Double_t        pfele_TrackptError[100];   //[ele_n]
   Double_t        pfele_Trackpt[100];   //[ele_n]
   Double_t        pfele_px[100];   //[pfele_n]
   Double_t        pfele_py[100];   //[pfele_n]
   Double_t        pfele_pz[100];   //[pfele_n]
   Double_t        pfele_eta[100];   //[pfele_n]
   Double_t        pfele_phi[100];   //[pfele_n]
   Double_t        pfele_charge[100];   //[pfele_n]
   Int_t           pfele_truth[100];   //[pfele_n]
   Int_t           pfele_trign[100];   //[pfele_n]
   Int_t           pfele_trig[100][500];   //[pfele_n]
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
   Double_t        muo_calocomp[100];   //[muo_n]
   Double_t        muo_calotower_e[100];   //[muo_n]
   Int_t           muo_prompttight[100];   //[muo_n]
   Int_t           muo_hitsTk[100];   //[muo_n]
   Int_t           muo_truth[100];   //[muo_n]
   Int_t           muo_trign[100];   //[muo_n]
   Int_t           muo_trig[100][500];   //[muo_n]
   Int_t           muo_ID[100][24];   //[muo_n]
   Int_t           muo_ChambersMatched[100];   //[muo_n]
   Double_t        muo_Valid_fraction[100];   //[muo_n]
   Double_t        muo_TrkChiNormCm[100];   //[muo_n]
   Int_t           muo_hitsCm[100];   //[muo_n]
   Double_t        muo_d0Cm[100];   //[muo_n]
   Double_t        muo_sd0Cm[100];   //[muo_n]
   Double_t        muo_d0OriginCm[100];   //[muo_n]
   Double_t        muo_d0bsCm[100];   //[muo_n]
   Double_t        muo_dzbsCm[100];   //[muo_n]
   Double_t        muo_vx[100];   //[muo_n]
   Double_t        muo_vy[100];   //[muo_n]
   Double_t        muo_vz[100];   //[muo_n]
   Int_t           muo_ValidMuonHitsCm[100];   //[muo_n]
   Int_t           muo_ValidTrackerHitsCm[100];   //[muo_n]
   Int_t           muo_ValidPixelHitsCm[100];   //[muo_n]
   Int_t           muo_TrackerLayersMeasCm[100];   //[muo_n]
   Int_t           muo_TrackerLayersNotMeasCm[100];   //[muo_n]
   Int_t           muo_LostHits[100];   //[muo_n]
   Double_t        muo_Cocktail_pt[100];   //[muo_n]
   Double_t        muo_Cocktail_phi[100];   //[muo_n]
   Double_t        muo_Cocktail_eta[100];   //[muo_n]
   Double_t        muo_TevReco_pt[100][7];   //[muo_n]
   Double_t        muo_TevReco_ptError[100][7];   //[muo_n]
   Double_t        muo_TevReco_eta[100][7];   //[muo_n]
   Double_t        muo_TevReco_phi[100][7];   //[muo_n]
   Double_t        muo_TevReco_chi2[100][7];   //[muo_n]
   Double_t        muo_TevReco_ndof[100][7];   //[muo_n]
   Double_t        muo_PFiso[100][9];      //[muo_n]
   Double_t        muo_PFCand_px[100];   //[muo_n]
   Double_t        muo_PFCand_py[100];   //[muo_n]
   Double_t        muo_PFCand_pz[100];   //[muo_n]
   Double_t        muo_PFCand_E[100];   //[muo_n]
   Double_t        muo_PFCand_eta[100];   //[muo_n]
   Double_t        muo_PFCand_phi[100];   //[muo_n]
   Int_t           muo_PFCand_pfid[100];   //[muo_n]
   Double_t        muo_PFCand_DeltaR[100];   //[muo_n]
   Int_t           PFmuo_n;
   Double_t        PFmuo_p[100];   //[PFmuo_n]
   Double_t        PFmuo_pt[100];   //[PFmuo_n]
   Double_t        PFmuo_E[100];   //[PFmuo_n]
   Double_t        PFmuo_Et[100];   //[PFmuo_n]
   Double_t        PFmuo_px[100];   //[PFmuo_n]
   Double_t        PFmuo_py[100];   //[PFmuo_n]
   Double_t        PFmuo_pz[100];   //[PFmuo_n]
   Double_t        PFmuo_eta[100];   //[PFmuo_n]
   Double_t        PFmuo_phi[100];   //[PFmuo_n]
   Double_t        PFmuo_Charge[100];   //[PFmuo_n]
   Double_t        PFmuo_particleIso[100];   //[PFmuo_n]
   Double_t        PFmuo_chadIso[100];   //[PFmuo_n]
   Double_t        PFmuo_nhadIso[100];   //[PFmuo_n]
   Double_t        PFmuo_gamIso[100];   //[PFmuo_n]
   Double_t        PFmuo_RelTrkIso[100];   //[PFmuo_n]
   Double_t        PFmuo_TrkIso[100];   //[PFmuo_n]
   Double_t        PFmuo_ECalIso[100];   //[PFmuo_n]
   Double_t        PFmuo_HCalIso[100];   //[PFmuo_n]
   Double_t        PFmuo_TrkIsoDep[100];   //[PFmuo_n]
   Double_t        PFmuo_ECalIsoDep[100];   //[PFmuo_n]
   Double_t        PFmuo_HCalIsoDep[100];   //[PFmuo_n]
   Double_t        PFmuo_AllIso[100];   //[PFmuo_n]
   Double_t        PFmuo_TrkChiNormCm[100];   //[PFmuo_n]
   Double_t        PFmuo_TrkChiNormTk[100];   //[PFmuo_n]
   Double_t        PFmuo_d0Cm[100];   //[PFmuo_n]
   Double_t        PFmuo_d0Tk[100];   //[PFmuo_n]
   Double_t        PFmuo_sd0Cm[100];   //[PFmuo_n]
   Double_t        PFmuo_sd0Tk[100];   //[PFmuo_n]
   Double_t        PFmuo_calocomp[100];   //[PFmuo_n]
   Double_t        PFmuo_calotower_e[100];   //[PFmuo_n]
   Int_t           PFmuo_hitsCm[100];   //[PFmuo_n]
   Int_t           PFmuo_hitsTk[100];   //[PFmuo_n]
   Int_t           PFmuo_ValidMuonHitsCm[100];   //[PFmuo_n]
   Int_t           PFmuo_ValidTrackerHitsCm[100];   //[PFmuo_n]
   Int_t           PFmuo_ValidPixelHitsCm[100];   //[PFmuo_n]
   Int_t           PFmuo_ChambersMatched[100];   //[PFmuo_n]
   Double_t        PFmuo_d0bsCm[100];   //[PFmuo_n]
   Double_t        PFmuo_d0OriginCm[100];   //[PFmuo_n]
   Double_t        PFmuo_dzbsCm[100];   //[PFmuo_n]
   Double_t        PFmuo_vx[100];   //[PFmuo_n]
   Double_t        PFmuo_vy[100];   //[PFmuo_n]
   Double_t        PFmuo_vz[100];   //[PFmuo_n]
   Int_t           PFmuo_TrackerLayersMeasCm[100];   //[PFmuo_n]
   Int_t           PFmuo_TrackerLayersNotMeasCm[100];   //[PFmuo_n]
   Double_t        PFmuo_Valid_fraction[100];   //[PFmuo_n]
   Int_t           tau_n;
   Double_t        tau_p[100];   //[tau_n]
   Double_t        tau_pt[100];   //[tau_n]
   Double_t        tau_E[100];   //[tau_n]
   Double_t        tau_Et[100];   //[tau_n]
   Double_t        tau_Px[100];   //[tau_n]
   Double_t        tau_Py[100];   //[tau_n]
   Double_t        tau_Pz[100];   //[tau_n]
   Double_t        tau_Eta[100];   //[tau_n]
   Double_t        tau_Phi[100];   //[tau_n]
   Int_t           tau_DecayMode[100];   //[tau_n]
   Double_t        tau_vx[100];   //[tau_n]
   Double_t        tau_vy[100];   //[tau_n]
   Double_t        tau_vz[100];   //[tau_n]
   Double_t        tau_vx2[100];   //[tau_n]
   Double_t        tau_vy2[100];   //[tau_n]
   Double_t        tau_vz2[100];   //[tau_n]
   Int_t           tau_trign[100];   //[tau_n]
   Int_t           tau_trig[100][500];   //[tau_n]
   Double_t        tau_ParticleIso[100];   //[tau_n]
   Double_t        tau_ChadIso[100];   //[tau_n]
   Double_t        tau_NhadIso[100];   //[tau_n]
   Double_t        tau_GamIso[100];   //[tau_n]
   Int_t           tau_PFChargedHadrCands[100];   //[tau_n]
   Int_t           tau_PFGammaCands[100];   //[tau_n]
   Double_t        tau_IsolationPFChargedHadrCandsPtSum[100];   //[tau_n]
   Double_t        tau_IsolationPFGammaCandsEtSum[100];   //[tau_n]
   Double_t        tau_EcalStripSumEOverPLead[100];   //[tau_n]
   Double_t        tau_EMfraction[100];   //[tau_n]
   Double_t        tau_Hcal3x3OverPLead[100];   //[tau_n]
   Double_t        tau_HcalMaxOverPLead[100];   //[tau_n]
   Double_t        tau_HcalTotOverPLead[100];   //[tau_n]
   Double_t        tau_LeadPFChargedHadrCandsignedSipt[100];   //[tau_n]
   Double_t        tau_PhiPhiMoment[100];   //[tau_n]
   Double_t        tau_EtaPhiMoment[100];   //[tau_n]
   Double_t        tau_EtaEtaMoment[100];   //[tau_n]
   Int_t           tau_NSignalTracks[100];   //[tau_n]
   Double_t        tau_ElectronPreIDOutput[100];   //[tau_n]
   Double_t        tau_PFLeadChargedPT[100];   //[tau_n]
   Double_t        tau_BremsRecoveryEOverPLead[100];   //[tau_n]
   Double_t        tau_id[100][10];   //[tau_n]
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
   TBranch        *b_global_pthat;   //!
   TBranch        *b_global_bfield;   //!
   TBranch        *b_global_store;   //!
   TBranch        *b_global_run;   //!
   TBranch        *b_global_event;   //!
   TBranch        *b_global_bx;   //!
   TBranch        *b_global_orbit;   //!
   TBranch        *b_global_exp;   //!
   TBranch        *b_global_isdata;   //!
   TBranch        *b_global_rho;   //!
   TBranch        *b_lumi_section;   //!
   TBranch        *b_lumi_del;   //!
   TBranch        *b_lumi_rec;   //!
   TBranch        *b_lumi_delerr;   //!
   TBranch        *b_lumi_recerr;   //!
   TBranch        *b_pu_n;   //!
   TBranch        *b_pu_vtxn;   //!
   TBranch        *b_pu_bunchx;   //!
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
   TBranch        *b_noise_hcal_eventChargeFraction;   //!
   TBranch        *b_noise_hcal_eventEMEnergy;   //!
   TBranch        *b_noise_hcal_eventEMFraction;   //!
   TBranch        *b_noise_hcal_eventHadEnergy;   //!
   TBranch        *b_noise_hcal_eventTrackEnergy;   //!
   TBranch        *b_noise_hcal_flatNoiseSumE;   //!
   TBranch        *b_noise_hcal_flatNoiseSumEt;   //!
   TBranch        *b_noise_hcal_HasBadRBXTS4TS5;   //!
   TBranch        *b_noise_hcal_isolatedNoiseSumE;   //!
   TBranch        *b_noise_hcal_isolatedNoiseSumEt;   //!
   TBranch        *b_noise_hcal_max10GeVHitTime;   //!
   TBranch        *b_noise_hcal_max25GeVHitTime;   //!
   TBranch        *b_noise_hcal_maxE10TS;   //!
   TBranch        *b_noise_hcal_maxE2Over10TS;   //!
   TBranch        *b_noise_hcal_maxE2TS;   //!
   TBranch        *b_noise_hcal_maxHPDHits;   //!
   TBranch        *b_noise_hcal_maxHPDNoOtherHits;   //!
   TBranch        *b_noise_hcal_maxRBXHits;   //!
   TBranch        *b_noise_hcal_maxZeros;   //!
   TBranch        *b_noise_hcal_min10GeVHitTime;   //!
   TBranch        *b_noise_hcal_min25GeVHitTime;   //!
   TBranch        *b_noise_hcal_minE10TS;   //!
   TBranch        *b_noise_hcal_minE2Over10TS;   //!
   TBranch        *b_noise_hcal_minE2TS;   //!
   TBranch        *b_noise_hcal_minHPDEMF;   //!
   TBranch        *b_noise_hcal_minRBXEMF;   //!
   TBranch        *b_noise_hcal_noiseFilterStatus;   //!
   TBranch        *b_noise_hcal_noiseType;   //!
   TBranch        *b_noise_hcal_num10GeVHits;   //!
   TBranch        *b_noise_hcal_num25GeVHits;   //!
   TBranch        *b_noise_hcal_numFlatNoiseChannels;   //!
   TBranch        *b_noise_hcal_numIsolatedNoiseChannels;   //!
   TBranch        *b_noise_hcal_numProblematicRBXs;   //!
   TBranch        *b_noise_hcal_numSpikeNoiseChannels;   //!
   TBranch        *b_noise_hcal_numTriangleNoiseChannels;   //!
   TBranch        *b_noise_hcal_numTS4TS5NoiseChannels;   //!
   TBranch        *b_noise_hcal_passHighLevelNoiseFilter;   //!
   TBranch        *b_noise_hcal_passLooseNoiseFilter;   //!
   TBranch        *b_noise_hcal_passTightNoiseFilter;   //!
   TBranch        *b_noise_hcal_rms10GeVHitTime;   //!
   TBranch        *b_noise_hcal_rms25GeVHitTime;   //!
   TBranch        *b_noise_hcal_spikeNoiseSumE;   //!
   TBranch        *b_noise_hcal_spikeNoiseSumEt;   //!
   TBranch        *b_noise_hcal_triangleNoiseSumE;   //!
   TBranch        *b_noise_hcal_triangleNoiseSumEt;   //!
   TBranch        *b_noise_hcal_TS4TS5NoiseSumE;   //!
   TBranch        *b_noise_hcal_TS4TS5NoiseSumEt;   //!
   TBranch        *b_noise_HBHE_filter_result;   //!
   TBranch        *b_trig_HLTName;   //!
   TBranch        *b_trig_n;   //!
   TBranch        *b_trig_L1prescale;   //!
   TBranch        *b_trig_HLTprescale;   //!
   TBranch        *b_trig_name;   //!
   TBranch        *b_trig_filter;   //!
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
   TBranch        *b_met_CaloMETInmHF;   //!
   TBranch        *b_met_CaloMETInpHF;   //!
   TBranch        *b_met_CaloMETPhiInmHF;   //!
   TBranch        *b_met_CaloMETPhiInpHF;   //!
   TBranch        *b_met_CaloSETInmHF;   //!
   TBranch        *b_met_CaloSETInpHF;   //!
   TBranch        *b_met_emEtFraction;   //!
   TBranch        *b_met_etFractionHadronic;   //!
   TBranch        *b_met_maxEtInEmTowers;   //!
   TBranch        *b_met_maxEtInHadTowers;   //!
   TBranch        *b_met_emEtInHF;   //!
   TBranch        *b_met_emEtInEE;   //!
   TBranch        *b_met_emEtInEB;   //!
   TBranch        *b_met_hadEtInHF;   //!
   TBranch        *b_met_hadEtInHE;   //!
   TBranch        *b_met_hadEtInHO;   //!
   TBranch        *b_met_hadEtInHB;   //!
   TBranch        *b_met_ChargedEMEtFraction;   //!
   TBranch        *b_met_ChargedHadEtFraction;   //!
   TBranch        *b_met_MuonEtFraction;   //!
   TBranch        *b_met_NeutralEMFraction;   //!
   TBranch        *b_met_NeutralHadEtFraction;   //!
   TBranch        *b_met_Type6EtFraction;   //!
   TBranch        *b_met_Type7EtFraction;   //!
   TBranch        *b_calojet_n;   //!
   TBranch        *b_calojet_E;   //!
   TBranch        *b_calojet_Et;   //!
   TBranch        *b_calojet_p;   //!
   TBranch        *b_calojet_pt;   //!
   TBranch        *b_calojet_pt_raw;   //!
   TBranch        *b_calojet_px;   //!
   TBranch        *b_calojet_py;   //!
   TBranch        *b_calojet_pz;   //!
   TBranch        *b_calojet_eta;   //!
   TBranch        *b_calojet_phi;   //!
   TBranch        *b_calojet_fem;   //!
   TBranch        *b_calojet_fhad;   //!
   TBranch        *b_calojet_btag;   //!
   TBranch        *b_calojet_charge;   //!
   TBranch        *b_calojet_fHPD;   //!
   TBranch        *b_calojet_fRBX;   //!
   TBranch        *b_calojet_n90hits;   //!
   TBranch        *b_calojet_n90;   //!
   TBranch        *b_calojet_flav;   //!
   TBranch        *b_calojet_truth;   //!
   TBranch        *b_calojet_const;   //!
   TBranch        *b_calojet_ID;   //!
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
   TBranch        *b_fatjet_n;   //!
   TBranch        *b_fatjet_nsub;   //!
   TBranch        *b_fatjet_pt;   //!
   TBranch        *b_fatjet_px;   //!
   TBranch        *b_fatjet_py;   //!
   TBranch        *b_fatjet_pz;   //!
   TBranch        *b_fatjet_E;   //!
   TBranch        *b_fatjet_eta;   //!
   TBranch        *b_fatjet_phi;   //!
   TBranch        *b_fatjet_sub_pt;   //!
   TBranch        *b_fatjet_sub_px;   //!
   TBranch        *b_fatjet_sub_py;   //!
   TBranch        *b_fatjet_sub_pz;   //!
   TBranch        *b_fatjet_sub_E;   //!
   TBranch        *b_fatjet_sub_eta;   //!
   TBranch        *b_fatjet_sub_phi;   //!
   TBranch        *b_fatjet_sub_fem;   //!
   TBranch        *b_fatjet_sub_fhad;   //!
   TBranch        *b_fatjet_sub_btag;   //!
   TBranch        *b_fatjet_sub_n90;   //!
   TBranch        *b_fatjet_sub_fHPD;   //!
   TBranch        *b_fatjet_sub_fRBX;   //!
   TBranch        *b_SC_n;   //!
   TBranch        *b_SC_truth;   //!
   TBranch        *b_SC_E;   //!
   TBranch        *b_SC_phi;   //!
   TBranch        *b_SC_eta;   //!
   TBranch        *b_SC_trign;   //!
   TBranch        *b_SC_trig;   //!
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
   TBranch        *b_ele_d0bs;   //!
   TBranch        *b_ele_sd0;   //!
   TBranch        *b_ele_hits;   //!
   TBranch        *b_ele_truth;   //!
   TBranch        *b_ele_isECal;   //!
   TBranch        *b_ele_isTracker;   //!
   TBranch        *b_ele_ValidHitFirstPxlB;   //!
   TBranch        *b_ele_TrkExpHitsInner;   //!
   TBranch        *b_ele_HCalOverEm;   //!
   TBranch        *b_ele_Dr03TkSumPt;   //!
   TBranch        *b_ele_Dr04HCalSumEt;   //!
   TBranch        *b_ele_Dr03HCalSumEt;   //!
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
   TBranch        *b_ele_trign;   //!
   TBranch        *b_ele_trig;   //!
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
   TBranch        *b_ele_SwissCross;   //!
   TBranch        *b_ele_EoverP;   //!
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
   TBranch        *b_pfele_trign;   //!
   TBranch        *b_pfele_trig;   //!
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
   TBranch        *b_muo_calocomp;   //!
   TBranch        *b_muo_calotower_e;   //!
   TBranch        *b_muo_prompttight;   //!
   TBranch        *b_muo_hitsTk;   //!
   TBranch        *b_muo_truth;   //!
   TBranch        *b_muo_trign;   //!
   TBranch        *b_muo_trig;   //!
   TBranch        *b_muo_ID;   //!
   TBranch        *b_muo_ChambersMatched;   //!
   TBranch        *b_muo_Valid_fraction;   //!
   TBranch        *b_muo_TrkChiNormCm;   //!
   TBranch        *b_muo_hitsCm;   //!
   TBranch        *b_muo_d0Cm;   //!
   TBranch        *b_muo_sd0Cm;   //!
   TBranch        *b_muo_d0OriginCm;   //!
   TBranch        *b_muo_d0bsCm;   //!
   TBranch        *b_muo_dzbsCm;   //!
   TBranch        *b_muo_vx;   //!
   TBranch        *b_muo_vy;   //!
   TBranch        *b_muo_vz;   //!
   TBranch        *b_muo_ValidMuonHitsCm;   //!
   TBranch        *b_muo_ValidTrackerHitsCm;   //!
   TBranch        *b_muo_ValidPixelHitsCm;   //!
   TBranch        *b_muo_TrackerLayersMeasCm;   //!
   TBranch        *b_muo_TrackerLayersNotMeasCm;   //!
   TBranch        *b_muo_LostHits;   //!
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
   TBranch        *b_PFmuo_n;   //!
   TBranch        *b_PFmuo_p;   //!
   TBranch        *b_PFmuo_pt;   //!
   TBranch        *b_PFmuo_E;   //!
   TBranch        *b_PFmuo_Et;   //!
   TBranch        *b_PFmuo_px;   //!
   TBranch        *b_PFmuo_py;   //!
   TBranch        *b_PFmuo_pz;   //!
   TBranch        *b_PFmuo_eta;   //!
   TBranch        *b_PFmuo_phi;   //!
   TBranch        *b_PFmuo_Charge;   //!
   TBranch        *b_PFmuo_particleIso;   //!
   TBranch        *b_PFmuo_chadIso;   //!
   TBranch        *b_PFmuo_nhadIso;   //!
   TBranch        *b_PFmuo_gamIso;   //!
   TBranch        *b_PFmuo_RelTrkIso;   //!
   TBranch        *b_PFmuo_TrkIso;   //!
   TBranch        *b_PFmuo_ECalIso;   //!
   TBranch        *b_PFmuo_HCalIso;   //!
   TBranch        *b_PFmuo_TrkIsoDep;   //!
   TBranch        *b_PFmuo_ECalIsoDep;   //!
   TBranch        *b_PFmuo_HCalIsoDep;   //!
   TBranch        *b_PFmuo_AllIso;   //!
   TBranch        *b_PFmuo_TrkChiNormCm;   //!
   TBranch        *b_PFmuo_TrkChiNormTk;   //!
   TBranch        *b_PFmuo_d0Cm;   //!
   TBranch        *b_PFmuo_d0Tk;   //!
   TBranch        *b_PFmuo_sd0Cm;   //!
   TBranch        *b_PFmuo_sd0Tk;   //!
   TBranch        *b_PFmuo_calocomp;   //!
   TBranch        *b_PFmuo_calotower_e;   //!
   TBranch        *b_PFmuo_hitsCm;   //!
   TBranch        *b_PFmuo_hitsTk;   //!
   TBranch        *b_PFmuo_ValidMuonHitsCm;   //!
   TBranch        *b_PFmuo_ValidTrackerHitsCm;   //!
   TBranch        *b_PFmuo_ValidPixelHitsCm;   //!
   TBranch        *b_PFmuo_ChambersMatched;   //!
   TBranch        *b_PFmuo_d0bsCm;   //!
   TBranch        *b_PFmuo_d0OriginCm;   //!
   TBranch        *b_PFmuo_dzbsCm;   //!
   TBranch        *b_PFmuo_vx;   //!
   TBranch        *b_PFmuo_vy;   //!
   TBranch        *b_PFmuo_vz;   //!
   TBranch        *b_PFmuo_TrackerLayersMeasCm;   //!
   TBranch        *b_PFmuo_TrackerLayersNotMeasCm;   //!
   TBranch        *b_PFmuo_Valid_fraction;   //!
   TBranch        *b_tau_n;   //!
   TBranch        *b_tau_p;   //!
   TBranch        *b_tau_pt;   //!
   TBranch        *b_tau_E;   //!
   TBranch        *b_tau_Et;   //!
   TBranch        *b_tau_Px;   //!
   TBranch        *b_tau_Py;   //!
   TBranch        *b_tau_Pz;   //!
   TBranch        *b_tau_Eta;   //!
   TBranch        *b_tau_Phi;   //!
   TBranch        *b_tau_DecayMode;   //!
   TBranch        *b_tau_vx;   //!
   TBranch        *b_tau_vy;   //!
   TBranch        *b_tau_vz;   //!
   TBranch        *b_tau_vx2;   //!
   TBranch        *b_tau_vy2;   //!
   TBranch        *b_tau_vz2;   //!
   TBranch        *b_tau_trign;   //!
   TBranch        *b_tau_trig;   //!
   TBranch        *b_tau_ParticleIso;   //!
   TBranch        *b_tau_ChadIso;   //!
   TBranch        *b_tau_NhadIso;   //!
   TBranch        *b_tau_GamIso;   //!
   TBranch        *b_tau_PFChargedHadrCands;   //!
   TBranch        *b_tau_PFGammaCands;   //!
   TBranch        *b_tau_IsolationPFChargedHadrCandsPtSum;   //!
   TBranch        *b_tau_IsolationPFGammaCandsEtSum;   //!
   TBranch        *b_tau_EcalStripSumEOverPLead;   //!
   TBranch        *b_tau_EMfraction;   //!
   TBranch        *b_tau_Hcal3x3OverPLead;   //!
   TBranch        *b_tau_HcalMaxOverPLead;   //!
   TBranch        *b_tau_HcalTotOverPLead;   //!
   TBranch        *b_tau_LeadPFChargedHadrCandsignedSipt;   //!
   TBranch        *b_tau_PhiPhiMoment;   //!
   TBranch        *b_tau_EtaPhiMoment;   //!
   TBranch        *b_tau_EtaEtaMoment;   //!
   TBranch        *b_tau_NSignalTracks;   //!
   TBranch        *b_tau_ElectronPreIDOutput;   //!
   TBranch        *b_tau_PFLeadChargedPT;   //!
   TBranch        *b_tau_BremsRecoveryEOverPLead;   //!
   TBranch        *b_tau_id;   //!
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

#endif /* TreeContent_h */
