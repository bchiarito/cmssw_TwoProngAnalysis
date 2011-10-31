//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug 25 13:11:30 2010 by ROOT version 5.26/00
// from TTree fTree/PhotonTree
// found on file: diphotonTree_PhotonJet_Pt30to50.root
//////////////////////////////////////////////////////////

#ifndef fTree_h
#define fTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"

#include <iomanip>
#include <iostream>
#include <math.h>
#include <string.h>


class fTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Event_run;
   Int_t           Event_LS;
   Int_t           Event_evnum;
   Int_t           GenEvent_signalProcessId;
   Double_t        GenEvent_binningValue;
   Double_t        Vtx_vx;
   Double_t        Vtx_vy;
   Double_t        Vtx_vz;
   Double_t        Vtx_sumPtTracks;
   Double_t        Vtx_ndof;
   Double_t        Vtx_d0;
   Int_t           Vtx_Nvtx;
   Int_t           Vtx_Ntracks;
   Bool_t          Vtx_isFake;
   Double_t        Vtx2_vx;
   Double_t        Vtx2_vy;
   Double_t        Vtx2_vz;
   Double_t        Vtx2_sumPtTracks;
   Double_t        Vtx2_ndof;
   Double_t        Vtx2_d0;
   Int_t           Vtx2_Nvtx;
   Int_t           Vtx2_Ntracks;
   Bool_t          Vtx2_isFake;
   Double_t        Vtx3_vx;
   Double_t        Vtx3_vy;
   Double_t        Vtx3_vz;
   Double_t        Vtx3_sumPtTracks;
   Double_t        Vtx3_ndof;
   Double_t        Vtx3_d0;
   Int_t           Vtx3_Nvtx;
   Int_t           Vtx3_Ntracks;
   Bool_t          Vtx3_isFake;
   Double_t        VtxGEN_vx;
   Double_t        VtxGEN_vy;
   Double_t        VtxGEN_vz;
   Double_t        VtxGEN_sumPtTracks;
   Double_t        VtxGEN_ndof;
   Double_t        VtxGEN_d0;
   Int_t           VtxGEN_Nvtx;
   Int_t           VtxGEN_Ntracks;
   Bool_t          VtxGEN_isFake;
   Double_t        rho;
   Double_t           pu_n;
   Double_t        BeamSpot_x0;
   Double_t        BeamSpot_y0;
   Double_t        BeamSpot_z0;
   Double_t        BeamSpot_sigmaZ;
   Double_t        BeamSpot_x0error;
   Double_t        BeamSpot_y0error;
   Double_t        BeamSpot_z0error;
   Double_t        BeamSpot_sigmaZ0error;
   Bool_t          L1trg_L1_Tech0;
   Bool_t          L1trg_L1_Tech36;
   Bool_t          L1trg_L1_Tech37;
   Bool_t          L1trg_L1_Tech38;
   Bool_t          L1trg_L1_Tech39;
   Bool_t          L1trg_L1_Tech40;
   Bool_t          L1trg_L1_Tech41;
   Bool_t          L1trg_L1_Tech42;
   Bool_t          L1trg_L1_Tech43;
   Bool_t          L1trg_L1_EG2;
   Bool_t          L1trg_L1_EG5;
   Bool_t          L1trg_L1_EG8;
   Int_t           TrigHLT_HLT_MinBiasBSC;
   Int_t           TrigHLT_HLT_MinBiasBSC_NoBPTX;
   Int_t           TrigHLT_HLT_MinBiasBSC_OR;
   Int_t           TrigHLT_HLT_L1_BscMinBiasOR_BptxPlusORMinus;
   Int_t           TrigHLT_HLT_L1SingleEG2;
   Int_t           TrigHLT_HLT_L1SingleEG5;
   Int_t           TrigHLT_HLT_L1SingleEG8;
   Int_t           TrigHLT_HLT_L1DoubleEG5;
   Int_t           TrigHLT_HLT_Photon10_L1R;
   Int_t           TrigHLT_HLT_Photon10_Cleaned_L1R;
   Int_t           TrigHLT_HLT_Photon15_L1R;
   Int_t           TrigHLT_HLT_Photon15_Cleaned_L1R;
   Int_t           TrigHLT_HLT_Photon15_LooseEcalIso_L1R;
   Int_t           TrigHLT_HLT_Photon15_LooseEcalIso_Cleaned_L1R;
   Int_t           TrigHLT_HLT_Photon15_TrackIso_L1R;
   Int_t           TrigHLT_HLT_Photon15_TrackIso_Cleaned_L1R;
   Int_t           TrigHLT_HLT_Photon17_Isol_SC17HE_L1R_v1;
   Int_t           TrigHLT_HLT_Photon17_SC17HE_L1R_v1;
   Int_t           TrigHLT_HLT_Photon20_L1R;
   Int_t           TrigHLT_HLT_Photon20_Cleaned_L1R;
   Int_t           TrigHLT_HLT_Photon20_NoHE_L1R;
   Int_t           TrigHLT_HLT_Photon22_SC22HE_L1R_v1;
   Int_t           TrigHLT_HLT_Photon25_Cleaned_L1R;
   Int_t           TrigHLT_HLT_Photon30_L1R;
   Int_t           TrigHLT_HLT_Photon30_Cleaned_L1R;
   Int_t           TrigHLT_HLT_Photon30_L1R_8E29;
   Int_t           TrigHLT_HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1;
   Int_t           TrigHLT_HLT_Photon35_Isol_Cleaned_L1R_v1;
   Int_t           TrigHLT_HLT_Photon40_CaloId_Cleaned_L1R_v1;
   Int_t           TrigHLT_HLT_Photon40_Isol_Cleaned_L1R_v1;
   Int_t           TrigHLT_HLT_Photon50_L1R;
   Int_t           TrigHLT_HLT_Photon50_Cleaned_L1R;
   Int_t           TrigHLT_HLT_Photon50_Cleaned_L1R_v1;
   Int_t           TrigHLT_HLT_Photon50_NoHE_L1R;
   Int_t           TrigHLT_HLT_Photon50_NoHE_Cleaned_L1R;
   Int_t           TrigHLT_HLT_Photon70_Cleaned_L1R_v1;
   Int_t           TrigHLT_HLT_Photon70_NoHE_Cleaned_L1R_v1;
   Int_t           TrigHLT_HLT_Photon100_NoHE_Cleaned_L1R_v1;
   Int_t           TrigHLT_HLT_Photon110_NoHE_Cleaned_L1R_v1;
   Int_t           TrigHLT_HLT_DoublePhoton5_L1R;
   Int_t           TrigHLT_HLT_DoublePhoton5_CEP_L1R;
   Int_t           TrigHLT_HLT_DoublePhoton5_CEP_L1R_v3;
   Int_t           TrigHLT_HLT_DoublePhoton5_Jpsi_L1R;
   Int_t           TrigHLT_HLT_DoublePhoton5_Upsilon_L1R;
   Int_t           TrigHLT_HLT_DoublePhoton10_L1R;
   Int_t           TrigHLT_HLT_DoublePhoton15_L1R;
   Int_t           TrigHLT_HLT_DoublePhoton17_L1R;
   Int_t           TrigHLT_HLT_DoublePhoton17_SingleIsol_L1R_v1;
   Int_t           TrigHLT_HLT_DoublePhoton20_L1R;
   Int_t           TrigHLT_HLT_DoublePhoton22_L1R_v1;
   Int_t           TrigHLT_HLT_DoublePhoton33_v1;
   Int_t           TrigHLT_HLT_DoublePhoton33_v2;
   Int_t           TrigHLT_HLT_DoublePhoton33_v3;
   Int_t           TrigHLT_HLT_DoublePhoton33_v5;
   Int_t           TrigHLT_HLT_DoublePhoton33_HEVT_v2;
   Int_t           TrigHLT_HLT_DoublePhoton5_IsoVL_CEP_v1;
   Int_t           TrigHLT_HLT_DoublePhoton5_IsoVL_CEP_v2;
   Int_t           TrigHLT_HLT_Photon125_NoSpikeFilter_v1;
   Int_t           TrigHLT_HLT_Photon125_NoSpikeFilter_v2;
   Int_t           TrigHLT_HLT_Photon125_NoSpikeFilter_v3;
   Int_t           TrigHLT_HLT_Photon20_CaloIdVL_IsoL_v1;
   Int_t           TrigHLT_HLT_Photon20_CaloIdVL_IsoL_v2;
   Int_t           TrigHLT_HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1;
   Int_t           TrigHLT_HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2;
   Int_t           TrigHLT_HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v3;
   Int_t           TrigHLT_HLT_Photon20_EBOnly_NoSpikeFilter_v1;
   Int_t           TrigHLT_HLT_Photon20_NoSpikeFilter_v1;
   Int_t           TrigHLT_HLT_Photon20_R9Id_Photon18_R9Id_v1;
   Int_t           TrigHLT_HLT_Photon20_R9Id_Photon18_R9Id_v2;
   Int_t           TrigHLT_HLT_Photon20_R9Id_Photon18_R9Id_v3;
   Int_t           TrigHLT_HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v1;
   Int_t           TrigHLT_HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v2;
   Int_t           TrigHLT_HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v3;
   Int_t           TrigHLT_HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v1;
   Int_t           TrigHLT_HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v2;
   Int_t           TrigHLT_HLT_Photon26_CaloIdL_IsoVL_Photon18_v1;
   Int_t           TrigHLT_HLT_Photon26_CaloIdL_IsoVL_Photon18_v2;
   Int_t           TrigHLT_HLT_Photon26_CaloIdL_IsoVL_Photon18_v3;
   Int_t           TrigHLT_HLT_Photon26_IsoVL_Photon18_IsoVL_v1;
   Int_t           TrigHLT_HLT_Photon26_IsoVL_Photon18_IsoVL_v2;
   Int_t           TrigHLT_HLT_Photon26_IsoVL_Photon18_IsoVL_v3;
   Int_t           TrigHLT_HLT_Photon26_IsoVL_Photon18_v1;
   Int_t           TrigHLT_HLT_Photon26_IsoVL_Photon18_v2;
   Int_t           TrigHLT_HLT_Photon26_IsoVL_Photon18_v3;
   Int_t           TrigHLT_HLT_Photon26_Photon18_v1;
   Int_t           TrigHLT_HLT_Photon26_Photon18_v2;
   Int_t           TrigHLT_HLT_Photon26_Photon18_v3;
   Int_t           TrigHLT_HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v1;
   Int_t           TrigHLT_HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v2;
   Int_t           TrigHLT_HLT_Photon30_CaloIdVL_IsoL_v1;
   Int_t           TrigHLT_HLT_Photon30_CaloIdVL_IsoL_v2;
   Int_t           TrigHLT_HLT_Photon30_CaloIdVL_IsoL_v3;
   Int_t           TrigHLT_HLT_Photon30_CaloIdVL_v1;
   Int_t           TrigHLT_HLT_Photon30_CaloIdVL_v2;
   Int_t           TrigHLT_HLT_Photon30_CaloIdVL_v3;
   Int_t           TrigHLT_HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1;
   Int_t           TrigHLT_HLT_Photon32_CaloIdL_Photon26_CaloIdL_v2;
   Int_t           TrigHLT_HLT_Photon32_CaloIdL_Photon26_CaloIdL_v3;
   Int_t           TrigHLT_HLT_Photon36_CaloIdL_Photon22_CaloIdL_v1;
   Int_t           TrigHLT_HLT_Photon36_CaloIdL_Photon22_CaloIdL_v2;
   Int_t           TrigHLT_HLT_Photon50_CaloIdVL_IsoL_v1;
   Int_t           TrigHLT_HLT_Photon50_CaloIdVL_IsoL_v2;
   Int_t           TrigHLT_HLT_Photon75_CaloIdVL_IsoL_v1;
   Int_t           TrigHLT_HLT_Photon75_CaloIdVL_IsoL_v2;
   Int_t           TrigHLT_HLT_Photon75_CaloIdVL_IsoL_v3;
   Int_t           TrigHLT_HLT_Photon75_CaloIdVL_v1;
   Int_t           TrigHLT_HLT_Photon75_CaloIdVL_v2;
   Int_t           TrigHLT_HLT_Photon75_CaloIdVL_v3;
   Int_t           TrigHLT_HLT_DoublePhoton40_MR150_v3;
   Int_t           TrigHLT_HLT_DoublePhoton40_R014_MR150_v3;
   Int_t           TrigHLT_HLT_DoublePhoton50_v2;
   Int_t           TrigHLT_HLT_DoublePhoton5_IsoVL_CEP_v4;
   Int_t           TrigHLT_HLT_DoublePhoton60_v2;
   Int_t           TrigHLT_HLT_Mu15_DoublePhoton15_CaloIdL_v6;
   Int_t           TrigHLT_HLT_Mu15_Photon20_CaloIdL_v6;
   Int_t           TrigHLT_HLT_Mu8_Photon20_CaloIdVT_IsoT_v5;
   Int_t           TrigHLT_HLT_Photon125_v2;
   Int_t           TrigHLT_HLT_Photon200_NoHE_v2;
   Int_t           TrigHLT_HLT_Photon20_CaloIdVL_IsoL_v4;
   Int_t           TrigHLT_HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v5;
   Int_t           TrigHLT_HLT_Photon20_R9Id_Photon18_R9Id_v5;
   Int_t           TrigHLT_HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v5;
   Int_t           TrigHLT_HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v4;
   Int_t           TrigHLT_HLT_Photon26_CaloIdL_IsoVL_Photon18_v5;
   Int_t           TrigHLT_HLT_Photon26_IsoVL_Photon18_IsoVL_v5;
   Int_t           TrigHLT_HLT_Photon26_IsoVL_Photon18_v5;
   Int_t           TrigHLT_HLT_Photon26_Photon18_v5;
   Int_t           TrigHLT_HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v4;
   Int_t           TrigHLT_HLT_Photon26_R9Id_Photon18_R9Id_v2;
   Int_t           TrigHLT_HLT_Photon30_CaloIdVL_IsoL_v5;
   Int_t           TrigHLT_HLT_Photon30_CaloIdVL_v5;
   Int_t           TrigHLT_HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v1;
   Int_t           TrigHLT_HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_v1;
   Int_t           TrigHLT_HLT_Photon36_CaloIdL_IsoVL_Photon22_v2;
   Int_t           TrigHLT_HLT_Photon36_CaloIdL_Photon22_CaloIdL_v4;
   Int_t           TrigHLT_HLT_Photon36_CaloId_IsoVL_Photon22_R9Id_v1;
   Int_t           TrigHLT_HLT_Photon36_IsoVL_Photon22_v2;
   Int_t           TrigHLT_HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v1;
   Int_t           TrigHLT_HLT_Photon36_R9Id_Photon22_R9Id_v1;
   Int_t           TrigHLT_HLT_Photon40_CaloIdL_Photon28_CaloIdL_v2;
   Int_t           TrigHLT_HLT_Photon40_R005_MR150_v3;
   Int_t           TrigHLT_HLT_Photon40_R014_MR450_v3;
   Int_t           TrigHLT_HLT_Photon40_R020_MR300_v3;
   Int_t           TrigHLT_HLT_Photon40_R025_MR200_v3;
   Int_t           TrigHLT_HLT_Photon40_R038_MR150_v3;
   Int_t           TrigHLT_HLT_Photon50_CaloIdVL_IsoL_v4;
   Int_t           TrigHLT_HLT_Photon50_CaloIdVL_v2;
   Int_t           TrigHLT_HLT_Photon70_CaloIdL_HT300_v6;
   Int_t           TrigHLT_HLT_Photon70_CaloIdL_HT350_v5;
   Int_t           TrigHLT_HLT_Photon70_CaloIdL_MHT50_v6;
   Int_t           TrigHLT_HLT_Photon70_CaloIdL_MHT70_v5;
   Int_t           TrigHLT_HLT_Photon75_CaloIdVL_IsoL_v5;
   Int_t           TrigHLT_HLT_Photon75_CaloIdVL_v5;
   Int_t           TrigHLT_HLT_Photon90_CaloIdVL_IsoL_v2;
   Int_t           TrigHLT_HLT_Photon90_CaloIdVL_v2;
   Int_t           TrigHLT_HLT_DoublePhoton33_HEVT_v1;
   Int_t           TrigHLT_HLT_DoublePhoton33_v4;
   Int_t           TrigHLT_HLT_DoublePhoton50_v1;
   Int_t           TrigHLT_HLT_DoublePhoton5_IsoVL_CEP_v3;
   Int_t           TrigHLT_HLT_DoublePhoton60_v1;
   Int_t           TrigHLT_HLT_Photon125_v1;
   Int_t           TrigHLT_HLT_Photon200_NoHE_v1;
   Int_t           TrigHLT_HLT_Photon20_CaloIdVL_IsoL_v3;
   Int_t           TrigHLT_HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v4;
   Int_t           TrigHLT_HLT_Photon20_R9Id_Photon18_R9Id_v4;
   Int_t           TrigHLT_HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v4;
   Int_t           TrigHLT_HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v3;
   Int_t           TrigHLT_HLT_Photon26_CaloIdL_IsoVL_Photon18_v4;
   Int_t           TrigHLT_HLT_Photon26_IsoVL_Photon18_IsoVL_v4;
   Int_t           TrigHLT_HLT_Photon26_IsoVL_Photon18_v4;
   Int_t           TrigHLT_HLT_Photon26_Photon18_v4;
   Int_t           TrigHLT_HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v3;
   Int_t           TrigHLT_HLT_Photon26_R9Id_Photon18_R9Id_v1;
   Int_t           TrigHLT_HLT_Photon30_CaloIdVL_IsoL_v4;
   Int_t           TrigHLT_HLT_Photon30_CaloIdVL_v4;
   Int_t           TrigHLT_HLT_Photon32_CaloIdL_Photon26_CaloIdL_v4;
   Int_t           TrigHLT_HLT_Photon36_CaloIdL_IsoVL_Photon22_v1;
   Int_t           TrigHLT_HLT_Photon36_CaloIdL_Photon22_CaloIdL_v3;
   Int_t           TrigHLT_HLT_Photon36_IsoVL_Photon22_v1;
   Int_t           TrigHLT_HLT_Photon40_CaloIdL_Photon28_CaloIdL_v1;
   Int_t           TrigHLT_HLT_Photon50_CaloIdVL_IsoL_v3;
   Int_t           TrigHLT_HLT_Photon50_CaloIdVL_v1;
   Int_t           TrigHLT_HLT_Photon75_CaloIdVL_IsoL_v4;
   Int_t           TrigHLT_HLT_Photon75_CaloIdVL_v4;
   Int_t           TrigHLT_HLT_Photon90_CaloIdVL_IsoL_v1;
   Int_t           TrigHLT_HLT_Photon90_CaloIdVL_v1;

   Int_t           nTightPhotons;
   Int_t           nFakeablePhotons;

   Int_t           GenPhoton1_status;
   Int_t           GenPhoton1_PdgId;
   Int_t           GenPhoton1_MotherPdgId;
   Int_t           GenPhoton1_GrandmotherPdgId;
   Double_t        GenPhoton1_pt;
   Double_t        GenPhoton1_eta;
   Double_t        GenPhoton1_phi;
   Int_t           GenPhoton2_status;
   Int_t           GenPhoton2_PdgId;
   Int_t           GenPhoton2_MotherPdgId;
   Int_t           GenPhoton2_GrandmotherPdgId;
   Double_t        GenPhoton2_pt;
   Double_t        GenPhoton2_eta;
   Double_t        GenPhoton2_phi;

   Double_t        Photon1_pt;
   Double_t        Photon1_eta;
   Double_t        Photon1_phi;
   Double_t        Photon1_detEta;
   Double_t        Photon1_detPhi;
   Double_t        Photon1_r9;
   Double_t        Photon1_sigmaIetaIeta;
   Double_t        Photon1_sigmaEtaEta;
   Double_t        Photon1_maxEnergyXtal;
   Double_t        Photon1_e1x5;
   Double_t        Photon1_e2x5;
   Double_t        Photon1_e3x3;
   Double_t        Photon1_e5x5;
   Double_t        Photon1_r1x5;
   Double_t        Photon1_r2x5;
   Double_t        Photon1_swisscross;
   Double_t        Photon1_eMax;
   Double_t        Photon1_eLeft;
   Double_t        Photon1_eRight;
   Double_t        Photon1_eTop;
   Double_t        Photon1_eBottom;
   Double_t        Photon1_eSecond;
   Double_t        Photon1_e2x2;
   Double_t        Photon1_e4x4;
   Double_t        Photon1_e2e9;
   Double_t        Photon1_maxRecHitTime;
   Double_t        Photon1_hadOverEm;
   Double_t        Photon1_hadDepth1OverEm;
   Double_t        Photon1_hadDepth2OverEm;
   Double_t        Photon1_hcalIso04;
   Double_t        Photon1_hcalIso03;
   Double_t        Photon1_ecalIso04;
   Double_t        Photon1_ecalIso03;
   Double_t        Photon1_trkIsoSumPtHollow04;
   Double_t        Photon1_trkIsoSumPtSolid04;
   Double_t        Photon1_trkIsoSumPtHollow03;
   Double_t        Photon1_trkIsoSumPtSolid03;
   Double_t        Photon1_esRatio;
   Double_t        Photon1_scRawEnergy;
   Double_t        Photon1_scPreshowerEnergy;
   Double_t        Photon1_scPhiWidth;
   Double_t        Photon1_scEtaWidth;
   Int_t           Photon1_scNumBasicClusters;
   Int_t           Photon1_trkIsoNtrksHollow04;
   Int_t           Photon1_trkIsoNtrksSolid04;
   Int_t           Photon1_trkIsoNtrksHollow03;
   Int_t           Photon1_trkIsoNtrksSolid03;
   Int_t           Photon1_severityLevel;
   Int_t           Photon1_recHitFlag;
   Int_t           Photon1_detId;
   Int_t           Photon1_iEtaY;
   Int_t           Photon1_iPhiX;
   Bool_t          Photon1_isEB;
   Bool_t          Photon1_isEE;
   Bool_t          Photon1_isEBEtaGap;
   Bool_t          Photon1_isEBPhiGap;
   Bool_t          Photon1_isEERingGap;
   Bool_t          Photon1_isEEDeeGap;
   Bool_t          Photon1_isEBEEGap;
   Bool_t          Photon1_hasPixelSeed;
   Bool_t          Photon1_isFakeable;
   Double_t        Photon2_pt;
   Double_t        Photon2_eta;
   Double_t        Photon2_phi;
   Double_t        Photon2_detEta;
   Double_t        Photon2_detPhi;
   Double_t        Photon2_r9;
   Double_t        Photon2_sigmaIetaIeta;
   Double_t        Photon2_sigmaEtaEta;
   Double_t        Photon2_maxEnergyXtal;
   Double_t        Photon2_e1x5;
   Double_t        Photon2_e2x5;
   Double_t        Photon2_e3x3;
   Double_t        Photon2_e5x5;
   Double_t        Photon2_r1x5;
   Double_t        Photon2_r2x5;
   Double_t        Photon2_swisscross;
   Double_t        Photon2_eMax;
   Double_t        Photon2_eLeft;
   Double_t        Photon2_eRight;
   Double_t        Photon2_eTop;
   Double_t        Photon2_eBottom;
   Double_t        Photon2_eSecond;
   Double_t        Photon2_e2x2;
   Double_t        Photon2_e4x4;
   Double_t        Photon2_e2e9;
   Double_t        Photon2_maxRecHitTime;
   Double_t        Photon2_hadOverEm;
   Double_t        Photon2_hadDepth1OverEm;
   Double_t        Photon2_hadDepth2OverEm;
   Double_t        Photon2_hcalIso04;
   Double_t        Photon2_hcalIso03;
   Double_t        Photon2_ecalIso04;
   Double_t        Photon2_ecalIso03;
   Double_t        Photon2_trkIsoSumPtHollow04;
   Double_t        Photon2_trkIsoSumPtSolid04;
   Double_t        Photon2_trkIsoSumPtHollow03;
   Double_t        Photon2_trkIsoSumPtSolid03;
   Double_t        Photon2_esRatio;
   Double_t        Photon2_scRawEnergy;
   Double_t        Photon2_scPreshowerEnergy;
   Double_t        Photon2_scPhiWidth;
   Double_t        Photon2_scEtaWidth;
   Int_t           Photon2_scNumBasicClusters;
   Int_t           Photon2_trkIsoNtrksHollow04;
   Int_t           Photon2_trkIsoNtrksSolid04;
   Int_t           Photon2_trkIsoNtrksHollow03;
   Int_t           Photon2_trkIsoNtrksSolid03;
   Int_t           Photon2_severityLevel;
   Int_t           Photon2_recHitFlag;
   Int_t           Photon2_detId;
   Int_t           Photon2_iEtaY;
   Int_t           Photon2_iPhiX;
   Bool_t          Photon2_isEB;
   Bool_t          Photon2_isEE;
   Bool_t          Photon2_isEBEtaGap;
   Bool_t          Photon2_isEBPhiGap;
   Bool_t          Photon2_isEERingGap;
   Bool_t          Photon2_isEEDeeGap;
   Bool_t          Photon2_isEBEEGap;
   Bool_t          Photon2_hasPixelSeed;
   Bool_t          Photon2_isFakeable;

   Int_t           MCMatchPhoton1_Status3_status;
   Int_t           MCMatchPhoton1_Status3_PdgId;
   Int_t           MCMatchPhoton1_Status3_MotherPdgId;
   Int_t           MCMatchPhoton1_Status3_GrandmotherPdgId;
   Double_t        MCMatchPhoton1_Status3_pt;
   Double_t        MCMatchPhoton1_Status3_eta;
   Double_t        MCMatchPhoton1_Status3_phi;
   Int_t           MCMatchPhoton2_Status3_status;
   Int_t           MCMatchPhoton2_Status3_PdgId;
   Int_t           MCMatchPhoton2_Status3_MotherPdgId;
   Int_t           MCMatchPhoton2_Status3_GrandmotherPdgId;
   Double_t        MCMatchPhoton2_Status3_pt;
   Double_t        MCMatchPhoton2_Status3_eta;
   Double_t        MCMatchPhoton2_Status3_phi;
   Int_t           MCMatchPhoton1_Status1_status;
   Int_t           MCMatchPhoton1_Status1_PdgId;
   Int_t           MCMatchPhoton1_Status1_MotherPdgId;
   Int_t           MCMatchPhoton1_Status1_GrandmotherPdgId;
   Double_t        MCMatchPhoton1_Status1_pt;
   Double_t        MCMatchPhoton1_Status1_eta;
   Double_t        MCMatchPhoton1_Status1_phi;
   Int_t           MCMatchPhoton2_Status1_status;
   Int_t           MCMatchPhoton2_Status1_PdgId;
   Int_t           MCMatchPhoton2_Status1_MotherPdgId;
   Int_t           MCMatchPhoton2_Status1_GrandmotherPdgId;
   Double_t        MCMatchPhoton2_Status1_pt;
   Double_t        MCMatchPhoton2_Status1_eta;
   Double_t        MCMatchPhoton2_Status1_phi;

   Double_t        Diphoton_Minv;
   Double_t        Diphoton_qt;
   Double_t        Diphoton_deltaPhi;
   Double_t        Diphoton_deltaEta;
   Double_t        Diphoton_deltaR;
   Double_t        Diphoton_cosThetaStar;
   Double_t        DiphotonVtx2_Minv;
   Double_t        DiphotonVtx2_qt;
   Double_t        DiphotonVtx2_deltaPhi;
   Double_t        DiphotonVtx2_deltaEta;
   Double_t        DiphotonVtx2_deltaR;
   Double_t        DiphotonVtx2_cosThetaStar;
   Double_t        DiphotonVtx3_Minv;
   Double_t        DiphotonVtx3_qt;
   Double_t        DiphotonVtx3_deltaPhi;
   Double_t        DiphotonVtx3_deltaEta;
   Double_t        DiphotonVtx3_deltaR;
   Double_t        DiphotonVtx3_cosThetaStar;

   Double_t        DiphotonGen_Minv;
   Double_t        DiphotonGen_qt;
   Double_t        DiphotonGen_deltaPhi;
   Double_t        DiphotonGen_deltaEta;
   Double_t        DiphotonGen_deltaR;
   Double_t        DiphotonGen_cosThetaStar;

   // List of branches
   TBranch        *b_Event;   //!
   TBranch        *b_GenEvent;   //! 
   TBranch        *b_Vtx;   //!
   TBranch        *b_Vtx2;   //!
   TBranch        *b_Vtx3;   //!
   TBranch        *b_VtxGEN;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_pu_n;   //!
   TBranch        *b_BeamSpot;   //!
   TBranch        *b_L1trg;   //!
   TBranch        *b_TrigHLT;   //!
   TBranch        *b_nTightPhotons;   //!
   TBranch        *b_nFakeablePhotons;   //!
   TBranch        *b_GenPhoton1;   //! 
   TBranch        *b_GenPhoton2;   //! 
   TBranch        *b_Photon1;   //!
   TBranch        *b_Photon2;   //!
   TBranch        *b_MCMatchPhoton1_Status3;   //! 
   TBranch        *b_MCMatchPhoton2_Status3;   //! 
   TBranch        *b_MCMatchPhoton1_Status1;   //! 
   TBranch        *b_MCMatchPhoton2_Status1;   //! 
   TBranch        *b_Diphoton;   //!
   TBranch        *b_DiphotonVtx2;   //!
   TBranch        *b_DiphotonVtx3;   //!
   TBranch        *b_DiphotonGen;   //! 


   // diphoton analysis
   Float_t _cutPhoton1Pt;
   Float_t _cutPhoton2Pt;
   Float_t _cutEta;

   Bool_t  _filterGen;
   TString _categoryEBEE;

   TString _fakeRateFile;

   TString _puFile;
   TH1F*   _puHistMC;
   TString _kFactorFile;
   TFile*  fInKfactor;

   Bool_t  _reweightPU;
   Bool_t  _Kfactor;
   Float_t _effScaleFactor;

   TString _fakeStatus;

   Double_t  _weight;

   TFile*  _outF;
   TTree*  _outT0;
   
   struct outputS {
     Double_t minv; 
     Double_t r9_photon1;
     Double_t r9_photon2;
     Double_t eta_photon1;
     Double_t eta_photon2;
     Double_t weight;
   };
   outputS outVar0; 


   TH1* h_mcMatch1;
   TH1* h_mcMatch2;
      
   TH1* h_Nvtx;
   TH1* h_TrigHLT;

   TH1* h_Diphoton_Minv;
   //   TH1* h_Diphoton_Minv_add;
   TH1* h_Diphoton_Minv_log;
   TH1* h_Diphoton_Minv_high;
   TH1* h_Diphoton_Minv_low;
   TH1* h_Diphoton_Minv_low_bin1GeV;
   TH1* h_Diphoton_qt;
   TH1* h_Diphoton_deltaPhi;
   TH1* h_Diphoton_deltaEta;
   TH1* h_Diphoton_deltaR;
   TH1* h_Diphoton_cosThetaStar;
   TH2* h_Diphoton_Minv_v_Photon1_pt; 
   TH2* h_Diphoton_Minv_v_Photon2_pt; 
   TH1* h_Diphoton_Minv_120to200;
   TH1* h_Diphoton_Minv_200to500;
   TH1* h_Diphoton_Minv_500to800;
   TH1* h_Diphoton_Minv_800toInf;
   TH1* h_Diphoton_Minv_Yousi5;
   TH1* h_Diphoton_Minv_Yousi10;
   TH1* h_Diphoton_Minv_Yousi40;
   TH1* h_Diphoton_Minv_Yousi100;

   TH1* h_Photon1_pt; 
   TH1* h_Photon1_pt_log; 
   TH1* h_Photon1_pt_zoom; 
   TH1* h_Photon1_eta; 
   TH1* h_Photon1_phi; 
   TH2* h_Photon1_occupancy; 
   
   TH1* h_Photon1_r9; 
   TH1* h_Photon1_sigmaIetaIeta; 
   TH1* h_Photon1_sigmaEtaEta; 
   
   TH1* h_Photon1_e2x2e4x4; 
   TH1* h_Photon1_e2e9; 
   TH1* h_Photon1_swisscross; 
   TH1* h_Photon1_severityLevel; 
   TH1* h_Photon1_recHitFlag; 
   TH1* h_Photon1_maxRecHitTime; 
   TH1* h_Photon1_maxRecHitTime_wide; 
   TH2* h_Photon1_e2e9_v_maxRecHitTime; 
   
   TH1* h_Photon1_hadOverEm; 
   TH1* h_Photon1_hadDepth1OverEm; 
   TH1* h_Photon1_hadDepth2OverEm; 
   
   TH1* h_Photon1_hcalIso04; 
   TH1* h_Photon1_hcalIso03; 
   
   TH1* h_Photon1_ecalIso04; 
   TH1* h_Photon1_ecalIso03; 
   
   TH1* h_Photon1_trkIsoSumPtHollow04; 
   TH1* h_Photon1_trkIsoSumPtSolid04; 
   TH1* h_Photon1_trkIsoNtrksHollow04; 
   TH1* h_Photon1_trkIsoNtrksSolid04; 
   
   TH1* h_Photon1_trkIsoSumPtHollow03; 
   TH1* h_Photon1_trkIsoSumPtSolid03; 
   TH1* h_Photon1_trkIsoNtrksHollow03; 
   TH1* h_Photon1_trkIsoNtrksSolid03; 
   
   TH1* h_Photon1_esRatio; 

   TH1* h_Photon2_pt; 
   TH1* h_Photon2_pt_log; 
   TH1* h_Photon2_pt_zoom; 
   TH1* h_Photon2_eta; 
   TH1* h_Photon2_phi; 
   TH2* h_Photon2_occupancy; 

   TH1* h_Photon2_r9; 
   TH1* h_Photon2_sigmaIetaIeta; 
   TH1* h_Photon2_sigmaEtaEta; 
   
   TH1* h_Photon2_e2x2e4x4; 
   TH1* h_Photon2_e2e9; 
   TH1* h_Photon2_swisscross; 
   TH1* h_Photon2_severityLevel; 
   TH1* h_Photon2_recHitFlag; 
   TH1* h_Photon2_maxRecHitTime; 
   TH1* h_Photon2_maxRecHitTime_wide; 
   TH2* h_Photon2_e2e9_v_maxRecHitTime; 
   
   TH1* h_Photon2_hadOverEm; 
   TH1* h_Photon2_hadDepth1OverEm; 
   TH1* h_Photon2_hadDepth2OverEm; 
   
   TH1* h_Photon2_hcalIso04; 
   TH1* h_Photon2_hcalIso03; 
   
   TH1* h_Photon2_ecalIso04; 
   TH1* h_Photon2_ecalIso03; 
   
   TH1* h_Photon2_trkIsoSumPtHollow04; 
   TH1* h_Photon2_trkIsoSumPtSolid04; 
   TH1* h_Photon2_trkIsoNtrksHollow04; 
   TH1* h_Photon2_trkIsoNtrksSolid04; 
   
   TH1* h_Photon2_trkIsoSumPtHollow03; 
   TH1* h_Photon2_trkIsoSumPtSolid03; 
   TH1* h_Photon2_trkIsoNtrksHollow03; 
   TH1* h_Photon2_trkIsoNtrksSolid03; 
   
   TH1* h_Photon2_esRatio; 

   // tight-tight ie signal sample

   TH1* h_FakeRate_tt_pt1;
   TH1* h_FakeRate_tt_pt1_zoom;
   TH1* h_FakeRate_tt_eta1;
   TH1* h_FakeRate_tt_phi1;
   TH1* h_FakeRate_tt_pt2;
   TH1* h_FakeRate_tt_pt2_zoom;
   TH1* h_FakeRate_tt_eta2;
   TH1* h_FakeRate_tt_phi2;
   TH1* h_FakeRate_tt_minv;
   TH1* h_FakeRate_tt_minv_high;
   TH1* h_FakeRate_tt_qt;
   TH1* h_FakeRate_tt_deltaPhi;
   TH1* h_FakeRate_tt_deltaEta;
   TH1* h_FakeRate_tt_deltaR;
   TH1* h_FakeRate_tt_cosThetaStar;
   
   // tight-fake
   TH1* h_FakeRate_tf_pt1;
   TH1* h_FakeRate_tf_pt1_zoom;
   TH1* h_FakeRate_tf_eta1;
   TH1* h_FakeRate_tf_phi1;
   TH1* h_FakeRate_tf_pt2;
   TH1* h_FakeRate_tf_pt2_zoom;
   TH1* h_FakeRate_tf_eta2;
   TH1* h_FakeRate_tf_phi2;
   TH1* h_FakeRate_tf_minv;
   TH1* h_FakeRate_tf_minv_high;
   TH1* h_FakeRate_tf_qt;
   TH1* h_FakeRate_tf_deltaPhi;
   TH1* h_FakeRate_tf_deltaEta;
   TH1* h_FakeRate_tf_deltaR;
   TH1* h_FakeRate_tf_cosThetaStar;
   
   TH1* h_FakeRate_tf_pt1_noweight;
   TH1* h_FakeRate_tf_pt2_noweight;

   TH1* h_FakeRate_ft_pt1_noweight;
   TH1* h_FakeRate_ft_pt2_noweight;

   // fake-fake
   TH1* h_FakeRate_ff_pt1;
   TH1* h_FakeRate_ff_pt1_zoom;
   TH1* h_FakeRate_ff_eta1;
   TH1* h_FakeRate_ff_phi1;
   TH1* h_FakeRate_ff_pt2;
   TH1* h_FakeRate_ff_pt2_zoom;
   TH1* h_FakeRate_ff_eta2;
   TH1* h_FakeRate_ff_phi2;
   TH1* h_FakeRate_ff_minv;   
   TH1* h_FakeRate_ff_minv_high;   
   TH1* h_FakeRate_ff_qt;
   TH1* h_FakeRate_ff_deltaPhi;
   TH1* h_FakeRate_ff_deltaEta;
   TH1* h_FakeRate_ff_deltaR;
   TH1* h_FakeRate_ff_cosThetaStar;

   TH1* h_FakeRate_ff_pt1_noweight;
   TH1* h_FakeRate_ff_pt2_noweight;

   TH2* h_trkIsoVRho1;
   TH2* h_ecalIsoVRho1;
   TH2* h_hcalIsoVRho1;
   TH2* h_trkIsoVRho2;
   TH2* h_ecalIsoVRho2;
   TH2* h_hcalIsoVRho2;

   TH1* h_FakeRate_tt_minv_120to200;
   TH1* h_FakeRate_tt_minv_200to500;
   TH1* h_FakeRate_tt_minv_500to800;
   TH1* h_FakeRate_tt_minv_800toInf;

   TH1* h_FakeRate_tf_minv_120to200;
   TH1* h_FakeRate_tf_minv_200to500;
   TH1* h_FakeRate_tf_minv_500to800;
   TH1* h_FakeRate_tf_minv_800toInf;

   TH1* h_FakeRate_ff_minv_120to200;
   TH1* h_FakeRate_ff_minv_200to500;
   TH1* h_FakeRate_ff_minv_500to800;
   TH1* h_FakeRate_ff_minv_800toInf;
   
   fTree(TTree *tree=0);
   virtual ~fTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   double kfactorF(float x) {

     //     return -6.08061 + 9.90359/TMath::Power(x,3.88349e-02); // old one
     //     return 2.91986557985310036-0.284606795636766441/TMath::Power(x, -2.32654814731609899e-01) ;     
     //     return 1.6;
     //     return 2.91986557985310036-0.284606795636766441/TMath::Power(x, -2.32654814731609899e-01) ;

     //     TH1D* hKfactor = (TH1D*)fInKfactor->Get("k-factorTmp") ;
     //     TGraph* g = new TGraph(hKfactor);
     //     return g->Eval(x);

     if (x<320) {
       return 1.73405 ;
     } else {
       return 1.07559 + 1.22200*TMath::Exp((-1.85612e-3)*x);
     }
   }

   vector<double> generate_flat10_weights(TH1D* data_npu_estimated, TH1F* mc_npu ){

     cout << " generating weights " << endl;
     // see SimGeneral/MixingModule/python/mix_E7TeV_FlatDist10_2011EarlyData_inTimeOnly_cfi.py; copy and paste from there:
     //     const double npu_probs[25] = {0.0698146584, 0.0698146584, 0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584 /* <-- 10*/,				   0.0630151648,0.0526654164,0.0402754482,0.0292988928,0.0194384503,0.0122016783,0.007207042,0.004003637,0.0020278322,
     //				   0.0010739954,0.0004595759,0.0002229748,0.0001028162,4.58337152809607E-05 /* <-- 24 */};
     
     //     const double npu_probs[25] = {  0.078138816,  0.076506664,  0.068656428,  0.075431185,  0.073189343,  0.073480934,  0.071769257,  0.067796803,  0.073704361,  0.073814181,  0.066293405,  0.057038231,  0.048566944,  0.034847019,   0.021276676,   0.01644649, 0.010870286,  0.007108004,  0.0035256,  0.001047077,  7.57379E-06,  3.97624E-05,  6.43772E-05,  0.000200705, 0.000179878 };


     // the latest
     //     const double npu_probs[25] = { 0.07813870, 0.07650480, 0.06865510, 0.07542560,0.07319340,0.07348120,0.07177530,0.06779940,0.07370080,0.07381820,0.06629420,0.05704160,0.04856530,0.03484440,0.02127510,0.01644530,0.01086950,0.00710748,0.00352534,0.00104700,0.00000757,0.00003976,0.00006437,0.00020069,0.00017986  };
     
     const double npu_probs[25] = {     0.106825,  0.0775458,  0.0778654,  0.0786776,  0.0779586,  0.0757749,  0.0692773,  0.0721932,  0.0653893,  0.0588784,  0.0558292,  0.0474808,  0.0405038,  0.0307707,  0.0219296,  0.0164705,  0.0101859,  0.00734981,  0.00410098,  0.0022369,  0.00135812,  0.000599169,  0.000346187,  0.000226353, 7.98892e-05 };

     double npu_mc[25];

     vector<double> result(25);
     //     cout << "mc_npu " << mc_npu << endl;
     double s = 0.0;
     double smc = 0.0;
     for(int npu=0; npu<25; ++npu){
       double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));
       if (mc_npu!=0) {
	 //	 cout << "we have MC histogram!" << endl;
	 cout << npu << endl;
	 cout << mc_npu->GetBinContent(mc_npu->GetXaxis()->FindBin(npu)) << endl;
       	 npu_mc[npu] = mc_npu->GetBinContent(mc_npu->GetXaxis()->FindBin(npu));                              
       } else {
	 npu_mc[npu] = npu_probs[npu];
       }
       //       cout << npu_estimated << " " << npu_mc << endl;
       result[npu] = npu_estimated / npu_mc[npu];
       s += npu_estimated;
       smc += npu_mc[npu];
     }
     //     cout << s << " " << smc << endl;
     // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
     for(int npu=0; npu<25; ++npu){
       result[npu] /= s;
       result[npu] *= smc;
     }
     return result;
   }
   
   double deltaPhi(double phi1, double phi2) { 
     double result = phi1 - phi2;
     while (result > M_PI) result -= 2*M_PI;
     while (result <= -M_PI) result += 2*M_PI;
     return result;
   }

};

#endif

#ifdef fTree_cxx

fTree::fTree(TTree *tree)
  :_cutPhoton1Pt(65)
  , _cutPhoton2Pt(65)
  , _cutEta(2.5)
  , _filterGen(kFALSE)
  , _categoryEBEE("ALL")
  , _fakeRateFile("fake_rate_fit_functions.root")
  //  , _puFile("pudist_160404-163869_Cert_JSON.root")
  , _puFile("Pileup_2011_EPS_8_jul.root")
  , _kFactorFile("kFactor_new.root")
  , _reweightPU(kTRUE)
  , _Kfactor(kFALSE)
  , _effScaleFactor(1.005)
  , _fakeStatus("TightTight")
  , _weight(1.0)
  , _outF(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  if (tree == 0) {
    std::cerr << "Give me a tree ! " << std::endl;
    return;
  }
  Init(tree);

}

fTree::~fTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t fTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t fTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void fTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event", &Event_run, &b_Event);
   fChain->SetBranchAddress("GenEvent", &GenEvent_signalProcessId, &b_GenEvent);
   fChain->SetBranchAddress("Vtx", &Vtx_vx, &b_Vtx);
   fChain->SetBranchAddress("Vtx2", &Vtx2_vx, &b_Vtx2);
   fChain->SetBranchAddress("Vtx3", &Vtx3_vx, &b_Vtx3);
   fChain->SetBranchAddress("VtxGEN", &VtxGEN_vx, &b_VtxGEN);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("pu_n", &pu_n, &b_pu_n);
   fChain->SetBranchAddress("BeamSpot", &BeamSpot_x0, &b_BeamSpot);
   fChain->SetBranchAddress("L1trg", &L1trg_L1_Tech0, &b_L1trg);
   fChain->SetBranchAddress("TrigHLT", &TrigHLT_HLT_MinBiasBSC, &b_TrigHLT);
   fChain->SetBranchAddress("nTightPhotons", &nTightPhotons, &b_nTightPhotons);
   fChain->SetBranchAddress("nFakeablePhotons", &nFakeablePhotons, &b_nFakeablePhotons);
   fChain->SetBranchAddress("GenPhoton1", &GenPhoton1_pt, &b_GenPhoton1);
   fChain->SetBranchAddress("GenPhoton2", &GenPhoton2_pt, &b_GenPhoton2);
   fChain->SetBranchAddress("Photon1", &Photon1_pt, &b_Photon1);
   fChain->SetBranchAddress("Photon2", &Photon2_pt, &b_Photon2);
   fChain->SetBranchAddress("MCMatchPhoton1_Status3", &MCMatchPhoton1_Status3_status, &b_MCMatchPhoton1_Status3);
   fChain->SetBranchAddress("MCMatchPhoton2_Status3", &MCMatchPhoton2_Status3_status, &b_MCMatchPhoton2_Status3);
   fChain->SetBranchAddress("MCMatchPhoton1_Status1", &MCMatchPhoton1_Status1_status, &b_MCMatchPhoton1_Status1);
   fChain->SetBranchAddress("MCMatchPhoton2_Status1", &MCMatchPhoton2_Status1_status, &b_MCMatchPhoton2_Status1);
   fChain->SetBranchAddress("Diphoton", &Diphoton_Minv, &b_Diphoton);
   fChain->SetBranchAddress("DiphotonVtx2", &DiphotonVtx2_Minv, &b_DiphotonVtx2);
   fChain->SetBranchAddress("DiphotonVtx3", &DiphotonVtx3_Minv, &b_DiphotonVtx3);
   fChain->SetBranchAddress("DiphotonGen", &DiphotonGen_Minv, &b_DiphotonGen);

   Notify();

   h_mcMatch1 = new TH1F("h_mcMatch1","match1",100,0,0.1);
   h_mcMatch2 = new TH1F("h_mcMatch2","match2",100,0,0.1);
   h_Nvtx = new TH1F("h_Nvtx","Num Vertices",20,0,20);
   h_TrigHLT = new TH1F("h_TrigHLT","HLT Results" , 186, 0., 186);
   //   h_TrigHLT->GetXaxis()->SetBinLabel(1,"HLT_MinBiasBSC");
   //   h_TrigHLT->GetXaxis()->SetBinLabel(2,"HLT_MinBiasBSC_NoBPTX");
   //   h_TrigHLT->GetXaxis()->SetBinLabel(3,"HLT_MinBiasBSC_OR");
   //   h_TrigHLT->GetXaxis()->SetBinLabel(4,"HLT_L1_BscMinBiasOR_BptxPlusORMinus");
   h_TrigHLT->GetXaxis()->SetBinLabel(1,"HLT_L1SingleEG2");
   h_TrigHLT->GetXaxis()->SetBinLabel(2,"HLT_L1SingleEG5");
   h_TrigHLT->GetXaxis()->SetBinLabel(3,"HLT_L1SingleEG8");
   h_TrigHLT->GetXaxis()->SetBinLabel(4,"HLT_L1DoubleEG5");
   h_TrigHLT->GetXaxis()->SetBinLabel(5,"HLT_Photon10_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(6,"HLT_Photon10_Cleaned_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(7,"HLT_Photon15_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(8,"HLT_Photon15_Cleaned_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(9,"HLT_Photon15_LooseEcalIso_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(10,"HLT_Photon15_LooseEcalIso_Cleaned_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(11,"HLT_Photon15_TrackIso_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(12,"HLT_Photon15_TrackIso_Cleaned_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(13,"HLT_Photon17_Isol_SC17HE_L1R_v1");
   h_TrigHLT->GetXaxis()->SetBinLabel(14,"HLT_Photon17_SC17HE_L1R_v1");
   h_TrigHLT->GetXaxis()->SetBinLabel(15,"HLT_Photon20_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(16,"HLT_Photon20_Cleaned_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(17,"HLT_Photon20_NoHE_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(18,"HLT_Photon22_SC22HE_L1R_v1");
   h_TrigHLT->GetXaxis()->SetBinLabel(19,"HLT_Photon25_Cleaned_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(20,"HLT_Photon30_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(21,"HLT_Photon30_Cleaned_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(22,"HLT_Photon30_L1R_8E29");
   h_TrigHLT->GetXaxis()->SetBinLabel(23,"HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1");
   h_TrigHLT->GetXaxis()->SetBinLabel(24,"HLT_Photon35_Isol_Cleaned_L1R_v1");
   h_TrigHLT->GetXaxis()->SetBinLabel(25,"HLT_Photon40_CaloId_Cleaned_L1R_v1");
   h_TrigHLT->GetXaxis()->SetBinLabel(26,"HLT_Photon40_Isol_Cleaned_L1R_v1");
   h_TrigHLT->GetXaxis()->SetBinLabel(27,"HLT_Photon50_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(28,"HLT_Photon50_Cleaned_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(29,"HLT_Photon50_Cleaned_L1R_v1");
   h_TrigHLT->GetXaxis()->SetBinLabel(30,"HLT_Photon50_NoHE_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(31,"HLT_Photon50_NoHE_Cleaned_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(32,"HLT_Photon70_Cleaned_L1R_v1");
   h_TrigHLT->GetXaxis()->SetBinLabel(33,"HLT_Photon70_NoHE_Cleaned_L1R_v1");
   h_TrigHLT->GetXaxis()->SetBinLabel(34,"HLT_Photon100_NoHE_Cleaned_L1R_v1");
   h_TrigHLT->GetXaxis()->SetBinLabel(35,"HLT_Photon110_NoHE_Cleaned_L1R_v1");
   h_TrigHLT->GetXaxis()->SetBinLabel(36,"HLT_DoublePhoton5_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(37,"HLT_DoublePhoton5_CEP_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(38,"HLT_DoublePhoton5_CEP_L1R_v3");
   h_TrigHLT->GetXaxis()->SetBinLabel(39,"HLT_DoublePhoton5_Jpsi_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(40,"HLT_DoublePhoton5_Upsilon_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(41,"HLT_DoublePhoton10_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(42,"HLT_DoublePhoton15_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(43,"HLT_DoublePhoton17_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(44,"HLT_DoublePhoton17_SingleIsol_L1R_v1");
   h_TrigHLT->GetXaxis()->SetBinLabel(45,"HLT_DoublePhoton20_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(46,"HLT_DoublePhoton22_L1R_v1");

   //   h_Diphoton_Minv = new TH1F("h_Diphoton_Minv","Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",125,0,2500); // 20 GeV bins
   //   h_Diphoton_Minv = new TH1F("h_Diphoton_Minv","Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",43,140,1000); // 20 GeV bins
   h_Diphoton_Minv = new TH1F("h_Diphoton_Minv","Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",1860,140,2000); // 20 GeV bins
   //
   //   h_Diphoton_Minv = new TH1F("h_Diphoton_Minv","Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",44,120,1000); // 20 GeV bins
   //   h_Diphoton_Minv = new TH1F("h_Diphoton_Minv","Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",1880,120,2000); // 40 GeV bins
   h_Diphoton_Minv_Yousi5 = new TH1F("h_Diphoton_Minv_Yousi5","Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",376,120,2000); 
   h_Diphoton_Minv_Yousi10 = new TH1F("h_Diphoton_Minv_Yousi10","Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",188,120,2000); 
   h_Diphoton_Minv_Yousi40 = new TH1F("h_Diphoton_Minv_Yousi40","Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",47,120,2000); 
   h_Diphoton_Minv_Yousi100 = new TH1F("h_Diphoton_Minv_Yousi100","Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",20,120,2000); 
   h_Diphoton_Minv_120to200 = new TH1F("h_Diphoton_Minv_120to200","Diphoton Invariant Mass (120-200 bin);M_{#gamma#gamma} [GeV/c^2]",1,120,200); 
   h_Diphoton_Minv_200to500 = new TH1F("h_Diphoton_Minv_200to500","Diphoton Invariant Mass (200-500 bin);M_{#gamma#gamma} [GeV/c^2]",1,200,500); 
   h_Diphoton_Minv_500to800 = new TH1F("h_Diphoton_Minv_500to800","Diphoton Invariant Mass (500-800 bin);M_{#gamma#gamma} [GeV/c^2]",1,500,800); 
   h_Diphoton_Minv_800toInf = new TH1F("h_Diphoton_Minv_800toInf","Diphoton Invariant Mass (800-10000 bin);M_{#gamma#gamma} [GeV/c^2]",1,800,100000); 

   h_Diphoton_Minv_log = new TH1F("h_Diphoton_Minv_log","Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",43,140,1000);  // 20 GeV bins
   h_Diphoton_Minv_high = new TH1F("h_Diphoton_Minv_high","Diphoton Invariant Mass (High mass events);M_{#gamma#gamma} [GeV/c^2]",35,500,1200); // 20 GeV bins
   h_Diphoton_Minv_low = new TH1F("h_Diphoton_Minv_low","Diphoton Invariant Mass (Low mass events);M_{#gamma#gamma} [GeV/c^2]",190,120,500); // 2 GeV bins
   h_Diphoton_Minv_low_bin1GeV = new TH1F("h_Diphoton_Minv_low_bin1GeV","Diphoton Invariant Mass (Low mass events);M_{#gamma#gamma} [GeV/c^2]",380,120,500); // 1 GeV bins
   h_Diphoton_qt = new TH1F("h_Diphoton_qt","Diphoton qt;#gamma#gamma qt [GeV]",20,0,50);
   h_Diphoton_deltaPhi = new TH1F("h_Diphoton_deltaPhi","Diphoton #Delta#phi;#gamma#gamma #Delta#phi",90,-3.14159,3.14159); 
   h_Diphoton_deltaEta = new TH1F("h_Diphoton_deltaEta","Diphoton #Delta#eta;#gamma#gamma #Delta#eta",40,-6.,6.); 
   h_Diphoton_deltaR = new TH1F("h_Diphoton_deltaR","Diphoton #DeltaR; #gamma#gamma #DeltaR",70,0,7.); 
   h_Diphoton_cosThetaStar = new TH1F("h_Diphoton_cosThetaStar","Diphoton |cos(#theta *)|; #gamma#gamma |cos#theta*|",20,0,1); 
   h_Diphoton_Minv_v_Photon1_pt = new TH2F("h_Diphoton_Minv_v_Photon1_pt","Diphoton Invariant Mass v Photon1 pt;M_{#gamma#gamma} [GeV/c^2];#gamma_{1} p_{T} [GeV]",50,0,500,100,0,500); 
   h_Diphoton_Minv_v_Photon2_pt = new TH2F("h_Diphoton_Minv_v_Photon2_pt","Diphoton Invariant Mass v Photon2 pt;M_{#gamma#gamma} [GeV/c^2];#gamma_{2} p_{T} [GeV]",50,0,500,100,0,500); 

   h_Photon1_pt = new TH1F("h_Photon1_pt","Photon1 pt;#gamma_{1} p_{T} [GeV]",44,60,500); 
   h_Photon1_pt_log = new TH1F("h_Photon1_pt_log","Photon1 pt;#gamma_{1} p_{T} [GeV]",44,60,500); 
   h_Photon1_pt_zoom = new TH1F("h_Photon1_pt_zoom","Photon1 pt;#gamma_{1} p_{T} [GeV]",30,0,150); 
   h_Photon1_eta = new TH1F("h_Photon1_eta","Photon1 #eta;#gamma_{1} #eta",25,-2.5,2.5); 
   h_Photon1_phi = new TH1F ("h_Photon1_phi","Photon1 #phi;#gamma_{1} #phi",20,-3.14159,3.14159); 
   h_Photon1_occupancy = new TH2F ("h_Photon1_occupancy","Photon1 (#eta,#phi);#gamma_{1} #eta; #gamma_{1} #phi",25,-2.5,2.5,20,-3.14159,3.14159); 

   h_Photon1_r9 = new TH1F("h_Photon1_r9","Photon1 R9;#gamma_{1} R9", 50, 0.1, 1.5);
   h_Photon1_sigmaIetaIeta = new TH1F("h_Photon1_sigmaIetaIeta","Photon1 #sigma_{i#etai#eta};#gamma_{1} #sigma_{i#etai#eta}", 25, 0, 0.05);
   h_Photon1_sigmaEtaEta = new TH1F("h_Photon1_sigmaEtaEta","Photon1 #sigma_{#eta#eta};#gamma_{1} #sigma_{#eta#eta}", 25, 0, 0.05);

   h_Photon1_e2x2e4x4 = new TH1F("h_Photon1_e2x2e4x4","Photon1 E2x2/E4x4;#gamma_{1} E_{2x2}/E_{4x4}", 20, 0.4, 1.2);
   h_Photon1_e2e9 = new TH1F("h_Photon1_e2e9","Photon1 E2/E9;#gamma_{1} E_{2}/E_{9}", 20, 0.4, 1.2);
   h_Photon1_swisscross = new TH1F("h_Photon1_swisscross","Photon1 swisscross;#gamma_{1} 1-E_{4}/E_{1}", 20, 0, 1.2);
   h_Photon1_severityLevel = new TH1F("h_Photon1_severityLevel","Photon1 severityLevel", 5, 0, 5);
   h_Photon1_severityLevel->GetXaxis()->SetBinLabel(1,"kGood");
   h_Photon1_severityLevel->GetXaxis()->SetBinLabel(2,"kProblematic");
   h_Photon1_severityLevel->GetXaxis()->SetBinLabel(3,"kRecovered");
   //   h_Photon1_severityLevel->GetXaxis()->SetBinLabel(4,"kTime");
   //   h_Photon1_severityLevel->GetXaxis()->SetBinLabel(5,"kWeird");
   //   h_Photon1_severityLevel->GetXaxis()->SetBinLabel(6,"kBad");
   h_Photon1_severityLevel->GetXaxis()->SetBinLabel(4,"kWeird");
   h_Photon1_severityLevel->GetXaxis()->SetBinLabel(5,"kBad");

   h_Photon1_recHitFlag = new TH1F("h_Photon1_recHitFlag","Photon1 recHitFlag", 15, 0, 15);

   // watch out, as this ordering changes with release
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(1,"kGood");                  // channel ok, the energy and time measurement are reliable
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(2,"kPoorReco");              // the energy is available from the UncalibRecHit, but approximate (bad shape, large chi2)
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(3,"kOutOfTime");             // the energy is available from the UncalibRecHit (sync reco), but the event is out of time
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(4,"kFaultyHardware");        // The energy is available from the UncalibRecHit, channel is faulty at some hardware level (e.g. noisy)
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(5,"kNoisy");                 // the channel is very noisy
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(6,"kPoorCalib");             // the energy is available from the UncalibRecHit, but the calibration of the channel is poor
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(7,"kSaturated");             // saturated channel (recovery not tried)
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(8,"kLeadingEdgeRecovered");  // saturated channel: energy estimated from the leading edge before saturation
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(9,"kNeighboursRecovered");   // saturated/isolated dead: energy estimated from neighbours
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(10,"kTowerRecovered");        // channel in TT with no data link, info retrieved from Trigger Primitive
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(11,"kFake");                  // the signal in the channel is a fake (e.g. a so-called spike)
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(12,"kFakeNeighbours");        // the signal in the channel is a fake and it is detected by looking at the neighbours
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(13,"kDead");                  // channel is dead and any recovery fails
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(14,"kKilled");                // MC only flag: the channel is killed in the real detector
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(15,"kUnknown");  

   h_Photon1_maxRecHitTime = new TH1F("h_Photon1_maxRecHitTime","Photon1 maxRecHitTime;#gamma_{1} time [ns]", 50, -5, 5); 
   h_Photon1_maxRecHitTime_wide = new TH1F("h_Photon1_maxRecHitTime_wide","Photon1 maxRecHitTime (zoom out);#gamma_{1} time [ns]", 200, -20, 20); 
   h_Photon1_e2e9_v_maxRecHitTime = new TH2F("h_Photon1_e2e9_v_maxRecHitTime","Photon1 E2/E9 v maxRecHitTime (zoom out);#gamma_{1} time [ns]; #gamma_{1} E_{2}/E_{9}", 200, -20, 20, 20, 0.4, 1.2); 
   
   h_Photon1_hadOverEm = new TH1F("h_Photon1_hadOverEm","Photon1 H/E;#gamma_{1} H/E", 40, 0, 0.06); 
   h_Photon1_hcalIso04 = new TH1F("h_Photon1_hcalIso04","Photon1 HCAL Iso #DeltaR 0.4;#gamma_{1} HCAL Iso #DeltaR 0.4", 50, 0, 4); 
   h_Photon1_hcalIso03 = new TH1F("h_Photon1_hcalIso03","Photon1 HCAL Iso #DeltaR 0.3;#gamma_{1} HCAL Iso #DeltaR 0.3", 50, 0, 4); 
   h_Photon1_ecalIso04 = new TH1F("h_Photon1_ecalIso04","Photon1 ECAL Iso #DeltaR 0.4;#gamma_{1} ECAL Iso #DeltaR 0.4", 30, -1, 4);
   h_Photon1_ecalIso03 = new TH1F("h_Photon1_ecalIso03","Photon1 ECAL Iso #DeltaR 0.3;#gamma_{1} ECAL Iso #DeltaR 0.3", 30, -1, 4); 
   
   h_Photon1_trkIsoSumPtHollow04 = new TH1F("h_Photon1_trkIsoSumPtHollow04","Photon1 Track Iso #Sigma p_{T} Hollow #DeltaR 0.4;#gamma_{1} Track Iso #Sigma p_{T} Hollow #DeltaR 0.4", 20, 0, 3); 
   h_Photon1_trkIsoSumPtSolid04 = new TH1F("h_Photon1_trkIsoSumPtSolid04","Photon1 Track Iso #Sigma p_{T} Solid #DeltaR 0.4;#gamma_{1} Track Iso #Sigma p_{T} Solid #DeltaR 0.4", 20, 0, 3); 
   h_Photon1_trkIsoSumPtHollow03 = new TH1F("h_Photon1_trkIsoSumPtHollow03","Photon1 Track Iso #Sigma p_{T} Hollow #DeltaR 0.3;#gamma_{1} Track Iso #Sigma p_{T} Hollow #DeltaR 0.3", 20, 0, 3); 
   h_Photon1_trkIsoSumPtSolid03 = new TH1F("h_Photon1_trkIsoSumPtSolid03","Photon1 Track Iso #Sigma p_{T} Solid #DeltaR 0.3;#gamma_{1} Track Iso #Sigma p_{T} Solid #DeltaR 0.3", 20, 0, 3); 

   h_Photon1_trkIsoNtrksHollow04 = new TH1F("h_Photon1_trkIsoNtrksHollow04","Photon1 Track Iso N_{tracks} Hollow #DeltaR 0.4;#gamma_{1} Track Iso N_{tracks} Hollow #DeltaR 0.4", 6, 0, 6); 
   h_Photon1_trkIsoNtrksSolid04 = new TH1F("h_Photon1_trkIsoNtrksSolid04","Photon1 Track Iso N_{tracks} Solid #DeltaR 0.4;#gamma_{1} Track Iso N_{tracks} Solid #DeltaR 0.4", 6, 0, 6); 
   h_Photon1_trkIsoNtrksHollow03 = new TH1F("h_Photon1_trkIsoNtrksHollow03","Photon1 Track Iso N_{tracks} Hollow #DeltaR 0.3;#gamma_{1} Track Iso N_{tracks} Hollow #DeltaR 0.3", 6, 0, 6); 
   h_Photon1_trkIsoNtrksSolid03 = new TH1F("h_Photon1_trkIsoNtrksSolid03","Photon1 Track Iso N_{tracks} Solid #DeltaR 0.3;#gamma_{1} Track Iso N_{tracks} Solid #DeltaR 0.3", 6, 0, 6); 

/*    h_Photon1_esRatio = new TH1F("h_Photon1_esRatio","Photon1_esRatio", , , ); */

   h_Photon2_pt = new TH1F("h_Photon2_pt","Photon2 pt;#gamma_{2} p_{T} [GeV]",44,60,500); 
   h_Photon2_pt_log = new TH1F("h_Photon2_pt_log","Photon2 pt;#gamma_{2} p_{T} [GeV]",44,60,500);
   h_Photon2_pt_zoom = new TH1F("h_Photon2_pt_zoom","Photon2 pt;#gamma_{2} p_{T} [GeV]",30,0,150); 
   h_Photon2_eta = new TH1F("h_Photon2_eta","Photon2 #eta;#gamma_{2} #eta",25,-2.5,2.5); 
   h_Photon2_phi = new TH1F ("h_Photon2_phi","Photon2 #phi;#gamma_{2} #phi",20,-3.14159,3.14159); 
   h_Photon2_occupancy = new TH2F ("h_Photon2_occupancy","Photon2 (#eta,#phi);#gamma_{2} #eta; #gamma_{2} #phi",25,-2.5,2.5,20,-3.14159,3.14159); 

   h_Photon2_r9 = new TH1F("h_Photon2_r9","Photon2 R9;#gamma_{2} R9", 50, 0.1, 1.5);
   h_Photon2_sigmaIetaIeta = new TH1F("h_Photon2_sigmaIetaIeta","Photon2 #sigma_{i#etai#eta};#gamma_{2} #sigma_{i#etai#eta}", 25, 0, 0.05);
   h_Photon2_sigmaEtaEta = new TH1F("h_Photon2_sigmaEtaEta","Photon2 #sigma_{#eta#eta};#gamma_{2} #sigma_{#eta#eta}", 25, 0, 0.05);

   h_Photon2_e2x2e4x4 = new TH1F("h_Photon2_e2x2e4x4","Photon2 E2x2/E4x4;#gamma_{2} E_{2x2}/E_{4x4}", 20, 0.4, 1.2);
   h_Photon2_e2e9 = new TH1F("h_Photon2_e2e9","Photon2 E2/E9;#gamma_{2} E_{2}/E_{9}", 20, 0.4, 1.2);
   h_Photon2_swisscross = new TH1F("h_Photon2_swisscross","Photon2 swisscross;#gamma_{2} 1-E_{4}/E_{1}", 20, 0, 1.2);
   h_Photon2_severityLevel = new TH1F("h_Photon2_severityLevel","Photon2 severityLevel", 5, 0, 5);
   h_Photon2_severityLevel->GetXaxis()->SetBinLabel(1,"kGood");
   h_Photon2_severityLevel->GetXaxis()->SetBinLabel(2,"kProblematic");
   h_Photon2_severityLevel->GetXaxis()->SetBinLabel(3,"kRecovered");
   //   h_Photon2_severityLevel->GetXaxis()->SetBinLabel(4,"kTime");
   //   h_Photon2_severityLevel->GetXaxis()->SetBinLabel(5,"kWeird");
   //   h_Photon2_severityLevel->GetXaxis()->SetBinLabel(6,"kBad");
   h_Photon2_severityLevel->GetXaxis()->SetBinLabel(4,"kWeird");
   h_Photon2_severityLevel->GetXaxis()->SetBinLabel(5,"kBad");

   h_Photon2_recHitFlag = new TH1F("h_Photon2_recHitFlag","Photon2 recHitFlag", 15, 0, 15);

   // watch out, as this ordering changes with release
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(1,"kGood");                  // channel ok, the energy and time measurement are reliable
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(2,"kPoorReco");              // the energy is available from the UncalibRecHit, but approximate (bad shape, large chi2)
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(3,"kOutOfTime");             // the energy is available from the UncalibRecHit (sync reco), but the event is out of time
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(4,"kFaultyHardware");        // The energy is available from the UncalibRecHit, channel is faulty at some hardware level (e.g. noisy)
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(5,"kNoisy");                 // the channel is very noisy
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(6,"kPoorCalib");             // the energy is available from the UncalibRecHit, but the calibration of the channel is poor
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(7,"kSaturated");             // saturated channel (recovery not tried)
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(8,"kLeadingEdgeRecovered");  // saturated channel: energy estimated from the leading edge before saturation
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(9,"kNeighboursRecovered");   // saturated/isolated dead: energy estimated from neighbours
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(10,"kTowerRecovered");        // channel in TT with no data link, info retrieved from Trigger Primitive
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(11,"kFake");                  // the signal in the channel is a fake (e.g. a so-called spike)
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(12,"kFakeNeighbours");        // the signal in the channel is a fake and it is detected by looking at the neighbours
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(13,"kDead");                  // channel is dead and any recovery fails
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(14,"kKilled");                // MC only flag: the channel is killed in the real detector
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(15,"kUnknown");  

   h_Photon2_maxRecHitTime = new TH1F("h_Photon2_maxRecHitTime","Photon2 maxRecHitTime;#gamma_{2} time [ns]", 50, -5, 5); 
   h_Photon2_maxRecHitTime_wide = new TH1F("h_Photon2_maxRecHitTime_wide","Photon2 maxRecHitTime (zoom out);#gamma_{2} time [ns]", 200, -20, 20); 
   h_Photon2_e2e9_v_maxRecHitTime = new TH2F("h_Photon2_e2e9_v_maxRecHitTime","Photon2 E2/E9 v maxRecHitTime (zoom out);#gamma_{2} time [ns]; #gamma_{2} E_{2}/E_{9}", 200, -20, 20, 20, 0.4, 1.2); 
   
   h_Photon2_hadOverEm = new TH1F("h_Photon2_hadOverEm","Photon2 H/E;#gamma_{2} H/E", 40, 0, 0.06); 
   h_Photon2_hcalIso04 = new TH1F("h_Photon2_hcalIso04","Photon2 HCAL Iso #DeltaR 0.4;#gamma_{2} HCAL Iso #DeltaR 0.4", 50, 0, 4); 
   h_Photon2_hcalIso03 = new TH1F("h_Photon2_hcalIso03","Photon2 HCAL Iso #DeltaR 0.3;#gamma_{2} HCAL Iso #DeltaR 0.3", 50, 0, 4); 
   h_Photon2_ecalIso04 = new TH1F("h_Photon2_ecalIso04","Photon2 ECAL Iso #DeltaR 0.4;#gamma_{2} ECAL Iso #DeltaR 0.4", 30, -1, 4);
   h_Photon2_ecalIso03 = new TH1F("h_Photon2_ecalIso03","Photon2 ECAL Iso #DeltaR 0.3;#gamma_{2} ECAL Iso #DeltaR 0.3", 30, -1, 4); 
   
   h_Photon2_trkIsoSumPtHollow04 = new TH1F("h_Photon2_trkIsoSumPtHollow04","Photon2 Track Iso #Sigma p_{T} Hollow #DeltaR 0.4;#gamma_{2} Track Iso #Sigma p_{T} Hollow #DeltaR 0.4", 20, 0, 3); 
   h_Photon2_trkIsoSumPtSolid04 = new TH1F("h_Photon2_trkIsoSumPtSolid04","Photon2 Track Iso #Sigma p_{T} Solid #DeltaR 0.4;#gamma_{2} Track Iso #Sigma p_{T} Solid #DeltaR 0.4", 20, 0, 3); 
   h_Photon2_trkIsoSumPtHollow03 = new TH1F("h_Photon2_trkIsoSumPtHollow03","Photon2 Track Iso #Sigma p_{T} Hollow #DeltaR 0.3;#gamma_{2} Track Iso #Sigma p_{T} Hollow #DeltaR 0.3", 20, 0, 3); 
   h_Photon2_trkIsoSumPtSolid03 = new TH1F("h_Photon2_trkIsoSumPtSolid03","Photon2 Track Iso #Sigma p_{T} Solid #DeltaR 0.3;#gamma_{2} Track Iso #Sigma p_{T} Solid #DeltaR 0.3", 20, 0, 3); 

   h_Photon2_trkIsoNtrksHollow04 = new TH1F("h_Photon2_trkIsoNtrksHollow04","Photon2 Track Iso N_{tracks} Hollow #DeltaR 0.4;#gamma_{2} Track Iso N_{tracks} Hollow #DeltaR 0.4", 6, 0, 6); 
   h_Photon2_trkIsoNtrksSolid04 = new TH1F("h_Photon2_trkIsoNtrksSolid04","Photon2 Track Iso N_{tracks} Solid #DeltaR 0.4;#gamma_{2} Track Iso N_{tracks} Solid #DeltaR 0.4", 6, 0, 6); 
   h_Photon2_trkIsoNtrksHollow03 = new TH1F("h_Photon2_trkIsoNtrksHollow03","Photon2 Track Iso N_{tracks} Hollow #DeltaR 0.3;#gamma_{2} Track Iso N_{tracks} Hollow #DeltaR 0.3", 6, 0, 6); 
   h_Photon2_trkIsoNtrksSolid03 = new TH1F("h_Photon2_trkIsoNtrksSolid03","Photon2 Track Iso N_{tracks} Solid #DeltaR 0.3;#gamma_{2} Track Iso N_{tracks} Solid #DeltaR 0.3", 6, 0, 6); 

/*    h_Photon2_esRatio = new TH1F("h_Photon2_esRatio","Photon2_esRatio", , , ); */

   // tight-tight ie signal sample
   h_FakeRate_tt_pt1  = new TH1F("h_FakeRate_tt_pt1","#gamma_{1} p_{T};#gamma_{1} p_{T}",44,60,500);
   h_FakeRate_tt_pt1_zoom  = new TH1F("h_FakeRate_tt_pt1_zoom","#gamma_{1} p_{T};#gamma_{1} p_{T}",30,0,150);
   h_FakeRate_tt_eta1 = new TH1F("h_FakeRate_tt_eta1","#gamma_{1} #eta;#gamma_{1} #eta",25,-2.5,2.5);
   h_FakeRate_tt_phi1 = new TH1F("h_FakeRate_tt_phi1","#gamma_{1} #phi;#gamma_{1} #phi",20,-3.14159,3.14159); 
   h_FakeRate_tt_pt2  = new TH1F("h_FakeRate_tt_pt2","#gamma_{2} p_{T};#gamma_{2} p_{T}",44,60,500);
   h_FakeRate_tt_pt2_zoom  = new TH1F("h_FakeRate_tt_pt2_zoom","#gamma_{2} p_{T};#gamma_{2} p_{T}",30,0,150);
   h_FakeRate_tt_eta2 = new TH1F("h_FakeRate_tt_eta2","#gamma_{2} #eta;#gamma_{2} #eta",25,-2.5,2.5);
   h_FakeRate_tt_phi2 = new TH1F("h_FakeRate_tt_phi2","#gamma_{2} #phi;#gamma_{2} #phi",20,-3.14159,3.14159); 
   //   h_FakeRate_tt_minv         = new TH1F("h_FakeRate_tt_minv",        "Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",43,140,1000);	    
   h_FakeRate_tt_minv         = new TH1F("h_FakeRate_tt_minv",        "Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",1860,140,2000);	    
h_FakeRate_tt_minv_high         = new TH1F("h_FakeRate_tt_minv_high",        "Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",35,500,1200);
   h_FakeRate_tt_qt           = new TH1F("h_FakeRate_tt_qt ",         "Diphoton qt;#gamma#gamma qt [GeV]",20,0,50);			    
   h_FakeRate_tt_deltaPhi     = new TH1F("h_FakeRate_tt_deltaPhi",    "Diphoton #Delta#phi;#gamma#gamma #Delta#phi",90,-3.14159,3.14159); 
   h_FakeRate_tt_deltaEta     = new TH1F("h_FakeRate_tt_deltaEta",    "Diphoton #Delta#eta;#gamma#gamma #Delta#eta",40,-6.,6.); 	    
   h_FakeRate_tt_deltaR       = new TH1F("h_FakeRate_tt_deltaR",      "Diphoton #DeltaR; #gamma#gamma #DeltaR",70,0,7.); 		    
   h_FakeRate_tt_cosThetaStar = new TH1F("h_FakeRate_tt_cosThetaStar","Diphoton |cos(#theta *)|; #gamma#gamma |cos#theta*|",20,0,1);      

   // tight-fake
   h_FakeRate_tf_pt1  = new TH1F("h_FakeRate_tf_pt1","#gamma_{1} p_{T};#gamma_{1} p_{T}",44,60,500);
   h_FakeRate_tf_pt1_zoom  = new TH1F("h_FakeRate_tf_pt1_zoom","#gamma_{1} p_{T};#gamma_{1} p_{T}",30,0,150);
   h_FakeRate_tf_eta1 = new TH1F("h_FakeRate_tf_eta1","#gamma_{1} #eta;#gamma_{1} #eta",25,-2.5,2.5);
   h_FakeRate_tf_phi1 = new TH1F("h_FakeRate_tf_phi1","#gamma_{1} #phi;#gamma_{1} #phi",20,-3.14159,3.14159); 
   h_FakeRate_tf_pt2  = new TH1F("h_FakeRate_tf_pt2","#gamma_{2} p_{T};#gamma_{2} p_{T}",44,60,500);
   h_FakeRate_tf_pt2_zoom  = new TH1F("h_FakeRate_tf_pt2_zoom","#gamma_{2} p_{T};#gamma_{2} p_{T}",30,0,150);
   h_FakeRate_tf_eta2 = new TH1F("h_FakeRate_tf_eta2","#gamma_{2} #eta;#gamma_{2} #eta",25,-2.5,2.5);
   h_FakeRate_tf_phi2 = new TH1F("h_FakeRate_tf_phi2","#gamma_{2} #phi;#gamma_{2} #phi",20,-3.14159,3.14159); 
   //   h_FakeRate_tf_minv         = new TH1F("h_FakeRate_tf_minv",        "Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",43,140,1000);	    
   h_FakeRate_tf_minv         = new TH1F("h_FakeRate_tf_minv",        "Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",1860,140,2000);	    
   h_FakeRate_tf_minv_high         = new TH1F("h_FakeRate_tf_minv_high",        "Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",35,500,1200);
   h_FakeRate_tf_qt           = new TH1F("h_FakeRate_tf_qt ",         "Diphoton qt;#gamma#gamma qt [GeV]",20,0,50);			    
   h_FakeRate_tf_deltaPhi     = new TH1F("h_FakeRate_tf_deltaPhi",    "Diphoton #Delta#phi;#gamma#gamma #Delta#phi",90,-3.14159,3.14159); 
   h_FakeRate_tf_deltaEta     = new TH1F("h_FakeRate_tf_deltaEta",    "Diphoton #Delta#eta;#gamma#gamma #Delta#eta",40,-6.,6.); 	    
   h_FakeRate_tf_deltaR       = new TH1F("h_FakeRate_tf_deltaR",      "Diphoton #DeltaR; #gamma#gamma #DeltaR",70,0,7.); 		    
   h_FakeRate_tf_cosThetaStar = new TH1F("h_FakeRate_tf_cosThetaStar","Diphoton |cos(#theta *)|; #gamma#gamma |cos#theta*|",20,0,1);      
   h_FakeRate_tf_pt1_noweight  = new TH1F("h_FakeRate_tf_pt1_noweight","TF #gamma_{1} (tight) p_{T}",100,0,150);
   h_FakeRate_tf_pt2_noweight  = new TH1F("h_FakeRate_tf_pt2_noweight","TF #gamma_{2} (fake) p_{T}",100,0,150);

   h_FakeRate_ft_pt1_noweight  = new TH1F("h_FakeRate_ft_pt1_noweight","FT #gamma_{1} (fake) p_{T}",100,0,150);
   h_FakeRate_ft_pt2_noweight  = new TH1F("h_FakeRate_ft_pt2_noweight","FT #gamma_{2} (tight) p_{T}",100,0,150);

   // fake-fake
   h_FakeRate_ff_pt1  = new TH1F("h_FakeRate_ff_pt1","#gamma_{1} p_{T};#gamma_{1} p_{T}",44,60,500);
   h_FakeRate_ff_pt1_zoom  = new TH1F("h_FakeRate_ff_pt1_zoom","#gamma_{1} p_{T};#gamma_{1} p_{T}",30,0,150);
   h_FakeRate_ff_eta1 = new TH1F("h_FakeRate_ff_eta1","#gamma_{1} #eta;#gamma_{1} #eta",25,-2.5,2.5);
   h_FakeRate_ff_phi1 = new TH1F("h_FakeRate_ff_phi1","#gamma_{1} #phi;#gamma_{1} #phi",20,-3.14159,3.14159); 
   h_FakeRate_ff_pt2  = new TH1F("h_FakeRate_ff_pt2","#gamma_{2} p_{T};#gamma_{2} p_{T}",44,60,500);
   h_FakeRate_ff_pt2_zoom  = new TH1F("h_FakeRate_ff_pt2_zoom","#gamma_{2} p_{T};#gamma_{2} p_{T}",30,0,150);
   h_FakeRate_ff_eta2 = new TH1F("h_FakeRate_ff_eta2","#gamma_{2} #eta;#gamma_{2} #eta",25,-2.5,2.5);
   h_FakeRate_ff_phi2 = new TH1F("h_FakeRate_ff_phi2","#gamma_{1} #phi;#gamma_{1} #phi",20,-3.14159,3.14159); 
   //   h_FakeRate_ff_minv         = new TH1F("h_FakeRate_ff_minv",        "Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",43,140,1000);	    
   h_FakeRate_ff_minv         = new TH1F("h_FakeRate_ff_minv",        "Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",1860,140,2000); 
   h_FakeRate_ff_minv_high    = new TH1F("h_FakeRate_ff_minv_high",        "Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",35,500,1200);
   h_FakeRate_ff_qt           = new TH1F("h_FakeRate_ff_qt ",         "Diphoton qt;#gamma#gamma qt [GeV]",20,0,50);			    
   h_FakeRate_ff_deltaPhi     = new TH1F("h_FakeRate_ff_deltaPhi",    "Diphoton #Delta#phi;#gamma#gamma #Delta#phi",90,-3.14159,3.14159); 
   h_FakeRate_ff_deltaEta     = new TH1F("h_FakeRate_ff_deltaEta",    "Diphoton #Delta#eta;#gamma#gamma #Delta#eta",40,-6.,6.); 	    
   h_FakeRate_ff_deltaR       = new TH1F("h_FakeRate_ff_deltaR",      "Diphoton #DeltaR; #gamma#gamma #DeltaR",70,0,7.); 		    
   h_FakeRate_ff_cosThetaStar = new TH1F("h_FakeRate_ff_cosThetaStar","Diphoton |cos(#theta *)|; #gamma#gamma |cos#theta*|",20,0,1);      
   h_FakeRate_ff_pt1_noweight  = new TH1F("h_FakeRate_ff_pt1_noweight","FF #gamma_{1} (fake) p_{T}",100,0,150);
   h_FakeRate_ff_pt2_noweight  = new TH1F("h_FakeRate_ff_pt2_noweight","FF #gamma_{2} (fake) p_{T}",100,0,150);   
   h_trkIsoVRho1 = new TH2F("h_trkIsoVRho1","Photon1 Track Iso;#rho; #gamma_{1} Track Iso", 10, 0, 10, 14, 0, 14); 
   h_ecalIsoVRho1 = new TH2F("h_ecalIsoVRho1","Photon1 Ecal Iso;#rho; #gamma_{1} Ecal Iso", 10, 0, 10, 14, 0, 14); 
   h_hcalIsoVRho1 = new TH2F("h_hcalIsoVRho1","Photon1 Hcal Iso;#rho; #gamma_{1} Hcal Iso", 10, 0, 10, 14, 0, 14); 

   h_trkIsoVRho2 = new TH2F("h_trkIsoVRho2","Photon2 Track Iso;#rho; #gamma_{2} Track Iso", 10, 0, 10, 14, 0, 14); 
   h_ecalIsoVRho2 = new TH2F("h_ecalIsoVRho2","Photon2 Ecal Iso;#rho; #gamma_{2} Ecal Iso", 10, 0, 10, 14, 0, 14); 
   h_hcalIsoVRho2 = new TH2F("h_hcalIsoVRho2","Photon2 Hcal Iso;#rho; #gamma_{2} Hcal Iso", 10, 0, 10, 14, 0, 14); 

   h_FakeRate_tt_minv_120to200 = new TH1F("h_FakeRate_tt_minv_120to200","Diphoton Invariant Mass (120-200 bin);M_{#gamma#gamma} [GeV/c^2]",1,120,200); 
   h_FakeRate_tt_minv_200to500 = new TH1F("h_FakeRate_tt_minv_200to500","Diphoton Invariant Mass (200-500 bin);M_{#gamma#gamma} [GeV/c^2]",1,200,500); 
   h_FakeRate_tt_minv_500to800 = new TH1F("h_FakeRate_tt_minv_500to800","Diphoton Invariant Mass (500-800 bin);M_{#gamma#gamma} [GeV/c^2]",1,500,800); 
   h_FakeRate_tt_minv_800toInf = new TH1F("h_FakeRate_tt_minv_800toInf","Diphoton Invariant Mass (800-Inf bin);M_{#gamma#gamma} [GeV/c^2]",1,800,100000); 

   h_FakeRate_tf_minv_120to200 = new TH1F("h_FakeRate_tf_minv_120to200","Diphoton Invariant Mass (120-200 bin);M_{#gamma#gamma} [GeV/c^2]",1,120,200); 
   h_FakeRate_tf_minv_200to500 = new TH1F("h_FakeRate_tf_minv_200to500","Diphoton Invariant Mass (200-500 bin);M_{#gamma#gamma} [GeV/c^2]",1,200,500); 
   h_FakeRate_tf_minv_500to800 = new TH1F("h_FakeRate_tf_minv_500to800","Diphoton Invariant Mass (500-800 bin);M_{#gamma#gamma} [GeV/c^2]",1,500,800);
   h_FakeRate_tf_minv_800toInf = new TH1F("h_FakeRate_tf_minv_800toInf","Diphoton Invariant Mass (800-Inf bin);M_{#gamma#gamma} [GeV/c^2]",1,800,100000); 

   h_FakeRate_ff_minv_120to200 = new TH1F("h_FakeRate_ff_minv_120to200","Diphoton Invariant Mass (120-200 bin);M_{#gamma#gamma} [GeV/c^2]",1,120,200); 
   h_FakeRate_ff_minv_200to500 = new TH1F("h_FakeRate_ff_minv_200to500","Diphoton Invariant Mass (200-500 bin);M_{#gamma#gamma} [GeV/c^2]",1,200,500); 
   h_FakeRate_ff_minv_500to800 = new TH1F("h_FakeRate_ff_minv_500to800","Diphoton Invariant Mass (500-800 bin);M_{#gamma#gamma} [GeV/c^2]",1,500,800); 
   h_FakeRate_ff_minv_800toInf = new TH1F("h_FakeRate_ff_minv_800toInf","Diphoton Invariant Mass (800-Inf bin);M_{#gamma#gamma} [GeV/c^2]",1,800,100000); 

   h_Diphoton_Minv->Sumw2();
   h_Diphoton_Minv_Yousi5->Sumw2();
   h_Diphoton_Minv_Yousi10->Sumw2();
   h_Diphoton_Minv_Yousi40->Sumw2();
   h_Diphoton_Minv_Yousi100->Sumw2();
   h_Diphoton_Minv_120to200->Sumw2();
   h_Diphoton_Minv_200to500->Sumw2();
   h_Diphoton_Minv_500to800->Sumw2();
   h_Diphoton_Minv_800toInf->Sumw2();
   h_Diphoton_Minv_log->Sumw2();
   h_Diphoton_Minv_high->Sumw2();
   h_Diphoton_Minv_low->Sumw2();
   h_FakeRate_tt_pt1->Sumw2();
   h_FakeRate_tt_pt1_zoom->Sumw2();
   h_FakeRate_tt_eta1->Sumw2();
   h_FakeRate_tt_phi1->Sumw2();
   h_FakeRate_tt_pt2->Sumw2();
   h_FakeRate_tt_pt2_zoom->Sumw2();
   h_FakeRate_tt_eta2->Sumw2();
   h_FakeRate_tt_phi2->Sumw2();
   h_FakeRate_tt_minv->Sumw2();
   h_FakeRate_tt_minv_high->Sumw2();
   h_FakeRate_tt_qt->Sumw2();
   h_FakeRate_tt_deltaPhi->Sumw2();
   h_FakeRate_tt_deltaEta->Sumw2();
   h_FakeRate_tt_deltaR->Sumw2();
   h_FakeRate_tt_cosThetaStar->Sumw2();
   
   h_FakeRate_tf_pt1->Sumw2();
   h_FakeRate_tf_pt1_zoom->Sumw2();
   h_FakeRate_tf_eta1->Sumw2();
   h_FakeRate_tf_phi1->Sumw2();
   h_FakeRate_tf_pt2->Sumw2();
   h_FakeRate_tf_pt2_zoom->Sumw2();
   h_FakeRate_tf_eta2->Sumw2();
   h_FakeRate_tf_phi2->Sumw2();
   h_FakeRate_tf_minv->Sumw2();
   h_FakeRate_tf_minv_high->Sumw2();
   h_FakeRate_tf_qt->Sumw2();
   h_FakeRate_tf_deltaPhi->Sumw2();
   h_FakeRate_tf_deltaEta->Sumw2();
   h_FakeRate_tf_deltaR->Sumw2();
   h_FakeRate_tf_cosThetaStar->Sumw2();
   
   h_FakeRate_tf_pt1_noweight->Sumw2();
   h_FakeRate_tf_pt2_noweight->Sumw2();

   h_FakeRate_ft_pt1_noweight->Sumw2();
   h_FakeRate_ft_pt2_noweight->Sumw2();

   h_FakeRate_ff_pt1->Sumw2();
   h_FakeRate_ff_pt1_zoom->Sumw2();
   h_FakeRate_ff_eta1->Sumw2();
   h_FakeRate_ff_phi1->Sumw2();
   h_FakeRate_ff_pt2->Sumw2();
   h_FakeRate_ff_pt2_zoom->Sumw2();
   h_FakeRate_ff_eta2->Sumw2();
   h_FakeRate_ff_phi2->Sumw2();
   h_FakeRate_ff_minv->Sumw2();   
   h_FakeRate_ff_minv_high->Sumw2();   
   h_FakeRate_ff_qt->Sumw2();
   h_FakeRate_ff_deltaPhi->Sumw2();
   h_FakeRate_ff_deltaEta->Sumw2();
   h_FakeRate_ff_deltaR->Sumw2();
   h_FakeRate_ff_cosThetaStar->Sumw2();

   h_FakeRate_ff_pt1_noweight->Sumw2();
   h_FakeRate_ff_pt2_noweight->Sumw2();

   h_FakeRate_tt_minv_120to200->Sumw2();
   h_FakeRate_tt_minv_200to500->Sumw2();
   h_FakeRate_tt_minv_500to800->Sumw2();
   h_FakeRate_tt_minv_800toInf->Sumw2();
   h_FakeRate_tf_minv_120to200->Sumw2();
   h_FakeRate_tf_minv_200to500->Sumw2();
   h_FakeRate_tf_minv_500to800->Sumw2();   
   h_FakeRate_tf_minv_800toInf->Sumw2();
   h_FakeRate_ff_minv_120to200->Sumw2();
   h_FakeRate_ff_minv_200to500->Sumw2();
   h_FakeRate_ff_minv_500to800->Sumw2();
   h_FakeRate_ff_minv_800toInf->Sumw2();

   _outT0 = new TTree("tDiphoton","Diphoton events passing selection");
   TBranch* b0 = _outT0->Branch("minv",&outVar0.minv,"minv/D");
   b0 = _outT0->Branch("weight",&outVar0.weight,"weight/D");
   b0 = _outT0->Branch("r9_photon1",&outVar0.r9_photon1,"r9_photon1/D");
   b0 = _outT0->Branch("r9_photon2",&outVar0.r9_photon2,"r9_photon2/D");
   b0 = _outT0->Branch("eta_photon1",&outVar0.eta_photon1,"eta_photon1/D");
   b0 = _outT0->Branch("eta_photon2",&outVar0.eta_photon2,"eta_photon2/D");

   fInKfactor = new TFile(_kFactorFile) ;

   return;
}

Bool_t fTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void fTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t fTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
  std::cout << entry << std::endl;
   return 1;
}
#endif // #ifdef fTree_cxx
