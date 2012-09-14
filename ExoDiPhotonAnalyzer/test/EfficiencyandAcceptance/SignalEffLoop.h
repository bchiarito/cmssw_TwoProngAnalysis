//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep 13 23:52:35 2012 by ROOT version 5.34/00
// from TTree fTree/PhotonTree
// found on file: diphoton_tree_RSGravToGG_kMpl-01_M-3250_TuneZ2star_8TeV-pythia6_Septh10thPileUp.root
//////////////////////////////////////////////////////////

#ifndef SignalEffLoop_h
#define SignalEffLoop_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class SignalEffLoop {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Declaration of leaf types
  Int_t           Event_run;
  Int_t           Event_LS;
  Int_t           Event_evnum;
  Double_t        Vtx_vx;
  Double_t        Vtx_vy;
  Double_t        Vtx_vz;
  Double_t        Vtx_sumPtTracks;
  Double_t        Vtx_ndof;
  Double_t        Vtx_d0;
  Int_t           Vtx_Nvtx;
  Int_t           Vtx_Ntracks;
  Bool_t          Vtx_isFake;
  Double_t        BeamSpot_x0;
  Double_t        BeamSpot_y0;
  Double_t        BeamSpot_z0;
  Double_t        BeamSpot_sigmaZ;
  Double_t        BeamSpot_x0error;
  Double_t        BeamSpot_y0error;
  Double_t        BeamSpot_z0error;
  Double_t        BeamSpot_sigmaZ0error;
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
  Int_t           TrigHLT_HLT_Photon20_CaloIdVL_IsoL_v5;
  Int_t           TrigHLT_HLT_Photon20_R9Id_Photon18_R9Id_v6;
  Int_t           TrigHLT_HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v6;
  Int_t           TrigHLT_HLT_Photon26_Photon18_v6;
  Int_t           TrigHLT_HLT_Photon26_IsoVL_Photon18_v6;
  Int_t           TrigHLT_HLT_Photon26_IsoVL_Photon18_IsoVL_v6;
  Int_t           TrigHLT_HLT_Photon26_CaloIdL_IsoVL_Photon18_v6;
  Int_t           TrigHLT_HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v5;
  Int_t           TrigHLT_HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v6;
  Int_t           TrigHLT_HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v5;
  Int_t           TrigHLT_HLT_Photon26_R9Id_Photon18_R9Id_v3;
  Int_t           TrigHLT_HLT_Photon30_CaloIdVL_v6;
  Int_t           TrigHLT_HLT_Photon30_CaloIdVL_IsoL_v6;
  Int_t           TrigHLT_HLT_Photon36_IsoVL_Photon22_v3;
  Int_t           TrigHLT_HLT_Photon36_CaloIdVL_Photon22_CaloIdVL_v1;
  Int_t           TrigHLT_HLT_Photon36_CaloIdL_Photon22_CaloIdL_v5;
  Int_t           TrigHLT_HLT_Photon36_CaloIdL_IsoVL_Photon22_v3;
  Int_t           TrigHLT_HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_v2;
  Int_t           TrigHLT_HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v2;
  Int_t           TrigHLT_HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v1;
  Int_t           TrigHLT_HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v2;
  Int_t           TrigHLT_HLT_Photon36_R9Id_Photon22_R9Id_v2;
  Int_t           TrigHLT_HLT_Photon40_CaloIdL_Photon28_CaloIdL_v3;
  Int_t           TrigHLT_HLT_Photon44_CaloIdL_Photon34_CaloIdL_v1;
  Int_t           TrigHLT_HLT_Photon48_CaloIdL_Photon38_CaloIdL_v1;
  Int_t           TrigHLT_HLT_Photon50_CaloIdVL_v3;
  Int_t           TrigHLT_HLT_Photon50_CaloIdVL_IsoL_v5;
  Int_t           TrigHLT_HLT_Photon70_CaloIdL_HT350_v6;
  Int_t           TrigHLT_HLT_Photon70_CaloIdL_HT400_v1;
  Int_t           TrigHLT_HLT_Photon70_CaloIdL_MHT70_v6;
  Int_t           TrigHLT_HLT_Photon70_CaloIdL_MHT90_v1;
  Int_t           TrigHLT_HLT_Photon75_CaloIdVL_v6;
  Int_t           TrigHLT_HLT_Photon75_CaloIdVL_IsoL_v6;
  Int_t           TrigHLT_HLT_Photon90_CaloIdVL_v3;
  Int_t           TrigHLT_HLT_Photon90_CaloIdVL_IsoL_v3;
  Int_t           TrigHLT_HLT_Photon135_v1;
  Int_t           TrigHLT_HLT_Photon225_NoHE_v1;
  Int_t           TrigHLT_HLT_Photon400_v1;
  Int_t           TrigHLT_HLT_Photon200_NoHE_v3;
  Int_t           TrigHLT_HLT_DoublePhoton33_HEVT_v3;
  Int_t           TrigHLT_HLT_DoublePhoton38_HEVT_v2;
  Int_t           TrigHLT_HLT_DoublePhoton60_v3;
  Int_t           TrigHLT_HLT_DoublePhoton80_v1;
  Int_t           TrigHLT_HLT_DoublePhoton5_IsoVL_CEP_v5;
  Int_t           TrigHLT_HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60_v5;
  Int_t           TrigHLT_HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60_v6;
  Int_t           TrigHLT_HLT_Photon26_CaloId10_Iso50_Photon18_R9Id85_Mass60_v5;
  Int_t           TrigHLT_HLT_Photon26_CaloId10_Iso50_Photon18_R9Id85_Mass60_v6;
  Int_t           TrigHLT_HLT_Photon26_Photon18_v11;
  Int_t           TrigHLT_HLT_Photon26_Photon18_v12;
  Int_t           TrigHLT_HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v5;
  Int_t           TrigHLT_HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v6;
  Int_t           TrigHLT_HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v1;
  Int_t           TrigHLT_HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v2;
  Int_t           TrigHLT_HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_v4;
  Int_t           TrigHLT_HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_v5;
  Int_t           TrigHLT_HLT_Photon26_R9Id85_Photon18_CaloId10_Iso50_Mass60_v5;
  Int_t           TrigHLT_HLT_Photon26_R9Id85_Photon18_CaloId10_Iso50_Mass60_v6;
  Int_t           TrigHLT_HLT_Photon26_R9Id85_Photon18_R9Id85_Mass60_v3;
  Int_t           TrigHLT_HLT_Photon26_R9Id85_Photon18_R9Id85_Mass60_v4;
  Int_t           TrigHLT_HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v5;
  Int_t           TrigHLT_HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v6;
  Int_t           TrigHLT_HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v5;
  Int_t           TrigHLT_HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v6;
  Int_t           TrigHLT_HLT_Photon36_Photon22_v5;
  Int_t           TrigHLT_HLT_Photon36_Photon22_v6;
  Int_t           TrigHLT_HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v5;
  Int_t           TrigHLT_HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v6;
  Int_t           TrigHLT_HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_v4;
  Int_t           TrigHLT_HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_v5;
  Int_t           TrigHLT_HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v5;
  Int_t           TrigHLT_HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v6;
  Int_t           TrigHLT_HLT_Photon36_R9Id85_Photon22_R9Id85_v3;
  Int_t           TrigHLT_HLT_Photon36_R9Id85_Photon22_R9Id85_v4;
  Int_t           TrigHLT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v6;
  Int_t           TrigHLT_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v7;
  Int_t           TrigHLT_HLT_DoubleEle33_CaloIdL_v13;
  Int_t           TrigHLT_HLT_DoubleEle33_CaloIdL_v14;
  Int_t           TrigHLT_HLT_DoubleEle33_CaloIdT_v10;
  Int_t           TrigHLT_HLT_DoubleEle33_CaloIdT_v9;
  Int_t           TrigHLT_HLT_DoublePhoton40_CaloIdL_Rsq0p035_v3;
  Int_t           TrigHLT_HLT_DoublePhoton40_CaloIdL_Rsq0p035_v4;
  Int_t           TrigHLT_HLT_DoublePhoton40_CaloIdL_Rsq0p035_v5;
  Int_t           TrigHLT_HLT_DoublePhoton40_CaloIdL_Rsq0p035_v6;
  Int_t           TrigHLT_HLT_DoublePhoton40_CaloIdL_Rsq0p06_v3;
  Int_t           TrigHLT_HLT_DoublePhoton40_CaloIdL_Rsq0p06_v4;
  Int_t           TrigHLT_HLT_DoublePhoton40_CaloIdL_Rsq0p06_v5;
  Int_t           TrigHLT_HLT_DoublePhoton40_CaloIdL_Rsq0p06_v6;
  Int_t           TrigHLT_HLT_DoublePhoton48_HEVT_v7;
  Int_t           TrigHLT_HLT_DoublePhoton48_HEVT_v8;
  Int_t           TrigHLT_HLT_DoublePhoton53_HEVT_v1;
  Int_t           TrigHLT_HLT_DoublePhoton53_HEVT_v2;
  Int_t           TrigHLT_HLT_DoublePhoton70_v5;
  Int_t           TrigHLT_HLT_DoublePhoton70_v6;
  Int_t           TrigHLT_HLT_DoublePhoton80_v6;
  Int_t           TrigHLT_HLT_DoublePhoton80_v7;
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
  Double_t        Photon1_sigmaIphiIphi;
  Double_t        Photon1_sigmaPhiPhi;
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
  Double_t        Photon2_sigmaIphiIphi;
  Double_t        Photon2_sigmaPhiPhi;
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
  Double_t        DiphotonGen_Minv;
  Double_t        DiphotonGen_qt;
  Double_t        DiphotonGen_deltaPhi;
  Double_t        DiphotonGen_deltaEta;
  Double_t        DiphotonGen_deltaR;
  Double_t        DiphotonGen_cosThetaStar;
  Double_t        Diphoton_Minv;
  Double_t        Diphoton_qt;
  Double_t        Diphoton_deltaPhi;
  Double_t        Diphoton_deltaEta;
  Double_t        Diphoton_deltaR;
  Double_t        Diphoton_cosThetaStar;
  Double_t        rho25;
  Int_t           pu_n;
  Double_t        MCPUWeight;

  // Histograms for Acceptance and Efficiency
  
  TH1F*  h_Diphoton_Minv;
  TH1F*  h_Diphoton_Minv_AccPassed;
  TH1F*  h_Diphoton_Minv_EffandAccPassed;
  
  TFile* _outfilename;
  
  Int_t ptmin;
  Int_t ptmax;
  Int_t nMinvbins;
  Int_t nptbins;
  
  Int_t Minvmin;
  Int_t Minvmax;
  
  
  // List of branches
  TBranch        *b_Event;   //!
  TBranch        *b_Vtx;   //!
  TBranch        *b_BeamSpot;   //!
  TBranch        *b_TrigHLT;   //!
  TBranch        *b_GenPhoton1;   //!
  TBranch        *b_GenPhoton2;   //!
  TBranch        *b_Photon1;   //!
  TBranch        *b_Photon2;   //!
  TBranch        *b_DiphotonGen;   //!
  TBranch        *b_Diphoton;   //!
  TBranch        *b_rho25;   //!
  TBranch        *b_pu_n;   //!
  TBranch        *b_MCPUWeight;   //!
  
  SignalEffLoop(TTree *tree=0);
  virtual ~SignalEffLoop();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef SignalEffLoop_cxx
SignalEffLoop::SignalEffLoop(TTree *tree) 
: fChain(0) 
, _outfilename(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("diphoton_tree_RSGravToGG_kMpl-01_M-3250_TuneZ2star_8TeV-pythia6_Septh10thPileUp.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("diphoton_tree_RSGravToGG_kMpl-01_M-3250_TuneZ2star_8TeV-pythia6_Septh10thPileUp.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("diphoton_tree_RSGravToGG_kMpl-01_M-3250_TuneZ2star_8TeV-pythia6_Septh10thPileUp.root:/diphotonSignalMCAnalyzer");
      dir->GetObject("fTree",tree);

   }
   Init(tree);
}

SignalEffLoop::~SignalEffLoop()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SignalEffLoop::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SignalEffLoop::LoadTree(Long64_t entry)
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

void SignalEffLoop::Init(TTree *tree)
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
   fChain->SetBranchAddress("Vtx", &Vtx_vx, &b_Vtx);
   fChain->SetBranchAddress("BeamSpot", &BeamSpot_x0, &b_BeamSpot);
   fChain->SetBranchAddress("TrigHLT", &TrigHLT_HLT_MinBiasBSC, &b_TrigHLT);
   fChain->SetBranchAddress("GenPhoton1", &GenPhoton1_status, &b_GenPhoton1);
   fChain->SetBranchAddress("GenPhoton2", &GenPhoton2_status, &b_GenPhoton2);
   fChain->SetBranchAddress("Photon1", &Photon1_pt, &b_Photon1);
   fChain->SetBranchAddress("Photon2", &Photon2_pt, &b_Photon2);
   fChain->SetBranchAddress("DiphotonGen", &DiphotonGen_Minv, &b_DiphotonGen);
   fChain->SetBranchAddress("Diphoton", &Diphoton_Minv, &b_Diphoton);
   fChain->SetBranchAddress("rho25", &rho25, &b_rho25);
   fChain->SetBranchAddress("pu_n", &pu_n, &b_pu_n);
   fChain->SetBranchAddress("MCPUWeight", &MCPUWeight, &b_MCPUWeight);
   Notify();
 
   // Set Axis Ranges here                                                                   

   ptmin = 20;
   ptmax= 3000;
   nptbins = (ptmax-ptmin)/20 ;
   Minvmin = 2800;
   Minvmax = 3800;
   nMinvbins = (Minvmax-Minvmin)/20;



   h_Diphoton_Minv = new TH1F("h_Diphoton_Minv","Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^{2}]",nMinvbins,Minvmin,Minvmax);


   h_Diphoton_Minv_AccPassed = new TH1F("h_Diphoton_Minv_AccPassed","Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^{2}]",nMinvbins,Minvmin,Minvmax);    
   
   
   h_Diphoton_Minv_EffandAccPassed = new TH1F("h_Diphoton_Minv_EffandAccPassed","Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^{2}]",nMinvbins,Minvmin,Minvmax);


}

Bool_t SignalEffLoop::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SignalEffLoop::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SignalEffLoop::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef SignalEffLoop_cxx
