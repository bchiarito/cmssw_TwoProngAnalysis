#ifndef TRIG_INFO_INC
#define TRIG_INFO_INC

//********************************************************************
// Definition of a struct that can be used for storing trig info
// in a tree, from different analysers
// Also includes a Fill function to fill the struct from the appropriate objects
// and a string that can be used to define the tree branch
// 
// $Id: TriggerInfo.h,v 1.11 2011/06/06 22:48:46 yma Exp $ 
// 
//********************************************************************

#include <string>

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

namespace ExoDiPhotons
{

  struct l1TrigInfo_t{

    // the following L1 tech bits are usually selected upon in GOODCOLL skim
    // so store the info here for later
    bool L1_Tech0; // BPTX AND
    bool L1_Tech36;  //L1Tech_BSC_halo_beam2_inner.v0
    bool L1_Tech37; //L1Tech_BSC_halo_beam2_outer.v0
    bool L1_Tech38; //L1Tech_BSC_halo_beam1_inner.v0
    bool L1_Tech39; //L1Tech_BSC_halo_beam1_outer.v0
    bool L1_Tech40; //L1Tech_BSC_minBias_threshold1.v0
    bool L1_Tech41; //L1Tech_BSC_minBias_threshold2.v0
    bool L1_Tech42; // L1Tech_BSC_splash_beam1.v0
    bool L1_Tech43; // L1Tech_BSC_splash_beam2.v0
    bool L1_EG2; // also L1 EG bits ?
    bool L1_EG5;
    bool L1_EG8;
  };

  // string for defining a tree branch for this struct
  std::string l1TrigBranchDefString("L1_Tech0/O:L1_Tech36:L1_Tech37:L1_Tech38:L1_Tech39:L1_Tech40:L1_Tech41:L1_Tech42:L1_Tech43:L1_EG2:L1_EG5:L1_EG8");


  // ideally would want a Fill function for L1 info too...



  // HLT info 
  // do far going with: MinBias, single and double photon, L1_EGs,
  // but not electron or photon+X
  // this list will need to be kept up to date since the HLT triggers evolve
  // eg with spike cleaning, prescales, etc
  // note that we use an int,, where 0=fail and 1=pass obviously, 
  // and -1=not present in menu for this event
  //   should we  store prescale value too?
  struct hltTrigInfo_t{
    int HLT_MinBiasBSC;
    int HLT_MinBiasBSC_NoBPTX;
    int HLT_MinBiasBSC_OR;
    int HLT_L1_BscMinBiasOR_BptxPlusORMinus;
    int HLT_L1SingleEG2;
    int HLT_L1SingleEG5;
    int HLT_L1SingleEG8;
    int HLT_L1DoubleEG5;
    int HLT_Photon10_L1R;
    int HLT_Photon10_Cleaned_L1R;
    int HLT_Photon15_L1R;
    int HLT_Photon15_Cleaned_L1R;
    int HLT_Photon15_LooseEcalIso_L1R;
    int HLT_Photon15_LooseEcalIso_Cleaned_L1R;
    int HLT_Photon15_TrackIso_L1R;
    int HLT_Photon15_TrackIso_Cleaned_L1R;
    int HLT_Photon17_Isol_SC17HE_L1R_v1;
    int HLT_Photon17_SC17HE_L1R_v1;
    int HLT_Photon20_L1R;
    int HLT_Photon20_Cleaned_L1R;
    int HLT_Photon20_NoHE_L1R;
    int HLT_Photon22_SC22HE_L1R_v1;
    int HLT_Photon25_Cleaned_L1R;
    int HLT_Photon30_L1R;
    int HLT_Photon30_Cleaned_L1R;
    int HLT_Photon30_L1R_8E29;
    int HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1;
    int HLT_Photon35_Isol_Cleaned_L1R_v1;
    int HLT_Photon40_CaloId_Cleaned_L1R_v1;
    int HLT_Photon40_Isol_Cleaned_L1R_v1;
    int HLT_Photon50_L1R;
    int HLT_Photon50_Cleaned_L1R;
    int HLT_Photon50_Cleaned_L1R_v1;
    int HLT_Photon50_NoHE_L1R;
    int HLT_Photon50_NoHE_Cleaned_L1R;
    int HLT_Photon70_Cleaned_L1R_v1;
    int HLT_Photon70_NoHE_Cleaned_L1R_v1;
    int HLT_Photon100_NoHE_Cleaned_L1R_v1;
    int HLT_Photon110_NoHE_Cleaned_L1R_v1;
    int HLT_DoublePhoton5_L1R;
    int HLT_DoublePhoton5_CEP_L1R;
    int HLT_DoublePhoton5_CEP_L1R_v3;
    int HLT_DoublePhoton5_Jpsi_L1R;
    int HLT_DoublePhoton5_Upsilon_L1R;
    int HLT_DoublePhoton10_L1R;
    int HLT_DoublePhoton15_L1R;
    int HLT_DoublePhoton17_L1R;
    int HLT_DoublePhoton17_SingleIsol_L1R_v1;
    int HLT_DoublePhoton20_L1R;
    int HLT_DoublePhoton22_L1R_v1;    
    //newly added, 2011A
    int HLT_DoublePhoton33_v1;
    int HLT_DoublePhoton33_v2;
    int HLT_DoublePhoton33_v3;
    int HLT_DoublePhoton33_v5;
    int HLT_DoublePhoton33_HEVT_v2;
    int HLT_DoublePhoton5_IsoVL_CEP_v1;
    int HLT_DoublePhoton5_IsoVL_CEP_v2;
    int HLT_Photon125_NoSpikeFilter_v1;
    int HLT_Photon125_NoSpikeFilter_v2;
    int HLT_Photon125_NoSpikeFilter_v3;
    int HLT_Photon20_CaloIdVL_IsoL_v1;
    int HLT_Photon20_CaloIdVL_IsoL_v2;
    int HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1;
    int HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2;
    int HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v3;
    int HLT_Photon20_EBOnly_NoSpikeFilter_v1;
    int HLT_Photon20_NoSpikeFilter_v1;
    int HLT_Photon20_R9Id_Photon18_R9Id_v1;
    int HLT_Photon20_R9Id_Photon18_R9Id_v2;
    int HLT_Photon20_R9Id_Photon18_R9Id_v3;
    int HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v1;
    int HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v2;
    int HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v3;
    int HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v1;
    int HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v2;
    int HLT_Photon26_CaloIdL_IsoVL_Photon18_v1;
    int HLT_Photon26_CaloIdL_IsoVL_Photon18_v2;
    int HLT_Photon26_CaloIdL_IsoVL_Photon18_v3;
    int HLT_Photon26_IsoVL_Photon18_IsoVL_v1;
    int HLT_Photon26_IsoVL_Photon18_IsoVL_v2;
    int HLT_Photon26_IsoVL_Photon18_IsoVL_v3;
    int HLT_Photon26_IsoVL_Photon18_v1;
    int HLT_Photon26_IsoVL_Photon18_v2;
    int HLT_Photon26_IsoVL_Photon18_v3;
    int HLT_Photon26_Photon18_v1;
    int HLT_Photon26_Photon18_v2;
    int HLT_Photon26_Photon18_v3;
    int HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v1;
    int HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v2;
    int HLT_Photon30_CaloIdVL_IsoL_v1;
    int HLT_Photon30_CaloIdVL_IsoL_v2;
    int HLT_Photon30_CaloIdVL_IsoL_v3;
    int HLT_Photon30_CaloIdVL_v1;
    int HLT_Photon30_CaloIdVL_v2;
    int HLT_Photon30_CaloIdVL_v3;
    int HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1;
    int HLT_Photon32_CaloIdL_Photon26_CaloIdL_v2;
    int HLT_Photon32_CaloIdL_Photon26_CaloIdL_v3;
    int HLT_Photon36_CaloIdL_Photon22_CaloIdL_v1;
    int HLT_Photon36_CaloIdL_Photon22_CaloIdL_v2;
    int HLT_Photon50_CaloIdVL_IsoL_v1;
    int HLT_Photon50_CaloIdVL_IsoL_v2;
    int HLT_Photon75_CaloIdVL_IsoL_v1;
    int HLT_Photon75_CaloIdVL_IsoL_v2;
    int HLT_Photon75_CaloIdVL_IsoL_v3;
    int HLT_Photon75_CaloIdVL_v1;
    int HLT_Photon75_CaloIdVL_v2;
    int HLT_Photon75_CaloIdVL_v3;
    int HLT_DoublePhoton40_MR150_v3;                                
    int HLT_DoublePhoton40_R014_MR150_v3;
    int HLT_DoublePhoton50_v2;
    int HLT_DoublePhoton5_IsoVL_CEP_v4;
    int HLT_DoublePhoton60_v2;
    int HLT_Mu15_DoublePhoton15_CaloIdL_v6;
    int HLT_Mu15_Photon20_CaloIdL_v6;
    int HLT_Mu8_Photon20_CaloIdVT_IsoT_v5;
    int HLT_Photon125_v2;
    int HLT_Photon200_NoHE_v2;
    int HLT_Photon20_CaloIdVL_IsoL_v4;
    int HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v5;
    int HLT_Photon20_R9Id_Photon18_R9Id_v5;
    int HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v5;
    int HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v4;
    int HLT_Photon26_CaloIdL_IsoVL_Photon18_v5;
    int HLT_Photon26_IsoVL_Photon18_IsoVL_v5;
    int HLT_Photon26_IsoVL_Photon18_v5;
    int HLT_Photon26_Photon18_v5;
    int HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v4;
    int HLT_Photon26_R9Id_Photon18_R9Id_v2;
    int HLT_Photon30_CaloIdVL_IsoL_v5;
    int HLT_Photon30_CaloIdVL_v5;
    int HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v1;
    int HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_v1;
    int HLT_Photon36_CaloIdL_IsoVL_Photon22_v2;
    int HLT_Photon36_CaloIdL_Photon22_CaloIdL_v4;
    int HLT_Photon36_CaloId_IsoVL_Photon22_R9Id_v1;
    int HLT_Photon36_IsoVL_Photon22_v2;
    int HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v1;
    int HLT_Photon36_R9Id_Photon22_R9Id_v1;
    int HLT_Photon40_CaloIdL_Photon28_CaloIdL_v2;
    int HLT_Photon40_R005_MR150_v3;
    int HLT_Photon40_R014_MR450_v3;
    int HLT_Photon40_R020_MR300_v3;
    int HLT_Photon40_R025_MR200_v3;
    int HLT_Photon40_R038_MR150_v3;
    int HLT_Photon50_CaloIdVL_IsoL_v4;
    int HLT_Photon50_CaloIdVL_v2;
    int HLT_Photon70_CaloIdL_HT300_v6;
    int HLT_Photon70_CaloIdL_HT350_v5;
    int HLT_Photon70_CaloIdL_MHT50_v6;
    int HLT_Photon70_CaloIdL_MHT70_v5;
    int HLT_Photon75_CaloIdVL_IsoL_v5;
    int HLT_Photon75_CaloIdVL_v5;
    int HLT_Photon90_CaloIdVL_IsoL_v2;
    int HLT_Photon90_CaloIdVL_v2;




    int HLT_DoublePhoton33_HEVT_v1;
    int HLT_DoublePhoton33_v4;
    int HLT_DoublePhoton50_v1;
    int HLT_DoublePhoton5_IsoVL_CEP_v3;
    int HLT_DoublePhoton60_v1;
    int HLT_Photon125_v1;
    int HLT_Photon200_NoHE_v1;
    int HLT_Photon20_CaloIdVL_IsoL_v3;
    int HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v4;
    int HLT_Photon20_R9Id_Photon18_R9Id_v4;
    int HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v4;
    int HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v3;
    int HLT_Photon26_CaloIdL_IsoVL_Photon18_v4;
    int HLT_Photon26_IsoVL_Photon18_IsoVL_v4;
    int HLT_Photon26_IsoVL_Photon18_v4;
    int HLT_Photon26_Photon18_v4;
    int HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v3;
    int HLT_Photon26_R9Id_Photon18_R9Id_v1;
    int HLT_Photon30_CaloIdVL_IsoL_v4;
    int HLT_Photon30_CaloIdVL_v4;
    int HLT_Photon32_CaloIdL_Photon26_CaloIdL_v4;
    int HLT_Photon36_CaloIdL_IsoVL_Photon22_v1;
    int HLT_Photon36_CaloIdL_Photon22_CaloIdL_v3;
    int HLT_Photon36_IsoVL_Photon22_v1;
    int HLT_Photon40_CaloIdL_Photon28_CaloIdL_v1;
    int HLT_Photon50_CaloIdVL_IsoL_v3;
    int HLT_Photon50_CaloIdVL_v1;
    int HLT_Photon75_CaloIdVL_IsoL_v4;
    int HLT_Photon75_CaloIdVL_v4;
    int HLT_Photon90_CaloIdVL_IsoL_v1;
    int HLT_Photon90_CaloIdVL_v1;
    
  };

  // string for defining a tree branch for this struct
  std::string hltTrigBranchDefString("HLT_MinBiasBSC/I:HLT_MinBiasBSC_NoBPTX:HLT_MinBiasBSC_OR:HLT_L1_BscMinBiasOR_BptxPlusORMinus:HLT_L1SingleEG2:HLT_L1SingleEG5:HLT_L1SingleEG8:HLT_L1DoubleEG5:HLT_Photon10_L1R:HLT_Photon10_Cleaned_L1R:HLT_Photon15_L1R:HLT_Photon15_Cleaned_L1R:HLT_Photon15_LooseEcalIso_L1R:HLT_Photon15_LooseEcalIso_Cleaned_L1R:HLT_Photon15_TrackIso_L1R:HLT_Photon15_TrackIso_Cleaned_L1R:HLT_Photon17_Isol_SC17HE_L1R_v1:HLT_Photon17_SC17HE_L1R_v1:HLT_Photon20_L1R:HLT_Photon20_Cleaned_L1R:HLT_Photon20_NoHE_L1R:HLT_Photon22_SC22HE_L1R_v1:HLT_Photon25_Cleaned_L1R:HLT_Photon30_L1R:HLT_Photon30_Cleaned_L1R:HLT_Photon30_L1R_8E29:HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1:HLT_Photon35_Isol_Cleaned_L1R_v1:HLT_Photon40_CaloId_Cleaned_L1R_v1:HLT_Photon40_Isol_Cleaned_L1R_v1:HLT_Photon50_L1R:HLT_Photon50_Cleaned_L1R:HLT_Photon50_Cleaned_L1R_v1:HLT_Photon50_NoHE_L1R:HLT_Photon50_NoHE_Cleaned_L1R:HLT_Photon70_Cleaned_L1R_v1:HLT_Photon70_NoHE_Cleaned_L1R_v1:HLT_Photon100_NoHE_Cleaned_L1R_v1:HLT_Photon110_NoHE_Cleaned_L1R_v1:HLT_DoublePhoton5_L1R:HLT_DoublePhoton5_CEP_L1R:HLT_DoublePhoton5_CEP_L1R_v3:HLT_DoublePhoton5_Jpsi_L1R:HLT_DoublePhoton5_Upsilon_L1R:HLT_DoublePhoton10_L1R:HLT_DoublePhoton15_L1R:HLT_DoublePhoton17_L1R:HLT_DoublePhoton17_SingleIsol_L1R_v1:HLT_DoublePhoton20_L1R:HLT_DoublePhoton22_L1R_v1:HLT_DoublePhoton33_v1:HLT_DoublePhoton33_v2:HLT_DoublePhoton33_v3:HLT_DoublePhoton33_v5:HLT_DoublePhoton33_HEVT_v2:HLT_DoublePhoton5_IsoVL_CEP_v1:HLT_DoublePhoton5_IsoVL_CEP_v2:HLT_Photon125_NoSpikeFilter_v1:HLT_Photon125_NoSpikeFilter_v2:HLT_Photon125_NoSpikeFilter_v3:HLT_Photon20_CaloIdVL_IsoL_v1:HLT_Photon20_CaloIdVL_IsoL_v2:HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1:HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2:HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v3:HLT_Photon20_EBOnly_NoSpikeFilter_v1:HLT_Photon20_NoSpikeFilter_v1:HLT_Photon20_R9Id_Photon18_R9Id_v1:HLT_Photon20_R9Id_Photon18_R9Id_v2:HLT_Photon20_R9Id_Photon18_R9Id_v3:HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v1:HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v2:HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v3:HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v1:HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v2:HLT_Photon26_CaloIdL_IsoVL_Photon18_v1:HLT_Photon26_CaloIdL_IsoVL_Photon18_v2:HLT_Photon26_CaloIdL_IsoVL_Photon18_v3:HLT_Photon26_IsoVL_Photon18_IsoVL_v1:HLT_Photon26_IsoVL_Photon18_IsoVL_v2:HLT_Photon26_IsoVL_Photon18_IsoVL_v3:HLT_Photon26_IsoVL_Photon18_v1:HLT_Photon26_IsoVL_Photon18_v2:HLT_Photon26_IsoVL_Photon18_v3:HLT_Photon26_Photon18_v1:HLT_Photon26_Photon18_v2:HLT_Photon26_Photon18_v3:HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v1:HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v2:HLT_Photon30_CaloIdVL_IsoL_v1:HLT_Photon30_CaloIdVL_IsoL_v2:HLT_Photon30_CaloIdVL_IsoL_v3:HLT_Photon30_CaloIdVL_v1:HLT_Photon30_CaloIdVL_v2:HLT_Photon30_CaloIdVL_v3:HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1:HLT_Photon32_CaloIdL_Photon26_CaloIdL_v2:HLT_Photon32_CaloIdL_Photon26_CaloIdL_v3:HLT_Photon36_CaloIdL_Photon22_CaloIdL_v1:HLT_Photon36_CaloIdL_Photon22_CaloIdL_v2:HLT_Photon50_CaloIdVL_IsoL_v1:HLT_Photon50_CaloIdVL_IsoL_v2:HLT_Photon75_CaloIdVL_IsoL_v1:HLT_Photon75_CaloIdVL_IsoL_v2:HLT_Photon75_CaloIdVL_IsoL_v3:HLT_Photon75_CaloIdVL_v1:HLT_Photon75_CaloIdVL_v2:HLT_Photon75_CaloIdVL_v3:HLT_DoublePhoton40_MR150_v3:HLT_DoublePhoton40_R014_MR150_v3:HLT_DoublePhoton50_v2:HLT_DoublePhoton5_IsoVL_CEP_v4:HLT_DoublePhoton60_v2:HLT_Mu15_DoublePhoton15_CaloIdL_v6:HLT_Mu15_Photon20_CaloIdL_v6:HLT_Mu8_Photon20_CaloIdVT_IsoT_v5:HLT_Photon125_v2:HLT_Photon200_NoHE_v2:HLT_Photon20_CaloIdVL_IsoL_v4:HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v5:HLT_Photon20_R9Id_Photon18_R9Id_v5:HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v5:HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v4:HLT_Photon26_CaloIdL_IsoVL_Photon18_v5:HLT_Photon26_IsoVL_Photon18_IsoVL_v5:HLT_Photon26_IsoVL_Photon18_v5:HLT_Photon26_Photon18_v5:HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v4:HLT_Photon26_R9Id_Photon18_R9Id_v2:HLT_Photon30_CaloIdVL_IsoL_v5:HLT_Photon30_CaloIdVL_v5:HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v1:HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_v1:HLT_Photon36_CaloIdL_IsoVL_Photon22_v2:HLT_Photon36_CaloIdL_Photon22_CaloIdL_v4:HLT_Photon36_CaloId_IsoVL_Photon22_R9Id_v1:HLT_Photon36_IsoVL_Photon22_v2:HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v1:HLT_Photon36_R9Id_Photon22_R9Id_v1:HLT_Photon40_CaloIdL_Photon28_CaloIdL_v2:HLT_Photon40_R005_MR150_v3:HLT_Photon40_R014_MR450_v3:HLT_Photon40_R020_MR300_v3:HLT_Photon40_R025_MR200_v3:HLT_Photon40_R038_MR150_v3:HLT_Photon50_CaloIdVL_IsoL_v4:HLT_Photon50_CaloIdVL_v2:HLT_Photon70_CaloIdL_HT300_v6:HLT_Photon70_CaloIdL_HT350_v5:HLT_Photon70_CaloIdL_MHT50_v6:HLT_Photon70_CaloIdL_MHT70_v5:HLT_Photon75_CaloIdVL_IsoL_v5:HLT_Photon75_CaloIdVL_v5:HLT_Photon90_CaloIdVL_IsoL_v2:HLT_Photon90_CaloIdVL_v2:HLT_DoublePhoton33_HEVT_v1:HLT_DoublePhoton33_v4:HLT_DoublePhoton50_v1:HLT_DoublePhoton5_IsoVL_CEP_v3:HLT_DoublePhoton60_v1:HLT_Photon125_v1:HLT_Photon200_NoHE_v1:HLT_Photon20_CaloIdVL_IsoL_v3:HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v4:HLT_Photon20_R9Id_Photon18_R9Id_v4:HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v4:HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v3:HLT_Photon26_CaloIdL_IsoVL_Photon18_v4:HLT_Photon26_IsoVL_Photon18_IsoVL_v4:HLT_Photon26_IsoVL_Photon18_v4:HLT_Photon26_Photon18_v4:HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v3:HLT_Photon26_R9Id_Photon18_R9Id_v1:HLT_Photon30_CaloIdVL_IsoL_v4:HLT_Photon30_CaloIdVL_v4:HLT_Photon32_CaloIdL_Photon26_CaloIdL_v4:HLT_Photon36_CaloIdL_IsoVL_Photon22_v1:HLT_Photon36_CaloIdL_Photon22_CaloIdL_v3:HLT_Photon36_IsoVL_Photon22_v1:HLT_Photon40_CaloIdL_Photon28_CaloIdL_v1:HLT_Photon50_CaloIdVL_IsoL_v3:HLT_Photon50_CaloIdVL_v1:HLT_Photon75_CaloIdVL_IsoL_v4:HLT_Photon75_CaloIdVL_v4:HLT_Photon90_CaloIdVL_IsoL_v1:HLT_Photon90_CaloIdVL_v1");

  // need also to have an Initialise() function?
  // No, we can do this inside the Fill function ...

  // Fill() function for HLT info
  // needs both a TriggerResults and TriggerNames object

  void FillHLTInfo(hltTrigInfo_t &hltInfo, const edm::TriggerResults *hltResults, const edm::TriggerNames &hltNames) {

    // first set all the triggers to -1
    // then if the path is not even present in the menu for this event
    // it will remain at -1, while regular pass/fail are 1/0

    hltInfo.HLT_MinBiasBSC = -1;
    hltInfo.HLT_MinBiasBSC_NoBPTX = -1;
    hltInfo.HLT_MinBiasBSC_OR = -1;
    hltInfo.HLT_L1_BscMinBiasOR_BptxPlusORMinus = -1;
    hltInfo.HLT_L1SingleEG2 = -1;
    hltInfo.HLT_L1SingleEG5 = -1;
    hltInfo.HLT_L1SingleEG8 = -1;
    hltInfo.HLT_L1DoubleEG5 = -1;
    hltInfo.HLT_Photon10_L1R = -1;
    hltInfo.HLT_Photon10_Cleaned_L1R = -1;
    hltInfo.HLT_Photon15_L1R = -1;
    hltInfo.HLT_Photon15_Cleaned_L1R = -1;
    hltInfo.HLT_Photon15_LooseEcalIso_L1R = -1;
    hltInfo.HLT_Photon15_LooseEcalIso_Cleaned_L1R = -1;
    hltInfo.HLT_Photon15_TrackIso_L1R = -1;
    hltInfo.HLT_Photon15_TrackIso_Cleaned_L1R = -1;
    hltInfo.HLT_Photon17_Isol_SC17HE_L1R_v1 = -1;
    hltInfo.HLT_Photon17_SC17HE_L1R_v1 = -1;
    hltInfo.HLT_Photon20_L1R = -1;
    hltInfo.HLT_Photon20_Cleaned_L1R = -1;
    hltInfo.HLT_Photon20_NoHE_L1R = -1;
    hltInfo.HLT_Photon22_SC22HE_L1R_v1 = -1;
    hltInfo.HLT_Photon25_Cleaned_L1R = -1;
    hltInfo.HLT_Photon30_L1R = -1;
    hltInfo.HLT_Photon30_Cleaned_L1R = -1;
    hltInfo.HLT_Photon30_L1R_8E29 = -1;
    hltInfo.HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1 = -1;
    hltInfo.HLT_Photon35_Isol_Cleaned_L1R_v1 = -1;
    hltInfo.HLT_Photon40_CaloId_Cleaned_L1R_v1 = -1;
    hltInfo.HLT_Photon40_Isol_Cleaned_L1R_v1 = -1;
    hltInfo.HLT_Photon50_L1R = -1;
    hltInfo.HLT_Photon50_Cleaned_L1R = -1;
    hltInfo.HLT_Photon50_Cleaned_L1R_v1 = -1;
    hltInfo.HLT_Photon50_NoHE_L1R = -1;
    hltInfo.HLT_Photon50_NoHE_Cleaned_L1R = -1;
    hltInfo.HLT_Photon70_Cleaned_L1R_v1 = -1;
    hltInfo.HLT_Photon70_NoHE_Cleaned_L1R_v1 = -1;
    hltInfo.HLT_Photon100_NoHE_Cleaned_L1R_v1 = -1;
    hltInfo.HLT_Photon110_NoHE_Cleaned_L1R_v1 = -1;
    hltInfo.HLT_DoublePhoton5_L1R = -1;
    hltInfo.HLT_DoublePhoton5_CEP_L1R = -1;
    hltInfo.HLT_DoublePhoton5_CEP_L1R_v3 = -1;
    hltInfo.HLT_DoublePhoton5_Jpsi_L1R = -1;
    hltInfo.HLT_DoublePhoton5_Upsilon_L1R = -1;
    hltInfo.HLT_DoublePhoton10_L1R = -1;
    hltInfo.HLT_DoublePhoton15_L1R = -1;
    hltInfo.HLT_DoublePhoton17_L1R = -1;
    hltInfo.HLT_DoublePhoton17_SingleIsol_L1R_v1 = -1;
    hltInfo.HLT_DoublePhoton20_L1R = -1;
    hltInfo.HLT_DoublePhoton22_L1R_v1 = -1;
    //new 2011A
    hltInfo.HLT_DoublePhoton33_v1 = -1;
    hltInfo.HLT_DoublePhoton33_v2 = -1;
    hltInfo.HLT_DoublePhoton33_v3 = -1;
    hltInfo.HLT_DoublePhoton33_v5 = -1;
    hltInfo.HLT_DoublePhoton33_HEVT_v2 = -1;
    hltInfo.HLT_DoublePhoton5_IsoVL_CEP_v1 = -1;
    hltInfo.HLT_DoublePhoton5_IsoVL_CEP_v2 = -1;
    hltInfo.HLT_Photon125_NoSpikeFilter_v1 = -1;
    hltInfo.HLT_Photon125_NoSpikeFilter_v2 = -1;
    hltInfo.HLT_Photon125_NoSpikeFilter_v3 = -1;
    hltInfo.HLT_Photon20_CaloIdVL_IsoL_v1 = -1;
    hltInfo.HLT_Photon20_CaloIdVL_IsoL_v2 = -1;
    hltInfo.HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1 = -1;
    hltInfo.HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2 = -1;
    hltInfo.HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v3 = -1;
    hltInfo.HLT_Photon20_EBOnly_NoSpikeFilter_v1 = -1;
    hltInfo.HLT_Photon20_NoSpikeFilter_v1 = -1;
    hltInfo.HLT_Photon20_R9Id_Photon18_R9Id_v1 = -1;
    hltInfo.HLT_Photon20_R9Id_Photon18_R9Id_v2 = -1;
    hltInfo.HLT_Photon20_R9Id_Photon18_R9Id_v3 = -1;
    hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v1 = -1;
    hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v2 = -1;
    hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v3 = -1;
    hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v1 = -1;
    hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v2 = -1;
    hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_v1 = -1;
    hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_v2 = -1;
    hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_v3 = -1;
    hltInfo.HLT_Photon26_IsoVL_Photon18_IsoVL_v1 = -1;
    hltInfo.HLT_Photon26_IsoVL_Photon18_IsoVL_v2 = -1;
    hltInfo.HLT_Photon26_IsoVL_Photon18_IsoVL_v3 = -1;
    hltInfo.HLT_Photon26_IsoVL_Photon18_v1 = -1;
    hltInfo.HLT_Photon26_IsoVL_Photon18_v2 = -1;
    hltInfo.HLT_Photon26_IsoVL_Photon18_v3 = -1;
    hltInfo.HLT_Photon26_Photon18_v1 = -1;
    hltInfo.HLT_Photon26_Photon18_v2 = -1;
    hltInfo.HLT_Photon26_Photon18_v3 = -1;
    hltInfo.HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v1 = -1;
    hltInfo.HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v2 = -1;
    hltInfo.HLT_Photon30_CaloIdVL_IsoL_v1 = -1;
    hltInfo.HLT_Photon30_CaloIdVL_IsoL_v2 = -1;
    hltInfo.HLT_Photon30_CaloIdVL_IsoL_v3 = -1;
    hltInfo.HLT_Photon30_CaloIdVL_v1 = -1;
    hltInfo.HLT_Photon30_CaloIdVL_v2 = -1;
    hltInfo.HLT_Photon30_CaloIdVL_v3 = -1;
    hltInfo.HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1 = -1;
    hltInfo.HLT_Photon32_CaloIdL_Photon26_CaloIdL_v2 = -1;
    hltInfo.HLT_Photon32_CaloIdL_Photon26_CaloIdL_v3 = -1;
    hltInfo.HLT_Photon36_CaloIdL_Photon22_CaloIdL_v1 = -1;
    hltInfo.HLT_Photon36_CaloIdL_Photon22_CaloIdL_v2 = -1;
    hltInfo.HLT_Photon50_CaloIdVL_IsoL_v1 = -1;
    hltInfo.HLT_Photon50_CaloIdVL_IsoL_v2 = -1;
    hltInfo.HLT_Photon75_CaloIdVL_IsoL_v1 = -1;
    hltInfo.HLT_Photon75_CaloIdVL_IsoL_v2 = -1;
    hltInfo.HLT_Photon75_CaloIdVL_IsoL_v3 = -1;
    hltInfo.HLT_Photon75_CaloIdVL_v1 = -1;
    hltInfo.HLT_Photon75_CaloIdVL_v2 = -1;
    hltInfo.HLT_Photon75_CaloIdVL_v3 = -1;
    hltInfo.HLT_DoublePhoton40_MR150_v3 = -1;
    hltInfo.HLT_DoublePhoton40_R014_MR150_v3 = -1;
    hltInfo.HLT_DoublePhoton50_v2 = -1;
    hltInfo.HLT_DoublePhoton5_IsoVL_CEP_v4 = -1;
    hltInfo.HLT_DoublePhoton60_v2 = -1;
    hltInfo.HLT_Mu15_DoublePhoton15_CaloIdL_v6 = -1;
    hltInfo.HLT_Mu15_Photon20_CaloIdL_v6 = -1;
    hltInfo.HLT_Mu8_Photon20_CaloIdVT_IsoT_v5 = -1;
    hltInfo.HLT_Photon125_v2 = -1;
    hltInfo.HLT_Photon200_NoHE_v2 = -1;
    hltInfo.HLT_Photon20_CaloIdVL_IsoL_v4 = -1;
    hltInfo.HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v5 = -1;
    hltInfo.HLT_Photon20_R9Id_Photon18_R9Id_v5 = -1;
    hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v5 = -1;
    hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v4 = -1;
    hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_v5 = -1;
    hltInfo.HLT_Photon26_IsoVL_Photon18_IsoVL_v5 = -1;
    hltInfo.HLT_Photon26_IsoVL_Photon18_v5 = -1;
    hltInfo.HLT_Photon26_Photon18_v5 = -1;
    hltInfo.HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v4 = -1;
    hltInfo.HLT_Photon26_R9Id_Photon18_R9Id_v2 = -1;
    hltInfo.HLT_Photon30_CaloIdVL_IsoL_v5 = -1;
    hltInfo.HLT_Photon30_CaloIdVL_v5 = -1;
    hltInfo.HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v1 = -1;
    hltInfo.HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_v1 = -1;
    hltInfo.HLT_Photon36_CaloIdL_IsoVL_Photon22_v2 = -1;
    hltInfo.HLT_Photon36_CaloIdL_Photon22_CaloIdL_v4 = -1;
    hltInfo.HLT_Photon36_CaloId_IsoVL_Photon22_R9Id_v1 = -1;
    hltInfo.HLT_Photon36_IsoVL_Photon22_v2 = -1;
    hltInfo.HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v1 = -1;
    hltInfo.HLT_Photon36_R9Id_Photon22_R9Id_v1 = -1;
    hltInfo.HLT_Photon40_CaloIdL_Photon28_CaloIdL_v2 = -1;
    hltInfo.HLT_Photon40_R005_MR150_v3 = -1;
    hltInfo.HLT_Photon40_R014_MR450_v3 = -1;
    hltInfo.HLT_Photon40_R020_MR300_v3 = -1;
    hltInfo.HLT_Photon40_R025_MR200_v3 = -1;
    hltInfo.HLT_Photon40_R038_MR150_v3 = -1;
    hltInfo.HLT_Photon50_CaloIdVL_IsoL_v4 = -1;
    hltInfo.HLT_Photon50_CaloIdVL_v2 = -1;
    hltInfo.HLT_Photon70_CaloIdL_HT300_v6 = -1;
    hltInfo.HLT_Photon70_CaloIdL_HT350_v5 = -1;
    hltInfo.HLT_Photon70_CaloIdL_MHT50_v6 = -1;
    hltInfo.HLT_Photon70_CaloIdL_MHT70_v5 = -1;
    hltInfo.HLT_Photon75_CaloIdVL_IsoL_v5 = -1;
    hltInfo.HLT_Photon75_CaloIdVL_v5 = -1;
    hltInfo.HLT_Photon90_CaloIdVL_IsoL_v2 = -1;
    hltInfo.HLT_Photon90_CaloIdVL_v2 = -1;
    hltInfo.HLT_DoublePhoton33_HEVT_v1 = -1;
    hltInfo.HLT_DoublePhoton33_v4 = -1;
    hltInfo.HLT_DoublePhoton50_v1 = -1;
    hltInfo.HLT_DoublePhoton5_IsoVL_CEP_v3 = -1;
    hltInfo.HLT_DoublePhoton60_v1 = -1;
    hltInfo.HLT_Photon125_v1 = -1;
    hltInfo.HLT_Photon200_NoHE_v1 = -1;
    hltInfo.HLT_Photon20_CaloIdVL_IsoL_v3 = -1;
    hltInfo.HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v4 = -1;
    hltInfo.HLT_Photon20_R9Id_Photon18_R9Id_v4 = -1;
    hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v4 = -1;
    hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v3 = -1;
    hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_v4 = -1;
    hltInfo.HLT_Photon26_IsoVL_Photon18_IsoVL_v4 = -1;
    hltInfo.HLT_Photon26_IsoVL_Photon18_v4 = -1;
    hltInfo.HLT_Photon26_Photon18_v4 = -1;
    hltInfo.HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v3 = -1;
    hltInfo.HLT_Photon26_R9Id_Photon18_R9Id_v1 = -1;
    hltInfo.HLT_Photon30_CaloIdVL_IsoL_v4 = -1;
    hltInfo.HLT_Photon30_CaloIdVL_v4 = -1;
    hltInfo.HLT_Photon32_CaloIdL_Photon26_CaloIdL_v4 = -1;
    hltInfo.HLT_Photon36_CaloIdL_IsoVL_Photon22_v1 = -1;
    hltInfo.HLT_Photon36_CaloIdL_Photon22_CaloIdL_v3 = -1;
    hltInfo.HLT_Photon36_IsoVL_Photon22_v1 = -1;
    hltInfo.HLT_Photon40_CaloIdL_Photon28_CaloIdL_v1 = -1;
    hltInfo.HLT_Photon50_CaloIdVL_IsoL_v3 = -1;
    hltInfo.HLT_Photon50_CaloIdVL_v1 = -1;
    hltInfo.HLT_Photon75_CaloIdVL_IsoL_v4 = -1;
    hltInfo.HLT_Photon75_CaloIdVL_v4 = -1;
    hltInfo.HLT_Photon90_CaloIdVL_IsoL_v1 = -1;
    hltInfo.HLT_Photon90_CaloIdVL_v1 = -1;
    
    // we'll just loop over all triggers in the current event
    // and check one-by-one the result for the triggers we are interested in

    for(int itrig=0;itrig<(int)hltNames.size();itrig++) {

      //      std::cout << hltNames.triggerName(itrig) << " " << hltResults->accept(itrig) << std::endl;

      if(hltNames.triggerName(itrig)=="HLT_MinBiasBSC")
	hltInfo.HLT_MinBiasBSC = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_MinBiasBSC_NoBPTX")
        hltInfo.HLT_MinBiasBSC_NoBPTX = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_MinBiasBSC_OR")
        hltInfo.HLT_MinBiasBSC_OR = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_L1_BscMinBiasOR_BptxPlusORMinus")
        hltInfo.HLT_L1_BscMinBiasOR_BptxPlusORMinus = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_L1SingleEG2")
        hltInfo.HLT_L1SingleEG2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_L1SingleEG5")
        hltInfo.HLT_L1SingleEG5 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_L1SingleEG8")
        hltInfo.HLT_L1SingleEG8 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_L1DoubleEG5")
        hltInfo.HLT_L1DoubleEG5 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon10_L1R")
        hltInfo.HLT_Photon10_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon10_Cleaned_L1R")
        hltInfo.HLT_Photon10_Cleaned_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon15_L1R")
        hltInfo.HLT_Photon15_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon15_Cleaned_L1R")
        hltInfo.HLT_Photon15_Cleaned_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon15_LooseEcalIso_L1R")
        hltInfo.HLT_Photon15_LooseEcalIso_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon15_LooseEcalIso_Cleaned_L1R")
        hltInfo.HLT_Photon15_LooseEcalIso_Cleaned_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon15_TrackIso_L1R")
        hltInfo.HLT_Photon15_TrackIso_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon15_TrackIso_Cleaned_L1R")
        hltInfo.HLT_Photon15_TrackIso_Cleaned_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon17_Isol_SC17HE_L1R_v1")
        hltInfo.HLT_Photon17_Isol_SC17HE_L1R_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon17_SC17HE_L1R_v1")
        hltInfo.HLT_Photon17_SC17HE_L1R_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon20_L1R")
        hltInfo.HLT_Photon20_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon20_Cleaned_L1R")
        hltInfo.HLT_Photon20_Cleaned_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon20_NoHE_L1R")
        hltInfo.HLT_Photon20_NoHE_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon22_SC22HE_L1R_v1")
        hltInfo.HLT_Photon22_SC22HE_L1R_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon25_Cleaned_L1R")
        hltInfo.HLT_Photon25_Cleaned_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon30_L1R")
        hltInfo.HLT_Photon30_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon30_Cleaned_L1R")
        hltInfo.HLT_Photon30_Cleaned_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon30_L1R_8E29")
        hltInfo.HLT_Photon30_L1R_8E29 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1")
        hltInfo.HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon35_Isol_Cleaned_L1R_v1")
        hltInfo.HLT_Photon35_Isol_Cleaned_L1R_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon40_CaloId_Cleaned_L1R_v1")
        hltInfo.HLT_Photon40_CaloId_Cleaned_L1R_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon40_Isol_Cleaned_L1R_v1")
        hltInfo.HLT_Photon40_Isol_Cleaned_L1R_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon50_L1R")
        hltInfo.HLT_Photon50_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon50_Cleaned_L1R")
        hltInfo.HLT_Photon50_Cleaned_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon50_Cleaned_L1R_v1")
        hltInfo.HLT_Photon50_Cleaned_L1R_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon50_NoHE_L1R")
        hltInfo.HLT_Photon50_NoHE_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon50_NoHE_Cleaned_L1R")
        hltInfo.HLT_Photon50_NoHE_Cleaned_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon70_Cleaned_L1R_v1")
        hltInfo.HLT_Photon70_Cleaned_L1R_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon70_NoHE_Cleaned_L1R_v1")
        hltInfo.HLT_Photon70_NoHE_Cleaned_L1R_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon100_NoHE_Cleaned_L1R_v1")
        hltInfo.HLT_Photon100_NoHE_Cleaned_L1R_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon110_NoHE_Cleaned_L1R_v1")
        hltInfo.HLT_Photon110_NoHE_Cleaned_L1R_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton5_L1R")
        hltInfo.HLT_DoublePhoton5_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton5_CEP_L1R")
        hltInfo.HLT_DoublePhoton5_CEP_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton5_CEP_L1R_v3")
        hltInfo.HLT_DoublePhoton5_CEP_L1R_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton5_Jpsi_L1R")
        hltInfo.HLT_DoublePhoton5_Jpsi_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton5_Upsilon_L1R")
        hltInfo.HLT_DoublePhoton5_Upsilon_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton10_L1R")
        hltInfo.HLT_DoublePhoton10_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton15_L1R")
        hltInfo.HLT_DoublePhoton15_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton17_L1R")
        hltInfo.HLT_DoublePhoton17_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton17_SingleIsol_L1R_v1")
        hltInfo.HLT_DoublePhoton17_SingleIsol_L1R_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton20_L1R")
        hltInfo.HLT_DoublePhoton20_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton22_L1R_v1")
        hltInfo.HLT_DoublePhoton22_L1R_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton33_v1")
	hltInfo.HLT_DoublePhoton33_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton33_v2")
	hltInfo.HLT_DoublePhoton33_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton33_v3")
	hltInfo.HLT_DoublePhoton33_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton33_v5")
	hltInfo.HLT_DoublePhoton33_v5 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton33_HEVT_v2")
	hltInfo.HLT_DoublePhoton33_HEVT_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton5_IsoVL_CEP_v1")
	hltInfo.HLT_DoublePhoton5_IsoVL_CEP_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton5_IsoVL_CEP_v2")
	hltInfo.HLT_DoublePhoton5_IsoVL_CEP_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon125_NoSpikeFilter_v1")
	hltInfo.HLT_Photon125_NoSpikeFilter_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon125_NoSpikeFilter_v2")
	hltInfo.HLT_Photon125_NoSpikeFilter_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon125_NoSpikeFilter_v3")
	hltInfo.HLT_Photon125_NoSpikeFilter_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon20_CaloIdVL_IsoL_v1")
	hltInfo.HLT_Photon20_CaloIdVL_IsoL_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon20_CaloIdVL_IsoL_v2")
	hltInfo.HLT_Photon20_CaloIdVL_IsoL_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1")
	hltInfo.HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2")
	hltInfo.HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v3")
	hltInfo.HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon20_EBOnly_NoSpikeFilter_v1")
	hltInfo.HLT_Photon20_EBOnly_NoSpikeFilter_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon20_NoSpikeFilter_v1")
	hltInfo.HLT_Photon20_NoSpikeFilter_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon20_R9Id_Photon18_R9Id_v1")
	hltInfo.HLT_Photon20_R9Id_Photon18_R9Id_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon20_R9Id_Photon18_R9Id_v2")
	hltInfo.HLT_Photon20_R9Id_Photon18_R9Id_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon20_R9Id_Photon18_R9Id_v3")
	hltInfo.HLT_Photon20_R9Id_Photon18_R9Id_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v1")
	hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v2")
	hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v3")
	hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v1")
	hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v2")
	hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_CaloIdL_IsoVL_Photon18_v1")
	hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_CaloIdL_IsoVL_Photon18_v2")
	hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_CaloIdL_IsoVL_Photon18_v3")
	hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_IsoVL_Photon18_IsoVL_v1")
	hltInfo.HLT_Photon26_IsoVL_Photon18_IsoVL_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_IsoVL_Photon18_IsoVL_v2")
	hltInfo.HLT_Photon26_IsoVL_Photon18_IsoVL_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_IsoVL_Photon18_IsoVL_v3")
	hltInfo.HLT_Photon26_IsoVL_Photon18_IsoVL_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_IsoVL_Photon18_v1")
	hltInfo.HLT_Photon26_IsoVL_Photon18_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_IsoVL_Photon18_v2")
	hltInfo.HLT_Photon26_IsoVL_Photon18_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_IsoVL_Photon18_v3")
	hltInfo.HLT_Photon26_IsoVL_Photon18_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_Photon18_v1")
	hltInfo.HLT_Photon26_Photon18_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_Photon18_v2")
	hltInfo.HLT_Photon26_Photon18_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_Photon18_v3")
	hltInfo.HLT_Photon26_Photon18_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v1")
	hltInfo.HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v2")
	hltInfo.HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon30_CaloIdVL_IsoL_v1")
	hltInfo.HLT_Photon30_CaloIdVL_IsoL_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon30_CaloIdVL_IsoL_v2")
	hltInfo.HLT_Photon30_CaloIdVL_IsoL_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon30_CaloIdVL_IsoL_v3")
	hltInfo.HLT_Photon30_CaloIdVL_IsoL_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon30_CaloIdVL_v1")
	hltInfo.HLT_Photon30_CaloIdVL_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon30_CaloIdVL_v2")
	hltInfo.HLT_Photon30_CaloIdVL_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon30_CaloIdVL_v3")
	hltInfo.HLT_Photon30_CaloIdVL_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1")
	hltInfo.HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon32_CaloIdL_Photon26_CaloIdL_v2")
	hltInfo.HLT_Photon32_CaloIdL_Photon26_CaloIdL_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon32_CaloIdL_Photon26_CaloIdL_v3")
	hltInfo.HLT_Photon32_CaloIdL_Photon26_CaloIdL_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon36_CaloIdL_Photon22_CaloIdL_v1")
	hltInfo.HLT_Photon36_CaloIdL_Photon22_CaloIdL_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon36_CaloIdL_Photon22_CaloIdL_v2")
	hltInfo.HLT_Photon36_CaloIdL_Photon22_CaloIdL_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon50_CaloIdVL_IsoL_v1")
	hltInfo.HLT_Photon50_CaloIdVL_IsoL_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon50_CaloIdVL_IsoL_v2")
	hltInfo.HLT_Photon50_CaloIdVL_IsoL_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon75_CaloIdVL_IsoL_v1")
	hltInfo.HLT_Photon75_CaloIdVL_IsoL_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon75_CaloIdVL_IsoL_v2")
	hltInfo.HLT_Photon75_CaloIdVL_IsoL_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon75_CaloIdVL_IsoL_v3")
	hltInfo.HLT_Photon75_CaloIdVL_IsoL_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon75_CaloIdVL_v1")
	hltInfo.HLT_Photon75_CaloIdVL_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon75_CaloIdVL_v2")
	hltInfo.HLT_Photon75_CaloIdVL_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon75_CaloIdVL_v3") 
	hltInfo.HLT_Photon75_CaloIdVL_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton40_MR150_v3") 
	hltInfo.HLT_DoublePhoton40_MR150_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton40_R014_MR150_v3") 
	hltInfo.HLT_DoublePhoton40_R014_MR150_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton50_v2") 
	hltInfo.HLT_DoublePhoton50_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton5_IsoVL_CEP_v4") 
	hltInfo.HLT_DoublePhoton5_IsoVL_CEP_v4 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton60_v2") 
	hltInfo.HLT_DoublePhoton60_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Mu15_DoublePhoton15_CaloIdL_v6") 
	hltInfo.HLT_Mu15_DoublePhoton15_CaloIdL_v6 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Mu15_Photon20_CaloIdL_v6") 
	hltInfo.HLT_Mu15_Photon20_CaloIdL_v6 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Mu8_Photon20_CaloIdVT_IsoT_v5") 
	hltInfo.HLT_Mu8_Photon20_CaloIdVT_IsoT_v5 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon125_v2") 
	hltInfo.HLT_Photon125_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon200_NoHE_v2") 
	hltInfo.HLT_Photon200_NoHE_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon20_CaloIdVL_IsoL_v4") 
	hltInfo.HLT_Photon20_CaloIdVL_IsoL_v4 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v5") 
	hltInfo.HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v5 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon20_R9Id_Photon18_R9Id_v5") 
	hltInfo.HLT_Photon20_R9Id_Photon18_R9Id_v5 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v5") 
	hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v5 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v4") 
	hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v4 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_CaloIdL_IsoVL_Photon18_v5") 
	hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_v5 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_IsoVL_Photon18_IsoVL_v5") 
	hltInfo.HLT_Photon26_IsoVL_Photon18_IsoVL_v5 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_IsoVL_Photon18_v5") 
	hltInfo.HLT_Photon26_IsoVL_Photon18_v5 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_Photon18_v5") 
	hltInfo.HLT_Photon26_Photon18_v5 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v4") 
	hltInfo.HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v4 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_R9Id_Photon18_R9Id_v2") 
	hltInfo.HLT_Photon26_R9Id_Photon18_R9Id_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon30_CaloIdVL_IsoL_v5") 
	hltInfo.HLT_Photon30_CaloIdVL_IsoL_v5 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon30_CaloIdVL_v5") 
	hltInfo.HLT_Photon30_CaloIdVL_v5 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v1") 
	hltInfo.HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_v1") 
	hltInfo.HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon36_CaloIdL_IsoVL_Photon22_v2") 
	hltInfo.HLT_Photon36_CaloIdL_IsoVL_Photon22_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon36_CaloIdL_Photon22_CaloIdL_v4") 
	hltInfo.HLT_Photon36_CaloIdL_Photon22_CaloIdL_v4 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon36_CaloId_IsoVL_Photon22_R9Id_v1") 
	hltInfo.HLT_Photon36_CaloId_IsoVL_Photon22_R9Id_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon36_IsoVL_Photon22_v2") 
	hltInfo.HLT_Photon36_IsoVL_Photon22_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v1") 
	hltInfo.HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon36_R9Id_Photon22_R9Id_v1") 
	hltInfo.HLT_Photon36_R9Id_Photon22_R9Id_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon40_CaloIdL_Photon28_CaloIdL_v2") 
	hltInfo.HLT_Photon40_CaloIdL_Photon28_CaloIdL_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon40_R005_MR150_v3") 
	hltInfo.HLT_Photon40_R005_MR150_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon40_R014_MR450_v3") 
	hltInfo.HLT_Photon40_R014_MR450_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon40_R020_MR300_v3") 
	hltInfo.HLT_Photon40_R020_MR300_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon40_R025_MR200_v3") 
	hltInfo.HLT_Photon40_R025_MR200_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon40_R038_MR150_v3") 
	hltInfo.HLT_Photon40_R038_MR150_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon50_CaloIdVL_IsoL_v4") 
	hltInfo.HLT_Photon50_CaloIdVL_IsoL_v4 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon50_CaloIdVL_v2") 
	hltInfo.HLT_Photon50_CaloIdVL_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon70_CaloIdL_HT300_v6")
	hltInfo.HLT_Photon70_CaloIdL_HT300_v6 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon70_CaloIdL_HT350_v5") 
	hltInfo.HLT_Photon70_CaloIdL_HT350_v5 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon70_CaloIdL_MHT50_v6") 
	hltInfo.HLT_Photon70_CaloIdL_MHT50_v6 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon70_CaloIdL_MHT70_v5") 
	hltInfo.HLT_Photon70_CaloIdL_MHT70_v5 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon75_CaloIdVL_IsoL_v5") 
	hltInfo.HLT_Photon75_CaloIdVL_IsoL_v5 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon75_CaloIdVL_v5") 
	hltInfo.HLT_Photon75_CaloIdVL_v5 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon90_CaloIdVL_IsoL_v2") 
	hltInfo.HLT_Photon90_CaloIdVL_IsoL_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon90_CaloIdVL_v2") 
	hltInfo.HLT_Photon90_CaloIdVL_v2 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton33_HEVT_v1") hltInfo.HLT_DoublePhoton33_HEVT_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton33_v4") hltInfo.HLT_DoublePhoton33_v4 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton50_v1") hltInfo.HLT_DoublePhoton50_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton5_IsoVL_CEP_v3") hltInfo.HLT_DoublePhoton5_IsoVL_CEP_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton60_v1") hltInfo.HLT_DoublePhoton60_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon125_v1") hltInfo.HLT_Photon125_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon200_NoHE_v1") hltInfo.HLT_Photon200_NoHE_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon20_CaloIdVL_IsoL_v3") hltInfo.HLT_Photon20_CaloIdVL_IsoL_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v4") hltInfo.HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v4 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon20_R9Id_Photon18_R9Id_v4") hltInfo.HLT_Photon20_R9Id_Photon18_R9Id_v4 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v4") hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v4 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v3") hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_CaloIdL_IsoVL_Photon18_v4") hltInfo.HLT_Photon26_CaloIdL_IsoVL_Photon18_v4 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_IsoVL_Photon18_IsoVL_v4") hltInfo.HLT_Photon26_IsoVL_Photon18_IsoVL_v4 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_IsoVL_Photon18_v4") hltInfo.HLT_Photon26_IsoVL_Photon18_v4 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_Photon18_v4") hltInfo.HLT_Photon26_Photon18_v4 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v3") hltInfo.HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon26_R9Id_Photon18_R9Id_v1") hltInfo.HLT_Photon26_R9Id_Photon18_R9Id_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon30_CaloIdVL_IsoL_v4") hltInfo.HLT_Photon30_CaloIdVL_IsoL_v4 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon30_CaloIdVL_v4") hltInfo.HLT_Photon30_CaloIdVL_v4 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon32_CaloIdL_Photon26_CaloIdL_v4") hltInfo.HLT_Photon32_CaloIdL_Photon26_CaloIdL_v4 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon36_CaloIdL_IsoVL_Photon22_v1") hltInfo.HLT_Photon36_CaloIdL_IsoVL_Photon22_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon36_CaloIdL_Photon22_CaloIdL_v3") hltInfo.HLT_Photon36_CaloIdL_Photon22_CaloIdL_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon36_IsoVL_Photon22_v1") hltInfo.HLT_Photon36_IsoVL_Photon22_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon40_CaloIdL_Photon28_CaloIdL_v1") hltInfo.HLT_Photon40_CaloIdL_Photon28_CaloIdL_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon50_CaloIdVL_IsoL_v3") hltInfo.HLT_Photon50_CaloIdVL_IsoL_v3 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon50_CaloIdVL_v1") hltInfo.HLT_Photon50_CaloIdVL_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon75_CaloIdVL_IsoL_v4") hltInfo.HLT_Photon75_CaloIdVL_IsoL_v4 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon75_CaloIdVL_v4") hltInfo.HLT_Photon75_CaloIdVL_v4 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon90_CaloIdVL_IsoL_v1") hltInfo.HLT_Photon90_CaloIdVL_IsoL_v1 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon90_CaloIdVL_v1") hltInfo.HLT_Photon90_CaloIdVL_v1 = (int) hltResults->accept(itrig);


    }

  }


}  //end of namespace


#endif
