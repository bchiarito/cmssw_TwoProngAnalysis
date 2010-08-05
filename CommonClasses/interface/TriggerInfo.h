#ifndef TRIG_INFO_INC
#define TRIG_INFO_INC

//********************************************************************
// Definition of a struct that can be used for storing trig info
// in a tree, from different analysers
// Also includes a Fill function to fill the struct from the appropriate objects
// and a string that can be used to define the tree branch
// 
// Conor, July 2010
// 
//********************************************************************

#include <string>

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



  // HLT info - MinBias trigger?  Definitely photon triggers!
  // this will need to be kept up to date since the HLT triggers evolve
  // eg with spike cleaning, prescales, etc
  // store prescale value too?
  struct hltTrigInfo_t{
    bool HLT_MinBiasBSC;
    bool HLT_MinBiasBSC_NoBPTX;
    bool HLT_MinBiasBSC_OR;
    bool HLT_L1_BscMinBiasOR_BptxPlusORMinus;
    bool HLT_Jet15U;
    bool HLT_Jet30U;
    bool HLT_Jet50U;
    bool HLT_L1SingleEG2;
    bool HLT_L1SingleEG2_NoBPTX;
    bool HLT_L1SingleEG5;
    bool HLT_L1SingleEG5_NoBPTX;
    bool HLT_L1SingleEG8;
    bool HLT_L1SingleEG20_NoBPTX;
    bool HLT_L1DoubleEG5;
    bool HLT_EgammaSuperClusterOnly_L1R;
    bool HLT_Photon10_L1R;
    bool HLT_Photon15_L1R;
    bool HLT_Photon15_TrackIso_L1R;
    bool HLT_Photon15_LooseEcalIso_L1R;
    bool HLT_Photon20_L1R;
    bool HLT_Photon30_L1R_8E29;
    bool HLT_DoublePhoton5_L1R;
    bool HLT_DoublePhoton10_L1R;

  };

  // string for defining a tree branch for this struct
  std::string hltTrigBranchDefString(":HLT_MinBiasBSC/O:HLT_MinBiasBSC_NoBPTX:HLT_MinBiasBSC_OR:HLT_L1_BscMinBiasOR_BptxPlusORMinus:HLT_Jet15U/O:HLT_Jet30U:HLT_Jet50U:HLT_L1SingleEG2:HLT_L1SingleEG2_NoBPTX:HLT_L1SingleEG5:HLT_L1SingleEG5_NoBPTX:HLT_L1SingleEG8:HLT_L1SingleEG20_NoBPTX:HLT_L1DoubleEG5:HLT_EgammaSuperClusterOnly_L1R:HLT_Photon10_L1R:HLT_Photon15_L1R:HLT_Photon15_TrackIso_L1R:HLT_Photon15_LooseEcalIso_L1R:HLT_Photon20_L1R:HLT_Photon30_L1R_8E29:HLT_DoublePhoton5_L1R:HLT_DoublePhoton10_L1R");
  
  // need also to have an Initialise() function?

  // Fill() function for HLT info to appear here ....




}  //end of namespace


#endif
