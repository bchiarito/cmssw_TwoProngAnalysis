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

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

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
    int HLT_Photon20_L1R;
    int HLT_Photon20_Cleaned_L1R;
    int HLT_Photon30_L1R;
    int HLT_Photon30_Cleaned_L1R;
    int HLT_Photon30_L1R_8E29;
    int HLT_Photon50_L1R;
    int HLT_Photon50_Cleaned_L1R;
    int HLT_DoublePhoton5_L1R;
    int HLT_DoublePhoton10_L1R;
    int HLT_DoublePhoton15_L1R;
    int HLT_DoublePhoton20_L1R;
  };

  // string for defining a tree branch for this struct
  std::string hltTrigBranchDefString("HLT_MinBiasBSC/I:HLT_MinBiasBSC_NoBPTX:HLT_MinBiasBSC_OR:HLT_L1_BscMinBiasOR_BptxPlusORMinus:HLT_L1SingleEG2:HLT_L1SingleEG5:HLT_L1SingleEG8:HLT_L1DoubleEG5:HLT_Photon10_L1R:HLT_Photon10_Cleaned_L1R:HLT_Photon15_L1R:HLT_Photon15_Cleaned_L1R:HLT_Photon15_LooseEcalIso_L1R:HLT_Photon15_LooseEcalIso_Cleaned_L1R:HLT_Photon15_TrackIso_L1R:HLT_Photon15_TrackIso_Cleaned_L1R:HLT_Photon20_L1R:HLT_Photon20_Cleaned_L1R:HLT_Photon30_L1R:HLT_Photon30_Cleaned_L1R:HLT_Photon30_L1R_8E29:HLT_Photon50_L1R:HLT_Photon50_Cleaned_L1R:HLT_DoublePhoton5_L1R:HLT_DoublePhoton10_L1R:HLT_DoublePhoton15_L1R:HLT_DoublePhoton20_L1R");


  
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
    hltInfo.HLT_Photon10_Cleaned_L1R = -1;
    hltInfo.HLT_Photon15_Cleaned_L1R = -1;
    hltInfo.HLT_Photon15_L1R = -1;
    hltInfo.HLT_Photon15_LooseEcalIso_L1R = -1;
    hltInfo.HLT_Photon15_TrackIso_L1R = -1;
    hltInfo.HLT_Photon20_Cleaned_L1R = -1;
    hltInfo.HLT_Photon20_L1R = -1;
    hltInfo.HLT_Photon30_Cleaned_L1R = -1;
    hltInfo.HLT_Photon30_L1R_8E29 = -1;
    hltInfo.HLT_Photon50_Cleaned_L1R = -1;
    hltInfo.HLT_Photon50_L1R = -1;
    hltInfo.HLT_DoublePhoton5_L1R = -1;
    hltInfo.HLT_DoublePhoton10_L1R = -1;
    hltInfo.HLT_DoublePhoton15_L1R = -1;
    hltInfo.HLT_DoublePhoton20_L1R = -1;
    hltInfo.HLT_Photon10_L1R = -1;
    hltInfo.HLT_Photon15_LooseEcalIso_Cleaned_L1R = -1;
    hltInfo.HLT_Photon15_TrackIso_Cleaned_L1R = -1;
    hltInfo.HLT_Photon30_L1R = -1;

    // we'll just loop over all triggers in the current event
    // and check one-by-one the result for the triggers we are interested in

    for(int itrig=0;itrig<hltResults->size();itrig++) {

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
      else if(hltNames.triggerName(itrig)=="HLT_Photon10_Cleaned_L1R")
	hltInfo.HLT_Photon10_Cleaned_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon15_Cleaned_L1R")
	hltInfo.HLT_Photon15_Cleaned_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon15_L1R")
	hltInfo.HLT_Photon15_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon15_LooseEcalIso_L1R")
	hltInfo.HLT_Photon15_LooseEcalIso_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon15_TrackIso_L1R")
	hltInfo.HLT_Photon15_TrackIso_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon20_Cleaned_L1R")
	hltInfo.HLT_Photon20_Cleaned_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon20_L1R")
	hltInfo.HLT_Photon20_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon30_Cleaned_L1R")
	hltInfo.HLT_Photon30_Cleaned_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon30_L1R_8E29")
	hltInfo.HLT_Photon30_L1R_8E29 = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon50_Cleaned_L1R")
	hltInfo.HLT_Photon50_Cleaned_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon50_L1R")
	hltInfo.HLT_Photon50_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton5_L1R")
	hltInfo.HLT_DoublePhoton5_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton10_L1R")
	hltInfo.HLT_DoublePhoton10_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton15_L1R")
	hltInfo.HLT_DoublePhoton15_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_DoublePhoton20_L1R")
	hltInfo.HLT_DoublePhoton20_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon10_L1R")
        hltInfo.HLT_Photon10_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon15_LooseEcalIso_Cleaned_L1R")
	hltInfo.HLT_Photon15_LooseEcalIso_Cleaned_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon15_TrackIso_Cleaned_L1R")
	hltInfo.HLT_Photon15_TrackIso_Cleaned_L1R = (int) hltResults->accept(itrig);
      else if(hltNames.triggerName(itrig)=="HLT_Photon30_L1R")
	hltInfo.HLT_Photon30_L1R = (int) hltResults->accept(itrig);
    }


  }


}  //end of namespace


#endif
