#ifndef MC_TRUE_OBJECT_INFO
#define MC_TRUE_OBJECT_INFO

//********************************************************************
// Definition of a struct that can be used for storing MC truth info
// for an object in a tree, from different analysers
// Also includes a Fill function to fill the struct from the appropriate objects
// and a string that can be used to define the tree branch
// 
// Conor, July 2010
// 
//********************************************************************

#include <string>

namespace ExoDiPhotons
{
  struct mcTrueObjectInfo_t {
    int status; // will be 3 or 1 
    int PdgId;
    int MotherPdgId;
    int GrandmotherPdgId; 
    double pt;
    double eta;
    double phi;
    //store deltaR of the match?
  };

  // string to define the tree branch
  std::string mcTrueObjectInfoBranchDefString("status/I:PdgId:MotherPdgId:GrandmotherPdgId:pt/D:eta/D:phi/D");

  // convenient function to fill the struct (pass by reference) from the desired object
  void FillMCTrueObjectInfo(mcTrueObjectInfo_t &mcTrueObjectInfo, const reco::GenParticle *genParticle) {
    
    mcTrueObjectInfo.status = genParticle->status();
    mcTrueObjectInfo.PdgId = genParticle->pdgId();

    // careful here when trying to get pointers to mother/grandmother
    if(genParticle->numberOfMothers()>0) {
      mcTrueObjectInfo.MotherPdgId = genParticle->mother()->pdgId();
      if(genParticle->mother()->numberOfMothers()>0) {
	mcTrueObjectInfo.GrandmotherPdgId = genParticle->mother()->mother()->pdgId();
      }
      else {
	mcTrueObjectInfo.GrandmotherPdgId = -9999999;
      }
    }
    else {
      mcTrueObjectInfo.MotherPdgId = -9999999;
      mcTrueObjectInfo.GrandmotherPdgId = -9999999;
    }



    mcTrueObjectInfo.pt = genParticle->pt();
    mcTrueObjectInfo.eta = genParticle->eta();
    mcTrueObjectInfo.phi = genParticle->phi();

  }
  





}// end of namespace

#endif

