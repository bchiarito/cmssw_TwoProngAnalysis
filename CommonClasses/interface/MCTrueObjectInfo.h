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
#include <vector>

// for MC
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "TMath.h"

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
    double isol04;
    double isol04ratio;
    double isol03;
    double isol03ratio;
    double isol02;
    double isol02ratio;
    //store deltaR of the match?
  };

  // string to define the tree branch
  std::string mcTrueObjectInfoBranchDefString("status/I:PdgId:MotherPdgId:GrandmotherPdgId:pt/D:eta/D:phi/D:isol04/D:isol04ratio/D:isol03/D:isol03ratio/D:isol02/D:isol02ratio/D");

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
  
  // Same FillMCTrueObjectInfo method but with reco::Candidates
  void FillMCTrueObjectInfo(mcTrueObjectInfo_t &mcTrueObjectInfo, const reco::Candidate *genParticle) {
    
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
  
  // also want to store MC truth event-level info
  // like signalProcess ID, and pthat value

  struct mcEventInfo_t {
    double binningValue;
    int SignalProcessId;
  };

  std::string mcEventInfoBranchDefString("binningValue/D:signalProcessId/I");

  void FillMCEventInfo(mcEventInfo_t &mcEventInfo, const GenEventInfoProduct *genEventInfo) {
    
    mcEventInfo.SignalProcessId = genEventInfo->signalProcessID();

    // get binning values 
    // Fabian says it's an array because in principle 
    // samples could be binned in more than 1 variable
    // generally this should be pthat value
    // for s^hat binned samples, will this be s^hat?
    const std::vector<double> genBinningValues = genEventInfo->binningValues();
    
    mcEventInfo.binningValue = genBinningValues[0];
  }


  bool compareGenPhotonsByPt(reco::GenParticle genphoton1, reco::GenParticle genphoton2) {return(genphoton1.pt() >= genphoton2.pt());}
  bool compareGenPhotonPointersByPt(const reco::GenParticle *genphoton1, const reco::GenParticle *genphoton2) {return(genphoton1->pt() >= genphoton2->pt());}

  //bool compareGenPhotonPairsByPt(std::pair<reco::GenParticle, float> genpair1, std::pair<reco::GenParticle, float> genpair2) {return(genpair1.first.pt() >= genpair2.first.pt());}
  //bool compareGenPhotonPairsByPt(std::pair<const reco::GenParticle*, const reco::Candidate*> genpair1, std::pair<const reco::GenParticle*, const reco::Candidate*> genpair2) {return(genpair1.first->pt() >= genpair2.first->pt());}


  bool compareGenPhotonPairsByPt(std::pair<reco::GenParticle, const reco::Candidate*> genpair1, std::pair<reco::GenParticle, const reco::Candidate*> genpair2) {return(genpair1.first.pt() >= genpair2.first.pt());}

  /*   bool compareGenPhotonPairsByDeltaR(std::pair<reco::GenParticle, const reco::Candidate*> genpair1, std::pair<reco::GenParticle, const reco::Candidate*> genpair2)  */
  /*   { */
  /*     float deltar1 = deltaR(genpair1.first.eta(),genpair1.first.phi(),genpair1.second->eta(),genpair1.second->phi()); */
  /*     float deltar2 = deltaR(genpair2.first.eta(),genpair2.first.phi(),genpair2.second->eta(),genpair2.second->phi()); */
  /*     return(deltar1 < deltar2); */
  /*   } */

}// end of namespace

#endif

