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
    int status;
    int motherStatus;
    int grandmotherStatus;
    int pdgId;
    int motherPdgId;
    int grandmotherPdgId;
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
  std::string mcTrueObjectInfoBranchDefString("status/I:motherStatus:grandmotherStatus:pdgId:motherPdgId:grandmotherPdgId:pt/D:eta:phi:isol04:isol04ratio:isol03:isol03ratio:isol02:isol02ratio");

  // convenient function to fill the struct (pass by reference) from the desired object
  void FillMCTrueObjectInfo(mcTrueObjectInfo_t &mcTrueObjectInfo, const reco::GenParticle *genParticle) {

    double minMotherDeltaR = 100000; // consider all mothers
    double minGrandmotherDeltaR = 100000;
    int motherIndex = 0;
    int grandmotherIndex = 0;
    
    mcTrueObjectInfo.status = genParticle->status();
    mcTrueObjectInfo.pdgId = genParticle->pdgId();

    // careful here when trying to get pointers to mother/grandmother
    if (genParticle->numberOfMothers() > 0) {
      // find best match in deltaR among all mothers
      for (unsigned int j = 0; j < genParticle->numberOfMothers(); j++) {
	double deltaR = reco::deltaR(genParticle->eta(),genParticle->phi(),genParticle->mother(j)->eta(),genParticle->mother(j)->phi());
	if (deltaR < minMotherDeltaR) {
	  minMotherDeltaR = deltaR;
	  motherIndex = j;
	}
      }
      mcTrueObjectInfo.motherPdgId = genParticle->mother(motherIndex)->status();
      mcTrueObjectInfo.motherPdgId = genParticle->mother(motherIndex)->pdgId();
      
      if (genParticle->mother(motherIndex)->numberOfMothers() > 0) {
	// find best match in deltaR among all mothers
	for (unsigned int j = 0; j < genParticle->mother(motherIndex)->numberOfMothers(); j++) {
	  double deltaR = reco::deltaR(genParticle->mother(motherIndex)->eta(),genParticle->mother(motherIndex)->phi(),
				       genParticle->mother(motherIndex)->mother(j)->eta(),genParticle->mother(motherIndex)->mother(j)->phi());
	  if (deltaR < minGrandmotherDeltaR) {
	    minGrandmotherDeltaR = deltaR;
	    grandmotherIndex = j;
	  }
	}
	mcTrueObjectInfo.grandmotherPdgId = genParticle->mother(motherIndex)->mother(grandmotherIndex)->status();
	mcTrueObjectInfo.grandmotherPdgId = genParticle->mother(motherIndex)->mother(grandmotherIndex)->pdgId();
      }
      else {
	mcTrueObjectInfo.grandmotherStatus = -9999999;
	mcTrueObjectInfo.grandmotherPdgId = -9999999;
      }
    }
    else {
      mcTrueObjectInfo.motherStatus = -9999999;
      mcTrueObjectInfo.grandmotherStatus = -9999999;
      mcTrueObjectInfo.motherPdgId = -9999999;
      mcTrueObjectInfo.grandmotherPdgId = -9999999;
    }
    
    mcTrueObjectInfo.pt = genParticle->pt();
    mcTrueObjectInfo.eta = genParticle->eta();
    mcTrueObjectInfo.phi = genParticle->phi();
  }
  
  // Same FillMCTrueObjectInfo method but with reco::Candidates
  void FillMCTrueObjectInfo(mcTrueObjectInfo_t &mcTrueObjectInfo, const reco::Candidate *genParticle) {

    double minMotherDeltaR = 100000; // consider all mothers
    double minGrandmotherDeltaR = 100000;
    int motherIndex = 0;
    int grandmotherIndex = 0;
    
    mcTrueObjectInfo.status = genParticle->status();
    mcTrueObjectInfo.pdgId = genParticle->pdgId();

    // careful here when trying to get pointers to mother/grandmother
    if (genParticle->numberOfMothers() > 0) {
      // find best match in deltaR among all mothers
      for (unsigned int j = 0; j < genParticle->numberOfMothers(); j++) {
	double deltaR = reco::deltaR(genParticle->eta(),genParticle->phi(),genParticle->mother(j)->eta(),genParticle->mother(j)->phi());
	if (deltaR < minMotherDeltaR) {
	  minMotherDeltaR = deltaR;
	  motherIndex = j;
	}
      }
      mcTrueObjectInfo.motherPdgId = genParticle->mother(motherIndex)->status();
      mcTrueObjectInfo.motherPdgId = genParticle->mother(motherIndex)->pdgId();

      if (genParticle->mother(motherIndex)->numberOfMothers() > 0) {
	// find best match in deltaR among all mothers
	for (unsigned int j = 0; j < genParticle->mother(motherIndex)->numberOfMothers(); j++) {
	  double deltaR = reco::deltaR(genParticle->mother(motherIndex)->eta(),genParticle->mother(motherIndex)->phi(),
				       genParticle->mother(motherIndex)->mother(j)->eta(),genParticle->mother(motherIndex)->mother(j)->phi());
	  if (deltaR < minGrandmotherDeltaR) {
	    minGrandmotherDeltaR = deltaR;
	    grandmotherIndex = j;
	  }
	}
	mcTrueObjectInfo.grandmotherStatus = genParticle->mother(motherIndex)->mother(grandmotherIndex)->status();
	mcTrueObjectInfo.grandmotherPdgId = genParticle->mother(motherIndex)->mother(grandmotherIndex)->pdgId();
      }
      else {
	mcTrueObjectInfo.grandmotherStatus = -9999999;
	mcTrueObjectInfo.grandmotherPdgId = -9999999;
      }
    }
    else {
      mcTrueObjectInfo.motherStatus = -9999999;
      mcTrueObjectInfo.grandmotherStatus = -9999999;
      mcTrueObjectInfo.motherPdgId = -9999999;
      mcTrueObjectInfo.grandmotherPdgId = -9999999;
    }
    
    mcTrueObjectInfo.pt = genParticle->pt();
    mcTrueObjectInfo.eta = genParticle->eta();
    mcTrueObjectInfo.phi = genParticle->phi();
  }
  
  // also want to store MC truth event-level info
  // like signalProcess ID, and pthat value

  void InitMCTrueObjectInfo(mcTrueObjectInfo_t &mcTrueObjectInfo) {
    mcTrueObjectInfo.status = -999999;
    mcTrueObjectInfo.motherStatus = -999999;
    mcTrueObjectInfo.grandmotherStatus = -999999;
    mcTrueObjectInfo.pdgId = -999999;
    mcTrueObjectInfo.motherPdgId = -999999;
    mcTrueObjectInfo.grandmotherPdgId = -999999;
    mcTrueObjectInfo.pt = -999999.99;
    mcTrueObjectInfo.eta = -999999.99;
    mcTrueObjectInfo.phi = -999999.99;
    mcTrueObjectInfo.isol04 = -999999.99;
    mcTrueObjectInfo.isol04ratio = -999999.99;
    mcTrueObjectInfo.isol03 = -999999.99;
    mcTrueObjectInfo.isol03ratio = -999999.99;
    mcTrueObjectInfo.isol02 = -999999.99;
    mcTrueObjectInfo.isol02ratio = -999999.99;
  }



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

