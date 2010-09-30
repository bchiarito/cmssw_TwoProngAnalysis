#ifndef PHOTON_ID_INC
#define PHOTON_ID_INC

//for photons
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"


namespace ExoDiPhotons{


  // now can apply each variable cut separately

  bool passesHadOverEmCut(const reco::Photon *photon) {
    
    bool result = false;

    if(photon->hadronicOverEm()<0.05)
      result = true;

    return result;
  }


  bool passesTrkIsoCut(const reco::Photon *photon) {
    
    bool result = false;

    double trkIsoCut = 3.5 + 0.001*photon->et();
    if(photon->trkSumPtHollowConeDR04()<trkIsoCut)
      result = true;

    return result;
  }


  bool passesEcalIsoCut(const reco::Photon *photon) {
    
    bool result = false;

    double ecalIsoCut = 4.2 + 0.006*photon->et();
    if(photon->ecalRecHitSumEtConeDR04()<ecalIsoCut)
      result = true;

    return result;
  }


  bool passesHcalIsoCut(const reco::Photon *photon) {
    
    bool result = false;

    double hcalIsoCut = 2.2 + 0.0025*photon->et();
    if(photon->hcalTowerSumEtConeDR04()<hcalIsoCut) 
      result = true;

    return result;
  }


  bool passesSigmaIetaIetaCut(const reco::Photon *photon) {
    
    bool result = false;

    // sigma ieta ieta cuts different in barrel and endcap
    if(photon->isEB()) {
      if(photon->sigmaIetaIeta()< 0.013)
	result = true;
      else
	result = false;
    }
    else if(photon->isEE()) {
      if(photon->sigmaIetaIeta()< 0.030)
	result = true;
      else
	result = false;
    }

    return result;
  }


  bool passesNoPixelSeedCut(const reco::Photon *photon) {
    
    bool result = false;

    if(photon->hasPixelSeed()==false)
      result = true;

    return result;
  }



  bool isTightPhoton(const reco::Photon *photon) {

    bool result = false;

    if(passesHadOverEmCut(photon) && passesTrkIsoCut(photon) && passesEcalIsoCut(photon) && passesHcalIsoCut(photon) && passesSigmaIetaIetaCut(photon) && passesNoPixelSeedCut(photon))
      result = true;

    return result; 

  }

  // for offline spike flagging in our analysers

  bool isASpike(const reco::Photon *photon) {

    bool thisPhotonIsASpike = false;

    // can do spike identification based on regular swiss cross
    // or on out-of-time
    // even add double-spike determination in here too

    // can be used on MC also
    
    return thisPhotonIsASpike;
  }

}


#endif
