#ifndef PHOTON_ID_INC
#define PHOTON_ID_INC

//for photons
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"


namespace ExoDiPhotons{

  bool isTightPhoton(const reco::Photon *photon) {

    bool result = false;

    // these cuts are just hardcoded for now ...

    bool hadOverEmResult = false;
    bool trkIsoResult =false;
    bool hcalIsoResult = false;
    bool ecalIsoResult = false;
    bool noPixelSeedResult = false;
    bool sigmaIetaIetaResult = false;
  
    if(photon->hadronicOverEm()<0.05)
      hadOverEmResult = true;


    double trkIsoCut = 2.0 + 0.001*photon->et();
    if(photon->trkSumPtHollowConeDR04()<trkIsoCut)
      trkIsoResult = true;

    double hcalIsoCut = 2.2 + 0.001*photon->et();
    if(photon->hcalTowerSumEtConeDR04()<hcalIsoCut) 
      hcalIsoResult = true;

    double ecalIsoCut = 4.2 + 0.003*photon->et();
    if(photon->ecalRecHitSumEtConeDR04()<ecalIsoCut)
      ecalIsoResult = true;

    if(photon->hasPixelSeed()==false)
      noPixelSeedResult = true; // ie it is true that it does NOT have a pixel seed!

    // sigmaIetaIeta should be included in tight photon ID too soon ...

    if(hadOverEmResult && trkIsoResult && ecalIsoResult && hcalIsoResult && noPixelSeedResult) 
      result = true;

    return result; 

  }

}


#endif
