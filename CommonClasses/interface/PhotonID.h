#ifndef PHOTON_ID_INC
#define PHOTON_ID_INC

//for photons
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "TMath.h"

namespace ExoDiPhotons{


  // now can apply each variable cut separately

  bool passesHadOverEmCut(const reco::Photon *photon) {
    
    bool result = false;

    if(photon->hadronicOverEm()<0.05)
      result = true;

    return result;
  }


  bool passesTrkIsoCut(const reco::Photon *photon,double rho25) {
    
    bool result = false;

    //    double trkIsoCut = 3.5 + 0.001*photon->et();
    //    double trkIsoCut = 2.0 + 0.001*photon->et();
    double trkIsoCut = 2.0+0.001*photon->et()+0.0167*rho25; 
    if(photon->trkSumPtHollowConeDR04()<trkIsoCut)
      result = true;

    return result;
  }


  bool passesEcalIsoCut(const reco::Photon *photon,double rho25) {
    
    bool result = false;

    //double ecalIsoCut = 4.2 + 0.006*photon->et();
    // double ecalIsoCut = 4.2 + 0.006*photon->et()+0.090*rho;
    double ecalIsoCut = 4.2 + 0.006*photon->et()+0.183*rho25;
    if(photon->ecalRecHitSumEtConeDR04()<ecalIsoCut)
      result = true;

    return result;
  }


  bool passesHcalIsoCut(const reco::Photon *photon,double rho25) {
    
    bool result = false;

    //   double hcalIsoCut = 2.2 + 0.0025*photon->et();
    double hcalIsoCut = 2.2 + 0.0025*photon->et()+0.062*rho25;
    if(photon->hcalTowerSumEtConeDR04()<hcalIsoCut) 
      result = true;

    return result;
  }


  bool passesSigmaIetaIetaCut(const reco::Photon *photon) {
    
    bool result = false;

    // sigma ieta ieta cuts different in barrel and endcap
    if(photon->isEB()) {
      if(photon->sigmaIetaIeta()< 0.011)
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
    //   return true;//keep Z->ee "photons" special

  }



  bool isTightPhoton(const reco::Photon *photon,double rho25) {

    bool result = false;

    if(passesHadOverEmCut(photon) && passesTrkIsoCut(photon,rho25) && passesEcalIsoCut(photon,rho25) && passesHcalIsoCut(photon,rho25) && passesSigmaIetaIetaCut(photon) && passesNoPixelSeedCut(photon))
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

  bool isBarrelPhoton(const reco::Photon *photon) {

    if ( fabs(photon->caloPosition().eta())<=1.4442)
      return true;
    else
      return false;

  }


  bool isGapPhoton(const reco::Photon *photon) {

    if ( fabs(photon->caloPosition().eta())>1.4442 && fabs(photon->caloPosition().eta())<1.566 ) return true; 
    else return false;
  }


  //bool isFakeableObject(const reco::Photon *photon) {
  bool isFakeableObject(const reco::Photon *photon,double rho25) {

    bool result = false;
    
    // for testing purposes:
    //    return true;

    double pt = photon->pt();

    double thisTrkIso = photon->trkSumPtHollowConeDR04();
    double thisEcalIso = photon->ecalRecHitSumEtConeDR04();
    double thisHcalIso = photon->hcalTowerSumEtConeDR04();

    // define fakeable objects as denominator for fake rate studies

    // we are defining them to be loose, but also exclude real photons (to large degree) 

    // these 'swing values' determine both loose limit and exclusion veto
    //double trkIsoSwingValue = 3.5 + 0.001*pt;
    double trkIsoSwingValue = 2.0 + 0.001*pt + 0.0167*rho25;
    double trkIsoLooseLimit = TMath::Min( 5.0*(trkIsoSwingValue), 0.2*pt );
    double trkIsoExclusion = trkIsoSwingValue;
    // for testing, to see if I can get a Fakeable event in small sample
    //    double trkIsoExclusion = 2.0 + 0.001*pt;


    //double ecalIsoSwingValue = 4.2 + 0.006*pt;
    double ecalIsoSwingValue = 4.2 + 0.006*pt + 0.183*rho25;;
    double ecalIsoLooseLimit = TMath::Min( 5.0*(ecalIsoSwingValue), 0.2*pt );
    double ecalIsoExclusion = ecalIsoSwingValue;

    //double hcalIsoSwingValue = 2.2 + 0.0025*pt;
    double hcalIsoSwingValue = 2.2 + 0.0025*pt + 0.062*rho25;
    double hcalIsoLooseLimit = TMath::Min( 5.0*(hcalIsoSwingValue), 0.2*pt );
    double hcalIsoExclusion = hcalIsoSwingValue;
    
    if( photon->hadronicOverEm()<0.05 && 
	thisTrkIso < trkIsoLooseLimit &&
	thisEcalIso < ecalIsoLooseLimit &&
	thisHcalIso < hcalIsoLooseLimit  && 
	( thisTrkIso > trkIsoExclusion || thisEcalIso > ecalIsoExclusion || thisHcalIso > hcalIsoExclusion || !passesSigmaIetaIetaCut(photon) )
	) {
      // then this passes our fakeable definition
      result = true;
    }
   
    return result;
  }

}


#endif
