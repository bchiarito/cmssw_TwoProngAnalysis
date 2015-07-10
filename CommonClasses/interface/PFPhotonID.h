#ifndef PFPHOTON_ID_INC
#define PFPHOTON_ID_INC

//for photons
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "TMath.h"

namespace ExoDiPhotons{


  //For the moment, the PF isolation variables cannot be accessed 
  //directly from the Photon class
  //Neither is the conversion safe electron veto
  //However the SingleTower H/E is
  //So we'll just provide the corresponding isolation values to each method

  //Now all of this will be for the medium ID (80%)
  //and for the moment only for the barrel

  //All PF isolation variables provided will already 
  //be rho corrected. So rho25 is useless as argument

  // now can apply each variable cut separately

  bool passesHadTowerOverEmCut(const reco::Photon *photon,TString CategoryID) {
    
    bool result = false;

    double hadTowerOverEmCut = 0.;

    if(photon->isEB()){
      if(CategoryID.Contains("Tight")) hadTowerOverEmCut = 0.05;
      if(CategoryID.Contains("Medium")) hadTowerOverEmCut = 0.05;
      if(CategoryID.Contains("Loose")) hadTowerOverEmCut = 0.05;
    }

    if(photon->isEE()){
      if(CategoryID.Contains("Tight")) hadTowerOverEmCut = 0.05;
      if(CategoryID.Contains("Medium")) hadTowerOverEmCut = 0.05;
      if(CategoryID.Contains("Loose")) hadTowerOverEmCut = 0.05;
    }

    if(photon->hadTowOverEm()<hadTowerOverEmCut)
      result = true;

    return result;
  }


  bool passesChargedHadronCut(const reco::Photon *photon,double rhocorPFChargedHadronIso,TString CategoryID) {
    
    bool result = false;

    std::pair<double,double> CHIsoCut (0.,0.);

    if(photon->isEB()){
      if(CategoryID.Contains("Tight")) {CHIsoCut.first = 0.7;CHIsoCut.second=0.;}
      if(CategoryID.Contains("Medium")) {CHIsoCut.first = 1.5;CHIsoCut.second=0.;}
      if(CategoryID.Contains("Loose")) {CHIsoCut.first = 2.6;CHIsoCut.second=0.;}
    }

    if(photon->isEE()){
      if(CategoryID.Contains("Tight")) {CHIsoCut.first = 0.5;CHIsoCut.second=0.;}
      if(CategoryID.Contains("Medium")) {CHIsoCut.first = 1.2;CHIsoCut.second=0.;}
      if(CategoryID.Contains("Loose")) {CHIsoCut.first = 2.3;CHIsoCut.second=0.;}
    }

    if( rhocorPFChargedHadronIso < (CHIsoCut.first + CHIsoCut.second * photon->et()))
      result = true;

    return result;
  }


  bool passesPhotonCut(const reco::Photon *photon,double rhocorPFPhotonIso,TString CategoryID) {
    
    bool result = false;

    std::pair<double,double> PHIsoCut (0.,0.);

    if(photon->isEB()){
      if(CategoryID.Contains("Tight")) {PHIsoCut.first = 0.5;PHIsoCut.second=0.005;}
      if(CategoryID.Contains("Medium")) {PHIsoCut.first = 0.7;PHIsoCut.second=0.005;}
      if(CategoryID.Contains("Loose")) {PHIsoCut.first = 1.3;PHIsoCut.second=0.005;}
    }

    if(photon->isEE()){
      if(CategoryID.Contains("Tight")) {PHIsoCut.first = 1.0;PHIsoCut.second=0.005;}
      if(CategoryID.Contains("Medium")) {PHIsoCut.first = 1.0;PHIsoCut.second=0.005;}
      if(CategoryID.Contains("Loose")) {PHIsoCut.first = 99999.;PHIsoCut.second=99999.;}
    }

    if( rhocorPFPhotonIso < (PHIsoCut.first + PHIsoCut.second * photon->et()))
      result = true;

    return result;
  }


  bool passesNeutralHadronCut(const reco::Photon *photon,double rhocorPFNeutralHadronIso,TString CategoryID) {
    
    bool result = false;

    std::pair<double,double> NHIsoCut (0.,0.);

    if(photon->isEB()){
      if(CategoryID.Contains("Tight")) {NHIsoCut.first = 0.4;NHIsoCut.second=0.04;}
      if(CategoryID.Contains("Medium")) {NHIsoCut.first = 1.0;NHIsoCut.second=0.04;}
      if(CategoryID.Contains("Loose")) {NHIsoCut.first = 3.5;NHIsoCut.second=0.04;}
    }

    if(photon->isEE()){
      if(CategoryID.Contains("Tight")) {NHIsoCut.first = 1.5;NHIsoCut.second=0.04;}
      if(CategoryID.Contains("Medium")) {NHIsoCut.first = 1.5;NHIsoCut.second=0.04;}
      if(CategoryID.Contains("Loose")) {NHIsoCut.first = 2.9;NHIsoCut.second=0.04;}
    }

    if( rhocorPFNeutralHadronIso < (NHIsoCut.first + NHIsoCut.second * photon->et()))
      result = true;

    return result;
  }


  bool passesPFSigmaIetaIetaCut(const reco::Photon *photon,TString CategoryID) {
    
    bool result = false;

    double sigmaIetaIetaCut = 0.;

    if(photon->isEB()){
      if(CategoryID.Contains("Tight")) sigmaIetaIetaCut = 0.011;
      if(CategoryID.Contains("Medium")) sigmaIetaIetaCut = 0.011;
      if(CategoryID.Contains("Loose")) sigmaIetaIetaCut = 0.012;
    }

    if(photon->isEE()){
      if(CategoryID.Contains("Tight")) sigmaIetaIetaCut = 0.031;
      if(CategoryID.Contains("Medium")) sigmaIetaIetaCut = 0.033;
      if(CategoryID.Contains("Loose")) sigmaIetaIetaCut = 0.034;
    }

    if(photon->sigmaIetaIeta()<sigmaIetaIetaCut)
      result = true;

    return result;
  }



  bool isPFTightPhoton(const reco::Photon *photon,double rhocorPFChargedHadronIso,double rhocorPFNeutralHadronIso,double rhocorPFPhotonIso,bool passesElectronVeto,TString CategoryID) {

    bool result = false;

    if(passesHadTowerOverEmCut(photon,CategoryID) && 
       passesChargedHadronCut(photon,rhocorPFChargedHadronIso,CategoryID) && 
       passesPhotonCut(photon,rhocorPFPhotonIso,CategoryID) && 
       passesNeutralHadronCut(photon,rhocorPFNeutralHadronIso,CategoryID) && 
       passesPFSigmaIetaIetaCut(photon,CategoryID) && 
       passesElectronVeto)
      result = true;

    return result; 

  }

  /*   // for offline spike flagging in our analysers */

  /*   bool isASpike(const reco::Photon *photon) { */

  /*     bool thisPhotonIsASpike = false; */

  /*     // can do spike identification based on regular swiss cross */
  /*     // or on out-of-time */
  /*     // even add double-spike determination in here too */

  /*     // can be used on MC also */
    
  /*     return thisPhotonIsASpike; */
  /*   } */

  /*   bool isBarrelPhoton(const reco::Photon *photon) { */

  /*     if ( fabs(photon->caloPosition().eta())<=1.4442) */
  /*       return true; */
  /*     else */
  /*       return false; */

  /*   } */


  /*   bool isGapPhoton(const reco::Photon *photon) { */

  /*     if ( fabs(photon->caloPosition().eta())>1.4442 && fabs(photon->caloPosition().eta())<1.566 ) return true;  */
  /*     else return false; */
  /*   } */


  bool isPFFakeableObject(const reco::Photon *photon,double thisCHIso,double thisNHIso,double thisPHIso,bool passesElectronVeto,TString CategoryID) {

    bool result = false;

    //Remember, these are already rho corrected isolation values that are provided
    
    // we are defining them to be loose, but also exclude real photons (to large degree) 
    // these 'swing values' determine both loose limit and exclusion veto

    double hadTowerOverEmCut = 0.;
    double sigmaIetaIetaCut = 0.;
    std::pair<double,double> CHIsoCut (0.,0.);
    std::pair<double,double> PHIsoCut (0.,0.);
    std::pair<double,double> NHIsoCut (0.,0.);

    //Set cut values for Single Tower H/E
    if(photon->isEB()){
      if(CategoryID.Contains("Tight")) hadTowerOverEmCut = 0.05;
      if(CategoryID.Contains("Medium")) hadTowerOverEmCut = 0.05;
      if(CategoryID.Contains("Loose")) hadTowerOverEmCut = 0.05;
    }

    if(photon->isEE()){
      if(CategoryID.Contains("Tight")) hadTowerOverEmCut = 0.05;
      if(CategoryID.Contains("Medium")) hadTowerOverEmCut = 0.05;
      if(CategoryID.Contains("Loose")) hadTowerOverEmCut = 0.05;
    }

    //Set cut values for Charged Hadron Isolation
    if(photon->isEB()){
      if(CategoryID.Contains("Tight")) {CHIsoCut.first = 0.7;CHIsoCut.second=0.;}
      if(CategoryID.Contains("Medium")) {CHIsoCut.first = 1.5;CHIsoCut.second=0.;}
      if(CategoryID.Contains("Loose")) {CHIsoCut.first = 2.6;CHIsoCut.second=0.;}
    }

    if(photon->isEE()){
      if(CategoryID.Contains("Tight")) {CHIsoCut.first = 0.5;CHIsoCut.second=0.;}
      if(CategoryID.Contains("Medium")) {CHIsoCut.first = 1.2;CHIsoCut.second=0.;}
      if(CategoryID.Contains("Loose")) {CHIsoCut.first = 2.3;CHIsoCut.second=0.;}
    }

    //Set cut values for Photon isolation
    if(photon->isEB()){
      if(CategoryID.Contains("Tight")) {PHIsoCut.first = 0.5;PHIsoCut.second=0.005;}
      if(CategoryID.Contains("Medium")) {PHIsoCut.first = 0.7;PHIsoCut.second=0.005;}
      if(CategoryID.Contains("Loose")) {PHIsoCut.first = 1.3;PHIsoCut.second=0.005;}
    }

    if(photon->isEE()){
      if(CategoryID.Contains("Tight")) {PHIsoCut.first = 1.0;PHIsoCut.second=0.005;}
      if(CategoryID.Contains("Medium")) {PHIsoCut.first = 1.0;PHIsoCut.second=0.005;}
      if(CategoryID.Contains("Loose")) {PHIsoCut.first = 99999.;PHIsoCut.second=99999.;}
    }

    //Set cut values for Neutral Hadron Isolation
    if(photon->isEB()){
      if(CategoryID.Contains("Tight")) {NHIsoCut.first = 0.4;NHIsoCut.second=0.04;}
      if(CategoryID.Contains("Medium")) {NHIsoCut.first = 1.0;NHIsoCut.second=0.04;}
      if(CategoryID.Contains("Loose")) {NHIsoCut.first = 3.5;NHIsoCut.second=0.04;}
    }

    if(photon->isEE()){
      if(CategoryID.Contains("Tight")) {NHIsoCut.first = 1.5;NHIsoCut.second=0.04;}
      if(CategoryID.Contains("Medium")) {NHIsoCut.first = 1.5;NHIsoCut.second=0.04;}
      if(CategoryID.Contains("Loose")) {NHIsoCut.first = 2.9;NHIsoCut.second=0.04;}
    }

    //Set cut values for sigmaIetaIeta 
    if(photon->isEB()){
      if(CategoryID.Contains("Tight")) sigmaIetaIetaCut = 0.011;
      if(CategoryID.Contains("Medium")) sigmaIetaIetaCut = 0.011;
      if(CategoryID.Contains("Loose")) sigmaIetaIetaCut = 0.012;
    }

    if(photon->isEE()){
      if(CategoryID.Contains("Tight")) sigmaIetaIetaCut = 0.031;
      if(CategoryID.Contains("Medium")) sigmaIetaIetaCut = 0.033;
      if(CategoryID.Contains("Loose")) sigmaIetaIetaCut = 0.034;
    }


    //For the moment we overwrite the parameters that define the Fakeable status
    //to be the always the ones that correspond to Loose 
    //Difference should not be so big, even in worst case
    //as only constant terms change


    if(photon->isEB()){
      hadTowerOverEmCut = 0.05;
      CHIsoCut.first = 2.6;CHIsoCut.second=0.;
      PHIsoCut.first = 1.3;PHIsoCut.second=0.005;
      NHIsoCut.first = 3.5;NHIsoCut.second=0.04;
      sigmaIetaIetaCut = 0.012;
    }



    if(photon->isEE()){
      hadTowerOverEmCut = 0.05;
      CHIsoCut.first = 2.3;CHIsoCut.second=0.;
      PHIsoCut.first = 99999.;PHIsoCut.second=99999.;
      NHIsoCut.first = 2.9;NHIsoCut.second=0.04;
      sigmaIetaIetaCut = 0.034;
    }


    double CHIsoSwingValue = CHIsoCut.first + CHIsoCut.second * photon->et();
    double CHIsoLooseLimit = TMath::Min( 5.0*(CHIsoSwingValue), 0.2*photon->et() );
    double CHIsoExclusion = CHIsoSwingValue;

    double PHIsoSwingValue = PHIsoCut.first + PHIsoCut.second * photon->et();
    double PHIsoLooseLimit = TMath::Min( 5.0*(PHIsoSwingValue), 0.2*photon->et() );
    double PHIsoExclusion = PHIsoSwingValue;

    double NHIsoSwingValue = NHIsoCut.first + NHIsoCut.second * photon->et();
    double NHIsoLooseLimit = TMath::Min( 5.0*(NHIsoSwingValue), 0.2*photon->et() );
    double NHIsoExclusion = NHIsoSwingValue;
    
    if( photon->hadTowOverEm() < hadTowerOverEmCut && 
	thisCHIso < CHIsoLooseLimit &&
	thisPHIso < PHIsoLooseLimit &&
	thisNHIso < NHIsoLooseLimit  && 
	//passesElectronVeto && 
	( thisCHIso > CHIsoExclusion || thisPHIso > PHIsoExclusion || thisNHIso > NHIsoExclusion || !passesPFSigmaIetaIetaCut(photon,CategoryID) )
	) {
      // then this passes our fakeable definition
      result = true;
    }
   
    return result;
  }


  std::vector<double> EffectiveAreas(const reco::Photon *photon) {

    std::vector<double> effarea;
    effarea.reserve(3);
    
    double effareaCH = 0.;
    double effareaNH = 0.;
    double effareaPH = 0.;


    //std::cout<<"photon abs eta in eff areas method in PFPhotonId.h "<<fabs(photon->eta())<<std::endl;

    //Setting effective areas for charged hadrons (CH), neutral hadrons (NH) and photons (PH)
    if(fabs(photon->eta()) < 1.){effareaCH = 0.012;effareaNH = 0.030;effareaPH = 0.148;}
    if( (fabs(photon->eta()) > 1.) && (fabs(photon->eta()) < 1.479) ){effareaCH = 0.010;effareaNH = 0.057;effareaPH = 0.130;}
    if( (fabs(photon->eta()) > 1.479) && (fabs(photon->eta()) < 2.) ){effareaCH = 0.014;effareaNH = 0.039;effareaPH = 0.112;}
    if( (fabs(photon->eta()) > 2.) && (fabs(photon->eta()) < 2.2) ){effareaCH = 0.012;effareaNH = 0.015;effareaPH = 0.216;}
    if( (fabs(photon->eta()) > 2.2) && (fabs(photon->eta()) < 2.3) ){effareaCH = 0.016;effareaNH = 0.024;effareaPH = 0.262;}
    if( (fabs(photon->eta()) > 2.3) && (fabs(photon->eta()) < 2.4) ){effareaCH = 0.020;effareaNH = 0.039;effareaPH = 0.260;}
    if(fabs(photon->eta()) > 2.4){effareaCH = 0.012;effareaNH = 0.072;effareaPH = 0.266;}
  
    effarea.push_back(effareaCH);
    effarea.push_back(effareaNH);
    effarea.push_back(effareaPH);

    return effarea;
  }

}

#endif

