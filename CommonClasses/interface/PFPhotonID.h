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

  bool passesHadTowerOverEmCut(const reco::Photon *photon,TString MethodID,TString CategoryID) {
    
    bool result = false;

    double hadTowerOverEmCut = 0.;

    if(MethodID.EqualTo("highpt")){
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
    }

    if(photon->hadTowOverEm()<hadTowerOverEmCut)
      result = true;

    return result;
  }


  bool passesChargedHadronCut(const reco::Photon *photon,double rhocorPFChargedHadronIso,TString MethodID,TString CategoryID) {
    
    bool result = false;

    std::pair<double,double> CHIsoCut (0.,0.);

    if(MethodID.EqualTo("highpt")){
      if(photon->isEB()){
	if(CategoryID.Contains("Tight")) {CHIsoCut.first = 5.;CHIsoCut.second=0.;}
	if(CategoryID.Contains("Medium")) {CHIsoCut.first = 5.;CHIsoCut.second=0.;}
	if(CategoryID.Contains("Loose")) {CHIsoCut.first = 5.;CHIsoCut.second=0.;}
      }

      if(photon->isEE()){
	if(CategoryID.Contains("Tight")) {CHIsoCut.first = 5.;CHIsoCut.second=0.;}
	if(CategoryID.Contains("Medium")) {CHIsoCut.first = 5.;CHIsoCut.second=0.;}
	if(CategoryID.Contains("Loose")) {CHIsoCut.first = 5.;CHIsoCut.second=0.;}
      }
    }

    if( rhocorPFChargedHadronIso < (CHIsoCut.first + CHIsoCut.second * photon->et()))
      result = true;

    return result;
  }


  bool passesPhotonCut(const reco::Photon *photon,double rhocorPFPhotonIso,TString MethodID,TString CategoryID) {
    
    bool result = false;

    std::pair<double,double> PHIsoCut (0.0, 0.0);

    if (MethodID.EqualTo("highpt")) {
      
      if (photon->isEB()) {
	if (CategoryID.Contains("Tight"))  {PHIsoCut.first = 0.25; PHIsoCut.second = 0.0045;}
	if (CategoryID.Contains("Medium")) {PHIsoCut.first = 0.25; PHIsoCut.second = 0.0045;}
	if (CategoryID.Contains("Loose"))  {PHIsoCut.first = 0.25; PHIsoCut.second = 0.0045;}
      }

      // 1.560 or 1.566?
      if (1.560 < fabs(photon->superCluster()->eta()) && fabs(photon->superCluster()->eta()) < 2.0) {
	if (CategoryID.Contains("Tight"))  {PHIsoCut.first = -0.5; PHIsoCut.second = 0.0045;}
	if (CategoryID.Contains("Medium")) {PHIsoCut.first = -0.5; PHIsoCut.second = 0.0045;}
	if (CategoryID.Contains("Loose"))  {PHIsoCut.first = -0.5; PHIsoCut.second = 0.0045;}
      }

      if (2.0 < fabs(photon->superCluster()->eta())) {
	if (CategoryID.Contains("Tight"))  {PHIsoCut.first = -0.5; PHIsoCut.second = 0.0030;}
	if (CategoryID.Contains("Medium")) {PHIsoCut.first = -0.5; PHIsoCut.second = 0.0030;}
	if (CategoryID.Contains("Loose"))  {PHIsoCut.first = -0.5; PHIsoCut.second = 0.0030;}
      }
      
    } // end "highpt"

    if (rhocorPFPhotonIso < (PHIsoCut.first + PHIsoCut.second * photon->pt()))
      result = true;

    return result;
  }


  bool passesNeutralHadronCut(const reco::Photon *photon,double rhocorPFNeutralHadronIso,TString MethodID,TString CategoryID) {
    
    bool result = false;

    std::pair<double,double> NHIsoCut (0.,0.);

    if(MethodID.EqualTo("highpt")){
      if(photon->isEB()){
	if(CategoryID.Contains("Tight")) {NHIsoCut.first = 99999.;NHIsoCut.second=0.;}
	if(CategoryID.Contains("Medium")) {NHIsoCut.first = 99999.;NHIsoCut.second=0.;}
	if(CategoryID.Contains("Loose")) {NHIsoCut.first = 99999.;NHIsoCut.second=0.;}
      }

      if(photon->isEE()){
	if(CategoryID.Contains("Tight")) {NHIsoCut.first = 99999.;NHIsoCut.second=0.;}
	if(CategoryID.Contains("Medium")) {NHIsoCut.first = 99999.;NHIsoCut.second=0.;}
	if(CategoryID.Contains("Loose")) {NHIsoCut.first = 99999.;NHIsoCut.second=0.;}
      }
    }

    if( rhocorPFNeutralHadronIso < (NHIsoCut.first + NHIsoCut.second * photon->et()))
      result = true;

    return result;
  }


  bool passesPFSigmaIetaIetaCut(const reco::Photon *photon, double full5x5SigmaIetaIeta, TString MethodID, TString CategoryID, bool isSaturated)
  {  
    bool result = false;

    double sigmaIetaIetaCut = 0.;

    if (MethodID.EqualTo("highpt")) {
      
      if (photon->isEB()) {
	if (CategoryID.Contains("Tight") || CategoryID.Contains("Medium") || CategoryID.Contains("Loose")) {
	  if (isSaturated) sigmaIetaIetaCut = 0.0112;
	  else             sigmaIetaIetaCut = 0.0105;
	}
      }
      
      if (photon->isEE()) {
	if (CategoryID.Contains("Tight") || CategoryID.Contains("Medium") || CategoryID.Contains("Loose")) {
	  if (isSaturated) sigmaIetaIetaCut = 0.0300;
	  else             sigmaIetaIetaCut = 0.0280;
	}
      }
      
    } // end "highpt"

    if (full5x5SigmaIetaIeta < sigmaIetaIetaCut)
      result = true;

    return result;
  }


  bool isPFTightPhoton(const reco::Photon *photon,double rhocorPFChargedHadronIso,double rhocorPFNeutralHadronIso,double rhocorPFPhotonIso,double full5x5SigmaIetaIeta,bool passesElectronVeto,TString MethodID,TString CategoryID,bool isSaturated) {

    bool result = false;

    if(passesHadTowerOverEmCut(photon,MethodID,CategoryID) &&
       passesChargedHadronCut(photon,rhocorPFChargedHadronIso,MethodID,CategoryID) &&
       passesPhotonCut(photon,rhocorPFPhotonIso,MethodID,CategoryID) &&
       passesNeutralHadronCut(photon,rhocorPFNeutralHadronIso,MethodID,CategoryID) &&
       passesPFSigmaIetaIetaCut(photon,full5x5SigmaIetaIeta,MethodID,CategoryID,isSaturated) &&
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


  bool isPFFakeableObject(const reco::Photon *photon, double thisCHIso, double thisNHIso, double thisPHIso, double thisSigmaIetaIeta,
			  bool passesElectronVeto, TString MethodID, TString CategoryID, bool isSaturated) {

    bool result = false;

    //Remember, these are already rho corrected isolation values that are provided
    
    // we are defining them to be loose, but also exclude real photons (to large degree) 
    // these 'swing values' determine both loose limit and exclusion veto

    double hadTowerOverEmCut = 0.;
    double sigmaIetaIetaCut = 0.;

    float arrCHIso[] = {0.,0.,0.};
    float arrNHIso[] = {0.,0.,0.};
    float arrPHIso[] = {0.,0.,0.};

    std::vector<float> CHIsoCut(arrCHIso, arrCHIso + sizeof(arrCHIso) / sizeof(arrCHIso[0]) );
    std::vector<float> NHIsoCut(arrNHIso, arrNHIso + sizeof(arrNHIso) / sizeof(arrNHIso[0]) );
    std::vector<float> PHIsoCut(arrPHIso, arrPHIso + sizeof(arrPHIso) / sizeof(arrPHIso[0]) );

    //Set cut values for Single Tower H/E
    if(MethodID.EqualTo("egamma")){
      if(photon->isEB()){
	if(CategoryID.Contains("Tight"))  hadTowerOverEmCut = 0.05;
	if(CategoryID.Contains("Medium")) hadTowerOverEmCut = 0.05;
	if(CategoryID.Contains("Loose"))  hadTowerOverEmCut = 0.05;
      }

      if(photon->isEE()){
	if(CategoryID.Contains("Tight"))  hadTowerOverEmCut = 0.05;
	if(CategoryID.Contains("Medium")) hadTowerOverEmCut = 0.05;
	if(CategoryID.Contains("Loose"))  hadTowerOverEmCut = 0.05;
      }
    }

    if(MethodID.EqualTo("highpt")){
      if(photon->isEB()){
	if(CategoryID.Contains("Tight"))  hadTowerOverEmCut = 0.05;
	if(CategoryID.Contains("Medium")) hadTowerOverEmCut = 0.05;
	if(CategoryID.Contains("Loose"))  hadTowerOverEmCut = 0.05;
      }

      if(photon->isEE()){
	if(CategoryID.Contains("Tight"))  hadTowerOverEmCut = 0.05;
	if(CategoryID.Contains("Medium")) hadTowerOverEmCut = 0.05;
	if(CategoryID.Contains("Loose"))  hadTowerOverEmCut = 0.05;
      }
    }


    //Set cut values for Charged Hadron Isolation
    if(MethodID.EqualTo("egamma")){
      if(photon->isEB()){
	if(CategoryID.Contains("Tight"))  {CHIsoCut[0] = 0.; CHIsoCut[1] = 0.91; CHIsoCut[2] = 0.;}
	if(CategoryID.Contains("Medium")) {CHIsoCut[0] = 0.; CHIsoCut[1] = 1.31; CHIsoCut[2] = 0.;}
	if(CategoryID.Contains("Loose"))  {CHIsoCut[0] = 0.; CHIsoCut[1] = 2.44; CHIsoCut[2] = 0.;}
      }

      if(photon->isEE()){
	if(CategoryID.Contains("Tight"))  {CHIsoCut[0] = 0.; CHIsoCut[1] = 0.65; CHIsoCut[2] = 0.;}
	if(CategoryID.Contains("Medium")) {CHIsoCut[0] = 0.; CHIsoCut[1] = 1.25; CHIsoCut[2] = 0.;}
	if(CategoryID.Contains("Loose"))  {CHIsoCut[0] = 0.; CHIsoCut[1] = 1.84; CHIsoCut[2] = 0.;}
      }
    }

    if(MethodID.EqualTo("highpt")){
      if(photon->isEB()){
	if(CategoryID.Contains("Tight"))  {CHIsoCut[0] = 0.; CHIsoCut[1] = 5.; CHIsoCut[2] = 0.;}
	if(CategoryID.Contains("Medium")) {CHIsoCut[0] = 0.; CHIsoCut[1] = 5.; CHIsoCut[2] = 0.;}
	if(CategoryID.Contains("Loose"))  {CHIsoCut[0] = 0.; CHIsoCut[1] = 5.; CHIsoCut[2] = 0.;}
      }

      if(photon->isEE()){
	if(CategoryID.Contains("Tight"))  {CHIsoCut[0] = 0.; CHIsoCut[1] = 5.; CHIsoCut[2] = 0.;}
	if(CategoryID.Contains("Medium")) {CHIsoCut[0] = 0.; CHIsoCut[1] = 5.; CHIsoCut[2] = 0.;}
	if(CategoryID.Contains("Loose"))  {CHIsoCut[0] = 0.; CHIsoCut[1] = 5.; CHIsoCut[2] = 0.;}
      }
    }

    //Set cut values for Photon isolation
    if(MethodID.EqualTo("egamma")){
      if(photon->isEB()){
	if(CategoryID.Contains("Tight"))  {PHIsoCut[0] = 0.; PHIsoCut[1] = 0.61; PHIsoCut[2] = 0.0043;}
	if(CategoryID.Contains("Medium")) {PHIsoCut[0] = 0.; PHIsoCut[1] = 1.33; PHIsoCut[2] = 0.0043;}
	if(CategoryID.Contains("Loose"))  {PHIsoCut[0] = 0.; PHIsoCut[1] = 1.92; PHIsoCut[2] = 0.0043;}
      }

      if(photon->isEE()){
	if(CategoryID.Contains("Tight"))  {PHIsoCut[0] = 0.; PHIsoCut[1] = 0.54; PHIsoCut[2] = 0.0041;}
	if(CategoryID.Contains("Medium")) {PHIsoCut[0] = 0.; PHIsoCut[1] = 1.02; PHIsoCut[2] = 0.0041;}
	if(CategoryID.Contains("Loose"))  {PHIsoCut[0] = 0.; PHIsoCut[1] = 2.15; PHIsoCut[2] = 0.0041;}
      }
    }

    if (MethodID.EqualTo("highpt")) {
      if (photon->isEB()){
	if(CategoryID.Contains("Tight"))  {PHIsoCut[0] = 0.; PHIsoCut[1] = 0.25; PHIsoCut[2] = 0.0045;}
	if(CategoryID.Contains("Medium")) {PHIsoCut[0] = 0.; PHIsoCut[1] = 0.25; PHIsoCut[2] = 0.0045;}
	if(CategoryID.Contains("Loose"))  {PHIsoCut[0] = 0.; PHIsoCut[1] = 0.25; PHIsoCut[2] = 0.0045;}
      }
      
      if (1.560 < fabs(photon->superCluster()->eta()) && fabs(photon->superCluster()->eta()) < 2.0) {
	if(CategoryID.Contains("Tight"))  {PHIsoCut[0] = 0.; PHIsoCut[1] = -0.5; PHIsoCut[2] = 0.0045;}
	if(CategoryID.Contains("Medium")) {PHIsoCut[0] = 0.; PHIsoCut[1] = -0.5; PHIsoCut[2] = 0.0045;}
	if(CategoryID.Contains("Loose"))  {PHIsoCut[0] = 0.; PHIsoCut[1] = -0.5; PHIsoCut[2] = 0.0045;}
      }
      
      if (2.0 < fabs(photon->superCluster()->eta())) {
	if(CategoryID.Contains("Tight"))  {PHIsoCut[0] = 0.; PHIsoCut[1] = -0.5; PHIsoCut[2] = 0.0030;}
	if(CategoryID.Contains("Medium")) {PHIsoCut[0] = 0.; PHIsoCut[1] = -0.5; PHIsoCut[2] = 0.0030;}
	if(CategoryID.Contains("Loose"))  {PHIsoCut[0] = 0.; PHIsoCut[1] = -0.5; PHIsoCut[2] = 0.0030;}
      }
    }

    //Set cut values for Neutral Hadron Isolation
    if(MethodID.EqualTo("egamma")){
      if(photon->isEB()){
	if(CategoryID.Contains("Tight"))  {NHIsoCut[0] = 0.33; NHIsoCut[1] = 0.5809; NHIsoCut[2] = 0.0044;}
	if(CategoryID.Contains("Medium")) {NHIsoCut[0] = 0.60; NHIsoCut[1] = 0.5809; NHIsoCut[2] = 0.0044;}
	if(CategoryID.Contains("Loose"))  {NHIsoCut[0] = 2.57; NHIsoCut[1] = 0.5809; NHIsoCut[2] = 0.0044;}
      }

      if(photon->isEE()){
	if(CategoryID.Contains("Tight"))  {NHIsoCut[0] = 0.93; NHIsoCut[1] = 0.9402; NHIsoCut[2] = 0.0040;}
	if(CategoryID.Contains("Medium")) {NHIsoCut[0] = 1.65; NHIsoCut[1] = 0.9402; NHIsoCut[2] = 0.0040;}
	if(CategoryID.Contains("Loose"))  {NHIsoCut[0] = 4.00; NHIsoCut[1] = 0.9402; NHIsoCut[2] = 0.0040;}
      }
    }

    if(MethodID.EqualTo("highpt")){
      if(photon->isEB()){
	if(CategoryID.Contains("Tight"))  {NHIsoCut[0] = 99999.; NHIsoCut[1] = 99999.; NHIsoCut[2] = 99999.;}
	if(CategoryID.Contains("Medium")) {NHIsoCut[0] = 99999.; NHIsoCut[1] = 99999.; NHIsoCut[2] = 99999.;}
	if(CategoryID.Contains("Loose"))  {NHIsoCut[0] = 99999.; NHIsoCut[1] = 99999.; NHIsoCut[2] = 99999.;}
      }

      if(photon->isEE()){
	if(CategoryID.Contains("Tight"))  {NHIsoCut[0] = 99999.; NHIsoCut[1] = 99999.; NHIsoCut[2] = 99999.;}
	if(CategoryID.Contains("Medium")) {NHIsoCut[0] = 99999.; NHIsoCut[1] = 99999.; NHIsoCut[2] = 99999.;}
	if(CategoryID.Contains("Loose"))  {NHIsoCut[0] = 99999.; NHIsoCut[1] = 99999.; NHIsoCut[2] = 99999.;}
      }
    }

    //Set cut values for sigmaIetaIeta
    if(MethodID.EqualTo("egamma")){
      if(photon->isEB()){
	if(CategoryID.Contains("Tight"))  sigmaIetaIetaCut = 0.0100;
	if(CategoryID.Contains("Medium")) sigmaIetaIetaCut = 0.0100;
	if(CategoryID.Contains("Loose"))  sigmaIetaIetaCut = 0.0103;
      }

      if(photon->isEE()){
	if(CategoryID.Contains("Tight"))  sigmaIetaIetaCut = 0.0267;
	if(CategoryID.Contains("Medium")) sigmaIetaIetaCut = 0.0267;
	if(CategoryID.Contains("Loose"))  sigmaIetaIetaCut = 0.0277;
      }
    }

    if(MethodID.EqualTo("highpt")){
      if (photon->isEB()) {
	if (CategoryID.Contains("Tight") || CategoryID.Contains("Medium") || CategoryID.Contains("Loose")) {
	  if (isSaturated) sigmaIetaIetaCut = 0.0112;
	  else             sigmaIetaIetaCut = 0.0105;
	}
      }

      if (photon->isEE()) {
	if (CategoryID.Contains("Tight") || CategoryID.Contains("Medium") || CategoryID.Contains("Loose")) {
	  if (isSaturated) sigmaIetaIetaCut = 0.0300;
	  else             sigmaIetaIetaCut = 0.0280;
	}
      }
    }


    double CHIsoSwingValue = CHIsoCut[1] + CHIsoCut[2] * photon->et();
    double CHIsoLooseLimit = TMath::Min( 5.0*(CHIsoSwingValue), 0.2*photon->et() );
    double CHIsoExclusion = CHIsoSwingValue;

    double PHIsoSwingValue = PHIsoCut[1] + PHIsoCut[2] * photon->et();
    double PHIsoLooseLimit = TMath::Min( 5.0*(PHIsoSwingValue), 0.2*photon->et() );
    double PHIsoExclusion = PHIsoSwingValue;

    //Special case for neutral hadron isolation
    //no cut in high pt photon ID method
    //only in egamma loose method
    double NHIsoSwingValue = 99999.;
    double NHIsoLooseLimit = 99999.;
    double NHIsoExclusion = 99999.;

    if(MethodID.EqualTo("egamma")){
      //if(photon->isEB())
      NHIsoSwingValue = NHIsoCut[0] + exp(NHIsoCut[1] + NHIsoCut[2] * photon->pt());
      //if(photon->isEE()) NHIsoSwingValue = NHIsoCut[1] + NHIsoCut[2] * photon->et();
      NHIsoLooseLimit = TMath::Min( 5.0*(NHIsoSwingValue), 0.2*photon->et() );
      NHIsoExclusion = NHIsoSwingValue;
    }

    std::cout << sigmaIetaIetaCut << std::endl;

    if( photon->hadTowOverEm() < hadTowerOverEmCut &&
    	thisCHIso < CHIsoLooseLimit &&
    	thisPHIso < PHIsoLooseLimit &&
    	thisNHIso < NHIsoLooseLimit  &&
    	passesElectronVeto &&
    	( thisCHIso > CHIsoExclusion || thisPHIso > PHIsoExclusion || thisNHIso > NHIsoExclusion || thisSigmaIetaIeta > sigmaIetaIetaCut)
    	) {
      // then this passes our fakeable definition
      result = true;
    }
   
    return result;
  }


  std::vector<double> EffectiveAreas(const reco::Photon *photon,TString MethodID,TString CategoryID){
    
    std::vector<double> effarea;
    effarea.reserve(3);
    
    double effareaCH = 0.;
    double effareaNH = 0.;
    double effareaPH = 0.;
    
    
    //std::cout<<"photon abs eta in eff areas method in PFPhotonId.h "<<fabs(photon->superCluster()->eta())<<std::endl;

    //Setting effective areas for charged hadrons (CH), neutral hadrons (NH) and photons (PH)

    // Spring15 50 ns (from RecoEgamma/PhotonIdentification/data/Spring15)
    if (MethodID.EqualTo("egamma")){
      if (fabs(photon->superCluster()->eta()) < 1.0)                                                  {effareaCH = 0.0157; effareaNH = 0.0143; effareaPH = 0.0725;}
      if (1.0   < fabs(photon->superCluster()->eta()) && fabs(photon->superCluster()->eta()) < 1.479) {effareaCH = 0.0143; effareaNH = 0.0210; effareaPH = 0.0604;}
      if (1.479 < fabs(photon->superCluster()->eta()) && fabs(photon->superCluster()->eta()) < 2.0)   {effareaCH = 0.0115; effareaNH = 0.0148; effareaPH = 0.0320;}
      if (2.0   < fabs(photon->superCluster()->eta()) && fabs(photon->superCluster()->eta()) < 2.2)   {effareaCH = 0.0094; effareaNH = 0.0082; effareaPH = 0.0512;}
      if (2.2   < fabs(photon->superCluster()->eta()) && fabs(photon->superCluster()->eta()) < 2.3)   {effareaCH = 0.0095; effareaNH = 0.0124; effareaPH = 0.0766;}
      if (2.3   < fabs(photon->superCluster()->eta()) && fabs(photon->superCluster()->eta()) < 2.4)   {effareaCH = 0.0068; effareaNH = 0.0186; effareaPH = 0.0949;}
      if (fabs(photon->superCluster()->eta()) > 2.4)                                                  {effareaCH = 0.0053; effareaNH = 0.0320; effareaPH = 0.1160;}
    }

    // 1.560 or 1.566?
    if (MethodID.EqualTo("highpt")){
      if (fabs(photon->superCluster()->eta()) < 0.9)                                                   {effareaPH = 0.17;}
      if (0.9   < fabs(photon->superCluster()->eta()) && fabs(photon->superCluster()->eta()) < 1.4442) {effareaPH = 0.14;}
      if (1.560 < fabs(photon->superCluster()->eta()) && fabs(photon->superCluster()->eta()) < 2.0)    {effareaPH = 0.11;}
      if (2.0   < fabs(photon->superCluster()->eta()) && fabs(photon->superCluster()->eta()) < 2.2)    {effareaPH = 0.14;}
      if (2.2   < fabs(photon->superCluster()->eta()) && fabs(photon->superCluster()->eta()) < 2.5)    {effareaPH = 0.22;}
    }

    effarea.push_back(effareaCH);
    effarea.push_back(effareaNH);
    effarea.push_back(effareaPH);

    return effarea;
  }

  double alphaPhotonHighPtID(const reco::Photon *photon)
  {

    double alpha = 0;

    // 1.560 or 1.566?
    if (fabs(photon->superCluster()->eta()) < 0.9)                                                    {alpha = 2.5;}
    if (0.9 < fabs(photon->superCluster()->eta())   && fabs(photon->superCluster()->eta()) < 1.4442 ) {alpha = 2.5;}
    if (1.560 < fabs(photon->superCluster()->eta()) && fabs(photon->superCluster()->eta()) < 2.0 )    {alpha = 2.5;}
    if (2.0 < fabs(photon->superCluster()->eta())   && fabs(photon->superCluster()->eta()) < 2.2 )    {alpha = 2.5;}
    if (2.2 < fabs(photon->superCluster()->eta())   && fabs(photon->superCluster()->eta()) < 2.5 )    {alpha = 2.5;}

    return alpha;
  }

  double kappaPhotonHighPtID(const reco::Photon *photon)
  {

    double kappa = 0;

    // 1.560 or 1.566?
    if (fabs(photon->superCluster()->eta()) < 0.9)                                                    {kappa = 0.0045;}
    if (0.9 < fabs(photon->superCluster()->eta())   && fabs(photon->superCluster()->eta()) < 1.4442 ) {kappa = 0.0045;}
    if (1.560 < fabs(photon->superCluster()->eta()) && fabs(photon->superCluster()->eta()) < 2.0 )    {kappa = 0.0045;}
    if (2.0 < fabs(photon->superCluster()->eta())   && fabs(photon->superCluster()->eta()) < 2.2 )    {kappa = 0.0030;}
    if (2.2 < fabs(photon->superCluster()->eta())   && fabs(photon->superCluster()->eta()) < 2.5 )    {kappa = 0.0030;}

    return kappa;
  }

  double corPhoIsoHighPtID(const reco::Photon *photon, TString MethodID, TString CategoryID, double phoIso, double rho)
  {
    std::vector<double> EA = EffectiveAreas(photon,MethodID,CategoryID);
    double corPhoIso =
      alphaPhotonHighPtID(photon) + phoIso - rho*EA[2] - kappaPhotonHighPtID(photon)*photon->pt();
    
    return corPhoIso;
  }

  bool passCorPhoIsoHighPtID(const reco::Photon *photon, TString MethodID, TString CategoryID, double phoIso, double rho)
  {
    bool passCorPhoIso = false;
    double corPhoIsoCut = 0;
    double corPhoIso = corPhoIsoHighPtID(photon, MethodID, CategoryID, phoIso, rho);

    if (photon->isEB()) corPhoIsoCut = 2.75;
    if (photon->isEE()) corPhoIsoCut = 2.00;

    if (corPhoIso < corPhoIsoCut) passCorPhoIso = true;
    
    return passCorPhoIso;
  }
  
  bool passHighPtID(const reco::Photon *photon, TString MethodID, TString CategoryID, double chIso, double phoIso, double sigIeIe, double rho, bool passCSEV, bool isSat)
  {
    return (passesHadTowerOverEmCut(photon, MethodID, CategoryID) &&
	    passesChargedHadronCut(photon, chIso, MethodID, CategoryID) &&
	    passCorPhoIsoHighPtID(photon, MethodID, CategoryID, phoIso, rho) &&
	    passesPFSigmaIetaIetaCut(photon, sigIeIe, MethodID, CategoryID, isSat) &&
	    passCSEV);
  }
  
}

#endif

