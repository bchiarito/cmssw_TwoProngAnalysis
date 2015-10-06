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

    std::pair<double,double> PHIsoCut (0.,0.);

    if(MethodID.EqualTo("highpt")){
      if(photon->isEB()){
	if(CategoryID.Contains("Tight")) {PHIsoCut.first = 1.;PHIsoCut.second=0.002;}
	if(CategoryID.Contains("Medium")) {PHIsoCut.first = 1.;PHIsoCut.second=0.002;}
	if(CategoryID.Contains("Loose")) {PHIsoCut.first = 1.;PHIsoCut.second=0.002;}
      }

      if(photon->isEE()){
	if(CategoryID.Contains("Tight")) {PHIsoCut.first = 0.;PHIsoCut.second=0.002;}
	if(CategoryID.Contains("Medium")) {PHIsoCut.first = 0.;PHIsoCut.second=0.002;}
	if(CategoryID.Contains("Loose")) {PHIsoCut.first = 0.;PHIsoCut.second=0.002;}
      }
    }

    if( rhocorPFPhotonIso < (PHIsoCut.first + PHIsoCut.second * photon->et()))
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


  bool passesPFSigmaIetaIetaCut(const reco::Photon *photon,double full5x5SigmaIetaIeta,TString MethodID,TString CategoryID) {
    
    bool result = false;

    double sigmaIetaIetaCut = 0.;

    if(MethodID.EqualTo("highpt")){
      if(photon->isEB()){
	if(CategoryID.Contains("Tight")) sigmaIetaIetaCut = 0.0105;
	if(CategoryID.Contains("Medium")) sigmaIetaIetaCut = 0.0105;
	if(CategoryID.Contains("Loose")) sigmaIetaIetaCut = 0.0105;
      }

      if(photon->isEE()){
	if(CategoryID.Contains("Tight")) sigmaIetaIetaCut = 0.028;
	if(CategoryID.Contains("Medium")) sigmaIetaIetaCut = 0.028;
	if(CategoryID.Contains("Loose")) sigmaIetaIetaCut = 0.028;
      }
    }

    if(full5x5SigmaIetaIeta<sigmaIetaIetaCut)
      result = true;

    return result;
  }



  bool isPFTightPhoton(const reco::Photon *photon,double rhocorPFChargedHadronIso,double rhocorPFNeutralHadronIso,double rhocorPFPhotonIso,double full5x5SigmaIetaIeta,bool passesElectronVeto,TString MethodID,TString CategoryID) {

    bool result = false;

    if(passesHadTowerOverEmCut(photon,MethodID,CategoryID) &&
       passesChargedHadronCut(photon,rhocorPFChargedHadronIso,MethodID,CategoryID) &&
       passesPhotonCut(photon,rhocorPFPhotonIso,MethodID,CategoryID) &&
       passesNeutralHadronCut(photon,rhocorPFNeutralHadronIso,MethodID,CategoryID) &&
       passesPFSigmaIetaIetaCut(photon,full5x5SigmaIetaIeta,MethodID,CategoryID) &&
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


  bool isPFFakeableObject(const reco::Photon *photon,double thisCHIso,double thisNHIso,double thisPHIso,double thisSigmaIetaIeta, bool passesElectronVeto,TString MethodID,TString CategoryID) {

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
	if(CategoryID.Contains("Tight")) hadTowerOverEmCut = 0.010;
	if(CategoryID.Contains("Medium")) hadTowerOverEmCut = 0.012;
	if(CategoryID.Contains("Loose")) hadTowerOverEmCut = 0.028;
      }

      if(photon->isEE()){
	if(CategoryID.Contains("Tight")) hadTowerOverEmCut = 0.015;
	if(CategoryID.Contains("Medium")) hadTowerOverEmCut = 0.023;
	if(CategoryID.Contains("Loose")) hadTowerOverEmCut = 0.093;
      }
    }

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


    //Set cut values for Charged Hadron Isolation
    if(MethodID.EqualTo("egamma")){
      if(photon->isEB()){
	if(CategoryID.Contains("Tight")) {CHIsoCut[0] = 0.;CHIsoCut[1]=1.66;CHIsoCut[2]=0.;}
	if(CategoryID.Contains("Medium")) {CHIsoCut[0] = 0.;CHIsoCut[1]=1.79;CHIsoCut[2]=0.;}
	if(CategoryID.Contains("Loose")) {CHIsoCut[0] = 0.;CHIsoCut[1]=2.67;CHIsoCut[2]=0.;}
      }

      if(photon->isEE()){
	if(CategoryID.Contains("Tight")) {CHIsoCut[0] = 0.;CHIsoCut[1]=1.04;CHIsoCut[2]=0.;}
	if(CategoryID.Contains("Medium")) {CHIsoCut[0] = 0.;CHIsoCut[1]=1.09;CHIsoCut[2]=0.;}
	if(CategoryID.Contains("Loose")) {CHIsoCut[0] = 0.;CHIsoCut[1]=1.79;CHIsoCut[2]=0.;}
      }
    }

    if(MethodID.EqualTo("highpt")){
      if(photon->isEB()){
	if(CategoryID.Contains("Tight")) {CHIsoCut[0] = 0.;CHIsoCut[1]=5.;CHIsoCut[2]=0.;}
	if(CategoryID.Contains("Medium")) {CHIsoCut[0] = 0.;CHIsoCut[1]=5.;CHIsoCut[2]=0.;}
	if(CategoryID.Contains("Loose")) {CHIsoCut[0] = 0.;CHIsoCut[1]=5.;CHIsoCut[2]=0.;}
      }

      if(photon->isEE()){
	if(CategoryID.Contains("Tight")) {CHIsoCut[0] = 0.;CHIsoCut[1]=5.;CHIsoCut[2]=0.;}
	if(CategoryID.Contains("Medium")) {CHIsoCut[0] = 0.;CHIsoCut[1]=5.;CHIsoCut[2]=0.;}
	if(CategoryID.Contains("Loose")) {CHIsoCut[0] = 0.;CHIsoCut[1]=5.;CHIsoCut[2]=0.;}
      }
    }

    //Set cut values for Photon isolation
    if(MethodID.EqualTo("egamma")){
      if(photon->isEB()){
	if(CategoryID.Contains("Tight")) {PHIsoCut[0] = 0.;PHIsoCut[1]=1.40;PHIsoCut[2]=0.0014;}
	if(CategoryID.Contains("Medium")) {PHIsoCut[0] = 0.;PHIsoCut[1]=1.90;PHIsoCut[2]=0.0014;}
	if(CategoryID.Contains("Loose")) {PHIsoCut[0] = 0.;PHIsoCut[1]=2.11;PHIsoCut[2]=0.0014;}
      }

      if(photon->isEE()){
	if(CategoryID.Contains("Tight")) {PHIsoCut[0] = 0.;PHIsoCut[1]=1.40;PHIsoCut[2]=0.0091;}
	if(CategoryID.Contains("Medium")) {PHIsoCut[0] = 0.;PHIsoCut[1]=1.90;PHIsoCut[2]=0.0091;}
	if(CategoryID.Contains("Loose")) {PHIsoCut[0] = 0.;PHIsoCut[1]=3.09;PHIsoCut[2]=0.0091;}
      }
    }

    if(MethodID.EqualTo("highpt")){
      if(photon->isEB()){
	if(CategoryID.Contains("Tight")) {PHIsoCut[0] = 0.;PHIsoCut[1]=1.;PHIsoCut[2]=0.002;}
	if(CategoryID.Contains("Medium")) {PHIsoCut[0] = 0.;PHIsoCut[1]=1.;PHIsoCut[2]=0.002;}
	if(CategoryID.Contains("Loose")) {PHIsoCut[0] = 0.;PHIsoCut[1]=1.;PHIsoCut[2]=0.002;}
      }

      if(photon->isEE()){
	if(CategoryID.Contains("Tight")) {PHIsoCut[0] = 0.;PHIsoCut[1]=0.;PHIsoCut[2]=0.002;}
	if(CategoryID.Contains("Medium")) {PHIsoCut[0] = 0.;PHIsoCut[1]=0.;PHIsoCut[2]=0.002;}
	if(CategoryID.Contains("Loose")) {PHIsoCut[0] = 0.;PHIsoCut[1]=0.;PHIsoCut[2]=0.002;}
      }
    }

    //Set cut values for Neutral Hadron Isolation
    if(MethodID.EqualTo("egamma")){
      if(photon->isEB()){
	if(CategoryID.Contains("Tight")) {NHIsoCut[0] = 0.14;NHIsoCut[1]=0.5408;NHIsoCut[2]=0.0028;}
	if(CategoryID.Contains("Medium")) {NHIsoCut[0] = 0.16;NHIsoCut[1]=0.5408;NHIsoCut[2]=0.0028;}
	if(CategoryID.Contains("Loose")) {NHIsoCut[0] = 7.23;NHIsoCut[1]=0.5408;NHIsoCut[2]=0.0028;}
      }

      if(photon->isEE()){
	if(CategoryID.Contains("Tight")) {NHIsoCut[0] = 0.;NHIsoCut[1]=3.89;NHIsoCut[2]=0.0172;}
	if(CategoryID.Contains("Medium")) {NHIsoCut[0] = 0.;NHIsoCut[1]=4.31;NHIsoCut[2]=0.0172;}
	if(CategoryID.Contains("Loose")) {NHIsoCut[0] = 0.;NHIsoCut[1]=8.89;NHIsoCut[2]=0.01725;}
      }
    }

    if(MethodID.EqualTo("highpt")){
      if(photon->isEB()){
	if(CategoryID.Contains("Tight")) {NHIsoCut[0] = 99999.;NHIsoCut[1]=99999.;NHIsoCut[2]=99999.;}
	if(CategoryID.Contains("Medium")) {NHIsoCut[0] = 99999.;NHIsoCut[1]=99999.;NHIsoCut[2]=99999.;}
	if(CategoryID.Contains("Loose")) {NHIsoCut[0] = 99999.;NHIsoCut[1]=99999.;NHIsoCut[2]=99999.;}
      }

      if(photon->isEE()){
	if(CategoryID.Contains("Tight")) {NHIsoCut[0] = 99999.;NHIsoCut[1]=99999.;NHIsoCut[2]=99999.;}
	if(CategoryID.Contains("Medium")) {NHIsoCut[0] = 99999.;NHIsoCut[1]=99999.;NHIsoCut[2]=99999.;}
	if(CategoryID.Contains("Loose")) {NHIsoCut[0] = 99999.;NHIsoCut[1]=99999.;NHIsoCut[2]=99999.;}
      }
    }

    //Set cut values for sigmaIetaIeta
    if(MethodID.EqualTo("egamma")){
      if(photon->isEB()){
	if(CategoryID.Contains("Tight")) sigmaIetaIetaCut = 0.0100;
	if(CategoryID.Contains("Medium")) sigmaIetaIetaCut = 0.0100;
	if(CategoryID.Contains("Loose")) sigmaIetaIetaCut = 0.0107;
      }

      if(photon->isEE()){
	if(CategoryID.Contains("Tight")) sigmaIetaIetaCut = 0.0265;
	if(CategoryID.Contains("Medium")) sigmaIetaIetaCut = 0.0267;
	if(CategoryID.Contains("Loose")) sigmaIetaIetaCut = 0.0272;
      }
    }

    if(MethodID.EqualTo("highpt")){
      if(photon->isEB()){
	if(CategoryID.Contains("Tight")) sigmaIetaIetaCut = 0.0105;
	if(CategoryID.Contains("Medium")) sigmaIetaIetaCut = 0.0105;
	if(CategoryID.Contains("Loose")) sigmaIetaIetaCut = 0.0105;
      }

      if(photon->isEE()){
	if(CategoryID.Contains("Tight")) sigmaIetaIetaCut = 0.028;
	if(CategoryID.Contains("Medium")) sigmaIetaIetaCut = 0.028;
	if(CategoryID.Contains("Loose")) sigmaIetaIetaCut = 0.028;
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
      if(photon->isEB()) NHIsoSwingValue = NHIsoCut[0] + exp(NHIsoCut[1] + NHIsoCut[2] * photon->et());
      if(photon->isEE()) NHIsoSwingValue = NHIsoCut[1] + NHIsoCut[2] * photon->et();
      NHIsoLooseLimit = TMath::Min( 5.0*(NHIsoSwingValue), 0.2*photon->et() );
      NHIsoExclusion = NHIsoSwingValue;
    }

    std::cout<<sigmaIetaIetaCut<<std::endl;

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
    if(MethodID.EqualTo("egamma")){
      if(fabs(photon->superCluster()->eta()) < 1.){effareaCH = 0.012;effareaNH = 0.030;effareaPH = 0.148;}
      if( (fabs(photon->superCluster()->eta()) > 1.) && (fabs(photon->superCluster()->eta()) < 1.479) ){effareaCH = 0.010;effareaNH = 0.057;effareaPH = 0.130;}
      if( (fabs(photon->superCluster()->eta()) > 1.479) && (fabs(photon->superCluster()->eta()) < 2.) ){effareaCH = 0.014;effareaNH = 0.039;effareaPH = 0.112;}
      if( (fabs(photon->superCluster()->eta()) > 2.) && (fabs(photon->superCluster()->eta()) < 2.2) ){effareaCH = 0.012;effareaNH = 0.015;effareaPH = 0.216;}
      if( (fabs(photon->superCluster()->eta()) > 2.2) && (fabs(photon->superCluster()->eta()) < 2.3) ){effareaCH = 0.016;effareaNH = 0.024;effareaPH = 0.262;}
      if( (fabs(photon->superCluster()->eta()) > 2.3) && (fabs(photon->superCluster()->eta()) < 2.4) ){effareaCH = 0.020;effareaNH = 0.039;effareaPH = 0.260;}
      if(fabs(photon->superCluster()->eta()) > 2.4){effareaCH = 0.012;effareaNH = 0.072;effareaPH = 0.266;}
    }

    if(MethodID.EqualTo("highpt")){
      if(fabs(photon->superCluster()->eta()) < 0.9){effareaPH = 0.21;}
      if( (fabs(photon->superCluster()->eta()) > 0.9) && (fabs(photon->superCluster()->eta()) < 1.4442) ){effareaPH = 0.2;}
      if( (fabs(photon->superCluster()->eta()) > 1.560) && (fabs(photon->superCluster()->eta()) < 2.) ){effareaPH = 0.14;}
      if( (fabs(photon->superCluster()->eta()) > 2.) && (fabs(photon->superCluster()->eta()) < 2.2) ){effareaPH = 0.22;}
      if( (fabs(photon->superCluster()->eta()) > 2.2) && (fabs(photon->superCluster()->eta()) < 2.5) ){effareaPH = 0.33;}
    }

    effarea.push_back(effareaCH);
    effarea.push_back(effareaNH);
    effarea.push_back(effareaPH);

    return effarea;
  }

}

#endif

