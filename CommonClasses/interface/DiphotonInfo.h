#ifndef DIPHOTON_INFO_INC
#define DIPHOTON_INFO_INC

//***********************************************
// Diphoton Info
//
//************************************************

#include <string>

#include "DataFormats/Math/interface/deltaPhi.h"
#include "TMath.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/Candidate/interface/LeafCandidate.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TLorentzVector.h"

namespace ExoDiPhotons{


  // diphoton info: Minv, q_T, delta phi, etc
  
  struct diphotonInfo_t{
    double Minv;
    double qt;
    double deltaPhi;
    double deltaEta;
    double deltaR;
    double deltaROld;
    double cosThetaStar;
    double cosThetaStarOld;
  };

  std::string diphotonInfoBranchDefString("Minv/D:qt:deltaPhi:deltaEta:deltaR:deltaROld:cosThetaStar:cosThetaStarOld");
  

  // the internal function which gets called by the others
  void FillDiphotonInfo(diphotonInfo_t &fDiphotonInfo, reco::LeafCandidate::LorentzVector photon_vector1, reco::LeafCandidate::LorentzVector photon_vector2) {  
    fDiphotonInfo.Minv = (photon_vector1+photon_vector2).M();
    // pt of the pair
    fDiphotonInfo.qt = (photon_vector1+photon_vector2).pt();
    // there's a CMS function for deltaPhi in DataFormats/Math
    fDiphotonInfo.deltaPhi = reco::deltaPhi(photon_vector1.phi(),photon_vector2.phi());
    fDiphotonInfo.deltaEta = photon_vector1.eta()-photon_vector2.eta(); // always highest pt - second
    // use CMS function for deltaR
    fDiphotonInfo.deltaR = reco::deltaR(photon_vector1.eta(),photon_vector1.phi(),photon_vector2.eta(),photon_vector2.phi());
    // old deltaR - leave in to check
    fDiphotonInfo.deltaROld = TMath::Sqrt(fDiphotonInfo.deltaPhi*fDiphotonInfo.deltaPhi+fDiphotonInfo.deltaEta*fDiphotonInfo.deltaEta);
    //    reco::LeafCandidate::LorentzVector diphoton = photon_vector1 + photon_vector2;
    //    double mymom = diphoton.P();

    // boost photon 1 into diphoton rest frame
    // the Boost() function changes the acutal Lorentz vector that its called
    // first lets boost it, then call CosTheta() after boosting

    //reco::LeafCandidate::LorentzVector photon1_clone = photon_vector1.Clone();
    TLorentzVector photon1_clone(photon_vector1.px(),photon_vector1.py(),photon_vector1.pz(),photon_vector1.E());
    TLorentzVector photon2_clone(photon_vector2.px(),photon_vector2.py(),photon_vector2.pz(),photon_vector2.E());
    TLorentzVector diphoton_system = photon1_clone+photon2_clone;
    photon1_clone.Boost(-diphoton_system.BoostVector());
    //photon1_clone.Boost(-(photon_vector1+photon_vector2).BoostVector());
    fDiphotonInfo.cosThetaStar = photon1_clone.CosTheta();
    // this older implementation was provided by Yousi Ma - should check it sometime
    fDiphotonInfo.cosThetaStarOld = fabs(photon_vector1.P() - photon_vector2.P())/(photon_vector1+photon_vector2).P();
  }

  // filling function for reco photons  
  void FillDiphotonInfo(diphotonInfo_t &fDiphotonInfo, const reco::Photon *photon1, const reco::Photon *photon2) {
    
    reco::LeafCandidate::LorentzVector photon_vector1 = photon1->p4();
    reco::LeafCandidate::LorentzVector photon_vector2 = photon2->p4();
    
    FillDiphotonInfo(fDiphotonInfo,photon_vector1,photon_vector2);
  }

  // same function, but for MC signal photons
  void FillDiphotonInfo(diphotonInfo_t &fDiphotonInfo, const reco::GenParticle *photon1, const reco::GenParticle *photon2) {
    
    reco::LeafCandidate::LorentzVector photon_vector1 = photon1->p4();
    reco::LeafCandidate::LorentzVector photon_vector2 = photon2->p4();
    
    FillDiphotonInfo(fDiphotonInfo,photon_vector1,photon_vector2);  
  }
  
  // same function, but with reco::Candidates
  void FillDiphotonInfo(diphotonInfo_t &fDiphotonInfo, const reco::Candidate *photon1, const reco::Candidate *photon2) {
    
    reco::LeafCandidate::LorentzVector photon_vector1 = photon1->p4();
    reco::LeafCandidate::LorentzVector photon_vector2 = photon2->p4();
    
    FillDiphotonInfo(fDiphotonInfo,photon_vector1,photon_vector2);
    
  }
  
  void InitDiphotonInfo(diphotonInfo_t &diphotonInfo) {
    diphotonInfo.Minv = -99999.99;
    diphotonInfo.qt = -99999.99;
    diphotonInfo.deltaPhi = -99999.99;
    diphotonInfo.deltaEta = -99999.99;
    diphotonInfo.deltaR = -99999.99;
    diphotonInfo.deltaROld = -99999.99;
    diphotonInfo.cosThetaStar = -99999.99;
    diphotonInfo.cosThetaStarOld = -99999.99;
  }

} //end of namespace


#endif
