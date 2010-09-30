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

namespace ExoDiPhotons{


  // diphoton info: Minv, q_T, delta phi, etc

  struct diphotonInfo_t{
    double Minv;
    double qt;
    double deltaPhi;
    double deltaEta;
    double deltaR;
    double cosThetaStar;
    
  };

  std::string diphotonInfoBranchDefString("Minv/D:qt:deltaPhi:deltaEta:deltaR:cosThetaStar");
  

    // the internal function which gets called by the others
  void FillDiphotonInfo(diphotonInfo_t &fDiphotonInfo, reco::LeafCandidate::LorentzVector photon_vector1, reco::LeafCandidate::LorentzVector photon_vector2) {

    fDiphotonInfo.Minv = (photon_vector1+photon_vector2).M();
    // pt of the pair
    fDiphotonInfo.qt = (photon_vector1+photon_vector2).pt();
    // there's a CMS function for deltaPhi in DataFormats/Math
    fDiphotonInfo.deltaPhi = reco::deltaPhi(photon_vector1.phi(),photon_vector2.phi());
    fDiphotonInfo.deltaEta = photon_vector1.eta()-photon_vector2.eta(); // always highest pt - second
    fDiphotonInfo.deltaR = TMath::Sqrt(fDiphotonInfo.deltaPhi*fDiphotonInfo.deltaPhi+fDiphotonInfo.deltaEta*fDiphotonInfo.deltaEta);
    fDiphotonInfo.cosThetaStar = -999.99; // need to calculate this!


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


  


} //end of namespace


#endif
