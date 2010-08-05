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

  void FillDiphotonInfo(diphotonInfo_t &fDiphotonInfo, const reco::Photon *photon1, const reco::Photon *photon2) {

    reco::LeafCandidate::LorentzVector photon_vector1 = photon1->p4();
    reco::LeafCandidate::LorentzVector photon_vector2 = photon2->p4();

    fDiphotonInfo.Minv = (photon_vector1+photon_vector2).M();
    // pt of the pair
    fDiphotonInfo.qt = (photon_vector1+photon_vector2).pt();
    // there's a CMS function for deltaPhi in DataFormats/Math
    fDiphotonInfo.deltaPhi = reco::deltaPhi(photon1->phi(),photon2->phi());
    fDiphotonInfo.deltaEta = photon1->eta()-photon2->eta(); // always highest pt - second
    fDiphotonInfo.deltaR = TMath::Sqrt(fDiphotonInfo.deltaPhi*fDiphotonInfo.deltaPhi+fDiphotonInfo.deltaEta*fDiphotonInfo.deltaEta);
    fDiphotonInfo.cosThetaStar = -999.99; // need to calculate this!
    
  }



} //end of namespace


#endif
