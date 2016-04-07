#ifndef CONVERSION_INFO_INC
#define CONVERSION_INFO_INC

//***********************************************
// Diphoton Info
//
//************************************************

#include <string>
#include <iostream>

#include "DataFormats/Math/interface/deltaPhi.h"
#include "TMath.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/Candidate/interface/LeafCandidate.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TLorentzVector.h"

// for conversions
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

namespace ExoDiPhotons{


  // diphoton info: Minv, q_T, delta phi, etc
  const Int_t maxConversions = 500;
  struct conversionInfo_t{

    Int_t nConversions;
    Double_t x[maxConversions];
    Double_t y[maxConversions];
    Double_t z[maxConversions];
    Double_t r[maxConversions];
    Double_t phi[maxConversions];
    Double_t dPhiTracksAtVtx[maxConversions];
    Double_t dPhiTracksAtEcal[maxConversions];
    Double_t dEtaTracksAtEcal[maxConversions];
    Double_t nTracks[maxConversions];
    Double_t dxy[maxConversions];
    Double_t dz[maxConversions];
    Double_t pairCotThetaSeparation[maxConversions];
    Double_t photonPt[maxConversions];
    Bool_t isConverted[maxConversions];

    // add dxy, dz, pairCotThetaSeparation , r = sqrt(x**2+y**2), phi also, link back to original photon

  };

  // std::string jetInfoBranchDefString("Minv/D:qt:deltaPhi:deltaEta:deltaR:deltaROld:cosThetaStar:cosThetaStarOld");
  // std::string jetInfoBranchDefString("nJets/I");
  // std::string jetInfoBranchDefString("nJets/I:pt[nJets]/D");
  

  void FillConversionInfo(conversionInfo_t &convInfo, const reco::ConversionRefVector convColl, double phoPt) {

    convInfo.nConversions = (Int_t)convColl.size();
    for (int i=0; i<convColl.nConversions; i++){
      reco::Conversion iConv = *(convsRefVec.at(i));
      x[i] = iConv.conversionVertex().position().X();
      y[i] = iConv.conversionVertex().position().Y();
      z[i] = iConv.conversionVertex().position().Z();
      dPhiTracksAtVtx[i] = iConv.dPhiTracksAtVtx();
      dPhiTracksAtEcal[i] = iConv.dPhiTracksAtEcal();
      dEtaTracksAtEcal[i] = iConv.dEtaTracksAtEcal();
      nTracks[i] = 1.*iConv.nTracks();
      isConverted[i] = iConv.isConverted();
      // std::cout << "jet pt " << jets->at(i).pt() << std::endl;
      pat::Jet jet = jets->at(i);
      convInfo.pt[i] = jet.pt();
      convInfo.eta[i] = jet.eta();
      convInfo.phi[i] = jet.phi();
      convInfo.mass[i] = jet.mass();
      convInfo.energy[i] = jet.energy();
      
      Bool_t loose;
      Bool_t tight;
      std::tie(loose,tight) = ExoDiPhotons::jetID(jet);
      convInfo.passLooseID[i] = loose;
      convInfo.passTightID[i] = tight;
    }


  }
  // void InitJetInfo(jetInfo_t &convInfo, const edm::View<pat::Jet>* jets) {

  //   convInfo.nJets = (Int_t)jets->size();
  //   // convInfo.pt = *(new Double_t [convInfo.nJets]);
    
  // }

  

} //end of namespace


#endif
