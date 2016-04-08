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
  

  void FillConversionInfo(conversionInfo_t &convInfo, const reco::ConversionRefVector& convColl, double phoPt, const reco::BeamSpot& bs) {
    using namespace std;
    cout << "in FillConversionInfo" << endl;
    convInfo.nConversions = (Int_t)convColl.size();
    for (int i=0; i<convInfo.nConversions; i++){
      cout << "in conversion collection loop 0" << endl;
      reco::Conversion iConv = *(convColl.at(i));
      cout << "in conversion collection loop 1" << endl;
      convInfo.x[i] = iConv.conversionVertex().position().X();
      convInfo.y[i] = iConv.conversionVertex().position().Y();
      convInfo.z[i] = iConv.conversionVertex().position().Z();
      cout << "in conversion collection loop 2" << endl;
      convInfo.r[i] = TMath::Sqrt((convInfo.x[i]*convInfo.x[i]) + (convInfo.y[i]*convInfo.y[i]));
      convInfo.phi[i] = TMath::ATan2( convInfo.y[i] , convInfo.x[i] );
      cout << "in conversion collection loop 3" << endl;
      convInfo.dPhiTracksAtVtx[i] = iConv.dPhiTracksAtVtx();
      convInfo.dPhiTracksAtEcal[i] = iConv.dPhiTracksAtEcal();
      convInfo.dEtaTracksAtEcal[i] = iConv.dEtaTracksAtEcal();
      cout << "in conversion collection loop 4" << endl;
      convInfo.nTracks[i] = 1.*iConv.nTracks();
      convInfo.dxy[i] = iConv.dxy(bs.position());
      convInfo.dz[i] = iConv.dz(bs.position());
      cout << "in conversion collection loop 5" << endl;
      convInfo.pairCotThetaSeparation[i] = iConv.pairCotThetaSeparation();
      cout << "in conversion collection loop 6" << endl;
      convInfo.photonPt[i] = phoPt;
      convInfo.isConverted[i] = iConv.isConverted();
      cout << "finish conversion collection loop" << endl;
    } // end loop over conversion collection


  }

} //end of namespace


#endif
