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
    Double_t dRToSc[maxConversions];
    Bool_t isConverted[maxConversions];

    // add dxy, dz, pairCotThetaSeparation , r = sqrt(x**2+y**2), phi also, link back to original photon

  };

  // std::string jetInfoBranchDefString("Minv/D:qt:deltaPhi:deltaEta:deltaR:deltaROld:cosThetaStar:cosThetaStarOld");
  // std::string jetInfoBranchDefString("nJets/I");
  // std::string jetInfoBranchDefString("nJets/I:pt[nJets]/D");
  
  // return dR between supercluster and conversion
  double scConvDr(const reco::SuperCluster &sc, const reco::Conversion &conv){

    // copied from RecoEgamma/EgammaTools/src/ConversionTools.cc

    math::XYZVector mom(conv.refittedPairMomentum());

    math::XYZPoint scpos(sc.position());
    math::XYZPoint cvtx(conv.conversionVertex().position());


    math::XYZVector cscvector = scpos - cvtx;
    double dR = (double) reco::deltaR(mom,cscvector);
    // float dEta = mom.eta() - cscvector.eta();
    // float dPhi = reco::deltaPhi(mom.phi(),cscvector.phi());

    // if (dR>dRMax) return false;
    // if (dEta>dEtaMax) return false;
    // if (dPhi>dPhiMax) return false;

    return dR;


  }
  void FillConversionInfo(conversionInfo_t &convInfo, const reco::SuperCluster& sc, const reco::ConversionCollection& convColl, double phoPt, const reco::BeamSpot& bs) {
    using namespace std;
    cout << "in FillConversionInfo" << endl;
    // convInfo.nConversions = (Int_t)convColl.size();
    int numMatched = 0;
    for (unsigned int i=0; i<convColl.size(); i++){
      cout << "in conversion collection loop 0" << endl;
      reco::Conversion iConv = convColl.at(i);
      bool matched = ConversionTools::matchesConversion(sc,iConv,0.2,999.,999.);
      if (matched){
        cout << "in conversion collection loop 1" << endl;
        convInfo.x[numMatched] = iConv.conversionVertex().position().X(); //numMatched is incremented at the end, so using it as the index is OK
        convInfo.y[numMatched] = iConv.conversionVertex().position().Y();
        convInfo.z[numMatched] = iConv.conversionVertex().position().Z();
        cout << "in conversion collection loop 2" << endl;
        convInfo.r[numMatched] = TMath::Sqrt((convInfo.x[numMatched]*convInfo.x[numMatched]) + (convInfo.y[numMatched]*convInfo.y[numMatched]));
        convInfo.phi[numMatched] = TMath::ATan2( convInfo.y[numMatched] , convInfo.x[numMatched] );
        cout << "in conversion collection loop 3" << endl;
        convInfo.dPhiTracksAtVtx[numMatched] = iConv.dPhiTracksAtVtx();
        cout << "in conversion collection loop 3.1" << endl;
        // convInfo.dPhiTracksAtEcal[numMatched] = iConv.dPhiTracksAtEcal();
        // cout << "in conversion collection loop 3.2" << endl;
        // convInfo.dEtaTracksAtEcal[numMatched] = iConv.dEtaTracksAtEcal();
        // cout << "in conversion collection loop 4" << endl;
        convInfo.dPhiTracksAtEcal[numMatched] = -100.;
        cout << "in conversion collection loop 3.2" << endl;
        convInfo.dEtaTracksAtEcal[numMatched] = -100.;
        convInfo.nTracks[numMatched] = 1.*iConv.nTracks();
        convInfo.dxy[numMatched] = iConv.dxy(bs.position());
        convInfo.dz[numMatched] = iConv.dz(bs.position());
        cout << "in conversion collection loop 5" << endl;
        convInfo.pairCotThetaSeparation[numMatched] = iConv.pairCotThetaSeparation();
        cout << "in conversion collection loop 6" << endl;
        convInfo.photonPt[numMatched] = phoPt;
        convInfo.dRToSc[numMatched] = scConvDr(sc,iConv);
        convInfo.isConverted[numMatched] = iConv.isConverted();
        numMatched++;
      }
      convInfo.nConversions = numMatched;
      cout << "finish conversion collection loop" << endl;
    } // end loop over conversion collection


  }

} //end of namespace


#endif
