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

  struct conversionInfo_t{

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> r;
    std::vector<double> phi;
    std::vector<double> dPhiTracksAtVtx;
    std::vector<double> nTracks;
    std::vector<double> dxy;
    std::vector<double> dz;
    std::vector<double> vtxChi2;
    std::vector<double> vtxNdof;
    std::vector<double> pairCotThetaSeparation;
    std::vector<double> photonPt;
    std::vector<double> dRToSc;
    std::vector<uint8_t> nSharedHits;
    std::vector<double> MVAout;
    std::vector<std::vector<float>> oneLegMVA;
    std::vector<std::vector<uint8_t>> nHitsBeforeVtx;
    std::vector<std::vector<int>> quality;

    // add dxy, dz, pairCotThetaSeparation , r = sqrt(x**2+y**2), phi also, link back to original photon

  };

  // std::string conversionInfoBranchDefString("nConversions/I:x[nConversions]/F:y[nConversions]:z[nConversions]:r[nConversions]:phi[nConversions]:dPhiTracksAtVtx[nConversions]:nTracks[nConversions]:dxy[nConversions]:dz[nConversions]:pairCotThetaSeparation[nConversions]:photonPt[nConversions]:dRToSc[nConversions]");
  
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
    
    // first need to clear vectors from the previous event!
    convInfo.x.clear();
    convInfo.y.clear();
    convInfo.z.clear();
    convInfo.r.clear();
    convInfo.phi.clear();
    convInfo.dPhiTracksAtVtx.clear();
    convInfo.nTracks.clear();
    convInfo.dxy.clear();
    convInfo.dz.clear();
    convInfo.vtxChi2.clear();
    convInfo.vtxNdof.clear();
    convInfo.pairCotThetaSeparation.clear();
    convInfo.photonPt.clear();
    convInfo.dRToSc.clear();
    convInfo.nSharedHits.clear();
    convInfo.MVAout.clear();
    convInfo.oneLegMVA.clear();
    convInfo.nHitsBeforeVtx.clear();
    convInfo.quality.clear();

    int numMatched = 0;
    for (unsigned int i=0; i<convColl.size(); i++){
      reco::Conversion iConv = convColl.at(i);
      bool matched = ConversionTools::matchesConversion(sc,iConv,0.2,999.,999.); // matched within dR of 0.2
      if (matched){
        convInfo.x.push_back(iConv.conversionVertex().position().X());
        convInfo.y.push_back(iConv.conversionVertex().position().Y());
        convInfo.z.push_back(iConv.conversionVertex().position().Z());
        convInfo.r.push_back(TMath::Sqrt((convInfo.x[numMatched]*convInfo.x[numMatched]) + (convInfo.y[numMatched]*convInfo.y[numMatched]))); //numMatched is incremented at the end, so using it as the index is OK
        convInfo.phi.push_back(TMath::ATan2( convInfo.y[numMatched] , convInfo.x[numMatched] ));
        convInfo.dPhiTracksAtVtx.push_back(iConv.dPhiTracksAtVtx());
        convInfo.nTracks.push_back(iConv.nTracks());
        convInfo.dxy.push_back(iConv.dxy(bs.position()));
        convInfo.dz.push_back(iConv.dz(bs.position()));
        convInfo.vtxChi2.push_back(iConv.conversionVertex().chi2());
        convInfo.vtxNdof.push_back(iConv.conversionVertex().ndof());
        convInfo.pairCotThetaSeparation.push_back(iConv.pairCotThetaSeparation());
        convInfo.photonPt.push_back(phoPt);
        convInfo.dRToSc.push_back(scConvDr(sc,iConv));
        convInfo.nSharedHits.push_back(iConv.nSharedHits());
        convInfo.MVAout.push_back(iConv.MVAout());
        convInfo.oneLegMVA.push_back(iConv.oneLegMVA());
        convInfo.nHitsBeforeVtx.push_back(iConv.nHitsBeforeVtx());

        // build quality vector
        std::vector<int> qualityVec;
        qualityVec.push_back(iConv.quality(reco::Conversion::ConversionQuality::generalTracksOnly));
        qualityVec.push_back(iConv.quality(reco::Conversion::ConversionQuality::arbitratedEcalSeeded));
        qualityVec.push_back(iConv.quality(reco::Conversion::ConversionQuality::arbitratedMerged));
        qualityVec.push_back(iConv.quality(reco::Conversion::ConversionQuality::arbitratedMergedEcalGeneral));
        qualityVec.push_back(iConv.quality(reco::Conversion::ConversionQuality::highPurity));
        qualityVec.push_back(iConv.quality(reco::Conversion::ConversionQuality::highEfficiency));
        qualityVec.push_back(iConv.quality(reco::Conversion::ConversionQuality::ecalMatched1Track));
        qualityVec.push_back(iConv.quality(reco::Conversion::ConversionQuality::ecalMatched2Track));

        convInfo.quality.push_back(qualityVec);
        numMatched++;

      }
    } // end loop over conversion collection


  }

} //end of namespace


#endif
