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
    std::vector<double> track1InnerPx;
    std::vector<double> track1InnerPy;
    std::vector<double> track1InnerPz;
    std::vector<double> track2InnerPx;
    std::vector<double> track2InnerPy;
    std::vector<double> track2InnerPz;
    std::vector<double> dRToSc;
    // std::vector<int> nHitsTrack1;
    // std::vector<int> nHitsTrack2;
    std::vector<uint8_t> nSharedHits;
    std::vector<double> MVAout;
    std::vector<std::vector<float>> oneLegMVA;
    std::vector<std::vector<uint8_t>> nHitsBeforeVtx;
    std::vector<int> isGeneralTracksOnly;
    std::vector<int> isArbitratedEcalSeeded;
    std::vector<int> isArbitratedMerged;
    std::vector<int> isArbitratedMergedEcalGeneral;
    std::vector<int> isHighPurity;
    std::vector<int> isHighEfficiency;
    std::vector<int> isEcalMatched1Track;
    std::vector<int> isEcalMatched2Track;

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
    convInfo.track1InnerPx.clear();
    convInfo.track1InnerPy.clear();
    convInfo.track1InnerPz.clear();
    convInfo.track2InnerPx.clear();
    convInfo.track2InnerPy.clear();
    convInfo.track2InnerPz.clear();
    convInfo.dRToSc.clear();
    // convInfo.nHitsTrack1.clear();
    // convInfo.nHitsTrack2.clear();
    convInfo.nSharedHits.clear();
    convInfo.MVAout.clear();
    convInfo.oneLegMVA.clear();
    convInfo.nHitsBeforeVtx.clear();
    convInfo.isGeneralTracksOnly.clear();
    convInfo.isArbitratedEcalSeeded.clear();
    convInfo.isArbitratedMerged.clear();
    convInfo.isArbitratedMergedEcalGeneral.clear();
    convInfo.isHighPurity.clear();
    convInfo.isHighEfficiency.clear();
    convInfo.isEcalMatched1Track.clear();
    convInfo.isEcalMatched2Track.clear();

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

        // quality flags
        convInfo.isGeneralTracksOnly.push_back(iConv.quality(reco::Conversion::ConversionQuality::generalTracksOnly));
        convInfo.isArbitratedEcalSeeded.push_back(iConv.quality(reco::Conversion::ConversionQuality::arbitratedEcalSeeded));
        convInfo.isArbitratedMerged.push_back(iConv.quality(reco::Conversion::ConversionQuality::arbitratedMerged));
        convInfo.isArbitratedMergedEcalGeneral.push_back(iConv.quality(reco::Conversion::ConversionQuality::arbitratedMergedEcalGeneral));
        convInfo.isHighPurity.push_back(iConv.quality(reco::Conversion::ConversionQuality::highPurity));
        convInfo.isHighEfficiency.push_back(iConv.quality(reco::Conversion::ConversionQuality::highEfficiency));
        convInfo.isEcalMatched1Track.push_back(iConv.quality(reco::Conversion::ConversionQuality::ecalMatched1Track));
        convInfo.isEcalMatched2Track.push_back(iConv.quality(reco::Conversion::ConversionQuality::ecalMatched2Track));

        //conversion momentum from innermost track
        const std::vector<math::XYZVectorF> momvecs = iConv.tracksPin();

        math::XYZVectorF track1mom = momvecs.at(0);
        convInfo.track1InnerPx.push_back( track1mom.X() );
        convInfo.track1InnerPy.push_back( track1mom.Y() );
        convInfo.track1InnerPz.push_back( track1mom.Z() );

        if (iConv.nTracks() == 2){ // track2 momentum vectors will be empty for one leg conversions
            math::XYZVectorF track2mom = momvecs.at(1);
            convInfo.track2InnerPx.push_back( track2mom.X() );
            convInfo.track2InnerPy.push_back( track2mom.Y() );
            convInfo.track2InnerPz.push_back( track2mom.Z() );
        }

        // fill track nHits
        // const reco::Track* track1 = iConv.tracks().at(0).get();
        // int nHits1 = 0;
        // int nHits2 = -99;
        // for (size_t i=0; i<track1->recHitsSize(); i++){
        //     TrackingRecHitRef iHit = track1->recHit(i);
        //     if ( iHit->isValid() ) nHits1++;
        // }
        // if (iConv.nTracks() == 2){ // i.e. only do this for two pronged conversions, nHits2=-99 for one pronged
        //     const reco::Track* track2 = iConv.tracks().at(1).get();
        //     nHits2 = 0;
        //     for (size_t i=0; i<track2->recHitsSize(); i++){
        //         TrackingRecHitRef iHit = track2->recHit(i);
        //         if ( iHit->isValid() ) nHits2++;
        //     }
        // }
        // convInfo.nHitsTrack1.push_back(nHits1);
        // convInfo.nHitsTrack2.push_back(nHits2);

        numMatched++;

      }
    } // end loop over conversion collection


  }

} //end of namespace


#endif
