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

  // std::string jetInfoBranchDefString("Minv/D:qt:deltaPhi:deltaEta:deltaR:deltaROld:cosThetaStar:cosThetaStarOld");
  // std::string jetInfoBranchDefString("nJets/I");
  // std::string jetInfoBranchDefString("nJets/I:pt[nJets]/D");
  std::string conversionInfoBranchDefString("nConversions/I:x[nConversions]/F:y[nConversions]:z[nConversions]:r[nConversions]:phi[nConversions]:dPhiTracksAtVtx[nConversions]:nTracks[nConversions]:dxy[nConversions]:dz[nConversions]:pairCotThetaSeparation[nConversions]:photonPt[nConversions]:dRToSc[nConversions]");
  
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
      bool matched = ConversionTools::matchesConversion(sc,iConv,0.2,999.,999.); // matched within dR of 0.2
      if (matched){
        cout << "in conversion collection loop 1" << endl;
        // convInfo.x[numMatched] = iConv.conversionVertex().position().X(); //numMatched is incremented at the end, so using it as the index is OK
        // convInfo.y[numMatched] = iConv.conversionVertex().position().Y();
        convInfo.x.push_back(iConv.conversionVertex().position().X());
        convInfo.y.push_back(iConv.conversionVertex().position().Y());
        convInfo.z.push_back(iConv.conversionVertex().position().Z());
        cout << "in conversion collection loop 2" << endl;
        convInfo.r.push_back(TMath::Sqrt((convInfo.x[numMatched]*convInfo.x[numMatched]) + (convInfo.y[numMatched]*convInfo.y[numMatched])));
        convInfo.phi.push_back(TMath::ATan2( convInfo.y[numMatched] , convInfo.x[numMatched] ));
        cout << "in conversion collection loop 3" << endl;
        convInfo.dPhiTracksAtVtx.push_back(iConv.dPhiTracksAtVtx());
        cout << "in conversion collection loop 3.1" << endl;
        convInfo.nTracks.push_back(iConv.nTracks());
        convInfo.dxy.push_back(iConv.dxy(bs.position()));
        convInfo.dz.push_back(iConv.dz(bs.position()));
        cout << "in conversion collection loop 5" << endl;
        convInfo.pairCotThetaSeparation.push_back(iConv.pairCotThetaSeparation());
        cout << "in conversion collection loop 6" << endl;
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
      // convInfo.nConversions = numMatched;
      cout << "finish conversion collection loop" << endl;
    } // end loop over conversion collection


  }

} //end of namespace


#endif
