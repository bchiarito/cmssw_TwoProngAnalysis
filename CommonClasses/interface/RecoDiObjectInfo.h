#ifndef RECO_DIOBJECT_INFO_INC
#define RECO_DIOBJECT_INFO_INC

//********************************************************************
// Definition of a struct that can be used for storing reco charged type decaying eta info
// in a tree, from different analysers
// Also includes a Fill function to fill the struct from the appropriate objects
// and a string that can be used to define the tree branch
// 
//  $Id: RecoTwoProngInfo.h,v 1.00 2016 16:26:48 charaf Exp $
// 
//********************************************************************

#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

// geometry
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

namespace ExoDiPhotons
{

  struct recoDiObjectInfo_t {
    Double_t pt;
    Double_t phi;
    Double_t eta;
    Double_t mass;
    Double_t px;
    Double_t py;
    Double_t pz;
    Double_t energy;

    Double_t dR;
    Double_t dPt;
    Double_t dPhi;
    Double_t dEta;
    Double_t dMass;
  };

  // also include a string that can be used to define the tree branch
  // obviously this needs to be kept up-to-date with the struct definition
  // but now at least this only needs to be done here in this file, 
  // rather than in each individual analyser 
  std::string recoDiObjectBranchDefString("pt/D:phi/D:eta/D:mass/D:px/D:py/D:pz/D:energy/D:dR/D:dR/D:dPt/D:dPhi/D:dEta/D:dMass/D");

  // also want a Fill function, that can fill the struct values from the appropriate objects
  // again, so that all editing only needs to be done here in this file
  void FillRecoDiObjectInfo(recoDiObjectInfo_t &recodiobjectinfo, TLorentzVector vec1, TLorentzVector vec2)
  {
    TLorentzVector comb = vec1 + vec2;
    recodiobjectinfo.pt = comb.Pt();
    recodiobjectinfo.phi = comb.Phi();
    recodiobjectinfo.eta = comb.Eta();
    recodiobjectinfo.mass = comb.M();
    recodiobjectinfo.px = comb.Px();
    recodiobjectinfo.py = comb.Py();
    recodiobjectinfo.pz = comb.Pz();
    recodiobjectinfo.energy = comb.E();

    recodiobjectinfo.dR = vec1.DeltaR(vec2);
    recodiobjectinfo.dPt = abs(vec1.Pt() - vec2.Pt());
    recodiobjectinfo.dPhi = vec1.DeltaPhi(vec2);
    recodiobjectinfo.dEta = abs(vec1.Eta() - vec2.Eta());
    recodiobjectinfo.dMass = abs(vec1.M() - vec2.M());
  }

  void InitRecoDiObjectInfo(recoDiObjectInfo_t &recodiobjectinfo)
  {
    recodiobjectinfo.pt = -99.9;
    recodiobjectinfo.phi = -99.9;
    recodiobjectinfo.eta = -99.9;
    recodiobjectinfo.mass = -99.9;
    recodiobjectinfo.px = -99.9;
    recodiobjectinfo.py = -99.9;
    recodiobjectinfo.pz = -99.9;
    recodiobjectinfo.energy = -99.9;

    recodiobjectinfo.dR = -99.9;
    recodiobjectinfo.dPt = -99.9;
    recodiobjectinfo.dPhi = -99.9;
    recodiobjectinfo.dEta = -99.9;
    recodiobjectinfo.dMass = -99.9;
  }

} //end of namespace

#endif
