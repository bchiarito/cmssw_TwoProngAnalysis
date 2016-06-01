#ifndef RECO_TWOPRONG_INFO_INC
#define RECO_TWOPRONG_INFO_INC

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

  struct recoTwoProngInfo_t {
    Double_t CHpos_pt;
    Double_t CHpos_phi;
    Double_t CHpos_eta;
    Double_t CHpos_mass;

    Double_t CHneg_pt;
    Double_t CHneg_phi;
    Double_t CHneg_eta;
    Double_t CHneg_mass;

    Double_t center_pt;
    Double_t center_phi;
    Double_t center_eta;
    Double_t center_mass;

    Double_t photon_pt;
    Double_t photon_phi;
    Double_t photon_eta;
    Double_t photon_mass;

    Double_t Eta_pt;
    Double_t Eta_phi;
    Double_t Eta_eta;
    Double_t Eta_mass;

    Double_t drNearestGen;

    Double_t chargedIso;
    Double_t neutralIso;
    Double_t egammaIso;

    Int_t photon_numGamma;
    Int_t photon_numE;
    Int_t genEtaIndex;
  };

  // also include a string that can be used to define the tree branch
  // obviously this needs to be kept up-to-date with the struct definition
  // but now at least this only needs to be done here in this file, 
  // rather than in each individual analyser 
  std::string recoTwoProngBranchDefString("CHpos_pt/D:CHpos_phi/D:CHpos_eta/D:CHpos_mass/D:CHneg_pt/D:CHneg_phi/D:CHneg_eta/D:CHneg_mass/D:center_pt/D:center_phi/D:center_eta/D:center_mass/D:photon_pt/D:photon_phi/D:photon_eta/D:photon_mass/D:Eta_pt/D:Eta_phi/D:Eta_eta/D:Eta_mass/D:drNearestGen/D:chargedIso/D:neutralIso/D:egammaIso/D:photon_numGamma/I:photon_numE/I:genEtaIndex/I");

  // also want a Fill function, that can fill the struct values from the appropriate objects
  // again, so that all editing only needs to be done here in this file
  void FillRecoTwoProngInfo(recoTwoProngInfo_t &recotwopronginfo, TLorentzVector &CHpos, TLorentzVector &CHneg, TLorentzVector &center,
                            TLorentzVector &photon, TLorentzVector &Eta, bool passed, bool matched, double genDR, int genIndex, 
                            double relchargediso, double relneutraliso, double relegammaiso, int numgamma, int nume) {

    recotwopronginfo.CHpos_pt = CHpos.Pt();
    recotwopronginfo.CHpos_phi = CHpos.Phi();
    recotwopronginfo.CHpos_eta = CHpos.Eta();
    recotwopronginfo.CHpos_mass = CHpos.M();

    recotwopronginfo.CHneg_pt = CHneg.Pt();
    recotwopronginfo.CHneg_phi = CHneg.Phi();
    recotwopronginfo.CHneg_eta = CHneg.Eta();
    recotwopronginfo.CHneg_mass = CHneg.M();

    recotwopronginfo.center_pt = center.Pt();
    recotwopronginfo.center_phi = center.Phi();
    recotwopronginfo.center_eta = center.Eta();
    recotwopronginfo.center_mass = center.M();

    recotwopronginfo.photon_pt = photon.Pt();
    recotwopronginfo.photon_phi = photon.Phi();
    recotwopronginfo.photon_eta = photon.Eta();
    recotwopronginfo.photon_mass = photon.M();

    recotwopronginfo.Eta_pt = Eta.Pt();
    recotwopronginfo.Eta_phi = Eta.Phi();
    recotwopronginfo.Eta_eta = Eta.Eta();
    recotwopronginfo.Eta_mass = Eta.M();

    recotwopronginfo.drNearestGen = genDR;

    recotwopronginfo.chargedIso = relchargediso;
    recotwopronginfo.neutralIso = relneutraliso;
    recotwopronginfo.egammaIso = relegammaiso;

    recotwopronginfo.photon_numGamma = numgamma;
    recotwopronginfo.photon_numE = nume;
    recotwopronginfo.genEtaIndex = genIndex;
  }

  void InitRecoTwoProngInfo(recoTwoProngInfo_t &recotwopronginfo) {
    recotwopronginfo.CHpos_pt = -9.9;
    recotwopronginfo.CHpos_phi = -9.9;
    recotwopronginfo.CHpos_eta = -9.9;
    recotwopronginfo.CHpos_mass = -9.9;

    recotwopronginfo.CHneg_pt = -9.9;
    recotwopronginfo.CHneg_phi = -9.9;
    recotwopronginfo.CHneg_eta = -9.9;
    recotwopronginfo.CHneg_mass = -9.9;

    recotwopronginfo.center_pt = -9.9;
    recotwopronginfo.center_phi = -9.9;
    recotwopronginfo.center_eta = -9.9;
    recotwopronginfo.center_mass = -9.9;

    recotwopronginfo.photon_pt = -9.9;
    recotwopronginfo.photon_phi = -9.9;
    recotwopronginfo.photon_eta = -9.9;
    recotwopronginfo.photon_mass = -9.9;

    recotwopronginfo.Eta_pt = -9.9;
    recotwopronginfo.Eta_phi = -9.9;
    recotwopronginfo.Eta_eta = -9.9;
    recotwopronginfo.Eta_mass = -9.9;

    recotwopronginfo.drNearestGen = -9.9;

    recotwopronginfo.chargedIso = -9.9;
    recotwopronginfo.neutralIso = -9.9;
    recotwopronginfo.egammaIso = -9.9;

    recotwopronginfo.photon_numGamma = -1;
    recotwopronginfo.photon_numE = -1;
    recotwopronginfo.genEtaIndex = -1;
  }

} //end of namespace

#endif
