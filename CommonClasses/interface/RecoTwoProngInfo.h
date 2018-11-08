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

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

namespace TwoProngAnalysis
{

  struct recoTwoProngInfo_t {
    Double_t CHpos_pt;
    Double_t CHpos_phi;
    Double_t CHpos_eta;
    Double_t CHpos_mass;
    Double_t CHpos_px;
    Double_t CHpos_py;
    Double_t CHpos_pz;
    Double_t CHpos_energy;

    Double_t CHneg_pt;
    Double_t CHneg_phi;
    Double_t CHneg_eta;
    Double_t CHneg_mass;
    Double_t CHneg_px;
    Double_t CHneg_py;
    Double_t CHneg_pz;
    Double_t CHneg_energy;

    Double_t center_pt;
    Double_t center_phi;
    Double_t center_eta;
    Double_t center_mass;
    Double_t center_px;
    Double_t center_py;
    Double_t center_pz;
    Double_t center_energy;

    Double_t photon_pt;
    Double_t photon_phi;
    Double_t photon_eta;
    Double_t photon_mass;
    Double_t photon_px;
    Double_t photon_py;
    Double_t photon_pz;
    Double_t photon_energy;

    Double_t pt;
    Double_t phi;
    Double_t eta;
    Double_t mass;
    Double_t px;
    Double_t py;
    Double_t pz;
    Double_t energy;
    Double_t grommedMass;

    Double_t genDR;

    Double_t chargedIso;
    Double_t neutralIso;
    Double_t egammaIso;

    Int_t photon_nPhotons;
    Int_t photon_nElectrons;
    Int_t genEtaIndex;
  };

  // also include a string that can be used to define the tree branch
  // obviously this needs to be kept up-to-date with the struct definition
  // but now at least this only needs to be done here in this file, 
  // rather than in each individual analyser 
  std::string recoTwoProngBranchDefString("CHpos_pt/D:CHpos_phi/D:CHpos_eta/D:CHpos_mass/D:CHpos_px/D:CHpos_py/D:CHpos_pz/D:CHpos_energy/D:CHneg_pt/D:CHneg_phi/D:CHneg_eta/D:CHneg_mass/D:CHneg_px/D:CHneg_py/D:CHneg_pz/D:CHneg_energy/D:center_pt/D:center_phi/D:center_eta/D:center_mass/D:center_px/D:center_py/D:center_pz/D:center_energy/D:photon_pt/D:photon_phi/D:photon_eta/D:photon_mass/D:photon_px/D:photon_py/D:photon_pz/D:photon_energy/D:pt/D:phi/D:eta/D:mass/D:px/D:py/D:pz/D:energy/D:grommedMass/D:genDR/D:chargedIso/D:neutralIso/D:egammaIso/D:photon_nPhotons/I:photon_nElectrons/I:genEtaIndex/I");

  // also want a Fill function, that can fill the struct values from the appropriate objects
  // again, so that all editing only needs to be done here in this file
  void FillRecoTwoProngInfo(recoTwoProngInfo_t &recotwopronginfo, TLorentzVector &CHpos, TLorentzVector &CHneg, TLorentzVector &center,
                            TLorentzVector &photon, TLorentzVector &Eta, double grommedmass, bool passed, bool matched, double genDR, int genIndex, 
                            double relchargediso, double relneutraliso, double relegammaiso, int numgamma, int nume) {

    recotwopronginfo.CHpos_pt = CHpos.Pt();
    recotwopronginfo.CHpos_phi = CHpos.Phi();
    recotwopronginfo.CHpos_eta = CHpos.Eta();
    recotwopronginfo.CHpos_mass = CHpos.M();
    recotwopronginfo.CHpos_px = CHpos.Px();
    recotwopronginfo.CHpos_py = CHpos.Py();
    recotwopronginfo.CHpos_pz = CHpos.Pz();
    recotwopronginfo.CHpos_energy = CHpos.E();

    recotwopronginfo.CHneg_pt = CHneg.Pt();
    recotwopronginfo.CHneg_phi = CHneg.Phi();
    recotwopronginfo.CHneg_eta = CHneg.Eta();
    recotwopronginfo.CHneg_mass = CHneg.M();
    recotwopronginfo.CHneg_px = CHneg.Px();
    recotwopronginfo.CHneg_py = CHneg.Py();
    recotwopronginfo.CHneg_pz = CHneg.Pz();
    recotwopronginfo.CHneg_energy = CHneg.E();

    recotwopronginfo.center_pt = center.Pt();
    recotwopronginfo.center_phi = center.Phi();
    recotwopronginfo.center_eta = center.Eta();
    recotwopronginfo.center_mass = center.M();
    recotwopronginfo.center_px = center.Px();
    recotwopronginfo.center_py = center.Py();
    recotwopronginfo.center_pz = center.Pz();
    recotwopronginfo.center_energy = center.E();

    recotwopronginfo.photon_pt = photon.Pt();
    recotwopronginfo.photon_phi = photon.Phi();
    recotwopronginfo.photon_eta = photon.Eta();
    recotwopronginfo.photon_mass = photon.M();
    recotwopronginfo.photon_px = photon.Px();
    recotwopronginfo.photon_py = photon.Py();
    recotwopronginfo.photon_pz = photon.Pz();
    recotwopronginfo.photon_energy = photon.E();

    recotwopronginfo.pt = Eta.Pt();
    recotwopronginfo.phi = Eta.Phi();
    recotwopronginfo.eta = Eta.Eta();
    recotwopronginfo.mass = Eta.M();
    recotwopronginfo.px = Eta.Px();
    recotwopronginfo.py = Eta.Py();
    recotwopronginfo.pz = Eta.Pz();
    recotwopronginfo.energy = Eta.E();
    recotwopronginfo.grommedMass = grommedmass;

    recotwopronginfo.genDR = genDR;

    recotwopronginfo.chargedIso = relchargediso;
    recotwopronginfo.neutralIso = relneutraliso;
    recotwopronginfo.egammaIso = relegammaiso;

    recotwopronginfo.photon_nPhotons = numgamma;
    recotwopronginfo.photon_nElectrons = nume;
    recotwopronginfo.genEtaIndex = genIndex;
  }

  void InitRecoTwoProngInfo(recoTwoProngInfo_t &recotwopronginfo) {
    recotwopronginfo.CHpos_pt = -99.9;
    recotwopronginfo.CHpos_phi = -99.9;
    recotwopronginfo.CHpos_eta = -99.9;
    recotwopronginfo.CHpos_mass = -99.9;
    recotwopronginfo.CHpos_px = -99.9;
    recotwopronginfo.CHpos_py = -99.9;
    recotwopronginfo.CHpos_pz = -99.9;
    recotwopronginfo.CHpos_energy = -99.9;

    recotwopronginfo.CHneg_pt = -99.9;
    recotwopronginfo.CHneg_phi = -99.9;
    recotwopronginfo.CHneg_eta = -99.9;
    recotwopronginfo.CHneg_mass = -99.9;
    recotwopronginfo.CHneg_px = -99.9;
    recotwopronginfo.CHneg_py = -99.9;
    recotwopronginfo.CHneg_pz = -99.9;
    recotwopronginfo.CHneg_energy = -99.9;

    recotwopronginfo.center_pt = -99.9;
    recotwopronginfo.center_phi = -99.9;
    recotwopronginfo.center_eta = -99.9;
    recotwopronginfo.center_mass = -99.9;
    recotwopronginfo.center_px = -99.9;
    recotwopronginfo.center_py = -99.9;
    recotwopronginfo.center_pz = -99.9;
    recotwopronginfo.center_energy = -99.9;

    recotwopronginfo.photon_pt = -99.9;
    recotwopronginfo.photon_phi = -99.9;
    recotwopronginfo.photon_eta = -99.9;
    recotwopronginfo.photon_mass = -99.9;
    recotwopronginfo.photon_px = -99.9;
    recotwopronginfo.photon_py = -99.9;
    recotwopronginfo.photon_pz = -99.9;
    recotwopronginfo.photon_energy = -99.9;

    recotwopronginfo.pt = -99.9;
    recotwopronginfo.phi = -99.9;
    recotwopronginfo.eta = -99.9;
    recotwopronginfo.mass = -99.9;
    recotwopronginfo.px = -99.9;
    recotwopronginfo.py = -99.9;
    recotwopronginfo.pz = -99.9;
    recotwopronginfo.energy = -99.9;
    recotwopronginfo.grommedMass = -99.9;

    recotwopronginfo.genDR = -99.9;

    recotwopronginfo.chargedIso = -99.9;
    recotwopronginfo.neutralIso = -99.9;
    recotwopronginfo.egammaIso = -99.9;

    recotwopronginfo.photon_nPhotons = -99;
    recotwopronginfo.photon_nElectrons = -99;
    recotwopronginfo.genEtaIndex = -99;
  }

} //end of namespace

#endif
