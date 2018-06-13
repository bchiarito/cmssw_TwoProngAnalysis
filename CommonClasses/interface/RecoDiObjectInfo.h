#ifndef RECO_DIOBJECT_INFO_INC
#define RECO_DIOBJECT_INFO_INC

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

namespace ExoDiPhotons
{
  struct recoDiObjectInfo_t {
    Double_t pt;
    Double_t phi;
    Double_t eta;
    Double_t mass;
    Double_t energy;
    Double_t part1_pt;
    Double_t part1_eta;
    Double_t part1_phi;
    Double_t part1_mass;
    Double_t part1_energy;
    Double_t part2_pt;
    Double_t part2_eta;
    Double_t part2_phi;
    Double_t part2_mass;
    Double_t part2_energy;
    Double_t dR;
  };

  // to set branch string
  std::string recoDiObjectBranchDefString("pt/D:phi/D:eta/D:mass/D:energy/D:part1_pt/D:part1_eta/D:part1_phi/D:part1_mass/D:part1_energy/D:part2_pt/D:part2_eta/D:part2_phi/D:part2_mass/D:part2_energy/D:dR/D");

  // fills the struct from two TLorentzVectors
  void FillRecoDiObjectInfo(recoDiObjectInfo_t &recodiobjectinfo, TLorentzVector vec1, TLorentzVector vec2)
  {
    TLorentzVector comb = vec1 + vec2;
    recodiobjectinfo.pt = comb.Pt();
    recodiobjectinfo.phi = comb.Phi();
    recodiobjectinfo.eta = comb.Eta();
    recodiobjectinfo.mass = comb.M();
    recodiobjectinfo.energy = comb.E();

    recodiobjectinfo.part1_pt = vec1.Pt();
    recodiobjectinfo.part1_eta = vec1.Eta();
    recodiobjectinfo.part1_phi = vec1.Phi();
    recodiobjectinfo.part1_mass = vec1.M();
    recodiobjectinfo.part1_energy = vec1.E();

    recodiobjectinfo.part2_pt = vec1.Pt();
    recodiobjectinfo.part2_eta = vec1.Eta();
    recodiobjectinfo.part2_phi = vec1.Phi();
    recodiobjectinfo.part2_mass = vec1.M();
    recodiobjectinfo.part2_energy = vec1.E();

    recodiobjectinfo.dR = vec1.DeltaR(vec2);
  }

  // initializes the struct
  void InitRecoDiObjectInfo(recoDiObjectInfo_t &recodiobjectinfo)
  {
    recodiobjectinfo.pt = -99.9;
    recodiobjectinfo.phi = -99.9;
    recodiobjectinfo.eta = -99.9;
    recodiobjectinfo.mass = -99.9;
    recodiobjectinfo.energy = -99.9;

    recodiobjectinfo.part1_pt = -99.9;
    recodiobjectinfo.part1_eta = -99.9;
    recodiobjectinfo.part1_phi = -99.9;
    recodiobjectinfo.part1_mass = -99.9;
    recodiobjectinfo.part1_energy = -99.9;

    recodiobjectinfo.part2_pt = -99.9;
    recodiobjectinfo.part2_eta = -99.9;
    recodiobjectinfo.part2_phi = -99.9;
    recodiobjectinfo.part2_mass = -99.9;
    recodiobjectinfo.part2_energy = -99.9;

    recodiobjectinfo.dR = -99.9;
  }

} //end of namespace

#endif
