#ifndef ZTOTAUHADPRESLECTION_H
#define ZTOTAUHADPRESLECTION_H

// cpp includes
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
// ROOT includes
#include "TVector2.h"
#include "TLorentzVector.h"
// cmssw framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
// trigger inlcudes
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
// pat includes
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
// edm utilities includes
#include "DataFormats/Math/interface/deltaR.h"

using std::string;
using std::vector;
using std::max;
using std::string;

namespace TauHadFilters
{
  const double Z_MASS = 91.1876;

  const string MUON_TRIGGER = "HLT_IsoMu24";
  const string MUON_TRIGGER_Tk = "HLT_IsoTkMu24";
  const double TRIGGEROBJ_MATCH_DR = 0.1;

  const double MUON_MIN_PT = 26;
  const double MUON_MAX_ETA = 2.1;
  const double MUON_MAX_RELISO = 0.15;
  const double MUON_MAX_DZ = 0.2;
  const double MUON_MAX_DXY = 0.045;

  const double DIMUON_MIN_DR = 0.05;
  const double DIMUON_Z_MASS_MIN = 60.0;
  const double DIMUON_Z_MASS_MAX = 120.0;

  const double LOOSEMUON_MAX_ETA = 2.4;

  const double EXTRAMUON_MIN_PT = 10.0;
  const double EXTRAMUON_MAX_ETA = 2.4;
  const double EXTRAMUON_MAX_RELISO = 0.25;
  const double EXTRAMUON_MAX_DZ = 0.2;
  const double EXTRAMUON_MAX_DXY = 0.045;

  const double EXTRAELECTRON_MIN_PT = 10.0;
  const double EXTRAELECTRON_MAX_ETA = 2.5;
  const double EXTRAELECTRON_MAX_RELISO = 0.3;
  const double EXTRAELECTRON_MAX_DZ = 0.2;
  const double EXTRAELECTRON_MAX_DXY = 0.045;
  const int    EXTRAELECTRON_MAX_LOSTTRACKS = 1;

  const double EXTRADIMUON_MAX_DR = 0.15;
  const double EXTRADIMUON_MIN_PT = 15;
  const double EXTRADIMUON_MAX_ETA = 2.4;
  const double EXTRADIMUON_MAX_RELISO = 0.25;
  const double EXTRADIMUON_MAX_DZ = 0.2;
  const double EXTRADIMUON_MAX_DXY = 0.045;

  const double TAUJET_MIN_LEADINGTRACKPT = 3;
  const double TAUJET_NEARBYGLOBALMUON_DR = 0.4;
  const double TAUJET_NEARBYGLOBALMUON_MIN_PT = 5;
  const double TAUJET_MIN_PT = 20;
  const double TAUJET_MAX_ETA = 2.4;

  const string BTAG_DISCRIMINATOR = "pfCombinedInclusiveSecondaryVertexV2BJetTags";
  const double JET_BTAG_WP = 0.890;
  const double MUONTAUPAIR_MIN_DR = 0.5;
  const double MAX_MT = 40;
  const double MIN_PZETA = -25;

  struct TauHadPreSelectionResult
  {
    // configuration
    bool usePatTau;
    // trigger
    bool passTrigger;
    bool passTriggerTk;
    string foundTrigger;
    string foundTriggerTk;
    // object counts
    int nTagMuons;
    int nProbeTaus;
    // event wide vetos
    bool passExtraMuonVeto;
    bool passDiMuonVeto;
    bool passExtraElectronVeto;
    bool passBtagVeto;
    double highestBtagDiscriminant;
    // muon-tau pair
    bool passMuonTauPair;
    const pat::Muon * tagMuon;
    const pat::Jet * probeTauJet;
    const pat::Tau * probeTau;
    bool pairAndPassMT;
    bool pairAndPassPzeta;
    double MT;
    double Pzeta;
    // full selection
    bool passPreSelection;
  };

  struct DiMuonPreSelectionResult
  {
    // configuration
    bool usePatTau;
    // trigger
    bool passTrigger;
    bool passTriggerTk;
    string foundTrigger;
    string foundTriggerTk;
    // object counts
    int nTagMuons;
    // event wide vetos
    bool passExtraMuonVeto;
    bool passDiMuonVeto;
    bool passExtraElectronVeto;
    // muon-muon pair
    const pat::Muon * tagMuon;
    const pat::Muon * tagMuon2;
    bool passDiMuon;
    bool passDiMuonOSN1;
    bool passDiMuonDRN1;
    bool passDiMuonMassWindowN1;
    // full selection
    bool passPreSelection;
  };

  typedef struct TauHadPreSelectionResult TauHadPreSelectionResult;
  typedef struct DiMuonPreSelectionResult DiMuonPreSelectionResult;

  // forward declare the preselection functions
  struct TauHadPreSelectionResult computePreSelectionResult(const edm::Event&, edm::Handle<edm::TriggerResults>&, edm::Handle<pat::TriggerObjectStandAloneCollection>, edm::Handle<pat::PackedTriggerPrescales>, edm::Handle<reco::VertexCollection>, edm::Handle<pat::TauCollection>, edm::Handle<pat::MuonCollection>, edm::Handle<pat::ElectronCollection>, edm::Handle<pat::JetCollection>, edm::Handle<pat::METCollection>, edm::Handle<double> , bool);

  struct DiMuonPreSelectionResult computeDiMuonPreSelectionResult(const edm::Event&, edm::Handle<edm::TriggerResults>&, edm::Handle<pat::TriggerObjectStandAloneCollection>, edm::Handle<pat::PackedTriggerPrescales>, edm::Handle<reco::VertexCollection>, edm::Handle<pat::TauCollection>, edm::Handle<pat::MuonCollection>, edm::Handle<pat::ElectronCollection>, edm::Handle<pat::JetCollection>, edm::Handle<pat::METCollection>, edm::Handle<double>);

  // helper functions
  double computeMuonIsolation(const pat::Muon * mu)
  {
    return (mu->pfIsolationR04().sumChargedHadronPt + max(0., mu->pfIsolationR04().sumNeutralHadronEt + mu->pfIsolationR04().sumPhotonEt - 0.5*mu->pfIsolationR04().sumPUPt))/mu->pt();
  }

  double computeElectronIsolation(const pat::Electron * el)
  {
    return (el->pfIsolationVariables().sumChargedHadronPt + max(0., el->pfIsolationVariables().sumNeutralHadronEt + el->pfIsolationVariables().sumPhotonEt - 0.5*el->pfIsolationVariables().sumPUPt))/el->pt();
  }
      
  double computeMT(const pat::Muon * mu, const pat::MET * met)
  {
    double dPhi = mu->phi() - met->phi();
    return sqrt(2 * mu->pt() * met->pt() * (1 - cos(dPhi)));    
  }

  double computePzeta(const pat::Muon * mu, const reco::Candidate * tau, const pat::MET * met)
  {
    TVector2 pTmuon;
    TVector2 pTtau;
    TVector2 pTmet;
    TVector2 zeta;
    pTmuon.SetMagPhi(mu->pt(), mu->phi());  
    pTtau.SetMagPhi(tau->pt(), tau->phi());  
    pTmet.SetMagPhi(met->pt(), met->phi());  
    zeta.SetMagPhi(1.0, (mu->phi()+tau->phi())/2.0);
    double PzetaAll = zeta * (pTmuon + pTtau + pTmet);
    double PzetaVis = zeta * (pTmuon + pTtau);
    return PzetaAll - 0.85 * PzetaVis;
  }

  bool betterPair(const pat::Muon * mu, const pat::Jet * tau, const pat::Muon * mu_old, const pat::Jet * tau_old)
  {
    if (mu->pt() + tau->pt() > mu_old->pt() + tau_old->pt()) return true;
    else return false; 
  }

  bool betterPair(const pat::Muon * mu, const pat::Tau * tau, const pat::Muon * mu_old, const pat::Tau * tau_old)
  {
    if (mu->pt() + tau->pt() > mu_old->pt() + tau_old->pt()) return true;
    else return false; 
  }

  // the tau_had tau_mu preselection
  struct TauHadPreSelectionResult computePreSelectionResult(const edm::Event& iEvent, edm::Handle<edm::TriggerResults>& triggerBits, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, edm::Handle<pat::PackedTriggerPrescales> triggerPrescales, edm::Handle<reco::VertexCollection> vertices, edm::Handle<pat::TauCollection> taus, edm::Handle<pat::MuonCollection> muons, edm::Handle<pat::ElectronCollection> electrons, edm::Handle<pat::JetCollection> jets, edm::Handle<pat::METCollection> mets, edm::Handle<double> rho, bool usePatTau = false)
  {
    struct TauHadPreSelectionResult result;

    const reco::Vertex &PV = vertices->front();
    pat::MET MET = (*mets)[0];

    // trigger 
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    string trigger_name = MUON_TRIGGER;
    bool trigger_bit = false;
    string trigger_found = "";
    string trigger_name_tk = MUON_TRIGGER_Tk;
    bool trigger_bit_tk = false;
    string trigger_found_tk = "";
    for (unsigned int i = 0, n = triggerBits->size(); i < n; i++)
    {
       string triggerName = names.triggerName(i);
       std::size_t pos = triggerName.find(trigger_name);
       if ( pos != std::string::npos ) {
         trigger_bit = triggerBits->accept(i);
         trigger_found = triggerName;
       }
       pos = triggerName.find(trigger_name_tk);
       if ( pos != std::string::npos ) {
         trigger_bit_tk = triggerBits->accept(i);
         trigger_found_tk = triggerName;
       }
    }

    // trigger objects
    vector<TLorentzVector> trigger_objects;
    vector<TLorentzVector> trigger_objects_muon;
    vector<TLorentzVector> trigger_objects_muontk;
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
       obj.unpackPathNames(names);
       TLorentzVector triggerobj;
       triggerobj.SetPtEtaPhiE(obj.pt(), obj.eta(), obj.phi(), obj.energy());
       std::vector<std::string> pathNamesAll = obj.pathNames(false);
       std::vector<std::string> pathNamesLast = obj.pathNames(true);
       for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
           bool isBoth = obj.hasPathName( pathNamesAll[h], true, true );
           if (!isBoth) continue;
           string pathName = pathNamesAll[h];
           std::size_t pos = pathName.find(trigger_name);
           if ( pos != std::string::npos ) {
             trigger_objects.push_back(triggerobj);
             trigger_objects_muon.push_back(triggerobj);
           }
           pos = pathName.find(trigger_name_tk);
           if ( pos != std::string::npos ) {
             trigger_objects.push_back(triggerobj);
             trigger_objects_muontk.push_back(triggerobj);
           }
       }
    }

    // muons
    vector<const pat::Muon *> passedMuons;
    for (const pat::Muon &muon : *muons) {
      if (muon.pt() > MUON_MIN_PT &&
          fabs(muon.eta()) < MUON_MAX_ETA &&
          computeMuonIsolation(&muon) < MUON_MAX_RELISO &&
          fabs(muon.muonBestTrack()->dz(PV.position())) < MUON_MAX_DZ &&
          fabs(muon.muonBestTrack()->dxy(PV.position())) < MUON_MAX_DXY &&
          muon.isTightMuon(PV) ) {
        TLorentzVector mu; mu.SetPtEtaPhiM(muon.pt(), muon.eta(), muon.phi(), muon.mass());
        bool matched_to_trigger_obj = false;
        for (TLorentzVector trigobj : trigger_objects) {
          if (trigobj.DeltaR(mu) < TRIGGEROBJ_MATCH_DR) matched_to_trigger_obj = true; }
        if (matched_to_trigger_obj) passedMuons.push_back(&muon);
        }
    }

    // extra lepton veto
    bool skipped_one_passing_muon = false;
    bool extraMuon = false;
    for (const pat::Muon &muon : *muons) {
      if (muon.pt() > EXTRAMUON_MIN_PT && 
          fabs(muon.eta()) < EXTRAMUON_MAX_ETA &&
          computeMuonIsolation(&muon) < EXTRAMUON_MAX_RELISO &&
          fabs(muon.muonBestTrack()->dz(PV.position())) < EXTRAMUON_MAX_DZ &&
          fabs(muon.muonBestTrack()->dxy(PV.position())) < EXTRAMUON_MAX_DXY &&
          muon.isMediumMuon() ) {
        if (muon.pt() > MUON_MIN_PT &&
            fabs(muon.eta()) < MUON_MAX_ETA &&
            computeMuonIsolation(&muon) < MUON_MAX_RELISO &&
            fabs(muon.muonBestTrack()->dz(PV.position())) < MUON_MAX_DZ &&
            fabs(muon.muonBestTrack()->dxy(PV.position())) < MUON_MAX_DXY &&
            muon.isMediumMuon() ) {
          if (skipped_one_passing_muon) extraMuon = true;
          if (!skipped_one_passing_muon) skipped_one_passing_muon = true;
        } else {
          extraMuon = true;
        }
      }
    }
    bool extraElectron = false;
    for (const pat::Electron &electron : *electrons) {
      if (electron.pt() > EXTRAELECTRON_MIN_PT &&
          fabs(electron.eta()) < EXTRAELECTRON_MAX_ETA &&
          fabs(electron.gsfTrack()->dz(PV.position())) < EXTRAELECTRON_MAX_DZ &&
          fabs(electron.gsfTrack()->dxy(PV.position())) < EXTRAELECTRON_MAX_DXY &&
          electron.passConversionVeto() &&
          computeElectronIsolation(&electron) < EXTRAELECTRON_MAX_RELISO &&
          electron.gsfTrack()->lost() <= EXTRAELECTRON_MAX_LOSTTRACKS)
        extraElectron = true;
    }
    bool diMuon = false;
    for (const pat::Muon &muon1 : *muons) {
      for (const pat::Muon &muon2 : *muons) {
        double DR = deltaR(muon1.eta(), muon1.phi(), muon2.eta(), muon2.phi());
        if (DR <= EXTRADIMUON_MAX_DR) continue;
        if (muon1.pt() > EXTRADIMUON_MIN_PT &&
            fabs(muon1.eta()) < EXTRADIMUON_MAX_ETA &&
            computeMuonIsolation(&muon1) < EXTRADIMUON_MAX_RELISO &&
            fabs(muon1.muonBestTrack()->dz(PV.position())) < EXTRADIMUON_MAX_DZ &&
            fabs(muon1.muonBestTrack()->dxy(PV.position())) < EXTRADIMUON_MAX_DXY &&
            muon1.isPFMuon() &&
            muon1.isGlobalMuon() &&
            muon1.isTrackerMuon() &&
            muon2.pt() > EXTRADIMUON_MIN_PT &&
            fabs(muon2.eta()) < EXTRADIMUON_MAX_ETA &&
            computeMuonIsolation(&muon2) < EXTRADIMUON_MAX_RELISO &&
            fabs(muon2.muonBestTrack()->dz(PV.position())) < EXTRADIMUON_MAX_DZ &&
            fabs(muon2.muonBestTrack()->dxy(PV.position())) < EXTRADIMUON_MAX_DXY &&
            muon2.isPFMuon() &&
            muon2.isGlobalMuon() &&
            muon2.isTrackerMuon() &&
            muon1.charge() * muon2.charge() < 0)
          diMuon = true;
      }
    }  

    // tau from jets and btag veto
    bool atLeastOneBTag = false;
    double highestBTagDiscriminant = -1.0;
    vector<const pat::Jet *> tauJetCands;
    vector<int> tauJetCandsCharge;
    for (const pat::Jet &jet : *jets) {
      if (jet.pt() > TAUJET_MIN_PT && fabs(jet.eta()) < TAUJET_MAX_ETA) {
        double btagdiscriminant = jet.bDiscriminator(BTAG_DISCRIMINATOR);
        if (btagdiscriminant > highestBTagDiscriminant) highestBTagDiscriminant = btagdiscriminant;
      }
      if (jet.pt() > TAUJET_MIN_PT && fabs(jet.eta()) < TAUJET_MAX_ETA && jet.bDiscriminator(BTAG_DISCRIMINATOR) > JET_BTAG_WP)
        atLeastOneBTag = true;
      bool noNearbyGlobalMuon = true;
      for (const pat::Muon &muon : *muons) {
        if (!muon.isGlobalMuon()) continue;
        if (muon.pt() < TAUJET_NEARBYGLOBALMUON_MIN_PT) continue;
        double deltaR = reco::deltaR(muon.eta(),muon.phi(),jet.eta(),jet.phi());
        if (deltaR < TAUJET_NEARBYGLOBALMUON_DR) noNearbyGlobalMuon = false;
      }
      if (jet.pt() > TAUJET_MIN_PT && fabs(jet.eta()) < TAUJET_MAX_ETA && noNearbyGlobalMuon) {
        double leading_track_pt = 0;
        int leading_track_charge = 0;
        for (unsigned int i = 0; i < jet.numberOfDaughters(); i++) {
          if (jet.daughter(i)->pdgId() == 22 || jet.daughter(i)->pdgId() == 111) continue;
          if (jet.daughter(i)->pt() > leading_track_pt) {
            leading_track_pt = jet.daughter(i)->pt();
            leading_track_charge = jet.daughter(i)->charge();
          }
        }
        if (leading_track_pt > TAUJET_MIN_LEADINGTRACKPT)
          tauJetCands.push_back(&jet);
          tauJetCandsCharge.push_back(leading_track_charge);
      }
    }

    // tau from pat::tau objects
    vector<const pat::Tau *> tauCands;
    for (const pat::Tau &tau : *taus) {
      bool noNearbyGlobalMuon = true;
      for (const pat::Muon &muon : *muons) {
        if (!muon.isGlobalMuon()) continue;
        if (muon.pt() < TAUJET_NEARBYGLOBALMUON_MIN_PT) continue;
        double deltaR = reco::deltaR(muon.eta(),muon.phi(),tau.eta(),tau.phi());
        if (deltaR < TAUJET_NEARBYGLOBALMUON_DR) noNearbyGlobalMuon = false;
      }
      if (!(tau.pt() > TAUJET_MIN_PT && fabs(tau.eta()) < TAUJET_MAX_ETA)) continue;
      if (!noNearbyGlobalMuon) continue;
      if (!(tau.leadChargedHadrCand()->pt() > TAUJET_MIN_LEADINGTRACKPT)) continue;
      if (!(tau.tauID("againstMuonTight3") > 0.5 && tau.tauID("againstElectronVLooseMVA6") > 0.5)) continue;
      if (tau.isTauIDAvailable("decayModeFindingOldDMs") && tau.tauID("decayModeFindingOldDMs") < 0.5) continue;
      tauCands.push_back(&tau);
    }

    // choose a muon-tau pair that passes DR and opp sign
    int muon_index = -1;
    int tau_index = -1;
    if (!usePatTau)
    { 
      for (unsigned int i = 0; i < passedMuons.size(); i++) {
        const pat::Muon & muon = *passedMuons[i];
        for (unsigned int j = 0; j < tauJetCands.size(); j++) {
          const pat::Jet & tau = *tauJetCands[j];
          int tau_charge = tauJetCandsCharge[j];
          double DR = reco::deltaR(muon.eta(),muon.phi(),tau.eta(),tau.phi());
          int charge_product = muon.charge() * tau_charge;
          if (DR < MUONTAUPAIR_MIN_DR || charge_product > 0) continue;
          if (muon_index == -1) { muon_index = i; tau_index = j; }
          else if ( betterPair(&muon, &tau, passedMuons[muon_index], tauJetCands[tau_index]) ) { muon_index = i; tau_index = j; }
        }
      }
      if (muon_index != -1) {
        result.tagMuon = passedMuons[muon_index];
        result.probeTauJet = tauJetCands[tau_index];
        result.probeTau = NULL;
        result.Pzeta = computePzeta(result.tagMuon, result.probeTauJet, &MET);
        result.MT = computeMT(result.tagMuon, &MET);
      }
    }
    else
    { 
      for (unsigned int i = 0; i < passedMuons.size(); i++) {
        const pat::Muon & muon = *passedMuons[i];
        for (unsigned int j = 0; j < tauCands.size(); j++) {
          const pat::Tau & tau = *tauCands[j];
          double DR = reco::deltaR(muon.eta(),muon.phi(),tau.eta(),tau.phi());
          int charge_product = muon.charge() * tau.charge();
          if (DR < MUONTAUPAIR_MIN_DR || charge_product > 0) continue;
          if (muon_index == -1) { muon_index = i; tau_index = j; }
          else if ( betterPair(&muon, &tau, passedMuons[muon_index], tauCands[tau_index]) ) { muon_index = i; tau_index = j; }
        }
      }
      if (muon_index != -1) {
        result.tagMuon = passedMuons[muon_index];
        result.probeTauJet = NULL;
        result.probeTau = tauCands[tau_index];
        result.Pzeta = computePzeta(result.tagMuon, result.probeTau, &MET);
        result.MT = computeMT(result.tagMuon, &MET);
      }
    }

    /// determine selection decisions
    result.usePatTau = usePatTau;
    // trigger
    result.passTrigger = trigger_bit;
    result.passTriggerTk = trigger_bit_tk;
    result.foundTrigger = trigger_found;
    result.foundTriggerTk = trigger_found_tk;
    // object counts
    result.nTagMuons = passedMuons.size();
    if (!usePatTau) result.nProbeTaus = tauJetCands.size();
    if (usePatTau)  result.nProbeTaus = tauCands.size();
    // event wide vetos
    result.passExtraMuonVeto = !extraMuon;
    result.passDiMuonVeto = !diMuon;
    result.passExtraElectronVeto = !extraElectron;
    result.passBtagVeto = !atLeastOneBTag;
    result.highestBtagDiscriminant = highestBTagDiscriminant;
    // muon tau pair
    result.passMuonTauPair = muon_index != -1;
    if (result.passMuonTauPair) {
      result.pairAndPassMT = result.MT < MAX_MT;
      result.pairAndPassPzeta = result.Pzeta > MIN_PZETA;
    } else {
      result.pairAndPassMT = false;
      result.pairAndPassPzeta = false;
      result.MT = -999.9;
      result.Pzeta = -999.9;
      result.tagMuon = NULL;
      result.probeTauJet = NULL;
      result.probeTau = NULL;
    }
    // full preselection
    result.passPreSelection = (result.passTrigger || result.passTriggerTk) &&
                              result.passMuonTauPair && result.pairAndPassMT && result.pairAndPassPzeta && 
                              result.passExtraMuonVeto && result.passDiMuonVeto && result.passExtraElectronVeto && 
                              result.passBtagVeto;

    return result;
  }


  // the dimuon preselection
  struct DiMuonPreSelectionResult computeDiMuonPreSelectionResult(const edm::Event& iEvent, edm::Handle<edm::TriggerResults>& triggerBits, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, edm::Handle<pat::PackedTriggerPrescales> triggerPrescales, edm::Handle<reco::VertexCollection> vertices, edm::Handle<pat::TauCollection> taus, edm::Handle<pat::MuonCollection> muons, edm::Handle<pat::ElectronCollection> electrons, edm::Handle<pat::JetCollection> jets, edm::Handle<pat::METCollection> mets, edm::Handle<double> rho)
  {
    struct DiMuonPreSelectionResult result;

    const reco::Vertex &PV = vertices->front();
    pat::MET MET = (*mets)[0];

    // trigger 
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    string trigger_name = MUON_TRIGGER;
    bool trigger_bit = false;
    string trigger_found = "";
    string trigger_name_tk = MUON_TRIGGER_Tk;
    bool trigger_bit_tk = false;
    string trigger_found_tk = "";
    for (unsigned int i = 0, n = triggerBits->size(); i < n; i++)
    {
       string triggerName = names.triggerName(i);
       std::size_t pos = triggerName.find(trigger_name);
       if ( pos != std::string::npos ) {
         trigger_bit = triggerBits->accept(i);
         trigger_found = triggerName;
       }
       pos = triggerName.find(trigger_name_tk);
       if ( pos != std::string::npos ) {
         trigger_bit_tk = triggerBits->accept(i);
         trigger_found_tk = triggerName;
       }
    }

    // trigger objects
    vector<TLorentzVector> trigger_objects;
    vector<TLorentzVector> trigger_objects_muon;
    vector<TLorentzVector> trigger_objects_muontk;
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
       obj.unpackPathNames(names);
       TLorentzVector triggerobj;
       triggerobj.SetPtEtaPhiE(obj.pt(), obj.eta(), obj.phi(), obj.energy());
       std::vector<std::string> pathNamesAll = obj.pathNames(false);
       std::vector<std::string> pathNamesLast = obj.pathNames(true);
       for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
           bool isBoth = obj.hasPathName( pathNamesAll[h], true, true );
           if (!isBoth) continue;
           string pathName = pathNamesAll[h];
           std::size_t pos = pathName.find(trigger_name);
           if ( pos != std::string::npos ) {
             trigger_objects.push_back(triggerobj);
             trigger_objects_muon.push_back(triggerobj);
           }
           pos = pathName.find(trigger_name_tk);
           if ( pos != std::string::npos ) {
             trigger_objects.push_back(triggerobj);
             trigger_objects_muontk.push_back(triggerobj);
           }
       }
    }
    
    // muons
    vector<const pat::Muon *> passedMuons;
    for (const pat::Muon &muon : *muons) {
      if (muon.pt() > MUON_MIN_PT &&
          fabs(muon.eta()) < MUON_MAX_ETA &&
          computeMuonIsolation(&muon) < MUON_MAX_RELISO &&
          fabs(muon.muonBestTrack()->dz(PV.position())) < MUON_MAX_DZ &&
          fabs(muon.muonBestTrack()->dxy(PV.position())) < MUON_MAX_DXY &&
          muon.isTightMuon(PV) ) {
        TLorentzVector mu; mu.SetPtEtaPhiM(muon.pt(), muon.eta(), muon.phi(), muon.mass());
        bool matched_to_trigger_obj = false;
        for (TLorentzVector trigobj : trigger_objects) {
          if (trigobj.DeltaR(mu) < TRIGGEROBJ_MATCH_DR) matched_to_trigger_obj = true; }
        if (matched_to_trigger_obj) passedMuons.push_back(&muon);
        }
    }

    // looser muons
    vector<const pat::Muon *> passedLooseMuons;
    for (const pat::Muon &muon : *muons) {
      if (muon.pt() > MUON_MIN_PT &&
          fabs(muon.eta()) < LOOSEMUON_MAX_ETA &&
          computeMuonIsolation(&muon) < MUON_MAX_RELISO &&
          fabs(muon.muonBestTrack()->dz(PV.position())) < MUON_MAX_DZ &&
          fabs(muon.muonBestTrack()->dxy(PV.position())) < MUON_MAX_DXY &&
          muon.isTightMuon(PV) ) {
        TLorentzVector mu; mu.SetPtEtaPhiM(muon.pt(), muon.eta(), muon.phi(), muon.mass());
        bool matched_to_trigger_obj = false;
        for (TLorentzVector trigobj : trigger_objects) {
          if (trigobj.DeltaR(mu) < TRIGGEROBJ_MATCH_DR) matched_to_trigger_obj = true; }
        if (matched_to_trigger_obj) passedLooseMuons.push_back(&muon);
      }
    }

    // extra lepton veto
    int skipped_muons = 0;
    bool extraMuon = false;
    for (const pat::Muon &muon : *muons) {
      if (!(muon.pt() > EXTRAMUON_MIN_PT && 
          fabs(muon.eta()) < EXTRAMUON_MAX_ETA &&
          computeMuonIsolation(&muon) < EXTRAMUON_MAX_RELISO &&
          fabs(muon.muonBestTrack()->dz(PV.position())) < EXTRAMUON_MAX_DZ &&
          fabs(muon.muonBestTrack()->dxy(PV.position())) < EXTRAMUON_MAX_DXY &&
          muon.isMediumMuon() ) )
        continue;
      if (skipped_muons < 2) {
        skipped_muons += 1;
      } else {
        extraMuon = true;
      }
    }
    bool extraElectron = false;
    for (const pat::Electron &electron : *electrons) {
      if (electron.pt() > EXTRAELECTRON_MIN_PT &&
          fabs(electron.eta()) < EXTRAELECTRON_MAX_ETA &&
          fabs(electron.gsfTrack()->dz(PV.position())) < EXTRAELECTRON_MAX_DZ &&
          fabs(electron.gsfTrack()->dxy(PV.position())) < EXTRAELECTRON_MAX_DXY &&
          electron.passConversionVeto() &&
          computeElectronIsolation(&electron) < EXTRAELECTRON_MAX_RELISO &&
          electron.gsfTrack()->lost() <= EXTRAELECTRON_MAX_LOSTTRACKS)
        extraElectron = true;
    }

    // choose the muon-muon system
    bool passDiMuon = false;
    result.tagMuon = NULL;
    result.tagMuon2 = NULL;
    double best_mmumu = -1;
    bool passDiMuonOS = false;
    bool passDiMuonDR = false;
    bool passDiMuonMassWindow = false;
    bool passDiMuonOSN1 = false;
    bool passDiMuonDRN1 = false;
    bool passDiMuonMassWindowN1 = false;
    if (passedMuons.size() >= 1) { // >=1 leading mu
      for (unsigned int i = 0; i < passedLooseMuons.size(); i++) {
        const pat::Muon & muon1 = *passedLooseMuons[i];
        for (unsigned int j = 0; j < passedLooseMuons.size(); j++) {
          const pat::Muon & muon2 = *passedLooseMuons[j];
          TLorentzVector mu1; mu1.SetPtEtaPhiM(muon1.pt(), muon1.eta(), muon1.phi(), muon1.mass());
          TLorentzVector mu2; mu2.SetPtEtaPhiM(muon2.pt(), muon2.eta(), muon2.phi(), muon2.mass());
          double mmumu = (mu1+mu2).M();
          if (reco::deltaR(muon1.eta(), muon1.phi(), muon2.eta(), muon2.phi()) > DIMUON_MIN_DR) passDiMuonDR = true; // dR > 0.5
          if (muon1.charge() * muon2.charge() < 0) passDiMuonOS = true; // OS
          if (mmumu >= DIMUON_Z_MASS_MIN && mmumu <= DIMUON_Z_MASS_MAX) passDiMuonMassWindow = true; // m_mumu in (60, 120)
          if (passDiMuonDR && passDiMuonOS && !passDiMuonMassWindow) passDiMuonMassWindowN1 = true;
          if (passDiMuonDR && !passDiMuonOS && passDiMuonMassWindow) passDiMuonOSN1 = true;
          if (!passDiMuonDR && passDiMuonOS && passDiMuonMassWindow) passDiMuonDRN1 = true;
          if (passDiMuonDR && passDiMuonOS && passDiMuonMassWindow) { passDiMuonDRN1 = true; passDiMuonOSN1 = true; passDiMuonMassWindowN1 = true; }
          if (!passDiMuonDR || !passDiMuonOS || !passDiMuonMassWindow) continue;
          passDiMuon = true;
          if ( fabs(mmumu - Z_MASS) < fabs(best_mmumu - Z_MASS) ) {
            if (muon1.pt() > muon2.pt()) {
              result.tagMuon = &muon1;
              result.tagMuon2 = &muon2;
            } else {
              result.tagMuon = &muon2;
              result.tagMuon2 = &muon1;
            }
            best_mmumu = mmumu;            
          }   
        }
      } 
    }

    
    /// determine selection decisions
    // trigger
    result.passTrigger = trigger_bit;
    result.passTriggerTk = trigger_bit_tk;
    result.foundTrigger = trigger_found;
    result.foundTriggerTk = trigger_found_tk;
    // object counts
    result.nTagMuons = passedLooseMuons.size();
    // event wide vetos
    result.passExtraMuonVeto = !extraMuon;
    result.passExtraElectronVeto = !extraElectron;
    // muon muon pair
    result.passDiMuon = passedMuons.size() >= 1 && passedLooseMuons.size() >= 2;
    result.passDiMuonOSN1 = passDiMuonOSN1;
    result.passDiMuonDRN1 = passDiMuonDRN1;
    result.passDiMuonMassWindowN1 = passDiMuonMassWindowN1;
    // full preselection
    result.passPreSelection = (result.passTrigger || result.passTriggerTk) &&
                              result.passExtraMuonVeto && result.passExtraElectronVeto && 
                              passDiMuon;

    return result;
  }
}

#endif
