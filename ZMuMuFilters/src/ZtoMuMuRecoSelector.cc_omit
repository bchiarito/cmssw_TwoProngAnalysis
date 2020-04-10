// -*- C++ -*-
//
// Package:    TauHadFilters/ZtoMuMuRecoSelector
// Class:      ZtoMuMuRecoSelector
// 
/**\class ZtoMuMuRecoSelector ZtoMuMuRecoSelector.cc TauHadFilters/ZtoMuMuRecoSelector/plugins/ZtoMuMuRecoSelector.cc

 Description: [one line class summary]

*/
//
// Original Author:  Brandon Chiarito
//         Created:  Wed, 05 Jul 2017 22:31:58 GMT
//
//

// cpp includes
#include <memory>
#include <math.h>
#include <algorithm>
#include <iostream>
// ROOT includes
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TTree.h"
// cmssw framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
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
#include "TwoProngAnalysis/ZTauTauFilters/interface/ZtoTauHadPreSelection.h"

// namesspace defaults
using std::vector;
using std::string;
using std::cout;
using std::endl;

//
// class declaration
//
class ZtoMuMuRecoSelector : public edm::EDFilter {
   public:
      explicit ZtoMuMuRecoSelector(const edm::ParameterSet&);
      ~ZtoMuMuRecoSelector();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;

      virtual void beginJob() override;
      virtual void endJob() override;
      virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // Configuration Parameters
      bool cfg_dumpCutflow;
      bool cfg_reducedSelection; // does not include MT or Pzeta cuts, or trigger requirement

      // EDM Collection lables
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      edm::EDGetTokenT<pat::TauCollection> tauToken_;
      edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
      edm::EDGetTokenT<pat::JetCollection> jetToken_;
      edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
      edm::EDGetTokenT<pat::METCollection> metToken_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfcandsToken_;
      edm::EDGetTokenT<reco::BeamSpot> beamToken_;
      edm::EDGetTokenT<double> rhoToken_;

      // cutflow counts
      int cutflow_total;
      int cutflow_foundTrigger;
      int cutflow_passTrigger;
      int cutflow_passMuon;
      int cutflow_passDiMuon;
      int cutflow_passExtraLeptonVeto;
      int cutflow_passReducedSelection;
      int cutflow_passPreSelection;
      // n-1
      int cutflow_N1_ExtraMuonVeto;
      int cutflow_N1_ExtraElectronVeto;
      int cutflow_N1_Trigger; 
      int cutflow_N1_DR;
      int cutflow_N1_OS;
      int cutflow_N1_MassWindow; 
      // trig eff
      int cutflow_passMuonAndTrigger;
};

//
// constructors
//
ZtoMuMuRecoSelector::ZtoMuMuRecoSelector(const edm::ParameterSet& iConfig) :
  cfg_dumpCutflow(iConfig.getUntrackedParameter<bool>("dumpCutflow")),
  cfg_reducedSelection(iConfig.getUntrackedParameter<bool>("reducedSelection"))
{
  triggerBits_ = consumes<edm::TriggerResults>( edm::InputTag("TriggerResults","","HLT") );
  triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>( edm::InputTag("selectedPatTrigger") );
  triggerPrescales_ = consumes<pat::PackedTriggerPrescales>( edm::InputTag("patTrigger") );
  vtxToken_ = consumes<vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));
  muonToken_ = consumes<std::vector<pat::Muon>>(edm::InputTag("slimmedMuons"));
  electronToken_ = consumes<std::vector<pat::Electron>>(edm::InputTag("slimmedElectrons"));
  tauToken_ = consumes<std::vector<pat::Tau>>(edm::InputTag("slimmedTaus"));
  photonToken_ = consumes<std::vector<pat::Photon>>(edm::InputTag("slimmedPhotons"));
  jetToken_ = consumes<std::vector<pat::Jet>>(edm::InputTag("slimmedJets"));
  fatjetToken_ = consumes<std::vector<pat::Jet>>(edm::InputTag("slimmedJetsAK8"));
  metToken_ = consumes<std::vector<pat::MET>>(edm::InputTag("slimmedMETs"));
  pfcandsToken_ = consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates"));
  beamToken_ = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  rhoToken_ = consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralChargedPileUp"));

  // cutflow counts
  cutflow_total = 0;
  cutflow_foundTrigger = 0;
  cutflow_passTrigger = 0;
  cutflow_passMuon = 0;
  cutflow_passDiMuon = 0;
  cutflow_passExtraLeptonVeto = 0;
  cutflow_passReducedSelection = 0;
  cutflow_passPreSelection = 0;
  // n-1
  cutflow_N1_ExtraMuonVeto = 0;
  cutflow_N1_ExtraElectronVeto = 0;
  cutflow_N1_Trigger = 0;
  cutflow_N1_DR = 0;
  cutflow_N1_OS = 0;
  cutflow_N1_MassWindow = 0;
  // trig eff
  cutflow_passMuonAndTrigger = 0;
}

//
// destructor
//
ZtoMuMuRecoSelector::~ZtoMuMuRecoSelector()
{
}

// ------------ method called on each new Event  ------------
bool
ZtoMuMuRecoSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // get event content
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);

  edm::Handle<pat::TauCollection> taus;
  iEvent.getByToken(tauToken_, taus);

  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);
 
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);

  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);
  pat::MET MET = (*mets)[0];

  edm::Handle<double> rho;
  iEvent.getByToken(rhoToken_, rho);

  // get preselection result
  TauHadFilters::DiMuonPreSelectionResult result;
  result = TauHadFilters::computeDiMuonPreSelectionResult(iEvent, triggerBits, triggerObjects, triggerPrescales, vertices, taus, muons, electrons, jets, mets, rho, TauHadFilters::mediumID, TauHadFilters::vvtightISO);
  bool passDiMuon = (result.tagMuon != NULL && result.tagMuon2 != NULL);

  // count
  cutflow_total += 1;
  if (result.foundTrigger != "") cutflow_foundTrigger += 1;
  if (result.passTrigger || result.passTriggerTk) cutflow_passTrigger += 1;
  if (result.nTagMuons > 0) cutflow_passMuon += 1;
  if (result.passDiMuon) cutflow_passDiMuon += 1;
  if (result.passExtraElectronVeto && result.passExtraMuonVeto) cutflow_passExtraLeptonVeto += 1;

  if (passDiMuon) cutflow_passReducedSelection += 1;
  if (result.passPreSelection) cutflow_passPreSelection += 1;

  if (passDiMuon && ( result.passTrigger || result.passTriggerTk ) && result.passExtraElectronVeto && true)                     cutflow_N1_ExtraMuonVeto += 1;
  if (passDiMuon && ( result.passTrigger || result.passTriggerTk ) && true                         && result.passExtraMuonVeto) cutflow_N1_ExtraElectronVeto += 1;
  if (passDiMuon && true               && result.passExtraElectronVeto && result.passExtraMuonVeto) cutflow_N1_Trigger += 1;

  if (result.passDiMuonDRN1) cutflow_N1_DR += 1;
  if (result.passDiMuonOSN1) cutflow_N1_OS += 1;
  if (result.passDiMuonMassWindowN1) cutflow_N1_MassWindow += 1;

  if ((result.passTrigger || result.passTriggerTk) && result.nTagMuons > 0) cutflow_passMuonAndTrigger += 1;

  // return
  if (!cfg_reducedSelection) return result.passPreSelection;
  else return passDiMuon;
}

void
ZtoMuMuRecoSelector::beginJob()
{
}

void
ZtoMuMuRecoSelector::endJob()
{
  if (cfg_reducedSelection) cout << "\nran filter in reduced mode" << endl;
  else cout << "\nran filter in full preselection mode" << endl;
  if (cfg_dumpCutflow) {
    // cutflow print
    cout << "total " << cutflow_total << endl;
    cout << "cutflow_foundTrigger " << cutflow_foundTrigger << endl;
    cout << "cutflow_passTrigger " << cutflow_passTrigger << endl;
    cout << "cutflow_passMuon " << cutflow_passMuon << endl;
    cout << "cutflow_passDiMuon " << cutflow_passDiMuon << endl;
    cout << "cutflow_passExtraLeptonVeto " << cutflow_passExtraLeptonVeto << endl;
    cout << "cutflow_passReducedSelection " << cutflow_passReducedSelection << endl;
    cout << "cutflow_passPreSelection " << cutflow_passPreSelection << endl;
    cout << "cutflow_N1_ExtraMuonVeto " << cutflow_N1_ExtraMuonVeto << endl;
    cout << "cutflow_N1_ExtraElectronVeto " << cutflow_N1_ExtraElectronVeto << endl;
    cout << "cutflow_N1_Trigger " << cutflow_N1_Trigger << endl; 
    cout << "cutflow_N1_DR " << cutflow_N1_DR << endl; 
    cout << "cutflow_N1_OS " << cutflow_N1_OS << endl; 
    cout << "cutflow_N1_MassWindow " << cutflow_N1_MassWindow << endl; 
    cout << "cutflow_passMuonAndTrigger " << cutflow_passMuonAndTrigger << endl;
  }
}

void
ZtoMuMuRecoSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

void
ZtoMuMuRecoSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}

void
ZtoMuMuRecoSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

void
ZtoMuMuRecoSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

void
ZtoMuMuRecoSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.addUntracked<bool>("dumpCutflow", true);
  desc.addUntracked<bool>("reducedSelection", false);
  descriptions.add("ZtoMuMuRecoSelectorFilter", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZtoMuMuRecoSelector);
