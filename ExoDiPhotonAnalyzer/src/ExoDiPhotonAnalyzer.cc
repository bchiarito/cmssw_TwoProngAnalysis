// -*- C++ -*-
//
// Package:    ExoDiPhotonAnalyzer
// Class:      ExoDiPhotonAnalyzer
// 
/**\class ExoDiPhotonAnalyzer ExoDiPhotonAnalyzer.cc DiPhotonAnalysis/ExoDiPhotonAnalyzer/src/ExoDiPhotonAnalyzer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Conor Henderson,40 1-B01,+41227671674,
//         Created:  Thu May  6 17:26:16 CEST 2010
// $Id: ExoDiPhotonAnalyzer.cc,v 1.32 2013/02/11 15:07:42 charaf Exp $
//
//


// system include files
#include <memory>

#include <algorithm>
#include <vector>
#include <utility>  // for std::pair
#include "TClonesArray.h"
#include "TVector3.h"
#include "TLorentzVector.h"

// for new photon code block
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"


// to use TfileService for histograms and trees
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TString.h"

// geometry
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
//#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
//#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
//#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
//#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"


//for vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// for beamspot
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// for ecal
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"


//for photons
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TMath.h"

//for trigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
// #include "CondFormats/L1TObjects/interface/L1GtTriggerMenuFwd.h"
// #include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h" 
// #include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerRecord.h" 
// #include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
// #include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
// #include "L1Trigger/GlobalTrigger/plugins/L1GlobalTrigger.h"
// #include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

// new CommonClasses approach
// these objects are all in the namespace 'ExoDiPhotons'
#include "DiPhotonAnalysis/CommonClasses/interface/RecoPhotonInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/RecoTwoProngInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/RecoDiObjectInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/GenParticleInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/TriggerInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/PhotonID.h"
#include "DiPhotonAnalysis/CommonClasses/interface/EventAndVertexInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/DiphotonInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/JetInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/ConversionInfo.h"

//new for PU gen
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// new for LumiReweighting 
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

//new for PF ID definition
#include "DiPhotonAnalysis/CommonClasses/interface/PFPhotonID.h"

//for conversion safe electron veto
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//new for PFIsolation code
//#include "EGamma/EGammaAnalysisTools/src/PFIsolationEstimator.cc"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

//-----------------taken from Ilya-----------------
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
//-----------------taken from Ilya-----------------

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

// for jets
#include "DataFormats/PatCandidates/interface/Jet.h"

// other miniaod objects
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

using namespace std;

  double phoKappaHighPtID(const pat::Photon *);
  double phoEAHighPtID(const pat::Photon* );
  double phoAlphaHighPtID(const pat::Photon *);
  bool passCorPhoIsoHighPtID(const pat::Photon* , double );
  bool passSigmaIetaIetaCut(const pat::Photon* , bool );
  bool passChargedHadronCut(const pat::Photon* );
  bool passHadTowerOverEmCut(const pat::Photon*);
  double corPhoIsoHighPtID(const pat::Photon*, double );
  bool compareCandsByPt(const edm::Ptr<const reco::Candidate> , const edm::Ptr<const reco::Candidate>);

//
// class declaration
//


class ExoDiPhotonAnalyzer : public edm::EDAnalyzer {
public:
  explicit ExoDiPhotonAnalyzer(const edm::ParameterSet&);
  ~ExoDiPhotonAnalyzer();


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // my functions
  bool isNeutral(int one, int two, int three);
  bool isCharged(int one, int two, int three);

  // ----------member data ---------------------------
  bool               fisMC;  // option to decide if MC or Data     
  bool               fisSignal;  // option to decide if Signal MC or Backround MC
  bool               fDebug;  // if set to False, mean to limit per event stdout output
  bool               fchargedDecayCutflow;  // option to print cutflow of the charged selection to stdout at the end of the job
  bool               fisAOD;

  float fCutflow_total;
  float fCutflow_oneCand;
  float fCutflow_twoCand;
  float fCutflow_onePass;
  float fCutflow_twoPass;
  float fCutflow_oneMatch;
  float fCutflow_twoMatch;
  float fCutflow_onePassMatch;
  float fCutflow_twoPassMatch;
  float fCutflow_passCharged;
  float fCutflow_passNeutral;
  float fCutflow_passEGamma;
  float fCutflow_passPhotonPt;
 
  // tools for clusters
  std::auto_ptr<noZS::EcalClusterLazyTools> lazyTools_;
  edm::InputTag recHitsEBTag_;
  edm::InputTag recHitsEETag_;
  edm::EDGetTokenT<EcalRecHitCollection> recHitsEBToken;
  edm::EDGetTokenT<EcalRecHitCollection> recHitsEEToken;

  // photon subroutines
  bool photon_isSaturated(const pat::Photon*, const EcalRecHitCollection *, const EcalRecHitCollection *,
       const CaloSubdetectorTopology*, const CaloSubdetectorTopology*);

  bool photon_passHighPtID(const pat::Photon*, double , bool ) ;

  TTree *fTree2;

  double fRho25;

  edm::EDGetTokenT<pat::PackedCandidateCollection> pfcandsToken_;
  edm::EDGetTokenT<vector<reco::GenParticle>> genToken_;
  edm::EDGetTokenT<vector<reco::Vertex>> pvToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamToken_;
  edm::EDGetTokenT<std::vector<pat::Jet>> ak4Token_;
  edm::EDGetTokenT<std::vector<pat::Photon>> photonToken_;
  edm::EDGetTokenT<std::vector<pat::Electron>> electronToken_;
  edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
  edm::EDGetTokenT<std::vector<pat::MET>> metToken_;
  double fCandidatePairDR;
  double fCandidatePairMinPt;
  double fCandidatePairIsolationDR;
  double fCandidatePairPhiBox;
  double fCandidatePairEtaBox;
  double fCandidatePairPhotonPtCut;
  double fCandidatePairChargedIsoCut;
  double fCandidatePairNeutralIsoCut;
  double fCandidatePairEGammaIsoCut;
  double fCandidatePairChargedIsoFakeCut;
  double fCandidatePairNeutralIsoFakeCut;
  double fCandidatePairEGammaIsoFakeCut;
  double fCandidatePairGenMatchDR;
  bool fOmitChargedDecayCode;
  bool fSkipPhotonMCCode;
  bool fTwoProngFakeRateCalcOnly;
  TH2F *fTwoProngFakeNume_even_pt;
  TH2F *fTwoProngFakeDeno_even_pt;
  TH2F *fTwoProngFakeRate_even_pt;
  TH2F *fTwoProngFakeNume_even_eta;
  TH2F *fTwoProngFakeDeno_even_eta;
  TH2F *fTwoProngFakeRate_even_eta;
  TH2F *fTwoProngFakeNume_even_phi;
  TH2F *fTwoProngFakeDeno_even_phi;
  TH2F *fTwoProngFakeRate_even_phi;
  TH2F *fTwoProngFakeNume_odd_pt;
  TH2F *fTwoProngFakeDeno_odd_pt;
  TH2F *fTwoProngFakeRate_odd_pt;
  TH2F *fTwoProngFakeNume_odd_eta;
  TH2F *fTwoProngFakeDeno_odd_eta;
  TH2F *fTwoProngFakeRate_odd_eta;
  TH2F *fTwoProngFakeNume_odd_phi;
  TH2F *fTwoProngFakeDeno_odd_phi;
  TH2F *fTwoProngFakeRate_odd_phi;

  edm::EDGetTokenT<double> rhoToken_;
  
  // MiniAOD case data members
  edm::EDGetToken photonsMiniAODToken_;
  edm::EDGetToken gedphotonsToken_;
  edm::EDGetToken patPhotonToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesMiniAODToken_;
  edm::EDGetToken genEvtInfoProdToken_;
  
  Float_t rho_;      // the rho variable

  // ** charged decay analysis **
  int fEventNum;
  int fRunNum;
  int fLumiNum;
  int fNumPVs;
  int fNumPF;
  int fNumPrunedPF;
  int fNumTwoProng;
  int fNumTwoProngPass;
  int fNumTwoProngMatched;
  int fNumTwoProngPassChargedIso;
  int fNumTwoProngPassNeutralIso;
  int fNumTwoProngPassEGammaIso;
  int fNumTwoProngPassPhotonPtIso;
  int fNumTwoProngFake;
  int fNumTightPhotons;
  int fNumTightPhotons_v2;
  int fNumOffset;
  double fHT;
  double fMET;
  double fMET_phi;
  int fNumElectrons;
  int fNumMuons;

  int fNumAK4jets;
  vector<Double_t> fAK4jet_pt;
  vector<Double_t> fAK4jet_eta;
  vector<Double_t> fAK4jet_phi;
  vector<Double_t> fAK4jet_mass;
  vector<Double_t> fAK4jet_px;
  vector<Double_t> fAK4jet_py;
  vector<Double_t> fAK4jet_pz;
  vector<Double_t> fAK4jet_energy;

  int fNumPhotons; 
  vector<Double_t> fPhoton_pt;
  vector<Double_t> fPhoton_eta;
  vector<Double_t> fPhoton_phi;
  vector<Double_t> fPhoton_mass;
  vector<Double_t> fPhoton_px;
  vector<Double_t> fPhoton_py;
  vector<Double_t> fPhoton_pz;
  vector<Double_t> fPhoton_energy;

  vector<Double_t> fCand_pt;
  vector<Double_t> fCand_eta;
  vector<Double_t> fCand_phi;
  vector<Double_t> fCand_mass;
  vector<Double_t> fCand_Mass;
  vector<Double_t> fCand_CHpos_pt;
  vector<Double_t> fCand_CHpos_eta;
  vector<Double_t> fCand_CHpos_phi;
  vector<Double_t> fCand_CHpos_mass;
  vector<Double_t> fCand_CHneg_pt;
  vector<Double_t> fCand_CHneg_eta;
  vector<Double_t> fCand_CHneg_phi;
  vector<Double_t> fCand_CHneg_mass;
  vector<Double_t> fCand_photon_pt;
  vector<Double_t> fCand_photon_eta;
  vector<Double_t> fCand_photon_phi;
  vector<Double_t> fCand_photon_mass;
  vector<Double_t> fCand_photon_nGamma;
  vector<Double_t> fCand_photon_nElectron;
  vector<Double_t> fCand_chargedIso;
  vector<Double_t> fCand_neutralIso;
  vector<Double_t> fCand_egammaIso;
  vector<Double_t> fCand_CHpos_vz;
  vector<Double_t> fCand_CHpos_vx;
  vector<Double_t> fCand_CHpos_vy;
  vector<Double_t> fCand_CHpos_dz;
  vector<Double_t> fCand_CHpos_dz_PV;
  vector<Double_t> fCand_CHpos_dz_beamspot;
  vector<Double_t> fCand_CHpos_dxy;
  vector<Double_t> fCand_CHpos_dxy_PV;
  vector<Double_t> fCand_CHpos_dxy_beamspot;
  vector<Double_t> fCand_CHneg_vz;
  vector<Double_t> fCand_CHneg_vx;
  vector<Double_t> fCand_CHneg_vy;
  vector<Double_t> fCand_CHneg_dz;
  vector<Double_t> fCand_CHneg_dz_PV;
  vector<Double_t> fCand_CHneg_dz_beamspot;
  vector<Double_t> fCand_CHneg_dxy;
  vector<Double_t> fCand_CHneg_dxy_PV;
  vector<Double_t> fCand_CHneg_dxy_beamspot;
  vector<Double_t> fCand_isoPF_vz;
  vector<Double_t> fCand_isoPF_vx;
  vector<Double_t> fCand_isoPF_vy;
  vector<Double_t> fCand_isoPF_dz;
  vector<Double_t> fCand_isoPF_dz_PV;
  vector<Double_t> fCand_isoPF_dz_beamspot;
  vector<Double_t> fCand_isoPF_dxy;
  vector<Double_t> fCand_isoPF_dxy_PV;
  vector<Double_t> fCand_isoPF_dxy_beamspot;
  vector<Int_t> fCand_nChargedIsoCone;
  vector<Int_t> fCand_nNeutralIsoCone;
  vector<Int_t> fCand_nEGammaIsoCone;
  vector<Double_t> fCand_genDR;
  vector<Bool_t> fCand_tight;
  vector<Bool_t> fCand_passChargedIso;
  vector<Bool_t> fCand_passNeutralIso;
  vector<Bool_t> fCand_passEGammaIso;
  vector<Bool_t> fCand_passPhotonPt;
  vector<Bool_t> fCand_loose;
  vector<Bool_t> fCand_match;

  vector<Double_t> fTwoProng_pt;
  vector<Double_t> fTwoProng_eta;
  vector<Double_t> fTwoProng_phi;
  vector<Double_t> fTwoProng_mass;
  vector<Double_t> fTwoProng_Mass;
  vector<Double_t> fTwoProng_px;
  vector<Double_t> fTwoProng_py;
  vector<Double_t> fTwoProng_pz;
  vector<Double_t> fTwoProng_energy;
  vector<Double_t> fTwoProng_CHpos_pt;
  vector<Double_t> fTwoProng_CHpos_eta;
  vector<Double_t> fTwoProng_CHpos_phi;
  vector<Double_t> fTwoProng_CHpos_mass;
  vector<Double_t> fTwoProng_CHneg_pt;
  vector<Double_t> fTwoProng_CHneg_eta;
  vector<Double_t> fTwoProng_CHneg_phi;
  vector<Double_t> fTwoProng_CHneg_mass;
  vector<Double_t> fTwoProng_CHpos_vz;
  vector<Double_t> fTwoProng_CHpos_vx;
  vector<Double_t> fTwoProng_CHpos_vy;
  vector<Double_t> fTwoProng_CHpos_dz;
  vector<Double_t> fTwoProng_CHpos_dz_PV;
  vector<Double_t> fTwoProng_CHpos_dz_beamspot;
  vector<Double_t> fTwoProng_CHpos_dxy;
  vector<Double_t> fTwoProng_CHpos_dxy_PV;
  vector<Double_t> fTwoProng_CHpos_dxy_beamspot;
  vector<Double_t> fTwoProng_CHneg_vz;
  vector<Double_t> fTwoProng_CHneg_vx;
  vector<Double_t> fTwoProng_CHneg_vy;
  vector<Double_t> fTwoProng_CHneg_dz;
  vector<Double_t> fTwoProng_CHneg_dz_PV;
  vector<Double_t> fTwoProng_CHneg_dz_beamspot;
  vector<Double_t> fTwoProng_CHneg_dxy;
  vector<Double_t> fTwoProng_CHneg_dxy_PV;
  vector<Double_t> fTwoProng_CHneg_dxy_beamspot;
  vector<Double_t> fTwoProng_isoPF_vz;
  vector<Double_t> fTwoProng_isoPF_vx;
  vector<Double_t> fTwoProng_isoPF_vy;
  vector<Double_t> fTwoProng_isoPF_dz;
  vector<Double_t> fTwoProng_isoPF_dz_PV;
  vector<Double_t> fTwoProng_isoPF_dz_beamspot;
  vector<Double_t> fTwoProng_isoPF_dxy;
  vector<Double_t> fTwoProng_isoPF_dxy_PV;
  vector<Double_t> fTwoProng_isoPF_dxy_beamspot;
  vector<Double_t> fTwoProng_photon_pt;
  vector<Double_t> fTwoProng_photon_eta;
  vector<Double_t> fTwoProng_photon_phi;
  vector<Double_t> fTwoProng_photon_mass;
  vector<Double_t> fTwoProng_photon_nGamma;
  vector<Double_t> fTwoProng_photon_nElectron;
  vector<Double_t> fTwoProng_chargedIso;
  vector<Double_t> fTwoProng_neutralIso;
  vector<Double_t> fTwoProng_egammaIso;
  vector<Int_t> fTwoProng_nChargedIsoCone;
  vector<Int_t> fTwoProng_nNeutralIsoCone;
  vector<Int_t> fTwoProng_nEGammaIsoCone;
  vector<Double_t> fTwoProng_genDR;
  vector<Bool_t> fTwoProng_match;

  vector<Double_t> fTwoProngLoose_pt;
  vector<Double_t> fTwoProngLoose_eta;
  vector<Double_t> fTwoProngLoose_phi;
  vector<Double_t> fTwoProngLoose_mass;
  vector<Double_t> fTwoProngLoose_Mass;
  vector<Double_t> fTwoProngLoose_px;
  vector<Double_t> fTwoProngLoose_py;
  vector<Double_t> fTwoProngLoose_pz;
  vector<Double_t> fTwoProngLoose_energy;
  vector<Double_t> fTwoProngLoose_CHpos_pt;
  vector<Double_t> fTwoProngLoose_CHpos_eta;
  vector<Double_t> fTwoProngLoose_CHpos_phi;
  vector<Double_t> fTwoProngLoose_CHpos_mass;
  vector<Double_t> fTwoProngLoose_CHneg_pt;
  vector<Double_t> fTwoProngLoose_CHneg_eta;
  vector<Double_t> fTwoProngLoose_CHneg_phi;
  vector<Double_t> fTwoProngLoose_CHneg_mass;
  vector<Double_t> fTwoProngLoose_CHpos_vz;
  vector<Double_t> fTwoProngLoose_CHpos_vx;
  vector<Double_t> fTwoProngLoose_CHpos_vy;
  vector<Double_t> fTwoProngLoose_CHpos_dz;
  vector<Double_t> fTwoProngLoose_CHpos_dz_PV;
  vector<Double_t> fTwoProngLoose_CHpos_dz_beamspot;
  vector<Double_t> fTwoProngLoose_CHpos_dxy;
  vector<Double_t> fTwoProngLoose_CHpos_dxy_PV;
  vector<Double_t> fTwoProngLoose_CHpos_dxy_beamspot;
  vector<Double_t> fTwoProngLoose_CHneg_vz;
  vector<Double_t> fTwoProngLoose_CHneg_vx;
  vector<Double_t> fTwoProngLoose_CHneg_vy;
  vector<Double_t> fTwoProngLoose_CHneg_dz;
  vector<Double_t> fTwoProngLoose_CHneg_dz_PV;
  vector<Double_t> fTwoProngLoose_CHneg_dz_beamspot;
  vector<Double_t> fTwoProngLoose_CHneg_dxy;
  vector<Double_t> fTwoProngLoose_CHneg_dxy_PV;
  vector<Double_t> fTwoProngLoose_CHneg_dxy_beamspot;
  vector<Double_t> fTwoProngLoose_isoPF_vz;
  vector<Double_t> fTwoProngLoose_isoPF_vx;
  vector<Double_t> fTwoProngLoose_isoPF_vy;
  vector<Double_t> fTwoProngLoose_isoPF_dz;
  vector<Double_t> fTwoProngLoose_isoPF_dz_PV;
  vector<Double_t> fTwoProngLoose_isoPF_dz_beamspot;
  vector<Double_t> fTwoProngLoose_isoPF_dxy;
  vector<Double_t> fTwoProngLoose_isoPF_dxy_PV;
  vector<Double_t> fTwoProngLoose_isoPF_dxy_beamspot;
  vector<Double_t> fTwoProngLoose_photon_pt;
  vector<Double_t> fTwoProngLoose_photon_eta;
  vector<Double_t> fTwoProngLoose_photon_phi;
  vector<Double_t> fTwoProngLoose_photon_mass;
  vector<Double_t> fTwoProngLoose_photon_nGamma;
  vector<Double_t> fTwoProngLoose_photon_nElectron;
  vector<Double_t> fTwoProngLoose_chargedIso;
  vector<Double_t> fTwoProngLoose_neutralIso;
  vector<Double_t> fTwoProngLoose_egammaIso;
  vector<Int_t> fTwoProngLoose_nChargedIsoCone;
  vector<Int_t> fTwoProngLoose_nNeutralIsoCone;
  vector<Int_t> fTwoProngLoose_nEGammaIsoCone;
  vector<Double_t> fTwoProngLoose_genDR;
  vector<Bool_t> fTwoProngLoose_match;

  vector<Double_t> fGenPhi_pt;
  vector<Double_t> fGenPhi_eta;
  vector<Double_t> fGenPhi_phi;
  vector<Double_t> fGenPhi_mass;
  vector<Double_t> fGenPhi_px;
  vector<Double_t> fGenPhi_py;
  vector<Double_t> fGenPhi_pz;
  vector<Double_t> fGenPhi_energy;
  vector<Double_t> fGenEta_pt;
  vector<Double_t> fGenEta_eta;
  vector<Double_t> fGenEta_phi;
  vector<Double_t> fGenEta_mass;
  vector<Double_t> fGenEta_px;
  vector<Double_t> fGenEta_py;
  vector<Double_t> fGenEta_pz;
  vector<Double_t> fGenEta_energy;
  vector<Double_t> fGenEta_candDR;
  vector<Double_t> fGenEta_passedCandDR;
  vector<Double_t> fGenEta_jetDR;
  int fGen_decayType;

  ExoDiPhotons::recoPhotonInfo_t fRecoTightPhotonInfo1;
  ExoDiPhotons::recoPhotonInfo_t fRecoTightPhotonInfo2;
  ExoDiPhotons::recoPhotonInfo_t fRecoTightPhotonInfo3;
  ExoDiPhotons::recoPhotonInfo_t fRecoTightPhotonInfo1_v2;
  ExoDiPhotons::recoPhotonInfo_t fRecoTightPhotonInfo2_v2;
  ExoDiPhotons::recoPhotonInfo_t fRecoTightPhotonInfo3_v2;

  ExoDiPhotons::recoDiObjectInfo_t fTwoProngTwoProngInfo; 
  ExoDiPhotons::recoDiObjectInfo_t fGammaTwoProngInfo; 
  ExoDiPhotons::recoDiObjectInfo_t fTwoProngFakeInfo; 
  ExoDiPhotons::recoDiObjectInfo_t fGammaFakeInfo; 
  ExoDiPhotons::recoDiObjectInfo_t fGammaGammaInfo; 
};

//
// constants, enums and typedefs
// 


//
// static data member definitions
//

//
// constructors and destructor
//
ExoDiPhotonAnalyzer::ExoDiPhotonAnalyzer(const edm::ParameterSet& iConfig)
  : fisMC(iConfig.getUntrackedParameter<bool>("isMC")),
    fisSignal(iConfig.getUntrackedParameter<bool>("isSignal")),
    fDebug(iConfig.getUntrackedParameter<bool>("debug")),
    fchargedDecayCutflow(iConfig.getUntrackedParameter<bool>("chargedDecayCutflow")),
    fCandidatePairDR(iConfig.getUntrackedParameter<double>("chargedHadronPairMinDeltaR")),
    fCandidatePairMinPt(iConfig.getUntrackedParameter<double>("chargedHadronMinPt")),
    fCandidatePairIsolationDR(iConfig.getUntrackedParameter<double>("isolationConeR")),
    fCandidatePairPhiBox(iConfig.getUntrackedParameter<double>("photonPhiBoxSize")),
    fCandidatePairEtaBox(iConfig.getUntrackedParameter<double>("photonEtaBoxSize")),
    fCandidatePairPhotonPtCut(iConfig.getUntrackedParameter<double>("photonPtCut")),
    fCandidatePairChargedIsoCut(iConfig.getUntrackedParameter<double>("chargedIsoCut")),
    fCandidatePairNeutralIsoCut(iConfig.getUntrackedParameter<double>("neutralIsoCut")),
    fCandidatePairEGammaIsoCut(iConfig.getUntrackedParameter<double>("egammaIsoCut")),
    fCandidatePairChargedIsoFakeCut(iConfig.getUntrackedParameter<double>("chargedIsoFakeMax")),
    fCandidatePairNeutralIsoFakeCut(iConfig.getUntrackedParameter<double>("neutralIsoFakeMax")),
    fCandidatePairEGammaIsoFakeCut(iConfig.getUntrackedParameter<double>("egammaIsoFakeMax")),
    fCandidatePairGenMatchDR(iConfig.getUntrackedParameter<double>("generatorEtaMatchDR")),
    fTwoProngFakeRateCalcOnly(iConfig.getUntrackedParameter<bool>("noTreeOnlyFakeRateHistos")),

    rhoToken_(consumes<double> (iConfig.getParameter<edm::InputTag>("rho")))
{
  // setting requested by Steve, omit all Brandon's code
  fOmitChargedDecayCode = ( iConfig.exists("omitChargedDecayCode") ? iConfig.getParameter<bool>("omitChargedDecayCode") : false );  
  fSkipPhotonMCCode = ( iConfig.exists("skipPhotonMCCode") ? iConfig.getParameter<bool>("skipPhotonMCCode") : true );  

  // Initialize cutflow variables
  fCutflow_total = 0;
  fCutflow_oneCand = 0;
  fCutflow_twoCand = 0;
  fCutflow_onePass = 0;
  fCutflow_twoPass = 0;
  fCutflow_oneMatch = 0;
  fCutflow_twoMatch = 0;
  fCutflow_onePassMatch = 0;
  fCutflow_twoPassMatch = 0;
  fCutflow_passCharged = 0;
  fCutflow_passNeutral = 0;
  fCutflow_passEGamma = 0;
  fCutflow_passPhotonPt = 0;

  pfcandsToken_ = consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates"));
  genToken_ = consumes<vector<reco::GenParticle>>(edm::InputTag("prunedGenParticles"));
  pvToken_ = consumes<vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));
  beamToken_ = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  ak4Token_ = consumes<std::vector<pat::Jet>>(edm::InputTag("slimmedJets"));
  photonToken_ = consumes<std::vector<pat::Photon>>(edm::InputTag("slimmedPhotons"));
  metToken_ = consumes<std::vector<pat::MET>>(edm::InputTag("slimmedMETs"));
  electronToken_ = consumes<std::vector<pat::Electron>>(edm::InputTag("slimmedElectrons"));
  muonToken_ = consumes<std::vector<pat::Muon>>(edm::InputTag("slimmedMuons"));

  edm::Service<TFileService> fs;
  // Branches for charged decay analysis
  if (!fOmitChargedDecayCode) {
  fTree2 = fs->make<TTree>("fTree2","ChargedDecayTree");
  // Event wide
  fTree2->Branch("eventNum",&fEventNum,"eventNum/I");
  fTree2->Branch("runNum",&fRunNum,"runNum/I");
  fTree2->Branch("lumiNum",&fLumiNum,"lumiNum/I");
  fTree2->Branch("nPV",&fNumPVs,"nPV/I");
  fTree2->Branch("nPF",&fNumPF,"nPF/I");
  fTree2->Branch("nPrunedPF",&fNumPrunedPF,"numPrunedPF/I");
  fTree2->Branch("HT",&fHT,"HT/D");
  fTree2->Branch("MET",&fMET,"MET/D");
  fTree2->Branch("MET_phi",&fMET_phi,"MET_phi/D");
  fTree2->Branch("nElectrons",&fNumElectrons,"nElectrons/I");
  fTree2->Branch("nMuons",&fNumMuons,"nMuons/I");
  // Jets
  fTree2->Branch("nJets",&fNumAK4jets,"nJets/I");
  fTree2->Branch("jet_pt",&fAK4jet_pt);
  fTree2->Branch("jet_eta",&fAK4jet_eta);
  fTree2->Branch("jet_phi",&fAK4jet_phi);
  fTree2->Branch("jet_mass",&fAK4jet_mass);
  fTree2->Branch("jet_px",&fAK4jet_px);
  fTree2->Branch("jet_py",&fAK4jet_py);
  fTree2->Branch("jet_pz",&fAK4jet_pz);
  fTree2->Branch("jet_energy",&fAK4jet_energy);
  // Photons
  fTree2->Branch("nPhotons",&fNumPhotons,"nPhotons/I");
  fTree2->Branch("photon_pt",&fPhoton_pt);
  fTree2->Branch("photon_eta",&fPhoton_eta);
  fTree2->Branch("photon_phi",&fPhoton_phi);
  fTree2->Branch("photon_mass",&fPhoton_mass);
  fTree2->Branch("photon_px",&fPhoton_px);
  fTree2->Branch("photon_py",&fPhoton_py);
  fTree2->Branch("photon_pz",&fPhoton_pz);
  fTree2->Branch("photon_energy",&fPhoton_energy);
  // Cutflow
  fTree2->Branch("nCands",&fNumTwoProng,"nCands/I");
  fTree2->Branch("nPass",&fNumTwoProngPass,"nPass/I");
  fTree2->Branch("nMatched",&fNumTwoProngMatched,"nMatched/I");
  fTree2->Branch("nPassChargedIso",&fNumTwoProngPassChargedIso,"nPassChargedIso/I");
  fTree2->Branch("nPassNeutralIso",&fNumTwoProngPassNeutralIso,"nPassNeutralIso/I");
  fTree2->Branch("nPassEGammaIso",&fNumTwoProngPassEGammaIso,"nPassEGammaIso/I");
  fTree2->Branch("nPassPhotonPtIso",&fNumTwoProngPassChargedIso,"nPassPhotonIso/I");
  fTree2->Branch("nFakes",&fNumTwoProngFake,"nFakes/I");
  fTree2->Branch("nTightPhotons",&fNumTightPhotons,"nTightPhotons/I");
  fTree2->Branch("nTightPhotons_v2",&fNumTightPhotons_v2,"nTightPhotons_v2/I");
  fTree2->Branch("nOffset",&fNumOffset,"nOffset/I");
  // Candidate information
  /*fTree2->Branch("Cand_pt",&fCand_pt);
  fTree2->Branch("Cand_eta",&fCand_eta);
  fTree2->Branch("Cand_phi",&fCand_phi);
  fTree2->Branch("Cand_mass",&fCand_mass);
  fTree2->Branch("Cand_Mass",&fCand_Mass);
  fTree2->Branch("Cand_CHpos_pt",&fCand_CHpos_pt);
  fTree2->Branch("Cand_CHpos_eta",&fCand_CHpos_eta);
  fTree2->Branch("Cand_CHpos_phi",&fCand_CHpos_phi);
  fTree2->Branch("Cand_CHpos_mass",&fCand_CHpos_mass);
  fTree2->Branch("Cand_CHneg_pt",&fCand_CHneg_pt);
  fTree2->Branch("Cand_CHneg_eta",&fCand_CHneg_eta);
  fTree2->Branch("Cand_CHneg_phi",&fCand_CHneg_phi);
  fTree2->Branch("Cand_CHneg_mass",&fCand_CHneg_mass);
  fTree2->Branch("Cand_photon_pt",&fCand_photon_pt);
  fTree2->Branch("Cand_photon_eta",&fCand_photon_eta);
  fTree2->Branch("Cand_photon_phi",&fCand_photon_phi);
  fTree2->Branch("Cand_photon_mass",&fCand_photon_mass);
  fTree2->Branch("Cand_photon_nGamma",&fCand_photon_nGamma);
  fTree2->Branch("Cand_photon_nElectron",&fCand_photon_nElectron);
  fTree2->Branch("Cand_chargedIso",&fCand_chargedIso);
  fTree2->Branch("Cand_neutralIso",&fCand_neutralIso);
  fTree2->Branch("Cand_egammaIso",&fCand_egammaIso);
  fTree2->Branch("Cand_CHpos_vz",&fCand_CHpos_vz);
  fTree2->Branch("Cand_CHpos_vx",&fCand_CHpos_vx);
  fTree2->Branch("Cand_CHpos_vy",&fCand_CHpos_vy);
  fTree2->Branch("Cand_CHpos_dz",&fCand_CHpos_dz);
  fTree2->Branch("Cand_CHpos_dz_PV",&fCand_CHpos_dz_PV);
  fTree2->Branch("Cand_CHpos_dz_beamspot",&fCand_CHpos_dz_beamspot);
  fTree2->Branch("Cand_CHpos_dxy",&fCand_CHpos_dxy);
  fTree2->Branch("Cand_CHpos_dxy_PV",&fCand_CHpos_dxy_PV);
  fTree2->Branch("Cand_CHpos_dxy_beamspot",&fCand_CHpos_dxy_beamspot);
  fTree2->Branch("Cand_CHneg_vz",&fCand_CHneg_vz);
  fTree2->Branch("Cand_CHneg_vx",&fCand_CHneg_vx);
  fTree2->Branch("Cand_CHneg_vy",&fCand_CHneg_vy);
  fTree2->Branch("Cand_CHneg_dz",&fCand_CHneg_dz);
  fTree2->Branch("Cand_CHneg_dz_PV",&fCand_CHneg_dz_PV);
  fTree2->Branch("Cand_CHneg_dz_beamspot",&fCand_CHneg_dz_beamspot);
  fTree2->Branch("Cand_CHneg_dxy",&fCand_CHneg_dxy);
  fTree2->Branch("Cand_CHneg_dxy_PV",&fCand_CHneg_dxy_PV);
  fTree2->Branch("Cand_CHneg_dxy_beamspot",&fCand_CHneg_dxy_beamspot);
  fTree2->Branch("Cand_isoPF_vz",&fCand_isoPF_vz);
  fTree2->Branch("Cand_isoPF_vx",&fCand_isoPF_vx);
  fTree2->Branch("Cand_isoPF_vy",&fCand_isoPF_vy);
  fTree2->Branch("Cand_isoPF_dz",&fCand_isoPF_dz);
  fTree2->Branch("Cand_isoPF_dz_PV",&fCand_isoPF_dz_PV);
  fTree2->Branch("Cand_isoPF_dz_beamspot",&fCand_isoPF_dz_beamspot);
  fTree2->Branch("Cand_isoPF_dxy",&fCand_isoPF_dxy);
  fTree2->Branch("Cand_isoPF_dxy_PV",&fCand_isoPF_dxy_PV);
  fTree2->Branch("Cand_isoPF_dxy_beamspot",&fCand_isoPF_dxy_beamspot);
  fTree2->Branch("Cand_nChargedIsoCone",&fCand_nChargedIsoCone);
  fTree2->Branch("Cand_nNeutralIsoCone",&fCand_nNeutralIsoCone);
  fTree2->Branch("Cand_nEGammaIsoCone",&fCand_nEGammaIsoCone);
  fTree2->Branch("Cand_genDR",&fCand_genDR);
  fTree2->Branch("Cand_pass",&fCand_tight);
  fTree2->Branch("Cand_passChargedIso",&fCand_passChargedIso);
  fTree2->Branch("Cand_passNeutralIso",&fCand_passNeutralIso);
  fTree2->Branch("Cand_passEGammaIso",&fCand_passEGammaIso);
  fTree2->Branch("Cand_passPhotonPt",&fCand_passPhotonPt);
  fTree2->Branch("Cand_fake",&fCand_loose);*/
  // Loose Candidate information, sorted by pt
  fTree2->Branch("TwoProngLoose_pt",&fTwoProngLoose_pt);
  fTree2->Branch("TwoProngLoose_eta",&fTwoProngLoose_eta);
  fTree2->Branch("TwoProngLoose_phi",&fTwoProngLoose_phi);
  fTree2->Branch("TwoProngLoose_mass",&fTwoProngLoose_mass);
  fTree2->Branch("TwoProngLoose_Mass",&fTwoProngLoose_Mass);
  fTree2->Branch("TwoProngLoose_px",&fTwoProngLoose_px);
  fTree2->Branch("TwoProngLoose_py",&fTwoProngLoose_py);
  fTree2->Branch("TwoProngLoose_pz",&fTwoProngLoose_pz);
  fTree2->Branch("TwoProngLoose_energy",&fTwoProngLoose_energy);
  fTree2->Branch("TwoProngLoose_CHpos_pt",&fTwoProngLoose_CHpos_pt);
  fTree2->Branch("TwoProngLoose_CHpos_eta",&fTwoProngLoose_CHpos_eta);
  fTree2->Branch("TwoProngLoose_CHpos_phi",&fTwoProngLoose_CHpos_phi);
  fTree2->Branch("TwoProngLoose_CHpos_mass",&fTwoProngLoose_CHpos_mass);
  fTree2->Branch("TwoProngLoose_CHneg_pt",&fTwoProngLoose_CHneg_pt);
  fTree2->Branch("TwoProngLoose_CHneg_eta",&fTwoProngLoose_CHneg_eta);
  fTree2->Branch("TwoProngLoose_CHneg_phi",&fTwoProngLoose_CHneg_phi);
  fTree2->Branch("TwoProngLoose_CHneg_mass",&fTwoProngLoose_CHneg_mass);
  fTree2->Branch("TwoProngLoose_photon_pt",&fTwoProngLoose_photon_pt);
  fTree2->Branch("TwoProngLoose_photon_eta",&fTwoProngLoose_photon_eta);
  fTree2->Branch("TwoProngLoose_photon_phi",&fTwoProngLoose_photon_phi);
  fTree2->Branch("TwoProngLoose_photon_mass",&fTwoProngLoose_photon_mass);
  fTree2->Branch("TwoProngLoose_photon_nGamma",&fTwoProngLoose_photon_nGamma);
  fTree2->Branch("TwoProngLoose_photon_nElectron",&fTwoProngLoose_photon_nElectron);
  fTree2->Branch("TwoProngLoose_chargedIso",&fTwoProngLoose_chargedIso);
  fTree2->Branch("TwoProngLoose_neutralIso",&fTwoProngLoose_neutralIso);
  fTree2->Branch("TwoProngLoose_egammaIso",&fTwoProngLoose_egammaIso);
  fTree2->Branch("TwoProngLoose_CHpos_vz",&fTwoProngLoose_CHpos_vz);
  fTree2->Branch("TwoProngLoose_CHpos_vx",&fTwoProngLoose_CHpos_vx);
  fTree2->Branch("TwoProngLoose_CHpos_vy",&fTwoProngLoose_CHpos_vy);
  fTree2->Branch("TwoProngLoose_CHpos_dz",&fTwoProngLoose_CHpos_dz);
  fTree2->Branch("TwoProngLoose_CHpos_dz_PV",&fTwoProngLoose_CHpos_dz_PV);
  fTree2->Branch("TwoProngLoose_CHpos_dz_beamspot",&fTwoProngLoose_CHpos_dz_beamspot);
  fTree2->Branch("TwoProngLoose_CHpos_dxy",&fTwoProngLoose_CHpos_dxy);
  fTree2->Branch("TwoProngLoose_CHpos_dxy_PV",&fTwoProngLoose_CHpos_dxy_PV);
  fTree2->Branch("TwoProngLoose_CHpos_dxy_beamspot",&fTwoProngLoose_CHpos_dxy_beamspot);
  fTree2->Branch("TwoProngLoose_CHneg_vz",&fTwoProngLoose_CHneg_vz);
  fTree2->Branch("TwoProngLoose_CHneg_vx",&fTwoProngLoose_CHneg_vx);
  fTree2->Branch("TwoProngLoose_CHneg_vy",&fTwoProngLoose_CHneg_vy);
  fTree2->Branch("TwoProngLoose_CHneg_dz",&fTwoProngLoose_CHneg_dz);
  fTree2->Branch("TwoProngLoose_CHneg_dz_PV",&fTwoProngLoose_CHneg_dz_PV);
  fTree2->Branch("TwoProngLoose_CHneg_dz_beamspot",&fTwoProngLoose_CHneg_dz_beamspot);
  fTree2->Branch("TwoProngLoose_CHneg_dxy",&fTwoProngLoose_CHneg_dxy);
  fTree2->Branch("TwoProngLoose_CHneg_dxy_PV",&fTwoProngLoose_CHneg_dxy_PV);
  fTree2->Branch("TwoProngLoose_CHneg_dxy_beamspot",&fTwoProngLoose_CHneg_dxy_beamspot);
  fTree2->Branch("TwoProngLoose_isoPF_vz",&fTwoProngLoose_isoPF_vz);
  fTree2->Branch("TwoProngLoose_isoPF_vx",&fTwoProngLoose_isoPF_vx);
  fTree2->Branch("TwoProngLoose_isoPF_vy",&fTwoProngLoose_isoPF_vy);
  fTree2->Branch("TwoProngLoose_isoPF_dz",&fTwoProngLoose_isoPF_dz);
  fTree2->Branch("TwoProngLoose_isoPF_dz_PV",&fTwoProngLoose_isoPF_dz_PV);
  fTree2->Branch("TwoProngLoose_isoPF_dz_beamspot",&fTwoProngLoose_isoPF_dz_beamspot);
  fTree2->Branch("TwoProngLoose_isoPF_dxy",&fTwoProngLoose_isoPF_dxy);
  fTree2->Branch("TwoProngLoose_isoPF_dxy_PV",&fTwoProngLoose_isoPF_dxy_PV);
  fTree2->Branch("TwoProngLoose_isoPF_dxy_beamspot",&fTwoProngLoose_isoPF_dxy_beamspot);
  fTree2->Branch("TwoProngLoose_genDR",&fTwoProngLoose_genDR);
  // Tight Candidate information, sorted by pt
  fTree2->Branch("TwoProng_pt",&fTwoProng_pt);
  fTree2->Branch("TwoProng_eta",&fTwoProng_eta);
  fTree2->Branch("TwoProng_phi",&fTwoProng_phi);
  fTree2->Branch("TwoProng_mass",&fTwoProng_mass);
  fTree2->Branch("TwoProng_Mass",&fTwoProng_Mass);
  fTree2->Branch("TwoProng_px",&fTwoProng_px);
  fTree2->Branch("TwoProng_py",&fTwoProng_py);
  fTree2->Branch("TwoProng_pz",&fTwoProng_pz);
  fTree2->Branch("TwoProng_energy",&fTwoProng_energy);
  fTree2->Branch("TwoProng_CHpos_pt",&fTwoProng_CHpos_pt);
  fTree2->Branch("TwoProng_CHpos_eta",&fTwoProng_CHpos_eta);
  fTree2->Branch("TwoProng_CHpos_phi",&fTwoProng_CHpos_phi);
  fTree2->Branch("TwoProng_CHpos_mass",&fTwoProng_CHpos_mass);
  fTree2->Branch("TwoProng_CHneg_pt",&fTwoProng_CHneg_pt);
  fTree2->Branch("TwoProng_CHneg_eta",&fTwoProng_CHneg_eta);
  fTree2->Branch("TwoProng_CHneg_phi",&fTwoProng_CHneg_phi);
  fTree2->Branch("TwoProng_CHneg_mass",&fTwoProng_CHneg_mass);
  fTree2->Branch("TwoProng_photon_pt",&fTwoProng_photon_pt);
  fTree2->Branch("TwoProng_photon_eta",&fTwoProng_photon_eta);
  fTree2->Branch("TwoProng_photon_phi",&fTwoProng_photon_phi);
  fTree2->Branch("TwoProng_photon_mass",&fTwoProng_photon_mass);
  fTree2->Branch("TwoProng_photon_nGamma",&fTwoProng_photon_nGamma);
  fTree2->Branch("TwoProng_photon_nElectron",&fTwoProng_photon_nElectron);
  fTree2->Branch("TwoProng_chargedIso",&fTwoProng_chargedIso);
  fTree2->Branch("TwoProng_neutralIso",&fTwoProng_neutralIso);
  fTree2->Branch("TwoProng_egammaIso",&fTwoProng_egammaIso);
  fTree2->Branch("TwoProng_CHpos_vz",&fTwoProng_CHpos_vz);
  fTree2->Branch("TwoProng_CHpos_vx",&fTwoProng_CHpos_vx);
  fTree2->Branch("TwoProng_CHpos_vy",&fTwoProng_CHpos_vy);
  fTree2->Branch("TwoProng_CHpos_dz",&fTwoProng_CHpos_dz);
  fTree2->Branch("TwoProng_CHpos_dz_PV",&fTwoProng_CHpos_dz_PV);
  fTree2->Branch("TwoProng_CHpos_dz_beamspot",&fTwoProng_CHpos_dz_beamspot);
  fTree2->Branch("TwoProng_CHpos_dxy",&fTwoProng_CHpos_dxy);
  fTree2->Branch("TwoProng_CHpos_dxy_PV",&fTwoProng_CHpos_dxy_PV);
  fTree2->Branch("TwoProng_CHpos_dxy_beamspot",&fTwoProng_CHpos_dxy_beamspot);
  fTree2->Branch("TwoProng_CHneg_vz",&fTwoProng_CHneg_vz);
  fTree2->Branch("TwoProng_CHneg_vx",&fTwoProng_CHneg_vx);
  fTree2->Branch("TwoProng_CHneg_vy",&fTwoProng_CHneg_vy);
  fTree2->Branch("TwoProng_CHneg_dz",&fTwoProng_CHneg_dz);
  fTree2->Branch("TwoProng_CHneg_dz_PV",&fTwoProng_CHneg_dz_PV);
  fTree2->Branch("TwoProng_CHneg_dz_beamspot",&fTwoProng_CHneg_dz_beamspot);
  fTree2->Branch("TwoProng_CHneg_dxy",&fTwoProng_CHneg_dxy);
  fTree2->Branch("TwoProng_CHneg_dxy_PV",&fTwoProng_CHneg_dxy_PV);
  fTree2->Branch("TwoProng_CHneg_dxy_beamspot",&fTwoProng_CHneg_dxy_beamspot);
  fTree2->Branch("TwoProng_isoPF_vz",&fTwoProng_isoPF_vz);
  fTree2->Branch("TwoProng_isoPF_vx",&fTwoProng_isoPF_vx);
  fTree2->Branch("TwoProng_isoPF_vy",&fTwoProng_isoPF_vy);
  fTree2->Branch("TwoProng_isoPF_dz",&fTwoProng_isoPF_dz);
  fTree2->Branch("TwoProng_isoPF_dz_PV",&fTwoProng_isoPF_dz_PV);
  fTree2->Branch("TwoProng_isoPF_dz_beamspot",&fTwoProng_isoPF_dz_beamspot);
  fTree2->Branch("TwoProng_isoPF_dxy",&fTwoProng_isoPF_dxy);
  fTree2->Branch("TwoProng_isoPF_dxy_PV",&fTwoProng_isoPF_dxy_PV);
  fTree2->Branch("TwoProng_isoPF_dxy_beamspot",&fTwoProng_isoPF_dxy_beamspot);
  fTree2->Branch("TwoProng_genDR",&fTwoProng_genDR);
  // Tight Photons, sorted by pt
  fTree2->Branch("Photon1",&fRecoTightPhotonInfo1,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  fTree2->Branch("Photon2",&fRecoTightPhotonInfo2,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  fTree2->Branch("Photon3",&fRecoTightPhotonInfo3,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  fTree2->Branch("Photon1_v2",&fRecoTightPhotonInfo1_v2,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  fTree2->Branch("Photon2_v2",&fRecoTightPhotonInfo2_v2,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  fTree2->Branch("Photon3_v2",&fRecoTightPhotonInfo3_v2,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  // Generator Objects
  fTree2->Branch("GenPhi_pt",&fGenPhi_pt);
  fTree2->Branch("GenPhi_eta",&fGenPhi_eta);
  fTree2->Branch("GenPhi_phi",&fGenPhi_phi);
  fTree2->Branch("GenPhi_mass",&fGenPhi_mass);
  fTree2->Branch("GenPhi_px",&fGenPhi_px);
  fTree2->Branch("GenPhi_py",&fGenPhi_py);
  fTree2->Branch("GenPhi_pz",&fGenPhi_pz);
  fTree2->Branch("GenPhi_energy",&fGenPhi_energy);
  fTree2->Branch("GenEta_pt",&fGenEta_pt);
  fTree2->Branch("GenEta_eta",&fGenEta_eta);
  fTree2->Branch("GenEta_phi",&fGenEta_phi);
  fTree2->Branch("GenEta_mass",&fGenEta_mass);
  fTree2->Branch("GenEta_px",&fGenEta_px);
  fTree2->Branch("GenEta_py",&fGenEta_py);
  fTree2->Branch("GenEta_pz",&fGenEta_pz);
  fTree2->Branch("GenEta_energy",&fGenEta_energy); 
  fTree2->Branch("GenEta_candDR",&fGenEta_candDR); 
  fTree2->Branch("GenEta_passedCandDR",&fGenEta_passedCandDR); 
  fTree2->Branch("GenEta_jetDR",&fGenEta_jetDR); 
  fTree2->Branch("Gen_decayType",&fGen_decayType);
  // Combined Objects
  fTree2->Branch("TwoProngTwoProng",&fTwoProngTwoProngInfo,ExoDiPhotons::recoDiObjectBranchDefString.c_str());
  fTree2->Branch("TwoProngFake",&fTwoProngFakeInfo,ExoDiPhotons::recoDiObjectBranchDefString.c_str());
  fTree2->Branch("GammaTwoProng",&fGammaTwoProngInfo,ExoDiPhotons::recoDiObjectBranchDefString.c_str());
  fTree2->Branch("GammaFake",&fGammaFakeInfo,ExoDiPhotons::recoDiObjectBranchDefString.c_str());
  fTree2->Branch("GammaGamma",&fGammaGammaInfo,ExoDiPhotons::recoDiObjectBranchDefString.c_str());
  // Fake rate histograms
  int num_even_mass_bins = 9;
  double even_mass_bins_array[num_even_mass_bins+1] = {0,.400,.600,.800,1.000,1.200,1.400,1.600,1.800,2.000};
  double *even_mass_bins = even_mass_bins_array;
  int num_odd_mass_bins = 9;
  double odd_mass_bins_array[num_odd_mass_bins+1] = {0,.300,.500,.700,.900,1.100,1.300,1.500,1.700,1.900};
  double *odd_mass_bins = odd_mass_bins_array;

  int num_pt_bins = 26;
  double pt_low = 0;
  double pt_high = 1300;
  int num_eta_bins = 20;
  double eta_low = -5;
  double eta_high = 5;
  int num_phi_bins = 24;
  double phi_low = -3.15;
  double phi_high = 3.15;

  fTwoProngFakeNume_even_pt = fs->make<TH2F>("twoprongfakenume_even_pt","Fake Numerator count, pt vs even mass bins",num_pt_bins,pt_low,pt_high,num_even_mass_bins,even_mass_bins);
  fTwoProngFakeDeno_even_pt = fs->make<TH2F>("twoprongfakedeno_even_pt","Fake Denominator count, pt vs even mass bins",num_pt_bins,pt_low,pt_high,num_even_mass_bins,even_mass_bins);
  fTwoProngFakeRate_even_pt = fs->make<TH2F>("twoprongfakerate_even_pt","Fake Rate, pt vs even mass bins",num_pt_bins,pt_low,pt_high,num_even_mass_bins,even_mass_bins);
  fTwoProngFakeNume_odd_pt = fs->make<TH2F>("twoprongfakenume_odd_pt","Fake Numerator count, pt vs odd mass bins",num_pt_bins,pt_low,pt_high,num_odd_mass_bins,odd_mass_bins);
  fTwoProngFakeDeno_odd_pt = fs->make<TH2F>("twoprongfakedeno_odd_pt","Fake Denominator count, pt vs odd mass bins",num_pt_bins,pt_low,pt_high,num_odd_mass_bins,odd_mass_bins);
  fTwoProngFakeRate_odd_pt = fs->make<TH2F>("twoprongfakerate_odd_pt","Fake Rate, pt vs odd mass bins",num_pt_bins,pt_low,pt_high,num_odd_mass_bins,odd_mass_bins);

  fTwoProngFakeNume_even_eta = fs->make<TH2F>("twoprongfakenume_even_eta","Fake Numerator count, eta vs even mass bins",num_eta_bins,eta_low,eta_high,num_even_mass_bins,even_mass_bins);
  fTwoProngFakeDeno_even_eta = fs->make<TH2F>("twoprongfakedeno_even_eta","Fake Denominator count, eta vs even mass bins",num_eta_bins,eta_low,eta_high,num_even_mass_bins,even_mass_bins);
  fTwoProngFakeRate_even_eta = fs->make<TH2F>("twoprongfakerate_even_eta","Fake Rate, eta vs even mass bins",num_eta_bins,eta_low,eta_high,num_even_mass_bins,even_mass_bins);
  fTwoProngFakeNume_odd_eta = fs->make<TH2F>("twoprongfakenume_odd_eta","Fake Numerator count, eta vs odd mass bins",num_eta_bins,eta_low,eta_high,num_odd_mass_bins,odd_mass_bins);
  fTwoProngFakeDeno_odd_eta = fs->make<TH2F>("twoprongfakedeno_odd_eta","Fake Denominator count, eta vs odd mass bins",num_eta_bins,eta_low,eta_high,num_odd_mass_bins,odd_mass_bins);
  fTwoProngFakeRate_odd_eta = fs->make<TH2F>("twoprongfakerate_odd_eta","Fake Rate, eta vs odd mass bins",num_eta_bins,eta_low,eta_high,num_odd_mass_bins,odd_mass_bins);

  fTwoProngFakeNume_even_phi = fs->make<TH2F>("twoprongfakenume_even_phi","Fake Numerator count, phi vs even mass bins",num_phi_bins,phi_low,phi_high,num_even_mass_bins,even_mass_bins);
  fTwoProngFakeDeno_even_phi = fs->make<TH2F>("twoprongfakedeno_even_phi","Fake Denominator count, phi vs even mass bins",num_phi_bins,phi_low,phi_high,num_even_mass_bins,even_mass_bins);
  fTwoProngFakeRate_even_phi = fs->make<TH2F>("twoprongfakerate_even_phi","Fake Rate, phi vs even mass bins",num_phi_bins,phi_low,phi_high,num_even_mass_bins,even_mass_bins);
  fTwoProngFakeNume_odd_phi = fs->make<TH2F>("twoprongfakenume_odd_phi","Fake Numerator count, phi vs odd mass bins",num_phi_bins,phi_low,phi_high,num_odd_mass_bins,odd_mass_bins);
  fTwoProngFakeDeno_odd_phi = fs->make<TH2F>("twoprongfakedeno_odd_phi","Fake Denominator count, phi vs odd mass bins",num_phi_bins,phi_low,phi_high,num_odd_mass_bins,odd_mass_bins);
  fTwoProngFakeRate_odd_phi = fs->make<TH2F>("twoprongfakerate_odd_phi","Fake Rate, phi vs odd mass bins",num_phi_bins,phi_low,phi_high,num_odd_mass_bins,odd_mass_bins);

  } // done with charged decay code conditional

    recHitsEBTag_ = iConfig.getUntrackedParameter<edm::InputTag>("RecHitsEBTag",edm::InputTag("reducedEgamma:reducedEBRecHits"));
    recHitsEETag_ = iConfig.getUntrackedParameter<edm::InputTag>("RecHitsEETag",edm::InputTag("reducedEgamma:reducedEERecHits"));
    recHitsEBToken = consumes < EcalRecHitCollection > (recHitsEBTag_);
    recHitsEEToken = consumes < EcalRecHitCollection > (recHitsEETag_);
}


ExoDiPhotonAnalyzer::~ExoDiPhotonAnalyzer()
{
 
}

// Helper methods
bool
ExoDiPhotonAnalyzer::isNeutral(int one, int two, int three)
{
  if (one == 22 && two == 22) return true;
  if (two == 22 && three == 22) return true;
  if (one == 22 && three == 22) return true;
  if (one == 111 && two == 111 && three == 111) return true;

  return false;
}

bool
ExoDiPhotonAnalyzer::isCharged(int one, int two, int three)
{
  if (one==211 && two==-211 && three==111) return true;
  if (one==-211 && two==211 && three==111) return true;
  if (one==111 && two==211 && three==-211) return true;
  if (one==111 && two==-211 && three==211) return true;
  if (one==211 && two==111 && three==-211) return true;
  if (one==-211 && two==111 && three==211) return true;
 
  if (one==211 && two==-211 && three==22) return true;
  if (one==-211 && two==211 && three==22) return true;
  if (one==22 && two==211 && three==-211) return true;
  if (one==22 && two==-211 && three==211) return true;
  if (one==211 && two==22 && three==-211) return true;
  if (one==-211 && two==22 && three==211) return true;
 
  return false;
}

void
ExoDiPhotonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;

  if (fDebug) cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;
  fEventNum = iEvent.id().event();
  fRunNum = iEvent.id().run();
  fLumiNum = iEvent.id().luminosityBlock();

  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho_ = *rhoH;

  // ecal information
  lazyTools_ = std::auto_ptr<noZS::EcalClusterLazyTools>( new noZS::EcalClusterLazyTools(iEvent, iSetup, recHitsEBToken, recHitsEEToken));   

  // get ecal barrel recHits for spike rejection
  edm::Handle<EcalRecHitCollection> recHitsEB_h;
  iEvent.getByToken(recHitsEBToken, recHitsEB_h );
  const EcalRecHitCollection * recHitsEB = 0;
  if ( ! recHitsEB_h.isValid() ) {
    LogError("ExoDiPhotonAnalyzer") << " ECAL Barrel RecHit Collection not available !"; return;
  } else {
    recHitsEB = recHitsEB_h.product();
  }

  edm::Handle<EcalRecHitCollection> recHitsEE_h;
  iEvent.getByToken(recHitsEEToken, recHitsEE_h );
  const EcalRecHitCollection * recHitsEE = 0;
  if ( ! recHitsEE_h.isValid() ) {
    LogError("ExoDiPhotonAnalyzer") << " ECAL Endcap RecHit Collection not available !"; return;
  } else {
    recHitsEE = recHitsEE_h.product();
  }

  edm::ESHandle<EcalChannelStatus> chStatus;
  iSetup.get<EcalChannelStatusRcd>().get(chStatus);
  const EcalChannelStatus *ch_status = chStatus.product(); 

  // Charged decay analysis
  if (fDebug) cout << ". start charged decay code" << endl;
  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(pfcandsToken_, pfcands);

  edm::Handle<vector<reco::GenParticle>> genparticles;
  if (fisMC) {
    iEvent.getByToken(genToken_, genparticles);
  }

  edm::Handle<vector<reco::Vertex>> primaryvertecies;
  iEvent.getByToken(pvToken_, primaryvertecies);

  edm::Handle<reco::BeamSpot> beamspot;
  iEvent.getByToken(beamToken_, beamspot);

  edm::Handle<std::vector<pat::Jet>> ak4jets;
  iEvent.getByToken(ak4Token_, ak4jets);

  edm::Handle<std::vector<pat::Photon>> miniaod_photons;
  iEvent.getByToken(photonToken_, miniaod_photons);

  edm::Handle<std::vector<pat::MET>> MET;
  iEvent.getByToken(metToken_, MET);

  edm::Handle<std::vector<pat::Electron>> electrons;
  iEvent.getByToken(electronToken_, electrons);

  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonToken_, muons);

  if (fDebug) cout << ". initializing branches" << endl;
  fCand_pt.clear();
  fCand_eta.clear();
  fCand_phi.clear();
  fCand_mass.clear();
  fCand_Mass.clear();
  fCand_CHpos_pt.clear();
  fCand_CHpos_eta.clear();
  fCand_CHpos_phi.clear();
  fCand_CHpos_mass.clear();
  fCand_CHneg_pt.clear();
  fCand_CHneg_eta.clear();
  fCand_CHneg_phi.clear();
  fCand_CHneg_mass.clear();
  fCand_photon_pt.clear();
  fCand_photon_eta.clear();
  fCand_photon_phi.clear();
  fCand_photon_mass.clear();
  fCand_photon_nGamma.clear();
  fCand_photon_nElectron.clear();
  fCand_chargedIso.clear();
  fCand_neutralIso.clear();
  fCand_egammaIso.clear();
  fCand_CHpos_vz.clear();
  fCand_CHpos_vx.clear();
  fCand_CHpos_vy.clear();
  fCand_CHpos_dz.clear();
  fCand_CHpos_dz_PV.clear();
  fCand_CHpos_dz_beamspot.clear();
  fCand_CHpos_dxy.clear();
  fCand_CHpos_dxy_PV.clear();
  fCand_CHpos_dxy_beamspot.clear();
  fCand_CHneg_vz.clear();
  fCand_CHneg_vx.clear();
  fCand_CHneg_vy.clear();
  fCand_CHneg_dz.clear();
  fCand_CHneg_dz_PV.clear();
  fCand_CHneg_dz_beamspot.clear();
  fCand_CHneg_dxy.clear();
  fCand_CHneg_dxy_PV.clear();
  fCand_CHneg_dxy_beamspot.clear();
  fCand_isoPF_vz.clear();
  fCand_isoPF_vx.clear();
  fCand_isoPF_vy.clear();
  fCand_isoPF_dz.clear();
  fCand_isoPF_dz_PV.clear();
  fCand_isoPF_dz_beamspot.clear();
  fCand_isoPF_dxy.clear();
  fCand_isoPF_dxy_PV.clear();
  fCand_isoPF_dxy_beamspot.clear();
  fCand_nChargedIsoCone.clear();
  fCand_nNeutralIsoCone.clear();
  fCand_nEGammaIsoCone.clear();
  fCand_genDR.clear();
  fCand_tight.clear();
  fCand_passChargedIso.clear();
  fCand_passNeutralIso.clear();
  fCand_passEGammaIso.clear();
  fCand_passPhotonPt.clear();
  fCand_loose.clear();
  fCand_match.clear();

  fTwoProngLoose_pt.clear();
  fTwoProngLoose_eta.clear();
  fTwoProngLoose_phi.clear();
  fTwoProngLoose_px.clear();
  fTwoProngLoose_py.clear();
  fTwoProngLoose_pz.clear();
  fTwoProngLoose_mass.clear();
  fTwoProngLoose_energy.clear();
  fTwoProngLoose_Mass.clear();
  fTwoProngLoose_CHpos_pt.clear();
  fTwoProngLoose_CHpos_eta.clear();
  fTwoProngLoose_CHpos_phi.clear();
  fTwoProngLoose_CHpos_mass.clear();
  fTwoProngLoose_CHneg_pt.clear();
  fTwoProngLoose_CHneg_eta.clear();
  fTwoProngLoose_CHneg_phi.clear();
  fTwoProngLoose_CHneg_mass.clear();
  fTwoProngLoose_photon_pt.clear();
  fTwoProngLoose_photon_eta.clear();
  fTwoProngLoose_photon_phi.clear();
  fTwoProngLoose_photon_mass.clear();
  fTwoProngLoose_photon_nGamma.clear();
  fTwoProngLoose_photon_nElectron.clear();
  fTwoProngLoose_chargedIso.clear();
  fTwoProngLoose_neutralIso.clear();
  fTwoProngLoose_egammaIso.clear();
  fTwoProngLoose_CHpos_vz.clear();
  fTwoProngLoose_CHpos_vx.clear();
  fTwoProngLoose_CHpos_vy.clear();
  fTwoProngLoose_CHpos_dz.clear();
  fTwoProngLoose_CHpos_dz_PV.clear();
  fTwoProngLoose_CHpos_dz_beamspot.clear();
  fTwoProngLoose_CHpos_dxy.clear();
  fTwoProngLoose_CHpos_dxy_PV.clear();
  fTwoProngLoose_CHpos_dxy_beamspot.clear();
  fTwoProngLoose_CHneg_vz.clear();
  fTwoProngLoose_CHneg_vx.clear();
  fTwoProngLoose_CHneg_vy.clear();
  fTwoProngLoose_CHneg_dz.clear();
  fTwoProngLoose_CHneg_dz_PV.clear();
  fTwoProngLoose_CHneg_dz_beamspot.clear();
  fTwoProngLoose_CHneg_dxy.clear();
  fTwoProngLoose_CHneg_dxy_PV.clear();
  fTwoProngLoose_CHneg_dxy_beamspot.clear();
  fTwoProngLoose_isoPF_vz.clear();
  fTwoProngLoose_isoPF_vx.clear();
  fTwoProngLoose_isoPF_vy.clear();
  fTwoProngLoose_isoPF_dz.clear();
  fTwoProngLoose_isoPF_dz_PV.clear();
  fTwoProngLoose_isoPF_dz_beamspot.clear();
  fTwoProngLoose_isoPF_dxy.clear();
  fTwoProngLoose_isoPF_dxy_PV.clear();
  fTwoProngLoose_isoPF_dxy_beamspot.clear();
  fTwoProngLoose_genDR.clear();
  fTwoProngLoose_match.clear();

  fTwoProng_pt.clear();
  fTwoProng_eta.clear();
  fTwoProng_phi.clear();
  fTwoProng_px.clear();
  fTwoProng_py.clear();
  fTwoProng_pz.clear();
  fTwoProng_mass.clear();
  fTwoProng_energy.clear();
  fTwoProng_Mass.clear();
  fTwoProng_CHpos_pt.clear();
  fTwoProng_CHpos_eta.clear();
  fTwoProng_CHpos_phi.clear();
  fTwoProng_CHpos_mass.clear();
  fTwoProng_CHneg_pt.clear();
  fTwoProng_CHneg_eta.clear();
  fTwoProng_CHneg_phi.clear();
  fTwoProng_CHneg_mass.clear();
  fTwoProng_photon_pt.clear();
  fTwoProng_photon_eta.clear();
  fTwoProng_photon_phi.clear();
  fTwoProng_photon_mass.clear();
  fTwoProng_photon_nGamma.clear();
  fTwoProng_photon_nElectron.clear();
  fTwoProng_chargedIso.clear();
  fTwoProng_neutralIso.clear();
  fTwoProng_egammaIso.clear();
  fTwoProng_CHpos_vz.clear();
  fTwoProng_CHpos_vx.clear();
  fTwoProng_CHpos_vy.clear();
  fTwoProng_CHpos_dz.clear();
  fTwoProng_CHpos_dz_PV.clear();
  fTwoProng_CHpos_dz_beamspot.clear();
  fTwoProng_CHpos_dxy.clear();
  fTwoProng_CHpos_dxy_PV.clear();
  fTwoProng_CHpos_dxy_beamspot.clear();
  fTwoProng_CHneg_vz.clear();
  fTwoProng_CHneg_vx.clear();
  fTwoProng_CHneg_vy.clear();
  fTwoProng_CHneg_dz.clear();
  fTwoProng_CHneg_dz_PV.clear();
  fTwoProng_CHneg_dz_beamspot.clear();
  fTwoProng_CHneg_dxy.clear();
  fTwoProng_CHneg_dxy_PV.clear();
  fTwoProng_CHneg_dxy_beamspot.clear();
  fTwoProng_isoPF_vz.clear();
  fTwoProng_isoPF_vx.clear();
  fTwoProng_isoPF_vy.clear();
  fTwoProng_isoPF_dz.clear();
  fTwoProng_isoPF_dz_PV.clear();
  fTwoProng_isoPF_dz_beamspot.clear();
  fTwoProng_isoPF_dxy.clear();
  fTwoProng_isoPF_dxy_PV.clear();
  fTwoProng_isoPF_dxy_beamspot.clear();
  fTwoProng_genDR.clear();
  fTwoProng_match.clear();

  fGenPhi_pt.clear();
  fGenPhi_eta.clear();
  fGenPhi_phi.clear();
  fGenPhi_energy.clear();
  fGenPhi_px.clear();
  fGenPhi_py.clear();
  fGenPhi_pz.clear();
  fGenPhi_energy.clear();
  fGenEta_pt.clear();
  fGenEta_eta.clear();
  fGenEta_phi.clear();
  fGenEta_energy.clear();
  fGenEta_candDR.clear();
  fGenEta_passedCandDR.clear();
  fGenEta_jetDR.clear();
  fGenEta_px.clear();
  fGenEta_py.clear();
  fGenEta_pz.clear();
  fGenEta_energy.clear();

  fAK4jet_pt.clear();
  fAK4jet_eta.clear();
  fAK4jet_phi.clear();
  fAK4jet_mass.clear();
  fAK4jet_px.clear();
  fAK4jet_py.clear();
  fAK4jet_pz.clear();
  fAK4jet_energy.clear();

  fPhoton_pt.clear();
  fPhoton_eta.clear();
  fPhoton_phi.clear();
  fPhoton_mass.clear();
  fPhoton_px.clear();
  fPhoton_py.clear();
  fPhoton_pz.clear();
  fPhoton_energy.clear();

  fNumPVs = primaryvertecies->size();
  fNumPF = pfcands->size();

  // Generator information if Signal
  if (fisSignal && fisMC) {
    for (unsigned int i = 0; i < genparticles->size(); i++) {
      const reco::GenParticle &genparticle = (*genparticles)[i];
      if (genparticle.pdgId() == 9000006 && genparticle.status() == 62) {
        fGenPhi_pt.push_back(genparticle.pt());
        fGenPhi_eta.push_back(genparticle.eta());
        fGenPhi_phi.push_back(genparticle.phi());
        fGenPhi_mass.push_back(genparticle.mass());
        fGenPhi_px.push_back(genparticle.px());
        fGenPhi_py.push_back(genparticle.py());
        fGenPhi_pz.push_back(genparticle.pz());
        fGenPhi_energy.push_back(genparticle.energy());
      }
      if ((genparticle.pdgId() == 221 || genparticle.pdgId() == 331) && genparticle.status() == 2) {
        fGenEta_pt.push_back(genparticle.pt());
        fGenEta_eta.push_back(genparticle.eta());
        fGenEta_phi.push_back(genparticle.phi());
        fGenEta_mass.push_back(genparticle.mass());
        fGenEta_px.push_back(genparticle.px());
        fGenEta_py.push_back(genparticle.py());
        fGenEta_pz.push_back(genparticle.pz());
        fGenEta_energy.push_back(genparticle.energy());
      }
    }
  }

  // Find all pairs of one CH pos and one CH neg within specified DR of each other
  int nPassCharged = 0;
  int nPassNeutral = 0;
  int nPassEGamma = 0;
  int nPassPhotonPt = 0;
  int nFake = 0;
  int pruned_count = 0;
  int nMatch = 0;
  int nPassMatch = 0;
  TLorentzVector LeadingFakeTwoProng; 
  if (fDebug) cout << ". starting candidate loop" << endl;
  for (unsigned int i = 0; i < pfcands->size(); i++) {
    const pat::PackedCandidate &pf1 = (*pfcands)[i];
    if (pf1.pt() < fCandidatePairMinPt) continue;
    if (pf1.fromPV()<=1) continue;
    pruned_count += 1;
    for (unsigned int j = i+1; j < pfcands->size(); j++) { // note loop starting with j=i+1, considers each pair exactly once
      const pat::PackedCandidate &pf2 = (*pfcands)[j];
      if (pf2.pt() < fCandidatePairMinPt) continue;
      if (pf2.fromPV()<=1) continue;
      if (!( ((pf1.pdgId() == 211) && (pf2.pdgId() == -211)) || ((pf1.pdgId() == -211) && (pf2.pdgId() == 211)) )) continue;
      TLorentzVector pfcand1;
      pfcand1.SetPtEtaPhiE(pf1.pt(), pf1.eta(), pf1.phiAtVtx(), pf1.energy());
      TLorentzVector pfcand2;
      pfcand2.SetPtEtaPhiE(pf2.pt(), pf2.eta(), pf2.phiAtVtx(), pf2.energy());
      // found CH pos and CH minus, now check DR
      double dr = pfcand1.DeltaR(pfcand2);
      if (dr < fCandidatePairDR) {
        // Pair is close together, now define photon
        TLorentzVector center;
        center = pfcand1 + pfcand2;
        TLorentzVector photon;
        TLorentzVector leading_pf_photon;
        int numgamma = 0;
        int nume = 0;
        int index_of_leading_pf_photon = -1;
        double pt_of_leading_pf_photon = -1;
        for (unsigned int k = 0; k < pfcands->size(); k++) {
          const pat::PackedCandidate &pf3 = (*pfcands)[k];
          if ((pf3.pdgId() != 22) && (abs(pf3.pdgId()) != 11)) continue; // only pf electron or pf photon contribute to definition of photon
          TLorentzVector pfcand3;
          pfcand3.SetPtEtaPhiE(pf3.pt(), pf3.eta(), pf3.phiAtVtx(), pf3.energy());
          if (fabs(pf3.phiAtVtx() - center.Phi()) < fCandidatePairPhiBox/2.0 &&
              fabs(pf3.eta() - center.Eta()) < fCandidatePairEtaBox/2.0) {
            photon = photon + pfcand3;
            if (pf3.pdgId() == 22) {
              numgamma += 1;
              // find leading photon pf
              if (pf3.pt() > pt_of_leading_pf_photon) {
                pt_of_leading_pf_photon = pf3.pt();
                index_of_leading_pf_photon = k;
              }
            }
            else if (abs(pf3.pdgId()) == 11) {
              nume += 1;
            }
          }
        } // end pf cand loop
        if (fDebug) cout << ". finished photon" << endl;
        int n = index_of_leading_pf_photon;
        if (n != -1) {
          leading_pf_photon.SetPtEtaPhiE((*pfcands)[n].pt(), (*pfcands)[n].eta(), (*pfcands)[n].phiAtVtx(), (*pfcands)[n].energy());
          TLorentzVector TwoProngObject;
          TwoProngObject = center + photon;
          TLorentzVector TwoProngObject_leadingPfPhoton;
          TwoProngObject_leadingPfPhoton = center + leading_pf_photon;
          double TwoProng_Mass = TwoProngObject_leadingPfPhoton.M();
          if (fabs(TwoProngObject.Eta()) > 2.5) continue;
          // CH pair: 
          // within dr 0.5
          // has at least one pf photon
          // |eta| < 2.5
          //   meets definition of candidate twoprong, fill vectors
          if (pf1.pdgId() > 0) {
            fCand_CHpos_pt.push_back(pfcand1.Pt());
            fCand_CHpos_eta.push_back(pfcand1.Eta());
            fCand_CHpos_phi.push_back(pfcand1.Phi());
            fCand_CHpos_mass.push_back(pfcand1.M());
            fCand_CHneg_pt.push_back(pfcand2.Pt());
            fCand_CHneg_eta.push_back(pfcand2.Eta());
            fCand_CHneg_phi.push_back(pfcand2.Phi());
            fCand_CHneg_mass.push_back(pfcand2.M());
            // impact parameters 
            fCand_CHpos_vz.push_back( pf1.vz() );
            fCand_CHpos_vx.push_back( pf1.vx() );
            fCand_CHpos_vy.push_back( pf1.vy() );
            fCand_CHpos_dz.push_back(          pf1.dzAssociatedPV() );
            fCand_CHpos_dz_PV.push_back(       pf1.dz(((*primaryvertecies)[0]).position()) );
            fCand_CHpos_dz_beamspot.push_back( pf1.dz(beamspot->position()));
            fCand_CHpos_dxy.push_back(          pf1.dxy() );
            fCand_CHpos_dxy_PV.push_back(       pf1.dxy(((*primaryvertecies)[0]).position()) );
            fCand_CHpos_dxy_beamspot.push_back( pf1.dxy(beamspot->position()));
            fCand_CHneg_vz.push_back( pf2.vz() );
            fCand_CHneg_vx.push_back( pf2.vx() );
            fCand_CHneg_vy.push_back( pf2.vy() );
            fCand_CHneg_dz.push_back(          pf2.dzAssociatedPV() );
            fCand_CHneg_dz_PV.push_back(       pf2.dz(((*primaryvertecies)[0]).position()) );
            fCand_CHneg_dz_beamspot.push_back( pf2.dz(beamspot->position()));
            fCand_CHneg_dxy.push_back(          pf2.dxy() );
            fCand_CHneg_dxy_PV.push_back(       pf2.dxy(((*primaryvertecies)[0]).position()) );
            fCand_CHneg_dxy_beamspot.push_back( pf2.dxy(beamspot->position()));
          } else {
            fCand_CHpos_pt.push_back(pfcand2.Pt());
            fCand_CHpos_eta.push_back(pfcand2.Eta());
            fCand_CHpos_phi.push_back(pfcand2.Phi());
            fCand_CHpos_mass.push_back(pfcand2.M());
            fCand_CHneg_pt.push_back(pfcand1.Pt());
            fCand_CHneg_eta.push_back(pfcand1.Eta());
            fCand_CHneg_phi.push_back(pfcand1.Phi());
            fCand_CHneg_mass.push_back(pfcand1.M());
            // impact parameters
            fCand_CHneg_vz.push_back( pf1.vz() );
            fCand_CHneg_vx.push_back( pf1.vx() );
            fCand_CHneg_vy.push_back( pf1.vy() );
            fCand_CHneg_dz.push_back(          pf1.dzAssociatedPV() );
            fCand_CHneg_dz_PV.push_back(       pf1.dz(((*primaryvertecies)[0]).position()) );
            fCand_CHneg_dz_beamspot.push_back( pf1.dz(beamspot->position()));
            fCand_CHneg_dxy.push_back(          pf1.dxy() );
            fCand_CHneg_dxy_PV.push_back(       pf1.dxy(((*primaryvertecies)[0]).position()) );
            fCand_CHneg_dxy_beamspot.push_back( pf1.dxy(beamspot->position()));
            fCand_CHpos_vz.push_back( pf2.vz() );
            fCand_CHpos_vx.push_back( pf2.vx() );
            fCand_CHpos_vy.push_back( pf2.vy() );
            fCand_CHpos_dz.push_back(          pf2.dzAssociatedPV() );
            fCand_CHpos_dz_PV.push_back(       pf2.dz(((*primaryvertecies)[0]).position()) );
            fCand_CHpos_dz_beamspot.push_back( pf2.dz(beamspot->position()));
            fCand_CHpos_dxy.push_back(          pf2.dxy() );
            fCand_CHpos_dxy_PV.push_back(       pf2.dxy(((*primaryvertecies)[0]).position()) );
            fCand_CHpos_dxy_beamspot.push_back( pf2.dxy(beamspot->position()));
          }
          fCand_photon_pt.push_back(photon.Pt());
          fCand_photon_eta.push_back(photon.Eta());
          fCand_photon_phi.push_back(photon.Phi());
          fCand_photon_mass.push_back(photon.M());
          fCand_photon_nGamma.push_back(numgamma);
          fCand_photon_nElectron.push_back(nume);
          fCand_pt.push_back(TwoProngObject.Pt());
          fCand_eta.push_back(TwoProngObject.Eta());
          fCand_phi.push_back(TwoProngObject.Phi());
          fCand_mass.push_back(TwoProngObject.M());
          fCand_Mass.push_back(TwoProng_Mass);
          // Now define isolations
          double chargedIso = 0;
          double neutralIso = 0;
          double egammaIso = 0;
          int chargedIsoCount = 0;
          int neutralIsoCount = 0;
          int egammaIsoCount = 0;
          for (unsigned int m = 0; m < pfcands->size(); m++) {
            const pat::PackedCandidate &pf4 = (*pfcands)[m];
            TLorentzVector pfcand4;
            pfcand4.SetPtEtaPhiE(pf4.pt(), pf4.eta(), pf4.phiAtVtx(), pf4.energy());
            if (pf4.fromPV() <= 1) continue;
            // charged (incl. muons)
            if (abs(pf4.pdgId()) == 13 || abs(pf4.pdgId()) == 211) {
              if ( center.DeltaR(pfcand4) < fCandidatePairIsolationDR && !(m == i || m == j) ) { // don't include one of CH from CH pair
                  chargedIso += pfcand4.Pt();
                  chargedIsoCount++;
                  fCand_isoPF_vz.push_back( pf4.vz() );
                  fCand_isoPF_vx.push_back( pf4.vx() );
                  fCand_isoPF_vy.push_back( pf4.vy() );
                  fCand_isoPF_dz.push_back(          pf4.dzAssociatedPV() );
                  fCand_isoPF_dz_PV.push_back(       pf4.dz(((*primaryvertecies)[0]).position()) );
                  fCand_isoPF_dz_beamspot.push_back( pf4.dz(beamspot->position()));
                  fCand_isoPF_dxy.push_back(          pf4.dxy() );
                  fCand_isoPF_dxy_PV.push_back(       pf4.dxy(((*primaryvertecies)[0]).position()) );
                  fCand_isoPF_dxy_beamspot.push_back( pf4.dxy(beamspot->position()));
               }
            // neutral
            } else if (pf4.pdgId() == 130) {
              if (center.DeltaR(pfcand4) < fCandidatePairIsolationDR) {
                neutralIso += pfcand4.Pt();
                  neutralIsoCount++;
              }
            // e gamma
            } else if (abs(pf4.pdgId()) == 11 || pf4.pdgId() == 22) {
              if ( (center.DeltaR(pfcand4) < fCandidatePairIsolationDR) &&
                   !(fabs(pf4.phiAtVtx() - center.Phi()) < fCandidatePairPhiBox/2.0 && fabs(pf4.eta() - center.Eta()) < fCandidatePairEtaBox/2.0)) {
                egammaIso += pfcand4.Pt();
                  egammaIsoCount++;
              }
            }
          } // end pf cand loop
          double relchargedIso = chargedIso / TwoProngObject.Pt();
          double relneutralIso = neutralIso / TwoProngObject.Pt();
          double relegammaIso = egammaIso / TwoProngObject.Pt();
          if (fDebug) cout << ". finished isolation" << endl;
          fCand_chargedIso.push_back(relchargedIso);
          fCand_neutralIso.push_back(relneutralIso);
          fCand_egammaIso.push_back(relegammaIso);
          fCand_nChargedIsoCone.push_back(chargedIsoCount);
          fCand_nNeutralIsoCone.push_back(neutralIsoCount);
          fCand_nEGammaIsoCone.push_back(egammaIsoCount);
          
          // Selection on Candidates
          bool passCharged = relchargedIso < fCandidatePairChargedIsoCut;
          bool passNeutral = relneutralIso < fCandidatePairNeutralIsoCut;
          bool passEGamma = relegammaIso < fCandidatePairEGammaIsoCut;
          bool passPhotonPt = photon.Pt() > fCandidatePairPhotonPtCut;
          bool tight = passCharged && passNeutral && passEGamma && passPhotonPt;
          bool loose = !tight && passPhotonPt &&
                       relchargedIso < fCandidatePairChargedIsoFakeCut &&
                       relneutralIso < fCandidatePairNeutralIsoFakeCut &&
                       relegammaIso < fCandidatePairEGammaIsoFakeCut;
          fCand_tight.push_back(tight);
          fCand_passChargedIso.push_back(passCharged);
          fCand_passNeutralIso.push_back(passNeutral);
          fCand_passEGammaIso.push_back(passEGamma);
          fCand_passPhotonPt.push_back(passPhotonPt);
          fCand_loose.push_back(loose);
          // Generator Matching
          bool match = false;
          double gen_dR = 99.9;
	        if (fisSignal && fisMC) {
	          for (unsigned int i = 0; i < genparticles->size(); i++) {
	            const reco::GenParticle &genparticle = (*genparticles)[i];
	            if ((genparticle.pdgId() == 221 || genparticle.pdgId() == 331) && genparticle.status() == 2) {
                TLorentzVector genEta;
                genEta.SetPtEtaPhiM(genparticle.pt(), genparticle.eta(), genparticle.phi(), genparticle.mass());
                double match_dR = genEta.DeltaR(TwoProngObject);
                if (match_dR < gen_dR) gen_dR = match_dR;
                if (match_dR < fCandidatePairGenMatchDR)
                  match = true;
              }
            }
          }
          fCand_match.push_back(match);
          fCand_genDR.push_back(gen_dR);
          // Cutflow variables
          if (passCharged) nPassCharged++;
          if (passNeutral) nPassNeutral++;
          if (passEGamma) nPassEGamma++;
          if (passPhotonPt) nPassPhotonPt++;
          if (match) nMatch++;
          if (match && tight) nPassMatch++;
          if (loose) nFake++;
          if (loose) { if (TwoProngObject.Pt() > LeadingFakeTwoProng.Pt()) LeadingFakeTwoProng = TwoProngObject; }
          // Fake rate histograms
          if (!fOmitChargedDecayCode) {
            if (tight) {
              fTwoProngFakeNume_even_pt->Fill(TwoProngObject.Pt(), TwoProng_Mass);
              fTwoProngFakeNume_even_eta->Fill(TwoProngObject.Eta(), TwoProng_Mass);
              fTwoProngFakeNume_even_phi->Fill(TwoProngObject.Phi(), TwoProng_Mass);
              fTwoProngFakeNume_odd_pt->Fill(TwoProngObject.Pt(), TwoProng_Mass);
              fTwoProngFakeNume_odd_eta->Fill(TwoProngObject.Eta(), TwoProng_Mass);
              fTwoProngFakeNume_odd_phi->Fill(TwoProngObject.Phi(), TwoProng_Mass);
            } if (loose) {
              fTwoProngFakeDeno_even_pt->Fill(TwoProngObject.Pt(), TwoProng_Mass);
              fTwoProngFakeDeno_even_eta->Fill(TwoProngObject.Eta(), TwoProng_Mass);
              fTwoProngFakeDeno_even_phi->Fill(TwoProngObject.Phi(), TwoProng_Mass);
              fTwoProngFakeDeno_odd_pt->Fill(TwoProngObject.Pt(), TwoProng_Mass);
              fTwoProngFakeDeno_odd_eta->Fill(TwoProngObject.Eta(), TwoProng_Mass);
              fTwoProngFakeDeno_odd_phi->Fill(TwoProngObject.Phi(), TwoProng_Mass);
            }
            if (fDebug) cout << ". finished fake rate filling" << endl;
          }
        }
      } // end conditionals on CH pair
    }
  } // end making candidates
  fNumPrunedPF = pruned_count;
  fNumTwoProngFake = nFake;

  // More matching, by gen Eta perspective now
  if (fisSignal && fisMC) {
    for (unsigned int i = 0; i < genparticles->size(); i++) {
      const reco::GenParticle &genparticle = (*genparticles)[i];
      TLorentzVector GenParticle;
      GenParticle.SetPtEtaPhiM(genparticle.pt(), genparticle.eta(), genparticle.phi(), genparticle.mass());
      if ((genparticle.pdgId() == 221 || genparticle.pdgId() == 331) && genparticle.status() == 2) {
        double candDR = 99.9;
        double passedCandDR = 99.9;
        double jetDR = 99.9;
        for (unsigned int j = 0; j < fCand_pt.size(); j++) {
          TLorentzVector Candidate;
          Candidate.SetPtEtaPhiM(fCand_pt[j], fCand_eta[j], fCand_phi[j], fCand_mass[j]);
          double dr = Candidate.DeltaR(GenParticle);
          if (dr < candDR) candDR = dr;
        }
        for (unsigned int j = 0; j < fCand_pt.size(); j++) {
          if (!fCand_tight[j]) continue;
          TLorentzVector PassedCandidate;
          PassedCandidate.SetPtEtaPhiM(fCand_pt[j], fCand_eta[j], fCand_phi[j], fCand_mass[j]);
          double dr = PassedCandidate.DeltaR(GenParticle);
          if (dr < passedCandDR) passedCandDR = dr;
        }
        for (unsigned int i = 0; i < ak4jets->size(); i++) {
          const pat::Jet &jet = (*ak4jets)[i];
          TLorentzVector Jet;
          Jet.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), jet.mass());
          double dr = Jet.DeltaR(GenParticle);
          if (dr < jetDR) jetDR = dr;
        }
        fGenEta_candDR.push_back(candDR);
        fGenEta_passedCandDR.push_back(passedCandDR);
        fGenEta_jetDR.push_back(jetDR);
      }
    }
  }

  // Create sorted-by-pt lists
  if (fDebug) cout << ". sorting" << endl;
  vector<unsigned int> sorted_indecies;
  for (unsigned int i = 0; i < fCand_pt.size(); i++) {
    double largestPtSoFar = -1.0;
    unsigned int largestPtSoFarIndex = 999;
    for (unsigned int j = 0; j < fCand_pt.size(); j++) {
      bool skip = false;
      for (unsigned int n = 0; n < sorted_indecies.size(); n++) {
        if (sorted_indecies[n] == j) {
          skip = true;
          break; } }
      if (skip) continue;
      else if (fCand_pt[j] > largestPtSoFar) {
        largestPtSoFar = fCand_pt[j];
        largestPtSoFarIndex = j; } }
    sorted_indecies.push_back(largestPtSoFarIndex);
  }
  if (fDebug) cout << ". finished sorting, filling passed collections" << endl;

  for (unsigned int i = 0; i < sorted_indecies.size(); i++) {
    unsigned int index = sorted_indecies[i];
    if (fCand_loose[index])
    {
      // Candidate is loose and is next leading, fill all loose candidate collections
      fTwoProngLoose_pt.push_back(fCand_pt[index]);
      fTwoProngLoose_eta.push_back(fCand_eta[index]);
      fTwoProngLoose_phi.push_back(fCand_phi[index]);
      TLorentzVector p;
      p.SetPtEtaPhiM(fCand_pt[index], fCand_eta[index], fCand_phi[index], fCand_mass[index]);
      fTwoProngLoose_px.push_back(p.Px());
      fTwoProngLoose_py.push_back(p.Py());
      fTwoProngLoose_pz.push_back(p.Pz());
      fTwoProngLoose_energy.push_back(p.E());
      fTwoProngLoose_mass.push_back(fCand_mass[index]);
      fTwoProngLoose_Mass.push_back(fCand_Mass[index]);
      fTwoProngLoose_CHpos_pt.push_back(fCand_CHpos_pt[index]);
      fTwoProngLoose_CHpos_eta.push_back(fCand_CHpos_eta[index]);
      fTwoProngLoose_CHpos_phi.push_back(fCand_CHpos_phi[index]);
      fTwoProngLoose_CHpos_mass.push_back(fCand_CHpos_mass[index]);
      fTwoProngLoose_CHpos_vz.push_back(fCand_CHpos_vz[index]);
      fTwoProngLoose_CHpos_vx.push_back(fCand_CHpos_vx[index]);
      fTwoProngLoose_CHpos_vy.push_back(fCand_CHpos_vy[index]);
      fTwoProngLoose_CHpos_dz.push_back(fCand_CHpos_dz[index]);
      fTwoProngLoose_CHpos_dz_PV.push_back(fCand_CHpos_dz_PV[index]);
      fTwoProngLoose_CHpos_dz_beamspot.push_back(fCand_CHpos_dz_beamspot[index]);
      fTwoProngLoose_CHpos_dxy.push_back(fCand_CHpos_dxy[index]);
      fTwoProngLoose_CHpos_dxy_PV.push_back(fCand_CHpos_dxy_PV[index]);
      fTwoProngLoose_CHpos_dxy_beamspot.push_back(fCand_CHpos_dxy_beamspot[index]);
      fTwoProngLoose_CHneg_pt.push_back(fCand_CHneg_pt[index]);
      fTwoProngLoose_CHneg_eta.push_back(fCand_CHneg_eta[index]);
      fTwoProngLoose_CHneg_phi.push_back(fCand_CHneg_phi[index]);
      fTwoProngLoose_CHneg_mass.push_back(fCand_CHneg_mass[index]);
      fTwoProngLoose_CHneg_vz.push_back(fCand_CHneg_vz[index]);
      fTwoProngLoose_CHneg_vx.push_back(fCand_CHneg_vx[index]);
      fTwoProngLoose_CHneg_vy.push_back(fCand_CHneg_vy[index]);
      fTwoProngLoose_CHneg_dz.push_back(fCand_CHneg_dz[index]);
      fTwoProngLoose_CHneg_dz_PV.push_back(fCand_CHneg_dz_PV[index]);
      fTwoProngLoose_CHneg_dz_beamspot.push_back(fCand_CHneg_dz_beamspot[index]);
      fTwoProngLoose_CHneg_dxy.push_back(fCand_CHneg_dxy[index]);
      fTwoProngLoose_CHneg_dxy_PV.push_back(fCand_CHneg_dxy_PV[index]);
      fTwoProngLoose_CHneg_dxy_beamspot.push_back(fCand_CHneg_dxy_beamspot[index]);
      fTwoProngLoose_photon_pt.push_back(fCand_photon_pt[index]);
      fTwoProngLoose_photon_eta.push_back(fCand_photon_eta[index]);
      fTwoProngLoose_photon_phi.push_back(fCand_photon_phi[index]);
      fTwoProngLoose_photon_mass.push_back(fCand_photon_mass[index]);
      fTwoProngLoose_photon_nGamma.push_back(fCand_photon_nGamma[index]);
      fTwoProngLoose_photon_nElectron.push_back(fCand_photon_nElectron[index]);
      fTwoProngLoose_chargedIso.push_back(fCand_chargedIso[index]);
      fTwoProngLoose_neutralIso.push_back(fCand_neutralIso[index]);
      fTwoProngLoose_egammaIso.push_back(fCand_egammaIso[index]);
      fTwoProngLoose_match.push_back(fCand_match[index]);
      fTwoProngLoose_genDR.push_back(fCand_genDR[index]);
    }
    if (fCand_tight[index])
    {
      // Candidate is tight and is next leading, fill all tight candidate collections
      fTwoProng_pt.push_back(fCand_pt[index]);
      fTwoProng_eta.push_back(fCand_eta[index]);
      fTwoProng_phi.push_back(fCand_phi[index]);
      TLorentzVector p;
      p.SetPtEtaPhiM(fCand_pt[index], fCand_eta[index], fCand_phi[index], fCand_mass[index]);
      fTwoProng_px.push_back(p.Px());
      fTwoProng_py.push_back(p.Py());
      fTwoProng_pz.push_back(p.Pz());
      fTwoProng_energy.push_back(p.E());
      fTwoProng_mass.push_back(fCand_mass[index]);
      fTwoProng_Mass.push_back(fCand_Mass[index]);
      fTwoProng_CHpos_pt.push_back(fCand_CHpos_pt[index]);
      fTwoProng_CHpos_eta.push_back(fCand_CHpos_eta[index]);
      fTwoProng_CHpos_phi.push_back(fCand_CHpos_phi[index]);
      fTwoProng_CHpos_mass.push_back(fCand_CHpos_mass[index]);
      fTwoProng_CHpos_vz.push_back(fCand_CHpos_vz[index]);
      fTwoProng_CHpos_vx.push_back(fCand_CHpos_vx[index]);
      fTwoProng_CHpos_vy.push_back(fCand_CHpos_vy[index]);
      fTwoProng_CHpos_dz.push_back(fCand_CHpos_dz[index]);
      fTwoProng_CHpos_dz_PV.push_back(fCand_CHpos_dz_PV[index]);
      fTwoProng_CHpos_dz_beamspot.push_back(fCand_CHpos_dz_beamspot[index]);
      fTwoProng_CHpos_dxy.push_back(fCand_CHpos_dxy[index]);
      fTwoProng_CHpos_dxy_PV.push_back(fCand_CHpos_dxy_PV[index]);
      fTwoProng_CHpos_dxy_beamspot.push_back(fCand_CHpos_dxy_beamspot[index]);
      fTwoProng_CHneg_pt.push_back(fCand_CHneg_pt[index]);
      fTwoProng_CHneg_eta.push_back(fCand_CHneg_eta[index]);
      fTwoProng_CHneg_phi.push_back(fCand_CHneg_phi[index]);
      fTwoProng_CHneg_mass.push_back(fCand_CHneg_mass[index]);
      fTwoProng_CHneg_vz.push_back(fCand_CHneg_vz[index]);
      fTwoProng_CHneg_vx.push_back(fCand_CHneg_vx[index]);
      fTwoProng_CHneg_vy.push_back(fCand_CHneg_vy[index]);
      fTwoProng_CHneg_dz.push_back(fCand_CHneg_dz[index]);
      fTwoProng_CHneg_dz_PV.push_back(fCand_CHneg_dz_PV[index]);
      fTwoProng_CHneg_dz_beamspot.push_back(fCand_CHneg_dz_beamspot[index]);
      fTwoProng_CHneg_dxy.push_back(fCand_CHneg_dxy[index]);
      fTwoProng_CHneg_dxy_PV.push_back(fCand_CHneg_dxy_PV[index]);
      fTwoProng_CHneg_dxy_beamspot.push_back(fCand_CHneg_dxy_beamspot[index]);
      fTwoProng_photon_pt.push_back(fCand_photon_pt[index]);
      fTwoProng_photon_eta.push_back(fCand_photon_eta[index]);
      fTwoProng_photon_phi.push_back(fCand_photon_phi[index]);
      fTwoProng_photon_mass.push_back(fCand_photon_mass[index]);
      fTwoProng_photon_nGamma.push_back(fCand_photon_nGamma[index]);
      fTwoProng_photon_nElectron.push_back(fCand_photon_nElectron[index]);
      fTwoProng_chargedIso.push_back(fCand_chargedIso[index]);
      fTwoProng_neutralIso.push_back(fCand_neutralIso[index]);
      fTwoProng_egammaIso.push_back(fCand_egammaIso[index]);
      fTwoProng_match.push_back(fCand_match[index]);
      fTwoProng_genDR.push_back(fCand_genDR[index]);
    }
  }
  if (fDebug) cout << ". finished passed collections" << endl;

  // Jets
  for (unsigned int i = 0; i < ak4jets->size(); i++) {
    const pat::Jet &jet = (*ak4jets)[i];
    fAK4jet_pt.push_back(jet.pt());
    fAK4jet_eta.push_back(jet.eta());
    fAK4jet_phi.push_back(jet.phi());
    fAK4jet_mass.push_back(jet.mass());
    fAK4jet_px.push_back(jet.px());
    fAK4jet_py.push_back(jet.py());
    fAK4jet_pz.push_back(jet.pz());
    fAK4jet_energy.push_back(jet.energy());
  }
  // Photons
  for (unsigned int i = 0; i < miniaod_photons->size(); i++) {
    const pat::Photon &photon = (*miniaod_photons)[i];
    fPhoton_pt.push_back(photon.pt());
    fPhoton_eta.push_back(photon.eta());
    fPhoton_phi.push_back(photon.phi());
    fPhoton_mass.push_back(photon.mass());
    fPhoton_px.push_back(photon.px());
    fPhoton_py.push_back(photon.py());
    fPhoton_pz.push_back(photon.pz());
    fPhoton_energy.push_back(photon.energy());
  }
  // Fill other event wide information
  fNumAK4jets = ak4jets->size();
  fNumPhotons = miniaod_photons->size();
  fNumElectrons = electrons->size();
  fNumMuons = muons->size();
  fNumTwoProng = fCand_pt.size();
  fNumTwoProngPass = fTwoProng_pt.size();
  fNumTwoProngMatched = nMatch;
  fNumTwoProngPassChargedIso = nPassCharged;
  fNumTwoProngPassNeutralIso = nPassNeutral;
  fNumTwoProngPassEGammaIso = nPassEGamma;
  fNumTwoProngPassChargedIso = nPassPhotonPt;
  fHT = 0.0;
  for (unsigned int i = 0; i < ak4jets->size(); i++) {
    const pat::Jet &jet = (*ak4jets)[i];
    if (jet.pt() < 30) continue;
    if (fabs(jet.eta()) > 2.5) continue;
    fHT += jet.pt();
  }
  fMET = (*MET)[0].pt();
  fMET_phi = (*MET)[0].phi();
  
  // Cutflow
  fCutflow_total++;
  if (fCand_pt.size() > 0) fCutflow_oneCand++;
  if (fCand_pt.size() > 1) fCutflow_twoCand++;
  if (fTwoProng_pt.size() > 0)  fCutflow_onePass++;
  if (fTwoProng_pt.size() > 1)  fCutflow_twoPass++;
  if (nPassCharged > 0) fCutflow_passCharged++;
  if (nPassNeutral > 0) fCutflow_passNeutral++;
  if (nPassEGamma > 0) fCutflow_passEGamma++;
  if (nPassPhotonPt > 0) fCutflow_passPhotonPt++;
  if (nMatch > 0) fCutflow_oneMatch++;
  if (nMatch > 1) fCutflow_twoMatch++;
  if (nPassMatch > 0) fCutflow_onePassMatch++;
  if (nPassMatch > 1) fCutflow_twoPassMatch++;
  if (fDebug) cout << ". finish charged decay code part one" << endl;

  if (fDebug) cout << ". starting charged decay code part two" << endl;

  // ***
  // New code block reimplementing photon id  
  // ***
  if (fDebug) cout << ". new block code enter" << endl;

  const CaloSubdetectorTopology* subDetTopologyEB_;
  const CaloSubdetectorTopology* subDetTopologyEE_;
  edm::ESHandle<CaloTopology> caloTopology;
  iSetup.get<CaloTopologyRecord>().get(caloTopology);
  subDetTopologyEB_ = caloTopology->getSubdetectorTopology(DetId::Ecal,EcalBarrel);
  subDetTopologyEE_ = caloTopology->getSubdetectorTopology(DetId::Ecal,EcalEndcap);

  edm::Handle<edm::View<pat::Photon> > ged_photons;
  iEvent.getByToken(gedphotonsToken_,ged_photons);

  if (fDebug) cout << ". got handles" << endl;

  bool isSat = false;
  std::vector<edm::Ptr<pat::Photon>> goodPhotons;
  //std::vector<std::pair<edm::Ptr<pat::Photon>, int> > realAndFakePhotons;

  for (size_t i = 0; i < ged_photons->size(); ++i) {
    const auto pho = ged_photons->ptrAt(i);
    
    isSat = photon_isSaturated(&(*pho), &(*recHitsEB), &(*recHitsEE), &(*subDetTopologyEB_), &(*subDetTopologyEE_));
    bool passID = photon_passHighPtID(&(*pho), rho_, isSat);
    //bool denominatorObject = ExoDiPhotons::passDenominatorCut(&(*pho), rho_, isSat);

    if(passID) {
      goodPhotons.push_back(pho);
      //realAndFakePhotons.push_back(std::pair<edm::Ptr<pat::Photon>, int>(pho, TRUE));
    }
    //if(denominatorObject) {
    //  realAndFakePhotons.push_back(std::pair<edm::Ptr<pat::Photon>, int>(pho, FAKE)); }
  }

  if (fDebug) cout << ". done photon loop" << endl;

  sort(goodPhotons.begin(),goodPhotons.end(),compareCandsByPt);

  if (fDebug) cout << ". done sorting" << endl;

  fNumTightPhotons_v2 = goodPhotons.size();

  InitRecoPhotonInfo(fRecoTightPhotonInfo1_v2);
  InitRecoPhotonInfo(fRecoTightPhotonInfo2_v2);
  InitRecoPhotonInfo(fRecoTightPhotonInfo3_v2);
  if (goodPhotons.size() > 0)
    ExoDiPhotons::FillRecoPhotonInfo(fRecoTightPhotonInfo1_v2,&(*goodPhotons[0]),lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent,iSetup); 
  if (goodPhotons.size() > 1)
    ExoDiPhotons::FillRecoPhotonInfo(fRecoTightPhotonInfo2_v2,&(*goodPhotons[1]),lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent,iSetup); 
  if (goodPhotons.size() > 2)
    ExoDiPhotons::FillRecoPhotonInfo(fRecoTightPhotonInfo3_v2,&(*goodPhotons[2]),lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent,iSetup); 

  if (fDebug) cout << ". new block code exit" << endl;
  // ***
  // end new code block
  // ***

  // Construct Di-Objects
  InitRecoDiObjectInfo(fTwoProngTwoProngInfo);
  // Passed Eta and Passed Eta
  if (fNumTwoProngPass >= 2)
  {
    TLorentzVector Eta1;
    Eta1.SetPtEtaPhiM(fTwoProng_pt[0], fTwoProng_eta[0], fTwoProng_phi[0], fTwoProng_mass[0]);
    TLorentzVector Eta2;
    Eta2.SetPtEtaPhiM(fTwoProng_pt[1], fTwoProng_eta[1], fTwoProng_phi[1], fTwoProng_mass[1]);
    FillRecoDiObjectInfo(fTwoProngTwoProngInfo, Eta1, Eta2);
    fTwoProngTwoProngInfo.dMass = fabs(fTwoProng_Mass[0] - fTwoProng_Mass[1]);
  }

  // generator decay type
  if (fisSignal)
  {
    int decayType = 0;
    for (unsigned int i = 0; i < genparticles->size(); i++) {
      const reco::GenParticle &genparticle = (*genparticles)[i];
      if (genparticle.pdgId() != 9000006 || genparticle.status() != 62) continue;
      if (fDebug) cout << genparticle.pdgId() << ", status=" << genparticle.status() << endl; 
      for (unsigned int j = 0; j < genparticle.numberOfDaughters(); j++) {
        const reco::Candidate * genparticle2 = genparticle.daughter(j);
        if (fDebug) cout << "-> " << genparticle2->pdgId() << ", status=" << genparticle2->status() << endl; 
        int decay1 = 0;
        int decay2 = 0;
        int decay3 = 0;
        for (unsigned int jj = 0; jj < genparticle2->numberOfDaughters(); jj++) {
          const reco::Candidate * genparticle3 = genparticle2->daughter(jj);
          if (fDebug) cout << "  -> " << genparticle3->pdgId() << ", status=" << genparticle3->status() << endl;
          if (jj == 0) decay1 = genparticle3->pdgId();
          if (jj == 1) decay2 = genparticle3->pdgId();
          if (jj == 2) decay3 = genparticle3->pdgId();
        }
        if (isNeutral(decay1, decay2, decay3)) decayType += 0;
        if (isCharged(decay1, decay2, decay3)) decayType += 1;
      }
    }
    fGen_decayType = decayType;
  }

  // Now fill fTree2, it's filled for every event 
  if (!fOmitChargedDecayCode && !fTwoProngFakeRateCalcOnly) fTree2->Fill();
}

void 
ExoDiPhotonAnalyzer::beginJob()
{
}

void 
ExoDiPhotonAnalyzer::endJob()
{
  // fake rate histograms
  if (!fOmitChargedDecayCode) {
  fTwoProngFakeNume_even_pt->Sumw2();
  fTwoProngFakeDeno_even_pt->Sumw2();
  fTwoProngFakeRate_even_pt->Add(fTwoProngFakeNume_even_pt);
  fTwoProngFakeRate_even_pt->Divide(fTwoProngFakeDeno_even_pt);
  fTwoProngFakeNume_even_eta->Sumw2();
  fTwoProngFakeDeno_even_eta->Sumw2();
  fTwoProngFakeRate_even_eta->Add(fTwoProngFakeNume_even_eta);
  fTwoProngFakeRate_even_eta->Divide(fTwoProngFakeDeno_even_eta);
  fTwoProngFakeNume_even_phi->Sumw2();
  fTwoProngFakeDeno_even_phi->Sumw2();
  fTwoProngFakeRate_even_phi->Add(fTwoProngFakeNume_even_phi);
  fTwoProngFakeRate_even_phi->Divide(fTwoProngFakeDeno_even_phi);
  fTwoProngFakeNume_odd_pt->Sumw2();
  fTwoProngFakeDeno_odd_pt->Sumw2();
  fTwoProngFakeRate_odd_pt->Add(fTwoProngFakeNume_odd_pt);
  fTwoProngFakeRate_odd_pt->Divide(fTwoProngFakeDeno_odd_pt);
  fTwoProngFakeNume_odd_eta->Sumw2();
  fTwoProngFakeDeno_odd_eta->Sumw2();
  fTwoProngFakeRate_odd_eta->Add(fTwoProngFakeNume_odd_eta);
  fTwoProngFakeRate_odd_eta->Divide(fTwoProngFakeDeno_odd_eta);
  fTwoProngFakeNume_odd_phi->Sumw2();
  fTwoProngFakeDeno_odd_phi->Sumw2();
  fTwoProngFakeRate_odd_phi->Add(fTwoProngFakeNume_odd_phi);
  fTwoProngFakeRate_odd_phi->Divide(fTwoProngFakeDeno_odd_phi);
  } // done with charged decay code conditional

  // Print Cutflow
  if (!fOmitChargedDecayCode && fchargedDecayCutflow) {
    cout << "Efficiencies for charged decay mode" << endl;
    cout << "Total number of events processed                        : " << fCutflow_total << " " << (fCutflow_total/fCutflow_total)*100 << "%" << endl;
    cout << "Events with at least one candidate Eta                  : " << fCutflow_oneCand << " " << (fCutflow_oneCand/fCutflow_total)*100 << "%" << endl;
    cout << "Events with at least two candidate Etas                 : " << fCutflow_twoCand << " " << (fCutflow_twoCand/fCutflow_total)*100 << "%" << endl;
    cout << "Events with at least one passing Eta                    : " << fCutflow_onePass << " " << (fCutflow_onePass/fCutflow_total)*100 << "%" << endl;
    cout << "Events with at least two passing Etas                   : " << fCutflow_twoPass << " " << (fCutflow_twoPass/fCutflow_total)*100 << "%" << endl;
    cout << "Events with at least one matched candidate              : " << fCutflow_oneMatch << " " << (fCutflow_oneMatch/fCutflow_total)*100 << "%" << endl;
    cout << "Events with at least two matched candidate              : " << fCutflow_twoMatch << " " << (fCutflow_twoMatch/fCutflow_total)*100 << "%" << endl;
    cout << "Events with at least one passing and matched candidate  : " << fCutflow_onePassMatch << " " << (fCutflow_onePassMatch/fCutflow_total)*100 << "%" << endl;
    cout << "Events with at least two passing and matched candidates : " << fCutflow_twoPassMatch << " " << (fCutflow_twoPassMatch/fCutflow_total)*100 << "%" << endl;
    cout << "Cutflow for two prong selection" << endl;
    cout << "Pass relative charged isolation        : " << fCutflow_passCharged << " " << (fCutflow_passCharged/fCutflow_oneCand)*100 << "%" << endl;
    cout << "Pass relative neutral isolation        : " << fCutflow_passNeutral << " " << (fCutflow_passNeutral/fCutflow_oneCand)*100 << "%" << endl;
    cout << "Pass relative egamma isolation         : " << fCutflow_passEGamma << " " << (fCutflow_passEGamma/fCutflow_oneCand)*100 << "%" << endl;
    cout << "Pass pt of combined photon requirement : " << fCutflow_passPhotonPt << " " << (fCutflow_passPhotonPt/fCutflow_oneCand)*100 << "%" << endl;
  }
}

bool ExoDiPhotonAnalyzer::photon_isSaturated(const pat::Photon *photon, const EcalRecHitCollection *recHitsEB, const EcalRecHitCollection *recHitsEE,
		   const CaloSubdetectorTopology* subDetTopologyEB_, const CaloSubdetectorTopology* subDetTopologyEE_) {
    using namespace std;
    
    bool isSat = false;
    DetId seedDetId = ((photon->superCluster())->seed())->seed();
    
    // check EB
    if (seedDetId.subdetId()==EcalBarrel) {
      CaloNavigator<DetId> cursor = CaloNavigator<DetId>(seedDetId,subDetTopologyEB_);
      for (int i = -2; i <= 2; ++i) {
      	for (int j = -2; j <= 2; ++j) {
      	  cursor.home();
      	  cursor.offsetBy(i,j);
      	  EcalRecHitCollection::const_iterator it = recHitsEB->find(*cursor);
      	  if(it != recHitsEB->end()) {
      	    /*cout << "Energy of (" << i << ", " << j << "): " << it-> energy()
	      << ", kSaturated: " << it->checkFlag(EcalRecHit::kSaturated)
	      << ", kDead: " << it->checkFlag(EcalRecHit::kDead)
	      << ", kKilled: " << it->checkFlag(EcalRecHit::kKilled)
	      << endl;*/
      	    if (it->checkFlag(EcalRecHit::kSaturated) && !it->checkFlag(EcalRecHit::kDead) && !it->checkFlag(EcalRecHit::kKilled)) {
      	      isSat = true;
      	    }
      	  }	  
      	}
      }
    }
    // check EE
    else if (seedDetId.subdetId()==EcalEndcap) {
      CaloNavigator<DetId> cursor = CaloNavigator<DetId>(seedDetId,subDetTopologyEE_);
      for (int i = -2; i <= 2; ++i) {
      	for (int j = -2; j <= 2; ++j) {
      	  cursor.home();
      	  cursor.offsetBy(i,j);
      	  EcalRecHitCollection::const_iterator it = recHitsEE->find(*cursor);
      	  if(it != recHitsEE->end()) {
      	    /*cout << "Energy of (" << i << ", " << j << "): " << it->energy()
	      << ", kSaturated: " << it->checkFlag(EcalRecHit::kSaturated)
	      << ", kDead: " << it->checkFlag(EcalRecHit::kDead)
	      << ", kKilled: " << it->checkFlag(EcalRecHit::kKilled)
	      << endl;*/
      	    if (it->checkFlag(EcalRecHit::kSaturated) && !it->checkFlag(EcalRecHit::kDead) && !it->checkFlag(EcalRecHit::kKilled)) {
      	      isSat = true;
      	    }
      	  }
      	}
      }
    }
    return isSat;
  }

bool ExoDiPhotonAnalyzer::photon_passHighPtID(const pat::Photon* photon, double rho, bool isSat) {
    if (
      passHadTowerOverEmCut(photon) &&
      passChargedHadronCut(photon) &&
      passSigmaIetaIetaCut(photon,isSat) &&
      passCorPhoIsoHighPtID(photon,rho) &&
      photon->passElectronVeto()
    ) return true;

    else return false;
  }
// H/E
  bool passHadTowerOverEmCut(const pat::Photon* photon) {
    double hOverE = photon->hadTowOverEm();
    if (hOverE < 0.05) return true;
    else return false;
  }
// CH ISO
  bool passChargedHadronCut(const pat::Photon* photon) {
    double chIsoCut = 5.;
    double chIso = photon->chargedHadronIso();
    if (chIso < chIsoCut) return true;
    else return false;
  }
// SIGMAiETAiETA
  bool passSigmaIetaIetaCut(const pat::Photon* photon, bool isSaturated) {
    double phoEta = fabs(photon->superCluster()->eta());
    double sIeIe = photon->full5x5_sigmaIetaIeta();
    double sIeIeCut = -1.;
    
    if (phoEta < 1.4442 && !isSaturated) sIeIeCut = 0.0105; 
    else if (phoEta < 1.4442 && isSaturated) sIeIeCut = 0.0112;
    else if (1.566 < phoEta && phoEta < 2.5 && !isSaturated) sIeIeCut = 0.0280; 
    else if (1.566 < phoEta && phoEta < 2.5 && isSaturated) sIeIeCut = 0.0300;

    if (sIeIe < sIeIeCut) return true;
    else return false;
  }
// COR ISO
  bool passCorPhoIsoHighPtID(const pat::Photon* photon, double rho) {
    double phoEta = fabs(photon->superCluster()->eta());
    double corPhoIsoCut = -999.9;
    double corPhoIso = corPhoIsoHighPtID(photon,rho);

    if (phoEta < 1.4442) corPhoIsoCut = 2.75;
    if (1.566 < phoEta && phoEta < 2.5) corPhoIsoCut = 2.00;

    if (corPhoIso < corPhoIsoCut) return true;
    else return false;
  }
  double corPhoIsoHighPtID(const pat::Photon* photon, double rho) {
    double phoIso = photon->photonIso();
    return (phoAlphaHighPtID(photon) + phoIso - rho*phoEAHighPtID(photon) - phoKappaHighPtID(photon)*photon->pt());
  }
  double phoAlphaHighPtID(const pat::Photon *photon) {
    double phoEta = fabs(photon->superCluster()->eta());
    if (phoEta < 1.4442) {
      if (phoEta < 0.9) {
	return 2.5;
      }
      else {
	return 2.5;
      }
    } // end EB
    else if (1.566 < phoEta && phoEta < 2.5) {
      if (phoEta < 2.0) {
	return 2.5;
      }
      else if (phoEta < 2.2) {
	return 2.5;
      }
      else {
	return 2.5;
      }
    } // end EE
    else {
      return 99999.99;
    }
  }
  double phoEAHighPtID(const pat::Photon* photon) {
    double phoEta = fabs(photon->superCluster()->eta());
    if (phoEta < 1.4442) {
      if (phoEta < 0.9) {
	return 0.17;
      }
      else {
	return 0.14;
      }
    } // end EB
    else if (1.566 < phoEta && phoEta < 2.5) {
      if (phoEta < 2.0) {
	return 0.11;
      }
      else if (phoEta < 2.2) {
	return 0.14;
      }
      else {
	return 0.22;
      }
    } // end EE
    else {
      return -99999.99;
    }
  }
  double phoKappaHighPtID(const pat::Photon *photon) {
    double phoEta = fabs(photon->superCluster()->eta());
    if (phoEta < 1.4442) {
      if (phoEta < 0.9) {
	return 0.0045;
      }
      else {
	return 0.0045;
      }
    } // end EB
    else if (1.566 < phoEta && phoEta < 2.5) {
      if (phoEta < 2.0) {
	return 0.003;
      }
      else if (phoEta < 2.2) {
	return 0.003;
      }
      else {
	return 0.003;
      }
    } // end EE
    else {
      return -99999.99;
    }
  }

// sorting by pt
bool compareCandsByPt(const edm::Ptr<const reco::Candidate> photon1, const edm::Ptr<const reco::Candidate> photon2) {
    return(photon1->pt()>=photon2->pt());
  }

//define this as a plug-in
DEFINE_FWK_MODULE(ExoDiPhotonAnalyzer);
