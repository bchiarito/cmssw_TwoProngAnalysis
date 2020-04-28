// -*- C++ -*-
//
// Package:    TwoProngAnalyzer
// Class:      TwoProngAnalyzer
// 
/**\class TwoProngAnalyzer TwoProngAnalyzer.cc TwoProngAnalysis/TwoProngAnalyzer/src/TwoProngAnalyzer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Conor Henderson,40 1-B01,+41227671674,
//         Created:  Thu May  6 17:26:16 CEST 2010
// $Id: TwoProngAnalyzer.cc,v 1.32 2013/02/11 15:07:42 charaf Exp $
//
//

// system include files
#include <vector>
#include <algorithm>
#include <cmath>

// ROOT includes
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"

// for new photon code block
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
  //#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

// fileservice
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// PAT objects
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

// for trigger
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

// Gen Event Info
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

// Pileup calculation
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

// TwoProngAnalsysis namespace
#include "TwoProngAnalysis/CommonClasses/interface/RecoPhotonInfo.h"
#include "TwoProngAnalysis/CommonClasses/interface/RecoDiObjectInfo.h"

// TauPreselection namespace
#include "TwoProngAnalysis/ZTauTauFilters/interface/ZtoTauHadPreSelection.h"
#include "TwoProngAnalysis/ZTauTauFilters/interface/ZtoTauHadTruthAlgorithms.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;



//
// class declaration
//


class TwoProngAnalyzer : public edm::EDAnalyzer {
public:
  explicit TwoProngAnalyzer(const edm::ParameterSet&);
  ~TwoProngAnalyzer();

  // constants
  const double PI0_MASS = 0.135;
  const double ETA_MASS = 0.548;

private:
  // EDAnalyzer methods
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // other subroutines
  static bool compareCandsByPt(const edm::Ptr<const reco::Candidate> , const edm::Ptr<const reco::Candidate>);
  bool isNeutral(int one, int two, int three);
  bool isCharged(int one, int two, int three);
  double iso_alpha(double, int);
  double iso_area(double, int);
  double iso_kappa(double, int);

  // photon subroutines
  bool photon_passHighPtID(const pat::Photon*, double , bool );
  bool photon_passHighPtID_loose1(const pat::Photon* , double , bool );
  bool photon_passHighPtID_loose2(const pat::Photon* , double , bool );
  bool photon_passHighPtID_loose3(const pat::Photon* , double , bool );
  bool photon_passHighPtID_loose4(const pat::Photon* , double , bool );
  bool photon_passHighPtID_loose5(const pat::Photon* , double , bool );
  bool photon_passHighPtID_base(const pat::Photon* , double , bool );
  bool photon_passHighPtID_coneHE(const pat::Photon* , double , bool );
  double photon_computeHE(const pat::Photon * photon);
  double photon_computeHE_coneBased(const pat::Photon * photon);
  double photon_computeIsoCh(const pat::Photon * photon);
  double photon_computeIsoGamma(const pat::Photon * photon, double rho);
  double photon_computeSigmaIetaIeta(const pat::Photon * photon);
  bool photon_passHE(const pat::Photon*);
  bool photon_passHE_coneBased(const pat::Photon*);
  bool photon_passIsoCh(const pat::Photon* );
  bool photon_passIsoGamma(const pat::Photon* , double );
  bool photon_passSigmaIetaIeta(const pat::Photon* , bool );
  bool photon_isSaturated(const pat::Photon*, const EcalRecHitCollection *, const EcalRecHitCollection *, const CaloSubdetectorTopology*, const CaloSubdetectorTopology*);
  double photon_scEta(const pat::Photon * photon);
  double photon_computeKappa(const pat::Photon *);
  double photon_computeEA(const pat::Photon* );
  double photon_computeAlpha(const pat::Photon *);

  // ----------member data ---------------------------

  // cmssw config file options
  bool               fMakeTrees;                   // include main Ttrees
  bool               fDebug;                       // turn on stdout for each event
  double             fMcXS;                        // set the value of the mc cross section branch
  double             fMcN;                         // set the value of the mc number generated branch
  bool               fFilterOnPhoton;              // save only events with one or more photons
  bool               fFilterOnTwoProng;            // save only events with one or more twoprongs
  bool               fFilterOnLepton;              // save only events with one or more tight muon (lepton+jets study)
  bool               fFilterForABCDStudy;          // save only events with at least one of: tight/loose twoprong, tight/loose photon
  bool               fincludeDalitzHistos;         // include dalitz plot histograms
  bool               fOldData;                     // runniing on miniAODv2 instead of miniAODv3 (80X vs 94X)
  
  // cmssw config file options, optional branches
  bool               fdontIncludeTwoProngs;        // don't include twoprong object branches
  bool               fincludeLooseTwoProngs;       // include loose twoprong branches, smaller isolation window
  bool               fincludeCandTwoProngs;        // include candidate twoprong branches, no isolation or asymmtery cuts
  bool               fincludeAsymTwoProngs;        // include asym sideband and loose asym sideband twoprong branches
  bool               fincludeMCInfo;               // include MC weight, pthat, and gen PU branches
  bool               fincludeSignalGenParticles;   // include gen particle branches for Phi and omega
  bool               fincludeOldPhotons;           // include old high-pt-id photon object branches, kept because they store more variables
  bool               fincludeBasePhotons;          // include base high-pt-id photon branches, only applies the electron veto
  bool               fincludeConeHEPhotons;        // include cone based high-pt-id photon branches
  bool               fincludeLoosePhotons;         // include loose high-pt-id photon branches

  // cmsssw config file options, Z study related
  bool               fincludeZDecayGenParticles;   // include gen particle branches for Z and its products
  bool               fincludeZTauHadBranches;      // include Z->tau_mu tau_had branches
  bool               fincludeLeptonBranches;       // include branches for lepton+jets control region for bkg estimate
  bool               fincludeZMuMuBranches;        // include Z->mu mu branches
  bool               fusePatTauForZPreBranches;    // use pat::tau instead of tau jet in preselection
  int                fmuonIDtype;                  // muon ID to use, 0, 1, 2 = loose, medium, tight
  int                fmuonISOtype;                 // muon ISO to use, 0, 1, 2, 3, 4, 5 = vloose, loose, medium, tight, vtight, vvtight

  // cmssw config file options, twoprong object
  double             ftwoprong_DR;
  double             ftwoprong_tracksMinPt;
  double             ftwoprong_IsolationDR;
  double             ftwoprong_PhiBox;
  double             ftwoprong_EtaBox;
  double             ftwoprong_PhotonPtCut;
  double             ftwoprong_ChargedIsoCut;
  double             ftwoprong_NeutralIsoCut;
  double             ftwoprong_EGammaIsoCut;
  double             ftwoprong_ChargedIsoFakeCut;
  double             ftwoprong_NeutralIsoFakeCut;
  double             ftwoprong_EGammaIsoFakeCut;
  double             ftwoprong_GenMatchDR;
  double             ftwoprong_AbsMaxEta;
  double             ftwoprong_MinPt;
  double             ftwoprong_TrackAsymmetryCut;
  double             ftwoprong_PhotonAsymmetryCut;
  bool               ftwoprong_OptionalExtraTrack;
  bool               ftwoprong_FlipAsymReq;
  bool               fincludeDalitzVariables;

  // counters for cutflow
  int cutflow_total;
  int cutflow_passFilter;

  // EDM Handles
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  edm::EDGetTokenT<double> rhoToken_; 
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfcandsToken_;
  edm::EDGetTokenT<vector<reco::GenParticle>> genToken_;
  edm::EDGetTokenT<vector<reco::GenJet>> genJetsToken_;
  edm::EDGetTokenT<vector<reco::Vertex>> pvToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamToken_;
  edm::EDGetTokenT<std::vector<pat::Jet>> ak4Token_;
  edm::EDGetTokenT<std::vector<pat::Photon>> photonToken_;
  edm::EDGetTokenT<std::vector<pat::Electron>> electronToken_;
  edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
  edm::EDGetTokenT<std::vector<pat::Tau>> tauToken_;
  edm::EDGetTokenT<std::vector<pat::MET>> metToken_;
  edm::EDGetToken gedphotonsToken_;
  edm::EDGetToken genEventInfoToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupInfoToken_;

  // tools for rechits collection
  std::auto_ptr<noZS::EcalClusterLazyTools> lazyTools_;
  edm::InputTag recHitsEBTag_;
  edm::InputTag recHitsEETag_;
  edm::EDGetTokenT<EcalRecHitCollection> recHitsEBToken;
  edm::EDGetTokenT<EcalRecHitCollection> recHitsEEToken;

  // dalitz histos
  TH2D *fHighvsMid;
  TH2D *fHighvsLow;
  TH2D *fMidvsLow;
  TH2D *fPhotonvsLarger;
  TH2D *fPhotonvsSmaller;
  TH2D *fDoubleStackedDalitz;
  TH2D *fTripleStackedDalitz;
  TH2D *fPhotonvsPositive;
  TH2D *fPhotonvsNegative;
  TH2D *fPositivevsNegative;
  TH1F *fPhotonFraction;
  TH1F *fPositiveFraction;
  TH1F *fNegativeFraction;
  TH1F *fHTverify;
 
  // Main Ntuple Ttree and branches
  TTree *fTree;
  double fzDecayType;
  double fGenZ_pt;
  double fGenZ_eta;
  double fGenZ_phi;
  double fGenZ_mass;
  double fpthat;
  double fHT_gen;
  double fnGenTaus;
  vector<Double_t> fGenTau_pt;
  vector<Double_t> fGenTau_eta;
  vector<Double_t> fGenTau_phi;
  vector<Double_t> fGenTau_mass;
  vector<Double_t> fGenTau_objDR;
  vector<Double_t> fGenTau_candobjDR;

  vector<Double_t> fGenPhi_pt;
  vector<Double_t> fGenPhi_eta;
  vector<Double_t> fGenPhi_phi;
  vector<Double_t> fGenPhi_mass;
  vector<Double_t> fGenPhi_vx;
  vector<Double_t> fGenPhi_vy;
  vector<Double_t> fGenPhi_vz;
  vector<Double_t> fGenPhi_vdiff_beamspot;
  vector<Double_t> fGenPhi_vdiff_PV;
  vector<Double_t> fGenOmega_pt;
  vector<Double_t> fGenOmega_eta;
  vector<Double_t> fGenOmega_phi;
  vector<Double_t> fGenOmega_mass;
  vector<Double_t> fGenOmega_vx;
  vector<Double_t> fGenOmega_vy;
  vector<Double_t> fGenOmega_vz;
  vector<Double_t> fGenOmega_vdiff_beamspot;
  vector<Double_t> fGenOmega_vdiff_PV;
  vector<Double_t> fGenOmega_decayMode;
  vector<Double_t> fGenOmega_neutral_pt;
  vector<Double_t> fGenOmega_neutral_eta;
  vector<Double_t> fGenOmega_neutral_phi;
  vector<Double_t> fGenOmega_neutral_mass;
  vector<Double_t> fGenOmega_positive_pt;
  vector<Double_t> fGenOmega_positive_eta;
  vector<Double_t> fGenOmega_positive_phi;
  vector<Double_t> fGenOmega_positive_mass;
  vector<Double_t> fGenOmega_negative_pt;
  vector<Double_t> fGenOmega_negative_eta;
  vector<Double_t> fGenOmega_negative_phi;
  vector<Double_t> fGenOmega_negative_mass;
  vector<Double_t> fGenOmega_posnegdr;
  vector<Double_t> fGenOmega_objDR;
  vector<Double_t> fGenOmega_candobjDR;
  vector<Double_t> fGenOmega_jetDR;
  vector<Double_t> fGenOmega_pattauDR;
  vector<Double_t> fGenOmega_pattau1DR;
  vector<Double_t> fGenOmega_pattau2DR;
  vector<Double_t> fGenOmega_pattau3DR;

  int fHLT_Photon175;
  int fHLT_Photon200;
  int fHLT_Photon22_Iso;
  int fHLT_Photon165_Iso;
  int fHLT_Photon110_Iso;
  int fHLT_Photon75_VBF;
  int fHLT_DoublePhoton60;
  int fHLT_DoublePhoton70;
  int fHLT_DoublePhoton85;
  int fHLT_Tau120;
  int fHLT_DoubleTau32;
  string fPhotonFoundTrigger;
  int fEventNum;
  int fRunNum;
  int fLumiNum;
  int fNumPVs;
  double fRho;
  int fNumPF;
  double fPV_x;
  double fPV_y;
  double fPV_z;
  double fBeamspot_x;
  double fBeamspot_y;
  double fBeamspot_z;
  double fHT_naive;    // include all jets past pt/eta cut
  double fHT_qcd;      // clean based on EnergyFraction(), try to only include QCD jets
  double fHT;          // clean based on DR to twoprong and photon
  double fST;          // add photon pt
  double fHT_l;        // clean based on DR to twoprong and muon
  double fST_l;        // add muon pt
  double fMET;
  double fMET_phi;
  double fMcW;
  double fMcWProd;
  double ftrueNpu;
  int fobsNpu;

  int fNumAK4jets;
  vector<Double_t> fAK4jet_pt;
  vector<Double_t> fAK4jet_eta;
  vector<Double_t> fAK4jet_phi;
  vector<Double_t> fAK4jet_mass;

  int fNumElectrons;
  vector<Double_t> fElectron_pt;
  vector<Double_t> fElectron_eta;
  vector<Double_t> fElectron_phi;
  vector<Double_t> fElectron_mass;

  int fNumMuons;
  vector<Double_t> fMuon_pt;
  vector<Double_t> fMuon_eta;
  vector<Double_t> fMuon_phi;
  vector<Double_t> fMuon_mass;

  int fNumTaus;
  vector<Double_t> fTau_pt;
  vector<Double_t> fTau_eta;
  vector<Double_t> fTau_phi;
  vector<Double_t> fTau_mass;

  int fNumPhotons; 
  vector<Double_t> fPhoton_pt;
  vector<Double_t> fPhoton_eta;
  vector<Double_t> fPhoton_scEta;
  vector<Double_t> fPhoton_phi;
  vector<Double_t> fPhoton_mass;
  vector<Double_t> fPhoton_isoGamma;
  vector<Double_t> fPhoton_isoCh;
  vector<Double_t> fPhoton_HE;
  vector<Double_t> fPhoton_coneHE;
  vector<Double_t> fPhoton_sigmaIetaIeta;
  vector<Int_t> fPhoton_passVeto;

  int fNumBaseIDPhotons;
  vector<Double_t> fBaseIDPhoton_pt;
  vector<Double_t> fBaseIDPhoton_eta;
  vector<Double_t> fBaseIDPhoton_scEta;
  vector<Double_t> fBaseIDPhoton_phi;
  vector<Double_t> fBaseIDPhoton_mass;
  vector<Double_t> fBaseIDPhoton_isoGamma;
  vector<Double_t> fBaseIDPhoton_isoCh;
  vector<Double_t> fBaseIDPhoton_HE;
  vector<Double_t> fBaseIDPhoton_coneHE;
  vector<Double_t> fBaseIDPhoton_sigmaIetaIeta;

  int fNumLoose1IDPhotons;
  vector<Double_t> fLoose1IDPhoton_pt;
  vector<Double_t> fLoose1IDPhoton_eta;
  vector<Double_t> fLoose1IDPhoton_scEta;
  vector<Double_t> fLoose1IDPhoton_phi;
  vector<Double_t> fLoose1IDPhoton_mass;
  vector<Double_t> fLoose1IDPhoton_isoGamma;
  vector<Double_t> fLoose1IDPhoton_isoCh;
  vector<Double_t> fLoose1IDPhoton_HE;
  vector<Double_t> fLoose1IDPhoton_coneHE;
  vector<Double_t> fLoose1IDPhoton_sigmaIetaIeta;

  int fNumLoose1IDPhotonsEndcap;
  vector<Double_t> fLoose1IDPhotonEndcap_pt;
  vector<Double_t> fLoose1IDPhotonEndcap_eta;
  vector<Double_t> fLoose1IDPhotonEndcap_scEta;
  vector<Double_t> fLoose1IDPhotonEndcap_phi;
  vector<Double_t> fLoose1IDPhotonEndcap_mass;
  vector<Double_t> fLoose1IDPhotonEndcap_isoGamma;
  vector<Double_t> fLoose1IDPhotonEndcap_isoCh;
  vector<Double_t> fLoose1IDPhotonEndcap_HE;
  vector<Double_t> fLoose1IDPhotonEndcap_coneHE;
  vector<Double_t> fLoose1IDPhotonEndcap_sigmaIetaIeta;

  int fNumLoose2IDPhotons;
  vector<Double_t> fLoose2IDPhoton_pt;
  vector<Double_t> fLoose2IDPhoton_eta;
  vector<Double_t> fLoose2IDPhoton_scEta;
  vector<Double_t> fLoose2IDPhoton_phi;
  vector<Double_t> fLoose2IDPhoton_mass;
  vector<Double_t> fLoose2IDPhoton_isoGamma;
  vector<Double_t> fLoose2IDPhoton_isoCh;
  vector<Double_t> fLoose2IDPhoton_HE;
  vector<Double_t> fLoose2IDPhoton_coneHE;
  vector<Double_t> fLoose2IDPhoton_sigmaIetaIeta;

  int fNumLoose2IDPhotonsEndcap;
  vector<Double_t> fLoose2IDPhotonEndcap_pt;
  vector<Double_t> fLoose2IDPhotonEndcap_eta;
  vector<Double_t> fLoose2IDPhotonEndcap_scEta;
  vector<Double_t> fLoose2IDPhotonEndcap_phi;
  vector<Double_t> fLoose2IDPhotonEndcap_mass;
  vector<Double_t> fLoose2IDPhotonEndcap_isoGamma;
  vector<Double_t> fLoose2IDPhotonEndcap_isoCh;
  vector<Double_t> fLoose2IDPhotonEndcap_HE;
  vector<Double_t> fLoose2IDPhotonEndcap_coneHE;
  vector<Double_t> fLoose2IDPhotonEndcap_sigmaIetaIeta;

  int fNumLoose3IDPhotons;
  vector<Double_t> fLoose3IDPhoton_pt;
  vector<Double_t> fLoose3IDPhoton_eta;
  vector<Double_t> fLoose3IDPhoton_scEta;
  vector<Double_t> fLoose3IDPhoton_phi;
  vector<Double_t> fLoose3IDPhoton_mass;
  vector<Double_t> fLoose3IDPhoton_isoGamma;
  vector<Double_t> fLoose3IDPhoton_isoCh;
  vector<Double_t> fLoose3IDPhoton_HE;
  vector<Double_t> fLoose3IDPhoton_coneHE;
  vector<Double_t> fLoose3IDPhoton_sigmaIetaIeta;

  int fNumLoose4IDPhotons;
  vector<Double_t> fLoose4IDPhoton_pt;
  vector<Double_t> fLoose4IDPhoton_eta;
  vector<Double_t> fLoose4IDPhoton_scEta;
  vector<Double_t> fLoose4IDPhoton_phi;
  vector<Double_t> fLoose4IDPhoton_mass;
  vector<Double_t> fLoose4IDPhoton_isoGamma;
  vector<Double_t> fLoose4IDPhoton_isoCh;
  vector<Double_t> fLoose4IDPhoton_HE;
  vector<Double_t> fLoose4IDPhoton_coneHE;
  vector<Double_t> fLoose4IDPhoton_sigmaIetaIeta;

  int fNumLoose5IDPhotons;
  vector<Double_t> fLoose5IDPhoton_pt;
  vector<Double_t> fLoose5IDPhoton_eta;
  vector<Double_t> fLoose5IDPhoton_scEta;
  vector<Double_t> fLoose5IDPhoton_phi;
  vector<Double_t> fLoose5IDPhoton_mass;
  vector<Double_t> fLoose5IDPhoton_isoGamma;
  vector<Double_t> fLoose5IDPhoton_isoCh;
  vector<Double_t> fLoose5IDPhoton_HE;
  vector<Double_t> fLoose5IDPhoton_coneHE;
  vector<Double_t> fLoose5IDPhoton_sigmaIetaIeta;
  vector<Int_t> fLoose5IDPhoton_passVeto;

  int fNumIDPhotons;
  vector<Double_t> fIDPhoton_pt;
  vector<Double_t> fIDPhoton_eta;
  vector<Double_t> fIDPhoton_scEta;
  vector<Double_t> fIDPhoton_phi;
  vector<Double_t> fIDPhoton_mass;
  vector<Double_t> fIDPhoton_isoGamma;
  vector<Double_t> fIDPhoton_isoCh;
  vector<Double_t> fIDPhoton_HE;
  vector<Double_t> fIDPhoton_coneHE;
  vector<Double_t> fIDPhoton_sigmaIetaIeta;

  int fNumIDPhotonsEndcap;
  vector<Double_t> fIDPhotonEndcap_pt;
  vector<Double_t> fIDPhotonEndcap_eta;
  vector<Double_t> fIDPhotonEndcap_scEta;
  vector<Double_t> fIDPhotonEndcap_phi;
  vector<Double_t> fIDPhotonEndcap_mass;
  vector<Double_t> fIDPhotonEndcap_isoGamma;
  vector<Double_t> fIDPhotonEndcap_isoCh;
  vector<Double_t> fIDPhotonEndcap_HE;
  vector<Double_t> fIDPhotonEndcap_coneHE;
  vector<Double_t> fIDPhotonEndcap_sigmaIetaIeta;

  int fNumConeHEIDPhotons;
  vector<Double_t> fConeHEIDPhoton_pt;
  vector<Double_t> fConeHEIDPhoton_eta;
  vector<Double_t> fConeHEIDPhoton_scEta;
  vector<Double_t> fConeHEIDPhoton_phi;
  vector<Double_t> fConeHEIDPhoton_mass;
  vector<Double_t> fConeHEIDPhoton_isoGamma;
  vector<Double_t> fConeHEIDPhoton_isoCh;
  vector<Double_t> fConeHEIDPhoton_HE;
  vector<Double_t> fConeHEIDPhoton_coneHE;
  vector<Double_t> fConeHEIDPhoton_sigmaIetaIeta;

  int fNumConeHEIDPhotonsEndcap;
  vector<Double_t> fConeHEIDPhotonEndcap_pt;
  vector<Double_t> fConeHEIDPhotonEndcap_eta;
  vector<Double_t> fConeHEIDPhotonEndcap_scEta;
  vector<Double_t> fConeHEIDPhotonEndcap_phi;
  vector<Double_t> fConeHEIDPhotonEndcap_mass;
  vector<Double_t> fConeHEIDPhotonEndcap_isoGamma;
  vector<Double_t> fConeHEIDPhotonEndcap_isoCh;
  vector<Double_t> fConeHEIDPhotonEndcap_HE;
  vector<Double_t> fConeHEIDPhotonEndcap_coneHE;
  vector<Double_t> fConeHEIDPhotonEndcap_sigmaIetaIeta;

  TwoProngAnalysis::recoPhotonInfo_t fRecoTightPhotonInfo1;
  TwoProngAnalysis::recoPhotonInfo_t fRecoTightPhotonInfo2;
  TwoProngAnalysis::recoPhotonInfo_t fRecoTightPhotonInfo3;

  int fnTwoProngCands;
  vector<Double_t> fTwoProngCand_pt;
  vector<Double_t> fTwoProngCand_eta;
  vector<Double_t> fTwoProngCand_phi;
  vector<Double_t> fTwoProngCand_mass;
  vector<Double_t> fTwoProngCand_mass_l;
  vector<Double_t> fTwoProngCand_Mass0;
  vector<Double_t> fTwoProngCand_MassPi0;
  vector<Double_t> fTwoProngCand_MassEta;
  vector<Double_t> fTwoProngCand_Mass300;
  vector<Int_t> fTwoProngCand_nExtraTracks;
  vector<Double_t> fTwoProngCand_CHpos_pt;
  vector<Double_t> fTwoProngCand_CHpos_eta;
  vector<Double_t> fTwoProngCand_CHpos_phi;
  vector<Double_t> fTwoProngCand_CHpos_mass;
  vector<Double_t> fTwoProngCand_CHpos_dz;
  vector<Double_t> fTwoProngCand_CHpos_dxy;
  vector<Double_t> fTwoProngCand_CHneg_pt;
  vector<Double_t> fTwoProngCand_CHneg_eta;
  vector<Double_t> fTwoProngCand_CHneg_phi;
  vector<Double_t> fTwoProngCand_CHneg_mass;
  vector<Double_t> fTwoProngCand_CHneg_dz;
  vector<Double_t> fTwoProngCand_CHneg_dxy;
  vector<Double_t> fTwoProngCand_photon_pt;
  vector<Double_t> fTwoProngCand_photon_eta;
  vector<Double_t> fTwoProngCand_photon_phi;
  vector<Double_t> fTwoProngCand_photon_mass;
  vector<Double_t> fTwoProngCand_photon_pt_l;
  vector<Double_t> fTwoProngCand_photon_eta_l;
  vector<Double_t> fTwoProngCand_photon_phi_l;
  vector<Double_t> fTwoProngCand_photon_mass_l;
  vector<Double_t> fTwoProngCand_photon_nGamma;
  vector<Double_t> fTwoProngCand_photon_nElectron;
  vector<Double_t> fTwoProngCand_chargedIso;
  vector<Double_t> fTwoProngCand_neutralIso;
  vector<Double_t> fTwoProngCand_egammaIso;
  vector<Double_t> fTwoProngCand_trackAsym;
  vector<Double_t> fTwoProngCand_photonAsym;
  vector<Int_t> fTwoProngCand_nChargedIsoCone;
  vector<Int_t> fTwoProngCand_nNeutralIsoCone;
  vector<Int_t> fTwoProngCand_nEGammaIsoCone;
  vector<Bool_t> fTwoProngCand_tight;
  vector<Bool_t> fTwoProngCand_loose;
  vector<Bool_t> fTwoProngCand_asym;
  vector<Bool_t> fTwoProngCand_asym_loose;
  vector<Double_t> fTwoProngCand_mPosPho;
  vector<Double_t> fTwoProngCand_mPosPho_l;
  vector<Double_t> fTwoProngCand_mPosPho_pi0;
  vector<Double_t> fTwoProngCand_mNegPho;
  vector<Double_t> fTwoProngCand_mNegPho_l;
  vector<Double_t> fTwoProngCand_mNegPho_pi0;
  vector<Double_t> fTwoProngCand_mPosNeg;
  vector<Double_t> fTwoProngCand_CHpos_p3;
  vector<Double_t> fTwoProngCand_CHneg_p3;
  vector<Double_t> fTwoProngCand_photon_p3;
  vector<Double_t> fTwoProngCand_genOmega_dR;
  vector<Double_t> fTwoProngCand_genTau_dR;

  int fnTwoProngs;
  vector<Double_t> fTwoProng_pt;
  vector<Double_t> fTwoProng_eta;
  vector<Double_t> fTwoProng_phi;
  vector<Double_t> fTwoProng_mass;
  vector<Double_t> fTwoProng_mass_l;
  vector<Double_t> fTwoProng_Mass0;
  vector<Double_t> fTwoProng_MassPi0;
  vector<Double_t> fTwoProng_MassEta;
  vector<Double_t> fTwoProng_Mass300;
  vector<Int_t> fTwoProng_nExtraTracks;
  vector<Double_t> fTwoProng_CHpos_pt;
  vector<Double_t> fTwoProng_CHpos_eta;
  vector<Double_t> fTwoProng_CHpos_phi;
  vector<Double_t> fTwoProng_CHpos_mass;
  vector<Double_t> fTwoProng_CHpos_dz;
  vector<Double_t> fTwoProng_CHpos_dxy;
  vector<Double_t> fTwoProng_CHneg_pt;
  vector<Double_t> fTwoProng_CHneg_eta;
  vector<Double_t> fTwoProng_CHneg_phi;
  vector<Double_t> fTwoProng_CHneg_mass;
  vector<Double_t> fTwoProng_CHneg_dz;
  vector<Double_t> fTwoProng_CHneg_dxy;
  vector<Double_t> fTwoProng_photon_pt;
  vector<Double_t> fTwoProng_photon_eta;
  vector<Double_t> fTwoProng_photon_phi;
  vector<Double_t> fTwoProng_photon_mass;
  vector<Double_t> fTwoProng_photon_pt_l;
  vector<Double_t> fTwoProng_photon_eta_l;
  vector<Double_t> fTwoProng_photon_phi_l;
  vector<Double_t> fTwoProng_photon_mass_l;
  vector<Double_t> fTwoProng_photon_nGamma;
  vector<Double_t> fTwoProng_photon_nElectron;
  vector<Double_t> fTwoProng_chargedIso;
  vector<Double_t> fTwoProng_neutralIso;
  vector<Double_t> fTwoProng_egammaIso;
  vector<Double_t> fTwoProng_trackAsym;
  vector<Double_t> fTwoProng_photonAsym;
  vector<Int_t> fTwoProng_nChargedIsoCone;
  vector<Int_t> fTwoProng_nNeutralIsoCone;
  vector<Int_t> fTwoProng_nEGammaIsoCone;
  vector<Bool_t> fTwoProng_tight;
  vector<Bool_t> fTwoProng_loose;
  vector<Double_t> fTwoProng_mPosPho;
  vector<Double_t> fTwoProng_mPosPho_l;
  vector<Double_t> fTwoProng_mPosPho_pi0;
  vector<Double_t> fTwoProng_mNegPho;
  vector<Double_t> fTwoProng_mNegPho_l;
  vector<Double_t> fTwoProng_mNegPho_pi0;
  vector<Double_t> fTwoProng_mPosNeg;
  vector<Double_t> fTwoProng_CHpos_p3;
  vector<Double_t> fTwoProng_CHneg_p3;
  vector<Double_t> fTwoProng_photon_p3;
  vector<Double_t> fTwoProng_genOmega_dR;
  vector<Double_t> fTwoProng_genTau_dR;

  int fnTwoProngsLoose;
  vector<Double_t> fTwoProngLoose_pt;
  vector<Double_t> fTwoProngLoose_eta;
  vector<Double_t> fTwoProngLoose_phi;
  vector<Double_t> fTwoProngLoose_mass;
  vector<Double_t> fTwoProngLoose_mass_l;
  vector<Double_t> fTwoProngLoose_Mass0;
  vector<Double_t> fTwoProngLoose_MassPi0;
  vector<Double_t> fTwoProngLoose_MassEta;
  vector<Double_t> fTwoProngLoose_Mass300;
  vector<Int_t> fTwoProngLoose_nExtraTracks;
  vector<Double_t> fTwoProngLoose_CHpos_pt;
  vector<Double_t> fTwoProngLoose_CHpos_eta;
  vector<Double_t> fTwoProngLoose_CHpos_phi;
  vector<Double_t> fTwoProngLoose_CHpos_mass;
  vector<Double_t> fTwoProngLoose_CHpos_dz;
  vector<Double_t> fTwoProngLoose_CHpos_dxy;
  vector<Double_t> fTwoProngLoose_CHneg_pt;
  vector<Double_t> fTwoProngLoose_CHneg_eta;
  vector<Double_t> fTwoProngLoose_CHneg_phi;
  vector<Double_t> fTwoProngLoose_CHneg_mass;
  vector<Double_t> fTwoProngLoose_CHneg_dz;
  vector<Double_t> fTwoProngLoose_CHneg_dxy;
  vector<Double_t> fTwoProngLoose_photon_pt;
  vector<Double_t> fTwoProngLoose_photon_eta;
  vector<Double_t> fTwoProngLoose_photon_phi;
  vector<Double_t> fTwoProngLoose_photon_mass;
  vector<Double_t> fTwoProngLoose_photon_pt_l;
  vector<Double_t> fTwoProngLoose_photon_eta_l;
  vector<Double_t> fTwoProngLoose_photon_phi_l;
  vector<Double_t> fTwoProngLoose_photon_mass_l;
  vector<Double_t> fTwoProngLoose_photon_nGamma;
  vector<Double_t> fTwoProngLoose_photon_nElectron;
  vector<Double_t> fTwoProngLoose_chargedIso;
  vector<Double_t> fTwoProngLoose_neutralIso;
  vector<Double_t> fTwoProngLoose_egammaIso;
  vector<Double_t> fTwoProngLoose_trackAsym;
  vector<Double_t> fTwoProngLoose_photonAsym;
  vector<Int_t> fTwoProngLoose_nChargedIsoCone;
  vector<Int_t> fTwoProngLoose_nNeutralIsoCone;
  vector<Int_t> fTwoProngLoose_nEGammaIsoCone;
  vector<Bool_t> fTwoProngLoose_tight;
  vector<Bool_t> fTwoProngLoose_loose;
  vector<Double_t> fTwoProngLoose_mPosPho;
  vector<Double_t> fTwoProngLoose_mPosPho_l;
  vector<Double_t> fTwoProngLoose_mPosPho_pi0;
  vector<Double_t> fTwoProngLoose_mNegPho;
  vector<Double_t> fTwoProngLoose_mNegPho_l;
  vector<Double_t> fTwoProngLoose_mNegPho_pi0;
  vector<Double_t> fTwoProngLoose_mPosNeg;
  vector<Double_t> fTwoProngLoose_CHpos_p3;
  vector<Double_t> fTwoProngLoose_CHneg_p3;
  vector<Double_t> fTwoProngLoose_photon_p3;
  vector<Double_t> fTwoProngLoose_genOmega_dR;
  vector<Double_t> fTwoProngLoose_genTau_dR;

  int fnTwoProngsAsym;
  vector<Double_t> fTwoProngAsym_pt;
  vector<Double_t> fTwoProngAsym_eta;
  vector<Double_t> fTwoProngAsym_phi;
  vector<Double_t> fTwoProngAsym_mass;
  vector<Double_t> fTwoProngAsym_mass_l;
  vector<Double_t> fTwoProngAsym_Mass0;
  vector<Double_t> fTwoProngAsym_MassPi0;
  vector<Double_t> fTwoProngAsym_MassEta;
  vector<Double_t> fTwoProngAsym_Mass300;
  vector<Int_t> fTwoProngAsym_nExtraTracks;
  vector<Double_t> fTwoProngAsym_CHpos_pt;
  vector<Double_t> fTwoProngAsym_CHpos_eta;
  vector<Double_t> fTwoProngAsym_CHpos_phi;
  vector<Double_t> fTwoProngAsym_CHpos_mass;
  vector<Double_t> fTwoProngAsym_CHpos_dz;
  vector<Double_t> fTwoProngAsym_CHpos_dxy;
  vector<Double_t> fTwoProngAsym_CHneg_pt;
  vector<Double_t> fTwoProngAsym_CHneg_eta;
  vector<Double_t> fTwoProngAsym_CHneg_phi;
  vector<Double_t> fTwoProngAsym_CHneg_mass;
  vector<Double_t> fTwoProngAsym_CHneg_dz;
  vector<Double_t> fTwoProngAsym_CHneg_dxy;
  vector<Double_t> fTwoProngAsym_photon_pt;
  vector<Double_t> fTwoProngAsym_photon_eta;
  vector<Double_t> fTwoProngAsym_photon_phi;
  vector<Double_t> fTwoProngAsym_photon_mass;
  vector<Double_t> fTwoProngAsym_photon_pt_l;
  vector<Double_t> fTwoProngAsym_photon_eta_l;
  vector<Double_t> fTwoProngAsym_photon_phi_l;
  vector<Double_t> fTwoProngAsym_photon_mass_l;
  vector<Double_t> fTwoProngAsym_photon_nGamma;
  vector<Double_t> fTwoProngAsym_photon_nElectron;
  vector<Double_t> fTwoProngAsym_chargedIso;
  vector<Double_t> fTwoProngAsym_neutralIso;
  vector<Double_t> fTwoProngAsym_egammaIso;
  vector<Double_t> fTwoProngAsym_trackAsym;
  vector<Double_t> fTwoProngAsym_photonAsym;
  vector<Int_t> fTwoProngAsym_nChargedIsoCone;
  vector<Int_t> fTwoProngAsym_nNeutralIsoCone;
  vector<Int_t> fTwoProngAsym_nEGammaIsoCone;
  vector<Bool_t> fTwoProngAsym_tight;
  vector<Bool_t> fTwoProngAsym_loose;
  vector<Double_t> fTwoProngAsym_mPosPho;
  vector<Double_t> fTwoProngAsym_mPosPho_l;
  vector<Double_t> fTwoProngAsym_mPosPho_pi0;
  vector<Double_t> fTwoProngAsym_mNegPho;
  vector<Double_t> fTwoProngAsym_mNegPho_l;
  vector<Double_t> fTwoProngAsym_mNegPho_pi0;
  vector<Double_t> fTwoProngAsym_mPosNeg;
  vector<Double_t> fTwoProngAsym_CHpos_p3;
  vector<Double_t> fTwoProngAsym_CHneg_p3;
  vector<Double_t> fTwoProngAsym_photon_p3;
  vector<Double_t> fTwoProngAsym_genOmega_dR;
  vector<Double_t> fTwoProngAsym_genTau_dR;

  int fnTwoProngsAsymLoose;
  vector<Double_t> fTwoProngAsymLoose_pt;
  vector<Double_t> fTwoProngAsymLoose_eta;
  vector<Double_t> fTwoProngAsymLoose_phi;
  vector<Double_t> fTwoProngAsymLoose_mass;
  vector<Double_t> fTwoProngAsymLoose_mass_l;
  vector<Double_t> fTwoProngAsymLoose_Mass0;
  vector<Double_t> fTwoProngAsymLoose_MassPi0;
  vector<Double_t> fTwoProngAsymLoose_MassEta;
  vector<Double_t> fTwoProngAsymLoose_Mass300;
  vector<Int_t> fTwoProngAsymLoose_nExtraTracks;
  vector<Double_t> fTwoProngAsymLoose_CHpos_pt;
  vector<Double_t> fTwoProngAsymLoose_CHpos_eta;
  vector<Double_t> fTwoProngAsymLoose_CHpos_phi;
  vector<Double_t> fTwoProngAsymLoose_CHpos_mass;
  vector<Double_t> fTwoProngAsymLoose_CHpos_dz;
  vector<Double_t> fTwoProngAsymLoose_CHpos_dxy;
  vector<Double_t> fTwoProngAsymLoose_CHneg_pt;
  vector<Double_t> fTwoProngAsymLoose_CHneg_eta;
  vector<Double_t> fTwoProngAsymLoose_CHneg_phi;
  vector<Double_t> fTwoProngAsymLoose_CHneg_mass;
  vector<Double_t> fTwoProngAsymLoose_CHneg_dz;
  vector<Double_t> fTwoProngAsymLoose_CHneg_dxy;
  vector<Double_t> fTwoProngAsymLoose_photon_pt;
  vector<Double_t> fTwoProngAsymLoose_photon_eta;
  vector<Double_t> fTwoProngAsymLoose_photon_phi;
  vector<Double_t> fTwoProngAsymLoose_photon_mass;
  vector<Double_t> fTwoProngAsymLoose_photon_pt_l;
  vector<Double_t> fTwoProngAsymLoose_photon_eta_l;
  vector<Double_t> fTwoProngAsymLoose_photon_phi_l;
  vector<Double_t> fTwoProngAsymLoose_photon_mass_l;
  vector<Double_t> fTwoProngAsymLoose_photon_nGamma;
  vector<Double_t> fTwoProngAsymLoose_photon_nElectron;
  vector<Double_t> fTwoProngAsymLoose_chargedIso;
  vector<Double_t> fTwoProngAsymLoose_neutralIso;
  vector<Double_t> fTwoProngAsymLoose_egammaIso;
  vector<Double_t> fTwoProngAsymLoose_trackAsym;
  vector<Double_t> fTwoProngAsymLoose_photonAsym;
  vector<Int_t> fTwoProngAsymLoose_nChargedIsoCone;
  vector<Int_t> fTwoProngAsymLoose_nNeutralIsoCone;
  vector<Int_t> fTwoProngAsymLoose_nEGammaIsoCone;
  vector<Bool_t> fTwoProngAsymLoose_tight;
  vector<Bool_t> fTwoProngAsymLoose_loose;
  vector<Double_t> fTwoProngAsymLoose_mPosPho;
  vector<Double_t> fTwoProngAsymLoose_mPosPho_l;
  vector<Double_t> fTwoProngAsymLoose_mPosPho_pi0;
  vector<Double_t> fTwoProngAsymLoose_mNegPho;
  vector<Double_t> fTwoProngAsymLoose_mNegPho_l;
  vector<Double_t> fTwoProngAsymLoose_mNegPho_pi0;
  vector<Double_t> fTwoProngAsymLoose_mPosNeg;
  vector<Double_t> fTwoProngAsymLoose_CHpos_p3;
  vector<Double_t> fTwoProngAsymLoose_CHneg_p3;
  vector<Double_t> fTwoProngAsymLoose_photon_p3;
  vector<Double_t> fTwoProngAsymLoose_genOmega_dR;
  vector<Double_t> fTwoProngAsymLoose_genTau_dR;

  Int_t fSCat;
  TwoProngAnalysis::recoDiObjectInfo_t fRecoPhiDiTwoProng;
  TwoProngAnalysis::recoDiObjectInfo_t fRecoPhiPhotonTwoProng;
  TwoProngAnalysis::recoDiObjectInfo_t fRecoPhiPhotonTwoProng_1;
  TwoProngAnalysis::recoDiObjectInfo_t fRecoPhiPhotonTwoProng_2;
  TwoProngAnalysis::recoDiObjectInfo_t fRecoPhiPhotonTwoProng_3;
  TwoProngAnalysis::recoDiObjectInfo_t fRecoPhiInclusive;

  Bool_t fpassMuonTrigger;
  Bool_t fpassMuonTriggerTk;
  string fMuonTrigger;
  string fMuonTriggerTk;
  Int_t fnTagMuons;
  Int_t fnProbeTaus;
  Bool_t fpassMuonTauPair;
  Bool_t fpassPreselection;
  Bool_t fpassReducedSelection;
  Bool_t fpassPzeta;
  Bool_t fpassMT;
  Bool_t fpassExtraElectronVeto;
  Bool_t fpassExtraMuonVeto;
  Bool_t fpassDiMuonVeto;
  Bool_t fpassExtraLeptonVeto;
  Bool_t fpassBtagVeto;
  Double_t fMT;
  Double_t fPzeta;
  Double_t fhighestBtag;
  Double_t fTagMuon_pt;
  Double_t fTagMuon_eta;
  Double_t fTagMuon_phi;
  Double_t fTagMuon_mass;
  Double_t fTagMuon_z;
  Double_t fTagMuon_dz;
  Double_t fTagMuon_dB;
  Double_t fTagMuon_dxy;
  Double_t fTagMuon_iso;
  Double_t fProbeTau_pt;
  Double_t fProbeTau_eta;
  Double_t fProbeTau_phi;
  Double_t fProbeTau_mass;
  Double_t fProbeTau_genDR;
  Bool_t fprobePassTauID;
  Double_t fProbeTau_tauID_pt;
  Double_t fProbeTau_tauID_eta;
  Double_t fProbeTau_tauID_phi;
  Double_t fProbeTau_tauID_mass;
  Double_t fProbeTau_tauID_probeDR;
  Bool_t fprobePassTwoProng;
  Double_t fProbeTau_twoprong_pt;
  Double_t fProbeTau_twoprong_eta;
  Double_t fProbeTau_twoprong_phi;
  Double_t fProbeTau_twoprong_mass;
  Double_t fProbeTau_twoprong_probeDR;
  TwoProngAnalysis::recoDiObjectInfo_t fMuonTwoProng;
  TwoProngAnalysis::recoDiObjectInfo_t fMuonTauID;
  TwoProngAnalysis::recoDiObjectInfo_t fMuonProbe;
  TwoProngAnalysis::recoDiObjectInfo_t fMuonTauID_pt;
  TwoProngAnalysis::recoDiObjectInfo_t fMuonTauID_zmass;
  TwoProngAnalysis::recoDiObjectInfo_t fMuonTwoProng_pt;
  TwoProngAnalysis::recoDiObjectInfo_t fMuonTwoProng_zmass;

  Double_t fTagMuon1_pt;
  Double_t fTagMuon1_eta;
  Double_t fTagMuon1_phi;
  Double_t fTagMuon1_mass;
  Double_t fTagMuon1_dz;
  Double_t fTagMuon1_iso;
  Double_t fTagMuon2_pt;
  Double_t fTagMuon2_eta;
  Double_t fTagMuon2_phi;
  Double_t fTagMuon2_mass;
  Double_t fTagMuon2_dz;
  Double_t fTagMuon2_iso;
  TwoProngAnalysis::recoDiObjectInfo_t fMuonMuon;

  Int_t fNTightMuons;
  vector<Double_t> fTightMuon_pt;
  vector<Double_t> fTightMuon_eta;
  vector<Double_t> fTightMuon_phi;
  vector<Double_t> fTightMuon_mass;
  Int_t fMuon_veto;
  Int_t fbtag_veto;
  Int_t fImag_W;
  vector<Double_t> fW_pt;
  vector<Double_t> fW_eta;
  vector<Double_t> fW_phi;
  vector<Double_t> fW_mass;
  vector<Double_t> fW_mT;
  Double_t fmT;
  Int_t fNum_Muons;

};

TwoProngAnalyzer::TwoProngAnalyzer(const edm::ParameterSet& iConfig)
  : 
    fMakeTrees(iConfig.getUntrackedParameter<bool>("makeTrees")),
    fDebug(iConfig.getUntrackedParameter<bool>("debug")),
    fMcXS(iConfig.getUntrackedParameter<double>("mcXS")),
    fMcN(iConfig.getUntrackedParameter<double>("mcN")),
    fFilterOnPhoton(iConfig.getUntrackedParameter<bool>("filterOnPhoton")),
    fFilterOnTwoProng(iConfig.getUntrackedParameter<bool>("filterOnTwoProng")),
    fFilterOnLepton(iConfig.getUntrackedParameter<bool>("filterOnLepton")),
    fFilterForABCDStudy(iConfig.getUntrackedParameter<bool>("filterForABCDStudy")),
    fincludeDalitzHistos(iConfig.getUntrackedParameter<bool>("includeDalitzHistos")),
    fOldData(iConfig.getUntrackedParameter<bool>("oldData")),
    fdontIncludeTwoProngs(iConfig.getUntrackedParameter<bool>("dontIncludeTwoProngs")),
    fincludeLooseTwoProngs(iConfig.getUntrackedParameter<bool>("includeLooseTwoProngs")),
    fincludeCandTwoProngs(iConfig.getUntrackedParameter<bool>("includeCandTwoProngs")),
    fincludeAsymTwoProngs(iConfig.getUntrackedParameter<bool>("includeAsymTwoProngs")),
    fincludeMCInfo(iConfig.getUntrackedParameter<bool>("includeMCInfo")),
    fincludeSignalGenParticles(iConfig.getUntrackedParameter<bool>("includeSignalGenParticles")),
    fincludeOldPhotons(iConfig.getUntrackedParameter<bool>("includeOldPhotons")),
    fincludeBasePhotons(iConfig.getUntrackedParameter<bool>("includeBasePhotons")),
    fincludeConeHEPhotons(iConfig.getUntrackedParameter<bool>("includeConeHEPhotons")),
    fincludeLoosePhotons(iConfig.getUntrackedParameter<bool>("includeLoosePhotons")),
    fincludeZDecayGenParticles(iConfig.getUntrackedParameter<bool>("includeZDecayGenParticles")),
    fincludeZTauHadBranches(iConfig.getUntrackedParameter<bool>("includeZTauHadBranches")),
    fincludeLeptonBranches(iConfig.getUntrackedParameter<bool>("includeLeptonBranches")),
    fincludeZMuMuBranches(iConfig.getUntrackedParameter<bool>("includeZMuMuBranches")),
    fusePatTauForZPreBranches(iConfig.getUntrackedParameter<bool>("usePatTauForZPreBranches")),
    fmuonIDtype(iConfig.getUntrackedParameter<int>("muonIDtype")),
    fmuonISOtype(iConfig.getUntrackedParameter<int>("muonISOtype")),
    ftwoprong_DR(iConfig.getUntrackedParameter<double>("twoprong_chargedHadronPairMinDR")),
    ftwoprong_tracksMinPt(iConfig.getUntrackedParameter<double>("twoprong_chargedHadronMinPt")),
    ftwoprong_IsolationDR(iConfig.getUntrackedParameter<double>("twoprong_isolationConeR")),
    ftwoprong_PhiBox(iConfig.getUntrackedParameter<double>("twoprong_photonPhiBoxSize")),
    ftwoprong_EtaBox(iConfig.getUntrackedParameter<double>("twoprong_photonEtaBoxSize")),
    ftwoprong_PhotonPtCut(iConfig.getUntrackedParameter<double>("twoprong_photonPtCut")),
    ftwoprong_ChargedIsoCut(iConfig.getUntrackedParameter<double>("twoprong_chargedIsoCut")),
    ftwoprong_NeutralIsoCut(iConfig.getUntrackedParameter<double>("twoprong_neutralIsoCut")),
    ftwoprong_EGammaIsoCut(iConfig.getUntrackedParameter<double>("twoprong_egammaIsoCut")),
    ftwoprong_ChargedIsoFakeCut(iConfig.getUntrackedParameter<double>("twoprong_chargedIsoLooseMax")),
    ftwoprong_NeutralIsoFakeCut(iConfig.getUntrackedParameter<double>("twoprong_neutralIsoLooseMax")),
    ftwoprong_EGammaIsoFakeCut(iConfig.getUntrackedParameter<double>("twoprong_egammaIsoLooseMax")),
    ftwoprong_GenMatchDR(iConfig.getUntrackedParameter<double>("twoprong_generatorMatchDR")),
    ftwoprong_AbsMaxEta(iConfig.getUntrackedParameter<double>("twoprong_AbsMaxEta")),
    ftwoprong_MinPt(iConfig.getUntrackedParameter<double>("twoprong_MinPt")),
    ftwoprong_TrackAsymmetryCut(iConfig.getUntrackedParameter<double>("twoprong_MinTrackAsymmetry")),
    ftwoprong_PhotonAsymmetryCut(iConfig.getUntrackedParameter<double>("twoprong_MinPhotonAsymmetry")),
    ftwoprong_OptionalExtraTrack(iConfig.getUntrackedParameter<bool>("twoprong_optionalExtraTrack")),
    ftwoprong_FlipAsymReq(iConfig.getUntrackedParameter<bool>("twoprong_flipAsymReq")),
    fincludeDalitzVariables(iConfig.getUntrackedParameter<bool>("twoprong_includeDalitzVariables")),
    rhoToken_(consumes<double> (iConfig.getParameter<edm::InputTag>("rho")))
{
  if (fDebug) cout << "Entering constructor... " << endl;

  // start cutflow counters
  cutflow_total = 0;
  cutflow_passFilter = 0;

  // MiniAOD event content
  triggerBits_ = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"));
  if (fOldData) {
    triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>(edm::InputTag("selectedPatTrigger"));
  } else {
    triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>(edm::InputTag("slimmedPatTrigger"));
  }
  triggerPrescales_ = consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"));
  pfcandsToken_ = consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates"));
  genToken_ = consumes<vector<reco::GenParticle>>(edm::InputTag("prunedGenParticles"));
  genJetsToken_ = consumes<vector<reco::GenJet>>(edm::InputTag("slimmedGenJets"));
  pvToken_ = consumes<vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));
  beamToken_ = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  ak4Token_ = consumes<std::vector<pat::Jet>>(edm::InputTag("slimmedJets"));
  photonToken_ = consumes<std::vector<pat::Photon>>(edm::InputTag("slimmedPhotons"));
  metToken_ = consumes<std::vector<pat::MET>>(edm::InputTag("slimmedMETs"));
  electronToken_ = consumes<std::vector<pat::Electron>>(edm::InputTag("slimmedElectrons"));
  muonToken_ = consumes<std::vector<pat::Muon>>(edm::InputTag("slimmedMuons"));
  tauToken_ = consumes<std::vector<pat::Tau>>(edm::InputTag("slimmedTaus"));
  gedphotonsToken_ = consumes<edm::View<pat::Photon>>( edm::InputTag("slimmedPhotons") );
  genEventInfoToken_ = mayConsume<GenEventInfoProduct>( edm::InputTag("generator") );
  vtxToken_ = consumes<vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));
  pileupInfoToken_ = consumes<vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"));

  recHitsEBTag_ = iConfig.getUntrackedParameter<edm::InputTag>("RecHitsEBTag",edm::InputTag("reducedEgamma:reducedEBRecHits"));
  recHitsEETag_ = iConfig.getUntrackedParameter<edm::InputTag>("RecHitsEETag",edm::InputTag("reducedEgamma:reducedEERecHits"));
  recHitsEBToken = consumes < EcalRecHitCollection > (recHitsEBTag_);
  recHitsEEToken = consumes < EcalRecHitCollection > (recHitsEETag_);

  // output setup
  edm::Service<TFileService> fs;
  // Branches for charged decay analysis
  if (fMakeTrees) {
  fTree = fs->make<TTree>("fTree","fTree");
  // Generator level
  if (fincludeMCInfo) {
  fTree->Branch("pthat",&fpthat,"pthat/D");
  fTree->Branch("HT_gen",&fHT_gen,"HT_gen/D");
  fTree->Branch("mcW",&fMcW,"mcW/D");
  fTree->Branch("mcWProd",&fMcWProd,"mcWProd/D");
  fTree->Branch("nPU_true",&ftrueNpu,"nPU_true/D");
  fTree->Branch("nPU_obs",&fobsNpu,"nPU_obs/I");
  }
  if (fincludeZDecayGenParticles) {
  fTree->Branch("zDecayType",&fzDecayType,"zDecayType/D");
  fTree->Branch("GenZ_pt",&fGenZ_pt,"GenZ_pt/D");
  fTree->Branch("GenZ_eta",&fGenZ_eta,"GenZ_eta/D");
  fTree->Branch("GenZ_phi",&fGenZ_phi,"GenZ_phi/D");
  fTree->Branch("GenZ_mass",&fGenZ_mass,"GenZ_mass/D");
  fTree->Branch("nGenTaus",&fnGenTaus,"nGenTaus/D");
  fTree->Branch("GenTau_pt",&fGenTau_pt);
  fTree->Branch("GenTau_eta",&fGenTau_eta);
  fTree->Branch("GenTau_phi",&fGenTau_phi);
  fTree->Branch("GenTau_mass",&fGenTau_mass);
  fTree->Branch("GenTau_objDR",&fGenTau_objDR);
  fTree->Branch("GenTau_candobjDR",&fGenTau_candobjDR);
  }
  if (fincludeSignalGenParticles) {
  fTree->Branch("GenPhi_pt",&fGenPhi_pt);
  fTree->Branch("GenPhi_eta",&fGenPhi_eta);
  fTree->Branch("GenPhi_phi",&fGenPhi_phi);
  fTree->Branch("GenPhi_mass",&fGenPhi_mass);
  fTree->Branch("GenPhi_vx",&fGenPhi_vx);
  fTree->Branch("GenPhi_vy",&fGenPhi_vy);
  fTree->Branch("GenPhi_vz",&fGenPhi_vz);
  fTree->Branch("GenPhi_vdiff_beamspot",&fGenPhi_vdiff_beamspot);
  fTree->Branch("GenPhi_vdiff_PV",&fGenPhi_vdiff_PV);
  fTree->Branch("GenOmega_pt",&fGenOmega_pt);
  fTree->Branch("GenOmega_eta",&fGenOmega_eta);
  fTree->Branch("GenOmega_phi",&fGenOmega_phi);
  fTree->Branch("GenOmega_mass",&fGenOmega_mass);
  fTree->Branch("GenOmega_vx",&fGenOmega_vx);
  fTree->Branch("GenOmega_vy",&fGenOmega_vy);
  fTree->Branch("GenOmega_vz",&fGenOmega_vz);
  fTree->Branch("GenOmega_vdiff_beamspot",&fGenOmega_vdiff_beamspot);
  fTree->Branch("GenOmega_vdiff_PV",&fGenOmega_vdiff_PV);
  fTree->Branch("GenOmega_decayMode",&fGenOmega_decayMode);
  fTree->Branch("GenOmega_neutral_pt",&fGenOmega_neutral_pt); 
  fTree->Branch("GenOmega_neutral_eta",&fGenOmega_neutral_eta); 
  fTree->Branch("GenOmega_neutral_phi",&fGenOmega_neutral_phi); 
  fTree->Branch("GenOmega_neutral_mass",&fGenOmega_neutral_mass); 
  fTree->Branch("GenOmega_positive_pt",&fGenOmega_positive_pt); 
  fTree->Branch("GenOmega_positive_eta",&fGenOmega_positive_eta); 
  fTree->Branch("GenOmega_positive_phi",&fGenOmega_positive_phi); 
  fTree->Branch("GenOmega_positive_mass",&fGenOmega_positive_mass); 
  fTree->Branch("GenOmega_negative_pt",&fGenOmega_negative_pt); 
  fTree->Branch("GenOmega_negative_eta",&fGenOmega_negative_eta); 
  fTree->Branch("GenOmega_negative_phi",&fGenOmega_negative_phi); 
  fTree->Branch("GenOmega_negative_mass",&fGenOmega_negative_mass); 
  fTree->Branch("GenOmega_posnegdr",&fGenOmega_posnegdr); 
  fTree->Branch("GenOmega_objDR",&fGenOmega_objDR); 
  fTree->Branch("GenOmega_candobjDR",&fGenOmega_candobjDR); 
  fTree->Branch("GenOmega_pattauDR",&fGenOmega_pattauDR); 
  fTree->Branch("GenOmega_pattau1DR",&fGenOmega_pattau1DR); 
  fTree->Branch("GenOmega_pattau2DR",&fGenOmega_pattau2DR); 
  fTree->Branch("GenOmega_pattau3DR",&fGenOmega_pattau3DR); 
  }
  // Trigger
  fTree->Branch("HLT_Photon175",&fHLT_Photon175,"HLT_Photon175/I");
  fTree->Branch("HLT_Photon200",&fHLT_Photon200,"HLT_Photon200/I");
  fTree->Branch("HLT_Photon22_Iso",&fHLT_Photon22_Iso,"HLT_Photon22_Iso/I");
  fTree->Branch("HLT_Photon165_Iso",&fHLT_Photon165_Iso,"HLT_Photon165_Iso/I");
  fTree->Branch("HLT_Photon110_Iso",&fHLT_Photon110_Iso,"HLT_Photon110_Iso/I");
  fTree->Branch("HLT_Photon75_VBF",&fHLT_Photon75_VBF,"HLT_Photon75_VBF/I");
  fTree->Branch("HLT_DoublePhoton60",&fHLT_DoublePhoton60,"HLT_DoublePhoton60/I");
  fTree->Branch("HLT_DoublePhoton70",&fHLT_DoublePhoton70,"HLT_DoublePhoton70/I");
  fTree->Branch("HLT_DoublePhoton75",&fHLT_DoublePhoton85,"HLT_DoublePhoton85/I");
  fTree->Branch("HLT_Tau120",&fHLT_Tau120,"HLT_Tau120/I");
  fTree->Branch("HLT_DoubleTau32",&fHLT_DoubleTau32,"HLT_DoubleTau32/I");
  fTree->Branch("photonFoundTrigger",&fPhotonFoundTrigger);
  // Event wide
  fTree->Branch("eventNum",&fEventNum,"eventNum/I");
  fTree->Branch("runNum",&fRunNum,"runNum/I");
  fTree->Branch("lumiNum",&fLumiNum,"lumiNum/I");
  fTree->Branch("PV_x",&fPV_x,"PV_x/D");
  fTree->Branch("PV_y",&fPV_y,"PV_y/D");
  fTree->Branch("PV_z",&fPV_z,"PV_z/D");
  fTree->Branch("beamspot_x",&fBeamspot_x,"beamspot_x/D");
  fTree->Branch("beamspot_y",&fBeamspot_y,"beamspot_y/D");
  fTree->Branch("beamspot_z",&fBeamspot_z,"beamspot_z/D");
  fTree->Branch("mcXS",&fMcXS,"mcXS/D");
  fTree->Branch("mcN",&fMcN,"mcN/D");
  fTree->Branch("nPV",&fNumPVs,"nPV/I");
  fTree->Branch("rho",&fRho,"rho/D");
  fTree->Branch("nPF",&fNumPF,"nPF/I");

  fTree->Branch("HT_naive",&fHT_naive,"HT_naive/D");
  fTree->Branch("HT_qcd",&fHT_qcd,"HT_qcd/D");
  fTree->Branch("HT",&fHT,"HT/D");
  fTree->Branch("ST",&fST,"ST/D");
  if (fincludeLeptonBranches) {
    fTree->Branch("HT_l",&fHT_l,"HT_l/D");
    fTree->Branch("ST_l",&fST_l,"ST_l/D");
  }
  fTree->Branch("MET",&fMET,"MET/D");
  fTree->Branch("MET_phi",&fMET_phi,"MET_phi/D");
  // Electrons
  fTree->Branch("nElectrons",&fNumElectrons,"nElectrons/I");
  fTree->Branch("Electron_pt",&fElectron_pt);
  fTree->Branch("Electron_eta",&fElectron_eta);
  fTree->Branch("Electron_phi",&fElectron_phi);
  fTree->Branch("Electron_mass",&fElectron_mass);
  // Muons
  fTree->Branch("nMuons",&fNumMuons,"nMuons/I");
  fTree->Branch("Muon_pt",&fMuon_pt);
  fTree->Branch("Muon_eta",&fMuon_eta);
  fTree->Branch("Muon_phi",&fMuon_phi);
  fTree->Branch("Muon_mass",&fMuon_mass);
  // Taus
  fTree->Branch("nTaus",&fNumTaus,"nTaus/I");
  fTree->Branch("Tau_pt",&fTau_pt);
  fTree->Branch("Tau_eta",&fTau_eta);
  fTree->Branch("Tau_phi",&fTau_phi);
  fTree->Branch("Tau_mass",&fTau_mass);
  // Jets
  fTree->Branch("nJets",&fNumAK4jets,"nJets/I");
  fTree->Branch("Jet_pt",&fAK4jet_pt);
  fTree->Branch("Jet_eta",&fAK4jet_eta);
  fTree->Branch("Jet_phi",&fAK4jet_phi);
  fTree->Branch("Jet_mass",&fAK4jet_mass);
  // Photons
  fTree->Branch("nPhotons",&fNumPhotons,"nPhotons/I");
  fTree->Branch("Photon_pt",&fPhoton_pt);
  fTree->Branch("Photon_eta",&fPhoton_eta);
  fTree->Branch("Photon_scEta",&fPhoton_scEta);
  fTree->Branch("Photon_phi",&fPhoton_phi);
  fTree->Branch("Photon_mass",&fPhoton_mass);
  fTree->Branch("Photon_isoGamma",&fPhoton_isoGamma);
  fTree->Branch("Photon_isoCh",&fPhoton_isoCh);
  fTree->Branch("Photon_HE",&fPhoton_HE);
  fTree->Branch("Photon_coneHE",&fPhoton_coneHE);
  fTree->Branch("Photon_sigmaIetaIeta",&fPhoton_sigmaIetaIeta);
  fTree->Branch("Photon_passVeto",&fPhoton_passVeto);
  // Base High-pt-id Photons, no cuts except electron veto
  if (fincludeBasePhotons) {
  fTree->Branch("nBaseIDPhotons",&fNumBaseIDPhotons,"nBaseIDPhotons/I");
  fTree->Branch("BaseIDPhoton_pt",&fBaseIDPhoton_pt);
  fTree->Branch("BaseIDPhoton_eta",&fBaseIDPhoton_eta);
  fTree->Branch("BaseIDPhoton_scEta",&fBaseIDPhoton_scEta);
  fTree->Branch("BaseIDPhoton_phi",&fBaseIDPhoton_phi);
  fTree->Branch("BaseIDPhoton_mass",&fBaseIDPhoton_mass);
  fTree->Branch("BaseIDPhoton_isoGamma",&fBaseIDPhoton_isoGamma);
  fTree->Branch("BaseIDPhoton_isoCh",&fBaseIDPhoton_isoCh);
  fTree->Branch("BaseIDPhoton_HE",&fBaseIDPhoton_HE);
  fTree->Branch("BaseIDPhoton_coneHE",&fBaseIDPhoton_coneHE);
  fTree->Branch("BaseIDPhoton_sigmaIetaIeta",&fBaseIDPhoton_sigmaIetaIeta);
  }
  // High-pt-id Photons
  fTree->Branch("nIDPhotons",&fNumIDPhotons,"nIDPhotons/I");
  fTree->Branch("IDPhoton_pt",&fIDPhoton_pt);
  fTree->Branch("IDPhoton_eta",&fIDPhoton_eta);
  fTree->Branch("IDPhoton_scEta",&fIDPhoton_scEta);
  fTree->Branch("IDPhoton_phi",&fIDPhoton_phi);
  fTree->Branch("IDPhoton_mass",&fIDPhoton_mass);
  fTree->Branch("IDPhoton_isoGamma",&fIDPhoton_isoGamma);
  fTree->Branch("IDPhoton_isoCh",&fIDPhoton_isoCh);
  fTree->Branch("IDPhoton_HE",&fIDPhoton_HE);
  fTree->Branch("IDPhoton_coneHE",&fIDPhoton_coneHE);
  fTree->Branch("IDPhoton_sigmaIetaIeta",&fIDPhoton_sigmaIetaIeta);
  // High-pt-id Photons, Endcap region
  fTree->Branch("nIDPhotonsEndcap",&fNumIDPhotonsEndcap,"nIDPhotonsEndcap/I");
  fTree->Branch("IDPhotonEndcap_pt",&fIDPhotonEndcap_pt);
  fTree->Branch("IDPhotonEndcap_eta",&fIDPhotonEndcap_eta);
  fTree->Branch("IDPhotonEndcap_scEta",&fIDPhotonEndcap_scEta);
  fTree->Branch("IDPhotonEndcap_phi",&fIDPhotonEndcap_phi);
  fTree->Branch("IDPhotonEndcap_mass",&fIDPhotonEndcap_mass);
  fTree->Branch("IDPhotonEndcap_isoGamma",&fIDPhotonEndcap_isoGamma);
  fTree->Branch("IDPhotonEndcap_isoCh",&fIDPhotonEndcap_isoCh);
  fTree->Branch("IDPhotonEndcap_HE",&fIDPhotonEndcap_HE);
  fTree->Branch("IDPhotonEndcap_coneHE",&fIDPhotonEndcap_coneHE);
  fTree->Branch("IDPhotonEndcap_sigmaIetaIeta",&fIDPhotonEndcap_sigmaIetaIeta);
  // High-pt-id Photons with conebased HE definition
  if (fincludeConeHEPhotons) {
  fTree->Branch("nConeHEIDPhotons",&fNumConeHEIDPhotons,"nConeHEIDPhotons/I");
  fTree->Branch("ConeHEIDPhoton_pt",&fConeHEIDPhoton_pt);
  fTree->Branch("ConeHEIDPhoton_eta",&fConeHEIDPhoton_eta);
  fTree->Branch("ConeHEIDPhoton_scEta",&fConeHEIDPhoton_scEta);
  fTree->Branch("ConeHEIDPhoton_phi",&fConeHEIDPhoton_phi);
  fTree->Branch("ConeHEIDPhoton_mass",&fConeHEIDPhoton_mass);
  fTree->Branch("ConeHEIDPhoton_isoGamma",&fConeHEIDPhoton_isoGamma);
  fTree->Branch("ConeHEIDPhoton_isoCh",&fConeHEIDPhoton_isoCh);
  fTree->Branch("ConeHEIDPhoton_HE",&fConeHEIDPhoton_HE);
  fTree->Branch("ConeHEIDPhoton_coneHE",&fConeHEIDPhoton_coneHE);
  fTree->Branch("ConeHEIDPhoton_sigmaIetaIeta",&fConeHEIDPhoton_sigmaIetaIeta);
  fTree->Branch("nConeHEIDPhotonsEndcap",&fNumConeHEIDPhotonsEndcap,"nConeHEIDPhotonsEndcap/I");
  fTree->Branch("ConeHEIDPhotonEndcap_pt",&fConeHEIDPhotonEndcap_pt);
  fTree->Branch("ConeHEIDPhotonEndcap_eta",&fConeHEIDPhotonEndcap_eta);
  fTree->Branch("ConeHEIDPhotonEndcap_scEta",&fConeHEIDPhotonEndcap_scEta);
  fTree->Branch("ConeHEIDPhotonEndcap_phi",&fConeHEIDPhotonEndcap_phi);
  fTree->Branch("ConeHEIDPhotonEndcap_mass",&fConeHEIDPhotonEndcap_mass);
  fTree->Branch("ConeHEIDPhotonEndcap_isoGamma",&fConeHEIDPhotonEndcap_isoGamma);
  fTree->Branch("ConeHEIDPhotonEndcap_isoCh",&fConeHEIDPhotonEndcap_isoCh);
  fTree->Branch("ConeHEIDPhotonEndcap_HE",&fConeHEIDPhotonEndcap_HE);
  fTree->Branch("ConeHEIDPhotonEndcap_coneHE",&fConeHEIDPhotonEndcap_coneHE);
  fTree->Branch("ConeHEIDPhotonEndcap_sigmaIetaIeta",&fConeHEIDPhotonEndcap_sigmaIetaIeta);
  }
  if (fincludeLoosePhotons) {
  // Loose1 High-pt-id Photons, every cut except iso gamma
  fTree->Branch("nLoose1IDPhotons",&fNumLoose1IDPhotons,"nLoose1IDPhotons/I");
  fTree->Branch("Loose1IDPhoton_pt",&fLoose1IDPhoton_pt);
  fTree->Branch("Loose1IDPhoton_eta",&fLoose1IDPhoton_eta);
  fTree->Branch("Loose1IDPhoton_scEta",&fLoose1IDPhoton_scEta);
  fTree->Branch("Loose1IDPhoton_phi",&fLoose1IDPhoton_phi);
  fTree->Branch("Loose1IDPhoton_mass",&fLoose1IDPhoton_mass);
  fTree->Branch("Loose1IDPhoton_isoGamma",&fLoose1IDPhoton_isoGamma);
  fTree->Branch("Loose1IDPhoton_isoCh",&fLoose1IDPhoton_isoCh);
  fTree->Branch("Loose1IDPhoton_HE",&fLoose1IDPhoton_HE);
  fTree->Branch("Loose1IDPhoton_coneHE",&fLoose1IDPhoton_coneHE);
  fTree->Branch("Loose1IDPhoton_sigmaIetaIeta",&fLoose1IDPhoton_sigmaIetaIeta);
  // Loose1 High-pt-id Photons, every cut except iso gamma, Endcap region
  fTree->Branch("nLoose1IDPhotonsEndcap",&fNumLoose1IDPhotonsEndcap,"nLoose1IDPhotonsEndcap/I");
  fTree->Branch("Loose1IDPhotonEndcap_pt",&fLoose1IDPhotonEndcap_pt);
  fTree->Branch("Loose1IDPhotonEndcap_eta",&fLoose1IDPhotonEndcap_eta);
  fTree->Branch("Loose1IDPhotonEndcap_scEta",&fLoose1IDPhotonEndcap_scEta);
  fTree->Branch("Loose1IDPhotonEndcap_phi",&fLoose1IDPhotonEndcap_phi);
  fTree->Branch("Loose1IDPhotonEndcap_mass",&fLoose1IDPhotonEndcap_mass);
  fTree->Branch("Loose1IDPhotonEndcap_isoGamma",&fLoose1IDPhotonEndcap_isoGamma);
  fTree->Branch("Loose1IDPhotonEndcap_isoCh",&fLoose1IDPhotonEndcap_isoCh);
  fTree->Branch("Loose1IDPhotonEndcap_HE",&fLoose1IDPhotonEndcap_HE);
  fTree->Branch("Loose1IDPhotonEndcap_coneHE",&fLoose1IDPhotonEndcap_coneHE);
  fTree->Branch("Loose1IDPhotonEndcap_sigmaIetaIeta",&fLoose1IDPhotonEndcap_sigmaIetaIeta);
  // Loose2 High-pt-id Photons, every cut except iso ch
  fTree->Branch("nLoose2IDPhotons",&fNumLoose2IDPhotons,"nLoose2IDPhotons/I");
  fTree->Branch("Loose2IDPhoton_pt",&fLoose2IDPhoton_pt);
  fTree->Branch("Loose2IDPhoton_eta",&fLoose2IDPhoton_eta);
  fTree->Branch("Loose2IDPhoton_scEta",&fLoose2IDPhoton_scEta);
  fTree->Branch("Loose2IDPhoton_phi",&fLoose2IDPhoton_phi);
  fTree->Branch("Loose2IDPhoton_mass",&fLoose2IDPhoton_mass);
  fTree->Branch("Loose2IDPhoton_isoGamma",&fLoose2IDPhoton_isoGamma);
  fTree->Branch("Loose2IDPhoton_isoCh",&fLoose2IDPhoton_isoCh);
  fTree->Branch("Loose2IDPhoton_HE",&fLoose2IDPhoton_HE);
  fTree->Branch("Loose2IDPhoton_coneHE",&fLoose2IDPhoton_coneHE);
  fTree->Branch("Loose2IDPhoton_sigmaIetaIeta",&fLoose2IDPhoton_sigmaIetaIeta);
  // Loose2 High-pt-id Photons, every cut except iso ch, Endcap region
  fTree->Branch("nLoose2IDPhotonsEndcap",&fNumLoose2IDPhotonsEndcap,"nLoose2IDPhotonsEndcap/I");
  fTree->Branch("Loose2IDPhotonEndcap_pt",&fLoose2IDPhotonEndcap_pt);
  fTree->Branch("Loose2IDPhotonEndcap_eta",&fLoose2IDPhotonEndcap_eta);
  fTree->Branch("Loose2IDPhotonEndcap_scEta",&fLoose2IDPhotonEndcap_scEta);
  fTree->Branch("Loose2IDPhotonEndcap_phi",&fLoose2IDPhotonEndcap_phi);
  fTree->Branch("Loose2IDPhotonEndcap_mass",&fLoose2IDPhotonEndcap_mass);
  fTree->Branch("Loose2IDPhotonEndcap_isoGamma",&fLoose2IDPhotonEndcap_isoGamma);
  fTree->Branch("Loose2IDPhotonEndcap_isoCh",&fLoose2IDPhotonEndcap_isoCh);
  fTree->Branch("Loose2IDPhotonEndcap_HE",&fLoose2IDPhotonEndcap_HE);
  fTree->Branch("Loose2IDPhotonEndcap_coneHE",&fLoose2IDPhotonEndcap_coneHE);
  fTree->Branch("Loose2IDPhotonEndcap_sigmaIetaIeta",&fLoose2IDPhotonEndcap_sigmaIetaIeta);
  // Loose3 High-pt-id Photons, every cut except H/E
  fTree->Branch("nLoose3IDPhotons",&fNumLoose3IDPhotons,"nLoose3IDPhotons/I");
  fTree->Branch("Loose3IDPhoton_pt",&fLoose3IDPhoton_pt);
  fTree->Branch("Loose3IDPhoton_eta",&fLoose3IDPhoton_eta);
  fTree->Branch("Loose3IDPhoton_scEta",&fLoose3IDPhoton_scEta);
  fTree->Branch("Loose3IDPhoton_phi",&fLoose3IDPhoton_phi);
  fTree->Branch("Loose3IDPhoton_mass",&fLoose3IDPhoton_mass);
  fTree->Branch("Loose3IDPhoton_isoGamma",&fLoose3IDPhoton_isoGamma);
  fTree->Branch("Loose3IDPhoton_isoCh",&fLoose3IDPhoton_isoCh);
  fTree->Branch("Loose3IDPhoton_HE",&fLoose3IDPhoton_HE);
  fTree->Branch("Loose3IDPhoton_coneHE",&fLoose3IDPhoton_coneHE);
  fTree->Branch("Loose3IDPhoton_sigmaIetaIeta",&fLoose3IDPhoton_sigmaIetaIeta);
  // Loose4 High-pt-id Photons, every cut except sigma_ietaieta
  fTree->Branch("nLoose4IDPhotons",&fNumLoose4IDPhotons,"nLoose4IDPhotons/I");
  fTree->Branch("Loose4IDPhoton_pt",&fLoose4IDPhoton_pt);
  fTree->Branch("Loose4IDPhoton_eta",&fLoose4IDPhoton_eta);
  fTree->Branch("Loose4IDPhoton_scEta",&fLoose4IDPhoton_scEta);
  fTree->Branch("Loose4IDPhoton_phi",&fLoose4IDPhoton_phi);
  fTree->Branch("Loose4IDPhoton_mass",&fLoose4IDPhoton_mass);
  fTree->Branch("Loose4IDPhoton_isoGamma",&fLoose4IDPhoton_isoGamma);
  fTree->Branch("Loose4IDPhoton_isoCh",&fLoose4IDPhoton_isoCh);
  fTree->Branch("Loose4IDPhoton_HE",&fLoose4IDPhoton_HE);
  fTree->Branch("Loose4IDPhoton_coneHE",&fLoose4IDPhoton_coneHE);
  fTree->Branch("Loose4IDPhoton_sigmaIetaIeta",&fLoose4IDPhoton_sigmaIetaIeta);
  // Loose5 High-pt-id Photons, every cut except electron veto
  fTree->Branch("nLoose5IDPhotons",&fNumLoose5IDPhotons,"nLoose5IDPhotons/I");
  fTree->Branch("Loose5IDPhoton_pt",&fLoose5IDPhoton_pt);
  fTree->Branch("Loose5IDPhoton_eta",&fLoose5IDPhoton_eta);
  fTree->Branch("Loose5IDPhoton_scEta",&fLoose5IDPhoton_scEta);
  fTree->Branch("Loose5IDPhoton_phi",&fLoose5IDPhoton_phi);
  fTree->Branch("Loose5IDPhoton_mass",&fLoose5IDPhoton_mass);
  fTree->Branch("Loose5IDPhoton_isoGamma",&fLoose5IDPhoton_isoGamma);
  fTree->Branch("Loose5IDPhoton_isoCh",&fLoose5IDPhoton_isoCh);
  fTree->Branch("Loose5IDPhoton_HE",&fLoose5IDPhoton_HE);
  fTree->Branch("Loose5IDPhoton_coneHE",&fLoose5IDPhoton_coneHE);
  fTree->Branch("Loose5IDPhoton_sigmaIetaIeta",&fLoose5IDPhoton_sigmaIetaIeta);
  fTree->Branch("Loose5IDPhoton_passVeto",&fLoose5IDPhoton_passVeto);
  }
  if(fincludeOldPhotons) {
  fTree->Branch("Photon1",&fRecoTightPhotonInfo1,TwoProngAnalysis::recoPhotonBranchDefString.c_str());
  fTree->Branch("Photon2",&fRecoTightPhotonInfo2,TwoProngAnalysis::recoPhotonBranchDefString.c_str());
  fTree->Branch("Photon3",&fRecoTightPhotonInfo3,TwoProngAnalysis::recoPhotonBranchDefString.c_str());
  }
  if (!fdontIncludeTwoProngs) {
  fTree->Branch("nTwoProngs",&fnTwoProngs,"nTwoProngs/I");
  // kinematics
  fTree->Branch("TwoProng_pt",&fTwoProng_pt);
  fTree->Branch("TwoProng_eta",&fTwoProng_eta);
  fTree->Branch("TwoProng_phi",&fTwoProng_phi);
  // mass
  fTree->Branch("TwoProng_mass",&fTwoProng_mass);
  fTree->Branch("TwoProng_mass_l",&fTwoProng_mass_l);
  fTree->Branch("TwoProng_Mass0",&fTwoProng_Mass0);
  fTree->Branch("TwoProng_MassPi0",&fTwoProng_MassPi0);
  fTree->Branch("TwoProng_MassEta",&fTwoProng_MassEta);
  fTree->Branch("TwoProng_Mass300",&fTwoProng_Mass300);
  // extra optional track 
  fTree->Branch("TwoProng_nExtraTracks",&fTwoProng_nExtraTracks);
  // CHpos constituent
  fTree->Branch("TwoProng_CHpos_pt",&fTwoProng_CHpos_pt);
  fTree->Branch("TwoProng_CHpos_eta",&fTwoProng_CHpos_eta);
  fTree->Branch("TwoProng_CHpos_phi",&fTwoProng_CHpos_phi);
  fTree->Branch("TwoProng_CHpos_mass",&fTwoProng_CHpos_mass);
  fTree->Branch("TwoProng_CHpos_dz",&fTwoProng_CHpos_dz);
  fTree->Branch("TwoProng_CHpos_dxy",&fTwoProng_CHpos_dxy);
  // CHneg constituent
  fTree->Branch("TwoProng_CHneg_pt",&fTwoProng_CHneg_pt);
  fTree->Branch("TwoProng_CHneg_eta",&fTwoProng_CHneg_eta);
  fTree->Branch("TwoProng_CHneg_phi",&fTwoProng_CHneg_phi);
  fTree->Branch("TwoProng_CHneg_mass",&fTwoProng_CHneg_mass);
  fTree->Branch("TwoProng_CHneg_dz",&fTwoProng_CHneg_dz);
  fTree->Branch("TwoProng_CHneg_dxy",&fTwoProng_CHneg_dxy);
  // photon constituent
  fTree->Branch("TwoProng_photon_pt",&fTwoProng_photon_pt);
  fTree->Branch("TwoProng_photon_eta",&fTwoProng_photon_eta);
  fTree->Branch("TwoProng_photon_phi",&fTwoProng_photon_phi);
  fTree->Branch("TwoProng_photon_mass",&fTwoProng_photon_mass);
  fTree->Branch("TwoProng_photon_pt_l",&fTwoProng_photon_pt_l);
  fTree->Branch("TwoProng_photon_eta_l",&fTwoProng_photon_eta_l);
  fTree->Branch("TwoProng_photon_phi_l",&fTwoProng_photon_phi_l);
  fTree->Branch("TwoProng_photon_mass_l",&fTwoProng_photon_mass_l);
  fTree->Branch("TwoProng_photon_nGamma",&fTwoProng_photon_nGamma);
  fTree->Branch("TwoProng_photon_nElectron",&fTwoProng_photon_nElectron);
  // cut variables
  fTree->Branch("TwoProng_chargedIso",&fTwoProng_chargedIso);
  fTree->Branch("TwoProng_neutralIso",&fTwoProng_neutralIso);
  fTree->Branch("TwoProng_egammaIso",&fTwoProng_egammaIso);
  fTree->Branch("TwoProng_trackAsym",&fTwoProng_trackAsym);
  fTree->Branch("TwoProng_photonAsym",&fTwoProng_photonAsym);
  fTree->Branch("TwoProng_nChargedIsoCone",&fTwoProng_nChargedIsoCone);
  fTree->Branch("TwoProng_nNeutralIsoCone",&fTwoProng_nNeutralIsoCone);
  fTree->Branch("TwoProng_nEGammaIsoCone",&fTwoProng_nEGammaIsoCone);
  // matching
  fTree->Branch("TwoProng_genOmega_dR",&fTwoProng_genOmega_dR);
  fTree->Branch("TwoProng_genTau_dR",&fTwoProng_genTau_dR);
  // dalitz
  if (fincludeDalitzVariables) {
  fTree->Branch("TwoProng_mPosPho",&fTwoProng_mPosPho);
  fTree->Branch("TwoProng_mPosPho_l",&fTwoProng_mPosPho_l);
  fTree->Branch("TwoProng_mPosPho_pi0",&fTwoProng_mPosPho_pi0);
  fTree->Branch("TwoProng_mNegPho",&fTwoProng_mNegPho);
  fTree->Branch("TwoProng_mNegPho_l",&fTwoProng_mNegPho_l);
  fTree->Branch("TwoProng_mNegPho_pi0",&fTwoProng_mNegPho_pi0);
  fTree->Branch("TwoProng_mPosNeg",&fTwoProng_mPosNeg);
  fTree->Branch("TwoProng_CHpos_p3",&fTwoProng_CHpos_p3);
  fTree->Branch("TwoProng_CHneg_p3",&fTwoProng_CHneg_p3);
  fTree->Branch("TwoProng_photon_p3",&fTwoProng_photon_p3);
  }
  // TwoProngs
  if(fincludeCandTwoProngs) {
    // Two Prong Candidate information, no isolation or asymmtery cuts
  fTree->Branch("nTwoProngCands",&fnTwoProngCands,"nTwoProngCands/I");
  // kinematics
  fTree->Branch("TwoProngCand_pt",&fTwoProngCand_pt);
  fTree->Branch("TwoProngCand_eta",&fTwoProngCand_eta);
  fTree->Branch("TwoProngCand_phi",&fTwoProngCand_phi);
  // mass
  fTree->Branch("TwoProngCand_mass",&fTwoProngCand_mass);
  fTree->Branch("TwoProngCand_mass_l",&fTwoProngCand_mass_l);
  fTree->Branch("TwoProngCand_Mass0",&fTwoProngCand_Mass0);
  fTree->Branch("TwoProngCand_MassPi0",&fTwoProngCand_MassPi0);
  fTree->Branch("TwoProngCand_MassEta",&fTwoProngCand_MassEta);
  fTree->Branch("TwoProngCand_Mass300",&fTwoProngCand_Mass300);
  // extra optional track 
  fTree->Branch("TwoProngCand_nExtraTracks",&fTwoProngCand_nExtraTracks);
  // CHpos constituent
  fTree->Branch("TwoProngCand_CHpos_pt",&fTwoProngCand_CHpos_pt);
  fTree->Branch("TwoProngCand_CHpos_eta",&fTwoProngCand_CHpos_eta);
  fTree->Branch("TwoProngCand_CHpos_phi",&fTwoProngCand_CHpos_phi);
  fTree->Branch("TwoProngCand_CHpos_mass",&fTwoProngCand_CHpos_mass);
  fTree->Branch("TwoProngCand_CHpos_dz",&fTwoProngCand_CHpos_dz);
  fTree->Branch("TwoProngCand_CHpos_dxy",&fTwoProngCand_CHpos_dxy);
  // CHneg constituent
  fTree->Branch("TwoProngCand_CHneg_pt",&fTwoProngCand_CHneg_pt);
  fTree->Branch("TwoProngCand_CHneg_eta",&fTwoProngCand_CHneg_eta);
  fTree->Branch("TwoProngCand_CHneg_phi",&fTwoProngCand_CHneg_phi);
  fTree->Branch("TwoProngCand_CHneg_mass",&fTwoProngCand_CHneg_mass);
  fTree->Branch("TwoProngCand_CHneg_dz",&fTwoProngCand_CHneg_dz);
  fTree->Branch("TwoProngCand_CHneg_dxy",&fTwoProngCand_CHneg_dxy);
  // photon constituent
  fTree->Branch("TwoProngCand_photon_pt",&fTwoProngCand_photon_pt);
  fTree->Branch("TwoProngCand_photon_eta",&fTwoProngCand_photon_eta);
  fTree->Branch("TwoProngCand_photon_phi",&fTwoProngCand_photon_phi);
  fTree->Branch("TwoProngCand_photon_mass",&fTwoProngCand_photon_mass);
  fTree->Branch("TwoProngCand_photon_pt_l",&fTwoProngCand_photon_pt_l);
  fTree->Branch("TwoProngCand_photon_eta_l",&fTwoProngCand_photon_eta_l);
  fTree->Branch("TwoProngCand_photon_phi_l",&fTwoProngCand_photon_phi_l);
  fTree->Branch("TwoProngCand_photon_mass_l",&fTwoProngCand_photon_mass_l);
  fTree->Branch("TwoProngCand_photon_nGamma",&fTwoProngCand_photon_nGamma);
  fTree->Branch("TwoProngCand_photon_nElectron",&fTwoProngCand_photon_nElectron);
  // cut variables
  fTree->Branch("TwoProngCand_chargedIso",&fTwoProngCand_chargedIso);
  fTree->Branch("TwoProngCand_neutralIso",&fTwoProngCand_neutralIso);
  fTree->Branch("TwoProngCand_egammaIso",&fTwoProngCand_egammaIso);
  fTree->Branch("TwoProngCand_trackAsym",&fTwoProngCand_trackAsym);
  fTree->Branch("TwoProngCand_photonAsym",&fTwoProngCand_photonAsym);
  fTree->Branch("TwoProngCand_nChargedIsoCone",&fTwoProngCand_nChargedIsoCone);
  fTree->Branch("TwoProngCand_nNeutralIsoCone",&fTwoProngCand_nNeutralIsoCone);
  fTree->Branch("TwoProngCand_nEGammaIsoCone",&fTwoProngCand_nEGammaIsoCone);
  fTree->Branch("TwoProngCand_tight",&fTwoProngCand_tight);
  fTree->Branch("TwoProngCand_loose",&fTwoProngCand_loose);
  fTree->Branch("TwoProngCand_asym",&fTwoProngCand_asym);
  fTree->Branch("TwoProngCand_asym_loose",&fTwoProngCand_asym_loose);
  // matching
  fTree->Branch("TwoProngCand_genOmega_dR",&fTwoProngCand_genOmega_dR);
  fTree->Branch("TwoProngCand_genTau_dR",&fTwoProngCand_genTau_dR);
  // dalitz
  if (fincludeDalitzVariables) {
  fTree->Branch("TwoProngCand_mPosPho",&fTwoProngCand_mPosPho);
  fTree->Branch("TwoProngCand_mPosPho_l",&fTwoProngCand_mPosPho_l);
  fTree->Branch("TwoProngCand_mPosPho_pi0",&fTwoProngCand_mPosPho_pi0);
  fTree->Branch("TwoProngCand_mNegPho",&fTwoProngCand_mNegPho);
  fTree->Branch("TwoProngCand_mNegPho_l",&fTwoProngCand_mNegPho_l);
  fTree->Branch("TwoProngCand_mNegPho_pi0",&fTwoProngCand_mNegPho_pi0);
  fTree->Branch("TwoProngCand_mPosNeg",&fTwoProngCand_mPosNeg);
  fTree->Branch("TwoProngCand_CHpos_p3",&fTwoProngCand_CHpos_p3);
  fTree->Branch("TwoProngCand_CHneg_p3",&fTwoProngCand_CHneg_p3);
  fTree->Branch("TwoProngCand_photon_p3",&fTwoProngCand_photon_p3);
  }
  }
  if(fincludeLooseTwoProngs) {
    // Loose TwoProng information
  fTree->Branch("nTwoProngsLoose",&fnTwoProngsLoose,"nTwoProngsLoose/I");
  // kinematics
  fTree->Branch("TwoProngLoose_pt",&fTwoProngLoose_pt);
  fTree->Branch("TwoProngLoose_eta",&fTwoProngLoose_eta);
  fTree->Branch("TwoProngLoose_phi",&fTwoProngLoose_phi);
  // mass
  fTree->Branch("TwoProngLoose_mass",&fTwoProngLoose_mass);
  fTree->Branch("TwoProngLoose_mass_l",&fTwoProngLoose_mass_l);
  fTree->Branch("TwoProngLoose_Mass0",&fTwoProngLoose_Mass0);
  fTree->Branch("TwoProngLoose_MassPi0",&fTwoProngLoose_MassPi0);
  fTree->Branch("TwoProngLoose_MassEta",&fTwoProngLoose_MassEta);
  fTree->Branch("TwoProngLoose_Mass300",&fTwoProngLoose_Mass300);
  // extra optional track 
  fTree->Branch("TwoProngLoose_nExtraTracks",&fTwoProngLoose_nExtraTracks);
  // CHpos constituent
  fTree->Branch("TwoProngLoose_CHpos_pt",&fTwoProngLoose_CHpos_pt);
  fTree->Branch("TwoProngLoose_CHpos_eta",&fTwoProngLoose_CHpos_eta);
  fTree->Branch("TwoProngLoose_CHpos_phi",&fTwoProngLoose_CHpos_phi);
  fTree->Branch("TwoProngLoose_CHpos_mass",&fTwoProngLoose_CHpos_mass);
  fTree->Branch("TwoProngLoose_CHpos_dz",&fTwoProngLoose_CHpos_dz);
  fTree->Branch("TwoProngLoose_CHpos_dxy",&fTwoProngLoose_CHpos_dxy);
  // CHneg constituent
  fTree->Branch("TwoProngLoose_CHneg_pt",&fTwoProngLoose_CHneg_pt);
  fTree->Branch("TwoProngLoose_CHneg_eta",&fTwoProngLoose_CHneg_eta);
  fTree->Branch("TwoProngLoose_CHneg_phi",&fTwoProngLoose_CHneg_phi);
  fTree->Branch("TwoProngLoose_CHneg_mass",&fTwoProngLoose_CHneg_mass);
  fTree->Branch("TwoProngLoose_CHneg_dz",&fTwoProngLoose_CHneg_dz);
  fTree->Branch("TwoProngLoose_CHneg_dxy",&fTwoProngLoose_CHneg_dxy);
  // photon constituent
  fTree->Branch("TwoProngLoose_photon_pt",&fTwoProngLoose_photon_pt);
  fTree->Branch("TwoProngLoose_photon_eta",&fTwoProngLoose_photon_eta);
  fTree->Branch("TwoProngLoose_photon_phi",&fTwoProngLoose_photon_phi);
  fTree->Branch("TwoProngLoose_photon_mass",&fTwoProngLoose_photon_mass);
  fTree->Branch("TwoProngLoose_photon_pt_l",&fTwoProngLoose_photon_pt_l);
  fTree->Branch("TwoProngLoose_photon_eta_l",&fTwoProngLoose_photon_eta_l);
  fTree->Branch("TwoProngLoose_photon_phi_l",&fTwoProngLoose_photon_phi_l);
  fTree->Branch("TwoProngLoose_photon_mass_l",&fTwoProngLoose_photon_mass_l);
  fTree->Branch("TwoProngLoose_photon_nGamma",&fTwoProngLoose_photon_nGamma);
  fTree->Branch("TwoProngLoose_photon_nElectron",&fTwoProngLoose_photon_nElectron);
  // cut variables
  fTree->Branch("TwoProngLoose_chargedIso",&fTwoProngLoose_chargedIso);
  fTree->Branch("TwoProngLoose_neutralIso",&fTwoProngLoose_neutralIso);
  fTree->Branch("TwoProngLoose_egammaIso",&fTwoProngLoose_egammaIso);
  fTree->Branch("TwoProngLoose_trackAsym",&fTwoProngLoose_trackAsym);
  fTree->Branch("TwoProngLoose_photonAsym",&fTwoProngLoose_photonAsym);
  fTree->Branch("TwoProngLoose_nChargedIsoCone",&fTwoProngLoose_nChargedIsoCone);
  fTree->Branch("TwoProngLoose_nNeutralIsoCone",&fTwoProngLoose_nNeutralIsoCone);
  fTree->Branch("TwoProngLoose_nEGammaIsoCone",&fTwoProngLoose_nEGammaIsoCone);
  // matching
  fTree->Branch("TwoProngLoose_genOmega_dR",&fTwoProngLoose_genOmega_dR);
  fTree->Branch("TwoProngLoose_genTau_dR",&fTwoProngLoose_genTau_dR);
  // dalitz
  if (fincludeDalitzVariables) {
  fTree->Branch("TwoProngLoose_mPosPho",&fTwoProngLoose_mPosPho);
  fTree->Branch("TwoProngLoose_mPosPho_l",&fTwoProngLoose_mPosPho_l);
  fTree->Branch("TwoProngLoose_mPosPho_pi0",&fTwoProngLoose_mPosPho_pi0);
  fTree->Branch("TwoProngLoose_mNegPho",&fTwoProngLoose_mNegPho);
  fTree->Branch("TwoProngLoose_mNegPho_l",&fTwoProngLoose_mNegPho_l);
  fTree->Branch("TwoProngLoose_mNegPho_pi0",&fTwoProngLoose_mNegPho_pi0);
  fTree->Branch("TwoProngLoose_mPosNeg",&fTwoProngLoose_mPosNeg);
  fTree->Branch("TwoProngLoose_CHpos_p3",&fTwoProngLoose_CHpos_p3);
  fTree->Branch("TwoProngLoose_CHneg_p3",&fTwoProngLoose_CHneg_p3);
  fTree->Branch("TwoProngLoose_photon_p3",&fTwoProngLoose_photon_p3);
  }
  }
  if(fincludeAsymTwoProngs) {
    // Asym TwoProng information
  fTree->Branch("nTwoProngsAsym",&fnTwoProngsAsym,"nTwoProngsAsym/I");
  // kinematics
  fTree->Branch("TwoProngAsym_pt",&fTwoProngAsym_pt);
  fTree->Branch("TwoProngAsym_eta",&fTwoProngAsym_eta);
  fTree->Branch("TwoProngAsym_phi",&fTwoProngAsym_phi);
  // mass
  fTree->Branch("TwoProngAsym_mass",&fTwoProngAsym_mass);
  fTree->Branch("TwoProngAsym_mass_l",&fTwoProngAsym_mass_l);
  fTree->Branch("TwoProngAsym_Mass0",&fTwoProngAsym_Mass0);
  fTree->Branch("TwoProngAsym_MassPi0",&fTwoProngAsym_MassPi0);
  fTree->Branch("TwoProngAsym_MassEta",&fTwoProngAsym_MassEta);
  fTree->Branch("TwoProngAsym_Mass300",&fTwoProngAsym_Mass300);
  // extra optional track 
  fTree->Branch("TwoProngAsym_nExtraTracks",&fTwoProngAsym_nExtraTracks);
  // CHpos constituent
  fTree->Branch("TwoProngAsym_CHpos_pt",&fTwoProngAsym_CHpos_pt);
  fTree->Branch("TwoProngAsym_CHpos_eta",&fTwoProngAsym_CHpos_eta);
  fTree->Branch("TwoProngAsym_CHpos_phi",&fTwoProngAsym_CHpos_phi);
  fTree->Branch("TwoProngAsym_CHpos_mass",&fTwoProngAsym_CHpos_mass);
  fTree->Branch("TwoProngAsym_CHpos_dz",&fTwoProngAsym_CHpos_dz);
  fTree->Branch("TwoProngAsym_CHpos_dxy",&fTwoProngAsym_CHpos_dxy);
  // CHneg constituent
  fTree->Branch("TwoProngAsym_CHneg_pt",&fTwoProngAsym_CHneg_pt);
  fTree->Branch("TwoProngAsym_CHneg_eta",&fTwoProngAsym_CHneg_eta);
  fTree->Branch("TwoProngAsym_CHneg_phi",&fTwoProngAsym_CHneg_phi);
  fTree->Branch("TwoProngAsym_CHneg_mass",&fTwoProngAsym_CHneg_mass);
  fTree->Branch("TwoProngAsym_CHneg_dz",&fTwoProngAsym_CHneg_dz);
  fTree->Branch("TwoProngAsym_CHneg_dxy",&fTwoProngAsym_CHneg_dxy);
  // photon constituent
  fTree->Branch("TwoProngAsym_photon_pt",&fTwoProngAsym_photon_pt);
  fTree->Branch("TwoProngAsym_photon_eta",&fTwoProngAsym_photon_eta);
  fTree->Branch("TwoProngAsym_photon_phi",&fTwoProngAsym_photon_phi);
  fTree->Branch("TwoProngAsym_photon_mass",&fTwoProngAsym_photon_mass);
  fTree->Branch("TwoProngAsym_photon_pt_l",&fTwoProngAsym_photon_pt_l);
  fTree->Branch("TwoProngAsym_photon_eta_l",&fTwoProngAsym_photon_eta_l);
  fTree->Branch("TwoProngAsym_photon_phi_l",&fTwoProngAsym_photon_phi_l);
  fTree->Branch("TwoProngAsym_photon_mass_l",&fTwoProngAsym_photon_mass_l);
  fTree->Branch("TwoProngAsym_photon_nGamma",&fTwoProngAsym_photon_nGamma);
  fTree->Branch("TwoProngAsym_photon_nElectron",&fTwoProngAsym_photon_nElectron);
  // cut variables
  fTree->Branch("TwoProngAsym_chargedIso",&fTwoProngAsym_chargedIso);
  fTree->Branch("TwoProngAsym_neutralIso",&fTwoProngAsym_neutralIso);
  fTree->Branch("TwoProngAsym_egammaIso",&fTwoProngAsym_egammaIso);
  fTree->Branch("TwoProngAsym_trackAsym",&fTwoProngAsym_trackAsym);
  fTree->Branch("TwoProngAsym_photonAsym",&fTwoProngAsym_photonAsym);
  fTree->Branch("TwoProngAsym_nChargedIsoCone",&fTwoProngAsym_nChargedIsoCone);
  fTree->Branch("TwoProngAsym_nNeutralIsoCone",&fTwoProngAsym_nNeutralIsoCone);
  fTree->Branch("TwoProngAsym_nEGammaIsoCone",&fTwoProngAsym_nEGammaIsoCone);
  // matching
  fTree->Branch("TwoProngAsym_genOmega_dR",&fTwoProngAsym_genOmega_dR);
  fTree->Branch("TwoProngAsym_genTau_dR",&fTwoProngAsym_genTau_dR);
  // dalitz
  if (fincludeDalitzVariables) {
  fTree->Branch("TwoProngAsym_mPosPho",&fTwoProngAsym_mPosPho);
  fTree->Branch("TwoProngAsym_mPosPho_l",&fTwoProngAsym_mPosPho_l);
  fTree->Branch("TwoProngAsym_mPosPho_pi0",&fTwoProngAsym_mPosPho_pi0);
  fTree->Branch("TwoProngAsym_mNegPho",&fTwoProngAsym_mNegPho);
  fTree->Branch("TwoProngAsym_mNegPho_l",&fTwoProngAsym_mNegPho_l);
  fTree->Branch("TwoProngAsym_mNegPho_pi0",&fTwoProngAsym_mNegPho_pi0);
  fTree->Branch("TwoProngAsym_mPosNeg",&fTwoProngAsym_mPosNeg);
  fTree->Branch("TwoProngAsym_CHpos_p3",&fTwoProngAsym_CHpos_p3);
  fTree->Branch("TwoProngAsym_CHneg_p3",&fTwoProngAsym_CHneg_p3);
  fTree->Branch("TwoProngAsym_photon_p3",&fTwoProngAsym_photon_p3);
  }
    // Asym Loose TwoProng information
  fTree->Branch("nTwoProngsAsymLoose",&fnTwoProngsAsymLoose,"nTwoProngsAsymLoose/I");
  // kinematics
  fTree->Branch("TwoProngAsymLoose_pt",&fTwoProngAsymLoose_pt);
  fTree->Branch("TwoProngAsymLoose_eta",&fTwoProngAsymLoose_eta);
  fTree->Branch("TwoProngAsymLoose_phi",&fTwoProngAsymLoose_phi);
  // mass
  fTree->Branch("TwoProngAsymLoose_mass",&fTwoProngAsymLoose_mass);
  fTree->Branch("TwoProngAsymLoose_mass_l",&fTwoProngAsymLoose_mass_l);
  fTree->Branch("TwoProngAsymLoose_Mass0",&fTwoProngAsymLoose_Mass0);
  fTree->Branch("TwoProngAsymLoose_MassPi0",&fTwoProngAsymLoose_MassPi0);
  fTree->Branch("TwoProngAsymLoose_MassEta",&fTwoProngAsymLoose_MassEta);
  fTree->Branch("TwoProngAsymLoose_Mass300",&fTwoProngAsymLoose_Mass300);
  // extra optional track 
  fTree->Branch("TwoProngAsymLoose_nExtraTracks",&fTwoProngAsymLoose_nExtraTracks);
  // CHpos constituent
  fTree->Branch("TwoProngAsymLoose_CHpos_pt",&fTwoProngAsymLoose_CHpos_pt);
  fTree->Branch("TwoProngAsymLoose_CHpos_eta",&fTwoProngAsymLoose_CHpos_eta);
  fTree->Branch("TwoProngAsymLoose_CHpos_phi",&fTwoProngAsymLoose_CHpos_phi);
  fTree->Branch("TwoProngAsymLoose_CHpos_mass",&fTwoProngAsymLoose_CHpos_mass);
  fTree->Branch("TwoProngAsymLoose_CHpos_dz",&fTwoProngAsymLoose_CHpos_dz);
  fTree->Branch("TwoProngAsymLoose_CHpos_dxy",&fTwoProngAsymLoose_CHpos_dxy);
  // CHneg constituent
  fTree->Branch("TwoProngAsymLoose_CHneg_pt",&fTwoProngAsymLoose_CHneg_pt);
  fTree->Branch("TwoProngAsymLoose_CHneg_eta",&fTwoProngAsymLoose_CHneg_eta);
  fTree->Branch("TwoProngAsymLoose_CHneg_phi",&fTwoProngAsymLoose_CHneg_phi);
  fTree->Branch("TwoProngAsymLoose_CHneg_mass",&fTwoProngAsymLoose_CHneg_mass);
  fTree->Branch("TwoProngAsymLoose_CHneg_dz",&fTwoProngAsymLoose_CHneg_dz);
  fTree->Branch("TwoProngAsymLoose_CHneg_dxy",&fTwoProngAsymLoose_CHneg_dxy);
  // photon constituent
  fTree->Branch("TwoProngAsymLoose_photon_pt",&fTwoProngAsymLoose_photon_pt);
  fTree->Branch("TwoProngAsymLoose_photon_eta",&fTwoProngAsymLoose_photon_eta);
  fTree->Branch("TwoProngAsymLoose_photon_phi",&fTwoProngAsymLoose_photon_phi);
  fTree->Branch("TwoProngAsymLoose_photon_mass",&fTwoProngAsymLoose_photon_mass);
  fTree->Branch("TwoProngAsymLoose_photon_pt_l",&fTwoProngAsymLoose_photon_pt_l);
  fTree->Branch("TwoProngAsymLoose_photon_eta_l",&fTwoProngAsymLoose_photon_eta_l);
  fTree->Branch("TwoProngAsymLoose_photon_phi_l",&fTwoProngAsymLoose_photon_phi_l);
  fTree->Branch("TwoProngAsymLoose_photon_mass_l",&fTwoProngAsymLoose_photon_mass_l);
  fTree->Branch("TwoProngAsymLoose_photon_nGamma",&fTwoProngAsymLoose_photon_nGamma);
  fTree->Branch("TwoProngAsymLoose_photon_nElectron",&fTwoProngAsymLoose_photon_nElectron);
  // cut variables
  fTree->Branch("TwoProngAsymLoose_chargedIso",&fTwoProngAsymLoose_chargedIso);
  fTree->Branch("TwoProngAsymLoose_neutralIso",&fTwoProngAsymLoose_neutralIso);
  fTree->Branch("TwoProngAsymLoose_egammaIso",&fTwoProngAsymLoose_egammaIso);
  fTree->Branch("TwoProngAsymLoose_trackAsym",&fTwoProngAsymLoose_trackAsym);
  fTree->Branch("TwoProngAsymLoose_photonAsym",&fTwoProngAsymLoose_photonAsym);
  fTree->Branch("TwoProngAsymLoose_nChargedIsoCone",&fTwoProngAsymLoose_nChargedIsoCone);
  fTree->Branch("TwoProngAsymLoose_nNeutralIsoCone",&fTwoProngAsymLoose_nNeutralIsoCone);
  fTree->Branch("TwoProngAsymLoose_nEGammaIsoCone",&fTwoProngAsymLoose_nEGammaIsoCone);
  // matching
  fTree->Branch("TwoProngAsymLoose_genOmega_dR",&fTwoProngAsymLoose_genOmega_dR);
  fTree->Branch("TwoProngAsymLoose_genTau_dR",&fTwoProngAsymLoose_genTau_dR);
  // dalitz
  if (fincludeDalitzVariables) {
  fTree->Branch("TwoProngAsymLoose_mPosPho",&fTwoProngAsymLoose_mPosPho);
  fTree->Branch("TwoProngAsymLoose_mPosPho_l",&fTwoProngAsymLoose_mPosPho_l);
  fTree->Branch("TwoProngAsymLoose_mPosPho_pi0",&fTwoProngAsymLoose_mPosPho_pi0);
  fTree->Branch("TwoProngAsymLoose_mNegPho",&fTwoProngAsymLoose_mNegPho);
  fTree->Branch("TwoProngAsymLoose_mNegPho_l",&fTwoProngAsymLoose_mNegPho_l);
  fTree->Branch("TwoProngAsymLoose_mNegPho_pi0",&fTwoProngAsymLoose_mNegPho_pi0);
  fTree->Branch("TwoProngAsymLoose_mPosNeg",&fTwoProngAsymLoose_mPosNeg);
  fTree->Branch("TwoProngAsymLoose_CHpos_p3",&fTwoProngAsymLoose_CHpos_p3);
  fTree->Branch("TwoProngAsymLoose_CHneg_p3",&fTwoProngAsymLoose_CHneg_p3);
  fTree->Branch("TwoProngAsymLoose_photon_p3",&fTwoProngAsymLoose_photon_p3);
  }
  }
  // Combined Objects
  fTree->Branch("sCat",&fSCat,"sCat/I");
  fTree->Branch("Obj_DiTwoProng",&fRecoPhiDiTwoProng,TwoProngAnalysis::recoDiObjectBranchDefString.c_str());
  fTree->Branch("Obj_PhotonTwoProng",&fRecoPhiPhotonTwoProng,TwoProngAnalysis::recoDiObjectBranchDefString.c_str());
  fTree->Branch("Obj_RecoPhiInclusive",&fRecoPhiInclusive,TwoProngAnalysis::recoDiObjectBranchDefString.c_str());
  }
  // Tau preseletion branches
  if (fincludeZTauHadBranches) {
  fTree->Branch("passMuonTrigger",&fpassMuonTrigger,"passMuonTrigger/O");
  fTree->Branch("passMuonTriggerTk",&fpassMuonTriggerTk,"passMuonTriggerTk/O");
  fTree->Branch("muonTrigger",&fMuonTrigger);
  fTree->Branch("muonTriggerTk",&fMuonTriggerTk);
  fTree->Branch("nTagMuons",&fnTagMuons,"nTagMuons/I");
  fTree->Branch("nProbeTaus",&fnProbeTaus,"nProbeTaus/I");
  fTree->Branch("passPzeta",&fpassPzeta,"passPzeta/O");
  fTree->Branch("passMT",&fpassMT,"passMT/O");
  fTree->Branch("passExtraElectronVeto",&fpassExtraElectronVeto,"passExtraElectronVeto/O");
  fTree->Branch("passExtraMuonVeto",&fpassExtraMuonVeto,"passExtraMuonVeto/O");
  fTree->Branch("passDiMuonVeto",&fpassDiMuonVeto,"passDiMuonVeto/O");
  fTree->Branch("passExtraLeptonVeto",&fpassExtraLeptonVeto,"passExtraLeptonVeto/O");
  fTree->Branch("passBtagVeto",&fpassBtagVeto,"passBtagVeto/O");
  fTree->Branch("MT",&fMT,"MT/D");
  fTree->Branch("Pzeta",&fPzeta,"Pzeta/D");
  fTree->Branch("highestBtag",&fhighestBtag,"highestBtag/D");
  fTree->Branch("passMuonTauPair",&fpassMuonTauPair,"passMuonTauPair/O");
  fTree->Branch("passPreselection",&fpassPreselection,"passPreselection/O");
  fTree->Branch("passReducedSelection",&fpassReducedSelection,"passReducedSelection/O");
  fTree->Branch("TagMuon_pt",&fTagMuon_pt,"TagMuon_pt/D");
  fTree->Branch("TagMuon_eta",&fTagMuon_eta,"TagMuon_eta/D");
  fTree->Branch("TagMuon_phi",&fTagMuon_phi,"TagMuon_phi/D");
  fTree->Branch("TagMuon_mass",&fTagMuon_mass,"TagMuon_mass/D");
  fTree->Branch("TagMuon_z",&fTagMuon_z,"TagMuon_z/D");
  fTree->Branch("TagMuon_dz",&fTagMuon_dz,"TagMuon_dz/D");
  fTree->Branch("TagMuon_dB",&fTagMuon_dB,"TagMuon_dB/D");
  fTree->Branch("TagMuon_dxy",&fTagMuon_dxy,"TagMuon_dxy/D");
  fTree->Branch("TagMuon_iso",&fTagMuon_iso,"TagMuon_iso/D");
  fTree->Branch("ProbeTau_pt",&fProbeTau_pt,"ProbeTau_pt/D");
  fTree->Branch("ProbeTau_eta",&fProbeTau_eta,"ProbeTau_eta/D");
  fTree->Branch("ProbeTau_phi",&fProbeTau_phi,"ProbeTau_phi/D");
  fTree->Branch("ProbeTau_mass",&fProbeTau_mass,"ProbeTau_mass/D");
  fTree->Branch("ProbeTau_genDR",&fProbeTau_genDR,"ProbeTau_genDR/D");
  fTree->Branch("probePassTauID",&fprobePassTauID,"probePassTauID/O");
  fTree->Branch("ProbeTau_tauID_pt",&fProbeTau_tauID_pt,"ProbeTau_tauID_pt/D");
  fTree->Branch("ProbeTau_tauID_eta",&fProbeTau_tauID_eta,"ProbeTau_tauID_eta/D");
  fTree->Branch("ProbeTau_tauID_phi",&fProbeTau_tauID_phi,"ProbeTau_tauID_phi/D");
  fTree->Branch("ProbeTau_tauID_mass",&fProbeTau_tauID_mass,"ProbeTau_tauID_mass/D");
  fTree->Branch("ProbeTau_tauID_probeDR",&fProbeTau_tauID_probeDR,"ProbeTau_tauID_probeDR/D");
  fTree->Branch("probePassTwoProng",&fprobePassTwoProng,"probePassTwoProng/O");
  fTree->Branch("ProbeTau_twoprong_pt",&fProbeTau_twoprong_pt,"ProbeTau_twoprong_pt/D");
  fTree->Branch("ProbeTau_twoprong_eta",&fProbeTau_twoprong_eta,"ProbeTau_twoprong_eta/D");
  fTree->Branch("ProbeTau_twoprong_phi",&fProbeTau_twoprong_phi,"ProbeTau_twoprong_phi/D");
  fTree->Branch("ProbeTau_twoprong_mass",&fProbeTau_twoprong_mass,"ProbeTau_twoprong_mass/D");
  fTree->Branch("ProbeTau_twoprong_probeDR",&fProbeTau_twoprong_probeDR,"ProbeTau_twoprong_probeDR/D");
  fTree->Branch("Obj_MuonProbe",&fMuonProbe,TwoProngAnalysis::recoDiObjectBranchDefString.c_str());
  fTree->Branch("Obj_MuonTauID",&fMuonTauID,TwoProngAnalysis::recoDiObjectBranchDefString.c_str());
  fTree->Branch("Obj_MuonTwoProng",&fMuonTwoProng,TwoProngAnalysis::recoDiObjectBranchDefString.c_str());
  fTree->Branch("Obj_MuonTauIDpt",&fMuonTauID_pt,TwoProngAnalysis::recoDiObjectBranchDefString.c_str());
  fTree->Branch("Obj_MuonTauIDzmass",&fMuonTauID_zmass,TwoProngAnalysis::recoDiObjectBranchDefString.c_str());
  fTree->Branch("Obj_MuonTwoProngpt",&fMuonTwoProng_pt,TwoProngAnalysis::recoDiObjectBranchDefString.c_str());
  fTree->Branch("Obj_MuonTwoProngzmass",&fMuonTwoProng_zmass,TwoProngAnalysis::recoDiObjectBranchDefString.c_str());
  }
  if (fincludeZMuMuBranches) {
  fTree->Branch("passMuonTrigger",&fpassMuonTrigger,"passMuonTrigger/O");
  fTree->Branch("passMuonTriggerTk",&fpassMuonTriggerTk,"passMuonTriggerTk/O");
  fTree->Branch("muonTrigger",&fMuonTrigger);
  fTree->Branch("muonTriggerTk",&fMuonTriggerTk);
  fTree->Branch("nTagMuons",&fnTagMuons,"nTagMuons/I");
  fTree->Branch("passExtraElectronVeto",&fpassExtraElectronVeto,"passExtraElectronVeto/O");
  fTree->Branch("passExtraMuonVeto",&fpassExtraMuonVeto,"passExtraMuonVeto/O");
  fTree->Branch("passPreselection",&fpassPreselection,"passPreselection/O");
  fTree->Branch("passReducedSelection",&fpassReducedSelection,"passReducedSelection/O");
  fTree->Branch("TagMuon1_pt",&fTagMuon1_pt,"TagMuon1_pt/D");
  fTree->Branch("TagMuon1_eta",&fTagMuon1_eta,"TagMuon1_eta/D");
  fTree->Branch("TagMuon1_phi",&fTagMuon1_phi,"TagMuon1_phi/D");
  fTree->Branch("TagMuon1_mass",&fTagMuon1_mass,"TagMuon1_mass/D");
  fTree->Branch("TagMuon1_dz",&fTagMuon1_dz,"TagMuon1_dz/D");
  fTree->Branch("TagMuon1_iso",&fTagMuon1_iso,"TagMuon1_iso/D");
  fTree->Branch("TagMuon2_pt",&fTagMuon2_pt,"TagMuon2_pt/D");
  fTree->Branch("TagMuon2_eta",&fTagMuon2_eta,"TagMuon2_eta/D");
  fTree->Branch("TagMuon2_phi",&fTagMuon2_phi,"TagMuon2_phi/D");
  fTree->Branch("TagMuon2_mass",&fTagMuon2_mass,"TagMuon2_mass/D");
  fTree->Branch("TagMuon2_dz",&fTagMuon2_dz,"TagMuon2_dz/D");
  fTree->Branch("TagMuon2_iso",&fTagMuon2_iso,"TagMuon2_iso/D");
  fTree->Branch("Obj_MuonMuon",&fMuonMuon,TwoProngAnalysis::recoDiObjectBranchDefString.c_str());
  }
  }
  if (fincludeLeptonBranches) {
   fTree->Branch("nTightMuons",&fNTightMuons,"nTightMuons/I");
   fTree->Branch("TightMuon_pt",&fTightMuon_pt);
   fTree->Branch("TightMuon_eta",&fTightMuon_eta);
   fTree->Branch("TightMuon_phi",&fTightMuon_phi);
   fTree->Branch("TightMuon_mass",&fTightMuon_mass);
   fTree->Branch("Muon_veto", &fMuon_veto, "Muon_veto/I");
   fTree->Branch("btag_veto",&fbtag_veto,"btag_veto/I");
   fTree->Branch("Imag_W", &fImag_W, "Imag_W/I");
   fTree->Branch("W_pt",&fW_pt);
   fTree->Branch("W_eta",&fW_eta);
   fTree->Branch("W_phi",&fW_phi);
   fTree->Branch("W_mass",&fW_mass);
   fTree->Branch("W_mT",&fW_mT);
   fTree->Branch("mT",&fmT,"mT/D");
   fTree->Branch("Num_Muons",&fNum_Muons,"Num_Muons/I"); // counts muons that pass pt and eta cuts, but not rest of muon ID, keeping in case stephen needs this
  }

  if(fincludeDalitzHistos) {
  fHighvsMid = fs->make<TH2D>("highvsmid","highvsmid",40,0,1,40,0,1);
  fHighvsLow = fs->make<TH2D>("highvslow","highvslow",40,0,1,40,0,1);
  fMidvsLow = fs->make<TH2D>("midvslow","midvslow",40,0,1,40,0,1);
  fPhotonvsLarger = fs->make<TH2D>("photonvslarger","photonvslarger",40,0,1,40,0,1);
  fPhotonvsSmaller = fs->make<TH2D>("photonvssmaller","photonvssmaller",40,0,1,40,0,1);
  fDoubleStackedDalitz = fs->make<TH2D>("dalitz_double","dalitz_double",40,0,1,40,0,1);
  fTripleStackedDalitz = fs->make<TH2D>("dalitz_triple","dalitz_triple",40,0,1,40,0,1);
  
  fPhotonvsPositive = fs->make<TH2D>("dalitz_photonvspos","dalitz_photonvspos",40,0,1,40,0,1);
  fPhotonvsNegative = fs->make<TH2D>("dalitz_photonvsneg","dalitz_photonvsneg",40,0,1,40,0,1);
  fPositivevsNegative = fs->make<TH2D>("dalitz_posvsneg","dalitz_posvsneg",40,0,1,40,0,1);

  fPhotonFraction = fs->make<TH1F>("1d_photon_frac","1d_pos_frac",40,0,1);
  fPositiveFraction = fs->make<TH1F>("1d_pos_frac","1d_pos_frac",40,0,1);
  fNegativeFraction = fs->make<TH1F>("1d_negative_frac","1d_pos_frac",40,0,1);
  fHTverify = fs->make<TH1F>("ht","ht",500,0,5000);
  }
  if (fDebug) cout << "Exiting constructor... " << endl;
}


TwoProngAnalyzer::~TwoProngAnalyzer()
{
 
}


void
TwoProngAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if (fDebug) cout << "Entering analyze()... " << endl;

  using namespace edm;
  using namespace std;
  using namespace reco;

  if (fDebug) cout << ". event " << iEvent.id().event() << " lumi " <<  iEvent.id().luminosityBlock() << " run " << iEvent.id().run() << endl;

  // clear member vectors
  fTwoProngCand_pt.clear();
  fTwoProngCand_eta.clear();
  fTwoProngCand_phi.clear();
  fTwoProngCand_mass.clear();
  fTwoProngCand_mass_l.clear();
  fTwoProngCand_Mass0.clear();
  fTwoProngCand_MassPi0.clear();
  fTwoProngCand_MassEta.clear();
  fTwoProngCand_Mass300.clear();
  fTwoProngCand_nExtraTracks.clear();
  fTwoProngCand_CHpos_pt.clear();
  fTwoProngCand_CHpos_eta.clear();
  fTwoProngCand_CHpos_phi.clear();
  fTwoProngCand_CHpos_mass.clear();
  fTwoProngCand_CHpos_dz.clear();
  fTwoProngCand_CHpos_dxy.clear();
  fTwoProngCand_CHneg_pt.clear();
  fTwoProngCand_CHneg_eta.clear();
  fTwoProngCand_CHneg_phi.clear();
  fTwoProngCand_CHneg_mass.clear();
  fTwoProngCand_CHneg_dz.clear();
  fTwoProngCand_CHneg_dxy.clear();
  fTwoProngCand_photon_pt.clear();
  fTwoProngCand_photon_eta.clear();
  fTwoProngCand_photon_phi.clear();
  fTwoProngCand_photon_mass.clear();
  fTwoProngCand_photon_pt_l.clear();
  fTwoProngCand_photon_eta_l.clear();
  fTwoProngCand_photon_phi_l.clear();
  fTwoProngCand_photon_mass_l.clear();
  fTwoProngCand_photon_nGamma.clear();
  fTwoProngCand_photon_nElectron.clear();
  fTwoProngCand_chargedIso.clear();
  fTwoProngCand_neutralIso.clear();
  fTwoProngCand_egammaIso.clear();
  fTwoProngCand_trackAsym.clear();
  fTwoProngCand_photonAsym.clear();
  fTwoProngCand_nChargedIsoCone.clear();
  fTwoProngCand_nNeutralIsoCone.clear();
  fTwoProngCand_nEGammaIsoCone.clear();
  fTwoProngCand_tight.clear();
  fTwoProngCand_loose.clear();
  fTwoProngCand_asym.clear();
  fTwoProngCand_asym_loose.clear();
  fTwoProngCand_genOmega_dR.clear();
  fTwoProngCand_genTau_dR.clear();
  fTwoProngCand_mPosPho.clear();
  fTwoProngCand_mPosPho_l.clear();
  fTwoProngCand_mPosPho_pi0.clear();
  fTwoProngCand_mNegPho.clear();
  fTwoProngCand_mNegPho_l.clear();
  fTwoProngCand_mNegPho_pi0.clear();
  fTwoProngCand_mPosNeg.clear();
  fTwoProngCand_CHpos_p3.clear();
  fTwoProngCand_CHneg_p3.clear();
  fTwoProngCand_photon_p3.clear();

  fTwoProngLoose_pt.clear();
  fTwoProngLoose_eta.clear();
  fTwoProngLoose_phi.clear();
  fTwoProngLoose_mass.clear();
  fTwoProngLoose_mass_l.clear();
  fTwoProngLoose_Mass0.clear();
  fTwoProngLoose_MassPi0.clear();
  fTwoProngLoose_MassEta.clear();
  fTwoProngLoose_Mass300.clear();
  fTwoProngLoose_nExtraTracks.clear();
  fTwoProngLoose_CHpos_pt.clear();
  fTwoProngLoose_CHpos_eta.clear();
  fTwoProngLoose_CHpos_phi.clear();
  fTwoProngLoose_CHpos_mass.clear();
  fTwoProngLoose_CHpos_dz.clear();
  fTwoProngLoose_CHpos_dxy.clear();
  fTwoProngLoose_CHneg_pt.clear();
  fTwoProngLoose_CHneg_eta.clear();
  fTwoProngLoose_CHneg_phi.clear();
  fTwoProngLoose_CHneg_mass.clear();
  fTwoProngLoose_CHneg_dz.clear();
  fTwoProngLoose_CHneg_dxy.clear();
  fTwoProngLoose_photon_pt.clear();
  fTwoProngLoose_photon_eta.clear();
  fTwoProngLoose_photon_phi.clear();
  fTwoProngLoose_photon_mass.clear();
  fTwoProngLoose_photon_pt_l.clear();
  fTwoProngLoose_photon_eta_l.clear();
  fTwoProngLoose_photon_phi_l.clear();
  fTwoProngLoose_photon_mass_l.clear();
  fTwoProngLoose_photon_nGamma.clear();
  fTwoProngLoose_photon_nElectron.clear();
  fTwoProngLoose_chargedIso.clear();
  fTwoProngLoose_neutralIso.clear();
  fTwoProngLoose_egammaIso.clear();
  fTwoProngLoose_trackAsym.clear();
  fTwoProngLoose_photonAsym.clear();
  fTwoProngLoose_nChargedIsoCone.clear();
  fTwoProngLoose_nNeutralIsoCone.clear();
  fTwoProngLoose_nEGammaIsoCone.clear();
  fTwoProngLoose_tight.clear();
  fTwoProngLoose_loose.clear();
  fTwoProngLoose_genOmega_dR.clear();
  fTwoProngLoose_genTau_dR.clear();
  fTwoProngLoose_mPosPho.clear();
  fTwoProngLoose_mPosPho_l.clear();
  fTwoProngLoose_mPosPho_pi0.clear();
  fTwoProngLoose_mNegPho.clear();
  fTwoProngLoose_mNegPho_l.clear();
  fTwoProngLoose_mNegPho_pi0.clear();
  fTwoProngLoose_mPosNeg.clear();
  fTwoProngLoose_CHpos_p3.clear();
  fTwoProngLoose_CHneg_p3.clear();
  fTwoProngLoose_photon_p3.clear();

  fTwoProngAsym_pt.clear();
  fTwoProngAsym_eta.clear();
  fTwoProngAsym_phi.clear();
  fTwoProngAsym_mass.clear();
  fTwoProngAsym_mass_l.clear();
  fTwoProngAsym_Mass0.clear();
  fTwoProngAsym_MassPi0.clear();
  fTwoProngAsym_MassEta.clear();
  fTwoProngAsym_Mass300.clear();
  fTwoProngAsym_nExtraTracks.clear();
  fTwoProngAsym_CHpos_pt.clear();
  fTwoProngAsym_CHpos_eta.clear();
  fTwoProngAsym_CHpos_phi.clear();
  fTwoProngAsym_CHpos_mass.clear();
  fTwoProngAsym_CHpos_dz.clear();
  fTwoProngAsym_CHpos_dxy.clear();
  fTwoProngAsym_CHneg_pt.clear();
  fTwoProngAsym_CHneg_eta.clear();
  fTwoProngAsym_CHneg_phi.clear();
  fTwoProngAsym_CHneg_mass.clear();
  fTwoProngAsym_CHneg_dz.clear();
  fTwoProngAsym_CHneg_dxy.clear();
  fTwoProngAsym_photon_pt.clear();
  fTwoProngAsym_photon_eta.clear();
  fTwoProngAsym_photon_phi.clear();
  fTwoProngAsym_photon_mass.clear();
  fTwoProngAsym_photon_pt_l.clear();
  fTwoProngAsym_photon_eta_l.clear();
  fTwoProngAsym_photon_phi_l.clear();
  fTwoProngAsym_photon_mass_l.clear();
  fTwoProngAsym_photon_nGamma.clear();
  fTwoProngAsym_photon_nElectron.clear();
  fTwoProngAsym_chargedIso.clear();
  fTwoProngAsym_neutralIso.clear();
  fTwoProngAsym_egammaIso.clear();
  fTwoProngAsym_trackAsym.clear();
  fTwoProngAsym_photonAsym.clear();
  fTwoProngAsym_nChargedIsoCone.clear();
  fTwoProngAsym_nNeutralIsoCone.clear();
  fTwoProngAsym_nEGammaIsoCone.clear();
  fTwoProngAsym_tight.clear();
  fTwoProngAsym_loose.clear();
  fTwoProngAsym_genOmega_dR.clear();
  fTwoProngAsym_genTau_dR.clear();
  fTwoProngAsym_mPosPho.clear();
  fTwoProngAsym_mPosPho_l.clear();
  fTwoProngAsym_mPosPho_pi0.clear();
  fTwoProngAsym_mNegPho.clear();
  fTwoProngAsym_mNegPho_l.clear();
  fTwoProngAsym_mNegPho_pi0.clear();
  fTwoProngAsym_mPosNeg.clear();
  fTwoProngAsym_CHpos_p3.clear();
  fTwoProngAsym_CHneg_p3.clear();
  fTwoProngAsym_photon_p3.clear();

  fTwoProngAsymLoose_pt.clear();
  fTwoProngAsymLoose_eta.clear();
  fTwoProngAsymLoose_phi.clear();
  fTwoProngAsymLoose_mass.clear();
  fTwoProngAsymLoose_mass_l.clear();
  fTwoProngAsymLoose_Mass0.clear();
  fTwoProngAsymLoose_MassPi0.clear();
  fTwoProngAsymLoose_MassEta.clear();
  fTwoProngAsymLoose_Mass300.clear();
  fTwoProngAsymLoose_nExtraTracks.clear();
  fTwoProngAsymLoose_CHpos_pt.clear();
  fTwoProngAsymLoose_CHpos_eta.clear();
  fTwoProngAsymLoose_CHpos_phi.clear();
  fTwoProngAsymLoose_CHpos_mass.clear();
  fTwoProngAsymLoose_CHpos_dz.clear();
  fTwoProngAsymLoose_CHpos_dxy.clear();
  fTwoProngAsymLoose_CHneg_pt.clear();
  fTwoProngAsymLoose_CHneg_eta.clear();
  fTwoProngAsymLoose_CHneg_phi.clear();
  fTwoProngAsymLoose_CHneg_mass.clear();
  fTwoProngAsymLoose_CHneg_dz.clear();
  fTwoProngAsymLoose_CHneg_dxy.clear();
  fTwoProngAsymLoose_photon_pt.clear();
  fTwoProngAsymLoose_photon_eta.clear();
  fTwoProngAsymLoose_photon_phi.clear();
  fTwoProngAsymLoose_photon_mass.clear();
  fTwoProngAsymLoose_photon_pt_l.clear();
  fTwoProngAsymLoose_photon_eta_l.clear();
  fTwoProngAsymLoose_photon_phi_l.clear();
  fTwoProngAsymLoose_photon_mass_l.clear();
  fTwoProngAsymLoose_photon_nGamma.clear();
  fTwoProngAsymLoose_photon_nElectron.clear();
  fTwoProngAsymLoose_chargedIso.clear();
  fTwoProngAsymLoose_neutralIso.clear();
  fTwoProngAsymLoose_egammaIso.clear();
  fTwoProngAsymLoose_trackAsym.clear();
  fTwoProngAsymLoose_photonAsym.clear();
  fTwoProngAsymLoose_nChargedIsoCone.clear();
  fTwoProngAsymLoose_nNeutralIsoCone.clear();
  fTwoProngAsymLoose_nEGammaIsoCone.clear();
  fTwoProngAsymLoose_tight.clear();
  fTwoProngAsymLoose_loose.clear();
  fTwoProngAsymLoose_genOmega_dR.clear();
  fTwoProngAsymLoose_genTau_dR.clear();
  fTwoProngAsymLoose_mPosPho.clear();
  fTwoProngAsymLoose_mPosPho_l.clear();
  fTwoProngAsymLoose_mPosPho_pi0.clear();
  fTwoProngAsymLoose_mNegPho.clear();
  fTwoProngAsymLoose_mNegPho_l.clear();
  fTwoProngAsymLoose_mNegPho_pi0.clear();
  fTwoProngAsymLoose_mPosNeg.clear();
  fTwoProngAsymLoose_CHpos_p3.clear();
  fTwoProngAsymLoose_CHneg_p3.clear();
  fTwoProngAsymLoose_photon_p3.clear();

  fTwoProng_pt.clear();
  fTwoProng_eta.clear();
  fTwoProng_phi.clear();
  fTwoProng_mass.clear();
  fTwoProng_mass_l.clear();
  fTwoProng_Mass0.clear();
  fTwoProng_MassPi0.clear();
  fTwoProng_MassEta.clear();
  fTwoProng_Mass300.clear();
  fTwoProng_nExtraTracks.clear();
  fTwoProng_CHpos_pt.clear();
  fTwoProng_CHpos_eta.clear();
  fTwoProng_CHpos_phi.clear();
  fTwoProng_CHpos_mass.clear();
  fTwoProng_CHpos_dz.clear();
  fTwoProng_CHpos_dxy.clear();
  fTwoProng_CHneg_pt.clear();
  fTwoProng_CHneg_eta.clear();
  fTwoProng_CHneg_phi.clear();
  fTwoProng_CHneg_mass.clear();
  fTwoProng_CHneg_dz.clear();
  fTwoProng_CHneg_dxy.clear();
  fTwoProng_photon_pt.clear();
  fTwoProng_photon_eta.clear();
  fTwoProng_photon_phi.clear();
  fTwoProng_photon_mass.clear();
  fTwoProng_photon_pt_l.clear();
  fTwoProng_photon_eta_l.clear();
  fTwoProng_photon_phi_l.clear();
  fTwoProng_photon_mass_l.clear();
  fTwoProng_photon_nGamma.clear();
  fTwoProng_photon_nElectron.clear();
  fTwoProng_chargedIso.clear();
  fTwoProng_neutralIso.clear();
  fTwoProng_egammaIso.clear();
  fTwoProng_trackAsym.clear();
  fTwoProng_photonAsym.clear();
  fTwoProng_nChargedIsoCone.clear();
  fTwoProng_nNeutralIsoCone.clear();
  fTwoProng_nEGammaIsoCone.clear();
  fTwoProng_tight.clear();
  fTwoProng_loose.clear();
  fTwoProng_genOmega_dR.clear();
  fTwoProng_genTau_dR.clear();
  fTwoProng_mPosPho.clear();
  fTwoProng_mPosPho_l.clear();
  fTwoProng_mPosPho_pi0.clear();
  fTwoProng_mNegPho.clear();
  fTwoProng_mNegPho_l.clear();
  fTwoProng_mNegPho_pi0.clear();
  fTwoProng_mPosNeg.clear();
  fTwoProng_CHpos_p3.clear();
  fTwoProng_CHneg_p3.clear();
  fTwoProng_photon_p3.clear();

  fAK4jet_pt.clear();
  fAK4jet_eta.clear();
  fAK4jet_phi.clear();
  fAK4jet_mass.clear();

  fElectron_pt.clear();
  fElectron_eta.clear();
  fElectron_phi.clear();
  fElectron_mass.clear();

  fMuon_pt.clear();
  fMuon_eta.clear();
  fMuon_phi.clear();
  fMuon_mass.clear();

  fTau_pt.clear();
  fTau_eta.clear();
  fTau_phi.clear();
  fTau_mass.clear();

  fPhoton_pt.clear();
  fPhoton_eta.clear();
  fPhoton_scEta.clear();
  fPhoton_phi.clear();
  fPhoton_mass.clear();
  fPhoton_isoGamma.clear();
  fPhoton_isoCh.clear();
  fPhoton_HE.clear();
  fPhoton_coneHE.clear();
  fPhoton_sigmaIetaIeta.clear();
  fPhoton_passVeto.clear();

  fBaseIDPhoton_pt.clear();
  fBaseIDPhoton_eta.clear();
  fBaseIDPhoton_scEta.clear();
  fBaseIDPhoton_phi.clear();
  fBaseIDPhoton_mass.clear();
  fBaseIDPhoton_isoGamma.clear();
  fBaseIDPhoton_isoCh.clear();
  fBaseIDPhoton_HE.clear();
  fBaseIDPhoton_coneHE.clear();
  fBaseIDPhoton_sigmaIetaIeta.clear();

  fLoose1IDPhoton_pt.clear();
  fLoose1IDPhoton_eta.clear();
  fLoose1IDPhoton_scEta.clear();
  fLoose1IDPhoton_phi.clear();
  fLoose1IDPhoton_mass.clear();
  fLoose1IDPhoton_isoGamma.clear();
  fLoose1IDPhoton_isoCh.clear();
  fLoose1IDPhoton_HE.clear();
  fLoose1IDPhoton_coneHE.clear();
  fLoose1IDPhoton_sigmaIetaIeta.clear();

  fLoose1IDPhotonEndcap_pt.clear();
  fLoose1IDPhotonEndcap_eta.clear();
  fLoose1IDPhotonEndcap_scEta.clear();
  fLoose1IDPhotonEndcap_phi.clear();
  fLoose1IDPhotonEndcap_mass.clear();
  fLoose1IDPhotonEndcap_isoGamma.clear();
  fLoose1IDPhotonEndcap_isoCh.clear();
  fLoose1IDPhotonEndcap_HE.clear();
  fLoose1IDPhotonEndcap_coneHE.clear();
  fLoose1IDPhotonEndcap_sigmaIetaIeta.clear();

  fLoose2IDPhoton_pt.clear();
  fLoose2IDPhoton_eta.clear();
  fLoose2IDPhoton_scEta.clear();
  fLoose2IDPhoton_phi.clear();
  fLoose2IDPhoton_mass.clear();
  fLoose2IDPhoton_isoGamma.clear();
  fLoose2IDPhoton_isoCh.clear();
  fLoose2IDPhoton_HE.clear();
  fLoose2IDPhoton_coneHE.clear();
  fLoose2IDPhoton_sigmaIetaIeta.clear();

  fLoose2IDPhotonEndcap_pt.clear();
  fLoose2IDPhotonEndcap_eta.clear();
  fLoose2IDPhotonEndcap_scEta.clear();
  fLoose2IDPhotonEndcap_phi.clear();
  fLoose2IDPhotonEndcap_mass.clear();
  fLoose2IDPhotonEndcap_isoGamma.clear();
  fLoose2IDPhotonEndcap_isoCh.clear();
  fLoose2IDPhotonEndcap_HE.clear();
  fLoose2IDPhotonEndcap_coneHE.clear();
  fLoose2IDPhotonEndcap_sigmaIetaIeta.clear();

  fLoose3IDPhoton_pt.clear();
  fLoose3IDPhoton_eta.clear();
  fLoose3IDPhoton_scEta.clear();
  fLoose3IDPhoton_phi.clear();
  fLoose3IDPhoton_mass.clear();
  fLoose3IDPhoton_isoGamma.clear();
  fLoose3IDPhoton_isoCh.clear();
  fLoose3IDPhoton_HE.clear();
  fLoose3IDPhoton_coneHE.clear();
  fLoose3IDPhoton_sigmaIetaIeta.clear();

  fLoose4IDPhoton_pt.clear();
  fLoose4IDPhoton_eta.clear();
  fLoose4IDPhoton_scEta.clear();
  fLoose4IDPhoton_phi.clear();
  fLoose4IDPhoton_mass.clear();
  fLoose4IDPhoton_isoGamma.clear();
  fLoose4IDPhoton_isoCh.clear();
  fLoose4IDPhoton_HE.clear();
  fLoose4IDPhoton_coneHE.clear();
  fLoose4IDPhoton_sigmaIetaIeta.clear();

  fLoose5IDPhoton_pt.clear();
  fLoose5IDPhoton_eta.clear();
  fLoose5IDPhoton_scEta.clear();
  fLoose5IDPhoton_phi.clear();
  fLoose5IDPhoton_mass.clear();
  fLoose5IDPhoton_isoGamma.clear();
  fLoose5IDPhoton_isoCh.clear();
  fLoose5IDPhoton_HE.clear();
  fLoose5IDPhoton_coneHE.clear();
  fLoose5IDPhoton_sigmaIetaIeta.clear();
  fLoose5IDPhoton_passVeto.clear();

  fIDPhoton_pt.clear();
  fIDPhoton_eta.clear();
  fIDPhoton_scEta.clear();
  fIDPhoton_phi.clear();
  fIDPhoton_mass.clear();
  fIDPhoton_isoGamma.clear();
  fIDPhoton_isoCh.clear();
  fIDPhoton_HE.clear();
  fIDPhoton_coneHE.clear();
  fIDPhoton_sigmaIetaIeta.clear();

  fIDPhotonEndcap_pt.clear();
  fIDPhotonEndcap_eta.clear();
  fIDPhotonEndcap_scEta.clear();
  fIDPhotonEndcap_phi.clear();
  fIDPhotonEndcap_mass.clear();
  fIDPhotonEndcap_isoGamma.clear();
  fIDPhotonEndcap_isoCh.clear();
  fIDPhotonEndcap_HE.clear();
  fIDPhotonEndcap_coneHE.clear();
  fIDPhotonEndcap_sigmaIetaIeta.clear();

  fConeHEIDPhoton_pt.clear();
  fConeHEIDPhoton_eta.clear();
  fConeHEIDPhoton_scEta.clear();
  fConeHEIDPhoton_phi.clear();
  fConeHEIDPhoton_mass.clear();
  fConeHEIDPhoton_isoGamma.clear();
  fConeHEIDPhoton_isoCh.clear();
  fConeHEIDPhoton_HE.clear();
  fConeHEIDPhoton_coneHE.clear();
  fConeHEIDPhoton_sigmaIetaIeta.clear();

  fConeHEIDPhotonEndcap_pt.clear();
  fConeHEIDPhotonEndcap_eta.clear();
  fConeHEIDPhotonEndcap_scEta.clear();
  fConeHEIDPhotonEndcap_phi.clear();
  fConeHEIDPhotonEndcap_mass.clear();
  fConeHEIDPhotonEndcap_isoGamma.clear();
  fConeHEIDPhotonEndcap_isoCh.clear();
  fConeHEIDPhotonEndcap_HE.clear();
  fConeHEIDPhotonEndcap_coneHE.clear();
  fConeHEIDPhotonEndcap_sigmaIetaIeta.clear();
  
  fzDecayType = -10;
  fGenZ_pt = -100.0;
  fGenZ_eta = -100.0;
  fGenZ_phi = -100.0;
  fGenZ_mass = -100.0;
  fGenTau_pt.clear();
  fGenTau_eta.clear();
  fGenTau_phi.clear();
  fGenTau_mass.clear();
  fGenTau_objDR.clear();
  fGenTau_candobjDR.clear();
  fGenPhi_pt.clear();
  fGenPhi_eta.clear();
  fGenPhi_phi.clear();
  fGenPhi_mass.clear();
  fGenPhi_vx.clear();
  fGenPhi_vy.clear();
  fGenPhi_vz.clear();
  fGenPhi_vdiff_beamspot.clear();
  fGenPhi_vdiff_PV.clear();
  fGenOmega_pt.clear();
  fGenOmega_eta.clear();
  fGenOmega_phi.clear();
  fGenOmega_mass.clear();
  fGenOmega_vx.clear();
  fGenOmega_vy.clear();
  fGenOmega_vz.clear();
  fGenOmega_vdiff_beamspot.clear();
  fGenOmega_vdiff_PV.clear();
  fGenOmega_decayMode.clear();
  fGenOmega_neutral_pt.clear();
  fGenOmega_neutral_eta.clear();
  fGenOmega_neutral_phi.clear();
  fGenOmega_neutral_mass.clear();
  fGenOmega_positive_pt.clear();
  fGenOmega_positive_eta.clear();
  fGenOmega_positive_phi.clear();
  fGenOmega_positive_mass.clear();
  fGenOmega_negative_pt.clear();
  fGenOmega_negative_eta.clear();
  fGenOmega_negative_phi.clear();
  fGenOmega_negative_mass.clear();
  fGenOmega_posnegdr.clear();
  fGenOmega_objDR.clear();
  fGenOmega_candobjDR.clear();
  fGenOmega_jetDR.clear();
  fGenOmega_pattauDR.clear();
  fGenOmega_pattau1DR.clear();
  fGenOmega_pattau2DR.clear();
  fGenOmega_pattau3DR.clear();

  fTightMuon_pt.clear();
  fTightMuon_eta.clear();
  fTightMuon_phi.clear();
  fTightMuon_mass.clear();
  fW_pt.clear();
  fW_eta.clear();
  fW_phi.clear();
  fW_mass.clear();
  fW_mT.clear();
  fMuon_veto = -1;
  fbtag_veto = -1;
  fImag_W = -1;
  fmT = -1.0;
  fNum_Muons = -1;

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);

  // trigger 
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects);
  iEvent.getByToken(triggerPrescales_, triggerPrescales);

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  string trigger_photon175  = "HLT_Photon175_v"; // 2016
  string trigger_photon200  = "HLT_Photon200_v"; // 2017, 2018
  string trigger_photon22_iso  = "HLT_Photon22_R9Id90_HE10_IsoM_v";
  string trigger_photon165_iso  = "HLT_Photon165_R9Id90_HE10_IsoM_v"; // 2016, 2017
  string trigger_photon110_iso  = "HLT_Photon110EB_TightID_TightIso_v"; // 2018
  string trigger_photon75_vbf  = "HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF_v"; // 2016
  string trigger_doublephoton60 = "HLT_DoublePhoton60_v"; // 2016, 2017
  string trigger_doublephoton70 = "HLT_DoublePhoton70_v"; // 2018
  string trigger_doublephoton85 = "HLT_DoublePhoton85_v"; // 2016, 2017, 2018
  string trigger_tau120 = "HLT_VLooseIsoPFTau120_Trk50_eta2p1_v"; // 2016
  string trigger_doubletau32 = "HLT_DoubleMediumIsoPFTau32_Trk1_eta2p1_Reg_v"; // 2016
  bool bit_photon175 = false;
  bool bit_photon200 = false;
  bool bit_photon22_iso = false;
  bool bit_photon165_iso = false;
  bool bit_photon110_iso = false;
  bool bit_photon75_vbf = false;
  bool bit_doublephoton60 = false;
  bool bit_doublephoton70 = false;
  bool bit_doublephoton85 = false;
  bool bit_tau120 = false;
  bool bit_doubletau32 = false;
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i)
  {
     string triggerName = names.triggerName(i);
     
     std::size_t pos_175 = triggerName.find(trigger_photon175);
     if ( pos_175 != std::string::npos ) {
       bit_photon175 = triggerBits->accept(i);
       fPhotonFoundTrigger = triggerName;
     }

     std::size_t pos_photon200 = triggerName.find(trigger_photon200);
     if ( pos_photon200 != std::string::npos ) { bit_photon200 = triggerBits->accept(i); }

     std::size_t pos_photon22_iso = triggerName.find(trigger_photon22_iso);
     if ( pos_photon22_iso != std::string::npos ) { bit_photon22_iso = triggerBits->accept(i); }

     std::size_t pos_photon165_iso = triggerName.find(trigger_photon165_iso);
     if ( pos_photon165_iso != std::string::npos ) { bit_photon165_iso = triggerBits->accept(i); }

     std::size_t pos_photon110_iso = triggerName.find(trigger_photon110_iso);
     if ( pos_photon110_iso != std::string::npos ) { bit_photon110_iso = triggerBits->accept(i); }

     std::size_t pos_photon75_vbf = triggerName.find(trigger_photon75_vbf);
     if ( pos_photon75_vbf != std::string::npos ) { bit_photon75_vbf = triggerBits->accept(i); }

     std::size_t pos_doublephoton60 = triggerName.find(trigger_doublephoton60);
     if ( pos_doublephoton60 != std::string::npos ) { bit_doublephoton60 = triggerBits->accept(i); }

     std::size_t pos_doublephoton70 = triggerName.find(trigger_doublephoton70);
     if ( pos_doublephoton70 != std::string::npos ) { bit_doublephoton70 = triggerBits->accept(i); }

     std::size_t pos_doublephoton85 = triggerName.find(trigger_doublephoton85);
     if ( pos_doublephoton85 != std::string::npos ) { bit_doublephoton85 = triggerBits->accept(i); }

     std::size_t pos_tau120 = triggerName.find(trigger_tau120);
     if ( pos_tau120 != std::string::npos ) { bit_tau120 = triggerBits->accept(i); }

     std::size_t pos_doubletau32 = triggerName.find(trigger_doubletau32);
     if ( pos_doubletau32 != std::string::npos ) { bit_doubletau32 = triggerBits->accept(i); }

  }
  fHLT_Photon175 = bit_photon175;
  fHLT_Photon200 = bit_photon200;
  fHLT_Photon22_Iso = bit_photon22_iso;
  fHLT_Photon165_Iso = bit_photon165_iso;
  fHLT_Photon110_Iso = bit_photon110_iso;
  fHLT_Photon75_VBF = bit_photon75_vbf;
  fHLT_DoublePhoton60 = bit_doublephoton60;
  fHLT_DoublePhoton70 = bit_doublephoton70;
  fHLT_DoublePhoton85 = bit_doublephoton85;
  fHLT_Tau120 = bit_tau120;
  fHLT_DoubleTau32 = bit_doubletau32;

  // ecal tool
  lazyTools_ = std::auto_ptr<noZS::EcalClusterLazyTools>( new noZS::EcalClusterLazyTools(iEvent, iSetup, recHitsEBToken, recHitsEEToken));   

  // get event products
  edm::Handle<EcalRecHitCollection> recHitsEB_h;
  iEvent.getByToken(recHitsEBToken, recHitsEB_h );
  const EcalRecHitCollection * recHitsEB = 0;
  if ( ! recHitsEB_h.isValid() ) {
    LogError("TwoProngAnalyzer") << " ECAL Barrel RecHit Collection not available !"; return;
  } else {
    recHitsEB = recHitsEB_h.product();
  }

  edm::Handle<EcalRecHitCollection> recHitsEE_h;
  iEvent.getByToken(recHitsEEToken, recHitsEE_h );
  const EcalRecHitCollection * recHitsEE = 0;
  if ( ! recHitsEE_h.isValid() ) {
    LogError("TwoProngAnalyzer") << " ECAL Endcap RecHit Collection not available !"; return;
  } else {
    recHitsEE = recHitsEE_h.product();
  }

  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);

  edm::ESHandle<EcalChannelStatus> chStatus;
  iSetup.get<EcalChannelStatusRcd>().get(chStatus);
  const EcalChannelStatus *ch_status = chStatus.product(); 

  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(pfcandsToken_, pfcands);

  edm::Handle<vector<reco::GenParticle>> genparticles;
  if (fincludeSignalGenParticles || fincludeZDecayGenParticles || fincludeMCInfo) {
    iEvent.getByToken(genToken_, genparticles);
  }

  edm::Handle<GenEventInfoProduct> genEventInfo;
  edm::Handle<vector<reco::GenJet>> genJets;
  if (fincludeMCInfo) {
    iEvent.getByToken(genEventInfoToken_, genEventInfo);
    iEvent.getByToken(genJetsToken_, genJets);
  }

  edm::Handle<vector<reco::Vertex>> primaryvertecies;
  iEvent.getByToken(pvToken_, primaryvertecies);
  const reco::Vertex & PV = (*primaryvertecies)[0];

  edm::Handle<reco::BeamSpot> beamspot;
  iEvent.getByToken(beamToken_, beamspot);

  fPV_x = PV.position().X();
  fPV_y = PV.position().Y();
  fPV_z = PV.position().Z();
  fBeamspot_x = beamspot->position().X();
  fBeamspot_y = beamspot->position().Y();
  fBeamspot_z = beamspot->position().Z();

  edm::Handle<std::vector<pat::Jet>> ak4jets;
  iEvent.getByToken(ak4Token_, ak4jets);

  edm::Handle<std::vector<pat::Photon>> miniaod_photons;
  iEvent.getByToken(photonToken_, miniaod_photons);

  edm::Handle<edm::View<pat::Photon> > ged_photons;
  iEvent.getByToken(gedphotonsToken_, ged_photons);

  edm::Handle<std::vector<pat::MET>> MET;
  iEvent.getByToken(metToken_, MET);
  pat::MET met = (*MET)[0];

  edm::Handle<std::vector<pat::Electron>> electrons;
  iEvent.getByToken(electronToken_, electrons);

  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonToken_, muons);

  edm::Handle<std::vector<pat::Tau>> taus;
  iEvent.getByToken(tauToken_, taus);

  edm::Handle<vector<PileupSummaryInfo>> pileupInfo;
  iEvent.getByToken(pileupInfoToken_, pileupInfo);

  if (fincludeMCInfo) {
    if (fDebug) cout << ". doing mc info" << endl;
    // MC weights
    fMcW = genEventInfo->weight();
    fMcWProd = genEventInfo->weightProduct();

    // pthat, pt of the leading quark (relevant for MC QCD dijet events)
    double pthat = -100.0;
    for (unsigned int i = 0; i < genparticles->size(); i++) {
      const reco::GenParticle &genparticle = (*genparticles)[i];
      int id = genparticle.pdgId();
      if (genparticle.status() == 23 && (id == 1 || id == 2 || id == 3 || id == 4 || id == 5 || id == 6 || id == 21)) {
        if (pthat < genparticle.pt()) {
          pthat = genparticle.pt();
        }
      }
    }
    fpthat = pthat;

    // HT_gen
    double ht_gen = 0;
    for (unsigned int i = 0; i < genJets->size(); i++) {
      const reco::GenJet &genjet = (*genJets)[i];
      if (genjet.pt() < 30) continue;
      if (fabs(genjet.eta()) > 2.5) continue;
      ht_gen += genjet.pt();
    }
    fHT_gen = ht_gen;

    // MC PU info
    float true_npu = -1;
    int obs_npu = -1;
    for (std::vector<PileupSummaryInfo>::const_iterator pvi = pileupInfo->begin(); pvi != pileupInfo->end(); ++pvi)
    {
      int bx = pvi->getBunchCrossing();
      if (bx == 0) {
        true_npu = pvi->getTrueNumInteractions();
        obs_npu = pvi->getPU_NumInteractions();
        break;
      }
    }
    ftrueNpu = true_npu;
    fobsNpu = obs_npu;
  } else {
    // use escape values for MC branches
    fMcW = -100.0;
    fMcWProd = -100.0;
    fpthat = -100.0;
    fHT_gen = -100.0;
    ftrueNpu = -100.0;
    fobsNpu = -100.0;
  }

  // Signal MC Generator information
  if (fincludeSignalGenParticles) {
    if (fDebug) cout << ". doing signal gen particles" << endl; 
    for (unsigned int i = 0; i < genparticles->size(); i++) {
      const reco::GenParticle &genparticle = (*genparticles)[i];
      if ((genparticle.pdgId() != 54 && genparticle.pdgId() != 9000006) || genparticle.status() != 62) continue;
      TLorentzVector resonance;
      resonance.SetPtEtaPhiM(genparticle.pt(),genparticle.eta(),genparticle.phi(),genparticle.mass());
      fGenPhi_pt.push_back(resonance.Pt());
      fGenPhi_eta.push_back(resonance.Eta());
      if (fabs((resonance).CosTheta()) == 1) fGenPhi_eta.push_back(1000.0);
      else fGenPhi_eta.push_back(resonance.Eta());
      fGenPhi_phi.push_back(resonance.Phi());
      fGenPhi_mass.push_back(resonance.M());
      double vx = genparticle.vx();
      double vy = genparticle.vy();
      double vz = genparticle.vz();
      double pvx = PV.position().X();
      double pvy = PV.position().Y();
      double bvx = beamspot->position().X();
      double bvy = beamspot->position().Y();
      double vdiff_beamspot = std::sqrt((vx-bvx)*(vx-bvx) + (vy-bvy)*(vy-bvy));
      double vdiff_PV = std::sqrt((vx-pvx)*(vx-pvx) + (vy-pvy)*(vy-pvy));
      fGenPhi_vx.push_back(vx);
      fGenPhi_vy.push_back(vy);
      fGenPhi_vz.push_back(vz);
      fGenPhi_vdiff_beamspot.push_back(vdiff_beamspot);
      fGenPhi_vdiff_PV.push_back(vdiff_PV);
      for (unsigned int j = 0; j < genparticle.numberOfDaughters(); j++) {
        const reco::GenParticle* genparticle2 = (const reco::GenParticle*) genparticle.daughter(j);
        if (!fOldData && genparticle2->pdgId() != 90000054 && genparticle2->pdgId() != 90000055) continue;
        TLorentzVector pseudoscalar;
        pseudoscalar.SetPtEtaPhiE(genparticle2->pt(),genparticle2->eta(),genparticle2->phi(),genparticle2->energy());
        TLorentzVector positivePion;
        TLorentzVector negativePion;
        TLorentzVector neutralContent;
        neutralContent.SetXYZT(0,0,0,0);
        int n_pos = 0;
        int n_neg = 0;
        int n_pi0 = 0;
        int n_gam = 0;
        int n_rho = 0;
        int n_ome = 0;
        int n_eta = 0;
        int n2_pos = 0;
        int n2_neg = 0;
        int n2_pi0 = 0;
        int n2_gam = 0;
        int n2_rho = 0;
        int n2_ome = 0;
        int n2_eta = 0;
        int decayMode = -1;
        if (fDebug) cout << ".. now looping on daughters of pomega, looking at daughter # " << j+1 << endl; 
        for (unsigned int jj = 0; jj < genparticle2->numberOfDaughters(); jj++) {
          const reco::Candidate* genparticle3 = genparticle2->daughter(jj);
          if (fDebug) cout << ".. got pdgId " << genparticle3->pdgId() << endl; 
          TLorentzVector genparticle3Vect;
          genparticle3Vect.SetPtEtaPhiE(genparticle3->pt(), genparticle3->eta(), genparticle3->phi(), genparticle3->energy());
          // count for decay mode
          if(genparticle3->pdgId()==211)  n_pos += 1;
          if(genparticle3->pdgId()==-211) n_neg += 1;
          if(genparticle3->pdgId()==111)  n_pi0 += 1;
          if(genparticle3->pdgId()==22)   n_gam += 1;
          if(genparticle3->pdgId()==113)  n_rho += 1;
          if(genparticle3->pdgId()==223)  n_ome += 1;
          if(genparticle3->pdgId()==221)  n_eta += 1;
          // add to four-vectors
          if(genparticle3->pdgId()==111 || genparticle3->pdgId()==22) neutralContent += genparticle3Vect;
          if(genparticle3->pdgId()==211) positivePion.SetPtEtaPhiE(genparticle3->pt(),genparticle3->eta(),genparticle3->phi(),genparticle3->energy());
          if(genparticle3->pdgId()==-211) negativePion.SetPtEtaPhiE(genparticle3->pt(),genparticle3->eta(),genparticle3->phi(),genparticle3->energy());
          for (unsigned int jjj = 0; jjj < genparticle3->numberOfDaughters(); jjj++) {
            const reco::Candidate* genparticle4 = genparticle3->daughter(jjj);
            if(!(genparticle3->pdgId()==221 || genparticle3->pdgId()==113 || genparticle3->pdgId()==223)) continue; // only decencd to decay of eta, omega, or rho
            TLorentzVector genparticle4Vect;
            genparticle4Vect.SetPtEtaPhiE(genparticle4->pt(), genparticle4->eta(), genparticle4->phi(), genparticle4->energy());
            // count for decay mode
            if(genparticle4->pdgId()==211)  n2_pos += 1;
            if(genparticle4->pdgId()==-211) n2_neg += 1;
            if(genparticle4->pdgId()==111)  n2_pi0 += 1;
            if(genparticle4->pdgId()==22)   n2_gam += 1;
            if(genparticle4->pdgId()==113)  n2_rho += 1;
            if(genparticle4->pdgId()==223)  n2_ome += 1;
            if(genparticle4->pdgId()==221)  n2_eta += 1;
            // add to four-vectors
            if(genparticle4->pdgId()==111 || genparticle4->pdgId()==22) neutralContent += genparticle4Vect;
            if(genparticle4->pdgId()==211) positivePion.SetPtEtaPhiE(genparticle4->pt(),genparticle4->eta(),genparticle4->phi(),genparticle4->energy());
            if(genparticle4->pdgId()==-211) negativePion.SetPtEtaPhiE(genparticle4->pt(),genparticle4->eta(),genparticle4->phi(),genparticle4->energy());
          }
        }
        if (n_gam==2) decayMode = 1;
        if (n_pi0==3) decayMode = 2;
        if (n_pos==1 && n_neg==1 && n_pi0==1) decayMode = 3;
        if (n_pos==1 && n_neg==1 && n_gam==1) decayMode = 4;
        if (n_pi0==2 && n_eta==1 && n2_gam==2) decayMode = 5;
        if (n_pi0==2 && n_eta==1 && n2_pi0==3) decayMode = 6;
        if (n_pi0==2 && n_eta==1 && n2_gam==2) decayMode = 7;
        //if (n_gam==2) decayMode = 8; // will get folded away into 1
        if (n_pos==1 && n_neg==1 && n_eta==1 && n2_gam==2) decayMode = 9;
        if (n_pos==1 && n_neg==1 && n_eta==1 && n2_pi0==3) decayMode = 10;
        if (n_pi0==2 && n_eta==1 && n2_pos==1 && n2_neg==1 && n2_pi0==1) decayMode = 11;
        if (n_pi0==2 && n_eta==1 && n2_pos==1 && n2_neg==1 && n2_gam==1) decayMode = 12;
        if (n_rho==1 && n_gam==1 && n2_pos==1 && n2_neg==1) decayMode = 13;
        if (n_ome==1 && n_gam==1 && n2_pos==1 && n2_neg==1 && n2_pi0==1) decayMode = 14;
        if (n_pos==1 && n_neg==1 && n_eta==1 && n2_pos==1 && n2_neg==1 && n2_pi0==1) decayMode = 15;
        if (n_pos==1 && n_neg==1 && n_eta==1 && n2_pos==1 && n2_neg==1 && n2_gam==1) decayMode = 16;
        if (decayMode==-1 && fDebug) cout << ".. decay mode finding failed!:" << " " << n_pos << " " << n_neg << " " << n_pi0 << " " << n_gam << " " << n_rho << " " << n_ome << " " << n_eta << endl;
        fGenOmega_pt.push_back(pseudoscalar.Pt());
        fGenOmega_eta.push_back(pseudoscalar.Eta());
        fGenOmega_phi.push_back(pseudoscalar.Phi());
        fGenOmega_mass.push_back(pseudoscalar.M());
        vx = genparticle2->vx();
        vy = genparticle2->vy();
        vz = genparticle2->vz();
        vdiff_beamspot = std::sqrt((vx-bvx)*(vx-bvx) + (vy-bvy)*(vy-bvy));
        vdiff_PV = std::sqrt((vx-pvx)*(vx-pvx) + (vy-pvy)*(vy-pvy));
        fGenOmega_vx.push_back(vx);
        fGenOmega_vy.push_back(vy);
        fGenOmega_vz.push_back(vz);
        fGenOmega_vdiff_beamspot.push_back(vdiff_beamspot);
        fGenOmega_vdiff_PV.push_back(vdiff_PV);
        fGenOmega_decayMode.push_back(decayMode);
        fGenOmega_neutral_pt.push_back(neutralContent.Pt());
        fGenOmega_neutral_eta.push_back(neutralContent.Eta());
        fGenOmega_neutral_phi.push_back(neutralContent.Phi());
        fGenOmega_neutral_mass.push_back(neutralContent.M());
        fGenOmega_positive_pt.push_back(positivePion.Pt());
        fGenOmega_positive_eta.push_back(positivePion.Eta());
        fGenOmega_positive_phi.push_back(positivePion.Phi());
        fGenOmega_positive_mass.push_back(positivePion.M());
        fGenOmega_negative_pt.push_back(negativePion.Pt());
        fGenOmega_negative_eta.push_back(negativePion.Eta());
        fGenOmega_negative_phi.push_back(negativePion.Phi());
        fGenOmega_negative_mass.push_back(negativePion.M());
        fGenOmega_posnegdr.push_back(positivePion.DeltaR(negativePion));
      } // end loop on daughters of Phi
    } // end loop on all gen particles
  }

  // Jets
  for (unsigned int i = 0; i < ak4jets->size(); i++) {
    const pat::Jet &jet = (*ak4jets)[i];
    if (jet.pt() < 30) continue;
    if (fabs(jet.eta()) > 2.5) continue;
    fAK4jet_pt.push_back(jet.pt());
    fAK4jet_eta.push_back(jet.eta());
    fAK4jet_phi.push_back(jet.phi());
    fAK4jet_mass.push_back(jet.mass());
  }
  fNumAK4jets = fAK4jet_pt.size();

  // Photons
  for (unsigned int i = 0; i < miniaod_photons->size(); i++) {
    const pat::Photon &photon = (*miniaod_photons)[i];
    fPhoton_pt.push_back(photon.pt());
    fPhoton_eta.push_back(photon.eta());
    fPhoton_scEta.push_back( photon_scEta(&photon) );
    fPhoton_phi.push_back(photon.phi());
    fPhoton_mass.push_back(photon.mass());
    fPhoton_isoGamma.push_back( photon_computeIsoGamma(&photon, *rhoH) );
    fPhoton_isoCh.push_back( photon_computeIsoCh(&photon) );
    fPhoton_HE.push_back( photon_computeHE(&photon) );
    fPhoton_coneHE.push_back( photon_computeHE_coneBased(&photon) );
    fPhoton_sigmaIetaIeta.push_back( photon_computeSigmaIetaIeta(&photon) );
    fPhoton_passVeto.push_back( photon.passElectronVeto() );
  }
  fNumPhotons = miniaod_photons->size();

  // Electrons
  for (unsigned int i = 0; i < electrons->size(); i++) {
    const pat::Electron &electron = (*electrons)[i];
    fElectron_pt.push_back(electron.pt());
    fElectron_eta.push_back(electron.eta());
    fElectron_phi.push_back(electron.phi());
    fElectron_mass.push_back(electron.mass());
  }
  fNumElectrons = electrons->size();

  // Muons
  for (unsigned int i = 0; i < muons->size(); i++) {
    const pat::Muon &muon = (*muons)[i];
    fMuon_pt.push_back(muon.pt());
    fMuon_eta.push_back(muon.eta());
    fMuon_phi.push_back(muon.phi());
    fMuon_mass.push_back(muon.mass());
  }
  fNumMuons = muons->size();

  // Taus
  for (unsigned int i = 0; i < taus->size(); i++) {
    const pat::Tau &tau = (*taus)[i];
    if (tau.pt() <= 20) continue;
    if (fabs(tau.eta()) >= 2.3) continue;
    if (tau.tauID("againstElectronVLooseMVA6") < 0.5) continue;
    if (tau.tauID("againstMuonTight3") < 0.5) continue;
    if (tau.tauID("byTightIsolationMVArun2v1DBnewDMwLT") < 0.5) continue;
    fTau_pt.push_back(tau.pt());
    fTau_eta.push_back(tau.eta());
    fTau_phi.push_back(tau.phi());
    fTau_mass.push_back(tau.mass());
  }
  fNumTaus = fTau_pt.size(); 

  // Event wide information
  fEventNum = iEvent.id().event();
  fRunNum = iEvent.id().run();
  fLumiNum = iEvent.id().luminosityBlock();
  fMET = (*MET)[0].pt();
  fMET_phi = (*MET)[0].phi();
  fNumPVs = primaryvertecies->size();
  fRho = *rhoH;
  fNumPF = pfcands->size();

  // Two prongs
  if (fDebug) cout << ". starting two prong code" << endl;
  // Find all pairs of one CH pos and one CH neg within specified DR of each other
  for (unsigned int i = 0; i < pfcands->size(); i++) {
    const pat::PackedCandidate &pf1 = (*pfcands)[i];
    if (pf1.pt() < ftwoprong_tracksMinPt) continue;
    if (pf1.fromPV()<=1) continue;
    for (unsigned int j = i+1; j < pfcands->size(); j++) { // note loop starting with j=i+1, considers each pair exactly once
      const pat::PackedCandidate &pf2 = (*pfcands)[j];
      if (pf2.pt() < ftwoprong_tracksMinPt) continue;
      if (pf2.fromPV()<=1) continue;
      if (!( ((pf1.pdgId() == 211) && (pf2.pdgId() == -211)) || ((pf1.pdgId() == -211) && (pf2.pdgId() == 211)) )) continue;
      TLorentzVector pfcand1;
      pfcand1.SetPtEtaPhiE(pf1.pt(), pf1.eta(), pf1.phiAtVtx(), pf1.energy());
      TLorentzVector pfcand2;
      pfcand2.SetPtEtaPhiE(pf2.pt(), pf2.eta(), pf2.phiAtVtx(), pf2.energy());
      // found CH pos and CH minus, now check DR
      double dr = pfcand1.DeltaR(pfcand2);
      if (dr < ftwoprong_DR) {
        // search for extra tracks
        vector<unsigned int> index_of_extras;
        for (unsigned int e = 0; e < pfcands->size(); e++) {
          if (e == i || e == j) continue;
          const pat::PackedCandidate &pfextra = (*pfcands)[e];
          if (pfextra.pt() < ftwoprong_tracksMinPt) continue;
          if (pfextra.fromPV()<=1) continue;
          if (abs(pfextra.pdgId()) != 211) continue;
          TLorentzVector pfcandextra;
          pfcandextra.SetPtEtaPhiE(pfextra.pt(), pfextra.eta(), pfextra.phiAtVtx(), pfextra.energy());
          double dr = max( pfcand1.DeltaR(pfcandextra), pfcand2.DeltaR(pfcandextra) );
          if (dr > ftwoprong_DR) continue;
          index_of_extras.push_back(e);
        }
        unsigned int index_of_extra = 99999;
        bool found_extra_track = false;
        TLorentzVector pfcandextra;       
        if (index_of_extras.size() == 1) {
          index_of_extra = index_of_extras[0];
          const pat::PackedCandidate &pfextra = (*pfcands)[index_of_extra];
          pfcandextra.SetPtEtaPhiE(pfextra.pt(), pfextra.eta(), pfextra.phiAtVtx(), pfextra.energy());
          found_extra_track = true;
        }
        if (index_of_extras.size() >= 2) {
          double largest_pt_extra = -1.0;
          for (unsigned int index : index_of_extras) {
            double pt_of_extra = ((*pfcands)[index]).pt();
            if (pt_of_extra > largest_pt_extra) {
              largest_pt_extra = pt_of_extra;
              index_of_extra = index; 
            }
          }
          const pat::PackedCandidate &pfextra = (*pfcands)[index_of_extra];
          pfcandextra.SetPtEtaPhiE(pfextra.pt(), pfextra.eta(), pfextra.phiAtVtx(), pfextra.energy());
          found_extra_track = true;
        }
        // now define photon
        TLorentzVector center;
        if (ftwoprong_OptionalExtraTrack && found_extra_track) {
          center = pfcand1 + pfcand2; }
        else {
          center = pfcand1 + pfcand2; }
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
          if (fabs(pf3.phiAtVtx() - center.Phi()) < ftwoprong_PhiBox/2.0 &&
              fabs(pf3.eta() - center.Eta()) < ftwoprong_EtaBox/2.0) {
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
        int n = index_of_leading_pf_photon;
        if (n != -1) {
          TLorentzVector TwoProngObject;
          TwoProngObject = center + photon;
          if (ftwoprong_OptionalExtraTrack && found_extra_track) {
            TwoProngObject = center + photon + pfcandextra;
          }

          leading_pf_photon.SetPtEtaPhiE((*pfcands)[n].pt(), (*pfcands)[n].eta(), (*pfcands)[n].phiAtVtx(), (*pfcands)[n].energy());

          TLorentzVector neutral_as_zeromass; neutral_as_zeromass.SetPtEtaPhiM(photon.Pt(), photon.Eta(), photon.Phi(), 0.0);
          TLorentzVector neutral_as_pi0mass; neutral_as_pi0mass.SetPtEtaPhiM(photon.Pt(), photon.Eta(), photon.Phi(), PI0_MASS);
          TLorentzVector neutral_as_etamass; neutral_as_etamass.SetPtEtaPhiM(photon.Pt(), photon.Eta(), photon.Phi(), ETA_MASS);
          TLorentzVector neutral_as_300mass; neutral_as_300mass.SetPtEtaPhiM(photon.Pt(), photon.Eta(), photon.Phi(), 0.3);

          if (fabs(TwoProngObject.Eta()) > ftwoprong_AbsMaxEta) continue;
          if (TwoProngObject.Pt() < ftwoprong_MinPt) continue;

          // CH pair within dr range
          // has at least one pf photon
          // |eta| <= eta requirement
          // pt >= min pt requirement
          //   meets definition of candidate twoprong, fill candidate vectors
          TLorentzVector poscand;
          TLorentzVector negcand;
          if (pf1.pdgId() > 0) {
            poscand = pfcand1;
            negcand = pfcand2;
            fTwoProngCand_CHpos_pt.push_back(pfcand1.Pt());
            fTwoProngCand_CHpos_eta.push_back(pfcand1.Eta());
            fTwoProngCand_CHpos_phi.push_back(pfcand1.Phi());
            fTwoProngCand_CHpos_mass.push_back(pfcand1.M());
            fTwoProngCand_CHneg_pt.push_back(pfcand2.Pt());
            fTwoProngCand_CHneg_eta.push_back(pfcand2.Eta());
            fTwoProngCand_CHneg_phi.push_back(pfcand2.Phi());
            fTwoProngCand_CHneg_mass.push_back(pfcand2.M());
            // impact parameters 
            fTwoProngCand_CHpos_dz.push_back( pf1.dz(PV.position()) );
            fTwoProngCand_CHpos_dxy.push_back( pf1.dxy(PV.position()));
            fTwoProngCand_CHneg_dz.push_back( pf2.dz( PV.position()) );
            fTwoProngCand_CHneg_dxy.push_back( pf2.dxy( PV.position()) );
          } else {
            poscand = pfcand2;
            negcand = pfcand1;
            fTwoProngCand_CHpos_pt.push_back(pfcand2.Pt());
            fTwoProngCand_CHpos_eta.push_back(pfcand2.Eta());
            fTwoProngCand_CHpos_phi.push_back(pfcand2.Phi());
            fTwoProngCand_CHpos_mass.push_back(pfcand2.M());
            fTwoProngCand_CHneg_pt.push_back(pfcand1.Pt());
            fTwoProngCand_CHneg_eta.push_back(pfcand1.Eta());
            fTwoProngCand_CHneg_phi.push_back(pfcand1.Phi());
            fTwoProngCand_CHneg_mass.push_back(pfcand1.M());
            // impact parameters
            fTwoProngCand_CHneg_dz.push_back( pf1.dz( PV.position()) );
            fTwoProngCand_CHneg_dxy.push_back( pf1.dxy( PV.position()) );
            fTwoProngCand_CHpos_dz.push_back( pf2.dz( PV.position()) );
            fTwoProngCand_CHpos_dxy.push_back( pf2.dxy( PV.position()) );
          }
          fTwoProngCand_photon_pt.push_back(photon.Pt());
          fTwoProngCand_photon_eta.push_back(photon.Eta());
          fTwoProngCand_photon_phi.push_back(photon.Phi());
          fTwoProngCand_photon_mass.push_back(photon.M());
          fTwoProngCand_photon_pt_l.push_back(leading_pf_photon.Pt());
          fTwoProngCand_photon_eta_l.push_back(leading_pf_photon.Eta());
          fTwoProngCand_photon_phi_l.push_back(leading_pf_photon.Phi());
          fTwoProngCand_photon_mass_l.push_back(leading_pf_photon.M());
          fTwoProngCand_photon_nGamma.push_back(numgamma);
          fTwoProngCand_photon_nElectron.push_back(nume);
          fTwoProngCand_pt.push_back(TwoProngObject.Pt());
          fTwoProngCand_eta.push_back(TwoProngObject.Eta());
          fTwoProngCand_phi.push_back(TwoProngObject.Phi());

          fTwoProngCand_mass.push_back(    TwoProngObject.M());
          fTwoProngCand_mass_l.push_back(  (center + leading_pf_photon).M() );
          fTwoProngCand_Mass0.push_back(   (center + neutral_as_zeromass).M() );
          fTwoProngCand_MassPi0.push_back( (center + neutral_as_pi0mass).M() );
          fTwoProngCand_MassEta.push_back( (center + neutral_as_etamass).M() );
          fTwoProngCand_Mass300.push_back( (center + neutral_as_300mass).M() );

          fTwoProngCand_mPosPho.push_back( (poscand + photon).M2() );
          fTwoProngCand_mPosPho_l.push_back( (poscand + leading_pf_photon).M2() );
          fTwoProngCand_mPosPho_pi0.push_back( (poscand + neutral_as_pi0mass).M2() );
          fTwoProngCand_mNegPho.push_back( (negcand + photon).M2() );
          fTwoProngCand_mNegPho_l.push_back( (negcand + leading_pf_photon).M2() );
          fTwoProngCand_mNegPho_pi0.push_back( (negcand + neutral_as_pi0mass).M2() );
          fTwoProngCand_mPosNeg.push_back( (poscand + negcand).M2() );
          fTwoProngCand_CHpos_p3.push_back(poscand.P());
          fTwoProngCand_CHneg_p3.push_back(negcand.P());
          fTwoProngCand_photon_p3.push_back(photon.P());

          fTwoProngCand_nExtraTracks.push_back(index_of_extras.size());

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
              if ( center.DeltaR(pfcand4) < ftwoprong_IsolationDR && !(m == i || m == j) ) { // don't include one of CH from CH pair 
                if ((ftwoprong_OptionalExtraTrack && found_extra_track && m != index_of_extra) || !ftwoprong_OptionalExtraTrack  || !found_extra_track) { // if including an extra track, skip its iso
                  chargedIso += pfcand4.Pt();
                  chargedIsoCount++;
                }
              }
            // neutral
            } else if (pf4.pdgId() == 130) {
              if (center.DeltaR(pfcand4) < ftwoprong_IsolationDR) {
                neutralIso += pfcand4.Pt();
                  neutralIsoCount++;
              }
            // e gamma
            } else if (abs(pf4.pdgId()) == 11 || pf4.pdgId() == 22) {
              if ( (center.DeltaR(pfcand4) < ftwoprong_IsolationDR) &&
                   !(fabs(pf4.phiAtVtx() - center.Phi()) < ftwoprong_PhiBox/2.0 && fabs(pf4.eta() - center.Eta()) < ftwoprong_EtaBox/2.0)) {
                egammaIso += pfcand4.Pt();
                  egammaIsoCount++;
              }
            }
          } // end pf cand loop
          double relchargedIso = chargedIso / TwoProngObject.Pt();
          double relneutralIso = neutralIso / TwoProngObject.Pt();
          double relegammaIso = egammaIso / TwoProngObject.Pt();
          fTwoProngCand_chargedIso.push_back(relchargedIso);
          fTwoProngCand_neutralIso.push_back(relneutralIso);
          fTwoProngCand_egammaIso.push_back(relegammaIso);
          fTwoProngCand_nChargedIsoCone.push_back(chargedIsoCount);
          fTwoProngCand_nNeutralIsoCone.push_back(neutralIsoCount);
          fTwoProngCand_nEGammaIsoCone.push_back(egammaIsoCount);
          
          // optimizing iso
          //double etaForIso = TwoProngObject.Eta();
          //double isoGammaCor1 = egammaIso + iso_alpha(etaForIso, 1) - fRho * iso_area(etaForIso, 1) - iso_kappa(etaForIso, 1) * TwoProngObject.Pt();
          //double isoGammaCor2 = egammaIso + iso_alpha(etaForIso, 2) - fRho * iso_area(etaForIso, 2) - iso_kappa(etaForIso, 2) * TwoProngObject.Pt();

          // Asymmetry variables
          double track_asymmetry = 1.0;
          double photon_asymmetry = 1.0;
          if (!ftwoprong_OptionalExtraTrack || !found_extra_track) {
            track_asymmetry = min(pfcand1.Pt(),pfcand2.Pt()) / max(pfcand1.Pt(),pfcand2.Pt());
            photon_asymmetry = min(pfcand1.Pt()+pfcand2.Pt(),photon.Pt()) / max(pfcand1.Pt()+pfcand2.Pt(),photon.Pt());
          }
          else {
            track_asymmetry = min( min(pfcand1.Pt(),pfcand2.Pt()), pfcandextra.Pt()) / max( max(pfcand1.Pt(),pfcand2.Pt()), pfcandextra.Pt() );
            photon_asymmetry = min(pfcand1.Pt()+pfcand2.Pt()+pfcandextra.Pt(),photon.Pt()) / max(pfcand1.Pt()+pfcand2.Pt()+pfcandextra.Pt(),photon.Pt());
          }
          bool passTrackAsymmetry = (track_asymmetry > ftwoprong_TrackAsymmetryCut);
          bool passPhotonAsymmetry = (photon_asymmetry > ftwoprong_PhotonAsymmetryCut);
          fTwoProngCand_trackAsym.push_back(track_asymmetry);
          fTwoProngCand_photonAsym.push_back(photon_asymmetry);
          if(ftwoprong_FlipAsymReq) {
            passTrackAsymmetry = !passTrackAsymmetry;
            passPhotonAsymmetry = !passPhotonAsymmetry;
          }

          // Selection on Candidates
          bool passCharged = relchargedIso < ftwoprong_ChargedIsoCut;
          bool passNeutral = relneutralIso < ftwoprong_NeutralIsoCut;
          bool passEGamma =  relegammaIso < ftwoprong_EGammaIsoCut;
          bool looseCharged = relchargedIso < ftwoprong_ChargedIsoFakeCut;
          bool looseNeutral = relneutralIso < ftwoprong_NeutralIsoFakeCut;
          bool looseEGamma =  relegammaIso  < ftwoprong_EGammaIsoFakeCut;
          bool passPhotonPt = photon.Pt() > ftwoprong_PhotonPtCut;

          bool tight = passCharged && passNeutral && passEGamma && passPhotonPt && passTrackAsymmetry && passPhotonAsymmetry;
          bool loose = !tight && looseCharged && looseNeutral && looseEGamma && passPhotonPt && passTrackAsymmetry && passPhotonAsymmetry;
          bool asym = passCharged && passNeutral && passEGamma && passPhotonPt && !(passTrackAsymmetry && passPhotonAsymmetry);
          bool asym_loose = !asym && looseCharged && looseNeutral && looseEGamma && passPhotonPt && !(passTrackAsymmetry && passPhotonAsymmetry);

          fTwoProngCand_tight.push_back(tight);
          fTwoProngCand_loose.push_back(loose);
          fTwoProngCand_asym.push_back(asym);
          fTwoProngCand_asym_loose.push_back(asym_loose);

          if (fDebug) cout << ". doing gen matching" << endl;
          // Generator matching to signal
          double genOmega_dR = 99.9;
	        if (fincludeSignalGenParticles) {
	          for (unsigned int i = 0; i < genparticles->size(); i++) {
	            const reco::GenParticle &genparticle = (*genparticles)[i];
	            if ((genparticle.pdgId() == 221) && genparticle.status() == 2) {
                TLorentzVector genEta;
                genEta.SetPtEtaPhiM(genparticle.pt(), genparticle.eta(), genparticle.phi(), genparticle.mass());
                double match_dR = genEta.DeltaR(TwoProngObject);
                if (match_dR < genOmega_dR) genOmega_dR = match_dR;
              }
            }
          }
          fTwoProngCand_genOmega_dR.push_back(genOmega_dR);
          // Generator matching to tau
          double genTau_dR = 99.9;
	        if (fincludeZDecayGenParticles) {
	          for (unsigned int i = 0; i < genparticles->size(); i++) {
	            const reco::GenParticle &genparticle = (*genparticles)[i];
	            if (abs(genparticle.pdgId()) == 15) {
                TLorentzVector genTau;
                genTau.SetPtEtaPhiM(genparticle.pt(), genparticle.eta(), genparticle.phi(), genparticle.mass());
                double match_dR = genTau.DeltaR(TwoProngObject);
                if (match_dR < genTau_dR) genTau_dR = match_dR;
              }
            }
          }
          fTwoProngCand_genTau_dR.push_back(genTau_dR);
        }
      } // end conditionals on CH pair
    }
  } // end making candidates
  // Create sorted-by-pt lists
  if (fDebug) cout << ". sorting" << endl;
  vector<unsigned int> sorted_indecies;
  for (unsigned int i = 0; i < fTwoProngCand_pt.size(); i++) {
    double largestPtSoFar = -1.0;
    unsigned int largestPtSoFarIndex = 999;
    for (unsigned int j = 0; j < fTwoProngCand_pt.size(); j++) {
      bool skip = false;
      for (unsigned int n = 0; n < sorted_indecies.size(); n++) {
        if (sorted_indecies[n] == j) {
          skip = true;
          break; } }
      if (skip) continue;
      else if (fTwoProngCand_pt[j] > largestPtSoFar) {
        largestPtSoFar = fTwoProngCand_pt[j];
        largestPtSoFarIndex = j; } }
    sorted_indecies.push_back(largestPtSoFarIndex);
  }
  for (unsigned int i = 0; i < sorted_indecies.size(); i++) {
    unsigned int index = sorted_indecies[i];
    if (fTwoProngCand_asym[index])
    {
      fTwoProngAsym_pt.push_back(fTwoProngCand_pt[index]);
      fTwoProngAsym_eta.push_back(fTwoProngCand_eta[index]);
      fTwoProngAsym_phi.push_back(fTwoProngCand_phi[index]);
      fTwoProngAsym_mass.push_back(fTwoProngCand_mass[index]);
      fTwoProngAsym_mass_l.push_back(fTwoProngCand_mass_l[index]);
      fTwoProngAsym_Mass0.push_back(fTwoProngCand_Mass0[index]);
      fTwoProngAsym_MassPi0.push_back(fTwoProngCand_MassPi0[index]);
      fTwoProngAsym_MassEta.push_back(fTwoProngCand_MassEta[index]);
      fTwoProngAsym_Mass300.push_back(fTwoProngCand_Mass300[index]);
      fTwoProngAsym_nExtraTracks.push_back(fTwoProngCand_nExtraTracks[index]);
      fTwoProngAsym_CHpos_pt.push_back(fTwoProngCand_CHpos_pt[index]);
      fTwoProngAsym_CHpos_eta.push_back(fTwoProngCand_CHpos_eta[index]);
      fTwoProngAsym_CHpos_phi.push_back(fTwoProngCand_CHpos_phi[index]);
      fTwoProngAsym_CHpos_mass.push_back(fTwoProngCand_CHpos_mass[index]);
      fTwoProngAsym_CHpos_dz.push_back(fTwoProngCand_CHpos_dz[index]);
      fTwoProngAsym_CHpos_dxy.push_back(fTwoProngCand_CHpos_dxy[index]);
      fTwoProngAsym_CHneg_pt.push_back(fTwoProngCand_CHneg_pt[index]);
      fTwoProngAsym_CHneg_eta.push_back(fTwoProngCand_CHneg_eta[index]);
      fTwoProngAsym_CHneg_phi.push_back(fTwoProngCand_CHneg_phi[index]);
      fTwoProngAsym_CHneg_mass.push_back(fTwoProngCand_CHneg_mass[index]);
      fTwoProngAsym_CHneg_dz.push_back(fTwoProngCand_CHneg_dz[index]);
      fTwoProngAsym_CHneg_dxy.push_back(fTwoProngCand_CHneg_dxy[index]);
      fTwoProngAsym_photon_pt.push_back(fTwoProngCand_photon_pt[index]);
      fTwoProngAsym_photon_eta.push_back(fTwoProngCand_photon_eta[index]);
      fTwoProngAsym_photon_phi.push_back(fTwoProngCand_photon_phi[index]);
      fTwoProngAsym_photon_mass.push_back(fTwoProngCand_photon_mass[index]);
      fTwoProngAsym_photon_pt_l.push_back(fTwoProngCand_photon_pt_l[index]);
      fTwoProngAsym_photon_eta_l.push_back(fTwoProngCand_photon_eta_l[index]);
      fTwoProngAsym_photon_phi_l.push_back(fTwoProngCand_photon_phi_l[index]);
      fTwoProngAsym_photon_mass_l.push_back(fTwoProngCand_photon_mass_l[index]);
      fTwoProngAsym_photon_nGamma.push_back(fTwoProngCand_photon_nGamma[index]);
      fTwoProngAsym_photon_nElectron.push_back(fTwoProngCand_photon_nElectron[index]);
      fTwoProngAsym_chargedIso.push_back(fTwoProngCand_chargedIso[index]);
      fTwoProngAsym_neutralIso.push_back(fTwoProngCand_neutralIso[index]);
      fTwoProngAsym_egammaIso.push_back(fTwoProngCand_egammaIso[index]);
      fTwoProngAsym_nChargedIsoCone.push_back(fTwoProngCand_nChargedIsoCone[index]);
      fTwoProngAsym_nNeutralIsoCone.push_back(fTwoProngCand_nNeutralIsoCone[index]);
      fTwoProngAsym_nEGammaIsoCone.push_back(fTwoProngCand_nChargedIsoCone[index]);
      fTwoProngAsym_genOmega_dR.push_back(fTwoProngCand_genOmega_dR[index]);
      fTwoProngAsym_genTau_dR.push_back(fTwoProngCand_genTau_dR[index]);
      fTwoProngAsym_trackAsym.push_back(fTwoProngCand_trackAsym[index]);
      fTwoProngAsym_photonAsym.push_back(fTwoProngCand_photonAsym[index]);
      fTwoProngAsym_mPosPho_l.push_back(fTwoProngCand_mPosPho_l[index]);
      fTwoProngAsym_mPosPho_pi0.push_back(fTwoProngCand_mPosPho_pi0[index]);
      fTwoProngAsym_mNegPho.push_back(fTwoProngCand_mNegPho[index]);
      fTwoProngAsym_mNegPho_l.push_back(fTwoProngCand_mNegPho_l[index]);
      fTwoProngAsym_mNegPho_pi0.push_back(fTwoProngCand_mNegPho_pi0[index]);
      fTwoProngAsym_mPosNeg.push_back(fTwoProngCand_mPosNeg[index]);
    }
    if (fTwoProngCand_asym_loose[index])
    {
      fTwoProngAsymLoose_pt.push_back(fTwoProngCand_pt[index]);
      fTwoProngAsymLoose_eta.push_back(fTwoProngCand_eta[index]);
      fTwoProngAsymLoose_phi.push_back(fTwoProngCand_phi[index]);
      fTwoProngAsymLoose_mass.push_back(fTwoProngCand_mass[index]);
      fTwoProngAsymLoose_mass_l.push_back(fTwoProngCand_mass_l[index]);
      fTwoProngAsymLoose_Mass0.push_back(fTwoProngCand_Mass0[index]);
      fTwoProngAsymLoose_MassPi0.push_back(fTwoProngCand_MassPi0[index]);
      fTwoProngAsymLoose_MassEta.push_back(fTwoProngCand_MassEta[index]);
      fTwoProngAsymLoose_Mass300.push_back(fTwoProngCand_Mass300[index]);
      fTwoProngAsymLoose_nExtraTracks.push_back(fTwoProngCand_nExtraTracks[index]);
      fTwoProngAsymLoose_CHpos_pt.push_back(fTwoProngCand_CHpos_pt[index]);
      fTwoProngAsymLoose_CHpos_eta.push_back(fTwoProngCand_CHpos_eta[index]);
      fTwoProngAsymLoose_CHpos_phi.push_back(fTwoProngCand_CHpos_phi[index]);
      fTwoProngAsymLoose_CHpos_mass.push_back(fTwoProngCand_CHpos_mass[index]);
      fTwoProngAsymLoose_CHpos_dz.push_back(fTwoProngCand_CHpos_dz[index]);
      fTwoProngAsymLoose_CHpos_dxy.push_back(fTwoProngCand_CHpos_dxy[index]);
      fTwoProngAsymLoose_CHneg_pt.push_back(fTwoProngCand_CHneg_pt[index]);
      fTwoProngAsymLoose_CHneg_eta.push_back(fTwoProngCand_CHneg_eta[index]);
      fTwoProngAsymLoose_CHneg_phi.push_back(fTwoProngCand_CHneg_phi[index]);
      fTwoProngAsymLoose_CHneg_mass.push_back(fTwoProngCand_CHneg_mass[index]);
      fTwoProngAsymLoose_CHneg_dz.push_back(fTwoProngCand_CHneg_dz[index]);
      fTwoProngAsymLoose_CHneg_dxy.push_back(fTwoProngCand_CHneg_dxy[index]);
      fTwoProngAsymLoose_photon_pt.push_back(fTwoProngCand_photon_pt[index]);
      fTwoProngAsymLoose_photon_eta.push_back(fTwoProngCand_photon_eta[index]);
      fTwoProngAsymLoose_photon_phi.push_back(fTwoProngCand_photon_phi[index]);
      fTwoProngAsymLoose_photon_mass.push_back(fTwoProngCand_photon_mass[index]);
      fTwoProngAsymLoose_photon_pt_l.push_back(fTwoProngCand_photon_pt_l[index]);
      fTwoProngAsymLoose_photon_eta_l.push_back(fTwoProngCand_photon_eta_l[index]);
      fTwoProngAsymLoose_photon_phi_l.push_back(fTwoProngCand_photon_phi_l[index]);
      fTwoProngAsymLoose_photon_mass_l.push_back(fTwoProngCand_photon_mass_l[index]);
      fTwoProngAsymLoose_photon_nGamma.push_back(fTwoProngCand_photon_nGamma[index]);
      fTwoProngAsymLoose_photon_nElectron.push_back(fTwoProngCand_photon_nElectron[index]);
      fTwoProngAsymLoose_chargedIso.push_back(fTwoProngCand_chargedIso[index]);
      fTwoProngAsymLoose_neutralIso.push_back(fTwoProngCand_neutralIso[index]);
      fTwoProngAsymLoose_egammaIso.push_back(fTwoProngCand_egammaIso[index]);
      fTwoProngAsymLoose_nChargedIsoCone.push_back(fTwoProngCand_nChargedIsoCone[index]);
      fTwoProngAsymLoose_nNeutralIsoCone.push_back(fTwoProngCand_nNeutralIsoCone[index]);
      fTwoProngAsymLoose_nEGammaIsoCone.push_back(fTwoProngCand_nChargedIsoCone[index]);
      fTwoProngAsymLoose_genOmega_dR.push_back(fTwoProngCand_genOmega_dR[index]);
      fTwoProngAsymLoose_genTau_dR.push_back(fTwoProngCand_genTau_dR[index]);
      fTwoProngAsymLoose_trackAsym.push_back(fTwoProngCand_trackAsym[index]);
      fTwoProngAsymLoose_photonAsym.push_back(fTwoProngCand_photonAsym[index]);
      fTwoProngAsymLoose_mPosPho_l.push_back(fTwoProngCand_mPosPho_l[index]);
      fTwoProngAsymLoose_mPosPho_pi0.push_back(fTwoProngCand_mPosPho_pi0[index]);
      fTwoProngAsymLoose_mNegPho.push_back(fTwoProngCand_mNegPho[index]);
      fTwoProngAsymLoose_mNegPho_l.push_back(fTwoProngCand_mNegPho_l[index]);
      fTwoProngAsymLoose_mNegPho_pi0.push_back(fTwoProngCand_mNegPho_pi0[index]);
      fTwoProngAsymLoose_mPosNeg.push_back(fTwoProngCand_mPosNeg[index]);
    }
    if (fTwoProngCand_loose[index])
    {
      // Candidate is loose and is next leading, fill all loose candidate collections
      fTwoProngLoose_pt.push_back(fTwoProngCand_pt[index]);
      fTwoProngLoose_eta.push_back(fTwoProngCand_eta[index]);
      fTwoProngLoose_phi.push_back(fTwoProngCand_phi[index]);
      fTwoProngLoose_mass.push_back(fTwoProngCand_mass[index]);
      fTwoProngLoose_mass_l.push_back(fTwoProngCand_mass_l[index]);
      fTwoProngLoose_Mass0.push_back(fTwoProngCand_Mass0[index]);
      fTwoProngLoose_MassPi0.push_back(fTwoProngCand_MassPi0[index]);
      fTwoProngLoose_MassEta.push_back(fTwoProngCand_MassEta[index]);
      fTwoProngLoose_Mass300.push_back(fTwoProngCand_Mass300[index]);
      fTwoProngLoose_nExtraTracks.push_back(fTwoProngCand_nExtraTracks[index]);
      fTwoProngLoose_CHpos_pt.push_back(fTwoProngCand_CHpos_pt[index]);
      fTwoProngLoose_CHpos_eta.push_back(fTwoProngCand_CHpos_eta[index]);
      fTwoProngLoose_CHpos_phi.push_back(fTwoProngCand_CHpos_phi[index]);
      fTwoProngLoose_CHpos_mass.push_back(fTwoProngCand_CHpos_mass[index]);
      fTwoProngLoose_CHpos_dz.push_back(fTwoProngCand_CHpos_dz[index]);
      fTwoProngLoose_CHpos_dxy.push_back(fTwoProngCand_CHpos_dxy[index]);
      fTwoProngLoose_CHneg_pt.push_back(fTwoProngCand_CHneg_pt[index]);
      fTwoProngLoose_CHneg_eta.push_back(fTwoProngCand_CHneg_eta[index]);
      fTwoProngLoose_CHneg_phi.push_back(fTwoProngCand_CHneg_phi[index]);
      fTwoProngLoose_CHneg_mass.push_back(fTwoProngCand_CHneg_mass[index]);
      fTwoProngLoose_CHneg_dz.push_back(fTwoProngCand_CHneg_dz[index]);
      fTwoProngLoose_CHneg_dxy.push_back(fTwoProngCand_CHneg_dxy[index]);
      fTwoProngLoose_photon_pt.push_back(fTwoProngCand_photon_pt[index]);
      fTwoProngLoose_photon_eta.push_back(fTwoProngCand_photon_eta[index]);
      fTwoProngLoose_photon_phi.push_back(fTwoProngCand_photon_phi[index]);
      fTwoProngLoose_photon_mass.push_back(fTwoProngCand_photon_mass[index]);
      fTwoProngLoose_photon_pt_l.push_back(fTwoProngCand_photon_pt_l[index]);
      fTwoProngLoose_photon_eta_l.push_back(fTwoProngCand_photon_eta_l[index]);
      fTwoProngLoose_photon_phi_l.push_back(fTwoProngCand_photon_phi_l[index]);
      fTwoProngLoose_photon_mass_l.push_back(fTwoProngCand_photon_mass_l[index]);
      fTwoProngLoose_photon_nGamma.push_back(fTwoProngCand_photon_nGamma[index]);
      fTwoProngLoose_photon_nElectron.push_back(fTwoProngCand_photon_nElectron[index]);
      fTwoProngLoose_chargedIso.push_back(fTwoProngCand_chargedIso[index]);
      fTwoProngLoose_neutralIso.push_back(fTwoProngCand_neutralIso[index]);
      fTwoProngLoose_egammaIso.push_back(fTwoProngCand_egammaIso[index]);
      fTwoProngLoose_nChargedIsoCone.push_back(fTwoProngCand_nChargedIsoCone[index]);
      fTwoProngLoose_nNeutralIsoCone.push_back(fTwoProngCand_nNeutralIsoCone[index]);
      fTwoProngLoose_nEGammaIsoCone.push_back(fTwoProngCand_nChargedIsoCone[index]);
      fTwoProngLoose_genOmega_dR.push_back(fTwoProngCand_genOmega_dR[index]);
      fTwoProngLoose_genTau_dR.push_back(fTwoProngCand_genTau_dR[index]);
      fTwoProngLoose_trackAsym.push_back(fTwoProngCand_trackAsym[index]);
      fTwoProngLoose_photonAsym.push_back(fTwoProngCand_photonAsym[index]);
      fTwoProngLoose_mPosPho_l.push_back(fTwoProngCand_mPosPho_l[index]);
      fTwoProngLoose_mPosPho_pi0.push_back(fTwoProngCand_mPosPho_pi0[index]);
      fTwoProngLoose_mNegPho.push_back(fTwoProngCand_mNegPho[index]);
      fTwoProngLoose_mNegPho_l.push_back(fTwoProngCand_mNegPho_l[index]);
      fTwoProngLoose_mNegPho_pi0.push_back(fTwoProngCand_mNegPho_pi0[index]);
      fTwoProngLoose_mPosNeg.push_back(fTwoProngCand_mPosNeg[index]);
    }
    if (fTwoProngCand_tight[index])
    {
      // Candidate is tight and is next leading, fill all tight candidate collections
      fTwoProng_pt.push_back(fTwoProngCand_pt[index]);
      fTwoProng_eta.push_back(fTwoProngCand_eta[index]);
      fTwoProng_phi.push_back(fTwoProngCand_phi[index]);
      fTwoProng_mass.push_back(fTwoProngCand_mass[index]);
      fTwoProng_mass_l.push_back(fTwoProngCand_mass_l[index]);
      fTwoProng_Mass0.push_back(fTwoProngCand_Mass0[index]);
      fTwoProng_MassPi0.push_back(fTwoProngCand_MassPi0[index]);
      fTwoProng_MassEta.push_back(fTwoProngCand_MassEta[index]);
      fTwoProng_Mass300.push_back(fTwoProngCand_Mass300[index]);
      fTwoProng_nExtraTracks.push_back(fTwoProngCand_nExtraTracks[index]);
      fTwoProng_CHpos_pt.push_back(fTwoProngCand_CHpos_pt[index]);
      fTwoProng_CHpos_eta.push_back(fTwoProngCand_CHpos_eta[index]);
      fTwoProng_CHpos_phi.push_back(fTwoProngCand_CHpos_phi[index]);
      fTwoProng_CHpos_mass.push_back(fTwoProngCand_CHpos_mass[index]);
      fTwoProng_CHpos_dz.push_back(fTwoProngCand_CHpos_dz[index]);
      fTwoProng_CHpos_dxy.push_back(fTwoProngCand_CHpos_dxy[index]);
      fTwoProng_CHneg_pt.push_back(fTwoProngCand_CHneg_pt[index]);
      fTwoProng_CHneg_eta.push_back(fTwoProngCand_CHneg_eta[index]);
      fTwoProng_CHneg_phi.push_back(fTwoProngCand_CHneg_phi[index]);
      fTwoProng_CHneg_mass.push_back(fTwoProngCand_CHneg_mass[index]);
      fTwoProng_CHneg_dz.push_back(fTwoProngCand_CHneg_dz[index]);
      fTwoProng_CHneg_dxy.push_back(fTwoProngCand_CHneg_dxy[index]);
      fTwoProng_photon_pt.push_back(fTwoProngCand_photon_pt[index]);
      fTwoProng_photon_eta.push_back(fTwoProngCand_photon_eta[index]);
      fTwoProng_photon_phi.push_back(fTwoProngCand_photon_phi[index]);
      fTwoProng_photon_mass.push_back(fTwoProngCand_photon_mass[index]);
      fTwoProng_photon_pt_l.push_back(fTwoProngCand_photon_pt_l[index]);
      fTwoProng_photon_eta_l.push_back(fTwoProngCand_photon_eta_l[index]);
      fTwoProng_photon_phi_l.push_back(fTwoProngCand_photon_phi_l[index]);
      fTwoProng_photon_mass_l.push_back(fTwoProngCand_photon_mass_l[index]);
      fTwoProng_photon_nGamma.push_back(fTwoProngCand_photon_nGamma[index]);
      fTwoProng_photon_nElectron.push_back(fTwoProngCand_photon_nElectron[index]);
      fTwoProng_chargedIso.push_back(fTwoProngCand_chargedIso[index]);
      fTwoProng_neutralIso.push_back(fTwoProngCand_neutralIso[index]);
      fTwoProng_egammaIso.push_back(fTwoProngCand_egammaIso[index]);
      fTwoProng_nChargedIsoCone.push_back(fTwoProngCand_nChargedIsoCone[index]);
      fTwoProng_nNeutralIsoCone.push_back(fTwoProngCand_nNeutralIsoCone[index]);
      fTwoProng_nEGammaIsoCone.push_back(fTwoProngCand_nChargedIsoCone[index]);
      fTwoProng_genOmega_dR.push_back(fTwoProngCand_genOmega_dR[index]);
      fTwoProng_genTau_dR.push_back(fTwoProngCand_genTau_dR[index]);
      fTwoProng_trackAsym.push_back(fTwoProngCand_trackAsym[index]);
      fTwoProng_photonAsym.push_back(fTwoProngCand_photonAsym[index]);
      fTwoProng_mPosPho.push_back(fTwoProngCand_mPosPho[index]);
      fTwoProng_mPosPho_l.push_back(fTwoProngCand_mPosPho_l[index]);
      fTwoProng_mPosPho_pi0.push_back(fTwoProngCand_mPosPho_pi0[index]);
      fTwoProng_mNegPho.push_back(fTwoProngCand_mNegPho[index]);
      fTwoProng_mNegPho_l.push_back(fTwoProngCand_mNegPho_l[index]);
      fTwoProng_mNegPho_pi0.push_back(fTwoProngCand_mNegPho_pi0[index]);
      fTwoProng_mPosNeg.push_back(fTwoProngCand_mPosNeg[index]);
      fTwoProng_CHpos_p3.push_back(fTwoProngCand_CHpos_p3[index]);
      fTwoProng_CHneg_p3.push_back(fTwoProngCand_CHneg_p3[index]);
      fTwoProng_photon_p3.push_back(fTwoProngCand_photon_p3[index]);
    }
  }
  fnTwoProngCands = fTwoProngCand_pt.size();
  fnTwoProngs = fTwoProng_pt.size();
  fnTwoProngsLoose = fTwoProngLoose_pt.size();
  fnTwoProngsAsym = fTwoProngAsym_pt.size();
  fnTwoProngsAsymLoose = fTwoProngAsymLoose_pt.size();
  if (fDebug) cout << ". finished twoprong collections" << endl;

  // Z decay type
  if (fincludeZDecayGenParticles) {
    if (fDebug) cout << ". doing z decay gen particles" << endl;
    fzDecayType = TauHadFilters::ZDecayType(genparticles);

    vector<const reco::Candidate *> leptons;
    for (unsigned int i = 0; i < genparticles->size(); i++) {
      const reco::GenParticle & genparticle = (*genparticles)[i];
      if (genparticle.status() != 21 && genparticle.status() != 22) continue;
      leptons = TauHadFilters::getLeptonObjects(genparticle);
      break;
    }
    if (leptons.size() == 2)
    {
      if (fDebug) cout << ". . found both Z daughters" << endl;
      TLorentzVector lepton1;
      lepton1.SetPtEtaPhiM(leptons[0]->pt(), leptons[0]->eta(), leptons[0]->phi(), leptons[0]->mass());
      TLorentzVector lepton2;
      lepton2.SetPtEtaPhiM(leptons[1]->pt(), leptons[1]->eta(), leptons[1]->phi(), leptons[1]->mass());
      fGenZ_pt = (lepton1+lepton2).Pt();
      fGenZ_phi = (lepton1+lepton2).Phi();
      fGenZ_mass = (lepton1+lepton2).M();
      if (fabs((lepton1+lepton2).CosTheta()) == 1) fGenZ_eta = 1000.0;
      else fGenZ_eta = (lepton1+lepton2).Eta();
      if (fDebug) cout << ". . done with Z daughters" << endl;
    }

    // matching, by gen tau perspective now
    for (unsigned int i = 0; i < genparticles->size(); i++) {
      if (fDebug) cout << ". . now matching hadronic taus" << endl;
      const reco::GenParticle & genparticle = (*genparticles)[i];
      if (!TauHadFilters::isHadronicTau(&genparticle)) continue; // only including hadronically decaying taus in gen tau collection
      TLorentzVector GenParticle;
      GenParticle.SetPtEtaPhiM(genparticle.pt(), genparticle.eta(), genparticle.phi(), genparticle.mass());
      fGenTau_pt.push_back(genparticle.pt());
      fGenTau_eta.push_back(genparticle.eta());
      fGenTau_phi.push_back(genparticle.phi());
      fGenTau_mass.push_back(genparticle.mass());
      double candDR = 99.9;
      double passedCandDR = 99.9;
      for (unsigned int j = 0; j < fTwoProngCand_pt.size(); j++) {
        TLorentzVector Candidate;
        Candidate.SetPtEtaPhiM(fTwoProngCand_pt[j], fTwoProngCand_eta[j], fTwoProngCand_phi[j], fTwoProngCand_mass[j]);
        double dr = Candidate.DeltaR(GenParticle);
        if (dr < candDR) candDR = dr;
      }
      for (unsigned int j = 0; j < fTwoProng_pt.size(); j++) {
        TLorentzVector PassedCandidate;
        PassedCandidate.SetPtEtaPhiM(fTwoProng_pt[j], fTwoProng_eta[j], fTwoProng_phi[j], fTwoProng_mass[j]);
        double dr = PassedCandidate.DeltaR(GenParticle);
        if (dr < passedCandDR) passedCandDR = dr;
      }
      fGenTau_objDR.push_back(passedCandDR);
      fGenTau_candobjDR.push_back(candDR);
    }
    fnGenTaus = fGenTau_pt.size();
    if (fDebug) cout << ". done z decay gen particles" << endl;
  }

  // matching, by gen omega perspective now
  if (fincludeSignalGenParticles) {
    if (fDebug) cout << ". doing phi signal gen particles" << endl;
    for (unsigned int i = 0; i < genparticles->size(); i++) {
      const reco::GenParticle &genparticle = (*genparticles)[i];
      if ((genparticle.pdgId() != 54 && genparticle.pdgId() != 9000006) || genparticle.status() != 62) continue;
      for (unsigned int j = 0; j < genparticle.numberOfDaughters(); j++) {
        const reco::Candidate* genparticle2 = genparticle.daughter(j);
        TLorentzVector GenParticle;
        GenParticle.SetPtEtaPhiM(genparticle2->pt(), genparticle2->eta(), genparticle2->phi(), genparticle2->mass());
        if (fOldData && !(genparticle2->pdgId() == 221 || genparticle2->pdgId() == 331)) continue;
        if (!fOldData && !(genparticle2->pdgId() == 90000054 || genparticle2->pdgId() == 90000055)) continue;
        double candDR = 99.9;
        double passedCandDR = 99.9;
        double jetDR = 99.9;
        double pattauDR = 99.9;
        double pattau1DR = 99.9;
        double pattau2DR = 99.9;
        double pattau3DR = 99.9;
        for (unsigned int j = 0; j < fTwoProngCand_pt.size(); j++) {
          TLorentzVector Candidate;
          Candidate.SetPtEtaPhiM(fTwoProngCand_pt[j], fTwoProngCand_eta[j], fTwoProngCand_phi[j], fTwoProngCand_mass[j]);
          double dr = Candidate.DeltaR(GenParticle);
          if (dr < candDR) candDR = dr;
        }
        for (unsigned int j = 0; j < fTwoProng_pt.size(); j++) {
          if (!fTwoProngCand_tight[j]) continue;
          TLorentzVector PassedCandidate;
          PassedCandidate.SetPtEtaPhiM(fTwoProng_pt[j], fTwoProng_eta[j], fTwoProng_phi[j], fTwoProng_mass[j]);
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
        for (unsigned int i = 0; i < taus->size(); i++) {
          const pat::Tau &tau = (*taus)[i];
          TLorentzVector Pattau;
          Pattau.SetPtEtaPhiM(tau.pt(), tau.eta(), tau.phi(), tau.mass());
          if (tau.pt() <= 20) continue;
          if (fabs(tau.eta()) >= 2.3) continue;
          if (tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT") >= 0.5) {
            if (Pattau.DeltaR(GenParticle) < pattauDR) pattauDR = Pattau.DeltaR(GenParticle);
            if ((tau.decayMode() == 0 || tau.decayMode() == 1 || tau.decayMode() == 2) && Pattau.DeltaR(GenParticle) < pattau1DR) pattau1DR = Pattau.DeltaR(GenParticle);
            if ((tau.decayMode() == 5 || tau.decayMode() == 6 || tau.decayMode() == 7) && Pattau.DeltaR(GenParticle) < pattau2DR) pattau2DR = Pattau.DeltaR(GenParticle);
            if ((tau.decayMode() == 10 || tau.decayMode() == 11) && Pattau.DeltaR(GenParticle) < pattau3DR) pattau3DR = Pattau.DeltaR(GenParticle);
          }
        }
        fGenOmega_objDR.push_back(passedCandDR);
        fGenOmega_candobjDR.push_back(candDR);
        fGenOmega_jetDR.push_back(jetDR);
        fGenOmega_pattauDR.push_back(pattauDR);
        fGenOmega_pattau1DR.push_back(pattau1DR);
        fGenOmega_pattau2DR.push_back(pattau2DR);
        fGenOmega_pattau3DR.push_back(pattau3DR);
      } // end loop on daughters of Phi gen particle
    } // end loop on all gen particles
    if (fDebug) cout << ". done phi signal gen particles" << endl;
  } // end conditional for signal mc branches

  // High-pT Photon id  
  if (fDebug) cout << ". high-pt-id" << endl;
  const CaloSubdetectorTopology* subDetTopologyEB_;
  const CaloSubdetectorTopology* subDetTopologyEE_;
  edm::ESHandle<CaloTopology> caloTopology;
  iSetup.get<CaloTopologyRecord>().get(caloTopology);
  subDetTopologyEB_ = caloTopology->getSubdetectorTopology(DetId::Ecal,EcalBarrel);
  subDetTopologyEE_ = caloTopology->getSubdetectorTopology(DetId::Ecal,EcalEndcap);
  bool isSat = false;
  std::vector<edm::Ptr<pat::Photon>> basePhotons;
  std::vector<edm::Ptr<pat::Photon>> loose1Photons;
  std::vector<edm::Ptr<pat::Photon>> loose2Photons;
  std::vector<edm::Ptr<pat::Photon>> loose3Photons;
  std::vector<edm::Ptr<pat::Photon>> loose4Photons;
  std::vector<edm::Ptr<pat::Photon>> loose5Photons;
  std::vector<edm::Ptr<pat::Photon>> goodPhotons;
  std::vector<edm::Ptr<pat::Photon>> coneHEphotons;
  for (size_t i = 0; i < ged_photons->size(); ++i) {
    const auto pho = ged_photons->ptrAt(i);
    isSat = photon_isSaturated(&(*pho), &(*recHitsEB), &(*recHitsEE), &(*subDetTopologyEB_), &(*subDetTopologyEE_));
    bool passID = photon_passHighPtID(&(*pho), fRho, isSat);
    if(passID) {
      goodPhotons.push_back(pho);
    }
    bool passIDConeHE = photon_passHighPtID_coneHE(&(*pho), fRho, isSat);
    if(passIDConeHE) {
      coneHEphotons.push_back(pho);
    }
    bool passLoose1ID = photon_passHighPtID_loose1(&(*pho), fRho, isSat);
    if(passLoose1ID) {
      loose1Photons.push_back(pho);
    }
    bool passLoose2ID = photon_passHighPtID_loose2(&(*pho), fRho, isSat);
    if(passLoose2ID) {
      loose2Photons.push_back(pho);
    }
    bool passLoose3ID = photon_passHighPtID_loose3(&(*pho), fRho, isSat);
    if(passLoose3ID) {
      loose3Photons.push_back(pho);
    }
    bool passLoose4ID = photon_passHighPtID_loose4(&(*pho), fRho, isSat);
    if(passLoose4ID) {
      loose4Photons.push_back(pho);
    }
    bool passLoose5ID = photon_passHighPtID_loose5(&(*pho), fRho, isSat);
    if(passLoose5ID) {
      loose5Photons.push_back(pho);
    }
    bool passBaseID = photon_passHighPtID_base(&(*pho), fRho, isSat);
    if(passBaseID) {
      basePhotons.push_back(pho);
    }
  }
  sort(loose1Photons.begin(),loose1Photons.end(),compareCandsByPt);
  sort(loose2Photons.begin(),loose2Photons.end(),compareCandsByPt);
  sort(loose3Photons.begin(),loose3Photons.end(),compareCandsByPt);
  sort(loose4Photons.begin(),loose4Photons.end(),compareCandsByPt);
  sort(loose5Photons.begin(),loose5Photons.end(),compareCandsByPt);
  sort(goodPhotons.begin(),goodPhotons.end(),compareCandsByPt);
  sort(coneHEphotons.begin(),coneHEphotons.end(),compareCandsByPt);
  for (unsigned int i = 0; i < basePhotons.size(); i++ )
  {
    fBaseIDPhoton_pt.push_back( (*basePhotons[i]).pt() );
    fBaseIDPhoton_eta.push_back( (*basePhotons[i]).eta() );
    fBaseIDPhoton_scEta.push_back( photon_scEta(&(*basePhotons[i])) );
    fBaseIDPhoton_phi.push_back( (*basePhotons[i]).phi() );
    fBaseIDPhoton_mass.push_back( (*basePhotons[i]).mass() );
    fBaseIDPhoton_isoGamma.push_back( photon_computeIsoGamma(&(*basePhotons[i]), *rhoH) );
    fBaseIDPhoton_isoCh.push_back( photon_computeIsoCh(&(*basePhotons[i])) );
    fBaseIDPhoton_HE.push_back( photon_computeHE(&(*basePhotons[i])) );
    fBaseIDPhoton_coneHE.push_back( photon_computeHE_coneBased(&(*basePhotons[i])) );
    fBaseIDPhoton_sigmaIetaIeta.push_back( photon_computeSigmaIetaIeta(&(*basePhotons[i])) );
  }
  for (unsigned int i = 0; i < loose1Photons.size(); i++ )
  {
    if (photon_scEta(&(*loose1Photons[i])) < 1.4442) {
      fLoose1IDPhoton_pt.push_back( (*loose1Photons[i]).pt() );
      fLoose1IDPhoton_eta.push_back( (*loose1Photons[i]).eta() );
      fLoose1IDPhoton_scEta.push_back( photon_scEta(&(*loose1Photons[i])) );
      fLoose1IDPhoton_phi.push_back( (*loose1Photons[i]).phi() );
      fLoose1IDPhoton_mass.push_back( (*loose1Photons[i]).mass() );
      fLoose1IDPhoton_isoGamma.push_back( photon_computeIsoGamma(&(*loose1Photons[i]), *rhoH) );
      fLoose1IDPhoton_isoCh.push_back( photon_computeIsoCh(&(*loose1Photons[i])) );
      fLoose1IDPhoton_HE.push_back( photon_computeHE(&(*loose1Photons[i])) );
      fLoose1IDPhoton_coneHE.push_back( photon_computeHE_coneBased(&(*loose1Photons[i])) );
      fLoose1IDPhoton_sigmaIetaIeta.push_back( photon_computeSigmaIetaIeta(&(*loose1Photons[i])) );
    } else if (photon_scEta(&(*loose1Photons[i])) > 1.566 && photon_scEta(&(*loose1Photons[i])) < 2.5) {
      fLoose1IDPhotonEndcap_pt.push_back( (*loose1Photons[i]).pt() );
      fLoose1IDPhotonEndcap_eta.push_back( (*loose1Photons[i]).eta() );
      fLoose1IDPhotonEndcap_scEta.push_back( photon_scEta(&(*loose1Photons[i])) );
      fLoose1IDPhotonEndcap_phi.push_back( (*loose1Photons[i]).phi() );
      fLoose1IDPhotonEndcap_mass.push_back( (*loose1Photons[i]).mass() );
      fLoose1IDPhotonEndcap_isoGamma.push_back( photon_computeIsoGamma(&(*loose1Photons[i]), *rhoH) );
      fLoose1IDPhotonEndcap_isoCh.push_back( photon_computeIsoCh(&(*loose1Photons[i])) );
      fLoose1IDPhotonEndcap_HE.push_back( photon_computeHE(&(*loose1Photons[i])) );
      fLoose1IDPhotonEndcap_coneHE.push_back( photon_computeHE_coneBased(&(*loose1Photons[i])) );
      fLoose1IDPhotonEndcap_sigmaIetaIeta.push_back( photon_computeSigmaIetaIeta(&(*loose1Photons[i])) );
    }
  }
  for (unsigned int i = 0; i < loose2Photons.size(); i++ )
  {
    if (photon_scEta(&(*loose2Photons[i])) < 1.4442) {
      fLoose2IDPhoton_pt.push_back( (*loose2Photons[i]).pt() );
      fLoose2IDPhoton_eta.push_back( (*loose2Photons[i]).eta() );
      fLoose2IDPhoton_scEta.push_back( photon_scEta(&(*loose2Photons[i])) );
      fLoose2IDPhoton_phi.push_back( (*loose2Photons[i]).phi() );
      fLoose2IDPhoton_mass.push_back( (*loose2Photons[i]).mass() );
      fLoose2IDPhoton_isoGamma.push_back( photon_computeIsoGamma(&(*loose2Photons[i]), *rhoH) );
      fLoose2IDPhoton_isoCh.push_back( photon_computeIsoCh(&(*loose2Photons[i])) );
      fLoose2IDPhoton_HE.push_back( photon_computeHE(&(*loose2Photons[i])) );
      fLoose2IDPhoton_coneHE.push_back( photon_computeHE_coneBased(&(*loose2Photons[i])) );
      fLoose2IDPhoton_sigmaIetaIeta.push_back( photon_computeSigmaIetaIeta(&(*loose2Photons[i])) );
    } else if (photon_scEta(&(*loose2Photons[i])) > 1.566 && photon_scEta(&(*loose2Photons[i])) < 2.5) {
      fLoose2IDPhotonEndcap_pt.push_back( (*loose2Photons[i]).pt() );
      fLoose2IDPhotonEndcap_eta.push_back( (*loose2Photons[i]).eta() );
      fLoose2IDPhotonEndcap_scEta.push_back( photon_scEta(&(*loose2Photons[i])) );
      fLoose2IDPhotonEndcap_phi.push_back( (*loose2Photons[i]).phi() );
      fLoose2IDPhotonEndcap_mass.push_back( (*loose2Photons[i]).mass() );
      fLoose2IDPhotonEndcap_isoGamma.push_back( photon_computeIsoGamma(&(*loose2Photons[i]), *rhoH) );
      fLoose2IDPhotonEndcap_isoCh.push_back( photon_computeIsoCh(&(*loose2Photons[i])) );
      fLoose2IDPhotonEndcap_HE.push_back( photon_computeHE(&(*loose2Photons[i])) );
      fLoose2IDPhotonEndcap_coneHE.push_back( photon_computeHE_coneBased(&(*loose2Photons[i])) );
      fLoose2IDPhotonEndcap_sigmaIetaIeta.push_back( photon_computeSigmaIetaIeta(&(*loose2Photons[i])) );
    }
  }
  for (unsigned int i = 0; i < loose3Photons.size(); i++ )
  {
    if (photon_scEta(&(*loose3Photons[i])) < 1.4442) {
      fLoose3IDPhoton_pt.push_back( (*loose3Photons[i]).pt() );
      fLoose3IDPhoton_eta.push_back( (*loose3Photons[i]).eta() );
      fLoose3IDPhoton_scEta.push_back( photon_scEta(&(*loose3Photons[i])) );
      fLoose3IDPhoton_phi.push_back( (*loose3Photons[i]).phi() );
      fLoose3IDPhoton_mass.push_back( (*loose3Photons[i]).mass() );
      fLoose3IDPhoton_isoGamma.push_back( photon_computeIsoGamma(&(*loose3Photons[i]), *rhoH) );
      fLoose3IDPhoton_isoCh.push_back( photon_computeIsoCh(&(*loose3Photons[i])) );
      fLoose3IDPhoton_HE.push_back( photon_computeHE(&(*loose3Photons[i])) );
      fLoose3IDPhoton_coneHE.push_back( photon_computeHE_coneBased(&(*loose3Photons[i])) );
      fLoose3IDPhoton_sigmaIetaIeta.push_back( photon_computeSigmaIetaIeta(&(*loose3Photons[i])) );
    }
  }
  for (unsigned int i = 0; i < loose4Photons.size(); i++ )
  {
    if (photon_scEta(&(*loose4Photons[i])) < 1.4442) {
      fLoose4IDPhoton_pt.push_back( (*loose4Photons[i]).pt() );
      fLoose4IDPhoton_eta.push_back( (*loose4Photons[i]).eta() );
      fLoose4IDPhoton_scEta.push_back( photon_scEta(&(*loose4Photons[i])) );
      fLoose4IDPhoton_phi.push_back( (*loose4Photons[i]).phi() );
      fLoose4IDPhoton_mass.push_back( (*loose4Photons[i]).mass() );
      fLoose4IDPhoton_isoGamma.push_back( photon_computeIsoGamma(&(*loose4Photons[i]), *rhoH) );
      fLoose4IDPhoton_isoCh.push_back( photon_computeIsoCh(&(*loose4Photons[i])) );
      fLoose4IDPhoton_HE.push_back( photon_computeHE(&(*loose4Photons[i])) );
      fLoose4IDPhoton_coneHE.push_back( photon_computeHE_coneBased(&(*loose4Photons[i])) );
      fLoose4IDPhoton_sigmaIetaIeta.push_back( photon_computeSigmaIetaIeta(&(*loose4Photons[i])) );
    }
  }
  for (unsigned int i = 0; i < loose5Photons.size(); i++ )
  {
    if (photon_scEta(&(*loose5Photons[i])) < 1.4442) {
      fLoose5IDPhoton_pt.push_back( (*loose5Photons[i]).pt() );
      fLoose5IDPhoton_eta.push_back( (*loose5Photons[i]).eta() );
      fLoose5IDPhoton_scEta.push_back( photon_scEta(&(*loose5Photons[i])) );
      fLoose5IDPhoton_phi.push_back( (*loose5Photons[i]).phi() );
      fLoose5IDPhoton_mass.push_back( (*loose5Photons[i]).mass() );
      fLoose5IDPhoton_isoGamma.push_back( photon_computeIsoGamma(&(*loose5Photons[i]), *rhoH) );
      fLoose5IDPhoton_isoCh.push_back( photon_computeIsoCh(&(*loose5Photons[i])) );
      fLoose5IDPhoton_HE.push_back( photon_computeHE(&(*loose5Photons[i])) );
      fLoose5IDPhoton_coneHE.push_back( photon_computeHE_coneBased(&(*loose5Photons[i])) );
      fLoose5IDPhoton_sigmaIetaIeta.push_back( photon_computeSigmaIetaIeta(&(*loose5Photons[i])) );
      fLoose5IDPhoton_passVeto.push_back( loose5Photons[i]->passElectronVeto() );
    }
  }
  for (unsigned int i = 0; i < goodPhotons.size(); i++ )
  {
    if (photon_scEta(&(*goodPhotons[i])) < 1.4442) {
      fIDPhoton_pt.push_back( (*goodPhotons[i]).pt() );
      fIDPhoton_eta.push_back( (*goodPhotons[i]).eta() );
      fIDPhoton_scEta.push_back( photon_scEta(&(*goodPhotons[i])) );
      fIDPhoton_phi.push_back( (*goodPhotons[i]).phi() );
      fIDPhoton_mass.push_back( (*goodPhotons[i]).mass() );
      fIDPhoton_isoGamma.push_back( photon_computeIsoGamma(&(*goodPhotons[i]), *rhoH) );
      fIDPhoton_isoCh.push_back( photon_computeIsoCh(&(*goodPhotons[i])) );
      fIDPhoton_HE.push_back( photon_computeHE(&(*goodPhotons[i])) );
      fIDPhoton_coneHE.push_back( photon_computeHE_coneBased(&(*goodPhotons[i])) );
      fIDPhoton_sigmaIetaIeta.push_back( photon_computeSigmaIetaIeta(&(*goodPhotons[i])) );
    } else if (photon_scEta(&(*goodPhotons[i])) > 1.566 && photon_scEta(&(*goodPhotons[i])) < 2.5) {
      fIDPhotonEndcap_pt.push_back( (*goodPhotons[i]).pt() );
      fIDPhotonEndcap_eta.push_back( (*goodPhotons[i]).eta() );
      fIDPhotonEndcap_scEta.push_back( photon_scEta(&(*goodPhotons[i])) );
      fIDPhotonEndcap_phi.push_back( (*goodPhotons[i]).phi() );
      fIDPhotonEndcap_mass.push_back( (*goodPhotons[i]).mass() );
      fIDPhotonEndcap_isoGamma.push_back( photon_computeIsoGamma(&(*goodPhotons[i]), *rhoH) );
      fIDPhotonEndcap_isoCh.push_back( photon_computeIsoCh(&(*goodPhotons[i])) );
      fIDPhotonEndcap_HE.push_back( photon_computeHE(&(*goodPhotons[i])) );
      fIDPhotonEndcap_coneHE.push_back( photon_computeHE_coneBased(&(*goodPhotons[i])) );
      fIDPhotonEndcap_sigmaIetaIeta.push_back( photon_computeSigmaIetaIeta(&(*goodPhotons[i])) );
    }
    
  }
  for (unsigned int i = 0; i < coneHEphotons.size(); i++ )
  {
    if (photon_scEta(&(*coneHEphotons[i])) < 1.4442) {
      fConeHEIDPhoton_pt.push_back( (*coneHEphotons[i]).pt() );
      fConeHEIDPhoton_eta.push_back( (*coneHEphotons[i]).eta() );
      fConeHEIDPhoton_scEta.push_back( photon_scEta(&(*coneHEphotons[i])) );
      fConeHEIDPhoton_phi.push_back( (*coneHEphotons[i]).phi() );
      fConeHEIDPhoton_mass.push_back( (*coneHEphotons[i]).mass() );
      fConeHEIDPhoton_isoGamma.push_back( photon_computeIsoGamma(&(*coneHEphotons[i]), *rhoH) );
      fConeHEIDPhoton_isoCh.push_back( photon_computeIsoCh(&(*coneHEphotons[i])) );
      fConeHEIDPhoton_HE.push_back( photon_computeHE(&(*coneHEphotons[i])) );
      fConeHEIDPhoton_coneHE.push_back( photon_computeHE_coneBased(&(*coneHEphotons[i])) );
      fConeHEIDPhoton_sigmaIetaIeta.push_back( photon_computeSigmaIetaIeta(&(*coneHEphotons[i])) );
    } else if (photon_scEta(&(*coneHEphotons[i])) > 1.566 && photon_scEta(&(*coneHEphotons[i])) < 2.5) {
      fConeHEIDPhotonEndcap_pt.push_back( (*coneHEphotons[i]).pt() );
      fConeHEIDPhotonEndcap_eta.push_back( (*coneHEphotons[i]).eta() );
      fConeHEIDPhotonEndcap_scEta.push_back( photon_scEta(&(*coneHEphotons[i])) );
      fConeHEIDPhotonEndcap_phi.push_back( (*coneHEphotons[i]).phi() );
      fConeHEIDPhotonEndcap_mass.push_back( (*coneHEphotons[i]).mass() );
      fConeHEIDPhotonEndcap_isoGamma.push_back( photon_computeIsoGamma(&(*coneHEphotons[i]), *rhoH) );
      fConeHEIDPhotonEndcap_isoCh.push_back( photon_computeIsoCh(&(*coneHEphotons[i])) );
      fConeHEIDPhotonEndcap_HE.push_back( photon_computeHE(&(*coneHEphotons[i])) );
      fConeHEIDPhotonEndcap_coneHE.push_back( photon_computeHE_coneBased(&(*coneHEphotons[i])) );
      fConeHEIDPhotonEndcap_sigmaIetaIeta.push_back( photon_computeSigmaIetaIeta(&(*coneHEphotons[i])) );
    }
  }
  fNumBaseIDPhotons = fBaseIDPhoton_pt.size();
  fNumIDPhotons = fIDPhoton_pt.size();
  fNumIDPhotonsEndcap = fIDPhotonEndcap_pt.size();
  fNumConeHEIDPhotons = fConeHEIDPhoton_pt.size();
  fNumConeHEIDPhotonsEndcap = fConeHEIDPhotonEndcap_pt.size();
  fNumLoose1IDPhotons = fLoose1IDPhoton_pt.size();
  fNumLoose1IDPhotonsEndcap = fLoose1IDPhotonEndcap_pt.size();
  fNumLoose2IDPhotons = fLoose2IDPhoton_pt.size();
  fNumLoose2IDPhotonsEndcap = fLoose2IDPhotonEndcap_pt.size();
  fNumLoose3IDPhotons = fLoose3IDPhoton_pt.size();
  fNumLoose4IDPhotons = fLoose4IDPhoton_pt.size();
  fNumLoose5IDPhotons = fLoose5IDPhoton_pt.size();
  // old style structs
  InitRecoPhotonInfo(fRecoTightPhotonInfo1);
  InitRecoPhotonInfo(fRecoTightPhotonInfo2);
  InitRecoPhotonInfo(fRecoTightPhotonInfo3);
  if (goodPhotons.size() > 0)
    TwoProngAnalysis::FillRecoPhotonInfo(fRecoTightPhotonInfo1,&(*goodPhotons[0]),lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent,iSetup); 
  if (goodPhotons.size() > 1)
    TwoProngAnalysis::FillRecoPhotonInfo(fRecoTightPhotonInfo2,&(*goodPhotons[1]),lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent,iSetup); 
  if (goodPhotons.size() > 2)
    TwoProngAnalysis::FillRecoPhotonInfo(fRecoTightPhotonInfo3,&(*goodPhotons[2]),lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent,iSetup); 
  if (fDebug) cout << ". done high-pt-id photons" << endl;

  // Sideband Category
  fSCat = 0;
  int T = fnTwoProngs;
  int L = fnTwoProngsLoose;
  int A = fnTwoProngsAsym;
  int B = fnTwoProngsAsymLoose;
  if (T==0 && L>=1 && A==0 && B==0) fSCat = 1; // Le
  if (T==0 && L==0 && A>=1 && B==0) fSCat = 2; // Ae
  if (T==0 && L==0 && A==0 && B>=1) fSCat = 3; // Be
  if (T==0 && L>=1 && A>=1 && B==0) fSCat = 4; // LAe
  if (T==0 && L==0 && A>=1 && B>=1) fSCat = 5; // ABe
  if (T==0 && L>=1 && A==0 && B>=1) fSCat = 6; // LBe
  if (T==0 && L>=1 && A>=1 && B>=1) fSCat = 7; // All three
  // as above plus has tight
  if (T>=1 && L==0 && A==0 && B==0) fSCat = 10; // T
  if (T>=1 && L>=1 && A==0 && B==0) fSCat = 11; // Le
  if (T>=1 && L==0 && A>=1 && B==0) fSCat = 12; // Ae
  if (T>=1 && L==0 && A==0 && B>=1) fSCat = 13; // Be
  if (T>=1 && L>=1 && A>=1 && B==0) fSCat = 14; // LAe
  if (T>=1 && L==0 && A>=1 && B>=1) fSCat = 15; // ABe
  if (T>=1 && L>=1 && A==0 && B>=1) fSCat = 16; // LBe
  if (T>=1 && L>=1 && A>=1 && B>=1) fSCat = 17; // All three

  // Construct Di-Objects
  // Di-TwoProng
  InitRecoDiObjectInfo(fRecoPhiDiTwoProng);
  if (fnTwoProngs >= 2)
  {
    TLorentzVector Eta1;
    Eta1.SetPtEtaPhiM(fTwoProng_pt[0], fTwoProng_eta[0], fTwoProng_phi[0], fTwoProng_mass[0]);
    TLorentzVector Eta2;
    Eta2.SetPtEtaPhiM(fTwoProng_pt[1], fTwoProng_eta[1], fTwoProng_phi[1], fTwoProng_mass[1]);
    FillRecoDiObjectInfo(fRecoPhiDiTwoProng, Eta1, Eta2);
  }

  // photon plus TwoProng
  InitRecoDiObjectInfo(fRecoPhiPhotonTwoProng);
  if (fnTwoProngs >= 1 && fNumIDPhotons >= 1)
  {
    TLorentzVector LeadingTwoProng;
    LeadingTwoProng.SetPtEtaPhiM(fTwoProng_pt[0], fTwoProng_eta[0], fTwoProng_phi[0], fTwoProng_mass[0]);
    TLorentzVector LeadingPhoton;
    LeadingPhoton.SetPtEtaPhiM(fIDPhoton_pt[0], fIDPhoton_eta[0], fIDPhoton_phi[0], fIDPhoton_mass[0]);
    FillRecoDiObjectInfo(fRecoPhiPhotonTwoProng, LeadingTwoProng, LeadingPhoton);
  }

  // TwoProng plus (photon or TwoProng) inclusive, use higher pt
  InitRecoDiObjectInfo(fRecoPhiInclusive);
  if (fnTwoProngs >=1 && (fnTwoProngs + fNumIDPhotons >= 2))
  {
    TLorentzVector LeadingTwoProng;
    LeadingTwoProng.SetPtEtaPhiM(fTwoProng_pt[0], fTwoProng_eta[0], fTwoProng_phi[0], fTwoProng_mass[0]);
    TLorentzVector SubLeadingTwoProng;
    TLorentzVector LeadingPhoton;
    TLorentzVector LeadingSecondary;
    if (fNumIDPhotons >= 1) LeadingPhoton.SetPtEtaPhiM(fIDPhoton_pt[0], fIDPhoton_eta[0], fIDPhoton_phi[0], fIDPhoton_mass[0]);
    if (fnTwoProngs >= 2) SubLeadingTwoProng.SetPtEtaPhiM(fTwoProng_pt[1], fTwoProng_eta[1], fTwoProng_phi[1], fTwoProng_mass[1]);
  
    if (fNumIDPhotons == 0) LeadingSecondary = SubLeadingTwoProng;
    if (fnTwoProngs == 1) LeadingSecondary = LeadingPhoton;
    if (fnTwoProngs >= 2 && fNumIDPhotons >=1) {
      if (SubLeadingTwoProng.Pt() > LeadingPhoton.Pt()) LeadingSecondary = SubLeadingTwoProng;
      else LeadingSecondary = LeadingPhoton;
    }
    FillRecoDiObjectInfo(fRecoPhiInclusive, LeadingTwoProng, LeadingSecondary);
  }

  // Z preselection tag and probe branches
  if (fDebug) cout << ". doing z preselection branches" << endl;
  TauHadFilters::TauHadPreSelectionResult result;
  result = TauHadFilters::computePreSelectionResult(iEvent, triggerBits, triggerObjects, triggerPrescales, vertices, taus, muons, electrons, ak4jets, MET, rhoH, fusePatTauForZPreBranches, static_cast<TauHadFilters::muonIDtype>(fmuonIDtype), static_cast<TauHadFilters::muonISOtype>(fmuonISOtype));

  fpassMuonTrigger = result.passTrigger;
  fpassMuonTriggerTk = result.passTriggerTk;
  fMuonTrigger = result.foundTrigger;
  fMuonTriggerTk = result.foundTriggerTk;
  fnTagMuons = result.nTagMuons;
  fnProbeTaus = result.nProbeTaus;
  fpassMuonTauPair = result.passMuonTauPair;
  fpassPzeta = result.pairAndPassPzeta;
  fpassMT = result.pairAndPassMT;
  fpassExtraElectronVeto = result.passExtraElectronVeto;
  fpassExtraMuonVeto = result.passExtraMuonVeto;
  fpassDiMuonVeto = result.passDiMuonVeto;
  fpassExtraLeptonVeto = result.passExtraElectronVeto && result.passExtraMuonVeto && result.passDiMuonVeto;
  fpassBtagVeto = result.passBtagVeto;
  fMT = result.MT;
  fPzeta = result.Pzeta;
  fhighestBtag = result.highestBtagDiscriminant;
  fpassPreselection = result.passPreSelection;
  
  fpassReducedSelection = result.passMuonTauPair && result.passExtraElectronVeto && result.passExtraMuonVeto && result.passDiMuonVeto && result.passBtagVeto;

  if (fpassMuonTauPair) {
    fTagMuon_pt = result.tagMuon->pt();
    fTagMuon_eta = result.tagMuon->eta();
    fTagMuon_phi = result.tagMuon->phi();
    fTagMuon_mass = result.tagMuon->mass();
    fTagMuon_z = (result.tagMuon->muonBestTrack())->vz();
    fTagMuon_dz = fabs( (result.tagMuon->muonBestTrack())->dz( PV.position() ) );
    fTagMuon_dB = fabs( result.tagMuon->dB() );
    fTagMuon_dxy = fabs( (result.tagMuon->muonBestTrack())->dxy( PV.position() ) );
    fTagMuon_iso = TauHadFilters::computeMuonIsolation(result.tagMuon);
    fProbeTau_pt = result.usePatTau ? result.probeTau->pt() : result.probeTauJet->pt();
    fProbeTau_eta = result.usePatTau ? result.probeTau->eta() : result.probeTauJet->eta();
    fProbeTau_phi = result.usePatTau ? result.probeTau->phi() : result.probeTauJet->phi();
    fProbeTau_mass = result.usePatTau ? result.probeTau->mass() : result.probeTauJet->mass();
    if(fincludeZDecayGenParticles) {
      TLorentzVector ProbeTau_p4;
      ProbeTau_p4.SetPtEtaPhiM(fProbeTau_pt, fProbeTau_eta, fProbeTau_phi, fProbeTau_mass);
      double dr_to_gen = 99.9;
      for (unsigned int i = 0; i < fGenTau_pt.size(); i++) {
        TLorentzVector GenTau_p4;
        GenTau_p4.SetPtEtaPhiM(fGenTau_pt[i], fGenTau_eta[i], fGenTau_phi[i], fGenTau_mass[i]);
        double dr_to_gen_i = ProbeTau_p4.DeltaR(GenTau_p4);
        if (dr_to_gen_i < 0.4 && dr_to_gen_i < dr_to_gen) {
          dr_to_gen = dr_to_gen_i;
        } 
      }
      fProbeTau_genDR = dr_to_gen;
    }
  } else {
    fTagMuon_pt = -999.9;
    fTagMuon_eta = -999.9;
    fTagMuon_phi = -999.9;
    fTagMuon_mass = -999.9;
    fTagMuon_z = -999.9;
    fTagMuon_dz = -999.9;
    fTagMuon_dB = -999.9;
    fTagMuon_dxy = -999.9;
    fTagMuon_iso = -999.9;
    fProbeTau_pt = -999.9;
    fProbeTau_eta = -999.9;
    fProbeTau_phi = -999.9;
    fProbeTau_mass = -999.9;
    fProbeTau_genDR = -999.9;
  }

  InitRecoDiObjectInfo(fMuonTwoProng);
  InitRecoDiObjectInfo(fMuonTauID);
  InitRecoDiObjectInfo(fMuonProbe);
  if (result.passMuonTauPair) {
    // visible Z with probe tau
    TLorentzVector z_muon;
    z_muon.SetPtEtaPhiM(result.tagMuon->pt(), result.tagMuon->eta(), result.tagMuon->phi(), result.tagMuon->mass());
    TLorentzVector z_probetau;
    if (result.usePatTau) z_probetau.SetPtEtaPhiM(result.probeTau->pt(), result.probeTau->eta(), result.probeTau->phi(), result.probeTau->mass());
    else z_probetau.SetPtEtaPhiM(result.probeTauJet->pt(), result.probeTauJet->eta(), result.probeTauJet->phi(), result.probeTauJet->mass());
    FillRecoDiObjectInfo(fMuonProbe, z_muon, z_probetau);

    // visible Z with twoprong
    double closest_dr_2prong = 100.0;
    unsigned int index_closest_dr_2prong = -1;
    for(unsigned int i = 0; i < fTwoProng_pt.size(); i++)
    {
      TLorentzVector temp_2prong;
      temp_2prong.SetPtEtaPhiM(fTwoProng_pt[i], fTwoProng_eta[i], fTwoProng_phi[i], fTwoProng_mass[i]);
      double temp_dr = z_probetau.DeltaR(temp_2prong);
      if (closest_dr_2prong > temp_dr && temp_dr < 0.4) {
        closest_dr_2prong = temp_dr;
        index_closest_dr_2prong = i;
      }
    }
    int ii = index_closest_dr_2prong;
    if (ii != -1) {
      fprobePassTwoProng = true;
      fProbeTau_twoprong_pt = fTwoProng_pt[ii];
      fProbeTau_twoprong_eta = fTwoProng_eta[ii];
      fProbeTau_twoprong_phi = fTwoProng_phi[ii];
      fProbeTau_twoprong_mass = fTwoProng_mass[ii];
      fProbeTau_twoprong_probeDR = closest_dr_2prong;
      TLorentzVector z_2prong;
      z_2prong.SetPtEtaPhiM(fTwoProng_pt[ii], fTwoProng_eta[ii], fTwoProng_phi[ii], fTwoProng_mass[ii]);
      FillRecoDiObjectInfo(fMuonTwoProng, z_muon, z_2prong);
    } else {
      fprobePassTwoProng = false;
      fProbeTau_twoprong_pt = -999.9;
      fProbeTau_twoprong_eta = -999.9;
      fProbeTau_twoprong_phi = -999.9;
      fProbeTau_twoprong_mass = -999.9;
      fProbeTau_twoprong_probeDR = -999.9;
    }
    
    // visible Z with pattau
    double closest_dr_tau = 100.0;
    unsigned int index_closest_dr_tau = -1;
    for(unsigned int i = 0; i < fTau_pt.size(); i++)
    {
      TLorentzVector temp_tau;
      temp_tau.SetPtEtaPhiM(fTau_pt[i], fTau_eta[i], fTau_phi[i], fTau_mass[i]);
      double temp_dr = z_probetau.DeltaR(temp_tau);
      if (closest_dr_tau > temp_dr && temp_dr < 0.4) {
        closest_dr_tau = temp_dr;
        index_closest_dr_tau = i;
      }
    }
    int jj = index_closest_dr_tau;
    if (jj != -1) {
      fprobePassTauID = true;
      fProbeTau_tauID_pt = fTau_pt[jj];
      fProbeTau_tauID_eta = fTau_eta[jj];
      fProbeTau_tauID_phi = fTau_phi[jj];
      fProbeTau_tauID_mass = fTau_mass[jj];
      fProbeTau_tauID_probeDR = closest_dr_tau;
      TLorentzVector z_pattau;
      z_pattau.SetPtEtaPhiM(fTau_pt[jj], fTau_eta[jj], fTau_phi[jj], fTau_mass[jj]);
      FillRecoDiObjectInfo(fMuonTauID, z_muon, z_pattau);
    } else {
      fprobePassTauID = false;
      fProbeTau_tauID_pt = -999.9;
      fProbeTau_tauID_eta = -999.9;
      fProbeTau_tauID_phi = -999.9;
      fProbeTau_tauID_mass = -999.9;
      fProbeTau_tauID_probeDR = -999.9;
    }
  }

  InitRecoDiObjectInfo(fMuonTwoProng_pt);
  InitRecoDiObjectInfo(fMuonTwoProng_zmass);
  InitRecoDiObjectInfo(fMuonTauID_pt);
  InitRecoDiObjectInfo(fMuonTauID_zmass);
  if (result.passMuonTauPair) {
    TLorentzVector the_muon;
    the_muon.SetPtEtaPhiM(result.tagMuon->pt(), result.tagMuon->eta(), result.tagMuon->phi(), result.tagMuon->mass());
    // form Z with pat::tau
    int best_tau_pt_i = -1;
    int best_tau_mass_i = -1;
    for(unsigned int i = 0; i < fTau_pt.size(); i++)
    {
      TLorentzVector temp_tau;
      temp_tau.SetPtEtaPhiM(fTau_pt[i], fTau_eta[i], fTau_phi[i], fTau_mass[i]);
      if (best_tau_pt_i == -1 || fTau_pt[i] > fTau_pt[best_tau_pt_i]) {
        best_tau_pt_i = i;
      }
      TLorentzVector temp_z;
      temp_z = temp_tau;
      temp_z += the_muon;
      if (best_tau_mass_i != -1) { 
        TLorentzVector other_tau;
        other_tau.SetPtEtaPhiM(fTau_pt[best_tau_mass_i], fTau_eta[best_tau_mass_i], fTau_phi[best_tau_mass_i], fTau_mass[best_tau_mass_i]);
        TLorentzVector other_z;
        other_z = other_tau;
        other_z += the_muon;
        if( fabs(temp_z.M()-TauHadFilters::Z_MASS) < fabs(other_z.M()-TauHadFilters::Z_MASS) ) {
          best_tau_mass_i = i;
        }
      }
      if (best_tau_mass_i == -1) { 
        best_tau_mass_i = i;
      }
    }
    // form Z with twoprong
    int best_2p_pt_i = -1;
    int best_2p_mass_i = -1;
    for(unsigned int i = 0; i < fTwoProng_pt.size(); i++)
    {
      TLorentzVector temp_2p;
      temp_2p.SetPtEtaPhiM(fTwoProng_pt[i], fTwoProng_eta[i], fTwoProng_phi[i], fTwoProng_mass[i]);
      if (best_2p_pt_i == -1 || fTwoProng_pt[i] > fTwoProng_pt[best_2p_pt_i]) {
        best_2p_pt_i = i;
      }
      TLorentzVector temp_z;
      temp_z = temp_2p;
      temp_z += the_muon;
      if (best_2p_mass_i != -1) { 
        TLorentzVector other_2p;
        other_2p.SetPtEtaPhiM(fTwoProng_pt[best_2p_mass_i], fTwoProng_eta[best_2p_mass_i], fTwoProng_phi[best_2p_mass_i], fTwoProng_mass[best_2p_mass_i]);
        TLorentzVector other_z;
        other_z = other_2p;
        other_z += the_muon;
        if( fabs(temp_z.M()-TauHadFilters::Z_MASS) < fabs(other_z.M()-TauHadFilters::Z_MASS) ) {
          best_2p_mass_i = i;
        }
      }
      if (best_2p_mass_i == -1) { 
        best_2p_mass_i = i;
      }
    }
    if (best_tau_pt_i != -1) {
      int t = best_tau_pt_i;
      TLorentzVector the_tau;
      the_tau.SetPtEtaPhiM(fTau_pt[t], fTau_eta[t], fTau_phi[t], fTau_mass[t]);
      FillRecoDiObjectInfo(fMuonTauID_pt, the_muon, the_tau);
    }
    if (best_tau_mass_i != -1) {
      int t = best_tau_mass_i;
      TLorentzVector the_tau;
      the_tau.SetPtEtaPhiM(fTau_pt[t], fTau_eta[t], fTau_phi[t], fTau_mass[t]);
      FillRecoDiObjectInfo(fMuonTauID_zmass, the_muon, the_tau);
    }
    if (best_2p_pt_i != -1) {
      int t = best_2p_pt_i;
      TLorentzVector the_tau;
      the_tau.SetPtEtaPhiM(fTwoProng_pt[t], fTwoProng_eta[t], fTwoProng_phi[t], fTwoProng_mass[t]);
      FillRecoDiObjectInfo(fMuonTwoProng_pt, the_muon, the_tau);
    }
    if (best_2p_mass_i != -1) {
      int t = best_2p_mass_i;
      TLorentzVector the_tau;
      the_tau.SetPtEtaPhiM(fTwoProng_pt[t], fTwoProng_eta[t], fTwoProng_phi[t], fTwoProng_mass[t]);
      FillRecoDiObjectInfo(fMuonTwoProng_zmass, the_muon, the_tau);
    }
  }

  if (fDebug) cout << ". doing mumu preselection branches" << endl;
  TauHadFilters::DiMuonPreSelectionResult dimuon_result;
  dimuon_result = TauHadFilters::computeDiMuonPreSelectionResult(iEvent, triggerBits, triggerObjects, triggerPrescales, vertices, taus, muons, electrons, ak4jets, MET, rhoH, static_cast<TauHadFilters::muonIDtype>(fmuonIDtype), static_cast<TauHadFilters::muonISOtype>(fmuonISOtype));
  bool passDiMuon = (dimuon_result.tagMuon != NULL && dimuon_result.tagMuon2 != NULL);

  fpassMuonTrigger = dimuon_result.passTrigger;
  fpassMuonTriggerTk = dimuon_result.passTriggerTk;
  fMuonTrigger = dimuon_result.foundTrigger;
  fMuonTriggerTk = dimuon_result.foundTriggerTk;
  fnTagMuons = dimuon_result.nTagMuons;
  fpassExtraElectronVeto = dimuon_result.passExtraElectronVeto;
  fpassExtraMuonVeto = dimuon_result.passExtraMuonVeto;
  fpassPreselection = dimuon_result.passPreSelection;
  fpassReducedSelection = passDiMuon;

  InitRecoDiObjectInfo(fMuonMuon);
  if (passDiMuon) {
    fTagMuon1_pt = dimuon_result.tagMuon->pt();
    fTagMuon1_eta = dimuon_result.tagMuon->eta();
    fTagMuon1_phi = dimuon_result.tagMuon->phi();
    fTagMuon1_mass = dimuon_result.tagMuon->mass();
    fTagMuon1_dz =  fabs( (dimuon_result.tagMuon->muonBestTrack())->dz( PV.position() ) );
    fTagMuon1_iso = TauHadFilters::computeMuonIsolation(dimuon_result.tagMuon);

    fTagMuon2_pt = dimuon_result.tagMuon2->pt();
    fTagMuon2_eta = dimuon_result.tagMuon2->eta();
    fTagMuon2_phi = dimuon_result.tagMuon2->phi();
    fTagMuon2_mass = dimuon_result.tagMuon2->mass();
    fTagMuon2_dz =  fabs( (dimuon_result.tagMuon2->muonBestTrack())->dz( PV.position() ) );
    fTagMuon2_iso = TauHadFilters::computeMuonIsolation(dimuon_result.tagMuon2);

    TLorentzVector muon1; muon1.SetPtEtaPhiM(dimuon_result.tagMuon->pt(), dimuon_result.tagMuon->eta(), dimuon_result.tagMuon->phi(), dimuon_result.tagMuon->mass());
    TLorentzVector muon2; muon2.SetPtEtaPhiM(dimuon_result.tagMuon2->pt(), dimuon_result.tagMuon2->eta(), dimuon_result.tagMuon2->phi(), dimuon_result.tagMuon2->mass());
    FillRecoDiObjectInfo(fMuonMuon, muon1, muon2);
  } else {
    fTagMuon1_pt = -999.9;
    fTagMuon1_eta = -999.9;
    fTagMuon1_phi = -999.9;
    fTagMuon1_mass = -999.9;
    fTagMuon1_dz = -999.9;
    fTagMuon1_iso = -999.9;

    fTagMuon2_pt = -999.9;
    fTagMuon2_eta = -999.9;
    fTagMuon2_phi = -999.9;
    fTagMuon2_mass = -999.9;
    fTagMuon2_dz = -999.9;
    fTagMuon2_iso = -999.9;
  }

  if (fincludeLeptonBranches) {
    const double MUON_MIN_PT = 26;
    const double MUON_MAX_ETA = 2.1;
    const double MUON_MAX_DZ = 0.2;
    const double MUON_MAX_DXY = 0.045;

    const double MUON_VLOOSE_RELISO = 0.4;
    const double MUON_LOOSE_RELISO = 0.25;
    const double MUON_MEDIUM_RELISO = 0.20;
    const double MUON_TIGHT_RELISO = 0.15;
    const double MUON_VTIGHT_RELISO = 0.10;
    const double MUON_VVTIGHT_RELISO = 0.05;

    // Tight Muons
    fNum_Muons = 0;
    TLorentzVector fMuon_vector;
    std::vector<const pat::Muon *> passedMuons;
    for (const pat::Muon &mu : *muons) {

      if (mu.pt() > MUON_MIN_PT && fabs(mu.eta()) < MUON_MAX_ETA){
        fNum_Muons++;
      }

      if (mu.pt() > MUON_MIN_PT &&
        fabs(mu.eta()) < MUON_MAX_ETA &&
        TauHadFilters::computeMuonIsolation(&mu) < MUON_TIGHT_RELISO &&
        fabs(mu.muonBestTrack()->dz(PV.position())) < MUON_MAX_DZ &&
        fabs(mu.muonBestTrack()->dxy(PV.position())) < MUON_MAX_DXY) {
          if (!mu.isTightMuon(PV)) continue;
          passedMuons.push_back(&mu);
          fTightMuon_pt.push_back(mu.pt());
          fTightMuon_eta.push_back(mu.eta());
          fTightMuon_phi.push_back(mu.phi());
          fTightMuon_mass.push_back(mu.mass());
        }
    }
    fMuon_veto = 1;
    if (fTightMuon_pt.size() > 0) fMuon_veto = 0;
    fNTightMuons = fTightMuon_pt.size();

    // B Veto
    fbtag_veto = 0;
    for (const pat::Jet &j : *ak4jets) {
        if (j.pt() < 20) continue;
        if (fabs(j.eta()) > 2.5) continue;
        if ( std::max(0.f,j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")) > 0.890) {
          fbtag_veto = 1;
          break;
        }
    }

    //CONSTRUCTING W
    fImag_W = 0;
    fmT = -1.0;
    double w_mass = 80.3790;
    if (fTightMuon_pt.size() > 0) {
      double v_pt, v_eta, v_phi;
      double v_pz = 0;
      double v_mass = 0.0;
      v_pt = met.pt();
      v_phi = met.phi();

      const pat::Muon &mu = *(passedMuons[0]);
      double u_pt, u_eta, u_phi, u_pz, u_mass;
      u_pt = mu.pt();
      u_eta = mu.eta();
      u_phi = mu.phi();
      u_mass = mu.mass();
      u_pz = u_pt * TMath::SinH(u_eta);

      double k = w_mass * w_mass / 2 + u_pt * v_pt * TMath::Cos(u_phi - v_phi);
      double u_p2 = u_pt * u_pt + u_pz * u_pz;

      //FROM TOP PAIR PRODUCTION ARTICLE
      double discr = 4 * k * k * u_pz * u_pz - 4 * u_pt * u_pt * (v_pt * v_pt * u_p2 - k * k);
      if (discr >= 0.0 && k * u_pz > 0) v_pz = (2 * k * u_pz - TMath::Sqrt(discr)) / (2 * u_pt * u_pt);
      else if (discr >= 0.0 && k * u_pz < 0) v_pz = (2 * k * u_pz + TMath::Sqrt(discr)) / (2 * u_pt * u_pt);
      else if (discr < 0.0) {        
        fImag_W = 1;
        v_pz = (2 * k * u_pz) / (2 * u_pt * u_pt);
      }
      double v_p2 = v_pz * v_pz + v_pt * v_pt;
      v_eta = TMath::ATanH( v_pz / TMath::Sqrt(v_p2));

      TLorentzVector w_vector, u_vector, v_vector;
      u_vector.SetPtEtaPhiM(u_pt, u_eta, u_phi, u_mass);
      v_vector.SetPtEtaPhiM(v_pt, v_eta, v_phi, v_mass);
      w_vector = u_vector + v_vector;

      //fmT = TMath::Sqrt(2 * fMuon_pt[0] * met.pt() * (1 - TMath::Cos((*passedMuons[0]).phi() - met.phi())));
      fmT = TMath::Sqrt(2 * u_pt * v_pt * (1 - TMath::Cos(u_phi - v_phi)));
      fW_mT.push_back(fmT);

      fW_pt.push_back(w_vector.Pt());
      fW_eta.push_back(w_vector.Eta());
      fW_phi.push_back(w_vector.Phi());
      fW_mass.push_back(w_vector.M());
    } // end conditional on number of tight muons
  } // end conditional to include lepton+jets branches

  // HT
  fHT = 0.0;
  fST = 0.0;
  fHT_naive = 0.0;
  fHT_qcd = 0.0;
  fHT_l = 0.0;
  fST_l = 0.0;
  for (unsigned int i = 0; i < ak4jets->size(); i++) {
    const pat::Jet &jet = (*ak4jets)[i];
    if (jet.pt() < 30) continue;
    if (fabs(jet.eta()) > 2.5) continue;
    TLorentzVector jet_vector; jet_vector.SetPtEtaPhiE(jet.pt(), jet.eta(), jet.phi(), jet.energy());
    TLorentzVector twoprong_vector;
    TLorentzVector photon_vector;
    if (fnTwoProngs>0) twoprong_vector.SetPtEtaPhiM(fTwoProng_pt[0], fTwoProng_eta[0], fTwoProng_phi[0], fTwoProng_mass[0]);
    if (fNumIDPhotons>0) photon_vector.SetPtEtaPhiM(fIDPhoton_pt[0], fIDPhoton_eta[0], fIDPhoton_phi[0], fIDPhoton_mass[0]);

    fHT_naive += jet.pt();
    fHT_qcd += jet.pt();
    fHT += jet.pt();
    fHT_l += jet.pt();

    // clean HT of twoprong and photon
    if (fnTwoProngs>0 && jet_vector.DeltaR(twoprong_vector) < 0.3) fHT -= jet.pt();
    else if (fNumIDPhotons>0 && jet_vector.DeltaR(photon_vector) < 0.3) fHT -= jet.pt();
    
    // clean old HT of all bad energyfraction() jets
    if (jet.muonEnergyFraction() > 0.7 || jet.electronEnergyFraction() > 0.6 || jet.photonEnergyFraction() > 0.6) fHT_qcd -= jet.pt();

    // add photon to ST
    fST = fHT;
    if (fNumIDPhotons>0) fST += fIDPhoton_pt[0];

    if (fincludeLeptonBranches) {

      TLorentzVector muon_vector;
      if (fNTightMuons>0) muon_vector.SetPtEtaPhiM(fTightMuon_pt[0], fTightMuon_eta[0], fTightMuon_phi[0], fTightMuon_mass[0]);

      // clean HT_l of twoprong and muon
      if (fnTwoProngs>0 && jet_vector.DeltaR(twoprong_vector) < 0.3) fHT_l -= jet.pt();
      else if (fNTightMuons>0 && jet_vector.DeltaR(muon_vector) < 0.3) fHT_l -= jet.pt();

      // add muon to ST_l
      fST_l = fHT_l;
      if (fNTightMuons>0) fST_l += fTightMuon_pt[0];
    }
  }

  // Now fill fTree
  cutflow_total++;
  if (fMakeTrees) {
    bool fill = true;
    if (fFilterOnPhoton && fNumIDPhotons==0) fill = false;
    if (fFilterOnTwoProng && fnTwoProngs==0) fill = false;
    if (fFilterOnLepton && fNTightMuons==0) fill = false;
    if (fFilterForABCDStudy && fNumIDPhotons + fNumLoose1IDPhotons + fNumLoose2IDPhotons + fnTwoProngs + fnTwoProngsLoose == 0) fill = false;
    if (fill) {
      cutflow_passFilter++;
      fTree->Fill();
    }
  }

  /* Histogram making */

  if (fincludeDalitzHistos) {
  for (unsigned int i = 0; i < fTwoProng_pt.size(); i++) {
    double norm = fTwoProng_pt[i];
    double ptg = fTwoProng_photon_pt[i] / norm;
    double pt1 = fTwoProng_CHpos_pt[i] / norm;
    double pt2 = fTwoProng_CHneg_pt[i] / norm;
    double high = max(max(ptg,pt1),pt2);
    double low = min(min(ptg,pt1),pt2);
    double mid;
    if(ptg > low && ptg < high) mid = ptg;
    else if(pt1 > low && pt1 < high) mid = pt1;
    else if(pt2 > low && pt2 < high) mid = pt2;
    else mid = low;
    double smaller = min(pt1,pt2);
    double larger = max(pt1,pt2);
    
    fHighvsMid->Fill(mid, high, fMcXS/fMcN);
    fHighvsLow->Fill(low, high, fMcXS/fMcN);
    fMidvsLow->Fill(low, mid, fMcXS/fMcN);
    fPhotonvsLarger->Fill(larger, ptg, fMcXS/fMcN);
    fPhotonvsSmaller->Fill(smaller, ptg, fMcXS/fMcN);

    fPhotonvsPositive->Fill(pt1, ptg, fMcXS/fMcN);
    fPhotonvsNegative->Fill(pt2, ptg, fMcXS/fMcN);
    fPositivevsNegative->Fill(pt2, pt1, fMcXS/fMcN);

    fPhotonFraction->Fill(ptg, fMcXS/fMcN);
    fPositiveFraction->Fill(pt1, fMcXS/fMcN);
    fNegativeFraction->Fill(pt2, fMcXS/fMcN);
    fHTverify->Fill(fHT, fMcXS/fMcN);
  }
  }
}

void 
TwoProngAnalyzer::beginJob()
{
  if (fMakeTrees) {
  cout << "\n===========================" << endl;
  cout << "= Ntuplizer Configuration =" << endl;
  cout << "===========================" << endl;
  cout << "makeTrees " << fMakeTrees << endl;
  cout << "debug " << fDebug << endl;
  cout << "mcXS " << fMcXS << endl;
  cout << "mcN " << fMcN << endl;
  cout << "filterOnPhoton " << fFilterOnPhoton << endl;
  cout << "filterOnTwoProng " << fFilterOnTwoProng << endl;
  cout << "filterOnLepton " << fFilterOnLepton << endl;
  cout << "filterForABCDStudy " << fFilterForABCDStudy << endl;
  cout << "stackedDalitzHistos " << fincludeDalitzHistos << endl;
  cout << "oldData " << fOldData << endl;
  cout << "===========================" << endl;
  cout << "noTwoProng " << fdontIncludeTwoProngs << endl;
  cout << "includeAllLooseObjects " << fincludeLooseTwoProngs << endl;
  cout << "includeAllCandObjects " << fincludeCandTwoProngs << endl;
  cout << "includeAsymTwoProngs " << fincludeAsymTwoProngs << endl;
  cout << "includeMCInfo " << fincludeMCInfo << endl;
  cout << "includeSignalGenParticles " << fincludeSignalGenParticles << endl;
  cout << "includeOldPhotons " << fincludeOldPhotons << endl;
  cout << "includeBasePhotons " << fincludeBasePhotons << endl;
  cout << "includeConeHEPhotons " << fincludeConeHEPhotons << endl;
  cout << "includeLoosePhotons " << fincludeLoosePhotons << endl;
  cout << "inlcudeZDecayGenParticles " << fincludeZDecayGenParticles << endl;
  cout << "includeTauTauBranches " << fincludeZTauHadBranches << endl;
  cout << "includeMuMuBranches " << fincludeZMuMuBranches << endl;
  cout << "includeLeptonBranches " << fincludeLeptonBranches << endl;
  cout << "usePatTauInPreselection " << fusePatTauForZPreBranches << endl;
  cout << "muonIDtype " << fmuonIDtype << endl;
  cout << "muonISOtype " << fmuonISOtype << endl;
  cout << "===========================" << endl;
  cout << "twoprong_DR " << ftwoprong_DR << endl;
  cout << "twoprong_tracksMinPt " << ftwoprong_tracksMinPt << endl;
  cout << "twoprong_IsolationDR " << ftwoprong_IsolationDR << endl;
  cout << "twoprong_PhiBox " << ftwoprong_PhiBox << endl;
  cout << "twoprong_EtaBox " << ftwoprong_EtaBox << endl;
  cout << "twoprong_PhotonPtCut " << ftwoprong_PhotonPtCut << endl;
  cout << "twoprong_ChargedIsoCut " << ftwoprong_ChargedIsoCut << endl;
  cout << "twoprong_NeutralIsoCut " << ftwoprong_NeutralIsoCut << endl;
  cout << "twoprong_EGammaIsoCut " << ftwoprong_EGammaIsoCut << endl;
  cout << "twoprong_ChargedIsoFakeCut " << ftwoprong_ChargedIsoFakeCut << endl;
  cout << "twoprong_NeutralIsoFakeCut " << ftwoprong_NeutralIsoFakeCut << endl;
  cout << "twoprong_EGammaIsoFakeCut " << ftwoprong_EGammaIsoFakeCut << endl;
  cout << "twoprong_GenMatchDR " << ftwoprong_GenMatchDR << endl;
  cout << "twoprong_AbsMaxEta " << ftwoprong_AbsMaxEta << endl;
  cout << "twoprong_MinPt " << ftwoprong_MinPt << endl;
  cout << "twoprong_TrackAsymmetryCut " << ftwoprong_TrackAsymmetryCut << endl;
  cout << "twoprong_PhotonAsymmetryCut " << ftwoprong_PhotonAsymmetryCut << endl;
  cout << "twoprong_OptionalExtraTrack " << ftwoprong_OptionalExtraTrack << endl;
  cout << "twoprong_FlipAsymReq " << ftwoprong_FlipAsymReq << endl;
  cout << "includeDalitzVariables " << fincludeDalitzVariables << endl;
  cout << "===========================" << endl;
  }
}

void 
TwoProngAnalyzer::endJob()
{
  // print a cutflow if using a filter 
  if ( fFilterOnPhoton || fFilterOnTwoProng || fFilterForABCDStudy || fFilterOnLepton) {
    cout << "\nCutflow report" << endl;
    cout << "==============" << endl;
    cout << "Total_Events " << cutflow_total << endl;
    cout << "Passed_Filter " << cutflow_passFilter << endl;
    cout << "==============" << endl;
  }

  // stacked dalitz histograms
  if (fincludeDalitzHistos) {

  fTripleStackedDalitz->Add(fHighvsMid);
  fTripleStackedDalitz->Add(fHighvsLow);
  fTripleStackedDalitz->Add(fMidvsLow);
  fDoubleStackedDalitz->Add(fPhotonvsLarger);
  fDoubleStackedDalitz->Add(fPhotonvsSmaller);

  }
}

// High-pt Id subroutines
bool TwoProngAnalyzer::photon_isSaturated(const pat::Photon *photon, const EcalRecHitCollection *recHitsEB, const EcalRecHitCollection *recHitsEE,
		                                         const CaloSubdetectorTopology* subDetTopologyEB_, const CaloSubdetectorTopology* subDetTopologyEE_)
{
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
          if (it->checkFlag(EcalRecHit::kSaturated) && !it->checkFlag(EcalRecHit::kDead) && !it->checkFlag(EcalRecHit::kKilled)) {
            isSat = true;
          }
        }
      }
    }
  }
  return isSat;
}

bool TwoProngAnalyzer::photon_passHighPtID(const pat::Photon* photon, double rho, bool isSat)
{
  if (
    photon_passHE(photon) &&
    photon_passIsoCh(photon) &&
    photon_passSigmaIetaIeta(photon,isSat) &&
    photon_passIsoGamma(photon,rho) &&
    photon->passElectronVeto()
  ) return true;

  else return false;
}

bool TwoProngAnalyzer::photon_passHighPtID_loose1(const pat::Photon* photon, double rho, bool isSat)
{
  if (
    photon_passHE(photon) &&
    photon_passIsoCh(photon) &&
    photon_passSigmaIetaIeta(photon,isSat) &&
    photon->passElectronVeto()
  ) return true;

  else return false;
}

bool TwoProngAnalyzer::photon_passHighPtID_loose2(const pat::Photon* photon, double rho, bool isSat)
{
  if (
    photon_passHE(photon) &&
    photon_passSigmaIetaIeta(photon,isSat) &&
    photon_passIsoGamma(photon,rho) &&
    photon->passElectronVeto()
  ) return true;

  else return false;
}

bool TwoProngAnalyzer::photon_passHighPtID_loose3(const pat::Photon* photon, double rho, bool isSat)
{
  if (
    photon_passIsoCh(photon) &&
    photon_passSigmaIetaIeta(photon,isSat) &&
    photon_passIsoGamma(photon,rho) &&
    photon->passElectronVeto()
  ) return true;

  else return false;
}

bool TwoProngAnalyzer::photon_passHighPtID_loose4(const pat::Photon* photon, double rho, bool isSat)
{
  if (
    photon_passHE(photon) &&
    photon_passIsoCh(photon) &&
    photon_passIsoGamma(photon,rho) &&
    photon->passElectronVeto()
  ) return true;

  else return false;
}

bool TwoProngAnalyzer::photon_passHighPtID_loose5(const pat::Photon* photon, double rho, bool isSat)
{
  if (
    photon_passHE(photon) &&
    photon_passIsoCh(photon) &&
    photon_passSigmaIetaIeta(photon,isSat) &&
    photon_passIsoGamma(photon,rho)
  ) return true;

  else return false;
}

bool TwoProngAnalyzer::photon_passHighPtID_base(const pat::Photon* photon, double rho, bool isSat)
{
  if (
    photon->passElectronVeto()
  ) return true;

  else return false;
}

bool TwoProngAnalyzer::photon_passHighPtID_coneHE(const pat::Photon* photon, double rho, bool isSat)
{
  if (
//    photon_passHE(photon) &&
    photon_passHE_coneBased(photon) &&
    photon_passIsoCh(photon) &&
    photon_passSigmaIetaIeta(photon,isSat) &&
    photon_passIsoGamma(photon,rho) &&
    photon->passElectronVeto()
  ) return true;

  else return false;
}

double TwoProngAnalyzer::photon_computeHE(const pat::Photon * photon)
{
  return photon->hadTowOverEm();
} 

double TwoProngAnalyzer::photon_computeHE_coneBased(const pat::Photon * photon)
{
  return photon->hadronicOverEm();
} 

double TwoProngAnalyzer::photon_computeIsoCh(const pat::Photon * photon)
{
  return photon->chargedHadronIso();
}

double TwoProngAnalyzer::photon_computeIsoGamma(const pat::Photon * photon, double rho)
{
  double phoIso = photon->photonIso();
  return (photon_computeAlpha(photon) + phoIso - rho*photon_computeEA(photon) - photon_computeKappa(photon)*photon->pt());
}

double TwoProngAnalyzer::photon_computeSigmaIetaIeta(const pat::Photon * photon)
{
  return photon->full5x5_sigmaIetaIeta();
}

double TwoProngAnalyzer::photon_scEta(const pat::Photon * photon)
{
  return fabs(photon->superCluster()->eta());
}

bool TwoProngAnalyzer::photon_passHE(const pat::Photon* photon)
{
  double hOverE = photon->hadTowOverEm();
  if (hOverE < 0.05) return true;
  else return false;
}

bool TwoProngAnalyzer::photon_passHE_coneBased(const pat::Photon* photon)
{
  double hOverE = photon->hadronicOverEm();
  if (hOverE < 0.05) return true;
  else return false;
}

bool TwoProngAnalyzer::photon_passIsoCh(const pat::Photon* photon)
{
  double chIsoCut = 5.;
  double chIso = photon->chargedHadronIso();
  if (chIso < chIsoCut) return true;
  else return false;
}

bool TwoProngAnalyzer::photon_passSigmaIetaIeta(const pat::Photon* photon, bool isSaturated)
{
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

bool TwoProngAnalyzer::photon_passIsoGamma(const pat::Photon* photon, double rho)
{
  double phoEta = fabs(photon->superCluster()->eta());
  double corPhoIsoCut = -999.9;
  double corPhoIso = photon_computeIsoGamma(photon,rho);

  if (phoEta < 1.4442) corPhoIsoCut = 2.75;
  if (1.566 < phoEta && phoEta < 2.5) corPhoIsoCut = 2.00;

  if (corPhoIso < corPhoIsoCut) return true;
  else return false;
}

double TwoProngAnalyzer::photon_computeAlpha(const pat::Photon *photon)
{
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

double TwoProngAnalyzer::photon_computeEA(const pat::Photon* photon)
{
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

double TwoProngAnalyzer::photon_computeKappa(const pat::Photon *photon)
{
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

// other routines
double TwoProngAnalyzer::iso_alpha(double eta, int ver)
{
  if (ver == 1)
  {
    if(eta<0.9) return 2.5;
    if(eta>=0.9 && eta<1.4442) return 2.5;
    if(eta>=1.4442 && eta<1.566) return -999.0;
    if(eta>1.566 && eta<2.0) return 2.5;
    if(eta>=2.0 && eta<2.2) return 2.5;
    if(eta>=2.2 && eta<2.5) return 2.5;
  }
  else if (ver == 2)
  {
    if(eta<0.9) return 0.99;
    if(eta>=0.9 && eta<1.4442) return 0.99;
    if(eta>=1.4442 && eta<1.566) return -999.0;
    if(eta>1.566 && eta<2.0) return 0.77;
    if(eta>=2.0 && eta<2.2) return 0.77;
    if(eta>=2.2 && eta<2.5) return 0.77;
  }
  return -999.0;
}

double TwoProngAnalyzer::iso_area(double eta, int ver)
{
  if (ver == 1)
  {
    if(eta<0.9) return 0.17;
    if(eta>=0.9 && eta<1.4442) return 0.14;
    if(eta>=1.4442 && eta<1.566) return -999.0;
    if(eta>1.566 && eta<2.0) return 0.11;
    if(eta>=2.0 && eta<2.2) return 0.14;
    if(eta>=2.2 && eta<2.5) return 0.22;
  }
  else if (ver == 2)
  {
    if(eta<0.9) return 0.15;
    if(eta>=0.9 && eta<1.4442) return 0.13;
    if(eta>=1.4442 && eta<1.566) return -999.0;
    if(eta>1.566 && eta<2.0) return 0.093;
    if(eta>=2.0 && eta<2.2) return 0.15;
    if(eta>=2.2 && eta<2.5) return 0.21;
  }
  return -999.0;
}

double TwoProngAnalyzer::iso_kappa(double eta, int ver)
{
  if (ver == 1)
  {
    if(eta<0.9) return 0.0045;
    if(eta>=0.9 && eta<1.4442) return 0.0045;
    if(eta>=1.4442 && eta<1.566) return -999.0;
    if(eta>1.566 && eta<2.0) return 0.003;
    if(eta>=2.0 && eta<2.2) return 0.003;
    if(eta>=2.2 && eta<2.5) return 0.003;
  }
  else if (ver == 2)
  {
    if(eta<0.9) return 0.0016;
    if(eta>=0.9 && eta<1.4442) return 0.0016;
    if(eta>=1.4442 && eta<1.566) return -999.0;
    if(eta>1.566 && eta<2.0) return 0.00075;
    if(eta>=2.0 && eta<2.2) return 0.00075;
    if(eta>=2.2 && eta<2.5) return 0.00075;
  }
  return -999.0;
}

bool TwoProngAnalyzer::isNeutral(int one, int two, int three)
{
  if (one == 22 && two == 22) return true;
  if (two == 22 && three == 22) return true;
  if (one == 22 && three == 22) return true;
  if (one == 111 && two == 111 && three == 111) return true;

  return false;
}

bool TwoProngAnalyzer::isCharged(int one, int two, int three)
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

bool TwoProngAnalyzer::compareCandsByPt(const edm::Ptr<const reco::Candidate> cand1, const edm::Ptr<const reco::Candidate> cand2)
{
  return(cand1->pt() >= cand2->pt());
}

//define this as a plug-in
DEFINE_FWK_MODULE(TwoProngAnalyzer);
