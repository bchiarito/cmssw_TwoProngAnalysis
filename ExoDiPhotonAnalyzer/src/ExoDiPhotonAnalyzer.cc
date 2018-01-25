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
#include <vector>
#include <algorithm>
#include <cmath>

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

// these objects are all in the namespace 'ExoDiPhotons'
#include "DiPhotonAnalysis/CommonClasses/interface/RecoPhotonInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/RecoDiObjectInfo.h"

// pat objects
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

// for trigger
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

// for gen event info
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

using namespace std;

// forward declare photon functions
double phoKappaHighPtID(const pat::Photon *);
double phoEAHighPtID(const pat::Photon* );
double phoAlphaHighPtID(const pat::Photon *);
bool passCorPhoIsoHighPtID(const pat::Photon* , double );
bool passSigmaIetaIetaCut(const pat::Photon* , bool );
bool passChargedHadronCut(const pat::Photon* );
bool passHadTowerOverEmCut(const pat::Photon*);
bool passHadDrConeOverEmCut(const pat::Photon*);
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
  double iso_alpha(double, int);
  double iso_area(double, int);
  double iso_kappa(double, int);
  vector<string> getDecay(const reco::Candidate &,int flag=0);

  // ----------member data ---------------------------

  // trigger
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;

  // global config file options
  bool               fDebug;                       // if set to False, mean to limit per event stdout output
  bool               fAddDrConePhotonCut;          // option needed for studying the Trigger ID, no longer needd
  bool               fincludeSignalGenParticles; // includes the GenPhi and other gen particles in ntuple
  bool               frunningOnTauTauMC;           // running on Z->tau tau MC, will do gen particle matching to hadronic taus
  bool               fincludeAllLooseObjects;    // include all loose twoprong objects in ntuple
  bool               fincludeAllCandObjects;     // include all twoprong candidate objects (no iso req) in ntuple
  bool               fincludeOldPhotons;         // include the Photon1,Photon2,Photon3 objects in ntuple
  bool               fincludeMCInfo;             // include MC weight in ntuple from GenEventInfo
  double             fMcXS;                      // the mc cross section for scaling purposes
  double             fMcN;                       // the mc number generated for scaling purposes
  bool               fMakeTrees;                // flag to include ttrees in output
  bool               fFakeRateHistos;       // flag to include histograms in output
  bool               fTriggerEffHistos;     // flag to include histograms in output
  bool               fTwoProngYieldHistos;  // flag to include histograms in output
  bool               fStackedDalitzHistos;  // flag to include stacked dalitz plots

  // ntuplizer config file options
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
  double fCandidateAbsMaxEta;
  double fCandidateMinPt;
  double fCandidateTrackAsymmetryCut;
  double fCandidatePhotonAsymmetryCut;
  bool fCandidateOptionalExtraTrack;

  // constants
  double PI0_MASS = 0.135;
  double ETA_MASS = 0.548;

  // EDM Handles
  edm::EDGetTokenT<double> rhoToken_; 
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfcandsToken_;
  edm::EDGetTokenT<vector<reco::GenParticle>> genToken_;
  edm::EDGetTokenT<vector<reco::Vertex>> pvToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamToken_;
  edm::EDGetTokenT<std::vector<pat::Jet>> ak4Token_;
  edm::EDGetTokenT<std::vector<pat::Photon>> photonToken_;
  edm::EDGetTokenT<std::vector<pat::Electron>> electronToken_;
  edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
  edm::EDGetTokenT<std::vector<pat::MET>> metToken_;
  edm::EDGetToken gedphotonsToken_;
  edm::EDGetToken genEventInfoToken_;

  // tools for rechits collection
  std::auto_ptr<noZS::EcalClusterLazyTools> lazyTools_;
  edm::InputTag recHitsEBTag_;
  edm::InputTag recHitsEETag_;
  edm::EDGetTokenT<EcalRecHitCollection> recHitsEBToken;
  edm::EDGetTokenT<EcalRecHitCollection> recHitsEEToken;

  // photon subroutines
  bool photon_isSaturated(const pat::Photon*, const EcalRecHitCollection *, const EcalRecHitCollection *,
                          const CaloSubdetectorTopology*, const CaloSubdetectorTopology*);
  bool photon_passHighPtID(const pat::Photon*, double , bool );
  bool photon_passHighPtID_DrCone(const pat::Photon* , double , bool );
  bool photon_passHighPtID_loose(const pat::Photon* , double , bool );
  bool photon_passHighPtID_base(const pat::Photon* , double , bool );

  // fake rate histos
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

  // trigger eff histos
  TH1F *fPhotonTriggerEff_all_Numerator;
  TH1F *fPhotonTriggerEff_all_Denominator;
  TH1F *fPhotonTriggerEff_all_Division;
  TH1F *fPhotonTriggerEff_Photon175_Numerator;
  TH1F *fPhotonTriggerEff_Photon175_Division;
  TH1F *fPhotonTriggerEff_Photon22_Iso_Numerator;
  TH1F *fPhotonTriggerEff_Photon22_Iso_Division;
  TH1F *fPhotonTriggerEff_ConeHE_all_Numerator;
  TH1F *fPhotonTriggerEff_ConeHE_all_Denominator;
  TH1F *fPhotonTriggerEff_ConeHE_all_Division;
  TH1F *fPhotonTriggerEff_ConeHE_Photon175_Numerator;
  TH1F *fPhotonTriggerEff_ConeHE_Photon175_Division;
  TH1F *fPhotonTriggerEff_ConeHE_Photon22_Iso_Numerator;
  TH1F *fPhotonTriggerEff_ConeHE_Photon22_Iso_Division;

  // twoprong yield histos
  TH1F *fTwoProngYield;

  // stacked dalitz histos
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
 
  // Main Ntuple Ttree and braches
  TTree *fTree2;
  int fHLT_Photon175;
  int fHLT_Photon22_Iso;
  int fEventNum;
  int fRunNum;
  int fLumiNum;
  int fNumPVs;
  double fRho;
  int fNumPF;
  int fNumPrunedPF;
  double fHT;
  double fMET;
  double fMET_phi;
  int fNumElectrons;
  int fNumMuons;
  double fMcW;
  double fMcWProd;

  double fTauDecayType;
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
  vector<Double_t> fGenPhi_px;
  vector<Double_t> fGenPhi_py;
  vector<Double_t> fGenPhi_pz;
  vector<Double_t> fGenPhi_energy;
  vector<Double_t> fGenOmega_pt;
  vector<Double_t> fGenOmega_eta;
  vector<Double_t> fGenOmega_phi;
  vector<Double_t> fGenOmega_mass;
  vector<Double_t> fGenOmega_px;
  vector<Double_t> fGenOmega_py;
  vector<Double_t> fGenOmega_pz;
  vector<Double_t> fGenOmega_energy;

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

  int fNumTightPhotons;
  int fNumTightPhotons_ConeHE;
  vector<Double_t> fBaseIDPhoton_pt;
  vector<Double_t> fBaseIDPhoton_eta;
  vector<Double_t> fBaseIDPhoton_phi;
  vector<Double_t> fBaseIDPhoton_mass;
  vector<Double_t> fBaseIDPhoton_iso_gamma;
  vector<Double_t> fBaseIDPhoton_iso_ch;
  vector<Double_t> fBaseIDPhoton_HE;
  vector<Double_t> fBaseIDPhoton_sigmaieie;
  vector<Double_t> fLooseIDPhoton_pt;
  vector<Double_t> fLooseIDPhoton_eta;
  vector<Double_t> fLooseIDPhoton_phi;
  vector<Double_t> fLooseIDPhoton_mass;
  vector<Double_t> fLooseIDPhoton_iso_gamma;
  vector<Double_t> fIDPhoton_pt;
  vector<Double_t> fIDPhoton_eta;
  vector<Double_t> fIDPhoton_phi;
  vector<Double_t> fIDPhoton_mass;
  vector<Double_t> fID2Photon_pt;
  vector<Double_t> fID2Photon_eta;
  vector<Double_t> fID2Photon_phi;
  vector<Double_t> fID2Photon_mass;
  ExoDiPhotons::recoPhotonInfo_t fRecoTightPhotonInfo1;
  ExoDiPhotons::recoPhotonInfo_t fRecoTightPhotonInfo2;
  ExoDiPhotons::recoPhotonInfo_t fRecoTightPhotonInfo3;

  int fNumTwoProng;
  int fNumTwoProngPass;
  int fNumTwoProngLoose;
  vector<Double_t> fCand_pt;
  vector<Double_t> fCand_eta;
  vector<Double_t> fCand_phi;
  vector<Double_t> fCand_mass;
  vector<Double_t> fCand_mass_l;
  vector<Double_t> fCand_Mass;
  vector<Double_t> fCand_Mass_l;
  vector<Double_t> fCand_MassEta;
  vector<Double_t> fCand_MassEta_l;
  vector<Double_t> fCand_Mass300;
  vector<Bool_t> fCand_FoundExtraTrack;
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
  vector<Double_t> fCand_photon_pt_l;
  vector<Double_t> fCand_photon_eta_l;
  vector<Double_t> fCand_photon_phi_l;
  vector<Double_t> fCand_photon_mass_l;
  vector<Double_t> fCand_photon_Mass;
  vector<Double_t> fCand_photon_nGamma;
  vector<Double_t> fCand_photon_nElectron;
  vector<Double_t> fCand_chargedIso;
  vector<Double_t> fCand_neutralIso;
  vector<Double_t> fCand_egammaIso;
  vector<Double_t> fCand_relchargedIso;
  vector<Double_t> fCand_relneutralIso;
  vector<Double_t> fCand_relegammaIso;
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
  vector<Double_t> fCand_trackAsym;
  vector<Double_t> fCand_photonAsym;
  vector<Double_t> fCand_mPosPho;
  vector<Double_t> fCand_mPosPho_l;
  vector<Double_t> fCand_mPosPho_pi0;
  vector<Double_t> fCand_mPosPho_lpi0;
  vector<Double_t> fCand_mNegPho;
  vector<Double_t> fCand_mNegPho_l;
  vector<Double_t> fCand_mNegPho_pi0;
  vector<Double_t> fCand_mNegPho_lpi0;
  vector<Double_t> fCand_mPosNeg;
  vector<Double_t> fCand_CHpos_p3;
  vector<Double_t> fCand_CHneg_p3;
  vector<Double_t> fCand_photon_p3;
  vector<Int_t> fCand_nChargedIsoCone;
  vector<Int_t> fCand_nNeutralIsoCone;
  vector<Int_t> fCand_nEGammaIsoCone;
  vector<Double_t> fCand_genOmega_dR;
  vector<Double_t> fCand_genTau_dR;
  vector<Bool_t> fCand_tight;
  vector<Bool_t> fCand_passChargedIso;
  vector<Bool_t> fCand_passNeutralIso;
  vector<Bool_t> fCand_passEGammaIso;
  vector<Bool_t> fCand_passPhotonPt;
  vector<Bool_t> fCand_loose;
  vector<Double_t> fCand_iso_gamma;
  vector<Double_t> fCand_iso_gamma_allPV;
  vector<Double_t> fCand_iso_gamma_rel;
  vector<Double_t> fCand_iso_ch;
  vector<Double_t> fCand_iso_ch_allPV;
  vector<Double_t> fCand_iso_ch_rel;
  vector<Double_t> fCand_iso_gammacor1;
  vector<Double_t> fCand_iso_gammacor2;

  vector<Double_t> fTwoProng_pt;
  vector<Double_t> fTwoProng_eta;
  vector<Double_t> fTwoProng_phi;
  vector<Double_t> fTwoProng_mass;
  vector<Double_t> fTwoProng_mass_l;
  vector<Double_t> fTwoProng_Mass;
  vector<Double_t> fTwoProng_Mass_l;
  vector<Double_t> fTwoProng_MassEta;
  vector<Double_t> fTwoProng_MassEta_l;
  vector<Double_t> fTwoProng_Mass300;
  vector<Bool_t> fTwoProng_FoundExtraTrack;
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
  vector<Double_t> fTwoProng_trackAsym;
  vector<Double_t> fTwoProng_photonAsym;
  vector<Double_t> fTwoProng_photon_pt;
  vector<Double_t> fTwoProng_photon_eta;
  vector<Double_t> fTwoProng_photon_phi;
  vector<Double_t> fTwoProng_photon_mass;
  vector<Double_t> fTwoProng_photon_pt_l;
  vector<Double_t> fTwoProng_photon_eta_l;
  vector<Double_t> fTwoProng_photon_phi_l;
  vector<Double_t> fTwoProng_photon_mass_l;
  vector<Double_t> fTwoProng_photon_Mass;
  vector<Double_t> fTwoProng_photon_nGamma;
  vector<Double_t> fTwoProng_photon_nElectron;
  vector<Double_t> fTwoProng_chargedIso;
  vector<Double_t> fTwoProng_neutralIso;
  vector<Double_t> fTwoProng_egammaIso;
  vector<Double_t> fTwoProng_mPosPho;
  vector<Double_t> fTwoProng_mPosPho_l;
  vector<Double_t> fTwoProng_mPosPho_pi0;
  vector<Double_t> fTwoProng_mPosPho_lpi0;
  vector<Double_t> fTwoProng_mNegPho;
  vector<Double_t> fTwoProng_mNegPho_l;
  vector<Double_t> fTwoProng_mNegPho_pi0;
  vector<Double_t> fTwoProng_mNegPho_lpi0;
  vector<Double_t> fTwoProng_mPosNeg;
  vector<Double_t> fTwoProng_CHpos_p3;
  vector<Double_t> fTwoProng_CHneg_p3;
  vector<Double_t> fTwoProng_photon_p3;
  vector<Int_t> fTwoProng_nChargedIsoCone;
  vector<Int_t> fTwoProng_nNeutralIsoCone;
  vector<Int_t> fTwoProng_nEGammaIsoCone;
  vector<Double_t> fTwoProng_genOmega_dR;
  vector<Double_t> fTwoProng_genTau_dR;

  vector<Double_t> fTwoProngLoose_pt;
  vector<Double_t> fTwoProngLoose_eta;
  vector<Double_t> fTwoProngLoose_phi;
  vector<Double_t> fTwoProngLoose_mass;
  vector<Double_t> fTwoProngLoose_mass_l;
  vector<Double_t> fTwoProngLoose_Mass;
  vector<Double_t> fTwoProngLoose_Mass_l;
  vector<Double_t> fTwoProngLoose_MassEta;
  vector<Double_t> fTwoProngLoose_MassEta_l;
  vector<Double_t> fTwoProngLoose_Mass300;
  vector<Bool_t> fTwoProngLoose_FoundExtraTrack;
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
  vector<Double_t> fTwoProngLoose_trackAsym;
  vector<Double_t> fTwoProngLoose_photonAsym;
  vector<Double_t> fTwoProngLoose_photon_pt;
  vector<Double_t> fTwoProngLoose_photon_eta;
  vector<Double_t> fTwoProngLoose_photon_phi;
  vector<Double_t> fTwoProngLoose_photon_mass;
  vector<Double_t> fTwoProngLoose_photon_pt_l;
  vector<Double_t> fTwoProngLoose_photon_eta_l;
  vector<Double_t> fTwoProngLoose_photon_phi_l;
  vector<Double_t> fTwoProngLoose_photon_mass_l;
  vector<Double_t> fTwoProngLoose_photon_Mass;
  vector<Double_t> fTwoProngLoose_photon_nGamma;
  vector<Double_t> fTwoProngLoose_photon_nElectron;
  vector<Double_t> fTwoProngLoose_chargedIso;
  vector<Double_t> fTwoProngLoose_neutralIso;
  vector<Double_t> fTwoProngLoose_egammaIso;
  vector<Double_t> fTwoProngLoose_mPosPho;
  vector<Double_t> fTwoProngLoose_mPosPho_l;
  vector<Double_t> fTwoProngLoose_mPosPho_pi0;
  vector<Double_t> fTwoProngLoose_mPosPho_lpi0;
  vector<Double_t> fTwoProngLoose_mNegPho;
  vector<Double_t> fTwoProngLoose_mNegPho_l;
  vector<Double_t> fTwoProngLoose_mNegPho_pi0;
  vector<Double_t> fTwoProngLoose_mNegPho_lpi0;
  vector<Double_t> fTwoProngLoose_mPosNeg;
  vector<Int_t> fTwoProngLoose_nChargedIsoCone;
  vector<Int_t> fTwoProngLoose_nNeutralIsoCone;
  vector<Int_t> fTwoProngLoose_nEGammaIsoCone;
  vector<Double_t> fTwoProngLoose_genOmega_dR;
  vector<Double_t> fTwoProngLoose_genTau_dR;

  ExoDiPhotons::recoDiObjectInfo_t fRecoPhiDiTwoProng; 
  ExoDiPhotons::recoDiObjectInfo_t fRecoPhiPhotonTwoProng; 
  ExoDiPhotons::recoDiObjectInfo_t fRecoPhiInclusive; 
};

//
// constructors and destructor
//
ExoDiPhotonAnalyzer::ExoDiPhotonAnalyzer(const edm::ParameterSet& iConfig)
  : 
    triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
    triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
    triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
    fDebug(iConfig.getUntrackedParameter<bool>("debug")),
    fAddDrConePhotonCut(iConfig.getUntrackedParameter<bool>("addPhotonCutDrConeHE")),
    fincludeSignalGenParticles(iConfig.getUntrackedParameter<bool>("includeSignalGenParticles")),
    frunningOnTauTauMC(iConfig.getUntrackedParameter<bool>("runningOnTauTauMC")),
    fincludeAllLooseObjects(iConfig.getUntrackedParameter<bool>("includeAllLooseObjects")),
    fincludeAllCandObjects(iConfig.getUntrackedParameter<bool>("includeAllCandObjects")),
    fincludeOldPhotons(iConfig.getUntrackedParameter<bool>("includeOldPhotons")),
    fincludeMCInfo(iConfig.getUntrackedParameter<bool>("includeMCInfo")),
    fMcXS(iConfig.getUntrackedParameter<double>("mcXS")),
    fMcN(iConfig.getUntrackedParameter<double>("mcN")),
    fMakeTrees(iConfig.getUntrackedParameter<bool>("makeTrees")),
    fFakeRateHistos(iConfig.getUntrackedParameter<bool>("fakeRateHistos")),
    fTriggerEffHistos(iConfig.getUntrackedParameter<bool>("triggerEffHistos")),
    fTwoProngYieldHistos(iConfig.getUntrackedParameter<bool>("twoprongYieldHistos")),
    fStackedDalitzHistos(iConfig.getUntrackedParameter<bool>("stackedDalitzHistos")),
    fCandidatePairDR(iConfig.getUntrackedParameter<double>("chargedHadronPairMinDR")),
    fCandidatePairMinPt(iConfig.getUntrackedParameter<double>("chargedHadronMinPt")),
    fCandidatePairIsolationDR(iConfig.getUntrackedParameter<double>("isolationConeR")),
    fCandidatePairPhiBox(iConfig.getUntrackedParameter<double>("photonPhiBoxSize")),
    fCandidatePairEtaBox(iConfig.getUntrackedParameter<double>("photonEtaBoxSize")),
    fCandidatePairPhotonPtCut(iConfig.getUntrackedParameter<double>("photonPtCut")),
    fCandidatePairChargedIsoCut(iConfig.getUntrackedParameter<double>("chargedIsoCut")),
    fCandidatePairNeutralIsoCut(iConfig.getUntrackedParameter<double>("neutralIsoCut")),
    fCandidatePairEGammaIsoCut(iConfig.getUntrackedParameter<double>("egammaIsoCut")),
    fCandidatePairChargedIsoFakeCut(iConfig.getUntrackedParameter<double>("chargedIsoLooseMax")),
    fCandidatePairNeutralIsoFakeCut(iConfig.getUntrackedParameter<double>("neutralIsoLooseMax")),
    fCandidatePairEGammaIsoFakeCut(iConfig.getUntrackedParameter<double>("egammaIsoLooseMax")),
    fCandidatePairGenMatchDR(iConfig.getUntrackedParameter<double>("generatorMatchDR")),
    fCandidateAbsMaxEta(iConfig.getUntrackedParameter<double>("candidateAbsMaxEta")),
    fCandidateMinPt(iConfig.getUntrackedParameter<double>("candidateMinPt")),
    fCandidateTrackAsymmetryCut(iConfig.getUntrackedParameter<double>("candidateTrackAsymmetryCut")),
    fCandidatePhotonAsymmetryCut(iConfig.getUntrackedParameter<double>("candidatePhotonAsymmetryCut")),
    fCandidateOptionalExtraTrack(iConfig.getUntrackedParameter<bool>("candidateOptionalExtraTrack")),
    rhoToken_(consumes<double> (iConfig.getParameter<edm::InputTag>("rho")))
{
  // MiniAOD event content
  pfcandsToken_ = consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates"));
  genToken_ = consumes<vector<reco::GenParticle>>(edm::InputTag("prunedGenParticles"));
  pvToken_ = consumes<vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));
  beamToken_ = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  ak4Token_ = consumes<std::vector<pat::Jet>>(edm::InputTag("slimmedJets"));
  photonToken_ = consumes<std::vector<pat::Photon>>(edm::InputTag("slimmedPhotons"));
  metToken_ = consumes<std::vector<pat::MET>>(edm::InputTag("slimmedMETs"));
  electronToken_ = consumes<std::vector<pat::Electron>>(edm::InputTag("slimmedElectrons"));
  muonToken_ = consumes<std::vector<pat::Muon>>(edm::InputTag("slimmedMuons"));
  gedphotonsToken_ = consumes<edm::View<pat::Photon>>( edm::InputTag("slimmedPhotons") );
  genEventInfoToken_ = mayConsume<GenEventInfoProduct>( edm::InputTag("generator") );

  recHitsEBTag_ = iConfig.getUntrackedParameter<edm::InputTag>("RecHitsEBTag",edm::InputTag("reducedEgamma:reducedEBRecHits"));
  recHitsEETag_ = iConfig.getUntrackedParameter<edm::InputTag>("RecHitsEETag",edm::InputTag("reducedEgamma:reducedEERecHits"));
  recHitsEBToken = consumes < EcalRecHitCollection > (recHitsEBTag_);
  recHitsEEToken = consumes < EcalRecHitCollection > (recHitsEETag_);

  // output setup
  edm::Service<TFileService> fs;
  // Branches for charged decay analysis
  if (fMakeTrees) {
  fTree2 = fs->make<TTree>("fTree2","ChargedDecayTree");
  // Trigger
  fTree2->Branch("HLT_Photon175",&fHLT_Photon175,"HLT_Photon175/I");
  fTree2->Branch("HLT_Photon22_Iso",&fHLT_Photon22_Iso,"HLT_Photon22_Iso/I");
  // Event wide
  fTree2->Branch("eventNum",&fEventNum,"eventNum/I");
  fTree2->Branch("runNum",&fRunNum,"runNum/I");
  fTree2->Branch("lumiNum",&fLumiNum,"lumiNum/I");
  fTree2->Branch("nPV",&fNumPVs,"nPV/I");
  fTree2->Branch("rho",&fRho,"rho/D");
  fTree2->Branch("nPF",&fNumPF,"nPF/I");
  fTree2->Branch("nPrunedPF",&fNumPrunedPF,"numPrunedPF/I");
  fTree2->Branch("HT",&fHT,"HT/D");
  fTree2->Branch("MET",&fMET,"MET/D");
  fTree2->Branch("MET_phi",&fMET_phi,"MET_phi/D");
  fTree2->Branch("nElectrons",&fNumElectrons,"nElectrons/I");
  fTree2->Branch("nMuons",&fNumMuons,"nMuons/I");
  fTree2->Branch("mcW",&fMcW,"mcW/D");
  fTree2->Branch("mcWProd",&fMcWProd,"mcWProd/D");
  fTree2->Branch("mcXS",&fMcXS,"mcXS/D");
  fTree2->Branch("mcN",&fMcN,"mcN/D");
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
  // TwoProngs
  fTree2->Branch("nTwoProngCands",&fNumTwoProng,"nTwoProngCands/I");
  fTree2->Branch("nTwoProngs",&fNumTwoProngPass,"nTwoProngs/I");
  fTree2->Branch("nTwoProngsLoose",&fNumTwoProngLoose,"nTwoProngsLoose/I");
  if(fincludeAllCandObjects) {
  // Candidate information
  fTree2->Branch("Cand_pt",&fCand_pt);
  fTree2->Branch("Cand_eta",&fCand_eta);
  fTree2->Branch("Cand_phi",&fCand_phi);
  fTree2->Branch("Cand_mass",&fCand_mass);
  fTree2->Branch("Cand_mass_l",&fCand_mass_l);
  fTree2->Branch("Cand_Mass",&fCand_Mass);
  fTree2->Branch("Cand_Mass_l",&fCand_Mass_l);
  fTree2->Branch("Cand_MassEta",&fCand_MassEta);
  fTree2->Branch("Cand_MassEta_l",&fCand_MassEta_l);
  fTree2->Branch("Cand_Mass300",&fCand_Mass300);
  fTree2->Branch("Cand_FoundExtraTrack",&fCand_FoundExtraTrack);
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
  fTree2->Branch("Cand_photon_pt_l",&fCand_photon_pt_l);
  fTree2->Branch("Cand_photon_eta_l",&fCand_photon_eta_l);
  fTree2->Branch("Cand_photon_phi_l",&fCand_photon_phi_l);
  fTree2->Branch("Cand_photon_mass_l",&fCand_photon_mass_l);
  fTree2->Branch("Cand_photon_Mass",&fCand_photon_Mass);
  fTree2->Branch("Cand_photon_nGamma",&fCand_photon_nGamma);
  fTree2->Branch("Cand_photon_nElectron",&fCand_photon_nElectron);
  fTree2->Branch("Cand_chargedIso",&fCand_chargedIso);
  fTree2->Branch("Cand_neutralIso",&fCand_neutralIso);
  fTree2->Branch("Cand_egammaIso",&fCand_egammaIso);
  fTree2->Branch("Cand_relchargedIso",&fCand_relchargedIso);
  fTree2->Branch("Cand_relneutralIso",&fCand_relneutralIso);
  fTree2->Branch("Cand_relegammaIso",&fCand_relegammaIso);
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
  fTree2->Branch("Cand_trackAsym",&fCand_trackAsym);
  fTree2->Branch("Cand_photonAsym",&fCand_photonAsym);
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
  fTree2->Branch("Cand_genOmega_dR",&fCand_genOmega_dR);
  fTree2->Branch("Cand_genTau_dR",&fCand_genTau_dR);
  fTree2->Branch("Cand_pass",&fCand_tight);
  fTree2->Branch("Cand_passChargedIso",&fCand_passChargedIso);
  fTree2->Branch("Cand_passNeutralIso",&fCand_passNeutralIso);
  fTree2->Branch("Cand_passEGammaIso",&fCand_passEGammaIso);
  fTree2->Branch("Cand_passPhotonPt",&fCand_passPhotonPt);
  fTree2->Branch("Cand_mPosPho",&fCand_mPosPho);
  fTree2->Branch("Cand_mPosPho_l",&fCand_mPosPho_l);
  fTree2->Branch("Cand_mPosPho_pi0",&fCand_mPosPho_pi0);
  fTree2->Branch("Cand_mPosPho_lpi0",&fCand_mPosPho_lpi0);
  fTree2->Branch("Cand_mNegPho",&fCand_mNegPho);
  fTree2->Branch("Cand_mNegPho_l",&fCand_mNegPho_l);
  fTree2->Branch("Cand_mNegPho_pi0",&fCand_mNegPho_pi0);
  fTree2->Branch("Cand_mNegPho_lpi0",&fCand_mNegPho_lpi0);
  fTree2->Branch("Cand_mPosNeg",&fCand_mPosNeg);
  fTree2->Branch("Cand_CHpos_p3",&fCand_CHpos_p3);
  fTree2->Branch("Cand_CHneg_p3",&fCand_CHneg_p3);
  fTree2->Branch("Cand_photon_p3",&fCand_photon_p3);
  fTree2->Branch("Cand_fake",&fCand_loose);
  fTree2->Branch("Cand_iso_gamma",&fCand_iso_gamma);
  fTree2->Branch("Cand_iso_gamma_allPV",&fCand_iso_gamma_allPV);
  fTree2->Branch("Cand_iso_gamma_rel",&fCand_iso_gamma_rel);
  fTree2->Branch("Cand_iso_ch",&fCand_iso_ch);
  fTree2->Branch("Cand_iso_ch_allPV",&fCand_iso_ch_allPV);
  fTree2->Branch("Cand_iso_ch_rel",&fCand_iso_ch_rel);
  fTree2->Branch("Cand_iso_gammacor1",&fCand_iso_gammacor1);
  fTree2->Branch("Cand_iso_gammacor2",&fCand_iso_gammacor2);
  }
  if(fincludeAllLooseObjects) {
  // Loose Candidate information, sorted by pt
  fTree2->Branch("TwoProngLoose_pt",&fTwoProngLoose_pt);
  fTree2->Branch("TwoProngLoose_eta",&fTwoProngLoose_eta);
  fTree2->Branch("TwoProngLoose_phi",&fTwoProngLoose_phi);
  fTree2->Branch("TwoProngLoose_mass",&fTwoProngLoose_mass);
  fTree2->Branch("TwoProngLoose_mass_l",&fTwoProngLoose_mass_l);
  fTree2->Branch("TwoProngLoose_Mass",&fTwoProngLoose_Mass);
  fTree2->Branch("TwoProngLoose_Mass_l",&fTwoProngLoose_Mass_l);
  fTree2->Branch("TwoProngLoose_MassEta",&fTwoProngLoose_MassEta);
  fTree2->Branch("TwoProngLoose_MassEta_l",&fTwoProngLoose_MassEta_l);
  fTree2->Branch("TwoProngLoose_Mass300",&fTwoProngLoose_Mass300);
  fTree2->Branch("TwoProngLoose_FoundExtraTrack",&fTwoProngLoose_FoundExtraTrack);
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
  fTree2->Branch("TwoProngLoose_photon_pt_l",&fTwoProngLoose_photon_pt_l);
  fTree2->Branch("TwoProngLoose_photon_eta_l",&fTwoProngLoose_photon_eta_l);
  fTree2->Branch("TwoProngLoose_photon_phi_l",&fTwoProngLoose_photon_phi_l);
  fTree2->Branch("TwoProngLoose_photon_mass",&fTwoProngLoose_photon_mass);
  fTree2->Branch("TwoProngLoose_photon_mass_l",&fTwoProngLoose_photon_mass_l);
  fTree2->Branch("TwoProngLoose_photon_Mass",&fTwoProngLoose_photon_Mass);
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
  fTree2->Branch("TwoProngLoose_trackAsym",&fTwoProngLoose_trackAsym);
  fTree2->Branch("TwoProngLoose_photonAsym",&fTwoProngLoose_photonAsym);
  fTree2->Branch("TwoProngLoose_mPosPho_l",&fTwoProngLoose_mPosPho_l);
  fTree2->Branch("TwoProngLoose_mPosPho_pi0",&fTwoProngLoose_mPosPho_pi0);
  fTree2->Branch("TwoProngLoose_mPosPho_lpi0",&fTwoProngLoose_mPosPho_lpi0);
  fTree2->Branch("TwoProngLoose_mNegPho",&fTwoProngLoose_mNegPho);
  fTree2->Branch("TwoProngLoose_mNegPho_l",&fTwoProngLoose_mNegPho_l);
  fTree2->Branch("TwoProngLoose_mNegPho_pi0",&fTwoProngLoose_mNegPho_pi0);
  fTree2->Branch("TwoProngLoose_mNegPho_lpi0",&fTwoProngLoose_mNegPho_lpi0);
  fTree2->Branch("TwoProngLoose_mPosNeg",&fTwoProngLoose_mPosNeg);
  fTree2->Branch("TwoProngLoose_genOmega_dR",&fTwoProngLoose_genOmega_dR);
  fTree2->Branch("TwoProngLoose_genTau_dR",&fTwoProngLoose_genTau_dR);
  }
  // Tight Candidate information, sorted by pt
  fTree2->Branch("TwoProng_pt",&fTwoProng_pt);
  fTree2->Branch("TwoProng_eta",&fTwoProng_eta);
  fTree2->Branch("TwoProng_phi",&fTwoProng_phi);
  fTree2->Branch("TwoProng_mass",&fTwoProng_mass);
  fTree2->Branch("TwoProng_mass_l",&fTwoProng_mass_l);
  fTree2->Branch("TwoProng_Mass",&fTwoProng_Mass);
  fTree2->Branch("TwoProng_Mass_l",&fTwoProng_Mass_l);
  fTree2->Branch("TwoProng_MassEta",&fTwoProng_MassEta);
  fTree2->Branch("TwoProng_MassEta_l",&fTwoProng_MassEta_l);
  fTree2->Branch("TwoProng_Mass300",&fTwoProng_Mass300);
  fTree2->Branch("TwoProng_FoundExtraTrack",&fTwoProng_FoundExtraTrack);
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
  fTree2->Branch("TwoProng_photon_pt_l",&fTwoProng_photon_pt_l);
  fTree2->Branch("TwoProng_photon_eta_l",&fTwoProng_photon_eta_l);
  fTree2->Branch("TwoProng_photon_phi_l",&fTwoProng_photon_phi_l);
  fTree2->Branch("TwoProng_photon_mass_l",&fTwoProng_photon_mass_l);
  fTree2->Branch("TwoProng_photon_Mass",&fTwoProng_photon_Mass);
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
  fTree2->Branch("TwoProng_trackAsym",&fTwoProng_trackAsym);
  fTree2->Branch("TwoProng_photonAsym",&fTwoProng_photonAsym);
  fTree2->Branch("TwoProng_genOmega_dR",&fTwoProng_genOmega_dR);
  fTree2->Branch("TwoProng_genTau_dR",&fTwoProng_genTau_dR);
  fTree2->Branch("TwoProng_mPosPho",&fTwoProng_mPosPho);
  fTree2->Branch("TwoProng_mPosPho_l",&fTwoProng_mPosPho_l);
  fTree2->Branch("TwoProng_mPosPho_pi0",&fTwoProng_mPosPho_pi0);
  fTree2->Branch("TwoProng_mPosPho_lpi0",&fTwoProng_mPosPho_lpi0);
  fTree2->Branch("TwoProng_mNegPho",&fTwoProng_mNegPho);
  fTree2->Branch("TwoProng_mNegPho_l",&fTwoProng_mNegPho_l);
  fTree2->Branch("TwoProng_mNegPho_pi0",&fTwoProng_mNegPho_pi0);
  fTree2->Branch("TwoProng_mNegPho_lpi0",&fTwoProng_mNegPho_lpi0);
  fTree2->Branch("TwoProng_mPosNeg",&fTwoProng_mPosNeg);
  fTree2->Branch("TwoProng_CHpos_p3",&fTwoProng_CHpos_p3);
  fTree2->Branch("TwoProng_CHneg_p3",&fTwoProng_CHneg_p3);
  fTree2->Branch("TwoProng_photon_p3",&fTwoProng_photon_p3);
  // EG Photons, no high-pt id cuts, except electron veto
  fTree2->Branch("BasePhoton_pt",&fBaseIDPhoton_pt);
  fTree2->Branch("BasePhoton_eta",&fBaseIDPhoton_eta);
  fTree2->Branch("BasePhoton_phi",&fBaseIDPhoton_phi);
  fTree2->Branch("BasePhoton_mass",&fBaseIDPhoton_mass);
  fTree2->Branch("BasePhoton_iso_gamma",&fBaseIDPhoton_iso_gamma);
  fTree2->Branch("BasePhoton_iso_ch",&fBaseIDPhoton_iso_ch);
  fTree2->Branch("BasePhoton_HE",&fBaseIDPhoton_HE);
  fTree2->Branch("BasePhoton_sigmaieie",&fBaseIDPhoton_sigmaieie);
  // Loose Photons, every cut except iso gamma
  fTree2->Branch("LoosePhoton_pt",&fLooseIDPhoton_pt);
  fTree2->Branch("LoosePhoton_eta",&fLooseIDPhoton_eta);
  fTree2->Branch("LoosePhoton_phi",&fLooseIDPhoton_phi);
  fTree2->Branch("LoosePhoton_mass",&fLooseIDPhoton_mass);
  fTree2->Branch("LoosePhoton_iso_gamma",&fLooseIDPhoton_iso_gamma);
  // Tight Photons, sorted by pt
  fTree2->Branch("nTightPhotons",&fNumTightPhotons,"nTightPhotons/I");
  if (fAddDrConePhotonCut) {
  fTree2->Branch("nTightPhotons_ConeHE",&fNumTightPhotons_ConeHE,"nTightPhotons_ConeHE/I");
  }
  fTree2->Branch("Photon_pt",&fIDPhoton_pt);
  fTree2->Branch("Photon_eta",&fIDPhoton_eta);
  fTree2->Branch("Photon_phi",&fIDPhoton_phi);
  fTree2->Branch("Photon_mass",&fIDPhoton_mass);
  if(fincludeOldPhotons) {
  fTree2->Branch("Photon1",&fRecoTightPhotonInfo1,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  fTree2->Branch("Photon2",&fRecoTightPhotonInfo2,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  fTree2->Branch("Photon3",&fRecoTightPhotonInfo3,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  }
  if (fAddDrConePhotonCut) {
  fTree2->Branch("Photon_ConeHE_pt",&fID2Photon_pt);
  fTree2->Branch("Photon_ConeHE_eta",&fID2Photon_eta);
  fTree2->Branch("Photon_ConeHE_phi",&fID2Photon_phi);
  fTree2->Branch("Photon_ConeHE_mass",&fID2Photon_mass);
  }
  // Generator Objects
  fTree2->Branch("tauDecayType",&fTauDecayType,"tauDecayType/D");
  fTree2->Branch("GenTau_pt",&fGenTau_pt);
  fTree2->Branch("GenTau_eta",&fGenTau_eta);
  fTree2->Branch("GenTau_phi",&fGenTau_phi);
  fTree2->Branch("GenTau_mass",&fGenTau_mass);
  fTree2->Branch("GenTau_objDR",&fGenTau_objDR);
  fTree2->Branch("GenTau_candobjDR",&fGenTau_candobjDR);
  fTree2->Branch("GenPhi_pt",&fGenPhi_pt);
  fTree2->Branch("GenPhi_eta",&fGenPhi_eta);
  fTree2->Branch("GenPhi_phi",&fGenPhi_phi);
  fTree2->Branch("GenPhi_mass",&fGenPhi_mass);
  fTree2->Branch("GenPhi_px",&fGenPhi_px);
  fTree2->Branch("GenPhi_py",&fGenPhi_py);
  fTree2->Branch("GenPhi_pz",&fGenPhi_pz);
  fTree2->Branch("GenPhi_energy",&fGenPhi_energy);
  fTree2->Branch("GenOmega_pt",&fGenOmega_pt);
  fTree2->Branch("GenOmega_eta",&fGenOmega_eta);
  fTree2->Branch("GenOmega_phi",&fGenOmega_phi);
  fTree2->Branch("GenOmega_mass",&fGenOmega_mass);
  fTree2->Branch("GenOmega_px",&fGenOmega_px);
  fTree2->Branch("GenOmega_py",&fGenOmega_py);
  fTree2->Branch("GenOmega_pz",&fGenOmega_pz);
  fTree2->Branch("GenOmega_energy",&fGenOmega_energy); 
  fTree2->Branch("GenOmega_neutral_pt",&fGenOmega_neutral_pt); 
  fTree2->Branch("GenOmega_neutral_eta",&fGenOmega_neutral_eta); 
  fTree2->Branch("GenOmega_neutral_phi",&fGenOmega_neutral_phi); 
  fTree2->Branch("GenOmega_neutral_mass",&fGenOmega_neutral_mass); 
  fTree2->Branch("GenOmega_positive_pt",&fGenOmega_positive_pt); 
  fTree2->Branch("GenOmega_positive_eta",&fGenOmega_positive_eta); 
  fTree2->Branch("GenOmega_positive_phi",&fGenOmega_positive_phi); 
  fTree2->Branch("GenOmega_positive_mass",&fGenOmega_positive_mass); 
  fTree2->Branch("GenOmega_negative_pt",&fGenOmega_negative_pt); 
  fTree2->Branch("GenOmega_negative_eta",&fGenOmega_negative_eta); 
  fTree2->Branch("GenOmega_negative_phi",&fGenOmega_negative_phi); 
  fTree2->Branch("GenOmega_negative_mass",&fGenOmega_negative_mass); 
  fTree2->Branch("GenOmega_posnegdr",&fGenOmega_posnegdr); 
  fTree2->Branch("GenOmega_objDR",&fGenOmega_objDR); 
  fTree2->Branch("GenOmega_candobjDR",&fGenOmega_candobjDR); 
  fTree2->Branch("GenOmega_jetDR",&fGenOmega_jetDR); 
  // Combined Objects
  fTree2->Branch("RecoPhiDiTwoProng",&fRecoPhiDiTwoProng,ExoDiPhotons::recoDiObjectBranchDefString.c_str());
  fTree2->Branch("RecoPhiPhotonTwoProng",&fRecoPhiPhotonTwoProng,ExoDiPhotons::recoDiObjectBranchDefString.c_str());
  fTree2->Branch("RecoPhiInclusive",&fRecoPhiInclusive,ExoDiPhotons::recoDiObjectBranchDefString.c_str());
  }
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

  if (fFakeRateHistos) {
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
  }

  // photon trigger eff histos
  if (fTriggerEffHistos) {
  int eff_num_bins = 23;
  float eff_bins[eff_num_bins+1] = {0,10,12,14,16,18,20,22,25,30,40,60,80,100,120,140,160,180,200,300,400,500,600,2000};

  fPhotonTriggerEff_all_Numerator = 
  fs->make<TH1F>("photon_trig_eff_nume","HLT_Photon175 or HLT_Photon22_R9Id90_HE10_IsoM;leading high-pt-id-photon p_{T};Numerator",eff_num_bins,&eff_bins[0]);
  fPhotonTriggerEff_all_Denominator = 
  fs->make<TH1F>("photon_trig_eff_deno","HLT_Photon175 or HLT_Photon22_R9Id90_HE10_IsoM;leading high-pt-id-photon p_{T};Demominator", eff_num_bins,&eff_bins[0]);
  fPhotonTriggerEff_all_Division = 
  fs->make<TH1F>("photon_trig_eff","HLT_Photon175 or HLT_Photon22_R9Id90_HE10_IsoM;leading high-pt-id-photon p_{T};Efficency",eff_num_bins,&eff_bins[0]);

  fPhotonTriggerEff_Photon175_Numerator = 
  fs->make<TH1F>("photon_trig_eff_175_nume","HLT_Photon175;leading high-pt-id-photon p_{T};Numerator",eff_num_bins,&eff_bins[0]);
  fPhotonTriggerEff_Photon175_Division = 
  fs->make<TH1F>("photon_trig_eff_175","HLT_Photon175;leading high-pt-id-photon p_{T};Numerator",eff_num_bins,&eff_bins[0]);

  fPhotonTriggerEff_Photon22_Iso_Numerator = 
  fs->make<TH1F>("photon_trig_eff_22iso_nume","HLT_Photon22_R9Id90_HE10_IsoM;leading high-pt-id-photon p_{T};Numerator",eff_num_bins,&eff_bins[0]);
  fPhotonTriggerEff_Photon22_Iso_Division = 
  fs->make<TH1F>("photon_trig_eff_22iso","HLT_Photon22_R9Id90_HE10_IsoM;leading high-pt-id-photon p_{T};Numerator",eff_num_bins,&eff_bins[0]);

  if (fAddDrConePhotonCut) {
  fPhotonTriggerEff_ConeHE_all_Numerator = 
  fs->make<TH1F>("photon_trig_eff_nume_coneHE","HLT_Photon175 or HLT_Photon22_R9Id90_HE10_IsoM;leading high-pt-id+ConeHE5-photon p_{T};Numerator",eff_num_bins,&eff_bins[0]);
  fPhotonTriggerEff_ConeHE_all_Denominator = 
  fs->make<TH1F>("photon_trig_eff_deno_coneHE","HLT_Photon175 or HLT_Photon22_R9Id90_HE10_IsoM;leading high-pt-id+ConeHE5-photon p_{T};Demominator", eff_num_bins,&eff_bins[0]);
  fPhotonTriggerEff_ConeHE_all_Division = 
  fs->make<TH1F>("photon_trig_eff_coneHE","HLT_Photon175 or HLT_Photon22_R9Id90_HE10_IsoM;leading high-pt-id+ConeHE5-photon p_{T};Efficency",eff_num_bins,&eff_bins[0]);

  fPhotonTriggerEff_ConeHE_Photon175_Numerator = 
  fs->make<TH1F>("photon_trig_eff_175_nume_coneHE","HLT_Photon175;leading high-pt-id+ConeHE5-photon p_{T};Numerator",eff_num_bins,&eff_bins[0]);
  fPhotonTriggerEff_ConeHE_Photon175_Division = 
  fs->make<TH1F>("photon_trig_eff_175_coneHE","HLT_Photon175;leading high-pt-id+ConeHE5-photon p_{T};Numerator",eff_num_bins,&eff_bins[0]);

  fPhotonTriggerEff_ConeHE_Photon22_Iso_Numerator = 
  fs->make<TH1F>("photon_trig_eff_22iso_nume_coneHE","HLT_Photon22_R9Id90_HE10_IsoM;leading high-pt-id+ConeHE5-photon p_{T};Numerator",eff_num_bins,&eff_bins[0]);
  fPhotonTriggerEff_ConeHE_Photon22_Iso_Division = 
  fs->make<TH1F>("photon_trig_eff_22iso_coneHE","HLT_Photon22_R9Id90_HE10_IsoM;leading high-pt-id+ConeHE5-photon p_{T};Numerator",eff_num_bins,&eff_bins[0]);
  } // end DrCone conditional
  } // end trigger eff histos conditional

  if (fTwoProngYieldHistos) {
  fTwoProngYield = fs->make<TH1F>("twoprongyield_pt","TwoProng Object Yield", 50, 0, 5000);
  }

  if(fStackedDalitzHistos) {
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
}


ExoDiPhotonAnalyzer::~ExoDiPhotonAnalyzer()
{
 
}

// Helper methods
double 
ExoDiPhotonAnalyzer::iso_alpha(double eta, int ver)
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

double 
ExoDiPhotonAnalyzer::iso_area(double eta, int ver)
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

double 
ExoDiPhotonAnalyzer::iso_kappa(double eta, int ver)
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

  if (fDebug) cout << "event " << iEvent.id().event() << " lumi " <<  iEvent.id().luminosityBlock() << " run " << iEvent.id().run() << endl;

  // clear member vectors
  fCand_pt.clear();
  fCand_eta.clear();
  fCand_phi.clear();
  fCand_mass.clear();
  fCand_mass_l.clear();
  fCand_Mass.clear();
  fCand_Mass_l.clear();
  fCand_MassEta.clear();
  fCand_MassEta_l.clear();
  fCand_Mass300.clear();
  fCand_FoundExtraTrack.clear();
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
  fCand_photon_pt_l.clear();
  fCand_photon_eta_l.clear();
  fCand_photon_phi_l.clear();
  fCand_photon_mass_l.clear();
  fCand_photon_Mass.clear();
  fCand_photon_nGamma.clear();
  fCand_photon_nElectron.clear();
  fCand_chargedIso.clear();
  fCand_neutralIso.clear();
  fCand_egammaIso.clear();
  fCand_relchargedIso.clear();
  fCand_relneutralIso.clear();
  fCand_relegammaIso.clear();
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
  fCand_trackAsym.clear();
  fCand_photonAsym.clear();
  fCand_nChargedIsoCone.clear();
  fCand_nNeutralIsoCone.clear();
  fCand_nEGammaIsoCone.clear();
  fCand_genOmega_dR.clear();
  fCand_genTau_dR.clear();
  fCand_mPosPho.clear();
  fCand_mPosPho_l.clear();
  fCand_mPosPho_pi0.clear();
  fCand_mPosPho_lpi0.clear();
  fCand_mNegPho.clear();
  fCand_mNegPho_l.clear();
  fCand_mNegPho_pi0.clear();
  fCand_mNegPho_lpi0.clear();
  fCand_mPosNeg.clear();
  fCand_CHpos_p3.clear();
  fCand_CHneg_p3.clear();
  fCand_photon_p3.clear();
  fCand_tight.clear();
  fCand_passChargedIso.clear();
  fCand_passNeutralIso.clear();
  fCand_passEGammaIso.clear();
  fCand_passPhotonPt.clear();
  fCand_loose.clear();
  fCand_tight.clear();
  fCand_passChargedIso.clear();
  fCand_passNeutralIso.clear();
  fCand_passEGammaIso.clear();
  fCand_passPhotonPt.clear();
  fCand_loose.clear();
  fCand_iso_gamma.clear();
  fCand_iso_gamma_allPV.clear();
  fCand_iso_gamma_rel.clear();
  fCand_iso_ch.clear();
  fCand_iso_ch_allPV.clear();
  fCand_iso_ch_rel.clear();
  fCand_iso_gammacor1.clear();
  fCand_iso_gammacor2.clear();

  fTwoProngLoose_pt.clear();
  fTwoProngLoose_eta.clear();
  fTwoProngLoose_phi.clear();
  fTwoProngLoose_px.clear();
  fTwoProngLoose_py.clear();
  fTwoProngLoose_pz.clear();
  fTwoProngLoose_energy.clear();
  fTwoProngLoose_mass.clear();
  fTwoProngLoose_mass_l.clear();
  fTwoProngLoose_Mass.clear();
  fTwoProngLoose_Mass_l.clear();
  fTwoProngLoose_MassEta.clear();
  fTwoProngLoose_MassEta_l.clear();
  fTwoProngLoose_Mass300.clear();
  fTwoProngLoose_FoundExtraTrack.clear();
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
  fTwoProngLoose_photon_pt_l.clear();
  fTwoProngLoose_photon_eta_l.clear();
  fTwoProngLoose_photon_phi_l.clear();
  fTwoProngLoose_photon_mass_l.clear();
  fTwoProngLoose_photon_Mass.clear();
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
  fTwoProngLoose_trackAsym.clear();
  fTwoProngLoose_photonAsym.clear();
  fTwoProngLoose_mPosPho_l.clear();
  fTwoProngLoose_mPosPho_pi0.clear();
  fTwoProngLoose_mPosPho_lpi0.clear();
  fTwoProngLoose_mNegPho.clear();
  fTwoProngLoose_mNegPho_l.clear();
  fTwoProngLoose_mNegPho_pi0.clear();
  fTwoProngLoose_mNegPho_lpi0.clear();
  fTwoProngLoose_mPosNeg.clear();
  fTwoProngLoose_genOmega_dR.clear();
  fTwoProngLoose_genTau_dR.clear();

  fTwoProng_pt.clear();
  fTwoProng_eta.clear();
  fTwoProng_phi.clear();
  fTwoProng_px.clear();
  fTwoProng_py.clear();
  fTwoProng_pz.clear();
  fTwoProng_energy.clear();
  fTwoProng_mass.clear();
  fTwoProng_mass_l.clear();
  fTwoProng_Mass.clear();
  fTwoProng_Mass_l.clear();
  fTwoProng_MassEta.clear();
  fTwoProng_MassEta_l.clear();
  fTwoProng_Mass300.clear();
  fTwoProng_FoundExtraTrack.clear();
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
  fTwoProng_photon_pt_l.clear();
  fTwoProng_photon_eta_l.clear();
  fTwoProng_photon_phi_l.clear();
  fTwoProng_photon_mass_l.clear();
  fTwoProng_photon_Mass.clear();
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
  fTwoProng_trackAsym.clear();
  fTwoProng_photonAsym.clear();
  fTwoProng_genOmega_dR.clear();
  fTwoProng_genTau_dR.clear();
  fTwoProng_mPosPho.clear();
  fTwoProng_mPosPho_l.clear();
  fTwoProng_mPosPho_pi0.clear();
  fTwoProng_mPosPho_lpi0.clear();
  fTwoProng_mNegPho.clear();
  fTwoProng_mNegPho_l.clear();
  fTwoProng_mNegPho_pi0.clear();
  fTwoProng_mNegPho_lpi0.clear();
  fTwoProng_mPosNeg.clear();
  fTwoProng_CHpos_p3.clear();
  fTwoProng_CHneg_p3.clear();
  fTwoProng_photon_p3.clear();

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

  fIDPhoton_pt.clear();
  fIDPhoton_eta.clear();
  fIDPhoton_phi.clear();
  fIDPhoton_mass.clear();

  fBaseIDPhoton_pt.clear();
  fBaseIDPhoton_eta.clear();
  fBaseIDPhoton_phi.clear();
  fBaseIDPhoton_mass.clear();
  fBaseIDPhoton_iso_gamma.clear();
  fBaseIDPhoton_iso_ch.clear();
  fBaseIDPhoton_HE.clear();
  fBaseIDPhoton_sigmaieie.clear();
  fLooseIDPhoton_pt.clear();
  fLooseIDPhoton_eta.clear();
  fLooseIDPhoton_phi.clear();
  fLooseIDPhoton_mass.clear();
  fLooseIDPhoton_iso_gamma.clear();
  fID2Photon_pt.clear();
  fID2Photon_eta.clear();
  fID2Photon_phi.clear();
  fID2Photon_mass.clear();

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
  fGenPhi_px.clear();
  fGenPhi_py.clear();
  fGenPhi_pz.clear();
  fGenPhi_energy.clear();
  fGenOmega_pt.clear();
  fGenOmega_eta.clear();
  fGenOmega_phi.clear();
  fGenOmega_mass.clear();
  fGenOmega_px.clear();
  fGenOmega_py.clear();
  fGenOmega_pz.clear();
  fGenOmega_energy.clear();
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

  // trigger 
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects);
  iEvent.getByToken(triggerPrescales_, triggerPrescales);

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  bool found_175 = false;
  bool found_22_iso = false;
  string trigger_photon175  = "HLT_Photon175_v";
  string trigger_photon22_iso  = "HLT_Photon22_R9Id90_HE10_IsoM_v";
  bool bit_photon175 = false;
  bool bit_photon22_iso = false;
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i)
  {
     string triggerName = names.triggerName(i);
     
     std::size_t pos_175 = triggerName.find(trigger_photon175);
     if ( pos_175 != std::string::npos ) {
       found_175 = true;
       bit_photon175 = triggerBits->accept(i);
     }

     std::size_t pos_22_iso = triggerName.find(trigger_photon22_iso);
     if ( pos_22_iso != std::string::npos ) {
       found_22_iso = true;
       bit_photon22_iso = triggerBits->accept(i);
     }
  }
  if(fDebug && !found_175) cout << "didn't find trigger! : " << trigger_photon175 << endl;
  if(fDebug && !found_22_iso) cout << "didn't find trigger! : " << trigger_photon22_iso << endl;
  if(fDebug && found_175) cout << "found : " << trigger_photon175 << ", bit: " << bit_photon175 << endl;
  if(fDebug && found_22_iso) cout << "found : " << trigger_photon22_iso << ", bit: " << bit_photon22_iso << endl;
  fHLT_Photon175 = bit_photon175;
  fHLT_Photon22_Iso = bit_photon22_iso;

  // ecal tool
  lazyTools_ = std::auto_ptr<noZS::EcalClusterLazyTools>( new noZS::EcalClusterLazyTools(iEvent, iSetup, recHitsEBToken, recHitsEEToken));   

  // get event products
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

  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);

  edm::ESHandle<EcalChannelStatus> chStatus;
  iSetup.get<EcalChannelStatusRcd>().get(chStatus);
  const EcalChannelStatus *ch_status = chStatus.product(); 

  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(pfcandsToken_, pfcands);

  edm::Handle<vector<reco::GenParticle>> genparticles;
  if (fincludeSignalGenParticles || frunningOnTauTauMC) {
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

  edm::Handle<edm::View<pat::Photon> > ged_photons;
  iEvent.getByToken(gedphotonsToken_,ged_photons);

  edm::Handle<std::vector<pat::MET>> MET;
  iEvent.getByToken(metToken_, MET);

  edm::Handle<std::vector<pat::Electron>> electrons;
  iEvent.getByToken(electronToken_, electrons);

  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonToken_, muons);

  edm::Handle<GenEventInfoProduct> genEventInfo;
  if (fincludeMCInfo) {
    iEvent.getByToken(genEventInfoToken_, genEventInfo);
    // fill MC weights
    fMcW = genEventInfo->weight();
    fMcWProd = genEventInfo->weightProduct();
  }

  // Signal MC Generator information
  if (fincludeSignalGenParticles) {
    for (unsigned int i = 0; i < genparticles->size(); i++) {
      const reco::GenParticle &genparticle = (*genparticles)[i];
      if (genparticle.pdgId() != 9000006 || genparticle.status() != 62) continue;
      TLorentzVector resonance;
      resonance.SetPtEtaPhiM(genparticle.pt(),genparticle.eta(),genparticle.phi(),genparticle.mass());
      fGenPhi_pt.push_back(resonance.Pt());
      fGenPhi_eta.push_back(resonance.Eta());
      fGenPhi_phi.push_back(resonance.Phi());
      fGenPhi_mass.push_back(resonance.M());
      fGenPhi_px.push_back(resonance.Px());
      fGenPhi_py.push_back(resonance.Py());
      fGenPhi_pz.push_back(resonance.Pz());
      fGenPhi_energy.push_back(resonance.E());
      if (fDebug) cout << genparticle.pdgId() << ", status=" << genparticle.status() << endl; 
      for (unsigned int j = 0; j < genparticle.numberOfDaughters(); j++) {
        const reco::Candidate* genparticle2 = genparticle.daughter(j);
        if (fDebug) cout << "-> " << genparticle2->pdgId() << ", status=" << genparticle2->status() << endl; 
        TLorentzVector pseudoscalar;
        pseudoscalar.SetPtEtaPhiE(genparticle2->pt(),genparticle2->eta(),genparticle2->phi(),genparticle2->energy());
        TLorentzVector positivePion;
        TLorentzVector negativePion;
        TLorentzVector neutralContent;
        neutralContent.SetXYZT(0,0,0,0);
        for (unsigned int jj = 0; jj < genparticle2->numberOfDaughters(); jj++) {
          const reco::Candidate* genparticle3 = genparticle2->daughter(jj);
          if (fDebug) cout << "  -> " << genparticle3->pdgId() << ", status=" << genparticle3->status() << endl; 
          TLorentzVector genparticle3Vect;
          genparticle3Vect.SetPtEtaPhiE(genparticle3->pt(), genparticle3->eta(), genparticle3->phi(), genparticle3->energy());
          if(genparticle3->pdgId()==111 || genparticle3->pdgId()==22) neutralContent += genparticle3Vect;
          if(fDebug && (genparticle3->pdgId()==111 || genparticle3->pdgId()==22)) cout << "  added neutral content: "<< genparticle3->pdgId() << endl;
          if(genparticle3->pdgId()==211) positivePion.SetPtEtaPhiE(genparticle3->pt(),genparticle3->eta(),genparticle3->phi(),genparticle3->energy());
          if(genparticle3->pdgId()==-211) negativePion.SetPtEtaPhiE(genparticle3->pt(),genparticle3->eta(),genparticle3->phi(),genparticle3->energy());
          for (unsigned int jjj = 0; jjj < genparticle3->numberOfDaughters(); jjj++) {
            const reco::Candidate* genparticle4 = genparticle3->daughter(jjj);
            if (fDebug) cout << "    -> " << genparticle4->pdgId() << ", status=" << genparticle4->status() << endl;
            if(genparticle3->pdgId()==111 || genparticle3->pdgId()==22) continue;
            TLorentzVector genparticle4Vect;
            genparticle4Vect.SetPtEtaPhiE(genparticle4->pt(), genparticle4->eta(), genparticle4->phi(), genparticle4->energy());
            if(genparticle4->pdgId()==111 || genparticle4->pdgId()==22) neutralContent += genparticle4Vect;
            if(fDebug && (genparticle4->pdgId()==111 || genparticle4->pdgId()==22)) cout << "    added neutral content: "<< genparticle4->pdgId() << endl;
            if(genparticle4->pdgId()==211) positivePion.SetPtEtaPhiE(genparticle4->pt(),genparticle4->eta(),genparticle4->phi(),genparticle4->energy());
            if(genparticle4->pdgId()==-211) negativePion.SetPtEtaPhiE(genparticle4->pt(),genparticle4->eta(),genparticle4->phi(),genparticle4->energy());
          }
        }
        fGenOmega_pt.push_back(pseudoscalar.Pt());
        fGenOmega_eta.push_back(pseudoscalar.Eta());
        fGenOmega_phi.push_back(pseudoscalar.Phi());
        fGenOmega_mass.push_back(pseudoscalar.M());
        fGenOmega_px.push_back(pseudoscalar.Px());
        fGenOmega_py.push_back(pseudoscalar.Py());
        fGenOmega_pz.push_back(pseudoscalar.Pz());
        fGenOmega_energy.push_back(pseudoscalar.E());
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
  if (fDebug) cout << ". done generator logic" << endl;

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
  fNumAK4jets = ak4jets->size();

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
  fNumPhotons = miniaod_photons->size();

  // Event wide information
  fEventNum = iEvent.id().event();
  fRunNum = iEvent.id().run();
  fLumiNum = iEvent.id().luminosityBlock();
  fNumElectrons = electrons->size();
  fNumMuons = muons->size();
  fHT = 0.0;
  for (unsigned int i = 0; i < ak4jets->size(); i++) {
    const pat::Jet &jet = (*ak4jets)[i];
    if (jet.pt() < 30) continue;
    if (fabs(jet.eta()) > 2.5) continue;
    fHT += jet.pt();
  }
  fMET = (*MET)[0].pt();
  fMET_phi = (*MET)[0].phi();
  fNumPVs = primaryvertecies->size();
  fRho = *rhoH;
  fNumPF = pfcands->size();

  // Two prongs
  if (fDebug) cout << ". starting two prong code" << endl;
  // Find all pairs of one CH pos and one CH neg within specified DR of each other
  int nLoose = 0;
  int pruned_count = 0;
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
        // search for extra tracks
        vector<unsigned int> extras;
        for (unsigned int e = 0; e < pfcands->size(); e++) {
          const pat::PackedCandidate &pfextra = (*pfcands)[e];
          if (pfextra.pt() < fCandidatePairMinPt) continue;
          if (pfextra.fromPV()<=1) continue;
          if (abs(pfextra.pdgId()) != 211) continue;
          TLorentzVector pfcandextra;
          pfcandextra.SetPtEtaPhiE(pfextra.pt(), pfextra.eta(), pfextra.phiAtVtx(), pfextra.energy());
          double dr = max( pfcand1.DeltaR(pfcandextra), pfcand2.DeltaR(pfcandextra) );
          if (dr > fCandidatePairDR) continue;
          extras.push_back(e);
        }
        unsigned int extra = 99999;
        bool found_extra_track = false;
        TLorentzVector pfcandextra;       
        if (extras.size() == 1) {
          const pat::PackedCandidate &pfextra = (*pfcands)[extras[0]];
          extra = extras[0];
          pfcandextra.SetPtEtaPhiE(pfextra.pt(), pfextra.eta(), pfextra.phiAtVtx(), pfextra.energy());
          found_extra_track = true;
        }
        // now define photon
        TLorentzVector center;
        if (fCandidateOptionalExtraTrack && found_extra_track) {
          center = pfcand1 + pfcand2 + pfcandextra; }
        else {
          center = pfcand1 + pfcand2;}
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
          TLorentzVector TwoProngObject;
          TwoProngObject = center + photon;

          leading_pf_photon.SetPtEtaPhiE((*pfcands)[n].pt(), (*pfcands)[n].eta(), (*pfcands)[n].phiAtVtx(), (*pfcands)[n].energy());

          TLorentzVector leading_pf_photon_as_pi0;
          leading_pf_photon_as_pi0.SetPtEtaPhiM((*pfcands)[n].pt(), (*pfcands)[n].eta(), (*pfcands)[n].phiAtVtx(), PI0_MASS);

          TLorentzVector leading_pf_photon_as_eta;
          leading_pf_photon_as_eta.SetPtEtaPhiM((*pfcands)[n].pt(), (*pfcands)[n].eta(), (*pfcands)[n].phiAtVtx(), ETA_MASS);

          TLorentzVector photon_summed_as_pi0;
          photon_summed_as_pi0.SetPtEtaPhiM(photon.Pt(), photon.Eta(), photon.Phi(), PI0_MASS);
           
          TLorentzVector photon_summed_as_eta;
          photon_summed_as_eta.SetPtEtaPhiM(photon.Pt(), photon.Eta(), photon.Phi(), ETA_MASS);

          TLorentzVector photon_summed_as_300;
          photon_summed_as_300.SetPtEtaPhiM(photon.Pt(), photon.Eta(), photon.Phi(), 0.3);

          if (fabs(TwoProngObject.Eta()) > fCandidateAbsMaxEta) continue;
          if (TwoProngObject.Pt() < fCandidateMinPt) continue;
          double track_asymmetry = min(pfcand1.Pt(),pfcand2.Pt()) / max(pfcand1.Pt(),pfcand2.Pt());
          double photon_asymmetry = min(pfcand1.Pt()+pfcand2.Pt(),photon.Pt()) / max(pfcand1.Pt()+pfcand2.Pt(),photon.Pt());
          if(track_asymmetry < fCandidateTrackAsymmetryCut) continue;
          if(photon_asymmetry < fCandidatePhotonAsymmetryCut) continue;
          // CH pair within dr range
          // has at least one pf photon
          // |eta| <= eta requirement
          // pt >= min pt requirement
          // passes asymmetry cuts
          //   meets definition of candidate twoprong, fill candidate vectors
          TLorentzVector poscand;
          TLorentzVector negcand;
          if (pf1.pdgId() > 0) {
            poscand = pfcand1;
            negcand = pfcand2;
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
            poscand = pfcand2;
            negcand = pfcand1;
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
          fCand_photon_pt_l.push_back(leading_pf_photon.Pt());
          fCand_photon_eta_l.push_back(leading_pf_photon.Eta());
          fCand_photon_phi_l.push_back(leading_pf_photon.Phi());
          fCand_photon_mass_l.push_back(leading_pf_photon.M());
          fCand_photon_Mass.push_back(PI0_MASS);
          fCand_photon_nGamma.push_back(numgamma);
          fCand_photon_nElectron.push_back(nume);
          fCand_pt.push_back(TwoProngObject.Pt());
          fCand_eta.push_back(TwoProngObject.Eta());
          fCand_phi.push_back(TwoProngObject.Phi());

          fCand_mass.push_back(TwoProngObject.M());
          fCand_Mass.push_back( (center + photon_summed_as_pi0).M() );
          fCand_mass_l.push_back( (center + leading_pf_photon).M() );
          fCand_Mass_l.push_back( (center + leading_pf_photon_as_pi0).M() );
          fCand_MassEta.push_back( (center + photon_summed_as_eta).M() );
          fCand_MassEta_l.push_back( (center + leading_pf_photon_as_eta).M() );
          fCand_Mass300.push_back( (center + photon_summed_as_300).M() );

          fCand_trackAsym.push_back(track_asymmetry);
          fCand_photonAsym.push_back(photon_asymmetry);

          fCand_mPosPho.push_back( (poscand + photon).M2() );
          fCand_mPosPho_l.push_back( (poscand + leading_pf_photon).M2() );
          fCand_mPosPho_pi0.push_back( (poscand + photon_summed_as_pi0).M2() );
          fCand_mPosPho_lpi0.push_back( (poscand + leading_pf_photon_as_pi0).M2() );
          fCand_mNegPho.push_back( (negcand + photon).M2() );
          fCand_mNegPho_l.push_back( (negcand + leading_pf_photon).M2() );
          fCand_mNegPho_pi0.push_back( (negcand + photon_summed_as_pi0).M2() );
          fCand_mNegPho_lpi0.push_back( (negcand + leading_pf_photon_as_pi0).M2() );
          fCand_mPosNeg.push_back( (poscand + negcand).M2() );

          fCand_CHpos_p3.push_back(poscand.P());
          fCand_CHneg_p3.push_back(negcand.P());
          fCand_photon_p3.push_back(photon.P());

          fCand_FoundExtraTrack.push_back(found_extra_track);

          // Now define isolations
          double chargedIso = 0;
          double chargedIso_allPV = 0;
          double neutralIso = 0;
          double egammaIso = 0;
          double egammaIso_allPV = 0;
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
                if ((found_extra_track && m != extra) || !found_extra_track) { // if we are including an extra track, skip it
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
          for (unsigned int m = 0; m < pfcands->size(); m++) {
            const pat::PackedCandidate &pf4 = (*pfcands)[m];
            TLorentzVector pfcand4;
            pfcand4.SetPtEtaPhiE(pf4.pt(), pf4.eta(), pf4.phiAtVtx(), pf4.energy());
            // charged (incl. muons)
            if (abs(pf4.pdgId()) == 13 || abs(pf4.pdgId()) == 211) {
              if ( center.DeltaR(pfcand4) < fCandidatePairIsolationDR && !(m == i || m == j) ) { // don't include one of CH from CH pair
                  chargedIso_allPV += pfcand4.Pt();
               }
            // e gamma
            } else if (abs(pf4.pdgId()) == 11 || pf4.pdgId() == 22) {
              if ( (center.DeltaR(pfcand4) < fCandidatePairIsolationDR) &&
                   !(fabs(pf4.phiAtVtx() - center.Phi()) < fCandidatePairPhiBox/2.0 && fabs(pf4.eta() - center.Eta()) < fCandidatePairEtaBox/2.0)) {
                egammaIso_allPV += pfcand4.Pt();
              }
            }
          } // end pf cand loop
          double relchargedIso = chargedIso / TwoProngObject.Pt();
          double relneutralIso = neutralIso / TwoProngObject.Pt();
          double relegammaIso = egammaIso / TwoProngObject.Pt();
          if (fDebug) cout << ". finished isolation" << endl;
          fCand_chargedIso.push_back(chargedIso);
          fCand_neutralIso.push_back(neutralIso);
          fCand_egammaIso.push_back(egammaIso);
          fCand_relchargedIso.push_back(relchargedIso);
          fCand_relneutralIso.push_back(relneutralIso);
          fCand_relegammaIso.push_back(relegammaIso);
          fCand_nChargedIsoCone.push_back(chargedIsoCount);
          fCand_nNeutralIsoCone.push_back(neutralIsoCount);
          fCand_nEGammaIsoCone.push_back(egammaIsoCount);
          
          double etaForIso = TwoProngObject.Eta();
          double isoGammaCor1 = egammaIso + iso_alpha(etaForIso, 1) - fRho * iso_area(etaForIso, 1) - iso_kappa(etaForIso, 1) * TwoProngObject.Pt();
          double isoGammaCor2 = egammaIso + iso_alpha(etaForIso, 2) - fRho * iso_area(etaForIso, 2) - iso_kappa(etaForIso, 2) * TwoProngObject.Pt();

          fCand_iso_gamma.push_back(egammaIso);
          fCand_iso_gamma_allPV.push_back(egammaIso_allPV);
          fCand_iso_gamma_rel.push_back(egammaIso/TwoProngObject.Pt());
          fCand_iso_ch.push_back(chargedIso);
          fCand_iso_ch_allPV.push_back(chargedIso_allPV);
          fCand_iso_ch_rel.push_back(chargedIso/TwoProngObject.Pt());
          fCand_iso_gammacor1.push_back(isoGammaCor1);
          fCand_iso_gammacor2.push_back(isoGammaCor2);
          
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
          fCand_genOmega_dR.push_back(genOmega_dR);
          if (fDebug) cout << ". finished gen omega matching" << endl;
          // Generator matching to tau
          double genTau_dR = 99.9;
	        if (frunningOnTauTauMC) {
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
          fCand_genTau_dR.push_back(genTau_dR);
          if (fDebug) cout << ". finished gen tau matching" << endl;

          if (loose) nLoose++;
          if (loose) { if (TwoProngObject.Pt() > LeadingFakeTwoProng.Pt()) LeadingFakeTwoProng = TwoProngObject; }
          // Fake rate Analysis histograms
          if (fFakeRateHistos) {
          if (tight) {
            fTwoProngFakeNume_even_pt->Fill(TwoProngObject.Pt(), TwoProngObject.M());
            fTwoProngFakeNume_even_eta->Fill(TwoProngObject.Eta(), TwoProngObject.M());
            fTwoProngFakeNume_even_phi->Fill(TwoProngObject.Phi(), TwoProngObject.M());
            fTwoProngFakeNume_odd_pt->Fill(TwoProngObject.Pt(), TwoProngObject.M());
            fTwoProngFakeNume_odd_eta->Fill(TwoProngObject.Eta(), TwoProngObject.M());
            fTwoProngFakeNume_odd_phi->Fill(TwoProngObject.Phi(), TwoProngObject.M());
          } if (loose) {
            fTwoProngFakeDeno_even_pt->Fill(TwoProngObject.Pt(), TwoProngObject.M());
            fTwoProngFakeDeno_even_eta->Fill(TwoProngObject.Eta(), TwoProngObject.M());
            fTwoProngFakeDeno_even_phi->Fill(TwoProngObject.Phi(), TwoProngObject.M());
            fTwoProngFakeDeno_odd_pt->Fill(TwoProngObject.Pt(), TwoProngObject.M());
            fTwoProngFakeDeno_odd_eta->Fill(TwoProngObject.Eta(), TwoProngObject.M());
            fTwoProngFakeDeno_odd_phi->Fill(TwoProngObject.Phi(), TwoProngObject.M());
          }
          if (fDebug) cout << ". finished fake rate filling" << endl;
          }
        }
      } // end conditionals on CH pair
    }
  } // end making candidates
  fNumPrunedPF = pruned_count;
  fNumTwoProngLoose = nLoose;

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
      fTwoProngLoose_mass_l.push_back(fCand_mass_l[index]);
      fTwoProngLoose_Mass.push_back(fCand_Mass[index]);
      fTwoProngLoose_Mass_l.push_back(fCand_Mass_l[index]);
      fTwoProngLoose_MassEta.push_back(fCand_MassEta[index]);
      fTwoProngLoose_MassEta_l.push_back(fCand_MassEta_l[index]);
      fTwoProngLoose_Mass300.push_back(fCand_Mass300[index]);
      fTwoProngLoose_FoundExtraTrack.push_back(fCand_FoundExtraTrack[index]);
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
      fTwoProngLoose_photon_pt_l.push_back(fCand_photon_pt_l[index]);
      fTwoProngLoose_photon_eta_l.push_back(fCand_photon_eta_l[index]);
      fTwoProngLoose_photon_phi_l.push_back(fCand_photon_phi_l[index]);
      fTwoProngLoose_photon_mass_l.push_back(fCand_photon_mass_l[index]);
      fTwoProngLoose_photon_Mass.push_back(fCand_photon_Mass[index]);
      fTwoProngLoose_photon_nGamma.push_back(fCand_photon_nGamma[index]);
      fTwoProngLoose_photon_nElectron.push_back(fCand_photon_nElectron[index]);
      fTwoProngLoose_chargedIso.push_back(fCand_chargedIso[index]);
      fTwoProngLoose_neutralIso.push_back(fCand_neutralIso[index]);
      fTwoProngLoose_egammaIso.push_back(fCand_egammaIso[index]);
      fTwoProngLoose_genOmega_dR.push_back(fCand_genOmega_dR[index]);
      fTwoProngLoose_genTau_dR.push_back(fCand_genTau_dR[index]);
      fTwoProngLoose_trackAsym.push_back(fCand_trackAsym[index]);
      fTwoProngLoose_photonAsym.push_back(fCand_photonAsym[index]);
      fTwoProngLoose_mPosPho_l.push_back(fCand_mPosPho_l[index]);
      fTwoProngLoose_mPosPho_pi0.push_back(fCand_mPosPho_pi0[index]);
      fTwoProngLoose_mPosPho_lpi0.push_back(fCand_mPosPho_lpi0[index]);
      fTwoProngLoose_mNegPho.push_back(fCand_mNegPho[index]);
      fTwoProngLoose_mNegPho_l.push_back(fCand_mNegPho_l[index]);
      fTwoProngLoose_mNegPho_pi0.push_back(fCand_mNegPho_pi0[index]);
      fTwoProngLoose_mNegPho_lpi0.push_back(fCand_mNegPho_lpi0[index]);
      fTwoProngLoose_mPosNeg.push_back(fCand_mPosNeg[index]);
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
      fTwoProng_mass_l.push_back(fCand_mass_l[index]);
      fTwoProng_Mass.push_back(fCand_Mass[index]);
      fTwoProng_Mass_l.push_back(fCand_Mass_l[index]);
      fTwoProng_MassEta.push_back(fCand_MassEta[index]);
      fTwoProng_MassEta_l.push_back(fCand_MassEta_l[index]);
      fTwoProng_Mass300.push_back(fCand_Mass300[index]);
      fTwoProng_FoundExtraTrack.push_back(fCand_FoundExtraTrack[index]);
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
      fTwoProng_photon_pt_l.push_back(fCand_photon_pt_l[index]);
      fTwoProng_photon_eta_l.push_back(fCand_photon_eta_l[index]);
      fTwoProng_photon_phi_l.push_back(fCand_photon_phi_l[index]);
      fTwoProng_photon_mass_l.push_back(fCand_photon_mass_l[index]);
      fTwoProng_photon_Mass.push_back(fCand_photon_Mass[index]);
      fTwoProng_photon_nGamma.push_back(fCand_photon_nGamma[index]);
      fTwoProng_photon_nElectron.push_back(fCand_photon_nElectron[index]);
      fTwoProng_chargedIso.push_back(fCand_chargedIso[index]);
      fTwoProng_neutralIso.push_back(fCand_neutralIso[index]);
      fTwoProng_egammaIso.push_back(fCand_egammaIso[index]);
      fTwoProng_genOmega_dR.push_back(fCand_genOmega_dR[index]);
      fTwoProng_genTau_dR.push_back(fCand_genTau_dR[index]);
      fTwoProng_trackAsym.push_back(fCand_trackAsym[index]);
      fTwoProng_photonAsym.push_back(fCand_photonAsym[index]);
      fTwoProng_mPosPho.push_back(fCand_mPosPho[index]);
      fTwoProng_mPosPho_l.push_back(fCand_mPosPho_l[index]);
      fTwoProng_mPosPho_pi0.push_back(fCand_mPosPho_pi0[index]);
      fTwoProng_mPosPho_lpi0.push_back(fCand_mPosPho_lpi0[index]);
      fTwoProng_mNegPho.push_back(fCand_mNegPho[index]);
      fTwoProng_mNegPho_l.push_back(fCand_mNegPho_l[index]);
      fTwoProng_mNegPho_pi0.push_back(fCand_mNegPho_pi0[index]);
      fTwoProng_mNegPho_lpi0.push_back(fCand_mNegPho_lpi0[index]);
      fTwoProng_mPosNeg.push_back(fCand_mPosNeg[index]);
      fTwoProng_CHpos_p3.push_back(fCand_CHpos_p3[index]);
      fTwoProng_CHneg_p3.push_back(fCand_CHneg_p3[index]);
      fTwoProng_photon_p3.push_back(fCand_photon_p3[index]);
    }
  }
  fNumTwoProng = fCand_pt.size();
  fNumTwoProngPass = fTwoProng_pt.size();
  if (fDebug) cout << ". finished passed collections" << endl;

  // tau decay type
  if (frunningOnTauTauMC) {
    vector<string> leptons;
    for (unsigned int i = 0; i < genparticles->size(); i++) {
      const reco::GenParticle & genparticle = (*genparticles)[i];
      if (genparticle.status() != 22) continue;
      leptons = getDecay(genparticle);
    }
    if (leptons.size() == 0) fTauDecayType = 0;
    if (leptons.size() == 1) fTauDecayType = -1;
    if (leptons.size() == 2) 
    {
      if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end()) fTauDecayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end()) fTauDecayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end()) fTauDecayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end()) fTauDecayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end()) fTauDecayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end()) fTauDecayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end()) fTauDecayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end()) fTauDecayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end()) fTauDecayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end()) fTauDecayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end()) fTauDecayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end()) fTauDecayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end()) fTauDecayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end()) fTauDecayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end()) fTauDecayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end()) fTauDecayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end()) fTauDecayType = 4;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end()) fTauDecayType = 4;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end()) fTauDecayType = 4;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end()) fTauDecayType = 4;
      else if (std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end()) fTauDecayType = 4;
      else if (std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end()) fTauDecayType = 4;
      else if (std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end()) fTauDecayType = 4;
      else if (std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end()) fTauDecayType = 4;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) fTauDecayType = 5.2;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) fTauDecayType = 5.1;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) fTauDecayType = 5.4;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) fTauDecayType = 5.3;
      else if (std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end()) fTauDecayType = 5.2;
      else if (std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end()) fTauDecayType = 5.1;
      else if (std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end()) fTauDecayType = 5.4;
      else if (std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end()) fTauDecayType = 5.3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end()) fTauDecayType = 6;
      else if (std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) fTauDecayType = 7;
      else if (std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) fTauDecayType = 8;
      else if (std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end()) fTauDecayType = 8;
      else if (std::find(leptons.begin(), leptons.end(), "e+") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "e-") != leptons.end()) fTauDecayType = 1;
      else if (std::find(leptons.begin(), leptons.end(), "mu+") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "mu-") != leptons.end()) fTauDecayType = 2;
      else
        fTauDecayType = 9;
    }
    if (leptons.size() > 2) fTauDecayType = 10;
  }

  // More matching, by gen tau perspective now
  if (frunningOnTauTauMC) {
    for (unsigned int i = 0; i < genparticles->size(); i++) {
      const reco::GenParticle & genparticle = (*genparticles)[i];
      if (abs(genparticle.pdgId()) != 15) continue;
      TLorentzVector GenParticle;
      GenParticle.SetPtEtaPhiM(genparticle.pt(), genparticle.eta(), genparticle.phi(), genparticle.mass());
      fGenTau_pt.push_back(genparticle.pt());
      fGenTau_eta.push_back(genparticle.eta());
      fGenTau_phi.push_back(genparticle.phi());
      fGenTau_mass.push_back(genparticle.mass());
      double candDR = 99.9;
      double passedCandDR = 99.9;
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
      fGenTau_objDR.push_back(passedCandDR);
      fGenTau_candobjDR.push_back(candDR);
    }
  }

  // More matching, by gen omega perspective now
  if (fincludeSignalGenParticles) {
    for (unsigned int i = 0; i < genparticles->size(); i++) {
      const reco::GenParticle &genparticle = (*genparticles)[i];
      if (genparticle.pdgId() != 9000006 || genparticle.status() != 62) continue;
      for (unsigned int j = 0; j < genparticle.numberOfDaughters(); j++) {
        const reco::Candidate* genparticle2 = genparticle.daughter(j);
        TLorentzVector GenParticle;
        GenParticle.SetPtEtaPhiM(genparticle2->pt(), genparticle2->eta(), genparticle2->phi(), genparticle2->mass());
        if (genparticle2->pdgId() == 221 || genparticle2->pdgId() == 331) {
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
          fGenOmega_objDR.push_back(passedCandDR);
          fGenOmega_candobjDR.push_back(candDR);
          fGenOmega_jetDR.push_back(jetDR);
        }
      }
    }
  }

  if (fDebug) cout << ". high-pt-id" << endl;

  // High-pT Photon id  
  const CaloSubdetectorTopology* subDetTopologyEB_;
  const CaloSubdetectorTopology* subDetTopologyEE_;
  edm::ESHandle<CaloTopology> caloTopology;
  iSetup.get<CaloTopologyRecord>().get(caloTopology);
  subDetTopologyEB_ = caloTopology->getSubdetectorTopology(DetId::Ecal,EcalBarrel);
  subDetTopologyEE_ = caloTopology->getSubdetectorTopology(DetId::Ecal,EcalEndcap);
  bool isSat = false;
  std::vector<edm::Ptr<pat::Photon>> basePhotons;
  std::vector<edm::Ptr<pat::Photon>> loosePhotons;
  std::vector<edm::Ptr<pat::Photon>> goodPhotons;
  std::vector<edm::Ptr<pat::Photon>> goodPhotons_ConeHE;
  for (size_t i = 0; i < ged_photons->size(); ++i) {
    const auto pho = ged_photons->ptrAt(i);
    isSat = photon_isSaturated(&(*pho), &(*recHitsEB), &(*recHitsEE), &(*subDetTopologyEB_), &(*subDetTopologyEE_));
    bool passID = photon_passHighPtID(&(*pho), fRho, isSat);
    if(passID) {
      goodPhotons.push_back(pho);
    }
    bool passIDConeHE = photon_passHighPtID_DrCone(&(*pho), fRho, isSat);
    if(passIDConeHE) {
      goodPhotons_ConeHE.push_back(pho);
    }
    bool passLooseID = photon_passHighPtID_loose(&(*pho), fRho, isSat);
    if(passLooseID) {
      loosePhotons.push_back(pho);
    }
    bool passBaseID = photon_passHighPtID_base(&(*pho), fRho, isSat);
    if(passBaseID) {
      basePhotons.push_back(pho);
    }
  }
  fNumTightPhotons = goodPhotons.size();
  if (fDebug) cout << ". done making photon collections" << endl;
  // sort and fill
  sort(loosePhotons.begin(),loosePhotons.end(),compareCandsByPt);
  sort(goodPhotons.begin(),goodPhotons.end(),compareCandsByPt);
  sort(goodPhotons_ConeHE.begin(),goodPhotons_ConeHE.end(),compareCandsByPt);
  if (fDebug) cout << ". done sorting photon collections" << endl;
  for (unsigned int i = 0; i < basePhotons.size(); i++ )
  {
    fBaseIDPhoton_pt.push_back( (*basePhotons[i]).pt() );
    fBaseIDPhoton_eta.push_back( (*basePhotons[i]).eta() );
    fBaseIDPhoton_phi.push_back( (*basePhotons[i]).phi() );
    fBaseIDPhoton_mass.push_back( (*basePhotons[i]).mass() );
    fBaseIDPhoton_iso_gamma.push_back( 0.0 );
    fBaseIDPhoton_iso_ch.push_back( 0.0 );
    fBaseIDPhoton_HE.push_back( 0.0 );
    fBaseIDPhoton_sigmaieie.push_back( 0.0 );
  }
  for (unsigned int i = 0; i < loosePhotons.size(); i++ )
  {
    fLooseIDPhoton_pt.push_back( (*loosePhotons[i]).pt() );
    fLooseIDPhoton_eta.push_back( (*loosePhotons[i]).eta() );
    fLooseIDPhoton_phi.push_back( (*loosePhotons[i]).phi() );
    fLooseIDPhoton_mass.push_back( (*loosePhotons[i]).mass() );
    fLooseIDPhoton_iso_gamma.push_back( 0.0 );
  }
  for (unsigned int i = 0; i < goodPhotons.size(); i++ )
  {
    fIDPhoton_pt.push_back( (*goodPhotons[i]).pt() );
    fIDPhoton_eta.push_back( (*goodPhotons[i]).eta() );
    fIDPhoton_phi.push_back( (*goodPhotons[i]).phi() );
    fIDPhoton_mass.push_back( (*goodPhotons[i]).mass() );
  }
  if (fAddDrConePhotonCut) {
    fNumTightPhotons_ConeHE = goodPhotons_ConeHE.size();
    for (unsigned int i = 0; i < goodPhotons_ConeHE.size(); i++ )
    {
      fID2Photon_pt.push_back( (*goodPhotons_ConeHE[i]).pt() );
      fID2Photon_eta.push_back( (*goodPhotons_ConeHE[i]).eta() );
      fID2Photon_phi.push_back( (*goodPhotons_ConeHE[i]).phi() );
      fID2Photon_mass.push_back( (*goodPhotons_ConeHE[i]).mass() );
    }
  }
  // fill old style structs, keeping because of extra info stored there in case needed
  InitRecoPhotonInfo(fRecoTightPhotonInfo1);
  InitRecoPhotonInfo(fRecoTightPhotonInfo2);
  InitRecoPhotonInfo(fRecoTightPhotonInfo3);
  if (goodPhotons.size() > 0)
    ExoDiPhotons::FillRecoPhotonInfo(fRecoTightPhotonInfo1,&(*goodPhotons[0]),lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent,iSetup); 
  if (goodPhotons.size() > 1)
    ExoDiPhotons::FillRecoPhotonInfo(fRecoTightPhotonInfo2,&(*goodPhotons[1]),lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent,iSetup); 
  if (goodPhotons.size() > 2)
    ExoDiPhotons::FillRecoPhotonInfo(fRecoTightPhotonInfo3,&(*goodPhotons[2]),lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent,iSetup); 
  if (fDebug) cout << ". done high-pt-id photons" << endl;

  // Construct Di-Objects
  InitRecoDiObjectInfo(fRecoPhiDiTwoProng);
  // Di-TwoProng
  if (fNumTwoProngPass >= 2)
  {
    TLorentzVector Eta1;
    Eta1.SetPtEtaPhiM(fTwoProng_pt[0], fTwoProng_eta[0], fTwoProng_phi[0], fTwoProng_mass[0]);
    TLorentzVector Eta2;
    Eta2.SetPtEtaPhiM(fTwoProng_pt[1], fTwoProng_eta[1], fTwoProng_phi[1], fTwoProng_mass[1]);
    FillRecoDiObjectInfo(fRecoPhiDiTwoProng, Eta1, Eta2);
  }
  InitRecoDiObjectInfo(fRecoPhiPhotonTwoProng);
  if (fDebug) cout << ". done making di-twoprong" << endl;
  // photon plus TwoProng
  if (fNumTwoProngPass >= 1 && fNumTightPhotons >= 1)
  {
    TLorentzVector LeadingTwoProng;
    LeadingTwoProng.SetPtEtaPhiM(fTwoProng_pt[0], fTwoProng_eta[0], fTwoProng_phi[0], fTwoProng_mass[0]);
    TLorentzVector LeadingPhoton;
    LeadingPhoton.SetPtEtaPhiM(fIDPhoton_pt[0], fIDPhoton_eta[0], fIDPhoton_phi[0], fIDPhoton_mass[0]);
    FillRecoDiObjectInfo(fRecoPhiPhotonTwoProng, LeadingTwoProng, LeadingPhoton);
  }
  InitRecoDiObjectInfo(fRecoPhiInclusive);
  if (fDebug) cout << ". done making photon-twoprong" << endl;
  // TwoProng plus (photon or TwoProng) inclusive
  if (fNumTwoProngPass >=1 && (fNumTwoProngPass + fNumTightPhotons >= 2))
  {
    TLorentzVector LeadingTwoProng;
    LeadingTwoProng.SetPtEtaPhiM(fTwoProng_pt[0], fTwoProng_eta[0], fTwoProng_phi[0], fTwoProng_mass[0]);
    TLorentzVector SubLeadingTwoProng;
    TLorentzVector LeadingPhoton;
    TLorentzVector LeadingSecondary;
    if (fNumTightPhotons >= 1) LeadingPhoton.SetPtEtaPhiM(fIDPhoton_pt[0], fIDPhoton_eta[0], fIDPhoton_phi[0], fIDPhoton_mass[0]);
    if (fNumTwoProngPass >= 2) SubLeadingTwoProng.SetPtEtaPhiM(fTwoProng_pt[1], fTwoProng_eta[1], fTwoProng_phi[1], fTwoProng_mass[1]);
  
    if (fNumTightPhotons == 0) LeadingSecondary = SubLeadingTwoProng;
    if (fNumTwoProngPass == 1) LeadingSecondary = LeadingPhoton;
    if (fNumTwoProngPass >= 2 && fNumTightPhotons >=1) {
      if (SubLeadingTwoProng.Pt() > LeadingPhoton.Pt()) LeadingSecondary = SubLeadingTwoProng;
      else LeadingSecondary = LeadingPhoton;
    }
    FillRecoDiObjectInfo(fRecoPhiInclusive, LeadingTwoProng, LeadingSecondary);
  }
  if (fDebug) cout << ". done making reco phi inclusive" << endl;

  // Now fill fTree2
  if (fMakeTrees) fTree2->Fill();

  // Optional Analysis code here

  // Photon Trigger Eff
  if (fTriggerEffHistos) {
  if (fNumTightPhotons > 0) {
    fPhotonTriggerEff_all_Denominator->Fill( (*goodPhotons[0]).pt() );     
    if(found_175 && found_22_iso && (fHLT_Photon175) )
        fPhotonTriggerEff_Photon175_Numerator->Fill( (*goodPhotons[0]).pt() );
    if(found_175 && found_22_iso && (fHLT_Photon22_Iso) )
        fPhotonTriggerEff_Photon22_Iso_Numerator->Fill( (*goodPhotons[0]).pt() );
    if(found_175 && found_22_iso && (fHLT_Photon175 || fHLT_Photon22_Iso) )
        fPhotonTriggerEff_all_Numerator->Fill( (*goodPhotons[0]).pt() );
  }
  if (fAddDrConePhotonCut) {
    if (fNumTightPhotons_ConeHE > 0) {
      fPhotonTriggerEff_ConeHE_all_Denominator->Fill( (*goodPhotons_ConeHE[0]).pt() );     
      if(found_175 && found_22_iso && (fHLT_Photon175) )
          fPhotonTriggerEff_ConeHE_Photon175_Numerator->Fill( (*goodPhotons_ConeHE[0]).pt() );
      if(found_175 && found_22_iso && (fHLT_Photon22_Iso) )
          fPhotonTriggerEff_ConeHE_Photon22_Iso_Numerator->Fill( (*goodPhotons_ConeHE[0]).pt() );
      if(found_175 && found_22_iso && (fHLT_Photon175 || fHLT_Photon22_Iso) )
          fPhotonTriggerEff_ConeHE_all_Numerator->Fill( (*goodPhotons_ConeHE[0]).pt() );
    }
  }
  if (fDebug) cout << ". done photon trigger efficiency histograms" << endl;
  }

  if (fTwoProngYieldHistos) {

  // twoprong yield analysis
  for (unsigned int i = 0; i < fTwoProng_pt.size(); i++) {
    fTwoProngYield->Fill(fTwoProng_pt[i]);
  }

  }

  if (fStackedDalitzHistos) {

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
ExoDiPhotonAnalyzer::beginJob()
{
  if (fMakeTrees) {
  cout << "===========================" << endl;
  cout << "= Ntuplizer Configuration =" << endl;
  cout << "===========================" << endl;
  cout << "Debug " << fDebug << endl;
  cout << "AddDrConePhotonCut " << fAddDrConePhotonCut << endl;
  cout << "includeSignalGenParticles " << fincludeSignalGenParticles << endl;
  cout << "includeAllLooseObjects " << fincludeAllLooseObjects << endl;
  cout << "includeAllCandObjects " << fincludeAllCandObjects << endl;
  cout << "includeOldPhotons " << fincludeOldPhotons << endl;
  cout << "includeMCInfo " << fincludeMCInfo << endl;
  cout << "McXS " << fMcXS << endl;
  cout << "McN " << fMcN << endl;
  cout << "MakeTrees " << fMakeTrees << endl;
  cout << "FakeRateHistos " << fFakeRateHistos << endl;
  cout << "TriggerEffHistos " << fTriggerEffHistos << endl;
  cout << "TwoProngYieldHistos " << fTwoProngYieldHistos << endl;
  cout << "StackedDalitzHistos " << fStackedDalitzHistos << endl;
  cout << "===========================" << endl;
  cout << "CandidatePairDR " << fCandidatePairDR << endl;
  cout << "CandidatePairMinPt " << fCandidatePairMinPt << endl;
  cout << "CandidatePairIsolationDR " << fCandidatePairIsolationDR << endl;
  cout << "CandidatePairPhiBox " << fCandidatePairPhiBox << endl;
  cout << "CandidatePairEtaBox " << fCandidatePairEtaBox << endl;
  cout << "CandidatePairPhotonPtCut " << fCandidatePairPhotonPtCut << endl;
  cout << "CandidatePairChargedIsoCut " << fCandidatePairChargedIsoCut << endl;
  cout << "CandidatePairNeutralIsoCut " << fCandidatePairNeutralIsoCut << endl;
  cout << "CandidatePairEGammaIsoCut " << fCandidatePairEGammaIsoCut << endl;
  cout << "CandidatePairChargedIsoFakeCut " << fCandidatePairChargedIsoFakeCut << endl;
  cout << "CandidatePairNeutralIsoFakeCut " << fCandidatePairNeutralIsoFakeCut << endl;
  cout << "CandidatePairEGammaIsoFakeCut " << fCandidatePairEGammaIsoFakeCut << endl;
  cout << "CandidatePairGenMatchDR " << fCandidatePairGenMatchDR << endl;
  cout << "CandidateAbsMaxEta " << fCandidateAbsMaxEta << endl;
  cout << "CandidateMinPt " << fCandidateMinPt << endl;
  cout << "CandidateTrackAsymmetryCut " << fCandidateTrackAsymmetryCut << endl;
  cout << "CandidatePhotonAsymmetryCut " << fCandidatePhotonAsymmetryCut << endl;
  cout << "CandidateOptionalExtraTrack " << fCandidateOptionalExtraTrack << endl;
  cout << "===========================" << endl;
  }
}

void 
ExoDiPhotonAnalyzer::endJob()
{
  // fake rate histograms
  if (fFakeRateHistos) {
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
  }

  // photon trigg eff histograms
  if (fTriggerEffHistos) {
  fPhotonTriggerEff_all_Denominator->Sumw2();
  fPhotonTriggerEff_all_Numerator->Sumw2();
  fPhotonTriggerEff_all_Division->Add(fPhotonTriggerEff_all_Numerator);
  fPhotonTriggerEff_all_Division->Divide(fPhotonTriggerEff_all_Denominator);

  fPhotonTriggerEff_Photon175_Numerator->Sumw2();
  fPhotonTriggerEff_Photon175_Division->Add(fPhotonTriggerEff_Photon175_Numerator);
  fPhotonTriggerEff_Photon175_Division->Divide(fPhotonTriggerEff_all_Denominator);

  fPhotonTriggerEff_Photon22_Iso_Numerator->Sumw2();
  fPhotonTriggerEff_Photon22_Iso_Division->Add(fPhotonTriggerEff_Photon22_Iso_Numerator);
  fPhotonTriggerEff_Photon22_Iso_Division->Divide(fPhotonTriggerEff_all_Denominator);

  if (fAddDrConePhotonCut) {
  fPhotonTriggerEff_ConeHE_all_Denominator->Sumw2();
  fPhotonTriggerEff_ConeHE_all_Numerator->Sumw2();
  fPhotonTriggerEff_ConeHE_all_Division->Add(fPhotonTriggerEff_ConeHE_all_Numerator);
  fPhotonTriggerEff_ConeHE_all_Division->Divide(fPhotonTriggerEff_ConeHE_all_Denominator);

  fPhotonTriggerEff_ConeHE_Photon175_Numerator->Sumw2();
  fPhotonTriggerEff_ConeHE_Photon175_Division->Add(fPhotonTriggerEff_ConeHE_Photon175_Numerator);
  fPhotonTriggerEff_ConeHE_Photon175_Division->Divide(fPhotonTriggerEff_ConeHE_all_Denominator);

  fPhotonTriggerEff_ConeHE_Photon22_Iso_Numerator->Sumw2();
  fPhotonTriggerEff_ConeHE_Photon22_Iso_Division->Add(fPhotonTriggerEff_ConeHE_Photon22_Iso_Numerator);
  fPhotonTriggerEff_ConeHE_Photon22_Iso_Division->Divide(fPhotonTriggerEff_ConeHE_all_Denominator);
  }
  }

  // stacke dalitz histograms
  if (fStackedDalitzHistos) {

  fTripleStackedDalitz->Add(fHighvsMid);
  fTripleStackedDalitz->Add(fHighvsLow);
  fTripleStackedDalitz->Add(fMidvsLow);
  fDoubleStackedDalitz->Add(fPhotonvsLarger);
  fDoubleStackedDalitz->Add(fPhotonvsSmaller);

  }
}

// High-pt Id subroutines

// determine if saturated, needed for high-pt-id
bool ExoDiPhotonAnalyzer::photon_isSaturated(const pat::Photon *photon, const EcalRecHitCollection *recHitsEB, const EcalRecHitCollection *recHitsEE,
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

// The high-pt-id
bool ExoDiPhotonAnalyzer::photon_passHighPtID(const pat::Photon* photon, double rho, bool isSat)
{
  if (
    passHadTowerOverEmCut(photon) &&
    passChargedHadronCut(photon) &&
    passSigmaIetaIetaCut(photon,isSat) &&
    passCorPhoIsoHighPtID(photon,rho) &&
    photon->passElectronVeto()
  ) return true;

  else return false;
}
// The high-pt-id without iso gamma
bool ExoDiPhotonAnalyzer::photon_passHighPtID_loose(const pat::Photon* photon, double rho, bool isSat)
{
  if (
    passHadTowerOverEmCut(photon) &&
    passChargedHadronCut(photon) &&
    passSigmaIetaIetaCut(photon,isSat) &&
    photon->passElectronVeto()
  ) return true;

  else return false;
}
// The minimal high-pt-id
bool ExoDiPhotonAnalyzer::photon_passHighPtID_base(const pat::Photon* photon, double rho, bool isSat)
{
  if (
    photon->passElectronVeto()
  ) return true;

  else return false;
}
// The modified high-pt-id
bool ExoDiPhotonAnalyzer::photon_passHighPtID_DrCone(const pat::Photon* photon, double rho, bool isSat)
{
  if (
    passHadTowerOverEmCut(photon) &&
    passHadDrConeOverEmCut(photon) && // new cut
    passChargedHadronCut(photon) &&
    passSigmaIetaIetaCut(photon,isSat) &&
    passCorPhoIsoHighPtID(photon,rho) &&
    photon->passElectronVeto()
  ) return true;

  else return false;
}

// photon subroutines as global functions, eventually move to private class functions

// H/E TOWER (part of high-pt-id)
bool passHadTowerOverEmCut(const pat::Photon* photon)
{
  double hOverE = photon->hadTowOverEm();
  if (hOverE < 0.05) return true;
  else return false;
}

// H/E DR CONE (method of trigger)
bool passHadDrConeOverEmCut(const pat::Photon* photon)
{
  double hOverE = photon->hadronicOverEm();
  if (hOverE < 0.05) return true;
  else return false;
}

// CH ISO
bool passChargedHadronCut(const pat::Photon* photon)
{
  double chIsoCut = 5.;
  double chIso = photon->chargedHadronIso();
  if (chIso < chIsoCut) return true;
  else return false;
}

// SIGMAiETAiETA
bool passSigmaIetaIetaCut(const pat::Photon* photon, bool isSaturated)
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

// COR ISO
bool passCorPhoIsoHighPtID(const pat::Photon* photon, double rho)
{
  double phoEta = fabs(photon->superCluster()->eta());
  double corPhoIsoCut = -999.9;
  double corPhoIso = corPhoIsoHighPtID(photon,rho);

  if (phoEta < 1.4442) corPhoIsoCut = 2.75;
  if (1.566 < phoEta && phoEta < 2.5) corPhoIsoCut = 2.00;

  if (corPhoIso < corPhoIsoCut) return true;
  else return false;
}

double corPhoIsoHighPtID(const pat::Photon* photon, double rho)
{
  double phoIso = photon->photonIso();
  return (phoAlphaHighPtID(photon) + phoIso - rho*phoEAHighPtID(photon) - phoKappaHighPtID(photon)*photon->pt());
}

double phoAlphaHighPtID(const pat::Photon *photon)
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

double phoEAHighPtID(const pat::Photon* photon)
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

double phoKappaHighPtID(const pat::Photon *photon)
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

// function for finding Ztoll decay type
vector<string>
ExoDiPhotonAnalyzer::getDecay(const reco::Candidate & genparticle, int flag)
{
  vector<string> products;

  if (flag == 1) { // parent is tau
    if (genparticle.pdgId() == 111) { products.push_back("npion"); return products; }
    if (abs(genparticle.pdgId()) == 211) { products.push_back("cpion"); return products; }
  }

  // ignore quarks and gluons and their daughters
  if (genparticle.pdgId() >= 1 && genparticle.pdgId() <= 8) { return products; }
  if (genparticle.pdgId() == 21) { return products; }

  if (genparticle.pdgId() == 11) { products.push_back("e-"); return products; }
  if (genparticle.pdgId() == -11) { products.push_back("e+"); return products; }
  if (genparticle.pdgId() == 13) { products.push_back("mu-"); return products; }
  if (genparticle.pdgId() == -13) { products.push_back("mu+"); return products; }

  if (abs(genparticle.pdgId()) == 15) {
    vector<string> tau_products; 
    for (unsigned int j = 0; j < genparticle.numberOfDaughters(); j++) {
      const reco::Candidate* daughter = genparticle.daughter(j);
      vector<string> daughter_products = getDecay(*daughter, 1);
      tau_products.insert(tau_products.end(), daughter_products.begin(), daughter_products.end());
    }
  
    if (tau_products.size() == 0) {
      if (genparticle.pdgId() == 15) { products.push_back("tau-unid"); return products; }
      if (genparticle.pdgId() == -15) { products.push_back("tau+unid"); return products; }
    }

    if (tau_products.size() == 1 && tau_products[0] == "e-") { products.push_back("tau-e"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "e+") { products.push_back("tau+e"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "mu-") { products.push_back("tau-mu"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "mu+") { products.push_back("tau+mu"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau-e") { products.push_back("tau-e"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau+e") { products.push_back("tau+e"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau-mu") { products.push_back("tau-mu"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau+mu") { products.push_back("tau+mu"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau+had10") { products.push_back("tau+had10"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau+had1") { products.push_back("tau+had1"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau+had30") { products.push_back("tau+had30"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau+had3") { products.push_back("tau+had3"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau-had10") { products.push_back("tau-had10"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau-had1") { products.push_back("tau-had1"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau-had30") { products.push_back("tau-had30"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau-had3") { products.push_back("tau-had3"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau+unid") { products.push_back("tau+unid"); return products; }
    if (tau_products.size() == 1 && tau_products[0] == "tau-unid") { products.push_back("tau-unid"); return products; }

    for (string prod : tau_products) {
      if (prod != "cpion" && prod != "npion") {  
        if (genparticle.pdgId() == 15) { products.push_back("tau-unid"); return products; }
        if (genparticle.pdgId() == -15) { products.push_back("tau+unid"); return products; }
      }
    }

    bool neutral = false;
    if (std::find(tau_products.begin(), tau_products.end(), "npion") != tau_products.end() ) neutral = true;

    int pion_count = 0;
    while( std::find(tau_products.begin(), tau_products.end(), "cpion") != tau_products.end() ) 
    {
      tau_products.erase( std::find(tau_products.begin(), tau_products.end(), "cpion") );
      pion_count += 1;
    }

    int charge = 0;
    if (genparticle.pdgId() == 15) charge = -1;
    if (genparticle.pdgId() == -15) charge = 1;

    if      (charge>0 && neutral  && pion_count == 1) products.push_back("tau+had10");
    else if (charge>0 && !neutral && pion_count == 1) products.push_back("tau+had1");
    else if (charge>0 && neutral  && pion_count == 3) products.push_back("tau+had30");
    else if (charge>0 && !neutral && pion_count == 3) products.push_back("tau+had3");
    else if (charge<0 && neutral  && pion_count == 1) products.push_back("tau-had10");
    else if (charge<0 && !neutral && pion_count == 1) products.push_back("tau-had1");
    else if (charge<0 && neutral  && pion_count == 3) products.push_back("tau-had30");
    else if (charge<0 && !neutral && pion_count == 3) products.push_back("tau-had3");
    else if (charge>0)                                products.push_back("tau+unid");
    else if (charge<0)                                products.push_back("tau-unid");

    return products;
  }

  for (unsigned int j = 0; j < genparticle.numberOfDaughters(); j++) {
    const reco::Candidate* daughter = genparticle.daughter(j);
    vector<string> daughter_products = getDecay(*daughter);
    products.insert(products.end(), daughter_products.begin(), daughter_products.end());
  } 
  return products;
}


// global sorting by pt function
bool compareCandsByPt(const edm::Ptr<const reco::Candidate> cand1, const edm::Ptr<const reco::Candidate> cand2)
{
  return(cand1->pt() >= cand2->pt());
}

//define this as a plug-in
DEFINE_FWK_MODULE(ExoDiPhotonAnalyzer);
