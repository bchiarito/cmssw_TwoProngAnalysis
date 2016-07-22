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

using namespace std;


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
  edm::InputTag      fPhotonTag;       //select photon collection 
  double             fMin_pt;          // min pt cut (photons)
  edm::InputTag      fHltInputTag;     // hltResults
  // edm::InputTag      fL1InputTag;      // L1 results
  edm::InputTag      fRho25Tag;  
  edm::InputTag      fpileupCollectionTag;
  edm::InputTag      fJetCollTag;         
  edm::LumiReWeighting    LumiWeights;      
  
  bool               fkRemoveSpikes;   // option to remove spikes before filling tree
  bool               fkRequireTightPhotons;  // option to require tight photon id in tree
  bool               fkRequireGenEventInfo;  // generated information for RS graviton files
  bool               fisMC;  // option to decide if MC or Data     
  bool               fisSignal;  // option to decide if Signal MC or Backround MC
  bool               fDebug;  // if set to False, mean to limit per event stdout output
  bool               fchargedDecayCutflow;  // option to print cutflow of the charged selection to stdout at the end of the job
  bool               fisAOD;
  string             fPUMCFileName;
  string             fPUDataFileName;
  string             fPUDataHistName;
  string             fPUMCHistName;   
  string             fPFIDCategory;   
  string             fIDMethod;   

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


  // to get L1 info, the L1 guide recommends to make this a member
  // this allows the event setup parts to be cached, rather than refetched every event
  // L1GtUtils m_l1GtUtils;

  // my Tree
  TTree *fTree;
  TTree *fTree2;

  // now also tight-fake and fake-fake trees
  TTree *fTightFakeTree;
  TTree *fFakeTightTree;
  TTree *fFakeFakeTree;

  ExoDiPhotons::eventInfo_t fEventInfo;
  ExoDiPhotons::vtxInfo_t fVtxInfo;
  // adding a second vertex
  ExoDiPhotons::vtxInfo_t fVtx2Info;
  // now even adding a third vtx!
  ExoDiPhotons::vtxInfo_t fVtx3Info;

  ExoDiPhotons::vtxInfo_t fVtxGENInfo;

      
  double fRho25;
  int fpu_n;
  int fold_pu_n; //Pileupbefore Bunch Crossing Correction
  int fBC;
  double fMCPUWeight;
  Int_t gv_n;
  
  TClonesArray* gv_pos;
  TClonesArray* gv_p3;
  
  Float_t gv_sumPtHi[100];
  Float_t gv_sumPtLo[100];
  Short_t gv_nTkHi[100];
  Short_t gv_nTkLo[100];

  ExoDiPhotons::beamSpotInfo_t fBeamSpotInfo;

  ExoDiPhotons::l1TrigInfo_t fL1TrigInfo;  
  ExoDiPhotons::hltTrigInfo_t fHLTInfo;

  int fNTightPhotons; // number of candidate photons in event (ie tight)
  int fNFakeablePhotons;  // number of 'fakeable objects' in event

  // //for PFIsolation Code
  // PFIsolationEstimator isolator04;
  // PFIsolationEstimator isolator03;
  // PFIsolationEstimator isolator02;

  ExoDiPhotons::recoPhotonInfo_t fRecoPhotonInfo1; // leading photon 
  ExoDiPhotons::recoPhotonInfo_t fRecoPhotonInfo2; // second photon

  ExoDiPhotons::jetInfo_t fJetInfo;
  ExoDiPhotons::conversionInfo_t fConvInfo1;
  ExoDiPhotons::conversionInfo_t fConvInfo2;
  ExoDiPhotons::conversionInfo_t fConvInfo_OneLeg1;
  ExoDiPhotons::conversionInfo_t fConvInfo_OneLeg2;
   
  ExoDiPhotons::diphotonInfo_t fDiphotonInfo;

  // diphoton info based on using hte second or third vtx in event
  ExoDiPhotons::diphotonInfo_t fDiphotonInfoVtx2; 
  ExoDiPhotons::diphotonInfo_t fDiphotonInfoVtx3; 


  TH1F* fpu_n_BeforeCuts; 
  TH1F* fpu_n_BeforeCutsAfterReWeight;
  TH1F *fNumTotalEvents;
  TH1F *fNumTotalWeightedEvents;

  // ** charged decay analysis **
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfcandsToken_;
  edm::EDGetTokenT<vector<reco::GenParticle>> genToken_;
  edm::EDGetTokenT<vector<reco::Vertex>> pvToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamToken_;
  edm::EDGetTokenT<std::vector<pat::Jet>> ak4Token_;
  double fJetPtCut;
  double fJetEtaCut;
  double fCandidatePairDR;
  double fCandidatePairMinPt;
  double fCandidatePairIsolationDR;
  double fCandidatePairPhiBox;
  double fCandidatePairEtaBox;
  double fCandidatePairPhotonPtCut;
  double fCandidatePairChargedIsoCut;
  double fCandidatePairNeutralIsoCut;
  double fCandidatePairEGammaIsoCut;
  double fCandidatePairGenMatchDR;
  // Fake rate histos
  TH1F *fTwoProngFakeRate_pt;
  TH1F *fTwoProngFakeRate_eta;
  TH1F *fTwoProngFakeRate_phi;
  TH1F *fTwoProngFakeNumer_pt;
  TH1F *fTwoProngFakeDenom_pt;
  TH1F *fTwoProngFakeNumer_eta;
  TH1F *fTwoProngFakeDenom_eta;
  TH1F *fTwoProngFakeNumer_phi;
  TH1F *fTwoProngFakeDenom_phi;

  //-----------------taken from Ilya-----------------
  // Format-independent data members
  edm::EDGetTokenT<double> rhoToken_;
  
  // AOD case data members
  edm::EDGetToken photonsToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
  
  // MiniAOD case data members
  edm::EDGetToken photonsMiniAODToken_;
  edm::EDGetToken patPhotonToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesMiniAODToken_;
  
  // Photon variables computed upstream in a special producer
  edm::EDGetTokenT<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMapToken_; 
  edm::EDGetTokenT<edm::ValueMap<float> > phoChargedIsolationToken_; 
  edm::EDGetTokenT<edm::ValueMap<float> > phoNeutralHadronIsolationToken_; 
  edm::EDGetTokenT<edm::ValueMap<float> > phoPhotonIsolationToken_; 

  // ID decision objects
  edm::EDGetTokenT<edm::ValueMap<bool> > phoLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > phoMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > phoTightIdMapToken_;

  // vertex token
  edm::EDGetToken vtxToken_;

  // beamspot
  edm::EDGetToken bsToken_;

  // trigger
  edm::EDGetToken trigToken_;

  // jets
  edm::EDGetToken jetsToken_;

  // conversions
  edm::EDGetToken twoLegToken_;
  edm::EDGetToken oneLegToken_;

  Float_t rho_;      // the rho variable

  // Effective area constants for all isolation types
  EffectiveAreas effAreaChHadrons_;
  EffectiveAreas effAreaNeuHadrons_;
  EffectiveAreas effAreaPhotons_;

  std::vector<Int_t> passLooseId_;
  std::vector<Int_t> passMediumId_;
  std::vector<Int_t> passTightId_;
  //-----------------taken from Ilya-----------------
  
  // ** charged decay analysis **
  int fEventNum;
  int fRunNum;
  int fLumiNum;
  int fNumPVs;
  int fNumPF;
  int fNumPrunedPF;
  int fNumCHpairs;
  int fNumCHpairsPass;
  int fNumCHpairsPassChargedIso;
  int fNumCHpairsPassNeutralIso;
  int fNumCHpairsPassEGammaIso;
  int fNumCHpairsPassPhotonPtIso;
  int fNumCHpairsFake;
  int fNumTightPhotons;
  double fHT;
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
  vector<Double_t> fCand_CHpos_dxy;
  vector<Double_t> fCand_CHpos_dxy_beamspot;
  vector<Double_t> fCand_CHpos_dxy_associated;
  vector<Double_t> fCand_CHneg_dxy;
  vector<Double_t> fCand_CHneg_dxy_beamspot;
  vector<Double_t> fCand_CHneg_dxy_associated;
  vector<Bool_t> fCand_pass;
  vector<Bool_t> fCand_passChargedIso;
  vector<Bool_t> fCand_passNeutralIso;
  vector<Bool_t> fCand_passEGammaIso;
  vector<Bool_t> fCand_passPhotonPt;
  vector<Bool_t> fCand_fake;
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
  vector<Double_t> fTwoProng_photon_pt;
  vector<Double_t> fTwoProng_photon_eta;
  vector<Double_t> fTwoProng_photon_phi;
  vector<Double_t> fTwoProng_photon_mass;
  vector<Double_t> fTwoProng_photon_nGamma;
  vector<Double_t> fTwoProng_photon_nElectron;
  vector<Double_t> fTwoProng_chargedIso;
  vector<Double_t> fTwoProng_neutralIso;
  vector<Double_t> fTwoProng_egammaIso;
  vector<Double_t> fTwoProng_CHpos_dxy;
  vector<Double_t> fTwoProng_CHpos_dxy_beamspot;
  vector<Double_t> fTwoProng_CHpos_dxy_associated;
  vector<Double_t> fTwoProng_CHneg_dxy;
  vector<Double_t> fTwoProng_CHneg_dxy_beamspot;
  vector<Double_t> fTwoProng_CHneg_dxy_associated;
  vector<Bool_t> fTwoProng_match;
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
  int fGen_decayType;
  ExoDiPhotons::recoPhotonInfo_t fRecoTightPhotonInfo1;
  ExoDiPhotons::recoPhotonInfo_t fRecoTightPhotonInfo2;
  ExoDiPhotons::recoPhotonInfo_t fRecoTightPhotonInfo3; 
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
  : fPhotonTag(iConfig.getUntrackedParameter<edm::InputTag>("photonCollection")),
    fMin_pt(iConfig.getUntrackedParameter<double>("ptMin")),
    fHltInputTag(iConfig.getUntrackedParameter<edm::InputTag>("hltResults")),
    // fL1InputTag(iConfig.getUntrackedParameter<edm::InputTag>("L1Results")),
    fRho25Tag(iConfig.getParameter<edm::InputTag>("rho25Correction")),
    fpileupCollectionTag(iConfig.getUntrackedParameter<edm::InputTag>("pileupCorrection")),
    fkRemoveSpikes(iConfig.getUntrackedParameter<bool>("removeSpikes")),
    fkRequireTightPhotons(iConfig.getUntrackedParameter<bool>("requireTightPhotons")),
    fkRequireGenEventInfo(iConfig.getUntrackedParameter<bool>("requireGenEventInfo")),
    fisMC(iConfig.getUntrackedParameter<bool>("isMC")),
    fisSignal(iConfig.getUntrackedParameter<bool>("isSignal")),
    fDebug(iConfig.getUntrackedParameter<bool>("debug")),
    fchargedDecayCutflow(iConfig.getUntrackedParameter<bool>("chargedDecayCutflow")),
    fisAOD(iConfig.getParameter<bool>("isAOD")),
    fPUMCFileName(iConfig.getUntrackedParameter<string>("PUMCFileName")),
    fPUDataFileName(iConfig.getUntrackedParameter<string>("PUDataFileName")), 
    fPUDataHistName(iConfig.getUntrackedParameter<string>("PUDataHistName")),
    fPUMCHistName(iConfig.getUntrackedParameter<string>("PUMCHistName")),
    fPFIDCategory(iConfig.getUntrackedParameter<string>("PFIDCategory")),
    fIDMethod(iConfig.getUntrackedParameter<string>("IDMethod")),
    fJetPtCut(iConfig.getUntrackedParameter<double>("jetPtCut")),
    fJetEtaCut(iConfig.getUntrackedParameter<double>("jetEtaCut")),
    fCandidatePairDR(iConfig.getUntrackedParameter<double>("chargedHadronPairMinDeltaR")),
    fCandidatePairMinPt(iConfig.getUntrackedParameter<double>("chargedHadronMinPt")),
    fCandidatePairIsolationDR(iConfig.getUntrackedParameter<double>("isolationConeR")),
    fCandidatePairPhiBox(iConfig.getUntrackedParameter<double>("photonPhiBoxSize")),
    fCandidatePairEtaBox(iConfig.getUntrackedParameter<double>("photonEtaBoxSize")),
    fCandidatePairPhotonPtCut(iConfig.getUntrackedParameter<double>("photonPtCut")),
    fCandidatePairChargedIsoCut(iConfig.getUntrackedParameter<double>("chargedIsoCut")),
    fCandidatePairNeutralIsoCut(iConfig.getUntrackedParameter<double>("neutralIsoCut")),
    fCandidatePairEGammaIsoCut(iConfig.getUntrackedParameter<double>("egammaIsoCut")),
    fCandidatePairGenMatchDR(iConfig.getUntrackedParameter<double>("generatorEtaMatchDR")),
    //-----------------taken from Ilya-----------------
    rhoToken_(consumes<double> (iConfig.getParameter<edm::InputTag>("rho"))),
    // Cluster shapes
    full5x5SigmaIEtaIEtaMapToken_(consumes <edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("full5x5SigmaIEtaIEtaMap"))),
    // Isolations
    phoChargedIsolationToken_(consumes <edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("phoChargedIsolation"))),
    phoNeutralHadronIsolationToken_(consumes <edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("phoNeutralHadronIsolation"))),
    phoPhotonIsolationToken_(consumes <edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("phoPhotonIsolation"))),
    phoLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("phoLooseIdMap"))),
    phoMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("phoMediumIdMap"))),
    phoTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("phoTightIdMap"))),
    // Objects containing effective area constants
    effAreaChHadrons_( (iConfig.getParameter<edm::FileInPath>("effAreaChHadFile")).fullPath() ),
    effAreaNeuHadrons_( (iConfig.getParameter<edm::FileInPath>("effAreaNeuHadFile")).fullPath() ),
    effAreaPhotons_( (iConfig.getParameter<edm::FileInPath>("effAreaPhoFile")).fullPath() )
    //-----------------taken from Ilya-----------------
{
  //now do what ever initialization is needed
  
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

  // LumiReweighting Tool

  std::cout<<"ExoDiPhotonAnalyzer: ID Method used "<<fIDMethod.c_str()
     <<"PF ID Category "<<fPFIDCategory.c_str()
     <<std::endl;

  //-----------------taken from Ilya-----------------
  //
  // Prepare tokens for all input collections and objects
  //
  // AOD tokens
  photonsToken_ = mayConsume<edm::View<reco::Photon> >
    (iConfig.getParameter<edm::InputTag>
     ("photons"));
  
  genParticlesToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticles"));
  
  // MiniAOD tokens
  photonsMiniAODToken_ = mayConsume<edm::View<reco::Photon> >
    (iConfig.getParameter<edm::InputTag>
     ("photonsMiniAOD"));
  patPhotonToken_ = mayConsume<edm::View<pat::Photon> >
    (iConfig.getParameter<edm::InputTag>
     ("photonsMiniAOD"));
  
  genParticlesMiniAODToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticlesMiniAOD"));

  if (fisAOD) vtxToken_ = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));
  else vtxToken_ = consumes<reco::VertexCollection>(edm::InputTag("offlineSlimmedPrimaryVertices"));

  if (!fisAOD) pfcandsToken_ = consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates"));
  if (!fisAOD && fisMC) genToken_ = consumes<vector<reco::GenParticle>>(edm::InputTag("prunedGenParticles"));
  if (!fisAOD) pvToken_ = consumes<vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));
  if (!fisAOD) beamToken_ = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  if (!fisAOD) ak4Token_ = consumes<std::vector<pat::Jet>>(edm::InputTag("slimmedJets"));

  bsToken_ = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  trigToken_ = consumes<edm::TriggerResults>(fHltInputTag);
  fJetCollTag = iConfig.getParameter<edm::InputTag>("jetCollection");
  if (!fisAOD) jetsToken_ = consumes< edm::View<pat::Jet> >(fJetCollTag);

  twoLegToken_ = consumes<reco::ConversionCollection>(edm::InputTag("reducedEgamma","reducedConversions"));
  oneLegToken_ = consumes<reco::ConversionCollection>(edm::InputTag("reducedEgamma","reducedSingleLegConversions"));
  //-----------------taken from Ilya-----------------
 
  edm::Service<TFileService> fs;
  fpu_n_BeforeCuts = fs->make<TH1F>("fpu_n_BeforeCuts","PileUpBeforeCuts",300,0,300);
  fpu_n_BeforeCutsAfterReWeight = fs->make<TH1F>("fpu_n_BeforeCutsAfterReWeight","PileUpBeforeCuts",300,0,300);
  fNumTotalEvents = fs->make<TH1F>("NumTotalEvents","Total number of events",4,0.,2.);
  fNumTotalWeightedEvents = fs->make<TH1F>("NumTotalWeightedEvents","Total weighted number of events",4,0.,2.);
  fTree = fs->make<TTree>("fTree","PhotonTree");

  fTree->Branch("Event",&fEventInfo,ExoDiPhotons::eventInfoBranchDefString.c_str());
  fTree->Branch("Vtx",&fVtxInfo,ExoDiPhotons::vtxInfoBranchDefString.c_str());
  //adding a second vtx
  fTree->Branch("Vtx2",&fVtx2Info,ExoDiPhotons::vtxInfoBranchDefString.c_str());
  fTree->Branch("Vtx3",&fVtx3Info,ExoDiPhotons::vtxInfoBranchDefString.c_str());
  fTree->Branch("VtxGEN",&fVtxGENInfo,ExoDiPhotons::vtxInfoBranchDefString.c_str());

  fTree->Branch("rho25",&fRho25,"rho25/D"); 
  fTree->Branch("pu_n", &fpu_n, "pu_n/I");
  fTree->Branch("old_pu_n", &fold_pu_n, "old_pu_n/I");
  fTree->Branch("MCPUWeight",&fMCPUWeight,"MCPUWeight/D");
  fTree->Branch("BeamSpot",&fBeamSpotInfo,ExoDiPhotons::beamSpotInfoBranchDefString.c_str());
  fTree->Branch("L1trg",&fL1TrigInfo,ExoDiPhotons::l1TrigBranchDefString.c_str());
  fTree->Branch("TrigHLT",&fHLTInfo,ExoDiPhotons::hltTrigBranchDefString.c_str());
  // add a branch for number of candidate photons in the event (tight and fakeable)
  fTree->Branch("nTightPhotons",&fNTightPhotons,"nTightPhotons/I");
  fTree->Branch("nFakeablePhotons",&fNFakeablePhotons,"nFakeablePhotons/I");
  // now with CommonClasses, use the string defined in the header
  fTree->Branch("Photon1",&fRecoPhotonInfo1,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  fTree->Branch("Photon2",&fRecoPhotonInfo2,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  fTree->Branch("JetInfo.pt","std::vector<double>",&fJetInfo.pt);
  fTree->Branch("JetInfo.eta","std::vector<double>",&fJetInfo.eta);
  fTree->Branch("JetInfo.phi","std::vector<double>",&fJetInfo.phi);
  fTree->Branch("JetInfo.mass","std::vector<double>",&fJetInfo.mass);
  fTree->Branch("JetInfo.energy","std::vector<double>",&fJetInfo.energy);
  fTree->Branch("JetInfo.passLooseID","std::vector<int>",&fJetInfo.passLooseID);
  fTree->Branch("JetInfo.passTightID","std::vector<int>",&fJetInfo.passTightID);
  fTree->Branch("JetInfo.HT",&fJetInfo.HT,"JetInfo.HT/D");
  fTree->Branch("JetInfo.missingHT",&fJetInfo.missingHT,"JetInfo.missingHT/D");

  fTree->Branch("ConvInfo1.x","std::vector<double>",&fConvInfo1.x);
  fTree->Branch("ConvInfo1.y","std::vector<double>",&fConvInfo1.y);
  fTree->Branch("ConvInfo1.z","std::vector<double>",&fConvInfo1.z);
  fTree->Branch("ConvInfo1.r","std::vector<double>",&fConvInfo1.r);
  fTree->Branch("ConvInfo1.phi","std::vector<double>",&fConvInfo1.phi);
  fTree->Branch("ConvInfo1.dPhiTracksAtVtx","std::vector<double>",&fConvInfo1.dPhiTracksAtVtx);
  fTree->Branch("ConvInfo1.nTracks","std::vector<double>",&fConvInfo1.nTracks);
  fTree->Branch("ConvInfo1.dxy","std::vector<double>",&fConvInfo1.dxy);
  fTree->Branch("ConvInfo1.dz","std::vector<double>",&fConvInfo1.dz);
  fTree->Branch("ConvInfo1.vtxChi2","std::vector<double>",&fConvInfo1.vtxChi2);
  fTree->Branch("ConvInfo1.vtxNdof","std::vector<double>",&fConvInfo1.vtxNdof);
  fTree->Branch("ConvInfo1.pairCotThetaSeparation","std::vector<double>",&fConvInfo1.pairCotThetaSeparation);
  fTree->Branch("ConvInfo1.photonPt","std::vector<double>",&fConvInfo1.photonPt);
  fTree->Branch("ConvInfo1.track1InnerPx","std::vector<double>",&fConvInfo1.track1InnerPx);
  fTree->Branch("ConvInfo1.track1InnerPy","std::vector<double>",&fConvInfo1.track1InnerPy);
  fTree->Branch("ConvInfo1.track1InnerPz","std::vector<double>",&fConvInfo1.track1InnerPz);
  fTree->Branch("ConvInfo1.track2InnerPx","std::vector<double>",&fConvInfo1.track2InnerPx);
  fTree->Branch("ConvInfo1.track2InnerPy","std::vector<double>",&fConvInfo1.track2InnerPy);
  fTree->Branch("ConvInfo1.track2InnerPz","std::vector<double>",&fConvInfo1.track2InnerPz);

  fTree->Branch("ConvInfo1.dRToSc","std::vector<double>",&fConvInfo1.dRToSc);
  // fTree->Branch("ConvInfo1.nHitsTrack1","std::vector<int>",&fConvInfo1.nHitsTrack1);
  // fTree->Branch("ConvInfo1.nHitsTrack2","std::vector<int>",&fConvInfo1.nHitsTrack2);
  fTree->Branch("ConvInfo1.nSharedHits","std::vector<uint8_t>",&fConvInfo1.nSharedHits);
  fTree->Branch("ConvInfo1.MVAout","std::vector<double>",&fConvInfo1.MVAout);
  fTree->Branch("ConvInfo1.oneLegMVA","std::vector<std::vector<float>>",&fConvInfo1.oneLegMVA);
  fTree->Branch("ConvInfo1.nHitsBeforeVtx","std::vector<std::vector<uint8_t>>",&fConvInfo1.nHitsBeforeVtx);
  fTree->Branch("ConvInfo1.isUndefinedAlgo","std::vector<int>",&fConvInfo1.isUndefinedAlgo);
  fTree->Branch("ConvInfo1.isEcalSeededAlgo","std::vector<int>",&fConvInfo1.isEcalSeededAlgo);
  fTree->Branch("ConvInfo1.isTrackerOnlyAlgo","std::vector<int>",&fConvInfo1.isTrackerOnlyAlgo);
  fTree->Branch("ConvInfo1.isMixedAlgo","std::vector<int>",&fConvInfo1.isMixedAlgo);
  fTree->Branch("ConvInfo1.isPflowAlgo","std::vector<int>",&fConvInfo1.isPflowAlgo);
  fTree->Branch("ConvInfo1.isGeneralTracksOnly","std::vector<int>",&fConvInfo1.isGeneralTracksOnly);
  fTree->Branch("ConvInfo1.isArbitratedEcalSeeded","std::vector<int>",&fConvInfo1.isArbitratedEcalSeeded);
  fTree->Branch("ConvInfo1.isArbitratedMerged","std::vector<int>",&fConvInfo1.isArbitratedMerged);
  fTree->Branch("ConvInfo1.isArbitratedMergedEcalGeneral","std::vector<int>",&fConvInfo1.isArbitratedMergedEcalGeneral);
  fTree->Branch("ConvInfo1.isHighPurity","std::vector<int>",&fConvInfo1.isHighPurity);
  fTree->Branch("ConvInfo1.isHighEfficiency","std::vector<int>",&fConvInfo1.isHighEfficiency);
  fTree->Branch("ConvInfo1.isEcalMatched1Track","std::vector<int>",&fConvInfo1.isEcalMatched1Track);
  fTree->Branch("ConvInfo1.isEcalMatched2Track","std::vector<int>",&fConvInfo1.isEcalMatched2Track);

  fTree->Branch("ConvInfo2.x","std::vector<double>",&fConvInfo2.x);
  fTree->Branch("ConvInfo2.y","std::vector<double>",&fConvInfo2.y);
  fTree->Branch("ConvInfo2.z","std::vector<double>",&fConvInfo2.z);
  fTree->Branch("ConvInfo2.r","std::vector<double>",&fConvInfo2.r);
  fTree->Branch("ConvInfo2.phi","std::vector<double>",&fConvInfo2.phi);
  fTree->Branch("ConvInfo2.dPhiTracksAtVtx","std::vector<double>",&fConvInfo2.dPhiTracksAtVtx);
  fTree->Branch("ConvInfo2.nTracks","std::vector<double>",&fConvInfo2.nTracks);
  fTree->Branch("ConvInfo2.dxy","std::vector<double>",&fConvInfo2.dxy);
  fTree->Branch("ConvInfo2.dz","std::vector<double>",&fConvInfo2.dz);
  fTree->Branch("ConvInfo2.vtxChi2","std::vector<double>",&fConvInfo2.vtxChi2);
  fTree->Branch("ConvInfo2.vtxNdof","std::vector<double>",&fConvInfo2.vtxNdof);
  fTree->Branch("ConvInfo2.pairCotThetaSeparation","std::vector<double>",&fConvInfo2.pairCotThetaSeparation);
  fTree->Branch("ConvInfo2.photonPt","std::vector<double>",&fConvInfo2.photonPt);
  fTree->Branch("ConvInfo2.track1InnerPx","std::vector<double>",&fConvInfo2.track1InnerPx);
  fTree->Branch("ConvInfo2.track1InnerPy","std::vector<double>",&fConvInfo2.track1InnerPy);
  fTree->Branch("ConvInfo2.track1InnerPz","std::vector<double>",&fConvInfo2.track1InnerPz);
  fTree->Branch("ConvInfo2.track2InnerPx","std::vector<double>",&fConvInfo2.track2InnerPx);
  fTree->Branch("ConvInfo2.track2InnerPy","std::vector<double>",&fConvInfo2.track2InnerPy);
  fTree->Branch("ConvInfo2.track2InnerPz","std::vector<double>",&fConvInfo2.track2InnerPz);

  fTree->Branch("ConvInfo2.dRToSc","std::vector<double>",&fConvInfo2.dRToSc);
  // fTree->Branch("ConvInfo2.nHitsTrack1","std::vector<int>",&fConvInfo2.nHitsTrack1);
  // fTree->Branch("ConvInfo2.nHitsTrack2","std::vector<int>",&fConvInfo2.nHitsTrack2);
  fTree->Branch("ConvInfo2.nSharedHits","std::vector<uint8_t>",&fConvInfo2.nSharedHits);
  fTree->Branch("ConvInfo2.MVAout","std::vector<double>",&fConvInfo2.MVAout);
  fTree->Branch("ConvInfo2.oneLegMVA","std::vector<std::vector<float>>",&fConvInfo2.oneLegMVA);
  fTree->Branch("ConvInfo2.nHitsBeforeVtx","std::vector<std::vector<uint8_t>>",&fConvInfo2.nHitsBeforeVtx);
  fTree->Branch("ConvInfo2.isUndefinedAlgo","std::vector<int>",&fConvInfo2.isUndefinedAlgo);
  fTree->Branch("ConvInfo2.isEcalSeededAlgo","std::vector<int>",&fConvInfo2.isEcalSeededAlgo);
  fTree->Branch("ConvInfo2.isTrackerOnlyAlgo","std::vector<int>",&fConvInfo2.isTrackerOnlyAlgo);
  fTree->Branch("ConvInfo2.isMixedAlgo","std::vector<int>",&fConvInfo2.isMixedAlgo);
  fTree->Branch("ConvInfo2.isPflowAlgo","std::vector<int>",&fConvInfo2.isPflowAlgo);
  fTree->Branch("ConvInfo2.isGeneralTracksOnly","std::vector<int>",&fConvInfo2.isGeneralTracksOnly);
  fTree->Branch("ConvInfo2.isArbitratedEcalSeeded","std::vector<int>",&fConvInfo2.isArbitratedEcalSeeded);
  fTree->Branch("ConvInfo2.isArbitratedMerged","std::vector<int>",&fConvInfo2.isArbitratedMerged);
  fTree->Branch("ConvInfo2.isArbitratedMergedEcalGeneral","std::vector<int>",&fConvInfo2.isArbitratedMergedEcalGeneral);
  fTree->Branch("ConvInfo2.isHighPurity","std::vector<int>",&fConvInfo2.isHighPurity);
  fTree->Branch("ConvInfo2.isHighEfficiency","std::vector<int>",&fConvInfo2.isHighEfficiency);
  fTree->Branch("ConvInfo2.isEcalMatched1Track","std::vector<int>",&fConvInfo2.isEcalMatched1Track);
  fTree->Branch("ConvInfo2.isEcalMatched2Track","std::vector<int>",&fConvInfo2.isEcalMatched2Track);
   
  fTree->Branch("ConvInfo_OneLeg1.x","std::vector<double>",&fConvInfo_OneLeg1.x);
  fTree->Branch("ConvInfo_OneLeg1.y","std::vector<double>",&fConvInfo_OneLeg1.y);
  fTree->Branch("ConvInfo_OneLeg1.z","std::vector<double>",&fConvInfo_OneLeg1.z);
  fTree->Branch("ConvInfo_OneLeg1.r","std::vector<double>",&fConvInfo_OneLeg1.r);
  fTree->Branch("ConvInfo_OneLeg1.phi","std::vector<double>",&fConvInfo_OneLeg1.phi);
  fTree->Branch("ConvInfo_OneLeg1.dPhiTracksAtVtx","std::vector<double>",&fConvInfo_OneLeg1.dPhiTracksAtVtx);
  fTree->Branch("ConvInfo_OneLeg1.nTracks","std::vector<double>",&fConvInfo_OneLeg1.nTracks);
  fTree->Branch("ConvInfo_OneLeg1.dxy","std::vector<double>",&fConvInfo_OneLeg1.dxy);
  fTree->Branch("ConvInfo_OneLeg1.dz","std::vector<double>",&fConvInfo_OneLeg1.dz);
  fTree->Branch("ConvInfo_OneLeg1.vtxChi2","std::vector<double>",&fConvInfo_OneLeg1.vtxChi2);
  fTree->Branch("ConvInfo_OneLeg1.vtxNdof","std::vector<double>",&fConvInfo_OneLeg1.vtxNdof);
  fTree->Branch("ConvInfo_OneLeg1.pairCotThetaSeparation","std::vector<double>",&fConvInfo_OneLeg1.pairCotThetaSeparation);
  fTree->Branch("ConvInfo_OneLeg1.photonPt","std::vector<double>",&fConvInfo_OneLeg1.photonPt);
  fTree->Branch("ConvInfo_OneLeg1.track1InnerPx","std::vector<double>",&fConvInfo_OneLeg1.track1InnerPx);
  fTree->Branch("ConvInfo_OneLeg1.track1InnerPy","std::vector<double>",&fConvInfo_OneLeg1.track1InnerPy);
  fTree->Branch("ConvInfo_OneLeg1.track1InnerPz","std::vector<double>",&fConvInfo_OneLeg1.track1InnerPz);
  fTree->Branch("ConvInfo_OneLeg1.track2InnerPx","std::vector<double>",&fConvInfo_OneLeg1.track2InnerPx);
  fTree->Branch("ConvInfo_OneLeg1.track2InnerPy","std::vector<double>",&fConvInfo_OneLeg1.track2InnerPy);
  fTree->Branch("ConvInfo_OneLeg1.track2InnerPz","std::vector<double>",&fConvInfo_OneLeg1.track2InnerPz);

  fTree->Branch("ConvInfo_OneLeg1.dRToSc","std::vector<double>",&fConvInfo_OneLeg1.dRToSc);
  // fTree->Branch("ConvInfo_OneLeg1.nHitsTrack1","std::vector<int>",&fConvInfo_OneLeg1.nHitsTrack1);
  // fTree->Branch("ConvInfo_OneLeg1.nHitsTrack2","std::vector<int>",&fConvInfo_OneLeg1.nHitsTrack2);
  fTree->Branch("ConvInfo_OneLeg1.nSharedHits","std::vector<uint8_t>",&fConvInfo_OneLeg1.nSharedHits);
  fTree->Branch("ConvInfo_OneLeg1.MVAout","std::vector<double>",&fConvInfo_OneLeg1.MVAout);
  fTree->Branch("ConvInfo_OneLeg1.oneLegMVA","std::vector<std::vector<float>>",&fConvInfo_OneLeg1.oneLegMVA);
  fTree->Branch("ConvInfo_OneLeg1.nHitsBeforeVtx","std::vector<std::vector<uint8_t>>",&fConvInfo_OneLeg1.nHitsBeforeVtx);
  fTree->Branch("ConvInfo_OneLeg1.isUndefinedAlgo","std::vector<int>",&fConvInfo_OneLeg1.isUndefinedAlgo);
  fTree->Branch("ConvInfo_OneLeg1.isEcalSeededAlgo","std::vector<int>",&fConvInfo_OneLeg1.isEcalSeededAlgo);
  fTree->Branch("ConvInfo_OneLeg1.isTrackerOnlyAlgo","std::vector<int>",&fConvInfo_OneLeg1.isTrackerOnlyAlgo);
  fTree->Branch("ConvInfo_OneLeg1.isMixedAlgo","std::vector<int>",&fConvInfo_OneLeg1.isMixedAlgo);
  fTree->Branch("ConvInfo_OneLeg1.isPflowAlgo","std::vector<int>",&fConvInfo_OneLeg1.isPflowAlgo);
  fTree->Branch("ConvInfo_OneLeg1.isGeneralTracksOnly","std::vector<int>",&fConvInfo_OneLeg1.isGeneralTracksOnly);
  fTree->Branch("ConvInfo_OneLeg1.isArbitratedEcalSeeded","std::vector<int>",&fConvInfo_OneLeg1.isArbitratedEcalSeeded);
  fTree->Branch("ConvInfo_OneLeg1.isArbitratedMerged","std::vector<int>",&fConvInfo_OneLeg1.isArbitratedMerged);
  fTree->Branch("ConvInfo_OneLeg1.isArbitratedMergedEcalGeneral","std::vector<int>",&fConvInfo_OneLeg1.isArbitratedMergedEcalGeneral);
  fTree->Branch("ConvInfo_OneLeg1.isHighPurity","std::vector<int>",&fConvInfo_OneLeg1.isHighPurity);
  fTree->Branch("ConvInfo_OneLeg1.isHighEfficiency","std::vector<int>",&fConvInfo_OneLeg1.isHighEfficiency);
  fTree->Branch("ConvInfo_OneLeg1.isEcalMatched1Track","std::vector<int>",&fConvInfo_OneLeg1.isEcalMatched1Track);
  fTree->Branch("ConvInfo_OneLeg1.isEcalMatched2Track","std::vector<int>",&fConvInfo_OneLeg1.isEcalMatched2Track);
   
  fTree->Branch("ConvInfo_OneLeg2.x","std::vector<double>",&fConvInfo_OneLeg2.x);
  fTree->Branch("ConvInfo_OneLeg2.y","std::vector<double>",&fConvInfo_OneLeg2.y);
  fTree->Branch("ConvInfo_OneLeg2.z","std::vector<double>",&fConvInfo_OneLeg2.z);
  fTree->Branch("ConvInfo_OneLeg2.r","std::vector<double>",&fConvInfo_OneLeg2.r);
  fTree->Branch("ConvInfo_OneLeg2.phi","std::vector<double>",&fConvInfo_OneLeg2.phi);
  fTree->Branch("ConvInfo_OneLeg2.dPhiTracksAtVtx","std::vector<double>",&fConvInfo_OneLeg2.dPhiTracksAtVtx);
  fTree->Branch("ConvInfo_OneLeg2.nTracks","std::vector<double>",&fConvInfo_OneLeg2.nTracks);
  fTree->Branch("ConvInfo_OneLeg2.dxy","std::vector<double>",&fConvInfo_OneLeg2.dxy);
  fTree->Branch("ConvInfo_OneLeg2.dz","std::vector<double>",&fConvInfo_OneLeg2.dz);
  fTree->Branch("ConvInfo_OneLeg2.vtxChi2","std::vector<double>",&fConvInfo_OneLeg2.vtxChi2);
  fTree->Branch("ConvInfo_OneLeg2.vtxNdof","std::vector<double>",&fConvInfo_OneLeg2.vtxNdof);
  fTree->Branch("ConvInfo_OneLeg2.pairCotThetaSeparation","std::vector<double>",&fConvInfo_OneLeg2.pairCotThetaSeparation);
  fTree->Branch("ConvInfo_OneLeg2.photonPt","std::vector<double>",&fConvInfo_OneLeg2.photonPt);
  fTree->Branch("ConvInfo_OneLeg2.track1InnerPx","std::vector<double>",&fConvInfo_OneLeg2.track1InnerPx);
  fTree->Branch("ConvInfo_OneLeg2.track1InnerPy","std::vector<double>",&fConvInfo_OneLeg2.track1InnerPy);
  fTree->Branch("ConvInfo_OneLeg2.track1InnerPz","std::vector<double>",&fConvInfo_OneLeg2.track1InnerPz);
  fTree->Branch("ConvInfo_OneLeg2.track2InnerPx","std::vector<double>",&fConvInfo_OneLeg2.track2InnerPx);
  fTree->Branch("ConvInfo_OneLeg2.track2InnerPy","std::vector<double>",&fConvInfo_OneLeg2.track2InnerPy);
  fTree->Branch("ConvInfo_OneLeg2.track2InnerPz","std::vector<double>",&fConvInfo_OneLeg2.track2InnerPz);

  fTree->Branch("ConvInfo_OneLeg2.dRToSc","std::vector<double>",&fConvInfo_OneLeg2.dRToSc);
  // fTree->Branch("ConvInfo_OneLeg2.nHitsTrack1","std::vector<int>",&fConvInfo_OneLeg2.nHitsTrack1);
  // fTree->Branch("ConvInfo_OneLeg2.nHitsTrack2","std::vector<int>",&fConvInfo_OneLeg2.nHitsTrack2);
  fTree->Branch("ConvInfo_OneLeg2.nSharedHits","std::vector<uint8_t>",&fConvInfo_OneLeg2.nSharedHits);
  fTree->Branch("ConvInfo_OneLeg2.MVAout","std::vector<double>",&fConvInfo_OneLeg2.MVAout);
  fTree->Branch("ConvInfo_OneLeg2.oneLegMVA","std::vector<std::vector<float>>",&fConvInfo_OneLeg2.oneLegMVA);
  fTree->Branch("ConvInfo_OneLeg2.nHitsBeforeVtx","std::vector<std::vector<uint8_t>>",&fConvInfo_OneLeg2.nHitsBeforeVtx);
  fTree->Branch("ConvInfo_OneLeg2.isUndefinedAlgo","std::vector<int>",&fConvInfo_OneLeg2.isUndefinedAlgo);
  fTree->Branch("ConvInfo_OneLeg2.isEcalSeededAlgo","std::vector<int>",&fConvInfo_OneLeg2.isEcalSeededAlgo);
  fTree->Branch("ConvInfo_OneLeg2.isTrackerOnlyAlgo","std::vector<int>",&fConvInfo_OneLeg2.isTrackerOnlyAlgo);
  fTree->Branch("ConvInfo_OneLeg2.isMixedAlgo","std::vector<int>",&fConvInfo_OneLeg2.isMixedAlgo);
  fTree->Branch("ConvInfo_OneLeg2.isPflowAlgo","std::vector<int>",&fConvInfo_OneLeg2.isPflowAlgo);
  fTree->Branch("ConvInfo_OneLeg2.isGeneralTracksOnly","std::vector<int>",&fConvInfo_OneLeg2.isGeneralTracksOnly);
  fTree->Branch("ConvInfo_OneLeg2.isArbitratedEcalSeeded","std::vector<int>",&fConvInfo_OneLeg2.isArbitratedEcalSeeded);
  fTree->Branch("ConvInfo_OneLeg2.isArbitratedMerged","std::vector<int>",&fConvInfo_OneLeg2.isArbitratedMerged);
  fTree->Branch("ConvInfo_OneLeg2.isArbitratedMergedEcalGeneral","std::vector<int>",&fConvInfo_OneLeg2.isArbitratedMergedEcalGeneral);
  fTree->Branch("ConvInfo_OneLeg2.isHighPurity","std::vector<int>",&fConvInfo_OneLeg2.isHighPurity);
  fTree->Branch("ConvInfo_OneLeg2.isHighEfficiency","std::vector<int>",&fConvInfo_OneLeg2.isHighEfficiency);
  fTree->Branch("ConvInfo_OneLeg2.isEcalMatched1Track","std::vector<int>",&fConvInfo_OneLeg2.isEcalMatched1Track);
  fTree->Branch("ConvInfo_OneLeg2.isEcalMatched2Track","std::vector<int>",&fConvInfo_OneLeg2.isEcalMatched2Track);


  fTree->Branch("Diphoton",&fDiphotonInfo,ExoDiPhotons::diphotonInfoBranchDefString.c_str());
  // diphoton info for second or thrid best vertex
  // only bothering to add this for tight-tight tree for now
  fTree->Branch("DiphotonVtx2",&fDiphotonInfoVtx2,ExoDiPhotons::diphotonInfoBranchDefString.c_str());
  fTree->Branch("DiphotonVtx3",&fDiphotonInfoVtx3,ExoDiPhotons::diphotonInfoBranchDefString.c_str());

  // Branches for charged decay analysis
  fTree2 = fs->make<TTree>("fTree2","ChargedDecayTree");
  // Event wide
  fTree2->Branch("eventNum",&fEventNum,"eventNum/I");
  fTree2->Branch("runNum",&fRunNum,"runNum/I");
  fTree2->Branch("lumiNum",&fLumiNum,"lumiNum/I");
  fTree2->Branch("nPV",&fNumPVs,"nPV/I");
  fTree2->Branch("nPF",&fNumPF,"nPF/I");
  fTree2->Branch("nPrunedPF",&fNumPrunedPF,"numPrunedPF/I");
  fTree2->Branch("HT",&fHT,"HT/D");
  // Cutflow
  fTree2->Branch("nCand",&fNumCHpairs,"nCand/I");
  fTree2->Branch("nPass",&fNumCHpairsPass,"nPass/I");
  fTree2->Branch("nPassChargedIso",&fNumCHpairsPassChargedIso,"nPassChargedIso/I");
  fTree2->Branch("nPassNeutralIso",&fNumCHpairsPassNeutralIso,"nPassNeutralIso/I");
  fTree2->Branch("nPassEGammaIso",&fNumCHpairsPassEGammaIso,"nPassEGammaIso/I");
  fTree2->Branch("nPassPhotonPtIso",&fNumCHpairsPassChargedIso,"nPassPhotonIso/I");
  fTree2->Branch("nFake",&fNumCHpairsFake,"nFake/I");
  fTree2->Branch("nTightPhoton",&fNumTightPhotons,"nTightPhoton/I");
  // Candidate information
  fTree2->Branch("Cand_pt",&fCand_pt);
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
  fTree2->Branch("Cand_CHpos_dxy",&fCand_CHpos_dxy);
  fTree2->Branch("Cand_CHpos_dxy_beamspot",&fCand_CHpos_dxy_beamspot);
  fTree2->Branch("Cand_CHpos_dxy_associated",&fCand_CHpos_dxy_associated);
  fTree2->Branch("Cand_CHneg_dxy",&fCand_CHneg_dxy);
  fTree2->Branch("Cand_CHneg_dxy_beamspot",&fCand_CHneg_dxy_beamspot);
  fTree2->Branch("Cand_CHneg_dxy_associated",&fCand_CHneg_dxy_associated);
  fTree2->Branch("Cand_pass",&fCand_pass);
  fTree2->Branch("Cand_passChargedIso",&fCand_passChargedIso);
  fTree2->Branch("Cand_passNeutralIso",&fCand_passNeutralIso);
  fTree2->Branch("Cand_passEGammaIso",&fCand_passEGammaIso);
  fTree2->Branch("Cand_passPhotonPt",&fCand_passPhotonPt);
  fTree2->Branch("Cand_fake",&fCand_fake);
  // Passing Candidate information, sorted by pt
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
  fTree2->Branch("TwoProng_CHpos_dxy",&fTwoProng_CHpos_dxy);
  fTree2->Branch("TwoProng_CHpos_dxy_beamspot",&fTwoProng_CHpos_dxy_beamspot);
  fTree2->Branch("TwoProng_CHpos_dxy_associated",&fTwoProng_CHpos_dxy_associated);
  fTree2->Branch("TwoProng_CHneg_dxy",&fTwoProng_CHneg_dxy);
  fTree2->Branch("TwoProng_CHneg_dxy_beamspot",&fTwoProng_CHneg_dxy_beamspot);
  fTree2->Branch("TwoProng_CHneg_dxy_associated",&fTwoProng_CHneg_dxy_associated);
  // Tight Photons, sorted by pt
  fTree2->Branch("Photon1",&fRecoTightPhotonInfo1,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  fTree2->Branch("Photon2",&fRecoTightPhotonInfo2,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  fTree2->Branch("Photon3",&fRecoTightPhotonInfo3,ExoDiPhotons::recoPhotonBranchDefString.c_str());
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
  fTree2->Branch("Gen_decayType",&fGen_decayType);
  // Combined Objects
  fTree2->Branch("TwoProngTwoProng",&fTwoProngTwoProngInfo,ExoDiPhotons::recoDiObjectBranchDefString.c_str());
  fTree2->Branch("TwoProngFake",&fTwoProngFakeInfo,ExoDiPhotons::recoDiObjectBranchDefString.c_str());
  fTree2->Branch("GammaTwoProng",&fGammaTwoProngInfo,ExoDiPhotons::recoDiObjectBranchDefString.c_str());
  fTree2->Branch("GammaFake",&fGammaFakeInfo,ExoDiPhotons::recoDiObjectBranchDefString.c_str());
  fTree2->Branch("GammaGamma",&fGammaGammaInfo,ExoDiPhotons::recoDiObjectBranchDefString.c_str());
  // Fake rate histograms
  fTwoProngFakeNumer_pt = fs->make<TH1F>("twoprongfakenumer_pt","Fake Numerator count for CH pairs, inverted charged iso cut, pt binned",26,0.0,1300.0);
  fTwoProngFakeDenom_pt = fs->make<TH1F>("twoprongfakedenom_pt","Fake Denominator count for CH pairs, inverted charged iso cut, pt binned",26,0.0,1300.0);
  fTwoProngFakeRate_pt = fs->make<TH1F>("twoprongfake_pt","Fake Rate for CH pairs, inverted charged iso cut, pt binned",26,0.,1300);
  fTwoProngFakeNumer_eta = fs->make<TH1F>("twoprongfakenumer_eta","Fake Numerator count for CH pairs, inverted charged iso cut, eta binned",24,-6.0,6.0);
  fTwoProngFakeDenom_eta = fs->make<TH1F>("twoprongfakedenom_eta","Fake Denominator count for CH pairs, inverted charged iso cut, eta binned",24,-6.0,6.0);
  fTwoProngFakeRate_eta = fs->make<TH1F>("twoprongfake_eta","Fake Rate for CH pairs, inverted charged iso cut, eta binned",24,-6.0,6.0);
  fTwoProngFakeNumer_phi = fs->make<TH1F>("twoprongfakenumer_phi","Fake Numerator count for CH pairs, inverted charged iso cut, phi binned",24,-6.0,6.0);
  fTwoProngFakeDenom_phi = fs->make<TH1F>("twoprongfakedenom_phi","Fake Denominator count for CH pairs, inverted charged iso cut, phi binned",24,-6.0,6.0);
  fTwoProngFakeRate_phi = fs->make<TH1F>("twoprongfake_phi","Fake Rate for CH pairs, inverted charged iso cut",24,-6.0,6.0);

  // repeating all this for each of tight-fake and fake-fake trees
  // basically they'll all point to the same structs, but the structs will contain
  // different values for the event, depending on the event category

  fTightFakeTree = fs->make<TTree>("fTightFakeTree","PhotonTightFakeTree");
  fTightFakeTree->Branch("Event",&fEventInfo,ExoDiPhotons::eventInfoBranchDefString.c_str());
  fTightFakeTree->Branch("Vtx",&fVtxInfo,ExoDiPhotons::vtxInfoBranchDefString.c_str());
  fTightFakeTree->Branch("Vtx2",&fVtx2Info,ExoDiPhotons::vtxInfoBranchDefString.c_str());
  fTightFakeTree->Branch("BeamSpot",&fBeamSpotInfo,ExoDiPhotons::beamSpotInfoBranchDefString.c_str());
  fTightFakeTree->Branch("L1trg",&fL1TrigInfo,ExoDiPhotons::l1TrigBranchDefString.c_str());
  fTightFakeTree->Branch("TrigHLT",&fHLTInfo,ExoDiPhotons::hltTrigBranchDefString.c_str());
  fTightFakeTree->Branch("nTightPhotons",&fNTightPhotons,"nTightPhotons/I");
  fTightFakeTree->Branch("nFakeablePhotons",&fNFakeablePhotons,"nFakeablePhotons/I");
  fTightFakeTree->Branch("Photon1",&fRecoPhotonInfo1,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  fTightFakeTree->Branch("Photon2",&fRecoPhotonInfo2,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  fTightFakeTree->Branch("Diphoton",&fDiphotonInfo,ExoDiPhotons::diphotonInfoBranchDefString.c_str());

  fFakeTightTree = fs->make<TTree>("fFakeTightTree","PhotonFakeTightTree");

  fFakeTightTree->Branch("Event",&fEventInfo,ExoDiPhotons::eventInfoBranchDefString.c_str());
  fFakeTightTree->Branch("Vtx",&fVtxInfo,ExoDiPhotons::vtxInfoBranchDefString.c_str());
  //adding a second vtx                                                                                                                  
  fFakeTightTree->Branch("Vtx2",&fVtx2Info,ExoDiPhotons::vtxInfoBranchDefString.c_str());
  fFakeTightTree->Branch("Vtx3",&fVtx3Info,ExoDiPhotons::vtxInfoBranchDefString.c_str());
  fFakeTightTree->Branch("BeamSpot",&fBeamSpotInfo,ExoDiPhotons::beamSpotInfoBranchDefString.c_str());
  fFakeTightTree->Branch("L1trg",&fL1TrigInfo,ExoDiPhotons::l1TrigBranchDefString.c_str());
  fFakeTightTree->Branch("TrigHLT",&fHLTInfo,ExoDiPhotons::hltTrigBranchDefString.c_str());
  // add a branch for number of candidate photons in the event (tight and fakeable)                                                      
  fFakeTightTree->Branch("nTightPhotons",&fNTightPhotons,"nTightPhotons/I");
  fFakeTightTree->Branch("nFakeablePhotons",&fNFakeablePhotons,"nFakeablePhotons/I");
  // now with CommonClasses, use the string defined in the header                                                                        
  fFakeTightTree->Branch("Photon1",&fRecoPhotonInfo1,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  fFakeTightTree->Branch("Photon2",&fRecoPhotonInfo2,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  fFakeTightTree->Branch("Diphoton",&fDiphotonInfo,ExoDiPhotons::diphotonInfoBranchDefString.c_str());    
  fFakeTightTree->Branch("DiphotonVtx2",&fDiphotonInfoVtx2,ExoDiPhotons::diphotonInfoBranchDefString.c_str());
  fFakeTightTree->Branch("DiphotonVtx3",&fDiphotonInfoVtx3,ExoDiPhotons::diphotonInfoBranchDefString.c_str());

  fFakeFakeTree = fs->make<TTree>("fFakeFakeTree","PhotonFakeFakeTree");
  fFakeFakeTree->Branch("Event",&fEventInfo,ExoDiPhotons::eventInfoBranchDefString.c_str());
  fFakeFakeTree->Branch("Vtx",&fVtxInfo,ExoDiPhotons::vtxInfoBranchDefString.c_str());
  fFakeFakeTree->Branch("Vtx2",&fVtx2Info,ExoDiPhotons::vtxInfoBranchDefString.c_str());
  fFakeFakeTree->Branch("BeamSpot",&fBeamSpotInfo,ExoDiPhotons::beamSpotInfoBranchDefString.c_str());
  fFakeFakeTree->Branch("L1trg",&fL1TrigInfo,ExoDiPhotons::l1TrigBranchDefString.c_str());
  fFakeFakeTree->Branch("TrigHLT",&fHLTInfo,ExoDiPhotons::hltTrigBranchDefString.c_str());
  fFakeFakeTree->Branch("nTightPhotons",&fNTightPhotons,"nTightPhotons/I");
  fFakeFakeTree->Branch("nFakeablePhotons",&fNFakeablePhotons,"nFakeablePhotons/I");
  fFakeFakeTree->Branch("Photon1",&fRecoPhotonInfo1,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  fFakeFakeTree->Branch("Photon2",&fRecoPhotonInfo2,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  fFakeFakeTree->Branch("Diphoton",&fDiphotonInfo,ExoDiPhotons::diphotonInfoBranchDefString.c_str());
  

  gv_pos = new TClonesArray("TVector3", 100);
  gv_p3 = new TClonesArray("TVector3", 100);

  // //new PFIsolation code
  // isolator04.initializePhotonIsolation(kTRUE);
  // isolator04.setConeSize(0.4);
  // isolator03.initializePhotonIsolation(kTRUE);
  // isolator03.setConeSize(0.3);
  // isolator02.initializePhotonIsolation(kTRUE);
  // isolator02.setConeSize(0.2);

  if (fisAOD){
    recHitsEBTag_ = iConfig.getUntrackedParameter<edm::InputTag>("RecHitsEBTag",edm::InputTag("reducedEcalRecHitsEB"));
    recHitsEETag_ = iConfig.getUntrackedParameter<edm::InputTag>("RecHitsEETag",edm::InputTag("reducedEcalRecHitsEE"));
    recHitsEBToken = consumes <edm::SortedCollection<EcalRecHit> > (recHitsEBTag_);
    recHitsEEToken = consumes <edm::SortedCollection<EcalRecHit> > (recHitsEETag_);
  }
  else{
    recHitsEBTag_ = iConfig.getUntrackedParameter<edm::InputTag>("RecHitsEBTag",edm::InputTag("reducedEgamma:reducedEBRecHits"));
    recHitsEETag_ = iConfig.getUntrackedParameter<edm::InputTag>("RecHitsEETag",edm::InputTag("reducedEgamma:reducedEERecHits"));
    recHitsEBToken = consumes < EcalRecHitCollection > (recHitsEBTag_);
    recHitsEEToken = consumes < EcalRecHitCollection > (recHitsEETag_);
  }

  // if (fisAOD){
  //   recHitsEBTag_ = iConfig.getParameter<edm::InputTag>("recHitsEB");
  //   recHitsEETag_ = iConfig.getParameter<edm::InputTag>("recHitsEE");
  //   recHitsEBToken = consumes< edm::SortedCollection<EcalRecHit> >(recHitsEBTag_);
  //   recHitsEEToken = consumes< edm::SortedCollection<EcalRecHit> >(recHitsEETag_);
  // }
  // else{
  //   recHitsEBTag_ = iConfig.getParameter<edm::InputTag>("recHitsEBMiniAOD");
  //   recHitsEETag_ = iConfig.getParameter<edm::InputTag>("recHitsEEMiniAOD");
  //   recHitsEBToken = consumes< edm::SortedCollection<EcalRecHit> >(recHitsEBTag_);
  //   recHitsEEToken = consumes< edm::SortedCollection<EcalRecHit> >(recHitsEETag_);
  // }
  //TO DISENTANGLE BETWEEN MINIAOD AND AOD

}


ExoDiPhotonAnalyzer::~ExoDiPhotonAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

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

// ------------ method called to for each event  ------------
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

  // basic event info
  ExoDiPhotons::InitEventInfo(fEventInfo,-5000.);
  ExoDiPhotons::FillEventInfo(fEventInfo,iEvent);
   
  //-----------------taken from Ilya-----------------
  // Retrieve the collection of photons from the event.
  // If we fail to retrieve the collection with the standard AOD
  // name, we next look for the one with the stndard miniAOD name. 
  //   We use exactly the same handle for AOD and miniAOD formats
  // since pat::Photon objects can be recast as reco::Photon objects.
  edm::Handle<edm::View<reco::Photon> > photons;
  // bool fisAOD = true;
  if (fisAOD) iEvent.getByToken(photonsToken_, photons);
  else iEvent.getByToken(photonsMiniAODToken_,photons);
  if (fDebug) cout << "Found photon collection with size = " << (*photons).size() << endl;
  // if( !photons.isValid() ){
  //   // fisAOD = false;
  //   iEvent.getByToken(photonsMiniAODToken_,photons);
  //   cout << "Found MiniAOD photon collection with size = " << (*photons).size() << endl;
  // }

  //-----------------taken from Ilya-----------------
  // Get generator level info
  edm::Handle<edm::View<reco::GenParticle> > genParticles;
  if( fisAOD )
    iEvent.getByToken(genParticlesToken_,genParticles);
  else
    iEvent.getByToken(genParticlesMiniAODToken_,genParticles);
  
  //Get rho
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho_ = *rhoH;
  
  // Get the full5x5 map
  edm::Handle<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMap;
  iEvent.getByToken(full5x5SigmaIEtaIEtaMapToken_, full5x5SigmaIEtaIEtaMap);

  // Get the isolation maps
  edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
  iEvent.getByToken(phoChargedIsolationToken_, phoChargedIsolationMap);
  edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
  iEvent.getByToken(phoNeutralHadronIsolationToken_, phoNeutralHadronIsolationMap);
  edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
  iEvent.getByToken(phoPhotonIsolationToken_, phoPhotonIsolationMap);

  // Get the photon ID data from the event stream.
  // Note: this implies that the VID ID modules have been run upstream.
  // If you need more info, check with the EGM group.
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  iEvent.getByToken(phoLooseIdMapToken_ ,loose_id_decisions);
  iEvent.getByToken(phoMediumIdMapToken_,medium_id_decisions);
  iEvent.getByToken(phoTightIdMapToken_ ,tight_id_decisions);
  //-----------------taken from Ilya-----------------




  edm::Handle<GenEventInfoProduct> GenInfoHandle;
  if(fkRequireGenEventInfo){
    iEvent.getByLabel("generator",GenInfoHandle);
    if(!GenInfoHandle.isValid()) {
      cout << "Gen Event Info Product collection empty! Bailing out!" <<endl;
      return;
    }
  }

  if(fkRequireGenEventInfo) {
    fEventInfo.pthat = GenInfoHandle->hasBinningValues() ? (GenInfoHandle->binningValues())[0] : 0.0 ;
    fEventInfo.alphaqcd = GenInfoHandle->alphaQCD();
    fEventInfo.alphaqed = GenInfoHandle->alphaQED();
    fEventInfo.qscale = GenInfoHandle->qScale();
    fEventInfo.processid = GenInfoHandle->signalProcessID();
    fEventInfo.weight = GenInfoHandle->weights()[0];
  }

  fNumTotalEvents->Fill(1.);
  fNumTotalWeightedEvents->Fill(1.,fEventInfo.weight);

  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  // get the vertex collection
  Handle<reco::VertexCollection> vertexColl;
  iEvent.getByToken(vtxToken_,vertexColl);
  // if (fisAOD) iEvent.getByLabel("offlinePrimaryVertices",vertexColl);
  // else iEvent.getByLabel("offlineSlimmedPrimaryVertices",vertexColl);
  //TO DISENTANGLE BETWEEN MINIAOD AND AOD
   
  if(!vertexColl.isValid()) {
    cout << "Vertex collection empty! Bailing out!" <<endl;
    return;
  }

  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  fpu_n = -99999.99; 
  fold_pu_n = -99999.99;
  fBC = -99999.99;
  fMCPUWeight = -99999.99;
    

  fVtxInfo.vx = -99999.99;
  fVtxInfo.vy = -99999.99;
  fVtxInfo.vz = -99999.99;
  fVtxInfo.isFake = true;   
  fVtxInfo.Ntracks = -99;
  fVtxInfo.sumPtTracks = -99999.99;
  fVtxInfo.ndof = -99999.99;
  fVtxInfo.d0 = -99999.99;

  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  fVtx2Info.vx = -99999.99;
  fVtx2Info.vy = -99999.99;
  fVtx2Info.vz = -99999.99;
  fVtx2Info.isFake = true;   
  fVtx2Info.Ntracks = -99;
  fVtx2Info.sumPtTracks = -99999.99;
  fVtx2Info.ndof = -99999.99;
  fVtx2Info.d0 = -99999.99;


  fVtx3Info.vx = -99999.99;
  fVtx3Info.vy = -99999.99;
  fVtx3Info.vz = -99999.99;
  fVtx3Info.isFake = true;   
  fVtx3Info.Ntracks = -99;
  fVtx3Info.sumPtTracks = -99999.99;
  fVtx3Info.ndof = -99999.99;
  fVtx3Info.d0 = -99999.99;

  fVtxGENInfo.vx = -99999.99;
  fVtxGENInfo.vy = -99999.99;
  fVtxGENInfo.vz = -99999.99;
  fVtxGENInfo.isFake = true;
  fVtxGENInfo.Ntracks = -99;
  fVtxGENInfo.sumPtTracks = -99999.99;
  fVtxGENInfo.ndof = -99999.99;
  fVtxGENInfo.d0 = -99999.99;


  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;


  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  if (fisMC){
    edm::Handle<std::vector<PileupSummaryInfo> > pileupHandle;
    iEvent.getByLabel(fpileupCollectionTag, pileupHandle);
    std::vector<PileupSummaryInfo>::const_iterator PUI;
   
    if (pileupHandle.isValid()){
    
      for (PUI = pileupHandle->begin();PUI != pileupHandle->end(); ++PUI){
      
	fBC = PUI->getBunchCrossing() ;
	if(fBC==0){ 
	  //Select only the in time bunch crossing with bunch crossing=0
	  PileupSummaryInfo oldpileup = (*pileupHandle.product())[0];
	  fpu_n = PUI->getTrueNumInteractions();
	  fold_pu_n = oldpileup.getPU_NumInteractions();
	  fpu_n_BeforeCuts->Fill(fpu_n);
         
	}
      }
    
   
      fMCPUWeight = LumiWeights.weight(fpu_n);  
      fpu_n_BeforeCutsAfterReWeight->Fill(fpu_n,fMCPUWeight);
  
    }
  } 
  //add rho correction

  //      double rho;

  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  // edm::Handle<double> rho25Handle;
  // iEvent.getByLabel(fRho25Tag, rho25Handle);

  // if (!rho25Handle.isValid()){
  //   cout<<"rho25 not found"<<endl;
  //   return;
  // }
        
  // fRho25 = *(rho25Handle.product());
  fRho25 = *rhoH;

  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  //for conversion safe electron veto, only needed in AOD
  edm::Handle<reco::ConversionCollection> hConversions;
  if (fisAOD) iEvent.getByLabel("allConversions", hConversions);
  
  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  edm::Handle<reco::GsfElectronCollection> hElectrons;
  if (fisAOD) iEvent.getByLabel("gedGsfElectrons", hElectrons); // don't need these for MiniAOD
  //edm::Handle<pat::ElectronCollection> hElectrons;
  //iEvent.getByLabel(edm::InputTag("slimmedElectrons"), hElectrons);
  //patElectrons_slimmedElectrons__PAT.obj.embeddedSuperCluster_
  //TO DISENTANGLE BETWEEN MINIAOD AND AOD
  if(!hElectrons.isValid() && fisAOD) {
    cout<<"no ged gsf electrons "<<endl;
    return;
  }
  
  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  // //for PFIsolation code
  // Handle<PFCandidateCollection> pfCandidatesColl;
  // iEvent.getByLabel("particleFlow",pfCandidatesColl);
  // const PFCandidateCollection * pfCandidates = pfCandidatesColl.product();

  // get offline beam spot
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(bsToken_, beamSpotHandle);
   
  fBeamSpotInfo.x0 = -99999999.99;
  fBeamSpotInfo.y0 = -99999999.99;
  fBeamSpotInfo.z0 = -99999999.99;
  fBeamSpotInfo.sigmaZ = -99999999.99;
  fBeamSpotInfo.x0error = -99999999.99;
  fBeamSpotInfo.y0error = -99999999.99;
  fBeamSpotInfo.z0error = -99999999.99;
  fBeamSpotInfo.sigmaZ0error = -99999999.99;

  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  if(beamSpotHandle.isValid()) {
    beamSpot = *beamSpotHandle.product();
    //   cout << beamSpot <<endl;
    ExoDiPhotons::FillBeamSpotInfo(fBeamSpotInfo,beamSpot);
  }


  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  // get jet information 
  // RIGHT NOW ONLY FOR MINIAOD

  edm::Handle<edm::View<pat::Jet>> jetsHandle;
  iEvent.getByToken(jetsToken_,jetsHandle);
  const edm::View<pat::Jet>* jets = jetsHandle.product();

  // ExoDiPhotons::InitJetInfo(fJetInfo,jets);
  ExoDiPhotons::FillJetInfo(fJetInfo,jets,fJetPtCut,fJetEtaCut);

  // for (unsigned int i=0; i< jets->size(); i++){
  //   pat::Jet iJet = jets->at(i);

  //   // jetID
  //   bool isLoose;
  //   bool isTight;
  //   std::tie(isLoose,isTight) = jetID(iJet);

  // }

  // get the trig info

  //trig results
  Handle<TriggerResults> hltResultsHandle;
  iEvent.getByToken(trigToken_,hltResultsHandle);
   
  if(!hltResultsHandle.isValid()) {
    cout << "HLT results not valid!" <<endl;
    cout << "Couldnt find TriggerResults with input tag " << fHltInputTag << endl;
    return;
  }

  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  const TriggerResults *hltResults = hltResultsHandle.product();
  //   cout << *hltResults <<endl;
  const TriggerNames & hltNames = iEvent.triggerNames(*hltResults);
  //   TriggerNames hltNames;
  //   hltNames.init(*hltResults);
  //   cout << "HLT Results" <<endl;


  // now we just use the FillHLTInfo() function from TrigInfo.h:
  ExoDiPhotons::FillHLTInfo(fHLTInfo,hltResults,hltNames);



  // L1 results
  // fL1TrigInfo.L1_Tech0 = false;
  // fL1TrigInfo.L1_Tech36 = false;
  // fL1TrigInfo.L1_Tech37 = false;
  // fL1TrigInfo.L1_Tech38 = false;
  // fL1TrigInfo.L1_Tech39 = false;
  // fL1TrigInfo.L1_Tech40 = false;
  // fL1TrigInfo.L1_Tech41 = false;
  // fL1TrigInfo.L1_Tech42 = false;
  // fL1TrigInfo.L1_Tech43 = false;
  // fL1TrigInfo.L1_EG2 = false;   
   
  // use the L1GtUtils class, following instructions in 
  // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideL1TriggerL1GtUtils
   
  // before accessing any result from L1GtUtils, one must retrieve and cache
  // the L1 trigger event setup
  // m_l1GtUtils.retrieveL1EventSetup(iSetup);

  // access L1 trigger results using public methods from L1GtUtils
  // always check on error code returned by that method
  // int iErrorCode = -1;
  // error 0 is okay; 1 means trig not found; else some other error
  // although I am probably just going to ignore it anyway...
   
  // since many of these L1 bits are actually masked, I think I should use the 
  // decision word BEFORE the mask, to be useful for selection
   
  // the L1 bits/names association for the tech triggers of interest is:
  // L1Tech_BPTX_plus_AND_minus.v0  	0
  // L1Tech_BSC_halo_beam2_inner.v0  	36
  // L1Tech_BSC_halo_beam2_outer.v0 	37
  // L1Tech_BSC_halo_beam1_inner.v0 	38
  // L1Tech_BSC_halo_beam1_outer.v0 	39
  // L1Tech_BSC_minBias_threshold1.v0 	40
  // L1Tech_BSC_minBias_threshold2.v0 	41
  // L1Tech_BSC_splash_beam1.v0 	42
  // L1Tech_BSC_splash_beam2.v0 	43 	
  

   
  // fL1TrigInfo.L1_Tech0 =    m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BPTX_plus_AND_minus.v0",iErrorCode);
  // fL1TrigInfo.L1_Tech36 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_halo_beam2_inner.v0",iErrorCode);
  // fL1TrigInfo.L1_Tech37 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_halo_beam2_outer.v0",iErrorCode);
  // fL1TrigInfo.L1_Tech38 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_halo_beam1_inner.v0",iErrorCode);
  // fL1TrigInfo.L1_Tech39 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_halo_beam1_outer.v0",iErrorCode);
  // fL1TrigInfo.L1_Tech40 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_minBias_threshold1.v0",iErrorCode);
  // fL1TrigInfo.L1_Tech41 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_minBias_threshold2.v0",iErrorCode);
  // fL1TrigInfo.L1_Tech42 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_splash_beam1.v0",iErrorCode);
  // fL1TrigInfo.L1_Tech43 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_splash_beam2.v0",iErrorCode);


  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  // ecal information
  lazyTools_ = std::auto_ptr<noZS::EcalClusterLazyTools>( new noZS::EcalClusterLazyTools(iEvent, iSetup, recHitsEBToken, recHitsEEToken));   

  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  // get ecal barrel recHits for spike rejection
  edm::Handle<EcalRecHitCollection> recHitsEB_h;
  iEvent.getByToken(recHitsEBToken, recHitsEB_h );
  const EcalRecHitCollection * recHitsEB = 0;
  if ( ! recHitsEB_h.isValid() ) {
    LogError("ExoDiPhotonAnalyzer") << " ECAL Barrel RecHit Collection not available !"; return;
  } else {
    recHitsEB = recHitsEB_h.product();
  }
  //TO DISENTANGLE BETWEEN MINIAOD AND AOD

  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  edm::Handle<EcalRecHitCollection> recHitsEE_h;
  iEvent.getByToken(recHitsEEToken, recHitsEE_h );
  const EcalRecHitCollection * recHitsEE = 0;
  if ( ! recHitsEE_h.isValid() ) {
    LogError("ExoDiPhotonAnalyzer") << " ECAL Endcap RecHit Collection not available !"; return;
  } else {
    recHitsEE = recHitsEE_h.product();
  }
  //TO DISENTANGLE BETWEEN MINIAOD AND AOD

  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  edm::ESHandle<EcalChannelStatus> chStatus;
  iSetup.get<EcalChannelStatusRcd>().get(chStatus);
  const EcalChannelStatus *ch_status = chStatus.product(); 

  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  //get the reference to 1st vertex for use in fGetIsolation
  //for PFIsolation calculation
  reco::VertexRef firstVtx(vertexColl,0);
   
  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  // get the photon collection
  // Handle<reco::PhotonCollection> photonColl;
  // iEvent.getByLabel(fPhotonTag,photonColl);
  // edm::Handle<pat::PhotonCollection> photonColl;
  // iEvent.getByLabel("slimmedPhotons",photonColl);

  // If photon collection is empty, exit
  // if (!photonColl.isValid()) {
  //   cout << "No Photons! Move along, there's nothing to see here .." <<endl;
  //   return;
  // }
  // const reco::PhotonCollection *myPhotonColl = photonColl.product();
  // cout<<"photoncoll size "<<myPhotonColl->size()<<endl;
  //   cout << "N reco photons = " << photonColl->size() <<endl;

  TString CategoryPFID(fPFIDCategory.c_str());
  TString MethodID(fIDMethod.c_str());

  // new approach - 
  // make vector of all selected Photons (tight ID, not spike, min pt, etc)
  // then sort at end by pt
  // also allows to count how often a third photon could be considered a candidate

  std::vector<reco::Photon> selectedPhotons; 

  // do a similar thing for 'Fakeable objects', for data-based fake rate approach
  std::vector<reco::Photon> fakeablePhotons; 

  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  int phoIndex = -1;

  std::map<float,bool> vetoMap;

  // if using MiniAOD, build the vetoMap using pat::Photons (ConversionTools::hasMatchedPromptElectron was not working for pat::Photons even though they inherit from reco::Photons)

  if (!fisAOD){
    edm::Handle<edm::View<pat::Photon> > patPhotons;
    iEvent.getByToken(patPhotonToken_, patPhotons);
    for (auto &p : *patPhotons){
      std::pair<float,bool> tempPair = std::make_pair<float,bool>(p.pt(),p.passElectronVeto());
      vetoMap.insert(tempPair);
    }
  }

  else{
    for(edm::View<reco::Photon>::const_iterator recoPhoton = photons->begin(); recoPhoton!=photons->end(); recoPhoton++) {
      // bool passelecveto = !ConversionTools::hasMatchedPromptElectron(recoPhoton->superCluster(), hElectrons, hConversions, beamSpot.position());
      std::pair<float,bool> tempPair = std::make_pair<float,bool>(recoPhoton->pt(),!ConversionTools::hasMatchedPromptElectron(recoPhoton->superCluster(), hElectrons, hConversions, beamSpot.position()));
      vetoMap.insert(tempPair);
    }

  }

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
  fCand_CHpos_dxy.clear();
  fCand_CHpos_dxy_beamspot.clear();
  fCand_CHpos_dxy_associated.clear();
  fCand_CHneg_dxy.clear();
  fCand_CHneg_dxy_beamspot.clear();
  fCand_CHneg_dxy_associated.clear();
  fCand_pass.clear();
  fCand_passChargedIso.clear();
  fCand_passNeutralIso.clear();
  fCand_passEGammaIso.clear();
  fCand_passPhotonPt.clear();
  fCand_fake.clear();
  fCand_match.clear();
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
  fTwoProng_CHpos_dxy.clear();
  fTwoProng_CHpos_dxy_beamspot.clear();
  fTwoProng_CHpos_dxy_associated.clear();
  fTwoProng_CHneg_dxy.clear();
  fTwoProng_CHneg_dxy_beamspot.clear();
  fTwoProng_CHneg_dxy_associated.clear();
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
  fGenEta_px.clear();
  fGenEta_py.clear();
  fGenEta_pz.clear();
  fGenEta_energy.clear();

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
  TLorentzVector LeadingFake; 
  if (fDebug) cout << ". starting candidate loop" << endl;
  for (unsigned int i = 0; i < pfcands->size(); i++) {
    const pat::PackedCandidate &pf1 = (*pfcands)[i];
    if (pf1.pt() < fCandidatePairMinPt) continue;
    pruned_count += 1;
    for (unsigned int j = i+1; j < pfcands->size(); j++) { // note loop starting with j=i+1, considers each pair exactly once
      const pat::PackedCandidate &pf2 = (*pfcands)[j];
      if (pf2.pt() < fCandidatePairMinPt) continue;
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
          // CH pair is close and has at least one pf photon, definition of candidate, fill vectors
          leading_pf_photon.SetPtEtaPhiE((*pfcands)[n].pt(), (*pfcands)[n].eta(), (*pfcands)[n].phiAtVtx(), (*pfcands)[n].energy());
          TLorentzVector EtaCandidate;
          EtaCandidate = center + photon;
          TLorentzVector EtaMassCandidate;
          EtaMassCandidate = center + leading_pf_photon;
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
            fCand_CHpos_dxy.push_back(pf1.dxy(((*primaryvertecies)[0]).position()));
            fCand_CHpos_dxy_beamspot.push_back(pf1.dxy(beamspot->position()));
            fCand_CHpos_dxy_associated.push_back(pf1.dxy());
            fCand_CHneg_dxy.push_back(pf2.dxy(((*primaryvertecies)[0]).position()));
            fCand_CHneg_dxy_beamspot.push_back(pf2.dxy(beamspot->position()));
            fCand_CHneg_dxy_associated.push_back(pf2.dxy());
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
            fCand_CHpos_dxy.push_back(pf2.dxy(((*primaryvertecies)[0]).position()));
            fCand_CHpos_dxy_beamspot.push_back(pf2.dxy(beamspot->position()));
            fCand_CHpos_dxy_associated.push_back(pf2.dxy());
            fCand_CHneg_dxy.push_back(pf1.dxy(((*primaryvertecies)[0]).position()));
            fCand_CHneg_dxy_beamspot.push_back(pf1.dxy(beamspot->position()));
            fCand_CHneg_dxy_associated.push_back(pf1.dxy());
          }
          fCand_photon_pt.push_back(photon.Pt());
          fCand_photon_eta.push_back(photon.Eta());
          fCand_photon_phi.push_back(photon.Phi());
          fCand_photon_mass.push_back(photon.M());
          fCand_photon_nGamma.push_back(numgamma);
          fCand_photon_nElectron.push_back(nume);
          fCand_pt.push_back(EtaCandidate.Pt());
          fCand_eta.push_back(EtaCandidate.Eta());
          fCand_phi.push_back(EtaCandidate.Phi());
          fCand_mass.push_back(EtaCandidate.M());
          fCand_Mass.push_back(EtaMassCandidate.M());
          // Now define isolations
          double chargedIso = 0;
          double neutralIso = 0;
          double egammaIso = 0;
          for (unsigned int m = 0; m < pfcands->size(); m++) {
            const pat::PackedCandidate &pf4 = (*pfcands)[m];
            TLorentzVector pfcand4;
            pfcand4.SetPtEtaPhiE(pf4.pt(), pf4.eta(), pf4.phiAtVtx(), pf4.energy());
            if (!fisSignal && pf4.fromPV() <= 1) continue;
            // charged (incl. muons)
            if (abs(pf4.pdgId()) == 13 || abs(pf4.pdgId()) == 211) {
              if ( center.DeltaR(pfcand4) < fCandidatePairIsolationDR && !(m == i || m == j) ) // iso cand can't be one of CH from CH pair
                  chargedIso += pfcand4.Pt();
            // neutral
            } else if (pf4.pdgId() == 130) {
              if (center.DeltaR(pfcand4) < fCandidatePairIsolationDR)
                neutralIso += pfcand4.Pt();
            // e gamma
            } else if (abs(pf4.pdgId()) == 11 || pf4.pdgId() == 22) {
              if ( (center.DeltaR(pfcand4) < fCandidatePairIsolationDR) &&
                   !(fabs(pf4.phiAtVtx() - center.Phi()) < fCandidatePairPhiBox/2.0 && fabs(pf4.eta() - center.Eta()) < fCandidatePairEtaBox/2.0))
                egammaIso += pfcand4.Pt();
            }
          } // end pf cand loop
          double relchargedIso = chargedIso / EtaCandidate.Pt();
          double relneutralIso = neutralIso / EtaCandidate.Pt();
          double relegammaIso = egammaIso / EtaCandidate.Pt();
          if (fDebug) cout << ". finished isolation" << endl;
          fCand_chargedIso.push_back(relchargedIso);
          fCand_neutralIso.push_back(relneutralIso);
          fCand_egammaIso.push_back(relegammaIso);
          
          // Selection on Candidates
          bool passCharged = relchargedIso < fCandidatePairChargedIsoCut;
          bool passNeutral = relneutralIso < fCandidatePairNeutralIsoCut;
          bool passEGamma = relegammaIso < fCandidatePairEGammaIsoCut;
          bool passPhotonPt = photon.Pt() > fCandidatePairPhotonPtCut;
          bool pass = passCharged && passNeutral && passEGamma && passPhotonPt;
          bool fake = !passCharged && passNeutral && passEGamma && passPhotonPt;
          fCand_pass.push_back(pass);
          fCand_passChargedIso.push_back(passCharged);
          fCand_passNeutralIso.push_back(passNeutral);
          fCand_passEGammaIso.push_back(passEGamma);
          fCand_passPhotonPt.push_back(passPhotonPt);
          fCand_fake.push_back(fake);
          // Generator Matching
          bool match = false;
	        if (fisSignal && fisMC) {
	          for (unsigned int i = 0; i < genparticles->size(); i++) {
	            const reco::GenParticle &genparticle = (*genparticles)[i];
	            if ((genparticle.pdgId() == 221 || genparticle.pdgId() == 331) && genparticle.status() == 2) {
                TLorentzVector genEta;
                genEta.SetPtEtaPhiM(genparticle.pt(), genparticle.eta(), genparticle.phi(), genparticle.mass());
                double match_dR = genEta.DeltaR(EtaCandidate);
                if (match_dR < fCandidatePairGenMatchDR)
                  match = true;
              }
            }
          }
          fCand_match.push_back(match);
          // Cutflow variables
          if (passCharged) nPassCharged++;
          if (passNeutral) nPassNeutral++;
          if (passEGamma) nPassEGamma++;
          if (passPhotonPt) nPassPhotonPt++;
          if (match) nMatch++;
          if (match && pass) nPassMatch++;
          // Fake rate histograms
          if (pass) {
            fTwoProngFakeNumer_pt->Fill(EtaCandidate.Pt());
            fTwoProngFakeNumer_eta->Fill(EtaCandidate.Eta());
            fTwoProngFakeNumer_phi->Fill(EtaCandidate.Phi());
          } if (fake) {
            fTwoProngFakeDenom_pt->Fill(EtaCandidate.Pt());
            fTwoProngFakeDenom_eta->Fill(EtaCandidate.Eta());
            fTwoProngFakeDenom_phi->Fill(EtaCandidate.Phi());
            nFake++;
            if (EtaCandidate.Pt() > LeadingFake.Pt()) LeadingFake = EtaCandidate;
          }
          if (fDebug) cout << ". finished fake rate filling" << endl;
        }
      } // end conditionals on CH pair
    }
  } // end making candidates
  fNumPrunedPF = pruned_count;
  fNumCHpairsFake = nFake;

  // Create sorted-by-pt list of passed candidates
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
    if (fCand_pass[index])
    {
      // Candidate passes and is next leading, fill all passed candidate collections
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
      fTwoProng_CHneg_pt.push_back(fCand_CHneg_pt[index]);
      fTwoProng_CHneg_eta.push_back(fCand_CHneg_eta[index]);
      fTwoProng_CHneg_phi.push_back(fCand_CHneg_phi[index]);
      fTwoProng_CHneg_mass.push_back(fCand_CHneg_mass[index]);
      fTwoProng_photon_pt.push_back(fCand_photon_pt[index]);
      fTwoProng_photon_eta.push_back(fCand_photon_eta[index]);
      fTwoProng_photon_phi.push_back(fCand_photon_phi[index]);
      fTwoProng_photon_mass.push_back(fCand_photon_mass[index]);
      fTwoProng_photon_nGamma.push_back(fCand_photon_nGamma[index]);
      fTwoProng_photon_nElectron.push_back(fCand_photon_nElectron[index]);
      fTwoProng_chargedIso.push_back(fCand_chargedIso[index]);
      fTwoProng_neutralIso.push_back(fCand_neutralIso[index]);
      fTwoProng_egammaIso.push_back(fCand_egammaIso[index]);
      fTwoProng_CHpos_dxy_beamspot.push_back(fCand_CHpos_dxy_beamspot[index]);
      fTwoProng_CHpos_dxy_associated.push_back(fCand_CHpos_dxy_associated[index]);
      fTwoProng_CHpos_dxy.push_back(fCand_CHpos_dxy[index]);
      fTwoProng_CHneg_dxy_beamspot.push_back(fCand_CHneg_dxy_beamspot[index]);
      fTwoProng_CHneg_dxy_associated.push_back(fCand_CHneg_dxy_associated[index]);
      fTwoProng_CHneg_dxy.push_back(fCand_CHneg_dxy[index]);
      fTwoProng_match.push_back(fCand_match[index]);
    }
  }
  if (fDebug) cout << ". finished passed collections" << endl;

  // Fill other event wide information
  fNumCHpairs = fCand_pt.size();
  fNumCHpairsPass = fTwoProng_pt.size();
  fNumCHpairsPassChargedIso = nPassCharged;
  fNumCHpairsPassNeutralIso = nPassNeutral;
  fNumCHpairsPassEGammaIso = nPassEGamma;
  fNumCHpairsPassChargedIso = nPassPhotonPt;
  fHT = 0;
  for (unsigned int i = 0; i < ak4jets->size(); i++) {
    const pat::Jet &jet = (*ak4jets)[i];
    if (jet.pt() < 30) continue;
    if (fabs(jet.eta()) > 2.5) continue;
    fHT += jet.pt();
  }
  
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

  // photon loop
  for(edm::View<reco::Photon>::const_iterator recoPhoton = photons->begin(); recoPhoton!=photons->end(); recoPhoton++) {
    
    phoIndex++;

    //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

    const reco::Photon testPho = *recoPhoton;
    edm::Ptr<reco::Photon> testPhoPtr(photons,phoIndex);

    //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

    //-----------------taken from Ilya-----------------
    float full5x5sigmaIetaIeta = (*full5x5SigmaIEtaIEtaMap)[testPhoPtr];
    float chIso =  (*phoChargedIsolationMap)[testPhoPtr];
    float nhIso =  (*phoNeutralHadronIsolationMap)[testPhoPtr];
    float phIso = (*phoPhotonIsolationMap)[testPhoPtr];

    // for MiniAOD, introduce chIso cut to avoid the case where a supercluster doesn't have rechit/cluster information that is counted on later
    // cout << "SK" << recoPhoton->r9() << " " << chIso << " " << recoPhoton->chargedHadronIso() << " " << 0.3*recoPhoton->pt() << endl;
    if (recoPhoton->chargedHadronIso() > 10.){
      if (fDebug) cout << "Photon with chargedHadronIso="<<recoPhoton->chargedHadronIso() << " pt=" << recoPhoton->pt() << "skipped!" << endl;
      continue;
    }
    
    float abseta = fabs( recoPhoton->superCluster()->eta());

    //we set the effective areas
    //based on whether we use egamma ID
    //or high pt photon ID
    float CHeffarea = 0.;
    float NHeffarea = 0.;
    float PHeffarea = 0.;

    //let's retrieve the effective areas
    //from exodiphotons (in case of high pt ID)
    //Remember effareaCH = 1st, effareaNH = 2nd, effareaPH = 3rd
    std::vector<double> effareas = ExoDiPhotons::EffectiveAreas(&(*recoPhoton),MethodID,CategoryPFID);

    if(MethodID.Contains("highpt")){
      CHeffarea = effareas[0];
      NHeffarea = effareas[1];
      PHeffarea = effareas[2];
    }

    if(MethodID.Contains("egamma")){
      CHeffarea = effAreaChHadrons_.getEffectiveArea(abseta);
      NHeffarea = effAreaNeuHadrons_.getEffectiveArea(abseta);
      PHeffarea = effAreaPhotons_.getEffectiveArea(abseta);
    }

    //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

    float isoChargedHadronsWithEA = std::max( (float)0.0, chIso - rho_*CHeffarea);
    float isoNeutralHadronsWithEA = std::max( (float)0.0, nhIso - rho_*NHeffarea);
    float isoPhotonsWithEA = std::max( (float)0.0, phIso - rho_*PHeffarea);
    
    //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

    bool isPassLoose  = (*loose_id_decisions)[testPhoPtr];
    bool isPassMedium = (*medium_id_decisions)[testPhoPtr];
    bool isPassTight  = (*tight_id_decisions)[testPhoPtr];
    
    if (fDebug) {
      cout<<full5x5sigmaIetaIeta<<" "
        <<chIso<<" "
        <<nhIso<<" "
        <<phIso<<" "
        <<isPassLoose<<" "
        <<isPassMedium<<" "
        <<isPassTight<<" "
        <<isoChargedHadronsWithEA<<" "
        <<isoNeutralHadronsWithEA<<" "
        <<isoPhotonsWithEA<<" "
        <<endl;
    }

    //-----------------taken from Ilya-----------------

    // now add selected photons to vector if:
    // tight ID
    // not in gap
    // not a spike
    // min pt cut  (same for all photons)
    // in EB only (use detector eta)

    // we will check both tight and fakeable photons
    // some cuts are common to both - EB and pt
    //apr 2011, remove EB cut - what should we do about the increased combinatorics?



    //Unfortunately, we have to compute PF isolation variables here
    //to check if the photon is tight or fakeable
    //But then, we have to recompute them later for each tight or fakeable photon
    //because once again the PF isolation variables are not accessible 
    //in the Photon class
    //We need to compute only 03 isol variables so far
    //because they are the official ones

    double rhocorPFIsoCH = isoChargedHadronsWithEA;
    double rhocorPFIsoNH = isoNeutralHadronsWithEA;
    double rhocorPFIsoPH = isoPhotonsWithEA;
    //double pfisoall = rhocorPFIsoCH + rhocorPFIsoNH + rhocorPFIsoPH;
    
    // conversion safe electron veto
    bool passelecveto = vetoMap.at(recoPhoton->pt());
    // bool passelecveto = !ConversionTools::hasMatchedPromptElectron(recoPhoton->superCluster(), hElectrons, hConversions, beamSpot.position());

    // if (passelecveto) cout << "Passed electron veto!" << endl;
    // else cout << "Failed electron veto!" << endl;
    //bool passelecveto = true;
    //TO DISENTANGLE BETWEEN MINIAOD AND AOD

    // cout << "Photon et, eta, phi = " << recoPhoton->et() <<", "<<recoPhoton->eta()<< ", "<< recoPhoton->phi();
    // cout << "; calo position eta = " << recoPhoton->caloPosition().eta();
    // cout << "; eMax/e3x3 = " << recoPhoton->maxEnergyXtal()/recoPhoton->e3x3();
    // cout << "; hadOverEm = " << recoPhoton->hadronicOverEm();
    // cout << "; pixelSeed = " << recoPhoton->hasPixelSeed();
    // cout << "; sigmaietaieta = " << recoPhoton->sigmaIetaIeta();
    // cout << "; CHiso = " << rhocorPFIsoCH;
    // cout << "; NHiso = " << rhocorPFIsoNH;
    // cout << "; PHiso = " << rhocorPFIsoPH;
    // cout<<""<<endl;

    // cout<<full5x5sigmaIetaIeta<<" "
    // 	<<chIso<<" "
    // 	<<nhIso<<" "
    // 	<<phIso<<" "
    // 	<<isPassLoose<<" "
    // 	<<isPassMedium<<" "
    // 	<<isPassTight<<" "
    // 	<<isoChargedHadronsWithEA<<" "
    // 	<<isoNeutralHadronsWithEA<<" "
    // 	<<isoPhotonsWithEA<<" "
    // 	<<endl;

    //trick to use existing saturation code
    // reco::SuperClusterRef sc1 = photon->superCluster();
    // if (*sc1).size() == 0{
      
    // }
    ExoDiPhotons::recoPhotonInfo_t tempInfo;
    ExoDiPhotons::FillRecoPhotonInfo(tempInfo,&(*recoPhoton),lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);
    Bool_t isSaturated = tempInfo.isSaturated;

    if (fDebug) {
      if(recoPhoton->pt() < fMin_pt) cout << "photon pt = " << recoPhoton->pt() << " which is less than " << fMin_pt << "!" << endl;
      if(ExoDiPhotons::isGapPhoton(&(*recoPhoton))) cout << "this is a photon in the gap!" << endl;
      if(ExoDiPhotons::isASpike(&(*recoPhoton))) cout << "this photon is a spike!" << endl;
    }

    if(recoPhoton->pt() < fMin_pt) continue;
    if(ExoDiPhotons::isGapPhoton(&(*recoPhoton))) continue;
    if(ExoDiPhotons::isASpike(&(*recoPhoton))) continue;

    if (fDebug) cout << "pt, gap, and spike cuts passed!" << endl;
    
    //Now we choose which ID to use (PF or Det)
    if(MethodID.Contains("highpt")){
      // if(ExoDiPhotons::isPFTightPhoton(&(*recoPhoton),rhocorPFIsoCH,rhocorPFIsoNH,rhocorPFIsoPH,full5x5sigmaIetaIeta,passelecveto,MethodID,CategoryPFID,isSaturated)){
      if(ExoDiPhotons::passHighPtID(&(*recoPhoton),MethodID,CategoryPFID,rhocorPFIsoCH,phIso,full5x5sigmaIetaIeta,rho_,passelecveto,isSaturated)){
        selectedPhotons.push_back(*recoPhoton);
        if (fDebug) cout << "photon passed the high pt id!" << endl;
      }
      else {
        if (fDebug) cout << "photon failed the high pt id!" << endl;
      }
    }
    else if(MethodID.Contains("egamma")){
      if(isPassLoose) selectedPhotons.push_back(*recoPhoton);
    }

    // also check for fakeable objects
    if(MethodID.Contains("highpt")){
      if(ExoDiPhotons::isPFFakeableObject(&(*recoPhoton),rhocorPFIsoCH,rhocorPFIsoNH,rhocorPFIsoPH,full5x5sigmaIetaIeta,passelecveto,MethodID,CategoryPFID,isSaturated) ) {

	//        cout << "Fakeable photon! ";
	//        cout << "Photon et, eta, phi = " << recoPhoton->et() <<", "<<recoPhoton->eta()<< ", "<< recoPhoton->phi();
	//      //     cout << "; calo position eta = " << recoPhoton->caloPosition().eta();
	// //      cout << "; eMax/e3x3 = " << recoPhoton->maxEnergyXtal()/recoPhoton->e3x3();
	//        cout << "; hadOverEm = " << recoPhoton->hadronicOverEm();
	//        cout << "; trkIso = " << recoPhoton->trkSumPtHollowConeDR04();
	//        cout << "; ecalIso = " << recoPhoton->ecalRecHitSumEtConeDR04();
	//        cout << "; hcalIso = " << recoPhoton->hcalTowerSumEtConeDR04();
	//        //      cout << "; pixelSeed = " << recoPhoton->hasPixelSeed();
	//        cout << "; sigmaietaieta = " << recoPhoton->sigmaIetaIeta();
	//        cout << endl;
	 
	 
	fakeablePhotons.push_back(*recoPhoton);
  if (fDebug) cout << "photon pased fakeable object cut" << endl;
      }
    }
    else if(MethodID.Contains("egamma")){
      if(ExoDiPhotons::isPFFakeableObject(&(*recoPhoton),rhocorPFIsoCH,rhocorPFIsoNH,rhocorPFIsoPH,full5x5sigmaIetaIeta,passelecveto,MethodID,CategoryPFID,isSaturated) ) {
	 
	//        cout << "Fakeable photon! ";
	//        cout << "Photon et, eta, phi = " << recoPhoton->et() <<", "<<recoPhoton->eta()<< ", "<< recoPhoton->phi();
	//      //     cout << "; calo position eta = " << recoPhoton->caloPosition().eta();
	// //      cout << "; eMax/e3x3 = " << recoPhoton->maxEnergyXtal()/recoPhoton->e3x3();
	//        cout << "; hadOverEm = " << recoPhoton->hadronicOverEm();
	//        cout << "; trkIso = " << recoPhoton->trkSumPtHollowConeDR04();
	//        cout << "; ecalIso = " << recoPhoton->ecalRecHitSumEtConeDR04();
	//        cout << "; hcalIso = " << recoPhoton->hcalTowerSumEtConeDR04();
	//        //      cout << "; pixelSeed = " << recoPhoton->hasPixelSeed();
	//        cout << "; sigmaietaieta = " << recoPhoton->sigmaIetaIeta();
	//        cout << endl;
	 	 
	fakeablePhotons.push_back(*recoPhoton);
      }
    }
       
  } //end reco photon loop


  // now sort the vector of selected photons by pt
  // (compare function is found in RecoPhotonInfo.h)
  sort(selectedPhotons.begin(),selectedPhotons.end(),ExoDiPhotons::comparePhotonsByPt);
  // check sorting
  //   cout << "After sorting photons" <<endl;
  //   for(std::vector<reco::Photon>::iterator myPhotonIter = selectedPhotons.begin();myPhotonIter<selectedPhotons.end();myPhotonIter++) {
  //     cout << "Photon et, eta, phi = " << myPhotonIter->et() <<", "<<myPhotonIter->eta()<< ", "<< myPhotonIter->phi()<<endl;
  //   }


  ////////////////////////////////////////////////////////////////////////////////
  //
  // Add conversion information for all tight photons (if there are 2)
  //
  ///////////////////////////////////////////////////////////////////////////////
  if (fDebug) cout << "got to conversions" << endl;
  if (selectedPhotons.size() >= 2){
    edm::Handle<reco::ConversionCollection> twoLegHandle;
    edm::Handle<reco::ConversionCollection> oneLegHandle;

    iEvent.getByToken(twoLegToken_,twoLegHandle);
    iEvent.getByToken(oneLegToken_,oneLegHandle);

    reco::ConversionCollection twoLegConversions = *(twoLegHandle.product());
    reco::ConversionCollection oneLegConversions = *(oneLegHandle.product());

    reco::Photon iPho1 = selectedPhotons.at(0);
    reco::Photon iPho2 = selectedPhotons.at(1);
    reco::SuperCluster sc1 = *(iPho1.superCluster());
    reco::SuperCluster sc2 = *(iPho2.superCluster());
    ExoDiPhotons::FillConversionInfo(fConvInfo1,sc1,twoLegConversions,iPho1.pt(), beamSpot);
    ExoDiPhotons::FillConversionInfo(fConvInfo2,sc2,twoLegConversions,iPho2.pt(), beamSpot);
    ExoDiPhotons::FillConversionInfo(fConvInfo_OneLeg1,sc1,oneLegConversions,iPho1.pt(), beamSpot);
    ExoDiPhotons::FillConversionInfo(fConvInfo_OneLeg2,sc2,oneLegConversions,iPho2.pt(), beamSpot);

  } // end conversion info block

  // now count many candidate photons we have in this event
  //   cout << "N candidate photons = " << selectedPhotons.size() <<endl;
  fNTightPhotons = selectedPhotons.size();
  if (fDebug) {
    if(fNTightPhotons >= 2) cout<<"great we have two tight photons"<<endl;
    else cout << "only " << fNTightPhotons << "photons passed the cuts!" << endl;
  }
  // now consider possible Fakeable objects

  // first sort by pt
  sort(fakeablePhotons.begin(),fakeablePhotons.end(),ExoDiPhotons::comparePhotonsByPt);

  //   cout << "N fakeable = " << fakeablePhotons.size() <<endl;
  fNFakeablePhotons = fakeablePhotons.size();

  // our new fake handling:
  // make decision on how to handle the event
  // based on the two highest pt objects, tight or fakeable

  // to do this, best to make a list of all the objects, tight or Fakeable
  // then sort this by photon pt
  // we pair the photon with a boolean to mark tight or fakeable status
  // our convention is that the output tree has an isFakeable leaf
  // hence tight=false and fakeable = true
  std::vector<std::pair<reco::Photon, bool> > allTightOrFakeableObjects;
   
  // actually, I realise that now we are using this vector of pairs,
  // we probably dont need anymore the original tight/fakeable separate vectors
  // because we could probably have filled this vector of pairs directly
  // but rather than change it now, we leave it the way it is 

  for(std::vector<reco::Photon>::iterator myPhotonIter = selectedPhotons.begin();myPhotonIter<selectedPhotons.end();myPhotonIter++) {
    //     cout << "Photon et, eta, phi = " << myPhotonIter->et() <<", "<<myPhotonIter->eta()<< ", "<< myPhotonIter->phi()<<endl;
    std::pair<reco::Photon, bool> myPair(*myPhotonIter,false); 
    // tight, so isFakeble=false

    allTightOrFakeableObjects.push_back(myPair);
  }
   
  for(std::vector<reco::Photon>::iterator myPhotonIter = fakeablePhotons.begin();myPhotonIter<fakeablePhotons.end();myPhotonIter++) {
    //     cout << "Photon et, eta, phi = " << myPhotonIter->et() <<", "<<myPhotonIter->eta()<< ", "<< myPhotonIter->phi()<<endl;
    std::pair<reco::Photon, bool> myPair(*myPhotonIter,true);
    // isFakeable = true
    allTightOrFakeableObjects.push_back(myPair);
  }
   
  // now sort according to photon pt
  sort(allTightOrFakeableObjects.begin(),allTightOrFakeableObjects.end(),ExoDiPhotons::comparePhotonPairsByPt);
  if (fDebug) cout << "allTightOrFakeableObjects.size(): " << allTightOrFakeableObjects.size() << endl;
  if(allTightOrFakeableObjects.size()>=2) {

    // now, we are always going to consider the top two objects
    // regardless of their 'nature'
    // so we can fill the structs now

    // order Photon1,2 by pt, no matter which is tight and which is fake
    // and then use isFakeable leaf to distinguish them

    // must specifically declare isFakeable status
    // the bool in hte pair uses the same convention, so:

    //since we want to have the PF isolation variables
    //we need edm::Ptr to the reco photons
    //and for that we need the index in the photon collection
    //but we lost track of it when sorting our vectors
    //so lets find it back

    if (fDebug) cout<<"here testing whether i get the correct pointer to my photon TorF"<<endl;
    reco::Photon* TorFObject1 = &(allTightOrFakeableObjects[0].first);
    reco::Photon* TorFObject2 = &(allTightOrFakeableObjects[1].first);
    if (fDebug) {
      cout<<TorFObject1->energy()<<" "
    <<TorFObject1->eta()<<" "
    <<TorFObject1->et()<<" "
    <<TorFObject1->phi()<<" "
    <<endl;
      cout<<TorFObject2->energy()<<" "
    <<TorFObject2->eta()<<" "
    <<TorFObject2->et()<<" "
    <<TorFObject2->phi()<<" "
    <<endl;
    }

    int indexTorFObject1 = -1;
    int indexTorFObject2 = -1;
    int myIndex = -1;

    for(edm::View<reco::Photon>::const_iterator recoPhoton = photons->begin(); recoPhoton!=photons->end(); recoPhoton++) {
      myIndex++;
      //const reco::Photon* myPhoton = &(*recoPhoton);
      //if(myPhoton == TorFObject1) {
      if(recoPhoton->pt() == TorFObject1->pt()) {
  if (fDebug)	cout<<"Great, I've found torf object 1 "<<endl;
	indexTorFObject1 = myIndex;
      }
      //if(myPhoton == TorFObject2) {
      if(recoPhoton->pt() == TorFObject2->pt()) {
  if (fDebug)	cout<<"Great, I've found torf object 2 "<<endl;
	indexTorFObject2 = myIndex;
      }
    }

    //in principle the indices should be always greater than 1
    //because we necessarily found the two objects in the 
    //photon collection


    edm::Ptr<reco::Photon> TorFObject1Ptr(photons,indexTorFObject1);
    if (fDebug) {
      cout<<TorFObject1Ptr->energy()<<" "
    <<TorFObject1Ptr->eta()<<" "
    <<TorFObject1Ptr->et()<<" "
    <<TorFObject1Ptr->phi()<<" "
    <<endl;
    }

    edm::Ptr<reco::Photon> TorFObject2Ptr(photons,indexTorFObject2);
    if (fDebug) {
      cout<<TorFObject2Ptr->energy()<<" "
    <<TorFObject2Ptr->eta()<<" "
    <<TorFObject2Ptr->et()<<" "
    <<TorFObject2Ptr->phi()<<" "
    <<endl;
    cout<<"here testing whether i get the correct pointer to my photon TorF"<<endl;
    }

    ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo1,&allTightOrFakeableObjects[0].first,lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);
    fRecoPhotonInfo1.isFakeable = allTightOrFakeableObjects[0].second;
    // fRecoPhotonInfo1.hasMatchedPromptElec = ConversionTools::hasMatchedPromptElectron((&allTightOrFakeableObjects[0].first)->superCluster(), hElectrons, hConversions, beamSpot.position());
    fRecoPhotonInfo1.hasMatchedPromptElec = !( vetoMap.at((&allTightOrFakeableObjects[0].first)->pt()) );
    
    fRecoPhotonInfo1.isTightPFPhoton = (*tight_id_decisions)[TorFObject1Ptr];
    fRecoPhotonInfo1.isMediumPFPhoton = (*medium_id_decisions)[TorFObject1Ptr];
    fRecoPhotonInfo1.isLoosePFPhoton = (*loose_id_decisions)[TorFObject1Ptr]; 

    //TO DISENTANGLE BETWEEN MINIAOD AND AOD

    fRecoPhotonInfo1.sigmaIetaIeta = (*full5x5SigmaIEtaIEtaMap)[TorFObject1Ptr];
    fRecoPhotonInfo1.PFIsoCharged03 = (*phoChargedIsolationMap)[TorFObject1Ptr];
    fRecoPhotonInfo1.PFIsoNeutral03 = (*phoNeutralHadronIsolationMap)[TorFObject1Ptr];
    fRecoPhotonInfo1.PFIsoPhoton03 = (*phoPhotonIsolationMap)[TorFObject1Ptr];
    fRecoPhotonInfo1.PFIsoAll03 = fRecoPhotonInfo1.PFIsoCharged03 + fRecoPhotonInfo1.PFIsoNeutral03 + fRecoPhotonInfo1.PFIsoPhoton03;

    float torfObject1Eta = abs(TorFObject1Ptr->superCluster()->eta());

    float CHeffarea = 0.;
    float NHeffarea = 0.;
    float PHeffarea = 0.;

    std::vector<double> effareas = ExoDiPhotons::EffectiveAreas((&allTightOrFakeableObjects[0].first),MethodID,CategoryPFID);
    if(MethodID.Contains("highpt")){
      CHeffarea = effareas[0];
      NHeffarea = effareas[1];
      PHeffarea = effareas[2];
    }
    if(MethodID.Contains("egamma")){
      CHeffarea = effAreaChHadrons_.getEffectiveArea(torfObject1Eta);
      NHeffarea = effAreaNeuHadrons_.getEffectiveArea(torfObject1Eta);
      PHeffarea = effAreaPhotons_.getEffectiveArea(torfObject1Eta);
    }

    //now the corrected PF isolation variables
    fRecoPhotonInfo1.rhocorPFIsoCharged03 = std::max((float)0.0,(float)fRecoPhotonInfo1.PFIsoCharged03-rho_*CHeffarea);
    fRecoPhotonInfo1.rhocorPFIsoNeutral03 = std::max((float)0.0,(float)fRecoPhotonInfo1.PFIsoNeutral03-rho_*NHeffarea);
    fRecoPhotonInfo1.rhocorPFIsoPhoton03 = std::max((float)0.0,(float)fRecoPhotonInfo1.PFIsoPhoton03-rho_*PHeffarea);
    fRecoPhotonInfo1.rhocorPFIsoAll03 = fRecoPhotonInfo1.rhocorPFIsoCharged03 + fRecoPhotonInfo1.rhocorPFIsoNeutral03 + fRecoPhotonInfo1.rhocorPFIsoPhoton03
      ;

    ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo2,&allTightOrFakeableObjects[1].first,lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);
    fRecoPhotonInfo2.isFakeable = allTightOrFakeableObjects[1].second;
    // fRecoPhotonInfo2.hasMatchedPromptElec = ConversionTools::hasMatchedPromptElectron((&allTightOrFakeableObjects[1].first)->superCluster(), hElectrons, hConversions, beamSpot.position());
    fRecoPhotonInfo2.hasMatchedPromptElec = !( vetoMap.at((&allTightOrFakeableObjects[1].first)->pt()) );

    fRecoPhotonInfo2.isTightPFPhoton = (*tight_id_decisions)[TorFObject2Ptr];
    fRecoPhotonInfo2.isMediumPFPhoton = (*medium_id_decisions)[TorFObject2Ptr];
    fRecoPhotonInfo2.isLoosePFPhoton = (*loose_id_decisions)[TorFObject2Ptr]; 

    //TO DISENTANGLE BETWEEN MINIAOD AND AOD


    fRecoPhotonInfo2.sigmaIetaIeta = (*full5x5SigmaIEtaIEtaMap)[TorFObject2Ptr];
    fRecoPhotonInfo2.PFIsoCharged03 = (*phoChargedIsolationMap)[TorFObject2Ptr];
    fRecoPhotonInfo2.PFIsoNeutral03 = (*phoNeutralHadronIsolationMap)[TorFObject2Ptr];
    fRecoPhotonInfo2.PFIsoPhoton03 = (*phoPhotonIsolationMap)[TorFObject2Ptr];
    fRecoPhotonInfo2.PFIsoAll03 = fRecoPhotonInfo2.PFIsoCharged03 + fRecoPhotonInfo2.PFIsoNeutral03 + fRecoPhotonInfo2.PFIsoPhoton03;

    float torfObject2Eta = abs(TorFObject2Ptr->superCluster()->eta());

    CHeffarea = 0.;
    NHeffarea = 0.;
    PHeffarea = 0.;

    effareas.clear();
    effareas = ExoDiPhotons::EffectiveAreas((&allTightOrFakeableObjects[1].first),MethodID,CategoryPFID);
    if(MethodID.Contains("highpt")){
      CHeffarea = effareas[0];
      NHeffarea = effareas[1];
      PHeffarea = effareas[2];
    }
    if(MethodID.Contains("egamma")){
      CHeffarea = effAreaChHadrons_.getEffectiveArea(torfObject2Eta);
      NHeffarea = effAreaNeuHadrons_.getEffectiveArea(torfObject2Eta);
      PHeffarea = effAreaPhotons_.getEffectiveArea(torfObject2Eta);
    }

    //now the corrected PF isolation variables
    fRecoPhotonInfo2.rhocorPFIsoCharged03 = std::max((float)0.0,(float)fRecoPhotonInfo2.PFIsoCharged03-rho_*CHeffarea);
    fRecoPhotonInfo2.rhocorPFIsoNeutral03 = std::max((float)0.0,(float)fRecoPhotonInfo2.PFIsoNeutral03-rho_*NHeffarea);
    fRecoPhotonInfo2.rhocorPFIsoPhoton03 = std::max((float)0.0,(float)fRecoPhotonInfo2.PFIsoPhoton03-rho_*PHeffarea);
    fRecoPhotonInfo2.rhocorPFIsoAll03 = fRecoPhotonInfo2.rhocorPFIsoCharged03 + fRecoPhotonInfo2.rhocorPFIsoNeutral03 + fRecoPhotonInfo2.rhocorPFIsoPhoton03
      ;

    ExoDiPhotons::InitDiphotonInfo(fDiphotonInfo);
    
    // fill diphoton info
    ExoDiPhotons::FillDiphotonInfo(fDiphotonInfo,&allTightOrFakeableObjects[0].first,&allTightOrFakeableObjects[1].first);

    //make an exception if there are 2 tight objects = top priority 

    if ( selectedPhotons.size() >=2 ){

      //since we want to have the PF isolation variables
      //we need edm::Ptr to the reco photons
      //and for that we need the index in the photon collection
      //but we lost track of it when sorting our vectors
      //so lets find it back

      if (fDebug) cout<<"here testing whether i get the correct pointer to my photon T"<<endl;
      reco::Photon* TObject1 = &(selectedPhotons[0]);
      reco::Photon* TObject2 = &(selectedPhotons[1]);
      if (fDebug) {
        cout<<TObject1->energy()<<" "
      <<TObject1->eta()<<" "
      <<TObject1->et()<<" "
      <<TObject1->phi()<<" "
      <<endl;
        cout<<TObject2->energy()<<" "
      <<TObject2->eta()<<" "
      <<TObject2->et()<<" "
      <<TObject2->phi()<<" "
      <<endl;
      }

      int indexTObject1 = -1;
      int indexTObject2 = -1;
      myIndex = -1;

      for(edm::View<reco::Photon>::const_iterator recoPhoton = photons->begin(); recoPhoton!=photons->end(); recoPhoton++) {
	myIndex++;
	//const reco::Photon* myPhoton = &(*recoPhoton);
	if(recoPhoton->pt() == TObject1->pt()) {
	  //if(myPhoton == TObject1) {
    if (fDebug) cout<<"Great, I've found t object 1 "<<endl;
	  indexTObject1 = myIndex;
	}
	if(recoPhoton->pt() == TObject2->pt()) {
	  //if(myPhoton == TObject2) {
    if (fDebug) cout<<"Great, I've found t object 2 "<<endl;
	  indexTObject2 = myIndex;
	}
      }

      //in principle the indices should be always greater than 1
      //because we necessarily found the two objects in the 
      //photon collection


      edm::Ptr<reco::Photon> TObject1Ptr(photons,indexTObject1);
      if (fDebug) {
        cout<<TObject1Ptr->energy()<<" "
      <<TObject1Ptr->eta()<<" "
      <<TObject1Ptr->et()<<" "
      <<TObject1Ptr->phi()<<" "
      <<endl;
      }

      edm::Ptr<reco::Photon> TObject2Ptr(photons,indexTObject2);
      if (fDebug) {
        cout<<TObject2Ptr->energy()<<" "
      <<TObject2Ptr->eta()<<" "
      <<TObject2Ptr->et()<<" "
      <<TObject2Ptr->phi()<<" "
      <<endl;
        cout<<"here testing whether i get the correct pointer to my photon T"<<endl;
      }

      // must specifically declare isFakeable status (should be Tight = not True = false                       
      ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo1,&selectedPhotons[0],lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);
      // fRecoPhotonInfo1.hasMatchedPromptElec = ConversionTools::hasMatchedPromptElectron((&selectedPhotons[0])->superCluster(), hElectrons, hConversions, beamSpot.position());
      fRecoPhotonInfo1.hasMatchedPromptElec = !( vetoMap.at((&selectedPhotons[0])->pt()) );

      fRecoPhotonInfo1.isTightPFPhoton = (*tight_id_decisions)[TObject1Ptr];
      fRecoPhotonInfo1.isMediumPFPhoton = (*medium_id_decisions)[TObject1Ptr];
      fRecoPhotonInfo1.isLoosePFPhoton = (*loose_id_decisions)[TObject1Ptr]; 

      //TO DISENTANGLE BETWEEN MINIAOD AND AOD

      fRecoPhotonInfo1.isFakeable = false;
      allTightOrFakeableObjects[0].second = false;//used in sorting, now it's faked if this 2 tight exception comes up

      //Now we store all PF isolation variables for the 1st photon (tight exception)

      fRecoPhotonInfo1.sigmaIetaIeta = (*full5x5SigmaIEtaIEtaMap)[TObject1Ptr];
      fRecoPhotonInfo1.PFIsoCharged03 = (*phoChargedIsolationMap)[TObject1Ptr];
      fRecoPhotonInfo1.PFIsoNeutral03 = (*phoNeutralHadronIsolationMap)[TObject1Ptr];
      fRecoPhotonInfo1.PFIsoPhoton03 = (*phoPhotonIsolationMap)[TObject1Ptr];
      fRecoPhotonInfo1.PFIsoAll03 = fRecoPhotonInfo1.PFIsoCharged03 + fRecoPhotonInfo1.PFIsoNeutral03 + fRecoPhotonInfo1.PFIsoPhoton03;

      float tObject1Eta = abs(TObject1Ptr->superCluster()->eta());

      CHeffarea = 0.;
      NHeffarea = 0.;
      PHeffarea = 0.;
      
      effareas.clear();
      effareas = ExoDiPhotons::EffectiveAreas((&selectedPhotons[0]),MethodID,CategoryPFID);
      if(MethodID.Contains("highpt")){
	CHeffarea = effareas[0];
	NHeffarea = effareas[1];
	PHeffarea = effareas[2];
      }
      if(MethodID.Contains("egamma")){
	CHeffarea = effAreaChHadrons_.getEffectiveArea(tObject1Eta);
	NHeffarea = effAreaNeuHadrons_.getEffectiveArea(tObject1Eta);
	PHeffarea = effAreaPhotons_.getEffectiveArea(tObject1Eta);
      }

      //now the corrected PF isolation variables
      fRecoPhotonInfo1.rhocorPFIsoCharged03 = std::max((float)0.0,(float)fRecoPhotonInfo1.PFIsoCharged03-rho_*CHeffarea);
      fRecoPhotonInfo1.rhocorPFIsoNeutral03 = std::max((float)0.0,(float)fRecoPhotonInfo1.PFIsoNeutral03-rho_*NHeffarea);
      fRecoPhotonInfo1.rhocorPFIsoPhoton03 = std::max((float)0.0,(float)fRecoPhotonInfo1.PFIsoPhoton03-rho_*PHeffarea);
      fRecoPhotonInfo1.rhocorPFIsoAll03 = fRecoPhotonInfo1.rhocorPFIsoCharged03 + fRecoPhotonInfo1.rhocorPFIsoNeutral03 + fRecoPhotonInfo1.rhocorPFIsoPhoton03;


      ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo2,&selectedPhotons[1],lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);
      // fRecoPhotonInfo2.hasMatchedPromptElec = ConversionTools::hasMatchedPromptElectron((&selectedPhotons[1])->superCluster(), hElectrons, hConversions, beamSpot.position());
      fRecoPhotonInfo2.hasMatchedPromptElec = !( vetoMap.at((&selectedPhotons[1])->pt()) );

      fRecoPhotonInfo2.isTightPFPhoton = (*tight_id_decisions)[TObject2Ptr];
      fRecoPhotonInfo2.isMediumPFPhoton = (*medium_id_decisions)[TObject2Ptr];
      fRecoPhotonInfo2.isLoosePFPhoton = (*loose_id_decisions)[TObject2Ptr]; 

      //TO DISENTANGLE BETWEEN MINIAOD AND AOD
      fRecoPhotonInfo2.isFakeable = false;
      allTightOrFakeableObjects[1].second = false;

      //Now we store all PF isolation variables for the 2st photon (tight exception)

      fRecoPhotonInfo2.sigmaIetaIeta = (*full5x5SigmaIEtaIEtaMap)[TObject2Ptr];
      fRecoPhotonInfo2.PFIsoCharged03 = (*phoChargedIsolationMap)[TObject2Ptr];
      fRecoPhotonInfo2.PFIsoNeutral03 = (*phoNeutralHadronIsolationMap)[TObject2Ptr];
      fRecoPhotonInfo2.PFIsoPhoton03 = (*phoPhotonIsolationMap)[TObject2Ptr];
      fRecoPhotonInfo2.PFIsoAll03 = fRecoPhotonInfo2.PFIsoCharged03 + fRecoPhotonInfo2.PFIsoNeutral03 + fRecoPhotonInfo2.PFIsoPhoton03;

      float tObject2Eta = abs(TObject2Ptr->superCluster()->eta());
      
      CHeffarea = 0.;
      NHeffarea = 0.;
      PHeffarea = 0.;
      
      effareas.clear();
      effareas = ExoDiPhotons::EffectiveAreas((&selectedPhotons[1]),MethodID,CategoryPFID);
      if(MethodID.Contains("highpt")){
	CHeffarea = effareas[0];
	NHeffarea = effareas[1];
	PHeffarea = effareas[2];
      }
      if(MethodID.Contains("egamma")){
	CHeffarea = effAreaChHadrons_.getEffectiveArea(tObject2Eta);
	NHeffarea = effAreaNeuHadrons_.getEffectiveArea(tObject2Eta);
	PHeffarea = effAreaPhotons_.getEffectiveArea(tObject2Eta);
      }

      //now the corrected PF isolation variables
      fRecoPhotonInfo2.rhocorPFIsoCharged03 = std::max((float)0.0,(float)fRecoPhotonInfo2.PFIsoCharged03-rho_*CHeffarea);
      fRecoPhotonInfo2.rhocorPFIsoNeutral03 = std::max((float)0.0,(float)fRecoPhotonInfo2.PFIsoNeutral03-rho_*NHeffarea);
      fRecoPhotonInfo2.rhocorPFIsoPhoton03 = std::max((float)0.0,(float)fRecoPhotonInfo2.PFIsoPhoton03-rho_*PHeffarea);
      fRecoPhotonInfo2.rhocorPFIsoAll03 = fRecoPhotonInfo2.rhocorPFIsoCharged03 + fRecoPhotonInfo2.rhocorPFIsoNeutral03 + fRecoPhotonInfo2.rhocorPFIsoPhoton03;

      ExoDiPhotons::InitDiphotonInfo(fDiphotonInfo);

      // fill diphoton info                                                                                                   
      ExoDiPhotons::FillDiphotonInfo(fDiphotonInfo,&selectedPhotons[0],&selectedPhotons[1]);

    }//end of 2 TT exception
     
    // for the pileup/vtx study.
    // lets try setting these photons to different vertices
    // and see how much these changes Mgg
     
    //      if(myVertices.size()>=2) {
    //        reco::Photon photon1_vtx2(allTightOrFakeableObjects[0].first);
    //        // keep calo poisiton same, but assign new vertex
    //        // this func also recalcs 4-vector after changing vertex!
    //        photon1_vtx2.setVertex(myVertices[1].position()); 

    //        reco::Photon photon2_vtx2(allTightOrFakeableObjects[1].first);
    //        photon2_vtx2.setVertex(myVertices[1].position()); 

    //        ExoDiPhotons::FillDiphotonInfo(fDiphotonInfoVtx2,&photon1_vtx2,&photon2_vtx2);

    //      }
     
    //      if(myVertices.size()>=3) {
    //        reco::Photon photon1_vtx3(allTightOrFakeableObjects[0].first);
    //        // keep calo poisiton same, but assign new vertex
    //        // this func also recalcs 4-vector after changing vertex!
    //        photon1_vtx3.setVertex(myVertices[2].position()); 

    //        reco::Photon photon2_vtx3(allTightOrFakeableObjects[1].first);
    //        photon2_vtx3.setVertex(myVertices[2].position()); 

    //        ExoDiPhotons::FillDiphotonInfo(fDiphotonInfoVtx3,&photon1_vtx3,&photon2_vtx3);

    //      }

    // now, we just decide which tree to fill: TT/TF/FF
    // remember the boolean is for isFakeable, so false is tight
    if(!allTightOrFakeableObjects[0].second && !allTightOrFakeableObjects[1].second) {
      // fill the tight-tight tree
      fTree->Fill();
      //       cout << "This event was TT" <<endl;
    }
    else if(!allTightOrFakeableObjects[0].second ) {
      // the BOTH tight option has been checked already
      // so now we just check the one-tight option
      // and fill the tight-fake tree instead
	   
      fTightFakeTree->Fill();
      //       cout << "This event was TF (or FT)" <<endl;
    }
    else if(!allTightOrFakeableObjects[1].second ) {
      fFakeTightTree->Fill();
    }
    else if(allTightOrFakeableObjects[0].second && allTightOrFakeableObjects[1].second) {
      // so if both isFakeable are true, then 
      // fill the fake-fake tree instead
      fRecoPhotonInfo1.isFakeable = true;
      fRecoPhotonInfo2.isFakeable = true;

      fFakeFakeTree->Fill();
      //       cout << "This event was FF" <<endl;
    }
    else {
      cout << "Neither TT, TF, FT nor FF?! Impossible!" <<endl;
    }

  } // end require two objects

  //
  //
  // The following is a duplicate of the above two-object block
  // However, it saves information for the leading three tight photons for every event, using fTree2
  //
  //
  if (fDebug) cout << ". starting charged decay code part two" << endl;
  InitRecoPhotonInfo(fRecoTightPhotonInfo1);
  InitRecoPhotonInfo(fRecoTightPhotonInfo2);
  InitRecoPhotonInfo(fRecoTightPhotonInfo3);
  int tightPhotonsCount = 0;
  for(unsigned int p = 0; p < allTightOrFakeableObjects.size(); p++)
  {
    if (!allTightOrFakeableObjects[p].second) tightPhotonsCount++;

    if (!allTightOrFakeableObjects[p].second && tightPhotonsCount == 1)
    {
      reco::Photon* TorFObject1 = &(allTightOrFakeableObjects[p].first);
      int indexTorFObject1 = -1;
      int myIndex = -1;
      for(edm::View<reco::Photon>::const_iterator recoPhoton = photons->begin(); recoPhoton!=photons->end(); recoPhoton++) {
        myIndex++;
        if(recoPhoton->pt() == TorFObject1->pt()) indexTorFObject1 = myIndex;
      }
      edm::Ptr<reco::Photon> TorFObject1Ptr(photons,indexTorFObject1);
      ExoDiPhotons::FillRecoPhotonInfo(fRecoTightPhotonInfo1,&allTightOrFakeableObjects[p].first,lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent,iSetup);
      fRecoTightPhotonInfo1.isFakeable = allTightOrFakeableObjects[p].second;
      fRecoTightPhotonInfo1.hasMatchedPromptElec = !(vetoMap.at((&allTightOrFakeableObjects[p].first)->pt())); 
      fRecoTightPhotonInfo1.isTightPFPhoton = (*tight_id_decisions)[TorFObject1Ptr];
      fRecoTightPhotonInfo1.isMediumPFPhoton = (*medium_id_decisions)[TorFObject1Ptr];
      fRecoTightPhotonInfo1.isLoosePFPhoton = (*loose_id_decisions)[TorFObject1Ptr]; 
      //TO DISENTANGLE BETWEEN MINIAOD AND AOD
      fRecoTightPhotonInfo1.sigmaIetaIeta = (*full5x5SigmaIEtaIEtaMap)[TorFObject1Ptr];
      fRecoTightPhotonInfo1.PFIsoCharged03 = (*phoChargedIsolationMap)[TorFObject1Ptr];
      fRecoTightPhotonInfo1.PFIsoNeutral03 = (*phoNeutralHadronIsolationMap)[TorFObject1Ptr];
      fRecoTightPhotonInfo1.PFIsoPhoton03 = (*phoPhotonIsolationMap)[TorFObject1Ptr];
      fRecoTightPhotonInfo1.PFIsoAll03 = fRecoTightPhotonInfo1.PFIsoCharged03 + fRecoTightPhotonInfo1.PFIsoNeutral03 + fRecoTightPhotonInfo1.PFIsoPhoton03;
      float torfObjectEta = abs(TorFObject1Ptr->superCluster()->eta());
      float CHeffarea = 0.;
      float NHeffarea = 0.;
      float PHeffarea = 0.;
      std::vector<double> effareas = ExoDiPhotons::EffectiveAreas((&allTightOrFakeableObjects[p].first),MethodID,CategoryPFID);
      if(MethodID.Contains("highpt")){
        CHeffarea = effareas[0];
        NHeffarea = effareas[1];
        PHeffarea = effareas[2];
      }
      if(MethodID.Contains("egamma")){
        CHeffarea = effAreaChHadrons_.getEffectiveArea(torfObjectEta);
        NHeffarea = effAreaNeuHadrons_.getEffectiveArea(torfObjectEta);
        PHeffarea = effAreaPhotons_.getEffectiveArea(torfObjectEta);
      }
      //now the corrected PF isolation variables
      fRecoTightPhotonInfo1.rhocorPFIsoCharged03 = std::max((float)0.0,(float)fRecoTightPhotonInfo1.PFIsoCharged03-rho_*CHeffarea);
      fRecoTightPhotonInfo1.rhocorPFIsoNeutral03 = std::max((float)0.0,(float)fRecoTightPhotonInfo1.PFIsoNeutral03-rho_*NHeffarea);
      fRecoTightPhotonInfo1.rhocorPFIsoPhoton03 = std::max((float)0.0,(float)fRecoTightPhotonInfo1.PFIsoPhoton03-rho_*PHeffarea);
      fRecoTightPhotonInfo1.rhocorPFIsoAll03 = fRecoTightPhotonInfo1.rhocorPFIsoCharged03 + fRecoTightPhotonInfo1.rhocorPFIsoNeutral03 + fRecoTightPhotonInfo1.rhocorPFIsoPhoton03;        
    }
    if (!allTightOrFakeableObjects[p].second && tightPhotonsCount == 2)
    {
      reco::Photon* TorFObject2 = &(allTightOrFakeableObjects[p].first);
      int indexTorFObject2 = -1;
      int myIndex = -1;
      for(edm::View<reco::Photon>::const_iterator recoPhoton = photons->begin(); recoPhoton!=photons->end(); recoPhoton++) {
        myIndex++;
        if(recoPhoton->pt() == TorFObject2->pt()) indexTorFObject2 = myIndex;
      }
      edm::Ptr<reco::Photon> TorFObject2Ptr(photons,indexTorFObject2);
      ExoDiPhotons::FillRecoPhotonInfo(fRecoTightPhotonInfo2,&allTightOrFakeableObjects[p].first,lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent,iSetup);
      fRecoTightPhotonInfo2.isFakeable = allTightOrFakeableObjects[p].second;
      fRecoTightPhotonInfo2.hasMatchedPromptElec = !(vetoMap.at((&allTightOrFakeableObjects[p].first)->pt())); 
      fRecoTightPhotonInfo2.isTightPFPhoton = (*tight_id_decisions)[TorFObject2Ptr];
      fRecoTightPhotonInfo2.isMediumPFPhoton = (*medium_id_decisions)[TorFObject2Ptr];
      fRecoTightPhotonInfo2.isLoosePFPhoton = (*loose_id_decisions)[TorFObject2Ptr]; 
      //TO DISENTANGLE BETWEEN MINIAOD AND AOD
      fRecoTightPhotonInfo2.sigmaIetaIeta = (*full5x5SigmaIEtaIEtaMap)[TorFObject2Ptr];
      fRecoTightPhotonInfo2.PFIsoCharged03 = (*phoChargedIsolationMap)[TorFObject2Ptr];
      fRecoTightPhotonInfo2.PFIsoNeutral03 = (*phoNeutralHadronIsolationMap)[TorFObject2Ptr];
      fRecoTightPhotonInfo2.PFIsoPhoton03 = (*phoPhotonIsolationMap)[TorFObject2Ptr];
      fRecoTightPhotonInfo2.PFIsoAll03 = fRecoTightPhotonInfo2.PFIsoCharged03 + fRecoTightPhotonInfo2.PFIsoNeutral03 + fRecoTightPhotonInfo2.PFIsoPhoton03;
      float torfObjectEta = abs(TorFObject2Ptr->superCluster()->eta());
      float CHeffarea = 0.;
      float NHeffarea = 0.;
      float PHeffarea = 0.;
      std::vector<double> effareas = ExoDiPhotons::EffectiveAreas((&allTightOrFakeableObjects[p].first),MethodID,CategoryPFID);
      if(MethodID.Contains("highpt")){
        CHeffarea = effareas[0];
        NHeffarea = effareas[1];
        PHeffarea = effareas[2];
      }
      if(MethodID.Contains("egamma")){
        CHeffarea = effAreaChHadrons_.getEffectiveArea(torfObjectEta);
        NHeffarea = effAreaNeuHadrons_.getEffectiveArea(torfObjectEta);
        PHeffarea = effAreaPhotons_.getEffectiveArea(torfObjectEta);
      }
      //now the corrected PF isolation variables
      fRecoTightPhotonInfo2.rhocorPFIsoCharged03 = std::max((float)0.0,(float)fRecoTightPhotonInfo2.PFIsoCharged03-rho_*CHeffarea);
      fRecoTightPhotonInfo2.rhocorPFIsoNeutral03 = std::max((float)0.0,(float)fRecoTightPhotonInfo2.PFIsoNeutral03-rho_*NHeffarea);
      fRecoTightPhotonInfo2.rhocorPFIsoPhoton03 = std::max((float)0.0,(float)fRecoTightPhotonInfo2.PFIsoPhoton03-rho_*PHeffarea);
      fRecoTightPhotonInfo2.rhocorPFIsoAll03 = fRecoTightPhotonInfo2.rhocorPFIsoCharged03 + fRecoTightPhotonInfo2.rhocorPFIsoNeutral03 + fRecoTightPhotonInfo2.rhocorPFIsoPhoton03;        
    }
    if (!allTightOrFakeableObjects[p].second && tightPhotonsCount == 3)
    {
      reco::Photon* TorFObject3 = &(allTightOrFakeableObjects[p].first);
      int indexTorFObject3 = -1;
      int myIndex = -1;
      for(edm::View<reco::Photon>::const_iterator recoPhoton = photons->begin(); recoPhoton!=photons->end(); recoPhoton++) {
        myIndex++;
        if(recoPhoton->pt() == TorFObject3->pt()) indexTorFObject3 = myIndex;
      }
      edm::Ptr<reco::Photon> TorFObject2Ptr(photons,indexTorFObject3);
      ExoDiPhotons::FillRecoPhotonInfo(fRecoTightPhotonInfo3,&allTightOrFakeableObjects[p].first,lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent,iSetup);
      fRecoTightPhotonInfo3.isFakeable = allTightOrFakeableObjects[p].second;
      fRecoTightPhotonInfo3.hasMatchedPromptElec = !(vetoMap.at((&allTightOrFakeableObjects[p].first)->pt())); 
      fRecoTightPhotonInfo3.isTightPFPhoton = (*tight_id_decisions)[TorFObject2Ptr];
      fRecoTightPhotonInfo3.isMediumPFPhoton = (*medium_id_decisions)[TorFObject2Ptr];
      fRecoTightPhotonInfo3.isLoosePFPhoton = (*loose_id_decisions)[TorFObject2Ptr]; 
      //TO DISENTANGLE BETWEEN MINIAOD AND AOD
      fRecoTightPhotonInfo3.sigmaIetaIeta = (*full5x5SigmaIEtaIEtaMap)[TorFObject2Ptr];
      fRecoTightPhotonInfo3.PFIsoCharged03 = (*phoChargedIsolationMap)[TorFObject2Ptr];
      fRecoTightPhotonInfo3.PFIsoNeutral03 = (*phoNeutralHadronIsolationMap)[TorFObject2Ptr];
      fRecoTightPhotonInfo3.PFIsoPhoton03 = (*phoPhotonIsolationMap)[TorFObject2Ptr];
      fRecoTightPhotonInfo3.PFIsoAll03 = fRecoTightPhotonInfo3.PFIsoCharged03 + fRecoTightPhotonInfo3.PFIsoNeutral03 + fRecoTightPhotonInfo3.PFIsoPhoton03;
      float torfObjectEta = abs(TorFObject2Ptr->superCluster()->eta());
      float CHeffarea = 0.;
      float NHeffarea = 0.;
      float PHeffarea = 0.;
      std::vector<double> effareas = ExoDiPhotons::EffectiveAreas((&allTightOrFakeableObjects[p].first),MethodID,CategoryPFID);
      if(MethodID.Contains("highpt")){
        CHeffarea = effareas[0];
        NHeffarea = effareas[1];
        PHeffarea = effareas[2];
      }
      if(MethodID.Contains("egamma")){
        CHeffarea = effAreaChHadrons_.getEffectiveArea(torfObjectEta);
        NHeffarea = effAreaNeuHadrons_.getEffectiveArea(torfObjectEta);
        PHeffarea = effAreaPhotons_.getEffectiveArea(torfObjectEta);
      }
      //now the corrected PF isolation variables
      fRecoTightPhotonInfo3.rhocorPFIsoCharged03 = std::max((float)0.0,(float)fRecoTightPhotonInfo3.PFIsoCharged03-rho_*CHeffarea);
      fRecoTightPhotonInfo3.rhocorPFIsoNeutral03 = std::max((float)0.0,(float)fRecoTightPhotonInfo3.PFIsoNeutral03-rho_*NHeffarea);
      fRecoTightPhotonInfo3.rhocorPFIsoPhoton03 = std::max((float)0.0,(float)fRecoTightPhotonInfo3.PFIsoPhoton03-rho_*PHeffarea);
      fRecoTightPhotonInfo3.rhocorPFIsoAll03 = fRecoTightPhotonInfo3.rhocorPFIsoCharged03 + fRecoTightPhotonInfo3.rhocorPFIsoNeutral03 + fRecoTightPhotonInfo3.rhocorPFIsoPhoton03;        
    }
  } // end loop over tight or fakable collection 
  if (fDebug) cout << ". done making tight photons" << endl;
  fNumTightPhotons = tightPhotonsCount;
  // Construct Di-Objects
  InitRecoDiObjectInfo(fTwoProngTwoProngInfo);
  InitRecoDiObjectInfo(fGammaTwoProngInfo);
  InitRecoDiObjectInfo(fTwoProngFakeInfo);
  InitRecoDiObjectInfo(fGammaFakeInfo);
  InitRecoDiObjectInfo(fGammaGammaInfo);
  // Passed Eta and Passed Eta
  if (fNumCHpairsPass >= 2)
  {
    TLorentzVector Eta1;
    Eta1.SetPtEtaPhiM(fTwoProng_pt[0], fTwoProng_eta[0], fTwoProng_phi[0], fTwoProng_mass[0]);
    TLorentzVector Eta2;
    Eta2.SetPtEtaPhiM(fTwoProng_pt[1], fTwoProng_eta[1], fTwoProng_phi[1], fTwoProng_mass[1]);
    FillRecoDiObjectInfo(fTwoProngTwoProngInfo, Eta1, Eta2);
    fTwoProngTwoProngInfo.dMass = fabs(fTwoProng_Mass[0] - fTwoProng_Mass[1]);
  }
  if (fDebug) cout << ". done making TwoProng TwoProng" << endl;
  // Passed Eta and Fake Eta
  if (fNumCHpairsPass >= 1 && fNumCHpairsFake >= 1)
  {
    TLorentzVector Eta;
    Eta.SetPtEtaPhiM(fTwoProng_pt[0], fTwoProng_eta[0], fTwoProng_phi[0], fTwoProng_mass[0]);
    FillRecoDiObjectInfo(fTwoProngFakeInfo, Eta, LeadingFake);
    fTwoProngFakeInfo.dMass = fabs(fTwoProng_Mass[0] - LeadingFake.M());
  }
  if (fDebug) cout << ". done making TwoProng Fake" << endl;
  // Tight Photon and Passed Eta
  if (fNumCHpairsPass >= 1 && fNumTightPhotons >= 1)
  {
    TLorentzVector Eta;
    Eta.SetPtEtaPhiM(fTwoProng_pt[0], fTwoProng_eta[0], fTwoProng_phi[0], fTwoProng_mass[0]);
    TLorentzVector Photon;
    Photon.SetPtEtaPhiM(fRecoTightPhotonInfo1.pt, fRecoTightPhotonInfo1.eta, fRecoTightPhotonInfo1.phi, 0);
    FillRecoDiObjectInfo(fGammaTwoProngInfo, Photon, Eta);
    fGammaTwoProngInfo.dMass = fabs(fTwoProng_Mass[0] - 0);
  }
  if (fDebug) cout << ". done making Gamma TwoProng" << endl;
  // Tight Photon and Fake Eta
  if (fNumCHpairsFake >= 1 && fNumTightPhotons >= 1)
  {
    TLorentzVector Photon;
    Photon.SetPtEtaPhiM(fRecoTightPhotonInfo1.pt, fRecoTightPhotonInfo1.eta, fRecoTightPhotonInfo1.phi, 0);
    FillRecoDiObjectInfo(fGammaFakeInfo, Photon, LeadingFake);
    fGammaFakeInfo.dMass = fabs(LeadingFake.M() - 0);
  }
  if (fDebug) cout << ". done making Gamma Fake" << endl;
  // Tight Photon and Tight Photon
  if (fNumTightPhotons >= 2)
  {
    TLorentzVector Photon1;
    Photon1.SetPtEtaPhiM(fRecoTightPhotonInfo1.pt, fRecoTightPhotonInfo1.eta, fRecoTightPhotonInfo1.phi, 0);
    TLorentzVector Photon2;
    Photon2.SetPtEtaPhiM(fRecoTightPhotonInfo2.pt, fRecoTightPhotonInfo2.eta, fRecoTightPhotonInfo2.phi, 0);
    FillRecoDiObjectInfo(fGammaGammaInfo, Photon1, Photon2);
  }
  if (fDebug) cout << ". done making Gamma Gamma" << endl;
  if (fDebug) cout << ". finished charged decay part two" << endl;

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

  // Cutflow for two prong selection


  // Now fill fTree2, it's filled for every event
  fTree2->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
ExoDiPhotonAnalyzer::beginJob()
{
  if (fisMC){
    LumiWeights = edm::LumiReWeighting(fPUMCFileName,fPUDataFileName,fPUMCHistName,fPUDataHistName);
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ExoDiPhotonAnalyzer::endJob()
{
  // fake rate histograms
  fTwoProngFakeNumer_pt->Sumw2();
  fTwoProngFakeDenom_pt->Sumw2();
  fTwoProngFakeRate_pt->Add(fTwoProngFakeNumer_pt);
  fTwoProngFakeRate_pt->Divide(fTwoProngFakeDenom_pt);
  fTwoProngFakeNumer_eta->Sumw2();
  fTwoProngFakeDenom_eta->Sumw2();
  fTwoProngFakeRate_eta->Add(fTwoProngFakeNumer_eta);
  fTwoProngFakeRate_eta->Divide(fTwoProngFakeDenom_eta);
  fTwoProngFakeNumer_phi->Sumw2();
  fTwoProngFakeDenom_phi->Sumw2();
  fTwoProngFakeRate_phi->Add(fTwoProngFakeNumer_phi);
  fTwoProngFakeRate_phi->Divide(fTwoProngFakeDenom_phi);

  // Print Cutflow
  if (fchargedDecayCutflow) {
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



//define this as a plug-in
DEFINE_FWK_MODULE(ExoDiPhotonAnalyzer);
