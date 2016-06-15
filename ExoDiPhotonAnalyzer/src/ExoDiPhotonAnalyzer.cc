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
  TH1F *fTwoProngFakeRate_pt;
  TH1F *fTwoProngFakeRate_eta;
  TH1F *fTwoProngFakeRate_phi;
  TH1F *fTwoProngFakeNumer_pt = new TH1F("twoprongfakenumer_eta","Fake Numerator count for CH pairs, inverted charged iso cut, eta binned",100,0.0,2000.0);
  TH1F *fTwoProngFakeDenom_pt = new TH1F("twoprongfakedenom_eta","Fake Denominator count for CH pairs, inverted charged iso cut, eta binned",100,0.0,2000.0);
  TH1F *fTwoProngFakeNumer_eta = new TH1F("twoprongfakenumer_eta","Fake Numerator count for CH pairs, inverted charged iso cut, eta binned",100,-10.0,10.0);
  TH1F *fTwoProngFakeDenom_eta = new TH1F("twoprongfakedenom_eta","Fake Denominator count for CH pairs, inverted charged iso cut, eta binned",100,-10.0,10.0);
  TH1F *fTwoProngFakeNumer_phi = new TH1F("twoprongfakenumer_phi","Fake Numerator count for CH pairs, inverted charged iso cut, phi binned",100,-10.0,10.0);
  TH1F *fTwoProngFakeDenom_phi = new TH1F("twoprongfakedenom_phi","Fake Denominator count for CH pairs, inverted charged iso cut, phi binned",100,-10.0,10.0);

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
  // Branch variables
  int fEventNum;
  int fRunNum;
  int fLumiNum;
  int fNumPVs;
  int fNumPF;
  int fNumPrunedPF;
  int fNumCHpairs;
  int fNumCHpairsPass;
  int fNumCHpairsFake;
  int fNumCHpairsMatch;
  int fNumCHpairsPassMatch;
  int fNumGenEta;
  ExoDiPhotons::recoTwoProngInfo_t fRecoTwoProngInfo1;
  ExoDiPhotons::recoTwoProngInfo_t fRecoTwoProngInfo2;
  ExoDiPhotons::genParticleInfo_t fGenEtaParticleInfo1;
  ExoDiPhotons::genParticleInfo_t fGenEtaParticleInfo2;
  double fMass_EtaEta;
  double fMass2_EtaEta;
  double fDR_EtaEta;
  double fDPhi_EtaEta;
  double fDEta_EtaEta;
  TLorentzVector fEtaEtaCand;
  TLorentzVector fEtaEtaCand2;
  vector<Double_t> fCands_pt;
  vector<Double_t> fCands_eta;
  vector<Double_t> fCands_phi;
  vector<Double_t> fCands_mass;
  vector<Double_t> fCands_px;
  vector<Double_t> fCands_py;
  vector<Double_t> fCands_pz;
  vector<Double_t> fCands_energy;
  vector<Double_t> fCands_chargediso;
  vector<Double_t> fCands_neutraliso;
  vector<Double_t> fCands_egammaiso;
  vector<Double_t> fPasses_pt;
  vector<Double_t> fPasses_eta;
  vector<Double_t> fPasses_phi;
  vector<Double_t> fPasses_mass;
  vector<Double_t> fPasses_px;
  vector<Double_t> fPasses_py;
  vector<Double_t> fPasses_pz;
  vector<Double_t> fPasses_energy;
  vector<Double_t> fFakes_pt;
  vector<Double_t> fFakes_eta;
  vector<Double_t> fFakes_phi;
  vector<Double_t> fFakes_mass;
  vector<Double_t> fFakes_px;
  vector<Double_t> fFakes_py;
  vector<Double_t> fFakes_pz;
  vector<Double_t> fFakes_energy;
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
  // Cutflow
  fTree2->Branch("nCands",&fNumCHpairs,"nCands/I");
  fTree2->Branch("nPass",&fNumCHpairsPass,"nPass/I");
  fTree2->Branch("nFake",&fNumCHpairsFake,"nFake/I");
  fTree2->Branch("nMatch",&fNumCHpairsMatch,"nMatch/I");
  fTree2->Branch("nPassMatch",&fNumCHpairsPassMatch,"nPassMatch/I");
  fTree2->Branch("nGenEta",&fNumGenEta,"nGenEta/I");
  // Kinematics and Particle Info
  fTree2->Branch("Eta1",&fRecoTwoProngInfo1,ExoDiPhotons::recoTwoProngBranchDefString.c_str());
  fTree2->Branch("Eta2",&fRecoTwoProngInfo2,ExoDiPhotons::recoTwoProngBranchDefString.c_str());
  fTree2->Branch("dREtaEta",&fDR_EtaEta,"dREtaEta/D");
  fTree2->Branch("dPhiEtaEta",&fDPhi_EtaEta,"dPhiEtaEta/D");
  fTree2->Branch("dEtaEtaEta",&fDEta_EtaEta,"dEtaEtaEta/D");
  fTree2->Branch("mEtaEta",&fMass_EtaEta,"mEtaEta/D");
  fTree2->Branch("mGrommedEtaEta",&fMass2_EtaEta,"mGrommedEtaEta/D");
  fTree2->Branch("vecEtaEta",&fEtaEtaCand);
  fTree2->Branch("vecGrommedEtaEta",&fEtaEtaCand2);
  fTree2->Branch("cands_pt",&fCands_pt);
  fTree2->Branch("cands_eta",&fCands_eta);
  fTree2->Branch("cands_phi",&fCands_phi);
  fTree2->Branch("cands_mass",&fCands_mass);
  fTree2->Branch("cands_px",&fCands_px);
  fTree2->Branch("cands_py",&fCands_py);
  fTree2->Branch("cands_pz",&fCands_pz);
  fTree2->Branch("cands_energy",&fCands_energy);
  fTree2->Branch("cands_chargediso",&fCands_chargediso);
  fTree2->Branch("cands_neutraliso",&fCands_neutraliso);
  fTree2->Branch("cands_egammaiso",&fCands_egammaiso);
  fTree2->Branch("passes_pt",&fPasses_pt);
  fTree2->Branch("passes_eta",&fPasses_eta);
  fTree2->Branch("passes_phi",&fPasses_phi);
  fTree2->Branch("passes_mass",&fPasses_mass);
  fTree2->Branch("passes_px",&fPasses_px);
  fTree2->Branch("passes_py",&fPasses_py);
  fTree2->Branch("passes_pz",&fPasses_pz);
  fTree2->Branch("passes_energy",&fPasses_energy);
  fTree2->Branch("fakes_pt",&fFakes_pt);
  fTree2->Branch("fakes_eta",&fFakes_eta);
  fTree2->Branch("fakes_phi",&fFakes_phi);
  fTree2->Branch("fakes_mass",&fFakes_mass);
  fTree2->Branch("fakes_px",&fFakes_px);
  fTree2->Branch("fakes_py",&fFakes_py);
  fTree2->Branch("fakes_pz",&fFakes_pz);
  fTree2->Branch("fakes_energy",&fFakes_energy);
  // Generator level
  fTree2->Branch("genEta1",&fGenEtaParticleInfo1,ExoDiPhotons::genParticleInfoBranchDefString.c_str());
  fTree2->Branch("genEta2",&fGenEtaParticleInfo2,ExoDiPhotons::genParticleInfoBranchDefString.c_str());
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

  // Fill branches for charged decay analysis
  if (fDebug) cout << "start charged decay code" << endl;
  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(pfcandsToken_, pfcands);

  edm::Handle<vector<reco::GenParticle>> genparticles;
  if (fisMC) {
    iEvent.getByToken(genToken_, genparticles);
  }

  edm::Handle<vector<reco::Vertex>> primaryvertecies;
  iEvent.getByToken(pvToken_, primaryvertecies);

  fCands_pt.clear();
  fCands_eta.clear();
  fCands_phi.clear();
  fCands_mass.clear();
  fCands_px.clear();
  fCands_py.clear();
  fCands_pz.clear();
  fCands_energy.clear();
  fCands_chargediso.clear();
  fCands_neutraliso.clear();
  fCands_egammaiso.clear();
  fPasses_pt.clear();
  fPasses_eta.clear();
  fPasses_phi.clear();
  fPasses_mass.clear();
  fPasses_px.clear();
  fPasses_py.clear();
  fPasses_pz.clear();
  fPasses_energy.clear();
  fFakes_pt.clear();
  fFakes_eta.clear();
  fFakes_phi.clear();
  fFakes_mass.clear();
  fFakes_px.clear();
  fFakes_py.clear();
  fFakes_pz.clear();
  fFakes_energy.clear();

  fNumPVs = primaryvertecies->size();
  fNumPF = pfcands->size();

  // Prepare generator Eta collection
  if (fisMC && fisSignal) {
    ExoDiPhotons::InitGenParticleInfo(fGenEtaParticleInfo1);
    ExoDiPhotons::InitGenParticleInfo(fGenEtaParticleInfo2);
    int genEtaCount = 0;
    for (unsigned int i = 0; i < genparticles->size(); ++i) {
      const reco::GenParticle &gen = (*genparticles)[i];
      if (gen.status() == 2 && abs(gen.pdgId()) == 221) {
        genEtaCount += 1;
        if (genEtaCount == 1)
          ExoDiPhotons::FillGenParticleInfo(fGenEtaParticleInfo1, gen);
        if (genEtaCount == 2)
          ExoDiPhotons::FillGenParticleInfo(fGenEtaParticleInfo2, gen);
      }
    } // end gen particle loop
    fNumGenEta = genEtaCount;
  } else {
    ExoDiPhotons::InitGenParticleInfo(fGenEtaParticleInfo1);
    ExoDiPhotons::InitGenParticleInfo(fGenEtaParticleInfo2);
    fNumGenEta = -1;
  }

  // Find all pairs of one CH pos and one CH neg within specified DR of each other
  vector<TLorentzVector> candidates_CHpos;
  vector<unsigned int> candidates_CHposIndex;
  vector<TLorentzVector> candidates_CHneg;
  vector<unsigned int> candidates_CHnegIndex;
  vector<TLorentzVector> candidates_center;
  vector<TLorentzVector> candidates_photon;
  vector<int> candidates_numgamma;
  vector<int> candidates_nume;
  vector<TLorentzVector> candidates_Eta;
  vector<TLorentzVector> candidates_EtaGrommed;
  int pruned_count = 0;
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
          if ((pf3.pdgId() != 22) && (pf3.pdgId() != 11)) continue; // electron or photon pf
          TLorentzVector pfcand3;
          pfcand3.SetPtEtaPhiE(pf3.pt(), pf3.eta(), pf3.phiAtVtx(), pf3.energy());
          if (abs(pf3.phiAtVtx() - center.Phi()) < fCandidatePairPhiBox/2.0 &&
              abs(pf3.eta() - center.Eta()) < fCandidatePairEtaBox/2.0) {
            if (pf3.pdgId() == 22) {
              numgamma += 1;
              // find leading photon pf
              if (pf3.pt() > pt_of_leading_pf_photon) {
                pt_of_leading_pf_photon = pf3.pt();
                index_of_leading_pf_photon = k;
              }
            }
            else nume += 1;
            photon = photon + pfcand3;
          }
        } // end pf cand loop
        if (fDebug) cout << ". finished photon" << endl;
        int n = index_of_leading_pf_photon;
        if (n != -1) {
          // CH pair is close and has at least one pf photon, definition of candidate, fill vectors
          leading_pf_photon.SetPtEtaPhiE((*pfcands)[n].pt(), (*pfcands)[n].eta(), (*pfcands)[n].phiAtVtx(), (*pfcands)[n].energy());
          TLorentzVector EtaCandidate;
          EtaCandidate = center + photon;
          TLorentzVector EtaGrommedCandidate;
          EtaGrommedCandidate = center + leading_pf_photon;
          if (pf1.pdgId() > 0) {
            candidates_CHpos.push_back(pfcand1);
            candidates_CHposIndex.push_back(i);
            candidates_CHneg.push_back(pfcand2);
            candidates_CHnegIndex.push_back(j);
          } else {
            candidates_CHpos.push_back(pfcand2);
            candidates_CHposIndex.push_back(j);
            candidates_CHneg.push_back(pfcand1);
            candidates_CHnegIndex.push_back(i);
          }
          candidates_center.push_back(center);
          candidates_photon.push_back(photon);
          candidates_numgamma.push_back(numgamma);
          candidates_nume.push_back(nume);
          candidates_Eta.push_back(EtaCandidate);
          candidates_EtaGrommed.push_back(EtaGrommedCandidate);
          // and fill branches
          fCands_pt.push_back(EtaCandidate.Pt());
          fCands_eta.push_back(EtaCandidate.Eta());
          fCands_phi.push_back(EtaCandidate.Phi());
          fCands_mass.push_back(EtaCandidate.M());
          fCands_px.push_back(EtaCandidate.Px());
          fCands_py.push_back(EtaCandidate.Py());
          fCands_pz.push_back(EtaCandidate.Pz());
          fCands_energy.push_back(EtaCandidate.E());
        }
      } // end conditionals on CH pair
    }
  } // end double loop over pf
  if (fDebug) cout << ". finish pairs finding" << endl;

  // Associate data to candidate pairs
  vector<bool> candidates_passed;
  vector<bool> candidates_matched;
  vector<double> candidates_nearestGenDR;
  vector<int> candidates_nearestGenIndex;
  vector<double> candidates_relchargediso;
  vector<double> candidates_relneutraliso;
  vector<double> candidates_relegammaiso;
  int cand_pairs_passed = 0;
  int cand_pairs_faked = 0;
  int cand_pairs_matched = 0;
  int cand_pairs_passed_matched = 0;
  int cand_pairs_passed_charged = 0;
  int cand_pairs_passed_neutral = 0;
  int cand_pairs_passed_egamma = 0;
  int cand_pairs_passed_photonpt = 0;
  for (unsigned int i = 0; i < candidates_center.size(); i++) {
    TLorentzVector center;
    center.SetPtEtaPhiM(candidates_center[i].Pt(), candidates_center[i].Eta(), candidates_center[i].Phi(), candidates_center[i].M());
    TLorentzVector EtaCandidate;
    EtaCandidate.SetPtEtaPhiM(candidates_Eta[i].Pt(), candidates_Eta[i].Eta(), candidates_Eta[i].Phi(), candidates_Eta[i].M());
    TLorentzVector photon;
    photon.SetPtEtaPhiM(candidates_photon[i].Pt(), candidates_photon[i].Eta(), candidates_photon[i].Phi(), candidates_photon[i].M());
    if (fDebug) cout << ". made Tvectors" << endl;

    // Define isolations
    double chargedIso = 0;
    double neutralIso = 0;
    double egammaIso = 0;
    for (unsigned int j = 0; j < pfcands->size(); ++j) {
      const pat::PackedCandidate &pf = (*pfcands)[j];
      TLorentzVector pfcand;
      pfcand.SetPtEtaPhiE(pf.pt(), pf.eta(), pf.phiAtVtx(), pf.energy());
      if (!fisSignal && pf.fromPV() <= 1) continue; // for signal, fromPV() is broken
      // charged (incl. muons)
      if (abs(pf.pdgId()) == 13 || abs(pf.pdgId()) == 211) {
        if ( center.DeltaR(pfcand) < fCandidatePairIsolationDR && !(j == candidates_CHposIndex[i] || j == candidates_CHnegIndex[i]) )
            chargedIso += pfcand.Pt();
      // neutral hadron
      } else if (pf.pdgId() == 130) {
        if (center.DeltaR(pfcand) < fCandidatePairIsolationDR)
          neutralIso += pfcand.Pt();
      // e gamma
      } else if (abs(pf.pdgId()) == 11 || pf.pdgId() == 22) {
        if ( (center.DeltaR(pfcand) < fCandidatePairIsolationDR) &&
             !(abs(pf.phiAtVtx() - center.Phi()) < fCandidatePairPhiBox/2.0 && abs(pf.eta() - center.Eta()) < fCandidatePairEtaBox/2.0))
          egammaIso += pfcand.Pt();
      }
    } // end pf cand loop
    if (fDebug) cout << ". finished isolation" << endl;

    // Selection on Eta Candidate
    double relchargedIso = chargedIso / EtaCandidate.Pt();
    double relneutralIso = neutralIso / EtaCandidate.Pt();
    double relegammaIso = egammaIso / EtaCandidate.Pt();
    fCands_chargediso.push_back(relchargedIso);
    fCands_neutraliso.push_back(relneutralIso);
    fCands_egammaiso.push_back(relegammaIso);
    bool passed = relchargedIso < fCandidatePairChargedIsoCut &&
                  relneutralIso < fCandidatePairNeutralIsoCut &&
                  relegammaIso < fCandidatePairEGammaIsoCut &&
                  photon.Pt() > fCandidatePairPhotonPtCut;
    bool passed_fake = relchargedIso > fCandidatePairChargedIsoCut &&
                       relneutralIso < fCandidatePairNeutralIsoCut &&
                       relegammaIso < fCandidatePairEGammaIsoCut &&
                       photon.Pt() > fCandidatePairPhotonPtCut;
    bool passedCharged = relchargedIso < fCandidatePairChargedIsoCut;
    bool passedNeutral = relneutralIso < fCandidatePairNeutralIsoCut;
    bool passedEGamma = relegammaIso < fCandidatePairEGammaIsoCut;
    bool passedPhotonPt = photon.Pt() > fCandidatePairPhotonPtCut;
    if (passed) {
      fTwoProngFakeNumer_pt->Fill(EtaCandidate.Pt());
      fTwoProngFakeNumer_eta->Fill(EtaCandidate.Eta());
      fTwoProngFakeNumer_phi->Fill(EtaCandidate.Phi());
      // store in vector of passes
      fPasses_pt.push_back(EtaCandidate.Pt());
      fPasses_eta.push_back(EtaCandidate.Eta());
      fPasses_phi.push_back(EtaCandidate.Phi());
      fPasses_mass.push_back(EtaCandidate.M());
      fPasses_px.push_back(EtaCandidate.Px());
      fPasses_py.push_back(EtaCandidate.Py());
      fPasses_pz.push_back(EtaCandidate.Pz());
      fPasses_energy.push_back(EtaCandidate.E());
    }
    if (passed_fake) {
      fTwoProngFakeDenom_pt->Fill(EtaCandidate.Pt());
      fTwoProngFakeDenom_eta->Fill(EtaCandidate.Eta());
      fTwoProngFakeDenom_phi->Fill(EtaCandidate.Phi());
      // store in vector of fakes
      fFakes_pt.push_back(EtaCandidate.Pt());
      fFakes_eta.push_back(EtaCandidate.Eta());
      fFakes_phi.push_back(EtaCandidate.Phi());
      fFakes_mass.push_back(EtaCandidate.M());
      fFakes_px.push_back(EtaCandidate.Px());
      fFakes_py.push_back(EtaCandidate.Py());
      fFakes_pz.push_back(EtaCandidate.Pz());
      fFakes_energy.push_back(EtaCandidate.E());
    }
    if (fDebug) cout << ". finished selection" << endl;

    // Generator Matching
    bool matched = false;
    double mingenDR = 99.9;
    int index = 99;
    if (fisMC && fisSignal && fNumGenEta>=2) {
      TLorentzVector genEta1;
      genEta1.SetPtEtaPhiE(fGenEtaParticleInfo1.pt, fGenEtaParticleInfo1.eta, fGenEtaParticleInfo1.phi, fGenEtaParticleInfo1.energy);
      TLorentzVector genEta2;
      genEta2.SetPtEtaPhiE(fGenEtaParticleInfo2.pt, fGenEtaParticleInfo2.eta, fGenEtaParticleInfo2.phi, fGenEtaParticleInfo2.energy);
      double DR1 = EtaCandidate.DeltaR(genEta1);
      double DR2 = EtaCandidate.DeltaR(genEta2);
      if (DR1 < DR2) {
        mingenDR = DR1;
        index = 1;
      } else {
        mingenDR = DR2;
        index = 2;
      }
      if (mingenDR < fCandidatePairGenMatchDR) matched = true;
    }
    if (fDebug) cout << ". finished matching" << endl;
    // Fill all vectors for each candidate all at once
    candidates_relchargediso.push_back(relchargedIso);
    candidates_relneutraliso.push_back(relneutralIso);
    candidates_relegammaiso.push_back(relegammaIso);
    candidates_passed.push_back(passed);
    candidates_nearestGenDR.push_back(mingenDR);
    candidates_nearestGenIndex.push_back(index);
    candidates_matched.push_back(matched);
    if (passed) cand_pairs_passed += 1;
    if (passed_fake) cand_pairs_faked += 1;
    if (matched) cand_pairs_matched += 1;
    if (passed && matched) cand_pairs_passed_matched += 1;
    if (passedCharged) cand_pairs_passed_charged += 1;
    if (passedNeutral) cand_pairs_passed_neutral += 1;
    if (passedEGamma) cand_pairs_passed_egamma += 1;
    if (passedPhotonPt) cand_pairs_passed_photonpt += 1;
  } // end candidate loop
  fNumPrunedPF = pruned_count;
  fNumCHpairs = candidates_center.size();
  fNumCHpairsPass = cand_pairs_passed;
  fNumCHpairsFake = cand_pairs_faked;
  fNumCHpairsMatch = cand_pairs_matched;
  fNumCHpairsPassMatch = cand_pairs_passed_matched;
  if (fDebug) cout << ". finished all cand data" << endl;
  fTwoProngFakeRate_pt->Add(fTwoProngFakeNumer_pt);
  fTwoProngFakeRate_pt->Divide(fTwoProngFakeDenom_pt);
  fTwoProngFakeRate_eta->Add(fTwoProngFakeNumer_eta);
  fTwoProngFakeRate_eta->Divide(fTwoProngFakeDenom_eta);
  fTwoProngFakeRate_phi->Add(fTwoProngFakeNumer_phi);
  fTwoProngFakeRate_phi->Divide(fTwoProngFakeDenom_phi);
  if (fDebug) cout << ". finished fake rate histos" << endl;

  // For leading two passing charged hadron pairs, save info into ttree
  int leading_index = -1;
  double leading_pt = -1;
  int subleading_index = -1;
  double subleading_pt = -1;
  for (unsigned int i = 0; i < candidates_passed.size(); i++) {
    if (!candidates_passed[i]) continue;
    if (candidates_Eta[i].Pt() > leading_pt) {
      subleading_index = leading_index;
      subleading_pt = leading_pt;
      leading_index = i;
      leading_pt = candidates_Eta[i].Pt();
    } else if (candidates_Eta[i].Pt() > subleading_pt) {
      subleading_index = i;
      subleading_pt = candidates_Eta[i].Pt();
    }
  } // end candidate loop
  if (fDebug) cout << ". finished sorting leading two" << endl;
  ExoDiPhotons::InitRecoTwoProngInfo(fRecoTwoProngInfo1);
  ExoDiPhotons::InitRecoTwoProngInfo(fRecoTwoProngInfo2);
  fMass_EtaEta = -99.9;
  fMass2_EtaEta = -99.9;
  fDR_EtaEta = -99.9;
  fDPhi_EtaEta = -99.9;
  fDEta_EtaEta = -99.9;
  int l = leading_index;
  int s = subleading_index;
  if (cand_pairs_passed == 1) {
    ExoDiPhotons::FillRecoTwoProngInfo(fRecoTwoProngInfo1, candidates_CHpos[l], candidates_CHneg[l], candidates_center[l],
                                       candidates_photon[l], candidates_Eta[l], candidates_EtaGrommed[l].M(), candidates_passed[l], candidates_matched[l], 
                                       candidates_nearestGenDR[l], candidates_nearestGenIndex[l], candidates_relchargediso[l],
                                       candidates_relneutraliso[l], candidates_relegammaiso[l], candidates_numgamma[l], candidates_nume[l]);
    ExoDiPhotons::InitRecoTwoProngInfo(fRecoTwoProngInfo2);
  } else if (cand_pairs_passed >= 2) {
    ExoDiPhotons::FillRecoTwoProngInfo(fRecoTwoProngInfo1, candidates_CHpos[l], candidates_CHneg[l], candidates_center[l],
                                       candidates_photon[l], candidates_Eta[l], candidates_EtaGrommed[l].M(), candidates_passed[l], candidates_matched[l], 
                                       candidates_nearestGenDR[l], candidates_nearestGenIndex[l], candidates_relchargediso[l],
                                       candidates_relneutraliso[l], candidates_relegammaiso[l], candidates_numgamma[l], candidates_nume[l]);
    ExoDiPhotons::FillRecoTwoProngInfo(fRecoTwoProngInfo2, candidates_CHpos[s], candidates_CHneg[s], candidates_center[s],
                                       candidates_photon[s], candidates_Eta[s], candidates_EtaGrommed[s].M(), candidates_passed[s], candidates_matched[s], 
                                       candidates_nearestGenDR[s], candidates_nearestGenIndex[s], candidates_relchargediso[s],
                                       candidates_relneutraliso[s], candidates_relegammaiso[s], candidates_numgamma[s], candidates_nume[s]);
    TLorentzVector phiCandidate = candidates_Eta[l] + candidates_Eta[s];
    fMass_EtaEta = phiCandidate.M();
    TLorentzVector phiCandidate2 = candidates_EtaGrommed[l] + candidates_EtaGrommed[s];
    fMass2_EtaEta = phiCandidate2.M();
    fEtaEtaCand2 = phiCandidate;
    fDR_EtaEta = candidates_Eta[l].DeltaR(candidates_Eta[s]);
    fDPhi_EtaEta = abs(candidates_Eta[l].DeltaPhi(candidates_Eta[s]));
    fDEta_EtaEta = abs(candidates_Eta[l].Eta() - candidates_Eta[s].Eta());
  } // end if block
  // Fill TTree
  fTree2->Fill();
  // Increment cutflow variables
  fCutflow_total += 1;
  if (fNumCHpairs>=1) fCutflow_oneCand += 1;
  if (fNumCHpairs>=2) fCutflow_twoCand += 1;
  if (fNumCHpairsPass>=1) fCutflow_onePass += 1;
  if (fNumCHpairsPass>=2) fCutflow_twoPass += 1;
  if (fNumCHpairsMatch>=1) fCutflow_oneMatch += 1;
  if (fNumCHpairsMatch>=2) fCutflow_twoMatch += 1;
  if (fNumCHpairsPassMatch>=1) fCutflow_onePassMatch += 1;
  if (fNumCHpairsPassMatch>=2) fCutflow_twoPassMatch += 1;
  if (cand_pairs_passed_charged>=1) fCutflow_passCharged +=1;
  if (cand_pairs_passed_neutral>=1) fCutflow_passNeutral +=1;
  if (cand_pairs_passed_egamma>=1) fCutflow_passEGamma +=1;
  if (cand_pairs_passed_photonpt>=1) fCutflow_passPhotonPt +=1;
  if (fDebug) cout << "finished charged decay code" << endl;

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
