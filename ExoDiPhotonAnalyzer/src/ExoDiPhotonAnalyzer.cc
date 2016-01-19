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
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenuFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h" 
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerRecord.h" 
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "L1Trigger/GlobalTrigger/plugins/L1GlobalTrigger.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

// new CommonClasses approach
// these objects are all in the namespace 'ExoDiPhotons'
#include "DiPhotonAnalysis/CommonClasses/interface/RecoPhotonInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/TriggerInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/PhotonID.h"
#include "DiPhotonAnalysis/CommonClasses/interface/EventAndVertexInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/DiphotonInfo.h"


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
  edm::InputTag      fL1InputTag;      // L1 results
  edm::InputTag      fRho25Tag;  
  edm::InputTag      fpileupCollectionTag;         
  edm::LumiReWeighting    LumiWeights;      
  
  bool               fkRemoveSpikes;   // option to remove spikes before filling tree
  bool               fkRequireTightPhotons;  // option to require tight photon id in tree
  bool               fkRequireGenEventInfo;  // generated information for RS graviton files
  bool               fisMC;  //option to decide if MC or Data     
  string             fPUMCFileName;
  string             fPUDataFileName;
  string             fPUDataHistName;
  string             fPUMCHistName;   
  string             fPFIDCategory;   
  string             fIDMethod;   

 
  // tools for clusters
  std::auto_ptr<noZS::EcalClusterLazyTools> lazyTools_;
  edm::InputTag recHitsEBTag_;
  edm::InputTag recHitsEETag_;
  edm::EDGetTokenT<EcalRecHitCollection> recHitsEBToken;
  edm::EDGetTokenT<EcalRecHitCollection> recHitsEEToken;




  // to get L1 info, the L1 guide recommends to make this a member
  // this allows the event setup parts to be cached, rather than refetched every event
  L1GtUtils m_l1GtUtils;

  // my Tree
  TTree *fTree;

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
   
  ExoDiPhotons::diphotonInfo_t fDiphotonInfo;

  // diphoton info based on using hte second or third vtx in event
  ExoDiPhotons::diphotonInfo_t fDiphotonInfoVtx2; 
  ExoDiPhotons::diphotonInfo_t fDiphotonInfoVtx3; 

  TH1F* fpu_n_BeforeCuts; 
  TH1F* fpu_n_BeforeCutsAfterReWeight;
  TH1F *fNumTotalEvents;
  TH1F *fNumTotalWeightedEvents;


  //-----------------taken from Ilya-----------------
  // Format-independent data members
  edm::EDGetTokenT<double> rhoToken_;
  
  // AOD case data members
  edm::EDGetToken photonsToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
  
  // MiniAOD case data members
  edm::EDGetToken photonsMiniAODToken_;
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

  Float_t rho_;      // the rho variable

  // Effective area constants for all isolation types
  EffectiveAreas effAreaChHadrons_;
  EffectiveAreas effAreaNeuHadrons_;
  EffectiveAreas effAreaPhotons_;

  std::vector<Int_t> passLooseId_;
  std::vector<Int_t> passMediumId_;
  std::vector<Int_t> passTightId_;
  //-----------------taken from Ilya-----------------

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
    fL1InputTag(iConfig.getUntrackedParameter<edm::InputTag>("L1Results")),
    fRho25Tag(iConfig.getParameter<edm::InputTag>("rho25Correction")),
    fpileupCollectionTag(iConfig.getUntrackedParameter<edm::InputTag>("pileupCorrection")),
    fkRemoveSpikes(iConfig.getUntrackedParameter<bool>("removeSpikes")),
    fkRequireTightPhotons(iConfig.getUntrackedParameter<bool>("requireTightPhotons")),
    fkRequireGenEventInfo(iConfig.getUntrackedParameter<bool>("requireGenEventInfo")),
    fisMC(iConfig.getUntrackedParameter<bool>("isMC")),
    fPUMCFileName(iConfig.getUntrackedParameter<string>("PUMCFileName")),
    fPUDataFileName(iConfig.getUntrackedParameter<string>("PUDataFileName")), 
    fPUDataHistName(iConfig.getUntrackedParameter<string>("PUDataHistName")),
    fPUMCHistName(iConfig.getUntrackedParameter<string>("PUMCHistName")),
    fPFIDCategory(iConfig.getUntrackedParameter<string>("PFIDCategory")),
    fIDMethod(iConfig.getUntrackedParameter<string>("IDMethod")),
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
  
  genParticlesMiniAODToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticlesMiniAOD"));
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
  fTree->Branch("Diphoton",&fDiphotonInfo,ExoDiPhotons::diphotonInfoBranchDefString.c_str());
  // diphoton info for second or thrid best vertex
  // only bothering to add this for tight-tight tree for now
  fTree->Branch("DiphotonVtx2",&fDiphotonInfoVtx2,ExoDiPhotons::diphotonInfoBranchDefString.c_str());
  fTree->Branch("DiphotonVtx3",&fDiphotonInfoVtx3,ExoDiPhotons::diphotonInfoBranchDefString.c_str());
  

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

  recHitsEBTag_ = iConfig.getUntrackedParameter<edm::InputTag>("RecHitsEBTag",edm::InputTag("reducedEcalRecHitsEB"));
  recHitsEETag_ = iConfig.getUntrackedParameter<edm::InputTag>("RecHitsEETag",edm::InputTag("reducedEcalRecHitsEE"));
  recHitsEBToken = consumes <edm::SortedCollection<EcalRecHit> > (recHitsEBTag_);
  recHitsEEToken = consumes <edm::SortedCollection<EcalRecHit> > (recHitsEETag_);
  // recHitsEBTag_ = iConfig.getUntrackedParameter<edm::InputTag>("RecHitsEBTag",edm::InputTag("reducedEgamma:reducedEBRecHits"));
  // recHitsEETag_ = iConfig.getUntrackedParameter<edm::InputTag>("RecHitsEETag",edm::InputTag("reducedEgamma:reducedEERecHits"));
  // recHitsEBToken = consumes < EcalRecHitCollection > (recHitsEBTag_);
  // recHitsEEToken = consumes < EcalRecHitCollection > (recHitsEETag_);
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


  cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

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
  bool isAOD = true;
  iEvent.getByToken(photonsToken_, photons);
  if( !photons.isValid() ){
    isAOD = false;
    iEvent.getByToken(photonsMiniAODToken_,photons);
  }

  //-----------------taken from Ilya-----------------
  // Get generator level info
  edm::Handle<edm::View<reco::GenParticle> > genParticles;
  if( isAOD )
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
  iEvent.getByLabel("offlinePrimaryVertices",vertexColl);
  //iEvent.getByLabel("offlineSlimmedPrimaryVertices",vertexColl);
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

  //for conversion safe electron veto
  edm::Handle<reco::ConversionCollection> hConversions;
  iEvent.getByLabel("allConversions", hConversions);
  
  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  edm::Handle<reco::GsfElectronCollection> hElectrons;
  iEvent.getByLabel("gedGsfElectrons", hElectrons);
  //edm::Handle<pat::ElectronCollection> hElectrons;
  //iEvent.getByLabel(edm::InputTag("slimmedElectrons"), hElectrons);
  //patElectrons_slimmedElectrons__PAT.obj.embeddedSuperCluster_
  //TO DISENTANGLE BETWEEN MINIAOD AND AOD
  if(!hElectrons.isValid()) {
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
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
   
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


  // get the trig info

  //trig results
  Handle<TriggerResults> hltResultsHandle;
  iEvent.getByLabel(fHltInputTag,hltResultsHandle);
   
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
  fL1TrigInfo.L1_Tech0 = false;
  fL1TrigInfo.L1_Tech36 = false;
  fL1TrigInfo.L1_Tech37 = false;
  fL1TrigInfo.L1_Tech38 = false;
  fL1TrigInfo.L1_Tech39 = false;
  fL1TrigInfo.L1_Tech40 = false;
  fL1TrigInfo.L1_Tech41 = false;
  fL1TrigInfo.L1_Tech42 = false;
  fL1TrigInfo.L1_Tech43 = false;
  fL1TrigInfo.L1_EG2 = false;   
   
  // use the L1GtUtils class, following instructions in 
  // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideL1TriggerL1GtUtils
   
  // before accessing any result from L1GtUtils, one must retrieve and cache
  // the L1 trigger event setup
  m_l1GtUtils.retrieveL1EventSetup(iSetup);

  // access L1 trigger results using public methods from L1GtUtils
  // always check on error code returned by that method
  int iErrorCode = -1;
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
  

   
  fL1TrigInfo.L1_Tech0 =    m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BPTX_plus_AND_minus.v0",iErrorCode);
  fL1TrigInfo.L1_Tech36 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_halo_beam2_inner.v0",iErrorCode);
  fL1TrigInfo.L1_Tech37 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_halo_beam2_outer.v0",iErrorCode);
  fL1TrigInfo.L1_Tech38 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_halo_beam1_inner.v0",iErrorCode);
  fL1TrigInfo.L1_Tech39 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_halo_beam1_outer.v0",iErrorCode);
  fL1TrigInfo.L1_Tech40 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_minBias_threshold1.v0",iErrorCode);
  fL1TrigInfo.L1_Tech41 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_minBias_threshold2.v0",iErrorCode);
  fL1TrigInfo.L1_Tech42 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_splash_beam1.v0",iErrorCode);
  fL1TrigInfo.L1_Tech43 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_splash_beam2.v0",iErrorCode);


  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  // ecal information
  lazyTools_ = std::auto_ptr<noZS::EcalClusterLazyTools>( new noZS::EcalClusterLazyTools(iEvent, iSetup, recHitsEBToken, recHitsEEToken));   

  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  // get ecal barrel recHits for spike rejection
  edm::Handle<EcalRecHitCollection> recHitsEB_h;
  iEvent.getByLabel(edm::InputTag("reducedEcalRecHitsEB"), recHitsEB_h );
  // iEvent.getByLabel(edm::InputTag("reducedEgamma:reducedEBRecHits"), recHitsEB_h );
  const EcalRecHitCollection * recHitsEB = 0;
  if ( ! recHitsEB_h.isValid() ) {
    LogError("ExoDiPhotonAnalyzer") << " ECAL Barrel RecHit Collection not available !"; return;
  } else {
    recHitsEB = recHitsEB_h.product();
  }
  //TO DISENTANGLE BETWEEN MINIAOD AND AOD

  //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  edm::Handle<EcalRecHitCollection> recHitsEE_h;
  iEvent.getByLabel(edm::InputTag("reducedEcalRecHitsEE"), recHitsEE_h );
  //iEvent.getByLabel(edm::InputTag("reducedEgamma:reducedEERecHits"), recHitsEE_h );
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
  Handle<reco::PhotonCollection> photonColl;
  iEvent.getByLabel(fPhotonTag,photonColl);
  // edm::Handle<pat::PhotonCollection> photonColl;
  // iEvent.getByLabel("slimmedPhotons",photonColl);

  // If photon collection is empty, exit
  if (!photonColl.isValid()) {
    cout << "No Photons! Move along, there's nothing to see here .." <<endl;
    return;
  }
  const reco::PhotonCollection *myPhotonColl = photonColl.product();
  cout<<"photoncoll size "<<myPhotonColl->size()<<endl;
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

  // photon loop
  for(reco::PhotonCollection::const_iterator recoPhoton = photonColl->begin(); recoPhoton!=photonColl->end(); recoPhoton++) {
    
    phoIndex++;

    //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

    const reco::Photon testPho = *recoPhoton;
    edm::Ptr<reco::Photon> testPhoPtr(photonColl,phoIndex);

    //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

    //-----------------taken from Ilya-----------------
    float full5x5sigmaIetaIeta = (*full5x5SigmaIEtaIEtaMap)[testPhoPtr];
    float chIso =  (*phoChargedIsolationMap)[testPhoPtr];
    float nhIso =  (*phoNeutralHadronIsolationMap)[testPhoPtr];
    float phIso = (*phoPhotonIsolationMap)[testPhoPtr];
    
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
    //and we also have to test the conversion safe electron veto
    bool passelecveto = !ConversionTools::hasMatchedPromptElectron(recoPhoton->superCluster(), hElectrons, hConversions, beamSpot.position());
    if (passelecveto) cout << "Passed electron veto!" << endl;
    else cout << "Failed electron veto!" << endl;
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
    ExoDiPhotons::recoPhotonInfo_t tempInfo;
    ExoDiPhotons::FillRecoPhotonInfo(tempInfo,&(*recoPhoton),lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);
    Bool_t isSaturated = tempInfo.isSaturated;

    if(recoPhoton->pt() < fMin_pt) cout << "photon pt = " << recoPhoton->pt() << " which is less than " << fMin_pt << "!" << endl;
    if(ExoDiPhotons::isGapPhoton(&(*recoPhoton))) cout << "this is a photon in the gap!" << endl;
    if(ExoDiPhotons::isASpike(&(*recoPhoton))) cout << "this photon is a spike!" << endl;

    if(recoPhoton->pt() < fMin_pt) continue;
    if(ExoDiPhotons::isGapPhoton(&(*recoPhoton))) continue;
    if(ExoDiPhotons::isASpike(&(*recoPhoton))) continue;

    cout << "pt, gap, and spike cuts passed!" << endl;
    
    //Now we choose which ID to use (PF or Det)
    if(MethodID.Contains("highpt")){
      if(ExoDiPhotons::isPFTightPhoton(&(*recoPhoton),rhocorPFIsoCH,rhocorPFIsoNH,rhocorPFIsoPH,full5x5sigmaIetaIeta,passelecveto,MethodID,CategoryPFID,isSaturated)){
        selectedPhotons.push_back(*recoPhoton);
        cout << "photon passed the high pt id!" << endl;
      }
      else cout << "photon failed the high pt id!" << endl;
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
  cout << "photon pased fakeable object cut" << endl;
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

  // now count many candidate photons we have in this event
  //   cout << "N candidate photons = " << selectedPhotons.size() <<endl;
  fNTightPhotons = selectedPhotons.size();
  if(fNTightPhotons >= 2) cout<<"great we have two tight photons"<<endl;
  else cout << "only " << fNTightPhotons << "photons passed the cuts!" << endl;

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
  cout << "allTightOrFakeableObjects.size(): " << allTightOrFakeableObjects.size() << endl;
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

    cout<<"here testing whether i get the correct pointer to my photon TorF"<<endl;
    reco::Photon* TorFObject1 = &(allTightOrFakeableObjects[0].first);
    reco::Photon* TorFObject2 = &(allTightOrFakeableObjects[1].first);
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

    int indexTorFObject1 = -1;
    int indexTorFObject2 = -1;
    int myIndex = -1;

    for(reco::PhotonCollection::const_iterator recoPhoton = photonColl->begin(); recoPhoton!=photonColl->end(); recoPhoton++) {
      myIndex++;
      //const reco::Photon* myPhoton = &(*recoPhoton);
      //if(myPhoton == TorFObject1) {
      if(recoPhoton->pt() == TorFObject1->pt()) {
	cout<<"Great, I've found torf object 1 "<<endl;
	indexTorFObject1 = myIndex;
      }
      //if(myPhoton == TorFObject2) {
      if(recoPhoton->pt() == TorFObject2->pt()) {
	cout<<"Great, I've found torf object 2 "<<endl;
	indexTorFObject2 = myIndex;
      }
    }

    //in principle the indices should be always greater than 1
    //because we necessarily found the two objects in the 
    //photon collection


    edm::Ptr<reco::Photon> TorFObject1Ptr(photonColl,indexTorFObject1);
    cout<<TorFObject1Ptr->energy()<<" "
	<<TorFObject1Ptr->eta()<<" "
	<<TorFObject1Ptr->et()<<" "
	<<TorFObject1Ptr->phi()<<" "
	<<endl;

    edm::Ptr<reco::Photon> TorFObject2Ptr(photonColl,indexTorFObject2);
    cout<<TorFObject2Ptr->energy()<<" "
	<<TorFObject2Ptr->eta()<<" "
	<<TorFObject2Ptr->et()<<" "
	<<TorFObject2Ptr->phi()<<" "
	<<endl;
    cout<<"here testing whether i get the correct pointer to my photon TorF"<<endl;

    ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo1,&allTightOrFakeableObjects[0].first,lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);
    fRecoPhotonInfo1.isFakeable = allTightOrFakeableObjects[0].second;
    fRecoPhotonInfo1.hasMatchedPromptElec = ConversionTools::hasMatchedPromptElectron((&allTightOrFakeableObjects[0].first)->superCluster(), hElectrons, hConversions, beamSpot.position());
    
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
    fRecoPhotonInfo2.hasMatchedPromptElec = ConversionTools::hasMatchedPromptElectron((&allTightOrFakeableObjects[1].first)->superCluster(), hElectrons, hConversions, beamSpot.position());

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

      cout<<"here testing whether i get the correct pointer to my photon T"<<endl;
      reco::Photon* TObject1 = &(selectedPhotons[0]);
      reco::Photon* TObject2 = &(selectedPhotons[1]);
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

      int indexTObject1 = -1;
      int indexTObject2 = -1;
      myIndex = -1;

      for(reco::PhotonCollection::const_iterator recoPhoton = photonColl->begin(); recoPhoton!=photonColl->end(); recoPhoton++) {
	myIndex++;
	//const reco::Photon* myPhoton = &(*recoPhoton);
	if(recoPhoton->pt() == TObject1->pt()) {
	  //if(myPhoton == TObject1) {
	  cout<<"Great, I've found t object 1 "<<endl;
	  indexTObject1 = myIndex;
	}
	if(recoPhoton->pt() == TObject2->pt()) {
	  //if(myPhoton == TObject2) {
	  cout<<"Great, I've found t object 2 "<<endl;
	  indexTObject2 = myIndex;
	}
      }

      //in principle the indices should be always greater than 1
      //because we necessarily found the two objects in the 
      //photon collection


      edm::Ptr<reco::Photon> TObject1Ptr(photonColl,indexTObject1);
      cout<<TObject1Ptr->energy()<<" "
	  <<TObject1Ptr->eta()<<" "
	  <<TObject1Ptr->et()<<" "
	  <<TObject1Ptr->phi()<<" "
	  <<endl;

      edm::Ptr<reco::Photon> TObject2Ptr(photonColl,indexTObject2);
      cout<<TObject2Ptr->energy()<<" "
	  <<TObject2Ptr->eta()<<" "
	  <<TObject2Ptr->et()<<" "
	  <<TObject2Ptr->phi()<<" "
	  <<endl;
      cout<<"here testing whether i get the correct pointer to my photon T"<<endl;

      // must specifically declare isFakeable status (should be Tight = not True = false                       
      ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo1,&selectedPhotons[0],lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);
      fRecoPhotonInfo1.hasMatchedPromptElec = ConversionTools::hasMatchedPromptElectron((&selectedPhotons[0])->superCluster(), hElectrons, hConversions, beamSpot.position());

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
      fRecoPhotonInfo2.hasMatchedPromptElec = ConversionTools::hasMatchedPromptElectron((&selectedPhotons[1])->superCluster(), hElectrons, hConversions, beamSpot.position());

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


  




#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
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
ExoDiPhotonAnalyzer::endJob() {
}



//define this as a plug-in
DEFINE_FWK_MODULE(ExoDiPhotonAnalyzer);
