// -*- C++ -*-
//
// Package:    ExoDiPhotonSignalMCAnalyzer
// Class:      ExoDiPhotonSignalMCAnalyzer
// 
/**\class ExoDiPhotonSignalMCAnalyzer ExoDiPhotonSignalMCAnalyzer.cc DiPhotonAnalysis/ExoDiPhotonSignalMCAnalyzer/src/ExoDiPhotonSignalMCAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Conor Henderson,40 1-B01,+41227671674,
//         Created:  Wed Jun 16 17:06:28 CEST 2010
// $Id: ExoDiPhotonSignalMCAnalyzer.cc,v 1.10 2013/02/12 14:01:15 scooper Exp $
//
//


// system include files
#include <memory>
#include <iomanip>

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

// geometry
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
//#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
//#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
//#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
//#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"

//for photons
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TMath.h"

#include "DataFormats/Math/interface/deltaR.h"

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

// for MC
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// new CommonClasses approach
// these objects are all in the namespace 'ExoDiPhotons'
#include "DiPhotonAnalysis/CommonClasses/interface/RecoPhotonInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/MCTrueObjectInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/TriggerInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/PhotonID.h"
#include "DiPhotonAnalysis/CommonClasses/interface/EventAndVertexInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/DiphotonInfo.h"
//new for PF ID definition
#include "DiPhotonAnalysis/CommonClasses/interface/PFPhotonID.h"

//new for PU gen
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// new for LumiReweighting
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

//for conversion safe electron veto
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//new for PFIsolation code
//#include "EGamma/EGammaAnalysisTools/src/PFIsolationEstimator.cc"

//-----------------taken from Ilya-----------------
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
//-----------------taken from Ilya-----------------

using namespace std;

//
// class declaration
//


class ExoDiPhotonSignalMCAnalyzer : public edm::EDAnalyzer {
   public:
      explicit ExoDiPhotonSignalMCAnalyzer(const edm::ParameterSet&);
      ~ExoDiPhotonSignalMCAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

  // ----------member data ---------------------------
  
  // input tags and parameters
  edm::InputTag           fPhotonTag;       // select photon collection 
  double                  fMin_pt;          // min pt cut (photons)
  edm::InputTag           fHltInputTag;     // hltResults
  edm::InputTag           fRho25Tag;
  edm::InputTag           pileupCollectionTag;
  edm::LumiReWeighting    LumiWeights;
  
  bool                    fkRemoveSpikes;         // option to remove spikes before filling tree
  bool                    fkRequireGenEventInfo;  // generated information for RS graviton files
  string                  PUMCFileName;
  string                  PUDataFileName;
  string                  PUDataHistName;
  string                  PUMCHistName;
  string                  fPFIDCategory;
  string                  fIDMethod;
  
  // tools for clusters
  //std::auto_ptr<EcalClusterLazyTools> lazyTools_;
  std::auto_ptr<noZS::EcalClusterLazyTools> lazyTools_;
  edm::InputTag recHitsEBTag_;
  edm::InputTag recHitsEETag_;
  edm::EDGetTokenT<EcalRecHitCollection> recHitsEBToken;
  edm::EDGetTokenT<EcalRecHitCollection> recHitsEEToken;

  // my Tree
  TTree *fTree;
  
  ExoDiPhotons::eventInfo_t fEventInfo;
  ExoDiPhotons::vtxInfo_t fVtxInfo;
  ExoDiPhotons::beamSpotInfo_t fBeamSpotInfo;
  
  ExoDiPhotons::hltTrigInfo_t fHLTInfo;
  
  ExoDiPhotons::mcTrueObjectInfo_t fUnstablePhoton1Info; // leading photon directly from RSG decay
  ExoDiPhotons::mcTrueObjectInfo_t fUnstablePhoton2Info;

  ExoDiPhotons::mcTrueObjectInfo_t fStablePhoton1Info; // leading final state photon coming from RSG decay (may pair produce)
  ExoDiPhotons::mcTrueObjectInfo_t fStablePhoton2Info;

  ExoDiPhotons::recoPhotonInfo_t fRecoPhoton1Info; // leading matched reco photon 
  ExoDiPhotons::recoPhotonInfo_t fRecoPhoton2Info; // subleading photon
  
  ExoDiPhotons::diphotonInfo_t fDiphotonGenUnstableInfo;
  ExoDiPhotons::diphotonInfo_t fDiphotonGenStableInfo;
  ExoDiPhotons::diphotonInfo_t fDiphotonRecoInfo;
 
  // Store PileUp Info
  double fRho25;
  int fpu_n;
  int fBC;
  double fMCPUWeight;

  // Events before selectoin
  TH1F* fpu_n_BeforeCuts;
  TH1F* fpu_n_BeforeCutsAfterReWeight;

  // SIC add
  // for PFIsolation Code
  //PFIsolationEstimator isolator04;
  //PFIsolationEstimator isolator03;
  //PFIsolationEstimator isolator02;

  //-----------------taken from Ilya-----------------
  // Format-independent data members
  edm::EDGetTokenT<double> rhoToken_;
  
  // AOD case data members
  edm::EDGetToken photonsToken_;
  //edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
  
  // MiniAOD case data members
  edm::EDGetToken photonsMiniAODToken_;
  //edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesMiniAODToken_;
  
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
ExoDiPhotonSignalMCAnalyzer::ExoDiPhotonSignalMCAnalyzer(const edm::ParameterSet& iConfig)
  : fPhotonTag(iConfig.getUntrackedParameter<edm::InputTag>("photonCollection")),
    fMin_pt(iConfig.getUntrackedParameter<double>("ptMin")),
    // note that the HLT process name can vary for different MC samples
    // so be sure to adjsut correctly in cfg
    fHltInputTag(iConfig.getUntrackedParameter<edm::InputTag>("hltResults")),
    fRho25Tag(iConfig.getParameter<edm::InputTag>("rho25Correction")),
    pileupCollectionTag(iConfig.getUntrackedParameter<edm::InputTag>("pileupCorrection")),
    fkRemoveSpikes(iConfig.getUntrackedParameter<bool>("removeSpikes")),
    fkRequireGenEventInfo(iConfig.getUntrackedParameter<bool>("requireGenEventInfo")),
    PUMCFileName(iConfig.getUntrackedParameter<string>("PUMCFileName")),
    PUDataFileName(iConfig.getUntrackedParameter<string>("PUDataFileName")),
    PUDataHistName(iConfig.getUntrackedParameter<string>("PUDataHistName")),
    PUMCHistName(iConfig.getUntrackedParameter<string>("PUMCHistName")),
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

  std::cout << "ExoDiPhotonAnalyzer: ID Method used " << fIDMethod.c_str() << ", "
	    << "PF ID Category " << fPFIDCategory.c_str()
	    << std::endl;

  //-----------------taken from Ilya-----------------
  //
  // Prepare tokens for all input collections and objects
  //
  // AOD tokens
  photonsToken_ = mayConsume<edm::View<reco::Photon> >
    (iConfig.getParameter<edm::InputTag>
     ("photons"));
  
  //  genParticlesToken_ = mayConsume<edm::View<reco::GenParticle> >
  //    (iConfig.getParameter<edm::InputTag>
  //     ("genParticles"));
  
  // MiniAOD tokens
  photonsMiniAODToken_ = mayConsume<edm::View<reco::Photon> >
    (iConfig.getParameter<edm::InputTag>
     ("photonsMiniAOD"));
  
  //  genParticlesMiniAODToken_ = mayConsume<edm::View<reco::GenParticle> >
  //    (iConfig.getParameter<edm::InputTag>
  //     ("genParticlesMiniAOD"));
  //-----------------taken from Ilya-----------------

  edm::Service<TFileService> fs;
  fTree = fs->make<TTree>("fTree","PhotonTree");
  fpu_n_BeforeCuts = fs->make<TH1F>("fpu_n_BeforeCuts","PileUpBeforeCuts",300,0,300);
  fpu_n_BeforeCutsAfterReWeight = fs->make<TH1F>("fpu_n_BeforeCutsAfterReWeight","PileUpBeforeCuts",300,0,300);
 
  // now with CommonClasses, use fthe string defined in the header

  fTree->Branch("Event",&fEventInfo,ExoDiPhotons::eventInfoBranchDefString.c_str());
  fTree->Branch("Vtx",&fVtxInfo,ExoDiPhotons::vtxInfoBranchDefString.c_str());
  fTree->Branch("BeamSpot",&fBeamSpotInfo,ExoDiPhotons::beamSpotInfoBranchDefString.c_str());
  fTree->Branch("TrigHLT",&fHLTInfo,ExoDiPhotons::hltTrigBranchDefString.c_str());

  fTree->Branch("UnstablePhoton1",&fUnstablePhoton1Info,ExoDiPhotons::mcTrueObjectInfoBranchDefString.c_str());
  fTree->Branch("UnstablePhoton2",&fUnstablePhoton2Info,ExoDiPhotons::mcTrueObjectInfoBranchDefString.c_str());

  fTree->Branch("StablePhoton1",&fStablePhoton1Info,ExoDiPhotons::mcTrueObjectInfoBranchDefString.c_str());
  fTree->Branch("StablePhoton2",&fStablePhoton2Info,ExoDiPhotons::mcTrueObjectInfoBranchDefString.c_str());

  fTree->Branch("Photon1",&fRecoPhoton1Info,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  fTree->Branch("Photon2",&fRecoPhoton2Info,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  
  // signal diphoton info? eg to probe true MC width?
  // reco diphoton info?
  
  fTree->Branch("DiphotonGenUnstable",&fDiphotonGenUnstableInfo,ExoDiPhotons::diphotonInfoBranchDefString.c_str());
  fTree->Branch("DiphotonGenStable",&fDiphotonGenStableInfo,ExoDiPhotons::diphotonInfoBranchDefString.c_str());
  fTree->Branch("Diphoton",&fDiphotonRecoInfo,ExoDiPhotons::diphotonInfoBranchDefString.c_str());
   
  
  //Pileup Info
  fTree->Branch("rho25",&fRho25,"rho25/D");
  fTree->Branch("pu_n", &fpu_n, "pu_n/I");
  fTree->Branch("MCPUWeight",&fMCPUWeight,"MCPUWeight/D");
  

  // SIC add
  //new PFIsolation code
  //isolator04.initializePhotonIsolation(kTRUE);
  //isolator04.setConeSize(0.4);
  //isolator03.initializePhotonIsolation(kTRUE);
  //isolator03.setConeSize(0.3);
  //isolator02.initializePhotonIsolation(kTRUE);
  //isolator02.setConeSize(0.2);
   
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


ExoDiPhotonSignalMCAnalyzer::~ExoDiPhotonSignalMCAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//




// ------------ method called to for each event  ------------
void
ExoDiPhotonSignalMCAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;

   //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;


   //cout<<"event"<<endl;
   //initialize
   fpu_n = -99999.99;
   fBC = -99999.99;
   fMCPUWeight = -99999.99;

   // basic event info
   ExoDiPhotons::FillEventInfo(fEventInfo,iEvent);

   //-----------------taken from Ilya-----------------
  // Retrieve the collection of photons from the event.
  // If we fail to retrieve the collection with the standard AOD
  // name, we next look for the one with the stndard miniAOD name. 
  //   We use exactly the same handle for AOD and miniAOD formats
  // since pat::Photon objects can be recast as reco::Photon objects.
  edm::Handle<edm::View<reco::Photon> > photons;
  //bool isAOD = true;
  iEvent.getByToken(photonsToken_, photons);
  if( !photons.isValid() ){
    //isAOD = false;
    iEvent.getByToken(photonsMiniAODToken_,photons);
  }
  //-----------------taken from Ilya-----------------
  // Get generator level info
  // edm::Handle<edm::View<reco::GenParticle> > genParticles;
  // if( isAOD )
  //  iEvent.getByToken(genParticlesToken_,genParticles);
  //  else
  //  iEvent.getByToken(genParticlesMiniAODToken_,genParticles);
  
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

  //fNumTotalEvents->Fill(1.);
  //fNumTotalWeightedEvents->Fill(1.,fEventInfo.weight);

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
  //   cout << "N vertices = " << vertexColl->size() <<endl;
  // fVtxInfo.Nvtx = vertexColl->size();

  fVtxInfo.vx = -99999.99;
  fVtxInfo.vy = -99999.99;
  fVtxInfo.vz = -99999.99;
  fVtxInfo.isFake = -99;   
  fVtxInfo.Ntracks = -99;
  fVtxInfo.sumPtTracks = -99999.99;
  fVtxInfo.ndof = -99999.99;
  fVtxInfo.d0 = -99999.99;
  
  //get the reference to 1st vertex for use in fGetIsolation
  //for PFIsolation calculation
  reco::VertexRef firstVtx(vertexColl,0);

  const reco::Vertex *vertex1 = NULL; // best vertex (by trk SumPt)
  // note for higher lumi, may want to also store second vertex, for pileup studies
  
  // even if the vertex has sumpt=0, its still enough to be the 'highest'
  double highestSumPtTracks1 = -1.0; 
  
  for(reco::VertexCollection::const_iterator vtx=vertexColl->begin(); vtx!=vertexColl->end(); vtx++) {
    
    // loop over assoc tracks to get sum pt
    //     fVtxInfo.sumPtTracks = 0.0;
    double sumPtTracks = 0.0;
    
    for(reco::Vertex::trackRef_iterator vtxTracks=vtx->tracks_begin(); vtxTracks!=vtx->tracks_end();vtxTracks++) {
      //       fVtxInfo.sumPtTracks += (**vtxTracks).pt();
      sumPtTracks += (**vtxTracks).pt();
    }
    
    //     cout << "Vtx x = "<< vtx->x()<<", y= "<< vtx->y()<<", z = " << vtx->z() << ";  N tracks = " << vtx->tracksSize() << "; isFake = " << vtx->isFake() <<", sumPt(tracks) = "<<sumPtTracks << "; ndof = " << vtx->ndof()<< "; d0 = " << vtx->position().rho() << endl;
     
    // and note that this vertex collection can contain vertices with Ntracks = 0
    // watch out for these!

    if(sumPtTracks > highestSumPtTracks1) {
      // then this new best
      vertex1 = &(*vtx);
      highestSumPtTracks1 = sumPtTracks;      
    }

  }// end vertex loop
  
  if(vertex1) {
    ExoDiPhotons::FillVertexInfo(fVtxInfo,vertex1);
    // fill the SumPt Tracks separately for now
    fVtxInfo.sumPtTracks = highestSumPtTracks1;  
  }

   // beam spot
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

   if(beamSpotHandle.isValid()) {
     beamSpot = *beamSpotHandle;
     ExoDiPhotons::FillBeamSpotInfo(fBeamSpotInfo,beamSpot);
   }
   
  //for PFIsolation code
  //Handle<PFCandidateCollection> pfCandidatesColl;
  //iEvent.getByLabel("particleFlow",pfCandidatesColl);
  //const PFCandidateCollection * pfCandidates = pfCandidatesColl.product();

  //for conversion safe electron veto
  edm::Handle<reco::ConversionCollection> hConversions;
  iEvent.getByLabel("allConversions", hConversions);
  
  edm::Handle<reco::GsfElectronCollection> hElectrons;
  //iEvent.getByLabel("gsfElectrons", hElectrons);
  iEvent.getByLabel("gedGsfElectrons", hElectrons);
  //edm::Handle<pat::ElectronCollection> hElectrons;
  //iEvent.getByLabel(edm::InputTag("slimmedElectrons"), hElectrons);
  //patElectrons_slimmedElectrons__PAT.obj.embeddedSuperCluster_
  //TO DISENTANGLE BETWEEN MINIAOD AND AOD
  if(!hElectrons.isValid()) {
    cout<<"no ged gsf electrons "<<endl;
    return;
  }


  // L1 info

  // HLT info
  Handle<TriggerResults> hltResultsHandle;
  iEvent.getByLabel(fHltInputTag,hltResultsHandle);
  
  if(!hltResultsHandle.isValid()) {
    cout << "HLT results not valid!" <<endl;
    cout << "Couldnt find TriggerResults with input tag " << fHltInputTag << endl;
    return;
  }

  const TriggerResults *hltResults = hltResultsHandle.product();
  //   cout << *hltResults <<endl;
  // this way of getting triggerNames should work, even when the 
  // trigger menu changes from one run to the next
  // alternatively, one could also use the HLTConfigPovider
  const TriggerNames & hltNames = iEvent.triggerNames(*hltResults);

  // now we just use the FillHLTInfo() function from TrigInfo.h:
  ExoDiPhotons::FillHLTInfo(fHLTInfo,hltResults,hltNames);
  
  //Add PileUp Information
  edm::Handle<std::vector<PileupSummaryInfo> > pileupHandle;
  iEvent.getByLabel(pileupCollectionTag, pileupHandle);
  std::vector<PileupSummaryInfo>::const_iterator PUI;
  
  if (pileupHandle.isValid()){
    for (PUI = pileupHandle->begin();PUI != pileupHandle->end(); ++PUI){
      fBC = PUI->getBunchCrossing() ;
      if(fBC==0){
	//Select only the in time bunch crossing with bunch crossing=0
	fpu_n = PUI->getTrueNumInteractions();
	fpu_n_BeforeCuts->Fill(fpu_n);
      }
    }
    fMCPUWeight = LumiWeights.weight(fpu_n);
    fpu_n_BeforeCutsAfterReWeight->Fill(fpu_n,fMCPUWeight);
  }
  
  //add rho correction
  //      double rho;
  // edm::Handle<double> rho25Handle;
  // iEvent.getByLabel(fRho25Tag, rho25Handle);

  //if (!rho25Handle.isValid()){
  //cout<<"rho25 not found"<<endl;
  //return;
  // }
  
  //fRho25 = *(rho25Handle.product());
  fRho25 = *rhoH;

  // ecal information
  // lazyTools_ = std::auto_ptr<EcalClusterLazyTools>( new  EcalClusterLazyTools(iEvent,iSetup,edm::InputTag("reducedEcalRecHitsEB"),edm::InputTag("reducedEcalRecHitsEE")));
  // //lazyTools_ = std::auto_ptr<EcalClusterLazyTools>( new EcalClusterLazyTools(iEvent,iSetup,edm::InputTag("ecalRecHit:EcalRecHitsEB"),edm::InputTag("ecalRecHit:EcalRecHitsEE")) );
  lazyTools_ = std::auto_ptr<noZS::EcalClusterLazyTools>( new noZS::EcalClusterLazyTools(iEvent, iSetup, recHitsEBToken, recHitsEEToken));

   // get ecal barrel recHits for spike rejection
   edm::Handle<EcalRecHitCollection> recHitsEB_h;
   iEvent.getByLabel(edm::InputTag("reducedEcalRecHitsEB"), recHitsEB_h );
   const EcalRecHitCollection * recHitsEB = 0;
   if ( ! recHitsEB_h.isValid() ) {
     LogError("ExoDiPhotonAnalyzer") << " ECAL Barrel RecHit Collection not available !"; return;
   } else {
     recHitsEB = recHitsEB_h.product();
   }
   //TO DISENTANGLE BETWEEN MINIAOD AND AOD

   edm::Handle<EcalRecHitCollection> recHitsEE_h;
   iEvent.getByLabel(edm::InputTag("reducedEcalRecHitsEE"), recHitsEE_h );
   const EcalRecHitCollection * recHitsEE = 0;
   if ( ! recHitsEE_h.isValid() ) {
    LogError("ExoDiPhotonAnalyzer") << " ECAL Endcap RecHit Collection not available !"; return;
   } else {
    recHitsEE = recHitsEE_h.product();
   }
   //TO DISENTANGLE BETWEEN MINIAOD AND AOD

   edm::ESHandle<EcalChannelStatus> chStatus;
   iSetup.get<EcalChannelStatusRcd>().get(chStatus);
   const EcalChannelStatus *ch_status = chStatus.product(); 
   
   // get the photon collection
   Handle<reco::PhotonCollection> photonColl;
   iEvent.getByLabel(fPhotonTag,photonColl);

   // If photon collection is empty, exit
   if (!photonColl.isValid()) {
     cout << "No Photons! Move along, there's nothing to see here .." <<endl;
     return;
   }

   //   const reco::PhotonCollection *myPhotonColl = photonColl.product();
   //cout<<"photoncoll size "<<myPhotonColl->size()<<endl;
   //   cout << "N reco photons = " << photonColl->size() <<endl;


   //cout << "N photons = " << photonColl->size() <<endl;

   //      //   photon loop
   //    for(reco::PhotonCollection::const_iterator recoPhoton = photonColl->begin(); recoPhoton!=photonColl->end(); recoPhoton++) {
   //
   //   //  cout << "Reco photon et, eta, phi = " << recoPhoton->et() <<", "<<recoPhoton->eta()<< ", "<< recoPhoton->phi();
   //   //  cout << "; eMax/e3x3 = " << recoPhoton->maxEnergyXtal()/recoPhoton->e3x3();
   //   //   cout << "; hadOverEm = " << recoPhoton->hadronicOverEm();
   //   //   cout << "; trkIso = " << recoPhoton->trkSumPtHollowConeDR04();
   //   //   cout << "; ecalIso = " << recoPhoton->ecalRecHitSumEtConeDR04();
   //   //   cout << "; hcalIso = " << recoPhoton->hcalTowerSumEtConeDR04();
   //   //   cout << "; pixelSeed = " << recoPhoton->hasPixelSeed();
   ////      cout << endl;
   //
   //
   //    } // end reco photon loop
   
   // need these??
   TString CategoryPFID(fPFIDCategory.c_str());
   TString MethodID(fIDMethod.c_str());
   
   //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

   /*        
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

     //we retrieve the effective areas
     //Remember effareaCH = 1st, effareaNH = 2nd, effareaPH = 3rd
     //std::vector<double> effareas = ExoDiPhotons::EffectiveAreas(&(*recoPhoton));
     // double rhocorPFIsoCH = isoChargedHadronsWithEA;
     //double rhocorPFIsoNH = isoNeutralHadronsWithEA;
     //double rhocorPFIsoPH = isoPhotonsWithEA;
     //double pfisoall = rhocorPFIsoCH + rhocorPFIsoNH + rhocorPFIsoPH;
     //and we also have to test the conversion safe electron veto
     //bool passelecveto = !ConversionTools::hasMatchedPromptElectron(recoPhoton->superCluster(), hElectrons, hConversions, beamSpot.position());
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

   
     //Now we choose which ID to use (PF or Det)
     if(MethodID.Contains("Detector")){
       if(ExoDiPhotons::isTightPhoton(&(*recoPhoton),fRho25) && !ExoDiPhotons::isGapPhoton(&(*recoPhoton)) && !ExoDiPhotons::isASpike(&(*recoPhoton))  ) {
	 //	    if( !ExoDiPhotons::isASpike(&(*recoPhoton))  ) {   
	 selectedPhotons.push_back(*recoPhoton);
       }
     }
     else if(MethodID.Contains("ParticleFlow")){
       //CAREFUL UNCOMMENT THAT WHEN DONE
       // if(!(ExoDiPhotons::isPFTightPhoton(&(*recoPhoton),rhocorPFIsoCH,rhocorPFIsoNH,rhocorPFIsoPH,full5x5sigmaIetaIeta,passelecveto,CategoryPFID))) continue;
       if(!isPassLoose) continue;
       if(ExoDiPhotons::isGapPhoton(&(*recoPhoton))) continue;
       if(ExoDiPhotons::isASpike(&(*recoPhoton))) continue;
       selectedPhotons.push_back(*recoPhoton);
     }

   
   */

   // =============================================
   // UNSTABLE AND STABLE PHOTON IDENTIFICATION
   // Find RSG in GenParticleCollection, and follow
   // its photon daughters.
   // =============================================

   // get GenParticleCollection   
   Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByLabel("genParticles",genParticles);

   if(!genParticles.isValid()) {
     cout << "No Gen Particles collection!" << endl;
     return;
   }

   ExoDiPhotons::InitMCTrueObjectInfo(fUnstablePhoton1Info);
   ExoDiPhotons::InitMCTrueObjectInfo(fUnstablePhoton2Info);

   ExoDiPhotons::InitMCTrueObjectInfo(fStablePhoton1Info);
   ExoDiPhotons::InitMCTrueObjectInfo(fStablePhoton2Info);

   const reco::GenParticle *unstablePhoton1 = NULL;
   const reco::GenParticle *unstablePhoton2 = NULL;

   const reco::GenParticle *stablePhoton1 = NULL;
   const reco::GenParticle *stablePhoton2 = NULL;

   using namespace reco;

   for(unsigned i = 0; i < genParticles->size(); ++ i) {
     
     const GenParticle & p = (*genParticles)[i];
     int id = p.pdgId();
     int st = p.status();  
     
     //const Candidate * mom = p.mother();
     //double pt = p.pt(), eta = p.eta(), phi = p.phi(), mass = p.mass();
     //double vx = p.vx(), vy = p.vy(), vz = p.vz();
     //int charge = p.charge();

     unsigned nDau = p.numberOfDaughters();
     
     if(id == 5100039 && p.isLastCopy()){     
       cout << "---RSG---" << endl;
       cout << "numberOfDaughters: " << nDau << endl;
       cout << "Status: " << st << endl;
       cout << "isLastCopy(): " << p.isLastCopy() << endl;
       cout << "GenParticle Pt, Eta, Phi: " << p.pt() << ", " << p.eta() << ", " << p.phi() << endl;
       for(unsigned j = 0; j < nDau; ++ j) {
	 
	 //const Candidate * dau = p.daughter( j );
	 const GenParticle * dau = (GenParticle *) p.daughter( j );
	 unsigned nDauDau = dau->numberOfDaughters();
	 int dauId = dau->pdgId();
	 
	 // since p is last copy, should be true, but will require it anyway
	 if(nDau == 2 && dauId == 22) {
	   //const GenParticle  *pho = dau;
	   cout << " -dauID: " << dauId << endl;
	   cout << "   isLastCopy: " << dau->isLastCopy() << endl;
	   cout << "   numberOfDaughters: " << nDauDau << endl;
	   cout << "   dauStatus: " << dau->status() << endl;
	   cout << "   dauPt, Eta, Phi: " << dau->pt() << ", " << dau->eta() << ", " << dau->phi() << endl;
	   
	   // assign unstable photons
	   if (!unstablePhoton1) unstablePhoton1 = dau;
	   else unstablePhoton2 = dau;

	   // assign stable photons
	   // CASE 1: Unstable photon is stable photon
	   if (nDauDau == 0 && dau->status() == 1) {
	     if (!stablePhoton1) stablePhoton1 = dau;
	     else stablePhoton2 = dau;
	   }
	   
	   // assign stable photons
	   // CASE 2: One unstable photon pair produces, while the other has a stable photon daughter,
	   // or they both pair produce
	   else if (nDauDau == 1 && dau->daughter(0)->pdgId() == 22) {
	     if (dau->daughter(0)->status() == 1) {
	       //const Candidate *dauDau = dau->daughter(0);
	       const GenParticle *dauDau = (GenParticle *) dau->daughter(0);
	       if (!stablePhoton1) stablePhoton1 = dauDau;
	       else stablePhoton2 = dauDau;
	     }
	     else if (dau->daughter(0)->numberOfDaughters()==2
		      && fabs(dau->daughter(0)->daughter(0)->pdgId()) == fabs(dau->daughter(0)->daughter(1)->pdgId())) {
	       //const Candidate *dauDauDau = dau->daughter(0)->daughter(0); // take either, choose 1st
	       const GenParticle *dauDauDau = (GenParticle *) dau->daughter(0)->daughter(0); // take either, choose 1st
	       cout << "  *Unstable photon pair produced.*" << endl;
	       if (!stablePhoton1) stablePhoton1 = dauDauDau;
	       else stablePhoton2 = dauDauDau;
	     }
	   }
	   else if (nDauDau == 2 && fabs(dau->daughter(0)->pdgId()) == fabs(dau->daughter(1)->pdgId())) {
	     //const Candidate *dauDau = dau->daughter(0); // take either, choose 1st
	     const GenParticle *dauDau = (GenParticle *) dau->daughter(0); // take either, choose 1st
	     cout << "  *Unstable photon pair produced.*" << endl;
	     if (!stablePhoton1) stablePhoton1 = dauDau;
	     else stablePhoton2 = dauDau;
	   }
	   else {
	     cout << "Didnt find Stable Photon!" << endl;
	   }
	   
	   // print dauDau info
	   for (unsigned j=0; j<nDauDau; ++j) {
	     //const Candidate *dauDau = dau->daughter( j );
	     const GenParticle *dauDau = (GenParticle *) dau->daughter(j);
	     cout << "  --dauID: " << dauDau->pdgId() << endl;
	     cout << "     isLastCopy: " << dauDau->isLastCopy() << endl;
	     cout << "     numberOfDaughters: " << dauDau->numberOfDaughters() << endl;
	     cout << "     dauStatus: " << dauDau->status() << endl;
	     cout << "     dauPt, Eta, Phi: " << dauDau->pt() << ", " << dauDau->eta() << ", " << dauDau->phi() << endl;
	     
	     //print dauDauDau info
	     unsigned nDauDauDau = dauDau->numberOfDaughters();
	     for (unsigned j=0; j<nDauDauDau; ++j) {
	       //const Candidate *dauDauDau = dauDau->daughter( j );
	       const GenParticle *dauDauDau = (GenParticle *) dauDau->daughter(j);
	       cout << "   ---dauID: " << dauDauDau->pdgId() << endl;
	       cout << "       isLastCopy: " << dauDauDau->isLastCopy() << endl;
	       cout << "       numberOfDaughters: " << dauDauDau->numberOfDaughters() << endl;
	       cout << "       dauStatus: " << dauDauDau->status() << endl;
	       cout << "       dauPt, Eta, Phi: " << dauDauDau->pt() << ", " << dauDauDau->eta() << ", " << dauDauDau->phi() << endl;
	     } // end dauDauDau loop
	   } // end dauDau loop
	 
	 } // end 2 photon check
	 else cout << "Didnt find 2 photons from RSG decay!" << endl;
       } // end dau loop
     } // end RSG check
   } // end genParticle loop


   // what if some of the unstable photons arent found? 
   if (!unstablePhoton1) {
     cout << "Couldnt find Unstable Photon 1 !" <<endl;
     fUnstablePhoton1Info.status = -99;
     fTree->Fill();
     return;
   }
   if (!unstablePhoton2) {
     cout << "Couldnt find Unstable Photon 2 !" <<endl;
     fUnstablePhoton2Info.status = -99;
     fTree->Fill();
     return;
   }

   // what if some of the stable photons arent found? 
   if (!stablePhoton1) {
     cout << "Couldnt find Stable Photon 1 !" <<endl;
     fStablePhoton1Info.status = -888.888;
   }
   if (!stablePhoton2) {
     cout << "Couldnt find Stable Photon 2 !" <<endl;
     fStablePhoton2Info.status = -888.888;
   }
  
   // reorder the unstable photons by pt,
   // which are correlted to the stable photons
   if (unstablePhoton2->pt()>unstablePhoton1->pt()) {
     const reco::GenParticle *tempUnstablePhoton = unstablePhoton1;
     unstablePhoton1 = unstablePhoton2;
     unstablePhoton2 = tempUnstablePhoton;
     //if (stablePhoton1 && stablePhoton2) {
     const reco::GenParticle *tempStablePhoton = stablePhoton1;
     stablePhoton1 = stablePhoton2; // could be NULL
     stablePhoton2 = tempStablePhoton; // could be NULL
     // }
   }

   // if stablePhoton pair produces, store pdgId of particle it pair produces to
   if (stablePhoton1 && stablePhoton1->pdgId() != 22) {
     cout << "StablePhoton1 is part of pair production." << endl;
     fStablePhoton1Info.status = -77;
     fStablePhoton1Info.PdgId = fabs(stablePhoton1->pdgId());
     fStablePhoton1Info.pt= -77.77;
     fStablePhoton1Info.eta= -77.77;
     fStablePhoton1Info.phi= -77.77;
   }
   if (stablePhoton2 && stablePhoton2->pdgId() != 22) {
     cout << "StablePhoton2 is part of pair production." << endl;
     fStablePhoton2Info.status = -77;
     fStablePhoton2Info.PdgId = fabs(stablePhoton2->pdgId());
     fStablePhoton2Info.pt= -77.77;
     fStablePhoton2Info.eta= -77.77;
     fStablePhoton2Info.phi= -77.77;
   }
      
   ExoDiPhotons::FillMCTrueObjectInfo(fUnstablePhoton1Info,unstablePhoton1);
   ExoDiPhotons::FillMCTrueObjectInfo(fUnstablePhoton2Info,unstablePhoton2);
   
   if(stablePhoton1 && stablePhoton1->pdgId() == 22) ExoDiPhotons::FillMCTrueObjectInfo(fStablePhoton1Info,stablePhoton1);
   if(stablePhoton2 && stablePhoton2->pdgId() == 22) ExoDiPhotons::FillMCTrueObjectInfo(fStablePhoton2Info,stablePhoton2);

   cout << endl;
   cout << "UnstablePhoton1 status, pdgId, pt: " << fUnstablePhoton1Info.status << ", " << fUnstablePhoton1Info.PdgId << ", " << fUnstablePhoton1Info.pt << endl;
   cout << "UnstablePhoton2 status, pdgId, pt: " << fUnstablePhoton2Info.status << ", " << fUnstablePhoton2Info.PdgId << ", " << fUnstablePhoton2Info.pt << endl;
   cout << endl;
   cout << "StablePhoton1 status, pdgId, pt: " << fStablePhoton1Info.status << ", " << fStablePhoton1Info.PdgId << ", " << fStablePhoton1Info.pt << endl;
   cout << "StablePhoton2 status, pdgId, pt: " << fStablePhoton2Info.status << ", " << fStablePhoton2Info.PdgId << ", " << fStablePhoton2Info.pt << endl;
   cout << endl;

   // ================================================================
   // RECO PHOTON IDENTIFICATION
   // Find best match in deltaR between stable photon and reco photon.
   // Also, check for match with unstable photon? If pair production?
   // ================================================================

   ExoDiPhotons::InitRecoPhotonInfo(fRecoPhoton1Info);
   ExoDiPhotons::InitRecoPhotonInfo(fRecoPhoton2Info);

   // for matching signal to reco photons
   const reco::Photon *matchPhoton1 = NULL;
   const reco::Photon *matchPhoton2 = NULL;
    
   //const reco::Photon *tempMatchPhoton1 = NULL;
   //const reco::Photon *tempMatchPhoton2 = NULL;   

   // use 0.2 as min deltaR cut
   double minDeltaR1 = 0.2;
   double minDeltaR2 = 0.2;

   // index needed to retrieve reco info 
   int matchPhoton1Index = -1;
   int matchPhoton2Index = -1;
   int phoIndex = -1;
   
   for (reco::PhotonCollection::const_iterator recoPhoton = photonColl->begin(); recoPhoton!=photonColl->end(); recoPhoton++) {	 
     phoIndex++;
     
     cout << "Reco Photons: " << endl;
     cout << "pt, eta, phi: " << recoPhoton->pt() << ", " << recoPhoton->eta() << ", " << recoPhoton->phi() << endl;
     
     // we have the stable photons
     // lets match each separately
     if (stablePhoton1 && stablePhoton1->pdgId()==22) {
       // there's a CMS function for deltaPhi in DataFormats/Math
       //double deltaPhi = reco::deltaPhi(stablePhoton1->phi(),recoPhoton->phi());
       //double deltaEta = stablePhoton1->eta()-recoPhoton->eta();
       //double deltaR = TMath::Sqrt(deltaPhi*deltaPhi+deltaEta*deltaEta);
       double deltaR = reco::deltaR(recoPhoton->eta(),recoPhoton->phi(),stablePhoton1->eta(),stablePhoton1->phi());

       cout << "StablePhoton1 deltaR   : " << deltaR << endl;
       cout << "StablePhoton1 minDeltaR: " << minDeltaR1 << endl;

       if (deltaR < minDeltaR1) {
	 // then this is the best match so far
	 minDeltaR1 = deltaR;
	 matchPhoton1 = &(*recoPhoton); //deref the iter to get what it points to
	 matchPhoton1Index = phoIndex;
       }       
     }
     if (stablePhoton2 && stablePhoton2->pdgId()==22) {
       // there's a CMS function for deltaPhi in DataFormats/Math
       //double deltaPhi = reco::deltaPhi(stablePhoton2->phi(),recoPhoton->phi());
       //double deltaEta = stablePhoton2->eta()-recoPhoton->eta();
       //double deltaR = TMath::Sqrt(deltaPhi*deltaPhi+deltaEta*deltaEta);
       double deltaR = reco::deltaR(recoPhoton->eta(),recoPhoton->phi(),stablePhoton2->eta(),stablePhoton2->phi());       

       cout << "StablePhoton2 deltaR   : " << deltaR << endl;
       cout << "StablePhoton2 minDeltaR: " << minDeltaR2 << endl;
		 
       if (deltaR < minDeltaR2) {
	 // then this is the best match so far
	 minDeltaR2 = deltaR;
	 matchPhoton2 = &(*recoPhoton); //deref the iter to get what it points to
	 matchPhoton2Index = phoIndex;
       }       
     }
     
   } //end recoPhoton loop to match to the present signal photon
	       

   cout << endl;

   
   // if no match found, then the match photon pointers are certainly NULL
   if(matchPhoton1) {		 
     //      cout << "Matched to signal photon!" <<endl;
     //      cout << "Reco photon et, eta, phi = " << matchPhoton1->et() <<", "<<matchPhoton1->eta()<< ", "<< matchPhoton1->phi();
     //      cout << "; eMax/e3x3 = " << matchPhoton1->maxEnergyXtal()/matchPhoton1->e3x3();
     //if(iEvent.id().event()==19486)
     //{
     //  cout << "MatchPhoton1: Event number: " << iEvent.id().event();
     // cout << "; et = " << matchPhoton1->et() << endl;
     // cout << "; hadOverEm = " << matchPhoton1->hadronicOverEm() << endl;
     // cout << "; trkIso = " << matchPhoton1->trkSumPtHollowConeDR04() << endl;
     // cout << "; trkIsoCut = " << (2.0+0.001*matchPhoton1->pt()+0.0167*fRho25) << endl;
     // cout << "; ecalIso = " << matchPhoton1->ecalRecHitSumEtConeDR04() << endl;
     // cout << "; ecalIsoCut = " << (4.2 + 0.006*matchPhoton1->pt()+0.183*fRho25) << endl;
     // cout << "; hcalIso = " << matchPhoton1->hcalTowerSumEtConeDR04() << endl;
     // cout << "; hcalIsoCut = " << (2.2 + 0.0025*matchPhoton1->pt()+0.062*fRho25) << endl;
     // cout << "; pixelSeed = " << matchPhoton1->hasPixelSeed() << endl;
     // cout << "; sigmaIetaIeta = " << matchPhoton1->sigmaIetaIeta() << endl;
     // cout << "; isGap = " <<  ExoDiPhotons::isGapPhoton(&(*matchPhoton1)) << endl;
     // cout << "; isSpike = " <<  ExoDiPhotons::isASpike(&(*matchPhoton1)) << endl; 
     // cout << "; passesHadOverEm = " << ExoDiPhotons::passesHadOverEmCut(matchPhoton1) << endl;
     // cout << ": passesTrkIso = " << ExoDiPhotons::passesTrkIsoCut(matchPhoton1,fRho25) << endl;
     // cout << "; passesEcalIsoCut = " << ExoDiPhotons::passesEcalIsoCut(matchPhoton1,fRho25) << endl;
     // cout << "; passesHcalIsoCut = " << ExoDiPhotons::passesHcalIsoCut(matchPhoton1,fRho25) << endl;
     // cout << "; passesSigmaIetaIetaCut = " << ExoDiPhotons::passesSigmaIetaIetaCut(matchPhoton1) << endl;
     // cout << endl << endl;
     //}
     
     
     // fill info into tree
     ExoDiPhotons::FillRecoPhotonInfo(fRecoPhoton1Info,matchPhoton1,lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent,iSetup);
     // SIC add all this to the function above?
     fRecoPhoton1Info.hasMatchedPromptElec = ConversionTools::hasMatchedPromptElectron(matchPhoton1->superCluster(), hElectrons, hConversions, beamSpot.position());


     edm::Ptr<reco::Photon> matchPho1Ptr(photonColl,matchPhoton1Index);

     fRecoPhoton1Info.sigmaIetaIeta = (*full5x5SigmaIEtaIEtaMap)[matchPho1Ptr];
     fRecoPhoton1Info.PFIsoCharged03 = (*phoChargedIsolationMap)[matchPho1Ptr];
     fRecoPhoton1Info.PFIsoNeutral03 = (*phoNeutralHadronIsolationMap)[matchPho1Ptr];
     fRecoPhoton1Info.PFIsoPhoton03 = (*phoPhotonIsolationMap)[matchPho1Ptr];
     fRecoPhoton1Info.PFIsoAll03 = fRecoPhoton1Info.PFIsoCharged03 + fRecoPhoton1Info.PFIsoNeutral03 + fRecoPhoton1Info.PFIsoPhoton03;

     
     float matchPho1Eta = abs(matchPho1Ptr->superCluster()->eta());
     
     //now the corrected PF isolation variables
     fRecoPhoton1Info.rhocorPFIsoCharged03 = std::max((float)0.0,(float)fRecoPhoton1Info.PFIsoCharged03-rho_*effAreaChHadrons_.getEffectiveArea(matchPho1Eta));
     fRecoPhoton1Info.rhocorPFIsoNeutral03 = std::max((float)0.0,(float)fRecoPhoton1Info.PFIsoNeutral03-rho_*effAreaNeuHadrons_.getEffectiveArea(matchPho1Eta));
     fRecoPhoton1Info.rhocorPFIsoPhoton03 = std::max((float)0.0,(float)fRecoPhoton1Info.PFIsoPhoton03-rho_*effAreaPhotons_.getEffectiveArea(matchPho1Eta));
     fRecoPhoton1Info.rhocorPFIsoAll03 = fRecoPhoton1Info.rhocorPFIsoCharged03 + fRecoPhoton1Info.rhocorPFIsoNeutral03 + fRecoPhoton1Info.rhocorPFIsoPhoton03;


     // highPt ID
     std::vector<double> matchPhoton1EA = ExoDiPhotons::EffectiveAreas(matchPhoton1,MethodID,CategoryPFID);
     fRecoPhoton1Info.EAPhotonHighPtID = matchPhoton1EA[2];
     fRecoPhoton1Info.alphaPhotonHighPtID = ExoDiPhotons::alphaPhotonHighPtID(matchPhoton1);
     fRecoPhoton1Info.kappaPhotonHighPtID = ExoDiPhotons::kappaPhotonHighPtID(matchPhoton1);;
     fRecoPhoton1Info.corPhotonIsoHighPtID = ExoDiPhotons::corPhoIsoHighPtID(matchPhoton1,MethodID,CategoryPFID,fRecoPhoton1Info.PFIsoPhoton03,rho_);
     fRecoPhoton1Info.isHighPtPFPhoton = ExoDiPhotons::passHighPtID(matchPhoton1,MethodID,CategoryPFID,fRecoPhoton1Info.PFIsoCharged03,
								      fRecoPhoton1Info.PFIsoPhoton03,fRecoPhoton1Info.sigmaIetaIeta,rho_,!fRecoPhoton1Info.hasMatchedPromptElec);

     /*
     cout << "fRecoPhoton1Info.EAPhotonHighPtID: " << fRecoPhoton1Info.EAPhotonHighPtID << endl;
     cout << "fRecoPhoton1Info.alphaPhotonHighPtID: " << fRecoPhoton1Info.alphaPhotonHighPtID << endl;
     cout << "fRecoPhoton1Info.kappaPhotonHighPtID: " << fRecoPhoton1Info.kappaPhotonHighPtID << endl;
     cout << "fRecoPhoton1Info.corPhotonIsoHighPtID : " << fRecoPhoton1Info.corPhotonIsoHighPtID << endl;
     cout << "fRecoPhoton1Info.isHighPtPFPhoton : " << fRecoPhoton1Info.isHighPtPFPhoton << endl;
     */
     
     //fill 02 and 04 cones?



     //we retrieve the effective areas
     //Remember effareaCH = 1st, effareaNH = 2nd, effareaPH = 3rd
     // std::vector<double> photon1TorFEffAreas = ExoDiPhotons::EffectiveAreas(matchPhoton1);
     // double pfisoall = isolator03.fGetIsolation(matchPhoton1,pfCandidates,firstVtx,vertexColl);
     // double rhocorPFIsoCH = max(isolator03.getIsolationCharged()-fRho25*photon1TorFEffAreas[0],0.);
     //double rhocorPFIsoNH = max(isolator03.getIsolationNeutral()-fRho25*photon1TorFEffAreas[1],0.);
     //double rhocorPFIsoPH = max(isolator03.getIsolationPhoton()-fRho25*photon1TorFEffAreas[2],0.);
     
     //fRecoPhoton1Info.PFIsoAll03 = pfisoall;
     //fRecoPhoton1Info.PFIsoCharged03 = isolator03.getIsolationCharged();
     //fRecoPhoton1Info.PFIsoNeutral03 = isolator03.getIsolationNeutral();
     //fRecoPhoton1Info.PFIsoPhoton03 = isolator03.getIsolationPhoton();      
/*
     fRecoPhoton1Info.PFIsoAll04 = isolator04.fGetIsolation(matchPhoton1,pfCandidates,firstVtx,vertexColl);
     fRecoPhoton1Info.PFIsoCharged04 = isolator04.getIsolationCharged();
     fRecoPhoton1Info.PFIsoNeutral04 = isolator04.getIsolationNeutral();
     fRecoPhoton1Info.PFIsoPhoton04 = isolator04.getIsolationPhoton();      

     fRecoPhoton1Info.PFIsoAll02 = isolator02.fGetIsolation(matchPhoton1,pfCandidates,firstVtx,vertexColl);
     fRecoPhoton1Info.PFIsoCharged02 = isolator02.getIsolationCharged();
     fRecoPhoton1Info.PFIsoNeutral02 = isolator02.getIsolationNeutral();
     fRecoPhoton1Info.PFIsoPhoton02 = isolator02.getIsolationPhoton();      

     //now the corrected PF isolation variables
     fRecoPhoton1Info.rhocorPFIsoCharged04 = max(fRecoPhoton1Info.PFIsoCharged04-fRho25*photon1TorFEffAreas[0],0.);
     fRecoPhoton1Info.rhocorPFIsoNeutral04 = max(fRecoPhoton1Info.PFIsoNeutral04-fRho25*photon1TorFEffAreas[1],0.);
     fRecoPhoton1Info.rhocorPFIsoPhoton04 = max(fRecoPhoton1Info.PFIsoPhoton04-fRho25*photon1TorFEffAreas[2],0.);
     fRecoPhoton1Info.rhocorPFIsoAll04 = fRecoPhoton1Info.rhocorPFIsoCharged04 + fRecoPhoton1Info.rhocorPFIsoNeutral04 + fRecoPhoton1Info.rhocorPFIsoPhoton04;
 */
     //fRecoPhoton1Info.rhocorPFIsoCharged03 = rhocorPFIsoCH;
     //fRecoPhoton1Info.rhocorPFIsoNeutral03 = rhocorPFIsoNH;
     //fRecoPhoton1Info.rhocorPFIsoPhoton03 = rhocorPFIsoPH;
     //fRecoPhoton1Info.rhocorPFIsoAll03 = fRecoPhoton1Info.rhocorPFIsoCharged03 + fRecoPhoton1Info.rhocorPFIsoNeutral03 + fRecoPhoton1Info.rhocorPFIsoPhoton03;
     /*
     fRecoPhoton1Info.rhocorPFIsoCharged02 = max(fRecoPhoton1Info.PFIsoCharged02-fRho25*photon1TorFEffAreas[0],0.);
     fRecoPhoton1Info.rhocorPFIsoNeutral02 = max(fRecoPhoton1Info.PFIsoNeutral02-fRho25*photon1TorFEffAreas[1],0.);
     fRecoPhoton1Info.rhocorPFIsoPhoton02 = max(fRecoPhoton1Info.PFIsoPhoton02-fRho25*photon1TorFEffAreas[2],0.);
     fRecoPhoton1Info.rhocorPFIsoAll02 = fRecoPhoton1Info.rhocorPFIsoCharged02 + fRecoPhoton1Info.rhocorPFIsoNeutral02 + fRecoPhoton1Info.rhocorPFIsoPhoton02;
     */

     
     fRecoPhoton1Info.isTightPFPhoton = (*tight_id_decisions)[matchPho1Ptr];
     fRecoPhoton1Info.isMediumPFPhoton = (*medium_id_decisions)[matchPho1Ptr];
     fRecoPhoton1Info.isLoosePFPhoton = (*loose_id_decisions)[matchPho1Ptr]; 


     // fill DET info?

     /*
     // now we fill the ID variables for det and PFID
     if(ExoDiPhotons::isTightPhoton(&(*matchPhoton1),fRho25) &&
        !ExoDiPhotons::isGapPhoton(&(*matchPhoton1)) &&
        !ExoDiPhotons::isASpike(&(*matchPhoton1))  )
       fRecoPhoton1Info.isTightDetPhoton = true;
     else
       fRecoPhoton1Info.isTightDetPhoton = false;
     */

     /*    
     // test the conversion safe electron veto
     bool passElecVeto = !ConversionTools::hasMatchedPromptElectron(matchPhoton1->superCluster(), hElectrons, hConversions, beamSpot.position());
     if(ExoDiPhotons::isPFTightPhoton(&(*matchPhoton1),rhocorPFIsoCH,rhocorPFIsoNH,rhocorPFIsoPH,passElecVeto,TString("Tight")) && 
	!ExoDiPhotons::isGapPhoton(&(*matchPhoton1)) && 
	!ExoDiPhotons::isASpike(&(*matchPhoton1))  )
       fRecoPhoton1Info.isTightPFPhoton = true;
     else
       fRecoPhoton1Info.isTightPFPhoton = false;
     if(ExoDiPhotons::isPFTightPhoton(&(*matchPhoton1),rhocorPFIsoCH,rhocorPFIsoNH,rhocorPFIsoPH,passElecVeto,TString("Medium")) && 
	!ExoDiPhotons::isGapPhoton(&(*matchPhoton1)) && 
	!ExoDiPhotons::isASpike(&(*matchPhoton1))  )
       fRecoPhoton1Info.isMediumPFPhoton = true;
     else
       fRecoPhoton1Info.isMediumPFPhoton = false;
     if(ExoDiPhotons::isPFTightPhoton(&(*matchPhoton1),rhocorPFIsoCH,rhocorPFIsoNH,rhocorPFIsoPH,passElecVeto,TString("Loose")) && 
	!ExoDiPhotons::isGapPhoton(&(*matchPhoton1)) && 
	!ExoDiPhotons::isASpike(&(*matchPhoton1))  )
       fRecoPhoton1Info.isLoosePFPhoton = true;
     else
       fRecoPhoton1Info.isLoosePFPhoton = false;
     */  
   }
   else {
     cout << "No match to stablePhoton1!" << endl;
     fRecoPhoton1Info.pt = -777.777;
   }
   
   if(matchPhoton2) {		 
     //      cout << "Matched to signal photon!" <<endl;
     //      cout << "Reco photon et, eta, phi = " << matchPhoton2->et() <<", "<<matchPhoton2->eta()<< ", "<< matchPhoton2->phi();

     // fill info into tree
     ExoDiPhotons::FillRecoPhotonInfo(fRecoPhoton2Info,matchPhoton2,lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);

     fRecoPhoton2Info.hasMatchedPromptElec = ConversionTools::hasMatchedPromptElectron(matchPhoton2->superCluster(), hElectrons, hConversions, beamSpot.position());

     edm::Ptr<reco::Photon> matchPho2Ptr(photonColl,matchPhoton2Index);
     
     fRecoPhoton2Info.sigmaIetaIeta = (*full5x5SigmaIEtaIEtaMap)[matchPho2Ptr];
     fRecoPhoton2Info.PFIsoCharged03 = (*phoChargedIsolationMap)[matchPho2Ptr];
     fRecoPhoton2Info.PFIsoNeutral03 = (*phoNeutralHadronIsolationMap)[matchPho2Ptr];
     fRecoPhoton2Info.PFIsoPhoton03 = (*phoPhotonIsolationMap)[matchPho2Ptr];
     fRecoPhoton2Info.PFIsoAll03 = fRecoPhoton2Info.PFIsoCharged03 + fRecoPhoton2Info.PFIsoNeutral03 + fRecoPhoton2Info.PFIsoPhoton03;
     
     float matchPho2Eta = abs(matchPho2Ptr->superCluster()->eta());
     
     //now the corrected PF isolation variables
     fRecoPhoton2Info.rhocorPFIsoCharged03 = std::max((float)0.0,(float)fRecoPhoton2Info.PFIsoCharged03-rho_*effAreaChHadrons_.getEffectiveArea(matchPho2Eta));
     fRecoPhoton2Info.rhocorPFIsoNeutral03 = std::max((float)0.0,(float)fRecoPhoton2Info.PFIsoNeutral03-rho_*effAreaNeuHadrons_.getEffectiveArea(matchPho2Eta));
     fRecoPhoton2Info.rhocorPFIsoPhoton03 = std::max((float)0.0,(float)fRecoPhoton2Info.PFIsoPhoton03-rho_*effAreaPhotons_.getEffectiveArea(matchPho2Eta));
     fRecoPhoton2Info.rhocorPFIsoAll03 = fRecoPhoton2Info.rhocorPFIsoCharged03 + fRecoPhoton2Info.rhocorPFIsoNeutral03 + fRecoPhoton2Info.rhocorPFIsoPhoton03;

     // highPt ID
     std::vector<double> matchPhoton2EA = ExoDiPhotons::EffectiveAreas(matchPhoton2,MethodID,CategoryPFID);
     fRecoPhoton2Info.EAPhotonHighPtID = matchPhoton2EA[2];
     fRecoPhoton2Info.alphaPhotonHighPtID = ExoDiPhotons::alphaPhotonHighPtID(matchPhoton2);
     fRecoPhoton2Info.kappaPhotonHighPtID = ExoDiPhotons::kappaPhotonHighPtID(matchPhoton2);;
     fRecoPhoton2Info.corPhotonIsoHighPtID = ExoDiPhotons::corPhoIsoHighPtID(matchPhoton2,MethodID,CategoryPFID,fRecoPhoton2Info.PFIsoPhoton03,rho_);
     fRecoPhoton2Info.isHighPtPFPhoton = ExoDiPhotons::passHighPtID(matchPhoton2,MethodID,CategoryPFID,fRecoPhoton2Info.PFIsoCharged03,
								      fRecoPhoton2Info.PFIsoPhoton03,fRecoPhoton2Info.sigmaIetaIeta,rho_,!fRecoPhoton2Info.hasMatchedPromptElec);
     
     fRecoPhoton2Info.isTightPFPhoton = (*tight_id_decisions)[matchPho2Ptr];
     fRecoPhoton2Info.isMediumPFPhoton = (*medium_id_decisions)[matchPho2Ptr];
     fRecoPhoton2Info.isLoosePFPhoton = (*loose_id_decisions)[matchPho2Ptr]; 
   }
   else {
     cout << "No match to stablePhoton2!" << endl;
     fRecoPhoton2Info.pt = -777.777;
   }

   cout << "RecoPhoton1 pt, eta, phi: " << fRecoPhoton1Info.pt << ", " << fRecoPhoton1Info.eta << ", " << fRecoPhoton1Info.phi << endl;
   cout << "RecoPhoton2 pt, eta, phi: " << fRecoPhoton2Info.pt << ", " << fRecoPhoton2Info.eta << ", " << fRecoPhoton2Info.phi << endl;
   cout << endl;
   cout << endl;

   // =================================
   // DIPHOTON IDENTIFICATION
   // If 2 photons, fill diphoton info.
   //  ================================

   ExoDiPhotons::InitDiphotonInfo(fDiphotonGenUnstableInfo);
   ExoDiPhotons::InitDiphotonInfo(fDiphotonGenStableInfo);
   ExoDiPhotons::InitDiphotonInfo(fDiphotonRecoInfo);
   
   if (unstablePhoton1 && unstablePhoton2) {
     ExoDiPhotons::FillDiphotonInfo(fDiphotonGenUnstableInfo,unstablePhoton1,unstablePhoton2);
   }

   if (stablePhoton1 && stablePhoton2 && stablePhoton1->pdgId()==22 && stablePhoton2->pdgId()==22) {
     ExoDiPhotons::FillDiphotonInfo(fDiphotonGenStableInfo,stablePhoton1,stablePhoton2);
   }

   if (matchPhoton1 && matchPhoton2) {
     ExoDiPhotons::FillDiphotonInfo(fDiphotonRecoInfo,matchPhoton1,matchPhoton2);
   }
   
   // Fill the tree every event.
   fTree->Fill();
   

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
ExoDiPhotonSignalMCAnalyzer::beginJob()
{
  LumiWeights = edm::LumiReWeighting(PUMCFileName,PUDataFileName,PUMCHistName,PUDataHistName);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ExoDiPhotonSignalMCAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ExoDiPhotonSignalMCAnalyzer);
