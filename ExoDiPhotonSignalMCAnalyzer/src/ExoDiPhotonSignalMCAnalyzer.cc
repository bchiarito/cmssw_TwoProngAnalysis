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
  edm::InputTag      fPhotonTag;       //select photon collection 
  double             fMin_pt;          // min pt cut (photons)
  edm::InputTag      fHltInputTag;     // hltResults
  edm::InputTag      fRho25Tag;
  edm::InputTag      pileupCollectionTag;
  edm::LumiReWeighting    LumiWeights;
  
  bool               fkRemoveSpikes;   // option to remove spikes before filling tree
  bool               fkRequireGenEventInfo;  // generated information for RS graviton files
  string             PUMCFileName;
  string             PUDataFileName;
  string             PUDataHistName;
  string             PUMCHistName;
  string             fPFIDCategory;
  string             fIDMethod;
  
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
  
  ExoDiPhotons::mcTrueObjectInfo_t fSignalPhoton1Info; // leading signal photon
  ExoDiPhotons::mcTrueObjectInfo_t fSignalPhoton2Info;
  ExoDiPhotons::recoPhotonInfo_t fRecoPhotonInfo1; // leading matched reco photon 
  ExoDiPhotons::recoPhotonInfo_t fRecoPhotonInfo2; // second photon
  
  ExoDiPhotons::diphotonInfo_t fDiphotonSignalInfo;
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

  std::cout << "ExoDiPhotonAnalyzer: ID Method used " << fIDMethod.c_str()
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
 
  // now with CommonClasses, use the string defined in the header

  fTree->Branch("Event",&fEventInfo,ExoDiPhotons::eventInfoBranchDefString.c_str());
  fTree->Branch("Vtx",&fVtxInfo,ExoDiPhotons::vtxInfoBranchDefString.c_str());
  fTree->Branch("BeamSpot",&fBeamSpotInfo,ExoDiPhotons::beamSpotInfoBranchDefString.c_str());
  fTree->Branch("TrigHLT",&fHLTInfo,ExoDiPhotons::hltTrigBranchDefString.c_str());

  fTree->Branch("GenPhoton1",&fSignalPhoton1Info,ExoDiPhotons::mcTrueObjectInfoBranchDefString.c_str());
  fTree->Branch("GenPhoton2",&fSignalPhoton2Info,ExoDiPhotons::mcTrueObjectInfoBranchDefString.c_str());

  fTree->Branch("Photon1",&fRecoPhotonInfo1,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  fTree->Branch("Photon2",&fRecoPhotonInfo2,ExoDiPhotons::recoPhotonBranchDefString.c_str());
  
  // signal diphoton info? eg to probe true MC width?
  // reco diphoton info?
  
  fTree->Branch("DiphotonGen",&fDiphotonSignalInfo,ExoDiPhotons::diphotonInfoBranchDefString.c_str());
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

   const reco::PhotonCollection *myPhotonColl = photonColl.product();
   cout<<"photoncoll size "<<myPhotonColl->size()<<endl;
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
   
   // new approach - 
   // make vector of all selected Photons (tight ID, not spike, min pt, etc)
   // then sort at end by pt
   // also allows to count how often a third photon could be considered a candidate

   std::vector<reco::Photon> selectedPhotons; 

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
    
     //cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;
     
     float abseta = fabs( recoPhoton->superCluster()->eta());
     float isoChargedHadronsWithEA = std::max( (float)0.0, chIso - rho_*effAreaChHadrons_.getEffectiveArea(abseta));
     float isoNeutralHadronsWithEA = std::max( (float)0.0, nhIso - rho_*effAreaNeuHadrons_.getEffectiveArea(abseta));
     float isoPhotonsWithEA = std::max( (float)0.0, phIso - rho_*effAreaPhotons_.getEffectiveArea(abseta));
    
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

     //CAREFUL UNCOMMENT THAT WHEN DONE
     if(recoPhoton->pt() < fMin_pt) continue;

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

     /*      
     // also check for fakeable objects
     //if(ExoDiPhotons::isFakeableObject(&(*recoPhoton)) ) {
     if(MethodID.Contains("Detector")){
       if(ExoDiPhotons::isFakeableObject(&(*recoPhoton),fRho25) ) {
	 
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
     else if(MethodID.Contains("ParticleFlow")){
       if(ExoDiPhotons::isPFFakeableObject(&(*recoPhoton),rhocorPFIsoCH,rhocorPFIsoNH,rhocorPFIsoPH,full5x5sigmaIetaIeta,passelecveto,CategoryPFID) ) {
	 
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
     */
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
     int myIndex = -1;

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
     //need??
     //allTightOrFakeableObjects[0].second = false;//used in sorting, now it's faked if this 2 tight exception comes up

     //Now we store all PF isolation variables for the 1st photon (tight exception)

     fRecoPhotonInfo1.sigmaIetaIeta = (*full5x5SigmaIEtaIEtaMap)[TObject1Ptr];
     fRecoPhotonInfo1.PFIsoCharged03 = (*phoChargedIsolationMap)[TObject1Ptr];
     fRecoPhotonInfo1.PFIsoNeutral03 = (*phoNeutralHadronIsolationMap)[TObject1Ptr];
     fRecoPhotonInfo1.PFIsoPhoton03 = (*phoPhotonIsolationMap)[TObject1Ptr];
     fRecoPhotonInfo1.PFIsoAll03 = fRecoPhotonInfo1.PFIsoCharged03 + fRecoPhotonInfo1.PFIsoNeutral03 + fRecoPhotonInfo1.PFIsoPhoton03;

     float tObject1Eta = abs(TObject1Ptr->superCluster()->eta());

     //now the corrected PF isolation variables
     fRecoPhotonInfo1.rhocorPFIsoCharged03 = std::max((float)0.0,(float)fRecoPhotonInfo1.PFIsoCharged03-rho_*effAreaChHadrons_.getEffectiveArea(tObject1Eta));
     fRecoPhotonInfo1.rhocorPFIsoNeutral03 = std::max((float)0.0,(float)fRecoPhotonInfo1.PFIsoNeutral03-rho_*effAreaChHadrons_.getEffectiveArea(tObject1Eta));
     fRecoPhotonInfo1.rhocorPFIsoPhoton03 = std::max((float)0.0,(float)fRecoPhotonInfo1.PFIsoPhoton03-rho_*effAreaChHadrons_.getEffectiveArea(tObject1Eta));
     fRecoPhotonInfo1.rhocorPFIsoAll03 = fRecoPhotonInfo1.rhocorPFIsoCharged03 + fRecoPhotonInfo1.rhocorPFIsoNeutral03 + fRecoPhotonInfo1.rhocorPFIsoPhoton03;

     ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo2,&selectedPhotons[1],lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);
     fRecoPhotonInfo2.hasMatchedPromptElec = ConversionTools::hasMatchedPromptElectron((&selectedPhotons[1])->superCluster(), hElectrons, hConversions, beamSpot.position());

     fRecoPhotonInfo2.isTightPFPhoton = (*tight_id_decisions)[TObject2Ptr];
     fRecoPhotonInfo2.isMediumPFPhoton = (*medium_id_decisions)[TObject2Ptr];
     fRecoPhotonInfo2.isLoosePFPhoton = (*loose_id_decisions)[TObject2Ptr]; 

     //TO DISENTANGLE BETWEEN MINIAOD AND AOD
     fRecoPhotonInfo2.isFakeable = false;
     //need??
     //allTightOrFakeableObjects[1].second = false;

     //Now we store all PF isolation variables for the 2st photon (tight exception)

     fRecoPhotonInfo2.sigmaIetaIeta = (*full5x5SigmaIEtaIEtaMap)[TObject2Ptr];
     fRecoPhotonInfo2.PFIsoCharged03 = (*phoChargedIsolationMap)[TObject2Ptr];
     fRecoPhotonInfo2.PFIsoNeutral03 = (*phoNeutralHadronIsolationMap)[TObject2Ptr];
     fRecoPhotonInfo2.PFIsoPhoton03 = (*phoPhotonIsolationMap)[TObject2Ptr];
     fRecoPhotonInfo2.PFIsoAll03 = fRecoPhotonInfo2.PFIsoCharged03 + fRecoPhotonInfo2.PFIsoNeutral03 + fRecoPhotonInfo2.PFIsoPhoton03;

     float tObject2Eta = abs(TObject2Ptr->superCluster()->eta());

     //now the corrected PF isolation variables
     fRecoPhotonInfo2.rhocorPFIsoCharged03 = std::max((float)0.0,(float)fRecoPhotonInfo2.PFIsoCharged03-rho_*effAreaChHadrons_.getEffectiveArea(tObject2Eta));
     fRecoPhotonInfo2.rhocorPFIsoNeutral03 = std::max((float)0.0,(float)fRecoPhotonInfo2.PFIsoNeutral03-rho_*effAreaChHadrons_.getEffectiveArea(tObject2Eta));
     fRecoPhotonInfo2.rhocorPFIsoPhoton03 = std::max((float)0.0,(float)fRecoPhotonInfo2.PFIsoPhoton03-rho_*effAreaChHadrons_.getEffectiveArea(tObject2Eta));
     fRecoPhotonInfo2.rhocorPFIsoAll03 = fRecoPhotonInfo2.rhocorPFIsoCharged03 + fRecoPhotonInfo2.rhocorPFIsoNeutral03 + fRecoPhotonInfo2.rhocorPFIsoPhoton03;

     // fill diphoton info                                                                                                   
     // add fDiphotonRecoInfo here -- change, eventually
     ExoDiPhotons::FillDiphotonInfo(fDiphotonRecoInfo,&selectedPhotons[0],&selectedPhotons[1]);

   } //end of 2 TT exception


   //   Handle<HepMCProduct> mcProduct;
   //   iEvent.getByLabel("generator",mcProduct);
   //   const HepMC::GenEvent *mcEvent = mcProduct->GetEvent();
   
   //     mcEvent->print( std::cout );

   // I used to use HepMCProduct, and access particles as HepMC::GenParticles
   // I now use the reco::GenParticles collection
   Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByLabel("genParticles",genParticles);

   if(!genParticles.isValid()) {
     cout << "No Gen Particles collection!" << endl;
     return;
   }

   fSignalPhoton1Info.status = -999999;
   fSignalPhoton1Info.PdgId = -999999;
   fSignalPhoton1Info.MotherPdgId = -999999;
   fSignalPhoton1Info.GrandmotherPdgId = -999999;
   fSignalPhoton1Info.pt = -999999.99;
   fSignalPhoton1Info.eta = -999999.99;
   fSignalPhoton1Info.phi = -999999.99;

   fSignalPhoton2Info.status = -999999;
   fSignalPhoton2Info.PdgId = -999999;
   fSignalPhoton2Info.MotherPdgId = -999999;
   fSignalPhoton2Info.GrandmotherPdgId = -999999;
   fSignalPhoton2Info.pt = -999999.99;
   fSignalPhoton2Info.eta = -999999.99;
   fSignalPhoton2Info.phi = -999999.99;

   const reco::GenParticle *signalPhoton1 = NULL;
   const reco::GenParticle *signalPhoton2 = NULL;

   // for matching signal to reco photons
   const reco::Photon *matchPhoton1 = NULL;
   const reco::Photon *matchPhoton2 = NULL;

   //cout<<"HereGenLoop"<<endl;
   for(reco::GenParticleCollection::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); ++genParticle) {
   //for(reco::GenParticle::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); ++genParticle) {
     
     // identify the status 1 (ie no further decay) photons
     // which came from the hard-scatt photons (status 3)
     // because I know they're there     
     //cout<<genParticle->status()<<endl;
     //cout<<genParticle->pdgId()<<endl;
     //cout<<"ROLLTIDE"<<endl;
   
     if(genParticle->status()==1&&genParticle->pdgId()==22) {
       //cout<<"1"<<endl;
       if(genParticle->numberOfMothers()>0) {
	 // cout<<"2"<<endl;
	 if(genParticle->mother()->status()==3&&genParticle->mother()->pdgId()==22) {
	   //   cout<<"3"<<endl;
           // further require that this status 3 photon itself came from the RS graviton 
	   if(genParticle->mother()->numberOfMothers()>0) {
	     //     cout<<"4"<<endl;
	     if(genParticle->mother()->mother()->pdgId()==5000039) {
	       //cout<<genParticle->mother()->mother()->pdgId()<<endl;
	       //cout << "MC particle: Status = "<< genParticle->status() << "; pdg id = "<< genParticle->pdgId() << "; pt, eta, phi = "
	       //     << genParticle->pt() << ", "<< genParticle->eta() << ", " << genParticle->phi() << endl;	   	       
	       
	       // now match this signal photon to best recoPhoton
	       const reco::Photon *tempMatchPhoton = NULL;
	       double minDeltaR = 1.0;
	       
	       for(reco::PhotonCollection::const_iterator recoPhoton = photonColl->begin(); recoPhoton!=photonColl->end(); recoPhoton++) {	 
		 // there's a CMS function for deltaPhi in DataFormats/Math
		 double deltaPhi = reco::deltaPhi(genParticle->phi(),recoPhoton->phi());
		 double deltaEta = genParticle->eta()-recoPhoton->eta();
		 double deltaR = TMath::Sqrt(deltaPhi*deltaPhi+deltaEta*deltaEta);
		 
		 if(deltaR<minDeltaR) {
		   // then this is the best match so far
		   minDeltaR = deltaR;
		   tempMatchPhoton = &(*recoPhoton); //deref the iter to get what it points to
		 }
	       } //end recoPhoton loop to match to the present signal photon
	       
	       // now assign our signal and matched photons to 1 or 2
	       if(!signalPhoton1) {
		 // then we havent found the first one yet, so this is it
		 signalPhoton1 = &(*genParticle);
                 //cout<<signalPhoton1<<endl;
		 matchPhoton1 = tempMatchPhoton;
	       }
	       else {
		 // we have already found one, so this is the second
		 signalPhoton2 = &(*genParticle);
		 matchPhoton2 = tempMatchPhoton;
	       }
	     }
	    
	   }
	 }
       }
       
     } //end status 1 req for  photons from RS graviton
     
     // identify other real photons in event
     // what about ISR photons? 
     
     //cout << "MC particle: Status = "<< genParticle->status() << "; pdg id = "<< genParticle->pdgId() << "; pt, eta, phi = " << genParticle->pt() << ", "<< genParticle->eta() << ", " << genParticle->phi() << endl;	   
     
     // what about a photon which converts late in detector?
     // (or, similarly, an electron which brems a photon)
     // how does this, which happens after pythia, get saved in event record?
     
     // or what if it is not even a photon at all, but say a jet?
     // try printing all final state particles above pt cut
     // but remember that jets can have lots of particles!
     // so maybe best not to cut on pt of each particle?
     //      if(genParticle->status()==1) {
     
     //        cout << "MC particle: Status = "<< genParticle->status() << "; pdg id = "<< genParticle->pdgId() << "; pt, eta, phi = " << genParticle->pt() << ", "<< genParticle->eta() << ", " << genParticle->phi();	   
     //        if(genParticle->numberOfMothers()>0) {
     // 	 cout << "; Mother pdg ID = " << genParticle->mother()->pdgId();
     //        }
     //        cout << endl;
     //      }

     
   } //end loop over gen particles


   // what if some of the signal Photons arent found? 
   if(!signalPhoton1) {
     cout << "Couldnt find signal Photon1 !" <<endl;
     fSignalPhoton1Info.status = -99;
     fTree->Fill();
     return;
   }
   if(!signalPhoton2) {
     cout << "Couldnt find signal Photon2 !" <<endl;
     fSignalPhoton2Info.status = -99;
     fTree->Fill();
     return;
   }
   
   // reorder the signalPhotons by pt
   if(signalPhoton2->pt()>signalPhoton1->pt()) {
     const reco::GenParticle *tempSignalPhoton = signalPhoton1;
     const reco::Photon *tempMatchPhoton = matchPhoton1;
     signalPhoton1 = signalPhoton2;
     signalPhoton2 = tempSignalPhoton;
     matchPhoton1 = matchPhoton2;
     matchPhoton2 = tempMatchPhoton;
   }
   
   if(signalPhoton1) {
     ExoDiPhotons::FillMCTrueObjectInfo(fSignalPhoton1Info,signalPhoton1);
     //      fSignalPhoton1Info.status = signalPhoton1->status();
     //      fSignalPhoton1Info.PdgId = signalPhoton1->pdgId();
     //      fSignalPhoton1Info.MotherPdgId = signalPhoton1->mother()->pdgId();
     //      fSignalPhoton1Info.GrandmotherPdgId = signalPhoton1->mother()->mother()->pdgId();
     //      fSignalPhoton1Info.pt = signalPhoton1->pt();
     //      fSignalPhoton1Info.eta = signalPhoton1->eta();
     //      fSignalPhoton1Info.phi = signalPhoton1->phi();
   }

   if(signalPhoton2) {
     ExoDiPhotons::FillMCTrueObjectInfo(fSignalPhoton2Info,signalPhoton2);
     //      fSignalPhoton2Info.status = signalPhoton2->status();
     //      fSignalPhoton2Info.PdgId = signalPhoton2->pdgId();
     //      fSignalPhoton2Info.MotherPdgId = signalPhoton2->mother()->pdgId();
     //      fSignalPhoton2Info.GrandmotherPdgId = signalPhoton2->mother()->mother()->pdgId();
     //      fSignalPhoton2Info.pt = signalPhoton2->pt();
     //      fSignalPhoton2Info.eta = signalPhoton2->eta();
     //      fSignalPhoton2Info.phi = signalPhoton2->phi();		 
   }
   
   /*
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
     ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo1,matchPhoton1,lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent,iSetup);
     // SIC add all this to the function above?
     fRecoPhotonInfo1.hasMatchedPromptElec = ConversionTools::hasMatchedPromptElectron(matchPhoton1->superCluster(), hElectrons, hConversions, beamSpot.position());
     //we retrieve the effective areas
     //Remember effareaCH = 1st, effareaNH = 2nd, effareaPH = 3rd
     std::vector<double> photon1TorFEffAreas = ExoDiPhotons::EffectiveAreas(matchPhoton1);
     double pfisoall = isolator03.fGetIsolation(matchPhoton1,pfCandidates,firstVtx,vertexColl);
     double rhocorPFIsoCH = max(isolator03.getIsolationCharged()-fRho25*photon1TorFEffAreas[0],0.);
     double rhocorPFIsoNH = max(isolator03.getIsolationNeutral()-fRho25*photon1TorFEffAreas[1],0.);
     double rhocorPFIsoPH = max(isolator03.getIsolationPhoton()-fRho25*photon1TorFEffAreas[2],0.);
     
     fRecoPhotonInfo1.PFIsoAll03 = pfisoall;
     fRecoPhotonInfo1.PFIsoCharged03 = isolator03.getIsolationCharged();
     fRecoPhotonInfo1.PFIsoNeutral03 = isolator03.getIsolationNeutral();
     fRecoPhotonInfo1.PFIsoPhoton03 = isolator03.getIsolationPhoton();      

     fRecoPhotonInfo1.PFIsoAll04 = isolator04.fGetIsolation(matchPhoton1,pfCandidates,firstVtx,vertexColl);
     fRecoPhotonInfo1.PFIsoCharged04 = isolator04.getIsolationCharged();
     fRecoPhotonInfo1.PFIsoNeutral04 = isolator04.getIsolationNeutral();
     fRecoPhotonInfo1.PFIsoPhoton04 = isolator04.getIsolationPhoton();      

     fRecoPhotonInfo1.PFIsoAll02 = isolator02.fGetIsolation(matchPhoton1,pfCandidates,firstVtx,vertexColl);
     fRecoPhotonInfo1.PFIsoCharged02 = isolator02.getIsolationCharged();
     fRecoPhotonInfo1.PFIsoNeutral02 = isolator02.getIsolationNeutral();
     fRecoPhotonInfo1.PFIsoPhoton02 = isolator02.getIsolationPhoton();      

     //now the corrected PF isolation variables
     fRecoPhotonInfo1.rhocorPFIsoCharged04 = max(fRecoPhotonInfo1.PFIsoCharged04-fRho25*photon1TorFEffAreas[0],0.);
     fRecoPhotonInfo1.rhocorPFIsoNeutral04 = max(fRecoPhotonInfo1.PFIsoNeutral04-fRho25*photon1TorFEffAreas[1],0.);
     fRecoPhotonInfo1.rhocorPFIsoPhoton04 = max(fRecoPhotonInfo1.PFIsoPhoton04-fRho25*photon1TorFEffAreas[2],0.);
     fRecoPhotonInfo1.rhocorPFIsoAll04 = fRecoPhotonInfo1.rhocorPFIsoCharged04 + fRecoPhotonInfo1.rhocorPFIsoNeutral04 + fRecoPhotonInfo1.rhocorPFIsoPhoton04;

     fRecoPhotonInfo1.rhocorPFIsoCharged03 = rhocorPFIsoCH;
     fRecoPhotonInfo1.rhocorPFIsoNeutral03 = rhocorPFIsoNH;
     fRecoPhotonInfo1.rhocorPFIsoPhoton03 = rhocorPFIsoPH;
     fRecoPhotonInfo1.rhocorPFIsoAll03 = fRecoPhotonInfo1.rhocorPFIsoCharged03 + fRecoPhotonInfo1.rhocorPFIsoNeutral03 + fRecoPhotonInfo1.rhocorPFIsoPhoton03;

     fRecoPhotonInfo1.rhocorPFIsoCharged02 = max(fRecoPhotonInfo1.PFIsoCharged02-fRho25*photon1TorFEffAreas[0],0.);
     fRecoPhotonInfo1.rhocorPFIsoNeutral02 = max(fRecoPhotonInfo1.PFIsoNeutral02-fRho25*photon1TorFEffAreas[1],0.);
     fRecoPhotonInfo1.rhocorPFIsoPhoton02 = max(fRecoPhotonInfo1.PFIsoPhoton02-fRho25*photon1TorFEffAreas[2],0.);
     fRecoPhotonInfo1.rhocorPFIsoAll02 = fRecoPhotonInfo1.rhocorPFIsoCharged02 + fRecoPhotonInfo1.rhocorPFIsoNeutral02 + fRecoPhotonInfo1.rhocorPFIsoPhoton02;

     // now we fill the ID variables for det and PFID
     if(ExoDiPhotons::isTightPhoton(&(*matchPhoton1),fRho25) &&
        !ExoDiPhotons::isGapPhoton(&(*matchPhoton1)) &&
        !ExoDiPhotons::isASpike(&(*matchPhoton1))  )
       fRecoPhotonInfo1.isTightDetPhoton = true;
     else
       fRecoPhotonInfo1.isTightDetPhoton = false;
    
     // test the conversion safe electron veto
     bool passElecVeto = !ConversionTools::hasMatchedPromptElectron(matchPhoton1->superCluster(), hElectrons, hConversions, beamSpot.position());
     if(ExoDiPhotons::isPFTightPhoton(&(*matchPhoton1),rhocorPFIsoCH,rhocorPFIsoNH,rhocorPFIsoPH,passElecVeto,TString("Tight")) && 
	!ExoDiPhotons::isGapPhoton(&(*matchPhoton1)) && 
	!ExoDiPhotons::isASpike(&(*matchPhoton1))  )
       fRecoPhotonInfo1.isTightPFPhoton = true;
     else
       fRecoPhotonInfo1.isTightPFPhoton = false;
     if(ExoDiPhotons::isPFTightPhoton(&(*matchPhoton1),rhocorPFIsoCH,rhocorPFIsoNH,rhocorPFIsoPH,passElecVeto,TString("Medium")) && 
	!ExoDiPhotons::isGapPhoton(&(*matchPhoton1)) && 
	!ExoDiPhotons::isASpike(&(*matchPhoton1))  )
       fRecoPhotonInfo1.isMediumPFPhoton = true;
     else
       fRecoPhotonInfo1.isMediumPFPhoton = false;
     if(ExoDiPhotons::isPFTightPhoton(&(*matchPhoton1),rhocorPFIsoCH,rhocorPFIsoNH,rhocorPFIsoPH,passElecVeto,TString("Loose")) && 
	!ExoDiPhotons::isGapPhoton(&(*matchPhoton1)) && 
	!ExoDiPhotons::isASpike(&(*matchPhoton1))  )
       fRecoPhotonInfo1.isLoosePFPhoton = true;
     else
       fRecoPhotonInfo1.isLoosePFPhoton = false;
   }
   else {
     //     cout << "No match to signal photon1!" <<endl;
     // as short cut for indicating this in tree
     // make sure the pt value is crazy
     fRecoPhotonInfo1.pt = -9999.99;
   }
   
   if(matchPhoton2) {		 
     //      cout << "Matched to signal photon!" <<endl;
     //      cout << "Reco photon et, eta, phi = " << matchPhoton2->et() <<", "<<matchPhoton2->eta()<< ", "<< matchPhoton2->phi();
     //      cout << "; eMax/e3x3 = " << matchPhoton2->maxEnergyXtal()/matchPhoton2->e3x3();
     //if(iEvent.id().event()==19486)
     //{
     //  cout << "MatchPhoton2: Event number: " << iEvent.id().event();
     // cout << "; et = " << matchPhoton2->et() << endl;
     // cout << "; hadOverEm = " << matchPhoton2->hadronicOverEm() << endl;
     // cout << "; trkIso = " << matchPhoton2->trkSumPtHollowConeDR04() << endl;
     // cout << "; trkIsoCut = " << (2.0+0.001*matchPhoton2->pt()+0.0167*fRho25) << endl;
     // cout << "; ecalIso = " << matchPhoton2->ecalRecHitSumEtConeDR04() << endl;
     // cout << "; ecalIsoCut = " << (4.2 + 0.006*matchPhoton2->pt()+0.183*fRho25) << endl;
     // cout << "; hcalIso = " << matchPhoton2->hcalTowerSumEtConeDR04() << endl;
     // cout << "; hcalIsoCut = " << (2.2 + 0.0025*matchPhoton2->pt()+0.062*fRho25) << endl;
     // cout << "; pixelSeed = " << matchPhoton2->hasPixelSeed() << endl;
     // cout << "; sigmaIetaIeta = " << matchPhoton2->sigmaIetaIeta() << endl;
     // cout << "; isGap = " <<  ExoDiPhotons::isGapPhoton(&(*matchPhoton2)) << endl;
     // cout << "; isSpike = " <<  ExoDiPhotons::isASpike(&(*matchPhoton2)) << endl; 
     // cout << "; passesHadOverEm = " << ExoDiPhotons::passesHadOverEmCut(matchPhoton2) << endl;
     // cout << ": passesTrkIso = " << ExoDiPhotons::passesTrkIsoCut(matchPhoton2,fRho25) << endl;
     // cout << "; passesEcalIsoCut = " << ExoDiPhotons::passesEcalIsoCut(matchPhoton2,fRho25) << endl;
     // cout << "; passesHcalIsoCut = " << ExoDiPhotons::passesHcalIsoCut(matchPhoton2,fRho25) << endl;
     // cout << "; passesSigmaIetaIetaCut = " << ExoDiPhotons::passesSigmaIetaIetaCut(matchPhoton2) << endl;
     // cout << endl << endl;
     //}


     // fill info into tree
     ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo2,matchPhoton2,lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);
     // SIC add all this to the function above?
     fRecoPhotonInfo2.hasMatchedPromptElec = ConversionTools::hasMatchedPromptElectron(matchPhoton2->superCluster(), hElectrons, hConversions, beamSpot.position());
     //we retrieve the effective areas
     //Remember effareaCH = 1st, effareaNH = 2nd, effareaPH = 3rd
     //Now we store all PF isolation variables for the 2nd photon
     std::vector<double> photon2TorFEffAreas = ExoDiPhotons::EffectiveAreas(matchPhoton2);
     double pfisoall = isolator03.fGetIsolation(matchPhoton2,pfCandidates,firstVtx,vertexColl);
     double rhocorPFIsoCH = max(isolator03.getIsolationCharged()-fRho25*photon2TorFEffAreas[0],0.);
     double rhocorPFIsoNH = max(isolator03.getIsolationNeutral()-fRho25*photon2TorFEffAreas[1],0.);
     double rhocorPFIsoPH = max(isolator03.getIsolationPhoton()-fRho25*photon2TorFEffAreas[2],0.);
     
     fRecoPhotonInfo2.PFIsoAll04 = isolator04.fGetIsolation(matchPhoton2,pfCandidates,firstVtx,vertexColl);
     fRecoPhotonInfo2.PFIsoCharged04 = isolator04.getIsolationCharged();
     fRecoPhotonInfo2.PFIsoNeutral04 = isolator04.getIsolationNeutral();
     fRecoPhotonInfo2.PFIsoPhoton04 = isolator04.getIsolationPhoton();      

     fRecoPhotonInfo2.PFIsoAll03 = pfisoall;
     fRecoPhotonInfo2.PFIsoCharged03 = isolator03.getIsolationCharged();
     fRecoPhotonInfo2.PFIsoNeutral03 = isolator03.getIsolationNeutral();
     fRecoPhotonInfo2.PFIsoPhoton03 = isolator03.getIsolationPhoton();      
    
     fRecoPhotonInfo2.PFIsoAll02 = isolator02.fGetIsolation(matchPhoton2,pfCandidates,firstVtx,vertexColl);
     fRecoPhotonInfo2.PFIsoCharged02 = isolator02.getIsolationCharged();
     fRecoPhotonInfo2.PFIsoNeutral02 = isolator02.getIsolationNeutral();
     fRecoPhotonInfo2.PFIsoPhoton02 = isolator02.getIsolationPhoton();      

     //now the corrected PF isolation variables
     fRecoPhotonInfo2.rhocorPFIsoCharged04 = max(fRecoPhotonInfo2.PFIsoCharged04-fRho25*photon2TorFEffAreas[0],0.);
     fRecoPhotonInfo2.rhocorPFIsoNeutral04 = max(fRecoPhotonInfo2.PFIsoNeutral04-fRho25*photon2TorFEffAreas[1],0.);
     fRecoPhotonInfo2.rhocorPFIsoPhoton04 = max(fRecoPhotonInfo2.PFIsoPhoton04-fRho25*photon2TorFEffAreas[2],0.);
     fRecoPhotonInfo2.rhocorPFIsoAll04 = fRecoPhotonInfo2.rhocorPFIsoCharged04 + fRecoPhotonInfo2.rhocorPFIsoNeutral04 + fRecoPhotonInfo2.rhocorPFIsoPhoton04;

     fRecoPhotonInfo2.rhocorPFIsoCharged03 = rhocorPFIsoCH;
     fRecoPhotonInfo2.rhocorPFIsoNeutral03 = rhocorPFIsoNH;
     fRecoPhotonInfo2.rhocorPFIsoPhoton03 = rhocorPFIsoPH;
     fRecoPhotonInfo2.rhocorPFIsoAll03 = fRecoPhotonInfo2.rhocorPFIsoCharged03 + fRecoPhotonInfo2.rhocorPFIsoNeutral03 + fRecoPhotonInfo2.rhocorPFIsoPhoton03;

     fRecoPhotonInfo2.rhocorPFIsoCharged02 = max(fRecoPhotonInfo2.PFIsoCharged02-fRho25*photon2TorFEffAreas[0],0.);
     fRecoPhotonInfo2.rhocorPFIsoNeutral02 = max(fRecoPhotonInfo2.PFIsoNeutral02-fRho25*photon2TorFEffAreas[1],0.);
     fRecoPhotonInfo2.rhocorPFIsoPhoton02 = max(fRecoPhotonInfo2.PFIsoPhoton02-fRho25*photon2TorFEffAreas[2],0.);
     fRecoPhotonInfo2.rhocorPFIsoAll02 = fRecoPhotonInfo2.rhocorPFIsoCharged02 + fRecoPhotonInfo2.rhocorPFIsoNeutral02 + fRecoPhotonInfo2.rhocorPFIsoPhoton02;

     // now we fill the ID variables for det and PFID
     if(ExoDiPhotons::isTightPhoton(matchPhoton2,fRho25) &&
        !ExoDiPhotons::isGapPhoton(matchPhoton2) &&
        !ExoDiPhotons::isASpike(matchPhoton2)  )
       fRecoPhotonInfo2.isTightDetPhoton = true;
     else
       fRecoPhotonInfo2.isTightDetPhoton = false;
    
     // test the conversion safe electron veto
     bool passElecVeto = !ConversionTools::hasMatchedPromptElectron(matchPhoton2->superCluster(), hElectrons, hConversions, beamSpot.position());
     if(ExoDiPhotons::isPFTightPhoton(&(*matchPhoton2),rhocorPFIsoCH,rhocorPFIsoNH,rhocorPFIsoPH,passElecVeto,TString("Tight")) && 
	!ExoDiPhotons::isGapPhoton(&(*matchPhoton2)) && 
	!ExoDiPhotons::isASpike(&(*matchPhoton2))  )
       fRecoPhotonInfo2.isTightPFPhoton = true;
     else
       fRecoPhotonInfo2.isTightPFPhoton = false;
     if(ExoDiPhotons::isPFTightPhoton(&(*matchPhoton2),rhocorPFIsoCH,rhocorPFIsoNH,rhocorPFIsoPH,passElecVeto,TString("Medium")) && 
	!ExoDiPhotons::isGapPhoton(&(*matchPhoton2)) && 
	!ExoDiPhotons::isASpike(&(*matchPhoton2))  )
       fRecoPhotonInfo2.isMediumPFPhoton = true;
     else
       fRecoPhotonInfo2.isMediumPFPhoton = false;
     if(ExoDiPhotons::isPFTightPhoton(&(*matchPhoton2),rhocorPFIsoCH,rhocorPFIsoNH,rhocorPFIsoPH,passElecVeto,TString("Loose")) && 
	!ExoDiPhotons::isGapPhoton(&(*matchPhoton2)) && 
	!ExoDiPhotons::isASpike(&(*matchPhoton2))  )
       fRecoPhotonInfo2.isLoosePFPhoton = true;
     else
       fRecoPhotonInfo2.isLoosePFPhoton = false;
   }
   else {
     //     cout << "No match to signal photon2!" <<endl;
     // as short cut for indicating this in tree
     // make sure the pt value is crazy
     fRecoPhotonInfo2.pt = -9999.99;
   }

   
   //   cout << endl;
   //   cout << endl;


   // put in diphoton info?
   // (both for signal and reco photons?)

   // fill diphoton info

   // start with silly values
   fDiphotonSignalInfo.Minv = -9999.99;
   fDiphotonRecoInfo.Minv = -9999.99;
   
   if(signalPhoton1&&signalPhoton2) {
     ExoDiPhotons::FillDiphotonInfo(fDiphotonSignalInfo,signalPhoton1,signalPhoton2);
   }

   if(matchPhoton1&&matchPhoton2) {
     ExoDiPhotons::FillDiphotonInfo(fDiphotonRecoInfo,matchPhoton1,matchPhoton2);
   }

   */

   // for this signal MC, want to fill the tree every event
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
