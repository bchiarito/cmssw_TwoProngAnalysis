// -*- C++ -*-
//
// Package:    ExoDiPhotonWithGenAnalyzer
// Class:      ExoDiPhotonWithGenAnalyzer
// 
/**\class ExoDiPhotonWithGenAnalyzer ExoDiPhotonWithGenAnalyzer.cc DiPhotonAnalysis/ExoDiPhotonWithGenAnalyzer/src/ExoDiPhotonWithGenAnalyzer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Conor Henderson,40 1-B01,+41227671674,
//         Created:  Thu May  6 17:26:16 CEST 2010
// $Id: ExoDiPhotonWithGenAnalyzer.cc,v 1.32 2013/02/11 15:07:42 charaf Exp $
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
#include "L1Trigger/GlobalTrigger/interface/L1GlobalTrigger.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

// new CommonClasses approach
// these objects are all in the namespace 'ExoDiPhotons'
#include "DiPhotonAnalysis/CommonClasses/interface/RecoPhotonInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/TriggerInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/PhotonID.h"
#include "DiPhotonAnalysis/CommonClasses/interface/EventAndVertexInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/DiphotonInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/MCTrueObjectInfo.h"


//new for PU gen
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// new for LumiReweighting 
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

//new for PF ID definition
#include "DiPhotonAnalysis/CommonClasses/interface/PFPhotonID.h"

//for conversion safe electron veto
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//new for PFIsolation code
#include "EGamma/EGammaAnalysisTools/src/PFIsolationEstimator.cc"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

using namespace std;


//
// class declaration
//


class ExoDiPhotonWithGenAnalyzer : public edm::EDAnalyzer {
public:
  explicit ExoDiPhotonWithGenAnalyzer(const edm::ParameterSet&);
  ~ExoDiPhotonWithGenAnalyzer();


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
  std::auto_ptr<EcalClusterLazyTools> lazyTools_;



  // to get L1 info, the L1 guide recommends to make this a member
  // this allows the event setup parts to be cached, rather than refetched every event
  L1GtUtils m_l1GtUtils;

  // my Tree
  TTree *fTree;

  //gen info
  ExoDiPhotons::mcTrueObjectInfo_t fSignalPhoton1Info; // leading signal photon
  ExoDiPhotons::mcTrueObjectInfo_t fSignalPhoton2Info;
  ExoDiPhotons::diphotonInfo_t fSignalDiphotonInfo;

  ExoDiPhotons::mcTrueObjectInfo_t fSignalUnstablePhoton1Info; // leading signal photon
  ExoDiPhotons::mcTrueObjectInfo_t fSignalUnstablePhoton2Info;
  ExoDiPhotons::diphotonInfo_t fSignalUnstableDiphotonInfo;

  float deltarPhoton1;
  float deltarPhoton2;
  float ptratioPhoton1;
  float ptratioPhoton2;
  //end of gen info

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

  //for PFIsolation Code
  PFIsolationEstimator isolator04;
  PFIsolationEstimator isolator03;
  PFIsolationEstimator isolator02;

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
ExoDiPhotonWithGenAnalyzer::ExoDiPhotonWithGenAnalyzer(const edm::ParameterSet& iConfig)
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
    fIDMethod(iConfig.getUntrackedParameter<string>("IDMethod"))
{
  //now do what ever initialization is needed
  // LumiReweighting Tool

  std::cout<<"ExoDiPhotonWithGenAnalyzer: ID Method used "<<fIDMethod.c_str()
	   <<"PF ID Category "<<fPFIDCategory.c_str()
	   <<std::endl;

  edm::Service<TFileService> fs;
  fpu_n_BeforeCuts = fs->make<TH1F>("fpu_n_BeforeCuts","PileUpBeforeCuts",300,0,300);
  fpu_n_BeforeCutsAfterReWeight = fs->make<TH1F>("fpu_n_BeforeCutsAfterReWeight","PileUpBeforeCuts",300,0,300);
  fNumTotalEvents = fs->make<TH1F>("NumTotalEvents","Total number of events",4,0.,2.);
  fNumTotalWeightedEvents = fs->make<TH1F>("NumTotalWeightedEvents","Total weighted number of events",4,0.,2.);
  fTree = fs->make<TTree>("fTree","PhotonTree");

  //gen info
  deltarPhoton1 = -20.;
  deltarPhoton2 = -20.;
  ptratioPhoton1 = -20.;
  ptratioPhoton2 = -20.;

  fTree->Branch("GenPhoton1",&fSignalPhoton1Info,ExoDiPhotons::mcTrueObjectInfoBranchDefString.c_str());
  fTree->Branch("GenPhoton2",&fSignalPhoton2Info,ExoDiPhotons::mcTrueObjectInfoBranchDefString.c_str());
  fTree->Branch("GenDiphoton",&fSignalDiphotonInfo,ExoDiPhotons::diphotonInfoBranchDefString.c_str());

  fTree->Branch("GenUnstablePhoton1",&fSignalUnstablePhoton1Info,ExoDiPhotons::mcTrueObjectInfoBranchDefString.c_str());
  fTree->Branch("GenUnstablePhoton2",&fSignalUnstablePhoton2Info,ExoDiPhotons::mcTrueObjectInfoBranchDefString.c_str());
  fTree->Branch("GenUnstableDiphoton",&fSignalUnstableDiphotonInfo,ExoDiPhotons::diphotonInfoBranchDefString.c_str());

  //   fTree->Branch("Deltar1",&deltarPhoton1,"Deltar1/F");
  //   fTree->Branch("Deltar2",&deltarPhoton2,"Deltar2/F");
  
  //   fTree->Branch("Ptratio1",&ptratioPhoton1,"Ptratio1/F");
  //   fTree->Branch("Ptratio2",&ptratioPhoton2,"Ptratio2/F");

  //end of gen info

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
  
  gv_pos = new TClonesArray("TVector3", 100);
  gv_p3 = new TClonesArray("TVector3", 100);

  //new PFIsolation code
  isolator04.initializePhotonIsolation(kTRUE);
  isolator04.setConeSize(0.4);
  isolator03.initializePhotonIsolation(kTRUE);
  isolator03.setConeSize(0.3);
  isolator02.initializePhotonIsolation(kTRUE);
  isolator02.setConeSize(0.2);

}


ExoDiPhotonWithGenAnalyzer::~ExoDiPhotonWithGenAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//



// ------------ method called to for each event  ------------
void
ExoDiPhotonWithGenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  //   cout <<  iEvent.id().run() << " " <<  iEvent.id().luminosityBlock() << " " << iEvent.id().event() << endl;

  // basic event info
  ExoDiPhotons::InitEventInfo(fEventInfo,-5000.);
  ExoDiPhotons::FillEventInfo(fEventInfo,iEvent);
   
  edm::Handle<GenEventInfoProduct> GenInfoHandle;
  if(fkRequireGenEventInfo) iEvent.getByLabel("generator",GenInfoHandle);

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

  // get the vertex collection
  Handle<reco::VertexCollection> vertexColl;
  iEvent.getByLabel("offlinePrimaryVertices",vertexColl);
   
  if(!vertexColl.isValid()) {
    cout << "Vertex collection empty! Bailing out!" <<endl;
    return;
  }
  //    //    cout << "N vertices = " << vertexColl->size() <<endl;
  //    //   fVtxInfo.Nvtx = vertexColl->size();
  //    // this just counts the collection size
  //    // may want to count N vtx with TrkPt> some cut ?

  
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


  //    // note for higher lumi, may want to also store second vertex, for pileup studies
  //    // to allow scalability for many vertices, use a vector and sort later
  std::vector<reco::Vertex> myVertices;
   
  for(reco::VertexCollection::const_iterator vtx=vertexColl->begin(); vtx!=vertexColl->end(); vtx++) {
    if(!vtx->isFake() && vtx->ndof()>4 && fabs(vtx->position().rho())<=2.0 && fabs(vtx->z())<=24.0  ) {
      myVertices.push_back(*vtx);
    }         
  }// end vertex loop
  sort(myVertices.begin(),myVertices.end(),ExoDiPhotons::sortVertices) ;

  //    // first count the number of good vertices
  fVtxInfo.Nvtx = myVertices.size();
  fVtx2Info.Nvtx = myVertices.size();
  fVtx3Info.Nvtx = myVertices.size();

  if(myVertices.size()>=1) {
    ExoDiPhotons::FillVertexInfo(fVtxInfo,&(*myVertices.begin()));
  }
  if(myVertices.size()>=2) {
    ExoDiPhotons::FillVertexInfo(fVtx2Info,&(*(myVertices.begin()+1)));
  }

  if(myVertices.size()>=3) {
    ExoDiPhotons::FillVertexInfo(fVtx3Info,&(*(myVertices.begin()+2)));
  }

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
  edm::Handle<double> rho25Handle;
  iEvent.getByLabel(fRho25Tag, rho25Handle);

  if (!rho25Handle.isValid()){
    cout<<"rho25 not found"<<endl;
    return;
  }
        
  fRho25 = *(rho25Handle.product());

  //for conversion safe electron veto
  edm::Handle<reco::ConversionCollection> hConversions;
  iEvent.getByLabel("allConversions", hConversions);
  
  edm::Handle<reco::GsfElectronCollection> hElectrons;
  //iEvent.getByLabel("gsfElectrons", hElectrons);
  iEvent.getByLabel("gedGsfElectrons", hElectrons);
  if(!hElectrons.isValid()) {
    cout<<"no ged gsf electrons "<<endl;
    return;
  }


  //for PFIsolation code
  Handle<PFCandidateCollection> pfCandidatesColl;
  iEvent.getByLabel("particleFlow",pfCandidatesColl);
  const PFCandidateCollection * pfCandidates = pfCandidatesColl.product();

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

  if(beamSpotHandle.isValid()) {
    beamSpot = *beamSpotHandle.product();
    //   cout << beamSpot <<endl;
    ExoDiPhotons::FillBeamSpotInfo(fBeamSpotInfo,beamSpot);
  }

  //trig results
  Handle<TriggerResults> hltResultsHandle;
  iEvent.getByLabel(fHltInputTag,hltResultsHandle);
   
  if(!hltResultsHandle.isValid()) {
    cout << "HLT results not valid!" <<endl;
    cout << "Couldnt find TriggerResults with input tag " << fHltInputTag << endl;
    return;
  }

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
   
  m_l1GtUtils.retrieveL1EventSetup(iSetup);

  int iErrorCode = -1;
   
  fL1TrigInfo.L1_Tech0 =    m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BPTX_plus_AND_minus.v0",iErrorCode);
  fL1TrigInfo.L1_Tech36 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_halo_beam2_inner.v0",iErrorCode);
  fL1TrigInfo.L1_Tech37 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_halo_beam2_outer.v0",iErrorCode);
  fL1TrigInfo.L1_Tech38 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_halo_beam1_inner.v0",iErrorCode);
  fL1TrigInfo.L1_Tech39 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_halo_beam1_outer.v0",iErrorCode);
  fL1TrigInfo.L1_Tech40 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_minBias_threshold1.v0",iErrorCode);
  fL1TrigInfo.L1_Tech41 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_minBias_threshold2.v0",iErrorCode);
  fL1TrigInfo.L1_Tech42 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_splash_beam1.v0",iErrorCode);
  fL1TrigInfo.L1_Tech43 =   m_l1GtUtils.decisionBeforeMask(iEvent,"L1Tech_BSC_splash_beam2.v0",iErrorCode);


  lazyTools_ = std::auto_ptr<EcalClusterLazyTools>( new 
						    EcalClusterLazyTools(iEvent,iSetup,edm::InputTag("reducedEcalRecHitsEB"),edm::InputTag("reducedEcalRecHitsEE")) 
						    );

  // get ecal barrel recHits for spike rejection
  edm::Handle<EcalRecHitCollection> recHitsEB_h;
  iEvent.getByLabel(edm::InputTag("reducedEcalRecHitsEB"), recHitsEB_h );
  const EcalRecHitCollection * recHitsEB = 0;
  if ( ! recHitsEB_h.isValid() ) {
    LogError("ExoDiPhotonWithGenAnalyzer") << " ECAL Barrel RecHit Collection not available !"; return;
  } else {
    recHitsEB = recHitsEB_h.product();
  }

  edm::Handle<EcalRecHitCollection> recHitsEE_h;
  iEvent.getByLabel(edm::InputTag("reducedEcalRecHitsEE"), recHitsEE_h );
  const EcalRecHitCollection * recHitsEE = 0;
  if ( ! recHitsEE_h.isValid() ) {
    LogError("ExoDiPhotonWithGenAnalyzer") << " ECAL Endcap RecHit Collection not available !"; return;
  } else {
    recHitsEE = recHitsEE_h.product();
  }

  edm::ESHandle<EcalChannelStatus> chStatus;
  iSetup.get<EcalChannelStatusRcd>().get(chStatus);
  const EcalChannelStatus *ch_status = chStatus.product(); 

  //get the reference to 1st vertex for use in fGetIsolation
  //for PFIsolation calculation
  reco::VertexRef firstVtx(vertexColl,0);
   
  //for gen info
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
  fSignalPhoton1Info.isol04 = -999999.99;
  fSignalPhoton1Info.isol04ratio = -999999.99;
  fSignalPhoton1Info.isol03 = -999999.99;
  fSignalPhoton1Info.isol03ratio = -999999.99;
  fSignalPhoton1Info.isol02 = -999999.99;
  fSignalPhoton1Info.isol02ratio = -999999.99;

  fSignalPhoton2Info.status = -999999;
  fSignalPhoton2Info.PdgId = -999999;
  fSignalPhoton2Info.MotherPdgId = -999999;
  fSignalPhoton2Info.GrandmotherPdgId = -999999;
  fSignalPhoton2Info.pt = -999999.99;
  fSignalPhoton2Info.eta = -999999.99;
  fSignalPhoton2Info.phi = -999999.99;
  fSignalPhoton2Info.isol04 = -999999.99;
  fSignalPhoton2Info.isol04ratio = -999999.99;
  fSignalPhoton2Info.isol03 = -999999.99;
  fSignalPhoton2Info.isol03ratio = -999999.99;
  fSignalPhoton2Info.isol02 = -999999.99;
  fSignalPhoton2Info.isol02ratio = -999999.99;

  fSignalUnstablePhoton1Info.status = -999999;
  fSignalUnstablePhoton1Info.PdgId = -999999;
  fSignalUnstablePhoton1Info.MotherPdgId = -999999;
  fSignalUnstablePhoton1Info.GrandmotherPdgId = -999999;
  fSignalUnstablePhoton1Info.pt = -999999.99;
  fSignalUnstablePhoton1Info.eta = -999999.99;
  fSignalUnstablePhoton1Info.phi = -999999.99;


  fSignalUnstablePhoton2Info.status = -999999;
  fSignalUnstablePhoton2Info.PdgId = -999999;
  fSignalUnstablePhoton2Info.MotherPdgId = -999999;
  fSignalUnstablePhoton2Info.GrandmotherPdgId = -999999;
  fSignalUnstablePhoton2Info.pt = -999999.99;
  fSignalUnstablePhoton2Info.eta = -999999.99;
  fSignalUnstablePhoton2Info.phi = -999999.99;

  const reco::Candidate *signalPhoton1 = NULL;
  const reco::Candidate *signalPhoton2 = NULL;

  const reco::GenParticle *signalUnstablePhoton1 = NULL;
  const reco::GenParticle *signalUnstablePhoton2 = NULL;

  std::vector<std::pair<reco::GenParticle, const reco::Candidate*> > genPhotonPairs;

  int indexgenpart = -1;
  int indexgenphoton = -1;

  for(reco::GenParticleCollection::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); ++genParticle) {

    indexgenpart++;

    if(genParticle->pdgId() != 22) continue;
    if(genParticle->status() != 3) continue;
    //if(genParticle->pt() <= 10.) continue;

    int indexDaughter = -1;
    int indexGrandDaughter = -1;
    float deltarMinDaughter = 1000.;

    if(genParticle->numberOfDaughters() == 0) continue;
    
    indexgenphoton++;
    for(unsigned int j=0;j<genParticle->numberOfDaughters();j++){
      if(genParticle->daughter(j)->pdgId() != 22) continue;
      for(unsigned int k=0;k<genParticle->daughter(j)->numberOfDaughters();k++){
	if(genParticle->daughter(j)->daughter(k)->pdgId() != 22) continue;	  
	float deltar = deltaR(genParticle->eta(),genParticle->phi(),genParticle->daughter(j)->daughter(k)->eta(),genParticle->daughter(j)->daughter(k)->phi());
	if(deltar < deltarMinDaughter){
	  indexDaughter = j;
	  indexGrandDaughter = k;
	  deltarMinDaughter = deltar;
	}	
      }//end of loop over grand daughters

      if(genParticle->daughter(j)->status() == 11) continue;
      float deltar = deltaR(genParticle->eta(),genParticle->phi(),genParticle->daughter(j)->eta(),genParticle->daughter(j)->phi());
      if(deltar < deltarMinDaughter){
	indexDaughter = j;
	deltarMinDaughter = deltar;
      }

    }//end of loop over daughters      
    
    const reco::Candidate *CorrectDaughter;

    //this is the case where an intermediate status 11 photon was found
    //so from status 3 to status 11 to status 1
    if(indexDaughter >=0 && indexGrandDaughter >= 0) CorrectDaughter = genParticle->daughter(indexDaughter)->daughter(indexGrandDaughter);

    //this is the case where no intermediate status 11 photon was found
    //so directly from status 3 to status 1
    if(indexDaughter >=0 && indexGrandDaughter < 0) CorrectDaughter = genParticle->daughter(indexDaughter);

    //this is the case where no branching is found
    //the photon produced (status 3) is the final one
    //we set CorrectDaughter to be the same as the photon
    if(indexDaughter < 0 && indexGrandDaughter < 0) CorrectDaughter = genParticle->daughter(0)->mother(indexgenphoton);

    std::pair<reco::GenParticle, const reco::Candidate*> genpair = std::make_pair(*genParticle, CorrectDaughter);
    genPhotonPairs.push_back(genpair);

    bool motherInHardProcess = true;
    //now we check that the mother is a parton only in the hard process
    //i.e. all its mother particles have pt = 0
    if(genParticle->numberOfMothers() == 0) continue;
    for(unsigned int l=0;l<genParticle->numberOfMothers();l++){
      if(genParticle->mother(l)->pt() != 0.) motherInHardProcess = false;
    }//end of loop over mothers
    cout<<""<<endl;

    if(!motherInHardProcess) {cout<<"not from hard process"<<endl;};

  }// ends gen particle loop
  
  cout<<"number of pairs found "<<genPhotonPairs.size()<<endl;

  if(genPhotonPairs.size() < 2) {cout<<"here problem"<<endl;return;}

  deltarPhoton1 = deltaR(genPhotonPairs[0].first.eta(),genPhotonPairs[0].first.phi(),genPhotonPairs[0].second->eta(),genPhotonPairs[0].second->phi());
  deltarPhoton2 = deltaR(genPhotonPairs[1].first.eta(),genPhotonPairs[1].first.phi(),genPhotonPairs[1].second->eta(),genPhotonPairs[1].second->phi());

  ptratioPhoton1 = genPhotonPairs[0].first.pt()/genPhotonPairs[0].second->pt();
  ptratioPhoton2 = genPhotonPairs[1].first.pt()/genPhotonPairs[1].second->pt();

  sort(genPhotonPairs.begin(),genPhotonPairs.end(),ExoDiPhotons::compareGenPhotonPairsByPt);

  const reco::Candidate *genStablePhoton1 = genPhotonPairs[0].second;
  const reco::Candidate *genStablePhoton2 = genPhotonPairs[1].second;

  fSignalPhoton1Info.isol04 = 0.;
  fSignalPhoton2Info.isol04 = 0.;
  fSignalPhoton1Info.isol04ratio = 0.;
  fSignalPhoton2Info.isol04ratio = 0.;

  fSignalPhoton1Info.isol03 = 0.;
  fSignalPhoton2Info.isol03 = 0.;
  fSignalPhoton1Info.isol03ratio = 0.;
  fSignalPhoton2Info.isol03ratio = 0.;

  fSignalPhoton1Info.isol02 = 0.;
  fSignalPhoton2Info.isol02 = 0.;
  fSignalPhoton1Info.isol02ratio = 0.;
  fSignalPhoton2Info.isol02ratio = 0.;

  //calculate the isolation
  for(reco::GenParticleCollection::const_iterator genParticle = genParticles->begin(); 
      genParticle != genParticles->end();++genParticle)
    {
      if(genParticle->status() != 1) continue;

      //around photon 1
      float deltar1 = deltaR(genStablePhoton1->eta(),genStablePhoton1->phi(),genParticle->eta(),genParticle->phi());
      if(deltar1 < 0.4 && deltar1 > 0.001){
	fSignalPhoton1Info.isol04 += genParticle->et();


// 	//--------------------TESTING ZONE----------------------------

// 	cout<< "number of pairs found "<<genPhotonPairs.size()<<endl;
	
// 	cout<< "MC particle (isol part): Status = "<< genParticle->status()
// 	    << "; pdg id = "<< genParticle->pdgId() 
// 	    << "; (px,py,pz,E) = ("<< genParticle->px() 
// 	    << "," << genParticle->py() 
// 	    << "," << genParticle->pz() 
// 	    << "," << genParticle->energy() << ")"
// 	    << "; pt, eta, phi = " << genParticle->pt() 
// 	    << ", "<< genParticle->eta() 
// 	    << ", " << genParticle->phi() 
// 	    << endl;

// 	for(unsigned int j=0;j<genParticle->numberOfMothers();j++){
// 	  const reco::Candidate *partmother = genParticle->mother(j);

// 	  cout << "   --> and its mother "<<j<<" : Status = "<< partmother->status() 
// 	       << "; pdg id = "<< partmother->pdgId() 
// 	       << "; (px,py,pz,E) = ("<< partmother->px() 
// 	       << "," << partmother->py() 
// 	       << "," << partmother->pz() 
// 	       << "," << partmother->energy() << ")"
// 	       << "; pt, eta, phi = " << partmother->pt() 
// 	       << ", "<< partmother->eta() 
// 	       << ", " << partmother->phi() 
// 	       << endl;



// 	  for(unsigned int k=0;k<genParticle->numberOfMothers();k++){
// 	    const reco::Candidate *partgrandmother = genParticle->mother(j)->mother(k);

// 	    cout << "   --> and its mother "<<k<<" : Status = "<< partgrandmother->status() 
// 		 << "; pdg id = "<< partgrandmother->pdgId() 
// 		 << "; (px,py,pz,E) = ("<< partgrandmother->px() 
// 		 << "," << partgrandmother->py() 
// 		 << "," << partgrandmother->pz() 
// 		 << "," << partgrandmother->energy() << ")"
// 		 << "; pt, eta, phi = " << partgrandmother->pt() 
// 		 << ", "<< partgrandmother->eta() 
// 		 << ", " << partgrandmother->phi() 
// 		 << endl;

// 	  }






// 	  //--------------------END OF TESTING ZONE----------------------------

	}
	if(deltar1 < 0.3 && deltar1 > 0.001){
	  fSignalPhoton1Info.isol03 += genParticle->et();
	}
	if(deltar1 < 0.2 && deltar1 > 0.001){
	  fSignalPhoton1Info.isol02 += genParticle->et();
	}

	//around photon 2
	float deltar2 = deltaR(genStablePhoton2->eta(),genStablePhoton2->phi(),genParticle->eta(),genParticle->phi());
	if(deltar2 < 0.4 && deltar2 > 0.001){
	  fSignalPhoton2Info.isol04 += genParticle->et();
	}
	if(deltar2 < 0.3 && deltar2 > 0.001){
	  fSignalPhoton2Info.isol03 += genParticle->et();
	}
	if(deltar2 < 0.2 && deltar2 > 0.001){
	  fSignalPhoton2Info.isol02 += genParticle->et();
	}
      }

      fSignalPhoton1Info.isol04ratio = fSignalPhoton1Info.isol04/genStablePhoton1->pt();;
      fSignalPhoton2Info.isol04ratio = fSignalPhoton2Info.isol04/genStablePhoton2->pt();;

      fSignalPhoton1Info.isol03ratio = fSignalPhoton1Info.isol03/genStablePhoton1->pt();;
      fSignalPhoton2Info.isol03ratio = fSignalPhoton2Info.isol03/genStablePhoton2->pt();;

      fSignalPhoton1Info.isol02ratio = fSignalPhoton1Info.isol02/genStablePhoton1->pt();;
      fSignalPhoton2Info.isol02ratio = fSignalPhoton2Info.isol02/genStablePhoton2->pt();;
  
      signalPhoton1 = genStablePhoton1;
      signalPhoton2 = genStablePhoton2;
  
      signalUnstablePhoton1 = &(genPhotonPairs[0].first);
      signalUnstablePhoton2 = &(genPhotonPairs[1].first);

      ExoDiPhotons::FillMCTrueObjectInfo(fSignalPhoton1Info,signalPhoton1);
      ExoDiPhotons::FillMCTrueObjectInfo(fSignalPhoton2Info,signalPhoton2);
      ExoDiPhotons::FillDiphotonInfo(fSignalDiphotonInfo,signalPhoton1,signalPhoton2);

      ExoDiPhotons::FillMCTrueObjectInfo(fSignalUnstablePhoton1Info,signalUnstablePhoton1);
      ExoDiPhotons::FillMCTrueObjectInfo(fSignalUnstablePhoton2Info,signalUnstablePhoton2);
      ExoDiPhotons::FillDiphotonInfo(fSignalUnstableDiphotonInfo,signalUnstablePhoton1,signalUnstablePhoton2);

      //end of gen info

      // get the photon collection
      Handle<reco::PhotonCollection> photonColl;
      iEvent.getByLabel(fPhotonTag,photonColl);

      // If photon collection is empty, exit
      if (!photonColl.isValid()) {
	cout << "No Photons! Move along, there's nothing to see here .." <<endl;
	return;
      }

      TString CategoryPFID(fPFIDCategory.c_str());
      TString MethodID(fIDMethod.c_str());

      std::vector<reco::Photon> selectedPhotons; 

      // photon loop
      for(reco::PhotonCollection::const_iterator recoPhoton = photonColl->begin(); recoPhoton!=photonColl->end(); recoPhoton++) {
	/*
	  cout << "Photon et, eta, phi = " << recoPhoton->et() <<", "<<recoPhoton->eta()<< ", "<< recoPhoton->phi();
	  cout << "; calo position eta = " << recoPhoton->caloPosition().eta();
	  cout << "; eMax/e3x3 = " << recoPhoton->maxEnergyXtal()/recoPhoton->e3x3();
	  cout << "; hadOverEm = " << recoPhoton->hadronicOverEm();
	  cout << "; trkIso = " << recoPhoton->trkSumPtHollowConeDR04();
	  cout << "; ecalIso = " << recoPhoton->ecalRecHitSumEtConeDR04();
	  cout << "; hcalIso = " << recoPhoton->hcalTowerSumEtConeDR04();
	  cout << "; pixelSeed = " << recoPhoton->hasPixelSeed();
	  cout << "; sigmaietaieta = " << recoPhoton->sigmaIetaIeta();      cout << endl;
	*/

	//we retrieve the effective areas
	//Remember effareaCH = 1st, effareaNH = 2nd, effareaPH = 3rd
	std::vector<double> effareas = ExoDiPhotons::EffectiveAreas(&(*recoPhoton));
	double pfisoall = isolator03.fGetIsolation(&(*recoPhoton),pfCandidates,firstVtx,vertexColl);
	double rhocorPFIsoCH = max(isolator03.getIsolationCharged()-fRho25*effareas[0],0.);
	double rhocorPFIsoNH = max(isolator03.getIsolationNeutral()-fRho25*effareas[1],0.);
	double rhocorPFIsoPH = max(isolator03.getIsolationPhoton()-fRho25*effareas[2],0.);
	//and we also have to test the conversion safe electron veto

	bool passelecveto = !ConversionTools::hasMatchedPromptElectron(recoPhoton->superCluster(), hElectrons, hConversions, beamSpot.position());

	if(ExoDiPhotons::isBarrelPhoton(&(*recoPhoton)) && (recoPhoton->pt()>=fMin_pt)) {

	  //Now we choose which ID to use (PF or Det)
	  if(MethodID.Contains("Detector")){
	    if(ExoDiPhotons::isTightPhoton(&(*recoPhoton),fRho25) && !ExoDiPhotons::isGapPhoton(&(*recoPhoton)) && !ExoDiPhotons::isASpike(&(*recoPhoton))  ) {
	      selectedPhotons.push_back(*recoPhoton);
	    }
	  }
	  else if(MethodID.Contains("ParticleFlow")){
	    if(ExoDiPhotons::isPFTightPhoton(&(*recoPhoton),rhocorPFIsoCH,rhocorPFIsoNH,rhocorPFIsoPH,passelecveto,CategoryPFID) && 
	       !ExoDiPhotons::isGapPhoton(&(*recoPhoton)) && 
	       !ExoDiPhotons::isASpike(&(*recoPhoton))  ) {
	      selectedPhotons.push_back(*recoPhoton);
	    }
	  }
      
	} //end first cuts on pt and (not applied, april2011) EB-only
       
      } //end reco photon loop

      sort(selectedPhotons.begin(),selectedPhotons.end(),ExoDiPhotons::comparePhotonsByPt);
      fNTightPhotons = selectedPhotons.size();

      ExoDiPhotons::InitRecoPhotonInfo(fRecoPhotonInfo1);
      ExoDiPhotons::InitRecoPhotonInfo(fRecoPhotonInfo2);

      if ( selectedPhotons.size() >=2 ){

	// must specifically declare isFakeable status (should be Tight = not True = false                       

	ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo1,&selectedPhotons[0],lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);
	fRecoPhotonInfo1.hasMatchedPromptElec = ConversionTools::hasMatchedPromptElectron((&selectedPhotons[0])->superCluster(), hElectrons, hConversions, beamSpot.position());
	fRecoPhotonInfo1.isFakeable = false;

	//Now we store all PF isolation variables for the 1st photon (tight exception)
	std::vector<double> photon1TEffAreas = ExoDiPhotons::EffectiveAreas((&selectedPhotons[0]));
     
	fRecoPhotonInfo1.PFIsoAll04 = isolator04.fGetIsolation((&selectedPhotons[0]),pfCandidates,firstVtx,vertexColl);
	fRecoPhotonInfo1.PFIsoCharged04 = isolator04.getIsolationCharged();
	fRecoPhotonInfo1.PFIsoNeutral04 = isolator04.getIsolationNeutral();
	fRecoPhotonInfo1.PFIsoPhoton04 = isolator04.getIsolationPhoton();      

	fRecoPhotonInfo1.PFIsoAll03 = isolator03.fGetIsolation((&selectedPhotons[0]),pfCandidates,firstVtx,vertexColl);
	fRecoPhotonInfo1.PFIsoCharged03 = isolator03.getIsolationCharged();
	fRecoPhotonInfo1.PFIsoNeutral03 = isolator03.getIsolationNeutral();
	fRecoPhotonInfo1.PFIsoPhoton03 = isolator03.getIsolationPhoton();      
    
	fRecoPhotonInfo1.PFIsoAll02 = isolator02.fGetIsolation((&selectedPhotons[0]),pfCandidates,firstVtx,vertexColl);
	fRecoPhotonInfo1.PFIsoCharged02 = isolator02.getIsolationCharged();
	fRecoPhotonInfo1.PFIsoNeutral02 = isolator02.getIsolationNeutral();
	fRecoPhotonInfo1.PFIsoPhoton02 = isolator02.getIsolationPhoton();      


	//now the corrected PF isolation variables
	fRecoPhotonInfo1.rhocorPFIsoCharged04 = max(fRecoPhotonInfo1.PFIsoCharged04-fRho25*photon1TEffAreas[0],0.);
	fRecoPhotonInfo1.rhocorPFIsoNeutral04 = max(fRecoPhotonInfo1.PFIsoNeutral04-fRho25*photon1TEffAreas[1],0.);
	fRecoPhotonInfo1.rhocorPFIsoPhoton04 = max(fRecoPhotonInfo1.PFIsoPhoton04-fRho25*photon1TEffAreas[2],0.);
	fRecoPhotonInfo1.rhocorPFIsoAll04 = fRecoPhotonInfo1.rhocorPFIsoCharged04 + fRecoPhotonInfo1.rhocorPFIsoNeutral04 + fRecoPhotonInfo1.rhocorPFIsoPhoton04;

	fRecoPhotonInfo1.rhocorPFIsoCharged03 = max(fRecoPhotonInfo1.PFIsoCharged03-fRho25*photon1TEffAreas[0],0.);
	fRecoPhotonInfo1.rhocorPFIsoNeutral03 = max(fRecoPhotonInfo1.PFIsoNeutral03-fRho25*photon1TEffAreas[1],0.);
	fRecoPhotonInfo1.rhocorPFIsoPhoton03 = max(fRecoPhotonInfo1.PFIsoPhoton03-fRho25*photon1TEffAreas[2],0.);
	fRecoPhotonInfo1.rhocorPFIsoAll03 = fRecoPhotonInfo1.rhocorPFIsoCharged03 + fRecoPhotonInfo1.rhocorPFIsoNeutral03 + fRecoPhotonInfo1.rhocorPFIsoPhoton03;

	fRecoPhotonInfo1.rhocorPFIsoCharged02 = max(fRecoPhotonInfo1.PFIsoCharged02-fRho25*photon1TEffAreas[0],0.);
	fRecoPhotonInfo1.rhocorPFIsoNeutral02 = max(fRecoPhotonInfo1.PFIsoNeutral02-fRho25*photon1TEffAreas[1],0.);
	fRecoPhotonInfo1.rhocorPFIsoPhoton02 = max(fRecoPhotonInfo1.PFIsoPhoton02-fRho25*photon1TEffAreas[2],0.);
	fRecoPhotonInfo1.rhocorPFIsoAll02 = fRecoPhotonInfo1.rhocorPFIsoCharged02 + fRecoPhotonInfo1.rhocorPFIsoNeutral02 + fRecoPhotonInfo1.rhocorPFIsoPhoton02;


	ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo2,&selectedPhotons[1],lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);
	fRecoPhotonInfo2.hasMatchedPromptElec = ConversionTools::hasMatchedPromptElectron((&selectedPhotons[1])->superCluster(), hElectrons, hConversions, beamSpot.position());
	fRecoPhotonInfo2.isFakeable = false;

	//Now we store all PF isolation variables for the 2st photon (tight exception)
	std::vector<double> photon2TEffAreas = ExoDiPhotons::EffectiveAreas((&selectedPhotons[1]));
     
	fRecoPhotonInfo2.PFIsoAll04 = isolator04.fGetIsolation((&selectedPhotons[1]),pfCandidates,firstVtx,vertexColl);
	fRecoPhotonInfo2.PFIsoCharged04 = isolator04.getIsolationCharged();
	fRecoPhotonInfo2.PFIsoNeutral04 = isolator04.getIsolationNeutral();
	fRecoPhotonInfo2.PFIsoPhoton04 = isolator04.getIsolationPhoton();      

	fRecoPhotonInfo2.PFIsoAll03 = isolator03.fGetIsolation((&selectedPhotons[1]),pfCandidates,firstVtx,vertexColl);
	fRecoPhotonInfo2.PFIsoCharged03 = isolator03.getIsolationCharged();
	fRecoPhotonInfo2.PFIsoNeutral03 = isolator03.getIsolationNeutral();
	fRecoPhotonInfo2.PFIsoPhoton03 = isolator03.getIsolationPhoton();      
    
	fRecoPhotonInfo2.PFIsoAll02 = isolator02.fGetIsolation((&selectedPhotons[1]),pfCandidates,firstVtx,vertexColl);
	fRecoPhotonInfo2.PFIsoCharged02 = isolator02.getIsolationCharged();
	fRecoPhotonInfo2.PFIsoNeutral02 = isolator02.getIsolationNeutral();
	fRecoPhotonInfo2.PFIsoPhoton02 = isolator02.getIsolationPhoton();      

	//now the corrected PF isolation variables
	fRecoPhotonInfo2.rhocorPFIsoCharged04 = max(fRecoPhotonInfo2.PFIsoCharged04-fRho25*photon2TEffAreas[0],0.);
	fRecoPhotonInfo2.rhocorPFIsoNeutral04 = max(fRecoPhotonInfo2.PFIsoNeutral04-fRho25*photon2TEffAreas[1],0.);
	fRecoPhotonInfo2.rhocorPFIsoPhoton04 = max(fRecoPhotonInfo2.PFIsoPhoton04-fRho25*photon2TEffAreas[2],0.);
	fRecoPhotonInfo2.rhocorPFIsoAll04 = fRecoPhotonInfo2.rhocorPFIsoCharged04 + fRecoPhotonInfo2.rhocorPFIsoNeutral04 + fRecoPhotonInfo2.rhocorPFIsoPhoton04;

	fRecoPhotonInfo2.rhocorPFIsoCharged03 = max(fRecoPhotonInfo2.PFIsoCharged03-fRho25*photon2TEffAreas[0],0.);
	fRecoPhotonInfo2.rhocorPFIsoNeutral03 = max(fRecoPhotonInfo2.PFIsoNeutral03-fRho25*photon2TEffAreas[1],0.);
	fRecoPhotonInfo2.rhocorPFIsoPhoton03 = max(fRecoPhotonInfo2.PFIsoPhoton03-fRho25*photon2TEffAreas[2],0.);
	fRecoPhotonInfo2.rhocorPFIsoAll03 = fRecoPhotonInfo2.rhocorPFIsoCharged03 + fRecoPhotonInfo2.rhocorPFIsoNeutral03 + fRecoPhotonInfo2.rhocorPFIsoPhoton03;

	fRecoPhotonInfo2.rhocorPFIsoCharged02 = max(fRecoPhotonInfo2.PFIsoCharged02-fRho25*photon2TEffAreas[0],0.);
	fRecoPhotonInfo2.rhocorPFIsoNeutral02 = max(fRecoPhotonInfo2.PFIsoNeutral02-fRho25*photon2TEffAreas[1],0.);
	fRecoPhotonInfo2.rhocorPFIsoPhoton02 = max(fRecoPhotonInfo2.PFIsoPhoton02-fRho25*photon2TEffAreas[2],0.);
	fRecoPhotonInfo2.rhocorPFIsoAll02 = fRecoPhotonInfo2.rhocorPFIsoCharged02 + fRecoPhotonInfo2.rhocorPFIsoNeutral02 + fRecoPhotonInfo2.rhocorPFIsoPhoton02;

	// fill diphoton info                                                                                                   
	ExoDiPhotons::FillDiphotonInfo(fDiphotonInfo,&selectedPhotons[0],&selectedPhotons[1]);

      }//end of 2 TT exception

      fTree->Fill();

    }//end of method


  // ------------ method called once each job just before starting event loop  ------------
  void 
    ExoDiPhotonWithGenAnalyzer::beginJob()
  {
    if (fisMC){
      LumiWeights = edm::LumiReWeighting(fPUMCFileName,fPUDataFileName,fPUMCHistName,fPUDataHistName);
    }
  }

  // ------------ method called once each job just after ending the event loop  ------------
  void 
    ExoDiPhotonWithGenAnalyzer::endJob() {
  }



  //define this as a plug-in
  DEFINE_FWK_MODULE(ExoDiPhotonWithGenAnalyzer);
