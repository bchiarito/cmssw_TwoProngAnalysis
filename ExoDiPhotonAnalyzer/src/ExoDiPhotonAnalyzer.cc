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
// $Id: ExoDiPhotonAnalyzer.cc,v 1.18 2010/11/22 17:24:04 chenders Exp $
//
//


// system include files
#include <memory>

#include <algorithm>
#include <vector>


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
      bool               fkRemoveSpikes;   // option to remove spikes before filling tree
      bool               fkRequireTightPhotons;  // option to require tight photon id in tree

      // tools for clusters
      std::auto_ptr<EcalClusterLazyTools> lazyTools_;

      // to get L1 info, the L1 guide recommends to make this a member
      // this allows the event setup parts to be cached, rather than refetched every event
      L1GtUtils m_l1GtUtils;

      // my Tree
      TTree *fTree;

  // now also tight-fake and fake-fake trees
      TTree *fTightFakeTree;
      TTree *fFakeFakeTree;

      ExoDiPhotons::eventInfo_t fEventInfo;
      ExoDiPhotons::vtxInfo_t fVtxInfo;
  // adding a second vertex
      ExoDiPhotons::vtxInfo_t fVtx2Info;
  // now even adding a third vtx!
      ExoDiPhotons::vtxInfo_t fVtx3Info;

      ExoDiPhotons::beamSpotInfo_t fBeamSpotInfo;

      ExoDiPhotons::l1TrigInfo_t fL1TrigInfo;  
      ExoDiPhotons::hltTrigInfo_t fHLTInfo;

      int fNTightPhotons; // number of candidate photons in event (ie tight)
      int fNFakeablePhotons;  // number of 'fakeable objects' in event

      ExoDiPhotons::recoPhotonInfo_t fRecoPhotonInfo1; // leading photon 
      ExoDiPhotons::recoPhotonInfo_t fRecoPhotonInfo2; // second photon

      ExoDiPhotons::diphotonInfo_t fDiphotonInfo;

  // diphoton info based on using hte second or third vtx in event
      ExoDiPhotons::diphotonInfo_t fDiphotonInfoVtx2; 
      ExoDiPhotons::diphotonInfo_t fDiphotonInfoVtx3; 

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
    fkRemoveSpikes(iConfig.getUntrackedParameter<bool>("removeSpikes")),
    fkRequireTightPhotons(iConfig.getUntrackedParameter<bool>("requireTightPhotons"))
{
   //now do what ever initialization is needed

  edm::Service<TFileService> fs;
  fTree = fs->make<TTree>("fTree","PhotonTree");

  fTree->Branch("Event",&fEventInfo,ExoDiPhotons::eventInfoBranchDefString.c_str());
  fTree->Branch("Vtx",&fVtxInfo,ExoDiPhotons::vtxInfoBranchDefString.c_str());

  //adding a second vtx
  fTree->Branch("Vtx2",&fVtx2Info,ExoDiPhotons::vtxInfoBranchDefString.c_str());
  fTree->Branch("Vtx3",&fVtx3Info,ExoDiPhotons::vtxInfoBranchDefString.c_str());

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


   // basic event info
   ExoDiPhotons::FillEventInfo(fEventInfo,iEvent);
   
   // get the vertex collection
   Handle<reco::VertexCollection> vertexColl;
   iEvent.getByLabel("offlinePrimaryVertices",vertexColl);

   if(!vertexColl.isValid()) {
     cout << "Vertex collection empty! Bailing out!" <<endl;
     return;
   }
   //   cout << "N vertices = " << vertexColl->size() <<endl;
   //   fVtxInfo.Nvtx = vertexColl->size();
   // this just counts the collection size
   // may want to count N vtx with TrkPt> some cut ?

   

   fVtxInfo.vx = -99999.99;
   fVtxInfo.vy = -99999.99;
   fVtxInfo.vz = -99999.99;
   fVtxInfo.isFake = -99;   
   fVtxInfo.Ntracks = -99;
   fVtxInfo.sumPtTracks = -99999.99;
   fVtxInfo.ndof = -99999.99;
   fVtxInfo.d0 = -99999.99;


   fVtx2Info.vx = -99999.99;
   fVtx2Info.vy = -99999.99;
   fVtx2Info.vz = -99999.99;
   fVtx2Info.isFake = -99;   
   fVtx2Info.Ntracks = -99;
   fVtx2Info.sumPtTracks = -99999.99;
   fVtx2Info.ndof = -99999.99;
   fVtx2Info.d0 = -99999.99;


   fVtx3Info.vx = -99999.99;
   fVtx3Info.vy = -99999.99;
   fVtx3Info.vz = -99999.99;
   fVtx3Info.isFake = -99;   
   fVtx3Info.Ntracks = -99;
   fVtx3Info.sumPtTracks = -99999.99;
   fVtx3Info.ndof = -99999.99;
   fVtx3Info.d0 = -99999.99;


   // note for higher lumi, may want to also store second vertex, for pileup studies
   // to allow scalability for many vertices, use a vector and sort later
   std::vector<reco::Vertex> myVertices;
   
   for(reco::VertexCollection::const_iterator vtx=vertexColl->begin(); vtx!=vertexColl->end(); vtx++) {
     
     // add to my vtx vector if not fake and ndof>4 and maxd0=2 and |vz|<24
     // ie default criteria
     if(!vtx->isFake() && vtx->ndof()>4 && fabs(vtx->position().rho())<=2.0 && fabs(vtx->z())<=24.0  ) {
       myVertices.push_back(*vtx);
     }      

     //     cout << "Vtx x = "<< vtx->x()<<", y= "<< vtx->y()<<", z = " << vtx->z() << ";  N tracks = " << vtx->tracksSize() << "; isFake = " << vtx->isFake() <<", sumPt(tracks) = "<< ExoDiPhotons::calcVtxSumPtTracks(&(*vtx)) << "; ndof = " << vtx->ndof()<< "; d0 = " << vtx->position().rho() << endl;
     
     // and note that this vertex collection can contain vertices with Ntracks = 0
     // watch out for these!
   
   }// end vertex loop


//    // loop over my selected vertices
//    for(std::vector<reco::Vertex>::iterator myVtxIter=myVertices.begin();myVtxIter<myVertices.end();myVtxIter++) {

//      cout << "MY MyVtxIter x = "<< myVtxIter->x()<<", y= "<< myVtxIter->y()<<", z = " << myVtxIter->z() << ";  N tracks = " << myVtxIter->tracksSize() << "; isFake = " << myVtxIter->isFake() <<", sumPt(tracks) = "<< ExoDiPhotons::calcVtxSumPtTracks(&(*myVtxIter))<< "; ndof = " << myVtxIter->ndof()<< "; d0 = " << myVtxIter->position().rho() << endl;

//    }

   // now sort the vertices after
   // can be either by Ntracks or TrackSumPt, depending how I write the function
   sort(myVertices.begin(),myVertices.end(),ExoDiPhotons::sortVertices) ;

   //   cout << "After sorting" << endl;

//    // check sorting
//    for(std::vector<reco::Vertex>::iterator myVtxIter=myVertices.begin();myVtxIter<myVertices.end();myVtxIter++) {

//      cout << "MY MyVtxIter x = "<< myVtxIter->x()<<", y= "<< myVtxIter->y()<<", z = " << myVtxIter->z() << ";  N tracks = " << myVtxIter->tracksSize() << "; isFake = " << myVtxIter->isFake() <<", sumPt(tracks) = "<< ExoDiPhotons::calcVtxSumPtTracks(&(*myVtxIter))<< "; ndof = " << myVtxIter->ndof()<< "; d0 = " << myVtxIter->position().rho() << endl;

//    }


   // first count the number of good vertices
   fVtxInfo.Nvtx = myVertices.size();
   fVtx2Info.Nvtx = myVertices.size();

   // now we will fill the vertex info structs from the sorted list
   if(myVertices.size()>=1) {
     ExoDiPhotons::FillVertexInfo(fVtxInfo,&(*myVertices.begin()));
   }
   if(myVertices.size()>=2) {
     ExoDiPhotons::FillVertexInfo(fVtx2Info,&(*(myVertices.begin()+1)));
   }

   if(myVertices.size()>=3) {
     ExoDiPhotons::FillVertexInfo(fVtx3Info,&(*(myVertices.begin()+2)));
   }
   




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
     beamSpot = *beamSpotHandle;
     //   cout << beamSpot <<endl;
     ExoDiPhotons::FillBeamSpotInfo(fBeamSpotInfo,beamSpot);
   }



   // get the trig info

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


   // ecal information
   
   lazyTools_ = std::auto_ptr<EcalClusterLazyTools>( new EcalClusterLazyTools(iEvent,iSetup,edm::InputTag("ecalRecHit:EcalRecHitsEB"),edm::InputTag("ecalRecHit:EcalRecHitsEE")) );

   // get ecal barrel recHits for spike rejection
   edm::Handle<EcalRecHitCollection> recHitsEB_h;
   iEvent.getByLabel(edm::InputTag("ecalRecHit:EcalRecHitsEB"), recHitsEB_h );
   const EcalRecHitCollection * recHitsEB = 0;
   if ( ! recHitsEB_h.isValid() ) {
     LogError("ExoDiPhotonAnalyzer") << " ECAL Barrel RecHit Collection not available !"; return;
   } else {
     recHitsEB = recHitsEB_h.product();
   }

   edm::Handle<EcalRecHitCollection> recHitsEE_h;
   iEvent.getByLabel(edm::InputTag("ecalRecHit:EcalRecHitsEE"), recHitsEE_h );
   const EcalRecHitCollection * recHitsEE = 0;
   if ( ! recHitsEE_h.isValid() ) {
     LogError("ExoDiPhotonAnalyzer") << " ECAL Endcap RecHit Collection not available !"; return;
   } else {
     recHitsEE = recHitsEE_h.product();
   }

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
   
   //   cout << "N reco photons = " << photonColl->size() <<endl;

   // new approach - 
   // make vector of all selected Photons (tight ID, not spike, min pt, etc)
   // then sort at end by pt
   // also allows to count how often a third photon could be considered a candidate

   std::vector<reco::Photon> selectedPhotons; 

   // do a similar thing for 'Fakeable objects', for data-based fake rate approach
   std::vector<reco::Photon> fakeablePhotons; 


   // photon loop
   for(reco::PhotonCollection::const_iterator recoPhoton = photonColl->begin(); recoPhoton!=photonColl->end(); recoPhoton++) {

     //          cout << "Photon et, eta, phi = " << recoPhoton->et() <<", "<<recoPhoton->eta()<< ", "<< recoPhoton->phi();
     //     cout << "; calo position eta = " << recoPhoton->caloPosition().eta();
//      cout << "; eMax/e3x3 = " << recoPhoton->maxEnergyXtal()/recoPhoton->e3x3();
//       cout << "; hadOverEm = " << recoPhoton->hadronicOverEm();
//      cout << "; trkIso = " << recoPhoton->trkSumPtHollowConeDR04();
//       cout << "; ecalIso = " << recoPhoton->ecalRecHitSumEtConeDR04();
//       cout << "; hcalIso = " << recoPhoton->hcalTowerSumEtConeDR04();
//       cout << "; pixelSeed = " << recoPhoton->hasPixelSeed();
//       cout << "; sigmaietaieta = " << recoPhoton->sigmaIetaIeta();
//      cout << endl;


     // now add selected photons to vector if:
     // tight ID
     // not in gap
     // not a spike
     // min pt cut  (for all photons)
     if(ExoDiPhotons::isTightPhoton(&(*recoPhoton)) && !ExoDiPhotons::isGapPhoton(&(*recoPhoton)) && !ExoDiPhotons::isASpike(&(*recoPhoton)) && (recoPhoton->pt()>=fMin_pt) ) {
       
       selectedPhotons.push_back(*recoPhoton);
     }

     // also check for fakeable objects
     if(ExoDiPhotons::isFakeableObject(&(*recoPhoton)) && !ExoDiPhotons::isGapPhoton(&(*recoPhoton)) && (recoPhoton->pt()>=fMin_pt) ) {
       
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


   // now consider possible Fakeable objects

   // first sort by pt
   sort(fakeablePhotons.begin(),fakeablePhotons.end(),ExoDiPhotons::comparePhotonsByPt);

   //   cout << "N fakeable = " << fakeablePhotons.size() <<endl;
   fNFakeablePhotons = fakeablePhotons.size();



   if(selectedPhotons.size()>=2) {
     //     cout << "then we have two photons - fill them to tree, you fool, dont just hang around here doing nothing!" <<endl;

     // can now use the Fill function specific to this recoPhoton struct
     ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo1,&selectedPhotons[0],lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);

     // specifically declare these photons NOT fake
     fRecoPhotonInfo1.isFakeable = false;

     ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo2,&selectedPhotons[1],lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);
      
     fRecoPhotonInfo2.isFakeable = false;

     // fill diphoton info
     ExoDiPhotons::FillDiphotonInfo(fDiphotonInfo,&selectedPhotons[0],&selectedPhotons[1]);
     
     // for the pileup/vtx study.
     // lets try setting these photons to different vertices
     // and see how much these changes Mgg
     
     if(myVertices.size()>=2) {
       reco::Photon photon1_vtx2(selectedPhotons[0]);
       // keep calo poisiton same, but assign new vertex
       // this func also recalcs 4-vector after changing vertex!
       photon1_vtx2.setVertex(myVertices[1].position()); 

       reco::Photon photon2_vtx2(selectedPhotons[1]);
       photon2_vtx2.setVertex(myVertices[1].position()); 

       ExoDiPhotons::FillDiphotonInfo(fDiphotonInfoVtx2,&photon1_vtx2,&photon2_vtx2);

     }
     
     if(myVertices.size()>=3) {
       reco::Photon photon1_vtx3(selectedPhotons[0]);
       // keep calo poisiton same, but assign new vertex
       // this func also recalcs 4-vector after changing vertex!
       photon1_vtx3.setVertex(myVertices[2].position()); 

       reco::Photon photon2_vtx3(selectedPhotons[1]);
       photon2_vtx3.setVertex(myVertices[2].position()); 

       ExoDiPhotons::FillDiphotonInfo(fDiphotonInfoVtx3,&photon1_vtx3,&photon2_vtx3);

     }



     // only fill the tree if there are two photons!
     fTree->Fill();


   }

   else if(selectedPhotons.size()==1 && fakeablePhotons.size()>=1) {
     
     // then is tight+fake category
     // will always take tight, + highest pt fake

     //selectedPhotons[0]
     // fakeablePhotons[0]
     

     // order Photon1,2 by pt, no matter which is tight and which is fake
     // and then use isFakeable leaf to distinguish them
     
     if(selectedPhotons[0].pt()>=fakeablePhotons[0].pt()) {

       // then higher pt photon is tight
       ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo1,&selectedPhotons[0],lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);
       ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo2,&fakeablePhotons[0],lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);
       ExoDiPhotons::FillDiphotonInfo(fDiphotonInfo,&selectedPhotons[0],&fakeablePhotons[0]);
       
       // specifically declare tight/fake nature
       fRecoPhotonInfo1.isFakeable = false;
       fRecoPhotonInfo2.isFakeable = true;
       
     }
     else {

       // then higher pt photon is fake
       ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo1,&fakeablePhotons[0],lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);
       ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo2,&selectedPhotons[0],lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);

       ExoDiPhotons::FillDiphotonInfo(fDiphotonInfo,&fakeablePhotons[0],&selectedPhotons[0]);
       
       // specifically declare tight/fake nature
       fRecoPhotonInfo1.isFakeable = true;
       fRecoPhotonInfo2.isFakeable = false;       
       
     }
     
     // and fill the tight-fake tree instead!
     fTightFakeTree->Fill();

   }
   
   else if(selectedPhotons.size()==0 && fakeablePhotons.size()>=2) {

     // then this is fake-fake category

     ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo1,&fakeablePhotons[0],lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);
     ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo2,&fakeablePhotons[1],lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);
     ExoDiPhotons::FillDiphotonInfo(fDiphotonInfo,&fakeablePhotons[0],&fakeablePhotons[1]);

     // specifically need to declare these photons as fake
     fRecoPhotonInfo1.isFakeable = true;
     fRecoPhotonInfo2.isFakeable = true;

     // and fill the fake-fake tree instead!
     fFakeFakeTree->Fill();
     
   }



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
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ExoDiPhotonAnalyzer::endJob() {
}



//define this as a plug-in
DEFINE_FWK_MODULE(ExoDiPhotonAnalyzer);
