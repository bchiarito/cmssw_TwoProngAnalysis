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
// $Id: ExoDiPhotonAnalyzer.cc,v 1.9 2010/08/24 20:29:13 chenders Exp $
//
//


// system include files
#include <memory>

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
  bool isASpike(const reco::Photon *photon);

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

      ExoDiPhotons::eventInfo_t fEventInfo;
      ExoDiPhotons::vtxInfo_t fVtxInfo;
      ExoDiPhotons::beamSpotInfo_t fBeamSpotInfo;

      ExoDiPhotons::l1TrigInfo_t fL1TrigInfo;  
      ExoDiPhotons::hltTrigInfo_t fHLTInfo;

      ExoDiPhotons::recoPhotonInfo_t fRecoPhotonInfo1; // leading photon 
      ExoDiPhotons::recoPhotonInfo_t fRecoPhotonInfo2; // second photon

      ExoDiPhotons::diphotonInfo_t fDiphotonInfo;

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
  fTree->Branch("BeamSpot",&fBeamSpotInfo,ExoDiPhotons::beamSpotInfoBranchDefString.c_str());


  fTree->Branch("L1trg",&fL1TrigInfo,ExoDiPhotons::l1TrigBranchDefString.c_str());
  
  fTree->Branch("TrigHLT",&fHLTInfo,ExoDiPhotons::hltTrigBranchDefString.c_str());
  
  // now with CommonClasses, use the string defined in the header
  fTree->Branch("Photon1",&fRecoPhotonInfo1,ExoDiPhotons::recoPhotonBranchDefString.c_str());

  fTree->Branch("Photon2",&fRecoPhotonInfo2,ExoDiPhotons::recoPhotonBranchDefString.c_str());

  fTree->Branch("Diphoton",&fDiphotonInfo,ExoDiPhotons::diphotonInfoBranchDefString.c_str());
}


ExoDiPhotonAnalyzer::~ExoDiPhotonAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

bool ExoDiPhotonAnalyzer::isASpike(const reco::Photon *photon)
{
  bool thisPhotonIsASpike = false;

  // for now, its simpler for me to use the older e1/e9 criteria to identify spikes
  /// but we need to move to the swiss cross cut ASAP
  
  if(photon->maxEnergyXtal()/photon->e3x3()>0.95)
    // then its a spike!
    thisPhotonIsASpike = true;

  return thisPhotonIsASpike;
}



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
   fVtxInfo.Nvtx = vertexColl->size();

   fVtxInfo.vx = -99999.99;
   fVtxInfo.vy = -99999.99;
   fVtxInfo.vz = -99999.99;
   fVtxInfo.isFake = -99;   
   fVtxInfo.Ntracks = -99;
   fVtxInfo.sumPtTracks = -99999.99;
   fVtxInfo.ndof = -99999.99;
   fVtxInfo.d0 = -99999.99;

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
   
   //   cout << "N photons = " << photonColl->size() <<endl;

   // we want the two highest Et photons for this analysis
   const reco::Photon *photon1 = NULL;
   const reco::Photon *photon2 = NULL;
   double highestEt1 = fMin_pt;
   double highestEt2 = fMin_pt;



   // photon loop
   for(reco::PhotonCollection::const_iterator recoPhoton = photonColl->begin(); recoPhoton!=photonColl->end(); recoPhoton++) {

//      cout << "Photon et, eta, phi = " << recoPhoton->et() <<", "<<recoPhoton->eta()<< ", "<< recoPhoton->phi();
//      cout << "; eMax/e3x3 = " << recoPhoton->maxEnergyXtal()/recoPhoton->e3x3();
//      cout << "; hadOverEm = " << recoPhoton->hadronicOverEm();
//      cout << "; trkIso = " << recoPhoton->trkSumPtHollowConeDR04();
//      cout << "; ecalIso = " << recoPhoton->ecalRecHitSumEtConeDR04();
//      cout << "; hcalIso = " << recoPhoton->hcalTowerSumEtConeDR04();
//      cout << "; pixelSeed = " << recoPhoton->hasPixelSeed();
//      cout << endl;

     // option to remove spikes, so only consider pt-ordering of non-spike photons
     // note that I have to do '&(*photon)' because I am using the iterator, 
     // while my function needs a pointer to the object (obtained by dereffing the iter) ...
     if(!fkRemoveSpikes || !isASpike(&(*recoPhoton))) {
       // ie either we dont care about spikes OR this photon is not a spike anyway
       // so in either case we continue to study this photon... 
       
       // similar logic for tight photon ID
       if(!fkRequireTightPhotons || ExoDiPhotons::isTightPhoton(&(*recoPhoton))) {
	 ///ie either we dont require tight photons OR this photon passes tight ID anyway

	 if(recoPhoton->et() >= fMin_pt && recoPhoton->et()>highestEt1) {
	   // formerly highest is now second highest
	   photon2 = photon1; // this can still be NULL at this point ...
	   if(photon2)
	     highestEt2 = photon2->et();
	   // and the new one is highest
	   photon1 =  &(*recoPhoton);
	   highestEt1 = photon1->et();
	   
	 }
	 else if(recoPhoton->et() >= fMin_pt && recoPhoton->et()>highestEt2) {
	   photon2 =  &(*recoPhoton);
	   highestEt2 = photon2->et();
	 }
	 
       } //ends tight photon ID check
     }//ends spike criteria check
       
   } //end reco photn loop

//    // now we have the two highest photons
    if(photon1) {
      //      cout << "Highest Et photon: et, eta, phi = " << photon1->et() <<", "<<photon1->eta()<< ", "<< photon1->phi()<<endl;


     // can now use the Fill function specific to this recoPhoton struct
      ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo1,photon1,lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);


    }
    if(photon2) {

      ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo2,photon2,lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);
      
      //      cout << "2nd Highest Et photon: et, eta, phi = " << photon2->et() <<", "<<photon2->eta()<< ", "<< photon2->phi()<<endl;


    }


    if(photon1&&photon2) {

      // fill diphoton info
      ExoDiPhotons::FillDiphotonInfo(fDiphotonInfo,photon1,photon2);

      // only fill the tree if there are two photons!
      fTree->Fill();
      
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
