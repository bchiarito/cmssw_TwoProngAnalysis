// -*- C++ -*-
//
// Package:    ExoDiPhotonBkgAnalyzer
// Class:      ExoDiPhotonBkgAnalyzer
// 
/**\class ExoDiPhotonBkgAnalyzer ExoDiPhotonBkgAnalyzer.cc DiPhotonAnalysis/ExoDiPhotonBkgAnalyzer/src/ExoDiPhotonBkgAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Conor Henderson,40 1-B01,+41227671674,
//         Created:  Mon Jun 28 12:37:19 CEST 2010
// $Id: ExoDiPhotonBkgAnalyzer.cc,v 1.5 2010/08/18 11:18:46 torimoto Exp $
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

#include "FWCore/ParameterSet/interface/InputTag.h"


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


// for MC
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


// new CommonClasses approach
// these objects all in the namespace 'ExoDiPhotons'
#include "DiPhotonAnalysis/CommonClasses/interface/RecoPhotonInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/MCTrueObjectInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/PhotonID.h"
#include "DiPhotonAnalysis/CommonClasses/interface/EventAndVertexInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/DiphotonInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/TriggerInfo.h"


using namespace std;


//
// class declaration
//




class ExoDiPhotonBkgAnalyzer : public edm::EDAnalyzer {
   public:
      explicit ExoDiPhotonBkgAnalyzer(const edm::ParameterSet&);
      ~ExoDiPhotonBkgAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      Float_t getESRatio(const reco::Photon *photon, const edm::Event& e, const edm::EventSetup& iSetup);

      // ----------member data ---------------------------
  // input tags and parameters
      edm::InputTag      fPhotonTag;       //select photon collection 
      double             fMin_pt;          // min pt cut (photons)
      edm::InputTag      fHltInputTag;     // hltResults

      // tools for clusters
      std::auto_ptr<EcalClusterLazyTools> lazyTools_;

      // my Tree
      TTree *fTree;
  
      ExoDiPhotons::eventInfo_t fEventInfo;
      ExoDiPhotons::vtxInfo_t fVtxInfo;
      ExoDiPhotons::beamSpotInfo_t fBeamSpotInfo;

      ExoDiPhotons::hltTrigInfo_t fHLTInfo;

      ExoDiPhotons::recoPhotonInfo_t fRecoPhotonInfo1; // leading photon 
      ExoDiPhotons::recoPhotonInfo_t fRecoPhotonInfo2; // second photon

      ExoDiPhotons::mcTrueObjectInfo_t fMCTrueObjectInfo1_Status3; 
      ExoDiPhotons::mcTrueObjectInfo_t fMCTrueObjectInfo2_Status3; 
      ExoDiPhotons::mcTrueObjectInfo_t fMCTrueObjectInfo1_Status1; 
      ExoDiPhotons::mcTrueObjectInfo_t fMCTrueObjectInfo2_Status1; 

      ExoDiPhotons::diphotonInfo_t fDiphotonInfo;

  // event normalisation histogram
  TH1F *fNorm_h;

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
ExoDiPhotonBkgAnalyzer::ExoDiPhotonBkgAnalyzer(const edm::ParameterSet& iConfig)
  : fPhotonTag(iConfig.getUntrackedParameter<edm::InputTag>("photonCollection")),
    fMin_pt(iConfig.getUntrackedParameter<double>("ptMin")),
    // note that the HLT process name can vary for different MC samples
    // so be sure to adjsut correctly in cfg
    fHltInputTag(iConfig.getUntrackedParameter<edm::InputTag>("hltResults"))
  
{
   //now do what ever initialization is needed

  edm::Service<TFileService> fs;
  fTree = fs->make<TTree>("fTree","PhotonTree");
  
  // now with CommonClasses, use the string defined in the header

  fTree->Branch("Event",&fEventInfo,ExoDiPhotons::eventInfoBranchDefString.c_str());
  fTree->Branch("Vtx",&fVtxInfo,ExoDiPhotons::vtxInfoBranchDefString.c_str());
  fTree->Branch("BeamSpot",&fBeamSpotInfo,ExoDiPhotons::beamSpotInfoBranchDefString.c_str());
  fTree->Branch("TrigHLT",&fHLTInfo,ExoDiPhotons::hltTrigBranchDefString.c_str());
  fTree->Branch("Photon1",&fRecoPhotonInfo1,ExoDiPhotons::recoPhotonBranchDefString.c_str());

  fTree->Branch("Photon2",&fRecoPhotonInfo2,ExoDiPhotons::recoPhotonBranchDefString.c_str());

  fTree->Branch("MCMatchPhoton1_Status3",&fMCTrueObjectInfo1_Status3,ExoDiPhotons::mcTrueObjectInfoBranchDefString.c_str());
  fTree->Branch("MCMatchPhoton2_Status3",&fMCTrueObjectInfo2_Status3,ExoDiPhotons::mcTrueObjectInfoBranchDefString.c_str());
  fTree->Branch("MCMatchPhoton1_Status1",&fMCTrueObjectInfo1_Status1,ExoDiPhotons::mcTrueObjectInfoBranchDefString.c_str());
  fTree->Branch("MCMatchPhoton2_Status1",&fMCTrueObjectInfo2_Status1,ExoDiPhotons::mcTrueObjectInfoBranchDefString.c_str());

  fTree->Branch("Diphoton",&fDiphotonInfo,ExoDiPhotons::diphotonInfoBranchDefString.c_str());


  // and the histogram for event normalisation purposes
  // Fill this with 0 if event fails, 1 if event passes
  fNorm_h = fs->make<TH1F>("fNorm_h","Normalisation Hist",4,0,2);
  

}


ExoDiPhotonBkgAnalyzer::~ExoDiPhotonBkgAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
ExoDiPhotonBkgAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   // basic event info
   ExoDiPhotons::FillEventInfo(fEventInfo,iEvent);

   // vertex

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
//      fVtxInfo.vx = vertex1->x();
//      fVtxInfo.vy = vertex1->y();
//      fVtxInfo.vz = vertex1->z();
//      fVtxInfo.isFake = vertex1->isFake(); 
//      fVtxInfo.Ntracks = vertex1->tracksSize();

//      fVtxInfo.ndof = vertex1->ndof();
//      fVtxInfo.d0 = vertex1->position().rho();


     ExoDiPhotons::FillVertexInfo(fVtxInfo,vertex1);
     // fill the SumPt Tracks separately for now
     fVtxInfo.sumPtTracks = highestSumPtTracks1;
     
   }



   //beam spot

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

   // get the photon collection
   Handle<reco::PhotonCollection> photonColl;
   iEvent.getByLabel(fPhotonTag,photonColl);

   // If photon collection is empty, exit
   if (!photonColl.isValid()) {
     cout << "No Photons! Move along, there's nothing to see here .." <<endl;
     return;
   }
   
   //   cout << "N reco photons = " << photonColl->size() <<endl;

   // we want the two highest Et photons for this analysis
   const reco::Photon *photon1 = NULL;
   const reco::Photon *photon2 = NULL;
   double highestEt1 = fMin_pt;
   double highestEt2 = fMin_pt;


   // photon loop
   for(reco::PhotonCollection::const_iterator recoPhoton = photonColl->begin(); recoPhoton!=photonColl->end(); recoPhoton++) {

//      cout << "Reco photon et, eta, phi = " << recoPhoton->et() <<", "<<recoPhoton->eta()<< ", "<< recoPhoton->phi();
//      cout << "; eMax/e3x3 = " << recoPhoton->maxEnergyXtal()/recoPhoton->e3x3();
//      cout << "; hadOverEm = " << recoPhoton->hadronicOverEm();
//      cout << "; trkIso = " << recoPhoton->trkSumPtHollowConeDR04();
//      cout << "; ecalIso = " << recoPhoton->ecalRecHitSumEtConeDR04();
//      cout << "; hcalIso = " << recoPhoton->hcalTowerSumEtConeDR04();
//      cout << "; pixelSeed = " << recoPhoton->hasPixelSeed();
//      cout << endl;

     //     if(!fkRequireTightPhotons || isTightPhoton(&(*recoPhoton))) {
     if(ExoDiPhotons::isTightPhoton(&(*recoPhoton))) {


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

       
     } // end if tight photon
     
   } //end of reco photon loop


   // now we have the two highest recoPhotons which pass our cuts
   
   if(photon1) {

     // can now use the Fill function specific to this recoPhoton struct
     ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo1,photon1);

     // anything using lazy tools still needs to be filled here though,
     // cause I had trouble getting it to work inside the Fill function

     reco::SuperClusterRef sc1 = photon1->superCluster();
     std::pair<DetId,float> maxc1 = lazyTools_->getMaximum(*sc1);

     //detId and ieta/iphi/iY/iX
     fRecoPhotonInfo1.detId = maxc1.first.rawId();
     
     if(maxc1.first.subdetId() == EcalBarrel) {
       EBDetId ebId( maxc1.first );
       fRecoPhotonInfo1.iEtaY = ebId.ieta(); // iEta in EB
       fRecoPhotonInfo1.iPhiX = ebId.iphi(); // iPhi in EB
     }
     else if (maxc1.first.subdetId() == EcalEndcap) {
       EEDetId eeId( maxc1.first );
       fRecoPhotonInfo1.iEtaY = eeId.iy(); // iY in EE
       fRecoPhotonInfo1.iPhiX = eeId.ix(); // iX in EE       
     }

     // swiss cross and other spike-related info
     if (maxc1.first.subdetId() == EcalBarrel) 
       fRecoPhotonInfo1.swisscross = EcalSeverityLevelAlgo::swissCross(maxc1.first, (*recHitsEB), 0);
     fRecoPhotonInfo1.eMax = lazyTools_->eMax(*sc1);
     fRecoPhotonInfo1.eLeft = lazyTools_->eLeft(*sc1);
     fRecoPhotonInfo1.eRight = lazyTools_->eRight(*sc1);
     fRecoPhotonInfo1.eTop = lazyTools_->eTop(*sc1);
     fRecoPhotonInfo1.eBottom = lazyTools_->eBottom(*sc1);
     fRecoPhotonInfo1.eSecond = lazyTools_->e2nd(*sc1);

     const reco::CaloClusterPtr  seed = sc1->seed();

     DetId id = lazyTools_->getMaximum(*seed).first; 
     float time  = -999., outOfTimeChi2 = -999., chi2 = -999.;
     int   flags=-1, severity = -1; 

     const EcalRecHitCollection & rechits = ( photon1->isEB() ? *recHitsEB : *recHitsEE); 
     EcalRecHitCollection::const_iterator it = rechits.find( id );
     if( it != rechits.end() ) { 
       time = it->time(); 
       outOfTimeChi2 = it->outOfTimeChi2();
       chi2 = it->chi2();
       flags = it->recoFlag();
       severity = EcalSeverityLevelAlgo::severityLevel( id, rechits, *chStatus );
     }

     const EcalChannelStatus *ch_status = chStatus.product(); 
     EcalChannelStatusMap::const_iterator chit;
     chit = ch_status->getMap().find(id.rawId());
     int mystatus = -99;
     if( chit != ch_status->getMap().end() ){
       EcalChannelStatusCode ch_code = (*chit);
       mystatus = ch_code.getStatusCode();
     }

     //     cout << "Photon1 seed " << lazyTools_->getMaximum(*seed).second << " " << maxc1.second << " " << lazyTools_->getMaximum(*seed).second -  maxc1.second << endl;     
     //     cout << "Photon1 time chStatus flags severity: " << time << " " << mystatus << " " << flags << " " << severity << endl;

     fRecoPhotonInfo1.severityLevel = severity;
     fRecoPhotonInfo1.recHitFlag = flags;
     fRecoPhotonInfo1.maxRecHitTime = time;

     fRecoPhotonInfo1.esRatio = getESRatio(photon1, iEvent, iSetup);

   }

   if(photon2) {
     ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo2,photon2);

     // anything using lazy tools still needs to be filled here though,
     // cause I had trouble getting it to work inside the Fill function

     reco::SuperClusterRef sc2 = photon2->superCluster();
     std::pair<DetId,float> maxc2 = lazyTools_->getMaximum(*sc2);

     //detId and ieta/iphi/iY/iX
     fRecoPhotonInfo2.detId = maxc2.first.rawId();
     
     if(maxc2.first.subdetId() == EcalBarrel) {
       EBDetId ebId( maxc2.first );
       fRecoPhotonInfo2.iEtaY = ebId.ieta(); // iEta in EB
       fRecoPhotonInfo2.iPhiX = ebId.iphi(); // iPhi in EB
     }
     else if (maxc2.first.subdetId() == EcalEndcap) {
       EEDetId eeId( maxc2.first );
       fRecoPhotonInfo2.iEtaY = eeId.iy(); // iY in EE
       fRecoPhotonInfo2.iPhiX = eeId.ix(); // iX in EE       
     }

     // swiss cross and other spike-related info
     if (maxc2.first.subdetId() == EcalBarrel) 
       fRecoPhotonInfo2.swisscross = EcalSeverityLevelAlgo::swissCross(maxc2.first, (*recHitsEB), 0);
     fRecoPhotonInfo2.eMax = lazyTools_->eMax(*sc2);
     fRecoPhotonInfo2.eLeft = lazyTools_->eLeft(*sc2);
     fRecoPhotonInfo2.eRight = lazyTools_->eRight(*sc2);
     fRecoPhotonInfo2.eTop = lazyTools_->eTop(*sc2);
     fRecoPhotonInfo2.eBottom = lazyTools_->eBottom(*sc2);
     fRecoPhotonInfo2.eSecond = lazyTools_->e2nd(*sc2);


     const reco::CaloClusterPtr  seed = sc2->seed();

     DetId id = lazyTools_->getMaximum(*seed).first;
     float time  = -999., outOfTimeChi2 = -999., chi2 = -999.;
     int   flags=-1, severity = -1;

     const EcalRecHitCollection & rechits = ( photon2->isEB() ? *recHitsEB : *recHitsEE);
     EcalRecHitCollection::const_iterator it = rechits.find( id );
     if( it != rechits.end() ) {
       time = it->time();
       outOfTimeChi2 = it->outOfTimeChi2();
       chi2 = it->chi2();
       flags = it->recoFlag();
       severity = EcalSeverityLevelAlgo::severityLevel( id, rechits, *chStatus );
     }

     const EcalChannelStatus *ch_status = chStatus.product();
     EcalChannelStatusMap::const_iterator chit;
     chit = ch_status->getMap().find(id.rawId());
     int mystatus = -99;
     if( chit != ch_status->getMap().end() ){
       EcalChannelStatusCode ch_code = (*chit);
       mystatus = ch_code.getStatusCode();
     }

     //     cout << "Photon2 seed " << lazyTools_->getMaximum(*seed).second << " " << maxc2.second << " " << lazyTools_->getMaximum(*seed).second -  maxc2.second << endl;
     //     cout << "Photon2 time chStatus flags severity: " << time << " " << mystatus << " " << flags << " " << severity << endl;

     fRecoPhotonInfo2.severityLevel = severity;
     fRecoPhotonInfo2.recHitFlag = flags;
     fRecoPhotonInfo2.maxRecHitTime = time;

     fRecoPhotonInfo2.esRatio = getESRatio(photon2, iEvent, iSetup);

   }


   // require that we have two photons passing our cuts
   if(photon1&&photon2) {

     // print both photons
//      cout << "Reco photon 1:  et, eta, phi = " << photon1->et() <<", "<<photon1->eta()<< ", "<< photon1->phi();
//      cout << "; eMax/e3x3 = " << photon1->maxEnergyXtal()/photon1->e3x3();
//      cout << "; hadOverEm = " << photon1->hadronicOverEm();
//      cout << "; trkIso = " << photon1->trkSumPtHollowConeDR04();
//      cout << "; ecalIso = " << photon1->ecalRecHitSumEtConeDR04();
//      cout << "; hcalIso = " << photon1->hcalTowerSumEtConeDR04();
//      cout << "; pixelSeed = " << photon1->hasPixelSeed();
//      cout << endl;

//      cout << "Reco photon 2:  et, eta, phi = " << photon2->et() <<", "<<photon2->eta()<< ", "<< photon2->phi();
//      cout << "; eMax/e3x3 = " << photon2->maxEnergyXtal()/photon2->e3x3();
//      cout << "; hadOverEm = " << photon2->hadronicOverEm();
//      cout << "; trkIso = " << photon2->trkSumPtHollowConeDR04();
//      cout << "; ecalIso = " << photon2->ecalRecHitSumEtConeDR04();
//      cout << "; hcalIso = " << photon2->hcalTowerSumEtConeDR04();
//      cout << "; pixelSeed = " << photon2->hasPixelSeed();
//      cout << endl;




     // now get MC truth info
     Handle<reco::GenParticleCollection> genParticles;
     iEvent.getByLabel("genParticles",genParticles);

     if(!genParticles.isValid()) {
       cout << "No Gen Particles collection!" << endl;
       return;
     }

     // we're going to match best status3 and status 1 MC particles
     const reco::GenParticle *genMatchPhoton1_Status3 = NULL;
     const reco::GenParticle *genMatchPhoton2_Status3 = NULL;
     const reco::GenParticle *genMatchPhoton1_Status1 = NULL;
     const reco::GenParticle *genMatchPhoton2_Status1 = NULL;

     double minDeltaR_pho1_status3 = 0.2;
     double minDeltaR_pho2_status3 = 0.2;
     double minDeltaR_pho1_status1 = 0.2;
     double minDeltaR_pho2_status1 = 0.2;

     for(reco::GenParticleCollection::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); ++genParticle) {

       //      if(genParticle->pt()>10) {
       //        cout << "MC particle: Status = "<< genParticle->status() << "; pdg id = "<< genParticle->pdgId() << "; pt, eta, phi = " << genParticle->pt() << ", "<< genParticle->eta() << ", " << genParticle->phi();
       //        if(genParticle->numberOfMothers()>0) {
       // 	 cout << "; Mother pid = " << genParticle->mother()->pdgId();

       // 	 if(genParticle->mother()->numberOfMothers()>0) {
       // 	   cout << "; GrandMother pid = " << genParticle->mother()->mother()->pdgId() << endl;
	   
       // 	 }
       //        }
       //      }

     
     
       // there's a CMS function for deltaPhi in DataFormats/Math
       double deltaPhi1 = reco::deltaPhi(genParticle->phi(),photon1->phi());
       double deltaEta1 = genParticle->eta()-photon1->eta();
       double deltaR1 = TMath::Sqrt(deltaPhi1*deltaPhi1+deltaEta1*deltaEta1);
     
       double deltaPhi2 = reco::deltaPhi(genParticle->phi(),photon2->phi());
       double deltaEta2 = genParticle->eta()-photon2->eta();
       double deltaR2 = TMath::Sqrt(deltaPhi2*deltaPhi2+deltaEta2*deltaEta2);



       if(genParticle->status()==3) {

	 if(deltaR1<minDeltaR_pho1_status3) {
	   // then this is best match so far
	   minDeltaR_pho1_status3 = deltaR1;
	   genMatchPhoton1_Status3 = &(*genParticle);
	 
	 }

	 // we can also allow this same gen particle to be matched to either reco photon
	 if(deltaR2<minDeltaR_pho2_status3) {
	   // then this is best match so far
	   minDeltaR_pho2_status3 = deltaR2;
	   genMatchPhoton2_Status3 = &(*genParticle);
	 
	 }

       } // end if status 3


       else if(genParticle->status()==1) {

	 if(deltaR1<minDeltaR_pho1_status1) {
	   // then this is best match so far
	   minDeltaR_pho1_status1 = deltaR1;
	   genMatchPhoton1_Status1 = &(*genParticle);
	 
	 }

	 // we can also allow this same gen particle to be matched to either reco photon
	 if(deltaR2<minDeltaR_pho2_status1) {
	   // then this is best match so far
	   minDeltaR_pho2_status1 = deltaR2;
	   genMatchPhoton2_Status1 = &(*genParticle);
	 
	 }
       } // end if status 1

     } //end loop over gen particles


     // its still possible that the matched photons are NULL, of course, so be careful
     if(genMatchPhoton1_Status3) {
       ExoDiPhotons::FillMCTrueObjectInfo(fMCTrueObjectInfo1_Status3,genMatchPhoton1_Status3);
       // print
       //       cout << "MC particle: Status = "<< genMatchPhoton1_Status3->status() << "; pdg id = "<< genMatchPhoton1_Status3->pdgId() << "; pt, eta, phi = " << genMatchPhoton1_Status3->pt() << ", "<< genMatchPhoton1_Status3->eta() << ", " << genMatchPhoton1_Status3->phi() <<endl;
       
     }
     else {
       // fill the default with nonsense values
       fMCTrueObjectInfo1_Status3.status = -999;
     }

     if(genMatchPhoton2_Status3) {
       ExoDiPhotons::FillMCTrueObjectInfo(fMCTrueObjectInfo2_Status3,genMatchPhoton2_Status3);
       // print
       //       cout << "MC particle: Status = "<< genMatchPhoton2_Status3->status() << "; pdg id = "<< genMatchPhoton2_Status3->pdgId() << "; pt, eta, phi = " << genMatchPhoton2_Status3->pt() << ", "<< genMatchPhoton2_Status3->eta() << ", " << genMatchPhoton2_Status3->phi() <<endl;
     }
     else {
       // fill the default with nonsense values
       fMCTrueObjectInfo2_Status3.status = -999;
     }


     if(genMatchPhoton1_Status1) {
       ExoDiPhotons::FillMCTrueObjectInfo(fMCTrueObjectInfo1_Status1,genMatchPhoton1_Status1);
       // print
       //       cout << "MC particle: Status = "<< genMatchPhoton1_Status1->status() << "; pdg id = "<< genMatchPhoton1_Status1->pdgId() << "; pt, eta, phi = " << genMatchPhoton1_Status1->pt() << ", "<< genMatchPhoton1_Status1->eta() << ", " << genMatchPhoton1_Status1->phi() <<endl;
     }
     else {
       // fill the default with nonsense values
       fMCTrueObjectInfo1_Status1.status = -999;
     }

     if(genMatchPhoton2_Status1) {
       ExoDiPhotons::FillMCTrueObjectInfo(fMCTrueObjectInfo2_Status1,genMatchPhoton2_Status1);
       // print
       //       cout << "MC particle: Status = "<< genMatchPhoton2_Status1->status() << "; pdg id = "<< genMatchPhoton2_Status1->pdgId() << "; pt, eta, phi = " << genMatchPhoton2_Status1->pt() << ", "<< genMatchPhoton2_Status1->eta() << ", " << genMatchPhoton2_Status1->phi() <<endl;
     }
     else {
       // fill the default with nonsense values
       fMCTrueObjectInfo2_Status1.status = -999;
     }


     // fill diphoton info
     ExoDiPhotons::FillDiphotonInfo(fDiphotonInfo,photon1,photon2);
     
     // worthwhile to also consider 'truth' diphoton info from genParticles?


     // only fill the tree if there are two photons!
     fTree->Fill();
     // and count this event for normalisation purposes
     fNorm_h->Fill(1);


   } //end of diphoton requirement
   else {

     // then there are not two tight photons in the event
     // hence count this in the norm hist
     fNorm_h->Fill(0);
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
ExoDiPhotonBkgAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ExoDiPhotonBkgAnalyzer::endJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
Float_t 
ExoDiPhotonBkgAnalyzer::getESRatio(const reco::Photon *photon, const edm::Event& e, const edm::EventSetup& iSetup){

  //get Geometry
  edm::ESHandle<CaloGeometry> caloGeometry;
  iSetup.get<CaloGeometryRecord>().get(caloGeometry);
  const CaloSubdetectorGeometry *geometry = caloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
  const CaloSubdetectorGeometry *& geometry_p = geometry;

  // Get ES rechits
  edm::Handle<EcalRecHitCollection> PreshowerRecHits;
  e.getByLabel(edm::InputTag("ecalPreshowerRecHit","EcalRecHitsES"), PreshowerRecHits);
  if( PreshowerRecHits.isValid() ) EcalRecHitCollection preshowerHits(*PreshowerRecHits);

  Float_t esratio=1.;

  if (fabs(photon->eta())>1.62) {

    const reco::CaloClusterPtr seed = (*photon).superCluster()->seed();    
    reco::CaloCluster cluster = (*seed);
    const GlobalPoint phopoint(cluster.x(), cluster.y(), cluster.z());
  
    DetId photmp1 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(phopoint, 1);
    DetId photmp2 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(phopoint, 2);
    ESDetId esfid = (photmp1 == DetId(0)) ? ESDetId(0) : ESDetId(photmp1);
    ESDetId esrid = (photmp2 == DetId(0)) ? ESDetId(0) : ESDetId(photmp2);

    int gs_esfid = -99;
    int gs_esrid = -99;
    gs_esfid = esfid.six()*32+esfid.strip();
    gs_esrid = esrid.siy()*32+esrid.strip();

    float esfe3 = 0.; 
    float esfe21 = 0.; 
    float esre3 = 0.; 
    float esre21 = 0.;

    const ESRecHitCollection *ESRH = PreshowerRecHits.product();
    EcalRecHitCollection::const_iterator esrh_it;
    for ( esrh_it = ESRH->begin(); esrh_it != ESRH->end(); esrh_it++) {
      ESDetId esdetid = ESDetId(esrh_it->id());
      if ( esdetid.plane()==1 ) {
	if ( esdetid.zside() == esfid.zside() &&
	     esdetid.siy() == esfid.siy() ) {
	  int gs_esid = esdetid.six()*32+esdetid.strip();
	  int ss = gs_esid-gs_esfid;
	  if ( TMath::Abs(ss)<=10) {
	    esfe21 += esrh_it->energy();
	  } 
	  if ( TMath::Abs(ss)<=1) {
	    esfe3 += esrh_it->energy();
	  } 
	}
      }
      if (esdetid.plane()==2 ){
	if ( esdetid.zside() == esrid.zside() &&
	     esdetid.six() == esrid.six() ) {
	  int gs_esid = esdetid.siy()*32+esdetid.strip();
	  int ss = gs_esid-gs_esrid;
	  if ( TMath::Abs(ss)<=10) {
	    esre21 += esrh_it->energy();
	  } 
	  if ( TMath::Abs(ss)<=1) {
	    esre3 += esrh_it->energy();
	  } 
	}
      }
    }
  
    if( (esfe21+esre21) == 0.) {
      esratio = 1.;
    }else{
      esratio = (esfe3+esre3) / (esfe21+esre21);
    }
    
  }
  return esratio;
  
}




//define this as a plug-in
DEFINE_FWK_MODULE(ExoDiPhotonBkgAnalyzer);
