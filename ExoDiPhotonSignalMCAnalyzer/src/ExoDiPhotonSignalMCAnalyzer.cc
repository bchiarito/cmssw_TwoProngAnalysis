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
// $Id: ExoDiPhotonSignalMCAnalyzer.cc,v 1.3 2010/09/29 17:08:12 chenders Exp $
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
#include "L1Trigger/GlobalTrigger/interface/L1GlobalTrigger.h"
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

      // tools for clusters
      std::auto_ptr<EcalClusterLazyTools> lazyTools_;


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


  fTree->Branch("SignalPhoton1",&fSignalPhoton1Info,ExoDiPhotons::mcTrueObjectInfoBranchDefString.c_str());
  fTree->Branch("SignalPhoton2",&fSignalPhoton2Info,ExoDiPhotons::mcTrueObjectInfoBranchDefString.c_str());


  fTree->Branch("MatchRecoPhoton1",&fRecoPhotonInfo1,ExoDiPhotons::recoPhotonBranchDefString.c_str());

  fTree->Branch("MatchRecoPhoton2",&fRecoPhotonInfo2,ExoDiPhotons::recoPhotonBranchDefString.c_str());

  // signal diphoton info? eg to probe true MC width?
  // reco diphoton info?

   fTree->Branch("DiphotonSignal",&fDiphotonSignalInfo,ExoDiPhotons::diphotonInfoBranchDefString.c_str());
   fTree->Branch("DiphotonReco",&fDiphotonRecoInfo,ExoDiPhotons::diphotonInfoBranchDefString.c_str());



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
   const EcalChannelStatus *ch_status = chStatus.product(); 



   // get the photon collection
   Handle<reco::PhotonCollection> photonColl;
   iEvent.getByLabel(fPhotonTag,photonColl);

   // If photon collection is empty, exit
   if (!photonColl.isValid()) {
     cout << "No Photons! Move along, there's nothing to see here .." <<endl;
     return;
   }
   
//      cout << "N photons = " << photonColl->size() <<endl;

//      //   photon loop
//    for(reco::PhotonCollection::const_iterator recoPhoton = photonColl->begin(); recoPhoton!=photonColl->end(); recoPhoton++) {

//      cout << "Reco photon et, eta, phi = " << recoPhoton->et() <<", "<<recoPhoton->eta()<< ", "<< recoPhoton->phi();
//      cout << "; eMax/e3x3 = " << recoPhoton->maxEnergyXtal()/recoPhoton->e3x3();
//      cout << "; hadOverEm = " << recoPhoton->hadronicOverEm();
//      cout << "; trkIso = " << recoPhoton->trkSumPtHollowConeDR04();
//      cout << "; ecalIso = " << recoPhoton->ecalRecHitSumEtConeDR04();
//      cout << "; hcalIso = " << recoPhoton->hcalTowerSumEtConeDR04();
//      cout << "; pixelSeed = " << recoPhoton->hasPixelSeed();
//      cout << endl;


//    } // end reco photon loop

   



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


   for(reco::GenParticleCollection::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); ++genParticle) {

     
     // identify the status 1 (ie no further decay) photons
     // which came from the hard-scatt photons (status 3)
     // because I know they're there     

     if(genParticle->status()==1&&genParticle->pdgId()==22) {
       if(genParticle->numberOfMothers()>0) {
	 if(genParticle->mother()->status()==3&&genParticle->mother()->pdgId()==22) {
	   // further require that this status 3 photon itself came from the RS graviton 
	   if(genParticle->mother()->numberOfMothers()>0) {
	     if(genParticle->mother()->mother()->pdgId()==5000039) {

	       
	       //	       cout << "MC particle: Status = "<< genParticle->status() << "; pdg id = "<< genParticle->pdgId() << "; pt, eta, phi = " << genParticle->pt() << ", "<< genParticle->eta() << ", " << genParticle->phi() << endl;	   	       
	       

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
     
     //     cout << "MC particle: Status = "<< genParticle->status() << "; pdg id = "<< genParticle->pdgId() << "; pt, eta, phi = " << genParticle->pt() << ", "<< genParticle->eta() << ", " << genParticle->phi() << endl;	   

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
   
      // if no match found, then the match photon pointers are certainly NULL
   if(matchPhoton1) {		 
//      cout << "Matched to signal photon!" <<endl;
//      cout << "Reco photon et, eta, phi = " << matchPhoton1->et() <<", "<<matchPhoton1->eta()<< ", "<< matchPhoton1->phi();
//      cout << "; eMax/e3x3 = " << matchPhoton1->maxEnergyXtal()/matchPhoton1->e3x3();
//      cout << "; hadOverEm = " << matchPhoton1->hadronicOverEm();
//      cout << "; trkIso = " << matchPhoton1->trkSumPtHollowConeDR04();
//      cout << "; ecalIso = " << matchPhoton1->ecalRecHitSumEtConeDR04();
//      cout << "; hcalIso = " << matchPhoton1->hcalTowerSumEtConeDR04();
//      cout << "; pixelSeed = " << matchPhoton1->hasPixelSeed();
//      cout << endl;

     // fill info into tree
     
     ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo1,matchPhoton1,lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);


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
//      cout << "; hadOverEm = " << matchPhoton2->hadronicOverEm();
//      cout << "; trkIso = " << matchPhoton2->trkSumPtHollowConeDR04();
//      cout << "; ecalIso = " << matchPhoton2->ecalRecHitSumEtConeDR04();
//      cout << "; hcalIso = " << matchPhoton2->hcalTowerSumEtConeDR04();
//      cout << "; pixelSeed = " << matchPhoton2->hasPixelSeed();
//      cout << endl;

     // fill info into tree
     ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo2,matchPhoton2,lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);

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
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ExoDiPhotonSignalMCAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ExoDiPhotonSignalMCAnalyzer);
