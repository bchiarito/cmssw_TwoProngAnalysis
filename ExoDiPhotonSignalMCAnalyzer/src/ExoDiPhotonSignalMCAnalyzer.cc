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
// $Id$
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
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"


//for photons
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TMath.h"



//for trigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"
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


using namespace std;

//
// class declaration
//

//structs for my output tree info


struct eventInfo_t{
  int run;
  int LS;
  int evnum;
};



struct vtxInfo_t{
  int Nvtx; // number of reco vertices in event
  // but I will only keep the info of the best two, I think
  double vx;
  double vy;
  double vz;
  int isFake;
  int Ntracks;
  double sumPtTracks;
  // ************* chi2 or some other quality criteria *****************
  // these are the two vars used in GoodVertexFilter module
  double ndof;
  double d0; 

};

// for mc truth info, want to match to both status 3 and status 1 objects

// what if there is more than one status 1 object in a reasonable dR cone?

struct mcTrueObjectInfo_t {
  int status; // will be 3 or 1 
  int PdgId;
  int MotherPdgId;
  int GrandmotherPdgId; 
  double pt;
  double eta;
  double phi;
  //store deltaR of the match?
};

struct recoPhotonInfo_t {
  double pt;
  double eta;
  double phi;
  // position in ECAL  - caloPosition;
  double detEta;
  double detPhi; // clearly should be identical to phi, so am using as sort of cross-check 
  // which detector channel, specifically?
  //useful to cross check channel masking and look for problem channels, etc, I think


  //double check the vertex assigned to the photon candidate
  double vx;
  double vy;
  double vz;
  

  //shower shape variables
  double r9;
  double sigmaIetaIeta;
  double sigmaEtaEta;
  double maxEnergyXtal;
  // eNxN ...
  double e1x5;
  double e2x5;
  double e3x3;
  double e5x5;
  double r1x5;
  double r2x5;

  // swiss-cross style info
  double eLeft;
  double eRight;
  double eBottom;
  double eTop;
  double eSwissCross;
  double eTwo; // this is eMax + highest of {eLeft,eRight,eBottom,eTop}
  
  double hadOverEm; 
  // note also that two hadronic depths are available
  double hadDepth1OverEm; 
  double hadDepth2OverEm; 


  //isolation variables
  // these have different options: cone size, hollowness, etc
  // must include them all!
  float hcalIso04;
  float hcalIso03;
  float ecalIso04;
  float ecalIso03;
  float trkIsoSumPtHollow04;
  float trkIsoSumPtSolid04;
  int trkIsoNtrksHollow04;
  int trkIsoNtrksSolid04;
  float trkIsoSumPtHollow03;
  float trkIsoSumPtSolid03;
  int trkIsoNtrksHollow03;
  int trkIsoNtrksSolid03;

  // supercluster info
  double scRawEnergy;
  double scPreshowerEnergy;
  double scPhiWidth;
  double scEtaWidth;
  int scNumBasicClusters; // number of basic clusters comprising superCluster

  // seed cluster info
  
  
  //fiducial flags
  bool isEB;//Photon is in EB
  bool isEE;//Photon is in EE
  bool isEBEtaGap;//Photon is in supermodule/supercrystal eta gap in EB
  bool isEBPhiGap;//Photon is in supermodule/supercrystal phi gap in EB
  bool isEERingGap;//Photon is in crystal ring gap in EE
  bool isEEDeeGap;//Photon is in crystal dee gap in EE
  bool isEBEEGap;//Photon is in border between EB and EE.

  // pixel seed match?
  bool hasPixelSeed;
  // note to self: weird problems with this var in middle of struct - try at end

  bool isTight; // convenient to evaluate this here, and store result in tree
  // clearly the tightID needs to be hardcoded here, in that case
  
};



class ExoDiPhotonSignalMCAnalyzer : public edm::EDAnalyzer {
   public:
      explicit ExoDiPhotonSignalMCAnalyzer(const edm::ParameterSet&);
      ~ExoDiPhotonSignalMCAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

  // my functions
  bool isTightPhoton(const reco::Photon *photon); // handy way to add tightID as bool in tree


      // ----------member data ---------------------------
  
  // input tags and parameters
      edm::InputTag      fPhotonTag;       //select photon collection 
      double             fMin_pt;          // min pt cut (photons)

      // tools for clusters
      std::auto_ptr<EcalClusterLazyTools> lazyTools_;


      // my Tree
      TTree *fTree;

      eventInfo_t fEventInfo;
      vtxInfo_t fVtxInfo;
      mcTrueObjectInfo_t fSignalPhoton1Info; // leading signal photon
      mcTrueObjectInfo_t fSignalPhoton2Info;
      recoPhotonInfo_t fRecoPhotonInfo1; // leading matched reco photon 
      recoPhotonInfo_t fRecoPhotonInfo2; // second photon
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
    fMin_pt(iConfig.getUntrackedParameter<double>("ptMin"))
{
   //now do what ever initialization is needed

  edm::Service<TFileService> fs;
  fTree = fs->make<TTree>("fTree","PhotonTree");

  fTree->Branch("Event",&fEventInfo,"run/I:LS:evnum");
  fTree->Branch("Vtx",&fVtxInfo,"Nvtx/I:vx/D:vy:vz:isFake/I:Ntracks/I:sumPtTracks/D:ndof:d0");
  fTree->Branch("SignalPhoton1",&fSignalPhoton1Info,"status/I:PdgId:MotherPdgId:GrandmotherPdgId:pt/D:eta/D:phi/D");
  fTree->Branch("SignalPhoton2",&fSignalPhoton2Info,"status/I:PdgId:MotherPdgId:GrandmotherPdgId:pt/D:eta/D:phi/D");

  //try pixel seed at end -seems to work better with all booleans at end of branch!
  fTree->Branch("MatchRecoPhoton1",&fRecoPhotonInfo1,"pt/D:eta:phi:detEta:detPhi:vx:vy:vz:r9:sigmaIetaIeta:sigmaEtaEta:maxEnergyXtal:e1x5:e2x5:e3x3:e5x5:r1x5:r2x5:eLeft:eRight:eBottom:eTop:eSwissCross:eTwo:hadOverEm:hadDepth1OverEm:hadDepth2OverEm:hcalIso04/f:hcalIso03/f:ecalIso04:ecalIso03:trkIsoSumPtHollow04:trkIsoSumPtSolid04:trkIsoNtrksHollow04/I:trkIsoNtrksSolid04/I:trkIsoSumPtHollow03/f:trkIsoSumPtSolid03/f:trkIsoNtrksHollow03/I:trkIsoNtrksSolid03/I:scRawEnergy/D:scPreshowerEnergy:scPhiWidth:scEtaWidth:scNumBasicClusters/I:isEB/O:isEE:isEBEtaGap:isEBPhiGap:isEERingGap:isEEDeeGap:isEBEEGap:hasPixelSeed:isTight");

  fTree->Branch("MatchRecoPhoton2",&fRecoPhotonInfo2,"pt/D:eta:phi:detEta:detPhi:vx:vy:vz:r9:sigmaIetaIeta:sigmaEtaEta:maxEnergyXtal:e1x5:e2x5:e3x3:e5x5:r1x5:r2x5:eLeft:eRight:eBottom:eTop:eSwissCross:eTwo:hadOverEm:hadDepth1OverEm:hadDepth2OverEm:hcalIso04/f:hcalIso03/f:ecalIso04:ecalIso03:trkIsoSumPtHollow04:trkIsoSumPtSolid04:trkIsoNtrksHollow04/I:trkIsoNtrksSolid04/I:trkIsoSumPtHollow03/f:trkIsoSumPtSolid03/f:trkIsoNtrksHollow03/I:trkIsoNtrksSolid03/I:scRawEnergy/D:scPreshowerEnergy:scPhiWidth:scEtaWidth:scNumBasicClusters/I:isEB/O:isEE:isEBEtaGap:isEBPhiGap:isEERingGap:isEEDeeGap:isEBEEGap:hasPixelSeed:isTight");

  // signal diphoton info? eg to probe true MC width?
  // reco diphoton info?

 

}


ExoDiPhotonSignalMCAnalyzer::~ExoDiPhotonSignalMCAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
bool ExoDiPhotonSignalMCAnalyzer::isTightPhoton(const reco::Photon *photon)
{
  bool result = false;

  // these cuts are just hardcoded for now ...

  bool hadOverEmResult = false;
  bool trkIsoResult =false;
  bool hcalIsoResult = false;
  bool ecalIsoResult = false;
  bool noPixelSeedResult = false;
  bool sigmaIetaIetaResult = false;
  
  if(photon->hadronicOverEm()<0.05)
    hadOverEmResult = true;


  double trkIsoCut = 2.0 + 0.001*photon->et();
  if(photon->trkSumPtHollowConeDR04()<trkIsoCut)
    trkIsoResult = true;

  double hcalIsoCut = 2.2 + 0.001*photon->et();
  if(photon->hcalTowerSumEtConeDR04()<hcalIsoCut) 
    hcalIsoResult = true;

  double ecalIsoCut = 4.2 + 0.003*photon->et();
  if(photon->ecalRecHitSumEtConeDR04()<ecalIsoCut)
    ecalIsoResult = true;

  if(photon->hasPixelSeed()==false)
    noPixelSeedResult = true; // ie it is true that it does NOT have a pixel seed!

  // sigmaIetaIeta should be included in tight photon ID too soon ...

  if(hadOverEmResult && trkIsoResult && ecalIsoResult && hcalIsoResult && noPixelSeedResult) 
    result = true;

  return result; 
}




// ------------ method called to for each event  ------------
void
ExoDiPhotonSignalMCAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // basic event info
   fEventInfo.run = iEvent.id().run();
   fEventInfo.LS = iEvent.id().luminosityBlock();
   fEventInfo.evnum = iEvent.id().event();

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
     fVtxInfo.vx = vertex1->x();
     fVtxInfo.vy = vertex1->y();
     fVtxInfo.vz = vertex1->z();
     fVtxInfo.isFake = vertex1->isFake(); 
     fVtxInfo.Ntracks = vertex1->tracksSize();
     fVtxInfo.sumPtTracks = highestSumPtTracks1;
     fVtxInfo.ndof = vertex1->ndof();
     fVtxInfo.d0 = vertex1->position().rho();
     
   }


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
     fSignalPhoton1Info.status = signalPhoton1->status();
     fSignalPhoton1Info.PdgId = signalPhoton1->pdgId();
     fSignalPhoton1Info.MotherPdgId = signalPhoton1->mother()->pdgId();
     fSignalPhoton1Info.GrandmotherPdgId = signalPhoton1->mother()->mother()->pdgId();
     fSignalPhoton1Info.pt = signalPhoton1->pt();
     fSignalPhoton1Info.eta = signalPhoton1->eta();
     fSignalPhoton1Info.phi = signalPhoton1->phi();
   }

   if(signalPhoton2) {
     fSignalPhoton2Info.status = signalPhoton2->status();
     fSignalPhoton2Info.PdgId = signalPhoton2->pdgId();
     fSignalPhoton2Info.MotherPdgId = signalPhoton2->mother()->pdgId();
     fSignalPhoton2Info.GrandmotherPdgId = signalPhoton2->mother()->mother()->pdgId();
     fSignalPhoton2Info.pt = signalPhoton2->pt();
     fSignalPhoton2Info.eta = signalPhoton2->eta();
     fSignalPhoton2Info.phi = signalPhoton2->phi();		 
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
     fRecoPhotonInfo1.pt = matchPhoton1->et();
     fRecoPhotonInfo1.eta = matchPhoton1->eta();
     fRecoPhotonInfo1.phi = matchPhoton1->phi();

     fRecoPhotonInfo1.detEta = matchPhoton1->caloPosition().eta();
     fRecoPhotonInfo1.detPhi = matchPhoton1->caloPosition().phi();     

     fRecoPhotonInfo1.vx = matchPhoton1->vx();
     fRecoPhotonInfo1.vy = matchPhoton1->vy();
     fRecoPhotonInfo1.vz = matchPhoton1->vz();


     fRecoPhotonInfo1.r9 = matchPhoton1->r9();
     fRecoPhotonInfo1.sigmaIetaIeta = matchPhoton1->sigmaIetaIeta();
     fRecoPhotonInfo1.sigmaEtaEta = matchPhoton1->sigmaEtaEta();
     fRecoPhotonInfo1.maxEnergyXtal = matchPhoton1->maxEnergyXtal();

     fRecoPhotonInfo1.e1x5 = matchPhoton1->e1x5();
     fRecoPhotonInfo1.e2x5 = matchPhoton1->e2x5();
     fRecoPhotonInfo1.e3x3 = matchPhoton1->e3x3();
     fRecoPhotonInfo1.e5x5 = matchPhoton1->e5x5();
     fRecoPhotonInfo1.r1x5 = matchPhoton1->r1x5();
     fRecoPhotonInfo1.r2x5 = matchPhoton1->r2x5();

     // swiss cross related
     reco::SuperClusterRef sc = matchPhoton1->superCluster();
     std::pair<DetId,float> maxc = lazyTools_->getMaximum(*sc);
     double scross = -999.99;
     if (maxc.first.subdetId() == EcalBarrel) 
       scross = EcalSeverityLevelAlgo::swissCross(maxc.first, (*recHitsEB), 0);
     fRecoPhotonInfo1.eSwissCross = scross;
     fRecoPhotonInfo1.eLeft = lazyTools_->eLeft(*sc);
     fRecoPhotonInfo1.eRight = lazyTools_->eRight(*sc);
     fRecoPhotonInfo1.eBottom = lazyTools_->eBottom(*sc);
     fRecoPhotonInfo1.eTop = lazyTools_->eTop(*sc);
     
     // sum of the max + highest of the {left,right,bottom,top}
     std::vector<double> eAround;
     eAround.push_back(lazyTools_->eLeft(*sc));
     eAround.push_back(lazyTools_->eRight(*sc));
     eAround.push_back(lazyTools_->eBottom(*sc));
     eAround.push_back(lazyTools_->eTop(*sc));
     sort(eAround.begin(),eAround.end());
     reverse(eAround.begin(),eAround.end());
     fRecoPhotonInfo1.eTwo = lazyTools_->eMax(*sc)+eAround[0];

     cout << "Swiss cross = " << fRecoPhotonInfo1.eSwissCross;
     cout << "; eTwo/e5x5 = " << fRecoPhotonInfo1.eTwo/fRecoPhotonInfo1.e5x5;
     cout<< endl;

     fRecoPhotonInfo1.hadOverEm = matchPhoton1->hadronicOverEm();
     fRecoPhotonInfo1.hadDepth1OverEm = matchPhoton1->hadronicDepth1OverEm();
     fRecoPhotonInfo1.hadDepth2OverEm = matchPhoton1->hadronicDepth2OverEm();
     
     fRecoPhotonInfo1.ecalIso04 = matchPhoton1->ecalRecHitSumEtConeDR04();
     fRecoPhotonInfo1.ecalIso03 = matchPhoton1->ecalRecHitSumEtConeDR03();
     fRecoPhotonInfo1.hcalIso04 = matchPhoton1->hcalTowerSumEtConeDR04();
     fRecoPhotonInfo1.hcalIso03 = matchPhoton1->hcalTowerSumEtConeDR03();

     fRecoPhotonInfo1.trkIsoSumPtHollow04 = matchPhoton1->trkSumPtHollowConeDR04();
     fRecoPhotonInfo1.trkIsoSumPtSolid04 = matchPhoton1->trkSumPtSolidConeDR04();
     fRecoPhotonInfo1.trkIsoNtrksHollow04 = matchPhoton1->nTrkHollowConeDR04();
     fRecoPhotonInfo1.trkIsoNtrksSolid04 = matchPhoton1->nTrkSolidConeDR04();
     
     fRecoPhotonInfo1.trkIsoSumPtHollow03 = matchPhoton1->trkSumPtHollowConeDR03();
     fRecoPhotonInfo1.trkIsoSumPtSolid03 = matchPhoton1->trkSumPtSolidConeDR03();
     fRecoPhotonInfo1.trkIsoNtrksHollow03 = matchPhoton1->nTrkHollowConeDR03();
     fRecoPhotonInfo1.trkIsoNtrksSolid03 = matchPhoton1->nTrkSolidConeDR03();

     fRecoPhotonInfo1.hasPixelSeed = matchPhoton1->hasPixelSeed();

     fRecoPhotonInfo1.isEB        = matchPhoton1->isEB();        
     fRecoPhotonInfo1.isEE	 = matchPhoton1->isEE();	 
     fRecoPhotonInfo1.isEBEtaGap	 = matchPhoton1->isEBEtaGap();	 
     fRecoPhotonInfo1.isEBPhiGap	 = matchPhoton1->isEBPhiGap();	 
     fRecoPhotonInfo1.isEERingGap = matchPhoton1->isEERingGap(); 
     fRecoPhotonInfo1.isEEDeeGap	 = matchPhoton1->isEEDeeGap();	 
     fRecoPhotonInfo1.isEBEEGap   = matchPhoton1->isEBEEGap();
     
     fRecoPhotonInfo1.scRawEnergy = matchPhoton1->superCluster()->rawEnergy();
     fRecoPhotonInfo1.scPreshowerEnergy = matchPhoton1->superCluster()->preshowerEnergy();
     fRecoPhotonInfo1.scPhiWidth = matchPhoton1->superCluster()->phiWidth();
     fRecoPhotonInfo1.scEtaWidth = matchPhoton1->superCluster()->etaWidth();
     fRecoPhotonInfo1.scNumBasicClusters = matchPhoton1->superCluster()->clustersSize();


     // also to add seed cluster info
     
     // check tight ID
     fRecoPhotonInfo1.isTight = isTightPhoton(matchPhoton1);

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
     fRecoPhotonInfo2.pt = matchPhoton2->et();
     fRecoPhotonInfo2.eta = matchPhoton2->eta();
     fRecoPhotonInfo2.phi = matchPhoton2->phi();

     fRecoPhotonInfo2.detEta = matchPhoton2->caloPosition().eta();
     fRecoPhotonInfo2.detPhi = matchPhoton2->caloPosition().phi();     

     fRecoPhotonInfo2.vx = matchPhoton2->vx();
     fRecoPhotonInfo2.vy = matchPhoton2->vy();
     fRecoPhotonInfo2.vz = matchPhoton2->vz();


     fRecoPhotonInfo2.r9 = matchPhoton2->r9();
     fRecoPhotonInfo2.sigmaIetaIeta = matchPhoton2->sigmaIetaIeta();
     fRecoPhotonInfo2.sigmaEtaEta = matchPhoton2->sigmaEtaEta();
     fRecoPhotonInfo2.maxEnergyXtal = matchPhoton2->maxEnergyXtal();

     fRecoPhotonInfo2.e1x5 = matchPhoton2->e1x5();
     fRecoPhotonInfo2.e2x5 = matchPhoton2->e2x5();
     fRecoPhotonInfo2.e3x3 = matchPhoton2->e3x3();
     fRecoPhotonInfo2.e5x5 = matchPhoton2->e5x5();
     fRecoPhotonInfo2.r1x5 = matchPhoton2->r1x5();
     fRecoPhotonInfo2.r2x5 = matchPhoton2->r2x5();

     // swiss cross related
     reco::SuperClusterRef sc = matchPhoton2->superCluster();
     std::pair<DetId,float> maxc = lazyTools_->getMaximum(*sc);
     double scross = -999.99;
     if (maxc.first.subdetId() == EcalBarrel) 
       scross = EcalSeverityLevelAlgo::swissCross(maxc.first, (*recHitsEB), 0);
     fRecoPhotonInfo2.eSwissCross = scross;
     fRecoPhotonInfo2.eLeft = lazyTools_->eLeft(*sc);
     fRecoPhotonInfo2.eRight = lazyTools_->eRight(*sc);
     fRecoPhotonInfo2.eBottom = lazyTools_->eBottom(*sc);
     fRecoPhotonInfo2.eTop = lazyTools_->eTop(*sc);
     
     // sum of the max + highest of the {left,right,bottom,top}
     std::vector<double> eAround;
     eAround.push_back(lazyTools_->eLeft(*sc));
     eAround.push_back(lazyTools_->eRight(*sc));
     eAround.push_back(lazyTools_->eBottom(*sc));
     eAround.push_back(lazyTools_->eTop(*sc));
     sort(eAround.begin(),eAround.end());
     reverse(eAround.begin(),eAround.end());
     fRecoPhotonInfo2.eTwo = lazyTools_->eMax(*sc)+eAround[0];

     cout << "Swiss cross = " << fRecoPhotonInfo2.eSwissCross;
     cout << "; eTwo/e5x5 = " << fRecoPhotonInfo2.eTwo/fRecoPhotonInfo2.e5x5;
     cout<< endl;

     fRecoPhotonInfo2.hadOverEm = matchPhoton2->hadronicOverEm();
     fRecoPhotonInfo2.hadDepth1OverEm = matchPhoton2->hadronicDepth1OverEm();
     fRecoPhotonInfo2.hadDepth2OverEm = matchPhoton2->hadronicDepth2OverEm();
     
     fRecoPhotonInfo2.ecalIso04 = matchPhoton2->ecalRecHitSumEtConeDR04();
     fRecoPhotonInfo2.ecalIso03 = matchPhoton2->ecalRecHitSumEtConeDR03();
     fRecoPhotonInfo2.hcalIso04 = matchPhoton2->hcalTowerSumEtConeDR04();
     fRecoPhotonInfo2.hcalIso03 = matchPhoton2->hcalTowerSumEtConeDR03();

     fRecoPhotonInfo2.trkIsoSumPtHollow04 = matchPhoton2->trkSumPtHollowConeDR04();
     fRecoPhotonInfo2.trkIsoSumPtSolid04 = matchPhoton2->trkSumPtSolidConeDR04();
     fRecoPhotonInfo2.trkIsoNtrksHollow04 = matchPhoton2->nTrkHollowConeDR04();
     fRecoPhotonInfo2.trkIsoNtrksSolid04 = matchPhoton2->nTrkSolidConeDR04();
     
     fRecoPhotonInfo2.trkIsoSumPtHollow03 = matchPhoton2->trkSumPtHollowConeDR03();
     fRecoPhotonInfo2.trkIsoSumPtSolid03 = matchPhoton2->trkSumPtSolidConeDR03();
     fRecoPhotonInfo2.trkIsoNtrksHollow03 = matchPhoton2->nTrkHollowConeDR03();
     fRecoPhotonInfo2.trkIsoNtrksSolid03 = matchPhoton2->nTrkSolidConeDR03();

     fRecoPhotonInfo2.hasPixelSeed = matchPhoton2->hasPixelSeed();

     fRecoPhotonInfo2.isEB        = matchPhoton2->isEB();        
     fRecoPhotonInfo2.isEE	 = matchPhoton2->isEE();	 
     fRecoPhotonInfo2.isEBEtaGap	 = matchPhoton2->isEBEtaGap();	 
     fRecoPhotonInfo2.isEBPhiGap	 = matchPhoton2->isEBPhiGap();	 
     fRecoPhotonInfo2.isEERingGap = matchPhoton2->isEERingGap(); 
     fRecoPhotonInfo2.isEEDeeGap	 = matchPhoton2->isEEDeeGap();	 
     fRecoPhotonInfo2.isEBEEGap   = matchPhoton2->isEBEEGap();
     
     fRecoPhotonInfo2.scRawEnergy = matchPhoton2->superCluster()->rawEnergy();
     fRecoPhotonInfo2.scPreshowerEnergy = matchPhoton2->superCluster()->preshowerEnergy();
     fRecoPhotonInfo2.scPhiWidth = matchPhoton2->superCluster()->phiWidth();
     fRecoPhotonInfo2.scEtaWidth = matchPhoton2->superCluster()->etaWidth();
     fRecoPhotonInfo2.scNumBasicClusters = matchPhoton2->superCluster()->clustersSize();

     // check tight ID
     fRecoPhotonInfo2.isTight = isTightPhoton(matchPhoton2);
   }
   else {
     //     cout << "No match to signal photon2!" <<endl;
     // as short cut for indicating this in tree
     // make sure the pt value is crazy
     fRecoPhotonInfo2.pt = -9999.99;
   }


   //   cout << endl;
   //   cout << endl;



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
