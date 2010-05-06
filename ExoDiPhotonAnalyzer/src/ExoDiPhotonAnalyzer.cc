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
// $Id$
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

//for photons
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TMath.h"

using namespace std;


//
// class declaration
//

//structs for my output tree info

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

// beam spot info
struct beamSpotInfo_t{
  double x0;
  double y0;
  double z0;
  double sigmaZ;
  double x0error;
  double y0error;
  double z0error;
  double sigmaZ0error;
  
};



//what other event info is relevant?
// at very least, identify the event: run, LS, event number
struct eventInfo_t{
  int run;
  int LS;
  int evnum;
};


// trigger info - MinBias trigger?  Definitely photon triggers!
struct trigInfo_t{
  bool HLT_Photon10; // for example
  bool L1_EG2; // also L1 bits
};

// single photon info
// store for each of the two photons in event


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
  
};


// also diphoton info: Minv, q_T, delta phi, etc

struct diphotonInfo_t{
  double Minv;
  double qt;
  double deltaPhi;
  double deltaEta;
  double deltaR;
  double cosThetaStar;
 
};




class ExoDiPhotonAnalyzer : public edm::EDAnalyzer {
   public:
      explicit ExoDiPhotonAnalyzer(const edm::ParameterSet&);
      ~ExoDiPhotonAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

  // my functions
  bool isTightPhoton(const reco::Photon *photon);
  bool isLoosePhoton(const reco::Photon *photon);


      // ----------member data ---------------------------
      edm::InputTag      fPhotonTag;       //select photon collection 
      double             fMin_pt;          // min pt cut (photons)

      // my Tree
      TTree *fTree;
      vtxInfo_t fVtxInfo;
      beamSpotInfo_t fBeamSpotInfo;
      eventInfo_t fEventInfo;
      trigInfo_t fTrigInfo;
      recoPhotonInfo_t fRecoPhotonInfo1; // leading photon 
      recoPhotonInfo_t fRecoPhotonInfo2; // second photon
      diphotonInfo_t fDiphotonInfo;
  

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
    fMin_pt(iConfig.getUntrackedParameter<double>("ptMin"))
{
   //now do what ever initialization is needed

  edm::Service<TFileService> fs;
  fTree = fs->make<TTree>("fTree","PhotonTree");
  fTree->Branch("Event",&fEventInfo,"run/I:LS:evnum");
  fTree->Branch("Vtx",&fVtxInfo,"Nvtx/I:vx/D:vy:vz:isFake/I:Ntracks/I:sumPtTracks/D:ndof:d0");
  fTree->Branch("BeamSpot",&fBeamSpotInfo,"x0/D:y0:z0:sigmaZ:x0error:y0error:z0error:sigmaZ0error");
  fTree->Branch("Trg",&fTrigInfo,"HLT_Photon10/O:L1_EG2");
  

  //try pixel seed at end
  fTree->Branch("Photon1",&fRecoPhotonInfo1,"pt/D:eta:phi:detEta:detPhi:vx:vy:vz:r9:sigmaIetaIeta:sigmaEtaEta:maxEnergyXtal:e1x5:e2x5:e3x3:e5x5:r1x5:r2x5:hadOverEm:hadDepth1OverEm:hadDepth2OverEm:hcalIso04/f:hcalIso03/f:ecalIso04:ecalIso03:trkIsoSumPtHollow04:trkIsoSumPtSolid04:trkIsoNtrksHollow04/I:trkIsoNtrksSolid04/I:trkIsoSumPtHollow03/f:trkIsoSumPtSolid03/f:trkIsoNtrksHollow03/I:trkIsoNtrksSolid03/I:scRawEnergy/D:scPreshowerEnergy:scPhiWidth:scEtaWidth:scNumBasicClusters/I:isEB/O:isEE:isEBEtaGap:isEBPhiGap:isEERingGap:isEEDeeGap:isEBEEGap:hasPixelSeed");

  fTree->Branch("Photon2",&fRecoPhotonInfo2,"pt/D:eta:phi:detEta:detPhi:vx:vy:vz:r9:sigmaIetaIeta:sigmaEtaEta:maxEnergyXtal:e1x5:e2x5:e3x3:e5x5:r1x5:r2x5:hadOverEm:hadDepth1OverEm:hadDepth2OverEm:hcalIso04/f:hcalIso03/f:ecalIso04:ecalIso03:trkIsoSumPtHollow04:trkIsoSumPtSolid04:trkIsoNtrksHollow04/I:trkIsoNtrksSolid04/I:trkIsoSumPtHollow03/f:trkIsoSumPtSolid03/f:trkIsoNtrksHollow03/I:trkIsoNtrksSolid03/I:scRawEnergy/D:scPreshowerEnergy:scPhiWidth:scEtaWidth:scNumBasicClusters/I:isEB/O:isEE:isEBEtaGap:isEBPhiGap:isEERingGap:isEEDeeGap:isEBEEGap:hasPixelSeed");


  fTree->Branch("Diphoton",&fDiphotonInfo,"Minv/D:qt:deltaPhi:deltaEta:deltaR:cosThetaStar");

}


ExoDiPhotonAnalyzer::~ExoDiPhotonAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
bool ExoDiPhotonAnalyzer::isTightPhoton(const reco::Photon *photon)
{
  bool result = false;

  return result;
}


bool ExoDiPhotonAnalyzer::isLoosePhoton(const reco::Photon *photon)
{
  bool result = false;

  return result;
}


// ------------ method called to for each event  ------------
void
ExoDiPhotonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
     fBeamSpotInfo.x0 = beamSpot.x0();
     fBeamSpotInfo.y0 = beamSpot.y0();
     fBeamSpotInfo.z0 = beamSpot.z0();
     fBeamSpotInfo.sigmaZ = beamSpot.sigmaZ();
     fBeamSpotInfo.x0error = beamSpot.x0Error();
     fBeamSpotInfo.y0error = beamSpot.y0Error();
     fBeamSpotInfo.z0error = beamSpot.z0Error();
     fBeamSpotInfo.sigmaZ0error = beamSpot.sigmaZ0Error();
   }



   // get the trig info
   fTrigInfo.HLT_Photon10 = false;
   fTrigInfo.L1_EG2 = false;



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

     //     cout << "Photon et, eta, phi = " << recoPhoton->et() <<", "<<recoPhoton->eta()<< ", "<< recoPhoton->phi()<<endl;


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


   } //end reco photn loop

//    // now we have the two highest photons
    if(photon1) {
      cout << "Highest Et photon: et, eta, phi = " << photon1->et() <<", "<<photon1->eta()<< ", "<< photon1->phi()<<endl;

// for photon tree
     fRecoPhotonInfo1.pt = photon1->et();
     fRecoPhotonInfo1.eta = photon1->eta();
     fRecoPhotonInfo1.phi = photon1->phi();

     fRecoPhotonInfo1.detEta = photon1->caloPosition().eta();
     fRecoPhotonInfo1.detPhi = photon1->caloPosition().phi();
     
     // since photon1 inherits from LeafCandidate, we can get the vertex position
     // that is associated with the photon:
     fRecoPhotonInfo1.vx = photon1->vx();
     fRecoPhotonInfo1.vy = photon1->vy();
     fRecoPhotonInfo1.vz = photon1->vz();


     fRecoPhotonInfo1.r9 = photon1->r9();
     fRecoPhotonInfo1.sigmaIetaIeta = photon1->sigmaIetaIeta();
     fRecoPhotonInfo1.sigmaEtaEta = photon1->sigmaEtaEta();
     fRecoPhotonInfo1.maxEnergyXtal = photon1->maxEnergyXtal();

     fRecoPhotonInfo1.e1x5 = photon1->e1x5();
     fRecoPhotonInfo1.e2x5 = photon1->e2x5();
     fRecoPhotonInfo1.e3x3 = photon1->e3x3();
     fRecoPhotonInfo1.e5x5 = photon1->e5x5();
     fRecoPhotonInfo1.r1x5 = photon1->r1x5();
     fRecoPhotonInfo1.r2x5 = photon1->r2x5();

     fRecoPhotonInfo1.hadOverEm = photon1->hadronicOverEm();
     fRecoPhotonInfo1.hadDepth1OverEm = photon1->hadronicDepth1OverEm();
     fRecoPhotonInfo1.hadDepth2OverEm = photon1->hadronicDepth2OverEm();

     
     fRecoPhotonInfo1.ecalIso04 = photon1->ecalRecHitSumEtConeDR04();
     fRecoPhotonInfo1.ecalIso03 = photon1->ecalRecHitSumEtConeDR03();
     fRecoPhotonInfo1.hcalIso04 = photon1->hcalTowerSumEtConeDR04();
     fRecoPhotonInfo1.hcalIso03 = photon1->hcalTowerSumEtConeDR03();

     fRecoPhotonInfo1.trkIsoSumPtHollow04 = photon1->trkSumPtHollowConeDR04();
     fRecoPhotonInfo1.trkIsoSumPtSolid04 = photon1->trkSumPtSolidConeDR04();
     fRecoPhotonInfo1.trkIsoNtrksHollow04 = photon1->nTrkHollowConeDR04();
     fRecoPhotonInfo1.trkIsoNtrksSolid04 = photon1->nTrkSolidConeDR04();
     
     fRecoPhotonInfo1.trkIsoSumPtHollow03 = photon1->trkSumPtHollowConeDR03();
     fRecoPhotonInfo1.trkIsoSumPtSolid03 = photon1->trkSumPtSolidConeDR03();
     fRecoPhotonInfo1.trkIsoNtrksHollow03 = photon1->nTrkHollowConeDR03();
     fRecoPhotonInfo1.trkIsoNtrksSolid03 = photon1->nTrkSolidConeDR03();


//      cout << "Trk iso cone 04: "<< photon1->trkSumPtHollowConeDR04();
//      cout << ", "<<  photon1->trkSumPtSolidConeDR04();
//      cout << ", "<<  photon1->nTrkHollowConeDR04();
//      cout << ", "<<  photon1->nTrkSolidConeDR04()<<endl;

//      cout << "In tree: " << fRecoPhotonInfo1.trkIsoSumPtHollow04;
//      cout << ", " << fRecoPhotonInfo1.trkIsoSumPtSolid04;
//      cout << ", " << fRecoPhotonInfo1.trkIsoNtrksSolid04;
     
     
//      cout << "Trk iso cone 03: "<< photon1->trkSumPtHollowConeDR03();
//      cout << ", "<<  photon1->trkSumPtSolidConeDR03();
//      cout << ", "<<  photon1->nTrkHollowConeDR03();
//      cout << ", "<<  photon1->nTrkSolidConeDR03();

//      cout << "In tree: " << fRecoPhotonInfo1.trkIsoSumPtHollow03;
//      cout << ", " << fRecoPhotonInfo1.trkIsoSumPtSolid03;
//      cout << ", " << fRecoPhotonInfo1.trkIsoNtrksSolid03;



     fRecoPhotonInfo1.hasPixelSeed = photon1->hasPixelSeed();

     //     cout << "Pixel seed: " << photon1->hasPixelSeed() <<endl;

     fRecoPhotonInfo1.isEB        = photon1->isEB();        
     fRecoPhotonInfo1.isEE	 = photon1->isEE();	 
     fRecoPhotonInfo1.isEBEtaGap	 = photon1->isEBEtaGap();	 
     fRecoPhotonInfo1.isEBPhiGap	 = photon1->isEBPhiGap();	 
     fRecoPhotonInfo1.isEERingGap = photon1->isEERingGap(); 
     fRecoPhotonInfo1.isEEDeeGap	 = photon1->isEEDeeGap();	 
     fRecoPhotonInfo1.isEBEEGap   = photon1->isEBEEGap();
     
     fRecoPhotonInfo1.scRawEnergy = photon1->superCluster()->rawEnergy();
     fRecoPhotonInfo1.scPreshowerEnergy = photon1->superCluster()->preshowerEnergy();
     fRecoPhotonInfo1.scPhiWidth = photon1->superCluster()->phiWidth();
     fRecoPhotonInfo1.scEtaWidth = photon1->superCluster()->etaWidth();
     fRecoPhotonInfo1.scNumBasicClusters = photon1->superCluster()->clustersSize();

     //     cout << "Nclusters = " << photon1->superCluster()->clustersSize()<<endl;

     // also to add seed cluster info
     


    }
    if(photon2) {
      cout << "2nd Highest Et photon: et, eta, phi = " << photon2->et() <<", "<<photon2->eta()<< ", "<< photon2->phi()<<endl;


// for photon tree
     fRecoPhotonInfo2.pt = photon2->et();
     fRecoPhotonInfo2.eta = photon2->eta();
     fRecoPhotonInfo2.phi = photon2->phi();

     fRecoPhotonInfo2.detEta = photon2->caloPosition().eta();
     fRecoPhotonInfo2.detPhi = photon2->caloPosition().phi();

     // channel number or something?
     
     // since photon2 inherits from LeafCandidate, we can get the vertex position
     // that is associated with the photon:
     fRecoPhotonInfo2.vx = photon2->vx();
     fRecoPhotonInfo2.vy = photon2->vy();
     fRecoPhotonInfo2.vz = photon2->vz();


     fRecoPhotonInfo2.r9 = photon2->r9();
     fRecoPhotonInfo2.sigmaIetaIeta = photon2->sigmaIetaIeta();
     fRecoPhotonInfo2.sigmaEtaEta = photon2->sigmaEtaEta();
     fRecoPhotonInfo2.maxEnergyXtal = photon2->maxEnergyXtal();

     fRecoPhotonInfo2.e1x5 = photon2->e1x5();
     fRecoPhotonInfo2.e2x5 = photon2->e2x5();
     fRecoPhotonInfo2.e3x3 = photon2->e3x3();
     fRecoPhotonInfo2.e5x5 = photon2->e5x5();
     fRecoPhotonInfo2.r1x5 = photon2->r1x5();
     fRecoPhotonInfo2.r2x5 = photon2->r2x5();

     fRecoPhotonInfo2.hadOverEm = photon2->hadronicOverEm();
     fRecoPhotonInfo2.hadDepth1OverEm = photon2->hadronicDepth1OverEm();
     fRecoPhotonInfo2.hadDepth2OverEm = photon2->hadronicDepth2OverEm();

     
     fRecoPhotonInfo2.ecalIso04 = photon2->ecalRecHitSumEtConeDR04();
     fRecoPhotonInfo2.ecalIso03 = photon2->ecalRecHitSumEtConeDR03();
     fRecoPhotonInfo2.hcalIso04 = photon2->hcalTowerSumEtConeDR04();
     fRecoPhotonInfo2.hcalIso03 = photon2->hcalTowerSumEtConeDR03();

     fRecoPhotonInfo2.trkIsoSumPtHollow04 = photon2->trkSumPtHollowConeDR04();
     fRecoPhotonInfo2.trkIsoSumPtSolid04 = photon2->trkSumPtSolidConeDR04();
     fRecoPhotonInfo2.trkIsoNtrksHollow04 = photon2->nTrkHollowConeDR04();
     fRecoPhotonInfo2.trkIsoNtrksSolid04 = photon2->nTrkSolidConeDR04();
     
     fRecoPhotonInfo2.trkIsoSumPtHollow03 = photon2->trkSumPtHollowConeDR03();
     fRecoPhotonInfo2.trkIsoSumPtSolid03 = photon2->trkSumPtSolidConeDR03();
     fRecoPhotonInfo2.trkIsoNtrksHollow03 = photon2->nTrkHollowConeDR03();
     fRecoPhotonInfo2.trkIsoNtrksSolid03 = photon2->nTrkSolidConeDR03();



//      cout << "Trk iso cone 04: "<< photon2->trkSumPtHollowConeDR04();
//      cout << ", "<<  photon2->trkSumPtSolidConeDR04();
//      cout << ", "<<  photon2->nTrkHollowConeDR04();
//      cout << ", "<<  photon2->nTrkSolidConeDR04()<<endl;

//      cout << "In tree: " << fRecoPhotonInfo2.trkIsoSumPtHollow04;
//      cout << ", " << fRecoPhotonInfo2.trkIsoSumPtSolid04;
//      cout << ", " << fRecoPhotonInfo2.trkIsoNtrksSolid04;
     
     
//      cout << "Trk iso cone 03: "<< photon2->trkSumPtHollowConeDR03();
//      cout << ", "<<  photon2->trkSumPtSolidConeDR03();
//      cout << ", "<<  photon2->nTrkHollowConeDR03();
//      cout << ", "<<  photon2->nTrkSolidConeDR03();

//      cout << "In tree: " << fRecoPhotonInfo2.trkIsoSumPtHollow03;
//      cout << ", " << fRecoPhotonInfo2.trkIsoSumPtSolid03;
//      cout << ", " << fRecoPhotonInfo2.trkIsoNtrksSolid03;



     fRecoPhotonInfo2.hasPixelSeed = photon2->hasPixelSeed();


     fRecoPhotonInfo2.isEB        = photon2->isEB();        
     fRecoPhotonInfo2.isEE	 = photon2->isEE();	 
     fRecoPhotonInfo2.isEBEtaGap	 = photon2->isEBEtaGap();	 
     fRecoPhotonInfo2.isEBPhiGap	 = photon2->isEBPhiGap();	 
     fRecoPhotonInfo2.isEERingGap = photon2->isEERingGap(); 
     fRecoPhotonInfo2.isEEDeeGap	 = photon2->isEEDeeGap();	 
     fRecoPhotonInfo2.isEBEEGap   = photon2->isEBEEGap();
     
     fRecoPhotonInfo2.scRawEnergy = photon2->superCluster()->rawEnergy();
     fRecoPhotonInfo2.scPreshowerEnergy = photon2->superCluster()->preshowerEnergy();
     fRecoPhotonInfo2.scPhiWidth = photon2->superCluster()->phiWidth();
     fRecoPhotonInfo2.scEtaWidth = photon2->superCluster()->etaWidth();
     fRecoPhotonInfo2.scNumBasicClusters = photon2->superCluster()->clustersSize();

     // also to add seed cluster info?
     

    }




    if(photon1&&photon2) {
      reco::LeafCandidate::LorentzVector photon_vector1 = photon1->p4();
      reco::LeafCandidate::LorentzVector photon_vector2 = photon2->p4();
      double invMass = (photon_vector1+photon_vector2).M();
      cout << "Inv Mass = " << invMass<<endl;
      fDiphotonInfo.Minv = invMass;
      // other diphotn variables
      //       deltaPhi, etc
      // pt of the pair
      fDiphotonInfo.qt = (photon_vector1+photon_vector2).pt();
      // there's a CMS function for deltaPhi in DataFormats/Math
      fDiphotonInfo.deltaPhi = reco::deltaPhi(photon1->phi(),photon2->phi());
      fDiphotonInfo.deltaEta = photon1->eta()-photon2->eta(); // always highest pt - second
      fDiphotonInfo.deltaR = TMath::Sqrt(fDiphotonInfo.deltaPhi*fDiphotonInfo.deltaPhi+fDiphotonInfo.deltaEta*fDiphotonInfo.deltaEta);
      fDiphotonInfo.cosThetaStar = -999.99; // need to calculate this!

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
