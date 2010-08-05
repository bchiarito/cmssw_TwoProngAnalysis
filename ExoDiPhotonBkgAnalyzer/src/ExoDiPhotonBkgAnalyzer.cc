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

      // ----------member data ---------------------------
  // input tags and parameters
      edm::InputTag      fPhotonTag;       //select photon collection 
      double             fMin_pt;          // min pt cut (photons)


      // my Tree
      TTree *fTree;
  
      ExoDiPhotons::eventInfo_t fEventInfo;
      ExoDiPhotons::vtxInfo_t fVtxInfo;
      ExoDiPhotons::beamSpotInfo_t fBeamSpotInfo;


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
    fMin_pt(iConfig.getUntrackedParameter<double>("ptMin"))
{
   //now do what ever initialization is needed

  edm::Service<TFileService> fs;
  fTree = fs->make<TTree>("fTree","PhotonTree");
  
  // now with CommonClasses, use the string defined in the header

  fTree->Branch("Event",&fEventInfo,ExoDiPhotons::eventInfoBranchDefString.c_str());
  fTree->Branch("Vtx",&fVtxInfo,ExoDiPhotons::vtxInfoBranchDefString.c_str());
  fTree->Branch("BeamSpot",&fBeamSpotInfo,ExoDiPhotons::beamSpotInfoBranchDefString.c_str());

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



   // trig info
   // TBC



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

     // need to fill swiss cross info by hand still

     // but for all other variables, 
     // can now use the Fill function specific to this recoPhoton struct
     ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo1,photon1);


   }

   if(photon2) {
     ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo2,photon2);
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

//define this as a plug-in
DEFINE_FWK_MODULE(ExoDiPhotonBkgAnalyzer);
