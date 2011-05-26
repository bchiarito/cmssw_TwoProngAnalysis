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
// $Id: ExoDiPhotonBkgAnalyzer.cc,v 1.10 2010/11/22 15:51:13 chenders Exp $
//
//


// system include files
#include <memory>
#include "TClonesArray.h"
#include "TVector3.h"
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Utilities/interface/InputTag.h"

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
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"




// new CommonClasses approach
// these objects all in the namespace 'ExoDiPhotons'
#include "DiPhotonAnalysis/CommonClasses/interface/RecoPhotonInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/MCTrueObjectInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/PhotonID.h"
#include "DiPhotonAnalysis/CommonClasses/interface/EventAndVertexInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/DiphotonInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/TriggerInfo.h"

//new for PU gen                                                                                                                             
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

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
  //      Float_t getESRatio(const reco::Photon *photon, const edm::Event& e, const edm::EventSetup& iSetup);

      // ----------member data ---------------------------
  // input tags and parameters
      edm::InputTag      fPhotonTag;       //select photon collection 
      double             fMin_pt;          // min pt cut (photons)
      edm::InputTag      fHltInputTag;     // hltResults

  edm::InputTag     fRhoTag;
  edm::InputTag     pileupCollectionTag;

      // tools for clusters
      std::auto_ptr<EcalClusterLazyTools> lazyTools_;

      // my Tree
      TTree *fTree;
  
      ExoDiPhotons::eventInfo_t fEventInfo;
      ExoDiPhotons::vtxInfo_t fVtxInfo;
      ExoDiPhotons::beamSpotInfo_t fBeamSpotInfo;


  ExoDiPhotons::vtxInfo_t fVtxGENInfo;
  double rho;
  int pu_n;

  Int_t gv_n;

  TClonesArray* gv_pos;
  TClonesArray* gv_p3;

  Float_t gv_sumPtHi[100];
  Float_t gv_sumPtLo[100];
  Short_t gv_nTkHi[100];
  Short_t gv_nTkLo[100];

      ExoDiPhotons::hltTrigInfo_t fHLTInfo;

      int fNTightPhotons; // number of candidate photons in event (ie tight)

      ExoDiPhotons::recoPhotonInfo_t fRecoPhotonInfo1; // leading photon 
      ExoDiPhotons::recoPhotonInfo_t fRecoPhotonInfo2; // second photon

      // MC truth event-level info
      ExoDiPhotons::mcEventInfo_t fMCEventInfo;

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
    fHltInputTag(iConfig.getUntrackedParameter<edm::InputTag>("hltResults")),
    fRhoTag(iConfig.getParameter<edm::InputTag>("rhoCorrection")),
    pileupCollectionTag(iConfig.getUntrackedParameter<edm::InputTag>("pileupCorrection"))  
{
   //now do what ever initialization is needed

  edm::Service<TFileService> fs;
  fTree = fs->make<TTree>("fTree","PhotonTree");
  
  // now with CommonClasses, use the string defined in the header

  fTree->Branch("Event",&fEventInfo,ExoDiPhotons::eventInfoBranchDefString.c_str());
  
  // MC truth event info
  fTree->Branch("GenEvent",&fMCEventInfo,ExoDiPhotons::mcEventInfoBranchDefString.c_str());

  fTree->Branch("Vtx",&fVtxInfo,ExoDiPhotons::vtxInfoBranchDefString.c_str());
  fTree->Branch("VtxGEN",&fVtxGENInfo,ExoDiPhotons::vtxInfoBranchDefString.c_str());

  fTree->Branch("rho",&rho,"rho/D");
  fTree->Branch("pu_n", &pu_n, "pu_n/I");

  fTree->Branch("BeamSpot",&fBeamSpotInfo,ExoDiPhotons::beamSpotInfoBranchDefString.c_str());
  fTree->Branch("TrigHLT",&fHLTInfo,ExoDiPhotons::hltTrigBranchDefString.c_str());

  fTree->Branch("nTightPhotons",&fNTightPhotons,"nTightPhotons/I");

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
   fVtxInfo.isFake = true;   
   fVtxInfo.Ntracks = -99;
   fVtxInfo.sumPtTracks = -99999.99;
   fVtxInfo.ndof = -99999.99;
   fVtxInfo.d0 = -99999.99;


   // note for higher lumi, may want to also store second vertex, for pileup studies
   // specially if we look at pileup MC

   // to allow scalability for many vertices, use a vector and sort later
   std::vector<reco::Vertex> myVertices;

   for(reco::VertexCollection::const_iterator vtx=vertexColl->begin(); vtx!=vertexColl->end(); vtx++) {

     // push back all MC vertices for now
     
     myVertices.push_back(*vtx);

   }// end vertex loop

   // now sort the vertices after
   // can be either by Ntracks or TrackSumPt, depending how I write the function
   sort(myVertices.begin(),myVertices.end(),ExoDiPhotons::sortVertices) ;

   fVtxInfo.Nvtx = myVertices.size();
   // then fill vtx info
   ExoDiPhotons::FillVertexInfo(fVtxInfo,&(*myVertices.begin()));


   //GEN info:
   edm::Handle<reco::GenParticleCollection> gpH;
   iEvent.getByLabel("genParticles", gpH);


   gv_n = 0;
   //      TClonesArray* gv_pos;                                                                                                           

   gv_pos->Clear();
   gv_p3->Clear();

   const float lowPtThrGenVtx = 0.1;
   const float highPtThrGenVtx = 0.5;
   if (gpH.isValid() ) {
     for(reco::GenParticleCollection::const_iterator it_gen =
	   gpH->begin(); it_gen!= gpH->end(); it_gen++){
       if( it_gen->status() != 3 || !(it_gen->vx()!=0. || it_gen->vy()!=0. || it_gen->vx()!=0.)  ) { continue; }

       // check for duplicate vertex                                                                                                             
       bool duplicate = false;
       for(Int_t itv = 0; itv < gv_n; itv++) {
	 TVector3 * checkVtx = (TVector3 *) gv_pos->At(itv);
	 if( (fabs(it_gen->vx()-checkVtx->X())<1e-5) &&  (fabs(it_gen->vy()-checkVtx->Y())<1e-5) && (fabs(it_gen->vz()-checkVtx->Z())<1e-5)) {
	   duplicate = true;
	   break;
	 }
       }

       if (duplicate) continue;

       new((*gv_pos)[gv_n]) TVector3();
       ((TVector3 *) gv_pos->At(gv_n))->SetXYZ(it_gen->vx(), it_gen->vy(), it_gen->vz());

       TVector3 * this_gv_pos = (TVector3 *) gv_pos->At(gv_n);
       TVector3 p3(0,0,0);

       gv_sumPtLo[gv_n] = 0;
       gv_nTkLo[gv_n] = 0;
       gv_sumPtHi[gv_n] = 0;
       gv_nTkHi[gv_n] = 0;

       for(reco::GenParticleCollection::const_iterator part = gpH->begin(); part!= gpH->end(); part++){
	 if( part->status() == 1 && part->charge() != 0 && fabs(part->eta())<2.5 &&
	     ( fabs(part->vx()-this_gv_pos->X())<1.e-5 && fabs(part->vy()-this_gv_pos->Y())<1.e-5 && fabs(part->vz()-this_gv_pos->Z())<1.e-5 ) )\
	   {

	     TVector3 m(part->px(),part->py(),part->pz());
	     p3 += m;
	     if( m.Pt() > lowPtThrGenVtx ) {
	       gv_sumPtLo[gv_n] += m.Pt();
	       gv_nTkLo[gv_n] += 1;
	       if( m.Pt() > highPtThrGenVtx ) {
		 gv_sumPtHi[gv_n] += m.Pt();
		 gv_nTkHi[gv_n] += 1;
	       }
	     }
	   }
       }
       new((*gv_p3)[gv_n]) TVector3();
       ((TVector3 *) gv_p3->At(gv_n))->SetXYZ(p3.X(),p3.Y(),p3.Z());

       gv_n++;
     }
   }



   fVtxGENInfo.Nvtx = gv_n;

   edm::Handle<std::vector<PileupSummaryInfo> > pileupHandle;
   iEvent.getByLabel(pileupCollectionTag, pileupHandle);

   if (pileupHandle.isValid()){
     PileupSummaryInfo pileup = (*pileupHandle.product())[0];

     pu_n = pileup.getPU_NumInteractions();
   }

   //add rho correction      

   edm::Handle<double> rhoHandle;
   iEvent.getByLabel(fRhoTag, rhoHandle);

   rho = *(rhoHandle.product());

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

   lazyTools_ = std::auto_ptr<EcalClusterLazyTools>( new EcalClusterLazyTools(iEvent,iSetup,edm::InputTag("reducedEcalRecHitsEB"),edm::InputTag("reducedEcalRecHitsEE")) );
   
   // get ecal barrel recHits for spike rejection
   edm::Handle<EcalRecHitCollection> recHitsEB_h;
   iEvent.getByLabel(edm::InputTag("reducedEcalRecHitsEB"), recHitsEB_h );
   const EcalRecHitCollection * recHitsEB = 0;
   if ( ! recHitsEB_h.isValid() ) {
     LogError("ExoDiPhotonAnalyzer") << " ECAL Barrel RecHit Collection not available !"; return;
   } else {
     recHitsEB = recHitsEB_h.product();
   }

   edm::Handle<EcalRecHitCollection> recHitsEE_h;
   iEvent.getByLabel(edm::InputTag("reducedEcalRecHitsEE"), recHitsEE_h );
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

     // now add selected photons to vector if:
     // tight ID
     // not in gap
     // not a spike
     // min pt cut  (for all photons)
     if(ExoDiPhotons::isTightPhoton(&(*recoPhoton)) && !ExoDiPhotons::isGapPhoton(&(*recoPhoton)) && (recoPhoton->pt()>=fMin_pt) ) {

     // fill all reco photons for bkg study
     //     if(true) {

       selectedPhotons.push_back(*recoPhoton);
       
     } // end if tight photon
     
   } //end of reco photon loop

   
   // now sort the vector of selected photons by pt
   // (compare function is found in RecoPhotonInfo.h)
   sort(selectedPhotons.begin(),selectedPhotons.end(),ExoDiPhotons::comparePhotonsByPt);

   // now we have the two highest recoPhotons which pass our cuts

   // now count many candidate photons we have in this event
   //   cout << "N candidate photons = " << selectedPhotons.size() <<endl;
   fNTightPhotons = selectedPhotons.size();

   // require that we have two photons passing our cuts
   if(selectedPhotons.size()>=2) {

     ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo1,&selectedPhotons[0],lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);

     ExoDiPhotons::FillRecoPhotonInfo(fRecoPhotonInfo2,&selectedPhotons[1],lazyTools_.get(),recHitsEB,recHitsEE,ch_status,iEvent, iSetup);
     
     fRecoPhotonInfo1.isFakeable = false;
     fRecoPhotonInfo2.isFakeable = false;

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

     // first start with MC event info
     Handle<GenEventInfoProduct> genEventInfo;
     iEvent.getByLabel("generator",genEventInfo);

     if(!genEventInfo.isValid()) {
       cout << "gen event info not valid!"<<endl;
       return;
     }
     
     ExoDiPhotons::FillMCEventInfo(fMCEventInfo,genEventInfo.product());
     //     cout << "Process id = " << genEventInfo->signalProcessID() <<endl;
     //     const std::vector<double> genBinningValues = genEventInfo->binningValues();
     //     cout << "Binning value = " << genBinningValues[0] <<endl;

     // now do gen particles

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
       double deltaPhi1 = reco::deltaPhi(genParticle->phi(),selectedPhotons[0].phi());
       double deltaEta1 = genParticle->eta()-selectedPhotons[0].eta();
       double deltaR1 = TMath::Sqrt(deltaPhi1*deltaPhi1+deltaEta1*deltaEta1);
     
       double deltaPhi2 = reco::deltaPhi(genParticle->phi(),selectedPhotons[1].phi());
       double deltaEta2 = genParticle->eta()-selectedPhotons[1].eta();
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
     ExoDiPhotons::FillDiphotonInfo(fDiphotonInfo,&selectedPhotons[0],&selectedPhotons[1]);
     
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
// Float_t 
// ExoDiPhotonBkgAnalyzer::getESRatio(const reco::Photon *photon, const edm::Event& e, const edm::EventSetup& iSetup){

//   //get Geometry
//   edm::ESHandle<CaloGeometry> caloGeometry;
//   iSetup.get<CaloGeometryRecord>().get(caloGeometry);
//   const CaloSubdetectorGeometry *geometry = caloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
//   const CaloSubdetectorGeometry *& geometry_p = geometry;

//   // Get ES rechits
//   edm::Handle<EcalRecHitCollection> PreshowerRecHits;
//   e.getByLabel(edm::InputTag("ecalPreshowerRecHit","EcalRecHitsES"), PreshowerRecHits);
//   if( PreshowerRecHits.isValid() ) EcalRecHitCollection preshowerHits(*PreshowerRecHits);

//   Float_t esratio=1.;

//   if (fabs(photon->eta())>1.62) {

//     const reco::CaloClusterPtr seed = (*photon).superCluster()->seed();    
//     reco::CaloCluster cluster = (*seed);
//     const GlobalPoint phopoint(cluster.x(), cluster.y(), cluster.z());
  
//     DetId photmp1 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(phopoint, 1);
//     DetId photmp2 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(phopoint, 2);
//     ESDetId esfid = (photmp1 == DetId(0)) ? ESDetId(0) : ESDetId(photmp1);
//     ESDetId esrid = (photmp2 == DetId(0)) ? ESDetId(0) : ESDetId(photmp2);

//     int gs_esfid = -99;
//     int gs_esrid = -99;
//     gs_esfid = esfid.six()*32+esfid.strip();
//     gs_esrid = esrid.siy()*32+esrid.strip();

//     float esfe3 = 0.; 
//     float esfe21 = 0.; 
//     float esre3 = 0.; 
//     float esre21 = 0.;

//     const ESRecHitCollection *ESRH = PreshowerRecHits.product();
//     EcalRecHitCollection::const_iterator esrh_it;
//     for ( esrh_it = ESRH->begin(); esrh_it != ESRH->end(); esrh_it++) {
//       ESDetId esdetid = ESDetId(esrh_it->id());
//       if ( esdetid.plane()==1 ) {
// 	if ( esdetid.zside() == esfid.zside() &&
// 	     esdetid.siy() == esfid.siy() ) {
// 	  int gs_esid = esdetid.six()*32+esdetid.strip();
// 	  int ss = gs_esid-gs_esfid;
// 	  if ( TMath::Abs(ss)<=10) {
// 	    esfe21 += esrh_it->energy();
// 	  } 
// 	  if ( TMath::Abs(ss)<=1) {
// 	    esfe3 += esrh_it->energy();
// 	  } 
// 	}
//       }
//       if (esdetid.plane()==2 ){
// 	if ( esdetid.zside() == esrid.zside() &&
// 	     esdetid.six() == esrid.six() ) {
// 	  int gs_esid = esdetid.siy()*32+esdetid.strip();
// 	  int ss = gs_esid-gs_esrid;
// 	  if ( TMath::Abs(ss)<=10) {
// 	    esre21 += esrh_it->energy();
// 	  } 
// 	  if ( TMath::Abs(ss)<=1) {
// 	    esre3 += esrh_it->energy();
// 	  } 
// 	}
//       }
//     }
  
//     if( (esfe21+esre21) == 0.) {
//       esratio = 1.;
//     }else{
//       esratio = (esfe3+esre3) / (esfe21+esre21);
//     }
    
//   }
//   return esratio;
  
// }




//define this as a plug-in
DEFINE_FWK_MODULE(ExoDiPhotonBkgAnalyzer);
