// -*- C++ -*-
//
// Package:    PhotonAnalyzer
// Class:      PhotonAnalyzer
// 
/**\class PhotonAnalyzer PhotonAnalyzer.cc PhotonAnalysis/PhotonAnalyzer/src/PhotonAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Conor Henderson
//         Created:  Mon Oct  5 13:21:34 CEST 2009
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

// to use TfileService for histograms
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"

//for photons
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

//#include "DataFormats/Candidate/interface/LeafCandidate.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "TMath.h"

using namespace std;

//
// class decleration
//

class PhotonAnalyzer : public edm::EDAnalyzer {
   public:
      explicit PhotonAnalyzer(const edm::ParameterSet&);
      ~PhotonAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
      edm::InputTag      fPhotonTag;       //select photon collection 
      double             fMin_pt;          // min pt cut (photons)

  // my histograms, using the TFileService
  TH1D *fPhoton1Et_h; // et of highest pt photon in event
  TH1D *fPhoton2Et_h; // et of 2nd highest pt photon in event
  TH1D *fPhoton1Eta_h; // eta of highest pt photon in event
  TH1D *fPhoton2Eta_h; // eta of 2nd highest pt photon in event
  TH1D *fPhoton1Phi_h; // phi of highest pt photon in event
  TH1D *fPhoton2Phi_h; // phi of 2nd highest pt photon in event

  // invariant mass of diphotons
  TH1D *fInvMass_h; 
  // delta phi/R or qT of photons?
  TH1D *fDeltaPhi_h; 


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
PhotonAnalyzer::PhotonAnalyzer(const edm::ParameterSet& iConfig)
  : fPhotonTag(iConfig.getUntrackedParameter<edm::InputTag>("photonCollection")),
    fMin_pt(iConfig.getUntrackedParameter<double>("ptMin"))
{
   //now do what ever initialization is needed

  // create histograms
  edm::Service<TFileService> fs;
  fPhoton1Et_h = fs->make<TH1D>("fPhoton1Et_h","Photon1 pt",100,0,200);
  fPhoton2Et_h = fs->make<TH1D>("fPhoton2Et_h","Photon2 pt",100,0,200);
  fPhoton1Eta_h = fs->make<TH1D>("fPhoton1Eta_h","Photon1 eta",80,-4,4);
  fPhoton2Eta_h = fs->make<TH1D>("fPhoton2Eta_h","Photon2 eta",80,-4,4);
  fPhoton1Phi_h = fs->make<TH1D>("fPhoton1Phi_h","Photon1 phi",100,-2*TMath::Pi(),2*TMath::Pi());
  fPhoton2Phi_h = fs->make<TH1D>("fPhoton2Phi_h","Photon2 phi",100,-2*TMath::Pi(),2*TMath::Pi());
  fInvMass_h = fs->make<TH1D>("fInvMass_h","Diphoton Invariant Mass",100,0,200);
  fDeltaPhi_h = fs->make<TH1D>("fDeltaPhi_h","Diphoton delta phi",100,-2*TMath::Pi(),2*TMath::Pi());
}


PhotonAnalyzer::~PhotonAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
PhotonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // how is vertex handled for photons?

   // get the photon collection
   Handle<reco::PhotonCollection> photonColl;
   iEvent.getByLabel(fPhotonTag,photonColl);

   // If photon collection is empty, exit
   if (!photonColl.isValid()) 
     return;

   // this is a diphoton analyzer, so:
   if(photonColl->size()<2) {
     return;
   }

   // NOTE: change couts to the Message Logger...
   cout << "N photons = " << photonColl->size() <<endl;

   // we want the two highest Et photons for this analysis
   const reco::Photon *photon1 = NULL;
   const reco::Photon *photon2 = NULL;
   double highestEt1 = fMin_pt;
   double highestEt2 = fMin_pt;

   // photon loop
   for(reco::PhotonCollection::const_iterator recoPhoton = photonColl->begin(); recoPhoton!=photonColl->end(); recoPhoton++) {

     cout << "Photon et, eta, phi = " << recoPhoton->et() <<", "<<recoPhoton->eta()<< ", "<< recoPhoton->phi()<<endl;
     
     //***********************************************************
     // require tight and loose photons!
     //***********************************************************


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
   }

   if(photon1) {
     fPhoton1Et_h->Fill(photon1->et());
     fPhoton1Eta_h->Fill(photon1->eta());
     fPhoton1Phi_h->Fill(photon1->phi());
     cout << "Hghest Et photon: et, eta, phi = " << photon1->et() <<", "<<photon1->eta()<< ", "<< photon1->phi()<<endl;
   }
   if(photon2) {
     fPhoton2Et_h->Fill(photon2->et());
     fPhoton2Eta_h->Fill(photon2->eta());
     fPhoton2Phi_h->Fill(photon2->phi());
     cout << "2nd Highest Et photon: et, eta, phi = " << photon2->et() <<", "<<photon2->eta()<< ", "<< photon2->phi()<<endl;
   }
   if(photon1&&photon2) {
     //     photon1
     reco::LeafCandidate::LorentzVector photon_vector1 = photon1->p4();
     reco::LeafCandidate::LorentzVector photon_vector2 = photon2->p4();
     double invMass = (photon_vector1+photon_vector2).M();
     cout << "Inv Mass = " << invMass<<endl;
     fInvMass_h->Fill(invMass);
     // there's a CMS function for deltaPhi in DataFormats/Math
     double deltaPhiPhotons = reco::deltaPhi(photon1->phi(),photon2->phi());
     cout << "delta Phi = " << deltaPhiPhotons;
     fDeltaPhi_h->Fill(deltaPhiPhotons);

   }


   cout << endl;
   
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
PhotonAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhotonAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhotonAnalyzer);
