// -*- C++ -*-
//
// Package:    ExoDiPhotonGENOnlyAnalyzer
// Class:      ExoDiPhotonGENOnlyAnalyzer
// 
/**\class ExoDiPhotonGENOnlyAnalyzer ExoDiPhotonGENOnlyAnalyzer.cc DiPhotonAnalysis/ExoDiPhotonGENOnlyAnalyzer/src/ExoDiPhotonGENOnlyAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Conor Henderson,40 1-B01,+41227671674,
//         Created:  Thu Aug 16 13:43:39 CEST 2012
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

#include "FWCore/Framework/interface/ESHandle.h"

// to use TfileService for histograms and trees
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TTree.h"

// for MC
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// new CommonClasses approach
// these objects are all in the namespace 'ExoDiPhotons'
//#include "DiPhotonAnalysis/CommonClasses/interface/RecoPhotonInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/MCTrueObjectInfo.h"
//#include "DiPhotonAnalysis/CommonClasses/interface/TriggerInfo.h"
//#include "DiPhotonAnalysis/CommonClasses/interface/PhotonID.h"
#include "DiPhotonAnalysis/CommonClasses/interface/EventAndVertexInfo.h"
//#include "DiPhotonAnalysis/CommonClasses/interface/DiphotonInfo.h"

using namespace std;

//
// class declaration
//

class ExoDiPhotonGENOnlyAnalyzer : public edm::EDAnalyzer {
   public:
      explicit ExoDiPhotonGENOnlyAnalyzer(const edm::ParameterSet&);
      ~ExoDiPhotonGENOnlyAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
      // my Tree
      TTree *fTree;

      ExoDiPhotons::eventInfo_t fEventInfo;


      ExoDiPhotons::mcTrueObjectInfo_t fSignalPhoton1Info; // leading signal photon
      ExoDiPhotons::mcTrueObjectInfo_t fSignalPhoton2Info;

  double fGravitonMass; // for graviton mass
  

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
ExoDiPhotonGENOnlyAnalyzer::ExoDiPhotonGENOnlyAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
edm::Service<TFileService> fs;
  fTree = fs->make<TTree>("fTree","PhotonTree");

  // now with CommonClasses, use the string defined in the header

  fTree->Branch("Event",&fEventInfo,ExoDiPhotons::eventInfoBranchDefString.c_str());

  // for graviton
  fTree->Branch("RSGraviton",&fGravitonMass,"Mass/D");

  fTree->Branch("GenPhoton1",&fSignalPhoton1Info,ExoDiPhotons::mcTrueObjectInfoBranchDefString.c_str());
  fTree->Branch("GenPhoton2",&fSignalPhoton2Info,ExoDiPhotons::mcTrueObjectInfoBranchDefString.c_str());


  // signal diphoton info? eg to probe true MC width?
}


ExoDiPhotonGENOnlyAnalyzer::~ExoDiPhotonGENOnlyAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ExoDiPhotonGENOnlyAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // basic event info
   ExoDiPhotons::FillEventInfo(fEventInfo,iEvent);


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

   for(reco::GenParticleCollection::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); ++genParticle) {

     // first lets pick out the RS graviton
        if(genParticle->pdgId()==5000039) {
	  //	  cout << "MC particle: Status = "<< genParticle->status() << "; pdg id = "<< genParticle->pdgId() << "; pt, eta, phi = " << genParticle->pt() << ", "<< genParticle->eta() << ", " << genParticle->phi() << endl;	   	       
	  //	  cout << "Mass = "<< genParticle->mass() <<endl;
	  
	  // note that there is both a status 2 and 3 graviton  here ...
	  
	  fGravitonMass = genParticle->mass();
	}


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


	       // now assign our signal photons to 1 or 2

	       if(!signalPhoton1) {
		 // then we havent found the first one yet, so this is it
		 signalPhoton1 = &(*genParticle);
	       }
	       else {
		 // we have already found one, so this is the second
		 signalPhoton2 = &(*genParticle);
	       }

	     }
	   }
	 }
       }
     } //end status 1 req for  photons from RS graviton

   }// ends gen particle loop

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
     signalPhoton1 = signalPhoton2;
     signalPhoton2 = tempSignalPhoton;
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


   // put in diphoton info?


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
ExoDiPhotonGENOnlyAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ExoDiPhotonGENOnlyAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
ExoDiPhotonGENOnlyAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ExoDiPhotonGENOnlyAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ExoDiPhotonGENOnlyAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ExoDiPhotonGENOnlyAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ExoDiPhotonGENOnlyAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ExoDiPhotonGENOnlyAnalyzer);
