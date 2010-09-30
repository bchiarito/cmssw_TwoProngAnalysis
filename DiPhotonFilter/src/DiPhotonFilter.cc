// -*- C++ -*-
//
// Package:    DiPhotonFilter
// Class:      DiPhotonFilter
// 
/**\class DiPhotonFilter DiPhotonFilter.cc DiPhotonAnalysis/DiPhotonFilter/src/DiPhotonFilter.cc

 Description: Simple filter for diphoton events

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Conor Henderson,40 1-B01,+41227671674,
//         Created:  Mon May 10 17:55:50 CEST 2010
// $Id: DiPhotonFilter.cc,v 1.1 2010/05/17 13:25:52 chenders Exp $
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//for photons
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

using namespace std;

//
// class declaration
//

class DiPhotonFilter : public edm::EDFilter {
   public:
      explicit DiPhotonFilter(const edm::ParameterSet&);
      ~DiPhotonFilter();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
      edm::InputTag      fPhotonTag;       //select photon collection 
      double             fMin_pt_photon1;  // min pt cut on leading photon
      double             fMin_pt_photon2;  // min pt cut on second photon

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
DiPhotonFilter::DiPhotonFilter(const edm::ParameterSet& iConfig)
  : fPhotonTag(iConfig.getUntrackedParameter<edm::InputTag>("photonCollection")),
    fMin_pt_photon1(iConfig.getUntrackedParameter<double>("ptMin_photon1")),
    fMin_pt_photon2(iConfig.getUntrackedParameter<double>("ptMin_photon2"))
{
   //now do what ever initialization is needed

}


DiPhotonFilter::~DiPhotonFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
DiPhotonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // the final result
   bool passEvent = false;

   // filter on two photons



   // get the photon collection
   Handle<reco::PhotonCollection> photonColl;
   iEvent.getByLabel(fPhotonTag,photonColl);

   // If photon collection is empty, exit
   if (!photonColl.isValid()) {
     cout << "No Photons! Move along, there's nothing to see here .." <<endl;
     return false;
   }

   // a quick short cut - if there are less than two reco photons in the collection,
   // just move on quickly
   if(photonColl->size()<2) {
     //     cout << "Less than two reco photons in event" <<endl;
     return false;
   }


   // can I just 'sort' the photonColl?
   // probably, if I write a simple compare function based on photon->pt()
   // but then how do I remove spikes?
   // so probably its just as simple to do things this way
   // in case I want to incorporate spike removal into the filter at a later point


   // we want the two highest Et photons for this analysis
   const reco::Photon *photon1 = NULL;
   const reco::Photon *photon2 = NULL;
   double highestEt1 = fMin_pt_photon1;
   double highestEt2 = fMin_pt_photon2;
   // by initialising these to the min thresholds for the two photons
   // then the outcome will be that both pointers will be only non-NULL
   // in the case where both thresholds passed
   
   // note that when photon2 cut is lower than photon1 cut, 
   // and the event contains a photon above the photon2 cut, but none above the photon1 cut
   // then what will happen is that the photon1 pointer will remain NULL, 
   // but the photon2 pointer will point to the highest pt photon 
   // (which is above the photon 2 cut)
   // this is apparently absurd
   // but its okay, the event wont pass anyway without both photons, 
   // so there is nothing to worry about 

   // photon loop
   for(reco::PhotonCollection::const_iterator recoPhoton = photonColl->begin(); recoPhoton!=photonColl->end(); recoPhoton++) {
     
     //     cout << "Photon et, eta, phi = " << recoPhoton->et() <<", "<<recoPhoton->eta()<< ", "<< recoPhoton->phi()<<endl;
     

     if(recoPhoton->et() >= highestEt1) {
       // then this is the new highest

       // but first make the old hihghest become the new second highest
       photon2 = photon1;
       // note that because of the way I have set this up, 
       // it is possible that photon1 is still NULL at this point
       // so be careful:
       if(photon2) {
	 highestEt2 = photon2->et();
       }
       // and now with this done, we can make the new one the highest
       photon1 =  &(*recoPhoton);
       highestEt1 = photon1->et();
       

     }
     else if(recoPhoton->et() >= highestEt2) {
       // then this is not highest, but is new second-highest
       photon2 =  &(*recoPhoton);
       highestEt2 = photon2->et();

     }


   } //end photon loop
   

   // now we have the two highest photons
//    if(photon1) {
//      cout << "Highest Et photon: et, eta, phi = " << photon1->et() <<", "<<photon1->eta()<< ", "<< photon1->phi()<<endl;
//    }
   
//     if(photon2) {
//       cout << "2nd Highest Et photon: et, eta, phi = " << photon2->et() <<", "<<photon2->eta()<< ", "<< photon2->phi()<<endl;      
//     }


    if((photon1)&&(photon1->et()>fMin_pt_photon1)&&(photon2)&&(photon2->et()>fMin_pt_photon2)) {
      //note that by my choice of initialisation, photon1 can never be non-NULL unless it is greater than fMin_pt_photon1, so the check should be redundant here
      // but I do it anyway
      
      // Anyway, if this is true, then we have passed this event
      
      //      cout << "Highest Et photon: et, eta, phi = " << photon1->et() <<", "<<photon1->eta()<< ", "<< photon1->phi()<<endl;
      //      cout << "2nd Highest Et photon: et, eta, phi = " << photon2->et() <<", "<<photon2->eta()<< ", "<< photon2->phi()<<endl;      
      //      cout << "This event passed" <<endl;

      passEvent = true;
    }

    return passEvent;


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
DiPhotonFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiPhotonFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiPhotonFilter);
