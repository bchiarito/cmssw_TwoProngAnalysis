// -*- C++ -*-
//
// Package:    PileupExtractor
// Class:      PileupExtractor
// 
/**\class PileupExtractor PileupExtractor.cc DiPhotonAnalysis/PileupExtractor/src/PileupExtractor.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Otman Charaf,42 1-015,+41227662353,
//         Created:  Thu Jan 17 14:57:25 CET 2013
// $Id$
//
//


// system include files
#include <memory>
#include <algorithm>
#include <vector>
#include <utility> 

// to use TfileService for histograms and trees
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TTree.h"
#include "TString.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

using namespace std;


//
// class declaration
//

class PileupExtractor : public edm::EDAnalyzer {
public:
  explicit PileupExtractor(const edm::ParameterSet&);
  ~PileupExtractor();

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
  edm::InputTag      fpileupCollectionTag;         
  string             fPUMCFileName;
  string             fPUMCHistName;   
  TH1F* fpu_n_BeforeCuts; 
  int fpu_n;
  int fBC;
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
PileupExtractor::PileupExtractor(const edm::ParameterSet& iConfig)
  :fpileupCollectionTag(iConfig.getUntrackedParameter<edm::InputTag>("pileupCollection"))
{
  //now do what ever initialization is needed

  edm::Service<TFileService> fs;
  fpu_n_BeforeCuts = fs->make<TH1F>("fpu_n_BeforeCuts","PileUpBeforeCuts",300,0,300);
}


PileupExtractor::~PileupExtractor()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PileupExtractor::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::Handle<std::vector<PileupSummaryInfo> > pileupHandle;
  iEvent.getByLabel(fpileupCollectionTag, pileupHandle);
  std::vector<PileupSummaryInfo>::const_iterator PUI;
   
  if (!pileupHandle.isValid()) cout<<"pu handle not valid"<<endl;

  if (pileupHandle.isValid()){
    for (PUI = pileupHandle->begin();PUI != pileupHandle->end(); ++PUI){      
      fBC = PUI->getBunchCrossing() ;
      //Select only the in time bunch crossing with bunch crossing=0
      if(fBC!=0) continue;      
      PileupSummaryInfo oldpileup = (*pileupHandle.product())[0];
      fpu_n = PUI->getTrueNumInteractions();
      fpu_n_BeforeCuts->Fill(fpu_n);         
    }//end of loop    
  }//end of condition of validity of pileup handle


}//end of analyze method


// ------------ method called once each job just before starting event loop  ------------
void 
PileupExtractor::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PileupExtractor::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
PileupExtractor::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
PileupExtractor::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
PileupExtractor::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
PileupExtractor::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PileupExtractor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PileupExtractor);
