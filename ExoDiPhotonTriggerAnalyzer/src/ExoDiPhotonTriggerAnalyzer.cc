// -*- C++ -*-
//
// Package:    ExoDiPhotonTriggerAnalyzer
// Class:      ExoDiPhotonTriggerAnalyzer
// 
/**\class ExoDiPhotonTriggerAnalyzer ExoDiPhotonTriggerAnalyzer.cc DiPhotonAnalysis/ExoDiPhotonTriggerAnalyzer/src/ExoDiPhotonTriggerAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Otman Charaf,42 1-015,+41227662353,
//         Created:  Sat Feb 23 18:03:42 CET 2013
// $Id$
//
//


// system include files
#include <memory>

#include <algorithm>
#include <vector>
#include <utility>  // for std::pair
#include "TClonesArray.h"
#include "TVector3.h"

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
#include "TString.h"

// geometry
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
//#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
//#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
//#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
//#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"


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

// new CommonClasses approach
// these objects are all in the namespace 'ExoDiPhotons'
#include "DiPhotonAnalysis/CommonClasses/interface/RecoPhotonInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/TriggerInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/PhotonID.h"
#include "DiPhotonAnalysis/CommonClasses/interface/EventAndVertexInfo.h"
#include "DiPhotonAnalysis/CommonClasses/interface/DiphotonInfo.h"


//new for PU gen
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// new for LumiReweighting 
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

//new for PF ID definition
#include "DiPhotonAnalysis/CommonClasses/interface/PFPhotonID.h"

//for conversion safe electron veto
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//new for PFIsolation code
#include "EGamma/EGammaAnalysisTools/src/PFIsolationEstimator.cc"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

using namespace std;


//
// class declaration
//

class ExoDiPhotonTriggerAnalyzer : public edm::EDAnalyzer {
public:
  explicit ExoDiPhotonTriggerAnalyzer(const edm::ParameterSet&);
  ~ExoDiPhotonTriggerAnalyzer();

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

  edm::InputTag      fPhotonTag;       //select photon collection 
  double             fMin_pt;          // min pt cut (photons)
  edm::InputTag      fHltInputTag;     // hltResults
  edm::InputTag      fL1InputTag;      // L1 results
  edm::InputTag      fRho25Tag;  
  edm::InputTag      fpileupCollectionTag;         
  edm::LumiReWeighting    LumiWeights;      
  
  bool               fkRemoveSpikes;   // option to remove spikes before filling tree
  bool               fkRequireTightPhotons;  // option to require tight photon id in tree
  bool               fkRequireGenEventInfo;  // generated information for RS graviton files
  bool               fisMC;  //option to decide if MC or Data     
  string             fPUMCFileName;
  string             fPUDataFileName;
  string             fPUDataHistName;
  string             fPUMCHistName;   
  string             fPFIDCategory;   
  string             fIDMethod;   

  double fRho25;

  //for PFIsolation Code
  PFIsolationEstimator isolator03;

  ExoDiPhotons::recoPhotonInfo_t fRecoPhotonInfo1; // leading photon 
  ExoDiPhotons::recoPhotonInfo_t fRecoPhotonInfo2; // second photon
  ExoDiPhotons::diphotonInfo_t fDiphotonInfo;

  ExoDiPhotons::hltTrigInfo_t fHLTInfo;

  TH1F *fPtPhotonRaw;
  TH1F *fNumTotalEvents;
  TH1F *fNumTotalEventsPassHLT;
  TH1F *fPtPhotonBarrelCut;
  TH1F *fPtPhotonBarrelCutPtCut;
  TH1F *fPtPhotonBarrelCutPtCutTightCut;
  TH1F *fPtTwoPhotonBarrelCutPtCutTightCutMass;
  TH1F *fPtTwoPhotonBarrelCutPtCutTightCutMassCut;

  TH1F *fMassAcceptance;
  TH1F *fMassBarrelBarrel;
  TH1F *fMassBarrelBarrelPtCut;
  TH1F *fMassBarrelBarrelPtCutTightCut;
  TH1F *fMassBarrelBarrelPtCutTightCutMassCut;

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
ExoDiPhotonTriggerAnalyzer::ExoDiPhotonTriggerAnalyzer(const edm::ParameterSet& iConfig)
  : fPhotonTag(iConfig.getUntrackedParameter<edm::InputTag>("photonCollection")),
    fMin_pt(iConfig.getUntrackedParameter<double>("ptMin")),
    fHltInputTag(iConfig.getUntrackedParameter<edm::InputTag>("hltResults")),
    fL1InputTag(iConfig.getUntrackedParameter<edm::InputTag>("L1Results")),
    fRho25Tag(iConfig.getParameter<edm::InputTag>("rho25Correction")),
    fpileupCollectionTag(iConfig.getUntrackedParameter<edm::InputTag>("pileupCorrection")),
    fkRemoveSpikes(iConfig.getUntrackedParameter<bool>("removeSpikes")),
    fkRequireTightPhotons(iConfig.getUntrackedParameter<bool>("requireTightPhotons")),
    fkRequireGenEventInfo(iConfig.getUntrackedParameter<bool>("requireGenEventInfo")),
    fisMC(iConfig.getUntrackedParameter<bool>("isMC")),
    fPUMCFileName(iConfig.getUntrackedParameter<string>("PUMCFileName")),
    fPUDataFileName(iConfig.getUntrackedParameter<string>("PUDataFileName")), 
    fPUDataHistName(iConfig.getUntrackedParameter<string>("PUDataHistName")),
    fPUMCHistName(iConfig.getUntrackedParameter<string>("PUMCHistName")),
    fPFIDCategory(iConfig.getUntrackedParameter<string>("PFIDCategory")),
    fIDMethod(iConfig.getUntrackedParameter<string>("IDMethod"))
{
  //now do what ever initialization is needed

  edm::Service<TFileService> fs;
  fNumTotalEvents = fs->make<TH1F>("fNumTotalEvents","fNumTotalEvents",2000,0.,2000.);
  fNumTotalEventsPassHLT = fs->make<TH1F>("fNumTotalEventsPassHLT","fNumTotalEventsPassHLT",2000,0.,2000.);
  fPtPhotonRaw = fs->make<TH1F>("fPtPhotonRaw","fPtPhotonRaw",2000,0.,2000.);
  fPtPhotonBarrelCut = fs->make<TH1F>("fPtPhotonBarrelCut","fPtPhotonBarrelCut",2000,0.,2000.);
  fPtPhotonBarrelCutPtCut = fs->make<TH1F>("fPtPhotonBarrelCutPtCut","fPtPhotonBarrelCutPtCut",2000,0.,2000.);
  fPtPhotonBarrelCutPtCutTightCut = fs->make<TH1F>("fPtPhotonBarrelCutPtCutTightCut","fPtPhotonBarrelCutPtCutTightCut",2000,0.,2000.);
  fPtTwoPhotonBarrelCutPtCutTightCutMass = fs->make<TH1F>("fPtTwoPhotonBarrelCutPtCutTightCutMass","fPtTwoPhotonBarrelCutPtCutTightCutMass",4000,0.,4000.);
  fPtTwoPhotonBarrelCutPtCutTightCutMassCut = fs->make<TH1F>("fPtTwoPhotonBarrelCutPtCutTightCutMassCut","fPtTwoPhotonBarrelCutPtCutTightCutMassCut",4000,0.,4000.);

  fMassAcceptance = fs->make<TH1F>("fMassAcceptance","fMassAcceptance",4000,0.,4000.);
  fMassBarrelBarrel = fs->make<TH1F>("fMassBarrelBarrel","fMassBarrelBarrel",4000,0.,4000.);
  fMassBarrelBarrelPtCut = fs->make<TH1F>("fMassBarrelBarrelPtCut","fMassBarrelBarrelPtCut",4000,0.,4000.);
  fMassBarrelBarrelPtCutTightCut = fs->make<TH1F>("fMassBarrelBarrelPtCutTightCut","fMassBarrelBarrelPtCutTightCut",4000,0.,4000.);
  fMassBarrelBarrelPtCutTightCutMassCut = fs->make<TH1F>("fMassBarrelBarrelPtCutTightCutMassCut","fMassBarrelBarrelPtCutTightCutMassCut",4000,0.,4000.);

  //new PFIsolation code
  isolator03.initializePhotonIsolation(kTRUE);
  isolator03.setConeSize(0.3);

}


ExoDiPhotonTriggerAnalyzer::~ExoDiPhotonTriggerAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ExoDiPhotonTriggerAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  TString CategoryPFID(fPFIDCategory.c_str());
  TString MethodID(fIDMethod.c_str());

  edm::Handle<double> rho25Handle;
  iEvent.getByLabel(fRho25Tag, rho25Handle);

  if (!rho25Handle.isValid()){
    cout<<"rho25 not found"<<endl;
    return;
  }
        
  fRho25 = *(rho25Handle.product());

  // get the vertex collection
  Handle<reco::VertexCollection> vertexColl;
  iEvent.getByLabel("offlinePrimaryVertices",vertexColl);
   
  if(!vertexColl.isValid()) {
    cout << "Vertex collection empty! Bailing out!" <<endl;
    return;
  }
  //get the reference to 1st vertex for use in fGetIsolation
  //for PFIsolation calculation
  reco::VertexRef firstVtx(vertexColl,0);
   
  //for conversion safe electron veto
  edm::Handle<reco::ConversionCollection> hConversions;
  iEvent.getByLabel("allConversions", hConversions);
  
  edm::Handle<reco::GsfElectronCollection> hElectrons;
  iEvent.getByLabel("gsfElectrons", hElectrons);
  
  //for PFIsolation code
  Handle<PFCandidateCollection> pfCandidatesColl;
  iEvent.getByLabel("particleFlow",pfCandidatesColl);
  const PFCandidateCollection * pfCandidates = pfCandidatesColl.product();

  // get offline beam spot
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
   
  // get the photon collection
  Handle<reco::PhotonCollection> photonColl;
  iEvent.getByLabel(fPhotonTag,photonColl);

  // If photon collection is empty, exit
  if (!photonColl.isValid()) {
    cout << "No Photons! Move along, there's nothing to see here .." <<endl;
    return;
  }
   
  // get the trig info

  //trig results
  Handle<TriggerResults> hltResultsHandle;
  iEvent.getByLabel(fHltInputTag,hltResultsHandle);
   
  if(!hltResultsHandle.isValid()) {
    cout << "HLT results not valid!" <<endl;
    cout << "Couldnt find TriggerResults with input tag " << fHltInputTag << endl;
    return;
  }

  const TriggerResults *hltResults = hltResultsHandle.product();
  const TriggerNames & hltNames = iEvent.triggerNames(*hltResults);
  // now we just use the FillHLTInfo() function from TrigInfo.h:
  ExoDiPhotons::FillHLTInfo(fHLTInfo,hltResults,hltNames);

  bool eventpasshlt = (fHLTInfo.HLT_DoublePhoton70_v3 == 1) || (fHLTInfo.HLT_DoublePhoton70_v4 == 1) || (fHLTInfo.HLT_DoublePhoton70_v5 == 1) || (fHLTInfo.HLT_DoublePhoton70_v6 == 1);
  //cout<<"eventpasshlt "<<eventpasshlt<<endl;

  //std::vector<reco::Photon> selectedPhotons; 
  std::vector<reco::Photon> AcceptancePhotons; 
  std::vector<reco::Photon> BarrelBarrelPhotons; 
  std::vector<reco::Photon> BarrelBarrelPtCutPhotons; 
  std::vector<reco::Photon> selectedPhotons; 

  fNumTotalEvents->Fill(-50.);

  if(!eventpasshlt) return;
  fNumTotalEventsPassHLT->Fill(-50.);

  for(reco::PhotonCollection::const_iterator recoPhoton = photonColl->begin(); recoPhoton!=photonColl->end(); recoPhoton++) {

    std::vector<double> effareas = ExoDiPhotons::EffectiveAreas(&(*recoPhoton));
    double pfisoall = isolator03.fGetIsolation(&(*recoPhoton),pfCandidates,firstVtx,vertexColl);
    double rhocorPFIsoCH = max(isolator03.getIsolationCharged()-fRho25*effareas[0],0.);
    double rhocorPFIsoNH = max(isolator03.getIsolationNeutral()-fRho25*effareas[1],0.);
    double rhocorPFIsoPH = max(isolator03.getIsolationPhoton()-fRho25*effareas[2],0.);
    //and we also have to test the conversion safe electron veto
    bool passelecveto = !ConversionTools::hasMatchedPromptElectron(recoPhoton->superCluster(), hElectrons, hConversions, beamSpot.position());

    fPtPhotonRaw->Fill(recoPhoton->pt());

    if( (fabs(recoPhoton->caloPosition().eta()) > 2.5) || ExoDiPhotons::isGapPhoton(&(*recoPhoton)) ) continue;
    AcceptancePhotons.push_back(*recoPhoton);

    if( !(ExoDiPhotons::isBarrelPhoton(&(*recoPhoton)) && !ExoDiPhotons::isGapPhoton(&(*recoPhoton)) ) ) continue;
    BarrelBarrelPhotons.push_back(*recoPhoton);
    fPtPhotonBarrelCut->Fill(recoPhoton->pt());

    if(recoPhoton->pt() < fMin_pt) continue;
    BarrelBarrelPtCutPhotons.push_back(*recoPhoton);
    fPtPhotonBarrelCutPtCut->Fill(recoPhoton->pt());

    //Now we choose which ID to use (PF or Det)
    if(MethodID.Contains("Detector")){
      if(ExoDiPhotons::isTightPhoton(&(*recoPhoton),fRho25) && !ExoDiPhotons::isGapPhoton(&(*recoPhoton)) && !ExoDiPhotons::isASpike(&(*recoPhoton))  ) {
	//	    if( !ExoDiPhotons::isASpike(&(*recoPhoton))  ) {   
	selectedPhotons.push_back(*recoPhoton);
	fPtPhotonBarrelCutPtCutTightCut->Fill(recoPhoton->pt());
      }
    }
    else if(MethodID.Contains("ParticleFlow")){
      if(ExoDiPhotons::isPFTightPhoton(&(*recoPhoton),rhocorPFIsoCH,rhocorPFIsoNH,rhocorPFIsoPH,passelecveto,CategoryPFID) && 
	 !ExoDiPhotons::isGapPhoton(&(*recoPhoton)) && 
	 !ExoDiPhotons::isASpike(&(*recoPhoton))  ) {
	//	    if( !ExoDiPhotons::isASpike(&(*recoPhoton))  ) {   
	selectedPhotons.push_back(*recoPhoton);
	fPtPhotonBarrelCutPtCutTightCut->Fill(recoPhoton->pt());
      }
    }

  }//end of loop over photon collection

  //sort all vector of photons
  sort(AcceptancePhotons.begin(),AcceptancePhotons.end(),ExoDiPhotons::comparePhotonsByPt);
  sort(BarrelBarrelPhotons.begin(),BarrelBarrelPhotons.end(),ExoDiPhotons::comparePhotonsByPt);
  sort(BarrelBarrelPtCutPhotons.begin(),BarrelBarrelPtCutPhotons.end(),ExoDiPhotons::comparePhotonsByPt);
  sort(selectedPhotons.begin(),selectedPhotons.end(),ExoDiPhotons::comparePhotonsByPt);

  if(AcceptancePhotons.size() >= 2){
    ExoDiPhotons::FillDiphotonInfo(fDiphotonInfo,&AcceptancePhotons[0],&AcceptancePhotons[1]);
    fMassAcceptance->Fill(fDiphotonInfo.Minv);
  }

  if(BarrelBarrelPhotons.size() >= 2){
    ExoDiPhotons::FillDiphotonInfo(fDiphotonInfo,&BarrelBarrelPhotons[0],&BarrelBarrelPhotons[1]);
    fMassBarrelBarrel->Fill(fDiphotonInfo.Minv);
  }

  if(BarrelBarrelPtCutPhotons.size() >= 2){
    ExoDiPhotons::FillDiphotonInfo(fDiphotonInfo,&BarrelBarrelPtCutPhotons[0],&BarrelBarrelPtCutPhotons[1]);
    fMassBarrelBarrelPtCut->Fill(fDiphotonInfo.Minv);
  }

  if(selectedPhotons.size() >= 2){
    ExoDiPhotons::FillDiphotonInfo(fDiphotonInfo,&selectedPhotons[0],&selectedPhotons[1]);
    fMassBarrelBarrelPtCutTightCut->Fill(fDiphotonInfo.Minv);
    if(fDiphotonInfo.Minv > 200.) fMassBarrelBarrelPtCutTightCutMassCut->Fill(fDiphotonInfo.Minv);
  }



}//end of method analyze


// ------------ method called once each job just before starting event loop  ------------
void 
ExoDiPhotonTriggerAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ExoDiPhotonTriggerAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
ExoDiPhotonTriggerAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ExoDiPhotonTriggerAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ExoDiPhotonTriggerAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ExoDiPhotonTriggerAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ExoDiPhotonTriggerAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ExoDiPhotonTriggerAnalyzer);
