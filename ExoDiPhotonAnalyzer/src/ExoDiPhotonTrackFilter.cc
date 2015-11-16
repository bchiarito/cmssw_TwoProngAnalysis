// -*- C++ -*-
//
// Package:   BeamSplash
// Class:     BeamSPlash
//
//
// Original Author:  Luca Malgeri

// Adapted for use with MiniAOD by Steven Kaplan (skaplan@cern.ch)

#include <memory>
#include <vector>
#include <map>
#include <set>

// user include files

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DiPhotonAnalysis/ExoDiPhotonAnalyzer/interface/ExoDiPhotonTrackFilter.h"


using namespace edm;
using namespace std;

ExoDiPhotonTrackFilter::ExoDiPhotonTrackFilter(const edm::ParameterSet& iConfig)
{
  
  applyfilter = iConfig.getUntrackedParameter<bool>("applyfilter",true);
  debugOn     = iConfig.getUntrackedParameter<bool>("debugOn",false);
  thresh =  iConfig.getUntrackedParameter<double>("thresh",0.2);
  numtrack = iConfig.getUntrackedParameter<unsigned int>("numtrack",10);
  cands_ = consumes<pat::PackedCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("src",edm::InputTag("packedPFCandidates")));

  // book histograms
  if (debugOn){
    hNumTracks = fs->make<TH1D>("hNumTracks","",701,-0.5,700.5);
    hNumHighPurity = fs->make<TH1D>("hNumHighPurity","",701,-0.5,700.5);
    hAccepted = fs->make<TH1D>("hAccepted","",2,-0.5,1.5);
    hFraction = fs->make<TH1D>("hFraction","",102,-0.01,1.01);
  }
}

ExoDiPhotonTrackFilter::~ExoDiPhotonTrackFilter()
{
}

bool ExoDiPhotonTrackFilter::filter( edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  bool accepted = false;
  //float fraction = 0;  
  // get GeneralTracks collection

  edm::Handle<pat::PackedCandidateCollection> packedCandsRef;
  iEvent.getByToken(cands_,packedCandsRef);    
  const pat::PackedCandidateCollection* packedCandsColl = packedCandsRef.product();

  //std::cout << "Total Number of Tracks " << tkColl->size() << std::endl;
  
  unsigned int numHighPurity = 0;
  unsigned int numTotalTracks = 0;
  for(pat::PackedCandidateCollection::const_iterator iCand = packedCandsColl->begin(); iCand != packedCandsColl->end(); ++iCand){

    // first see if the candidate is charged
    int candCharge = iCand->charge();
    if (candCharge == 0) continue;
    else{
      numTotalTracks++;
      bool isHighPurity = iCand->trackHighPurity();
      if (isHighPurity) numHighPurity++;

    }

  }
  double fraction = ( (float)numHighPurity ) / ( (float)numTotalTracks );

  // in the scenario of numTotalTracks < numtrack, acccept the event anyway
  if (numTotalTracks < numtrack) accepted = true;
  else if (fraction > thresh) accepted = true;

  if (debugOn) {
    // int ievt = iEvent.id().event();
    // int irun = iEvent.id().run();
    // int ils = iEvent.luminosityBlock();
    // int bx = iEvent.bunchCrossing();
    
    // std::cout << "ExoDiPhotonTrackFilter_debug: Run " << irun << " Event " << ievt << " Lumi Block " << ils << " Bunch Crossing " << bx << " Fraction " << fraction << " NTracks " << numTotalTracks << " Accepted " << accepted << std::endl;
  
  //fill debug histograms
  hNumTracks->Fill(numTotalTracks);
  hNumHighPurity->Fill(numHighPurity);
  hAccepted->Fill(accepted);
  hFraction->Fill(fraction);
  }
 
  if (applyfilter)
    return accepted;
  else
    return true;

}

//define this as a plug-in
DEFINE_FWK_MODULE(ExoDiPhotonTrackFilter);
