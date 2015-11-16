// -*- C++ -*-
//
// Package:   FilterOutScraping
// Class:     FilterOutScraping
//
// Original Author:  Luca Malgeri

// Adapted for use with MiniAOD by Steven Kaplan (skaplan@cern.ch)

#ifndef ExoDiphotonTrackFilter_H
#define ExoDiphotonTrackFilter_H

// system include files
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
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

// for plotting
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//
// class declaration
//


class ExoDiPhotonTrackFilter : public edm::EDFilter {
public:
  explicit ExoDiPhotonTrackFilter( const edm::ParameterSet & );
  ~ExoDiPhotonTrackFilter();
  
private:
  virtual bool filter ( edm::Event &, const edm::EventSetup&) override;
  
  bool applyfilter;
  bool debugOn;
  double thresh;
  unsigned int numtrack;
  edm::EDGetTokenT<pat::PackedCandidateCollection> cands_;

  // plots for debugging
  edm::Service<TFileService> fs;
  TH1D* hNumTracks;
  TH1D* hNumHighPurity;
  TH1D* hAccepted;
  TH1D* hFraction;
};

#endif