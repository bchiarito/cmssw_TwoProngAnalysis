// Event and Vertex info

#ifndef EVENT_AND_VERTEX_INFO_INC
#define EVENT_AND_VERTEX_INFO_INC

#include <string>

#include "FWCore/Framework/interface/Event.h"

//for vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// for beamspot
#include "DataFormats/BeamSpot/interface/BeamSpot.h"


namespace ExoDiPhotons
{

  // event info 

  struct eventInfo_t{
    Long64_t run;
    Long64_t LS;
    Long64_t evnum;
    Long64_t processid;
    Float_t pthat;
    Float_t alphaqcd;
    Float_t alphaqed;
    Float_t qscale;
    Float_t weight;

  };

  std::string eventInfoBranchDefString("run/L:LS:evnum:processid:pthat/F:alphaqcd:alphaqed:qscale:weight");

  void FillEventInfo(eventInfo_t &eventInfo,const edm::Event& iEvent) {
    
    eventInfo.run = (Long64_t)iEvent.id().run();
    eventInfo.LS = (Long64_t)iEvent.id().luminosityBlock();
    eventInfo.evnum = (Long64_t)iEvent.id().event();
    
  }

  void InitEventInfo(eventInfo_t &eventInfo, float value) {
    
    eventInfo.run = (Long64_t)value;
    eventInfo.LS = (Long64_t)value;
    eventInfo.evnum = (Long64_t)value;
    
    eventInfo.processid = (Long64_t)value;
    eventInfo.pthat = value;
    eventInfo.alphaqcd = value;
    eventInfo.alphaqed = value;
    eventInfo.qscale = value;
    eventInfo.weight = value;

  }

  // vertex info

  struct vtxInfo_t{
    Double_t vx;
    Double_t vy;
    Double_t vz;

    Double_t sumPtTracks;
    // ************* chi2 or some other quality criteria *****************
    // these are the two vars used in GoodVertexFilter module
    Double_t ndof;
    Double_t d0; 

    Int_t Nvtx; // number of reco vertices in event
    // but I will only keep the info of the best two, I think
    Int_t Ntracks;
    Bool_t isFake;

  };

  std::string vtxInfoBranchDefString("vx/D:vy:vz:sumPtTracks/D:ndof/D:d0/D:Nvtx/I:Ntracks/I:isFake/O");


  // simple function to get sumPtTracks for vtx
  double calcVtxSumPtTracks(const reco::Vertex *vtx) {

    // loop over assoc tracks to get sum pt
    double sumPtTracks = 0.0;
    for(reco::Vertex::trackRef_iterator vtxTracks=vtx->tracks_begin(); vtxTracks!=vtx->tracks_end();vtxTracks++) {

      sumPtTracks += (**vtxTracks).pt();
    }
    return sumPtTracks;
  }


  void FillVertexInfo(vtxInfo_t &vtxInfo, const reco::Vertex *vertex) {

    vtxInfo.vx = (double) vertex->x();
    vtxInfo.vy = (double)  vertex->y();
    vtxInfo.vz = (double) vertex->z();
    vtxInfo.isFake =  vertex->isFake(); 
    vtxInfo.Ntracks =  vertex->tracksSize();
    vtxInfo.ndof = (double) vertex->ndof();
    vtxInfo.d0 = (double) vertex->position().rho();

    // now I have a function to get sumPtTracks
    vtxInfo.sumPtTracks = calcVtxSumPtTracks(vertex);

  }


  // I also use my own compare function for sorting vertices in the analyser
  bool sortVertices(const reco::Vertex &vtx1, const reco::Vertex &vtx2) 
  {
    // sort by Ntracks, with highest first
    //    return(vtx1.tracksSize()>=vtx2.tracksSize());
    
    // or by TrackSumPt
    return(calcVtxSumPtTracks(&vtx1)>=calcVtxSumPtTracks(&vtx2));
    
  }





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
  
  std::string beamSpotInfoBranchDefString("x0/D:y0:z0:sigmaZ:x0error:y0error:z0error:sigmaZ0error");


  void FillBeamSpotInfo(beamSpotInfo_t &beamSpotInfo, reco::BeamSpot beamSpot) {

    beamSpotInfo.x0 = beamSpot.x0();
    beamSpotInfo.y0 = beamSpot.y0();
    beamSpotInfo.z0 = beamSpot.z0();
    beamSpotInfo.sigmaZ = beamSpot.sigmaZ();
    beamSpotInfo.x0error = beamSpot.x0Error();
    beamSpotInfo.y0error = beamSpot.y0Error();
    beamSpotInfo.z0error = beamSpot.z0Error();
    beamSpotInfo.sigmaZ0error = beamSpot.sigmaZ0Error();
  }


  // initialise function also?




} // end of namespace


#endif
