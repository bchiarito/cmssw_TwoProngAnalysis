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
    int run;
    int LS;
    int evnum;
  };

  std::string eventInfoBranchDefString("run/I:LS:evnum");

  void FillEventInfo(eventInfo_t &eventInfo,const edm::Event& iEvent) {
    
    eventInfo.run = iEvent.id().run();
    eventInfo.LS = iEvent.id().luminosityBlock();
    eventInfo.evnum = iEvent.id().event();
    
  }


  // vertex info

  struct vtxInfo_t{
    int Nvtx; // number of reco vertices in event
    // but I will only keep the info of the best two, I think
    double vx;
    double vy;
    double vz;
    int isFake;
    int Ntracks;
    double sumPtTracks;
    // ************* chi2 or some other quality criteria *****************
    // these are the two vars used in GoodVertexFilter module
    double ndof;
    double d0; 

  };

  std::string vtxInfoBranchDefString("Nvtx/I:vx/D:vy:vz:isFake/I:Ntracks/I:sumPtTracks/D:ndof:d0");

  void FillVertexInfo(vtxInfo_t &vtxInfo, const reco::Vertex *vertex) {

     vtxInfo.vx = vertex->x();
     vtxInfo.vy = vertex->y();
     vtxInfo.vz = vertex->z();
     vtxInfo.isFake = vertex->isFake(); 
     vtxInfo.Ntracks = vertex->tracksSize();
     vtxInfo.ndof = vertex->ndof();
     vtxInfo.d0 = vertex->position().rho();

     // special note: SumPtTracks is not a member of vertex class
     // so it needs to be calculated separately by looping over track collection
     // for now, since this is already done in main analyser, just make
     // sure to fill this in the analyser code 
     vtxInfo.sumPtTracks = -9999.99;

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
