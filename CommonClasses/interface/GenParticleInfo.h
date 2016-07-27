#ifndef GEN_PARTICLE_INFO_INC
#define GEN_PARTICLE_INFO_INC

//********************************************************************
// Definition of a struct that can be used for storing reco charged type decaying eta info
// in a tree, from different analysers
// Also includes a Fill function to fill the struct from the appropriate objects
// and a string that can be used to define the tree branch
// 
//  $Id: RecoTwoProngInfo.h,v 1.00 2016 16:26:48 charaf Exp $
// 
//********************************************************************

#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

// geometry
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

namespace ExoDiPhotons
{

  struct genParticleInfo_t {
    Double_t pt;
    Double_t phi;
    Double_t eta;
    Double_t mass;
    Double_t px;
    Double_t py;
    Double_t pz;
    Double_t energy;
  };

  // also include a string that can be used to define the tree branch
  // obviously this needs to be kept up-to-date with the struct definition
  // but now at least this only needs to be done here in this file, 
  // rather than in each individual analyser 
  std::string genParticleInfoBranchDefString("pt/D:phi/D:eta/D:mass/D:px/D:py/D:pz/D:energy/D");

  // also want a Fill function, that can fill the struct values from the appropriate objects
  // again, so that all editing only needs to be done here in this file
  void FillGenParticleInfo(genParticleInfo_t &genparticleinfo, reco::GenParticle genparticle) {
    genparticleinfo.pt = genparticle.pt();
    genparticleinfo.eta = genparticle.eta();
    genparticleinfo.phi = genparticle.phi();
    genparticleinfo.mass = genparticle.mass();

    genparticleinfo.px = genparticle.px();
    genparticleinfo.py = genparticle.py();
    genparticleinfo.pz = genparticle.pz();
    genparticleinfo.energy = genparticle.energy();
  }

  void InitGenParticleInfo(genParticleInfo_t &genparticleinfo) {
    genparticleinfo.pt = -99.9;
    genparticleinfo.eta = -99.9;
    genparticleinfo.phi = -99.9;
    genparticleinfo.mass = -99.9;

    genparticleinfo.px = -99.9;
    genparticleinfo.py = -99.9;
    genparticleinfo.pz = -99.9;
    genparticleinfo.energy = -99.9;
  }
  bool compareGenParticlesByPt(const reco::GenParticle &photon1, const reco::GenParticle &photon2) {

    // sorts such that highest pt photon first
    return(photon1.pt()>=photon2.pt());
  }

} //end of namespace

#endif
