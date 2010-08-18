#ifndef RECO_PHOTON_INFO_INC
#define RECO_PHOTON_INFO_INC

//********************************************************************
// Definition of a struct that can be used for storing reco photon info
// in a tree, from different analysers
// Also includes a Fill function to fill the struct from the appropriate objects
// and a string that can be used to define the tree branch
// 
//  $Id: RecoPhotonInfo.h,v 1.3 2010/08/18 12:06:16 torimoto Exp $
// 
//********************************************************************

#include <string>

// for ecal
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"


//for photons
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"


namespace ExoDiPhotons
{

  struct recoPhotonInfo_t {
    double pt;
    double eta;
    double phi;
    // position in ECAL  - caloPosition;
    double detEta;
    double detPhi; // clearly should be identical to phi, so am using as sort of cross-check 
    // which detector channel, specifically?
    //useful to cross check channel masking and look for problem channels, etc I expect
    int detId;
    int iEtaY; // iEta if EB, iY if EE
    int iPhiX; // iPhi if EB, iX if EE


    //double check the vertex assigned to the photon candidate
    double vx;
    double vy;
    double vz;
  

    //shower shape variables
    double r9;
    double sigmaIetaIeta;
    double sigmaEtaEta;
    double maxEnergyXtal;
    // eNxN ...
    double e1x5;
    double e2x5;
    double e3x3;
    double e5x5;
    double r1x5;
    double r2x5;
  
    // swiss cross and other spike-related flags, eg kOutOfTime
    double swisscross;
    double eMax; // believe this is same as maxEnergyCrystal, but its a different getter, so let's just put it in to be on the safe side
    double eLeft;
    double eRight;
    double eTop;
    double eBottom;
    double eSecond; // second highest rec hit in the cluster (not sure if required to be within 3x3 or not?)

    // ecal severity level
    int severityLevel;
    int recHitFlag;

    // rec hit timing
    double maxRecHitTime;

    double hadOverEm;  
    // note also that two hadronic depths are available
    double hadDepth1OverEm; 
    double hadDepth2OverEm; 
 

    //isolation variables
    // these have different options: cone size, hollowness, etc
    // must include them all!
    float hcalIso04;
    float hcalIso03;
    float ecalIso04;
    float ecalIso03;
    float trkIsoSumPtHollow04;
    float trkIsoSumPtSolid04;
    int trkIsoNtrksHollow04;
    int trkIsoNtrksSolid04;
    float trkIsoSumPtHollow03;
    float trkIsoSumPtSolid03;
    int trkIsoNtrksHollow03;
    int trkIsoNtrksSolid03;

    // es ratio
    float esRatio;

    // supercluster info
    double scRawEnergy;
    double scPreshowerEnergy;
    double scPhiWidth;
    double scEtaWidth;
    int scNumBasicClusters; // number of basic clusters comprising superCluster

    // seed cluster info
  
  
    //fiducial flags
    bool isEB;//Photon is in EB
    bool isEE;//Photon is in EE
    bool isEBEtaGap;//Photon is in supermodule/supercrystal eta gap in EB
    bool isEBPhiGap;//Photon is in supermodule/supercrystal phi gap in EB
    bool isEERingGap;//Photon is in crystal ring gap in EE
    bool isEEDeeGap;//Photon is in crystal dee gap in EE
    bool isEBEEGap;//Photon is in border between EB and EE.

    // pixel seed match?
    bool hasPixelSeed;
    // note to self: weird problems with this var in middle of struct - try at end
    
  };


  // also include a string that can be used to define the tree branch
  // obviously this needs to be kept up-to-date with the struct definition
  // but now at least this only needs to be done here in this file, 
  // rather than in each individual analyser 
  std::string recoPhotonBranchDefString("pt/D:eta:phi:detEta:detPhi:detId/I:iEtaY/I:iPhiX/I:vx/D:vy:vz:r9:sigmaIetaIeta:sigmaEtaEta:maxEnergyXtal:e1x5:e2x5:e3x3:e5x5:r1x5:r2x5:swisscross:eMax:eLeft:eRight:eTop:eBottom:eSecond:severityLevel/I:recHitFlag/I:maxRecHitTime/D:hadOverEm:hadDepth1OverEm:hadDepth2OverEm:hcalIso04/f:hcalIso03/f:ecalIso04:ecalIso03:trkIsoSumPtHollow04:trkIsoSumPtSolid04:trkIsoNtrksHollow04/I:trkIsoNtrksSolid04/I:trkIsoSumPtHollow03/f:trkIsoSumPtSolid03/f:trkIsoNtrksHollow03/I:trkIsoNtrksSolid03/I:esRatio/f:scRawEnergy/D:scPreshowerEnergy:scPhiWidth:scEtaWidth:scNumBasicClusters/I:isEB/O:isEE:isEBEtaGap:isEBPhiGap:isEERingGap:isEEDeeGap:isEBEEGap:hasPixelSeed");

  

  // also want a Fill function, that can fill the struct values from the appropriate objects
  // again, so that all editing only needs to be done here in this file

  void FillRecoPhotonInfo(recoPhotonInfo_t &recoPhotonInfo, const reco::Photon *photon) {
    
    recoPhotonInfo.pt = photon->et();
    recoPhotonInfo.eta = photon->eta();
    recoPhotonInfo.phi = photon->phi();
    
    recoPhotonInfo.detEta = photon->caloPosition().eta();
    recoPhotonInfo.detPhi = photon->caloPosition().phi();
    
    //   detId and crystal eta/phi
    // these use lazyTools and so need to be filled inside the analyser itself still
    recoPhotonInfo.detId = -99999;
    recoPhotonInfo.iEtaY = -99999;
    recoPhotonInfo.iPhiX = -99999;


    // since photon inherits from LeafCandidate, we can get the vertex position
    // that is associated with the photon:
    recoPhotonInfo.vx = photon->vx();
    recoPhotonInfo.vy = photon->vy();
    recoPhotonInfo.vz = photon->vz();


     recoPhotonInfo.r9 = photon->r9();
     recoPhotonInfo.sigmaIetaIeta = photon->sigmaIetaIeta();
     recoPhotonInfo.sigmaEtaEta = photon->sigmaEtaEta();
     recoPhotonInfo.maxEnergyXtal = photon->maxEnergyXtal();

     recoPhotonInfo.e1x5 = photon->e1x5();
     recoPhotonInfo.e2x5 = photon->e2x5();
     recoPhotonInfo.e3x3 = photon->e3x3();
     recoPhotonInfo.e5x5 = photon->e5x5();
     recoPhotonInfo.r1x5 = photon->r1x5();
     recoPhotonInfo.r2x5 = photon->r2x5();

     // swiss cross and related use lazyTools, so need to be filled in analyser itself
     recoPhotonInfo.swisscross = -999999.99;
     recoPhotonInfo.eMax = -999999.99; 
     recoPhotonInfo.eLeft = -999999.99;
     recoPhotonInfo.eRight = -999999.99;
     recoPhotonInfo.eTop = -999999.99;
     recoPhotonInfo.eBottom = -999999.99;
     recoPhotonInfo.eSecond = -999999.99;

     // official ecal rec hit severity level
     recoPhotonInfo.severityLevel = -999;
     recoPhotonInfo.recHitFlag = -999;
     // time of highest energy rec hit
     recoPhotonInfo.maxRecHitTime = -9999999.99;


     recoPhotonInfo.hadOverEm = photon->hadronicOverEm();
     recoPhotonInfo.hadDepth1OverEm = photon->hadronicDepth1OverEm();
     recoPhotonInfo.hadDepth2OverEm = photon->hadronicDepth2OverEm();

     
     recoPhotonInfo.ecalIso04 = photon->ecalRecHitSumEtConeDR04();
     recoPhotonInfo.ecalIso03 = photon->ecalRecHitSumEtConeDR03();
     recoPhotonInfo.hcalIso04 = photon->hcalTowerSumEtConeDR04();
     recoPhotonInfo.hcalIso03 = photon->hcalTowerSumEtConeDR03();

     recoPhotonInfo.trkIsoSumPtHollow04 = photon->trkSumPtHollowConeDR04();
     recoPhotonInfo.trkIsoSumPtSolid04 = photon->trkSumPtSolidConeDR04();
     recoPhotonInfo.trkIsoNtrksHollow04 = photon->nTrkHollowConeDR04();
     recoPhotonInfo.trkIsoNtrksSolid04 = photon->nTrkSolidConeDR04();
     
     recoPhotonInfo.trkIsoSumPtHollow03 = photon->trkSumPtHollowConeDR03();
     recoPhotonInfo.trkIsoSumPtSolid03 = photon->trkSumPtSolidConeDR03();
     recoPhotonInfo.trkIsoNtrksHollow03 = photon->nTrkHollowConeDR03();
     recoPhotonInfo.trkIsoNtrksSolid03 = photon->nTrkSolidConeDR03();

     recoPhotonInfo.esRatio = -999.;

     recoPhotonInfo.hasPixelSeed = photon->hasPixelSeed();


     recoPhotonInfo.isEB        = photon->isEB();        
     recoPhotonInfo.isEE	 = photon->isEE();	 
     recoPhotonInfo.isEBEtaGap	 = photon->isEBEtaGap();	 
     recoPhotonInfo.isEBPhiGap	 = photon->isEBPhiGap();	 
     recoPhotonInfo.isEERingGap = photon->isEERingGap(); 
     recoPhotonInfo.isEEDeeGap	 = photon->isEEDeeGap();	 
     recoPhotonInfo.isEBEEGap   = photon->isEBEEGap();
     
     recoPhotonInfo.scRawEnergy = photon->superCluster()->rawEnergy();
     recoPhotonInfo.scPreshowerEnergy = photon->superCluster()->preshowerEnergy();
     recoPhotonInfo.scPhiWidth = photon->superCluster()->phiWidth();
     recoPhotonInfo.scEtaWidth = photon->superCluster()->etaWidth();
     recoPhotonInfo.scNumBasicClusters = photon->superCluster()->clustersSize();


     // FIXME: swiss cross info too!

     // FIXME: ecal cell id of seed crystal?





  }// end of fill reco photon info

  

} //end of namespace

#endif
