#ifndef RECO_PHOTON_INFO_INC
#define RECO_PHOTON_INFO_INC

//********************************************************************
// Definition of a struct that can be used for storing reco photon info
// in a tree, from different analysers
// Also includes a Fill function to fill the struct from the appropriate objects
// and a string that can be used to define the tree branch
// 
//  $Id: RecoPhotonInfo.h,v 1.9 2010/11/22 14:44:27 chenders Exp $
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
//#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
//#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
//#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
//#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"


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

    double e2x2;
    double e4x4;
    double e2e9;

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

    // since we will now store also 'Fakeable objects' from data
    bool isFakeable;
  };


  // also include a string that can be used to define the tree branch
  // obviously this needs to be kept up-to-date with the struct definition
  // but now at least this only needs to be done here in this file, 
  // rather than in each individual analyser 
  std::string recoPhotonBranchDefString("pt/D:eta:phi:detEta:detPhi:detId/I:iEtaY/I:iPhiX/I:vx/D:vy:vz:r9:sigmaIetaIeta:sigmaEtaEta:maxEnergyXtal:e1x5:e2x5:e3x3:e5x5:r1x5:r2x5:swisscross:eMax:eLeft:eRight:eTop:eBottom:eSecond:e2x2:e4x4:e2e9:severityLevel/I:recHitFlag/I:maxRecHitTime/D:hadOverEm:hadDepth1OverEm:hadDepth2OverEm:hcalIso04/f:hcalIso03/f:ecalIso04:ecalIso03:trkIsoSumPtHollow04:trkIsoSumPtSolid04:trkIsoNtrksHollow04/I:trkIsoNtrksSolid04/I:trkIsoSumPtHollow03/f:trkIsoSumPtSolid03/f:trkIsoNtrksHollow03/I:trkIsoNtrksSolid03/I:esRatio/f:scRawEnergy/D:scPreshowerEnergy:scPhiWidth:scEtaWidth:scNumBasicClusters/I:isEB/O:isEE:isEBEtaGap:isEBPhiGap:isEERingGap:isEEDeeGap:isEBEEGap:hasPixelSeed:isFakeable");


  // useful function for ESratio
  // was originally in BkgAnalyser - now move to here
  double getESRatio(const reco::Photon *photon, const edm::Event& e, const edm::EventSetup& iSetup) {

    //get Geometry
    edm::ESHandle<CaloGeometry> caloGeometry;
    iSetup.get<CaloGeometryRecord>().get(caloGeometry);
    const CaloSubdetectorGeometry *geometry = caloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
    const CaloSubdetectorGeometry *& geometry_p = geometry;

    // Get ES rechits
    edm::Handle<EcalRecHitCollection> PreshowerRecHits;
    e.getByLabel(edm::InputTag("ecalPreshowerRecHit","EcalRecHitsES"), PreshowerRecHits);
    if( PreshowerRecHits.isValid() ) EcalRecHitCollection preshowerHits(*PreshowerRecHits);

    Float_t esratio=1.;

    if (fabs(photon->eta())>1.62) {

      const reco::CaloClusterPtr seed = (*photon).superCluster()->seed();    
      reco::CaloCluster cluster = (*seed);
      const GlobalPoint phopoint(cluster.x(), cluster.y(), cluster.z());
  
      DetId photmp1 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(phopoint, 1);
      DetId photmp2 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(phopoint, 2);
      ESDetId esfid = (photmp1 == DetId(0)) ? ESDetId(0) : ESDetId(photmp1);
      ESDetId esrid = (photmp2 == DetId(0)) ? ESDetId(0) : ESDetId(photmp2);

      int gs_esfid = -99;
      int gs_esrid = -99;
      gs_esfid = esfid.six()*32+esfid.strip();
      gs_esrid = esrid.siy()*32+esrid.strip();

      float esfe3 = 0.; 
      float esfe21 = 0.; 
      float esre3 = 0.; 
      float esre21 = 0.;

      const ESRecHitCollection *ESRH = PreshowerRecHits.product();
      EcalRecHitCollection::const_iterator esrh_it;
      for ( esrh_it = ESRH->begin(); esrh_it != ESRH->end(); esrh_it++) {
	ESDetId esdetid = ESDetId(esrh_it->id());
	if ( esdetid.plane()==1 ) {
	  if ( esdetid.zside() == esfid.zside() &&
	       esdetid.siy() == esfid.siy() ) {
	    int gs_esid = esdetid.six()*32+esdetid.strip();
	    int ss = gs_esid-gs_esfid;
	    if ( TMath::Abs(ss)<=10) {
	      esfe21 += esrh_it->energy();
	    } 
	    if ( TMath::Abs(ss)<=1) {
	      esfe3 += esrh_it->energy();
	    } 
	  }
	}
	if (esdetid.plane()==2 ){
	  if ( esdetid.zside() == esrid.zside() &&
	       esdetid.six() == esrid.six() ) {
	    int gs_esid = esdetid.siy()*32+esdetid.strip();
	    int ss = gs_esid-gs_esrid;
	    if ( TMath::Abs(ss)<=10) {
	      esre21 += esrh_it->energy();
	    } 
	    if ( TMath::Abs(ss)<=1) {
	      esre3 += esrh_it->energy();
	    } 
	  }
	}
      }
  
      if( (esfe21+esre21) == 0.) {
	esratio = 1.;
      }else{
	esratio = (esfe3+esre3) / (esfe21+esre21);
      }
    
    }
    return esratio;
    
  }

  float recHitE( const DetId id, const EcalRecHitCollection &recHits )
  {
    if ( id == DetId(0) ) {
      return 0;
    } else {
      EcalRecHitCollection::const_iterator it = recHits.find( id );
      if ( it != recHits.end() ) return (*it).energy();
    }
    return 0;
  }

  float recHitE( const DetId id, const EcalRecHitCollection & recHits,
					int di, int dj )
  {
    // in the barrel:   di = dEta   dj = dPhi
    // in the endcap:   di = dX     dj = dY
    
    DetId nid;
    if( id.subdetId() == EcalBarrel) nid = EBDetId::offsetBy( id, di, dj );
    else if( id.subdetId() == EcalEndcap) nid = EEDetId::offsetBy( id, di, dj );

    return ( nid == DetId(0) ? 0 : recHitE( nid, recHits ) );
  }

  float recHitApproxEt( const DetId id, const EcalRecHitCollection &recHits )
  {
    // for the time being works only for the barrel
    if ( id.subdetId() == EcalBarrel ) {
      return recHitE( id, recHits ) / cosh( EBDetId::approxEta( id ) );
    }
    return 0;
  }



  float E2overE9( const DetId id, const EcalRecHitCollection & recHits, float recHitEtThreshold = 10.0 , float recHitEtThreshold2 = 1.0 , bool avoidIeta85=false, bool KillSecondHit=true)
  {
    
    // compute e2overe9
    //  
    //   | | | |
    //   +-+-+-+
    //   | |1|2|
    //   +-+-+-+
    //   | | | |
    //
    //   1 - input hit,  2 - highest energy hit in a 3x3 around 1
    // 
    //   rechit 1 must have E_t > recHitEtThreshold
    //   rechit 2 must have E_t > recHitEtThreshold2
    //
    //   function returns value of E2/E9 centered around 1 (E2=energy of hits 1+2) if energy of 1>2
    //
    //   if energy of 2>1 and KillSecondHit is set to true, function returns value of E2/E9 centered around 2
    //   *provided* that 1 is the highest energy hit in a 3x3 centered around 2, otherwise, function returns 0
    
    
    if ( id.subdetId() == EcalBarrel ) {
      
      EBDetId ebId( id );
      
      // avoid recHits at |eta|=85 where one side of the neighbours is missing
      if ( abs(ebId.ieta())==85 && avoidIeta85) return 0;
      
      // select recHits with Et above recHitEtThreshold      
      float e1 = recHitE( id, recHits );
      float ete1=recHitApproxEt( id, recHits );     
      
      // check that rechit E_t is above threshold      
      if (ete1 < std::min(recHitEtThreshold,recHitEtThreshold2) ) return 0;      
      if (ete1 < recHitEtThreshold && !KillSecondHit ) return 0;		

      float e2=-1;
      float ete2=0;
      float s9 = 0;
      
      // coordinates of 2nd hit relative to central hit
      int e2eta=0;
      int e2phi=0;
      
      // LOOP OVER 3x3 ARRAY CENTERED AROUND HIT 1
      
      for ( int deta = -1; deta <= +1; ++deta ) {
	for ( int dphi = -1; dphi <= +1; ++dphi ) {
	  
	  // compute 3x3 energy	  
	  float etmp=recHitE( id, recHits, deta, dphi );
	  s9 += etmp;
	  
	  EBDetId idtmp=EBDetId::offsetBy(id,deta,dphi);
	  float eapproxet=recHitApproxEt( idtmp, recHits );
	  
	  // remember 2nd highest energy deposit (above threshold) in 3x3 array 
	  if (etmp>e2 && eapproxet>recHitEtThreshold2 && !(deta==0 && dphi==0)) {
	    e2=etmp;
	    ete2=eapproxet;
	    e2eta=deta;
	    e2phi=dphi;	    
	  }
	}
      }
      
      if ( e1 == 0 )  return 0;
      
      // return 0 if 2nd hit is below threshold
      if ( e2 == -1 ) return 0;
      
      // compute e2/e9 centered around 1st hit      
      float e2nd=e1+e2;
      float e2e9=0;
      
      if (s9!=0) e2e9=e2nd/s9;
      
      // if central hit has higher energy than 2nd hit
      //  return e2/e9 if 1st hit is above E_t threshold      
      if (e1 > e2 && ete1>recHitEtThreshold) return e2e9;
      
      // if second hit has higher energy than 1st hit      
      if ( e2 > e1 ) { 
	
	// return 0 if user does not want to flag 2nd hit, or
	// hits are below E_t thresholds - note here we
	// now assume the 2nd hit to be the leading hit.
	
	if (!KillSecondHit || ete2<recHitEtThreshold || ete1<recHitEtThreshold2) {
	  return 0;
	}	
	else {
	  
	  // LOOP OVER 3x3 ARRAY CENTERED AROUND HIT 2
	  
	  float s92nd=0;
          
	  float e2nd_prime=0;
	  int e2prime_eta=0;
	  int e2prime_phi=0;
	  
	  EBDetId secondid=EBDetId::offsetBy(id,e2eta,e2phi);
	  
	  for ( int deta = -1; deta <= +1; ++deta ) {
	    for ( int dphi = -1; dphi <= +1; ++dphi ) {
	      
	      // compute 3x3 energy
	      float etmp=recHitE( secondid, recHits, deta, dphi );
	      s92nd += etmp;
	      
	      if (etmp>e2nd_prime && !(deta==0 && dphi==0)) {
		e2nd_prime=etmp;
		e2prime_eta=deta;
		e2prime_phi=dphi;
	      }	      
	    }
	  }
	  
	  // if highest energy hit around E2 is not the same as the input hit, return 0;
	  
	  if (!(e2prime_eta==-e2eta && e2prime_phi==-e2phi)) 
	    { 
	      return 0;
	    }
	  
	  // compute E2/E9 around second hit 
	  float e2e9_2=0;
	  if (s92nd!=0) e2e9_2=e2nd/s92nd;
          
	  //   return the value of E2/E9 calculated around 2nd hit          
	  return e2e9_2;
	  
	}
      }
    } else if ( id.subdetId() == EcalEndcap ) {
      // only used for EB at the moment
      return 0;
    }
    return 0;
  }
  

  // also want a Fill function, that can fill the struct values from the appropriate objects
  // again, so that all editing only needs to be done here in this file

  // now trying with LazyTools etc...
  // lazyTools_ is a std::auto_ptr<EcalClusterLazyTools> in our analysers 
  // but you cannot declare std::auto_ptr in this fucntion arguments
  // because then calling the function transfers ownership of the auto_ptr
  // which then gets deleted when the function exits
  // and seg vios result
  // My solution: make the func prototype just EcalClusterLazyTools* 
  // and in the analyser, pass just the underlying pointer via lazyTools.get()
  // This seems to work, miraculously enough

  void FillRecoPhotonInfo(recoPhotonInfo_t &recoPhotonInfo, const reco::Photon *photon, EcalClusterLazyTools* lazyTools_,const EcalRecHitCollection * recHitsEB, const EcalRecHitCollection * recHitsEE, const EcalChannelStatus *ch_status,const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    

    recoPhotonInfo.pt = photon->et();
    recoPhotonInfo.eta = photon->eta();
    recoPhotonInfo.phi = photon->phi();
    
    recoPhotonInfo.detEta = photon->caloPosition().eta();
    recoPhotonInfo.detPhi = photon->caloPosition().phi();
    
    //   detId and crystal eta/phi

    reco::SuperClusterRef sc1 = photon->superCluster();
    std::pair<DetId,float> maxc1 = lazyTools_->getMaximum(*sc1);

    recoPhotonInfo.detId = maxc1.first.rawId();

    recoPhotonInfo.iEtaY = -99999;
    recoPhotonInfo.iPhiX = -99999;

     if(maxc1.first.subdetId() == EcalBarrel) {
       EBDetId ebId( maxc1.first );
       recoPhotonInfo.iEtaY = ebId.ieta(); // iEta in EB
       recoPhotonInfo.iPhiX = ebId.iphi(); // iPhi in EB
     }
     else if (maxc1.first.subdetId() == EcalEndcap) {
       EEDetId eeId( maxc1.first );
       recoPhotonInfo.iEtaY = eeId.iy(); // iY in EE
       recoPhotonInfo.iPhiX = eeId.ix(); // iX in EE       
     }


     // swiss cross and other spike-related info
     // we also need EB and EE rechits for some of this

     if (maxc1.first.subdetId() == EcalBarrel) {
       recoPhotonInfo.swisscross = EcalSeverityLevelAlgo::swissCross(maxc1.first, (*recHitsEB), 0);
     }
     else {
       recoPhotonInfo.swisscross = -999.99;
     }
     //     cout << "(Internal) Swiss Cross = "<< recoPhotonInfo.swisscross <<endl;
     
     recoPhotonInfo.eMax = lazyTools_->eMax(*sc1);
     recoPhotonInfo.eLeft = lazyTools_->eLeft(*sc1);
     recoPhotonInfo.eRight = lazyTools_->eRight(*sc1);
     recoPhotonInfo.eTop = lazyTools_->eTop(*sc1);
     recoPhotonInfo.eBottom = lazyTools_->eBottom(*sc1);
     recoPhotonInfo.eSecond = lazyTools_->e2nd(*sc1);

     recoPhotonInfo.e2x2 = lazyTools_->e2x2(*sc1);
     recoPhotonInfo.e4x4 = lazyTools_->e4x4(*sc1);

     recoPhotonInfo.e2e9 = E2overE9(maxc1.first, (*recHitsEB));

     //     cout << "(Internal) eMax = "<< recoPhotonInfo.eMax <<endl;
     //     cout << "(Internal) eSecond = "<< recoPhotonInfo.eSecond <<endl;
     
     // official ecal rec hit severity level & rec hit flag
     // also, time of highest energy rec hit

     const reco::CaloClusterPtr  seed = sc1->seed();

     DetId id = lazyTools_->getMaximum(*seed).first; 
     float time  = -999., outOfTimeChi2 = -999., chi2 = -999.;
     int   flags=-1, severity = -1; 

     const EcalRecHitCollection & rechits = ( photon->isEB() ? *recHitsEB : *recHitsEE); 
     EcalRecHitCollection::const_iterator it = rechits.find( id );
     if( it != rechits.end() ) { 
       time = it->time(); 
       outOfTimeChi2 = it->outOfTimeChi2();
       chi2 = it->chi2();
       flags = it->recoFlag();
       severity = EcalSeverityLevelAlgo::severityLevel( id, rechits, (*ch_status) );
     }

     
     EcalChannelStatusMap::const_iterator chit;
     chit = ch_status->getMap().find(id.rawId());
     int mystatus = -99;
     if( chit != ch_status->getMap().end() ){
       EcalChannelStatusCode ch_code = (*chit);
       mystatus = ch_code.getStatusCode();
     }

     recoPhotonInfo.severityLevel = severity;
     recoPhotonInfo.recHitFlag = flags;
     recoPhotonInfo.maxRecHitTime = time;

     //     cout << "(Internal) severity = " <<      recoPhotonInfo.severityLevel<<endl;
     //     cout << "(Internal) rechit flag = " << recoPhotonInfo.recHitFlag <<endl;
     //     cout << "(Internal)  rechit time = " << recoPhotonInfo.maxRecHitTime <<endl;

     //ES ratio - use the helper function
     recoPhotonInfo.esRatio = getESRatio(photon, iEvent, iSetup);


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


     // by default I'm going to set the new isFakeable to false here
     // if it is to be true, it needs to be filled inside hte analyser
     recoPhotonInfo.isFakeable = false;


  }// end of fill reco photon info

  // new compare function for sorting reco photon vectors
  bool comparePhotonsByPt(const reco::Photon &photon1, const reco::Photon &photon2) {

    // sorts such that highest pt photon first
    return(photon1.pt()>=photon2.pt());
  }

  // and also for std::pairs of reco photon with an int for tight/fakeable status
  bool comparePhotonPairsByPt(const std::pair<reco::Photon, bool> &pair1, const std::pair<reco::Photon, bool> &pair2) {

    // sort so that highest pt photon first
    return(pair1.first.pt()>=pair2.first.pt());
  }

  

} //end of namespace

#endif
