// -*- C++ -*-
//
// Package:    ADDrootTreeMaker
// Class:      ADDrootTreeMaker
// 
/**\class ADDrootTreeMaker ADDrootTreeMaker.cc MyAnaSpace/ADDrootTreeMaker/src/ADDrootTreeMaker.cc
   
Description: <one line class summary>

Implementation:
     <Notes on implementation>
*/
//
// Original Author:  "Duong Nguyen"
//         Created:  Thu Oct  29 11:50:24 CDT 2009
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "MyAnaSpace/ADDrootTreeMaker/interface/ADDrootTreeMaker.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h" 
#include "DataFormats/JetReco/interface/GenJet.h"

#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

// root include files
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include <iostream>
#include "../interface/ADDrootTreeMaker.h"

#include <sstream>

//
// constructors and destructor
//i
using namespace std ;
const bool isMC=true ;
const bool isReco = false ;
ADDrootTreeMaker::ADDrootTreeMaker(const edm::ParameterSet& iConfig)
{
  //==============================================
	// get input parameter
	//==============================================
	
	nPhotons_ = 0 ;
	nEvent_ = 0 ;
	
	nPhotonEvent_ = 0 ;
	nIsoPho40_ = 0 ;
	nPho25_ = 0 ;
	nDoublePho20_ = 0 ;
	
	fOutName_ = iConfig.getParameter<std::string>("outFileName") ;
	cout << "\n OutFileName: " << fOutName_ ;
	////////////////////////////////////////
  // initialized the TTree
  ////////////////////////////////////////

  rootfile_ = new TFile(fOutName_.c_str(),"RECREATE") ;
  tree_ = new TTree("ADDdiphotonTree","");

  tree_->Branch("evt_run",&evt_run_,"evt_run/I");
  tree_->Branch("evt_event",&evt_event_,"evt_event/I");

	tree_->Branch("nJets", &nJets_, "nJets/I") ;
	tree_->Branch("eta_recJet", eta_recJet_, "eta_recJet[nJets]/F") ;
	tree_->Branch("phi_recJet", phi_recJet_, "phi_recJet[nJets]/F") ;
	tree_->Branch("e_recJet", e_recJet_, "e_recJet[nJets]/F") ;
	tree_->Branch("et_recJet", et_recJet_, "et_recJet[nJets]/F") ;
	tree_->Branch("p_recJet", p_recJet_, "p_recJet[nJets]/F") ;
	tree_->Branch("emf_recJet", emf_recJet_, "emf_recJet[nJets]/F") ;
	tree_->Branch("genE_recJet", genE_recJet_, "genE_recJet[nJets]/F") ;
	tree_->Branch("genP_recJet", genP_recJet_, "genP_recJet[nJets]/F") ;
	tree_->Branch("genEt_recJet", genEt_recJet_, "genEt_recJet[nJets]/F") ;
	tree_->Branch("genEta_recJet", genEta_recJet_, "genEta_recJet[nJets]/F") ;
	tree_->Branch("genPhi_recJet", genPhi_recJet_, "genPhi_recJet[nJets]/F") ;
	
	tree_->Branch("nPhos", &nPhos_, "nPhos/I") ;
	tree_->Branch("eta_recPho", eta_recPho_, "eta_recPho[nPhos]/F") ;
	tree_->Branch("phi_recPho", phi_recPho_, "phi_recPho[nPhos]/F") ;
	tree_->Branch("e_recPho", e_recPho_, "e_recPho[nPhos]/F") ;
	tree_->Branch("pt_recPho", pt_recPho_, "pt_recPho[nPhos]/F") ;
	tree_->Branch("p_recPho", p_recPho_, "p_recPho[nPhos]/F") ;
	tree_->Branch("e5x5_recPho", e5x5_recPho_, "e5x5_recPho[nPhos]/F") ;
	tree_->Branch("e3x3_recPho", e3x3_recPho_, "e3x3_recPho[nPhos]/F") ;
	tree_->Branch("genE_recPho", genE_recPho_, "genE_recPho[nPhos]/F") ;
	tree_->Branch("genPt_recPho", genPt_recPho_, "genPt_recPho[nPhos]/F") ;
	tree_->Branch("genEta_recPho", genEta_recPho_, "genEta_recPho[nPhos]/F") ;
	tree_->Branch("genPhi_recPho", genPhi_recPho_, "genPhi_recPho[nPhos]/F") ;	
	tree_->Branch("genPdgId_recPho", genPdgId_recPho_, "genPdgId_recPho[nPhos]/I") ;
	tree_->Branch("motherGenPdgId_recPho", motherGenPdgId_recPho_, "motherGenPdgId_recPho[nPhos]/I") ;
	
	tree_->Branch("nRecoPhos", &nRecoPhos_, "nRecoPhos/I") ;
	tree_->Branch("eta_recoPho", eta_recoPho_, "eta_recoPho[nRecoPhos]/F") ;
	tree_->Branch("phi_recoPho", phi_recoPho_, "phi_recoPho[nRecoPhos]/F") ;
	tree_->Branch("e_recoPho", e_recoPho_, "e_recoPho[nRecoPhos]/F") ;
	tree_->Branch("pt_recoPho", pt_recoPho_, "pt_recoPho[nRecoPhos]/F") ;
	
	tree_->Branch("nGenPho", &nGenPho_, "nGenPho/I") ;
	tree_->Branch("genE_pho", genE_pho_, "genE_pho[nGenPho]/F") ;
	tree_->Branch("genPt_pho", genPt_pho_, "genPt_pho[nGenPho]/F") ;
	tree_->Branch("genEta_pho", genEta_pho_, "genEta_pho[nGenPho]/F") ;
	tree_->Branch("genPhi_pho", genPhi_pho_, "genPhi_pho[nGenPho]/F") ;	
	tree_->Branch("genPdgId_pho", genPdgId_pho_, "genPdgId_pho[nGenPho]/I") ;
	tree_->Branch("motherGenPdgId_pho", motherGenPdgId_pho_, "motherGenPdgId_pho[nGenPho]/I") ;
	tree_->Branch("motherGenStatus_pho", motherGenStatus_pho_, "motherGenStatus_pho[nGenPho]/I") ;

	tree_->Branch("isAlsoElectron", isAlsoElectron_, "isAlsoElectron[nPhos]/I") ;
	tree_->Branch("doPassDefaultPhotonIDcut", doPassDefaultPhotonIDcut_, "doPassDefaultPhotonIDcut[nPhos]/I") ;
	tree_->Branch("isolationEcalRecHit", isolationEcalRecHit_, "isolationEcalRecHit[nPhos]/F") ;
	tree_->Branch("isolationHcalRecHit", isolationHcalRecHit_, "isolationHcalRecHit[nPhos]/F") ;
	tree_->Branch("isolationSolidTrkCone", isolationSolidTrkCone_, "isolationSolidTrkCone[nPhos]/F") ;
	tree_->Branch("isolationHollowTrkCone", isolationHollowTrkCone_, "isolationHollowTrkCone[nPhos]/F") ;
	tree_->Branch("nTrkSolidCone", nTrkSolidCone_, "nTrkSolidCone[nPhos]/I") ;
	tree_->Branch("nTrkHollowCone", nTrkHollowCone_, "nTrkHollowCone[nPhos]/I") ;
	tree_->Branch("r9", r9_, "r9[nPhos]/F") ;
	tree_->Branch("hadronicOverEm", hadronicOverEm_, "hadronicOverEm[nPhos]/F") ;
	tree_->Branch("hasPixelSeed", hasPixelSeed_, "hasPixelSeed[nPhos]/I") ;
	
	tree_->Branch("isolationEcalRecHit_recoPho", isolationEcalRecHit_recoPho_, "isolationEcalRecHit_recoPho[nRecoPhos]/F") ;
	tree_->Branch("isolationHcalRecHit_recoPho", isolationHcalRecHit_recoPho_, "isolationHcalRecHit_recoPho[nRecoPhos]/F") ;
	tree_->Branch("isolationSolidTrkCone_recoPho", isolationSolidTrkCone_recoPho_, "isolationSolidTrkCone_recoPho[nRecoPhos]/F") ;
	tree_->Branch("isolationHollowTrkCone_recoPho", isolationHollowTrkCone_recoPho_, "isolationHollowTrkCone_recoPho[nRecoPhos]/F") ;
	tree_->Branch("nTrkSolidCone_recoPho", nTrkSolidCone_recoPho_, "nTrkSolidCone_recoPho[nRecoPhos]/I") ;
	tree_->Branch("nTrkHollowCone_recoPho", nTrkHollowCone_recoPho_, "nTrkHollowCone_recoPho[nRecoPhos]/I") ;
	tree_->Branch("r9_recoPho", r9_recoPho_, "r9_recoPho[nRecoPhos]/F") ;
	tree_->Branch("hadronicOverEm_recoPho", hadronicOverEm_recoPho_, "hadronicOverEm_recoPho[nRecoPhos]/F") ;
	tree_->Branch("hasPixelSeed_recoPho", hasPixelSeed_recoPho_, "hasPixelSeed_recoPho[nRecoPhos]/I") ;
	
	tree_->Branch("nPhotonPair", &nPhotonPair_, "nPhotonPair/I") ;
	tree_->Branch("indPho1_phoPair", indPho1_phoPair_, "indPho1_phoPair[nPhotonPair]/I") ;
	tree_->Branch("indPho2_phoPair", indPho2_phoPair_, "indPho2_phoPair[nPhotonPair]/I") ;
	
	tree_->Branch("invMass_phoPair", invMass_phoPair_, "invMass_phoPair[nPhotonPair]/F") ;
	tree_->Branch("genInvMass_phoPair", genInvMass_phoPair_, "genInvMass_phoPair[nPhotonPair]/F") ;
	tree_->Branch("dR_phoPair", dR_phoPair_, "dR_phoPair[nPhotonPair]/F") ;
	tree_->Branch("dPhi_phoPair", dPhi_phoPair_, "dPhi_phoPair[nPhotonPair]/F") ;
	
	tree_->Branch("nRecoPhotonPair", &nRecoPhotonPair_, "nRecoPhotonPair/I") ;
	tree_->Branch("indPho1_phoPair_recoPho", indPho1_phoPair_recoPho_, "indPho1_phoPair_recoPho[nRecoPhotonPair]/I") ;
	tree_->Branch("indPho2_phoPair_recoPho", indPho2_phoPair_recoPho_, "indPho2_phoPair_recoPho[nRecoPhotonPair]/I") ;
	
	tree_->Branch("invMass_phoPair_recoPho", invMass_phoPair_recoPho_, "invMass_phoPair_recoPho[nRecoPhotonPair]/F") ;
	tree_->Branch("dR_phoPair_recoPho", dR_phoPair_recoPho_, "dR_phoPair_recoPho[nRecoPhotonPair]/F") ;
	tree_->Branch("dPhi_phoPair_recoPho", dPhi_phoPair_recoPho_, "dPhi_phoPair_recoPho[nRecoPhotonPair]/F") ;
	
	tree_->Branch("nPhoJetPair", &nPhoJetPair_, "nPhoJetPair/I") ;
	tree_->Branch("minDrJet_phoJetPair", minDrJet_phoJetPair_, "minDrJet_phoJetPair[nPhoJetPair]/F") ;
	tree_->Branch("indJet_phoJetPair", indJet_phoJetPair_, "indJet_phoJetPair[nPhoJetPair]/I") ;
	tree_->Branch("minDrPho_phoJetPair", minDrPho_phoJetPair_, "minDrPho_phoJetPair[nPhoJetPair]/F") ;
	tree_->Branch("indPho_phoJetPair", indPho_phoJetPair_, "indPho_phoJetPair[nPhoJetPair]/I") ;
	tree_->Branch("dR_phoJetPair", dR_phoJetPair_, "dR_phoJetPair[nPhoJetPair]/F") ;
	tree_->Branch("dPhi_phoJetPair", dPhi_phoJetPair_, "dPhi_phoJetPair[nPhoJetPair]/F") ;
	tree_->Branch("invMass_phoJetPair", invMass_phoJetPair_, "invMass_phoJetPair[nPhoJetPair]/F") ;
	
	tree_->Branch("nMCpart", &nMCpart_, "nMCpart/I") ;
	tree_->Branch("eta_genPart", eta_genPart_, "eta_genPart[nMCpart]/F") ;
	tree_->Branch("phi_genPart", phi_genPart_, "phi_genPart[nMCpart]/F") ;
	tree_->Branch("pt_genPart", pt_genPart_, "pt_genPart[nMCpart]/F") ;
	tree_->Branch("e_genPart", e_genPart_, "e_genPart[nMCpart]/F") ;
	tree_->Branch("pdgId_genPart", pdgId_genPart_, "pdgId_genPart[nMCpart]/I") ;
	tree_->Branch("status_genPart", status_genPart_, "status_genPart[nMCpart]/I") ;
	tree_->Branch("motherPdgId_genPart", motherPdgId_genPart_, "motherPdgId_genPart[nMCpart]/I") ;
	tree_->Branch("motherStatus_genPart", motherStatus_genPart_, "motherStatus_genPart[nMCpart]/I") ;
	tree_->Branch("grandMotherPdgId_genPart", grandMotherPdgId_genPart_, "grandMotherPdgId_genPart[nMCpart]/I") ;
	tree_->Branch("grandMotherStatus_genPart", grandMotherStatus_genPart_, "grandMotherStatus_genPart[nMCpart]/I") ;
	
	tree_->Branch("nGenJet", &nGenJet_, "nGenJet/I") ;
	tree_->Branch("eta_genJet", eta_genJet_, "eta_genJet[nGenJet]/F") ;
	tree_->Branch("phi_genJet", phi_genJet_, "phi_genJet[nGenJet]/F") ;
	tree_->Branch("et_genJet", et_genJet_, "et_genJet[nGenJet]/F") ;
	tree_->Branch("e_genJet", e_genJet_, "e_genJet[nGenJet]/F") ;

	
	tree_->Branch("pThat", &pThat_, "pThat/F") ;
	tree_->Branch("processId", &processId_, "processId/I") ;
	
	tree_->Branch("HLT_IsoPhoton40_L1R", &HLT_IsoPhoton40_L1R_, "HLT_IsoPhoton40_L1R/I") ;
  tree_->Branch("HLT_Photon25_L1R", &HLT_Photon25_L1R_, "HLT_Photon25_L1/I") ;
  tree_->Branch("HLT_DoubleIsoPhoton20_L1R", &HLT_DoubleIsoPhoton20_L1R_, "HLT_DoubleIsoPhoton20_L1R/I") ;

	tree_->Branch("nMCpart3", &nMCpart3_, "nMCpart3/I") ;
	tree_->Branch("genEta_part3", genEta_part3_, "genEta_part3[nMCpart3]/F") ;
	tree_->Branch("genPhi_part3", genPhi_part3_, "genPhi_part3[nMCpart3]/F") ;
	tree_->Branch("genPt_part3", genPt_part3_, "genPt_part3[nMCpart3]/F") ;
	tree_->Branch("pdgId_part3", pdgId_part3_, "pdgId_part3[nMCpart3]/I") ;
	tree_->Branch("nPi0", &nPi0_, "nPi0/I") ;
	tree_->Branch("eta_pi0", eta_pi0_, "eta_pi0[nPi0]/F") ;
	tree_->Branch("phi_pi0", phi_pi0_, "phi_pi0[nPi0]/F") ;
	tree_->Branch("pt_pi0", pt_pi0_, "pt_pi0[nPi0]/F") ;
	
	hNevent_ = new TH1F("Nevent", "", 2, 0, 2) ;
}

ADDrootTreeMaker::~ADDrootTreeMaker()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

  rootfile_->cd();
  tree_->Write();
  hNevent_->Write() ;
	rootfile_->Close();
  
}


//
// member functions
//

// ------------ method called to for each event  ------------
void ADDrootTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
	
	nEvent_ += 1 ;
	hNevent_->Fill(1.5) ;
	//////////////////////////////////////
	//Input branches
	//////////////////////////////////////
	
	edm::Handle< std::vector<pat::Photon> > myPhos ;
	iEvent.getByLabel("selectedLayer1Photons", myPhos) ;
//	iEvent.getByLabel("cleanLayer1Photons", myPhos) ;

	edm::Handle< std::vector<pat::Jet> > myJets ;
	iEvent.getByLabel("selectedLayer1Jets", myJets) ;
//	iEvent.getByLabel("cleanLayer1Jets", myJets) ;
	
	edm::Handle<double> pThat ;
	iEvent.getByLabel("genEventScale", pThat) ;

	edm::Handle<std::vector<reco::GenParticle> > myGenParts ;
	iEvent.getByLabel("genParticles", myGenParts) ;
	
	edm::Handle<std::vector<reco::GenJet> > myGenJets ;
	iEvent.getByLabel("iterativeCone5GenJets", myGenJets) ;
	
	edm::Handle<std::vector<reco::Photon> > myRecoPhos ;
	iEvent.getByLabel("photons", myRecoPhos) ;

	//////////////////////////////////////////
	//Fill process Id
	//////////////////////////////////////////
	Handle<edm::HepMCProduct> genEvt;
	iEvent.getByLabel("generator", genEvt);
	if(!(genEvt.isValid())) {
	  LogError("ADDdiphotonTreeAnalyzer") << "genParticles not found";
	  return;
	}

	/////////////////////////////////////////
	//Fill trigger information
	/////////////////////////////////////////
	
	edm::Handle<edm::TriggerResults>  hltresults ;
	edm::InputTag trigResult("TriggerResults::HLT") ;
	iEvent.getByLabel(trigResult, hltresults) ;

	edm::TriggerNames triggerNames ;

	HLT_IsoPhoton40_L1R_ = 0 ;
  HLT_Photon25_L1R_ = 0 ;
  HLT_DoubleIsoPhoton20_L1R_ = 0 ;

  if (hltresults.isValid()) {
    int ntrigs = hltresults->size();
    if (ntrigs==0){std::cout << "%HLTInfo -- No trigger name given in TriggerResults of the input " << std::endl;}

    triggerNames.init(*hltresults);

    for (int itrig = 0; itrig != ntrigs; ++itrig){

      string trigName=triggerNames.triggerName(itrig);
      bool accept = hltresults->accept(itrig);

      if (trigName == "HLT_IsoPhoton40_L1R") HLT_IsoPhoton40_L1R_ = accept ;
      if (trigName == "HLT_Photon25_L1R") HLT_Photon25_L1R_ = accept ;
      if (trigName == "HLT_DoubleIsoPhoton20_L1R") HLT_DoubleIsoPhoton20_L1R_ = accept ;
      //std::cout << "%HLTInfo --  Number of HLT Triggers: " << ntrigs << std::endl;
      //std::cout << "%HLTInfo --  HLTTrigger(" << itrig << "): " << trigName << " = " << accept << std::endl;
    } //end loop over trigger
  } //end if (hltresults.isValid())

	//////////////////////////////////////////
	//Fill signal process
	//////////////////////////////////////////
	//////////////////////////////////////////
	
	const HepMC::GenEvent* myGenEvent = genEvt->GetEvent();
	processId_ = myGenEvent->signal_process_id() ;
	
	/*
	if (nEvent_ <= 10) {
		cout << "\n ===============================" ;
		cout << "\n Event: " << nEvent_ ;
		myGenEvent->print() ;
	}
	*/
	/*
	if (processId_ >= 0) {
		//cout << "\n ============diphoton event=============" ;
		for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin();   p != myGenEvent->particles_end(); ++p ) {
			HepMC::GenVertex* parentVertex = (*p)->production_vertex() ;
			if (parentVertex !=0 ) {
			  if ((*p)->status() == 1 && (*p)->pdg_id() == 22 && parentVertex->particles_in_size() >= 1) {
				  //cout << "\n Number of mother: " << parentVertex->particles_in_size();
				  HepMC::GenVertex::particles_in_const_iterator parentItr = parentVertex->particles_in_const_begin() ;
				  if ((*parentItr)->pdg_id() == 21) {
						//cout << "\n Gamma: " << (*p)->pdg_id() << "  " << (*p)->status() ;
						//for ( ; parentItr != parentVertex->particles_in_const_end(); parentItr++) {
						//	cout << "\n Mother: " << (*parentItr)->pdg_id() << "  " << (*parentItr)->status() ;
						 cout << "\n Gamma has parent as gluon: " << std::endl ;
						 myGenEvent->print() ;
					}
				}
			}
		}
	}*/
	/*
	int nPromptPho = 0 ;
	for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin();   p != myGenEvent->particles_end(); ++p ) {
		HepMC::GenVertex* parentVertex = (*p)->production_vertex() ;
		if ((*p)->status() == 1 && (*p)->pdg_id() == 22 && parentVertex->particles_in_size() >= 1) {
			if (parentVertex !=0 ) {
			  HepMC::GenVertex::particles_in_const_iterator parentItr = parentVertex->particles_in_const_begin() ;
				if ((*parentItr)->pdg_id() == 22 && (*parentItr)->status() == 3) {
					if (parentVertex->particles_in_size() > 1) cout << "\n Number of mother: " << parentVertex->particles_in_size();
					nPromptPho++ ;
				}
			}
		}
	}
	if(nPromptPho >= 2) { cout << "\n diphoton event: " << nPromptPho ; }
	*/
	////////////////////////////////////////
	//Fill pThat
	///////////////////////////////////////
	pThat_ = (*pThat) ;
	////////////////////////////////////////
  // Fill the Event Info
  ////////////////////////////////////////
  
  evt_run_ = iEvent.id().run();
  evt_event_ = iEvent.id().event();
	//////////////////////////////////////////
	//Fill jet
	/////////////////////////////////////////
	nJets_ = 0 ;
	for (unsigned int i = 0 ; (i < (*myJets).size()) && (i < maxJets) ; ++i) {
		
		eta_recJet_[i] = (*myJets)[i].eta() ;
		phi_recJet_[i] = (*myJets)[i].phi() ;
		e_recJet_[i] = (*myJets)[i].energy() ;
		et_recJet_[i] = (*myJets)[i].et() ;
		p_recJet_[i] = (*myJets)[i].p() ;
		emf_recJet_[i] = (*myJets)[i].emEnergyFraction() ;
		
		const reco::GenJet* genJetPtr = (*myJets)[i].genJet() ;
		if (genJetPtr != 0) {
			genEta_recJet_[i] = genJetPtr->eta() ;
			genPhi_recJet_[i] = genJetPtr->phi() ;
			genEt_recJet_[i] = genJetPtr->et() ;
			genE_recJet_[i] = genJetPtr->energy() ;
			genP_recJet_[i] = genJetPtr->p() ;
		} //end if (genJetPtr)
		else {
			genEta_recJet_[i] = -100 ;
			genPhi_recJet_[i] = -100 ;
			genEt_recJet_[i] = -100 ;
			genE_recJet_[i] = -100 ;
			genP_recJet_[i] = -100 ;
		}
		
		nJets_++ ;
		
	}
	
	////////////////////////////////////////
	//Fill genJet
	////////////////////////////////////////
	nGenJet_ = 0 ;

	for (unsigned int i = 0 ; (i < (*myGenJets).size()) && (i < maxGenJet) ; ++i) {
		
		eta_genJet_[i] = (*myGenJets)[i].eta() ;
		phi_genJet_[i] = (*myGenJets)[i].phi() ;
		et_genJet_[i] = (*myGenJets)[i].et() ;
		e_genJet_[i] = (*myGenJets)[i].energy() ;

		nGenJet_++ ;
	
	}

	
	/////////////////////////////////////////
	//Fill photon
	////////////////////////////////////////
	nPhos_ = 0 ;
	
	bool doPassEtaPtCut = false ;
	
	for (unsigned int i = 0 ; (i < (*myPhos).size()) && (i < maxPhos) ; ++i) {
		eta_recPho_[i] = (*myPhos)[i].eta() ;
		phi_recPho_[i] = (*myPhos)[i].phi() ;
		pt_recPho_[i] = (*myPhos)[i].pt() ;
		e_recPho_[i] = (*myPhos)[i].energy() ;
		p_recPho_[i] = (*myPhos)[i].p() ;
		e5x5_recPho_[i] = (*myPhos)[i].e5x5() ;
		e3x3_recPho_[i] = (*myPhos)[i].e3x3() ;
		isAlsoElectron_[i] = (*myPhos)[i].isElectron() ;
		nTrkSolidCone_[i] = (*myPhos)[i].nTrkSolidConeDR04() ;
		nTrkHollowCone_[i] = (*myPhos)[i].nTrkHollowConeDR04() ;
//		isolationEcalRecHit_[i] = (*myPhos)[i].isolationEcalRecHit() ;
//		isolationHcalRecHit_[i] = (*myPhos)[i].isolationHcalRecHit() ;
//		isolationSolidTrkCone_[i] = (*myPhos)[i].isolationSolidTrkCone() ;
//		isolationHollowTrkCone_[i] = (*myPhos)[i].isolationHollowTrkCone() ;
		isolationEcalRecHit_[i] = (*myPhos)[i].ecalRecHitSumEtConeDR04() ;
		isolationHcalRecHit_[i] = (*myPhos)[i].hcalTowerSumEtConeDR04() ;
		isolationSolidTrkCone_[i] = (*myPhos)[i].trkSumPtSolidConeDR04() ;
		isolationHollowTrkCone_[i] = (*myPhos)[i].trkSumPtHollowConeDR04() ;
		r9_[i] = (*myPhos)[i].r9() ;
		hadronicOverEm_[i] = (*myPhos)[i].hadronicOverEm() ;
		hasPixelSeed_[i] = (*myPhos)[i].hasPixelSeed() ;
//		cout << "\n Pixel seed: " << hasPixelSeed_[i] ;	
		doPassDefaultPhotonIDcut_[i] = DoPassDefaultPhotonIDcut((*myPhos)[i]) ;
		
		if (fabs(eta_recPho_[i]) < 1.5 && pt_recPho_[i] > 200) doPassEtaPtCut = true ;
				
		const reco::Candidate* genPhoPtr = (*myPhos)[i].genParticle() ;
		if (genPhoPtr != 0) {
			genEta_recPho_[i] = genPhoPtr->eta() ;
			genPhi_recPho_[i] = genPhoPtr->phi() ;
			genPt_recPho_[i] = genPhoPtr->et() ;
			genE_recPho_[i] = genPhoPtr->energy() ;
			genP_recPho_[i] = genPhoPtr->p() ;
			genPdgId_recPho_[i] = genPhoPtr->pdgId() ;
			motherGenPdgId_recPho_[i] = genPhoPtr->mother(0)->pdgId() ;
		}
		else {
			genEta_recPho_[i] = -100 ;
			genPhi_recPho_[i] = -100 ;
			genPt_recPho_[i] = -100 ;
			genE_recPho_[i] = -100 ;
			genP_recPho_[i] = -100 ;
			genPdgId_recPho_[i] = -1000 ;
			motherGenPdgId_recPho_[i] = -1000 ;
		}
		
		
		nPhos_++ ;
	
	}
  
  /////////////////////////////////////////
  //fill reco::Photon
	/////////////////////////////////////////
	nRecoPhos_ = 0 ;
	for (unsigned int i = 0 ; (i < (*myRecoPhos).size()) && (i < maxPhos) ; ++i) {
		eta_recoPho_[i] = (*myRecoPhos)[i].eta() ;
		phi_recoPho_[i] = (*myRecoPhos)[i].phi() ;
		pt_recoPho_[i] = (*myRecoPhos)[i].pt() ;
		e_recoPho_[i] = (*myRecoPhos)[i].energy() ;
		nTrkSolidCone_recoPho_[i] = (*myRecoPhos)[i].nTrkSolidConeDR04() ;
		nTrkHollowCone_recoPho_[i] = (*myRecoPhos)[i].nTrkHollowConeDR04() ;
		isolationEcalRecHit_recoPho_[i] = (*myRecoPhos)[i].ecalRecHitSumEtConeDR04() ;
		isolationHcalRecHit_recoPho_[i] = (*myRecoPhos)[i].hcalTowerSumEtConeDR04() ;
		isolationSolidTrkCone_recoPho_[i] = (*myRecoPhos)[i].trkSumPtSolidConeDR04() ;
		isolationHollowTrkCone_recoPho_[i] = (*myRecoPhos)[i].trkSumPtHollowConeDR04() ;
		r9_recoPho_[i] = (*myRecoPhos)[i].r9() ;
		hadronicOverEm_recoPho_[i] = (*myRecoPhos)[i].hadronicOverEm() ;
		hasPixelSeed_recoPho_[i] = (*myRecoPhos)[i].hasPixelSeed() ;
		
		nRecoPhos_++ ;
	
	}
	/////////////////////////////////////////
	/////////////////////////////////////////
	//check trigger
	/////////////////////////////////////////
	if (doPassEtaPtCut) {
		nPhotonEvent_++ ;
		if (HLT_IsoPhoton40_L1R_) nIsoPho40_++ ;
		if (HLT_Photon25_L1R_) nPho25_++ ;
		if (HLT_DoubleIsoPhoton20_L1R_) nDoublePho20_++ ;
	}
	/////////////////////////////////////////
	//Fill gen particle info
	////////////////////////////////////////
	
	nMCpart_ = 0 ;
	nGenPho_ = 0 ;
	nMCpart3_ = 0 ;
	nPi0_ = 0 ;
	
	for (unsigned int i = 0 ; i < (*myGenParts).size() ; ++i) {
		//fill gen pho
		if ((*myGenParts)[i].pdgId() == 22 && (*myGenParts)[i].status() == 1) {
			
			genPdgId_pho_[nGenPho_] = (*myGenParts)[i].pdgId() ;
			genE_pho_[nGenPho_] = (*myGenParts)[i].energy() ;
			genPt_pho_[nGenPho_] = (*myGenParts)[i].et() ;
			genEta_pho_[nGenPho_] = (*myGenParts)[i].eta() ;
			genPhi_pho_[nGenPho_] = (*myGenParts)[i].phi() ;
			motherGenPdgId_pho_[nGenPho_] = (*myGenParts)[i].mother(0)->pdgId() ;
			motherGenStatus_pho_[nGenPho_] = (*myGenParts)[i].mother(0)->status() ;

			nGenPho_++ ;
		
		}
		//store particle status 3
		if ((*myGenParts)[i].status()==3) {
			bool hasDauSta3 = false ;
			size_t nDau = (*myGenParts)[i].numberOfDaughters() ;
			for (size_t iDau = 0 ; iDau < nDau ; ++iDau) {
				if ((*myGenParts)[i].daughter(iDau)->status()==3) hasDauSta3 = true ;
			}
			if (!hasDauSta3) {
				genEta_part3_[nMCpart3_] = (*myGenParts)[i].eta() ;
				genPhi_part3_[nMCpart3_] = (*myGenParts)[i].phi() ;
				genPt_part3_[nMCpart3_] = (*myGenParts)[i].et() ;
				pdgId_part3_[nMCpart3_] = (*myGenParts)[i].pdgId() ;
				
				nMCpart3_++ ;
				/*
				size_t nMo = (*myGenParts)[i].numberOfMothers() ;
				size_t nDau = (*myGenParts)[i].numberOfDaughters() ;
				cout << "\n Particle: " << (*myGenParts)[i].pdgId() << "  " << (*myGenParts)[i].eta() << "  " << (*myGenParts)[i].phi() << "  " << (*myGenParts)[i].et() ;
				cout << " | " ;
				for (size_t iMo = 0 ; iMo < nMo ; ++iMo) {
					cout << "    " << (*myGenParts)[i].mother(iMo)->pdgId() ;
				}
				cout << " | " ;
				for (size_t iDau = 0 ; iDau < nDau ; ++iDau) {
					cout << "    " << (*myGenParts)[i].daughter(iDau)->pdgId() ;
				}
				cout << "\n" ;
				*/
			}
		}
		
		//fill pi0
		if ((*myGenParts)[i].pdgId()==111) {
			eta_pi0_[nPi0_] = (*myGenParts)[i].eta() ;
			phi_pi0_[nPi0_] = (*myGenParts)[i].phi() ;
			pt_pi0_[nPi0_] = (*myGenParts)[i].et() ;
			nPi0_++ ;
		}
		
		/*
		//fill gen particle
		eta_genPart_[nMCpart_] = (*myGenParts)[i].eta() ;
		phi_genPart_[nMCpart_] = (*myGenParts)[i].phi() ;
		pt_genPart_[nMCpart_] = (*myGenParts)[i].et() ;
		e_genPart_[nMCpart_] = (*myGenParts)[i].energy() ;
		pdgId_genPart_[nMCpart_] = (*myGenParts)[i].pdgId() ;
		status_genPart_[nMCpart_] = (*myGenParts)[i].status() ;
		if ((*myGenParts)[i].mother(0) != 0) {
			motherPdgId_genPart_[nMCpart_] = (*myGenParts)[i].mother(0)->pdgId() ;
			motherStatus_genPart_[nMCpart_] = (*myGenParts)[i].mother(0)->status() ;
			if ((*myGenParts)[i].mother(0)->mother(0) != 0) {
				grandMotherPdgId_genPart_[nMCpart_] = (*myGenParts)[i].mother(0)->mother(0)->pdgId() ;
				grandMotherStatus_genPart_[nMCpart_] = (*myGenParts)[i].mother(0)->mother(0)->status() ;
			}
			else {
				grandMotherPdgId_genPart_[nMCpart_] = -30000 ;
				grandMotherStatus_genPart_[nMCpart_] = -1 ;
			}
		}
		else {
			motherPdgId_genPart_[nMCpart_] = -30000 ;
			motherStatus_genPart_[nMCpart_] = -1 ;
			grandMotherPdgId_genPart_[nMCpart_] = -1000 ;
			grandMotherStatus_genPart_[nMCpart_] = -1 ;
		}
		
		nMCpart_++ ;
		*/
	}
	
	/////////////////////////////////////////
	// Fill diphoton object
	/////////////////////////////////////////
	
	nPhotonPair_ = 0 ;
	for (unsigned int i = 0 ; i < (*myPhos).size() ; ++i) {
	  for (unsigned int j = i + 1 ; j < (*myPhos).size() ; ++j) { 
			if ((nPhotonPair_ <= maxPhotonPair) && (*myPhos)[i].pt() > 10 && (*myPhos)[j].pt() > 10) {// && fabs((*myPhos)[i].eta()) < 2.5 && fabs((*myPhos)[j].eta()) < 2.5) {
			indPho1_phoPair_[nPhotonPair_] = i ;
			indPho2_phoPair_[nPhotonPair_] = j ;
			
			invMass_phoPair_[nPhotonPair_] =  InvariantMass((*myPhos)[i].energy(), (*myPhos)[i].px(), (*myPhos)[i].py(), (*myPhos)[i].pz(), (*myPhos)[j].energy(), (*myPhos)[j].px(), (*myPhos)[j].py(), (*myPhos)[j].pz());
			dR_phoPair_[nPhotonPair_] = dR((*myPhos)[i].eta(), (*myPhos)[i].phi(), (*myPhos)[j].eta(), (*myPhos)[j].phi()) ;
			dPhi_phoPair_[nPhotonPair_] = dPhi((*myPhos)[i].phi(), (*myPhos)[j].phi()) ;
			const reco::Candidate* genPho1 = (*myPhos)[i].genParticle() ;
			const reco::Candidate* genPho2 = (*myPhos)[j].genParticle() ;
			if (genPho1 != 0 && genPho2 != 0) genInvMass_phoPair_[nPhotonPair_] = InvariantMass(genPho1->energy(), genPho1->px(), genPho1->py(), genPho1->pz(), genPho2->energy(), genPho2->px(), genPho2->py(), genPho2->pz()) ;
			else genInvMass_phoPair_[nPhotonPair_] = -1 ;

			nPhotonPair_ += 1 ;
			
			} 
		}
	}

////////////////////////////////////////
//Fill reco diphoton pair
////////////////////////////////////////
nRecoPhotonPair_ = 0 ;
for (unsigned int i = 0 ; i < (*myRecoPhos).size() ; ++i) {
	for (unsigned int j = i + 1  ; j < (*myRecoPhos).size() ; ++j) {	
			if ((nRecoPhotonPair_ <= maxPhotonPair) && (*myRecoPhos)[i].pt() > 10 && (*myRecoPhos)[j].pt() > 10) {// && fabs((*myRecoPhos)[i].eta()) < 2.5 && fabs((*myRecoPhos)[j].eta()) < 2.5) {
			indPho1_phoPair_recoPho_[nRecoPhotonPair_] = i ;
			indPho2_phoPair_recoPho_[nRecoPhotonPair_] = j ;
			
			invMass_phoPair_recoPho_[nRecoPhotonPair_] =  InvariantMass((*myRecoPhos)[i].energy(), (*myRecoPhos)[i].px(), (*myRecoPhos)[i].py(), (*myRecoPhos)[i].pz(), (*myRecoPhos)[j].energy(), (*myRecoPhos)[j].px(), (*myRecoPhos)[j].py(), (*myRecoPhos)[j].pz());
			dR_phoPair_recoPho_[nRecoPhotonPair_] = dR((*myRecoPhos)[i].eta(), (*myRecoPhos)[i].phi(), (*myRecoPhos)[j].eta(), (*myRecoPhos)[j].phi()) ;
			dPhi_phoPair_recoPho_ [nRecoPhotonPair_] = dPhi((*myRecoPhos)[i].phi(), (*myRecoPhos)[j].phi()) ;

			nRecoPhotonPair_ += 1 ;
			}
	}
}
	///////////////////////////////////////
	//Fill photonJet Pair
	///////////////////////////////////////

	nPhoJetPair_ = 0 ;
	for (unsigned int iPho = 0 ; iPho < (*myPhos).size() ; ++iPho) {
		for (unsigned int iJet = 0 ; iJet < (*myJets).size() ; ++iJet) {
			float dRphoJet = dR((*myPhos)[iPho].eta(), (*myPhos)[iPho].phi(), (*myJets)[iJet].eta(), (*myJets)[iJet].phi()) ;
			//cut on photon and jet pt and eta
			if (fabs((*myPhos)[iPho].eta()) < 2.4 && (*myPhos)[iPho].et() > 50 && fabs((*myJets)[iJet].eta()) < 2.4 && (*myJets)[iJet].et() > 50 && dRphoJet > 0.7) {
				
				minDrPho_phoJetPair_[nPhoJetPair_] = minDr((*myPhos)[iPho].eta(), (*myPhos)[iPho].phi(), (*myJets)) ;
				indPho_phoJetPair_[nPhoJetPair_] = iPho ;
					
				minDrJet_phoJetPair_[nPhoJetPair_] = minDr((*myJets)[iJet].eta(), (*myJets)[iJet].phi(), (*myPhos)) ;
				indJet_phoJetPair_[nPhoJetPair_] = iJet ;
				
			//	invMass_phoJetPair_[nPhoJetPair_] = InvariantMass((*myPhos)[iPho].energy(), (*myPhos)[iPho].px(), (*myPhos)[iPho].py(), (*myPhos)[iPho].pz(), (*myJets)[iJet].energy(), (*myJets)[iJet].px(), (*myJets)[iJet].py(), (*myJets)[iJet].pz()) ;
				invMass_phoJetPair_[nPhoJetPair_] = InvariantMass((*myPhos)[iPho].p4(), (*myJets)[iJet].p4()) ;
				dR_phoJetPair_[nPhoJetPair_] = dR((*myPhos)[iPho].eta(), (*myPhos)[iPho].phi(), (*myJets)[iJet].eta(), (*myJets)[iJet].phi()) ;
				dPhi_phoJetPair_[nPhoJetPair_] = dPhi((*myPhos)[iPho].phi(), (*myJets)[iJet].phi()) ; 

				nPhoJetPair_++ ;

			} //end cut on photon and jet pt and eta
		} //end loop over jet
	} //end loop over photons
	
	////////////////////////////////////////
	///filter and fill the tree
	///////////////////////////////////////
	//if (nPhotonPair_ > 0 || nPhoJetPair_ > 0) {
		rootfile_->cd();
		tree_->Fill();
		hNevent_->Fill(0.5) ;
	//}
}

float ADDrootTreeMaker::InvariantMass(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& v1, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& v2) {
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > v = v1 + v2 ;
  return v.mass() ;
}

float ADDrootTreeMaker::InvariantMass(double energy1, double px1, double py1, double pz1, double energy2, double px2, double py2, double pz2) {
  //float mass2 = pow((energy1 + energy2),2) - pow((px1+px2),2) - pow((py1+py2),2) - pow((pz1+pz2),2) ;
	//if (mass2 < 0) return -sqrt(-mass2) ;
	float mass2 = 2*energy1*energy2*(1 - (px1*px2+py1*py2+pz1*pz2)/(energy1*energy2)) ;
	return sqrt(mass2) ;
}

float ADDrootTreeMaker::dPhi(double phi1, double phi2) {
  float dPhiTmp = phi1 - phi2 ;
 	if (dPhiTmp <= -TMath::Pi()) return (dPhiTmp + 2.0*TMath::Pi()) ;
  if (dPhiTmp > TMath::Pi()) return (2.0*TMath::Pi()-dPhiTmp) ;
  return fabs(dPhiTmp) ;
}

float ADDrootTreeMaker::dR(double eta1, double phi1, double eta2, double phi2) {
  return sqrt(pow(eta1-eta2, 2) + pow(dPhi(phi1,phi2), 2)) ;
}

float ADDrootTreeMaker::minDr(double eta, double phi, const std::vector<pat::Jet>& jets) {
	float minDrVal = 100 ;
	for (int i = 0 ; i < jets.size() ; ++i)
		if (jets[i].et() > 40 && dR(eta, phi, jets[i].eta(), jets[i].phi()) < minDrVal) minDrVal = dR(eta, phi, jets[i].eta(), jets[i].phi()) ;
	return minDrVal ;
}

float ADDrootTreeMaker::minDr(double eta, double phi, const std::vector<pat::Photon>& phos) {
	float minDrVal = 100 ;
	for (int i = 0 ; i < phos.size() ; ++i)
		if (phos[i].et() > 40 && dR(eta, phi, phos[i].eta(), phos[i].phi()) < minDrVal) minDrVal = dR(eta, phi, phos[i].eta(), phos[i].phi()) ;
	return minDrVal ;
}

int ADDrootTreeMaker::DoPassDefaultPhotonIDcut(const pat::Photon& pho) {
	if (pho.hadronicOverEm() < 0.05 &&
			pho.trkSumPtHollowConeDR04() < 5. &&
			pho.ecalRecHitSumEtConeDR04() < 10. &&
			pho.hcalTowerSumEtConeDR04() < 5. ) return 1 ;
	return 0 ;
}

// ------------ method called once each job just before starting event loop  ------------
void 
ADDrootTreeMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ADDrootTreeMaker::endJob() {
	std::cout << "\n Total number of Event: " << nEvent_ ;
	std::string s;
	std::stringstream out;
	out << nEvent_;
	s = out.str();
	std::cout << "\n Total number of Event (string): " << s ;
	tree_->SetTitle(s.c_str()) ;
	cout << "\n Trigger: " << "  " << nPhotonEvent_ << "  " << nIsoPho40_ << "  " << nPho25_ << "  " << nDoublePho20_ ;
}



//define this as a plug-in
DEFINE_FWK_MODULE(ADDrootTreeMaker);
