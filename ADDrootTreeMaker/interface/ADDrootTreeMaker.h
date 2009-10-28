#ifndef __CALO_CALIB_H__
#define __CALO_CALIB_H__

#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "Math/LorentzVector.h"
#include "TH1F.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"


//
// foward definitions
//

class TFile;
class TTree;

#include "TObject.h"


//
// class decleration
//

class ADDrootTreeMaker : public edm::EDAnalyzer {
 public:
  explicit ADDrootTreeMaker(const edm::ParameterSet&);
  ~ADDrootTreeMaker();
//  float InvariantMass(const LorentzVector&, const LorentzVector&) ;
  float InvariantMass(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >&, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >&) ;
	float InvariantMass(double, double, double, double, double, double, double, double) ;
	float dPhi(double, double) ;
	float dR(double, double, double, double) ;
	float minDr(double, double, const std::vector<pat::Jet>&) ;
	float minDr(double, double, const std::vector<pat::Photon>&) ;
	int DoPassDefaultPhotonIDcut(const pat::Photon&) ;
 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;



  // ----------member data ---------------------------
  

  TFile* rootfile_;
  TTree* tree_;
  int evt_run_, evt_event_;
	
	static const int maxpgens=1000;
  int npgens_;
  int pgen_id_[maxpgens], pgen_energy_[maxpgens]  ;
	int nEvent_ ;
	
	static const int maxPhotonPair=1000 ;
	int nPhotons_ ;
	int nPhotonPair_ ;
	int indPho1_phoPair_[maxPhotonPair], indPho2_phoPair_[maxPhotonPair] ;
	float invMass_phoPair_[maxPhotonPair] ;
	float dR_phoPair_[maxPhotonPair] ;
	float dPhi_phoPair_[maxPhotonPair] ;
	float genInvMass_phoPair_[maxPhotonPair] ;
	
	int nRecoPhotonPair_ ;
	int indPho1_phoPair_recoPho_[maxPhotonPair], indPho2_phoPair_recoPho_[maxPhotonPair] ;
	float invMass_phoPair_recoPho_[maxPhotonPair] ;
	float dR_phoPair_recoPho_[maxPhotonPair] ;
	float dPhi_phoPair_recoPho_[maxPhotonPair] ;

	static const int maxPhoJetPair = 1000 ;
	int nPhoJetPair_ ;
	float minDrJet_phoJetPair_[maxPhoJetPair];
	int indJet_phoJetPair_[maxPhoJetPair];
	float minDrPho_phoJetPair_[maxPhoJetPair];
	int indPho_phoJetPair_[maxPhoJetPair] ;
	
	float invMass_phoJetPair_[maxPhoJetPair] ;
	float dR_phoJetPair_[maxPhoJetPair] ;
	float dPhi_phoJetPair_[maxPhoJetPair] ;
	float genInvMass_phoJetPair_[maxPhoJetPair] ;

	static const int maxGenPho = 10000 ;
	int nGenPho_ ;
	float genE_pho_[maxGenPho], genPt_pho_[maxGenPho], genEta_pho_[maxGenPho], genPhi_pho_[maxGenPho] ;
	int genPdgId_pho_[maxGenPho], motherGenPdgId_pho_[maxGenPho], motherGenStatus_pho_[maxGenPho];
	static const int maxPhos = 100 ;
	int nPhos_ ;
	float eta_recPho_[maxPhos], phi_recPho_[maxPhos], pt_recPho_[maxPhos], e_recPho_[maxPhos], p_recPho_[maxPhos], e5x5_recPho_[maxPhos], e3x3_recPho_[maxPhos] ;
	int nRecoPhos_ ;
	float eta_recoPho_[maxPhos], phi_recoPho_[maxPhos], pt_recoPho_[maxPhos], e_recoPho_[maxPhos];
	int doPassDefaultPhotonIDcut_[maxPhos] ;
	int isAlsoElectron_[maxPhos], nTrkSolidCone_[maxPhos], nTrkHollowCone_[maxPhos] ; 
	float isolationEcalRecHit_[maxPhos], isolationHcalRecHit_[maxPhos], isolationSolidTrkCone_[maxPhos], isolationHollowTrkCone_[maxPhos], r9_[maxPhos], hadronicOverEm_[maxPhos] ;
	int hasPixelSeed_[maxPhos] ;
	
	int nTrkSolidCone_recoPho_[maxPhos], nTrkHollowCone_recoPho_[maxPhos] ; 
	float isolationEcalRecHit_recoPho_[maxPhos], isolationHcalRecHit_recoPho_[maxPhos], isolationSolidTrkCone_recoPho_[maxPhos], isolationHollowTrkCone_recoPho_[maxPhos], r9_recoPho_[maxPhos], hadronicOverEm_recoPho_[maxPhos] ;
	int hasPixelSeed_recoPho_[maxPhos] ;
	
	float genEta_recPho_[maxPhos], genPhi_recPho_[maxPhos], genPt_recPho_[maxPhos], genE_recPho_[maxPhos], genP_recPho_[maxPhos], genMotherPdg_recPho_[maxPhos] ;
	int genPdgId_recPho_[maxPhos], motherGenPdgId_recPho_[maxPhos] ;
	//float invMass_[maxPhotonPair] ;
	//static const int maxjets = 100 ;
	static const int maxJets = 100 ;
	int nJets_ ;
	float eta_recJet_[maxJets], phi_recJet_[maxJets], et_recJet_[maxJets], e_recJet_[maxJets], p_recJet_[maxJets], emf_recJet_[maxJets] ;
	float genEta_recJet_[maxJets], genPhi_recJet_[maxJets], genEt_recJet_[maxJets], genE_recJet_[maxJets], genP_recJet_[maxJets] ;
	
	float pThat_ ;
	int processId_ ;
	
	static const int maxMCpart = 3000 ;
	int nMCpart_ ;
	float eta_genPart_[maxMCpart], phi_genPart_[maxMCpart], pt_genPart_[maxMCpart], e_genPart_[maxMCpart] ;
	int pdgId_genPart_[maxMCpart], status_genPart_[maxMCpart], motherPdgId_genPart_[maxMCpart], motherStatus_genPart_[maxMCpart], grandMotherPdgId_genPart_[maxMCpart], grandMotherStatus_genPart_[maxMCpart] ;
	
	static const int maxGenJet = 1000 ;
	int nGenJet_ ;
	float eta_genJet_[maxGenJet], phi_genJet_[maxGenJet], et_genJet_[maxGenJet], e_genJet_[maxGenJet] ; 
	
	int HLT_IsoPhoton40_L1R_;
	int HLT_Photon25_L1R_ ;
	int HLT_DoubleIsoPhoton20_L1R_ ;

	int nPhotonEvent_, nIsoPho40_, nPho25_, nDoublePho20_ ;

	static const int maxMCpart3 = 20 ;
	int nMCpart3_ ;
	float genEta_part3_[maxMCpart3], genPhi_part3_[maxMCpart3], genPt_part3_[maxMCpart3] ;
	int pdgId_part3_[maxMCpart3] ;
		
	static const int nMaxPi0 = 1000 ;
	int nPi0_ ;
	float eta_pi0_[nMaxPi0], phi_pi0_[nMaxPi0], pt_pi0_[nMaxPi0] ;
	
	TH1F* hNevent_ ;
	std::string fOutName_ ;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//


#endif
