#ifndef JET_INFO_INC
#define JET_INFO_INC

//***********************************************
// Diphoton Info
//
//************************************************

#include <string>
#include <iostream>

#include "DataFormats/Math/interface/deltaPhi.h"
#include "TMath.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/Candidate/interface/LeafCandidate.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TLorentzVector.h"

// for jets
#include "DataFormats/PatCandidates/interface/Jet.h"

namespace ExoDiPhotons{


  // diphoton info: Minv, q_T, delta phi, etc
  const Int_t maxJets = 500;
  struct jetInfo_t{

    Int_t nJets;
    Double_t pt[maxJets];
    Double_t eta[maxJets];
    Double_t phi[maxJets];
    Double_t mass[maxJets];
    Double_t energy[maxJets];
    Bool_t passLooseID[maxJets];
    Bool_t passTightID[maxJets];


  };

  // std::string jetInfoBranchDefString("Minv/D:qt:deltaPhi:deltaEta:deltaR:deltaROld:cosThetaStar:cosThetaStarOld");
  // std::string jetInfoBranchDefString("nJets/I");
  // std::string jetInfoBranchDefString("nJets/I:pt[nJets]/D");
  

  std::pair<bool,bool> jetID(const pat::Jet& pfjet){

    float NHF = pfjet.neutralHadronEnergyFraction();
    float NEMF = pfjet.neutralEmEnergyFraction();
    float CHF = pfjet.chargedHadronEnergyFraction();
    // float MUF = pfjet.muonEnergyFraction();
    float CEMF = pfjet.chargedEmEnergyFraction();
    int NumConst = pfjet.chargedMultiplicity()+pfjet.neutralMultiplicity();
    int NumNeutralParticles =pfjet.neutralMultiplicity();
    int CHM = pfjet.chargedMultiplicity();

    bool looseJetID = false;
    bool tightJetID = false;

    if ( fabs(pfjet.eta()) <=3. ){
      looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((fabs(pfjet.eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || fabs(pfjet.eta())>2.4) && fabs(pfjet.eta())<=3.0;

      tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((fabs(pfjet.eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || fabs(pfjet.eta())>2.4) && fabs(pfjet.eta())<=3.0;

    }
    else if ( fabs(pfjet.eta()) > 3.){
      looseJetID = (NEMF<0.90 && NumNeutralParticles>10 && fabs(pfjet.eta())>3.0 );
      tightJetID = (NEMF<0.90 && NumNeutralParticles>10 && fabs(pfjet.eta())>3.0 );
    }

    return std::make_pair(looseJetID,tightJetID);

  }

  void FillJetInfo(jetInfo_t &fJetInfo, const edm::View<pat::Jet>* jets) {

    fJetInfo.nJets = (Int_t)jets->size();
    for (int i=0; i<fJetInfo.nJets; i++){
      // std::cout << "jet pt " << jets->at(i).pt() << std::endl;
      pat::Jet jet = jets->at(i);
      fJetInfo.pt[i] = jet.pt();
      fJetInfo.eta[i] = jet.eta();
      fJetInfo.phi[i] = jet.phi();
      fJetInfo.mass[i] = jet.mass();
      fJetInfo.energy[i] = jet.energy();
      
      Bool_t loose;
      Bool_t tight;
      std::tie(loose,tight) = ExoDiPhotons::jetID(jet);
      fJetInfo.passLooseID[i] = loose;
      fJetInfo.passTightID[i] = tight;
    }


  }
  // void InitJetInfo(jetInfo_t &fJetInfo, const edm::View<pat::Jet>* jets) {

  //   fJetInfo.nJets = (Int_t)jets->size();
  //   // fJetInfo.pt = *(new Double_t [fJetInfo.nJets]);
    
  // }

  

} //end of namespace


#endif