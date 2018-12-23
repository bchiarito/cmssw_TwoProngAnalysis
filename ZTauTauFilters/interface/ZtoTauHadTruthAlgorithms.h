#ifndef ZTOTAUHADTRUTHALGORITHMS_H
#define ZTOTAUHADTRUTHALGORITHMS_H

// cpp includes
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
// ROOT includes
#include "TVector2.h"
#include "TLorentzVector.h"
// cmssw framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
// trigger inlcudes
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
// pat includes
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
// edm utilities includes
#include "DataFormats/Math/interface/deltaR.h"

using std::string;
using std::vector;
using std::max;
using std::string;

namespace TauHadFilters
{
  // forward declare all functions
  bool isHadronicTau(const reco::GenParticle *);
  bool isAncestorOfZ(const reco::Candidate *);
  bool isAncestorOfTau(const reco::Candidate *);
  bool notTerminalTau(const reco::Candidate *);
  const reco::Candidate * tauDaughter(const reco::Candidate *);
  double ZDecayType(edm::Handle<vector<reco::GenParticle>>&);
  vector<string> getDecay(const reco::Candidate &, int flag=0);
  vector<const reco::Candidate *> getLeptonObjects(const reco::Candidate &);

  bool isHadronicTau(const reco::GenParticle * genparticle)
  {
    vector<string> leptons;
    if (abs(genparticle->pdgId()) != 15) return false;
    //for(size_t i=0; i<genparticle->numberOfMothers(); i++)
    //{
    //  if(isAncestorOfTau(genparticle->mother(i))) return false;
    //}
    if (notTerminalTau(genparticle)) return false;
    leptons = getDecay(*genparticle);
    if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() ||
        std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() ||
        std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() ||
        std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() ||
        std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end() ||
        std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end() ||
        std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end() ||
        std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end()) return true;
    else return false;
  }

  bool isAncestorOfZ(const reco::Candidate * particle)
  {
    if(particle->pdgId() == 23) return true;
    for(size_t i=0; i<particle->numberOfMothers(); i++)
    {
      if(isAncestorOfZ(particle->mother(i))) return true;
    }
    return false;
  }

  bool isAncestorOfTau(const reco::Candidate * particle)
  {
    if(abs(particle->pdgId()) == 15) return true;
    for(size_t i=0; i<particle->numberOfMothers(); i++)
    {
      if(isAncestorOfTau(particle->mother(i))) return true;
    }
    return false;
  }

  bool notTerminalTau(const reco::Candidate * particle)
  {
    for (unsigned int i = 0; i < particle->numberOfDaughters(); i++) {
      const reco::Candidate * d = particle->daughter(i);
      if (abs(d->pdgId()) == 15) return true;  
    }
    return false;
  }

  const reco::Candidate * tauDaughter(const reco::Candidate * particle)
  {
    for (unsigned int i = 0; i < particle->numberOfDaughters(); i++) {
      const reco::Candidate * d = particle->daughter(i);
      if (abs(d->pdgId()) == 15) return d;
    }
    return NULL;
  }

  double ZDecayType(edm::Handle<vector<reco::GenParticle>>& genparticles)
  {
    double decayType = -10.0;
    vector<string> leptons;
    for (unsigned int i = 0; i < genparticles->size(); i++) {
      const reco::GenParticle & genparticle = (*genparticles)[i];
      if (genparticle.status() != 21 && genparticle.status() != 22) continue;
      leptons = getDecay(genparticle);
    }
    if (leptons.size() == 0) 
      decayType = 0;
    if (leptons.size() == 1) 
      decayType = -1;
    if (leptons.size() == 2) 
    {
      if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end()) decayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end()) decayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end()) decayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end()) decayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end()) decayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end()) decayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end()) decayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end()) decayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end()) decayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end()) decayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end()) decayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end()) decayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end()) decayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end()) decayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end()) decayType = 3;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end()) decayType = 3;

      else if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end()) decayType = 4;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end()) decayType = 4;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end()) decayType = 4;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end()) decayType = 4;
      else if (std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end()) decayType = 4;
      else if (std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end()) decayType = 4;
      else if (std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end()) decayType = 4;
      else if (std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end()) decayType = 4;

      else if (std::find(leptons.begin(), leptons.end(), "tau+had10") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) decayType = 5.2;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had1") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) decayType = 5.1;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had30") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) decayType = 5.4;
      else if (std::find(leptons.begin(), leptons.end(), "tau+had3") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) decayType = 5.3;
      else if (std::find(leptons.begin(), leptons.end(), "tau-had10") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end()) decayType = 5.2;
      else if (std::find(leptons.begin(), leptons.end(), "tau-had1") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end()) decayType = 5.1;
      else if (std::find(leptons.begin(), leptons.end(), "tau-had30") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end()) decayType = 5.4;
      else if (std::find(leptons.begin(), leptons.end(), "tau-had3") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end()) decayType = 5.3;

      else if (std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end()) decayType = 6;

      else if (std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) decayType = 7;

      else if (std::find(leptons.begin(), leptons.end(), "tau+e") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau-mu") != leptons.end()) decayType = 8;
      else if (std::find(leptons.begin(), leptons.end(), "tau-e") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "tau+mu") != leptons.end()) decayType = 8;

      else if (std::find(leptons.begin(), leptons.end(), "e+") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "e-") != leptons.end()) decayType = 1;

      else if (std::find(leptons.begin(), leptons.end(), "mu+") != leptons.end() &&
          std::find(leptons.begin(), leptons.end(), "mu-") != leptons.end()) decayType = 2;

      else
        decayType = 9;
    }
    if (leptons.size() > 2)
    {
      decayType = 10;
    }
    return decayType;
  }

  vector<string> getDecay(const reco::Candidate & genparticle, int flag)
  {
    vector<string> products;

    if (flag == 1) { // parent is tau
      if (genparticle.pdgId() == 111) { products.push_back("npion"); return products; }
      if (abs(genparticle.pdgId()) == 211) { products.push_back("cpion"); return products; }
    }

    // ignore quarks and gluons and their daughters, unless status 21
    if (genparticle.status() != 21 && genparticle.pdgId() >= 1 && genparticle.pdgId() <= 8) { return products; }
    if (genparticle.status() != 21 && genparticle.pdgId() == 21) { return products; }

    if (genparticle.pdgId() == 11) { products.push_back("e-"); return products; }
    if (genparticle.pdgId() == -11) { products.push_back("e+"); return products; }
    if (genparticle.pdgId() == 13) { products.push_back("mu-"); return products; }
    if (genparticle.pdgId() == -13) { products.push_back("mu+"); return products; }

    if (abs(genparticle.pdgId()) == 15) {
      vector<string> tau_products; 
      for (unsigned int j = 0; j < genparticle.numberOfDaughters(); j++) {
        const reco::Candidate* daughter = genparticle.daughter(j);
        vector<string> daughter_products = getDecay(*daughter, 1);
        tau_products.insert(tau_products.end(), daughter_products.begin(), daughter_products.end());
      }
    
      if (tau_products.size() == 0) {
        if (genparticle.pdgId() == 15) { products.push_back("tau-unid"); return products; }
        if (genparticle.pdgId() == -15) { products.push_back("tau+unid"); return products; }
      }

      if (tau_products.size() == 1 && tau_products[0] == "e-") { products.push_back("tau-e"); return products; }
      if (tau_products.size() == 1 && tau_products[0] == "e+") { products.push_back("tau+e"); return products; }
      if (tau_products.size() == 1 && tau_products[0] == "mu-") { products.push_back("tau-mu"); return products; }
      if (tau_products.size() == 1 && tau_products[0] == "mu+") { products.push_back("tau+mu"); return products; }
      if (tau_products.size() == 1 && tau_products[0] == "tau-e") { products.push_back("tau-e"); return products; }
      if (tau_products.size() == 1 && tau_products[0] == "tau+e") { products.push_back("tau+e"); return products; }
      if (tau_products.size() == 1 && tau_products[0] == "tau-mu") { products.push_back("tau-mu"); return products; }
      if (tau_products.size() == 1 && tau_products[0] == "tau+mu") { products.push_back("tau+mu"); return products; }
      if (tau_products.size() == 1 && tau_products[0] == "tau+had10") { products.push_back("tau+had10"); return products; }
      if (tau_products.size() == 1 && tau_products[0] == "tau+had1") { products.push_back("tau+had1"); return products; }
      if (tau_products.size() == 1 && tau_products[0] == "tau+had30") { products.push_back("tau+had30"); return products; }
      if (tau_products.size() == 1 && tau_products[0] == "tau+had3") { products.push_back("tau+had3"); return products; }
      if (tau_products.size() == 1 && tau_products[0] == "tau-had10") { products.push_back("tau-had10"); return products; }
      if (tau_products.size() == 1 && tau_products[0] == "tau-had1") { products.push_back("tau-had1"); return products; }
      if (tau_products.size() == 1 && tau_products[0] == "tau-had30") { products.push_back("tau-had30"); return products; }
      if (tau_products.size() == 1 && tau_products[0] == "tau-had3") { products.push_back("tau-had3"); return products; }
      if (tau_products.size() == 1 && tau_products[0] == "tau+unid") { products.push_back("tau+unid"); return products; }
      if (tau_products.size() == 1 && tau_products[0] == "tau-unid") { products.push_back("tau-unid"); return products; }

      for (string prod : tau_products) {
        if (prod != "cpion" && prod != "npion") {  
          if (genparticle.pdgId() == 15) { products.push_back("tau-unid"); return products; }
          if (genparticle.pdgId() == -15) { products.push_back("tau+unid"); return products; }
        }
      }

      bool neutral = false;
      if (std::find(tau_products.begin(), tau_products.end(), "npion") != tau_products.end() ) neutral = true;

      int pion_count = 0;
      while( std::find(tau_products.begin(), tau_products.end(), "cpion") != tau_products.end() ) 
      {
        tau_products.erase( std::find(tau_products.begin(), tau_products.end(), "cpion") );
        pion_count += 1;
      }

      int charge = 0;
      if (genparticle.pdgId() == 15) charge = -1;
      if (genparticle.pdgId() == -15) charge = 1;

      if      (charge>0 && neutral  && pion_count == 1) products.push_back("tau+had10");
      else if (charge>0 && !neutral && pion_count == 1) products.push_back("tau+had1");
      else if (charge>0 && neutral  && pion_count == 3) products.push_back("tau+had30");
      else if (charge>0 && !neutral && pion_count == 3) products.push_back("tau+had3");
      else if (charge<0 && neutral  && pion_count == 1) products.push_back("tau-had10");
      else if (charge<0 && !neutral && pion_count == 1) products.push_back("tau-had1");
      else if (charge<0 && neutral  && pion_count == 3) products.push_back("tau-had30");
      else if (charge<0 && !neutral && pion_count == 3) products.push_back("tau-had3");
      else if (charge>0)                                products.push_back("tau+unid");
      else if (charge<0)                                products.push_back("tau-unid");

      return products;
      
    }

    for (unsigned int j = 0; j < genparticle.numberOfDaughters(); j++) {
      const reco::Candidate* daughter = genparticle.daughter(j);
      vector<string> daughter_products = getDecay(*daughter);
      products.insert(products.end(), daughter_products.begin(), daughter_products.end());
    }
    
    return products;
  }

  vector<const reco::Candidate *> getLeptonObjects(const reco::Candidate & genparticle)
  {
    vector<const reco::Candidate *> products;

    // ignore quarks and gluons and their daughters, unless status 21
    if (genparticle.status() != 21 && genparticle.pdgId() >= 1 && genparticle.pdgId() <= 8) { return products; }
    if (genparticle.status() != 21 && genparticle.pdgId() == 21) { return products; }

    // if a lepton, return this lepton
    if (abs(genparticle.pdgId()) == 11) { products.push_back(&genparticle); return products; }
    if (abs(genparticle.pdgId()) == 13) { products.push_back(&genparticle); return products; }
    if (abs(genparticle.pdgId()) == 15) { products.push_back(&genparticle); return products; }

    // otherwise, look for leptons in daughters and return those
    for (unsigned int j = 0; j < genparticle.numberOfDaughters(); j++) {
      const reco::Candidate* daughter = genparticle.daughter(j);
      vector<const reco::Candidate *> daughter_products = getLeptonObjects(*daughter);
      products.insert(products.end(), daughter_products.begin(), daughter_products.end());
    }
    return products;
  }
}

#endif
