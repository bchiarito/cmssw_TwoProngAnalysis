#define CalcSignalEfficiency_cxx
#include "CalcSignalEfficiency.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>

using namespace std;

void CalcSignalEfficiency::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L CalcSignalEfficiency.C
//      Root > CalcSignalEfficiency t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   //std::cout<<nentries<<endl;
   int neventspass = 0;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if(fabs(Diphoton_deltaPhi) < 2.8) continue;

      //Start the cuts

      Float_t dphi1 = GenPhoton1_phi - Photon1_phi;
      if(dphi1 > 3.1415926) dphi1-=2.*3.1415926;
      if(dphi1 <= -3.1415926) dphi1+=2.*3.1415926;
      Float_t deta1 = GenPhoton1_eta - Photon1_eta;
      Float_t deltar1 = sqrt(deta1*deta1 + dphi1*dphi1);
      if(deltar1 > 0.2) continue;

      //Event cut
      if(Diphoton_Minv < 200.) continue;

      //Photon1 fiducial cuts
      if(fabs(Photon1_eta) > 1.442) continue;
      if(Photon1_pt < 80.) continue;

      //Photon1 id and isolation cuts
      if(Photon1_hadOverEm > 0.05) continue;
      if(Photon1_trkIsoSumPtHollow04 > (2.0+0.001*Photon1_pt+0.0167*rho25)) continue;
      if(Photon1_ecalIso04 > (4.2 + 0.006*Photon1_pt+0.183*rho25)) continue;
      if(Photon1_hcalIso04 > (2.2 + 0.0025*Photon1_pt+0.062*rho25)) continue;
      if(Photon1_sigmaIetaIeta > 0.011) continue;
      if(Photon1_hasPixelSeed) continue;

      Float_t dphi2 = GenPhoton2_phi - Photon2_phi;
      if(dphi2 > 3.1415926) dphi2-=2.*3.1415926;
      if(dphi2 <= -3.1415926) dphi2+=2.*3.1415926;
      Float_t deta2 = GenPhoton2_eta - Photon2_eta;
      Float_t deltar2 = sqrt(deta2*deta2 + dphi2*dphi2);
      if(deltar2 > 0.2) continue;

      //Photon2 fiducial cuts
      if(fabs(Photon2_eta) > 1.442) continue;
      if(Photon2_pt < 80.) continue;

      //Photon2 id and isolation cuts
      if(Photon2_hadOverEm > 0.05) continue;
      if(Photon2_trkIsoSumPtHollow04 > (2.0+0.001*Photon2_pt+0.0167*rho25)) continue;
      if(Photon2_ecalIso04 > (4.2 + 0.006*Photon2_pt+0.183*rho25)) continue;
      if(Photon2_hcalIso04 > (2.2 + 0.0025*Photon2_pt+0.062*rho25)) continue;
      if(Photon2_sigmaIetaIeta > 0.011) continue;
      if(Photon2_hasPixelSeed) continue;

      neventspass++;
   }

   Float_t effi = (neventspass*1.)/nentries;
   Float_t effierror = sqrt((effi*(1.-effi))/nentries);

   //std::cout<<"#entries "<<nentries<<std::endl;
   //std::cout<<"neventspass "<<neventspass<<std::endl;
   std::cout<<"effi "<<effi<<" +- "<<effierror<<std::endl;

}
