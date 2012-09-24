#define CutEfficiencyLoop_cxx
#include "CutEfficiencyLoop.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "AcceptanceCuts.C"
#include "ReconEfficiency.C"
#include  <iostream>

void CutEfficiencyLoop::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L CutEfficiencyLoop.C
//      Root > CutEfficiencyLoop t
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

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      h_Diphoton_Minv->Fill(Diphoton_Minv,MCPUWeight);
 
     if (passesAcptCuts(Diphoton_Minv,Photon1_detEta,Photon1_pt) && passesAcptCuts(Diphoton_Minv,Photon2_detEta,Photon2_pt)){

      h_Diphoton_Minv_Acceptance->Fill(Diphoton_Minv,MCPUWeight);

      if (PassEffCutsminusPixel(Photon1_pt,Photon1_hasPixelSeed,Photon1_sigmaIetaIeta,rho25,Photon1_trkIsoSumPtHollow04,Photon1_hcalIso04,Photon1_ecalIso04,Photon1_hadOverEm) && PassEffCutsminusPixel(Photon2_pt,Photon2_hasPixelSeed,Photon2_sigmaIetaIeta,rho25,Photon2_trkIsoSumPtHollow04,Photon2_hcalIso04,Photon2_ecalIso04,Photon2_hadOverEm)){
	h_Diphoton_NoPixelSeed_Before->Fill(Diphoton_Minv,MCPUWeight); 
	if (!PassesPixelSeed(Photon1_hasPixelSeed) && !PassesPixelSeed(Photon2_hasPixelSeed)) {
	  h_Diphoton_NoPixelSeed_After->Fill(Diphoton_Minv,MCPUWeight);
	  	  
	}
      }
     
   
      if (PassesEffCutsminusHoverE(Photon1_pt,Photon1_hasPixelSeed,Photon1_sigmaIetaIeta,rho25,Photon1_trkIsoSumPtHollow04,Photon1_hcalIso04,Photon1_ecalIso04,Photon1_hadOverEm) && PassesEffCutsminusHoverE(Photon2_pt,Photon2_hasPixelSeed,Photon2_sigmaIetaIeta,rho25,Photon2_trkIsoSumPtHollow04,Photon2_hcalIso04,Photon2_ecalIso04,Photon2_hadOverEm)){
	h_Diphoton_Minv_HoverE_Before->Fill(Diphoton_Minv,MCPUWeight);
	if ( passesHoverECut(Photon1_hadOverEm) && passesHoverECut(Photon2_hadOverEm)) {
	  h_Diphoton_Minv_HoverE_After->Fill(Diphoton_Minv,MCPUWeight);
	  
        }
      }
      

      if (PassesEffCutsminusSigmaIetaIeta(Photon1_pt,Photon1_hasPixelSeed,Photon1_sigmaIetaIeta,rho25,Photon1_trkIsoSumPtHollow04,Photon1_hcalIso04,Photon1_ecalIso04,Photon1_hadOverEm) && PassesEffCutsminusSigmaIetaIeta(Photon2_pt,Photon2_hasPixelSeed,Photon2_sigmaIetaIeta,rho25,Photon2_trkIsoSumPtHollow04,Photon2_hcalIso04,Photon2_ecalIso04,Photon2_hadOverEm)){
        h_Diphoton_Minv_SigmaIetaIeta_Before->Fill(Diphoton_Minv,MCPUWeight);
        if ( passesHoverECut(Photon1_sigmaIetaIeta) &&  passesHoverECut(Photon2_sigmaIetaIeta)) {
	  h_Diphoton_Minv_SigmaIetaIeta_After->Fill(Diphoton_Minv,MCPUWeight);	 
        }
      }
     
      if (PassesEffCutsminstrckIso(Photon1_pt,Photon1_hasPixelSeed,Photon1_sigmaIetaIeta,rho25,Photon1_trkIsoSumPtHollow04,Photon1_hcalIso04,Photon1_ecalIso04,Photon1_hadOverEm) && PassesEffCutsminstrckIso(Photon2_pt,Photon2_hasPixelSeed,Photon2_sigmaIetaIeta,rho25,Photon2_trkIsoSumPtHollow04,Photon2_hcalIso04,Photon2_ecalIso04,Photon2_hadOverEm)){
        h_Diphoton_Minv_TrackIso_Before->Fill(Diphoton_Minv,MCPUWeight);
        if (trckSumPtHollow04Cut(Photon1_trkIsoSumPtHollow04,Photon1_pt,rho25) && trckSumPtHollow04Cut(Photon2_trkIsoSumPtHollow04,Photon2_pt,rho25)) {
          h_Diphoton_Minv_TrackIso_After->Fill(Diphoton_Minv,MCPUWeight);
        }
      }
   
      if (PassesEffCutsminusEcal(Photon1_pt,Photon1_hasPixelSeed,Photon1_sigmaIetaIeta,rho25,Photon1_trkIsoSumPtHollow04,Photon1_hcalIso04,Photon1_ecalIso04,Photon1_hadOverEm) &&  PassesEffCutsminusEcal(Photon2_pt,Photon2_hasPixelSeed,Photon2_sigmaIetaIeta,rho25,Photon2_trkIsoSumPtHollow04,Photon2_hcalIso04,Photon2_ecalIso04,Photon2_hadOverEm)){
        h_Diphoton_Minv_EcalIso_Before->Fill(Diphoton_Minv,MCPUWeight);
        if (EcalIsoCut(Photon1_ecalIso04,Photon1_pt,rho25) && EcalIsoCut(Photon2_ecalIso04,Photon2_pt,rho25)) {
          h_Diphoton_Minv_EcalIso_After->Fill(Diphoton_Minv,MCPUWeight);
        }
      }

        
      if (PassesEffCutsminusHcal(Photon1_pt,Photon1_hasPixelSeed,Photon1_sigmaIetaIeta,rho25,Photon1_trkIsoSumPtHollow04,Photon1_hcalIso04,Photon1_ecalIso04,Photon1_hadOverEm) &&  PassesEffCutsminusHcal(Photon2_pt,Photon2_hasPixelSeed,Photon2_sigmaIetaIeta,rho25,Photon2_trkIsoSumPtHollow04,Photon2_hcalIso04,Photon2_ecalIso04,Photon2_hadOverEm)){
        h_Diphoton_Minv_HcalIso_Before->Fill(Diphoton_Minv,MCPUWeight);
        if (HcalIsoCut(Photon1_hcalIso04,Photon1_pt,rho25) && HcalIsoCut(Photon2_hcalIso04,Photon2_pt,rho25)) {
          h_Diphoton_Minv_HcalIso_After->Fill(Diphoton_Minv,MCPUWeight);
        }
      }
     }
   }

    
     _outfilename->cd();
     h_Diphoton_Minv->Write();
     h_Diphoton_Minv_Acceptance->Write();
     h_Diphoton_Minv_HoverE_Before->Write();
     h_Diphoton_Minv_TrackIso_Before->Write();
     h_Diphoton_Minv_HcalIso_Before->Write();
     h_Diphoton_NoPixelSeed_Before->Write();
     h_Diphoton_Minv_EcalIso_Before->Write();
     h_Diphoton_Minv_SigmaIetaIeta_Before->Write();
     h_Diphoton_Minv_HoverE_After->Write();
     h_Diphoton_Minv_TrackIso_After->Write();
     h_Diphoton_Minv_HcalIso_After->Write();
     h_Diphoton_NoPixelSeed_After->Write();
     h_Diphoton_Minv_EcalIso_After->Write();
     h_Diphoton_Minv_SigmaIetaIeta_After->Write();
 
}
