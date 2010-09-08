#define fTreeData_cxx
#include "fTreeData.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


void fTreeData::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L fTreeData.C
//      Root > fTreeData t
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
      if (ientry==0) {
	std::cout << "Photon pT cuts: " << _cutPhoton1Pt << " " << _cutPhoton2Pt << std::endl;
      }

      if ( Photon1_pt>_cutPhoton1Pt && Photon2_pt>_cutPhoton2Pt ) {

	//	if (TrigHLT_HLT_MinBiasBSC>0) h_TrigHLT->Fill(0);
	//	if (TrigHLT_HLT_MinBiasBSC_NoBPTX>0) h_TrigHLT->Fill(1);
	//	if (TrigHLT_HLT_MinBiasBSC_OR>0) h_TrigHLT->Fill(2);
	//	if (TrigHLT_HLT_L1_BscMinBiasOR_BptxPlusORMinus>0) h_TrigHLT->Fill(3);
	if (TrigHLT_HLT_L1SingleEG2>0) h_TrigHLT->Fill(0);
	if (TrigHLT_HLT_L1SingleEG5>0) h_TrigHLT->Fill(1);
	if (TrigHLT_HLT_L1SingleEG8>0) h_TrigHLT->Fill(2);
	if (TrigHLT_HLT_L1DoubleEG5>0) h_TrigHLT->Fill(3);
	if (TrigHLT_HLT_Photon10_L1R>0) h_TrigHLT->Fill(4);
	if (TrigHLT_HLT_Photon10_Cleaned_L1R>0) h_TrigHLT->Fill(5);
	if (TrigHLT_HLT_Photon15_L1R>0) h_TrigHLT->Fill(6);
	if (TrigHLT_HLT_Photon15_Cleaned_L1R>0) h_TrigHLT->Fill(7);
	if (TrigHLT_HLT_Photon15_LooseEcalIso_L1R>0) h_TrigHLT->Fill(8);
	if (TrigHLT_HLT_Photon15_LooseEcalIso_Cleaned_L1R>0) h_TrigHLT->Fill(9);
	if (TrigHLT_HLT_Photon15_TrackIso_L1R>0) h_TrigHLT->Fill(10);
	if (TrigHLT_HLT_Photon15_TrackIso_Cleaned_L1R>0) h_TrigHLT->Fill(11);
	if (TrigHLT_HLT_Photon20_L1R>0) h_TrigHLT->Fill(12);
	if (TrigHLT_HLT_Photon20_Cleaned_L1R>0) h_TrigHLT->Fill(13);
	if (TrigHLT_HLT_Photon30_L1R>0) h_TrigHLT->Fill(14);
	if (TrigHLT_HLT_Photon30_Cleaned_L1R>0) h_TrigHLT->Fill(15);
	if (TrigHLT_HLT_Photon30_L1R_8E29>0) h_TrigHLT->Fill(16);
	if (TrigHLT_HLT_Photon50_L1R>0) h_TrigHLT->Fill(17);
	if (TrigHLT_HLT_Photon50_Cleaned_L1R>0) h_TrigHLT->Fill(18);
	if (TrigHLT_HLT_DoublePhoton5_L1R>0) h_TrigHLT->Fill(19);
	if (TrigHLT_HLT_DoublePhoton10_L1R>0) h_TrigHLT->Fill(20);
	if (TrigHLT_HLT_DoublePhoton15_L1R>0) h_TrigHLT->Fill(21);
	if (TrigHLT_HLT_DoublePhoton20_L1R>0) h_TrigHLT->Fill(22);


	h_Photon1_pt->Fill(Photon1_pt);
	h_Photon1_eta->Fill(Photon1_eta);
	h_Photon1_phi->Fill(Photon1_phi);

   	h_Photon1_r9->Fill(Photon1_r9);
 	h_Photon1_sigmaIetaIeta->Fill(Photon1_sigmaIetaIeta);
 	h_Photon1_sigmaEtaEta->Fill(Photon1_sigmaEtaEta);

 	if (Photon1_swisscross>-100) h_Photon1_swisscross->Fill(Photon1_swisscross);
 	h_Photon1_severityLevel->Fill(Photon1_severityLevel);
 	h_Photon1_recHitFlag->Fill(Photon1_recHitFlag);
 	h_Photon1_maxRecHitTime->Fill(Photon1_maxRecHitTime);
	
 	h_Photon1_hadOverEm->Fill(Photon1_hadOverEm);
 	h_Photon1_hcalIso04->Fill(Photon1_hcalIso04);
 	h_Photon1_hcalIso03->Fill(Photon1_hcalIso03);
 	h_Photon1_ecalIso04->Fill(Photon1_ecalIso04);
 	h_Photon1_ecalIso03->Fill(Photon1_ecalIso03);

 	h_Photon1_trkIsoSumPtHollow04->Fill(Photon1_trkIsoSumPtHollow04);
 	h_Photon1_trkIsoSumPtSolid04->Fill(Photon1_trkIsoSumPtSolid04);
 	h_Photon1_trkIsoSumPtHollow03->Fill(Photon1_trkIsoSumPtHollow03);
 	h_Photon1_trkIsoSumPtSolid03->Fill(Photon1_trkIsoSumPtSolid03);

 	h_Photon1_trkIsoNtrksHollow04->Fill(Photon1_trkIsoNtrksHollow04);
 	h_Photon1_trkIsoNtrksSolid04->Fill(Photon1_trkIsoNtrksSolid04);
 	h_Photon1_trkIsoNtrksHollow03->Fill(Photon1_trkIsoNtrksHollow03);
 	h_Photon1_trkIsoNtrksSolid03->Fill(Photon1_trkIsoNtrksSolid03);

// 	h_Photon1_esRatio->Fill(Photon1_esRatio);

	h_Photon2_pt->Fill(Photon2_pt);
	h_Photon2_eta->Fill(Photon2_eta);
	h_Photon2_phi->Fill(Photon2_phi);

   	h_Photon2_r9->Fill(Photon2_r9);
 	h_Photon2_sigmaIetaIeta->Fill(Photon2_sigmaIetaIeta);
 	h_Photon2_sigmaEtaEta->Fill(Photon2_sigmaEtaEta);

 	if (Photon2_swisscross>-100) h_Photon2_swisscross->Fill(Photon2_swisscross);
 	h_Photon2_severityLevel->Fill(Photon2_severityLevel);
 	h_Photon2_recHitFlag->Fill(Photon2_recHitFlag);
 	h_Photon2_maxRecHitTime->Fill(Photon2_maxRecHitTime);
	
 	h_Photon2_hadOverEm->Fill(Photon2_hadOverEm);
 	h_Photon2_hcalIso04->Fill(Photon2_hcalIso04);
 	h_Photon2_hcalIso03->Fill(Photon2_hcalIso03);
 	h_Photon2_ecalIso04->Fill(Photon2_ecalIso04);
 	h_Photon2_ecalIso03->Fill(Photon2_ecalIso03);

 	h_Photon2_trkIsoSumPtHollow04->Fill(Photon2_trkIsoSumPtHollow04);
 	h_Photon2_trkIsoSumPtSolid04->Fill(Photon2_trkIsoSumPtSolid04);
 	h_Photon2_trkIsoSumPtHollow03->Fill(Photon2_trkIsoSumPtHollow03);
 	h_Photon2_trkIsoSumPtSolid03->Fill(Photon2_trkIsoSumPtSolid03);

 	h_Photon2_trkIsoNtrksHollow04->Fill(Photon2_trkIsoNtrksHollow04);
 	h_Photon2_trkIsoNtrksSolid04->Fill(Photon2_trkIsoNtrksSolid04);
 	h_Photon2_trkIsoNtrksHollow03->Fill(Photon2_trkIsoNtrksHollow03);
 	h_Photon2_trkIsoNtrksSolid03->Fill(Photon2_trkIsoNtrksSolid03);

// 	h_Photon2_esRatio->Fill(Photon2_esRatio);

        h_Diphoton_Minv->Fill(Diphoton_Minv);
        h_Diphoton_qt->Fill(Diphoton_qt);
        h_Diphoton_deltaPhi->Fill(Diphoton_deltaPhi);
        h_Diphoton_deltaEta->Fill(Diphoton_deltaEta);
        h_Diphoton_deltaR->Fill(Diphoton_deltaR);

      }

   }

   if ( _outF != 0 ) {
     _outF->cd();
     h_TrigHLT->Write();
     h_Diphoton_Minv->Write();
     h_Diphoton_qt->Write();
     h_Diphoton_deltaPhi->Write();
     h_Diphoton_deltaEta->Write();
     h_Diphoton_deltaR->Write();
     h_Photon1_pt->Write();
     h_Photon1_eta->Write();
     h_Photon1_phi->Write();
     h_Photon1_r9->Write();
     h_Photon1_sigmaIetaIeta->Write();
     h_Photon1_sigmaEtaEta->Write();
     h_Photon1_swisscross->Write();
     h_Photon1_severityLevel->Write();
     h_Photon1_recHitFlag->Write();
     h_Photon1_maxRecHitTime->Write();
     h_Photon1_hadOverEm->Write();
     h_Photon1_hcalIso04->Write();
     h_Photon1_hcalIso03->Write();
     h_Photon1_ecalIso04->Write();
     h_Photon1_ecalIso03->Write();
     h_Photon1_trkIsoSumPtHollow04->Write();
     h_Photon1_trkIsoSumPtSolid04->Write();
     h_Photon1_trkIsoSumPtHollow03->Write();
     h_Photon1_trkIsoSumPtSolid03->Write();
     h_Photon1_trkIsoNtrksHollow04->Write();
     h_Photon1_trkIsoNtrksSolid04->Write();
     h_Photon1_trkIsoNtrksHollow03->Write();
     h_Photon1_trkIsoNtrksSolid03->Write();
     h_Photon2_pt->Write();
     h_Photon2_eta->Write();
     h_Photon2_phi->Write();
     h_Photon2_r9->Write();
     h_Photon2_sigmaIetaIeta->Write();
     h_Photon2_sigmaEtaEta->Write();
     h_Photon2_swisscross->Write();
     h_Photon2_severityLevel->Write();
     h_Photon2_recHitFlag->Write();
     h_Photon2_maxRecHitTime->Write();
     h_Photon2_hadOverEm->Write();
     h_Photon2_hcalIso04->Write();
     h_Photon2_hcalIso03->Write();
     h_Photon2_ecalIso04->Write();
     h_Photon2_ecalIso03->Write();
     h_Photon2_trkIsoSumPtHollow04->Write();
     h_Photon2_trkIsoSumPtSolid04->Write();
     h_Photon2_trkIsoSumPtHollow03->Write();
     h_Photon2_trkIsoSumPtSolid03->Write();
     h_Photon2_trkIsoNtrksHollow04->Write();
     h_Photon2_trkIsoNtrksSolid04->Write();
     h_Photon2_trkIsoNtrksHollow03->Write();
     h_Photon2_trkIsoNtrksSolid03->Write();
     _outF->Write();
     _outF->Close();
   }
   return;
}
