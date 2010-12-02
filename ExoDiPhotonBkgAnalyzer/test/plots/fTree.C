#define fTree_cxx
#include "fTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


void fTree::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L fTree.C
//      Root > fTree t
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

   Int_t nPass=0;
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if (ientry==0) {
	std::cout << "TOTAL ENTRIES: " << nentries << std::endl;
	std::cout << "Photon pT cuts: " << _cutPhoton1Pt << " " << _cutPhoton2Pt << " " << _cutEta << std::endl;
      }

      // filter out diphoton events if PhotonJet
      if (_filterGen && (GenEvent_signalProcessId == 18)) continue;
      
      if ( Photon1_pt>_cutPhoton1Pt && Photon2_pt>_cutPhoton2Pt 
	   && fabs(Photon1_detEta)<_cutEta && fabs(Photon2_detEta)<_cutEta 
	   && (fabs(Photon1_detEta)<1.4442 || fabs(Photon1_detEta)>1.566)
	   && (fabs(Photon2_detEta)<1.4442 || fabs(Photon2_detEta)>1.566)
	   ) {

	if (_onlyEB) {	  
	  if (!Photon1_isEB || !Photon2_isEB) continue;
	}

	nPass++;

	// 	if (Diphoton_Minv>200.0) {
	// 	  std::cout << "High Mass Event: "
	// 		    << Diphoton_Minv 
	// 		    << " : " << Photon1_isEB << Photon1_isEE
	// 		    << " : " << Photon2_isEB << Photon2_isEE
	
	// 	    << " : " << Event_run << ":" << Event_LS << ":" << Event_evnum 
	// 	    << std::endl;
	// 	}
	
	//	std::cout //<< "Run:Lumi:Event " 
	//	  << Event_run << ":" << Event_LS << ":" << Event_evnum 
	//	  << "GenSignalProcessID \t"  << GenEvent_signalProcessId 
	//	  << std::endl;

	// 	std::cout 
	// 	  << "\tMinv:pT1:pT2:eta1:eta2           \t"
	// 	  << Diphoton_Minv << ":"
	// 	  << Photon1_pt << ":"
	// 	  << Photon2_pt << ":"
	// 	  << Photon1_detEta << ":"
	// 	  << Photon2_detEta 
	// 	  << std::endl; std::cout
	// 	  << "\thoe1:hoe2:sieie1:sieie2          \t"
	// 	  << Photon1_hadOverEm << ":"
	// 	  << Photon2_hadOverEm << ":"
	// 	  << Photon1_sigmaIetaIeta << ":"
	// 	  << Photon2_sigmaIetaIeta 
	// 	  << std::endl;	std::cout 
	// 	  << "\tecal1:ecal2:hcal1:hcal2:trk1:trk2\t"
	// 	  << Photon1_ecalIso04 << ":" 
	// 	  << Photon2_ecalIso04 << ":" 
	// 	  << Photon1_hcalIso04 << ":" 
	// 	  << Photon2_hcalIso04 << ":" 
	// 	  << Photon1_trkIsoSumPtHollow04 << ":" 
	// 	  << Photon2_trkIsoSumPtHollow04  
	// 	  << std::endl;	std::cout 
	// 	  << "\tVtx_vz                           \t"
	// 	  << Vtx_vz << ":"
	// 	  << std::endl;
	
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
        if (TrigHLT_HLT_Photon17_Isol_SC17HE_L1R_v1>0) h_TrigHLT->Fill(12);
        if (TrigHLT_HLT_Photon17_SC17HE_L1R_v1>0) h_TrigHLT->Fill(13);
        if (TrigHLT_HLT_Photon20_L1R>0) h_TrigHLT->Fill(14);
        if (TrigHLT_HLT_Photon20_Cleaned_L1R>0) h_TrigHLT->Fill(15);
        if (TrigHLT_HLT_Photon20_NoHE_L1R>0) h_TrigHLT->Fill(16);
        if (TrigHLT_HLT_Photon22_SC22HE_L1R_v1>0) h_TrigHLT->Fill(17);
        if (TrigHLT_HLT_Photon25_Cleaned_L1R>0) h_TrigHLT->Fill(18);
        if (TrigHLT_HLT_Photon30_L1R>0) h_TrigHLT->Fill(19);
        if (TrigHLT_HLT_Photon30_Cleaned_L1R>0) h_TrigHLT->Fill(20);
        if (TrigHLT_HLT_Photon30_L1R_8E29>0) h_TrigHLT->Fill(21);
        if (TrigHLT_HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1>0) h_TrigHLT->Fill(22);
        if (TrigHLT_HLT_Photon35_Isol_Cleaned_L1R_v1>0) h_TrigHLT->Fill(23);
        if (TrigHLT_HLT_Photon40_CaloId_Cleaned_L1R_v1>0) h_TrigHLT->Fill(24);
        if (TrigHLT_HLT_Photon40_Isol_Cleaned_L1R_v1>0) h_TrigHLT->Fill(25);
        if (TrigHLT_HLT_Photon50_L1R>0) h_TrigHLT->Fill(26);
        if (TrigHLT_HLT_Photon50_Cleaned_L1R>0) h_TrigHLT->Fill(27);
        if (TrigHLT_HLT_Photon50_Cleaned_L1R_v1>0) h_TrigHLT->Fill(28);
        if (TrigHLT_HLT_Photon50_NoHE_L1R>0) h_TrigHLT->Fill(29);
        if (TrigHLT_HLT_Photon50_NoHE_Cleaned_L1R>0) h_TrigHLT->Fill(30);
        if (TrigHLT_HLT_Photon70_Cleaned_L1R_v1>0) h_TrigHLT->Fill(31);
        if (TrigHLT_HLT_Photon70_NoHE_Cleaned_L1R_v1>0) h_TrigHLT->Fill(32);
        if (TrigHLT_HLT_Photon100_NoHE_Cleaned_L1R_v1>0) h_TrigHLT->Fill(33);
        if (TrigHLT_HLT_Photon110_NoHE_Cleaned_L1R_v1>0) h_TrigHLT->Fill(34);
        if (TrigHLT_HLT_DoublePhoton5_L1R>0) h_TrigHLT->Fill(35);
        if (TrigHLT_HLT_DoublePhoton5_CEP_L1R>0) h_TrigHLT->Fill(36);
        if (TrigHLT_HLT_DoublePhoton5_CEP_L1R_v3>0) h_TrigHLT->Fill(37);
        if (TrigHLT_HLT_DoublePhoton5_Jpsi_L1R>0) h_TrigHLT->Fill(38);
        if (TrigHLT_HLT_DoublePhoton5_Upsilon_L1R>0) h_TrigHLT->Fill(39);
        if (TrigHLT_HLT_DoublePhoton10_L1R>0) h_TrigHLT->Fill(40);
        if (TrigHLT_HLT_DoublePhoton15_L1R>0) h_TrigHLT->Fill(41);
        if (TrigHLT_HLT_DoublePhoton17_L1R>0) h_TrigHLT->Fill(42);
        if (TrigHLT_HLT_DoublePhoton17_SingleIsol_L1R_v1>0) h_TrigHLT->Fill(43);
        if (TrigHLT_HLT_DoublePhoton20_L1R>0) h_TrigHLT->Fill(44);
        if (TrigHLT_HLT_DoublePhoton22_L1R_v1>0) h_TrigHLT->Fill(45);

	h_Photon1_pt->Fill(Photon1_pt);
	h_Photon1_eta->Fill(Photon1_eta);
	h_Photon1_phi->Fill(Photon1_phi);

   	h_Photon1_r9->Fill(Photon1_r9);
 	h_Photon1_sigmaIetaIeta->Fill(Photon1_sigmaIetaIeta);
 	h_Photon1_sigmaEtaEta->Fill(Photon1_sigmaEtaEta);

 	if (Photon1_swisscross>-100) h_Photon1_swisscross->Fill(Photon1_swisscross);
	if (Photon1_e2e9>0) h_Photon1_e2e9->Fill(Photon1_e2e9);
	if (Photon1_e4x4!=0 && Photon1_e4x4>-100 && Photon1_e2x2>-100) h_Photon1_e2x2e4x4->Fill(Photon1_e2x2/Photon1_e4x4);
 	h_Photon1_severityLevel->Fill(Photon1_severityLevel);
 	h_Photon1_recHitFlag->Fill(Photon1_recHitFlag);
 	h_Photon1_maxRecHitTime->Fill(Photon1_maxRecHitTime);
 	h_Photon1_maxRecHitTime_wide->Fill(Photon1_maxRecHitTime);
	
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
	if (Photon2_e2e9>0) h_Photon2_e2e9->Fill(Photon2_e2e9);
	if (Photon2_e4x4!=0 && Photon2_e4x4>-100 && Photon2_e2x2>-100) h_Photon2_e2x2e4x4->Fill(Photon2_e2x2/Photon2_e4x4);
 	h_Photon2_severityLevel->Fill(Photon2_severityLevel);
 	h_Photon2_recHitFlag->Fill(Photon2_recHitFlag);
 	h_Photon2_maxRecHitTime->Fill(Photon2_maxRecHitTime);
 	h_Photon2_maxRecHitTime_wide->Fill(Photon2_maxRecHitTime);
	
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

        if (Diphoton_Minv>200) h_Diphoton_Minv_high->Fill(Diphoton_Minv);

        h_Diphoton_qt->Fill(Diphoton_qt);
        h_Diphoton_deltaPhi->Fill(Diphoton_deltaPhi);
        h_Diphoton_deltaEta->Fill(Diphoton_deltaEta);
        h_Diphoton_deltaR->Fill(Diphoton_deltaR);
        h_Diphoton_cosThetaStar->Fill(Diphoton_cosThetaStar);

      }

   }

   if ( _outF != 0 ) {
     _outF->cd();
     h_TrigHLT->Write();
     h_Diphoton_Minv->Write();
     h_Diphoton_Minv_high->Write();
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
     h_Photon1_maxRecHitTime_wide->Write();
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
     h_Photon2_maxRecHitTime_wide->Write();
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
   std::cout << "PASSING " << nPass << std::endl;
   return;
}
