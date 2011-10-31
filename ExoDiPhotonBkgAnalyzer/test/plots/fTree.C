#define fTree_cxx
#include "fTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>

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

   //
   // fake rate function for TF/FT/FF weighting events
   //
   TFile f_fakeratefunc(_fakeRateFile.Data());
   //   //   TF1 *fake_rate_fn = (TF1*) f_fakeratefunc.Get("fit_pow2_pixelveto");
   //   TF1 *fake_rate_fn = (TF1*) f_fakeratefunc.Get("fit_pow_pv");
   //   TF1 *fake_rate_fn = (TF1*) f_fakeratefunc.Get("fit_pow_pv");
   //   TF1*  fR = new TF1("fR", "0.01598+2431.92/x^2.6771", 30, 1000) ;
   //   TF1* fake_rate_fn  = new TF1("fake_rate_fn", "0.01598+2431.92/x^2.6771", 30, 200) ;

   //   TF1* fake_rate_fn    = new TF1("fake_rate_fn", "(-2.49781e-02+2.38899e+00/x^8.24583e-01)", 30, 1000) ;
   TF1* fake_rate_fn  = new TF1("fake_rate_fn", "(3.90649e-02+1.70945e+03/x^2.85695e+00)", 30, 1000) ;
   TF1* fake_rate_fn_EE = new TF1("fake_rate_fn_EE", "(4.50000e-02+9.22729e+07/x^6.37759e+00)", 30, 1000) ;

   //
   // kfactor
   //
   //   TF1* fKfactor = new TF1("kfactor", "-6.08061+9.90359/x^(3.88349e-02)", 100, 2000) ;
   //   TF1* fKfactor = new TF1("kfactor", "2.91986557985310036-0.284606795636766441/TMath::Power(x, -2.32654814731609899e-01)");

   //
   // PU distribution for re-weighting MC
   //

   vector<double> weights;

   if (_reweightPU) {
     TFile f_pu(_puFile);
     TH1D* histo = (TH1D*)f_pu.Get("pileup");
     //   vector<double> weights = generate_flat10_weights(histo,_puHistMC);
     weights = generate_flat10_weights(histo,_puHistMC);

     // this is the one
     //   Double_t weights[25] =  { 0.17923132, 0.35685511, 0.92815860, 1.40599997, 1.89257346, 2.04946601, 1.97097349, 1.73835990, 1.20440035, 0.83068139, 0.59277453, 0.41333785, 0.27472192, 0.20563187, 0.17254676, 0.10959478, 0.07831745, 0.05460577, 0.04859357, 0.07010237, 4.03994868, 0.31272140, 0.07666291, 0.00954727, 0.00405157     };
     
     //   Double_t weights[25] =  {    0.23894408, 0.38916675, 0.98219276, 1.45306598, 1.91620181, 2.03977899, 1.93553725, 1.68969438, 1.16193943, 0.79736259, 0.56706082, 0.39465926, 0.26210803, 0.19630832, 0.16509289, 0.10531444, 0.07579143, 0.05340644, 0.04824095, 0.07100313, 4.19888260, 0.33559630, 0.08548517, 0.01113009, 0.00496592 };
     
     //   Double_t weights[25] =  { 0.12166539, 0.41845364, 0.95846996, 1.56405252, 2.01196932, 2.16243274, 2.03053270, 1.98707804, 1.44884805, 1.06033535, 0.74564233, 0.45954123, 0.30539659, 0.20309582, 0.12108392, 0.07919775, 0.04615414, 0.02982704, 0.01629000, 0.01079210, 0.00588811, 0.00345084, 0.00214830, 0.00103641, 0.00071530   };
     
     //   Double_t weights[25] =  {  0          ,0.181628314,0.666156665,1.255546193,1.51942196 ,1.45958309 ,1.230185811,0.994068323,0.766019388,0.670261417,0.622132937,0.640068892,0.825382922,1.139633547,2.32635392 ,2.812910793,3.375494193,	   ,0	   ,0	   ,0	   ,0	   ,0	   ,0          ,0             }; 
     
     //   Double_t weights[25] =  {    0.23894408, 0.38916675, 0.98219276, 1.45306598, 1.91620181, 2.03977899, 1.93553725, 1.68969438, 1.16193943, 0.79736259, 0.56706082, 0.39465926, 0.26210803, 0.19630832, 0.16509289, 0.10531444, 0.07579143, 0.05340644, 0.04824095, 0.07100313, 4.19888260, 0.33559630, 0.08548517, 0.01113009, 0.00496592 };
     
     float blah=0;
     //   for (unsigned int i=0; i< weights.size(); i++) {
     for (unsigned int i=0; i< 25; i++) {
       cout << "PU weight " << i << " " <<  weights[i] << endl;
       blah+=weights[i];
     }
     //   cout << "Average " <<  blah/25 << endl;
   }

   cout << " now to actually loooop..." << endl;
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //   for (Long64_t jentry=0; jentry<500;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (ientry==0) {
	std::cout << "TOTAL ENTRIES: " << nentries << std::endl;
	std::cout << "Photon pT cuts: " << _cutPhoton1Pt << " " << _cutPhoton2Pt << " " << _cutEta << std::endl;
	std::cout << "Diphoton Kfactor ? " << _Kfactor << endl;
	std::cout << "PU reweight ? " << _reweightPU << endl;
      }

      // filter out diphoton events if PhotonJet
      if (_filterGen && (GenEvent_signalProcessId == 18)) continue;

      // PU event weight
      Double_t puWeight = 1.0;
      //      cout << pu_n << endl;
      if (_reweightPU) {
	if (pu_n<0) {
	  puWeight = weights[0];
	}
	else if (pu_n<25) {
	  puWeight = weights[pu_n];
	  //      	puWeight = weights[Vtx_Nvtx];
	  //	cout << pu_n << " " << puWeight << endl;
	} else {
	  puWeight = weights[24];
	}
      }
      
      if (_categoryEBEE=="EBEB") {	  
	//        if ((fabs(Photon1_detEta)>1.4442) || (fabs(Photon2_detEta)>1.4442)) continue;
        if ((!Photon1_isEB) || (!Photon2_isEB)) continue;
      } else if (_categoryEBEE=="EEEE") {
	if (!Photon1_isEE || !Photon2_isEE) continue;
      } else if (_categoryEBEE=="EBEE") {
	if (( Photon1_isEB && Photon2_isEB ) || (Photon1_isEE && Photon2_isEE)) continue;
      } else if (_categoryEBEE=="NOTEE") {
	if ( Photon1_isEE && Photon2_isEE) continue;
      }

      if ( Photon1_pt>_cutPhoton1Pt && Photon2_pt>_cutPhoton2Pt	

	   // EBEE Gap rejection
	   && (!Photon1_isEBEEGap) && (!Photon2_isEBEEGap)

	   && fabs(Photon1_detEta)<2.5 && fabs(Photon2_detEta)<2.5 

	   //	   && (Diphoton_Minv>1000.0)  
	   && (Diphoton_Minv>140.0)	   

	   // halo rejection
	   //	   && fabs(Diphoton_deltaPhi)>0.05
	   
	   //	   && (Photon1_maxRecHitTime>-10)&&(Photon2_maxRecHitTime>-10)
	   
	   ) {

	Double_t weightFake = 1.0;
	Double_t weightKFactor = 1.0;
	//	if (_Kfactor) weightKFactor=fKfactor->Eval(Diphoton_Minv);
	if (_Kfactor) {
	  weightKFactor=kfactorF(Diphoton_Minv);
	  //	  cout << "blah " << 	  weightKFactor << endl;
	  weightKFactor=weightKFactor*_effScaleFactor*_effScaleFactor;
	}

	nPass++;	  
	  
	//	cout << _fakeStatus << " " << nPass << " " << Diphoton_Minv << " " << Photon1_isEB << " " << Photon2_isEB << endl;

	//	if (!Photon1_isFakeable && !Photon2_isFakeable) {
	if (_fakeStatus=="TightTight") {
	  
	   // Tight ID
	  if (! ( 
		 (Photon1_ecalIso04 < 4.2+0.006*Photon1_pt) && (Photon2_ecalIso04 < 4.2+0.006*Photon2_pt) 
		 && (Photon1_hcalIso04 < 2.2+0.0025*Photon1_pt) && (Photon2_hcalIso04 < 2.2+0.0025*Photon2_pt) 
		 && (Photon1_hadOverEm<0.05) && (Photon2_hadOverEm<0.05)
		 && (Photon1_trkIsoSumPtHollow04< 2+0.001*Photon1_pt) && (Photon2_trkIsoSumPtHollow04< 2+0.001*Photon2_pt)
		 && ( (Photon1_isEB && Photon1_sigmaIetaIeta<0.013) || (Photon1_isEE && Photon1_sigmaIetaIeta<0.030) )
		 && ( (Photon2_isEB && Photon2_sigmaIetaIeta<0.013) || (Photon2_isEE && Photon2_sigmaIetaIeta<0.030) )
		 && (!Photon1_hasPixelSeed) && (!Photon2_hasPixelSeed) 
		  )) continue;	  

	  //	  std::cout  << Event_run << ":" << Event_evnum << ":" << Event_LS << ":" << Diphoton_Minv << ":" << Photon1_pt << ":" << Photon2_pt << ":" << Photon1_detEta << ":" << Photon2_detEta << std::endl;

	  //	  std::cout << "Toyoko : " << Event_run << "\t" << Event_evnum << "\t" << Event_LS// << "\t" << Diphoton_Minv << "\t" << Photon1_pt << "\t " << Photon2_pt << "\t" << Photon1_detEta << "\t" << Photon2_detEta << "\t" << Photon1_isEB << "\t" << Photon2_isEB << std::endl;

// 	  if (_filterGen) {
// 	    if (Diphoton_Minv>1000.0) continue;
// 	  }

//           	  if (Diphoton_Minv>800.0) {
//         	  std::cout << "High Mass Event: " << std::endl;
//          	  std::cout << " Run, LS, Event# : " << Event_run << "  " << Event_LS << " " << Event_evnum  << std::endl;
//          	  std::cout << " Diphoton minv, pt, cosTheta*, Deta ,Dphi :" << Diphoton_Minv  << " " << Diphoton_qt << " " << Diphoton_cosThetaStar << " " << Diphoton_deltaEta << " " << Diphoton_deltaPhi << std::endl;
//          	  std::cout << " Photon1, 2 isEB : " << Photon1_isEB << " " << Photon2_isEB << endl;
//          	  std::cout << " Photon1, 2 pt : " << Photon1_pt << " " << Photon2_pt << endl;
//          	  std::cout << " Photon1, 2 eta : " << Photon1_eta << " " << Photon2_eta << endl;
//          	  std::cout << " Photon1, 2 phi : " << Photon1_phi << " " << Photon2_phi << endl;
//          	  std::cout << " Photon1, 2 ecalIso : " << Photon1_ecalIso04 << " " << Photon2_ecalIso04 << endl;
//          	  std::cout << " Photon1, 2 hcalIso : " << Photon1_hcalIso04 << " " << Photon2_hcalIso04 << endl;
//          	  std::cout << " Photon1, 2 trackIso : " << Photon1_trkIsoSumPtHollow04 << " " << Photon2_trkIsoSumPtHollow04 << endl;
//          	  std::cout << " Photon1, 2 R9 : " << Photon1_r9 << " " << Photon2_r9 << endl;
//         	  std::cout << " Photon1, 2 maxRecHitTime : " << Photon1_maxRecHitTime << " " << Photon2_maxRecHitTime << endl;
//         	  std::cout << " Photon1, 2 maxRecHitTime : " << Photon1_maxRecHitTime << " " << Photon2_maxRecHitTime << endl;

//       //  	    //	    std::cout << " Photon1, 2 swisscross : " << Photon1_swisscross << " " << Photon2_swisscross << endl;
//       //  	    //	    std::cout << " Photon1, 2 e2/e9 : " << Photon1_e2e9 << " " << Photon2_e2e9 << endl;

//        	  }

// 	// MC Truth matching
	
// 	// status 3 = from hard scattering
// 	// status 1 = final state photon
	
// 	std::cout << Event_run << "\t" << Event_evnum << "\t" << Event_LS << endl;
	
// 	cout << "Photon1 " << Photon1_pt << " " <<   Photon1_eta << " "   << Photon1_phi << endl;
// 	cout << "Photon2 " << Photon2_pt << " " <<   Photon2_eta << " "   << Photon2_phi << endl;      
	
// 	cout << "MCMatchPhoton1_Status3 " << MCMatchPhoton1_Status3_status << " " <<     MCMatchPhoton1_Status3_PdgId << " " <<   MCMatchPhoton1_Status3_MotherPdgId << " " <<   MCMatchPhoton1_Status3_GrandmotherPdgId << " " <<    MCMatchPhoton1_Status3_pt << " " <<   MCMatchPhoton1_Status3_eta << " "   << MCMatchPhoton1_Status3_phi << endl;
// 	cout << "MCMatchPhoton2_Status3 " << MCMatchPhoton2_Status3_status << " " <<     MCMatchPhoton2_Status3_PdgId << " " <<   MCMatchPhoton2_Status3_MotherPdgId << " " <<   MCMatchPhoton2_Status3_GrandmotherPdgId << " " <<    MCMatchPhoton2_Status3_pt << " " <<   MCMatchPhoton2_Status3_eta << " "   << MCMatchPhoton2_Status3_phi << endl;
// 	cout << "MCMatchPhoton1_Status1 " << MCMatchPhoton1_Status1_status << " " <<     MCMatchPhoton1_Status1_PdgId << " " <<   MCMatchPhoton1_Status1_MotherPdgId << " " <<   MCMatchPhoton1_Status1_GrandmotherPdgId << " " <<    MCMatchPhoton1_Status1_pt << " " <<   MCMatchPhoton1_Status1_eta << " "   << MCMatchPhoton1_Status1_phi << endl;
// 	cout << "MCMatchPhoton2_Status1 " << MCMatchPhoton2_Status1_status << " " <<     MCMatchPhoton2_Status1_PdgId << " " <<   MCMatchPhoton2_Status1_MotherPdgId << " " <<   MCMatchPhoton2_Status1_GrandmotherPdgId << " " <<  MCMatchPhoton2_Status1_pt << " " <<   MCMatchPhoton2_Status1_eta << " "   << MCMatchPhoton2_Status1_phi << endl;


// 	bool mcMatchPhoton1Status3 = kFALSE;
// 	bool mcMatchPhoton2Status3 = kFALSE;
	
// 	if (fabs(MCMatchPhoton1_Status3_phi)<100&&fabs(MCMatchPhoton2_Status3_phi<100)) {

// 	  double dPhi1 = deltaPhi(Photon1_phi,MCMatchPhoton1_Status3_phi);
// 	  double dEta1 = Photon1_eta-MCMatchPhoton1_Status3_eta;
// 	  double dR1 = TMath::Sqrt(dPhi1*dPhi1 + dEta1*dEta1);
	  
// 	  double dPhi2 = deltaPhi(Photon2_phi,MCMatchPhoton2_Status3_phi);
// 	  double dEta2 = Photon2_eta-MCMatchPhoton2_Status3_eta;
// 	  double dR2 = TMath::Sqrt(dPhi2*dPhi2 + dEta2*dEta2);
	  
// 	  //      cout << "dR1 dR2 " <<  dR1 << " " << dR2 << endl;
	  
// 	  if (dR1<0.1&&MCMatchPhoton1_Status3_PdgId==22) { mcMatchPhoton1Status3=kTRUE; }
// 	  if (dR2<0.1&&MCMatchPhoton2_Status3_PdgId==22) { mcMatchPhoton2Status3=kTRUE; }
// 	  h_mcMatch1->Fill(dR1);
// 	  h_mcMatch2->Fill(dR2);
// 	  if (mcMatchPhoton1Status3) cout << "Photon1 Status3 TRUTH MATCHED!! " << endl;
// 	  if (mcMatchPhoton2Status3) cout << "Photon2 Status3 TRUTH MATCHED!! " << endl;

// 	  if (!mcMatchPhoton1Status3) cout << "Photon1 Status3 NOT TRUTH MATCHED!! " << endl;
// 	  if (!mcMatchPhoton2Status3) cout << "Photon2 Status3 NOT TRUTH MATCHED!! " << endl;

// 	  if (mcMatchPhoton1Status3&&mcMatchPhoton2Status3) cout << "Both Status3 TRUTH MATCHED!! " << endl;
// 	}

// 	bool mcMatchPhoton1Status1 = kFALSE;
// 	bool mcMatchPhoton2Status1 = kFALSE;

// 	if (fabs(MCMatchPhoton1_Status1_phi)<100&&fabs(MCMatchPhoton2_Status1_phi<100)) {

// 	  double dPhi1 = deltaPhi(Photon1_phi,MCMatchPhoton1_Status1_phi);
// 	  double dEta1 = Photon1_eta-MCMatchPhoton1_Status1_eta;
// 	  double dR1 = TMath::Sqrt(dPhi1*dPhi1 + dEta1*dEta1);
	  
// 	  double dPhi2 = deltaPhi(Photon2_phi,MCMatchPhoton2_Status1_phi);
// 	  double dEta2 = Photon2_eta-MCMatchPhoton2_Status1_eta;
// 	  double dR2 = TMath::Sqrt(dPhi2*dPhi2 + dEta2*dEta2);
	  
// 	  //      cout << "dR1 dR2 " <<  dR1 << " " << dR2 << endl;
	  
// 	  if (dR1<0.1&&MCMatchPhoton1_Status1_PdgId==22) { mcMatchPhoton1Status1=kTRUE; }
// 	  if (dR2<0.1&&MCMatchPhoton2_Status1_PdgId==22) { mcMatchPhoton2Status1=kTRUE; }
// 	  h_mcMatch1->Fill(dR1);
// 	  h_mcMatch2->Fill(dR2);
// 	  if (mcMatchPhoton1Status1) cout << "Photon1 Status1 TRUTH MATCHED!! " << endl;
// 	  if (mcMatchPhoton2Status1) cout << "Photon2 Status1 TRUTH MATCHED!! " << endl;

// 	  if (!mcMatchPhoton1Status1) cout << "Photon1 Status1 NOT TRUTH MATCHED!! " << endl;
// 	  if (!mcMatchPhoton2Status1) cout << "Photon2 Status1 NOT TRUTH MATCHED!! " << endl;

// 	  if (mcMatchPhoton1Status1&&mcMatchPhoton2Status1) cout << "Both Status1 TRUTH MATCHED!! " << endl;
// 	}
// 	cout << "Signal Process ID " << GenEvent_signalProcessId << " " << mcMatchPhoton1Status3 << " " << mcMatchPhoton1Status1 << " " << mcMatchPhoton2Status3 << " " << mcMatchPhoton2Status1 << endl;

// 	//	if (mcMatchPhoton1Status1&&mcMatchPhoton1Status3&&!mcMatchPhoton2Status3&&mcMatchPhoton2Status1) {

// 	//	}

// 	// Photon1 status3-matched
// 	// if ( (mcMatchPhoton1!=kTRUE)||(mcMatchPhoton2==kTRUE) ) continue;
// 	// Photon2 status3-matched
// 	//	if ( (mcMatchPhoton2!=kTRUE)||(mcMatchPhoton1==kTRUE) ) continue;

// 	//	if ( (mcMatchPhoton1==kTRUE&&mcMatchPhoton2==kTRUE) || (mcMatchPhoton1==kFALSE&&mcMatchPhoton2==kFALSE)) continue;

	 	  
	  //	if (TrigHLT_HLT_MinBiasBSC>0) h_TrigHLT->Fill(0);
	  //	if (TrigHLT_HLT_MinBiasBSC_NoBPTX>0) h_TrigHLT->Fill(1);
	  //	if (TrigHLT_HLT_MinBiasBSC_OR>0) h_TrigHLT->Fill(2);
	  //	if (TrigHLT_HLT_L1_BscMinBiasOR_BptxPlusORMinus>0) h_TrigHLT->Fill(3);
	  if (TrigHLT_HLT_L1SingleEG2>0) h_TrigHLT->Fill(0.0,puWeight*weightKFactor);
	  if (TrigHLT_HLT_L1SingleEG5>0) h_TrigHLT->Fill(1,puWeight*weightKFactor);
	  if (TrigHLT_HLT_L1SingleEG8>0) h_TrigHLT->Fill(2,puWeight*weightKFactor);
	  if (TrigHLT_HLT_L1DoubleEG5>0) h_TrigHLT->Fill(3,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon10_L1R>0) h_TrigHLT->Fill(4,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon10_Cleaned_L1R>0) h_TrigHLT->Fill(5,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon15_L1R>0) h_TrigHLT->Fill(6,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon15_Cleaned_L1R>0) h_TrigHLT->Fill(7,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon15_LooseEcalIso_L1R>0) h_TrigHLT->Fill(8,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon15_LooseEcalIso_Cleaned_L1R>0) h_TrigHLT->Fill(9,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon15_TrackIso_L1R>0) h_TrigHLT->Fill(10,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon15_TrackIso_Cleaned_L1R>0) h_TrigHLT->Fill(11,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon17_Isol_SC17HE_L1R_v1>0) h_TrigHLT->Fill(12,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon17_SC17HE_L1R_v1>0) h_TrigHLT->Fill(13,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon20_L1R>0) h_TrigHLT->Fill(14,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon20_Cleaned_L1R>0) h_TrigHLT->Fill(15,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon20_NoHE_L1R>0) h_TrigHLT->Fill(16,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon22_SC22HE_L1R_v1>0) h_TrigHLT->Fill(17,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon25_Cleaned_L1R>0) h_TrigHLT->Fill(18,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon30_L1R>0) h_TrigHLT->Fill(19,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon30_Cleaned_L1R>0) h_TrigHLT->Fill(20,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon30_L1R_8E29>0) h_TrigHLT->Fill(21,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1>0) h_TrigHLT->Fill(22,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon35_Isol_Cleaned_L1R_v1>0) h_TrigHLT->Fill(23,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon40_CaloId_Cleaned_L1R_v1>0) h_TrigHLT->Fill(24,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon40_Isol_Cleaned_L1R_v1>0) h_TrigHLT->Fill(25,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon50_L1R>0) h_TrigHLT->Fill(26,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon50_Cleaned_L1R>0) h_TrigHLT->Fill(27,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon50_Cleaned_L1R_v1>0) h_TrigHLT->Fill(28,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon50_NoHE_L1R>0) h_TrigHLT->Fill(29,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon50_NoHE_Cleaned_L1R>0) h_TrigHLT->Fill(30,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon70_Cleaned_L1R_v1>0) h_TrigHLT->Fill(31,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon70_NoHE_Cleaned_L1R_v1>0) h_TrigHLT->Fill(32,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon100_NoHE_Cleaned_L1R_v1>0) h_TrigHLT->Fill(33,puWeight*weightKFactor);
	  if (TrigHLT_HLT_Photon110_NoHE_Cleaned_L1R_v1>0) h_TrigHLT->Fill(34,puWeight*weightKFactor);
	  if (TrigHLT_HLT_DoublePhoton5_L1R>0) h_TrigHLT->Fill(35,puWeight*weightKFactor);
	  if (TrigHLT_HLT_DoublePhoton5_CEP_L1R>0) h_TrigHLT->Fill(36,puWeight*weightKFactor);
	  if (TrigHLT_HLT_DoublePhoton5_CEP_L1R_v3>0) h_TrigHLT->Fill(37,puWeight*weightKFactor);
	  if (TrigHLT_HLT_DoublePhoton5_Jpsi_L1R>0) h_TrigHLT->Fill(38,puWeight*weightKFactor);
	  if (TrigHLT_HLT_DoublePhoton5_Upsilon_L1R>0) h_TrigHLT->Fill(39,puWeight*weightKFactor);
	  if (TrigHLT_HLT_DoublePhoton10_L1R>0) h_TrigHLT->Fill(40,puWeight*weightKFactor);
	  if (TrigHLT_HLT_DoublePhoton15_L1R>0) h_TrigHLT->Fill(41,puWeight*weightKFactor);
	  if (TrigHLT_HLT_DoublePhoton17_L1R>0) h_TrigHLT->Fill(42,puWeight*weightKFactor);
	  if (TrigHLT_HLT_DoublePhoton17_SingleIsol_L1R_v1>0) h_TrigHLT->Fill(43,puWeight*weightKFactor);
	  if (TrigHLT_HLT_DoublePhoton20_L1R>0) h_TrigHLT->Fill(44,puWeight*weightKFactor);
	  if (TrigHLT_HLT_DoublePhoton22_L1R_v1>0) h_TrigHLT->Fill(45,puWeight*weightKFactor);

	  h_Nvtx->Fill(Vtx_Nvtx,puWeight*weightKFactor);
	  
	  h_Photon1_pt->Fill(Photon1_pt,puWeight*weightKFactor);
	  h_Photon1_pt_log->Fill(Photon1_pt,puWeight*weightKFactor);
	  h_Photon1_pt_zoom->Fill(Photon1_pt,puWeight*weightKFactor);
	  h_Photon1_eta->Fill(Photon1_eta,puWeight*weightKFactor);
	  h_Photon1_phi->Fill(Photon1_phi,puWeight*weightKFactor);
	  h_Photon1_occupancy->Fill(Photon1_eta, Photon1_phi,puWeight*weightKFactor);
	  
	  h_Photon1_r9->Fill(Photon1_r9,puWeight*weightKFactor);
	  h_Photon1_sigmaIetaIeta->Fill(Photon1_sigmaIetaIeta,puWeight*weightKFactor);
	  h_Photon1_sigmaEtaEta->Fill(Photon1_sigmaEtaEta,puWeight*weightKFactor);
	  
	  if (Photon1_swisscross>-100) h_Photon1_swisscross->Fill(Photon1_swisscross,puWeight*weightKFactor);
	  if (Photon1_e2e9>0) h_Photon1_e2e9->Fill(Photon1_e2e9,puWeight*weightKFactor);
	  if (Photon1_e4x4!=0 && Photon1_e4x4>-100 && Photon1_e2x2>-100) h_Photon1_e2x2e4x4->Fill(Photon1_e2x2/Photon1_e4x4,puWeight*weightKFactor);
	  h_Photon1_severityLevel->Fill(Photon1_severityLevel,puWeight*weightKFactor);
	  h_Photon1_recHitFlag->Fill(Photon1_recHitFlag,puWeight*weightKFactor);
	  h_Photon1_maxRecHitTime->Fill(Photon1_maxRecHitTime,puWeight*weightKFactor);
	  h_Photon1_maxRecHitTime_wide->Fill(Photon1_maxRecHitTime,puWeight*weightKFactor);
	  h_Photon1_e2e9_v_maxRecHitTime->Fill(Photon1_maxRecHitTime,Photon1_e2e9,puWeight*weightKFactor);
	  	  
	  h_Photon1_hadOverEm->Fill(Photon1_hadOverEm,puWeight*weightKFactor);
	  h_Photon1_hcalIso04->Fill(Photon1_hcalIso04,puWeight*weightKFactor);
	  h_Photon1_hcalIso03->Fill(Photon1_hcalIso03,puWeight*weightKFactor);
	  h_Photon1_ecalIso04->Fill(Photon1_ecalIso04,puWeight*weightKFactor);
	  h_Photon1_ecalIso03->Fill(Photon1_ecalIso03,puWeight*weightKFactor);
	  
	  h_Photon1_trkIsoSumPtHollow04->Fill(Photon1_trkIsoSumPtHollow04,puWeight*weightKFactor);
	  h_Photon1_trkIsoSumPtSolid04->Fill(Photon1_trkIsoSumPtSolid04,puWeight*weightKFactor);
	  h_Photon1_trkIsoSumPtHollow03->Fill(Photon1_trkIsoSumPtHollow03,puWeight*weightKFactor);
	  h_Photon1_trkIsoSumPtSolid03->Fill(Photon1_trkIsoSumPtSolid03,puWeight*weightKFactor);
	  
	  h_Photon1_trkIsoNtrksHollow04->Fill(Photon1_trkIsoNtrksHollow04,puWeight*weightKFactor);
	  h_Photon1_trkIsoNtrksSolid04->Fill(Photon1_trkIsoNtrksSolid04,puWeight*weightKFactor);
	  h_Photon1_trkIsoNtrksHollow03->Fill(Photon1_trkIsoNtrksHollow03,puWeight*weightKFactor);
	  h_Photon1_trkIsoNtrksSolid03->Fill(Photon1_trkIsoNtrksSolid03,puWeight*weightKFactor);
	  
	  // 	h_Photon1_esRatio->Fill(Photon1_esRatio,puWeight*weightKFactor);
	  
	  h_Photon2_pt->Fill(Photon2_pt,puWeight*weightKFactor);
	  h_Photon2_pt_log->Fill(Photon2_pt,puWeight*weightKFactor);
	  h_Photon2_pt_zoom->Fill(Photon2_pt,puWeight*weightKFactor);
	  h_Photon2_eta->Fill(Photon2_eta,puWeight*weightKFactor);
	  h_Photon2_phi->Fill(Photon2_phi,puWeight*weightKFactor);
	  h_Photon2_occupancy->Fill(Photon2_eta, Photon2_phi,puWeight*weightKFactor);

	  h_Photon2_r9->Fill(Photon2_r9,puWeight*weightKFactor);
	  h_Photon2_sigmaIetaIeta->Fill(Photon2_sigmaIetaIeta,puWeight*weightKFactor);
	  h_Photon2_sigmaEtaEta->Fill(Photon2_sigmaEtaEta,puWeight*weightKFactor);
	  
	  if (Photon2_swisscross>-100) h_Photon2_swisscross->Fill(Photon2_swisscross,puWeight*weightKFactor);
	  if (Photon2_e2e9>0) h_Photon2_e2e9->Fill(Photon2_e2e9,puWeight*weightKFactor);
	  if (Photon2_e4x4!=0 && Photon2_e4x4>-100 && Photon2_e2x2>-100) h_Photon2_e2x2e4x4->Fill(Photon2_e2x2/Photon2_e4x4,puWeight*weightKFactor);
	  h_Photon2_severityLevel->Fill(Photon2_severityLevel,puWeight*weightKFactor);
	  h_Photon2_recHitFlag->Fill(Photon2_recHitFlag,puWeight*weightKFactor);
	  h_Photon2_maxRecHitTime->Fill(Photon2_maxRecHitTime,puWeight*weightKFactor);
	  h_Photon2_maxRecHitTime_wide->Fill(Photon2_maxRecHitTime,puWeight*weightKFactor);
	  h_Photon2_e2e9_v_maxRecHitTime->Fill(Photon2_maxRecHitTime,Photon1_e2e9,puWeight*weightKFactor);
	  
	  h_Photon2_hadOverEm->Fill(Photon2_hadOverEm,puWeight*weightKFactor);
	  h_Photon2_hcalIso04->Fill(Photon2_hcalIso04,puWeight*weightKFactor);
	  h_Photon2_hcalIso03->Fill(Photon2_hcalIso03,puWeight*weightKFactor);
	  h_Photon2_ecalIso04->Fill(Photon2_ecalIso04,puWeight*weightKFactor);
	  h_Photon2_ecalIso03->Fill(Photon2_ecalIso03,puWeight*weightKFactor);
	  
	  h_Photon2_trkIsoSumPtHollow04->Fill(Photon2_trkIsoSumPtHollow04,puWeight*weightKFactor);
	  h_Photon2_trkIsoSumPtSolid04->Fill(Photon2_trkIsoSumPtSolid04,puWeight*weightKFactor);
	  h_Photon2_trkIsoSumPtHollow03->Fill(Photon2_trkIsoSumPtHollow03,puWeight*weightKFactor);
	  h_Photon2_trkIsoSumPtSolid03->Fill(Photon2_trkIsoSumPtSolid03,puWeight*weightKFactor);
	  
	  h_Photon2_trkIsoNtrksHollow04->Fill(Photon2_trkIsoNtrksHollow04,puWeight*weightKFactor);
	  h_Photon2_trkIsoNtrksSolid04->Fill(Photon2_trkIsoNtrksSolid04,puWeight*weightKFactor);
	  h_Photon2_trkIsoNtrksHollow03->Fill(Photon2_trkIsoNtrksHollow03,puWeight*weightKFactor);
	  h_Photon2_trkIsoNtrksSolid03->Fill(Photon2_trkIsoNtrksSolid03,puWeight*weightKFactor);
	  
	  // 	h_Photon2_esRatio->Fill(Photon2_esRatio,puWeight*weightKFactor);
	  
	  h_Diphoton_Minv->Fill(Diphoton_Minv,puWeight*weightKFactor);

	  h_Diphoton_Minv_Yousi5->Fill(Diphoton_Minv,puWeight*weightKFactor);

 	  if (_filterGen&&Diphoton_Minv>1000.0) {
 	  } else {	    
	    h_Diphoton_Minv_Yousi10->Fill(Diphoton_Minv,puWeight*weightKFactor);
	    h_Diphoton_Minv_Yousi40->Fill(Diphoton_Minv,puWeight*weightKFactor);
	    h_Diphoton_Minv_Yousi100->Fill(Diphoton_Minv,puWeight*weightKFactor);
	  }

	  if (Diphoton_Minv>120&&Diphoton_Minv<200) h_Diphoton_Minv_120to200->Fill(Diphoton_Minv,puWeight*weightKFactor);
	  if (Diphoton_Minv>200&&Diphoton_Minv<500) h_Diphoton_Minv_200to500->Fill(Diphoton_Minv,puWeight*weightKFactor);
	  if (Diphoton_Minv>500&&Diphoton_Minv<800) h_Diphoton_Minv_500to800->Fill(Diphoton_Minv,puWeight*weightKFactor);
	  if (Diphoton_Minv>800) h_Diphoton_Minv_800toInf->Fill(Diphoton_Minv,puWeight*weightKFactor);

	  h_Diphoton_Minv_log->Fill(Diphoton_Minv,puWeight*weightKFactor);
	  //	  h_Diphoton_Minv_add->Fill(Diphoton_Minv,puWeight*weightKFactor);
	  h_Diphoton_Minv_v_Photon1_pt->Fill(Diphoton_Minv,Photon1_pt,puWeight*weightKFactor);
	  h_Diphoton_Minv_v_Photon2_pt->Fill(Diphoton_Minv,Photon2_pt,puWeight*weightKFactor);	 

	  if (Diphoton_Minv>500) h_Diphoton_Minv_high->Fill(Diphoton_Minv,puWeight*weightKFactor);
	  if (Diphoton_Minv<500) h_Diphoton_Minv_low->Fill(Diphoton_Minv,puWeight*weightKFactor);
	  if (Diphoton_Minv<500) h_Diphoton_Minv_low_bin1GeV->Fill(Diphoton_Minv,puWeight*weightKFactor);

	  h_Diphoton_qt->Fill(Diphoton_qt,puWeight*weightKFactor);
	  h_Diphoton_deltaPhi->Fill(Diphoton_deltaPhi,puWeight*weightKFactor);
	  h_Diphoton_deltaEta->Fill(Diphoton_deltaEta,puWeight*weightKFactor);
	  h_Diphoton_deltaR->Fill(Diphoton_deltaR,puWeight*weightKFactor);
	  h_Diphoton_cosThetaStar->Fill(Diphoton_cosThetaStar,puWeight*weightKFactor);
	  
	  // Now for fake rate stuff
	  
	  // then this is signal sample = tight-tight
	  weightFake = 1;
	  // fill all histograms with appropriate weight for this event
	  h_FakeRate_tt_pt1->Fill(Photon1_pt,weightFake);
	  h_FakeRate_tt_pt1_zoom->Fill(Photon1_pt,weightFake);
	  h_FakeRate_tt_eta1->Fill(Photon1_eta,weightFake);
	  h_FakeRate_tt_phi1->Fill(Photon1_phi,weightFake);
	  h_FakeRate_tt_pt2->Fill(Photon2_pt,weightFake);
	  h_FakeRate_tt_pt2_zoom->Fill(Photon2_pt,weightFake);
	  h_FakeRate_tt_eta2->Fill(Photon2_eta,weightFake);	  
	  h_FakeRate_tt_phi2->Fill(Photon2_phi,weightFake);	  
	  h_FakeRate_tt_minv->Fill(Diphoton_Minv,weightFake);	  
	  if (Diphoton_Minv>200) h_FakeRate_tt_minv_high->Fill(Diphoton_Minv,weightFake);	  
	  h_FakeRate_tt_qt->Fill(Diphoton_qt,weightFake);
	  h_FakeRate_tt_deltaPhi->Fill(Diphoton_deltaPhi,weightFake);
	  h_FakeRate_tt_deltaEta->Fill(Diphoton_deltaEta,weightFake);
	  h_FakeRate_tt_deltaR->Fill(Diphoton_deltaR,weightFake);
	  h_FakeRate_tt_cosThetaStar->Fill(Diphoton_cosThetaStar,weightFake);

	  h_trkIsoVRho1->Fill(Photon1_trkIsoSumPtHollow04,rho);
	  h_ecalIsoVRho1->Fill(Photon1_ecalIso04,rho);
	  h_hcalIsoVRho1->Fill(Photon1_hcalIso04,rho);
	  h_trkIsoVRho2->Fill(Photon2_trkIsoSumPtHollow04,rho);
	  h_ecalIsoVRho2->Fill(Photon2_ecalIso04,rho);
	  h_hcalIsoVRho2->Fill(Photon2_hcalIso04,rho);

	  if (Diphoton_Minv>120&&Diphoton_Minv<200) h_FakeRate_tt_minv_120to200->Fill(Diphoton_Minv,weightFake);
	  if (Diphoton_Minv>200&&Diphoton_Minv<500) h_FakeRate_tt_minv_200to500->Fill(Diphoton_Minv,weightFake);
	  if (Diphoton_Minv>500&&Diphoton_Minv<800) h_FakeRate_tt_minv_500to800->Fill(Diphoton_Minv,weightFake);
	  if (Diphoton_Minv>800) h_FakeRate_tt_minv_800toInf->Fill(Diphoton_Minv,weightFake);
	  
// 	  if (_outT0!=0) {    
// 	    outVar0.weight = _weight;
// 	    outVar0.minv = Diphoton_Minv;
// 	    outVar0.r9_photon1 = Photon1_r9;
// 	    outVar0.r9_photon2 = Photon2_r9;
// 	    outVar0.eta_photon1 = Photon1_detEta;
// 	    outVar0.eta_photon2 = Photon2_detEta;
// 	    _outT0->Fill();
// 	  }       
	  
	  //	} else if (Photon1_isFakeable && !Photon2_isFakeable) {
	} else if (_fakeStatus=="FakeTight") {	  

	  // this will be the appropriate weight for this event
	  // depending one whichever of Photon1/Photon2 (or both) is Fakeable
	  
	  // obtain a weight for this event from the fake rate function
	  // for this photon pt	
	  if ( Photon1_isEB ) {
	    weightFake = fake_rate_fn->Eval(Photon1_pt);
	  } else {
	    weightFake = fake_rate_fn_EE->Eval(Photon1_pt);
	  }

	  // fill all histograms with appropriate weight for this event
	  h_FakeRate_tf_pt1->Fill(Photon1_pt,weightFake);
	  h_FakeRate_tf_pt1_zoom->Fill(Photon1_pt,weightFake);
	  h_FakeRate_tf_eta1->Fill(Photon1_eta,weightFake);
	  h_FakeRate_tf_phi1->Fill(Photon1_phi,weightFake);
	  h_FakeRate_tf_pt2->Fill(Photon2_pt,weightFake);
	  h_FakeRate_tf_pt2_zoom->Fill(Photon2_pt,weightFake);
	  h_FakeRate_tf_eta2->Fill(Photon2_eta,weightFake);	  
	  h_FakeRate_tf_phi2->Fill(Photon2_phi,weightFake);	  
	  h_FakeRate_tf_minv->Fill(Diphoton_Minv,weightFake);
	  if (Diphoton_Minv>200) h_FakeRate_tf_minv_high->Fill(Diphoton_Minv,weightFake);	  
	  h_FakeRate_tf_qt->Fill(Diphoton_qt,weightFake);
	  h_FakeRate_tf_deltaPhi->Fill(Diphoton_deltaPhi,weightFake);
	  h_FakeRate_tf_deltaEta->Fill(Diphoton_deltaEta,weightFake);
	  h_FakeRate_tf_deltaR->Fill(Diphoton_deltaR,weightFake);
	  h_FakeRate_tf_cosThetaStar->Fill(Diphoton_cosThetaStar,weightFake);

	  h_FakeRate_ft_pt1_noweight->Fill(Photon1_pt);
	  h_FakeRate_ft_pt2_noweight->Fill(Photon2_pt);

	  if (Diphoton_Minv>120&&Diphoton_Minv<200) h_FakeRate_tf_minv_120to200->Fill(Diphoton_Minv,weightFake);
	  if (Diphoton_Minv>200&&Diphoton_Minv<500) h_FakeRate_tf_minv_200to500->Fill(Diphoton_Minv,weightFake);
	  if (Diphoton_Minv>500&&Diphoton_Minv<800) h_FakeRate_tf_minv_500to800->Fill(Diphoton_Minv,weightFake);
	  if (Diphoton_Minv>800) h_FakeRate_tf_minv_800toInf->Fill(Diphoton_Minv,weightFake);

	  //	} else if (!Photon1_isFakeable && Photon2_isFakeable) {	
	} else if (_fakeStatus=="TightFake") {	  

	  // tight-fake, but other way round
	  if ( Photon2_isEB ) {
	    weightFake = fake_rate_fn->Eval(Photon2_pt);
	  } else {
	    weightFake = fake_rate_fn_EE->Eval(Photon2_pt);
	  }

	  // fill all histograms with appropriate weight for this event
	  h_FakeRate_tf_pt1->Fill(Photon1_pt,weightFake);
	  h_FakeRate_tf_pt1_zoom->Fill(Photon1_pt,weightFake);
	  h_FakeRate_tf_eta1->Fill(Photon1_eta,weightFake);
	  h_FakeRate_tf_phi1->Fill(Photon1_phi,weightFake);
	  h_FakeRate_tf_pt2->Fill(Photon2_pt,weightFake);
	  h_FakeRate_tf_pt2_zoom->Fill(Photon2_pt,weightFake);
	  h_FakeRate_tf_eta2->Fill(Photon2_eta,weightFake);	  
	  h_FakeRate_tf_phi2->Fill(Photon2_phi,weightFake);	  
	  h_FakeRate_tf_minv->Fill(Diphoton_Minv,weightFake);
	  if (Diphoton_Minv>200) h_FakeRate_tf_minv_high->Fill(Diphoton_Minv,weightFake);	  
	  h_FakeRate_tf_qt->Fill(Diphoton_qt,weightFake);
	  h_FakeRate_tf_deltaPhi->Fill(Diphoton_deltaPhi,weightFake);
	  h_FakeRate_tf_deltaEta->Fill(Diphoton_deltaEta,weightFake);
	  h_FakeRate_tf_deltaR->Fill(Diphoton_deltaR,weightFake);
	  h_FakeRate_tf_cosThetaStar->Fill(Diphoton_cosThetaStar,weightFake);

	  h_FakeRate_tf_pt1_noweight->Fill(Photon1_pt);
	  h_FakeRate_tf_pt2_noweight->Fill(Photon2_pt);

	  if (Diphoton_Minv>120&&Diphoton_Minv<200) h_FakeRate_tf_minv_120to200->Fill(Diphoton_Minv,weightFake);
	  if (Diphoton_Minv>200&&Diphoton_Minv<500) h_FakeRate_tf_minv_200to500->Fill(Diphoton_Minv,weightFake);
	  if (Diphoton_Minv>500&&Diphoton_Minv<800) h_FakeRate_tf_minv_500to800->Fill(Diphoton_Minv,weightFake);
	  if (Diphoton_Minv>800) h_FakeRate_tf_minv_800toInf->Fill(Diphoton_Minv,weightFake);

	  //	} else if (Photon1_isFakeable && Photon2_isFakeable) {

	} else if (_fakeStatus=="FakeFake") {	  

	  // then both are fake!

	  double tempPhoton1 = 1;
	  double tempPhoton2 = 1;

	  if ( Photon2_isEB ) {
	    tempPhoton2 = fake_rate_fn->Eval(Photon2_pt);
	  } else {
	    tempPhoton2 = fake_rate_fn_EE->Eval(Photon2_pt);
	  }
	  if ( Photon1_isEB ) {
	    tempPhoton1 = fake_rate_fn->Eval(Photon1_pt);
	  } else {
	    tempPhoton1 = fake_rate_fn_EE->Eval(Photon1_pt);
	  }

	  weightFake = tempPhoton1*tempPhoton1;
	  //	  weightFake = fake_rate_fn->Eval(Photon1_pt)*fake_rate_fn->Eval(Photon2_pt);

	  // fill all histograms with appropriate weight for this event

	  h_FakeRate_ff_pt1->Fill(Photon1_pt,weightFake);
	  h_FakeRate_ff_pt1_zoom->Fill(Photon1_pt,weightFake);
	  h_FakeRate_ff_eta1->Fill(Photon1_eta,weightFake);
	  h_FakeRate_ff_phi1->Fill(Photon1_phi,weightFake);
	  h_FakeRate_ff_pt2->Fill(Photon2_pt,weightFake);
	  h_FakeRate_ff_pt2_zoom->Fill(Photon2_pt,weightFake);
	  h_FakeRate_ff_eta2->Fill(Photon2_eta,weightFake);	  
	  h_FakeRate_ff_phi2->Fill(Photon2_phi,weightFake);	  
	  h_FakeRate_ff_minv->Fill(Diphoton_Minv,weightFake);
	  if (Diphoton_Minv>200) h_FakeRate_ff_minv_high->Fill(Diphoton_Minv,weightFake);	  
	  h_FakeRate_ff_qt->Fill(Diphoton_qt,weightFake);
	  h_FakeRate_ff_deltaPhi->Fill(Diphoton_deltaPhi,weightFake);
	  h_FakeRate_ff_deltaEta->Fill(Diphoton_deltaEta,weightFake);
	  h_FakeRate_ff_deltaR->Fill(Diphoton_deltaR,weightFake);
	  h_FakeRate_ff_cosThetaStar->Fill(Diphoton_cosThetaStar,weightFake);

	  h_FakeRate_ff_pt1_noweight->Fill(Photon1_pt);
	  h_FakeRate_ff_pt2_noweight->Fill(Photon2_pt);

	  if (Diphoton_Minv>120&&Diphoton_Minv<200) h_FakeRate_ff_minv_120to200->Fill(Diphoton_Minv,weightFake);
	  if (Diphoton_Minv>200&&Diphoton_Minv<500) h_FakeRate_ff_minv_200to500->Fill(Diphoton_Minv,weightFake);
	  if (Diphoton_Minv>500&&Diphoton_Minv<800) h_FakeRate_ff_minv_500to800->Fill(Diphoton_Minv,weightFake);
	  if (Diphoton_Minv>800) h_FakeRate_ff_minv_800toInf->Fill(Diphoton_Minv,weightFake);

	} else {
	  cout << "You should not be in here..." << endl;
	}	

      } // cuts

   }

   if ( _outF != 0 ) {
     _outF->cd();
     h_mcMatch1->Write();
     h_mcMatch2->Write();
     h_Nvtx->Write();
     h_TrigHLT->Write();
     h_Diphoton_Minv->Write();
     h_Diphoton_Minv_Yousi5->Write();
     h_Diphoton_Minv_Yousi10->Write();
     h_Diphoton_Minv_Yousi40->Write();
     h_Diphoton_Minv_Yousi100->Write();
     h_Diphoton_Minv_120to200->Write();
     h_Diphoton_Minv_200to500->Write();
     h_Diphoton_Minv_500to800->Write();
     h_Diphoton_Minv_800toInf->Write();
     h_Diphoton_Minv_log->Write();
     //     h_Diphoton_Minv_add->Write();
     h_Diphoton_Minv_high->Write();
     h_Diphoton_Minv_low->Write();
     h_Diphoton_Minv_low_bin1GeV->Write();
     h_Diphoton_qt->Write();
     h_Diphoton_deltaPhi->Write();
     h_Diphoton_deltaEta->Write();
     h_Diphoton_deltaR->Write();
     h_Diphoton_cosThetaStar->Write();
     h_Diphoton_Minv_v_Photon1_pt->Write();
     h_Diphoton_Minv_v_Photon2_pt->Write();
     h_Photon1_pt->Write();
     h_Photon1_pt_log->Write();
     h_Photon1_pt_zoom->Write();
     h_Photon1_eta->Write();
     h_Photon1_phi->Write();
     h_Photon1_occupancy->Write();
     h_Photon1_r9->Write();
     h_Photon1_sigmaIetaIeta->Write();
     h_Photon1_sigmaEtaEta->Write();
     h_Photon1_e2e9->Write();
     h_Photon1_e2x2e4x4->Write();
     h_Photon1_swisscross->Write();
     h_Photon1_swisscross->Write();
     h_Photon1_severityLevel->Write();
     h_Photon1_recHitFlag->Write();
     h_Photon1_maxRecHitTime->Write();
     h_Photon1_maxRecHitTime_wide->Write();
     h_Photon1_e2e9_v_maxRecHitTime->Write();
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
     h_Photon2_pt_log->Write();
     h_Photon2_pt_zoom->Write();
     h_Photon2_eta->Write();
     h_Photon2_phi->Write();
     h_Photon2_occupancy->Write();
     h_Photon2_r9->Write();
     h_Photon2_sigmaIetaIeta->Write();
     h_Photon2_sigmaEtaEta->Write();
     h_Photon2_e2e9->Write();
     h_Photon2_e2x2e4x4->Write();
     h_Photon2_swisscross->Write();
     h_Photon2_severityLevel->Write();
     h_Photon2_recHitFlag->Write();
     h_Photon2_maxRecHitTime->Write();
     h_Photon2_maxRecHitTime_wide->Write();
     h_Photon2_e2e9_v_maxRecHitTime->Write();
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
     h_FakeRate_tt_pt1->Write();
     h_FakeRate_tt_pt1_zoom->Write();
     h_FakeRate_tt_eta1->Write();
     h_FakeRate_tt_phi1->Write();
     h_FakeRate_tt_pt2->Write();
     h_FakeRate_tt_pt2_zoom->Write();
     h_FakeRate_tt_eta2->Write();
     h_FakeRate_tt_phi2->Write();
     h_FakeRate_tt_minv->Write();   
     h_FakeRate_tt_minv_high->Write();   
     h_FakeRate_tt_qt->Write();
     h_FakeRate_tt_deltaPhi->Write();
     h_FakeRate_tt_deltaEta->Write();
     h_FakeRate_tt_deltaR->Write();
     h_FakeRate_tt_cosThetaStar->Write();
     h_FakeRate_tf_pt1->Write();
     h_FakeRate_tf_pt1_zoom->Write();
     h_FakeRate_tf_eta1->Write();
     h_FakeRate_tf_phi1->Write();
     h_FakeRate_tf_pt2->Write();
     h_FakeRate_tf_pt2_zoom->Write();
     h_FakeRate_tf_eta2->Write();
     h_FakeRate_tf_phi2->Write();
     h_FakeRate_tf_minv->Write();
     h_FakeRate_tf_minv_high->Write();
     h_FakeRate_tf_qt->Write();
     h_FakeRate_tf_deltaPhi->Write();
     h_FakeRate_tf_deltaEta->Write();
     h_FakeRate_tf_deltaR->Write();
     h_FakeRate_tf_cosThetaStar->Write();
     h_FakeRate_tf_pt1_noweight->Write();
     h_FakeRate_tf_pt2_noweight->Write();
     h_FakeRate_ft_pt1_noweight->Write();
     h_FakeRate_ft_pt2_noweight->Write();
     h_FakeRate_ff_pt1->Write();
     h_FakeRate_ff_pt1_zoom->Write();
     h_FakeRate_ff_eta1->Write();
     h_FakeRate_ff_phi1->Write();
     h_FakeRate_ff_pt2->Write();
     h_FakeRate_ff_pt2_zoom->Write();
     h_FakeRate_ff_eta2->Write();
     h_FakeRate_ff_phi2->Write();
     h_FakeRate_ff_minv->Write();
     h_FakeRate_ff_minv_high->Write();
     h_FakeRate_ff_qt->Write();
     h_FakeRate_ff_deltaPhi->Write();
     h_FakeRate_ff_deltaEta->Write();
     h_FakeRate_ff_deltaR->Write();
     h_FakeRate_ff_cosThetaStar->Write();
     h_FakeRate_ff_pt1_noweight->Write();
     h_FakeRate_ff_pt2_noweight->Write();
     h_trkIsoVRho1->Write();
     h_ecalIsoVRho1->Write();
     h_hcalIsoVRho1->Write();
     h_trkIsoVRho2->Write();
     h_ecalIsoVRho2->Write();
     h_hcalIsoVRho2->Write();
     h_FakeRate_tt_minv_120to200->Write();
     h_FakeRate_tt_minv_200to500->Write();
     h_FakeRate_tt_minv_500to800->Write();
     h_FakeRate_tt_minv_800toInf->Write();
     h_FakeRate_tf_minv_120to200->Write();
     h_FakeRate_tf_minv_200to500->Write();
     h_FakeRate_tf_minv_500to800->Write();
     h_FakeRate_tf_minv_800toInf->Write();   
     h_FakeRate_ff_minv_120to200->Write();
     h_FakeRate_ff_minv_200to500->Write();
     h_FakeRate_ff_minv_500to800->Write();
     h_FakeRate_ff_minv_800toInf->Write();

     _outT0->Write();  
     _outF->Write();
     _outF->Close();
   }
   std::cout << _fakeStatus << " PASSING " << nPass << std::endl;
   return;
}
