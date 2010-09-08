//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug 25 13:11:30 2010 by ROOT version 5.26/00
// from TTree fTreeData/PhotonTree
// found on file: diphotonTree_PhotonJet_Pt30to50.root
//////////////////////////////////////////////////////////

#ifndef fTreeData_h
#define fTreeData_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TH1.h"
#include "TH2.h"

#include <iomanip>
#include <iostream>
#include <math.h>
#include <string.h>


class fTreeData {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Event_run;
   Int_t           Event_LS;
   Int_t           Event_evnum;
   Int_t           Vtx_Nvtx;
   Double_t        Vtx_vx;
   Double_t        Vtx_vy;
   Double_t        Vtx_vz;
   Int_t           Vtx_isFake;
   Int_t           Vtx_Ntracks;
   Double_t        Vtx_sumPtTracks;
   Double_t        Vtx_ndof;
   Double_t        Vtx_d0;
   Double_t        BeamSpot_x0;
   Double_t        BeamSpot_y0;
   Double_t        BeamSpot_z0;
   Double_t        BeamSpot_sigmaZ;
   Double_t        BeamSpot_x0error;
   Double_t        BeamSpot_y0error;
   Double_t        BeamSpot_z0error;
   Double_t        BeamSpot_sigmaZ0error;
   Int_t           TrigHLT_HLT_MinBiasBSC;
   Int_t           TrigHLT_HLT_MinBiasBSC_NoBPTX;
   Int_t           TrigHLT_HLT_MinBiasBSC_OR;
   Int_t           TrigHLT_HLT_L1_BscMinBiasOR_BptxPlusORMinus;
   Int_t           TrigHLT_HLT_L1SingleEG2;
   Int_t           TrigHLT_HLT_L1SingleEG5;
   Int_t           TrigHLT_HLT_L1SingleEG8;
   Int_t           TrigHLT_HLT_L1DoubleEG5;
   Int_t           TrigHLT_HLT_Photon10_L1R;
   Int_t           TrigHLT_HLT_Photon10_Cleaned_L1R;
   Int_t           TrigHLT_HLT_Photon15_L1R;
   Int_t           TrigHLT_HLT_Photon15_Cleaned_L1R;
   Int_t           TrigHLT_HLT_Photon15_LooseEcalIso_L1R;
   Int_t           TrigHLT_HLT_Photon15_LooseEcalIso_Cleaned_L1R;
   Int_t           TrigHLT_HLT_Photon15_TrackIso_L1R;
   Int_t           TrigHLT_HLT_Photon15_TrackIso_Cleaned_L1R;
   Int_t           TrigHLT_HLT_Photon20_L1R;
   Int_t           TrigHLT_HLT_Photon20_Cleaned_L1R;
   Int_t           TrigHLT_HLT_Photon30_L1R;
   Int_t           TrigHLT_HLT_Photon30_Cleaned_L1R;
   Int_t           TrigHLT_HLT_Photon30_L1R_8E29;
   Int_t           TrigHLT_HLT_Photon50_L1R;
   Int_t           TrigHLT_HLT_Photon50_Cleaned_L1R;
   Int_t           TrigHLT_HLT_DoublePhoton5_L1R;
   Int_t           TrigHLT_HLT_DoublePhoton10_L1R;
   Int_t           TrigHLT_HLT_DoublePhoton15_L1R;
   Int_t           TrigHLT_HLT_DoublePhoton20_L1R;
   Double_t        Photon1_pt;
   Double_t        Photon1_eta;
   Double_t        Photon1_phi;
   Double_t        Photon1_detEta;
   Double_t        Photon1_detPhi;
   Int_t           Photon1_detId;
   Int_t           Photon1_iEtaY;
   Int_t           Photon1_iPhiX;
   Double_t        Photon1_vx;
   Double_t        Photon1_vy;
   Double_t        Photon1_vz;
   Double_t        Photon1_r9;
   Double_t        Photon1_sigmaIetaIeta;
   Double_t        Photon1_sigmaEtaEta;
   Double_t        Photon1_maxEnergyXtal;
   Double_t        Photon1_e1x5;
   Double_t        Photon1_e2x5;
   Double_t        Photon1_e3x3;
   Double_t        Photon1_e5x5;
   Double_t        Photon1_r1x5;
   Double_t        Photon1_r2x5;
   Double_t        Photon1_swisscross;
   Double_t        Photon1_eMax;
   Double_t        Photon1_eLeft;
   Double_t        Photon1_eRight;
   Double_t        Photon1_eTop;
   Double_t        Photon1_eBottom;
   Double_t        Photon1_eSecond;
   Int_t           Photon1_severityLevel;
   Int_t           Photon1_recHitFlag;
   Double_t        Photon1_maxRecHitTime;
   Double_t        Photon1_hadOverEm;
   Double_t        Photon1_hadDepth1OverEm;
   Double_t        Photon1_hadDepth2OverEm;
   Float_t         Photon1_hcalIso04;
   Float_t         Photon1_hcalIso03;
   Float_t         Photon1_ecalIso04;
   Float_t         Photon1_ecalIso03;
   Float_t         Photon1_trkIsoSumPtHollow04;
   Float_t         Photon1_trkIsoSumPtSolid04;
   Int_t           Photon1_trkIsoNtrksHollow04;
   Int_t           Photon1_trkIsoNtrksSolid04;
   Float_t         Photon1_trkIsoSumPtHollow03;
   Float_t         Photon1_trkIsoSumPtSolid03;
   Int_t           Photon1_trkIsoNtrksHollow03;
   Int_t           Photon1_trkIsoNtrksSolid03;
   Float_t         Photon1_esRatio;
   Double_t        Photon1_scRawEnergy;
   Double_t        Photon1_scPreshowerEnergy;
   Double_t        Photon1_scPhiWidth;
   Double_t        Photon1_scEtaWidth;
   Int_t           Photon1_scNumBasicClusters;
   Bool_t          Photon1_isEB;
   Bool_t          Photon1_isEE;
   Bool_t          Photon1_isEBEtaGap;
   Bool_t          Photon1_isEBPhiGap;
   Bool_t          Photon1_isEERingGap;
   Bool_t          Photon1_isEEDeeGap;
   Bool_t          Photon1_isEBEEGap;
   Bool_t          Photon1_hasPixelSeed;
   Double_t        Photon2_pt;
   Double_t        Photon2_eta;
   Double_t        Photon2_phi;
   Double_t        Photon2_detEta;
   Double_t        Photon2_detPhi;
   Int_t           Photon2_detId;
   Int_t           Photon2_iEtaY;
   Int_t           Photon2_iPhiX;
   Double_t        Photon2_vx;
   Double_t        Photon2_vy;
   Double_t        Photon2_vz;
   Double_t        Photon2_r9;
   Double_t        Photon2_sigmaIetaIeta;
   Double_t        Photon2_sigmaEtaEta;
   Double_t        Photon2_maxEnergyXtal;
   Double_t        Photon2_e1x5;
   Double_t        Photon2_e2x5;
   Double_t        Photon2_e3x3;
   Double_t        Photon2_e5x5;
   Double_t        Photon2_r1x5;
   Double_t        Photon2_r2x5;
   Double_t        Photon2_swisscross;
   Double_t        Photon2_eMax;
   Double_t        Photon2_eLeft;
   Double_t        Photon2_eRight;
   Double_t        Photon2_eTop;
   Double_t        Photon2_eBottom;
   Double_t        Photon2_eSecond;
   Int_t           Photon2_severityLevel;
   Int_t           Photon2_recHitFlag;
   Double_t        Photon2_maxRecHitTime;
   Double_t        Photon2_hadOverEm;
   Double_t        Photon2_hadDepth1OverEm;
   Double_t        Photon2_hadDepth2OverEm;
   Float_t         Photon2_hcalIso04;
   Float_t         Photon2_hcalIso03;
   Float_t         Photon2_ecalIso04;
   Float_t         Photon2_ecalIso03;
   Float_t         Photon2_trkIsoSumPtHollow04;
   Float_t         Photon2_trkIsoSumPtSolid04;
   Int_t           Photon2_trkIsoNtrksHollow04;
   Int_t           Photon2_trkIsoNtrksSolid04;
   Float_t         Photon2_trkIsoSumPtHollow03;
   Float_t         Photon2_trkIsoSumPtSolid03;
   Int_t           Photon2_trkIsoNtrksHollow03;
   Int_t           Photon2_trkIsoNtrksSolid03;
   Float_t         Photon2_esRatio;
   Double_t        Photon2_scRawEnergy;
   Double_t        Photon2_scPreshowerEnergy;
   Double_t        Photon2_scPhiWidth;
   Double_t        Photon2_scEtaWidth;
   Int_t           Photon2_scNumBasicClusters;
   Bool_t          Photon2_isEB;
   Bool_t          Photon2_isEE;
   Bool_t          Photon2_isEBEtaGap;
   Bool_t          Photon2_isEBPhiGap;
   Bool_t          Photon2_isEERingGap;
   Bool_t          Photon2_isEEDeeGap;
   Bool_t          Photon2_isEBEEGap;
   Bool_t          Photon2_hasPixelSeed;
//   Int_t           MCMatchPhoton1_Status3_status;
//   Int_t           MCMatchPhoton1_Status3_PdgId;
//   Int_t           MCMatchPhoton1_Status3_MotherPdgId;
//   Int_t           MCMatchPhoton1_Status3_GrandmotherPdgId;
//   Double_t        MCMatchPhoton1_Status3_pt;
//   Double_t        MCMatchPhoton1_Status3_eta;
//   Double_t        MCMatchPhoton1_Status3_phi;
//   Int_t           MCMatchPhoton2_Status3_status;
//   Int_t           MCMatchPhoton2_Status3_PdgId;
//   Int_t           MCMatchPhoton2_Status3_MotherPdgId;
//   Int_t           MCMatchPhoton2_Status3_GrandmotherPdgId;
//   Double_t        MCMatchPhoton2_Status3_pt;
//   Double_t        MCMatchPhoton2_Status3_eta;
//   Double_t        MCMatchPhoton2_Status3_phi;
//   Int_t           MCMatchPhoton1_Status1_status;
//   Int_t           MCMatchPhoton1_Status1_PdgId;
//   Int_t           MCMatchPhoton1_Status1_MotherPdgId;
//   Int_t           MCMatchPhoton1_Status1_GrandmotherPdgId;
//   Double_t        MCMatchPhoton1_Status1_pt;
//   Double_t        MCMatchPhoton1_Status1_eta;
//   Double_t        MCMatchPhoton1_Status1_phi;
//   Int_t           MCMatchPhoton2_Status1_status;
//   Int_t           MCMatchPhoton2_Status1_PdgId;
//   Int_t           MCMatchPhoton2_Status1_MotherPdgId;
//   Int_t           MCMatchPhoton2_Status1_GrandmotherPdgId;
//   Double_t        MCMatchPhoton2_Status1_pt;
//   Double_t        MCMatchPhoton2_Status1_eta;
//   Double_t        MCMatchPhoton2_Status1_phi;
   Double_t        Diphoton_Minv;
   Double_t        Diphoton_qt;
   Double_t        Diphoton_deltaPhi;
   Double_t        Diphoton_deltaEta;
   Double_t        Diphoton_deltaR;
   Double_t        Diphoton_cosThetaStar;

   // List of branches
   TBranch        *b_Event;   //!
   TBranch        *b_Vtx;   //!
   TBranch        *b_BeamSpot;   //!
   TBranch        *b_TrigHLT;   //!
   TBranch        *b_Photon1;   //!
   TBranch        *b_Photon2;   //!
//   TBranch        *b_MCMatchPhoton1_Status3;   //!
//   TBranch        *b_MCMatchPhoton2_Status3;   //!
//   TBranch        *b_MCMatchPhoton1_Status1;   //!
//   TBranch        *b_MCMatchPhoton2_Status1;   //!
   TBranch        *b_Diphoton;   //!

   // diphoton analysis
   Float_t _cutPhoton1Pt;
   Float_t _cutPhoton2Pt;

   TFile*  _outF;
   
   TH1* h_TrigHLT;

   TH1* h_Diphoton_Minv;
   TH1* h_Diphoton_qt;
   TH1* h_Diphoton_deltaPhi;
   TH1* h_Diphoton_deltaEta;
   TH1* h_Diphoton_deltaR;

   TH1* h_Photon1_pt; 
   TH1* h_Photon1_eta; 
   TH1* h_Photon1_phi; 
   
   TH1* h_Photon1_r9; 
   TH1* h_Photon1_sigmaIetaIeta; 
   TH1* h_Photon1_sigmaEtaEta; 
   
   TH1* h_Photon1_swisscross; 
   TH1* h_Photon1_severityLevel; 
   TH1* h_Photon1_recHitFlag; 
   TH1* h_Photon1_maxRecHitTime; 
   
   TH1* h_Photon1_hadOverEm; 
   TH1* h_Photon1_hadDepth1OverEm; 
   TH1* h_Photon1_hadDepth2OverEm; 
   
   TH1* h_Photon1_hcalIso04; 
   TH1* h_Photon1_hcalIso03; 
   
   TH1* h_Photon1_ecalIso04; 
   TH1* h_Photon1_ecalIso03; 
   
   TH1* h_Photon1_trkIsoSumPtHollow04; 
   TH1* h_Photon1_trkIsoSumPtSolid04; 
   TH1* h_Photon1_trkIsoNtrksHollow04; 
   TH1* h_Photon1_trkIsoNtrksSolid04; 
   
   TH1* h_Photon1_trkIsoSumPtHollow03; 
   TH1* h_Photon1_trkIsoSumPtSolid03; 
   TH1* h_Photon1_trkIsoNtrksHollow03; 
   TH1* h_Photon1_trkIsoNtrksSolid03; 
   
   TH1* h_Photon1_esRatio; 


   TH1* h_Photon2_pt; 
   TH1* h_Photon2_eta; 
   TH1* h_Photon2_phi; 
   
   TH1* h_Photon2_r9; 
   TH1* h_Photon2_sigmaIetaIeta; 
   TH1* h_Photon2_sigmaEtaEta; 
   
   TH1* h_Photon2_swisscross; 
   TH1* h_Photon2_severityLevel; 
   TH1* h_Photon2_recHitFlag; 
   TH1* h_Photon2_maxRecHitTime; 
   
   TH1* h_Photon2_hadOverEm; 
   TH1* h_Photon2_hadDepth1OverEm; 
   TH1* h_Photon2_hadDepth2OverEm; 
   
   TH1* h_Photon2_hcalIso04; 
   TH1* h_Photon2_hcalIso03; 
   
   TH1* h_Photon2_ecalIso04; 
   TH1* h_Photon2_ecalIso03; 
   
   TH1* h_Photon2_trkIsoSumPtHollow04; 
   TH1* h_Photon2_trkIsoSumPtSolid04; 
   TH1* h_Photon2_trkIsoNtrksHollow04; 
   TH1* h_Photon2_trkIsoNtrksSolid04; 
   
   TH1* h_Photon2_trkIsoSumPtHollow03; 
   TH1* h_Photon2_trkIsoSumPtSolid03; 
   TH1* h_Photon2_trkIsoNtrksHollow03; 
   TH1* h_Photon2_trkIsoNtrksSolid03; 
   
   TH1* h_Photon2_esRatio; 

   fTreeData(TTree *tree=0);
   virtual ~fTreeData();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef fTreeData_cxx

fTreeData::fTreeData(TTree *tree)
  :_cutPhoton1Pt(10)
  , _cutPhoton2Pt(10)
  , _outF(0)
	    //  , h_TrigHLT(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  if (tree == 0) {
    std::cerr << "Give me a tree ! " << std::endl;
    return;
  }
  Init(tree);

}

fTreeData::~fTreeData()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t fTreeData::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t fTreeData::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void fTreeData::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event", &Event_run, &b_Event);
   fChain->SetBranchAddress("Vtx", &Vtx_Nvtx, &b_Vtx);
   fChain->SetBranchAddress("BeamSpot", &BeamSpot_x0, &b_BeamSpot);
   fChain->SetBranchAddress("TrigHLT", &TrigHLT_HLT_MinBiasBSC, &b_TrigHLT);
   fChain->SetBranchAddress("Photon1", &Photon1_pt, &b_Photon1);
   fChain->SetBranchAddress("Photon2", &Photon2_pt, &b_Photon2);
//   fChain->SetBranchAddress("MCMatchPhoton1_Status3", &MCMatchPhoton1_Status3_status, &b_MCMatchPhoton1_Status3);
//   fChain->SetBranchAddress("MCMatchPhoton2_Status3", &MCMatchPhoton2_Status3_status, &b_MCMatchPhoton2_Status3);
//   fChain->SetBranchAddress("MCMatchPhoton1_Status1", &MCMatchPhoton1_Status1_status, &b_MCMatchPhoton1_Status1);
//   fChain->SetBranchAddress("MCMatchPhoton2_Status1", &MCMatchPhoton2_Status1_status, &b_MCMatchPhoton2_Status1);
   fChain->SetBranchAddress("Diphoton", &Diphoton_Minv, &b_Diphoton);
   Notify();

   h_TrigHLT = new TH1F("h_TrigHLT","HLT Results" , 23, 0., 23);
   //   h_TrigHLT->GetXaxis()->SetBinLabel(1,"HLT_MinBiasBSC");
   //   h_TrigHLT->GetXaxis()->SetBinLabel(2,"HLT_MinBiasBSC_NoBPTX");
   //   h_TrigHLT->GetXaxis()->SetBinLabel(3,"HLT_MinBiasBSC_OR");
   //   h_TrigHLT->GetXaxis()->SetBinLabel(4,"HLT_L1_BscMinBiasOR_BptxPlusORMinus");
   h_TrigHLT->GetXaxis()->SetBinLabel(1,"HLT_L1SingleEG2");
   h_TrigHLT->GetXaxis()->SetBinLabel(2,"HLT_L1SingleEG5");
   h_TrigHLT->GetXaxis()->SetBinLabel(3,"HLT_L1SingleEG8");
   h_TrigHLT->GetXaxis()->SetBinLabel(4,"HLT_L1DoubleEG5");
   h_TrigHLT->GetXaxis()->SetBinLabel(5,"HLT_Photon10_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(6,"HLT_Photon10_Cleaned_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(7,"HLT_Photon15_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(8,"HLT_Photon15_Cleaned_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(9,"HLT_Photon15_LooseEcalIso_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(10,"HLT_Photon15_LooseEcalIso_Cleaned_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(11,"HLT_Photon15_TrackIso_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(12,"HLT_Photon15_TrackIso_Cleaned_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(13,"HLT_Photon20_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(14,"HLT_Photon20_Cleaned_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(15,"HLT_Photon30_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(16,"HLT_Photon30_Cleaned_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(17,"HLT_Photon30_L1R_8E29");
   h_TrigHLT->GetXaxis()->SetBinLabel(18,"HLT_Photon50_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(19,"HLT_Photon50_Cleaned_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(20,"HLT_DoublePhoton5_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(21,"HLT_DoublePhoton10_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(22,"HLT_DoublePhoton15_L1R");
   h_TrigHLT->GetXaxis()->SetBinLabel(23,"HLT_DoublePhoton20_L1R");   

   h_Diphoton_Minv = new TH1F("h_Diphoton_Minv","Diphoton Invariant Mass;m(#gamma#gamma) [GeV]",100,0,500);
   h_Diphoton_qt = new TH1F("h_Diphoton_qt","Diphoton qt;#gamma#gamma qt [GeV]",30,0,30);
   h_Diphoton_deltaPhi = new TH1F("h_Diphoton_deltaPhi","Diphoton #Delta#phi;#gamma#gamma #Delta#phi",180,-3.14159,3.14159); 
   h_Diphoton_deltaEta = new TH1F("h_Diphoton_deltaEta","Diphoton #Delta#eta;#gamma#gamma #Delta#eta",86,-6.,6.); 
   h_Diphoton_deltaR = new TH1F("h_Diphoton_deltaR","Diphoton #DeltaR; #gamma#gamma #DeltaR",70,0,7.); 

   h_Photon1_pt = new TH1F("h_Photon1_pt","Photon1 pt;#gamma_{1} p_{T} [GeV]",100,0,1000); 
   h_Photon1_eta = new TH1F("h_Photon1_eta","Photon1 #eta;#gamma_{1} #eta",86,-1.5,1.5); 
   h_Photon1_phi = new TH1F ("h_Photon1_phi","Photon1 #phi;#gamma_{1} #phi",180,-3.14159,3.14159); 

   h_Photon1_r9 = new TH1F("h_Photon1_r9","Photon1 R9;#gamma_{1} R9", 50, 0.1, 1.5);
   h_Photon1_sigmaIetaIeta = new TH1F("h_Photon1_sigmaIetaIeta","Photon1 #sigma_{i#etai#eta};#gamma_{1} #sigma_{i#etai#eta}", 50, 0, 0.10);
   h_Photon1_sigmaEtaEta = new TH1F("h_Photon1_sigmaEtaEta","Photon1 #sigma_{#eta#eta};#gamma_{1} #sigma_{#eta#eta}", 50, 0, 0.10);

   h_Photon1_swisscross = new TH1F("h_Photon1_swisscross","Photon1 swisscross;#gamma_{1} 1-E_{4}/E_{1}", 60, 0, 1.2);
   h_Photon1_severityLevel = new TH1F("h_Photon1_severityLevel","Photon1 severityLevel", 5, 0, 5);
   h_Photon1_severityLevel->GetXaxis()->SetBinLabel(1,"kGood");
   h_Photon1_severityLevel->GetXaxis()->SetBinLabel(2,"kProblematic");
   h_Photon1_severityLevel->GetXaxis()->SetBinLabel(3,"kRecovered");
   //   h_Photon1_severityLevel->GetXaxis()->SetBinLabel(4,"kTime");
   //   h_Photon1_severityLevel->GetXaxis()->SetBinLabel(5,"kWeird");
   //   h_Photon1_severityLevel->GetXaxis()->SetBinLabel(6,"kBad");
   h_Photon1_severityLevel->GetXaxis()->SetBinLabel(4,"kWeird");
   h_Photon1_severityLevel->GetXaxis()->SetBinLabel(5,"kBad");

   h_Photon1_recHitFlag = new TH1F("h_Photon1_recHitFlag","Photon1 recHitFlag", 15, 0, 15);

   // watch out, as this ordering changes with release
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(1,"kGood");                  // channel ok, the energy and time measurement are reliable
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(2,"kPoorReco");              // the energy is available from the UncalibRecHit, but approximate (bad shape, large chi2)
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(3,"kOutOfTime");             // the energy is available from the UncalibRecHit (sync reco), but the event is out of time
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(4,"kFaultyHardware");        // The energy is available from the UncalibRecHit, channel is faulty at some hardware level (e.g. noisy)
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(5,"kNoisy");                 // the channel is very noisy
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(6,"kPoorCalib");             // the energy is available from the UncalibRecHit, but the calibration of the channel is poor
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(7,"kSaturated");             // saturated channel (recovery not tried)
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(8,"kLeadingEdgeRecovered");  // saturated channel: energy estimated from the leading edge before saturation
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(9,"kNeighboursRecovered");   // saturated/isolated dead: energy estimated from neighbours
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(10,"kTowerRecovered");        // channel in TT with no data link, info retrieved from Trigger Primitive
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(11,"kFake");                  // the signal in the channel is a fake (e.g. a so-called spike)
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(12,"kFakeNeighbours");        // the signal in the channel is a fake and it is detected by looking at the neighbours
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(13,"kDead");                  // channel is dead and any recovery fails
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(14,"kKilled");                // MC only flag: the channel is killed in the real detector
   h_Photon1_recHitFlag->GetXaxis()->SetBinLabel(15,"kUnknown");  

   h_Photon1_maxRecHitTime = new TH1F("h_Photon1_maxRecHitTime","Photon1 maxRecHitTime;#gamma_{1} time [ns]", 100, -20, 20); 
   
   h_Photon1_hadOverEm = new TH1F("h_Photon1_hadOverEm","Photon1 H/E;#gamma_{1} H/E", 50, 0, 0.06); 
   h_Photon1_hcalIso04 = new TH1F("h_Photon1_hcalIso04","Photon1 HCAL Iso #DeltaR 0.4;#gamma_{1} HCAL Iso #DeltaR 0.4", 100, 0, 4); 
   h_Photon1_hcalIso03 = new TH1F("h_Photon1_hcalIso03","Photon1 HCAL Iso #DeltaR 0.3;#gamma_{1} HCAL Iso #DeltaR 0.3", 100, 0, 4); 
   h_Photon1_ecalIso04 = new TH1F("h_Photon1_ecalIso04","Photon1 ECAL Iso #DeltaR 0.4;#gamma_{1} ECAL Iso #DeltaR 0.4", 50, -1, 4);
   h_Photon1_ecalIso03 = new TH1F("h_Photon1_ecalIso03","Photon1 ECAL Iso #DeltaR 0.3;#gamma_{1} ECAL Iso #DeltaR 0.3", 50, -1, 4); 
   
   h_Photon1_trkIsoSumPtHollow04 = new TH1F("h_Photon1_trkIsoSumPtHollow04","Photon1 Track Iso #Sigma p_{T} Hollow #DeltaR 0.4;#gamma_{1} Track Iso #Sigma p_{T} Hollow #DeltaR 0.4", 30, 0, 3); 
   h_Photon1_trkIsoSumPtSolid04 = new TH1F("h_Photon1_trkIsoSumPtSolid04","Photon1 Track Iso #Sigma p_{T} Solid #DeltaR 0.4;#gamma_{1} Track Iso #Sigma p_{T} Solid #DeltaR 0.4", 30, 0, 3); 
   h_Photon1_trkIsoSumPtHollow03 = new TH1F("h_Photon1_trkIsoSumPtHollow03","Photon1 Track Iso #Sigma p_{T} Hollow #DeltaR 0.3;#gamma_{1} Track Iso #Sigma p_{T} Hollow #DeltaR 0.3", 30, 0, 3); 
   h_Photon1_trkIsoSumPtSolid03 = new TH1F("h_Photon1_trkIsoSumPtSolid03","Photon1 Track Iso #Sigma p_{T} Solid #DeltaR 0.3;#gamma_{1} Track Iso #Sigma p_{T} Solid #DeltaR 0.3", 30, 0, 3); 

   h_Photon1_trkIsoNtrksHollow04 = new TH1F("h_Photon1_trkIsoNtrksHollow04","Photon1 Track Iso N_{tracks} Hollow #DeltaR 0.4;#gamma_{1} Track Iso N_{tracks} Hollow #DeltaR 0.4", 6, 0, 6); 
   h_Photon1_trkIsoNtrksSolid04 = new TH1F("h_Photon1_trkIsoNtrksSolid04","Photon1 Track Iso N_{tracks} Solid #DeltaR 0.4;#gamma_{1} Track Iso N_{tracks} Solid #DeltaR 0.4", 6, 0, 6); 
   h_Photon1_trkIsoNtrksHollow03 = new TH1F("h_Photon1_trkIsoNtrksHollow03","Photon1 Track Iso N_{tracks} Hollow #DeltaR 0.3;#gamma_{1} Track Iso N_{tracks} Hollow #DeltaR 0.3", 6, 0, 6); 
   h_Photon1_trkIsoNtrksSolid03 = new TH1F("h_Photon1_trkIsoNtrksSolid03","Photon1 Track Iso N_{tracks} Solid #DeltaR 0.3;#gamma_{1} Track Iso N_{tracks} Solid #DeltaR 0.3", 6, 0, 6); 

/*    h_Photon1_esRatio = new TH1F("h_Photon1_esRatio","Photon1_esRatio", , , ); */


   h_Photon2_pt = new TH1F("h_Photon2_pt","Photon2 pt;#gamma_{2} p_{T} [GeV]",100,0,1000); 
   h_Photon2_eta = new TH1F("h_Photon2_eta","Photon2 #eta;#gamma_{2} #eta",86,-1.5,1.5); 
   h_Photon2_phi = new TH1F ("h_Photon2_phi","Photon2 #phi;#gamma_{2} #phi",180,-3.14159,3.14159); 

   h_Photon2_r9 = new TH1F("h_Photon2_r9","Photon2 R9;#gamma_{2} R9", 50, 0.1, 1.5);
   h_Photon2_sigmaIetaIeta = new TH1F("h_Photon2_sigmaIetaIeta","Photon2 #sigma_{i#etai#eta};#gamma_{2} #sigma_{i#etai#eta}", 50, 0, 0.10);
   h_Photon2_sigmaEtaEta = new TH1F("h_Photon2_sigmaEtaEta","Photon2 #sigma_{#eta#eta};#gamma_{2} #sigma_{#eta#eta}", 50, 0, 0.10);

   h_Photon2_swisscross = new TH1F("h_Photon2_swisscross","Photon2 swisscross;#gamma_{2} 1-E_{4}/E_{1}", 60, 0, 1.2);
   h_Photon2_severityLevel = new TH1F("h_Photon2_severityLevel","Photon2 severityLevel", 5, 0, 5);
   h_Photon2_severityLevel->GetXaxis()->SetBinLabel(1,"kGood");
   h_Photon2_severityLevel->GetXaxis()->SetBinLabel(2,"kProblematic");
   h_Photon2_severityLevel->GetXaxis()->SetBinLabel(3,"kRecovered");
   //   h_Photon2_severityLevel->GetXaxis()->SetBinLabel(4,"kTime");
   //   h_Photon2_severityLevel->GetXaxis()->SetBinLabel(5,"kWeird");
   //   h_Photon2_severityLevel->GetXaxis()->SetBinLabel(6,"kBad");
   h_Photon2_severityLevel->GetXaxis()->SetBinLabel(4,"kWeird");
   h_Photon2_severityLevel->GetXaxis()->SetBinLabel(5,"kBad");

   h_Photon2_recHitFlag = new TH1F("h_Photon2_recHitFlag","Photon2 recHitFlag", 15, 0, 15);

   // watch out, as this ordering changes with release
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(1,"kGood");                  // channel ok, the energy and time measurement are reliable
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(2,"kPoorReco");              // the energy is available from the UncalibRecHit, but approximate (bad shape, large chi2)
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(3,"kOutOfTime");             // the energy is available from the UncalibRecHit (sync reco), but the event is out of time
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(4,"kFaultyHardware");        // The energy is available from the UncalibRecHit, channel is faulty at some hardware level (e.g. noisy)
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(5,"kNoisy");                 // the channel is very noisy
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(6,"kPoorCalib");             // the energy is available from the UncalibRecHit, but the calibration of the channel is poor
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(7,"kSaturated");             // saturated channel (recovery not tried)
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(8,"kLeadingEdgeRecovered");  // saturated channel: energy estimated from the leading edge before saturation
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(9,"kNeighboursRecovered");   // saturated/isolated dead: energy estimated from neighbours
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(10,"kTowerRecovered");        // channel in TT with no data link, info retrieved from Trigger Primitive
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(11,"kFake");                  // the signal in the channel is a fake (e.g. a so-called spike)
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(12,"kFakeNeighbours");        // the signal in the channel is a fake and it is detected by looking at the neighbours
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(13,"kDead");                  // channel is dead and any recovery fails
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(14,"kKilled");                // MC only flag: the channel is killed in the real detector
   h_Photon2_recHitFlag->GetXaxis()->SetBinLabel(15,"kUnknown");  

   h_Photon2_maxRecHitTime = new TH1F("h_Photon2_maxRecHitTime","Photon2 maxRecHitTime;#gamma_{2} time [ns]", 100, -20, 20); 
   
   h_Photon2_hadOverEm = new TH1F("h_Photon2_hadOverEm","Photon2 H/E;#gamma_{2} H/E", 50, 0, 0.06); 
   h_Photon2_hcalIso04 = new TH1F("h_Photon2_hcalIso04","Photon2 HCAL Iso #DeltaR 0.4;#gamma_{2} HCAL Iso #DeltaR 0.4", 100, 0, 4); 
   h_Photon2_hcalIso03 = new TH1F("h_Photon2_hcalIso03","Photon2 HCAL Iso #DeltaR 0.3;#gamma_{2} HCAL Iso #DeltaR 0.3", 100, 0, 4); 
   h_Photon2_ecalIso04 = new TH1F("h_Photon2_ecalIso04","Photon2 ECAL Iso #DeltaR 0.4;#gamma_{2} ECAL Iso #DeltaR 0.4", 50, -1, 4);
   h_Photon2_ecalIso03 = new TH1F("h_Photon2_ecalIso03","Photon2 ECAL Iso #DeltaR 0.3;#gamma_{2} ECAL Iso #DeltaR 0.3", 50, -1, 4); 
   
   h_Photon2_trkIsoSumPtHollow04 = new TH1F("h_Photon2_trkIsoSumPtHollow04","Photon2 Track Iso #Sigma p_{T} Hollow #DeltaR 0.4;#gamma_{2} Track Iso #Sigma p_{T} Hollow #DeltaR 0.4", 30, 0, 3); 
   h_Photon2_trkIsoSumPtSolid04 = new TH1F("h_Photon2_trkIsoSumPtSolid04","Photon2 Track Iso #Sigma p_{T} Solid #DeltaR 0.4;#gamma_{2} Track Iso #Sigma p_{T} Solid #DeltaR 0.4", 30, 0, 3); 
   h_Photon2_trkIsoSumPtHollow03 = new TH1F("h_Photon2_trkIsoSumPtHollow03","Photon2 Track Iso #Sigma p_{T} Hollow #DeltaR 0.3;#gamma_{2} Track Iso #Sigma p_{T} Hollow #DeltaR 0.3", 30, 0, 3); 
   h_Photon2_trkIsoSumPtSolid03 = new TH1F("h_Photon2_trkIsoSumPtSolid03","Photon2 Track Iso #Sigma p_{T} Solid #DeltaR 0.3;#gamma_{2} Track Iso #Sigma p_{T} Solid #DeltaR 0.3", 30, 0, 3); 

   h_Photon2_trkIsoNtrksHollow04 = new TH1F("h_Photon2_trkIsoNtrksHollow04","Photon2 Track Iso N_{tracks} Hollow #DeltaR 0.4;#gamma_{2} Track Iso N_{tracks} Hollow #DeltaR 0.4", 6, 0, 6); 
   h_Photon2_trkIsoNtrksSolid04 = new TH1F("h_Photon2_trkIsoNtrksSolid04","Photon2 Track Iso N_{tracks} Solid #DeltaR 0.4;#gamma_{2} Track Iso N_{tracks} Solid #DeltaR 0.4", 6, 0, 6); 
   h_Photon2_trkIsoNtrksHollow03 = new TH1F("h_Photon2_trkIsoNtrksHollow03","Photon2 Track Iso N_{tracks} Hollow #DeltaR 0.3;#gamma_{2} Track Iso N_{tracks} Hollow #DeltaR 0.3", 6, 0, 6); 
   h_Photon2_trkIsoNtrksSolid03 = new TH1F("h_Photon2_trkIsoNtrksSolid03","Photon2 Track Iso N_{tracks} Solid #DeltaR 0.3;#gamma_{2} Track Iso N_{tracks} Solid #DeltaR 0.3", 6, 0, 6); 

/*    h_Photon2_esRatio = new TH1F("h_Photon2_esRatio","Photon2_esRatio", , , ); */


}

Bool_t fTreeData::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void fTreeData::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t fTreeData::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
  std::cout << entry << std::endl;
   return 1;
}
#endif // #ifdef fTreeData_cxx
