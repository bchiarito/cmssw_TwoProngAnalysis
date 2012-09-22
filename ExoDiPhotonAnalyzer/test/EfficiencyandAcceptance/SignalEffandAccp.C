#define  SignalEffandAcc
#include "SignalEffLoop.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "AcceptanceCuts.C"
#include "ReconEfficiency.C"
#include  <iostream>
#include "SignalEffLoop.C"
#include <fstream>
#include <iomanip>
#include <TMath.h>
TString TreeName[18] = {"RSGravToGG_kMpl-001_M-750","RSGravToGG_kMpl-001_M-1000","RSGravToGG_kMpl-001_M-1250","RSGravToGG_kMpl-001_M-1500","RSGravToGG_kMpl-001_M-1750","RSGravToGG_kMpl-001_M-2000","RSGravToGG_kMpl-001_M-3000","RSGravToGG_kMpl-005_M-1750","RSGravToGG_kMpl-005_M-2000","RSGravToGG_kMpl-005_M-2500","RSGravToGG_kMpl-005_M-2750","RSGravToGG_kMpl-005_M-3000","RSGravToGG_kMpl-01_M-2250","RSGravToGG_kMpl-01_M-2500","RSGravToGG_kMpl-01_M-2750","RSGravToGG_kMpl-01_M-3000","RSGravToGG_kMpl-01_M-3250","RSGravToGG_kMpl-01_M-3500"};

void CreateAcceptanceandEffHists( )
{
  TString outputfilename[18];
  TString inputfilename[18];
  TFile*  OutFile[18];
  TFile* infile[18];
  for (Int_t i=0;i<18;i++){
    outputfilename[i] = "histograms_"+TreeName[i];
    OutFile[i] = new TFile(outputfilename[i],"RECREATE");
    inputfilename[i] = "~/diphoton/DiPhotonTrees/SignalPointsSept10th/diphoton_tree_"+TreeName[i]+"_TuneZ2star_8TeV-pythia6_private_Reco_PUSept10th.root";
    infile[i]=TFile::Open(inputfilename[i]);
    
    TTree* Tree = (TTree*)infile[i]->Get("diphotonSignalMCAnalyzer/fTree");            
    SignalEffLoop* treereader = new SignalEffLoop(Tree);
    treereader->_outfilename = OutFile[i];
    treereader->Loop();
  }
}  
void GetEffandAccNum() 
{  
  ofstream filename("EffPrivateOutputRSTabletst.txt",ios::app);
   Double_t Efficiency[18];
  Double_t Acceptance[18];
  Double_t AcceptancetimesEfficiency[18];   
  TString inputfilename[18];
  TFile* inputfile[18];
  TH1F* h_Diphoton_Minv[18];
  TH1F* h_Diphoton_Minv_AccPassed[18];
  TH1F* h_Diphoton_Minv_EffandAccPassed[18];
  Double_t nbinsDiphoton_Minv[18];
  Double_t nbinsAccPassed[18];
  Double_t nbinsEffandAccPassed[18];
  Double_t HistIntMinvError[18];
  Double_t HistIntFinalMinvError[18];
  Double_t HistIntAccPassedError[18];
  Double_t HistIntFinalAccPassedError[18];
  Double_t HistIntEffandAccPassedError[18];
  Double_t HistIntFinalEffandAccPassedError[18];
  Double_t TotalError[18]; 
  Double_t AcceptanceError[18];
  Double_t AcceptancetimesEffError[18];
  Double_t EffError[18];
  
   for(Int_t i=0;i<18;i++){
        

     inputfilename[i] = "histograms_"+TreeName[i];
     inputfile[i] = TFile::Open(inputfilename[i]);
     h_Diphoton_Minv[i]=(TH1F*)inputfile[i]->Get("h_Diphoton_Minv");
     h_Diphoton_Minv_AccPassed[i] = (TH1F*)inputfile[i]->Get("h_Diphoton_Minv_AccPassed");
     h_Diphoton_Minv_EffandAccPassed[i] = (TH1F*)inputfile[i]->Get("h_Diphoton_Minv_EffandAccPassed"); 
     
     h_Diphoton_Minv[i]->Sumw2();
     h_Diphoton_Minv_AccPassed[i]->Sumw2();
     h_Diphoton_Minv_EffandAccPassed[i]->Sumw2();
     
     Acceptance[i] = (h_Diphoton_Minv_AccPassed[i]->GetEntries())/(h_Diphoton_Minv[i]->GetEntries());
     AcceptancetimesEfficiency[i] = (h_Diphoton_Minv_EffandAccPassed[i]->GetEntries())/(h_Diphoton_Minv[i]->GetEntries());
     Efficiency[i] = AcceptancetimesEfficiency[i]/Acceptance[i];
     
     cout<<"Here"<<endl; 
     
     AcceptanceError[i] = TMath::Sqrt((Acceptance[i]*(1-Acceptance[i])/(h_Diphoton_Minv[i]->GetEntries())));
     EffError[i] = TMath::Sqrt((Efficiency[i]*(1-Efficiency[i]))/(h_Diphoton_Minv[i]->GetEntries()));
     AcceptancetimesEffError[i] = TMath::Sqrt((TMath::Power(AcceptanceError[i],2)+ TMath::Power(EffError[i],2)));
     filename<<Acceptance[i]<<" "<<AcceptanceError[i]<<" "<<Efficiency[i]<<" "<<EffError[i]<<" "<<AcceptancetimesEfficiency[i]<<" "<<AcceptancetimesEffError[i]<<endl;
     
   }
   filename.close(); 
}
  
