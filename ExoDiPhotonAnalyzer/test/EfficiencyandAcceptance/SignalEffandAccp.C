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

TString TreeName[2] = {"diphoton_tree_RSGravToGG_kMpl-001_M-750_TuneZ2star_8TeV-pythia6_Sept10thPileUp.root","diphoton_tree_RSGravToGG_kMpl-01_M-3250_TuneZ2star_8TeV-pythia6_Septh10thPileUp.root"};

void CreateAcceptanceandEffHists( )
{
  TString outputfilename[2];
  TString inputfilename[2];
  TFile*  OutFile[2];
  TFile* infile[2];
  for (Int_t i=0;i<2;i++){
    outputfilename[i] = "histograms_"+TreeName[i];
    OutFile[i] = new TFile(outputfilename[i],"RECREATE");
    inputfilename[i] = TreeName[i];
    infile[i]=TFile::Open(inputfilename[i]);
    
    TTree* Tree = (TTree*)infile[i]->Get("diphotonSignalMCAnalyzer/fTree");            
    SignalEffLoop* treereader = new SignalEffLoop(Tree);
    treereader->_outfilename = OutFile[i];
    treereader->Loop();
  }
}  
void GetEffandAccNum() 
{  
  ofstream filename("EffOutput.txt",ios::app);
  filename<<"Sample"<<setw(40)<<"Acceptance"<<setw(40)<<"Efficiency"<<setw(40)<<"AcceptancetimesEfficiency"<<endl<<endl;
  Double_t Efficiency[2];
  Double_t Acceptance[2];
  Double_t AcceptancetimesEfficiency[2];   
  TString inputfilename[2];
  TFile* inputfile[2];
  TH1F* h_Diphoton_Minv[2];
  TH1F* h_Diphoton_Minv_AccPassed[2];
  TH1F* h_Diphoton_Minv_EffandAccPassed[2];
  
  for(Int_t i=0;i<2;i++){
    inputfilename[i] = "histograms_"+TreeName[i];
    inputfile[i] = TFile::Open(inputfilename[i]);
    h_Diphoton_Minv[i]=(TH1F*)inputfile[i]->Get("h_Diphoton_Minv");
    h_Diphoton_Minv_AccPassed[i] = (TH1F*)inputfile[i]->Get("h_Diphoton_Minv_AccPassed");
    h_Diphoton_Minv_EffandAccPassed[i] = (TH1F*)inputfile[i]->Get("h_Diphoton_Minv_EffandAccPassed"); 
    
    
    Acceptance[i] = (h_Diphoton_Minv_AccPassed[i]->GetEntries())/(h_Diphoton_Minv[i]->GetEntries());
    AcceptancetimesEfficiency[i] = (h_Diphoton_Minv_EffandAccPassed[i]->GetEntries())/(h_Diphoton_Minv[i]->GetEntries());
    
    Efficiency[i] = AcceptancetimesEfficiency[i]/Acceptance[i];
    
    filename<<TreeName[i]<<setw(20)<<Acceptance[i]<<setw(20)<<Efficiency[i]<<setw(20)<<AcceptancetimesEfficiency[i]<<endl<<endl;
    
  }
  filename.close(); 
}
  
