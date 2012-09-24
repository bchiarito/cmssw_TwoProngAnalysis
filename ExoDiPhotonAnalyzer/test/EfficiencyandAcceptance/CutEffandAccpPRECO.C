#define  SignalEffandAcc
#include "CutEfficiencyLoop.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "AcceptanceCuts.C"
#include "ReconEfficiency.C"
#include  <iostream>
#include "CutEfficiencyLoop.C"
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
    outputfilename[i] = "histograms_Cuts_"+TreeName[i]+".root";
    OutFile[i] = new TFile(outputfilename[i],"RECREATE");
    inputfilename[i] = "~/diphoton/DiPhotonTrees/SignalPointsSept10th/diphoton_tree_"+TreeName[i]+"_TuneZ2star_8TeV-pythia6_private_Reco_PUSept10th.root";
    infile[i]=TFile::Open(inputfilename[i]);    
    TTree* Tree = (TTree*)infile[i]->Get("diphotonSignalMCAnalyzer/fTree");            
    CutEfficiencyLoop* treereader = new CutEfficiencyLoop(Tree);
    treereader->_outfilename = OutFile[i];
    treereader->Loop();
  }
}  

void GetEffandAccNum() 
{  
  ofstream filename("EffCutsPrivateOutputTESTTS.txt",ios::app);
    filename<<"AcceptanceEff"<<setw(15)<<"HoverEff"<<setw(15)<<"EcalIsoEff"<<setw(15)<<"HcalIsoEff"<<setw(15)<<"TrackIsoEff"<<setw(15)<<"SigmaIetaIetaEff"<<setw(15)<<"NoPixelSeed"<<endl<<endl;
  TString inputfilename[18];
  TFile* inputfile[18];
  TH1F*  h_Diphoton_Minv[18];
  TH1F*  h_Diphoton_Minv_Acceptance[18];
  TH1F*  h_Diphoton_Minv_HoverE_After[18];
  TH1F*  h_Diphoton_Minv_TrackIso_After[18];
  TH1F*  h_Diphoton_Minv_HcalIso_After[18];
  TH1F*  h_Diphoton_NoPixelSeed_After[18];
  TH1F*  h_Diphoton_Minv_EcalIso_After[18];
  TH1F*  h_Diphoton_Minv_SigmaIetaIeta_After[18];
  TH1F*  h_Diphoton_Minv_HoverE_Before[18];
  TH1F*  h_Diphoton_Minv_TrackIso_Before[18];
  TH1F*  h_Diphoton_Minv_HcalIso_Before[18];
  TH1F*  h_Diphoton_NoPixelSeed_Before[18];
  TH1F*  h_Diphoton_Minv_EcalIso_Before[18];
  TH1F*  h_Diphoton_Minv_SigmaIetaIeta_Before[18];

  
  Double_t AcceptanceEff[18];
  Double_t TrackIsoEff[18];
  Double_t HoverEff[18];
  Double_t HcalIsoEff[18];
  Double_t EcalIsoEff[18];
  Double_t SigmaIetaIetaEff[18];
  Double_t NoPixelSeedEff[18];


  for(Int_t i=0;i<18;i++){
  
    cout<<"In Loop"<<endl;
    inputfilename[i] = "histograms_Cuts_"+TreeName[i]+".root";
    inputfile[i] = TFile::Open(inputfilename[i]);

      
    h_Diphoton_Minv[i]=(TH1F*)inputfile[i]->Get("h_Diphoton_Minv");
    h_Diphoton_Minv_Acceptance[i]=(TH1F*)inputfile[i]->Get("h_Diphoton_Minv_Acceptance");
    h_Diphoton_Minv_HoverE_After[i]=(TH1F*)inputfile[i]->Get("h_Diphoton_Minv_HoverE_After");  
    h_Diphoton_Minv_TrackIso_After[i]=(TH1F*)inputfile[i]->Get("h_Diphoton_Minv_TrackIso_After");
    h_Diphoton_Minv_HcalIso_After[i]=(TH1F*)inputfile[i]->Get("h_Diphoton_Minv_HcalIso_After");
    h_Diphoton_NoPixelSeed_After[i]=(TH1F*)inputfile[i]->Get("h_Diphoton_NoPixelSeed_After");
    h_Diphoton_Minv_EcalIso_After[i]=(TH1F*)inputfile[i]->Get("h_Diphoton_Minv_EcalIso_After");
    h_Diphoton_Minv_SigmaIetaIeta_After[i]=(TH1F*)inputfile[i]->Get("h_Diphoton_Minv_SigmaIetaIeta_After");
    h_Diphoton_Minv_HoverE_Before[i]=(TH1F*)inputfile[i]->Get("h_Diphoton_Minv_HoverE_Before");
    h_Diphoton_Minv_TrackIso_Before[i]=(TH1F*)inputfile[i]->Get("h_Diphoton_Minv_TrackIso_Before");
    h_Diphoton_Minv_HcalIso_Before[i]=(TH1F*)inputfile[i]->Get("h_Diphoton_Minv_HcalIso_Before");
    cout<<"HERE"<<endl;
    h_Diphoton_NoPixelSeed_Before[i]=(TH1F*)inputfile[i]->Get("h_Diphoton_NoPixelSeed_Before");
    h_Diphoton_Minv_EcalIso_Before[i]=(TH1F*)inputfile[i]->Get("h_Diphoton_Minv_EcalIso_Before");
    h_Diphoton_Minv_SigmaIetaIeta_Before[i]=(TH1F*)inputfile[i]->Get("h_Diphoton_Minv_SigmaIetaIeta_Before");

    AcceptanceEff[i]=((h_Diphoton_Minv_Acceptance[i]->GetEntries())/(h_Diphoton_Minv[i]->GetEntries()));
    TrackIsoEff[i]=((h_Diphoton_Minv_TrackIso_After[i]->GetEntries())/(h_Diphoton_Minv_TrackIso_Before[i]->GetEntries()));
    cout<<"AfterTRCK"<<endl;
    HoverEff[i]=((h_Diphoton_Minv_HoverE_After[i]->GetEntries())/(( h_Diphoton_Minv_HoverE_Before[i]->GetEntries())));
    cout<<"AfterHOverE"<<endl;
    HcalIsoEff[i]=((h_Diphoton_Minv_HcalIso_After[i]->GetEntries())/(h_Diphoton_Minv_HcalIso_Before[i]->GetEntries()));
    cout<<"AfterHCal"<<endl;
    cout<<h_Diphoton_Minv_EcalIso_Before[i]->GetEntries()<<endl;
    EcalIsoEff[i]=((h_Diphoton_Minv_EcalIso_After[i]->GetEntries())/(h_Diphoton_Minv_EcalIso_Before[i]->GetEntries()));
    cout<<"AfterCal"<<endl;
    SigmaIetaIetaEff[i]=((h_Diphoton_Minv_SigmaIetaIeta_After[i]->GetEntries())/(h_Diphoton_Minv_SigmaIetaIeta_Before[i]->GetEntries()));
    cout<<"AfterSigmaIetaIeta"<<endl;
    cout<<h_Diphoton_NoPixelSeed_After[i]->GetEntries()<<endl;
    NoPixelSeedEff[i]=((h_Diphoton_NoPixelSeed_After[i]->GetEntries())/(h_Diphoton_NoPixelSeed_Before[i]->GetEntries()));
    cout<<"AfterNoPixelSeed"<<endl;
     filename<<AcceptanceEff[i]<<setw(15)<<HoverEff[i]<<setw(15)<<EcalIsoEff[i]<<setw(15)<<HcalIsoEff[i]<<setw(15)<<TrackIsoEff[i]<<setw(15)<<SigmaIetaIetaEff[i]<<setw(15)<<NoPixelSeedEff[i]<<endl<<endl;     
   }  
 
 filename.close(); 
}
 
