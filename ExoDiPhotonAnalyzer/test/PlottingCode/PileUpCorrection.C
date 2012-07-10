                                   
#include <TObjArray.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TText.h>
#include "TH1.h"
#include "TH2.h"
#include <TStyle.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TF1.h>
#include <iomanip>
#include <math.h>
#include <string.h>
#include <TSystem.h>
#include "PlottingCodeLoop.C"

using namespace std;


void PileUpCorrection(TString DataPileUp="Data", TString JSON2 ="Blac",Int_t Lumi =5050)
{

  TString  MCPileUpFileName = "/afs/cern.ch/work/j/jcarson/private/DiPhotonTrees/diphoton_tree_DiPhotonBorn_Pt25to250_Summer12.root";
  TString  DataPileUpFileName = "/afs/cern.ch/work/j/jcarson/private/DiPhotonTrees/"+DataPileUp+".root";

  cout<<MCPileUpFileName.Data()<<endl;
  cout<<DataPileUpFileName.Data()<<endl;
 
  TFile* DataPileUpFile = TFile::Open(DataPileUpFileName.Data());
  TH1D* DataPileUpHisto;
 DataPileUpHisto = (TH1D*)DataPileUpFile->Get("pileup");

TChain *chain =  new TChain("diphotonAnalyzer/fTree");
 chain->Add(MCPileUpFileName.Data());
 cout<<chain<<endl;
 TH1F* MCPileUpHisto = new TH1F("MCPileUpHisto","MCPileUp",100,0,100); 
chain->Draw("pu_n >> MCPileUpHisto");
 MCPileUpHisto->Print();
Double_t normMC = MCPileUpHisto->GetEntries();
 MCPileUpHisto->Scale(1/normMC);
 Double_t normData = (DataPileUpHisto->Integral());
 cout<<normData<<endl;
 DataPileUpHisto->Scale(1/normData); 
 TH1D* PileUpCorrection = new TH1D("pileupcor","pileupcor",100,0,100);
 PileUpCorrection->Divide(DataPileUpHisto,MCPileUpHisto,1,1);
TFile* PileUpCorrectionFile = new TFile("PileUpCorrection"+JSON2+".root","RECREATE");


PileUpCorrectionFile->cd();
PileUpCorrection->Write();
 DataPileUpHisto->Write();
 MCPileUpHisto->Write();
PileUpCorrectionFile->Close();
}
