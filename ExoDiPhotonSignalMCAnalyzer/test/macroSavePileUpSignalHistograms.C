#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"

#include "TString.h"

#include <iostream>

void storePUHisto(TString inputfilepath, TString outputfilepath, TString histoname){

  TFile *inputFile = TFile::Open(inputfilepath.Data());
  //This should have always the same name
  TH1F *histo = (TH1F*)inputFile->Get("diphotonSignalMCAnalyzer/fpu_n_BeforeCuts");
  //set the name according to the name provided in cfg for crab
  histo->SetName(histoname.Data());
  histo->SetBins(80,0,80);
  //Open a new file to store this histo to
  TFile *outputFile = new TFile(outputfilepath.Data(),"recreate");
  outputFile->cd();
  histo->Write();
  outputFile->cd();
  outputFile->Close();
}
