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
#include "CalcSignalEfficiency.C"

using namespace std;

void CalcEfficiencyForFile(TString inputFile)
{
 
  //cout << "Entering CalcEfficiencies with parameters:" <<endl;
  cout<<"input file: "<<inputFile.Data()
      <<endl;

  TChain *chain = new TChain("diphotonSignalMCAnalyzer/fTree");
  chain->Add(inputFile.Data());
  //cout << "#entries = " << chain->GetEntries() <<endl;

  CalcSignalEfficiency* efficalculator = new CalcSignalEfficiency(chain);
  efficalculator->Loop();
  
}


void CalculateAllEfficiencies()
{

CalcEfficiencyForFile("/afs/cern.ch/work/c/charaf/public/diphoton_tree_RSGravToGG_kMpl-001_M-750_TuneZ2star_8TeV-pythia6_private_Reco_PUSept10th.root");
CalcEfficiencyForFile("/afs/cern.ch/work/c/charaf/public/diphoton_tree_RSGravToGG_kMpl-001_M-1000_TuneZ2star_8TeV-pythia6_private_Reco_PUSept10th.root");
CalcEfficiencyForFile("/afs/cern.ch/work/c/charaf/public/diphoton_tree_RSGravToGG_kMpl-001_M-1250_TuneZ2star_8TeV-pythia6_private_Reco_PUSept10th.root");
CalcEfficiencyForFile("/afs/cern.ch/work/c/charaf/public/diphoton_tree_RSGravToGG_kMpl-001_M-1500_TuneZ2star_8TeV-pythia6_private_Reco_PUSept10th.root");
CalcEfficiencyForFile("/afs/cern.ch/work/c/charaf/public/diphoton_tree_RSGravToGG_kMpl-001_M-1750_TuneZ2star_8TeV-pythia6_private_Reco_PUSept10th.root");
CalcEfficiencyForFile("/afs/cern.ch/work/c/charaf/public/diphoton_tree_RSGravToGG_kMpl-001_M-2000_TuneZ2star_8TeV-pythia6_private_Reco_PUSept10th.root");
CalcEfficiencyForFile("/afs/cern.ch/work/c/charaf/public/diphoton_tree_RSGravToGG_kMpl-001_M-3000_TuneZ2star_8TeV-pythia6_private_Reco_PUSept10th.root");

CalcEfficiencyForFile("/afs/cern.ch/work/c/charaf/public/diphoton_tree_RSGravToGG_kMpl-005_M-1750_TuneZ2star_8TeV-pythia6_private_Reco_PUSept10th.root");
CalcEfficiencyForFile("/afs/cern.ch/work/c/charaf/public/diphoton_tree_RSGravToGG_kMpl-005_M-2000_TuneZ2star_8TeV-pythia6_private_Reco_PUSept10th.root");
CalcEfficiencyForFile("/afs/cern.ch/work/c/charaf/public/diphoton_tree_RSGravToGG_kMpl-005_M-2500_TuneZ2star_8TeV-pythia6_private_Reco_PUSept10th.root");
CalcEfficiencyForFile("/afs/cern.ch/work/c/charaf/public/diphoton_tree_RSGravToGG_kMpl-005_M-2750_TuneZ2star_8TeV-pythia6_private_Reco_PUSept10th.root");
CalcEfficiencyForFile("/afs/cern.ch/work/c/charaf/public/diphoton_tree_RSGravToGG_kMpl-005_M-3000_TuneZ2star_8TeV-pythia6_private_Reco_PUSept10th.root");

CalcEfficiencyForFile("/afs/cern.ch/work/c/charaf/public/diphoton_tree_RSGravToGG_kMpl-01_M-1500_TuneZ2star_8TeV-pythia6_private_Reco_PUSept10th.root");
CalcEfficiencyForFile("/afs/cern.ch/work/c/charaf/public/diphoton_tree_RSGravToGG_kMpl-01_M-2250_TuneZ2star_8TeV-pythia6_private_Reco_PUSept10th.root");
CalcEfficiencyForFile("/afs/cern.ch/work/c/charaf/public/diphoton_tree_RSGravToGG_kMpl-01_M-2500_TuneZ2star_8TeV-pythia6_private_Reco_PUSept10th.root");
CalcEfficiencyForFile("/afs/cern.ch/work/c/charaf/public/diphoton_tree_RSGravToGG_kMpl-01_M-2750_TuneZ2star_8TeV-pythia6_private_Reco_PUSept10th.root");
CalcEfficiencyForFile("/afs/cern.ch/work/c/charaf/public/diphoton_tree_RSGravToGG_kMpl-01_M-3000_TuneZ2star_8TeV-pythia6_private_Reco_PUSept10th.root");
CalcEfficiencyForFile("/afs/cern.ch/work/c/charaf/public/diphoton_tree_RSGravToGG_kMpl-01_M-3250_TuneZ2star_8TeV-pythia6_private_Reco_PUSept10th.root");
CalcEfficiencyForFile("/afs/cern.ch/work/c/charaf/public/diphoton_tree_RSGravToGG_kMpl-01_M-3500_TuneZ2star_8TeV-pythia6_private_Reco_PUSept10th.root");


}
