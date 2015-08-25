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

TString XTitles[12] = {"M_{#gamma #gamma} (GeV)","q_{t} (GeV)","#Delta #phi","#Delta #eta","#Delta R","cos #theta^{*}","#gamma_{1} p_{t} (GeV)","#gamma_{1} #eta","#gamma_{1} #phi","#gamma_{2} p_{t} (GeV)","#gamma_{2} #eta","#gamma_{2} #phi"};
TString YTitles[12] = {"Events / 20 GeV","Events / 10 GeV","","","","","Entries / 20 GeV","","","Entries / 20 GeV","",""};

TString tfHists[13] = {"h_FakeRate_tf_minv","h_FakeRate_tf_qt","h_FakeRate_tf_deltaPhi","h_FakeRate_tf_deltaEta","h_FakeRate_tf_deltaR","h_FakeRate_tf_cosThetaStar","h_FakeRate_tf_pt1","h_FakeRate_tf_eta1","h_FakeRate_tf_phi1","h_FakeRate_tf_pt2","h_FakeRate_tf_eta2","h_FakeRate_tf_phi2","h_FakeRate_tf_minv_FineBinning"};

TString ffHists[13] = {"h_FakeRate_ff_minv","h_FakeRate_ff_qt","h_FakeRate_ff_deltaPhi","h_FakeRate_ff_deltaEta","h_FakeRate_ff_deltaR","h_FakeRate_ff_cosThetaStar","h_FakeRate_ff_pt1","h_FakeRate_ff_eta1","h_FakeRate_ff_phi1","h_FakeRate_ff_pt2","h_FakeRate_ff_eta2","h_FakeRate_ff_phi2","h_FakeRate_ff_minv_FineBinning"};

TString ftHists[13] = {"h_FakeRate_ft_minv","h_FakeRate_ft_qt","h_FakeRate_ft_deltaPhi","h_FakeRate_ft_deltaEta","h_FakeRate_ft_deltaR","h_FakeRate_ft_cosThetaStar","h_FakeRate_ft_pt1","h_FakeRate_ft_eta1","h_FakeRate_ft_phi1","h_FakeRate_ft_pt2","h_FakeRate_ft_eta2","h_FakeRate_ft_phi2","h_FakeRate_ft_minv_FineBinning"};


TString GammaJetHists[13] = {"h_GammaJet_minv","h_GammaJet_qt","h_GammaJet_deltaPhi","h_GammaJet_deltaEta","h_GammaJet_deltaR","h_GammaJet_cosThetaStar","h_GammaJet_pt1","h_GammaJet_eta1", "h_GammaJet_phi1","h_GammaJet_pt2","h_GammaJet_eta2","h_GammaJet_phi2","h_GammaJet_minv_FineBinning"};

TString JetJetHists[13] = {"h_JetJet_minv","h_JetJet_qt","h_JetJet_deltaPhi","h_JetJet_deltaEta","h_JetJet_deltaR","h_JetJet_cosThetaStar","h_JetJet_pt1","h_JetJet_eta1", "h_JetJet_phi1","h_JetJet_pt2","h_JetJet_eta2","h_JetJet_phi2","h_JetJet_minv_FineBinning"};

TString nameHists[44] = {"h_Photon1_pt","h_Photon1_pt_log","h_Photon1_eta","h_Photon1_phi","h_Photon2_pt","h_Photon2_pt_log","h_Photon2_eta","h_Photon2_phi","h_Diphoton_Minv","h_Diphoton_Minv_log","h_Diphoton_qt","h_Diphoton_deltaR","h_Diphoton_deltaEta","h_Diphoton_cosThetaStar","h_Diphoton_deltaPhi","h_Vtx_Nvtx","h_Vtx_vx","h_Vtx_vy","h_Vtx_vz","h_Photon1_sigmaIetaIeta","h_Photon1_sigmaEtaEta","h_Photon1_hadOverEm","h_Photon1_trkIsoSumPtHollow04","h_Photon1_trkIsoNtrksHollow04","h_Pho\
ton1_hcalIso04","h_Photon1_ecalIso04","h_Photon1_detEta","h_Photon2_sigmaIetaIeta","h_Photon2_sigmaEtaEta","h_Photon2_hadOverEm","h_Photon2_trkIsoSumPtHollow04","h_Photon2_trkIsoNtrksHollow04","h_Photon2_hcalIso04","h_Photon2_ecalIso04","h_Photon2_detEta","h_Diphoton_qt_log","h_Photon1_hadOverEm_log","h_Photon2_hadOverEm_log","h_Photon2_trkIsoSumPtHollow04_log","h_Photon1_trkIsoSumPtHollow04_log","h_Photon2_hcalIso04_log","h_Photon1_hcalIso04_log","h_Diphoton_Minv_FineBinning","h_Cumalitive_DiphotonMinv"};


TString fNames[12] = {"h_Diphoton_Minv","h_Diphoton_qt","h_Diphoton_deltaPhi","h_Diphoton_deltaEta","h_Diphoton_deltaR","h_Diphoton_cosThetaStar","h_Photon1_pt","h_Photon1_eta","h_Photon1_phi","h_Photon2_pt","h_Photon2_eta","h_Photon2_phi"};

//TString TreeFileLocation = "/afs/cern.ch/user/s/scooper/work/private/data/diPhotonTrees/";
//TString TreeFileLocation = "/afs/cern.ch/work/s/scooper/public/4Otman/DiPhotonTrees/";
TString TreeFileLocation = "/data2/scooper/DiPhotons/Trees/AddMediumLoosePFID/";
//TString HistogramFileLocation = "/afs/cern.ch/user/s/scooper/work/private/results/diPhotonHistogramsPF/";
//TString HistogramFileLocation = "/data2/scooper/DiPhotons/HistogramFiles/AddMediumLoosePFID/";
//TString HistogramFileLocation = "/afs/cern.ch/user/s/scooper/work/public/DiPhotonHistograms/PFID_19p5invFb_looseID_photon1p5pctScaleShifts/";
TString HistogramFileLocation = "/afs/cern.ch/user/s/scooper/work/public/DiPhotonHistograms/PFID_19p5invFb_looseID_photon1p5pctScaleShifts_300minv/";

int inputfiles=4;
double xsec[4]={25.41,1.079e-2,15.53,3.202e-4};
int ngenevents[4]={500254,500038,500050,500352};

//int inputfiles=3;
//double xsec[4]={15.53,3.202e-4,75.39};
//int ngenevents[4]={500050,500352,1154970};

TFile* ftemp[4];

void CreateHistogramFiles(TString Sample = "Diphoton", TString SampleType = "data", TString JSON = "NOJSON")
{
 
  cout << "Entering CreateHistogramFiles method with parameters:" <<endl;
  cout<<"Sample: "<<Sample.Data()
      <<" Sample type: "<<SampleType.Data()
      <<" JSON: "<<JSON.Data()
      <<endl;

  TString outName;
  TFile *outfilename;

  TString  inputfile= TreeFileLocation+Sample+".root";
   
  std::string treePath = SampleType=="signal" ? "diphotonSignalMCAnalyzer/fTree" : "diphotonAnalyzer/fTree";
  TChain *chain_tt = new TChain(treePath.c_str());

  cout << "inputfile: " << inputfile.Data() << endl;
  chain_tt->Add(inputfile.Data());
  cout << "TT entries = " << chain_tt->GetEntries() <<endl;
  
  cout << "CreateHistogramFiles: Start processing main loop entries" <<endl;

  PlottingCodeLoop* treereader = new PlottingCodeLoop(chain_tt);
  treereader->_cutPhoton1pt = 80.;
  treereader->_cutPhoton2pt = 80.;
  treereader->_fakeStatus="TightTight";
  treereader->_JSON=JSON.Data();
  treereader->_SampleType=SampleType.Data();
  outName = TString::Format("histograms_%s.root", Sample.Data());
  outfilename = new TFile(HistogramFileLocation+Sample.Data()+"/"+outName.Data(),"RECREATE");
  treereader->_outputfile = outfilename;
  treereader->Loop();
  
  cout << "CreateHistogramFiles: Finished processing main loop entries" <<endl;


  if (SampleType=="data"){

    cout << "CreateHistogramFiles: Start processing TF, FF, FT trees. This is data." <<endl;

    TChain *chain_tf = new TChain("diphotonAnalyzer/fTightFakeTree");
    chain_tf->Add(inputfile.Data());
    cout << "TF entries = " << chain_tf->GetEntries() <<endl;
    
    TChain *chain_ft = new TChain("diphotonAnalyzer/fFakeTightTree");
    chain_ft->Add(inputfile.Data());
    cout << "FT entries = " << chain_ft->GetEntries() <<endl;

    TChain *chain_ff = new TChain("diphotonAnalyzer/fFakeFakeTree");
    chain_ff->Add(inputfile.Data());
    cout << "FF entries = " << chain_ff->GetEntries() <<endl;

    PlottingCodeLoop* treereaderTF = new PlottingCodeLoop(chain_tf);
    treereaderTF->_fakeStatus = "TightFake";    
    treereaderTF->_cutPhoton1pt = 80.;
    treereaderTF->_cutPhoton2pt = 80.;
    treereader->_JSON=JSON.Data();
    treereader->_SampleType=SampleType.Data();
    outName = TString::Format("histograms_%s_TF.root",Sample.Data());
    outfilename = new TFile(HistogramFileLocation+Sample+"/"+outName,"RECREATE");
    treereaderTF->_outputfile = outfilename;
    treereaderTF->Loop(); 
    
    PlottingCodeLoop* treereaderFT = new PlottingCodeLoop(chain_ft);
    treereaderFT->_fakeStatus = "FakeTight";    
    treereaderFT->_cutPhoton1pt = 80.;
    treereaderFT->_cutPhoton2pt = 80.;
    treereader->_JSON=JSON.Data();
    treereader->_SampleType=SampleType.Data();
    outName = TString::Format("histograms_%s_FT.root",Sample.Data());
    outfilename = new TFile(HistogramFileLocation+Sample+"/"+outName,"RECREATE");
    treereaderFT->_outputfile = outfilename;
    treereaderFT->Loop();     
    
    PlottingCodeLoop* treereaderFF = new PlottingCodeLoop(chain_ff);
    treereaderFF->_fakeStatus = "FakeFake";    
    treereaderFF->_cutPhoton1pt = 80.;
    treereaderFF->_cutPhoton2pt = 80.;
    treereader->_JSON=JSON.Data();
    treereader->_SampleType=SampleType.Data();
    outName = TString::Format("histograms_%s_FF.root",Sample.Data());
    outfilename= new TFile(HistogramFileLocation+Sample+"/"+outName.Data(),"RECREATE");
    treereaderFF->_outputfile = outfilename;
    treereaderFF->Loop();  
     
    cout << "TT entries = " << chain_tt->GetEntries() <<endl;
    cout<<" TF entries = " <<chain_tf->GetEntries()<<endl;
    cout<<" FT entries = " <<chain_ft->GetEntries()<<endl;
    cout<<" FF entries = " <<chain_ff->GetEntries()<<endl;
 
  }  // ends if data

  outfilename->cd();
  outfilename->Close();

} 



void makeplots(TString Sample = "Diphoton",TString lumi = "1", TString JSONFile="NOJSON" ,TString SampleType = "data")
{

  cout<<"Entering makeplots method with parameters: "<<endl; 
  cout<<"Sample: "<<Sample.Data()
      <<" Sample type: "<<SampleType.Data()
      <<" JSON: "<<JSONFile.Data()
      <<" lumi: "<<lumi.Data()<<" /pb"
      <<endl;

  Double_t crosssection = 1.;
  Int_t numevents = 1.;
  Double_t luminumber = lumi.Atof();

  char canvasname[10];
  TCanvas*  c[42];
  TString histogramfileinput = HistogramFileLocation+Sample+"/histograms_"+Sample+".root";

  cout<<"Opening histogram root file "<<histogramfileinput.Data()<<endl;
  TFile* fhists = TFile::Open(histogramfileinput.Data());
  fhists->cd();

  cout<<"Successfully Opened histogram root file "<<histogramfileinput.Data()<<endl;

  TPaveText *LumiLabel = new TPaveText(.6,.7,.8,.85,"NDC");
  LumiLabel->SetTextSize(0.03);
  LumiLabel->SetFillStyle(0);
  LumiLabel->SetBorderSize(0);
  LumiLabel->AddText("CMS Internal");

  if (SampleType=="data" || SampleType=="datamc" || SampleType=="stitchdatamc"){
    TString pico(" pb^{-1}");
    TString Luminosity = lumi.Append(pico);
    LumiLabel->AddText(Luminosity.Data());
    LumiLabel->AddText(JSONFile.Data());
  }

  cout<<"Plotting histograms "<<endl;

  //   if(Sample.Contains(("diphoton_tree_DiPhotonBorn_Pt25To250_Summer12_Sept10th"))) {crosssection = xsec[0];numevents = ngenevents[0];}
  //   if(Sample.Contains(("diphoton_tree_DiPhotonBorn_Pt250ToInf_Summer12_Sept10th"))) {crosssection = xsec[1];numevents = ngenevents[1];}
  //   if(Sample.Contains(("diphoton_tree_DiPhotonBox_Pt25To250_Summer12_Sept10th"))) {crosssection = xsec[2];numevents = ngenevents[2];}
  //   if(Sample.Contains(("diphoton_tree_DiPhotonBox_Pt250ToInf_Summer12_Sept10th"))) {crosssection = xsec[3];numevents = ngenevents[3];}

  gStyle->SetOptStat("ourmei")  ;

  TH1F* histos[42];
  for(int i=0;i<42;i++)
    {
      sprintf(canvasname,"c%i",i);
      c[i] = new TCanvas(canvasname, canvasname, 800., 600.);
      cout<<"Canvas "<<canvasname<<endl;
      c[i]->cd();
      histos[i] = (TH1F*)fhists->Get(nameHists[i].Data());
      histos[i]->Sumw2();
      cout<<"Histogram "<<histos[i]->GetName()<<endl;

      //histos[i]->Scale((luminumber*crosssection)/numevents);

      histos[i]->GetXaxis()->CenterTitle();
      histos[i]->GetYaxis()->CenterTitle();
      histos[i]->SetMarkerStyle(20);
      histos[i]->SetMarkerColor(1);

      histos[i]->SetMinimum(0.);

      if (nameHists[i].Contains("_log"))
	{
	  histos[i]->SetMinimum(1.e-2);
	  c[i]->SetLogy(1);
	}

      histos[i]->Draw("HIST EP");

      LumiLabel->Draw();
      c[i]->SaveAs(nameHists[i]+"_"+Sample+".C");
      c[i]->SaveAs(nameHists[i]+"_"+Sample+".png");
      c[i]->SaveAs(nameHists[i]+"_"+Sample+".pdf");
      c[i]->SaveAs(nameHists[i]+"_"+Sample+".eps");

      c[i]->cd();
      c[i]->Close();
    }//end of loop over all histograms

  cout<<"All histograms plotted and saved()"<<endl;

  fhists->cd();
  fhists->Close();

}//end of makeplots method






void MakeCombinedMCHistos()
{

  cout << "Entering MakeCombinedMCHistos method" <<endl;


  //FIXME-------------------------
  //DANGEROUS HERE
  //TRY NOT TO HARD CODE
  //THE ROOT FILES
  //------------------------------


  gROOT->SetStyle("Plain");
  gStyle->SetOptStat("ourme");

  //--------------------------------
  //ALWAYS HAS TO BE SCALED TO 1 pb-1

  const double lumi = 1.0; // units of pb-1
  //const

  TString inputmcfiles[inputfiles];

  //     inputmcfiles[0]= HistogramFileLocation+"diphoton_tree_DiPhotonBorn_Pt25To250_Summer12_Sept10th/histograms_diphoton_tree_DiPhotonBorn_Pt25To250_Summer12_Sept10th.root";
  //     inputmcfiles[1]= HistogramFileLocation+"diphoton_tree_DiPhotonBorn_Pt250ToInf_Summer12_Sept10th/histograms_diphoton_tree_DiPhotonBorn_Pt250ToInf_Summer12_Sept10th.root";
  //     inputmcfiles[2]= HistogramFileLocation+"diphoton_tree_DiPhotonBox_Pt25To250_Summer12_Sept10th/histograms_diphoton_tree_DiPhotonBox_Pt25To250_Summer12_Sept10th.root";	
  //     inputmcfiles[3]= HistogramFileLocation+"diphoton_tree_DiPhotonBox_Pt250ToInf_Summer12_Sept10th/histograms_diphoton_tree_DiPhotonBox_Pt250ToInf_Summer12_Sept10th.root";

//     inputmcfiles[0]= HistogramFileLocation+"ExoDiPhotonAnalyzer_PFDec14th_53X_BornPt25To250/histograms_ExoDiPhotonAnalyzer_PFDec14th_53X_BornPt25To250.root";
//     inputmcfiles[1]= HistogramFileLocation+"ExoDiPhotonAnalyzer_PFDec14th_53X_BornPt250ToInf/histograms_ExoDiPhotonAnalyzer_PFDec14th_53X_BornPt250ToInf.root";
//     inputmcfiles[2]= HistogramFileLocation+"ExoDiPhotonAnalyzer_PFDec14th_53X_BoxPt25To250/histograms_ExoDiPhotonAnalyzer_PFDec14th_53X_BoxPt25To250.root";
//     inputmcfiles[3]= HistogramFileLocation+"ExoDiPhotonAnalyzer_PFDec14th_53X_BoxPt250ToInf/histograms_ExoDiPhotonAnalyzer_PFDec14th_53X_BoxPt250ToInf.root";
  
//   inputmcfiles[0]= HistogramFileLocation+"ExoDiPhotonAnalyzer_DETDec14th_53X_BornPt25To250/histograms_ExoDiPhotonAnalyzer_DETDec14th_53X_BornPt25To250.root";
//   inputmcfiles[1]= HistogramFileLocation+"ExoDiPhotonAnalyzer_DETDec14th_53X_BornPt250ToInf/histograms_ExoDiPhotonAnalyzer_DETDec14th_53X_BornPt250ToInf.root";
//   inputmcfiles[2]= HistogramFileLocation+"ExoDiPhotonAnalyzer_DETDec14th_53X_BoxPt25To250/histograms_ExoDiPhotonAnalyzer_DETDec14th_53X_BoxPt25To250.root";
//   inputmcfiles[3]= HistogramFileLocation+"ExoDiPhotonAnalyzer_DETDec14th_53X_BoxPt250ToInf/histograms_ExoDiPhotonAnalyzer_DETDec14th_53X_BoxPt250ToInf.root";
  
//   inputmcfiles[0]= HistogramFileLocation+"ExoDiPhotonAnalyzer_PFAug10th_52X_BornPt25To250/histograms_ExoDiPhotonAnalyzer_PFAug10th_52X_BornPt25To250.root";
//   inputmcfiles[1]= HistogramFileLocation+"ExoDiPhotonAnalyzer_PFAug10th_52X_BornPt250ToInf/histograms_ExoDiPhotonAnalyzer_PFAug10th_52X_BornPt250ToInf.root";
//   inputmcfiles[2]= HistogramFileLocation+"ExoDiPhotonAnalyzer_PFAug10th_52X_BoxPt25To250/histograms_ExoDiPhotonAnalyzer_PFAug10th_52X_BoxPt25To250.root";
//   inputmcfiles[3]= HistogramFileLocation+"ExoDiPhotonAnalyzer_PFAug10th_52X_BoxPt250ToInf/histograms_ExoDiPhotonAnalyzer_PFAug10th_52X_BoxPt250ToInf.root";
  
    inputmcfiles[0]= HistogramFileLocation+"ExoDiPhotonAnalyzer_PFDec14th_53X_BornPt25To250/histograms_ExoDiPhotonAnalyzer_PFDec14th_53X_BornPt25To250.root";
    inputmcfiles[1]= HistogramFileLocation+"ExoDiPhotonAnalyzer_PFDec14th_53X_BornPt250ToInf/histograms_ExoDiPhotonAnalyzer_PFDec14th_53X_BornPt250ToInf.root";
    inputmcfiles[2]= HistogramFileLocation+"ExoDiPhotonAnalyzer_PFDec14th_53X_BoxPt25To250/histograms_ExoDiPhotonAnalyzer_PFDec14th_53X_BoxPt25To250.root";
    inputmcfiles[3]= HistogramFileLocation+"ExoDiPhotonAnalyzer_PFDec14th_53X_BoxPt250ToInf/histograms_ExoDiPhotonAnalyzer_PFDec14th_53X_BoxPt250ToInf.root";
    //inputmcfiles[2]= HistogramFileLocation+"ExoDiPhotonAnalyzer_PFDec14th_53X_DiPhotonJetsMadGraph/histograms_ExoDiPhotonAnalyzer_PFDec14th_53X_DiPhotonJetsMadGraph.root";
  


  cout << "Input MC ROOT files: " <<endl;

  //------------------
  //ABSOLUTELY USELESS
  //ALREADY HAVE INPUTMCFILES
  //------------------

  TString filenames[inputfiles];
  for(int ifile=0;ifile<inputfiles;ifile++) {
    filenames[ifile] = inputmcfiles[ifile];
    cout << filenames[ifile]<<endl;
  }

  cout << "Opening MC ROOT files: " <<endl;

  for(int ifile=0;ifile<inputfiles;ifile++) {
    ftemp[ifile] = TFile::Open(filenames[ifile].Data());
  }

  int nHists = 43;
  TH1F* indHistos[nHists][inputfiles];
  Double_t totalError = 0.;

  for(int ifile=0;ifile<inputfiles;ifile++) {

    cout << "file "<<filenames[ifile].Data()<<endl;
    cout << "xsec " << xsec[ifile] << " #events " << ngenevents[ifile] << " lumi " << lumi << endl;

    for (int i=0; i<nHists; i++) {
      indHistos[i][ifile] = (TH1F*)ftemp[ifile]->Get(nameHists[i].Data());
      indHistos[i][ifile]->Sumw2();

      //if(i == 29) cout<<indHistos[i][ifile]->GetName()<<" entries "<<indHistos[i][ifile]->GetEntries()<<" integral "<<indHistos[i][ifile]->Integral()<<endl;

      //       //------------------------------
      //       //WHAT IS THE POINTS OF ALL THIS
      //       if (i==1) {
      // 	//  error N_tot = sqrt(C1^2 * n1 + C2^2 * n2 + ...)
      // 	Double_t ci = xsec[ifile]/ngenevents[ifile];
      // 	//      totalError = totalError + ci*ci*(indHistos[i][ifile]->Integral());
      // 	totalError = totalError + ci*TMath::Sqrt(indHistos[i][ifile]->Integral());
      //       }
      //       //------------------------------

      indHistos[i][ifile]->Scale((lumi/ngenevents[ifile])*xsec[ifile]);

      //if(i == 29) cout<<indHistos[i][ifile]->GetName()<<" entries "<<indHistos[i][ifile]->GetEntries()<<" integral "<<indHistos[i][ifile]->Integral()<<endl;

    }//end of loop over all histograms
  }//end of loop over all MC background files

  TH1F* allHistos[nHists];

  // instantiate the sum histos with histo in first file
  for (int i=0; i<nHists; i++) {
    allHistos[i] = indHistos[i][0];
    if(i == 29) cout<<allHistos[i]->GetName()<<" entries "<<allHistos[i]->GetEntries()<<" integral "<<allHistos[i]->Integral()<<endl;
  }

  // note we exclude first file (we start at ifile=1) since we used it to instantiate the sum histo
  for(int ifile=1;ifile<inputfiles;ifile++) {
    for (int i=0; i<nHists; i++) {
      allHistos[i]->Add(indHistos[i][ifile]);
      //if(i == 29) cout<<allHistos[i]->GetName()<<" entries "<<allHistos[i]->GetEntries()<<" integral "<<allHistos[i]->Integral()<<endl;
    }
  }

  //------------------------------
  //AGAIN HARDCODED FOR MC TREE ALL
  //------------------------------
  TString MCALLOut = HistogramFileLocation+"diphoton_tree_MC_all/";
  TFile *fout = TFile::Open(TString(MCALLOut+"histograms_diphoton_tree_MC_all.root").Data(),"RECREATE");
  fout->cd();

  for (int i=0; i<nHists; i++) {
    allHistos[i]->Write();
  }
  fout->Close();

  //makeplots("diphoton_tree_MC_all","1","MC","mc");
  //makeplots("diphoton_tree_MC_all","10252","MC","mc");

}//enf of method MakeCombinedMCHistos


TH1F* MakeChists(TH1F* hist){
  TH1F* cHist = (TH1F*) hist->Clone("cHist");
  cHist->SetDirectory(0);

  int maxBin = hist->GetNbinsX()+1;
  for(int binNr=0;binNr<=hist->GetNbinsX();binNr++){
    float nrEntries =hist->Integral(binNr,maxBin);
    cHist->SetBinContent(binNr,nrEntries);
  }
  return cHist;
}



void OverlayMCandData( TString Sample  = "Diphoton", TString lumi = "1", TString JSONFile = "NOJSON")
{

  cout<<"Entering OverlayMCandData method with parameters: "<<endl; 
  cout<<"Sample: "<<Sample.Data()
      <<" JSON: "<<JSONFile.Data()
      <<" lumi: "<<lumi.Data()
      <<" /pb"
      <<endl;

  //---------------------
  //AGAIN HARDCODED FOR MC 
  //---------------------
  TString histogramdata =HistogramFileLocation+Sample+"/histograms_"+Sample+".root";
  TString histogrammc=HistogramFileLocation+"diphoton_tree_MC_all/histograms_diphoton_tree_MC_all.root";
  TFile* fmchists = TFile::Open(histogrammc.Data());
  TFile* fdatahists = TFile::Open(histogramdata.Data());

  char canvasname[10];
  TCanvas*  c[42];

  TString pico(" pb^{-1}");
  TPaveText *LumiLabel = new TPaveText(.55,.75,.75,.9,"NDC");
  LumiLabel->SetTextSize(0.03);
  LumiLabel->SetFillStyle(0);
  LumiLabel->SetBorderSize(0);
  LumiLabel->AddText("CMS Internal");

  TString Luminosity = lumi.Append(pico);
  LumiLabel->AddText(Luminosity.Data());
  LumiLabel->AddText(JSONFile.Data());

  TH1F* histosdata[42];
  TH1F* histosmc[42];
  TLegend* DataMCLegend[42];
  float lumiNumber = lumi.Atof();
  cout<<"Lumi number "<<lumiNumber<<endl;

  gStyle->SetOptStat("ourmei");

  for(int i=0;i<42;i++)
    {

      sprintf(canvasname,"c%i",i);
      c[i] = new TCanvas(canvasname, canvasname, 800., 600.);
      cout<<"Canvas "<<canvasname<<endl;
      c[i]->cd();

      histosdata[i] = (TH1F*)fdatahists->Get(nameHists[i].Data());
      histosdata[i]->SetMinimum(.01);
      cout<<"Data histogram name "<<histosdata[i]->GetName()<<endl;

      histosmc[i] = (TH1F*)fmchists->Get(nameHists[i].Data());
      histosmc[i]->Scale(lumiNumber);
      histosmc[i]->SetMinimum(.01);
      cout<<"MC histogram name "<<histosmc[i]->GetName()<<endl;

      DataMCLegend[i]= new TLegend(0.70,0.6,0.88,0.75,"","NDC");
      DataMCLegend[i]->SetFillStyle(0);
      DataMCLegend[i]->SetTextSize(0.03);

      histosmc[i]->GetXaxis()->CenterTitle();
      histosmc[i]->GetYaxis()->CenterTitle();
      histosmc[i]->SetFillColor(kBlue);
      histosmc[i]->SetMaximum(1.3*histosdata[i]->GetMaximum());

      histosdata[i]->SetMarkerColor(1);
      histosdata[i]->SetMaximum(1.3*histosdata[i]->GetMaximum());
      //histosdata[i]->SetMinimum(.8*histosmc[i]->GetMinimum());
      //histosdata[i]->Draw("HIST EP");

      //histosmc[i]->Draw("SAME HIST");
      histosmc[i]->Draw("HIST");
      histosmc[i]->Draw("sameaxis");
      //histosdata[i]->Draw("HIST SAME EP");
      histosdata[i]->Draw("EP,SAMES");

      DataMCLegend[i]->AddEntry(histosdata[i],"Data","LEP");
      DataMCLegend[i]->AddEntry(histosmc[i],"SM Diphoton","F");
      //gPad->RedrawAxis();
      LumiLabel->Draw();
      DataMCLegend[i]->Draw();

      if (nameHists[i].Contains("_log"))
        {
	  c[i]->SetLogy(1);
	}

      c[i]->SaveAs(nameHists[i]+"_"+Sample.Data()+"Overlay"+".C");
      c[i]->SaveAs(nameHists[i]+"_"+Sample.Data()+"Overlay"+".png");
      c[i]->SaveAs(nameHists[i]+"_"+Sample.Data()+"Overlay"+".pdf");
      c[i]->SaveAs(nameHists[i]+"_"+Sample.Data()+"Overlay"+".eps");
    }//end of loop over histograms


  //AGAIN HERE INV MASS INDEX 9 IS HARDCODED
  //SHOULD CHANGE THAT
  //ALSO UNDERSTAND WHY CUMUL PLOT HAS NO SAME NUMBER OF ENTRIES

  TCanvas* CumulativeCanvas = new TCanvas("CumulativeCanvas", "CumulativeCanvas", 800., 600.);
  CumulativeCanvas->cd();
  CumulativeCanvas->SetLogy(1);
  TH1F* h_Diphoton_Minv_log_DataCumm =  MakeChists(histosdata[9]);
  TH1F* h_Diphoton_Minv_log_MCCumm = MakeChists(histosmc[9]);

  h_Diphoton_Minv_log_MCCumm->SetFillColor(kBlue);

  TPaveText *LumiLabelCumm = new TPaveText(.55,.75,.75,.9,"NDC");
  LumiLabelCumm->SetTextSize(0.03);
  LumiLabelCumm->SetFillStyle(0);
  LumiLabelCumm->SetBorderSize(0);
  LumiLabelCumm->AddText("CMS Internal");
  LumiLabelCumm->AddText(TString(lumi+" pb^{-1}").Data());
  LumiLabelCumm->AddText(JSONFile.Data());

  TLegend* DataMCLegendCumm = new TLegend(0.70,0.6,0.88,0.75,"","NDC");
  DataMCLegendCumm->AddEntry( h_Diphoton_Minv_log_DataCumm,"Data","LEP");
  DataMCLegendCumm->AddEntry(h_Diphoton_Minv_log_MCCumm,"SM Diphoton","F");
  DataMCLegendCumm->SetFillColor(0);


  h_Diphoton_Minv_log_DataCumm->SetMaximum(2.5*h_Diphoton_Minv_log_DataCumm->GetMaximum());
  h_Diphoton_Minv_log_MCCumm->SetMaximum(2.5*h_Diphoton_Minv_log_DataCumm->GetMaximum());

  h_Diphoton_Minv_log_DataCumm->SetMinimum(.01);
  h_Diphoton_Minv_log_MCCumm->SetMinimum(.01);

  h_Diphoton_Minv_log_MCCumm->Draw("HIST");
  h_Diphoton_Minv_log_MCCumm->Draw("SAMEAXIS");
  h_Diphoton_Minv_log_DataCumm->Draw("EP,SAME");

  //gPad->RedrawAxis();
  LumiLabelCumm->Draw();
  DataMCLegendCumm->Draw();
  //CumulativeCanvas->SaveAs(nameHists[43]+Sample.Data()+"Overlay"+".png");
  CumulativeCanvas->SaveAs(nameHists[42]+"_"+Sample.Data()+"Overlay"+".png");
  CumulativeCanvas->SaveAs(nameHists[42]+"_"+Sample.Data()+"Overlay"+".eps");
  CumulativeCanvas->SaveAs(nameHists[42]+"_"+Sample.Data()+"Overlay"+".pdf");
  CumulativeCanvas->SaveAs(nameHists[42]+"_"+Sample.Data()+"Overlay"+".C");

  //   //Closing all files

  //   fmchists->cd();
  //   fmchists->Close();

  //   fdatahists->cd();
  //   fdatahists->Close();


}//end of method OverlayMCAndData



void StitchBackgroundandMC(TString Sample = "Diphoton", TString lumi = "1", TString JSONFile = "NOJSON")
{

  cout<<"Entering StitchBackgroundandMC method with parameters: "<<endl; 
  cout<<"Sample: "<<Sample.Data()
      <<" JSON: "<<JSONFile.Data()
      <<" lumi: "<<lumi.Data()
      <<" /pb"
      <<endl;

  TString histogramdata = HistogramFileLocation+Sample+"/histograms_"+Sample+".root";
  TFile* fdatahists = TFile::Open(histogramdata.Data());
  TH1F* histosmcdminv[42][4];
  TH1F* histosdatadminv[42];

  char canvasname[10];
  TCanvas* c[42];

  TPaveText *LumiLabel = new TPaveText(.6,.7,.8,.85,"NDC");
  LumiLabel->SetTextSize(0.03);
  LumiLabel->SetFillStyle(0);
  LumiLabel->SetBorderSize(0);
  LumiLabel->AddText("CMS Internal");

  TString pico(" pb^{-1}");
  TString Luminosity = lumi.Append(pico);
  LumiLabel->AddText(Luminosity.Data());
  LumiLabel->AddText(JSONFile.Data());
  float lumiNumber = lumi.Atof();
  cout<<"Lumi number "<<lumiNumber<<endl;

  Color_t fillColors[inputfiles];
  fillColors[0] = kGreen;
  fillColors[1] = kRed;
  fillColors[2] = kBlue;
  fillColors[3] = 6;

  //AGAIN HERE HARDCODED

  TString inputmcfiles[inputfiles];
  inputmcfiles[0]= HistogramFileLocation+"diphoton_tree_DiPhotonBorn_Pt25To250_Summer12_Sept10th/histograms_diphoton_tree_DiPhotonBorn_Pt25To250_Summer12_Sept10th.root";
  inputmcfiles[1]= HistogramFileLocation+"diphoton_tree_DiPhotonBorn_Pt250ToInf_Summer12_Sept10th/histograms_diphoton_tree_DiPhotonBorn_Pt250ToInf_Summer12_Sept10th.root";
  inputmcfiles[2]= HistogramFileLocation+"diphoton_tree_DiPhotonBox_Pt25To250_Summer12_Sept10th/histograms_diphoton_tree_DiPhotonBox_Pt25To250_Summer12_Sept10th.root";
  inputmcfiles[3]= HistogramFileLocation+"diphoton_tree_DiPhotonBox_Pt250ToInf_Summer12_Sept10th/histograms_diphotn_tree_DiPhotonBox_Pt250ToInf_Summer12_Sept10th.root";

  TString filenames[inputfiles];
  for(int ifile=0;ifile<inputfiles;ifile++) {
    filenames[ifile] = inputmcfiles[ifile];
  }

  TString SampleNames[4]={"Born_Pt_25To250","Born_Pt_250ToInf","Box_Pt_25To250","Box_Pt_250ToInf"};
  TLegend *StackLegend[42];
  THStack *StackMC[42];

  for(int i=0;i<42;i++){

    sprintf(canvasname,"c%i",i);
    c[i] = new TCanvas(canvasname, canvasname, 800., 600.);
    cout<<"Canvas "<<canvasname<<endl;
    c[i]->cd();

    StackLegend[i]= new TLegend(0.70,0.6,0.88,0.75,"","NDC");
    StackLegend[i]->SetFillStyle(0);
    StackLegend[i]->SetTextSize(0.03);

    TString StackName(nameHists[i]+"_Stack");
    StackMC[i]= new THStack(StackName.Data(),StackName.Data());
    histosdatadminv[i]=(TH1F*)fdatahists->Get(nameHists[i].Data());

    for(int ifile=0;ifile<inputfiles;ifile++) {
      ftemp[ifile] = TFile::Open(filenames[ifile].Data());
      histosmcdminv[ifile][i]= (TH1F*)ftemp[ifile]->Get(nameHists[i].Data());
      histosmcdminv[ifile][i]->SetFillColor(fillColors[ifile]);
      histosmcdminv[ifile][i]->Scale((lumiNumber/ngenevents[ifile])*xsec[ifile]);

      StackMC[i]->Add(histosmcdminv[ifile][i]);

      if (nameHists[i].Contains("_log"))
	{
	  histosmcdminv[ifile][i]->SetMinimum(0.01);
	  c[i]->SetLogy(1);
	}

      //histosdatadminv[i]->SetMinimum(.8*(StackMC[i]->GetMinimum()));
      StackLegend[i]->AddEntry(histosdatadminv[i],"Data","LEP");
      StackLegend[i]->AddEntry(histosmcdminv[ifile][i],SampleNames[ifile].Data(),"f");

      StackMC[i]->Draw("HIST") ;
      histosdatadminv[i]->Draw("EP,SAME");

      LumiLabel->Draw();
      StackLegend[i]->Draw();
      gPad->RedrawAxis();

      c[i]->SaveAs(nameHists[i]+"_"+Sample.Data()+"StitchOverlay"+".C");
      c[i]->SaveAs(nameHists[i]+"_"+Sample.Data()+"StitchOverlay"+".png");
      c[i]->SaveAs(nameHists[i]+"_"+Sample.Data()+"StitchOverlay"+".pdf");
      c[i]->SaveAs(nameHists[i]+"_"+Sample.Data()+"StitchOverlay"+".eps");

    }//end of loop over all MC samples
  }//end of loop over all histograms

  fdatahists->cd();
  fdatahists->Close();

  for(int ifile=0;ifile<inputfiles;ifile++) {
    ftemp[ifile]->cd();
    ftemp[ifile]->Close();
  }

}//end of method StitchBackgroundandMC


void fakeratehistos(TString Sample = "Diphoton",TString lumi = "1", TString JSONFile="NOJSON" )
{
  TH1F* h_JetJet_pt1;
  TH1F* h_JetJet_eta1;
  TH1F* h_JetJet_phi1;
  TH1F* h_JetJet_pt2;
  TH1F* h_JetJet_eta2;
  TH1F* h_JetJet_phi2;
  TH1F* h_JetJet_minv;
  TH1F* h_JetJet_minv_FineBinning;
  TH1F* h_JetJet_qt;
  TH1F* h_JetJet_deltaPhi;
  TH1F* h_JetJet_deltaEta;
  TH1F* h_JetJet_deltaR;
  TH1F* h_JetJet_cosThetaStar;
  TH1F* h_JetJet_minv_UpperBand;
  TH1F* h_JetJet_minv_LowerBand;
  TH1F* h_JetJet_minv_UpperError;
  TH1F* h_JetJet_minv_LowerError;

  TH1F* h_GammaJet_pt1;
  TH1F* h_GammaJet_eta1;
  TH1F* h_GammaJet_phi1;
  TH1F* h_GammaJet_pt2;
  TH1F* h_GammaJet_eta2;
  TH1F* h_GammaJet_phi2;
  TH1F* h_GammaJet_minv;
  TH1F* h_GammaJet_minv_FineBinning;
  TH1F* h_GammaJet_qt;
  TH1F* h_GammaJet_deltaPhi;
  TH1F* h_GammaJet_deltaEta;
  TH1F* h_GammaJet_deltaR;
  TH1F* h_GammaJet_cosThetaStar;
  TH1F* h_GammaJet_minv_UpperBand;
  TH1F* h_GammaJet_minv_LowerBand;
  TH1F* h_GammaJet_minv_UpperError;
  TH1F* h_GammaJet_minv_LowerError;

  h_GammaJet_pt1  = new TH1F("h_GammaJet_pt1","#gamma_{1} p_{T};#gamma_{1} p_{T}",42,60.,900.);
  h_GammaJet_eta1 = new TH1F("h_GammaJet_eta1","#gamma_{1} #eta;#gamma_{1} #eta",60,-3.,3.);
  h_GammaJet_phi1 = new TH1F("h_GammaJet_phi1","#gamma_{1} #phi;#gamma_{1} #phi",36,-3.14159,3.14159);
  h_GammaJet_pt2  = new TH1F("h_GammaJet_pt2","#gamma_{2} p_{T};#gamma_{2} p_{T}",42,60.,900.);
  h_GammaJet_eta2 = new TH1F("h_GammaJet_eta2","#gamma_{2} #eta;#gamma_{2} #eta",60,-3.,3.);
  h_GammaJet_phi2 = new TH1F("h_GammaJet_phi2","#gamma_{1} #phi;#gamma_{1} #phi",36,-3.14159,3.14159);
  //   h_GammaJet_minv         = new TH1F("h_GammaJet_minv",        "Diphoton Invariant Mass;M_{#gamma#gamma} [GeV]",43,140,1000);
  h_GammaJet_minv         = new TH1F("h_GammaJet_minv",        "Diphoton Invariant Mass;M_{#gamma#gamma} [GeV]",199,20.,4000.);
  h_GammaJet_minv_FineBinning         = new TH1F("h_GammaJet_minv_FineBinning",        "Diphoton Invariant Mass;M_{#gamma#gamma} [GeV]",3980,20.,4000.);
  h_GammaJet_qt           = new TH1F("h_GammaJet_qt ",         "Diphoton qt;#gamma#gamma qt [GeV]",50,0.,600.);
  h_GammaJet_deltaPhi     = new TH1F("h_GammaJet_deltaPhi",    "Diphoton #Delta#phi;#gamma#gamma #Delta#phi",36,-3.14159,3.14159);
  h_GammaJet_deltaEta     = new TH1F("h_GammaJet_deltaEta",    "Diphoton #Delta#eta;#gamma#gamma #Delta#eta",60,-6.,6.);
  h_GammaJet_deltaR       = new TH1F("h_GammaJet_deltaR",      "Diphoton #DeltaR; #gamma#gamma #DeltaR",70,0.,7.);
  h_GammaJet_cosThetaStar = new TH1F("h_GammaJet_cosThetaStar","Diphoton |cos(#theta *)|; #gamma#gamma |cos#theta*|",20,0.,1.);
  h_GammaJet_minv_UpperBand         = new TH1F("h_GammaJet_minv_UpperBand",        "Diphoton Invariant Mass Upper Syst. Band;M_{#gamma#gamma} [GeV]",3980,20.,4000.);
  h_GammaJet_minv_LowerBand         = new TH1F("h_GammaJet_minv_LowerBand",        "Diphoton Invariant Mass Lower Syst. Band;M_{#gamma#gamma} [GeV]",3980,20.,4000.);
  h_GammaJet_minv_UpperError         = new TH1F("h_GammaJet_minv_UpperError",        "Diphoton Invariant Mass Upper Syst. Error;M_{#gamma#gamma} [GeV]",3980,20.,4000.);
  h_GammaJet_minv_LowerError         = new TH1F("h_GammaJet_minv_LowerError",        "Diphoton Invariant Mass Lower Syst. Error;M_{#gamma#gamma} [GeV]",3980,20.,4000.);


  h_JetJet_pt1  = new TH1F("h_JetJet_pt1","#gamma_{1} p_{T};#gamma_{1} p_{T}",42,60.,900.);
  h_JetJet_eta1 = new TH1F("h_JetJet_eta1","#gamma_{1} #eta;#gamma_{1} #eta",60,-3.,3.);
  h_JetJet_phi1 = new TH1F("h_JetJet_phi1","#gamma_{1} #phi;#gamma_{1} #phi",36,-3.14159,3.14159);
  h_JetJet_pt2  = new TH1F("h_JetJet_pt2","#gamma_{2} p_{T};#gamma_{2} p_{T}",42,60.,900.);
  h_JetJet_eta2 = new TH1F("h_JetJet_eta2","#gamma_{2} #eta;#gamma_{2} #eta",60,-3.,3.);
  h_JetJet_phi2 = new TH1F("h_JetJet_phi2","#gamma_{1} #phi;#gamma_{1} #phi",36,-3.14159,3.14159);
  //   h_JetJet_minv         = new TH1F("h_JetJet_minv",        "Diphoton Invariant Mass;M_{#gamma#gamma} [GeV]",43,140,1000);
  h_JetJet_minv         = new TH1F("h_JetJet_minv",        "Diphoton Invariant Mass;M_{#gamma#gamma} [GeV]",199,20.,4000.);
  h_JetJet_minv_FineBinning         = new TH1F("h_JetJet_minv_FineBinning",        "Diphoton Invariant Mass;M_{#gamma#gamma} [GeV]",3980,20.,4000.);
  h_JetJet_qt           = new TH1F("h_JetJet_qt ",         "Diphoton qt;#gamma#gamma qt [GeV]",50,0.,600.);
  h_JetJet_deltaPhi     = new TH1F("h_JetJet_deltaPhi",    "Diphoton #Delta#phi;#gamma#gamma #Delta#phi",36,-3.14159,3.14159);
  h_JetJet_deltaEta     = new TH1F("h_JetJet_deltaEta",    "Diphoton #Delta#eta;#gamma#gamma #Delta#eta",60,-6.,6.);
  h_JetJet_deltaR       = new TH1F("h_JetJet_deltaR",      "Diphoton #DeltaR; #gamma#gamma #DeltaR",70,0.,7.);
  h_JetJet_cosThetaStar = new TH1F("h_JetJet_cosThetaStar","Diphoton |cos(#theta *)|; #gamma#gamma |cos#theta*|",20,0.,1.);
  h_JetJet_minv_UpperBand         = new TH1F("h_JetJet_minv_UpperBand",        "Diphoton Invariant Mass Upper Syst. Band;M_{#gamma#gamma} [GeV]",3980,20.,4000.);
  h_JetJet_minv_LowerBand         = new TH1F("h_JetJet_minv_LowerBand",        "Diphoton Invariant Mass Lower Syst. Band;M_{#gamma#gamma} [GeV]",3980,20.,4000.);
  h_JetJet_minv_UpperError         = new TH1F("h_JetJet_minv_UpperError",        "Diphoton Invariant Mass Upper Syst. Error;M_{#gamma#gamma} [GeV]",3980,20.,4000.);
  h_JetJet_minv_LowerError         = new TH1F("h_JetJet_minv_LowerError",        "Diphoton Invariant Mass Lower Syst. Error;M_{#gamma#gamma} [GeV]",3980,20.,4000.);

  TString histoTFlocation=HistogramFileLocation.Data()+Sample+"/histograms_"+Sample+"_TF.root";
  TString histoFFlocation=HistogramFileLocation.Data()+Sample+"/histograms_"+Sample+"_FF.root";
  TString histoFTlocation=HistogramFileLocation.Data()+Sample+"/histograms_"+Sample+"_FT.root";
    
  cout<<"Tight Fake Histogram Location: "<<histoTFlocation<<endl;
  cout<<"Fake Tight Histogram Location: "<<histoFTlocation<<endl;
  cout<<"Tight Tight Histogram Location: "<<histoFFlocation<<endl;

  TFile* histoFF = TFile::Open(histoFFlocation.Data());
  TFile* histoFT = TFile::Open(histoFTlocation.Data()); 
  TFile* histoTF = TFile::Open(histoTFlocation.Data());  

  cout<<"Opened histogram root files"<<endl;
    
  TString GammaJetoutName = TString::Format("histograms_%s_GammaJet.root",Sample.Data());
  TFile* histoFileGammaJet = new TFile(HistogramFileLocation.Data()+Sample+"/"+GammaJetoutName,"RECREATE");       
    
  TString JetJetoutName = TString::Format("histograms_%s_JetJet.root",Sample.Data());
  TFile* histoFileJetJet = new TFile(HistogramFileLocation.Data()+Sample+"/"+JetJetoutName,"RECREATE");

  int nFakeHists=13;

  TH1F* hFF[13];
  TH1F* hFT[13];
  TH1F* hTF[13];
  TH1F* histosGammaJet[13];
  TH1F* histosJetJet[13];
    
  TList* GammaJetHistList = new TList;    
  TList* JetJetHistList = new TList;

  GammaJetHistList->Add(h_GammaJet_minv);
  GammaJetHistList->Add(h_GammaJet_qt);
  GammaJetHistList->Add(h_GammaJet_deltaPhi);
  GammaJetHistList->Add(h_GammaJet_deltaEta);
  GammaJetHistList->Add(h_GammaJet_deltaR);
  GammaJetHistList->Add(h_GammaJet_cosThetaStar);
  GammaJetHistList->Add(h_GammaJet_pt1);
  GammaJetHistList->Add(h_GammaJet_eta1);
  GammaJetHistList->Add(h_GammaJet_phi1);
  GammaJetHistList->Add(h_GammaJet_pt2);
  GammaJetHistList->Add(h_GammaJet_eta2);
  GammaJetHistList->Add(h_GammaJet_phi2);
  GammaJetHistList->Add(h_GammaJet_minv_FineBinning);

  JetJetHistList->Add(h_JetJet_minv);
  JetJetHistList->Add(h_JetJet_qt);
  JetJetHistList->Add(h_JetJet_deltaPhi);
  JetJetHistList->Add(h_JetJet_deltaEta);
  JetJetHistList->Add(h_JetJet_deltaR);
  JetJetHistList->Add(h_JetJet_cosThetaStar);
  JetJetHistList->Add(h_JetJet_pt1);
  JetJetHistList->Add(h_JetJet_eta1);
  JetJetHistList->Add(h_JetJet_phi1);
  JetJetHistList->Add(h_JetJet_pt2);
  JetJetHistList->Add(h_JetJet_eta2);
  JetJetHistList->Add(h_JetJet_phi2);
  JetJetHistList->Add(h_JetJet_minv_FineBinning);


  for(int i=0;i<13;i++) {

    cout<<"Getting fake rate histograms"<<endl;
    cout<<ffHists[i].Data()<<endl;
    cout<<ftHists[i].Data()<<endl;
    cout<<tfHists[i].Data()<<endl;

    hFF[i]=(TH1F*)histoFF->Get(ffHists[i].Data());
    hFT[i]=(TH1F*)histoFT->Get(ftHists[i].Data());
    hTF[i]=(TH1F*)histoTF->Get(tfHists[i].Data());

    cout<<"Fake rate histograms: #entries "<<endl;
    cout<<"FF "<<hFF[i]->GetEntries()<<endl;;  
    cout<<"FT "<<hFT[i]->GetEntries()<<endl;;  
    cout<<"TF "<<hTF[i]->GetEntries()<<endl;;  

    histosGammaJet[i]= (TH1F*)GammaJetHistList->At(i);
    histosGammaJet[i]->Sumw2();
    histosGammaJet[i]->Add(hTF[i],hFT[i],1.,1.);

    //why is it subtracted twice ???
    //1 for TF and 1 for FT, I guess !!!
    histosGammaJet[i]->Add(hFF[i],-2.);

    cout<<"Gamma+Jet total contribution "<<histosGammaJet[i]->Integral()<<endl;

    histosJetJet[i]= (TH1F*)JetJetHistList->At(i);
    histosJetJet[i]->Sumw2();
    histosJetJet[i]->Add(hFF[i],1.); 

    cout<<"Jet+Jet total contribution "<<histosJetJet[i]->Integral()<<endl;

    //if(i == 0){
    if(i == 12){
      //Special treatment for upper bands
      h_GammaJet_minv_UpperBand->Add(hTF[i],hFT[i],1.05,1.05);
      h_GammaJet_minv_UpperBand->Add(hFF[i],-2.*1.05*1.05);
      h_GammaJet_minv_UpperBand->Add(histosGammaJet[i],-1.);
      h_GammaJet_minv_UpperBand->Scale(sqrt(3.));
      h_GammaJet_minv_UpperBand->Add(histosGammaJet[i],1.);

      h_JetJet_minv_UpperBand->Add(hFF[i],1.05*1.05); 
      h_JetJet_minv_UpperBand->Add(histosJetJet[i],-1.); 
      h_JetJet_minv_UpperBand->Scale(sqrt(3.));
      h_JetJet_minv_UpperBand->Add(histosJetJet[i],1.); 

      //And special treatment for lower bands
      h_GammaJet_minv_LowerBand->Add(hTF[i],hFT[i],0.95,0.95);
      h_GammaJet_minv_LowerBand->Add(hFF[i],-2.*0.95*0.95);
      h_GammaJet_minv_LowerBand->Add(histosGammaJet[i],-1.);
      h_GammaJet_minv_LowerBand->Scale(sqrt(3.));
      h_GammaJet_minv_LowerBand->Add(histosGammaJet[i],1.);

      h_JetJet_minv_LowerBand->Add(hFF[i],0.95*0.95); 
      h_JetJet_minv_LowerBand->Add(histosJetJet[i],-1.); 
      h_JetJet_minv_LowerBand->Scale(sqrt(3.));
      h_JetJet_minv_LowerBand->Add(histosJetJet[i],1.); 

      //Now we also compute the syst errors as 
      //a function of the mass and store it in a TH1
      h_GammaJet_minv_UpperError->Add(h_GammaJet_minv_UpperBand,histosGammaJet[i],1.,-1.);
      //h_GammaJet_minv_UpperError->Scale(sqrt(3.));
      h_GammaJet_minv_LowerError->Add(h_GammaJet_minv_LowerBand,histosGammaJet[i],-1.,1.);
      //h_GammaJet_minv_LowerError->Scale(sqrt(3.));

      h_JetJet_minv_UpperError->Add(h_JetJet_minv_UpperBand,histosJetJet[i],1.,-1.);
      //h_JetJet_minv_UpperError->Scale(sqrt(3.));
      h_JetJet_minv_LowerError->Add(h_JetJet_minv_LowerBand,histosJetJet[i],-1.,1.);
      //h_JetJet_minv_LowerError->Scale(sqrt(3.));
    }

    histoFileGammaJet->cd();
    histosGammaJet[i]->Write(); 

    histoFileJetJet->cd();
    histosJetJet[i]->Write();
      
  }//end of loop over histograms
  
  histoFileGammaJet->cd();
  h_GammaJet_minv_UpperBand->Write();
  h_GammaJet_minv_LowerBand->Write();
  h_GammaJet_minv_UpperError->Write();
  h_GammaJet_minv_LowerError->Write();
  histoFileGammaJet->Close();
  
  histoFileJetJet->cd();
  h_JetJet_minv_UpperBand->Write();
  h_JetJet_minv_LowerBand->Write();
  h_JetJet_minv_UpperError->Write();
  h_JetJet_minv_LowerError->Write();
  histoFileJetJet->Close();
  
  histoFF->cd();
  histoFF->Close();

  histoFT->cd();
  histoFT->Close();

  histoTF->cd();
  histoTF->Close();

  char canvasnamegj[10];
  char canvasnamejj[10];
  TCanvas*  cgj[13];
  TCanvas*  cjj[13];
      
  TPaveText *LumiLabel = new TPaveText(.6,.7,.8,.85,"NDC");
  LumiLabel->SetTextSize(0.03);
  LumiLabel->SetFillStyle(0);
  LumiLabel->SetBorderSize(0);
  LumiLabel->AddText("CMS Internal");
  TString pico(" pb^{-1}");
  TString Luminosity = lumi.Append(pico);
  LumiLabel->AddText(Luminosity.Data());
  LumiLabel->AddText(JSONFile.Data());

  for(int in=0;in<12;in++)
    {
      sprintf(canvasnamegj,"cgj%i",in);
      cgj[in] = new TCanvas(canvasnamegj, canvasnamegj, 800., 600.);

      histosGammaJet[in]->GetXaxis()->CenterTitle();
      histosGammaJet[in]->GetYaxis()->CenterTitle();
      histosGammaJet[in]->SetMarkerStyle(20);
      histosGammaJet[in]->SetMarkerColor(1);
      histosGammaJet[in]->Draw();

      sprintf(canvasnamejj,"cjj%i",in);
      cjj[in] = new TCanvas(canvasnamejj, canvasnamejj, 800., 600.);

      histosJetJet[in]->GetXaxis()->CenterTitle();
      histosJetJet[in]->GetYaxis()->CenterTitle();
      histosJetJet[in]->SetMarkerStyle(20);
      histosJetJet[in]->SetMarkerColor(1);
      histosJetJet[in]->Draw();

      if ((strcmp(histosGammaJet[in]->GetName(),"h_GammaJet_minv")==0 )|| (strcmp(histosGammaJet[in]->GetName(),"h_GammaJet_pt1")==0 ) || strcmp(histosGammaJet[in]->GetName(),"h_GammaJet_pt2")==0 ){
	cgj[in]->SetLogy(1);
	histosGammaJet[in]->Draw();
      }

      if ((strcmp(histosJetJet[in]->GetName(),"h_JetJet_minv")==0 )|| (strcmp(histosJetJet[in]->GetName(),"h_JetJet_pt1")==0 ) || strcmp(histosGammaJet[in]->GetName(),"h_JetJet_pt2")==0 ){
	cjj[in]->SetLogy(1);
	histosJetJet[in]->Draw();
      }

      LumiLabel->Draw();

      //cgj[in]->SaveAs(GammaJetHists[in]+"_"+Sample+".C");
      //cgj[in]->SaveAs(GammaJetHists[in]+"_"+Sample+".png");
      //cgj[in]->SaveAs(GammaJetHists[in]+"_"+Sample+".pdf");
      //cgj[in]->SaveAs(GammaJetHists[in]+"_"+Sample+".eps");
      //cjj[in]->SaveAs(JetJetHists[in]+"_"+Sample+".C");
      //cjj[in]->SaveAs(JetJetHists[in]+"_"+Sample+".png");
      //cjj[in]->SaveAs(JetJetHists[in]+"_"+Sample+".pdf");
      //cjj[in]->SaveAs(JetJetHists[in]+"_"+Sample+".eps");

    }//end of loop over histograms

}//end of method fakeratehistos



void OverlayMCFake(TString Sample = "Diphoton", TString lumi = "1", TString JSONFile = "N0JSON")
{

  gStyle->SetOptStat("ourmei");

  TCanvas* CumulativeCanvas = new TCanvas("CumulativeCanvas", "CumulativeCanvas", 800., 600.);
  TCanvas* dataoverlayedCumulative = new TCanvas("dataoverlayedCumulative", "Data Overlayed onto Cumulative Plot",800.,600.);
  THStack *StackCumm;

  TPaveText *LumiLabelCumm = new TPaveText(.55,.75,.75,.9,"NDC");
  LumiLabelCumm->SetTextSize(0.03);
  LumiLabelCumm->SetFillStyle(0);
  LumiLabelCumm->SetBorderSize(0);
  LumiLabelCumm->AddText("CMS Internal");
  TString LuminosityCumm = lumi;
  LumiLabelCumm->AddText(LuminosityCumm.Data());
  LumiLabelCumm->AddText(JSONFile.Data());

  TH1F* h_DatadivBack = new TH1F("h_DatadivBack","",199,20.,4000.); 
  TH1F* h_DataMinusBack = new TH1F("h_DataMinusBack","",199,20.,4000.); 
  TH1F* h_DataMinusBackdivBack = new TH1F("h_DataMinusBackdivBack","",199,20.,4000.); 
  TH1F* h_TotalBackground = new TH1F("h_TotalBackground","",199,20.,4000.);

  h_DatadivBack->Sumw2();
  h_DataMinusBack->Sumw2();
  h_DataMinusBackdivBack->Sumw2();
  h_TotalBackground->Sumw2();

  TH1F* h_DatadivBackZoom = new TH1F("h_DatadivBackZoom","",39,20.,800.); 
  TH1F* h_DataMinusBackdivBackZoom = new TH1F("h_DataMinusBackdivBackZoom","",39,20.,800.); 

  h_DatadivBackZoom->Sumw2();
  h_DataMinusBackdivBackZoom->Sumw2();

  TString histogramdata = HistogramFileLocation+Sample+"/histograms_"+Sample+".root";
  TFile* fdatahists = TFile::Open(histogramdata.Data());

  TString histogramJetJet =HistogramFileLocation+Sample+"/histograms_"+Sample+"_JetJet.root";
  TString histogramGammaJet =HistogramFileLocation+Sample+"/histograms_"+Sample+"_GammaJet.root";	
  TFile* fJetJethists=TFile::Open(histogramJetJet.Data());
  TFile* fGammaJethists=TFile::Open(histogramGammaJet.Data());

  //HARDCODED
  TString histogramMC =HistogramFileLocation+"/diphoton_tree_MC_all/histograms_diphoton_tree_MC_all.root";
  TFile* fMChists=TFile::Open(histogramMC.Data());

  TH1F* histosmc[12];
  TH1F* histosdata[12];
  TH1F* histosJetJet[12];
  TH1F* histosGammaJet[12];
    
  char canvasname[10];
  TCanvas*  c[12];
  
  TPaveText *LumiLabel = new TPaveText(.62,.75,.72,.90,"NDC");
  LumiLabel->SetTextSize(0.03);
  LumiLabel->SetFillStyle(0);
  LumiLabel->SetBorderSize(0);
  LumiLabel->SetTextSize(0.03);
  LumiLabel->SetFillStyle(0);
  LumiLabel->SetBorderSize(0);
  LumiLabel->AddText("CMS Internal");
  TString pico(" pb^{-1}");
  TString Luminosity = lumi.Append(pico);
  LumiLabel->AddText(Luminosity.Data());
  LumiLabel->AddText(JSONFile.Data());
  float lumiNumber = lumi.Atof();
  cout<<"lumi float value (pb^{-1})"<<lumiNumber<<endl;
  
  TLegend *StackLegend[12];
  THStack *StackMC[12];
  for(int i=0;i<12;i++){
   
    sprintf(canvasname,"c%i",i);
    c[i] = new TCanvas(canvasname, canvasname, 800., 600.);
    c[i]->cd();

    StackLegend[i]= new TLegend(0.23,0.71,0.41,0.86,"","NDC");
    StackLegend[i]->SetFillStyle(0);
    StackLegend[i]->SetTextSize(0.03);

    TString StackName(fNames[i]+"_Stack");
    //StackMC[i]= new THStack(StackName.Data(),StackName.Data());
    StackMC[i]= new THStack(StackName.Data(),"");

    cout<<"Getting data histogram "<<endl;
    histosdata[i]=(TH1F*)fdatahists->Get(fNames[i].Data());
    histosdata[i]->SetMarkerColor(1);
    cout<<"histo "<<histosdata[i]->GetName()<<" "<<histosdata[i]->GetEntries()<<" entries"<<histosdata[i]->Integral()<<" (integral)"<<endl;

    cout<<"Getting MC histogram "<<endl;
    histosmc[i]= (TH1F*)fMChists->Get(fNames[i].Data());
    histosmc[i]->SetFillColor(33);
    histosmc[i]->Scale(lumiNumber);
    cout<<"histo "<<histosmc[i]->GetName()<<" "<<histosmc[i]->GetEntries()<<" entries "<<histosmc[i]->Integral()<<" (integral)"<<endl;

    cout<<"Getting Jet+Jet histogram "<<endl;
    histosJetJet[i]=(TH1F*)fJetJethists->Get(JetJetHists[i].Data());
    histosJetJet[i]->SetFillColor(36); 
    cout<<"histo "<<histosJetJet[i]->GetName()<<" "<<histosJetJet[i]->GetEntries()<<" entries "<<histosJetJet[i]->Integral()<<" (integral)"<<endl;

    cout<<"Getting Gamma+Jet histogram "<<endl;
    histosGammaJet[i]=(TH1F*)fGammaJethists->Get(GammaJetHists[i].Data());
    histosGammaJet[i]->SetFillColor(38);
    cout<<"histo "<<histosGammaJet[i]->GetName()<<" "<<histosGammaJet[i]->GetEntries()<<" entries "<<histosGammaJet[i]->Integral()<<" (integral)"<<endl;


    if (i==0){
      h_TotalBackground->Add(histosJetJet[i],histosGammaJet[i],1.,1.);
      h_TotalBackground->Add(histosmc[i],1.);

      //We should have Data/Bckgd and (Data-Bckgd)/Bckgd
      h_DatadivBack->Divide(histosdata[i],h_TotalBackground,1.,1.);
      h_DataMinusBack->Add(histosdata[i],h_TotalBackground,1.,-1.);
      h_DataMinusBackdivBack->Divide(h_DataMinusBack,h_TotalBackground,1.,1.);

      //Zooming parts
      h_DatadivBackZoom->Divide(histosdata[i],h_TotalBackground,1.,1.);
      h_DataMinusBackdivBackZoom->Divide(h_DataMinusBack,h_TotalBackground,1.,1.);

      TCanvas* CanvasDatadivBack = new TCanvas("CDatadivBack","CDatadivBack", 800., 600.);
      CanvasDatadivBack->cd();

      //h_DatadivBack->SetMaximum(2.*h_DatadivBack->GetMaximum());	
      //h_DatadivBack->SetMinimum(0.);	
      h_DatadivBack->SetMaximum(8.2);	
      h_DatadivBack->SetMinimum(0.001);	
      h_DatadivBack->SetEntries(histosdata[i]->GetEntries());	
      h_DatadivBack->Draw();
      
      h_DatadivBack->GetXaxis()->SetTitle("M_{#gamma #gamma} (GeV)");
      h_DatadivBack->GetYaxis()->SetTitle("Data/Background");

      CanvasDatadivBack->SaveAs("h_DatadivBack_"+Sample+"CompleteOverlay"+".png");
      CanvasDatadivBack->SaveAs("h_DatadivBack_"+Sample+"CompleteOverlay"+".eps");
      CanvasDatadivBack->SaveAs("h_DatadivBack_"+Sample+"CompleteOverlay"+".C");
      CanvasDatadivBack->SaveAs("h_DatadivBack_"+Sample+"CompleteOverlay"+".pdf");

      Float_t maxvalue = max(h_DataMinusBackdivBack->GetMaximum(),fabs(h_DataMinusBackdivBack->GetMinimum()));

      //h_DataMinusBackdivBack->SetMaximum(2.*maxvalue);	
      //h_DataMinusBackdivBack->SetMinimum(-2.*maxvalue);	
      h_DataMinusBackdivBack->SetMaximum(6.);	
      h_DataMinusBackdivBack->SetMinimum(-6.);	
      h_DataMinusBackdivBack->SetEntries(histosdata[i]->GetEntries());	
      h_DataMinusBackdivBack->Draw();
	      
      h_DataMinusBackdivBack->GetXaxis()->SetTitle("M_{#gamma #gamma} (GeV)");
      h_DataMinusBackdivBack->GetYaxis()->SetTitle("(Data - Background)/Background");

      CanvasDatadivBack->SaveAs("h_DataMinusBackdivBack_"+Sample+"CompleteOverlay"+".png");
      CanvasDatadivBack->SaveAs("h_DataMinusBackdivBack_"+Sample+"CompleteOverlay"+".eps");
      CanvasDatadivBack->SaveAs("h_DataMinusBackdivBack_"+Sample+"CompleteOverlay"+".C");
      CanvasDatadivBack->SaveAs("h_DataMinusBackdivBack_"+Sample+"CompleteOverlay"+".pdf");

      //Zooming parts
      h_DatadivBack->SetMaximum(2.);	
      h_DatadivBack->SetMinimum(0.001);	
      h_DatadivBack->GetXaxis()->SetRangeUser(0.,800.);
      h_DatadivBack->SetEntries(histosdata[i]->GetEntries());	
      h_DatadivBack->Draw();
      

      CanvasDatadivBack->SaveAs("h_DatadivBack_Zoom_"+Sample+"CompleteOverlay"+".png");
      CanvasDatadivBack->SaveAs("h_DatadivBack_Zoom_"+Sample+"CompleteOverlay"+".eps");
      CanvasDatadivBack->SaveAs("h_DatadivBack_Zoom_"+Sample+"CompleteOverlay"+".C");
      CanvasDatadivBack->SaveAs("h_DatadivBack_Zoom_"+Sample+"CompleteOverlay"+".pdf");

      h_DataMinusBackdivBack->SetMaximum(1.5);	
      h_DataMinusBackdivBack->SetMinimum(-1.5);	
      h_DataMinusBackdivBack->GetXaxis()->SetRangeUser(0.,800.);
      h_DataMinusBackdivBack->SetEntries(histosdata[i]->GetEntries());	
      h_DataMinusBackdivBack->Draw();
	      
      CanvasDatadivBack->SaveAs("h_DataMinusBackdivBack_Zoom_"+Sample+"CompleteOverlay"+".png");
      CanvasDatadivBack->SaveAs("h_DataMinusBackdivBack_Zoom_"+Sample+"CompleteOverlay"+".eps");
      CanvasDatadivBack->SaveAs("h_DataMinusBackdivBack_Zoom_"+Sample+"CompleteOverlay"+".C");
      CanvasDatadivBack->SaveAs("h_DataMinusBackdivBack_Zoom_"+Sample+"CompleteOverlay"+".pdf");

    }

    c[i]->cd();    

    StackMC[i]->Add(histosJetJet[i]);
    StackMC[i]->Add(histosGammaJet[i]);
    StackMC[i]->Add(histosmc[i]);

    StackLegend[i]->AddEntry(histosdata[0],"Data", "LEP");
    StackLegend[i]->AddEntry(histosmc[0],"SM Diphoton","f");
    StackLegend[i]->AddEntry(histosGammaJet[0],"Photon+Jet","f");
    StackLegend[i]->AddEntry(histosJetJet[0],"Jet+Jet","f");

    if ( (strcmp(fNames[i],"h_Diphoton_Minv_high")==0 ) || (strcmp(fNames[i],"h_Photon1_pt")==0 ) || (strcmp(fNames[i],"h_Diphoton_Minv")==0) || (strcmp(fNames[i],"h_Photon2_pt")==0) || (strcmp(fNames[i],"h_Diphoton_qt")==0) || (strcmp(fNames[i],"h_Diphoton_deltaR")==0) ){
      StackMC[i]->SetMaximum(10.*histosdata[i]->GetMaximum());
      histosdata[i]->SetMaximum(10.*histosdata[i]->GetMaximum());
      c[i]->SetLogy(1);
    }
    
    histosdata[i]->SetMinimum(1.e-2);
    StackMC[i]->SetMinimum(1.e-2);

    StackMC[i]->SetMaximum(1.7*histosdata[i]->GetMaximum());
    histosdata[i]->SetMaximum(1.7*histosdata[i]->GetMaximum());

    StackMC[i]->Draw("HIST");
    histosdata[i]->Draw("EP,SAMES");

    histosdata[i]->GetXaxis()->SetTitle(XTitles[i].Data());
    histosdata[i]->GetYaxis()->SetTitle(YTitles[i].Data());
    StackMC[i]->GetHistogram()->GetXaxis()->SetTitle(XTitles[i].Data());
    StackMC[i]->GetHistogram()->GetYaxis()->SetTitle(YTitles[i].Data());

    if ((strcmp(fNames[i],"h_Diphoton_cosThetaStar")==0) || (strcmp(fNames[i],"h_Photon1_phi")==0) ||  (strcmp(fNames[i],"h_Photon2_phi")==0)){
      histosdata[i]->SetMinimum(0.);
    }

    LumiLabel->Draw("sames");
    StackLegend[i]->Draw("sames");
    gPad->RedrawAxis();
    gPad->Update();
    
    c[i]->SaveAs(fNames[i]+"_"+Sample.Data()+"CompleteOverlay"+".C");
    c[i]->SaveAs(fNames[i]+"_"+Sample.Data()+"CompleteOverlay"+".png");
    c[i]->SaveAs(fNames[i]+"_"+Sample.Data()+"CompleteOverlay"+".pdf");
    c[i]->SaveAs(fNames[i]+"_"+Sample.Data()+"CompleteOverlay"+".eps");

    if (i==2){
      c[i]->SetLogy(1);

      StackMC[i]->SetMaximum(27.*histosdata[i]->GetMaximum());
      histosdata[i]->SetMaximum(27.*histosdata[i]->GetMaximum());

      c[i]->SaveAs(fNames[i]+"_log_"+Sample.Data()+"CompleteOverlay"+".png"); 
      c[i]->SaveAs(fNames[i]+"_log_"+Sample.Data()+"CompleteOverlay"+".eps"); 
      c[i]->SaveAs(fNames[i]+"_log_"+Sample.Data()+"CompleteOverlay"+".pdf"); 
      c[i]->SaveAs(fNames[i]+"_log_"+Sample.Data()+"CompleteOverlay"+".C"); 
    }      

    if (i==0){

      CumulativeCanvas->cd();
      CumulativeCanvas->SetLogy(1);

      TString StackCummName(fNames[i]+"_Stack");
      StackCumm = new THStack(StackCummName.Data(),"");

      TH1F* h_Diphoton_Minv_log_DataCumm =  MakeChists(histosdata[0]);
      TH1F* h_Diphoton_Minv_log_MCCumm = MakeChists(histosmc[0]); 
      TH1F* h_Diphoton_Minv_log_JetJetCumm = MakeChists(histosJetJet[0]);
      TH1F* h_Diphoton_Minv_log_GammaJetCumm = MakeChists(histosGammaJet[0]);
      
      h_Diphoton_Minv_log_MCCumm->SetFillColor(33);
      h_Diphoton_Minv_log_JetJetCumm->SetFillColor(36);
      h_Diphoton_Minv_log_GammaJetCumm->SetFillColor(38);
      
      StackCumm->Add(h_Diphoton_Minv_log_JetJetCumm);
      StackCumm->Add(h_Diphoton_Minv_log_GammaJetCumm);
      StackCumm->Add(h_Diphoton_Minv_log_MCCumm);
      
      TLegend* DataMCLegendCumm = new TLegend(0.23,0.71,0.41,0.86,"","NDC");
      DataMCLegendCumm->AddEntry( h_Diphoton_Minv_log_DataCumm,"Data","LEP");
      DataMCLegendCumm->AddEntry(h_Diphoton_Minv_log_MCCumm,"SM Diphoton","F");
      DataMCLegendCumm->AddEntry(h_Diphoton_Minv_log_GammaJetCumm,"Gamma+Jet","F");
      DataMCLegendCumm->AddEntry(h_Diphoton_Minv_log_JetJetCumm,"Jet+Jet","F");
      DataMCLegendCumm->SetFillColor(0);

      StackCumm->SetMaximum(25*h_Diphoton_Minv_log_DataCumm->GetMaximum());
      h_Diphoton_Minv_log_DataCumm->SetMaximum(25*h_Diphoton_Minv_log_DataCumm->GetMaximum());

      StackCumm->Draw("HIST");
      h_Diphoton_Minv_log_DataCumm->Draw("EP,SAMES");
      gPad->RedrawAxis();

      h_Diphoton_Minv_log_DataCumm->GetXaxis()->SetTitle(XTitles[i].Data());
      h_Diphoton_Minv_log_DataCumm->GetYaxis()->SetTitle(YTitles[i].Data());
      StackCumm->GetHistogram()->GetXaxis()->SetTitle(XTitles[i].Data());
      StackCumm->GetHistogram()->GetYaxis()->SetTitle(YTitles[i].Data());

      LumiLabelCumm->Draw();
      DataMCLegendCumm->Draw();
      CumulativeCanvas->SaveAs(fNames[i]+"_"+Sample.Data()+"Cumulative"+".png");
      CumulativeCanvas->SaveAs(fNames[i]+"_"+Sample.Data()+"Cumulative"+".eps");
      CumulativeCanvas->SaveAs(fNames[i]+"_"+Sample.Data()+"Cumulative"+".pdf");
      CumulativeCanvas->SaveAs(fNames[i]+"_"+Sample.Data()+"Cumulative"+".C");
      
      dataoverlayedCumulative->cd();
      dataoverlayedCumulative->SetLogy(1);
      StackCumm->Draw("HIST");
      h_Diphoton_Minv_log_DataCumm->Draw("EP,SAMES");
      gPad->RedrawAxis();
      LumiLabelCumm->Draw();
      DataMCLegendCumm->Draw();
      histosdata[i]->SetMarkerColor(kRed);   
      
      histosdata[i]->Draw("EP,SAMES");
      dataoverlayedCumulative->SaveAs(fNames[i]+Sample.Data()+"DataOverlayedCumulative"+".png");
      dataoverlayedCumulative->SaveAs(fNames[i]+Sample.Data()+"DataOverlayedCumulative"+".eps");
      dataoverlayedCumulative->SaveAs(fNames[i]+Sample.Data()+"DataOverlayedCumulative"+".pdf");
      dataoverlayedCumulative->SaveAs(fNames[i]+Sample.Data()+"DataOverlayedCumulative"+".C");
             
    }
    gPad->Update();

  }//end of loop over histos

}//end of method OverlayMCFake
















































void CreateValidationPlotFiles(TString Sample = "Validation", TString SampleType = "mc", TString JSON = "NOJSON")
{
 
  cout << "Entering CreateHistogramFiles method with parameters:" <<endl;
  cout<<"Sample: "<<Sample.Data()<<" Sample type: "<<SampleType.Data()<<" JSON: "<<JSON.Data()<<endl;

  TString outName;
  TFile *outfilename;

  TString  inputfile= TreeFileLocation+Sample+".root";
   
  TChain *chain_tt = new TChain("mcprodValidator/fTree");
  chain_tt->Add(inputfile.Data());
  cout << "TT entries = " << chain_tt->GetEntries() <<endl;
  
  cout << "CreateHistogramFiles: Start processing main loop entries" <<endl;

  PlottingCodeLoop* treereader = new PlottingCodeLoop(chain_tt);
  treereader->_cutPhoton1pt = 70.;
  treereader->_cutPhoton2pt = 70.;
  treereader->_fakeStatus="TightTight";
  treereader->_JSON=JSON.Data();
  treereader->_SampleType=SampleType.Data();
  outName = TString::Format("histograms_%s.root", Sample.Data());
  outfilename = new TFile(HistogramFileLocation+Sample.Data()+"/"+outName.Data(),"RECREATE");
  treereader->_outputfile = outfilename;
  treereader->Loop();
  
  cout << "CreateHistogramFiles: Finished processing main loop entries" <<endl;

  outfilename->cd();
  outfilename->Close();

}
























void MakeValidationPlots( TString Sample  = "Validation", TString lumi = "10000", TString JSONFile = "NOJSON")
{

  cout<<"Entering MakeValidationPlots method with parameters: "<<endl; 
  cout<<"Sample: "<<Sample.Data()<<" JSON: "<<JSONFile.Data()<<" lumi: "<<lumi.Data()<<" /pb"<<endl;

  TString histogramCentralProduction =HistogramFileLocation+"diphoton_tree_PrivateProduction_RSGravToGG_kMpl01_M3250/histograms_diphoton_PrivateProduction_RSGravToGG_kMpl01_M3250.root";
  TString histogramPrivateProduction=HistogramFileLocation+"diphoton_tree_PrivateProduction_RSGravToGG_kMpl01_M3250/histograms_diphoton_PrivateProduction_RSGravToGG_kMpl01_M3250.root";
  TFile* fPrivateProductionhists = TFile::Open(histogramPrivateProduction.Data());
  TFile* fCentralProductionhists = TFile::Open(histogramCentralProduction.Data());

  Int_t nGenCentralProduction = 25080;
  Int_t nGenPrivateProduction = 25080;
  char canvasname[10];
  TCanvas* c[42];

  TString pico(" pb^{-1}");
  TPaveText *LumiLabel = new TPaveText(.55,.75,.75,.9,"NDC");
  LumiLabel->SetTextSize(0.03);
  LumiLabel->SetFillStyle(0);
  LumiLabel->SetBorderSize(0);
  LumiLabel->AddText("CMS Internal");

  TString Luminosity = lumi.Append(pico);
  LumiLabel->AddText(Luminosity.Data());
  LumiLabel->AddText(JSONFile.Data());

  TH1F* histosCentralProduction[42];
  TH1F* histosPrivateProduction[42];
  TLegend* Legend[42];
  float lumiNumber = lumi.Atof();
  cout<<"Lumi number "<<lumiNumber<<endl;

  for(int i=0;i<42;i++)
    {

      sprintf(canvasname,"c%i",i);
      c[i] = new TCanvas(canvasname, canvasname, 800., 600.);
      cout<<"Canvas "<<canvasname<<endl;
      c[i]->cd();

      histosCentralProduction[i] = (TH1F*)fPrivateProductionhists->Get(nameHists[i].Data());
      histosCentralProduction[i]->Scale((lumiNumber*1.91e-5)/nGenCentralProduction);
      histosCentralProduction[i]->SetMinimum(.01);
      cout<<"CentralProduction histogram name"<<histosCentralProduction[i]->GetName()<<endl;

      histosPrivateProduction[i] = (TH1F*)fPrivateProductionhists->Get(nameHists[i].Data());
      histosPrivateProduction[i]->Scale((lumiNumber*1.91e-5)/nGenPrivateProduction);
      histosPrivateProduction[i]->SetMinimum(.01);
      cout<<"PrivateProduction histogram name"<<histosPrivateProduction[i]->GetName()<<endl;

      Legend[i]= new TLegend(0.70,0.6,0.88,0.75,"","NDC");
      Legend[i]->SetFillStyle(0);
      Legend[i]->SetTextSize(0.03);

      histosPrivateProduction[i]->GetXaxis()->CenterTitle();
      histosPrivateProduction[i]->GetYaxis()->CenterTitle();
      histosPrivateProduction[i]->SetFillColor(kBlue);
      histosPrivateProduction[i]->SetMaximum(1.3*histosCentralProduction[i]->GetMaximum());

      histosCentralProduction[i]->SetMarkerColor(1);
      histosCentralProduction[i]->SetMaximum(1.3*histosCentralProduction[i]->GetMaximum());
      //histosdata[i]->SetMinimum(.8*histosmc[i]->GetMinimum());
      //histosdata[i]->Draw("HIST EP");

      //histosmc[i]->Draw("SAME HIST");
      histosPrivateProduction[i]->Draw("HIST");
      //histosdata[i]->Draw("HIST SAME EP");
      histosCentralProduction[i]->Draw("EP,SAME");

      Legend[i]->AddEntry(histosCentralProduction[i],"CentralProduction","LEP");
      Legend[i]->AddEntry(histosPrivateProduction[i],"SM Diphoton","F");
      gPad->RedrawAxis();
      //histosdata[i]->Draw("EPSAME");
      LumiLabel->Draw();
      Legend[i]->Draw();

      if (nameHists[i].Contains("_log"))
        {
	  c[i]->SetLogy(1);
	}

      c[i]->SaveAs(nameHists[i]+"_"+Sample.Data()+"Overlay"+".C");
      c[i]->SaveAs(nameHists[i]+"_"+Sample.Data()+"Overlay"+".png");
      c[i]->SaveAs(nameHists[i]+"_"+Sample.Data()+"Overlay"+".pdf");
      c[i]->SaveAs(nameHists[i]+"_"+Sample.Data()+"Overlay"+".eps");
    }//end of loop over histograms


  TCanvas* CumulativeCanvas = new TCanvas("CumulativeCanvas", "CumulativeCanvas", 800., 600.);
  CumulativeCanvas->cd();
  CumulativeCanvas->SetLogy(1);
  TH1F* h_Diphoton_Minv_log_CentralProductionCumm =  MakeChists(histosCentralProduction[9]);
  TH1F* h_Diphoton_Minv_log_PRIVATEPRODUCTIONCumm = MakeChists(histosPrivateProduction[9]);

  h_Diphoton_Minv_log_PRIVATEPRODUCTIONCumm->SetFillColor(kBlue);

  TPaveText *LumiLabelCumm = new TPaveText(.55,.75,.75,.9,"NDC");
  LumiLabelCumm->SetTextSize(0.03);
  LumiLabelCumm->SetFillStyle(0);
  LumiLabelCumm->SetBorderSize(0);
  LumiLabelCumm->AddText("CMS Internal");
  LumiLabelCumm->AddText(TString(lumi+" pb^{-1}").Data());
  LumiLabelCumm->AddText(JSONFile.Data());

  TLegend* LegendCumm = new TLegend(0.70,0.6,0.88,0.75,"","NDC");
  LegendCumm->AddEntry( h_Diphoton_Minv_log_CentralProductionCumm,"CentralProduction","LEP");
  LegendCumm->AddEntry(h_Diphoton_Minv_log_PRIVATEPRODUCTIONCumm,"SM Diphoton","F");
  LegendCumm->SetFillColor(0);


  h_Diphoton_Minv_log_CentralProductionCumm->SetMaximum(2.5*h_Diphoton_Minv_log_CentralProductionCumm->GetMaximum());
  h_Diphoton_Minv_log_PRIVATEPRODUCTIONCumm->SetMaximum(2.5*h_Diphoton_Minv_log_CentralProductionCumm->GetMaximum());

  h_Diphoton_Minv_log_CentralProductionCumm->SetMinimum(.01);
  h_Diphoton_Minv_log_PRIVATEPRODUCTIONCumm->SetMinimum(.01);

  h_Diphoton_Minv_log_PRIVATEPRODUCTIONCumm->Draw("HIST");
  h_Diphoton_Minv_log_CentralProductionCumm->Draw("EPSAME");

  gPad->RedrawAxis();
  LumiLabelCumm->Draw();
  LegendCumm->Draw();
  CumulativeCanvas->SaveAs(nameHists[43]+Sample.Data()+"Overlay"+".png");
  CumulativeCanvas->SaveAs(nameHists[43]+Sample.Data()+"Overlay"+".eps");
  CumulativeCanvas->SaveAs(nameHists[43]+Sample.Data()+"Overlay"+".pdf");
  CumulativeCanvas->SaveAs(nameHists[43]+Sample.Data()+"Overlay"+".C");

}//end of method MakeValidationPlots


























































void MakeYieldsTable(TString Sample = "Diphoton", Float_t lumiNumber = 1.)
{

  gStyle->SetOptStat("ourmei");

  TString histogramdata = HistogramFileLocation+Sample+"/histograms_"+Sample+".root";
  TFile* fdatahists = TFile::Open(histogramdata.Data());

  TString histogramJetJet =HistogramFileLocation+Sample+"/histograms_"+Sample+"_JetJet.root";
  TString histogramGammaJet =HistogramFileLocation+Sample+"/histograms_"+Sample+"_GammaJet.root";	
  TFile* fJetJethists=TFile::Open(histogramJetJet.Data());
  TFile* fGammaJethists=TFile::Open(histogramGammaJet.Data());

  //HARDCODED
  TString histogramMC =HistogramFileLocation+"/diphoton_tree_MC_all/histograms_diphoton_tree_MC_all.root";
  TFile* fMChists=TFile::Open(histogramMC.Data());

  TH1F* histosmc;
  TH1F* histosdata;
  TH1F* histosJetJet;
  TH1F* histosJetJetUpperError;
  TH1F* histosJetJetLowerError;
  TH1F* histosGammaJet;
  TH1F* histosGammaJetUpperError;
  TH1F* histosGammaJetLowerError;
  

  cout<<"Getting data histogram "<<endl;
  histosdata=(TH1F*)fdatahists->Get("h_Diphoton_Minv_FineBinning");
  cout<<"histo "<<histosdata->GetName()<<" "<<histosdata->GetEntries()<<" entries"<<histosdata->Integral()<<" (integral)"<<endl;

  cout<<"Getting MC histogram "<<endl;
  histosmc= (TH1F*)fMChists->Get("h_Diphoton_Minv_FineBinning");
  histosmc->Scale(lumiNumber);
  cout<<"histo "<<histosmc->GetName()<<" "<<histosmc->GetEntries()<<" entries "<<histosmc->Integral()<<" (integral)"<<endl;

  cout<<"Getting Jet+Jet histograms "<<endl;
  histosJetJet=(TH1F*)fJetJethists->Get("h_JetJet_minv_FineBinning");
  histosJetJetUpperError=(TH1F*)fJetJethists->Get("h_JetJet_minv_UpperError");
  histosJetJetLowerError=(TH1F*)fJetJethists->Get("h_JetJet_minv_LowerError");
  cout<<"histo "<<histosJetJet->GetName()<<" "<<histosJetJet->GetEntries()<<" entries "<<histosJetJet->Integral()<<" (integral)"<<endl;
  cout<<"histo "<<histosJetJetUpperError->GetName()<<" "<<histosJetJetUpperError->GetEntries()<<" entries "<<histosJetJetUpperError->Integral()<<" (integral)"<<endl;
  cout<<"histo "<<histosJetJetLowerError->GetName()<<" "<<histosJetJetLowerError->GetEntries()<<" entries "<<histosJetJetLowerError->Integral()<<" (integral)"<<endl;

  cout<<"Getting Gamma+Jet histograms "<<endl;
  histosGammaJet=(TH1F*)fGammaJethists->Get("h_GammaJet_minv_FineBinning");
  histosGammaJetUpperError=(TH1F*)fGammaJethists->Get("h_GammaJet_minv_UpperError");
  histosGammaJetLowerError=(TH1F*)fGammaJethists->Get("h_GammaJet_minv_LowerError");
  cout<<"histo "<<histosGammaJet->GetName()<<" "<<histosGammaJet->GetEntries()<<" entries "<<histosGammaJet->Integral()<<" (integral)"<<endl;
  cout<<"histo "<<histosGammaJetUpperError->GetName()<<" "<<histosGammaJetUpperError->GetEntries()<<" entries "<<histosGammaJetUpperError->Integral()<<" (integral)"<<endl;
  cout<<"histo "<<histosGammaJetLowerError->GetName()<<" "<<histosGammaJetLowerError->GetEntries()<<" entries "<<histosGammaJetLowerError->Integral()<<" (integral)"<<endl;

  //special part to calculate integrals in different mass ranges

  //They all have the same binning
  //Let's take the data one
  Int_t binnr200 = -1;
  Int_t binnr500 = -1;
  Int_t binnr750 = -1;
  Int_t binnr1000 = -1;
  Int_t binnr1250 = -1;

  Int_t nbinsX = histosdata->GetNbinsX();

  for(int nbin = 0; nbin <  histosdata->GetNbinsX(); nbin++)
    {
      Float_t binlowedge = histosdata->GetBinLowEdge(nbin);
      //cout<<"binlowedge "<<binlowedge<<endl;
      if(binlowedge >= 200. && binnr200 == -1 ) {binnr200 = nbin;}
      if(binlowedge >= 500. && binnr500 == -1 ) {binnr500 = nbin;}
      if(binlowedge >= 750. && binnr750 == -1 ) {binnr750 = nbin;}
      if(binlowedge >= 1000. && binnr1000 == -1 ) {binnr1000 = nbin;}
      if(binlowedge >= 1250. && binnr1250 == -1 ) {binnr1250 = nbin;}
    }
  cout<<" binnr200 "<<binnr200<<endl;
  cout<<" binnr500 "<<binnr500<<endl;
  cout<<" binnr750 "<<binnr750<<endl;
  cout<<" binnr1000 "<<binnr1000<<endl;
  cout<<" binnr1250 "<<binnr1250<<endl;

  //Once we have the bin numbers
  //we compute the integrals
  //and their stat uncertainties
  //and their syst uncertainties

  //For each contribution, we take nbinsX+1 to have the overflow

  Float_t entriesJetJet = histosJetJet->Integral()+histosJetJet->GetBinContent(nbinsX+1);
  Float_t entries200JetJet = histosJetJet->Integral(binnr200,nbinsX+1);
  Float_t entries500JetJet = histosJetJet->Integral(binnr500,nbinsX+1);
  Float_t entries750JetJet = histosJetJet->Integral(binnr750,nbinsX+1);
  Float_t entries1000JetJet = histosJetJet->Integral(binnr1000,nbinsX+1);
  Float_t entries1250JetJet = histosJetJet->Integral(binnr1250,nbinsX+1);
	
  Float_t errorsJetJet = histosJetJetUpperError->Integral()+histosJetJetUpperError->GetBinContent(nbinsX+1);
  Float_t errors200JetJet = histosJetJetUpperError->Integral(binnr200,nbinsX+1);
  Float_t errors500JetJet = histosJetJetUpperError->Integral(binnr500,nbinsX+1);
  Float_t errors750JetJet = histosJetJetUpperError->Integral(binnr750,nbinsX+1);
  Float_t errors1000JetJet = histosJetJetUpperError->Integral(binnr1000,nbinsX+1);
  Float_t errors1250JetJet = histosJetJetUpperError->Integral(binnr1250,nbinsX+1);
	




  printf("%s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s \n","Jet+Jet & ",entriesJetJet," & "
	 ,entries200JetJet," & "
	 ,entries500JetJet," & "
	 ,entries750JetJet," & "
	 ,entries1000JetJet," & "
	 ,entries1250JetJet," \\ "
	 );

  printf("%s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s \n","stat & ",sqrt(entriesJetJet)," & "
	 ,sqrt(entries200JetJet)," & "
	 ,sqrt(entries500JetJet)," & "
	 ,sqrt(entries750JetJet)," & "
	 ,sqrt(entries1000JetJet)," & "
	 ,sqrt(entries1250JetJet)," \\ "
	 );

  printf("%s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s \n","syst & ",errorsJetJet," & "
	 ,errors200JetJet," & "
	 ,errors500JetJet," & "
	 ,errors750JetJet," & "
	 ,errors1000JetJet," & "
	 ,errors1250JetJet," \\ "
	 );


	
  Float_t entriesGammaJet = histosGammaJet->Integral()+histosGammaJet->GetBinContent(nbinsX+1);
  Float_t entries200GammaJet = histosGammaJet->Integral(binnr200,nbinsX+1);
  Float_t entries500GammaJet = histosGammaJet->Integral(binnr500,nbinsX+1);
  Float_t entries750GammaJet = histosGammaJet->Integral(binnr750,nbinsX+1);
  Float_t entries1000GammaJet = histosGammaJet->Integral(binnr1000,nbinsX+1);
  Float_t entries1250GammaJet = histosGammaJet->Integral(binnr1250,nbinsX+1);
	
  Float_t errorsGammaJet = histosGammaJetUpperError->Integral()+histosGammaJetUpperError->GetBinContent(nbinsX+1);
  Float_t errors200GammaJet = histosGammaJetUpperError->Integral(binnr200,nbinsX+1);
  Float_t errors500GammaJet = histosGammaJetUpperError->Integral(binnr500,nbinsX+1);
  Float_t errors750GammaJet = histosGammaJetUpperError->Integral(binnr750,nbinsX+1);
  Float_t errors1000GammaJet = histosGammaJetUpperError->Integral(binnr1000,nbinsX+1);
  Float_t errors1250GammaJet = histosGammaJetUpperError->Integral(binnr1250,nbinsX+1);
	
  printf("%s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s \n","Gamma+Jet & ",entriesGammaJet," & "
	 ,entries200GammaJet," & "
	 ,entries500GammaJet," & "
	 ,entries750GammaJet," & "
	 ,entries1000GammaJet," & "
	 ,entries1250GammaJet," \\ "
	 );

  printf("%s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s \n","stat & ",sqrt(entriesGammaJet)," & "
	 ,sqrt(entries200GammaJet)," & "
	 ,sqrt(entries500GammaJet)," & "
	 ,sqrt(entries750GammaJet)," & "
	 ,sqrt(entries1000GammaJet)," & "
	 ,sqrt(entries1250GammaJet)," \\ "
	 );

  printf("%s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s \n","syst & ",errorsGammaJet," & "
	 ,errors200GammaJet," & "
	 ,errors500GammaJet," & "
	 ,errors750GammaJet," & "
	 ,errors1000GammaJet," & "
	 ,errors1250GammaJet," \\ "
	 );


  Float_t entriesFake = entriesGammaJet + entriesJetJet;
  Float_t entries200Fake = entries200GammaJet + entries200JetJet;
  Float_t entries500Fake = entries500GammaJet + entries500JetJet;
  Float_t entries750Fake = entries750GammaJet + entries750JetJet;
  Float_t entries1000Fake = entries1000GammaJet + entries1000JetJet;
  Float_t entries1250Fake = entries1250GammaJet + entries1250JetJet;
	
  Float_t errorsFake = sqrt( (errorsGammaJet * errorsGammaJet) + (errorsJetJet * errorsJetJet) ); 
  Float_t errors200Fake = sqrt( (errors200GammaJet * errors200GammaJet) + (errors200JetJet * errors200JetJet) ); 
  Float_t errors500Fake = sqrt( (errors500GammaJet * errors500GammaJet) + (errors500JetJet * errors500JetJet) ); 
  Float_t errors750Fake = sqrt( (errors750GammaJet * errors750GammaJet) + (errors750JetJet * errors750JetJet) ); 
  Float_t errors1000Fake = sqrt( (errors1000GammaJet * errors1000GammaJet) + (errors1000JetJet * errors1000JetJet) ); 
  Float_t errors1250Fake = sqrt( (errors1250GammaJet * errors1250GammaJet) + (errors1250JetJet * errors1250JetJet) ); 
	
  printf("%s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s \n","Total Fake & ",entriesFake," & "
	 ,entries200Fake," & "
	 ,entries500Fake," & "
	 ,entries750Fake," & "
	 ,entries1000Fake," & "
	 ,entries1250Fake," \\ "
	 );

  printf("%s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s \n","stat & ",sqrt(entriesFake)," & "
	 ,sqrt(entries200Fake)," & "
	 ,sqrt(entries500Fake)," & "
	 ,sqrt(entries750Fake)," & "
	 ,sqrt(entries1000Fake)," & "
	 ,sqrt(entries1250Fake)," \\ "
	 );

  printf("%s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s \n","syst & ",errorsFake," & "
	 ,errors200Fake," & "
	 ,errors500Fake," & "
	 ,errors750Fake," & "
	 ,errors1000Fake," & "
	 ,errors1250Fake," \\ "
	 );


  Float_t entriesMC = histosmc->Integral()+histosmc->GetBinContent(nbinsX+1);
  Float_t entries200MC = histosmc->Integral(binnr200,nbinsX+1);
  Float_t entries500MC = histosmc->Integral(binnr500,nbinsX+1);
  Float_t entries750MC = histosmc->Integral(binnr750,nbinsX+1);
  Float_t entries1000MC = histosmc->Integral(binnr1000,nbinsX+1);
  Float_t entries1250MC = histosmc->Integral(binnr1250,nbinsX+1);
	
  //No syst uncertainties for MC yet
  Float_t errorsMC = 0.;
  Float_t errors200MC = 0.;
  Float_t errors500MC = 0.;
  Float_t errors750MC = 0.;
  Float_t errors1000MC = 0.;
  Float_t errors1250MC = 0.;
	
  printf("%s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s \n","MC & ",entriesMC," & "
	 ,entries200MC," & "
	 ,entries500MC," & "
	 ,entries750MC," & "
	 ,entries1000MC," & "
	 ,entries1250MC," \\ "
	 );

  printf("%s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s \n","stat & ",sqrt(entriesMC)," & "
	 ,sqrt(entries200MC)," & "
	 ,sqrt(entries500MC)," & "
	 ,sqrt(entries750MC)," & "
	 ,sqrt(entries1000MC)," & "
	 ,sqrt(entries1250MC)," \\ "
	 );

  printf("%s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s \n","syst & ",errorsMC," & "
	 ,errors200MC," & "
	 ,errors500MC," & "
	 ,errors750MC," & "
	 ,errors1000MC," & "
	 ,errors1250MC," \\ "
	 );



  Float_t entriesBackground = entriesGammaJet + entriesJetJet + entriesMC;
  Float_t entries200Background = entries200GammaJet + entries200JetJet + entries200MC;
  Float_t entries500Background = entries500GammaJet + entries500JetJet + entries500MC;
  Float_t entries750Background = entries750GammaJet + entries750JetJet + entries750MC;
  Float_t entries1000Background = entries1000GammaJet + entries1000JetJet + entries1000MC;
  Float_t entries1250Background = entries1250GammaJet + entries1250JetJet + entries1250MC;
	
  Float_t errorsBackground = sqrt( (errorsGammaJet * errorsGammaJet) + (errorsJetJet * errorsJetJet) + (errorsMC * errorsMC) ); 
  Float_t errors200Background = sqrt( (errors200GammaJet * errors200GammaJet) + (errors200JetJet * errors200JetJet) + (errors200MC * errors200MC) ); 
  Float_t errors500Background = sqrt( (errors500GammaJet * errors500GammaJet) + (errors500JetJet * errors500JetJet) + (errors500MC * errors500MC) ); 
  Float_t errors750Background = sqrt( (errors750GammaJet * errors750GammaJet) + (errors750JetJet * errors750JetJet) + (errors750MC * errors750MC) ); 
  Float_t errors1000Background = sqrt( (errors1000GammaJet * errors1000GammaJet) + (errors1000JetJet * errors1000JetJet) + (errors1000MC * errors1000MC) ); 
  Float_t errors1250Background = sqrt( (errors1250GammaJet * errors1250GammaJet) + (errors1250JetJet * errors1250JetJet) + (errors1250MC * errors1250MC) ); 
	
  printf("%s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s \n","Total Background & ",entriesBackground," & "
	 ,entries200Background," & "
	 ,entries500Background," & "
	 ,entries750Background," & "
	 ,entries1000Background," & "
	 ,entries1250Background," \\ "
	 );

  printf("%s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s \n","stat & ",sqrt(entriesBackground)," & "
	 ,sqrt(entries200Background)," & "
	 ,sqrt(entries500Background)," & "
	 ,sqrt(entries750Background)," & "
	 ,sqrt(entries1000Background)," & "
	 ,sqrt(entries1250Background)," \\ "
	 );

  printf("%s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s \n","syst & ",errorsBackground," & "
	 ,errors200Background," & "
	 ,errors500Background," & "
	 ,errors750Background," & "
	 ,errors1000Background," & "
	 ,errors1250Background," \\ "
	 );


  Float_t entriesData = histosdata->Integral()+histosdata->GetBinContent(nbinsX+1);
  Float_t entries200Data = histosdata->Integral(binnr200,nbinsX+1);
  Float_t entries500Data = histosdata->Integral(binnr500,nbinsX+1);
  Float_t entries750Data = histosdata->Integral(binnr750,nbinsX+1);
  Float_t entries1000Data = histosdata->Integral(binnr1000,nbinsX+1);
  Float_t entries1250Data = histosdata->Integral(binnr1250,nbinsX+1);
	
  Float_t errorsData = 0.;
  Float_t errors200Data = 0.;
  Float_t errors500Data = 0.;
  Float_t errors750Data = 0.;
  Float_t errors1000Data = 0.;
  Float_t errors1250Data = 0.;
	
  printf("%s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s \n","Data & ",entriesData," & "
	 ,entries200Data," & "
	 ,entries500Data," & "
	 ,entries750Data," & "
	 ,entries1000Data," & "
	 ,entries1250Data," \\ "
	 );
	
  printf("%s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s \n","stat & ",sqrt(entriesData)," & "
	 ,sqrt(entries200Data)," & "
	 ,sqrt(entries500Data)," & "
	 ,sqrt(entries750Data)," & "
	 ,sqrt(entries1000Data)," & "
	 ,sqrt(entries1250Data)," \\ "
	 );

	



  Float_t entriesDataOverBackground = entriesData/entriesBackground;
  Float_t entries200DataOverBackground = entries200Data/entries200Background;
  Float_t entries500DataOverBackground = entries500Data/entries500Background;
  Float_t entries750DataOverBackground = entries750Data/entries750Background;
  Float_t entries1000DataOverBackground = entries1000Data/entries1000Background;
  Float_t entries1250DataOverBackground = entries1250Data/entries1250Background;
	
  Float_t errorsDataOverBackground = entriesDataOverBackground * sqrt( (errorsData/entriesData)*(errorsData/entriesData) + (errorsBackground/entriesBackground)*(errorsBackground/entriesBackground) );
  Float_t errors200DataOverBackground = entries200DataOverBackground * sqrt( (errors200Data/entries200Data)*(errors200Data/entries200Data) + (errors200Background/entries200Background)*(errors200Background/entries200Background) );
  Float_t errors500DataOverBackground = entries500DataOverBackground * sqrt( (errors500Data/entries500Data)*(errors500Data/entries500Data) + (errors500Background/entries500Background)*(errors500Background/entries500Background) );
  Float_t errors750DataOverBackground = entries750DataOverBackground * sqrt( (errors750Data/entries750Data)*(errors750Data/entries750Data) + (errors750Background/entries750Background)*(errors750Background/entries750Background) );
  Float_t errors1000DataOverBackground = entries1000DataOverBackground * sqrt( (errors1000Data/entries1000Data)*(errors1000Data/entries1000Data) + (errors1000Background/entries1000Background)*(errors1000Background/entries1000Background) );
  Float_t errors1250DataOverBackground = entries1250DataOverBackground * sqrt( (errors1250Data/entries1250Data)*(errors1250Data/entries1250Data) + (errors1250Background/entries1250Background)*(errors1250Background/entries1250Background) );

  printf("%s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s \n","DataOverBackground & ",entriesDataOverBackground," & "
	 ,entries200DataOverBackground," & "
	 ,entries500DataOverBackground," & "
	 ,entries750DataOverBackground," & "
	 ,entries1000DataOverBackground," & "
	 ,entries1250DataOverBackground," \\ "
	 );	

  printf("%s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s \n","stat & ",entriesDataOverBackground * sqrt(1./entriesData + 1./entriesBackground)," & "
	 ,entries200DataOverBackground * sqrt(1./entries200Data + 1./entries200Background)," & "
	 ,entries500DataOverBackground * sqrt(1./entries500Data + 1./entries500Background)," & "
	 ,entries750DataOverBackground * sqrt(1./entries750Data + 1./entries750Background)," & "
	 ,entries1000DataOverBackground * sqrt(1./entries1000Data + 1./entries1000Background)," & "
	 ,entries1250DataOverBackground * sqrt(1./entries1250Data + 1./entries1250Background)," \\ "
	 );

  printf("%s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s %.2f %s \n","syst & ",errorsDataOverBackground," & "
	 ,errors200DataOverBackground," & "
	 ,errors500DataOverBackground," & "
	 ,errors750DataOverBackground," & "
	 ,errors1000DataOverBackground," & "
	 ,errors1250DataOverBackground," \\ "
	 );



}//end of method MakeYieldsTable




































































void MakeYieldsTableForMassRanges(TString Sample, Float_t lumiNumber, std::vector<Float_t> minmasses, std::vector<Float_t> maxmasses)
{

  gStyle->SetOptStat("ourmei");

  TString histogramdata = HistogramFileLocation+Sample+"/histograms_"+Sample+".root";
  TFile* fdatahists = TFile::Open(histogramdata.Data());

  TString histogramJetJet =HistogramFileLocation+Sample+"/histograms_"+Sample+"_JetJet.root";
  TString histogramGammaJet =HistogramFileLocation+Sample+"/histograms_"+Sample+"_GammaJet.root";	
  TFile* fJetJethists=TFile::Open(histogramJetJet.Data());
  TFile* fGammaJethists=TFile::Open(histogramGammaJet.Data());

  //HARDCODED
  TString histogramMC =HistogramFileLocation+"/diphoton_tree_MC_all/histograms_diphoton_tree_MC_all.root";
  TFile* fMChists=TFile::Open(histogramMC.Data());

  TH1F* histosmc;
  TH1F* histosdata;
  TH1F* histosJetJet;
  TH1F* histosJetJetUpperError;
  TH1F* histosJetJetLowerError;
  TH1F* histosGammaJet;
  TH1F* histosGammaJetUpperError;
  TH1F* histosGammaJetLowerError;
  

  cout<<"Getting data histogram "<<endl;
  histosdata=(TH1F*)fdatahists->Get("h_Diphoton_Minv_FineBinning");
  cout<<"histo "<<histosdata->GetName()<<" "<<histosdata->GetEntries()<<" entries"<<histosdata->Integral()<<" (integral)"<<endl;

  cout<<"Getting MC histogram "<<endl;
  histosmc= (TH1F*)fMChists->Get("h_Diphoton_Minv_FineBinning");
  histosmc->Scale(lumiNumber);
  cout<<"histo "<<histosmc->GetName()<<" "<<histosmc->GetEntries()<<" entries "<<histosmc->Integral()<<" (integral)"<<endl;

  cout<<"Getting Jet+Jet histograms "<<endl;
  histosJetJet=(TH1F*)fJetJethists->Get("h_JetJet_minv_FineBinning");
  histosJetJetUpperError=(TH1F*)fJetJethists->Get("h_JetJet_minv_UpperError");
  histosJetJetLowerError=(TH1F*)fJetJethists->Get("h_JetJet_minv_LowerError");
  cout<<"histo "<<histosJetJet->GetName()<<" "<<histosJetJet->GetEntries()<<" entries "<<histosJetJet->Integral()<<" (integral)"<<endl;
  cout<<"histo "<<histosJetJetUpperError->GetName()<<" "<<histosJetJetUpperError->GetEntries()<<" entries "<<histosJetJetUpperError->Integral()<<" (integral)"<<endl;
  cout<<"histo "<<histosJetJetLowerError->GetName()<<" "<<histosJetJetLowerError->GetEntries()<<" entries "<<histosJetJetLowerError->Integral()<<" (integral)"<<endl;

  cout<<"Getting Gamma+Jet histograms "<<endl;
  histosGammaJet=(TH1F*)fGammaJethists->Get("h_GammaJet_minv_FineBinning");
  histosGammaJetUpperError=(TH1F*)fGammaJethists->Get("h_GammaJet_minv_UpperError");
  histosGammaJetLowerError=(TH1F*)fGammaJethists->Get("h_GammaJet_minv_LowerError");
  cout<<"histo "<<histosGammaJet->GetName()<<" "<<histosGammaJet->GetEntries()<<" entries "<<histosGammaJet->Integral()<<" (integral)"<<endl;
  cout<<"histo "<<histosGammaJetUpperError->GetName()<<" "<<histosGammaJetUpperError->GetEntries()<<" entries "<<histosGammaJetUpperError->Integral()<<" (integral)"<<endl;
  cout<<"histo "<<histosGammaJetLowerError->GetName()<<" "<<histosGammaJetLowerError->GetEntries()<<" entries "<<histosGammaJetLowerError->Integral()<<" (integral)"<<endl;

  //special part to calculate integrals in different mass ranges

  //They all have the same binning
  //Let's take the data one
  std::vector<Int_t> binminnumbers;
  std::vector<Int_t> binmaxnumbers;

  Int_t nbinsX = histosdata->GetNbinsX();
  Int_t binnr = nbinsX;

  cout<<"---------------Mass intervals to seek----------------"<<endl;
  for(int index = 0;index < minmasses.size();index++){
    cout<<"interval ["<<minmasses[index]<<","<<maxmasses[index]<<"]"<<endl;
  }

  for(unsigned int index = 0;index < minmasses.size();index++) {
    binnr = nbinsX;
    for(int nbin = 0; nbin <  histosdata->GetNbinsX(); nbin++) {
      Float_t binlowedge = histosdata->GetBinLowEdge(nbin);
      if(binlowedge >= minmasses[index] && binnr == nbinsX ) {binnr = nbin;}
    }
    //cout<<" binnr "<<binnr<<endl;
    binminnumbers.push_back(binnr);
  }

  for(unsigned int index = 0;index < maxmasses.size();index++) {
    binnr = nbinsX;
    for(int nbin = 0; nbin <  histosdata->GetNbinsX(); nbin++) {
      Float_t binlowedge = histosdata->GetBinLowEdge(nbin);
      if(binlowedge >= maxmasses[index] && binnr == nbinsX ) {binnr = nbin;}
    }
    //cout<<" binnr-1 "<<binnr-1<<endl;
    binmaxnumbers.push_back(binnr-1);
  }

  //Once we have the bin numbers
  //we compute the integrals
  //and their stat uncertainties
  //and their syst uncertainties

  for(unsigned int index = 0;index < binminnumbers.size();index++) {
    cout<<binminnumbers[index]<<"---"<<binmaxnumbers[index]<<endl;
  }

  //For each contribution, we take nbinsX+1 to have the overflow

  std::vector<Float_t> entriesJetJet;
  std::vector<Float_t> errorsJetJet;

  for(unsigned int index=0;index<binminnumbers.size();index++){
    entriesJetJet.push_back(histosJetJet->Integral(binminnumbers[index],binmaxnumbers[index]));
    errorsJetJet.push_back(histosJetJetUpperError->Integral(binminnumbers[index],binmaxnumbers[index]));
  }

  printf("%s","JetJet #entries = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",entriesJetJet[index]);
  }
  printf("}; \n");
	
  printf("%s","JetJet stat = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",sqrt(entriesJetJet[index]));
  }
  printf("}; \n");
	
  printf("%s","JetJet syst = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",errorsJetJet[index]);
  }
  printf("}; \n");


  std::vector<Float_t> entriesGammaJet;
  std::vector<Float_t> errorsGammaJet;

  for(unsigned int index=0;index<binminnumbers.size();index++){
    entriesGammaJet.push_back(histosGammaJet->Integral(binminnumbers[index],binmaxnumbers[index]));
    errorsGammaJet.push_back(histosGammaJetUpperError->Integral(binminnumbers[index],binmaxnumbers[index]));
  }

  printf("%s","GammaJet #entries = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",entriesGammaJet[index]);
  }
  printf("}; \n");
	
  printf("%s","GammaJet stat = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",sqrt(entriesGammaJet[index]));
  }
  printf("}; \n");
	
  printf("%s","GammaJet syst = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",errorsGammaJet[index]);
  }
  printf("}; \n");



  std::vector<Float_t> entriesFake;
  std::vector<Float_t> errorsFake;

  for(unsigned int index=0;index<binminnumbers.size();index++){
    entriesFake.push_back(entriesGammaJet[index] + entriesJetJet[index]);
    errorsFake.push_back(sqrt( (errorsGammaJet[index] * errorsGammaJet[index]) + (errorsJetJet[index] * errorsJetJet[index]) ));
  }

  printf("%s","Fake #entries = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",entriesFake[index]);
  }
  printf("}; \n");
	
  printf("%s","Fake stat = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",sqrt(entriesFake[index]));
  }
  printf("}; \n");
	
  printf("%s","Fake syst = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",errorsFake[index]);
  }
  printf("}; \n");



  std::vector<Float_t> entriesMC;
  std::vector<Float_t> errorsMC;

  for(unsigned int index=0;index<binminnumbers.size();index++){
    entriesMC.push_back(histosmc->Integral(binminnumbers[index],binmaxnumbers[index]));
    errorsMC.push_back(0.);
  }

  printf("%s","MC #entries = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",entriesMC[index]);
  }
  printf("}; \n");
	
  printf("%s","MC stat = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",sqrt(entriesMC[index]));
  }
  printf("}; \n");
	
  printf("%s","MC syst = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",errorsMC[index]);
  }
  printf("}; \n");



  std::vector<Float_t> entriesBackground;
  std::vector<Float_t> errorsBackground;

  for(unsigned int index=0;index<binminnumbers.size();index++){
    entriesBackground.push_back(entriesGammaJet[index] + entriesJetJet[index] + entriesMC[index]);
    errorsBackground.push_back(sqrt( (errorsGammaJet[index] * errorsGammaJet[index]) + (errorsJetJet[index] * errorsJetJet[index]) + (errorsMC[index] * errorsMC[index]) ));
  }

  printf("%s","Background #entries = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",entriesBackground[index]);
  }
  printf("}; \n");
	
  printf("%s","Background stat = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",sqrt(entriesBackground[index]));
  }
  printf("}; \n");
	
  printf("%s","Background syst = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",errorsBackground[index]);
  }
  printf("}; \n");




  std::vector<Float_t> entriesData;
  std::vector<Float_t> errorsData;

  for(unsigned int index=0;index<binminnumbers.size();index++){
    entriesData.push_back(histosdata->Integral(binminnumbers[index],binmaxnumbers[index]));
    errorsData.push_back(0.);
  }

  printf("%s","Data #entries = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",entriesData[index]);
  }
  printf("}; \n");
	
  printf("%s","Data stat = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",sqrt(entriesData[index]));
  }
  printf("}; \n");
	
  printf("%s","Data syst = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",errorsData[index]);
  }
  printf("}; \n");

  std::vector<Float_t> entriesDataOverBackground;
  std::vector<Float_t> errorsDataOverBackground;

  for(unsigned int index=0;index<binminnumbers.size();index++){
    entriesDataOverBackground.push_back(entriesData[index]/entriesBackground[index]);
    errorsDataOverBackground.push_back(entriesDataOverBackground[index] * sqrt( (errorsData[index]/entriesData[index])*(errorsData[index]/entriesData[index]) + (errorsBackground[index]/entriesBackground[index])*(errorsBackground[index]/entriesBackground[index]) ));
  }

  printf("%s","Data/Background #entries = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",entriesDataOverBackground[index]);
  }
  printf("}; \n");
	
  printf("%s","Data/Background stat = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",sqrt(entriesDataOverBackground[index]));
  }
  printf("}; \n");
	
  printf("%s","Data/Background syst = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",errorsDataOverBackground[index]);
  }
  printf("}; \n");



}//end of method MakeYieldsTable













































void MakeYieldsTableForMassRangesForSignal(std::vector<TString> Samples, std::vector<Float_t> minmasses, std::vector<Float_t> maxmasses)
{

  gStyle->SetOptStat("ourmei");

  std::vector<TFile*> fsamplehists;
  for(int index=0;index<Samples.size();index++){
    TString histogramsample = HistogramFileLocation+Samples[index]+"/histograms_"+Samples[index]+".root";
    TFile *filetempo = TFile::Open(histogramsample.Data());
    fsamplehists.push_back(filetempo);
  }

  std::vector<TH1F*> histossample;
  for(int index=0;index<Samples.size();index++){
    cout<<"Getting sample histogram "<<endl;
    TH1F *histotempo = (TH1F*)fsamplehists[index]->Get("h_Diphoton_Minv_FineBinning");
    histossample.push_back(histotempo);
    cout<<"histo "<<histossample[index]->GetName()<<" "
	<<histossample[index]->GetEntries()<<" entries"
	<<histossample[index]->Integral()<<" (integral)"
	<<endl;
  }

  //special part to calculate integrals in different mass ranges

  //They all have the same binning
  //Let's take the sample one
  std::vector<Int_t> binminnumbers;
  std::vector<Int_t> binmaxnumbers;

  Int_t nbinsX = histossample[0]->GetNbinsX();
  Int_t binnr = nbinsX;

  cout<<"---------------Mass intervals to seek----------------"<<endl;
  for(int index = 0;index < minmasses.size();index++){
    cout<<"interval ["<<minmasses[index]<<","<<maxmasses[index]<<"]"<<endl;
  }

  for(unsigned int index = 0;index < minmasses.size();index++) {
    binnr = nbinsX;
    for(int nbin = 0; nbin <  histossample[index]->GetNbinsX(); nbin++) {
      Float_t binlowedge = histossample[index]->GetBinLowEdge(nbin);
      if(binlowedge >= minmasses[index] && binnr == nbinsX ) {binnr = nbin;}
    }
    //cout<<" binnr "<<binnr<<endl;
    binminnumbers.push_back(binnr);
  }

  for(unsigned int index = 0;index < maxmasses.size();index++) {
    binnr = nbinsX;
    for(int nbin = 0; nbin <  histossample[index]->GetNbinsX(); nbin++) {
      Float_t binlowedge = histossample[index]->GetBinLowEdge(nbin);
      if(binlowedge >= maxmasses[index] && binnr == nbinsX ) {binnr = nbin;}
    }
    //cout<<" binnr-1 "<<binnr-1<<endl;
    binmaxnumbers.push_back(binnr-1);
  }

  //Once we have the bin numbers
  //we compute the integrals
  //and their stat uncertainties
  //and their syst uncertainties

  for(unsigned int index = 0;index < binminnumbers.size();index++) {
    cout<<binminnumbers[index]<<"---"<<binmaxnumbers[index]<<endl;
  }

  //For each contribution, we take nbinsX+1 to have the overflow

  std::vector<Float_t> entriesSamples;
  std::vector<Float_t> errorsSamples;

  for(unsigned int index=0;index<binminnumbers.size();index++){
    entriesSamples.push_back(histossample[index]->Integral(binminnumbers[index],binmaxnumbers[index]));
    errorsSamples.push_back(0.);
  }

  printf("%s","Samples #entries = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",entriesSamples[index]);
  }
  printf("}; \n");
	
  printf("%s","Samples stat = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",sqrt(entriesSamples[index]));
  }
  printf("}; \n");
	
  printf("%s","Samples syst = {");
  for(unsigned int index=0; index < binmaxnumbers.size(); index ++) {
    printf("%.5f,",errorsSamples[index]);
  }
  printf("}; \n");


}//end of method MakeYieldsTable
































void CalculateYieldsInMassRanges(Int_t couplingValue = 1, Int_t numsigmas = 2)
{

  cout<<"Entering CalculateYieldsInMassRanges with parameters: "
      <<" couplingValue (x 100) "<<couplingValue
      <<" numsigmas "<<numsigmas
      <<endl;

  const Int_t n1=7;
  const Int_t n5=5;
  const Int_t n10=6;

  std::vector<TString> sampleNames;

  Float_t mass1[n1] = {750, 1000, 1250, 1500, 1750, 2000, 3000 };
  Float_t mass5[n5] = { 1750, 2000, 2500, 2750, 3000 };
  Float_t mass10[n10] = { 2250, 2500, 2750, 3000, 3250, 3500 };
  
  Float_t sigma1[n1] = { 5.2041, 6.50326, 8.53237, 10.0917, 11.5976, 13.6845, 22.4852 };
  Float_t sigma5[n5] = { 13.7476, 16.795, 20.4539, 22.7374, 25.2982 };    
  Float_t sigma10[n10] = { 26.3803, 30.8038, 35.9283, 38.7399, 41.8178, 40.2991 };

  std::vector<TString> sampleNames1;
  sampleNames1.push_back("diphoton_tree_RSGravToGG_kMpl-001_M-750_TuneZ2star_8TeV-pythia6_merged");
  sampleNames1.push_back("diphoton_tree_RSGravToGG_kMpl-001_M-1000_TuneZ2star_8TeV-pythia6_merged");
  sampleNames1.push_back("diphoton_tree_RSGravToGG_kMpl-001_M-1250_TuneZ2star_8TeV-pythia6_merged");
  sampleNames1.push_back("diphoton_tree_RSGravToGG_kMpl-001_M-1500_TuneZ2star_8TeV-pythia6_merged");
  sampleNames1.push_back("diphoton_tree_RSGravToGG_kMpl-001_M-1750_TuneZ2star_8TeV-pythia6_merged");
  sampleNames1.push_back("diphoton_tree_RSGravToGG_kMpl-001_M-2000_TuneZ2star_8TeV-pythia6_merged");
  //sampleNames1.push_back("diphoton_tree_RSGravToGG_kMpl-001_M-2250_TuneZ2star_8TeV-pythia6_merged");
  //sampleNames1.push_back("diphoton_tree_RSGravToGG_kMpl-001_M-2500_TuneZ2star_8TeV-pythia6_merged");
  sampleNames1.push_back("diphoton_tree_RSGravToGG_kMpl-001_M-3000_TuneZ2star_8TeV-pythia6_merged");

  std::vector<TString> sampleNames5;
  //sampleNames5.push_back("diphoton_tree_RSGravToGG_kMpl-005_M-1250_TuneZ2star_8TeV-pythia6_merged");
  sampleNames5.push_back("diphoton_tree_RSGravToGG_kMpl-005_M-1750_TuneZ2star_8TeV-pythia6_merged");
  sampleNames5.push_back("diphoton_tree_RSGravToGG_kMpl-005_M-2000_TuneZ2star_8TeV-pythia6_merged");
  //sampleNames5.push_back("diphoton_tree_RSGravToGG_kMpl-005_M-2250_TuneZ2star_8TeV-pythia6_merged");
  sampleNames5.push_back("diphoton_tree_RSGravToGG_kMpl-005_M-2500_TuneZ2star_8TeV-pythia6_merged");
  sampleNames5.push_back("diphoton_tree_RSGravToGG_kMpl-005_M-2750_TuneZ2star_8TeV-pythia6_merged");
  sampleNames5.push_back("diphoton_tree_RSGravToGG_kMpl-005_M-3000_TuneZ2star_8TeV-pythia6_merged");

  std::vector<TString> sampleNames10;
  //sampleNames10.push_back("diphoton_tree_RSGravToGG_kMpl-01_M-1500_TuneZ2star_8TeV-pythia6_merged");
  //sampleNames10.push_back("diphoton_tree_RSGravToGG_kMpl-01_M-1750_TuneZ2star_8TeV-pythia6_merged");
  //sampleNames10.push_back("diphoton_tree_RSGravToGG_kMpl-01_M-2000_TuneZ2star_8TeV-pythia6_merged");
  sampleNames10.push_back("diphoton_tree_RSGravToGG_kMpl-01_M-2250_TuneZ2star_8TeV-pythia6_merged");
  sampleNames10.push_back("diphoton_tree_RSGravToGG_kMpl-01_M-2500_TuneZ2star_8TeV-pythia6_merged");
  sampleNames10.push_back("diphoton_tree_RSGravToGG_kMpl-01_M-2750_TuneZ2star_8TeV-pythia6_merged");
  sampleNames10.push_back("diphoton_tree_RSGravToGG_kMpl-01_M-3000_TuneZ2star_8TeV-pythia6_merged");
  sampleNames10.push_back("diphoton_tree_RSGravToGG_kMpl-01_M-3250_TuneZ2star_8TeV-pythia6_merged");
  sampleNames10.push_back("diphoton_tree_RSGravToGG_kMpl-01_M-3500_TuneZ2star_8TeV-pythia6_merged");

  std::vector<Float_t> masses;
  std::vector<Float_t> sigmas;

  if(couplingValue == 1){
    for(int index = 0;index < n1;index++){
      masses.push_back(mass1[index]);
      sigmas.push_back(sigma1[index]);
      sampleNames.push_back(sampleNames1[index]);
    }
  }

  if(couplingValue == 5){
    for(int index = 0;index < n5;index++){
      masses.push_back(mass5[index]);
      sigmas.push_back(sigma5[index]);
      sampleNames.push_back(sampleNames5[index]);
    }
  }

  if(couplingValue == 10){
    for(int index = 0;index < n10;index++){
      masses.push_back(mass10[index]);
      sigmas.push_back(sigma10[index]);
      sampleNames.push_back(sampleNames10[index]);
    }
  }

  std::vector<Float_t> minmasses;
  std::vector<Float_t> maxmasses;

  for(unsigned int index=0;index<masses.size();index++){
    minmasses.push_back(masses[index] - numsigmas * sigmas[index]);
    maxmasses.push_back(masses[index] + numsigmas * sigmas[index]);
  }


  //MakeYieldsTableForMassRanges("ExoDiPhotonAnalyzer_DataABC",10252,minmasses,maxmasses;)
  //MakeYieldsTableForMassRanges("ExoDiPhotonAnalyzer_DataABC",19620,minmasses,maxmasses);

  //For the signal points
  MakeYieldsTableForMassRangesForSignal(sampleNames,minmasses,maxmasses);
}





