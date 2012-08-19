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

TString tfHists[12] = {"h_FakeRate_tf_minv","h_FakeRate_tf_qt","h_FakeRate_tf_deltaPhi","h_FakeRate_tf_deltaEta","h_FakeRate_tf_deltaR","h_FakeRate_tf_cosThetaStar","h_FakeRate_tf_pt1","h_FakeRate_tf_eta1","h_FakeRate_tf_phi1","h_FakeRate_tf_pt2","h_FakeRate_tf_eta2","h_FakeRate_tf_phi2"};

TString ffHists[12] = {"h_FakeRate_ff_minv","h_FakeRate_ff_qt","h_FakeRate_ff_deltaPhi","h_FakeRate_ff_deltaEta","h_FakeRate_ff_deltaR","h_FakeRate_ff_cosThetaStar","h_FakeRate_ff_pt1","h_FakeRate_ff_eta1","h_FakeRate_ff_phi1","h_FakeRate_ff_pt2","h_FakeRate_ff_eta2","h_FakeRate_ff_phi2"};

TString ftHists[12] = {"h_FakeRate_ft_minv","h_FakeRate_ft_qt","h_FakeRate_ft_deltaPhi","h_FakeRate_ft_deltaEta","h_FakeRate_ft_deltaR","h_FakeRate_ft_cosThetaStar","h_FakeRate_ft_pt1","h_FakeRate_ft_eta1","h_FakeRate_ft_phi1","h_FakeRate_ft_pt2","h_FakeRate_ft_eta2","h_FakeRate_ft_phi2"};


TString GammaJetHists[12] = {"h_GammaJet_minv","h_GammaJet_qt","h_GammaJet_deltaPhi","h_GammaJet_deltaEta","h_GammaJet_deltaR","h_GammaJet_cosThetaStar","h_GammaJet_pt1","h_GammaJet_eta1", "h_GammaJet_phi1","h_GammaJet_pt2","h_GammaJet_eta2","h_GammaJet_phi2"};

TString JetJetHists[12] = {"h_JetJet_minv","h_JetJet_qt","h_JetJet_deltaPhi","h_JetJet_deltaEta","h_JetJet_deltaR","h_JetJet_cosThetaStar","h_JetJet_pt1","h_JetJet_eta1", "h_JetJet_phi1","h_JetJet_pt2","h_JetJet_eta2","h_JetJet_phi2"};

TString nameHists[43] = {"h_Photon1_pt","h_Photon1_pt_log","h_Photon1_eta","h_Photon1_phi","h_Photon2_pt","h_Photon2_pt_log","h_Photon2_eta","h_Photon2_phi","h_Diphoton_Minv","h_Diphoton_Minv_log", "h_Diphoton_qt","h_Diphoton_deltaR","h_Diphoton_deltaEta","h_Diphoton_cosThetaStar","h_Diphoton_deltaPhi","h_Vtx_Nvtx","h_Vtx_vx","h_Vtx_vy","h_Vtx_vz","h_Photon1_sigmaIetaIeta","h_Photon1_sigmaEtaEta","h_Photon1_hadOverEm","h_Photon1_trkIsoSumPtHollow04","h_Photon1_trkIsoNtrksHollow04","h_Pho\
ton1_hcalIso04","h_Photon1_ecalIso04","h_Photon1_detEta","h_Photon2_sigmaIetaIeta","h_Photon2_sigmaEtaEta","h_Photon2_hadOverEm","h_Photon2_trkIsoSumPtHollow04","h_Photon2_trkIsoNtrksHollow04","h_Photon2_hcalIso04","h_Photon2_ecalIso04","h_Photon2_detEta","h_Diphoton_qt_log","h_Photon1_hadOverEm_log","h_Photon2_hadOverEm","h_Photon2_trkIsoSumPtHollow04_log","h_Photon1_trkIsoSumPtHollow04_log","h_Photon2_hcalIso04_log","h_Photon1_hcalIso04_log","h_Cumalitive_DiphotonMinv"};


TString fNames[12] = {"h_Diphoton_Minv","h_Diphoton_qt","h_Diphoton_deltaPhi","h_Diphoton_deltaEta","h_Diphoton_deltaR","h_Diphoton_cosThetaStar","h_Photon1_pt","h_Photon1_eta","h_Photon1_phi","h_Photon2_pt","h_Photon2_eta","h_Photon2_phi"};

TString TreeFileLocation = "/afs/cern.ch/work/j/jcarson/private/DiPhotonTrees/";
TString HistogramFileLocation = "/afs/cern.ch/work/j/jcarson/private/DiPhotonTrees/histograms/";

double xsec[4]={25.48,0.029038,15.54,0.0011805};
//in pb
int inputfiles=4 ;

int ngenevents[4]={500254,500038,500050,500352};


TFile* ftemp[4];


void CreateHistogramFiles(TString Sample = "Diphoton", TString SampleType = "data", TString JSON = "BOPB")
{
 
  TString outName;
  TFile *outfilename;

  TString  inputfile= TreeFileLocation+"/"+Sample+".root";
  
 
  TChain *chain_tt = new TChain("diphotonAnalyzer/fTree");
    chain_tt->Add(inputfile.Data());

  cout << "TT entries = " << chain_tt->GetEntries() <<endl;
  
  
  PlottingCodeLoop* treereader = new PlottingCodeLoop(chain_tt);
  treereader->_cutPhoton1pt = 70;
  treereader->_cutPhoton2pt = 70;
  treereader->_fakeStatus="TightTight";
  treereader->_JSON=JSON.Data();
  cout<<JSON.Data()<<endl;
  treereader->_SampleType=SampleType.Data();
  cout<<SampleType.Data()<<endl;
  outName = TString::Format("histograms_%s.root", Sample.Data());
  outfilename = new TFile(HistogramFileLocation+Sample.Data()+"/"+outName.Data(),"RECREATE");
  treereader->_outputfile = outfilename;
 treereader->Loop();
  
  cout << "Here 1" <<endl;


  if (SampleType=="data"){
    cout << "Here 2" <<endl;
    TChain *chain_tf = new TChain("diphotonAnalyzer/fTightFakeTree");
      chain_tf->Add(inputfile.Data());
    cout << "TF entries = " << chain_tf->GetEntries() <<endl;
    
    TChain *chain_ft = new TChain("diphotonAnalyzer/fFakeTightTree");
      chain_ft->Add(inputfile.Data());

    
    
    TChain *chain_ff = new TChain("diphotonAnalyzer/fFakeFakeTree");
      chain_ff->Add(inputfile.Data());

    
    
    PlottingCodeLoop* treereaderTF = new PlottingCodeLoop(chain_tf);
    treereaderTF->_fakeStatus = "TightFake";    
    treereaderTF->_cutPhoton1pt = 70;
    treereaderTF->_cutPhoton2pt = 70;
    treereader->_JSON=JSON.Data();
    treereader->_SampleType=SampleType.Data();
    outName = TString::Format("histograms_%s_TF.root",Sample.Data());
    outfilename = new TFile(HistogramFileLocation+Sample+"/"+outName,"RECREATE");
    treereaderTF->_outputfile = outfilename;
    
    treereaderTF->Loop(); 
    
    PlottingCodeLoop* treereaderFT = new PlottingCodeLoop(chain_ft);
    treereaderFT->_fakeStatus = "FakeTight";    
    treereaderFT->_cutPhoton1pt = 70;
    treereaderFT->_cutPhoton2pt = 70;
    treereader->_JSON=JSON.Data();
    treereader->_SampleType=SampleType.Data();
    outName = TString::Format("histograms_%s_FT.root",Sample.Data());
    outfilename = new TFile(HistogramFileLocation+Sample+"/"+outName,"RECREATE");
    treereaderFT->_outputfile = outfilename;
    treereaderFT->Loop();     
    
    PlottingCodeLoop* treereaderFF = new PlottingCodeLoop(chain_ff);
    treereaderFF->_fakeStatus = "FakeFake";    
    treereaderFF->_cutPhoton1pt = 70;
    treereaderFF->_cutPhoton2pt = 70;
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

} 



void makeplots(TString Sample = "Diphoton",TString lumi = "300", TString JSONFile="April20.json" ,TString SampleType = "data")
{
  char canvasname[42];
  TCanvas*  c[42];
  TString histogramfileinput = HistogramFileLocation+"histograms_"+Sample+".root";
  TFile* fhists = TFile::Open(histogramfileinput.Data());
  fhists->cd();
  cout<<"ROLLTIDE"<<endl;
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
  TH1F* histos[42];
  for(int i=0;i<42;i++)
    {
      cout<<i<<endl;
      sprintf(canvasname,"c%i",i);
      cout<<canvasname<<endl;
      histos[i] = (TH1F*)fhists->Get(nameHists[i].Data());
      c[i] = new TCanvas(canvasname, canvasname, 800, 600);
      cout<<histos[i]<<endl;
      c[i]->cd();
      histos[i]->GetXaxis()->CenterTitle();
      histos[i]->GetYaxis()->CenterTitle();
      histos[i]->SetMarkerStyle(20);
      histos[i]->SetMarkerColor(1);
      histos[i]->Draw("HIST EP");




      //      if (((nameHists[i]=="h_Photon1_pt_log") || (nameHists[i]=="h_Photon2_pt_log") || (nameHists[i]=="h_Diphoton_Minv_log")))
      if (nameHists[i].Contains("_log"))
	{
	  c[i]->SetLogy(1);
	}
      LumiLabel->Draw();
      c[i]->SaveAs(nameHists[i]+"_"+Sample+".C");
      c[i]->SaveAs(nameHists[i]+"_"+Sample+".png");
      c[i]->SaveAs(nameHists[i]+"_"+Sample+".pdf");
    }
}






void MakeCombinedMCHistos()
{

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat("ourme");
  const double lumi = 1.0; // units of pb-1
  //const

  TString inputmcfiles[inputfiles];
  inputmcfiles[0]= HistogramFileLocation+"diphoton_tree_DiPhotonBorn_Pt25to250_Summer12/histograms_diphoton_tree_DiPhotonBorn_Pt25to250_Summer12.root";
  inputmcfiles[1]= HistogramFileLocation+"diphoton_tree_DiPhotonBorn_Pt250toInf_Summer12/histograms_diphoton_tree_DiPhotonBorn_Pt250toInf_Summer12.root";
  inputmcfiles[2]= HistogramFileLocation+"diphoton_tree_DiPhotonBox_Pt25to250_Summer12/histograms_diphoton_tree_DiPhotonBox_Pt25to250_Summer12.root";	
  inputmcfiles[3]= HistogramFileLocation+"diphoton_tree_DiPhotonBox_Pt250toInf_Summer12/histograms_diphoton_tree_DiPhotonBox_Pt250toInf_Summer12.root";

    TString filenames[inputfiles];
    for(int ifile=0;ifile<inputfiles;ifile++) {
      filenames[ifile] = inputmcfiles[ifile];
      cout << filenames[ifile]<<endl;
    }

    for(int ifile=0;ifile<inputfiles;ifile++) {
      ftemp[ifile] = TFile::Open(filenames[ifile].Data());
    }
    int nHists = 42;
    TH1F* indHistos[nHists][inputfiles];
    Double_t totalError = 0;

    for(int ifile=0;ifile<inputfiles;ifile++) {
      for (int i=0; i<nHists; i++) {
	indHistos[i][ifile] = (TH1F*)ftemp[ifile]->Get(nameHists[i].Data());

	if (i==1) {
	  //  error N_tot = sqrt(C1^2 * n1 + C2^2 * n2 + ...)
	  Double_t ci = xsec[ifile]/ngenevents[ifile];
	  //      totalError = totalError + ci*ci*(indHistos[i][ifile]->Integral());
	  totalError = totalError + ci*TMath::Sqrt(indHistos[i][ifile]->Integral());
	}

	indHistos[i][ifile]->Sumw2();
	indHistos[i][ifile]->Scale((1.0/ngenevents[ifile])*xsec[ifile]*lumi);
	cout << "xsec orig lumi " << xsec[ifile] << " " << ngenevents[ifile] << " " << lumi << endl;
	//      cout << indHistos[i][ifile]->GetEntries() << endl;
      }
    }

    TH1F* allHistos[nHists];

    // instantiate the sum histos with histo in first file
    for (int i=0; i<nHists; i++) {
      allHistos[i] = indHistos[i][0];
    }

    for(int ifile=0;ifile<inputfiles;ifile++) {
      for (int i=0; i<nHists; i++) {
	if (ifile>0){ allHistos[i]->Add(indHistos[i][ifile]);} // note we exclude first file since we used it to instantiate the sum histo
      }
    }
    TString MCALLOut = HistogramFileLocation+"/diphoton_tree_MC_all/";
      TFile fout(MCALLOut+"histograms_diphoton_tree_MC_all.root","RECREATE");

    for (int i=0; i<nHists; i++) {
      allHistos[i]->Write();
    }
    fout.Close();

    makeplots("diphoton_tree_MC_all","1","MC","mc");

}


TH1F* MakeChists(TH1F* hist){
  TH1F* cHist = (TH1F*) hist->Clone("cHist");
  cHist->SetDirectory(0);

  int maxBin = hist->GetNbinsX()+1;
  //int maxBin = hist->GetNbinsX();
  for(int binNr=0;binNr<=hist->GetNbinsX();binNr++){
    float nrEntries =hist->Integral(binNr,maxBin);
    cHist->SetBinContent(binNr,nrEntries);
  }
  return cHist;
}



void OverlayMCandData( TString Sample  = "", TString lumi = "1", TString JSONFile = "April20")
{
  TString histogramdata =HistogramFileLocation+Sample+"/histograms_"+Sample+".root";
  TString histogrammc=HistogramFileLocation+"diphoton_tree_MC_all/histograms_diphoton_tree_MC_all.root";
  TFile* fmchists = TFile::Open(histogrammc.Data());
  TFile* fdatahists = TFile::Open(histogramdata.Data());

  char canvasname[42];
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
  cout<<lumiNumber<<endl;
  for(int i=0;i<42;i++)
    {
      cout<<i<<endl;
      sprintf(canvasname,"c%i",i);
      cout<<canvasname<<endl;
      histosdata[i] = (TH1F*)fdatahists->Get(nameHists[i].Data());
      histosdata[i]->SetMinimum(.01);
      histosmc[i] = (TH1F*)fmchists->Get(nameHists[i].Data());
      histosmc[i]->Scale(lumiNumber);
      histosmc[i]->SetMinimum(.01);
      c[i] = new TCanvas(canvasname, canvasname, 800, 600);
      cout<<histosmc[i]<<endl;
      c[i]->cd();
      DataMCLegend[i]= new TLegend(0.70,0.6,0.88,0.75,"","NDC");
      DataMCLegend[i]->SetFillStyle(0);
      DataMCLegend[i]->SetTextSize(0.03);
      histosmc[i]->GetXaxis()->CenterTitle();
      histosmc[i]->GetYaxis()->CenterTitle();
      histosmc[i]->SetFillColor(kBlue);
      histosdata[i]->SetMarkerColor(1);
      histosdata[i]->SetMaximum(1.3*histosdata[i]->GetMaximum());
      histosdata[i]->SetMinimum(.8*histosmc[i]->GetMinimum());
      histosdata[i]->Draw("HIST EP");
      histosmc[i]->Draw("SAME HIST");
      histosdata[i]->Draw("HIST SAME EP");
      DataMCLegend[i]->AddEntry(histosdata[i],"Data","LEP");
      DataMCLegend[i]->AddEntry(histosmc[i],"SM Diphoton","F");
      gPad->RedrawAxis();
      histosdata[i]->Draw("EPSAME");
      LumiLabel->Draw();
      DataMCLegend[i]->Draw();
      //if (((nameHists[i]=="h_Photon1_pt_log") || (nameHists[i]=="h_Photon2_pt_log") || (nameHists[i]=="h_Diphoton_Minv_log")))
      if (nameHists[i].Contains("_log"))
        {
	  c[i]->SetLogy(1);
	}

      c[i]->SaveAs(nameHists[i]+"_"+Sample.Data()+"Overlay"+".C");
      c[i]->SaveAs(nameHists[i]+"_"+Sample.Data()+"Overlay"+".png");
      c[i]->SaveAs(nameHists[i]+"_"+Sample.Data()+"Overlay"+".pdf");
    }
  TCanvas* CummalitiveCanvas = new TCanvas("CummalitiveCanvas", "CummalitiveCanvas", 800, 600);
  CummalitiveCanvas->cd();
  CummalitiveCanvas->SetLogy(1);
  TH1F* h_Diphoton_Minv_log_DataCumm =  MakeChists(histosdata[9]);
  TH1F* h_Diphoton_Minv_log_MCCumm = MakeChists(histosmc[9]);

  h_Diphoton_Minv_log_MCCumm->SetFillColor(kBlue);
  DataMCLegend[1]->Draw();

  TPaveText *LumiLabelCumm = new TPaveText(.55,.75,.75,.9,"NDC");
  LumiLabelCumm->SetTextSize(0.03);
  LumiLabelCumm->SetFillStyle(0);
  LumiLabelCumm->SetBorderSize(0);
  LumiLabelCumm->AddText("CMS Internal");
  TString LuminosityCumm = lumi;
  LumiLabelCumm->AddText(LuminosityCumm.Data());
  LumiLabelCumm->AddText(JSONFile.Data());

  TLegend* DataMCLegendCumm = new TLegend(0.70,0.6,0.88,0.75,"","NDC");
  DataMCLegendCumm->AddEntry( h_Diphoton_Minv_log_DataCumm,"Data","LEP");
  DataMCLegendCumm->AddEntry(h_Diphoton_Minv_log_MCCumm,"SM Diphoton","F");
  DataMCLegendCumm->SetFillColor(0);
  h_Diphoton_Minv_log_DataCumm->SetMaximum(2.5*h_Diphoton_Minv_log_DataCumm->GetMaximum());
  h_Diphoton_Minv_log_DataCumm->SetMinimum(.01);
  h_Diphoton_Minv_log_DataCumm->Draw("EP");
  h_Diphoton_Minv_log_MCCumm->Draw("HISTSAME");
  h_Diphoton_Minv_log_MCCumm->SetMinimum(.01);
  h_Diphoton_Minv_log_DataCumm->SetMaximum(2.5*h_Diphoton_Minv_log_DataCumm->GetMaximum());
  h_Diphoton_Minv_log_DataCumm->Draw("EPSAME");
  gPad->RedrawAxis();
  h_Diphoton_Minv_log_DataCumm->SetMaximum(2.5*h_Diphoton_Minv_log_DataCumm->GetMaximum());
  LumiLabelCumm->Draw();
  DataMCLegendCumm->Draw();
  CummalitiveCanvas->SaveAs(nameHists[43]+Sample.Data()+"Overlay"+".png");
  return;
}



void StitchBackgroundandMC(TString Sample = "2", TString lumi = "1", TString JSONFile = "April20")
{
  TString histogramdata = HistogramFileLocation.Data()+Sample+"/histograms_"+Sample+".root";
  TFile* fdatahists = TFile::Open(histogramdata.Data());
  TH1F* histosmcdminv[42][4];
  TH1F* histosdatadminv[42];

  char canvasname[41];
  TCanvas*  c[42];

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
  cout<<lumiNumber<<endl;

  Color_t fillColors[inputfiles];
  fillColors[0] = kGreen;
  fillColors[1] = kRed;
  fillColors[2] = kBlue;
  fillColors[3] = 6;

  TString inputmcfiles[inputfiles];
  inputmcfiles[0]= HistogramFileLocation+"diphoton_tree_DiPhotonBorn_Pt25to250_Summer12/histograms_diphoton_tree_DiPhotonBorn_Pt25to250_Summer12.root";
  inputmcfiles[1]= HistogramFileLocation+"diphoton_tree_DiPhotonBorn_Pt250toInf_Summer12/histograms_diphoton_tree_DiPhotonBorn_Pt250toInf_Summer12.root";
  inputmcfiles[2]= HistogramFileLocation+"diphoton_tree_DiPhotonBox_Pt25to250_Summer12/histograms_diphoton_tree_DiPhotonBox_Pt25to250_Summer12.root";
  inputmcfiles[3]= HistogramFileLocation+"diphoton_tree_DiPhotonBox_Pt250toInf_Summer12/histograms_diphotn_tree_DiPhotonBox_Pt250toInf_Summer12.root";

    TString filenames[inputfiles];
  for(int ifile=0;ifile<inputfiles;ifile++) {
    filenames[ifile] = inputmcfiles[ifile];
  }
  TString SampleNames[4]={"Born_Pt_25to250","Born_Pt_250toInf","Box_Pt_25to250","Box_Pt_250toInf"};
  TLegend *StackLegend[42];
  THStack *StackMC[42];
  for(int i=0;i<42;i++){
    c[i] = new TCanvas(canvasname, canvasname, 800, 600);
    c[i]->cd();
    StackLegend[i]= new TLegend(0.70,0.6,0.88,0.75,"","NDC");
    StackLegend[i]->SetFillStyle(0);
    StackLegend[i]->SetTextSize(0.03);
    StackMC[i]= new THStack(nameHists[i].Data(),nameHists[i].Data());
    histosdatadminv[i]=(TH1F*)fdatahists->Get(nameHists[i].Data());
    for(int ifile=0;ifile<inputfiles;ifile++) {
      ftemp[ifile] = TFile::Open(filenames[ifile].Data());
      histosmcdminv[ifile][i]= (TH1F*)ftemp[ifile]->Get(nameHists[i].Data());
      histosmcdminv[ifile][i]->SetFillColor(fillColors[ifile]);
      histosmcdminv[ifile][i]->Scale((1.0/ngenevents[ifile])*xsec[ifile]*1);
      histosmcdminv[ifile][i]->Scale(lumiNumber);
      StackMC[i]->Add(histosmcdminv[ifile][i]);
      StackLegend[i]->AddEntry(histosmcdminv[ifile][i],SampleNames[ifile].Data(),"f");
      //    if (((nameHists[i]=="h_Photon1_pt_log") || (nameHists[i]=="h_Photon2_pt_log") || (nameHists[i]=="h_Diphoton_Minv_log")))
      if (nameHists[i].Contains("_log"))
      {
	c[i]->SetLogy(1);
      }
    histosdatadminv[i]->SetMinimum(.8*(StackMC[i]->GetMinimum()));
    StackLegend[i]->AddEntry(histosdatadminv[i],"Data","LEP");
    histosdatadminv[i]->Draw("HIST EP");
    StackMC[i]->Draw("SAME") ;
    histosdatadminv[i]->Draw("SAME EP");
    LumiLabel->Draw();
    StackLegend[i]->Draw();
    gPad->RedrawAxis();


    c[i]->SaveAs(nameHists[i]+"_"+Sample.Data()+"StitchOverlay"+".C");
    c[i]->SaveAs(nameHists[i]+"_"+Sample.Data()+"StitchOverlay"+".png");
    c[i]->SaveAs(nameHists[i]+"_"+Sample.Data()+"StitchOverlay"+".pdf");

  }
}
}

 void fakeratehistos(TString Sample = "Diphoton",TString lumi = "300", TString JSONFile="April20.json" )
  {
    TH1F* h_JetJet_pt1;
    TH1F* h_JetJet_eta1;
    TH1F* h_JetJet_phi1;
    TH1F* h_JetJet_pt2;
    TH1F* h_JetJet_eta2;
    TH1F* h_JetJet_phi2;
    TH1F* h_JetJet_minv;
    TH1F* h_JetJet_qt;
    TH1F* h_JetJet_deltaPhi;
    TH1F* h_JetJet_deltaEta;
    TH1F* h_JetJet_deltaR;
    TH1F* h_JetJet_cosThetaStar;

  
    TH1F* h_GammaJet_pt1;
    TH1F* h_GammaJet_eta1;
    TH1F* h_GammaJet_phi1;
    TH1F* h_GammaJet_pt2;
    TH1F* h_GammaJet_eta2;
    TH1F* h_GammaJet_phi2;
    TH1F* h_GammaJet_minv;
    TH1F* h_GammaJet_qt;
    TH1F* h_GammaJet_deltaPhi;
    TH1F* h_GammaJet_deltaEta;
    TH1F* h_GammaJet_deltaR;
    TH1F* h_GammaJet_cosThetaStar;


    h_GammaJet_pt1  = new TH1F("h_GammaJet_pt1","#gamma_{1} p_{T};#gamma_{1} p_{T}",42,60,900);
    h_GammaJet_eta1 = new TH1F("h_GammaJet_eta1","#gamma_{1} #eta;#gamma_{1} #eta",40,-3.14159,3.14159);
    h_GammaJet_phi1 = new TH1F("h_GammaJet_phi1","#gamma_{1} #phi;#gamma_{1} #phi",40,-3.14159,3.14159);
    h_GammaJet_pt2  = new TH1F("h_GammaJet_pt2","#gamma_{2} p_{T};#gamma_{2} p_{T}",42,60,900);
    h_GammaJet_eta2 = new TH1F("h_GammaJet_eta2","#gamma_{2} #eta;#gamma_{2} #eta",40,-3.14159,3.14159);
    h_GammaJet_phi2 = new TH1F("h_GammaJet_phi2","#gamma_{1} #phi;#gamma_{1} #phi",40,-3.14159,3.14159);
    //   h_GammaJet_minv         = new TH1F("h_GammaJet_minv",        "Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",43,140,1000);
    h_GammaJet_minv         = new TH1F("h_GammaJet_minv",        "Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",89,20,1800);
        h_GammaJet_qt           = new TH1F("h_GammaJet_qt ",         "Diphoton qt;#gamma#gamma qt [GeV]",50,0,600);
    h_GammaJet_deltaPhi     = new TH1F("h_GammaJet_deltaPhi",    "Diphoton #Delta#phi;#gamma#gamma #Delta#phi",90,-3.14159,3.14159);
    h_GammaJet_deltaEta     = new TH1F("h_GammaJet_deltaEta",    "Diphoton #Delta#eta;#gamma#gamma #Delta#eta",40,-6.,6.);
    h_GammaJet_deltaR       = new TH1F("h_GammaJet_deltaR",      "Diphoton #DeltaR; #gamma#gamma #DeltaR",70,0,7.);
    h_GammaJet_cosThetaStar = new TH1F("h_GammaJet_cosThetaStar","Diphoton |cos(#theta *)|; #gamma#gamma |cos#theta*|",20,0,1);



    h_JetJet_pt1  = new TH1F("h_JetJet_pt1","#gamma_{1} p_{T};#gamma_{1} p_{T}",42,60,900);
    h_JetJet_eta1 = new TH1F("h_JetJet_eta1","#gamma_{1} #eta;#gamma_{1} #eta",40,-3.14159,3.14159);
    h_JetJet_phi1 = new TH1F("h_JetJet_phi1","#gamma_{1} #phi;#gamma_{1} #phi",40,-3.14159,3.14159);
    h_JetJet_pt2  = new TH1F("h_JetJet_pt2","#gamma_{2} p_{T};#gamma_{2} p_{T}",42,60,900);
    h_JetJet_eta2 = new TH1F("h_JetJet_eta2","#gamma_{2} #eta;#gamma_{2} #eta",40,-3.14159,3.14159);
    h_JetJet_phi2 = new TH1F("h_JetJet_phi2","#gamma_{1} #phi;#gamma_{1} #phi",40,-3.14159,3.14159);
    //   h_JetJet_minv         = new TH1F("h_JetJet_minv",        "Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",43,140,1000);
    h_JetJet_minv         = new TH1F("h_JetJet_minv",        "Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",89,20,1800);
    h_JetJet_qt           = new TH1F("h_JetJet_qt ",         "Diphoton qt;#gamma#gamma qt [GeV]",50,0,600);
    h_JetJet_deltaPhi     = new TH1F("h_JetJet_deltaPhi",    "Diphoton #Delta#phi;#gamma#gamma #Delta#phi",90,-3.14159,3.14159);
    h_JetJet_deltaEta     = new TH1F("h_JetJet_deltaEta",    "Diphoton #Delta#eta;#gamma#gamma #Delta#eta",40,-6.,6.);
    h_JetJet_deltaR       = new TH1F("h_JetJet_deltaR",      "Diphoton #DeltaR; #gamma#gamma #DeltaR",70,0,7.);
    h_JetJet_cosThetaStar = new TH1F("h_JetJet_cosThetaStar","Diphoton |cos(#theta *)|; #gamma#gamma |cos#theta*|",20,0,1);

    TString histoTFlocation=HistogramFileLocation.Data()+Sample+"/histograms_"+Sample+"_TF.root";
    cout<<histoTFlocation<<endl;
    TString histoFFlocation=HistogramFileLocation.Data()+Sample+"/histograms_"+Sample+"_FF.root";
    TString histoFTlocation=HistogramFileLocation.Data()+Sample+"/histograms_"+Sample+"_FT.root";
    
    TFile* histoFF = TFile::Open(histoFFlocation.Data());
    TFile* histoFT = TFile::Open(histoFTlocation.Data()); 
    TFile* histoTF = TFile::Open(histoTFlocation.Data());  
    
    TString GammaJetoutName = TString::Format("histograms_%s_GammaJet.root",Sample.Data());
    TFile* histoFileGammaJet = new TFile(HistogramFileLocation.Data()+Sample+"/"+GammaJetoutName,"RECREATE");       
    
    TString JetJetoutName = TString::Format("histograms_%s_JetJet.root",Sample.Data());
    TFile* histoFileJetJet = new TFile(HistogramFileLocation.Data()+Sample+"/"+JetJetoutName,"RECREATE");

    int nFakeHists=12;


    
  TH1F* hFF[12];
  TH1F* hFT[12];
  TH1F* hTF[12];
  TH1F* histosGammaJet[12];
  TH1F* histosJetJet[12];
    
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


    TList *list[12];
    TList *listFF[12];
  for(int i=0;i<12;i++) {
    hFF[i]=(TH1F*)histoFF->Get(ffHists[i].Data());
    hFT[i]=(TH1F*)histoFT->Get(ftHists[i].Data());
    hTF[i]=(TH1F*)histoTF->Get(tfHists[i].Data());
    cout<<hFF[i]<<endl;
    histosGammaJet[i]= (TH1F*)GammaJetHistList->At(i);
    histosGammaJet[i]->Sumw2();
    cout<<hFF[i]->GetEntries()<<endl;;  
  histosGammaJet[i]->Add(hTF[i],hFT[i],1,1);
    cout<<histosGammaJet[i]->GetEntries()<<endl;
    histosGammaJet[i]->Add(hFF[i],-2);
     cout<<histosGammaJet[i]<<endl;
    histosJetJet[i]= (TH1F*)JetJetHistList->At(i);
    histosJetJet[i]->Add(hFF[i],1); 

      
      histoFileGammaJet->cd();
      histosGammaJet[i]->Write(); 
      histoFileJetJet->cd();
      histosJetJet[i]->Write();
      
      }
  
  histoFileGammaJet->Close();
  histoFileJetJet->Close();
  

      char canvasnamegj[13];
      char canvasnamejj[13];
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
	  sprintf(canvasnamegj,"cgj%in",in);
	  cgj[in] = new TCanvas(canvasnamegj, canvasnamegj, 800, 600);
	  histosGammaJet[in]->GetXaxis()->CenterTitle();
	  histosGammaJet[in]->GetYaxis()->CenterTitle();
	  histosGammaJet[in]->SetMarkerStyle(20);
	  histosGammaJet[in]->SetMarkerColor(1);
	  histosGammaJet[in]->Draw();
	  sprintf(canvasnamejj,"cjj%in",in);
          cjj[in] = new TCanvas(canvasnamejj, canvasnamejj, 800, 600);
          histosJetJet[in]->GetXaxis()->CenterTitle();
          histosJetJet[in]->GetYaxis()->CenterTitle();
          histosJetJet[in]->SetMarkerStyle(20);
          histosJetJet[in]->SetMarkerColor(1);
          histosJetJet[in]->Draw();
	  if ((strcmp(histosGammaJet[in]->GetName(),"h_GammaJet_minv")==0 )|| (strcmp(histosGammaJet[in]->GetName(),"h_GammaJet_pt1")==0 ) || strcmp(histosGammaJet[in]->GetName(),"h_GammaJet_pt2")==0 ){
	    cgj[in]->SetLogy(1);
          }
          if ((strcmp(histosJetJet[in]->GetName(),"h_JetJet_minv")==0 )|| (strcmp(histosJetJet[in]->GetName(),"h_JetJet_pt1")==0 ) || strcmp(histosGammaJet[in]->GetName(),"h_JetJet_pt2")==0 ){

          cjj[in]->SetLogy(1);
          }
	  LumiLabel->Draw();
	  cgj[in]->SaveAs(GammaJetHists[in]+"_"+Sample+".C");
	  cgj[in]->SaveAs(GammaJetHists[in]+"_"+Sample+".png");
	  cgj[in]->SaveAs(GammaJetHists[in]+"_"+Sample+".pdf");
	  cjj[in]->SaveAs(JetJetHists[in]+"_"+Sample+".C");
          cjj[in]->SaveAs(JetJetHists[in]+"_"+Sample+".png");
          cjj[in]->SaveAs(JetJetHists[in]+"_"+Sample+".pdf");


}
  }

void OverlayMCFake(TString Sample = "2", TString lumi = "1", TString JSONFile = "April20")
{





  TCanvas* CummalitiveCanvas = new TCanvas("CummalitiveCanvas", "CummalitiveCanvas", 800, 600);
  TCanvas* dataoverlayedCumulative = new TCanvas("dataoverlayedCumulative", "Data Overlayed onto Cumulative Plot",800,600);
  THStack *StackCumm;

  TPaveText *LumiLabelCumm = new TPaveText(.55,.75,.75,.9,"NDC");
  LumiLabelCumm->SetTextSize(0.03);
  LumiLabelCumm->SetFillStyle(0);
  LumiLabelCumm->SetBorderSize(0);
  LumiLabelCumm->AddText("CMS Internal");
  TString LuminosityCumm = lumi;
  LumiLabelCumm->AddText(LuminosityCumm.Data());
  LumiLabelCumm->AddText(JSONFile.Data());
  TH1F* h_DatadivBack = new TH1F("h_DatadivBack","Data/Background",89,20,1800); 
  TH1F* h_TotalBackground = new TH1F("h_TotalBackground","Data/Background",89,20,1800);
  TH1F* dataoncumulative = new TH1F("h_TotalBackground"," (Cumulative) Diphoton Invariant Mass;M_{#gamma#gamma} [GeV/c^2]",89,20,1800);
  TString histogramdata = HistogramFileLocation+Sample+"/histograms_"+Sample+".root";
  TFile* fdatahists = TFile::Open(histogramdata.Data());
  TString histogramJetJet =HistogramFileLocation+Sample+"/histograms_"+Sample+"_JetJet.root";
  TString histogramGammaJet =HistogramFileLocation+Sample+"/histograms_"+Sample+"_GammaJet.root";	
  TFile* fJetJethists=TFile::Open(histogramJetJet.Data());
  cout<<fJetJethists<<endl; 
  TFile* fGammaJethists=TFile::Open(histogramGammaJet.Data());
  TString histogramMC =HistogramFileLocation+"/diphoton_tree_MC_all/histograms_diphoton_tree_MC_all.root";
  TFile* fMChists=TFile::Open(histogramMC.Data());
  TH1F* histosmc[12];
  TH1F* histosdata[12];
  TH1F* histosJetJet[12];
  TH1F* histosGammaJet[12];
    
  char canvasname[12];
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
  cout<<lumiNumber<<endl;
  
  
  TLegend *StackLegend[12];
  THStack *StackMC[12];
  for(int i=0;i<12;i++){
   
    c[i] = new TCanvas(canvasname, canvasname, 800, 600);
    c[i]->cd();
    StackLegend[i]= new TLegend(0.70,0.60,0.88,0.75,"","NDC");
    StackLegend[i]->SetFillStyle(0);
    StackLegend[i]->SetTextSize(0.03);
    StackMC[i]= new THStack(fNames[i].Data(),fNames[i].Data());
    cout<<StackMC[i]<<endl;
    histosdata[i]=(TH1F*)fdatahists->Get(fNames[i].Data());
    histosdata[i]->SetMarkerColor(1);
    histosmc[i]= (TH1F*)fMChists->Get(fNames[i].Data());
    cout<<histosmc[i]->GetEntries()<<endl;
    histosmc[i]->SetFillColor(33);
    histosmc[i]->Scale(lumiNumber);
    cout<<"before"<<endl;
    histosJetJet[i]=(TH1F*)fJetJethists->Get(JetJetHists[i].Data());
    cout<<"RAW"<<endl;
    histosGammaJet[i]=(TH1F*)fGammaJethists->Get(GammaJetHists[i].Data());
    cout<<histosJetJet[i]<<endl;
    if (i==0){

      h_TotalBackground->Add(histosJetJet[i],histosGammaJet[i],1,1);
      h_TotalBackground->Add(histosmc[i],1);
      h_DatadivBack->Divide(histosdata[i],h_TotalBackground,1,1); 
      TCanvas* CanvasDatadivBack = new TCanvas("CDatadivBack","CDatadivBack", 800, 600);
      CanvasDatadivBack->cd();
      h_DatadivBack->SetMaximum(3);	
      h_DatadivBack->Draw();
      
	      CanvasDatadivBack->SaveAs("h_DatadivBack_"+Sample+"CompleteOverlay"+".png");
	      
    }   

    c[i]->cd();    
    histosJetJet[i]->SetFillColor(36); 
    histosGammaJet[i]->SetFillColor(38);
    
    
    StackMC[i]->Add(histosJetJet[i]);
    StackMC[i]->Add(histosGammaJet[i]);
    StackMC[i]->Add(histosmc[i]);
    StackLegend[i]->AddEntry(histosdata[0],"Data", "LEP");
    StackLegend[i]->AddEntry(histosmc[0],"SM Diphoton","f");
    StackLegend[i]->AddEntry(histosGammaJet[0],"Photon+Jet","f");
    StackLegend[i]->AddEntry(histosJetJet[0],"Jet+Jet","f");
    if ((strcmp(fNames[i],"h_Diphoton_Minv_high")==0 )|| (strcmp(fNames[i],"h_Photon1_pt")==0 ) || (strcmp(fNames[i],"h_Diphoton_Minv")==0) ||( strcmp(fNames[i],"h_Photon2_pt"))==0){
       c[i]->SetLogy(1);
    }
    
    histosdata[i]->Draw("HIST EP");
    histosdata[i]->SetMinimum(.8*StackMC[i]->GetMinimum());
    histosdata[i]->SetMaximum(1.7*histosdata[i]->GetMaximum());
    StackMC[i]->Draw("HIST SAME") ;
    histosdata[i]->Draw("EPSAME");
    if ((strcmp(fNames[i],"h_Diphoton_cosThetaStar")==0) || (strcmp(fNames[i],"h_Photon1_phi")==0) ||  (strcmp(fNames[i],"h_Photon2_phi")==0)){
      histosdata[i]->SetMinimum(0);
    }
    LumiLabel->Draw();
    StackLegend[i]->Draw();
    gPad->RedrawAxis();
    gPad->Update();
    
    c[i]->SaveAs(fNames[i]+"_"+Sample.Data()+"CompleteOverlay"+".C");
    c[i]->SaveAs(fNames[i]+"_"+Sample.Data()+"CompleteOverlay"+".png");
    c[i]->SaveAs(fNames[i]+"_"+Sample.Data()+"CompleteOverlay"+".pdf");
    if (i==2){
      c[i]->SetLogy(1);
      c[i]->SaveAs(fNames[i]+"_log_"+Sample.Data()+"CompleteOverlay"+".png"); 
 }      
    if (i==0){
      StackCumm = new THStack(fNames[i].Data(),fNames[i].Data());
      CummalitiveCanvas->cd();
      CummalitiveCanvas->SetLogy(1);
      TH1F* h_Diphoton_Minv_log_DataCumm =  MakeChists(histosdata[0]);
      TH1F* h_Diphoton_Minv_log_MCCumm = MakeChists(histosmc[0]); 
      TH1F* h_Diphoton_Minv_log_JetJetCumm = MakeChists(histosJetJet[0]);
      TH1F* h_Diphoton_Minv_log_GammaJetCumm = MakeChists(histosJetJet[0]);
      
      h_Diphoton_Minv_log_MCCumm->SetFillColor(33);
      h_Diphoton_Minv_log_JetJetCumm->SetFillColor(36);
      h_Diphoton_Minv_log_GammaJetCumm->SetFillColor(38);
      
      StackCumm->Add(h_Diphoton_Minv_log_JetJetCumm);
      StackCumm->Add(h_Diphoton_Minv_log_GammaJetCumm);
      StackCumm->Add(h_Diphoton_Minv_log_MCCumm);
      
      
      TLegend* DataMCLegendCumm = new TLegend(0.70,0.6,0.88,0.75,"","NDC");
      DataMCLegendCumm->AddEntry( h_Diphoton_Minv_log_DataCumm,"Data","LEP");
      DataMCLegendCumm->AddEntry(h_Diphoton_Minv_log_MCCumm,"SM Diphoton","F");
      DataMCLegendCumm->AddEntry(h_Diphoton_Minv_log_GammaJetCumm,"Gamma+Jet","F");
      DataMCLegendCumm->AddEntry(h_Diphoton_Minv_log_JetJetCumm,"Jet+Jet","F");
      DataMCLegendCumm->SetFillColor(0);
      h_Diphoton_Minv_log_DataCumm->SetMaximum(2.5*h_Diphoton_Minv_log_DataCumm->GetMaximum());
      h_Diphoton_Minv_log_DataCumm->Draw("EP");
      StackCumm->Draw("HISTSAME");
      h_Diphoton_Minv_log_DataCumm->SetMaximum(2.5*h_Diphoton_Minv_log_DataCumm->GetMaximum());
      h_Diphoton_Minv_log_DataCumm->Draw("EPSAME");
      gPad->RedrawAxis();
      h_Diphoton_Minv_log_DataCumm->SetMaximum(2.5*h_Diphoton_Minv_log_DataCumm->GetMaximum());
      LumiLabelCumm->Draw();
      DataMCLegendCumm->Draw();
      CummalitiveCanvas->SaveAs(fNames[i]+Sample.Data()+"Cumulative"+".png");
      
      dataoverlayedCumulative->cd();
      dataoverlayedCumulative->SetLogy(1);
      h_Diphoton_Minv_log_DataCumm->Draw("EP");
      StackCumm->Draw("HISTSAME");
      h_Diphoton_Minv_log_DataCumm->Draw("EPSAME");
      gPad->RedrawAxis();
      LumiLabelCumm->Draw();
      DataMCLegendCumm->Draw();
      histosdata[i]->SetMarkerColor(kRed);   
      
       histosdata[i]->Draw("EPSAME");
        dataoverlayedCumulative->SaveAs(fNames[i]+Sample.Data()+"DataOverlayedCumulative"+".png");
      
  
      
       
     }
    gPad->Update();

    }

  }


      


 






 
