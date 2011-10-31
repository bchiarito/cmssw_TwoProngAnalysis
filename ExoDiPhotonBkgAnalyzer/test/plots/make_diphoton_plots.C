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
#include <iostream>
#include <math.h>
#include <string.h>
#include <tdrstyle.C>
//#include "fTree.h"

//const int nHists=77+15+1+2+2+3+2+1+1+1+5+1+9+1;
const int nHists=125;
TString nameHists[nHists] = { "h_TrigHLT",  "h_Nvtx", "h_Diphoton_Minv", "h_Diphoton_Minv_log", "h_Diphoton_Minv_high", "h_Diphoton_Minv_low", "h_Diphoton_Minv_low_bin1GeV", "h_Diphoton_qt", "h_Diphoton_deltaPhi", "h_Diphoton_deltaEta", "h_Diphoton_deltaR", "h_Diphoton_cosThetaStar", "h_Photon1_pt", "h_Photon1_pt_log", "h_Photon1_pt_zoom", "h_Photon1_eta", "h_Photon1_phi", "h_Photon1_r9", "h_Photon1_sigmaIetaIeta", "h_Photon1_sigmaEtaEta", "h_Photon1_e2x2e4x4", "h_Photon1_e2e9", "h_Photon1_swisscross", "h_Photon1_severityLevel", "h_Photon1_recHitFlag", "h_Photon1_maxRecHitTime", "h_Photon1_maxRecHitTime_wide", "h_Photon1_e2e9_v_maxRecHitTime", "h_Photon1_hadOverEm", "h_Photon1_hcalIso04", "h_Photon1_hcalIso03", "h_Photon1_ecalIso04", "h_Photon1_ecalIso03", "h_Photon1_trkIsoSumPtHollow04", "h_Photon1_trkIsoSumPtSolid04", "h_Photon1_trkIsoNtrksHollow04", "h_Photon1_trkIsoNtrksSolid04", "h_Photon1_trkIsoSumPtHollow03", "h_Photon1_trkIsoSumPtSolid03", "h_Photon1_trkIsoNtrksHollow03", "h_Photon1_trkIsoNtrksSolid03", "h_Photon2_pt", "h_Photon2_pt_log", "h_Photon2_pt_zoom", "h_Photon2_eta", "h_Photon2_phi", "h_Photon2_r9", "h_Photon2_sigmaIetaIeta", "h_Photon2_sigmaEtaEta", "h_Photon2_e2x2e4x4", "h_Photon2_e2e9", "h_Photon2_swisscross", "h_Photon2_severityLevel", "h_Photon2_recHitFlag", "h_Photon2_maxRecHitTime", "h_Photon2_maxRecHitTime_wide", "h_Photon2_e2e9_v_maxRecHitTime", "h_Photon2_hadOverEm", "h_Photon2_hcalIso04", "h_Photon2_hcalIso03", "h_Photon2_ecalIso04", "h_Photon2_ecalIso03", "h_Photon2_trkIsoSumPtHollow04", "h_Photon2_trkIsoSumPtSolid04", "h_Photon2_trkIsoNtrksHollow04", "h_Photon2_trkIsoNtrksSolid04", "h_Photon2_trkIsoSumPtHollow03", "h_Photon2_trkIsoSumPtSolid03", "h_Photon2_trkIsoNtrksHollow03", "h_Photon2_trkIsoNtrksSolid03","h_FakeRate_tt_pt1","h_FakeRate_tt_eta1","h_FakeRate_tt_pt2","h_FakeRate_tt_eta2","h_FakeRate_tt_minv","h_FakeRate_tt_minv_high","h_FakeRate_tt_qt","h_FakeRate_tt_deltaPhi","h_FakeRate_tt_deltaEta","h_FakeRate_tt_deltaR","h_FakeRate_tt_cosThetaStar","h_FakeRate_tf_pt1","h_FakeRate_tf_eta1","h_FakeRate_tf_pt2","h_FakeRate_tf_eta2","h_FakeRate_tf_minv","h_FakeRate_tf_minv_high","h_FakeRate_tf_qt","h_FakeRate_tf_deltaPhi","h_FakeRate_tf_deltaEta","h_FakeRate_tf_deltaR","h_FakeRate_tf_cosThetaStar","h_FakeRate_ff_pt1","h_FakeRate_ff_eta1","h_FakeRate_ff_pt2","h_FakeRate_ff_eta2","h_FakeRate_ff_minv","h_FakeRate_ff_minv_high","h_FakeRate_ff_qt","h_FakeRate_ff_deltaPhi","h_FakeRate_ff_deltaEta","h_FakeRate_ff_deltaR","h_FakeRate_ff_cosThetaStar","h_Photon1_occupancy","h_Photon2_occupancy","h_Diphoton_Minv_120to200","h_Diphoton_Minv_200to500","h_Diphoton_Minv_500to800", "h_Diphoton_Minv_800toInf","h_Diphoton_Minv_Yousi5", "h_Diphoton_Minv_Yousi10","h_Diphoton_Minv_Yousi40","h_Diphoton_Minv_Yousi100","h_FakeRate_tt_minv_120to200","h_FakeRate_tt_minv_200to500","h_FakeRate_tt_minv_500to800","h_FakeRate_tt_minv_800toInf","h_FakeRate_tf_minv_120to200","h_FakeRate_tf_minv_200to500","h_FakeRate_tf_minv_500to800","h_FakeRate_tf_minv_800toInf","h_FakeRate_ff_minv_120to200","h_FakeRate_ff_minv_200to500","h_FakeRate_ff_minv_500to800","h_FakeRate_ff_minv_800toInf"

};

// prints histos on individual canvases
void draw_individual_histos(TString sample, Bool_t kPrint=kTRUE, Bool_t kStacked=kFALSE, TString categoryEBEE="ALL") {
  
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1,0);
  gStyle->SetOptStat("ourme");
  //  gStyle->SetOptStat(10);
  
  TString outName = TString::Format("histograms_%s_%s.root",sample.Data(),categoryEBEE.Data());

  TFile* fhists = new TFile(outName.Data());

  TCanvas* c[nHists];
  TH1F* histos[nHists];
  char cname[100]; 

  for (int i=0; i<nHists; i++) {
    //    cout << "******" << endl;
    //    cout << i << " " << nameHists[i] << endl;
    sprintf(cname,"c%i",i);
    histos[i] = (TH1F*)fhists->Get(nameHists[i].Data());
    TString histoTitle = histos[i]->GetTitle();
    histoTitle = "(" + sample + ") " + histoTitle;
    histos[i]->SetTitle(histoTitle.Data());

    if (i==0) { 
            continue;
      c[i] = new TCanvas(cname, nameHists[i].Data(), 1000, 60); 
      c[i]->cd();
      c[i]->SetBottomMargin(0.4);
      //      histos[i]->GetXaxis()->SetTitleSize(0.01);
      histos[i]->GetXaxis()->SetBit(TAxis::kLabelsVert);
      if (kStacked) { histos[i]->Draw("hist"); } else { histos[i]->Draw(); }
    } else { 
      c[i] = new TCanvas(cname, nameHists[i].Data(), 600, 600);
      c[i]->cd();
      if (kStacked) { histos[i]->Draw("hist"); } else { histos[i]->Draw(); }
      if ((nameHists[i]=="h_Photon1_e2e9_v_maxRecHitTime")||(nameHists[i]=="h_Photon2_e2e9_v_maxRecHitTime")) { histos[i]->Draw("colz"); c[i]->SetLogz(); }
      if ((nameHists[i]=="h_Photon1_occupancy")||(nameHists[i]=="h_Photon2_occupancy")) { histos[i]->Draw("colz"); }
      if ((nameHists[i]=="h_Photon1_hadOverEm")||(nameHists[i]=="h_Photon1_hcalIso04")||(nameHists[i]=="h_Photon1_hcalIso03")||(nameHists[i]=="h_Photon1_trkIsoSumPtHollow04")||(nameHists[i]=="h_Photon1_trkIsoSumPtSolid04")||(nameHists[i]=="h_Photon1_trkIsoSumPtHollow03")||(nameHists[i]=="h_Photon1_trkIsoSumPtSolid03")||(nameHists[i]=="h_Photon1_maxRecHitTime_wide")||(nameHists[i]=="h_Photon1_e2x2e4x4")||(nameHists[i]=="h_Photon1_e2e9")||(nameHists[i]=="h_Photon2_hadOverEm")||(nameHists[i]=="h_Photon2_hcalIso04")||(nameHists[i]=="h_Photon2_hcalIso03")||(nameHists[i]=="h_Photon2_trkIsoSumPtHollow04")||(nameHists[i]=="h_Photon2_trkIsoSumPtSolid04")||(nameHists[i]=="h_Photon2_trkIsoSumPtHollow03")||(nameHists[i]=="h_Photon2_trkIsoSumPtSolid03")||(nameHists[i]=="h_Photon2_maxRecHitTime_wide")||(nameHists[i]=="h_Photon2_e2x2e4x4")||(nameHists[i]=="h_Photon2_e2e9")||(nameHists[i]=="h_Diphoton_Minv_log")||(nameHists[i]=="h_Photon1_pt_log")||(nameHists[i]=="h_Photon2_pt_log")||(nameHists[i]=="h_Photon1_pt_zoom")||(nameHists[i]=="h_Photon2_pt_zoom")) { 
	//      c[i]->SetLogy();    // this one
      }
      //      if ((nameHists[i]=="h_Photon1_hadOverEm")||(nameHists[i]=="h_Photon1_hcalIso04")||(nameHists[i]=="h_Photon1_hcalIso03")||(nameHists[i]=="h_Photon1_trkIsoSumPtHollow04")||(nameHists[i]=="h_Photon1_trkIsoSumPtSolid04")||(nameHists[i]=="h_Photon1_trkIsoSumPtHollow03")||(nameHists[i]=="h_Photon1_trkIsoSumPtSolid03")||(nameHists[i]=="h_Photon1_maxRecHitTime_wide")||(nameHists[i]=="h_Photon1_e2x2e4x4")||(nameHists[i]=="h_Photon1_e2e9")||(nameHists[i]=="h_Photon2_hadOverEm")||(nameHists[i]=="h_Photon2_hcalIso04")||(nameHists[i]=="h_Photon2_hcalIso03")||(nameHists[i]=="h_Photon2_trkIsoSumPtHollow04")||(nameHists[i]=="h_Photon2_trkIsoSumPtSolid04")||(nameHists[i]=="h_Photon2_trkIsoSumPtHollow03")||(nameHists[i]=="h_Photon2_trkIsoSumPtSolid03")||(nameHists[i]=="h_Photon2_maxRecHitTime_wide")||(nameHists[i]=="h_Photon2_e2x2e4x4")||(nameHists[i]=="h_Photon2_e2e9")) { 
      //	c[i]->SetLogy();   
      //      }
    }
    char fname[100];     
    sprintf(fname,"%s_%s.png",nameHists[i].Data(),sample.Data());  
    if (kPrint) c[i]->Print(fname);
    sprintf(fname,"%s_%s.pdf",nameHists[i].Data(),sample.Data());  
    if (kPrint) c[i]->Print(fname);
    //    sprintf(fname,"%s_%s.C",nameHists[i].Data(),sample.Data());  
    //    if (kPrint) c[i]->Print(fname);
  }

  return;

}


TCanvas* cumulativeBackground(TH1F* histoData =0, THStack* stackMC = 0, int nsamp=4) {

  TList* mcHists = stackMC->GetHists();
  //  mcHists->Print();
  TIter next (mcHists);

  const int nsamples=nsamp;
  TH1F* h[nsamples];
  for (int i=0; i<nsamples;i++) {
    h[i] = (TH1F*)next();
  }  
  //  h[0] = (TH1F*)next();
  //  h[1] = (TH1F*)next();
  //  h[2] = (TH1F*)next();
  //  h[3] = (TH1F*)next();
    
  float countData = 0;
  float countMC[nsamples];
  for (int i=0; i<nsamples;i++) {
    countMC[i]=0;
  }

  TH1F* hCumulData = (TH1F*)histoData->Clone("hCumulData") ;
  TH1F* hCumulMC = (TH1F*)h[0]->Clone("hCumulMC") ;

  int nbins =  h[0]->GetNbinsX(); 
  float binwidth = h[0]->GetBinWidth(1);

  if (hCumulData->GetNbinsX() != hCumulMC->GetNbinsX()) {
    cout << "Data and MC histos have different bin sizes! " << hCumulData->GetNbinsX() << " " << hCumulMC->GetNbinsX() << endl;
    return 0;
  }
  const int ncumulbins = nbins;
  float countCumulData[ncumulbins];
  float countCumulMC[ncumulbins];
  for (int i=0;i<ncumulbins; i++) {
    countCumulData[i]=0;
    countCumulMC[i]=0;
  }

  for (int i=nbins; i>0; i--) {
    //  for (int i=nbins; i>0; i--) {
    if (i==nbins) {
      countData += histoData->GetBinContent(i);
      countData += histoData->GetBinContent(i+1);
    } else {
      countData += histoData->GetBinContent(i);
    }
    countCumulData[i-1] += countData; 
    for (int j=0; j<nsamples;j++) {
      countMC[j] += h[j]->GetBinContent(i);
      countCumulMC[i-1] += countMC[j];
      //      cout << i << " " << j << " " << 100+(i-1)*binwidth << " " << h[j]->GetBinContent(i) << " " << countMC[j] << " " << countCumulMC[i-1]  << endl;
    }
    if (i==21) cout << "Cumulative Background Above " << 100+(i-1)*binwidth << " : " << countCumulMC[i-1]  << " (Expected) " << countCumulData[i-1] << " (Actual) " << endl;
    if (i==36) cout << "Cumulative Background Above " << 100+(i-1)*binwidth << " : " << countCumulMC[i-1]  << " (Expected) " << countCumulData[i-1] << " (Actual) " << endl;


  }

  for(int i=1; i<nbins+1; i++) {
    hCumulData->SetBinContent(i,countCumulData[i-1]);
    hCumulMC->SetBinContent(i,countCumulMC[i-1]);
  }
  float totalCount = 0;
  for (int j=0; j<nsamples;j++) {
    totalCount+= countMC[j];
  }
  //  cout << "Total Count " << totalCount << endl;
  //  cout << countCumulMC[0] << " " << countCumulData[0] << endl;

  cout << hCumulData->GetMaximum() << " " << hCumulMC->GetMaximum() << endl;
  hCumulData->SetMaximum(1.05*TMath::Max(hCumulData->GetMaximum(),hCumulMC->GetMaximum()));
  cout << hCumulData->GetMaximum() << " " << hCumulMC->GetMaximum() << endl;

  TCanvas* cCumul = new TCanvas("cCumul","Cumulative Background", 600, 600);
  cCumul->cd();
  hCumulMC->SetFillColor(38);
  hCumulData->Draw("e");
  hCumulMC->Draw("samehist");
  hCumulData->Draw("esame");

  //  f.cd();
  //  cCumul->Write();
  return cCumul;

}


// prints fake rate histos
void fake_rate_histos3(TString sample, Bool_t kPrint=kTRUE, Float_t lumi=0.0, TString categoryEBEE="ALL") {

  cout << " FAKE RATE METHOD " << endl;

  gROOT->SetStyle("Plain");
  //  gStyle->SetOptStat("oue");
  gStyle->SetOptStat(000000);
  gStyle->SetMarkerStyle(20);
  
  TString outName = TString::Format("histograms_%s_TF_%s.root",sample.Data(),categoryEBEE.Data());
  TFile* fileTF = new TFile(outName.Data());
  
  outName = TString::Format("histograms_%s_FT_%s.root",sample.Data(),categoryEBEE.Data());
  TFile* fileFT = new TFile(outName.Data());
  
  outName = TString::Format("histograms_%s_FF_%s.root",sample.Data(),categoryEBEE.Data());
  TFile* fileFF = new TFile(outName.Data());

  outName = TString::Format("histograms_%s_%s.root",sample.Data(),categoryEBEE.Data());
  TFile* fhists = new TFile(outName.Data(),"UPDATE");
  fhists->cd();

  // color scheme
  //  Color_t kRealDiphotonMCColor = kBlue;//-9;
  //  Color_t kTightFakeColor = kGreen-6;
  //  Color_t kFakeFakeColor = kMagenta-7;
  int kRealDiphotonMCColor = TColor::GetColor("#66ccff");
  int kTightFakeColor = TColor::GetColor("#3399ff");
  int kFakeFakeColor = TColor::GetColor("#0066cc");
  
  outName = TString::Format("histograms_DiPhoton_all_%s.root",categoryEBEE.Data());
  TFile* fmc = new TFile(outName.Data());

  const int nFakeHists = 13+3+1;

  TString diphotonHists[nFakeHists] = {"h_Diphoton_Minv","h_Diphoton_Minv_high","h_Diphoton_qt","h_Diphoton_deltaPhi","h_Diphoton_deltaEta","h_Diphoton_deltaR","h_Diphoton_cosThetaStar","h_Photon1_pt","h_Photon1_eta","h_Photon1_phi","h_Photon2_pt","h_Photon2_eta","h_Photon2_phi","h_Diphoton_Minv_120to200","h_Diphoton_Minv_200to500","h_Diphoton_Minv_500to800","h_Diphoton_Minv_800toInf"};

  TString ttHists[nFakeHists] = {"h_FakeRate_tt_minv","h_FakeRate_tt_minv_high","h_FakeRate_tt_qt","h_FakeRate_tt_deltaPhi","h_FakeRate_tt_deltaEta","h_FakeRate_tt_deltaR","h_FakeRate_tt_cosThetaStar","h_FakeRate_tt_pt1","h_FakeRate_tt_eta1", "h_FakeRate_tt_phi1","h_FakeRate_tt_pt2","h_FakeRate_tt_eta2","h_FakeRate_tt_phi2","h_FakeRate_tt_minv_120to200","h_FakeRate_tt_minv_200to500","h_FakeRate_tt_minv_500to800","h_FakeRate_tt_minv_800toInf" };

  TString tfHists[nFakeHists] = {"h_FakeRate_tf_minv","h_FakeRate_tf_minv_high","h_FakeRate_tf_qt","h_FakeRate_tf_deltaPhi","h_FakeRate_tf_deltaEta","h_FakeRate_tf_deltaR","h_FakeRate_tf_cosThetaStar","h_FakeRate_tf_pt1","h_FakeRate_tf_eta1", "h_FakeRate_tf_phi1","h_FakeRate_tf_pt2","h_FakeRate_tf_eta2","h_FakeRate_tf_phi2","h_FakeRate_tf_minv_120to200","h_FakeRate_tf_minv_200to500","h_FakeRate_tf_minv_500to800","h_FakeRate_tf_minv_800toInf"};

  TString ffHists[nFakeHists] = {"h_FakeRate_ff_minv","h_FakeRate_ff_minv_high","h_FakeRate_ff_qt","h_FakeRate_ff_deltaPhi","h_FakeRate_ff_deltaEta","h_FakeRate_ff_deltaR","h_FakeRate_ff_cosThetaStar","h_FakeRate_ff_pt1","h_FakeRate_ff_eta1", "h_FakeRate_ff_phi1","h_FakeRate_ff_pt2","h_FakeRate_ff_eta2","h_FakeRate_ff_phi2","h_FakeRate_ff_minv_120to200","h_FakeRate_ff_minv_200to500","h_FakeRate_ff_minv_500to800","h_FakeRate_ff_minv_800toInf"};
  
  TString fNames[nFakeHists] = { "minv", "minv_high", "qt", "deltaPhi", "deltaEta", "deltaR", "cosThetaStar", "pt1", "eta1", "phi1", "pt2", "eta2", "phi2", "minv_120to200", "minv_200to500","minv500to800","minv800toInf" };

  Double_t kFactor = 1.0;

  setTDRStyle();
 
 //  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  TCanvas *c1 = new TCanvas("c1", "",4,30,600,600);
  gStyle->SetOptStat(0);
  c1->Range(-52.28915,-4.771183,1273.012,4.098671);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  //  c1->SetLogy();
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetLeftMargin(0.13);
  c1->SetRightMargin(0.04);
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.13);
  c1->SetFrameLineWidth(2);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderSize(2);
  c1->SetFrameLineWidth(2);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderSize(2);

   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetLogy();
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.13);
   c1->SetRightMargin(0.04);
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.13);
   c1->SetFrameLineWidth(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderSize(2);
   c1->SetFrameLineWidth(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderSize(2);

       
  for (int i=0; i<nFakeHists; i++) {
    //    cout << i << " " << diphotonHists[i] << endl;
    TH1F* h_Diphoton = (TH1F*) fmc->Get(diphotonHists[i].Data());
    h_Diphoton->SetFillColor(kRealDiphotonMCColor);

//     if (i>1) {

//       // flat scale factor for all but minv plots
//       h_Diphoton->Scale(kFactor*lumi);

//     } else {
      
//       // mgg dependent k-factor for 
        
//       TF1* fKfactor = new TF1("kfactor", "-6.08061+9.90359/x^(3.88349e-02)", 100, 2000) ;

//       int nbins      = h_Diphoton->GetNbinsX();
//       float binwidth = h_Diphoton->GetBinWidth(1);
//       float start    = h_Diphoton->GetBinLowEdge(1);
      
//       for (int j=1; j<nbins+1; j++) {
// 	float mass =  start+(j-1)*binwidth;
// 	//	cout << j << " " << mass << " " << h_Diphoton->GetBinContent(j) << " " <<  fKfactor->Eval(mass) << " " << h_Diphoton->GetBinContent(j)*fKfactor->Eval(mass) << endl;
	
// 	h_Diphoton->SetBinContent(j,h_Diphoton->GetBinContent(j)*fKfactor->Eval(mass));
//       }
//       h_Diphoton->Scale(lumi);          

//     }

    h_Diphoton->Scale(lumi);          
    
    //    cout << "Scaling Fake Rate Estimates To: " << lumi << endl;
    
    TH1F *h_FakeRate_tt = (TH1F*) fhists->Get(ttHists[i].Data()); 
    TH1F *h_FakeRate_tf = (TH1F*) fileTF->Get(tfHists[i].Data());  
    TH1F *h_FakeRate_ft = (TH1F*) fileFT->Get(tfHists[i].Data());  
    TH1F *h_FakeRate_ff = (TH1F*) fileFF->Get(ffHists[i].Data());
    
    //     cout << h_FakeRate_tt->GetEntries() << " " << h_FakeRate_tf->GetEntries() + h_FakeRate_ft->GetEntries() << " " << h_FakeRate_ff->GetEntries() << endl;
    
    // add TF and FT here
    
    //    h_FakeRate_tf->Add(h_FakeRate_ft);
    
    h_FakeRate_tf->SetFillColor(kTightFakeColor);	
    h_FakeRate_ft->SetFillColor(kTightFakeColor);	
    h_FakeRate_ff->SetFillColor(kFakeFakeColor);
    
    
    TH1F* h_FakeRate_gJet = (TH1F*)h_FakeRate_tf->Clone();
    h_FakeRate_gJet->Add(h_FakeRate_ft);
    h_FakeRate_gJet->Add(h_FakeRate_ff,-2);
    

    TH1F* hErr = (TH1F*)h_FakeRate_gJet->Clone("hErr");
    hErr->Add(h_FakeRate_ff);
    hErr->Add(h_Diphoton);
    
    // error
    //     TH1F* hErr = (TH1F*)h_FakeRate_gJet->Clone("hErr") ;
    //     hErr->Add(h_FakeRate_ff);
    //     hErr->Add(h_Diphoton);
    
    TH1F* hAll = (TH1F*)hErr->Clone("hAll");
    
    for (int iBin = 1 ; iBin <= hErr->GetNbinsX() ; ++iBin) {
      
      //       cout << iBin << endl;
      float x = h_FakeRate_tf->GetBinContent(iBin) +  h_FakeRate_ft->GetBinContent(iBin);
      float y = h_FakeRate_ff->GetBinContent(iBin) ;
      float z = h_Diphoton->GetBinContent(iBin) ;
      //      float gamjet = x - 2*y;

      // error on total 
      float value = hErr->GetBinContent(iBin) ;
      //      float err = sqrt( (0.2*x)*(0.2*x) + (0.2*sqrt(2)*y)*(0.2*sqrt(2)*y) + (0.13*z)*(0.13*z));
      //      float err = sqrt( ((x - 2*y)*0.2)*((x - 2*y)*0.2) + (0.13*z)*(0.13*z) ); //0.2*gamjet;
      float err = sqrt( ((x - 2*y)*0.2)*((x - 2*y)*0.2) + (0.18*z)*(0.18*z) ); //0.2*gamjet;
      float errStat = 0; //hErr->GetBinError(iBin) ;
      float errTotal = 0;
      if (value!=0) errTotal = value*TMath::Sqrt((err/value)*(err/value) + (errStat/value)*(errStat/value)); 
      hErr->SetBinError(iBin, errTotal) ;
      
      //       cout << "Total " << value << " +- " << errStat << " +- " << err << endl;
      
      // error on FF
      value = h_FakeRate_ff->GetBinContent(iBin) ;
      //      err = 0.2*sqrt(2)*y;
      err = 0.2*2*y;
      errStat = 0; //h_FakeRate_ff->GetBinError(iBin) ;
      float errTotal = 0;
      if (value!=0) errTotal = value*TMath::Sqrt((err/value)*(err/value) + (errStat/value)*(errStat/value)); 
      h_FakeRate_ff->SetBinError(iBin, errTotal) ;              
      
      //       cout << "QCD   " << value << " +- " << errStat << " +- " << err << endl;
      
      // error on gammaJet
      value = h_FakeRate_gJet->GetBinContent(iBin) ;
      //      err = sqrt( (0.2*x)*(0.2*x) + 4*(0.2*sqrt(2)*y)*(0.2*sqrt(2)*y) );
      err = (x - 4*y)*0.2; //(gamjet - 2*y)*0.2;
      errStat = 0; //h_FakeRate_gJet->GetBinError(iBin) ; 
      float errTotal = 0;
      if (value!=0) errTotal = value*TMath::Sqrt((err/value)*(err/value) + (errStat/value)*(errStat/value)); 
      h_FakeRate_gJet->SetBinError(iBin, errTotal) ;              
      
      //       cout << "GammaJet " << value << " +- " << errStat << " +- " << err << endl;
      
      // error on Diphoton
      value = h_Diphoton->GetBinContent(iBin) ;
      //      err = sqrt( (0.13*z)*(0.13*z) );
      err = sqrt( (0.18*z)*(0.18*z) );
      errStat = 0; //h_Diphoton->GetBinError(iBin) ;
      float errTotal = 0;
      if (value!=0) errTotal = value*TMath::Sqrt((err/value)*(err/value) + (errStat/value)*(errStat/value)); 
      h_Diphoton->SetBinError(iBin, errTotal) ;                     
      
      //       cout << "Diphoton " << value << " +- " << errStat << " +- " << err << endl;
      
    }
    
    //      if (i<1) {
    
    //        // compute yields and errors
    
//        const int nmbins = 3; // 120to200,200to500,500to1000
//        TString namebins[nmbins] = {"120to200", "200to500","500to1000"}

//        //       const int nsamp = 5;  //  QCD, GammaJet, Diphoton, Total, Data
//        const int nsamp = 1;  //  QCD, GammaJet, Diphoton, Total, Data
//        float num[nmbins][nsamp];
//        float error[nmbins][nsamp];
//        float errorSyst[nmbins][nsamp];
       
//        for (int i=0; i<nmbins; i++) {
// 	 for (int j=0; j<nsamp; j++) {
// 	   num[i][j]=0;
// 	   error[i][j]=0;
// 	 }
//        }
       
//        const int nbins = h_Diphoton_Minv->GetNbinsX();       
       
//        for (int i=1; i<nbins+1; i++) {
	 
// 	 double N[nsamp];	
// 	 N[0] = h_FakeRate_ff->GetBinContent(i); // qcd
// 	 //	 N[1] = h_FakeRate_gJet->GetBinContent(i);      // gammajet
// 	 //	 N[2] = h_Diphoton->GetBinContent(i);    //diphoton
// 	 //	 N[3] = h_FakeRate_ff->GetBinContent(i)+h_FakeRate_gJet->GetBinContent(i)+h_Diphoton->GetBinContent(i);
// 	 //	 N[4] = h_FakeRate_tt->GetBinContent(i); //data
	 
// 	 double eStat[nsamp];	
// 	 eStat[0] = h_FakeRate_ff->GetBinError(i);    
// 	 //	 eStat[1] = h_FakeRate_gJet->GetBinError(i);    
// 	 //	 eStat[2] = h_Diphoton->GetBinError(i);
// 	 //	 eStat[3] = hAll->GetBinError(i);    
// 	 //	 eStat[4] = h_FakeRate_tt->GetBinError(i);

// 	 cout << "stat " << eStat[0]/N[0] << endl;
// 	 //	 cout << "stat " << eStat[0]/N[0] << " " << eStat[1]/N[1] << " " << eStat[2]/N[2] << " " << eStat[3]/N[3] << " " << eStat[4]/N[4] << endl;

// 	 //	 cout << i << " " << 120 + (i-1)*20 << " " << x << " " << y << " " << z << endl;

// 	 //	 cout << i << " " << 120 + (i-1)*20  << " : " << N[0] << " +- " << e[0] << " : " << N[1] << " +- " << e[1] << endl;

// 	 if (i<5) {
// 	   for (int j=0; j<nsamp; j++)  {
// 	     num[0][j] += N[j];
// 	     if (N[j]!=0) {
// 	       cout << " before " << j << " " << N[j] << " " << (eStat[j]/N[j]) << " " <<  eStat[j] << endl;
// 	       //	       error[0][j] += (eStat[j]/N[j])*(eStat[j]/N[j]);	      
// 	       error[0][j] += (eStat[j])*(eStat[j]);	      
// 	     }
// 	   }	   

// 	   cout << " blah " << i << " " << N[0] << " " << N[0]*sqrt(error[0][0]) << endl;

// 	 } else if (i>=5 && i<20) {
// 	   for (int j=0; j<nsamp; j++)  {
// 	     num[1][j] += N[j];
// 	     if (N[j]!=0) {
// 	       //	       error[1][j] += (eStat[j]/N[j])*(eStat[j]/N[j]);
// 	       error[1][j] += (eStat[j])*(eStat[j]);
// 	     }
// 	   }
// 	 } else {
// 	   for (int j=0; j<nsamp; j++)  {
// 	     num[2][j] += N[j];
// 	     if (N[j]!=0) {
// 	       //	       error[2][j] += (eStat[j]/N[j])*(eStat[j]/N[j]);
// 	       error[2][j] += (eStat[j])*(eStat[j]);
// 	     }
// 	   }
// 	 }		 
	 
//        }

//        cout << " errors before systematic " << endl;
//        for (int i=0; i<nmbins; i++) {
// 	 //	 cout << namebins[i].Data() << " " << num[i][0]*sqrt(error[i][0]) << " " << num[i][1]*sqrt(error[i][1]) << " " << num[i][2]*sqrt(error[i][2]) << " " << num[i][3]*sqrt(error[i][3]) << " " << num[i][4]*sqrt(error[i][4]) << endl;
// 	 //	 cout << namebins[i].Data() << " " << num[i][0]*sqrt(error[i][0]) << endl;
// 	 cout << namebins[i].Data() << " " << sqrt(error[i][0]) << endl;
//        }

//        return;

//        cout << " now for systematic " << endl;              
//        for (int i=0; i<nmbins; i++)  {

// 	 float x = num[i][1]-2*num[i][0]; //h_FakeRate_tf->GetBinContent(i) +  h_FakeRate_ft->GetBinContent(i);
// 	 float y = num[i][0];             //h_FakeRate_ff->GetBinContent(i) ;
// 	 float z = num[i][2];             //h_Diphoton->GetBinContent(i) ;
	 
// 	 cout << "xyz " <<x << " " << y << " " << z << endl;
// 	 double eSyst[nsamp];	
// 	 eSyst[0] = 0.2*sqrt(2)*y;
// 	 eSyst[1] = sqrt( (0.2*x)*(0.2*x) + (4*0.2*sqrt(2)*y)*(4*0.2*sqrt(2)*y) );
// 	 eSyst[2] = sqrt( (0.08*z)*(0.08*z) );
// 	 eSyst[3] = sqrt( (0.2*x)*(0.2*x) + (0.2*sqrt(2)*y)*(0.2*sqrt(2)*y) + (0.08*z)*(0.08*z));
// 	 eSyst[4] = 0;

// 	 cout << "syst " << eSyst[0] << " " << eSyst[1] << " " << eSyst[2] << " " << eSyst[3] << " " << eSyst[4] << endl;
	 	 	 
// 	 errorSyst[i][0] += (eSyst[0]/num[i][0])*(eSyst[0]/num[i][0]);
// 	 errorSyst[i][1] += (eSyst[1]/num[i][1])*(eSyst[1]/num[i][1]);
// 	 errorSyst[i][2] += (eSyst[2]/num[i][2])*(eSyst[2]/num[i][2]);
// 	 errorSyst[i][3] += (eSyst[3]/num[i][3])*(eSyst[3]/num[i][3]);
// 	 errorSyst[i][4] += (eSyst[4]/num[i][4])*(eSyst[4]/num[i][4]);

//        }
      

//        for (int i=0; i<nmbins; i++)  {

// 	 for (int j=0; j<nsamp; j++)  {
// 	   error[i][j] = num[i][j]*sqrt(error[i][j]);
// 	   errorSyst[i][j] = num[i][j]*sqrt(errorSyst[i][j]);
// 	 }

// 	 cout << "total " << error[i][0] << " " << error[i][1] << " " << error[i][2] << " " << error[i][3] << " " << error[i][4] << endl;
       
//        }

//        for (int i=0; i<nmbins; i++)  {

// 	 printf(">>>> BACKGROUND ESTIMATE FOR %s\n",namebins[i].Data());
// 	 printf("Diphoton \t%0.1f\t $\pm$ %0.1f\t $\pm$ %0.1f\n", num[i][2],error[i][2],errorSyst[i][2]);
// 	 printf("GammaJet \t%0.1f\t $\pm$ %0.1f\t $\pm$ %0.1f\n",num[i][1],error[i][1],errorSyst[i][1]);
// 	 printf("DiJet    \t%0.1f\t $\pm$ %0.1f\t $\pm$ %0.1f\n",num[i][0],error[i][0],errorSyst[i][0]);
// 	 printf("Total    \t%0.1f\t $\pm$ %0.1f\t $\pm$ %0.1f\n",num[i][3],error[i][3],errorSyst[i][3]);
// 	 printf("Data     \t%i\n",num[i][4],error[i][4]);
	
	 
//        }
       
       
// //        cout << ">>>> BACKGROUND ESTIMATE FOR " << fNames[i] << endl;
// //        cout << "Diphoton \t" << h_Diphoton->GetBinContent(1)      << " +- \t" << h_Diphoton->GetBinError(1) << endl;
// //        cout << "GammaJet \t" << h_FakeRate_gJet->GetBinContent(1) << " +- \t" << h_FakeRate_gJet->GetBinError(1) << endl;
// //        cout << "DiJet \t"    << h_FakeRate_ff->GetBinContent(1)   << " +- \t" << h_FakeRate_ff->GetBinError(1) << endl; 
// //        cout << "Total \t"    << hErr->GetBinContent(1)            << " +- \t" << hErr->GetBinError(1) << " " << endl;
// //        cout << "Data  \t"    << h_FakeRate_tt->GetBinContent(1)   << endl;      


//      }
     
     char fname[100];     
     sprintf(fname,"h_FakeRate_%s_stack",fNames[i].Data());

     THStack *h_FakeRate_stack = new THStack(fname,"");
  
     h_FakeRate_stack->Add(h_FakeRate_ff);
     h_FakeRate_stack->Add(h_FakeRate_gJet); 
     h_FakeRate_stack->Add(h_Diphoton);

//      // printout
//      if (i==0) {
//        printf(" Diphoton  Estimate %0.1f\n",h_Diphoton->Integral());
//        printf(" GammaJet  Estimate %0.1f\n",h_FakeRate_gJet->Integral());
//        printf(" DiJet     Estimate %0.1f\n",h_FakeRate_ff->Integral());
//        printf(" Total Bkg Estimate %0.1f\n",h_Diphoton->Integral()+h_FakeRate_gJet->Integral()+h_FakeRate_ff->Integral());
//        printf(" Total Data         %0.1f\n",h_FakeRate_tt->Integral());
//        //      cout << " Diphoton Estimate " << h_Diphoton->Integral() << endl;
//        //      cout << " GammaJet  Estimate " << h_FakeRate_gJet->Integral() << endl;
//        //      cout << " DiJet  Estimate " << h_FakeRate_ff->Integral() << endl;
//        //      cout << " Total Bkg   Estimate " << h_Diphoton->Integral() + h_FakeRate_gJet->Integral() + h_FakeRate_ff->Integral()  << endl;
//        //      cout << " Total Data           " << h_FakeRate_tt->Integral() << endl;
//      }
     

      if (fNames[i]=="minv_120to200"||fNames[i]=="minv_200to500"||fNames[i]=="minv500to800"||fNames[i]=="minv800toInf") {

	cout << ">>>> BACKGROUND ESTIMATE FOR " << fNames[i] << endl;
	cout << "Diphoton \t" << h_Diphoton->GetBinContent(1)      << " +- \t" << h_Diphoton->GetBinError(1) << endl;
	cout << "GammaJet \t" << h_FakeRate_gJet->GetBinContent(1) << " +- \t" << h_FakeRate_gJet->GetBinError(1) << endl;
	cout << "DiJet \t"    << h_FakeRate_ff->GetBinContent(1)   << " +- \t" << h_FakeRate_ff->GetBinError(1) << endl; 
	cout << "Total \t"    << hErr->GetBinContent(1)            << " +- \t" << hErr->GetBinError(1) << " " << endl;
	cout << "Data  \t"    << h_FakeRate_tt->GetBinContent(1)   << endl;
     
	printf(">>>> BACKGROUND ESTIMATE FOR %s\n",fNames[i].Data());
	printf("Diphoton &\t%0.1f $~\pm $ \t%0.1f\n", h_Diphoton->GetBinContent(1),h_Diphoton->GetBinError(1));
	printf("$\gamma$Jet &\t%0.1f $~\pm $ \t%0.1f\n",h_FakeRate_gJet->GetBinContent(1),h_FakeRate_gJet->GetBinError(1));
	printf("DiJet    &\t%0.1f $~\pm $ \t%0.1f\n",h_FakeRate_ff->GetBinContent(1) , h_FakeRate_ff->GetBinError(1)); 
	printf("Total    &\t%0.1f $~\pm $ \t%0.1f\n",hErr->GetBinContent(1)          , hErr->GetBinError(1));
	printf("Data     &\t%i\n",h_FakeRate_tt->GetBinContent(1)  );
	

	//	cout << " TF Total " << h_FakeRate_tf->GetBinContent(1) +  h_FakeRate_ft->GetBinContent(1) << endl;
	//	cout << " FF Total " << h_FakeRate_ff->GetBinContent(1)  << endl;

     
      }


     // legend
      TLegend *leg = new TLegend(0.6761745,0.701049,0.8020134,0.9125874,NULL,"brNDC");
   //      TLegend *leg = new TLegend(0.6426174,0.5804196,0.7818792,0.791958,NULL,"brNDC");
      //     TLegend *leg = new TLegend(0.7332215,0.5804196,0.8724832,0.791958,NULL,"brNDC");
     //     TLegend *leg = new TLegend(0.5234899,0.5961538,0.6627517,0.8076923,NULL,"brNDC");
     leg->SetBorderSize(1);
     leg->SetTextFont(62);
     leg->SetTextSize(0.03225806);
     leg->SetLineColor(0);
     leg->SetFillColor(0);
     leg->AddEntry(h_FakeRate_tt,"Observed");
     leg->AddEntry(h_Diphoton,"Diphoton");
     leg->AddEntry(h_FakeRate_gJet,"#gamma+jet");
     leg->AddEntry(h_FakeRate_ff,"Dijet");
     leg->AddEntry(hErr,"Syst. Uncertainty");

     // now draw
     c1->cd();
     c1->SetLogy(0);
     //     h_FakeRate_stack->SetMaximum(1.1*TMath::Max(h_FakeRate_stack->GetMaximum(),h_FakeRate_tt->GetMaximum()));
     h_FakeRate_tt->SetMaximum(1.3*TMath::Max(h_FakeRate_stack->GetMaximum(),h_FakeRate_tt->GetMaximum()));
     //     h_FakeRate_tt->SetMaximum(70);
     h_FakeRate_tt->SetMinimum(0.001);

     if ( fNames[i]=="phi1" || fNames[i]=="phi2" ) h_FakeRate_tt->SetMinimum(0);

     h_FakeRate_tt->SetTitle("");
     h_FakeRate_tt->Draw();
     h_FakeRate_tt->GetXaxis()->SetTitle(h_Diphoton->GetXaxis()->GetTitle());
     h_FakeRate_tt->GetYaxis()->SetTitle("Entries/20 GeV/c^{2}");
     h_FakeRate_tt->GetYaxis()->SetTitleOffset(1.5O);
     h_FakeRate_stack->Draw("histsame");
     h_FakeRate_tt->SetMarkerStyle(20);
     h_FakeRate_tt->SetMarkerSize(1.5);

     //     hErr->SetLineColor(0) ;
     //     hErr->SetFillColor(2) ;
     //     hErr->SetFillStyle(3244) ;

     ci = TColor::GetColor("#9900cc");
     hErr->SetFillColor(ci);
     hErr->SetFillStyle(3001);

     // hErr->Draw("e") ;
     //     hErr->SetMarkerColor(kRed);
     //     hErr->SetMarkerStyle(24);
     //     hErr->SetLineColor(kRed);
     hErr->Draw("same E2") ;
     h_FakeRate_tt->Draw("epsame");
    
     if ( fNames[i]!="deltaEta" && fNames[i]!="eta1" && fNames[i]!="eta2" && fNames[i]!="phi1" || fNames[i]!="phi2") leg->Draw();

     if ( fNames[i]=="minv") {
       //       TPaveText *pt = new TPaveText(0.2718121,0.7828467,0.4916107,0.8905109,"blNDC");
       TPaveText *pt = new TPaveText(0.4127517,0.8583916,0.5738255,0.9353147,"blNDC");
       //       TPaveText *pt = new TPaveText(0.3053691,0.8601399,0.466443,0.9370629,"blNDC");
   //       TPaveText *pt = new TPaveText(0.3271812,0.8146853,0.5469799,0.9230769,"blNDC");

       pt->SetName("CMS Preliminary");
       pt->SetBorderSize(1);
       pt->SetFillColor(0);
       pt->SetLineColor(0);
       pt->SetTextSize(0.0454545);
       TText *text = pt->AddText("CMS Preliminary");
       pt->Draw();
    
       //       pt = new TPaveText(0.6442953,0.7828467,0.8842282,0.8905109,"brNDC");
       pt = new TPaveText(0.3741611,0.7814685,0.614094,0.8618881,"brNDC");
   //       pt = new TPaveText(0.6442953,0.8094406,0.8842282,0.9178322,"brNDC");

       //      pt = new TPaveText(0.2600671,0.6788321,0.5,0.7992701,"brNDC");
       char pname[100];     
       sprintf(pname,"%0.0f pb^{-1} at 7 TeV",lumi);      
       pt->SetName(pname);
       pt->SetFillColor(0);
       pt->SetBorderSize(1);
       pt->SetLineColor(0);
       pt->SetTextSize(0.0454545);
       text = pt->AddText(pname);
       pt->Draw();
    
     }

     sprintf(fname,"h_FakeRate_%s_%s.png",fNames[i].Data(),sample.Data());  
     if (kPrint) c1->Print(fname);
     sprintf(fname,"h_FakeRate_%s_%s.pdf",fNames[i].Data(),sample.Data());  
     if (kPrint) c1->Print(fname);
     sprintf(fname,"h_FakeRate_%s_%s.C",fNames[i].Data(),sample.Data());  
     if (kPrint) c1->Print(fname);    
  
     if ( fNames[i]=="minv" || fNames[i]=="pt1" || fNames[i]=="pt2") {

       c1->SetLogy(1); c1->Update();
       sprintf(fname,"h_FakeRate_%s_log_%s.png",fNames[i].Data(),sample.Data());  
       if (kPrint) c1->Print(fname);
       sprintf(fname,"h_FakeRate_%s_log_%s.pdf",fNames[i].Data(),sample.Data());  
       if (kPrint) c1->Print(fname);
       //      sprintf(fname,"h_FakeRate_%s_log_%s.C",fNames[i].Data(),sample.Data());  
       //      if (kPrint) c1->Print(fname);    
       c1->SetLogy(0);

     }

     fhists->cd();

     if ( fNames[i]=="minv" ) {  
       //       cout << " WRITING HERR " << endl;
       hErr->Write(); 
     }

     h_FakeRate_stack->Write();

     // Cumulative Background
     if (i==0) {
       TCanvas* cCumul = cumulativeBackground(h_FakeRate_tt,h_FakeRate_stack,3);
       cCumul->SetLogy();
       cCumul->Write();
       char fname[100];     
       sprintf(fname,"%s_%s.png","h_CumulativeBackground",sample.Data());  
       cCumul->Print(fname);
       sprintf(fname,"%s_%s.pdf","h_CumulativeBackground",sample.Data());  
       cCumul->Print(fname);
     }

  }

  fhists->ls();
  fhists->Close();

  return;

}



void make_diphoton_plots(TString sample = "PhotonJet_Pt30to50", Bool_t kPrint=kTRUE, Bool_t isData=kFALSE, Float_t lumi=0.0, TString categoryEBEE="ALL", Bool_t isSignal=kFALSE, TString version="MC_36X_V2") {

  TString infile;
  if (isData) {
    //    TString infile = "/Users/toyoko/Work/cms/physics/diphoton/2011/ntuples/DATA_42X_V1/diphoton_tree_"+sample+".root";    
    TString infile = "diphoton_tree_"+sample+".root";
  } else {
    TString infile = "/Users/toyoko/Work/cms/physics/diphoton/2011/ntuples/MC_42X_V1/diphoton_tree_"+sample+".root";
  }
  cout << "Making quick check plots for: " << infile << endl;

  gSystem->Load("fTree_C.so");

  // input file
  TFile* f = TFile::Open(infile.Data());  

  TString tempName;

  TTree* fChain;
  fTree* reader;

  if (isData) { 
          
    tempName = "diphotonAnalyzer/fTree"; 
    fChain = (TTree*)f->Get(tempName.Data()); 
    cout << "All entries = " << fChain->GetEntries() <<endl;

  } else if (isSignal) { 

    tempName = "diphotonAnalyzer/fTree"; 
    fChain = (TTree*)f->Get(tempName.Data()); 

  } else { 

    tempName = "diphotonAnalyzer/fTree"; 
    fChain = (TTree*)f->Get(tempName.Data()); 

  }
    
  // set up output
  if (isData) sample = "data_"+sample;
  cout << sample << endl;
  TString outName = TString::Format("histograms_%s_%s.root",sample.Data(),categoryEBEE.Data());
  cout << outName << endl;
  TFile* outF = new TFile(outName.Data(),"RECREATE");

  // loop
  reader = new fTree(fChain);
  reader->_fakeStatus = "TightTight";

  reader->_reweightPU=kTRUE;
  reader->_puHistMC = (TH1F*)f->Get("diphotonAnalyzer/histPUMC");
  reader->_cutPhoton1Pt = 70;
  reader->_cutPhoton2Pt = 70;
  //  reader->_cutEta = 2.5;

  if (isData) reader->_reweightPU=kFALSE;

  if (categoryEBEE!="ALL" && categoryEBEE!="EBEB" && categoryEBEE!="EEEE"  && categoryEBEE!="EBEE" && categoryEBEE!="NOTEE") {
    cout << "Wrong Category" << endl;
    return;
  }
  reader->_categoryEBEE = categoryEBEE;
  cout << "RUNNING ON CATEGORY " << categoryEBEE << endl;

  if ((sample=="PhotonJet_Pt0to15")||(sample=="PhotonJet_Pt15to20")||(sample=="PhotonJet_Pt20to30")||(sample=="PhotonJet_Pt30to50")||(sample=="PhotonJet_Pt50to80")||(sample=="PhotonJet_Pt80to120")||(sample=="PhotonJet_Pt120to170")||(sample=="PhotonJet_Pt170to300")||(sample=="PhotonJet_Pt300to500")||(sample=="PhotonJet_Pt500toInf")||(sample=="G_Pt_0to15_TuneZ2_7TeV_pythia6")||(sample=="G_Pt_15to30_TuneZ2_7TeV_pythia6")||(sample=="G_Pt_30to50_TuneZ2_7TeV_pythia6")||(sample=="G_Pt_50to80_TuneZ2_7TeV_pythia6")||(sample=="G_Pt_80to120_TuneZ2_7TeV_pythia6")||(sample=="G_Pt_120to170_TuneZ2_7TeV_pythia6")||(sample=="G_Pt_170to300_TuneZ2_7TeV_pythia6")||(sample=="G_Pt_300to470_TuneZ2_7TeV_pythia6")||(sample=="G_Pt_470to800_TuneZ2_7TeV_pythia6")||(sample=="G_Pt_800to1400_TuneZ2_7TeV_pythia6")||(sample=="G_Pt_1400to1800_TuneZ2_7TeV_pythia6")||(sample=="G_Pt_1800_TuneZ2_7TeV_pythia6")||(sample=="G_Pt_0to15")||(sample=="G_Pt_15to30")||(sample=="G_Pt_30to50")||(sample=="G_Pt_50to80")||(sample=="G_Pt_80to120")||(sample=="G_Pt_120to170")||(sample=="G_Pt_170to300")||(sample=="G_Pt_300to470")||(sample=="G_Pt_470to800")||(sample=="G_Pt_800to1400")||(sample=="G_Pt_1400to1800")||(sample=="G_Pt_1800") || (sample=="G_pt30to50")|| (sample=="G_pt50to80")|| (sample=="G_pt80to120")|| (sample=="G_pt120to170")|| (sample=="G_pt170to300")|| (sample=="G_pt300to470")|| (sample=="G_pt470to800")|| (sample=="G_pt800to1400")|| (sample=="G_pt1400to1800")|| (sample=="G_pt1800toInf") || (sample=="G_2EM_Pt20toInf"))
    {
      cout << "Filtering on GEN quantities" << endl;
      reader->_filterGen = kTRUE;
    }

  if ( (sample=="DiPhotonBorn_Pt10to25")||(sample=="DiPhotonBorn_Pt25to250")||(sample=="DiPhotonBorn_Pt250toInf")||(sample=="DiPhotonBox_Pt10to25")||(sample=="DiPhotonBox_Pt25to250")||(sample=="DiPhotonBox_Pt250toInf") )
    {
      cout << "Applying Diphoton KFactor" << endl;
      reader->_Kfactor = kTRUE;
    }
  
  reader->_outF = outF;
  reader->Loop();
  
  if (isData) {

    TChain *chain_tt = new TChain("diphotonAnalyzer/fTree");
    chain_tt->Add(infile.Data());
    cout << "TT entries = " << chain_tt->GetEntries() <<endl;
    
    TChain *chain_tf = new TChain("diphotonAnalyzer/fTightFakeTree");
    chain_tf->Add(infile.Data());
    cout << "TF entries = " << chain_tf->GetEntries() <<endl;

    TChain *chain_ft = new TChain("diphotonAnalyzer/fFakeTightTree");
    chain_ft->Add(infile.Data());
    cout << "FT entries = " << chain_ft->GetEntries() <<endl;
    
    TChain *chain_ff = new TChain("diphotonAnalyzer/fFakeFakeTree");
    chain_ff->Add(infile.Data());
    cout << "FF entries = " << chain_ff->GetEntries() <<endl;

    fTree*  readerTF = new fTree(chain_tf);
    readerTF->_fakeStatus = "TightFake";    
    readerTF->_cutPhoton1Pt = 70;
    readerTF->_cutPhoton2Pt = 70;
    readerTF->_reweightPU=kFALSE;
    readerTF->_categoryEBEE = categoryEBEE;
    outName = TString::Format("histograms_%s_TF_%s.root",sample.Data(),categoryEBEE.Data());
    outF = new TFile(outName.Data(),"RECREATE");
    readerTF->_outF = outF;

    readerTF->Loop(); 

    fTree*  readerFT = new fTree(chain_ft);
    readerFT->_fakeStatus = "FakeTight";    
    readerFT->_cutPhoton1Pt = 70;
    readerFT->_cutPhoton2Pt = 70;
    readerFT->_reweightPU=kFALSE;
    readerFT->_categoryEBEE = categoryEBEE;
    outName = TString::Format("histograms_%s_FT_%s.root",sample.Data(),categoryEBEE.Data());
    outF = new TFile(outName.Data(),"RECREATE");
    readerFT->_outF = outF;
    readerFT->Loop();     

    fTree*  readerFF = new fTree(chain_ff);
    readerFF->_fakeStatus = "FakeFake";    
    readerFF->_cutPhoton1Pt = 70;
    readerFF->_cutPhoton2Pt = 70;
    readerFF->_reweightPU=kFALSE;
    readerFF->_categoryEBEE = categoryEBEE;
    outName = TString::Format("histograms_%s_FF_%s.root",sample.Data(),categoryEBEE.Data());
    outF = new TFile(outName.Data(),"RECREATE");
    readerFF->_outF = outF;
    readerFF->Loop();     
    
  }

  // now for making pretty plots

  cout << "making pretty plots" << endl;
  TFile* fhists = new TFile(outName.Data(),"UPDATE");
  fhists->cd();

  // want to store the normalization histo
  if (!isData&&!isSignal) {
    //     tempName = "diphotonBkgAnalyzer/fNorm_h"; 
    tempName = "diphotonAnalyzer/fNorm_h"; 
    TH1F* h_norm = (TH1F*) f->Get(tempName.Data());
    if (h_norm) {
      cout << h_norm->GetEntries() << endl;
      h_norm->Write();
      fhists->Write();
    }
  }
  fhists->Close();
  
  cout << "draw individual histos" << endl;
  draw_individual_histos(sample,kPrint,kFALSE,categoryEBEE);

  cout << "fake rate histos" << endl;
  if (isData && sample!="DataMC" ) fake_rate_histos3(sample,kPrint,lumi,categoryEBEE);

  return;

}

void merge(TString sample = "PhotonJet",TString categoryEBEE="ALL"){

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat("ourme");

  //number of files in this set
  int temp_nfiles = 0;
   if (sample=="PhotonJet_filtered" || sample=="PhotonJet") {
     temp_nfiles = 10; 
   } else if (sample=="DiPhoton") {
     //     temp_nfiles = 6; 
     temp_nfiles = 4; 
   } else if (sample=="QCD_2EM") {
     temp_nfiles = 2; 
   } else if (sample=="QCDDiJet") {
     temp_nfiles = 20; 
   } else if (sample=="QCD_EMEnriched") {
     temp_nfiles = 3; 
   } else  if (sample=="G") {
     temp_nfiles = 10; 
   } else  if (sample=="G_2EM") {
     temp_nfiles = 1; 
   } else if (sample=="QCD") {
     temp_nfiles = 12; 
   } else if (sample=="DYtoEE") {
     temp_nfiles = 1; 
   } else {
     cout << "Wrong sample name!" << endl;
     return;
   }

  const int nfiles = temp_nfiles;

  //number of gen events
  int ngen[nfiles];
  if (sample=="DiPhoton") {
    
    //        ngen[0] = 537445;
    //        ngen[1] = 546355;
    //        ngen[2] = 777725;
    //        ngen[3] = 789470;

    ngen[0] = 532864;
    ngen[1] = 526240;
    ngen[2] = 518288;
    ngen[3] = 514514;
    
  } else if (sample=="QCD_2EM") {
    
    ngen[0] = 2700408;
    ngen[1] = 21276029;
    
  } else  if (sample=="G_2EM") {
    
    ngen[0] = 1182075;
    
  } else if (sample=="DYtoEE") {
    
    ngen[0] = 2059987;
    
  } else {
    cout << "Wrong sample name!" << endl;
    return;
  }
  
  const double lumi = 1.0; // units of pb-1

  //  TString ntupleDir = "rfio:/castor/cern.ch/user/t/torimoto/physics/diphoton/ntuples/mc/";
  TString ntupleDir = ".";
  
  TString labels[nfiles];

  //cross-sections (units of pb)
  // from  https://twiki.cern.ch/twiki/bin/viewauth/CMS/ProductionReProcessingSpring10
  double xsecs[nfiles];  

  // handy labels for each file, and for the overall set
   if (sample=="PhotonJet_filtered" || sample=="PhotonJet") {

     labels[0] = "Pt0to15";
     labels[1] = "Pt15to20";
     labels[2] = "Pt20to30";
     labels[3] = "Pt30to50";
     labels[4] = "Pt50to80";
     labels[5] = "Pt80to120";
     labels[6] = "Pt120to170";
     labels[7] = "Pt170to300";
     labels[8] = "Pt300to500";
     labels[9] = "Pt500toInf";

     xsecs[0] = 8.446e+07;
     xsecs[1] = 1.147e+05 ;
     xsecs[2] = 5.718e+04 ;
     xsecs[3] = 1.652e+04 ;
     xsecs[4] = 2.723e+03 ;
     xsecs[5] = 4.462e+02 ;
     xsecs[6] = 8.443e+01 ;
     xsecs[7] = 2.255e+01 ;
     xsecs[8] = 1.545e+00 ;
     xsecs[9] = 9.230e-02 ;    

     // assuming all the files follow a consistent naming scheme
     // where the only thing that changes is the 'label' as defined above
   } else if (sample=="DiPhoton") {

//      labels[0] = "Born_Pt10to25"; // Born
//      labels[1] = "Born_Pt25to250"; // Born
//      labels[2] = "Born_Pt250toInf"; // Born
//      labels[3] = "Box_Pt10to25"; // Box
//      labels[4] = "Box_Pt25to250"; // Box
//      labels[5] = "Box_Pt250toInf"; // Box

//      xsecs[0] = 236.4 ;
//      xsecs[1] = 22.37 ;
//      xsecs[2] = 8.072e-03 ;
//      xsecs[3] = 358.2 ;
//      xsecs[4] = 12.37 ;
//      xsecs[5] = 2.08e-04 ;

     labels[0] = "Born_Pt25to250"; // Born
     labels[1] = "Born_Pt250toInf"; // Born
     labels[2] = "Box_Pt25to250"; // Box
     labels[3] = "Box_Pt250toInf"; // Box

     xsecs[0] = 22.37 ;
     xsecs[1] = 8.072e-03 ;
     xsecs[2] = 12.37 ;
     xsecs[3] = 2.08e-04 ;

   } else if (sample=="QCDDiJet") {

     labels[0] = "Pt0to15";
     labels[1] = "Pt15to20";
     labels[2] = "Pt20to30";
     labels[3] = "Pt30to50";
     labels[4] = "Pt50to80";
     labels[5] = "Pt80to120";
     labels[6] = "Pt120to170";
     labels[7] = "Pt170to230";
     labels[8] = "Pt230to300";
     labels[9] = "Pt300to380";
     labels[10] = "Pt380to470";
     labels[11] = "Pt470to600";
     labels[12] = "Pt600to800";
     labels[13] = "Pt800to1000";
     labels[14] = "Pt1000to1400";
     labels[15] = "Pt1400to1800";
     labels[16] = "Pt1800to2200";
     labels[17] = "Pt2200to2600";
     labels[18] = "Pt2600to3000";
     labels[19] = "Pt3000to3500";

     xsecs[0] = 4.844e+10;
     xsecs[1] = 5.794e+08;
     xsecs[2] = 2.361e+08;
     xsecs[3] = 5.311e+07;
     xsecs[4] = 6.358e+06;
     xsecs[5] = 7.849e+05;
     xsecs[6] = 1.151e+05;
     xsecs[7] = 2.014e+04;
     xsecs[8] = 4.094e+03;
     xsecs[9] = 9.346e+02;
     xsecs[10] = 2.338e+02;
     xsecs[11] = 7.021e+01;
     xsecs[12] = 1.557e+01;
     xsecs[13] = 1.843e+00;
     xsecs[14] = 3.318e-01;
     xsecs[15] = 1.086e-02;
     xsecs[16] = 3.499e-04;
     xsecs[17] = 7.549e-06;
     xsecs[18] = 6.465e-08;
     xsecs[19] = 6.295e-11;

   } else if (sample=="QCD_EMEnriched") {

     labels[0] = "Pt20to30";
     labels[1] = "Pt30to80";
     labels[2] = "Pt80to170";

     xsecs[0] = 0.0073*(235.5e+06);
     xsecs[1] = 0.059*(59.3e+06);
     xsecs[2] = 0.148*(906000.0);

   } else if (sample=="QCD_2EM") {

     labels[0] = "pt30to40";
     labels[1] = "pt40toInf";

     xsecs[0] = 0.00023*(41800000.);
     xsecs[1] = 0.00216*(18700000.);

   } else if (sample=="G") {
  
    labels[0] = "pt30to50";	   
    labels[1] = "pt50to80";	   
    labels[2] = "pt80to120";   
    labels[3] = "pt120to170";  
    labels[4] = "pt170to300";  
    labels[5] = "pt300to470";  
    labels[6] = "pt470to800";  
    labels[7] = "pt800to1400"; 
    labels[8] = "pt1400to1800";
    labels[9] = "pt1800toInf"; 

    xsecs[0] = 16690;	     
    xsecs[1] = 2722;	     
    xsecs[2] = 447.2;	     
    xsecs[3] = 84.17;	     
    xsecs[4] = 22.64;	     
    xsecs[5] = 1.493;	     
    xsecs[6] = 0.1323;	     
    xsecs[7] = 0.003481;    
    xsecs[8] = 0.0000127;   
    xsecs[9] = 0.0000002936;

   } else if (sample=="G_2EM") {

     labels[0] = "Pt20toInf";

     xsecs[0] = 0.0064*(77100);

  } else if (sample=="QCD") {

    labels[0]  = "Pt30to50";	
    labels[1]  = "Pt50to80";	
    labels[2]  = "Pt80to120";	
    labels[3]  = "Pt120to170";	
    labels[4]  = "Pt170to300";	
    labels[5]  = "pt300to470";	
    labels[6]  = "Pt470to600";	
    labels[7]  = "Pt600to800";	
    labels[8]  = "Pt800to1000";	
    labels[9]  = "Pt1000to1400";	
    labels[10] = "Pt1400to1800";	
    labels[11] = "Pt1800toInf";  

    xsecs[0]  = 53120000; 
    xsecs[1]  = 6359000;  
    xsecs[2]  = 784300;	  
    xsecs[3]  = 115100;	  
    xsecs[4]  = 24260;	  
    xsecs[5]  = 1168;	  
    xsecs[6]  = 70.22;	  
    xsecs[7]  = 15.55;	  
    xsecs[8]  = 1.844;	  
    xsecs[9]  = 0.3321;	  
    xsecs[10] = 0.01087;  
    xsecs[11] = 0.0003575;

   } else if (sample=="DYtoEE") {

     labels[0] = "M20"; 

     xsecs[0] = 1300.0;

 }


  TString fileNames[nfiles];
  
  for(int ifile=0;ifile<nfiles;ifile++) {
    if (sample=="PhotonJet" || sample=="QCDDiJet" || sample=="QCD_EMEnriched" || sample=="G" || sample=="QCD" ) {
      //    fileNames[ifile] = TString::Format("%s/%s_%s/histograms_%s_%s.root",ntupleDir.Data(),sample.Data(),labels[ifile].Data(),sample.Data(),labels[ifile].Data());
      fileNames[ifile] = TString::Format("%s/histograms_%s_%s_%s.root",ntupleDir.Data(),sample.Data(),labels[ifile].Data(),categoryEBEE.Data());
    } else if (sample=="DiPhoton") {
      fileNames[ifile] = TString::Format("%s/histograms_%s%s_%s.root",ntupleDir.Data(),sample.Data(),labels[ifile].Data(),categoryEBEE.Data());
    } else {
      fileNames[ifile] = TString::Format("%s/histograms_%s_%s_%s.root",ntupleDir.Data(),sample.Data(),labels[ifile].Data(),categoryEBEE.Data());
      
    }
        cout << "File name = " << fileNames[ifile] <<endl;
  }
  
  //get the normalisation numbers from the norm hist
  int orig_events[nfiles];
  int npass[nfiles];

  TFile *ftemp[nfiles];
  TH1F *norm_h[nfiles];
  for(int ifile=0;ifile<nfiles;ifile++) { 
    ftemp[ifile] = TFile::Open(fileNames[ifile].Data());
    norm_h[ifile] = (TH1F*) ftemp[ifile]->Get("fNorm_h");

    //    if (norm_h[ifile]) {
    //      orig_events[ifile] = norm_h[ifile]->GetEntries();
    //      npass[ifile] = norm_h[ifile]->GetBinContent(norm_h[ifile]->FindBin(1));    
    //    } else {
    orig_events[ifile] = ngen[ifile];
    //    }
    //    orig_events[ifile] = 1.0;

    //    cout << "File name = " << fileNames[ifile];
    //    cout << "; Orig = " << orig_events[ifile]<< "; npass = " <<npass[ifile] <<endl;
    cout << fileNames[ifile] << " " << xsecs[ifile] << " " << orig_events[ifile] << " " << orig_events[ifile]/xsecs[ifile] << endl;
    //    ftemp[ifile]->Close();  

  }
  
  // instantiate the sum histos with histo in first file
  //  for (int i=0; i<nHists; i++) {
  //    allHistos[i] = (TH1F*)ftemp[0]->Get(nameHists[i].Data());
  //    allHistos[i]->Sumw2();
  //    allHistos[i]->Scale((1.0/orig_events[0])*xsecs[0]*lumi); 
  //  }

  // each individual one
  TH1F* indHistos[nHists][nfiles];

  Double_t totalError = 0;

  for(int ifile=0;ifile<nfiles;ifile++) { 
    for (int i=0; i<nHists; i++) {
      indHistos[i][ifile] = (TH1F*)ftemp[ifile]->Get(nameHists[i].Data());

      if (i==1) {
	//  error N_tot = sqrt(C1^2 * n1 + C2^2 * n2 + ...)        
	Double_t ci = xsecs[ifile]/orig_events[ifile];
	//	totalError = totalError + ci*ci*(indHistos[i][ifile]->Integral());
	totalError = totalError + ci*TMath::Sqrt(indHistos[i][ifile]->Integral());
	cout << sample.Data() << "_" << labels[ifile].Data() << " "  << ifile << " " << xsecs[ifile] << " " << orig_events[ifile] << " " 
	     << indHistos[i][ifile]->Integral() << " => " << ci << " " << ci*ci << " " << ci*ci*indHistos[i][ifile]->Integral() << " " << totalError << endl;	
      }

      indHistos[i][ifile]->Sumw2();
      indHistos[i][ifile]->Scale((1.0/orig_events[ifile])*xsecs[ifile]*lumi); 
      cout << "xsec orig lumi " << xsecs[ifile] << " " << orig_events[ifile] << " " << lumi << endl;
      //      cout << indHistos[i][ifile]->GetEntries() << endl;
    }
  }
  //  totalError = TMath::Sqrt(totalError);
  cout << sample.Data() << " Total Error = " << totalError << endl;
  
  TH1F* allHistos[nHists];

  // instantiate the sum histos with histo in first file
  for (int i=0; i<nHists; i++) {
    allHistos[i] = indHistos[i][0];
  }
  
  for(int ifile=0;ifile<nfiles;ifile++) {  
    for (int i=0; i<nHists; i++) {
      if (ifile>0) allHistos[i]->Add(indHistos[i][ifile]); // note we exclude first file since we used it to instantiate the sum histo
    }
  }

  for (int i=0; i<nHists; i++) {
    if (i==1) {
      cout <<  sample.Data() << " " << allHistos[i]->Integral()  << " +- " << totalError << endl ;
    }
  }
  
  TString outFileName = TString::Format("histograms_%s_all_%s.root",sample.Data(),categoryEBEE.Data());
  TFile fout(outFileName.Data(),"RECREATE");
  
  for (int i=0; i<nHists; i++) {
    allHistos[i]->Write();
  }
  //  fout.ls();
  cout << "Results written to: " << outFileName.Data() << endl;
  fout.Close();

  TString printLabel =  sample + "_all";

  draw_individual_histos(printLabel,kTRUE,kFALSE,categoryEBEE);
  
  return;

}



void mergeAllMC( Float_t lumiScaleFactor=1.0, TString categoryEBEE="ALL" ) {

  // lumi = lumiScaleFactor * 1/pb

  Double_t kFactor = 1.0;

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat("ourme");

  TString ntupleDir = ".";

  //  const int nfiles = 4;
  const int nfiles = 3;
  TString fileNames[nfiles];
  TString labels[nfiles];

  labels[0] = "G_2EM";
  //  labels[0] = "PhotonJet";
  //  labels[1] = "PhotonJet_filtered";
  //  labels[2] = "QCDDiJet";
  //  labels[1] = "QCD_EMEnriched";
  labels[1] = "QCD_2EM";
  labels[2] = "DiPhoton";
  //  labels[3] = "DYtoEE";

  Color_t fillColors[nfiles];
  fillColors[0] = kGreen+2;
  //  fillColors[2] = kRed+1;
  fillColors[1] = 6;
  fillColors[2] = kBlue; 
  //  fillColors[3] = 94; 

  TFile *ftemp[nfiles];

  for(int ifile=0;ifile<nfiles;ifile++) {    
    //    fileNames[ifile] = TString::Format("%s/%s_all/histograms_%s_all.root",ntupleDir.Data(),labels[ifile].Data(),labels[ifile].Data());
    fileNames[ifile] = TString::Format("histograms_%s_all_%s.root",labels[ifile].Data(),categoryEBEE.Data());
    ftemp[ifile] = TFile::Open(fileNames[ifile].Data());
    cout << "File name = " << fileNames[ifile] << endl;
  }

  // each individual one
  TH1F* indHistos[nHists][nfiles];

  for(int ifile=0;ifile<nfiles;ifile++) { 
    for (int i=0; i<nHists; i++) {
      indHistos[i][ifile] = (TH1F*)ftemp[ifile]->Get(nameHists[i].Data());
      indHistos[i][ifile]->Sumw2();
      indHistos[i][ifile]->SetFillColor(fillColors[ifile]);
      if (ifile==2) { 
	cout << "SCALING BY K-FACTOR: " << kFactor << " " << lumiScaleFactor << " " << kFactor*lumiScaleFactor << endl;
	indHistos[i][ifile]->Scale(lumiScaleFactor*kFactor);
	//	indHistos[i][ifile]->Scale(lumiScaleFactor*kFactor/0.556355);
	//	indHistos[i][ifile]->Scale(lumiScaleFactor/0.443645);
      } else {
	indHistos[i][ifile]->Scale(lumiScaleFactor);
	//	indHistos[i][ifile]->Scale(lumiScaleFactor/0.556355);
	//	indHistos[i][ifile]->Scale(lumiScaleFactor/0.443645);
      }
    }
  }

  //  TH1F* allHistos[nHists];

  THStack *allHistos[nHists];

  // instantiate the sum histos with histo in first file
  for (int i=0; i<nHists; i++) {
    //    allHistos[i] = indHistos[i][0];
    //    indHistos[i][ifile] = (TH1F*)ftemp[ifile]->Get(nameHists[i].Data());    
    TString hname1 = indHistos[i][0]->GetTitle();
    TString hname2 = indHistos[i][0]->GetXaxis()->GetTitle();
    TString hname = hname1 + ";" + hname2;
    allHistos[i] = new THStack(nameHists[i].Data(),hname.Data());
  }

  float sumMC = 0;
  for(int ifile=0;ifile<nfiles;ifile++) {  
    for (int i=0; i<nHists; i++) {
      //      if (ifile>0) allHistos[i]->Add(indHistos[i][ifile]); // note we exclude first file since we used it to instantiate the sum histo
      allHistos[i]->Add(indHistos[i][ifile]); 
      if (i==1) {
	cout << ifile << " " << labels[ifile].Data() << " " << indHistos[i][ifile]->Integral() << endl;	
	sumMC = sumMC + indHistos[i][ifile]->Integral();
      }
    }
  }
  cout << "Total " << sumMC << endl;

  TString sample = TString::Format("allMC_%0.1fpb",lumiScaleFactor);

  cout << sample << endl;

  TString outFileName = TString::Format("histograms_%s_%s.root",sample.Data(),categoryEBEE.Data());
  TFile fout(outFileName.Data(),"RECREATE");

  cout << outFileName << endl;

  for (int i=0; i<nHists; i++) {
    allHistos[i]->Write();
  }
  //  fout.ls();
  cout << "Results written to: " << outFileName.Data() << endl;
  fout.Close();

  draw_individual_histos(sample,kTRUE,kTRUE,categoryEBEE);

  return;

}


void overlayDataMC(  TString datalabel = "Apr29_43pb",  TString mclabel = "allMC_43.0pb", TString categoryEBEE="ALL", Bool_t scaleMC2DATA = kFALSE){
  
  gROOT->SetStyle("Plain");
  //  gStyle->SetOptStat("ourme");
  gStyle->SetOptStat("ourme");

  //  TString ntupleDir = "rfio:/castor/cern.ch/user/t/torimoto/physics/diphoton/ntuples/mc/";
  TString ntupleDir = ".";
  
  const int nfiles = 2;
  TString fileNames[nfiles];

  TString labels[nfiles];

  labels[1]=mclabel;
  labels[0] = "data_"+datalabel;

  TFile *ftemp[nfiles];

  for(int ifile=0;ifile<nfiles;ifile++) {    
    fileNames[ifile] = TString::Format("histograms_%s_%s.root",labels[ifile].Data(),categoryEBEE.Data());
    ftemp[ifile] = TFile::Open(fileNames[ifile].Data());
    cout << "File name = " << fileNames[ifile] << endl;
  }
  TString sample = "DataMC";

  TCanvas* c[nHists];
  char cname[100];   

  TH1F* histosData[nHists];
  THStack* histosMC[nHists];
  float maxData[nHists];
  float maxMC[nHists];

  for (int i=0; i<nHists; i++) {
    //    if (i==99) continue;
    if (i==55) continue;
    if (i==56) continue;
    //    cout << i << " " << nameHists[i].Data() << endl;
    // Data
    histosData[i] = (TH1F*)ftemp[0]->Get(nameHists[i].Data());      
    TString histoTitle = histosData[i]->GetTitle();
    //    histoTitle = histoTitle + " (" + sample + ")";
    //    histosData[i]->SetTitle(histoTitle.Data());	       
    histosData[i]->SetTitle("");
    histosData[i]->SetMarkerStyle(20);

    if ((nameHists[i]=="h_Diphoton_Minv_log")||(nameHists[i]=="h_Diphoton_Minv")||(nameHists[i]=="h_Diphoton_Minv_low")||(nameHists[i]=="h_Diphoton_Minv_high")) {
      histosData[i]->SetMarkerSize(0.8);
    } else {
      histosData[i]->SetMarkerSize(1.0);
    }
    maxData[i] = histosData[i]->GetMaximum();
    // MC
    histosMC[i] = (THStack*)ftemp[1]->Get(nameHists[i].Data());      
    histosMC[i]->SetTitle("");
    maxMC[i] = histosMC[i]->GetMaximum();  
  }

  THStack *histosMC_scaledToData[nHists];

  for (int i=0; i<nHists; i++) {
    if (i==54) continue;
    if (i==25) continue;
    if (i==55) continue;
    if (i==56) continue;
    //    cout << i << " " << nameHists[i].Data() << endl;
    sprintf(cname,"c%i",i);
    // cout << i << " " << cname << endl;
    if (i==0) {
      c[i] = new TCanvas(cname, nameHists[i].Data(), 1000, 600); 
      c[i]->SetBottomMargin(0.4);
      histosData[i]->GetXaxis()->SetTitleSize(0.01);
      histosData[i]->GetXaxis()->SetBit(TAxis::kLabelsVert);
    } else {
      c[i] = new TCanvas(cname, nameHists[i].Data(), 600, 600);
    }
    c[i]->cd();

    if (scaleMC2DATA) {
      
      int nBinsData = histosData[i]->GetNbinsX();
      int nData = histosData[i]->Integral(1,nBinsData); //histosData[i]->GetEntries();
      
      TList* mcHists = histosMC[i]->GetHists();
      //	mcHists->Print();
      TIter next (mcHists);

      //      const int nsamp=4; 
      const int nsamp=3; 
      float nMC = 0;
      TH1F* histMC[nsamp];
      
      for (int ij=0; ij<nsamp; ij++) {
	histMC[ij] = (TH1F*)next();
	int nBinsMC = histMC[ij]->GetNbinsX();
	nMC += histMC[ij]->Integral(1,nBinsMC);
      }
      
      //	cout << nData << " " << nMC << " " << nameHists[i].Data() << endl;
      float scaleMC = 1;
      if (nMC!=0) {
	scaleMC = 1.0*nData/nMC;
      } else {
	scaleMC = 0;
      }
      cout << i << nameHists[i].Data() << " DATA / MC = SCALE " << nData << " " << nMC << " " << scaleMC << endl;
      
      char fname[100];     	
      histosMC_scaledToData[i] = new THStack(fname,"");
      sprintf(fname,"%s_scaledToData",histosMC[i]->GetName());
      
      for (int ij=0; ij<nsamp; ij++) {
	histMC[ij]->Scale(scaleMC);
	histosMC_scaledToData[i]->Add(histMC[ij]);
      }		

    }

    if (scaleMC2DATA) {

      histosMC_scaledToData[i]->SetMaximum(TMath::Max(histosMC_scaledToData[i]->GetMaximum(),histosData[i]->GetMaximum()));
      histosData[i]->SetMaximum(TMath::Max(histosMC_scaledToData[i]->GetMaximum(),histosData[i]->GetMaximum()));
      
      histosMC_scaledToData[i]->Draw("hist");        
      histosMC_scaledToData[i]->GetXaxis()->SetTitle(histosData[i]->GetXaxis()->GetTitle());
      histosData[i]->Draw("e same");              
      c[i]->Update();

    } else {
      histosMC[i]->SetMaximum(TMath::Max(histosMC[i]->GetMaximum(),histosData[i]->GetMaximum()));
      histosData[i]->SetMaximum(TMath::Max(histosMC[i]->GetMaximum(),histosData[i]->GetMaximum()));

      //    histosData[i]->Draw("e");              
      
      histosMC[i]->Draw("hist");        
      histosMC[i]->GetXaxis()->SetTitle(histosData[i]->GetXaxis()->GetTitle());
      histosData[i]->Draw("e same");              
      c[i]->Update();
    }

    //    }    
    if ((nameHists[i]=="h_Photon1_hadOverEm")||(nameHists[i]=="h_Photon1_hcalIso04")||(nameHists[i]=="h_Photon1_hcalIso03")||(nameHists[i]=="h_Photon1_trkIsoSumPtHollow04")||(nameHists[i]=="h_Photon1_trkIsoSumPtSolid04")||(nameHists[i]=="h_Photon1_trkIsoSumPtHollow03")||(nameHists[i]=="h_Photon1_trkIsoSumPtSolid03")||(nameHists[i]=="h_Photon1_maxRecHitTime_wide")||(nameHists[i]=="h_Photon1_e2x2e4x4")||(nameHists[i]=="h_Photon1_e2e9")||(nameHists[i]=="h_Photon2_hadOverEm")||(nameHists[i]=="h_Photon2_hcalIso04")||(nameHists[i]=="h_Photon2_hcalIso03")||(nameHists[i]=="h_Photon2_trkIsoSumPtHollow04")||(nameHists[i]=="h_Photon2_trkIsoSumPtSolid04")||(nameHists[i]=="h_Photon2_trkIsoSumPtHollow03")||(nameHists[i]=="h_Photon2_trkIsoSumPtSolid03")||(nameHists[i]=="h_Photon2_maxRecHitTime_wide")||(nameHists[i]=="h_Photon2_e2x2e4x4")||(nameHists[i]=="h_Photon2_e2e9")||(nameHists[i]=="h_Diphoton_Minv_log")||(nameHists[i]=="h_Photon1_pt_log")||(nameHists[i]=="h_Photon2_pt_log")||(nameHists[i]=="h_Photon1_pt_zoom")||(nameHists[i]=="h_Photon2_pt_zoom")) { 
      c[i]->SetLogy();   
    }
        
    TList* toyo = histosMC[i]->GetHists();
    TIter next (toyo);
    TH1F* h_TF = (TH1F*)next(); 
    TH1F* h_FF = (TH1F*)next(); 
    TH1F* h_Diphoton = (TH1F*)next(); 
    //    TH1F* h_DY = (TH1F*)next(); 

    // legend
    TLegend *leg = new TLegend(0.5234899,0.5961538,0.6627517,0.8076923,NULL,"brNDC");
    leg->SetBorderSize(1);
    leg->SetTextFont(62);
    leg->SetTextSize(0.03225806);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    leg->AddEntry(histosData[i],"data");
    //    leg->AddEntry(h_DY,"Drell-Yan e+e- MC");
    leg->AddEntry(h_Diphoton,"Diphoton MC");
    leg->AddEntry(h_FF,"Dijets EM-Enriched MC");
    leg->AddEntry(h_TF,"#gamma jets MC");
    if (nameHists[i]=="h_Photon1_pt"||nameHists[i]=="h_Photon2_pt"||nameHists[i]=="h_Diphoton_Minv"||nameHists[i]=="h_Diphoton_Minv_log"||nameHists[i]=="h_Diphoton_Minv_low"||nameHists[i]=="h_Diphoton_Minv_high") {
      leg->Draw();
    }
    char fname[100];     
    sprintf(fname,"%s_%s.png",nameHists[i].Data(),sample.Data());  
    c[i]->Print(fname);
    sprintf(fname,"%s_%s.pdf",nameHists[i].Data(),sample.Data());  
    c[i]->Print(fname);
    //    sprintf(fname,"%s_%s.C",nameHists[i].Data(),sample.Data());  
    //    c[i]->Print(fname);
  }

  TString outFileName = TString::Format("histograms_%s_%s_%s.root",sample.Data(),datalabel.Data(),categoryEBEE.Data());
  TFile fout(outFileName.Data(),"RECREATE");

  // cumulative background plot
  TH1F* histoData = (TH1F*)ftemp[0]->Get("h_Diphoton_Minv");
  
  THStack* stackMC;

  if (scaleMC2DATA) {
    stackMC = histosMC_scaledToData[2];
  } else {
    stackMC = (THStack*) ftemp[1]->Get("h_Diphoton_Minv");
  }  
  TCanvas* cCumul = cumulativeBackground(histoData, stackMC,3);
  //  cCumul->SetLogy();
  cCumul->Write();
  char fname[100];     
  sprintf(fname,"%s_%s.png","h_CumulativeBackground",sample.Data());  
  cCumul->Print(fname);
  sprintf(fname,"%s_%s.pdf","h_CumulativeBackground",sample.Data());  
  cCumul->Print(fname);
  
  //  cout << "writing file.." << endl;
  for (int i=0; i<nHists; i++) {
    cout << i << " " << nameHists[i].Data() << endl;
    if (i==54) continue;
    if (i==25) continue;
    if (i==55) continue;
    if (i==56) continue;
    c[i]->Write();
  }
  
  //  fout.ls();
  cout << "Results written to: " << outFileName.Data() << endl;
  fout.Close();

  return;

}




void mergeAllData(TString categoryEBEE="ALL") {

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat("ourme");

  TString ntupleDir = ".";

  const int nfiles = 8;
  TString fileNames[nfiles];
  TString labels[nfiles];

  labels[0] = "38X_sept17rereco";
  labels[1] = "EG_Run2010A_Sep17ReReco_Oct5update";
  labels[2] = "Photon_Run2010B_PromptReco_146428_146644_Oct4";
  labels[3] = "Photon_Run2010B_146729_147116_Oct8th";

  labels[4] = "Photon_Run2010B_147115_147454_3978nb";
  labels[5] = "Photon_Run2010B_Oct22_4150nb";
  labels[6] = "Photon_Run2010B_promptreco_oct29_6807nb";
  labels[7] = "Nov5_12831nb";

  TFile *ftemp[nfiles];

  for(int ifile=0;ifile<nfiles;ifile++) {    
    //    fileNames[ifile] = TString::Format("%s/%s_all/histograms_%s_all.root",ntupleDir.Data(),labels[ifile].Data(),labels[ifile].Data());
    fileNames[ifile] = TString::Format("histograms_data_%s_%s.root",labels[ifile].Data(),categoryEBEE.Data());
    ftemp[ifile] = TFile::Open(fileNames[ifile].Data());
    cout << "File name = " << fileNames[ifile] << endl;
  }

  // each individual one
  TH1F* indHistos[nHists][nfiles];

  for(int ifile=0;ifile<nfiles;ifile++) { 
    for (int i=0; i<nHists; i++) {
      indHistos[i][ifile] = (TH1F*)ftemp[ifile]->Get(nameHists[i].Data());
    }
  }

  TH1F* allHistos[nHists];

  // instantiate the sum histos with histo in first file
  for (int i=0; i<nHists; i++) {
    allHistos[i] = indHistos[i][0];
  }

  for(int ifile=0;ifile<nfiles;ifile++) {  
    for (int i=0; i<nHists; i++) {
      if (ifile>0) allHistos[i]->Add(indHistos[i][ifile]); // note we exclude first file since we used it to instantiate the sum histo
    }
  }

  TString sample = TString::Format("allData");

  cout << sample << endl;

  TString outFileName = TString::Format("histograms_%s_%s.root",sample.Data(),categoryEBEE.Data());
  TFile fout(outFileName.Data(),"RECREATE");

  cout << outFileName << endl;

  for (int i=0; i<nHists; i++) {
    allHistos[i]->Write();
  }
  //  fout.ls();
  cout << "Results written to: " << outFileName.Data() << endl;
  fout.Close();

  //  draw(sample);
  draw_individual_histos(sample,kTRUE,kFALSE);

  return;

}



