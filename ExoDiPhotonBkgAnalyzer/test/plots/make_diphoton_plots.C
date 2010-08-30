
void make_diphoton_plots(TString sample = "PhotonJet_Pt30to50", Bool_t kPrint=kTRUE, Bool_t isData=kFALSE) {

  TString infile;
  if (isData) {
    TString infile = "rfio:/castor/cern.ch/user/y/yma/RSGravitons/Aug2010/diphoton_tree_132440-141961.root";
  } else {
    TString infile = "rfio:/castor/cern.ch/user/t/torimoto/physics/diphoton/ntuples/mc/"+sample+"/diphotonTree_"+sample+".root";
  }
  cout << "Making quick check plots for: " << infile << endl;

  if (isData) {
    gSystem->Load("fTreeData_C.so");
  } else {
    gSystem->Load("fTree_C.so");
  }

  // input file
  TFile* f = TFile::Open(infile.Data());
  
  TString tempName;
  if (isData) { tempName = "diphotonAnalyzer/fTree";  } 
  else        { tempName = "diphotonBkgAnalyzer/fTree"; }
  TTree* fChain = (TTree*)f->Get(tempName.Data());

  // set up output
  TString outName = TString::Format("histograms_%s.root",sample.Data());
  TFile* outF = new TFile(outName.Data(),"RECREATE");

  // loop
  if (isData) {
    fTreeData* reader = new fTreeData(fChain);
    reader->_outF = outF;
    reader->Loop();
  } else {
    fTree* reader = new fTree(fChain);
    reader->_outF = outF;
    reader->Loop();
  }

  // now for making pretty plots

  TFile* fhists = new TFile(outName.Data(),"UPDATE");
  fhists->cd();

  // want to store the normalization histo
  if (!isData) {
    tempName = "diphotonBkgAnalyzer/fNorm_h"; 
    h_norm = (TH1F*) f->Get(tempName.Data());
    cout << h_norm->GetEntries() << endl;
    h_norm->Write();
    fhists->Write();
  }
  fhists->Close();

  draw(sample,kPrint);
  draw_individual_histos(sample,kPrint);

  return;

}


void merge(TString sample = "PhotonJet"){

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat("ourme");

  //number of files in this set
  int temp_nfiles = 0;
  if (sample=="PhotonJet") {
    temp_nfiles = 10; 
  } else if (sample=="DiPhoton") {
    temp_nfiles = 6; 
  } else {
    cout << "Wrong sample name!" << endl;
    return;
  }

  const double lumi = 1; // units of pb-1

  //  TString ntupleDir = "rfio:/castor/cern.ch/user/t/torimoto/physics/diphoton/ntuples/mc/";
  TString ntupleDir = ".";
  
  const int nfiles = temp_nfiles;
  TString labels[nfiles];
  double xsecs[nfiles];  


  // handy labels for each file, and for the overall set
  if (sample=="PhotonJet") {
  
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

    //cross-sections (units of pb)
    // from  https://twiki.cern.ch/twiki/bin/viewauth/CMS/ProductionReProcessingSpring10
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

    labels[0] = "Born_Pt10to25"; // Born
    labels[1] = "Born_Pt25to250"; // Born
    labels[2] = "Born_Pt250toInf"; // Born
    labels[3] = "Box_Pt10to25"; // Box
    labels[4] = "Box_Pt25to250"; // Box
    labels[5] = "Box_Pt250toInf"; // Box

    xsecs[0] = 236.4 ;
    xsecs[1] = 22.37 ;
    xsecs[2] = 8.072e-03 ;
    xsecs[3] = 358.2 ;
    xsecs[4] = 12.37 ;
    xsecs[5] = 2.08e-04 ;
  }

  TString fileNames[nfiles];
  
  for(int ifile=0;ifile<nfiles;ifile++) {
    if (sample=="PhotonJet") {
      //    fileNames[ifile] = TString::Format("%s/%s_%s/histograms_%s_%s.root",ntupleDir.Data(),sample.Data(),labels[ifile].Data(),sample.Data(),labels[ifile].Data());
      fileNames[ifile] = TString::Format("%s/histograms_%s_%s.root",ntupleDir.Data(),sample.Data(),labels[ifile].Data(),sample.Data(),labels[ifile].Data());
    } else if (sample=="DiPhoton") {
      fileNames[ifile] = TString::Format("%s/histograms_%s%s.root",ntupleDir.Data(),sample.Data(),labels[ifile].Data(),sample.Data(),labels[ifile].Data());
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

    orig_events[ifile] = norm_h[ifile]->GetEntries();
    npass[ifile] = norm_h[ifile]->GetBinContent(norm_h[ifile]->FindBin(1));

    cout << "File name = " << fileNames[ifile];
    cout << "; Orig = " << orig_events[ifile]<< "; npass = " <<npass[ifile] <<endl;
    //    ftemp[ifile]->Close();
  }

  
  const int nHists=52;
  TString nameHists[nHists] = { "h_TrigHLT", "h_Diphoton_Minv", "h_Diphoton_qt", "h_Diphoton_deltaPhi", "h_Diphoton_deltaEta", "h_Diphoton_deltaR", "h_Photon1_pt", "h_Photon1_eta", "h_Photon1_phi", "h_Photon1_r9", "h_Photon1_sigmaIetaIeta", "h_Photon1_sigmaEtaEta", "h_Photon1_swisscross", "h_Photon1_severityLevel", "h_Photon1_recHitFlag", "h_Photon1_maxRecHitTime", "h_Photon1_hadOverEm", "h_Photon1_hcalIso04", "h_Photon1_hcalIso03", "h_Photon1_ecalIso04", "h_Photon1_ecalIso03", "h_Photon1_trkIsoSumPtHollow04", "h_Photon1_trkIsoSumPtSolid04", "h_Photon1_trkIsoNtrksHollow04", "h_Photon1_trkIsoNtrksSolid04", "h_Photon1_trkIsoSumPtHollow03", "h_Photon1_trkIsoSumPtSolid03", "h_Photon1_trkIsoNtrksHollow03", "h_Photon1_trkIsoNtrksSolid03", "h_Photon2_pt", "h_Photon2_eta", "h_Photon2_phi", "h_Photon2_r9", "h_Photon2_sigmaIetaIeta", "h_Photon2_sigmaEtaEta", "h_Photon2_swisscross", "h_Photon2_severityLevel", "h_Photon2_recHitFlag", "h_Photon2_maxRecHitTime", "h_Photon2_hadOverEm", "h_Photon2_hcalIso04", "h_Photon2_hcalIso03", "h_Photon2_ecalIso04", "h_Photon2_ecalIso03", "h_Photon2_trkIsoSumPtHollow04", "h_Photon2_trkIsoSumPtSolid04", "h_Photon2_trkIsoNtrksHollow04", "h_Photon2_trkIsoNtrksSolid04", "h_Photon2_trkIsoSumPtHollow03", "h_Photon2_trkIsoSumPtSolid03", "h_Photon2_trkIsoNtrksHollow03", "h_Photon2_trkIsoNtrksSolid03"  };


  // instantiate the sum histos with histo in first file
  //  for (int i=0; i<nHists; i++) {
  //    allHistos[i] = (TH1F*)ftemp[0]->Get(nameHists[i].Data());
  //    allHistos[i]->Sumw2();
  //    allHistos[i]->Scale((1.0/orig_events[0])*xsecs[0]*lumi); 
  //  }

  // each individual one
  TH1F* indHistos[nHists][nfiles];

  for(int ifile=0;ifile<nfiles;ifile++) { 
    for (int i=0; i<nHists; i++) {
      indHistos[i][ifile] = (TH1F*)ftemp[ifile]->Get(nameHists[i].Data());
      indHistos[i][ifile]->Sumw2();
      indHistos[i][ifile]->Scale((1.0/orig_events[ifile])*xsecs[ifile]*lumi); 
      //      cout << indHistos[i][ifile]->GetEntries() << endl;
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

  TString outFileName = TString::Format("histograms_%s_all.root",sample.Data());
  TFile fout(outFileName.Data(),"RECREATE");
  
  for (int i=0; i<nHists; i++) {
    allHistos[i]->Write();
  }
  //  fout.ls();
  cout << "Results written to: " << outFileName.Data() << endl;
  fout.Close();

  TString printLabel =  sample + "_all";

  draw(printLabel);
  draw_individual_histos(printLabel);
  
  return;

}

// draws histograms organized into a few canvases
void draw(TString sample, Bool_t kPrint=kTRUE) {

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(10);
  
  TString outName = TString::Format("histograms_%s.root",sample.Data());
  TFile* fhists = new TFile(outName.Data());

  const int nHists=10;
  TCanvas* c[nHists];
  char cname[100]; 

  for (int i=0; i<nHists; i++) {
    sprintf(cname,"c%i",i);
    //    int x = (i%3)*600;     //int x = (i%3)*600;
    //    int y = (i/3)*100;     //int y = (i/3)*200;
    //    c[i] =  new TCanvas(cname,cname,x,y,600,400);
    c[i] = new TCanvas(cname, cname,1000,600);
    if (i>0) c[i]->Divide(3,2);
  }

  char fname[100]; 

  // trigger
  c[0]->cd();
  c0->SetBottomMargin(0.4);
  h_TrigHLT->GetXaxis()->SetTitleSize(0.01);
  h_TrigHLT->GetXaxis()->SetBit(TAxis::kLabelsVert);
  h_TrigHLT->Draw();
  sprintf(fname,"trigger_%s.png",sample.Data());  
  if (kPrint) c[0]->Print(fname);

  // photon 1
  c[1]->cd(1);
  h_Photon1_pt->Draw();
  c[1]->cd(2);
  h_Photon1_eta->Draw();
  c[1]->cd(3);
  h_Photon1_phi->Draw();
  c[1]->cd(4);
  h_Photon1_r9->Draw();
  c[1]->cd(5);
  h_Photon1_sigmaIetaIeta->Draw();
  c[1]->cd(6);
  h_Photon1_sigmaEtaEta->Draw();
  sprintf(fname,"kinematics_photon1_%s.png",sample.Data());  
  if (kPrint) c[1]->Print(fname);

  c[2]->cd(1);
  h_Photon1_swisscross->Draw();
  c[2]->cd(2);
  h_Photon1_severityLevel->Draw();
  c[2]->cd(3);
  h_Photon1_recHitFlag->Draw();
  c[2]->cd(4);
  h_Photon1_maxRecHitTime->Draw();
  sprintf(fname,"spike_photon1_%s.png",sample.Data());  
  if (kPrint) c[2]->Print(fname);

  c[3]->cd(1);	
  //  c[3]->GetPad(1)->SetLogy(1);
  c3_1->SetLogy(1);
  h_Photon1_hadOverEm->Draw();
  c[3]->cd(2);	
  c3_2->SetLogy(1);
  h_Photon1_hcalIso04->Draw();
  c[3]->cd(3);	
  c3_3->SetLogy(1);
  h_Photon1_hcalIso03->Draw();
  c[3]->cd(4);	
  h_Photon1_ecalIso04->Draw();
  c[3]->cd(5);	
  h_Photon1_ecalIso03->Draw();
  sprintf(fname,"caloIso_photon1_%s.png",sample.Data());  
  if (kPrint) c[3]->Print(fname);

  c[4] = new TCanvas("c4", "c4",1200,600);
  c[4]->Divide(4,2);
  c[4]->cd(1);	
  c4_1->SetLogy(1);
  h_Photon1_trkIsoSumPtHollow04->Draw();
  c[4]->cd(2);	
  c4_2->SetLogy(1);
  h_Photon1_trkIsoSumPtSolid04->Draw();
  c[4]->cd(3);	
  c4_3->SetLogy(1);
  h_Photon1_trkIsoSumPtHollow03->Draw();
  c[4]->cd(4);	
  c4_4->SetLogy(1);
  h_Photon1_trkIsoSumPtSolid03->Draw();
  c[4]->cd(5);	
  h_Photon1_trkIsoNtrksHollow04->Draw();
  c[4]->cd(6);	
  h_Photon1_trkIsoNtrksSolid04->Draw();
  c[4]->cd(7);	
  h_Photon1_trkIsoNtrksHollow03->Draw();
  c[4]->cd(8);	
  h_Photon1_trkIsoNtrksSolid03->Draw();
  sprintf(fname,"trackIso_photon1_%s.png",sample.Data());  
  if (kPrint) c[4]->Print(fname);

  // photon2

  c[5]->cd(1);
  h_Photon2_pt->Draw();
  c[5]->cd(2);
  h_Photon2_eta->Draw();
  c[5]->cd(3);
  h_Photon2_phi->Draw();
  c[5]->cd(4);
  h_Photon2_r9->Draw();
  c[5]->cd(5);
  h_Photon2_sigmaIetaIeta->Draw();
  c[5]->cd(6);
  h_Photon2_sigmaEtaEta->Draw();
  sprintf(fname,"kinematics_photon2_%s.png",sample.Data());  
  if (kPrint) c[5]->Print(fname);

  c[6]->cd(1);
  h_Photon2_swisscross->Draw();
  c[6]->cd(2);
  h_Photon2_severityLevel->Draw();
  c[6]->cd(3);
  h_Photon2_recHitFlag->Draw();
  c[6]->cd(4);
  h_Photon2_maxRecHitTime->Draw();
  sprintf(fname,"spike_photon2_%s.png",sample.Data());  
  if (kPrint) c[6]->Print(fname);

  c[7]->cd(1);	
  //  c[3]->GetPad(1)->SetLogy(1);
  c7_1->SetLogy(1);
  h_Photon2_hadOverEm->Draw();
  c[7]->cd(2);	
  c7_2->SetLogy(1);
  h_Photon2_hcalIso04->Draw();
  c[7]->cd(3);	
  c7_3->SetLogy(1);
  h_Photon2_hcalIso03->Draw();
  c[7]->cd(4);	
  h_Photon2_ecalIso04->Draw();
  c[7]->cd(5);	
  h_Photon2_ecalIso03->Draw();
  sprintf(fname,"caloIso_photon2_%s.png",sample.Data());  
  if (kPrint) c[7]->Print(fname);

  c[8] = new TCanvas("c8", "c8",1200,600);
  c[8]->Divide(4,2);
  c[8]->cd(1);	
  c8_1->SetLogy(1);
  h_Photon2_trkIsoSumPtHollow04->Draw();
  c[8]->cd(2);	
  c8_2->SetLogy(1);
  h_Photon2_trkIsoSumPtSolid04->Draw();
  c[8]->cd(3);	
  c8_3->SetLogy(1);
  h_Photon2_trkIsoSumPtHollow03->Draw();
  c[8]->cd(4);	
  c8_4->SetLogy(1);
  h_Photon2_trkIsoSumPtSolid03->Draw();
  c[8]->cd(5);	
  h_Photon2_trkIsoNtrksHollow04->Draw();
  c[8]->cd(6);	
  h_Photon2_trkIsoNtrksSolid04->Draw();
  c[8]->cd(7);	
  h_Photon2_trkIsoNtrksHollow03->Draw();
  c[8]->cd(8);	
  h_Photon2_trkIsoNtrksSolid03->Draw();
  sprintf(fname,"trackIso_photon2_%s.png",sample.Data());  
  if (kPrint) c[8]->Print(fname);

  // diphoton
  c[9] = new TCanvas("c9", "c9",1200,600);
  c[9]->Divide(4,2);
  c[9]->cd(1);	
  h_Diphoton_Minv->Draw();
  c[9]->cd(2);	
  h_Diphoton_qt->Draw();
  c[9]->cd(3);	
  h_Diphoton_deltaPhi->Draw();
  c[9]->cd(4);	
  h_Diphoton_deltaEta->Draw();
  c[9]->cd(5);	
  h_Diphoton_deltaR->Draw();
  sprintf(fname,"diphoton_%s.png",sample.Data());  
  if (kPrint) c[9]->Print(fname);

  return;
}


// prints histos on individual canvases
void draw_individual_histos(TString sample, Bool_t kPrint=kTRUE) {

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(10);
  
  TString outName = TString::Format("histograms_%s.root",sample.Data());
  TFile* fhists = new TFile(outName.Data());

  const int nHists=52;
  TString nameHists[nHists] = { "h_TrigHLT", "h_Diphoton_Minv", "h_Diphoton_qt", "h_Diphoton_deltaPhi", "h_Diphoton_deltaEta", "h_Diphoton_deltaR", "h_Photon1_pt", "h_Photon1_eta", "h_Photon1_phi", "h_Photon1_r9", "h_Photon1_sigmaIetaIeta", "h_Photon1_sigmaEtaEta", "h_Photon1_swisscross", "h_Photon1_severityLevel", "h_Photon1_recHitFlag", "h_Photon1_maxRecHitTime", "h_Photon1_hadOverEm", "h_Photon1_hcalIso04", "h_Photon1_hcalIso03", "h_Photon1_ecalIso04", "h_Photon1_ecalIso03", "h_Photon1_trkIsoSumPtHollow04", "h_Photon1_trkIsoSumPtSolid04", "h_Photon1_trkIsoNtrksHollow04", "h_Photon1_trkIsoNtrksSolid04", "h_Photon1_trkIsoSumPtHollow03", "h_Photon1_trkIsoSumPtSolid03", "h_Photon1_trkIsoNtrksHollow03", "h_Photon1_trkIsoNtrksSolid03", "h_Photon2_pt", "h_Photon2_eta", "h_Photon2_phi", "h_Photon2_r9", "h_Photon2_sigmaIetaIeta", "h_Photon2_sigmaEtaEta", "h_Photon2_swisscross", "h_Photon2_severityLevel", "h_Photon2_recHitFlag", "h_Photon2_maxRecHitTime", "h_Photon2_hadOverEm", "h_Photon2_hcalIso04", "h_Photon2_hcalIso03", "h_Photon2_ecalIso04", "h_Photon2_ecalIso03", "h_Photon2_trkIsoSumPtHollow04", "h_Photon2_trkIsoSumPtSolid04", "h_Photon2_trkIsoNtrksHollow04", "h_Photon2_trkIsoNtrksSolid04", "h_Photon2_trkIsoSumPtHollow03", "h_Photon2_trkIsoSumPtSolid03", "h_Photon2_trkIsoNtrksHollow03", "h_Photon2_trkIsoNtrksSolid03"  };


  TCanvas* c[nHists];
  TH1F* histos[nHists];
  char cname[100]; 

  for (int i=0; i<nHists; i++) {
    sprintf(cname,"c%i",i);
    histos[i] = (TH1F*)fhists->Get(nameHists[i].Data());
    TString histoTitle = histos[i]->GetTitle();
    histoTitle = histoTitle + " (" + sample + ")";
    histos[i]->SetTitle(histoTitle.Data());
    //    cout << histoTitle << endl;

    if (i==0) { 
      c[i] = new TCanvas(cname, nameHists[i].Data(), 1000, 600); 
      c[i]->cd();
      c[i]->SetBottomMargin(0.4);
      histos[i]->Draw();
      histos[i]->GetXaxis()->SetTitleSize(0.01);
      histos[i]->GetXaxis()->SetBit(TAxis::kLabelsVert);
    } else { 
      c[i] = new TCanvas(cname, nameHists[i].Data(), 600, 600);
      c[i]->cd();
      histos[i]->Draw();
      if ((nameHists[i]=="h_Photon1_hadOverEm")||(nameHists[i]=="h_Photon1_hcalIso04")||(nameHists[i]=="h_Photon1_hcalIso03")||(nameHists[i]=="h_Photon1_trkIsoSumPtHollow04")||(nameHists[i]=="h_Photon1_trkIsoSumPtSolid04")||(nameHists[i]=="h_Photon1_trkIsoSumPtHollow03")||(nameHists[i]=="h_Photon1_trkIsoSumPtSolid03")) { 
	c[i]->SetLogy();   
      }
    }
    
    char fname[100];     
    sprintf(fname,"%s_%s.png",nameHists[i].Data(),sample.Data());  
    if (kPrint) c[i]->Print(fname);
  }

  return;

}
