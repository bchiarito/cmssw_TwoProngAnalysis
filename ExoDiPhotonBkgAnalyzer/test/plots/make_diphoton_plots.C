
void make_diphoton_plots(TString sample = "PhotonJet_Pt30to50") {

  TString infile = "rfio:/castor/cern.ch/user/t/torimoto/physics/diphoton/ntuples/mc/"+sample+"/diphotonTree_"+sample+".root";
  
  //  if (!infile) {
  //    cout << " No input file specified !" << endl;
  //    return;
  //  }

  cout << "Making quick check plots for: " << infile << endl;

  gSystem->Load("fTree_C.so");

  TFile* f = TFile::Open(infile.Data());

  TTree* fChain = (TTree*)f->Get("diphotonBkgAnalyzer/fTree");

  fTree* reader = new fTree(fChain);

  reader->Loop();
 
  gROOT->SetStyle("Plain");

  gStyle->SetOptStat(10);

  const int nHists=9;
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

  c[0]->cd();
  c0->SetBottomMargin(0.4);
  h_TrigHLT->GetXaxis()->SetTitleSize(0.01);
  h_TrigHLT->GetXaxis()->SetBit(TAxis::kLabelsVert);
  h_TrigHLT->Draw();
  sprintf(fname,"trigger_%s.png",sample.Data());  
  c[0]->Print(fname);

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
  c[1]->Print(fname);

  c[2]->cd(1);
  h_Photon1_swisscross->Draw();
  c[2]->cd(2);
  h_Photon1_severityLevel->Draw();
  c[2]->cd(3);
  h_Photon1_recHitFlag->Draw();
  c[2]->cd(4);
  h_Photon1_maxRecHitTime->Draw();
  sprintf(fname,"spike_photon1_%s.png",sample.Data());  
  c[2]->Print(fname);

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
  c[3]->Print(fname);

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
  c[4]->Print(fname);

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
  c[5]->Print(fname);

  c[6]->cd(1);
  h_Photon2_swisscross->Draw();
  c[6]->cd(2);
  h_Photon2_severityLevel->Draw();
  c[6]->cd(3);
  h_Photon2_recHitFlag->Draw();
  c[6]->cd(4);
  h_Photon2_maxRecHitTime->Draw();
  sprintf(fname,"spike_photon2_%s.png",sample.Data());  
  c[6]->Print(fname);

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
  c[7]->Print(fname);


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
  c[8]->Print(fname);
  
  //  delete reader;
  return;

}
