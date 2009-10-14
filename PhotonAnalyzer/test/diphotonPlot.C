void diphotonPlot(Char_t* sampleType = "zee"){

  gROOT->SetStyle("Plain");

  TFile* f;

  if (sampleType=="zee") {
    f = new TFile("diphoton_zee_hists.root");
  } else if (sampleType=="qcd15") {
    f = new TFile("diphoton_qcd15_hists.root");
  } else if (sampleType=="qcd30") {
    f = new TFile("diphoton_qcd30_hists.root");
  } else {
    cout << "wrong sampleType!" << endl;
    return;
  }

  f->cd("photonAnalyzer");

  c = new TCanvas("c","photon plots",1200,600);
  c->Divide(4,2);
  c->cd(1); 
  fPhoton1Et_h->Draw();
  fPhoton1Et_h->GetXaxis()->SetTitle("#gamma_{1} p_{T} (GeV)");
  c->cd(2); 
  fPhoton2Et_h->Draw();
  fPhoton2Et_h->GetXaxis()->SetTitle("#gamma_{2} p_{T} (GeV)");
  c->cd(3); 
  fPhoton1Eta_h->Draw();
  fPhoton1Eta_h->GetXaxis()->SetTitle("#gamma_{1} #eta");
  c->cd(4); 
  fPhoton2Eta_h->Draw();
  fPhoton2Eta_h->GetXaxis()->SetTitle("#gamma_{2} #eta");
  c->cd(5); 
  fPhoton1Phi_h->Draw();
  fPhoton1Phi_h->GetXaxis()->SetTitle("#gamma_{1} #phi");
  c->cd(6); 
  fPhoton2Phi_h->Draw(); 
  fPhoton2Phi_h->GetXaxis()->SetTitle("#gamma_{2} #phi");
  c->cd(7); 
  fDeltaPhi_h->Draw();
  fDeltaPhi_h->GetXaxis()->SetTitle("#Delta #phi");
  c->cd(8); 
  fInvMass_h->Draw();
  fInvMass_h->GetXaxis()->SetTitle("m_{#gamma #gamma} (Gev)");

  char fileName[100];
  sprintf( fileName,"diphoton_%s.png", sampleType);
  c->Print(fileName);

  return;

}
