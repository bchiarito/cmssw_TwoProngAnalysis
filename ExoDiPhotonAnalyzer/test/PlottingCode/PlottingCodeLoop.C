
#define PlottingCodeLoop_cxx
#include "PlottingCodeLoop.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>

using namespace std;

void PlottingCodeLoop::Loop()
{

  if (fChain == 0) return;
 
  cout<<"PlottingCodeLoop::Loop(): "<<fChain->GetEntries()<<" entries in Tree Chain"<<endl;

  Double_t PUweight = 1.;
  Double_t weightFake = 1.;      
  Double_t KFactorweight = 1.; 

  std::vector<Float_t> FRparams;
  //   //Jet dataset(sinin)
  //   FRparams.push_back(0.03957);
  //   FRparams.push_back(3105.);
  //   FRparams.push_back(2.897);
  
  //Jet dataset(sinin with pixel seed veto)
  FRparams.push_back(0.0213);
  FRparams.push_back(344.3);
  FRparams.push_back(2.336);
  
  //     //Photon dataset
  //     FRparams.push_back(0.06529);
  //     FRparams.push_back(2.749E5);
  //     FRparams.push_back(3.917);

  ////Old version of the FR
  //   FRparams.push_back(0.103);
  //   FRparams.push_back(2.597E5);
  //   FRparams.push_back(3.834);

  TString FRformula = TString::Format("(%g+(%g/x^(%g)))",FRparams[0],FRparams[1],FRparams[2]);
  cout<<"Fake rate function used: "<<FRformula.Data()<<endl;
  TF1* fake_rate_fn  = new TF1("fake_rate_fn",FRformula.Data(),10.,2000.) ;

  std::vector<Float_t> KFactorparams;
  ////2011
  //   KFactorparams.push_back(1.07559);
  //   KFactorparams.push_back(1.22200);
  //   KFactorparams.push_back(-0.00185612);

  //2012
  KFactorparams.push_back(0.955490);
  KFactorparams.push_back(1.36034);
  KFactorparams.push_back(-0.0019794);

  TString KFactorformula = TString::Format("%g+(%g*exp(%g*x))",KFactorparams[0],KFactorparams[1],KFactorparams[2]);
  cout<<"K factor function used: "<<KFactorformula.Data()<<endl;
  TF1* KFactorFunction = new TF1("KFactorFunction",KFactorformula.Data(),10.,2000.);   

  Long64_t nentries = fChain->GetEntriesFast();
   
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    if((Photon1_isEB) && (Photon2_isEB) && (Photon1_pt>=80.) && (Photon2_pt>=80.) ){
       
      if(Diphoton_Minv < 200.) continue;

      //set it to true by default
      Bool_t triggerPass = true;
      //Then check the path only if this is data
      if(_SampleType == "data")
      	triggerPass = ( (TrigHLT_HLT_DoublePhoton70_v3 == 1) || (TrigHLT_HLT_DoublePhoton70_v4 == 1) || (TrigHLT_HLT_DoublePhoton70_v5 == 1) || (TrigHLT_HLT_DoublePhoton70_v6 == 1) );      
      if(!triggerPass) continue;
      
      if (_fakeStatus=="TightTight"){
		
	weightFake = 1.;   
	PUweight = 1.;
	KFactorweight = 1.; 

	if (_SampleType == "mc"){                      
	  PUweight = MCPUWeight;
	  // 	  if(Diphoton_Minv < 320.) KFactorweight = 1.73405;
	  // 	  if(Diphoton_Minv >= 320.) KFactorweight = KFactorFunction->Eval(Diphoton_Minv);	  
	  KFactorweight = KFactorFunction->Eval(Diphoton_Minv);
	}//end of if sampletype == mc					

	//cout<<"Pileup weight "<<PUweight<<" K factor "<<KFactorweight<<endl;	
	//cout<<"Filling TT histograms "<<endl;

	h_Photon1_pt->Fill(Photon1_pt,PUweight*KFactorweight);		
	h_Photon1_pt_log->Fill(Photon1_pt,PUweight*KFactorweight);
	h_Photon1_eta->Fill(Photon1_eta,PUweight*KFactorweight);
	h_Photon1_phi->Fill(Photon1_phi,PUweight*KFactorweight);
	h_Photon2_pt->Fill(Photon2_pt,PUweight*KFactorweight);
	h_Photon2_pt_log->Fill(Photon2_pt,PUweight*KFactorweight);
	h_Photon2_phi->Fill(Photon2_phi,PUweight*KFactorweight);
	h_Photon2_eta->Fill(Photon2_eta,PUweight*KFactorweight);
	h_Diphoton_Minv->Fill(Diphoton_Minv,PUweight*KFactorweight);
	h_Diphoton_Minv_log->Fill(Diphoton_Minv,PUweight*KFactorweight);
	h_Diphoton_qt->Fill(Diphoton_qt,PUweight*KFactorweight);
	h_Diphoton_qt_log->Fill(Diphoton_qt,PUweight*KFactorweight);
	h_Diphoton_deltaR->Fill(Diphoton_deltaR,PUweight*KFactorweight);
	h_Diphoton_deltaPhi->Fill(Diphoton_deltaPhi,PUweight*KFactorweight);
	h_Diphoton_deltaEta->Fill(Diphoton_deltaEta,PUweight*KFactorweight);
	h_Diphoton_cosThetaStar->Fill(Diphoton_cosThetaStar,PUweight*KFactorweight);
	h_Vtx_vx->Fill(Vtx_vx,PUweight*KFactorweight);
	h_Vtx_vy->Fill(Vtx_vy,PUweight*KFactorweight);
	h_Vtx_vz->Fill(Vtx_vz,PUweight*KFactorweight);
	h_Vtx_Nvtx->Fill(Vtx_Nvtx,PUweight*KFactorweight);
	h_Photon1_sigmaIetaIeta->Fill(Photon1_sigmaIetaIeta,PUweight*KFactorweight);
	h_Photon1_sigmaEtaEta->Fill(Photon1_sigmaEtaEta,PUweight*KFactorweight);
	h_Photon1_hadOverEm->Fill(Photon1_hadOverEm,PUweight*KFactorweight);
	h_Photon1_hadOverEm_log->Fill(Photon1_hadOverEm,PUweight*KFactorweight);
	h_Photon1_trkIsoSumPtHollow04->Fill(Photon1_trkIsoSumPtHollow04,PUweight*KFactorweight);
	h_Photon1_trkIsoSumPtHollow04_log->Fill(Photon1_trkIsoSumPtHollow04,PUweight*KFactorweight);
	h_Photon1_trkIsoNtrksHollow04->Fill(Photon1_trkIsoNtrksHollow04,PUweight*KFactorweight);
	h_Photon1_hcalIso04->Fill(Photon1_hcalIso04,PUweight*KFactorweight);
	h_Photon1_hcalIso04_log->Fill(Photon1_hcalIso04,PUweight*KFactorweight);
	h_Photon1_ecalIso04->Fill(Photon1_ecalIso04,PUweight*KFactorweight);
	h_Photon1_detEta->Fill(Photon1_detEta,PUweight*KFactorweight);
	h_Photon2_sigmaIetaIeta->Fill(Photon2_sigmaIetaIeta,PUweight*KFactorweight);
	h_Photon2_sigmaEtaEta->Fill(Photon2_sigmaEtaEta,PUweight*KFactorweight);
	h_Photon2_hadOverEm->Fill(Photon2_hadOverEm,PUweight*KFactorweight);
	h_Photon2_hadOverEm_log->Fill(Photon2_hadOverEm,PUweight*KFactorweight);
	h_Photon2_trkIsoSumPtHollow04->Fill(Photon2_trkIsoSumPtHollow04,PUweight*KFactorweight);
	h_Photon2_trkIsoSumPtHollow04_log->Fill(Photon2_trkIsoSumPtHollow04,PUweight*KFactorweight);
	h_Photon2_trkIsoNtrksHollow04->Fill(Photon2_trkIsoNtrksHollow04,PUweight*KFactorweight);
	h_Photon2_hcalIso04->Fill(Photon2_hcalIso04,PUweight*KFactorweight);
	h_Photon2_hcalIso04_log->Fill(Photon2_hcalIso04,PUweight*KFactorweight);
	h_Photon2_ecalIso04->Fill(Photon2_ecalIso04,PUweight*KFactorweight);
	h_Photon2_detEta->Fill(Photon2_detEta,PUweight*KFactorweight);

	//Why is the fake rate applied here ???

	h_FakeRate_tt_pt1->Fill(Photon1_pt,weightFake*PUweight*KFactorweight);
	h_FakeRate_tt_pt1_zoom->Fill(Photon1_pt,weightFake*PUweight*KFactorweight);
	h_FakeRate_tt_eta1->Fill(Photon1_eta,weightFake*PUweight*KFactorweight);
	h_FakeRate_tt_phi1->Fill(Photon1_phi,weightFake*PUweight*KFactorweight);
	h_FakeRate_tt_pt2->Fill(Photon2_pt,weightFake*PUweight*KFactorweight);
	h_FakeRate_tt_pt2_zoom->Fill(Photon2_pt,weightFake*PUweight*KFactorweight);
	h_FakeRate_tt_eta2->Fill(Photon2_eta,weightFake*PUweight*KFactorweight);
	h_FakeRate_tt_phi2->Fill(Photon2_phi,weightFake*PUweight*KFactorweight);
	h_FakeRate_tt_minv->Fill(Diphoton_Minv,weightFake*PUweight*KFactorweight);
	if (Diphoton_Minv>200) h_FakeRate_tt_minv_high->Fill(Diphoton_Minv,weightFake*PUweight*KFactorweight);
	h_FakeRate_tt_qt->Fill(Diphoton_qt,weightFake*PUweight*KFactorweight);
	h_FakeRate_tt_deltaPhi->Fill(Diphoton_deltaPhi,weightFake*PUweight*KFactorweight);
	h_FakeRate_tt_deltaEta->Fill(Diphoton_deltaEta,weightFake*PUweight*KFactorweight);
	h_FakeRate_tt_deltaR->Fill(Diphoton_deltaR,weightFake*PUweight*KFactorweight);
	h_FakeRate_tt_cosThetaStar->Fill(Diphoton_cosThetaStar,weightFake*PUweight*KFactorweight);		
		
      }//end of if status == TightTight

      //TightFake
	      
      if ((abs(Diphoton_deltaPhi)>0.05)){

	if (_fakeStatus=="FakeTight") {  
	  weightFake = fake_rate_fn->Eval(Photon1_pt);

	  h_FakeRate_ft_pt1->Fill(Photon1_pt,weightFake);
	  h_FakeRate_ft_pt1_zoom->Fill(Photon1_pt,weightFake);
	  h_FakeRate_ft_eta1->Fill(Photon1_eta,weightFake);
	  h_FakeRate_ft_phi1->Fill(Photon1_phi,weightFake);
	  h_FakeRate_ft_pt2->Fill(Photon2_pt,weightFake);
	  h_FakeRate_ft_pt2_zoom->Fill(Photon2_pt,weightFake);
	  h_FakeRate_ft_eta2->Fill(Photon2_eta,weightFake);  
	  h_FakeRate_ft_phi2->Fill(Photon2_phi,weightFake);  
	  h_FakeRate_ft_minv->Fill(Diphoton_Minv,weightFake);
	  if (Diphoton_Minv>200) h_FakeRate_ft_minv_high->Fill(Diphoton_Minv,weightFake);  
	  h_FakeRate_ft_qt->Fill(Diphoton_qt,weightFake);
	  h_FakeRate_ft_deltaPhi->Fill(Diphoton_deltaPhi,weightFake);
	  h_FakeRate_ft_deltaEta->Fill(Diphoton_deltaEta,weightFake);
	  h_FakeRate_ft_deltaR->Fill(Diphoton_deltaR,weightFake);
	  h_FakeRate_ft_cosThetaStar->Fill(Diphoton_cosThetaStar,weightFake);
		  
	  h_FakeRate_ft_pt1_noweight->Fill(Photon1_pt);
	  h_FakeRate_ft_pt2_noweight->Fill(Photon2_pt);
	}//end of FT status cut
		
	if (_fakeStatus=="TightFake") {
	  weightFake = fake_rate_fn->Eval(Photon2_pt);

	  h_FakeRate_tf_pt1->Fill(Photon1_pt,weightFake);		  
	  h_FakeRate_tf_pt1_zoom->Fill(Photon1_pt,weightFake) ;		  
	  h_FakeRate_tf_eta1->Fill(Photon1_eta,weightFake);
	  h_FakeRate_tf_phi1->Fill(Photon1_phi,weightFake);
	  h_FakeRate_tf_pt2->Fill(Photon2_pt,weightFake);
	  h_FakeRate_tf_pt2_zoom->Fill(Photon2_pt,weightFake);
	  h_FakeRate_tf_eta2->Fill(Photon2_eta,weightFake);  
	  h_FakeRate_tf_phi2->Fill(Photon2_phi,weightFake);  
	  h_FakeRate_tf_minv->Fill(Diphoton_Minv,weightFake);
	  if (Diphoton_Minv>200) h_FakeRate_tf_minv_high->Fill(Diphoton_Minv,weightFake);  
	  h_FakeRate_tf_qt->Fill(Diphoton_qt,weightFake);
	  h_FakeRate_tf_deltaPhi->Fill(Diphoton_deltaPhi,weightFake);
	  h_FakeRate_tf_deltaEta->Fill(Diphoton_deltaEta,weightFake);
	  h_FakeRate_tf_deltaR->Fill(Diphoton_deltaR,weightFake);
	  h_FakeRate_tf_cosThetaStar->Fill(Diphoton_cosThetaStar,weightFake);
		  
	  h_FakeRate_tf_pt1_noweight->Fill(Photon1_pt);
	  h_FakeRate_tf_pt2_noweight->Fill(Photon2_pt);
	}//end of TF status cut
		
	if (_fakeStatus=="FakeFake") {
	  weightFake = (fake_rate_fn->Eval(Photon1_pt))*(fake_rate_fn->Eval(Photon2_pt));

	  h_FakeRate_ff_pt1->Fill(Photon1_pt,weightFake);
	  h_FakeRate_ff_pt1_zoom->Fill(Photon1_pt,weightFake);
	  h_FakeRate_ff_eta1->Fill(Photon1_eta,weightFake);
	  h_FakeRate_ff_phi1->Fill(Photon1_phi,weightFake);
	  h_FakeRate_ff_pt2->Fill(Photon2_pt,weightFake);
	  h_FakeRate_ff_pt2_zoom->Fill(Photon2_pt,weightFake);
	  h_FakeRate_ff_eta2->Fill(Photon2_eta,weightFake);  
	  h_FakeRate_ff_phi2->Fill(Photon2_phi,weightFake);  
	  h_FakeRate_ff_minv->Fill(Diphoton_Minv,weightFake);
	  if (Diphoton_Minv>200) h_FakeRate_ff_minv_high->Fill(Diphoton_Minv,weightFake);  
	  h_FakeRate_ff_qt->Fill(Diphoton_qt,weightFake);
	  h_FakeRate_ff_deltaPhi->Fill(Diphoton_deltaPhi,weightFake);
	  h_FakeRate_ff_deltaEta->Fill(Diphoton_deltaEta,weightFake);
	  h_FakeRate_ff_deltaR->Fill(Diphoton_deltaR,weightFake);
	  h_FakeRate_ff_cosThetaStar->Fill(Diphoton_cosThetaStar,weightFake);
		  
	  h_FakeRate_ff_pt1_noweight->Fill(Photon1_pt);
	  h_FakeRate_ff_pt2_noweight->Fill(Photon2_pt);
	}//end of FF status cut

      }//end of deltaphi cut

    }//end of barrel and pt cut
     
  }//end of loop over entries    
   
  _outputfile->cd();

  if (_fakeStatus=="TightTight"){
    h_Photon1_pt->Write();
    h_Photon1_pt_log->Write();
    h_Photon1_eta->Write();
    h_Photon1_phi->Write();
    h_Photon2_pt->Write();
    h_Photon2_pt_log->Write();
    h_Photon2_phi->Write();
    h_Photon2_eta->Write();
    h_Diphoton_Minv->Write();
    h_Diphoton_Minv_log->Write();
    h_Diphoton_qt->Write();
    h_Diphoton_qt_log->Write();
    h_Diphoton_deltaR->Write();
    h_Diphoton_deltaPhi->Write();
    h_Diphoton_deltaEta->Write();
    h_Diphoton_cosThetaStar->Write(); 
    h_Vtx_vx->Write();
    h_Vtx_Nvtx->Write();
    h_Vtx_vy->Write();
    h_Vtx_vz->Write();
    h_Photon1_sigmaIetaIeta->Write();
    h_Photon1_sigmaEtaEta->Write();
    h_Photon1_hadOverEm->Write();
    h_Photon1_hadOverEm_log->Write();
    h_Photon1_trkIsoSumPtHollow04->Write();
    h_Photon1_trkIsoSumPtHollow04_log->Write();
    h_Photon1_trkIsoNtrksHollow04->Write();
    h_Photon1_hcalIso04->Write();
    h_Photon1_hcalIso04_log->Write();
    h_Photon1_ecalIso04->Write();
    h_Photon1_detEta->Write();
    h_Photon2_sigmaIetaIeta->Write();
    h_Photon2_sigmaEtaEta->Write();
    h_Photon2_hadOverEm->Write();
    h_Photon2_hadOverEm_log->Write();
    h_Photon2_trkIsoSumPtHollow04->Write();
    h_Photon2_trkIsoSumPtHollow04_log->Write();
    h_Photon2_trkIsoNtrksHollow04->Write();
    h_Photon2_hcalIso04->Write();
    h_Photon2_hcalIso04_log->Write();  
    h_Photon2_ecalIso04->Write();
    h_Photon2_detEta->Write();
    h_FakeRate_tt_pt1->Write();
    h_FakeRate_tt_pt1_zoom->Write();
    h_FakeRate_tt_eta1->Write();
    h_FakeRate_tt_phi1->Write();
    h_FakeRate_tt_pt2->Write();
    h_FakeRate_tt_pt2_zoom->Write();
    h_FakeRate_tt_eta2->Write();
    h_FakeRate_tt_phi2->Write();
    h_FakeRate_tt_minv->Write();   
    h_FakeRate_tt_minv_high->Write();   
    h_FakeRate_tt_qt->Write();
    h_FakeRate_tt_deltaPhi->Write();
    h_FakeRate_tt_deltaEta->Write();
    h_FakeRate_tt_deltaR->Write();
    h_FakeRate_tt_cosThetaStar->Write();
  }

  if (_fakeStatus=="TightFake"){
    h_FakeRate_tf_pt1->Write();
    h_FakeRate_tf_pt1_zoom->Write();
    h_FakeRate_tf_eta1->Write();
    h_FakeRate_tf_phi1->Write();
    h_FakeRate_tf_pt2->Write();
    h_FakeRate_tf_pt2_zoom->Write();
    h_FakeRate_tf_eta2->Write();
    h_FakeRate_tf_phi2->Write();
    h_FakeRate_tf_minv->Write();
    h_FakeRate_tf_minv_high->Write();
    h_FakeRate_tf_qt->Write();
    h_FakeRate_tf_deltaPhi->Write();
    h_FakeRate_tf_deltaEta->Write();
    h_FakeRate_tf_deltaR->Write();
    h_FakeRate_tf_cosThetaStar->Write();
    h_FakeRate_tf_pt1_noweight->Write();
    h_FakeRate_tf_pt2_noweight->Write();
  }

  if (_fakeStatus=="FakeTight"){
    h_FakeRate_ft_pt1->Write();
    h_FakeRate_ft_pt1_zoom->Write();
    h_FakeRate_ft_eta1->Write();
    h_FakeRate_ft_phi1->Write();
    h_FakeRate_ft_pt2->Write();
    h_FakeRate_ft_pt2_zoom->Write();
    h_FakeRate_ft_eta2->Write();
    h_FakeRate_ft_phi2->Write();
    h_FakeRate_ft_minv->Write();
    h_FakeRate_ft_minv_high->Write();
    h_FakeRate_ft_qt->Write();
    h_FakeRate_ft_deltaPhi->Write();
    h_FakeRate_ft_deltaEta->Write();
    h_FakeRate_ft_deltaR->Write();
    h_FakeRate_ft_cosThetaStar->Write();
    h_FakeRate_ft_pt1_noweight->Write();
    h_FakeRate_ft_pt2_noweight->Write();
  }

  if (_fakeStatus=="FakeFake"){
    h_FakeRate_ff_pt1->Write();
    h_FakeRate_ff_pt1_zoom->Write();
    h_FakeRate_ff_eta1->Write();
    h_FakeRate_ff_phi1->Write();
    h_FakeRate_ff_pt2->Write();
    h_FakeRate_ff_pt2_zoom->Write();
    h_FakeRate_ff_eta2->Write();
    h_FakeRate_ff_phi2->Write();
    h_FakeRate_ff_minv->Write();
    h_FakeRate_ff_minv_high->Write();
    h_FakeRate_ff_qt->Write();
    h_FakeRate_ff_deltaPhi->Write();
    h_FakeRate_ff_deltaEta->Write();
    h_FakeRate_ff_deltaR->Write();
    h_FakeRate_ff_cosThetaStar->Write();
    h_FakeRate_ff_pt1_noweight->Write();
    h_FakeRate_ff_pt2_noweight->Write();
  }
  _outputfile->cd();
  _outputfile->Close();   

}//end of Loop method
