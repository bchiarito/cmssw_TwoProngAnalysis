#ifndef ReconEfficiency_INC
#define ReconEfficiency_INC

#include <string>

bool passesHoverECut(double HoverE){
  bool resultHoverE = false;
  if (HoverE<0.05){
    resultHoverE = true;
  }
  return resultHoverE;
}

bool EcalIsoCut(double ecalIso04, double PhotonPt,double rho){ 
  bool resultEcalIso = false;
  if (ecalIso04 < (4.2+0.006*PhotonPt+0.183*rho)) {
    resultEcalIso = true;
  }
  return resultEcalIso;
}

bool HcalIsoCut(double hcalIso04,double PhotonPt,double rho){      
  bool resultHcalIso = false;
  if (hcalIso04<(2.2+0.0025*PhotonPt+0.062*rho)){
    resultHcalIso = true;
  }
  return resultHcalIso;
}


bool trckSumPtHollow04Cut(double trckPthollow04, double PhotonPt,double rho){
  bool resultrckSumPt = false;
  if (trckPthollow04<(2+0.001*PhotonPt+0.0167*rho)){
    resultrckSumPt = true;
  }
  return resultrckSumPt; 
}

bool  SigmaIetaIetaCut(double PhotonSigmaIetaIeta){
  bool resultSigmaIetaIeta=false;
  if (PhotonSigmaIetaIeta<0.011){
    resultSigmaIetaIeta=true;
  }
  return resultSigmaIetaIeta;
}


bool PassesPixelSeed(bool PhotonPixelSeed){
  bool result = false;
  if (PhotonPixelSeed){
    result = true;
  }
  return result;
  
}

bool PassesEffCuts(double PhotonPt,bool PixelSeed,double PhotonSigmaIetaIeta,double rho,double  trckPthollow04, double hcalIso04, double ecalIso04,double HoverE){
  bool resultPassesCuts = false;   
  if ( passesHoverECut(HoverE) && EcalIsoCut(ecalIso04,PhotonPt,rho) && HcalIsoCut(hcalIso04,PhotonPt,rho) && SigmaIetaIetaCut(PhotonSigmaIetaIeta) && !PassesPixelSeed(PixelSeed) && trckSumPtHollow04Cut(trckPthollow04,PhotonPt,rho)) {
    resultPassesCuts = true;
  }
  return resultPassesCuts;

}

bool PassEffCutsminusPixel(double PhotonPt,bool PixelSeed,double PhotonSigmaIetaIeta,double rho,double  trckPthollow04, double hcalIso04, double ecalIso04,double HoverE){
  bool resultPassesCuts = false;
  if ( passesHoverECut(HoverE) && EcalIsoCut(ecalIso04,PhotonPt,rho) && HcalIsoCut(hcalIso04,PhotonPt,rho) && SigmaIetaIetaCut(PhotonSigmaIetaIeta) && trckSumPtHollow04Cut(trckPthollow04,PhotonPt,rho)) {
    resultPassesCuts = true;
  }
  return resultPassesCuts;
}

bool PassesEffCutsminusHoverE(double PhotonPt,bool PixelSeed,double PhotonSigmaIetaIeta,double rho,double  trckPthollow04, double hcalIso04, double ecalIso04,double HoverE){
  bool resultPassesCuts = false;
  if ( EcalIsoCut(ecalIso04,PhotonPt,rho) && HcalIsoCut(hcalIso04,PhotonPt,rho) && SigmaIetaIetaCut(PhotonSigmaIetaIeta) && !PassesPixelSeed(PixelSeed) && trckSumPtHollow04Cut(trckPthollow04,PhotonPt,rho)) {
    resultPassesCuts = true;
  }
  return resultPassesCuts;
}

bool PassesEffCutsminusEcal(double PhotonPt,bool PixelSeed,double PhotonSigmaIetaIeta,double rho,double  trckPthollow04, double hcalIso04, double ecalIso04,double HoverE){
  bool resultPassesCuts = false;
  if ( passesHoverECut(HoverE) && HcalIsoCut(hcalIso04,PhotonPt,rho) && SigmaIetaIetaCut(PhotonSigmaIetaIeta) && !PassesPixelSeed(PixelSeed) && trckSumPtHollow04Cut(trckPthollow04,PhotonPt,rho)) {
    resultPassesCuts = true;
  }
  return resultPassesCuts;
}

bool PassesEffCutsminusHcal(double PhotonPt,bool PixelSeed,double PhotonSigmaIetaIeta,double rho,double  trckPthollow04, double hcalIso04, double ecalIso04,double HoverE){
  bool resultPassesCuts = false;
  if ( passesHoverECut(HoverE) && EcalIsoCut(ecalIso04,PhotonPt,rho) && SigmaIetaIetaCut(PhotonSigmaIetaIeta) && !PassesPixelSeed(PixelSeed) &&  trckSumPtHollow04Cut(trckPthollow04,PhotonPt,rho)) {
    resultPassesCuts = true;
  }
  return resultPassesCuts;

}

bool PassesEffCutsminusSigmaIetaIeta(double PhotonPt,bool PixelSeed,double PhotonSigmaIetaIeta,double rho,double  trckPthollow04, double hcalIso04, double ecalIso04,double HoverE){
  bool resultPassesCuts = false;
  if ( passesHoverECut(HoverE) && EcalIsoCut(ecalIso04,PhotonPt,rho) && HcalIsoCut(hcalIso04,PhotonPt,rho) && !PassesPixelSeed(PixelSeed) &&  trckSumPtHollow04Cut(trckPthollow04,PhotonPt,rho)) {
    resultPassesCuts = true;
  }
  return resultPassesCuts;

}

bool PassesEffCutsminstrckIso(double PhotonPt,bool PixelSeed,double PhotonSigmaIetaIeta,double rho,double  trckPthollow04, double hcalIso04, double ecalIso04,double HoverE){
  bool resultPassesCuts = false;
  if ( passesHoverECut(HoverE) && EcalIsoCut(ecalIso04,PhotonPt,rho) && HcalIsoCut(hcalIso04,PhotonPt,rho) && SigmaIetaIetaCut(PhotonSigmaIetaIeta) && !PassesPixelSeed(PixelSeed) ){
    resultPassesCuts = true;
  }
  return resultPassesCuts;

}


#endif
			   
