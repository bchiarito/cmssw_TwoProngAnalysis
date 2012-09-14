#ifndef AcceptanceCuts_INC
#define AcceptanceCuts_INC

#include <string>

bool passesPhotonPtCut(double photonpt){
  bool result = false;
  if (photonpt>70){
    result = true;
  }
  return result;
}

bool inBarrell(double PhotonDetEta){ 
  bool inEB = false;
  if ( fabs(PhotonDetEta)<=1.442) {
    inEB= true;
  }
  return inEB;
}

bool MinvCut(double Minv){
  bool resultMinv = false;
  if (Minv >140) {
    resultMinv = true;
  }
  return resultMinv;
}

bool passesAcptCuts(double Minv, double PhotonDetEta, double Photonpt){
  bool resultpassesAcptCuts = false;
  if (MinvCut(Minv) && inBarrell(PhotonDetEta) && passesPhotonPtCut(Photonpt)){
    resultpassesAcptCuts = true;
  }
  return resultpassesAcptCuts;
}
#endif
