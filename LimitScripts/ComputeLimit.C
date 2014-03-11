#include <string>
#include <stdlib.h>
#include <sstream>
#include <fstream>

#include <TROOT.h>


std::string floatToString(float num)
{
  using namespace std;
  ostringstream myStream;
  myStream << num << flush;
  return(myStream.str()); //returns the string form of the stringstream object
}


std::string intToString(int num)
{
  using namespace std;
  ostringstream myStream;
  myStream << num << flush;
  return(myStream.str()); //returns the string form of the stringstream object
}


void ComputeLimit(float lumi, float lumiError, float totalEff, float totalEffErrStat, float totalEffErrSyst,
    float totalEffMScaleSystUp, float totalEffMScaleSystDown, float totalEffMResSystUp, float totalEffMResSystDown,
    float totalEffPileupSystUp, float totalEffPileupSystDown,
    float nBackground, float nBackgroundErrStat, float nBackgroundErrSyst,
    float nDataObs, int mass, float coupling, float halfWidth,
    float totalXSec, int massWindowLow, int massWindowHigh, std::string fileName)
{
  using namespace std;
  
  float totalEffErr = sqrt(pow(totalEffErrStat,2)+pow(totalEffErrSyst,2));
  float nBackgroundErr = sqrt(pow(nBackgroundErrStat,2)+pow(nBackgroundErrSyst,2));

  cout << "Compute limit for coupling= " << coupling << " mass= " << mass << endl;
  cout << "Efficiency*Acceptance= " << totalEff << " +/- " << totalEffErr << " nBackground= " << nBackground << " +/- " << nBackgroundErr << endl;
  cout << "nDataObs= " << nDataObs << " halfWidth= " << halfWidth << " totalXSec= " << totalXSec << endl;

  string cmsswBase;
  cmsswBase = string(getenv("CMSSW_BASE"));
  string cl95MacroPath = cmsswBase;
  cl95MacroPath+="/src/StatisticalTools/RooStatsRoutines/root/roostats_cl95.C";

  string roofitsys;
  roofitsys = string(getenv("ROOFITSYS"));
  gROOT->Reset();
  string includePath = "-I"+roofitsys;
  gSystem->SetIncludePath((includePath+"/include").c_str());
  string command = ".L ";
  command+=cl95MacroPath;
  command+="+g";
  gROOT->ProcessLine(command.c_str());
  SetParameter("Optimize",false);
  SetParameter("NToys",10000);
  SetParameter("NClsSteps",25);
  SetParameter("MakePlot",true);
  SetParameter("WriteResult",true);
  SetParameter("PlotHypoTestResult",true);
  
  // new interface
  //limit = ROOT.GetClsLimit(lumi, lumiErr, modelPoint.GetTotalEff(), 
  //                         math.sqrt((modelPoint.GetTotalEff()*(1.-modelPoint.GetTotalEff()))/25000.),
  //                         modelPoint.GetNBackground(), modelPoint.GetNBackgroundErr(),
  //                         modelPoint.GetNDataObs());
  
  // legacy interface
  string plotFileName = "plots_cl95_";
  plotFileName+=floatToString(coupling);
  plotFileName+="_";
  plotFileName+=intToString(mass);
  plotFileName+=+".pdf";
  LimitResult limit = roostats_limit(lumi, lumiError, totalEff,
                                     totalEffErr,
                                     nBackground, nBackgroundErr,
                                     nDataObs, false, 0, "cls",
                                     plotFileName.c_str(),0);
                                     //"",123456); // seed for testing

 // open file
 ofstream myfile(fileName.c_str());
 if(myfile.is_open())
 {
   // write file
   myfile << "Coupling: " << coupling << "\n";
   myfile << "Mass: " << mass << "\n";
   myfile << "TotalXSection: " << totalXSec << "\n";
   myfile << "TotalEff: " << totalEff << "\n";
   myfile << "TotalEffErrStat: " << totalEffErrStat << "\n";
   myfile << "TotalEffErrSyst: " << totalEffErrSyst << "\n";
   myfile << "TotalEffMScaleSystUp: " << totalEffMScaleSystUp << "\n";
   myfile << "TotalEffMScaleSystDown: " << totalEffMScaleSystDown << "\n";
   myfile << "TotalEffMResSystUp: " << totalEffMResSystUp << "\n";
   myfile << "TotalEffMResSystDown: " << totalEffMResSystDown << "\n";
   myfile << "TotalEffPileupSystUp: " << totalEffPileupSystUp << "\n";
   myfile << "TotalEffPileupSystDown: " << totalEffPileupSystDown << "\n";
   myfile << "HalfWidth: " << halfWidth << "\n";
   myfile << "OptMassWindowLow: " << massWindowLow << "\n";
   myfile << "OptMassWindowHigh: " << massWindowHigh << "\n";
   myfile << "NDataObs: " << nDataObs << "\n";
   myfile << "NBackground: " << nBackground << "\n";
   myfile << "NBackgroundErrStat: " << nBackgroundErrStat << "\n";
   myfile << "NBackgroundErrSyst: " << nBackgroundErrSyst << "\n";
   myfile << "ExpectedLimit: " << limit.GetExpectedLimit() << "\n";
   myfile << "ExpectedLimitOneSigmaHigh: " << limit.GetOneSigmaHighRange() << "\n";
   myfile << "ExpectedLimitOneSigmaLow: " << limit.GetOneSigmaLowRange() << "\n";
   myfile << "ExpectedLimitTwoSigmaHigh: " << limit.GetTwoSigmaHighRange() << "\n";
   myfile << "ExpectedLimitTwoSigmaLow: " << limit.GetTwoSigmaLowRange() << "\n";
   myfile << "ObservedLimit: " << limit.GetObservedLimit() << "\n\n";
   myfile << flush;
 }
 myfile.close();

 gApplication->Terminate();

}



