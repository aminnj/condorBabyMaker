////////////////////////////////////////////////////////////////////////////////
//
// JetResolution_t
// ---------------
//
//            11/05/2010 Philipp Schieferdecker <philipp.schieferdecker@cern.ch>
////////////////////////////////////////////////////////////////////////////////


#include "JetResolution.h"


#include <TROOT.h>
#include <TApplication.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TMath.h>

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>


using namespace std;


//______________________________________________________________________________
int main(int argc,char**argv)
{
  if (argc>1&&string(argv[1])=="--help") {
    cout<<"USAGE: JetResolution_t --era <era> --alg <alg> --nevts <n> --gaussian"<<endl;
    return 0;
  }
  
  string   era("Spring10");
//   string   alg("AK5Calo");
//   string   alg("AK5JPT");
  string   alg("AK5PF");
  unsigned nevts(100);
  bool     doGaussian(false);
  
  for (int i=1;i<argc;i++) {
    if      (string(argv[i])=="--era")      { era = argv[i+1]; i++; }
    else if (string(argv[i])=="--alg")      { alg = argv[i+1]; i++; }
    else if (string(argv[i])=="--nevts")    { stringstream ss; ss<<argv[i+1]; ss>>nevts; i++; }
    else if (string(argv[i])=="--gaussian") { doGaussian = true; }
    else {
      cerr<<"ERROR: unknown option '"<<argv[i]<<"'"<<endl;
      return 0;
    }
  }
  
  cout<<"era:      "<<era<<endl;
  cout<<"alg:      "<<alg<<endl;
  cout<<"nevts:    "<<nevts<<endl;
  cout<<"gaussian: "<<doGaussian<<endl<<endl;
  
  string path = "txts";
  
  string ptFileName  = path + "/" + era + "_PtResolution_" +alg+".txt";
  string etaFileName = path + "/" + era + "_EtaResolution_"+alg+".txt";
  string phiFileName = path + "/" + era + "_PhiResolution_"+alg+".txt";

  cout<<ptFileName<<endl;
  cout<<etaFileName<<endl;
  cout<<phiFileName<<endl;
  cout<<endl;

  JetResolution ptResol(ptFileName,doGaussian);
  JetResolution etaResol(etaFileName,doGaussian);
  JetResolution phiResol(phiFileName,doGaussian);
  
  // SIMPLE TEST
  float pt =  70.0;
  float eta = 3.0;
  float phi = 3.0;

  cout<<"pT="<<pt<<" eta="<<eta<<" phi="<<phi<<endl;

  TF1* fPtResol  = ptResol.resolutionEtaPt(eta,pt);
  TF1* fEtaResol = etaResol.resolutionEtaPt(eta,pt);
  TF1* fPhiResol = phiResol.resolutionEtaPt(eta,pt);

  fPtResol ->Print(); cout<<endl;
  fEtaResol->Print(); cout<<endl;
  fPhiResol->Print(); cout<<endl;

  for (int i = 0; i < 20; i++){

    float rndpt   = fPtResol ->GetRandom();
    float rndeta  = fEtaResol->GetRandom();
    float rndphi  = fPhiResol->GetRandom();
   
    float jetpt   = rndpt*pt;
    float jeteta  = rndeta+eta;
    float jetphi  = rndphi+phi;
    
    cout<<"pT="<<jetpt<<" eta="<<jeteta<<" phi="<<jetphi<<endl;       
  }
 
  return 0;
}

