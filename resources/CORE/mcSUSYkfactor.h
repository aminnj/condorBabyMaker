#ifndef MCSUSYKFACTOR_H
#define MCSUSYKFACTOR_H
//----------------------------------------------------
// k-factors for LMx samples and SUSY Scan samples -Sanjay.Padhi@cern.ch
//
// Usage:
//       kfactorSUSY("lm0");
//   Or  kfactorSUSY(m0, m12, "tanbeta3");    
//
//--------------------------------------------------------

#include "CMS2.h"
#include "TH2F.h"

float kfactorSUSY(string sample);
float kfactorSUSY(float m0, float m12, string sample);
double lmdata(int ipart1, int ipart2, string prefix);
Float_t GetValueTH2FS(Float_t x, Float_t y, TH2F* h);
int sfinalState(int ipart1, int ipart2);
float cmssm_loxsec(float m0, float m12, string filename="../data/loxsec_m02TeV.root");

#endif

