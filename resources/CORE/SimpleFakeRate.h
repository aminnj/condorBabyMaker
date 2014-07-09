#ifndef SimpleFakeRate_H
#define SimpleFakeRate_H
//----------------------------------------------------
// Code to read a fake rate from a 2D histogram.
// Eta must be on the X axis, pt on the Y axis
//
// Usage:
// -- Initialize once (assumes that last pt bin is overflow)
// -- Do this for each fake rate you want to use (V1, V2, ...)  
// SimpleFakeRate fr("file.root","histname");
//
// -- same as above but takes FR=0 if pt > max in histogram;
// SimpleFakeRate fr("file.root","histname", false);
//
// -- get fake rate and error for given pt, eta
// float thisFrValue = fr.getFR(pt, eta);
// float thisFrError = fr.getFRerr(pt, eta);
//
//--------------------------------------------------------
//
#include "TSystem.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TAxis.h"
#include <iostream>

using namespace std;

class SimpleFakeRate {
 public:
  // constructor
  SimpleFakeRate(const char* filename, const char* histname , 
                 bool lastBinIsOverflow=true);

  // destructor
  ~SimpleFakeRate() {};

  // return fake rate and error
  float getFR( float pt, float eta) {return getValue(pt,eta,1);};
  float getFRerr( float pt, float eta) {return getValue(pt,eta,2);};

  // a dump routine
  void print();

  // get the histograms (you should not need it ever)
  TH2F* getHist() {return hist_;};

 private:

  float getValue(float pt, float eta, int iflag);
  TH2F* hist_;
  bool overflowFlag_;
  float ptmax_;
  float ptmin_;
  float etamax_;
  float etamin_;
  bool useAbsEta_;
};
#endif
