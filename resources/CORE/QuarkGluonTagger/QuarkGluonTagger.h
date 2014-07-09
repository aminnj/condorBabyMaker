#ifndef QuarksGluonTagger_h
#define QuarksGluonTagger_h

// C++
#include <stdint.h>
#include <vector>

// ROOT
#include "TMath.h"
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TROOT.h"
#include "TFormula.h"

// CMS2
#include "../CMS2.h"
#include "../trackSelections.h"

// Header
#include "QGLikelihoodCalculator.h"

#include "../utilities.h"

using namespace std;
using namespace tas;

//struct indP4_{
//  LorentzVector p4obj;
//  int p4ind;
//};

float getLRM (int ijet , int power);

float constituentPtDistribution(int ijet);

float QGtagger(LorentzVector p4 ,int ijet, QGLikelihoodCalculator * );
 

#endif
