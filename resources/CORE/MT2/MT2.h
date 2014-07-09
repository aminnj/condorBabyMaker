//---------------------------------------------------------------------------------
// Filename      MT2.cc
// Author:       dbarge
// Created:      February 22 2010
// Modified:     February 09 2011
// Description:  Library to calculate MT2 and MT2J
// References:   http://arxiv.org/PS_cache/hep-ph/pdf/0304/0304226v1.pdf
//               http://arxiv.org/PS_cache/arxiv/pdf/0810/0810.5178v2.pdf
//---------------------------------------------------------------------------------

#ifndef MT2_H
#define MT2_H

#include "MT2Utility.h"
#include "Math/LorentzVector.h"
#include "TMath.h"
#include <vector>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

enum enum_mt2_method { BISECT, GRID };

///////////////////////////////////////////////////////////////
// MT2 Calculated with the Bisection method from Cheng & Han //
///////////////////////////////////////////////////////////////

// MT2 declaration
double MT2(
  const float,
  const float,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
  float = 0.0,
  bool = false
);

// MT2J declaration
double MT2J(
  const float,
  const float,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
  const std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >,
  float = 0.0,
  enum_mt2_method = BISECT,
  bool = false
);

///////////////////////////////////////////////
// MT2 Calculated Using a Simple Grid Method //
///////////////////////////////////////////////

class TMt2 {

  private:
    
    int           grid_size_;
    int           grid_spacing_;
    float         mt2_;
    LorentzVector p4_nu1_;
    LorentzVector p4_nu2_;

  public:

    TMt2    (void);
    ~TMt2   (void);

    float GetMt2  ( const float , const float , const LorentzVector , const LorentzVector ,                               const float = 0.0, bool = false );
    //float GetMt2j ( const float , const float , const LorentzVector , const LorentzVector , const vector<LorentzVector> , const float = 0.0, bool = false );

    inline float         Mt2   (void) { return mt2_;    }
    inline LorentzVector Nu1p4 (void) { return p4_nu1_; }
    inline LorentzVector Nu2p4 (void) { return p4_nu2_; }

};

// Constructor & Destructor
TMt2::TMt2  (void) { 
  grid_size_    = 500;
  grid_spacing_ = 2;
  mt2_          = -999.0;
};
TMt2::~TMt2 (void) { };


#endif

