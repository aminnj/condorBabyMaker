#ifndef UTILITIES_H
#define UTILITIES_H

#include "TMath.h"
 #include "Math/VectorUtil.h"
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

// This is meant to be passed as the third argument, the predicate, of the standard library sort algorithm
inline bool sortByPt(const LorentzVector &vec1, const LorentzVector &vec2 ) {
    return vec1.pt() > vec2.pt();
}

// Delta Phi
inline float deltaPhi( float phi1 , float phi2 ) {
  float dphi = fabs( phi1 - phi2 );
  if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
  return dphi;
}

// Delta R
inline double deltaR(float eta1 , float phi1 , float eta2 , float phi2) {
  double dphi = std::min(::fabs(phi1 - phi2), 2 * M_PI - fabs(phi1 - phi2));
  double deta = eta1 - eta2;
  return sqrt(dphi*dphi + deta*deta);
}


bool   hypsOverlap(int, int );

//double trkIsolation(int trk_index);

#endif
