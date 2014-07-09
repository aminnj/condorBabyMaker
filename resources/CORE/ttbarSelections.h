#ifndef TTBARSELECTIONS_H
#define TTBARSELECTIONS_H

#include <vector>
#include "Math/LorentzVector.h"


typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;


///******************************************************************************************/     
////return the MET and the MET phi, correcting for mus that are not corrected for by default
///******************************************************************************************/     
//std::pair<float,float> getMet(const string algo, unsigned int hypIdx);

/*****************************************************************************************/
//hypothesis disambiguation. Returns the hypothesis that has the highest sum Pt
/*****************************************************************************************/
unsigned int selectHypByHighestSumPt(const vector<unsigned int> &v_goodHyps);


#endif

