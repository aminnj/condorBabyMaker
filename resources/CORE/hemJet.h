#include <vector>
#include "Math/LorentzVector.h"
#include "Math/Vector3D.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

bool sortByMass(LorentzVector, LorentzVector);
float LundDist(LorentzVector&, LorentzVector&);
std::vector<LorentzVector> getHemJets(std::vector<LorentzVector>&);
