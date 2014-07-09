#include "hemJet.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

using namespace std;


float LundDist(LorentzVector& v1, LorentzVector& v2){
  //v1 is the hemisphere axis
  float p2 = sqrt( v2.Px()*v2.Px() + v2.Py()*v2.Py() + v2.Pz()*v2.Pz() );
  return ( (v1.E() - (v1.Vect()).Dot(v2.Vect())/p2) * ( v1.E()/( (v1.E() + v2.E())*(v1.E() + v2.E()) ) ) );
}

std::vector<LorentzVector> getHemJets(std::vector<LorentzVector>& jets){

  if(jets.size() < 2) throw std::runtime_error("less than 2 jets used for hemisphere reconstruction");

  LorentzVector axis1;
  LorentzVector axis2;
  LorentzVector tempAxis1;
  LorentzVector tempAxis2;
  unsigned int axis1Index = 0;
  unsigned int axis2Index = 0;

  std::vector<int> hemAssoc (jets.size(), 0);

  //seed the axes with the pair of jets that form the highest invariant mass
  float maxMass = -1.0;
  for(unsigned int i=0; i<jets.size(); i++) {
    for(unsigned int j=0; j<jets.size(); j++) {

      if(i==j) continue;
      if( (jets.at(i) + jets.at(j)).M() > maxMass ){
        maxMass = (jets.at(i) + jets.at(j)).M();
        axis1Index = i;
        axis2Index = j;
      }

    }
  }

  //set initial axes 
  axis1 = jets.at(axis1Index);
  axis2 = jets.at(axis2Index);

//debug
//cout << "axis1Index = " << axis1Index << endl;
//cout << "axis2Index = " << axis2Index << endl;
  
  //store hemisphere associations for jets that serve as axis seeds
  hemAssoc.at(axis1Index) = 1;
  hemAssoc.at(axis2Index) = 2;

  int counter = 0;
  bool move = true;

  while (move && counter<100) { 

    move = false;
    ++counter;
    if (counter==200) cout << "hemisphere reconstruction not converging" << endl;

    tempAxis1 = LorentzVector(0,0,0,0);
    tempAxis2 = LorentzVector(0,0,0,0);

    //assign each jet to a hemisphere, max of one can be moved per iteration
    for(unsigned int jet=0; jet<jets.size(); jet++){
      if (jet==axis1Index) {
        tempAxis1 += jets.at(jet);
      }
      else if (jet==axis2Index) {
        tempAxis2 += jets.at(jet);
      }
      else {
        if(!move) {
          LorentzVector newAxis1 = hemAssoc.at(jet)==1 ? axis1 - jets.at(jet) : axis1;
          LorentzVector newAxis2 = hemAssoc.at(jet)==2 ? axis2 - jets.at(jet) : axis2;
          //cout << "newAxis1: " << newAxis1.E() << endl;
          //cout << "newAxis2: " << newAxis2.E() << endl;
          //cout << "lund dists: " << LundDist(newAxis1, jets.at(jet)) << ", " << LundDist(newAxis2, jets.at(jet)) << endl;//debug
          if (LundDist(newAxis1, jets.at(jet)) <= LundDist(newAxis2, jets.at(jet))) {
            if (hemAssoc.at(jet) != 1) move = true;
            tempAxis1 += jets.at(jet);
            hemAssoc.at(jet) = 1;
          }
          else {
            if (hemAssoc.at(jet) != 2) move = true;
            tempAxis2 += jets.at(jet);
            hemAssoc.at(jet) = 2;
          }
        }
        else {
          if (hemAssoc.at(jet)==1) tempAxis1 += jets.at(jet);
          else if (hemAssoc.at(jet)==2) tempAxis2 += jets.at(jet);
        }
      }
    }
    
    //new axes for each hemisphere composed of the momenta of the jets in that hemisphere
    axis1 = tempAxis1;
    axis2 = tempAxis2;

  }

  std::vector<LorentzVector> result;
  float E1 = sqrt(axis1.Px()*axis1.Px() + axis1.Py()*axis1.Py() + axis1.Pz()*axis1.Pz()); //used to force hemisphere to be massless
  LorentzVector hem1 = LorentzVector(axis1.Px(), axis1.Py(), axis1.Pz(), E1);
  float E2 = sqrt(axis2.Px()*axis2.Px() + axis2.Py()*axis2.Py() + axis2.Pz()*axis2.Pz()); //used to force hemisphere to be massless
  LorentzVector hem2 = LorentzVector(axis2.Px(), axis2.Py(), axis2.Pz(), E2);
  result.push_back(hem1);
  result.push_back(hem2);

/*
  //for debugging
  cout << endl;
  for(unsigned int i=0; i<hemAssoc.size(); i++){

    cout << i << ", " << hemAssoc.at(i) << endl;

  }
*/

  return(result);  

}
