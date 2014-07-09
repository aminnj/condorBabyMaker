///////// 
/// code below based on 
/// http://cmslxr.fnal.gov/lxr/source/PhysicsTools/CandUtils/interface/EventShapeVariables.h
/// still need to add the 
/// http://cmslxr.fnal.gov/lxr/source/PhysicsTools/CandUtils/interface/Thrust.h
////////

#ifndef EventShape_h
#define EventShape_h

//ROOT
#include "TMath.h"
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "Math/Vector3D.h"
#include "Math/Point3D.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"

// CMS2
#include "./CMS2.h"

//#include "Utils.h"

#include "./utilities.h"
#include <vector>

using namespace std; 
//using namespace tas;

typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double> > XYZPoint;
typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> > XYZVector; 

class EventShape {

 public:
  
  //  explicit EventShape(const std::vector<XYZVector>& inputVectors);
  EventShape(const std::vector<LorentzVector> & inputVectors);
  
  //  EventShape();
  // default destructor
  ~EventShape();
  
  double sphericity(double = 2.)  const;
  double aplanarity(double = 2.)  const;
 
  double circularity(const unsigned int& numberOfSteps = 1000) const;

 private:
  
  TMatrixDSym compMomentumTensor(double = 2.) const;
  TVectorD compEigenValues(double = 2.) const;
  
  /// cashing of input vectors
  std::vector<XYZVector> inputVectors_;
  
};

#endif
