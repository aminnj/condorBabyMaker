//
// taken from http://cmslxr.fnal.gov/lxr/source/PhysicsTools/CandUtils/src/Thrust.cc
//

#ifndef Thrust_h
#define Thrust_h

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

class Thrust {
 public:
  /// spatial vector
  //  typedef math::XYZVector Vector;
  typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double> > XYZPoint;
  typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> > XYZVector;

  //  const std::vector<LorentzVector> & inputVectors;

  /// constructor from first and last iterators
  //  template<typename const_iterator>
  Thrust(const std::vector<LorentzVector> & inputVec) :
    thrust_(0), axis_(0, 0, 0), pSum_(0), 
    n_(inputVec.size()), p_(n_) {
    if (n_ == 0) return;
    
    init(inputVec);
    
  } 

    /// thrust value (in the range [0.5, 1.0])
    double thrust() const { return thrust_; } 
    /// thrust axis (with magnitude = 1)
    const XYZVector& axis() const { return axis_; } 
    
 private:
    double thrust_;
    XYZVector axis_;
    double pSum_;
    const unsigned int n_;
    std::vector<XYZVector> p_;

    struct ThetaPhi {
      ThetaPhi(double t, double p) : theta( t ), phi( p ) { }
      double theta, phi;
    };
    double thrust(const XYZVector & theAxis) const; 
    ThetaPhi initialAxis() const;
    ThetaPhi finalAxis(ThetaPhi) const;
    XYZVector axis(double theta, double phi) const;
    XYZVector axis(const ThetaPhi & tp) const {
      return axis(tp.theta, tp.phi);
    }
    void parabola(double & a, double & b, double & c, 
		  const XYZVector &, const XYZVector &, const XYZVector &) const;
    void init(const std::vector<LorentzVector> & inputVectors);
};

#endif
