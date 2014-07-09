#include "MT2.h"
#include "math.h"
#include "Math/VectorUtil.h"
#include "iostream"

using namespace std;

///////////////////////////////////////////////////////////////
// MT2 Calculated with the Bisection method from Cheng & Han //
///////////////////////////////////////////////////////////////

// MT2( MET_MAGNITUDE, MET_PHI, P4_LEPTON_1, P4_LEPTON_2, MASS_INVISIBLE_PARTICLE )
double MT2(
  const float met,
  const float metPhi,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > v1,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > v2,
  float invisible_particle_mass,
  bool verbose
){

//--- code here follows documentation in MT2Utility.cc ---//

  // initialize arrays for lepton 1, lepton 2, MET
  double pa[3]    = {0,0,0};
  double pb[3]    = {0,0,0};
  double pmiss[3] = {0,0,0};

  // Set the masses
  if( v1.M2() >= 0 ){
    pa[0] = v1.M();
  } else {
    pa[0] = 0.0;
    if(verbose) cout << "ERROR: v1.M2 < 0 ... Setting v1.M() = 0" << endl;
  }
  if( v2.M2() >= 0 ){
    pb[0] = v2.M();
  } else {
    pb[0] = 0.0;
    if(verbose) cout << "ERROR: v2.M2 < 0 ... Setting v2.M() = 0" << endl;
  }

  // set the transverse momenta for the leptons & MET
  pa[1]     = (double) v1.Px();
  pa[2]     = (double) v1.Py();
  pb[1]     = (double) v2.Px();
  pb[2]     = (double) v2.Py();
  pmiss[0]  = 0.0;                      // not used
  pmiss[1]  = (double) met*cos(metPhi);
  pmiss[2]  = (double) met*sin(metPhi);

  // instantiate mt2 class, set momenta and mass of invisible particle
  mt2_bisect::mt2 mt2_event;
  mt2_event.set_momenta( pa, pb, pmiss );
  mt2_event.set_mn( invisible_particle_mass );

  //
  return mt2_event.get_mt2();
}

// MT2J( MET_MAGNITUDE, MET_PHI, P4_LEPTON_1, P4_LEPTON_2, VECT_P4_Jets, MASS_INVISIBLE_PARTICLE, MT2_CALCULATION_METHOD )
double MT2J(
  const float met,
  const float metPhi,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > p4_lepton_1,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > p4_lepton_2,
  const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > vect_p4_jets,
  float,
  enum enum_mt2_method method_mt2,
  bool
){
  if( vect_p4_jets.size() < 2 ){
    cout << "MT2.cc: error MT2J called with less than 2 jets... returning mt2 value of -1.0" << endl;
    return -1.0; 
  }
  double mt2_min = std::numeric_limits<double>::max();
  for(unsigned int j1=0; j1 < vect_p4_jets.size(); j1++){
  for(unsigned int j2=0; j2 < vect_p4_jets.size(); j2++){
    if(j1==j2) continue;
    double mt2_j1_j2;
    if( method_mt2 == BISECT){
      mt2_j1_j2 = MT2( met, metPhi, p4_lepton_1 + vect_p4_jets.at(j1), p4_lepton_2 + vect_p4_jets.at(j2) );
    } else {
      cout << "ERROR: Undefined calculation method... returning -1" << endl;
      return -1.0;
    }
    if(mt2_j1_j2 < mt2_min){
      mt2_min = mt2_j1_j2;
    }
  }}
  return mt2_min;
}

//////////
// GRID //
//////////

//float TMt2::GetMt2j (
//  const float                 met,
//  const float                 metPhi,
//  const LorentzVector         p4_lep1,
//  const LorentzVector         p4_lep2,
//  const vector<LorentzVector> v_p4_jets,
//  const float                 invisible_particle_mass,
//  const bool                  verbose
//){
//  mt2_ = -999.0;
//  return mt2_;
//}

float TMt2::GetMt2 (
  const float         met,
  const float         metPhi,
  const LorentzVector p4_lep1,
  const LorentzVector p4_lep2,
  const float         invisible_particle_mass,
  const bool          verbose
){

  //
  mt2_ = -999.0;

  // Check Input
  if( invisible_particle_mass != 0 ){
    cout << "ERROR: Non zero invisible particle masses are not implemented... Exiting." << endl;
    exit(1);
  }
  if( grid_spacing_ > grid_size_ ){
    cout << endl;
    cout << "ERROR: Grid size must be larger than grid spacing... Exiting." << endl;
    cout << endl;
    exit(1);
  }

  // Particle masses... M^2 < 0 sometimes due to floating poing imprecision, leptons are in the massless limit anyway, set M = 0 when M^2 < 0
  double mass1;
  double mass2;
  if( p4_lep1.M2() >= 0 ){
    mass1 = p4_lep1.M();
  }
  else {
    mass1 = 0.0;
    if(verbose) cout << "p4_lep1.M() < 0... Setting p4_lep1.M() = 0" << endl;
  }
  if( p4_lep2.M2() >= 0 ){
    mass2 = p4_lep2.M();
  }
  else{
    mass2 = 0.0;
    if(verbose) cout << "p4_lep2.M() < 0... Setting p4_lep2.M() = 0" << endl;
  }

  // Set transverse quantities 
  double Et1    = p4_lep1.Et();
  double Et2    = p4_lep2.Et();
  double px1    = p4_lep1.Px();
  double py1    = p4_lep1.Py();
  double px2    = p4_lep2.Px();
  double py2    = p4_lep2.Py();
  double metx   = met*cos(metPhi);
  double mety   = met*sin(metPhi);

  // Initialize grid variables
  double EtNu1 = 0;
  double EtNu2 = 0;
  double mtsq1 = 0;
  double mtsq2 = 0;
  double nu1Px = 0;
  double nu1Py = 0;
  double nu2Px = 0;
  double nu2Py = 0;
  double Max   = 0;
  double Min   = std::numeric_limits<double>::max();

  // Grid Minimization
  for(int x = ( -1*grid_size_ ); x <= grid_size_; x += grid_spacing_ ){
  for(int y = ( -1*grid_size_ ); y <= grid_size_; y += grid_spacing_ ){
      double trial1_Px = x;
      double trial1_Py = y;
      double trial2_Px = ( metx - x );
      double trial2_Py = ( mety - y );
      EtNu1            = sqrt( pow( trial1_Px, 2 ) + pow( trial1_Py, 2 ) );
      EtNu2            = sqrt( pow( trial2_Px, 2 ) + pow( trial2_Py, 2 ) );
      mtsq1            = pow( mass1, 2 ) + 2 * ( Et1 * EtNu1 - px1 * trial1_Px - py1 * trial1_Py );
      mtsq2            = pow( mass2, 2 ) + 2 * ( Et2 * EtNu2 - px2 * trial2_Px - py2 * trial2_Py );
      Max              = max( sqrt( mtsq1 ), sqrt( mtsq2 ) );
      if ( Max < Min ){
        Min = Max;
        nu1Px = trial1_Px;
        nu1Py = trial1_Py;
        nu2Px = trial2_Px;
        nu2Py = trial2_Py;
      }
  } }
  if( Min == std::numeric_limits<double>::max() ){
    cout << "MT2_GRID: Error could not find minimum." << endl;
    exit(1);
  }

  // Get the ( transverse ) 4-vectors of the neutrinos chosen by MT2
  p4_nu1_.SetPxPyPzE( nu1Px, nu1Py, 0.0, sqrt( pow( nu1Px, 2) + pow( nu1Py, 2) ) );
  p4_nu2_.SetPxPyPzE( nu2Px, nu2Py, 0.0, sqrt( pow( nu2Px, 2) + pow( nu2Py, 2) ) );

  //
  mt2_ = Min;
  return mt2_;

}
