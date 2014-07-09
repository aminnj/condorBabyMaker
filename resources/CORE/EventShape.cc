#include "./EventShape.h"

EventShape::~EventShape(){};

EventShape::EventShape(const std::vector<LorentzVector> & inputVectors) {

  inputVectors_.reserve( inputVectors.size() );
  
  for(unsigned int i_jet = 0; i_jet < inputVectors.size(); i_jet++) { 
    
    std::vector<XYZVector> p;
    
    LorentzVector p4 = inputVectors.at(i_jet);
    
    XYZVector Vjet=XYZVector(p4.px(), p4.py(), p4.pz());
    
    inputVectors_.push_back(Vjet);
    
  }
  
}

/// helper function to fill the 3 dimensional momentum tensor from the inputVecotrs where needed
TMatrixDSym  EventShape::compMomentumTensor(double r) const {
   TMatrixDSym momentumTensor(3);
   momentumTensor.Zero();
 
   if ( inputVectors_.size() < 2 ){
     return momentumTensor;
   }
 
   // fill momentumTensor from inputVectors
   double norm = 1.;
   for ( int i = 0; i < (int)inputVectors_.size(); ++i ){
     double p2 = inputVectors_[i].Dot(inputVectors_[i]);
     double pR = ( r == 2. ) ? p2 : TMath::Power(p2, 0.5*r);
     norm += pR;
     double pRminus2 = ( r == 2. ) ? 1. : TMath::Power(p2, 0.5*r - 1.);
     momentumTensor(0,0) += pRminus2*inputVectors_[i].x()*inputVectors_[i].x();
     momentumTensor(0,1) += pRminus2*inputVectors_[i].x()*inputVectors_[i].y();
     momentumTensor(0,2) += pRminus2*inputVectors_[i].x()*inputVectors_[i].z();
     momentumTensor(1,0) += pRminus2*inputVectors_[i].y()*inputVectors_[i].x();
     momentumTensor(1,1) += pRminus2*inputVectors_[i].y()*inputVectors_[i].y();
     momentumTensor(1,2) += pRminus2*inputVectors_[i].y()*inputVectors_[i].z();
     momentumTensor(2,0) += pRminus2*inputVectors_[i].z()*inputVectors_[i].x();
     momentumTensor(2,1) += pRminus2*inputVectors_[i].z()*inputVectors_[i].y();
     momentumTensor(2,2) += pRminus2*inputVectors_[i].z()*inputVectors_[i].z();
   }
 
   //std::cout << "momentumTensor:" << std::endl;
   //std::cout << momentumTensor(0,0) << " " << momentumTensor(0,1) << " " << momentumTensor(0,2) 
   //          << momentumTensor(1,0) << " " << momentumTensor(1,1) << " " << momentumTensor(1,2) 
   //          << momentumTensor(2,0) << " " << momentumTensor(2,1) << " " << momentumTensor(2,2) << std::endl;
 
   // return momentumTensor normalized to determinant 1
   return (1./norm)*momentumTensor;
}
 
/// helper function to fill the 3 dimensional vector of eigen-values;
/// the largest (smallest) eigen-value is stored at index position 0 (2)
TVectorD EventShape::compEigenValues(double r) const {

   TVectorD eigenValues(3);
   TMatrixDSym myTensor = compMomentumTensor(r);
   if( myTensor.IsSymmetric() ){
     if( myTensor.NonZeros() != 0 ) myTensor.EigenVectors(eigenValues);
   }
 
   // CV: TMatrixDSym::EigenVectors returns eigen-values and eigen-vectors
   //     ordered by descending eigen-values, so no need to do any sorting here...
   // std::cout << "eigenValues(0) = " << eigenValues(0) << ","
   //           << " eigenValues(1) = " << eigenValues(1) << ","
   //           << " eigenValues(2) = " << eigenValues(2) << std::endl;
 
   return eigenValues;
}

double  EventShape::sphericity(double r) const {

  TVectorD eigenValues = compEigenValues(r);
  return 1.5*(eigenValues(1) + eigenValues(2));

}

double  EventShape::circularity(const unsigned int& numberOfSteps) const {

   const double deltaPhi=2*TMath::Pi()/numberOfSteps;
   double circ=-1, phi=0, area = 0;
   for(unsigned int i=0;i<inputVectors_.size();i++) {
     area+=TMath::Sqrt(inputVectors_[i].x()*inputVectors_[i].x()+inputVectors_[i].y()*inputVectors_[i].y());
   }
   for(unsigned int i=0; i<numberOfSteps; ++i){
     phi+=deltaPhi;
     double sum=0, tmp=0.;
     for(unsigned int j=0; j<inputVectors_.size(); ++j){
       sum+=TMath::Abs(TMath::Cos(phi)*inputVectors_[j].x()+TMath::Sin(phi)*inputVectors_[j].y());
     }
     tmp=TMath::Pi()/2*sum/area;
     if( circ<0 || tmp<circ ){
       circ=tmp;
     }
   }
   return circ;
}

double 
EventShape::aplanarity(double r) const
{
   TVectorD eigenValues = compEigenValues(r);
   return 1.5*eigenValues(2);
}
