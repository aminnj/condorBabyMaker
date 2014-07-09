#ifndef triggerUtils_h
#define triggerUtils_h

#include "Math/LorentzVector.h"
#include "Math/Point3D.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

bool goodEGTrigger5July2010 (bool);
int nHLTObjects(const char*);
LorentzVector p4HLTObject(const char*, int) ;
int idHLTObject(const char*, int) ;
void PrintTriggers();

bool passUnprescaledHLTTrigger(const char* arg);
bool passUnprescaledHLTTrigger(const char* arg, const LorentzVector &obj);

int HLT_prescale( const char* arg );
int L1_prescale ( const char* arg );

#endif
