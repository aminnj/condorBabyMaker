//--------------------------------------------------
// Trigger utilities stolen from Derek and adapted
//-------------------------------------------------

#include "TSystem.h"
#include "triggerUtils.h"
#include "CMS2.h"
#include "Math/VectorUtil.h"

//--------------------------------------------------
// EG trigger selection from 5 July 2010
// data: Photon10 OR Electron 10 OR Photon 15
// MC  : E10 OR Photon15
// (Note: in data Photon10 is prescaled at some point)
// (Note: in MC Photon15 has a different L1 requirement, so we 
//        use Photon10 and we tighten the pt threshold)
//--------------------------------------------------
bool goodEGTrigger5July2010 (bool mc) {
  
  if (mc) {
    int e10 = nHLTObjects("HLT_Ele10_LW_L1R");
    if (e10 != 0) return true;

    int p10 = nHLTObjects("HLT_Photon10_L1R");
    for (int i=0; i<p10; i++) {
      LorentzVector v = p4HLTObject("HLT_Photon10_L1R", i);
      if (v.pt() > 15.) return true;
    }
 
  } else {  // data now
    int e10 = nHLTObjects("HLT_Ele10_LW_L1R");
    if (e10 != 0) return true;

    int p10 = nHLTObjects("HLT_Photon10_L1R");
    if (p10 != 0) return true;

    p10 = nHLTObjects("HLT_Photon10_Cleaned_L1R");
    if (p10 != 0) return true;

    int p15 = nHLTObjects("HLT_Photon15_L1R");
    if (p15 != 0) return true;

    p15 = nHLTObjects("HLT_Photon15_Cleaned_L1R");
    if (p15 != 0) return true;

 
  }
  return false;
}



//-----------------------------------------------------
// Returns the nth object that passes a given trigger
// (n starts from 0)
///----------------------------------------------------
LorentzVector p4HLTObject(const char* arg, int objNumber){
 
  TString HLTTrigger( arg );
  int trigIndx = -1;
  vector<TString>::const_iterator begin_it = cms2.hlt_trigNames().begin();
  vector<TString>::const_iterator end_it = cms2.hlt_trigNames().end();
  vector<TString>::const_iterator found_it = find(begin_it, end_it, HLTTrigger );
  if( (found_it != end_it) ){
    trigIndx = found_it - begin_it;
    //cout << "p4HLTObject: Found Trigger: " << arg << endl;
  }
  else {
    cout << "p4HLTObject: Cannot find Trigger: " << arg << endl;
    gSystem->Exit(1);
  }

  int nobj = cms2.hlt_trigObjs_p4().at(trigIndx).size();
  if (nobj == 0 ) {
    cout << "ERROR: nobj == 0" << endl;
    gSystem->Exit(1);
  }

  if (objNumber > (nobj-1)) {
    cout << "ERROR: requested object number " << objNumber << " but we only have " << nobj <<endl;
    gSystem->Exit(1);
  }

  return cms2.hlt_trigObjs_p4().at(trigIndx).at(objNumber);

}

// trigger id
int idHLTObject(const char* arg, int objNumber){

  TString HLTTrigger( arg );
  int trigIndx = -1;
  vector<TString>::const_iterator begin_it = cms2.hlt_trigNames().begin();
  vector<TString>::const_iterator end_it = cms2.hlt_trigNames().end();
  vector<TString>::const_iterator found_it = find(begin_it, end_it, HLTTrigger );
  if( (found_it != end_it) ){
    trigIndx = found_it - begin_it;
    //cout << "idHLTObject: Found Trigger: " << arg << endl;
  }
  else {
    cout << "idHLTObject: Cannot find Trigger: " << arg << endl;
    gSystem->Exit(1);
  }

  int nobj = cms2.hlt_trigObjs_id().at(trigIndx).size();
  if (nobj == 0 ) {
    cout << "ERROR: nobj == 0" << endl;
    gSystem->Exit(1);
  }

  if (objNumber > (nobj-1)) {
    cout << "ERROR: requested object number " << objNumber << " but we only have " << nobj <<endl;
    gSystem->Exit(1);
  }

  return cms2.hlt_trigObjs_id().at(trigIndx).at(objNumber);

}



//--------------------------------------------------------
// Returns the number of objects passing a given trigger
// Returns zero if the trigger failed
// Returns -1 if the trigger passed but no onjects were found
//--------------------------------------------------------
int nHLTObjects(const char* arg ){

  // put the trigger name into a string
  TString HLTTrigger( arg );

  // Did the trigger pass?
  if ( !(cms2.passHLTTrigger(HLTTrigger)) ) return 0;

  // The trigger passed, see how many associated objects there are
  int trigIndx = -1;
  vector<TString>::const_iterator begin_it = cms2.hlt_trigNames().begin();
  vector<TString>::const_iterator end_it = cms2.hlt_trigNames().end();
  vector<TString>::const_iterator found_it = find(begin_it, end_it, HLTTrigger );
  if( (found_it != end_it) ){
    trigIndx = found_it - begin_it;
    //cout << "nHLTObjects: Found Trigger: " << arg << endl;
  }
  else {
    cout << "nHLTObjects: Cannot find Trigger " << arg << endl;
    return 0;
  }

  int nobj = cms2.hlt_trigObjs_p4().at(trigIndx).size();
  /*
  for( unsigned int i=0; i < nobj; i++ ){
    cout << "\t" << i << ", (pt, eta, phi): " << cms2.hlt_trigObjs_p4().at(trigIndx).at(i).pt() << " "
                  << cms2.hlt_trigObjs_p4().at(trigIndx).at(i).eta() << " " << cms2.hlt_trigObjs_p4().at(trigIndx).at(i).phi() << endl;
  }
  */

  // cout << " Number of jets = " << njets << endl;

  if (nobj == 0) return -1;
  return nobj;
}
//---------------------------------------------
// Utility to print trigger names in the file
//---------------------------------------------
void PrintTriggers(){
  for( unsigned int i = 0; i < cms2.hlt_trigNames().size(); i++ ){
    cout << cms2.passHLTTrigger(cms2.hlt_trigNames().at(i).Data()) << "\t"
         << cms2.hlt_prescales().at(i) << "\t" 
         << cms2.hlt_trigNames().at(i).Data() << endl;

  } 
  cout << endl;
}

//---------------------------------------------
// Check if trigger is unprescaled and passes
//---------------------------------------------
bool passUnprescaledHLTTrigger(const char* arg){

  // put the trigger name into a string
  TString HLTTrigger( arg );

  // Did the trigger pass?
  if ( !(cms2.passHLTTrigger(HLTTrigger)) ) return false;

  // The trigger passed, check the pre-scale
  int trigIndx = -1;
  vector<TString>::const_iterator begin_it = cms2.hlt_trigNames().begin();
  vector<TString>::const_iterator end_it   = cms2.hlt_trigNames().end();
  vector<TString>::const_iterator found_it = find(begin_it, end_it, HLTTrigger );
  if( (found_it != end_it) ){
    trigIndx = found_it - begin_it;
  }
  else {
    //this should not happen
    cout << "passUnprescaledTrigger: Cannot find Trigger " << arg << endl;
    return false;
  }

  //sanity check (this should not happen)
  if( strcmp( arg , cms2.hlt_trigNames().at(trigIndx) ) != 0 ){
    cout << "Error! trig names don't match" << endl;
    cout << "Found trig name " << cms2.hlt_trigNames().at(trigIndx) << endl;
    cout << "Prescale        " << cms2.hlt_prescales().at(trigIndx) << endl;
    exit(0);
  }

  //return true only if pre-scale = 1
  if( cms2.hlt_prescales().at(trigIndx) == 1 ) return true;

  return false;

}

//---------------------------------------------
// Check if trigger is unprescaled and passes
// for a specific object, specified by a p4
//---------------------------------------------
bool passUnprescaledHLTTrigger(const char* arg, const LorentzVector &obj){

  // put the trigger name into a string
  TString HLTTrigger( arg );

  // find the index of this trigger
  int trigIdx = -1;
  vector<TString>::const_iterator begin_it = cms2.hlt_trigNames().begin();
  vector<TString>::const_iterator end_it = cms2.hlt_trigNames().end();
  vector<TString>::const_iterator found_it = find(begin_it, end_it, HLTTrigger);
  if(found_it != end_it) trigIdx = found_it - begin_it;
  else return false; // trigger was not found

  // get the vector of p4 passing this trigger
  std::vector<LorentzVector> trigObjs = cms2.hlt_trigObjs_p4()[trigIdx];

  // if no trigger objects then fail
  if (trigObjs.size() == 0) return false; 

  // does the trigger match this lepton
  float drMin = 999.99;
  for (size_t i = 0; i < trigObjs.size(); ++i)
  {
    float dr = ROOT::Math::VectorUtil::DeltaR(trigObjs[i], obj);
    if (dr < drMin) drMin = dr;
  }

  // if the closest trigger object
  // is further than 0.1 then fail
  if (drMin > 0.1) return false;

  // if we got to here then
  // the trigger passed, check the pre-scale

  //sanity check (this should not happen)
  if( strcmp( arg , cms2.hlt_trigNames().at(trigIdx) ) != 0 ){
    cout << "Error! trig names don't match" << endl;
    cout << "Found trig name " << cms2.hlt_trigNames().at(trigIdx) << endl;
    cout << "Prescale        " << cms2.hlt_prescales().at(trigIdx) << endl;
    exit(0);
  }

  //return true only if pre-scale = 1
  if( cms2.hlt_prescales().at(trigIdx) == 1 ) return true;

  return false;

}


//this function returns the HLT pre-scale for a given trigger name

int HLT_prescale( const char* arg ){

 // put the trigger name into a string
  TString HLTTrigger( arg );

  // Did the trigger pass?
  if ( !(cms2.passHLTTrigger(HLTTrigger)) ) return -1;

  // The trigger passed, check the pre-scale
  int trigIndx = -1;
  vector<TString>::const_iterator begin_it = cms2.hlt_trigNames().begin();
  vector<TString>::const_iterator end_it   = cms2.hlt_trigNames().end();
  vector<TString>::const_iterator found_it = find(begin_it, end_it, HLTTrigger );
  if( (found_it != end_it) ){
    trigIndx = found_it - begin_it;
  }
  else {
    //this should not happen
    cout << "passUnprescaledTrigger: Cannot find Trigger " << arg << endl;
    return false;
  }

  //sanity check (this should not happen)
  if( strcmp( arg , cms2.hlt_trigNames().at(trigIndx) ) != 0 ){
    cout << "Error! trig names don't match" << endl;
    cout << "Found trig name " << cms2.hlt_trigNames().at(trigIndx) << endl;
    cout << "Prescale        " << cms2.hlt_prescales().at(trigIndx) << endl;
    exit(0);
  }

  //return prescale
  return cms2.hlt_prescales().at(trigIndx);


}


//this function returns the L1 prescale for a given trigger name

int L1_prescale( const char* arg ){

  // put the trigger name into a string
  TString trigger( arg );

  // Did the trigger pass?
  if ( !(cms2.passL1Trigger(trigger)) ) return -1;

  // The trigger passed, check the pre-scale
  int trigIndx = -1;
  vector<TString>::const_iterator begin_it = cms2.l1_trigNames().begin();
  vector<TString>::const_iterator end_it   = cms2.l1_trigNames().end();
  vector<TString>::const_iterator found_it = find(begin_it, end_it, trigger );
  if( (found_it != end_it) ){
    trigIndx = found_it - begin_it;
  }
  else {
    //this should not happen
    cout << "L1_prescale: Cannot find Trigger " << arg << endl;
    return false;
  }

  //sanity check (this should not happen)
  if( strcmp( arg , cms2.l1_trigNames().at(trigIndx) ) != 0 ){
    cout << "Error! trig names don't match" << endl;
    cout << "Found trig name " << cms2.l1_trigNames().at(trigIndx) << endl;
    cout << "Prescale        " << cms2.l1_prescales().at(trigIndx) << endl;
    exit(0);
  }

  return cms2.l1_prescales().at(trigIndx);

}
