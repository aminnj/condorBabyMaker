// Header
#include "utilities.h"

// C++ includes
#include <vector>

// ROOT includes
#include "Math/VectorUtil.h"
#include "CMS2.h"
using std::vector;

//return true if one of the leptons is the same in both hyps
bool hypsOverlap(int idxa, int idxb){
  int idlta   = cms2.hyp_lt_id()[idxa];
  int idlla   = cms2.hyp_ll_id()[idxa];
  int ilta    = cms2.hyp_lt_index()[idxa];
  int illa    = cms2.hyp_ll_index()[idxa];
  int idltb   = cms2.hyp_lt_id()[idxb];
  int idllb   = cms2.hyp_ll_id()[idxb];
  int iltb    = cms2.hyp_lt_index()[idxb];
  int illb    = cms2.hyp_ll_index()[idxb];
  int matches = (idlta == idltb && ilta == iltb) + (idlla == idllb && illa == illb) + (idlta == idllb && ilta == illb) + (idlla == idltb && illa == iltb);
  return matches>0;
}

int match4vector(const LorentzVector &lvec, const vector<LorentzVector> &vec, double cut=10.0 ){

  if( vec.size() == 0 ) return -1;
  //cout << "size of vec = " << vec.size() << endl;
  double dR = cut; 
  double x;
  int iret = -1;
  for ( unsigned int i=0; i < vec.size();++i) {
    x = ROOT::Math::VectorUtil::DeltaR( lvec, vec[i] );
    if (x < dR ) {dR = x; iret = i;}
  }
  return iret;
}

std::vector<LorentzVector> p4sInCone(const LorentzVector &refvec, const std::vector<LorentzVector> &invec, double coneSize=0.5 ) {
  vector<LorentzVector> result;
  if ( invec.size() == 0 ) return result;
  double dR = coneSize; 
  double x = 0.0;
  for ( unsigned int i=0; i < invec.size();++i) {
    x = ROOT::Math::VectorUtil::DeltaR( refvec, invec[i] );
    if (x < dR ) {result.push_back(invec[i]);}
  }
  return result;
}

std::vector<unsigned int> idxInCone(const LorentzVector &refvec, const std::vector<LorentzVector> &invec, double coneSize=0.5 ) {
  vector<unsigned int > result;
  if ( invec.size() == 0 ) return result;
  double dR = coneSize; 
  for ( unsigned int i=0; i < invec.size();++i) {
    if ( ROOT::Math::VectorUtil::DeltaR( refvec, invec[i] ) < dR ) {result.push_back(i);}
  }
  return result;
}

/*
  
- Depricated

double trkIsolation(int trk_index) {
  //
  // calculate track isolation following electron isolation definition in ElectronMaker.cc
  //
  // sum up all track.pt around track if track fulfills:
  // dR < 0.3
  // dR > 0.01
  // d0 < 0.1
  // dZ < 0.5
  // pT >= 1.0
  // nHits > 7

  float dRConeMin = 0.01;
  float dRConeMax = 0.3;
  float vtxDiffZMax = 0.5;
  float tkVtxDMax = 0.1;
  float ptMin = 1.0;
  int nHits = 7;

  double isoResult = -10.;
  if( cms2.trks_trk_p4().size() == 0 ) {
    std::cout << "Configuration Error: track collection is not set!" <<std::endl;
    return isoResult;
  }

  double sumPt = 0;

  for ( unsigned int trk = 0;
	trk < cms2.trks_trk_p4().size();
	++trk ) {
    if ( cms2.trks_z0()[trk] > tkVtxDMax ) continue;
    double pt = cms2.trks_trk_p4()[trk].pt();
    if (  pt < ptMin ) continue;
    if ( cms2.trks_validHits()[trk] <= nHits ) continue;
    double dR = ROOT::Math::VectorUtil::DeltaR(cms2.trks_trk_p4()[trk_index], cms2.trks_trk_p4()[trk]);
    if (dR < dRConeMin) continue;
    if ( dR > dRConeMax ) continue;
    double dZ = TMath::Abs(cms2.trks_z0()[trk_index] - cms2.trks_z0()[trk]);
    if ( dZ >= vtxDiffZMax) continue;
    sumPt += pt;
  }

  return sumPt;
}
*/

/*

- This is not the right place for this ( utilities is intended for possible inclusion by other selections, this belongs in Tools, or maybe eventSelections )
- Classes need to be defined in the header file

//this function reads in a list of events from a text file
//and checks if the current event is present in the list
//USAGE:
//EventInListIdentified eid;
//bool found = eid.event_in_list( EventInList() );

class EventInList{
public:
  EventInList(){
    run_      = cms2.evt_run();
    lumi_     = cms2.evt_lumiBlock();
    event_    = cms2.evt_event();
    dataset_  = cms2.evt_dataset();
  }
  ~EventInList() {}

  int run          () const { return run_;      }
  int lumi         () const { return lumi_;     }
  int event        () const { return event_;    }
  TString dataset  () const { return dataset_;  }
  
private:
  int run_, lumi_, event_;
  TString dataset_;
  
};

class EventInListIdentifier{
public:
  EventInListIdentifier(char* inputFile = "eventList.txt"){
    runs.clear();
    lumis.clear();
    events.clear();

    ifstream ifile( inputFile );
    if( !ifile.is_open() ){
      cout << "ERROR CANNOT FIND EVENT LIST " << inputFile << endl;
      exit(0);
    }
    cout << "Opened event list " << inputFile << endl;
    
    while( ifile >> temprun >> templumi >> tempevent ){
      runs.push_back( temprun );
      lumis.push_back( templumi );
      events.push_back( tempevent );
      cout << "Added: " << temprun << " " << templumi << " " << tempevent << endl;
    }
  }
  ~EventInListIdentifier() {}
  
  bool event_in_list(const EventInList &id){
    
    for( unsigned int i = 0 ; i < runs.size() ; ++i ){
    
      if( id.run() == runs.at(i) && id.lumi() == lumis.at(i) && id.event() == events.at(i) ){
        cout << "Found event in list! " << id.dataset() << " " << id.run() << " " << id.lumi() << " " << id.event() << endl;
        return true;
      }
    }
    
    return false;
  }
  
private:
  vector<int> runs;
  vector<int> lumis;
  vector<int> events;
  int temprun;
  int templumi;
  int tempevent;
};

*/
