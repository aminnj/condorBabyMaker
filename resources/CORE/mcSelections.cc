#include <math.h>
#include <stdlib.h>
#include <set>
#include "TDatabasePDG.h"
#include "Math/VectorUtil.h"
#include "CMS2.h"
#include "mcSelections.h"
#include "Math/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

//-------------------------------------------------
// Auxiliary function to scan the doc line and 
// identify DY-> ee vs mm vs tt
//-------------------------------------------------
int getDrellYanType() {
  bool foundEP = false;
  bool foundEM = false;
  bool foundMP = false;
  bool foundMM = false;
  bool foundTP = false;
  bool foundTM = false;
  for (unsigned int i = 0; i < cms2.genps_id().size(); ++i) {
    if ( cms2.genps_id_mother().at(i) == 23 ){
      switch ( TMath::Abs(cms2.genps_id().at(i)) ){
      case 11:
	return 0;
	break;
      case 13:
	return 1;
	break;
      case 15:
	return 2;
	break;
      default:
	break;
      }
    }
    switch ( cms2.genps_id().at(i) ){
    case 11:
      foundEM = true;
      break;
    case -11:
      foundEP = true;
      break;
    case 13:
      foundMM = true;
      break;
    case -13:
      foundMP = true;
      break;
    case 15:
      foundTM = true;
      break;
    case -15:
      foundTP = true;
      break;
    default:
      break;
    }
  }
  
  if ( foundEP && foundEM ) return 0;  //DY->ee
  if ( foundMP && foundMM ) return 1;  //DY->mm
  if ( foundTP && foundTM ) return 2;  //DY->tautau
  std::cout << "Does not look like a DY event" << std::endl;
  return 999;
}

int getVVType() {
  // types:
  //   0 - WW
  //   1 - WZ
  //   2 - ZZ
  unsigned int nZ(0);
  unsigned int nW(0);
  std::vector<std::vector<int> > leptons;
  std::vector<int> mothers;

  bool verbose = false;

  for (unsigned int i = 0; i < cms2.genps_id().size(); ++i) {
    int pid = cms2.genps_id().at(i);
    int mid = cms2.genps_id_mother().at(i);
    if ( verbose ) std::cout << "Gen particle id: " << pid << ",\t mother id: " << mid <<std::endl;
    if ( abs(pid)<11 || abs(pid)>16 ) continue;
    if ( mid == 23 ) ++nZ;
    if ( abs(mid) == 24 ) ++nW;
    // now we need to really understand the pattern.
    unsigned int mIndex = 0;
    while ( mIndex < mothers.size() && mid != mothers[mIndex] ) ++mIndex;
    if ( mIndex == mothers.size() ) {
      mothers.push_back(mid);
      leptons.push_back(std::vector<int>());
    }
    leptons[mIndex].push_back(pid);
    if (mothers.size()>3){
      if (verbose) std::cout << "WARNING: failed to identify event (too many mothers)" << std::endl;
      return 999;
    }
  }

  if ( nZ == 4 ) {
    if ( verbose ) std::cout << "Event type ZZ" << std::endl;
    return 2;
  }
  if ( nW == 4 ) {
    if ( verbose ) std::cout << "Event type WW" << std::endl;
    return 0;
  }
  if ( nW == 2 && nZ == 2 ) {
    if ( verbose ) std::cout << "Event type WZ" << std::endl;
    return 1;
  }
  unsigned int nNus(0);
  for ( unsigned int i=0; i<mothers.size(); ++i ){
      nNus += leptons[i].size();
  }
  if ( mothers.size() < 3 && nNus == 4){
    for ( unsigned int i=0; i<mothers.size(); ++i ){
      if ( mothers[i] != 23 && abs(mothers[i]) != 24 ){
	if( leptons[i].size() != 2 && leptons[i].size() != 4){
	  if (verbose) std::cout << "WARNING: failed to identify event (unexpected number of daughters)" << std::endl;
	  if (verbose) std::cout << "\tnumber of daughters for first mother: " <<  leptons[0].size() << std::endl;
	  if (verbose) std::cout << "\tnumber of daughters for second mother: " <<  leptons[1].size() << std::endl;
	  return 999;
	}
	if ( abs(leptons[i][0]) == abs(leptons[i][1]) )
	  nZ += 2;
	else
	  nW += 2;
	if ( leptons[i].size()==4 ){
	  // now it's a wild guess, it's fraction should be small
	  if ( abs(leptons[i][2]) == abs(leptons[i][3]) )
	    nZ += 2;
	  else
	    nW += 2;
	}
      }
    }
  } else {
    // here be dragons
    
    // if we have 2 leptons and 3 neutrinos and they all of the same
    // generation, we assume it's ZZ (can be WZ also), but if
    // neutrinos are from different generations, than we conclude it's
    // WZ. 
    
    std::set<int> nus;
    for ( unsigned int i=0; i<mothers.size(); ++i )
      for ( unsigned int j=0; j<leptons[i].size(); ++j ) 
	if ( abs(leptons[i][j]) == 12 ||
	     abs(leptons[i][j]) == 14 ||
	     abs(leptons[i][j]) == 16 )
	  nus.insert(abs(leptons[i][j]));
    
    if ( nNus == 5 ){
      if ( nus.size() == 1 ) return 2;
      if ( nus.size() == 2 ) return 1;
    }
    
    if ( verbose ) std::cout << "WARNING: failed to identify event" << std::endl;
    return 999;
  }

  if ( nZ+nW != 4 ){
    if (verbose) std::cout << "WARNING: failed to identify event (wrong number of bosons)" << std::endl;
    if (verbose) std::cout << "\tfirst mother id: " << mothers[0] << std::endl;
    if (verbose) std::cout << "\tsecond mother id: " << mothers[1] << std::endl;
    if (verbose) std::cout << "\tnumber of daughters for first mother: " << leptons[0].size() << std::endl;
    if (verbose) std::cout << "\tnumber of daughters for second mother: " << leptons[1].size() << std::endl;
    if (verbose) std::cout << "\tnumber of Zs: " << nZ << std::endl;
    if (verbose) std::cout << "\tnumber of Ws: " << nW << std::endl;
    return 999;
  }

  if ( nZ == 4 ) {
    if ( verbose ) std::cout << "Event type ZZ" << std::endl;
    return 2;
  }
  if ( nW == 4 ) {
    if ( verbose ) std::cout << "Event type WW" << std::endl;
    return 0;
  }
  // this covers screws in logic, i.e. most hard to identify events end up being WZ
  if ( verbose ) std::cout << "Event type WZ (can be wrong)" << std::endl;
  return 1;
}

//-------------------------------------------------
// Auxiliary function to scan the doc line and 
// identify DY-> ee vs mm vs tt
//-------------------------------------------------
int getWType() 
{
     bool foundE = false;
     bool foundNuE = false;
     bool foundM = false;
     bool foundNuM = false;
     bool foundT = false;
     bool foundNuT = false;
     for (unsigned int i = 0; i < cms2.genps_id().size(); ++i) {
	  if ( abs(cms2.genps_id_mother().at(i)) == 24 ){
	       switch ( TMath::Abs(cms2.genps_id().at(i)) ){
	       case 11:
		    return 0;
		    break;
	       case 13:
		    return 1;
		    break;
	       case 15:
		    return 2;
		    break;
	       default:
		    break;
	       }
	  }
	  switch ( abs(cms2.genps_id().at(i)) ){
	  case 11:
	       foundE = true;
	       break;
	  case 12:
	       foundNuE = true;
	       break;
	  case 13:
	       foundM = true;
	       break;
	  case 14:
	       foundNuM = true;
	       break;
	  case 15:
	       foundT = true;
	       break;
	  case 16:
	       foundNuT = true;
	       break;
	  default:
	       break;
	  }
     }
     
     if ( foundE && foundNuE ) return 0;  //W->e
     if ( foundM && foundNuM ) return 1;  //W->m
     if ( foundT && foundNuT ) return 2;  //W->t
     std::cout << "Does not look like a W event" << std::endl;
     return 999;
}

//--------------------------------------------
// Booleans for DY
//------------------------------------------
bool isDYee() {
  if (getDrellYanType() == 0) return true;
  return false;
}
bool isDYmm() {
  if (getDrellYanType() == 1) return true;
  return false;
}
bool isDYtt() {
  if (getDrellYanType() == 2) return true;
  return false;
}

//--------------------------------------------
// Booleans for Wjets
//------------------------------------------
bool isWe() 
{
     if (getWType() == 0) return true;
     return false;
}
bool isWm() 
{
     if (getWType() == 1) return true;
     return false;
}
bool isWt() 
{
     if (getWType() == 2) return true;
     return false;
}

//--------------------------------------------
// Booleans for VVjets
//------------------------------------------

bool isWW() {
  if (getVVType() == 0) return true;
  return false;
}

bool isWZ() {
  if (getVVType() == 1) return true;
  return false;
}

bool isZZ() {
  if (getVVType() == 2) return true;
  return false;
}

//--------------------------------------------
// ZZ type:
// 0 for Z1 --> ee, mm; Z2 --> ee, mm
// 1 for Z1 --> ee, mm; Z2 --> tt (and v.v.)
// 2 for Z1 --> tt; Z2 --> tt
// 995 to 999 otherwise
//------------------------------------------
int getZZType() 
{
     int foundEP = 0;
     int foundEM = 0;
     int foundMP = 0;
     int foundMM = 0;
     int foundTP = 0;
     int foundTM = 0;
     for (unsigned int i = 0; i < cms2.genps_id().size(); ++i) {
	  switch ( cms2.genps_id().at(i) ){
	  case 11:
	       foundEM++;
	       break;
	  case -11:
	       foundEP++;
	       break;
	  case 13:
	       foundMM++;
	       break;
	  case -13:
	       foundMP++;
	       break;
	  case 15:
	       foundTM++;
	       break;
	  case -15:
	       foundTP++;
	       break;
	  default:
	       break;
	  }
     }
  
     if (foundEM == foundEP && foundMM == foundMP && (foundEM != 0 || foundMM != 0)) {
	  // both Zs decay to e or mu
	  if (foundEM + foundMM == 2)
	       return 0;
	  // one Z decays to e or mu
	  else if (foundEM + foundMM == 1) 
	       // other Z decays to tau
	       if (foundTP == 1 && foundTM == 1)
		    return 1;
	       else return 995;
	  else return 996;
     } else if (foundEM == 0 && foundEP == 0 && foundMM == 0 && foundMP == 0) {
	  // both Zs decay to tau
	  if (foundTP == 2 && foundTM == 2)
	       return 2;
	  else return 997;
     } else return 998;
     return 999;
}

//------------------------------------------------------------
// Not a selection function per se, but useful nonetheless:
// dumps the documentation lines for this event
//------------------------------------------------------------
int dumpDocLines() {

  //////////////////////////////////
  // Initialize particle database //
  //////////////////////////////////
  TDatabasePDG *pdg = new TDatabasePDG();

  //////////////////
  // Print Header //
  ////////////////// 
  cout << "                " << "   pt    " << "  phi  " << "      eta   " << "    mass  " << "status " << "Mother  " << endl;     
  std::cout << "---------------------------------------------------------------------" << std::endl;

  /////////////////////////////
  // Loop over gen particles //
  /////////////////////////////
  int size = cms2.genps_id().size();
  for (int j=0; j<size; j++) {

    // mass
    float m2 = cms2.genps_p4().at(j).M2();
    float m  = m2 >= 0 ? sqrt(m2) : 0.0;

    //////////////////////////////////////////////////
    // Print information about the jth gen particle //
    //////////////////////////////////////////////////
    cout << setw(4)  << left  <<                    j << " "
         << setw(10) << left  <<                    pdg->GetParticle(cms2.genps_id().at(j))->GetName()        << " "
	       << setw(7)  << right << setprecision(4) << cms2.genps_p4().at(j).pt()                                << "  "
	       << setw(7)  << right << setprecision(4) << cms2.genps_p4().at(j).phi()                               << "  "
	       << setw(10) << right << setprecision(4) << cms2.genps_p4().at(j).eta()                               << "  "
	       << setw(7)  << right << setprecision(4) << m                                                         << "  "
         << setw(4)  << right <<                    cms2.genps_status().at(j)                                 << " "
         << setw(10) << left  <<                    pdg->GetParticle(cms2.genps_id_mother().at(j))->GetName() << " " 
         << endl;

    //////////////////////////////////////////////
    // Lepton daughters of the jth gen particle //
    //////////////////////////////////////////////
    if(cms2.genps_lepdaughter_id()[j].size() > 0) {
      cout << endl;
      cout << "  Daughters:" << endl;
      for(unsigned int i = 0; i < cms2.genps_lepdaughter_id()[j].size(); i++) {

        // mass
        float m2_daught = cms2.genps_lepdaughter_p4().at(j).at(i).M2();
	      float m_daught  = m2_daught >= 0 ? sqrt(m2_daught) : 0.0;

        ///////////////////////////////////////////
        // Print information about the daughters //
        ///////////////////////////////////////////
	      cout << setw(2)  << left  <<                    "    " << i << " "
	           << setw(10) << left  <<                    pdg->GetParticle(cms2.genps_lepdaughter_id().at(j).at(i))->GetName() << " "
	           << setw(7)  << right << setprecision(4) << cms2.genps_lepdaughter_p4().at(j).at(i).pt()                         << "  "
	           << setw(7)  << right << setprecision(4) << cms2.genps_lepdaughter_p4().at(j).at(i).phi()                        << "  "
	           << setw(10) << right << setprecision(4) << cms2.genps_lepdaughter_p4().at(j).at(i).eta()                        << "  "
	           << setw(7)  << right << setprecision(4) << m_daught                                                             << "  " 
             << endl;
      
      }
      cout << endl;
    }
  
  }
  delete pdg;
  return 0;
}

int ttbarconstituents(int i_hyp){
  // Categories:
  //WW = both leptons from W = 1
  //WO = one of the two leptons from W and the other not = 2
  //OO = neither of the two leptons is from a W = 3

  int lttype = leptonIsFromW(cms2.hyp_lt_index()[i_hyp],cms2.hyp_lt_id()[i_hyp]);
  int lltype = leptonIsFromW(cms2.hyp_ll_index()[i_hyp],cms2.hyp_ll_id()[i_hyp]);
  if (lltype > 0 && lttype > 0) return 1;
  else if( (lltype >0 && lttype <= 0) || (lttype >0 && lltype <=0) ) return 2;
  else if( (lltype <=0 && lttype <=0) )return 3;
  else { cout << "bug in ttbarconstituents"; return -999;}
}

//--------------------------------------------------------
// Determines if the lepton in question is from W/Z
// and if its charge is correct
//
// Note that if we have
//     W->lepton->lepton gamma
// where the gamma is at large angles and it is the
// gamma that gives the lepton signature in the detector,
// then this returns "not from W/Z".  This is by design.
//
// Note W->tau->lepton is tagged as "from W"
//
// Inputs:  idx   = index in the els or mus block
//          id    = lepton ID (11 or 13 or -11 or -13)
//
// Output:  2 = from W/Z incorrect charge
//          1 = from W/Z   correct charge
//          0 = not matched to a lepton (= fake)
//         -1 = lepton from b decay
//         -2 = lepton from c decay
//         -3 = lepton from some other source
//        
// Authors: Claudio in consultation with fkw 22-july-09    
//
// 22-nov-2010...by "accident" this code returns "fromW" (1 or 2)
// in many of the cases where the lepton is from a susy particle,
// because it looks at whether or not it is matched to a status==3
// lepton.  By this procedure is not 100% efficient.  We have now
// added a flag to try to do this more systematically.
//    Claudio & Derek
//---------------------------------------------------------
int leptonIsFromW(int idx, int id, bool alsoSusy) {

  // get the matches to status=1 and status=3
  int st1_id = 0;
  int st3_id = 0;
  int st1_motherid = 0;
  if (abs(id) == 11) {
    st1_id = cms2.els_mc_id()[idx];
    st3_id = cms2.els_mc3_id()[idx];
    st1_motherid = cms2.els_mc_motherid()[idx];
    //a true lepton from a W will have it's motherid==24
    //if the lepton comes from a tau decay that comes from a W, 
    //we have to do some work to trace the parentage
    //to do this, we have to go to the status==3 block because 
    //the daughter info is not in status==1
    if(abs(st1_motherid)==15) {
      bool foundelectronneutrino = false; //ensures that the matched electron from a tau really came from a W
      for(unsigned int i = 0; i < cms2.genps_id().size(); i++) {//status 3 loop
	if(abs(cms2.genps_id()[i]) == 15 ) { //make sure we get the right tau!
	  cms2.genps_lepdaughter_id()[i].size(); 
	  for(unsigned int j = 0; j < cms2.genps_lepdaughter_id()[i].size(); j++) { //loop over the tau's status1 daughters
	    if(abs(cms2.genps_lepdaughter_id()[i][j]) == 12)
	      foundelectronneutrino = true;
	    float dr = ROOT::Math::VectorUtil::DeltaR(cms2.els_mc_p4()[idx], cms2.genps_lepdaughter_p4()[i][j]);
	    if (dr < 0.0001) { //should be the same exact status==1 gen particle!
	      st1_motherid = cms2.genps_id_mother()[i];
	      continue;
	    }//if (dr < 0.0001)
	  }//loop over the tau's daughters
	  if(!foundelectronneutrino)
	    st1_motherid = -9999;
	}//tau
      }//status 3 loop
    }//if(abs(st1_motherid)==15) {
  } else if (abs(id) == 13) {
    st1_id = cms2.mus_mc_id()[idx];
    st3_id = cms2.mus_mc3_id()[idx];
    st1_motherid = cms2.mus_mc_motherid()[idx];
    //a true lepton from a W will have it's motherid==24
    //if the lepton comes from a tau decay that comes from a W, 
    //we have to do some work to trace the parentage
    //to do this, we have to go to the status==3 block because 
    //the daughter info is not in status==1
    if(abs(st1_motherid)==15) {
      bool foundmuonneutrino = false;
      for(unsigned int i = 0; i < cms2.genps_id().size(); i++) {//status 3 loop
	if(abs(cms2.genps_id()[i]) == 15 ) { //make sure we get the right tau!
	  cms2.genps_lepdaughter_id()[i].size();
	  for(unsigned int j = 0; j < cms2.genps_lepdaughter_id()[i].size(); j++) {//loop over the tau's status1 daughters
	    if(abs(cms2.genps_lepdaughter_id()[i][j]) == 14)
	      foundmuonneutrino = true;
	    float dr = ROOT::Math::VectorUtil::DeltaR(cms2.mus_mc_p4()[idx], cms2.genps_lepdaughter_p4()[i][j]);
	    if (dr < 0.0001) { //should be the same exact status==1 gen particle!
 	      st1_motherid = cms2.genps_id_mother()[i];
	      continue;
	    }//if (dr < 0.0001)
	  }//loop over the tau's daughters
	  if(!foundmuonneutrino)
	    st1_motherid = -9999;
	}//tau
      }//status 3 loop
    }//if(abs(st1_motherid)==15) {
  } else {
    std::cout << "You fool.  Give me +/- 11 or +/- 13 please" << std::endl;
    return false;
  }


  // debug
  // std::cout << "id=" << id << " st1_id=" << st1_id;
  // std::cout << " st3_id=" << st3_id;
  // std::cout << " st1_motherid=" << st1_motherid << std::endl;

  // Step 1
  // Look at status 1 match, it should be either a lepton or
  // a photon if it comes from W/Z.
  // The photon case takes care of collinear FSR
  if ( !(abs(st1_id) == abs(id) || st1_id == 22)) return 0;

  // Step 2
  // If the status 1 match is a photon, its mother must be
  // a lepton.  Otherwise it is not FSR
  if (st1_id == 22) {
    if (abs(st1_motherid) != abs(id)) return 0;
  }

  // At this point we are matched (perhaps via FSR) to
  // a status 1 lepton.  This means that we are left with
  // leptons from W, taus, bottom, charm, as well as dalitz decays

  // Step 3
  // A no-brainer: pick up vanilla W->lepton decays
  // (should probably add Higgs, SUSY, W' etc...not for now)
  if (st1_id ==  id && abs(st1_motherid) == 24) return 1; // W
  if (st1_id == -id && abs(st1_motherid) == 24) return 2; // W
  if (st1_id ==  id &&   st1_motherid    == 23) return 1; // Z
  if (st1_id == -id &&   st1_motherid    == 23) return 2; // Z
  if ( alsoSusy ) {
    if (st1_id ==  id && abs(st1_motherid) > 1e6) return 1; // exotica
    if (st1_id == -id && abs(st1_motherid) > 1e6) return 2; // exotica
  }

  // Step 4
  // Another no-brainer: pick up leptons matched to status=3
  // leptons.  This should take care of collinear FSR
  // This also picks up a bunch of SUSY decays
  if (st3_id ==  id) return 1;
  if (st3_id == -id) return 2;
 
  // Step 5
  // Now we need to go after the W->tau->lepton.  
  // We exploit the fact that in t->W->tau the tau shows up
  // at status=3.  We also use the fact that the tau decay products
  // are highly boosted, so the direction of the status=3 tau and
  // the lepton from tau decays are about the same
  //
  // We do not use the status=1 links because there is not
  // enough information to distinguish
  // W->tau->lepton  or W->tau->lepton gamma
  //  from
  // B->tau->lepton or B->tau->lepton gamma
  //if (abs(st3_id) == 15 && id*st3_id > 0) return 1;
  //if (abs(st3_id) == 15 && id*st3_id < 0) return 2;
  if(abs(st3_id) == 15) {
    //have to find the index of the status3 particle by dR
    //because the indices are buggy
    unsigned int mc3idx = 999999;
    LorentzVector lepp4 =  abs(id)==11 ? cms2.els_p4()[idx] : cms2.mus_p4()[idx];
    double mindR = 9999;
    for(unsigned int i = 0; i < cms2.genps_id().size(); i++) {
      float dr = ROOT::Math::VectorUtil::DeltaR(lepp4, cms2.genps_p4()[i]);
      if(dr < mindR) {
	mindR = dr;
	mc3idx = i;
      }
    }
    bool foundElOrMuNu = false;    
    for(unsigned int i = 0; i < cms2.genps_lepdaughter_p4()[mc3idx].size(); i++) {
      if(abs(cms2.genps_lepdaughter_id()[mc3idx][i]) == 12 || abs(cms2.genps_lepdaughter_id()[mc3idx][i]) == 14)
	foundElOrMuNu = true;
    }
    if(!foundElOrMuNu) //comes from a hadronic decay of the tau
      return -3;
    if(id*st3_id > 0) 
      return 1;
    else return 2;
  }
  
  
  // Step 6
  // If we get here, we have a non-W lepton
  // Now we figure out if it is from b, c, or "other"
  // There are a couple of caveats
  // (a) b/c --> lepton --> lepton gamma (ie FSR) is labelled as "other"
  // (b) b   --> tau --> lepton is labelled as "other"
  if ( abs(st1_id) == abs(id) && idIsBeauty(st1_motherid)) return -1;
  if ( abs(st1_id) == abs(id) && idIsCharm(st1_motherid))  return -2;
  return -3;
}
//---------------------------------------------------------

bool trueGammaFromMuon(int electron) {
  // true gamma reconstructed as electron 
  // gamma coming from muon
  if(TMath::Abs(cms2.els_mc_id()[electron]) == 22 && TMath::Abs(cms2.els_mc_motherid()[electron]) == 13) { // ok, abs of photon makes no sense ;)
    //    std::cout<<"Gamma from muon event - r: " << cms2.evt_run() << " e: " << cms2.evt_event() << " l: " << cms2.evt_lumiBlock() << std::endl;
    return true;
  }
  return false;
}

// count genp leptons
//-------------------------------------------------- 
// Returns the number of e,mu, and tau in the doc lines 
//----------------------------------------------------- 
int leptonGenpCount(int& nele, int& nmuon, int& ntau) { 
  nele=0; 
  nmuon=0; 
  ntau=0; 
  int size = cms2.genps_id().size(); 
  for (int jj=0; jj<size; jj++) { 
    if (abs(cms2.genps_id().at(jj)) == 11) nele++; 
    if (abs(cms2.genps_id().at(jj)) == 13) nmuon++; 
    if (abs(cms2.genps_id().at(jj)) == 15) ntau++; 
  }
  
  return nele + nmuon + ntau;
 
} 

int leptonGenpCount_lepTauDecays(int& nele, int& nmuon, int& ntau) { 
  nele=0; 
  nmuon=0; 
  ntau=0; 
  int size = cms2.genps_id().size(); 
  for (int jj=0; jj<size; jj++) { 
    if (abs(cms2.genps_id().at(jj)) == 11) nele++; 
    if (abs(cms2.genps_id().at(jj)) == 13) nmuon++; 
    if (abs(cms2.genps_id().at(jj)) == 15) {
      for(unsigned int kk = 0; kk < cms2.genps_lepdaughter_id()[jj].size(); kk++) {
	int daughter = abs(cms2.genps_lepdaughter_id()[jj][kk]);
	//we count neutrino's because that guarantees that 
	//there is a corresponding lepton and that it comes from
	// a leptonic tau decay. You can get electrons from converted photons
	//which are radiated by charged pions from the tau decay but thats
	//hadronic and we don't care for those 
	if( daughter == 12 || daughter == 14)
	  ntau++; 
      }//daughter loop
    }//if tau
  }//genps loop
  
  return nele + nmuon + ntau;
}


//---------------------------------------------------------
int genpDileptonType(){
  //0 mumu; 1 emu; 2 ee
  
  unsigned int nmus = 0;
  unsigned int nels = 0;
  int size = cms2.genps_id().size();
  for (int jj=0; jj<size; jj++) {
    if (abs(cms2.genps_id().at(jj)) == 11) nels++;
    if (abs(cms2.genps_id().at(jj)) == 13) nmus++;
  }

  if ((nels + nmus) != 2){
    return -1;
  }

  int dilType = -1;
  if (nmus == 2) dilType = 0;
  if (nels == 2) dilType = 2;
  if (nels == 1 && nmus == 1) dilType = 1;
  return dilType;
}

// -----------------------------------------
// MC helper functions for fakerate tests:
// -----------------------------------------
int elFakeMCCategory(int i_el) {
  int category = -1;
  if(
     (abs(cms2.els_mc_id()[i_el])       == 11  && 
      abs(cms2.els_mc_motherid()[i_el]) == 22) ||
     (abs(cms2.els_mc_id()[i_el])       == 22) ||
     (abs(cms2.els_mc_id()[i_el])        > 100 && 
      abs(cms2.els_mc_id()[i_el])        < 200)||
     (abs(cms2.els_mc_id()[i_el])       == 11  && 
      abs(cms2.els_mc_motherid()[i_el]) == 111)
     ) {
    // electrons from gamma (conversion)
    category = 1;
  }
  else if(
          (abs(cms2.els_mc_id()[i_el]) > 200     && 
           abs(cms2.els_mc_id()[i_el]) < 400  )  ||
          (abs(cms2.els_mc_id()[i_el]) > 2000    && 
           abs(cms2.els_mc_id()[i_el]) < 4000 )  ||
          (abs(cms2.els_mc_id()[i_el]) == 11 && 
           abs(cms2.els_mc_motherid()[i_el]) > 200     && 
           abs(cms2.els_mc_motherid()[i_el]) < 400  )  || 
          (abs(cms2.els_mc_id()[i_el]) == 11 && 
           abs(cms2.els_mc_motherid()[i_el]) > 2000    && 
           abs(cms2.els_mc_motherid()[i_el]) < 4000 )  
          ) {
    // electron candidate or its mother is a light hadron
    category = 2;
  }
  else if( ( abs(cms2.els_mc_id()[i_el]) == 11 
             && abs(cms2.els_mc_motherid()[i_el]) >=400
             && abs(cms2.els_mc_motherid()[i_el]) <=600 )  || 
           ( abs(cms2.els_mc_id()[i_el]) == 11 
             && abs(cms2.els_mc_motherid()[i_el]) >=4000
             && abs(cms2.els_mc_motherid()[i_el]) <=6000 )
           ) {
    // heavy hadrons
    category = 3;
  }
  else {
    // the rest
    category = 4;
  }
  return category;
}
int muFakeMCCategory(int i_mu) {
  int category = -1;
  
  if( // punchthrough / sailthrough
     (abs(cms2.mus_mc_id()[i_mu]) != 13 )
     ) {
    category = 1;
  }
  else if( 
          abs(cms2.mus_mc_id()[i_mu]) == 13 && 
          abs(cms2.mus_mc_motherid()[i_mu]) < 400 
          ) {
    // light hadrons
    category = 2;
  }
  else if(
          ( abs(cms2.mus_mc_id()[i_mu]) == 13          &&
            abs(cms2.mus_mc_motherid()[i_mu]) >=400   &&
            abs(cms2.mus_mc_motherid()[i_mu]) <=600 ) || 
          ( abs(cms2.mus_mc_id()[i_mu]) == 13          &&
            abs(cms2.mus_mc_motherid()[i_mu]) >=4000   &&
            abs(cms2.mus_mc_motherid()[i_mu]) <=6000 )
          ) {
    // heavy hadrons
    category = 3;
  }
  else {
    // the rest
    category = 4;
  }
  return category;
}

bool idIsCharm(int id) {
  id = abs(id);
  if (
      id == 4       ||
      id == 411     ||
      id == 421     ||
      id == 10411   ||
      id == 10421   ||
      id == 413     ||
      id == 423     ||
      id == 10413   ||
      id == 10423   ||
      id == 20413   ||
      id == 20423   ||
      id == 415     ||
      id == 425     ||
      id == 431     ||
      id == 10431   ||
      id == 433     ||
      id == 10433   ||
      id == 20433   ||
      id == 435     ||
      id == 441     ||
      id == 10441   ||
      id == 100441  ||
      id == 443     ||
      id == 10443   ||
      id == 20443   ||
      id == 100443  ||
      id == 30443   ||
      id == 9000443 ||
      id == 9010443 ||
      id == 9020443 ||
      id == 445     ||
      id == 9000445 ||
      id == 4122    ||
      id == 4222    ||
      id == 4212    ||
      id == 4112    ||
      id == 4224    ||
      id == 4214    ||
      id == 4114    ||
      id == 4232    ||
      id == 4132    ||
      id == 4322    ||
      id == 4312    ||
      id == 4324    ||
      id == 4314    ||
      id == 4332    ||
      id == 4334    ||
      id == 4412    ||
      id == 4422    ||
      id == 4414    ||
      id == 4424    ||
      id == 4432    ||
      id == 4434    ||
      id == 4444
      ) {
    return true;
  }
  else return false;
}

bool idIsBeauty(int id) {
  id = abs(id);
  if (
      id == 5       ||
      id == 511     ||
      id == 521     ||
      id == 10511   ||
      id == 10521   ||
      id == 513     ||
      id == 523     ||
      id == 10513   ||
      id == 10523   ||
      id == 20513   ||
      id == 20523   ||
      id == 515     ||
      id == 525     ||
      id == 531     ||
      id == 10531   ||
      id == 533     ||
      id == 10533   ||
      id == 20533   ||
      id == 535     ||
      id == 541     ||
      id == 10541   ||
      id == 543     ||
      id == 10543   ||
      id == 20543   ||
      id == 545     ||
      id == 551     ||
      id == 10551   ||
      id == 100551  ||
      id == 110551  ||
      id == 200551  ||
      id == 210551  ||
      id == 553     ||
      id == 10553   ||
      id == 20553   ||
      id == 30553   ||
      id == 100553  ||
      id == 110553  ||
      id == 120553  ||
      id == 130553  ||
      id == 200553  ||
      id == 210553  ||
      id == 220553  ||
      id == 300553  ||
      id == 9000553 ||
      id == 9010553 ||
      id == 555     ||
      id == 10555   ||
      id == 20555   ||
      id == 100555  ||
      id == 110555  ||
      id == 120555  ||
      id == 200555  ||
      id == 557     ||
      id == 100557  ||
      id == 5122    || 
      id == 5112    ||
      id == 5212    ||
      id == 5222    ||
      id == 5114    ||
      id == 5214    ||
      id == 5224    ||
      id == 5132    ||
      id == 5232    ||
      id == 5312    ||
      id == 5322    ||
      id == 5314    ||
      id == 5324    ||
      id == 5332    ||
      id == 5334    ||
      id == 5142    ||
      id == 5242    ||
      id == 5412    ||
      id == 5422    ||
      id == 5414    ||
      id == 5424    ||
      id == 5342    ||
      id == 5432    ||
      id == 5434    ||
      id == 5442    ||
      id == 5444    ||
      id == 5512    ||
      id == 5522    ||
      id == 5514    ||
      id == 5524    ||
      id == 5532    ||
      id == 5534    ||
      id == 5542    ||
      id == 5544    ||
      id == 5554 
      ) {
    return true;
  }
  else return false;
}

// -------------------------------------
// quick and dirty function to see if
// lepton is prompt (approximately)
// 
// Note: As I re-read this with a fresher
// eye this certainly isn't perfect as 
// there will be real leptons I miss
// as well as fake leptons I pick up when
// the matching messes up.  I suppose there
// are a couple different ways to implement
// this equivalently or more better. 
// -------------------------------------
bool isNotPromptSusyLeptonHyp(int idx) {

	 // require lepton matched to lepton of same flavor at status 1
	 if (abs(cms2.hyp_lt_id()[idx]) == 11 && abs(cms2.els_mc_id()[cms2.hyp_lt_index()[idx]]) != 11)
		  return true;
	 if (abs(cms2.hyp_ll_id()[idx]) == 11 && abs(cms2.els_mc_id()[cms2.hyp_ll_index()[idx]]) != 11)
		  return true;
	 if (abs(cms2.hyp_lt_id()[idx]) == 13 && abs(cms2.mus_mc_id()[cms2.hyp_lt_index()[idx]]) != 13)
		  return true;
	 if (abs(cms2.hyp_ll_id()[idx]) == 13 && abs(cms2.mus_mc_id()[cms2.hyp_ll_index()[idx]]) != 13)
		  return true;

	 // require status 3 mother to be a SUSY particle or a W or Z
	 if (abs(cms2.hyp_lt_id()[idx]) == 11)
		  if (abs(cms2.els_mc3_motherid()[cms2.hyp_lt_index()[idx]]) < 1000000 && abs(cms2.els_mc3_motherid()[cms2.hyp_lt_index()[idx]]) != 24 && abs(cms2.els_mc3_motherid()[cms2.hyp_lt_index()[idx]]) != 23)
			   return true;						 
	 if (abs(cms2.hyp_ll_id()[idx]) == 11)
		  if (abs(cms2.els_mc3_motherid()[cms2.hyp_ll_index()[idx]]) < 1000000 && abs(cms2.els_mc3_motherid()[cms2.hyp_ll_index()[idx]]) != 24 && abs(cms2.els_mc3_motherid()[cms2.hyp_ll_index()[idx]]) != 23)
			   return true;
	 if (abs(cms2.hyp_lt_id()[idx]) == 13)
		  if (abs(cms2.mus_mc3_motherid()[cms2.hyp_lt_index()[idx]]) < 1000000 && abs(cms2.mus_mc3_motherid()[cms2.hyp_lt_index()[idx]]) != 24 && abs(cms2.mus_mc3_motherid()[cms2.hyp_lt_index()[idx]]) != 23)
			   return true;
	 if (abs(cms2.hyp_ll_id()[idx]) == 13)
		  if (abs(cms2.mus_mc3_motherid()[cms2.hyp_ll_index()[idx]]) < 1000000 && abs(cms2.mus_mc3_motherid()[cms2.hyp_ll_index()[idx]]) != 24 && abs(cms2.mus_mc3_motherid()[cms2.hyp_ll_index()[idx]]) != 23)
			   return true;

	 // require status 1 mother to be a SUSY particle, W, Z or tau
	 if (abs(cms2.hyp_lt_id()[idx]) == 11)
		  if (abs(cms2.els_mc_motherid()[cms2.hyp_lt_index()[idx]]) < 1000000 && abs(cms2.els_mc_motherid()[cms2.hyp_lt_index()[idx]]) != 24 && abs(cms2.els_mc_motherid()[cms2.hyp_lt_index()[idx]]) != 15 && abs(cms2.els_mc_motherid()[cms2.hyp_lt_index()[idx]]) != 23)
			   return true;						 
	 if (abs(cms2.hyp_ll_id()[idx]) == 11)
		  if (abs(cms2.els_mc_motherid()[cms2.hyp_ll_index()[idx]]) < 1000000 && abs(cms2.els_mc_motherid()[cms2.hyp_ll_index()[idx]]) != 24 && abs(cms2.els_mc_motherid()[cms2.hyp_ll_index()[idx]]) != 15 && abs(cms2.els_mc_motherid()[cms2.hyp_ll_index()[idx]]) != 23)
			   return true;
	 if (abs(cms2.hyp_lt_id()[idx]) == 13)
		  if (abs(cms2.mus_mc_motherid()[cms2.hyp_lt_index()[idx]]) < 1000000 && abs(cms2.mus_mc_motherid()[cms2.hyp_lt_index()[idx]]) != 24 && abs(cms2.mus_mc_motherid()[cms2.hyp_lt_index()[idx]]) != 15 && abs(cms2.mus_mc_motherid()[cms2.hyp_lt_index()[idx]]) != 23)
			   return true;
	 if (abs(cms2.hyp_ll_id()[idx]) == 13)
		  if (abs(cms2.mus_mc_motherid()[cms2.hyp_ll_index()[idx]]) < 1000000 && abs(cms2.mus_mc_motherid()[cms2.hyp_ll_index()[idx]]) != 24 && abs(cms2.mus_mc_motherid()[cms2.hyp_ll_index()[idx]]) != 15 && abs(cms2.mus_mc_motherid()[cms2.hyp_ll_index()[idx]]) != 23)
			   return true;


	 return false;
}


int mc3idx_eormu(int id, int idx, float maxDr, float minPt){
  LorentzVector lepp4 =  abs(id)==11 ? cms2.els_p4()[idx] : cms2.mus_p4()[idx];
  
  float dRMin = 999999;
  int nGens = cms2.genps_id().size();
  int genidx = -1;
  for (int iG = 0; iG < nGens; ++iG){
    if (cms2.genps_status()[iG]==3){
      float dr = ROOT::Math::VectorUtil::DeltaR(lepp4, cms2.genps_p4()[iG]);
      if (dr < maxDr && cms2.genps_p4()[iG].pt() > minPt && dr < dRMin){
	genidx = iG;
	dRMin = dr;
      }
    }
  }
  return genidx;
}

float mc3dr_eormu(int id, int idx, float maxDr, float minPt){
  int genidx = mc3idx_eormu(id,idx,maxDr,minPt);
  if (genidx < 0) return 999999;
  LorentzVector lepp4 =  abs(id)==11 ? cms2.els_p4()[idx] : cms2.mus_p4()[idx];
  return ROOT::Math::VectorUtil::DeltaR(lepp4, cms2.genps_p4()[genidx]);
}
