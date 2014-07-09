#include <assert.h>
#include <algorithm>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "triggerUtils.h"
#include "CMS2.h"

#include "susySelections.h"
#include "triggerUtils.h"
#include "electronSelections.h"
#include "muonSelections.h"
#include "eventSelections.h"
#include "mcSelections.h"

using namespace tas;
using namespace wp2012;

//--------------------------------------------------------------------------------
// loose electron WP of:
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification
//--------------------------------------------------------------------------------

bool passElectronSelection_ZMet2012_v1_NoIso(int index, bool vetoTransition, bool eta24){

  if( vetoTransition && fabs(cms2.els_etaSC()[index]) > 1.4442 && fabs(cms2.els_etaSC()[index]) < 1.566 ) return false;
  if( eta24 && fabs(cms2.els_p4()[index].eta()) > 2.4 )                                                   return false;


  electronIdComponent_t answer_loose_2012 = electronId_WP2012(index, LOOSE);
  if ((answer_loose_2012 & PassWP2012CutsNoIso) == PassWP2012CutsNoIso) return true;
  
  return false;
}

bool passElectronSelection_ZMet2012_v1_Iso(int index, bool vetoTransition, bool eta24){

  if( vetoTransition && fabs(cms2.els_etaSC()[index]) > 1.4442 && fabs(cms2.els_etaSC()[index]) < 1.566 ) return false;
  if( eta24 && fabs(cms2.els_p4()[index].eta()) > 2.4 )                                                   return false;

  electronIdComponent_t answer_loose_2012 = electronId_WP2012(index, LOOSE);
  if ((answer_loose_2012 & PassWP2012CutsIso) == PassWP2012CutsIso) return true;
  
  return false;
}

bool passElectronSelection_ZMet2012_v1(int index, bool vetoTransition, bool eta24){

  if( vetoTransition && fabs(cms2.els_etaSC()[index]) > 1.4442 && fabs(cms2.els_etaSC()[index]) < 1.566 ) return false;
  if( eta24 && fabs(cms2.els_p4()[index].eta()) > 2.4 )                                                   return false;

  electronIdComponent_t answer_loose_2012 = electronId_WP2012(index, LOOSE);
  if ((answer_loose_2012 & PassAllWP2012Cuts) == PassAllWP2012Cuts)  return true;
  
  return false;
}

//------------------------------------------------------------------------------------------------
// loose electron WP of:
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification
// v2: update electron isolation branches and use d0/dz of GSF track w.r.t. 1st GOOD vertex
//------------------------------------------------------------------------------------------------

bool overlapMuon_ZMet2012_v1(int index , float ptcut = 10.0 ){

  for( unsigned int imu = 0 ; imu < mus_p4().size(); ++imu ){

    float dr = ROOT::Math::VectorUtil::DeltaR( cms2.els_p4().at(index) , cms2.mus_p4().at(imu) );
    
    if( dr > 0.1                           ) continue;
    if( cms2.mus_p4().at(imu).pt() < ptcut ) continue;
    if( !muonId( imu , ZMet2012_v1 )       ) continue;
    
    return true;
  }

  return false;

}

bool overlapMuon_ZMet2012_v2(int index , float ptcut = 10.0 ){

  for( unsigned int imu = 0 ; imu < mus_p4().size(); ++imu ){

    float dr = ROOT::Math::VectorUtil::DeltaR( cms2.els_p4().at(index) , cms2.mus_p4().at(imu) );
    
    if( dr > 0.1                           ) continue;
    if( cms2.mus_p4().at(imu).pt() < ptcut ) continue;
    if ( (((cms2.mus_type().at(imu)) & (1<<1)) == 0) && (((cms2.mus_type().at(imu)) & (1<<2)) == 0) ) continue;
    
    return true;
  }

  return false;

}

bool passElectronSelection_ZMet2012_v2_NoIso(int index, bool vetoTransition, bool eta24, bool useOldIsolation ){

  if( vetoTransition && fabs(cms2.els_etaSC()[index]) > 1.4442 && fabs(cms2.els_etaSC()[index]) < 1.566 ) return false;
  if( eta24 && fabs(cms2.els_p4()[index].eta()) > 2.4 )                                                   return false;
  if( overlapMuon_ZMet2012_v1(index,10.0) )                                                               return false;

  electronIdComponent_t answer_loose_2012 = electronId_WP2012_v2(index, LOOSE, useOldIsolation);
  if ((answer_loose_2012 & PassWP2012CutsNoIso) == PassWP2012CutsNoIso) return true;
  
  return false;
}

bool passElectronSelection_ZMet2012_v2_Iso(int index, bool vetoTransition, bool eta24, bool useOldIsolation ){

  if( vetoTransition && fabs(cms2.els_etaSC()[index]) > 1.4442 && fabs(cms2.els_etaSC()[index]) < 1.566 ) return false;
  if( eta24 && fabs(cms2.els_p4()[index].eta()) > 2.4 )                                                   return false;
  if( overlapMuon_ZMet2012_v1(index,10.0) )                                                               return false;

  electronIdComponent_t answer_loose_2012 = electronId_WP2012_v2(index, LOOSE, useOldIsolation);
  if ((answer_loose_2012 & PassWP2012CutsIso) == PassWP2012CutsIso) return true;
  
  return false;
}

bool passElectronSelection_ZMet2012_v2(int index, bool vetoTransition, bool eta24, bool useOldIsolation ){

  if( vetoTransition && fabs(cms2.els_etaSC()[index]) > 1.4442 && fabs(cms2.els_etaSC()[index]) < 1.566 ) return false;
  if( eta24 && fabs(cms2.els_p4()[index].eta()) > 2.4 )                                                   return false;
  if( overlapMuon_ZMet2012_v1(index,10.0) )                                                               return false;

  electronIdComponent_t answer_loose_2012 = electronId_WP2012_v2(index, LOOSE, useOldIsolation);
  if ((answer_loose_2012 & PassAllWP2012Cuts) == PassAllWP2012Cuts)  return true;
  
  return false;
}

bool passElectronSelection_ZMet2012_v2_DetIso(int index, bool vetoTransition, bool eta24, bool useOldIsolation ){

  if( vetoTransition && fabs(cms2.els_etaSC()[index]) > 1.4442 && fabs(cms2.els_etaSC()[index]) < 1.566 ) return false;
  if( eta24 && fabs(cms2.els_p4()[index].eta()) > 2.4 )                                                   return false;
  if( overlapMuon_ZMet2012_v1(index,10.0) )                                                               return false;
  if( electronIsolation_rel_v1(index, true ) > 0.15 )                                                     return false;              

  electronIdComponent_t answer_loose_2012 = electronId_WP2012_v2(index, LOOSE, useOldIsolation);
  if ((answer_loose_2012 & PassWP2012CutsNoIso) == PassWP2012CutsNoIso) return true;
  
  return false;
}

//------------------------------------------------------------------------------------------------
// loose electron WP of:
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification
// v3: same as v2, but update electron effective areas
//------------------------------------------------------------------------------------------------

bool passElectronSelection_ZMet2012_v3_NoIso(int index, bool vetoTransition, bool eta24, bool useOldIsolation ){

  if( vetoTransition && fabs(cms2.els_etaSC()[index]) > 1.4442 && fabs(cms2.els_etaSC()[index]) < 1.566 ) return false;
  if( eta24 && fabs(cms2.els_p4()[index].eta()) > 2.4 )                                                   return false;
  if( overlapMuon_ZMet2012_v1(index,10.0) )                                                               return false;

  electronIdComponent_t answer_loose_2012 = electronId_WP2012_v3(index, LOOSE, useOldIsolation);
  if ((answer_loose_2012 & PassWP2012CutsNoIso) == PassWP2012CutsNoIso) return true;
  
  return false;
}

bool passElectronSelection_ZMet2012_v3_Iso(int index, bool vetoTransition, bool eta24, bool useOldIsolation ){

  if( vetoTransition && fabs(cms2.els_etaSC()[index]) > 1.4442 && fabs(cms2.els_etaSC()[index]) < 1.566 ) return false;
  if( eta24 && fabs(cms2.els_p4()[index].eta()) > 2.4 )                                                   return false;
  if( overlapMuon_ZMet2012_v1(index,10.0) )                                                               return false;

  electronIdComponent_t answer_loose_2012 = electronId_WP2012_v3(index, LOOSE, useOldIsolation);
  if ((answer_loose_2012 & PassWP2012CutsIso) == PassWP2012CutsIso) return true;
  
  return false;
}

bool passElectronSelection_ZMet2012_v3(int index, bool vetoTransition, bool eta24, bool useOldIsolation ){

  if( vetoTransition && fabs(cms2.els_etaSC()[index]) > 1.4442 && fabs(cms2.els_etaSC()[index]) < 1.566 ) return false;
  if( eta24 && fabs(cms2.els_p4()[index].eta()) > 2.4 )                                                   return false;
  if( overlapMuon_ZMet2012_v1(index,10.0) )                                                               return false;

  electronIdComponent_t answer_loose_2012 = electronId_WP2012_v3(index, LOOSE, useOldIsolation);
  if ((answer_loose_2012 & PassAllWP2012Cuts) == PassAllWP2012Cuts)  return true;
  
  return false;
}

//--------------------------------------------------------------------------------
// medium electron WP of:
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification
//--------------------------------------------------------------------------------

bool passElectronSelection_Stop2012_v1_NoIso(int index, bool vetoTransition, bool eta24){

  if( vetoTransition && fabs(cms2.els_etaSC()[index]) > 1.4442 && fabs(cms2.els_etaSC()[index]) < 1.566 ) return false;
  if( eta24 && fabs(cms2.els_p4()[index].eta()) > 2.4 )                                                   return false;


  electronIdComponent_t answer_loose_2012 = electronId_WP2012(index, MEDIUM);
  if ((answer_loose_2012 & PassWP2012CutsNoIso) == PassWP2012CutsNoIso) return true;
  
  return false;
}

bool passElectronSelection_Stop2012_v1_Iso(int index, bool vetoTransition, bool eta24){

  if( vetoTransition && fabs(cms2.els_etaSC()[index]) > 1.4442 && fabs(cms2.els_etaSC()[index]) < 1.566 ) return false;
  if( eta24 && fabs(cms2.els_p4()[index].eta()) > 2.4 )                                                   return false;

  electronIdComponent_t answer_loose_2012 = electronId_WP2012(index, MEDIUM);
  if ((answer_loose_2012 & PassWP2012CutsIso) == PassWP2012CutsIso) return true;
  
  return false;
}

bool passElectronSelection_Stop2012_v1(int index, bool vetoTransition, bool eta24){

  if( vetoTransition && fabs(cms2.els_etaSC()[index]) > 1.4442 && fabs(cms2.els_etaSC()[index]) < 1.566 ) return false;
  if( eta24 && fabs(cms2.els_p4()[index].eta()) > 2.4 )                                                   return false;

  electronIdComponent_t answer_loose_2012 = electronId_WP2012(index, MEDIUM);
  if ((answer_loose_2012 & PassAllWP2012Cuts) == PassAllWP2012Cuts)  return true;
  
  return false;
}

//------------------------------------------------------------------------------------------------
// medium electron WP of:
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification
// v2: update electron isolation branches and use d0/dz of GSF track w.r.t. 1st GOOD vertex
//------------------------------------------------------------------------------------------------

bool passElectronSelection_Stop2012_v2_NoIso(int index, bool vetoTransition, bool eta24, bool useOldIsolation ){

  if( vetoTransition && fabs(cms2.els_etaSC()[index]) > 1.4442 && fabs(cms2.els_etaSC()[index]) < 1.566 ) return false;
  if( eta24 && fabs(cms2.els_p4()[index].eta()) > 2.4 )                                                   return false;
  if( overlapMuon_ZMet2012_v1(index,10.0) )                                                               return false;

  electronIdComponent_t answer_loose_2012 = electronId_WP2012_v2(index, MEDIUM, useOldIsolation);
  if ((answer_loose_2012 & PassWP2012CutsNoIso) == PassWP2012CutsNoIso) return true;
  
  return false;
}

bool passElectronSelection_Stop2012_v2_Iso(int index, bool vetoTransition, bool eta24, bool useOldIsolation ){

  if( vetoTransition && fabs(cms2.els_etaSC()[index]) > 1.4442 && fabs(cms2.els_etaSC()[index]) < 1.566 ) return false;
  if( eta24 && fabs(cms2.els_p4()[index].eta()) > 2.4 )                                                   return false;
  if( overlapMuon_ZMet2012_v1(index,10.0) )                                                               return false;

  electronIdComponent_t answer_loose_2012 = electronId_WP2012_v2(index, MEDIUM, useOldIsolation);
  if ((answer_loose_2012 & PassWP2012CutsIso) == PassWP2012CutsIso) return true;
  
  return false;
}

bool passElectronSelection_Stop2012_v2(int index, bool vetoTransition, bool eta24, bool useOldIsolation ){

  if( vetoTransition && fabs(cms2.els_etaSC()[index]) > 1.4442 && fabs(cms2.els_etaSC()[index]) < 1.566 ) return false;
  if( eta24 && fabs(cms2.els_p4()[index].eta()) > 2.4 )                                                   return false;
  if( overlapMuon_ZMet2012_v1(index,10.0) )                                                               return false;

  electronIdComponent_t answer_loose_2012 = electronId_WP2012_v2(index, MEDIUM, useOldIsolation);
  if ((answer_loose_2012 & PassAllWP2012Cuts) == PassAllWP2012Cuts)  return true;
  
  return false;
}

//------------------------------------------------------------------------------------------------
// medium electron WP of:
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification
// v3: relax selection applied to muons used for overlap removal and update effective areas
//------------------------------------------------------------------------------------------------

bool passElectronSelection_Stop2012_v3_NoIso(int index, bool vetoTransition, bool eta24, bool useOldIsolation ){

  if( vetoTransition && fabs(cms2.els_etaSC()[index]) > 1.4442 && fabs(cms2.els_etaSC()[index]) < 1.566 ) return false;
  if( eta24 && fabs(cms2.els_p4()[index].eta()) > 2.4 )                                                   return false;
  if( overlapMuon_ZMet2012_v2(index,10.0) )                                                               return false;

  electronIdComponent_t answer_loose_2012 = electronId_WP2012_v3(index, MEDIUM, useOldIsolation);
  if ((answer_loose_2012 & PassWP2012CutsNoIso) == PassWP2012CutsNoIso) return true;
  
  return false;
}

bool passElectronSelection_Stop2012_v3_Iso(int index, bool vetoTransition, bool eta24, bool useOldIsolation ){

  if( vetoTransition && fabs(cms2.els_etaSC()[index]) > 1.4442 && fabs(cms2.els_etaSC()[index]) < 1.566 ) return false;
  if( eta24 && fabs(cms2.els_p4()[index].eta()) > 2.4 )                                                   return false;
  if( overlapMuon_ZMet2012_v2(index,10.0) )                                                               return false;

  electronIdComponent_t answer_loose_2012 = electronId_WP2012_v3(index, MEDIUM, useOldIsolation);
  if ((answer_loose_2012 & PassWP2012CutsIso) == PassWP2012CutsIso) return true;
  
  return false;
}

bool passElectronSelection_Stop2012_v3(int index, bool vetoTransition, bool eta24, bool useOldIsolation ){

  if( vetoTransition && fabs(cms2.els_etaSC()[index]) > 1.4442 && fabs(cms2.els_etaSC()[index]) < 1.566 ) return false;
  if( eta24 && fabs(cms2.els_p4()[index].eta()) > 2.4 )                                                   return false;
  if( overlapMuon_ZMet2012_v2(index,10.0) )                                                               return false;

  electronIdComponent_t answer_loose_2012 = electronId_WP2012_v3(index, MEDIUM, useOldIsolation);
  if ((answer_loose_2012 & PassAllWP2012Cuts) == PassAllWP2012Cuts)  return true;
  
  return false;
}



//--------------------------------------------------------
// leptonOrTauIsFromW(int idx, int id, bool alsoSusy)
// this function is the same as leptonIsFromW in mcSelections
// except that it distinguishes between W->e/mu (result 1 or 2)
// vs. W->tau->e/mu (result 3 or 4)
//--------------------------------------------------------

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
// Output:  4 = from W/Z->tau->l incorrect charge  // <<<--- Hooberman 04/29/11: distinguish 
//          3 = from W/Z->tau->l   correct charge  // <<<--- W->l (1 or 2) vs. W->tau->l (3 or 4)
//          2 = from W/Z incorrect charge
//          1 = from W/Z   correct charge
//          0 = not matched to a lepton (= fake)
//         -1 = lepton from b decay
//         -2 = lepton from c decay
//         -3 = lepton from some other source
//
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
int leptonOrTauIsFromW(int idx, int id, bool alsoSusy) {

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
  //std::cout << "id=" << id << " st1_id=" << st1_id;
  //std::cout << " st3_id=" << st3_id;
  //std::cout << " st1_motherid=" << st1_motherid << std::endl;

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
    //if(id*st3_id > 0) 
    //  return 1;       
    //else return 2;
    if(id*st3_id > 0) 
      return 3;       // <<<--- Hooberman 04/29/11: distinguish btw W->l (1 or 2) vs. W->tau->l (3 or 4)
    else return 4;
  }
  

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


/*****************************************************************************************/
//print event info
/*****************************************************************************************/
void printEventInfo(){
  //cout << cms2.evt_dataset() << endl;
  cout << cms2.evt_run() << " " << cms2.evt_lumiBlock() << " " << cms2.evt_event() << endl;
}


/*****************************************************************************************/
//veto Z->mumugamma events                                                                
/*****************************************************************************************/
bool vetoZmumuGamma( unsigned int hypIdx , float emax , float minmass , float maxmass ){
  
  //we only care about mumu hyp types
  if( cms2.hyp_type().at(hypIdx) != 0 ) return false; 

  //we only want to veto events that started *outside* of Z window 
  //and moved inside after including ecal deposit
  if( cms2.hyp_p4().at(hypIdx).mass() > minmass && cms2.hyp_p4().at(hypIdx).mass() < maxmass )
    return false;

  //get ecal deposits for each muon
  int ill = cms2.hyp_ll_index().at( hypIdx );
  int ilt = cms2.hyp_lt_index().at( hypIdx );

  //energy = ET * cosh(eta)
  float ell = cms2.mus_iso_ecalvetoDep().at(ill) * cosh( cms2.mus_p4().at(ill).eta() ) ;
  float elt = cms2.mus_iso_ecalvetoDep().at(ilt) * cosh( cms2.mus_p4().at(ilt).eta() ) ;

  //don't veto event if ecal deposits are both less than emax
  if( ell < emax && elt < emax ) return false;

  LorentzVector vll = cms2.hyp_ll_p4().at(hypIdx);
  LorentzVector vlt = cms2.hyp_lt_p4().at(hypIdx);

  //if ecal deposit is greater than emax, scale muon momentum up
  if( ell > emax ) vll = vll * ( 1 + ell / vll.energy() );
  if( elt > emax ) vlt = vlt * ( 1 + elt / vlt.energy() );

  float mass = ( vll + vlt ).mass();

  //check if dimuon mass, including extra ecal energy, is inside Z mass window
  if( mass > minmass && mass < maxmass ){
    //cout << "Veto Z->mumugamma event!" << endl;
    //cout << cms2.evt_run() << " " << cms2.evt_lumiBlock() << " " << cms2.evt_event() << endl;
    return true;
  }

  return false;
}

/*****************************************************************************************/
//generalized Z veto
/*****************************************************************************************/
bool ZVetoGeneral( float ptcut , float minmass ,  float maxmass , SelectionType mutype ){

  for(unsigned int i = 0; i < hyp_p4().size(); ++i) {

    //check that hyp leptons come from same vertex
    if( !hypsFromSameVtx( i ) )    continue;

    //opposite-sign, same-flavor
    if( cms2.hyp_lt_id().at(i)     *  cms2.hyp_ll_id().at(i) > 0 ) continue;
    if( cms2.hyp_type().at(i) == 1 || cms2.hyp_type().at(i) == 2 ) continue;

    //check that both lepton pt's > ptcut
    if( cms2.hyp_ll_p4().at(i).pt() < ptcut ) continue; 
    if( cms2.hyp_lt_p4().at(i).pt() < ptcut ) continue; 
          
    //muon ID
    if ( abs( cms2.hyp_ll_id().at(i) ) == 13  && !( muonId( cms2.hyp_ll_index().at(i) , mutype ) ) )   continue;
    if ( abs( cms2.hyp_lt_id().at(i) ) == 13  && !( muonId( cms2.hyp_lt_index().at(i) , mutype ) ) )   continue;
          
    //electron ID
    if ( abs( cms2.hyp_ll_id().at(i) ) == 11  && !( pass_electronSelection( cms2.hyp_ll_index().at(i) , electronSelection_el_OSV3 , false , false ))) continue;
    if ( abs( cms2.hyp_lt_id().at(i) ) == 11  && !( pass_electronSelection( cms2.hyp_lt_index().at(i) , electronSelection_el_OSV3 , false , false ))) continue;
          
    if( cms2.hyp_p4().at(i).mass() > minmass && cms2.hyp_p4().at(i).mass() < maxmass ){
      //cout << "General Z veto: mass " << cms2.hyp_p4().at(i).mass() << endl;
      //printEventInfo();
      return true; 
    }
  }

  return false;
}


//-------------------------------------------------------
// get exact trigger name corresponding to given pattern
//-------------------------------------------------------
TString triggerName(TString triggerPattern){

  bool    foundTrigger  = false;
  TString exact_hltname = "";

  for( unsigned int itrig = 0 ; itrig < hlt_trigNames().size() ; ++itrig ){
    if( TString( hlt_trigNames().at(itrig) ).Contains( triggerPattern ) ){
      foundTrigger  = true;
      exact_hltname = hlt_trigNames().at(itrig);
      break;
    }
  }

  if( !foundTrigger) return "TRIGGER_NOT_FOUND";

  return exact_hltname;

}


//---------------------------------------------
// Check if trigger passes
//---------------------------------------------

bool passHLTTriggerPattern(const char* arg){

  TString HLTTriggerPattern(arg);
  TString HLTTrigger = triggerName( HLTTriggerPattern );

  if( HLTTrigger.Contains("TRIGGER_NOT_FOUND")){
    return false;
  }
  return passHLTTrigger( HLTTrigger );
}


//---------------------------------------------
// Check if trigger is unprescaled and passes
//---------------------------------------------
bool passUnprescaledHLTTriggerPattern(const char* arg){

  //cout << "Checking for pattern " << arg << endl;

  TString HLTTriggerPattern(arg);

  TString HLTTrigger = triggerName( HLTTriggerPattern );

  //cout << "Found trigger " << HLTTrigger << endl;

  if( HLTTrigger.Contains("TRIGGER_NOT_FOUND")){
    //cout << "Didn't find trigger!" << endl;
    return false;
  }
  return passUnprescaledHLTTrigger( HLTTrigger );

}

//---------------------------------------------
// function returns:
// -1: no matching trigger found
//  0: trigger not passed
//  1: trigger passed, un-prescaled
//  N: trigger passed, prescale N
//---------------------------------------------

int passTriggerPrescale(const char* arg){

  //Find exact trigger name
  TString HLTTriggerPattern(arg);
  TString HLTTrigger = triggerName( HLTTriggerPattern );

  //Return -1 if no matching trigger found
  if( HLTTrigger.Contains("TRIGGER_NOT_FOUND") )  return -1;
 
  //Return 0 if trigger didn't pass
  if( !passHLTTrigger( HLTTrigger ) ) return 0;

  //Return 1 if trigger passes and is unprescaled
  if( passUnprescaledHLTTrigger( HLTTrigger ) ) return 1;

  //Return prescale if prescaled trigger passes
  return HLT_prescale( HLTTrigger );

}

//---------------------------------------------
// single muon triggers for lljj bump search
//---------------------------------------------

bool passMuMuJJTrigger_v1( bool isData ) {

  if( isData ){
    
    //-----------------------------------------------------------------------------
    if (evt_run() >= 160329 && evt_run() <= 163261){
      if( passUnprescaledHLTTriggerPattern("HLT_IsoMu15_v5") )          return true;
    }
    //-----------------------------------------------------------------------------
    else if (evt_run() >= 163269 && evt_run() <= 164236){
      if( passUnprescaledHLTTriggerPattern("HLT_IsoMu24_v2") )          return true;
    }
    //-----------------------------------------------------------------------------
    else if (evt_run() >= 165088 && evt_run() <= 165887){
      if( passUnprescaledHLTTriggerPattern("HLT_IsoMu24_v4") )          return true;
    }
    //-----------------------------------------------------------------------------
    else if (evt_run() == 166346 ){
      if( passUnprescaledHLTTriggerPattern("HLT_IsoMu24_v6") )          return true;
    }
    //-----------------------------------------------------------------------------
    else if (evt_run() >= 165922 && evt_run() <= 167043){
      if( passUnprescaledHLTTriggerPattern("HLT_IsoMu24_v5") )          return true;
    }
    //-----------------------------------------------------------------------------
    else if (evt_run() >= 167078 && evt_run() <= 170053){
      if( passUnprescaledHLTTriggerPattern("HLT_IsoMu24_v7") )          return true;
    }
    //-----------------------------------------------------------------------------
    else if (evt_run() >= 170071 && evt_run() <= 173198){
      if( passUnprescaledHLTTriggerPattern("HLT_IsoMu24_v8") )          return true;
    }
    //-----------------------------------------------------------------------------
    else if (evt_run() >= 173212 && evt_run() <= 178380){
      if( passUnprescaledHLTTriggerPattern("HLT_IsoMu30_eta2p1_v3") )   return true;
    }
    //-----------------------------------------------------------------------------
    else if (evt_run() >= 178420 && evt_run() <= 179889){
      if( passUnprescaledHLTTriggerPattern("HLT_IsoMu30_eta2p1_v6") )   return true;
    }
    //-----------------------------------------------------------------------------
    else if (evt_run() >= 179959 && evt_run() <= 180093){
      if( passUnprescaledHLTTriggerPattern("HLT_IsoMu30_eta2p1_v7") )   return true;
    }
    //-----------------------------------------------------------------------------
  }

  else{
    if( passUnprescaledHLTTriggerPattern("HLT_IsoMu24_v") )  return true;
  }

  return false;
}


//-----------------------------------------------
// single lepton + jets triggers for stop search
//-----------------------------------------------

bool passSingleLepSUSYTrigger2011_v1( bool isData , int lepType ) {

  //------------------------------------------------------------------------
  // These are the all the triggers considered for single lepton+jets 
  //------------------------------------------------------------------------

  // no triggers required for MC
  if( !isData ) return true;

  if( passSingleLep2JetSUSYTrigger2011( isData , lepType ) )   return true; // l+dijet+MHT triggers
  if( passSingleLep3JetSUSYTrigger2011( isData , lepType ) )   return true; // l+trijet triggers
  if( passSingleMuTrigger2011(          isData , lepType ) )   return true; // single muon triggers

  return false;
}

bool passSingleLep2JetSUSYTrigger2011( bool isData , int lepType ) {

  //-------------------------------------------------------
  // These are the trigger options for lepton+2jets+MET 
  //-------------------------------------------------------

  // no triggers required for MC
  if( !isData ) return true;

  // electron channel
  if( lepType == 0 ){
    if( passUnprescaledHLTTriggerPattern("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v") )                                   return true; // 160329-164236
    if( passUnprescaledHLTTriggerPattern("HLT_Ele17_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_CentralJet25_PFMHT15_v") ) return true; // 165088-165887
    if( passUnprescaledHLTTriggerPattern("HLT_Ele22_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_CentralJet25_PFMHT20_v") ) return true; // 166979-173198
    if( passUnprescaledHLTTriggerPattern("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_CentralJet25_PFMHT20_v") ) return true; // 170826-176309
    if( passUnprescaledHLTTriggerPattern("HLT_Ele30_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_DiCentralJet30_PFMHT25_v") )            return true; // 173212-178380
    if( passUnprescaledHLTTriggerPattern("HLT_Ele27_WP80_DiCentralPFJet25_PFMHT15_v") )                                      return true; // 178420-180291
  }

  // muon channel
  else if( lepType == 1 ){    
    if( passUnprescaledHLTTriggerPattern("HLT_IsoMu17_v") )          return true; // 160329-165887
    if( passUnprescaledHLTTriggerPattern("HLT_IsoMu24_v") )          return true; // 160329-173198
    if( passUnprescaledHLTTriggerPattern("HLT_IsoMu30_eta2p1_v") )   return true; // 173212-180291
  }

  else{
    cout << "susySelections.cc:: ERROR unrecognized lepType " << lepType << ", quitting" << endl;
    exit(0);
  }

  return false;
}

bool passSingleMuTrigger2011( bool isData , int lepType ) {
  
  //----------------------------
  // single muon triggers
  //----------------------------

  // no triggers required for MC
  if( !isData ) return true;

  // false for electron channel
  if( lepType == 0 ){
    return false;
  }

  // muon channel
  else if( lepType == 1 ){    
    if( passUnprescaledHLTTriggerPattern("HLT_IsoMu17_v") )          return true; // 160329-165887
    if( passUnprescaledHLTTriggerPattern("HLT_IsoMu24_v") )          return true; // 160329-173198
    if( passUnprescaledHLTTriggerPattern("HLT_IsoMu30_eta2p1_v") )   return true; // 173212-180291
  }

  else{
    cout << "susySelections.cc:: " << __LINE__ << " ERROR unrecognized lepType " << lepType << ", quitting" << endl;
    exit(0);
  }

  return false;
}

bool passSingleLep3JetSUSYTrigger2011( bool isData , int lepType ) {

  //-------------------------------------------------------
  //These are the triggers for lepton+3jets
  //-------------------------------------------------------

  // no triggers required for MC
  if( !isData ) return true;

  // electron channel
  if( lepType == 0 ){
    if( passUnprescaledHLTTriggerPattern("HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v") )                     return true; // 160329-164236
    if( passUnprescaledHLTTriggerPattern("HLT_Ele25_CaloIdVT_TrkIdT_TriCentralJet30_v") )                     return true; // 165088-165887
    if( passUnprescaledHLTTriggerPattern("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30_v") )    return true; // 165922-178380
    if( passUnprescaledHLTTriggerPattern("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v") )  return true; // 178420-180291
  }

  // muon channel
  else if( lepType == 1 ){
    if( passUnprescaledHLTTriggerPattern("HLT_Mu17_TriCentralJet30_v") )                   return true;   // 160329-165887
    if( passUnprescaledHLTTriggerPattern("HLT_IsoMu17_TriCentralJet30_v") )                return true;   // 165922-173198
    if( passUnprescaledHLTTriggerPattern("HLT_IsoMu17_eta2p1_TriCentralJet30_v") )         return true;   // 173212-178380
    if( passUnprescaledHLTTriggerPattern("HLT_IsoMu17_eta2p1_TriCentralPFJet30_v") )       return true;   // 178420-180291
  }

  else{
    cout << "susySelections.cc:: ERROR unrecognized lepType " << lepType << ", quitting" << endl;
    exit(0);
  }

  return false;
}


/*****************************************************************************************/
//passes the OS SUSY trigger selection 2011
/*****************************************************************************************/

bool passSUSYTrigger2011_v1( bool isData , int hypType , bool highpt ) {
  
  //----------------------------------------
  // no trigger requirements applied to MC
  //----------------------------------------
  
  if( !isData ) return true; 
  
  //---------------------------------
  // triggers for lepton-HT datasets
  //---------------------------------
  
  if( !highpt ) {
  
    //mm
    if( hypType == 0 ){
      if( passUnprescaledHLTTriggerPattern("HLT_DoubleMu3_HT150_v") )       return true;
      if( passUnprescaledHLTTriggerPattern("HLT_DoubleMu3_HT160_v") )       return true;
      if( passUnprescaledHLTTriggerPattern("HLT_DoubleMu3_Mass4_HT150_v") ) return true;
      if( passUnprescaledHLTTriggerPattern("HLT_DoubleMu5_Mass4_HT150_v") ) return true;
      if( passUnprescaledHLTTriggerPattern("HLT_DoubleMu5_Mass8_HT150_v") ) return true;
    }
    
    //em
    else if( hypType == 1 || hypType == 2 ){
      if( passUnprescaledHLTTriggerPattern("HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT150_v") )       return true; 
      if( passUnprescaledHLTTriggerPattern("HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT150_v") )       return true;
      if( passUnprescaledHLTTriggerPattern("HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT160_v") )       return true; 
      if( passUnprescaledHLTTriggerPattern("HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT160_v") )       return true;
      if( passUnprescaledHLTTriggerPattern("HLT_Mu5_Ele8_CaloIdT_TrkIdVL_Mass4_HT150_v") ) return true;
      if( passUnprescaledHLTTriggerPattern("HLT_Mu5_Ele8_CaloIdT_TrkIdVL_Mass8_HT150_v") ) return true;
    }
    
    //ee
    else if( hypType == 3 ){
      if( passUnprescaledHLTTriggerPattern("HLT_DoubleEle8_CaloIdL_TrkIdVL_HT150_v") )       return true;
      if( passUnprescaledHLTTriggerPattern("HLT_DoubleEle8_CaloIdT_TrkIdVL_HT150_v") )       return true;
      if( passUnprescaledHLTTriggerPattern("HLT_DoubleEle8_CaloIdL_TrkIdVL_HT160_v") )       return true;
      if( passUnprescaledHLTTriggerPattern("HLT_DoubleEle8_CaloIdT_TrkIdVL_HT160_v") )       return true;
      if( passUnprescaledHLTTriggerPattern("HLT_DoubleEle8_CaloIdT_TrkIdVL_Mass8_HT150_v") ) return true;
    }
  }
  
  //---------------------------------
  // triggers for dilepton datasets
  //---------------------------------
  
  else{
  
    //mm
    if( hypType == 0 ){
      if( passUnprescaledHLTTriggerPattern("HLT_DoubleMu7_v") )   return true;
      if( passUnprescaledHLTTriggerPattern("HLT_Mu13_Mu7_v" ) )   return true;
      if( passUnprescaledHLTTriggerPattern("HLT_Mu13_Mu8_v" ) )   return true;
      if( passUnprescaledHLTTriggerPattern("HLT_Mu17_Mu8_v" ) )   return true;
    }
    
    //em
    else if( hypType == 1 || hypType == 2 ){
      if( passUnprescaledHLTTriggerPattern("HLT_Mu17_Ele8_CaloIdL_v") )   return true;
      if( passUnprescaledHLTTriggerPattern("HLT_Mu8_Ele17_CaloIdL_v") )   return true;
      if( passUnprescaledHLTTriggerPattern("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v") )   return true;
      if( passUnprescaledHLTTriggerPattern("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v") )   return true;

    }
    
    //ee
    else if( hypType == 3 ){
      if( passUnprescaledHLTTriggerPattern("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v") )                                   return true;
      if( passUnprescaledHLTTriggerPattern("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v") ) return true;
      if( passUnprescaledHLTTriggerPattern("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") ) return true;
    }                                     
  }        
  
  return false;
    
}


bool passSUSYTrigger2012_v1( int hypType ) {

  if (!cms2.evt_isRealData()) return true;

  //mm
  if( hypType == 0 ){
    if( passUnprescaledHLTTriggerPattern("HLT_Mu17_Mu8_v" ) )     return true;
    if( passUnprescaledHLTTriggerPattern("HLT_Mu17_TkMu8_v" ) )   return true;
  }
  
  //em
  else if( hypType == 1 || hypType == 2 ){
    if( passUnprescaledHLTTriggerPattern("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") ) return true;
    if( passUnprescaledHLTTriggerPattern("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") ) return true;
 }
  
  //ee
  else if( hypType == 3 ){
    if( passUnprescaledHLTTriggerPattern("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") ) return true;
  }

  return false;
}

bool passSUSYTrigger2012_v2( bool isData ) {

  if( !isData ) return true;

  if( passUnprescaledHLTTriggerPattern("HLT_Mu17_Mu8_v" ) )                                           return true;
  if( passUnprescaledHLTTriggerPattern("HLT_Mu17_TkMu8_v" ) )                                         return true;
  if( passUnprescaledHLTTriggerPattern("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") )        return true;
  if( passUnprescaledHLTTriggerPattern("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") )        return true;
  if( passUnprescaledHLTTriggerPattern("HLT_IsoMu24_v" ) )                                            return true;
  if( passUnprescaledHLTTriggerPattern("HLT_IsoMu24_eta2p1_v" ) )                                     return true;
  if( passUnprescaledHLTTriggerPattern("HLT_Ele27_WP80_v" ) )                                         return true;
  if( passUnprescaledHLTTriggerPattern("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") ) return true;

  return false;
}


/*****************************************************************************************/
//passes the OS SUSY trigger selection
/*****************************************************************************************/
bool passSUSYTrigger_v1( bool isData , int hypType ) {

  //int run = cms2.evt_run();

  //currently do not require any triggers for MC
  if( !isData ) return true;

  //mumu
  if( hypType == 0 ){

    // This is overkill, as Mu15 should be a susbset of Mu11 (for example)
    // But why not.
  
    if( passUnprescaledHLTTrigger("HLT_DoubleMu3_v2") )   return true;
    if( passUnprescaledHLTTrigger("HLT_DoubleMu5_v1") )   return true;
    if( passUnprescaledHLTTrigger("HLT_Mu11") )           return true;
    if( passUnprescaledHLTTrigger("HLT_Mu13_v1") )        return true;
    if( passUnprescaledHLTTrigger("HLT_Mu15_v1") )        return true;
    if( passUnprescaledHLTTrigger("HLT_Mu17_v1") )        return true; //<<<---Added 2e32
    if( passUnprescaledHLTTrigger("HLT_Mu19_v1") )        return true; //<<<---Added 2e32
    
    //if( run <= 147116 ){
    if( passUnprescaledHLTTrigger("HLT_DoubleMu3") )    return true; //136033-147116
    if( passUnprescaledHLTTrigger("HLT_Mu9") )          return true; //136033-147116
    //}
    //if( run <= 144114 ){
    if( passUnprescaledHLTTrigger("HLT_Mu7") )          return true; //140116-144114
    //}
    //if( run <= 141882 ){
    if( passUnprescaledHLTTrigger("HLT_Mu5") )          return true; //136033-141882
    //}
    
  }
 


  //ee
  else if( hypType == 3 ){

    //Added for 2e32 single ele trigs:
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v1")) return true;
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2")) return true;
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3")) return true;

    if( passUnprescaledHLTTrigger("HLT_Ele22_SW_TighterEleId_L1R_v2"))     return true;
    if( passUnprescaledHLTTrigger("HLT_Ele22_SW_TighterEleId_L1R_v3")) return true;
    if( passUnprescaledHLTTrigger("HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2")) return true;

    if( passUnprescaledHLTTrigger("HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1")) return true;

    if( passUnprescaledHLTTrigger("HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1")) return true;
    if( passUnprescaledHLTTrigger("HLT_Ele32_SW_TighterEleId_L1R_v2")) return true;

   // This is a family of never prescaled Ele17 triggers
    // Some of them are probably prescaled at 2e32?
    //got prescaled at 148058...
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TightEleId_L1R") )          return true; 
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TighterEleId_L1R_v1") )     return true;
    
    // These are unprescaled double triggers
    if( passUnprescaledHLTTrigger("HLT_DoubleEle15_SW_L1R_v1") )                 return true;
    if( passUnprescaledHLTTrigger("HLT_DoubleEle17_SW_L1R_v1") )                 return true; //    
           
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1") ) return true; // 147390-->
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v2") ) return true; // 147390-->
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1") )  return true; // 147196-148058 

    // These are double triggers that became prescaled at some point
    //if( run <= 147116 ){
    if( passUnprescaledHLTTrigger("HLT_DoubleEle10_SW_L1R") ) return true; //141956-147116
    //}
    //if( run <= 141882 ){
    if( passUnprescaledHLTTrigger("HLT_DoubleEle5_SW_L1R") )  return true; //136033-141882
    //}
      
      
    //if( run <= 147116 ){
    // This is a family of Ele17 triggers that are not
    // there anymore or are not unprescaled anymore.
    // The loosest one is indicated as master
    // In principle, the other two could be skipped
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_CaloEleId_L1R") )                 return true; //146428-147116 <---- master
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_EleId_L1R") )                     return true; //146428-147116
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_LooseEleId_L1R") )                return true; //146428-147116
    //}

    //if( run <= 144114 ){
    // This is a family of Ele15 triggers that are not
    // there anymore or are not unprescaled anymore.
    // The loosest one is indicated as master
    // In principle, the other one could be skipped
    if( passUnprescaledHLTTrigger("HLT_Ele15_SW_CaloEleId_L1R") ) return true; //141956-144114 <---- master
    if( passUnprescaledHLTTrigger("HLT_Ele15_SW_EleId_L1R") )     return true; //141956-144114
    //}
    
    // These are Ele15 that got prescaled in a different
    // time frame as the previous family...use both
    //if( run <= 143962 ){
    if( passUnprescaledHLTTrigger("HLT_Ele15_SW_L1R") ) return true; //140058-143962
    //}
    //if( run <= 141882 ){
    if( passUnprescaledHLTTrigger("HLT_Ele15_LW_L1R") ) return true; //136033-141882
    //}
    
    // This is a bonus trigger, should not be needed, but why not?
    //if( run <= 144114 ){
    if( passUnprescaledHLTTrigger("HLT_Ele20_SW_L1R") ) return true; //140058-144114
    //}

    // This is the Ele10 family that came and went at ramdom times
    // Use all (although in principle the one with EleID could be
    // skipped in the subset of runs where triggers with no EleId
    // were present)
    //if( run <= 144114 ){
    if( passUnprescaledHLTTrigger("HLT_Ele10_SW_EleId_L1R") ) return true;      //141956-144114
    //}
    //if( run <= 141882 ){
    if( passUnprescaledHLTTrigger("HLT_Ele10_LW_EleId_L1R") ) return true;      //136033-141882
    //} 
    //if( run <= 139980 ){
    if( passUnprescaledHLTTrigger("HLT_Ele10_LW_L1R") )       return true;      //136033-139980
    if( passUnprescaledHLTTrigger("HLT_Ele10_SW_L1R") )       return true;      //139195-139980
    //}

  }

  //emu
  else if( hypType == 1 || hypType == 2 ){
    
    //---------------------------
    // single muon triggers
    //---------------------------
    if( passUnprescaledHLTTrigger("HLT_Mu11") )           return true;
    if( passUnprescaledHLTTrigger("HLT_Mu13_v1") )        return true;
    if( passUnprescaledHLTTrigger("HLT_Mu15_v1") )        return true;
    
    //if( run <= 147116 ){
    if( passUnprescaledHLTTrigger("HLT_Mu9") )          return true; //136033-147116
    //}
    //if( run <= 144114 ){
    if( passUnprescaledHLTTrigger("HLT_Mu7") )          return true; //140116-144114
    //}
    //if( run <= 141882 ){
    if( passUnprescaledHLTTrigger("HLT_Mu5") )          return true; //136033-141882
    //}

    if( passUnprescaledHLTTrigger("HLT_Mu17_v1") )        return true; //<<<---Added 2e32
    if( passUnprescaledHLTTrigger("HLT_Mu19_v1") )        return true; //<<<---Added 2e32
 
    
    //---------------------------
    // single electron triggers
    //---------------------------
    
    //Added 2e32
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v1")) return true;
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2")) return true;
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3")) return true;

    if( passUnprescaledHLTTrigger("HLT_Ele22_SW_TighterEleId_L1R_v2"))     return true;
    if( passUnprescaledHLTTrigger("HLT_Ele22_SW_TighterEleId_L1R_v3")) return true;
    if( passUnprescaledHLTTrigger("HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2")) return true;

    if( passUnprescaledHLTTrigger("HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1")) return true;

    if( passUnprescaledHLTTrigger("HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1")) return true;
    if( passUnprescaledHLTTrigger("HLT_Ele32_SW_TighterEleId_L1R_v2")) return true;

    
    // This is a family of never prescaled Ele17 triggers
    // Some of them are probably prescaled at 2e32?
    // they did get prescaled 148058... never say never...
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TightEleId_L1R") )          return true;
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TighterEleId_L1R_v1") )     return true;
    
    //if( run <= 147116 ){
    // This is a family of Ele17 triggers that are not
    // there anymore or are not unprescaled anymore.
    // The loosest one is indicated as master
    // In principle, the other two could be skipped
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_CaloEleId_L1R") )                 return true; //146428-147116
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_EleId_L1R") )                     return true; //146428-147116
    if( passUnprescaledHLTTrigger("HLT_Ele17_SW_LooseEleId_L1R") )                return true; //146428-147116
    //}
    
    //if( run <= 144114 ){
    // This is a family of Ele15 triggers that are not
    // there anymore or are not unprescaled anymore.
    // The loosest one is indicated as master
    // In principle, the other one could be skipped
    if( passUnprescaledHLTTrigger("HLT_Ele15_SW_CaloEleId_L1R") ) return true; //141956-144114 
    if( passUnprescaledHLTTrigger("HLT_Ele15_SW_EleId_L1R") )     return true; //141956-144114
    //}
    
    // These are Ele15 that got prescaled in a different
    // time frame as the previous family...use both
    //if( run <= 143962 ){
    if( passUnprescaledHLTTrigger("HLT_Ele15_SW_L1R") ) return true; //140058-143962
    //}
    //if( run <= 141882 ){
    if( passUnprescaledHLTTrigger("HLT_Ele15_LW_L1R") ) return true; //136033-141882
    //}
    
    // This is a bonus trigger, should not be needed, but why not?
    //if( run <= 144114 ){
    if( passUnprescaledHLTTrigger("HLT_Ele20_SW_L1R") ) return true; //140058-144114
    //}
    
    // This is the Ele10 family that came and went at ramdom times
    // Use all (although in principle the one with EleID could be
    // skipped in the subset of runs where triggers with no EleId
    // were present)
    //if( run <= 144114 ){
    if( passUnprescaledHLTTrigger("HLT_Ele10_SW_EleId_L1R") ) return true;      //141956-144114
    //}
    //if( run <= 141882 ){
    if( passUnprescaledHLTTrigger("HLT_Ele10_LW_EleId_L1R") ) return true;      //136033-141882
    //} 
    //if( run <= 139980 ){
    if( passUnprescaledHLTTrigger("HLT_Ele10_LW_L1R") )       return true;      //136033-139980
    if( passUnprescaledHLTTrigger("HLT_Ele10_SW_L1R") )       return true;      //139195-139980
    //}
    
    //---------------------------
    // e-mu cross triggers
    //---------------------------
    
    if( passUnprescaledHLTTrigger("HLT_Mu5_Ele5_v1") )          return true;
    if( passUnprescaledHLTTrigger("HLT_Mu5_Ele9_v1") )          return true;
    
    //Added 2e32
    if( passUnprescaledHLTTrigger("HLT_Mu11_Ele8_v1") )     return true;    // Mu3 L1 Seed
    if( passUnprescaledHLTTrigger("HLT_Mu8_Ele8_v1") )      return true;    // Mu3 L1 Seed
    if( passUnprescaledHLTTrigger("HLT_Mu5_Ele13_v2") )     return true;    // EG8 L1 Seed
    if( passUnprescaledHLTTrigger("HLT_Mu5_Ele17_v1") )     return true;    // EG8 L1 Seed & open mu Seed
    if( passUnprescaledHLTTrigger("HLT_Mu5_Ele17_v2") )     return true;    // EG8 L1 Seed & open mu Seed

  }
  
  else{
    cout << "ERROR: unrecognized hyptype " << hypType << ", quitting" << endl;
    exit(0);
  }

  return false;
}


//----------------------------------------------------------
//this is a simplified version of the SUSY triggers
//----------------------------------------------------------


bool passSimpleSUSYTrigger_v1( bool isData ) {

  //currently do not require any triggers for MC
  if( !isData ) return true;

  if( passUnprescaledHLTTrigger("HLT_Mu9") )                                    return true;            
  if( passUnprescaledHLTTrigger("HLT_Mu15_v1 147196-148862") )                  return true;
  if( passUnprescaledHLTTrigger("HLT_Ele10_LW_L1R") )                           return true;             
  if( passUnprescaledHLTTrigger("HLT_Ele10_LW_EleId_L1R") )                     return true;      
  if( passUnprescaledHLTTrigger("HLT_Ele15_SW_L1R") )                           return true; 
  if( passUnprescaledHLTTrigger("HLT_Ele15_SW_CaloEleId_L1R") )                 return true;  
  if( passUnprescaledHLTTrigger("HLT_Ele17_SW_CaloEleId_L1R") )                 return true; 
  if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TightEleId_L1R") )                return true; 
  if( passUnprescaledHLTTrigger("HLT_DoubleEle10_SW_L1R") )                     return true;     
  if( passUnprescaledHLTTrigger("HLT_DoubleEle15_SW_L1R_v1") )                  return true; 
  if( passUnprescaledHLTTrigger("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1") )  return true; 
  if( passUnprescaledHLTTrigger("HLT_Mu5_Ele5_v1") )                            return true; 
  if( passUnprescaledHLTTrigger("HLT_Mu8_Ele8_v1") )                            return true; 
  if( passUnprescaledHLTTrigger("HLT_Mu5_Ele17_v1") )                           return true;

  return false;
}
/*****************************************************************************************/
//hypothesis disambiguation. Returns the hypothesis that has mass closest to MZ
/*****************************************************************************************/
unsigned int selectBestZHyp(const vector<unsigned int> &v_goodHyps) {
  
  float mindeltam         = 100000;
  unsigned int bestHypIdx = 0;
  for(unsigned int i = 0; i < v_goodHyps.size(); i++) {
    unsigned int index = v_goodHyps.at(i);
    float deltam = fabs( hyp_p4()[index].mass() - 91. );
    if( deltam < mindeltam ) {
      mindeltam = deltam;
      bestHypIdx = index;
    }
  }

  return bestHypIdx;

}

/*
int sfinalState(int ipart1, int ipart2) {
//#                   ng     neutralino/chargino + gluino                     #
//#                   ns     neutralino/chargino + squark                     #
//#                   nn     neutralino/chargino pair combinations            #
//#                   ll     slepton pair combinations                        #
//#                   sb     squark-antisquark                                #
//#                   ss     squark-squark                                    #
//#                   tb     stop-antistop                                    #
//#                   bb     sbottom-antisbottom                              #
//#                   gg     gluino pair                                      #
//#                   sg     squark + gluino                                  #

 int index = -1;

// Ugly can be done in a better way using list/objects ... oh well
 // ng case
 if ((abs(ipart1) > 1000021 && abs(ipart1) < 1000039 && abs(ipart2) == 1000021) ||
      (abs(ipart2) > 1000021 && abs(ipart2) < 1000039 && abs(ipart1) == 1000021)) index = 0;
// ns case 
 else if ((abs(ipart1) > 1000021 && abs(ipart1) < 1000039 && (abs(ipart2) < 1000007 || (abs(ipart2) > 2000000 && abs(ipart2) < 2000007))) ||
    (abs(ipart2) > 1000021 && abs(ipart2) < 1000039 && (abs(ipart1) < 1000007 || (abs(ipart1) > 2000000 && abs(ipart1) < 2000007)))) index = 1;
// nn case 
 else if (abs(ipart1) > 1000021 && abs(ipart1) < 1000039 &&  abs(ipart2) > 1000021 && abs(ipart2) < 1000039) index = 2;
// ll case
 else if (((abs(ipart1) > 1000010 && abs(ipart1) < 1000017) || (abs(ipart1) > 2000010 && abs(ipart1) < 2000016)) &&
         ((abs(ipart2) > 1000010 && abs(ipart2) < 1000017) || (abs(ipart2) > 2000010 && abs(ipart2) < 2000016))) index = 3;
// tb case
 else if ((abs(ipart1) == 1000006 || abs(ipart1) == 2000006) && (abs(ipart2) == 1000006 || abs(ipart2) == 2000006)) index = 6;
// bb case
 else if ((abs(ipart1) == 1000005 || abs(ipart1) == 2000005) && (abs(ipart2) == 1000005 || abs(ipart2) == 2000005)) index = 7;
// sb case
 else if ((ipart1 * ipart2 < 0) &&
         (abs(ipart1) < 1000007 || (abs(ipart1) > 2000000 && abs(ipart1) < 2000007)) &&
         (abs(ipart2) < 1000007 || (abs(ipart2) > 2000000 && abs(ipart2) < 2000007))) index = 4;
// ss case
 else if ((ipart1 * ipart2 > 0) &&
         (abs(ipart1) < 1000007 || (abs(ipart1) > 2000000 && abs(ipart1) < 2000007)) &&
         (abs(ipart2) < 1000007 || (abs(ipart2) > 2000000 && abs(ipart2) < 2000007))) index = 5;
// gg case 
  else if (ipart1 == 1000021 && ipart2 == 1000021) index = 8;
// sg case 
 else if ((abs(ipart1) == 1000021 && (abs(ipart2) < 1000007 || (abs(ipart2) > 2000000 && abs(ipart2) < 2000007))) ||
    (abs(ipart2) == 1000021 && (abs(ipart1) < 1000007 || (abs(ipart1) > 2000000 && abs(ipart1) < 2000007)))) index = 9;
 else
   cout << "Warning mcSUSYkfactor: index not defined " << ipart1 << "  " << ipart2 << endl;

  return index;
}
*/

/*
double lmdata(int ipart1, int ipart2, string prefix) {

  static const Double_t lm0[10] = {1.06604, 1.00369, 1.27186, 1.19103, 1.44681, 1.22883, 1.5649, 1.70195, 1.99721, 1.33951};
  static const Double_t lm0scale[10] = {1.04522, 1.0001, 1.27058, 1.19107, 1.39667, 1.18765, 1.53141, 1.66519, 1.88761, 1.27273 };

  static const Double_t lm1[10] = {0.988762, 1.02694, 1.22581, 1.24139, 1.45174, 1.18939, 1.72769, 1.7931, 2.33333, 1.39205};
  static const Double_t lm2[10] = {0.971508, 1.07009, 1.18325, 1.1879, 1.47269, 1.16915, 1.86207, 1.94417, 2.96117, 1.57895};
  static const Double_t lm3[10] = {1.03489, 1.03535, 1.2383, 1.12839, 1.42798, 1.20694, 1.74914, 1.83868, 2.46354, 1.46457};
  static const Double_t lm4[10] = {0.978699, 1.04381, 1.21136, 1.16334, 1.47239, 1.18924, 1.779, 1.87219, 2.54733, 1.47039};
  static const Double_t lm5[10] = {0.969489, 1.07524, 1.17991, 1.14481, 1.48387, 1.16456, 1.88094, 2.00574, 3.05195, 1.6087};
  static const Double_t lm6[10] = {0.971497, 1.08695, 1.16543, 1.1852, 1.47826, 1.15094, 1.9199, 2.04659, 3.32237, 1.66159};
  static const Double_t lm8[10] = {1.02828, 1.07339, 1.2137, 1.08113, 1.43871, 1.20667, 1.82724, 1.97961, 2.91919, 1.63415};
  static const Double_t lm9[10] = {1.08102, 1.21649, 3.31807, 1.09983, 2.16036, 1.30526, 2.18081, 2.2218, 2.17647, 2.14039};

 
  int index = sfinalState(ipart1, ipart2);
  if (index < 0) return 1.0;

  if (prefix == "lm0") return lm0[index];
  else if (prefix == "lm1") return lm1[index];
  else if (prefix == "lm2") return lm2[index];
  else if (prefix == "lm3") return lm3[index];
  else if (prefix == "lm4") return lm4[index];
  else if (prefix == "lm5") return lm5[index];
  else if (prefix == "lm6") return lm6[index];
  else if (prefix == "lm8") return lm8[index];
  else if (prefix == "lm9") return lm9[index];
  else if (prefix == "lm0scale") return lm0scale[index];

  return 1.0;
}




float kfactorSUSY(string sample)
{
  float kfactor = 1.0;
  TDatabasePDG *pdg = new TDatabasePDG();
  std::vector<int> interactions;
  
  for (int j=0; j<cms2.genps_id().size(); j++) {
    if (cms2.genps_status().at(j) != 3) continue;
    int ID = abs(cms2.genps_id().at(j));
    int mID = abs(cms2.genps_id_mother().at(j));
    if (ID > 1000000 && ID < 2000016 && (mID < 7 || mID ==21 || mID == 22 || mID == 23 || mID == 24)) { // kept the bosons in case of screw ups

     // Check the mother of the SM is a proton
       bool isProton = false;       
       for (int k=0; k < j; k++) {
            int pID = abs(cms2.genps_id().at(k));
            int mpID = abs(cms2.genps_id_mother().at(k));
         if ((pID < 7 || pID ==21 || pID == 22 || pID == 23 || pID == 24)&&(mpID == 2212)) isProton = true; 
       }
      if (!isProton) continue;  
      interactions.push_back(cms2.genps_id().at(j));

    if (interactions.size() > 2) 
         cout << setw(4) << left << j << " WARNING mcSUSYkfactor: Something is wrong with "
         << setw(10) << left << pdg->GetParticle(cms2.genps_id().at(j))->GetName() << " "
         << setw(10) << left << cms2.genps_id().at(j) << " "
         << setw(7) << right << setprecision(4) << cms2.genps_p4().at(j).pt() << "  "
         << setw(7) << right << setprecision(4) << cms2.genps_p4().at(j).phi() << "  "
         << setw(10) << right << setprecision(4) << cms2.genps_p4().at(j).eta() << "  "
         << setw(4) << right << cms2.genps_status().at(j) << " "
         << setw(10) << left << pdg->GetParticle(cms2.genps_id_mother().at(j))->GetName()
         << " using k=1 " << endl;
   }
  }

  delete pdg;

  if (interactions.size() == 2 ) {
     kfactor = lmdata(interactions[0], interactions[1], sample);
   } else {
    cout << "WARNING mcSUSYkfactor: Number of sparticles found: " << interactions.size() << " using kfactor=1" <<endl;
    kfactor = 1.0;
//    dumpDocLines();
  }

  return kfactor;
}
*/
