
# ifndef MCSELECTIONS_H
# define MCSELECTIONS_H

#include "CMS2.h"

int getDrellYanType();
bool isDYee();
bool isDYmm();
bool isDYtt();
bool isWe();
bool isWm();
bool isWt();
bool isWW();
bool isWZ();
bool isZZ();
int getZZType();
int dumpDocLines();

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
//          v     = 4-vector of reco lepton
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
int leptonIsFromW(int idx, int id, bool alsoSusy=false);
int ttbarconstituents(int i_hyp);

bool trueGammaFromMuon(int electron);
// count genp leptons
//-------------------------------------------------- 
// Returns the number of e,mu, and tau in the doc lines 
//----------------------------------------------------- 
int leptonGenpCount(int& nele, int& nmuon, int& ntau);
int leptonGenpCount_lepTauDecays(int& nele, int& nmuon, int& ntau);
int genpDileptonType();
  //0 mumu; 1 emu; 2 ee
  
// -----------------------------------------
// MC helper functions for fakerate tests:
// -----------------------------------------
int elFakeMCCategory(int i_el);
int muFakeMCCategory(int i_mu);
bool idIsCharm(int id);
bool idIsBeauty(int id);

// -------------------------------------
// quick and dirty function to see if
// lepton is prompt (approximately)
// -------------------------------------
bool isNotPromptSusyLeptonHyp(int idx);

int mc3idx_eormu(int id, int idx, float maxDr = 0.5, float minPt = 1);
float mc3dr_eormu(int id, int idx, float maxDr = 0.5, float minPt = 1);

#endif

