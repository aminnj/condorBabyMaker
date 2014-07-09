#ifndef SUSYSELECTIONS_H
#define SUSYSELECTIONS_H

#include <vector>
#include "Math/LorentzVector.h"
#include "muonSelections.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

int leptonOrTauIsFromW(int idx, int id, bool alsoSusy = false);

//---------------------------------------------
// single muon triggers for lljj bump search
//---------------------------------------------

bool passMuMuJJTrigger_v1( bool isData );

/*****************************************************************************************/
//print event info
/*****************************************************************************************/
void printEventInfo();

/*****************************************************************************************/
//veto Z->mumugamma events
/*****************************************************************************************/
bool vetoZmumuGamma( unsigned int hypIdx , float emax = 6. , 
                     float minmass = 76. , float maxmass = 106.);

/*****************************************************************************************/
//passes the SUSY 2011 trigger selections
/*****************************************************************************************/
bool passSUSYTrigger2011_v1( bool isData , int hypType , bool highpt );

/*****************************************************************************************/
//passes the SUSY 2012 trigger selections
/*****************************************************************************************/
bool passSUSYTrigger2012_v1( int hypType );
bool passSUSYTrigger2012_v2( bool isData );

/*****************************************************************************************/
//passes the single-lepton SUSY 2011 trigger selections
/*****************************************************************************************/
bool passSingleLepSUSYTrigger2011_v1( bool isData , int lepType );
bool passSingleLep2JetSUSYTrigger2011( bool isData , int lepType );
bool passSingleLep3JetSUSYTrigger2011( bool isData , int lepType );
bool passSingleMuTrigger2011( bool isData , int lepType );

/*****************************************************************************************/
//passes the SUSY trigger selections
/*****************************************************************************************/
bool passSUSYTrigger_v1( bool isData , int hypType );

/*****************************************************************************************/
//passes the simplified version of the SUSY trigger selections
/*****************************************************************************************/
bool passSimpleSUSYTrigger_v1( bool isData );


/*****************************************************************************************/
//hypothesis disambiguation. Returns the hypothesis that has mass closest to MZ
/*****************************************************************************************/
unsigned int selectBestZHyp(const vector<unsigned int> &v_goodHyps);

/*****************************************************************************************/
//generalized Z veto
/*****************************************************************************************/
bool ZVetoGeneral( float ptcut = 20 , float minmass = 76 ,  float maxmass = 106 , SelectionType = OSGeneric_v3 );

//---------------------------------------------
// Check if trigger is unprescaled and passes
//---------------------------------------------
bool passUnprescaledHLTTriggerPattern(const char* arg);
bool passHLTTriggerPattern(const char* arg);
int passTriggerPrescale(const char* arg);
TString triggerName(TString triggerPattern);

bool passElectronSelection_ZMet2012_v1_NoIso(int index , bool vetoTransition=false, bool eta24=false);
bool passElectronSelection_ZMet2012_v1_Iso(int index   , bool vetoTransition=false, bool eta24=false);
bool passElectronSelection_ZMet2012_v1(int index       , bool vetoTransition=false, bool eta24=false);

bool passElectronSelection_ZMet2012_v2_NoIso(int index , bool vetoTransition=false, bool eta24=false, bool useOldIsolation=false);
bool passElectronSelection_ZMet2012_v2_Iso(int index   , bool vetoTransition=false, bool eta24=false, bool useOldIsolation=false);
bool passElectronSelection_ZMet2012_v2(int index       , bool vetoTransition=false, bool eta24=false, bool useOldIsolation=false);

bool passElectronSelection_ZMet2012_v3_NoIso(int index , bool vetoTransition=false, bool eta24=false, bool useOldIsolation=false);
bool passElectronSelection_ZMet2012_v3_Iso(int index   , bool vetoTransition=false, bool eta24=false, bool useOldIsolation=false);
bool passElectronSelection_ZMet2012_v3(int index       , bool vetoTransition=false, bool eta24=false, bool useOldIsolation=false);

bool passElectronSelection_Stop2012_v1_NoIso(int index , bool vetoTransition=false, bool eta24=false);
bool passElectronSelection_Stop2012_v1_Iso(int index   , bool vetoTransition=false, bool eta24=false);
bool passElectronSelection_Stop2012_v1(int index       , bool vetoTransition=false, bool eta24=false);

bool passElectronSelection_Stop2012_v2_NoIso(int index , bool vetoTransition=false, bool eta24=false, bool useOldIsolation=false);
bool passElectronSelection_Stop2012_v2_Iso(int index   , bool vetoTransition=false, bool eta24=false, bool useOldIsolation=false);
bool passElectronSelection_Stop2012_v2(int index       , bool vetoTransition=false, bool eta24=false, bool useOldIsolation=false);

bool passElectronSelection_Stop2012_v3_NoIso(int index , bool vetoTransition=false, bool eta24=false, bool useOldIsolation=false);
bool passElectronSelection_Stop2012_v3_Iso(int index   , bool vetoTransition=false, bool eta24=false, bool useOldIsolation=false);
bool passElectronSelection_Stop2012_v3(int index       , bool vetoTransition=false, bool eta24=false, bool useOldIsolation=false);


#endif

