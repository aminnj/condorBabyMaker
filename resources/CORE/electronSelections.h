#ifndef ELECTRONSELECTIONS_H
#define ELECTRONSELECTIONS_H

// C++
#include <stdint.h>
#include <vector>

// ROOT
#include "TMath.h"

// Header
#include "electronSelectionsParameters.h"

// CMS2
#include "CMS2.h"

//
typedef ULong64_t   uint64;
typedef uint64      cuts_t;
typedef uint64      electronIdComponent_t;

/////////////////////////////////////////////
// This is the menu of electron selections //
/////////////////////////////////////////////

enum EleSelectionType {

///////////////
// Isolation //
///////////////
 
  ELEISO_REL010,             // rel iso (fixed below 20 GeV) < 0.10, 0.3 cone size for all, 1 GeV pedestal sub in EB
  ELEISO_REL015,             // rel iso (fixed below 20 GeV) < 0.15, 0.3 cone size for all, 1 GeV pedestal sub in EB
  ELEISO_REL040,             // rel iso (fixed below 20 GeV) < 0.40, 0.3 cone size for all, 1 GeV pedestal sub in EB
  ELEISO_REL100,             // rel iso (fixed below 20 GeV) < 1.00, 0.3 cone size for all, 1 GeV pedestal sub in EB
  ELEISO_REL010_WW,          // rel iso (fixed below 20 GeV) < 0.10, 0.3 cone size for all, 1 GeV pedestal sub in EB/EE
  ELEISO_REL040_WW,          // rel iso (fixed below 20 GeV) < 0.40, 0.3 cone size for all, 1 GeV pedestal sub in EB/EE
  ELEISO_REL100_WW,          // rel iso (fixed below 20 GeV) < 1.00, 0.3 cone size for all, 1 GeV pedestal sub in EB/EE
  ELEISO_SMURFV4,            // non-truncated relative pf iso with cut [0.15,0.09] for [barrel,endcap]
  ELEISO_SMURFV5,            // non-truncated relative pf iso with cut [0.13,0.09] for [barrel,endcap]
  ELEISO_RELNT010,           // non-truncated relative iso < 0.10, 0.3 cone size for all, 1 GeV pedestal sub in EB
  ELEISO_RELNT015,           // non-truncated relative iso < 0.15, 0.3 cone size for all, 1 GeV pedestal sub in EB
  ELEISO_RELNT040,           // non-truncated relative iso < 0.40, 0.3 cone size for all, 1 GeV pedestal sub in EB
  ELEISO_TRK_RELNT020,       // non-truncated relative Tracker iso < 0.20, 0.3 cone size for all
  ELEISO_ECAL_RELNT020,      // non-truncated relative ECAL    iso < 0.20, 0.3 cone size for all, 1 GeV pedestal sub in EB
  ELEISO_ECAL_RELNT020_NPS,  // non-truncated relative ECAL    iso < 0.20, 0.3 cone size for all, no pedestal sub in EB
  ELEISO_HCAL_RELNT020,      // non-truncated relative HCAL    iso < 0.20, 0.3 cone size for all
  ELEISO_ECAL_REL020,        // truncated relative ECAL    iso < 0.20, 0.3 cone size for all, 1 GeV pedestal sub in EB
  ELEISO_HCAL_REL020,        // truncated relative HCAL    iso < 0.20, 0.3 cone size for all
  ELEISO_FASTJET_REL005,     // truncated reliso < 0.05, 0.3 cone size, 1 GeV pedestal subtraction, fastjet-corrected
  ELEISO_FASTJET_REL010,     // truncated reliso < 0.05, 0.3 cone size, 1 GeV pedestal subtraction, fastjet-corrected
  ELEISO_FASTJET_REL015,     // truncated reliso < 0.05, 0.3 cone size, 1 GeV pedestal subtraction, fastjet-corrected
  ELEISO_COR_RELNT010,       // truncated correction reliso < 0.1, cor_iso = ntiso - (ln(pt) * nvtx)/(30+pt)

//////////////////////
// Impact Parameter //
//////////////////////

  ELEIP_200,         // d0 corrected for beamspot < 0.02
  ELEIP_400,         // d0 corrected for beamspot < 0.04
  ELEIP_PV_200,      // d0 corrected for primary vertex < 0.02
  ELEIP_PV_wwV1,     // d0 (PV) < 0.02 and dz (PV) < 1.0
  ELEIP_PV_SMURFV3,  // d0 (PV) < 0.02 and dz (PV) < 0.2, using first PV
  ELEIP_PV_DZ_1MM,   // dz (PV) < 0.1, using first PV
  ELEIP_PV_OSV2,     // d0 (PV) < 0.04 and dz (PV) < 1, using first PV
  ELEIP_PV_OSV2_FO,  // d0 (PV) < 0.2 and dz (PV) < 1, using first PV
  ELEIP_SS200,       // 2011 SS numerator d0 cut


/////////////////////////////
// Electron Identification //
/////////////////////////////

  ELEID_SMURFV1_EXTRA,                // pass smurf v1 electron ID
  ELEID_SMURFV2_EXTRA,                // pass smurf v2 electron ID
  ELEID_SMURFV3_EXTRA,                // pass smurf v3 electron ID
  ELEID_SMURFV1SS_EXTRA,              // electron ID with VBTF80
  ELEID_SMURFV2SS_EXTRA,              // electron ID with VBTF80
  ELEID_VBTF_35X_95,                  // VBTF95 electron ID (35X)
  ELEID_VBTF_35X_90,                  // VBTF90 electron ID (35X)
  ELEID_VBTF_35X_80,                  // VBTF80 electron ID (35X)
  ELEID_VBTF_35X_70,                  // VBTF70 electron ID (35X)
  ELEID_VBTF_80_NOHOEEND,             // VBTF80 electron ID no HoE in endcap
  ELEID_VBTF_85_NOHOEEND,             // VBTF85 electron ID no HoE in endcap
  ELEID_VBTF_85,                      // VBTF85 electron ID
  ELEID_VBTF_70_NOHOEEND,             // VBTF70 electron ID no HoE in endcap
  ELEID_VBTF_90_HLT,                  // VBTF90 electron ID with HoE and dPhiIn cuts tuned to represent HLT requirements for CaloIdL_TrkIdVL
  ELEID_VBTF_90_HLT_CALOIDT_TRKIDVL,  // VBTF90 electron ID with HoE and dPhiIn cuts tuned to represent HLT requirements for CaloIdT_TrkIdVL
  ELEID_CIC_V03_MEDIUM,               // CIC_MEDIUM electron ID (V03)
  ELEID_VBTF_95_NOHOEEND,             // VBTF80 electron ID no HoE in endcap
  // ELEID_WP2012_LOOSE_NOISO,          // WP2012 LOOSE ELECTRON ID, NO ISO
  // ELEID_WP2012_LOOSE_NOISO_NOIP,     // WP2012 LOOSE ELECTRON ID, NO ISO, NO IP
  ELEID_WP2012_MEDIUM_NOISO,          // WP2012 MEDIUM ELECTRON ID, NO ISO
  ELEID_WP2012_MEDIUM_NOISO_NOIP,     // WP2012 MEDIUM ELECTRON ID, NO ISO, NO IP

//////////////////////////
// Conversion Rejection //
//////////////////////////

  ELENOTCONV_MIT,                // mit conversion rejection v11 
  ELENOTCONV_DISTDCOT002,        // dist < 0.02 && dcot(theta) < 0.02
  ELENOTCONV_HITPATTERN,         // < 2 missing hits
  ELENOTCONV_HITPATTERN_0MHITS,  // < 1 missing hits
  ELENOTCONV_DISTDCOT002_OLD,    // dist < 0.02 && dcot(theta) < 0.02 using dist/dcot taken from the GsfElectron object

//////////////////////
// Basic Selections //
//////////////////////

  ELEETA_250,               // |eta| < 2.50
  ELEETA_240,               // |eta| < 2.40 

  ELEPT_010,                // Pt > 10

  ELENOMUON_010,            // no muon within dR < 0.1
  ELENOMUON_010_SS,         // no muon passing same sign numerator selection within dR < 0.1

  ELESEED_ECAL,             // seed must have been found by at least the ecal algo

  ELECHARGE_NOTFLIP3AGREE,  // Not a charge flip and CTF, GSF, and Pixel-SuperCluster charges agree

  ELE_NOT_TRANSITION,       // SC |eta| < 1.4442 OR SC |eta| > 1.556 (veto transition region)
 
  //*************************************************************************************//
  // This is the last enumerator element                                                 //
  // The electron selection bitmask will not work with more than 63 element in this enum //
  // DO NOT ADD ANY ELEMENTS AFTER THIS!                                                 //
  //*************************************************************************************//

  ELE_LAST,

};

// Assuming the constants in EleSelectionType have default numeric values ( 0, 1, 2, ... N ),
// the last enum constant will have integer value N
// For the bitmasks to work, N must be <= 63
static bool shown = true;
inline void checkElectronSelections(void){
  int n    = (int) EleSelectionType(ELE_LAST);
  int nMax = (int) 8*sizeof(1ll) - 1;
  if( n > nMax ){
    cout << endl << "ERROR at line " << __LINE__ << " in " << __FILE__ << ":" << endl;
    cout << "enum \"EleSelectionType\" has " << n << " elements but cannot have more than " << nMax << " elements... Exiting." << endl << endl;
    exit(1);
  }
  else{
    if( !shown ){
      cout << endl << "There are " << ( nMax - n ) << " available selectors left in enum EleSelectionType" << endl;
      cout << "\t( " << __FILE__ << " )" << endl << endl;
    }
    shown = true;
  }
}

///////////////////
// Opposite Sign //
///////////////////

static const cuts_t electronSelection_el_OSV2_noiso = 
  (1ll<<ELEID_VBTF_90_HLT_CALOIDT_TRKIDVL) | // VBTF90, tightened to match CaloIdT+TrkIdVL
  (1ll<<ELEIP_PV_OSV2)                     | // d0(PV) < 0.04 cm, dz(PV) < 1.0 cm
  (1ll<<ELENOMUON_010)                     | // no muon dR < 0.1
  (1ll<<ELENOTCONV_HITPATTERN)             | // <=1 missing hits
  (1ll<<ELENOTCONV_MIT)                    | // MIT conversion rejection
  (1ll<<ELEPT_010)                         | // electron p_T > 10 GeV
  (1ll<<ELEETA_250);                         // |eta| < 2.5

static const cuts_t electronSelection_el_OSV2_iso = 
  (1ll<<ELEISO_ECAL_RELNT020_NPS)          | // ecal/pt < 0.2 (matches HLT requirement)
  (1ll<<ELEISO_REL015);                      // reliso < 0.15, truncated, 1 GeV EB PS

static const cuts_t electronSelection_el_OSV2 = 
  electronSelection_el_OSV2_iso | electronSelection_el_OSV2_noiso;

static const cuts_t electronSelection_el_OSV2_FO = 
  (1ll<<ELEID_VBTF_90_HLT_CALOIDT_TRKIDVL) | // VBTF90, tightened to match CaloIdT+TrkIdVL
  (1ll<<ELEIP_PV_OSV2_FO)                  | // d0(PV) < 0.2 cm, dz(PV) < 1.0 cm
  (1ll<<ELENOMUON_010)                     | // no muon dR < 0.1
  (1ll<<ELENOTCONV_HITPATTERN)             | // <=1 missing hits
  (1ll<<ELENOTCONV_MIT)                    | // MIT conversion rejection
  (1ll<<ELEPT_010)                         | // electron p_T > 10 GeV
  (1ll<<ELEISO_ECAL_RELNT020_NPS)          | // ecal/pt < 0.2 (matches HLT requirement)
  (1ll<<ELEISO_REL040)                     | // reliso < 0.4, truncated, 1 GeV EB PS
  (1ll<<ELEETA_250);                         // |eta| < 2.5

static const cuts_t electronSelection_el_VBTF90 = 
  (1ll<<ELEID_VBTF_90_HLT);                  // VBTF90, tightened to match CaloIdT+TrkIdVL

static const cuts_t electronSelection_el_VBTF95_NOHOEEND = 
  (1ll<<ELEID_VBTF_95_NOHOEEND);             // VBTF95, no H/E endcap requirement

static const cuts_t electronSelection_el_OSV3_noiso = 
  (1ll<<ELEID_VBTF_90_HLT_CALOIDT_TRKIDVL) | // VBTF90, tightened to match CaloIdT+TrkIdVL
  (1ll<<ELEIP_PV_OSV2)                     | // d0(PV) < 0.04 cm, dz(PV) < 1.0 cm
  (1ll<<ELENOMUON_010)                     | // no muon dR < 0.1
  (1ll<<ELENOTCONV_HITPATTERN)             | // <=1 missing hits
  (1ll<<ELENOTCONV_DISTDCOT002)            | // dist/dcot(theta) conversion rejection
  (1ll<<ELEPT_010)                         | // electron p_T > 10 GeV
  (1ll<<ELE_NOT_TRANSITION)                | // veto electrons with SC in transition region 
  (1ll<<ELEETA_250);                         // |eta| < 2.5

static const cuts_t electronSelection_el_OSV3_iso = 
  (1ll<<ELEISO_ECAL_RELNT020_NPS)          | // ecal/pt < 0.2 (matches HLT requirement)
  (1ll<<ELEISO_RELNT015);                    // reliso < 0.15, non-truncated, 1 GeV EB PS

static const cuts_t electronSelection_el_OSV3 = 
  electronSelection_el_OSV3_iso | electronSelection_el_OSV3_noiso;

static const cuts_t electronSelection_el_OSV3_FO = 
  (1ll<<ELEID_VBTF_90_HLT_CALOIDT_TRKIDVL) | // VBTF90, tightened to match CaloIdT+TrkIdVL
  (1ll<<ELEIP_PV_OSV2_FO)                  | // d0(PV) < 0.2 cm, dz(PV) < 1.0 cm
  (1ll<<ELENOMUON_010)                     | // no muon dR < 0.1
  (1ll<<ELENOTCONV_HITPATTERN)             | // <=1 missing hits
  (1ll<<ELENOTCONV_DISTDCOT002)            | // dist/dcot(theta) conversion rejection
  (1ll<<ELEPT_010)                         | // electron p_T > 10 GeV
  (1ll<<ELEISO_ECAL_RELNT020_NPS)          | // ecal/pt < 0.2 (matches HLT requirement)
  (1ll<<ELEISO_RELNT040)                   | // reliso < 0.4, non-truncated, 1 GeV EB PS
  (1ll<<ELE_NOT_TRANSITION)                | // veto electrons with SC in transition region 
  (1ll<<ELEETA_250);                         // |eta| < 2.5


///////////////
// Higgs, WW //
///////////////

// ======================== WW ============================
//
// The standard WW impact parameter cut
//
//---------------------------------------------------------
static const cuts_t electronSelection_ww_ip =
	 (1ll<<ELEIP_200);
//---------------------------------------------------------
//
// The standard WW isolation cut
//
//---------------------------------------------------------
static const cuts_t electronSelection_ww_iso =
	 (1ll<<ELEISO_REL010);
//---------------------------------------------------------

//--------beginning of WW V0 cuts--------------------------

//---------------------------------------------------------
// WWV0 base cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV0_base  = 
	 (1ll<<ELEETA_250) |
	 (1ll<<ELENOTCONV_DISTDCOT002);

//---------------------------------------------------------
// WWV0 impact parameter cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV0_ip  = 
	 (1ll<<ELEIP_PV_200);

//---------------------------------------------------------
// WWV0 id cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV0_id  = 
	 (1ll<<ELEID_VBTF_35X_80) ;

//---------------------------------------------------------
// WWV0 isolation cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV0_iso  = 
	 (1ll<<ELEISO_REL010_WW); 

//---------------------------------------------------------
// WWV0 selection
//--------------------------------------------------------
static const cuts_t electronSelection_wwV0  = 
	 electronSelection_wwV0_base |
	 electronSelection_wwV0_ip   |
	 electronSelection_wwV0_id   |
	 electronSelection_wwV0_iso;

//---------------------------------------------------------
// WWV0 fakeable object baseline definition
//---------------------------------------------------------
static const cuts_t electronSelectionFO_wwV0_baseline =
	 electronSelection_wwV0_base;

//---------------------------------------------------------
// WWV0 fakeable object definition v1
// extrapolating in isolation and id
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV0_v1 =
	 electronSelectionFO_wwV0_baseline |
	 (1ll<<ELEISO_REL040_WW); 

//---------------------------------------------------------
// WWV0 fakeable object definition v2
// extrapolating in id
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV0_v2 =
	 electronSelectionFO_wwV0_baseline |
	 electronSelection_wwV0_iso;

//---------------------------------------------------------
// WWV0 fakeable object definition v3
// extrapolating in iso
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV0_v3 =
	 electronSelectionFO_wwV0_baseline |
	 electronSelection_wwV0_id |
	 (1ll<<ELEISO_REL100_WW); 

//---------------------------------------------------------
// WWV0 fakeable object definition v4
// extrapolating in iso
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV0_v4 =
	 electronSelectionFO_wwV0_baseline |
	 (1ll<<ELEID_VBTF_35X_90) |
	 (1ll<<ELEISO_REL100_WW); 

//--------end of WW V0 cuts--------------------------------

//--------beginning of WW V0b cuts--------------------------
// At the time of this writing the only difference between V0
// and V0b is that WW V0b additionally requires zero missing
// hits conversion rejection in the electron selection. Conversion
// rejection is in the baseline selection, so the following
// is really an exact copy of the above, except that the
// baseline includes ELENOTCONV_HITPATTERN_0MHITS, and V0 ->
// V0b everywhere. However, I will write the selection out in
// its full glory, rather than just set things equal to the
// above, in recognition that we may choose to evolve V0b into
// a real second generation selection that is more than just
// a copy of V0 plus the 0mhits conversion rejection ;)

//---------------------------------------------------------
// WWV0b base cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV0b_base  = 
	 (1ll<<ELEETA_250) |
	 (1ll<<ELENOTCONV_DISTDCOT002) |
	 (1ll<<ELENOTCONV_HITPATTERN_0MHITS);

//---------------------------------------------------------
// WWV0b impact parameter cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV0b_ip  = 
	 (1ll<<ELEIP_PV_200);

//---------------------------------------------------------
// WWV0b id cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV0b_id  = 
	 (1ll<<ELEID_VBTF_35X_80) ;

//---------------------------------------------------------
// WWV0b isolation cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV0b_iso  = 
	 (1ll<<ELEISO_REL010); 

//---------------------------------------------------------
// WWV0b selection
//--------------------------------------------------------
static const cuts_t electronSelection_wwV0b  = 
	 electronSelection_wwV0b_base |
	 electronSelection_wwV0b_ip   |
	 electronSelection_wwV0b_id   |
	 electronSelection_wwV0b_iso;

//---------------------------------------------------------
// WWV0b fakeable object baseline definition
//---------------------------------------------------------
static const cuts_t electronSelectionFO_wwV0b_baseline =
	 electronSelection_wwV0b_base;

//---------------------------------------------------------
// WWV0b fakeable object definition v1
// extrapolating in isolation and id
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV0b_v1 =
	 electronSelectionFO_wwV0b_baseline |
	 (1ll<<ELEISO_REL040); 

//---------------------------------------------------------
// WWV0b fakeable object definition v2
// extrapolating in id
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV0b_v2 =
	 electronSelectionFO_wwV0b_baseline |
	 electronSelection_wwV0b_iso;

//---------------------------------------------------------
// WWV0b fakeable object definition v3
// extrapolating in iso
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV0b_v3 =
	 electronSelectionFO_wwV0b_baseline |
	 electronSelection_wwV0b_id |
	 (1ll<<ELEISO_REL100); 

//---------------------------------------------------------
// WWV0b fakeable object definition v4
// extrapolating in iso
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV0b_v4 =
	 electronSelectionFO_wwV0b_baseline |
	 (1ll<<ELEID_VBTF_35X_90) |
	 (1ll<<ELEISO_REL100); 

//--------end of WW V0b cuts--------------------------------

//---------------------------------------------------------
// WWV1 base cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV1_base  = 
         (1ll<<ELEETA_250) ;

//---------------------------------------------------------
// WWV1 convrej
//---------------------------------------------------------
static const cuts_t electronSelection_wwV1_convrej  = 
	 (1ll<<ELENOTCONV_DISTDCOT002) | 
       //(1ll<<ELENOTCONV_HITPATTERN39X_0MHITS);
	 (1ll<<ELENOTCONV_HITPATTERN_0MHITS);

//---------------------------------------------------------
// WWV1 impact parameter cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV1_ip  = 
	 (1ll<<ELEIP_PV_wwV1);

//---------------------------------------------------------
// WWV1 id cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV1_id  = 
	 (1ll<<ELEID_VBTF_35X_80) ;

//---------------------------------------------------------
// WWV1 isolation cut
//---------------------------------------------------------
static const cuts_t electronSelection_wwV1_iso  = 
	 (1ll<<ELEISO_REL010); 

//---------------------------------------------------------
// WWV1 selection
//--------------------------------------------------------
static const cuts_t electronSelection_wwV1  = 
	 electronSelection_wwV1_base |
	 electronSelection_wwV1_convrej |
	 electronSelection_wwV1_ip   |
	 electronSelection_wwV1_id   |
	 electronSelection_wwV1_iso;

//---------------------------------------------------------
// WWV1 fakeable object baseline definition
//---------------------------------------------------------
static const cuts_t electronSelectionFO_wwV1_baseline =
	 electronSelection_wwV1_base | 
         electronSelection_wwV1_convrej;


//---------------------------------------------------------
// WWV1 fakeable object definition v1
// extrapolating in isolation and id
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV1_v1 =
	 electronSelectionFO_wwV1_baseline |
	 (1ll<<ELEISO_REL040); 

//---------------------------------------------------------
// WWV1 fakeable object definition v2
// extrapolating in id
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV1_v2 =
	 electronSelectionFO_wwV1_baseline |
	 electronSelection_wwV1_iso;

//---------------------------------------------------------
// WWV1 fakeable object definition v3
// extrapolating in iso
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV1_v3 =
         electronSelectionFO_wwV1_baseline |
	 electronSelection_wwV1_id |
	 (1ll<<ELEISO_REL100); 

//---------------------------------------------------------
// WWV1 fakeable object definition v4
// extrapolating in iso
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_wwV1_v4 =
	 electronSelectionFO_wwV1_baseline |
	 (1ll<<ELEID_VBTF_35X_90) |
	 (1ll<<ELEISO_REL100); 

//--------end of WW V1 cuts--------------------------------

//--------SMURF V1 cuts--------------------------------
/*
static const cuts_t electronSelection_smurfV1_baseline  = 
	 electronSelection_wwV1_base |
	 electronSelection_wwV1_convrej |
	 electronSelection_wwV1_ip;
static const cuts_t electronSelection_smurfV1_iso  = 
         (1ll<<ELEISO_SMURFV1);
static const cuts_t electronSelection_smurfV1_id  = 
	 electronSelection_wwV1_id |
         (1ll<<ELEID_SMURFV1_EXTRA);
static const cuts_t electronSelection_smurfV1  = 
         electronSelection_smurfV1_baseline |
         electronSelection_smurfV1_iso |
         electronSelection_smurfV1_id;
*/
//--------end of SMURF V1 cuts--------------------------------

//--------SMURF V2 cuts--------------------------------
static const cuts_t electronSelection_smurfV2_baseline  = 
	 electronSelection_wwV1_base |
	 electronSelection_wwV1_convrej |
	 electronSelection_wwV1_ip;
static const cuts_t electronSelection_smurfV2_iso  = 
         (1ll<<ELEISO_RELNT010);
static const cuts_t electronSelection_smurfV2_id  = 
	 electronSelection_wwV1_id |
         (1ll<<ELEID_SMURFV2_EXTRA);
static const cuts_t electronSelection_smurfV2  = 
         electronSelection_smurfV2_baseline |
         electronSelection_smurfV2_iso |
         electronSelection_smurfV2_id;
//--------end of SMURF V2 cuts--------------------------------

//--------SMURF V3 cuts--------------------------------
static const cuts_t electronSelection_smurfV3_ip  = 
         (1ll<<ELEIP_PV_SMURFV3);
static const cuts_t electronSelection_smurfV3_baseline  = 
	 electronSelection_wwV1_base |
	 electronSelection_smurfV3_ip;
static const cuts_t electronSelection_smurfV3_convrej  = 
	 (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
         (1ll<<ELENOTCONV_MIT);
static const cuts_t electronSelection_smurfV3_iso  = 
         (1ll<<ELEISO_RELNT010);
static const cuts_t electronSelection_smurfV3_id  = 
	 (1ll<<ELEID_VBTF_80_NOHOEEND) |
         (1ll<<ELEID_SMURFV3_EXTRA);

static const cuts_t electronSelection_smurfV1ss_id  =
         (1ll<<ELEID_VBTF_80_NOHOEEND) |
         (1ll<<ELEID_SMURFV1SS_EXTRA);

static const cuts_t electronSelection_smurfV2ss_id  =
         (1ll<<ELEID_VBTF_80_NOHOEEND) |
         (1ll<<ELEID_SMURFV2SS_EXTRA);

static const cuts_t electronSelection_smurfV3  = 
         electronSelection_smurfV3_baseline |
         electronSelection_smurfV3_convrej |
         electronSelection_smurfV3_iso |
         electronSelection_smurfV3_id;
//--------end of SMURF V3 cuts--------------------------------

//--------SMURF V4 cuts--------------------------------
static const cuts_t electronSelection_smurfV4_ip  = 
         (1ll<<ELEIP_PV_SMURFV3) | (1ll<<ELEIP_PV_DZ_1MM);
static const cuts_t electronSelection_smurfV4_baseline  = 
	 electronSelection_wwV1_base |
	 electronSelection_smurfV4_ip;
static const cuts_t electronSelection_smurfV4_convrej  = 
	 (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
         (1ll<<ELENOTCONV_MIT);
static const cuts_t electronSelection_smurfV4_iso  = 
         (1ll<<ELEISO_SMURFV4);
static const cuts_t electronSelection_smurfV4_id  = 
	 (1ll<<ELEID_VBTF_80_NOHOEEND) |
         (1ll<<ELEID_SMURFV3_EXTRA);
static const cuts_t electronSelection_smurfV4  = 
         electronSelection_smurfV4_baseline |
         electronSelection_smurfV4_convrej |
         electronSelection_smurfV4_iso |
         electronSelection_smurfV4_id;
//--------end of SMURF V4 cuts--------------------------------

//--------SMURF V5 cuts--------------------------------
static const cuts_t electronSelection_smurfV5_ip  = 
         (1ll<<ELEIP_PV_SMURFV3) | (1ll<<ELEIP_PV_DZ_1MM);
static const cuts_t electronSelection_smurfV5_baseline  = 
	 electronSelection_wwV1_base |
	 electronSelection_smurfV5_ip;
static const cuts_t electronSelection_smurfV5_convrej  = 
	 (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
         (1ll<<ELENOTCONV_MIT);
static const cuts_t electronSelection_smurfV5_iso  = 
         (1ll<<ELEISO_SMURFV5);
static const cuts_t electronSelection_smurfV5_id  = 
	 (1ll<<ELEID_VBTF_80_NOHOEEND) |
         (1ll<<ELEID_SMURFV3_EXTRA);
static const cuts_t electronSelection_smurfV5  = 
         electronSelection_smurfV5_baseline |
         electronSelection_smurfV5_convrej |
         electronSelection_smurfV5_iso |
         electronSelection_smurfV5_id;
static const cuts_t electronSelection_smurfV6  = 
         electronSelection_smurfV5;
//--------end of SMURF V5 cuts--------------------------------

//--------SMURF FakableObject cuts--------------------------------
static const cuts_t electronSelectionFO_el_smurf_base =
  (1ll<<ELEETA_250) | (1ll<<ELEIP_PV_DZ_1MM) |
  electronSelection_smurfV3_convrej;
//---------------------------------------------------------
// Fakeable object definition (option V3)
// extrapolating in isolation as much as the trigger allows
// *! USE WITH CARE !*
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_smurf_v3 =
  electronSelectionFO_el_smurf_base |
  (1ll<<ELEID_VBTF_80_NOHOEEND);
//---------------------------------------------------------
// Fakeable object definition (option V1)
// extrapolating in isolation as much as the trigger allows
// and in partial id
// *! USE WITH CARE !*
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_smurf_v1 =
  electronSelectionFO_el_smurf_base |
  (1ll<<ELEID_VBTF_90_HLT);
//---------------------------------------------------------
// Fakeable object definition (option V4)
// extrapolating in partial id and partial isolation
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_smurf_v4 =
  electronSelectionFO_el_smurf_base |
  (1ll<<ELEID_VBTF_90_HLT) |
  (1ll<<ELEISO_TRK_RELNT020) |
  (1ll<<ELEISO_ECAL_RELNT020) |
  (1ll<<ELEISO_HCAL_RELNT020);
//---------------------------------------------------------
// Fakeable object definition (option V2)
// extrapolating in partial id
//---------------------------------------------------------
static const cuts_t electronSelectionFO_el_smurf_v2 =
  electronSelectionFO_el_smurf_base |
  (1ll<<ELEID_VBTF_90_HLT) |
  (1ll<<ELEISO_RELNT010); 
//--------end of SMURF FakableObject cuts------------------

// ======================== SS ============================
//
// The standard SS selection
//---------------------------------------------------------
// SS NoIso selections


/////////////////////////////////////
// 2011 SS Selections              //
/////////////////////////////////////

//baseline SS FO selections
static const cuts_t electronSelectionFO_SS_baseline =
           (1ll<<ELEETA_240) |
           (1ll<<ELENOMUON_010);

// Analysis Selection (fake rate numerator)
static const cuts_t electronSelection_ssV3_noIso = 
           electronSelectionFO_SS_baseline    |
           electronSelection_smurfV3_convrej  |
           electronSelection_smurfV1ss_id     |
           (1ll<<ELEIP_SS200)                 |
           (1ll<<ELECHARGE_NOTFLIP3AGREE);

static const cuts_t electronSelection_ssV3_iso =
                 (1ll<<ELEISO_RELNT015) |
                 (1ll<<ELEISO_ECAL_RELNT020_NPS);

static const cuts_t electronSelection_ssV3 = 
                       electronSelection_ssV3_noIso |
                       electronSelection_ssV3_iso;


// Loose "Fakeable Object" Selection (fake rate denominators)

static const cuts_t electronSelectionFOV3_ssVBTF80_v1 =       // V1 - relaxed Id & Isolation
                 electronSelectionFO_SS_baseline    |
                 electronSelection_smurfV3_convrej  |                
                 (1ll<<ELECHARGE_NOTFLIP3AGREE)     |
                 (1ll<<ELEISO_ECAL_RELNT020_NPS)         |
                 (1ll<<ELEISO_HCAL_RELNT020);

static const cuts_t electronSelectionFOV3_ssVBTF80_v2 =       // V2 - relaxed Id
                 electronSelectionFO_SS_baseline    |
                 electronSelection_smurfV3_convrej  |                
                 (1ll<<ELECHARGE_NOTFLIP3AGREE)     |
                 (1ll<<ELEISO_RELNT015)             |
                 (1ll<<ELEISO_ECAL_RELNT020_NPS);

static const cuts_t electronSelectionFOV3_ssVBTF80_v3 =       // V3 - relaxed isolation 
                 electronSelectionFO_SS_baseline    |
                 electronSelection_smurfV3_convrej  |
                 electronSelection_smurfV1ss_id     |              
                 (1ll<<ELECHARGE_NOTFLIP3AGREE)     |
                 (1ll<<ELEISO_ECAL_RELNT020_NPS)    |
                 (1ll<<ELEISO_HCAL_RELNT020);


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// Analysis Selection (fake rate numerator)
static const cuts_t electronSelection_ssV4_noIso = 
           electronSelectionFO_SS_baseline     |
           (1ll<<ELENOTCONV_DISTDCOT002)       |
           (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
           electronSelection_smurfV1ss_id      |
           (1ll<<ELEIP_SS200)                  |
           (1ll<<ELESEED_ECAL)                 |
           (1ll<<ELE_NOT_TRANSITION)           |
           (1ll<<ELECHARGE_NOTFLIP3AGREE);

static const cuts_t electronSelection_ssV4_iso =
                 (1ll<<ELEISO_RELNT015);

static const cuts_t electronSelection_ssV4 = 
                       electronSelection_ssV4_noIso |
                       electronSelection_ssV4_iso;

// Loose "Fakeable Object" Selection (fake rate denominators)

static const cuts_t electronSelectionFOV4_ssVBTF80_v1 =       // V1 - relaxed Id & Isolation
                 electronSelectionFO_SS_baseline     |
                 (1ll<<ELE_NOT_TRANSITION)           |
                 (1ll<<ELENOTCONV_DISTDCOT002)       |
                 (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
                 (1ll<<ELECHARGE_NOTFLIP3AGREE)      |
                 (1ll<<ELESEED_ECAL);

static const cuts_t electronSelectionFOV4_ssVBTF80_v2 =       // V2 - relaxed Id
                 electronSelectionFO_SS_baseline     |
                 (1ll<<ELE_NOT_TRANSITION)           |
                 (1ll<<ELENOTCONV_DISTDCOT002)       |
                 (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
                 (1ll<<ELECHARGE_NOTFLIP3AGREE)      |
                 (1ll<<ELESEED_ECAL)                 |
                 (1ll<<ELEISO_RELNT015);

static const cuts_t electronSelectionFOV4_ssVBTF80_v3 =       // V3 - relaxed isolation (relaxed all the way; we store the relIso and can cut on it separately in the babies or elsewhere)
                 electronSelectionFO_SS_baseline     |
                 (1ll<<ELE_NOT_TRANSITION)           |
                 (1ll<<ELENOTCONV_DISTDCOT002)       |
                 (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
                 electronSelection_smurfV1ss_id      |              
                 (1ll<<ELECHARGE_NOTFLIP3AGREE)      |
                 (1ll<<ELESEED_ECAL);

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

//baseline SS FO selections
static const cuts_t electronSelectionFO_SS_baselineV2 =
           (1ll<<ELEETA_240) |
           (1ll<<ELENOMUON_010_SS);

// Analysis Selection (fake rate numerator)
static const cuts_t electronSelection_ssV5_noIso = 
           electronSelectionFO_SS_baselineV2   |
           (1ll<<ELENOTCONV_DISTDCOT002)       |
           (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
           electronSelection_smurfV1ss_id      |
           (1ll<<ELEIP_SS200)                  |
           (1ll<<ELESEED_ECAL)                 |
           (1ll<<ELE_NOT_TRANSITION)           |
           (1ll<<ELECHARGE_NOTFLIP3AGREE);

static const cuts_t electronSelection_ssV5_iso =
                 (1ll<<ELEISO_RELNT015);

static const cuts_t electronSelection_ssV5 = 
                       electronSelection_ssV5_noIso |
                       electronSelection_ssV5_iso;

// Loose "Fakeable Object" Selection (fake rate denominators)

static const cuts_t electronSelectionFOV5_ssVBTF80_v1 =       // V1 - relaxed Id & Isolation
                 electronSelectionFO_SS_baselineV2   |
                 (1ll<<ELE_NOT_TRANSITION)           |
                 (1ll<<ELENOTCONV_DISTDCOT002)       |
                 (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
                 (1ll<<ELECHARGE_NOTFLIP3AGREE)      |
                 (1ll<<ELESEED_ECAL);

static const cuts_t electronSelectionFOV5_ssVBTF80_v2 =       // V2 - relaxed Id
                 electronSelectionFO_SS_baselineV2   |
                 (1ll<<ELE_NOT_TRANSITION)           |
                 (1ll<<ELENOTCONV_DISTDCOT002)       |
                 (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
                 (1ll<<ELECHARGE_NOTFLIP3AGREE)      |
                 (1ll<<ELESEED_ECAL)                 |
                 (1ll<<ELEISO_RELNT015);

static const cuts_t electronSelectionFOV5_ssVBTF80_v3 =       // V3 - relaxed isolation (relaxed all the way; we store the relIso and can cut on it separately in the babies or elsewhere)
                 electronSelectionFO_SS_baselineV2   |
                 (1ll<<ELE_NOT_TRANSITION)           |
                 (1ll<<ELENOTCONV_DISTDCOT002)       |
                 (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
                 electronSelection_smurfV1ss_id      |              
                 (1ll<<ELECHARGE_NOTFLIP3AGREE)      |
                 (1ll<<ELESEED_ECAL);

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// Analysis Selection (fake rate numerator)
static const cuts_t electronSelection_ssV5_noIso_noConvCuts = 
           electronSelectionFO_SS_baselineV2   |
           electronSelection_smurfV1ss_id      |
           (1ll<<ELEIP_SS200)                  |
           (1ll<<ELESEED_ECAL)                 |
           (1ll<<ELE_NOT_TRANSITION)           |
           (1ll<<ELECHARGE_NOTFLIP3AGREE);

static const cuts_t electronSelection_ssV5_iso_noConvCuts =
                 (1ll<<ELEISO_RELNT015);

static const cuts_t electronSelection_ssV5_noConvCuts = 
                       electronSelection_ssV5_noIso_noConvCuts |
                       electronSelection_ssV5_iso_noConvCuts;

// Loose "Fakeable Object" Selection (fake rate denominators)

static const cuts_t electronSelectionFOV5_ssVBTF80_noConvCuts_v1 =       // V1 - relaxed Id & Isolation
                 electronSelectionFO_SS_baselineV2   |
                 (1ll<<ELE_NOT_TRANSITION)           |
                 (1ll<<ELECHARGE_NOTFLIP3AGREE)      |
                 (1ll<<ELESEED_ECAL);

static const cuts_t electronSelectionFOV5_ssVBTF80_noConvCuts_v2 =       // V2 - relaxed Id
                 electronSelectionFO_SS_baselineV2   |
                 (1ll<<ELE_NOT_TRANSITION)           |
                 (1ll<<ELECHARGE_NOTFLIP3AGREE)      |
                 (1ll<<ELESEED_ECAL)                 |
                 (1ll<<ELEISO_RELNT015);

static const cuts_t electronSelectionFOV5_ssVBTF80_noConvCuts_v3 =       // V3 - relaxed isolation (relaxed all the way; we store the relIso and can cut on it separately in the babies or elsewhere)
                 electronSelectionFO_SS_baselineV2   |
                 (1ll<<ELE_NOT_TRANSITION)           |
                 electronSelection_smurfV1ss_id      |              
                 (1ll<<ELECHARGE_NOTFLIP3AGREE)      |
                 (1ll<<ELESEED_ECAL);

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// Analysis Selection (fake rate numerator)
static const cuts_t electronSelection_ssV6_noIso = 
           electronSelectionFO_SS_baselineV2   |
           electronSelection_smurfV1ss_id      |
           (1ll<<ELENOTCONV_DISTDCOT002_OLD)   |
           (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
           (1ll<<ELEIP_SS200)                  |
           (1ll<<ELESEED_ECAL)                 |
           (1ll<<ELE_NOT_TRANSITION)           |
           (1ll<<ELECHARGE_NOTFLIP3AGREE);

static const cuts_t electronSelection_ssV6_iso =
                 (1ll<<ELEISO_RELNT015);

static const cuts_t electronSelection_ssV6 = 
                       electronSelection_ssV6_noIso |
                       electronSelection_ssV6_iso;

// Loose "Fakeable Object" Selection (fake rate denominators)

static const cuts_t electronSelectionFOV6_ssVBTF80_v1 =       // V1 - relaxed Id & Isolation
                 electronSelectionFO_SS_baselineV2   |
                 (1ll<<ELENOTCONV_DISTDCOT002_OLD)   |
                 (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
                 (1ll<<ELE_NOT_TRANSITION)           |
                 (1ll<<ELECHARGE_NOTFLIP3AGREE)      |
                 (1ll<<ELESEED_ECAL);

static const cuts_t electronSelectionFOV6_ssVBTF80_v2 =       // V2 - relaxed Id
                 electronSelectionFO_SS_baselineV2   |
                 (1ll<<ELENOTCONV_DISTDCOT002_OLD)   |
                 (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
                 (1ll<<ELE_NOT_TRANSITION)           |
                 (1ll<<ELECHARGE_NOTFLIP3AGREE)      |
                 (1ll<<ELESEED_ECAL)                 |
                 (1ll<<ELEISO_RELNT015);

static const cuts_t electronSelectionFOV6_ssVBTF80_v3 =       // V3 - relaxed isolation (relaxed all the way; we store the relIso and can cut on it separately in the babies or elsewhere)
                 electronSelectionFO_SS_baselineV2   |
                 electronSelection_smurfV1ss_id      |
                 (1ll<<ELENOTCONV_DISTDCOT002_OLD)   |
                 (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
                 (1ll<<ELE_NOT_TRANSITION)           |
                 (1ll<<ELECHARGE_NOTFLIP3AGREE)      |
                 (1ll<<ELESEED_ECAL);

/////////////////////////////////////
// End 2011 SS Selections          //
/////////////////////////////////////


/////////////////////////////////////
// 2012 SS Selections              //
/////////////////////////////////////

// Analysis Selection (fake rate numerator)
static const cuts_t electronSelection_ssV7_noIso = 
                       electronSelectionFO_SS_baselineV2   |
                       (1ll<<ELEID_WP2012_MEDIUM_NOISO)    |
                       (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
                       (1ll<<ELE_NOT_TRANSITION)           |
                       (1ll<<ELECHARGE_NOTFLIP3AGREE);

static const cuts_t electronSelection_ssV7_iso =
                       (1ll<<ELEISO_RELNT015);

static const cuts_t electronSelection_ssV7 = 
                       electronSelection_ssV7_noIso |
                       electronSelection_ssV7_iso;

// Loose "Fakeable Object" Selection (fake rate denominators)

static const cuts_t electronSelectionFOV7_v1 =       // V1 - relaxed Id & Isolation
                 electronSelectionFO_SS_baselineV2   |
                 (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
                 (1ll<<ELE_NOT_TRANSITION)           |
                 (1ll<<ELECHARGE_NOTFLIP3AGREE);

static const cuts_t electronSelectionFOV7_v2 =       // V2 - relaxed Id
                 electronSelectionFO_SS_baselineV2   |
                 electronSelection_ssV7_iso          |
                 (1ll<<ELENOTCONV_HITPATTERN_0MHITS) |
                 (1ll<<ELE_NOT_TRANSITION)           |
                 (1ll<<ELECHARGE_NOTFLIP3AGREE);

static const cuts_t electronSelectionFOV7_v3 =       // V3 - relaxed isolation (relaxed all the way; we store the relIso and can cut on it separately in the babies or elsewhere)
                 electronSelectionFO_SS_baselineV2     |
                 (1ll<<ELEID_WP2012_MEDIUM_NOISO_NOIP) |
                 (1ll<<ELENOTCONV_HITPATTERN_0MHITS)   |
                 (1ll<<ELE_NOT_TRANSITION)             |
                 (1ll<<ELECHARGE_NOTFLIP3AGREE);

/////////////////////////////////////
// End 2012 SS Selections          //
/////////////////////////////////////

// /////////////////////////////////////
// // 2012 TTV Selections             //
// /////////////////////////////////////

// note, may need to change dz and d0 cuts
static const cuts_t electronSelection_TTVTightv1_noIso =
					   //					   (1ll<<ELEETA_240)                      |
					   (1ll<<ELEID_WP2012_MEDIUM_NOISO)       |
                       (1ll<<ELENOTCONV_HITPATTERN_0MHITS)    |
					   (1ll<<ELE_NOT_TRANSITION);
					   

//fakeable object definition
static const cuts_t electronSelection_TTVTightFOv1 =
							 //					   (1ll<<ELEETA_240)                     |
					   (1ll<<ELEID_WP2012_MEDIUM_NOISO_NOIP)  |
					   (1ll<<ELENOTCONV_HITPATTERN_0MHITS)    |
					   (1ll<<ELE_NOT_TRANSITION);              



// /////////////////////////////////////
// // End 2012 TTV Selections         //
// /////////////////////////////////////











////////////////////////////
// enums for internal use //
////////////////////////////

enum EgammaFiduciality {
  ISEB,
  ISEBEEGAP,
  ISEE,
  ISEEGAP,
  ISEBETAGAP,
  ISEBPHIGAP,
  ISEEDEEGAP,
  ISEERINGGAP,
  ISGAP
};

enum EgammaElectronType {
  ISECALENERGYCORRECTED,  // if false, the electron "ecalEnergy" is just the supercluster energy 
  ISMOMENTUMCORRECTED,    // has E-p combination been applied
  ISECALDRIVEN,
  ISTRACKERDRIVEN
};

enum ElectronIDComponent {
  ELEID_ID,
  ELEID_ISO,
  ELEID_CONV,
  ELEID_IP,
};

namespace wp2012 {
enum ElectronIDComponentWP2012 {
    DETAIN          = (1<<0),
    DPHIIN          = (1<<1),
    SIGMAIETAIETA   = (1<<2),
    HOE             = (1<<3),
    OOEMOOP         = (1<<4),
    D0VTX           = (1<<5),
    DZVTX           = (1<<6),
    ISO             = (1<<7),
    VTXFIT          = (1<<8),
    MHITS           = (1<<9)
};
static const electronIdComponent_t PassAllWP2012Cuts = DETAIN | DPHIIN | SIGMAIETAIETA | HOE 
                                | OOEMOOP | D0VTX | DZVTX | ISO | VTXFIT | MHITS;

static const electronIdComponent_t PassWP2012CutsNoIso = DETAIN | DPHIIN | SIGMAIETAIETA | HOE 
                                | OOEMOOP | D0VTX | DZVTX | VTXFIT | MHITS;

static const electronIdComponent_t PassWP2012CutsIso = ISO;

static const electronIdComponent_t PassWP2012CutsNoIsoNoIP = DETAIN | DPHIIN | SIGMAIETAIETA | HOE 
                                | OOEMOOP | DZVTX | VTXFIT | MHITS;
};


// master selection function
bool pass_electronSelectionCompareMask(const cuts_t cuts_passed, const cuts_t selectionType);
bool pass_electronSelection(const unsigned int index, const cuts_t selectionType, bool applyAlignmentCorrection = false, bool removedEtaCutInEndcap = false, bool useGsfTrack = true);
cuts_t electronSelection(const unsigned int index, bool applyAlignmentCorrection = false, bool removedEtaCutInEndcap = false, bool useGsfTrack = true);

// "smurf" electron id
// WARNING!!! this is not the full smurf selection, just the additional ID on top of VBTF80 for low pt guys
bool electronId_smurf_v1(const unsigned int index);
bool electronId_smurf_v2(const unsigned int index);
bool electronId_smurf_v3(const unsigned int index);
bool electronId_smurf_v1ss(const unsigned int index);
bool electronId_smurf_v2ss(const unsigned int index);

// "cand" electron id
//bool electronId_cand(const unsigned int index, const cand_tightness tightness, bool applyAlignementCorrection = false, bool removedEtaCutInEndcap = false);
//bool electronId_extra(const unsigned int index);

// WP2012
electronIdComponent_t electronId_WP2012(const unsigned int index, const wp2012_tightness tightness);
electronIdComponent_t electronId_WP2012_v2(const unsigned int index, const wp2012_tightness tightness, bool useOldIsolation = false);
electronIdComponent_t electronId_WP2012_v3(const unsigned int index, const wp2012_tightness tightness, bool useOldIsolation = false);
electronIdComponent_t electronId_WP2012_noIso_useElEtaForIsEB(const unsigned int index, const wp2012_tightness tightness);  // same as v2 except uses el->eta() to determine isEB and no iso

// "VBTF" id
electronIdComponent_t electronId_VBTF(const unsigned int index, const vbtf_tightness tightness,  bool applyAlignementCorrection = false, bool removedEtaCutInEndcap = false);

// "CIC" id
electronIdComponent_t electronId_CIC(const unsigned int index, const unsigned int version, const cic_tightness tightness, bool applyAlignementCorrection = false, bool removedEtaCutInEndcap =false);

unsigned int eidClassify(const unsigned int version, const unsigned int index);
bool eidComputeCut(double x, double et, double cut_min, double cut_max, bool gtn=false);

//electron ID WP as of June 8th 2011
electronIdComponent_t passLikelihoodId(unsigned int index, float lhValue, int workingPoint);
//electron ID WP as of September 8th 2011
bool passLikelihoodId_v2(unsigned int index, float lhValue, int workingPoint);

// relative isolation 
// - standard track isolation from CMSSW
// --- CMSSW >= 3_5_X track jurassic strip half width 0.015
// --- CMSSW < 3_5_X no jurassic iso for tracks
// - ecal pedestal of 1 GeV subtracted in EB as Max(0, ecalIso - 1.0)
// - hcal iso as usual
float electronIsolation_rel(const unsigned int index, bool use_calo_iso);
//non-truncated relative iso
//float electronIsolation_rel_v1Original(const unsigned int index, bool use_calo_iso);
float electronIsolation_rel_v1(const unsigned int, bool);
float electronIsolation_ECAL_rel_v1(const unsigned int, bool useEBps = true);
float electronIsolation_HCAL_rel_v1(const unsigned int);
float electronIsolation_ECAL_rel(const unsigned int);
float electronIsolation_HCAL_rel(const unsigned int);
float electronIsolation_rel_FastJet(const unsigned int, bool);
float electronIsolation_rel_v1_FastJet(const unsigned int, bool);
float el_fastjet_rel_offset(const unsigned int);
float electronIsolation_rel_ww(const unsigned int index, bool use_calo_iso);  // the difference from above is that the pedestal sub is applied on both EB/EE
float electronIsoValuePF(const unsigned int iel, unsigned int ivtx, float coner=0.4, float minptn=1.0, float dzcut=0.1, 
			 float footprintdr=0.07, float gammastripveto=0.025, float elestripveto=-999., int filterId = 0);
float electronIsolation_cor_rel_v1(const unsigned int, bool);

// remove electrons that are overlapping with a muon
bool electronId_noMuon(const unsigned int index);
bool electronId_noMuon_SS(const unsigned int index);

// conversion rejection
bool isFromConversionHitPattern(const unsigned int index);
bool isFromConversionPartnerTrack(const unsigned int index);
bool isFromConversionPartnerTrack_v2(const unsigned int index);
bool isFromConversionMIT(const unsigned int index);

//electron charge using the majority logic of the egamma group
int getChargeUsingMajorityLogic(int elIdx, float minFracSharedHits = 0.45);

//charge flip rejection
bool isChargeFlip(int elIndex);
bool isChargeFlip3agree(int elIndex); 

// spike rejection for electrons
bool isSpikeElectron(const unsigned int index);

// position correction for electrons
void electronCorrection_pos(const unsigned int index, float &dEtaIn, float &dPhiIn);

// d0 corrected by the primary vertex
double electron_d0PV(unsigned int index);
double electron_d0PV_wwV1(unsigned int index);
double electron_d0PV_mindz(unsigned int index);
double electron_dzPV_wwV1(unsigned int index);
double electron_d0PV_smurfV3(unsigned int index);
double electron_dzPV_smurfV3(unsigned int index);

// dz
double dzPV(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv);

//
// 2012 PF Iso
//

void electronIsoValuePF2012(float &pfiso_ch, float &pfiso_em, float &pfiso_nh, const float R, const unsigned int iel, const int ivtx, bool barrelVetoes = false);
void electronIsoValuePF2012reco(float &pfiso_ch, float &pfiso_em, float &pfiso_nh, const float R, const unsigned int iel, const int ivtx, float neutral_threshold = 0.5);
float electronIsoValuePF2012_FastJetEffArea( int index , float conesize = 0.3, int ivtx = 0);
float electronIsoValuePF2012_FastJetEffArea_v2(int index, float conesize = 0.3, int ivtx = 0, bool useOldIsolation = false);
float electronIsoValuePF2012_FastJetEffArea_v3(int index, float conesize = 0.3, int ivtx = 0, bool useOldIsolation = false);
float electronRadialIsolation(int index, float &chiso, float &nhiso, float &emiso, float neutral_et_threshold = 1.0, float cone_size = 0.3, bool barrelVetoes = false, bool verbose = false);

float electronIsoValuePF2012_FastJetEffArea_HWW( int index );

// calculate Effective area (updated to value from Egamma)
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaEARhoCorrection
// Topic revision: r12 - 28-Nov-2012
float fastJetEffArea03_v2(const float eta);
float fastJetEffArea04_v2(const float eta);

float fastJetEffArea03_v1(const float eta);
float fastJetEffArea04_v1(const float eta);

//
// 2012 cut based ID
//





#endif

