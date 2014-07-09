#ifndef MUON_SELECTIONS_H
#define MUON_SELECTIONS_H

///////////////
// Selectors //
///////////////

enum SelectionType { 

    ///////////////////
    // Opposite Sign //
    ///////////////////

      // OS Analysis
      OSGeneric_v4,
      OSGeneric_v3,

      // OS Fakes
      OSGeneric_v3_FO,

      // OSZ
      OSZ_v4,
      OSZ_v3,
      OSZ_v2,

    ////////////////////
    // Same Sign 2011 //
    ///////////////////

      // Analysis
      NominalSSv4,
   
      // Fakes
      muonSelectionFO_ssV4,

    ////////////////////
    // Same Sign 2012 //
    ///////////////////

      // Analysis
      NominalSSv5,
   
      // Fakes
      muonSelectionFO_ssV5,

    ///////////////
    // Higgs, WW //
    ///////////////

      // WW
      NominalWWV0,
      muonSelectionFO_mu_ww,
      muonSelectionFO_mu_ww_iso10,
      NominalWWV1,
      muonSelectionFO_mu_wwV1,
      muonSelectionFO_mu_wwV1_iso10,
      muonSelectionFO_mu_wwV1_iso10_d0,

      // Analysis
      NominalSmurfV3,
      NominalSmurfV4,
      NominalSmurfV5,
      NominalSmurfV6,

      // Fakes
      muonSelectionFO_mu_smurf_04,
      muonSelectionFO_mu_smurf_10,
      ZMet2012_v1,
      ZMet2012_detiso_v1,

    ///////////////
    // TTV 2012  //
    ///////////////
      
      // Analysis
      NominalTTZ_loose_v1,
      NominalTTZ_tight_v1,

      // Fakes
      NominalTTZ_looseFO_v1,
      NominalTTZ_tightFO_v1
}; 

////////////////////
// Identification //
////////////////////

bool muonId           ( unsigned int index, SelectionType type);
bool muonIdNotIsolated( unsigned int index, SelectionType type);
bool isGoodStandardMuon( unsigned int index );

///////////////
// Isolation //
///////////////
double muonIsoValue          ( unsigned int , bool = true );
double muonIsoValue_TRK      ( unsigned int , bool = true );
double muonIsoValue_ECAL     ( unsigned int , bool = true );
double muonIsoValue_HCAL     ( unsigned int , bool = true );
double muonIsoValue_FastJet  ( unsigned int , bool = true );
double mu_fastjet_rel_offset ( unsigned int , bool = true );

double muonIsoValuePF        ( unsigned int imu, unsigned int ivtx, float coner=0.4, float minptn=1.0, float dzcut=0.1, int filterId = 0);
void muonIsoValuePF2012  (float &pfiso_ch, float &pfiso_em, float &pfiso_nh, const float R, const unsigned int imu, const int ivtx, float neutral_et_threshold = 0.5);
float muonIsoValuePF2012_FastJetEffArea( int index , float conesize, float effective_area, int ivtx);
double muonCorIsoValue (unsigned int , bool = true);
float muonRadialIsolation (unsigned int imu, float &chiso, float &nhiso, float &emiso, float neutral_et_threshold = 1.0, float cone_size = 0.3, bool verbose = false);
float muonIsoValuePF2012_deltaBeta(unsigned int imu);

///////////////////////
// Cosmics Rejection //
///////////////////////

bool isCosmics(unsigned int index);

/////////////////////////////
// Muon d0 corrected by PV //
/////////////////////////////

double mud0PV         (unsigned int index);
double mud0PV_wwV1    (unsigned int index);
double mudzPV_wwV1    (unsigned int index);
double mud0PV_smurfV3 (unsigned int index);
double mudzPV_smurfV3 (unsigned int index);

/////////////////////////////////////////////////////////////////
// checks if muon is also pfmuon, and pfmuon pt = reco muon pt //
/////////////////////////////////////////////////////////////////

bool isPFMuon( int index , bool requireSamePt = true , float dpt_max = 1.0 );

struct mu2012_tightness 
{
    enum value_type 
    {
        LOOSE,
        TIGHT,
        static_size
    };
};

bool passes_muid_wp2012(const unsigned int index, const mu2012_tightness::value_type tightness);

#endif

