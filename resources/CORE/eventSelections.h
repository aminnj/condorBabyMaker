#ifndef EVENTSELECTIONS_H
#define EVENTSELECTIONS_H
#include <sys/types.h>
//----------------------------------------------------------------
// A ridicolusly simple function, but since the Z veto is used 
// in two places, might as well centralize it to keep consistency
//----------------------------------------------------------------
bool inZmassWindow (float mass);

//----------------------------------------------------------------
// standard event cleaning
// for data
//----------------------------------------------------------------
bool cleaning_standard(bool isData);

//----------------------------------------------------------------
// 20 October 2010
// standard event cleaning used for SS analysis
//----------------------------------------------------------------
bool cleaning_standardOctober2010();


//----------------------------------------------------------------
// 26 April 2011
// standard event cleaning used for OS analysis
//----------------------------------------------------------------
bool cleaning_standardApril2011();

//----------------------------------------------------------------
// 04 November 2011
// standard event cleaning used for SS analysis
//----------------------------------------------------------------
bool cleaning_standardNovember2011();

////----------------------------------------------------------------
//// 5 August 2010
//// standard event cleaning
//// for low pt dilepton / fake rate data studies
////----------------------------------------------------------------
//bool cleaning_standardAugust2010(bool isdata);

//----------------------------------------------------------------
// standard event cleaning
// for low pt dilepton / fake rate data studies
//----------------------------------------------------------------
bool cleaning_standardNoBSC(bool isData);

//----------------------------------------------------------------
// require BPTX
//----------------------------------------------------------------
bool cleaning_BPTX(bool isData);

//----------------------------------------------------------------
// require bit 40 or 41 passed
//----------------------------------------------------------------
bool cleaning_BSC();

//----------------------------------------------------------------
// require bits 36-39 DIDN't pass
// to reject beam halo
//----------------------------------------------------------------
bool cleaning_beamHalo();

// ----------------------------------------------------------------
// 5 August 2010
// at least 1 good vertex
// z position increased from 15 to 24 cm
// ----------------------------------------------------------------
bool cleaning_goodVertexAugust2010();

// ----------------------------------------------------------------
// 26 April 2011
// at least 1 good vertex
// ----------------------------------------------------------------
bool cleaning_goodVertexApril2011();

//----------------------------------------------------------------
// if >= 10 tracks, require at least 25% high purity
//----------------------------------------------------------------
bool cleaning_goodTracks();

//----------------------------------------------------------------
// checks whether a vertex is good or not
//----------------------------------------------------------------
bool isGoodVertex(size_t ivtx);

//----------------------------------------------------------------
// checks whether the leptons of a given
// hypothesis come from the same good vertex
// by checking if both leptons are within dz
// of 1cm of the same PV
//----------------------------------------------------------------
bool hypsFromSameVtx(size_t hypIdx);

//----------------------------------------------------------------
// checks whether the leptons of a given
// hypothesis come from the same good vertex
// by checking if both leptons are within dz
// of 0.2 cm of the same PV and if that PV is
// the closest vertex to each lepton
//----------------------------------------------------------------
int hypsFromSameVtx2011(size_t hypIdx, float dz = 0.2, bool requireClosest = false);

// find first good vertex
int firstGoodVertex ();

//----------------------------------------------------------------
// checks whether the leptons of a given
// hypothesis come from the same good vertex
// by checking if both leptons are within dz
// of 1cm of the same PV
//----------------------------------------------------------------
bool hypsFromFirstGoodVertex(size_t hypIdx, float dz_cut = 1.0);

/*****************************************************************************************/
// number of good vertices in the event
/*****************************************************************************************/
int numberOfGoodVertices(void);

//

int chargedHadronVertex( const unsigned int );

//----------------------------------------------------------------
// These are the MET filters described here:
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters
//----------------------------------------------------------------
bool passCSCBeamHaloFilter();
bool passHBHEFilter();
bool passHCALLaserFilter();
bool passECALDeadCellFilter();
bool passTrackingFailureFilter();
bool passeeBadScFilter();
bool passECALLaserFilter();
bool passMETFilters();


#endif

