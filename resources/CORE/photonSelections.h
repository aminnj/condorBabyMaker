#ifndef PHOTON_SELECTIONS_H
#define PHOTON_SELECTIONS_H

// photon id
// see https://twiki.cern.ch/twiki/bin/viewauth/CMS/PhotonID

// photon selection choices
enum PhotonSelectionType { 
  Yuri = 0
}; 

// declarations
bool photonId( const unsigned int, PhotonSelectionType = Yuri );
bool isSpikePhoton( const unsigned int );
int  isGoodEMObject( const unsigned int index );
bool isGoodEMObject2012( const unsigned int index );
bool photon_VGamma_2011(const int index);

#endif

