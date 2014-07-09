#include "Math/VectorUtil.h"
#include "IsoTrackVeto.h"
#include "../eventSelections.h"
#include "../trackSelections.h"
#include "../CMS2.h"

float TrackIso(int thisPf){

  const bool dovtxcut=false;

  // no cuts  
  float iso_dr03_dz01_pt00 = 0.0;
  
  for (int ipf = 0; ipf < (int)cms2.pfcands_p4().size(); ipf++) {

    if( ipf == thisPf ) continue; // skip this PFCandidate
    
    if(cms2.pfcands_charge().at(ipf) == 0 ) continue; // skip neutrals                                                                                                                          

    // when we want to find the second isolated muon and electron
    // we do not use the electron and muon in the isolation sum,                                                                                              
    // to avoid overlap with the other lepton in the event  
    //if((abs(cms2.pfcands_particleId().at(thisPf))==13 || abs(cms2.pfcands_particleId().at(thisPf))==11) && abs(cms2.pfcands_particleId().at(ipf))==13) continue;    
    //if((abs(cms2.pfcands_particleId().at(thisPf))==13 || abs(cms2.pfcands_particleId().at(thisPf))==11) && abs(cms2.pfcands_particleId().at(ipf))==11) continue;    


    //----------------------------------------                                                                                                                                                   
    // find closest PV and dz w.r.t. that PV                                                                                                                                                     
   //----------------------------------------                                                                                                                                                   

    float mindz = 999.;
    int vtxi    = -1; 

    if (dovtxcut) {
      for (int ivtx = 0; ivtx < (int)cms2.vtxs_position().size(); ivtx++) {

        int itrk = cms2.pfcands_trkidx().at(ipf);

        if( itrk >= (int)cms2.trks_trk_p4().size() || itrk < 0 ){
          //note: this should only happen for electrons which do not have a matched track                                                                                                            
          //currently we are just ignoring these guys                                                                                                                                                
          continue;
        }

        if(!isGoodVertex(ivtx)) continue;

        float mydz = trks_dz_pv(itrk,ivtx).first;

        if (fabs(mydz) < fabs(mindz)) {
          mindz = mydz;
          vtxi = ivtx;
        }

      }

      //----------------------------------------------------------------------------                                                                                                               
      // require closest PV is signal PV, dz cut, exclude tracks near hyp leptons                                                                                                                  
      //----------------------------------------------------------------------------                                                                                                               

      if ( vtxi != 0 )     continue;
    } else {

      int itrk = -1;

      if (abs(cms2.pfcands_particleId().at(ipf))!=11) {
        itrk = cms2.pfcands_trkidx().at(ipf);
        if( itrk >= (int)cms2.trks_trk_p4().size() || itrk < 0 ) continue;
        mindz=trks_dz_pv(itrk,0).first;
      }

      if (abs(cms2.pfcands_particleId().at(ipf))==11 && cms2.pfcands_pfelsidx().at(ipf)>=0) {
        itrk = cms2.els_gsftrkidx().at(cms2.pfcands_pfelsidx().at(ipf));
        if( itrk >= (int)cms2.gsftrks_p4().size() || itrk < 0 ) continue;
        mindz=gsftrks_dz_pv(itrk,0).first;
      }

    }

    //---------------------------------------                                                                                                                                     
    // passes cuts, add up isolation value                                                                                                                                        
    //---------------------------------------                                                                                                                                                    

    double coneR=0.3;
    double dr=ROOT::Math::VectorUtil::DeltaR( cms2.pfcands_p4().at(ipf) , cms2.pfcands_p4().at(thisPf) );
    if( dr > coneR ) continue; // skip pfcands outside the cone                                     

    if( cms2.pfcands_p4().at(ipf).pt()>=0.0 && fabs(mindz) <= 0.1) iso_dr03_dz01_pt00 += cms2.pfcands_p4().at(ipf).pt();

  } // end loop of the cand 

  return iso_dr03_dz01_pt00;

}
