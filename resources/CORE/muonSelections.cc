// Header
#include "muonSelections.h"
#include "eventSelections.h"

// C++ includes
#include <iostream>

// ROOT includes
#include "Math/VectorUtil.h"

// CMS2 Includes
#include "eventSelections.h"
#include "trackSelections.h"
#include "CMS2.h"
//#include "ssSelections.h"

using namespace tas;

////////////////////
// Identification //
////////////////////

bool muonId(unsigned int index, SelectionType type){

    float isovalue;
    bool  truncated = true;

    switch(type) {

        ///////////////////
        // Opposite Sign //
        ///////////////////

    case OSGeneric_v4:
        if (!muonIdNotIsolated( index, type )) return false;      
        return muonIsoValuePF(index,0,0.3) < 0.15;
        break;
    case OSGeneric_v3:
        truncated = false;
        isovalue  = 0.15;
        break;
    case OSGeneric_v3_FO:
        truncated = false;
        isovalue  = 0.4;
        break;
    case OSZ_v4:
        isovalue = 0.15;
        break;
    case OSZ_v3:
        if (!muonIdNotIsolated( index, type )) return false;      
        return muonIsoValuePF(index,0,0.3) < 0.15;
        break;
    case OSZ_v2:
        isovalue = 0.15;
        break;

        ////////////////////
        // Same Sign 2011 //
        ////////////////////
  
    case NominalSSv4:
        if (!muonIdNotIsolated(index, type)) return false;
        return (muonIsoValue(index, false) < 0.15);
        break;
    case muonSelectionFO_ssV4:
        if (!muonIdNotIsolated(index, type)) return false;
        return (muonIsoValue(index, false) < 0.40);
        break;

    case ZMet2012_v1:
        if (!muonIdNotIsolated(index, type)) return false;
        return (muonIsoValuePF2012_deltaBeta(index) < 0.15);
        break;

    case ZMet2012_detiso_v1:
        if (!muonIdNotIsolated(index, ZMet2012_v1)) return false;
        return (muonIsoValue(index,false)<0.15);
        break;

        ////////////////////
        // Same Sign 2012 //
        ////////////////////
  
    case NominalSSv5:
        if (!muonIdNotIsolated(index, type)) return false;
        return (muonIsoValuePF2012_deltaBeta(index) < 0.1);
        break;
    case muonSelectionFO_ssV5:
        if (!muonIdNotIsolated(index, type)) return false;
        return (muonIsoValuePF2012_deltaBeta(index) < 0.4);
        break;

        ///////////////
        // Higgs, WW //
        ///////////////

        // WW
    case NominalWWV0:
    case NominalWWV1:
        isovalue = 0.15;
        break;
    case muonSelectionFO_mu_wwV1:
    case muonSelectionFO_mu_ww:
        isovalue = 0.40;
        break;
    case muonSelectionFO_mu_smurf_04:
        if (!muonIdNotIsolated( index, type )) return false;
        return muonIsoValuePF(index,0,0.3) < 0.40;
        break;
    case muonSelectionFO_mu_wwV1_iso10_d0:
    case muonSelectionFO_mu_wwV1_iso10:
    case muonSelectionFO_mu_ww_iso10:
        isovalue = 1.0;
        break;

        // SMURF
    case muonSelectionFO_mu_smurf_10:
        if (!muonIdNotIsolated( index, type )) return false;
        return muonIsoValuePF(index,0,0.3) < 1.0;
        break;
    case NominalSmurfV3:
        if (!muonIdNotIsolated( index, type )) return false;
        if (cms2.mus_p4().at(index).pt()<20) 
            return muonIsoValue(index,false) < 0.1;
        else
            return muonIsoValue(index,false) < 0.15;
        break;
    case NominalSmurfV4:
        if (!muonIdNotIsolated( index, type )) return false;
        if (cms2.mus_p4().at(index).pt()>20) {
            if (TMath::Abs(cms2.mus_p4()[index].eta())<1.479) return muonIsoValuePF(index,0) < 0.22;
            else return muonIsoValuePF(index,0) < 0.20;
        } else {
            return muonIsoValuePF(index,0) < 0.11;
        }
        break;
    case NominalSmurfV5:
    case NominalSmurfV6:
        if (!muonIdNotIsolated( index, type )) return false;
        if (cms2.mus_p4().at(index).pt()>20) {
            if (TMath::Abs(cms2.mus_p4()[index].eta())<1.479) return muonIsoValuePF(index,0,0.3) < 0.13;
            else return muonIsoValuePF(index,0,0.3) < 0.09;
        } else {
            if (TMath::Abs(cms2.mus_p4()[index].eta())<1.479) return muonIsoValuePF(index,0,0.3) < 0.06;
            else return muonIsoValuePF(index,0,0.3) < 0.05;
        }
        break;

        ///////////////
        // TTV 2012  //
        ///////////////
        
        // Analysis
    case NominalTTZ_loose_v1:
        if (!muonIdNotIsolated(index, type)) return false;
        return (muonIsoValuePF2012_deltaBeta(index) < 0.15);
    case NominalTTZ_tight_v1:
        if (!muonIdNotIsolated(index, type)) return false;
        return (muonIsoValuePF2012_deltaBeta(index) < 0.10);

      // Fakes
    case NominalTTZ_looseFO_v1:
        if (!muonIdNotIsolated(index, type)) return false;
        return (muonIsoValuePF2012_deltaBeta(index) < 0.40);
    case NominalTTZ_tightFO_v1:
        if (!muonIdNotIsolated(index, type)) return false;
        return (muonIsoValuePF2012_deltaBeta(index) < 0.40);


        /////////////
        // Default //
        /////////////
    default:
        std::cout << "muonID ERROR: requested muon type is not defined. Abort." << std::endl;
        exit(1);
        return false;
    } 
    return 
        muonIdNotIsolated( index, type ) &&   // Id
        muonIsoValue(index,truncated) < isovalue;           // Isolation cut
}


bool isGoodStandardMuon( unsigned int index ){
    if ( TMath::Abs( mus_p4()[index].eta() ) > 2.4 )              return false;
    if ( mus_gfit_chi2()[index] / mus_gfit_ndof()[index] >= 50 )  return false;
    if ( ( ( mus_type()[index] ) & (1<<1) ) == 0 )                return false;
    if ( ( ( mus_type()[index] ) & (1<<2) ) == 0 )                return false;
    if ( mus_validHits()[index] < 11 )                            return false;
    if ( mus_gfit_validSTAHits()[index] == 0)                     return false;
    return true;
}

////////////////////
// Identification //
////////////////////

bool muonIdNotIsolated(unsigned int index, SelectionType type) {

    if ( cms2.mus_p4()[index].pt() < 5.0) {
        // std::cout << "muonID ERROR: requested muon is too low pt,  Abort." << std::endl;
        return false;
    }

    int vtxidx = firstGoodVertex();
    int trkidx = cms2.mus_trkidx().at(index);
    
    // Muon Selections that are standard for Analysis & Fake Selections
    bool standardMuon = true;
    if ( TMath::Abs( mus_p4()[index].eta() ) > 2.4 )              standardMuon = false;
    if ( mus_gfit_chi2()[index] / mus_gfit_ndof()[index] >= 50 )  standardMuon = false;
    if ( ( ( mus_type()[index] ) & (1<<1) ) == 0 )                standardMuon = false;
    if ( ( ( mus_type()[index] ) & (1<<2) ) == 0 )                standardMuon = false;
    if ( mus_validHits()[index] < 11 )                            standardMuon = false;
    if ( mus_gfit_validSTAHits()[index] == 0)                     standardMuon = false;
    if ( mus_ptErr()[index] / mus_p4()[index].pt() > 0.1 )        standardMuon = false;


    //
    switch (type) {

    case NominalWWV0:
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits  
        if (cms2.mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mud0PV(index)) >= 0.02)              return false; // d0 from pvtx
        return true;
        break;

    case muonSelectionFO_mu_wwV1:
    case muonSelectionFO_mu_wwV1_iso10:
    case NominalWWV1:
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits  
        if (cms2.mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mud0PV_wwV1(index)) >= 0.02)         return false; // d0 from pvtx
        if (TMath::Abs(mudzPV_wwV1(index)) >= 1.0)          return false; // dz from pvtx
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt()>0.1) return false;
        if (cms2.trks_valid_pixelhits().at(cms2.mus_trkidx().at(index))==0) return false;
        if (cms2.mus_nmatches().at(index)<2) return false;
        return true;
        break;

    case muonSelectionFO_mu_wwV1_iso10_d0: // same as muonSelectionFO_mu_wwV1_iso10 but with looser d0 cut
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits  
        if (cms2.mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mud0PV_wwV1(index)) >= 0.2)         return false; // d0 from pvtx
        if (TMath::Abs(mudzPV_wwV1(index)) >= 1.0)          return false; // dz from pvtx
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt()>0.1) return false;
        if (cms2.trks_valid_pixelhits().at(cms2.mus_trkidx().at(index))==0) return false;
        if (cms2.mus_nmatches().at(index)<2) return false;
        return true;
        break;

    case muonSelectionFO_mu_ww:
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits  
        if (cms2.mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mud0PV(index)) >= 0.02)              return false; // d0 from pvtx
        return true;

    case muonSelectionFO_mu_ww_iso10:
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits  
        if (cms2.mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mud0PV(index)) >= 0.02)              return false; // d0 from pvtx
        return true;

    case OSGeneric_v4:
        //baseline selector for 2011 OS analysis
        if( !standardMuon )                                                      return false; // |eta| < 2.4, chisq/ndof < 50, tracker & global muon, 11 or more TRK hits, glb fit mu hits, dpt/pt < 0.1
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; // glb fit chisq
        if (TMath::Abs(mud0PV_smurfV3(index)) > 0.02)                            return false; // d0(PV) < 0.02 cm
        if (TMath::Abs(mudzPV_smurfV3(index)) > 1  )                             return false; // dz(PV) < 1 cm
        return true;

    case OSGeneric_v3:
        //baseline selector for 2011 OS analysis
        if( !standardMuon )                                                      return false; // |eta| < 2.4, chisq/ndof < 50, tracker & global muon, 11 or more TRK hits, glb fit mu hits, dpt/pt < 0.1
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; // glb fit chisq
        if (TMath::Abs(mud0PV_smurfV3(index)) > 0.02)                            return false; // d0(PV) < 0.02 cm
        if (TMath::Abs(mudzPV_smurfV3(index)) > 1  )                             return false; // dz(PV) < 1 cm
        return true;

    case OSGeneric_v3_FO:
	    // Fakes for 2011: reliso < 0.4, d0 < 0.2, chisq/ndof < 50
        if( !standardMuon )                                                      return false; // |eta| < 2.4, chisq/ndof < 50, tracker & global muon, 11 or more TRK hits, glb fit mu hits, dpt/pt < 0.1
        if (TMath::Abs(mud0PV_smurfV3(index)) > 0.2)                             return false; // d0(PV) < 0.2 cm
        if (TMath::Abs(mudzPV_smurfV3(index)) > 1  )                             return false; // dz(PV) < 1 cm
        return true;

    case OSZ_v4:
        // baseline selector for 2011 Z+MET analysis
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)                       return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; // glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)                         return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)                                 return false; // # of tracker hits
        if (cms2.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mud0PV_smurfV3(index)) > 0.02)                            return false; // d0(PV) < 0.02 cm
        if (TMath::Abs(mudzPV_smurfV3(index)) > 1  )                             return false; // dz(PV) < 1 cm
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt()>0.1)         return false; // dpt/pt 
        if (cms2.mus_iso_ecalvetoDep().at(index) > 4)                            return false; // ECalE < 4 
        if (cms2.mus_iso_hcalvetoDep().at(index) > 6)                            return false; // HCalE < 6 
        if (!isPFMuon(index,true,1.0))                                           return false; // require muon is pfmuon with same pt
        return true;
        break;

    case OSZ_v3:
        // baseline selector for 2011 Z+MET analysis
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)                       return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; // glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)                         return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)                                 return false; // # of tracker hits
        if (cms2.mus_iso_ecalvetoDep().at(index) > 4)                            return false; // ECalE < 4 
        if (cms2.mus_iso_hcalvetoDep().at(index) > 6)                            return false; // HCalE < 6 
        if (cms2.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(cms2.mus_d0corr().at(index)) > 0.02)                      return false; // d0 from beamspot
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt()>0.1)         return false; // dpt/pt 
        if (cms2.mus_iso_ecalvetoDep().at(index) > 4)                            return false; // ECalE < 4 
        if (cms2.mus_iso_hcalvetoDep().at(index) > 6)                            return false; // HCalE < 6 
        if (!isPFMuon(index,true,1.0))                                           return false; // require muon is pfmuon with same pt
        return true;
        break;

    case OSZ_v2:
        // baseline selector for 2011 Z+MET analysis
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)                       return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; // glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)                         return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)                                 return false; // # of tracker hits
        if (cms2.mus_iso_ecalvetoDep().at(index) > 4)                            return false; // ECalE < 4 
        if (cms2.mus_iso_hcalvetoDep().at(index) > 6)                            return false; // HCalE < 6 
        if (cms2.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(cms2.mus_d0corr().at(index)) > 0.02)                      return false; // d0 from beamspot
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt()>0.1)         return false; // dpt/pt 
        if (cms2.mus_iso_ecalvetoDep().at(index) > 4)                            return false; // ECalE < 4 
        if (cms2.mus_iso_hcalvetoDep().at(index) > 6)                            return false; // HCalE < 6 
        if (!isPFMuon(index,true,1.0))                                           return false; // require muon is pfmuon with same pt
        return true;
        break;

    case NominalSmurfV3:
    case NominalSmurfV4:
    case NominalSmurfV5:
        if (type == NominalSmurfV3 || type == NominalSmurfV4 || type == NominalSmurfV5){
            if (cms2.mus_p4().at(index).pt()<20){
                if (TMath::Abs(mud0PV_smurfV3(index)) >= 0.01)    return false; // d0 from pvtx
            } else {
                if (TMath::Abs(mud0PV_smurfV3(index)) >= 0.02)    return false; // d0 from pvtx
            }
        } else {
            if (TMath::Abs(mud0PV_smurfV3(index)) >= 0.2)    return false; // d0 from pvtx
        }
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits  
        if (cms2.mus_gfit_validSTAHits().at(index)==0 )     return false; // Glb fit must have hits in mu chambers
        if (TMath::Abs(mudzPV_smurfV3(index)) >= 0.1)       return false; // dz from pvtx
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt()>0.1) return false;
        if (cms2.trks_valid_pixelhits().at(cms2.mus_trkidx().at(index))==0) return false;
        if (cms2.mus_nmatches().at(index)<2) return false;
        return true;
        break;

        //baseline selector for 2011 SS analysis
    case NominalSSv4:
        if (fabs(cms2.mus_p4().at(index).eta()) > 2.4)                           return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; // glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)                         return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)                                 return false; // # of tracker hits
        if (cms2.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt() > 0.1)       return false; // dpt/pt < 0.1
        if (cms2.mus_iso_ecalvetoDep().at(index) > 4)                            return false; // ECalE < 4 
        if (cms2.mus_iso_hcalvetoDep().at(index) > 6)                            return false; // HCalE < 6
        // cut on d0, dz using first good vertex
        // if there isn't a good vertex, use the beamSpot
        if (vtxidx < 0 || cms2.mus_trkidx().at(index) < 0) {
            if (fabs(cms2.mus_d0corr().at(index)) > 0.02)
                return false;
        }
        else if (vtxidx >= 0) {
            if (fabs(trks_d0_pv(cms2.mus_trkidx().at(index), vtxidx).first) > 0.02)
                return false;            
        }
        return true;
        break;
        //baseline selector for 2011 SS analysis
    case muonSelectionFO_ssV4:
        if (fabs(cms2.mus_p4().at(index).eta()) > 2.4)                           return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 50) return false; // glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (((cms2.mus_type().at(index)) & (1<<2)) == 0)                         return false; // tracker muon
        if (cms2.mus_validHits().at(index) < 11)                                 return false; // # of tracker hits
        if (cms2.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt() > 0.1)       return false; // dpt/pt < 0.1
        // cut on d0, dz using first good vertex
        // if there isn't a good vertex, use the beamSpot
        if (vtxidx < 0 || cms2.mus_trkidx().at(index) < 0) {
            if (fabs(cms2.mus_d0corr().at(index)) > 0.2)
                return false;
        }
        else if (vtxidx >= 0) {
            if (fabs(trks_d0_pv(cms2.mus_trkidx().at(index), vtxidx).first) > 0.2)
                return false;            
        }
        return true;
        break;

        // muon POG tight muon requirements, with d0 cut tightened to 0.02 cm
	// see: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId

    case ZMet2012_v1:
        if (fabs(cms2.mus_p4().at(index).eta()) > 2.4)                           return false; // eta cut
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (cms2.mus_pid_PFMuon().at(index) == 0)                                return false; // pf muon
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; // glb fit chisq
        if (cms2.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
	if (cms2.mus_numberOfMatchedStations().at(index) < 2)                    return false; // require muon segements in at least two muon stations

        // cut on d0, dz using first good vertex
        // if there isn't a good vertex, use the beamSpot
        if (trkidx < 0)                                                          return false; // require a matching track
        if (vtxidx < 0 || trkidx < 0) {
	  cout << __FILE__ << " " << __LINE__ << endl;
	  cout << "WARNING: didn't find any good vertices, should never get here" << endl;

            if (fabs(cms2.mus_d0corr().at(index)) > 0.02) return false;
            if (fabs(cms2.mus_z0corr().at(index)) > 0.5)  return false;
        }
        else if (vtxidx >= 0) {
            if (fabs(trks_d0_pv(trkidx, vtxidx).first) > 0.02) return false;
            if (fabs(trks_dz_pv(trkidx, vtxidx).first) > 0.5)  return false;
        }
        else return false;     

        if (cms2.trks_valid_pixelhits().at(trkidx) == 0)                         return false; // require at least 1 valid pixel hit
        if (cms2.trks_nlayers().at(trkidx) < 6)                                  return false; // require at least 6 tracker layers with hits

        return true;
        break;

        //baseline selector for 2012 SS analysis
    case NominalSSv5:
        if (fabs(cms2.mus_p4().at(index).eta()) > 2.4)                           return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; // glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (cms2.mus_pid_PFMuon().at(index) == 0)                                return false; // pf muon
        if (cms2.mus_numberOfMatchedStations().at(index) < 2)                    return false; // require muon segements in at least two muon stations

        if (trkidx < 0)                                                          return false; // require a matching track
        if (cms2.trks_nlayers().at(trkidx) < 6)                                  return false; // require at least 6 tracker layers with hits
        if (cms2.trks_valid_pixelhits().at(trkidx) == 0)                         return false; // require at least 1 valid pixel hit

        if (cms2.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        if (cms2.mus_iso_ecalvetoDep().at(index) > 4)                            return false; // ECalE < 4 
        if (cms2.mus_iso_hcalvetoDep().at(index) > 6)                            return false; // HCalE < 6
        // cut on d0, dz using first good vertex
        // if there isn't a good vertex, use the beamSpot
        if (vtxidx < 0 || trkidx < 0) {
            if (fabs(cms2.mus_d0corr().at(index)) > 0.02)
                return false;
            if (fabs(cms2.mus_z0corr().at(index)) > 0.1)
                return false;
        }
        else if (vtxidx >= 0) {
            if (fabs(trks_d0_pv(trkidx, vtxidx).first) > 0.02)
                return false;
            if (fabs(trks_dz_pv(trkidx, vtxidx).first) > 0.1)
                return false;
        }
        else return false;     
        return true;
        break;
        //baseline FO selector for 2012 SS analysis
    case muonSelectionFO_ssV5:
        if (fabs(cms2.mus_p4().at(index).eta()) > 2.4)                           return false; // eta cut
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 50) return false; // glb fit chisq
        if (((cms2.mus_type().at(index)) & (1<<1)) == 0)                         return false; // global muon
        if (cms2.mus_pid_PFMuon().at(index) == 0)                                return false; // pf muon
        if (cms2.mus_numberOfMatchedStations().at(index) < 2)                    return false; // require muon segements in at least two muon stations

        if (trkidx < 0)                                                          return false; // require a matching track
        if (cms2.trks_nlayers().at(trkidx) < 6)                                  return false; // require at least 6 tracker layers with hits
        if (cms2.trks_valid_pixelhits().at(trkidx) == 0)                         return false; // require at least 1 valid pixel hit

        if (cms2.mus_gfit_validSTAHits().at(index) == 0)                         return false; // Glb fit must have hits in mu chambers
        // cut on d0, dz using first good vertex
        // if there isn't a good vertex, use the beamSpot
        if (vtxidx < 0 || trkidx < 0) {
            if (fabs(cms2.mus_d0corr().at(index)) > 0.2)
                return false;
            if (fabs(cms2.mus_z0corr().at(index)) > 0.1)
                return false;
        }
        else if (vtxidx >= 0) {
            if (fabs(trks_d0_pv(cms2.mus_trkidx().at(index), vtxidx).first) > 0.2)
                return false;            
            if (fabs(trks_dz_pv(trkidx, vtxidx).first) > 0.1)
                return false;
        }
        else return false;
        return true;
        break;

    case muonSelectionFO_mu_smurf_04:
    case muonSelectionFO_mu_smurf_10:
    case NominalSmurfV6:
    {
        if (type == NominalSmurfV6){
            if (cms2.mus_p4().at(index).pt()<20){
                if (TMath::Abs(mud0PV_smurfV3(index)) >= 0.01)    return false; // d0 from pvtx
            } else {
                if (TMath::Abs(mud0PV_smurfV3(index)) >= 0.02)    return false; // d0 from pvtx
            }
        } else {
            if (TMath::Abs(mud0PV_smurfV3(index)) >= 0.2)    return false; // d0 from pvtx
        }
        if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4)  return false; // eta cut
        if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits  
        if (TMath::Abs(mudzPV_smurfV3(index)) >= 0.1)       return false; // dz from pvtx
        if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt()>0.1) return false;
        if (cms2.trks_valid_pixelhits().at(cms2.mus_trkidx().at(index))==0) return false;
        bool goodMuonGlobalMuon = false;
        if (((cms2.mus_type().at(index)) & (1<<1)) != 0) { // global muon
            goodMuonGlobalMuon = true;
            if (cms2.mus_nmatches().at(index)<2) goodMuonGlobalMuon = false;
            if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) goodMuonGlobalMuon = false; //glb fit chisq
            if (cms2.mus_gfit_validSTAHits().at(index)==0 ) goodMuonGlobalMuon = false; // Glb fit must have hits in mu chambers
        } 
        bool goodMuonTrackerMuon = false;
        if (((cms2.mus_type().at(index)) & (1<<2)) != 0) { // tracker muon
            goodMuonTrackerMuon = true;
            if (cms2.mus_pid_TMLastStationTight().at(index) == 0 ) goodMuonTrackerMuon = false; // last station tight
        }
        return goodMuonGlobalMuon || goodMuonTrackerMuon;
        break;
    }
    case NominalTTZ_loose_v1:
        if (!passes_muid_wp2012(index, mu2012_tightness::TIGHT)) return false;
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false;
        if (trkidx < 0 || vtxidx < 0) return false;
        if (fabs(trks_d0_pv(trkidx, vtxidx).first) > 0.02) return false;
		if (fabs(trks_dz_pv(trkidx, vtxidx).first) > 0.10) return false;
        return true;
        break;
    case NominalTTZ_tight_v1:
        if (!passes_muid_wp2012(index, mu2012_tightness::TIGHT)) return false;
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false;
        if (trkidx < 0 || vtxidx < 0) return false;
        if (fabs(trks_d0_pv(trkidx, vtxidx).first) > 0.02) return false;
        if (fabs(trks_dz_pv(trkidx, vtxidx).first) > 0.10) return false;
        if (cms2.mus_iso_ecalvetoDep().at(index) > 4) return false; // ECalE < 4 
        if (cms2.mus_iso_hcalvetoDep().at(index) > 6) return false; // HCalE < 6
        return true;
        break;
    case NominalTTZ_looseFO_v1:
        if (!passes_muid_wp2012(index, mu2012_tightness::TIGHT)) return false;
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 50) return false;
        if (trkidx < 0 || vtxidx < 0) return false;
        if (fabs(trks_d0_pv(trkidx, vtxidx).first) > 0.20) return false;
		if (fabs(trks_dz_pv(trkidx, vtxidx).first) > 0.10) return false;
        return true;
        break;
    case NominalTTZ_tightFO_v1:
        if (!passes_muid_wp2012(index, mu2012_tightness::TIGHT)) return false;
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 50) return false;
        if (trkidx < 0 || vtxidx < 0) return false;
        if (fabs(trks_d0_pv(trkidx, vtxidx).first) > 0.20) return false;
        if (fabs(trks_dz_pv(trkidx, vtxidx).first) > 0.10) return false;
        return true;
        break;        

    default:
        std::cout << "muonID ERROR: requested muon type is not defined. Abort." << std::endl;
        return false;
    }
}




////////////////////////////
// Isolation Calculations //
////////////////////////////

double muonIsoValue(unsigned int index, bool truncated ){
    return ( muonIsoValue_TRK( index, truncated ) + muonIsoValue_ECAL( index, truncated ) + muonIsoValue_HCAL( index, truncated ) );
}
double muonIsoValue_FastJet(unsigned int index, bool truncated ){
    return ( muonIsoValue_TRK( index, truncated ) + TMath::Max( muonIsoValue_ECAL( index, truncated ) + muonIsoValue_HCAL( index, truncated ) - mu_fastjet_rel_offset(index,truncated) , 0.0 ) );
}
double mu_fastjet_rel_offset(unsigned int index, bool truncated ){
    double pt        = cms2.mus_p4().at(index).pt();
    if(truncated) pt = max( pt, 20.0 );
    double offset = TMath::Pi() * pow( 0.3 , 2 ) * cms2.evt_rho();
    return offset / pt;
}
double muonIsoValue_TRK(unsigned int index, bool truncated ){
    double pt        = cms2.mus_p4().at(index).pt();
    if(truncated) pt = max( pt, 20.0 );
    return cms2.mus_iso03_sumPt().at(index) / pt;
}
double muonIsoValue_ECAL(unsigned int index, bool truncated ){
    double pt  = cms2.mus_p4().at(index).pt();
    if(truncated) pt = max( pt, 20.0 );
    return cms2.mus_iso03_emEt().at(index) / pt;
}
double muonIsoValue_HCAL(unsigned int index, bool truncated){
    double pt  = cms2.mus_p4().at(index).pt();
    if(truncated) pt = max( pt, 20.0 );
    return cms2.mus_iso03_hadEt().at(index) / pt;
}
double muonCorIsoValue (unsigned int index, bool truncated) {
    double ntiso  = muonIsoValue(index, truncated);
    double pt     = cms2.mus_p4().at(index).pt();
    int nvtxs     = numberOfGoodVertices();
    double coriso = ntiso - ((TMath::Log(pt)*nvtxs)/(30*pt));
    return coriso;
}

#ifdef PFISOFROMNTUPLE
double muonIsoValuePF( unsigned int imu, unsigned int ivtx, float coner, float minptn, float dzcut, int filterId){
    if (fabs(coner-0.3)<0.0001) {
        if (cms2.mus_iso03_pf().at(imu)<-99.) return 9999.;
        return cms2.mus_iso03_pf().at(imu)/cms2.mus_p4().at(imu).pt();
    } else if (fabs(coner-0.4)<0.0001) {
        if (cms2.mus_iso04_pf().at(imu)<-99.) return 9999.;
        return cms2.mus_iso04_pf().at(imu)/cms2.mus_p4().at(imu).pt();
    } else {
        cout << "muonIsoValuePF: CONE SIZE NOT SUPPORTED" << endl;
        return 9999.;
    }
}
#else
double muonIsoValuePF( unsigned int imu, unsigned int ivtx, float coner, float minptn, float dzcut, int filterId){
    float pfciso = 0;
    float pfniso = 0;
    int mutkid = cms2.mus_trkidx().at(imu);
    float mudz = mutkid>=0 ? trks_dz_pv(mutkid,ivtx).first : cms2.mus_sta_z0corr().at(imu);
    for (unsigned int ipf=0; ipf<cms2.pfcands_p4().size(); ++ipf){
        float dR = ROOT::Math::VectorUtil::DeltaR( pfcands_p4().at(ipf), mus_p4().at(imu) );
        if (dR>coner) continue;
        float pfpt = cms2.pfcands_p4().at(ipf).pt();
        int pfid = abs(cms2.pfcands_particleId().at(ipf));
        if (filterId!=0 && filterId!=pfid) continue;
        if (cms2.pfcands_charge().at(ipf)==0) {
            //neutrals
            if (pfpt>minptn) pfniso+=pfpt;
        } else {
            //charged
            //avoid double counting of muon itself
            int pftkid = cms2.pfcands_trkidx().at(ipf);
            if (mutkid>=0 && pftkid>=0 && mutkid==pftkid) continue;
            //first check electrons with gsf track
            if (abs(cms2.pfcands_particleId().at(ipf))==11 && cms2.pfcands_pfelsidx().at(ipf)>=0 && cms2.pfels_elsidx().at(cms2.pfcands_pfelsidx().at(ipf))>=0) {
                int gsfid = cms2.els_gsftrkidx().at(cms2.pfels_elsidx().at(cms2.pfcands_pfelsidx().at(ipf))); 
                if (gsfid>=0) { 
                    if(fabs(gsftrks_dz_pv( gsfid,ivtx ).first - mudz )<dzcut) {//dz cut
                        pfciso+=pfpt;
                    }   
                    continue;//and avoid double counting
                }
            }
            //then check anything that has a ctf track
            if (cms2.pfcands_trkidx().at(ipf)>=0) {//charged (with a ctf track)
                if(fabs( trks_dz_pv(cms2.pfcands_trkidx().at(ipf),ivtx).first - mudz )<dzcut) {//dz cut
                    pfciso+=pfpt;
                }
            } 
        }
    } 
    return (pfciso+pfniso)/cms2.mus_p4().at(imu).pt();
}
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Remove cosmics by looking for back-to-back muon-track pairs ( http://indico.cern.ch/contributionDisplay.py?contribId=2&confId=86834 ) //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool isCosmics(unsigned int index){
    for (int itrk=0; itrk < int(cms2.trks_trk_p4().size()); ++itrk) {
        const LorentzVector& mu_p4  = cms2.mus_trk_p4().at(index);
        const LorentzVector& trk_p4 = cms2.trks_trk_p4().at(itrk);
        double sprod = mu_p4.px()*trk_p4.px()+mu_p4.py()*trk_p4.py()+mu_p4.pz()*trk_p4.pz();
        if ( acos( -(sprod/trk_p4.P()/mu_p4.P()) ) < 0.01 &&
             fabs(trk_p4.pt()-mu_p4.pt())/mu_p4.pt() < 0.05 )
            return true;
    }
    return false;
}

/////////////////////////////
// Muon d0 corrected by PV //
////////////////////////////

double mud0PV(unsigned int index){
    if ( cms2.vtxs_sumpt().empty() ) return 9999.;
    unsigned int iMax = 0;
    double sumPtMax = cms2.vtxs_sumpt().at(0);
    for ( unsigned int i = iMax+1; i < cms2.vtxs_sumpt().size(); ++i )
        if ( cms2.vtxs_sumpt().at(i) > sumPtMax ){
            iMax = i;
            sumPtMax = cms2.vtxs_sumpt().at(i);
        }
    double dxyPV = cms2.mus_d0()[index]-
        cms2.vtxs_position()[iMax].x()*sin(cms2.mus_trk_p4()[index].phi())+
        cms2.vtxs_position()[iMax].y()*cos(cms2.mus_trk_p4()[index].phi());
    return dxyPV;
}

double mud0PV_wwV1(unsigned int index){
    if ( cms2.vtxs_sumpt().empty() ) return 9999.;
    double sumPtMax = -1;
    int iMax = -1;
    for ( unsigned int i = 0; i < cms2.vtxs_sumpt().size(); ++i ){
        // if (!isGoodVertex(i)) continue;
        // Copied from eventSelections.cc 
        if (cms2.vtxs_isFake()[i]) continue;
        if (cms2.vtxs_ndof()[i] < 4.) continue;
        if (cms2.vtxs_position()[i].Rho() > 2.0) continue;
        if (fabs(cms2.vtxs_position()[i].Z()) > 24.0) continue;
        if ( cms2.vtxs_sumpt().at(i) > sumPtMax ){
            iMax = i;
            sumPtMax = cms2.vtxs_sumpt().at(i);
        }
    }
    if (iMax<0) return 9999.;
    double dxyPV = cms2.mus_d0()[index]-
        cms2.vtxs_position()[iMax].x()*sin(cms2.mus_trk_p4()[index].phi())+
        cms2.vtxs_position()[iMax].y()*cos(cms2.mus_trk_p4()[index].phi());
    return dxyPV;
}

double mud0PV_smurfV3(unsigned int index){
    int vtxIndex = 0;
    double dxyPV = cms2.mus_d0()[index]-
        cms2.vtxs_position()[vtxIndex].x()*sin(cms2.mus_trk_p4()[index].phi())+
        cms2.vtxs_position()[vtxIndex].y()*cos(cms2.mus_trk_p4()[index].phi());
    return dxyPV;
}

double dzPV_mu(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv){
    return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt();
}

double mudzPV_smurfV3(unsigned int index){
    int vtxIndex = 0;
    double dzpv = dzPV_mu(cms2.mus_vertex_p4()[index], cms2.mus_trk_p4()[index], cms2.vtxs_position()[vtxIndex]);
    return dzpv;
}

double mudzPV_wwV1(unsigned int index){
    if ( cms2.vtxs_sumpt().empty() ) return 9999.;
    double sumPtMax = -1;
    int iMax = -1;
    for ( unsigned int i = 0; i < cms2.vtxs_sumpt().size(); ++i ){
        // if (!isGoodVertex(i)) continue;
        // Copied from eventSelections.cc 
        if (cms2.vtxs_isFake()[i]) continue;
        if (cms2.vtxs_ndof()[i] < 4.) continue;
        if (cms2.vtxs_position()[i].Rho() > 2.0) continue;
        if (fabs(cms2.vtxs_position()[i].Z()) > 24.0) continue;
        if ( cms2.vtxs_sumpt().at(i) > sumPtMax ){
            iMax = i;
            sumPtMax = cms2.vtxs_sumpt().at(i);
        }
    }
    if (iMax<0) return 9999.;
    // double dzpv = cms2.mus_z0corr()[index]-cms2.vtxs_position()[iMax].z();
    const LorentzVector& vtx = cms2.mus_vertex_p4()[index];
    const LorentzVector& p4 = cms2.mus_trk_p4()[index];
    const LorentzVector& pv = cms2.vtxs_position()[iMax];
    return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt(); 
    /* directly from NtupleMacros/WW/doAnalysis.cc
       double dzpv = dzPV(cms2.mus_vertex_p4()[index], cms2.mus_trk_p4()[index], cms2.vtxs_position()[iMax]);
       double dzPV(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv){
       return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt();
       }*/
}

bool isPFMuon( int index , bool requireSamePt , float dpt_max ){

    int ipf = cms2.mus_pfmusidx().at( index );

    //--------------------------
    // require matched pfmuon
    //--------------------------

    if( ipf >= int(cms2.pfmus_p4().size()) || ipf < 0 ) return false;

    //----------------------------------------------------
    // require PFMuon pt = reco muon pt (within dpt_max)
    //----------------------------------------------------

    if( requireSamePt ){

        float pt_pf = cms2.pfmus_p4().at(ipf).pt();
        float pt    = cms2.mus_p4().at(index).pt();

        if( fabs( pt_pf - pt ) > dpt_max ) return false;

    }

    return true;

}

void muonIsoValuePF2012 (float &pfiso_ch, float &pfiso_em, float &pfiso_nh, const float R, const unsigned int imu, const int ivtx, float neutral_et_threshold)
{

    // isolation sums
    pfiso_ch = 0.0;
    pfiso_em = 0.0; 
    pfiso_nh = 0.0;
       
    // loop on pfcandidates
    for (unsigned int ipf = 0; ipf < cms2.pfcands_p4().size(); ++ipf) {
            
        // skip electrons and muons
        const int particleId = abs(cms2.pfcands_particleId()[ipf]);
        if (particleId == 11)    continue;
        if (particleId == 13)    continue;
    
        // deltaR between electron and cadidate
        const float dR = ROOT::Math::VectorUtil::DeltaR(cms2.pfcands_p4()[ipf], cms2.mus_p4()[imu]);
        if (dR > R)              continue;

        // charged hadrons closest vertex
        // should be the primary vertex
        if (particleId == 211 || particleId == 321 || particleId == 2212 || particleId == 999211) {
            if (cms2.pfcands_vtxidx().at(ipf) != ivtx) continue;
            if (dR < 0.0001)
                continue;
        }
        if (particleId == 22 || particleId == 130 || particleId == 111 || particleId == 310 || particleId == 2112) {
            if (cms2.pfcands_p4().at(ipf).pt() < neutral_et_threshold)
                continue;
            if (dR < 0.01)
                continue;
        }

        // add to isolation sum
        if (particleId == 211 || particleId == 321 || particleId == 2212 || particleId == 999211)      pfiso_ch += cms2.pfcands_p4()[ipf].pt();
        if (particleId == 22)                                                                          pfiso_em += cms2.pfcands_p4()[ipf].pt();
        if (particleId == 130 || particleId == 111 || particleId == 310 || particleId == 2112)         pfiso_nh += cms2.pfcands_p4()[ipf].pt();
    }
}

float muonIsoValuePF2012_FastJetEffArea(int index, float conesize, float effective_area, int ivtx)
{
    float pt     = cms2.mus_p4()[index].pt();

    // pf iso
    // calculate from the ntuple for now...
    float pfiso_ch = 0.0;
    float pfiso_em = 0.0;
    float pfiso_nh = 0.0;
    muonIsoValuePF2012(pfiso_ch, pfiso_em, pfiso_nh, conesize, index, ivtx);

    // rho
    float rhoPrime = std::max(cms2.evt_rho(), float(0.0));
    float pfiso_n = std::max(pfiso_em + pfiso_nh - rhoPrime * effective_area, float(0.0));
    float pfiso = (pfiso_ch + pfiso_n) / pt;   

    return pfiso;    
}

float muonRadialIsolation (unsigned int imu, float &chiso, float &nhiso, float &emiso, float neutral_et_threshold, float cone_size, bool verbose)
{
    float radial_iso = 0.;
    chiso = 0.;
    nhiso = 0.;
    emiso = 0.;

    int ivtx = firstGoodVertex();

    LorentzVector p4 = (cms2.mus_trk_p4().at(imu).pt() > 0.01) ? cms2.mus_trk_p4().at(imu) : cms2.mus_sta_p4().at(imu);
    if (p4.pt() < 0.01)
        return -9999.;

    for (unsigned int ipf = 0; ipf < cms2.pfcands_p4().size(); ipf++) {

        // skip electrons and muons
        const int particleId = abs(cms2.pfcands_particleId().at(ipf));
        if (particleId == 11) {
            if (verbose)
                std::cout << "Skipping electron with id, pt, eta = " << cms2.pfcands_particleId().at(ipf) << ", " << cms2.pfcands_p4().at(ipf).pt() << ", " << cms2.pfcands_p4().at(ipf).eta() << std::endl;
            continue;
        }
        if (particleId == 13) {
            if (verbose)
                std::cout << "Skipping muon with id, pt, eta = " << cms2.pfcands_particleId().at(ipf) << ", " << cms2.pfcands_p4().at(ipf).pt() << ", " << cms2.pfcands_p4().at(ipf).eta() << std::endl;
            continue;
        }

        // in the event that the muon is not a PF muon, need to remove any other PF cand reconstructed using the same track as the muon
        if (!cms2.mus_pid_PFMuon().at(imu) && cms2.mus_trkidx().at(imu) >= 0 && cms2.mus_trkidx().at(imu) == cms2.pfcands_trkidx().at(ipf)) {
            if (verbose)
                std::cout << "Skipping PF cand with same track as muon with id, pt, eta = " << cms2.pfcands_particleId().at(ipf) << ", " << cms2.pfcands_p4().at(ipf).pt() << ", " << cms2.pfcands_p4().at(ipf).eta() << std::endl;
            continue;
        }

        const float dr = ROOT::Math::VectorUtil::DeltaR(cms2.pfcands_p4().at(ipf), cms2.mus_p4().at(imu));
        if (dr > cone_size) {
            if (verbose)
                std::cout << "Skipping PF candidate outside of cone with id, pt, eta = " << cms2.pfcands_particleId().at(ipf) << ", " << cms2.pfcands_p4().at(ipf).pt() << ", " << cms2.pfcands_p4().at(ipf).eta() << std::endl;            
            continue;
        }
        if (dr < 0.01) {
            if (verbose)
                std::cout << "Skipping PF candidate in veto cone with id, pt, eta = " << cms2.pfcands_particleId().at(ipf) << ", " << cms2.pfcands_p4().at(ipf).pt() << ", " << cms2.pfcands_p4().at(ipf).eta() << std::endl;            
            continue;
        }

        // deal with charged
        if (cms2.pfcands_charge().at(ipf) != 0) {
            if (cms2.pfcands_vtxidx().at(ipf) != ivtx) {
                if (verbose)
                    std::cout << "Skipping PF candidate from other vertex  with id, pt, eta, ivtx = " << cms2.pfcands_particleId().at(ipf) << ", " << cms2.pfcands_p4().at(ipf).pt() << ", " 
                              << cms2.pfcands_p4().at(ipf).eta() << ", " << cms2.pfcands_vtxidx().at(ipf) << std::endl;
                continue;
            }
            radial_iso += cms2.pfcands_p4().at(ipf).pt() * (1 - 3*dr) / cms2.mus_p4().at(imu).pt();
            chiso += cms2.pfcands_p4().at(ipf).pt() * (1 - 3*dr) / cms2.mus_p4().at(imu).pt();
            if (verbose)
                std::cout << "Summing CH with id, pt, eta = " << cms2.pfcands_particleId().at(ipf) << ", " << cms2.pfcands_p4().at(ipf).pt() << ", " << cms2.pfcands_p4().at(ipf).eta() << std::endl;            
        }
        else if (cms2.pfcands_p4().at(ipf).pt() > neutral_et_threshold) {
            radial_iso += cms2.pfcands_p4().at(ipf).pt() * (1 - 3*dr) / cms2.mus_p4().at(imu).pt();
            if (particleId == 22) {
                emiso += cms2.pfcands_p4().at(ipf).pt() * (1 - 3*dr) / cms2.mus_p4().at(imu).pt();
                if (verbose)
                    std::cout << "Summing EM with id, pt, eta = " << cms2.pfcands_particleId().at(ipf) << ", " << cms2.pfcands_p4().at(ipf).pt() << ", " << cms2.pfcands_p4().at(ipf).eta() << std::endl;            
            }
            else {
                nhiso += cms2.pfcands_p4().at(ipf).pt() * (1 - 3*dr) / cms2.mus_p4().at(imu).pt();
                if (verbose)
                    std::cout << "Summing NH with id, pt, eta = " << cms2.pfcands_particleId().at(ipf) << ", " << cms2.pfcands_p4().at(ipf).pt() << ", " << cms2.pfcands_p4().at(ipf).eta() << std::endl;            
            }
        }
    } // loop over pfcands

    return radial_iso;
}

float muonIsoValuePF2012_deltaBeta(unsigned int imu)
{
    const float chiso = cms2.mus_isoR03_pf_ChargedHadronPt().at(imu);
    const float nhiso = cms2.mus_isoR03_pf_NeutralHadronEt().at(imu);
    const float emiso = cms2.mus_isoR03_pf_PhotonEt().at(imu);
    const float deltaBeta = cms2.mus_isoR03_pf_PUPt().at(imu);
    const float pt = cms2.mus_p4().at(imu).pt();

    const float absiso = chiso + max(0.0, nhiso + emiso - 0.5 * deltaBeta);
    return (absiso / pt);
}

bool passes_muid_wp2012 (const unsigned int index, const mu2012_tightness::value_type tightness)
{
    const bool is_global  = ((cms2.mus_type().at(index) & (1<<1)) != 0);
    const bool is_tracker = ((cms2.mus_type().at(index) & (1<<2)) != 0);
    const bool is_pfmu    = ((cms2.mus_type().at(index) & (1<<5)) != 0);

    const int vtxidx = firstGoodVertex();

    switch (tightness) {
        
    case mu2012_tightness::LOOSE: {
        if (!is_pfmu) return false;
        if (!is_global && !is_tracker) return false;      
        return true;
		}
		break;

    case mu2012_tightness::TIGHT: {
        if (!is_global) return false;
        if (!is_pfmu) return false;
        if (cms2.mus_gfit_validSTAHits().at(index) < 1) return false;
        if (cms2.mus_numberOfMatchedStations().at(index) < 2) return false;

        const int ctfidx = cms2.mus_trkidx().at(index);
        if (ctfidx < 0 || vtxidx < 0) return false;
        const std::pair<double, double> cord0 = trks_d0_pv(ctfidx, vtxidx);
        const std::pair<double, double> cordz = trks_dz_pv(ctfidx, vtxidx);
        if (fabs(cord0.first) > 0.2) return false;
        if (fabs(cordz.first) > 0.5) return false;
        if (cms2.trks_valid_pixelhits().at(ctfidx) < 1) return false;
        if (cms2.trks_nlayers().at(ctfidx) < 6) return false;
        return true;
		}
		break;

	default: {/*do nothing*/}
    } // end switch block

    return false;
}

