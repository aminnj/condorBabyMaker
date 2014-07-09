//
// met selections
//

#include <math.h>
#include <algorithm>
#include "TMath.h"
#include "TVector2.h"

#include "CMS2.h"
#include "trackSelections.h"
#include "metSelections.h"
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"

#include "jetcorr/FactorizedJetCorrector.h"
#include "jetSelections.h"
#include "eventSelections.h"
#include <string>
#include "Math/PtEtaPhiE4D.h"
#include "Math/LorentzVector.h"

//---------------------------------------------
// function to calculate latest tcMET
//---------------------------------------------
#include "tcmet/getTcmetFromCaloMet.icc"
#include "tcmet/getResponseFunction_fit.icc"

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > LorentzVector2;

//----------------------------------------------------------
//this function takes met, and performs a type1 correction
//using the given jet collection and L2L3-corrections
//----------------------------------------------------------

struct mu_jet_dr {
    bool operator () (const std::pair<int, int> &v1, const std::pair<int, int>  &v2) 
    {
        float dr1  = ROOT::Math::VectorUtil::DeltaR(cms2.pfjets_p4().at(v1.first), cms2.mus_p4().at(v1.second));
        float dr2  = ROOT::Math::VectorUtil::DeltaR(cms2.pfjets_p4().at(v2.first), cms2.mus_p4().at(v2.second));
        return dr1 < dr2;
    }
};

metStruct customType1Met( float metx , float mety , float sumet  , VofP4 jets , vector<float> cors )
{

    for( unsigned int i = 0 ; i < jets.size() ; ++i ){
        LorentzVector vdiff = jets.at(i) * cors.at(i) - jets.at(i);
        metx  -= vdiff.x();
        mety  -= vdiff.y();
        sumet += vdiff.pt();
    }

    metStruct myStruct;
    myStruct.metx   = metx;
    myStruct.mety   = mety;
    myStruct.met    = sqrt( metx * metx + mety * mety );
    myStruct.sumet  = sumet;
    myStruct.metphi = atan2( mety , metx );

    return myStruct;

}


metStruct correctedTCMET(bool printout, ostream& ostr) 
{
    // static because we only want to get the response function once
    static TH2F* rf = getResponseFunction_fit();
    return getTcmetFromCaloMet(rf, printout,ostr);
}

//---------------------------------------------------
// Function that checks whether met (or tcmet) was 
// corrected for a given muon.  This uses the value maps
//---------------------------------------------------
bool wasMetCorrectedForThisMuon(int imu, whichMetType type) {
    bool answer=true;
    switch(type) {
        case usingTcMet:
            if (cms2.mus_tcmet_flag().at(imu) == 0 || 
                    cms2.mus_tcmet_flag().at(imu) == 4) answer = false;
            break;
        case usingTcMet_looper:
            if (cms2.mus_tcmet_flag().at(imu) == 0 || cms2.mus_ptErr()[imu] / cms2.mus_p4()[imu].pt() > 1 ||
                    cms2.mus_tcmet_flag().at(imu) == 4) answer = false;
            break;
        case usingTcMet35X:
            //    if (cms2.evt35X_mus_tcmet_flag().at(imu) == 0 || 
            // 	cms2.evt35X_mus_tcmet_flag().at(imu) == 4) answer = false;
            cout << "metSelections.cc: Now using 38x MC, usingTcMet35X does not exist" << endl;
            break;
        case usingCaloMet:
            if (cms2.mus_met_flag().at(imu) == 0) answer = false;
            break;
        default:
            std::cout << "Illegal call to wasMetCorrectedForThisMuon" <<std::endl;
    }

    return answer;
}

//-----------------------------------------------------------
// Function that corrects the met (or tcmet) for a given
// muon in case this was not done in reco.  Uses value maps
//-------------------------------------------------------------
void fixMetForThisMuon(int imu, float& metX, float& metY, whichMetType type) {
    bool wasItCorrected = wasMetCorrectedForThisMuon(imu, type);
    if (!wasItCorrected) {
        switch(type) {

            case usingTcMet:
                if (cms2.mus_tcmet_flag()[imu] == 0) {//not corrected
                    metX += cms2.mus_met_deltax()[imu] - cms2.mus_p4()[imu].x();
                    metY += cms2.mus_met_deltay()[imu] - cms2.mus_p4()[imu].y();
                } else if (cms2.mus_tcmet_flag()[imu] == 4) {
                    metX += - cms2.mus_tcmet_deltax()[imu] + cms2.trks_trk_p4()[cms2.mus_trkidx()[imu]].px() // undo the pion correction
                        + cms2.mus_met_deltax()[imu] - cms2.mus_p4()[imu].x(); // perform the muon correction
                    metY += - cms2.mus_tcmet_deltay()[imu] + cms2.trks_trk_p4()[cms2.mus_trkidx()[imu]].py() // undo the pion correction
                        + cms2.mus_met_deltay()[imu] - cms2.mus_p4()[imu].y(); // perform the muon correction
                }
                break;

            case usingTcMet_looper:
                if (cms2.mus_tcmet_flag()[imu] == 0 || cms2.mus_ptErr()[imu] / cms2.mus_p4()[imu].pt() > 1) {//not corrected
                    metX += cms2.mus_met_deltax()[imu] - cms2.mus_p4()[imu].x();
                    metY += cms2.mus_met_deltay()[imu] - cms2.mus_p4()[imu].y();
                } else if (cms2.mus_tcmet_flag()[imu] == 4) {
                    metX += - cms2.mus_tcmet_deltax()[imu] + cms2.trks_trk_p4()[cms2.mus_trkidx()[imu]].px() // undo the pion correction
                        + cms2.mus_met_deltax()[imu] - cms2.mus_p4()[imu].x(); // perform the muon correction
                    metY += - cms2.mus_tcmet_deltay()[imu] + cms2.trks_trk_p4()[cms2.mus_trkidx()[imu]].py() // undo the pion correction
                        + cms2.mus_met_deltay()[imu] - cms2.mus_p4()[imu].y(); // perform the muon correction
                }
                break;


            case usingTcMet35X:
                //      if (cms2.evt35X_mus_tcmet_flag()[imu] == 0) {//not corrected
                // 	metX += cms2.mus_met_deltax()[imu] - cms2.mus_p4()[imu].x();
                // 	metY += cms2.mus_met_deltay()[imu] - cms2.mus_p4()[imu].y();
                //      } else if (cms2.evt35X_mus_tcmet_flag()[imu] == 4) {
                // 	metX += - cms2.evt35X_mus_tcmet_deltax()[imu] + cms2.trks_trk_p4()[cms2.mus_trkidx()[imu]].px() // undo the pion correction
                // 	  + cms2.mus_met_deltax()[imu] - cms2.mus_p4()[imu].x(); // perform the muon correction
                // 	metY += - cms2.evt35X_mus_tcmet_deltay()[imu] + cms2.trks_trk_p4()[cms2.mus_trkidx()[imu]].py() // undo the pion correction
                // 	  + cms2.mus_met_deltay()[imu] - cms2.mus_p4()[imu].y(); // perform the muon correction
                //      }
                cout << "metSelections.cc: Now using 38x MC, usingTcMet35X does not exist" << endl;
                break;


            case usingCaloMet:
                metX += cms2.mus_met_deltax()[imu] - cms2.mus_p4()[imu].x();
                metY += cms2.mus_met_deltay()[imu] - cms2.mus_p4()[imu].y();
                break;
        }
    }
}

//-----------------------------------------------------------
// Function that corrects the met (or tcmet) for a given
// muon in case this was not done in reco.  Uses value maps
//-------------------------------------------------------------
void fixMetForThisMuon(int imu, float& metX, float& metY, float& sumET, whichMetType type) {
    bool wasItCorrected = wasMetCorrectedForThisMuon(imu, type);
    if (!wasItCorrected) {
        switch(type) {

            case usingTcMet:
                if (cms2.mus_tcmet_flag()[imu] == 0) {
                    metX += cms2.mus_met_deltax()[imu] - cms2.mus_p4()[imu].x();
                    metY += cms2.mus_met_deltay()[imu] - cms2.mus_p4()[imu].y();
                    sumET -= sqrt(cms2.mus_met_deltax()[imu] * cms2.mus_met_deltax()[imu] + cms2.mus_met_deltay()[imu] * cms2.mus_met_deltay()[imu]) - cms2.mus_p4()[imu].pt(); 
                } else if (cms2.mus_tcmet_flag()[imu] == 4) {
                    metX += - cms2.mus_tcmet_deltax()[imu] + cms2.trks_trk_p4()[cms2.mus_trkidx()[imu]].px() // undo the pion correction
                        + cms2.mus_met_deltax()[imu] - cms2.mus_p4()[imu].x(); // perform the muon correction
                    metY += - cms2.mus_tcmet_deltay()[imu] + cms2.trks_trk_p4()[cms2.mus_trkidx()[imu]].py() // undo the pion correction
                        + cms2.mus_met_deltay()[imu] - cms2.mus_p4()[imu].y(); // perform the muon correction
                    sumET -= sqrt(cms2.mus_met_deltax()[imu] * cms2.mus_met_deltax()[imu] + cms2.mus_met_deltay()[imu] * cms2.mus_met_deltay()[imu]) - cms2.mus_p4()[imu].pt()
                        + sqrt(cms2.mus_tcmet_deltax()[imu] * cms2.mus_tcmet_deltax()[imu] + cms2.mus_tcmet_deltay()[imu] * cms2.mus_tcmet_deltay()[imu]) + cms2.trks_trk_p4()[cms2.mus_trkidx()[imu]].pt();
                }
                break;

            case usingTcMet_looper:
                if (cms2.mus_tcmet_flag()[imu] == 0 || cms2.mus_ptErr()[imu] / cms2.mus_p4()[imu].pt() > 1) {
                    metX += cms2.mus_met_deltax()[imu] - cms2.mus_p4()[imu].x();
                    metY += cms2.mus_met_deltay()[imu] - cms2.mus_p4()[imu].y();
                    sumET -= sqrt(cms2.mus_met_deltax()[imu] * cms2.mus_met_deltax()[imu] + cms2.mus_met_deltay()[imu] * cms2.mus_met_deltay()[imu]) - cms2.mus_p4()[imu].pt(); 
                } else if (cms2.mus_tcmet_flag()[imu] == 4) {
                    metX += - cms2.mus_tcmet_deltax()[imu] + cms2.trks_trk_p4()[cms2.mus_trkidx()[imu]].px() // undo the pion correction
                        + cms2.mus_met_deltax()[imu] - cms2.mus_p4()[imu].x(); // perform the muon correction
                    metY += - cms2.mus_tcmet_deltay()[imu] + cms2.trks_trk_p4()[cms2.mus_trkidx()[imu]].py() // undo the pion correction
                        + cms2.mus_met_deltay()[imu] - cms2.mus_p4()[imu].y(); // perform the muon correction
                    sumET -= sqrt(cms2.mus_met_deltax()[imu] * cms2.mus_met_deltax()[imu] + cms2.mus_met_deltay()[imu] * cms2.mus_met_deltay()[imu]) - cms2.mus_p4()[imu].pt()
                        + sqrt(cms2.mus_tcmet_deltax()[imu] * cms2.mus_tcmet_deltax()[imu] + cms2.mus_tcmet_deltay()[imu] * cms2.mus_tcmet_deltay()[imu]) + cms2.trks_trk_p4()[cms2.mus_trkidx()[imu]].pt();
                }
                break;


            case usingTcMet35X:
                //      if (cms2.evt35X_mus_tcmet_flag()[imu] == 0) {
                //        metX += cms2.mus_met_deltax()[imu] - cms2.mus_p4()[imu].x();
                //        metY += cms2.mus_met_deltay()[imu] - cms2.mus_p4()[imu].y();
                //        sumET -= sqrt(cms2.mus_met_deltax()[imu] * cms2.mus_met_deltax()[imu] + cms2.mus_met_deltay()[imu] * cms2.mus_met_deltay()[imu]) - cms2.mus_p4()[imu].pt(); 
                //      } else if (cms2.mus_tcmet_flag()[imu] == 4) {
                //        metX += - cms2.evt35X_mus_tcmet_deltax()[imu] + cms2.trks_trk_p4()[cms2.mus_trkidx()[imu]].px() // undo the pion correction
                //          + cms2.mus_met_deltax()[imu] - cms2.mus_p4()[imu].x(); // perform the muon correction
                //        metY += - cms2.evt35X_mus_tcmet_deltay()[imu] + cms2.trks_trk_p4()[cms2.mus_trkidx()[imu]].py() // undo the pion correction
                //          + cms2.mus_met_deltay()[imu] - cms2.mus_p4()[imu].y(); // perform the muon correction
                //        sumET -= sqrt(cms2.mus_met_deltax()[imu] * cms2.mus_met_deltax()[imu] + cms2.mus_met_deltay()[imu] * cms2.mus_met_deltay()[imu]) - cms2.mus_p4()[imu].pt()
                //          + sqrt(cms2.evt35X_mus_tcmet_deltax()[imu] * cms2.evt35X_mus_tcmet_deltax()[imu] + cms2.evt35X_mus_tcmet_deltay()[imu] * cms2.evt35X_mus_tcmet_deltay()[imu]) + cms2.trks_trk_p4()[cms2.mus_trkidx()[imu]].pt();
                //      }
                cout << "metSelections.cc: Now using 38x MC, usingTcMet35X does not exist" << endl;
                break;



            case usingCaloMet:
                metX += cms2.mus_met_deltax()[imu] - cms2.mus_p4()[imu].x();
                metY += cms2.mus_met_deltay()[imu] - cms2.mus_p4()[imu].y();
                break;
        }
    }
}

//---------------------------------------------
// function to calculate projected MET.
// takes three parameters as input:
//
// 1. met
// 2. met phi
// 3. hypothesis index
//---------------------------------------------
float projectedMET( float met, float metPhi, int hyp_index ) {

    float deltaPhi = nearestHypLeptonPhi(metPhi, hyp_index);
    return ((deltaPhi < TMath::Pi() / 2.) ? met * sin(deltaPhi) : met);
}

//---------------------------------------------
// as above but simpler for single lepton events
//---------------------------------------------
float projectedMETW( float met, float metPhi, float leptonPhi) {
    float deltaPhi = acos(cos(metPhi - leptonPhi));
    return ((deltaPhi < TMath::Pi() / 2.) ? met * sin(deltaPhi) : met);
}


//---------------------------------------------
// utility function find deltaPhi between met
// and nearest hypothesis lepton
//---------------------------------------------
float nearestHypLeptonPhi( float metPhi, int hyp_index ) {

    //WARNING!  This was designed to work in a dilepton environment - NOT a trilepton 
    float tightDPhi = min(fabs(cms2.hyp_lt_p4()[hyp_index].phi() - metPhi), (float)(2 * TMath::Pi()) - fabs(cms2.hyp_lt_p4()[hyp_index].phi() - metPhi));
    float looseDPhi = min(fabs(cms2.hyp_ll_p4()[hyp_index].phi() - metPhi), (float)(2 * TMath::Pi()) - fabs(cms2.hyp_ll_p4()[hyp_index].phi() - metPhi));

    return min(tightDPhi, looseDPhi);

}

//---------------------------------------------
// correct tcMET for any hypothesis muons
// that have not been corrected for
//---------------------------------------------
metStruct correctTCMETforHypMuons (int hyp_index, float met_x, float met_y, float sumet)
{
    metStruct tcmetStruct;
    tcmetStruct.met     = sqrt(met_x * met_x + met_y * met_y);
    tcmetStruct.metphi  = atan2(met_y, met_x);
    tcmetStruct.metx    = met_x;
    tcmetStruct.mety    = met_y;
    tcmetStruct.sumet   = sumet; 

    if (cms2.hyp_type()[hyp_index] ==3)
        return tcmetStruct;

    unsigned int i_lt = cms2.hyp_lt_index()[hyp_index];
    unsigned int i_ll = cms2.hyp_ll_index()[hyp_index];

    if (abs(cms2.hyp_lt_id()[hyp_index]) == 13)
    {
        if(cms2.mus_tcmet_flag()[i_lt] == 0)
        {
            met_x += cms2.mus_met_deltax()[i_lt] - cms2.mus_p4()[i_lt].x();
            met_y += cms2.mus_met_deltay()[i_lt] - cms2.mus_p4()[i_lt].y();
        }
        else if (cms2.mus_tcmet_flag()[i_lt] == 4)
        {
            met_x += -cms2.mus_tcmet_deltax()[i_lt] + cms2.mus_met_deltax()[i_lt] - cms2.mus_p4()[i_lt].x() + cms2.trks_trk_p4()[cms2.mus_trkidx()[i_lt]].x(); 
            met_y += -cms2.mus_tcmet_deltay()[i_lt] + cms2.mus_met_deltay()[i_lt] - cms2.mus_p4()[i_lt].y() + cms2.trks_trk_p4()[cms2.mus_trkidx()[i_lt]].y(); 
        }
    }
    if (abs(cms2.hyp_ll_id()[hyp_index]) == 13)
    {
        if(cms2.mus_tcmet_flag()[i_ll] == 0)
        { 
            met_x += cms2.mus_met_deltax()[i_ll] - cms2.mus_p4()[i_ll].x(); 
            met_y += cms2.mus_met_deltay()[i_ll] - cms2.mus_p4()[i_ll].y(); 
        }
        else if (cms2.mus_tcmet_flag()[i_ll] == 4)
        { 
            met_x += -cms2.mus_tcmet_deltax()[i_ll] + cms2.mus_met_deltax()[i_ll] - cms2.mus_p4()[i_ll].x() + cms2.trks_trk_p4()[cms2.mus_trkidx()[i_ll]].x();  
            met_y += -cms2.mus_tcmet_deltay()[i_ll] + cms2.mus_met_deltay()[i_ll] - cms2.mus_p4()[i_ll].y() + cms2.trks_trk_p4()[cms2.mus_trkidx()[i_ll]].y();  
        } 
    }

    tcmetStruct.met     = sqrt(met_x * met_x + met_y * met_y);
    tcmetStruct.metphi  = atan2(met_y, met_x);
    tcmetStruct.metx    = met_x;
    tcmetStruct.mety    = met_y;
    tcmetStruct.sumet   = sumet; 

    return tcmetStruct;
}

/*
   double dzPV(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv){
   return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt();
   }
 */

metStruct trackerMET( int hyp_index, double deltaZCut,
        const std::vector<LorentzVector>* jets )
{
    if ( cms2.vtxs_sumpt().empty() ) return metStruct();
    double pX(0), pY(0);
    pX -= cms2.hyp_lt_p4().at(hyp_index).px();
    pY -= cms2.hyp_lt_p4().at(hyp_index).py();
    pX -= cms2.hyp_ll_p4().at(hyp_index).px();
    pY -= cms2.hyp_ll_p4().at(hyp_index).py();

    for (unsigned int i=0; i<cms2.pfcands_particleId().size(); ++i){
        if ( cms2.pfcands_charge().at(i)==0 ) continue;
        if ( fabs(ROOT::Math::VectorUtil::DeltaR(cms2.pfcands_p4().at(i),cms2.hyp_lt_p4().at(hyp_index)))<0.1 ) continue;
        if ( fabs(ROOT::Math::VectorUtil::DeltaR(cms2.pfcands_p4().at(i),cms2.hyp_ll_p4().at(hyp_index)))<0.1 ) continue;
        if ( jets ){
            bool matched = false;
            for ( std::vector<LorentzVector>::const_iterator jet = jets->begin(); jet != jets->end(); ++jet )
                if ( fabs(ROOT::Math::VectorUtil::DeltaR(cms2.pfcands_p4().at(i),*jet))<0.5 ) matched=true;
            if (matched) continue;
        }

        int trkIndex = cms2.pfcands_trkidx().at(i);
        if (trkIndex<0) continue;
        double dzpv = dzPV(cms2.trks_vertex_p4()[trkIndex], cms2.trks_trk_p4()[trkIndex], cms2.vtxs_position().front());

        if ( fabs(dzpv) > deltaZCut) continue;

        pX -= cms2.pfcands_p4().at(i).px();
        pY -= cms2.pfcands_p4().at(i).py();
    }

    if (jets){
        for ( std::vector<LorentzVector>::const_iterator jet = jets->begin(); jet != jets->end(); ++jet ){
            pX -= jet->px();
            pY -= jet->py();
        }
    }
    metStruct met;
    met.met     = sqrt(pX * pX + pY * pY);
    met.metphi  = atan2(pY, pX);
    met.metx = pX;
    met.mety = pY;
    return met;
}

LorentzVector cmsReducedMET(LorentzVector sumJet, LorentzVector lep1, LorentzVector lep2, LorentzVector metP4, int version) 
{
    LorentzVector Q  = lep1 + lep2;
    float bisectorPhi = min(lep1.Phi(), lep2.Phi()) + ROOT::Math::VectorUtil::DeltaPhi(lep1, lep2)/2;

    TVector2 projDilepton(Q.Px(), Q.Py());
    TVector2 projSumJet(sumJet.Px(), sumJet.Py());
    TVector2 projMET(metP4.Pt()*cos(metP4.Phi()), metP4.Pt()*sin(metP4.Phi()));

    if (version == 1) {
        projDilepton = projDilepton.Rotate(-bisectorPhi);
        projSumJet   = projSumJet.Rotate(-bisectorPhi);
        projMET      = projMET.Rotate(-bisectorPhi);
    } else if (version == 2) {
        projDilepton = projDilepton.Rotate(-Q.Phi());
        projSumJet   = projSumJet.Rotate(-Q.Phi());
        projMET      = projMET.Rotate(-Q.Phi());
    }

    TVector2 unclustered = -1*projMET;                     // projDilepton - 1.*(projMET + projDilepton) + delta;
    TVector2 clustered   = projDilepton + 1.*projSumJet;   // + delta;

    TVector2 reducedMET = TVector2((fabs(unclustered.Px()) < fabs(clustered.Px()) ? unclustered.Px() : clustered.Px()),
            (fabs(unclustered.Py()) < fabs(clustered.Py()) ? unclustered.Py() : clustered.Py()));

    return LorentzVector(reducedMET.Px(), reducedMET.Py(), 0, reducedMET.Mod());
}

std::pair<float, float> cmsReducedMET_v2(LorentzVector lep1, LorentzVector lep2, const std::vector<LorentzVector> &jets)
{

    //define the leading and sub-leading lepton
    TVector2 lead, trail;
    if (lep1.Pt() > lep2.Pt()) {
        lead = TVector2(lep1.px(), lep1.py());
        trail = TVector2(lep2.px(), lep2.py());
    } else {
        lead = TVector2(lep2.px(), lep2.py());
        trail = TVector2(lep1.px(), lep1.py());
    }

    //define the thrust and dilepton
    TVector2 dil = lead+trail;
    TVector2 thr = lead-trail;
    float dphill = fabs(ROOT::Math::VectorUtil::DeltaPhi(lep1, lep2));

    //define the longitudinal and perpendicular axis
    TVector2 a_l, a_t;
    if (dphill >= TMath::Pi()/2) {
      a_l = thr.Unit();
      a_t = a_l.Rotate(TMath::Pi()/2);
      if(a_t * lead < 0) a_t *= -1;
    } else {
      a_t = dil.Unit();
      a_l = a_t.Rotate(TMath::Pi()/2);
      if(a_l * lead < 0) a_l *= -1;
    }

    //project the dilepton
    float dileptonProj_long = dil * a_l;
    float dileptonProj_perp = dil * a_t;

    //project the jet sum
    float sumJetProj_long = 0.;
    float sumJetProj_perp = 0.;
    for (unsigned int j = 0; j < jets.size(); ++j) {
        TVector2 jet(jets[j].Px(), jets[j].Py());
        sumJetProj_long += jet*a_l;
        sumJetProj_perp += jet*a_t;
    }

    //project the met
    TVector2 pfMET(cms2.evt_pfmet() * cos(cms2.evt_pfmetPhi()), cms2.evt_pfmet() * sin(cms2.evt_pfmetPhi()));
    //float metProj_long = pfMET * a_l;
    //float metProj_perp = pfMET * a_t;

    TVector2 uncl = pfMET + dil;
    float unclProj_long = uncl * a_l;
    float unclProj_perp = uncl * a_t;

    //take the minimum recoil possible depending on the event category type
    float recoilProj_long = min(sumJetProj_long, float(-1. * (unclProj_long)));
    recoilProj_long = min(recoilProj_long, float(0.));
    float recoilProj_perp = min(sumJetProj_perp, float(-1. * (unclProj_perp)));
    recoilProj_perp = min(recoilProj_perp, float(0.));  

    //
    // CMS MINIMIZED VERSION
    //

    // unclustered
    float unclRedMet_long = dileptonProj_long - 1.0 * unclProj_long;
    float unclRedMet_perp = dileptonProj_perp - 1.0 * unclProj_perp;
    // clustered
    float cluRedMet_long  = dileptonProj_long + 1.0 * sumJetProj_long;
    float cluRedMet_perp  = dileptonProj_perp + 1.0 * sumJetProj_perp;
    //float cluRedMet       = sqrt(pow(cluRedMet_long, 2) + pow(cluRedMet_perp, 2));   

    //
    // CMS INDEPEDENT MINIMIZATION VERSION
    //

    float reducedMETIndminRmet_long = (fabs(unclRedMet_long) < fabs(cluRedMet_long) ? unclRedMet_long : cluRedMet_long); 
    float reducedMETIndminRmet_perp = (fabs(unclRedMet_perp) < fabs(cluRedMet_perp) ? unclRedMet_perp : cluRedMet_perp); 
    //float redMETIndminRmet          = sqrt(pow(reducedMETIndminRmet_long, 2) + pow(reducedMETIndminRmet_perp, 2));
    TVector2 redMETIndminRmetxy        = reducedMETIndminRmet_long * a_l + reducedMETIndminRmet_perp * a_t;

    return std::make_pair<float, float> (reducedMETIndminRmet_long, reducedMETIndminRmet_perp);

}

//-----------------------------------------------------
// function to scale the hadronic component of the MET
//-----------------------------------------------------
std::pair<float, float> scaleMET(std::pair<float, float> p_met, LorentzVector p4_dilep, double rescale) {
    
    float met    = p_met.first;
    float metPhi = p_met.second;
    float metx   = met * cos(metPhi);
    float mety   = met * sin(metPhi);

    float lepx = p4_dilep.Px();
    float lepy = p4_dilep.Py();
      
    //hadronic component of MET (well, mostly), scaled
    float metHx     = (metx + lepx) * rescale;
    float metHy     = (mety + lepy) * rescale;
    float metNewx   = metHx - lepx;
    float metNewy   = metHy - lepy;
    float metNewPhi = atan2(metNewy, metNewx);

    return ( std::make_pair(sqrt(metNewx * metNewx + metNewy * metNewy), metNewPhi) );
}

MetCorrector::MetCorrector(std::vector<std::string> &list_of_files)
{
    //
    // expect a vector of strings in the following format:
    // 
    // for MC, expect three strings
    //
    // L1Fast
    // L2
    // L3
    //
    // for Data, expect four strings
    //
    // L1Fast
    // L2
    // L3
    // L2L3residual
    //
    size_t number_of_files = list_of_files.size();
    if (number_of_files < 3 || number_of_files > 4) {
        std::cout << "Invalid list of files.  Expect three files for MC, four for data." << std::endl;
        std::cout << "Note: The order of the files is important."                        << std::endl;
        std::cout << "For MC, please provide, in order: L1Fast, L2, L3."                 << std::endl;
        std::cout << "For data, please provide, in order: L1Fast, L2, L3, L2L3residual." << std::endl;
        return;
    }

    std::vector<std::string> offset_corrector_filenames;
    offset_corrector_filenames.push_back(list_of_files.at(0));

    std::vector<std::string> full_corrector_filenames;
    full_corrector_filenames = list_of_files;

    offset_corrector = makeJetCorrector(offset_corrector_filenames);
    full_corrector   = makeJetCorrector(full_corrector_filenames);
}

MetCorrector::~MetCorrector ()
{
    delete offset_corrector;
    delete full_corrector;
}

std::pair<float, float> MetCorrector::getCorrectedMET(std::pair<float, float> &uncorr_met)
{
    return correctMETforJES(uncorr_met);
}

std::pair<float, float> MetCorrector::getCorrectedMET()
{
    std::pair<float, float> uncorr_met = make_pair(cms2.evt_pfmet(), cms2.evt_pfmetPhi());
    return correctMETforJES(uncorr_met);
}

std::pair<float, float> MetCorrector::correctMETforJES(std::pair<float, float> uncorr_met)
{
    float met    = uncorr_met.first;
    float metPhi = uncorr_met.second;
    float metx   = met * cos(metPhi);
    float mety   = met * sin(metPhi);

    for (unsigned int idx = 0; idx < cms2.pfjets_p4().size(); idx++) {
        
        LorentzVector jetp4 = cms2.pfjets_p4().at(idx);

        //
        // veto events with EM fraction > 0.9
        float emfrac = (cms2.pfjets_chargedEmE().at(idx) + cms2.pfjets_neutralEmE().at(idx)) / cms2.pfjets_p4().at(idx).E();
        if (emfrac > 0.9) continue;

        LorentzVector2 tmpjetp4 = LorentzVector2(jetp4);
        if (fabs(tmpjetp4.eta()) > 9.9) continue;
        // if (tmpjetp4.eta() > 4.7) tmpjetp4.SetEta(4.7);
        // if (tmpjetp4.eta() < -4.7) tmpjetp4.SetEta(-4.7);   

        offset_corrector->setRho(cms2.evt_ww_rho_vor());
        offset_corrector->setJetA(cms2.pfjets_area().at(idx));
        offset_corrector->setJetEta(tmpjetp4.eta());
        offset_corrector->setJetPt(tmpjetp4.pt());
        float offset_corr = offset_corrector->getCorrection();        

        full_corrector->setRho(cms2.evt_ww_rho_vor());
        full_corrector->setJetA(cms2.pfjets_area().at(idx));
        full_corrector->setJetEta(tmpjetp4.eta());
        full_corrector->setJetPt(tmpjetp4.pt());
        float full_corr = full_corrector->getCorrection();

        for (unsigned int ipfc = 0; ipfc < cms2.pfjets_pfcandIndicies().at(idx).size(); ipfc++) {

            int index = cms2.pfjets_pfcandIndicies().at(idx).at(ipfc);
            if (index < 0) {
                std::cout << "Found a PF candidate in a PF jet with bad index..." << std::endl;
                continue;
            }

            //
            // remove SA or global muons from jets before correcting
            //
            if (abs(cms2.pfcands_particleId().at(index)) == 13) {

                int pfmusidx = cms2.pfcands_pfmusidx().at(index);
                if (pfmusidx < 0) {
                    std::cout << "Found a PF candidate with |id| == 13 with a bad pfmus index..." << std::endl;
                    continue;
                }
            
                int musidx = cms2.pfmus_musidx().at(pfmusidx);
                if (musidx < 0) {
                    std::cout << "Found a PF muon with a bad reco muon index..." << std::endl;
                    continue;
                }

                bool is_global     = !(((cms2.mus_type().at(musidx)) & (1<<1)) == 0);
                bool is_standalone = !(((cms2.mus_type().at(musidx)) & (1<<3)) == 0);            
                if (!is_global && !is_standalone) continue;
                jetp4 -= cms2.pfcands_p4().at(index);
            }
            else {
                //
                // now we need to look for pf cands that aren't IDed as muons but have non-null muonrefs, sigh...
                //
                int trkidx = cms2.pfcands_trkidx().at(index);
                for (unsigned int imu = 0; imu < cms2.mus_p4().size(); imu++) {
                    if (cms2.mus_trkidx().at(imu) != trkidx) continue;
                    bool is_global     = !(((cms2.mus_type().at(imu)) & (1<<1)) == 0);
                    bool is_standalone = !(((cms2.mus_type().at(imu)) & (1<<3)) == 0);
                    if (!is_global && !is_standalone) continue;
                    jetp4 -= cms2.pfcands_p4().at(index);
                }
            }
        } // end loop over jet constituents

        //
        // only correct MET for jets with corrected pt > 10 GeV
        //
        if (full_corr * jetp4.pt() < 10.) continue;

        metx += jetp4.px() * (offset_corr -full_corr);
        mety += jetp4.py() * (offset_corr - full_corr);
    }

    met    = sqrt(pow(metx, 2) + pow(mety, 2));
    metPhi = atan2(mety, metx);

    return make_pair(met, metPhi);    
}
