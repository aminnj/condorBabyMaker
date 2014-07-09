#include "Math/VectorUtil.h"
#include <math.h>
#include "CMS2.h"
#include "trackSelections.h"

// return a pair of d0, d0err of a ctf track with respect to a primary vertex
std::pair<double, double> trks_d0_pv (int itrk, int ipv)
{
    // assume the layout of the covariance matrix is (Vxx, Vxy, Vxz)
    //						      (Vyx, Vyy, ...)
    const double bx  = cms2.vtxs_position().at(ipv).x()   ;
    const double by  = cms2.vtxs_position().at(ipv).y()   ;
    const double vxx = cms2.vtxs_covMatrix().at(ipv).at(0);
    const double vxy = cms2.vtxs_covMatrix().at(ipv).at(1);
    const double vyy = cms2.vtxs_covMatrix().at(ipv).at(4);

    const double phi      = cms2.trks_trk_p4().at(itrk).phi();
    const double d0vtx    = cms2.trks_d0().at(itrk) - bx * sin(phi) + by * cos(phi);
    const double d0err    = cms2.trks_d0Err().at(itrk);
    const double phierr   = cms2.trks_phiErr().at(itrk);
    const double d0phicov = cms2.trks_d0phiCov().at(itrk);

    // we will let the optimizer take care of subexpression
    // elimination for this one...
    const double d0err2vtx = d0err * d0err 
        - 2 * (bx * cos(phi) + by * sin(phi)) * d0phicov
        + (bx * cos(phi) + by * sin(phi)) * (bx * cos(phi) + by * sin(phi)) * phierr * phierr
        + sin(phi) * sin(phi) * vxx + cos(phi) * cos(phi) * vyy
        - 2 * sin(phi) * cos(phi) * vxy;
    if (d0err2vtx >= 0) 
        return std::pair<double, double>(d0vtx, sqrt(d0err2vtx));

    std::cerr << "Oh no!  sigma^2(d0corr) < 0!" << std::endl;
    return std::pair<double, double>(d0vtx, -sqrt(-d0err2vtx));
}

// return a pair of d0, d0err of a gsf track with respect to a primary vertex
std::pair<double , double> gsftrks_d0_pv (int itrk, int ipv)
{
    // assume the layout of the covariance matrix is (Vxx, Vxy, Vxz)
    //						      (Vyx, Vyy, ...)
    const double bx  = cms2.vtxs_position().at(ipv).x()   ;
    const double by  = cms2.vtxs_position().at(ipv).y()   ;
    const double vxx = cms2.vtxs_covMatrix().at(ipv).at(0);
    const double vxy = cms2.vtxs_covMatrix().at(ipv).at(1);
    const double vyy = cms2.vtxs_covMatrix().at(ipv).at(4);

    const double phi      = cms2.gsftrks_p4().at(itrk).phi();
    const double d0vtx    = cms2.gsftrks_d0().at(itrk) - bx * sin(phi) + by * cos(phi);
    const double d0err    = cms2.gsftrks_d0Err().at(itrk);
    const double phierr   = cms2.gsftrks_phiErr().at(itrk);
    const double d0phicov = cms2.gsftrks_d0phiCov().at(itrk);

    // we will let the optimizer take care of subexpression
    // elimination for this one...
    const double d0err2vtx = d0err * d0err 
        - 2 * (bx * cos(phi) + by * sin(phi)) * d0phicov
        + (bx * cos(phi) + by * sin(phi)) * (bx * cos(phi) + by * sin(phi)) * phierr * phierr
        + sin(phi) * sin(phi) * vxx + cos(phi) * cos(phi) * vyy
        - 2 * sin(phi) * cos(phi) * vxy;
    if (d0err2vtx >= 0) 
        return std::pair<double, double>(d0vtx, sqrt(d0err2vtx));

    std::cerr << "Oh no!  sigma^2(d0corr) < 0!" << std::endl;
    return std::pair<double, double>(d0vtx, -sqrt(-d0err2vtx));
}

// return a pair of dz, dzerr of a ctf track with respect to a primary vertex
std::pair<double, double> trks_dz_pv (int itrk, int ipv)
{
    LorentzVector pv = cms2.vtxs_position().at(ipv);
    double pvxErr    = cms2.vtxs_xError().at(ipv)  ;
    double pvyErr    = cms2.vtxs_yError().at(ipv)  ;
    double pvzErr    = cms2.vtxs_zError().at(ipv)  ;

    //LorentzVector tkp = cms2.trks_trk_p4().at(itrk);
    //LorentzVector tkv = cms2.trks_vertex_p4().at(itrk);

    double phi        = cms2.trks_trk_p4().at(itrk).phi();
    double theta      = cms2.trks_trk_p4().at(itrk).theta();

    double ddzdpvx    = cos(phi)*1./tan(theta);
    double ddzdpvy    = sin(phi)*1./tan(theta);
    double ddzdphi    = -1*pv.x()*sin(phi)*1./tan(theta) + pv.y()*cos(phi)*1./tan(theta);
    double ddzdtheta  = -1*1/sin(theta)*1/sin(theta) * (pv.x()*cos(phi) + pv.y()*sin(phi));

    ddzdpvx   *= ddzdpvx;
    ddzdpvy   *= ddzdpvy;
    ddzdphi   *= ddzdphi;
    ddzdtheta *= ddzdtheta;

    double z0Err    = cms2.trks_z0Err().at(itrk);
    double phiErr   = cms2.trks_phiErr().at(itrk);
    double thetaErr = cms2.trks_etaErr().at(itrk)*sin(theta);

    z0Err    *= z0Err;
    phiErr   *= phiErr;
    thetaErr *= thetaErr;
    pvxErr   *= pvxErr;
    pvyErr   *= pvyErr;
    pvzErr   *= pvzErr;

    double value = cms2.trks_z0().at(itrk) - pv.z() + (pv.x()*cos(phi) + pv.y()*sin(phi) )*1./tan(theta);

    //note that the error does not account for correlations since we do not store the track covariance matrix
    double error = sqrt(z0Err + pvzErr + ddzdpvx*pvxErr + ddzdpvy*pvyErr + ddzdphi*phiErr + ddzdtheta*thetaErr);

    return std::pair<double, double>(value, error);
}

std::pair<double, double> gsftrks_dz_pv (int itrk, int ipv)
{
    LorentzVector pv = cms2.vtxs_position().at(ipv);
    double pvxErr    = cms2.vtxs_xError().at(ipv)  ;
    double pvyErr    = cms2.vtxs_yError().at(ipv)  ;
    double pvzErr    = cms2.vtxs_zError().at(ipv)  ;

    //LorentzVector tkp = cms2.gsftrks_p4().at(itrk);
    //LorentzVector tkv = cms2.gsftrks_vertex_p4().at(itrk);

    double phi   = cms2.gsftrks_p4().at(itrk).phi();
    double theta = cms2.gsftrks_p4().at(itrk).theta();

    double ddzdpvx   = cos(phi)*1./tan(theta);
    double ddzdpvy   = sin(phi)*1./tan(theta);
    double ddzdphi   = -1*pv.x()*sin(phi)*1./tan(theta) + pv.y()*cos(phi)*1./tan(theta);
    double ddzdtheta = -1*1/sin(theta)*1/sin(theta) * (pv.x()*cos(phi) + pv.y()*sin(phi));

    ddzdpvx   *= ddzdpvx;
    ddzdpvy   *= ddzdpvy;
    ddzdphi   *= ddzdphi;
    ddzdtheta *= ddzdtheta;

    double z0Err    = cms2.gsftrks_z0Err().at(itrk);
    double phiErr   = cms2.gsftrks_phiErr().at(itrk);
    double thetaErr = cms2.gsftrks_etaErr().at(itrk)*sin(theta);

    z0Err    *= z0Err;
    phiErr   *= phiErr;
    thetaErr *= thetaErr;
    pvxErr   *= pvxErr;
    pvyErr   *= pvyErr;
    pvzErr   *= pvzErr;

    double value = cms2.gsftrks_z0().at(itrk) - pv.z() + (pv.x()*cos(phi) + pv.y()*sin(phi) )*1./tan(theta);

    //note that the error does not account for correlations since we do not store the track covariance matrix
    double error = sqrt(z0Err + pvzErr + ddzdpvx*pvxErr + ddzdpvy*pvyErr + ddzdphi*phiErr + ddzdtheta*thetaErr);

    return std::pair<double, double>(value, error);
}


//----------------------------------------------------------------
// Simple function that tells you whether or not a track passed 
// a particular quality flag.
//----------------------------------------------------------------
bool isTrackQuality( int index, int cuts ) {
    return ( ( cms2.trks_qualityMask().at(index) & cuts ) == cuts );
}


float ctfIsoValuePF(const unsigned int itrk, unsigned int ivtx, float coner, float minptn, float dzcut) {

    float trkdz = trks_dz_pv(itrk,ivtx).first;
    //float trketa = cms2.trks_trk_p4().at(itrk).eta();

    float pfciso = 0.;
    float pfniso = 0.;
    //float pffootprint = 0.;
    for (unsigned int ipf = 0; ipf < cms2.pfcands_p4().size(); ++ipf) {

        float dR = ROOT::Math::VectorUtil::DeltaR(cms2.pfcands_p4().at(ipf), cms2.trks_trk_p4().at(itrk));

        if (dR > coner) continue;

        float pfpt  = cms2.pfcands_p4().at(ipf).pt();
        //float pfeta = cms2.pfcands_p4().at(ipf).eta();    
        //float deta  = fabs(pfeta - trketa);
        //int pfid    = abs(cms2.pfcands_particleId().at(ipf));
        if (cms2.pfcands_charge().at(ipf) == 0) { 
            if (pfpt > minptn)
                pfniso += pfpt;           
        }
        else {
            int pftkid = cms2.pfcands_trkidx().at(ipf);
            if (pftkid >= 0 && static_cast<int>(itrk) == pftkid)
                continue;

            if (pftkid >= 0) {
                if ( fabs(trks_dz_pv(cms2.pfcands_trkidx().at(ipf), ivtx).first - trkdz) < dzcut) 
                    pfciso += pfpt;
            }
        }
    } // end loop over particle flow candidates

    return (pfciso+pfniso)/cms2.trks_trk_p4().at(itrk).pt();
}

int associateTrackToVertex(const unsigned int itrk) {
    
    int best_vtx = -1;
    float dz = -999.;
    for (unsigned int ivtx = 0; ivtx < cms2.vtxs_position().size(); ivtx++) {
        float tmp_dz = trks_dz_pv(itrk, ivtx).first;
        if (fabs(tmp_dz) > fabs(dz))
            continue;

        dz = tmp_dz;
        best_vtx = ivtx;
    }

    return best_vtx;
}
