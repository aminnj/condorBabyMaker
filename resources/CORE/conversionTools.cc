#include "CMS2.h"
#include "conversionTools.h"


using namespace std;

//-----------------------------------------------------------------------------
// construct
ConversionInfo::ConversionInfo()
    : dist_              ( -9999.0)
    , dcot_              ( -9999.0)
    , radiusOfConversion_( -9999.0)
    , pointOfConversion_ ( -9999.0,-9999.0,-9999.0)
    , ctfPartnerIndex_   ( -9999)
    , gsfPartnerIndex_   ( -9999)
    , deltaMissingHits_  ( -9999)
    , flag_              ( -9999)
{
}

ConversionInfo::ConversionInfo
( 
    double p_dist,
    double p_dcot,
    double p_radiusOfConversion,
    const XYZPoint& p_pointOfConversion,      
    int p_ctfPartnerIndex,
    int p_gsfPartnerIndex,
    int p_deltaMissingHits,
    int p_flag
)
    : dist_              (p_dist)
    , dcot_              (p_dcot)
    , radiusOfConversion_(p_radiusOfConversion)
    , pointOfConversion_ (p_pointOfConversion)
    , ctfPartnerIndex_   (p_ctfPartnerIndex)
    , gsfPartnerIndex_   (p_gsfPartnerIndex)
    , deltaMissingHits_  (p_deltaMissingHits)
    , flag_              (p_flag)
{
}

// destory
ConversionInfo::~ConversionInfo()
{
}

//-----------------------------------------------------------------------------
//Code which does the work of looping over the tracks and searching for the 
//partners
std::vector<ConversionInfo> getConversionInfos(const int gsfElectronIdx,                
                           const double bFieldAtOrigin,
                           const double minFracSharedHits) {


    //the electron's CTF track must share at least 45% of the inner hits
    //with the electron's GSF track
    int elctfidx = -999;
    int elgsfidx = cms2.els_gsftrkidx()[gsfElectronIdx];
    if(cms2.els_trkshFrac()[gsfElectronIdx] > minFracSharedHits)
        elctfidx = cms2.els_trkidx()[gsfElectronIdx];

    //these vectors are for those candidate partner tracks that pass our cuts
    vector<ConversionInfo> v_candidatePartners;
    //track indices required to make references
    int ctftk_i = 0;
    int gsftk_i = 0;


    //loop over the CTF tracks and try to find the partner track
    for( ; ctftk_i < static_cast<int>(cms2.trks_trk_p4().size()); ctftk_i++) {

        if(ctftk_i == elctfidx)
            continue;

        //apply quality cuts to remove bad tracks
        if(cms2.trks_ptErr()[ctftk_i]/cms2.trks_trk_p4()[ctftk_i].Pt() > 0.05)
            continue;
        if(cms2.trks_validHits()[ctftk_i] < 5)
            continue;

        if(elctfidx > -1) {      
            if(fabs(cms2.trks_trk_p4()[ctftk_i].Pt() - cms2.trks_trk_p4()[elctfidx].Pt())/cms2.trks_trk_p4()[elctfidx].Pt() < 0.2)
                continue;
        }

        //use the electron's CTF track, if not null, to search for the partner track
        //look only in a cone of 0.5 to save time, and require that the track is opp. sign
        if(elctfidx > -1 ) {
            if(ROOT::Math::VectorUtil::DeltaR(cms2.trks_trk_p4()[elctfidx], cms2.trks_trk_p4()[ctftk_i]) < 0.5 &&
                    (cms2.trks_charge()[elctfidx] + cms2.trks_charge()[ctftk_i] == 0) ) {

                ConversionInfo convInfo = getConversionInfo(cms2.trks_trk_p4()[elctfidx], 
                        cms2.trks_charge()[elctfidx],
                        cms2.trks_d0()[elctfidx],
                        cms2.trks_z0()[elctfidx],
                        cms2.trks_trk_p4()[ctftk_i], 
                        cms2.trks_charge()[ctftk_i],
                        cms2.trks_d0()[ctftk_i],
                        bFieldAtOrigin);

                //need to add the track reference information for completeness
                //because the overloaded fnc above does not make a trackRef
                int deltaMissingHits = cms2.trks_exp_innerlayers()[ctftk_i] - cms2.trks_exp_innerlayers()[elctfidx];
                convInfo = ConversionInfo(convInfo.dist(),
                        convInfo.dcot(),
                        convInfo.radiusOfConversion(),
                        convInfo.pointOfConversion(),
                        ctftk_i,
                        -9999,
                        deltaMissingHits,
                        0);

                v_candidatePartners.push_back(convInfo);
            }
        }//using the electron's CTF track



        //now we check using the electron's gsf track
        if(ROOT::Math::VectorUtil::DeltaR(cms2.gsftrks_p4()[elgsfidx], cms2.trks_trk_p4()[ctftk_i]) < 0.5 &&
                (cms2.gsftrks_charge()[elgsfidx] + cms2.trks_charge()[ctftk_i] == 0) &&
                cms2.gsftrks_ptErr()[elgsfidx]/cms2.gsftrks_p4()[elgsfidx].Pt() < 0.25) {

            int deltaMissingHits    = cms2.trks_exp_innerlayers()[ctftk_i] - cms2.gsftrks_exp_innerlayers()[elgsfidx];
            ConversionInfo convInfo = getConversionInfo(cms2.gsftrks_p4()[elgsfidx], 
                    cms2.gsftrks_charge()[elgsfidx],
                    cms2.gsftrks_d0()[elgsfidx],
                    cms2.gsftrks_z0()[elgsfidx],
                    cms2.trks_trk_p4()[ctftk_i],
                    cms2.trks_charge()[ctftk_i],
                    cms2.trks_d0()[ctftk_i],
                    bFieldAtOrigin);
            convInfo = ConversionInfo(convInfo.dist(),
                    convInfo.dcot(),
                    convInfo.radiusOfConversion(),
                    convInfo.pointOfConversion(),
                    ctftk_i,
                    -9999,
                    deltaMissingHits,
                    1);

            v_candidatePartners.push_back(convInfo);
        }//using the electron's GSF track

    }//loop over the CTF track collection


    //------------------------------------------------------ Loop over GSF collection ----------------------------------//
    for( ; gsftk_i < static_cast<int>(cms2.gsftrks_p4().size()); gsftk_i++) {

        //reject the electron's own gsfTrack
        if(elgsfidx == gsftk_i)
            continue;

        //apply quality cuts to remove bad tracks
        if(cms2.gsftrks_ptErr()[gsftk_i]/(cms2.gsftrks_p4()[gsftk_i].Pt()) > 0.5)
            continue;
        if(cms2.gsftrks_validHits()[gsftk_i] < 5)
            continue;

        if(fabs(cms2.gsftrks_p4()[elgsfidx].Pt() - cms2.gsftrks_p4()[gsftk_i].Pt())/cms2.gsftrks_p4()[elgsfidx].Pt() < 0.25)
            continue;

        //try using the electron's CTF track first if it exists
        //look only in a cone of 0.5 around the electron's track
        //require opposite sign
        if(elctfidx > -1 ) {
            if(ROOT::Math::VectorUtil::DeltaR(cms2.trks_trk_p4()[elctfidx], cms2.gsftrks_p4()[gsftk_i]) < 0.5 &&
                    (cms2.trks_charge()[elctfidx] + cms2.gsftrks_charge()[gsftk_i] == 0)) {

                int deltaMissingHits = cms2.gsftrks_exp_innerlayers()[gsftk_i] - cms2.trks_exp_innerlayers()[elctfidx];
                ConversionInfo convInfo = getConversionInfo(cms2.trks_trk_p4()[elctfidx],
                        cms2.trks_charge()[elctfidx],
                        cms2.trks_d0()[elctfidx],
                        cms2.trks_z0()[elctfidx],
                        cms2.gsftrks_p4()[gsftk_i], 
                        cms2.gsftrks_charge()[gsftk_i],
                        cms2.gsftrks_d0()[gsftk_i],
                        bFieldAtOrigin);
                //fill the Ref info
                convInfo = ConversionInfo(convInfo.dist(),
                        convInfo.dcot(),
                        convInfo.radiusOfConversion(),
                        convInfo.pointOfConversion(),
                        -9999,
                        gsftk_i,
                        deltaMissingHits,
                        2);
                v_candidatePartners.push_back(convInfo);
            }
        }

        //use the electron's gsf track
        if(ROOT::Math::VectorUtil::DeltaR(cms2.gsftrks_p4()[elgsfidx], cms2.gsftrks_p4()[gsftk_i]) < 0.5 &&
                (cms2.gsftrks_charge()[elgsfidx] + cms2.gsftrks_charge()[gsftk_i] == 0) &&
                (cms2.gsftrks_ptErr()[elgsfidx]/cms2.gsftrks_p4()[elgsfidx].Pt() < 0.5)) {
            ConversionInfo convInfo = getConversionInfo(cms2.gsftrks_p4()[elgsfidx], 
                    cms2.gsftrks_charge()[elgsfidx],
                    cms2.gsftrks_d0()[elgsfidx],
                    cms2.gsftrks_z0()[elgsfidx],
                    cms2.gsftrks_p4()[gsftk_i], 
                    cms2.gsftrks_charge()[gsftk_i],
                    cms2.gsftrks_d0()[gsftk_i],
                    bFieldAtOrigin); 

            int deltaMissingHits = cms2.gsftrks_exp_innerlayers()[gsftk_i] - cms2.gsftrks_exp_innerlayers()[elgsfidx];
            convInfo = ConversionInfo(convInfo.dist(),
                    convInfo.dcot(),
                    convInfo.radiusOfConversion(),
                    convInfo.pointOfConversion(),
                    -9999,
                    gsftk_i,
                    deltaMissingHits,
                    3);

            v_candidatePartners.push_back(convInfo);
        }
    }//loop over the gsf track collection

    return v_candidatePartners;

}


//------------------------------------------------------------------------------------
//Function which calculates dist, dcot, radius of conversion and the 
//point of conversion
ConversionInfo getConversionInfo(const LorentzVector el_tk_p4, 
                 const int el_charge,
                 const float el_d0,
                 const float el_dz,
                 const LorentzVector cand_p4,
                 const int cand_charge,
                 const float cand_d0,
                 const double bFieldAtOrigin) {

  

  //now calculate the conversion related information
  double elCurvature = -0.3*bFieldAtOrigin*(el_charge/el_tk_p4.pt())/100.;
  double rEl = fabs(1./elCurvature);
  double xEl = -1*(1./elCurvature - el_d0)*sin(el_tk_p4.phi());
  double yEl = (1./elCurvature - el_d0)*cos(el_tk_p4.phi());


  double candCurvature = -0.3*bFieldAtOrigin*(cand_charge/cand_p4.pt())/100.;
  double rCand = fabs(1./candCurvature);
  double xCand = -1*(1./candCurvature - cand_d0)*sin(cand_p4.phi());
  double yCand = (1./candCurvature - cand_d0)*cos(cand_p4.phi());

  double d = sqrt(pow(xEl-xCand, 2) + pow(yEl-yCand , 2));
  double dist = d - (rEl + rCand);
  double dcot = 1./tan(el_tk_p4.theta()) - 1./tan(cand_p4.theta());

  //get the point of conversion
  double xa1 = xEl   + (xCand-xEl) * rEl/d;
  double xa2 = xCand + (xEl-xCand) * rCand/d;
  double ya1 = yEl   + (yCand-yEl) * rEl/d;
  double ya2 = yCand + (yEl-yCand) * rCand/d;

  double x=.5*(xa1+xa2);
  double y=.5*(ya1+ya2);
  double rconv = sqrt(pow(x,2) + pow(y,2));
  double z = el_dz + rEl*el_tk_p4.pz()*TMath::ACos(1-pow(rconv,2)/(2.*pow(rEl,2)))/el_tk_p4.Pt();

  XYZPoint convPoint(x, y, z);

  //now assign a sign to the radius of conversion
  float tempsign = el_tk_p4.px()*x + el_tk_p4.py()*y;
  tempsign = tempsign/fabs(tempsign);
  rconv = tempsign*rconv;

  //return an instance of ConversionInfo, but with a NULL track refs
  return ConversionInfo(dist, dcot, rconv, convPoint, -9999, -9999, -9999, -9999);

}



//------------------------------------------------------------------------------------
//ranks partners by kind, then passes to arbitration function
ConversionInfo findBestConversionMatch(const std::vector<ConversionInfo>& v_convCandidates) {
  using namespace std;

  if(v_convCandidates.size() == 0)
    return   ConversionInfo(-9999.,-9999.,-9999.,
                XYZPoint(-9999.,-9999.,-9999),
                -9999, -9999,
                -9999, -9999);


  if(v_convCandidates.size() == 1)
    return v_convCandidates.at(0);

  vector<ConversionInfo> v_0;
  vector<ConversionInfo> v_1;
  vector<ConversionInfo> v_2;
  vector<ConversionInfo> v_3;
  //loop over the candidates
  for(unsigned int i = 1; i < v_convCandidates.size(); i++) {
    ConversionInfo temp = v_convCandidates.at(i);

    if(temp.flag() == 0) {
      bool isConv = false;
      if(fabs(temp.dist()) < 0.02 &&
     fabs(temp.dcot()) < 0.02 &&
     temp.deltaMissingHits() < 3 &&
     temp.radiusOfConversion() > -2)
    isConv = true;
      if(sqrt(pow(temp.dist(),2) + pow(temp.dcot(),2)) < 0.05  &&
     temp.deltaMissingHits() < 2 &&
     temp.radiusOfConversion() > -2)
    isConv = true;

      if(isConv)
    v_0.push_back(temp);
    }

    if(temp.flag() == 1) {

      if(sqrt(pow(temp.dist(),2) + pow(temp.dcot(),2)) < 0.05  &&
     temp.deltaMissingHits() < 2 &&
     temp.radiusOfConversion() > -2)
    v_1.push_back(temp);
    }
    if(temp.flag() == 2) {

      if(sqrt(pow(temp.dist(),2) + pow(temp.dcot()*temp.dcot(),2)) < 0.05 &&
     temp.deltaMissingHits() < 2 &&
     temp.radiusOfConversion() > -2)
    v_2.push_back(temp);

    }
    if(temp.flag() == 3) {

      if(sqrt(temp.dist()*temp.dist() + temp.dcot()*temp.dcot()) < 0.05
     && temp.deltaMissingHits() < 2
     && temp.radiusOfConversion() > -2)
    v_3.push_back(temp);

    }

  }//candidate conversion loop

  //now do some arbitration

  //give preference to conversion partners found in the CTF collection
  //using the electron's CTF track
  if(v_0.size() > 0)
    return arbitrateConversionPartnersbyR(v_0);

  if(v_1.size() > 0)
    return arbitrateConversionPartnersbyR(v_1);

  if(v_2.size() > 0)
    return arbitrateConversionPartnersbyR(v_2);

  if(v_3.size() > 0)
    return arbitrateConversionPartnersbyR(v_3);


  //if we get here, we didn't find a candidate conversion partner that
  //satisfied even the loose selections
  //return the the closest partner by R
  return arbitrateConversionPartnersbyR(v_convCandidates);

}

//--------------------------------------------------------------------------------------------------------
//takes in a vector of candidate conversion partners
//and arbitrates between them returning the one with the
//smallest R=sqrt(dist*dist + dcot*dcot)
ConversionInfo arbitrateConversionPartnersbyR(const std::vector<ConversionInfo>& v_convCandidates) {

  if(v_convCandidates.size() == 1)
    return v_convCandidates.at(0);

  ConversionInfo arbitratedConvInfo = v_convCandidates.at(0);
  double R = sqrt(pow(arbitratedConvInfo.dist(),2) + pow(arbitratedConvInfo.dcot(),2));

  for(unsigned int i = 1; i < v_convCandidates.size(); i++) {
    ConversionInfo temp = v_convCandidates.at(i);
    double temp_R = sqrt(pow(temp.dist(),2) + pow(temp.dcot(),2));
    if(temp_R < R) {
      R = temp_R;
      arbitratedConvInfo = temp;
    }

  }

  return arbitratedConvInfo;

}
