#include <assert.h>
#include <algorithm>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "electronSelections.h"
#include "electronSelectionsParameters.h"
#include "muonSelections.h"
#include "metSelections.h"




#include "triggerUtils.h"
#include "CMS2.h"

#include "ttbarSelections.h"

using namespace tas;


///******************************************************************************************/     
////return the MET and the MET phi, correcting for mus that are not corrected for by default
///******************************************************************************************/     
//std::pair<float,float> getMet(const string algo, unsigned int hypIdx) {
//  
//  if(algo != "tcMET" && algo != "muCorMET" && algo != "pfMET" && algo != "tcMET35X" && algo != "tcMET_looper") {
//    cout << algo << "IS NOT A RECOGNIZED MET ALGORITHM!!!!! PLEASE CHECK YOUR CODE!!!";
//    return make_pair(-99999., -99999.);
//  }
//
//  
//  if(algo == "tcMET") {
//
//    float tcmetX = evt_tcmet()*cos(evt_tcmetPhi());
//    float tcmetY = evt_tcmet()*sin(evt_tcmetPhi());
//    
//    if(abs(hyp_lt_id()[hypIdx]) == 13)
//      fixMetForThisMuon(hyp_lt_index().at(hypIdx), tcmetX, tcmetY, usingTcMet);
//    if(abs(hyp_ll_id()[hypIdx]) == 13)
//      fixMetForThisMuon(hyp_ll_index().at(hypIdx), tcmetX, tcmetY, usingTcMet);
//    
//    return make_pair(sqrt(tcmetX*tcmetX + tcmetY*tcmetY), atan2(tcmetY, tcmetX));
//  }
///*
//  if(algo == "tcMET_looper") {
//
//    metStruct tcmetStruct = correctedTCMET();
//    float tcmetX = tcmetStruct.metx;
//    float tcmetY = tcmetStruct.mety;
//    
//    if(abs(hyp_lt_id()[hypIdx]) == 13)
//      fixMetForThisMuon(hyp_lt_index().at(hypIdx), tcmetX, tcmetY, usingTcMet_looper);
//    if(abs(hyp_ll_id()[hypIdx]) == 13)
//      fixMetForThisMuon(hyp_ll_index().at(hypIdx), tcmetX, tcmetY, usingTcMet_looper);
//    
//    return make_pair(sqrt(tcmetX*tcmetX + tcmetY*tcmetY), atan2(tcmetY, tcmetX));
//  }
//
//
//  if(algo == "tcMET35X") {
//
//    float tcmetX = evt35X_tcmet()*cos(evt35X_tcmetPhi());
//    float tcmetY = evt35X_tcmet()*sin(evt35X_tcmetPhi());
//    
//    if(abs(hyp_lt_id()[hypIdx]) == 13)
//      fixMetForThisMuon(hyp_lt_index().at(hypIdx), tcmetX, tcmetY, usingTcMet35X);
//    if(abs(hyp_ll_id()[hypIdx]) == 13)
//      fixMetForThisMuon(hyp_ll_index().at(hypIdx), tcmetX, tcmetY, usingTcMet35X);
//    
//    return make_pair(sqrt(tcmetX*tcmetX + tcmetY*tcmetY), atan2(tcmetY, tcmetX));
//  }
//*/
//
//  
//  if(algo == "muCorMET") {
//
//    float metX = evt_metMuonCorr()*cos(evt_metMuonCorrPhi());
//    float metY = evt_metMuonCorr()*sin(evt_metMuonCorrPhi());
//    
//    if(abs(hyp_lt_id()[hypIdx]) == 13)
//      fixMetForThisMuon(hyp_lt_index().at(hypIdx), metX, metY, usingCaloMet);
//    if(abs(hyp_ll_id()[hypIdx]) == 13)
//      fixMetForThisMuon(hyp_ll_index().at(hypIdx), metX, metY, usingCaloMet);
//
//    return make_pair(sqrt(metX*metX + metY*metY), atan2(metY, metX));
//  }
//  
//  //nothing to do here because they're perfect
//  if(algo == "pfMET") 
//    return make_pair(evt_pfmet(), evt_pfmetPhi());
//  
//  
//  return make_pair(-99999., -99999);
//  
//}


/*****************************************************************************************/
//hypothesis disambiguation. Returns the hypothesis that has the highest sum Pt
/*****************************************************************************************/
unsigned int selectHypByHighestSumPt(const vector<unsigned int> &v_goodHyps) {
  
  float maxSumPt = 0;
  unsigned int bestHypIdx = 0;
  for(unsigned int i = 0; i < v_goodHyps.size(); i++) {
    
    unsigned int index = v_goodHyps.at(i);
    float sumPt = hyp_lt_p4()[index].Pt() + hyp_ll_p4()[index].Pt();
    if( sumPt > maxSumPt) {
      maxSumPt = sumPt;
      bestHypIdx = index;
    }
  }

  return bestHypIdx;

}



