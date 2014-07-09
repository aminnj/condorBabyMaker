#include "./QuarkGluonTagger.h"

// Liner Radial Momentum taken from http://arxiv.org/abs/1106.3076v2

float getLRM(int ijet , int power) {

  float sum_pt_ri=0;
  float sum_pt=0;

  vector<int> cands = pfjets_pfcandIndicies().at(ijet);

  for (unsigned int ivc=0; ivc<cands.size(); ivc++) {
    
    int ican = cands[ivc];

    if ( pfcands_charge().at(ican)==0 ) continue;

    /////// dz constraint

    int itrk = -1;
    float dz=9999.;

    if (abs(pfcands_particleId().at(ican))!=11) {
      itrk = pfcands_trkidx().at(ican);
      if( itrk >= (int)trks_trk_p4().size() || itrk < 0 ) continue;
      dz=trks_dz_pv(itrk,0).first;
    }

    if (abs(pfcands_particleId().at(ican))==11 && pfcands_pfelsidx().at(ican)>=0) {
      itrk = els_gsftrkidx().at(pfcands_pfelsidx().at(ican));
      if( itrk >= (int)gsftrks_p4().size() || itrk < 0 ) continue;
      dz=gsftrks_dz_pv(itrk,0).first;
    }

    if(fabs(dz)>0.05) continue; 

    /////// STORE pt and dr

    float pt = pfcands_p4().at(ican).pt();
    float ri=deltaR(pfcands_p4().at(ican).Rapidity(), pfcands_p4().at(ican).Phi(), pfjets_p4().at(ijet).Rapidity(), pfjets_p4().at(ijet).Phi());

    sum_pt_ri += pow(pt,power) * pow(ri,power);
    sum_pt += pow(pt,power) ;
    
    //    int ipfel = cms2.pfcands_pfelsidx().at(ican);

  }    

  float lrm = (sum_pt > 0.) ? (sum_pt_ri / sum_pt ) : 9999.;

  return lrm;

}

float constituentPtDistribution(int ijet ) {
 
  float sum_pt2 = 0.;
  float sum_pt  = 0.;  

  vector<int> cands = pfjets_pfcandIndicies().at(ijet);

  for (unsigned int ivc=0;ivc<cands.size();ivc++) {
    
    int ican = cands[ivc];

    float pt = pfcands_p4().at(ican).pt();
    float pt2 = pt*pt;
    
    sum_pt += pt;
    sum_pt2 += pt2;
    
    //    if (cms2.pfcands_charge().at(ican)==0) continue;
    //    int ipfel = cms2.pfcands_pfelsidx().at(ican);

  }    
    
  float ptD_value = (sum_pt>0.) ? sqrt( sum_pt2 / (sum_pt*sum_pt) ) : 0.;
 
  return ptD_value;
 
}

// this function port
// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/pandolf/QuarkGluonTagger/src/QuarkGluonTagger.cc?revision=1.9&view=markup
// the tagger is defined only for jets with pt>20 and eta < 2.5 and the design performance is achieved for jets with |eta| < 2.0
// here we need use the rho with EtaMax = 2.5
// taks as input the p4 of the corrected (L1,L2,L3) jet

float QGtagger(LorentzVector p4 , int ijet, QGLikelihoodCalculator * qglikeli_) {

  //  bool isData=true;
  
  double corPt=p4.pt();
  float ptD = constituentPtDistribution(ijet); 
  double rho = evt_kt6pf_foregiso_rho();

  int nCharged = pfjets_chargedHadronMultiplicity().at(ijet);
  int nNeutral = pfjets_neutralHadronMultiplicity().at(ijet)+pfjets_photonMultiplicity().at(ijet);
  
  ////  cout << "corPt " << corPt << " ptD " << ptD << " neu " << nNeutral  << " ch " << nCharged << endl;
  
  //// compute the LD:
  float qgl(-1.0);
  
  if (nCharged + nNeutral > 0 ) {
    if( fabs(p4.eta())<2.5  && corPt>20)
      qgl = qglikeli_->computeQGLikelihood(corPt,rho,nCharged,nNeutral,ptD);
    else
      qgl = -1.;
  }
  
  return qgl;
  

}


