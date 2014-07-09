// $Id: jetSelections.cc,v 1.40 2013/01/28 23:44:26 kelley Exp $

#include <algorithm>
#include <utility>
#include "Math/VectorUtil.h"
#include "jetSelections.h"
#include "trackSelections.h"
#include "jetcorr/JetCorrectorParameters.icc"
#include "jetcorr/FactorizedJetCorrector.icc"
#include "jetcorr/SimpleJetCorrector.icc"
#include "jetcorr/JetCorrectionUncertainty.icc"
#include "jetcorr/SimpleJetCorrectionUncertainty.icc"
#include "utilities.h"

using std::vector;
using std::pair;

// define this type for speed: allows us to get a vector of selected
// jets, potentially with a correction factor, without having to make
// copies
typedef vector<pair<const LorentzVector *, double> > jets_with_corr_t;

//-----------------------------------------------------------------------------
// These are the jetMet rescaling factors for the jet (+1 sigma or -1 sigma)
// flag == 0:  return 1
// flag == 1:  scale up   by 1 sigma
// flag == -1: scale down by 1 sigma
// Uses the ICHE2012 prescription from Filip
//------------------------------------------------------------------------------
float getJetMetSyst(int flag, float pt, float eta) {

  if (flag == 0) return 1.;
  if (flag > 1 || flag < -1) {
    std::cout << "Illegal call to getJetMetSyst" << std::endl;
    return 1.0;
  }

  float result;
  const int nBinsEta = 2 ;
  const int nBinsPt = 10 ;
  float uncert[nBinsEta][nBinsPt] =
    {
      { 12.0,  7.5,  6.3,  4.7,  3.7,  2.7,  2.6,  2.5,  2.4,  2.3  },
      { 20.0,  17.2, 14.5,  10.7,  7.9,  6.9,  6.0,  5.2,  4.7,  4.4  }
    };			
  int ipt = int(pt/10.);
  float rest = pt/10. - ipt;
  int ieta = 0;
  if (eta < -2.5 || eta > 2.5) {ieta = 1;}
  
  if (pt > 10. && pt < 100.) { 
    result = uncert[ieta][ipt-1] + rest*(uncert[ieta][ipt] - uncert[ieta][ipt-1]);
  } else if (pt >= 100.){
    result = uncert[ieta][9];
  } else {
    result = uncert[ieta][0];
  }		
  result = 1.0 + flag*result*0.01;
  return result;
}
  



// function to give us the indices of jets passing kinematic and cleaning cuts
static jets_with_corr_t getJets_fast (unsigned int i_hyp, enum JetType type, enum CleaningType cleaning,
                                      double deltaR, double min_pt, double max_eta, double rescale = 1.0, int systFlag = 0)
{
    // JPT, PF or calo jets?  Introduce this variable so we only have to decide once
    const vector<LorentzVector> *jets = 0;
    switch (type) {
    case JETS_TYPE_CALO_CORR: case JETS_TYPE_CALO_UNCORR:
        jets = &cms2.jets_p4();
        break;
    case JETS_TYPE_PF_CORR: case JETS_TYPE_PF_UNCORR: case JETS_TYPE_PF_FAST_CORR: case JETS_TYPE_PF_FAST_CORR_RESIDUAL:
        jets = &cms2.pfjets_p4();
        break;
#if haveGEN	  
    case JETS_TYPE_GEN:
        jets = &cms2.genjets_p4();
        break;
#endif
    }
    jets_with_corr_t ret;
    ret.reserve(jets->size()); // reserve so we don't have to realloc later, which is slow
    for (unsigned int i = 0; i < jets->size(); ++i) {
        //------------------------------------------------------------
        // min pt cut
        //------------------------------------------------------------
        double corr = 1;
        // CALO_CORR and PF_CORR need to be pt-corrected
        switch (type) {
        case JETS_TYPE_CALO_CORR:
	  corr = cms2.jets_cor().at(i) * rescale * getJetMetSyst(systFlag, cms2.jets_p4().at(i).pt(), cms2.jets_p4().at(i).eta());
            break;
        case JETS_TYPE_PF_CORR:
//            corr = cms2.pfjets_cor().at(i) * rescale * getJetMetSyst(systFlag, cms2.pfjets_p4().at(i).pt(), cms2.pfjets_p4().at(i).eta());
            corr = cms2.pfjets_cor().at(i) * rescale * getJetMetSyst(systFlag, cms2.pfjets_p4().at(i).pt() * cms2.pfjets_cor().at(i), cms2.pfjets_p4().at(i).eta());
            break;
        case JETS_TYPE_PF_FAST_CORR:
//            corr = cms2.pfjets_corL1FastL2L3().at(i) * rescale * getJetMetSyst(systFlag, cms2.pfjets_p4().at(i).pt(), cms2.pfjets_p4().at(i).eta());
            corr = cms2.pfjets_corL1FastL2L3().at(i) * rescale * getJetMetSyst(systFlag, cms2.pfjets_p4().at(i).pt() * cms2.pfjets_corL1FastL2L3().at(i), cms2.pfjets_p4().at(i).eta());
            break;
        case JETS_TYPE_PF_FAST_CORR_RESIDUAL:
//            corr = cms2.pfjets_corL1FastL2L3residual().at(i) * rescale * getJetMetSyst(systFlag, cms2.pfjets_p4().at(i).pt(), cms2.pfjets_p4().at(i).eta());
            corr = cms2.pfjets_corL1FastL2L3residual().at(i) * rescale * getJetMetSyst(systFlag, cms2.pfjets_p4().at(i).pt() * cms2.pfjets_corL1FastL2L3residual().at(i), cms2.pfjets_p4().at(i).eta());
            break;
        case JETS_TYPE_CALO_UNCORR: 
        case JETS_TYPE_PF_UNCORR:
#if haveGEN	  
        case JETS_TYPE_GEN:
#endif
            corr = rescale;
            break;
        }
        const double pt = jets->at(i).pt() * corr;
        if (pt < min_pt) 
            goto conti;
        //------------------------------------------------------------
        // max |eta| cut
        //------------------------------------------------------------
        if (fabs(jets->at(i).eta()) > max_eta)
            goto conti;
        //------------------------------------------------------------
        // lepton cleaning
        //------------------------------------------------------------
        switch (cleaning) {
        case JETS_CLEAN_NONE:
            break;
        case JETS_CLEAN_HYP_E_MU:
        {
            const LorentzVector &lt = cms2.hyp_lt_p4().at(i_hyp);
            if (ROOT::Math::VectorUtil::DeltaR(jets->at(i), lt) < deltaR)
                goto conti;
            const LorentzVector &ll = cms2.hyp_ll_p4().at(i_hyp);
            if (ROOT::Math::VectorUtil::DeltaR(jets->at(i), ll) < deltaR)
                goto conti;
            break;
        }
        case JETS_CLEAN_HYP_E:
        {
            const LorentzVector &lt = cms2.hyp_lt_p4().at(i_hyp);
            const int lt_id = cms2.hyp_lt_id().at(i_hyp);
            if (abs(lt_id) == 11 && ROOT::Math::VectorUtil::DeltaR(jets->at(i), lt) < deltaR)
                goto conti;
            const LorentzVector &ll = cms2.hyp_ll_p4().at(i_hyp);
            const int ll_id = cms2.hyp_ll_id().at(i_hyp);
            if (abs(ll_id) == 11 && ROOT::Math::VectorUtil::DeltaR(jets->at(i), ll) < deltaR)
                goto conti;
            break;
        }
        case JETS_CLEAN_SINGLE_E:
        {
            const LorentzVector &e = cms2.els_p4().at(i_hyp);
            if (ROOT::Math::VectorUtil::DeltaR(jets->at(i), e) < deltaR)
                goto conti;
            break;
        }
        case JETS_CLEAN_SINGLE_MU:
        {
            const LorentzVector &e = cms2.mus_p4().at(i_hyp);
            if (ROOT::Math::VectorUtil::DeltaR(jets->at(i), e) < deltaR)
                goto conti;
            break;
        }
        default:
            assert(false);
        }
        //------------------------------------------------------------
        // jet ID cuts
        //------------------------------------------------------------
        if (type == JETS_TYPE_CALO_CORR || type == JETS_TYPE_CALO_UNCORR){
            if (not passesCaloJetID(jets->at(i)))
                goto conti;
        }
        if (type == JETS_TYPE_PF_FAST_CORR || type == JETS_TYPE_PF_CORR || type == JETS_TYPE_PF_UNCORR || type == JETS_TYPE_PF_FAST_CORR_RESIDUAL){
            if (not passesPFJetID(i)) 
                goto conti;
        }
        //------------------------------------------------------------
        // jet passed all cuts
        //------------------------------------------------------------
        ret.push_back(pair<const LorentzVector *, double>(&jets->at(i), corr));
        continue;
        //------------------------------------------------------------
        // jet failed
        //------------------------------------------------------------
    conti:
        ret.push_back(pair<const LorentzVector *, double>(0, corr));
    }
    return ret;
}

struct jets_pt_gt {
    bool operator () (const LorentzVector &v1, const LorentzVector &v2) 
    {
        return v1.pt() > v2.pt();
    }
};

// functions that we let other people use
vector<LorentzVector> getJets (unsigned int i_hyp, bool sort_, 
                               enum JetType type, enum CleaningType cleaning,
                               double deltaR, double min_pt, double max_eta, double rescale, int systFlag)

{
  jets_with_corr_t jets = getJets_fast(i_hyp, type, cleaning, deltaR, min_pt, max_eta, rescale, systFlag);
    vector<LorentzVector> ret;
    ret.reserve(jets.size());
    for (unsigned int i = 0; i < jets.size(); ++i) {
        // correct the jet momentum if a corrected jet type was requested
        if (jets[i].first != 0)
            ret.push_back(*jets[i].first * jets[i].second);
    }
    if (sort_)
        sort(ret.begin(), ret.end(), jets_pt_gt());
    return ret;
}

std::vector<bool> getJetFlags (unsigned int i_hyp, enum JetType type, enum CleaningType cleaning,
                               double deltaR, double min_pt, double max_eta, double rescale,  int systFlag)
{
  jets_with_corr_t jets = getJets_fast(i_hyp, type, cleaning, deltaR, min_pt, max_eta, rescale, systFlag);
    vector<bool> ret;
    ret.reserve(jets.size());
    for (unsigned int i = 0; i < jets.size(); ++i) {
        // correct the jet momentum if a corrected jet type was requested
        ret.push_back(jets[i].first != 0);
    }
    return ret;
}


int nJets (unsigned int i_hyp, enum JetType type, enum CleaningType cleaning,
           double deltaR, double min_pt, double max_eta, double rescale,  int systFlag)
{
  jets_with_corr_t jets = getJets_fast(i_hyp, type, cleaning, deltaR, min_pt, max_eta, rescale, systFlag);
    int ret = 0;
    for (unsigned int i = 0; i < jets.size(); ++i) {
        // correct the jet momentum if a corrected jet type was requested
        if (jets[i].first != 0)
            ret++;
    }
    return ret;
}

double sumPt (unsigned int i_hyp, enum JetType type, enum CleaningType cleaning,
              double deltaR, double min_pt, double max_eta, double rescale,  int systFlag)
{
  jets_with_corr_t jets = getJets_fast(i_hyp, type, cleaning, deltaR, min_pt, max_eta, rescale, systFlag);
    double ret = 0;
    for (unsigned int i = 0; i < jets.size(); ++i) {
        // correct the jet momentum if a corrected jet type was requested
        if (jets[i].first != 0)
            ret += jets[i].first->pt() * jets[i].second;
    }
    return ret;
}

static FactorizedJetCorrector *jetCorrector = 0;

void setJetCorrector (FactorizedJetCorrector *jc) 
{
    if (jetCorrector != 0)
        delete jetCorrector;
    jetCorrector = jc;
}

class FactorizedJetCorrector *makeJetCorrector (const char *l2corr, 
                                                const char *l3corr, 
                                                const char *l2l3_residual_corr)
{
    std::vector<std::string> corrs;
    corrs.reserve(3);
    corrs.push_back(l2corr);
    corrs.push_back(l3corr);
    corrs.push_back(l2l3_residual_corr);
    return makeJetCorrector(corrs);
}

class FactorizedJetCorrector *makeJetCorrector (const std::vector<std::string> &corrs)
{
    vector<JetCorrectorParameters> vParam;
    for (std::vector<std::string>::const_iterator i = corrs.begin(), i_end = corrs.end();
         i != i_end; ++i) {
        // do some rigmarole to evaluate env variables in the strings
        const std::string cmd = "echo ";
        FILE *f = popen((cmd + *i).c_str(), "r");
        if (!f) {
            perror((std::string("Opening pipe to execute ") + cmd + *i).c_str());
            return 0;
        }
        char corr_name[1024];
        int s = fscanf(f, " %1024s\n", corr_name);
        if (s != 1) {
            perror("reading file list");
        }
        assert(s == 1);
        JetCorrectorParameters JetCorPar(corr_name);
        // printf("%s\n", corr_name);
        vParam.push_back(JetCorrectorParameters(corr_name));
    }
    return new FactorizedJetCorrector(vParam);
}

double jetCorrection (const LorentzVector &jet, FactorizedJetCorrector *jc)
{
    jc->setJetPt(jet.pt());
    jc->setJetEta(jet.eta());
    return jc->getCorrection();
}

double jetCorrection (const LorentzVector &jet)
{
    assert(jetCorrector != 0);
    jetCorrector->setJetPt(jet.pt());
    jetCorrector->setJetEta(jet.eta());
    return jetCorrector->getCorrection();
}

double jetCorrection (int ijet) 
{ 
    return jetCorrection(cms2.jets_p4()[ijet]); 
}

bool jetPassesLooseJetID(int ijet)
{
    if (fabs(cms2.jets_p4()[ijet].eta()) < 3.)
    {
        if (cms2.jets_fHPD()[ijet] > 0.98)
            return false;
        if (cms2.jets_emFrac()[ijet] < 0.01)
            return false;
        if (cms2.jets_n90Hits()[ijet] < 2)
            return false;
    }

    return true;
}

bool passesCaloJetID (const LorentzVector &jetp4)
{
    int jet_idx = -1;
    double minDR = 999;

    for (unsigned int i = 0; i < cms2.jets_p4().size(); i++)
    {
        double deltaR = ROOT::Math::VectorUtil::DeltaR(jetp4, cms2.jets_p4()[i]);

        if (deltaR < minDR)
        {
            minDR = deltaR;
            jet_idx = i;
        }
    }

    if (jet_idx < 0)
        return false;

    if (cms2.jets_emFrac()[jet_idx] < 0.01 || cms2.jets_fHPD()[jet_idx] > 0.98 || cms2.jets_n90Hits()[jet_idx] < 2)
        return false;

    return true;
}

bool passesPFJetID(unsigned int pfJetIdx) {

    float pfjet_chf_  = cms2.pfjets_chargedHadronE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();
    float pfjet_nhf_  = cms2.pfjets_neutralHadronE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();
    float pfjet_cef_  = cms2.pfjets_chargedEmE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();
    float pfjet_nef_  = cms2.pfjets_neutralEmE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();
    int   pfjet_cm_   = cms2.pfjets_chargedMultiplicity()[pfJetIdx];
    int   pfjet_mult_ = pfjet_cm_ + cms2.pfjets_neutralMultiplicity()[pfJetIdx];

    if (pfjet_nef_ >= 0.99)
        return false;
    if (pfjet_nhf_ >= 0.99)
        return false;
    if (pfjet_mult_ < 2)
        return false;

    if (fabs(cms2.pfjets_p4()[pfJetIdx].eta()) < 2.4)
    {
        if (pfjet_chf_ < 1e-6)
			return false;
        if (pfjet_cm_ < 1)
			return false;
        if (pfjet_cef_ >= 0.99)
			return false;
    }

    return true;
}  

//added by jgran March 19, 2014.
//used for all hadronic MT2
bool passesLoosePFJetID(unsigned int pfJetIdx) {

    float pfjet_chf_  = cms2.pfjets_chargedHadronE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();
    float pfjet_nhf_  = cms2.pfjets_neutralHadronE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();
    float pfjet_cef_  = cms2.pfjets_chargedEmE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();
    float pfjet_nef_  = cms2.pfjets_neutralEmE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();
    int   pfjet_cm_   = cms2.pfjets_chargedMultiplicity()[pfJetIdx];

    if (cms2.pfjets_pfcandIndicies().size() < 2)
        return false;
    if (pfjet_nef_ >= 0.99)
        return false;
    if (pfjet_nhf_ >= 0.99)
        return false;
    if (pfjet_cm_ < 1)
        return false;
    if (pfjet_chf_ < 1e-6)
        return false;
    if (pfjet_cef_ >= 0.99)
        return false;

    return true;
}  



float randomConeEventDensity(){

    //define the eta bins
    vector<float> etabins;
    for (int i=0;i<8;++i) etabins.push_back(-2.1+0.6*i);
    //define the phi bins
    vector<float> phibins;
    for (int i=0;i<10;i++) phibins.push_back(-TMath::Pi()+(2*i+1)*TMath::TwoPi()/20.);
  
    float etadist = etabins[1]-etabins[0];
    float phidist = phibins[1]-phibins[0];
    float etahalfdist = (etabins[1]-etabins[0])/2.;
    float phihalfdist = (phibins[1]-phibins[0])/2.;

    vector<float> sumPFNallSMDQ;
    sumPFNallSMDQ.reserve(80);
    for (unsigned int ieta=0;ieta<etabins.size();++ieta) {
        for (unsigned int iphi=0;iphi<phibins.size();++iphi) {
            float pfniso_ieta_iphi = 0;
            for (unsigned int ipf=0; ipf<cms2.pfcands_p4().size(); ++ipf){
                if (fabs(etabins[ieta]-cms2.pfcands_p4().at(ipf).eta())>etahalfdist) continue;
                if (fabs(deltaPhi(phibins[iphi],cms2.pfcands_p4().at(ipf).phi()))>phihalfdist) continue;
                pfniso_ieta_iphi+=cms2.pfcands_p4().at(ipf).pt();
            }
            sumPFNallSMDQ.push_back(pfniso_ieta_iphi);
        }
    }

    float evt_smdq = 0;
    sort(sumPFNallSMDQ.begin(),sumPFNallSMDQ.end());
  
    if (sumPFNallSMDQ.size()%2) evt_smdq = sumPFNallSMDQ[(sumPFNallSMDQ.size()-1)/2];
    else evt_smdq = (sumPFNallSMDQ[sumPFNallSMDQ.size()/2]+sumPFNallSMDQ[(sumPFNallSMDQ.size()-2)/2])/2.;
    return evt_smdq/(etadist*phidist);
}


float jetDz(int ijet, int ivtx) {
    //add protection against non-pf jets
    vector<int> cands = cms2.pfjets_pfcandIndicies().at(ijet);
    float jptsq   = 0;
    float jptsqdz = 0;
    for (unsigned int ivc=0;ivc<cands.size();ivc++) {
        int ican = cands[ivc];
        if (cms2.pfcands_charge().at(ican)==0) continue;
        int ipfel = cms2.pfcands_pfelsidx().at(ican);
        int iel=-1, igsf=-1;
        if (ipfel>=0) iel   = cms2.pfels_elsidx().at(ipfel);
        if (iel>=0) igsf  = cms2.els_gsftrkidx().at(iel);
        float ptsq = pow(cms2.pfcands_p4().at(ican).pt(),2);
        if (igsf>=0) {
            float dzc = gsftrks_dz_pv(igsf,ivtx).first;
            jptsq+=ptsq;
            jptsqdz+=ptsq*dzc;
        } else {
            int itrk = cms2.pfcands_trkidx().at(ican);
            if (itrk<0) continue;
            float dzc = trks_dz_pv(itrk,ivtx).first;
            jptsq+=ptsq;
            jptsqdz+=ptsq*dzc;
        }
    }
    if (jptsq>1E-6) return jptsqdz/jptsq;
    else return 99999.;
}

float dz_trk_vtx( const unsigned int trkidx, const unsigned int vtxidx ){
  
  //return ((cms2.trks_vertex_p4()[trkidx].z()-cms2.vtxs_position()[vtxidx].z()) - ((cms2.trks_vertex_p4()[trkidx].x()-cms2.vtxs_position()[vtxidx].x()) * cms2.trks_trk_p4()[trkidx].px() + (cms2.trks_vertex_p4()[trkidx].y() - cms2.vtxs_position()[vtxidx].y()) * cms2.trks_trk_p4()[trkidx].py())/cms2.trks_trk_p4()[trkidx].pt() * cms2.trks_trk_p4()[trkidx].pz()/cms2.trks_trk_p4()[trkidx].pt());
  return trks_dz_pv(trkidx,vtxidx).first;

}

float pfjet_beta(int ijet, int power , float dzcut , int ivtx , bool verbose ) {

  //---------------------------------------------------------------------------------
  // cands is the vector of PFCandidate indices matched to the pfjet with index ijet
  //---------------------------------------------------------------------------------

  vector<int> cands = cms2.pfjets_pfcandIndicies().at(ijet);

  float pt_tot = 0.0;
  float pt_vtx = 0.0;

  //-------------------------------------
  // loop over PFCandidates in jet
  //-------------------------------------
  
  for (unsigned int ivc=0;ivc<cands.size();ivc++) {

    int ican = cands[ivc];

    if( ican > static_cast<int>(cms2.pfcands_charge().size()) ){
      cout << __FILE__ << " " << __LINE__ << " ERROR! ican, numCandidates " << ican << ", " << cms2.pfcands_charge().size() << endl;
      return -3.0;
    }

    //--------------------
    // skip neutrals 
    //--------------------
    
    if (cms2.pfcands_charge().at(ican)==0) continue;

    //-------------------------------------------------------------
    // get track from PFCandidate, this is used to calculate dz
    //-------------------------------------------------------------
    
    int itrk = cms2.pfcands_trkidx().at(ican);
    
    if( itrk >= (int) cms2.trks_trk_p4().size() || itrk < 0 ){
      if( verbose ){
	std::cout << __FILE__ << " " << __LINE__ << " WARNING! skipping electron with pt " << cms2.pfcands_p4().at(ican).pt() << endl;
      }
      //note: this should only happen for electrons which do not have a matched track
      //currently we are just ignoring these guys
      continue;
    }
    
    //-------------------------------------------------------------------
    // calculate scalar sum of  pT(track)^(power) for all tracks in jet
    //-------------------------------------------------------------------

    pt_tot += pow( cms2.pfcands_p4().at(ican).pt() , power );

    //-------------------------------------------------------------------    
    // require dz(trk,vtx) < dzcut
    //-------------------------------------------------------------------

    float dz = trks_dz_pv(itrk,ivtx).first;  
    if ( fabs(dz) > dzcut )     continue;

    //-------------------------------------------------------------------
    // calculate scalar sum of  pT(track)^(power) for all tracks in jet
    // which are within dz < dzcut of the given vtx
    //-------------------------------------------------------------------
    
    pt_vtx += pow( cms2.pfcands_p4().at(ican).pt() , power );
    
  }

  //-------------------------------------------------------------------    
  // calculate beta, the ratio of pt_vtx to pt_tot
  //-------------------------------------------------------------------

  float beta = -2.0;
  if( pt_tot > 1e-5 ) beta = pt_vtx / pt_tot;
  return beta;

}

vector<LorentzVector> getBtaggedJets (unsigned int i_hyp, bool sort_, enum JetType type, enum CleaningType cleaning,
                                      enum BtagType btag_type, double deltaR, double min_pt, double max_eta, 
                                      double rescale,  int systFlag)
{
  jets_with_corr_t jets = getJets_fast(i_hyp, type, cleaning, deltaR, min_pt, max_eta, rescale, systFlag);
    const vector<float> *btags = 0;
    const float btag_wp = BtagWP[btag_type];

    switch (btag_type) {
    case JETS_BTAG_NONE:
        btags = 0;
    case JETS_BTAG_TCHEL: case JETS_BTAG_TCHEM:
        switch (type) {
        case JETS_TYPE_CALO_CORR: case JETS_TYPE_CALO_UNCORR:
            btags = &cms2.jets_trackCountingHighEffBJetTag();
            break;
        case JETS_TYPE_PF_UNCORR: case JETS_TYPE_PF_CORR: case JETS_TYPE_PF_FAST_CORR: case JETS_TYPE_PF_FAST_CORR_RESIDUAL:
            btags = &cms2.pfjets_trackCountingHighEffBJetTag();
            break;
        default:
            assert(false);
        }
        break;
    case JETS_BTAG_TCHPM: case JETS_BTAG_TCHPT:
        switch (type) {
        case JETS_TYPE_CALO_CORR: case JETS_TYPE_CALO_UNCORR:
            btags = &cms2.jets_trackCountingHighPurBJetTag();
            break;
        case JETS_TYPE_PF_UNCORR: case JETS_TYPE_PF_CORR: case JETS_TYPE_PF_FAST_CORR: case JETS_TYPE_PF_FAST_CORR_RESIDUAL:
            btags = &cms2.pfjets_trackCountingHighPurBJetTag();
            break;
        default:
            assert(false);
        }
        break;
    case JETS_BTAG_SSVHEM:
        switch (type) {
        case JETS_TYPE_CALO_CORR: case JETS_TYPE_CALO_UNCORR:
            btags = &cms2.jets_simpleSecondaryVertexHighEffBJetTag();
            break;
        case JETS_TYPE_PF_UNCORR: case JETS_TYPE_PF_CORR: case JETS_TYPE_PF_FAST_CORR: case JETS_TYPE_PF_FAST_CORR_RESIDUAL:
            btags = &cms2.pfjets_simpleSecondaryVertexHighEffBJetTag();
            break;
        default:
            assert(false);
        }
        break;
    case JETS_BTAG_SSVHPT:
        switch (type) {
        case JETS_TYPE_CALO_CORR: case JETS_TYPE_CALO_UNCORR:
            btags = &cms2.jets_simpleSecondaryVertexHighPurBJetTags();
            break;
        case JETS_TYPE_PF_UNCORR: case JETS_TYPE_PF_CORR: case JETS_TYPE_PF_FAST_CORR: case JETS_TYPE_PF_FAST_CORR_RESIDUAL:
            btags = &cms2.pfjets_simpleSecondaryVertexHighPurBJetTags();
            break;
        default:
            assert(false);
        }
        break;
    case JETS_BTAG_CSVL: case JETS_BTAG_CSVM: case JETS_BTAG_CSVT:
        switch (type) {
        case JETS_TYPE_CALO_CORR: case JETS_TYPE_CALO_UNCORR:
            btags = &cms2.jets_combinedSecondaryVertexBJetTag();
            break;
        case JETS_TYPE_PF_UNCORR: case JETS_TYPE_PF_CORR: case JETS_TYPE_PF_FAST_CORR: case JETS_TYPE_PF_FAST_CORR_RESIDUAL:
            btags = &cms2.pfjets_combinedSecondaryVertexBJetTag();
            break;
        default:
            assert(false);
        }
        break;
    default:
        assert(false);
    }

    vector<LorentzVector> ret;
    ret.reserve(jets.size());
    if (jets.size() == 0 || btags == 0 )
        return ret;

    for (unsigned int i = 0; i < jets.size(); ++i) {
        if (jets[i].first != 0 && btags->at(i) > btag_wp )
            ret.push_back(*jets[i].first * jets[i].second);         
    }

    if (sort_)
        sort(ret.begin(), ret.end(), jets_pt_gt());
    return ret;
}

std::vector<bool> getBtaggedJetFlags (unsigned int i_hyp, enum JetType type, enum CleaningType cleaning,
                                      enum BtagType btag_type, double deltaR, double min_pt, double max_eta, 
                                      double rescale,  int systFlag)
{
  jets_with_corr_t jets = getJets_fast(i_hyp, type, cleaning, deltaR, min_pt, max_eta, rescale, systFlag);
    const vector<float> *btags = 0;
    const float btag_wp = BtagWP[btag_type];

    switch (btag_type) {
    case JETS_BTAG_NONE:
        btags = 0;
    case JETS_BTAG_TCHEL: case JETS_BTAG_TCHEM:
        switch (type) {
        case JETS_TYPE_CALO_CORR: case JETS_TYPE_CALO_UNCORR:
            btags = &cms2.jets_trackCountingHighEffBJetTag();
            break;
        case JETS_TYPE_PF_UNCORR: case JETS_TYPE_PF_CORR: case JETS_TYPE_PF_FAST_CORR: case JETS_TYPE_PF_FAST_CORR_RESIDUAL:
            btags = &cms2.pfjets_trackCountingHighEffBJetTag();
            break;
        default:
            assert(false);
        }
        break;
    case JETS_BTAG_TCHPM: case JETS_BTAG_TCHPT:
        switch (type) {
        case JETS_TYPE_CALO_CORR: case JETS_TYPE_CALO_UNCORR:
            btags = &cms2.jets_trackCountingHighPurBJetTag();
            break;
        case JETS_TYPE_PF_UNCORR: case JETS_TYPE_PF_CORR: case JETS_TYPE_PF_FAST_CORR: case JETS_TYPE_PF_FAST_CORR_RESIDUAL:
            btags = &cms2.pfjets_trackCountingHighPurBJetTag();
            break;
        default:
            assert(false);
        }
        break;
    case JETS_BTAG_SSVHEM:
        switch (type) {
        case JETS_TYPE_CALO_CORR: case JETS_TYPE_CALO_UNCORR:
            btags = &cms2.jets_simpleSecondaryVertexHighEffBJetTag();
            break;
        case JETS_TYPE_PF_UNCORR: case JETS_TYPE_PF_CORR: case JETS_TYPE_PF_FAST_CORR: case JETS_TYPE_PF_FAST_CORR_RESIDUAL:
            btags = &cms2.pfjets_simpleSecondaryVertexHighEffBJetTag();
            break;
        default:
            assert(false);
        }
        break;
    case JETS_BTAG_SSVHPT:
        switch (type) {
        case JETS_TYPE_CALO_CORR: case JETS_TYPE_CALO_UNCORR:
            btags = &cms2.jets_simpleSecondaryVertexHighPurBJetTags();
            break;
        case JETS_TYPE_PF_UNCORR: case JETS_TYPE_PF_CORR: case JETS_TYPE_PF_FAST_CORR: case JETS_TYPE_PF_FAST_CORR_RESIDUAL:
            btags = &cms2.pfjets_simpleSecondaryVertexHighPurBJetTags();
            break;
        default:
            assert(false);
        }
        break;
    case JETS_BTAG_CSVL: case JETS_BTAG_CSVM: case JETS_BTAG_CSVT:
        switch (type) {
        case JETS_TYPE_CALO_CORR: case JETS_TYPE_CALO_UNCORR:
            btags = &cms2.jets_combinedSecondaryVertexBJetTag();
            break;
        case JETS_TYPE_PF_UNCORR: case JETS_TYPE_PF_CORR: case JETS_TYPE_PF_FAST_CORR: case JETS_TYPE_PF_FAST_CORR_RESIDUAL:
            btags = &cms2.pfjets_combinedSecondaryVertexBJetTag();
            break;
        default:
            assert(false);
        }
        break;
    default:
        assert(false);
    }

    vector<bool> ret;
    ret.reserve(jets.size());
    if (jets.size() == 0 || btags == 0 )
        return ret;

    for (unsigned int i = 0; i < jets.size(); ++i) {
        ret.push_back(jets[i].first != 0 && btags->at(i) > btag_wp );
    }

    return ret;
}

int nBtaggedJets (unsigned int i_hyp, enum JetType type, enum CleaningType cleaning,
                  enum BtagType btag_type, double deltaR, double min_pt, double max_eta, double rescale,  int systFlag)
{
  std::vector<LorentzVector> jets = getBtaggedJets(i_hyp, false, type, cleaning, btag_type, deltaR, min_pt, max_eta, rescale, systFlag);
    return jets.size();
}
