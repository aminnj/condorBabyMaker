#ifndef SSSELECTIONS_H
#define SSSELECTIONS_H

#include <vector>
#include <string>
#include "Math/LorentzVector.h"
#include "CMS2.h"
#include "jetSelections.h"
#include "electronSelections.h"
#include "jetcorr/JetCorrectionUncertainty.h"

/////////////////////////////////////////////////////////////////
///                                                           ///
///                                                           ///
///                                                           ///
///          2012 Selections                                  ///
///                                                           ///
///                                                           ///
///                                                           ///
/////////////////////////////////////////////////////////////////

namespace samesign 
{
    ////////////////////////////////////////////////////////////////////////////////////////////     
    // 2012 good lepton (passes ID)
    ////////////////////////////////////////////////////////////////////////////////////////////     
    bool isGoodLepton(int id, int idx, bool use_el_eta = true);


    ////////////////////////////////////////////////////////////////////////////////////////////     
    // 2012 isolated lepton
    ////////////////////////////////////////////////////////////////////////////////////////////     
    bool isIsolatedLepton(int id, int idx);


    ////////////////////////////////////////////////////////////////////////////////////////////     
    // 2012 lepton isolation value
    ////////////////////////////////////////////////////////////////////////////////////////////     
    double leptonIsolation(int id, int idx);


    ////////////////////////////////////////////////////////////////////////////////////////////     
    // 2012 effective area 
    ////////////////////////////////////////////////////////////////////////////////////////////     
    float EffectiveArea03(int id, int idx);
    float EffectiveArea04(int id, int idx);


    ////////////////////////////////////////////////////////////////////////////////////////////     
    // 2012 numerator lepton (passes ID and isolation)
    ////////////////////////////////////////////////////////////////////////////////////////////     
    bool isNumeratorLepton(int id, int idx, bool use_el_eta = true);


    ////////////////////////////////////////////////////////////////////////////////////////////     
    // 2012 numerator hypothesis (passes ID and isolation)
    ////////////////////////////////////////////////////////////////////////////////////////////     
    bool isNumeratorHypothesis(int idx, bool use_el_eta = true);


    ////////////////////////////////////////////////////////////////////////////////////////////     
    // 2012 denominator lepton (relaxed ID and Isolation)
    ////////////////////////////////////////////////////////////////////////////////////////////     
    bool isDenominatorLepton(int id, int idx, bool use_el_eta = true);


    ////////////////////////////////////////////////////////////////////////////////////////////     
    // 2012 denominator hypothesis (relaxed ID and Isolation)
    ////////////////////////////////////////////////////////////////////////////////////////////     
    bool isDenominatorHypothesis(int idx, bool use_el_eta = true);


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 require electron GSF, CTF and SC charges agree
    ///////////////////////////////////////////////////////////////////////////////////////////
    bool passThreeChargeRequirement(int elIdx);


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 calculate PF-based isolation for electrons with rho*Aeff correction
    ///////////////////////////////////////////////////////////////////////////////////////////
    float electronIsolationPF2012_cone03(int idx);  // uses ∆R < 0.3
    float electronIsolationPF2012_cone04(int idx);  // uses ∆R < 0.4
    float electronIsolationPF2012(int idx);         // wrapper ∆R < 0.3 version which is used in the analysis


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 passes dilepton trigger
    ///////////////////////////////////////////////////////////////////////////////////////////

    // analysis type:
    //   0 --> use high pT analysis triggers
    //   1 --> use low pT analysis triggers
    //   2 --> use very low pT analysis triggers
    //   anything else will return false

    bool passesTrigger(int hyp_type, int analysis_type);
    bool passesTriggerHighPt(int hyp_type);
    bool passesTriggerLowPt (int hyp_type);
    bool passesTriggerVeryLowPt(int hyp_type);


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 extra Z veto for b-tagged same sign analysis
    ///////////////////////////////////////////////////////////////////////////////////////////
    bool makesExtraZ(int idx, bool apply_id_iso = true); 


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 extra gamma* veto for b-tagged same sign analysis
    ///////////////////////////////////////////////////////////////////////////////////////////
    bool makesExtraGammaStar(int idx, bool apply_id_iso = true); 


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 get jets and perform overlap removal with numerator e/mu with pt > x (defaults are 20/20 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////

    // JEC taken from ntuple
    std::vector<LorentzVector> getJets(int idx, enum JetType type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0, float rescale = 1.0, int systFlag = 0);

    // JEC applied otf
    std::vector<LorentzVector> getJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0, float rescale = 1.0, int systFlag = 0);

    // JEC uncertainty applied otf
    std::vector<LorentzVector> getJets(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0);	 

    // JEC AND JEC uncertainty applied otf
    std::vector<LorentzVector> getJets(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type,  enum JetType type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0);	 


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 get jet flags and perform overlap removal with numerator e/mu with pt > x (defaults are 20/20 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////

    // JEC taken from ntuple
    std::vector<bool> getJetFlags(int idx, enum JetType type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0, float rescale = 1.0, int systFlag = 0);

    // JEC applied otf
    std::vector<bool> getJetFlags(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0, float rescale = 1.0, int systFlag = 0);

    // JEC uncertainty applied otf
    std::vector<bool> getJetFlags(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0);	 

    // JEC AND JEC uncertainty applied otf
    std::vector<bool> getJetFlags(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type,  enum JetType type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0);	 


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 get sumpt, skip jets overlapping with numerator e/mu with pt>x (defaults are 20/20 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////

    // JEC taken from ntuple
    float sumJetPt(int idx, enum JetType type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0, float rescale = 1.0, int systFlag = 0);

    // JEC applied otf
    float sumJetPt(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0, float rescale = 1.0, int systFlag = 0);

    // JEC uncertainty applied otf
    float sumJetPt(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0);	 

    // JEC AND JEC uncertainty applied otf
    float sumJetPt(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0);	 


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 get njets, skip jets overlapping with numerator e/mu with pt>x (defaults are 20/20 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////

    // JEC taken from ntuple
    int nJets(int idx, enum JetType type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0, float rescale = 1.0, int systFlag = 0);

    // JEC applied otf
    int nJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0, float rescale = 1.0, int systFlag = 0);

    // JEC uncertainty applied otf
    int nJets(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0);	 

    // JEC AND JEC uncertainty applied otf
    int nJets(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0);	 


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 get b-tagged jets and perform overlap removal with numerator e/mu with pt > x (defaults are 20/20 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////

    // JEC taken from ntuple
    std::vector<LorentzVector> getBtaggedJets(int idx, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0, float rescale = 1.0, int systFlag = 0);

    // JEC applied otf
    std::vector<LorentzVector> getBtaggedJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0, float rescale = 1.0, int systFlag = 0);

    // JEC uncertainty applied otf
    std::vector<LorentzVector> getBtaggedJets(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0);	 

    // JEC AND JEC uncertainty applied otf
    std::vector<LorentzVector> getBtaggedJets(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0); 


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 get b-tagged jet flags and perform overlap removal with numerator e/mu with pt > x (defaults are 20/20 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////

    // JEC taken from ntuple
    std::vector<bool> getBtaggedJetFlags(int idx, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0, float rescale = 1.0, int systFlag = 0);

    // JEC applied otf
    std::vector<bool> getBtaggedJetFlags(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0, float rescale = 1.0, int systFlag = 0);

    // JEC uncertainty applied otf
    std::vector<bool> getBtaggedJetFlags(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0);	 

    // JEC AND JEC uncertainty applied otf
    std::vector<bool> getBtaggedJetFlags(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0); 


    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 get sumpt, skip jets overlapping with numerator e/mu with pt>x (defaults are 20/20 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////

    // JEC taken from ntuple
    int nBtaggedJets(int idx, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0, float rescale = 1.0, int systFlag = 0);

    // JEC applied otf
    int nBtaggedJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0, float rescale = 1.0, int systFlag = 0);

    // JEC uncertainty applied otf
    int nBtaggedJets(int idx, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0);	 

    // JEC AND JEC uncertainty applied otf
    int nBtaggedJets(int idx, FactorizedJetCorrector* jet_corrector, JetCorrectionUncertainty *jet_unc, enum JetScaleType scale_type, enum JetType type, enum BtagType btag_type, float deltaR = 0.4, float min_pt = 40.0, float max_eta = 2.4, float mu_minpt = 20.0, float ele_minpt = 20.0);	 

    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 rescale the jet energy resolution (JER) 
    ///////////////////////////////////////////////////////////////////////////////////////////
    void smearJETScaleJetsMetHt(std::vector<LorentzVector>& vjets_p4, float& met, float& met_phi, float& ht, const unsigned int seed);
    void smearJETScaleJets(std::vector<LorentzVector>& vjets_p4, const unsigned int seed);
    void smearJETScaleJetsMetHt
    (
        std::vector<LorentzVector>& vjets_p4, 
        float& met,
        float& met_phi,
        float& ht, 
        int idx,
        enum JetType type,
        const unsigned int seed,
        float deltaR = 0.4,
        float min_pt = 40.0,
        float max_eta = 2.4,
        float mu_minpt = 20.0,
        float ele_minpt = 20.0
    );

    ///////////////////////////////////////////////////////////////////////////////////////////
    // 2012 rescale the MET by scaling up/down the unclustered erngy 
    ///////////////////////////////////////////////////////////////////////////////////////////
    float scaleMET
    (
        const float met,
        const float met_phi,
        int idx,
        enum JetType type,
        float deltaR = 0.4,
        float min_pt = 40.0,
        float max_eta = 2.4,
        float mu_minpt = 20.0,
        float ele_minpt = 20.0,
        const int scale_type = 0,
        const float scale = 0.1
    );

    ///////////////////////////////////////////////////////////////////////////////////////////	 
    // 2012 get vector of good els p4s	 
    ///////////////////////////////////////////////////////////////////////////////////////////	 
    std::vector<LorentzVector> getGoodElectrons(const float ptcut = 20.0f);	 
    std::vector<std::pair<LorentzVector, unsigned int> > getNumeratorElectrons(const float ptcut = 20.0f);	 


    ///////////////////////////////////////////////////////////////////////////////////////////	 
    // 2012 get vector of good mus p4s	 
    ///////////////////////////////////////////////////////////////////////////////////////////	 
    std::vector<LorentzVector> getGoodMuons(const float ptcut = 20.0f);
    std::vector<std::pair<LorentzVector, unsigned int> > getNumeratorMuons(const float ptcut = 20.0f);	 

} // namespace samesign


/////////////////////////////////////////////////////////////////
///                                                           ///
///                                                           ///
///                                                           ///
///          2011 Selections                                  ///
///                                                           ///
///                                                           ///
///                                                           ///
/////////////////////////////////////////////////////////////////

namespace samesign2011 
{

    enum IsolationType { DET_ISO, TIGHT_DET_ISO, COR_DET_ISO };

    ////////////////////////////////////////////////////////////////////////////////////////////     
    // 2011 good lepton
    ////////////////////////////////////////////////////////////////////////////////////////////     
    bool isGoodLepton(int id, int idx);


    ////////////////////////////////////////////////////////////////////////////////////////////     
    // 2011 isolated lepton
    ////////////////////////////////////////////////////////////////////////////////////////////     
    bool isIsolatedLepton(int id, int idx, enum IsolationType iso_type = DET_ISO);


    ////////////////////////////////////////////////////////////////////////////////////////////     
    // 2011 numerator lepton
    ////////////////////////////////////////////////////////////////////////////////////////////     
    bool isNumeratorLepton(int id, int idx, enum IsolationType iso_type = DET_ISO);


    ////////////////////////////////////////////////////////////////////////////////////////////     
    // 2011 numerator hypothesis
    ////////////////////////////////////////////////////////////////////////////////////////////     
    bool isNumeratorHypothesis(int idx, enum IsolationType iso_type = DET_ISO);


    ////////////////////////////////////////////////////////////////////////////////////////////     
    // 2011 denominator lepton
    ////////////////////////////////////////////////////////////////////////////////////////////     
    bool isDenominatorLepton(int id, int idx, enum IsolationType iso_type = DET_ISO);


    ///////////////////////////////////////////////////////////////////////////////////////////
    // require electron GSF, CTF and SC charges agree
    ///////////////////////////////////////////////////////////////////////////////////////////
    bool passThreeChargeRequirement(int elIdx);


    ///////////////////////////////////////////////////////////////////////////////////////////
    // get jets and perform overlap removal with numerator e/mu with pt > x (defaults are 10/5 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////
    std::vector<LorentzVector> getJets(int idx, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt = 5, double ele_minpt = 10, enum IsolationType iso_type = DET_ISO, double rescale = 1.0);


    ///////////////////////////////////////////////////////////////////////////////////////////
    // get jets and apply an on-the-fly JEC and perform overlap removal with numerator
    // e/mu with pt > x (defaults are 10/5 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////
    std::vector<LorentzVector> getJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt = 5, double ele_minpt = 10, enum IsolationType iso_type = DET_ISO, double rescale = 1.0);


    ///////////////////////////////////////////////////////////////////////////////////////////
    // get jets and perform overlap removal with numerator e/mu with pt > x (defaults are 10/5 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////
    std::vector<bool> getJetFlags(int idx, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt = 5, double ele_minpt = 10, enum IsolationType iso_type = DET_ISO, double rescale = 1.0);


    ///////////////////////////////////////////////////////////////////////////////////////////
    // get jets and apply an on-the-fly JEC and perform overlap removal with numerator
    // e/mu with pt > x (defaults are 10/5 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////
    std::vector<bool> getJetFlags(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt = 5, double ele_minpt = 10, enum IsolationType iso_type = DET_ISO, double rescale = 1.0);


    ///////////////////////////////////////////////////////////////////////////////////////////
    // get sumpt, skip jets overlapping with numerator e/mu with pt>x (defaults are 10/5 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////
    float sumJetPt(int idx, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt = 5, double ele_minpt = 10, enum IsolationType iso_type = DET_ISO, double rescale = 1.0);


    ///////////////////////////////////////////////////////////////////////////////////////////
    // same as above, but allowing use of on-the-fly JEC corrections
    ///////////////////////////////////////////////////////////////////////////////////////////
    float sumJetPt(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt = 5, double ele_minpt = 10, enum IsolationType iso_type = DET_ISO, double rescale = 1.0);


    ///////////////////////////////////////////////////////////////////////////////////////////
    // get sumpt, skip jets overlapping with numerator e/mu with pt>x (defaults are 10/5 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////
    int nJets(int idx, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt = 5, double ele_minpt = 10, enum IsolationType iso_type = DET_ISO, double rescale = 1.0);


    ///////////////////////////////////////////////////////////////////////////////////////////
    // same as above, but allowing use of on-the-fly JEC corrections
    ///////////////////////////////////////////////////////////////////////////////////////////
    int nJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, double deltaR, double min_pt, double max_eta, double mu_minpt = 5, double ele_minpt = 10, enum IsolationType iso_type = DET_ISO, double rescale = 1.0);


    ///////////////////////////////////////////////////////////////////////////////////////////
    // extra Z veto for generic same sign analysis
    ///////////////////////////////////////////////////////////////////////////////////////////
    bool overlapsOtherNNHypInZ(int idx, enum IsolationType iso_type = DET_ISO); //similar to makesExtraZ, uses hyps only


    ///////////////////////////////////////////////////////////////////////////////////////////
    // extra Z veto for b-tagged same sign analysis
    ///////////////////////////////////////////////////////////////////////////////////////////
    bool makesExtraZ(int idx, enum IsolationType iso_type = COR_DET_ISO, bool apply_id_iso = false); //similar to makesExtraZ, uses hyps only


    ///////////////////////////////////////////////////////////////////////////////////////////
    // passes dilepton trigger
    ///////////////////////////////////////////////////////////////////////////////////////////
    bool passesTrigger(bool is_data, int hyp_type, bool is_high_pt);


    ///////////////////////////////////////////////////////////////////////////////////////////
    // get jets and perform overlap removal with numerator e/mu with pt > x (defaults are 10/5 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////
    std::vector<LorentzVector> getBtaggedJets(int idx, enum JetType type, enum BtagType btag_type, double deltaR, double min_pt, double max_eta, double mu_minpt = 5, double ele_minpt = 10, enum IsolationType iso_type = DET_ISO, double rescale = 1.0);

    ///////////////////////////////////////////////////////////////////////////////////////////
    // get jets and apply an on-the-fly JEC and perform overlap removal with numerator
    // e/mu with pt > x (defaults are 10/5 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////
    std::vector<LorentzVector> getBtaggedJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, double deltaR, double min_pt, double max_eta, double mu_minpt = 5, double ele_minpt = 10, enum IsolationType iso_type = DET_ISO, double rescale = 1.0);


    ///////////////////////////////////////////////////////////////////////////////////////////
    // get sumpt, skip jets overlapping with numerator e/mu with pt>x (defaults are 10/5 GeV)
    ///////////////////////////////////////////////////////////////////////////////////////////
    int nBtaggedJets(int idx, enum JetType type, enum BtagType btag_type, double deltaR, double min_pt, double max_eta, double mu_minpt = 5, double ele_minpt = 10, enum IsolationType iso_type = DET_ISO, double rescale = 1.0);


    ///////////////////////////////////////////////////////////////////////////////////////////
    // same as above, but allowing use of on-the-fly JEC corrections
    ///////////////////////////////////////////////////////////////////////////////////////////
    int nBtaggedJets(int idx, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, double deltaR, double min_pt, double max_eta, double mu_minpt = 5, double ele_minpt = 10, enum IsolationType iso_type = DET_ISO, double rescale = 1.0);

} // namespace 2011

#endif

