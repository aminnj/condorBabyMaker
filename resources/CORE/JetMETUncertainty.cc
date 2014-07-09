#include "JetMETUncertainty.h"
#include "jetcorr/JetCorrectionUncertainty.h"
#include "jetSmearingTools.h"
#include <math.h>

using namespace std;

//
// default constructor
//
JetMETUncertainty::JetMETUncertainty ()
    : ele_unc_b_(0.006), ele_unc_e_(0.015), mu_unc_(0.01), uncl_unc_(0.1), have_jetcorr_unc_(false), have_jetsmear_(false)
{
    jetCorrectionUncertainty_ = 0;
    jetSmearer_               = 0;
}

//
// constructor taking jet parameter file names
//
JetMETUncertainty::JetMETUncertainty(string jetcorr_unc_file_name,
                                     vector<string> &jetsmear_file_names)
    : ele_unc_b_(0.006), ele_unc_e_(0.015), mu_unc_(0.01), uncl_unc_(0.1)
{
    have_jetcorr_unc_ = false;
    jetcorr_unc_file_name_ = jetcorr_unc_file_name;
    if (!jetcorr_unc_file_name.empty()) {
        jetCorrectionUncertainty_ = new JetCorrectionUncertainty(jetcorr_unc_file_name_);
        if (jetCorrectionUncertainty_ != 0)
            have_jetcorr_unc_ = true;
        else
            cout << "Trouble getting jet correction uncertainty service." << endl;
    }

    have_jetsmear_ = false;
    jetsmear_file_names_ = vector<string>(jetsmear_file_names);
    if (jetsmear_file_names_.size() == 3) {
        jetSmearer_ = makeJetSmearer(jetsmear_file_names_);
        if (jetSmearer_ != 0)
            have_jetsmear_ = true;
        else
            cout << "Trouble getting jet resolution smearing service." << endl;
    }
}

//
// return total met uncertainty
// return a pair<double, double>
// first entry is the uncertainty on the UP scaling
// second entry is the uncertainty on the DOWN scaling
//
dpair JetMETUncertainty::GetTotalUncertainty(dpair good_met)
{
    good_met_ = dpair(good_met);
    double base_met = good_met_.first;

    vector<dpair> vmet_up;
    vector<dpair> vmet_down;

    SmearMETForElectronUncertainty ();
    vmet_up.push_back(met_ele_up_);
    vmet_down.push_back(met_ele_down_);

    SmearMETForMuonUncertainty ();
    vmet_up.push_back(met_mu_up_);
    vmet_down.push_back(met_mu_down_);

    SmearMETForUnclusteredEnergyUncertainty ();
    vmet_up.push_back(met_uncl_up_);
    vmet_down.push_back(met_uncl_down_);

    if (have_jetcorr_unc_) {
        SmearMETForJESUncertainty ();
        vmet_up.push_back(met_jes_up_);
        vmet_down.push_back(met_jes_down_);
    }
    
    if (have_jetsmear_) {
        SmearMETForJERUncertainty ();
        vmet_up.push_back(met_jer_up_);
        vmet_down.push_back(met_jer_down_);
    }

    assert(vmet_up.size() == vmet_down.size());

    double unc_up   = 0.;
    double unc_down = 0.;
    for (unsigned int idx = 0; idx < vmet_up.size(); idx++) {
        double tmp_up   = (vmet_up.at(idx).first - base_met) / base_met;
        unc_up += pow(tmp_up, 2);

        double tmp_down = (vmet_down.at(idx).first - base_met) / base_met;
        unc_down += pow(tmp_down, 2);
    }

    return make_pair(sqrt(unc_up), sqrt(unc_down));
}

//
// get scaled MET
// component is after scaling either electrons (ELE), muons (MU), unclustered energy (UNCL)
// jet energy scale (JES), or jet energy resolution (JER)
// scale_type is either UP or DOWN
// returns a pair<double, double>
// first entry is the MET
// second entry is the MET phi
//
dpair JetMETUncertainty::GetScaledMET (dpair good_met, Component component, ScaleType scale_type)
{
    good_met_ = dpair(good_met);

    if (component == ELE) {
        SmearMETForElectronUncertainty ();
        
        return (scale_type == DOWN) ? met_ele_down_ : met_ele_up_;
    }
    else if (component == MU) {
        SmearMETForMuonUncertainty ();
        
        return (scale_type == DOWN) ? met_mu_down_ : met_mu_up_;
    }
    else if (component == UNCL) {
        SmearMETForUnclusteredEnergyUncertainty ();
        
        return (scale_type == DOWN) ? met_uncl_down_ : met_uncl_up_;
    }
    else if (component == JES) {
        SmearMETForJESUncertainty ();

        return (scale_type == DOWN) ? met_jes_down_ : met_jes_up_;
    }
    else if (component == JER) {
        SmearMETForJERUncertainty ();

        return (scale_type == DOWN) ? met_jer_down_ : met_jer_up_;
    }
    else {
        cout << "Do not recognize the component you asked for." << endl;
        return make_pair(-1., 0.);
    }
}

//
// smear met for electron energy scale uncertainty
// default uncertainty is 0.6% in barrel
// and 1.5% in endcap
//
void JetMETUncertainty::SmearMETForElectronUncertainty ()
{
    double dmet_x = good_met_.first * cos(good_met_.second);
    double dmet_y = good_met_.first * sin(good_met_.second);
    double umet_x = good_met_.first * cos(good_met_.second);
    double umet_y = good_met_.first * sin(good_met_.second);
    for (unsigned int idx = 0; idx < good_els_.size(); idx++) {
        double eta   = fabs(good_els_.at(idx).eta());
        double px    = good_els_.at(idx).px();
        double py    = good_els_.at(idx).py();
        double scale = (eta < 1.479) ? ele_unc_b_ : ele_unc_e_;

        dmet_x -= DOWN * scale * px;
        dmet_y -= DOWN * scale * py;

        umet_x -= UP * scale * px;
        umet_y -= UP * scale * py;
    }

    met_ele_up_   = make_pair(sqrt( pow(umet_x, 2) + pow(umet_y, 2)), atan2(umet_y, umet_x));
    met_ele_down_ = make_pair(sqrt( pow(dmet_x, 2) + pow(dmet_y, 2)), atan2(dmet_y, dmet_x));
}

//
// smear met for muon momentum scale uncertainty
// default uncertainty is 1.0%
//
void JetMETUncertainty::SmearMETForMuonUncertainty ()
{
    double dmet_x = good_met_.first * cos(good_met_.second);
    double dmet_y = good_met_.first * sin(good_met_.second);
    double umet_x = good_met_.first * cos(good_met_.second);
    double umet_y = good_met_.first * sin(good_met_.second);
    for (unsigned int idx = 0; idx < good_mus_.size(); idx++) {
        double px    = good_mus_.at(idx).px();
        double py    = good_mus_.at(idx).py();        

        dmet_x -= DOWN * mu_unc_ * px;
        dmet_y -= DOWN * mu_unc_ * py;

        umet_x -= UP * mu_unc_ * px;
        umet_y -= UP * mu_unc_ * py;
    }

    met_mu_up_   = make_pair(sqrt( pow(umet_x, 2) + pow(umet_y, 2)), atan2(umet_y, umet_x));
    met_mu_down_ = make_pair(sqrt( pow(dmet_x, 2) + pow(dmet_y, 2)), atan2(dmet_y, dmet_x));
}

//
// smear met for unclustered energy scale uncertainty
// default uncertainty is 10%
//
void JetMETUncertainty::SmearMETForUnclusteredEnergyUncertainty ()
{  
    dpair tmp_met = GetUnclusteredMET ();

    double dmet_x = good_met_.first * cos(good_met_.second);
    double dmet_y = good_met_.first * sin(good_met_.second);

    dmet_x -= DOWN * uncl_unc_ * tmp_met.first;
    dmet_y -= DOWN * uncl_unc_ * tmp_met.second;

    double umet_x = good_met_.first * cos(good_met_.second);
    double umet_y = good_met_.first * sin(good_met_.second);

    umet_x -= DOWN * uncl_unc_ * tmp_met.first;
    umet_y -= DOWN * uncl_unc_ * tmp_met.second;

    met_uncl_up_   = make_pair(sqrt( pow(umet_x, 2) + pow(umet_y, 2)), atan2(umet_y, umet_x));
    met_uncl_down_ = make_pair(sqrt( pow(dmet_x, 2) + pow(dmet_y, 2)), atan2(dmet_y, dmet_x));
}

//
// calculate component of MET from unclustered energy
// start with input MET and remove the following
// 1) good electrons
// 2) good muons
// 3) good jets
//
dpair JetMETUncertainty::GetUnclusteredMET ()
{
    double met_x = good_met_.first * cos(good_met_.second);
    double met_y = good_met_.first * sin(good_met_.second);
    for (unsigned int idx = 0; idx < good_els_.size(); idx++) {
        met_x += good_els_.at(idx).px();
        met_y += good_els_.at(idx).py();
    }

    for (unsigned int idx = 0; idx < good_mus_.size(); idx++) {
        met_x += good_mus_.at(idx).px();
        met_y += good_mus_.at(idx).py();
    }

    for (unsigned int idx = 0; idx < good_jets_.size(); idx++) {
        met_x += good_jets_.at(idx).px();
        met_y += good_jets_.at(idx).py();
    }

    return make_pair(met_x, met_y);
}

//
// smear met for JES uncertainty
// uncertainty determined on jet-by-jet basis
// using text file from jetcorr/data
// the good_jets input to this method are the final corrected jets used by the analysis
// after all selections on the jets
// the default selections should be
// |eta| < 4.7 and corrected pt > 10 GeV
//
void JetMETUncertainty::SmearMETForJESUncertainty ()
{
    if (!have_jetcorr_unc_) {
        cout << "Attempting to smear MET for JES but do not have a valid jet correction uncertainty service." << endl;
        met_jes_up_   = make_pair(-1., 0.);
        met_jes_down_ = make_pair(-1., 0.);
        return;
    }

    double dmet_x = good_met_.first * cos(good_met_.second);
    double dmet_y = good_met_.first * sin(good_met_.second);
    double umet_x = good_met_.first * cos(good_met_.second);
    double umet_y = good_met_.first * sin(good_met_.second);

    jetcorr_uncertainties_.clear();
    for (unsigned int idx = 0; idx < good_jets_.size(); idx++) {
        double px = good_jets_.at(idx).px();
        double py = good_jets_.at(idx).py();

        jetCorrectionUncertainty_->setJetEta(good_jets_.at(idx).eta());
        jetCorrectionUncertainty_->setJetPt(good_jets_.at(idx).pt());
        double scale = jetCorrectionUncertainty_->getUncertainty(true);
        jetcorr_uncertainties_.push_back(scale);

        dmet_x -= DOWN * scale * px;
        dmet_y -= DOWN * scale * py;

        umet_x -= UP * scale * px;
        umet_y -= UP * scale * py;
    }

    met_jes_up_   = make_pair(sqrt( pow(umet_x, 2) + pow(umet_y, 2)), atan2(umet_y, umet_x));
    met_jes_down_ = make_pair(sqrt( pow(dmet_x, 2) + pow(dmet_y, 2)), atan2(dmet_y, dmet_x));
}

//
// smear met for JER uncertainty
// uncertainty determined on jet-by-jet basis
// using text files jetsmear/data
// the good_jets input to this method are the final corrected jets used by the analysis
// after all selections on the jets
// the default selections should be
// |eta| < 4.7 and corrected pt > 10 GeV
//
void JetMETUncertainty::SmearMETForJERUncertainty ()
{
    if (!have_jetsmear_) {
        cout << "Attempting to smear MET for JER but do not have valid jet smearing service." << endl;
        met_jer_up_   = make_pair(-1., 0.);
        met_jer_down_ = make_pair(-1., 0.);
        return;
    }

    smeared_jets_ = smearJets(good_jets_, jetSmearer_);

    //
    // this is all just to make the up/down variation easier to deal with
    //
    smeared_jet_sf_.clear();
    smeared_jet_sf_.reserve(smeared_jets_.size());
    for (unsigned int idx = 0; idx < smeared_jets_.size(); idx++) {
        double diff = (smeared_jets_.at(idx).pt() - good_jets_.at(idx).pt()) / good_jets_.at(idx).pt();
        smeared_jet_sf_.push_back(diff - 1.);
    }
    assert(smeared_jet_sf_.size() == good_jets_.size());

    double dmet_x = good_met_.first * cos(good_met_.second);
    double dmet_y = good_met_.first * sin(good_met_.second);
    double umet_x = good_met_.first * cos(good_met_.second);
    double umet_y = good_met_.first * sin(good_met_.second);

    for (unsigned int idx = 0; idx < good_jets_.size(); idx++) {
        double px    = good_jets_.at(idx).px();
        double py    = good_jets_.at(idx).py();
        double scale = smeared_jet_sf_.at(idx);

        dmet_x -= DOWN * scale * px;
        dmet_y -= DOWN * scale * py;

        umet_x -= UP * scale * px;
        umet_y -= UP * scale * py;        
    }

    met_jer_up_   = make_pair(sqrt( pow(umet_x, 2) + pow(umet_y, 2) ), atan2(umet_y, umet_x));
    met_jer_down_ = make_pair(sqrt( pow(dmet_x, 2) + pow(dmet_y, 2) ), atan2(dmet_y, dmet_x));
}

//
// set input parameters
// 1) vector of good electrons as defined by analysis selections
// 2) vector of good muons as defined by analysis selections
// 3) vector of good jets - default should be |eta| < 4.7 and corrected pt > 10 GeV, but it can be defined by the user
//
void JetMETUncertainty::SetInputParameters (VofP4s &good_els, VofP4s &good_mus, VofP4s &good_jets)
{
    good_els_  = VofP4s(good_els);
    good_mus_  = VofP4s(good_mus);
    good_jets_ = VofP4s(good_jets);
}

//
// set jet correction uncertainty pararmeter file name
//
void JetMETUncertainty::SetJetCorrUncFileName(string jetcorr_unc_file_name)
{
    have_jetcorr_unc_ = false;
    jetcorr_unc_file_name_ = std::string(jetcorr_unc_file_name);
    if (!jetcorr_unc_file_name.empty()) {
        jetCorrectionUncertainty_ = new JetCorrectionUncertainty(jetcorr_unc_file_name_);
        if (jetCorrectionUncertainty_ != 0)
            have_jetcorr_unc_ = true;
    }
}

//
// set jet smearing parameter file names
//
void JetMETUncertainty::SetJetSmearFileNames(vector<string> &jetsmear_file_names)
{
    have_jetsmear_ = false;
    jetsmear_file_names_ = std::vector<std::string>(jetsmear_file_names);
    if (jetsmear_file_names_.size() == 3) {
        jetSmearer_ = makeJetSmearer(jetsmear_file_names_);
        if (jetSmearer_ != 0)
            have_jetsmear_ = true;
    }
}

//
// set electron energy scale uncertainty
// first argument is the barrel uncertainty
// second argument is the endcap uncertainty
//
void JetMETUncertainty::SetElectronUncertainty(double barrel_uncertainty, double endcap_uncertainty)
{
    ele_unc_b_ = barrel_uncertainty;
    ele_unc_e_ = endcap_uncertainty;
}

//
// set muon momentum scale uncertainty
//
void JetMETUncertainty::SetMuonUncertainty(double muon_uncertainty)
{
    mu_unc_ = muon_uncertainty;
}

//
// set unclustered energy scale uncertainty
//
void JetMETUncertainty::SetUnclusteredEnergyUncertainty(double unclustered_uncertainty)
{
    uncl_unc_ = unclustered_uncertainty;
}

//
// get electron energy scale uncertainties
// first entry is the barrel uncertainty
// second entry is the endcap uncertainty
//
dpair JetMETUncertainty::GetElectronUncertainty()
{
    return make_pair(ele_unc_b_, ele_unc_e_);
}

//
// get muon momentum scale uncertaintya
//
double JetMETUncertainty::GetMuonUncertainty()
{
    return mu_unc_;
}

//
// get unclustered energy scale uncertainty
//
double JetMETUncertainty::GetUnclusteredEnergyUncertainty()
{
    return uncl_unc_;
}

//
// get SF used to smear jets for JER uncertainty
//
vector<double> JetMETUncertainty::GetSmearedJetScaleFactors ()
{
    return smeared_jet_sf_;
}


//
// get JES uncertainty for jets
//
vector<double> JetMETUncertainty::GetJESUncertainties ()
{
    return jetcorr_uncertainties_;
}
