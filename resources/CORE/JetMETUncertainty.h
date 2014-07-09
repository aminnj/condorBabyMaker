//-*- C++ -*-
//
// Package:    CMS2/NtupleMacros/CORE
// Class:      JetMETUncertainty
// 
/**\class JetMETUncertainty JetMETUcertainty.cc CMS2/NtupleMacros/CORE/JetMETUcertainty.cc

   Description: Class for determing jet and met uncertainties.


   Implementation: Following recommendations of JetMET POG:

   https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription
   https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#MET_Systematics_Tools
   https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution



   Usage: Use in looper (not necessarily cms2)
   
   1) Create list of input files, instance of class at beginning of looper.
   std::string jes_uncertainty_file;
   if (is_data)
      jes_uncertainty_file = "jetcorr/data/GR_R_42_V23_AK5PF_Uncertainty.txt"
   else
      jes_uncertainty_file = "jetcorr/data/DESIGN42_V17_AK5PF_Uncertainty.txt"

   std::vector<std::string> vector_of_jet_smearing_files;
   vector_of_jet_smearing_files.push_back("jetsmear/data/Spring10_PtResolution_AK5PF.txt");
   vector_of_jet_smearing_files.push_back("jetsmear/data/Spring10_PhiResolution_AK5PF.txt");
   vector_of_jet_smearing_files.push_back("jetsmear/data/jet_resolutions.txt");

   JetMETUncertainty *jetMETUncertainty = new Jetmetuncertainty(jes_uncertainty_file, vector_of_jet_smearing_files);   


   2) Modify any other uncertainties you want

   By default, uses the following uncertainies:
      a) electron energy scale uncertainty: 0.6% (1.5%) in the barrel (encap)
      b) muon momentum scale uncertainty: 1.0%
      c) unclustered energy scale uncertainty: 10%

   Modify whichever you want using:
      jetMETUncertainty->SetElectronUncertainty(double barrel_uncertainty, double endcap_uncertainty);
      jetMETUncertainty->SetMuonUncertainty(double muon_uncertainty);
      jetMETUncertainty->SetUnclusteredEnergyUncertainty(double unclustered_uncertainty);


   3) For each event, hyp, whatever:

   Set input parameters: jetMETUncertainty->SetInputParameters(VofP4s &good_els, VofP4s &good_mus, VofP4s &good_jets);

   Here, the inputs are defined as:
      a) good electrons as defined by analysis selection
      b) good muons as defined by analysis selection
      c) good jets: you should use |eta| < 4.7 && L1FastL2L3(Res) corrected pt > 10 GeV unless you know what you're doing

   Get the uncertainty (up, down): std::pair<double, double> met_unc = jetMETUncertainty->GetTotalUncertainty(std::pair<double, double> input_met);
      a) input met is pair of (met, met_phi)
      b) returns uncertainty on scaling met (up, down)


   4) How to use the uncertainties:
      a) Simultaneously vary MET and jets --UP-- and determine impact on yields passing full selection.
      b) Simultaneously vary MET and jets --DOWN-- and determine impact on yields passing full selection.
      c) Total JetMET uncertainty is max(UP, DOWN) variation in yields.

*/
//
// Original Author:  Frank Golf
//         Created:  Fri April  6 18:14 PDT 2012
//
//


#ifndef JETMETUNCERTAINTY_H
#define JETMETUNCERTAINTY_H

#include "Math/LorentzVector.h"
#include <vector>
#include <string>
#include <utility>
#include "jetcorr/JetCorrectionUncertainty.h"
#include "jetsmear/JetSmearer.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
typedef std::vector<LorentzVector> VofP4s;
typedef std::pair<double, double> dpair;

class JetMETUncertainty {
public:
    JetMETUncertainty  ();
    JetMETUncertainty  (std::string jetcorr_unc_file_name, std::vector<std::string> &jetsmear_file_names);
    ~JetMETUncertainty ()
    {
        delete jetCorrectionUncertainty_;
        delete jetSmearer_;
    }

    enum ScaleType { DOWN = -1, UP = 1 };
    enum Component { ELE = 0, MU, UNCL, JES, JER, N_COMPONENTS };

    void SetJetCorrUncFileName (std::string jetcorr_unc_file_name);
    void SetJetSmearFileNames (std::vector<std::string> &jetsmear_file_names);
    void SetElectronUncertainty(double barrel_uncertainty, double endcap_uncertainty);
    void SetMuonUncertainty(double muon_uncertainty);
    void SetUnclusteredEnergyUncertainty(double unclustered_uncertainty);
    void SetInputParameters (VofP4s &good_els, VofP4s &good_mus, VofP4s &good_jets);

    dpair  GetElectronUncertainty();
    double GetMuonUncertainty();
    double GetUnclusteredEnergyUncertainty();
    dpair  GetScaledMET (dpair good_met, Component component, ScaleType scale_type);
    dpair  GetTotalUncertainty(dpair good_met);

    std::vector<double> GetSmearedJetScaleFactors ();
    std::vector<double> GetJESUncertainties ();

private:
    
    double ele_unc_b_;     // electron energy scale uncertainty, barrel
    double ele_unc_e_;     // electron energy scale uncertainty, endcap
    double mu_unc_;        // muon energy scale uncertainty
    double uncl_unc_;      // unclustered energy uncertainty        

    std::string jetcorr_unc_file_name_;
    std::vector<std::string> jetsmear_file_names_;

    bool have_jetcorr_unc_;
    bool have_jetsmear_;

    JetCorrectionUncertainty *jetCorrectionUncertainty_;
    JetSmearer *jetSmearer_;
    
    dpair met_ele_up_;
    dpair met_ele_down_;
    dpair met_mu_up_;
    dpair met_mu_down_;
    dpair met_uncl_up_;
    dpair met_uncl_down_;
    dpair met_jes_up_;
    dpair met_jes_down_;
    dpair met_jer_up_;
    dpair met_jer_down_;
    dpair met_up_;
    dpair met_down_;

    VofP4s good_els_;
    VofP4s good_mus_;
    VofP4s good_jets_;
    dpair  good_met_;

    VofP4s smeared_jets_;
    std::vector<double> smeared_jet_sf_;
    std::vector<double> jetcorr_uncertainties_;

    void SmearMETForElectronUncertainty          ();
    void SmearMETForMuonUncertainty              ();
    void SmearMETForUnclusteredEnergyUncertainty ();
    void SmearMETForJESUncertainty               ();
    void SmearMETForJERUncertainty               ();

    dpair GetUnclusteredMET ();
};

#endif
