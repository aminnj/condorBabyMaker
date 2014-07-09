// -*- C++ -*-

// $Id: jetSelections.h,v 1.24 2012/08/25 00:32:31 fgolf Exp $

#ifndef JETSELECTIONS_H
#define JETSELECTIONS_H

#include <vector>
#include "CMS2.h"

//set to 1 if you want to use gen jets
#define haveGEN 0

enum JetType {
    JETS_TYPE_PF_FAST_CORR_RESIDUAL,
    JETS_TYPE_PF_FAST_CORR,
    JETS_TYPE_PF_CORR,
    JETS_TYPE_PF_UNCORR,
    JETS_TYPE_CALO_CORR,
    JETS_TYPE_CALO_UNCORR,
#if haveGEN
    JETS_TYPE_GEN
#endif
};

enum CleaningType {
    JETS_CLEAN_NONE,			// dirty, dirty jets
    JETS_CLEAN_HYP_E_MU,		// e or mu from hypothesis
    JETS_CLEAN_HYP_E,			// e only from hypothesis
    JETS_CLEAN_SINGLE_E,		// e in single-lepton final state (QCD)
    JETS_CLEAN_SINGLE_MU,      // mu in single-lepton final state (QCD)
    JETS_CLEAN_SS_E_MU         // mu in single-lepton final state (QCD)
};

enum BtagType {
    JETS_BTAG_NONE,
    JETS_BTAG_TCHEL,
    JETS_BTAG_TCHEM,
    JETS_BTAG_TCHPM,
    JETS_BTAG_TCHPT,
    JETS_BTAG_SSVHEM,
    JETS_BTAG_SSVHPT,
    JETS_BTAG_CSVL,
    JETS_BTAG_CSVM,
    JETS_BTAG_CSVT
};

enum JetScaleType {
    JETS_SCALE_DOWN = -1,
    JETS_SCALE_UP = 1
};

static const float BtagWP[] = {-999999., 1.7, 3.3, 1.93, 3.41, 1.74, 2.00, 0.244, 0.679, 0.898};

#define JET_DEFAULT_TYPE 	JETS_TYPE_PF_UNCORR
#define JET_DEFAULT_CLEANING 	JETS_CLEAN_HYP_E_MU
#define JET_DEFAULT_DR		0.4
#define JET_DEFAULT_PT		30
#define JET_DEFAULT_ETA		2.5
#define JETS_DEFAULT_BTAG  JETS_BTAG_NONE

// vector of p4's of the jets passing selections
std::vector<LorentzVector> getJets (unsigned int i_hyp,  // hyp or single-e to use for cleaning
                                    bool sort = false,
                                    enum JetType = JET_DEFAULT_TYPE,
                                    enum CleaningType = JET_DEFAULT_CLEANING,
                                    double deltaR = JET_DEFAULT_DR,
                                    double min_pt = JET_DEFAULT_PT,
                                    double max_eta = JET_DEFAULT_ETA,
                                    double rescale = 1.0,
				    int systFlag = 0);

// vector of bools aligned with the jet collection telling you which
// jets passed the selections
std::vector<bool> getJetFlags (unsigned int i_hyp,  // hyp or single-e to use for cleaning
                               enum JetType = JET_DEFAULT_TYPE,
                               enum CleaningType = JET_DEFAULT_CLEANING,
                               double deltaR = JET_DEFAULT_DR,
                               double min_pt = JET_DEFAULT_PT,
                               double max_eta = JET_DEFAULT_ETA,
                               double rescale = 1.0,
			       int systFlag = 0);

// number of jets passing selections
int nJets (unsigned int i_hyp,  // hyp or single-e to use for cleaning
           enum JetType = JET_DEFAULT_TYPE,
           enum CleaningType = JET_DEFAULT_CLEANING,
           double deltaR = JET_DEFAULT_DR,
           double min_pt = JET_DEFAULT_PT,
           double max_eta = JET_DEFAULT_ETA,
           double rescale = 1.0,
	   int systFlag = 0);

// scalar sum pt of jets passing selections
double sumPt (unsigned int i_hyp,  // hyp or single-e to use for cleaning
              enum JetType = JET_DEFAULT_TYPE,
              enum CleaningType = JET_DEFAULT_CLEANING,
              double deltaR = JET_DEFAULT_DR,
              double min_pt = JET_DEFAULT_PT,
              double max_eta = JET_DEFAULT_ETA,
              double rescale = 1.0,
	      int systFlag = 0);

// code to retrieve jet corrections from jet-correction text files
class FactorizedJetCorrector;
FactorizedJetCorrector *makeJetCorrector (const char *l2corr 		 = "$CMSSW_BASE/src/CondFormats/JetMETObjects/data/Spring10_L2Relative_AK5Calo.txt", 
                                          const char *l3corr 		 = "$CMSSW_BASE/src/CondFormats/JetMETObjects/data/Spring10_L3Absolute_AK5Calo.txt", 
                                          const char *l2l3_residual_corr = "$CMSSW_BASE/src/CondFormats/JetMETObjects/data/Spring10DataV1_L2L3Residual_AK5Calo.txt");
FactorizedJetCorrector *makeJetCorrector (const std::vector<std::string> &vector_of_file_names);
// either specify a jet corrector in the call to jetCorrection 
double jetCorrection (const LorentzVector &jet, 
                      FactorizedJetCorrector *jetCorrector);
// or set it once and have it be used it in all subsequent calls to jetCorrection
void setJetCorrector (FactorizedJetCorrector *);
double jetCorrection (const LorentzVector &jet);
double jetCorrection (int ijet);
bool jetPassesLooseJetID(int ijet);
bool passesCaloJetID (const LorentzVector &jetp4);
bool passesPFJetID(unsigned int pfJetIdx);
bool passesLoosePFJetID(unsigned int pfJetIdx);

float randomConeEventDensity();

float jetDz(int ijet, int ivtx);

// vector of p4's of the jets passing selections
std::vector<LorentzVector> getBtaggedJets (unsigned int i_hyp,  // hyp or single-e to use for cleaning
                                           bool sort = false,
                                           enum JetType = JET_DEFAULT_TYPE,
                                           enum CleaningType = JET_DEFAULT_CLEANING,
                                           enum BtagType = JETS_DEFAULT_BTAG,
                                           double deltaR = JET_DEFAULT_DR,
                                           double min_pt = JET_DEFAULT_PT,
                                           double max_eta = JET_DEFAULT_ETA,
                                           double rescale = 1.0,
					   int systFlag = 0);

// vector of bools aligned with the jet collection telling you which
// jets passed the selections
std::vector<bool> getBtaggedJetFlags (unsigned int i_hyp,  // hyp or single-e to use for cleaning
                                      enum JetType = JET_DEFAULT_TYPE,
                                      enum CleaningType = JET_DEFAULT_CLEANING,
                                      enum BtagType = JETS_DEFAULT_BTAG,
                                      double deltaR = JET_DEFAULT_DR,
                                      double min_pt = JET_DEFAULT_PT,
                                      double max_eta = JET_DEFAULT_ETA,
                                      double rescale = 1.0,
				      int systFlag = 0);

// number of jets passing selections
int nBtaggedJets (unsigned int i_hyp,  // hyp or single-e to use for cleaning
                  enum JetType = JET_DEFAULT_TYPE,
                  enum CleaningType = JET_DEFAULT_CLEANING,
                  enum BtagType = JETS_DEFAULT_BTAG,
                  double deltaR = JET_DEFAULT_DR,
                  double min_pt = JET_DEFAULT_PT,
                  double max_eta = JET_DEFAULT_ETA,
                  double rescale = 1.0,
		  int systFlag = 0);

// the jet met systematic...see the cc file for details
float getJetMetSyst(int flag, float pt, float eta);


// this function calculates the fraction of the pt of charged particles in a jet associated to the vertex ivtx
float pfjet_beta(int ijet, int power = 1, float dzcut = 0.05 , int ivtx = 0, bool verbose = false );

#endif // SEL_JETS_H
