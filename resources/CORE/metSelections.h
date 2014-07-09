#ifndef METSELECTIONS_H
#define METSELECTIONS_H

#include "electronSelections.h"

//typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
typedef std::vector<LorentzVector> VofP4;

class FactorizedJetCorrector;

//
// met selections
//

struct metStruct{
  metStruct() : met(-999.), metphi(-999.), metx(-999.), mety(-999.), sumet(-999.)  {}
  float met;
  float metphi;
  float metx;
  float mety;
  float sumet;
};

enum whichMetType {
  usingTcMet = 1,
  usingCaloMet = 2,
  usingTcMet35X = 3,
  usingTcMet_looper = 4
};

//---------------------------------------------------
// Function that checks whether met (or tcmet) was 
// corrected for a given muon.  This uses the value maps
//---------------------------------------------------
bool wasMetCorrectedForThisMuon(int imu, whichMetType type);

//-----------------------------------------------------------
// Function that corrects the met (or tcmet) for a given
// muon in case this was not done in reco.  Uses value maps
//-------------------------------------------------------------
void fixMetForThisMuon(int imu, float& metX, float& metY, whichMetType type);

//-----------------------------------------------------------
// Function that corrects the met (or tcmet) for a given
// muon in case this was not done in reco.  Uses value maps
// same things as above but also corrects sumET
//-------------------------------------------------------------
void fixMetForThisMuon(int imu, float& metX, float& metY, float& sumET, whichMetType type);


//---------------------------------------------
// function to calculate projected MET.
// takes three parameters as input:
//
// 1. met
// 2. met phi
// 3. hypothesis index
//---------------------------------------------
float projectedMET( float met, float metPhi, int hyp_index );

//---------------------------------------------
// as above but simpler for single lepton events
//---------------------------------------------
float projectedMETW( float met, float metPhi, float leptonPhi);

//---------------------------------------------
// utility function find deltaPhi between met
// and nearest hypothesis lepton
//---------------------------------------------
float nearestHypLeptonPhi( float metPhi, int hyp_index);

//---------------------------------------------
// correct tcMET for any hypothesis muons
// that have not been corrected for
//---------------------------------------------
metStruct correctTCMETforHypMuons (int hyp_index, float met_x, float met_y, float sumet);

//---------------------------------------------
// function to calculate latest tcMET
//---------------------------------------------
metStruct correctedTCMET(bool printout = false, ostream& ostr = std::cout);

//---------------------------------------------
//function to perform custom Type1 correction
//---------------------------------------------
metStruct customType1Met( float metx , float mety , float sumet  , VofP4 jets , std::vector<float> cors );

//
metStruct trackerMET( int hyp_index, double deltaZCut = 0.2, 
                      const VofP4* jets = 0 );


//---------------------------------------------
// function to calculate CMS reduced MET
// (simplified implementation from Nate Odell)
//---------------------------------------------
LorentzVector cmsReducedMET(LorentzVector sumJet, LorentzVector lep1, LorentzVector lep2, LorentzVector metP4, int version);
std::pair<float, float> cmsReducedMET_v2(LorentzVector lep1, LorentzVector lep2, const std::vector<LorentzVector> &jets);

//-----------------------------------------------------
// function to scale the hadronic component of the MET
//-----------------------------------------------------
std::pair<float, float> scaleMET(std::pair<float, float> p_met, LorentzVector p4_dilep, double rescale = 1.0);


class MetCorrector
{
public:
    MetCorrector (std::vector<std::string> &list_of_files);
    ~MetCorrector ();
    std::pair<float, float> getCorrectedMET(std::pair<float, float> &uncorr_met);
    std::pair<float, float> getCorrectedMET();

private:
    std::pair<float, float> correctMETforJES(std::pair<float, float>);

    FactorizedJetCorrector *offset_corrector;
    FactorizedJetCorrector *full_corrector;
};

#endif

