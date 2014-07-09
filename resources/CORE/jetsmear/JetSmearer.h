#ifndef JETSMEARER_H
#define JETSMEARER_H

#include <string>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "JetResolution.h"
#include <vector>
#include <map>
#include <utility>
#include "SigInputObj.h"
#include <fstream>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

enum resolutionType { caloEE, caloEB, caloHE, caloHO, caloHF, caloHB, jet, electron, tau, muon,PFtype1,PFtype2, PFtype3, PFtype4, PFtype5, PFtype6, PFtype7 };
enum resolutionFunc { ET, PHI,TRACKP,CONSTPHI };

class JetSmearer 
{
public:
    JetSmearer ();
    JetSmearer (const std::string ptFileName, const std::string phiFileName, const std::string resFileName);
    ~JetSmearer ();
    void setResFileNames (const std::string ptFileName, const std::string phiFileName, const std::string resFileName);
    LorentzVector smearJet (LorentzVector p4);
    double getJetResolution(LorentzVector p4);
    void setDeltaR (double dr);
    double getJetPtThreshold ();
    void setDelimiter (std::string);
    std::string getDelimiter ();

private:
    double deltaR_;
    TRandom3* rand_;
    std::pair<double, double> getKjet(LorentzVector p4);
    double getRjet(LorentzVector p4);
    int matchRecoJetToGenJet(LorentzVector p4);
    double getRandom (double sigma, double mean = 0.);
    JetResolution* ptResol_;
    JetResolution* phiResol_;
    void initializeJetResolutions (const std::string ptFileName, const std::string phiFileName);
    void addResolutions (const std::string resFileName);
    double ptResolThreshold_;
    //temporary fix for low pT jet resolutions
    //First index, eta bins, from 0 to 5;
    //Second index, pt bins, from 3 to 23 GeV;
    std::vector<double> jdpt[10];
    std::vector<double> jdphi[10];
    typedef std::pair<resolutionType, resolutionFunc> functionCombo;
    typedef std::vector<double> functionPars;
    std::map<functionCombo,functionPars> functionmap_;
    void addfunction(const resolutionType type, const resolutionFunc func, std::vector<double> parameters);
    std::string getLine(std::ifstream *filep, std::string identifier);
    SigInputObj evalPFJet(const LorentzVector p4)  const;
    std::string res_delim_;
};

#endif
