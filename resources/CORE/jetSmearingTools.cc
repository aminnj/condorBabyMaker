#include "jetSmearingTools.h"
#include <iostream>
#include "jetsmear/JetSmearer.h"
#include "jetsmear/SigInputObj.h"
#include "jetsmear/JetResolution.h"
#include "jetSelections.h"

// function to smear jet energy to account for differences in data-MC JER
JetSmearer* makeJetSmearer(const char* ptFileName, const char* phiFileName, const char* resFileName)
{
    std::vector<std::string> vector_of_file_names;
    vector_of_file_names.reserve(3);
    vector_of_file_names.push_back(ptFileName);
    vector_of_file_names.push_back(phiFileName);
    vector_of_file_names.push_back(resFileName);
    return makeJetSmearer(vector_of_file_names);
}

JetSmearer* makeJetSmearer(std::vector<std::string> &vector_of_file_names)
{
    vector<std::string> vParam;
    for (std::vector<std::string>::const_iterator i = vector_of_file_names.begin(), i_end = vector_of_file_names.end(); i != i_end; ++i) {
        // do some rigmarole to evaluate env variables in the strings
        // std::cout << "file name: " << *i << std::endl;
        const std::string cmd = "echo ";
        FILE *f = popen((cmd + *i).c_str(), "r");
        if (!f) {
            perror((std::string("Opening pipe to execute ") + cmd + *i).c_str());
            return 0;
        }
        char corr_name[1024];
        int s = fscanf(f, " %1024s\n", corr_name);
        // std::cout << "s = " << s << std::endl;
        if (s != 1) {
            perror("reading file list");
        }
        assert(s == 1);
        // JetCorrectorParameters JetCorPar(corr_name);
        // printf("%s\n", corr_name);
        vParam.push_back(std::string(corr_name));
        fclose(f);
    }

    assert(vParam.size() == 3);
    // for (unsigned int idx = 0; idx < vParam.size(); idx++)
    //     std::cout << "file" << idx << ": " << vParam.at(idx) << std::endl;
    return (new JetSmearer(vParam[0], vParam[1], vParam[2]));
}

LorentzVector smearJet(LorentzVector p4, JetSmearer* jetSmearer)
{
    return jetSmearer->smearJet(p4);
}

std::vector<LorentzVector> smearJets(std::vector<LorentzVector>& vp4s, JetSmearer* jetSmearer)
{
    std::vector<LorentzVector> smeared_jets;
    for (unsigned int idx = 0; idx < vp4s.size(); idx++) {
        LorentzVector smeared_jet = jetSmearer->smearJet(vp4s.at(idx));
        smeared_jets.push_back(smeared_jet);
    }

    return smeared_jets;
}

//-----------------------------------------------------
// function to smear the MET for differences in jet
// resolution between data and MC
//-----------------------------------------------------
std::pair<double, double> smearMETforJER(std::pair<double, double> in_met, JetSmearer *jetSmearer, std::vector<LorentzVector>& vp4s) {

    double met     = in_met.first;
    double met_phi = in_met.second;
    double met_x   = met * cos(met_phi);
    double met_y   = met * sin(met_phi);

    std::vector<LorentzVector> vsmeared_jets = smearJets(vp4s, jetSmearer);
    assert(vp4s.size() == vsmeared_jets.size());
    for (unsigned int idx = 0; idx < vp4s.size(); idx++) {
        met_x += (vp4s.at(idx).px() - vsmeared_jets.at(idx).px());
        met_y += (vp4s.at(idx).py() - vsmeared_jets.at(idx).py());
    }
    
    delete jetSmearer;

    return (std::make_pair(sqrt(met_x * met_x + met_y * met_y), atan2(met_y, met_x)));
}


//-----------------------------------------------------
// get (relative) resolution of a jet
// NOTE: this gets the resolution of a MC jet only
//-----------------------------------------------------
double getJetResolution(LorentzVector p4, JetSmearer* jetSmearer)
{
    return jetSmearer->getJetResolution(p4);
}

std::vector<double> getJetResolutions(std::vector<LorentzVector>& vp4s, JetSmearer* jetSmearer)
{
    std::vector<double> jet_rel_res;
    for (unsigned int idx = 0; idx < vp4s.size(); idx++)
    {
        jet_rel_res.push_back(getJetResolution(vp4s.at(idx), jetSmearer));
    }

    return jet_rel_res;
}
