
#include "JetSmearer.h"
#include "JetResolution.h"
#include "../CMS2.h"
#include "TString.h"
#include "TObjArray.h"
#include "TF1.h"

void JetSmearer::setDeltaR (double dr) {
    
    deltaR_ = dr;
    return;
}

LorentzVector JetSmearer::smearJet (LorentzVector p4)
{
    double scale = 1.;

    int midx = matchRecoJetToGenJet(p4);
    double kjet = getKjet(p4).first;
    double jet_energy = p4.energy();
    double rjet = getRjet(p4);
    if (midx < 0) {
        double sigma = sqrt(kjet*kjet - 1.) * rjet;
        double mean = 0.;
        double rand = getRandom(mean, sigma);
        scale = 1. + (rand / jet_energy);
    }
    else {
        double genjet_energy = cms2.genjets_p4().at(midx).energy();
        scale = 1. + (kjet - 1.) * ((jet_energy - genjet_energy) / jet_energy);
    }

    return p4 * scale;
}

JetSmearer::JetSmearer() : deltaR_(0.5), res_delim_(",")
{
    rand_ = new TRandom3();
    ptResol_ = 0;
    phiResol_ = 0;
}

JetSmearer::JetSmearer(const std::string ptFileName, const std::string phiFileName, const std::string resFileName) : deltaR_(0.5), res_delim_(",")
{
    rand_ = new TRandom3();
    ptResol_ = 0;
    phiResol_ = 0;
    // std::cout << "ptFileName: " << ptFileName << std::endl;
    // std::cout << "phiFileName: " << phiFileName << std::endl;
    // std::cout << "resFileName: " << resFileName << std::endl;
    initializeJetResolutions(ptFileName, phiFileName);
    addResolutions(resFileName);
}

std::pair<double, double> JetSmearer::getKjet (LorentzVector p4) {

    double eta = fabs(p4.eta());

    if (eta < 1.1)
        return make_pair(1.062, 0.028);
    else if (eta < 1.7)
        return make_pair(1.084, 0.034);
    else if (eta < 2.3)
        return make_pair(1.029, 0.048);
    else if (eta < 4.7)
        return make_pair(1.153, 0.076);
    else
        return make_pair(1., 0.);
}

double JetSmearer::getRjet (LorentzVector p4) {

    SigInputObj jet_res_obj = evalPFJet(p4);
    if (jet_res_obj.get_energy() > 0)
        return (p4.energy() * (jet_res_obj.get_sigma_e() / jet_res_obj.get_energy()));
    else
        return 0.;
}

int JetSmearer::matchRecoJetToGenJet(LorentzVector p4) {
    
    double mindr = deltaR_;
    int midx = -1;
    for (unsigned int idx = 0; idx < cms2.genjets_p4().size(); idx++) {
        
        double tmp_dr = ROOT::Math::VectorUtil::DeltaR(p4, cms2.genjets_p4().at(idx));
        if (tmp_dr < mindr) {
            mindr = tmp_dr;
            midx = idx;
        }        
    }

    return midx;
}

double JetSmearer::getRandom (double sigma, double mean) {
    
    return rand_->Gaus(mean, sigma);
}

void JetSmearer::initializeJetResolutions (const std::string ptFileName, const std::string phiFileName) {

    // std::cout << "ptFileName: " << ptFileName << std::endl;
    // std::cout << "phiFileName: " << phiFileName << std::endl;
  
    if (ptResol_ == 0) {
        ptResol_  = new JetResolution(ptFileName, false);
        phiResol_ = new JetResolution(phiFileName, false);     
    }

    return;
}

void JetSmearer::setResFileNames (const std::string ptFileName, const std::string phiFileName, const std::string resFileName) {

    initializeJetResolutions(ptFileName, phiFileName);
    addResolutions(resFileName);

    return;
}

void JetSmearer::addResolutions (const std::string resFileName) {

    TString line;
    std::ifstream* res_file = new std::ifstream();
    res_file->open(resFileName.c_str());
    TObjArray *tokens = 0;
    TObjString *objtoken = 0;
    int num_tokens = -1;

    // perform some quick sanity checks that we're using the file we think we are
    line = TString(getLine(res_file, "resolutionsEra"));
    tokens = line.Tokenize(TString(res_delim_));
    num_tokens = tokens->GetEntries();
    if (num_tokens == 2) {
        objtoken = (TObjString*)tokens->At(1);
        TString token = objtoken->GetString();
        if (!token.Contains("Spring10")) {
            std::cout << "didn't fine resolutionsEra.  Exiting." << std::endl;
            res_file->close();
            delete res_file;
            delete tokens;
            return;
        }
    }
    else {
        std::cout << "found " << num_tokens << " tokens when 2 expected.  exiting." << std::endl;
        res_file->close();
        delete res_file;
        delete tokens;

        return;
    }

    line = TString(getLine(res_file, "resolutionsAlgo"));
    tokens = line.Tokenize(TString(res_delim_));
    num_tokens = tokens->GetEntries();
    if (num_tokens == 2) {
        objtoken = (TObjString*)tokens->At(1);
        TString token = objtoken->GetString();
        if (!token.Contains("AK5PF")) {
            std::cout << "didn't fine resolutionsAlgo.  Exiting." << std::endl;
            res_file->close();
            delete res_file;
            delete tokens;
            return;   
        }
    }
    else {
        std::cout << "found " << num_tokens << " tokens when 2 expected.  exiting." << std::endl;
        res_file->close();
        delete res_file;
        delete tokens;
        return;   
    }

    // get the pt threshold
    line = TString(getLine(res_file, "ptresolthreshold"));
    tokens = line.Tokenize(TString(res_delim_));
    num_tokens = tokens->GetEntries();
    if (num_tokens == 2) {
        objtoken = (TObjString*)tokens->At(1);
        TString token = objtoken->GetString();
        ptResolThreshold_ = token.Atof();
    }

    // get pf numbers
    std::vector<double> et_params;
    std::vector<double> phi_params;

    for (unsigned int nparams = 1; nparams <= 7; nparams++) {
        et_params.clear();
        phi_params.clear();

        line = TString(getLine(res_file, Form("PF_EtResType%d", nparams)));
        tokens = line.Tokenize(TString(res_delim_));
        num_tokens = tokens->GetEntries();
        for (int idx = 1; idx < num_tokens; idx++) {
            objtoken = (TObjString*)tokens->At(idx);
            TString token = objtoken->GetString();
            et_params.push_back(token.Atof());
        }

        line = TString(getLine(res_file, Form("PF_PhiResType%d", nparams)));
        tokens = line.Tokenize(TString(res_delim_));
        num_tokens = tokens->GetEntries();
        for (int idx = 1; idx < num_tokens; idx++) {
            objtoken = (TObjString*)tokens->At(idx);
            TString token = objtoken->GetString();
            phi_params.push_back(token.Atof());
        }

        addfunction(resolutionType(PFtype1+nparams-1), ET, et_params);
        addfunction(resolutionType(PFtype1+nparams-1), PHI, phi_params);
    } 

    // get temporary low pT pfjet resolutions
    for (int ieta = 0; ieta < 10; ieta++) {
        et_params.clear();
        phi_params.clear();

        line = TString(getLine(res_file, Form("jdpt%d", ieta)));
        tokens = line.Tokenize(TString(res_delim_));
        num_tokens = tokens->GetEntries();
        for (int idx = 1; idx < num_tokens; idx++) {
            objtoken = (TObjString*)tokens->At(idx);
            TString token = objtoken->GetString();
            et_params.push_back(token.Atof());
        }

        line = TString(getLine(res_file, Form("jdphi%d", ieta)));
        tokens = line.Tokenize(TString(res_delim_));
        num_tokens = tokens->GetEntries();
        for (int idx = 1; idx < num_tokens; idx++) {
            objtoken = (TObjString*)tokens->At(idx);
            TString token = objtoken->GetString();
            phi_params.push_back(token.Atof());
        }
        
        jdpt[ieta]  = et_params;
        jdphi[ieta] = phi_params;
    }

    // clean-up after ourselves    
    res_file->close();
    delete res_file;
    delete tokens;

    return;
}

void JetSmearer::addfunction(const resolutionType type, const resolutionFunc func, std::vector<double> parameters) {

    functionCombo mypair(type,func);
    functionmap_[mypair]=parameters;
    return;
}

std::string JetSmearer::getLine(std::ifstream *filep, std::string identifier) {

    // need to re-set the file for reading
    filep->clear();
    filep->seekg(0, std::ios::beg);

    std::string matched_line = std::string();
    bool found_match = false;
    TObjArray *tokens = new TObjArray();
    TObjString *objtoken = new TObjString();
    if (filep->is_open()) {
        std::string line;
        while (filep->good()) {
            getline(*filep, line);
            TString tline(line);
            if (tline.BeginsWith("#"))
                continue;
            tokens = tline.Tokenize(TString(res_delim_));
            if (tokens->GetEntries() > 0) {
                objtoken = (TObjString*)tokens->At(0);
                TString token = objtoken->GetString();
                if (token.Contains(identifier.c_str())) {
                    found_match = true;
                    matched_line = line;
                }
            }
        }
    }    

    if (!found_match)
        std::cout << "didn't find a line with identifier " << identifier << std::endl;

    // clean-up after ourselves
    delete tokens;

    return matched_line;
}

SigInputObj JetSmearer::evalPFJet(const LorentzVector p4) const {
    
    double jpt  = p4.pt();
    double jphi = p4.phi();
    double jeta = p4.eta();
    double jdeltapt = 999.;
    double jdeltapphi = 999.;

    if(jpt < ptResolThreshold_ && jpt < 20.){ //use temporary fix for low pT jets
        double feta = TMath::Abs(jeta);
        int ieta = feta<5.? int(feta/0.5) : 9; //bin size = 0.5 
        int ipt  = jpt>3. ? int(jpt-3./2) : 0; //bin size =2, starting from ptmin=3GeV
        jdeltapt   = jdpt[ieta][ipt];
        jdeltapphi = jpt*jdphi[ieta][ipt];
    }
    else{
        //use the resolution functions at |eta|=5 to avoid crash for jets with large eta.
        if(jeta>5) jeta=5;
        if(jeta<-5) jeta=-5;
        TF1* fPtEta  = ptResol_->parameterEta("sigma",jeta);
        TF1* fPhiEta = phiResol_->parameterEta("sigma",jeta);
        jdeltapt   = jpt>ptResolThreshold_ ? jpt*fPtEta->Eval(jpt)  : jpt*fPtEta->Eval(ptResolThreshold_);
        jdeltapphi = jpt>ptResolThreshold_ ? jpt*fPhiEta->Eval(jpt) : jpt*fPhiEta->Eval(ptResolThreshold_);
        delete fPtEta;
        delete fPhiEta;
    }

    std::string inputtype = "jet";
    SigInputObj obj_jet(inputtype,jpt,jphi,jdeltapt,jdeltapphi);
    //std::cout << "RESOLUTIONS JET: " << jpt << "   " << jphi<< "   " <<jdeltapt << "   " << jdeltapphi << std::endl;
    return obj_jet;
}

double JetSmearer::getJetPtThreshold ()
{
    return ptResolThreshold_;
}

void JetSmearer::setDelimiter (std::string delimiter)
{
    res_delim_ = delimiter;
}

std::string JetSmearer::getDelimiter ()
{
    return res_delim_;
}

JetSmearer::~JetSmearer ()
{
    delete rand_;
    delete ptResol_;
    delete phiResol_;
}

//-----------------------------------------------------
// get (relative) resolution of a jet
// NOTE: this gets the resolution of a MC jet only
//-----------------------------------------------------
double JetSmearer::getJetResolution(LorentzVector p4)
{
    double rjet = getRjet(p4);
    return rjet/p4.energy();
}
