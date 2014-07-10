#ifndef ScanChain_h
#define ScanChain_h

// C++ includes
#include <string>
#include <vector>

// ROOT includes
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "Math/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

class babyMaker {

 public:

  babyMaker() {};
  ~babyMaker() {
    delete BabyFile_;
    delete BabyTree_;
  };

  void ScanChain(TChain*, std::string = "testSample", bool isMC = true);

  void MakeBabyNtuple(const char *);
  void InitBabyNtuple();
  void FillBabyNtuple();
  void CloseBabyNtuple();

 private:

  TFile *BabyFile_;
  TTree *BabyTree_;

  //baby ntuple variables
  float pfmet;
  float pfmet_type1cor;
  float genmet;
  float scale1fb;
  int isRealData;
  unsigned int event;
  unsigned int lumiBlock;
  unsigned int run;

  bool dielectronTrigger;
  bool dimuonTrigger;
  bool electronmuonTrigger;

  std::vector<LorentzVector> els_p4_b;
  std::vector<LorentzVector> mus_p4_b;
  std::vector<LorentzVector> pfjets_p4_b;

  std::vector<int> genps_id_b;
  std::vector<int> genps_id_mother_b;
  std::vector<LorentzVector> genps_p4_b;
  std::vector<LorentzVector> genjets_p4_b;

  std::vector<bool> looseEl;
  std::vector<bool> looseMu;
  std::vector<bool> tightEl;
  std::vector<bool> tightMu;
  std::vector<bool> passesLoosePFJetID_b;

  std::vector<float> corL1FastL2L3;

};

#endif

void babyMaker::MakeBabyNtuple(const char *BabyFilename){

  //
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  BabyFile_ = new TFile(Form("%s", BabyFilename), "RECREATE");
  BabyFile_->cd();
  BabyTree_ = new TTree("tree", "A Baby Ntuple");

  BabyTree_->Branch("pfmet", &pfmet );
  BabyTree_->Branch("pfmet_type1cor", &pfmet_type1cor );
  BabyTree_->Branch("gen_met", &genmet );
  BabyTree_->Branch("evt_scale1fb", &scale1fb );
  BabyTree_->Branch("evt_isRealData", &isRealData );
  BabyTree_->Branch("evt_event", &event );
  BabyTree_->Branch("evt_lumiBlock", &lumiBlock );
  BabyTree_->Branch("evt_run", &run );

  BabyTree_->Branch("dielectronTrigger", &dielectronTrigger );
  BabyTree_->Branch("dimuonTrigger", &dimuonTrigger );
  BabyTree_->Branch("electronmuonTrigger", &electronmuonTrigger );
  
  BabyTree_->Branch("els_p4", &els_p4_b );
  BabyTree_->Branch("mus_p4", &mus_p4_b );
  BabyTree_->Branch("pfjets_p4", &pfjets_p4_b );
  BabyTree_->Branch("genps_id", &genps_id_b );
  BabyTree_->Branch("genps_id_mother", &genps_id_mother_b );
  BabyTree_->Branch("genps_p4", &genps_p4_b );
  BabyTree_->Branch("genjets_p4", &genjets_p4_b );
  BabyTree_->Branch("looseEl", &looseEl );
  BabyTree_->Branch("looseMu", &looseMu );
  BabyTree_->Branch("tightEl", &tightEl );
  BabyTree_->Branch("tightMu", &tightMu );
  BabyTree_->Branch("passesLoosePFJetID", &passesLoosePFJetID_b );
  BabyTree_->Branch("pfjets_corL1FastL2L3", &corL1FastL2L3 );
   
  return;
}

void babyMaker::InitBabyNtuple () {

  pfmet = -999.0;
  pfmet_type1cor = -999.0;
  genmet = -999.0;
  scale1fb = -999.0;
  isRealData = -999;
  event = 999;
  lumiBlock = 999;
  run = 999;

  dielectronTrigger = false;
  dimuonTrigger = false;
  electronmuonTrigger = false;

  els_p4_b.clear();
  mus_p4_b.clear();
  pfjets_p4_b.clear();

  genps_id_b.clear();
  genps_id_mother_b.clear();
  genps_p4_b.clear();
  genjets_p4_b.clear();

  looseEl.clear();
  looseMu.clear();
  tightEl.clear();
  tightMu.clear();
  passesLoosePFJetID_b.clear();

  corL1FastL2L3.clear();
  
  return;
}

void babyMaker::FillBabyNtuple(){
  BabyTree_->Fill();
  return;
}

void babyMaker::CloseBabyNtuple(){
  BabyFile_->cd();
  BabyTree_->Write();
  BabyFile_->Close();
  return;
}
