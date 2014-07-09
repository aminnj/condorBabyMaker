#ifndef Resolutions_cxx
#define Resolutions_cxx
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TRandom3.h"

#include <vector>
#include <iostream>

#ifndef __CINT__
#include "JR_Standalone/JetResolution.h"
#include "JR_Standalone/JetResolution.cc"
#include "JR_Standalone/JetCorrectorParameters_tm.h"
#include "JR_Standalone/JetCorrectorParameters_tm.cc"
#include "JR_Standalone/Utilities.cc"
#endif

using namespace std;

void smear_JetMET(vector <TLorentzVector> & orig_jets, const TVector2 & orig_met, vector <TLorentzVector> & smear_jets, TVector2 & smear_met, TRandom3* rand3, vector <TF1*> vPtRes, vector <TF1*> vEtaRes, vector <TF1*> vPhiRes, TLorentzVector lep_sum);
void Dont_smear_JetMET(vector <TLorentzVector> & orig_jets, const TVector2 & orig_met, vector <TLorentzVector> & smear_jets, TVector2 & smear_met);

void smear_JetMET(vector <TLorentzVector> & orig_jets, const TVector2 & orig_met, vector <TLorentzVector> & smear_jets, TVector2 & smear_met, TRandom3* rand3, vector <TF1*> vPtRes, vector <TF1*> vEtaRes, vector <TF1*> vPhiRes, TLorentzVector lep_sum){

  smear_jets.clear();
  
  double sum_jpx = 0;
  double sum_jpy = 0;
  
  double sum_jpx_sm = 0;
  double sum_jpy_sm = 0;

  double Pt_sm, Eta_sm, Phi_sm;

  TLorentzVector v_temp;

//   double a, b, c;
//   double Et;

  for (unsigned int sui = 0; sui < orig_jets.size(); sui++){

//     a = 1.84; b = 1.14; c = 0.027;
    
//     if (orig_jets.at(sui).Eta() < 1.4){
//       a = 1.78; b = 1.3; c = 0.053;
//     }

//     Et = orig_jets.at(sui).Pt();
    
//     Pt_sm = Et + sqrt(a*a + b*b*Et + c*c*Et*Et) * rand3->Gaus();

//     Eta_sm = orig_jets.at(sui).Eta();
//     Phi_sm = orig_jets.at(sui).Phi();

    Long64_t iseed = (Long64_t) 10E10;
    
    gRandom->SetSeed(rand3->Integer(iseed));
//     Pt_sm  = orig_jets.at(sui).Pt()  * (1 + (vPtRes.at(sui)->GetRandom() - 1) * 1.1);
    Pt_sm  = orig_jets.at(sui).Pt()  * vPtRes.at(sui)->GetRandom();

    gRandom->SetSeed(rand3->Integer(iseed));
    Eta_sm = orig_jets.at(sui).Eta() + vEtaRes.at(sui)->GetRandom();
    
    gRandom->SetSeed(rand3->Integer(iseed));
    Phi_sm = orig_jets.at(sui).Phi() + vPhiRes.at(sui)->GetRandom();
    
    v_temp.SetPtEtaPhiM(Pt_sm, Eta_sm, Phi_sm, orig_jets.at(sui).M());
   
    sum_jpx += orig_jets.at(sui).Px();
    sum_jpy += orig_jets.at(sui).Py();

    sum_jpx_sm += v_temp.Px();
    sum_jpy_sm += v_temp.Py();

    smear_jets.push_back(v_temp);
  }

  double unclust_metx = orig_met.Px() + sum_jpx + lep_sum.Px();
  double unclust_mety = orig_met.Py() + sum_jpy + lep_sum.Py();

  //10% resolution
  double unclust_metx_sm = unclust_metx * (1 + 0.1*rand3->Gaus());
  double unclust_mety_sm = unclust_mety * (1 + 0.1*rand3->Gaus());

  smear_met.Set(orig_met.Px() + sum_jpx - unclust_metx - sum_jpx_sm + unclust_metx_sm, orig_met.Py() + sum_jpy - unclust_mety - sum_jpy_sm + unclust_mety_sm);
}

void Dont_smear_JetMET(vector <TLorentzVector> & orig_jets, const TVector2 & orig_met, vector <TLorentzVector> & smear_jets, TVector2 & smear_met){ 
  smear_met  = orig_met;
  smear_jets = orig_jets;
}
#endif
