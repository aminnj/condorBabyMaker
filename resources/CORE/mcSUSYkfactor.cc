#include <math.h>
#include <stdlib.h>
#include <set>
#include "TDatabasePDG.h"
#include "Math/VectorUtil.h"
#include "CMS2.h"
#include "mcSUSYkfactor.h"
#include "Math/LorentzVector.h"
#include "mcSelections.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TAxis.h"
#include <iostream>

#include "mcSelections.h"
//#include "CMS2.cc"


typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

int sfinalState(int ipart1, int ipart2) {
// #                   ng     neutralino/chargino + gluino                     #
// #                   ns     neutralino/chargino + squark                     #
// #                   nn     neutralino/chargino pair combinations            #
// #                   ll     slepton pair combinations                        #
// #                   sb     squark-antisquark                                #
// #                   ss     squark-squark                                    #
// #                   tb     stop-antistop                                    #
// #                   bb     sbottom-antisbottom                              #
// #                   gg     gluino pair                                      #
// #                   sg     squark + gluino                                  #

 int index = -1;

// Ugly can be done in a better way using list/objects ... oh well
 // ng case
 if ((abs(ipart1) > 1000021 && abs(ipart1) < 1000039 && abs(ipart2) == 1000021) ||
      (abs(ipart2) > 1000021 && abs(ipart2) < 1000039 && abs(ipart1) == 1000021)) index = 0;
// ns case 
 else if ((abs(ipart1) > 1000021 && abs(ipart1) < 1000039 && (abs(ipart2) < 1000007 || (abs(ipart2) > 2000000 && abs(ipart2) < 2000007))) ||
    (abs(ipart2) > 1000021 && abs(ipart2) < 1000039 && (abs(ipart1) < 1000007 || (abs(ipart1) > 2000000 && abs(ipart1) < 2000007)))) index = 1;
// nn case 
 else if (abs(ipart1) > 1000021 && abs(ipart1) < 1000039 &&  abs(ipart2) > 1000021 && abs(ipart2) < 1000039) index = 2;
// ll case
 else if (((abs(ipart1) > 1000010 && abs(ipart1) < 1000017) || (abs(ipart1) > 2000010 && abs(ipart1) < 2000016)) &&
         ((abs(ipart2) > 1000010 && abs(ipart2) < 1000017) || (abs(ipart2) > 2000010 && abs(ipart2) < 2000016))) index = 3;
// tb case
 else if ((abs(ipart1) == 1000006 || abs(ipart1) == 2000006) && (abs(ipart2) == 1000006 || abs(ipart2) == 2000006)) index = 6;
// bb case
 else if ((abs(ipart1) == 1000005 || abs(ipart1) == 2000005) && (abs(ipart2) == 1000005 || abs(ipart2) == 2000005)) index = 7;
// sb case
 else if ((ipart1 * ipart2 < 0) &&
         (abs(ipart1) < 1000007 || (abs(ipart1) > 2000000 && abs(ipart1) < 2000007)) &&
         (abs(ipart2) < 1000007 || (abs(ipart2) > 2000000 && abs(ipart2) < 2000007))) index = 4;
// ss case
 else if ((ipart1 * ipart2 > 0) &&
         (abs(ipart1) < 1000007 || (abs(ipart1) > 2000000 && abs(ipart1) < 2000007)) &&
         (abs(ipart2) < 1000007 || (abs(ipart2) > 2000000 && abs(ipart2) < 2000007))) index = 5;
// gg case 
  else if (ipart1 == 1000021 && ipart2 == 1000021) index = 8;
// sg case 
 else if ((abs(ipart1) == 1000021 && (abs(ipart2) < 1000007 || (abs(ipart2) > 2000000 && abs(ipart2) < 2000007))) ||
    (abs(ipart2) == 1000021 && (abs(ipart1) < 1000007 || (abs(ipart1) > 2000000 && abs(ipart1) < 2000007)))) index = 9;
 else
   cout << "Warning mcSUSYkfactor: index not defined " << ipart1 << "  " << ipart2 << endl;

  return index;
}


Float_t GetValueTH2FS(Float_t x, Float_t y, TH2F* h)
{
  Int_t binx = h->GetXaxis()->FindBin(x);
  Int_t biny = h->GetYaxis()->FindBin(y);
  return h->GetBinContent(binx, biny);

}

float kfactorSUSY(float m0, float m12, string sample){

 char* filename;
 float kfactor = 1.0; 
 char *flavor[10] = {(char*)"ng", (char*)"ns", (char*)"nn", (char*)"ll", (char*)"sb", (char*)"ss", (char*)"tb", (char*)"bb", (char*)"gg", (char*)"sg"};

 if      ( sample == "tanbeta3"              ) filename = (char*)"../data/nlotanbeta3.root";
 else if ( sample == "tanbeta3Scale05"       ) filename = (char*)"../data/nlotanbeta305.root"; 
 else if ( sample == "tanbeta3Scale20"       ) filename = (char*)"../data/nlotanbeta320.root"; 
 else if ( sample == "tanbeta10"             ) filename = (char*)"../data/nlotanbeta10.root";
 else if ( sample == "tanbeta10Scale05"      ) filename = (char*)"../data/nlotanbeta1005.root";
 else if ( sample == "tanbeta10Scale20"      ) filename = (char*)"../data/nlotanbeta1020.root";
 else if ( sample == "tanbeta10_2012"        ) filename = (char*)"../data/msugra_tanb10_2012_kfactor_nominal.root";
 else if ( sample == "tanbeta10Scale05_2012" ) filename = (char*)"../data/msugra_tanb10_2012_kfactor_scale05.root";
 else if ( sample == "tanbeta10Scale20_2012" ) filename = (char*)"../data/msugra_tanb10_2012_kfactor_scale20.root";
 else if ( sample == "tanbeta10_2012final"   ) filename = (char*)"../data/combined_kfactor_msugra_m0_m12_10_0_1.root";
 else if ( sample == "tanbeta10Up_2012final" ) filename = (char*)"../data/combined_kfactorUp_msugra_m0_m12_10_0_1.root";
 else if ( sample == "tanbeta10Dn_2012final" ) filename = (char*)"../data/combined_kfactorDn_msugra_m0_m12_10_0_1.root";
 else {
  return kfactor;
  cout << "WARNING:: CORE/kfactorSUSY wrong string given, using kfactor=1" << endl; 
 }

  TFile infile(filename);
  if (infile.IsZombie()) {
    cout << "Error opening file " << filename << endl;
    gSystem->Exit(1);
  }

// determine the process

  TDatabasePDG *pdg = new TDatabasePDG(); 
  std::vector<int> interactions; 
 
   for (unsigned int j=0; j<cms2.genps_id().size(); j++) { 
    if (cms2.genps_status().at(j) != 3) continue; 
    int ID = abs(cms2.genps_id().at(j)); 
    int mID = abs(cms2.genps_id_mother().at(j)); 
    if (ID > 1000000 && ID < 2000016 && (mID < 7 || mID ==21 || mID == 22 || mID == 23 || mID == 24)) { // kept the bosons in case of screw ups 
 
     // Check the mother of the SM is a proton 
       bool isProton = false;        
       for (unsigned int k=0; k < j; k++) { 
            int pID = abs(cms2.genps_id().at(k)); 
            int mpID = abs(cms2.genps_id_mother().at(k)); 
         if ((pID < 7 || pID ==21 || pID == 22 || pID == 23 || pID == 24)&&(mpID == 2212)) isProton = true;  
       } 
//      if (!isProton) continue;   
      if (!isProton) { 
          dumpDocLines(); 
          continue; 
       } 
      interactions.push_back(cms2.genps_id().at(j)); 
 
    if (interactions.size() > 2)  
         cout << setw(4) << left << j << " WARNING mcSUSYkfactor: Something is wrong with " 
         << setw(10) << left << pdg->GetParticle(cms2.genps_id().at(j))->GetName() << " " 
         << setw(10) << left << cms2.genps_id().at(j) << " " 
         << setw(7) << right << setprecision(4) << cms2.genps_p4().at(j).pt() << "  " 
         << setw(7) << right << setprecision(4) << cms2.genps_p4().at(j).phi() << "  " 
         << setw(10) << right << setprecision(4) << cms2.genps_p4().at(j).eta() << "  " 
         << setw(4) << right << cms2.genps_status().at(j) << " " 
         << setw(10) << left << pdg->GetParticle(cms2.genps_id_mother().at(j))->GetName() 
         << " using k=1 " << endl; 
   } 
  } 
 
  delete pdg; 
 
  if (interactions.size() == 2 ) { 
      int index = sfinalState(interactions[0], interactions[1]); 
      if (index < 0) return 1.0; 
      TH2F* thisHist = dynamic_cast<TH2F *>( infile.Get( Form("h%s", flavor[index])) );
      if ( thisHist == 0) {
        cout << "Cannot find histogram in file " << filename << endl;
        gSystem->Exit(1);
      }
      kfactor = GetValueTH2FS(m0, m12, thisHist);
    delete thisHist;
    infile.Close();

   } else { 
    cout << "WARNING mcSUSYkfactor: Number of sparticles found: " << interactions.size() << " using kfactor=1" <<endl; 
    kfactor = 1.0; 
   } 
 
  return kfactor; 
}

/*
Float_t GetValueTH2FS(Float_t x, Float_t y, TH2F* h) 
{ 
  Int_t binx = h->GetXaxis()->FindBin(x);
  Int_t biny = h->GetYaxis()->FindBin(y);
  return h->GetBinContent(binx, biny);
 
} 
*/

float cmssm_loxsec(float m0, float m12, string filename){

  //string filename = "../data/loxsec_m02TeV.root";

  TFile infile(filename.c_str());
  if (infile.IsZombie()) {
    cout << "Error opening file " << filename << endl;
    gSystem->Exit(1);
  }

  TH2F* hist = (TH2F*) infile.Get("xsec");
  float xsec = GetValueTH2FS(m0,m12,hist);

  infile.Close();

  return xsec;
}

float kfactorSUSY(string sample)
{
  float kfactor = 1.0;
  TDatabasePDG *pdg = new TDatabasePDG();
  std::vector<int> interactions;

   for (unsigned int j=0; j<cms2.genps_id().size(); j++) {
    if (cms2.genps_status().at(j) != 3) continue;
    int ID = abs(cms2.genps_id().at(j));
    int mID = abs(cms2.genps_id_mother().at(j));
    if (ID > 1000000 && ID < 2000016 && (mID < 7 || mID ==21 || mID == 22 || mID == 23 || mID == 24)) { // kept the bosons in case of screw ups

     // Check the mother of the SM is a proton
       bool isProton = false;       
       for (unsigned int k=0; k < j; k++) {
            int pID = abs(cms2.genps_id().at(k));
            int mpID = abs(cms2.genps_id_mother().at(k));
         if ((pID < 7 || pID ==21 || pID == 22 || pID == 23 || pID == 24)&&(mpID == 2212)) isProton = true; 
       }
//      if (!isProton) continue;  
      if (!isProton) {
          dumpDocLines();
          continue;
       }
      interactions.push_back(cms2.genps_id().at(j));

    if (interactions.size() > 2) 
         cout << setw(4) << left << j << " WARNING mcSUSYkfactor: Something is wrong with "
         << setw(10) << left << pdg->GetParticle(cms2.genps_id().at(j))->GetName() << " "
         << setw(10) << left << cms2.genps_id().at(j) << " "
         << setw(7) << right << setprecision(4) << cms2.genps_p4().at(j).pt() << "  "
         << setw(7) << right << setprecision(4) << cms2.genps_p4().at(j).phi() << "  "
         << setw(10) << right << setprecision(4) << cms2.genps_p4().at(j).eta() << "  "
         << setw(4) << right << cms2.genps_status().at(j) << " "
         << setw(10) << left << pdg->GetParticle(cms2.genps_id_mother().at(j))->GetName()
         << " using k=1 " << endl;
   }
  }

  delete pdg;

  if (interactions.size() == 2 ) {
     kfactor = lmdata(interactions[0], interactions[1], sample);
   } else {
    cout << "WARNING mcSUSYkfactor: Number of sparticles found: " << interactions.size() << " using kfactor=1" <<endl;
    kfactor = 1.0;
//    dumpDocLines();
   }

  return kfactor;
}

double lmdata(int ipart1, int ipart2, string prefix) {

  static const Double_t lm0[10] = {1.06604, 1.00369, 1.27186, 1.19103, 1.44681, 1.22883, 1.5649, 1.70195, 1.99721, 1.33951};
  static const Double_t lm0scale[10] = {1.04522, 1.0001, 1.27058, 1.19107, 1.39667, 1.18765, 1.53141, 1.66519, 1.88761, 1.27273 };

  static const Double_t lm1[10] = {0.988762, 1.02694, 1.22581, 1.24139, 1.45174, 1.18939, 1.72769, 1.7931, 2.33333, 1.39205};
  static const Double_t lm2[10] = {0.971508, 1.07009, 1.18325, 1.1879, 1.47269, 1.16915, 1.86207, 1.94417, 2.96117, 1.57895};
  static const Double_t lm3[10] = {1.03489, 1.03535, 1.2383, 1.12839, 1.42798, 1.20694, 1.74914, 1.83868, 2.46354, 1.46457};
  static const Double_t lm4[10] = {0.978699, 1.04381, 1.21136, 1.16334, 1.47239, 1.18924, 1.779, 1.87219, 2.54733, 1.47039};
  static const Double_t lm5[10] = {0.969489, 1.07524, 1.17991, 1.14481, 1.48387, 1.16456, 1.88094, 2.00574, 3.05195, 1.6087};
  static const Double_t lm6[10] = {0.971497, 1.08695, 1.16543, 1.1852, 1.47826, 1.15094, 1.9199, 2.04659, 3.32237, 1.66159};
  static const Double_t lm8[10] = {1.02828, 1.07339, 1.2137, 1.08113, 1.43871, 1.20667, 1.82724, 1.97961, 2.91919, 1.63415};
  static const Double_t lm9[10] = {1.08102, 1.21649, 3.31807, 1.09983, 2.16036, 1.30526, 2.18081, 2.2218, 2.17647, 2.14039};

 
  int index = sfinalState(ipart1, ipart2);
  if (index < 0) return 1.0;

  if (prefix == "lm0") return lm0[index];
  else if (prefix == "lm1") return lm1[index];
  else if (prefix == "lm2") return lm2[index];
  else if (prefix == "lm3") return lm3[index];
  else if (prefix == "lm4") return lm4[index];
  else if (prefix == "lm5") return lm5[index];
  else if (prefix == "lm6") return lm6[index];
  else if (prefix == "lm8") return lm8[index];
  else if (prefix == "lm9") return lm9[index];
  else if (prefix == "lm0scale") return lm0scale[index];

  return 1.0;
}

