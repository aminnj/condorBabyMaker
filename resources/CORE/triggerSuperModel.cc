#include "triggerSuperModel.h"
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "CMS2.h"

//----------------------------------------
// To be used on Monte Carlo events.
// Input is hypothesis index.
// Output is the 2010 trigger efficiency for the 
// given hypothesis.
//
// Claudio, Jae, Avi  3 Nov 2010
//----------------------------------------

float triggerSuperModelEffic(int hyp) {

  //----------------------------------------------
  // Some important inputs
  //-----------------------------------------------
  // Plateau efficiency of muon trigger
  // (For now we ignore the fact that in the mu9 period 
  //  it was not quite as good)
  float effmu   = 0.93;   // eta<2.1  
  float effmu24 = 0.40;   // 2.1<eta<2.4

  // Fraction of luminosity where mu9 was not prescaled
  // run <= 147116 
  float f9=0.215;

  // Fraction of luminosity where mu9 was prescaled and mu11
  // was unprescaled
  // 147196 <= run <= 148058
  float f11=0.273;

  // Fraction of luminosity where mu9 and mu11 were prescaled 
  // and mu15 was not
  // 148058 <= run <= 149442
  //float f15 = 1.0 - f9 - f11;

  // Fraction of luminosity with 100% efficient ele10 trigger
  // 136033 <= run <= 139980
  float e10=0.002;

  //Fraction of luminosity with 100% efficienct ele15 GeV trigger
  // 139980<run<=144114
  float e15=0.086;

  //Fraction of luminosity with 100% efficient ele17 trigger
  // 144114<run<=147116
  float e17=0.127;

  //Fraction of luminosity with somewhat inefficient (~93%) ele17 trigger
  // 147116<run<=148058
  float e17b=0.273;

  //Efficiency of 17 GeV trigger in 147116<run<=148058
  float eff17b=0.93; 

  //Fraction of luminosity with messy electron trigger
  // run>148058
  float emess=0.512;

  //Efficiency in pt bins for messy electron trigger period
  // run>148058  (rounded to the nearest 1%)
  float eff17to22 =0.99;
  float eff22to27 =0.97;
  float eff27to32 =0.98;
  float eff32andUp=1.00;
    
  //-------------------------------------------------
  // Algorithm applies only to events with one
  // lepton above 20 and one letpon above 10.
  // If this is not the case, return zero and complain.
  //-------------------------------------------------
  bool badHyp = false;
  if( TMath::Max(cms2.hyp_ll_p4()[hyp].pt(),cms2.hyp_lt_p4()[hyp].pt()) < 20.)
    badHyp=true;
  if( TMath::Min(cms2.hyp_ll_p4()[hyp].pt(),cms2.hyp_lt_p4()[hyp].pt()) < 10.)
    badHyp=true;
  if (badHyp) {
    std::cout << "Bad inputs to trigger SuperModelEffic" << std::endl;
    return 0.0;
  }

  //-------------------------------------------------
  // If it is a dielectron event, we return 100%
  //-------------------------------------------------
  if (cms2.hyp_type()[hyp] == 3) return 1.0;

  //-------------------------------------------------
  // Dimuon events.  
  //-------------------------------------------------

  float pt1;
  float pt2;
  float eta1;
  float eta2;

  // Convention: 1 is the highest pt guy, 2 is the lowest pt guy
  if (cms2.hyp_type()[hyp] == 0) {
    pt1 = cms2.hyp_ll_p4()[hyp].pt();
    pt2 = cms2.hyp_lt_p4()[hyp].pt();
    eta1= TMath::Abs(cms2.hyp_ll_p4()[hyp].eta());
    eta2= TMath::Abs(cms2.hyp_lt_p4()[hyp].eta());
    if (pt2 > pt1) {
      pt1 = cms2.hyp_lt_p4()[hyp].pt();
      pt2 = cms2.hyp_ll_p4()[hyp].pt();
      eta1= TMath::Abs(cms2.hyp_lt_p4()[hyp].eta());
      eta2= TMath::Abs(cms2.hyp_ll_p4()[hyp].eta());
    }

    // both above 15 and both in eta<2.1
    if (pt2>=15 && eta1<=2.1 && eta2<=2.1) {
      float eff = 1 - (1-effmu)*(1-effmu);
      return eff;
    }

    // the 2nd one between 11 and 15, both in eta<2.1
    if (pt2>11 && pt2<15 && eta1<=2.1 && eta2<=2.1) {
      float eff = effmu + (f9+f11)*effmu*(1-effmu);
      return eff;
    }

    // the 2nd one between 10 and 11, both in eta<2.1
    if (pt2<11 && eta1<=2.1 && eta2<=2.1) {
      float eff = effmu + f9*effmu*(1-effmu);
      return eff;
    }

    // both at high eta
    if (eta1>2.1 && eta2>2.1) {
      float eff = effmu*effmu;  // this neglects the trigger at high eta
      if (pt2>=15) {
	eff = eff + 2*(1.-effmu)*effmu24;
      } else if (pt2>=11) {
	eff = eff + (1.+f9+f11)*(1-effmu)*effmu24;
      } else {
	eff = eff + (1.+f9)*(1-effmu)*effmu24;
      }
      return eff;
    }

    // One with pt>15 eta<2.1.  The other with eta>2.1
    if ( (pt1>=15 && eta1<=2.1 && eta2>2.1) || 
         (pt2>=15 && eta2<=2.1 && eta1>2.1)    ) {
      float eff = effmu;   // this neglects the trigger at high eta
      if ( (pt1>=15 && eta1>2.1) || (pt2>=15 && eta2>2.1) ) {
	eff = eff + (1-effmu)*effmu24;
      } else if ( (pt1>=11 && eta1>2.1) || (pt2>=11 && eta2>2.1) ) {
	eff = eff + (f9+f11)*(1-effmu)*effmu24;
      } else {
	eff = eff + f9*(1-effmu)*effmu24;
      }
      return eff;
    }
    

    // First with 11<pt<15 eta<2.1.  Second one with eta>2.1
    if ( (pt1>=11 && pt1<15 && eta1<=2.1 && eta2>2.1) || 
	 (pt2>=11 && pt2<15 && eta2<=2.1 && eta1>2.1)   ) {
      float eff = (f9+f11)*effmu + (1-f9-f11)*effmu*effmu; 
      eff = eff + (1-effmu)*effmu24; // allow trigger at high eta
      return eff;
    }

    // First with 10<pt<11 eta<2.1.  Second one with eta>2.1
    if ( (pt1<11 && eta1<=2.1 && eta2>2.1) || 
	 (pt2<11 && eta2<=2.1 && eta1>2.1)   ) {
      float eff = f9*effmu + (1-f9)*effmu*effmu;
      eff = eff + (1-effmu)*effmu24;  // allow trigger at high eta
      return eff;
    }

    // We should never get here!
    std::cout << "----------" <<std::endl;
    std::cout << "Logic failure for mu-mu events in triggerSuperModel" << std::endl;
    std::cout << "This should never happen -- do not ignore" << std::endl;
    std::cout << "----------" << std::endl;
    return 0.0;

  } // close mumu code block

  //-------------------------------------------------
  // emu events
  //-------------------------------------------------
  if (cms2.hyp_type()[hyp] == 1 || cms2.hyp_type()[hyp] == 2) {

    float ptmu;
    float ptele;
    float etamu;

    if (TMath::Abs(cms2.hyp_ll_id()[hyp]) == 13) {      
      ptmu  = cms2.hyp_ll_p4()[hyp].pt();
      ptele = cms2.hyp_lt_p4()[hyp].pt();
      etamu=  TMath::Abs(cms2.hyp_ll_p4()[hyp].eta());
    } else {
      ptmu  = cms2.hyp_lt_p4()[hyp].pt();
      ptele = cms2.hyp_ll_p4()[hyp].pt();
      etamu=  TMath::Abs(cms2.hyp_lt_p4()[hyp].eta());
    }


    // muon in eta<2.1 and pt>15;
    if (ptmu>= 15 && etamu<= 2.1) {      
      float delta;
      if (ptele<=15.) {
	delta=e10;
      } else if (ptele<17) {
	  delta=e10+e15;
      } else if (ptele>=17) {
	delta = e10+e15+e17+e17b*eff17b;
      } if (ptele<=22) {
	delta = e10+e15+e17+e17b*eff17b+emess*eff17to22;
      } else if (ptele<=27) {
	delta = e10+e15+e17+e17b*eff17b+emess*eff22to27;
      } else if (ptele<=32) {
	delta = e10+e15+e17+e17b*eff17b+emess*eff27to32;
      } else {
	delta = e10+e15+e17+e17b*eff17b+emess*eff32andUp;
      }
      float eff = effmu + (1-effmu)*delta;
      return eff;
    }

    // muon in eta<2.1 and pt beween 11 and 15 or 10 and 11
    // (note: here the electron has pt>20 for sure)
    if (ptmu< 15 && etamu<= 2.1) {      
      float eleff;
      float f=f9; 
      if (ptmu>=11)  f=f11+f9;
      if (ptele <= 22) eleff=eff17to22;
      if (ptele <= 27) eleff=eff22to27;
      if (ptele <= 32) eleff=eff27to32;
      if (ptele >  32) eleff=eff32andUp;
      float delta2 = (1-f)*(effmu + (1-effmu)*eleff);
      float delta3 = f*(1-effmu)*( (e10+e15+e17)/f + eff17b*(f-e10-e15-e17)/f);
      float eff = f*effmu + delta2 +delta3;
      return eff;
    }

    // muon in eta>2.1
    if (etamu>2.1) {
      float eff;
      if (ptele<=15) {
	eff = e10;
      } else if (ptele<=17) {
	eff = e10 + e15;
      } else if (ptele <= 22) {
	eff = e10 + e15 + e17 + e17b*eff17b + emess*eff17to22;
      } else if (ptele <= 27) {
	eff = e10 + e15 + e17 + e17b*eff17b + emess*eff22to27;
      } else if (ptele <= 32) {
	eff = e10 + e15 + e17 + e17b*eff17b + emess*eff27to32;
      } else {
	eff = e10 + e15 + e17 + e17b*eff17b + emess*eff32andUp;
      }
      // now allow for the possibility of a muon trigger
      if (ptmu>=15) {
	eff = eff + (1-eff)*effmu24;
      } else if (ptmu>=11) {
	eff = eff + (f9+f11)*(1-eff)*effmu24;
      } else {
	eff = eff + f9*(1-eff)*effmu24;
      }
      return eff;
    }

    // We should never get here!
    std::cout << "----------" <<std::endl;
    std::cout << "Logic failure for e-mu events in triggerSuperModel" << std::endl;
    std::cout << "This should never happen -- do not ignore" << std::endl;
    std::cout << "----------" << std::endl;
    return 0.0;

  } // Close emu code block

    // We should never get here!
    std::cout << "----------" <<std::endl;
    std::cout << "Logic failure in triggerSuperModel" << std::endl;
    std::cout << "This should never happen -- do not ignore" << std::endl;
    std::cout << "----------" << std::endl;
    return 0.0;

}
