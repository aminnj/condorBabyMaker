
#include <assert.h>
#include <iostream>
#include <algorithm>

//core headers
#include "ttvSelections.h"
#include "electronSelections.h"
#include "electronSelectionsParameters.h"
#include "muonSelections.h"
#include "metSelections.h"
#include "ssSelections.h"
#include "jetSelections.h"
#include "triggerUtils.h"
#include "eventSelections.h"
#include "utilities.h"
#include "susySelections.h"
#include "jetcorr/FactorizedJetCorrector.h"

//Root headers
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"

using namespace ttv;

struct jet_pt_gt_ttv {
  bool operator () (const LorentzVector &v1, const LorentzVector &v2) 
  {
	return v1.pt() > v2.pt();
  }
};

bool ttv::passesTrigger (int hyp_type)
{
  //----------------------------------------
  // no trigger requirements applied to MC
  //----------------------------------------
  if (!cms2.evt_isRealData())
  {
	  return true; 
  }

  //---------------------------------
  // triggers for dilepton datasets
  //---------------------------------
  //mm
  if (hyp_type == 0) 
  {
	  if( passUnprescaledHLTTriggerPattern("HLT_Mu17_Mu8_v" ) )   { return true; }
  }
  //em
  else if ((hyp_type == 1 || hyp_type == 2)) 
  {
	  if( passUnprescaledHLTTriggerPattern("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") )  { return true; }
	  if( passUnprescaledHLTTriggerPattern("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") )  { return true; }
  }
  //ee
  else if (hyp_type == 3) 
  {
	  if( passUnprescaledHLTTriggerPattern("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") ) { return true; }
  }
    
  return false;
}

bool ttv::passesSingleLeptonTrigger (TrileptonInfo& trilep_info)
{
  return passesSingleLeptonTrigger( trilep_info.z.lep1.id, trilep_info.z.lep2.id, trilep_info.w.id );
}
bool ttv::passesSingleLeptonTrigger (int id1, int id2, int id3)
{
  if (!cms2.evt_isRealData())
  {
	  return true; 
  }
  //mxx
  if(abs(id1) == 13 || abs(id2) == 13 || abs(id3) == 13){
	if( passUnprescaledHLTTriggerPattern("HLT_IsoMu24_eta2p1_v") ) { return true; }
	if( passUnprescaledHLTTriggerPattern("HLT_IsoMu24_v")        ) { return true; }
  }
  return false;
}

bool ttv::passesDiLeptonTrigger (TrileptonInfo& trilep_info)
{
  return passesDiLeptonTrigger( trilep_info.z.lep1.id, trilep_info.z.lep2.id, trilep_info.w.id );
}

bool ttv::passesDiLeptonTrigger (int id1, int id2, int id3)
{
  if (!cms2.evt_isRealData())
  {
	  return true; 
  }

  //mmx
  if (abs(id1*id2) == 13*13 || abs(id1*id3) == 13*13 || abs(id2*id3) == 13*13) 
  {
	  if( passUnprescaledHLTTriggerPattern("HLT_Mu17_Mu8_v")         ) { return true; }
	  if( passUnprescaledHLTTriggerPattern("HLT_Mu17_TkMu8_v")       ) { return true; }
  }   
  //eex
  if (abs(id1*id2) == 11*11 || abs(id1*id3) == 11*11 || abs(id2*id3) == 11*11)
  {
	  if( passUnprescaledHLTTriggerPattern("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") ) { return true; }
	  if( passUnprescaledHLTTriggerPattern("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")                                         ) { return true; }
  }
  //emx
  if (abs(id1*id2) == 11*13 || abs(id1*id3) == 11*13 || abs(id2*id3) == 11*13) 
  {
	  if( passUnprescaledHLTTriggerPattern("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") )  { return true; }
	  if( passUnprescaledHLTTriggerPattern("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") )  { return true; }
  }
  return false;
}

bool ttv::passesTriLeptonTrigger (TrileptonInfo& trilep_info)
{
  return passesTriLeptonTrigger( trilep_info.z.lep1.id, trilep_info.z.lep2.id, trilep_info.w.id );
}
bool ttv::passesTriLeptonTrigger (int id1, int id2, int id3)
{
  if (!cms2.evt_isRealData())
  {
	  return true; 
  }

  //mmm
  if (abs(id1*id2*id3) == 13*13*13) 
  {
	  if( passUnprescaledHLTTriggerPattern("HLT_DoubleMu5_IsoMu5_v") ) { return true; }
	  if( passUnprescaledHLTTriggerPattern("HLT_TripleMu5_v")        ) { return true; }
  }
  //eee
  if (abs(id1*id2*id3) == 11*11*11) 
  {
	  if( passUnprescaledHLTTriggerPattern("HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_CaloIdT_TrkIdVL_v") ) { return true; }
	  if( passUnprescaledHLTTriggerPattern("HLT_Ele15_Ele8_Ele5_CaloIdL_TrkIdVL_v")                   ) { return true; }
  }
  //mme
  if (abs(id1*id2*id3) == 11*13*13) 
  {
	  if( passUnprescaledHLTTriggerPattern("HLT_DoubleMu5_Ele8_CaloIdT_TrkIdVL_v")                ){ return true; }
	  if( passUnprescaledHLTTriggerPattern("HLT_DoubleMu8_Ele8_CaloIdT_TrkIdVL_v")                ){ return true; }
  }
  //eem
  if (abs(id1*id2*id3) == 11*11*13)
  {
	  if( passUnprescaledHLTTriggerPattern("HLT_Mu8_DoubleEle8_CaloIdT_TrkIdVL_v")                ){ return true; }
	  if( passUnprescaledHLTTriggerPattern("HLT_Mu8_Ele8_CaloIdT_TrkIdVL_Ele8_CaloIdL_TrkIdVL_v") ){ return true; }
  }
  return false;
}

bool ttv::passesTrigger (TrileptonInfo& trilep_info)
{
  return passesTrigger( trilep_info.z.lep1.id, trilep_info.z.lep2.id, trilep_info.w.id );
}

bool ttv::passesTrigger (int id1, int id2, int id3)
{
  //----------------------------------------
  // no trigger requirements applied to MC
  //----------------------------------------  
  if (!cms2.evt_isRealData())
  {
	  return true; 
  }

  //---------------------------------
  // triggers for [mono,di,tri]lepton datasets
  //---------------------------------
  if( passesSingleLeptonTrigger(id1, id2, id3) ){ return true; }
  if( passesDiLeptonTrigger(id1, id2, id3) ){ return true; }
  if( passesTriLeptonTrigger(id1, id2, id3) ){ return true; }

  return false;
}

bool ttv::isGoodLepton (int id, int idx, LeptonType::value_type lep_type)
{
  if (abs(id) == 11)
  {
	  if (lep_type == ttv::LeptonType::LOOSE)
      {
		  if( fabs(cms2.els_etaSC()[idx]) > 1.4442 && fabs(cms2.els_etaSC()[idx]) < 1.566 )     { return false; }
		  if( fabs(cms2.els_p4()[idx].eta()) > 2.5 )                                            { return false; }
		  if( ttv::overlapMuon(idx, ttv::LeptonType::LOOSE, (float)20.0) )                      { return false; }
		  electronIdComponent_t answer_loose_2012 = electronId_WP2012_v2(idx, LOOSE);
		  if ((answer_loose_2012 & wp2012::PassWP2012CutsNoIso) != wp2012::PassWP2012CutsNoIso) { return false; } 
		  return true;
      }
	  else if (lep_type == ttv::LeptonType::TIGHT)
      {
		  if( ttv::overlapMuon(idx, ttv::LeptonType::LOOSE, 20.0) )                              { return false; }
		  return pass_electronSelection(idx, electronSelection_TTVTightv1_noIso);
      }
	  else if (lep_type == ttv::LeptonType::LOOSEDILEPMVA || lep_type == ttv::LeptonType::TIGHTDILEPMVA)
      {
		  if( fabs(cms2.els_etaSC()[idx]) > 1.4442 && fabs(cms2.els_etaSC()[idx]) < 1.566 )      { return false; }
		  if( fabs(cms2.els_p4()[idx].eta()) > 2.5 )                                             { return false; }
		  if( ttv::overlapMuon(idx, ttv::LeptonType::LOOSE, (float)20.0) )                       { return false; }
		  electronIdComponent_t answer_medium_2012 = electronId_WP2012_v2(idx, MEDIUM);
		  electronIdComponent_t cutmask = wp2012::MHITS | wp2012::VTXFIT | wp2012::D0VTX | wp2012::DZVTX;
		  if( (answer_medium_2012 & cutmask) != cutmask )                                        { return false; }
		  electronIdComponent_t answer_loose_2012 = electronId_WP2012_v2(idx, LOOSE);
		  electronIdComponent_t fiducialMask = wp2012::DETAIN | wp2012::DPHIIN | wp2012::SIGMAIETAIETA | wp2012::HOE;
		  if( (answer_loose_2012 & fiducialMask) != fiducialMask )                               { return false; } 
		  if (abs(cms2.els_etaSC().at(idx)) < 1.479) 
		  {
			if( cms2.els_tkIso().at(idx)/cms2.els_p4().at(idx).pt() > 0.2 )                      { return false; }
			if( max((cms2.els_ecalIso().at(idx) - 1.0)/cms2.els_p4().at(idx).pt() , 0. ) > 0.2 ) { return false; }
			if( cms2.els_hcalIso().at(idx)/cms2.els_p4().at(idx).pt() > 0.2 )                    { return false; } 
		  }
		  else
		  {
			if( cms2.els_tkIso().at(idx)/cms2.els_p4().at(idx).pt() > 0.2 )                      { return false; }			
			if( cms2.els_ecalIso().at(idx)/cms2.els_p4().at(idx).pt() > 0.2 )                    { return false; } 		
			if( cms2.els_hcalIso().at(idx)/cms2.els_p4().at(idx).pt() > 0.2 )                    { return false; } 
		  }
		  return true;
      }
	  else if (lep_type == ttv::LeptonType::LOOSETRILEPMVA || lep_type == ttv::LeptonType::TIGHTTRILEPMVA)
      {
		  if( fabs(cms2.els_etaSC()[idx]) > 1.4442 && fabs(cms2.els_etaSC()[idx]) < 1.566 )      { return false; }
		  if( fabs(cms2.els_p4()[idx].eta()) > 2.5 )                                             { return false; }
		  if( ttv::overlapMuon(idx, ttv::LeptonType::LOOSE, (float)20.0) )                       { return false; }
		  electronIdComponent_t answer_medium_2012 = electronId_WP2012_v2(idx, MEDIUM);
		  electronIdComponent_t cutmask = wp2012::MHITS | wp2012::VTXFIT | wp2012::D0VTX | wp2012::DZVTX;
		  if( (answer_medium_2012 & cutmask) != cutmask )                                        { return false; } 
		  electronIdComponent_t answer_loose_2012 = electronId_WP2012_v2(idx, LOOSE);
		  electronIdComponent_t fiducialMask = wp2012::DETAIN | wp2012::DPHIIN | wp2012::SIGMAIETAIETA | wp2012::HOE;
		  if( (answer_loose_2012 & fiducialMask) != fiducialMask )                               { return false; }
		  return true;
      }
	  else
      {
		  std::cout<<"Warning: in ttvSelections bool ttv::isGoodLepton, unsupported enum "<<lep_type<<" for electrons."<<std::endl;
      }
  }
  else if (abs(id) == 13)
  {
	  if (lep_type == ttv::LeptonType::LOOSE)
      {
		  return (muonIdNotIsolated(idx, NominalTTZ_loose_v1));
      }
	  else if (lep_type == ttv::LeptonType::TIGHT)
      {
		  return (muonIdNotIsolated(idx, NominalTTZ_tight_v1));
      }
	  else
      {
		  std::cout<<"Warning: in ttvSelections bool ttv::isGoodLepton, unsupported enum "<<lep_type<<" for muons."<<endl;
      }
  }

  return false;
}

bool ttv::isGoodLepton (int id, int idx, LeptonType::value_type lep_type, float mvaValue)
{
  if( ! ttv::isGoodLepton (id, idx, lep_type) )              { return false; }
  if( !(mvaValue > ttv::getTrigMVAThreshold(idx, lep_type)) ){ return false; }

  return true;
}

bool ttv::isIsolatedLepton (int id, int idx, LeptonType::value_type lep_type)
{
  if (abs(id) == 11)
  {
	  if (lep_type == ttv::LeptonType::LOOSE || lep_type == ttv::LeptonType::LOOSEDILEPMVA || lep_type == ttv::LeptonType::LOOSETRILEPMVA)
      {
		  return (electronIsoValuePF2012_FastJetEffArea_v2( idx ) < 0.15);
      }
	  else if (lep_type == ttv::LeptonType::TIGHT || lep_type == ttv::LeptonType::TIGHTDILEPMVA || lep_type == ttv::LeptonType::TIGHTTRILEPMVA)
      {
		  return (electronIsoValuePF2012_FastJetEffArea_v2( idx ) < 0.09);
      }
	  else
	  {
		  cout<<"Warning: in ttvSelections bool ttv::isIsolatedLepton, unsupported enum "<<lep_type<<" for electrons."<<endl;
	  }
  }
  else if (abs(id) == 13)
  {
	  if (lep_type == ttv::LeptonType::LOOSE)
      {
		  return (muonIsoValuePF2012_deltaBeta(idx) < 0.15);
      }
	  else if (lep_type == ttv::LeptonType::TIGHT)
      {
		  return (muonIsoValuePF2012_deltaBeta(idx) < 0.10);
      }
	  else
	  {
		  cout<<"Warning: in ttvSelections bool ttv::isIsolatedLepton, unsupported enum "<<lep_type<<" for muons."<<endl;
	  }
  }

  return false;
}

bool ttv::isNumeratorLepton (int id, int idx, LeptonType::value_type lep_type)
{
  return (ttv::isGoodLepton(id, idx, lep_type) && ttv::isIsolatedLepton(id, idx, lep_type));
}

bool ttv::isNumeratorLepton (int id, int idx, LeptonType::value_type lep_type, float mvaValue)
{
  if(! ttv::isGoodLepton(id, idx, lep_type, mvaValue) ){ return false; }
  if(! ttv::isIsolatedLepton(id, idx, lep_type) )      { return false; }	
  return true;
}

bool ttv::isDenominatorLepton (int id, int idx, LeptonType::value_type lep_type)
{
  if (abs(id) == 11)
  {
	  if (lep_type == ttv::LeptonType::LOOSE)
      {
		  if( fabs(cms2.els_etaSC()[idx]) > 1.4442 && fabs(cms2.els_etaSC()[idx]) < 1.566 )             { return false; }
		  if( fabs(cms2.els_p4()[idx].eta()) > 2.5 )                                                    { return false; }
		  if( ttv::overlapMuon(idx, ttv::LeptonType::LOOSE, (float)20.0) )                              { return false; }
		  electronIdComponent_t answer_loose_2012 = electronId_WP2012_v2(idx, LOOSE);
		  if ((answer_loose_2012 & wp2012::PassWP2012CutsNoIsoNoIP) != wp2012::PassWP2012CutsNoIsoNoIP) { return false; }  // if we add the d0 back in, change from NoIsoNoIP to NoIso
		  return true;
      }
	  else if (lep_type == ttv::LeptonType::TIGHT)
      {
		  if( ttv::overlapMuon(idx, ttv::LeptonType::LOOSE, 20.0) )                                     { return false; }
		  return (pass_electronSelection(idx, electronSelection_TTVTightFOv1));
      }
  }
  else if (abs(id) == 13)
  {
	  if (lep_type == ttv::LeptonType::LOOSE)
      {
		  return (muonId(idx, NominalTTZ_looseFO_v1));
      }
	  else if (lep_type == ttv::LeptonType::TIGHT)
      {
		  return (muonId(idx, NominalTTZ_tightFO_v1));
      }
  }

  return false;
}

bool ttv::isDenominatorLepton (int id, int idx, LeptonType::value_type lep_type, float mvaValue)
{
  if( ! ttv::isDenominatorLepton (id, idx, lep_type) )       { return false; }
  if( !(mvaValue > ttv::getTrigMVAThreshold(idx, lep_type)) ){ return false; }
  return true;
}

bool ttv::overlapMuon(int idx, LeptonType::value_type lep_type, float pt, float eta, float deltaR)
{
  for( unsigned int muidx = 0; muidx < cms2.mus_p4().size(); muidx++)
  {
	  if (cms2.mus_p4()[muidx].pt() > pt &&
		  abs(cms2.mus_p4()[muidx].eta()) < eta &&
		  ttv::isNumeratorLepton(13, muidx, lep_type) &&
		  ROOT::Math::VectorUtil::DeltaR(cms2.mus_p4()[muidx], cms2.els_p4()[idx]) < deltaR)
	  { 
		  return true; 
	  }
  }
  return false;
}

std::vector<LorentzVector> ttv::getGoodLeptons(LeptonType::value_type lep_type, float pt, float mu_eta, float ele_eta)
{
  std::vector<LorentzVector> good_leps;

  for (unsigned int midx = 0; midx < cms2.mus_p4().size(); midx++)
  {
	  if (cms2.mus_p4().at(midx).pt() < pt)            { continue; }
	  if (fabs(cms2.mus_p4().at(midx).eta()) > mu_eta) { continue; }
	  if (!isGoodLepton(13, midx, lep_type))           { continue; }
	  good_leps.push_back(cms2.mus_p4().at(midx));
  }

  for (unsigned int eidx = 0; eidx < cms2.els_p4().size(); eidx++)
  {
	  if (cms2.els_p4().at(eidx).pt() < pt)             { continue; }
	  if (fabs(cms2.els_p4().at(eidx).eta()) > ele_eta) { continue; }
	  if (!isGoodLepton(11, eidx, lep_type))            { continue; }
	  good_leps.push_back(cms2.els_p4().at(eidx));
  }

  return good_leps;
}

float ttv::getTrigMVAThreshold(int idx, LeptonType::value_type lep_type)
{
  float looselowptValues[3] = { -0.756, -0.396, -0.276 };
  float looseValues[3]      = { 0.949,  0.921,  0.956  };
  float tightlowptValues[3] = { -0.524, 0.204,  0.292  };
  float tightValues[3]      = { 0.956,  0.949,  0.968  };

  if( lep_type == LeptonType::TIGHTDILEPMVA || lep_type == LeptonType::TIGHTTRILEPMVA )
  {
	  if( cms2.els_p4()[idx].pt() > 10. && cms2.els_p4()[idx].pt() <= 20. )
	  {
		  if( abs(cms2.els_p4()[idx].eta()) <= 0.8 )
		  {
			  return tightlowptValues[0];
		  }
		  else if( abs(cms2.els_p4()[idx].eta()) > 0.8 && abs(cms2.els_p4()[idx].eta()) <= 1.479 )
		  {
			  return tightlowptValues[1];
		  }
		  else if( abs(cms2.els_p4()[idx].eta()) > 1.479 )
		  {
			  return tightlowptValues[2];
		  }
	  } 
	  else if( cms2.els_p4()[idx].pt() > 20. )
	  {
		  if( abs(cms2.els_p4()[idx].eta()) <= 0.8 )
		  {
			  return tightValues[0];
		  }
		  else if( abs(cms2.els_p4()[idx].eta()) > 0.8 && abs(cms2.els_p4()[idx].eta()) <= 1.479 )
		  {
			  return tightValues[1];
		  }
		  else if( abs(cms2.els_p4()[idx].eta()) > 1.479 )
		  {
			  return tightValues[2];
		  }
	  }  
  }
  else if( lep_type == LeptonType::LOOSEDILEPMVA || lep_type == LeptonType::LOOSETRILEPMVA )
  {
	  if( cms2.els_p4()[idx].pt() > 10. && cms2.els_p4()[idx].pt() <= 20. )
	  {
		  if( abs(cms2.els_p4()[idx].eta()) <= 0.8 )
		  {
			  return looselowptValues[0];
		  }
		  else if( abs(cms2.els_p4()[idx].eta()) > 0.8 && abs(cms2.els_p4()[idx].eta()) <= 1.479 )
		  {
			  return looselowptValues[1];
		  }
		  else if( abs(cms2.els_p4()[idx].eta()) > 1.479 )
		  {
			  return looselowptValues[2];
		  }
	  } 
	  else if( cms2.els_p4()[idx].pt() > 20. )
	  {
		  if( abs(cms2.els_p4()[idx].eta()) <= 0.8 )
		  {
			  return looseValues[0];
		  }
		  else if( abs(cms2.els_p4()[idx].eta()) > 0.8 && abs(cms2.els_p4()[idx].eta()) <= 1.479 )
		  {
			  return looseValues[1];
		  }
		  else if( abs(cms2.els_p4()[idx].eta()) > 1.479 )
		  {
			  return looseValues[2];
		  }
	  } 
  }
  else
  {
	  std::cout<<"Warning: In ttvSelections::getTrigMVAThreshold, unsupported enum "<<lep_type<<" was used as argument."<<endl;
  }
  return 2.; // or should this be a garbage value? this will guarentee that a comparison to a legitimate value will fail.
}
float ttv::getNonTrigMVAThreshold(int idx, LeptonType::value_type lep_type)
{
  // this one is not implemented yet and likely won't be used
  // dummy assignment to suppress warning
  idx = -9999;  
  lep_type = LeptonType::static_size;
  return 2.; // or should this be a garbage value? this will guarentee that a comparison to a legitimate value will fail.
}





std::vector<LorentzVector> ttv::getJets(std::vector<LorentzVector>& leps, enum JetType type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{

  std::vector<LorentzVector> tmp_jets = getJets(999999, true, type, JETS_CLEAN_NONE, deltaR, min_pt, max_eta, (double) rescale, systFlag);

  // ok, now perform the rest of the lepton overlap removal
  // and the impose the pt requirement after applying the
  // extra corrections
  std::vector<LorentzVector> final_jets;

  for (unsigned int jidx = 0; jidx < tmp_jets.size(); jidx++)
  {
	  LorentzVector vjet = tmp_jets.at(jidx);
	  bool jetIsLep = false;
	  for (unsigned int lidx = 0; lidx < leps.size(); lidx++)
      {
		  if (ROOT::Math::VectorUtil::DeltaR(vjet, leps.at(lidx)) < deltaR)
          {
			  jetIsLep = true;
			  break;
          }
      }

	  if (jetIsLep) { continue; }        

	  final_jets.push_back(vjet);
  }

  sort(final_jets.begin(), final_jets.end(), jet_pt_gt_ttv());
  return final_jets;    
}

std::vector<LorentzVector> ttv::getJets(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
  std::vector<bool> tmp_jet_flags = ttv::getJetFlags(leps, type, deltaR, 0., max_eta, 1., 0);

  // now impose the pt requirement after applying the extra corrections
  std::vector<LorentzVector> final_jets;
  for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) 
  {

	if (!tmp_jet_flags.at(jidx)) { continue; }

	LorentzVector vjet = cms2.pfjets_p4().at(jidx);
	jet_corrector->setRho(cms2.evt_ww_rho_vor());
	jet_corrector->setJetA(cms2.pfjets_area().at(jidx));
	jet_corrector->setJetPt(cms2.pfjets_p4().at(jidx).pt());
	jet_corrector->setJetEta(cms2.pfjets_p4().at(jidx).eta());        
	double jet_cor = jet_corrector->getCorrection();

	double unc = (systFlag == 0) ? 1. : getJetMetSyst(systFlag, vjet.pt(), vjet.eta());        
	vjet *= jet_cor * rescale * unc;
	if (vjet.pt() < min_pt)      { continue; }

	final_jets.push_back(vjet);
  }

  sort(final_jets.begin(), final_jets.end(), jet_pt_gt_ttv());
  return final_jets;
}

std::vector<bool> ttv::getJetFlags(std::vector<LorentzVector>& leps, enum JetType type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
  std::vector<bool> tmp_jet_flags = getJetFlags(999999, type, JETS_CLEAN_NONE, (double)deltaR, min_pt, (double)max_eta, (double)rescale, systFlag);
  std::vector<LorentzVector> tmp_jets = getJets(999999, false, type, JETS_CLEAN_NONE, deltaR, min_pt, max_eta, (double) rescale, systFlag);

  // ok, now perform the rest of the lepton overlap removal
  // and the impose the pt requirement after applying the
  // extra corrections
  std::vector<bool> final_jets;
  for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) 
  {

	if (!tmp_jet_flags.at(jidx)) 
	{
	  final_jets.push_back(tmp_jet_flags.at(jidx));
	  continue;
	}
        
	LorentzVector vjet = tmp_jets.at(jidx);
	bool jetIsLep = false;
	for (unsigned int lidx = 0; lidx < leps.size(); lidx++)
	{
		if (ROOT::Math::VectorUtil::DeltaR(vjet, leps.at(lidx)) < deltaR)
		{
			jetIsLep = true;
			break;
		}
	}
            
	final_jets.push_back(!jetIsLep);
  }

  return final_jets;
}

std::vector<bool> ttv::getJetFlags(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
  std::vector<bool> tmp_jet_flags = ttv::getJetFlags(leps, type, deltaR, 0., max_eta, 1., 0);
  std::vector<LorentzVector> const * tmp_jets = NULL;
  if (type <= JETS_TYPE_PF_UNCORR)
	{ tmp_jets = &cms2.pfjets_p4(); }
  else if (type <= JETS_TYPE_CALO_UNCORR)
	{ tmp_jets = &cms2.jets_p4(); }
  else
	{ tmp_jets = &cms2.genjets_p4(); }

  assert(tmp_jets->size() == tmp_jet_flags.size());

  // ok, now perform the rest of the lepton overlap removal
  // and the impose the pt requirement after applying the
  // extra corrections
  std::vector<bool> final_jets;
  for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) 
  {

	if (!tmp_jet_flags.at(jidx)) 
	{
	  final_jets.push_back(tmp_jet_flags.at(jidx));
	  continue;
	}
        
	LorentzVector vjet = tmp_jets->at(jidx);
	jet_corrector->setRho(cms2.evt_ww_rho_vor());
	jet_corrector->setJetA(cms2.pfjets_area().at(jidx));
	jet_corrector->setJetPt(cms2.pfjets_p4().at(jidx).pt());
	jet_corrector->setJetEta(cms2.pfjets_p4().at(jidx).eta());        
	double jet_cor = jet_corrector->getCorrection();

	double unc = (systFlag == 0) ? 1. : getJetMetSyst(systFlag, vjet.pt(), vjet.eta());        
	vjet *= jet_cor * rescale * unc;
	bool is_good_jet = (vjet.pt() >= min_pt);

	final_jets.push_back(is_good_jet);
  }

  return final_jets;
}
    
float ttv::sumJetPt(std::vector<LorentzVector>& leps, enum JetType type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
  std::vector<LorentzVector> tmp_jets = ttv::getJets(leps, type, deltaR, min_pt, max_eta, rescale, systFlag);
  float ht = 0.;
  for (unsigned int jidx = 0; jidx < tmp_jets.size(); jidx++)
	{ ht += tmp_jets.at(jidx).pt(); }

  return ht;
}


float ttv::sumJetPt(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
  std::vector<LorentzVector> tmp_jets = ttv::getJets(leps, jet_corrector, type, deltaR, min_pt, max_eta, rescale, systFlag);
  float ht = 0.;
  for (unsigned int jidx = 0; jidx < tmp_jets.size(); jidx++)
	{ ht += tmp_jets.at(jidx).pt(); }

  return ht;
}

int ttv::nJets(std::vector<LorentzVector>& leps, enum JetType type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
  std::vector<LorentzVector> tmp_jets = ttv::getJets(leps, type, deltaR, min_pt, max_eta, rescale, systFlag);
  return tmp_jets.size();
}

int ttv::nJets(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
  std::vector<LorentzVector> tmp_jets = ttv::getJets(leps, jet_corrector, type, deltaR, min_pt, max_eta, rescale, systFlag);
  return tmp_jets.size();
}

std::vector<LorentzVector> ttv::getBtaggedJets(std::vector<LorentzVector>& leps, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
  std::vector<LorentzVector> tmp_jets = getBtaggedJets(999999, true, type, JETS_CLEAN_NONE, btag_type, deltaR, min_pt, max_eta, (double) rescale, systFlag);

  // ok, now perform the rest of the lepton overlap removal
  // and the impose the pt requirement after applying the
  // extra corrections
  std::vector<LorentzVector> final_jets;
  for (unsigned int jidx = 0; jidx < tmp_jets.size(); jidx++)
  {
	  LorentzVector vjet = tmp_jets.at(jidx);
	  bool jetIsLep = false;
	  for (unsigned int lidx = 0; lidx < leps.size(); lidx++)
      {
		  if (ROOT::Math::VectorUtil::DeltaR(vjet, leps.at(lidx)) < deltaR)
          {
			  jetIsLep = true;
			  break;
          }
      }

	  if (jetIsLep) { continue; }       

	  final_jets.push_back(vjet);
    }

  sort(final_jets.begin(), final_jets.end(), jet_pt_gt_ttv());
  return final_jets;    
}

std::vector<LorentzVector> ttv::getBtaggedJets(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
  std::vector<bool> tmp_jet_flags = ttv::getBtaggedJetFlags(leps, type, btag_type, deltaR, 0., max_eta, 1., 0);

  // now impose the pt requirement after applying the extra corrections
  std::vector<LorentzVector> final_jets;
  for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) 
  {

	if (!tmp_jet_flags.at(jidx))
	  { continue; }

	LorentzVector vjet = cms2.pfjets_p4().at(jidx);
	jet_corrector->setRho(cms2.evt_ww_rho_vor());
	jet_corrector->setJetA(cms2.pfjets_area().at(jidx));
	jet_corrector->setJetPt(cms2.pfjets_p4().at(jidx).pt());
	jet_corrector->setJetEta(cms2.pfjets_p4().at(jidx).eta());        
	double jet_cor = jet_corrector->getCorrection();

	double unc = (systFlag == 0) ? 1. : getJetMetSyst(systFlag, vjet.pt(), vjet.eta());        
	vjet *= jet_cor * rescale * unc;
	if (vjet.pt() < min_pt) { continue; }

	final_jets.push_back(vjet);
  }

  sort(final_jets.begin(), final_jets.end(), jet_pt_gt_ttv());
  return final_jets;
}

std::vector<bool> ttv::getBtaggedJetFlags(std::vector<LorentzVector>& leps, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
  std::vector<bool> tmp_jet_flags = getBtaggedJetFlags(999999, type, JETS_CLEAN_NONE, btag_type, (double)deltaR, (double)min_pt, (double)max_eta, (double)rescale, systFlag);
  std::vector<LorentzVector> tmp_jets = getBtaggedJets(999999, false, type, JETS_CLEAN_NONE, btag_type, deltaR, min_pt, max_eta, (double) rescale, systFlag);

  // ok, now perform the rest of the lepton overlap removal
  // and the impose the pt requirement after applying the
  // extra corrections
  std::vector<bool> final_jets;
  for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

	if (!tmp_jet_flags.at(jidx)) {
	  final_jets.push_back(tmp_jet_flags.at(jidx));
	  continue;
	}
        
	LorentzVector vjet = tmp_jets.at(jidx);
	bool jetIsLep = false;
	for (unsigned int lidx = 0; lidx < leps.size(); lidx++)
	  {
		if (ROOT::Math::VectorUtil::DeltaR(vjet, leps.at(lidx)) < deltaR)
		  {
			jetIsLep = true;
			break;
		  }
	  }
            
	final_jets.push_back(!jetIsLep);
  }

  return final_jets;
}

std::vector<bool> ttv::getBtaggedJetFlags(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
  std::vector<bool> tmp_jet_flags = ttv::getBtaggedJetFlags(leps, type, btag_type, deltaR, 0., max_eta, 1., 0);
  std::vector<LorentzVector> const *tmp_jets = NULL;
  if (type <= JETS_TYPE_PF_UNCORR)
	{ tmp_jets = &cms2.pfjets_p4(); }
  else if (type <= JETS_TYPE_CALO_UNCORR)
	{ tmp_jets = &cms2.jets_p4(); }
  else
	{ tmp_jets = &cms2.genjets_p4(); }

  assert(tmp_jets->size() == tmp_jet_flags.size());

  // ok, now perform the rest of the lepton overlap removal
  // and the impose the pt requirement after applying the
  // extra corrections
  std::vector<bool> final_jets;
  for (unsigned int jidx = 0; jidx < tmp_jet_flags.size(); jidx++) {

	if (!tmp_jet_flags.at(jidx)) {
	  final_jets.push_back(tmp_jet_flags.at(jidx));
	  continue;
	}
        
	LorentzVector vjet = tmp_jets->at(jidx);
	jet_corrector->setRho(cms2.evt_ww_rho_vor());
	jet_corrector->setJetA(cms2.pfjets_area().at(jidx));
	jet_corrector->setJetPt(cms2.pfjets_p4().at(jidx).pt());
	jet_corrector->setJetEta(cms2.pfjets_p4().at(jidx).eta());        
	double jet_cor = jet_corrector->getCorrection();

	double unc = (systFlag == 0) ? 1. : getJetMetSyst(systFlag, vjet.pt(), vjet.eta());        
	vjet *= jet_cor * rescale * unc;
	bool is_good_jet = (vjet.pt() >= min_pt);

	final_jets.push_back(is_good_jet);
  }

  return final_jets;
}

int ttv::nBtaggedJets(std::vector<LorentzVector>& leps, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
  std::vector<LorentzVector> tmp_jets = ttv::getBtaggedJets(leps, type, btag_type, deltaR, min_pt, max_eta, rescale, systFlag);
  return tmp_jets.size();
}

int ttv::nBtaggedJets(std::vector<LorentzVector>& leps, FactorizedJetCorrector* jet_corrector, enum JetType type, enum BtagType btag_type, float deltaR, float min_pt, float max_eta, float rescale, int systFlag)
{
  std::vector<LorentzVector> tmp_jets = ttv::getBtaggedJets(leps, jet_corrector, type, btag_type, deltaR, min_pt, max_eta, rescale, systFlag);
  return tmp_jets.size();
}

LorentzVector ttv::LeptonInfo::p4()
{
  return (abs(id) == 11) ? cms2.els_p4().at(idx) : cms2.mus_p4().at(idx);
}

ttv::DileptonInfo::DileptonInfo (int id1_, int idx1_, int id2_, int idx2_)
{
  lep1.id  = id1_;
  lep1.idx = idx1_;
  lep2.id  = id2_;
  lep2.idx = idx2_;
  hyp_index = -1;
}

ttv::DileptonInfo::DileptonInfo (int hyp_index_)
{
  hyp_index = hyp_index_;
  lep1.id   = cms2.hyp_lt_id().at(hyp_index);
  lep1.idx  = cms2.hyp_lt_index().at(hyp_index);
  lep2.id   = cms2.hyp_ll_id().at(hyp_index);
  lep2.idx  = cms2.hyp_ll_index().at(hyp_index);
}

float ttv::DileptonInfo::sumpt ()
{
  return lep1.p4().pt() + lep2.p4().pt();
}

LorentzVector ttv::DileptonInfo::p4 ()
{
  return (lep1.p4() + lep2.p4());
}

float ttv::DileptonInfo::mass ()
{
  return p4().mass();
}

ttv::TrileptonInfo::TrileptonInfo (int id1_, int idx1_, int id2_, int idx2_, int id3_, int idx3_)
{
  z.lep1.id  = id1_;
  z.lep1.idx = idx1_;
  z.lep2.id  = id2_;
  z.lep2.idx = idx2_;
  w.id       = id3_;
  w.idx      = idx3_;
}

ttv::TrileptonInfo::TrileptonInfo (int hyp_index_, int id3_, int idx3_) :
  z(hyp_index_),
  w(id3_, idx3_)
{}

float ttv::TrileptonInfo::sumpt ()
{
  return z.sumpt() + w.p4().pt();
}
