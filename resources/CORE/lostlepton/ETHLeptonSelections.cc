#include "ETHLeptonSelections.h"
#include "../electronSelections.h"
#include "../muonSelections.h"
#include "../trackSelections.h"
#include "../eventSelections.h"

bool ETHSelection::goodElectron(int index){

  if (cms2.els_p4()[index].pt() < 10.0) return false;
  if (fabs(cms2.els_p4()[index].eta()) > 2.4) return false;
  if (fabs(cms2.els_p4()[index].eta()) > 1.442 && fabs(cms2.els_p4()[index].eta()) < 1.566) return false;
  if (fabs(electron_d0PV(index)) > 0.04) return false;
  if (fabs(electron_dzPV_smurfV3(index)) > 0.2 ) return false;
  //if (electronIsoValuePF(index, 0) > 0.15) return false;
  if(electronIsoValuePF2012_FastJetEffArea_v3(index) > 0.15) return false;

  if (fabs(cms2.els_p4()[index].eta()) < 1.479) { //barrel
    if (cms2.els_dEtaIn()[index] > 0.007) return false;
    if (cms2.els_dPhiIn()[index] > 0.8) return false;
    if (cms2.els_sigmaIEtaIEta()[index] > 0.01) return false;
    if (cms2.els_hOverE()[index] > 0.15) return false;
  }
  else { //endcap
    if (cms2.els_dEtaIn()[index] > 0.01) return false;
    if (cms2.els_dPhiIn()[index] > 0.7) return false;
    if (cms2.els_sigmaIEtaIEta()[index] > 0.03) return false;
  }

  return true;

}

bool ETHSelection::goodMuon(int index){

  if (cms2.mus_p4()[index].pt() < 10.0) return false;
  if (fabs(cms2.mus_p4()[index].eta()) > 2.4) return false;

  const bool is_global  = ((cms2.mus_type().at(index) & (1<<1)) != 0);
  const bool is_pfmu    = ((cms2.mus_type().at(index) & (1<<5)) != 0);

  const int vtxidx = firstGoodVertex();

  if (!is_global) return false;
  if (!is_pfmu) return false;
  if (cms2.mus_gfit_validSTAHits().at(index) < 1) return false;
  if (cms2.mus_numberOfMatchedStations().at(index) < 2) return false;

  const int ctfidx = cms2.mus_trkidx().at(index);
  if (ctfidx < 0 || vtxidx < 0) return false;
  const std::pair<double, double> cord0 = trks_d0_pv(ctfidx, vtxidx);
  const std::pair<double, double> cordz = trks_dz_pv(ctfidx, vtxidx);
  if (fabs(cord0.first) > 0.2) return false;
  if (fabs(cordz.first) > 0.5) return false;
  if (cms2.trks_valid_pixelhits().at(ctfidx) < 1) return false;
  if (cms2.trks_nlayers().at(ctfidx) < 6) return false;

  if ( cms2.mus_gfit_chi2()[index] / cms2.mus_gfit_ndof()[index] >= 10 ) return false;
  if (muonIsoValuePF2012_deltaBeta(index) > 0.2) return false;

  return true;

}

bool ETHSelection::goodTau(int index){

  if(cms2.taus_pf_p4()[index].pt() < 20) return false;
  if(fabs(cms2.taus_pf_p4()[index].eta()) > 2.3) return false;
  if(cms2.taus_pf_byDecayModeFinding()[index] < 0.5) return false;
  if(cms2.taus_pf_againstElectronLoose()[index] < 0.5) return false;
  if(cms2.taus_pf_againstMuonTight2()[index] < 0.5) return false;
  if(cms2.taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits()[index] < 0.5) return false;
  
  return true;
}
