// C++
#include <iostream>
#include <vector>
#include <set>

// ROOT
#include "TDirectory.h"
#include "TTreeCache.h"
#include "Math/VectorUtil.h"

// CMS2
#include "CORE/CMS2.h"
#include "CORE/mcSelections.h"
#include "CORE/susySelections.h"
#include "CORE/ssSelections.h"
#include "CORE/jetSelections.h"
#include "CORE/hemJet.h"

// header
#include "ScanChain.h"

using namespace std;
using namespace tas;


void babyMaker::ScanChain(TChain* chain, std::string baby_name, bool isMC) {

    if(isMC)
        MakeBabyNtuple( Form("%s_MC.root", baby_name.c_str()) );
    else
        MakeBabyNtuple( Form("%s_Data.root", baby_name.c_str()) );

  // File Loop
  int nDuplicates = 0;
  int nEvents = chain->GetEntries();
  unsigned int nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  while ( (currentFile = (TFile*)fileIter.Next()) ) {

    // Get File Content
    TFile f( currentFile->GetTitle() );
    TTree *tree = (TTree*)f.Get("Events");
    TTreeCache::SetLearnEntries(10);
    tree->SetCacheSize(128*1024*1024);
    cms2.Init(tree);
    
    // Event Loop
    unsigned int nEventsTree = tree->GetEntriesFast();
    for( unsigned int iEvent = 0; iEvent < nEventsTree; ++iEvent) {
    

        //if(nEventsTotal > 1000) break;
      // Get Event Content
      tree->LoadTree(iEvent);
      cms2.GetEntry(iEvent);
      ++nEventsTotal;
    
      // Progress
      CMS2::progress( nEventsTotal, nEventsChain );

      InitBabyNtuple();

      //analysis
     
     // loose cuts to make baby smaller 
      if(evt_pfmet_type1cor() < 20.0) continue;
      if(els_p4().size() + mus_p4().size() < 3) continue;

      pfmet = evt_pfmet();
      pfmet_type1cor = evt_pfmet_type1cor();
      scale1fb = evt_scale1fb();
      isRealData = evt_isRealData();
      event = evt_event();
      lumiBlock = evt_lumiBlock();
      run = evt_run();

      els_p4_b = els_p4();
      mus_p4_b = mus_p4();
      pfjets_p4_b = pfjets_p4();

     
      if(isMC) { 
          genmet = gen_met();
          genps_id_b = genps_id();
          genps_p4_b = genps_p4();
          genps_id_mother_b = genps_id_mother();
          genjets_p4_b = genjets_p4();  
      } else {
          dielectronTrigger = passHLTTrigger(TString("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v19"));
          dimuonTrigger = passHLTTrigger(TString("HLT_Mu17_Mu8_v22"));
          electronmuonTrigger = passHLTTrigger(TString("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9")) ||
                                passHLTTrigger(TString("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9"));
      }

      for(int iEl = 0; iEl < els_p4().size(); iEl++) {
        looseEl.push_back( passElectronSelection_Stop2012_v3(iEl,true,true,false) );
        tightEl.push_back( samesign::isNumeratorLepton(11,iEl) );
      }

      for(int iMu = 0; iMu < mus_p4().size(); iMu++) {
        looseMu.push_back( muonId(iMu, ZMet2012_v1) );
        tightMu.push_back( samesign::isNumeratorLepton(13,iMu) );
      }

      for(int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
          passesLoosePFJetID_b.push_back( passesLoosePFJetID(iJet) );
          corL1FastL2L3.push_back( pfjets_corL1FastL2L3().at(iJet) );
      }
      // evt_scale1fb()
      // els_p4()
      // mus_p4()
      //
      //bool looseEl = passElectronSelection_Stop2012_v3(iEl,true,true,false);
      //bool looseMu = muonId(iMu, ZMet2012_v1);
      //bool tightEl = samesign::isNumeratorLepton(11,goodToP4MapEl[pair[4]]);
      //bool tightMu = samesign::isNumeratorLepton(13,goodToP4MapMu[pair[4]]);
      //
      // pfjets_p4()
      // pfjets_corL1FastL2L3(iJet)
      // passesLoosePFJetID(iJet)
      // evt_pfmet_type1cor() (cut on MET at 20? here)
      // 
      // genps_id()
      // genps_p4()
      // genps_id_mother()
      // genjets_p4()
      // gen_met()
      // 
      // evt_isRealData()
      // evt_event()
      // evt_lumiBlock()
      // evt_run()
      //

      FillBabyNtuple();

   }//end loop on events in a file
  
    delete tree;
    f.Close();
  }//end loop on files
  
  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }

  cout << nDuplicates << " duplicate events were skipped." << endl;

  CloseBabyNtuple();

  return;
}
