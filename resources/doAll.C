{
    gSystem->AddIncludePath("-I./CORE");
    gSystem->Load("./CORE/libCMS2NtupleMacrosCORE.so");  

    // following 4 lines only work if CMSSW_BASE is set (so cmsenv required)
    gSystem->Load("libFWCoreFWLite.so");
    AutoLibraryLoader::enable();
    gSystem->Load("libDataFormatsFWLite.so");
    gSystem->Load("libDataFormatsPatCandidates.so");

    gROOT->ProcessLine(".L ScanChain.C++");
    //gSystem->Load("ScanChain_C.so");

    TChain *ch = new TChain("Events"); 
    ch->Add("INPUT_NAMES");
    babyMaker *looper = new babyMaker();
    //                            isMC
    looper->ScanChain(ch, "BABY_NAME", IS_MC); 

}
