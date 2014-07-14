## Instructions

# Chunk submissions

* Modify ScanChain.C and ScanChain.h in the `resources` folder.
* Then do something like this, for example:
```
. createChunks.sh mc "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/*.root" baby_DYJetsToLL_V1 10
```
(Note that the files were placed in quotes to prevent bash from expanding the `*`.)

This tells the baby maker that it is data (so doesn't fill in gen information, depending on if you kept that feature in ScanChain.C), and tells the TChain that the files are all `*.root` files inside `/hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-24/DoubleMu_Run2012C-PromptReco-v2_AOD/merged/`. A folder named `baby_DoubleMu_Run2012C_VN1` will be created and a condor job with 10 processes submitted. Finally, the output baby `.root` file will have the name `baby_DoubleMu_Run2012C_VN1_Data.root`.
(The python helper script splits up the `.root` files into approximately equal chunks where the number of chunks is specified by user.)
