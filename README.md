## Instructions

# Single job submissions

* Check out repository (obviously)
* Modify ScanChain.C and ScanChain.h in the `resources` folder.
* Then do:
```. create.sh [MC/data] [input folder] [baby name]```

For example, consider 
```. create.sh Data "/hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-24/DoubleMu_Run2012C-PromptReco-v2_AOD/merged/*.root" baby_DoubleMu_Run2012C_VN1```
(Note that the files were placed in quotes to prevent bash from expanding the `*`.)

This tells the baby maker that it is data (so doesn't fill in gen information, depending on if you kept that feature in ScanChain.C), and tells the TChain that the files are all `*.root` files inside `/hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-24/DoubleMu_Run2012C-PromptReco-v2_AOD/merged/`. A folder named `baby_DoubleMu_Run2012C_VN1` will be created and a condor job submitted. Finally, the output baby `.root` file will have the name `baby_DoubleMu_Run2012C_VN1_Data.root`.


# Chunk submissions

* Check out repository (obviously)
* Modify ScanChain.C and ScanChain.h in the `resources` folder.
* Then do:
```. createChunks.sh ["MC"/"data"] [input file names] [baby/folder name] [num chunks]"``

For example, consider
```
. createChunks.sh mc "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/*.root" baby_DYJetsToLL_V1 10
```

This creates 10 jobs that will run over the given DY sample, splitting up the root files into 10 approximately equal chunks to run in parallel. Thus, 10 folders will be created with the names baby_DYJetsToLL_V1_000 to baby_DYJetsToLL_V1_009, each containing the corresponding `.root` file for that job. Afterwards, they can be combined (if needed) using the `hadd` command.
