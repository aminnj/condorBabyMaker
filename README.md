## Instructions

* Check out repository (obviously)
* Modify ScanChain.C and ScanChain.h in the `resources` folder.
* Then do:
```. create.sh [MC/data] [input folder] [baby name]```

For example, consider 
```. create.sh Data "/hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-24/DoubleMu_Run2012C-PromptReco-v2_AOD/merged/*.root" baby_DoubleMu_Run2012C_VN1```
(Note that the files were placed in quotes to prevent bash from expanding the `*`.)

This tells the baby maker that it is data (so doesn't fill in gen information, depending on if you kept that feature in ScanChain.C), and tells the TChain that the files are all `*.root` files inside `/hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-24/DoubleMu_Run2012C-PromptReco-v2_AOD/merged/`. A folder named `baby_DoubleMu_Run2012C_VN1` will be created and a condor job submitted. Finally, the output baby `.root` file will have the name `baby_DoubleMu_Run2012C_VN1_Data.root`.
