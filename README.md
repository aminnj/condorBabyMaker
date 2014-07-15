## Instructions

# Chunk submissions

* Modify `ScanChain.C` and `ScanChain.h` in the `resources` folder. Keep those filenames the same. `CORE` can also be swapped out for a newer version.
* Source `createChunks.sh` with no parameters for argument hints.

```
$ . createChunks.sh
Usage: . createChunks.sh ["MC"/"data"] [input file names] [baby/folder name] [num chunks]
Note: quotes around the files prevents bash from expanding wildcards!

```
* Then do something like this, for example:

```
. createChunks.sh mc "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/*.root" baby_DYJetsToLL_V1 10
```
(Note that the files were placed in quotes to prevent bash from expanding the `*`.)

This tells the baby maker that it is data (so doesn't fill in gen information,
depending on if you kept that feature in ScanChain.C), and tells the TChain
that the files are all `*.root` files inside
`/hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-24/DoubleMu_Run2012C-PromptReco-v2_AOD/merged/`.
A folder named `baby_DoubleMu_Run2012C_VN1` will be created and a condor job
with 10 processes submitted. Finally, the output baby `.root` files will have
the names `baby_DoubleMu_Run2012C_VN1_#_Data.root` where `#` corresponds to the
process number. Afterwards, you can use the `hadd` command to combine them into
one `.root` file. In thise case, you could do:

```
hadd baby_DoubleMu_Run2012C_VN1.root baby_DoubleMu_Run2012C_VN1_*.root
```

(The python helper script splits up the `.root` files into approximately equal
chunks where the number of chunks is specified by user.)
