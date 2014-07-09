#!/bin/bash

# only works on el6 comps
# submitting from uaf-7 works fine
uname -a
echo "initial files"
ls -l

#This stuff to get it to run
export CMS_PATH=/code/osgcode/cmssoft/cms
source /code/osgcode/cmssoft/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc472
source /code/osgcode/fgolf/5.30-patches/bin/thisroot.sh
export ROOTSYS=$ROOTSYS:$HOME/CORE
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH:$HOME/CORE
export PATH=$HOME/bin:$ROOTSYS/bin:$HOME/CORE:/opt/d-cache/dcap/bin:$PATH
export PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH


echo $SCRAM_ARCH

scramv1 project CMSSW CMSSW_5_3_19
cd CMSSW_5_3_19/src
eval `scramv1 runtime -sh`

echo $CMSSW_RELEASE_BASE

# copy files into CMSSW_.../src directory
cp -rp ../../CORE .
cp -p ../../ScanChain.{C,h} .
cp -p ../../doAll.C .
root -b -q -n -l doAll.C
echo "files after running root"
ls -l

# move root file back and delete CMSSW folder -- performance improvement
mv *.root ../../
cd ../../
rm -rf CMSSW*/

echo "final files"
ls -l

