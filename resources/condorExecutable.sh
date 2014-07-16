#!/bin/bash

# only works on el6 comps
# submitting from uaf-7 or uaf-9 works fine
uname -a
echo "initial files"
ls -l

job=$1
echo $job

#This stuff to get it to run
export CMS_PATH=/code/osgcode/cmssoft/cms
source /code/osgcode/cmssoft/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc472
source /code/osgcode/fgolf/5.30-patches/bin/thisroot.sh
export ROOTSYS=$ROOTSYS:$HOME/CORE
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH:$HOME/CORE
export PATH=$HOME/bin:$ROOTSYS/bin:$HOME/CORE:/opt/d-cache/dcap/bin:$PATH
export PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH

export HADOOPDIR=/hadoop/cms/store/user/namin/wz/
export BABYNAME=$(grep "looper->" doAll_${job}.C | cut -d '"' -f2 | sed "s/_[0-9]*$//")

echo $SCRAM_ARCH

scramv1 project CMSSW CMSSW_5_3_19
cd CMSSW_5_3_19/src
eval `scramv1 runtime -sh`

echo $CMSSW_RELEASE_BASE

# copy files into CMSSW_.../src directory
cp -rp ../../CORE .
cp -p ../../ScanChain.{C,h} .
cp -p ../../doAll_${job}.C .
root -b -q -n -l doAll_${job}.C
echo "files after running root"
ls -l

echo ${HADOOPDIR}/${BABYNAME}
output=$(ls -1 *.root)
echo $output
lcg-cp -b -D srmv2 --vo cms --connect-timeout 2400 --verbose file://`pwd`/${output} srm://bsrm-3.t2.ucsd.edu:8443/srm/v2/server?SFN=${HADOOPDIR}/${BABYNAME}/${job}.root
cd ../../
rm -rf CMSSW*/


