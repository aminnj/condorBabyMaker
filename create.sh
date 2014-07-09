if [ $# -ne 3 ]; then
    echo "Usage: . create.sh [MC/data] [folder name] [baby name]"
    return
fi

isMC=$1
inputFolder=$2
babyName=$3
resources="resources"


echo ">>> MC or Data: $isMC"
echo ">>> Baby name (and folder name): $babyName"
echo ">>> Input folder: $inputFolder"

# die if no voms proxy
vomsInfo=$(voms-proxy-info | grep "path")
if [ -z "$vomsInfo" ]; then
    echo ">>> Output of command \"voms-proxy-info\" is malformed"
    echo ">>> Create a voms proxy before proceeding"
    echo ">>>   \" voms-proxy-init -valid 120:00 \""
    return
fi

mkdir $babyName
cd $babyName


cp -p ../$resources/doAll.C .
cp -p ../$resources/ScanChain.* .
cp -p ../$resources/condor* .
cp -rp ../$resources/CORE .

echo ">>> Copied files from $resources"

# give ROOT script info on data/MC
if [ "$isMC" == "MC" ] || [ "$isMC" == "mc" ]; then
    echo ">>> Telling script this is MC"
    sed -i "s/IS_MC/true/g" doAll.C
else
    echo ">>> Telling script this is not MC"
    sed -i "s/IS_MC/false/g" doAll.C
fi

# substitute baby name
sed -i "s/BABY_NAME/$babyName/g" doAll.C

# substitute input file names
sed -i "s,INPUT_NAMES,$inputFolder/*.root,g" doAll.C

# substitute voms proxy into condor submission file
proxyFile=$(echo $vomsInfo | awk '{print $3}')
echo ">>> Telling condor to use proxy file $proxyFile"
cp -p ../$resources/condorFile .
sed -i "s,PROXY_FILE,$proxyFile," condorFile

echo ">>> Submitting $babyName"
condor_submit condorFile

