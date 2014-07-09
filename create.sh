if [ $# -ne 3 ]; then
    echo "Usage: . create.sh [MC/data] [folder name] [baby name]"
    return
fi

isMC=$1
inputFolder=$2
babyName=$3
resources="resources"

mkdir $babyName
cd $babyName

echo "MC or Data: $isMC"
echo "Baby name (and folder name): $babyName"
echo "Input folder: $inputFolder"

cp -p ../$resources/doAll.C .
cp -p ../$resources/ScanChain.* .
cp -p ../$resources/condor* .
cp -rp ../$resources/CORE .

# give ROOT script info on data/MC
if [ "$isMC" == "MC" ] || [ "$isMC" == "mc" ]; then
    sed -i "s/IS_MC/true/g" doAll.C
else
    sed -i "s/IS_MC/false/g" doAll.C
fi

# substitute baby name
sed -i "s/BABY_NAME/$babyName/g" doAll.C

# substitute input file names
sed -i "s,INPUT_NAMES,$inputFolder/*.root,g" doAll.C

# substitute voms proxy into condor submission file
vomsInfo=$(voms-proxy-info | grep "path")
if [ -z "$vomsInfo" ]; then
    echo "Output of command \"voms-proxy-info\" is malformed"
    echo "Create a voms proxy before proceeding"
    return
fi
proxyFile=$(echo $vomsInfo | awk '{print $3}')
cp -p ../$resources/condorFile .
sed -i "s,PROXY_FILE,$proxyFile," condorFile

echo "Submitting $babyName"
condor_submit condorFile


# -rwxr-xr-x  1 namin namin  986 Jul  8 23:04 condorExecutable.sh
# -rw-rw-r--  1 namin namin  344 Jul  8 23:08 condorFile
# drwxr-xr-x 11 namin namin 4096 Jun 30 23:20 CORE
# -rw-r--r--  1 namin namin  750 Jul  8 22:51 doAll.C
# -rw-r--r--  1 namin namin 3990 Jul  8 18:26 ScanChain.C
# -rw-r--r--  1 namin namin 3332 Jul  8 16:24 ScanChain.h

