if [ $# -ne 4 ]; then
    echo "Usage: . createChunks.sh [\"MC\"/\"data\"] [input file names] [baby/folder name] [num chunks]"
    echo "Note: quotes around the files prevents bash from expanding wildcards!"
    return
fi

# die if no voms proxy
vomsInfo=$(voms-proxy-info | grep "path")
if [ -z "$vomsInfo" ]; then
    echo ">>> Output of command \"voms-proxy-info\" is malformed"
    echo ">>> Create a voms proxy before proceeding"
    echo ">>>   \" voms-proxy-init -valid 120:00 \""
    return
fi

# substitute voms proxy into condor submission file
proxyFile=$(echo $vomsInfo | awk '{print $3}')

# prevent * from expanding
set -f
python split.py $1 "$2" $3 $4 $proxyFile
set +f
