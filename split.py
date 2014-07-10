import sys, glob, os

def chunks(ls, n):
    return [ ls[i::n] for i in xrange(n) ]

if(len(sys.argv) < 6): 
    print "python script dying :("
    sys.exit()

print sys.argv
# XXX pass in isMC wildcard babyName, numChunks,proxyFile
#wildcard = "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-23/*1.root"
#numChunks = 3
#isMC=(sys.argv[1].lower() == "mc")
#babyName="baby_TEST"
#proxyFile="/tmp/x509up_u31567"

isMC=(sys.argv[1].lower() == "mc")
wildcard = sys.argv[2]
babyName = sys.argv[3]
numChunks = int(sys.argv[4])
proxyFile = sys.argv[5]
resources="resources"

files = glob.glob(wildcard)

print ">>> Is this MC?: %s" % isMC
print ">>> Baby name (and folder name): %s" % babyName
print ">>> Input files: %s" % resources
print ">>> Number of chunks: %i" % numChunks

condorFileStr = """
universe = vanilla
Executable = CONDOR_EXECUTABLE
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = ScanChain.C,ScanChain.h,doAll.C,CORE
Output = out.stdout
Error = out.stderr
Log = out.log
+DESIRED_Sites="UCSD"
+Owner = undefined
notify_user = dummymail000@FNAL.GOV
x509userproxy = PROXY_FILE
Queue 1
"""
doAllStrBegin = """
{
    gSystem->AddIncludePath("-I./CORE");
    gSystem->Load("./CORE/libCMS2NtupleMacrosCORE.so");
    gSystem->Load("libFWCoreFWLite.so");
    AutoLibraryLoader::enable();
    gSystem->Load("libDataFormatsFWLite.so");
    gSystem->Load("libDataFormatsPatCandidates.so");
    gROOT->ProcessLine(".L ScanChain.C++");
    TChain *ch = new TChain("Events");
"""

doAllStrEnd = """
    babyMaker *looper = new babyMaker();
    looper->ScanChain(ch, "BABY_NAME", IS_MC);
}
"""
doAllStrEnd = doAllStrEnd.replace("IS_MC", str(isMC).lower())
condorFileStr = condorFileStr.replace("PROXY_FILE", proxyFile)

#for chunk in range(3,4): # FIXME range(numChunks):
for chunk in range(numChunks):
    folderName = babyName + "_" + str(chunk).zfill(3)
    condorName = babyName + "_" + str(chunk).zfill(3) + ".sh"

    print ">>> Making", folderName
    os.system("mkdir -p %s" % folderName)
    os.system("cp -p %s/ScanChain.* %s/" % (resources, folderName))
    os.system("cp -rp %s/CORE %s/" % (resources, folderName))
    os.system("cp -p %s/condorExecutable.sh %s/%s" % (resources, folderName, condorName))

    fdoAll = open("%s/doAll.C" % folderName, "w")
    fcondorFile = open("%s/condorFile" % folderName, "w")

    #print doAllStrBegin
    fdoAll.write( doAllStrBegin )
    for fileName in chunks(files, numChunks)[chunk]:
        #print "    ch->Add(\"%s\");" % fileName
        fdoAll.write( "    ch->Add(\"%s\");\n" % fileName )
    #print doAllStrEnd.replace("BABY_NAME", folderName)
    fdoAll.write( doAllStrEnd.replace("BABY_NAME", folderName) )

    #print condorFileStr.replace("CONDOR_EXECUTABLE", condorName)
    fcondorFile.write( condorFileStr.replace("CONDOR_EXECUTABLE", condorName) )
    print ">>> Made", folderName

    fdoAll.close()
    fcondorFile.close()
    print ">>> Submitting", folderName
    os.system("(cd %s && condor_submit condorFile)" % folderName)


