import sys, glob, os

def chunks(ls, n):
    return [ ls[i::n] for i in xrange(n) ]

if(len(sys.argv) < 6): 
    print "python script dying :("
    sys.exit()

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
print ">>> Number of files: %i" % len(files)
print ">>> Number of chunks: %i" % numChunks
print ">>> Proxy file: %s" % proxyFile

condorFileStr = """
universe = vanilla
Executable = CONDOR_EXECUTABLE
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Output = out.stdout.$(Process)
Error = out.stderr.$(Process)
Log = out.log.$(Process)
+DESIRED_Sites="UCSD"
+Owner = undefined
notify_user = dummymail000@FNAL.GOV
x509userproxy = PROXY_FILE
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

folderName = babyName
condorName = babyName + ".sh"
os.system("mkdir -p %s" % folderName)

fcondorFile = open("%s/condorFile" % folderName, "w")
fcondorFile.write( condorFileStr.replace("CONDOR_EXECUTABLE", condorName) )

os.system("cp -p %s/ScanChain.* %s/" % (resources, folderName))
os.system("cp -rp %s/CORE %s/" % (resources, folderName))
os.system("cp -p %s/condorExecutable.sh %s/%s" % (resources, folderName, condorName))

for chunk in range(numChunks):

    doAllName = "doAll_%i.C" % chunk
    print ">>> Making", doAllName, "for chunk",chunk,"with",len(chunks(files,numChunks)[chunk]),"files"

    fdoAll = open("%s/%s" % (folderName, doAllName), "w")

    fdoAll.write( doAllStrBegin )
    for fileName in chunks(files, numChunks)[chunk]:
        fdoAll.write( "    ch->Add(\"%s\");\n" % fileName )
    fdoAll.write( doAllStrEnd.replace("BABY_NAME", folderName+"_"+str(chunk)) )

    fcondorFile.write( "\n\n" )
    fcondorFile.write( "Transfer_Input_Files = ScanChain.C,ScanChain.h,%s,CORE\n" % doAllName )
    fcondorFile.write( "Arguments = %i\n" % chunk )
    fcondorFile.write( "Queue 1\n" )

    fdoAll.close()

fcondorFile.close()

print ">>> Submitting", folderName
os.system("(cd %s && condor_submit condorFile)" % folderName)

print ">>> Done"
