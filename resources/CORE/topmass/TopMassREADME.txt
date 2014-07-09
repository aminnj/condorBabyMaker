//------------------------
Base requirement: 
//------------------------
You need a version of LHAPDF installed. On our tas machines you can
find it at: /tas07/disk00/jribnik/lhapdf-5.8.3/

LHAPDF source code is 32b, to compile it at the UAF ( 64b arch ) use:
  > CFLAGS="-m32" CXXFLAGS="-m32" LDFLAGS="-m32" ./configure --prefix=/home/users/dbarge/LHAPD --disable-pyext
  > make
  > make install

Download your favourite PDF and direct the code to use it:
 - download from http://projects.hepforge.org/lhapdf/
 (e.g. 
wget http://svn.hepforge.org/lhapdf/pdfsets/5/cteq61.LHgrid

and point LHAPDF::initPDFSet in ttdilepsolve.cpp to where you saved
the file. The default location is CORE/topmass/pdfs

Also maybe adjust the resolution file path:
string JR_Standalone_Path = "CORE/topmass/JR_Standalone/txts/";
in getTopMassEstimate.icc

//------------------------
The following additions are needed:
//------------------------
-- in your root setup:
// Load the LHA library
gSystem->Load("/tas07/disk00/jribnik/lhapdf-5.8.3/lib/.libs/libLHAPDF.so"); //on tas
gSystem->Load("/nfs-3/userdata/kalavase/lhapdf-5.8.3/lib/.libs/libLHAPDF.so"); //on uaf
// your root Include Path should have at leas this - look e.g. at
setup.C in OSSusy:
gSystem->AddIncludePath(" -w -I../CORE/topmass -I/tas07/disk00/jribnik/lhapdf/include"); //tas
gSystem->AddIncludePath(" -w -I../CORE/topmass -I/nfs-3/userdata/kalavase/lhapdf/include"); //uaf
-- in your macro, before looper->ScanChain(...
	 and after you set the stuff above:


//Load the ttbar mass solver (incl. LHA lib)
 gROOT->ProcessLine(".L ../CORE/topmass/ttdilepsolve.cpp+");

-- in your code:
//add these includes:
#include "../CORE/topmass/ttdilepsolve.h"
#include "../CORE/topmass/getTopMassEstimate.icc"

-- instanciate at top of ScanChain:

ttdilepsolve * d_llsol = new ttdilepsolve;

-- and where you want the top mass use:

float topMassEst = getTopMassEstimate(d_llsol, hypIdx, Yourjets_p4(), tcmet, tcmetphi);

-- to enable jet smearing:
Call the Estimator this way:
float topMassEst = getTopMassEstimate(d_llsol, hypIdx, Yourjets_p4(), tcmet, tcmetphi, 100);
The "100" at the end is the number of smear iterations. Aram recommends
100 for MC or 1000 for data 
The smearing code is horribly slow (seconds per ttbar event) - so run it only on the final set of
a few skimmed events... or compile it properly.
