CXX=g++

ROOTCFLAGS  = $(shell root-config --cflags)
ROOTLIBS    = $(shell root-config --libs)

ROOTIL      = $(ROOTCFLAGS) $(ROOTLIBS)

MyTest:
	$(CXX) $(ROOTIL)  -o MyTestJetRes JetResolution_t.cc JetResolution.cc JetCorrectorParameters_tm.cc

Standalone:
	$(CXX) $(ROOTIL)  -o standalone Standalone.cc JetResolution.cc JetCorrectorParameters_tm.cc


