In a simple looper ( made with makeCMS2ClassFiles.C ) it should be enough to:

	#include "MT2.cc"

provided your enivronment is set up so that CMS2/NtupleMacros/CORE/MT2 is in your ROOT path. This can be done with:
 
	gSystem->AddIncludePath()

Then you can access the MT2 functions

	MT2( MET_MAGNITUDE, MET_PHI, P4_LEPTON_1, P4_LEPTON_2 )
	MT2J( MET_MAGNITUDE, MET_PHI, P4_LEPTON_1, P4_LEPTON_2, VECT_P4_Jets )

in your looper.
