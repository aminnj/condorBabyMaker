#include "MITConversionUtilities.h"
#include "CMS2.h"
#include "TMath.h"


bool isMITConversion(unsigned int elidx, 
		     int nWrongHitsMax, 
		     float probMin,
		     float dlMin,
		     bool matchCTF,
		     bool requireArbitratedMerged) {

  unsigned int nconvs = cms2.convs_isConverted().size();
  if(nconvs == 0) 
    return false;
  bool isGoodConversion = false;

  for(unsigned int iconv = 0; iconv < nconvs; iconv++) {
    
    bool conversionMatchFound = false;
    for(unsigned int itk = 0; itk < cms2.convs_tkidx()[iconv].size(); itk++) {

      if(cms2.convs_tkalgo()[iconv][itk] == 29 && cms2.convs_tkidx()[iconv][itk] == cms2.els_gsftrkidx()[elidx])
	conversionMatchFound = true;
      if(matchCTF) {
	if(cms2.convs_tkalgo()[iconv][itk] > 3 && cms2.convs_tkalgo()[iconv][itk] < 14 && cms2.convs_tkalgo()[iconv][itk] != 12 && cms2.convs_tkidx()[iconv][itk] == cms2.els_trkidx()[elidx])
	  conversionMatchFound = true;
      }
    
      if(conversionMatchFound)
	break;
    }
    
    
    if(conversionMatchFound==false)
      continue;
    
    if( TMath::Prob( cms2.convs_chi2()[iconv], (Int_t)cms2.convs_ndof()[iconv] )  > probMin && cms2.convs_dl()[iconv] > dlMin ) isGoodConversion = true;
    if(requireArbitratedMerged) {
      if(cms2.convs_quality()[iconv] & 4)
	isGoodConversion = true;
      else 
	isGoodConversion = false;
    }

    for(unsigned int j = 0; j < cms2.convs_nHitsBeforeVtx()[iconv].size(); j++) {
      if(cms2.convs_nHitsBeforeVtx()[iconv][j] > nWrongHitsMax)
	isGoodConversion = false;
    }
      
    if(isGoodConversion)
      break;
      
      
  }//loop over convserions


  return isGoodConversion;
}
