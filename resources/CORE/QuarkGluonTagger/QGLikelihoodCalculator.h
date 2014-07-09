// ------------------------------------------------------------
//  
//    QGLikelihoodCalculator - Class
//    for the computation of the QG likelihood.
//
// ------------------------------------------------------------

#ifndef QGLikelihoodCalculator_h
#define QGLikelihoodCalculator_h

#include <string>
#include <vector>
//#include "TFormula.h"
//#include "../CMS2.h"
//#include "Math/VectorUtil.h"

#include <vector>
#include <string>
#include <utility>

#include "../jetcorr/JetCorrectorParameters.h"
#include "../jetcorr/SimpleJetCorrector.h"

class QGLikelihoodCalculator {

 public:

  QGLikelihoodCalculator( const std::string& fileName_nCharged="UserCode/pandolf/QuarkGluonTagger/data/QGTaggerConfig_nCharged_AK5PF.txt", const std::string& fileName_nNeutral="UserCode/pandolf/QuarkGluonTagger/data/QGTaggerConfig_nNeutral_AK5PF.txt", const std::string& fileName_ptD="UserCode/pandolf/QuarkGluonTagger/data/QGTaggerConfig_ptD_AK5PF.txt");

   ~QGLikelihoodCalculator();

   float computeQGLikelihood( float pt, float rhoPF, int nCharged, int nNeutral, float ptD );
  
 private:

  JetCorrectorParameters *jcp_nCharged_quark_;
  JetCorrectorParameters *jcp_nCharged_gluon_;
  JetCorrectorParameters *jcp_nNeutral_quark_;
  JetCorrectorParameters *jcp_nNeutral_gluon_;
  JetCorrectorParameters *jcp_ptD_quark_;
  JetCorrectorParameters *jcp_ptD_gluon_;

  SimpleJetCorrector *sjc_nCharged_quark_;
  SimpleJetCorrector *sjc_nCharged_gluon_;
  SimpleJetCorrector *sjc_nNeutral_quark_;
  SimpleJetCorrector *sjc_nNeutral_gluon_;
  SimpleJetCorrector *sjc_ptD_quark_;
  SimpleJetCorrector *sjc_ptD_gluon_;

};


#endif
