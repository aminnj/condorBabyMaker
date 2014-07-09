//#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "QGLikelihoodCalculator.h"

using namespace std;

QGLikelihoodCalculator::QGLikelihoodCalculator( const std::string& fileName_nCharged, const std::string& fileName_nNeutral, const std::string& fileName_ptD ) {

  /*
  // those input files are hard coded for now
  std::string path_nCharged = "CORE/QuarkGluonTagger/QGTaggerConfig_nCharged_AK5PF.txt";
  std::string path_nNeutral = "CORE/QuarkGluonTagger/QGTaggerConfig_nNeutral_AK5PF.txt";
  std::string path_ptD = "CORE/QuarkGluonTagger/QGTaggerConfig_ptD_AK5PF.txt";
  */

  jcp_nCharged_quark_ = new JetCorrectorParameters(fileName_nCharged, "quark");
  jcp_nCharged_gluon_ = new JetCorrectorParameters(fileName_nCharged, "gluon");

  jcp_nNeutral_quark_ = new JetCorrectorParameters(fileName_nNeutral, "quark");
  jcp_nNeutral_gluon_ = new JetCorrectorParameters(fileName_nNeutral, "gluon");

  jcp_ptD_quark_ = new JetCorrectorParameters(fileName_ptD, "quark");
  jcp_ptD_gluon_ = new JetCorrectorParameters(fileName_ptD, "gluon");

  /*
  //
  // COMMENT for NOW
  //

  //check that provided files are for correct variables:
  if( jcp_nCharged_quark_->definitions().level() != "QGL_nCharged_quark" )
    throw cms::Exception("QuarkGluonTagger Config File Error") << "quark section of file \'" << fileName_nCharged << "\' is not of the proper format. Check your input files.";
  if( jcp_nCharged_gluon_->definitions().level() != "QGL_nCharged_gluon" )
    throw cms::Exception("QuarkGluonTagger Config File Error") << "gluon section of file \'" << fileName_nCharged << "\' is not of the proper format. Check your input files.";
  if( jcp_nNeutral_quark_->definitions().level() != "QGL_nNeutral_quark" )
    throw cms::Exception("QuarkGluonTagger Config File Error") << "quark section of file \'" << fileName_nNeutral << "\' is not of the proper format. Check your input files.";
  if( jcp_nNeutral_gluon_->definitions().level() != "QGL_nNeutral_gluon" )
    throw cms::Exception("QuarkGluonTagger Config File Error") << "gluon section of file \'" << fileName_nNeutral << "\' is not of the proper format. Check your input files.";
  if( jcp_ptD_quark_->definitions().level() != "QGL_ptD_quark" )
    throw cms::Exception("QuarkGluonTagger Config File Error") << "quark section of file \'" << fileName_ptD << "\' is not of the proper format. Check your input files.";
  if( jcp_ptD_gluon_->definitions().level() != "QGL_ptD_gluon" )
    throw cms::Exception("QuarkGluonTagger Config File Error") << "gluon section of file \'" << fileName_ptD << "\' is not of the proper format. Check your input files.";
  */

  sjc_nCharged_quark_ = new SimpleJetCorrector(*jcp_nCharged_quark_);
  sjc_nCharged_gluon_ = new SimpleJetCorrector(*jcp_nCharged_gluon_);

  sjc_nNeutral_quark_ = new SimpleJetCorrector(*jcp_nNeutral_quark_);
  sjc_nNeutral_gluon_ = new SimpleJetCorrector(*jcp_nNeutral_gluon_);

  sjc_ptD_quark_ = new SimpleJetCorrector(*jcp_ptD_quark_);
  sjc_ptD_gluon_ = new SimpleJetCorrector(*jcp_ptD_gluon_);

}


QGLikelihoodCalculator::~QGLikelihoodCalculator() {

  delete jcp_nCharged_quark_;
  delete jcp_nCharged_gluon_;
  delete jcp_nNeutral_quark_;
  delete jcp_nNeutral_gluon_;
  delete jcp_ptD_quark_;
  delete jcp_ptD_gluon_;

  delete sjc_nCharged_quark_;
  delete sjc_nCharged_gluon_;
  delete sjc_nNeutral_quark_;
  delete sjc_nNeutral_gluon_;
  delete sjc_ptD_quark_;
  delete sjc_ptD_gluon_;

}




float QGLikelihoodCalculator::computeQGLikelihood( float pt, float rhoPF, int nCharged, int nNeutral, float ptD ) {

  //  cout << "pt " << pt << " rhoPF " << rhoPF << " nNeutral " << nNeutral << " nCharged " << nCharged << endl;

  std::vector<float> v_pt_rho;
  v_pt_rho.push_back( pt );
  v_pt_rho.push_back( rhoPF );

  std::vector<float> v_nCharged;
  v_nCharged.push_back( (float)nCharged );

  std::vector<float> v_nNeutral;
  v_nNeutral.push_back( (float)nNeutral );

  std::vector<float> v_ptD;
  v_ptD.push_back( ptD );

  //  cout << "== pt " << v_pt_rho.at(0) << " rhoPF " << v_pt_rho.at(1) << " nNeutral " << nNeutral << " nCharged " << nCharged << endl;
  //  cout << " ciao " << sjc_nCharged_quark_<< endl;

  float quarkProb_nCharged = sjc_nCharged_quark_->correction(v_pt_rho, v_nCharged); 


  float gluonProb_nCharged = sjc_nCharged_gluon_->correction(v_pt_rho, v_nCharged);

  float quarkProb_nNeutral = sjc_nNeutral_quark_->correction(v_pt_rho, v_nNeutral);
  float gluonProb_nNeutral = sjc_nNeutral_gluon_->correction(v_pt_rho, v_nNeutral);

  float quarkProb_ptD = sjc_ptD_quark_->correction(v_pt_rho, v_ptD);
  float gluonProb_ptD = sjc_ptD_gluon_->correction(v_pt_rho, v_ptD);


  float quarkProb = quarkProb_nCharged*quarkProb_nNeutral*quarkProb_ptD;
  float gluonProb = gluonProb_nCharged*gluonProb_nNeutral*gluonProb_ptD;

  float QGLikelihood = (gluonProb+quarkProb>0.) ? quarkProb / (gluonProb + quarkProb ) : -1.;

  return QGLikelihood;

}
