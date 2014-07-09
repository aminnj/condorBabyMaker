#include "SimpleFakeRate.h"
//*******************************************************
//
//  Usage of this code is documented in the header file
//
//*******************************************************

//-------------------------------------------------------------
// Class Constructor.  Includes some sanity checks.
//-------------------------------------------------------------
SimpleFakeRate::SimpleFakeRate(const char* filename, const char* histname, 
                               bool lastBinIsOverflow) {

  // fill the overflow flag
  overflowFlag_ = lastBinIsOverflow;

  // open the histogram file
  TFile infile(filename);
  if (infile.IsZombie()) {
    cout << "Error opening file " << filename << endl;
    gSystem->Exit(1);
  }
 
  // let's find the histogram
  // TH2F* thisHist = dynamic_cast<TH2F *>( infile.Get( histname) );
  hist_ = dynamic_cast<TH2F *>( infile.Get( histname) );
  if ( hist_ == 0) {
    cout << "Cannot find histogram " << histname << 
            " in file " << filename << endl;
    gSystem->Exit(1);
  }

  //prevent this histogram from being deleted once the file is closed
  hist_->SetDirectory(0);

  // now let's close the histogram file
  infile.Close();

  // now let's find the maximum and minimum pt and eta
  ptmax_  = hist_->GetYaxis()->GetBinUpEdge(hist_->GetYaxis()->GetNbins());
  etamax_ = hist_->GetXaxis()->GetBinUpEdge(hist_->GetXaxis()->GetNbins());
  ptmin_  = hist_->GetYaxis()->GetBinLowEdge(1);
  etamin_ = hist_->GetXaxis()->GetBinLowEdge(1);

  // if the maximum eta is negative complain
  if (etamax_ < 0.) {
    cout << "Histogram " << histname << " in file " << filename <<
      " has -ve maximum eta = " << etamax_ << " something is fishy" <<endl; 
  }


  // now find out if this looks like a fake rate for positive eta only
  if ( etamin_ > -0.001 && etamax_ > 0.) {
    //cout << "------" << endl;
    //cout << "The fake rate in histogram " << histname << endl;
    //cout << " in file " << filename << " is only for +ve eta" << endl;
    //cout << " .. we will assume that it is the same at -ve eta" << endl;
    //cout << "------" << endl;
    useAbsEta_ = true;
  } else {
    useAbsEta_ = false;
  }
}
//-------------------------------------------------------------------
// Real code that gets the fake rate (iflag=1) or the error (iflag=2)
//--------------------------------------------------------------------
float SimpleFakeRate::getValue(float pt, float eta, int iflag) {

  // handling of pt overflow
  float thispt = pt;
  if (thispt > ptmax_) {
    if (overflowFlag_) {
      thispt = ptmax_ - 0.001;
    } else {
      return 0.0;
    }
  }

  // handling of negative eta
  float thiseta = eta;
  if (useAbsEta_ && eta<0.) thiseta= -eta;

  // handling of overflow/underflow
  if (thispt < ptmin_ || thiseta > etamax_ || thiseta < etamin_) {
    cout << "------" << endl;
    cout << "Fake rate requested for eta = " << thiseta << " and pt = " << thispt <<endl;
    cout << "Parameters are out of range.  Return FR/FRerr=0.  Limits are:"<< endl;
    cout << "Pt:  " << ptmin_ << " to " << ptmax_  << endl;
    cout << "Eta: " << etamin_<< " to " << etamax_ << endl;
    cout << "------" << endl;
    return 0.0;
  }

  // get it
  float f = -1.;
  if (iflag == 1) {
    f  = hist_->GetBinContent(hist_->FindBin(thiseta,thispt));
  } else if (iflag == 2) {
    f  = hist_->GetBinError(hist_->FindBin(thiseta,thispt));
  }
  return f;
}
//-------------------------------------------
// Print it out
//--------------------------------------------
void SimpleFakeRate::print(){
  cout << "Fake rate print not implemented yet" << endl;
}
