#include <iostream>
#include <vector>

#include "Math/Point3D.h"
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"


/*
  Example of how to use this code:

//electron loop
for(unsigned int elidx = 0; elidx < els_p4().size(); elidx++) {

  vector<ConversionInfo> v_convInfos = getConversionInfos(elidx, evt_bField(), 0.45);	
  ConversionInfo blah = findBestConversionMatch(v_convInfos);
    cout << "Electron Pt: " << els_p4()[elidx].Pt() << endl;
    cout << "dist: " << blah.dist() << endl;
    cout << "dcot: " << blah.dcot() << endl;
    cout << "flag: " << blah.flag() << endl;
    cout << "rad. of conv" << blah.radiusOfConversion() << endl;
  }
  
*/


typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double> > XYZPoint;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

class ConversionInfo
{
	public:
		// constuct
		ConversionInfo();
		ConversionInfo
		( 
			double p_dist,
		 	double p_dcot,
		 	double p_radiusOfConversion,
		 	const XYZPoint& p_pointOfConversion,      
		 	int p_ctfPartnerIndex,
		 	int p_gsfPartnerIndex,
		 	int p_deltaMissingHits,
		 	int p_flag
		);

		// destroy
		~ConversionInfo();

		// members
		double dist() const {return dist_ ; }
		double dcot() const {return dcot_ ; }
		double radiusOfConversion() const { return radiusOfConversion_ ; }
		int ctfPartnerIndex() const { return ctfPartnerIndex_; } 
		int gsfPartnerIndex() const { return gsfPartnerIndex_; }
		XYZPoint pointOfConversion() const { return pointOfConversion_ ; }

		//if the partner track is found in the  CTF track collection,
		//we return a ref to the CTF partner track
		int deltaMissingHits() const { return deltaMissingHits_ ; }

		//if(flag == 0) //Partner track found in the CTF collection using the electron's CTF track
		//if(flag == 1) //Partner track found in the CTF collection using the electron's GSF track
		//if(flag == 2) //Partner track found in the GSF collection using the electron's CTF track
		//if(flag == 3) //Partner track found in the GSF collection using the electron's GSF track

		int flag() const { return flag_ ; }
	private:

		double dist_;
		double dcot_;
		double radiusOfConversion_;
		XYZPoint pointOfConversion_;
		int ctfPartnerIndex_;
		int gsfPartnerIndex_;
		int deltaMissingHits_;
		int flag_;

} ;





ConversionInfo getConversionInfo(const LorentzVector el_tk_p4, 
				 const int el_charge,
				 const float el_d0,
				 const float el_dz,
				 const LorentzVector cand_p4,
				 const int cand_charge,
				 const float cand_d0,
				 const double bFieldAtOrigin);


ConversionInfo arbitrateConversionPartnersbyR(const std::vector<ConversionInfo>& v_convCandidates);


ConversionInfo findBestConversionMatch(const std::vector<ConversionInfo>& v_convCandidates);

std::vector<ConversionInfo> getConversionInfos(const int gsfElectronIdx,				
					       const double bFieldAtOrigin,
					       const double minFracSharedHits);
