
#include <iostream>
#include "electronSelectionsParameters.h"

/*
   Useful regexp incase you feel a little bit sick
   241,807s/cut\(.*\) = cms.vdouble(/double cut\1_tmp\[18\] = {/g

 */



void eidGetWP2012(const wp2012_tightness tightness, std::vector<double> &cutdeta, std::vector<double> &cutdphi, std::vector<double> &cuthoe, std::vector<double> &cutsee, std::vector<double> &cutooemoop, std::vector<double> &cutd0vtx, std::vector<double> &cutdzvtx, std::vector<bool> &cutvtxfit, std::vector<int> &cutmhit, std::vector<double> &cutrelisohighpt, std::vector<double> &cutrelisolowpt)
{

    switch (tightness) {
        case VETO:
            {
                double dEtaIn_tmp[2]        = {0.007, 0.010};
                double dPhiIn_tmp[2]        = {0.800, 0.700};
                double sigmaIEtaIEta_tmp[2] = {0.010, 0.030};
                double hoe_tmp[2]           = {0.150, 999.9};
                double ooemoop_tmp[2]       = {999.9, 999.9};
                double d0Vtx_tmp[2]         = {0.040, 0.040};
                double dzVtx_tmp[2]         = {0.200, 0.200};
                bool vtxFit_tmp[2]          = {false, false};
                int mHits_tmp[2]            = {999, 999};
                double isoHi_tmp[2]         = {0.150, 0.150};
                double isoLo_tmp[2]         = {0.150, 0.150};
                eidAssign(cutdeta,          dEtaIn_tmp, 2);
                eidAssign(cutdphi,          dPhiIn_tmp, 2);
                eidAssign(cutsee,           sigmaIEtaIEta_tmp, 2);
                eidAssign(cuthoe,           hoe_tmp, 2);
                eidAssign(cutooemoop,       ooemoop_tmp, 2);
                eidAssign(cutd0vtx,         d0Vtx_tmp, 2);
                eidAssign(cutdzvtx,         dzVtx_tmp, 2);
                eidAssign(cutvtxfit,        vtxFit_tmp, 2);
                eidAssign(cutmhit,          mHits_tmp, 2);
                eidAssign(cutrelisohighpt,  isoHi_tmp, 2);
                eidAssign(cutrelisolowpt,   isoLo_tmp, 2);
                return;
            }
        case LOOSE:
            {
                double dEtaIn_tmp[2]        = {0.007, 0.009};
                double dPhiIn_tmp[2]        = {0.150, 0.100};
                double sigmaIEtaIEta_tmp[2] = {0.010, 0.030};
                double hoe_tmp[2]           = {0.120, 0.100};
                double ooemoop_tmp[2]       = {0.050, 0.050};
                double d0Vtx_tmp[2]         = {0.020, 0.020};
                double dzVtx_tmp[2]         = {0.200, 0.200};
                bool vtxFit_tmp[2]          = {true, true};
                int mHits_tmp[2]            = {1, 1};
                double isoHi_tmp[2]         = {0.150, 0.150};
                double isoLo_tmp[2]         = {0.150, 0.100};
                eidAssign(cutdeta,          dEtaIn_tmp, 2);
                eidAssign(cutdphi,          dPhiIn_tmp, 2);
                eidAssign(cutsee,           sigmaIEtaIEta_tmp, 2);
                eidAssign(cuthoe,           hoe_tmp, 2);
                eidAssign(cutooemoop,       ooemoop_tmp, 2);
                eidAssign(cutd0vtx,         d0Vtx_tmp, 2);
                eidAssign(cutdzvtx,         dzVtx_tmp, 2);
                eidAssign(cutvtxfit,        vtxFit_tmp, 2);
                eidAssign(cutmhit,          mHits_tmp, 2);
                eidAssign(cutrelisohighpt,  isoHi_tmp, 2);
                eidAssign(cutrelisolowpt,   isoLo_tmp, 2);
                return;
            }
        case MEDIUM:
            {
                double dEtaIn_tmp[2]        = {0.004, 0.007};
                double dPhiIn_tmp[2]        = {0.060, 0.030};
                double sigmaIEtaIEta_tmp[2] = {0.010, 0.030};
                double hoe_tmp[2]           = {0.120, 0.100};
                double ooemoop_tmp[2]       = {0.050, 0.050};
                double d0Vtx_tmp[2]         = {0.020, 0.020};
                double dzVtx_tmp[2]         = {0.100, 0.100};
                bool vtxFit_tmp[2]          = {true, true};
                int mHits_tmp[2]            = {1, 1};
                double isoHi_tmp[2]         = {0.150, 0.150};
                double isoLo_tmp[2]         = {0.150, 0.100};
                eidAssign(cutdeta,          dEtaIn_tmp, 2);
                eidAssign(cutdphi,          dPhiIn_tmp, 2);
                eidAssign(cutsee,           sigmaIEtaIEta_tmp, 2);
                eidAssign(cuthoe,           hoe_tmp, 2);
                eidAssign(cutooemoop,       ooemoop_tmp, 2);
                eidAssign(cutd0vtx,         d0Vtx_tmp, 2);
                eidAssign(cutdzvtx,         dzVtx_tmp, 2);
                eidAssign(cutvtxfit,        vtxFit_tmp, 2);
                eidAssign(cutmhit,          mHits_tmp, 2);
                eidAssign(cutrelisohighpt,  isoHi_tmp, 2);
                eidAssign(cutrelisolowpt,   isoLo_tmp, 2);
                return;
            }
        case TIGHT:
            {
                double dEtaIn_tmp[2]        = {0.004, 0.005};
                double dPhiIn_tmp[2]        = {0.030, 0.020};
                double sigmaIEtaIEta_tmp[2] = {0.010, 0.030};
                double hoe_tmp[2]           = {0.120, 0.100};
                double ooemoop_tmp[2]       = {0.050, 0.050};
                double d0Vtx_tmp[2]         = {0.020, 0.020};
                double dzVtx_tmp[2]         = {0.100, 0.100};
                bool vtxFit_tmp[2]          = {true, true};
                int mHits_tmp[2]            = {0, 0};
                double isoHi_tmp[2]         = {0.100, 0.100};
                double isoLo_tmp[2]         = {0.100, 0.070};
                eidAssign(cutdeta,          dEtaIn_tmp, 2);
                eidAssign(cutdphi,          dPhiIn_tmp, 2);
                eidAssign(cutsee,           sigmaIEtaIEta_tmp, 2);
                eidAssign(cuthoe,           hoe_tmp, 2);
                eidAssign(cutooemoop,       ooemoop_tmp, 2);
                eidAssign(cutd0vtx,         d0Vtx_tmp, 2);
                eidAssign(cutdzvtx,         dzVtx_tmp, 2);
                eidAssign(cutvtxfit,        vtxFit_tmp, 2);
                eidAssign(cutmhit,          mHits_tmp, 2);
                eidAssign(cutrelisohighpt,  isoHi_tmp, 2);
                eidAssign(cutrelisolowpt,   isoLo_tmp, 2);
                return;
            }

        default:
            std::cout << "[eidGetWP2012] ERROR! Invalid tightness level" << std::endl;

    }

    return;

}

void eidGetCand(const cand_tightness tightness, std::vector<double> &cutdeta, std::vector<double> &cutdphi, std::vector<double> &cuthoe, std::vector<double> &cutslat)
{

    switch (tightness) {
        case CAND_01:
            {
                double dEtaInThresholds_tmp[2]               = {0.007, 0.010};
                double dPhiInThresholds_tmp[2]               = {0.020, 0.025};
                double hoeThresholds_tmp[2]                  = {0.01, 0.01};
                double latThresholds_tmp[2]                  = {0.90, 0.03};
                eidAssign(cutdeta, dEtaInThresholds_tmp, 2);
                eidAssign(cutdphi, dPhiInThresholds_tmp, 2);
                eidAssign(cuthoe, hoeThresholds_tmp, 2);
                eidAssign(cutslat, latThresholds_tmp, 2);
                return;
            }
        case CAND_02:
            {
                double dEtaInThresholds_tmp[2]               = {0.005, 0.007};
                double dPhiInThresholds_tmp[2]               = {0.020, 0.025};
                double hoeThresholds_tmp[2]                  = {0.01, 0.01};
                double latThresholds_tmp[2]                  = {0.94, 0.03};
                eidAssign(cutdeta, dEtaInThresholds_tmp, 2);
                eidAssign(cutdphi, dPhiInThresholds_tmp, 2);
                eidAssign(cuthoe, hoeThresholds_tmp, 2);
                eidAssign(cutslat, latThresholds_tmp, 2);
                return;
            }
        default:
            std::cout << "[eidGetCand] ERROR! Invalid tightness level" << std::endl;
    }

    return;
}


void eidGetVBTF(const vbtf_tightness tightness, std::vector<double> &cutdeta, std::vector<double> &cutdphi, std::vector<double> &cuthoe, std::vector<double> &cutsee, std::vector<double> &cutreliso)
{

    switch (tightness) {
        case VBTF_35X_95:
            {
                double isoThresholds_tmp[2]                 = {0.15,     0.1};
                double sigmaIEtaIEtaThresholds_tmp[2]       = {0.01,    0.03};
                double dPhiInThresholds_tmp[2]              = {0.8,     0.7};
                double dEtaInThresholds_tmp[2]              = {0.007,   0.01};
                double hoeThresholds_tmp[2]                 = {0.5,     0.07};
                eidAssign(cutreliso, isoThresholds_tmp, 2);
                eidAssign(cutdeta, dEtaInThresholds_tmp, 2);
                eidAssign(cutdphi, dPhiInThresholds_tmp, 2);
                eidAssign(cuthoe, hoeThresholds_tmp, 2);
                eidAssign(cutsee, sigmaIEtaIEtaThresholds_tmp, 2);
                return;
            }

        case VBTF_35X_90:
            {
                double isoThresholds_tmp[2]                 = {0.1,     0.07};
                double sigmaIEtaIEtaThresholds_tmp[2]       = {0.01,    0.03};
                double dPhiInThresholds_tmp[2]              = {0.8,     0.7};
                double dEtaInThresholds_tmp[2]              = {0.007,   0.009};
                double hoeThresholds_tmp[2]                 = {0.12,     0.05};
                eidAssign(cutreliso, isoThresholds_tmp, 2);
                eidAssign(cutdeta, dEtaInThresholds_tmp, 2);
                eidAssign(cutdphi, dPhiInThresholds_tmp, 2);
                eidAssign(cuthoe, hoeThresholds_tmp, 2);
                eidAssign(cutsee, sigmaIEtaIEtaThresholds_tmp, 2);
                return;
            }   

        case VBTF_35X_85:
            {
                double isoThresholds_tmp[2]                 = {0.09,     0.06};
                double sigmaIEtaIEtaThresholds_tmp[2]       = {0.01,    0.03};
                double dPhiInThresholds_tmp[2]              = {0.06,     0.04};
                double dEtaInThresholds_tmp[2]              = {0.006,   0.007};            
                double hoeThresholds_tmp[2]                 = {0.04,     0.025};
                eidAssign(cutreliso, isoThresholds_tmp, 2);
                eidAssign(cutdeta, dEtaInThresholds_tmp, 2);
                eidAssign(cutdphi, dPhiInThresholds_tmp, 2);
                eidAssign(cuthoe, hoeThresholds_tmp, 2);
                eidAssign(cutsee, sigmaIEtaIEtaThresholds_tmp, 2);
                return;
            }

        case VBTF_35X_80:
            {
                double isoThresholds_tmp[2]                 = {0.07,     0.06};
                double sigmaIEtaIEtaThresholds_tmp[2]       = {0.01,    0.03};
                double dPhiInThresholds_tmp[2]              = {0.06,     0.03};
                double dEtaInThresholds_tmp[2]              = {0.004,   0.007};            
                double hoeThresholds_tmp[2]                 = {0.04,     0.025};
                eidAssign(cutreliso, isoThresholds_tmp, 2);
                eidAssign(cutdeta, dEtaInThresholds_tmp, 2);
                eidAssign(cutdphi, dPhiInThresholds_tmp, 2);
                eidAssign(cuthoe, hoeThresholds_tmp, 2);
                eidAssign(cutsee, sigmaIEtaIEtaThresholds_tmp, 2);
                return;
            }

        case VBTF_35X_70:
            {
                double isoThresholds_tmp[2]                 = {0.05,     0.04};
                double sigmaIEtaIEtaThresholds_tmp[2]       = {0.01,    0.03};
                double dPhiInThresholds_tmp[2]              = {0.03,     0.02};
                double dEtaInThresholds_tmp[2]              = {0.003,   0.005};
                double hoeThresholds_tmp[2]                 = {0.025,     0.012};
                eidAssign(cutreliso, isoThresholds_tmp, 2);
                eidAssign(cutdeta, dEtaInThresholds_tmp, 2);
                eidAssign(cutdphi, dPhiInThresholds_tmp, 2);
                eidAssign(cuthoe, hoeThresholds_tmp, 2);
                eidAssign(cutsee, sigmaIEtaIEtaThresholds_tmp, 2);
                return;
            }

        case VBTF_35X_60:
            {
                double isoThresholds_tmp[2]                 = {0.04,     0.03};
                double sigmaIEtaIEtaThresholds_tmp[2]       = {0.01,    0.03};
                double dPhiInThresholds_tmp[2]              = {0.02,     0.02};
                double dEtaInThresholds_tmp[2]              = {0.0025,   0.003};
                double hoeThresholds_tmp[2]                 = {0.025,     0.009};
                eidAssign(cutreliso, isoThresholds_tmp, 2);
                eidAssign(cutdeta, dEtaInThresholds_tmp, 2);
                eidAssign(cutdphi, dPhiInThresholds_tmp, 2);
                eidAssign(cuthoe, hoeThresholds_tmp, 2);
                eidAssign(cutsee, sigmaIEtaIEtaThresholds_tmp, 2);
                return;
            }

        case VBTF_35Xr2_70:
            {
                double isoThresholds_tmp[2]                 = {0.04,     0.03};
                double sigmaIEtaIEtaThresholds_tmp[2]       = {0.01,    0.03};
                double dPhiInThresholds_tmp[2]              = {0.03,     0.02};
                double dEtaInThresholds_tmp[2]              = {0.004,   0.005};
                double hoeThresholds_tmp[2]                 = {0.025,     0.025};
                eidAssign(cutreliso, isoThresholds_tmp, 2);
                eidAssign(cutdeta, dEtaInThresholds_tmp, 2);
                eidAssign(cutdphi, dPhiInThresholds_tmp, 2);
                eidAssign(cuthoe, hoeThresholds_tmp, 2);
                eidAssign(cutsee, sigmaIEtaIEtaThresholds_tmp, 2);
                return;
            }

        case VBTF_35Xr2_60:
            {
                double isoThresholds_tmp[2]                 = {0.03,     0.02};
                double sigmaIEtaIEtaThresholds_tmp[2]       = {0.01,    0.03};
                double dPhiInThresholds_tmp[2]              = {0.025,     0.02};
                double dEtaInThresholds_tmp[2]              = {0.004,   0.005};
                double hoeThresholds_tmp[2]                 = {0.025,     0.025};
                eidAssign(cutreliso, isoThresholds_tmp, 2);
                eidAssign(cutdeta, dEtaInThresholds_tmp, 2);
                eidAssign(cutdphi, dPhiInThresholds_tmp, 2);
                eidAssign(cuthoe, hoeThresholds_tmp, 2);
                eidAssign(cutsee, sigmaIEtaIEtaThresholds_tmp, 2);
                return;
            }

        case VBTF_80_NOHOEEND:
            {
                double isoThresholds_tmp[2]                 = {0.07,     0.06};
                double sigmaIEtaIEtaThresholds_tmp[2]       = {0.01,    0.03};
                double dPhiInThresholds_tmp[2]              = {0.06,     0.03};
                double dEtaInThresholds_tmp[2]              = {0.004,   0.007};            
                double hoeThresholds_tmp[2]                 = {0.04,     9999.};
                eidAssign(cutreliso, isoThresholds_tmp, 2);
                eidAssign(cutdeta, dEtaInThresholds_tmp, 2);
                eidAssign(cutdphi, dPhiInThresholds_tmp, 2);
                eidAssign(cuthoe, hoeThresholds_tmp, 2);
                eidAssign(cutsee, sigmaIEtaIEtaThresholds_tmp, 2);
                return;
            }

        case VBTF_85_NOHOEEND:
            {
                double isoThresholds_tmp[2]                 = { 0.09  ,  0.06  };
                double sigmaIEtaIEtaThresholds_tmp[2]       = { 0.01  ,  0.03  };
                double dPhiInThresholds_tmp[2]              = { 0.06  ,  0.04  };
                double dEtaInThresholds_tmp[2]              = { 0.006 , 0.007  };
                double hoeThresholds_tmp[2]                 = { 0.04  , 9999.  };
                eidAssign(cutreliso, isoThresholds_tmp, 2);
                eidAssign(cutdeta, dEtaInThresholds_tmp, 2);
                eidAssign(cutdphi, dPhiInThresholds_tmp, 2);
                eidAssign(cuthoe, hoeThresholds_tmp, 2);
                eidAssign(cutsee, sigmaIEtaIEtaThresholds_tmp, 2);
                return;
            }

        case VBTF_85:
            {
                double isoThresholds_tmp[2]                 = { 0.09  ,  0.06  };
                double sigmaIEtaIEtaThresholds_tmp[2]       = { 0.01  ,  0.03  };
                double dPhiInThresholds_tmp[2]              = { 0.06  ,  0.04  };
                double dEtaInThresholds_tmp[2]              = { 0.006 , 0.007  };
                double hoeThresholds_tmp[2]                 = { 0.04  , 0.025  };
                eidAssign(cutreliso, isoThresholds_tmp, 2);
                eidAssign(cutdeta, dEtaInThresholds_tmp, 2);
                eidAssign(cutdphi, dPhiInThresholds_tmp, 2);
                eidAssign(cuthoe, hoeThresholds_tmp, 2);
                eidAssign(cutsee, sigmaIEtaIEtaThresholds_tmp, 2);
                return;
            }

        case VBTF_70_NOHOEEND:
            {
                double isoThresholds_tmp[2]                 = {0.04,     0.03};
                double sigmaIEtaIEtaThresholds_tmp[2]       = {0.01,    0.03};
                double dPhiInThresholds_tmp[2]              = {0.03,     0.02};
                double dEtaInThresholds_tmp[2]              = {0.004,   0.005};
                double hoeThresholds_tmp[2]                 = {0.025,     9999.};
                eidAssign(cutreliso, isoThresholds_tmp, 2);
                eidAssign(cutdeta, dEtaInThresholds_tmp, 2);
                eidAssign(cutdphi, dPhiInThresholds_tmp, 2);
                eidAssign(cuthoe, hoeThresholds_tmp, 2);
                eidAssign(cutsee, sigmaIEtaIEtaThresholds_tmp, 2);
                return;
            }
        case VBTF_90_HLT:
            {
                double isoThresholds_tmp[2]                 = {0.1,     0.07};
                double sigmaIEtaIEtaThresholds_tmp[2]       = {0.01,    0.03};
                double dPhiInThresholds_tmp[2]              = {0.15,    0.10};
                double dEtaInThresholds_tmp[2]              = {0.007,   0.009};
                double hoeThresholds_tmp[2]                 = {0.12,    0.10};
                eidAssign(cutreliso, isoThresholds_tmp, 2);
                eidAssign(cutdeta, dEtaInThresholds_tmp, 2);
                eidAssign(cutdphi, dPhiInThresholds_tmp, 2);
                eidAssign(cuthoe, hoeThresholds_tmp, 2);
                eidAssign(cutsee, sigmaIEtaIEtaThresholds_tmp, 2);
                return;
            }   

        case VBTF_90_HLT_CALOIDT_TRKIDVL:
            {
                double isoThresholds_tmp[2]                 = {0.1,     0.07};
                double sigmaIEtaIEtaThresholds_tmp[2]       = {0.01,    0.03};
                double dPhiInThresholds_tmp[2]              = {0.15,    0.10};
                double dEtaInThresholds_tmp[2]              = {0.007,   0.009};
                double hoeThresholds_tmp[2]                 = {0.10,    0.075};
                eidAssign(cutreliso, isoThresholds_tmp, 2);
                eidAssign(cutdeta, dEtaInThresholds_tmp, 2);
                eidAssign(cutdphi, dPhiInThresholds_tmp, 2);
                eidAssign(cuthoe, hoeThresholds_tmp, 2);
                eidAssign(cutsee, sigmaIEtaIEtaThresholds_tmp, 2);
                return;
            }   

        case VBTF_95_NOHOEEND:
            {
                double isoThresholds_tmp[2]                 = {0.15,    0.10};
                double sigmaIEtaIEtaThresholds_tmp[2]       = {0.01,    0.03};
                double dPhiInThresholds_tmp[2]              = {0.80,    0.70};
                double dEtaInThresholds_tmp[2]              = {0.007,   0.01};
                double hoeThresholds_tmp[2]                 = {0.15,    999.};
                eidAssign(cutreliso, isoThresholds_tmp, 2);
                eidAssign(cutdeta, dEtaInThresholds_tmp, 2);
                eidAssign(cutdphi, dPhiInThresholds_tmp, 2);
                eidAssign(cuthoe, hoeThresholds_tmp, 2);
                eidAssign(cutsee, sigmaIEtaIEtaThresholds_tmp, 2);
                return;
            }   

        default:
            std::cout << "[eidGetVBTF] ERROR! Invalid tightness level" << std::endl;
    }

    return;
}

void eidGetCIC_V06(const cic_tightness tightness, std::vector<double>& cutiso_sum, std::vector<double>& cutiso_sumoet, std::vector<double>& cuthoe, std::vector<double>& cutsee, std::vector<double>& cutdphiin, std::vector<double>& cutdetain, std::vector<double>& cuteseedopcor, std::vector<double>& cutfmishits, std::vector<double>& cutdcotdist, std::vector<double>& cutip_gsf, std::vector<double>& cutiso_sumoetl, std::vector<double> & cuthoel, std::vector<double>& cutseel, std::vector<double>& cutdphiinl, std::vector<double>& cutdetainl, std::vector<double>& cutip_gsfl)
{

    switch (tightness) {
        case CIC_VERYLOOSE:
            {   

                double cutdcotdist_tmp[9] = {3.89e-02, 3.90e-02, 3.96e-02, 3.92e-02, 3.95e-02, 3.97e-02, 3.92e-02, 3.95e-02, 2.98e-02};
                double cutdetain_tmp[9] = {1.33e-02, 6.97e-03, 2.43e-02, 2.47e-02, 5.50e-02, 2.20e-02, 4.31e-02, 3.84e-02, 3.13e-02};
                double cutdetainl_tmp[9] = { 1.28e-02, 5.93e-03, 2.64e-02, 2.72e-02, 6.72e-02, 2.06e-02, 1.92e-02, 1.97e-01, 2.91e-02};
                double cutdphiin_tmp[9] = {9.71e-02, 2.70e-01, 3.59e-01, 8.36e-02, 4.42e-01, 3.34e-01, 3.63e-01, 4.04e-01, 9.20e-01};
                double cutdphiinl_tmp[9] = {7.93e-02, 2.66e-01, 3.60e-01, 9.12e-02, 4.42e-01, 3.33e-01, 3.39e-01, 6.61e-01, 2.92e-01};
                double cuteseedopcor_tmp[9] = {6.35e-01, 3.27e-01, 4.00e-01, 7.31e-01, 3.50e-01, 4.54e-01, 1.27e-01, 2.91e-01, 6.28e-02};
                double cutfmishits_tmp[9] = {4.50e+00, 1.50e+00, 1.50e+00, 4.50e+00, 2.50e+00, 1.50e+00, 3.50e+00, 4.50e+00, 3.50e+00};
                double cuthoe_tmp[9] = {2.30e-01, 1.16e-01, 1.48e-01, 3.66e-01, 1.01e-01, 1.46e-01, 4.29e-01, 4.42e-01, 4.00e-01};
                double cuthoel_tmp[9] = {2.44e-01, 1.17e-01, 1.48e-01, 3.60e-01, 7.69e-02, 1.46e-01, 3.26e-01, 3.83e-01, 3.93e-01};
                double cutip_gsf_tmp[9] = {8.48e-02, 1.05e-01, 1.78e-01, 8.78e-02, 7.13e-01, 4.77e-01, 4.30e-01, 5.69e+00, 4.76e-01};
                double cutip_gsfl_tmp[9] = {8.70e-02, 1.09e-01, 1.79e-01, 7.55e-02, 7.14e-01, 5.24e-01, 9.01e-01, 1.84e+00, 3.01e-01};
                double cutiso_sum_tmp[9] = {2.56e+01, 1.70e+01, 1.76e+01, 1.86e+01, 8.79e+00, 1.25e+01, 2.14e+01, 2.34e+01, 3.23e+00};
                double cutiso_sumoet_tmp[9] = {5.38e+01, 1.07e+01, 1.03e+01, 4.02e+01, 5.81e+00, 8.01e+00, 9.27e+00, 1.15e+01, 8.86e+02};
                double cutiso_sumoetl_tmp[9] = {1.76e+01, 1.10e+01, 1.14e+01, 1.37e+01, 6.28e+00, 8.27e+00, 1.59e+01, 1.58e+01, 8.08e+00};
                double cutsee_tmp[9] = {1.57e-02, 1.20e-02, 1.84e-02, 3.98e-02, 3.24e-02, 3.81e-02, 1.25e-02, 6.42e-02, 6.69e-02};
                double cutseel_tmp[9] = {1.77e-02, 1.23e-02, 1.92e-02, 4.73e-02, 3.54e-02, 4.87e-02, 1.59e-02, 6.17e-02, 1.19e-01};

                eidAssign(cutdcotdist, cutdcotdist_tmp, 9);
                eidAssign(cutdetain, cutdetain_tmp, 9);
                eidAssign(cutdetainl, cutdetainl_tmp, 9);
                eidAssign(cutdphiin, cutdphiin_tmp, 9);
                eidAssign(cutdphiinl, cutdphiinl_tmp, 9);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 9);
                eidAssign(cutfmishits, cutfmishits_tmp, 9);
                eidAssign(cuthoe, cuthoe_tmp, 9);
                eidAssign(cuthoel, cuthoel_tmp, 9);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 9);
                eidAssign(cutip_gsfl, cutip_gsfl_tmp, 9);
                eidAssign(cutiso_sum, cutiso_sum_tmp, 9);
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 9);
                eidAssign(cutiso_sumoetl, cutiso_sumoetl_tmp, 9);
                eidAssign(cutsee, cutsee_tmp, 9);
                eidAssign(cutseel, cutseel_tmp, 9);

                return;
            }

        case CIC_LOOSE:
            {

                double cutdcotdist_tmp[9] = {3.87e-02, 3.50e-02, 3.18e-02, 3.92e-02, 3.94e-02, 3.97e-02, 3.10e-02, 3.95e-02, 1.10e-02};
                double cutdetain_tmp[9] = {1.33e-02, 5.28e-03, 1.44e-02, 2.19e-02, 1.25e-02, 1.37e-02, 2.18e-02, 3.84e-02, 2.75e-02};
                double cutdetainl_tmp[9] = {1.26e-02, 4.88e-03, 1.68e-02, 2.67e-02, 1.21e-02, 1.31e-02, 1.92e-02, 1.97e-01, 2.84e-02};
                double cutdphiin_tmp[9] = {9.36e-02, 2.46e-01, 3.25e-01, 8.18e-02, 3.22e-01, 2.83e-01, 3.54e-01, 4.04e-01, 6.80e-01};
                double cutdphiinl_tmp[9] = {7.93e-02, 2.44e-01, 3.11e-01, 9.12e-02, 3.04e-01, 2.82e-01, 3.39e-01, 6.61e-01, 2.91e-01};
                double cuteseedopcor_tmp[9] = {6.37e-01, 8.79e-01, 4.02e-01, 7.45e-01, 3.67e-01, 4.88e-01, 1.27e-01, 7.19e-01, 6.28e-02};
                double cutfmishits_tmp[9] = {4.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 3.50e+00, 3.50e+00, 3.50e+00};
                double cuthoe_tmp[9] = {1.96e-01, 7.92e-02, 1.48e-01, 3.66e-01, 6.88e-02, 1.45e-01, 4.29e-01, 4.42e-01, 4.00e-01};
                double cuthoel_tmp[9] = {2.26e-01, 7.95e-02, 1.48e-01, 3.60e-01, 6.23e-02, 1.46e-01, 3.26e-01, 3.83e-01, 3.92e-01};
                double cutip_gsf_tmp[9] = {8.48e-02, 9.95e-02, 1.75e-01, 6.97e-02, 5.65e-01, 4.77e-01, 4.30e-01, 3.32e+00, 1.61e-01};
                double cutip_gsfl_tmp[9] = {7.58e-02, 9.81e-02, 1.76e-01, 6.66e-02, 5.65e-01, 5.16e-01, 9.01e-01, 1.12e+00, 8.42e-02};
                double cutiso_sum_tmp[9] = {2.02e+01, 1.31e+01, 1.56e+01, 1.61e+01, 8.61e+00, 1.10e+01, 1.31e+01, 1.63e+01, 2.37e+00};
                double cutiso_sumoet_tmp[9] = {1.49e+01, 8.33e+00, 7.64e+00, 1.22e+01, 4.49e+00, 5.59e+00, 7.44e+00, 7.31e+00, 2.74e+01};
                double cutiso_sumoetl_tmp[9] = {1.21e+01, 9.18e+00, 8.66e+00, 9.43e+00, 4.34e+00, 5.73e+00, 1.08e+01, 1.05e+01, 8.08e+00};
                double cutsee_tmp[9] = {1.57e-02, 1.12e-02, 1.40e-02, 3.95e-02, 3.10e-02, 3.37e-02, 1.11e-02, 6.13e-02, 6.69e-02};
                double cutseel_tmp[9] = {1.77e-02, 1.15e-02, 1.50e-02, 4.55e-02, 3.24e-02, 4.46e-02, 1.22e-02, 6.17e-02, 1.19e-01};

                eidAssign(cutdcotdist, cutdcotdist_tmp, 9);
                eidAssign(cutdetain, cutdetain_tmp, 9);
                eidAssign(cutdetainl, cutdetainl_tmp, 9);
                eidAssign(cutdphiin, cutdphiin_tmp, 9);
                eidAssign(cutdphiinl, cutdphiinl_tmp, 9);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 9);
                eidAssign(cutfmishits, cutfmishits_tmp, 9);
                eidAssign(cuthoe, cuthoe_tmp, 9);
                eidAssign(cuthoel, cuthoel_tmp, 9);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 9);
                eidAssign(cutip_gsfl, cutip_gsfl_tmp, 9);
                eidAssign(cutiso_sum, cutiso_sum_tmp, 9);
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 9);
                eidAssign(cutiso_sumoetl, cutiso_sumoetl_tmp, 9);
                eidAssign(cutsee, cutsee_tmp, 9);
                eidAssign(cutseel, cutseel_tmp, 9);

                return;
            }

        case CIC_MEDIUM:
            {

                double cutdcotdist_tmp[9] = {3.32e-02, 2.92e-02, 2.49e-02, 3.92e-02, 3.41e-02, 3.96e-02, 2.91e-02, 3.95e-02, 7.71e-03};
                double cutdetain_tmp[9] = {1.33e-02, 4.48e-03, 9.22e-03, 1.54e-02, 7.26e-03, 1.24e-02, 1.29e-02, 3.84e-02, 1.88e-02};
                double cutdetainl_tmp[9] = {1.21e-02, 4.22e-03, 9.18e-03, 1.61e-02, 6.45e-03, 1.16e-02, 1.23e-02, 6.20e-02, 2.43e-02};
                double cutdphiin_tmp[9] = {7.09e-02, 2.43e-01, 2.96e-01, 7.98e-02, 2.35e-01, 2.76e-01, 3.42e-01, 4.04e-01, 2.99e-01};
                double cutdphiinl_tmp[9] = {7.42e-02, 2.43e-01, 2.97e-01, 9.12e-02, 2.26e-01, 2.76e-01, 3.34e-01, 5.58e-01, 2.91e-01};
                double cuteseedopcor_tmp[9] = {6.42e-01, 9.44e-01, 4.53e-01, 7.62e-01, 3.67e-01, 5.57e-01, 1.98e-01, 9.15e-01, 6.28e-02};
                double cutfmishits_tmp[9] = {4.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01};
                double cuthoe_tmp[9] = {1.96e-01, 6.30e-02, 1.48e-01, 3.66e-01, 5.66e-02, 1.45e-01, 4.29e-01, 4.28e-01, 3.99e-01};
                double cuthoel_tmp[9] = {2.19e-01, 6.19e-02, 1.47e-01, 3.58e-01, 4.61e-02, 1.46e-01, 3.26e-01, 3.81e-01, 3.89e-01};
                double cutip_gsf_tmp[9] = {2.45e-02, 9.74e-02, 1.48e-01, 5.49e-02, 5.65e-01, 3.33e-01, 2.04e-01, 5.41e-01, 1.21e-01};
                double cutip_gsfl_tmp[9] = {1.92e-02, 9.81e-02, 1.33e-01, 4.34e-02, 5.65e-01, 3.24e-01, 2.33e-01, 4.30e-01, 6.44e-02};
                double cutiso_sum_tmp[9] = {1.44e+01, 1.12e+01, 1.09e+01, 1.08e+01, 6.35e+00, 9.78e+00, 1.30e+01, 1.62e+01, 1.96e+00};
                double cutiso_sumoet_tmp[9] = {1.01e+01, 6.41e+00, 6.00e+00, 8.14e+00, 3.90e+00, 4.76e+00, 6.86e+00, 6.48e+00, 1.74e+01};
                double cutiso_sumoetl_tmp[9] = {9.44e+00, 7.67e+00, 7.15e+00, 7.34e+00, 3.35e+00, 4.70e+00, 8.32e+00, 7.55e+00, 6.25e+00};
                double cutsee_tmp[9] = {1.30e-02, 1.09e-02, 1.18e-02, 3.94e-02, 3.04e-02, 3.28e-02, 1.00e-02, 3.73e-02, 6.69e-02};
                double cutseel_tmp[9] = {1.42e-02, 1.11e-02, 1.29e-02, 4.32e-02, 2.96e-02, 3.82e-02, 1.01e-02, 4.45e-02, 1.19e-01};


                eidAssign(cutdcotdist, cutdcotdist_tmp, 9);
                eidAssign(cutdetain, cutdetain_tmp, 9);
                eidAssign(cutdetainl, cutdetainl_tmp, 9);
                eidAssign(cutdphiin, cutdphiin_tmp, 9);
                eidAssign(cutdphiinl, cutdphiinl_tmp, 9);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 9);
                eidAssign(cutfmishits, cutfmishits_tmp, 9);
                eidAssign(cuthoe, cuthoe_tmp, 9);
                eidAssign(cuthoel, cuthoel_tmp, 9);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 9);
                eidAssign(cutip_gsfl, cutip_gsfl_tmp, 9);
                eidAssign(cutiso_sum, cutiso_sum_tmp, 9);
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 9);
                eidAssign(cutiso_sumoetl, cutiso_sumoetl_tmp, 9);
                eidAssign(cutsee, cutsee_tmp, 9);
                eidAssign(cutseel, cutseel_tmp, 9);

                return;
            }

        case CIC_TIGHT:
            {

                double cutdcotdist_tmp[9] = {2.68e-02, 2.36e-02, 2.21e-02, 3.72e-02, 3.17e-02, 3.61e-02, 2.55e-02, 3.75e-02, 2.16e-04};
                double cutdetain_tmp[9] = {8.92e-03, 3.96e-03, 8.50e-03, 1.34e-02, 6.27e-03, 1.05e-02, 1.12e-02, 3.09e-02, 1.88e-02};
                double cutdetainl_tmp[9] = {9.23e-03, 3.77e-03, 8.70e-03, 1.39e-02, 5.60e-03, 9.40e-03, 1.07e-02, 6.20e-02, 4.10e-03};
                double cutdphiin_tmp[9] = {6.37e-02, 1.53e-01, 2.90e-01, 7.69e-02, 1.81e-01, 2.34e-01, 3.42e-01, 3.93e-01, 2.84e-01};
                double cutdphiinl_tmp[9] = {6.92e-02, 2.33e-01, 2.96e-01, 8.65e-02, 1.85e-01, 2.76e-01, 3.34e-01, 3.53e-01, 2.90e-01};
                double cuteseedopcor_tmp[9] = {6.52e-01, 9.69e-01, 9.12e-01, 7.79e-01, 3.67e-01, 6.99e-01, 3.28e-01, 9.67e-01, 5.89e-01};
                double cutfmishits_tmp[9] = {4.50e+00, 1.50e+00, 5.00e-01, 1.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
                double cuthoe_tmp[9] = {1.74e-01, 4.88e-02, 1.46e-01, 3.64e-01, 4.93e-02, 1.45e-01, 4.29e-01, 4.20e-01, 3.99e-01};
                double cuthoel_tmp[9] = {2.19e-01, 5.25e-02, 1.47e-01, 3.57e-01, 4.25e-02, 1.45e-01, 3.26e-01, 3.80e-01, 1.32e-01};
                double cutip_gsf_tmp[9] = {1.58e-02, 8.25e-02, 1.15e-01, 4.05e-02, 5.40e-01, 1.51e-01, 7.74e-02, 4.17e-01, 7.80e-02};
                double cutip_gsfl_tmp[9] = {1.27e-02, 6.26e-02, 9.68e-02, 3.02e-02, 5.65e-01, 1.46e-01, 7.90e-02, 4.10e-01, 4.79e-02};
                double cutiso_sum_tmp[9] = {1.23e+01, 9.77e+00, 1.01e+01, 9.77e+00, 6.13e+00, 7.55e+00, 1.30e+01, 1.62e+01, 1.78e+00};
                double cutiso_sumoet_tmp[9] = {7.75e+00, 5.45e+00, 5.67e+00, 5.97e+00, 3.17e+00, 3.86e+00, 6.06e+00, 5.31e+00, 1.05e+01};
                double cutiso_sumoetl_tmp[9] = {7.56e+00, 5.08e+00, 5.77e+00, 5.74e+00, 2.37e+00, 3.32e+00, 4.97e+00, 5.46e+00, 3.82e+00};
                double cutsee_tmp[9] = {1.16e-02, 1.07e-02, 1.08e-02, 3.49e-02, 2.89e-02, 3.08e-02, 9.87e-03, 3.37e-02, 4.40e-02};
                double cutseel_tmp[9] = {1.27e-02, 1.08e-02, 1.13e-02, 4.19e-02, 2.81e-02, 3.02e-02, 9.76e-03, 4.28e-02, 2.98e-02};

                eidAssign(cutdcotdist, cutdcotdist_tmp, 9);
                eidAssign(cutdetain, cutdetain_tmp, 9);
                eidAssign(cutdetainl, cutdetainl_tmp, 9);
                eidAssign(cutdphiin, cutdphiin_tmp, 9);
                eidAssign(cutdphiinl, cutdphiinl_tmp, 9);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 9);
                eidAssign(cutfmishits, cutfmishits_tmp, 9);
                eidAssign(cuthoe, cuthoe_tmp, 9);
                eidAssign(cuthoel, cuthoel_tmp, 9);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 9);
                eidAssign(cutip_gsfl, cutip_gsfl_tmp, 9);
                eidAssign(cutiso_sum, cutiso_sum_tmp, 9);
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 9);
                eidAssign(cutiso_sumoetl, cutiso_sumoetl_tmp, 9);
                eidAssign(cutsee, cutsee_tmp, 9);
                eidAssign(cutseel, cutseel_tmp, 9);

                return;
            }

        case CIC_SUPERTIGHT:
            {

                double cutdcotdist_tmp[9] = {2.11e-02, 1.86e-02, 1.55e-02, 3.40e-02, 2.85e-02, 3.32e-02, 1.64e-02, 3.75e-02, 1.30e-04};
                double cutdetain_tmp[9] = {7.84e-03, 3.67e-03, 7.00e-03, 1.28e-02, 5.65e-03, 9.53e-03, 1.08e-02, 2.97e-02, 7.24e-03};
                double cutdetainl_tmp[9] = {7.61e-03, 3.28e-03, 6.57e-03, 1.03e-02, 5.05e-03, 8.55e-03, 1.07e-02, 2.94e-02, 4.10e-03};
                double cutdphiin_tmp[9] = {4.83e-02, 7.39e-02, 2.38e-01, 5.74e-02, 1.29e-01, 2.13e-01, 3.31e-01, 3.93e-01, 2.84e-01};
                double cutdphiinl_tmp[9] = {5.79e-02, 7.21e-02, 2.18e-01, 7.70e-02, 1.41e-01, 2.11e-01, 2.43e-01, 3.53e-01, 2.89e-01};
                double cuteseedopcor_tmp[9] = {7.32e-01, 9.77e-01, 9.83e-01, 8.55e-01, 4.31e-01, 7.35e-01, 4.18e-01, 9.99e-01, 5.89e-01};
                double cutfmishits_tmp[9] = {3.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
                double cuthoe_tmp[9] = {9.19e-02, 4.11e-02, 1.42e-01, 3.35e-01, 3.82e-02, 1.41e-01, 4.29e-01, 4.01e-01, 3.99e-01};
                double cuthoel_tmp[9] = {7.51e-02, 3.81e-02, 1.41e-01, 3.32e-01, 3.10e-02, 1.43e-01, 2.35e-01, 3.80e-01, 1.32e-01};
                double cutip_gsf_tmp[9] = {1.42e-02, 2.66e-02, 1.06e-01, 3.38e-02, 3.23e-01, 1.07e-01, 7.74e-02, 2.32e-01, 7.80e-02};
                double cutip_gsfl_tmp[9] = {1.15e-02, 2.72e-02, 8.41e-02, 2.49e-02, 4.17e-01, 1.02e-01, 7.90e-02, 1.69e-01, 4.79e-02};
                double cutiso_sum_tmp[9] = {8.95e+00, 8.18e+00, 8.75e+00, 7.47e+00, 5.43e+00, 5.87e+00, 8.16e+00, 1.02e+01, 1.78e+00};
                double cutiso_sumoet_tmp[9] = {6.45e+00, 5.14e+00, 4.99e+00, 5.21e+00, 2.65e+00, 3.12e+00, 4.52e+00, 4.72e+00, 3.68e+00};
                double cutiso_sumoetl_tmp[9] = {6.02e+00, 3.96e+00, 4.23e+00, 4.73e+00, 1.99e+00, 2.64e+00, 3.72e+00, 3.81e+00, 1.44e+00};
                double cutsee_tmp[9] = {1.09e-02, 1.05e-02, 1.05e-02, 3.24e-02, 2.81e-02, 2.95e-02, 9.77e-03, 2.75e-02, 2.95e-02};
                double cutseel_tmp[9] = {1.12e-02, 1.05e-02, 1.07e-02, 3.51e-02, 2.75e-02, 2.87e-02, 9.59e-03, 2.67e-02, 2.98e-02};

                eidAssign(cutdcotdist, cutdcotdist_tmp, 9);
                eidAssign(cutdetain, cutdetain_tmp, 9);
                eidAssign(cutdetainl, cutdetainl_tmp, 9);
                eidAssign(cutdphiin, cutdphiin_tmp, 9);
                eidAssign(cutdphiinl, cutdphiinl_tmp, 9);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 9);
                eidAssign(cutfmishits, cutfmishits_tmp, 9);
                eidAssign(cuthoe, cuthoe_tmp, 9);
                eidAssign(cuthoel, cuthoel_tmp, 9);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 9);
                eidAssign(cutip_gsfl, cutip_gsfl_tmp, 9);
                eidAssign(cutiso_sum, cutiso_sum_tmp, 9);
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 9);
                eidAssign(cutiso_sumoetl, cutiso_sumoetl_tmp, 9);
                eidAssign(cutsee, cutsee_tmp, 9);
                eidAssign(cutseel, cutseel_tmp, 9);

                return;
            }

        case CIC_HYPERTIGHT1:
            {

                double cutdcotdist_tmp[9] = {1.48e-02, 1.50e-02, 8.25e-03, 3.16e-02, 2.85e-02, 3.15e-02, 6.62e-03, 3.48e-02, 3.63e-06};
                double cutdetain_tmp[9] = {6.51e-03, 3.51e-03, 5.53e-03, 9.16e-03, 5.30e-03, 8.28e-03, 1.08e-02, 2.97e-02, 7.24e-03};
                double cutdetainl_tmp[9] = {6.05e-03, 3.23e-03, 4.93e-03, 8.01e-03, 4.93e-03, 7.91e-03, 1.03e-02, 2.94e-02, 4.10e-03};
                double cutdphiin_tmp[9] = {4.83e-02, 4.91e-02, 2.30e-01, 3.48e-02, 7.44e-02, 2.04e-01, 9.95e-02, 3.93e-01, 2.84e-01};
                double cutdphiinl_tmp[9] = {4.74e-02, 4.51e-02, 2.18e-01, 2.99e-02, 7.37e-02, 2.11e-01, 9.99e-02, 3.53e-01, 2.89e-01};
                double cuteseedopcor_tmp[9] = {7.72e-01, 9.90e-01, 1.01e+00, 8.55e-01, 9.11e-01, 7.72e-01, 9.17e-01, 1.06e+00, 7.63e-01};
                double cutfmishits_tmp[9] = {3.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
                double cuthoe_tmp[9] = {6.17e-02, 3.70e-02, 1.41e-01, 2.91e-01, 3.82e-02, 1.34e-01, 4.19e-01, 3.87e-01, 3.93e-01};
                double cuthoel_tmp[9] = {4.43e-02, 3.57e-02, 1.41e-01, 2.81e-01, 3.07e-02, 1.28e-01, 2.27e-01, 3.80e-01, 1.32e-01};
                double cutip_gsf_tmp[9] = {1.21e-02, 1.76e-02, 6.01e-02, 2.96e-02, 1.74e-01, 9.70e-02, 7.74e-02, 1.33e-01, 7.80e-02};
                double cutip_gsfl_tmp[9] = {1.01e-02, 1.56e-02, 6.87e-02, 2.13e-02, 1.25e-01, 8.16e-02, 7.90e-02, 1.30e-01, 4.79e-02};
                double cutiso_sum_tmp[9] = {7.92e+00, 6.85e+00, 7.87e+00, 6.77e+00, 4.47e+00, 5.28e+00, 6.57e+00, 1.02e+01, 1.78e+00};
                double cutiso_sumoet_tmp[9] = {5.20e+00, 3.93e+00, 3.88e+00, 4.10e+00, 2.40e+00, 2.43e+00, 3.49e+00, 3.94e+00, 3.01e+00};
                double cutiso_sumoetl_tmp[9] = {4.18e+00, 3.12e+00, 3.44e+00, 3.25e+00, 1.77e+00, 2.06e+00, 2.83e+00, 3.12e+00, 1.43e+00};
                double cutsee_tmp[9] = {1.05e-02, 1.04e-02, 1.01e-02, 3.24e-02, 2.80e-02, 2.85e-02, 9.67e-03, 2.61e-02, 2.95e-02};
                double cutseel_tmp[9] = {1.04e-02, 1.03e-02, 1.01e-02, 3.04e-02, 2.74e-02, 2.78e-02, 9.58e-03, 2.54e-02, 2.83e-02};

                eidAssign(cutdcotdist, cutdcotdist_tmp, 9);
                eidAssign(cutdetain, cutdetain_tmp, 9);
                eidAssign(cutdetainl, cutdetainl_tmp, 9);
                eidAssign(cutdphiin, cutdphiin_tmp, 9);
                eidAssign(cutdphiinl, cutdphiinl_tmp, 9);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 9);
                eidAssign(cutfmishits, cutfmishits_tmp, 9);
                eidAssign(cuthoe, cuthoe_tmp, 9);
                eidAssign(cuthoel, cuthoel_tmp, 9);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 9);
                eidAssign(cutip_gsfl, cutip_gsfl_tmp, 9);
                eidAssign(cutiso_sum, cutiso_sum_tmp, 9);
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 9);
                eidAssign(cutiso_sumoetl, cutiso_sumoetl_tmp, 9);
                eidAssign(cutsee, cutsee_tmp, 9);
                eidAssign(cutseel, cutseel_tmp, 9);

                return;
            }

        case CIC_HYPERTIGHT2:
            {

                double cutdcotdist_tmp[9] = {1.15e-02, 1.07e-02, 4.01e-03, 2.97e-02, 2.85e-02, 3.10e-02, 9.34e-04, 3.40e-02, 2.82e-07};
                double cutdetain_tmp[9] = {5.29e-03, 2.56e-03, 4.89e-03, 7.89e-03, 5.30e-03, 7.37e-03, 8.91e-03, 9.36e-03, 5.94e-03};
                double cutdetainl_tmp[9] = {4.48e-03, 2.59e-03, 4.42e-03, 6.54e-03, 4.93e-03, 6.98e-03, 8.49e-03, 9.06e-03, -4.81e-03};
                double cutdphiin_tmp[9] = {2.41e-02, 3.83e-02, 1.48e-01, 2.91e-02, 3.15e-02, 1.57e-01, 8.90e-02, 1.02e-01, 2.81e-01};
                double cutdphiinl_tmp[9] = {2.13e-02, 3.79e-02, 1.25e-01, 2.24e-02, 3.69e-02, 1.64e-01, 9.99e-02, 9.23e-02, 2.37e-01};
                double cuteseedopcor_tmp[9] = {1.03e+00, 9.95e-01, 1.03e+00, 1.01e+00, 9.46e-01, 9.03e-01, 9.97e-01, 1.14e+00, 8.00e-01};
                double cutfmishits_tmp[9] = {1.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01};
                double cuthoe_tmp[9] = {4.94e-02, 3.45e-02, 1.40e-01, 2.02e-01, 3.82e-02, 1.19e-01, 1.23e-01, 3.82e-01, 2.50e-01};
                double cuthoel_tmp[9] = {4.04e-02, 3.42e-02, 1.31e-01, 1.85e-01, 3.01e-02, 1.27e-01, 2.27e-01, 3.80e-01, 1.32e-01};
                double cutip_gsf_tmp[9] = {1.14e-02, 1.38e-02, 5.29e-02, 1.87e-02, 1.31e-01, 8.63e-02, 7.74e-02, 1.04e-01, 2.42e-02};
                double cutip_gsfl_tmp[9] = {9.83e-03, 1.35e-02, 4.27e-02, 1.72e-02, 1.25e-01, 7.92e-02, 7.90e-02, 1.30e-01, 3.40e-02};
                double cutiso_sum_tmp[9] = {6.40e+00, 5.77e+00, 6.54e+00, 5.22e+00, 3.86e+00, 4.63e+00, 6.31e+00, 1.02e+01, 1.78e+00};
                double cutiso_sumoet_tmp[9] = {4.03e+00, 3.03e+00, 3.24e+00, 3.13e+00, 2.05e+00, 2.01e+00, 2.99e+00, 3.44e+00, 2.76e+00};
                double cutiso_sumoetl_tmp[9] = {3.08e+00, 2.31e+00, 2.84e+00, 2.53e+00, 1.65e+00, 1.72e+00, 2.34e+00, 3.11e+00, 1.35e+00};
                double cutsee_tmp[9] = {1.03e-02, 1.03e-02, 9.88e-03, 3.03e-02, 2.79e-02, 2.79e-02, 9.67e-03, 2.52e-02, 2.58e-02};
                double cutseel_tmp[9] = {1.02e-02, 1.02e-02, 9.80e-03, 2.90e-02, 2.74e-02, 2.75e-02, 9.58e-03, 2.49e-02, 2.50e-02};

                eidAssign(cutdcotdist, cutdcotdist_tmp, 9);
                eidAssign(cutdetain, cutdetain_tmp, 9);
                eidAssign(cutdetainl, cutdetainl_tmp, 9);
                eidAssign(cutdphiin, cutdphiin_tmp, 9);
                eidAssign(cutdphiinl, cutdphiinl_tmp, 9);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 9);
                eidAssign(cutfmishits, cutfmishits_tmp, 9);
                eidAssign(cuthoe, cuthoe_tmp, 9);
                eidAssign(cuthoel, cuthoel_tmp, 9);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 9);
                eidAssign(cutip_gsfl, cutip_gsfl_tmp, 9);
                eidAssign(cutiso_sum, cutiso_sum_tmp, 9);
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 9);
                eidAssign(cutiso_sumoetl, cutiso_sumoetl_tmp, 9);
                eidAssign(cutsee, cutsee_tmp, 9);
                eidAssign(cutseel, cutseel_tmp, 9);

                return;
            }

        case CIC_HYPERTIGHT3:
            {

                double cutdcotdist_tmp[9] = {9.63e-03, 5.11e-03, 1.95e-04, 2.97e-02, 2.85e-02, 2.18e-02, 2.61e-05, 2.57e-02, 2.82e-07};
                double cutdetain_tmp[9] = {4.86e-03, 2.29e-03, 4.40e-03, 7.79e-03, 4.07e-03, 6.33e-03, 7.70e-03, 7.93e-03, 5.94e-03};
                double cutdetainl_tmp[9] = {4.48e-03, 2.30e-03, 4.14e-03, 6.04e-03, 3.87e-03, 6.09e-03, 7.97e-03, 8.04e-03, -4.81e-03};
                double cutdphiin_tmp[9] = {2.41e-02, 2.88e-02, 7.39e-02, 2.91e-02, 1.91e-02, 1.14e-01, 3.61e-02, 8.92e-02, 2.81e-01};
                double cutdphiinl_tmp[9] = {1.95e-02, 3.42e-02, 8.06e-02, 2.22e-02, 2.26e-02, 9.73e-02, 4.51e-02, 9.23e-02, 2.37e-01};
                double cuteseedopcor_tmp[9] = {1.07e+00, 1.01e+00, 1.08e+00, 1.01e+00, 9.69e-01, 9.10e-01, 1.04e+00, 1.20e+00, 8.00e-01};
                double cutfmishits_tmp[9] = {5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01};
                double cuthoe_tmp[9] = {3.52e-02, 3.45e-02, 1.33e-01, 1.88e-01, 2.72e-02, 1.19e-01, 9.28e-02, 2.46e-01, 2.50e-01};
                double cuthoel_tmp[9] = {4.04e-02, 3.40e-02, 1.31e-01, 1.84e-01, 2.64e-02, 1.18e-01, 9.76e-02, 2.53e-01, 1.32e-01};
                double cutip_gsf_tmp[9] = {1.14e-02, 1.26e-02, 3.79e-02, 1.68e-02, 1.21e-01, 5.29e-02, 7.74e-02, 3.35e-02, 2.42e-02};
                double cutip_gsfl_tmp[9] = {9.83e-03, 1.18e-02, 3.59e-02, 1.56e-02, 1.20e-01, 5.36e-02, 7.90e-02, 2.88e-02, 3.40e-02};
                double cutiso_sum_tmp[9] = {5.40e+00, 5.41e+00, 5.88e+00, 4.32e+00, 3.86e+00, 4.33e+00, 5.87e+00, 9.05e+00, 1.78e+00};
                double cutiso_sumoet_tmp[9] = {3.03e+00, 2.50e+00, 2.58e+00, 2.44e+00, 1.91e+00, 1.76e+00, 2.92e+00, 3.13e+00, 2.76e+00};
                double cutiso_sumoetl_tmp[9] = {2.36e+00, 2.02e+00, 2.29e+00, 1.89e+00, 1.65e+00, 1.69e+00, 2.03e+00, 2.79e+00, 1.35e+00};
                double cutsee_tmp[9] = {1.03e-02, 1.01e-02, 9.84e-03, 2.89e-02, 2.74e-02, 2.73e-02, 9.47e-03, 2.44e-02, 2.58e-02};
                double cutseel_tmp[9] = {1.02e-02, 1.00e-02, 9.73e-03, 2.79e-02, 2.73e-02, 2.69e-02, 9.40e-03, 2.46e-02, 2.50e-02};

                eidAssign(cutdcotdist, cutdcotdist_tmp, 9);
                eidAssign(cutdetain, cutdetain_tmp, 9);
                eidAssign(cutdetainl, cutdetainl_tmp, 9);
                eidAssign(cutdphiin, cutdphiin_tmp, 9);
                eidAssign(cutdphiinl, cutdphiinl_tmp, 9);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 9);
                eidAssign(cutfmishits, cutfmishits_tmp, 9);
                eidAssign(cuthoe, cuthoe_tmp, 9);
                eidAssign(cuthoel, cuthoel_tmp, 9);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 9);
                eidAssign(cutip_gsfl, cutip_gsfl_tmp, 9);
                eidAssign(cutiso_sum, cutiso_sum_tmp, 9);
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 9);
                eidAssign(cutiso_sumoetl, cutiso_sumoetl_tmp, 9);
                eidAssign(cutsee, cutsee_tmp, 9);
                eidAssign(cutseel, cutseel_tmp, 9);

                return;
            }

        case CIC_HYPERTIGHT4:
            {

                double cutdcotdist_tmp[9] = {2.70e-04, 1.43e-04, 1.95e-04, 2.64e-03, 2.82e-02, 1.64e-02, 2.61e-05, 2.57e-02, 2.82e-07};
                double cutdetain_tmp[9] = {2.44e-03, 1.67e-03, 2.26e-03, 3.43e-03, 3.51e-03, 3.52e-03, 2.98e-03, 4.79e-03, 5.94e-03};
                double cutdetainl_tmp[9] = {2.34e-03, 1.29e-03, 2.30e-03, 3.30e-03, 3.61e-03, 3.84e-03, 2.53e-03, 3.66e-03, -4.81e-03};
                double cutdphiin_tmp[9] = {8.44e-03, 5.21e-03, 2.18e-02, 1.39e-02, 7.82e-03, 1.52e-02, 2.59e-02, 3.87e-02, 2.81e-01};
                double cutdphiinl_tmp[9] = {5.77e-03, 3.20e-03, 2.85e-02, 2.22e-02, 7.00e-03, 1.84e-02, 2.91e-02, 4.40e-02, 2.37e-01};
                double cuteseedopcor_tmp[9] = {1.15e+00, 1.01e+00, 1.21e+00, 1.07e+00, 9.69e-01, 9.10e-01, 1.08e+00, 1.36e+00, 8.00e-01};
                double cutfmishits_tmp[9] = {5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01};
                double cuthoe_tmp[9] = {2.39e-02, 2.68e-02, 2.12e-02, 1.03e-01, 9.92e-03, 7.07e-02, 7.12e-02, 1.48e-01, 2.50e-01};
                double cuthoel_tmp[9] = {2.87e-02, 1.94e-02, 2.16e-02, 5.68e-02, 1.35e-02, 4.04e-02, 7.98e-02, 1.50e-01, 1.32e-01};
                double cutip_gsf_tmp[9] = {7.61e-03, 5.22e-03, 3.79e-02, 1.02e-02, 4.62e-02, 1.82e-02, 7.74e-02, 3.35e-02, 2.42e-02};
                double cutip_gsfl_tmp[9] = {7.81e-03, 4.25e-03, 3.08e-02, 1.04e-02, 2.35e-02, 2.45e-02, 7.90e-02, 2.88e-02, 3.40e-02};
                double cutiso_sum_tmp[9] = {5.40e+00, 5.41e+00, 5.88e+00, 4.32e+00, 3.86e+00, 4.33e+00, 5.86e+00, 9.05e+00, 1.78e+00};
                double cutiso_sumoet_tmp[9] = {2.53e+00, 2.10e+00, 1.87e+00, 1.84e+00, 1.79e+00, 1.61e+00, 2.53e+00, 1.98e+00, 2.76e+00};
                double cutiso_sumoetl_tmp[9] = {2.28e+00, 2.02e+00, 2.04e+00, 1.69e+00, 1.65e+00, 1.61e+00, 2.03e+00, 1.82e+00, 1.35e+00};
                double cutsee_tmp[9] = {9.99e-03, 9.61e-03, 9.65e-03, 2.75e-02, 2.61e-02, 2.64e-02, 9.18e-03, 2.44e-02, 2.58e-02};
                double cutseel_tmp[9] = {9.66e-03, 9.69e-03, 9.58e-03, 2.73e-02, 2.66e-02, 2.66e-02, 8.64e-03, 2.46e-02, 2.50e-02};

                eidAssign(cutdcotdist, cutdcotdist_tmp, 9);
                eidAssign(cutdetain, cutdetain_tmp, 9);
                eidAssign(cutdetainl, cutdetainl_tmp, 9);
                eidAssign(cutdphiin, cutdphiin_tmp, 9);
                eidAssign(cutdphiinl, cutdphiinl_tmp, 9);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 9);
                eidAssign(cutfmishits, cutfmishits_tmp, 9);
                eidAssign(cuthoe, cuthoe_tmp, 9);
                eidAssign(cuthoel, cuthoel_tmp, 9);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 9);
                eidAssign(cutip_gsfl, cutip_gsfl_tmp, 9);
                eidAssign(cutiso_sum, cutiso_sum_tmp, 9);
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 9);
                eidAssign(cutiso_sumoetl, cutiso_sumoetl_tmp, 9);
                eidAssign(cutsee, cutsee_tmp, 9);
                eidAssign(cutseel, cutseel_tmp, 9);

                return;
            }

        default:
            std::cout << "[eidGetCIC] ERROR! Invalid tightness level" << std::endl;

    }

    return;

}

void eidGetCIC_V04(const cic_tightness tightness, std::vector<double> &cutdcotdist, std::vector<double> &cutdetain, std::vector<double> &cutdphiin, std::vector<double> &cuteseedopcor, std::vector<double> &cutet, std::vector<double> &cutfmishits, std::vector<double> &cuthoe, std::vector<double> &cutip_gsf, std::vector<double> &cutiso_sum, std::vector<double> &cutiso_sumoet, std::vector<double> &cutsee)
{

    switch (tightness) {
        case CIC_VERYLOOSE:
            {
                double cutdcotdist_tmp[9] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutet_tmp[9] = {
                    0., 0., 0., 0., 0., 0., 0., 0., 0.
                };
                double cutip_gsf_tmp[9] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutiso_sum_tmp[9] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutiso_sumoet_tmp[9] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutdetain_tmp[9] = {
                    1.30e-02, 1.00e-02, 3.10e-02, 2.80e-02, 1.26e-02, 2.89e-02, 2.42e-02, 5.13e-02, 2.21e-02
                };
                double cutdphiin_tmp[9] = {
                    7.51e-02, 3.71e-01, 4.20e-01, 9.86e-02, 2.88e-01, 3.89e-01, 3.77e-01, 4.32e-01, 4.57e-01
                };
                double cuteseedopcor_tmp[9] = {
                    6.31e-01, 2.37e-01, 3.04e-01, 8.04e-01, 1.63e-01, 5.03e-01, 2.78e-01, 3.10e-01, 1.31e-01
                };
                double cutfmishits_tmp[9] = {
                    4.50e+00, 1.50e+00, 1.50e+00, 7.50e+00, 2.50e+00, 2.50e+00, 3.50e+00, 4.50e+00, 3.50e+00
                };
                double cuthoe_tmp[9] = {
                    2.47e-01, 1.11e-01, 1.49e-01, 3.82e-01, 7.17e-02, 1.47e-01, 1.16e+00, 5.04e+00, 3.33e+00
                };
                double cutsee_tmp[9] = {
                    1.92e-02, 1.98e-02, 2.53e-02, 5.28e-02, 3.91e-02, 4.61e-02, 2.66e-02, 6.58e-02, 3.20e+00
                };
                eidAssign(cutdcotdist, cutdcotdist_tmp, 9);
                eidAssign(cutdetain, cutdetain_tmp, 9);
                eidAssign(cutdphiin, cutdphiin_tmp, 9);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 9);
                eidAssign(cutet, cutet_tmp, 9);
                eidAssign(cutfmishits, cutfmishits_tmp, 9);
                eidAssign(cuthoe, cuthoe_tmp, 9);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 9);
                eidAssign(cutiso_sum, cutiso_sum_tmp, 9);
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 9);
                eidAssign(cutsee, cutsee_tmp, 9);
                return;
            }

        case CIC_LOOSE:
            {
                double cutdcotdist_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutet_tmp[18] = {
                    0., 0., 0., 0., 0., 0., 0., 0., 0.
                };
                double cutip_gsf_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutiso_sum_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutiso_sumoet_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutdetain_tmp[18] = {
                    1.30e-02, 5.95e-03, 3.10e-02, 1.68e-02, 8.44e-03, 1.70e-02, 1.55e-02, 5.13e-02, 1.61e-02
                };
                double cutdphiin_tmp[18] = {
                    7.51e-02, 3.30e-01, 4.20e-01, 9.86e-02, 2.84e-01, 3.28e-01, 3.77e-01, 4.32e-01, 3.74e-01
                };
                double cuteseedopcor_tmp[18] = {
                    6.31e-01, 3.02e-01, 3.04e-01, 8.10e-01, 2.23e-01, 5.03e-01, 2.78e-01, 3.10e-01, 4.69e-01
                };
                double cutfmishits_tmp[18] = {
                    4.50e+00, 1.50e+00, 1.50e+00, 2.50e+00, 2.50e+00, 1.50e+00, 2.50e+00, 4.50e+00, 5.00e-01
                };
                double cuthoe_tmp[18] = {
                    2.47e-01, 7.78e-02, 1.49e-01, 3.82e-01, 4.70e-02, 1.12e-01, 1.16e+00, 5.04e+00, 1.35e+00
                };
                double cutsee_tmp[18] = {
                    1.92e-02, 1.31e-02, 2.53e-02, 5.27e-02, 3.29e-02, 4.19e-02, 2.65e-02, 6.58e-02, 1.38e-01
                };
                eidAssign(cutdcotdist, cutdcotdist_tmp, 9);
                eidAssign(cutdetain, cutdetain_tmp, 9);
                eidAssign(cutdphiin, cutdphiin_tmp, 9);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 9);
                eidAssign(cutet, cutet_tmp, 9);
                eidAssign(cutfmishits, cutfmishits_tmp, 9);
                eidAssign(cuthoe, cuthoe_tmp, 9);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 9);
                eidAssign(cutiso_sum, cutiso_sum_tmp, 9);
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 9);
                eidAssign(cutsee, cutsee_tmp, 9);
                return;
            } 

        case CIC_MEDIUM:
            {
                double cutdcotdist_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutet_tmp[18] = {
                    0., 0., 0., 0., 0., 0., 0., 0., 0.
                };
                double cutip_gsf_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutiso_sum_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutiso_sumoet_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutdetain_tmp[18] = {
                    1.19e-02, 4.20e-03, 1.07e-02, 1.49e-02, 6.56e-03, 1.19e-02, 1.16e-02, 5.13e-02, 6.37e-03
                };
                double cutdphiin_tmp[18] = {
                    7.51e-02, 2.93e-01, 3.58e-01, 9.53e-02, 1.62e-01, 2.99e-01, 2.76e-01, 4.32e-01, 2.57e-01
                };
                double cuteseedopcor_tmp[18] = {
                    6.31e-01, 8.14e-01, 7.60e-01, 8.18e-01, 7.56e-01, 5.35e-01, 6.20e-01, 7.88e-01, 8.85e-01
                };
                double cutfmishits_tmp[18] = {
                    1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 2.50e+00, 1.50e+00, 5.00e-01
                };
                double cuthoe_tmp[18] = {
                    2.46e-01, 6.80e-02, 1.35e-01, 3.73e-01, 2.33e-02, 5.58e-02, 8.80e-01, 5.04e+00, 3.78e-02
                };
                double cutsee_tmp[18] = {
                    1.92e-02, 1.13e-02, 1.47e-02, 3.84e-02, 3.05e-02, 3.36e-02, 1.35e-02, 5.05e-02, 2.79e-02
                };
                eidAssign(cutdcotdist, cutdcotdist_tmp, 9);
                eidAssign(cutdetain, cutdetain_tmp, 9);
                eidAssign(cutdphiin, cutdphiin_tmp, 9);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 9);
                eidAssign(cutet, cutet_tmp, 9);
                eidAssign(cutfmishits, cutfmishits_tmp, 9);
                eidAssign(cuthoe, cuthoe_tmp, 9);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 9);
                eidAssign(cutiso_sum, cutiso_sum_tmp, 9);
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 9);
                eidAssign(cutsee, cutsee_tmp, 9);
                return;
            } 

        case CIC_TIGHT:
            {

                double cutdcotdist_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutet_tmp[18] = {
                    0., 0., 0., 0., 0., 0., 0., 0., 0.
                };
                double cutip_gsf_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutiso_sum_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutiso_sumoet_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutdetain_tmp[18] = {
                    9.28e-03, 3.56e-03, 7.16e-03, 1.31e-02, 5.81e-03, 9.79e-03, 1.15e-02, 1.66e-02, 3.19e-03
                };
                double cutdphiin_tmp[18] = {
                    4.66e-02, 7.80e-02, 2.64e-01, 4.42e-02, 3.20e-02, 2.37e-01, 8.25e-02, 2.07e-01, 5.39e-02
                };
                double cuteseedopcor_tmp[18] = {
                    6.48e-01, 8.97e-01, 8.91e-01, 8.39e-01, 8.35e-01, 6.49e-01, 6.76e-01, 8.70e-01, 9.91e-01
                };
                double cutfmishits_tmp[18] = {
                    1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 5.00e-01, 2.50e+00, 5.00e-01, 5.00e-01
                };
                double cuthoe_tmp[18] = {
                    9.94e-02, 5.61e-02, 1.05e-01, 9.73e-02, 1.81e-02, 3.06e-02, 5.57e-01, 5.04e+00, 1.06e-03
                };
                double cutsee_tmp[18] = {
                    1.56e-02, 1.07e-02, 1.23e-02, 3.35e-02, 2.98e-02, 3.06e-02, 1.07e-02, 3.79e-02, 1.01e-02
                }; 
                eidAssign(cutdcotdist, cutdcotdist_tmp, 9);
                eidAssign(cutdetain, cutdetain_tmp, 9);
                eidAssign(cutdphiin, cutdphiin_tmp, 9);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 9);
                eidAssign(cutet, cutet_tmp, 9);
                eidAssign(cutfmishits, cutfmishits_tmp, 9);
                eidAssign(cuthoe, cuthoe_tmp, 9);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 9);
                eidAssign(cutiso_sum, cutiso_sum_tmp, 9);                
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 9);
                eidAssign(cutsee, cutsee_tmp, 9); 
                return; 
            } 

        case CIC_SUPERTIGHT:
            {
                double cutdcotdist_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutet_tmp[18] = {
                    0., 0., 0., 0., 0., 0., 0., 0., 0.
                };
                double cutip_gsf_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutiso_sum_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutiso_sumoet_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutdetain_tmp[18] = {
                    9.28e-03, 3.41e-03, 5.60e-03, 9.00e-03, 5.15e-03, 8.03e-03, 1.06e-02, 1.51e-02, 3.15e-03
                };
                double cutdphiin_tmp[18] = {
                    3.23e-02, 3.51e-02, 1.61e-01, 3.06e-02, 2.07e-02, 6.15e-02, 5.82e-02, 6.05e-02, 3.66e-02
                };
                double cuteseedopcor_tmp[18] = {
                    7.35e-01, 9.41e-01, 9.53e-01, 8.86e-01, 8.85e-01, 9.38e-01, 7.98e-01, 9.26e-01, 1.02e+00
                };
                double cutfmishits_tmp[18] = {
                    1.50e+00, 1.50e+00, 5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01, 2.50e+00, 5.00e-01, 5.00e-01
                };
                double cuthoe_tmp[18] = {
                    5.25e-02, 4.62e-02, 4.98e-02, 5.05e-02, 1.60e-02, 2.04e-02, 2.18e-01, 5.02e+00, 2.96e-05
                };
                double cutsee_tmp[18] = {
                    1.22e-02, 1.04e-02, 1.12e-02, 3.01e-02, 2.82e-02, 2.88e-02, 9.95e-03, 2.74e-02, 9.07e-03
                };         
                eidAssign(cutdcotdist, cutdcotdist_tmp, 9);
                eidAssign(cutdetain, cutdetain_tmp, 9);
                eidAssign(cutdphiin, cutdphiin_tmp, 9);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 9);
                eidAssign(cutet, cutet_tmp, 9);
                eidAssign(cutfmishits, cutfmishits_tmp, 9);
                eidAssign(cuthoe, cuthoe_tmp, 9);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 9); 
                eidAssign(cutiso_sum, cutiso_sum_tmp, 9);                
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 9);
                eidAssign(cutsee, cutsee_tmp, 9); 
                return;  
            } 

        case CIC_HYPERTIGHT1:
            {
                double cutdcotdist_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutet_tmp[18] = {
                    0., 0., 0., 0., 0., 0., 0., 0., 0.
                };
                double cutip_gsf_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutiso_sum_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutiso_sumoet_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutdetain_tmp[18] = {
                    6.70e-03, 3.18e-03, 4.15e-03, 8.05e-03, 3.95e-03, 6.92e-03, 1.06e-02, 1.46e-02, 1.15e-03
                };
                double cutdphiin_tmp[18] = {
                    2.07e-02, 2.22e-02, 9.49e-02, 2.43e-02, 1.58e-02, 2.70e-02, 3.84e-02, 3.84e-02, 3.47e-02
                };
                double cuteseedopcor_tmp[18] = {
                    1.04e+00, 9.63e-01, 9.91e-01, 9.69e-01, 8.98e-01, 9.55e-01, 8.16e-01, 9.76e-01, 1.04e+00
                };
                double cutfmishits_tmp[18] = {
                    1.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01
                };
                double cuthoe_tmp[18] = {
                    3.83e-02, 3.84e-02, 3.50e-02, 2.33e-02, 1.26e-02, 1.64e-02, 1.12e-01, 2.72e+00, 2.30e-06
                };
                double cutsee_tmp[18] = {
                    1.05e-02, 1.01e-02, 1.02e-02, 2.86e-02, 2.72e-02, 2.75e-02, 9.67e-03, 2.54e-02, 8.85e-03
                };
                eidAssign(cutdcotdist, cutdcotdist_tmp, 9);
                eidAssign(cutdetain, cutdetain_tmp, 9);
                eidAssign(cutdphiin, cutdphiin_tmp, 9);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 9);
                eidAssign(cutet, cutet_tmp, 9);
                eidAssign(cutfmishits, cutfmishits_tmp, 9);
                eidAssign(cuthoe, cuthoe_tmp, 9);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 9);  
                eidAssign(cutiso_sum, cutiso_sum_tmp, 9);                
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 9);
                eidAssign(cutsee, cutsee_tmp, 9);          
                return; 
            } 

        case CIC_HYPERTIGHT2:
            {
                double cutdcotdist_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutet_tmp[18] = {
                    0., 0., 0., 0., 0., 0., 0., 0., 0.
                };
                double cutip_gsf_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutiso_sum_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutiso_sumoet_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutdetain_tmp[18] = {
                    3.12e-03, 2.09e-03, 3.22e-03, 4.93e-03, 2.99e-03, 5.23e-03, 8.85e-03, 1.21e-02, 1.15e-03
                };
                double cutdphiin_tmp[18] = {
                    1.45e-02, 7.75e-03, 3.30e-02, 1.17e-02, 7.93e-03, 1.64e-02, 2.10e-02, 2.40e-02, 3.47e-02
                };
                double cuteseedopcor_tmp[18] = {
                    1.09e+00, 9.90e-01, 1.07e+00, 9.71e-01, 9.60e-01, 9.55e-01, 9.87e-01, 1.04e+00, 1.04e+00
                };
                double cutfmishits_tmp[18] = {
                    5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01
                };
                double cuthoe_tmp[18] = {
                    3.25e-02, 3.29e-02, 2.96e-02, 1.52e-02, 1.22e-02, 1.31e-02, 6.02e-02, 1.22e-01, 2.30e-06
                };
                double cutsee_tmp[18] = {
                    9.85e-03, 9.79e-03, 9.64e-03, 2.72e-02, 2.64e-02, 2.65e-02, 9.37e-03, 2.37e-02, 8.85e-03
                };
                eidAssign(cutdcotdist, cutdcotdist_tmp, 9);
                eidAssign(cutdetain, cutdetain_tmp, 9);
                eidAssign(cutdphiin, cutdphiin_tmp, 9);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 9);
                eidAssign(cutet, cutet_tmp, 9);
                eidAssign(cutfmishits, cutfmishits_tmp, 9);
                eidAssign(cuthoe, cuthoe_tmp, 9);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 9); 
                eidAssign(cutiso_sum, cutiso_sum_tmp, 9);                
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 9);
                eidAssign(cutsee, cutsee_tmp, 9);   
                return; 
            } 

        case CIC_HYPERTIGHT3:
            {
                double cutdcotdist_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutet_tmp[18] = {
                    0., 0., 0., 0., 0., 0., 0., 0., 0.
                };
                double cutip_gsf_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutiso_sum_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutiso_sumoet_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutdetain_tmp[18] = {
                    1.73e-03, 1.23e-03, 2.47e-03, 4.18e-03, 2.10e-03, 3.28e-03, 8.85e-03, 7.63e-03, 1.15e-03
                };
                double cutdphiin_tmp[18] = {
                    5.32e-03, 2.99e-03, 1.72e-02, 6.27e-03, 7.93e-03, 9.69e-03, 7.60e-03, 1.67e-02, 3.47e-02
                };
                double cuteseedopcor_tmp[18] = {
                    1.14e+00, 1.00e+00, 1.14e+00, 1.02e+00, 9.64e-01, 1.15e+00, 1.05e+00, 1.22e+00, 1.04e+00
                };
                double cutfmishits_tmp[18] = {
                    5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01
                };
                double cuthoe_tmp[18] = {
                    2.93e-02, 2.80e-02, 2.81e-02, 1.35e-02, 9.55e-03, 1.21e-02, 3.27e-02, 3.45e-02, 2.30e-06
                };
                double cutsee_tmp[18] = {
                    9.41e-03, 9.54e-03, 9.34e-03, 2.62e-02, 2.56e-02, 2.56e-02, 9.08e-03, 2.35e-02, 8.85e-03
                };
                eidAssign(cutdcotdist, cutdcotdist_tmp, 9);
                eidAssign(cutdetain, cutdetain_tmp, 9);
                eidAssign(cutdphiin, cutdphiin_tmp, 9);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 9);
                eidAssign(cutet, cutet_tmp, 9);
                eidAssign(cutfmishits, cutfmishits_tmp, 9);
                eidAssign(cuthoe, cuthoe_tmp, 9);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 9);  
                eidAssign(cutiso_sum, cutiso_sum_tmp, 9);                
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 9);
                eidAssign(cutsee, cutsee_tmp, 9);                
                return; 
            } 

        case CIC_HYPERTIGHT4:
            {
                double cutdcotdist_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutet_tmp[18] = {
                    0., 0., 0., 0., 0., 0., 0., 0., 0.
                };
                double cutip_gsf_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutiso_sum_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutiso_sumoet_tmp[18] = {
                    9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
                };
                double cutdetain_tmp[18] = {
                    1.47e-03, 9.67e-04, 1.87e-03, 2.27e-03, 2.10e-03, 1.77e-03, 5.30e-03, 3.60e-03, 1.15e-03
                };
                double cutdphiin_tmp[18] = {
                    5.03e-03, 2.96e-03, 9.84e-03, 5.85e-03, 7.93e-03, 9.08e-03, 7.10e-03, 1.24e-02, 3.47e-02
                };
                double cuteseedopcor_tmp[18] = {
                    1.15e+00, 1.01e+00, 1.15e+00, 1.13e+00, 9.64e-01, 1.30e+00, 1.07e+00, 1.29e+00, 1.04e+00
                };
                double cutfmishits_tmp[18] = {
                    -5.00e-01, -5.00e-01, 5.00e-01, -5.00e-01, -5.00e-01, -5.00e-01, -5.00e-01, 5.00e-01, -5.00e-01
                };
                double cuthoe_tmp[18] = {
                    2.28e-03, 2.18e-03, 2.64e-02, 1.05e-02, 9.55e-03, 7.95e-03, 3.27e-02, 2.63e-02, 2.30e-06
                };
                double cutsee_tmp[18] = {
                    9.15e-03, 9.35e-03, 9.16e-03, 2.54e-02, 2.56e-02, 2.49e-02, 8.94e-03, 2.33e-02, 8.85e-03
                };
                eidAssign(cutdcotdist, cutdcotdist_tmp, 9);
                eidAssign(cutdetain, cutdetain_tmp, 9);
                eidAssign(cutdphiin, cutdphiin_tmp, 9);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 9);
                eidAssign(cutet, cutet_tmp, 9);
                eidAssign(cutfmishits, cutfmishits_tmp, 9);
                eidAssign(cuthoe, cuthoe_tmp, 9);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 9);   
                eidAssign(cutiso_sum, cutiso_sum_tmp, 9);                
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 9);
                eidAssign(cutsee, cutsee_tmp, 9);           
                return;  
            } 

        default:    
            std::cout << "[eidGetCIC] ERROR! Invalid tightness level" << std::endl;

    }           

}


void eidGetCIC_V03(const cic_tightness tightness, std::vector<double> &cutdcotdist, std::vector<double> &cutdetain, std::vector<double> &cutdphiin, std::vector<double> &cuteseedopcor, std::vector<double> &cutet, std::vector<double> &cutfmishits, std::vector<double> &cuthoe, std::vector<double> &cutip_gsf, std::vector<double> &cutiso_sum, std::vector<double> &cutiso_sumoet, std::vector<double> &cutsee)
{  
    switch (tightness) {
        case CIC_VERYLOOSE:
            {  
                double cutdcotdist_tmp[27] = {
                    3.95e-02, 3.94e-02, 3.98e-02, 3.96e-02, 3.94e-02, 3.92e-02, 3.80e-02, 3.90e-02, 3.98e-02,
                    3.93e-02, 3.92e-02, 3.84e-02, 3.94e-02, 3.89e-02, 3.90e-02, 3.96e-02, 3.94e-02, 2.79e-02,
                    3.85e-02, 3.00e-02, 1.02e-02, 3.93e-02, 3.48e-02, 3.88e-02, 1.51e-02, 3.67e-02, 3.99e-02
                };
                double cutdetain_tmp[27] = {
                    1.08e-02, 6.70e-03, 2.72e-02, 1.59e-02, 1.19e-02, 1.94e-02, 1.52e-02, 4.77e-02, 2.75e-02,
                    1.21e-02, 4.20e-03, 1.26e-02, 1.38e-02, 7.97e-03, 1.41e-02, 1.26e-02, 2.51e-02, 2.34e-02,
                    1.38e-02, 3.89e-03, 1.24e-02, 1.97e-02, 8.08e-03, 1.55e-02, 1.27e-02, 2.14e-02, 8.62e-03
                };
                double cutdphiin_tmp[27] = {
                    4.10e-02, 2.78e-01, 3.95e-01, 4.70e-02, 2.74e-01, 3.28e-01, 3.44e-01, 4.65e-01, 6.57e-01,
                    6.10e-02, 2.71e-01, 3.43e-01, 7.31e-02, 2.27e-01, 2.79e-01, 2.85e-01, 3.48e-01, 2.71e-01,
                    9.35e-02, 7.33e-02, 3.58e-01, 1.19e-01, 1.96e-01, 3.19e-01, 1.21e-01, 3.39e-01, 2.86e-01
                };
                double cuteseedopcor_tmp[27] = {
                    7.79e-01, 2.64e-01, 4.83e-01, 9.04e-01, 1.65e-01, 6.36e-01, 1.08e-01, 2.84e-01, 7.61e-02,
                    5.81e-01, 1.97e-01, 4.74e-01, 8.02e-01, 1.22e-01, 6.07e-01, 3.98e-01, 7.05e-01, 5.89e-01,
                    4.91e-01, 8.61e-01, 2.66e-01, 7.95e-01, 8.17e-01, 5.07e-01, 2.73e-01, 7.56e-01, 2.46e-01
                };
                double cutet_tmp[27] = {
                    -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
                    -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
                    1.20e+01, 1.20e+01, 1.20e+01, 1.20e+01, 1.20e+01, 1.20e+01, 1.20e+01, 1.20e+01, 1.20e+01
                };
                double cutfmishits_tmp[27] = {
                    4.50e+00, 1.50e+00, 1.50e+00, 5.50e+00, 2.50e+00, 2.50e+00, 3.50e+00, 5.50e+00, 9.50e+00,
                    2.50e+00, 1.50e+00, 1.50e+00, 2.50e+00, 1.50e+00, 1.50e+00, 2.50e+00, 2.50e+00, 1.50e+00,
                    2.50e+00, 1.50e+00, 2.50e+00, 1.50e+00, 1.50e+00, 5.00e-01, 2.50e+00, 2.50e+00, 5.00e-01
                };
                double cuthoe_tmp[27] = {
                    1.68e-01, 9.29e-02, 1.45e-01, 3.70e-01, 8.32e-02, 1.47e-01, 4.05e-01, 2.68e+00, 2.99e+00,
                    2.52e-01, 8.32e-02, 1.49e-01, 3.69e-01, 4.34e-02, 1.47e-01, 6.02e-01, 2.77e+00, 1.20e+00,
                    2.60e-01, 7.18e-02, 1.18e-01, 3.94e-01, 2.81e-02, 5.31e-02, 1.14e+00, 4.65e+00, 2.86e+00
                };
                double cutip_gsf_tmp[27] = {
                    4.31e-02, 8.76e-02, 1.39e-01, 1.71e-01, 1.98e-01, 3.58e-01, 9.43e-01, 2.29e-01, 1.24e-01,
                    3.25e-02, 8.48e-02, 1.69e-01, 9.04e-02, 1.47e-01, 1.40e-01, 6.47e-02, 1.97e-01, 2.53e-01,
                    1.40e-02, 1.44e-02, 6.93e-02, 3.48e-02, 1.31e-01, 1.34e-01, 8.59e-02, 7.91e-02, 8.23e-02
                };
                double cutiso_sum_tmp[27] = {
                    4.25e+01, 1.27e+01, 1.34e+01, 2.40e+01, 6.38e+00, 9.62e+00, 1.42e+01, 3.30e+01, 4.26e+00,
                    1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05,
                    1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05
                };
                double cutiso_sumoet_tmp[27] = {
                    4.54e+01, 2.00e+01, 1.54e+01, 3.12e+01, 9.96e+00, 1.55e+01, 1.70e+01, 1.12e+01, 5.12e+00,
                    3.98e+01, 2.15e+01, 1.69e+01, 2.77e+01, 1.04e+01, 1.16e+01, 2.37e+01, 2.00e+01, 7.58e+00,
                    3.20e+01, 3.02e+01, 2.38e+01, 2.16e+01, 1.32e+01, 1.80e+01, 3.04e+01, 3.10e+01, 1.65e+01
                };
                double cutsee_tmp[27] = {
                    1.90e-02, 1.74e-02, 2.53e-02, 4.71e-02, 3.30e-02, 3.75e-02, 1.60e-02, 4.16e-02, 3.19e+00,
                    1.88e-02, 1.12e-02, 1.55e-02, 3.91e-02, 3.25e-02, 3.71e-02, 1.29e-02, 5.49e-02, 1.01e+03,
                    1.74e-02, 1.19e-02, 1.58e-02, 4.99e-02, 3.48e-02, 4.48e-02, 1.20e-02, 7.56e-02, 2.96e-02
                };
                eidAssign(cutdcotdist, cutdcotdist_tmp, 27);
                eidAssign(cutdetain, cutdetain_tmp, 27);
                eidAssign(cutdphiin, cutdphiin_tmp, 27);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 27);
                eidAssign(cutet, cutet_tmp, 27);
                eidAssign(cutfmishits, cutfmishits_tmp, 27);
                eidAssign(cuthoe, cuthoe_tmp, 27);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 27);
                eidAssign(cutiso_sum, cutiso_sum_tmp, 27);
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 27);
                eidAssign(cutsee, cutsee_tmp, 27);
                return;
            }

        case CIC_LOOSE:
            {
                double cutdcotdist_tmp[27] = {
                    3.93e-02, 3.92e-02, 3.97e-02, 3.94e-02, 3.93e-02, 3.90e-02, 3.78e-02, 3.88e-02, 3.82e-02,
                    3.85e-02, 1.67e-02, 3.25e-03, 3.94e-02, 3.87e-02, 3.88e-02, 2.27e-02, 2.58e-02, 1.27e-02,
                    2.98e-02, 3.00e-02, 9.46e-03, 3.90e-02, 2.31e-02, 2.78e-02, 1.62e-03, 3.67e-02, 1.99e-02
                };
                double cutdetain_tmp[27] = {
                    9.89e-03, 4.84e-03, 1.46e-02, 1.46e-02, 9.02e-03, 1.72e-02, 1.37e-02, 4.77e-02, 2.75e-02,
                    9.67e-03, 3.77e-03, 9.24e-03, 1.30e-02, 6.66e-03, 1.23e-02, 1.25e-02, 2.28e-02, 1.12e-02,
                    1.06e-02, 3.80e-03, 8.97e-03, 1.39e-02, 6.67e-03, 1.22e-02, 1.22e-02, 1.93e-02, 2.39e-03
                };
                double cutdphiin_tmp[27] = {
                    4.10e-02, 2.75e-01, 3.65e-01, 4.70e-02, 2.73e-01, 2.96e-01, 3.29e-01, 4.65e-01, 6.27e-01,
                    5.81e-02, 9.54e-02, 3.27e-01, 7.02e-02, 5.82e-02, 2.79e-01, 1.17e-01, 3.18e-01, 2.46e-01,
                    8.21e-02, 5.20e-02, 2.92e-01, 1.16e-01, 4.35e-02, 3.12e-01, 1.18e-01, 2.96e-01, 4.59e-02
                };
                double cuteseedopcor_tmp[27] = {
                    7.80e-01, 3.02e-01, 4.83e-01, 9.04e-01, 1.68e-01, 6.45e-01, 1.08e-01, 2.84e-01, 3.24e-01,
                    5.91e-01, 2.86e-01, 4.88e-01, 8.13e-01, 7.91e-01, 6.72e-01, 3.98e-01, 8.34e-01, 8.78e-01,
                    5.15e-01, 9.37e-01, 8.06e-01, 8.16e-01, 8.50e-01, 5.07e-01, 3.67e-01, 8.30e-01, 6.48e-01
                };
                double cutet_tmp[27] = {
                    -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
                    -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
                    1.20e+01, 1.20e+01, 1.20e+01, 1.20e+01, 1.20e+01, 1.20e+01, 1.20e+01, 1.20e+01, 1.25e+01
                };
                double cutfmishits_tmp[27] = {
                    4.50e+00, 1.50e+00, 1.50e+00, 2.50e+00, 2.50e+00, 1.50e+00, 2.50e+00, 2.50e+00, 1.50e+00,
                    2.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 5.00e-01, 2.50e+00, 2.50e+00, 5.00e-01,
                    2.50e+00, 1.50e+00, 5.00e-01, 1.50e+00, 1.50e+00, 5.00e-01, 2.50e+00, 5.00e-01, 5.00e-01
                };
                double cuthoe_tmp[27] = {
                    1.66e-01, 7.71e-02, 1.44e-01, 3.70e-01, 4.97e-02, 1.39e-01, 4.01e-01, 2.68e+00, 5.16e-01,
                    2.34e-01, 5.56e-02, 1.44e-01, 3.68e-01, 3.10e-02, 1.20e-01, 6.02e-01, 2.01e+00, 1.05e+00,
                    1.04e-01, 6.30e-02, 5.65e-02, 3.80e-01, 1.92e-02, 2.94e-02, 5.37e-01, 4.65e+00, 1.87e+00
                };
                double cutip_gsf_tmp[27] = {
                    4.31e-02, 7.67e-02, 1.39e-01, 1.01e-01, 1.49e-01, 1.54e-01, 9.32e-01, 1.50e-01, 1.24e-01,
                    2.38e-02, 4.67e-02, 7.59e-02, 3.69e-02, 1.47e-01, 9.86e-02, 6.26e-02, 1.95e-01, 1.16e-01,
                    1.22e-02, 1.25e-02, 6.93e-02, 1.62e-02, 8.90e-02, 6.73e-02, 4.67e-02, 6.51e-02, 2.21e-02
                };
                double cutiso_sum_tmp[27] = {
                    3.15e+01, 1.03e+01, 8.80e+00, 1.10e+01, 6.13e+00, 6.94e+00, 7.52e+00, 9.00e+00, 3.50e+00,
                    1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05,
                    1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05
                };
                double cutiso_sumoet_tmp[27] = {
                    2.89e+01, 1.53e+01, 1.20e+01, 1.83e+01, 7.17e+00, 9.42e+00, 1.10e+01, 9.81e+00, 3.94e+00,
                    2.27e+01, 1.59e+01, 1.23e+01, 1.70e+01, 7.58e+00, 8.89e+00, 1.52e+01, 1.27e+01, 6.17e+00,
                    2.08e+01, 2.12e+01, 1.72e+01, 1.55e+01, 9.37e+00, 1.06e+01, 1.98e+01, 2.21e+01, 1.56e+01
                };
                double cutsee_tmp[27] = {
                    1.75e-02, 1.27e-02, 1.77e-02, 3.73e-02, 3.14e-02, 3.29e-02, 1.57e-02, 4.09e-02, 1.40e-01,
                    1.69e-02, 1.06e-02, 1.42e-02, 3.63e-02, 3.22e-02, 3.54e-02, 1.17e-02, 3.72e-02, 2.82e+01,
                    1.71e-02, 1.13e-02, 1.40e-02, 4.03e-02, 3.23e-02, 4.11e-02, 1.04e-02, 4.36e-02, 1.14e-02
                };
                eidAssign(cutdcotdist, cutdcotdist_tmp, 27);
                eidAssign(cutdetain, cutdetain_tmp, 27);
                eidAssign(cutdphiin, cutdphiin_tmp, 27);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 27);
                eidAssign(cutet, cutet_tmp, 27);
                eidAssign(cutfmishits, cutfmishits_tmp, 27);
                eidAssign(cuthoe, cuthoe_tmp, 27);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 27);
                eidAssign(cutiso_sum, cutiso_sum_tmp, 27);
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 27);
                eidAssign(cutsee, cutsee_tmp, 27);
                return;
            }

        case CIC_MEDIUM:
            {
                double cutdcotdist_tmp[27] = {
                    3.93e-02, 2.91e-02, 3.23e-02, 3.94e-02, 3.90e-02, 3.90e-02, 3.75e-02, 3.88e-02, 3.82e-02,
                    3.18e-02, 4.69e-04, 9.09e-05, 3.81e-02, 2.32e-02, 1.12e-02, 1.78e-02, 1.93e-02, 9.75e-03,
                    8.35e-04, 2.63e-02, 4.58e-04, 3.26e-02, 5.50e-03, 1.82e-02, 1.62e-03, 2.67e-02, 1.55e-03
                };
                double cutdetain_tmp[27] = {
                    9.50e-03, 4.01e-03, 8.81e-03, 1.31e-02, 6.99e-03, 1.09e-02, 1.08e-02, 4.58e-02, 1.68e-02,
                    8.32e-03, 3.53e-03, 7.05e-03, 1.29e-02, 5.47e-03, 9.19e-03, 1.22e-02, 1.69e-02, 2.24e-03,
                    1.06e-02, 3.80e-03, 6.88e-03, 1.08e-02, 6.13e-03, 8.99e-03, 1.19e-02, 1.69e-02, 1.44e-03
                };
                double cutdphiin_tmp[27] = {
                    4.04e-02, 9.19e-02, 3.23e-01, 4.70e-02, 2.50e-01, 2.71e-01, 2.77e-01, 4.52e-01, 6.27e-01,
                    5.57e-02, 5.21e-02, 3.25e-01, 6.23e-02, 2.52e-02, 1.23e-01, 6.29e-02, 2.53e-01, 1.46e-02,
                    7.62e-02, 2.44e-02, 1.71e-01, 5.60e-02, 2.27e-02, 6.54e-02, 7.80e-02, 8.19e-02, 3.46e-02
                };
                double cuteseedopcor_tmp[27] = {
                    7.84e-01, 3.48e-01, 5.70e-01, 9.04e-01, 1.68e-01, 6.45e-01, 2.00e-01, 2.92e-01, 8.94e-01,
                    8.12e-01, 9.21e-01, 8.64e-01, 8.68e-01, 8.73e-01, 9.55e-01, 4.88e-01, 9.21e-01, 9.62e-01,
                    5.15e-01, 9.53e-01, 9.52e-01, 8.23e-01, 8.79e-01, 9.74e-01, 8.03e-01, 9.23e-01, 1.01e+00
                };
                double cutet_tmp[27] = {
                    -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
                    -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
                    1.20e+01, 1.20e+01, 1.20e+01, 1.22e+01, 1.21e+01, 1.22e+01, 1.20e+01, 1.22e+01, 1.77e+01
                };
                double cutfmishits_tmp[27] = {
                    2.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 2.50e+00, 2.50e+00, 1.50e+00,
                    2.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 2.50e+00, 5.00e-01, 5.00e-01,
                    2.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01
                };
                double cuthoe_tmp[27] = {
                    9.70e-02, 5.25e-02, 1.12e-01, 3.65e-01, 3.72e-02, 7.75e-02, 2.67e-01, 2.65e+00, 4.15e-01,
                    1.08e-01, 5.03e-02, 1.41e-01, 3.63e-01, 2.37e-02, 6.00e-02, 6.02e-01, 1.30e+00, 6.53e-01,
                    5.99e-02, 5.81e-02, 3.10e-02, 1.06e-01, 1.45e-02, 2.28e-02, 1.24e-01, 4.65e+00, 6.82e-01
                };
                double cutip_gsf_tmp[27] = {
                    4.09e-02, 5.19e-02, 1.16e-01, 4.47e-02, 1.28e-01, 1.26e-01, 2.42e-01, 1.19e-01, 1.24e-01,
                    1.40e-02, 2.04e-02, 5.64e-02, 3.04e-02, 8.69e-02, 4.84e-02, 5.71e-02, 1.89e-01, 2.04e-02,
                    1.02e-02, 9.96e-03, 4.00e-02, 1.15e-02, 1.90e-02, 3.41e-02, 1.88e-02, 2.35e-02, 4.20e-03
                };
                double cutiso_sum_tmp[27] = {
                    1.80e+01, 8.66e+00, 6.93e+00, 8.76e+00, 4.11e+00, 4.94e+00, 5.38e+00, 5.71e+00, 1.88e+00,
                    1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05,
                    1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05
                };
                double cutiso_sumoet_tmp[27] = {
                    1.76e+01, 1.16e+01, 8.63e+00, 1.42e+01, 4.70e+00, 5.97e+00, 7.31e+00, 8.88e+00, 2.65e+00,
                    1.48e+01, 1.29e+01, 9.81e+00, 1.11e+01, 6.59e+00, 6.67e+00, 1.16e+01, 9.95e+00, 6.04e+00,
                    1.39e+01, 1.46e+01, 1.26e+01, 1.01e+01, 5.91e+00, 7.28e+00, 1.56e+01, 1.62e+01, 1.26e+01
                };
                double cutsee_tmp[27] = {
                    1.51e-02, 1.12e-02, 1.40e-02, 3.31e-02, 3.08e-02, 3.17e-02, 1.27e-02, 3.37e-02, 5.23e-02,
                    1.47e-02, 1.04e-02, 1.29e-02, 3.51e-02, 3.02e-02, 3.18e-02, 1.14e-02, 3.60e-02, 8.02e-01,
                    1.58e-02, 1.06e-02, 1.27e-02, 3.77e-02, 3.15e-02, 3.68e-02, 9.68e-03, 3.39e-02, 1.01e-02
                };
                eidAssign(cutdcotdist, cutdcotdist_tmp, 27);
                eidAssign(cutdetain, cutdetain_tmp, 27);
                eidAssign(cutdphiin, cutdphiin_tmp, 27);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 27);
                eidAssign(cutet, cutet_tmp, 27);
                eidAssign(cutfmishits, cutfmishits_tmp, 27);
                eidAssign(cuthoe, cuthoe_tmp, 27);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 27);
                eidAssign(cutiso_sum, cutiso_sum_tmp, 27);
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 27);
                eidAssign(cutsee, cutsee_tmp, 27);
                return;
            }

        case CIC_TIGHT:
            {
                double cutdcotdist_tmp[27] = {
                    3.93e-02, 2.56e-02, 6.91e-03, 3.94e-02, 3.86e-02, 3.90e-02, 3.25e-02, 3.84e-02, 3.82e-02,
                    2.45e-02, 2.81e-04, 5.46e-05, 3.42e-02, 2.32e-02, 1.07e-03, 1.78e-02, 1.93e-02, 7.58e-04,
                    1.08e-04, 2.48e-02, 4.58e-04, 1.29e-02, 1.19e-03, 1.82e-02, 4.53e-05, 1.89e-02, 9.28e-04
                };
                double cutdetain_tmp[27] = {
                    8.11e-03, 3.41e-03, 6.33e-03, 1.03e-02, 6.67e-03, 1.00e-02, 1.06e-02, 1.45e-02, 1.63e-02,
                    7.60e-03, 2.59e-03, 5.11e-03, 9.41e-03, 4.30e-03, 8.57e-03, 1.20e-02, 1.69e-02, 1.72e-03,
                    8.61e-03, 3.62e-03, 6.01e-03, 9.25e-03, 4.89e-03, 8.32e-03, 1.19e-02, 1.69e-02, 9.96e-04
                };
                double cutdphiin_tmp[27] = {
                    4.04e-02, 4.99e-02, 2.63e-01, 4.20e-02, 4.84e-02, 2.41e-01, 2.42e-01, 2.31e-01, 2.86e-01,
                    5.52e-02, 3.38e-02, 1.54e-01, 6.23e-02, 1.83e-02, 3.92e-02, 5.47e-02, 5.88e-02, 6.54e-03,
                    4.20e-02, 2.17e-02, 8.85e-02, 4.45e-02, 1.41e-02, 2.34e-02, 6.50e-02, 2.58e-02, 3.46e-02
                };
                double cuteseedopcor_tmp[27] = {
                    7.84e-01, 3.66e-01, 5.70e-01, 9.11e-01, 2.98e-01, 6.45e-01, 5.10e-01, 4.97e-01, 9.32e-01,
                    8.35e-01, 9.68e-01, 9.69e-01, 9.23e-01, 8.98e-01, 9.80e-01, 6.30e-01, 9.71e-01, 1.00e+00,
                    5.15e-01, 9.63e-01, 9.86e-01, 8.23e-01, 8.79e-01, 1.01e+00, 9.31e-01, 9.37e-01, 1.05e+00
                };
                double cutet_tmp[27] = {
                    -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
                    -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
                    1.37e+01, 1.32e+01, 1.36e+01, 1.42e+01, 1.41e+01, 1.39e+01, 1.29e+01, 1.49e+01, 1.77e+01
                };
                double cutfmishits_tmp[27] = {
                    2.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 5.00e-01, 2.50e+00, 5.00e-01, 5.00e-01,
                    2.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01,
                    2.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01
                };
                double cuthoe_tmp[27] = {
                    7.83e-02, 3.87e-02, 1.05e-01, 1.18e-01, 2.27e-02, 6.20e-02, 1.30e-01, 2.47e+00, 3.80e-01,
                    8.88e-02, 5.03e-02, 9.55e-02, 7.41e-02, 1.50e-02, 3.00e-02, 5.89e-01, 1.13e+00, 6.12e-01,
                    4.94e-02, 4.61e-02, 2.92e-02, 3.69e-02, 1.13e-02, 1.45e-02, 1.24e-01, 2.05e+00, 6.10e-01
                };
                double cutip_gsf_tmp[27] = {
                    2.13e-02, 4.22e-02, 6.32e-02, 3.61e-02, 7.30e-02, 1.26e-01, 1.71e-01, 1.19e-01, 3.72e-02,
                    1.31e-02, 1.46e-02, 5.64e-02, 1.52e-02, 2.22e-02, 2.68e-02, 3.14e-02, 8.84e-02, 3.74e-03,
                    8.52e-03, 7.61e-03, 1.43e-02, 1.06e-02, 1.27e-02, 1.19e-02, 1.23e-02, 2.35e-02, 3.63e-03
                };
                double cutiso_sum_tmp[27] = {
                    1.18e+01, 8.31e+00, 6.26e+00, 6.18e+00, 3.28e+00, 4.38e+00, 4.17e+00, 5.40e+00, 1.57e+00,
                    1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05,
                    1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05
                };
                double cutiso_sumoet_tmp[27] = {
                    1.37e+01, 1.16e+01, 7.14e+00, 9.98e+00, 3.52e+00, 4.87e+00, 6.24e+00, 7.96e+00, 2.53e+00,
                    1.12e+01, 1.19e+01, 7.88e+00, 8.16e+00, 5.58e+00, 5.03e+00, 1.14e+01, 8.15e+00, 5.79e+00,
                    1.04e+01, 1.11e+01, 1.04e+01, 7.47e+00, 5.08e+00, 5.90e+00, 1.18e+01, 1.41e+01, 1.17e+01
                };
                double cutsee_tmp[27] = {
                    1.43e-02, 1.05e-02, 1.23e-02, 3.24e-02, 3.07e-02, 3.01e-02, 1.09e-02, 2.70e-02, 2.92e-02,
                    1.33e-02, 1.04e-02, 1.16e-02, 3.32e-02, 2.96e-02, 3.10e-02, 9.81e-03, 3.07e-02, 7.20e-02,
                    1.49e-02, 1.05e-02, 1.10e-02, 3.42e-02, 3.07e-02, 3.03e-02, 9.54e-03, 2.65e-02, 1.01e-02
                };
                eidAssign(cutdcotdist, cutdcotdist_tmp, 27);
                eidAssign(cutdetain, cutdetain_tmp, 27);
                eidAssign(cutdphiin, cutdphiin_tmp, 27);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 27);
                eidAssign(cutet, cutet_tmp, 27);
                eidAssign(cutfmishits, cutfmishits_tmp, 27);
                eidAssign(cuthoe, cuthoe_tmp, 27);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 27);
                eidAssign(cutiso_sum, cutiso_sum_tmp, 27);
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 27);
                eidAssign(cutsee, cutsee_tmp, 27);
                return;
            }

        case CIC_SUPERTIGHT:
            {
                double cutdcotdist_tmp[27] = {
                    3.44e-02, 2.37e-02, 6.91e-03, 3.93e-02, 3.64e-02, 2.44e-02, 2.85e-02, 3.84e-02, 3.82e-02,
                    2.00e-02, 2.81e-04, 1.53e-06, 2.19e-02, 2.20e-02, 1.07e-03, 1.72e-02, 1.93e-02, 7.58e-04,
                    2.34e-05, 2.29e-02, 4.58e-04, 1.29e-02, 3.33e-05, 1.82e-02, 1.27e-06, 4.07e-03, 9.28e-04
                };
                double cutdetain_tmp[27] = {
                    8.11e-03, 3.41e-03, 5.58e-03, 8.58e-03, 6.18e-03, 8.69e-03, 1.04e-02, 1.21e-02, 1.63e-02,
                    6.43e-03, 2.40e-03, 4.08e-03, 7.25e-03, 3.28e-03, 6.70e-03, 1.20e-02, 1.47e-02, 1.72e-03,
                    8.61e-03, 2.19e-03, 1.62e-03, 9.05e-03, 4.67e-03, 3.40e-03, 1.19e-02, 4.73e-03, 6.40e-04
                };
                double cutdphiin_tmp[27] = {
                    3.78e-02, 3.63e-02, 2.43e-01, 3.40e-02, 2.16e-02, 1.41e-01, 1.20e-01, 1.65e-01, 1.08e-01,
                    4.91e-02, 2.30e-02, 1.18e-01, 4.77e-02, 9.59e-03, 2.26e-02, 4.63e-02, 2.87e-02, 6.54e-03,
                    4.01e-02, 7.23e-03, 8.85e-02, 3.19e-02, 1.07e-02, 1.44e-02, 6.50e-02, 2.04e-02, 3.46e-02
                };
                double cuteseedopcor_tmp[27] = {
                    8.17e-01, 8.86e-01, 8.88e-01, 9.22e-01, 3.16e-01, 7.12e-01, 6.88e-01, 8.72e-01, 9.60e-01,
                    8.86e-01, 9.78e-01, 9.96e-01, 9.55e-01, 9.49e-01, 9.80e-01, 8.48e-01, 1.03e+00, 1.00e+00,
                    6.86e-01, 9.68e-01, 9.86e-01, 9.28e-01, 9.02e-01, 1.19e+00, 1.01e+00, 9.68e-01, 1.07e+00
                };
                double cutet_tmp[27] = {
                    -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
                    -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
                    1.58e+01, 1.52e+01, 1.73e+01, 1.61e+01, 1.56e+01, 1.70e+01, 1.42e+01, 1.83e+01, 1.86e+01
                };
                double cutfmishits_tmp[27] = {
                    2.50e+00, 1.50e+00, 5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01,
                    5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01,
                    5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01
                };
                double cuthoe_tmp[27] = {
                    7.07e-02, 3.00e-02, 4.77e-02, 1.18e-01, 2.27e-02, 4.25e-02, 1.30e-01, 5.90e-01, 1.66e-01,
                    7.18e-02, 4.83e-02, 8.90e-02, 4.99e-02, 1.24e-02, 2.46e-02, 5.55e-01, 7.76e-01, 6.12e-01,
                    2.42e-02, 3.01e-02, 1.07e-02, 1.73e-02, 1.09e-02, 1.45e-02, 1.12e-01, 7.45e-01, 6.10e-01
                };
                double cutip_gsf_tmp[27] = {
                    1.61e-02, 2.78e-02, 6.32e-02, 2.85e-02, 7.30e-02, 7.33e-02, 7.21e-02, 4.57e-02, 2.49e-02,
                    1.19e-02, 1.04e-02, 4.10e-02, 1.31e-02, 2.14e-02, 1.46e-02, 2.07e-02, 5.33e-02, 3.74e-03,
                    8.31e-03, 6.58e-03, 8.53e-03, 8.37e-03, 9.24e-03, 7.92e-03, 1.11e-02, 7.78e-03, 2.61e-03
                };
                double cutiso_sum_tmp[27] = {
                    8.92e+00, 6.66e+00, 5.59e+00, 4.74e+00, 3.21e+00, 3.25e+00, 4.06e+00, 3.79e+00, 1.17e+00,
                    1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05,
                    1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05
                };
                double cutiso_sumoet_tmp[27] = {
                    9.52e+00, 9.10e+00, 7.14e+00, 6.60e+00, 3.33e+00, 3.53e+00, 3.33e+00, 5.54e+00, 1.78e+00,
                    8.78e+00, 8.13e+00, 6.90e+00, 6.14e+00, 3.58e+00, 4.11e+00, 1.00e+01, 6.19e+00, 5.79e+00,
                    9.13e+00, 9.62e+00, 8.46e+00, 6.22e+00, 4.80e+00, 5.70e+00, 1.06e+01, 1.31e+01, 6.75e+00
                };
                double cutsee_tmp[27] = {
                    1.36e-02, 1.04e-02, 1.16e-02, 3.21e-02, 2.87e-02, 2.95e-02, 1.02e-02, 2.67e-02, 2.69e-02,
                    1.27e-02, 1.02e-02, 1.11e-02, 3.24e-02, 2.90e-02, 3.05e-02, 9.56e-03, 2.70e-02, 7.20e-02,
                    1.07e-02, 1.05e-02, 9.49e-03, 3.42e-02, 2.95e-02, 2.65e-02, 9.21e-03, 2.61e-02, 9.67e-03
                };
                eidAssign(cutdcotdist, cutdcotdist_tmp, 27);
                eidAssign(cutdetain, cutdetain_tmp, 27);
                eidAssign(cutdphiin, cutdphiin_tmp, 27);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 27);
                eidAssign(cutet, cutet_tmp, 27);
                eidAssign(cutfmishits, cutfmishits_tmp, 27);
                eidAssign(cuthoe, cuthoe_tmp, 27);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 27);
                eidAssign(cutiso_sum, cutiso_sum_tmp, 27);
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 27);
                eidAssign(cutsee, cutsee_tmp, 27);
                return;
            }

        case CIC_HYPERTIGHT1:
            {
                double cutdcotdist_tmp[27] = {
                    2.85e-02, 2.37e-02, 1.93e-04, 3.81e-02, 3.26e-02, 2.14e-02, 2.85e-02, 3.82e-02, 3.82e-02,
                    1.61e-02, 2.81e-04, 1.53e-06, 1.26e-02, 2.20e-02, 2.99e-05, 6.79e-03, 7.70e-03, 7.58e-04,
                    6.54e-07, 1.30e-02, 9.88e-05, 1.29e-02, 7.19e-06, 3.93e-03, 9.86e-08, 4.07e-03, 9.28e-04
                };
                double cutdetain_tmp[27] = {
                    7.12e-03, 3.41e-03, 5.47e-03, 8.16e-03, 6.09e-03, 7.79e-03, 9.22e-03, 1.07e-02, 1.14e-02,
                    6.04e-03, 2.19e-03, 3.87e-03, 6.40e-03, 2.70e-03, 4.58e-03, 1.20e-02, 1.15e-02, 1.72e-03,
                    3.21e-03, 1.93e-03, 5.76e-04, 6.87e-03, 1.45e-03, 1.93e-03, 4.18e-03, 4.73e-03, 6.40e-04
                };
                double cutdphiin_tmp[27] = {
                    3.64e-02, 2.63e-02, 1.37e-01, 3.03e-02, 1.50e-02, 5.37e-02, 8.45e-02, 7.06e-02, 1.08e-01,
                    4.58e-02, 2.30e-02, 5.98e-02, 4.28e-02, 6.43e-03, 1.44e-02, 3.39e-02, 2.55e-02, 6.54e-03,
                    4.01e-02, 6.10e-03, 4.02e-02, 2.59e-02, 3.77e-03, 9.49e-03, 3.83e-02, 1.69e-02, 3.46e-02
                };
                double cuteseedopcor_tmp[27] = {
                    8.17e-01, 9.57e-01, 8.91e-01, 9.54e-01, 9.04e-01, 8.30e-01, 6.88e-01, 8.72e-01, 9.61e-01,
                    1.01e+00, 9.84e-01, 1.00e+00, 9.55e-01, 9.52e-01, 1.07e+00, 9.82e-01, 1.07e+00, 1.00e+00,
                    9.79e-01, 9.90e-01, 9.86e-01, 9.61e-01, 9.02e-01, 1.19e+00, 1.01e+00, 9.68e-01, 1.07e+00
                };
                double cutet_tmp[27] = {
                    -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
                    -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
                    1.73e+01, 1.71e+01, 1.84e+01, 1.87e+01, 1.67e+01, 1.76e+01, 1.58e+01, 1.97e+01, 1.86e+01
                };
                double cutfmishits_tmp[27] = {
                    2.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01,
                    5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01,
                    5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01, 5.00e-01, 5.00e-01
                };
                double cuthoe_tmp[27] = {
                    6.66e-02, 2.91e-02, 4.27e-02, 5.06e-02, 2.17e-02, 2.96e-02, 6.72e-02, 4.55e-01, 1.29e-02,
                    7.18e-02, 3.58e-02, 8.26e-02, 2.53e-02, 1.24e-02, 2.46e-02, 5.17e-01, 3.77e-01, 6.12e-01,
                    6.78e-04, 3.01e-02, 1.07e-02, 1.73e-02, 2.36e-03, 3.13e-03, 1.12e-01, 7.45e-01, 6.10e-01
                };
                double cutip_gsf_tmp[27] = {
                    1.61e-02, 1.96e-02, 4.76e-02, 2.30e-02, 4.98e-02, 5.48e-02, 6.85e-02, 4.57e-02, 2.36e-02,
                    1.14e-02, 9.25e-03, 2.25e-02, 1.21e-02, 2.14e-02, 1.46e-02, 2.07e-02, 3.86e-02, 3.74e-03,
                    6.71e-03, 5.56e-03, 6.48e-03, 5.90e-03, 3.78e-03, 6.80e-03, 8.10e-03, 7.78e-03, 2.61e-03
                };
                double cutiso_sum_tmp[27] = {
                    7.98e+00, 5.90e+00, 4.72e+00, 4.12e+00, 3.21e+00, 2.75e+00, 4.06e+00, 3.47e+00, 1.14e+00,
                    1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05,
                    1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05
                };
                double cutiso_sumoet_tmp[27] = {
                    8.49e+00, 7.29e+00, 6.17e+00, 6.04e+00, 3.03e+00, 3.02e+00, 2.77e+00, 5.22e+00, 1.78e+00,
                    7.20e+00, 6.71e+00, 6.29e+00, 5.04e+00, 2.51e+00, 3.74e+00, 9.68e+00, 5.91e+00, 5.79e+00,
                    8.00e+00, 8.97e+00, 6.88e+00, 5.13e+00, 4.13e+00, 4.92e+00, 9.78e+00, 1.13e+01, 6.56e+00
                };
                double cutsee_tmp[27] = {
                    1.26e-02, 1.03e-02, 1.11e-02, 3.21e-02, 2.87e-02, 2.93e-02, 9.38e-03, 2.49e-02, 2.69e-02,
                    1.26e-02, 9.97e-03, 1.02e-02, 3.24e-02, 2.71e-02, 2.59e-02, 9.51e-03, 2.47e-02, 7.20e-02,
                    1.00e-02, 1.00e-02, 9.49e-03, 3.42e-02, 2.95e-02, 2.64e-02, 8.96e-03, 2.60e-02, 9.67e-03
                };
                eidAssign(cutdcotdist, cutdcotdist_tmp, 27);
                eidAssign(cutdetain, cutdetain_tmp, 27);
                eidAssign(cutdphiin, cutdphiin_tmp, 27);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 27);
                eidAssign(cutet, cutet_tmp, 27);
                eidAssign(cutfmishits, cutfmishits_tmp, 27);
                eidAssign(cuthoe, cuthoe_tmp, 27);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 27);
                eidAssign(cutiso_sum, cutiso_sum_tmp, 27);
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 27);
                eidAssign(cutsee, cutsee_tmp, 27);
                return;
            }

        case CIC_HYPERTIGHT2:
            {
                double cutdcotdist_tmp[27] = {
                    2.71e-02, 2.37e-02, 6.96e-05, 1.13e-02, 2.00e-02, 1.07e-02, 2.85e-02, 3.81e-02, 3.82e-02,
                    1.61e-02, 2.81e-04, 1.53e-06, 3.51e-04, 1.32e-02, 2.99e-05, 6.79e-03, 2.16e-04, 7.58e-04,
                    6.54e-07, 1.56e-03, 5.93e-05, 1.67e-03, 4.31e-06, 2.36e-03, 9.86e-08, 4.07e-03, 9.28e-04
                };
                double cutdetain_tmp[27] = {
                    4.91e-03, 3.27e-03, 5.37e-03, 7.39e-03, 5.81e-03, 5.75e-03, 9.22e-03, 7.79e-03, 1.14e-02,
                    4.60e-03, 2.19e-03, 2.92e-03, 4.38e-03, 2.10e-03, 4.58e-03, 8.40e-03, 1.15e-02, 1.72e-03,
                    2.17e-03, 5.96e-04, 5.76e-04, 5.37e-03, 9.71e-04, 3.88e-04, 4.18e-03, 4.73e-03, 6.40e-04
                };
                double cutdphiin_tmp[27] = {
                    3.21e-02, 1.90e-02, 1.05e-01, 2.71e-02, 1.20e-02, 2.52e-02, 8.45e-02, 2.77e-02, 1.08e-01,
                    4.15e-02, 2.30e-02, 2.54e-02, 1.68e-02, 6.31e-03, 1.05e-02, 2.46e-02, 1.83e-02, 6.54e-03,
                    4.42e-03, 2.18e-03, 4.02e-02, 7.00e-03, 3.19e-03, 8.76e-03, 3.83e-02, 1.69e-02, 3.46e-02
                };
                double cuteseedopcor_tmp[27] = {
                    8.95e-01, 9.61e-01, 8.95e-01, 9.85e-01, 9.04e-01, 9.56e-01, 6.88e-01, 8.72e-01, 9.61e-01,
                    1.06e+00, 9.93e-01, 1.00e+00, 1.03e+00, 9.52e-01, 1.17e+00, 1.02e+00, 1.07e+00, 1.00e+00,
                    1.05e+00, 9.90e-01, 9.86e-01, 9.65e-01, 9.02e-01, 1.19e+00, 1.01e+00, 9.68e-01, 1.07e+00
                };
                double cutet_tmp[27] = {
                    -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
                    -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
                    1.92e+01, 1.83e+01, 1.88e+01, 1.90e+01, 1.79e+01, 1.83e+01, 1.58e+01, 1.99e+01, 1.86e+01
                };
                double cutfmishits_tmp[27] = {
                    2.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01,
                    5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01,
                    5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01, 5.00e-01, 5.00e-01
                };
                double cuthoe_tmp[27] = {
                    4.06e-02, 2.91e-02, 2.82e-02, 4.01e-02, 2.06e-02, 2.02e-02, 6.72e-02, 3.29e-01, 1.29e-02,
                    3.71e-02, 2.98e-02, 8.26e-02, 2.33e-02, 1.16e-02, 1.23e-02, 2.37e-01, 1.30e-01, 6.12e-01,
                    4.07e-04, 3.01e-02, 1.07e-02, 1.73e-02, 1.42e-03, 1.88e-03, 1.12e-01, 7.45e-01, 6.10e-01
                };
                double cutip_gsf_tmp[27] = {
                    1.27e-02, 1.11e-02, 3.52e-02, 2.30e-02, 1.86e-02, 5.48e-02, 6.85e-02, 4.57e-02, 2.36e-02,
                    9.67e-03, 8.83e-03, 1.50e-02, 1.21e-02, 2.14e-02, 1.46e-02, 1.89e-02, 3.59e-02, 3.74e-03,
                    3.63e-03, 4.39e-03, 5.17e-03, 3.81e-03, 2.68e-03, 3.45e-03, 8.10e-03, 7.78e-03, 2.61e-03
                };
                double cutiso_sum_tmp[27] = {
                    7.36e+00, 5.40e+00, 3.37e+00, 3.83e+00, 2.13e+00, 2.75e+00, 4.06e+00, 2.97e+00, 1.14e+00,
                    1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05,
                    1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05
                };
                double cutiso_sumoet_tmp[27] = {
                    6.63e+00, 5.05e+00, 4.52e+00, 5.81e+00, 2.46e+00, 2.75e+00, 2.61e+00, 5.22e+00, 1.78e+00,
                    5.96e+00, 5.61e+00, 5.68e+00, 3.48e+00, 2.47e+00, 2.84e+00, 5.64e+00, 5.48e+00, 5.79e+00,
                    7.32e+00, 8.27e+00, 6.37e+00, 4.70e+00, 4.13e+00, 4.92e+00, 9.78e+00, 1.13e+01, 6.56e+00
                };
                double cutsee_tmp[27] = {
                    1.25e-02, 1.03e-02, 1.09e-02, 3.21e-02, 2.72e-02, 2.88e-02, 9.38e-03, 2.40e-02, 2.69e-02,
                    1.15e-02, 9.68e-03, 9.95e-03, 3.24e-02, 2.67e-02, 2.55e-02, 9.51e-03, 2.45e-02, 7.20e-02,
                    9.61e-03, 1.00e-02, 9.38e-03, 3.36e-02, 2.91e-02, 2.63e-02, 8.96e-03, 2.60e-02, 9.67e-03
                };
                eidAssign(cutdcotdist, cutdcotdist_tmp, 27);
                eidAssign(cutdetain, cutdetain_tmp, 27);
                eidAssign(cutdphiin, cutdphiin_tmp, 27);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 27);
                eidAssign(cutet, cutet_tmp, 27);
                eidAssign(cutfmishits, cutfmishits_tmp, 27);
                eidAssign(cuthoe, cuthoe_tmp, 27);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 27);
                eidAssign(cutiso_sum, cutiso_sum_tmp, 27);
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 27);
                eidAssign(cutsee, cutsee_tmp, 27);
                return;
            }

        case CIC_HYPERTIGHT3:
            {
                double cutdcotdist_tmp[27] = {
                    2.28e-02, 2.37e-02, 6.96e-05, 3.15e-04, 2.00e-02, 1.07e-02, 2.85e-02, 3.81e-02, 3.82e-02,
                    7.71e-03, 2.81e-04, 1.53e-06, 3.51e-04, 1.03e-03, 2.33e-06, 6.79e-03, 7.76e-05, 7.58e-04,
                    3.93e-07, 9.36e-04, 5.93e-05, 1.00e-03, 4.31e-06, 2.36e-03, 9.86e-08, 4.07e-03, 9.28e-04
                };
                double cutdetain_tmp[27] = {
                    4.12e-03, 3.27e-03, 3.51e-03, 7.16e-03, 4.52e-03, 5.16e-03, 9.22e-03, 7.79e-03, 1.14e-02,
                    4.43e-03, 1.60e-03, 2.44e-03, 3.82e-03, 2.04e-03, 2.79e-03, 8.40e-03, 1.15e-02, 1.72e-03,
                    1.15e-03, 4.69e-04, 5.76e-04, 5.37e-03, 9.71e-04, 2.54e-04, 4.18e-03, 4.73e-03, 6.40e-04
                };
                double cutdphiin_tmp[27] = {
                    2.32e-02, 1.90e-02, 6.71e-02, 2.71e-02, 1.02e-02, 1.66e-02, 8.45e-02, 2.45e-02, 1.08e-01,
                    3.46e-02, 2.30e-02, 1.89e-02, 1.68e-02, 6.31e-03, 9.39e-03, 2.38e-02, 1.68e-02, 6.54e-03,
                    2.91e-03, 2.18e-03, 4.02e-02, 3.56e-03, 3.19e-03, 8.76e-03, 3.83e-02, 1.69e-02, 3.46e-02
                };
                double cuteseedopcor_tmp[27] = {
                    1.03e+00, 9.61e-01, 9.60e-01, 9.97e-01, 9.25e-01, 9.56e-01, 6.88e-01, 8.72e-01, 9.61e-01,
                    1.10e+00, 9.94e-01, 1.01e+00, 1.04e+00, 9.52e-01, 1.17e+00, 1.03e+00, 1.07e+00, 1.00e+00,
                    1.06e+00, 9.94e-01, 9.86e-01, 9.66e-01, 9.02e-01, 1.19e+00, 1.01e+00, 9.68e-01, 1.07e+00
                };
                double cutet_tmp[27] = {
                    -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
                    -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
                    1.92e+01, 1.83e+01, 1.90e+01, 1.92e+01, 1.88e+01, 1.83e+01, 1.58e+01, 1.99e+01, 1.86e+01
                };
                double cutfmishits_tmp[27] = {
                    5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01,
                    5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01, -5.00e-01, 5.00e-01, -5.00e-01, -5.00e-01,
                    5.00e-01, 1.50e+00, -5.00e-01, -5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01, 5.00e-01, 5.00e-01
                };
                double cuthoe_tmp[27] = {
                    2.99e-02, 2.91e-02, 2.38e-02, 3.98e-02, 2.06e-02, 1.33e-02, 6.72e-02, 3.29e-01, 1.29e-02,
                    3.61e-02, 1.92e-02, 8.26e-02, 2.33e-02, 9.04e-04, 9.57e-04, 2.37e-01, 1.16e-01, 6.12e-01,
                    2.44e-04, 3.01e-02, 1.07e-02, 1.04e-02, 1.42e-03, 1.88e-03, 1.12e-01, 7.45e-01, 6.10e-01
                };
                double cutip_gsf_tmp[27] = {
                    1.13e-02, 9.84e-03, 3.47e-02, 1.77e-02, 1.86e-02, 1.31e-02, 6.85e-02, 4.57e-02, 2.36e-02,
                    8.59e-03, 8.79e-03, 8.97e-03, 5.28e-03, 2.14e-02, 4.82e-03, 1.89e-02, 3.59e-02, 3.74e-03,
                    3.63e-03, 2.03e-03, 5.17e-03, 3.56e-03, 2.68e-03, 3.45e-03, 8.10e-03, 7.78e-03, 2.61e-03
                };
                double cutiso_sum_tmp[27] = {
                    6.70e+00, 4.26e+00, 2.67e+00, 3.68e+00, 2.13e+00, 2.05e+00, 4.06e+00, 2.75e+00, 1.14e+00,
                    1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05,
                    1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05
                };
                double cutiso_sumoet_tmp[27] = {
                    6.04e+00, 3.21e+00, 3.95e+00, 3.46e+00, 2.46e+00, 2.23e+00, 2.61e+00, 5.22e+00, 1.78e+00,
                    5.78e+00, 5.24e+00, 4.69e+00, 3.17e+00, 2.45e+00, 2.80e+00, 5.41e+00, 5.48e+00, 5.79e+00,
                    7.32e+00, 7.66e+00, 6.37e+00, 4.70e+00, 4.13e+00, 4.92e+00, 9.78e+00, 1.13e+01, 6.56e+00
                };
                double cutsee_tmp[27] = {
                    1.18e-02, 1.03e-02, 1.09e-02, 3.09e-02, 2.72e-02, 2.84e-02, 9.14e-03, 2.36e-02, 2.69e-02,
                    9.99e-03, 9.64e-03, 9.91e-03, 3.10e-02, 2.60e-02, 2.53e-02, 9.51e-03, 2.45e-02, 7.20e-02,
                    9.40e-03, 9.98e-03, 9.38e-03, 3.17e-02, 2.91e-02, 2.63e-02, 8.96e-03, 2.60e-02, 9.67e-03
                };
                eidAssign(cutdcotdist, cutdcotdist_tmp, 27);
                eidAssign(cutdetain, cutdetain_tmp, 27);
                eidAssign(cutdphiin, cutdphiin_tmp, 27);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 27);
                eidAssign(cutet, cutet_tmp, 27);
                eidAssign(cutfmishits, cutfmishits_tmp, 27);
                eidAssign(cuthoe, cuthoe_tmp, 27);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 27);
                eidAssign(cutiso_sum, cutiso_sum_tmp, 27);
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 27);
                eidAssign(cutsee, cutsee_tmp, 27);
                return;
            }

        case CIC_HYPERTIGHT4:
            {
                double cutdcotdist_tmp[27] = {
                    6.38e-04, 6.64e-04, 6.96e-05, 8.83e-06, 2.00e-02, 3.05e-03, 2.76e-02, 3.81e-02, 3.82e-02,
                    2.16e-04, 7.87e-06, 1.19e-07, 3.51e-04, 1.03e-03, 2.33e-06, 6.79e-03, 7.76e-05, 7.58e-04,
                    3.93e-07, 9.36e-04, 5.93e-05, 1.00e-03, 4.31e-06, 2.36e-03, 9.86e-08, 4.07e-03, 9.28e-04
                };
                double cutdetain_tmp[27] = {
                    4.12e-03, 2.80e-03, 3.51e-03, 4.65e-03, 2.38e-03, 4.34e-03, 9.22e-03, 7.79e-03, 1.14e-02,
                    2.85e-03, 1.13e-03, 1.12e-03, 3.39e-03, 2.04e-03, 2.79e-03, 8.40e-03, 1.15e-02, 1.72e-03,
                    2.87e-04, 3.54e-04, 5.76e-04, 5.37e-03, 9.71e-04, 2.50e-04, 4.18e-03, 4.73e-03, 6.40e-04
                };
                double cutdphiin_tmp[27] = {
                    2.32e-02, 1.90e-02, 4.38e-02, 2.71e-02, 6.28e-03, 1.07e-02, 8.45e-02, 2.42e-02, 1.08e-01,
                    1.19e-02, 5.36e-03, 1.14e-02, 1.53e-02, 6.31e-03, 9.39e-03, 2.38e-02, 1.68e-02, 6.54e-03,
                    2.44e-03, 2.18e-03, 4.02e-02, 3.56e-03, 3.19e-03, 8.76e-03, 3.83e-02, 1.69e-02, 3.46e-02
                };
                double cuteseedopcor_tmp[27] = {
                    1.03e+00, 9.68e-01, 9.60e-01, 9.97e-01, 9.25e-01, 9.56e-01, 6.88e-01, 8.72e-01, 9.61e-01,
                    1.13e+00, 1.00e+00, 1.01e+00, 1.04e+00, 9.52e-01, 1.17e+00, 1.03e+00, 1.07e+00, 1.00e+00,
                    1.06e+00, 9.94e-01, 9.86e-01, 9.66e-01, 9.02e-01, 1.19e+00, 1.01e+00, 9.68e-01, 1.07e+00
                };
                double cutet_tmp[27] = {
                    -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
                    -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
                    1.93e+01, 1.83e+01, 1.90e+01, 1.92e+01, 1.93e+01, 1.83e+01, 1.58e+01, 1.99e+01, 1.86e+01
                };
                double cutfmishits_tmp[27] = {
                    5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01,
                    5.00e-01, 5.00e-01, -5.00e-01, 5.00e-01, -5.00e-01, -5.00e-01, 5.00e-01, -5.00e-01, -5.00e-01,
                    -5.00e-01, 1.50e+00, -5.00e-01, -5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01, 5.00e-01, 5.00e-01
                };
                double cuthoe_tmp[27] = {
                    2.99e-02, 2.64e-02, 1.35e-02, 2.16e-02, 1.14e-02, 7.32e-03, 6.72e-02, 3.29e-01, 1.29e-02,
                    1.43e-02, 1.49e-03, 6.42e-03, 5.03e-03, 9.04e-04, 9.57e-04, 2.37e-01, 1.16e-01, 6.12e-01,
                    1.46e-04, 3.01e-02, 1.07e-02, 1.04e-02, 1.42e-03, 1.88e-03, 1.12e-01, 7.45e-01, 6.10e-01
                };
                double cutip_gsf_tmp[27] = {
                    1.13e-02, 9.84e-03, 2.42e-02, 1.43e-02, 1.86e-02, 8.63e-03, 6.85e-02, 4.57e-02, 2.36e-02,
                    7.46e-03, 8.79e-03, 3.93e-03, 5.08e-03, 2.14e-02, 4.82e-03, 1.89e-02, 3.59e-02, 3.74e-03,
                    3.03e-03, 8.37e-04, 5.17e-03, 3.56e-03, 2.68e-03, 3.45e-03, 8.10e-03, 7.78e-03, 2.61e-03
                };
                double cutiso_sum_tmp[27] = {
                    4.61e+00, 3.55e+00, 2.39e+00, 3.41e+00, 1.49e+00, 1.87e+00, 4.01e+00, 2.75e+00, 1.14e+00,
                    1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05,
                    1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05
                };
                double cutiso_sumoet_tmp[27] = {
                    4.31e+00, 2.39e+00, 3.41e+00, 3.10e+00, 2.39e+00, 2.23e+00, 2.56e+00, 5.22e+00, 1.78e+00,
                    5.70e+00, 5.24e+00, 4.47e+00, 3.13e+00, 2.45e+00, 2.80e+00, 5.20e+00, 5.48e+00, 5.79e+00,
                    7.32e+00, 7.08e+00, 6.37e+00, 4.70e+00, 4.13e+00, 4.92e+00, 9.78e+00, 1.13e+01, 6.56e+00
                };
                double cutsee_tmp[27] = {
                    1.09e-02, 1.01e-02, 1.02e-02, 2.99e-02, 2.64e-02, 2.84e-02, 9.11e-03, 2.36e-02, 2.69e-02,
                    9.27e-03, 9.51e-03, 9.18e-03, 2.71e-02, 2.60e-02, 2.53e-02, 9.51e-03, 2.45e-02, 7.20e-02,
                    9.40e-03, 9.98e-03, 9.38e-03, 3.17e-02, 2.91e-02, 2.63e-02, 8.96e-03, 2.60e-02, 9.67e-03
                };
                eidAssign(cutdcotdist, cutdcotdist_tmp, 18);
                eidAssign(cutdetain, cutdetain_tmp, 18);
                eidAssign(cutdphiin, cutdphiin_tmp, 18);
                eidAssign(cuteseedopcor, cuteseedopcor_tmp, 18);
                eidAssign(cutet, cutet_tmp, 18);
                eidAssign(cutfmishits, cutfmishits_tmp, 18);
                eidAssign(cuthoe, cuthoe_tmp, 18);
                eidAssign(cutip_gsf, cutip_gsf_tmp, 18);
                eidAssign(cutiso_sum, cutiso_sum_tmp, 18);
                eidAssign(cutiso_sumoet, cutiso_sumoet_tmp, 18);
                eidAssign(cutsee, cutsee_tmp, 18);
                return;
            }

        default:
            std::cout << "[eidGetCIC] ERROR! Invalid tightness level" << std::endl;

    }

}

void eidGetCIC_V02(const cic_tightness tightness, std::vector<double> &cutdeta, std::vector<double> &cutdphi, std::vector<double> &cuteopin, std::vector<double> &cutet, std::vector<double> &cuthoe, std::vector<double> &cutip, std::vector<double> &cutisoecal, std::vector<double> &cutisohcal, std::vector<double> &cutisotk, std::vector<double> &cutmishits, std::vector<double> &cutsee)
{

    switch (tightness) {
        case CIC_VERYLOOSE:
            {
                double cutdeta_tmp[18] = {
                    9.59e-03, 5.11e-03, 1.46e-02, 1.37e-02, 9.65e-03, 1.48e-02,
                    1.14e-02, 4.41e-03, 1.06e-02, 1.50e-02, 7.62e-03, 1.33e-02,
                    1.39e-02, 4.14e-03, 1.38e-02, 1.49e-02, 6.86e-03, 1.32e-02};
                double cutdphi_tmp[18] = {
                    3.75e-02, 1.16e-01, 1.19e-01, 4.88e-02, 1.19e-01, 1.19e-01,
                    6.69e-02, 7.52e-02, 1.19e-01, 7.62e-02, 9.95e-02, 1.20e-01,
                    9.06e-02, 5.08e-02, 1.19e-01, 6.71e-02, 3.39e-02, 1.19e-01};
                double cuteopin_tmp[18] = {
                    8.78e-01, 8.02e-01, 8.14e-01, 9.42e-01, 7.35e-01, 7.74e-01,
                    8.29e-01, 8.30e-01, 8.05e-01, 7.99e-01, 7.02e-01, 7.88e-01,
                    8.18e-01, 8.21e-01, 8.01e-01, 8.47e-01, 7.16e-01, 7.87e-01};
                double cutet_tmp[18] = {
                    0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
                    0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
                    1.20e+01, 1.20e+01, 1.20e+01, 1.20e+01, 1.21e+01, 1.20e+01};
                double cuthoe_tmp[18] = {
                    8.96e-02, 9.63e-02, 9.82e-02, 1.01e-01, 6.65e-02, 9.21e-02,
                    9.98e-02, 6.15e-02, 1.02e-01, 1.01e-01, 4.31e-02, 9.94e-02,
                    7.88e-02, 7.14e-02, 9.97e-02, 9.84e-02, 1.89e-02, 9.84e-02};
                double cutip_tmp[18] = {
                    2.53e-02, 7.64e-02, 9.69e-02, 8.91e-02, 4.45e-01, 2.71e-01,
                    3.71e-02, 3.92e-02, 8.79e-02, 4.67e-02, 5.00e-01, 3.08e-01,
                    7.33e+00, 1.88e-02, 6.04e+00, 7.31e+00, 2.02e+00, 5.77e+00};
                double cutisoecal_tmp[18] = {
                    3.34e+01, 2.87e+01, 1.16e+01, 2.74e+01, 1.86e+01, 2.47e+01,
                    1.04e+02, 1.02e+02, 2.28e+01, 7.07e+01, 3.84e+01, 5.37e+01,
                    2.52e+01, 2.70e+02, 2.97e+01, 5.30e+01, 3.20e+01, 1.59e+01};
                double cutisohcal_tmp[18] = {
                    1.35e+01, 1.27e+01, 1.14e+01, 1.48e+01, 9.86e+00, 1.76e+01,
                    4.87e+01, 2.35e+01, 2.15e+01, 5.05e+01, 6.27e+00, 7.36e+00,
                    1.89e+01, 5.36e+01, 3.87e+01, 1.74e+01, 5.00e+01, 3.45e+00};
                double cutisotk_tmp[18] = {
                    2.44e+01, 1.89e+01, 1.94e+01, 2.78e+01, 9.63e+00, 1.74e+01,
                    4.11e+01, 1.46e+01, 2.18e+01, 6.36e+01, 9.11e+00, 1.78e+01,
                    1.68e+01, 1.19e+01, 1.42e+01, 2.32e+01, 1.36e+01, 2.37e+01};
                double cutmishits_tmp[18] = {
                    5.50e+00, 1.50e+00, 5.50e+00, 4.50e+00, 2.50e+00, 2.50e+00,
                    3.50e+00, 5.50e+00, 7.50e+00, 4.50e+00, 2.50e+00, 1.50e+00,
                    1.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
                double cutsee_tmp[18] = {
                    1.73e-02, 1.34e-02, 2.32e-02, 3.45e-02, 3.15e-02, 3.45e-02,
                    1.83e-02, 1.16e-02, 1.40e-02, 3.49e-02, 3.09e-02, 3.38e-02,
                    1.78e-02, 1.14e-02, 1.42e-02, 3.50e-02, 3.18e-02, 3.47e-02};
                eidAssign(cutdeta, cutdeta_tmp, 18);
                eidAssign(cutdphi, cutdphi_tmp, 18);
                eidAssign(cuteopin, cuteopin_tmp, 18);
                eidAssign(cutet, cutet_tmp, 18);
                eidAssign(cuthoe, cuthoe_tmp, 18);
                eidAssign(cutip, cutip_tmp, 18);
                eidAssign(cutisoecal, cutisoecal_tmp, 18);
                eidAssign(cutisohcal, cutisohcal_tmp, 18);
                eidAssign(cutisotk, cutisotk_tmp, 18);
                eidAssign(cutmishits, cutmishits_tmp, 18);
                eidAssign(cutsee, cutsee_tmp, 18);
                return;
            }
        case CIC_LOOSE:
            {   
                double cutdeta_tmp[18] = {
                    9.58e-03, 4.06e-03, 1.22e-02, 1.37e-02, 8.37e-03, 1.27e-02,
                    1.10e-02, 3.36e-03, 9.77e-03, 1.50e-02, 6.75e-03, 1.09e-02,
                    1.40e-02, 5.08e-03, 1.09e-02, 1.46e-02, 5.06e-03, 1.27e-02};
                double cutdphi_tmp[18] = {
                    3.72e-02, 1.14e-01, 1.18e-01, 4.88e-02, 1.17e-01, 1.19e-01,
                    6.06e-02, 5.48e-02, 1.17e-01, 7.00e-02, 3.55e-02, 1.17e-01,
                    8.80e-02, 4.50e-02, 1.18e-01, 9.19e-02, 2.36e-02, 5.15e-02};
                double cuteopin_tmp[18] = {
                    8.78e-01, 8.02e-01, 8.14e-01, 9.42e-01, 7.35e-01, 7.74e-01,
                    8.29e-01, 9.09e-01, 8.29e-01, 8.13e-01, 8.60e-01, 8.97e-01,
                    8.17e-01, 8.31e-01, 8.18e-01, 8.61e-01, 7.87e-01, 7.89e-01};
                double cutet_tmp[18] = {
                    0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
                    0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
                    1.20e+01, 1.20e+01, 1.29e+01, 1.21e+01, 1.25e+01, 1.29e+01};
                double cuthoe_tmp[18] = {
                    8.87e-02, 9.34e-02, 9.49e-02, 9.86e-02, 4.31e-02, 8.78e-02,
                    9.70e-02, 5.09e-02, 9.80e-02, 9.91e-02, 3.21e-02, 9.28e-02,
                    6.63e-02, 7.17e-02, 9.66e-02, 7.58e-02, 1.49e-02, 1.31e-02};
                double cutip_tmp[18] = {
                    2.46e-02, 7.60e-02, 9.66e-02, 8.85e-02, 4.41e-01, 2.05e-01,
                    2.92e-02, 2.93e-02, 6.19e-02, 2.51e-02, 1.59e-01, 8.15e-02,
                    7.29e+00, 1.06e-02, 5.76e+00, 6.89e+00, 1.27e+00, 5.89e+00};
                double cutisoecal_tmp[18] = {
                    3.34e+01, 2.81e+01, 7.32e+00, 2.74e+01, 7.33e+00, 2.17e+01,
                    9.38e+01, 1.02e+02, 1.21e+01, 2.60e+01, 8.91e+00, 1.00e+01,
                    1.61e+01, 3.13e+01, 1.69e+01, 1.54e+01, 1.33e+01, 3.77e+01};
                double cutisohcal_tmp[18] = {
                    1.35e+01, 9.93e+00, 7.56e+00, 1.48e+01, 8.10e+00, 1.08e+01,
                    4.27e+01, 2.01e+01, 9.11e+00, 1.04e+01, 6.89e+00, 5.59e+00,
                    8.53e+00, 9.59e+00, 2.42e+01, 2.78e+00, 8.67e+00, 2.88e-01};
                double cutisotk_tmp[18] = {
                    2.43e+01, 8.45e+00, 1.44e+01, 2.78e+01, 6.02e+00, 1.05e+01,
                    1.41e+01, 1.02e+01, 1.45e+01, 1.91e+01, 6.10e+00, 1.41e+01,
                    8.59e+00, 8.33e+00, 8.30e+00, 8.93e+00, 8.60e+00, 1.60e+01};
                double cutmishits_tmp[18] = {
                    5.50e+00, 1.50e+00, 5.50e+00, 2.50e+00, 2.50e+00, 2.50e+00,
                    3.50e+00, 5.50e+00, 5.00e-01, 1.50e+00, 2.50e+00, 5.00e-01,
                    1.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
                double cutsee_tmp[18] = {
                    1.72e-02, 1.15e-02, 1.43e-02, 3.44e-02, 2.95e-02, 3.04e-02,
                    1.45e-02, 1.08e-02, 1.28e-02, 3.47e-02, 3.07e-02, 3.16e-02,
                    1.80e-02, 1.10e-02, 1.32e-02, 3.49e-02, 3.10e-02, 3.27e-02};
                eidAssign(cutdeta, cutdeta_tmp, 18);
                eidAssign(cutdphi, cutdphi_tmp, 18);
                eidAssign(cuteopin, cuteopin_tmp, 18);
                eidAssign(cutet, cutet_tmp, 18);
                eidAssign(cuthoe, cuthoe_tmp, 18);
                eidAssign(cutip, cutip_tmp, 18);
                eidAssign(cutisoecal, cutisoecal_tmp, 18);
                eidAssign(cutisohcal, cutisohcal_tmp, 18);
                eidAssign(cutisotk, cutisotk_tmp, 18);
                eidAssign(cutmishits, cutmishits_tmp, 18);
                eidAssign(cutsee, cutsee_tmp, 18);
                return;
            }
        case CIC_MEDIUM:
            {   
                double cutdeta_tmp[18] = {
                    9.37e-03, 3.28e-03, 7.70e-03, 1.37e-02, 6.57e-03, 8.60e-03,
                    1.02e-02, 2.86e-03, 8.26e-03, 1.18e-02, 7.03e-03, 8.67e-03,
                    1.23e-02, 2.67e-03, 1.01e-02, 1.45e-02, 4.98e-03, 1.28e-02};
                double cutdphi_tmp[18] = {
                    3.72e-02, 5.32e-02, 1.18e-01, 4.87e-02, 6.60e-02, 1.19e-01,
                    5.72e-02, 3.40e-02, 1.12e-01, 6.91e-02, 1.82e-02, 1.13e-01,
                    8.00e-02, 8.89e-02, 7.75e-02, 4.54e-02, 1.61e-02, 4.15e-02};
                double cuteopin_tmp[18] = {
                    8.78e-01, 8.07e-01, 8.43e-01, 9.42e-01, 7.35e-01, 7.74e-01,
                    8.27e-01, 9.54e-01, 8.50e-01, 8.34e-01, 9.04e-01, 9.71e-01,
                    8.23e-01, 9.44e-01, 8.23e-01, 8.46e-01, 7.71e-01, 1.09e+00};
                double cutet_tmp[18] = {
                    0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
                    0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
                    1.32e+01, 1.39e+01, 1.50e+01, 1.45e+01, 1.52e+01, 1.53e+01};
                double cuthoe_tmp[18] = {
                    8.86e-02, 9.02e-02, 9.39e-02, 9.81e-02, 3.43e-02, 5.64e-02,
                    7.45e-02, 4.83e-02, 8.09e-02, 9.60e-02, 2.42e-02, 8.72e-02,
                    6.44e-03, 8.80e-02, 4.03e-02, 2.51e-02, 1.51e-02, 2.07e-02};
                double cutip_tmp[18] = {
                    2.45e-02, 5.64e-02, 8.79e-02, 8.59e-02, 3.91e-01, 1.61e-01,        
                    1.22e-02, 1.95e-02, 3.77e-02, 2.10e-02, 1.26e-01, 4.50e-02,        
                    7.33e+00, 9.38e-03, 5.28e+00, 6.52e+00, 7.69e-01, 5.90e+00};
                double cutisoecal_tmp[18] = {
                    2.87e+01, 2.81e+01, 5.29e+00, 2.10e+01, 4.11e+00, 4.86e+00,
                    2.79e+01, 2.67e+01, 8.86e+00, 1.47e+01, 4.94e+00, 6.09e+00,
                    1.12e+01, 1.49e+01, 1.36e+01, 1.31e+01, 1.06e+01, 1.34e+01};
                double cutisohcal_tmp[18] = {        
                    1.15e+01, 7.62e+00, 8.80e+00, 5.23e+00, 7.17e+00, 2.55e+00,        
                    3.17e+01, 2.12e+01, 6.10e+00, 6.03e+00, 8.84e+00, 3.71e+00,
                    3.52e+00, 5.94e+00, 2.00e+01, 1.29e-01, 4.10e+00, 1.35e-02};
                double cutisotk_tmp[18] = {
                    1.00e+01, 6.04e+00, 8.94e+00, 1.29e+01, 4.14e+00, 8.03e+00, 
                    8.51e+00, 6.73e+00, 7.49e+00, 1.19e+01, 3.59e+00, 6.70e+00,
                    5.88e+00, 5.95e+00, 6.60e+00, 5.87e+00, 6.21e+00, 7.44e-01};
                double cutmishits_tmp[18] = {
                    5.50e+00, 1.50e+00, 5.50e+00, 2.50e+00, 2.50e+00, 1.50e+00,
                    3.50e+00, 5.50e+00, 5.00e-01, 1.50e+00, 1.50e+00, 5.00e-01,
                    5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
                double cutsee_tmp[18] = {
                    1.69e-02, 1.07e-02, 1.20e-02, 3.18e-02, 2.83e-02, 2.89e-02,
                    1.32e-02, 1.06e-02, 1.19e-02, 3.29e-02, 2.96e-02, 3.00e-02, 
                    1.66e-02, 1.12e-02, 1.15e-02, 3.48e-02, 2.96e-02, 3.33e-02};
                eidAssign(cutdeta, cutdeta_tmp, 18);
                eidAssign(cutdphi, cutdphi_tmp, 18);
                eidAssign(cuteopin, cuteopin_tmp, 18);
                eidAssign(cutet, cutet_tmp, 18);
                eidAssign(cuthoe, cuthoe_tmp, 18);
                eidAssign(cutip, cutip_tmp, 18);
                eidAssign(cutisoecal, cutisoecal_tmp, 18);
                eidAssign(cutisohcal, cutisohcal_tmp, 18);
                eidAssign(cutisotk, cutisotk_tmp, 18);
                eidAssign(cutmishits, cutmishits_tmp, 18);
                eidAssign(cutsee, cutsee_tmp, 18);
                return;
            }
        case CIC_TIGHT:
            {   
                double cutdeta_tmp[18] = {
                    9.15e-03, 3.02e-03, 6.10e-03, 1.35e-02, 5.65e-03, 7.93e-03,
                    1.02e-02, 2.66e-03, 1.06e-02, 9.03e-03, 7.66e-03, 7.23e-03,
                    1.16e-02, 2.03e-03, 6.59e-03, 1.48e-02, 5.55e-03, 1.28e-02};
                double cutdphi_tmp[18] = {
                    3.69e-02, 3.07e-02, 1.17e-01, 4.75e-02, 2.16e-02, 1.17e-01,
                    3.72e-02, 2.46e-02, 4.26e-02, 6.12e-02, 1.42e-02, 3.90e-02,
                    7.37e-02, 5.66e-02, 3.59e-02, 1.87e-02, 1.20e-02, 3.58e-02};
                double cuteopin_tmp[18] = {
                    8.78e-01, 8.59e-01, 8.74e-01, 9.44e-01, 7.37e-01, 7.73e-01,
                    8.60e-01, 9.67e-01, 9.17e-01, 8.12e-01, 9.15e-01, 1.01e+00,
                    8.47e-01, 9.53e-01, 9.79e-01, 8.41e-01, 7.71e-01, 1.09e+00};
                double cutet_tmp[18] = {
                    0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
                    0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
                    1.46e+01, 1.44e+01, 1.61e+01, 1.52e+01, 1.53e+01, 1.61e+01};
                double cuthoe_tmp[18] = {
                    8.71e-02, 2.89e-02, 7.83e-02, 9.46e-02, 2.45e-02, 3.63e-02,
                    6.71e-02, 4.80e-02, 6.14e-02, 9.24e-02, 1.58e-02, 4.90e-02,
                    3.82e-02, 9.15e-02, 4.51e-02, 4.52e-02, 1.96e-03, 4.30e-03};
                double cutip_tmp[18] = {
                    2.39e-02, 2.70e-02, 7.68e-02, 2.31e-02, 1.78e-01, 9.57e-02,
                    1.02e-02, 1.68e-02, 4.30e-02, 1.66e-02, 5.94e-02, 3.08e-02,
                    2.10e+00, 5.27e-03, 3.17e+00, 4.91e+00, 7.69e-01, 5.90e+00};
                double cutisoecal_tmp[18] = {
                    2.00e+01, 2.72e+01, 4.48e+00, 1.35e+01, 4.56e+00, 3.19e+00,
                    1.22e+01, 1.31e+01, 7.42e+00, 7.67e+00, 4.12e+00, 4.85e+00,
                    1.01e+01, 1.24e+01, 1.11e+01, 1.10e+01, 1.06e+01, 1.34e+01};
                double cutisohcal_tmp[18] = {
                    1.09e+01, 7.01e+00, 8.75e+00, 3.51e+00, 7.75e+00, 1.62e+00,
                    1.16e+01, 9.90e+00, 4.97e+00, 5.33e+00, 3.18e+00, 2.32e+00,
                    1.64e-01, 5.46e+00, 1.20e+01, 6.04e-03, 4.10e+00, 6.28e-04};
                double cutisotk_tmp[18] = {
                    6.53e+00, 4.60e+00, 6.00e+00, 8.63e+00, 3.11e+00, 7.77e+00,
                    5.42e+00, 4.81e+00, 4.06e+00, 6.47e+00, 2.80e+00, 3.45e+00,
                    5.29e+00, 5.18e+00, 1.54e+01, 5.38e+00, 4.47e+00, 3.47e-02};
                double cutmishits_tmp[18] = {
                    5.50e+00, 1.50e+00, 5.00e-01, 1.50e+00, 2.50e+00, 5.00e-01,
                    3.50e+00, 5.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01,
                    5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
                double cutsee_tmp[18] = {
                    1.31e-02, 1.06e-02, 1.15e-02, 3.06e-02, 2.80e-02, 2.93e-02,
                    1.31e-02, 1.06e-02, 1.15e-02, 3.17e-02, 2.90e-02, 2.89e-02,
                    1.42e-02, 1.06e-02, 1.03e-02, 3.50e-02, 2.96e-02, 3.33e-02};
                eidAssign(cutdeta, cutdeta_tmp, 18);
                eidAssign(cutdphi, cutdphi_tmp, 18);
                eidAssign(cuteopin, cuteopin_tmp, 18);
                eidAssign(cutet, cutet_tmp, 18);
                eidAssign(cuthoe, cuthoe_tmp, 18);
                eidAssign(cutip, cutip_tmp, 18);
                eidAssign(cutisoecal, cutisoecal_tmp, 18);
                eidAssign(cutisohcal, cutisohcal_tmp, 18);
                eidAssign(cutisotk, cutisotk_tmp, 18);
                eidAssign(cutmishits, cutmishits_tmp, 18);
                eidAssign(cutsee, cutsee_tmp, 18);
                return;
            }
        case CIC_SUPERTIGHT:
            {   
                double cutdeta_tmp[18] = {
                    8.92e-03, 2.77e-03, 5.35e-03, 1.19e-02, 5.21e-03, 6.88e-03,
                    9.66e-03, 2.66e-03, 1.02e-02, 7.68e-03, 5.22e-03, 5.91e-03,
                    1.16e-02, 2.03e-03, 6.59e-03, 1.48e-02, 5.55e-03, 1.28e-02};
                double cutdphi_tmp[18] = {
                    3.70e-02, 2.06e-02, 1.13e-01, 4.59e-02, 1.46e-02, 1.12e-01,
                    3.35e-02, 1.99e-02, 2.98e-02, 3.07e-02, 1.02e-02, 1.94e-02,
                    7.37e-02, 5.66e-02, 3.59e-02, 1.87e-02, 1.20e-02, 3.58e-02};
                double cuteopin_tmp[18] = {
                    8.78e-01, 9.35e-01, 8.87e-01, 9.50e-01, 8.61e-01, 7.65e-01,
                    8.99e-01, 9.75e-01, 9.69e-01, 7.99e-01, 9.47e-01, 1.00e+00,
                    8.47e-01, 9.53e-01, 9.79e-01, 8.41e-01, 7.71e-01, 1.09e+00};
                double cutet_tmp[18] = {
                    0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
                    0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
                    1.46e+01, 1.44e+01, 1.61e+01, 1.52e+01, 1.53e+01, 1.61e+01};
                double cuthoe_tmp[18] = {
                    8.13e-02, 2.19e-02, 2.82e-02, 5.87e-02, 1.94e-02, 1.66e-02,
                    5.51e-02, 3.56e-02, 7.75e-02, 7.83e-02, 2.14e-02, 3.76e-02,
                    3.82e-02, 9.15e-02, 4.51e-02, 4.52e-02, 1.96e-03, 4.30e-03};
                double cutip_tmp[18] = {
                    2.37e-02, 1.97e-02, 7.91e-02, 1.38e-02, 6.88e-02, 4.15e-02,
                    9.36e-03, 9.66e-03, 1.45e-02, 1.37e-02, 6.01e-02, 1.32e-02,
                    2.10e+00, 5.27e-03, 3.17e+00, 4.91e+00, 7.69e-01, 5.90e+00};
                double cutisoecal_tmp[18] = {
                    9.09e+00, 1.18e+01, 3.92e+00, 4.71e+00, 1.34e+01, 2.58e+00,
                    8.75e+00, 8.10e+00, 6.41e+00, 5.70e+00, 3.47e+00, 4.17e+00,
                    1.01e+01, 1.24e+01, 1.11e+01, 1.10e+01, 1.06e+01, 1.34e+01};
                double cutisohcal_tmp[18] = {
                    9.93e+00, 8.32e+00, 2.76e+00, 2.98e+00, 3.81e+00, 1.08e+00,
                    6.66e+00, 6.67e+00, 2.79e+00, 3.10e+00, 1.56e+00, 2.00e+00,
                    1.64e-01, 5.46e+00, 1.20e+01, 6.04e-03, 4.10e+00, 6.28e-04};
                double cutisotk_tmp[18] = {
                    5.02e+00, 3.71e+00, 3.80e+00, 7.73e+00, 2.13e+00, 8.76e+00,
                    4.46e+00, 3.73e+00, 2.79e+00, 3.25e+00, 2.51e+00, 2.18e+00,
                    5.29e+00, 5.18e+00, 1.54e+01, 5.38e+00, 4.47e+00, 3.47e-02};
                double cutmishits_tmp[18] = {
                    5.50e+00, 1.50e+00, 5.00e-01, 1.50e+00, 1.50e+00, 5.00e-01,
                    5.00e-01, 5.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01,
                    5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
                double cutsee_tmp[18] = {
                    1.18e-02, 1.05e-02, 1.12e-02, 3.05e-02, 2.76e-02, 2.86e-02,
                    1.25e-02, 1.01e-02, 1.08e-02, 3.37e-02, 2.77e-02, 2.93e-02,
                    1.42e-02, 1.06e-02, 1.03e-02, 3.50e-02, 2.96e-02, 3.33e-02};
                eidAssign(cutdeta, cutdeta_tmp, 18);
                eidAssign(cutdphi, cutdphi_tmp, 18);
                eidAssign(cuteopin, cuteopin_tmp, 18);
                eidAssign(cutet, cutet_tmp, 18);
                eidAssign(cuthoe, cuthoe_tmp, 18);
                eidAssign(cutip, cutip_tmp, 18);
                eidAssign(cutisoecal, cutisoecal_tmp, 18);
                eidAssign(cutisohcal, cutisohcal_tmp, 18);
                eidAssign(cutisotk, cutisotk_tmp, 18);
                eidAssign(cutmishits, cutmishits_tmp, 18);
                eidAssign(cutsee, cutsee_tmp, 18);
                return;
            }
        case CIC_HYPERTIGHT1:
            {   
                double cutdeta_tmp[18] = {
                    8.81e-03, 2.39e-03, 5.60e-03, 8.00e-03, 5.30e-03, 6.18e-03,
                    9.96e-03, 2.43e-03, 1.06e-02, 7.83e-03, 6.12e-03, 4.94e-03,
                    1.16e-02, 2.03e-03, 6.59e-03, 1.48e-02, 5.55e-03, 1.28e-02};
                double cutdphi_tmp[18] = {
                    3.64e-02, 1.61e-02, 4.33e-02, 4.19e-02, 1.34e-02, 4.40e-02,
                    2.75e-02, 2.22e-02, 2.17e-02, 2.04e-02, 9.14e-03, 1.13e-02,
                    7.37e-02, 5.66e-02, 3.59e-02, 1.87e-02, 1.20e-02, 3.58e-02};
                double cuteopin_tmp[18] = {
                    8.79e-01, 9.65e-01, 9.04e-01, 9.50e-01, 7.91e-01, 7.74e-01,
                    9.20e-01, 9.83e-01, 9.89e-01, 7.98e-01, 9.70e-01, 1.02e+00,
                    8.47e-01, 9.53e-01, 9.79e-01, 8.41e-01, 7.71e-01, 1.09e+00};
                double cutet_tmp[18] = {
                    0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
                    0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
                    1.46e+01, 1.44e+01, 1.61e+01, 1.52e+01, 1.53e+01, 1.61e+01};
                double cuthoe_tmp[18] = {        
                    5.29e-02, 1.71e-02, 1.95e-02, 3.93e-02, 1.38e-02, 1.13e-02,        
                    5.01e-02, 2.99e-02, 5.38e-02, 5.44e-02, 2.05e-02, 1.35e-02,
                    3.82e-02, 9.15e-02, 4.51e-02, 4.52e-02, 1.96e-03, 4.30e-03};
                double cutip_tmp[18] = {
                    2.38e-02, 2.02e-02, 4.23e-02, 1.39e-02, 2.00e-02, 3.00e-02,
                    8.61e-03, 1.50e-02, 1.08e-02, 1.07e-02, 1.65e-02, 1.06e-02,
                    2.10e+00, 5.27e-03, 3.17e+00, 4.91e+00, 7.69e-01, 5.90e+00};    
                double cutisoecal_tmp[18] = {        
                    5.42e+00, 9.14e+00, 3.65e+00, 3.92e+00, 2.90e+00, 2.45e+00,
                    6.92e+00, 6.49e+00, 5.08e+00, 5.26e+00, 2.71e+00, 3.57e+00,
                    1.01e+01, 1.24e+01, 1.11e+01, 1.10e+01, 1.06e+01, 1.34e+01};
                double cutisohcal_tmp[18] = {
                    9.95e+00, 8.15e+00, 4.13e+00, 2.98e+00, 3.22e+00, 9.25e-01,
                    5.35e+00, 5.67e+00, 2.45e+00, 1.91e+00, 1.01e+00, 1.78e+00,
                    1.64e-01, 5.46e+00, 1.20e+01, 6.04e-03, 4.10e+00, 6.28e-04};
                double cutisotk_tmp[18] = {
                    4.26e+00, 3.26e+00, 2.34e+00, 8.33e+00, 1.40e+00, 3.45e+00,
                    3.35e+00, 3.16e+00, 2.37e+00, 2.22e+00, 2.14e+00, 1.58e+00,
                    5.29e+00, 5.18e+00, 1.54e+01, 5.38e+00, 4.47e+00, 3.47e-02};
                double cutmishits_tmp[18] = {
                    5.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01,
                    5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01,
                    5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
                double cutsee_tmp[18] = {
                    1.13e-02, 1.06e-02, 1.08e-02, 2.96e-02, 2.73e-02, 2.86e-02,
                    1.25e-02, 1.00e-02, 1.05e-02, 3.40e-02, 2.69e-02, 2.87e-02, 
                    1.42e-02, 1.06e-02, 1.03e-02, 3.50e-02, 2.96e-02, 3.33e-02};
                eidAssign(cutdeta, cutdeta_tmp, 18);
                eidAssign(cutdphi, cutdphi_tmp, 18);
                eidAssign(cuteopin, cuteopin_tmp, 18);
                eidAssign(cutet, cutet_tmp, 18);
                eidAssign(cuthoe, cuthoe_tmp, 18);
                eidAssign(cutip, cutip_tmp, 18);
                eidAssign(cutisoecal, cutisoecal_tmp, 18);
                eidAssign(cutisohcal, cutisohcal_tmp, 18);
                eidAssign(cutisotk, cutisotk_tmp, 18);
                eidAssign(cutmishits, cutmishits_tmp, 18);
                eidAssign(cutsee, cutsee_tmp, 18);
                return;
            }
        case CIC_HYPERTIGHT2:
            {   
                double cutdeta_tmp[18] = {
                    8.53e-03, 2.22e-03, 5.15e-03, 6.94e-03, 6.10e-03, 6.05e-03,
                    8.34e-03, 1.88e-03, 2.46e-03, 4.61e-03, 3.20e-03, 3.28e-03,
                    1.16e-02, 2.03e-03, 6.59e-03, 1.48e-02, 5.55e-03, 1.28e-02};
                double cutdphi_tmp[18] = {
                    2.01e-02, 1.21e-02, 2.13e-02, 4.45e-02, 1.15e-02, 1.70e-02,
                    3.62e-02, 1.93e-02, 1.37e-02, 1.81e-02, 6.43e-03, 9.47e-03,
                    7.37e-02, 5.66e-02, 3.59e-02, 1.87e-02, 1.20e-02, 3.58e-02};
                double cuteopin_tmp[18] = {
                    8.79e-01, 9.81e-01, 8.77e-01, 9.61e-01, 9.29e-01, 7.84e-01,
                    8.88e-01, 1.01e+00, 1.02e+00, 9.12e-01, 1.01e+00, 1.03e+00,
                    8.47e-01, 9.53e-01, 9.79e-01, 8.41e-01, 7.71e-01, 1.09e+00};
                double cutet_tmp[18] = {
                    0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
                    0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
                    1.46e+01, 1.44e+01, 1.61e+01, 1.52e+01, 1.53e+01, 1.61e+01};
                double cuthoe_tmp[18] = {
                    2.91e-02, 2.35e-02, 1.84e-02, 2.58e-02, 1.19e-02, 1.02e-02,
                    4.38e-02, 2.61e-03, 5.73e-02, 4.52e-02, 5.93e-02, 6.31e-04,
                    3.82e-02, 9.15e-02, 4.51e-02, 4.52e-02, 1.96e-03, 4.30e-03};
                double cutip_tmp[18] = {
                    1.95e-02, 1.97e-02, 1.58e-02, 1.48e-02, 1.31e-02, 1.80e-02,
                    7.56e-03, 1.81e-02, 1.11e-02, 9.29e-03, 2.39e-02, 6.64e-03,
                    2.10e+00, 5.27e-03, 3.17e+00, 4.91e+00, 7.69e-01, 5.90e+00};
                double cutisoecal_tmp[18] = {
                    6.72e+00, 5.83e+00, 3.59e+00, 2.91e+00, 2.63e+00, 2.54e+00,
                    5.28e+00, 5.43e+00, 4.82e+00, 3.68e+00, 2.54e+00, 3.25e+00,
                    1.01e+01, 1.24e+01, 1.11e+01, 1.10e+01, 1.06e+01, 1.34e+01};
                double cutisohcal_tmp[18] = {
                    8.30e+00, 8.23e+00, 4.99e+00, 2.16e+00, 3.45e+00, 7.47e+00,
                    3.12e+00, 4.67e+00, 2.82e+00, 1.12e+00, 9.91e-01, 1.59e+00,
                    1.64e-01, 5.46e+00, 1.20e+01, 6.04e-03, 4.10e+00, 6.28e-04};
                double cutisotk_tmp[18] = {
                    3.30e+00, 2.34e+00, 1.74e+00, 7.26e+00, 9.64e-01, 1.35e+00,
                    2.45e+00, 2.66e+00, 1.89e+00, 1.83e+00, 5.53e-01, 7.39e-02,
                    5.29e+00, 5.18e+00, 1.54e+01, 5.38e+00, 4.47e+00, 3.47e-02};
                double cutmishits_tmp[18] = {
                    5.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01,
                    5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01,
                    5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
                double cutsee_tmp[18] = {
                    1.10e-02, 1.01e-02, 1.09e-02, 2.99e-02, 2.72e-02, 2.77e-02,
                    1.07e-02, 9.80e-03, 1.01e-02, 3.01e-02, 2.71e-02, 2.62e-02,
                    1.42e-02, 1.06e-02, 1.03e-02, 3.50e-02, 2.96e-02, 3.33e-02};
                eidAssign(cutdeta, cutdeta_tmp, 18);
                eidAssign(cutdphi, cutdphi_tmp, 18);
                eidAssign(cuteopin, cuteopin_tmp, 18);
                eidAssign(cutet, cutet_tmp, 18);
                eidAssign(cuthoe, cuthoe_tmp, 18);
                eidAssign(cutip, cutip_tmp, 18);
                eidAssign(cutisoecal, cutisoecal_tmp, 18);
                eidAssign(cutisohcal, cutisohcal_tmp, 18);
                eidAssign(cutisotk, cutisotk_tmp, 18);
                eidAssign(cutmishits, cutmishits_tmp, 18);
                eidAssign(cutsee, cutsee_tmp, 18);
                return;
            }
        case CIC_HYPERTIGHT3:
            {   
                double cutdeta_tmp[18] = {
                    8.19e-03, 2.41e-03, 5.34e-03, 5.93e-03, 5.79e-03, 5.99e-03,
                    3.80e-03, 1.66e-03, 1.71e-03, 4.59e-03, 2.41e-03, 3.31e-03,
                    1.16e-02, 2.03e-03, 6.59e-03, 1.48e-02, 5.55e-03, 1.28e-02};
                double cutdphi_tmp[18] = {
                    1.72e-02, 1.02e-02, 1.57e-02, 3.54e-02, 8.80e-03, 1.14e-02,
                    3.97e-02, 5.64e-03, 1.32e-02, 1.08e-02, 8.53e-03, 9.25e-03,
                    7.37e-02, 5.66e-02, 3.59e-02, 1.87e-02, 1.20e-02, 3.58e-02};
                double cuteopin_tmp[18] = {
                    8.91e-01, 9.84e-01, 9.22e-01, 9.61e-01, 9.53e-01, 7.97e-01,
                    9.51e-01, 1.01e+00, 1.12e+00, 1.26e+00, 1.03e+00, 1.12e+00,
                    8.47e-01, 9.53e-01, 9.79e-01, 8.41e-01, 7.71e-01, 1.09e+00};
                double cutet_tmp[18] = {
                    0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
                    0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
                    1.46e+01, 1.44e+01, 1.61e+01, 1.52e+01, 1.53e+01, 1.61e+01};
                double cuthoe_tmp[18] = {
                    2.49e-02, 1.63e-02, 1.94e-02, 1.72e-02, 1.26e-02, 1.12e-02,
                    5.70e-02, 7.18e-02, 4.03e-02, 8.01e-02, 2.77e-03, 4.80e-02,
                    3.82e-02, 9.15e-02, 4.51e-02, 4.52e-02, 1.96e-03, 4.30e-03};
                double cutip_tmp[18] = {
                    1.91e-02, 1.03e-02, 2.19e-02, 1.44e-02, 1.36e-02, 1.25e-02,
                    7.44e-03, 7.45e-03, 8.02e-03, 5.16e-03, 5.57e-03, 3.86e-03,
                    2.10e+00, 5.27e-03, 3.17e+00, 4.91e+00, 7.69e-01, 5.90e+00};
                double cutisoecal_tmp[18] = {
                    4.03e+00, 9.81e+00, 3.34e+00, 3.74e+00, 2.05e+00, 2.30e+00,
                    4.60e+00, 4.69e+00, 4.80e+00, 3.65e+00, 2.65e+00, 3.96e+00,
                    1.01e+01, 1.24e+01, 1.11e+01, 1.10e+01, 1.06e+01, 1.34e+01};
                double cutisohcal_tmp[18] = {
                    5.00e+00, 3.99e+00, 6.63e+00, 1.94e+00, 2.62e+00, 4.54e+00,
                    2.59e+00, 3.93e+00, 1.21e+00, 1.44e+00, 3.90e+00, 3.83e-01,
                    1.64e-01, 5.46e+00, 1.20e+01, 6.04e-03, 4.10e+00, 6.28e-04};
                double cutisotk_tmp[18] = {
                    2.66e+00, 1.30e+00, 1.07e+00, 2.39e+00, 9.35e-01, 8.17e-01,
                    2.07e+00, 2.89e+00, 1.92e+00, 6.11e-01, 2.58e-02, 3.45e-03,
                    5.29e+00, 5.18e+00, 1.54e+01, 5.38e+00, 4.47e+00, 3.47e-02};
                double cutmishits_tmp[18] = {
                    5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01,
                    5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01,
                    5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
                double cutsee_tmp[18] = {
                    1.07e-02, 1.01e-02, 1.02e-02, 3.13e-02, 2.72e-02, 2.67e-02,
                    9.46e-03, 9.61e-03, 9.85e-03, 2.69e-02, 2.55e-02, 2.53e-02,
                    1.42e-02, 1.06e-02, 1.03e-02, 3.50e-02, 2.96e-02, 3.33e-02};
                eidAssign(cutdeta, cutdeta_tmp, 18);
                eidAssign(cutdphi, cutdphi_tmp, 18);
                eidAssign(cuteopin, cuteopin_tmp, 18);
                eidAssign(cutet, cutet_tmp, 18);
                eidAssign(cuthoe, cuthoe_tmp, 18);
                eidAssign(cutip, cutip_tmp, 18);
                eidAssign(cutisoecal, cutisoecal_tmp, 18);
                eidAssign(cutisohcal, cutisohcal_tmp, 18);
                eidAssign(cutisotk, cutisotk_tmp, 18);
                eidAssign(cutmishits, cutmishits_tmp, 18);
                eidAssign(cutsee, cutsee_tmp, 18);
                return;
            }
        case CIC_HYPERTIGHT4:
            {   
                double cutdeta_tmp[18] = {
                    8.22e-03, 2.04e-03, 8.82e-03, 4.98e-03, 5.35e-03, 4.54e-03,
                    3.41e-03, 3.92e-03, 1.64e-03, 1.00e-02, 3.32e-03, 6.44e-03,
                    1.16e-02, 2.03e-03, 6.59e-03, 1.48e-02, 5.55e-03, 1.28e-02};
                double cutdphi_tmp[18] = {
                    1.82e-02, 5.88e-03, 1.11e-02, 1.61e-02, 5.32e-03, 9.21e-03,
                    8.33e-03, 3.82e-03, 1.61e-02, 7.76e-03, 7.18e-03, 9.03e-03,
                    7.37e-02, 5.66e-02, 3.59e-02, 1.87e-02, 1.20e-02, 3.58e-02};
                double cuteopin_tmp[18] = {
                    9.37e-01, 9.84e-01, 9.45e-01, 9.65e-01, 9.67e-01, 1.13e+00,
                    8.97e-01, 1.01e+00, 1.26e+00, 1.29e+00, 1.03e+00, 1.33e+00,
                    8.47e-01, 9.53e-01, 9.79e-01, 8.41e-01, 7.71e-01, 1.09e+00};
                double cutet_tmp[18] = {
                    0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
                    0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
                    1.46e+01, 1.44e+01, 1.61e+01, 1.52e+01, 1.53e+01, 1.61e+01};
                double cuthoe_tmp[18] = {
                    2.54e-02, 2.34e-02, 1.60e-03, 2.63e-02, 1.54e-02, 1.80e-02,
                    1.47e-02, 5.92e-02, 6.22e-02, 1.18e-02, 1.29e-04, 4.94e-02,
                    3.82e-02, 9.15e-02, 4.51e-02, 4.52e-02, 1.96e-03, 4.30e-03};
                double cutip_tmp[18] = {
                    1.72e-02, 1.38e-02, 1.45e-02, 1.07e-02, 5.34e-02, 1.01e-02,
                    5.37e-03, 5.31e-03, 3.90e-03, 4.58e-03, 4.94e-03, 3.90e-03,
                    2.10e+00, 5.27e-03, 3.17e+00, 4.91e+00, 7.69e-01, 5.90e+00};
                double cutisoecal_tmp[18] = {
                    3.49e+00, 2.68e+00, 2.44e+00, 8.86e+00, 1.01e+00, 1.87e+00,
                    4.46e+00, 4.85e+00, 4.74e+00, 3.29e+00, 2.50e+00, 4.44e+00,
                    1.01e+01, 1.24e+01, 1.11e+01, 1.10e+01, 1.06e+01, 1.34e+01};
                double cutisohcal_tmp[18] = {
                    1.56e+00, 6.13e+00, 6.24e+00, 1.03e+00, 2.49e+00, 6.78e-01,
                    2.66e+00, 2.84e+00, 2.55e+00, 6.70e-02, 4.76e+00, 4.80e-01,
                    1.64e-01, 5.46e+00, 1.20e+01, 6.04e-03, 4.10e+00, 6.28e-04};
                double cutisotk_tmp[18] = {
                    2.00e+00, 7.84e-01, 7.16e-01, 1.53e+00, 9.23e-01, 2.91e-01,
                    4.32e-01, 7.54e-01, 2.05e+00, 2.85e-02, 1.20e-03, 1.61e-04,
                    5.29e+00, 5.18e+00, 1.54e+01, 5.38e+00, 4.47e+00, 3.47e-02};
                double cutmishits_tmp[18] = {
                    5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01,
                    5.00e-01, 5.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01,
                    5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
                double cutsee_tmp[18] = {
                    1.04e-02, 1.06e-02, 9.82e-03, 3.29e-02, 2.69e-02, 2.65e-02,
                    9.21e-03, 9.35e-03, 9.99e-03, 2.62e-02, 2.55e-02, 2.52e-02,
                    1.42e-02, 1.06e-02, 1.03e-02, 3.50e-02, 2.96e-02, 3.33e-02};
                eidAssign(cutdeta, cutdeta_tmp, 18);
                eidAssign(cutdphi, cutdphi_tmp, 18);
                eidAssign(cuteopin, cuteopin_tmp, 18);
                eidAssign(cutet, cutet_tmp, 18);
                eidAssign(cuthoe, cuthoe_tmp, 18);
                eidAssign(cutip, cutip_tmp, 18);
                eidAssign(cutisoecal, cutisoecal_tmp, 18);
                eidAssign(cutisohcal, cutisohcal_tmp, 18);
                eidAssign(cutisotk, cutisotk_tmp, 18);
                eidAssign(cutmishits, cutmishits_tmp, 18);
                eidAssign(cutsee, cutsee_tmp, 18);
                return;
            }
        default:
            std::cout << "[eidGetCIC] ERROR! Invalid tightness level" << std::endl;
    }

    return;

}

void eidAssign(std::vector<double> &cutarr, double cutvals[], unsigned int size)
{
    cutarr.clear();
    for (unsigned int i = 0; i < size; ++i) {
        cutarr.push_back(cutvals[i]);
    }
}


void eidAssign(std::vector<bool> &cutarr, bool cutvals[], unsigned int size)
{
    cutarr.clear();
    for (unsigned int i = 0; i < size; ++i) {
        cutarr.push_back(cutvals[i]);
    }
}

void eidAssign(std::vector<int> &cutarr, int cutvals[], unsigned int size)
{
    cutarr.clear();
    for (unsigned int i = 0; i < size; ++i) {
        cutarr.push_back(cutvals[i]);
    }
}


