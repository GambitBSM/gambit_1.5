//==========================================================================
// This file has been automatically generated for Pythia 8
// MadGraph5_aMC@NLO v. 2.3.0, 2015-07-01
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef Pythia8_parameters_MC4BSM_2012_UFO_H
#define Pythia8_parameters_MC4BSM_2012_UFO_H

#include <complex> 

#include "Pythia8/ParticleData.h"
#include "Pythia8/StandardModel.h"
#include "Pythia8/SusyLesHouches.h"

using namespace std; 

namespace Pythia8 
{

class Parameters_MC4BSM_2012_UFO
{
  public:

    static Parameters_MC4BSM_2012_UFO * getInstance(); 

    // Model parameters independent of aS
    double mdl_Wev, mdl_Wuv, mdl_Wpe2, mdl_Wpe1, mdl_WH, mdl_WW, mdl_WZ,
        mdl_WT, mdl_Mev, mdl_Muv, mdl_MH, mdl_MZ, mdl_MB, mdl_MT, mdl_MTA,
        mdl_lam2p, mdl_lam1p, mdl_lam2, mdl_lam1, mdl_MM12, mdl_MM2, mdl_MM1,
        mdl_ymtau, mdl_ymt, mdl_ymb, mdl_Gf, aEWM1, ZERO, mdl_MM1__exp__2,
        mdl_MM2__exp__2, mdl_MM1__exp__4, mdl_MM12__exp__4, mdl_MM2__exp__4,
        mdl_sqrt__2, mdl_MPe1, mdl_MPe2, mdl_th, mdl_MZ__exp__2,
        mdl_MZ__exp__4, mdl_MH__exp__2, mdl_cos__th, mdl_sin__th, mdl_aEW,
        mdl_MW, mdl_sqrt__aEW, mdl_ee, mdl_MW__exp__2, mdl_sw2, mdl_cw,
        mdl_sqrt__sw2, mdl_sw, mdl_g1, mdl_gw, mdl_v, mdl_v__exp__2, mdl_lam,
        mdl_yb, mdl_yt, mdl_ytau, mdl_muH, mdl_ee__exp__2, mdl_gw__exp__2,
        mdl_cw__exp__2, mdl_sw__exp__2;
    std::complex<double> mdl_complexi; 
    // Model parameters dependent on aS
    double aS, mdl_sqrt__aS, G, mdl_G__exp__2; 
    // Model couplings independent of aS
    std::complex<double> GC_1, GC_2, GC_3, GC_4, GC_5, GC_6, GC_7, GC_8, GC_12,
        GC_13, GC_14, GC_15, GC_16, GC_17, GC_18, GC_19, GC_20, GC_21, GC_22,
        GC_23, GC_24, GC_25, GC_26, GC_27, GC_28, GC_29, GC_30, GC_31, GC_32,
        GC_33, GC_34, GC_35, GC_36, GC_37, GC_38, GC_39, GC_40, GC_41, GC_42,
        GC_43, GC_44, GC_45, GC_46, GC_47, GC_48, GC_49, GC_50, GC_51, GC_52,
        GC_53, GC_54, GC_55, GC_56, GC_57, GC_58, GC_59, GC_60, GC_61, GC_62,
        GC_63, GC_64, GC_65, GC_66, GC_67, GC_68, GC_69, GC_70, GC_71, GC_72,
        GC_73, GC_74, GC_75, GC_76, GC_77, GC_78;
    // Model couplings dependent on aS
    std::complex<double> GC_9, GC_11, GC_10; 

    // Set parameters that are unchanged during the run
    void setIndependentParameters(ParticleData * & pd, Couplings * & csm,
        SusyLesHouches * & slhaPtr);
    // Set couplings that are unchanged during the run
    void setIndependentCouplings(); 
    // Set parameters that are changed event by event
    void setDependentParameters(ParticleData * & pd, Couplings * & csm,
        SusyLesHouches * & slhaPtr, double alpS);
    // Set couplings that are changed event by event
    void setDependentCouplings(); 

    // Print parameters that are unchanged during the run
    void printIndependentParameters(); 
    // Print couplings that are unchanged during the run
    void printIndependentCouplings(); 
    // Print parameters that are changed event by event
    void printDependentParameters(); 
    // Print couplings that are changed event by event
    void printDependentCouplings(); 


  private:
    static Parameters_MC4BSM_2012_UFO * instance; 
}; 

}  // end namespace Pythia8
#endif  // Pythia8_parameters_MC4BSM_2012_UFO_H

