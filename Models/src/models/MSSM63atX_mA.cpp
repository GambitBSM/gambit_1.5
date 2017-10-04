///  GAMBIT: Global and Modular BSM Inference Tool
///  *********************************************
///
///  Translation function definitions for all
///  'mA' versions of 63 parameter MSSM, back to
///  the corresponding mhu2 mhd2 parameterisations
///
///  Contains translation functions for:
///   MSSM63atQ_mA     --> MSSM63atQ
///   MSSM63atMGUT_mA  --> MSSM63atMGUT
///   MSSM63atMSUSY_mA --> MSSM63atMSUSY
///
///  *********************************************
///
///  Authors
///  =======
///
///  (add name and date if you modify)
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2017 Aug
///
///  \author Ben Farmer
///          (benjamin.farmer@fysik.su.se)
///  \date 2017 Sep, Oct
///
///  *********************************************

#include <string>
#include <vector>

#include "gambit/Models/model_macros.hpp"
#include "gambit/Models/model_helpers.hpp"
#include "gambit/Logs/logger.hpp"
#include "gambit/Utils/util_functions.hpp"

#include "gambit/Models/models/MSSM63atQ_mA.hpp"
#include "gambit/Models/models/MSSM63atMSUSY_mA.hpp"
#include "gambit/Models/models/MSSM63atMGUT_mA.hpp"
#include "gambit/Elements/spectrum.hpp"

// Activate debug output
//#define MSSM63_mA_DBUG

// General helper translation function definition
namespace Gambit { 
  void MSSM_mA_to_MSSM_mhud(const ModelParameters &myP, ModelParameters &targetP, const SubSpectrum& HE)
  {
     // Copy all the common parameters of MSSM63at<X>_mA into MSSM63at<X>
     targetP.setValues(myP,false); // Set "missing_is_error" flag to false since e.g. mA parameter from MSSM63atQ_mA does not exist in MSSM63atQ. Similar for variants at other scales.
  
     // Set the sign of mu
     targetP.setValue("SignMu", Gambit::sgn(myP["mu"]));
  
     // Now only the "mHu2" and "mHd2" parameters are left unset. Extract these from the Spectrum object dependency.
     // Make sure the high-scale value was correctly added to the spectrum wrapper object
     if (HE.has(Par::mass2,"mHu2") and HE.has(Par::mass2,"mHd2"))
     {
        targetP.setValue("mHu2", HE.get(Par::mass2,"mHu2"));
        targetP.setValue("mHd2", HE.get(Par::mass2,"mHd2"));
     }
     else
     {
        model_error().raise(LOCAL_INFO,"Parameter with name 'mHu2' or 'mHd2' (type Par::mass2) not found in Spectrum object! "
                                       "Translation from MSSM<X>_mA to MSSM<X> is not possible without this value. "
                                       "Please use a Spectrum wrapper that provides it.");
     }
  
  }
}
/// @{ Actual translation function definitions

// Need to define MODEL and PARENT in order for helper macros to work correctly
#define MODEL  MSSM63atQ_mA
#define PARENT MSSM63atQ
void MODEL_NAMESPACE::MSSM63atQ_mA_to_MSSM63atQ (const ModelParameters &myP, ModelParameters &targetP)
{
   USE_MODEL_PIPE(PARENT)
   logger()<<"Running interpret_as_parent calculations for MSSM63atQ_mA --> MSSM63atQ..."<<LogTags::info<<EOM;
   const SubSpectrum& HE = Dep::unimproved_MSSM_spectrum->get_HE();
   MSSM_mA_to_MSSM_mhud(myP, targetP, HE);
}
#undef PARENT
#undef MODEL

#define MODEL  MSSM63atMSUSY_mA
#define PARENT MSSM63atMSUSY
void MODEL_NAMESPACE::MSSM63atMSUSY_mA_to_MSSM63atMSUSY (const ModelParameters &myP, ModelParameters &targetP)
{
   USE_MODEL_PIPE(PARENT)
   logger()<<"Running interpret_as_parent calculations for MSSM63atMSUSY_mA --> MSSM63atMSUSY..."<<LogTags::info<<EOM;
   const SubSpectrum& HE = Dep::unimproved_MSSM_spectrum->get_HE();
   MSSM_mA_to_MSSM_mhud(myP, targetP, HE);
}
#undef PARENT
#undef MODEL

#define MODEL  MSSM63atMGUT_mA
#define PARENT MSSM63atMGUT
void MODEL_NAMESPACE::MSSM63atMGUT_mA_to_MSSM63atMGUT (const ModelParameters &myP, ModelParameters &targetP)
{
   USE_MODEL_PIPE(PARENT)
   logger()<<"Running interpret_as_parent calculations for MSSM63atMGUT_mA --> MSSM63atMGUT..."<<LogTags::info<<EOM;
   const SubSpectrum& HE = Dep::unimproved_MSSM_spectrum->get_HE();
   MSSM_mA_to_MSSM_mhud(myP, targetP, HE);
}
#undef PARENT
#undef MODEL
