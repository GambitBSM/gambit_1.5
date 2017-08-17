///  GAMBIT: Global and Modular BSM Inference Tool
///  *********************************************
///
///  MSSM63atQ_mA translation function definitions
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
///  *********************************************

#include <string>
#include <vector>

#include "gambit/Models/model_macros.hpp"
#include "gambit/Models/model_helpers.hpp"
#include "gambit/Logs/logger.hpp"
#include "gambit/Utils/util_functions.hpp"

#include "gambit/Models/models/MSSM63atQ_mA.hpp"
#include "gambit/Elements/spectrum.hpp"

// Activate debug output
//#define MSSM63atQ_mA_DBUG

// Need to define MODEL and PARENT in order for helper macros to work correctly
#define MODEL  MSSM63atQ_mA
#define PARENT MSSM63atQ

// Translation function definition
void MODEL_NAMESPACE::MSSM63atQ_mA_to_MSSM63atQ (const ModelParameters &myP, ModelParameters &targetP)
{
   USE_MODEL_PIPE(PARENT) // get pipe for "interpret as PARENT" function
   logger()<<"Running interpret_as_parent calculations for MSSM63atQ_mA --> MSSM63atQ..."<<LogTags::info<<EOM;

   // Copy all the common parameters of MSSM63atMGUT into MSSM63atQ
   targetP.setValues(myP);

   // Set the sign of mu
   targetP.setValue("SignMu", Gambit::sgn(myP["mu"]));

   // Now only the "mHu2" and "mHd2" parameters are left unset. Extract these from the Spectrum object dependency.
   const SubSpectrum& HE = Dep::unimproved_MSSM_spectrum->get_HE();

   // Make sure the high-scale value was correctly added to the spectrum wrapper object
   if (HE.has(Par::mass2,"mHu2") and HE.has(Par::mass2,"mHd2"))
   {
      targetP.setValue("mHu2", HE.get(Par::mass2,"mHu2"));
      targetP.setValue("mHd2", HE.get(Par::mass2,"mHd2"));
   }
   else
   {
      model_error().raise(LOCAL_INFO,"Parameter with name 'mHu2' or 'mHd2' (type Par::mass2) not found in Spectrum object! "
                                     "Translation from MSSM63atQ_mA to MSSM63atQ is not possible without this value. "
                                     "Please use a Spectrum wrapper that provides it.");
   }

   // Done! Check that everything is ok if desired.
   #ifdef MSSM63atQ_mA_DBUG
     std::cout << "MSSM63atQ_mA parameters:" << myP << std::endl;
     std::cout << "MSSM63atQ parameters   :" << targetP << std::endl;
   #endif
}

#undef PARENT
#undef MODEL
