///  GAMBIT: Global and Modular BSM Inference Tool
///  *********************************************
///
///  Translation function definitions for
///  63 parameter MSSM variants back to the 
///  'master' MSSM63atQ generic model
///  
///  Contains interpret-as-parent definitions for
///   MSSM63atMGUT
///   MSSM63atMSUSY
///  back to their common parent:
///   MSSM63atQ
///
///  *********************************************
///
///  Authors
///  =======
///
///  (add name and date if you modify)
///
///  \author Ben Farmer
///          (benjamin.farmer@fysik.su.se)
///  \date 2015 Aug, 2017 Oct
///
///  *********************************************

#include <string>
#include <vector>

#include "gambit/Models/model_macros.hpp"
#include "gambit/Models/model_helpers.hpp"
#include "gambit/Logs/logger.hpp"
#include "gambit/Utils/util_functions.hpp"

#include "gambit/Models/models/MSSM63atQ.hpp"
#include "gambit/Models/models/MSSM63atMGUT.hpp"
#include "gambit/Models/models/MSSM63atMSUSY.hpp"
#include "gambit/Elements/spectrum.hpp"

using namespace Gambit::Utils;

// General helper translation function definition
namespace Gambit { 
  void MSSMatX_to_MSSMatQ(const ModelParameters &myP, ModelParameters &targetP, const SubSpectrum& HE)
  {
    // Copy all the parameters of MSSM63atMGUT into MSSM63atQ
    targetP.setValues(myP);

    // Now only the "Qin" parameter is left unset. Need to extract this from the Spectrum object dependency.
    // Make sure the high-scale value was correctly added to the spectrum wrapper object
    if( HE.has(Par::mass1,"high_scale") )
    {
       targetP.setValue("Qin", HE.get(Par::mass1,"high_scale") );
    }
    else
    {
       model_error().raise(LOCAL_INFO,"Parameter with name 'high_scale' (type Par::mass1) not found in Spectrum object! Translation from MSSM63at<X> to MSSM63atQ is not possible without this value. Please use a Spectrum wrapper which provides it.");
    }
    // Done!
  }
}


/// @{ Translation function definitions
#define PARENT MSSM63atQ
#define MODEL  MSSM63atMGUT
void MODEL_NAMESPACE::MSSM63atMGUT_to_MSSM63atQ (const ModelParameters &myP, ModelParameters &targetP)
{
   USE_MODEL_PIPE(PARENT) // get pipe for "interpret as PARENT" function
   logger()<<"Running interpret_as_parent calculations for MSSM63atMGUT --> MSSM63atQ..."<<LogTags::info<<EOM;
   const SubSpectrum& HE = Dep::unimproved_MSSM_spectrum->get_HE();
   MSSMatX_to_MSSMatQ(myP, targetP, HE);
}
#undef MODEL

#define MODEL  MSSM63atMSUSY
void MODEL_NAMESPACE::MSSM63atMSUSY_to_MSSM63atQ (const ModelParameters &myP, ModelParameters &targetP)
{
   USE_MODEL_PIPE(PARENT) // get pipe for "interpret as PARENT" function
   logger()<<"Running interpret_as_parent calculations for MSSM63atMSUSY --> MSSM63atQ..."<<LogTags::info<<EOM;
   const SubSpectrum& HE = Dep::unimproved_MSSM_spectrum->get_HE();
   MSSMatX_to_MSSMatQ(myP, targetP, HE);
}
#undef MODEL
#undef PARENT
/// @}
