///  GAMBIT: Global and Modular BSM Inference Tool
///  *********************************************
///
///  Model translation functions for cosmology models
///
///  *********************************************
///
///  Authors
///  =======
///
///  (add name and date if you modify)
///
///  \author Patrick Stöcker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2019 Feb
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///   \date 2019 Feb, Jun
///
///  *********************************************

#include <string>
#include <vector>

#include "gambit/Models/model_macros.hpp"
#include "gambit/Models/model_helpers.hpp"
#include "gambit/Logs/logger.hpp"
#include "gambit/Utils/util_functions.hpp"

#include "gambit/Models/models/CosmoNuisanceModels.hpp"

/////////////// translation functions for cosmological nuisance parameter models ///////////////////////////


#define MODEL cosmo_nuisance_params_Pantheon
  #define PARENT cosmo_nuisance_params_JLA
    // Translation function definition
    void MODEL_NAMESPACE::cosmo_nuisance_params_Pantheon_to_cosmo_nuisance_params_JLA (const ModelParameters &myP, ModelParameters &targetP)
    {
      USE_MODEL_PIPE(PARENT) // get pipe for "interpret as PARENT" function
      logger()<<"Running interpret_as_parent calculations for cosmo_nuisance_params_Pantheon --> cosmo_nuisance_params_JLA ..."<<LogTags::info<<EOM;
    
      targetP.setValue("M", myP.getValue("M"));
      // the rest (alpha, beta & Delta_M) automatically defaults to 0
    }
  #undef PARENT
#undef MODEL

