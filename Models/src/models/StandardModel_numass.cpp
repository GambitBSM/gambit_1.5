///  GAMBIT: Global and Modular BSM Inference Tool
///  *********************************************
///
///  Model translation functions for neutrino
///  masses in cosmology
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
///  \date 2019 Feb, Jul
///
///  *********************************************

#include <string>
#include <vector>

#include "gambit/Models/model_macros.hpp"
#include "gambit/Models/model_helpers.hpp"
#include "gambit/Logs/logger.hpp"
#include "gambit/Utils/util_functions.hpp"
#include "gambit/Utils/numerical_constants.hpp"

#include "gambit/Models/models/StandardModel_numass.hpp"

/////////////// translation functions ///////////////////////////

#define MODEL StandardModel_numass_single
#define PARENT StandardModel_SLHA2

// Translation function definition
void MODEL_NAMESPACE::StandardModel_numass_single_to_StandardModel_SLHA2 (const ModelParameters &myP, ModelParameters &targetP)
{
  USE_MODEL_PIPE(PARENT) // get pipe for "interpret as PARENT" function
  logger()<<"Running interpret_as_parent calculations for StandardModel_numass_single --> StandardModel_SLHA2..."<<LogTags::info<<EOM;

  targetP.setValues( myP, false);

  targetP.setValue("mNu1", 1e-9*myP.getValue("mNu") );
  targetP.setValue("mNu2", 0.0 );
  targetP.setValue("mNu3", 0.0 );
}

#undef PARENT
#undef MODEL

#define MODEL StandardModel_numass_degenerate
#define PARENT StandardModel_SLHA2

// Translation function definition
void MODEL_NAMESPACE::StandardModel_numass_degenerate_to_StandardModel_SLHA2 (const ModelParameters &myP, ModelParameters &targetP)
{
  USE_MODEL_PIPE(PARENT) // get pipe for "interpret as PARENT" function
  logger()<<"Running interpret_as_parent calculations for StandardModel_numass_degenerate --> StandardModel_SLHA2..."<<LogTags::info<<EOM;

  targetP.setValues( myP, false);

  double mNu = 1e-9*myP.getValue("Smu")/3.;

  targetP.setValue("mNu1", mNu );
  targetP.setValue("mNu2", mNu );
  targetP.setValue("mNu3", mNu );
}

#undef PARENT
#undef MODEL
