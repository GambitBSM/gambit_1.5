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

// (PS:) Please note that in this implementation we assume the parameter dmNu3l
// to be the one which has the higher absolute value in the given hierarchy.
// This agrres with the convention of NuFit (cf. 1811.05487)
//
// !-> This inverted to the definition of StandardModel_mNudiff <-!
#define MODEL StandardModel_numass_split
#define PARENT StandardModel_SLHA2

void MODEL_NAMESPACE::StandardModel_numass_split_to_StandardModel_SLHA2 (const ModelParameters &myP, ModelParameters &targetP)
{

  logger()<<"Running interpret_as_parent calculations for StandardModel_numass_split --> StandardModel_SLHA2."<<LogTags::info<<EOM;

  // Copy all the "non-neutrino" paramters of the SLHA2 block.
  targetP.setValues(myP,false);

  if (myP["dmNu3l"] > 0.)  // normal hierarchy, l = 1
  {
   targetP.setValue("mNu1", myP["mNu_light"]*1e-9);
   targetP.setValue("mNu2",
       pow(myP["mNu_light"]*myP["mNu_light"]+myP["dmNu21"], 0.5)*1e-9);
   targetP.setValue("mNu3",
       pow(myP["mNu_light"]*myP["mNu_light"]+myP["dmNu3l"], 0.5)*1e-9);
  }
  else // inverted hierarchy, l = 2
  {
   targetP.setValue("mNu3", myP["mNu_light"]*1e-9);
   targetP.setValue("mNu2",
       pow(myP["mNu_light"]*myP["mNu_light"]-myP["dmNu3l"], 0.5)*1e-9);
   targetP.setValue("mNu1",
       pow(myP["mNu_light"]*myP["mNu_light"]-myP["dmNu3l"]-myP["dmNu21"], 0.5)*1e-9);
  }
}

#undef PARENT
#undef MODEL
