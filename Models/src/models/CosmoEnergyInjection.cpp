///  GAMBIT: Global and Modular BSM Inference Tool
///  *********************************************
///
///  Model translation functions for Cosmological Energy injection
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
///  \date 2019 Sep
///
///  *********************************************

#include <string>
#include <vector>

#include "gambit/Models/model_macros.hpp"
#include "gambit/Models/model_helpers.hpp"
#include "gambit/Logs/logger.hpp"
#include "gambit/Utils/util_functions.hpp"
#include "gambit/Utils/numerical_constants.hpp"

#include "gambit/Models/models/CosmoEnergyInjection.hpp"

#define MODEL DecayingDM_photon
void MODEL_NAMESPACE::DecayingDM_photon_to_DecayingDM_general (const ModelParameters &myP, ModelParameters &targetP)
{
  USE_MODEL_PIPE(DecayingDM_general)
  logger()<<"Running interpret_as_parent calculations for DecayingDM_photon --> DecayingDM_general ..."<<LogTags::info<<EOM;

  targetP.setValue("mass", myP.getValue("mass"));
  targetP.setValue("lifetime", myP.getValue("lifetime"));
  targetP.setValue("fraction", myP.getValue("fraction"));

  // All decays into photons (BR into electrons is 1.0)
  targetP.setValue("BR", 1.0);
}
#undef MODEL

#define MODEL DecayingDM_electron
void MODEL_NAMESPACE::DecayingDM_electron_to_DecayingDM_general (const ModelParameters &myP, ModelParameters &targetP)
{
  USE_MODEL_PIPE(DecayingDM_general)
  logger()<<"Running interpret_as_parent calculations for DecayingDM_electron --> DecayingDM_general ..."<<LogTags::info<<EOM;

  targetP.setValue("mass", myP.getValue("mass"));
  targetP.setValue("lifetime", myP.getValue("lifetime"));
  targetP.setValue("fraction", myP.getValue("fraction"));

  // All decays into electrons
  targetP.setValue("BR", 0.0);
}
#undef MODEL
