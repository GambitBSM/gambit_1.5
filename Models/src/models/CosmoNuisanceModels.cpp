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


#define MODEL cosmo_nuisance_Pantheon
// Translation function definition
void MODEL_NAMESPACE::cosmo_nuisance_Pantheon_to_cosmo_nuisance_JLA (const ModelParameters &myP, ModelParameters &targetP)
{
  logger()<<"Running interpret_as_parent calculations for cosmo_nuisance_Pantheon --> cosmo_nuisance_JLA ..."<<LogTags::info<<EOM;

  targetP.setValue("M", myP.getValue("M"));
  // the rest (alpha, beta & Delta_M) automatically defaults to 0
}
#undef MODEL

#define MODEL cosmo_nuisance_BK14priors
// Translation function definition
void MODEL_NAMESPACE::cosmo_nuisance_BK14priors_to_cosmo_nuisance_BK14 (const ModelParameters &myP, ModelParameters &targetP)
{
  logger()<<"Running interpret_as_parent calculations for cosmo_nuisance_BK14priors --> cosmo_nuisance_BK14 ..."<<LogTags::info<<EOM;

  targetP.setValue("BBbetadust", myP.getValue("BBbetadust"));
  targetP.setValue("BBbetasync", myP.getValue("BBbetasync"));
  // the rest BBdust,BBsync,BBalphadust,BBTdust,BBalphasync,BBdustsynccorr,EEtoBB_dust,EEtoBB_sync automatically defaults to 0
}
#undef MODEL

