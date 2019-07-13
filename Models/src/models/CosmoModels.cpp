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
///  \date 2019 Feb, Jun
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///   \date 2019 Feb
///
///  *********************************************

#include <string>
#include <vector>

#include "gambit/Models/model_macros.hpp"
#include "gambit/Models/model_helpers.hpp"
#include "gambit/Logs/logger.hpp"
#include "gambit/Utils/util_functions.hpp"
#include "gambit/Utils/numerical_constants.hpp"

#include "gambit/Models/models/CosmoModels.hpp"

/////////////// translation functions for Cosmology models ///////////////////////////

#define MODEL etaBBN
void MODEL_NAMESPACE::etaBBN_to_etaBBN_rBBN_rCMB_dNeffBBN_dNeffCMB (const ModelParameters &myP, ModelParameters &targetP)
{
  logger()<<"Running interpret_as_parent calculations for etaBBN --> etaBBN_rBBN_rCMB_dNeffBBN_dNeffCMB ..."<<LogTags::info<<EOM;

  // Set eta_BBN
  targetP.setValue("eta_BBN", myP.getValue("eta_BBN"));

  // No dNeff due to an altered neutrino temperature
  targetP.setValue("r_BBN", 1.);
  targetP.setValue("r_CMB", 1.);

  // No dNeff due to additional radiation
  targetP.setValue("dNeff_BBN", 0.);
  targetP.setValue("dNeff_CMB", 0.);
}
#undef MODEL

#define MODEL rBBN_rCMB_dNeffBBN_dNeffCMB
#define PARENT etaBBN_rBBN_rCMB_dNeffBBN_dNeffCMB
void MODEL_NAMESPACE::rBBN_rCMB_dNeffBBN_dNeffCMB_to_etaBBN_rBBN_rCMB_dNeffBBN_dNeffCMB (const ModelParameters &myP, ModelParameters &targetP)
{
  USE_MODEL_PIPE(PARENT) // get pipe for "INTERPRET_AS_PARENT_DEPENDENCY"
  logger()<<"Running interpret_as_parent calculations for rBBN_rCMB_dNeffBBN_dNeffCMB --> etaBBN_rBBN_rCMB_dNeffBBN_dNeffCMB ..."<<LogTags::info<<EOM;

  // Set eta_BBN (This requires that eta0 / omega_b is known, i.e. LCDM is in use)
  // -- The dependency_resolver will figure that out --
  targetP.setValue("eta_BBN", *Dep::eta0);

  // Set the respective values of r
  targetP.setValue("r_BBN", myP.getValue("r_BBN"));
  targetP.setValue("r_CMB", myP.getValue("r_CMB"));

  // Set the respective values of dNeff
  targetP.setValue("dNeff_BBN", myP.getValue("dNeff_BBN"));
  targetP.setValue("dNeff_CMB", myP.getValue("dNeff_CMB"));
}
#undef PARENT
#undef MODEL

#define MODEL rBBN_rCMB
void MODEL_NAMESPACE::rBBN_rCMB_to_rBBN_rCMB_dNeffBBN_dNeffCMB (const ModelParameters &myP, ModelParameters &targetP)
{
  logger()<<"Running interpret_as_parent calculations for rBBN_rCMB --> rBBN_rCMB_dNeffBBN_dNeffCMB ..."<<LogTags::info<<EOM;

  // Set the respective values of r
  targetP.setValue("r_BBN", myP.getValue("r_BBN"));
  targetP.setValue("r_CMB", myP.getValue("r_CMB"));

  // No dNeff due to additional radiation
  targetP.setValue("dNeff_BBN", 0.);
  targetP.setValue("dNeff_CMB", 0.);
}
#undef MODEL

#define MODEL rCMB
void MODEL_NAMESPACE::rCMB_to_rBBN_rCMB (const ModelParameters &myP, ModelParameters &targetP)
{
  logger()<<"Running interpret_as_parent calculations for rCMB --> rBBN_rCMB ..."<<LogTags::info<<EOM;

  // Set the respective values of r
  targetP.setValue("r_BBN", myP.getValue("r_CMB"));
  targetP.setValue("r_CMB", myP.getValue("r_CMB"));
}
#undef MODEL

#define MODEL dNeffBBN_dNeffCMB
void MODEL_NAMESPACE::dNeffBBN_dNeffCMB_to_rBBN_rCMB_dNeffBBN_dNeffCMB (const ModelParameters &myP, ModelParameters &targetP)
{
  logger()<<"Running interpret_as_parent calculations for dNeffBBN_dNeffCMB --> rBBN_rCMB_dNeffBBN_dNeffCMB ..."<<LogTags::info<<EOM;

  // No dNeff due to an altered neutrino temperature
  targetP.setValue("r_BBN", 1.);
  targetP.setValue("r_CMB", 1.);

  // Set the respective values of dNeff
  targetP.setValue("dNeff_BBN", myP.getValue("dNeff_BBN"));
  targetP.setValue("dNeff_CMB", myP.getValue("dNeff_CMB"));
}
#undef MODEL

#define MODEL dNeffCMB
void MODEL_NAMESPACE::dNeffCMB_to_dNeffBBN_dNeffCMB (const ModelParameters &myP, ModelParameters &targetP)
{
  logger()<<"Running interpret_as_parent calculations for dNeffCMB --> dNeffBBN_dNeffCMB ..."<<LogTags::info<<EOM;

  // Set the respective values of dNeff
  targetP.setValue("dNeff_BBN", myP.getValue("dNeff_CMB"));
  targetP.setValue("dNeff_CMB", myP.getValue("dNeff_CMB"));
}
#undef MODEL
