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
#include "gambit/Utils/numerical_constants.hpp"

#include "gambit/Models/models/CosmoModels.hpp"

/////////////// translation functions for Cosmology models ///////////////////////////

#define MODEL  LCDM_dNeffCMB_dNeffBBN
#define PARENT LCDM_dNeffCMB_dNeffBBN_etaBBN

// Translation function definition
void MODEL_NAMESPACE::LCDM_dNeffCMB_dNeffBBN_to_LCDM_dNeffCMB_dNeffBBN_etaBBN (const ModelParameters &myP, ModelParameters &targetP)
{
  USE_MODEL_PIPE(PARENT) // get pipe for "interpret as PARENT" function
  logger()<<"Running interpret_as_parent calculations for LCDM_dNeffCMB_dNeffBBN --> LCDM_dNeffCMB_dNeffBBN_etaBBN ..."<<LogTags::info<<EOM;

  targetP.setValue("omega_b", myP.getValue("omega_b"));
  targetP.setValue("omega_cdm",myP.getValue("omega_cdm"));
  targetP.setValue("H0", myP.getValue("H0") );
  targetP.setValue("ln10A_s", myP.getValue("ln10A_s") );
  targetP.setValue("n_s", myP.getValue("n_s") );
  targetP.setValue("tau_reio", myP.getValue("tau_reio") );
  targetP.setValue("dNeff", myP.getValue("dNeff") );
  targetP.setValue("dNeff_BBN", myP.getValue("dNeff_BBN") );

  targetP.setValue("eta_BBN", *Dep::etaBBN );
}

#undef PARENT
#undef MODEL

#define MODEL  LCDM_dNeffCMB
#define PARENT LCDM_dNeffCMB_dNeffBBN

// Translation function definition
void MODEL_NAMESPACE::LCDM_dNeffCMB_to_LCDM_dNeffCMB_dNeffBBN (const ModelParameters &myP, ModelParameters &targetP)
{
  USE_MODEL_PIPE(PARENT) // get pipe for "interpret as PARENT" function
  logger()<<"Running interpret_as_parent calculations for LCDM_dNeffCMB --> LCDM_dNeffCMB_dNeffBBN ..."<<LogTags::info<<EOM;

  targetP.setValue("omega_b", myP.getValue("omega_b"));
  targetP.setValue("omega_cdm",myP.getValue("omega_cdm"));
  targetP.setValue("H0", myP.getValue("H0") );
  targetP.setValue("ln10A_s", myP.getValue("ln10A_s") );
  targetP.setValue("n_s", myP.getValue("n_s") );
  targetP.setValue("tau_reio", myP.getValue("tau_reio") );
  targetP.setValue("dNeff", myP.getValue("dNeff") );
  targetP.setValue("dNeff_BBN", myP.getValue("dNeff") );
}

#undef PARENT
#undef MODEL

#define MODEL  LCDM_ExtdNeffCMB_ExtetaBBN
#define PARENT LCDM_dNeffCMB_dNeffBBN_etaBBN

// Translation function definition
void MODEL_NAMESPACE::LCDM_ExtdNeffCMB_ExtetaBBN_to_LCDM_dNeffCMB_dNeffBBN_etaBBN (const ModelParameters &myP, ModelParameters &targetP)
{
  USE_MODEL_PIPE(PARENT) // get pipe for "interpret as PARENT" function
  logger()<<"Running interpret_as_parent calculations for LCDM_ExtdNeffCMB_ExtetaBBN --> LCDM_dNeffCMB_dNeffBBN_etaBBN ..."<<LogTags::info<<EOM;

  targetP.setValue("omega_b", myP.getValue("omega_b"));
  targetP.setValue("omega_cdm",myP.getValue("omega_cdm"));
  targetP.setValue("H0", myP.getValue("H0") );
  targetP.setValue("ln10A_s", myP.getValue("ln10A_s") );
  targetP.setValue("n_s", myP.getValue("n_s") );
  targetP.setValue("tau_reio", myP.getValue("tau_reio") );

  targetP.setValue("eta_BBN", *Dep::etaBBN);
  targetP.setValue("dNeff_BBN", 0. );
  targetP.setValue("dNeff", *Dep::ExtdNeffCMB );
}

#undef PARENT
#undef MODEL

#define MODEL  LCDM
#define PARENT LCDM_dNeffCMB

// Translation function definition
void MODEL_NAMESPACE::LCDM_to_LCDM_dNeffCMB (const ModelParameters &myP, ModelParameters &targetP)
{
  USE_MODEL_PIPE(PARENT) // get pipe for "interpret as PARENT" function
  logger()<<"Running interpret_as_parent calculations for LCDM --> LCDM_dNeffCMB ..."<<LogTags::info<<EOM;

  targetP.setValue("omega_b", myP.getValue("omega_b"));
  targetP.setValue("omega_cdm",myP.getValue("omega_cdm"));
  targetP.setValue("H0", myP.getValue("H0") );
  targetP.setValue("ln10A_s", myP.getValue("ln10A_s") );
  targetP.setValue("n_s", myP.getValue("n_s") );
  targetP.setValue("tau_reio", myP.getValue("tau_reio") );
  targetP.setValue("dNeff", 0. );
}

#undef PARENT
#undef MODEL


#define MODEL  cosmo_nuisance_params_Pantheon
#define PARENT cosmo_nuisance_params_JLA

// Translation function definition
void MODEL_NAMESPACE::cosmo_nuisance_params_Pantheon_to_cosmo_nuisance_params_JLA (const ModelParameters &myP, ModelParameters &targetP)
{
  USE_MODEL_PIPE(PARENT) // get pipe for "interpret as PARENT" function
  logger()<<"Running interpret_as_parent calculations for cosmo_nuisance_params_Pantheon --> cosmo_nuisance_params_JLA ..."<<LogTags::info<<EOM;

  targetP.setValue("M", myP.getValue("M"));
  targetP.setValue("alpha",0.);
  targetP.setValue("beta", 0. );
  targetP.setValue("Delta_M", 0. );
}

#undef PARENT
#undef MODEL

