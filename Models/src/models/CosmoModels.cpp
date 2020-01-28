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
///   \date 2019 Feb, Jun
///
///  \author Sanjay Bloor
///          (sanjay.bloor12@imperial.ac.uk)
///   \date 2019 Nov
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

// ΛCDM parameters without those relating to the primordial power spectrum (A_s, n_s)
// This model should be scanned alongside an inflationary model able to provide
// a primordial power spectrum. 
#define MODEL LCDM_no_primordial
void MODEL_NAMESPACE::LCDM_to_LCDM_no_primordial(const ModelParameters &myP, ModelParameters &targetP)
{
  logger()<<"Running interpret_as_parent calculations for LCDM --> LCDM_no_primordial ..." << LogTags::info << EOM;

  // Same non-primordial parameters as vanilla ΛCDM.
  targetP.setValue("H0",        myP.getValue("H0"));
  targetP.setValue("omega_b",   myP.getValue("omega_b"));
  targetP.setValue("omega_cdm", myP.getValue("omega_cdm"));
  targetP.setValue("tau_reio",  myP.getValue("tau_reio"));

  // Don't take any of the model parameters
  targetP.setValue("n_s", 0.);
  targetP.setValue("ln10A_s", 0.); 
  // *Technically* log of 0 is undefined but I think this is okay for the translation functions...
}
#undef MODEL

#define MODEL etaBBN
void MODEL_NAMESPACE::etaBBN_to_etaBBN_rBBN_rCMB_dNurBBN_dNurCMB (const ModelParameters &myP, ModelParameters &targetP)
{
  logger()<<"Running interpret_as_parent calculations for etaBBN --> etaBBN_rBBN_rCMB_dNurBBN_dNurCMB ..."<<LogTags::info<<EOM;

  // Set eta_BBN
  targetP.setValue("eta_BBN", myP.getValue("eta_BBN"));

  // No dNeff due to an altered neutrino temperature
  targetP.setValue("r_BBN", 1.);
  targetP.setValue("r_CMB", 1.);

  // No dNeff due to additional radiation
  targetP.setValue("dNur_BBN", 0.);
  targetP.setValue("dNur_CMB", 0.);
}
#undef MODEL

#define MODEL rBBN_rCMB_dNurBBN_dNurCMB
#define PARENT etaBBN_rBBN_rCMB_dNurBBN_dNurCMB
void MODEL_NAMESPACE::rBBN_rCMB_dNurBBN_dNurCMB_to_etaBBN_rBBN_rCMB_dNurBBN_dNurCMB (const ModelParameters &myP, ModelParameters &targetP)
{
  USE_MODEL_PIPE(PARENT) // get pipe for "INTERPRET_AS_PARENT_DEPENDENCY"
  logger()<<"Running interpret_as_parent calculations for rBBN_rCMB_dNurBBN_dNurCMB --> etaBBN_rBBN_rCMB_dNurBBN_dNurCMB ..."<<LogTags::info<<EOM;

  // Set eta_BBN (This requires that eta0 / omega_b is known, i.e. LCDM is in use)
  // -- The dependency_resolver will figure that out --
  targetP.setValue("eta_BBN", *Dep::eta0);

  // Set the respective values of r
  targetP.setValue("r_BBN", myP.getValue("r_BBN"));
  targetP.setValue("r_CMB", myP.getValue("r_CMB"));

  // Set the respective values of dNeff
  targetP.setValue("dNur_BBN", myP.getValue("dNur_BBN"));
  targetP.setValue("dNur_CMB", myP.getValue("dNur_CMB"));
}
#undef PARENT
#undef MODEL

#define MODEL rBBN_rCMB
void MODEL_NAMESPACE::rBBN_rCMB_to_rBBN_rCMB_dNurBBN_dNurCMB (const ModelParameters &myP, ModelParameters &targetP)
{
  logger()<<"Running interpret_as_parent calculations for rBBN_rCMB --> rBBN_rCMB_dNurBBN_dNurCMB ..."<<LogTags::info<<EOM;

  // Set the respective values of r
  targetP.setValue("r_BBN", myP.getValue("r_BBN"));
  targetP.setValue("r_CMB", myP.getValue("r_CMB"));

  // No dNeff due to additional radiation
  targetP.setValue("dNur_BBN", 0.);
  targetP.setValue("dNur_CMB", 0.);
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

#define MODEL dNurBBN_dNurCMB
void MODEL_NAMESPACE::dNurBBN_dNurCMB_to_rBBN_rCMB_dNurBBN_dNurCMB (const ModelParameters &myP, ModelParameters &targetP)
{
  logger()<<"Running interpret_as_parent calculations for dNurBBN_dNurCMB --> rBBN_rCMB_dNurBBN_dNurCMB ..."<<LogTags::info<<EOM;

  // No dNeff due to an altered neutrino temperature
  targetP.setValue("r_BBN", 1.);
  targetP.setValue("r_CMB", 1.);

  // Set the respective values of dNeff
  targetP.setValue("dNur_BBN", myP.getValue("dNur_BBN"));
  targetP.setValue("dNur_CMB", myP.getValue("dNur_CMB"));
}
#undef MODEL

#define MODEL dNurCMB
void MODEL_NAMESPACE::dNurCMB_to_dNurBBN_dNurCMB (const ModelParameters &myP, ModelParameters &targetP)
{
  logger()<<"Running interpret_as_parent calculations for dNurCMB --> dNurBBN_dNurCMB ..."<<LogTags::info<<EOM;

  // Set the respective values of dNeff
  targetP.setValue("dNur_BBN", myP.getValue("dNur_CMB"));
  targetP.setValue("dNur_CMB", myP.getValue("dNur_CMB"));
}
#undef MODEL
