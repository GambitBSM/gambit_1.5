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

#define MODEL  LCDM_dNeffCMB_dNeffBBN
#define PARENT LCDM_dNeffCMB_dNeffBBN_etaBBN

// Translation function definition
void MODEL_NAMESPACE::LCDM_dNeffCMB_dNeffBBN_to_LCDM_dNeffCMB_dNeffBBN_etaBBN (const ModelParameters &myP, ModelParameters &targetP)
{
  USE_MODEL_PIPE(PARENT) // get pipe for "interpret as PARENT" function
  logger()<<"Running interpret_as_parent calculations for LCDM_dNeffCMB_dNeffBBN --> LCDM_dNeffCMB_dNeffBBN_etaBBN..."<<LogTags::info<<EOM;

  targetP.setValue("omega_b", myP.getValue("omega_b"));
  targetP.setValue("omega_cdm",myP.getValue("omega_cdm"));
  targetP.setValue("H0", myP.getValue("H0") );
  targetP.setValue("ln10A_s", myP.getValue("ln10A_s") );
  targetP.setValue("n_s", myP.getValue("n_s") );
  targetP.setValue("tau_reio", myP.getValue("tau_reio") );
  targetP.setValue("dNeff", myP.getValue("dNeff") );
  targetP.setValue("dNeff_BBN", myP.getValue("dNeff_BBN") );

  // eta_BBN = eta_CMB, computed from omega_b
  double ngamma, nb,eta_CMB;
  ngamma = 16*pi*zeta3*pow(*Dep::T_cmb*kb/hc,3); // photon number density today
  nb = myP.getValue("omega_b")*3*100*1e3*100*1e3/Mpc/Mpc/(8*pi*Gn*m_proton_g); // baryon number density today
  eta_CMB =  nb/ngamma;

  targetP.setValue("eta_BBN", eta_CMB );
}

#undef PARENT
#undef MODEL

#define MODEL  LCDM_dNeffCMB
#define PARENT LCDM_dNeffCMB_dNeffBBN

// Translation function definition
void MODEL_NAMESPACE::LCDM_dNeffCMB_to_LCDM_dNeffCMB_dNeffBBN (const ModelParameters &myP, ModelParameters &targetP)
{
  USE_MODEL_PIPE(PARENT) // get pipe for "interpret as PARENT" function
  logger()<<"Running interpret_as_parent calculations for LCDM_dNeffCMB --> LCDM_dNeffCMB_dNeffBBN..."<<LogTags::info<<EOM;

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

#define MODEL  LCDM_dNeffExt
#define PARENT LCDM_dNeffCMB_dNeffBBN_etaBBN

// Translation function definition
void MODEL_NAMESPACE::LCDM_dNeffExt_to_LCDM_dNeffCMB_dNeffBBN_etaBBN (const ModelParameters &myP, ModelParameters &targetP)
{
  USE_MODEL_PIPE(PARENT) // get pipe for "interpret as PARENT" function
  logger()<<"Running interpret_as_parent calculations for LCDM_dNeffExt --> LCDM_dNeffCMB_dNeffBBN_etaBBN..."<<LogTags::info<<EOM;

  targetP.setValue("omega_b", myP.getValue("omega_b"));
  targetP.setValue("omega_cdm",myP.getValue("omega_cdm"));
  targetP.setValue("H0", myP.getValue("H0") );
  targetP.setValue("ln10A_s", myP.getValue("ln10A_s") );
  targetP.setValue("n_s", myP.getValue("n_s") );
  targetP.setValue("tau_reio", myP.getValue("tau_reio") );

  // (PS) Needs fix. In the most genral case dNeffCMB and dNeffBBN are different.
  // The dependency should be either on std::vector<double> or on two distinct
  // capabilities.
  double ngamma, nb,eta_CMB;
  ngamma = 16*pi*zeta3*pow(*Dep::T_cmb*kb/hc,3); // photon number density today
  nb = myP.getValue("omega_b")*3*100*1e3*100*1e3/Mpc/Mpc/(8*pi*Gn*m_proton_g); // baryon number density today
  eta_CMB =  nb/ngamma;
  std::cout<< "Expected eta_CMB " << eta_CMB << std::endl;
  std::cout<< "Result Ext calc " << *Dep::external_dNeff_etaBBN << std::endl;
  //std::cout<< "Gotten eta_CMB "<< myP.getValue("eta_BBN")<<endl;
  double ratio = 1; 
  targetP.setValue("eta_BBN", eta_CMB*ratio );
  targetP.setValue("dNeff_BBN", 0. );
  targetP.setValue("dNeff", *Dep::external_dNeff_etaBBN );
}

#undef PARENT
#undef MODEL

#define MODEL  LCDM
#define PARENT LCDM_dNeffCMB

// Translation function definition
void MODEL_NAMESPACE::LCDM_to_LCDM_dNeffCMB (const ModelParameters &myP, ModelParameters &targetP)
{
  USE_MODEL_PIPE(PARENT) // get pipe for "interpret as PARENT" function
  logger()<<"Running interpret_as_parent calculations for LCDM --> LCDM_dNeffCMB..."<<LogTags::info<<EOM;

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

