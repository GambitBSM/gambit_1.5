//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Rollcall header for the developement
///  version of the CosmoBit module.
///
///  Compile-time registration of available
///  functions and variables from this backend.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Selim C. Hotinli
///          (selim.hotinli14@imperial.ac.uk)
///  \date 2017 Jul
///  \date 2017 Nov
///
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2017 Nov
///  \date 2018 Jan,Feb, Mar
///  \date 2019 Jan, Feb, June
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2018 Jun
///  \date 2019 Mar
///
///  \author Sanjay Bloor
///          (sanjay.bloor12@imperial.ac.uk)
///  \date 2019 June
///
///  *********************************************

#ifndef __CosmoBit_rollcall_hpp__
#define __CosmoBit_rollcall_hpp__

#include "gambit/CosmoBit/CosmoBit_types.hpp"

//#include <valarray>

#define MODULE CosmoBit
START_MODULE

  #define CAPABILITY injection_spectrum
  START_CAPABILITY
    #define FUNCTION injection_spectrum_ToyModel
    START_FUNCTION(DarkAges::injectionSpectrum)
    ALLOW_MODELS(TestDecayingDM)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY DM_fraction
  START_CAPABILITY
    #define FUNCTION DM_fraction_CosmoALP
    START_FUNCTION(double)
    DEPENDENCY(T_cmb, double)
    DEPENDENCY(minimum_abundance,double)
    DEPENDENCY(lifetime,double)
    ALLOW_MODEL_DEPENDENCE(CosmoALP,LCDM_dNeffCMB_dNeffBBN_etaBBN)
    MODEL_GROUP(alp,(CosmoALP))
    MODEL_GROUP(cosmo,(LCDM_dNeffCMB_dNeffBBN_etaBBN))
    ALLOW_MODEL_COMBINATION(cosmo,alp)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY lifetime
  START_CAPABILITY
    #define FUNCTION lifetime_CosmoALP
    START_FUNCTION(double)
    ALLOW_MODELS(CosmoALP)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY minimum_abundance
  START_CAPABILITY
    #define FUNCTION minimum_abundance_CosmoALP
    START_FUNCTION(double)
    ALLOW_MODELS(CosmoALP)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY energy_injection_efficiency
  START_CAPABILITY
    #define FUNCTION energy_injection_efficiency_func
    START_FUNCTION(DarkAges::fz_table)
    ALLOW_MODELS(TestDecayingDM)
    BACKEND_REQ(DA_efficiency_function, (DarkAges_tag), DarkAges::fz_table,())
    #undef FUNCTION
  #undef CAPABILITY
  
  #define CAPABILITY f_effective
  START_CAPABILITY
    #define FUNCTION f_effective_func
    START_FUNCTION(double)
    ALLOW_MODELS(TestDecayingDM)
    DEPENDENCY(energy_injection_efficiency,DarkAges::fz_table)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY lnL_A_planck
  START_CAPABILITY
    #define FUNCTION lnL_A_planck_gaussian
    START_FUNCTION(double)
    ALLOW_MODELS(Planck_TTTEEE,Planck_TT,Planck_lite)
    #undef FUNCTION
  #undef CAPABILITY
  
  #define CAPABILITY NuMasses_SM
  START_CAPABILITY
    #define FUNCTION set_NuMasses_SM
    START_FUNCTION(map_str_dbl)
    ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN,TestDecayingDM,StandardModel_SLHA2)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY class_Nur
  START_CAPABILITY
    #define FUNCTION set_class_Nur_LCDM_family
    START_FUNCTION(double)
    ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN,TestDecayingDM,StandardModel_SLHA2)
    DEPENDENCY(NuMasses_SM, map_str_dbl)
    #undef FUNCTION
    #define FUNCTION set_class_Nur_CosmoALP
    START_FUNCTION(double)
    MODEL_GROUP(SM,(StandardModel_SLHA2))
    MODEL_GROUP(cosmo,(LCDM))
    MODEL_GROUP(dark,(CosmoALP))
    ALLOW_MODEL_COMBINATION(cosmo,dark,SM)
    //ALLOW_MODELS(CosmoALP)
    //ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN,TestDecayingDM,StandardModel_SLHA2)
    DEPENDENCY(external_dNeff_etaBBN, map_str_dbl)
    DEPENDENCY(NuMasses_SM, map_str_dbl)
    #undef FUNCTION
  
  #undef CAPABILITY

  #define CAPABILITY class_set_parameter
  START_CAPABILITY

    #define FUNCTION class_set_parameter_LCDM_family
    START_FUNCTION(CosmoBit::Class_container)
    DEPENDENCY(Helium_abundance,std::vector<double>)
    DEPENDENCY(T_cmb, double)
    DEPENDENCY(T_ncdm, double)
    DEPENDENCY(class_Nur, double)
    DEPENDENCY(NuMasses_SM, map_str_dbl )
    ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN,TestDecayingDM,StandardModel_SLHA2)
    MODEL_CONDITIONAL_DEPENDENCY(lifetime,double,TestDecayingDM)
    MODEL_CONDITIONAL_DEPENDENCY(DM_fraction,double,TestDecayingDM)
    #undef FUNCTION


    #define FUNCTION class_set_parameter_LCDM_SingletDM
    START_FUNCTION(CosmoBit::Class_container)
    DEPENDENCY(T_cmb, double)
    DEPENDENCY(mwimp,double)
    DEPENDENCY(sigmav,double)
    ALLOW_MODEL_DEPENDENCE(LCDM,SingletDM)
    MODEL_GROUP(cosmo,(LCDM))
    MODEL_GROUP(dark,(SingletDM))
    ALLOW_MODEL_COMBINATION(cosmo,dark)
    #undef FUNCTION

    #define FUNCTION class_set_parameter_LCDMtensor
    START_FUNCTION(CosmoBit::Class_container)
    DEPENDENCY(T_cmb, double)
    ALLOW_MODELS(LCDMtensor)
    #undef FUNCTION

    #define FUNCTION class_set_parameter_inf_SR1quad_LCDMt
    START_FUNCTION(CosmoBit::Class_container)
    ALLOW_MODELS(inf_SR1quad_LCDMt)
    DEPENDENCY(T_cmb, double)
    BACKEND_REQ(multimodecode_gambit_driver,(modecode_tag), void, (gambit_inflation_observables*,int& ,int& ,  int& ,  int& ,int& ,  int& ,  int& ,  int& ,  int& ,int& ,double& ,int& ,  int& ,double& ,int& ,double*,double*,  int& ,  int& ,double*,double*,double*,double& ,double& ,double& ,  int& ,int& ,double& ,double*,double*,double*,double*,double& ,double&))
    #undef FUNCTION

    #define FUNCTION class_set_parameter_inf_1quarInf_LCDMt
    START_FUNCTION(CosmoBit::Class_container)
    ALLOW_MODELS(inf_1quarInf_LCDMt)
    DEPENDENCY(T_cmb, double)
    BACKEND_REQ(multimodecode_gambit_driver,(modecode_tag), void, (gambit_inflation_observables*,int& ,int& ,  int& ,  int& ,int& ,  int& ,  int& ,  int& ,  int& ,int& ,double& ,int& ,  int& ,double& ,int& ,double*,double*,  int& ,  int& ,double*,double*,double*,double& ,double& ,double& ,  int& ,int& ,double& ,double*,double*,double*,double*,double& ,double&))
    #undef FUNCTION

    #define FUNCTION class_set_parameter_inf_1mono32Inf_LCDMt
    START_FUNCTION(CosmoBit::Class_container)
    ALLOW_MODELS(inf_1mono32Inf_LCDMt)
    DEPENDENCY(T_cmb, double)
    BACKEND_REQ(multimodecode_gambit_driver,(modecode_tag), void, (gambit_inflation_observables*,int& ,int& ,  int& ,  int& ,int& ,  int& ,  int& ,  int& ,  int& ,int& ,double& ,int& ,  int& ,double& ,int& ,double*,double*,  int& ,  int& ,double*,double*,double*,double& ,double& ,double& ,  int& ,int& ,double& ,double*,double*,double*,double*,double& ,double&))
    #undef FUNCTION

    #define FUNCTION class_set_parameter_inf_1linearInf_LCDMt
    START_FUNCTION(CosmoBit::Class_container)
    ALLOW_MODELS(inf_1linearInf_LCDMt)
    DEPENDENCY(T_cmb, double)
    BACKEND_REQ(multimodecode_gambit_driver,(modecode_tag), void, (gambit_inflation_observables*,int& ,int& ,  int& ,  int& ,int& ,  int& ,  int& ,  int& ,  int& ,int& ,double& ,int& ,  int& ,double& ,int& ,double*,double*,  int& ,  int& ,double*,double*,double*,double& ,double& ,double& ,  int& ,int& ,double& ,double*,double*,double*,double*,double& ,double&))
    #undef FUNCTION

    #define FUNCTION class_set_parameter_inf_1naturalInf_LCDMt
    START_FUNCTION(CosmoBit::Class_container)
    ALLOW_MODELS(inf_1naturalInf_LCDMt)
    DEPENDENCY(T_cmb, double)
    BACKEND_REQ(multimodecode_gambit_driver,(modecode_tag), void, (gambit_inflation_observables*,int& ,int& ,  int& ,  int& ,int& ,  int& ,  int& ,  int& ,  int& ,int& ,double& ,int& ,  int& ,double& ,int& ,double*,double*,  int& ,  int& ,double*,double*,double*,double& ,double& ,double& ,  int& ,int& ,double& ,double*,double*,double*,double*,double& ,double&))
    #undef FUNCTION


    #define FUNCTION class_set_parameter_inf_1hilltopInf_LCDMt
    START_FUNCTION(CosmoBit::Class_container)
    ALLOW_MODELS(inf_1hilltopInf_LCDMt)
    DEPENDENCY(T_cmb, double)
    BACKEND_REQ(multimodecode_gambit_driver,(modecode_tag), void, (gambit_inflation_observables*,int& ,int& ,  int& ,  int& ,int& ,  int& ,  int& ,  int& ,  int& ,int& ,double& ,int& ,  int& ,double& ,int& ,double*,double*,  int& ,  int& ,double*,double*,double*,double& ,double& ,double& ,  int& ,int& ,double& ,double*,double*,double*,double*,double& ,double&))
    #undef FUNCTION


    #define FUNCTION class_set_parameter_inf_smashInf_LCDMt
    START_FUNCTION(CosmoBit::Class_container)
    ALLOW_MODELS(inf_smashInf_LCDMt)
    DEPENDENCY(T_cmb, double)
    BACKEND_REQ(multimodecode_gambit_driver,(modecode_tag), void, (gambit_inflation_observables*,int& ,int& ,  int& ,  int& ,int& ,  int& ,  int& ,  int& ,  int& ,int& ,double& ,int& ,  int& ,double& ,int& ,double*,double*,  int& ,  int& ,double*,double*,double*,double& ,double& ,double& ,  int& ,int& ,double& ,double*,double*,double*,double*,double& ,double&))
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY Cl_TT
  START_CAPABILITY
    #define FUNCTION class_get_Cl_TT
    START_FUNCTION(std::vector<double>)
    BACKEND_REQ(class_get_cl,(class_tag),std::vector<double>, (str))
    FORCE_SAME_BACKEND(class_tag)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY Cl_TE
  START_CAPABILITY
    #define FUNCTION class_get_Cl_TE
    START_FUNCTION(std::vector<double>)
    BACKEND_REQ(class_get_cl,(class_tag),std::vector<double>, (str))
    FORCE_SAME_BACKEND(class_tag)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY Cl_EE
  START_CAPABILITY
    #define FUNCTION class_get_Cl_EE
    START_FUNCTION(std::vector<double>)
    BACKEND_REQ(class_get_cl,(class_tag),std::vector<double>, (str))
    FORCE_SAME_BACKEND(class_tag)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY Cl_BB
  START_CAPABILITY
    #define FUNCTION class_get_Cl_BB
    START_FUNCTION(std::vector<double>)
    BACKEND_REQ(class_get_cl,(class_tag),std::vector<double>, (str))
    FORCE_SAME_BACKEND(class_tag)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY Cl_PhiPhi
  START_CAPABILITY
    #define FUNCTION class_get_Cl_PhiPhi
    START_FUNCTION(std::vector<double>)
    BACKEND_REQ(class_get_cl,(class_tag),std::vector<double>, (str))
    FORCE_SAME_BACKEND(class_tag)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY compute_Planck_lowp_TT_loglike
  START_CAPABILITY
    #define FUNCTION function_Planck_lowp_TT_loglike
    START_FUNCTION(double)
    DEPENDENCY(Cl_TT,std::vector<double>)
    DEPENDENCY(Cl_TE,std::vector<double>)
    DEPENDENCY(Cl_EE,std::vector<double>)
    DEPENDENCY(Cl_BB,std::vector<double>)
    ALLOW_MODELS(Planck_TTTEEE,Planck_TT,Planck_lite)
    BACKEND_REQ(clik_compute_loglike, (clik_tag), double ,(clik_object*,double*,clik_error**))
    BACKEND_REQ(clik_initialize_error, (clik_tag), clik_error* ,())
    BACKEND_REQ(return_lowp_TT,(clik_tag),clik_object*,())
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY compute_Planck_high_TT_loglike
  START_CAPABILITY
    #define FUNCTION function_Planck_high_TT_loglike
    START_FUNCTION(double)
    DEPENDENCY(Cl_TT,std::vector<double>)
    ALLOW_MODELS(Planck_TT)
    BACKEND_REQ(clik_compute_loglike, (clik_tag), double ,(clik_object*,double*,clik_error**))
    BACKEND_REQ(clik_initialize_error, (clik_tag), clik_error* ,())
    BACKEND_REQ(return_high_TT,(clik_tag),clik_object*,())
    #undef FUNCTION

    #define FUNCTION function_Planck_high_TT_lite_loglike
    START_FUNCTION(double)
    DEPENDENCY(Cl_TT,std::vector<double>)
    ALLOW_MODELS(Planck_lite)
    BACKEND_REQ(clik_compute_loglike, (clik_tag), double ,(clik_object*,double*,clik_error**))
    BACKEND_REQ(clik_initialize_error, (clik_tag), clik_error* ,())
    BACKEND_REQ(return_high_TT_lite,(clik_tag),clik_object*,())
    #undef FUNCTION

    #define FUNCTION function_Planck_high_TTTEEE_loglike
    START_FUNCTION(double)
    DEPENDENCY(Cl_TT,std::vector<double>)
    DEPENDENCY(Cl_TE,std::vector<double>)
    DEPENDENCY(Cl_EE,std::vector<double>)
    ALLOW_MODELS(Planck_TTTEEE)
    BACKEND_REQ(clik_compute_loglike, (clik_tag), double ,(clik_object*,double*,clik_error**))
    BACKEND_REQ(clik_initialize_error, (clik_tag), clik_error* ,())
    BACKEND_REQ(return_high_TTTEEE,(clik_tag),clik_object*,())
    #undef FUNCTION

    #define FUNCTION function_Planck_high_TTTEEE_lite_loglike
    START_FUNCTION(double)
    DEPENDENCY(Cl_TT,std::vector<double>)
    DEPENDENCY(Cl_TE,std::vector<double>)
    DEPENDENCY(Cl_EE,std::vector<double>)
    ALLOW_MODELS(Planck_lite)
    BACKEND_REQ(clik_compute_loglike, (clik_tag), double ,(clik_object*,double*,clik_error**))
    BACKEND_REQ(clik_initialize_error, (clik_tag), clik_error* ,())
    BACKEND_REQ(return_high_TTTEEE_lite,(clik_tag),clik_object*,())
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY compute_Planck_lensing_loglike
  START_CAPABILITY
    #define FUNCTION function_Planck_lensing_loglike
    START_FUNCTION(double)
    DEPENDENCY(Cl_TT,std::vector<double>)
    DEPENDENCY(Cl_TE,std::vector<double>)
    DEPENDENCY(Cl_EE,std::vector<double>)
    DEPENDENCY(Cl_PhiPhi,std::vector<double>)
    ALLOW_MODELS(Planck_TTTEEE,Planck_TT,Planck_lite)
    BACKEND_REQ(clik_lensing_compute_loglike, (clik_tag), double ,(clik_lensing_object*,double*,clik_error**))
    BACKEND_REQ(clik_initialize_error, (clik_tag), clik_error* ,())
    BACKEND_REQ(return_lensing,(clik_tag),clik_lensing_object*,())
    #undef FUNCTION
  #undef CAPABILITY


  #define CAPABILITY T_cmb
    START_CAPABILITY
    #define FUNCTION set_T_cmb
      START_FUNCTION(double)
    #undef FUNCTION
  #undef CAPABILITY
  
#define CAPABILITY T_ncdm_SM // needed in addition to T_ncdm since T_ncdm of non-SM models assume a fiducial value to base calculation on 
    START_CAPABILITY
    #define FUNCTION set_T_ncdm_SM
      START_FUNCTION(double)
    #undef FUNCTION
  #undef CAPABILITY
  
  #define CAPABILITY T_ncdm
    START_CAPABILITY
    #define FUNCTION set_T_ncdm
      START_FUNCTION(double)
      DEPENDENCY(T_ncdm_SM,double)
    #undef FUNCTION
    #define FUNCTION set_T_ncdm_CosmoALP
      START_FUNCTION(double)
      ALLOW_MODELS(CosmoALP)
      DEPENDENCY(external_dNeff_etaBBN, map_str_dbl)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY n0_g
    START_CAPABILITY

    #define FUNCTION compute_n0_g
      START_FUNCTION(double)
      DEPENDENCY(T_cmb, double)
      ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN)
    #undef FUNCTION

    // #define FUNCTION get_n0_g_classy
    //   START_FUNCTION(double)
    //   DEPENDENCY(get_Classy_cosmo_container, CosmoBit::Classy_cosmo_container)
    //   ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN)
    // #undef FUNCTION

  #undef CAPABILITY
  
  #define CAPABILITY Omega0_m
    START_CAPABILITY

    #define FUNCTION compute_Omega0_m
      START_FUNCTION(double)
      DEPENDENCY(Omega0_b, double)
      DEPENDENCY(Omega0_cdm, double)
      DEPENDENCY(Omega0_ncdm, double)
      ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN)
    #undef FUNCTION

    #define FUNCTION get_Omega0_m_classy
      START_FUNCTION(double)
      BACKEND_REQ(class_get_Omega0_m,(classy),double,())
      ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN)
    #undef FUNCTION

  #undef CAPABILITY


  #define CAPABILITY Omega0_b
    START_CAPABILITY
    
    #define FUNCTION compute_Omega0_b
      START_FUNCTION(double)
      ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN)
    #undef FUNCTION

  #undef CAPABILITY

  
  #define CAPABILITY Omega0_cdm
    START_CAPABILITY

    #define FUNCTION compute_Omega0_cdm
      START_FUNCTION(double)
      ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN)
    #undef FUNCTION

  #undef CAPABILITY

  
  #define CAPABILITY Omega0_r
    START_CAPABILITY

    #define FUNCTION compute_Omega0_r
      START_FUNCTION(double)
      DEPENDENCY(Omega0_g, double)
      DEPENDENCY(Omega0_ur, double)
      ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN)
    #undef FUNCTION

    #define FUNCTION get_Omega0_r_classy
      START_FUNCTION(double)
      BACKEND_REQ(class_get_Omega0_r,(classy),double,())
    #undef FUNCTION

  #undef CAPABILITY

  
  #define CAPABILITY Omega0_g
    START_CAPABILITY

    #define FUNCTION compute_Omega0_g
      START_FUNCTION(double)
      DEPENDENCY(T_cmb, double)
      ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN)
    #undef FUNCTION

  #undef CAPABILITY
  

  #define CAPABILITY Omega0_ur
    START_CAPABILITY

    #define FUNCTION compute_Omega0_ur
      START_FUNCTION(double)
        DEPENDENCY(Omega0_g, double)
        DEPENDENCY(class_Nur, double)
    #undef FUNCTION

    #define FUNCTION get_Omega0_ur_classy  
      START_FUNCTION(double)
        BACKEND_REQ(class_get_Omega0_ur,(classy),double,())
    #undef FUNCTION

  #undef CAPABILITY


  #define CAPABILITY Omega0_ncdm_tot
    START_CAPABILITY

    #define FUNCTION get_Omega0_ncdm_tot_classy 
      START_FUNCTION(double)
        BACKEND_REQ(class_get_Omega0_ncdm_tot,(classy),double,())
    #undef FUNCTION

  #undef CAPABILITY
  

  #define CAPABILITY Omega0_ncdm
    START_CAPABILITY

    #define FUNCTION compute_Omega0_ncdm
      START_FUNCTION(double)
      DEPENDENCY(T_cmb, double)
      ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN,StandardModel_SLHA2)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY eta0
    START_CAPABILITY

    #define FUNCTION calculate_eta0
      START_FUNCTION(double)
      DEPENDENCY(T_cmb, double)
      ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN) // To get etaCMB for LCDM_dNeffCMB_dNeffBBN_etaBBN
      ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN) // Allow for all direct childs of LCDM_dNeffCMB_dNeffBBN_etaBBN. Needed for the translation into LCDM_dNeffCMB_dNeffBBN_etaBBN
      ALLOW_MODELS(LCDM_ExtdNeffCMB_ExtetaBBN)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY rs_drag
  START_CAPABILITY

    #define FUNCTION get_rs_drag_classy
      START_FUNCTION(double)
      BACKEND_REQ(class_get_rs,(classy),double,())
    #undef FUNCTION

  #undef CAPABILITY

  /// Good for cross-checks, innit.
  #define CAPABILITY Neff
  START_CAPABILITY

    #define FUNCTION get_Neff_classy
      START_FUNCTION(double)
      BACKEND_REQ(class_get_Neff,(classy),double,())
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY etaCMB
    START_CAPABILITY

    #define FUNCTION calculate_etaCMB_SM // here: eta0 = etaCMB
      START_FUNCTION(double)
      DEPENDENCY(eta0, double)
      ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN) // To get etaCMB for LCDM_dNeffCMB_dNeffBBN_etaBBN
      ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN) // Allow for all direct childs of LCDM_dNeffCMB_dNeffBBN_etaBBN. Needed for the translation into LCDM_dNeffCMB_dNeffBBN_etaBBN
      ALLOW_MODELS(LCDM_ExtdNeffCMB_ExtetaBBN)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY etaBBN
   START_CAPABILITY

   #define FUNCTION set_etaBBN // etaBBN is model parameter_dict_for_MPLike
      START_FUNCTION(double)
      ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN) // To get etaCMB for LCDM_dNeffCMB_dNeffBBN_etaBBN
    #undef FUNCTION
   #define FUNCTION calculate_etaBBN_SM // etaBBN = etaCMB
      START_FUNCTION(double)
      DEPENDENCY(etaCMB, double)
      ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN) // Allow for all direct childs of LCDM_dNeffCMB_dNeffBBN_etaBBN. Needed for the translation into LCDM_dNeffCMB_dNeffBBN_etaBBN
    #undef FUNCTION
   #define FUNCTION calculate_etaBBN_ALP
     START_FUNCTION(double)
        DEPENDENCY(etaCMB, double)
        DEPENDENCY(external_dNeff_etaBBN, map_str_dbl)
        ALLOW_MODELS(LCDM_ExtdNeffCMB_ExtetaBBN) // To make sure this function is used to fulfill etaBBN capability in model translation function 
   #undef FUNCTION
  #undef CAPABILITY

  // compute dNeff AND etaBBN for non-standard models
  #define CAPABILITY external_dNeff_etaBBN
    START_CAPABILITY
    #define FUNCTION compute_dNeff_etaBBN_ALP
     START_FUNCTION(map_str_dbl)
     ALLOW_MODELS(CosmoALP)
     DEPENDENCY(lifetime, double)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY ExtdNeffCMB
   START_CAPABILITY
   #define FUNCTION calculate_dNeffCMB_ALP
     START_FUNCTION(double)
     DEPENDENCY(external_dNeff_etaBBN, map_str_dbl)
   #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY Sigma8
    START_CAPABILITY

    #define FUNCTION compute_Sigma8
      START_FUNCTION(double)
      DEPENDENCY(Omega0_m, double)
      BACKEND_REQ(class_get_sigma8,(class_tag),double,(double))
    #undef FUNCTION
  
    #define FUNCTION get_Sigma8_classy
      START_FUNCTION(double)
      DEPENDENCY(Omega0_m, double)
      BACKEND_REQ(class_get_sigma8,(class_tag),double,())
    #undef FUNCTION

  #undef CAPABILITY

// AlterBBN related functions & capabilities
  #define CAPABILITY AlterBBN_setInput
    START_CAPABILITY
    #define FUNCTION AlterBBN_Input_LCDM_dNeffCMB_dNeffBBN_etaBBN
     START_FUNCTION(map_str_dbl)
     ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN)
    #undef FUNCTION
  #undef CAPABILITY


  #define CAPABILITY Helium_abundance
   START_CAPABILITY
    #define FUNCTION get_Helium_abundance
      START_FUNCTION(std::vector<double>)
      DEPENDENCY(BBN_abundances, CosmoBit::BBN_container)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY Deuterium_abundance
   START_CAPABILITY
    #define FUNCTION get_Deuterium_abundance
      START_FUNCTION(std::vector<double>)
      DEPENDENCY(BBN_abundances, CosmoBit::BBN_container)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY Helium3_abundance
   START_CAPABILITY
    #define FUNCTION get_Helium3_abundance
      START_FUNCTION(std::vector<double>)
      DEPENDENCY(BBN_abundances, CosmoBit::BBN_container)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY Lithium7_abundance
   START_CAPABILITY
    #define FUNCTION get_Lithium7_abundance
      START_FUNCTION(std::vector<double>)
      DEPENDENCY(BBN_abundances, CosmoBit::BBN_container)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY Beryllium7_abundance
   START_CAPABILITY
    #define FUNCTION get_Beryllium7_abundance
      START_FUNCTION(std::vector<double>)
      DEPENDENCY(BBN_abundances, CosmoBit::BBN_container)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY BBN_abundances
  START_CAPABILITY
    #define FUNCTION compute_BBN_abundances
      START_FUNCTION(CosmoBit::BBN_container)
      DEPENDENCY(AlterBBN_setInput, map_str_dbl)
      BACKEND_REQ(call_nucl_err, (libbbn), int, (map_str_dbl&,double*,double*))
      BACKEND_REQ(get_NNUC, (libbbn), int, ())
      BACKEND_OPTION( (AlterBBN), (libbbn) )
    #undef FUNCTION
  #undef CAPABILITY

 #define CAPABILITY BBN_LogLike
   START_CAPABILITY
   #define FUNCTION compute_BBN_LogLike
   START_FUNCTION(double)
   DEPENDENCY(BBN_abundances, CosmoBit::BBN_container)
   DEPENDENCY(AlterBBN_setInput, map_str_dbl)
  #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY H0_LogLike
   START_CAPABILITY
    #define FUNCTION compute_H0_LogLike
    START_FUNCTION(double)
    ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY Pantheon_LogLike
   START_CAPABILITY
   #define FUNCTION compute_Pantheon_LogLike
    START_FUNCTION(double)
    ALLOW_MODELS(cosmo_nuisance_params_Pantheon)
    BACKEND_REQ(class_get_Dl,(class_tag),double,(double))
   #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY BAO_LogLike
   START_CAPABILITY
   #define FUNCTION compute_BAO_LogLike
    START_FUNCTION(double)
    BACKEND_REQ(class_get_Da,(class_tag),double,(double))
    BACKEND_REQ(class_get_Hz,(class_tag),double,(double))
    BACKEND_REQ(class_get_rs,(class_tag),double,())
    FORCE_SAME_BACKEND(class_tag)
   #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY sigma8_LogLike
     START_CAPABILITY
     #define FUNCTION compute_sigma8_LogLike
      START_FUNCTION(double)
      DEPENDENCY(Omega0_m, double)
      ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN)
      BACKEND_REQ(class_get_sigma8,(class_tag),double,(double))
     #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY parameter_dict_for_MPLike
     START_CAPABILITY
     #define FUNCTION set_parameter_dict_for_MPLike
      START_FUNCTION(pybind11::dict)
      ALLOW_MODELS(cosmo_nuisance_params_JLA,cosmo_nuisance_params_BK14,cosmo_nuisance_params_CFHTLens_correlation)
      ALLOW_MODELS(cosmo_nuisance_params_euclid_lensing,cosmo_nuisance_params_euclid_pk,cosmo_nuisance_params_ISW)
      ALLOW_MODELS(cosmo_nuisance_params_kids450_qe_likelihood_public,cosmo_nuisance_params_ska,cosmo_nuisance_params_wmap)
      // if you implement new MontePython likelihoods with new nuisance parameters add the name of your new
      // nuisance parameter model (to be defined in Models/include/gambit/Models/models/CosmoNuisanceModels.hpp)
      //ALLOW_MODELS(cosmo_nuisance_params_FOR_YOUR_NEW_LIKE) 
     #undef FUNCTION
     #define FUNCTION pass_empty_parameter_dict_for_MPLike
       START_FUNCTION(pybind11::dict)
     #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY MP_experiment_names
     START_CAPABILITY
     #define FUNCTION set_MP_experiment_names
      START_FUNCTION(map_str_str)
      BACKEND_REQ(get_MP_availible_likelihoods, (libmontepythonlike), std::vector<str>, ())
     #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY cosmo_args_from_MPLike
     START_CAPABILITY
     #define FUNCTION init_cosmo_args_from_MPLike
      START_FUNCTION(pybind11::dict)
      DEPENDENCY(MP_experiment_names, map_str_str)
      BACKEND_REQ(create_MP_likelihood_objects, (libmontepythonlike), map_str_pyobj,    (pybind11::object&, map_str_str&))
      BACKEND_REQ(create_MP_data_object,        (libmontepythonlike), pybind11::object, (map_str_str&))
     #undef FUNCTION
  #undef CAPABILITY

    #define CAPABILITY get_classy_cosmo_container
     START_CAPABILITY
     #define FUNCTION init_classy_cosmo_container_with_MPLike
      START_FUNCTION(CosmoBit::ClassyInput)
      ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN)
      DEPENDENCY(cosmo_args_from_MPLike,  pybind11::dict)
      DEPENDENCY(set_classy_parameters,   CosmoBit::ClassyInput)
     #undef FUNCTION
     #define FUNCTION init_classy_cosmo_container
      START_FUNCTION(CosmoBit::ClassyInput)
      ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN)
      DEPENDENCY(set_classy_parameters, CosmoBit::ClassyInput)
     #undef FUNCTION
  #undef CAPABILITY

    #define CAPABILITY set_classy_parameters
     START_CAPABILITY
     #define FUNCTION set_classy_parameters_LCDM
      START_FUNCTION(CosmoBit::ClassyInput)
      ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN)
      DEPENDENCY(NuMasses_SM, map_str_dbl)
      DEPENDENCY(T_cmb,             double)
      DEPENDENCY(T_ncdm,            double)
      DEPENDENCY(class_Nur,         double)
      DEPENDENCY(Helium_abundance,  std::vector<double>)
      MODEL_CONDITIONAL_DEPENDENCY(model_dep_classy_parameters, pybind11::dict,  TestDecayingDM)
      //MODEL_CONDITIONAL_DEPENDENCY(lifetime,                    double,             TestDecayingDM)
      //MODEL_CONDITIONAL_DEPENDENCY(DM_fraction,                 double,             TestDecayingDM)
      //MODEL_CONDITIONAL_DEPENDENCY(energy_injection_efficiency, DarkAges::fz_table, TestDecayingDM)
     #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY model_dep_classy_parameters
     START_CAPABILITY
     #define FUNCTION model_dep_classy_parameters_TestDecayingDM
      START_FUNCTION(pybind11::dict)
      ALLOW_MODELS(TestDecayingDM)
      DEPENDENCY(lifetime,                    double)
      DEPENDENCY(DM_fraction,                 double)
      DEPENDENCY(energy_injection_efficiency, DarkAges::fz_table)
     #undef FUNCTION
  #undef CAPABILITY
  
  /// Calculates lnL for each experiment using the experiment names
  #define CAPABILITY MP_LogLikes
    START_CAPABILITY
    #define FUNCTION calc_MP_LogLikes
      START_FUNCTION(map_str_dbl) 
      ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN)
      DEPENDENCY(MP_experiment_names,         map_str_str)
      DEPENDENCY(parameter_dict_for_MPLike,   pybind11::dict)
      BACKEND_REQ(get_classy_cosmo_object,               (classy),             pybind11::object,      ())
      BACKEND_REQ(get_MP_loglike,               (libmontepythonlike), double,           (const CosmoBit::MPLike_data_container&, pybind11::object&, std::string&))
      BACKEND_REQ(create_MP_data_object,        (libmontepythonlike), pybind11::object, (map_str_str&))
      BACKEND_REQ(create_MP_likelihood_objects, (libmontepythonlike), map_str_pyobj,    (pybind11::object&, map_str_str&))
    #undef FUNCTION
  #undef CAPABILITY
      
  /// Calculates the total lnL from MontePython 
  #define CAPABILITY MP_Combined_LogLike
    START_CAPABILITY
    #define FUNCTION calc_MP_combined_LogLike
      START_FUNCTION(double)
      DEPENDENCY(MP_LogLikes, map_str_dbl)
      ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN)
    #undef FUNCTION
  #undef CAPABILITY

#undef MODULE

#endif /* defined __CosmoBit_rollcall_hpp__ */
