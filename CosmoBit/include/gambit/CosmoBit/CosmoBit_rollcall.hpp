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
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2018 Jun
///  \date 2019 Mar
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
    BACKEND_REQ(DA_efficiency_function, (DarkAges_tag) ,DarkAges::fz_table,())
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

  #define CAPABILITY class_get_spectra
  START_CAPABILITY
    #define FUNCTION class_get_spectra_func
    START_FUNCTION(CosmoBit::Class_container)
    BACKEND_REQ(get_ptr_to_class,(class_tag),CosmoBit::Class_container,())
    BACKEND_REQ(class_get_cl,(class_tag),std::vector< double >, (str))
    FORCE_SAME_BACKEND(class_tag)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY compute_Planck_lowp_TT_loglike
  START_CAPABILITY
    #define FUNCTION function_Planck_lowp_TT_loglike
    START_FUNCTION(double)
    DEPENDENCY(class_get_spectra,CosmoBit::Class_container)
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
    DEPENDENCY(class_get_spectra,CosmoBit::Class_container)
    ALLOW_MODELS(Planck_TT)
    BACKEND_REQ(clik_compute_loglike, (clik_tag), double ,(clik_object*,double*,clik_error**))
    BACKEND_REQ(clik_initialize_error, (clik_tag), clik_error* ,())
    BACKEND_REQ(return_high_TT,(clik_tag),clik_object*,())
    #undef FUNCTION

    #define FUNCTION function_Planck_high_TTTEEE_loglike
    START_FUNCTION(double)
    DEPENDENCY(class_get_spectra,CosmoBit::Class_container)
    ALLOW_MODELS(Planck_TTTEEE)
    BACKEND_REQ(clik_compute_loglike, (clik_tag), double ,(clik_object*,double*,clik_error**))
    BACKEND_REQ(clik_initialize_error, (clik_tag), clik_error* ,())
    BACKEND_REQ(return_high_TTTEEE,(clik_tag),clik_object*,())
    #undef FUNCTION

    #define FUNCTION function_Planck_high_TT_lite_loglike
    START_FUNCTION(double)
    DEPENDENCY(class_get_spectra,CosmoBit::Class_container)
    ALLOW_MODELS(Planck_lite)
    BACKEND_REQ(clik_compute_loglike, (clik_tag), double ,(clik_object*,double*,clik_error**))
    BACKEND_REQ(clik_initialize_error, (clik_tag), clik_error* ,())
    BACKEND_REQ(return_high_TT_lite,(clik_tag),clik_object*,())
    #undef FUNCTION

    #define FUNCTION function_Planck_high_TTTEEE_lite_loglike
    START_FUNCTION(double)
    DEPENDENCY(class_get_spectra,CosmoBit::Class_container)
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
    DEPENDENCY(class_get_spectra,CosmoBit::Class_container)
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
  
#define CAPABILITY T_ncdm_SM // needed in additino to T_ncdm since T_ncdm of non-SM models assume a fiducial value to base calculation on 
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
        //ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN)
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
   #define FUNCTION set_etaBBN // etaBBN is model parameter
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
  #undef CAPABILITY

// AlterBBN related functions & capabilities
  #define CAPABILITY AlterBBN_modelinfo
    START_CAPABILITY
    #define FUNCTION AlterBBN_fill_LCDM_dNeffCMB_dNeffBBN_etaBBN
     START_FUNCTION(relicparam)
     ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN)
     BACKEND_OPTION( (AlterBBN, 2.0), (libbbn) )
     BACKEND_REQ(Init_cosmomodel, (libbbn), void, (relicparam*))
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
      DEPENDENCY(AlterBBN_modelinfo, relicparam)
      BACKEND_REQ(nucl_err, (libbbn), int, (const relicparam*,double*,double*))
      BACKEND_REQ(get_NNUC, (libbbn), int, ())
    #undef FUNCTION
  #undef CAPABILITY

 #define CAPABILITY BBN_LogLike
   START_CAPABILITY
   #define FUNCTION compute_BBN_LogLike
   START_FUNCTION(double)
   DEPENDENCY(BBN_abundances, CosmoBit::BBN_container)
   BACKEND_OPTION( (AlterBBN, 2.0), (libbbn) )
   DEPENDENCY(AlterBBN_modelinfo, relicparam)
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
    ALLOW_MODEL_DEPENDENCE(LCDM_dNeffCMB_dNeffBBN_etaBBN,cosmo_nuisance_params)
    MODEL_GROUP(cosmology, (LCDM_dNeffCMB_dNeffBBN_etaBBN))
    MODEL_GROUP(nuisance, (cosmo_nuisance_params))
    ALLOW_MODEL_COMBINATION(cosmology,nuisance)
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

#undef MODULE

#endif /* defined __CosmoBit_rollcall_hpp__ */
