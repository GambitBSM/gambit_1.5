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
///  *********************************************

#ifndef __CosmoBit_rollcall_hpp__
#define __CosmoBit_rollcall_hpp__

#include "gambit/CosmoBit/CosmoBit_types.hpp"

#define MODULE CosmoBit
START_MODULE

  #define CAPABILITY injection_spectrum
  START_CAPABILITY
    #define FUNCTION injection_spectrum_ToyModel
    START_FUNCTION(DarkAges::injectionSpectrum)
    ALLOW_MODELS(TestDecayingDM)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY DM_mass
  START_CAPABILITY
    #define FUNCTION DM_mass_ToyModel
    START_FUNCTION(double)
    ALLOW_MODELS(TestDecayingDM)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY DM_fraction
  START_CAPABILITY
    #define FUNCTION DM_fraction_ToyModel
    START_FUNCTION(double)
    ALLOW_MODELS(TestDecayingDM)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY lifetime
  START_CAPABILITY
    #define FUNCTION lifetime_ToyModel
    START_FUNCTION(double)
    ALLOW_MODELS(TestDecayingDM)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY f_effective
  START_CAPABILITY
    #define FUNCTION f_effective_func
    START_FUNCTION(double)
    ALLOW_MODELS(TestDecayingDM)
    BACKEND_REQ(DA_efficiency_function, (DarkAges_tag),DarkAges::fz_table,())
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY lnL_A_planck
  START_CAPABILITY
    #define FUNCTION lnL_A_planck_gaussian
    START_FUNCTION(double)
    ALLOW_MODELS(Planck_TTTEEE,Planck_TT,Planck_lite)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY class_set_Smu
  START_CAPABILITY
    #define FUNCTION class_set_Smu_LCDM
    START_FUNCTION(CosmoBit::Class_container)
    ALLOW_MODELS(LCDM)
    #undef FUNCTION

    #define FUNCTION class_set_Smu_LCDM_Smu_dNeffCMB_dNeffBBN_etaBBN
    START_FUNCTION(CosmoBit::Class_container)
    ALLOW_MODELS(LCDM_Smu_dNeffCMB_dNeffBBN_etaBBN)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY class_set_parameter
  START_CAPABILITY

    #define FUNCTION class_set_parameter_LCDM_family
    START_FUNCTION(CosmoBit::Class_container)
    DEPENDENCY(Helium_abundance,std::vector<double>)
    DEPENDENCY(T_cmb, double)
    DEPENDENCY(class_set_Smu, CosmoBit::Class_container)
    ALLOW_MODELS(LCDM,LCDM_Smu_dNeffCMB_dNeffBBN_etaBBN)
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



// ------------------------


  #define CAPABILITY T_cmb
    START_CAPABILITY
    #define FUNCTION set_T_cmb
      START_FUNCTION(double)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY eta
    START_CAPABILITY
    #define FUNCTION calculate_eta
      START_FUNCTION(double)
      DEPENDENCY(T_cmb, double)
      ALLOW_MODELS(LCDM)
      // TODO: atm calculation of eta implemented twice: once in CosmoModels, once here. put default LCDM into model tree, hand nu masses!
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY dNeffExt
    START_CAPABILITY
    #define FUNCTION compute_dNeffExt_ALP
      START_FUNCTION(double)
      //ALLOW_MODELS(ALP) # TODO: refer to correct model name
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY etaBBN_ALP
    START_CAPABILITY
    #define FUNCTION compute_etaBBN_ALP
      START_FUNCTION(double)
      // TODO: refer to correct model name
      //MODEL_GROUP(cosmology, (LCDM, LCDM_Smu_dNeffCMB_dNeffBBN_etaBBN))
      //MODEL_GROUP(particle, (ALP))
      //ALLOW_MODEL_COMBINATION(cosmology,particle)
    #undef FUNCTION
  #undef CAPABILITY

// AlterBBN related functions & capabilities
  #define CAPABILITY AlterBBN_modelinfo
    START_CAPABILITY
    #define FUNCTION AlterBBN_fill_LCDM
      START_FUNCTION(relicparam)
      ALLOW_MODELS(LCDM)
      DEPENDENCY(eta, double)
      BACKEND_OPTION( (AlterBBN, 2.0), (libbbn) )
      BACKEND_REQ(Init_cosmomodel, (libbbn), void, (relicparam*))
    #undef FUNCTION

    #define FUNCTION AlterBBN_fill_LCDM_Smu_dNeffCMB_dNeffBBN_etaBBN
     START_FUNCTION(relicparam)
     ALLOW_MODELS(LCDM_Smu_dNeffCMB_dNeffBBN_etaBBN)
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
   //BACKEND_REQ(bbn_excluded_chi2, (libbbn), int, (const relicparam*))
   ALLOW_MODELS(LCDM,LCDM_Smu_dNeffCMB_dNeffBBN_etaBBN)
  #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY H0_LogLike
   START_CAPABILITY
    #define FUNCTION compute_H0_LogLike
    START_FUNCTION(double)
    ALLOW_MODELS(LCDM, LCDM_Smu_dNeffCMB_dNeffBBN_etaBBN)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY Pantheon_LogLike
   START_CAPABILITY
   #define FUNCTION compute_Pantheon_LogLike
    START_FUNCTION(double)
    ALLOW_MODEL_DEPENDENCE(LCDM,cosmo_nuisance_params)
    MODEL_GROUP(cosmology, (LCDM, LCDM_Smu_dNeffCMB_dNeffBBN_etaBBN))
    MODEL_GROUP(nuisance, (cosmo_nuisance_params))
    ALLOW_MODEL_COMBINATION(cosmology,nuisance)
    BACKEND_REQ(class_get_Dl,(class_tag),double,(double))
   #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY BAO_LogLike
   START_CAPABILITY
   #define FUNCTION compute_BAO_LogLike
    START_FUNCTION(double)
    ALLOW_MODELS(LCDM, LCDM_Smu_dNeffCMB_dNeffBBN_etaBBN)
    BACKEND_REQ(class_get_Da,(class_tag),double,(double))
    BACKEND_REQ(class_get_Hz,(class_tag),double,(double))
    BACKEND_REQ(class_get_rs,(class_tag),double,())
   #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY sigma8_LogLike
     START_CAPABILITY
     #define FUNCTION compute_sigma8_LogLike
      START_FUNCTION(double)
      ALLOW_MODELS(LCDM, LCDM_Smu_dNeffCMB_dNeffBBN_etaBBN)
      BACKEND_REQ(class_get_sigma8,(class_tag),double,(double))
      BACKEND_REQ(class_get_Omega_m,(class_tag),double,())
     #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY Sigma8
     START_CAPABILITY
     #define FUNCTION compute_Sigma8
      START_FUNCTION(double)
      ALLOW_MODELS(LCDM, LCDM_Smu_dNeffCMB_dNeffBBN_etaBBN)
      BACKEND_REQ(class_get_sigma8,(class_tag),double,(double))
      BACKEND_REQ(class_get_Omega_m,(class_tag),double,())
     #undef FUNCTION
  #undef CAPABILITY

#undef MODULE

#endif /* defined __CosmoBit_rollcall_hpp__ */
