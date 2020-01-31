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
///  \date 2019 June, Nov
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
    #define FUNCTION injection_spectrum_annihilatingDM
    START_FUNCTION(DarkAges::injectionSpectrum)
    ALLOW_MODELS(AnnihilatingDM_general)
    #undef FUNCTION

    #define FUNCTION injection_spectrum_decayingDM
    START_FUNCTION(DarkAges::injectionSpectrum)
    ALLOW_MODELS(DecayingDM_general)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY DM_fraction
  START_CAPABILITY
    #define FUNCTION DM_fraction_ALP
    START_FUNCTION(double)
    DEPENDENCY(T_cmb, double)
    DEPENDENCY(minimum_fraction,double)
    DEPENDENCY(lifetime,double)
    DEPENDENCY(RD_oh2, double)
    ALLOW_MODEL_DEPENDENCE(GeneralCosmoALP,LCDM)
    MODEL_GROUP(alp,(GeneralCosmoALP))
    MODEL_GROUP(cosmo,(LCDM))
    ALLOW_MODEL_COMBINATION(cosmo,alp)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY total_DM_abundance
  START_CAPABILITY
    #define FUNCTION total_DM_abundance_ALP
    START_FUNCTION(double)
    DEPENDENCY(T_cmb, double)
    DEPENDENCY(DM_fraction,double)
    ALLOW_MODEL_DEPENDENCE(GeneralCosmoALP,LCDM)
    MODEL_GROUP(alp,(GeneralCosmoALP))
    MODEL_GROUP(cosmo,(LCDM))
    ALLOW_MODEL_COMBINATION(cosmo,alp)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY lifetime
  START_CAPABILITY
    #define FUNCTION lifetime_ALP_agg
    START_FUNCTION(double)
    ALLOW_MODELS(GeneralCosmoALP)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY minimum_abundance
  START_CAPABILITY
    #define FUNCTION minimum_abundance_ALP
    START_FUNCTION(double)
    ALLOW_MODELS(GeneralCosmoALP)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY minimum_fraction
  START_CAPABILITY
    #define FUNCTION minimum_fraction_ALP
    START_FUNCTION(double)
    DEPENDENCY(T_cmb, double)
    DEPENDENCY(minimum_abundance,double)
    ALLOW_MODEL_DEPENDENCE(GeneralCosmoALP,LCDM)
    MODEL_GROUP(alp,(GeneralCosmoALP))
    MODEL_GROUP(cosmo,(LCDM))
    ALLOW_MODEL_COMBINATION(cosmo,alp)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY energy_injection_efficiency
  START_CAPABILITY
    #define FUNCTION energy_injection_efficiency_func
    START_FUNCTION(DarkAges::fz_table)
    ALLOW_MODELS(AnnihilatingDM_general,DecayingDM_general)
    BACKEND_REQ(DA_efficiency_function, (DarkAges_tag), DarkAges::fz_table,())
    #undef FUNCTION
  #undef CAPABILITY
  
  #define CAPABILITY f_effective
  START_CAPABILITY
    #define FUNCTION f_effective_func
    START_FUNCTION(double)
    ALLOW_MODELS(AnnihilatingDM_general,DecayingDM_general)
    DEPENDENCY(energy_injection_efficiency,DarkAges::fz_table)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY Planck_nuisance_prior_loglike
  START_CAPABILITY
    #define FUNCTION compute_Planck_nuisance_prior_loglike
    START_FUNCTION(double)
    ALLOW_MODELS(Planck_lite,Planck_TT,Planck_TTTEEE)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY Planck_sz_prior_loglike
  START_CAPABILITY
    #define FUNCTION compute_Planck_sz_prior
    START_FUNCTION(double)
    ALLOW_MODELS(Planck_TT,Planck_TTTEEE)
    #undef FUNCTION
  #undef CAPABILITY


  /// capabilities related to setting neutrino masses, 
  /// temperature, ncdm components & number of ultra-relativistic species Nur
  /// ----------------------
  #define CAPABILITY NuMasses_SM
  START_CAPABILITY
    #define FUNCTION set_NuMasses_SM_baseline
    START_FUNCTION(map_str_dbl)
    #undef FUNCTION

    #define FUNCTION set_NuMasses_SM
    START_FUNCTION(map_str_dbl)
    ALLOW_MODELS(StandardModel_SLHA2)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY mNu_tot
  START_CAPABILITY
    #define FUNCTION get_mNu_tot
    START_FUNCTION(double)
    DEPENDENCY(NuMasses_SM, map_str_dbl)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY N_ur
  START_CAPABILITY
    #define FUNCTION get_N_ur
    START_FUNCTION(double)
    MODEL_CONDITIONAL_DEPENDENCY(etaBBN_rBBN_rCMB_dNurBBN_dNurCMB_parameters,ModelParameters,etaBBN_rBBN_rCMB_dNurBBN_dNurCMB)
    DEPENDENCY(NuMasses_SM, map_str_dbl)
    #undef FUNCTION
  #undef CAPABILITY
 /// ------------------------
    

  /// capabilities related to setting input options for CLASS right 
  /// cosmo parameters, temperature and number of ultra-relativistic species Nur
  /// ----------------------
  #define CAPABILITY baseline_classy_input
  START_CAPABILITY
    #define FUNCTION set_baseline_classy_input
    START_FUNCTION(pybind11::dict)
    
    ALLOW_MODELS(LCDM,LCDM_no_primordial)
    //ALLOW_MODELS(LCDM)
    MODEL_CONDITIONAL_DEPENDENCY(classy_parameters_EnergyInjection, pybind11::dict, AnnihilatingDM_general, DecayingDM_general)
    MODEL_CONDITIONAL_DEPENDENCY(classy_parameters_PlanckLike, pybind11::dict, Planck_lite,Planck_TTTEEE,Planck_TT,plik_dx11dr2_HM_v18_TT)

    DEPENDENCY(T_cmb,                 double)
    DEPENDENCY(Helium_abundance,      std::vector<double>)
    DEPENDENCY(NuMasses_classy_input, pybind11::dict)
    #undef FUNCTION
  #undef CAPABILITY

    // needed to be able to initialise CLASS either with the runoptions 
    // asked for by MontePython Likelihoods (t modes, Pk at specific z,..) 
    // or without if MontePython is not in use
   #define CAPABILITY get_classy_cosmo_container
     START_CAPABILITY
     #define FUNCTION init_classy_cosmo_container_with_MPLike
      START_FUNCTION(CosmoBit::ClassyInput)
      DEPENDENCY(cosmo_args_from_MPLike,  pybind11::dict)
      DEPENDENCY(set_classy_parameters,   CosmoBit::ClassyInput)
     #undef FUNCTION
     #define FUNCTION init_classy_cosmo_container
      START_FUNCTION(CosmoBit::ClassyInput)
      DEPENDENCY(set_classy_parameters, CosmoBit::ClassyInput)
     #undef FUNCTION
  #undef CAPABILITY

    // set CLASS input parameters for ..
    #define CAPABILITY set_classy_parameters
     START_CAPABILITY
    
     // H0, tau_reio, Omega_m, Omega_b plus an external primordial power spectrum
     #define FUNCTION set_classy_parameters_primordial_ps
      START_FUNCTION(CosmoBit::ClassyInput)
         ALLOW_MODELS(LCDM_no_primordial)
         DEPENDENCY(baseline_classy_input, pybind11::dict)
         DEPENDENCY(primordial_power_spectrum, primordial_ps)
     #undef FUNCTION

     // H0, tau_reio, Omega_m, Omega_b plus an external *parametrised* primordial power spectrum
     #define FUNCTION set_classy_parameters_parametrised_ps
      START_FUNCTION(CosmoBit::ClassyInput)
         ALLOW_MODELS(LCDM_no_primordial, LCDM)
         DEPENDENCY(baseline_classy_input, pybind11::dict)
         DEPENDENCY(parametrised_power_spectrum,   parametrised_ps)
     #undef FUNCTION

  #undef CAPABILITY

  // set extra parameters for CLASS run for energy injection by non-standard cosmological models
  #define CAPABILITY classy_parameters_EnergyInjection
     START_CAPABILITY
     #define FUNCTION set_classy_parameters_EnergyInjection_AnnihilatingDM
      START_FUNCTION(pybind11::dict)
      ALLOW_MODELS(AnnihilatingDM_general)
      DEPENDENCY(energy_injection_efficiency, DarkAges::fz_table)
     #undef FUNCTION

     #define FUNCTION set_classy_parameters_EnergyInjection_DecayingDM
      START_FUNCTION(pybind11::dict)
      ALLOW_MODELS(DecayingDM_general)
      DEPENDENCY(energy_injection_efficiency, DarkAges::fz_table)
     #undef FUNCTION
  #undef CAPABILITY

  // set extra parameters for CLASS run if Planck CMB likelihoods are included
  #define CAPABILITY classy_parameters_PlanckLike
     START_CAPABILITY
     #define FUNCTION set_classy_parameters_PlanckLike
      START_FUNCTION(pybind11::dict)
      ALLOW_MODELS(Planck_lite,Planck_TTTEEE,Planck_TT,plik_dx11dr2_HM_v18_TT)
     #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY NuMasses_classy_input
  START_CAPABILITY
    #define FUNCTION set_NuMasses_classy_input
    START_FUNCTION(pybind11::dict)
    DEPENDENCY(T_ncdm,            double)
    DEPENDENCY(NuMasses_SM, map_str_dbl)
    DEPENDENCY(N_ur,double)
    #undef FUNCTION
  #undef CAPABILITY

 /// -----------

 /// primodial power spectrum related capabilities (MultiMode)
 ///
  /* MultiModeCode and power spectra */

  // Initialise settings for MultiModeCode
  #define CAPABILITY multimode_input_parameters
    START_CAPABILITY
    #define FUNCTION set_multimode_inputs
      START_FUNCTION(multimode_inputs)
      ALLOW_MODELS(Inflation_SR1quad,Inflation_1quar,Inflation_1mono23,Inflation_1linear,Inflation_1natural,Inflation_smash)
    #undef FUNCTION
  #undef CAPABILITY

  // // pass settings to multimode, run it and return the structure containing the results 
  // #define CAPABILITY multimode_results
  //   START_CAPABILITY
  //   #define FUNCTION get_multimode_results
  //     START_FUNCTION(gambit_inflation_observables)
  //     ALLOW_MODELS(Inflation_SR1quad,Inflation_1quar,Inflation_1mono23,Inflation_1linear,Inflation_1natural,Inflation_smash)
  //     DEPENDENCY(multimode_pk_setting,int)
  //     DEPENDENCY(inf_inputs, multimode_inputs)
  //     BACKEND_REQ(multimodecode_gambit_driver,(modecode_tag), void, (gambit_inflation_observables*,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,double&,int&,int&,double&,int&,double*,double*,int&,int&,double*,double*,double*,double&,double&,double&,int&,int&,double&,double*,double*,double*,double*,double&,double&))
  //   #undef FUNCTION
  // #undef CAPABILITY 

  #define CAPABILITY parametrised_power_spectrum
    START_CAPABILITY
    
    #define FUNCTION get_multimode_parametrised_ps
      START_FUNCTION(parametrised_ps)
      ALLOW_MODELS(Inflation_SR1quad,Inflation_1quar,Inflation_1mono23,Inflation_1linear,Inflation_1natural,Inflation_smash)
      //MODEL_GROUP(inflation,(Inflation_SR1quad,Inflation_1quar,Inflation_1mono23,Inflation_1linear,Inflation_1natural))
      //MODEL_GROUP(cosmo,(LCDM_no_primordial))
      //ALLOW_MODEL_COMBINATION(cosmo,inflation)
      DEPENDENCY(multimode_input_parameters, multimode_inputs)
      BACKEND_REQ(multimodecode_gambit_driver, (), gambit_inflation_observables, (int&  ,int& ,  int& , int& , double* , double* ,
																																					      	double* , double&, double&, double& , int& ,
																																				      		double& , double& , int& , int& ))
    #undef FUNCTION

    #define FUNCTION get_parametrised_ps_LCDM
      START_FUNCTION(parametrised_ps)
      ALLOW_MODELS(LCDM)
    #undef FUNCTION
    
    #define FUNCTION get_parametrised_ps_SMASH
      START_FUNCTION(parametrised_ps)
      ALLOW_MODELS(Inflation_smash)
    #undef FUNCTION

  #undef CAPABILITY 

  // pass settings to multimode, run it and return the structure containing the results 
  #define CAPABILITY primordial_power_spectrum
    START_CAPABILITY

    #define FUNCTION get_multimode_primordial_ps
      START_FUNCTION(primordial_ps)
      //ALLOW_MODELS(LCDM_no_primordial) // todo check models...
      ALLOW_MODELS(Inflation_SR1quad,Inflation_1quar,Inflation_1mono23,Inflation_1linear,Inflation_1natural,Inflation_smash)
      DEPENDENCY(multimode_input_parameters, multimode_inputs)
      BACKEND_REQ(multimodecode_gambit_driver, (), gambit_inflation_observables,(int&  ,int& ,  int& , int& , double* , double* ,
																																								 double* , double&, double&, double& , int& ,
																																								 double& , double& , int& , int& ))
    #undef FUNCTION

    /*
    #define FUNCTION get_LCDM_primordial_ps
    START_FUNCTION(primordial_ps)
    ALLOW_MODELS(LCDM)
    #undef FUNCTION
    */

    /*
    #define FUNCTION get_SMASH_primordial_ps
    START_FUNCTION(primordial_ps)
    ALLOW_MODELS(Inflation_smash)
    #undef FUNCTION
    */

  #undef CAPABILITY

  #define CAPABILITY unlensed_Cl_TT
  START_CAPABILITY
    #define FUNCTION class_get_unlensed_Cl_TT
    START_FUNCTION(std::vector<double>)
    BACKEND_REQ(class_get_unlensed_cl,(class_tag),std::vector<double>, (str))
    FORCE_SAME_BACKEND(class_tag)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY unlensed_Cl_TE
  START_CAPABILITY
    #define FUNCTION class_get_unlensed_Cl_TE
    START_FUNCTION(std::vector<double>)
    BACKEND_REQ(class_get_unlensed_cl,(class_tag),std::vector<double>, (str))
    FORCE_SAME_BACKEND(class_tag)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY unlensed_Cl_EE
  START_CAPABILITY
    #define FUNCTION class_get_unlensed_Cl_EE
    START_FUNCTION(std::vector<double>)
    BACKEND_REQ(class_get_unlensed_cl,(class_tag),std::vector<double>, (str))
    FORCE_SAME_BACKEND(class_tag)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY unlensed_Cl_BB
  START_CAPABILITY
    #define FUNCTION class_get_unlensed_Cl_BB
    START_FUNCTION(std::vector<double>)
    BACKEND_REQ(class_get_unlensed_cl,(class_tag),std::vector<double>, (str))
    FORCE_SAME_BACKEND(class_tag)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY unlensed_Cl_PhiPhi
  START_CAPABILITY
    #define FUNCTION class_get_unlensed_Cl_PhiPhi
    START_FUNCTION(std::vector<double>)
    BACKEND_REQ(class_get_unlensed_cl,(class_tag),std::vector<double>, (str))
    FORCE_SAME_BACKEND(class_tag)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY lensed_Cl_TT
  START_CAPABILITY
    #define FUNCTION class_get_lensed_Cl_TT
    START_FUNCTION(std::vector<double>)
    BACKEND_REQ(class_get_lensed_cl,(class_tag),std::vector<double>, (str))
    FORCE_SAME_BACKEND(class_tag)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY lensed_Cl_TE
  START_CAPABILITY
    #define FUNCTION class_get_lensed_Cl_TE
    START_FUNCTION(std::vector<double>)
    BACKEND_REQ(class_get_lensed_cl,(class_tag),std::vector<double>, (str))
    FORCE_SAME_BACKEND(class_tag)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY lensed_Cl_EE
  START_CAPABILITY
    #define FUNCTION class_get_lensed_Cl_EE
    START_FUNCTION(std::vector<double>)
    BACKEND_REQ(class_get_lensed_cl,(class_tag),std::vector<double>, (str))
    FORCE_SAME_BACKEND(class_tag)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY lensed_Cl_BB
  START_CAPABILITY
    #define FUNCTION class_get_lensed_Cl_BB
    START_FUNCTION(std::vector<double>)
    BACKEND_REQ(class_get_lensed_cl,(class_tag),std::vector<double>, (str))
    FORCE_SAME_BACKEND(class_tag)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY lensed_Cl_PhiPhi
  START_CAPABILITY
    #define FUNCTION class_get_lensed_Cl_PhiPhi
    START_FUNCTION(std::vector<double>)
    BACKEND_REQ(class_get_lensed_cl,(class_tag),std::vector<double>, (str))
    FORCE_SAME_BACKEND(class_tag)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY Planck_lowl_loglike
  START_CAPABILITY
    #define FUNCTION function_Planck_lowl_TT_2015_loglike
    START_FUNCTION(double)
    DEPENDENCY(lensed_Cl_TT,std::vector<double>)
    DEPENDENCY(T_cmb,double)
    ALLOW_MODELS(Planck_TTTEEE,Planck_TT,Planck_lite)
    BACKEND_REQ(plc_loglike_lowl_TT_2015,(),double,(double*))
    #undef FUNCTION

    #define FUNCTION function_Planck_lowl_TEB_2015_loglike
    START_FUNCTION(double)
    DEPENDENCY(lensed_Cl_TT,std::vector<double>)
    DEPENDENCY(lensed_Cl_TE,std::vector<double>)
    DEPENDENCY(lensed_Cl_EE,std::vector<double>)
    DEPENDENCY(lensed_Cl_BB,std::vector<double>)
    DEPENDENCY(T_cmb,double)
    ALLOW_MODELS(Planck_TTTEEE,Planck_TT,Planck_lite)
    BACKEND_REQ(plc_loglike_lowl_TEB_2015,(),double,(double*))
    #undef FUNCTION

    #define FUNCTION function_Planck_lowl_TT_2018_loglike
    START_FUNCTION(double)
    DEPENDENCY(lensed_Cl_TT,std::vector<double>)
    DEPENDENCY(T_cmb,double)
    ALLOW_MODELS(Planck_TTTEEE,Planck_TT,Planck_lite)
    BACKEND_REQ(plc_loglike_lowl_TT_2018,(),double,(double*))
    #undef FUNCTION

    #define FUNCTION function_Planck_lowl_EE_2018_loglike
    START_FUNCTION(double)
    DEPENDENCY(lensed_Cl_EE,std::vector<double>)
    DEPENDENCY(T_cmb,double)
    ALLOW_MODELS(Planck_TTTEEE,Planck_TT,Planck_lite)
    BACKEND_REQ(plc_loglike_lowl_EE_2018,(),double,(double*))
    #undef FUNCTION

    #define FUNCTION function_Planck_lowl_TTEE_2018_loglike
    START_FUNCTION(double)
    DEPENDENCY(lensed_Cl_TT,std::vector<double>)
    DEPENDENCY(lensed_Cl_EE,std::vector<double>)
    DEPENDENCY(T_cmb,double)
    ALLOW_MODELS(Planck_TTTEEE,Planck_TT,Planck_lite)
    BACKEND_REQ(plc_loglike_lowl_TT_2018,(plc_tag),double,(double*))
    BACKEND_REQ(plc_loglike_lowl_EE_2018,(plc_tag),double,(double*))
    FORCE_SAME_BACKEND(plc_tag)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY Planck_highl_loglike
  START_CAPABILITY
    #define FUNCTION function_Planck_highl_TT_2015_loglike
    START_FUNCTION(double)
    DEPENDENCY(lensed_Cl_TT,std::vector<double>)
    DEPENDENCY(T_cmb,double)
    ALLOW_MODELS(Planck_TT)
    BACKEND_REQ(plc_loglike_highl_TT_2015,(),double,(double*))
    #undef FUNCTION

    #define FUNCTION function_Planck_highl_TT_lite_2015_loglike
    START_FUNCTION(double)
    DEPENDENCY(lensed_Cl_TT,std::vector<double>)
    DEPENDENCY(T_cmb,double)
    ALLOW_MODELS(Planck_lite)
    BACKEND_REQ(plc_loglike_highl_TT_lite_2015,(),double,(double*))
    #undef FUNCTION

    #define FUNCTION function_Planck_highl_TTTEEE_2015_loglike
    START_FUNCTION(double)
    DEPENDENCY(lensed_Cl_TT,std::vector<double>)
    DEPENDENCY(lensed_Cl_TE,std::vector<double>)
    DEPENDENCY(lensed_Cl_EE,std::vector<double>)
    DEPENDENCY(T_cmb,double)
    ALLOW_MODELS(Planck_TTTEEE)
    BACKEND_REQ(plc_loglike_highl_TTTEEE_2015,(),double,(double*))
    #undef FUNCTION

    #define FUNCTION function_Planck_highl_TTTEEE_lite_2015_loglike
    START_FUNCTION(double)
    DEPENDENCY(lensed_Cl_TT,std::vector<double>)
    DEPENDENCY(lensed_Cl_TE,std::vector<double>)
    DEPENDENCY(lensed_Cl_EE,std::vector<double>)
    DEPENDENCY(T_cmb,double)
    ALLOW_MODELS(Planck_lite)
    BACKEND_REQ(plc_loglike_highl_TTTEEE_lite_2015,(),double,(double*))
    #undef FUNCTION

    #define FUNCTION function_Planck_highl_TT_2018_loglike
    START_FUNCTION(double)
    DEPENDENCY(lensed_Cl_TT,std::vector<double>)
    DEPENDENCY(T_cmb,double)
    ALLOW_MODELS(Planck_TT)
    BACKEND_REQ(plc_loglike_highl_TT_2018,(),double,(double*))
    #undef FUNCTION

    #define FUNCTION function_Planck_highl_TT_lite_2018_loglike
    START_FUNCTION(double)
    DEPENDENCY(lensed_Cl_TT,std::vector<double>)
    DEPENDENCY(T_cmb,double)
    ALLOW_MODELS(Planck_lite)
    BACKEND_REQ(plc_loglike_highl_TT_lite_2018,(),double,(double*))
    #undef FUNCTION

    #define FUNCTION function_Planck_highl_TTTEEE_2018_loglike
    START_FUNCTION(double)
    DEPENDENCY(lensed_Cl_TT,std::vector<double>)
    DEPENDENCY(lensed_Cl_TE,std::vector<double>)
    DEPENDENCY(lensed_Cl_EE,std::vector<double>)
    DEPENDENCY(T_cmb,double)
    ALLOW_MODELS(Planck_TTTEEE)
    BACKEND_REQ(plc_loglike_highl_TTTEEE_2018,(),double,(double*))
    #undef FUNCTION

    #define FUNCTION function_Planck_highl_TTTEEE_lite_2018_loglike
    START_FUNCTION(double)
    DEPENDENCY(lensed_Cl_TT,std::vector<double>)
    DEPENDENCY(lensed_Cl_TE,std::vector<double>)
    DEPENDENCY(lensed_Cl_EE,std::vector<double>)
    DEPENDENCY(T_cmb,double)
    ALLOW_MODELS(Planck_lite)
    BACKEND_REQ(plc_loglike_highl_TTTEEE_lite_2018,(),double,(double*))
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY Planck_lensing_loglike
  START_CAPABILITY
    #define FUNCTION function_Planck_lensing_2015_loglike
    START_FUNCTION(double)
    DEPENDENCY(lensed_Cl_TT,std::vector<double>)
    DEPENDENCY(lensed_Cl_TE,std::vector<double>)
    DEPENDENCY(lensed_Cl_EE,std::vector<double>)
    DEPENDENCY(lensed_Cl_PhiPhi,std::vector<double>)
    DEPENDENCY(T_cmb,double)
    ALLOW_MODELS(Planck_TTTEEE,Planck_TT,Planck_lite)
    BACKEND_REQ(plc_loglike_lensing_2015,(),double,(double*))
    #undef FUNCTION

    #define FUNCTION function_Planck_lensing_2018_loglike
    START_FUNCTION(double)
    DEPENDENCY(lensed_Cl_TT,std::vector<double>)
    DEPENDENCY(lensed_Cl_TE,std::vector<double>)
    DEPENDENCY(lensed_Cl_EE,std::vector<double>)
    DEPENDENCY(lensed_Cl_PhiPhi,std::vector<double>)
    DEPENDENCY(T_cmb,double)
    ALLOW_MODELS(Planck_TTTEEE,Planck_TT,Planck_lite)
    BACKEND_REQ(plc_loglike_lensing_2018,(),double,(double*))
    #undef FUNCTION

    #define FUNCTION function_Planck_lensing_marged_2018_loglike
    START_FUNCTION(double)
    DEPENDENCY(lensed_Cl_PhiPhi,std::vector<double>)
    ALLOW_MODELS(Planck_TTTEEE,Planck_TT,Planck_lite)
    BACKEND_REQ(plc_loglike_lensing_marged_2018,(),double,(double*))
    #undef FUNCTION
  #undef CAPABILITY


  #define CAPABILITY T_cmb
    START_CAPABILITY
    #define FUNCTION set_T_cmb
      START_FUNCTION(double)
    #undef FUNCTION
  #undef CAPABILITY
  
  // needed in addition to T_ncdm since T_ncdm of non-SM models 
  // assume a fiducial value to base calculation on 
  #define CAPABILITY T_ncdm_SM 
    START_CAPABILITY
    #define FUNCTION T_ncdm_SM
      START_FUNCTION(double)
      ALLOW_MODELS(LCDM)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY T_ncdm
    START_CAPABILITY
    #define FUNCTION set_T_ncdm
      START_FUNCTION(double)
      ALLOW_MODELS(etaBBN_rBBN_rCMB_dNurBBN_dNurCMB)
      DEPENDENCY(T_ncdm_SM,double)
    #undef FUNCTION
    
    // needed again here otherwise dependency resolution won't work
    // for SM and dependencies on T_ncdm
    #define FUNCTION set_T_ncdm_SM
      START_FUNCTION(double)
      DEPENDENCY(T_ncdm_SM, double)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY n0_g
    START_CAPABILITY

    #define FUNCTION compute_n0_g
      START_FUNCTION(double)
      DEPENDENCY(T_cmb, double)
    #undef FUNCTION

    // #define FUNCTION get_n0_g_classy
    //   START_FUNCTION(double)
    //   DEPENDENCY(get_Classy_cosmo_container, CosmoBit::Classy_cosmo_container)
    // #undef FUNCTION

  #undef CAPABILITY
  
  #define CAPABILITY Omega0_m
    START_CAPABILITY

    #define FUNCTION compute_Omega0_m
      START_FUNCTION(double)
      DEPENDENCY(Omega0_b, double)
      DEPENDENCY(Omega0_cdm, double)
      DEPENDENCY(Omega0_ncdm, double)
    #undef FUNCTION

    #define FUNCTION get_Omega0_m_classy
      START_FUNCTION(double)
      BACKEND_REQ(class_get_Omega0_m,(class_tag),double,())
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY Omega0_b
    START_CAPABILITY

    #define FUNCTION compute_Omega0_b
      START_FUNCTION(double)
      ALLOW_MODELS(LCDM)
    #undef FUNCTION

  #undef CAPABILITY

  
  #define CAPABILITY Omega0_cdm
    START_CAPABILITY

    #define FUNCTION compute_Omega0_cdm
      START_FUNCTION(double)
      ALLOW_MODELS(LCDM)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY Omega0_r
    START_CAPABILITY

    #define FUNCTION compute_Omega0_r
      START_FUNCTION(double)
      DEPENDENCY(Omega0_g, double)
      DEPENDENCY(Omega0_ur, double)
    #undef FUNCTION

    #define FUNCTION get_Omega0_r_classy
      START_FUNCTION(double)
      BACKEND_REQ(class_get_Omega0_r,(class_tag),double,())
    #undef FUNCTION

  #undef CAPABILITY


  #define CAPABILITY Omega0_g
    START_CAPABILITY

    #define FUNCTION compute_Omega0_g
      START_FUNCTION(double)
      DEPENDENCY(T_cmb, double)
      ALLOW_MODELS(LCDM)
    #undef FUNCTION

  #undef CAPABILITY
  

  #define CAPABILITY Omega0_ur
    START_CAPABILITY

    #define FUNCTION compute_Omega0_ur
      START_FUNCTION(double)
      DEPENDENCY(Omega0_g, double)
      DEPENDENCY(N_ur, double)
    #undef FUNCTION
  
    #define FUNCTION get_Omega0_ur_classy  
      START_FUNCTION(double)
      BACKEND_REQ(class_get_Omega0_ur,(class_tag),double,())
    #undef FUNCTION

  #undef CAPABILITY


  #define CAPABILITY Omega0_ncdm_tot
    START_CAPABILITY

    #define FUNCTION get_Omega0_ncdm_tot_classy 
      START_FUNCTION(double)
      BACKEND_REQ(class_get_Omega0_ncdm_tot,(class_tag),double,())
    #undef FUNCTION

  #undef CAPABILITY
  

  #define CAPABILITY Omega0_ncdm
    START_CAPABILITY

    #define FUNCTION compute_Omega0_ncdm
      START_FUNCTION(double)
      DEPENDENCY(T_cmb, double)
      DEPENDENCY(mNu_tot,double)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY eta0
    START_CAPABILITY

    #define FUNCTION calculate_eta0
      START_FUNCTION(double)
      DEPENDENCY(T_cmb, double)
      ALLOW_MODELS(LCDM)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY rs_drag
  START_CAPABILITY

    #define FUNCTION get_rs_drag_classy
      START_FUNCTION(double)
      BACKEND_REQ(class_get_rs,(class_tag),double,())
    #undef FUNCTION

  #undef CAPABILITY

  /// Good for cross-checks, innit.
  #define CAPABILITY Neff
  START_CAPABILITY

    #define FUNCTION get_Neff_classy
      START_FUNCTION(double)
      BACKEND_REQ(class_get_Neff,(class_tag),double,())
    #undef FUNCTION

  #undef CAPABILITY

/* Is there even a possible distinction between etaCMB and eta0?
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
*/

/* This capability lives now in the etaBBN_rBBN_rCMB_dNurBBN_dNurCMB model
   by using the MAP_TO_CAPABILITY macro. It might be celaner to have this definition here.

   #define CAPABILITY etaBBN
     START_CAPABILITY
     #define FUNCTION set_etaBBN // etaBBN is model parameter
       START_FUNCTION(double)
       ALLOW_MODELS(etaBBN_rBBN_rCMB_dNurBBN_dNurCMB) // To get etaCMB for etaBBN_rBBN_rCMB_dNurBBN_dNurCMB
     #undef FUNCTION
   #undef CAPABILITY
*/
  // compute dNeff AND etaBBN for non-standard models
  #define CAPABILITY external_dNeff_etaBBN
    START_CAPABILITY
    #define FUNCTION compute_dNeff_etaBBN_ALP
     START_FUNCTION(map_str_dbl)
     ALLOW_MODELS(GeneralCosmoALP)
     DEPENDENCY(total_DM_abundance, double)
     DEPENDENCY(lifetime, double)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY Sigma8
    START_CAPABILITY

    #define FUNCTION get_Sigma8_classy
      START_FUNCTION(double)
      DEPENDENCY(Omega0_m, double)
      BACKEND_REQ(class_get_sigma8,(class_tag),double,())
    #undef FUNCTION

  #undef CAPABILITY

  // AlterBBN related functions & capabilities
  #define CAPABILITY AlterBBN_setInput
    START_CAPABILITY
    #define FUNCTION AlterBBN_Input
      START_FUNCTION(map_str_dbl)
      DEPENDENCY(eta0, double)
      MODEL_CONDITIONAL_DEPENDENCY(etaBBN_rBBN_rCMB_dNurBBN_dNurCMB_parameters,ModelParameters,etaBBN_rBBN_rCMB_dNurBBN_dNurCMB)
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
    ALLOW_MODELS(LCDM)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY sigma8_LogLike
    START_CAPABILITY
    #define FUNCTION compute_sigma8_LogLike
      START_FUNCTION(double)
      DEPENDENCY(Omega0_m, double)
      BACKEND_REQ(class_get_sigma8,(class_tag),double,())
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
      START_FUNCTION(map_str_map_str_str)
      BACKEND_REQ(get_MP_availible_likelihoods, (libmontepythonlike), std::vector<str>, ())
     #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY cosmo_args_from_MPLike
     START_CAPABILITY
     #define FUNCTION init_cosmo_args_from_MPLike
      START_FUNCTION(pybind11::dict)
      DEPENDENCY(MP_experiment_names, map_str_map_str_str)
      BACKEND_REQ(create_MP_likelihood_objects, (libmontepythonlike), map_str_pyobj,    (pybind11::object&, map_str_str&))
      BACKEND_REQ(create_MP_data_object,        (libmontepythonlike), pybind11::object, (map_str_str&))
     #undef FUNCTION
  #undef CAPABILITY


  /* MontePython */
  
  /// Calculates lnL for each experiment using the experiment names
  #define CAPABILITY MP_LogLikes
    START_CAPABILITY
    #define FUNCTION calc_MP_LogLikes
      START_FUNCTION(map_str_dbl) 
      DEPENDENCY(MP_experiment_names,         map_str_map_str_str)
      DEPENDENCY(parameter_dict_for_MPLike,   pybind11::dict)
      BACKEND_REQ(get_classy_cosmo_object,               (class_tag),             pybind11::object,      ())
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
    #undef FUNCTION
  #undef CAPABILITY

  /// Calculates lnL for each experiment using the experiment names
  /// but DOES NOT add them to the lnL used to steer GAMBIT scans!
  #define CAPABILITY MP_Observables
    START_CAPABILITY
    #define FUNCTION calc_MP_observables
      START_FUNCTION(map_str_dbl) 
      DEPENDENCY(MP_experiment_names,         map_str_map_str_str)
      DEPENDENCY(parameter_dict_for_MPLike,   pybind11::dict)
      BACKEND_REQ(get_classy_cosmo_object,               (class_tag),             pybind11::object,      ())
      BACKEND_REQ(get_MP_loglike,               (libmontepythonlike), double,           (const CosmoBit::MPLike_data_container&, pybind11::object&, std::string&))
      BACKEND_REQ(create_MP_data_object,        (libmontepythonlike), pybind11::object, (map_str_str&))
      BACKEND_REQ(create_MP_likelihood_objects, (libmontepythonlike), map_str_pyobj,    (pybind11::object&, map_str_str&))
    #undef FUNCTION
  #undef CAPABILITY

#undef MODULE

#endif /* defined __CosmoBit_rollcall_hpp__ */
