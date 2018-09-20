//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Rollcall declarations for module functions
///  contained in SpecBit_VS.cpp
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author James McKay
///          (j.mckay14@imperial.ac.uk)
///    \date Nov 2015
///
///
///  *********************************************

#ifndef __SpecBit_VS_rollcall_hpp__
#define __SpecBit_VS_rollcall_hpp__

  #define CAPABILITY check_EW_stability_ScalarSingletDM_Z3
  START_CAPABILITY
    #define FUNCTION check_EW_stability_ScalarSingletDM_Z3
    START_FUNCTION(double)
    DEPENDENCY(ScalarSingletDM_Z3_spectrum, Spectrum)
    ALLOW_MODEL_DEPENDENCE(ScalarSingletDM_Z3)
    MODEL_GROUP(higgs,   (StandardModel_Higgs_running))
    MODEL_GROUP(singlet, (ScalarSingletDM_Z3_running))
    ALLOW_MODEL_COMBINATION(higgs, singlet)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY vacuum_stability
  START_CAPABILITY

    #define FUNCTION find_min_lambda_ScalarSingletDM_Z2
    START_FUNCTION(dbl_dbl_bool)
    DEPENDENCY(SMINPUTS, SMInputs)
    DEPENDENCY(ScalarSingletDM_Z2_spectrum, Spectrum)
    ALLOW_MODEL_DEPENDENCE(StandardModel_Higgs_running, ScalarSingletDM_Z2_running)
    MODEL_GROUP(higgs,   (StandardModel_Higgs_running))
    MODEL_GROUP(singlet, (ScalarSingletDM_Z2_running))
    ALLOW_MODEL_COMBINATION(higgs, singlet)
    #undef FUNCTION

		#define FUNCTION find_min_lambda_ScalarSingletDM_Z3
    START_FUNCTION(dbl_dbl_bool)
    DEPENDENCY(SMINPUTS, SMInputs)
    DEPENDENCY(ScalarSingletDM_Z3_spectrum, Spectrum)
    ALLOW_MODEL_DEPENDENCE(StandardModel_Higgs_running, ScalarSingletDM_Z3)
    MODEL_GROUP(higgs,   (StandardModel_Higgs_running))
    MODEL_GROUP(singlet, (ScalarSingletDM_Z3_running))
    ALLOW_MODEL_COMBINATION(higgs, singlet)
    #undef FUNCTION

	  #define FUNCTION find_min_lambda_MDM
    START_FUNCTION(dbl_dbl_bool)
    DEPENDENCY(SMINPUTS, SMInputs)
    DEPENDENCY(MDM_spectrum, Spectrum)
    ALLOW_MODEL_DEPENDENCE(StandardModel_Higgs_running, MDM)
    MODEL_GROUP(higgs,   (StandardModel_Higgs_running))
    MODEL_GROUP(mdm, (MDM))
    ALLOW_MODEL_COMBINATION(higgs, mdm)
    #undef FUNCTION

  #undef CAPABILITY


  #define CAPABILITY VS_likelihood
  START_CAPABILITY

    #define FUNCTION get_likelihood
    START_FUNCTION(double)
    DEPENDENCY(vacuum_stability, dbl_dbl_bool)
    #undef FUNCTION

  #undef CAPABILITY


  #define CAPABILITY expected_lifetime
    START_CAPABILITY
    #define FUNCTION get_expected_lifetime
    START_FUNCTION(double)
    DEPENDENCY(vacuum_stability, dbl_dbl_bool)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY check_perturb_min_lambda
    START_CAPABILITY
    #define FUNCTION get_check_perturb_min_lambda
    START_FUNCTION(double)
    DEPENDENCY(vacuum_stability, dbl_dbl_bool)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY get_lambdaB
    START_CAPABILITY
    #define FUNCTION get_lambdaB
    START_FUNCTION(double)
    DEPENDENCY(vacuum_stability, dbl_dbl_bool)
    #undef FUNCTION
  #undef CAPABILITY


#endif

