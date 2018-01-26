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

#ifndef __SpecBit_SingletDMZ3_hpp__
#define __SpecBit_SingletDMZ3_hpp__

  #define CAPABILITY check_EW_stability_SingletDMZ3
  START_CAPABILITY
    #define FUNCTION check_EW_stability_SingletDMZ3
    START_FUNCTION(double)
    DEPENDENCY(SingletDMZ3_spectrum, Spectrum)
    ALLOW_MODEL_DEPENDENCE(SingletDMZ3)
    MODEL_GROUP(higgs,   (StandardModel_Higgs_running))
    MODEL_GROUP(singlet, (SingletDMZ3))
    ALLOW_MODEL_COMBINATION(higgs, singlet)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY vacuum_stability
  START_CAPABILITY

    #define FUNCTION find_min_lambda_SingletDM
    START_FUNCTION(dbl_dbl_bool)
    DEPENDENCY(SMINPUTS, SMInputs)
    DEPENDENCY(SingletDM_spectrum, Spectrum)
    ALLOW_MODEL_DEPENDENCE(StandardModel_Higgs_running, SingletDM_running)
    MODEL_GROUP(higgs,   (StandardModel_Higgs_running))
    MODEL_GROUP(singlet, (SingletDM_running,SingletDMZ3))
    ALLOW_MODEL_COMBINATION(higgs, singlet)
    #undef FUNCTION

		#define FUNCTION find_min_lambda_SingletDMZ3
    START_FUNCTION(dbl_dbl_bool)
    DEPENDENCY(SMINPUTS, SMInputs)
    DEPENDENCY(SingletDMZ3_spectrum, Spectrum)
    ALLOW_MODEL_DEPENDENCE(StandardModel_Higgs_running, SingletDMZ3)
    MODEL_GROUP(higgs,   (StandardModel_Higgs_running))
    MODEL_GROUP(singlet, (SingletDMZ3))
    ALLOW_MODEL_COMBINATION(higgs, singlet)
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

