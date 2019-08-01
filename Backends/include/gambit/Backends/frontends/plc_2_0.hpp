//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend header for the plc backend.
///
///  Compile-time registration of available
///  functions and variables from this backend.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Selim C. Hotinli
///          (selim.hotinli14@pimperial.ac.uk)
///  \date 2016 Sep
///
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2017 Nov
///  \date 2019 July, Aug
///
///  *********************************************

#define BACKENDNAME plc
#define BACKENDLANG CC
#define VERSION 2.0
#define SAFE_VERSION 2_0

#include "gambit/Utils/util_functions.hpp"

LOAD_LIBRARY

// Internal function of plc (Not intended to be called directly)
BE_FUNCTION(initError, clik_error* , (),"initError","clik_initialize_error")
BE_FUNCTION(isError, int , (clik_error*),"_isError","clik_is_error")
BE_FUNCTION(stringError, void , (char*, clik_error*),"stringError","clik_string_error")
BE_FUNCTION(clik_init, clik_object*, (char*,clik_error**),"clik_init","clik_initialize")
BE_FUNCTION(clik_lensing_init, clik_lensing_object*, (char*,clik_error**),"clik_lensing_init","clik_lensing_initialize")
BE_FUNCTION(clik_compute, double, (clik_object*,double*,clik_error**), "clik_compute","clik_compute_loglike")
BE_FUNCTION(clik_lensing_compute, double, (clik_lensing_object*,double*,clik_error**), "clik_lensing_compute","clik_lensing_compute_loglike")

// All relevant data and variables will be kept within the fronted.
// Define convenience functions for the communication with the outside world
BE_CONV_FUNCTION(plc_loglike_highl_TTTEEE_2015,double,(double*),"plc_loglike_highl_TTTEEE_2015",(Planck_TTTEEE))
BE_CONV_FUNCTION(plc_loglike_highl_TT_2015,double,(double*),"plc_loglike_highl_TT_2015",(Planck_TT))
BE_CONV_FUNCTION(plc_loglike_highl_TTTEEE_lite_2015,double,(double*),"plc_loglike_highl_TTTEEE_lite_2015",(Planck_lite))
BE_CONV_FUNCTION(plc_loglike_highl_TT_lite_2015,double,(double*),"plc_loglike_highl_TT_lite_2015",(Planck_lite))
BE_CONV_FUNCTION(plc_loglike_lowl_TEB_2015,double,(double*),"plc_loglike_lowl_TEB_2015",(Planck_TTTEEE, Planck_TT,Planck_lite))
BE_CONV_FUNCTION(plc_loglike_lowl_TT_2015,double,(double*),"plc_loglike_lowl_TT_2015",(Planck_TTTEEE, Planck_TT,Planck_lite))
BE_CONV_FUNCTION(plc_loglike_lensing_2015,double,(double*),"plc_loglike_lensing_2015",(Planck_TTTEEE, Planck_TT, Planck_lite))

// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
