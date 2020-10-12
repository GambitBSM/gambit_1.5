///  Rollcall header for the NuclearBit module.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Anders Kvellestad
///          (anders.kvellestad@fys.uio.no)
///  \date 2020 Sep
///
///  *********************************************

#ifndef __NuclearBit_rollcall_hpp__
#define __NuclearBit_rollcall_hpp__

#include "gambit/cmake/cmake_variables.hpp"
#include "gambit/Elements/shared_types.hpp"


#define MODULE NuclearBit
START_MODULE

  // Run the gledeli backend and get all results in a map<str,double>
  #define CAPABILITY gledeliResults
  START_CAPABILITY
    #define FUNCTION getGledeliResults
    START_FUNCTION(map_str_dbl)
    ALLOW_MODELS(GSFModel2, GSFModel3, GSFModel4, GSFModel5, GSFModel10, GSFModel20, NLDModelCT_and_discretes, NLDModel10, NLDModel15, NLDModel20)
    BACKEND_REQ(gledeliBE_set_model_pars, (), void, (pybind11::dict&))
    BACKEND_REQ(gledeliBE_set_model_names, (), void, (pybind11::list&))
    BACKEND_REQ(gledeliBE_run, (), pybind11::dict, (pybind11::dict&))
    BACKEND_REQ(gledeliBE_get_results, (), pybind11::dict, ())
    #undef FUNCTION
  #undef CAPABILITY

  // Get the total loglike returned by the gledeli backend.
  // Depends on the result from gledeli (capability gledeliResults above).
  #define CAPABILITY gledeliLogLike
  START_CAPABILITY
    #define FUNCTION getGledeliLogLike
    START_FUNCTION(double)
    DEPENDENCY(gledeliResults, map_str_dbl)
    #undef FUNCTION
  #undef CAPABILITY


#undef MODULE

#endif /* defined __NuclearBit_rollcall_hpp__ */
