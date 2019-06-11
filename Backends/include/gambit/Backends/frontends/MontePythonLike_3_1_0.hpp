//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend header for the DirectDM backend.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2019 June
///
///  *********************************************

#define BACKENDNAME MontePythonLike
#define BACKENDLANG Python2
#define VERSION 3.1.0
#define SAFE_VERSION 3_1_0

LOAD_LIBRARY

BE_CONV_FUNCTION(get_MP_loglike, double, (const CosmoBit::Classy_cosmo_container&, pybind11::object&), "get_MP_loglike")

BE_CONV_FUNCTION(create_data_object, pybind11::object, (), "create_data_object")
BE_CONV_FUNCTION(create_likelihood_objects, map_str_dbl, (pybind11::object &), "create_likelihood_objects")

//BE_INI_DEPENDENCY(classy_python_obj,CosmoBit::Classy_cosmo_container)

// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
