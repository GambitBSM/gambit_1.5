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

#define BACKENDNAME classy
#define BACKENDLANG Python2
#define VERSION 2.6.3
#define SAFE_VERSION 2_6_3

LOAD_LIBRARY

//BE_FUNCTION(init, void, (), "init", "MontePythonLike_init")

BE_CONV_FUNCTION(classy_2_6_3_create_python_obj, void, (pybind11::object&), "classy_2_6_3_create_python_obj")
BE_CONV_FUNCTION(classy_2_6_3_set_parameter, void, (pybind11::object&,pybind11::dict&), "classy_2_6_3_set_parameter")
//BE_CONV_FUNCTION(classy_2_6_3_set_parameter, void, (pybind11::dict, pybind11::dict), "classy_2_6_3_set_parameter")

// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
