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

//BE_FUNCTION(init, void, (), "init", "MontePythonLike_init")

BE_CONV_FUNCTION(test_MontePythonLike, void, (pybind11::object&), "test_MontePythonLike")

// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
