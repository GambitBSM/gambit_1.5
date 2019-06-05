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
#define BACKENDLANG Python
#define VERSION 3.1.0
#define SAFE_VERSION 3_1_0

LOAD_LIBRARY

//BE_FUNCTION(init, void, (), "init", "MontePythonLike_init")

BE_CONV_FUNCTION(test_MontePythonLike, void, (), "test_MontePythonLike")
//BE_CONV_FUNCTION(get_NR_WCs_flav, map_str_dbl, (map_str_dbl&, double&, int&, std::string&), "get_NR_WCs_flav")
//BE_CONV_FUNCTION(get_NR_WCs_EW, map_str_dbl, (map_str_dbl&, double&, double&, double&, double&, std::string&), "get_NR_WCs_EW")

// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
