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

BE_CONV_FUNCTION(get_MP_loglike, double, (std::string), "get_MP_loglike")

BE_INI_DEPENDENCY(classy_python_obj,CosmoBit::Classy_cosmo_container)

// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
