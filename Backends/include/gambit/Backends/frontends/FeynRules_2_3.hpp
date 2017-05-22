//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for FeynRules 2.3.
///
///  *********************************************
///
///  Authors (add name and sate if you modify):
///
///  \author Tomas Gonzalo
///  \date 2017 Mar
///
///  *********************************************

#define BACKENDNAME FeynRules
#define BACKENDLANG MATHEMATICA
#define VERSION 2.3
#define SAFE_VERSION 2_3

LOAD_LIBRARY


/* Convenience functions (declarations) */

// Undefine macros toa void conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
