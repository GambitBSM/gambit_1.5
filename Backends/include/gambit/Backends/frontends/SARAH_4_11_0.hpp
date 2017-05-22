//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for the SARAH backend.
///
///  *********************************************
///
///  Authors (add name and sate if you modify):
///
///  \author Tomas Gonzalo
///  \date 2017 Apr
///
///  *********************************************

#define BACKENDNAME SARAH
#define BACKENDLANG MATHEMATICA
#define VERSION 4.11.0
#define SAFE_VERSION 4_11_0

LOAD_LIBRARY

/* Convenience functions (declarations) */

// Undefine macros toa void conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
