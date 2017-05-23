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

BE_ALLOW_MODELS(CMSSM)

BE_FUNCTION(Start, void, (const MString&), "Start", "SARAH_Start")
//BE_FUNCTION(Make_SPheno, void, (), "Make_SPheno", "SARAH_Make_SPheno")

/* Convenience functions (declarations) */

// Undefine macros toa void conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
