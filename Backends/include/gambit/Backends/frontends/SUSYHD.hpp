//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for SUSY.
///
///  *********************************************
///
///  Authors (add name and sate if you modify):
///
///  \author Tomas Gonzalo
///  \date 2016 Dec
///
///  *********************************************

#define BACKENDNAME SUSYHD
#define BACKENDLANG MATHEMATICA
#define VERSION 1.0.2
#define SAFE_VERSION 1_0_2

LOAD_LIBRARY

//BE_FUNCTION(MHiggs, double, (MList), "MHiggs", "")
/*BE_FUNCTION(DeltaMHiggs, double, (MList), "\[Delta]MHiggs", "")
BE_FUNCTION(SetSMParameters, void, (MReal, MReal), "SetSMParameters", "")
*/
/* Convenience functions (declarations) */

BE_INI_FUNCTION {}
END_BE_INI_FUNCTION

// Undefine macros toa void conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
