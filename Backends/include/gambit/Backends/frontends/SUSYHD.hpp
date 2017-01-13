//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for SUSYHD.
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

BE_FUNCTION(MHiggs, MReal, (const MList<MReal>&), "MHiggs", "SUSYHD_mh")
BE_FUNCTION(DeltaMHiggs, MReal, (MList<MReal>&), "\[Delta]MHiggs", "SUSYHD_mh")
BE_FUNCTION(SetSMparameters, MVoid, (const MReal&, const MReal&), "SetSMparameters", "SUSYHD_mh")

/* Convenience functions (declarations) */

BE_INI_FUNCTION {}
END_BE_INI_FUNCTION

// Undefine macros toa void conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
