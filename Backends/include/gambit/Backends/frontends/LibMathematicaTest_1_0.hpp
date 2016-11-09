//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  A Mathematica backend example.
///
///  *********************************************
///
///  Authors (add name and sate if you modify):
///
///  \author Tomas Gonzalo
///  \date 2016 Sept
///
///  *********************************************

#define BACKENDNAME libMathematicaTest
#define BACKENDLANG MATHEMATICA
#define VERSION 1.0
#define SAFE_VERSION 1_0

LOAD_LIBRARY

BE_ALLOW_MODELS(CMSSM)

BE_FUNCTION(CalculateSquare, double, (const int&), "CalculateSquare","MathematicaTest")

//BE_VARIABLE(Number, int, "Number", "MathematicaTest")

namespace Gambit
{
  namespace Backends
  {
    namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)
    {

      /* Convenience functions go here */
      //BE_CONV_FUNCTION(run_Mathematica_Test, double, (const int&), "MathematicaTest")

    } /* end namespace BACKENDNAME_SAFE_VERSION */
  } /*end namespace Backends */
} /* end namespace Gambit */

BE_INI_FUNCTION{}
END_BE_INI_FUNCTION

// Undefine macros toa void conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
