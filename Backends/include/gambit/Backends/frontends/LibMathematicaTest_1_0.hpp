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
BE_FUNCTION(CalculateSum, double, (const double&, const double&), "CalculateSum", "MathematicaTest")
BE_FUNCTION(PrintVar, int, (), "PrintVar", "MathematicaTest")
BE_FUNCTION(PrintVarorVar2, double, (const bool&), "PrintVarorVar2", "MathematicaTest")
BE_FUNCTION(VarEqualVar2, bool, (), "VarEqualVar2", "MathematicaTest")
BE_FUNCTION(StringTest, str, (const str&), "StringTest", "MathematicaTest")
BE_FUNCTION(VoidTest, void, (), "VoidTest", "MathematicaTest")

BE_VARIABLE(Var, int, "Var", "MathematicaTest")
BE_VARIABLE(Var2, double, "Var2", "MathematicaTest")

/* Convenience functions (declarations) */
BE_CONV_FUNCTION(run_Mathematica_Test, double, (const int&, const int&), "MathematicaTest")

// Undefine macros toa void conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
