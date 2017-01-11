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

BE_FUNCTION(CalculateSquare, MReal, (const MInteger&), "CalculateSquare","MathematicaTest")
BE_FUNCTION(CalculateSum, MReal, (const MReal&, const MReal&), "CalculateSum", "MathematicaTest")
BE_FUNCTION(PrintVar, MInteger, (), "PrintVar", "MathematicaTest")
BE_FUNCTION(PrintVarorVar2, MReal, (const MBool&), "PrintVarorVar2", "MathematicaTest")
BE_FUNCTION(VarEqualVar2, MBool, (), "VarEqualVar2", "MathematicaTest")
BE_FUNCTION(StringTest, MString, (const MString&), "StringTest", "MathematicaTest")
BE_FUNCTION(VoidTest, MVoid, (), "VoidTest", "MathematicaTest")
BE_FUNCTION(ExtractElement, MInteger, (const MList<MInteger>&, const MInteger&), "ExtractElement", "MathematicaTest")
BE_FUNCTION(SquareList, MList<MInteger>, (const MList<MInteger>&), "SquareList", "MathematicaTest")

BE_VARIABLE(Var, MInteger, "Var", "MathematicaTest")
BE_VARIABLE(Var2, MReal, "Var2", "MathematicaTest")

/* Convenience functions (declarations) */
BE_CONV_FUNCTION(run_Mathematica_Test, double, (const int&, const int&), "MathematicaTest")

// Undefine macros toa void conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
