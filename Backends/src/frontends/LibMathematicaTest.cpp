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
///  \date 2016 Nov
///
///  *********************************************

#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/LibMathematicaTest_1_0.hpp"

// Convenience functions (definitions)
BE_NAMESPACE
{

  double run_Mathematica_Test(const int& val1, const int& val2)
  {

     double result1 = CalculateSquare(val1);
     double result2 = CalculateSquare(val2);

     cout << "Var = " << *Var << endl;
     //cout << "Var = " << PrintVar() << endl;

//     *Var = 121;

  //   cout << "Var = " << PrintVar() << endl;

     return result1 + result2;

  }

}
END_BE_NAMESPACE

// Initialisation function (definition)
BE_INI_FUNCTION
{

}
END_BE_INI_FUNCTION
