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

     cout << "result1 = " << result1 << endl;
     cout << "restul2 = " << result2 << endl;

     cout << "Var = " << PrintVar() << endl;

     *Var = 121;

     cout << "Var = " << PrintVar() << endl;

     *Var = *Var + 7;

     cout << "Var = " << PrintVar() << endl;

     cout << "Var2 = " << *Var2 << endl;

     cout << "Sum of Var and Var2 = " << CalculateSum(*Var, *Var2) << endl;
     cout << PrintVarorVar2(true) << endl;
     cout << PrintVarorVar2(false) << endl;

     cout << "Var == Var2? " << VarEqualVar2() << endl;
     *Var2 = *Var;
     cout << "Var == Var2? " << VarEqualVar2() << endl;

     cout << StringTest("You") << endl;

     VoidTest();

     cout << PrintVar() << endl;
     cout << *Var << endl;

     MList<MInteger> list = {1,2,3,4,5,6};

     cout << list << endl;
     cout << ExtractElement(list, 3) << endl;
     cout << SquareList(list) << endl;

     cout << "deltavar = " << *deltavar << endl;
     *deltavar = 1001;
     cout << "deltavar = " << *deltavar << endl;
     cout << "deltavar + 2 = " << *deltavar + 2 << endl;

     cout << Gammafunc(*Var) << endl;

     return CalculateSum(result1,result2);

  }

}
END_BE_NAMESPACE

// Initialisation function (definition)
BE_INI_FUNCTION
{

}
END_BE_INI_FUNCTION
