//  GAMBIT: Global and Modular BSM Inference Tool
//
//  *********************************************
//
//  Test function
//
//  *********************************************
//
//  Authors
//
//  Suraj Krishnamurthy
//  2017 January
//
//  *********************************************

#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/DarkBit/DarkBit_rollcall.hpp"

namespace Gambit
{
  namespace DarkBit
  {
    void lHood(double& result)
    {
      using namespace Pipes::lHood;
      double m = 0.5;
      double del = 0.01;
      double x = *Param["x"];
      double y = *Param["y"];
//      if(y<0.4 || y>0.6 || x<0.4 || x>0.6)
//      {
//        result = 2*(-pow((y-m),2)/(2*del) - pow((x-m),2)/(2*del));
//      }
//      else
//      {
        result = -pow((x-m),2)/(2*del) -pow((y-m),2)/(2*del);
//      }
    }
  }
}
