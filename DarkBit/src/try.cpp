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
      double y_mean = 0.5;
      //double x = *Param["x"];
      double y = *Param["y"];
      //ln L = (-(y[i]-y0)^2)/(2*(del_y[i]^2))
      //ln L = (-(y[i]-y0)^2)/(2*(del_y[i]^2)) + step function
      result = pow((y-y_mean),2);
    }
  }
}
