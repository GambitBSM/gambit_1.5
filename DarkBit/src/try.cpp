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
    double getSE(double arr[])
    {
     int i, size = 40000;
      double mean = 0.5;
      double SE = 0, SE_num = 0;
      for(i=0; i<size; i++)
      {
          SE_num += pow((arr[i]-mean),2);
      }
      SE = SE_num/(size*(size-1));
      return SE;
    }
    void lHood(double& result)
    {
      using namespace Pipes::lHood;
      double x = *Param["x"];
      double y = *Param["y"];
      double y_arr[];
      double m = 0.5;
      int j, s = 40000;
      for(j=0; j<s; j++)
      {
          y_arr[j] = y;
      }
      double SE;
      SE = getSE(y_arr);
//      ln L = (-(y[i]-y0)^2)/(2*(del_y[i]^2))
//      ln L = (-(y[i]-y0)^2)/(2*(del_y[i]^2)) + step function
      if(y<m)
      {
        result = 0.0;
      }
      else
      {
        result = -pow((y-m),2)/(2*SE);
      }
    }
  }
}
