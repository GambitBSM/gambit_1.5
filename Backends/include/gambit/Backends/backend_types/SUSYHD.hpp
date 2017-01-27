//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Definitions of types
///  for the SUSYHD backend.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///   
///  \author Tomas Gonzalo
///          (t.e.gonzalo@fys.uio.no)
///  \date 2015 Apr
///
///  *********************************************

#include "gambit/Utils/util_types.hpp"

#ifndef __SUSYHD_types_hpp__
#define __SUSYHD_types_hpp__

namespace Gambit
{
  // Container for the mass of the Higgs
  struct shd_HiggsMassObs
  {
    MReal MH;
    MReal deltaMH;
  };
}

#endif
