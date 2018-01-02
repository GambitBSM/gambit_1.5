//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Class python_variable, needed to overload
///  constructor and assignment operators to send
///  messages throught pybind11.
///
///  ***********************************************
///
///  Authors
///  =======
///
///  (add name and date if you modify)
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2017 Dec
///
///  ***********************************************

#ifndef __python_variable_hpp__
#define __python_variable_hpp__

#include "gambit/Elements/ini_catch.hpp"
#include "gambit/Backends/backend_singleton.hpp"
#include "gambit/cmake/cmake_variables.hpp"

#include <pybind11/pybind11.h>

namespace Gambit
{
  namespace Backends
  {
    // Class python_variable
    template <typename TYPE>
    class python_variable {};
  }

}

#endif /* __python_variable_hpp__ */
