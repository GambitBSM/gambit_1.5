//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Templated class for holding and executing
///  pointers to Python backends.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2017 Dec
///
///  *********************************************

#ifndef __python_function_hpp__
#define __python_function_hpp__

#include "gambit/Elements/ini_catch.hpp"
#include "gambit/Backends/backend_singleton.hpp"
#include "gambit/cmake/cmake_variables.hpp"

//#include some pybind header

namespace Gambit
{
  namespace Backends
  {
    // Class python_function
    template <typename TYPE>
    class python_function {};
  }

}

#endif /* __function_hpp__ */
