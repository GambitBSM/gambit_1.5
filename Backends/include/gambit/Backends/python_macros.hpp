//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Macros for creating Python backend functions
///  and variables.
///
///  *********************************************
///
///  Authos (add name and date if you modify):
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2017 Dec
///
///  *********************************************

#ifndef __python_macros_hpp__
#define __python_macros_hpp__

#include "gambit/cmake/cmake_variables.hpp"

#ifdef HAVE_PYBIND11 //FIXME

  // Just use the dummy functions if pybind11 is missing
  #define BE_FUNCTION_I_PY BE_FUNCTION_I_DUMMY
  #define BE_VARIABLE_I_PY BE_VARIABLE_I_DUMMY

#else

  #include "gambit/Utils/util_macros.hpp"
  #include "gambit/Backends/python_function.hpp"
  #include "gambit/Backends/python_variable.hpp"

  /// Backend function macro for Python backends
  #define BE_FUNCTION_I_PY(NAME, TYPE, ARGLIST, SYMBOLNAME, CAPABILITY, MODELS)                 \
  namespace Gambit                                                                              \
  {                                                                                             \
    namespace Backends                                                                          \
    {                                                                                           \
      namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)                                               \
      {                                                                                         \
          /*FIXME*/                                                                             \
      }                                                                                         \
    }                                                                                           \
  }

  /// Backend variable macro for Python
  #define BE_VARIABLE_I_PY(NAME, TYPE, SYMBOLNAME, CAPABILITY, MODELS)                          \
  namespace Gambit                                                                              \
  {                                                                                             \
    namespace Backends                                                                          \
    {                                                                                           \
      namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)                                               \
      {                                                                                         \
         /*FIXME*/                                                                              \
      }                                                                                         \
    }                                                                                           \
  }

#endif // HAVE_PYBIND11

#endif // #defined __python_macros_hpp
