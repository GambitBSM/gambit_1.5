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

#include "gambit/Utils/util_macros.hpp"
#ifdef HAVE_PYBIND11
  #include "pybind11/pybind11.h"
#endif

/// Backend function macro for Python backends
#ifndef HAVE_PYBIND11 //FIXME
  #define BE_FUNCTION_I_PY(NAME, TYPE, ARGLIST, SYMBOLNAME, CAPABILITY, MODELS)                 \
  namespace Gambit                                                                              \
  {                                                                                             \
    namespace Backends                                                                          \
    {                                                                                           \
      namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)                                               \
      {                                                                                         \
                                                                                                \
      }                                                                                         \
    }                                                                                           \
  }
#else
  /* If Python is not available in the system return a dummy funtion */
  #define BE_FUNCTION_I_PY(NAME, TYPE, ARGLIST, SYMBOLNAME, CAPABILITY, MODELS)                 \
  namespace Gambit                                                                              \
  {                                                                                             \
    namespace Backends                                                                          \
    {                                                                                           \
      namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)                                               \
      {                                                                                         \
                                                                                                \
        /* Define a type NAME_type to be a suitable function pointer. */                        \
        typedef TYPE (*NAME##_type) CONVERT_VARIADIC_ARG(ARGLIST);                              \
                                                                                                \
        TYPE NAME##_function FUNCTION_ARGS(ARGLIST)                                             \
        {                                                                                       \
                                                                                                \
          /* Do something inconsequential with the args to skip compiler warnings. */           \
          BOOST_PP_IF(ISEMPTY(ARGLIST),,                                                        \
            BOOST_PP_SEQ_FOR_EACH_I(VOIDARG, , BOOST_PP_TUPLE_TO_SEQ(ARGLIST)))                 \
                                                                                                \
          /* Return a dummy value, unless the function type is void */                          \
          BOOST_PP_IF(IS_TYPE(void, STRIP(TYPE)), DUMMY, return new TYPE();)                    \
        }      /* FIXME Should this be new?? */                                                 \
                                                                                                \
        extern const NAME##_type NAME = NAME##_function;                                        \
                                                                                                \
      }                                                                                         \
    }                                                                                           \
  }
#endif

/// Backend variable macro for Python
#ifndef HAVE_PYBIND11 //FIXME
  #define BE_VARIABLE_I_MATH(NAME, TYPE, SYMBOLNAME, CAPABILITY, MODELS)                        \
  namespace Gambit                                                                              \
  {                                                                                             \
    namespace Backends                                                                          \
    {                                                                                           \
      namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)                                               \
      {                                                                                         \
                                                                                                \
      }                                                                                         \
    }                                                                                           \
  }
#else
  /* If pybind11 is not available, define a dummy variable */
  #define BE_VARIABLE_I_PY(NAME, TYPE, SYMBOLNAME, CAPABILITY, MODELS)                          \
  namespace Gambit                                                                              \
  {                                                                                             \
    namespace Backends                                                                          \
    {                                                                                           \
      namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)                                               \
      {                                                                                         \
                                                                                                \
        extern TYPE* const NAME = new TYPE();                                                   \
                                                                                                \
        TYPE* CAT(getptr,NAME)() { return NAME; }                                               \
                                                                                                \
      }                                                                                         \
    }                                                                                           \
  }
#endif


namespace Gambit
{
  namespace Backends
  {
    // Class python_variable
    template <typename TYPE>
    class python_variable {};
  }

}


#endif // #defined __python_macros_hpp
