//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Macros for creating mathematica functions and
///  sending and receiving packets through WSTP
///
///  *********************************************
///
///  Authos (add name and date if you modify):
///
///  \author Tomas Gonzalo
///          (t.e.gonzalo@fys.uio.no)
///  \date 2016 Oct
///
///  *********************************************

#ifndef __MATHEMATICA_MACROS_HPP__
#define __MATHEMATICA_MACROS_HPP__

#include "gambit/Utils/util_macros.hpp"
#include "gambit/Backends/ini_functions.hpp"

/// If not defined already, define Mathematica
#ifndef MATHEMATICA
#define MATHEMATICA 3
#endif

/// Macro to help identifying the language of the backend
#ifndef DEFINED_BACKENDLANG 
#define DEFINED_BACKENDLANG ()
#endif

/// Macro that determines whether the language of the backend is mathematica
#define USING_MATHEMATICA IF_ELSE_TOKEN_DEFINED(BACKENDLANG,                                    \
        BOOST_PP_EQUAL(BACKENDLANG, MATHEMATICA), 0)

/// Macro for testing stuff
#define TEST(NAME, STUFF) int NAME##_stuff = print_stuff(STRINGIFY(STUFF));

/// Macros to give names to an argument list
#define ARGS_WITH_NAMES(ARGLIST) CONVERT_VARIADIC_ARG(ARGLIST)

/// Backend function macro for mathematica
#define BE_FUNCTION_I_MATH(NAME, TYPE, ARGLIST, SYMBOLNAME, CAPABILITY, MODELS)                 \
namespace Gambit                                                                                \
{                                                                                               \
  namespace Backends                                                                            \
  {                                                                                             \
    namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)                                                 \
    {                                                                                           \
                                                                                                \
      /* Define a type NAME_type to be a suitable function pointer. */                          \
      typedef TYPE (*NAME##_type) CONVERT_VARIADIC_ARG(ARGLIST);                                \
                                                                                                \
      TYPE NAME_function (const int& val)                                                       \
      {                                                                                         \
        TYPE return_value;                                                                      \
                                                                                                \
        if(!WSPutFunction((WSLINK)pHandle, "CalculateSquare", 1)                                \
            or !WSPutInteger((WSLINK)pHandle, 7)                                                \
            or !WSEndPacket((WSLINK)pHandle))                                                   \
          cout << "Error sending packet" << endl;                                               \
                                                                                                \
        int pkt;                                                                                \
        while( (pkt = WSNextPacket((WSLINK)pHandle), pkt) && pkt != RETURNPKT)                  \
        {                                                                                       \
          WSNewPacket((WSLINK)pHandle);                                                         \
          if (WSError((WSLINK)pHandle)) cout << "Error reading package" << endl;                \
        }                                                                                       \
                                                                                                \
        double square;                                                                          \
        if (WSGetReal((WSLINK)pHandle, &square))                                                \
        {                                                                                       \
          cout << "Calculate square of " << 7 << " is " << square << endl;                      \
        }                                                                                       \
        else                                                                                    \
        {                                                                                       \
          cout << "Error" << endl;                                                              \
        }                                                                                       \
                                                                                                \
        return return_value;                                                                    \
      }                                                                                         \
                                                                                                \
      extern const NAME##_type NAME = NAME_function;                                            \
                                                                                                \
    }                                                                                           \
  }                                                                                             \
}

#endif // __MATHEMATICA_MACROS_HPP__
