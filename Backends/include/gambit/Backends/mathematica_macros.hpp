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
#define TEST(NAME, STUFF) int NAME##_stuff = print_stuff(STUFF);

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
      TYPE NAME##_function (const int& val)                                                     \
      {                                                                                         \
        TYPE return_value;                                                                      \
                                                                                                \
        /* If TYPE is a numeric type, send N first */                                           \
        if(IS_TYPE(TYPE, int) or IS_TYPE(TYPE, float) or IS_TYPE(TYPE, double))                 \
          if(!WSPutFunction((WSLINK)pHandle, "N", 1))                                           \
            backend_warning().raise(LOCAL_INFO, "Error sending packet through WSTP");           \
                                                                                                \
        /* Send the symbol name next */                                                         \
        if(!WSPutFunction((WSLINK)pHandle, SYMBOLNAME, 1))                                      \
          backend_warning().raise(LOCAL_INFO, "Error sending packet through WSTP");             \
                                                                                                \
        /* Now send all the arguments */                                                        \
        if(!WSPutInteger((WSLINK)pHandle, val))                                                 \
          backend_warning().raise(LOCAL_INFO, "Error sending packet throug WSTP");              \
                                                                                                \
        /* Last, mark the end of the message */                                                 \
        if(!WSEndPacket((WSLINK)pHandle))                                                       \
          backend_warning().raise(LOCAL_INFO, "Error sending packet through WSTP");             \
                                                                                                \
        /* Wait to receive a packet from the kernel */                                          \
        int pkt;                                                                                \
        while( (pkt = WSNextPacket((WSLINK)pHandle), pkt) && pkt != RETURNPKT)                  \
        {                                                                                       \
          WSNewPacket((WSLINK)pHandle);                                                         \
          if (WSError((WSLINK)pHandle))                                                         \
            backend_warning().raise(LOCAL_INFO, "Error reading packet from WSTP");              \
        }                                                                                       \
                                                                                                \
        /* Read the received packet into the return value */                                    \
        if (!WSGetReal((WSLINK)pHandle, &return_value))                                         \
          backend_warning().raise(LOCAL_INFO, "Error reading packet from WSTP");                \
                                                                                                \
        return return_value;                                                                    \
      }                                                                                         \
                                                                                                \
      extern const NAME##_type NAME = NAME##_function;                                          \
      TYPE NAME##_thingy = (*NAME)(12);                                                         \
      TEST(NAME, NAME##_thingy)                                                                 \
                                                                                                \
    }                                                                                           \
  }                                                                                             \
}

#endif // __MATHEMATICA_MACROS_HPP__
