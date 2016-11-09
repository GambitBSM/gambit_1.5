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

#include <boost/preprocessor/seq/seq.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/seq/to_tuple.hpp>
#include <boost/preprocessor/tuple/to_seq.hpp>
#include <boost/preprocessor/tuple/elem.hpp>
#include <boost/preprocessor/tuple/size.hpp>

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
#define ARG_NAME(R,DATA,INDEX,ELEM) (ELEM arg##INDEX)
#define FUNCTION_ARGS_SEQ(ARGLIST) BOOST_PP_IF(ISEMPTY(ARGLIST), (),                            \
        BOOST_PP_SEQ_FOR_EACH_I(ARG_NAME, , BOOST_PP_TUPLE_TO_SEQ(ARGLIST)))
#define FUNCTION_ARGS(ARGLIST) BOOST_PP_SEQ_TO_TUPLE(FUNCTION_ARGS_SEQ(ARGLIST))

/// Macros for identifying WSTP types
#define WSPutVOID(...) 1
#define WSPutUNKNOWN(...) 0
#define WSTPTYPE(TYPE) BOOST_PP_IF(IS_TYPE(void, TYPE), VOID,                                   \
                       BOOST_PP_IF(IS_TYPE(int, TYPE), Integer32,                               \
                       BOOST_PP_IF(IS_TYPE(float, TYPE), Real32,                                \
                       BOOST_PP_IF(IS_TYPE(double, TYPE), Real64,                               \
                       BOOST_PP_IF(IS_TYPE(bool, TYPE), Integer8,                               \
                       BOOST_PP_IF(IS_TYPE(char, TYPE), Integer8,                               \
                       BOOST_PP_IF(IS_TYPE(string, TYPE), String, UNKNOWN)))))))

// Macros for stripping to basic types
#define STRIP_void void
#define STRIP_int int
#define STRIP_float float
#define STRIP_double double
#define STRIP_char char
#define STRIP_string string
#define STRIP_const DUMMY
#define STRIP_CONST(TYPE) STRIP_CONST_I(TYPE)
#define STRIP_CONST_I(TYPE) CAT(STRIP_,TYPE)

/// Macros for sending data through WSTP
#define WSPUTARG(R, DATA, INDEX, ELEM)                                                          \
  if(!CAT(WSPut,WSTPTYPE(STRIP_CONST(ELEM))) ((WSLINK)pHandle, CAT(arg,INDEX)))                 \
    backend_warning().raise(LOCAL_INFO, "Error sending packet through WSTP");

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
      TYPE NAME##_function FUNCTION_ARGS(ARGLIST)                                               \
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
        BOOST_PP_IF(ISEMPTY(ARGLIST),,                                                          \
          BOOST_PP_SEQ_FOR_EACH_I(WSPUTARG, , BOOST_PP_TUPLE_TO_SEQ(ARGLIST)))                  \
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
        /* Read the received packet into the return value, unless it's void */                  \
        BOOST_PP_IF(IS_TYPE(void, STRIP_CONST(TYPE)), DUMMY,                                    \
          if (!CAT(WSGet,WSTPTYPE(STRIP_CONST(TYPE)))((WSLINK)pHandle, &return_value))          \
            backend_warning().raise(LOCAL_INFO, "Error reading packet from WSTP");              \
                                                                                                \
          return return_value;                                                                  \
        )                                                                                       \
      }                                                                                         \
                                                                                                \
      extern const NAME##_type NAME = NAME##_function;                                          \
                                                                                                \
    }                                                                                           \
  }                                                                                             \
}

#endif // __MATHEMATICA_MACROS_HPP__
