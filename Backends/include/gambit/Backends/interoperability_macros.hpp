//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Macros for interfacing with backends written
///  in other languages.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2017 Dec
///
///  \author Tomas Gonzalo
///          (t.e.gonzalo@fys.uio.no)
///  \date 2016 Oct
///
///  *********************************************

#ifndef __interoperability_macros_hpp__
#define __interoperability_macros_hpp__

/// If not defined already, define the backend languages
#ifndef UNKNOWN
  #define UNKNOWN 0
#endif
#ifndef CC
  #define CC 1
#endif
#ifndef CXX
  #define CXX 2
#endif
#ifndef FORTRAN
  #define FORTRAN 3
#endif
#ifndef MATHEMATICA
  #define MATHEMATICA 4
#endif
#ifndef PYTHON
  #define PYTHON 5
#endif

/// Macro to help identifying the language of the backend
#ifndef DEFINED_BACKENDLANG
#define DEFINED_BACKENDLANG ()
#endif

/// Macro to choose between mathematica types, python types and normal types
#ifdef HAVE_MATHEMATICA
  #ifdef HAVE_PYBIND11
    #define MATH_TYPE(TYPE)                                                                     \
     IF_ELSEIF(USING_MATHEMATICA, mathematica_variable<TYPE>,                                   \
               USING_PYTHON, python_variable<TYPE>,                                             \
               /*USING NONE OF THE ABOVE*/ TYPE)
  #else
    #define MATH_TYPE(TYPE) BOOST_PP_IF(USING_MATHEMATICA,mathematica_variable<TYPE>,TYPE)
  #endif
#else
  #ifdef HAVE_PYBIND11
    #define MATH_TYPE(TYPE) BOOST_PP_IF(USING_PYTHON,python_variable<TYPE>,TYPE)
  #else
    #define MATH_TYPE(TYPE) TYPE
  #endif
#endif

/// Macro that determines whether the language of the backend is C
#define USING_CC IF_ELSE_TOKEN_DEFINED(BACKENDLANG,                                             \
        BOOST_PP_EQUAL(BACKENDLANG, CC), 0)

/// Macro that determines whether the language of the backend is C++
#define USING_CXX IF_ELSE_TOKEN_DEFINED(BACKENDLANG,                                            \
        BOOST_PP_EQUAL(BACKENDLANG, CXX), 0)

/// Macro that determines whether the language of the backend is Fortran
#define USING_FORTRAN IF_ELSE_TOKEN_DEFINED(BACKENDLANG,                                        \
        BOOST_PP_EQUAL(BACKENDLANG, FORTRAN), 0)

/// Macro that determines whether the language of the backend is Mathematica
#define USING_MATHEMATICA IF_ELSE_TOKEN_DEFINED(BACKENDLANG,                                    \
        BOOST_PP_EQUAL(BACKENDLANG, MATHEMATICA), 0)

/// Macro that determines whether the language of the backend is Mathematica
#define USING_PYTHON IF_ELSE_TOKEN_DEFINED(BACKENDLANG,                                         \
        BOOST_PP_EQUAL(BACKENDLANG, PYTHON), 0)

#endif // #defined __interoperability_macros_hpp__