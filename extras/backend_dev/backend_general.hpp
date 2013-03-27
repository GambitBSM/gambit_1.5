/*
 * General macros for loading a shared library
 * and constructing pointers to the variables and functions 
 * within the library. 
 * 
 * \author Anders Kvellestad
 * \date 2013-03-26
 * 
 */

#ifndef __BACKEND_GENERAL_HPP__
#define __BACKEND_GENERAL_HPP__

#include <iostream>
#include <string>
#include "dlfcn.h"


// A container struct for code that needs to be executed as 
// initialization code at startup. When an instance of the 'ini_code'
// struct is created, the constructor will execute the void function 
// pointed to by the constructor argument.
namespace GAMBIT
{
  struct ini_code 
  {
    ini_code(void (*unroll)()) { (*unroll)(); } 
  };
} // end namespace GAMBIT



//
// Handy macros
//
#define CAT(X,Y) X##Y
#define TYPENAME(X) CAT(X,_type)



//
// Macro containing initialization code
//
#define LOAD_LIBRARY                                                        \
namespace GAMBIT                                                          \
{                                                                           \
  namespace Backend                                                       \
  {                                                                         \
    namespace BACKENDNAME                                                 \
    {                                                                        \
      void * pHandle;                                                      \
      void * pSym;                                                         \
                                                                             \
      void loadLibrary()                                                    \
      {                                                                      \
        pHandle = dlmopen(LM_ID_NEWLM, LIBPATH, RTLD_LAZY);                  \
        /* pHandle = dlopen(LIBPATH, RTLD_LAZY); */                          \
                                                                             \
        if(not pHandle) { std::cout << dlerror() << std::endl; }             \
      }                                                                        \
                                                                               \
      /*The code within the void function 'loadLibrary' is executed             \  
        when we create the following instance of the 'ini_code' struct. */     \
      namespace ini                                                           \
      {                                                                        \
        ini_code BACKENDNAME(&loadLibrary);                                  \
      }                                                                        \
                                                                                \
    } /* end namespace BACKENDNAME */                                          \
  } /* end namespace Backend */                                                \
} /* end namespace GAMBIT */                                                   \



//
// Macro for constructing pointers to library variables
//
#define BE_VARIABLE(NAME, TYPE, SYMBOLNAME, POINTERNAME)                    \
namespace GAMBIT                                                           \
{                                                                            \
  namespace Backend                                                        \
  {                                                                          \
    namespace BACKENDNAME                                                  \
    {                                                                        \
      TYPE * POINTERNAME;                                                     \
                                                                             \
      void CAT(constructVarPointer_,NAME)()                                   \
      {                                                                       \
        pSym = dlsym(pHandle, SYMBOLNAME);                                    \
        POINTERNAME = (TYPE*) pSym;                                            \
      }                                                                        \
                                                                               \
      /* The code within the void function 'constructVarPointer'              \
         is executed when we create the following instance of                 \
         the 'ini_code' struct. */                                            \
      namespace ini                                                          \
      {                                                                        \
        ini_code NAME(&CAT(constructVarPointer_,NAME));                        \
      }                                                                        \
                                                                               \
    } /* end namespace BACKENDNAME */                                          \
  } /* end namespace Backend */                                                \
} /* end namespace GAMBIT */                                                   \



//
// Macro for constructing pointers to library functions
//
#define BE_FUNCTION(NAME, TYPE, ARGTYPES, SYMBOLNAME, POINTERNAME)           \
namespace GAMBIT                                                           \
{                                                                            \
  namespace Backend                                                        \
  {                                                                          \
    namespace BACKENDNAME                                                  \
    {                                                                        \
      typedef TYPE (*TYPENAME(NAME))ARGTYPES;                              \
      TYPENAME(NAME) POINTERNAME;                                             \
                                                                             \
      void CAT(constructFuncPointer_,NAME)()                                   \
      {                                                                       \
        pSym = dlsym(pHandle, SYMBOLNAME);                                    \
        POINTERNAME = (TYPENAME(NAME)) pSym;                                  \
      }                                                                       \
                                                                               \
      /* The code within the void function 'constructVarPointer_NAME'         \
         is executed when we create the following instance of                 \
         the 'ini_code' struct. */                                             \
      namespace ini                                                          \
      {                                                                        \
        ini_code NAME(&CAT(constructFuncPointer_,NAME));                        \
      }                                                                        \
                                                                               \
    } /* end namespace BACKENDNAME */                                          \
  } /* end namespace Backend */                                                \
} /* end namespace GAMBIT */                                                   \


#endif // __BACKEND_GENERAL_HPP__
