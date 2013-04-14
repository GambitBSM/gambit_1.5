/*
 * General macros for loading a shared library
 * and constructing pointers to the variables and functions 
 * within the library. 
 * 
 * \author Anders Kvellestad
 * \date 2013-03-26
 *  
 * Modified: 2013-04-05
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
#define STRINGIFY(X) STRINGIFY2(X)
#define STRINGIFY2(X) #X



//
// Macro containing initialization code
//
#define LOAD_LIBRARY                                                        \
namespace GAMBIT                                                          \
{                                                                           \
  namespace Backend                                                        \
  {                                                                         \
    namespace BACKENDNAME                                                  \
    {                                                                       \
                                                                            \
      void * pHandle;                                                      \
      void * pSym;                                                         \
                                                                            \
      void loadLibrary()                                                   \
      {                                                                     \
        pHandle = dlmopen(LM_ID_NEWLM, LIBPATH, RTLD_LAZY);                 \
        if(not pHandle) { std::cout << dlerror() << std::endl; }           \
      }                                                                     \
                                                                            \
      /*The code within the void function 'loadLibrary' is executed         \  
        when we create the following instance of the 'ini_code' struct. */  \
      namespace ini                                                        \
      {                                                                     \
        ini_code BACKENDNAME(&loadLibrary);                                \
      }                                                                     \
                                                                            \
    } /* end namespace BACKENDNAME */                                       \
  } /* end namespace Backend */                                             \
} /* end namespace GAMBIT */                                                \



//
// Macro for constructing pointers to library variables
// and defining simple 'get/set' functions
//
#define BE_VARIABLE(NAME, TYPE, SYMBOLNAME, POINTERNAME)                     \
namespace GAMBIT                                                           \
{                                                                            \
  namespace Backend                                                         \
  {                                                                          \
    namespace BACKENDNAME                                                   \
    {                                                                        \
                                                                             \
      TYPE * POINTERNAME;                                                    \
                                                                             \
      void CAT(constructVarPointer_,NAME)()                                 \
      {                                                                      \
        pSym = dlsym(pHandle, SYMBOLNAME);                                   \
        POINTERNAME = (TYPE*) pSym;                                          \
      }                                                                      \
                                                                             \
      /* The code within the void function 'constructVarPointer_NAME'        \
         is executed when we create the following instance of                \
         the 'ini_code' struct. */                                           \
      namespace ini                                                         \
      {                                                                      \
        ini_code NAME(&CAT(constructVarPointer_,NAME));                     \
      }                                                                      \
                                                                             \
      /* Construct 'get' function */                                         \
      TYPE CAT(get,NAME)() { return *POINTERNAME; }                         \
                                                                             \
      /* Construct 'set' function */                                         \
      void CAT(set,NAME)(TYPE a) { *POINTERNAME = a; }                      \
                                                                             \
    } /* end namespace BACKENDNAME */                                        \
  } /* end namespace Backend */                                              \
} /* end namespace GAMBIT */                                                 \



//
// Macro for constructing pointers to library functions
//
#define BE_FUNCTION(NAME, TYPE, ARGTYPES, SYMBOLNAME)                         \
namespace GAMBIT                                                            \
{                                                                             \
  namespace Backend                                                         \
  {                                                                           \
    namespace BACKENDNAME                                                  \
    {                                                                         \
                                                                              \
      /* Define a type NAME_type to be a suitable function pointer. */        \
      typedef TYPE (*TYPENAME(NAME))ARGTYPES;                                \
      /* Declare a pointer NAME of type NAME_type */                          \
      TYPENAME(NAME) NAME;                                                    \
                                                                              \
      void CAT(constructFuncPointer_,NAME)()                                 \
      {                                                                       \
        /* Obtain a void pointer (pSym) to the library symbol. */             \
        pSym = dlsym(pHandle, SYMBOLNAME);                                    \
        /* Convert it to type (NAME_type) and assign it to pointer NAME. */  \
        NAME = (TYPENAME(NAME)) pSym;                                         \
      }                                                                       \
                                                                              \
      /* The code within the void function 'constructVarPointer_NAME'         \
         is executed when we create the following instance of                 \
         the 'ini_code' struct. */                                            \
      namespace ini                                                          \
      {                                                                       \
        ini_code NAME(&CAT(constructFuncPointer_,NAME));                     \
      }                                                                       \
                                                                              \
    } /* end namespace BACKENDNAME */                                         \
  } /* end namespace Backend */                                               \
} /* end namespace GAMBIT */                                                  \


#endif // __BACKEND_GENERAL_HPP__
