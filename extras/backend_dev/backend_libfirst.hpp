/*
 * Example of how to use the macros in 'backend_general.hpp' 
 * to set up a backend for a specific library
 * 
 * \author Anders Kvellestad
 * \date 2013-03-26
 * 
 */

#include "backend_general.hpp"


//
// ------ Specific backend code ----------
//


// Define some macros needed by backend_general.hpp
#define LIBPATH      "./libfirst.so"
#ifdef BACKENDRENAME
  #define BACKENDNAME BACKENDRENAME
#else
  #define BACKENDNAME LibFirst
#endif

// The following macro loads the library in LIBPATH 
// when this header file is included somewhere.

LOAD_LIBRARY

// We now have a void pointer GAMBIT::Backend::BACKENDNAME::pHandle
// that points to the library. 
//
// This pointer is now used to get specific typed pointers 
// to the various symbols in the library.
// 
// First we use a macro to set up pointers to library variables:
//
// syntax: BE_VARIABLE([choose variable name], [type], "[exact symbol name]", [choose pointer name])

BE_VARIABLE(someInt, int, "someInt", pSomeInt)
BE_VARIABLE(someDouble, double, "someDouble", pSomeDouble)

// We now have the following pointers:
// GAMBIT::Backend::BACKENDNAME::pSomeInt     (int *)
// GAMBIT::Backend::BACKENDNAME::pSomeDouble  (double *)
//
// We use a similar macro to get pointers to library functions:
//
// syntax: BE_VARIABLE([choose variable name], [type], [arguement types], "[exact symbol name]", [choose pointer name])

BE_FUNCTION(initialize, void, (int), "_Z10initializei", pInitialize)
BE_FUNCTION(someFunction, void, (), "_Z12someFunctionv", pSomeFunction)
BE_FUNCTION(returnResult, double, (), "_Z12returnResultv", pReturnResult)

// We now have the following pointers:
// GAMBIT::Backend::BACKENDNAME::pInitialize       (void *)(int)
// GAMBIT::Backend::BACKENDNAME::pSomeFunction     (void *)()
// GAMBIT::Backend::BACKENDNAME::pReturnResult     (double *)()


//
// At this point we have access to a set of pointers to the library symbols, 
// and can construct any kind of interface we want. 
//
// The simplest example would be to simply construct get/set functions for the
// variables and simple wrappers for the functions:
//

namespace GAMBIT                                                           
{                                                                            
  namespace Backend                                                        
  {                                                                          
    namespace BACKENDNAME                                                  
    {                                                                        

      // variable: someInt
      int getSomeInt() { return *pSomeInt; }

      void setSomeInt(int a) { *pSomeInt = a; }

      // variable: someDouble
      double getSomeDouble() { return *pSomeDouble; }

      void setSomeDouble(double a) { *pSomeDouble = a; }


      // function: initialize()
      void initialize(int a) { pInitialize(a); }

      // function: someFunction()
      void someFunction() { pSomeFunction(); }

      // function: returnResult()
      double returnResult() { return pReturnResult(); }


    } /* end namespace BACKENDNAME */                                          
  } /* end namespace Backend */                                                
} /* end namespace GAMBIT */                                                   


// Undefine macros to avoid conflict with other backends
#undef LIBPATH 
#undef BACKENDNAME

