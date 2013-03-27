/*
 * Example of how to use the macros in 'backend_general.hpp' 
 * to set up a backend for a specific library
 * 
 * \author Anders Kvellestad
 * \date 2013-03-26
 * 
 */

#ifndef __BACKEND_LIBFIRST_COPY_HPP__
#define __BACKEND_LIBFIRST_COPY_HPP__

#include "backend_general.hpp"


//
// ------ Specific backend code ----------
//


#define LIBPATH      "./libfirst.so"
#define BACKENDNAME LibFirstCopy

LOAD_LIBRARY

BE_VARIABLE(someInt, int, "someInt", pSomeInt)
BE_VARIABLE(someDouble, double, "someDouble", pSomeDouble)

BE_FUNCTION(initialize, void, (int), "_Z10initializei", pInitialize)
BE_FUNCTION(someFunction, void, (), "_Z12someFunctionv", pSomeFunction)
BE_FUNCTION(returnResult, double, (), "_Z12returnResultv", pReturnResult)


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



#undef LIBPATH 
#undef BACKENDNAME

#endif /* __BACKEND_LIBFIRST_COPY_HPP__ */
