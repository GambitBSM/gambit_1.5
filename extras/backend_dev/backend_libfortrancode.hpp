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


#define LIBPATH "./libfortrancode.so"
#ifdef BACKENDRENAME
  #define BACKENDNAME BACKENDRENAME
#else
  #define BACKENDNAME LibFortranCode
#endif


/* Load library */
LOAD_LIBRARY


/* Struct to be used as type for a specific common block
 * in the Fortran library */
struct commonBlock_type
{
  double dummyDouble;
  int dummyArray[3];
};

/* Creating pointers to variables, functions and subroutines 
 * (treated as void functions) in the Fortran library.
 * Note that 'single' arguments must be passed by reference (&), 
 * as this is the Fortran default. */

BE_VARIABLE(commonBlock, commonBlock_type, "commonblock_", pCommonBlock)
BE_VARIABLE(testInt, int, "__fortrancode_MOD_testint", pTestInt)

BE_FUNCTION(printMe, void, (double[3]), "__fortrancode_MOD_printme", pPrintMe)
BE_FUNCTION(total, int, (int &, int &), "__fortrancode_MOD_total", pTotal)


/* Setting up the backend interface */

namespace GAMBIT                                                           
{                                                                            
  namespace Backend                                                        
  {                                                                          
    namespace BACKENDNAME                                                  
    {                                                                        

      // variable: testInt
      int getTestInt() { return *pTestInt; }
      void setTestInt(int a) { *pTestInt = a; }
      
      // variable: commonBlock
      commonBlock_type getCommonBlock() { return *pCommonBlock; }
      void setCommonBlock(commonBlock_type a) { *pCommonBlock = a; }

      
      // function: printMe(int)
      void printMe(double a[3]) { pPrintMe(a); }
      
      // function: total(double, double)
      int total(int a, int b) { return pTotal(a,b); }


    } /* end namespace BACKENDNAME */                                          
  } /* end namespace Backend */                                                
} /* end namespace GAMBIT */                                                   


// Undefine macros to avoid conflict with other backends
#undef LIBPATH 
#undef BACKENDNAME


