/*
 * Example of how to use the macros in 'backend_general.hpp' 
 * to set up a backend for a Fortran library
 * 
 * \author Anders Kvellestad
 * \date 2013-03-26
 * 
 * Modified: 2013-04-05
 */

#include "backend_general.hpp"


/* Start stuff inside include braces; this may need to go in another separate header. */

#ifndef __BACKEND_LIBFORTRANCODE_HPP__
#define __BACKEND_LIBFORTRANCODE_HPP__

  /* Struct to be used as type for a specific common block
   * in the Fortran library. */
  struct commonBlock_type
  {
    double dummyDouble;
    int dummyArray[3];
  };

#endif // end stuff in include braces


/* Specify the path to the shared library along with a backend name. */

#define LIBPATH "./libfortrancode.so"
#ifdef BACKENDRENAME
  #define BACKENDNAME BACKENDRENAME
#else
  #define BACKENDNAME LibFortranCode
#endif


/* Load library, obtain pointers to the various library symbols
 * and set up a minimal interface consisting of get/set functions
 * for the variables and function pointers for the functions.
 * 
 * Note 1: Fortran subroutines are treated as void functions.
 * Note 2: The Fortran default is to pass arguments by reference, so 
 *         (non-array) arguments must be specified with 
 *         a trailing ampersand (&). If appropriate, also add 
 *         a 'const' declaration (as seen below). */

LOAD_LIBRARY

BE_VARIABLE(CommonBlock, commonBlock_type, "commonblock_", pCommonBlock)
BE_VARIABLE(TestInt, int, "__fortrancode_MOD_testint", pTestInt)

BE_FUNCTION(printMe, void, (const double[3]), "__fortrancode_MOD_printme")
BE_FUNCTION(total, int, (const int &, const int &), "__fortrancode_MOD_total")

/* We now have access to the following:
 * 
 * Get/set functions:
 * 
 * commonBlock_type GAMBIT::Backend::BACKENDNAME::getCommonBlock()
 * void GAMBIT::Backend::BACKENDNAME::setCommonBlock(commonBlock_type)
 * 
 * Pointers:
 * 
 * GAMBIT::Backend::BACKENDNAME::pCommonBlock   (commonBlock_type *)
 * GAMBIT::Backend::BACKENDNAME::pTestInt       (int *)
 * 
 * Function pointers:
 * 
 * GAMBIT::Backend::BACKENDNAME::printMe   (* void)(double[3] &)
 * GAMBIT::Backend::BACKENDNAME::total     (* int)(int &, int &)     */


/* At this point we have a minimal interface to the loaded library.
 * Any additional convenince functions could be constructed below 
 * using the available pointers. */


namespace GAMBIT                                                           
{                                                                            
  namespace Backend                                                        
  {                                                                          
    namespace BACKENDNAME                                                  
    {                                                                        


      /* Convenience functions go here */


    } // end namespace BACKENDNAME
  } // end namespace Backend
} // end namespace GAMBIT


/* Undefine macros to avoid conflict with other backends */
#undef LIBPATH 
#undef BACKENDNAME


