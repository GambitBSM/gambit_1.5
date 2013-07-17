/* 
 * Another backend test...
 *
 * \author Anders Kvellestad
 */

/* Hack to bypass core */
#define IN_CORE

#include <backend_macros.hpp>

#define LIBPATH      "./lib_fortran.so"
#ifdef BACKENDRENAME
  #define BACKENDNAME BACKENDRENAME
#else
  #define BACKENDNAME LibFortran
#endif
#define VERSION 1.0


LOAD_LIBRARY


//typedef void (*SubroutineType)(const int&);
//typedef float (*FunctionType)(const int&);

BE_FUNCTION(runMe, void, ( float (*)(const int&), const int&), "runme_", "runMe")
BE_FUNCTION(externalFunction, float, (const int&), "externalfunction_", "externalFunction")





namespace GAMBIT
{
  namespace Backends
  {
    namespace BACKENDNAME
    {

      /* Convenience functions go here */

    } /* end namespace BACKENDNAME */                                          
  } /* end namespace Backends */                                                
} /* end namespace GAMBIT */                                                   


//BE_CONV_FUNCTION(awesomenessByAnders, double, "awesomeness")

// Undefine macros to avoid conflict with other backends
#undef LIBPATH 
#undef BACKENDNAME

