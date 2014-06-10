/* 
 * Example of how to use the macros in 'backend_general.hpp' 
 * to set up a backend for a specific library.
 * 
 * \author Anders Kvellestad
 * \date 2013-03-26  
 * 
 * Modified: 2013-04-05, Pat Scott 2013-06-23
 */


/* A small hack to avoid having to run everything with the Core */
#define IN_CORE



#include "backend_macros.hpp"

/* Specify the path to the shared library along with a backend name. */

#define LIBPATH      "./libfcode.so"
#ifdef BACKENDRENAME
  #define BACKENDNAME BACKENDRENAME
#else
  #define BACKENDNAME LibFcode
#endif


/* The following macro loads the library (using dlmopen) in LIBPATH 
 * when this header file is included somewhere. */

LOAD_LIBRARY


/* Syntax for BE_FUNCTION:
 * BE_FUNCTION([choose function name], [type], [arguement types], "[exact symbol name]", "[choose capability name]")
 * 
 * The last argument (capability name) is optional. 
 * If left out (as done below) it will default to "[backend name]_[function name]_capability"
 * (e.g. "LibFirst_initialize_capability")
 *
 * Note: Fortran uses pass-by-reference as default, so the argument list in BE_FUNCTION must contain '&' 
 * when referring to Fortran codes.
 */

BE_FUNCTION(addOneSubr, void, (int&), "fcode_mp_addonesubr_")
BE_FUNCTION(addOneFunc, int, (int&), "fcode_mp_addonefunc_")


/* Syntax for BE_VARIABLE:
 * BE_VARIABLE([choose variable name], [type], "[exact symbol name]")  */

//BE_VARIABLE(SomeInt, int, "someInt")




namespace GAMBIT
{
  namespace Backends
  {
    namespace BACKENDNAME
    {

      /* Convenience functions go here */

    } /* end namespace BACKENDNAME */                                          
  } /* end namespace Backend */                                                
} /* end namespace GAMBIT */                                                   


// Undefine macros to avoid conflict with other backends
#undef LIBPATH 
#undef BACKENDNAME

