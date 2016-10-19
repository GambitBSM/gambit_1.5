//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Definitions of the class for mathematica wrapper functions
///
///  ***********************************************
///
///  Authors
///  =======
///
///  (add name and date if you modify)
///
///  \author Tomas Gonzalo
///          (t.e.gonzalo@fys.uio.no)
///  \date 2016 Oct
///
///  ***********************************************

#ifdef HAVE_MATHEMATICA
#include "mathematica_functions.hpp"

namespace Gambit
{
   namespace Backends
   {
  
     // Constructor
     mathematica_function::mathematica_function(void *pHandle, str symbol_name) :
       WSlink((WSLINK) pHandle),
       function_name(symbol_name);

     // Actual function that takes care of the wrapping
     void mathematica_function::wrapper_function()
     {
       cout << "do something" << endl;

     }
  }

}

#endif /* HAVE_MATHEMATICA */
