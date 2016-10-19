//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Declarations of the class for mathematica wrapper functions
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
#include MATHEMATICA_WSTP_H

#ifndef __mathematica_functions_hpp__
#define __mathematica_functions_hpp__

namespace Gambit
{
  namespace Backends
  {
  
    class mathematica_function
    {
      public:

       // Constructor
       mathematica_function(void *, str);

       // Actual function that takes care of the wrapping
       void wrapper_function();

      protected:

        // Pointer to the WSTP link stablished during loading 
        WSLINK WSlink;

        // Name of the function to the called through WSTP
        str function_name;

    };

  }
}

#endif /* __mathematica_functions_hpp__ */

#endif /* HAVE_MATHEMATICA */
